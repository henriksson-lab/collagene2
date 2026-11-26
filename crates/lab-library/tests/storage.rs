use collagene_lab_library::*;
use std::{
    path::PathBuf,
    sync::{
        Arc, Barrier,
        atomic::{AtomicU64, Ordering},
    },
};

// Deliberately opaque fixtures: native decoding/sharing approval belongs to the
// trusted host. These tests exercise database integrity, never molecular authority.
fn payload(text: &str) -> Payload {
    let bytes = text.as_bytes().to_vec();
    Payload {
        kind: PayloadKind::NativeDocument,
        archive_hash: Digest::parse(blake3::hash(&bytes).to_hex().to_string()).unwrap(),
        bytes,
        selection: None,
        sharing: None,
    }
}
fn workspace_payload(text: &str) -> Payload {
    let mut p = payload(text);
    p.kind = PayloadKind::WorkspaceArchive;
    p.selection = Some(NativeSelection {
        document_id: "local-doc-1".into(),
        revision: Digest::parse(
            blake3::hash(b"selected native revision")
                .to_hex()
                .to_string(),
        )
        .unwrap(),
        head_event: Some("event-1234".into()),
    });
    p.sharing = Some(WorkspaceSharing::DocumentClosure);
    p
}
fn fixture() -> (Library, LabAccess, LabAccess) {
    let mut db = Library::open_in_memory().unwrap();
    let a = db.create_lab("Lab A").unwrap();
    let b = db.create_lab("Lab B").unwrap();
    let a = db.trusted_host_access(&a.id).unwrap();
    let b = db.trusted_host_access(&b.id).unwrap();
    (db, a, b)
}
struct TempDb(PathBuf);
impl TempDb {
    fn new() -> Self {
        static NEXT: AtomicU64 = AtomicU64::new(0);
        let dir = std::env::temp_dir().join(format!(
            "collagene-lab-storage-{}-{}-{}",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos(),
            NEXT.fetch_add(1, Ordering::Relaxed)
        ));
        std::fs::create_dir(&dir).unwrap();
        Self(dir)
    }
    fn file(&self) -> PathBuf {
        self.0.join("library.sqlite")
    }
}
impl Drop for TempDb {
    fn drop(&mut self) {
        let _ = std::fs::remove_dir_all(&self.0);
    }
}

#[test]
fn folders_are_lab_scoped_acyclic_and_cannot_delete_data() {
    let (mut db, a, b) = fixture();
    let root = db.create_folder(&a, None, "Plasmids").unwrap();
    assert_eq!(root.generation, 1);
    let child = db.create_folder(&a, Some(&root.id), "Vectors").unwrap();
    let doc = db
        .create_document(&a, Some(&child.id), "pExample", &payload("archive-v1"))
        .unwrap();
    let before = db.lab(&a).unwrap().catalog_generation;
    assert!(matches!(db.folder(&b, &root.id), Err(Error::NotFound)));
    assert!(matches!(db.document(&b, &doc.id), Err(Error::NotFound)));
    assert!(matches!(
        db.create_folder(&b, Some(&root.id), "foreign"),
        Err(Error::NotFound)
    ));
    assert!(matches!(
        db.move_folder(&a, &root.id, root.generation, Some(&child.id)),
        Err(Error::Invalid("folder cycle"))
    ));
    assert!(matches!(
        db.move_folder(&b, &root.id, root.generation, None),
        Err(Error::NotFound)
    ));
    assert!(matches!(
        db.delete_empty_folder(&a, &root.id, root.generation),
        Err(Error::NotEmpty)
    ));
    assert!(matches!(
        db.delete_empty_folder(&a, &child.id, child.generation),
        Err(Error::NotEmpty)
    ));
    assert_eq!(db.lab(&a).unwrap().catalog_generation, before);
    let moved = db.move_document(&a, &doc.id, doc.generation, None).unwrap();
    assert_eq!(moved.head, doc.head);
    assert!(moved.folder.is_none());
    db.delete_empty_folder(&a, &child.id, child.generation)
        .unwrap();
    db.delete_empty_folder(&a, &root.id, root.generation)
        .unwrap();
    assert_eq!(
        db.payload(&a, &doc.id, &doc.head.id).unwrap().bytes,
        b"archive-v1"
    );
}
#[test]
fn names_metadata_cas_and_cross_type_sibling_uniqueness_are_atomic() {
    let (mut db, a, b) = fixture();
    let folder = db.create_folder(&a, None, "Vectors").unwrap();
    assert!(matches!(
        db.create_document(&a, None, "Vectors", &payload("a")),
        Err(Error::AlreadyExists)
    ));
    db.create_folder(&b, None, "Vectors").unwrap();
    let renamed = db.rename_folder(&a, &folder.id, 1, "Backbones").unwrap();
    assert_eq!(renamed.generation, 2);
    assert!(matches!(
        db.rename_folder(&a, &folder.id, 1, "Stale"),
        Err(Error::MetadataConflict { actual: 2, .. })
    ));
    let doc = db
        .create_document(&a, None, "Insert", &payload("a"))
        .unwrap();
    let renamed = db
        .rename_document(&a, &doc.id, doc.generation, "Insert 2")
        .unwrap();
    assert_eq!(renamed.head.id, doc.head.id);
    assert!(matches!(
        db.move_document(&a, &doc.id, doc.generation, Some(&folder.id)),
        Err(Error::MetadataConflict { .. })
    ));
    for name in ["", " ", "bad\nname", "bad\0name"] {
        assert!(db.create_folder(&a, None, name).is_err());
    }
    assert!(db.create_folder(&a, None, &"é".repeat(257)).is_err());
    let foreign = db.create_folder(&b, None, "foreign").unwrap();
    assert!(matches!(
        db.move_document(&a, &doc.id, renamed.generation, Some(&foreign.id)),
        Err(Error::NotFound)
    ));
    assert_eq!(db.document(&a, &doc.id).unwrap(), renamed);
}
#[test]
fn immutable_cas_history_and_payload_pins_survive_restart() {
    let temp = TempDb::new();
    let (lab, access, id, first, second, metadata) = {
        let mut db = Library::open(temp.file()).unwrap();
        let lab = db.create_lab("Shared lab").unwrap();
        let access = db.trusted_host_access(&lab.id).unwrap();
        let doc = db
            .create_document(
                &access,
                None,
                "pArchive",
                &workspace_payload("opaque native workspace closure"),
            )
            .unwrap();
        let first = doc.head.id.clone();
        let next = workspace_payload("opaque updated native workspace closure");
        let doc = db.save(&access, &doc.id, &first, &next).unwrap();
        let second = doc.head.id.clone();
        assert_ne!(first, second);
        assert_eq!(doc.head.parent.as_ref(), Some(&first));
        assert_eq!(doc.head.number, 2);
        let generation = db.lab(&access).unwrap().catalog_generation;
        assert!(
            matches!(db.save(&access,&doc.id,&first,&payload("stale")),Err(Error::Conflict{actual,..}) if actual==second)
        );
        assert_eq!(db.lab(&access).unwrap().catalog_generation, generation);
        assert_eq!(db.revisions(&access, &doc.id, 0, 20).unwrap().len(), 2);
        (lab, access, doc.id.clone(), first, second, doc)
    };
    let mut db = Library::open(temp.file()).unwrap();
    let reopened = db.trusted_host_access(&lab.id).unwrap();
    assert_eq!(reopened.lab_id(), access.lab_id());
    assert_eq!(db.document(&reopened, &id).unwrap(), metadata);
    let original = db.payload(&reopened, &id, &first).unwrap();
    let current = db.payload(&reopened, &id, &second).unwrap();
    assert_eq!(original.bytes, b"opaque native workspace closure");
    assert_eq!(current.bytes, b"opaque updated native workspace closure");
    assert_eq!(current.selection, workspace_payload("x").selection);
    assert_eq!(current.sharing, Some(WorkspaceSharing::DocumentClosure));
    let page = db.revision_page(&reopened, &id, 0, 1, None).unwrap();
    assert_eq!(page.revisions[0].id, second);
    assert_eq!(page.next_offset, Some(1));
    let older = db
        .revision_page(&reopened, &id, 1, 1, Some(page.catalog_generation))
        .unwrap();
    assert_eq!(older.revisions[0].id, first);
}
#[test]
fn concurrent_writers_do_not_overwrite_or_append_losing_history() {
    let temp = TempDb::new();
    let path = temp.file();
    let (access, doc) = {
        let mut db = Library::open(&path).unwrap();
        let lab = db.create_lab("Lab").unwrap();
        let access = db.trusted_host_access(&lab.id).unwrap();
        let doc = db
            .create_document(&access, None, "pCAS", &payload("original"))
            .unwrap();
        (access, doc)
    };
    let barrier = Arc::new(Barrier::new(2));
    let handles = ["writer-a", "writer-b"].map(|text| {
        let path = path.clone();
        let access = access.clone();
        let doc = doc.clone();
        let barrier = barrier.clone();
        std::thread::spawn(move || {
            let mut db = Library::open(path).unwrap();
            barrier.wait();
            db.save(&access, &doc.id, &doc.head.id, &payload(text))
        })
    });
    let results = handles.map(|handle| handle.join().unwrap());
    assert_eq!(results.iter().filter(|r| r.is_ok()).count(), 1);
    assert_eq!(
        results
            .iter()
            .filter(|r| matches!(r, Err(Error::Conflict { .. })))
            .count(),
        1
    );
    let db = Library::open(&path).unwrap();
    let revisions = db.revisions(&access, &doc.id, 0, 20).unwrap();
    assert_eq!(revisions.len(), 2);
    assert_eq!(revisions[0].parent.as_ref(), Some(&doc.head.id));
    assert_eq!(
        db.payload(&access, &doc.id, &doc.head.id).unwrap().bytes,
        b"original"
    );
}
#[test]
fn bounded_snapshot_pages_and_scoped_literal_name_search() {
    let (mut db, a, b) = fixture();
    assert_eq!(
        db.browse(&a, None, 0, 20, None).unwrap().catalog_generation,
        0
    );
    let folder = db.create_folder(&a, None, "Alpha folder").unwrap();
    let child = db.create_folder(&a, Some(&folder.id), "Nested").unwrap();
    let root = db
        .create_document(&a, None, "Alpha root", &payload("a"))
        .unwrap();
    let nested = db
        .create_document(&a, Some(&child.id), "Alpha 100%_nested", &payload("b"))
        .unwrap();
    db.create_document(&b, None, "Alpha foreign", &payload("c"))
        .unwrap();
    let first = db.browse(&a, None, 0, 1, None).unwrap();
    assert_eq!(first.entries.len(), 1);
    assert_eq!(first.next_offset, Some(1));
    assert!(matches!(
        db.browse(&a, None, 1, 1, None),
        Err(Error::Invalid("continuation snapshot"))
    ));
    assert_eq!(
        db.browse(&a, None, 1, 1, Some(first.catalog_generation))
            .unwrap()
            .entries,
        vec![Entry::Document(root.clone())]
    );
    let all = db
        .search(&a, "ALPHA", &SearchScope::AllFolders, 0, 20, None)
        .unwrap();
    assert_eq!(all.entries.len(), 3);
    let direct = db
        .search(
            &a,
            "Alpha",
            &SearchScope::Folder {
                folder: Some(folder.id.clone()),
                recursive: false,
            },
            0,
            20,
            None,
        )
        .unwrap();
    assert!(direct.entries.is_empty());
    let recursive = db
        .search(
            &a,
            "100%_",
            &SearchScope::Folder {
                folder: Some(folder.id),
                recursive: true,
            },
            0,
            20,
            None,
        )
        .unwrap();
    assert_eq!(recursive.entries, vec![Entry::Document(nested)]);
    db.rename_document(&a, &root.id, root.generation, "Changed")
        .unwrap();
    assert!(matches!(
        db.browse(&a, None, 1, 1, Some(first.catalog_generation)),
        Err(Error::CatalogConflict { .. })
    ));
    assert!(db.browse(&a, None, 0, MAX_PAGE_SIZE + 1, None).is_err());
}
#[test]
fn quotas_bad_hashes_and_ambiguous_workspace_sharing_fail_without_writes() {
    let mut db = Library::open_in_memory_with_options(StorageOptions {
        max_lab_archive_bytes: 8,
    })
    .unwrap();
    let lab = db.create_lab("Lab").unwrap();
    let access = db.trusted_host_access(&lab.id).unwrap();
    let doc = db
        .create_document(&access, None, "p", &payload("1234"))
        .unwrap();
    let next = db
        .save(&access, &doc.id, &doc.head.id, &payload("5678"))
        .unwrap();
    let before = db.lab(&access).unwrap().catalog_generation;
    assert!(matches!(
        db.save(&access, &doc.id, &next.head.id, &payload("x")),
        Err(Error::Limit("lab archive bytes"))
    ));
    let mut bad = payload("wrong");
    bad.archive_hash = payload("different").archive_hash;
    assert!(bad.validate().is_err());
    let mut workspace = workspace_payload("a");
    workspace.sharing = None;
    assert!(workspace.validate().is_err());
    workspace.sharing = Some(WorkspaceSharing::ExplicitWholeWorkspace);
    workspace.selection = None;
    assert!(workspace.validate().is_err());
    assert!(db.create_document(&access, None, "bad", &bad).is_err());
    assert_eq!(db.lab(&access).unwrap().catalog_generation, before);
    assert_eq!(db.revisions(&access, &doc.id, 0, 20).unwrap().len(), 2);
}
#[test]
fn oversized_folder_depth_and_subtree_moves_leave_the_tree_unchanged() {
    let (mut db, a, _) = fixture();
    let mut parent = None;
    for i in 0..MAX_FOLDER_DEPTH {
        let folder = db
            .create_folder(&a, parent.as_ref(), &format!("level-{i}"))
            .unwrap();
        parent = Some(folder.id);
    }
    let before = db.lab(&a).unwrap().catalog_generation;
    assert!(matches!(
        db.create_folder(&a, parent.as_ref(), "too deep"),
        Err(Error::Limit("folder depth"))
    ));
    assert_eq!(db.lab(&a).unwrap().catalog_generation, before);
    let branch = db.create_folder(&a, None, "separate").unwrap();
    db.create_folder(&a, Some(&branch.id), "child").unwrap();
    assert!(matches!(
        db.move_folder(&a, &branch.id, branch.generation, parent.as_ref()),
        Err(Error::Limit("folder depth"))
    ));
    assert!(db.folder(&a, &branch.id).unwrap().parent.is_none());
}
#[test]
fn unrelated_database_is_not_reinterpreted_and_revisions_cannot_be_mutated() {
    let temp = TempDb::new();
    let path = temp.file();
    let raw = rusqlite::Connection::open(&path).unwrap();
    raw.execute("CREATE TABLE unrelated(value TEXT)", [])
        .unwrap();
    drop(raw);
    assert!(matches!(
        Library::open(&path),
        Err(Error::UnsupportedDatabase)
    ));
    let temp = TempDb::new();
    let path = temp.file();
    let (access, doc) = {
        let mut db = Library::open(&path).unwrap();
        let lab = db.create_lab("Lab").unwrap();
        let access = db.trusted_host_access(&lab.id).unwrap();
        let doc = db
            .create_document(&access, None, "p", &payload("protected"))
            .unwrap();
        (access, doc)
    };
    let raw = rusqlite::Connection::open(&path).unwrap();
    assert!(
        raw.execute(
            "UPDATE revisions SET archive=?1 WHERE id=?2",
            rusqlite::params![b"tampered".as_slice(), doc.head.id.as_str()]
        )
        .is_err()
    );
    assert!(
        raw.execute("DELETE FROM revisions WHERE id=?1", [doc.head.id.as_str()])
            .is_err()
    );
    drop(raw);
    let db = Library::open(&path).unwrap();
    assert_eq!(
        db.payload(&access, &doc.id, &doc.head.id).unwrap().bytes,
        b"protected"
    );
}
