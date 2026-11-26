//! Real SQLite/WAL and durable receipt fixtures. Payloads intentionally remain
//! opaque native/archive bytes: a database backup does not grant molecular state.
use collagene_lab_library::{
    backup::{self, Policy},
    mutations::*,
    *,
};
use std::{
    fs,
    path::PathBuf,
    sync::atomic::{AtomicU64, Ordering},
};
struct Temp(PathBuf);
impl Temp {
    fn new() -> Self {
        static N: AtomicU64 = AtomicU64::new(0);
        let p = std::env::temp_dir().join(format!(
            "collagene-backup-fixture-{}-{}-{}",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos(),
            N.fetch_add(1, Ordering::Relaxed)
        ));
        let mut b = fs::DirBuilder::new();
        #[cfg(unix)]
        {
            use std::os::unix::fs::DirBuilderExt;
            b.mode(0o700);
        }
        b.create(&p).unwrap();
        Self(p)
    }
    fn path(&self, name: &str) -> PathBuf {
        self.0.join(name)
    }
}
impl Drop for Temp {
    fn drop(&mut self) {
        let _ = fs::remove_dir_all(&self.0);
    }
}
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
fn mutation(db: &Library, operation: MutationOperation) -> Mutation {
    Mutation {
        id: db.new_mutation_id().unwrap(),
        operation,
    }
}
#[test]
fn online_wal_snapshot_restores_every_revision_folder_receipt_and_private_checkpoint_after_restart()
{
    let t = Temp::new();
    let path = t.path("live.sqlite");
    let mut db = Library::open(&path).unwrap();
    let lab = db.create_lab("Lab fixture").unwrap();
    let access = db.trusted_host_access(&lab.id).unwrap();
    let folder = db.create_folder(&access, None, "Plasmids").unwrap();
    let nested = db
        .create_folder(&access, Some(&folder.id), "Reviewed")
        .unwrap();
    let import = mutation(
        &db,
        MutationOperation::ImportDocument {
            folder: Some(nested.id.clone()),
            name: "Actual stored record".into(),
            payload: payload("opaque native revision one"),
        },
    );
    let key = CheckpointKey::repository("durable-outbox");
    let checkpoint = CheckpointValue {
        version: 1,
        bytes: b"private local pending mutation fixture".to_vec(),
    };
    let applied = db
        .apply_reliable_with(&access, &import, |_| {
            Ok(vec![CheckpointWrite {
                key: key.clone(),
                expected: None,
                value: Some(checkpoint.clone()),
            }])
        })
        .unwrap();
    let MutationOutcome::Applied(Applied::Imported(doc)) = &applied.receipt.outcome else {
        panic!("import")
    };
    let original = doc.clone();
    let second = db
        .save(
            &access,
            &doc.id,
            &doc.head.id,
            &payload("opaque native revision two"),
        )
        .unwrap();
    let conflict = mutation(
        &db,
        MutationOperation::SaveRevision {
            id: doc.id.clone(),
            base_revision: doc.head.id.clone(),
            payload: payload("stale write"),
        },
    );
    let rejected = db.apply_reliable(&access, &conflict).unwrap();
    assert!(matches!(
        rejected.receipt.outcome,
        MutationOutcome::Rejected(_)
    ));
    let mut writer = Library::open(&path).unwrap();
    let other = writer
        .create_folder(&access, None, "Second live connection")
        .unwrap();
    assert!(fs::metadata(t.path("live.sqlite-wal")).unwrap().len() > 0);
    let policy = Policy::single_lab(lab.id.clone());
    let output = t.path("online.cglbak");
    let manifest = backup::create(&db, &output, &policy).unwrap();
    assert_eq!(backup::verify(&output, &policy).unwrap(), manifest);
    assert_eq!(manifest.table_rows["revisions"], 2);
    assert_eq!(manifest.table_rows["cg_reliable_mutations"], 2);
    let root = t.path("restored");
    assert_eq!(
        backup::restore_new(&output, &root, &policy).unwrap(),
        manifest
    );
    let mut restored = Library::open(root.join("library.sqlite")).unwrap();
    let restored_access = restored.trusted_host_access(&lab.id).unwrap();
    assert_eq!(
        restored.folder(&restored_access, &folder.id).unwrap(),
        folder
    );
    assert_eq!(
        restored.folder(&restored_access, &nested.id).unwrap(),
        nested
    );
    assert_eq!(restored.folder(&restored_access, &other.id).unwrap(), other);
    assert_eq!(
        restored.document(&restored_access, &original.id).unwrap(),
        second
    );
    assert_eq!(
        restored
            .payload(&restored_access, &original.id, &original.head.id)
            .unwrap()
            .bytes,
        b"opaque native revision one"
    );
    assert_eq!(
        restored
            .payload(&restored_access, &original.id, &second.head.id)
            .unwrap()
            .bytes,
        b"opaque native revision two"
    );
    let replay = restored.apply_reliable(&restored_access, &import).unwrap();
    assert!(replay.replayed);
    assert_eq!(replay.receipt, applied.receipt);
    let replay = restored
        .apply_reliable(&restored_access, &conflict)
        .unwrap();
    assert!(replay.replayed);
    assert_eq!(replay.receipt, rejected.receipt);
    assert_eq!(
        restored
            .read_checkpoint(&restored_access, &key)
            .unwrap()
            .unwrap()
            .value,
        checkpoint
    );
    // No global live database is replaced: restoring and editing a new copy cannot
    // advance the source head or modify the immutable backup.
    restored
        .save(
            &restored_access,
            &second.id,
            &second.head.id,
            &payload("new private restored head"),
        )
        .unwrap();
    assert_eq!(db.document(&access, &second.id).unwrap(), second);
    assert_eq!(backup::verify(&output, &policy).unwrap(), manifest);
}
#[test]
fn corrupt_truncated_trailing_or_foreign_lab_backup_never_creates_restore_destination() {
    let t = Temp::new();
    let mut db = Library::open_in_memory().unwrap();
    let lab = db.create_lab("One").unwrap();
    let a = db.trusted_host_access(&lab.id).unwrap();
    db.create_document(&a, None, "p", &payload("archive bytes"))
        .unwrap();
    let p = Policy::single_lab(lab.id.clone());
    let output = t.path("source.cglbak");
    backup::create(&db, &output, &p).unwrap();
    let original = fs::read(&output).unwrap();
    let mut bad = original.clone();
    let last = bad.len() - 1;
    bad[last] ^= 1;
    let mut bad_manifest = original.clone();
    bad_manifest[44] ^= 1;
    let mut trailing = original.clone();
    trailing.push(0);
    let mut future = original.clone();
    future[7] = 2;
    for (i, bytes) in [
        bad,
        bad_manifest,
        original[..original.len() - 1].to_vec(),
        trailing,
        future,
    ]
    .into_iter()
    .enumerate()
    {
        let path = t.path(&format!("broken-{i}"));
        fs::write(&path, bytes).unwrap();
        #[cfg(unix)]
        {
            use std::os::unix::fs::PermissionsExt;
            fs::set_permissions(&path, fs::Permissions::from_mode(0o600)).unwrap();
        }
        let dest = t.path(&format!("reject-{i}"));
        assert!(backup::verify(&path, &p).is_err());
        assert!(backup::restore_new(&path, &dest, &p).is_err());
        assert!(!dest.exists());
    }
    let wrong = Policy::single_lab(LabId::parse("lab-11111111111111111111111111111111").unwrap());
    assert!(backup::verify(&output, &wrong).is_err());
    assert!(backup::restore_new(&output, &t.path("wrong-lab"), &wrong).is_err());
    assert_eq!(fs::read(&output).unwrap(), original);
    let other = db
        .create_lab("Other lab must not be silently bundled")
        .unwrap();
    assert!(backup::create(&db, &t.path("unexpected.cglbak"), &p).is_err());
    let full = Policy {
        expected_labs: vec![lab.id, other.id],
        max_database_bytes: backup::DEFAULT_MAX_DATABASE_BYTES,
    };
    assert_eq!(
        backup::create(&db, &t.path("explicit-multiple.cglbak"), &full)
            .unwrap()
            .labs
            .len(),
        2
    );
}
#[test]
fn backup_and_restore_are_immutable_no_overwrites_and_enforce_database_budget() {
    let t = Temp::new();
    let mut db = Library::open_in_memory().unwrap();
    let lab = db.create_lab("Bounded").unwrap();
    let p = Policy::single_lab(lab.id);
    let file = t.path("snapshot");
    backup::create(&db, &file, &p).unwrap();
    let original = fs::read(&file).unwrap();
    assert!(backup::create(&db, &file, &p).is_err());
    assert_eq!(fs::read(&file).unwrap(), original);
    let tiny = Policy {
        max_database_bytes: 100,
        ..p.clone()
    };
    assert!(backup::create(&db, &t.path("too-big"), &tiny).is_err());
    assert!(!t.path("too-big").exists());
    assert!(backup::verify(&file, &tiny).is_err());
    let existing = t.path("existing");
    fs::create_dir(&existing).unwrap();
    assert!(backup::restore_new(&file, &existing, &p).is_err());
    assert!(!existing.join("library.sqlite").exists());
    let occupied = t.path("occupied");
    fs::create_dir(&occupied).unwrap();
    fs::write(occupied.join("library.sqlite"), b"must not change").unwrap();
    assert!(backup::restore_new(&file, &occupied, &p).is_err());
    assert_eq!(
        fs::read(occupied.join("library.sqlite")).unwrap(),
        b"must not change"
    );
    assert!(!fs::read_dir(&t.0).unwrap().any(|r| {
        r.unwrap()
            .file_name()
            .to_string_lossy()
            .starts_with(".collagene-backup-")
    }));
}
#[test]
fn unsupported_schema_and_corrupted_blob_bindings_refuse_before_publication() {
    let t = Temp::new();
    let source = t.path("live.sqlite");
    let mut db = Library::open(&source).unwrap();
    let lab = db.create_lab("Schema").unwrap();
    let access = db.trusted_host_access(&lab.id).unwrap();
    db.create_document(&access, None, "p", &payload("original opaque bytes"))
        .unwrap();
    let p = Policy::single_lab(lab.id);
    let raw = rusqlite::Connection::open(&source).unwrap();
    raw.execute_batch("CREATE TABLE unexpected_private_contract(x TEXT)")
        .unwrap();
    assert!(backup::create(&db, &t.path("foreign"), &p).is_err());
    assert!(!t.path("foreign").exists());
    let trigger: String = raw
        .query_row(
            "SELECT sql FROM sqlite_master WHERE name='immutable_revision_update'",
            [],
            |r| r.get(0),
        )
        .unwrap();
    raw.execute_batch("DROP TABLE unexpected_private_contract; DROP TRIGGER immutable_revision_update; UPDATE revisions SET archive=zeroblob(archive_size)").unwrap();
    raw.execute_batch(&trigger).unwrap();
    assert!(backup::create(&db, &t.path("corrupt"), &p).is_err());
    assert!(!t.path("corrupt").exists());
}
#[cfg(unix)]
#[test]
fn unix_private_modes_and_symlink_refusal_do_not_follow_or_overwrite_foreign_paths() {
    use std::os::unix::fs::{PermissionsExt, symlink};
    let t = Temp::new();
    let mut db = Library::open_in_memory().unwrap();
    let lab = db.create_lab("Private").unwrap();
    let p = Policy::single_lab(lab.id);
    let file = t.path("good");
    backup::create(&db, &file, &p).unwrap();
    assert_eq!(
        fs::metadata(&file).unwrap().permissions().mode() & 0o777,
        0o600
    );
    let dest = t.path("restored");
    backup::restore_new(&file, &dest, &p).unwrap();
    assert_eq!(
        fs::metadata(&dest).unwrap().permissions().mode() & 0o777,
        0o700
    );
    assert_eq!(
        fs::metadata(dest.join("library.sqlite"))
            .unwrap()
            .permissions()
            .mode()
            & 0o777,
        0o600
    );
    let link = t.path("link");
    symlink(&file, &link).unwrap();
    assert!(backup::verify(&link, &p).is_err());
    assert!(backup::create(&db, &link, &p).is_err());
    fs::set_permissions(&file, fs::Permissions::from_mode(0o644)).unwrap();
    assert!(backup::verify(&file, &p).is_err());
}
