use collagene_lab_library::{mutations::*, *};
use std::{
    path::PathBuf,
    sync::{
        Arc, Barrier,
        atomic::{AtomicUsize, Ordering},
    },
};
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
fn folder(receipt: &Receipt) -> Folder {
    let MutationOutcome::Applied(Applied::CreatedFolder(folder)) = &receipt.outcome else {
        panic!("expected created folder")
    };
    folder.clone()
}
fn document(receipt: &Receipt) -> Document {
    let MutationOutcome::Applied(Applied::Imported(document)) = &receipt.outcome else {
        panic!("expected imported document")
    };
    document.clone()
}
fn value(text: &str) -> CheckpointValue {
    CheckpointValue {
        version: 1,
        bytes: text.as_bytes().to_vec(),
    }
}
fn checkpoint(name: &str, expected: Option<Digest>, text: &str) -> CheckpointWrite {
    CheckpointWrite {
        key: CheckpointKey::repository(name),
        expected,
        value: Some(value(text)),
    }
}
fn fixture() -> (Library, LabAccess, LabAccess) {
    let mut db = Library::open_in_memory().unwrap();
    let a = db.create_lab("A").unwrap();
    let b = db.create_lab("B").unwrap();
    let a = db.trusted_host_access(&a.id).unwrap();
    let b = db.trusted_host_access(&b.id).unwrap();
    (db, a, b)
}
struct Temp(PathBuf);
impl Temp {
    fn new() -> Self {
        static NEXT: AtomicUsize = AtomicUsize::new(0);
        let path = std::env::temp_dir().join(format!(
            "collagene-reliable-{}-{}-{}",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos(),
            NEXT.fetch_add(1, Ordering::Relaxed)
        ));
        std::fs::create_dir(&path).unwrap();
        Self(path)
    }
    fn file(&self) -> PathBuf {
        self.0.join("library.sqlite")
    }
}
impl Drop for Temp {
    fn drop(&mut self) {
        let _ = std::fs::remove_dir_all(&self.0);
    }
}

#[test]
fn lost_import_reply_replays_original_receipt_after_restart_and_later_edits() {
    let temp = Temp::new();
    let (mut db, access) = {
        let mut db = Library::open(temp.file()).unwrap();
        let access = db.local_repository_access("local", "Local").unwrap();
        (db, access)
    };
    let request = mutation(
        &db,
        MutationOperation::ImportDocument {
            folder: None,
            name: "pOne".into(),
            payload: payload("original"),
        },
    );
    let first = db.apply_reliable(&access, &request).unwrap();
    assert!(!first.replayed);
    let doc = document(&first.receipt);
    db.rename_document(&access, &doc.id, doc.generation, "Later display name")
        .unwrap();
    let generation = db.lab(&access).unwrap().catalog_generation;
    drop(db);
    let mut db = Library::open(temp.file()).unwrap();
    let reopened = db
        .local_repository_access("local", "Ignored restart label")
        .unwrap();
    assert_eq!(access.lab_id(), reopened.lab_id());
    let replay = db.apply_reliable(&reopened, &request).unwrap();
    assert!(replay.replayed);
    assert_eq!(replay.receipt, first.receipt);
    assert_eq!(db.lab(&reopened).unwrap().catalog_generation, generation);
    assert_eq!(
        db.browse(&reopened, None, 0, 20, None)
            .unwrap()
            .entries
            .len(),
        1
    );
    assert_eq!(
        db.lookup_reliable(&reopened, &request.id, &request.intent_digest().unwrap())
            .unwrap(),
        Some(first.receipt)
    );
    assert_eq!(
        db.payload(&reopened, &doc.id, &doc.head.id).unwrap().bytes,
        b"original"
    );
}
#[test]
fn mutation_ids_bind_canonical_intent_and_are_scoped_to_one_lab() {
    let (mut db, a, b) = fixture();
    let request = mutation(
        &db,
        MutationOperation::CreateFolder {
            parent: None,
            name: "Vectors".into(),
        },
    );
    let canonical = br#"{"CreateFolder":{"parent":null,"name":"Vectors"}}"#;
    let mut hash = blake3::Hasher::new();
    hash.update(b"collagene-library-mutation-v1\0");
    hash.update(canonical);
    assert_eq!(
        request.intent_digest().unwrap().as_str(),
        hash.finalize().to_hex().as_str()
    );
    let original = db.apply_reliable(&a, &request).unwrap();
    let mut altered = request.clone();
    altered.operation = MutationOperation::CreateFolder {
        parent: None,
        name: "Other".into(),
    };
    assert!(matches!(
        db.apply_reliable(&a, &altered),
        Err(MutationError::IntentMismatch)
    ));
    assert!(matches!(
        db.lookup_reliable(&a, &request.id, &altered.intent_digest().unwrap()),
        Err(MutationError::IntentMismatch)
    ));
    assert!(
        db.lookup_reliable(&b, &request.id, &request.intent_digest().unwrap())
            .unwrap()
            .is_none()
    );
    let other = db.apply_reliable(&b, &request).unwrap();
    assert_ne!(folder(&original.receipt).id, folder(&other.receipt).id);
    let imported = mutation(
        &db,
        MutationOperation::ImportDocument {
            folder: None,
            name: "p".into(),
            payload: payload("one"),
        },
    );
    db.apply_reliable(&a, &imported).unwrap();
    let mut changed = imported.clone();
    if let MutationOperation::ImportDocument { payload, .. } = &mut changed.operation {
        *payload = crate::payload("two");
    }
    assert_ne!(
        imported.intent_digest().unwrap(),
        changed.intent_digest().unwrap()
    );
    assert!(matches!(
        db.apply_reliable(&a, &changed),
        Err(MutationError::IntentMismatch)
    ));
}
#[test]
fn conflict_receipts_replay_the_original_snapshot_instead_of_a_later_head() {
    let (mut db, a, _) = fixture();
    let doc = db.create_document(&a, None, "p", &payload("one")).unwrap();
    let second = db.save(&a, &doc.id, &doc.head.id, &payload("two")).unwrap();
    let request = mutation(
        &db,
        MutationOperation::SaveRevision {
            id: doc.id.clone(),
            base_revision: doc.head.id.clone(),
            payload: payload("stale"),
        },
    );
    let conflict = db.apply_reliable(&a, &request).unwrap();
    assert!(
        matches!(&conflict.receipt.outcome,MutationOutcome::Rejected(Rejection{code:RejectionCode::Conflict,conflict:Some(MutationConflict::Document{current}),..})if current.head==second.head)
    );
    let third = db
        .save(&a, &doc.id, &second.head.id, &payload("three"))
        .unwrap();
    let replay = db.apply_reliable(&a, &request).unwrap();
    assert_eq!(replay.receipt, conflict.receipt);
    assert!(replay.replayed);
    assert_eq!(db.document(&a, &doc.id).unwrap().head, third.head);
    assert_eq!(db.revisions(&a, &doc.id, 0, 20).unwrap().len(), 3);
    let new_request = Mutation {
        id: db.new_mutation_id().unwrap(),
        operation: request.operation,
    };
    let latest = db.apply_reliable(&a, &new_request).unwrap();
    assert!(
        matches!(latest.receipt.outcome,MutationOutcome::Rejected(Rejection{conflict:Some(MutationConflict::Document{current}),..})if current.head==third.head)
    );
}
#[test]
fn finalizer_and_checkpoint_cas_failures_roll_back_content_quota_and_receipt_together() {
    let mut db = Library::open_in_memory_with_options(StorageOptions {
        max_lab_archive_bytes: 4,
    })
    .unwrap();
    let access = db.local_repository_access("local", "Local").unwrap();
    let request = mutation(
        &db,
        MutationOperation::ImportDocument {
            folder: None,
            name: "p".into(),
            payload: payload("1234"),
        },
    );
    assert!(matches!(
        db.apply_reliable_with(&access, &request, |_| Err(MutationError::Aborted(
            "simulated repository failure".into()
        ))),
        Err(MutationError::Aborted(_))
    ));
    assert!(
        db.lookup_reliable(&access, &request.id, &request.intent_digest().unwrap())
            .unwrap()
            .is_none()
    );
    assert!(
        db.browse(&access, None, 0, 20, None)
            .unwrap()
            .entries
            .is_empty()
    );
    db.commit_checkpoints(&access, &[checkpoint("state", None, "initial")])
        .unwrap();
    let initial = db
        .read_checkpoint(&access, &CheckpointKey::repository("state"))
        .unwrap()
        .unwrap();
    assert!(matches!(
        db.apply_reliable_with(&access, &request, |_| Ok(vec![checkpoint(
            "state",
            None,
            "wrong stale CAS"
        )])),
        Err(MutationError::CheckpointConflict { .. })
    ));
    assert!(
        db.lookup_reliable(&access, &request.id, &request.intent_digest().unwrap())
            .unwrap()
            .is_none()
    );
    let calls = Arc::new(AtomicUsize::new(0));
    let called = calls.clone();
    let committed = db
        .apply_reliable_with(&access, &request, |receipt| {
            called.fetch_add(1, Ordering::SeqCst);
            let doc = document(receipt);
            Ok(vec![checkpoint(
                "state",
                Some(initial.digest.clone()),
                doc.id.as_str(),
            )])
        })
        .unwrap();
    let replay = db
        .apply_reliable_with(&access, &request, |_| {
            calls.fetch_add(1, Ordering::SeqCst);
            panic!("replay must not rerun finalizer")
        })
        .unwrap();
    assert_eq!(committed.receipt, replay.receipt);
    assert_eq!(calls.load(Ordering::SeqCst), 1);
    let retained = db
        .read_checkpoint(&access, &CheckpointKey::repository("state"))
        .unwrap()
        .unwrap();
    assert_eq!(
        retained.value.bytes,
        document(&committed.receipt).id.as_str().as_bytes()
    );
    assert_eq!(db.lab(&access).unwrap().catalog_generation, 1);
}
#[test]
fn private_checkpoint_batches_are_atomic_bounded_and_persistent_without_catalog_mutations() {
    let temp = Temp::new();
    let mut db = Library::open(temp.file()).unwrap();
    let access = db.local_repository_access("local", "Local").unwrap();
    let other = db.local_repository_access("second", "Second").unwrap();
    db.commit_checkpoints(
        &access,
        &[
            checkpoint("state", None, "one"),
            checkpoint("outbox-mut-001", None, "queued"),
        ],
    )
    .unwrap();
    let first = db
        .read_checkpoint(&access, &CheckpointKey::repository("state"))
        .unwrap()
        .unwrap();
    assert!(matches!(
        db.commit_checkpoints(
            &access,
            &[
                checkpoint("state", Some(first.digest.clone()), "must roll back"),
                checkpoint("outbox-mut-001", None, "stale")
            ]
        ),
        Err(MutationError::CheckpointConflict { .. })
    ));
    assert_eq!(
        db.read_checkpoint(&access, &first.key).unwrap(),
        Some(first.clone())
    );
    assert!(db.read_checkpoint(&other, &first.key).unwrap().is_none());
    assert_eq!(db.lab(&access).unwrap().catalog_generation, 0);
    assert!(
        db.commit_checkpoints(
            &access,
            &[CheckpointWrite {
                key: CheckpointKey {
                    namespace: "credentials".into(),
                    key: "secret".into()
                },
                expected: None,
                value: Some(value("not accepted"))
            }]
        )
        .is_err()
    );
    assert!(
        db.commit_checkpoints(
            &access,
            &[checkpoint("same", None, "a"), checkpoint("same", None, "b")]
        )
        .is_err()
    );
    let mut invalid = checkpoint("zero-version", None, "x");
    invalid.value.as_mut().unwrap().version = 0;
    assert!(db.commit_checkpoints(&access, &[invalid]).is_err());
    drop(db);
    let mut db = Library::open(temp.file()).unwrap();
    let access = db.local_repository_access("local", "Local").unwrap();
    assert_eq!(
        db.read_checkpoint(&access, &first.key).unwrap(),
        Some(first)
    );
}
#[test]
fn concurrent_duplicate_creation_returns_one_identity_and_one_receipt() {
    let temp = Temp::new();
    let path = temp.file();
    let mut db = Library::open(&path).unwrap();
    let access = db.local_repository_access("local", "Local").unwrap();
    let request = mutation(
        &db,
        MutationOperation::ImportDocument {
            folder: None,
            name: "Only one".into(),
            payload: payload("data"),
        },
    );
    drop(db);
    let barrier = Arc::new(Barrier::new(2));
    let handles = [0, 1].map(|_| {
        let barrier = barrier.clone();
        let access = access.clone();
        let request = request.clone();
        let path = path.clone();
        std::thread::spawn(move || {
            let mut db = Library::open(path).unwrap();
            barrier.wait();
            db.apply_reliable(&access, &request).unwrap()
        })
    });
    let [a, b] = handles.map(|h| h.join().unwrap());
    assert_eq!(a.receipt, b.receipt);
    assert_ne!(a.replayed, b.replayed);
    let mut db = Library::open(&path).unwrap();
    assert_eq!(
        db.browse(&access, None, 0, 20, None).unwrap().entries.len(),
        1
    );
    assert_eq!(db.lab(&access).unwrap().catalog_generation, 1);
}
#[test]
fn reliable_tree_and_metadata_operations_keep_existing_invariants_and_apis() {
    let (mut db, a, b) = fixture();
    let create = mutation(
        &db,
        MutationOperation::CreateFolder {
            parent: None,
            name: "Parent".into(),
        },
    );
    let parent = folder(&db.apply_reliable(&a, &create).unwrap().receipt);
    let child = db.create_folder(&a, Some(&parent.id), "Child").unwrap();
    for operation in [
        MutationOperation::MoveFolder {
            id: parent.id.clone(),
            expected_generation: parent.generation,
            parent: Some(child.id.clone()),
        },
        MutationOperation::DeleteEmptyFolder {
            id: parent.id.clone(),
            expected_generation: parent.generation,
        },
    ] {
        let request = mutation(&db, operation);
        assert!(matches!(
            db.apply_reliable(&a, &request).unwrap().receipt.outcome,
            MutationOutcome::Rejected(_)
        ));
    }
    let request = mutation(
        &db,
        MutationOperation::RenameFolder {
            id: parent.id.clone(),
            expected_generation: parent.generation,
            name: "Renamed".into(),
        },
    );
    let result = db.apply_reliable(&a, &request).unwrap();
    assert!(matches!(
        result.receipt.outcome,
        MutationOutcome::Applied(Applied::RenamedFolder(_))
    ));
    let foreign = mutation(
        &db,
        MutationOperation::CreateFolder {
            parent: Some(parent.id.clone()),
            name: "Forbidden".into(),
        },
    );
    assert!(matches!(
        db.apply_reliable(&b, &foreign).unwrap().receipt.outcome,
        MutationOutcome::Rejected(Rejection {
            code: RejectionCode::NotFound,
            ..
        })
    ));
    let doc = db.create_document(&a, None, "p", &payload("a")).unwrap();
    let rename = mutation(
        &db,
        MutationOperation::RenameDocument {
            id: doc.id.clone(),
            expected_generation: doc.generation,
            name: "pRenamed".into(),
        },
    );
    let MutationOutcome::Applied(Applied::RenamedDocument(renamed)) =
        db.apply_reliable(&a, &rename).unwrap().receipt.outcome
    else {
        panic!("rename failed")
    };
    let move_doc = mutation(
        &db,
        MutationOperation::MoveDocument {
            id: doc.id.clone(),
            expected_generation: renamed.generation,
            folder: Some(child.id.clone()),
        },
    );
    assert!(matches!(
        db.apply_reliable(&a, &move_doc).unwrap().receipt.outcome,
        MutationOutcome::Applied(Applied::MovedDocument(_))
    ));
    assert_eq!(db.payload(&a, &doc.id, &doc.head.id).unwrap().bytes, b"a");
    let empty = db.create_folder(&a, None, "Empty").unwrap();
    let delete = mutation(
        &db,
        MutationOperation::DeleteEmptyFolder {
            id: empty.id.clone(),
            expected_generation: empty.generation,
        },
    );
    let first = db.apply_reliable(&a, &delete).unwrap();
    let again = db.apply_reliable(&a, &delete).unwrap();
    assert_eq!(first.receipt, again.receipt);
    assert!(matches!(db.folder(&a, &empty.id), Err(Error::NotFound)));
}

#[test]
fn concurrent_distinct_save_nonces_preserve_one_winner_and_replayable_loser() {
    let temp = Temp::new();
    let path = temp.file();
    let mut db = Library::open(&path).unwrap();
    let access = db.local_repository_access("local", "Local").unwrap();
    let original = db
        .create_document(&access, None, "pCAS", &payload("original"))
        .unwrap();
    let requests = ["left edit", "right edit"].map(|text| {
        mutation(
            &db,
            MutationOperation::SaveRevision {
                id: original.id.clone(),
                base_revision: original.head.id.clone(),
                payload: payload(text),
            },
        )
    });
    drop(db);
    let barrier = Arc::new(Barrier::new(2));
    let handles = requests.clone().map(|request| {
        let barrier = barrier.clone();
        let access = access.clone();
        let path = path.clone();
        std::thread::spawn(move || {
            let mut db = Library::open(path).unwrap();
            barrier.wait();
            db.apply_reliable(&access, &request).unwrap()
        })
    });
    let receipts = handles.map(|h| h.join().unwrap().receipt);
    assert_eq!(
        receipts
            .iter()
            .filter(|r| matches!(r.outcome, MutationOutcome::Applied(Applied::Saved { .. })))
            .count(),
        1
    );
    assert_eq!(
        receipts
            .iter()
            .filter(|r| matches!(
                r.outcome,
                MutationOutcome::Rejected(Rejection {
                    code: RejectionCode::Conflict,
                    ..
                })
            ))
            .count(),
        1
    );
    let mut db = Library::open(&path).unwrap();
    assert_eq!(db.revisions(&access, &original.id, 0, 20).unwrap().len(), 2);
    for (request, receipt) in requests.iter().zip(receipts) {
        let replay = db.apply_reliable(&access, request).unwrap();
        assert!(replay.replayed);
        assert_eq!(replay.receipt, receipt);
    }
    assert_eq!(
        db.payload(&access, &original.id, &original.head.id)
            .unwrap()
            .bytes,
        b"original"
    );
}
