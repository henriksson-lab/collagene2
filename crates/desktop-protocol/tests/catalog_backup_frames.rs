use collagene_desktop_protocol::*;
use library::{Digest, LabId, Text};
fn op(n: u8) -> OperationId {
    OperationId::new(format!("op-{n:032x}")).unwrap()
}
fn receipt() -> CatalogBackupReceipt {
    CatalogBackupReceipt {
        id: op(1),
        basename: Text::new(format!("{}.cglbak", op(1).as_str())).unwrap(),
        archive_hash: Digest::of(b"independent archive"),
        archive_bytes: 5000,
        database_hash: Digest::of(b"independent SQLite payload"),
        database_bytes: 4096,
        schema_hash: Digest::of(b"schema"),
        schema_version: library::DATABASE_SCHEMA_VERSION,
        lab: LabId::new("lab-0123456789abcdef0123456789abcdef").unwrap(),
        lab_count: 1,
        documents: 2,
        folders: 1,
        revisions: 3,
        mutation_receipts: 4,
        created_unix_seconds: 123,
        request_checkpoint: None,
    }
}
#[test]
fn backup_receipts_bind_selected_archive_hash_and_distinct_recovery_operation() {
    let r = receipt();
    r.validate().unwrap();
    for version in [0, 1, library::DATABASE_SCHEMA_VERSION + 1] {
        let mut unsupported = r.clone();
        unsupported.schema_version = version;
        assert!(unsupported.validate().is_err());
    }
    let verify = CatalogBackupCommand::Verify {
        id: op(1),
        expected_hash: Some(r.archive_hash.clone()),
    };
    CatalogBackupSuccess::Verified(r.clone())
        .validate_for(&verify, &op(3))
        .unwrap();
    let mut wrong = r.clone();
    wrong.archive_hash = Digest::of(b"changed after selection");
    assert!(
        CatalogBackupSuccess::Verified(wrong)
            .validate_for(&verify, &op(3))
            .is_err()
    );
    assert!(
        CatalogBackupSuccess::Created(r.clone())
            .validate_for(&CatalogBackupCommand::Create, &op(3))
            .is_err()
    );
    let restore = CatalogBackupCommand::RestoreCopy {
        id: op(1),
        expected_hash: r.archive_hash.clone(),
    };
    let copy = CatalogRecoveryReceipt {
        operation_id: op(3),
        backup: r,
        directory_basename: Text::new(op(3).as_str()).unwrap(),
        database_basename: Text::new("library.sqlite").unwrap(),
    };
    CatalogBackupSuccess::Restored(copy.clone())
        .validate_for(&restore, &op(3))
        .unwrap();
    assert!(
        CatalogBackupSuccess::Restored(copy.clone())
            .validate_for(&restore, &op(4))
            .is_err()
    );
    let mut wrong_path = copy;
    wrong_path.database_basename = Text::new("../library.sqlite").unwrap();
    assert!(
        CatalogBackupSuccess::Restored(wrong_path)
            .validate()
            .is_err()
    );
}
