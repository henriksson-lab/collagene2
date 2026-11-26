use collagene_desktop_protocol::*;
use library::{Digest, LabId, Text};
fn op(n: u8) -> OperationId {
    OperationId::new(format!("op-{n:032x}")).unwrap()
}
fn request(command: Command) -> RequestEnvelope {
    RequestEnvelope {
        request_id: 17,
        operation_id: op(3),
        command,
    }
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
        schema_version: 1,
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
fn backup_frames_require_version_two_even_for_generic_failure_acknowledgments() {
    let old = request(Command::WorkingMetadata);
    let legacy = encode_request(&old).unwrap();
    assert_eq!(&legacy[4..6], &1u16.to_le_bytes());
    assert!(!catalog_backup_request_hint(&legacy));
    let q = request(Command::CatalogBackup(CatalogBackupCommand::Verify {
        id: op(1),
        expected_hash: None,
    }));
    let bytes = encode_request(&q).unwrap();
    assert_eq!(&bytes[4..6], &2u16.to_le_bytes());
    assert!(catalog_backup_request_hint(&bytes));
    assert_eq!(decode_request(&bytes).unwrap(), q);
    let mut obsolete = bytes.clone();
    obsolete[4..6].copy_from_slice(&1u16.to_le_bytes());
    assert!(decode_request(&obsolete).is_err());
    let r = ResponseEnvelope {
        request_id: q.request_id,
        operation_id: q.operation_id.clone(),
        outcome: Outcome::Failure(Failure::Unavailable),
    };
    let old_ack = encode_response(&r).unwrap();
    assert!(decode_response_for(&old_ack, &q).is_err());
    let ack = encode_response_for_frame(&r, &bytes).unwrap();
    assert_eq!(&ack[4..6], &2u16.to_le_bytes());
    assert_eq!(decode_response_for(&ack, &q).unwrap(), r);
    let old_r = ResponseEnvelope {
        request_id: old.request_id,
        operation_id: old.operation_id.clone(),
        outcome: Outcome::Failure(Failure::Unavailable),
    };
    assert_eq!(
        &encode_response_for_frame(&old_r, &legacy).unwrap()[4..6],
        &1u16.to_le_bytes()
    );
}
#[test]
fn backup_receipts_bind_selected_archive_hash_and_distinct_recovery_operation() {
    let r = receipt();
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
