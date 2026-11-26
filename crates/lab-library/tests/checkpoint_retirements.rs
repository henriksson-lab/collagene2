use collagene_lab_library::{
    Library,
    mutations::{CheckpointKey, CheckpointValue, CheckpointWrite, MutationError},
};
fn value(n: u8) -> CheckpointValue {
    CheckpointValue {
        version: 1,
        bytes: vec![n],
    }
}
#[test]
fn all_131_retirements_are_atomic_and_late_cas_failure_rolls_back_every_write() {
    let mut db = Library::open_in_memory().unwrap();
    let access = db
        .local_repository_access("retirement-test", "Private native")
        .unwrap();
    for batch in 0..16 {
        let rows: Vec<_> = (0..8)
            .map(|i| CheckpointWrite {
                key: CheckpointKey::repository(format!("job-{}", batch * 8 + i)),
                expected: None,
                value: Some(value(i as u8)),
            })
            .collect();
        db.commit_checkpoints(&access, &rows).unwrap();
    }
    let mut writes = Vec::new();
    for n in 0..128 {
        let key = CheckpointKey::repository(format!("job-{n}"));
        let c = db.read_checkpoint(&access, &key).unwrap().unwrap();
        writes.push(CheckpointWrite {
            key,
            expected: Some(c.digest),
            value: None,
        });
    }
    for key in ["state", "tree-jobs", "remote-removal-test"] {
        writes.push(CheckpointWrite {
            key: CheckpointKey::repository(key),
            expected: None,
            value: Some(value(9)),
        });
    }
    assert_eq!(writes.len(), 131);
    assert!(
        db.commit_checkpoints(&access, &writes).is_err(),
        "ordinary8-writeAPI must remain bounded"
    );
    let mut wrong = writes.clone();
    wrong[127].expected = None;
    assert!(matches!(
        db.commit_checkpoint_retirements(&access, &wrong),
        Err(MutationError::CheckpointConflict { .. })
    ));
    for row in &writes[..128] {
        assert!(db.read_checkpoint(&access, &row.key).unwrap().is_some());
    }
    for row in &writes[128..] {
        assert!(db.read_checkpoint(&access, &row.key).unwrap().is_none());
    }
    db.commit_checkpoint_retirements(&access, &writes).unwrap();
    for row in &writes[..128] {
        assert!(db.read_checkpoint(&access, &row.key).unwrap().is_none());
    }
    for row in &writes[128..] {
        assert_eq!(
            db.read_checkpoint(&access, &row.key)
                .unwrap()
                .unwrap()
                .value,
            value(9)
        );
    }
}
#[test]
fn retirement_bounds_and_lab_scope_never_delete_foreign_metadata() {
    let mut db = Library::open_in_memory().unwrap();
    let first = db.local_repository_access("one", "One").unwrap();
    let second = db.local_repository_access("two", "Two").unwrap();
    let key = CheckpointKey::repository("preserve");
    db.commit_checkpoints(
        &first,
        &[CheckpointWrite {
            key: key.clone(),
            expected: None,
            value: Some(value(1)),
        }],
    )
    .unwrap();
    let pin = db.read_checkpoint(&first, &key).unwrap().unwrap();
    assert!(
        db.commit_checkpoint_retirements(
            &second,
            &[CheckpointWrite {
                key: key.clone(),
                expected: Some(pin.digest.clone()),
                value: None
            }]
        )
        .is_err()
    );
    assert_eq!(db.read_checkpoint(&first, &key).unwrap().unwrap(), pin);
    let too_many: Vec<_> = (0..129)
        .map(|n| CheckpointWrite {
            key: CheckpointKey::repository(format!("absent-{n}")),
            expected: None,
            value: None,
        })
        .collect();
    assert!(db.commit_checkpoint_retirements(&first, &too_many).is_err());
    let replacements: Vec<_> = (0..4)
        .map(|n| CheckpointWrite {
            key: CheckpointKey::repository(format!("new-{n}")),
            expected: None,
            value: Some(value(2)),
        })
        .collect();
    assert!(
        db.commit_checkpoint_retirements(&first, &replacements)
            .is_err()
    );
    assert_eq!(db.read_checkpoint(&first, &key).unwrap().unwrap(), pin);
}
