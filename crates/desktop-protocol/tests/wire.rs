use collagene_desktop_protocol::*;
use library::{Blob, Digest, Text};
fn op(n: u8) -> OperationId {
    OperationId::new(format!("op-{n:032x}")).unwrap()
}
fn request(command: Command) -> RequestEnvelope {
    RequestEnvelope {
        request_id: 7,
        operation_id: op(7),
        command,
    }
}
fn response(request: &RequestEnvelope, outcome: Outcome) -> ResponseEnvelope {
    ResponseEnvelope {
        request_id: request.request_id,
        operation_id: request.operation_id.clone(),
        outcome,
    }
}
fn frame(raw: &[u8]) -> Vec<u8> {
    let mut b = vec![];
    b.extend_from_slice(b"DTQ\0");
    b.extend_from_slice(&SCHEMA_VERSION.to_le_bytes());
    b.extend_from_slice(&(raw.len() as u32).to_le_bytes());
    b.extend_from_slice(blake3::hash(raw).as_bytes());
    b.extend_from_slice(raw);
    b
}
fn metadata(current: &[u8], previous: Option<&[u8]>) -> WorkingMetadata {
    WorkingMetadata {
        version: Some(Digest::of(b"independent pair version")),
        current: Some(SlotMetadata {
            archive_hash: Digest::of(current),
            bytes: current.len() as u64,
        }),
        previous: previous.map(|p| SlotMetadata {
            archive_hash: Digest::of(p),
            bytes: p.len() as u64,
        }),
    }
}

#[test]
fn desktop_frames_reject_wrong_namespace_versions_corruption_and_length_before_decode() {
    let r = request(Command::WorkingMetadata);
    let b = encode_request(&r).unwrap();
    assert_eq!(&b[..4], b"DTQ\0");
    assert_eq!(decode_request(&b).unwrap(), r);
    assert!(decode_response(&b).is_err());
    assert!(library::decode_request(&b).is_err());
    let mut bad = b.clone();
    *bad.last_mut().unwrap() ^= 1;
    assert_eq!(decode_request(&bad), Err(Error::Integrity));
    let mut newer = b.clone();
    let unsupported = SCHEMA_VERSION + 1;
    newer[4..6].copy_from_slice(&unsupported.to_le_bytes());
    assert_eq!(decode_request(&newer), Err(Error::Version(unsupported)));
    let mut trailing = b[HEADER_BYTES..].to_vec();
    trailing.push(0);
    assert_eq!(decode_request(&frame(&trailing)), Err(Error::Malformed));
    let mut oversized = vec![0; HEADER_BYTES];
    oversized[..4].copy_from_slice(b"DTQ\0");
    oversized[4..6].copy_from_slice(&SCHEMA_VERSION.to_le_bytes());
    oversized[6..10].copy_from_slice(&((MAX_REQUEST_BYTES + 1) as u32).to_le_bytes());
    assert_eq!(decode_request(&oversized), Err(Error::TooLarge));
    oversized[..4].copy_from_slice(b"DTR\0");
    oversized[6..10].copy_from_slice(&((MAX_RESPONSE_BYTES + 1) as u32).to_le_bytes());
    assert_eq!(decode_response(&oversized), Err(Error::TooLarge));
    assert!(OperationId::new("../../profile").is_err());
}

#[test]
fn caller_previous_working_archive_is_bound_independently_and_read_cas_is_exact() {
    // These are deliberately opaque byte fixtures: native host/loaders, not IPC,
    // establish valid workspace contents and scientific ancestry.
    let current = b"new current native archive";
    let previous = b"caller's independent previous archive";
    let r = request(Command::ReplaceWorking {
        expected_version: None,
        current: Blob::new(current.to_vec()).unwrap(),
        previous: Blob::new(previous.to_vec()).unwrap(),
    });
    let reply = response(
        &r,
        Outcome::Success(Success::WorkingMetadata(metadata(current, Some(previous)))),
    );
    reply.validate_for(&r).unwrap();
    assert_eq!(
        decode_response_for(&encode_response(&reply).unwrap(), &r).unwrap(),
        reply
    );
    let mut duplicated = reply.clone();
    let Outcome::Success(Success::WorkingMetadata(m)) = &mut duplicated.outcome else {
        unreachable!()
    };
    m.previous = metadata(current, None).current;
    assert_eq!(duplicated.validate_for(&r), Err(Error::Mismatch));
    let version = Digest::of(b"pinned pair");
    let read = request(Command::ReadWorking {
        slot: WorkingSlot::Previous,
        expected_version: Some(version.clone()),
    });
    let good = response(
        &read,
        Outcome::Success(Success::WorkingArchive {
            slot: WorkingSlot::Previous,
            version: version.clone(),
            archive: Blob::new(previous.to_vec()).unwrap(),
        }),
    );
    good.validate_for(&read).unwrap();
    let wrong = response(
        &read,
        Outcome::Success(Success::WorkingArchive {
            slot: WorkingSlot::Current,
            version,
            archive: Blob::new(current.to_vec()).unwrap(),
        }),
    );
    assert_eq!(wrong.validate_for(&read), Err(Error::Mismatch));
    assert!(
        request(Command::SaveCurrent {
            expected_version: None,
            archive: Blob::new(vec![]).unwrap()
        })
        .validate()
        .is_err()
    );
    assert!(
        SlotMetadata {
            archive_hash: Digest::of(b"x"),
            bytes: (MAX_ARCHIVE_BYTES + 1) as u64
        }
        .validate()
        .is_err()
    );
}

#[test]
fn working_and_preference_conflicts_cannot_cross_scopes_or_claim_unchanged_versions() {
    let old = Digest::of(b"old version");
    let bytes = Blob::new(b"small view prefs".to_vec()).unwrap();
    let r = request(Command::SavePreferences {
        expected_version: Some(old.clone()),
        bytes: bytes.clone(),
    });
    let p = Preferences {
        version: Some(Digest::of(b"new preference version")),
        bytes: Blob::new(b"other view prefs".to_vec()).unwrap(),
    };
    response(
        &r,
        Outcome::Failure(Failure::Conflict {
            current: Conflict::Preferences(p),
        }),
    )
    .validate_for(&r)
    .unwrap();
    let wrong = response(
        &r,
        Outcome::Failure(Failure::Conflict {
            current: Conflict::Working(metadata(b"workspace", None)),
        }),
    );
    assert_eq!(wrong.validate_for(&r), Err(Error::Mismatch));
    let stale = response(
        &r,
        Outcome::Failure(Failure::Conflict {
            current: Conflict::Preferences(Preferences {
                version: Some(old),
                bytes: bytes.clone(),
            }),
        }),
    );
    assert_eq!(stale.validate_for(&r), Err(Error::Mismatch));
    assert!(Blob::<MAX_PREFERENCE_BYTES>::new(vec![0; MAX_PREFERENCE_BYTES + 1]).is_err());
    let saved = response(
        &r,
        Outcome::Success(Success::Preferences(Preferences {
            version: Some(Digest::of(b"saved preferences")),
            bytes,
        })),
    );
    saved.validate_for(&r).unwrap();
}

#[test]
fn small_control_path_checks_whole_control_and_reports_commit_too_late_honestly() {
    let target = op(1);
    let cancel = request(Command::Cancel {
        operation_id: target.clone(),
    });
    let bytes = encode_request(&cancel).unwrap();
    let control = decode_control(&bytes).unwrap().unwrap();
    assert_eq!(
        control,
        ControlRequest {
            request_id: 7,
            operation_id: op(7),
            control: Control::Cancel {
                operation_id: target.clone()
            }
        }
    );
    let too_late = response(
        &cancel,
        Outcome::Success(Success::Cancelled {
            operation_id: target.clone(),
            status: OperationState::Committing,
        }),
    );
    too_late.validate_for(&cancel).unwrap();
    let falsely_accepted = response(
        &cancel,
        Outcome::Success(Success::Cancelled {
            operation_id: target,
            status: OperationState::Running,
        }),
    );
    assert!(falsely_accepted.validate().is_err());
    let mut corrupt = bytes.clone();
    *corrupt.last_mut().unwrap() ^= 1;
    assert!(decode_control(&corrupt).is_err());
    let mut oversized = bytes[HEADER_BYTES..].to_vec();
    oversized.resize(MAX_CONTROL_BYTES + 1, 0);
    assert_eq!(decode_control(&frame(&oversized)), Err(Error::TooLarge));
    let normal = request(Command::SaveCurrent {
        expected_version: None,
        archive: Blob::new(vec![1; 4096]).unwrap(),
    });
    let mut raw = encode_request(&normal).unwrap();
    assert_eq!(decode_control(&raw).unwrap(), None);
    // None is deliberately not a grant to execute a normal frame. The full
    // decoder still rejects corruption before the host can run any mutation.
    *raw.last_mut().unwrap() ^= 1;
    assert_eq!(decode_control(&raw).unwrap(), None);
    assert_eq!(decode_request(&raw), Err(Error::Integrity));
}

#[test]
fn hostile_prefixes_and_unrelated_library_conflicts_do_not_bind_to_an_operation() {
    #[derive(serde::Serialize)]
    struct Raw {
        request_id: u64,
        operation_id: String,
        command: Command,
    }
    let raw = postcard::to_allocvec(&Raw {
        request_id: 7,
        operation_id: "x".repeat(1024 * 1024),
        command: Command::Bootstrap,
    })
    .unwrap();
    assert!(decode_control(&frame(&raw)).is_err());
    assert!(decode_request(&frame(&raw)).is_err());
    let folder = |n: u8| library::FolderId::new(format!("fld-{n:032x}")).unwrap();
    let r = request(Command::Library(library::Request::RenameFolder {
        id: folder(1),
        expected_generation: 1,
        name: Text::new("new name").unwrap(),
    }));
    let wrong = response(
        &r,
        Outcome::Failure(Failure::Library(library::Failure {
            code: library::ErrorCode::Conflict,
            message: Text::new("Changed").unwrap(),
            conflict: Some(library::Conflict::Folder {
                current: library::Folder {
                    id: folder(2),
                    parent: None,
                    name: Text::new("unrelated").unwrap(),
                    generation: 2,
                },
            }),
        })),
    );
    assert_eq!(wrong.validate_for(&r), Err(Error::Mismatch));
    let mut wrong_operation = response(&r, Outcome::Failure(Failure::Cancelled));
    wrong_operation.operation_id = op(8);
    assert_eq!(wrong_operation.validate_for(&r), Err(Error::Mismatch));
}
