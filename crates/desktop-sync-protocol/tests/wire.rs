use collagene_desktop_sync_protocol::{library, *};

fn remote() -> Remote {
    Remote {
        id: RemoteId::new(format!("rem-{}", "1".repeat(32))).unwrap(),
        name: library::Text::new("Lab library").unwrap(),
        endpoint: library::Text::new("https://lab.example").unwrap(),
        expected_lab: library::LabId::new(format!("lab-{}", "2".repeat(32))).unwrap(),
    }
}
fn request(command: Command) -> RequestEnvelope {
    RequestEnvelope {
        request_id: 17,
        operation_id: OperationId::new(format!("op-{}", "3".repeat(32))).unwrap(),
        command,
    }
}
fn context() -> RemoteContext {
    RemoteContext {
        remote: remote(),
        lab: library::LabContext {
            lab_id: remote().expected_lab,
            name: library::Text::new("Lab").unwrap(),
        },
    }
}
fn answer(request: &RequestEnvelope, success: Success) -> ResponseEnvelope {
    ResponseEnvelope {
        request_id: request.request_id,
        operation_id: request.operation_id.clone(),
        context: Some(context()),
        outcome: Outcome::Success(success),
    }
}
fn checked_frame(magic: &[u8; 4], raw: &[u8]) -> Vec<u8> {
    let mut out = magic.to_vec();
    out.extend_from_slice(&SCHEMA_VERSION.to_le_bytes());
    out.extend_from_slice(&(raw.len() as u32).to_le_bytes());
    out.extend_from_slice(blake3::hash(raw).as_bytes());
    out.extend_from_slice(raw);
    out
}

#[test]
fn local_bookmarks_bind_checkpoint_and_page_without_remote_authority() {
    assert!(
        request(Command::ListBookmarks {
            page: LocalPageRequest {
                offset: 1,
                limit: 1
            },
            expected_checkpoint: None
        })
        .validate()
        .is_err()
    );
    let checkpoint = library::Digest::of(b"repository checkpoint");
    let request = request(Command::ListBookmarks {
        page: LocalPageRequest {
            offset: 0,
            limit: 1,
        },
        expected_checkpoint: Some(checkpoint.clone()),
    });
    let mut reply = answer(
        &request,
        Success::Bookmarks(BookmarkPage {
            items: library::Rows::new(vec![Bookmark::NativeBranch {
                document: library::Text::new("document-1").unwrap(),
                revision: library::Digest::of(b"old sequence"),
                event: Some(library::Text::new("event-1").unwrap()),
                undo_count: 2,
                redo_count: 1,
                incoming_revision: library::Digest::of(b"new sequence"),
                incoming_event: None,
            }])
            .unwrap(),
            next_offset: Some(1),
            checkpoint: Some(checkpoint),
        }),
    );
    assert!(
        reply.validate_for(&request).is_err(),
        "a local bookmark reply cannot claim a remote context"
    );
    reply.context = None;
    assert_eq!(
        decode_response_for(&encode_response(&reply).unwrap(), &request).unwrap(),
        reply
    );
    let Outcome::Success(Success::Bookmarks(page)) = &mut reply.outcome else {
        unreachable!()
    };
    page.checkpoint = Some(library::Digest::of(b"different checkpoint"));
    assert!(reply.validate_for(&request).is_err());
    let Outcome::Success(Success::Bookmarks(page)) = &mut reply.outcome else {
        unreachable!()
    };
    page.checkpoint = match &request.command {
        Command::ListBookmarks {
            expected_checkpoint,
            ..
        } => expected_checkpoint.clone(),
        _ => unreachable!(),
    };
    page.next_offset = Some(2);
    assert!(
        reply.validate_for(&request).is_err(),
        "continuation must start after exactly the returned rows"
    );
}

#[test]
fn secret_is_bounded_on_decode_and_redacted_in_enclosing_debug() {
    let token = "private-native-access-token";
    let req = request(Command::SetCredential {
        remote: remote().id,
        secret: Secret::new(token).unwrap(),
        retention: CredentialRetention::Session,
    });
    let debug = format!("{req:?}");
    assert!(!debug.contains(token));
    assert!(debug.contains("REDACTED"));
    let decoded = decode_request(&encode_request(&req).unwrap()).unwrap();
    let Command::SetCredential { secret, .. } = decoded.command else {
        panic!()
    };
    assert_eq!(secret.expose_secret(), token);
    assert!(Secret::new("").is_err());
    assert!(Secret::new("bad\r\nAuthorization").is_err());
    assert!(Secret::new("a".repeat(MAX_SECRET_BYTES + 1)).is_err());
    // Text is wire-compatible with Secret, but may have been constructed by an
    // untrusted caller using a larger bound. The receiver must reject it.
    let raw = postcard::to_allocvec(
        &library::Text::<4096>::new("x".repeat(MAX_SECRET_BYTES + 1)).unwrap(),
    )
    .unwrap();
    assert!(postcard::from_bytes::<Secret>(&raw).is_err());
    let reply = answer(
        &req,
        Success::CredentialSet {
            remote: remote().id,
            retention: CredentialRetention::Session,
        },
    );
    let response = encode_response(&reply).unwrap();
    assert!(!response.windows(token.len()).any(|w| w == token.as_bytes()));
    decode_response_for(&response, &req).unwrap();
}

#[test]
fn checked_frames_reject_damage_trailing_data_wrong_namespace_and_bounds() {
    let req = request(Command::Remotes);
    let frame = encode_request(&req).unwrap();
    assert_eq!(decode_request(&frame).unwrap(), req);
    assert!(frame.len() <= MAX_REQUEST_BYTES);
    let mut damaged = frame.clone();
    *damaged.last_mut().unwrap() ^= 1;
    assert_eq!(decode_request(&damaged), Err(Error::Integrity));
    let mut wrong = frame.clone();
    wrong[1] = b'T';
    assert_eq!(decode_request(&wrong), Err(Error::Malformed));
    let mut raw = postcard::to_allocvec(&req).unwrap();
    raw.push(0);
    assert_eq!(
        decode_request(&checked_frame(b"DSQ\0", &raw)),
        Err(Error::Malformed)
    );
    let mut version = frame.clone();
    let unsupported = SCHEMA_VERSION + 1;
    version[4..6].copy_from_slice(&unsupported.to_le_bytes());
    assert_eq!(decode_request(&version), Err(Error::Version(unsupported)));
    let oversized = checked_frame(b"DSQ\0", &vec![0; MAX_REQUEST_BYTES - HEADER_BYTES + 1]);
    assert_eq!(oversized.len(), MAX_REQUEST_BYTES + 1);
    assert_eq!(decode_request(&oversized), Err(Error::TooLarge));
    let reply = answer(
        &req,
        Success::Remotes(library::Rows::new(vec![remote()]).unwrap()),
    );
    assert!(decode_request(&encode_response(&reply).unwrap()).is_err());
}

#[test]
fn browse_replies_bind_lab_locator_parent_page_and_operation() {
    let parent = library::FolderId::new(format!("fld-{}", "4".repeat(32))).unwrap();
    let req = request(Command::Browse {
        remote: remote().id,
        parent: Some(parent.clone()),
        page: library::PageRequest {
            offset: 0,
            limit: 1,
            expected_catalog_generation: Some(7),
        },
    });
    let child = library::Folder {
        id: library::FolderId::new(format!("fld-{}", "5".repeat(32))).unwrap(),
        parent: Some(parent),
        name: library::Text::new("Constructs").unwrap(),
        generation: 1,
    };
    let good = answer(
        &req,
        Success::Browsed(library::Page {
            catalog_generation: 7,
            entries: library::Rows::new(vec![library::Entry::Folder(child)]).unwrap(),
            next_offset: Some(1),
        }),
    );
    good.validate_in_remote(&req, &remote()).unwrap();
    decode_response_for(&encode_response(&good).unwrap(), &req).unwrap();
    let mut bad = good.clone();
    bad.context = None;
    assert!(bad.validate_for(&req).is_err());
    let mut bad = good.clone();
    bad.operation_id = OperationId::new(format!("op-{}", "6".repeat(32))).unwrap();
    assert!(bad.validate_for(&req).is_err());
    let mut bad = good.clone();
    bad.context.as_mut().unwrap().remote.endpoint =
        library::Text::new("https://other.example").unwrap();
    assert!(bad.validate_in_remote(&req, &remote()).is_err());
    let mut bad = good.clone();
    bad.context.as_mut().unwrap().lab.lab_id =
        library::LabId::new(format!("lab-{}", "7".repeat(32))).unwrap();
    assert!(bad.validate_for(&req).is_err());
    let mut bad = good.clone();
    let Outcome::Success(Success::Browsed(p)) = &mut bad.outcome else {
        panic!()
    };
    p.catalog_generation = 8;
    assert!(bad.validate_for(&req).is_err());
    let mut bad = good;
    let Outcome::Success(Success::Browsed(p)) = &mut bad.outcome else {
        panic!()
    };
    p.next_offset = Some(2);
    assert!(bad.validate_for(&req).is_err());
}

#[test]
fn queue_receipt_pins_saved_local_revision_and_does_not_claim_attempted_work() {
    let document = library::DocumentId::new(format!("doc-{}", "8".repeat(32))).unwrap();
    let revision = library::RevisionId::new("a".repeat(64)).unwrap();
    let req = request(Command::QueuePush {
        local_document: document.clone(),
        expected_local_revision: revision.clone(),
        remote: remote().id,
    });
    let pending = PendingPush {
        job: PushJobId::new(format!("mut-{}", "9".repeat(32))).unwrap(),
        remote: remote().id,
        target: LocalTarget::Document(document),
        local_revision: Some(revision),
        attempted: false,
        rejection: None,
    };
    let good = answer(&req, Success::Queued(pending.clone()));
    good.validate_for(&req).unwrap();
    let mut stale = pending.clone();
    stale.local_revision = None;
    assert!(
        answer(&req, Success::Queued(stale))
            .validate_for(&req)
            .is_err()
    );
    let mut already_attempted = pending;
    already_attempted.attempted = true;
    assert!(
        answer(&req, Success::Queued(already_attempted))
            .validate_for(&req)
            .is_err()
    );
}

#[test]
fn tree_controls_preserve_mount_scope_and_never_claim_partial_clone_complete() {
    let destination = library::FolderId::new(format!("fld-{}", "a".repeat(32))).unwrap();
    let req = request(Command::BeginTreeClone {
        remote: remote().id,
        root: TreeRoot::LabRoot,
        destination: Some(destination.clone()),
        limits: TreeLimits::default(),
    });
    let status = TreeStatus {
        job: TreeJobId::new(format!("tree-{}", "b".repeat(32))).unwrap(),
        remote: remote().id,
        root: TreeRoot::LabRoot,
        destination: Some(destination),
        phase: TreePhase::Scanning,
        catalog_generation: None,
        scanned_nodes: 0,
        folders_created: 0,
        documents_created: 0,
        issue: None,
    };
    answer(&req, Success::TreeStarted(status.clone()))
        .validate_for(&req)
        .unwrap();
    let mut wrong_mount = status.clone();
    wrong_mount.destination = None;
    assert!(
        answer(&req, Success::TreeStarted(wrong_mount))
            .validate_for(&req)
            .is_err()
    );
    let mut incomplete = status;
    incomplete.phase = TreePhase::Complete;
    incomplete.scanned_nodes = 2;
    incomplete.documents_created = 1;
    assert!(incomplete.validate().is_err());
    let too_many = request(Command::ResumeTreeClone {
        job: incomplete.job,
        max_steps: 1001,
    });
    assert!(encode_request(&too_many).is_err());
    let bad_page = request(Command::Pending {
        page: LocalPageRequest {
            offset: 0,
            limit: 201,
        },
    });
    assert!(encode_request(&bad_page).is_err());
    assert!(RemoteId::new(format!("rem-{}", "A".repeat(32))).is_err());
    let unsafe_remote = request(Command::AddRemote {
        name: library::Text::new("Lab").unwrap(),
        endpoint: library::Text::new("https://token@lab.example/private").unwrap(),
        expected_lab: remote().expected_lab,
    });
    assert!(unsafe_remote.validate().is_err());
}

#[test]
fn retained_conflict_receipt_summary_validates_without_fabricating_live_lab_context() {
    let job = PushJobId::new(format!("mut-{}", "a".repeat(32))).unwrap();
    let request = request(Command::ResumePush { job: job.clone() });
    let current = library::Folder {
        id: library::FolderId::new(format!("fld-{}", "b".repeat(32))).unwrap(),
        parent: None,
        name: library::Text::new("New remote name").unwrap(),
        generation: 2,
    };
    let rejection = Rejection {
        code: library::ErrorCode::Conflict,
        message: library::Text::new("The item changed; review the retained current state.")
            .unwrap(),
        conflict: Some(library::Conflict::Folder { current }),
    };
    let mut reply = answer(
        &request,
        Success::Resumed(PushReceipt {
            job,
            intent: library::Digest::of(b"original durable mutation"),
            catalog_generation: 7,
            outcome: PushOutcome::Rejected(rejection),
        }),
    );
    reply.context = None; // Historical summary, not a live remote read reply.
    decode_response_for(&encode_response(&reply).unwrap(), &request).unwrap();
    let Outcome::Success(Success::Resumed(receipt)) = &mut reply.outcome else {
        panic!()
    };
    let PushOutcome::Rejected(rejection) = &mut receipt.outcome else {
        panic!()
    };
    rejection.code = library::ErrorCode::NotFound;
    assert!(reply.validate_for(&request).is_err());
    let Outcome::Success(Success::Resumed(receipt)) = &mut reply.outcome else {
        panic!()
    };
    let PushOutcome::Rejected(rejection) = &mut receipt.outcome else {
        panic!()
    };
    rejection.code = library::ErrorCode::Conflict;
    let Some(library::Conflict::Folder { current }) = &mut rejection.conflict else {
        panic!()
    };
    current.generation = 0;
    assert!(reply.validate_for(&request).is_err());
}

#[test]
fn tree_parent_review_and_restart_bind_exact_scope_without_remote_authority() {
    let folder = library::Folder {
        id: library::FolderId::new(format!("fld-{}", "a".repeat(32))).unwrap(),
        parent: None,
        name: library::Name::new("Reviewed local parent").unwrap(),
        generation: 7,
    };
    let tree = TreeStatus {
        job: TreeJobId::new(format!("tree-{}", "b".repeat(32))).unwrap(),
        remote: remote().id,
        root: TreeRoot::LabRoot,
        destination: Some(folder.id.clone()),
        phase: TreePhase::Blocked,
        catalog_generation: Some(3),
        scanned_nodes: 2,
        folders_created: 1,
        documents_created: 0,
        issue: Some(TreeIssue::LocalParentChanged {
            folder: folder.id.clone(),
        }),
    };
    let review = TreeParentReview {
        tree: tree.clone(),
        folder,
    };
    review.validate().unwrap();
    let inspect = request(Command::InspectTreeParentChange {
        job: tree.job.clone(),
    });
    let mut reply = answer(&inspect, Success::TreeParentReview(review.clone()));
    assert!(
        reply.validate_for(&inspect).is_err(),
        "local parent metadata cannot claim remote authority"
    );
    reply.context = None;
    assert_eq!(
        decode_response_for(&encode_response(&reply).unwrap(), &inspect).unwrap(),
        reply
    );
    let accept = request(Command::AcceptTreeParentChange {
        review: review.clone(),
    });
    let mut acknowledged = tree.clone();
    acknowledged.issue = None;
    let mut reply = answer(&accept, Success::TreeParentAccepted(acknowledged.clone()));
    reply.context = None;
    assert_eq!(
        decode_response_for(&encode_response(&reply).unwrap(), &accept).unwrap(),
        reply
    );
    acknowledged.folders_created = 2;
    reply.outcome = Outcome::Success(Success::TreeParentAccepted(acknowledged));
    assert!(
        reply.validate_for(&accept).is_err(),
        "acknowledgment cannot import extra nodes"
    );
    let mut forged = review;
    forged.folder.id = library::FolderId::new(format!("fld-{}", "c".repeat(32))).unwrap();
    assert!(
        request(Command::AcceptTreeParentChange { review: forged })
            .validate()
            .is_err()
    );
    let restart = request(Command::RestartTreeClone {
        job: tree.job.clone(),
    });
    let mut restarted = tree;
    restarted.phase = TreePhase::Scanning;
    restarted.issue = None;
    restarted.catalog_generation = None;
    restarted.scanned_nodes = 0;
    restarted.folders_created = 0;
    let mut reply = answer(&restart, Success::TreeUpdated(restarted.clone()));
    reply.context = None;
    reply.validate_for(&restart).unwrap();
    restarted.folders_created = 1;
    restarted.scanned_nodes = 1;
    reply.outcome = Outcome::Success(Success::TreeUpdated(restarted));
    assert!(
        reply.validate_for(&restart).is_err(),
        "a restart receipt cannot hide applied nodes"
    );
    let mut previous = encode_request(&restart).unwrap();
    previous[4..6].copy_from_slice(&(SCHEMA_VERSION - 1).to_le_bytes());
    assert_eq!(
        decode_request(&previous),
        Err(Error::Version(SCHEMA_VERSION - 1))
    );
}
