use collagene_desktop_sync_protocol as wire;
use wire::library;
fn remote() -> wire::RemoteId {
    wire::RemoteId::new(format!("rem-{:032x}", 1)).unwrap()
}
fn id() -> wire::ReviewId {
    wire::ReviewId::new(format!("review-{:032x}", 1)).unwrap()
}
fn description() -> wire::RemovalDescription {
    wire::RemovalDescription {
        local_lab: library::LabId::new(format!("lab-{:032x}", 2)).unwrap(),
        remote: wire::PublicationRemote {
            id: remote(),
            name: library::Text::new("Lab").unwrap(),
            endpoint: library::Text::new("https://lab.example").unwrap(),
            lab: library::LabId::new(format!("lab-{:032x}", 3)).unwrap(),
        },
        expected_checkpoint: Some(library::Digest::of(b"checkpoint")),
        document_bindings: 2,
        folder_bindings: 1,
        completed_tree_jobs: 1,
    }
}
fn request(command: wire::Command) -> wire::RequestEnvelope {
    wire::RequestEnvelope {
        request_id: 1,
        operation_id: wire::OperationId::new(format!("op-{:032x}", 1)).unwrap(),
        command,
    }
}
fn response(q: &wire::RequestEnvelope, success: wire::Success) -> wire::ResponseEnvelope {
    wire::ResponseEnvelope {
        request_id: q.request_id,
        operation_id: q.operation_id.clone(),
        context: None,
        outcome: wire::Outcome::Success(success),
    }
}
#[test]
fn current_removal_and_metadata_commands_reject_noncurrent_headers() {
    let commands = [
        wire::Command::PrepareRemoteRemoval { remote: remote() },
        wire::Command::ConfirmRemoteRemoval {
            remote: remote(),
            review: id(),
        },
        wire::Command::DiscardRemoteRemovalReview { review: id() },
        wire::Command::CheckRemoteRemovalOutcome { remote: remote() },
    ];
    for c in commands {
        let q = request(c);
        let bytes = wire::encode_request(&q).unwrap();
        assert_eq!(&bytes[4..6], &3u16.to_le_bytes());
        assert_eq!(wire::decode_request(&bytes).unwrap(), q);
        for old in [1u16, 2] {
            let mut bad = bytes.clone();
            bad[4..6].copy_from_slice(&old.to_le_bytes());
            assert!(wire::decode_request(&bad).is_err());
        }
    }
    let ordinary = wire::encode_request(&request(wire::Command::Remotes)).unwrap();
    assert_eq!(&ordinary[4..6], &wire::SCHEMA_VERSION.to_le_bytes());
    let prior = wire::encode_request(&request(wire::Command::PrepareMetadataPublication {
        target: wire::PublicationTarget::Folder(
            library::FolderId::new(format!("fld-{:032x}", 1)).unwrap(),
        ),
        remote: remote(),
    }))
    .unwrap();
    assert_eq!(&prior[4..6], &wire::SCHEMA_VERSION.to_le_bytes());
}
#[test]
fn removal_reply_is_exactly_request_bound_and_unknown_receipt_versions_or_counts_fail() {
    let q = request(wire::Command::PrepareRemoteRemoval { remote: remote() });
    let review = wire::RemoteRemovalReview {
        id: id(),
        description: description(),
        digest: library::Digest::of(b"review"),
    };
    let r = response(&q, wire::Success::RemoteRemovalPrepared { review });
    let bytes = wire::encode_response_for_frame(&r, &wire::encode_request(&q).unwrap()).unwrap();
    wire::decode_response_for(&bytes, &q).unwrap();
    for old in [1u16, 2] {
        let mut b = bytes.clone();
        b[4..6].copy_from_slice(&old.to_le_bytes());
        assert!(wire::decode_response(&b).is_err());
    }
    let mut wrong = q.clone();
    wrong.command = wire::Command::PrepareRemoteRemoval {
        remote: wire::RemoteId::new(format!("rem-{:032x}", 9)).unwrap(),
    };
    assert!(r.validate_for(&wrong).is_err());
    let d = description();
    let receipt = wire::RemovalReceipt {
        schema_version: 1,
        local_lab: d.local_lab,
        remote: d.remote,
        review_digest: library::Digest::of(b"review"),
        removed_document_bindings: 2,
        removed_folder_bindings: 1,
        removed_completed_tree_jobs: 1,
    };
    let mut bad = receipt.clone();
    bad.schema_version = 2;
    assert!(bad.validate().is_err());
    bad = receipt.clone();
    bad.removed_completed_tree_jobs = 129;
    assert!(bad.validate().is_err());
    let q = request(wire::Command::CheckRemoteRemovalOutcome { remote: remote() });
    let r = response(
        &q,
        wire::Success::RemoteRemovalOutcome {
            remote: remote(),
            receipt: Some(receipt),
        },
    );
    wire::decode_response_for(&wire::encode_response(&r).unwrap(), &q).unwrap();
}
