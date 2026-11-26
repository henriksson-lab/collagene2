use collagene_library_sync::{self as sync, library as wire, *};
use wire::{Blob, Digest, Payload, PayloadKind, Text};
fn id(n: u8) -> MutationId {
    MutationId::new(format!("mut-{n:032x}")).unwrap()
}
fn payload(bytes: &[u8]) -> Payload {
    Payload {
        kind: PayloadKind::NativeDocument,
        bytes: Blob::new(bytes.to_vec()).unwrap(),
        archive_hash: Digest::of(bytes),
        selection: None,
        sharing: None,
    }
}
fn request(n: u8, expected_lab: wire::LabId, operation: MutationOperation) -> RequestEnvelope {
    RequestEnvelope {
        request_id: u64::from(n),
        expected_lab,
        request: Request::Apply {
            mutation: Mutation {
                id: id(n),
                operation,
            },
        },
    }
}
fn frame(raw: &[u8]) -> Vec<u8> {
    let mut b = vec![];
    b.extend_from_slice(b"CSQ\0");
    b.extend_from_slice(&1u16.to_le_bytes());
    b.extend_from_slice(&(raw.len() as u32).to_le_bytes());
    b.extend_from_slice(blake3::hash(raw).as_bytes());
    b.extend_from_slice(raw);
    b
}

#[test]
fn sync_namespace_bounds_and_malformed_frames_do_not_alias_library_v1() {
    let request = RequestEnvelope {
        request_id: 1,
        expected_lab: wire::LabId::new("lab-00000000000000000000000000000001").unwrap(),
        request: Request::Lookup {
            mutation_id: id(1),
            expected_intent: Digest::of(b"intent"),
        },
    };
    let bytes = encode_request(&request).unwrap();
    let prefix = request_prefix(&bytes).unwrap();
    assert_eq!(prefix.request_id, request.request_id);
    assert_eq!(prefix.expected_lab, request.expected_lab.as_str());
    assert!(
        (bytes.as_ptr() as usize..bytes.as_ptr() as usize + bytes.len())
            .contains(&(prefix.expected_lab.as_ptr() as usize)),
        "lab prefix must borrow the frame, not allocate an archive"
    );
    assert_eq!(&bytes[..4], b"CSQ\0");
    assert_eq!(decode_request(&bytes).unwrap(), request);
    assert!(wire::decode_request(&bytes).is_err());
    assert!(decode_response(&bytes).is_err());
    let mut corrupted = bytes.clone();
    *corrupted.last_mut().unwrap() ^= 1;
    assert_eq!(decode_request(&corrupted), Err(Error::Integrity));
    let mut future = bytes.clone();
    future[4..6].copy_from_slice(&2u16.to_le_bytes());
    assert_eq!(decode_request(&future), Err(Error::Version(2)));
    let mut trailing = bytes[HEADER_BYTES..].to_vec();
    trailing.push(0);
    assert_eq!(decode_request(&frame(&trailing)), Err(Error::Malformed));
    let mut huge = vec![0; HEADER_BYTES];
    huge[..4].copy_from_slice(b"CSQ\0");
    huge[4..6].copy_from_slice(&1u16.to_le_bytes());
    huge[6..10].copy_from_slice(&((MAX_MESSAGE_BYTES + 1) as u32).to_le_bytes());
    assert_eq!(decode_request(&huge), Err(Error::TooLarge));
    let mut huge_reply = huge;
    huge_reply[..4].copy_from_slice(b"CSR\0");
    huge_reply[6..10].copy_from_slice(&((MAX_RESPONSE_BYTES + 1) as u32).to_le_bytes());
    assert_eq!(decode_response(&huge_reply), Err(Error::TooLarge));
    #[derive(serde::Serialize)]
    #[allow(dead_code)]
    enum Raw {
        Apply {
            mutation: (),
        },
        Lookup {
            mutation_id: String,
            expected_intent: Digest,
        },
    }
    #[derive(serde::Serialize)]
    struct RawEnvelope {
        request_id: u64,
        expected_lab: wire::LabId,
        request: Raw,
    }
    let oversized_name = postcard::to_allocvec(&RawEnvelope {
        request_id: 1,
        expected_lab: wire::LabId::new("lab-00000000000000000000000000000001").unwrap(),
        request: Raw::Lookup {
            mutation_id: "x".repeat(1024 * 1024),
            expected_intent: Digest::of(b"intent"),
        },
    })
    .unwrap();
    assert_eq!(
        decode_request(&frame(&oversized_name)),
        Err(Error::Malformed)
    );
    assert!(MutationId::new("mut-../../file").is_err());
}

#[test]
fn intent_pins_operation_cas_native_bytes_and_sharing_but_not_retry_nonce() {
    let operation = MutationOperation::ImportDocument {
        folder: None,
        name: Text::new("p1").unwrap(),
        payload: payload(b"opaque native bytes"),
    };
    let a = Mutation {
        id: id(1),
        operation: operation.clone(),
    };
    let b = Mutation {
        id: id(2),
        operation: operation.clone(),
    };
    assert_eq!(a.intent_digest().unwrap(), b.intent_digest().unwrap());
    let changed = Mutation {
        id: id(1),
        operation: MutationOperation::ImportDocument {
            folder: None,
            name: Text::new("p2").unwrap(),
            payload: payload(b"opaque native bytes"),
        },
    };
    assert_ne!(a.intent_digest().unwrap(), changed.intent_digest().unwrap());
    let different_bytes = Mutation {
        id: id(1),
        operation: MutationOperation::ImportDocument {
            folder: None,
            name: Text::new("p1").unwrap(),
            payload: payload(b"different native bytes"),
        },
    };
    assert_ne!(
        a.intent_digest().unwrap(),
        different_bytes.intent_digest().unwrap()
    );
    let mut shared = payload(b"opaque archive");
    shared.kind = PayloadKind::WorkspaceArchive;
    shared.selection = Some(wire::NativeSelection {
        document_id: Text::new("native-doc").unwrap(),
        revision: Digest::of(b"molecule"),
        head_event: None,
    });
    shared.sharing = Some(wire::WorkspaceSharing::DocumentClosure);
    let x = MutationOperation::ImportDocument {
        folder: None,
        name: Text::new("p").unwrap(),
        payload: shared.clone(),
    };
    shared.sharing = Some(wire::WorkspaceSharing::ExplicitWholeWorkspace);
    let y = MutationOperation::ImportDocument {
        folder: None,
        name: Text::new("p").unwrap(),
        payload: shared,
    };
    assert_ne!(x.intent_digest().unwrap(), y.intent_digest().unwrap());
    #[cfg(not(target_arch = "wasm32"))]
    {
        for mutation in [a, b, changed, different_bytes] {
            let native = sync::native::to_storage_mutation(&mutation).unwrap();
            assert_eq!(
                mutation.intent_digest().unwrap().as_str(),
                native.intent_digest().unwrap().as_str()
            );
        }
    }
}

#[cfg(not(target_arch = "wasm32"))]
mod native_tests {
    use super::*;
    use collagene_core::{Molecule, Topology, protocol};
    use collagene_lab_library::{LabAccess, Library};
    fn native_payload(name: &str, dna: &str) -> Payload {
        let mol = Molecule::dna(name, dna, Topology::Circular).unwrap();
        payload(&protocol::save_document(&mol).unwrap())
    }
    fn setup() -> (Library, LabAccess) {
        let mut db = Library::open_in_memory().unwrap();
        let lab = db.create_lab("Shared lab").unwrap();
        let access = db.trusted_host_access(&lab.id).unwrap();
        (db, access)
    }
    fn apply(db: &mut Library, a: &LabAccess, r: &RequestEnvelope) -> ResponseEnvelope {
        let reply = sync::native::dispatch(db, a, r).unwrap();
        reply.validate_for(r).unwrap();
        decode_response_for(&encode_response(&reply).unwrap(), r).unwrap()
    }
    fn receipt(reply: ResponseEnvelope) -> Receipt {
        match reply.outcome {
            ResponseOutcome::Receipt(r) => r,
            other => panic!("expected durable receipt: {other:?}"),
        }
    }
    fn imported(r: &Receipt) -> wire::Document {
        let MutationOutcome::Applied(Applied::Imported(d)) = &r.outcome else {
            panic!("expected imported receipt")
        };
        d.clone()
    }
    fn saved(r: &Receipt) -> wire::Document {
        let MutationOutcome::Applied(Applied::Saved { document, .. }) = &r.outcome else {
            panic!("expected saved receipt")
        };
        document.clone()
    }

    #[test]
    fn lost_import_response_replays_exact_created_identity_after_other_metadata_changes() {
        let (mut db, access) = setup();
        let import = request(
            1,
            wire::LabId::new(access.lab_id().as_str()).unwrap(),
            MutationOperation::ImportDocument {
                folder: None,
                name: Text::new("pReporter").unwrap(),
                payload: native_payload("pReporter", "AGTCCGAT"),
            },
        );
        let first = receipt(apply(&mut db, &access, &import));
        let original = imported(&first);
        let native_id = collagene_lab_library::DocumentId::parse(original.id.as_str()).unwrap();
        db.rename_document(
            &access,
            &native_id,
            original.generation,
            "Renamed by another member",
        )
        .unwrap();
        let retry = receipt(apply(&mut db, &access, &import));
        assert_eq!(retry, first);
        assert_eq!(
            db.document(&access, &native_id).unwrap().name,
            "Renamed by another member"
        );
        assert_eq!(
            db.browse(&access, None, 0, 200, None)
                .unwrap()
                .entries
                .len(),
            1
        );
        let Request::Apply { mutation } = &import.request else {
            unreachable!()
        };
        let lookup = RequestEnvelope {
            request_id: 99,
            expected_lab: wire::LabId::new(access.lab_id().as_str()).unwrap(),
            request: Request::Lookup {
                mutation_id: mutation.id.clone(),
                expected_intent: mutation.intent_digest().unwrap(),
            },
        };
        assert_eq!(receipt(apply(&mut db, &access, &lookup)), first);
        let changed = Mutation {
            id: mutation.id.clone(),
            operation: MutationOperation::ImportDocument {
                folder: None,
                name: Text::new("do not create duplicate").unwrap(),
                payload: native_payload("other", "AGTCCGAT"),
            },
        };
        let mismatch = RequestEnvelope {
            request_id: 100,
            expected_lab: wire::LabId::new(access.lab_id().as_str()).unwrap(),
            request: Request::Apply { mutation: changed },
        };
        assert!(matches!(
            apply(&mut db, &access, &mismatch).outcome,
            ResponseOutcome::IntentMismatch { .. }
        ));
        assert_eq!(
            db.browse(&access, None, 0, 200, None)
                .unwrap()
                .entries
                .len(),
            1
        );
        let other = db.create_lab("Other lab").unwrap();
        let other = db.trusted_host_access(&other.id).unwrap();
        let mut other_lookup = lookup.clone();
        other_lookup.expected_lab = wire::LabId::new(other.lab_id().as_str()).unwrap();
        let absent = apply(&mut db, &other, &other_lookup);
        assert!(matches!(absent.outcome, ResponseOutcome::Unknown { .. }));
        assert!(
            absent
                .validate_for_in_context(
                    &lookup,
                    &wire::LabId::new(access.lab_id().as_str()).unwrap()
                )
                .is_err()
        );
    }

    #[test]
    fn lost_save_never_duplicates_revision_and_rejected_conflict_keeps_original_snapshot() {
        let (mut db, access) = setup();
        let import = request(
            11,
            wire::LabId::new(access.lab_id().as_str()).unwrap(),
            MutationOperation::ImportDocument {
                folder: None,
                name: Text::new("p1").unwrap(),
                payload: native_payload("A", "AGTCCGAT"),
            },
        );
        let d = imported(&receipt(apply(&mut db, &access, &import)));
        let save = request(
            12,
            wire::LabId::new(access.lab_id().as_str()).unwrap(),
            MutationOperation::SaveRevision {
                id: d.id.clone(),
                base_revision: d.head.id.clone(),
                payload: native_payload("B", "AGTCCGAA"),
            },
        );
        let first = receipt(apply(&mut db, &access, &save));
        let b = saved(&first);
        assert_eq!(receipt(apply(&mut db, &access, &save)), first);
        let native_id = collagene_lab_library::DocumentId::parse(d.id.as_str()).unwrap();
        assert_eq!(db.revisions(&access, &native_id, 0, 200).unwrap().len(), 2);
        let third = request(
            13,
            wire::LabId::new(access.lab_id().as_str()).unwrap(),
            MutationOperation::SaveRevision {
                id: d.id.clone(),
                base_revision: b.head.id,
                payload: native_payload("C", "AGTCCGAC"),
            },
        );
        let c = saved(&receipt(apply(&mut db, &access, &third)));
        let stale = request(
            14,
            wire::LabId::new(access.lab_id().as_str()).unwrap(),
            MutationOperation::SaveRevision {
                id: d.id.clone(),
                base_revision: d.head.id,
                payload: native_payload("unpublished", "AGTCCGAG"),
            },
        );
        let conflict = receipt(apply(&mut db, &access, &stale));
        let MutationOutcome::Rejected(Rejection {
            code: RejectionCode::Conflict,
            conflict: Some(MutationConflict::Document { current }),
            ..
        }) = &conflict.outcome
        else {
            panic!("expected durable conflict")
        };
        assert_eq!(current.head.id, c.head.id);
        let fourth = request(
            15,
            wire::LabId::new(access.lab_id().as_str()).unwrap(),
            MutationOperation::SaveRevision {
                id: d.id,
                base_revision: c.head.id,
                payload: native_payload("D", "AGTCCGTT"),
            },
        );
        let latest = saved(&receipt(apply(&mut db, &access, &fourth)));
        assert_eq!(receipt(apply(&mut db, &access, &stale)), conflict);
        assert_eq!(
            db.document(&access, &native_id).unwrap().head.id.as_str(),
            latest.head.id.as_str()
        );
        assert_eq!(db.revisions(&access, &native_id, 0, 200).unwrap().len(), 4);
    }

    #[test]
    fn malformed_native_payload_has_no_effect_or_false_durable_success() {
        let (mut db, access) = setup();
        let invalid = request(
            20,
            wire::LabId::new(access.lab_id().as_str()).unwrap(),
            MutationOperation::ImportDocument {
                folder: None,
                name: Text::new("bad").unwrap(),
                payload: payload(b"not a native document"),
            },
        );
        let reply = apply(&mut db, &access, &invalid);
        assert!(matches!(
            reply.outcome,
            ResponseOutcome::Failure(wire::Failure {
                code: wire::ErrorCode::Invalid,
                ..
            })
        ));
        let Request::Apply { mutation } = &invalid.request else {
            unreachable!()
        };
        let lookup = RequestEnvelope {
            request_id: 21,
            expected_lab: wire::LabId::new(access.lab_id().as_str()).unwrap(),
            request: Request::Lookup {
                mutation_id: mutation.id.clone(),
                expected_intent: mutation.intent_digest().unwrap(),
            },
        };
        assert!(matches!(
            apply(&mut db, &access, &lookup).outcome,
            ResponseOutcome::Unknown { .. }
        ));
        assert!(
            db.browse(&access, None, 0, 200, None)
                .unwrap()
                .entries
                .is_empty()
        );
    }

    #[test]
    fn expected_lab_is_a_precondition_before_apply_lookup_or_native_review() {
        let (mut db, a) = setup();
        let b_lab = db
            .create_lab("Reconfigured same-origin server lab")
            .unwrap();
        let b = db.trusted_host_access(&b_lab.id).unwrap();
        let intended = request(
            30,
            wire::LabId::new(a.lab_id().as_str()).unwrap(),
            MutationOperation::ImportDocument {
                folder: None,
                name: Text::new("pPinnedToA").unwrap(),
                payload: native_payload("Pinned native payload", "AGTCCGAT"),
            },
        );
        let committed_a = receipt(apply(&mut db, &a, &intended));
        assert!(matches!(
            sync::native::dispatch(&mut db, &b, &intended),
            Err(sync::native::DispatchError::WrongLab)
        ));
        assert_eq!(db.lab(&b).unwrap().catalog_generation, 0);
        assert!(
            db.browse(&b, None, 0, 200, None)
                .unwrap()
                .entries
                .is_empty()
        );
        let Request::Apply { mutation } = &intended.request else {
            unreachable!()
        };
        let native_id =
            collagene_lab_library::mutations::MutationId::parse(mutation.id.as_str()).unwrap();
        let intent =
            collagene_lab_library::Digest::parse(mutation.intent_digest().unwrap().as_str())
                .unwrap();
        assert!(
            db.lookup_reliable(&b, &native_id, &intent)
                .unwrap()
                .is_none()
        );
        let lookup = RequestEnvelope {
            request_id: 31,
            expected_lab: intended.expected_lab.clone(),
            request: Request::Lookup {
                mutation_id: mutation.id.clone(),
                expected_intent: mutation.intent_digest().unwrap(),
            },
        };
        assert!(matches!(
            sync::native::dispatch(&mut db, &b, &lookup),
            Err(sync::native::DispatchError::WrongLab)
        ));
        let mut invalid = intended.clone();
        let Request::Apply { mutation } = &mut invalid.request else {
            unreachable!()
        };
        let MutationOperation::ImportDocument { payload, .. } = &mut mutation.operation else {
            unreachable!()
        };
        payload.archive_hash = Digest::of(b"intentionally wrong native checksum");
        assert!(
            matches!(
                sync::native::dispatch(&mut db, &b, &invalid),
                Err(sync::native::DispatchError::WrongLab)
            ),
            "lab precondition must run before full payload/native validation"
        );
        // Explicitly pinning B is not an authority grant: the host already
        // independently granted B's access. The same nonce is scoped by lab.
        let mut intended_b = intended.clone();
        intended_b.expected_lab = wire::LabId::new(b.lab_id().as_str()).unwrap();
        let committed_b = receipt(apply(&mut db, &b, &intended_b));
        assert_ne!(imported(&committed_a).id, imported(&committed_b).id);
        assert!(
            matches!(
                sync::native::dispatch(&mut db, &b, &lookup),
                Err(sync::native::DispatchError::WrongLab)
            ),
            "an existing matching B receipt must not bypass the A pin"
        );
        let mut lookup_b = lookup;
        lookup_b.expected_lab = intended_b.expected_lab;
        assert_eq!(receipt(apply(&mut db, &b, &lookup_b)), committed_b);
        assert_eq!(db.browse(&b, None, 0, 200, None).unwrap().entries.len(), 1);
    }
}
