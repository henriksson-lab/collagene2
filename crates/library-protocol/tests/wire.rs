use collagene_library_protocol::*;

fn folder(n: u8) -> FolderId {
    FolderId::new(format!("fld-{n:032x}")).unwrap()
}
fn document(n: u8) -> DocumentId {
    DocumentId::new(format!("doc-{n:032x}")).unwrap()
}
fn revision(n: u8) -> RevisionId {
    RevisionId::new(format!("{n:064x}")).unwrap()
}
fn context() -> Option<LabContext> {
    Some(LabContext {
        lab_id: LabId::new(format!("lab-{:032x}", 1)).unwrap(),
        name: Text::new("Our lab").unwrap(),
    })
}
fn payload() -> Payload {
    let bytes = b"host native-loader validation is a separate boundary".to_vec();
    Payload {
        archive_hash: Digest::of(&bytes),
        bytes: Blob::new(bytes).unwrap(),
        kind: PayloadKind::NativeDocument,
        selection: None,
        sharing: None,
    }
}
fn summary(n: u8, parent: Option<RevisionId>, p: &Payload) -> RevisionSummary {
    RevisionSummary {
        id: revision(n),
        parent,
        number: u64::from(n),
        archive_hash: p.archive_hash.clone(),
        archive_size: p.bytes.len() as u64,
        kind: p.kind,
        selection: p.selection.clone(),
        sharing: p.sharing,
        created_at: 42,
    }
}
fn doc(p: &Payload) -> Document {
    Document {
        id: document(1),
        folder: None,
        name: Text::new("pReporter").unwrap(),
        generation: 1,
        head: summary(1, None, p),
    }
}
fn repair_hash(bytes: &mut [u8]) {
    let hash = blake3::hash(&bytes[HEADER_BYTES..]);
    bytes[10..42].copy_from_slice(hash.as_bytes());
}

#[test]
fn versioned_binary_envelope_rejects_corruption_trailing_payload_and_other_direction() {
    let request = RequestEnvelope {
        request_id: 7,
        request: Request::Browse {
            parent: None,
            page: PageRequest::default(),
        },
    };
    let bytes = encode_request(&request).unwrap();
    assert_eq!(&bytes[..4], b"CLQ\0");
    assert_eq!(&bytes[4..6], &1u16.to_le_bytes());
    assert_eq!(
        u32::from_le_bytes(bytes[6..10].try_into().unwrap()) as usize,
        bytes.len() - HEADER_BYTES
    );
    assert_eq!(decode_request(&bytes).unwrap(), request);
    assert!(decode_response(&bytes).is_err());
    let mut corrupt = bytes.clone();
    *corrupt.last_mut().unwrap() ^= 1;
    assert_eq!(decode_request(&corrupt), Err(Error::Integrity));
    let mut future = bytes.clone();
    future[4..6].copy_from_slice(&2u16.to_le_bytes());
    assert_eq!(decode_request(&future), Err(Error::Version(2)));
    let mut trailing = bytes;
    trailing.push(0);
    let size = (trailing.len() - HEADER_BYTES) as u32;
    trailing[6..10].copy_from_slice(&size.to_le_bytes());
    repair_hash(&mut trailing);
    assert_eq!(decode_request(&trailing), Err(Error::Malformed));
    let mut oversized = vec![0; HEADER_BYTES];
    oversized[..4].copy_from_slice(b"CLQ\0");
    oversized[4..6].copy_from_slice(&1u16.to_le_bytes());
    oversized[6..10].copy_from_slice(&((MAX_MESSAGE_BYTES + 1) as u32).to_le_bytes());
    assert_eq!(decode_request(&oversized), Err(Error::TooLarge));
}

#[test]
fn bounded_binary_visitors_reject_lengths_before_owned_collection_allocation() {
    let text = postcard::to_allocvec(&"12345").unwrap();
    assert!(postcard::from_bytes::<Text<4>>(&text).is_err());
    let bytes = postcard::to_allocvec(&vec![1u8; 17]).unwrap();
    assert!(postcard::from_bytes::<Blob<16>>(&bytes).is_err());
    let rows = postcard::to_allocvec(&vec![1u8, 2, 3]).unwrap();
    assert!(postcard::from_bytes::<Rows<u8, 2>>(&rows).is_err());
    // 2^32-1 announced entries, no backing data. The length visitor rejects it.
    assert!(postcard::from_bytes::<Rows<u8, 2>>(&[255, 255, 255, 255, 15]).is_err());
    assert!(FolderId::new("../../outside").is_err());
    assert!(DocumentId::new("fld-00000000000000000000000000000001").is_err());
    let malicious = postcard::to_allocvec(&"fld-0000000000000000000000000000000Z").unwrap();
    assert!(postcard::from_bytes::<FolderId>(&malicious).is_err());
    let unicode = "λ".repeat(256);
    validate_name(&unicode).unwrap();
    assert!(validate_name(&format!("{unicode}λ")).is_err());
    assert!(validate_name("\0hidden").is_err());
}

#[test]
fn save_reply_binds_document_parent_payload_and_host_context_without_accepting_identity_claims() {
    let p = payload();
    let request = RequestEnvelope {
        request_id: 4,
        request: Request::SaveRevision {
            id: document(1),
            base_revision: revision(1),
            payload: p.clone(),
        },
    };
    let mut d = doc(&p);
    let r = summary(2, Some(revision(1)), &p);
    d.head = r.clone();
    d.generation = 2;
    let reply = ResponseEnvelope {
        request_id: 4,
        context: context(),
        outcome: Outcome::Success(Success::Saved {
            document: d.clone(),
            revision: r.clone(),
        }),
    };
    reply.validate_for(&request).unwrap();
    let bytes = encode_response(&reply).unwrap();
    assert_eq!(decode_response_for(&bytes, &request).unwrap(), reply);
    let wrong_lab = LabId::new(format!("lab-{:032x}", 2)).unwrap();
    assert_eq!(
        reply.validate_for_in_context(&request, &wrong_lab),
        Err(Error::ResponseMismatch)
    );
    let mut unrelated = reply.clone();
    {
        let Outcome::Success(Success::Saved {
            revision: child,
            document: changed,
        }) = &mut unrelated.outcome
        else {
            unreachable!()
        };
        child.parent = Some(revision(3));
        changed.head = child.clone();
    }
    assert_eq!(
        unrelated.validate_for(&request),
        Err(Error::ResponseMismatch)
    );
    let mut unscoped = reply;
    unscoped.context = None;
    assert!(unscoped.validate().is_err());
    let conflict = ResponseEnvelope {
        request_id: 4,
        context: context(),
        outcome: Outcome::Failure(Failure {
            code: ErrorCode::Conflict,
            message: Text::new("The saved head changed; reopen or explicitly save a copy.")
                .unwrap(),
            conflict: Some(Conflict::Document { current: d.clone() }),
        }),
    };
    // Conflict is returned explicitly; it is never converted into successful overwrite.
    let _ = r;
    assert!(conflict.validate_for(&request).is_ok());
}

#[test]
fn snapshots_scope_and_folder_membership_are_bound_on_page_replies() {
    let request = RequestEnvelope {
        request_id: 9,
        request: Request::Browse {
            parent: Some(folder(1)),
            page: PageRequest {
                offset: 0,
                limit: 2,
                expected_catalog_generation: Some(3),
            },
        },
    };
    let f = Folder {
        id: folder(2),
        parent: Some(folder(1)),
        name: Text::new("Vectors").unwrap(),
        generation: 1,
    };
    let reply = ResponseEnvelope {
        request_id: 9,
        context: context(),
        outcome: Outcome::Success(Success::Browsed(Page {
            catalog_generation: 3,
            entries: Rows::new(vec![Entry::Folder(f)]).unwrap(),
            next_offset: None,
        })),
    };
    reply.validate_for(&request).unwrap();
    let mut stale = reply.clone();
    let Outcome::Success(Success::Browsed(page)) = &mut stale.outcome else {
        unreachable!()
    };
    page.catalog_generation = 4;
    assert_eq!(stale.validate_for(&request), Err(Error::ResponseMismatch));
    let mut wrong_parent = reply;
    let Outcome::Success(Success::Browsed(page)) = &mut wrong_parent.outcome else {
        unreachable!()
    };
    let mut entries = page.entries.clone().into_vec();
    let Entry::Folder(f) = &mut entries[0] else {
        unreachable!()
    };
    f.parent = None;
    page.entries = Rows::new(entries).unwrap();
    assert_eq!(
        wrong_parent.validate_for(&request),
        Err(Error::ResponseMismatch)
    );
    assert!(
        PageRequest {
            offset: 1,
            limit: 20,
            expected_catalog_generation: None
        }
        .validate()
        .is_err()
    );
    Page {
        catalog_generation: 0,
        entries: Rows::new(vec![]).unwrap(),
        next_offset: None,
    }
    .validate()
    .unwrap();
}

#[test]
fn workspace_payload_requires_explicit_sharing_selection_and_exact_archive_integrity() {
    let mut p = payload();
    p.kind = PayloadKind::WorkspaceArchive;
    assert!(p.validate().is_err());
    p.sharing = Some(WorkspaceSharing::DocumentClosure);
    p.selection = Some(NativeSelection {
        document_id: Text::new("local-doc").unwrap(),
        revision: Digest::of(b"selected sequence model"),
        head_event: Some(Text::new("event-selected").unwrap()),
    });
    p.validate().unwrap();
    let r = summary(1, None, &p);
    let mut d = doc(&p);
    d.head = r.clone();
    let reply = ResponseEnvelope {
        request_id: 3,
        context: context(),
        outcome: Outcome::Success(Success::Opened {
            document: d,
            revision: r,
            payload: p.clone(),
        }),
    };
    reply.validate().unwrap();
    let mut wrong = reply;
    let Outcome::Success(Success::Opened {
        payload: opened_payload,
        ..
    }) = &mut wrong.outcome
    else {
        unreachable!()
    };
    opened_payload.sharing = Some(WorkspaceSharing::ExplicitWholeWorkspace);
    assert!(wrong.validate().is_err());
    let mut bad_hash = p.clone();
    bad_hash.archive_hash = Digest::of(b"different archive");
    assert_eq!(bad_hash.validate(), Err(Error::Integrity));
    let mut overlarge = summary(1, None, &payload());
    overlarge.archive_size = (MAX_DOCUMENT_BYTES + 1) as u64;
    assert_eq!(overlarge.validate(), Err(Error::TooLarge));
    // These arbitrary bytes pass only transport integrity. A host must reject
    // them through the native loader before admitting anything to the database.
    assert!(!p.bytes.starts_with(b"PK"));
}

#[test]
fn mutation_variants_are_bounded_and_errors_cannot_pretend_to_succeed() {
    let requests = vec![
        Request::CreateFolder {
            parent: None,
            name: Text::new("Plasmids").unwrap(),
        },
        Request::RenameFolder {
            id: folder(1),
            expected_generation: 1,
            name: Text::new("Vectors").unwrap(),
        },
        Request::MoveFolder {
            id: folder(1),
            expected_generation: 1,
            parent: Some(folder(2)),
        },
        Request::DeleteEmptyFolder {
            id: folder(1),
            expected_generation: 1,
        },
        Request::ImportDocument {
            folder: None,
            name: Text::new("p1").unwrap(),
            payload: payload(),
        },
        Request::OpenDocument {
            id: document(1),
            revision: Some(revision(1)),
        },
        Request::RenameDocument {
            id: document(1),
            expected_generation: 1,
            name: Text::new("p2").unwrap(),
        },
        Request::MoveDocument {
            id: document(1),
            expected_generation: 1,
            folder: Some(folder(1)),
        },
        Request::Search {
            query: Text::new("reporter").unwrap(),
            scope: SearchScope::AllFolders,
            page: PageRequest::default(),
        },
        Request::ListRevisions {
            id: document(1),
            page: PageRequest::default(),
        },
    ];
    for (i, r) in requests.into_iter().enumerate() {
        let r = RequestEnvelope {
            request_id: (i + 1) as u64,
            request: r,
        };
        assert_eq!(decode_request(&encode_request(&r).unwrap()).unwrap(), r);
    }
    assert!(
        Request::MoveFolder {
            id: folder(1),
            expected_generation: 1,
            parent: Some(folder(1))
        }
        .validate()
        .is_err()
    );
    let unauthed = ResponseEnvelope {
        request_id: 1,
        context: None,
        outcome: Outcome::Failure(Failure {
            code: ErrorCode::Unauthenticated,
            message: Text::new("Sign in to the lab server.").unwrap(),
            conflict: None,
        }),
    };
    unauthed.validate().unwrap();
    let malformed = ResponseEnvelope {
        request_id: 1,
        context: context(),
        outcome: Outcome::Failure(Failure {
            code: ErrorCode::Conflict,
            message: Text::new("changed").unwrap(),
            conflict: None,
        }),
    };
    assert!(malformed.validate().is_err());
}
