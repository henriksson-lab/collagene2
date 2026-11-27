//! Native binary SDK boundary only; no HTTP listener or application probe.
use collagene_lab_library::{self as store, collaboration as native};
use collagene_library_protocol::{self as wire, discussion as api};
fn member(
    db: &store::Library,
    lab: &store::LabAccess,
    id: &str,
    role: native::Role,
) -> native::MemberAccess {
    db.trusted_member_access(
        lab,
        native::Author::new(native::ActorId::parse(id).unwrap(), id).unwrap(),
        role,
    )
    .unwrap()
}
fn call(
    db: &mut store::Library,
    m: &native::MemberAccess,
    request: api::Request,
) -> api::ResponseEnvelope {
    let request = api::RequestEnvelope {
        request_id: 17,
        request,
    };
    let frame = api::encode_request(&request).unwrap();
    assert!(frame.starts_with(b"CDQ\0\x02\x00"));
    assert!(wire::decode_request(&frame).is_err());
    let decoded = api::decode_request(&frame).unwrap();
    let response = collagene_library_host::discussion::dispatch(db, m, &decoded).unwrap();
    let response = api::decode_response(&api::encode_response(&response).unwrap()).unwrap();
    response
        .validate_for(&request, &wire::LabId::new(m.lab_id().as_str()).unwrap())
        .unwrap();
    response
}
#[test]
fn sdk_roundtrip_uses_authenticated_author_and_exact_revision_with_typed_conflicts() {
    let mut db = store::Library::open_in_memory().unwrap();
    let lab = db.create_lab("Lab").unwrap();
    let access = db.trusted_host_access(&lab.id).unwrap();
    let a = member(&db, &access, "scientist-a", native::Role::Contributor);
    let reader = member(&db, &access, "reader", native::Role::Reader);
    let bytes = b"opaque stored archive".to_vec();
    let payload = store::Payload {
        kind: store::PayloadKind::NativeDocument,
        archive_hash: store::Digest::parse(wire::Digest::of(&bytes).as_str()).unwrap(),
        bytes,
        selection: None,
        sharing: None,
    };
    let doc = db
        .create_document(&access, None, "pExample", &payload)
        .unwrap();
    let target = api::Target {
        document: wire::DocumentId::new(doc.id.as_str()).unwrap(),
        revision: wire::RevisionId::new(doc.head.id.as_str()).unwrap(),
    };
    let identified = call(&mut db, &a, api::Request::Identify);
    assert_eq!(identified.member.author.id.as_str(), "scientist-a");
    let request = api::Request::Append {
        target: target.clone(),
        expected_generation: 0,
        expected_actor: wire::Text::new("scientist-a").unwrap(),
        action: api::Action::OpenThread {
            body: wire::Text::new("<b>Plain comment</b>").unwrap(),
        },
    };
    let appended = call(&mut db, &a, request.clone());
    assert!(
        matches!(appended.outcome,api::Outcome::Success(api::Success::Appended{event,..}) if event.author.id.as_str()=="scientist-a"&&event.sequence==1)
    );
    assert!(matches!(
        call(&mut db, &a, request.clone()).outcome,
        api::Outcome::Failure(api::Failure {
            code: wire::ErrorCode::Conflict,
            actual_generation: Some(1)
        })
    ));
    assert!(matches!(
        call(&mut db, &reader, request).outcome,
        api::Outcome::Failure(api::Failure {
            code: wire::ErrorCode::Forbidden,
            actual_generation: None
        })
    ));
    let other_person = member(&db, &access, "scientist-b", native::Role::Contributor);
    let wrong_identity = api::Request::Append {
        target: target.clone(),
        expected_generation: 1,
        expected_actor: wire::Text::new("scientist-a").unwrap(),
        action: api::Action::Review {
            decision: api::ReviewDecision::Approved,
            body: wire::Text::new("").unwrap(),
        },
    };
    assert!(matches!(
        call(&mut db, &other_person, wrong_identity).outcome,
        api::Outcome::Failure(api::Failure {
            code: wire::ErrorCode::Forbidden,
            ..
        })
    ));
    let page = call(
        &mut db,
        &reader,
        api::Request::List {
            target: target.clone(),
            offset: 0,
            limit: 1,
            expected_generation: None,
        },
    );
    let api::Outcome::Success(api::Success::Listed(page)) = page.outcome else {
        panic!("list")
    };
    assert_eq!(page.events.len(), 1);
    assert_eq!(page.events[0].author.id.as_str(), "scientist-a");
    let other = db.create_lab("Other").unwrap();
    let other = db.trusted_host_access(&other.id).unwrap();
    let foreign = member(&db, &other, "scientist-a", native::Role::Maintainer);
    assert!(matches!(
        call(
            &mut db,
            &foreign,
            api::Request::List {
                target,
                offset: 0,
                limit: 1,
                expected_generation: None
            }
        )
        .outcome,
        api::Outcome::Failure(api::Failure {
            code: wire::ErrorCode::NotFound,
            ..
        })
    ));
}
#[test]
fn bounded_frames_and_response_binding_refuse_changed_author_target_action_and_paging() {
    let lab = wire::LabId::new(format!("lab-{}", "a".repeat(32))).unwrap();
    let target = api::Target {
        document: wire::DocumentId::new(format!("doc-{}", "b".repeat(32))).unwrap(),
        revision: wire::RevisionId::new("c".repeat(64)).unwrap(),
    };
    let member = api::Member {
        author: api::Author {
            id: wire::Text::new("a").unwrap(),
            display_name: wire::Text::new("Scientist").unwrap(),
        },
        role: api::Role::Contributor,
    };
    let action = api::Action::Review {
        decision: api::ReviewDecision::Approved,
        body: wire::Text::new("").unwrap(),
    };
    let request = api::RequestEnvelope {
        request_id: 9,
        request: api::Request::Append {
            target: target.clone(),
            expected_generation: 0,
            expected_actor: member.author.id.clone(),
            action: action.clone(),
        },
    };
    let response = api::ResponseEnvelope {
        request_id: 9,
        context: wire::LabContext {
            lab_id: lab.clone(),
            name: wire::Text::new("Lab").unwrap(),
        },
        member: member.clone(),
        outcome: api::Outcome::Success(api::Success::Appended {
            target: target.clone(),
            event: api::Event {
                sequence: 1,
                author: member.author.clone(),
                role: member.role,
                created_at: 1,
                action,
            },
        }),
    };
    response.validate_for(&request, &lab).unwrap();
    let mut changed = response.clone();
    if let api::Outcome::Success(api::Success::Appended { event, .. }) = &mut changed.outcome {
        event.author.id = wire::Text::new("b").unwrap();
    }
    assert!(changed.validate_for(&request, &lab).is_err());
    let mut changed = response.clone();
    if let api::Outcome::Success(api::Success::Appended { target, .. }) = &mut changed.outcome {
        target.revision = wire::RevisionId::new("d".repeat(64)).unwrap();
    }
    assert!(changed.validate_for(&request, &lab).is_err());
    let mut changed = response.clone();
    if let api::Outcome::Success(api::Success::Appended { event, .. }) = &mut changed.outcome {
        event.action = api::Action::OpenThread {
            body: wire::Text::new("different").unwrap(),
        };
    }
    assert!(changed.validate_for(&request, &lab).is_err());
    let mut bytes = api::encode_request(&request).unwrap();
    bytes[4] = 1;
    assert!(matches!(
        api::decode_request(&bytes),
        Err(wire::Error::Version(1))
    ));
    assert!(api::decode_request(&vec![0; api::MAX_REQUEST_BYTES + 1]).is_err());
    assert!(
        api::Action::OpenThread {
            body: wire::Text::new("  ").unwrap()
        }
        .validate()
        .is_err()
    );
    assert!(wire::Text::<8192>::new("é".repeat(4097)).is_err());
    let request = api::RequestEnvelope {
        request_id: 9,
        request: api::Request::List {
            target: target.clone(),
            offset: 0,
            limit: 1,
            expected_generation: None,
        },
    };
    let mut paged = response;
    let api::Outcome::Success(api::Success::Appended { event, .. }) = paged.outcome else {
        unreachable!()
    };
    paged.outcome = api::Outcome::Success(api::Success::Listed(api::Page {
        target,
        generation: 2,
        events: wire::Rows::new(vec![event]).unwrap(),
        next_offset: None,
    }));
    assert!(paged.validate_for(&request, &lab).is_err());
    if let api::Outcome::Success(api::Success::Listed(p)) = &mut paged.outcome {
        p.next_offset = Some(1);
    }
    paged.validate_for(&request, &lab).unwrap();
}
