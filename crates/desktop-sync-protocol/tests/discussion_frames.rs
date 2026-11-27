use collagene_desktop_sync_protocol::{
    self as wire,
    library::{self, discussion as api},
};
#[test]
fn full_discussion_page_binds_remote_actor_revision_and_current_desktop_schema() {
    let remote = wire::Remote {
        id: wire::RemoteId::new(format!("rem-{}", "a".repeat(32))).unwrap(),
        name: library::Text::new("Remote").unwrap(),
        endpoint: library::Text::new("https://lab.example").unwrap(),
        expected_lab: library::LabId::new(format!("lab-{}", "b".repeat(32))).unwrap(),
    };
    let target = api::Target {
        document: library::DocumentId::new(format!("doc-{}", "c".repeat(32))).unwrap(),
        revision: library::RevisionId::new("d".repeat(64)).unwrap(),
    };
    let member = api::Member {
        author: api::Author {
            id: library::Text::new("scientist").unwrap(),
            display_name: library::Text::new("Scientist").unwrap(),
        },
        role: api::Role::Contributor,
    };
    let request = wire::RequestEnvelope {
        request_id: 9,
        operation_id: wire::OperationId::new(format!("op-{}", "e".repeat(32))).unwrap(),
        command: wire::Command::Discussion {
            remote: remote.id.clone(),
            request: api::Request::List {
                target: target.clone(),
                offset: 0,
                limit: 200,
                expected_generation: None,
            },
        },
    };
    let context = wire::RemoteContext {
        remote: remote.clone(),
        lab: library::LabContext {
            lab_id: remote.expected_lab.clone(),
            name: library::Text::new("Lab").unwrap(),
        },
    };
    let rows = (1..=200)
        .map(|sequence| api::Event {
            sequence,
            author: member.author.clone(),
            role: member.role,
            created_at: 1,
            action: api::Action::OpenThread {
                body: library::Text::new("x".repeat(8192)).unwrap(),
            },
        })
        .collect();
    let response = wire::ResponseEnvelope {
        request_id: 9,
        operation_id: request.operation_id.clone(),
        context: Some(context.clone()),
        outcome: wire::Outcome::Success(wire::Success::Discussion(api::ResponseEnvelope {
            request_id: 9,
            context: context.lab,
            member,
            outcome: api::Outcome::Success(api::Success::Listed(api::Page {
                target,
                generation: 200,
                events: library::Rows::new(rows).unwrap(),
                next_offset: None,
            })),
        })),
    };
    let bytes = wire::encode_response(&response).unwrap();
    assert!(bytes.len() > 1 << 20);
    assert!(bytes.len() < wire::MAX_RESPONSE_BYTES);
    wire::decode_response_for(&bytes, &request)
        .unwrap()
        .validate_in_remote(&request, &remote)
        .unwrap();
    let mut changed = response.clone();
    if let wire::Outcome::Success(wire::Success::Discussion(reply)) = &mut changed.outcome {
        reply.context.lab_id = library::LabId::new(format!("lab-{}", "f".repeat(32))).unwrap();
    }
    assert!(wire::encode_response(&changed).is_ok());
    assert!(changed.validate_for(&request).is_err());
    let bytes = wire::encode_request(&request).unwrap();
    assert_eq!(&bytes[4..6], &4u16.to_le_bytes());
    let mut old = bytes;
    old[4..6].copy_from_slice(&3u16.to_le_bytes());
    assert!(matches!(
        wire::decode_request(&old),
        Err(wire::Error::Version(3))
    ));
}
