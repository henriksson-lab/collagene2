//! Stateful batch messages retain one private source frame. Caller-supplied
//! historical summaries cannot populate that private computation state.
use collagene_core::{
    assembly_recipe::Preparation,
    chemistry::Terminus,
    combinatorial_batch::{self as batch, RowRequest, RowState},
    combinatorial_design::{
        Bin, Choice, ChoiceSelection, Design, Method, SelectedConstruct, TargetDeclaration,
    },
    protocol::{
        self,
        request::Request,
        response::Response,
        schema::{AnalysisJob, AnalysisResult},
        session::AnalysisSession,
    },
    run_lineage::EventBinding,
    verification::Status,
    workspace::Workspace,
};
use collagene_model::{Error, Molecule, Topology};

fn binding(workspace: &Workspace, index: usize) -> EventBinding {
    let doc = &workspace.documents[index];
    EventBinding {
        document_id: doc.id.clone(),
        event_id: doc.head_event.clone().unwrap(),
        revision: doc.head.clone(),
    }
}
fn fixture() -> (Workspace, batch::Request) {
    let workspace = Workspace::from_documents(
        "sequential batch",
        vec![
            Molecule::dna("backbone", "AACC", Topology::Linear).unwrap(),
            Molecule::dna("insert", "TTAA", Topology::Linear).unwrap(),
            Molecule::dna("independent wrong target", "GGGG", Topology::Linear).unwrap(),
        ],
        0,
    )
    .unwrap();
    let choices = vec![
        ChoiceSelection {
            bin_id: "backbone".into(),
            choice_id: "source-a".into(),
        },
        ChoiceSelection {
            bin_id: "insert".into(),
            choice_id: "source-b".into(),
        },
    ];
    let ready = |id: &str, target| {
        RowRequest::Ready(SelectedConstruct {
            id: id.into(),
            name: id.into(),
            choices: choices.clone(),
            target,
            overrides: None,
        })
    };
    let matrix = batch::Request {
        design: Some(Design {
            material_bindings: vec![],
            id: "wire-library".into(),
            name: "Declared duplex alternatives".into(),
            circular: false,
            bins: [("backbone", "source-a", 0), ("insert", "source-b", 1)]
                .into_iter()
                .map(|(bin, id, index)| Bin {
                    id: bin.into(),
                    label: bin.into(),
                    choices: vec![Choice {
                        generation: None,
                        id: id.into(),
                        label: id.into(),
                        source: binding(&workspace, index),
                        selector: collagene_core::assembly_intent::Selector::ReferenceRange {
                            start: 0,
                            end: 4,
                            reverse: false,
                        },
                        preparation_reverse: false,
                        preparation: Some(Preparation::DeclaredDuplex {
                            left: Terminus::CLEAVED,
                            right: Terminus::CLEAVED,
                        }),
                        treatment: None,
                    }],
                })
                .collect(),
            method: Method::Ligation {
                allow_nicked: false,
            },
        }),
        rows: vec![
            ready(
                "construct-good",
                TargetDeclaration::OriginalRegions {
                    substitutions: vec![],
                    edits: vec![],
                    junction_edits: vec![],
                },
            ),
            RowRequest::DraftError {
                id: "construct-unfinished".into(),
                name: "Unfinished primer".into(),
                choices: choices.clone(),
                message: "Primer draft is not valid DNA.".into(),
            },
            ready(
                "construct-wrong",
                TargetDeclaration::Existing(binding(&workspace, 2)),
            ),
        ],
    };
    matrix.validate().unwrap();
    (workspace, matrix)
}
fn begin(workspace: &Workspace, matrix: &batch::Request) -> Request {
    Request {
        id: 128,
        revision: batch::begin_revision(workspace, matrix).unwrap(),
        job: AnalysisJob::BeginCombinatorialBatch {
            workspace: workspace.clone(),
            matrix: matrix.clone(),
        },
    }
}
fn advance(context: &str, index: usize) -> Request {
    Request {
        id: 129 + index as u32,
        revision: context.into(),
        job: AnalysisJob::AdvanceCombinatorialBatch {
            context: context.into(),
            index,
        },
    }
}
fn capture(workspace: &Workspace, context: &str) -> Request {
    Request {
        id: 130,
        revision: batch::capture_revision(workspace, context).unwrap(),
        job: AnalysisJob::CaptureCombinatorialBatch {
            workspace: workspace.clone(),
            current_context: context.into(),
        },
    }
}
fn selector(wire: &[u8]) -> u32 {
    let ((_, _, selector), _) = postcard::take_from_bytes::<(u32, &str, u32)>(&wire[42..]).unwrap();
    selector
}
fn ask(session: &mut AnalysisSession, request: &Request) -> (Vec<u8>, Response) {
    let wire = session
        .process_request(&protocol::request::Request::encode(request).unwrap())
        .unwrap();
    let response = protocol::response::Response::decode(&wire).unwrap();
    assert_eq!(response.id, request.id);
    assert_eq!(response.revision, request.revision);
    (wire, response)
}
fn response_wire(response: &Response) -> Vec<u8> {
    let body = postcard::to_allocvec(response).unwrap();
    let mut wire = b"CGR\0".to_vec();
    wire.extend_from_slice(&protocol::schema::WORKER_SCHEMA_VERSION.to_le_bytes());
    wire.extend_from_slice(&(body.len() as u32).to_le_bytes());
    wire.extend_from_slice(blake3::hash(&body).as_bytes());
    wire.extend(body);
    wire
}

#[test]
fn one_private_source_frame_runs_each_row_and_captures_historical_evidence_only() {
    let (workspace, matrix) = fixture();
    let request = begin(&workspace, &matrix);
    let context = request.revision.clone();
    let mut session = AnalysisSession::default();
    assert_eq!(
        selector(&protocol::request::Request::encode(&request).unwrap()),
        126
    );
    let (wire, response) = ask(&mut session, &request);
    assert_eq!(selector(&wire), 118);
    let AnalysisResult::CombinatorialBatchStarted(start) = response.result else {
        panic!("missing batch start")
    };
    assert_eq!(
        (start.next_index, start.total_rows, start.finished),
        (0, 3, false)
    );
    assert!(start.row.is_none());
    assert!(
        session
            .process_request(&protocol::request::Request::encode(&advance(&context, 2)).unwrap())
            .is_err()
    );
    assert!(
        session
            .process_request(
                &protocol::request::Request::encode(&capture(&workspace, &context)).unwrap()
            )
            .is_err()
    );
    let mut evidence = vec![];
    for index in 0..3 {
        let (wire, response) = ask(&mut session, &advance(&context, index));
        assert_eq!(selector(&wire), 119);
        let AnalysisResult::CombinatorialBatchAdvanced(progress) = response.result else {
            panic!("missing row progress")
        };
        assert_eq!(progress.next_index, index + 1);
        assert_eq!(progress.total_rows, 3);
        assert_eq!(progress.finished, index == 2);
        let row = progress.row.unwrap();
        assert_eq!(row.index, index);
        evidence.push(row);
    }
    assert_eq!(evidence[0].construct_id, "construct-good");
    assert_eq!(evidence[0].state, RowState::Evaluated);
    assert_eq!(evidence[0].status, Some(Status::ConsistentUnderAssumptions));
    assert!(evidence[0].final_product_revision.is_some());
    assert_eq!(evidence[0].input_bindings[0].event, binding(&workspace, 0));
    assert_eq!(evidence[0].input_bindings[1].event, binding(&workspace, 1));
    assert_eq!(evidence[1].state, RowState::DraftError);
    assert!(evidence[1].status.is_none());
    assert!(evidence[1].final_product_revision.is_none());
    assert_eq!(evidence[2].state, RowState::Evaluated);
    assert_eq!(evidence[2].status, Some(Status::Invalid));
    assert!(evidence[2].final_product_revision.is_some());
    let (_, response) = ask(&mut session, &advance(&context, 0));
    let AnalysisResult::CombinatorialBatchAdvanced(retry) = response.result else {
        panic!("missing retry")
    };
    assert_eq!(retry.next_index, 3);
    assert_eq!(retry.row.unwrap().construct_id, "construct-good");
    let (wire, response) = ask(&mut session, &capture(&workspace, &context));
    assert_eq!(selector(&wire), 120);
    let AnalysisResult::CombinatorialBatchCaptured(commit) = response.result else {
        panic!("missing private capture")
    };
    commit
        .validate_binding(&batch::capture_revision(&workspace, &context).unwrap())
        .unwrap();
    for (before, after) in workspace.documents.iter().zip(&commit.workspace.documents) {
        assert_eq!(before.head, after.head);
        assert_eq!(before.head_event, after.head_event);
    }
    let record = batch::load(&commit.workspace, &commit.attachment_hash).unwrap();
    assert_eq!(record.rows.len(), 3);
    assert_eq!(record.rows[0].construct_id, "construct-good");
    assert_eq!(record.rows[1].state, RowState::DraftError);
    assert_eq!(record.rows[2].status, Some(Status::Invalid));
    let mut corrupted = commit.workspace;
    let asset = corrupted
        .attachments
        .get_mut(&commit.attachment_hash)
        .unwrap();
    asset.bytes[0] ^= 1;
    assert!(batch::load(&corrupted, &commit.attachment_hash).is_err());
}

#[test]
fn private_session_missing_stale_or_replaced_context_cannot_supply_receipts() {
    let (workspace, matrix) = fixture();
    let request = begin(&workspace, &matrix);
    let context = request.revision.clone();
    let mut session = AnalysisSession::default();
    assert!(
        protocol::execution::process_request(
            &protocol::request::Request::encode(&request).unwrap()
        )
        .is_err()
    );
    assert!(
        session
            .process_request(&protocol::request::Request::encode(&advance(&context, 0)).unwrap())
            .is_err()
    );
    ask(&mut session, &request);
    let mut bad = request.clone();
    bad.revision = "0".repeat(64);
    assert!(
        session
            .process_request(&protocol::request::Request::encode(&bad).unwrap())
            .is_err()
    );
    ask(&mut session, &advance(&context, 0));
    let mut stale = advance(&context, 1);
    stale.revision = "0".repeat(64);
    assert!(
        session
            .process_request(&protocol::request::Request::encode(&stale).unwrap())
            .is_err()
    );
    assert!(
        session
            .process_request(
                &protocol::request::Request::encode(&advance(&"0".repeat(64), 1)).unwrap()
            )
            .is_err()
    );
    let mut other = matrix.clone();
    if let RowRequest::Ready(row) = &mut other.rows[0] {
        row.name = "Changed declaration".into();
    }
    let replacement = begin(&workspace, &other);
    assert_ne!(replacement.revision, context);
    ask(&mut session, &replacement);
    assert!(
        session
            .process_request(&protocol::request::Request::encode(&advance(&context, 0)).unwrap())
            .is_err()
    );
    ask(&mut session, &advance(&replacement.revision, 0));
}

#[test]
fn schema56_full_and_borrowed_prefixes_cannot_hide_variable_length_batch_selectors() {
    let (workspace, matrix) = fixture();
    let request = begin(&workspace, &matrix);
    let context = request.revision.clone();
    let mut session = AnalysisSession::default();
    let (begin_reply, _) = ask(&mut session, &request);
    let (advance_reply, _) = ask(&mut session, &advance(&context, 0));
    ask(&mut session, &advance(&context, 1));
    ask(&mut session, &advance(&context, 2));
    let (capture_reply, _) = ask(&mut session, &capture(&workspace, &context));
    let pairs = [
        (
            protocol::request::Request::encode(&request).unwrap(),
            false,
            126u32,
        ),
        (
            protocol::request::Request::encode(&advance(&context, 0)).unwrap(),
            false,
            127,
        ),
        (
            protocol::request::Request::encode(&capture(&workspace, &context)).unwrap(),
            false,
            128,
        ),
        (begin_reply, true, 118),
        (advance_reply, true, 119),
        (capture_reply, true, 120),
    ];
    for (mut wire, response, tag) in pairs {
        assert_eq!(selector(&wire), tag);
        wire[4..6].copy_from_slice(&56u16.to_le_bytes());
        let rejected = if response {
            protocol::response::Response::decode(&wire).map(|_| ())
        } else {
            session.process_request(&wire).map(|_| ())
        };
        assert!(matches!(rejected, Err(Error::Unsupported(_))));
        let body = postcard::to_allocvec(&(7u32, "context", tag)).unwrap();
        if tag >= 128 {
            assert!(body.ends_with(&[0x80 | (tag as u8 & 0x7f), 1]));
        }
        wire.truncate(42);
        wire[6..10].copy_from_slice(&(body.len() as u32).to_le_bytes());
        wire.extend(body);
        let rejected = if response {
            protocol::response::Response::decode(&wire).map(|_| ())
        } else {
            session.process_request(&wire).map(|_| ())
        };
        assert!(matches!(rejected, Err(Error::Unsupported(_))));
    }
}

#[test]
fn response_shape_and_source_tokens_are_checked_without_restoring_private_state() {
    let (workspace, matrix) = fixture();
    let request = begin(&workspace, &matrix);
    let context = request.revision.clone();
    let mut session = AnalysisSession::default();
    let (_, start) = ask(&mut session, &request);
    let mut bad = start;
    bad.revision = "0".repeat(64);
    assert!(protocol::response::Response::decode(&response_wire(&bad)).is_err());
    let (_, progress) = ask(&mut session, &advance(&context, 0));
    let mut bad = progress;
    let AnalysisResult::CombinatorialBatchAdvanced(row) = &mut bad.result else {
        unreachable!()
    };
    row.next_index = usize::MAX;
    assert!(protocol::response::Response::decode(&response_wire(&bad)).is_err());
    ask(&mut session, &advance(&context, 1));
    ask(&mut session, &advance(&context, 2));
    let (_, commit) = ask(&mut session, &capture(&workspace, &context));
    let mut bad = commit.clone();
    bad.revision = "0".repeat(64);
    assert!(protocol::response::Response::decode(&response_wire(&bad)).is_err());
    let mut fresh = AnalysisSession::default();
    assert!(
        fresh
            .process_request(
                &protocol::request::Request::encode(&capture(&workspace, &context)).unwrap()
            )
            .is_err()
    );
    let AnalysisResult::CombinatorialBatchCaptured(saved) = commit.result else {
        unreachable!()
    };
    assert!(batch::load(&saved.workspace, &saved.attachment_hash).is_ok());
}

#[test]
fn ordinary_progress_limits_and_explicit_workspace_budget_whitelist_remain_distinct() {
    let maximum = collagene_core::formats::MAX_IMPORT_BYTES;
    let mut wire = vec![0; 42 + maximum + 1];
    wire[4..6].copy_from_slice(&protocol::schema::WORKER_SCHEMA_VERSION.to_le_bytes());
    wire[6..10].copy_from_slice(&((maximum + 1) as u32).to_le_bytes());
    for (magic, tag, response) in [
        (b"CGQ\0", 127u32, false),
        (b"CGR\0", 118u32, true),
        (b"CGR\0", 119u32, true),
    ] {
        wire[..4].copy_from_slice(magic);
        let body = postcard::to_allocvec(&(7u32, "context", tag)).unwrap();
        wire[42..42 + body.len()].copy_from_slice(&body);
        let rejected = if response {
            protocol::response::Response::decode(&wire).map(|_| ())
        } else {
            AnalysisSession::default()
                .process_request(&wire)
                .map(|_| ())
        };
        assert!(matches!(rejected, Err(Error::TooLarge)));
    }
    for (magic, tag, response) in [
        (b"CGQ\0", 126u32, false),
        (b"CGQ\0", 128u32, false),
        (b"CGR\0", 120u32, true),
    ] {
        wire[..4].copy_from_slice(magic);
        let body = postcard::to_allocvec(&(7u32, "context", tag)).unwrap();
        wire[42..42 + body.len()].copy_from_slice(&body);
        // These selectors may use workspace budget but still reject the invalid
        // checksum/body. They do not accidentally inherit ordinary-job routing.
        let rejected = if response {
            protocol::response::Response::decode(&wire).map(|_| ())
        } else {
            AnalysisSession::default()
                .process_request(&wire)
                .map(|_| ())
        };
        assert!(matches!(rejected, Err(Error::Format(_))));
    }
}
