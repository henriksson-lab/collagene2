//! Retained runs cross the binary boundary as data, with exact workspace bindings.
use collagene_core::{
    protocol::{
        self,
        request::Request,
        response::Response,
        schema::{AnalysisJob, AnalysisResult},
    },
    run_lineage::*,
    verification::*,
    workspace::*,
};
use collagene_model::{Error, Molecule, Topology};
fn fixture() -> (Workspace, DesignPlan, EventBinding, Vec<InputEventBinding>) {
    let source = Molecule::dna("source", "AACCGG", Topology::Linear).unwrap();
    let target = Molecule::dna("target", "AACCGG", Topology::Linear).unwrap();
    let w =
        Workspace::from_documents("run fixture", vec![source.clone(), target.clone()], 1).unwrap();
    let binding = |i: usize| EventBinding {
        document_id: w.documents[i].id.clone(),
        event_id: w.documents[i].head_event.clone().unwrap(),
        revision: w.documents[i].head.clone(),
    };
    let plan = DesignPlan {
        name: "reverse".into(),
        target_revision: target.revision(),
        inputs: vec![PlanInput {
            id: "s".into(),
            revision: source.revision(),
            molecule: source,
        }],
        steps: vec![Operation::Overlap {
            fragments: vec![Source::Input("s".into())],
            overlaps: vec![],
            circular: false,
        }],
    };
    let target = binding(1);
    let inputs = vec![InputEventBinding {
        input_id: "s".into(),
        event: binding(0),
    }];
    (w, plan, target, inputs)
}
fn envelope(magic: &[u8; 4], version: u16, value: &impl serde::Serialize) -> Vec<u8> {
    let p = postcard::to_allocvec(value).unwrap();
    raw(magic, version, &p)
}
fn raw(magic: &[u8; 4], version: u16, p: &[u8]) -> Vec<u8> {
    let mut b = magic.to_vec();
    b.extend(version.to_le_bytes());
    b.extend((p.len() as u32).to_le_bytes());
    b.extend(blake3::hash(p).as_bytes());
    b.extend(p);
    b
}
fn record(
    w: Workspace,
    plan: DesignPlan,
    target: EventBinding,
    inputs: Vec<InputEventBinding>,
) -> Request {
    Request {
        id: 103,
        revision: workspace_fingerprint(&w).unwrap(),
        job: AnalysisJob::RecordCloningRun {
            workspace: w,
            plan,
            target_binding: target,
            input_bindings: inputs,
            allow_incomplete: false,
        },
    }
}
#[test]
fn record_replay_source_binding_and_current_schema_roundtrip() {
    let (w, p, t, i) = fixture();
    let original = workspace_fingerprint(&w).unwrap();
    let req = record(w, p, t, i);
    let encoded = protocol::request::Request::encode(&req).unwrap();
    assert_eq!(
        &encoded[4..6],
        &collagene_core::protocol::schema::WORKER_SCHEMA_VERSION.to_le_bytes()
    );
    for v in 1..12 {
        assert!(protocol::execution::process_request(&envelope(b"CGQ\0", v, &req)).is_err());
    }
    let result = protocol::response::Response::decode(
        &protocol::execution::process_request(&encoded).unwrap(),
    )
    .unwrap();
    assert_eq!(result.revision, original);
    let AnalysisResult::RunCommit(commit) = result.result else {
        panic!("not a run commit")
    };
    assert_eq!(commit.receipt.products.len(), 1);
    let response = Response {
        id: 103,
        revision: original.clone(),
        result: AnalysisResult::RunCommit(commit.clone()),
    };
    for v in 1..12 {
        assert!(protocol::response::Response::decode(&envelope(b"CGR\0", v, &response)).is_err());
    }
    let save = Request {
        id: 103,
        revision: commit
            .workspace
            .head(&commit.workspace.active_document)
            .unwrap()
            .revision(),
        job: AnalysisJob::SaveWorkspace {
            workspace: commit.workspace.clone(),
        },
    };
    assert!(protocol::execution::process_request(&envelope(b"CGQ\0", 11, &save)).is_err());
    let loaded = Response {
        id: 103,
        revision: original.clone(),
        result: AnalysisResult::WorkspaceLoaded(commit.workspace.clone()),
    };
    assert!(protocol::response::Response::decode(&envelope(b"CGR\0", 11, &loaded)).is_err());
    let mut wrong = response;
    wrong.revision = "0".repeat(64);
    assert!(
        protocol::response::Response::decode(&envelope(
            b"CGR\0",
            protocol::schema::WORKER_SCHEMA_VERSION,
            &wrong
        ))
        .is_err()
    );
    let replay = Request {
        id: 104,
        revision: workspace_fingerprint(&commit.workspace).unwrap(),
        job: AnalysisJob::ReplayCloningRun {
            workspace: commit.workspace,
            run_id: commit.receipt.run_id,
        },
    };
    let r = protocol::response::Response::decode(
        &protocol::execution::process_request(
            &protocol::request::Request::encode(&replay).unwrap(),
        )
        .unwrap(),
    )
    .unwrap();
    let AnalysisResult::RunReplay(report) = r.result else {
        panic!()
    };
    assert!(report.products_match && report.status_matches && report.engine_matches);
    let response = Response {
        id: 104,
        revision: r.revision.clone(),
        result: AnalysisResult::RunReplay(report),
    };
    assert!(protocol::response::Response::decode(&envelope(b"CGR\0", 11, &response)).is_err());
    let mut substituted = response;
    substituted.revision = "0".repeat(64);
    assert!(
        protocol::response::Response::decode(&envelope(
            b"CGR\0",
            protocol::schema::WORKER_SCHEMA_VERSION,
            &substituted
        ))
        .is_err()
    );
    let mut wrong = replay;
    wrong.revision = original;
    assert!(
        protocol::execution::process_request(&protocol::request::Request::encode(&wrong).unwrap())
            .is_err()
    );
}
#[test]
fn large_run_workspace_payloads_are_whitelisted_but_ordinary_current_payloads_remain_bounded() {
    let (mut w, p, t, i) = fixture();
    for value in [7u8, 9] {
        let bytes = vec![value; 16 * 1024 * 1024 + 1];
        w.attachments.insert(
            blake3::hash(&bytes).to_hex().to_string(),
            WorkspaceAttachment {
                kind: "evidence".into(),
                name: format!("large{value}"),
                media_type: "application/octet-stream".into(),
                bytes,
            },
        );
    }
    let request = record(w, p, t, i);
    let encoded = protocol::request::Request::encode(&request).unwrap();
    assert!(encoded.len() > 32 * 1024 * 1024);
    let bytes = protocol::execution::process_request(&encoded).unwrap();
    assert!(bytes.len() > 32 * 1024 * 1024);
    let AnalysisResult::RunCommit(commit) =
        protocol::response::Response::decode(&bytes).unwrap().result
    else {
        panic!()
    };
    let replay = Request {
        id: 7,
        revision: workspace_fingerprint(&commit.workspace).unwrap(),
        job: AnalysisJob::ReplayCloningRun {
            workspace: commit.workspace,
            run_id: commit.receipt.run_id,
        },
    };
    assert!(protocol::request::Request::encode(&replay).unwrap().len() > 32 * 1024 * 1024);
    assert!(
        protocol::execution::process_request(&protocol::request::Request::encode(&replay).unwrap())
            .is_ok()
    );
    // Current ordinary variants remain bounded before their intentionally invalid
    // oversized bodies are deserialized. Run/workspace variants above are larger.
    for (magic, tags) in [(b"CGQ\0", vec![0]), (b"CGR\0", vec![2])] {
        for tag in tags {
            let mut p = vec![1, 1, b'x', tag];
            p.resize(32 * 1024 * 1024 + 1, 0);
            let b = raw(magic, protocol::schema::WORKER_SCHEMA_VERSION, &p);
            let r = if magic == b"CGQ\0" {
                protocol::execution::process_request(&b).map(|_| ())
            } else {
                protocol::response::Response::decode(&b).map(|_| ())
            };
            assert!(matches!(r, Err(Error::TooLarge)));
        }
    }
}
