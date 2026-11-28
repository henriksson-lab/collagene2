//! Semantic binary selectors retain exact independent source roles and path evidence.
use collagene_core::{
    alignment::{ComparisonOptions, design_consensus as native},
    documents::workspace::Workspace,
    material::lineage::EventBinding,
    protocol::{
        self,
        request::Request,
        response::Response,
        schema::{AnalysisJob, AnalysisResult, WORKER_SCHEMA_VERSION},
    },
};
use collagene_model::{Error, Molecule, Topology};

fn fixture() -> (Workspace, native::Request) {
    let workspace = Workspace::from_documents(
        "Independent supplied consensus",
        vec![
            Molecule::dna("separately bound design", "ACGTT", Topology::Linear).unwrap(),
            Molecule::dna("supplied consensus", "ACATT", Topology::Linear).unwrap(),
        ],
        1,
    )
    .unwrap();
    let event = |index: usize| {
        let document = &workspace.documents[index];
        EventBinding {
            document_id: document.id.clone(),
            event_id: document.head_event.clone().unwrap(),
            revision: document.head.clone(),
        }
    };
    let declaration = native::Request {
        schema_version: 1,
        design: event(0),
        consensus: event(1),
        comparison: ComparisonOptions {
            both_query_strands: false,
            reference_rotations: false,
            max_origins: 1,
            max_total_cells: 100_000,
            ..Default::default()
        },
        consensus_axis: native::ConsensusAxis::StoredTopology,
        cds: vec![],
        limits: Default::default(),
    };
    (workspace, declaration)
}
fn request(workspace: Workspace, declaration: native::Request) -> Request {
    Request {
        id: 74,
        revision: native::revision(&workspace, &declaration).unwrap(),
        job: AnalysisJob::AssessDesignConsensus {
            workspace,
            request: declaration,
        },
    }
}
fn raw<T: serde::Serialize>(magic: &[u8; 4], value: &T) -> Vec<u8> {
    let payload = postcard::to_allocvec(value).unwrap();
    let mut bytes = magic.to_vec();
    bytes.extend_from_slice(&WORKER_SCHEMA_VERSION.to_le_bytes());
    bytes.extend_from_slice(&(payload.len() as u32).to_le_bytes());
    bytes.extend_from_slice(blake3::hash(&payload).as_bytes());
    bytes.extend_from_slice(&payload);
    bytes
}

#[test]
fn stateless_and_session_semantic_dispatch_preserve_literal_path_and_exact_roles() {
    let (workspace, declaration) = fixture();
    let request = request(workspace.clone(), declaration.clone());
    let wire = request.encode().unwrap();
    let ((_, _, selector), _) = postcard::take_from_bytes::<(u32, &str, u32)>(&wire[42..]).unwrap();
    assert_eq!(selector, 214);
    let mut session = protocol::session::AnalysisSession::default();
    for bytes in [
        protocol::execution::process_request(&wire).unwrap(),
        session.process_request(&wire).unwrap(),
    ] {
        let ((_, _, selector), _) =
            postcard::take_from_bytes::<(u32, &str, u32)>(&bytes[42..]).unwrap();
        assert_eq!(selector, 204);
        let response = Response::decode(&bytes).unwrap();
        assert_eq!(response.id, request.id);
        assert_eq!(response.revision, request.revision);
        let AnalysisResult::DesignConsensus(report) = &response.result else {
            panic!("expected semantic design-consensus report")
        };
        report.validate_request(&workspace, &declaration).unwrap();
        report
            .evidence
            .validate_sources(&workspace, &report.pairwise)
            .unwrap();
        assert_eq!(report.evidence.request.design, declaration.design);
        assert_eq!(report.evidence.request.consensus, declaration.consensus);
        assert!(report.evidence.selected_path_conditional);
        let difference = &report.pairwise.comparison.alignment.differences[0];
        assert_eq!(
            (difference.reference_start, difference.reference_end),
            (2, 3)
        );
        assert_eq!((difference.query_start, difference.query_end), (2, 3));
        assert_eq!(difference.reference_bases, "G");
        assert_eq!(difference.query_bases, "A");
        let mut forged = response;
        let AnalysisResult::DesignConsensus(report) = &mut forged.result else {
            unreachable!()
        };
        report.pairwise.comparison.alignment.differences[0].query_bases = "T".into();
        // Correct framing/checksum does not authenticate modified path evidence.
        assert!(Response::decode(&raw(b"CGR\0", &forged)).is_err());
    }
}

#[test]
fn current_frames_refuse_stale_context_wrong_event_and_obsolete_schema() {
    let (workspace, declaration) = fixture();
    let request = request(workspace, declaration);
    let mut stale = request.clone();
    stale.revision = "0".repeat(64);
    assert!(stale.encode().is_err());
    assert!(protocol::execution::process_request(&raw(b"CGQ\0", &stale)).is_err());
    let mut wrong_source = request.clone();
    let AnalysisJob::AssessDesignConsensus {
        request: declaration,
        ..
    } = &mut wrong_source.job
    else {
        unreachable!()
    };
    declaration.consensus.document_id = declaration.design.document_id.clone();
    assert!(wrong_source.encode().is_err());
    assert!(protocol::execution::process_request(&raw(b"CGQ\0", &wrong_source)).is_err());
    let wire = request.encode().unwrap();
    let response = protocol::execution::process_request(&wire).unwrap();
    for (mut bytes, reply) in [(wire, false), (response, true)] {
        bytes[4..6].copy_from_slice(&(WORKER_SCHEMA_VERSION - 1).to_le_bytes());
        let result = if reply {
            Response::decode(&bytes).map(|_| ())
        } else {
            protocol::execution::process_request(&bytes).map(|_| ())
        };
        assert!(matches!(result, Err(Error::Unsupported(_))));
    }
}

#[test]
fn borrowed_semantic_reply_prefix_has_explicit_large_budget_without_bulk_allocation() {
    let payload_length = 32 * 1024 * 1024 + 1;
    let claimed_frame_length = 42 + payload_length;
    let prefix = |selector: u32| {
        let mut prefix = b"CGR\0".to_vec();
        prefix.extend_from_slice(&WORKER_SCHEMA_VERSION.to_le_bytes());
        prefix.extend_from_slice(&(payload_length as u32).to_le_bytes());
        prefix.extend_from_slice(&[0; 32]);
        prefix.extend(postcard::to_allocvec(&(74u32, "a".repeat(64), selector)).unwrap());
        prefix
    };
    let semantic = prefix(204);
    assert!(semantic.len() < 256);
    assert!(
        Response::admitted_frame_limit(&semantic, claimed_frame_length).unwrap() > payload_length
    );
    assert!(matches!(
        Response::admitted_frame_limit(&prefix(0), claimed_frame_length),
        Err(Error::TooLarge)
    ));
    // Prefix admission only chooses a resource cap; it cannot admit a truncated body.
    assert!(Response::decode(&semantic).is_err());
}
