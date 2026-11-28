//! Binary transport preserves finite declared states and immutable source identity.
use collagene_core::{
    protocol::{
        self,
        request::Request,
        response::Response,
        schema::{AnalysisJob, AnalysisResult},
    },
    sequence::crispr::{self, base_editing},
};
use collagene_model::{Molecule, Topology};

fn inputs(limit: usize) -> (Molecule, base_editing::Request) {
    let source = Molecule::dna(
        "Literal base-editor source",
        "AAACCCCAAAAAAAAAAAAAAGG",
        Topology::Linear,
    )
    .unwrap();
    let guide_options = crispr::GuideOptions::default();
    let guide = crispr::design_guides(&source, &guide_options)
        .unwrap()
        .candidates
        .into_iter()
        .find(|candidate| !candidate.reverse && candidate.guide_sequence == "AAACCCCAAAAAAAAAAAAA")
        .unwrap();
    let request = base_editing::Request {
        guide_options,
        guide,
        profile: base_editing::be1_2016().unwrap(),
        intended_positions: vec![4],
        max_outcomes: limit,
        cds: None,
    };
    (source, request)
}

#[test]
fn base_editor_worker_keeps_finite_caps_and_biological_nonexhaustiveness_separate() {
    let mut session = protocol::session::AnalysisSession::default();
    for retained in [3, 16] {
        let (source, declaration) = inputs(retained);
        let before = source.clone();
        let request = Request {
            id: 91,
            revision: base_editing::revision(&source, &declaration).unwrap(),
            job: AnalysisJob::AnalyzeBaseEditing {
                molecule: source.clone(),
                request: declaration.clone(),
            },
        };
        let bytes = request.encode().unwrap();
        for response in [
            protocol::execution::process_request(&bytes).unwrap(),
            session.process_request(&bytes).unwrap(),
        ] {
            let response = Response::decode(&response).unwrap();
            assert_eq!(
                (response.id, &response.revision),
                (request.id, &request.revision)
            );
            let AnalysisResult::BaseEditing(report) = response.result else {
                panic!("base-editor report");
            };
            report.validate_for(&source, &declaration).unwrap();
            assert_eq!(report.finite_state_count, 16);
            assert_eq!(report.states.len(), retained);
            assert_eq!(report.finite_enumeration_complete, retained == 16);
            assert!(!report.biological_outcomes_exhaustive);
            assert_eq!(report.states[0].reference_guide, "AAACCCCAAAAAAAAAAAAA");
            assert!(report.states[0].edited_positions.is_empty());
            assert_eq!(report.states[0].satisfies_intent, Some(false));
            assert!(
                report.cds.is_none(),
                "missing explicit CDS assumptions never fabricate consequences"
            );
        }
        assert_eq!(source, before);
    }
}

#[test]
fn base_editor_worker_refuses_foreign_source_context_and_altered_retained_states() {
    let (source, declaration) = inputs(16);
    let request = Request {
        id: 92,
        revision: base_editing::revision(&source, &declaration).unwrap(),
        job: AnalysisJob::AnalyzeBaseEditing {
            molecule: source.clone(),
            request: declaration,
        },
    };
    let mut foreign = request.clone();
    let AnalysisJob::AnalyzeBaseEditing { molecule, .. } = &mut foreign.job else {
        unreachable!()
    };
    molecule.name = "Another exact source, identical DNA letters".into();
    assert!(foreign.encode().is_err());
    let mut wrong_context = request.clone();
    wrong_context.revision = "0".repeat(64);
    assert!(wrong_context.encode().is_err());
    let mut response = Response::decode(
        &protocol::execution::process_request(&request.encode().unwrap()).unwrap(),
    )
    .unwrap();
    let AnalysisResult::BaseEditing(report) = &mut response.result else {
        unreachable!()
    };
    report.states[0].reference_guide.replace_range(0..1, "T");
    assert!(Response::decode(&response.encode().unwrap()).is_err());
}
