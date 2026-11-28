use collagene_core::{
    protocol::{
        self,
        request::Request,
        schema::{AnalysisJob, AnalysisResult},
    },
    sequence_similarity as native,
};
use collagene_model::{Molecule, Topology};
fn inputs() -> (native::Source, Vec<native::Source>, native::Options) {
    (
        native::Source {
            key: "query".into(),
            molecule: Molecule::dna("query", "CGAA", Topology::Linear).unwrap(),
        },
        vec![native::Source {
            key: "ring".into(),
            molecule: Molecule::dna("ring", "AATTCCG", Topology::Circular).unwrap(),
        }],
        native::Options::default(),
    )
}
#[test]
fn current_worker_job_and_retained_session_bind_sources_and_original_coordinates() {
    let (query, subjects, options) = inputs();
    let revision = native::context_revision(&query, &subjects, &options).unwrap();
    let request = Request {
        id: 7,
        revision: revision.clone(),
        job: AnalysisJob::FindSimilarSequences {
            query,
            subjects,
            options,
        },
    };
    let bytes = protocol::request::Request::encode(&request).unwrap();
    assert_eq!(
        &bytes[4..6],
        &protocol::schema::WORKER_SCHEMA_VERSION.to_le_bytes()
    );
    let ((_, _, tag), _) = postcard::take_from_bytes::<(u32, &str, u32)>(&bytes[42..]).unwrap();
    assert_eq!(tag, 191);
    let mut session = protocol::session::AnalysisSession::default();
    for reply in [
        protocol::execution::process_request(&bytes).unwrap(),
        session.process_request(&bytes).unwrap(),
    ] {
        let ((_, _, tag), _) = postcard::take_from_bytes::<(u32, &str, u32)>(&reply[42..]).unwrap();
        assert_eq!(tag, 181);
        let response = protocol::response::Response::decode(&reply).unwrap();
        let AnalysisResult::SequenceSimilarity(report) = response.result else {
            panic!("search result")
        };
        report.validate_binding(&revision).unwrap();
        assert!(report.complete);
        assert_eq!(report.hits[0].score, 8);
        assert_eq!(report.hits[0].reference_spans.len(), 2);
        assert_eq!(report.subjects[0].length, 7);
    }
    let mut old = bytes;
    old[4..6].copy_from_slice(&114u16.to_le_bytes());
    assert!(protocol::execution::process_request(&old).is_err());
    assert!(session.process_request(&old).is_err());
}
#[test]
fn report_context_engine_rank_statistics_and_completeness_tampering_are_refused() {
    let (query, subjects, options) = inputs();
    let report = native::search(&query, &subjects, &options, &mut || false).unwrap();
    let revision = report.context_revision.clone();
    let mut changed = report.clone();
    changed.engine_version = "sequence-similarity-1".into();
    assert!(changed.validate_binding(&revision).is_err());
    let mut changed = report.clone();
    changed.query_revision = "a".repeat(64);
    assert!(changed.validate_binding(&revision).is_err());
    let mut changed = report.clone();
    changed.complete = false;
    assert!(changed.validate_binding(&revision).is_err());
    let mut changed = report.clone();
    changed.hits[0].statistics.certain_matches += 1;
    assert!(changed.validate_binding(&revision).is_err());
    let mut changed = report.clone();
    changed.subjects[0].best.as_mut().unwrap().query_span.end += 1;
    assert!(changed.validate_binding(&revision).is_err());
    let mut changed = report.clone();
    changed.subjects[0].score_upper_bound_reached = false;
    changed.subjects[0].best_score_proven = true;
    assert!(changed.validate_binding(&revision).is_err());
    let mut changed = report;
    changed.work_used = options.max_work + 1;
    assert!(changed.validate_binding(&revision).is_err());
    let request = Request {
        id: 8,
        revision: "b".repeat(64),
        job: AnalysisJob::FindSimilarSequences {
            query,
            subjects,
            options,
        },
    };
    let bytes = protocol::request::Request::encode(&request).unwrap();
    assert!(protocol::execution::process_request(&bytes).is_err());
}
