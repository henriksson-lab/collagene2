use collagene_core::{
    protocol::{
        execution::process_request,
        request::Request,
        response::Response,
        schema::{AnalysisJob, AnalysisResult, WORKER_SCHEMA_VERSION},
        session::AnalysisSession,
    },
    sequence::search::*,
    tracks::view::*,
};
use collagene_model::{Molecule, Topology};
fn source() -> Molecule {
    Molecule::dna("source", "AAAAGAATTCCCC", Topology::Circular).unwrap()
}
fn options() -> SearchOptions {
    SearchOptions {
        query: "EcoRI".into(),
        mode: SearchMode::Enzymes,
        ..Default::default()
    }
}
fn wire(magic: &[u8; 4], version: u16, value: &impl serde::Serialize) -> Vec<u8> {
    let payload = postcard::to_allocvec(value).unwrap();
    let mut b = magic.to_vec();
    b.extend(version.to_le_bytes());
    b.extend((payload.len() as u32).to_le_bytes());
    b.extend(blake3::hash(&payload).as_bytes());
    b.extend(payload);
    b
}
#[test]
fn search_and_track_jobs_roundtrip_bound_sources() {
    let m = source();
    for job in [
        AnalysisJob::SearchDocument {
            molecule: m.clone(),
            options: options(),
        },
        AnalysisJob::ViewTracks {
            molecule: m.clone(),
            window: TrackWindow {
                start: 0,
                end: m.sequence.len(),
            },
            options: TrackOptions::default(),
        },
    ] {
        let request = Request {
            id: 7,
            revision: m.revision(),
            job,
        };
        let encoded = collagene_core::protocol::request::Request::encode(&request).unwrap();
        assert_eq!(
            &encoded[4..6],
            &collagene_core::protocol::schema::WORKER_SCHEMA_VERSION.to_le_bytes()
        );
        let response = collagene_core::protocol::response::Response::decode(
            &process_request(&encoded).unwrap(),
        )
        .unwrap();
        assert_eq!(response.revision, m.revision());
        assert!(!matches!(response.result, AnalysisResult::Error(_)));
        for version in 1..9 {
            assert!(process_request(&wire(b"CGQ\0", version, &request)).is_err());
        }
    }
}

#[test]
fn new_results_cannot_claim_old_versions_and_corruption_fails() {
    let m = source();
    let report = search_document(&m, &options()).unwrap();
    let response = Response {
        id: 2,
        revision: m.revision(),
        result: AnalysisResult::DocumentSearch(report.clone()),
    };
    for version in 1..9 {
        assert!(
            collagene_core::protocol::response::Response::decode(&wire(
                b"CGR\0", version, &response
            ))
            .is_err()
        );
    }
    assert!(
        collagene_core::protocol::response::Response::decode(&wire(
            b"CGR\0",
            collagene_core::protocol::schema::WORKER_SCHEMA_VERSION,
            &response
        ))
        .is_ok()
    );
    let mut corrupt = wire(
        b"CGR\0",
        collagene_core::protocol::schema::WORKER_SCHEMA_VERSION,
        &response,
    );
    *corrupt.last_mut().unwrap() ^= 1;
    assert!(collagene_core::protocol::response::Response::decode(&corrupt).is_err());
    let mut invalid = report;
    invalid.hits[0].parts[0].end = m.sequence.len() + 1;
    let response = Response {
        id: 2,
        revision: m.revision(),
        result: AnalysisResult::DocumentSearch(invalid),
    };
    assert!(
        collagene_core::protocol::response::Response::decode(&wire(
            b"CGR\0",
            collagene_core::protocol::schema::WORKER_SCHEMA_VERSION,
            &response
        ))
        .is_err()
    );
    let report = view_tracks(
        &m,
        TrackWindow { start: 0, end: 10 },
        &TrackOptions::default(),
    )
    .unwrap();
    let response = Response {
        id: 2,
        revision: m.revision(),
        result: AnalysisResult::ViewTracks(report),
    };
    for version in 1..9 {
        assert!(
            collagene_core::protocol::response::Response::decode(&wire(
                b"CGR\0", version, &response
            ))
            .is_err()
        );
    }
    assert!(
        collagene_core::protocol::response::Response::decode(&wire(
            b"CGR\0",
            collagene_core::protocol::schema::WORKER_SCHEMA_VERSION,
            &response
        ))
        .is_ok()
    );
}
#[test]
fn stale_input_and_result_revisions_are_rejected() {
    let m = source();
    let r = Request {
        id: 0,
        revision: "stale".into(),
        job: AnalysisJob::SearchDocument {
            molecule: m.clone(),
            options: options(),
        },
    };
    assert!(
        process_request(&collagene_core::protocol::request::Request::encode(&r).unwrap()).is_err()
    );
    let response = Response {
        id: 0,
        revision: "stale".into(),
        result: AnalysisResult::DocumentSearch(search_document(&m, &options()).unwrap()),
    };
    assert!(
        collagene_core::protocol::response::Response::decode(&wire(
            b"CGR\0",
            WORKER_SCHEMA_VERSION,
            &response
        ))
        .is_err()
    );
}
#[test]
fn oversized_new_variant_is_rejected_before_payload_deserialization() {
    let mut payload = postcard::to_allocvec(&(0u32, "", 30u32)).unwrap();
    payload.resize(32 * 1024 * 1024 + 1, 0);
    let mut b = b"CGQ\0".to_vec();
    b.extend(collagene_core::protocol::schema::WORKER_SCHEMA_VERSION.to_le_bytes());
    b.extend((payload.len() as u32).to_le_bytes());
    b.extend(blake3::hash(&payload).as_bytes());
    b.extend(payload);
    assert!(matches!(
        process_request(&b),
        Err(collagene_model::Error::TooLarge)
    ));
}
#[test]
fn retained_session_rejects_missing_stale_and_replaced_source() {
    let a = source();
    let b = Molecule::dna("other", "ACGTACGT", Topology::Linear).unwrap();
    let mut session = AnalysisSession::default();
    let ordinary = collagene_core::protocol::request::Request::encode(&Request {
        id: 77,
        revision: a.revision(),
        job: AnalysisJob::Search {
            molecule: a.clone(),
            query: "AAAA".into(),
        },
    })
    .unwrap();
    assert!(matches!(
        collagene_core::protocol::response::Response::decode(
            &session.process_request(&ordinary).unwrap()
        )
        .unwrap()
        .result,
        AnalysisResult::Search(_)
    ));

    let cached = |revision: String| {
        collagene_core::protocol::request::Request::encode(&Request {
            id: 4,
            revision,
            job: AnalysisJob::CachedViewTracks {
                window: TrackWindow { start: 0, end: 8 },
                options: TrackOptions::default(),
            },
        })
        .unwrap()
    };
    assert!(session.process_request(&cached(a.revision())).is_err());
    let prep = |m: Molecule| {
        collagene_core::protocol::request::Request::encode(&Request {
            id: 1,
            revision: m.revision(),
            job: AnalysisJob::PrepareTracks { molecule: m },
        })
        .unwrap()
    };
    let r = collagene_core::protocol::response::Response::decode(
        &session.process_request(&prep(a.clone())).unwrap(),
    )
    .unwrap();
    assert!(matches!(r.result, AnalysisResult::TracksPrepared(_)));
    assert!(matches!(
        collagene_core::protocol::response::Response::decode(
            &session.process_request(&cached(a.revision())).unwrap()
        )
        .unwrap()
        .result,
        AnalysisResult::ViewTracks(_)
    ));
    assert!(session.process_request(&cached("stale".into())).is_err());
    session.process_request(&prep(b.clone())).unwrap();
    assert!(session.process_request(&cached(a.revision())).is_err());
    assert!(session.process_request(&cached(b.revision())).is_ok());
    let malformed = Request {
        id: 1,
        revision: "stale".into(),
        job: AnalysisJob::PrepareTracks { molecule: a },
    };
    assert!(
        session
            .process_request(
                &collagene_core::protocol::request::Request::encode(&malformed).unwrap()
            )
            .is_err()
    );
    assert!(session.process_request(&cached(b.revision())).is_ok());
    let mut old = prep(b);
    old[4..6].copy_from_slice(&8u16.to_le_bytes());
    assert!(session.process_request(&old).is_err());
}
