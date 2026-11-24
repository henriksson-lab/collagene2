use collagene_core::{Molecule, Topology, protocol::*, search::*, view_tracks::*};
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
fn new_jobs_roundtrip_bound_source_and_keep_old_search_tag() {
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
        let encoded = encode_request(&request).unwrap();
        assert_eq!(&encoded[4..6], &22u16.to_le_bytes());
        let response = decode_response(&process_request(&encoded).unwrap()).unwrap();
        assert_eq!(response.revision, m.revision());
        assert!(!matches!(response.result, AnalysisResult::Error(_)));
        for version in 1..9 {
            assert!(process_request(&wire(b"CGQ\0", version, &request)).is_err());
        }
    }
    let old = Request {
        id: 8,
        revision: m.revision(),
        job: AnalysisJob::Search {
            molecule: m.clone(),
            query: "AAAA".into(),
        },
    };
    for version in 1..=9 {
        let response =
            decode_response(&process_request(&wire(b"CGQ\0", version, &old)).unwrap()).unwrap();
        assert!(matches!(response.result, AnalysisResult::Search(_)));
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
        assert!(decode_response(&wire(b"CGR\0", version, &response)).is_err());
    }
    assert!(decode_response(&wire(b"CGR\0", 9, &response)).is_ok());
    let mut corrupt = wire(b"CGR\0", 9, &response);
    *corrupt.last_mut().unwrap() ^= 1;
    assert!(decode_response(&corrupt).is_err());
    let mut invalid = report;
    invalid.hits[0].parts[0].end = m.sequence.len() + 1;
    let response = Response {
        id: 2,
        revision: m.revision(),
        result: AnalysisResult::DocumentSearch(invalid),
    };
    assert!(decode_response(&wire(b"CGR\0", 9, &response)).is_err());
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
        assert!(decode_response(&wire(b"CGR\0", version, &response)).is_err());
    }
    assert!(decode_response(&wire(b"CGR\0", 9, &response)).is_ok());
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
    assert!(process_request(&encode_request(&r).unwrap()).is_err());
    let response = Response {
        id: 0,
        revision: "stale".into(),
        result: AnalysisResult::DocumentSearch(search_document(&m, &options()).unwrap()),
    };
    assert!(decode_response(&wire(b"CGR\0", 9, &response)).is_err());
}
#[test]
fn oversized_new_variant_is_rejected_before_payload_deserialization() {
    let mut payload = postcard::to_allocvec(&(0u32, "", 30u32)).unwrap();
    payload.resize(32 * 1024 * 1024 + 1, 0);
    let mut b = b"CGQ\0".to_vec();
    b.extend(9u16.to_le_bytes());
    b.extend((payload.len() as u32).to_le_bytes());
    b.extend(blake3::hash(&payload).as_bytes());
    b.extend(payload);
    assert!(matches!(
        process_request(&b),
        Err(collagene_core::Error::TooLarge)
    ));
}
#[test]
fn retained_session_rejects_missing_stale_and_replaced_source() {
    let a = source();
    let b = Molecule::dna("other", "ACGTACGT", Topology::Linear).unwrap();
    let mut session = AnalysisSession::default();
    let ordinary = encode_request(&Request {
        id: 77,
        revision: a.revision(),
        job: AnalysisJob::Search {
            molecule: a.clone(),
            query: "AAAA".into(),
        },
    })
    .unwrap();
    assert!(matches!(
        decode_response(&session.process_request(&ordinary).unwrap())
            .unwrap()
            .result,
        AnalysisResult::Search(_)
    ));

    let cached = |revision: String| {
        encode_request(&Request {
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
        encode_request(&Request {
            id: 1,
            revision: m.revision(),
            job: AnalysisJob::PrepareTracks { molecule: m },
        })
        .unwrap()
    };
    let r = decode_response(&session.process_request(&prep(a.clone())).unwrap()).unwrap();
    assert!(matches!(r.result, AnalysisResult::TracksPrepared(_)));
    assert!(matches!(
        decode_response(&session.process_request(&cached(a.revision())).unwrap())
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
            .process_request(&encode_request(&malformed).unwrap())
            .is_err()
    );
    assert!(session.process_request(&cached(b.revision())).is_ok());
    let mut old = prep(b);
    old[4..6].copy_from_slice(&8u16.to_le_bytes());
    assert!(session.process_request(&old).is_err());
}
