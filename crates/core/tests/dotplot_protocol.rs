//! New dot-plot selectors are ordinary bounded scientific messages, not workspaces.
use collagene_core::{
    Molecule, Topology,
    dotplot::*,
    protocol::{self, AnalysisJob, AnalysisResult, Request, Response},
};
fn fixture() -> (Molecule, Molecule, DotPlotOptions) {
    (
        Molecule::dna("reference", "ACGTTT", Topology::Linear).unwrap(),
        Molecule::dna("query", "ACG", Topology::Linear).unwrap(),
        DotPlotOptions {
            window: 3,
            min_matches: 3,
            ..Default::default()
        },
    )
}
fn envelope(magic: &[u8; 4], version: u16, payload: &[u8]) -> Vec<u8> {
    let mut bytes = magic.to_vec();
    bytes.extend(version.to_le_bytes());
    bytes.extend((payload.len() as u32).to_le_bytes());
    bytes.extend(blake3::hash(payload).as_bytes());
    bytes.extend(payload);
    bytes
}
fn raw(magic: &[u8; 4], version: u16, value: &impl serde::Serialize) -> Vec<u8> {
    envelope(magic, version, &postcard::to_allocvec(value).unwrap())
}
#[test]
fn dotplot_native_worker_roundtrip_binds_both_sources_options_and_original_reverse_coordinates() {
    let (reference, query, options) = fixture();
    let req = Request {
        id: 91,
        revision: reference.revision(),
        job: AnalysisJob::DotPlot {
            molecule: reference.clone(),
            query: query.clone(),
            options: options.clone(),
        },
    };
    let bytes = protocol::encode_request(&req).unwrap();
    assert_eq!(&bytes[4..6], &21u16.to_le_bytes());
    let response = protocol::decode_response(&protocol::process_request(&bytes).unwrap()).unwrap();
    assert_eq!(response.id, 91);
    assert_eq!(response.revision, req.revision);
    let AnalysisResult::DotPlot(report) = response.result else {
        panic!("wrong result")
    };
    report
        .validate_bindings_for(&reference, &query, &options)
        .unwrap();
    assert!(report.complete);
    assert!(!report.sampled);
    assert!(report.hits.iter().any(|h| h.reference_start == 0
        && h.query_start == 0
        && !h.reverse
        && h.certain_matches == 3));
    assert!(report.hits.iter().any(|h| h.reference_start == 1
        && h.query_start == 2
        && h.reverse
        && h.certain_matches == 3));
    let wrong = Molecule::dna("otherquery", "AAA", Topology::Linear).unwrap();
    assert!(
        report
            .validate_bindings_for(&reference, &wrong, &options)
            .is_err()
    );
    let mut changed = options.clone();
    changed.stride = 2;
    assert!(
        report
            .validate_bindings_for(&reference, &query, &changed)
            .is_err()
    );
    let mismatched = Response {
        id: 91,
        revision: "0".repeat(64),
        result: AnalysisResult::DotPlot(report),
    };
    assert!(protocol::decode_response(&raw(b"CGR\0", 14, &mismatched)).is_err());
    let stale = Request {
        revision: "0".repeat(64),
        ..req
    };
    assert!(protocol::process_request(&protocol::encode_request(&stale).unwrap()).is_err());
}
#[test]
fn dotplot_tag_floors_apply_before_decoding_and_never_receive_the_workspace_message_budget() {
    let (reference, query, options) = fixture();
    let req = Request {
        id: 91,
        revision: reference.revision(),
        job: AnalysisJob::DotPlot {
            molecule: reference.clone(),
            query: query.clone(),
            options: options.clone(),
        },
    };
    let report = analyze(&reference, &query, &options).unwrap();
    let res = Response {
        id: 91,
        revision: reference.revision(),
        result: AnalysisResult::DotPlot(report),
    };
    for old in 1..14 {
        assert!(protocol::process_request(&raw(b"CGQ\0", old, &req)).is_err());
        assert!(protocol::decode_response(&raw(b"CGR\0", old, &res)).is_err());
    }
    for (magic, tag) in [(b"CGQ\0", 46u32), (b"CGR\0", 44u32)] {
        let prefix = postcard::to_allocvec(&(91u32, "reference", tag)).unwrap();
        let bytes = envelope(magic, 13, &prefix);
        let err = if magic == b"CGQ\0" {
            protocol::process_request(&bytes).unwrap_err()
        } else {
            protocol::decode_response(&bytes).unwrap_err()
        };
        assert!(err.to_string().contains("schema 14"), "{err}");
        let mut oversized = prefix;
        oversized.resize(32 * 1024 * 1024 + 1, 0);
        let bytes = envelope(magic, 14, &oversized);
        let err = if magic == b"CGQ\0" {
            protocol::process_request(&bytes).unwrap_err()
        } else {
            protocol::decode_response(&bytes).unwrap_err()
        };
        assert!(matches!(err, collagene_core::Error::TooLarge));
    }
}
#[test]
fn source_and_query_evidence_preflight_precedes_revision_hashing() {
    let (reference, query, options) = fixture();
    for query_oversized in [false, true] {
        let mut r = reference.clone();
        let mut q = query.clone();
        if query_oversized {
            q.metadata
                .insert("oversized".into(), "x".repeat(32 * 1024 * 1024));
        } else {
            r.metadata
                .insert("oversized".into(), "x".repeat(32 * 1024 * 1024));
        }
        // Direct helper avoids a premature giant hash even when no wire is involved.
        assert!(matches!(
            source_revisions(&r, &q, &options),
            Err(collagene_core::Error::TooLarge)
        ));
    }
}
