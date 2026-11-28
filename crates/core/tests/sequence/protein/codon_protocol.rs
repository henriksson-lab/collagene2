//! Codon design is a separate source-bound model, never a physical cloning operation.
use collagene_core::{
    formats::sequence as formats,
    protocol::{
        self,
        request::Request,
        response::Response,
        schema::{AnalysisJob, AnalysisResult},
    },
    sequence::codon_design::*,
};
use collagene_model::{Molecule, Topology};
fn source() -> Molecule {
    formats::import_text(">required [molecule=protein]\nMK*\n")
        .unwrap()
        .molecules
        .remove(0)
}
fn envelope(magic: &[u8; 4], version: u16, value: &impl serde::Serialize) -> Vec<u8> {
    raw(magic, version, &postcard::to_allocvec(value).unwrap())
}
fn raw(magic: &[u8; 4], version: u16, payload: &[u8]) -> Vec<u8> {
    let mut out = magic.to_vec();
    out.extend(version.to_le_bytes());
    out.extend((payload.len() as u32).to_le_bytes());
    out.extend(blake3::hash(payload).as_bytes());
    out.extend(payload);
    out
}
fn roundtrip(request: &Request) -> Response {
    let bytes = protocol::request::Request::encode(request).unwrap();
    assert_eq!(
        &bytes[4..6],
        &collagene_core::protocol::schema::WORKER_SCHEMA_VERSION.to_le_bytes()
    );
    let response = protocol::response::Response::decode(
        &protocol::execution::process_request(&bytes).unwrap(),
    )
    .unwrap();
    assert_eq!(response.id, request.id);
    assert_eq!(response.revision, request.revision);
    response
}
#[test]
fn design_and_independent_check_are_source_options_and_candidate_bound() {
    let source = source();
    let options = DesignOptions {
        preferences: collagene_core::sequence::codon_design::CodonPreferences::Explicit {
            weights: vec![CodonWeight {
                codon: "AAG".into(),
                weight: 100,
            }],
        },
        ..Default::default()
    };
    let request = Request {
        id: 91,
        revision: source.revision(),
        job: AnalysisJob::DesignCodons {
            molecule: source.clone(),
            options: options.clone(),
        },
    };
    let result = roundtrip(&request);
    let AnalysisResult::CodonDesign(report) = result.result else {
        panic!("wrong design result")
    };
    report.validate_for_source(&source, &options).unwrap();
    let candidate = report.candidate.clone().unwrap();
    assert_eq!(candidate.sequence.as_str(), "ATGAAGTAA");
    let checked = roundtrip(&Request {
        id: 92,
        revision: source.revision(),
        job: AnalysisJob::CheckCodonDesign {
            source: source.clone(),
            candidate: candidate.clone(),
            options: options.clone(),
        },
    });
    let AnalysisResult::CodonCheck(check) = checked.result else {
        panic!("wrong check result")
    };
    check.validate_for(&source, &candidate, &options).unwrap();
    assert!(check.complete);
    assert_eq!(check.satisfies, Some(true));
    let wrong = Molecule::dna("wrong", "ATGGAGTAA", Topology::Linear).unwrap();
    let rejected = roundtrip(&Request {
        id: 93,
        revision: source.revision(),
        job: AnalysisJob::CheckCodonDesign {
            source: source.clone(),
            candidate: wrong,
            options: options.clone(),
        },
    });
    let AnalysisResult::CodonCheck(failed) = rejected.result else {
        panic!()
    };
    assert_eq!(failed.satisfies, Some(false));
    let mut changed = options;
    changed.genetic_code = 2;
    assert!(report.validate_for_source(&source, &changed).is_err());
    let mut stale = request;
    stale.revision = "0".repeat(64);
    assert!(
        protocol::execution::process_request(&protocol::request::Request::encode(&stale).unwrap())
            .is_err()
    );
    let mut response = Response {
        id: 91,
        revision: "0".repeat(64),
        result: AnalysisResult::CodonDesign(report),
    };
    assert!(
        protocol::response::Response::decode(&envelope(
            b"CGR\0",
            protocol::schema::WORKER_SCHEMA_VERSION,
            &response
        ))
        .is_err()
    );
    response.result = AnalysisResult::CodonCheck(check);
    assert!(
        protocol::response::Response::decode(&envelope(
            b"CGR\0",
            protocol::schema::WORKER_SCHEMA_VERSION,
            &response
        ))
        .is_err()
    );
}
#[test]
fn new_tags_reject_old_versions_before_payload_and_retain_ordinary_limits() {
    let source = source();
    let options = DesignOptions::default();
    for job in [
        AnalysisJob::DesignCodons {
            molecule: source.clone(),
            options: options.clone(),
        },
        AnalysisJob::CheckCodonDesign {
            source: source.clone(),
            candidate: Molecule::dna("candidate", "ATGAAATAA", Topology::Linear).unwrap(),
            options: options.clone(),
        },
    ] {
        let req = Request {
            id: 1,
            revision: source.revision(),
            job,
        };
        for old in 1..13 {
            assert!(protocol::execution::process_request(&envelope(b"CGQ\0", old, &req)).is_err());
        }
    }
    let design = design(&source, &options).unwrap();
    let candidate = design.candidate.clone().unwrap();
    let check = check(&source, &candidate, &options).unwrap();
    for result in [
        AnalysisResult::CodonDesign(design),
        AnalysisResult::CodonCheck(check),
    ] {
        let response = Response {
            id: 1,
            revision: source.revision(),
            result,
        };
        for old in 1..13 {
            assert!(
                protocol::response::Response::decode(&envelope(b"CGR\0", old, &response)).is_err()
            );
        }
    }
    // Independent literal prefix: malformed missing variant bodies must hit version
    // gates without attempting to deserialize the new schema under an older label.
    for tag in [44u32, 45] {
        let prefix = postcard::to_allocvec(&(1u32, "source", tag)).unwrap();
        let err = protocol::execution::process_request(&raw(b"CGQ\0", 12, &prefix)).unwrap_err();
        assert!(matches!(err, collagene_model::Error::Unsupported(_)));
    }
    for tag in [42u32, 43] {
        let prefix = postcard::to_allocvec(&(1u32, "source", tag)).unwrap();
        let err = protocol::response::Response::decode(&raw(b"CGR\0", 12, &prefix)).unwrap_err();
        assert!(matches!(err, collagene_model::Error::Unsupported(_)));
    }
    let mut payload = postcard::to_allocvec(&(1u32, "source", 44u32)).unwrap();
    payload.resize(32 * 1024 * 1024 + 1, 0);
    assert!(protocol::execution::process_request(&raw(b"CGQ\0", 13, &payload)).is_err());
    let mut payload = postcard::to_allocvec(&(1u32, "source", 42u32)).unwrap();
    payload.resize(32 * 1024 * 1024 + 1, 0);
    assert!(protocol::response::Response::decode(&raw(b"CGR\0", 13, &payload)).is_err());
}

#[test]
fn codon_source_evidence_is_bounded_before_revision_hashing() {
    let mut source = source();
    source.metadata.insert(
        "source:oversized_evidence".into(),
        "x".repeat(4 * 1024 * 1024 + 1),
    );
    let options = DesignOptions::default();
    for job in [
        AnalysisJob::DesignCodons {
            molecule: source.clone(),
            options: options.clone(),
        },
        AnalysisJob::CheckCodonDesign {
            source: source.clone(),
            candidate: Molecule::dna("candidate", "ATGAAATAA", Topology::Linear).unwrap(),
            options: options.clone(),
        },
    ] {
        // Deliberately stale revision distinguishes preflight from a revision
        // comparison that needlessly hashes the oversized retained evidence.
        let request = Request {
            id: 94,
            revision: "0".repeat(64),
            job,
        };
        let bytes = protocol::request::Request::encode(&request).unwrap();
        assert!(matches!(
            protocol::execution::process_request(&bytes),
            Err(collagene_model::Error::TooLarge)
        ));
    }
}
