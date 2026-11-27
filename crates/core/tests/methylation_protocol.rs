//! Append-only worker gates and exact source/options/profile assessment bindings.
use collagene_core::{
    Error, Molecule, Topology,
    methylation::*,
    protocol::{self, AnalysisJob, AnalysisResult, Request, Response},
};
fn fixture() -> (Molecule, MethylationOptions) {
    (
        Molecule::dna("assessment", "TCTAGATCGAATTCGGATCCCTCGAG", Topology::Linear).unwrap(),
        MethylationOptions::default(),
    )
}
fn raw(magic: &[u8; 4], version: u16, payload: &[u8]) -> Vec<u8> {
    let mut out = magic.to_vec();
    out.extend(version.to_le_bytes());
    out.extend((payload.len() as u32).to_le_bytes());
    out.extend(blake3::hash(payload).as_bytes());
    out.extend(payload);
    out
}
#[test]
fn worker_assessment_binds_source_states_profiles_and_options_without_new_product_authority() {
    let (molecule, options) = fixture();
    let revision = batch_revision(&molecule, &options).unwrap();
    let request = Request {
        id: 401,
        revision: revision.clone(),
        job: AnalysisJob::AssessMethylation {
            molecule: molecule.clone(),
            options: options.clone(),
        },
    };
    let bytes = protocol::encode_request(&request).unwrap();
    assert_eq!(
        &bytes[4..6],
        &collagene_core::protocol::WORKER_SCHEMA_VERSION.to_le_bytes()
    );
    let response = protocol::decode_response(&protocol::process_request(&bytes).unwrap()).unwrap();
    assert_eq!(response.revision, revision);
    let AnalysisResult::Methylation(report) = response.result else {
        panic!("wrong result")
    };
    report.validate_bindings_for(&molecule, &options).unwrap();
    replay(&molecule, &options, &report).unwrap();
    assert_eq!(report.profiles.len(), 12);
    assert!(!report.sites.is_empty());
    let mut changed = options.clone();
    changed.dam = MethylationState::FullyModified;
    assert_ne!(batch_revision(&molecule, &changed).unwrap(), revision);
    assert!(report.validate_bindings_for(&molecule, &changed).is_err());
    let wrong = Molecule::dna("other", "TCTAGATCGAATTCGGATCCCTCGAA", Topology::Linear).unwrap();
    assert!(report.validate_bindings_for(&wrong, &options).is_err());
    let mut forged = report.clone();
    forged.profile_revision = "0".repeat(64);
    assert!(forged.validate_bindings_for(&molecule, &options).is_err());
    let mut stale = request;
    stale.revision = "0".repeat(64);
    assert!(protocol::process_request(&protocol::encode_request(&stale).unwrap()).is_err());
    let response = Response {
        id: 401,
        revision: "0".repeat(64),
        result: AnalysisResult::Methylation(report),
    };
    assert!(
        protocol::decode_response(&raw(
            b"CGR\0",
            protocol::WORKER_SCHEMA_VERSION,
            &postcard::to_allocvec(&response).unwrap()
        ))
        .is_err()
    );
}
#[test]
fn methylation_tags_require_schema_sixteen_before_body_decode_and_keep_ordinary_caps() {
    for (magic, tag) in [(b"CGQ\0", 48u32), (b"CGR\0", 46u32)] {
        let prefix = postcard::to_allocvec(&(402u32, "assessment", tag)).unwrap();
        for version in 1..16 {
            let bytes = raw(magic, version, &prefix);
            let error = if magic == b"CGQ\0" {
                protocol::process_request(&bytes).unwrap_err()
            } else {
                protocol::decode_response(&bytes).unwrap_err()
            };
            assert!(
                matches!(error, collagene_core::Error::Unsupported(_)),
                "{error}"
            );
        }
        let mut large = prefix;
        large.resize(collagene_core::formats::MAX_IMPORT_BYTES + 1, 0);
        let bytes = raw(magic, protocol::WORKER_SCHEMA_VERSION, &large);
        let error = if magic == b"CGQ\0" {
            protocol::process_request(&bytes).unwrap_err()
        } else {
            protocol::decode_response(&bytes).unwrap_err()
        };
        assert!(matches!(error, Error::TooLarge));
    }
}
#[test]
fn bounded_source_evidence_precedes_revision_hashing() {
    let (mut molecule, options) = fixture();
    molecule
        .metadata
        .insert("oversized".into(), "x".repeat(9 * 1024 * 1024));
    let request = Request {
        id: 403,
        revision: "stale".into(),
        job: AnalysisJob::AssessMethylation { molecule, options },
    };
    assert!(matches!(
        protocol::process_request(&protocol::encode_request(&request).unwrap()),
        Err(Error::TooLarge)
    ));
}
