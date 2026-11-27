use collagene_core::{
    Molecule, Topology, mutation_intent,
    protocol::{self, AnalysisJob, AnalysisResult, Request},
};
fn request() -> Request {
    let molecule = Molecule::dna("original", "ACGTCAGT", Topology::Circular).unwrap();
    let request = mutation_intent::Request {
        name: "declared mutant".into(),
        edit_start: 6,
        edit_end: 10,
        inserted: "TT".into(),
    };
    Request {
        id: 410,
        revision: mutation_intent::intent_revision(&molecule, &request).unwrap(),
        job: AnalysisJob::BuildMutationIntent { molecule, request },
    }
}
#[test]
fn source_edit_roundtrip_has_literal_target_and_rejects_pre_41_envelopes() {
    let req = request();
    let bytes = protocol::encode_request(&req).unwrap();
    let ((_, _, tag), _) = postcard::take_from_bytes::<(u32, &str, u32)>(&bytes[42..]).unwrap();
    assert_eq!(tag, 102);
    let mut old = bytes.clone();
    old[4..6].copy_from_slice(&40u16.to_le_bytes());
    assert!(protocol::process_request(&old).is_err());
    let bytes = protocol::process_request(&bytes).unwrap();
    let ((_, _, tag), _) = postcard::take_from_bytes::<(u32, &str, u32)>(&bytes[42..]).unwrap();
    assert_eq!(tag, 94);
    let response = protocol::decode_response(&bytes).unwrap();
    let AnalysisResult::MutationIntent(report) = response.result else {
        panic!("expected declared reference")
    };
    report.validate_binding(&req.revision).unwrap();
    assert_eq!(report.target.sequence.as_str(), "TTGTCA");
    assert_eq!(report.deleted, "GTAC");
    assert_eq!(report.origin_shift, 6);
    let mut old = bytes;
    old[4..6].copy_from_slice(&40u16.to_le_bytes());
    assert!(protocol::decode_response(&old).is_err());
}
#[test]
fn an_edit_or_source_change_cannot_reuse_the_previous_context() {
    let mut req = request();
    if let AnalysisJob::BuildMutationIntent { request, .. } = &mut req.job {
        request.inserted = "A".into();
    }
    assert!(protocol::process_request(&protocol::encode_request(&req).unwrap()).is_err());
    let mut req = request();
    if let AnalysisJob::BuildMutationIntent { molecule, .. } = &mut req.job {
        molecule.name = "another source".into();
    }
    assert!(protocol::process_request(&protocol::encode_request(&req).unwrap()).is_err());
}
