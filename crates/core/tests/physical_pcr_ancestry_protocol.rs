//! The physical-template semantic fix is versioned; old reports stay documentary.
use collagene_core::{
    Molecule, Topology,
    analysis::Primer,
    protocol::{self, AnalysisJob, AnalysisResult, Request, Response},
    verification::{self, DesignPlan, Operation, PlanInput, Source, Status},
};
fn wire(version: u16, response: &Response) -> Vec<u8> {
    let payload = postcard::to_allocvec(response).unwrap();
    let mut bytes = b"CGR\0".to_vec();
    bytes.extend(version.to_le_bytes());
    bytes.extend((payload.len() as u32).to_le_bytes());
    bytes.extend(blake3::hash(&payload).as_bytes());
    bytes.extend(payload);
    bytes
}
#[test]
fn current_scientific_receipts_use_current_worker_and_historical20_stays_readable() {
    let original =
        Molecule::dna("Original reference", "ACGATTCGAACTGCTA", Topology::Linear).unwrap();
    let target = Molecule::dna(
        "Independent target",
        original.sequence.as_str(),
        Topology::Linear,
    )
    .unwrap();
    let plan = DesignPlan {
        name: "Original template hypothesis PCR".into(),
        target_revision: target.revision(),
        inputs: vec![PlanInput {
            id: "original".into(),
            revision: original.revision(),
            molecule: original,
        }],
        steps: vec![Operation::Pcr {
            template: Source::Input("original".into()),
            forward: Primer::new("Forward", "ACGATTCG", "").unwrap(),
            reverse: Primer::new("Reverse", "TAGCAGTT", "").unwrap(),
        }],
    };
    let request = Request {
        id: 80,
        revision: target.revision(),
        job: AnalysisJob::Verify { target, plan },
    };
    let encoded = protocol::process_request(&protocol::encode_request(&request).unwrap()).unwrap();
    assert_eq!(
        &encoded[4..6],
        &protocol::WORKER_SCHEMA_VERSION.to_le_bytes()
    );
    let response = protocol::decode_response(&encoded).unwrap();
    let AnalysisResult::Verification(report) = &response.result else {
        panic!("verification reply");
    };
    assert_eq!(
        report.status,
        Status::ConsistentUnderAssumptions,
        "{:?}",
        report.findings
    );
    assert_eq!(report.engine_version, verification::ENGINE_VERSION);
    assert!(report.engine_version.ends_with("/science-25"));
    assert!(protocol::decode_response(&wire(91, &response)).is_ok());
    assert!(protocol::decode_response(&wire(90, &response)).is_err());
    assert!(protocol::decode_response(&wire(89, &response)).is_err());
    assert!(protocol::decode_response(&wire(88, &response)).is_err());
    assert!(protocol::decode_response(&wire(87, &response)).is_err());
    assert!(protocol::decode_response(&wire(81, &response)).is_err());
    assert!(protocol::decode_response(&wire(80, &response)).is_err());
    assert!(protocol::decode_response(&wire(79, &response)).is_err());
    assert!(protocol::decode_response(&wire(73, &response)).is_err());

    // Synthetic documentary fixture of the preceding unchanged original-input
    // PCR model. Lowering its stamp never makes it a current-engine receipt.
    let mut historical = response;
    let AnalysisResult::Verification(report) = &mut historical.result else {
        unreachable!()
    };
    let current = verification::ENGINE_VERSION;
    let old = current.replace("/science-25", "/science-20");
    report.engine_version = old.clone();
    for product in &mut report.products {
        for value in product.metadata.values_mut() {
            if value == current {
                *value = old.clone();
            }
        }
    }
    let decoded = protocol::decode_response(&wire(79, &historical)).unwrap();
    let AnalysisResult::Verification(report) = decoded.result else {
        unreachable!()
    };
    assert!(report.engine_version.ends_with("/science-20"));
    assert_ne!(report.engine_version, verification::ENGINE_VERSION);
    assert!(protocol::decode_response(&wire(72, &historical)).is_err());
}
