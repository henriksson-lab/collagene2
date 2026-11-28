//! The physical-template semantic fix admits only the current report contract.
use collagene_core::{
    planning::plan::{DesignPlan, Operation, PlanInput, Source},
    protocol::{
        self,
        request::Request,
        response::Response,
        schema::{AnalysisJob, AnalysisResult},
    },
    sequence::analysis::Primer,
    verification::{self, evidence::Status},
};
use collagene_model::{Molecule, Topology};
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
fn current_encoding_roundtrips_scientific_and_documentary_receipts() {
    let original =
        Molecule::dna("Original reference", "ACGATTCGAACTGCTA", Topology::Linear).unwrap();
    let target = Molecule::dna(
        "Independent target",
        original.sequence.as_str(),
        Topology::Linear,
    )
    .unwrap();
    let plan = DesignPlan {
        standard: None,
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
            template_strand: collagene_core::chemistry::Strand::Top,
        }],
    };
    let request = Request {
        id: 80,
        revision: target.revision(),
        job: AnalysisJob::Verify { target, plan },
    };
    let encoded = protocol::execution::process_request(
        &protocol::request::Request::encode(&request).unwrap(),
    )
    .unwrap();
    assert_eq!(
        &encoded[4..6],
        &protocol::schema::WORKER_SCHEMA_VERSION.to_le_bytes()
    );
    let response = protocol::response::Response::decode(&encoded).unwrap();
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
    assert!(
        report
            .engine_version
            .eq(collagene_core::verification::ENGINE_VERSION)
    );
    assert!(
        protocol::response::Response::decode(&wire(
            protocol::schema::WORKER_SCHEMA_VERSION,
            &response
        ))
        .is_ok()
    );
    assert!(protocol::response::Response::decode(&wire(92, &response)).is_err());
    assert!(protocol::response::Response::decode(&wire(91, &response)).is_err());
    assert!(protocol::response::Response::decode(&wire(90, &response)).is_err());
    assert!(protocol::response::Response::decode(&wire(89, &response)).is_err());
    assert!(protocol::response::Response::decode(&wire(88, &response)).is_err());
    assert!(protocol::response::Response::decode(&wire(87, &response)).is_err());
    assert!(protocol::response::Response::decode(&wire(81, &response)).is_err());
    assert!(protocol::response::Response::decode(&wire(80, &response)).is_err());
    assert!(protocol::response::Response::decode(&wire(79, &response)).is_err());
    assert!(protocol::response::Response::decode(&wire(73, &response)).is_err());

    // A superseded engine stamp is refused even inside a current frame.
    let mut historical = response;
    let AnalysisResult::Verification(report) = &mut historical.result else {
        unreachable!()
    };
    let current = verification::ENGINE_VERSION;
    let old = current.replace(collagene_core::verification::ENGINE_VERSION, "/science-20");
    report.engine_version = old.clone();
    for step in &mut report.steps {
        let product = step
            .product
            .as_mut()
            .expect("each declared step must reconstruct");
        for value in product.metadata.values_mut() {
            if value == current {
                *value = old.clone();
            }
        }
    }
    assert!(
        protocol::response::Response::decode(&wire(
            protocol::schema::WORKER_SCHEMA_VERSION,
            &historical,
        ))
        .is_err()
    );
    assert!(protocol::response::Response::decode(&wire(72, &historical)).is_err());
}
