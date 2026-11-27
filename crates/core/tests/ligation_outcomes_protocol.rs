//! Current ligation outcomes retain their intended path, alternatives and geometry.
use collagene_core::{
    Molecule, Topology,
    chemistry::Terminus,
    protocol::{self, AnalysisJob, AnalysisResult},
    transfer,
    verification::{
        self, DesignPlan, Operation, OrientedSource, PhysicalProcessingEvidence, PlanInput, Source,
        Status,
    },
    workspace,
};

fn fixture() -> (DesignPlan, Molecule) {
    // The intended product is independently written, never selected from the
    // generated alternatives. Reversing the second operand instead makes GGGG.
    let target = Molecule::dna("Independent intended DNA", "AAAACCCC", Topology::Linear).unwrap();
    let inputs = [("a", "AAAA"), ("b", "CCCC")]
        .into_iter()
        .map(|(id, text)| {
            let molecule = Molecule::dna(id, text, Topology::Linear).unwrap();
            PlanInput {
                id: id.into(),
                revision: molecule.revision(),
                molecule,
            }
        })
        .collect();
    let mut steps = ["a", "b"]
        .into_iter()
        .map(|id| Operation::DeclareSourceTermini {
            template: Source::Input(id.into()),
            left: Terminus::CLEAVED,
            right: Terminus::CLEAVED,
        })
        .collect::<Vec<_>>();
    steps.push(Operation::Ligate {
        fragments: (0..2)
            .map(|i| OrientedSource {
                source: Source::Step(i),
                reverse: false,
            })
            .collect(),
        circular: false,
        allow_nicked: false,
        options: Default::default(),
    });
    (
        DesignPlan {
            name: "Intended path plus alternatives".into(),
            target_revision: target.revision(),
            inputs,
            steps,
        },
        target,
    )
}

fn tag(value: &impl serde::Serialize) -> u32 {
    postcard::take_from_bytes::<u32>(&postcard::to_allocvec(value).unwrap())
        .unwrap()
        .0
}

#[test]
fn ligation_and_evidence_keep_intended_and_extractable_alternatives() {
    let (plan, target) = fixture();
    assert_eq!(tag(&plan.steps[2]), 4);
    let report = verification::verify(&plan, &target);
    assert_eq!(
        report.status,
        Status::ConsistentUnderAssumptions,
        "{:?}",
        report.findings
    );
    assert_eq!(report.products[2].sequence.as_str(), "AAAACCCC");
    let evidence = &report.physical_processing_reports[0].evidence;
    assert_eq!(tag(evidence), 7);
    let PhysicalProcessingEvidence::LigationOutcomes(evidence) = evidence else {
        panic!("outcomes")
    };
    assert_eq!(evidence.report.intended.sequence.as_str(), "AAAACCCC");
    assert!(
        evidence
            .report
            .alternatives
            .iter()
            .any(|p| p.sequence.as_str() == "AAAAGGGG")
    );
    evidence.validate_documentary().unwrap();
    let serialized = serde_json::to_value(evidence).unwrap();
    assert!(
        serialized
            .get("report")
            .unwrap()
            .get("alternatives")
            .is_some()
    );
    assert!(!report.exhaustive_alternative_search);
}

#[test]
fn current_worker_and_candidate_roundtrips_validate_full_outcome_geometry() {
    let (plan, target) = fixture();
    let candidate = transfer::Candidate::new(plan.clone());
    let bytes = transfer::save_candidate(&candidate).unwrap();
    assert_eq!(
        &bytes[4..6],
        &protocol::CANDIDATE_SCHEMA_VERSION.to_le_bytes()
    );
    assert_eq!(
        transfer::load_candidate(&bytes).unwrap().plan_revision,
        candidate.plan_revision
    );
    assert_eq!(workspace::PLAN_SCHEMA_VERSION, 27);
    let request = protocol::Request {
        id: 4,
        revision: target.revision(),
        job: AnalysisJob::Verify { plan, target },
    };
    let bytes = protocol::encode_request(&request).unwrap();
    let bytes = protocol::process_request(&bytes).unwrap();
    let mut response = protocol::decode_response(&bytes).unwrap();
    let AnalysisResult::Verification(report) = &mut response.result else {
        panic!("verification")
    };
    let PhysicalProcessingEvidence::LigationOutcomes(evidence) =
        &mut report.physical_processing_reports[0].evidence
    else {
        panic!("ligation outcomes")
    };
    assert_eq!(evidence.report.intended.sequence.as_str(), "AAAACCCC");
    assert!(!evidence.report.alternatives.is_empty());
    // An edited topology cannot be paired with the unchanged physical duplex.
    evidence.report.intended.topology = Topology::Circular;
    let payload = postcard::to_allocvec(&response).unwrap();
    let mut forged = b"CGR\0".to_vec();
    forged.extend(protocol::WORKER_SCHEMA_VERSION.to_le_bytes());
    forged.extend((payload.len() as u32).to_le_bytes());
    forged.extend(blake3::hash(&payload).as_bytes());
    forged.extend(payload);
    assert!(protocol::decode_response(&forged).is_err());
}
