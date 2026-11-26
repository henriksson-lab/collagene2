//! Append-only outcome evidence must not acquire meaning under older envelopes.
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
    steps.push(Operation::LigateWithOutcomes {
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
fn operation29_and_evidence8_keep_intended_and_extractable_alternatives() {
    let (plan, target) = fixture();
    assert_eq!(tag(&plan.steps[2]), 29);
    let report = verification::verify(&plan, &target);
    assert_eq!(
        report.status,
        Status::ConsistentUnderAssumptions,
        "{:?}",
        report.findings
    );
    assert_eq!(report.products[2].sequence.as_str(), "AAAACCCC");
    let evidence = &report.physical_processing_reports[0].evidence;
    assert_eq!(tag(evidence), 8);
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
fn new_semantics_cannot_hide_in_old_worker_or_candidate_envelopes() {
    let (plan, target) = fixture();
    let request = protocol::Request {
        id: 29,
        revision: target.revision(),
        job: AnalysisJob::Verify {
            plan: plan.clone(),
            target: target.clone(),
        },
    };
    let bytes = protocol::encode_request(&request).unwrap();
    let mut old = bytes.clone();
    old[4..6].copy_from_slice(&91u16.to_le_bytes());
    assert_eq!(&old[6..], &bytes[6..]);
    assert!(protocol::process_request(&old).is_err());
    let response_bytes = protocol::process_request(&bytes).unwrap();
    let mut response = protocol::decode_response(&response_bytes).unwrap();
    // Removing only the engine-version fence must still leave the new evidence
    // capability fence intact.
    let AnalysisResult::Verification(report) = &mut response.result else {
        panic!("verification")
    };
    report.engine_version = report.engine_version.replace("/science-26", "/science-25");
    let payload = postcard::to_allocvec(&response).unwrap();
    let mut bytes = b"CGR\0".to_vec();
    bytes.extend(protocol::WORKER_SCHEMA_VERSION.to_le_bytes());
    bytes.extend((payload.len() as u32).to_le_bytes());
    bytes.extend(blake3::hash(&payload).as_bytes());
    bytes.extend(payload);
    protocol::decode_response(&bytes).unwrap();
    let mut old = bytes.clone();
    old[4..6].copy_from_slice(&91u16.to_le_bytes());
    assert_eq!(&old[6..], &bytes[6..]);
    assert!(protocol::decode_response(&old).is_err());
    let bytes = transfer::save_candidate(&transfer::Candidate::new(plan)).unwrap();
    assert_eq!(&bytes[4..6], &23u16.to_le_bytes());
    let mut old = bytes.clone();
    old[4..6].copy_from_slice(&22u16.to_le_bytes());
    assert_eq!(&old[6..], &bytes[6..]);
    assert!(transfer::load_candidate(&old).is_err());
    assert_eq!(workspace::PLAN_SCHEMA_VERSION, 21);
}
