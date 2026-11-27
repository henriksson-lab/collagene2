//! Literal contribution edits are independent of primer/product DNA across binary transport.
use collagene_core::{
    Molecule, Topology,
    assembly_edit_intent::{self, Edit, Request as IntentRequest},
    assembly_intent,
    assembly_intent::{Contribution, Selector},
    assembly_substitution_intent,
    protocol::{self, AnalysisJob, AnalysisResult, Request},
    verification::PlanInput,
};
fn wire(magic: &[u8; 4], version: u16, value: &impl serde::Serialize) -> Vec<u8> {
    let payload = postcard::to_allocvec(value).unwrap();
    let mut bytes = magic.to_vec();
    bytes.extend(version.to_le_bytes());
    bytes.extend((payload.len() as u32).to_le_bytes());
    bytes.extend(blake3::hash(&payload).as_bytes());
    bytes.extend(payload);
    bytes
}
fn declaration() -> (Vec<PlanInput>, IntentRequest) {
    let molecule = Molecule::dna("Unchanged original", "AACCGGTT", Topology::Linear).unwrap();
    let inputs = vec![PlanInput {
        id: "original".into(),
        revision: molecule.revision(),
        molecule,
    }];
    let request = IntentRequest {
        base: assembly_substitution_intent::Request {
            base: assembly_intent::Request {
                name: "Independent edited requirement".into(),
                fragments: vec![
                    Contribution {
                        role: "left".into(),
                        input_id: "original".into(),
                        selector: Selector::ReferenceRange {
                            start: 0,
                            end: 4,
                            reverse: false,
                        },
                    },
                    Contribution {
                        role: "right".into(),
                        input_id: "original".into(),
                        selector: Selector::ReferenceRange {
                            start: 0,
                            end: 4,
                            reverse: true,
                        },
                    },
                ],
                circular: false,
                material_bindings: vec![],
            },
            junction_edits: vec![],
            substitutions: vec![],
        },
        edits: vec![
            Edit {
                role: "left".into(),
                start: 2,
                end: 2,
                expected: "".into(),
                replacement: "TT".into(),
            },
            Edit {
                role: "right".into(),
                start: 1,
                end: 3,
                expected: "GT".into(),
                replacement: "A".into(),
            },
        ],
    };
    (inputs, request)
}
#[test]
fn contribution_edits_append_163_155_and_reject_old_or_stale_bindings() {
    assert!(protocol::WORKER_SCHEMA_VERSION >= 79);
    let (inputs, declaration) = declaration();
    let originals = inputs.clone();
    let revision = assembly_edit_intent::intent_revision(&inputs, &declaration, None).unwrap();
    let request = Request {
        id: 163,
        revision: revision.clone(),
        job: AnalysisJob::BuildContributionEditIntent {
            inputs,
            request: declaration,
            workspace: None,
        },
    };
    let encoded = protocol::encode_request(&request).unwrap();
    let ((_, _, tag), _) = postcard::take_from_bytes::<(u32, &str, u32)>(&encoded[42..]).unwrap();
    assert_eq!(tag, 163);
    for version in [1, 76, 77, 78] {
        assert!(protocol::process_request(&wire(b"CGQ\0", version, &request)).is_err());
    }
    let encoded = protocol::process_request(&encoded).unwrap();
    let ((_, _, tag), _) = postcard::take_from_bytes::<(u32, &str, u32)>(&encoded[42..]).unwrap();
    assert_eq!(tag, 155);
    let response = protocol::decode_response(&encoded).unwrap();
    let AnalysisResult::ContributionEditIntent(report) = &response.result else {
        panic!("literal intent")
    };
    report.validate_binding(&revision).unwrap();
    report
        .validate_sources(&originals, &revision, None)
        .unwrap();
    assert_eq!(report.target.sequence.as_str(), "AATTCCGAT");
    assert_eq!(report.applied.len(), 2);
    assert_eq!(report.fragments.len(), 2);
    assert!(collagene_core::chemistry::has_unmapped_physical_evidence(
        &report.target
    ));
    for version in [1, 76, 77, 78] {
        assert!(protocol::decode_response(&wire(b"CGR\0", version, &response)).is_err());
    }
    let mut stale = request;
    stale.revision = "0".repeat(64);
    assert!(protocol::process_request(&protocol::encode_request(&stale).unwrap()).is_err());
    let mut stale = response.clone();
    stale.revision = "0".repeat(64);
    assert!(
        protocol::decode_response(&wire(b"CGR\0", protocol::WORKER_SCHEMA_VERSION, &stale))
            .is_err()
    );
    let mut forged = response;
    let AnalysisResult::ContributionEditIntent(report) = &mut forged.result else {
        unreachable!()
    };
    report.target.sequence = collagene_core::Sequence::dna("AAAAAAAAA").unwrap();
    assert!(
        protocol::decode_response(&wire(b"CGR\0", protocol::WORKER_SCHEMA_VERSION, &forged))
            .is_err()
    );
}
#[test]
fn wrong_original_literal_and_revision_are_rejected_before_worker_execution() {
    let (mut inputs, mut request) = declaration();
    request.edits[1].expected = "AC".into();
    assert!(assembly_edit_intent::intent_revision(&inputs, &request, None).is_err());
    request.edits[1].expected = "GT".into();
    let revision = assembly_edit_intent::intent_revision(&inputs, &request, None).unwrap();
    inputs[0].revision = "0".repeat(64);
    let request = Request {
        id: 7,
        revision,
        job: AnalysisJob::BuildContributionEditIntent {
            inputs,
            request,
            workspace: None,
        },
    };
    assert!(protocol::process_request(&protocol::encode_request(&request).unwrap()).is_err());
}
