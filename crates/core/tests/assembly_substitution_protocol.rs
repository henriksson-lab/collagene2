//! Independently declared base changes preserve original source bindings across the worker.
use collagene_core::{
    Molecule, Topology, assembly_intent,
    assembly_intent::{Contribution, Selector},
    assembly_junction_edit::Edit,
    assembly_substitution_intent::{self, Request as IntentRequest, Substitution},
    protocol::{self, AnalysisJob, AnalysisResult, Request},
    verification::PlanInput,
};
fn wire(magic: &[u8; 4], v: u16, value: &impl serde::Serialize) -> Vec<u8> {
    let p = postcard::to_allocvec(value).unwrap();
    let mut b = magic.to_vec();
    b.extend(v.to_le_bytes());
    b.extend((p.len() as u32).to_le_bytes());
    b.extend(blake3::hash(&p).as_bytes());
    b.extend(p);
    b
}
fn declaration(
    reverse: bool,
    cross_origin: bool,
) -> (Vec<PlanInput>, IntentRequest, &'static str, usize) {
    let molecule = Molecule::dna(
        "Original",
        if cross_origin { "AAGTCC" } else { "AACCGGTT" },
        if cross_origin {
            Topology::Circular
        } else {
            Topology::Linear
        },
    )
    .unwrap();
    let input = PlanInput {
        id: "source".into(),
        revision: molecule.revision(),
        molecule,
    };
    let r = IntentRequest {
        base: assembly_intent::Request {
            name: "Independent substitution intent".into(),
            fragments: vec![Contribution {
                role: "one".into(),
                input_id: "source".into(),
                selector: Selector::ReferenceRange {
                    start: if cross_origin { 4 } else { 0 },
                    end: 8,
                    reverse,
                },
            }],
            circular: cross_origin,
            material_bindings: vec![],
        },
        junction_edits: vec![],
        substitutions: vec![Substitution {
            role: "one".into(),
            offset: 2,
            expected: if cross_origin { "A" } else { "C" }.into(),
            replacement: if cross_origin { "G" } else { "A" }.into(),
        }],
    };
    (
        vec![input],
        r,
        if cross_origin { "CCGA" } else { "AAACGGTT" },
        if cross_origin {
            0
        } else if reverse {
            5
        } else {
            2
        },
    )
}
#[test]
fn substituted_intent_161_153_roundtrips_forward_reverse_and_origin_with_literal_source_positions()
{
    assert!(protocol::WORKER_SCHEMA_VERSION >= 77);
    for (reverse, cross) in [(false, false), (true, false), (false, true)] {
        let (inputs, request, expected, source_position) = declaration(reverse, cross);
        let before = inputs[0].molecule.clone();
        let revision =
            assembly_substitution_intent::intent_revision(&inputs, &request, None).unwrap();
        let r = Request {
            id: 161,
            revision: revision.clone(),
            job: AnalysisJob::BuildSubstitutedAssemblyIntent {
                inputs,
                request,
                workspace: None,
            },
        };
        let bytes = protocol::encode_request(&r).unwrap();
        let ((_, _, tag), _) = postcard::take_from_bytes::<(u32, &str, u32)>(&bytes[42..]).unwrap();
        assert_eq!(tag, 159);
        assert!(protocol::process_request(&wire(b"CGQ\0", 76, &r)).is_err());
        let bytes = protocol::process_request(&bytes).unwrap();
        let ((_, _, tag), _) = postcard::take_from_bytes::<(u32, &str, u32)>(&bytes[42..]).unwrap();
        assert_eq!(tag, 151);
        let response = protocol::decode_response(&bytes).unwrap();
        let AnalysisResult::SubstitutedAssemblyIntent(report) = &response.result else {
            panic!("intent")
        };
        report.validate_binding(&revision).unwrap();
        assert_eq!(report.target.sequence.as_str(), expected);
        assert_eq!(report.applied.len(), 1);
        assert_eq!(
            report.applied[0].source_position,
            collagene_core::assembly_intent::ContributionPosition::Reference {
                position: source_position,
                reverse
            }
        );
        assert_eq!(report.applied[0].input_revision, before.revision());
        assert!(protocol::decode_response(&wire(b"CGR\0", 76, &response)).is_err());
        let mut stale = r;
        stale.revision = "0".repeat(64);
        assert!(protocol::process_request(&protocol::encode_request(&stale).unwrap()).is_err());
        let mut stale = response;
        stale.revision = "0".repeat(64);
        assert!(
            protocol::decode_response(&wire(b"CGR\0", protocol::WORKER_SCHEMA_VERSION, &stale))
                .is_err()
        );
    }
}
#[test]
fn trims_and_literal_inserts_precede_retained_core_substitution_and_expected_base_is_checked() {
    let (inputs, mut request, _, _) = declaration(false, false);
    request.base.fragments = vec![
        Contribution {
            role: "left".into(),
            input_id: "source".into(),
            selector: Selector::ReferenceRange {
                start: 0,
                end: 4,
                reverse: false,
            },
        },
        Contribution {
            role: "right".into(),
            input_id: "source".into(),
            selector: Selector::ReferenceRange {
                start: 4,
                end: 8,
                reverse: false,
            },
        },
    ];
    request.junction_edits = vec![Edit {
        junction: 0,
        remove_left: 1,
        remove_right: 1,
        inserted: "AAA".into(),
    }];
    request.substitutions = vec![Substitution {
        role: "left".into(),
        offset: 1,
        expected: "A".into(),
        replacement: "T".into(),
    }];
    let revision = assembly_substitution_intent::intent_revision(&inputs, &request, None).unwrap();
    let r = Request {
        id: 2,
        revision,
        job: AnalysisJob::BuildSubstitutedAssemblyIntent {
            inputs: inputs.clone(),
            request: request.clone(),
            workspace: None,
        },
    };
    let response = protocol::decode_response(
        &protocol::process_request(&protocol::encode_request(&r).unwrap()).unwrap(),
    )
    .unwrap();
    let AnalysisResult::SubstitutedAssemblyIntent(report) = response.result else {
        panic!("intent")
    };
    assert_eq!(report.target.sequence.as_str(), "ATCAAAGTT");
    assert_eq!(
        report.applied[0].source_position,
        collagene_core::assembly_intent::ContributionPosition::Reference {
            position: 1,
            reverse: false
        }
    );
    assert_eq!(report.applied[0].output_position, 1);
    request.substitutions[0].expected = "C".into();
    assert!(assembly_substitution_intent::intent_revision(&inputs, &request, None).is_err());
}
