use collagene_core::{
    cloning::{
        assembly::{
            intent::{self as assembly_intent, Contribution, Selector},
            junction_edit::Edit,
        },
        golden_gate::{GoldenGateEnzyme, junction_edit as edit},
    },
    planning::plan::PlanInput,
    protocol::{
        self,
        request::Request,
        schema::{AnalysisJob, AnalysisResult},
    },
};
use collagene_model::{Molecule, Topology};

fn fixture() -> (Vec<PlanInput>, edit::Request) {
    let inputs = [("a", "ATGACCCCGGGGTTTCCAGA"), ("b", "ATCCGTTAGCGACCTTGGAA")]
        .into_iter()
        .map(|(id, dna)| {
            let molecule = Molecule::dna(id, dna, Topology::Linear).unwrap();
            PlanInput {
                id: id.into(),
                revision: molecule.revision(),
                molecule,
            }
        })
        .collect::<Vec<_>>();
    let base = assembly_intent::Request {
        name: "source-declared edited circle".into(),
        fragments: inputs
            .iter()
            .map(|input| Contribution {
                role: input.id.clone(),
                input_id: input.id.clone(),
                selector: Selector::ReferenceRange {
                    start: 0,
                    end: 20,
                    reverse: false,
                },
            })
            .collect(),
        circular: true,
        material_bindings: vec![],
    };
    let rows = base
        .fragments
        .iter()
        .map(|f| edit::Row::DesignPcr {
            role: f.role.clone(),
            input_id: f.input_id.clone(),
            reverse: false,
            forward_length: 8,
            reverse_length: 8,
        })
        .collect();
    (
        inputs,
        edit::Request {
            base,
            rows,
            edits: vec![Edit {
                junction: 0,
                remove_left: 0,
                remove_right: 0,
                inserted: "GC".into(),
            }],
            enzyme: GoldenGateEnzyme::BsaI,
            clamp: "ATGC".into(),
            spacer: "A".into(),
            terminal_right_overhang: None,
            allow_nicked: false,
        },
    )
}

#[test]
fn golden_gate_junction_messages_bind_edits_and_reject_pre_46_selectors() {
    let (inputs, request) = fixture();
    let revision = edit::design_revision(&inputs, &request).unwrap();
    let bytes = protocol::request::Request::encode(&Request {
        id: 460,
        revision: revision.clone(),
        job: AnalysisJob::DesignGoldenGateJunctionPrimers {
            inputs: inputs.clone(),
            request: request.clone(),
        },
    })
    .unwrap();
    let ((_, _, tag), _) = postcard::take_from_bytes::<(u32, &str, u32)>(&bytes[42..]).unwrap();
    assert_eq!(tag, 110);
    let mut old = bytes.clone();
    old[4..6].copy_from_slice(&45u16.to_le_bytes());
    assert!(protocol::execution::process_request(&old).is_err());
    let response = protocol::execution::process_request(&bytes).unwrap();
    let ((_, _, tag), _) = postcard::take_from_bytes::<(u32, &str, u32)>(&response[42..]).unwrap();
    assert_eq!(tag, 102);
    let AnalysisResult::GoldenGateJunctionPrimers(report) =
        protocol::response::Response::decode(&response)
            .unwrap()
            .result
    else {
        panic!("expected Golden Gate junction proposal")
    };
    report.validate_binding(&revision).unwrap();
    assert_eq!(
        report.intent.target.sequence.as_str(),
        "ATGACCCCGGGGTTTCCAGAGCATCCGTTAGCGACCTTGGAA"
    );
    assert_eq!(report.intent.target.topology, Topology::Circular);
    let mut old = response;
    old[4..6].copy_from_slice(&45u16.to_le_bytes());
    assert!(protocol::response::Response::decode(&old).is_err());
    let mut changed = request;
    changed.edits[0].inserted = "TT".into();
    let stale = protocol::request::Request::encode(&Request {
        id: 461,
        revision,
        job: AnalysisJob::DesignGoldenGateJunctionPrimers {
            inputs,
            request: changed,
        },
    })
    .unwrap();
    assert!(protocol::execution::process_request(&stale).is_err());
}
