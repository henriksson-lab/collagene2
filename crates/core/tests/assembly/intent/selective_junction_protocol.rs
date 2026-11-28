use collagene_core::{
    cloning::{
        assembly::{
            intent::{self as assembly_intent, Contribution, Selector},
            junction_edit as edit, selective_edit as selective,
        },
        pcr::primers::DesignOptions,
    },
    planning::plan::PlanInput,
    protocol::{
        self,
        request::Request,
        schema::{AnalysisJob, AnalysisResult},
    },
};
use collagene_model::{Molecule, Topology};
fn fixture() -> (Vec<PlanInput>, selective::Request) {
    let inputs = [("l", "ATGACCCCGGGGTTTCCAGA"), ("r", "ATCCGTTAGCGACCTTGGAA")]
        .into_iter()
        .map(|(id, dna)| {
            let molecule = Molecule::dna(id, dna, Topology::Linear).unwrap();
            PlanInput {
                id: id.into(),
                revision: molecule.revision(),
                molecule,
            }
        })
        .collect();
    let base = assembly_intent::Request {
        name: "declared edited assembly".into(),
        fragments: vec![
            Contribution {
                role: "left".into(),
                input_id: "l".into(),
                selector: Selector::ReferenceRange {
                    start: 0,
                    end: 20,
                    reverse: false,
                },
            },
            Contribution {
                role: "right".into(),
                input_id: "r".into(),
                selector: Selector::ReferenceRange {
                    start: 0,
                    end: 20,
                    reverse: false,
                },
            },
        ],
        circular: false,
        material_bindings: vec![],
    };
    let mut options = DesignOptions::default();
    options.min_length = 8;
    options.max_length = 8;
    options.max_tm_difference_c = 100.;
    options.min_gc_fraction = 0.;
    options.max_gc_fraction = 1.;
    options.max_homopolymer = 200;
    options.max_three_prime_complementarity = 200;
    options.max_pairs = 1;
    (
        inputs,
        selective::Request {
            base: base.clone(),
            edits: vec![edit::Edit {
                junction: 0,
                remove_left: 1,
                remove_right: 2,
                inserted: "GC".into(),
            }],
            rows: base
                .fragments
                .iter()
                .map(|f| selective::Row::DesignPcr {
                    role: f.role.clone(),
                    input_id: f.input_id.clone(),
                    reverse: false,
                })
                .collect(),
            overlaps: vec![15],
            primer_options: options,
            profile_options: Default::default(),
        },
    )
}

#[test]
fn selected_row_protocol_is_append_only_and_context_bound() {
    let (inputs, request) = fixture();
    let revision = selective::design_revision(&inputs, &request).unwrap();
    let job = Request {
        id: 451,
        revision: revision.clone(),
        job: AnalysisJob::DesignSelectiveAssemblyPrimers {
            inputs: inputs.clone(),
            request: request.clone(),
        },
    };
    let bytes = protocol::request::Request::encode(&job).unwrap();
    let ((_, _, tag), _) = postcard::take_from_bytes::<(u32, &str, u32)>(&bytes[42..]).unwrap();
    assert_eq!(tag, 109);
    let mut old = bytes.clone();
    old[4..6].copy_from_slice(&44u16.to_le_bytes());
    assert!(protocol::execution::process_request(&old).is_err());
    let response = protocol::execution::process_request(&bytes).unwrap();
    let ((_, _, tag), _) = postcard::take_from_bytes::<(u32, &str, u32)>(&response[42..]).unwrap();
    assert_eq!(tag, 101);
    let AnalysisResult::SelectiveAssemblyPrimers(report) =
        protocol::response::Response::decode(&response)
            .unwrap()
            .result
    else {
        panic!("wrong result")
    };
    report.validate_binding(&revision).unwrap();
    assert_eq!(report.designs.len(), 2);
    let mut old = response;
    old[4..6].copy_from_slice(&44u16.to_le_bytes());
    assert!(protocol::response::Response::decode(&old).is_err());
    let mut changed = request;
    changed.overlaps[0] = 16;
    let stale = protocol::request::Request::encode(&Request {
        id: 452,
        revision,
        job: AnalysisJob::DesignSelectiveAssemblyPrimers {
            inputs,
            request: changed,
        },
    })
    .unwrap();
    assert!(protocol::execution::process_request(&stale).is_err());
}
