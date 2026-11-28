use collagene_core::{
    assembly_intent::{self, Contribution, Selector},
    assembly_junction_edit as edit,
    primers::DesignOptions,
    protocol::{
        self,
        request::Request,
        schema::{AnalysisJob, AnalysisResult},
    },
    verification::PlanInput,
};
use collagene_model::{Molecule, Topology};
fn fixture() -> (Vec<PlanInput>, edit::Request) {
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
        edit::Request {
            base,
            edits: vec![edit::Edit {
                junction: 0,
                remove_left: 1,
                remove_right: 2,
                inserted: "GC".into(),
            }],
            overlap: 15,
            primer_options: options,
            profile_options: Default::default(),
        },
    )
}
#[test]
fn edited_intent_and_primer_proposal_have_stable_tags_and_reject_pre_44_headers() {
    let (inputs, request) = fixture();
    let revision = edit::intent_revision(&inputs, &request.base, &request.edits, None).unwrap();
    let job = Request {
        id: 440,
        revision: revision.clone(),
        job: AnalysisJob::BuildEditedAssemblyIntent {
            inputs: inputs.clone(),
            base: request.base.clone(),
            edits: request.edits.clone(),
            workspace: None,
        },
    };
    let bytes = protocol::request::Request::encode(&job).unwrap();
    let ((_, _, tag), _) = postcard::take_from_bytes::<(u32, &str, u32)>(&bytes[42..]).unwrap();
    assert_eq!(tag, 107);
    let mut old = bytes.clone();
    old[4..6].copy_from_slice(&43u16.to_le_bytes());
    assert!(protocol::execution::process_request(&old).is_err());
    let bytes = protocol::execution::process_request(&bytes).unwrap();
    let ((_, _, tag), _) = postcard::take_from_bytes::<(u32, &str, u32)>(&bytes[42..]).unwrap();
    assert_eq!(tag, 99);
    let AnalysisResult::EditedAssemblyIntent(intent) =
        protocol::response::Response::decode(&bytes).unwrap().result
    else {
        panic!("expected independent intended reference")
    };
    intent.validate_binding(&revision).unwrap();
    assert_eq!(
        intent.target.sequence.as_str(),
        "ATGACCCCGGGGTTTCCAGGCCCGTTAGCGACCTTGGAA"
    );
    let mut old = bytes;
    old[4..6].copy_from_slice(&43u16.to_le_bytes());
    assert!(protocol::response::Response::decode(&old).is_err());
    let revision = edit::design_revision(&inputs, &request).unwrap();
    let bytes = protocol::request::Request::encode(&Request {
        id: 441,
        revision: revision.clone(),
        job: AnalysisJob::DesignEditedAssemblyPrimers { inputs, request },
    })
    .unwrap();
    let ((_, _, tag), _) = postcard::take_from_bytes::<(u32, &str, u32)>(&bytes[42..]).unwrap();
    assert_eq!(tag, 108);
    let bytes = protocol::execution::process_request(&bytes).unwrap();
    let ((_, _, tag), _) = postcard::take_from_bytes::<(u32, &str, u32)>(&bytes[42..]).unwrap();
    assert_eq!(tag, 100);
    let AnalysisResult::EditedAssemblyPrimers(report) =
        protocol::response::Response::decode(&bytes).unwrap().result
    else {
        panic!("expected original-source primer proposal")
    };
    report.validate_binding(&revision).unwrap();
    assert_eq!(report.intent.target.revision(), intent.target.revision());
    assert_eq!(report.recipe.fragments.len(), 2);
}
#[test]
fn changing_a_declared_junction_cannot_reuse_the_old_intent_context() {
    let (inputs, mut request) = fixture();
    let revision = edit::intent_revision(&inputs, &request.base, &request.edits, None).unwrap();
    request.edits[0].inserted = "TT".into();
    let bytes = protocol::request::Request::encode(&Request {
        id: 442,
        revision,
        job: AnalysisJob::BuildEditedAssemblyIntent {
            inputs,
            base: request.base,
            edits: request.edits,
            workspace: None,
        },
    })
    .unwrap();
    assert!(protocol::execution::process_request(&bytes).is_err());
}
