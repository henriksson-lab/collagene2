use collagene_core::{
    cloning::{
        golden_gate::{GoldenGateEnzyme, GoldenGateOptions, prepared as pool, primers as design},
        ligation::recipe as ligation_recipe,
    },
    planning::plan::PlanInput,
    protocol::{
        self,
        request::Request,
        schema::{AnalysisJob, AnalysisResult},
    },
    verification,
};
use collagene_model::{Molecule, Topology};

fn fixture() -> (Vec<PlanInput>, Molecule, pool::Request) {
    // Literal target is fixed before the independent preparation designer runs.
    let target = Molecule::dna(
        "desired",
        "ATGACCCCGGGGTTTCCAGAATGAGTTAGCGACCTTGGAA",
        Topology::Circular,
    )
    .unwrap();
    let inputs = [("a", "ATGACCCCGGGGTTTCCAGA"), ("b", "ATGAGTTAGCGACCTTGGAA")]
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
    let proposal = design::design(
        &inputs,
        &target,
        &design::Request {
            name: "original PCR and BsaI preparations".into(),
            fragments: inputs
                .iter()
                .map(|i| design::Fragment {
                    role: i.id.clone(),
                    input_id: i.id.clone(),
                    start: 0,
                    end: 20,
                    reverse: false,
                    forward_length: 8,
                    reverse_length: 8,
                    right_overhang: "ATGA".into(),
                })
                .collect(),
            enzyme: GoldenGateEnzyme::BsaI,
            clamp: "ATGC".into(),
            spacer: "A".into(),
        },
    )
    .unwrap();
    let request = pool::Request {
        recipe: ligation_recipe::Recipe {
            standard: None,
            name: "two-role pool".into(),
            fragments: proposal
                .fragments
                .into_iter()
                .map(|source| ligation_recipe::Fragment {
                    source,
                    treatment: None,
                })
                .collect(),
            circular: true,
            allow_nicked: false,
        },
        options: GoldenGateOptions {
            max_fragments: 2,
            ..Default::default()
        },
    };
    (inputs, target, request)
}

#[test]
fn prepared_pool_and_unverified_selection_have_append_only_bound_47_messages() {
    let (inputs, target, request) = fixture();
    let revision = pool::inspect_revision(&inputs, &target, &request).unwrap();
    let bytes = protocol::request::Request::encode(&Request {
        id: 470,
        revision: revision.clone(),
        job: AnalysisJob::InspectPreparedGoldenGatePool {
            inputs: inputs.clone(),
            target: target.clone(),
            request: request.clone(),
        },
    })
    .unwrap();
    let ((_, _, tag), _) = postcard::take_from_bytes::<(u32, &str, u32)>(&bytes[42..]).unwrap();
    assert_eq!(tag, 111);
    let mut old = bytes.clone();
    old[4..6].copy_from_slice(&46u16.to_le_bytes());
    assert!(protocol::execution::process_request(&old).is_err());
    let bytes = protocol::execution::process_request(&bytes).unwrap();
    let ((_, _, tag), _) = postcard::take_from_bytes::<(u32, &str, u32)>(&bytes[42..]).unwrap();
    assert_eq!(tag, 103);
    let AnalysisResult::PreparedGoldenGatePool(report) =
        protocol::response::Response::decode(&bytes).unwrap().result
    else {
        panic!("expected actual prepared pool")
    };
    report.validate_binding(&revision).unwrap();
    let mut old = bytes;
    old[4..6].copy_from_slice(&46u16.to_le_bytes());
    assert!(protocol::response::Response::decode(&old).is_err());
    let a = report
        .assessments
        .iter()
        .find(|a| a.target_matches)
        .unwrap();
    let selection = pool::Selection {
        context_revision: report.context_revision.clone(),
        candidate_index: a.candidate_index,
        directed_sequence_key: a.directed_sequence_key.clone(),
        physical_sequence_key: a.physical_sequence_key.clone(),
    };
    let selected_revision = pool::selection_revision(&report, &selection).unwrap();
    let bytes = protocol::request::Request::encode(&Request {
        id: 471,
        revision: selected_revision.clone(),
        job: AnalysisJob::SelectPreparedGoldenGateCandidate {
            report: Box::new(report),
            selection,
        },
    })
    .unwrap();
    let ((_, _, tag), _) = postcard::take_from_bytes::<(u32, &str, u32)>(&bytes[42..]).unwrap();
    assert_eq!(tag, 112);
    let mut old = bytes.clone();
    old[4..6].copy_from_slice(&46u16.to_le_bytes());
    assert!(protocol::execution::process_request(&old).is_err());
    let bytes = protocol::execution::process_request(&bytes).unwrap();
    let ((_, _, tag), _) = postcard::take_from_bytes::<(u32, &str, u32)>(&bytes[42..]).unwrap();
    assert_eq!(tag, 104);
    let AnalysisResult::PreparedGoldenGateCandidate(selected) =
        protocol::response::Response::decode(&bytes).unwrap().result
    else {
        panic!("expected unverified candidate")
    };
    selected.validate_binding(&selected_revision).unwrap();
    assert_eq!(
        verification::engine::verify(&selected.plan, &target).status,
        verification::evidence::Status::ConsistentUnderAssumptions
    );
    let mut old = bytes;
    old[4..6].copy_from_slice(&46u16.to_le_bytes());
    assert!(protocol::response::Response::decode(&old).is_err());
    let mut changed = request;
    changed.options.max_states = 1;
    let bytes = protocol::request::Request::encode(&Request {
        id: 472,
        revision,
        job: AnalysisJob::InspectPreparedGoldenGatePool {
            inputs,
            target,
            request: changed,
        },
    })
    .unwrap();
    assert!(protocol::execution::process_request(&bytes).is_err());
}
