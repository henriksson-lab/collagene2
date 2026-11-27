use collagene_core::{
    Molecule, Topology,
    golden_gate::{GoldenGateEnzyme, GoldenGateOptions},
    golden_gate_primer_design as design, ligation_recipe, prepared_golden_gate as pool,
    protocol::{self, AnalysisJob, AnalysisResult, Request},
    verification::{self, PlanInput},
};

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
    let bytes = protocol::encode_request(&Request {
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
    assert!(protocol::process_request(&old).is_err());
    let bytes = protocol::process_request(&bytes).unwrap();
    let ((_, _, tag), _) = postcard::take_from_bytes::<(u32, &str, u32)>(&bytes[42..]).unwrap();
    assert_eq!(tag, 103);
    let AnalysisResult::PreparedGoldenGatePool(report) =
        protocol::decode_response(&bytes).unwrap().result
    else {
        panic!("expected actual prepared pool")
    };
    report.validate_binding(&revision).unwrap();
    let mut old = bytes;
    old[4..6].copy_from_slice(&46u16.to_le_bytes());
    assert!(protocol::decode_response(&old).is_err());
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
    let bytes = protocol::encode_request(&Request {
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
    assert!(protocol::process_request(&old).is_err());
    let bytes = protocol::process_request(&bytes).unwrap();
    let ((_, _, tag), _) = postcard::take_from_bytes::<(u32, &str, u32)>(&bytes[42..]).unwrap();
    assert_eq!(tag, 104);
    let AnalysisResult::PreparedGoldenGateCandidate(selected) =
        protocol::decode_response(&bytes).unwrap().result
    else {
        panic!("expected unverified candidate")
    };
    selected.validate_binding(&selected_revision).unwrap();
    assert_eq!(
        verification::verify(&selected.plan, &target).status,
        verification::Status::ConsistentUnderAssumptions
    );
    let mut old = bytes;
    old[4..6].copy_from_slice(&46u16.to_le_bytes());
    assert!(protocol::decode_response(&old).is_err());
    let mut changed = request;
    changed.options.max_states = 1;
    let bytes = protocol::encode_request(&Request {
        id: 472,
        revision,
        job: AnalysisJob::InspectPreparedGoldenGatePool {
            inputs,
            target,
            request: changed,
        },
    })
    .unwrap();
    assert!(protocol::process_request(&bytes).is_err());
}
