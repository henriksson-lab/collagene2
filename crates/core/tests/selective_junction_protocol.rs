use collagene_core::{
    Molecule, Topology, assembly_intent, assembly_junction_edit as edit,
    overlap_primer_design::Fragment,
    primers::DesignOptions,
    protocol::{self, AnalysisJob, AnalysisResult, Request},
    selective_junction_edit as selective,
    verification::PlanInput,
};
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
            Fragment {
                role: "left".into(),
                input_id: "l".into(),
                start: 0,
                end: 20,
                reverse: false,
            },
            Fragment {
                role: "right".into(),
                input_id: "r".into(),
                start: 0,
                end: 20,
                reverse: false,
            },
        ],
        circular: false,
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
                    reverse: f.reverse,
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
    let bytes = protocol::encode_request(&job).unwrap();
    let ((_, _, tag), _) = postcard::take_from_bytes::<(u32, &str, u32)>(&bytes[42..]).unwrap();
    assert_eq!(tag, 111);
    let mut old = bytes.clone();
    old[4..6].copy_from_slice(&44u16.to_le_bytes());
    assert!(protocol::process_request(&old).is_err());
    let response = protocol::process_request(&bytes).unwrap();
    let ((_, _, tag), _) = postcard::take_from_bytes::<(u32, &str, u32)>(&response[42..]).unwrap();
    assert_eq!(tag, 103);
    let AnalysisResult::SelectiveAssemblyPrimers(report) =
        protocol::decode_response(&response).unwrap().result
    else {
        panic!("wrong result")
    };
    report.validate_binding(&revision).unwrap();
    assert_eq!(report.designs.len(), 2);
    let mut old = response;
    old[4..6].copy_from_slice(&44u16.to_le_bytes());
    assert!(protocol::decode_response(&old).is_err());
    let mut changed = request;
    changed.overlaps[0] = 16;
    let stale = protocol::encode_request(&Request {
        id: 452,
        revision,
        job: AnalysisJob::DesignSelectiveAssemblyPrimers {
            inputs,
            request: changed,
        },
    })
    .unwrap();
    assert!(protocol::process_request(&stale).is_err());
}
