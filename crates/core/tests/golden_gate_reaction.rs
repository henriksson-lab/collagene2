use collagene_core::{
    Molecule, Topology,
    analysis::Primer,
    assembly_recipe::Preparation,
    golden_gate::GoldenGateOptions,
    golden_gate_reaction as reaction,
    prepared_golden_gate::Selection,
    protocol::{self, AnalysisJob, AnalysisResult, Request},
    verification::{self, Operation, PlanInput, Status},
};
const A: &str = "ATGACCCCGGGGTTTCCAGA";
const B: &str = "ATGAGTTAGCGACCTTGGAA";
const AB: &str = "ATGACCCCGGGGTTTCCAGAATGAGTTAGCGACCTTGGAA";
fn fixture() -> (Vec<PlanInput>, Molecule, reaction::Request) {
    let target = Molecule::dna("independent insert construct", AB, Topology::Circular).unwrap();
    let sources = [
        (
            "vector",
            format!("ATGCGGTCTCA{A}ATGATGAGACCGCAT"),
            Topology::Circular,
        ),
        ("insert", B.into(), Topology::Linear),
        ("carry", "ACGTACGTACGT".into(), Topology::Circular),
    ];
    let inputs = sources
        .into_iter()
        .map(|(id, dna, topology)| {
            let molecule = Molecule::dna(id, &dna, topology).unwrap();
            PlanInput {
                id: id.into(),
                revision: molecule.revision(),
                molecule,
            }
        })
        .collect();
    let request = reaction::Request {
        name: "all material entering the reaction".into(),
        rows: vec![
            reaction::Row {
                role: "whole vector".into(),
                input_id: "vector".into(),
                substrate: reaction::Substrate::WholeSource,
            },
            reaction::Row {
                role: "complete amplicon".into(),
                input_id: "insert".into(),
                substrate: reaction::Substrate::Prepared {
                    preparation: Preparation::ExactPcr {
                        forward: Primer::new("manual forward", "ATGAGTTA", "ATGCGGTCTCA").unwrap(),
                        reverse: Primer::new("manual reverse", "TTCCAAGG", "ATGCGGTCTCATCAT")
                            .unwrap(),
                    },
                    treatment: None,
                },
            },
            reaction::Row {
                role: "uncut contaminant".into(),
                input_id: "carry".into(),
                substrate: reaction::Substrate::WholeSource,
            },
        ],
        options: GoldenGateOptions {
            max_fragments: 3,
            ..Default::default()
        },
    };
    (inputs, target, request)
}
#[test]
fn whole_circle_and_whole_amplicon_include_every_digest_piece_and_uncut_carryover() {
    let (inputs, target, request) = fixture();
    let before = serde_json::to_vec(&inputs).unwrap();
    let report = reaction::inspect(&inputs, &target, &request).unwrap();
    let pool = report.pool.as_ref().unwrap();
    // Two circular-vector fragments, three linear-amplicon fragments, one uncut circle.
    assert_eq!(pool.fragments.len(), 6);
    assert_eq!(pool.carryover_sources, vec![2]);
    assert_eq!(pool.sources.len(), 3);
    assert_eq!(report.probe_plan.steps.len(), 2);
    assert!(matches!(report.probe_plan.steps[0], Operation::Pcr { .. }));
    assert!(
        pool.candidates.iter().any(|c| !c.restriction_site_free),
        "discarded vector material must remain available for recutting backgrounds"
    );
    let matching = report
        .assessments
        .iter()
        .find(|a| a.target_matches && pool.candidates[a.candidate_index].restriction_site_free)
        .unwrap();
    assert_eq!(matching.source_roles, vec![0, 1]);
    let s = Selection {
        context_revision: report.context_revision.clone(),
        candidate_index: matching.candidate_index,
        directed_sequence_key: matching.directed_sequence_key.clone(),
        physical_sequence_key: matching.physical_sequence_key.clone(),
    };
    let selected = reaction::select(&report, &s).unwrap();
    assert_eq!(
        serde_json::to_value(&selected.plan.steps[..1]).unwrap(),
        serde_json::to_value(&report.probe_plan.steps[..1]).unwrap()
    );
    assert_eq!(
        verification::verify(&selected.plan, &target).status,
        Status::ConsistentUnderAssumptions
    );
    let self_a = report
        .assessments
        .iter()
        .find(|a| pool.candidates[a.candidate_index].duplex.top.as_str() == A)
        .unwrap();
    let self_plan = reaction::select(
        &report,
        &Selection {
            context_revision: report.context_revision.clone(),
            candidate_index: self_a.candidate_index,
            directed_sequence_key: self_a.directed_sequence_key.clone(),
            physical_sequence_key: self_a.physical_sequence_key.clone(),
        },
    )
    .unwrap();
    assert_eq!(
        verification::verify(&self_plan.plan, &target).status,
        Status::Invalid
    );
    assert_eq!(serde_json::to_vec(&inputs).unwrap(), before);
}
#[test]
fn empty_capped_or_failed_probe_keeps_its_real_enumeration_without_selecting_a_winner() {
    let (inputs, target, mut request) = fixture();
    request.options.max_states = 1;
    request.options.min_fragments = 2;
    let capped = reaction::inspect(&inputs, &target, &request).unwrap();
    let pool = capped
        .pool
        .as_ref()
        .expect("enumeration must survive an unavailable index0 probe");
    assert!(!pool.enumeration_complete);
    assert!(pool.candidates.is_empty());
    assert!(!pool.truncation_reasons.is_empty());
    assert_eq!(capped.verification.products.len(), 1);
    let mut missing = request.clone();
    missing.options.max_states = 50_000;
    if let reaction::Substrate::Prepared {
        preparation: Preparation::ExactPcr { forward, .. },
        ..
    } = &mut missing.rows[1].substrate
    {
        *forward = Primer::new("bad", "AAAAAAAA", "").unwrap();
    }
    assert!(
        reaction::inspect(&inputs, &target, &missing)
            .unwrap()
            .pool
            .is_none()
    );
    let mut poisoned = inputs.clone();
    poisoned[0]
        .molecule
        .metadata
        .insert("duplex_model".into(), "unverified-reference-only".into());
    poisoned[0].revision = poisoned[0].molecule.revision();
    let mut valid = request.clone();
    valid.options.max_states = 50_000;
    assert!(
        reaction::inspect(&poisoned, &target, &valid)
            .unwrap()
            .pool
            .is_none()
    );
}
#[test]
fn reaction_protocol_is_append_only_and_changed_keys_or_substrates_cannot_reuse_context() {
    let (inputs, target, request) = fixture();
    let revision = reaction::inspect_revision(&inputs, &target, &request).unwrap();
    let bytes = protocol::encode_request(&Request {
        id: 490,
        revision: revision.clone(),
        job: AnalysisJob::InspectGoldenGateReaction {
            inputs: inputs.clone(),
            target: target.clone(),
            request: request.clone(),
        },
    })
    .unwrap();
    let ((_, _, tag), _) = postcard::take_from_bytes::<(u32, &str, u32)>(&bytes[42..]).unwrap();
    assert_eq!(tag, 115);
    let mut old = bytes.clone();
    old[4..6].copy_from_slice(&48u16.to_le_bytes());
    assert!(protocol::process_request(&old).is_err());
    let bytes = protocol::process_request(&bytes).unwrap();
    let ((_, _, tag), _) = postcard::take_from_bytes::<(u32, &str, u32)>(&bytes[42..]).unwrap();
    assert_eq!(tag, 107);
    let AnalysisResult::GoldenGateReaction(report) =
        protocol::decode_response(&bytes).unwrap().result
    else {
        panic!("reaction result required")
    };
    let mut old = bytes;
    old[4..6].copy_from_slice(&48u16.to_le_bytes());
    assert!(protocol::decode_response(&old).is_err());
    let a = report
        .assessments
        .iter()
        .find(|a| a.target_matches)
        .unwrap();
    let mut s = Selection {
        context_revision: revision.clone(),
        candidate_index: a.candidate_index,
        directed_sequence_key: a.directed_sequence_key.clone(),
        physical_sequence_key: a.physical_sequence_key.clone(),
    };
    let expected = reaction::selection_revision(&report, &s).unwrap();
    let encoded = protocol::encode_request(&Request {
        id: 491,
        revision: expected.clone(),
        job: AnalysisJob::SelectGoldenGateReaction {
            report: Box::new(report.clone()),
            selection: s.clone(),
        },
    })
    .unwrap();
    let ((_, _, tag), _) = postcard::take_from_bytes::<(u32, &str, u32)>(&encoded[42..]).unwrap();
    assert_eq!(tag, 116);
    let response = protocol::process_request(&encoded).unwrap();
    let ((_, _, tag), _) = postcard::take_from_bytes::<(u32, &str, u32)>(&response[42..]).unwrap();
    assert_eq!(tag, 108);
    let AnalysisResult::GoldenGateReactionCandidate(selected) =
        protocol::decode_response(&response).unwrap().result
    else {
        panic!("candidate required")
    };
    selected.validate_binding(&expected).unwrap();
    s.directed_sequence_key = "0".repeat(64);
    assert!(reaction::select(&report, &s).is_err());
    let mut changed = report.clone();
    changed.request.rows[0].substrate = reaction::Substrate::Prepared {
        preparation: Preparation::WholeLinear,
        treatment: None,
    };
    assert!(changed.validate_binding(&revision).is_err());
    let mut stale = inputs;
    stale[0].revision = "0".repeat(64);
    assert!(reaction::inspect(&stale, &target, &request).is_err());
}
