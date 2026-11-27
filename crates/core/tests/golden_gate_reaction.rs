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
                generation: Some(collagene_core::modified_duplex::GenerationAssumptions {
                    dam: collagene_core::methylation::MethylationState::Unmodified,
                    dcm: collagene_core::methylation::MethylationState::Unmodified,
                    cpg: collagene_core::methylation::MethylationState::Unmodified,
                }),
                role: "whole vector".into(),
                input_id: "vector".into(),
                substrate: reaction::Substrate::WholeSource { origin: None },
            },
            reaction::Row {
                generation: Some(collagene_core::modified_duplex::GenerationAssumptions {
                    dam: collagene_core::methylation::MethylationState::Unmodified,
                    dcm: collagene_core::methylation::MethylationState::Unmodified,
                    cpg: collagene_core::methylation::MethylationState::Unmodified,
                }),
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
                generation: Some(collagene_core::modified_duplex::GenerationAssumptions {
                    dam: collagene_core::methylation::MethylationState::Unmodified,
                    dcm: collagene_core::methylation::MethylationState::Unmodified,
                    cpg: collagene_core::methylation::MethylationState::Unmodified,
                }),
                role: "uncut contaminant".into(),
                input_id: "carry".into(),
                substrate: reaction::Substrate::WholeSource { origin: None },
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
    assert_eq!(report.probe_plan.steps.len(), 5);
    assert!(matches!(
        report.probe_plan.steps[0],
        Operation::DeclareSourceModifications { .. }
    ));
    assert!(matches!(report.probe_plan.steps[1], Operation::Pcr { .. }));
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
        serde_json::to_value(&selected.plan.steps[..4]).unwrap(),
        serde_json::to_value(&report.probe_plan.steps[..4]).unwrap()
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
    assert_eq!(
        capped.verification.products.len(),
        4,
        "three declared-generation products and PCR remain; no capped candidate was fabricated"
    );
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

#[test]
fn key_bound_native_operation_reconstructs_the_selected_path_without_target_selection() {
    let (inputs, target, request) = fixture();
    let discovery = reaction::inspect(&inputs, &target, &request).unwrap();
    let pool = discovery.pool.as_ref().unwrap();
    let index = pool
        .candidates
        .iter()
        .position(|candidate| {
            candidate.duplex.top.as_str() == AB && candidate.restriction_site_free
        })
        .unwrap();
    let candidate = &pool.candidates[index];
    let selection = Selection {
        context_revision: discovery.context_revision.clone(),
        candidate_index: index,
        directed_sequence_key: candidate.directed_sequence_key.clone(),
        physical_sequence_key: candidate.physical_sequence_key.clone(),
    };
    let plan = reaction::compile_selected(&inputs, &target, &request, &selection).unwrap();
    assert!(matches!(
        plan.steps.last(),
        Some(Operation::SelectedGoldenGate { .. })
    ));
    assert_eq!(
        serde_json::to_value(&plan.steps[..4]).unwrap(),
        serde_json::to_value(&discovery.probe_plan.steps[..4]).unwrap()
    );
    let fresh = verification::verify(&plan, &target);
    assert_eq!(fresh.status, Status::ConsistentUnderAssumptions);
    assert_eq!(fresh.products.last().unwrap().sequence.as_str(), AB);
    assert!(
        fresh
            .findings
            .iter()
            .any(|finding| finding.contains("1 uncut carryover"))
    );
    assert_eq!(
        fresh.engine_version,
        collagene_core::verification::ENGINE_VERSION
    );

    // Matching the target cannot rescue a changed selection key, including when
    // an otherwise valid imported/deserialized plan supplies it.
    for directed in [true, false] {
        let mut forged: verification::DesignPlan =
            serde_json::from_slice(&serde_json::to_vec(&plan).unwrap()).unwrap();
        let Operation::SelectedGoldenGate {
            directed_sequence_key,
            physical_sequence_key,
            ..
        } = forged.steps.last_mut().unwrap()
        else {
            unreachable!()
        };
        if directed {
            *directed_sequence_key = "0".repeat(64);
        } else {
            *physical_sequence_key = "0".repeat(64);
        }
        let rejected = verification::verify(&forged, &target);
        assert_eq!(rejected.status, Status::Invalid);
        assert_eq!(
            rejected.products.len(),
            4,
            "only three declared-generation products and PCR survive a rejected final candidate"
        );
        assert!(
            rejected
                .findings
                .iter()
                .any(|finding| finding.contains("candidate keys disagree"))
        );
    }

    // A separately declared wrong target still rejects the exact selected AB path.
    let wrong_target = Molecule::dna("independent wrong target", A, Topology::Circular).unwrap();
    let mut wrong_selection = selection.clone();
    wrong_selection.context_revision =
        reaction::inspect_revision(&inputs, &wrong_target, &request).unwrap();
    let wrong_plan =
        reaction::compile_selected(&inputs, &wrong_target, &request, &wrong_selection).unwrap();
    let wrong = verification::verify(&wrong_plan, &wrong_target);
    assert_eq!(wrong.status, Status::Invalid);
    assert_eq!(wrong.products.last().unwrap().sequence.as_str(), AB);
}

#[test]
fn key_bound_compilation_preserves_explicit_choice_and_refuses_stale_context() {
    let (inputs, target, request) = fixture();
    let report = reaction::inspect(&inputs, &target, &request).unwrap();
    let pool = report.pool.as_ref().unwrap();
    let index = pool
        .candidates
        .iter()
        .position(|candidate| candidate.duplex.top.as_str() == A && candidate.restriction_site_free)
        .unwrap();
    let candidate = &pool.candidates[index];
    let mut selection = Selection {
        context_revision: report.context_revision.clone(),
        candidate_index: index,
        directed_sequence_key: candidate.directed_sequence_key.clone(),
        physical_sequence_key: candidate.physical_sequence_key.clone(),
    };
    let plan = reaction::compile_selected(&inputs, &target, &request, &selection).unwrap();
    let native = verification::verify(&plan, &target);
    assert_eq!(native.status, Status::Invalid);
    assert_eq!(
        native.products.last().unwrap().sequence.as_str(),
        A,
        "target does not choose the AB alternative"
    );
    selection.context_revision = "0".repeat(64);
    assert!(reaction::compile_selected(&inputs, &target, &request, &selection).is_err());
    selection.context_revision = report.context_revision;
    selection.candidate_index = request.options.max_products;
    assert!(reaction::compile_selected(&inputs, &target, &request, &selection).is_err());
    selection.candidate_index = index;
    selection.physical_sequence_key = "x".repeat(64);
    assert!(reaction::compile_selected(&inputs, &target, &request, &selection).is_err());
}

#[test]
fn keyed_reaction_protocol_and_junctions_keep_the_native_path_identity() {
    let (inputs, target, request) = fixture();
    let report = reaction::inspect(&inputs, &target, &request).unwrap();
    let pool = report.pool.as_ref().unwrap();
    let index = pool
        .candidates
        .iter()
        .position(|c| c.duplex.top.as_str() == AB && c.restriction_site_free)
        .unwrap();
    let candidate = &pool.candidates[index];
    let selection = Selection {
        context_revision: report.context_revision,
        candidate_index: index,
        directed_sequence_key: candidate.directed_sequence_key.clone(),
        physical_sequence_key: candidate.physical_sequence_key.clone(),
    };
    let plan = reaction::compile_selected(&inputs, &target, &request, &selection).unwrap();
    let html = verification::export_protocol_html(&plan, &target);
    assert!(html.contains(&selection.directed_sequence_key));
    assert!(html.contains(&selection.physical_sequence_key));
    let pdf = collagene_core::pdf::export_protocol_pdf(&plan, &target).unwrap();
    assert!(pdf.bytes.starts_with(b"%PDF-"));
    let junctions =
        collagene_core::assembly_junctions::inspect_plan(&plan, &target, &Default::default())
            .unwrap();
    assert!(junctions.complete);
    assert_eq!(
        junctions.verification.status,
        Status::ConsistentUnderAssumptions
    );
    assert!(!junctions.junctions.is_empty());
    let encoded = postcard::to_allocvec(plan.steps.last().unwrap()).unwrap();
    assert_eq!(
        encoded[0], 18,
        "new op is appended; existing discriminants stay unchanged"
    );
}
