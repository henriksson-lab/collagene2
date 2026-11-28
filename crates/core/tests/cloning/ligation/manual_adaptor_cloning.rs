//! A hand-specified adaptor/backbone construct, built through manual graph edits.
use collagene_core::{
    chemistry::{self, AnnealOptions, Terminus},
    documents::{
        transfer::Candidate,
        workspace::{self, Workspace},
    },
    material::lineage::{self as run_lineage, EventBinding, InputEventBinding},
    planning::{
        editing as edit,
        plan::{DesignPlan, Operation, OrientedSource, PlanInput, Source},
    },
    verification::{self as v, evidence::Status},
};
use collagene_model::{Molecule, Topology};
const BACKBONE: &str = "GGAATTCTTACG";
const TOP: &str = "AATTGCTACG";
const BOTTOM: &str = "AATTCGTAGC";
const TARGET: &str = "AATTCTTACGGGAATTGCTACG";
fn input(id: &str, dna: &str, topology: Topology) -> PlanInput {
    let molecule = Molecule::dna(id, dna, topology).unwrap();
    PlanInput {
        id: id.into(),
        revision: molecule.revision(),
        molecule,
    }
}
fn fixture() -> (DesignPlan, Molecule) {
    let inputs = vec![
        input("adaptor-top", TOP, Topology::Linear),
        input("adaptor-bottom", BOTTOM, Topology::Linear),
        input("original-backbone", BACKBONE, Topology::Circular),
    ];
    let ids = inputs.iter().map(|i| i.id.clone()).collect::<Vec<_>>();
    let initial = vec![Operation::AnnealOligos {
        top: Source::Input("adaptor-top".into()),
        bottom: Source::Input("adaptor-bottom".into()),
        options: AnnealOptions {
            min_paired_bases: 6,
            top_chemistry: Terminus::CLEAVED,
            bottom_chemistry: Terminus::CLEAVED,
            ..Default::default()
        },
        candidate_index: 0,
    }];
    let steps = edit::inserted(
        &ids,
        &initial,
        1,
        Operation::Digest {
            template: Source::Input("original-backbone".into()),
            enzymes: vec!["EcoRI".into()],
            fragment_index: 0,
        },
    )
    .unwrap();
    let steps = edit::inserted(
        &ids,
        &steps,
        2,
        Operation::Ligate {
            fragments: vec![
                OrientedSource {
                    source: Source::Step(1),
                    reverse: false,
                },
                OrientedSource {
                    source: Source::Step(0),
                    reverse: false,
                },
            ],
            circular: true,
            allow_nicked: false,
            options: Default::default(),
        },
    )
    .unwrap();
    let target = Molecule::dna(
        "Independently specified22bp circle",
        TARGET,
        Topology::Circular,
    )
    .unwrap();
    (
        DesignPlan {
            standard: None,
            name: "Annealed adaptor into a separately digested backbone".into(),
            target_revision: target.revision(),
            inputs,
            steps,
        },
        target,
    )
}
#[test]
fn literal_adaptor_annealing_and_actual_backbone_digest_form_the_frozen_circle() {
    let (plan, target) = fixture();
    let candidate = Candidate::new(plan.clone());
    candidate.validate().unwrap();
    let r = v::engine::verify(&plan, &target);
    assert_eq!(
        r.status,
        Status::ConsistentUnderAssumptions,
        "{:?}",
        r.findings
    );
    assert_eq!(r.reconstructed_count(), 3);
    assert_eq!(
        r.product(2)
            .expect("selected step must reconstruct")
            .sequence
            .as_str(),
        TARGET
    );
    let annealed = chemistry::anneal_oligos(
        &plan.inputs[0].molecule.sequence,
        &plan.inputs[1].molecule.sequence,
        &match &plan.steps[0] {
            Operation::AnnealOligos { options, .. } => options.clone(),
            _ => unreachable!(),
        },
    )
    .unwrap();
    assert!(annealed.complete);
    assert_eq!(annealed.candidates.len(), 1);
    let insert = &annealed.candidates[0].duplex;
    assert_eq!(insert.top.as_str(), TOP);
    assert_eq!(insert.bottom.as_str(), BOTTOM);
    assert_eq!(insert.bottom_start, 4);
    assert_eq!(annealed.candidates[0].paired_bases, 6);
    let (left, right) = insert.ends().unwrap();
    assert_eq!(left.overhang, "AATT");
    assert_eq!(right.overhang, "AATT");
    let restriction = chemistry::digest(
        &chemistry::Duplex::from_molecule(&plan.inputs[2].molecule).unwrap(),
        &[collagene_core::sequence::analysis::Enzyme {
            name: "Literal EcoRI".into(),
            recognition: "GAATTC".into(),
            top_cut: 1,
            bottom_cut: 5,
        }],
    )
    .unwrap();
    assert_eq!(restriction.len(), 1);
    assert_eq!(restriction[0].top_start, 2);
    assert_eq!(restriction[0].bottom_start, 6);
    assert_eq!(restriction[0].duplex.top.as_str(), "AATTCTTACGGG");
    assert_eq!(
        TARGET.len() * 2,
        BACKBONE.len() * 2 + TOP.len() + BOTTOM.len()
    );
    let summary: chemistry::DuplexSummary = serde_json::from_str(
        &r.product(2)
            .expect("selected step must reconstruct")
            .metadata["duplex_model"],
    )
    .unwrap();
    assert_eq!(summary.topology, Topology::Circular);
    assert!(summary.nicks.is_empty());
    let html = v::output::export_protocol_html(&plan, &target);
    assert!(html.contains(TOP));
    assert!(html.contains(BOTTOM));
    assert!(html.contains("EcoRI"));
    let pdf = collagene_core::output::pdf::export_protocol_pdf(&plan, &target).unwrap();
    assert_eq!(pdf.status, Status::ConsistentUnderAssumptions);
    assert!(pdf.bytes.starts_with(b"%PDF-"));
}
#[test]
fn dependency_preserving_edits_reverify_while_invalid_order_or_target_never_inherits_success() {
    let (mut plan, target) = fixture();
    let ids = plan.inputs.iter().map(|i| i.id.clone()).collect::<Vec<_>>();
    let before = serde_json::to_vec(&plan.steps).unwrap();
    assert!(edit::removed(&ids, &plan.steps, 0).is_err());
    assert!(edit::moved(&ids, &plan.steps, 0, 2).is_err());
    assert_eq!(serde_json::to_vec(&plan.steps).unwrap(), before);
    plan.steps = edit::moved(&ids, &plan.steps, 1, 0).unwrap();
    let r = v::engine::verify(&plan, &target);
    assert_eq!(
        r.status,
        Status::ConsistentUnderAssumptions,
        "{:?}",
        r.findings
    );
    assert_eq!(
        r.product(2)
            .expect("selected step must reconstruct")
            .sequence
            .as_str(),
        TARGET
    );
    let wrong = Molecule::dna(
        "Wrong independent target",
        "AATTCTTACGGGAATTGCTACC",
        Topology::Circular,
    )
    .unwrap();
    plan.target_revision = wrong.revision();
    assert_eq!(v::engine::verify(&plan, &wrong).status, Status::Invalid);
    plan.target_revision = target.revision();
    if let Operation::Digest { enzymes, .. } = &mut plan.steps[0] {
        *enzymes = vec!["HindIII".into()];
    }
    assert_ne!(
        v::engine::verify(&plan, &target).status,
        Status::ConsistentUnderAssumptions
    );
}
fn event(w: &Workspace, i: usize) -> EventBinding {
    let d = &w.documents[i];
    EventBinding {
        document_id: d.id.clone(),
        event_id: d.head_event.clone().unwrap(),
        revision: d.head.clone(),
    }
}
#[test]
fn edited_candidate_roundtrip_and_retained_run_preserve_all_original_source_events() {
    let (plan, target) = fixture();
    let value = Candidate::new(plan.clone());
    let bytes = collagene_core::documents::transfer::save_candidate(&value).unwrap();
    let read = collagene_core::documents::transfer::load_candidate(&bytes).unwrap();
    assert_eq!(read.plan_revision, value.plan_revision);
    let mut docs = plan
        .inputs
        .iter()
        .map(|i| i.molecule.clone())
        .collect::<Vec<_>>();
    docs.push(target.clone());
    let original = Workspace::from_documents("Manual adaptor design", docs, 3).unwrap();
    let bindings = plan
        .inputs
        .iter()
        .enumerate()
        .map(|(i, input)| InputEventBinding {
            input_id: input.id.clone(),
            event: event(&original, i),
        })
        .collect();
    let saved =
        run_lineage::record_run(&original, &plan, event(&original, 3), bindings, false).unwrap();
    let archived =
        workspace::load_workspace(&workspace::save_workspace(&saved.workspace).unwrap()).unwrap();
    let record = run_lineage::load_run(&archived, &saved.receipt.run_id).unwrap();
    assert_eq!(record.inputs.len(), 3);
    assert_eq!(record.products.len(), 3);
    let final_event = archived
        .history
        .iter()
        .find(|e| e.id == record.products[2].event_id)
        .unwrap();
    assert!(final_event.parents.contains(&record.products[0].event_id));
    assert!(final_event.parents.contains(&record.products[1].event_id));
    assert!(!final_event.parents.contains(&record.target.event_id));
    let replay = run_lineage::replay_run(&archived, &record.run_id).unwrap();
    assert_eq!(
        replay.current_report.status,
        Status::ConsistentUnderAssumptions
    );
    assert!(replay.products_match);
    // The original documents remain input records; no reference preview is
    // substituted for either of the fresh private duplex substrates.
    for i in 0..3 {
        assert_eq!(archived.documents[i].head, original.documents[i].head);
    }
}
