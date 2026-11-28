//! Supplied oligos prepare a changed insert from unchanged stock before assembly.
use collagene_core::{
    analysis::Primer,
    assembly_diagnostics::{self, Stage},
    assembly_recipe::{self as a, FragmentSpec, Preparation, RunPreparation},
    chemistry::EndTreatment,
    cloning_outputs, ligation_recipe as lig, overlap_extension,
    overlap_extension_pcr::{PrimerPair, Request},
    overlap_profiles::ProfileOptions,
    run_lineage::{EventBinding, InputEventBinding},
    run_material::{self, Selection},
    verification::{self, DesignPlan, Operation, PlanInput, Source, Status},
    workspace::Workspace,
};
use collagene_model::{Molecule, Topology};
const ORIGINAL: &str = "AGTCCGATGCTACGTAGGACTTCAGACCTGTAGATCGTCA";
const EDITED: &str = "AGTCCGATGCTACGTATTTTCAGACCTGTAGATCGTCA";
const FORWARD_TARGET: &str = "CACGTTCGAGTCCGATGCTACGTATTTTCAGACCTGTAGATCGTCA";
const REVERSE_TARGET: &str = "CACGTTCGTGACGATCTACAGGTCTGAAAATACGTAGCATCGGACT";
fn input(id: &str, text: &str) -> PlanInput {
    let molecule = Molecule::dna(id, text, Topology::Linear).unwrap();
    PlanInput {
        id: id.into(),
        revision: molecule.revision(),
        molecule,
    }
}
fn pair(f: &str, ft: &str, r: &str, rt: &str) -> PrimerPair {
    PrimerPair {
        forward: Primer::new("Forward", f, ft).unwrap(),
        reverse: Primer::new("Reverse", r, rt).unwrap(),
        template_strand: collagene_core::chemistry::Strand::Top,
    }
}
fn request() -> Request {
    Request {
        left: pair("AGTCCGAT", "", "TACGTAGC", "GGTCTGAAAA"),
        right: pair("TTCAGACC", "GCTACGTATT", "TGACGATC", ""),
        outer: pair("AGTCCGAT", "", "TGACGATC", ""),
        extension: overlap_extension::Options {
            overlap: 18,
            min_paired_bases: 8,
            ..Default::default()
        },
    }
}
fn row(id: &str, preparation: Preparation, reverse: bool) -> FragmentSpec {
    FragmentSpec {
        generation: None,
        role: format!("{id}-role"),
        input_id: id.into(),
        preparation,
        reverse,
    }
}
fn overlap(reverse: bool) -> (Vec<PlanInput>, Molecule, a::Recipe) {
    let (vector, reverse_primer, target) = if reverse {
        ("CACGTTCGTGACGATC", "GATCGTCA", REVERSE_TARGET)
    } else {
        ("CACGTTCGAGTCCGAT", "ATCGGACT", FORWARD_TARGET)
    };
    (
        vec![input("vector", vector), input("insert", ORIGINAL)],
        Molecule::dna(
            "Independent multipart GGAC to TT intent",
            target,
            Topology::Linear,
        )
        .unwrap(),
        a::Recipe {
            name: "Edited insert and independent vector".into(),
            fragments: vec![
                row(
                    "vector",
                    Preparation::ExactPcr {
                        forward: Primer::new("Vector forward", "CACGTTCG", "").unwrap(),
                        reverse: Primer::new("Vector reverse", reverse_primer, "").unwrap(),
                        template_strand: collagene_core::chemistry::Strand::Top,
                    },
                    false,
                ),
                row(
                    "insert",
                    Preparation::OverlapExtensionPcr { request: request() },
                    reverse,
                ),
            ],
            overlaps: vec![8],
            circular: false,
            options: ProfileOptions::default(),
        },
    )
}
#[test]
fn literal_edited_insert_then_multipart_overlap_retains_both_branches_and_outer_pcr() {
    for reverse in [false, true] {
        let (inputs, target, recipe) = overlap(reverse);
        let originals = serde_json::to_value(&inputs).unwrap();
        let preview = a::preview(&inputs, &target, &recipe).unwrap();
        assert_eq!(
            preview.verification.status,
            Status::ConsistentUnderAssumptions,
            "{:?}",
            preview.verification.findings
        );
        assert_eq!(preview.compiled.assembly_step, 5);
        assert_eq!(
            preview
                .verification
                .products
                .iter()
                .map(|p| p.sequence.len())
                .collect::<Vec<_>>(),
            vec![16, 26, 30, 38, 38, 46]
        );
        assert_eq!(
            preview.verification.products[1].sequence.as_str(),
            "AGTCCGATGCTACGTATTTTCAGACC"
        );
        assert_eq!(
            preview.verification.products[2].sequence.as_str(),
            "GCTACGTATTTTCAGACCTGTAGATCGTCA"
        );
        assert_eq!(preview.prepared[1].sequence.as_str(), EDITED);
        assert_eq!(
            preview.verification.products[5].sequence.as_str(),
            target.sequence.as_str()
        );
        let row = &preview.compiled.rows[1];
        assert_eq!(row.preparation_step, Some(1));
        assert!(row.repair_step.is_none());
        assert!(matches!(row.prepared_source, Source::Step(4)));
        assert_eq!(row.reverse, reverse);
        assert!(
            matches!(&preview.compiled.plan.steps[1], Operation::Pcr {template: Source::Input(id), ..} if id == "insert")
        );
        assert!(
            matches!(&preview.compiled.plan.steps[2], Operation::Pcr {template: Source::Input(id), ..} if id == "insert")
        );
        assert!(matches!(
            &preview.compiled.plan.steps[4],
            Operation::Pcr {
                template: Source::Step(3),
                ..
            }
        ));
        let diagnostics = assembly_diagnostics::overlap(&preview).unwrap();
        assert_eq!(diagnostics.rows[1].total_steps, 4);
        assert_eq!(diagnostics.rows[1].completed_steps, 4);
        assert!(diagnostics.failure.is_none());
        preview.validate_binding(&preview.context_revision).unwrap();
        assert_eq!(serde_json::to_value(&inputs).unwrap(), originals);
    }
}
#[test]
fn actual_prepared_duplex_can_be_phosphorylated_and_ligated_without_reference_materialization() {
    let (mut inputs, target, overlap) = overlap(false);
    inputs[0] = input("vector", "CACGTTCG");
    let recipe = lig::Recipe {
        name: "Physical edited insert ligation".into(),
        fragments: vec![
            lig::Fragment {
                source: row(
                    "vector",
                    Preparation::ExactPcr {
                        forward: Primer::new("Vector F", "CACG", "").unwrap(),
                        reverse: Primer::new("Vector R", "CGAA", "").unwrap(),
                        template_strand: collagene_core::chemistry::Strand::Top,
                    },
                    false,
                ),
                treatment: Some(EndTreatment::Phosphorylate),
            },
            lig::Fragment {
                source: overlap.fragments[1].clone(),
                treatment: Some(EndTreatment::Phosphorylate),
            },
        ],
        circular: false,
        allow_nicked: false,
    };
    let preview = lig::preview(
        &inputs,
        &target,
        &recipe,
        collagene_core::assembly_recipe::CleavagePolicy::AssumedComplete,
    )
    .unwrap();
    assert_eq!(
        preview.verification.status,
        Status::ConsistentUnderAssumptions,
        "{:?}",
        preview.verification.findings
    );
    assert_eq!(preview.compiled.ligation_step, 7);
    assert_eq!(
        preview
            .verification
            .products
            .last()
            .unwrap()
            .sequence
            .as_str(),
        FORWARD_TARGET
    );
    assert!(preview.physical.iter().all(Option::is_some));
    assert!(preview.product_physical.as_ref().unwrap().nicks.is_empty());
    let diagnostics = assembly_diagnostics::ligation(&preview).unwrap();
    assert_eq!(diagnostics.rows[1].total_steps, 5);
    assert!(diagnostics.failure.is_none());
}
fn binding(w: &Workspace, index: usize) -> EventBinding {
    let d = &w.documents[index];
    EventBinding {
        document_id: d.id.clone(),
        event_id: d.head_event.clone().unwrap(),
        revision: d.head.clone(),
    }
}
#[test]
fn saved_original_is_replayed_before_both_branches_without_changing_its_origin_target() {
    let source = input("original", ORIGINAL);
    let original_target = Molecule::dna(
        "Independent original PCR target",
        ORIGINAL,
        Topology::Linear,
    )
    .unwrap();
    let plan = DesignPlan {
        name: "Original physical template".into(),
        target_revision: original_target.revision(),
        inputs: vec![source.clone()],
        steps: vec![Operation::Pcr {
            template: Source::Input("original".into()),
            forward: Primer::new("F", "AGTCCGAT", "").unwrap(),
            reverse: Primer::new("R", "TGACGATC", "").unwrap(),
            template_strand: collagene_core::chemistry::Strand::Top,
        }],
    };
    let w = Workspace::from_documents("Origin", vec![source.molecule, original_target], 1).unwrap();
    let commit = cloning_outputs::finalize(
        &w,
        &plan,
        binding(&w, 1),
        vec![InputEventBinding {
            input_id: "original".into(),
            event: binding(&w, 0),
        }],
        &cloning_outputs::Options {
            name: "Saved original stock PCR".into(),
            show_prepared: false,
            save_primers: false,
            allow_incomplete: false,
        },
    )
    .unwrap();
    let event = commit.receipt.construct.clone();
    let origin = run_material::capsule(
        &commit.workspace,
        &Selection {
            run_id: commit.receipt.run.run_id.clone(),
            step_index: 0,
            event: event.clone(),
        },
    )
    .unwrap();
    let immutable_origin = serde_json::to_value(&origin).unwrap();
    let (mut inputs, target, mut recipe) = overlap(false);
    inputs[1].molecule = commit.workspace.revisions[&event.revision].clone();
    inputs[1].revision = event.revision;
    recipe.fragments[1].preparation = Preparation::FromRun {
        origin,
        preparation: RunPreparation::OverlapExtensionPcr { request: request() },
    };
    let preview = a::preview(&inputs, &target, &recipe).unwrap();
    assert_eq!(
        preview.verification.status,
        Status::ConsistentUnderAssumptions,
        "{:?}",
        preview.verification.findings
    );
    assert_eq!(preview.compiled.assembly_step, 6);
    assert!(matches!(
        &preview.compiled.plan.steps[2],
        Operation::Pcr {
            template: Source::Step(1),
            ..
        }
    ));
    assert!(matches!(
        &preview.compiled.plan.steps[3],
        Operation::Pcr {
            template: Source::Step(1),
            ..
        }
    ));
    let Operation::ReconstructRunProduct { origin, .. } = &preview.compiled.plan.steps[1] else {
        panic!("saved reconstruction prefix");
    };
    assert_eq!(serde_json::to_value(origin).unwrap(), immutable_origin);
    assert_eq!(origin.target.sequence.as_str(), ORIGINAL);
    assert_ne!(origin.target.revision(), target.revision());
    assert_eq!(
        assembly_diagnostics::overlap(&preview).unwrap().rows[1].total_steps,
        5
    );
    assert_eq!(
        preview
            .verification
            .products
            .last()
            .unwrap()
            .sequence
            .as_str(),
        FORWARD_TARGET
    );
}
fn rebind_plan(preview: &mut a::Preview) {
    preview.verification.plan_revision = verification::plan_revision(&preview.compiled.plan);
}
#[test]
fn graph_correspondence_rejects_wrong_branch_pair_outer_source_skips_and_shared_ownership() {
    let (inputs, target, recipe) = overlap(false);
    let preview = a::preview(&inputs, &target, &recipe).unwrap();
    for kind in 0..6 {
        let mut bad = preview.clone();
        match kind {
            0 => {
                if let Operation::Pcr { template, .. } = &mut bad.compiled.plan.steps[2] {
                    *template = Source::Step(1);
                }
            }
            1 => {
                if let Operation::OverlapExtension { right, .. } = &mut bad.compiled.plan.steps[3] {
                    right.source = Source::Step(1);
                }
            }
            2 => {
                if let Operation::OverlapExtension { left, .. } = &mut bad.compiled.plan.steps[3] {
                    left.reverse = true;
                }
            }
            3 => {
                if let Operation::Pcr { template, .. } = &mut bad.compiled.plan.steps[4] {
                    *template = Source::Step(1);
                }
            }
            4 => bad.compiled.rows[1].prepared_source = Source::Step(3),
            _ => bad.compiled.rows[1].preparation_step = Some(0),
        }
        rebind_plan(&mut bad);
        assert!(bad.validate().is_err(), "tampered graph {kind}");
    }
    let mut wrong = target;
    wrong.sequence =
        collagene_model::Sequence::dna("CACGTTCGAGTCCGATGCTACGTAGGACTTCAGACCTGTAGATCGTCA").unwrap();
    let mismatch = a::preview(&inputs, &wrong, &recipe).unwrap();
    assert_eq!(mismatch.verification.status, Status::Invalid);
    assert_eq!(
        mismatch
            .verification
            .products
            .last()
            .unwrap()
            .sequence
            .as_str(),
        FORWARD_TARGET
    );
}
#[test]
fn invalid_oligos_and_search_limits_fail_before_source_copy_and_failure_points_to_extension() {
    let (inputs, target, mut recipe) = overlap(false);
    for kind in 0..5 {
        let mut request = request();
        match kind {
            0 => request.outer.forward.tail = "A".repeat(1001),
            1 => request.left.reverse.name = "bad\nname".into(),
            2 => request.right.forward.tail = "NNNN".into(),
            3 => request.extension.max_candidates = 0,
            _ => request.extension.max_comparisons = usize::MAX,
        }
        assert!(request.validate_shape().is_err());
        recipe.fragments[1].preparation = Preparation::OverlapExtensionPcr { request };
        assert!(
            a::source_preflight(
                &inputs
                    .iter()
                    .map(|i| (i.id.as_str(), &i.molecule))
                    .collect::<Vec<_>>(),
                &target,
                &recipe
            )
            .is_err()
        );
    }
    let mut bounded = request();
    bounded.extension.max_comparisons = bounded.extension.overlap;
    recipe.fragments[1].preparation = Preparation::OverlapExtensionPcr { request: bounded };
    let uncertain = a::preview(&inputs, &target, &recipe).unwrap();
    assert_eq!(uncertain.verification.status, Status::Incomplete);
    assert_eq!(uncertain.verification.products.len(), 6);
    assert_eq!(
        uncertain.verification.products[5].sequence.as_str(),
        FORWARD_TARGET
    );
    assert!(
        assembly_diagnostics::overlap(&uncertain)
            .unwrap()
            .failure
            .is_none()
    );
    let mut request = request();
    request.extension.overlap = 17;
    recipe.fragments[1].preparation = Preparation::OverlapExtensionPcr { request };
    let failed = a::preview(&inputs, &target, &recipe).unwrap();
    assert_ne!(
        failed.verification.status,
        Status::ConsistentUnderAssumptions
    );
    let location = assembly_diagnostics::overlap(&failed)
        .unwrap()
        .failure
        .unwrap();
    assert_eq!(location.stage, Stage::OverlapExtension);
    assert_eq!(location.row, Some(1));
    assert_eq!(location.step, 3);
}
