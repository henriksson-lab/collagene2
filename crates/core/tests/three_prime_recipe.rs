//! Independent ordinary-cloning declarations; joined diagnostic DNA never becomes material.
use collagene_core::{
    Feature, Location, Molecule, Topology,
    analysis::Primer,
    assembly_recipe::{self as recipe, FragmentSpec, Preparation},
    chemistry::{self, EndRepairMethod, FivePrime, Terminus, ThreePrime},
    homology_annealing::{Options, Report},
    run_lineage::{self, EventBinding, InputEventBinding},
    run_material,
    three_prime_resection::Cuts,
    verification::{
        self, DesignPlan, Operation, OrientedSource, PhysicalProcessingEvidence, PlanInput, Source,
        Status,
    },
    workspace::{self, Workspace},
};
const LEFT: &str = "AGTCCGATGCTACGTAGGAC";
const RIGHT: &str = "TTCAGACCTGTAGATCGTCA";
const CIRCLE: &str = "AGTCCGATGCTACGTAGGACTTCAGACCTGTAGATCGTCA";
fn input(id: &str, molecule: Molecule) -> PlanInput {
    PlanInput {
        id: id.into(),
        revision: molecule.revision(),
        molecule,
    }
}
fn fixture(circular: bool) -> (Vec<PlanInput>, recipe::Recipe, Options) {
    let mut left = Molecule::dna("Original left template", LEFT, Topology::Linear).unwrap();
    left.features.push(Feature {
        kind: "misc_feature".into(),
        location: Location::range(0, 8),
        qualifiers: vec![("label".into(), Some("Original left calls".into()))],
        needs_review: false,
    });
    let right = Molecule::dna("Original right template", RIGHT, Topology::Linear).unwrap();
    let inputs = vec![input("left", left), input("right", right)];
    let fragments = vec![
        FragmentSpec {
            generation: None,
            role: "left-role".into(),
            input_id: "left".into(),
            reverse: false,
            preparation: Preparation::ExactPcr {
                forward: Primer::new("Left F", "AGTCCGAT", "").unwrap(),
                reverse: Primer::new("Left R + downstream homology", "GTCCTACG", "GGTCTGAA")
                    .unwrap(),
                template_strand: collagene_core::chemistry::Strand::Top,
            },
        },
        FragmentSpec {
            generation: None,
            role: "right-role".into(),
            input_id: "right".into(),
            reverse: false,
            preparation: Preparation::ExactPcr {
                forward: Primer::new("Right F", "TTCAGACC", "").unwrap(),
                reverse: Primer::new(
                    "Right R + circle closure",
                    "TGACGATC",
                    if circular { "ATCGGACT" } else { "" },
                )
                .unwrap(),
                template_strand: collagene_core::chemistry::Strand::Top,
            },
        },
    ];
    let overlaps = if circular { vec![8, 8] } else { vec![8] };
    let mut profile = collagene_core::overlap_profiles::ProfileOptions::default();
    profile.scan_repeated_overlaps = false;
    profile.method = collagene_core::overlap_profiles::OverlapMethod::InFusion;
    let request = recipe::Recipe {
        name: "Declared physical homology route".into(),
        fragments,
        overlaps: overlaps.clone(),
        circular,
        options: profile,
    };
    let policy = Options {
        overlaps,
        topology: if circular {
            Topology::Circular
        } else {
            Topology::Linear
        },
        cuts: if circular {
            vec![Cuts { left: 8, right: 8 }, Cuts { left: 8, right: 8 }]
        } else {
            vec![Cuts { left: 0, right: 8 }, Cuts { left: 8, right: 0 }]
        },
        min_paired_bases: 4,
        max_comparisons: 100_000,
        max_candidates: 20,
    };
    (inputs, request, policy)
}
fn binding(w: &Workspace, index: usize) -> EventBinding {
    let d = &w.documents[index];
    EventBinding {
        document_id: d.id.clone(),
        event_id: d.head_event.clone().unwrap(),
        revision: d.head.clone(),
    }
}
fn annealed(p: &recipe::Preview) -> &Report {
    p.verification
        .physical_processing_reports
        .iter()
        .find_map(|s| match &s.evidence {
            PhysicalProcessingEvidence::HomologyAnnealing(r) => Some(r),
            _ => None,
        })
        .unwrap()
}
#[test]
fn pcr_hydroxyl_fragments_anneal_to_independent_circle_but_every_bond_stays_unsealed() {
    let (inputs, request, policy) = fixture(true);
    let target = Molecule::dna("Independent circle40", CIRCLE, Topology::Circular).unwrap();
    let p = recipe::preview_annealing(&inputs, &target, &request, &policy).unwrap();
    p.validate_binding(
        &recipe::annealing_preview_revision(&inputs, &target, &request, &policy).unwrap(),
    )
    .unwrap();
    assert_eq!(
        p.verification.status,
        Status::Incomplete,
        "{:?}",
        p.verification.findings
    );
    assert!(matches!(
        p.compiled.plan.steps[2],
        Operation::HomologyAnnealing { .. }
    ));
    assert_eq!(
        p.verification.products.last().unwrap().sequence.as_str(),
        CIRCLE
    );
    let r = annealed(&p);
    r.validate_documentary().unwrap();
    assert_eq!(r.segments.len(), 4);
    assert_eq!(r.paired_diagnostic.as_ref().unwrap().nicks.len(), 4);
    assert!(
        r.junctions
            .iter()
            .all(|j| j.top.five_prime == FivePrime::Hydroxyl
                && j.bottom.five_prime == FivePrime::Hydroxyl)
    );
    assert_eq!(r.junctions[0].top.start, 20);
    assert_eq!(r.junctions[0].bottom.start, 28);
    assert_eq!(r.junctions[1].top.start, 40);
    assert_eq!(r.junctions[1].bottom.start, 48);
    let product = p.verification.products.last().unwrap();
    assert!(chemistry::has_unmapped_physical_evidence(product));
    assert!(chemistry::Duplex::from_molecule(product).is_err());
    let junctions = collagene_core::assembly_junctions::inspect_plan(
        &p.compiled.plan,
        &target,
        &Default::default(),
    )
    .unwrap();
    assert!(matches!(
        junctions.method,
        Some(collagene_core::assembly_junctions::Method::HomologyAnnealing)
    ));
}
#[test]
fn independent_wrong_target_does_not_select_or_repair_a_different_physical_path() {
    let (inputs, request, policy) = fixture(true);
    let right = Molecule::dna("Independent circle40", CIRCLE, Topology::Circular).unwrap();
    let wrong = Molecule::dna(
        "Different declared target",
        "AGTCCGATGCTACGTAGGACTTCAGACCTGTAGATCGTCC",
        Topology::Circular,
    )
    .unwrap();
    let a = recipe::preview_annealing(&inputs, &right, &request, &policy).unwrap();
    let b = recipe::preview_annealing(&inputs, &wrong, &request, &policy).unwrap();
    assert_eq!(b.verification.status, Status::Invalid);
    assert_eq!(annealed(&a), annealed(&b));
    assert_eq!(
        a.verification.products.last().unwrap().sequence,
        b.verification.products.last().unwrap().sequence
    );
    assert_ne!(a.context_revision, b.context_revision);
}
#[test]
fn excess_resection_reports_actual_gaps_and_comparison_limits_remain_incomplete() {
    let (inputs, request, mut policy) = fixture(true);
    let target = Molecule::dna("Independent circle40", CIRCLE, Topology::Circular).unwrap();
    policy.cuts[0].right = 10;
    policy.cuts[1].left = 11;
    let p = recipe::preview_annealing(&inputs, &target, &request, &policy).unwrap();
    let r = annealed(&p);
    assert_eq!(p.verification.status, Status::Incomplete);
    assert!(r.paired_diagnostic.is_none());
    assert_eq!((r.junctions[0].top.start, r.junctions[0].top.end), (18, 20));
    assert_eq!(
        (r.junctions[0].bottom.start, r.junctions[0].bottom.end),
        (28, 31)
    );
    assert_eq!(
        p.verification.products.last().unwrap().sequence.as_str(),
        CIRCLE
    );
    let (_, _, mut policy) = fixture(true);
    policy.max_comparisons = 16;
    let p = recipe::preview_annealing(&inputs, &target, &request, &policy).unwrap();
    let r = annealed(&p);
    assert_eq!(r.comparisons, 16);
    assert!(!r.terminal_search_complete);
    assert_eq!(p.verification.status, Status::Incomplete);
    let mut forged = p.clone();
    let step = forged
        .verification
        .physical_processing_reports
        .iter_mut()
        .find_map(|s| match &mut s.evidence {
            PhysicalProcessingEvidence::HomologyAnnealing(r) => Some(r),
            _ => None,
        })
        .unwrap();
    step.junctions[0].top.five_prime = FivePrime::Phosphate;
    assert!(forged.validate().is_err());
}
#[test]
fn saved_annealed_reference_cannot_supply_pcr_resection_or_saved_material_recovery() {
    let (inputs, request, policy) = fixture(true);
    let target = Molecule::dna("Independent circle40", CIRCLE, Topology::Circular).unwrap();
    let p = recipe::preview_annealing(&inputs, &target, &request, &policy).unwrap();
    for operation in [
        Operation::Pcr {
            template: Source::Step(2),
            forward: Primer::new("F", "AGTCCGAT", "").unwrap(),
            reverse: Primer::new("R", "TGACGATC", "").unwrap(),
            template_strand: collagene_core::chemistry::Strand::Top,
        },
        Operation::ResectThreePrime {
            source: Source::Step(2),
            cuts: Cuts { left: 1, right: 1 },
        },
    ] {
        let mut plan = p.compiled.plan.clone();
        plan.steps.push(operation);
        let result = verification::verify(&plan, &target);
        assert_ne!(result.status, Status::ConsistentUnderAssumptions);
        assert_eq!(result.products.len(), 3);
    }
    let w = Workspace::from_documents(
        "Exact original events",
        vec![
            inputs[0].molecule.clone(),
            inputs[1].molecule.clone(),
            target,
        ],
        2,
    )
    .unwrap();
    let saved = run_lineage::record_run(
        &w,
        &p.compiled.plan,
        binding(&w, 2),
        vec![
            InputEventBinding {
                input_id: "left".into(),
                event: binding(&w, 0),
            },
            InputEventBinding {
                input_id: "right".into(),
                event: binding(&w, 1),
            },
        ],
        true,
    )
    .unwrap();
    let output = saved
        .receipt
        .products
        .iter()
        .find(|p| p.step_index == 2)
        .unwrap();
    let selection = run_material::Selection {
        run_id: saved.receipt.run_id.clone(),
        step_index: 2,
        event: EventBinding {
            document_id: output.document_id.clone(),
            event_id: output.event_id.clone(),
            revision: output.revision.clone(),
        },
    };
    let archive = workspace::save_workspace(&saved.workspace).unwrap();
    let loaded = workspace::load_workspace(&archive).unwrap();
    assert!(run_material::capsule(&loaded, &selection).is_err());
    let evidence = collagene_core::physical_run_evidence::load(&loaded, &saved.receipt.run_id)
        .unwrap()
        .unwrap();
    assert!(
        evidence
            .steps
            .iter()
            .any(|s| matches!(s.evidence, PhysicalProcessingEvidence::HomologyAnnealing(_)))
    );
}
fn resection_plan() -> (DesignPlan, Molecule, Molecule) {
    let mut source = Molecule::dna("Original forty", CIRCLE, Topology::Linear).unwrap();
    source.features.push(Feature {
        kind: "misc_feature".into(),
        location: Location::range(35, 40),
        qualifiers: vec![(
            "label".into(),
            Some("Partly resected source annotation".into()),
        )],
        needs_review: false,
    });
    let target = Molecule::dna("Independent shortened37", &CIRCLE[..37], Topology::Linear).unwrap();
    let plan = DesignPlan {
        name: "Explicit actual3′ shortening".into(),
        target_revision: target.revision(),
        inputs: vec![input("stock", source.clone())],
        steps: vec![
            Operation::DeclareSourceTermini {
                template: Source::Input("stock".into()),
                left: Terminus::PCR,
                right: Terminus::PCR,
            },
            Operation::ResectThreePrime {
                source: Source::Step(0),
                cuts: Cuts { left: 2, right: 3 },
            },
        ],
    };
    (plan, source, target)
}
#[test]
fn actual_resected_operand_retains_extruded_bottom_and_clips_top_annotation_then_saved_origin_can_be_filled()
 {
    let (plan, source, target) = resection_plan();
    let report = verification::verify(&plan, &target);
    assert_eq!(
        report.status,
        Status::ConsistentUnderAssumptions,
        "{:?}",
        report.findings
    );
    let actual = &report.products[1];
    assert_eq!(actual.sequence.as_str(), &CIRCLE[..37]);
    assert!(actual.features.iter().any(|f| {
        f.needs_review
            && f.location
                .spans()
                .iter()
                .any(|(a, b, _)| *a == 35 && *b == 37)
    }));
    let r = report
        .physical_processing_reports
        .iter()
        .find_map(|s| match &s.evidence {
            PhysicalProcessingEvidence::ThreePrimeResection(r) => Some(r),
            _ => None,
        })
        .unwrap();
    assert_eq!(
        (
            r.product.top.len(),
            r.product.bottom.len(),
            r.product.bottom_start
        ),
        (37, 38, 2)
    );
    assert_eq!(r.product.left.five_prime, FivePrime::Hydroxyl);
    assert_eq!(r.product.right.three_prime, ThreePrime::Hydroxyl);
    r.validate_documentary().unwrap();
    let w = Workspace::from_documents(
        "Literal original source and target",
        vec![source, target],
        1,
    )
    .unwrap();
    let saved = run_lineage::record_run(
        &w,
        &plan,
        binding(&w, 1),
        vec![InputEventBinding {
            input_id: "stock".into(),
            event: binding(&w, 0),
        }],
        false,
    )
    .unwrap();
    let product = saved
        .receipt
        .products
        .iter()
        .find(|p| p.step_index == 1)
        .unwrap();
    let restored =
        workspace::load_workspace(&workspace::save_workspace(&saved.workspace).unwrap()).unwrap();
    let origin = run_material::capsule(
        &restored,
        &run_material::Selection {
            run_id: saved.receipt.run_id,
            step_index: 1,
            event: EventBinding {
                document_id: product.document_id.clone(),
                event_id: product.event_id.clone(),
                revision: product.revision.clone(),
            },
        },
    )
    .unwrap();
    let independent = Molecule::dna("Independent filled forty", CIRCLE, Topology::Linear).unwrap();
    let downstream = DesignPlan {
        name: "Reuse retained5′ extensions, fill only templated positions".into(),
        target_revision: independent.revision(),
        inputs: vec![input(
            "prepared",
            restored.revisions[&product.revision].clone(),
        )],
        steps: vec![
            Operation::ReconstructRunProduct {
                template: Source::Input("prepared".into()),
                origin,
            },
            Operation::RepairEnds {
                template: Source::Step(0),
                method: EndRepairMethod::FillFivePrimeOverhangs,
            },
        ],
    };
    let result = verification::verify(&downstream, &independent);
    assert_eq!(
        result.status,
        Status::ConsistentUnderAssumptions,
        "{:?}",
        result.findings
    );
    assert_eq!(result.products.last().unwrap().sequence.as_str(), CIRCLE);
}
#[test]
fn source_accessibility_and_single_linear_vs_self_circle_policies_are_not_implicitly_changed() {
    let (plan, source, target) = resection_plan();
    let mut unknown = plan.clone();
    unknown.steps = vec![Operation::ResectThreePrime {
        source: Source::Input("stock".into()),
        cuts: Cuts { left: 2, right: 3 },
    }];
    let report = verification::verify(&unknown, &target);
    assert_ne!(report.status, Status::ConsistentUnderAssumptions);
    assert!(report.products.is_empty());
    let mut blocked = plan;
    let Operation::DeclareSourceTermini { left, .. } = &mut blocked.steps[0] else {
        panic!()
    };
    left.three_prime = ThreePrime::Blocked;
    assert_ne!(
        verification::verify(&blocked, &target).status,
        Status::ConsistentUnderAssumptions
    );
    let mut single = DesignPlan {
        name: "No phantom single linear assembly".into(),
        target_revision: target.revision(),
        inputs: vec![input("stock", source)],
        steps: vec![Operation::HomologyAnnealing {
            fragments: vec![OrientedSource {
                source: Source::Input("stock".into()),
                reverse: false,
            }],
            options: Default::default(),
            policy: Options {
                overlaps: vec![],
                topology: Topology::Linear,
                cuts: vec![Cuts { left: 2, right: 3 }],
                min_paired_bases: 4,
                max_comparisons: 1000,
                max_candidates: 10,
            },
        }],
    };
    assert_eq!(
        verification::verify(&single, &target).status,
        Status::Invalid
    );
    let literal = format!("{}{}", CIRCLE, &CIRCLE[..8]);
    let source =
        Molecule::dna("Original self-circle precursor", &literal, Topology::Linear).unwrap();
    let target =
        Molecule::dna("Independent self-circle forty", CIRCLE, Topology::Circular).unwrap();
    single.inputs = vec![input("stock", source)];
    single.target_revision = target.revision();
    single.steps = vec![
        Operation::DeclareSourceTermini {
            template: Source::Input("stock".into()),
            left: Terminus::PCR,
            right: Terminus::PCR,
        },
        Operation::HomologyAnnealing {
            fragments: vec![OrientedSource {
                source: Source::Step(0),
                reverse: false,
            }],
            options: Default::default(),
            policy: Options {
                overlaps: vec![8],
                topology: Topology::Circular,
                cuts: vec![Cuts { left: 8, right: 8 }],
                min_paired_bases: 4,
                max_comparisons: 1000,
                max_candidates: 10,
            },
        },
    ];
    assert_eq!(
        verification::verify(&single, &target).status,
        Status::Incomplete
    );
}
#[test]
fn literal_recipe_binding_and_previous_five_prime_route_remain_distinct() {
    let (inputs, request, policy) = fixture(true);
    let target = Molecule::dna("Independent circle40", CIRCLE, Topology::Circular).unwrap();
    let mut wrong = policy.clone();
    wrong.overlaps[0] = 7;
    assert!(recipe::compile_annealing(&inputs, &target, &request, &wrong).is_err());
    let mut wrong = policy.clone();
    wrong.topology = Topology::Linear;
    wrong.overlaps.pop();
    assert!(recipe::compile_annealing(&inputs, &target, &request, &wrong).is_err());
    let old_policy = collagene_core::homology_assembly::Options {
        overlaps: request.overlaps.clone(),
        ..Default::default()
    };
    assert!(recipe::compile_homology(&inputs, &target, &request, &old_policy).is_err());
    let legacy = recipe::preview(&inputs, &target, &request).unwrap();
    assert!(matches!(
        legacy.compiled.plan.steps.last(),
        Some(Operation::ProfiledOverlap { .. })
    ));
    let (mut plan, _, target) = resection_plan();
    let before = serde_json::to_value(&plan).unwrap();
    let steps = collagene_core::plan_editing::moved(&["stock".into()], &plan.steps, 1, 1).unwrap();
    plan.steps = steps;
    assert_eq!(serde_json::to_value(&plan).unwrap(), before);
    assert_eq!(
        verification::verify(&plan, &target).status,
        Status::ConsistentUnderAssumptions
    );
}
#[test]
fn printable_protocols_describe_explicit_cuts_gaps_and_unsealed_joining_bonds() {
    let (inputs, request, mut policy) = fixture(true);
    policy.cuts[0].right = 10;
    let target = Molecule::dna("Independent circle40", CIRCLE, Topology::Circular).unwrap();
    let p = recipe::compile_annealing(&inputs, &target, &request, &policy).unwrap();
    let html = verification::export_protocol_html(&p.plan, &target);
    assert!(html.contains("exposed5′tail annealing"));
    assert!(html.contains("unfilled gaps"));
    assert!(html.contains("no host repair"));
    let pdf = collagene_core::pdf::export_protocol_pdf(&p.plan, &target).unwrap();
    assert!(pdf.bytes.starts_with(b"%PDF-"));
    assert_eq!(pdf.status, Status::Incomplete);
}
