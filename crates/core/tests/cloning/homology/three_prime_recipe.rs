//! Independent ordinary-cloning declarations; joined diagnostic DNA never becomes material.
use collagene_core::{
    chemistry::{self, EndRepairMethod, FivePrime, Terminus, ThreePrime},
    cloning::{
        assembly::recipe::{self, FragmentSpec, Preparation},
        overlap::{
            annealing::{Options, Report},
            resection::Cuts,
        },
    },
    documents::workspace::{self, Workspace},
    material::{
        lineage::{self as run_lineage, EventBinding, InputEventBinding},
        run as run_material,
    },
    planning::plan::{DesignPlan, Operation, OrientedSource, PlanInput, Source},
    sequence::analysis::Primer,
    verification::{
        self,
        evidence::{PhysicalProcessingEvidence, Status},
    },
};
use collagene_model::{Feature, Location, Molecule, Topology};
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
    let mut profile = collagene_core::cloning::overlap::profiles::ProfileOptions::default();
    profile.scan_repeated_overlaps = false;
    profile.method = collagene_core::cloning::overlap::profiles::OverlapMethod::InFusion;
    let request = recipe::Recipe {
        name: "Declared physical homology route".into(),
        fragments,
        overlaps: overlaps.clone(),
        circular,
        options: profile,
    };
    let policy = Options {
        planned_rescue: Some(
            collagene_core::cloning::overlap::annealing::PlannedRescue::EcoliTransformation,
        ),
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
        p.verification
            .final_product()
            .expect("final declared step must reconstruct")
            .sequence
            .as_str(),
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
    let product = p
        .verification
        .final_product()
        .expect("final declared step must reconstruct");
    assert!(chemistry::has_unmapped_physical_evidence(product));
    assert!(chemistry::Duplex::from_molecule(product).is_err());
    let junctions = collagene_core::cloning::assembly::junctions::inspect_plan(
        &p.compiled.plan,
        &target,
        &Default::default(),
    )
    .unwrap();
    assert!(matches!(
        junctions.method,
        Some(collagene_core::cloning::assembly::junctions::Method::HomologyAnnealing)
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
        a.verification
            .final_product()
            .expect("final declared step must reconstruct")
            .sequence,
        b.verification
            .final_product()
            .expect("final declared step must reconstruct")
            .sequence
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
        p.verification
            .final_product()
            .expect("final declared step must reconstruct")
            .sequence
            .as_str(),
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
        let result = verification::engine::verify(&plan, &target);
        assert_ne!(result.status, Status::ConsistentUnderAssumptions);
        assert_eq!(result.reconstructed_count(), 3);
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
    assert!(collagene_core::material::capsule::capsule(&loaded, &selection).is_err());
    let evidence =
        collagene_core::material::physical_evidence::load(&loaded, &saved.receipt.run_id)
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
    let report = verification::engine::verify(&plan, &target);
    assert_eq!(
        report.status,
        Status::ConsistentUnderAssumptions,
        "{:?}",
        report.findings
    );
    let actual = report.product(1).expect("selected step must reconstruct");
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
    let origin = collagene_core::material::capsule::capsule(
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
    let result = verification::engine::verify(&downstream, &independent);
    assert_eq!(
        result.status,
        Status::ConsistentUnderAssumptions,
        "{:?}",
        result.findings
    );
    assert_eq!(
        result
            .final_product()
            .expect("final declared step must reconstruct")
            .sequence
            .as_str(),
        CIRCLE
    );
}
#[test]
fn source_accessibility_and_single_linear_vs_self_circle_policies_are_not_implicitly_changed() {
    let (plan, source, target) = resection_plan();
    let mut unknown = plan.clone();
    unknown.steps = vec![Operation::ResectThreePrime {
        source: Source::Input("stock".into()),
        cuts: Cuts { left: 2, right: 3 },
    }];
    let report = verification::engine::verify(&unknown, &target);
    assert_ne!(report.status, Status::ConsistentUnderAssumptions);
    assert!((report.reconstructed_count() == 0));
    let mut blocked = plan;
    let Operation::DeclareSourceTermini { left, .. } = &mut blocked.steps[0] else {
        panic!()
    };
    left.three_prime = ThreePrime::Blocked;
    assert_ne!(
        verification::engine::verify(&blocked, &target).status,
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
                planned_rescue: None,
                overlaps: vec![],
                topology: Topology::Linear,
                cuts: vec![Cuts { left: 2, right: 3 }],
                min_paired_bases: 4,
                max_comparisons: 1000,
                max_candidates: 10,
            },
        }],
    };
    let refused = verification::engine::verify(&single, &target);
    assert_eq!(refused.status, Status::Incomplete);
    assert!(refused.product(0).is_none());
    assert_eq!(
        refused.steps[0].state,
        verification::evidence::StepState::Withheld
    );
    assert!(refused.findings.iter().any(|finding| finding.rule
        == verification::checks::RuleId::HomologyPolicy
        && finding.verdict == verification::checks::CheckVerdict::Indeterminate));
    refused.validate_documentary().unwrap();
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
                planned_rescue: None,
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
        verification::engine::verify(&single, &target).status,
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
    let old_policy = collagene_core::cloning::overlap::assembly::Options {
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
    let steps =
        collagene_core::planning::editing::moved(&["stock".into()], &plan.steps, 1, 1).unwrap();
    plan.steps = steps;
    assert_eq!(serde_json::to_value(&plan).unwrap(), before);
    assert_eq!(
        verification::engine::verify(&plan, &target).status,
        Status::ConsistentUnderAssumptions
    );
}
#[test]
fn printable_protocols_describe_explicit_cuts_gaps_and_unsealed_joining_bonds() {
    let (inputs, request, mut policy) = fixture(true);
    policy.cuts[0].right = 10;
    let target = Molecule::dna("Independent circle40", CIRCLE, Topology::Circular).unwrap();
    let p = recipe::compile_annealing(&inputs, &target, &request, &policy).unwrap();
    let html = verification::output::export_protocol_html(&p.plan, &target);
    assert!(html.contains("exposed5′tail annealing"));
    assert!(html.contains("unfilled gaps"));
    assert!(html.contains("no host repair"));
    let pdf = collagene_core::output::pdf::export_protocol_pdf(&p.plan, &target).unwrap();
    assert!(pdf.bytes.starts_with(b"%PDF-"));
    assert_eq!(pdf.status, Status::Incomplete);
}

#[test]
fn infusion_rescue_is_explicit_documentary_and_cannot_change_physical_completion() {
    let (inputs, request, mut policy) = fixture(true);
    let target = Molecule::dna("Independent circle", CIRCLE, Topology::Circular).unwrap();
    policy.planned_rescue = None;
    assert!(recipe::compile_annealing(&inputs, &target, &request, &policy).is_err());
    policy.planned_rescue =
        Some(collagene_core::cloning::overlap::annealing::PlannedRescue::EcoliTransformation);
    let preview = recipe::preview_annealing(&inputs, &target, &request, &policy).unwrap();
    assert_eq!(preview.verification.status, Status::Incomplete);
    assert_eq!(
        annealed(&preview).options.planned_rescue,
        policy.planned_rescue
    );
    let mut generic = request.clone();
    generic.options.method = collagene_core::cloning::overlap::profiles::OverlapMethod::HiFi;
    policy.planned_rescue = None;
    let no_rescue = recipe::preview_annealing(&inputs, &target, &generic, &policy).unwrap();
    assert_eq!(annealed(&preview).reference, annealed(&no_rescue).reference);
    assert_eq!(annealed(&preview).segments, annealed(&no_rescue).segments);
    let html = verification::output::export_protocol_html(&preview.compiled.plan, &target);
    assert!(html.contains("Planned next step"));
    assert!(html.contains("not simulated"));
}

#[test]
fn infusion_native_ecori_vector_reverse_insert_and_linker_keep_bottom_extension_calls() {
    use collagene_core::cloning::overlap::{
        annealing::PlannedRescue,
        profiles::{OverlapMethod, ProfileOptions},
    };
    let mut vector = Molecule::dna(
        "EcoRI vector",
        "AATTCCTACGATTCGACGTACCGGATACCGTACCTGAGCTAACG",
        Topology::Circular,
    )
    .unwrap();
    vector.features.push(Feature {
        kind: "misc_feature".into(),
        location: Location::range(0, 4),
        qualifiers: vec![("label".into(), Some("Existing vector calls".into()))],
        needs_review: false,
    });
    let mut first =
        Molecule::dna("Reverse insert", "TTCAGACCTGTAGATCGTCA", Topology::Linear).unwrap();
    first.features.push(Feature {
        kind: "misc_feature".into(),
        location: Location::range(0, 8),
        qualifiers: vec![("label".into(), Some("Reversed insert calls".into()))],
        needs_review: false,
    });
    let second = Molecule::dna(
        "Second insert",
        "TGCAGATCCGTAACGCTTGAAGTCCGATGCTACGTAGGAC",
        Topology::Linear,
    )
    .unwrap();
    let inputs = vec![
        input("vector", vector),
        input("first", first),
        input("second", second),
    ];
    // The independent intended genotype retains the four real bottom-extension
    // AATT calls and the explicitly declared GGC linker; neither is target repair.
    let target = Molecule::dna("Independent111", "AATTCCTACGATTCGACGTACCGGATACCGTACCTGAGCTAACGAATTTGACGATCTACAGGTCTGAAGGCTGCAGATCCGTAACGCTTGAAGTCCGATGCTACGTAGGAC", Topology::Circular).unwrap();
    assert_eq!(target.sequence.len(), 111);
    let request = recipe::Recipe {
        name: "Three-part native-end In-Fusion declaration".into(),
        circular: true,
        overlaps: vec![20, 20, 20],
        options: ProfileOptions {
            method: OverlapMethod::InFusion,
            scan_repeated_overlaps: false,
            ..Default::default()
        },
        fragments: vec![
            FragmentSpec {
                role: "vector".into(),
                input_id: "vector".into(),
                generation: None,
                reverse: false,
                preparation: Preparation::Digest {
                    enzymes: vec!["EcoRI".into()],
                    fragment_index: 0,
                    repair: None,
                },
            },
            FragmentSpec {
                role: "first".into(),
                input_id: "first".into(),
                generation: None,
                reverse: true,
                preparation: Preparation::ExactPcr {
                    forward: Primer::new(
                        "First F + linker/overlap",
                        "TTCAGACC",
                        "TCAAGCGTTACGGATCTGCAGCC",
                    )
                    .unwrap(),
                    reverse: Primer::new(
                        "First R + vector overlap",
                        "TGACGATC",
                        "CGTACCTGAGCTAACGAATT",
                    )
                    .unwrap(),
                    template_strand: chemistry::Strand::Top,
                },
            },
            FragmentSpec {
                role: "second".into(),
                input_id: "second".into(),
                generation: None,
                reverse: false,
                preparation: Preparation::ExactPcr {
                    forward: Primer::new("Second F", "TGCAGATC", "").unwrap(),
                    reverse: Primer::new(
                        "Second R + vector overlap",
                        "GTCCTACG",
                        "TACGTCGAATCGTAGGAATT",
                    )
                    .unwrap(),
                    template_strand: chemistry::Strand::Top,
                },
            },
        ],
    };
    let policy = Options {
        overlaps: vec![20, 20, 20],
        topology: Topology::Circular,
        cuts: vec![
            Cuts {
                left: 16,
                right: 16,
            },
            Cuts {
                left: 20,
                right: 20,
            },
            Cuts {
                left: 20,
                right: 20,
            },
        ],
        min_paired_bases: 12,
        max_comparisons: 100_000,
        max_candidates: 50,
        planned_rescue: Some(PlannedRescue::EcoliTransformation),
    };
    let preview = recipe::preview_annealing(&inputs, &target, &request, &policy).unwrap();
    assert_eq!(
        preview.verification.status,
        Status::Incomplete,
        "{:?}",
        preview.verification
    );
    let product = preview
        .verification
        .final_product()
        .expect("final declared step must reconstruct");
    let reversed_feature = product
        .features
        .iter()
        .find(|feature| {
            feature.qualifiers.iter().any(|(key, value)| {
                key == "label" && value.as_deref() == Some("Reversed insert calls")
            })
        })
        .unwrap();
    assert_eq!(
        reversed_feature.location,
        Location::Complement(Box::new(Location::range(60, 68)))
    );
    let evidence = annealed(&preview);
    assert_eq!(evidence.reference.as_str(), target.sequence.as_str());
    assert_eq!(evidence.segments.len(), 6);
    assert_eq!(
        evidence
            .segments
            .iter()
            .map(|s| s.sequence.len())
            .sum::<usize>(),
        222
    );
    assert_eq!(evidence.junctions.len(), 3);
    assert!(
        evidence
            .junctions
            .iter()
            .all(|j| j.top.start == j.top.end && j.bottom.start == j.bottom.end)
    );
    assert_eq!(evidence.preparations[0].product.top.len(), 28);
    assert_eq!(evidence.preparations[0].product.bottom.len(), 28);
    let html = verification::output::export_protocol_html(&preview.compiled.plan, &target);
    assert!(html.contains("Planned next step"));
    assert!(
        collagene_core::output::pdf::export_protocol_pdf(&preview.compiled.plan, &target)
            .unwrap()
            .bytes
            .starts_with(b"%PDF-")
    );
    let mut wrong = target.clone();
    wrong.sequence =
        collagene_model::Sequence::dna(&format!("{}A", &target.sequence.as_str()[..110])).unwrap();
    let mismatch = recipe::preview_annealing(&inputs, &wrong, &request, &policy).unwrap();
    assert_eq!(mismatch.verification.status, Status::Invalid);
    assert_eq!(annealed(&mismatch).reference, evidence.reference);
}
