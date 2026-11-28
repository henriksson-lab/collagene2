//! Native saved-fragment reuse through real PCR and physical overlap assembly.
use collagene_core::{
    cloning::{
        assembly::recipe::{self as a, FragmentSpec, Preparation, RunPreparation},
        overlap::assembly as homology_assembly,
        restriction::catalog as catalog_digest,
    },
    documents::workspace::{self, Workspace},
    material::{
        derivation as digest,
        lineage::{EventBinding, InputEventBinding},
        outputs as cloning_outputs, run as run_material,
    },
    planning::{
        oligos as plan_oligos,
        plan::{DesignPlan, Operation, OrientedSource, PlanInput, Source},
    },
    sequence::analysis::Primer,
    verification::{self, evidence::Status},
};
use collagene_model::{Molecule, Topology};
fn input(id: &str, m: Molecule) -> PlanInput {
    PlanInput {
        id: id.into(),
        revision: m.revision(),
        molecule: m,
    }
}
fn dna(name: &str, text: &str, topology: Topology) -> Molecule {
    Molecule::dna(name, text, topology).unwrap()
}
fn event(w: &Workspace, index: usize) -> EventBinding {
    let d = &w.documents[index];
    EventBinding {
        document_id: d.id.clone(),
        event_id: d.head_event.clone().unwrap(),
        revision: d.head.clone(),
    }
}
fn primers(f: &str, r: &str) -> RunPreparation {
    RunPreparation::ExactPcr {
        forward: Primer::new("New forward", f, "").unwrap(),
        reverse: Primer::new("New reverse", r, "").unwrap(),
        template_strand: collagene_core::chemistry::Strand::Top,
    }
}
struct Fixture {
    w: Workspace,
    fragment: EventBinding,
    capsule: digest::Capsule,
    vector: usize,
    target: usize,
    incomplete: bool,
}
fn fixture(incomplete: bool) -> Fixture {
    let (
        mut docs,
        mut plan,
        original_target,
        fragment_text,
        vector_text,
        target_text,
        fragment_index,
    ) = if incomplete {
        let left = dna("Left stock", "AGTCCGATGCTA", Topology::Linear);
        let right = dna("Right stock", "GCTATTCGAATTCCAGTC", Topology::Linear);
        let target = dna(
            "Independent closed circle",
            "AGTCCGATGCTATTCGAATTCC",
            Topology::Circular,
        );
        let policy = homology_assembly::Options {
            overlaps: vec![4, 4],
            topology: Topology::Circular,
            min_paired_bases: 4,
            max_comparisons: 3 * (12 + 18) + 2 * 22 + 8,
            ..Default::default()
        };
        let plan = DesignPlan {
            standard: None,
            name: "Original finite-search physical circle".into(),
            target_revision: target.revision(),
            inputs: vec![input("left", left.clone()), input("right", right.clone())],
            steps: vec![
                Operation::Pcr {
                    template: Source::Input("left".into()),
                    forward: Primer::new("Historical LF", "AGTC", "").unwrap(),
                    reverse: Primer::new("Historical LR", "TAGC", "").unwrap(),
                    template_strand: collagene_core::chemistry::Strand::Top,
                },
                Operation::Pcr {
                    template: Source::Input("right".into()),
                    forward: Primer::new("Historical RF", "GCTA", "").unwrap(),
                    reverse: Primer::new("Historical RR", "GACT", "").unwrap(),
                    template_strand: collagene_core::chemistry::Strand::Top,
                },
                Operation::HomologyAssembly {
                    fragments: vec![
                        OrientedSource {
                            source: Source::Step(0),
                            reverse: false,
                        },
                        OrientedSource {
                            source: Source::Step(1),
                            reverse: false,
                        },
                    ],
                    options: Default::default(),
                    policy,
                },
            ],
        };
        (
            vec![left, right, target],
            plan,
            2,
            "AATTCCAGTCCGATGCTATTCG",
            "CACGTTCGAATTCCAG",
            "CACGTTCGAATTCCAGTCCGATGCTATTCG",
            0,
        )
    } else {
        let stock = dna(
            "Original stock",
            "AGTCCGATGAATTCCGTAGGACTTCAGACCTGTAGATCGTCA",
            Topology::Linear,
        );
        let target = dna(
            "Independent original PCR",
            "AGTCCGATGAATTCCGTAGGACTTCAGACCTGTAGATCGTCA",
            Topology::Linear,
        );
        let plan = DesignPlan {
            standard: None,
            name: "Original physical PCR".into(),
            target_revision: target.revision(),
            inputs: vec![input("stock", stock.clone())],
            steps: vec![Operation::Pcr {
                template: Source::Input("stock".into()),
                forward: Primer::new("Historical forward", "AGTCCGAT", "").unwrap(),
                reverse: Primer::new("Historical reverse", "TGACGATC", "").unwrap(),
                template_strand: collagene_core::chemistry::Strand::Top,
            }],
        };
        (
            vec![stock, target],
            plan,
            1,
            "AATTCCGTAGGACTTCAGACCTGTAGATCGTCA",
            "CACGTTCGAATTCCGT",
            "CACGTTCGAATTCCGTAGGACTTCAGACCTGTAGATCGTCA",
            1,
        )
    };
    // Repeated-overlap profiling is not relevant to the chosen physical kernel.
    if let Some(Operation::HomologyAssembly { options, .. }) = plan.steps.last_mut() {
        options.scan_repeated_overlaps = false;
    }
    let vector = docs.len();
    docs.push(dna("Fresh vector stock", vector_text, Topology::Linear));
    let target = docs.len();
    docs.push(dna(
        "Independent final literal",
        target_text,
        Topology::Linear,
    ));
    let w = Workspace::from_documents("Saved fragment downstream cloning", docs, target).unwrap();
    let bindings = plan
        .inputs
        .iter()
        .enumerate()
        .map(|(i, p)| InputEventBinding {
            input_id: p.id.clone(),
            event: event(&w, i),
        })
        .collect();
    let report =
        verification::engine::verify(&plan, w.head(&w.documents[original_target].id).unwrap());
    assert_eq!(
        report.reconstructed_count(),
        plan.steps.len(),
        "{:?}",
        report.findings
    );
    assert_eq!(
        report.status,
        if incomplete {
            Status::Incomplete
        } else {
            Status::ConsistentUnderAssumptions
        },
        "{:?}",
        report.findings
    );
    let saved = cloning_outputs::finalize(
        &w,
        &plan,
        event(&w, original_target),
        bindings,
        &cloning_outputs::Options {
            allow_incomplete: incomplete,
            ..Default::default()
        },
    )
    .unwrap();
    let origin = collagene_core::material::capsule::capsule(
        &saved.workspace,
        &run_material::Selection {
            run_id: saved.receipt.run.run_id.clone(),
            step_index: plan.steps.len() - 1,
            event: saved.receipt.construct.clone(),
        },
    )
    .unwrap();
    let c = catalog_digest::builtin_catalog().unwrap();
    let id = c
        .edition()
        .enzymes
        .iter()
        .find(|e| e.name == "EcoRI")
        .unwrap()
        .id
        .clone();
    let derived = digest::prepare_install(
        &saved.workspace,
        &digest::Request {
            sources: vec![collagene_core::material::derivation::BoundSource {
                role: collagene_core::material::derivation::SourceRole::Template,
                event: saved.receipt.construct,
                origin: digest::SourceOrigin::Saved(
                    collagene_core::material::run::MaterialOrigin::Run(origin),
                ),
            }],
            derivation: digest::Derivation::CatalogDigest {
                conditions: None,
                selection: catalog_digest::select(&c, &[id]).unwrap(),
                options: Default::default(),
                fragment_index: fragment_index,
            },
            name: "Native selected digest fragment".into(),
            save_primers: false,
        },
    )
    .unwrap();
    let fragment = derived.receipt.product;
    assert_eq!(
        derived.workspace.revisions[&fragment.revision]
            .sequence
            .as_str(),
        fragment_text
    );
    let capsule = digest::capsule(&derived.workspace, &fragment).unwrap();
    Fixture {
        w: derived.workspace,
        fragment,
        capsule,
        vector,
        target,
        incomplete,
    }
}
fn compile(f: &Fixture) -> (Vec<PlanInput>, Molecule, a::Recipe, a::Compiled) {
    let fragment = f.w.revisions[&f.fragment.revision].clone();
    let vector = f.w.head(&f.w.documents[f.vector].id).unwrap().clone();
    let target = f.w.head(&f.w.documents[f.target].id).unwrap().clone();
    let reverse = if f.incomplete { "CGAATAGC" } else { "TGACGATC" };
    let vector_reverse = if f.incomplete { "CTGGAATT" } else { "ACGGAATT" };
    let forward = if f.incomplete { "AATTCCAG" } else { "AATTCCGT" };
    let inputs = vec![input("fragment", fragment), input("vector", vector)];
    let recipe = a::Recipe {
        standard: None,
        name: "Saved digest then new PCR and physical assembly".into(),
        fragments: vec![
            FragmentSpec {
                generation: None,
                role: "vector-role".into(),
                input_id: "vector".into(),
                preparation: match primers("CACGTTCG", vector_reverse) {
                    RunPreparation::ExactPcr {
                        forward, reverse, ..
                    } => Preparation::ExactPcr {
                        forward,
                        reverse,
                        template_strand: collagene_core::chemistry::Strand::Top,
                    },
                    _ => unreachable!(),
                },
                reverse: false,
            },
            FragmentSpec {
                generation: None,
                role: "fragment-role".into(),
                input_id: "fragment".into(),
                preparation: Preparation::FromDerived {
                    origin: f.capsule.clone(),
                    preparation: primers(forward, reverse),
                },
                reverse: false,
            },
        ],
        overlaps: vec![8],
        circular: false,
        options: Default::default(),
    };
    let policy = homology_assembly::Options {
        overlaps: vec![8],
        topology: Topology::Linear,
        min_paired_bases: 8,
        ..Default::default()
    };
    let compiled = a::compile_homology(&inputs, &target, &recipe, &policy).unwrap();
    (inputs, target, recipe, compiled)
}
fn finish(f: Fixture) {
    let (_, target, _, compiled) = compile(&f);
    assert_eq!(compiled.rows[1].preparation_step, Some(1));
    assert!(matches!(compiled.rows[1].prepared_source, Source::Step(2)));
    assert!(
        matches!(&compiled.plan.steps[1],Operation::ReconstructMaterial{template:Source::Input(id),..} if id=="fragment")
    );
    assert!(matches!(
        &compiled.plan.steps[2],
        Operation::Pcr {
            template: Source::Step(1),
            ..
        }
    ));
    let report = verification::engine::verify(&compiled.plan, &target);
    assert_eq!(
        report.reconstructed_count(),
        compiled.plan.steps.len(),
        "{:?}",
        report.findings
    );
    assert!(verification::engine::plan_target_matches(
        &compiled.plan,
        report
            .final_product()
            .expect("final declared step must reconstruct"),
        &target
    ));
    assert_eq!(
        report.status,
        if f.incomplete {
            Status::Incomplete
        } else {
            Status::ConsistentUnderAssumptions
        },
        "{:?}",
        report.findings
    );
    let order = plan_oligos::build(&compiled.plan, &target, &Default::default()).unwrap();
    assert_eq!(order.uses.len(), 4);
    assert!(
        order
            .uses
            .iter()
            .all(|u| !u.primer_name.starts_with("Historical"))
    );
    let saved = cloning_outputs::finalize(
        &f.w,
        &compiled.plan,
        event(&f.w, f.target),
        vec![
            InputEventBinding {
                input_id: "fragment".into(),
                event: f.fragment.clone(),
            },
            InputEventBinding {
                input_id: "vector".into(),
                event: event(&f.w, f.vector),
            },
        ],
        &cloning_outputs::Options {
            allow_incomplete: f.incomplete,
            ..Default::default()
        },
    )
    .unwrap();
    let archive = workspace::save_workspace(&saved.workspace).unwrap();
    let restored = workspace::load_workspace(&archive).unwrap();
    let fresh = run_material::reconstruct(
        &restored,
        &run_material::Selection {
            run_id: saved.receipt.run.run_id,
            step_index: compiled.assembly_step,
            event: saved.receipt.construct,
        },
    )
    .unwrap();
    assert_eq!(fresh.reference().sequence, target.sequence);
    assert_eq!(fresh.report().status == Status::Incomplete, f.incomplete);
    assert!(fresh.tracked().duplex.nicks.is_empty());
}
#[test]
fn actual_saved_pcr_digest_then_new_pcr_and_physical_assembly_records_and_replays() {
    finish(fixture(false));
}
#[test]
fn actual_saved_circle_digest_keeps_incomplete_origin_through_new_pcr_and_saved_assembly() {
    finish(fixture(true));
}
#[test]
fn changed_saved_digest_event_and_revision_cannot_be_recorded_as_matching_input() {
    let f = fixture(false);
    let (inputs, target, mut recipe, compiled) = compile(&f);
    if let Preparation::FromDerived { origin, .. } = &mut recipe.fragments[1].preparation {
        origin.product.revision = "00".repeat(32);
    }
    assert!(a::compile(&inputs, &target, &recipe).is_err());
    let mut wrong = f.fragment.clone();
    wrong.event_id = event(&f.w, 0).event_id;
    assert!(
        cloning_outputs::finalize(
            &f.w,
            &compiled.plan,
            event(&f.w, f.target),
            vec![
                InputEventBinding {
                    input_id: "fragment".into(),
                    event: wrong
                },
                InputEventBinding {
                    input_id: "vector".into(),
                    event: event(&f.w, f.vector)
                }
            ],
            &Default::default()
        )
        .is_err()
    );
}

#[test]
fn saved_digest_preparation_preserves_actual_source_modifications_before_any_new_pcr() {
    use collagene_core::chemistry::{
        methylation::MethylationState, tracked_duplex::GenerationAssumptions,
    };
    let source = dna(
        "Explicit original methylated circle",
        "AGTCCGATGAATTCCGTAGGACTTCAGACCTGTAGATCGTCA",
        Topology::Circular,
    );
    let target = dna(
        "Independent intact circle",
        "AGTCCGATGAATTCCGTAGGACTTCAGACCTGTAGATCGTCA",
        Topology::Circular,
    );
    let cut_target = dna(
        "Independent linear cut literal",
        "AATTCCGTAGGACTTCAGACCTGTAGATCGTCAAGTCCGATG",
        Topology::Linear,
    );
    let plan = DesignPlan {
        standard: None,
        name: "Original source assumptions".into(),
        target_revision: target.revision(),
        inputs: vec![input("stock", source.clone())],
        steps: vec![Operation::DeclareSourceModifications {
            template: Source::Input("stock".into()),
            assumptions: GenerationAssumptions {
                host: None,
                dam: MethylationState::FullyModified,
                dcm: MethylationState::Unknown,
                cpg: MethylationState::Unknown,
            },
        }],
    };
    let w = Workspace::from_documents(
        "Native generation ancestry",
        vec![source, target, cut_target.clone()],
        1,
    )
    .unwrap();
    let saved = cloning_outputs::finalize(
        &w,
        &plan,
        event(&w, 1),
        vec![InputEventBinding {
            input_id: "stock".into(),
            event: event(&w, 0),
        }],
        &Default::default(),
    )
    .unwrap();
    let origin = collagene_core::material::capsule::capsule(
        &saved.workspace,
        &run_material::Selection {
            run_id: saved.receipt.run.run_id,
            step_index: 0,
            event: saved.receipt.construct.clone(),
        },
    )
    .unwrap();
    let c = catalog_digest::builtin_catalog().unwrap();
    let id = c
        .edition()
        .enzymes
        .iter()
        .find(|e| e.name == "EcoRI")
        .unwrap()
        .id
        .clone();
    let derived = digest::prepare_install(
        &saved.workspace,
        &digest::Request {
            sources: vec![collagene_core::material::derivation::BoundSource {
                role: collagene_core::material::derivation::SourceRole::Template,
                event: saved.receipt.construct,
                origin: digest::SourceOrigin::Saved(
                    collagene_core::material::run::MaterialOrigin::Run(origin),
                ),
            }],
            derivation: digest::Derivation::CatalogDigest {
                conditions: None,
                selection: catalog_digest::select(&c, &[id]).unwrap(),
                options: Default::default(),
                fragment_index: 0,
            },
            name: "Native selected digest fragment".into(),
            save_primers: false,
        },
    )
    .unwrap();
    assert!(derived.receipt.modifications.known_modified > 0);
    let product = derived.receipt.product.clone();
    let prepared = a::compile_preparations(
        &[input(
            "cut",
            derived.workspace.revisions[&product.revision].clone(),
        )],
        &cut_target,
        "Reuse actual cut material",
        &[FragmentSpec {
            generation: None,
            role: "cut-role".into(),
            input_id: "cut".into(),
            reverse: false,
            preparation: Preparation::FromDerived {
                origin: digest::capsule(&derived.workspace, &product).unwrap(),
                preparation: RunPreparation::WholeLinear,
            },
        }],
        collagene_core::cloning::assembly::recipe::CleavagePolicy::AssumedComplete,
    )
    .unwrap();
    let plan = DesignPlan {
        standard: None,
        name: "Actual cut ancestry".into(),
        target_revision: prepared.target_revision,
        inputs: prepared.inputs,
        steps: prepared.steps,
    };
    let report = verification::engine::verify(&plan, &cut_target);
    assert_eq!(
        report.status,
        Status::ConsistentUnderAssumptions,
        "{:?}",
        report.findings
    );
    let evidence = &report.physical_processing_reports[0].evidence;
    let verification::evidence::PhysicalProcessingEvidence::DerivedOrigin(evidence) = evidence
    else {
        panic!("missing actual digest reconstruction evidence")
    };
    assert_eq!(
        evidence.material.modifications.known_modified,
        derived.receipt.modifications.known_modified
    );
    assert_eq!(
        evidence.material.modifications.unknown_eligible,
        derived.receipt.modifications.unknown_eligible
    );
    assert_eq!(
        serde_json::to_value(&evidence.material.modifications.generation).unwrap(),
        serde_json::to_value(&derived.receipt.modifications.generation).unwrap()
    );
    assert_eq!(
        evidence.material.duplex.top_length,
        derived.receipt.duplex.top_length
    );
}
