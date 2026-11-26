//! Independent literal strand-coverage fixtures for ordinary exact PCR.
use collagene_core::{
    Molecule, Topology,
    analysis::{self, Primer},
    chemistry::{AnnealOptions, EndTreatment},
    cloning_outputs, material_pcr,
    run_lineage::{EventBinding, InputEventBinding},
    run_material::{self, MaterialOrigin},
    verification::{DesignPlan, Operation, OrientedSource, PlanInput, Source},
    workspace::Workspace,
};
fn dna(name: &str, text: &str, topology: Topology) -> Molecule {
    Molecule::dna(name, text, topology).unwrap()
}
fn input(id: &str, molecule: Molecule) -> PlanInput {
    PlanInput {
        id: id.into(),
        revision: molecule.revision(),
        molecule,
    }
}
fn event(w: &Workspace, index: usize) -> EventBinding {
    let d = &w.documents[index];
    EventBinding {
        document_id: d.id.clone(),
        event_id: d.head_event.clone().unwrap(),
        revision: d.head.clone(),
    }
}
fn saved(docs: Vec<Molecule>, plan: DesignPlan, target: usize) -> (Molecule, MaterialOrigin) {
    let w = Workspace::from_documents("Actual strand source", docs, target).unwrap();
    let inputs = plan
        .inputs
        .iter()
        .enumerate()
        .map(|(i, p)| InputEventBinding {
            input_id: p.id.clone(),
            event: event(&w, i),
        })
        .collect();
    let result = cloning_outputs::finalize(
        &w,
        &plan,
        event(&w, target),
        inputs,
        &cloning_outputs::Options {
            allow_incomplete: true,
            ..Default::default()
        },
    )
    .unwrap();
    let source = result.workspace.revisions[&result.receipt.construct.revision].clone();
    let origin = run_material::capsule(
        &result.workspace,
        &run_material::Selection {
            run_id: result.receipt.run.run_id,
            step_index: plan.steps.len() - 1,
            event: result.receipt.construct,
        },
    )
    .unwrap();
    (source, MaterialOrigin::Run(origin))
}
fn staggered() -> (Molecule, MaterialOrigin) {
    // Actual bottom is an independently ordered 17-mer. Its plus-oriented
    // coverage is AGTCCGATGCTACCTAC, not a complement invented from top13.
    let top = dna("Top thirteen", "AGTCCGATGCTAC", Topology::Linear);
    let bottom = dna("Bottom seventeen", "GTAGGTAGCATCGGACT", Topology::Linear);
    let target = dna(
        "Independent top reference",
        "AGTCCGATGCTAC",
        Topology::Linear,
    );
    let plan = DesignPlan {
        name: "Two literal unequal annealed oligos".into(),
        target_revision: target.revision(),
        inputs: vec![input("top", top.clone()), input("bottom", bottom.clone())],
        steps: vec![Operation::AnnealOligos {
            top: Source::Input("top".into()),
            bottom: Source::Input("bottom".into()),
            options: AnnealOptions {
                min_paired_bases: 12,
                ..Default::default()
            },
            candidate_index: 0,
        }],
    };
    saved(vec![top, bottom, target], plan, 2)
}
#[test]
fn actual_thirteen_seventeen_strands_support_top_contained_pcr_without_silent_repair() {
    let (source, origin) = staggered();
    let before = source.clone();
    let forward = Primer::new("F", "AGTC", "").unwrap();
    let reverse = Primer::new("R", "GTAGC", "").unwrap();
    let previews = material_pcr::simulate(&source, &origin, &forward, &reverse).unwrap();
    assert_eq!(previews.len(), 1);
    assert_eq!(previews[0].molecule.sequence.as_str(), "AGTCCGATGCTAC");
    assert!(
        previews[0]
            .findings
            .iter()
            .any(|s| s.contains("initially staggered"))
    );
    assert_eq!(source, before);
}
#[test]
fn bottom_only_extra_bases_creating_second_geometry_are_explicit_unsupported_not_false_unique() {
    let (source, origin) = staggered();
    let before = source.clone();
    let forward = Primer::new("F", "AGTC", "").unwrap();
    let reverse = Primer::new("R", "GTAG", "").unwrap();
    assert_eq!(
        analysis::simulate_pcr_geometry(&source, &forward, &reverse)
            .unwrap()
            .len(),
        1
    );
    let union = dna(
        "Independently known strand coverage",
        "AGTCCGATGCTACCTAC",
        Topology::Linear,
    );
    assert_eq!(
        analysis::simulate_pcr_geometry(&union, &forward, &reverse)
            .unwrap()
            .len(),
        2
    );
    let error = material_pcr::simulate(&source, &origin, &forward, &reverse)
        .unwrap_err()
        .to_string();
    assert!(
        error.contains("outside the top-reference-contained"),
        "{error}"
    );
    assert_eq!(source, before);
}
#[test]
fn freshly_replayed_unsealed_nicks_remain_unsupported_for_pcr_and_source_is_unchanged() {
    let left = dna("Left original", "TTGCAATCGTCGATGCTAGC", Topology::Linear);
    let right = dna("Right original", "ATCGTCGATGCTAGCGGTAC", Topology::Linear);
    let target = dna(
        "Independent nicked circle",
        "TTGCAATCGTCGATGCTAGCATCGTCGATGCTAGCGGTAC",
        Topology::Circular,
    );
    let plan = DesignPlan {
        name: "One phosphorylated and one hydroxyl fragment".into(),
        target_revision: target.revision(),
        inputs: vec![input("left", left.clone()), input("right", right.clone())],
        steps: vec![
            Operation::Pcr {
                template: Source::Input("left".into()),
                forward: Primer::new("LF", "TTGCAATC", "").unwrap(),
                reverse: Primer::new("LR", "GCTAGCAT", "").unwrap(),
            },
            Operation::TreatEnds {
                template: Source::Step(0),
                treatment: EndTreatment::Phosphorylate,
            },
            Operation::Pcr {
                template: Source::Input("right".into()),
                forward: Primer::new("RF", "ATCGTCGA", "").unwrap(),
                reverse: Primer::new("RR", "GTACCGCT", "").unwrap(),
            },
            Operation::Ligate {
                fragments: vec![
                    OrientedSource {
                        source: Source::Step(1),
                        reverse: false,
                    },
                    OrientedSource {
                        source: Source::Step(2),
                        reverse: false,
                    },
                ],
                circular: true,
                allow_nicked: true,
            },
        ],
    };
    let (source, origin) = saved(vec![left, right, target], plan, 2);
    let before = source.clone();
    let forward = Primer::new("F", "TTGCAATC", "").unwrap();
    let reverse = Primer::new("R", "GTACCGCT", "").unwrap();
    let error = material_pcr::simulate(&source, &origin, &forward, &reverse)
        .unwrap_err()
        .to_string();
    assert!(error.contains("nick propagation"), "{error}");
    assert_eq!(source, before);
}
