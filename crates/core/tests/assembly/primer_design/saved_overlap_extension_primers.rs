//! Literal actual-strand OE primer proposals from exact native Run/Derived events.
use collagene_core::{
    assembly_edit_intent::Edit,
    chemistry::{AnnealOptions, EndRepairMethod, Strand},
    cloning_outputs, material_derivation,
    overlap_extension_primer_design::{self as oe, BoundOrigin, Request},
    overlap_primer_design::Fragment,
    run_lineage::{EventBinding, InputEventBinding},
    run_material::{self, MaterialOrigin},
    verification::{DesignPlan, Operation, PlanInput, Source, Status},
    workspace::{self, Workspace},
};
use collagene_model::{Molecule, Topology};
const TOP: &str = "AGTCCGATGCTACGTAGGACTTCAGACCTGTAGATCGTCA";
const BOTTOM: &str = "GGGGTGACGATCTACAGGTCTGAAGTCCTACGTAGCATCGGACT";
const EDITED_BOTTOM: &str = "GGGGTGACGATCTACAAATGAAGTCCTACGTAGCATCGGACT";
fn pin(w: &Workspace, i: usize) -> EventBinding {
    let d = &w.documents[i];
    EventBinding {
        document_id: d.id.clone(),
        event_id: d.head_event.clone().unwrap(),
        revision: d.head.clone(),
    }
}
fn input(id: &str, molecule: Molecule) -> PlanInput {
    PlanInput {
        id: id.into(),
        revision: molecule.revision(),
        molecule,
    }
}
fn saved() -> (Workspace, EventBinding, MaterialOrigin) {
    let top = Molecule::dna("Actual top40", TOP, Topology::Linear).unwrap();
    let bottom = Molecule::dna("Actual bottom44", BOTTOM, Topology::Linear).unwrap();
    let target = Molecule::dna("Independent annealing top", TOP, Topology::Linear).unwrap();
    let w = Workspace::from_documents(
        "Unequal actual native source",
        vec![top.clone(), bottom.clone(), target.clone()],
        2,
    )
    .unwrap();
    let plan = DesignPlan {
        name: "Explicit40/44 annealing".into(),
        target_revision: target.revision(),
        inputs: vec![input("top", top), input("bottom", bottom)],
        steps: vec![Operation::AnnealOligos {
            top: Source::Input("top".into()),
            bottom: Source::Input("bottom".into()),
            options: AnnealOptions {
                // Stop after the real40-call candidate before the last permitted offset.
                // This deliberately retains search incompleteness through reuse.
                min_paired_bases: 39,
                max_candidates: 1,
                ..Default::default()
            },
            candidate_index: 0,
        }],
    };
    let installed = cloning_outputs::finalize(
        &w,
        &plan,
        pin(&w, 2),
        vec![
            InputEventBinding {
                input_id: "top".into(),
                event: pin(&w, 0),
            },
            InputEventBinding {
                input_id: "bottom".into(),
                event: pin(&w, 1),
            },
        ],
        &cloning_outputs::Options {
            name: "Saved actual40/44".into(),
            allow_incomplete: true,
            ..Default::default()
        },
    )
    .unwrap();
    let event = installed.receipt.construct;
    let origin = MaterialOrigin::Run(
        run_material::capsule(
            &installed.workspace,
            &run_material::Selection {
                run_id: installed.receipt.run.run_id,
                step_index: 0,
                event: event.clone(),
            },
        )
        .unwrap(),
    );
    (installed.workspace, event, origin)
}
fn request(event: EventBinding, origin: MaterialOrigin) -> Request {
    let mut o = oe::Options::default();
    o.homology_flank = 4;
    o.min_paired_bases = 8;
    o.primer_options.min_length = 8;
    o.primer_options.max_length = 8;
    o.primer_options.max_pairs = 1;
    o.primer_options.min_gc_fraction = 0.;
    o.primer_options.max_gc_fraction = 1.;
    o.primer_options.max_tm_difference_c = 100.;
    o.primer_options.max_homopolymer = 200;
    o.primer_options.max_three_prime_complementarity = 200;
    Request {
        fragment: Fragment {
            role: "insert".into(),
            input_id: "source".into(),
            start: 0,
            end: 44,
            reverse: false,
        },
        edits: vec![Edit {
            role: "insert".into(),
            start: 16,
            end: 20,
            expected: "GGTC".into(),
            replacement: "AA".into(),
        }],
        options: o,
        template_strand: Strand::Bottom,
        material_binding: Some(BoundOrigin { event, origin }),
    }
}
fn check(w: &Workspace, event: EventBinding, origin: MaterialOrigin) {
    let source = &w.revisions[&event.revision];
    let before = workspace::save_workspace(w).unwrap();
    let r = request(event, origin);
    let context = oe::design_revision(source, &r, Some(w)).unwrap();
    let report = oe::design(source, &r, Some(w)).unwrap();
    report.validate_sources(source, &context, Some(w)).unwrap();
    assert_eq!(report.source_frame.sequence.as_str(), BOTTOM);
    assert_eq!(report.source_frame.template_strand, Strand::Bottom);
    assert_eq!(report.expected.sequence.as_str(), EDITED_BOTTOM);
    assert_eq!(
        report
            .verification
            .products
            .last()
            .unwrap()
            .sequence
            .as_str(),
        EDITED_BOTTOM
    );
    assert_eq!(
        report.verification.status,
        Status::Incomplete,
        "{:?}",
        report.verification.findings
    );
    assert_eq!(report.plan.inputs[0].molecule, *source);
    assert_eq!(report.plan.inputs.len(), 1);
    assert_eq!(report.plan.steps.len(), 5);
    assert_eq!(report.preparation.left.template_strand, Strand::Bottom);
    assert_eq!(report.preparation.right.template_strand, Strand::Bottom);
    assert_eq!(report.preparation.outer.template_strand, Strand::Top);
    assert!(
        report.search.iter().all(|s| s.unique_geometry),
        "{:?}",
        report.search
    );
    assert_eq!(workspace::save_workspace(w).unwrap(), before);
    assert!(oe::design(source, &r, None).is_err());
    let mut unbound = r.clone();
    unbound.material_binding = None;
    assert!(oe::design(source, &unbound, Some(w)).is_err());
    let mut other = r.clone();
    other.template_strand = Strand::Top;
    assert!(oe::design(source, &other, Some(w)).is_err()); // 44 calls are absent from actual top40
    let mut wrong = r.clone();
    wrong.material_binding.as_mut().unwrap().event = pin(w, 0);
    assert!(oe::preflight(source, &wrong, Some(w)).is_err());
    assert!(
        report
            .validate_sources(
                source,
                &oe::design_revision(
                    source,
                    &Request {
                        fragment: Fragment {
                            reverse: true,
                            ..r.fragment.clone()
                        },
                        ..r.clone()
                    },
                    Some(w)
                )
                .unwrap(),
                Some(w)
            )
            .is_err()
    );
}
#[test]
fn exact_run_source_uses_actual_bottom_extension_and_retains_unsealed_origin_findings() {
    let (w, event, origin) = saved();
    let w = workspace::load_workspace(&workspace::save_workspace(&w).unwrap()).unwrap();
    check(&w, event, origin);
}
#[test]
fn derived_source_remains_a_native_recipe_and_never_becomes_an_original_bottom_stock() {
    let (w, event, origin) = saved();
    let next = material_derivation::prepare_install(
        &w,
        &material_derivation::Request {
            sources: vec![collagene_core::material_derivation::BoundSource {
                role: collagene_core::material_derivation::SourceRole::Template,
                event: event,
                origin: material_derivation::SourceOrigin::Saved(origin),
            }],
            derivation: material_derivation::Derivation::EndRepair {
                method: EndRepairMethod::FillFivePrimeOverhangs,
            },
            name: "Explicit filled native derivative".into(),
            save_primers: false,
        },
    )
    .unwrap();
    let event = next.receipt.product;
    let origin =
        MaterialOrigin::Derived(material_derivation::capsule(&next.workspace, &event).unwrap());
    let source = &next.workspace.revisions[&event.revision];
    let r = request(event, origin);
    let p = oe::design(source, &r, Some(&next.workspace)).unwrap();
    assert!(matches!(
        p.plan.steps[0],
        Operation::ReconstructMaterial { .. }
    ));
    assert_eq!(p.expected.sequence.as_str(), EDITED_BOTTOM);
    assert_eq!(
        p.verification.products.last().unwrap().sequence.as_str(),
        EDITED_BOTTOM
    );
    assert!(
        p.verification
            .findings
            .iter()
            .any(|s| s.contains("incomplete")
                || s.contains("unsealed")
                || s.contains("conditional")),
        "{:?}",
        p.verification.findings
    );
    p.validate_sources(source, &p.context_revision, Some(&next.workspace))
        .unwrap();
}
