//! Literal independently ordered 13/17-base strands; no repair or missing DNA.
use collagene_core::{
    chemistry::AnnealOptions,
    cloning::restriction::{
        catalog as catalog_digest,
        cut_selection::{Boundary, Choice},
    },
    documents::workspace::Workspace,
    material::{
        derivation::{self as native, Derivation, SourceOrigin},
        lineage::{EventBinding, InputEventBinding},
        outputs as cloning_outputs,
        run::{self as run_material, MaterialOrigin},
    },
    planning::plan::{DesignPlan, Operation, PlanInput, Source},
    verification,
};
use collagene_model::{Molecule, Topology};
fn event(w: &Workspace, index: usize) -> EventBinding {
    let d = &w.documents[index];
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
#[test]
fn unequal_native_source_conserves_both_actual_strands_in_repeated_no_cut_inspection_and_capture() {
    let top = Molecule::dna("Literal top13", "AGTCCGATGCTAC", Topology::Linear).unwrap();
    let bottom = Molecule::dna("Literal bottom17", "GTAGGTAGCATCGGACT", Topology::Linear).unwrap();
    let target = Molecule::dna(
        "Independent intended top reference",
        "AGTCCGATGCTAC",
        Topology::Linear,
    )
    .unwrap();
    let plan = DesignPlan {
        name: "Explicit unequal oligo annealing".into(),
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
    let w = Workspace::from_documents("Actual13/17 source", vec![top, bottom, target.clone()], 2)
        .unwrap();
    let report = verification::engine::verify(&plan, &target);
    assert_eq!(report.products.len(), 1, "{:?}", report.findings);
    let saved = cloning_outputs::finalize(
        &w,
        &plan,
        event(&w, 2),
        vec![
            InputEventBinding {
                input_id: "top".into(),
                event: event(&w, 0),
            },
            InputEventBinding {
                input_id: "bottom".into(),
                event: event(&w, 1),
            },
        ],
        &cloning_outputs::Options {
            allow_incomplete: true,
            ..Default::default()
        },
    )
    .unwrap();
    let original = saved.receipt.construct;
    let origin = MaterialOrigin::Run(
        collagene_core::material::capsule::capsule(
            &saved.workspace,
            &run_material::Selection {
                run_id: saved.receipt.run.run_id,
                step_index: 0,
                event: original.clone(),
            },
        )
        .unwrap(),
    );
    let c = catalog_digest::builtin_catalog().unwrap();
    let id = c
        .edition()
        .enzymes
        .iter()
        .find(|e| e.name == "EcoRI")
        .unwrap()
        .id
        .clone();
    let selection = catalog_digest::select(&c, &[id]).unwrap();
    let options = collagene_core::cloning::restriction::inspection::InspectionOptions::default();
    let source = saved.workspace.revisions[&original.revision].clone();
    let before = source.clone();
    let preview = native::digest::inspect_digest(
        "actual13/17",
        &source,
        &SourceOrigin::Saved(origin.clone()),
        &selection,
        &options,
    )
    .unwrap();
    assert!(preview.inspection.sites.is_empty());
    assert!(preview.inspection.products_determined);
    assert_eq!(
        preview.source_material.as_ref().unwrap().duplex.top_length,
        13
    );
    assert_eq!(
        preview
            .source_material
            .as_ref()
            .unwrap()
            .duplex
            .bottom_length,
        17
    );
    assert!(
        preview.inspection.validate().is_err(),
        "Legacy DTO-only validation still requires full equal-length source accounting"
    );
    preview
        .inspection
        .validate_with_source(&preview.source_material.as_ref().unwrap().duplex)
        .unwrap();
    let selected = preview
        .selected_fragment(&Choice {
            first: Boundary::SourceStart,
            second: Boundary::SourceEnd,
            complement: false,
        })
        .unwrap();
    assert_eq!(selected.duplex.top_length, 13);
    assert_eq!(selected.duplex.bottom_length, 17);
    assert_eq!(selected.duplex.bottom_start, 0);
    let mut missing_bottom = preview.inspection.clone();
    missing_bottom.fragments[0].duplex.bottom_length = 16;
    assert!(
        missing_bottom
            .validate_with_source(&preview.source_material.as_ref().unwrap().duplex)
            .is_err()
    );
    let first = native::prepare_install(
        &saved.workspace,
        &native::Request {
            sources: vec![collagene_core::material::derivation::BoundSource {
                role: collagene_core::material::derivation::SourceRole::Template,
                event: original,
                origin: SourceOrigin::Saved(origin),
            }],
            derivation: Derivation::CatalogDigest {
                conditions: None,
                selection: selection.clone(),
                options: options.clone(),
                fragment_index: 0,
            },
            name: "Retained actual13/17 fragment".into(),
            save_primers: false,
        },
    )
    .unwrap();
    let origin =
        MaterialOrigin::Derived(native::capsule(&first.workspace, &first.receipt.product).unwrap());
    let second = native::prepare_install(
        &first.workspace,
        &native::Request {
            sources: vec![collagene_core::material::derivation::BoundSource {
                role: collagene_core::material::derivation::SourceRole::Template,
                event: first.receipt.product.clone(),
                origin: SourceOrigin::Saved(origin),
            }],
            derivation: Derivation::CatalogDigest {
                conditions: None,
                selection,
                options,
                fragment_index: 0,
            },
            name: "Same actual unequal strands again".into(),
            save_primers: false,
        },
    )
    .unwrap();
    let fresh = native::reconstruct(&second.workspace, &second.receipt.product).unwrap();
    assert_eq!(fresh.tracked().duplex.top.as_str(), "AGTCCGATGCTAC");
    assert_eq!(fresh.tracked().duplex.bottom.as_str(), "GTAGGTAGCATCGGACT");
    assert_eq!(fresh.tracked().duplex.bottom_start, 0);
    assert!(fresh.tracked().duplex.nicks.is_empty());
    assert_eq!(source, before);
}
