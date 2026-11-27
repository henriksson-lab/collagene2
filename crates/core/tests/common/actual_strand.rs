//! Shared literal source fixture, not a production origin resolver.
use collagene_core::{
    Molecule, Topology,
    assembly_intent::{self, Contribution, MaterialBinding, Selector},
    chemistry::{AnnealOptions, Strand, Terminus},
    cloning_outputs,
    run_lineage::{EventBinding, InputEventBinding},
    run_material::{self, MaterialSelection},
    verification::{DesignPlan, Operation, PlanInput, Source},
    workspace::Workspace,
};
pub fn event(w: &Workspace, i: usize) -> EventBinding {
    let d = &w.documents[i];
    EventBinding {
        document_id: d.id.clone(),
        event_id: d.head_event.clone().unwrap(),
        revision: d.head.clone(),
    }
}
pub fn input(id: &str, m: Molecule) -> PlanInput {
    PlanInput {
        id: id.into(),
        revision: m.revision(),
        molecule: m,
    }
}
pub fn fixture() -> (
    Workspace,
    MaterialSelection,
    PlanInput,
    assembly_intent::Request,
) {
    let top = Molecule::dna("Top21", "AATTGCTACGGAATTCTTACG", Topology::Linear).unwrap();
    let bottom = Molecule::dna("Bottom21", "AATTCGTAAGAATTCCGTAGC", Topology::Linear).unwrap();
    let target = Molecule::dna(
        "Independent top21",
        "AATTGCTACGGAATTCTTACG",
        Topology::Linear,
    )
    .unwrap();
    let w = Workspace::from_documents(
        "Literal native strand origin",
        vec![top.clone(), bottom.clone(), target.clone()],
        2,
    )
    .unwrap();
    let plan = DesignPlan {
        name: "Actual staggered annealing".into(),
        target_revision: target.revision(),
        inputs: vec![input("top", top), input("bottom", bottom)],
        steps: vec![Operation::AnnealOligos {
            top: Source::Input("top".into()),
            bottom: Source::Input("bottom".into()),
            options: AnnealOptions {
                min_paired_bases: 17,
                top_chemistry: Terminus::CLEAVED,
                bottom_chemistry: Terminus::CLEAVED,
                ..Default::default()
            },
            candidate_index: 0,
        }],
    };
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
            name: "Saved actual strands".into(),
            ..Default::default()
        },
    )
    .unwrap();
    let selection = MaterialSelection::Run(run_material::Selection {
        run_id: saved.receipt.run.run_id,
        step_index: 0,
        event: saved.receipt.construct.clone(),
    });
    let origins =
        run_material::documentary_material_origins(&saved.workspace, &[selection.clone()]).unwrap();
    let reference = saved.workspace.revisions[&selection.event().revision].clone();
    let fresh = origins[0].replay(&reference).unwrap();
    assert_eq!(fresh.tracked().duplex.bottom_start, 4);
    assert_eq!(
        fresh.tracked().duplex.bottom.as_str(),
        "AATTCGTAAGAATTCCGTAGC"
    );
    let state = collagene_core::modified_duplex::summary(fresh.tracked())
        .unwrap()
        .state_revision;
    let request = assembly_intent::Request {
        name: "Independent actual-bottom window".into(),
        fragments: vec![Contribution {
            role: "reverse digest requirement".into(),
            input_id: "saved".into(),
            selector: Selector::MaterialStrandRange {
                strand: Strand::Bottom,
                start: 0,
                end: 10,
                expected_state_revision: state,
            },
        }],
        circular: true,
        material_bindings: vec![MaterialBinding {
            input_id: "saved".into(),
            event: selection.event().clone(),
            origin: origins.into_iter().next().unwrap(),
        }],
    };
    (
        saved.workspace,
        selection,
        input("saved", reference),
        request,
    )
}
