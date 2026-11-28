use collagene_core::{
    assembly_intent::Selector,
    assembly_recipe::{Preparation, RunPreparation},
    catalog_digest,
    chemistry::{AnnealOptions, Strand, Terminus},
    cloning_outputs, combinatorial_batch as batch, combinatorial_design as library,
    combinatorial_release as release, modified_duplex,
    run_lineage::{EventBinding, InputEventBinding},
    run_material::{self, MaterialOrigin},
    verification::{DesignPlan, Operation, PlanInput, Source, Status},
    workspace::{self, Workspace},
};
use collagene_model::{Feature, Location, Molecule, Topology};
pub const EXPECTED: &str = "AATTCGTAAG";
pub fn event(w: &Workspace, index: usize) -> EventBinding {
    let d = &w.documents[index];
    EventBinding {
        document_id: d.id.clone(),
        event_id: d.head_event.clone().unwrap(),
        revision: d.head.clone(),
    }
}
pub fn native_source(
    top_text: &str,
    bottom_text: &str,
    min_pair: usize,
) -> (Workspace, EventBinding, MaterialOrigin) {
    let mut top = Molecule::dna("Actual top", top_text, Topology::Linear).unwrap();
    if top_text == "AATTGCTACGGAATTCTTACG" {
        top.features.push(Feature {
            kind: "misc_feature".into(),
            location: Location::range(15, 21),
            qualifiers: vec![("label".into(), Some("Exact top mapping".into()))],
            needs_review: false,
        });
    }
    let bottom = Molecule::dna("Actual bottom", bottom_text, Topology::Linear).unwrap();
    let target = Molecule::dna("Independent source top", top_text, Topology::Linear).unwrap();
    let private = Molecule::dna(
        "Unrelated private head",
        "CCCCAAAAGGGGTTTT",
        Topology::Circular,
    )
    .unwrap();
    let w = Workspace::from_documents(
        "Actual library sources",
        vec![top.clone(), bottom.clone(), target.clone(), private],
        0,
    )
    .unwrap();
    let plan = DesignPlan {
        name: "Actual unequal-offset native strands".into(),
        target_revision: target.revision(),
        inputs: vec![
            PlanInput {
                id: "top".into(),
                revision: top.revision(),
                molecule: top,
            },
            PlanInput {
                id: "bottom".into(),
                revision: bottom.revision(),
                molecule: bottom,
            },
        ],
        steps: vec![Operation::AnnealOligos {
            top: Source::Input("top".into()),
            bottom: Source::Input("bottom".into()),
            options: AnnealOptions {
                min_paired_bases: min_pair,
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
            name: "Saved actual source".into(),
            show_prepared: false,
            save_primers: false,
            allow_incomplete: false,
        },
    )
    .unwrap();
    let source = saved.receipt.construct;
    let w =
        workspace::load_workspace(&workspace::save_workspace(&saved.workspace).unwrap()).unwrap();
    let capsule = run_material::capsule(
        &w,
        &run_material::Selection {
            run_id: saved.receipt.run.run_id,
            step_index: 0,
            event: source.clone(),
        },
    )
    .unwrap();
    (w, source, MaterialOrigin::Run(capsule))
}
pub fn fixture() -> (Workspace, library::Design, library::SelectedConstruct) {
    let (w, source, origin) = native_source("AATTGCTACGGAATTCTTACG", "AATTCGTAAGAATTCCGTAGC", 17);
    let MaterialOrigin::Run(capsule) = origin.clone() else {
        unreachable!()
    };
    let fresh = origin.replay(&w.revisions[&source.revision]).unwrap();
    assert_eq!(fresh.tracked().duplex.bottom_start, 4);
    let state = modified_duplex::summary(fresh.tracked())
        .unwrap()
        .state_revision;
    let catalog = catalog_digest::builtin_catalog().unwrap();
    let eco = catalog
        .edition()
        .enzymes
        .iter()
        .find(|e| e.name == "EcoRI")
        .unwrap();
    let selection = catalog_digest::select(&catalog, &[eco.id.clone()]).unwrap();
    let design = library::Design {
        id: "strand-library".into(),
        name: "Actual bottom library".into(),
        circular: true,
        method: library::Method::Ligation {
            allow_nicked: false,
        },
        material_bindings: vec![library::LibraryMaterialBinding {
            event: source.clone(),
            origin,
        }],
        bins: vec![library::Bin {
            id: "insert".into(),
            label: "Sticky native insert".into(),
            choices: vec![library::Choice {
                id: "actual-bottom".into(),
                label: "Actual stored bottom0..10".into(),
                source,
                selector: Selector::MaterialStrandRange {
                    strand: Strand::Bottom,
                    start: 0,
                    end: 10,
                    expected_state_revision: state,
                },
                preparation_reverse: true,
                preparation: Some(Preparation::FromRun {
                    origin: capsule,
                    preparation: RunPreparation::CatalogDigest {
                        selection,
                        fragment_index: 1,
                        repair: None,
                    },
                }),
                treatment: None,
                generation: None,
            }],
        }],
    };
    let row = library::SelectedConstruct {
        id: "native-circle".into(),
        name: "Independent bottom10 circle".into(),
        choices: vec![library::ChoiceSelection {
            bin_id: "insert".into(),
            choice_id: "actual-bottom".into(),
        }],
        target: library::TargetDeclaration::OriginalRegions {
            junction_edits: vec![],
            substitutions: vec![],
            edits: vec![],
        },
        overrides: None,
    };
    (w, design, row)
}
pub fn captured() -> (batch::Commit, String) {
    let (w, design, row) = fixture();
    let mut session = batch::Session::begin(
        w,
        batch::Request {
            design: Some(design),
            rows: vec![batch::RowRequest::Ready(row)],
        },
    )
    .unwrap();
    let context = session.summary().context_revision;
    let evidence = session.advance(&context, 0).unwrap().row.unwrap();
    assert_eq!(evidence.status, Some(Status::ConsistentUnderAssumptions));
    let declaration = evidence.declaration_revision.unwrap();
    (session.into_capture(&context).unwrap(), declaration)
}
pub fn released() -> (release::Commit, EventBinding, String) {
    let (saved, _) = captured();
    let mut session = release::ReleaseSession::begin(
        saved.workspace,
        release::Request {
            batch_hash: saved.attachment_hash,
            rows: vec![release::Selection {
                construct_id: "native-circle".into(),
                options: cloning_outputs::Options {
                    name: "Released bottom10".into(),
                    show_prepared: false,
                    save_primers: false,
                    allow_incomplete: false,
                },
            }],
        },
    )
    .unwrap();
    let context = session.summary().context_revision;
    session.advance(&context, 0).unwrap();
    let row = session.inspect(0).unwrap().row;
    (
        session.into_commit(&context).unwrap(),
        row.output.construct,
        row.output.run.run_id,
    )
}
