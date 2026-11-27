//! The batch ledger reserves the same declared PCR/OE work as native recipes.
//! Oversized declarations are rejected by preflight; no costly PCR is executed.
use collagene_core::{
    Molecule, Topology,
    analysis::Primer,
    assembly_recipe::{Preparation, RunPreparation},
    chemistry::Terminus,
    cloning_outputs, combinatorial_batch as batch, combinatorial_design as library,
    material_derivation as derivation, overlap_extension, overlap_extension_pcr,
    run_lineage::{EventBinding, InputEventBinding},
    run_material,
    verification::{self, DesignPlan, Operation, PlanInput, Source},
    workspace::Workspace,
};
fn event(w: &Workspace, index: usize) -> EventBinding {
    let d = &w.documents[index];
    EventBinding {
        document_id: d.id.clone(),
        event_id: d.head_event.clone().unwrap(),
        revision: d.head.clone(),
    }
}
fn oe() -> overlap_extension_pcr::Request {
    let pair = overlap_extension_pcr::PrimerPair {
        forward: Primer::new("Declared F", &"A".repeat(1000), "").unwrap(),
        reverse: Primer::new("Declared R", &"T".repeat(1000), "").unwrap(),
        template_strand: collagene_core::chemistry::Strand::Top,
    };
    overlap_extension_pcr::Request {
        left: pair.clone(),
        right: pair.clone(),
        outer: pair,
        extension: overlap_extension::Options {
            overlap: 8,
            min_paired_bases: 8,
            ..Default::default()
        },
    }
}
fn stock() -> String {
    format!("AGTCCGAT{}GATCGTCA", "AC".repeat(1992))
}
fn source(family: usize) -> (Workspace, EventBinding, Preparation) {
    let dna = stock();
    assert_eq!(dna.len(), 4000);
    let m = Molecule::dna("Actual original stock", &dna, Topology::Linear).unwrap();
    let w = Workspace::from_documents(
        "Native cost sources",
        vec![
            m.clone(),
            Molecule::dna("Independent stock reference", &dna, Topology::Linear).unwrap(),
            Molecule::dna("Small other source", "CGTA", Topology::Linear).unwrap(),
        ],
        0,
    )
    .unwrap();
    let (w, pin, preparation) = match family {
        0 => {
            let pin = event(&w, 0);
            (w, pin, Preparation::OverlapExtensionPcr { request: oe() })
        }
        1 => {
            let input = PlanInput {
                id: "stock".into(),
                revision: m.revision(),
                molecule: m,
            };
            let plan = DesignPlan {
                name: "Actual native source termini".into(),
                target_revision: event(&w, 1).revision.clone(),
                inputs: vec![input],
                steps: vec![Operation::DeclareSourceTermini {
                    template: Source::Input("stock".into()),
                    left: Terminus::CLEAVED,
                    right: Terminus::CLEAVED,
                }],
            };
            let commit = cloning_outputs::finalize(
                &w,
                &plan,
                event(&w, 1),
                vec![InputEventBinding {
                    input_id: "stock".into(),
                    event: event(&w, 0),
                }],
                &cloning_outputs::Options {
                    name: "Retained native stock".into(),
                    show_prepared: false,
                    save_primers: false,
                    allow_incomplete: false,
                },
            )
            .unwrap();
            let pin = commit.receipt.construct;
            let origin = run_material::capsule(
                &commit.workspace,
                &run_material::Selection {
                    run_id: commit.receipt.run.run_id,
                    step_index: 0,
                    event: pin.clone(),
                },
            )
            .unwrap();
            (
                commit.workspace,
                pin,
                Preparation::FromRun {
                    origin,
                    preparation: RunPreparation::OverlapExtensionPcr { request: oe() },
                },
            )
        }
        2 => {
            let forward = Primer::new("Unique F", "AGTCCGAT", "").unwrap();
            let reverse = Primer::new("Unique R", "TGACGATC", "").unwrap();
            let previews = verification::simulate_pcr_documents(
                &m,
                &forward,
                &reverse,
                collagene_core::chemistry::Strand::Top,
            )
            .unwrap();
            assert_eq!(previews.len(), 1);
            let request = derivation::Request {
                source: event(&w, 0),
                origin: derivation::SourceOrigin::OriginalReference,
                derivation: derivation::Derivation::ExactPcr {
                    forward,
                    reverse,
                    selection: derivation::PcrSelection {
                        selected_index: 0,
                        expected_preview_revision: previews[0].molecule.revision(),
                        geometry: previews[0].geometry.clone(),
                    },
                    template_strand: collagene_core::chemistry::Strand::Top,
                },
                name: "Actual native PCR stock".into(),
                save_primers: false,
            };
            let installed = derivation::prepare_install(&w, &request).unwrap();
            let pin = installed.receipt.product;
            let origin = derivation::capsule(&installed.workspace, &pin).unwrap();
            (
                installed.workspace,
                pin,
                Preparation::FromDerived {
                    origin,
                    preparation: RunPreparation::OverlapExtensionPcr { request: oe() },
                },
            )
        }
        _ => unreachable!(),
    };
    (w, pin, preparation)
}
#[test]
fn original_run_and_derived_oe_share_native_cost_and_stop_the_batch_before_unreserved_work() {
    for family in 0..3 {
        let (w, pin, preparation) = source(family);
        let before = serde_json::to_vec(&w).unwrap();
        let material_bindings = match &preparation {
            Preparation::FromRun { origin, .. } => vec![library::LibraryMaterialBinding {
                event: pin.clone(),
                origin: run_material::MaterialOrigin::Run(origin.clone()),
            }],
            Preparation::FromDerived { origin, .. } => vec![library::LibraryMaterialBinding {
                event: pin.clone(),
                origin: run_material::MaterialOrigin::Derived(origin.clone()),
            }],
            _ => vec![],
        };
        let d = library::Design {
            material_bindings,
            id: "cost-library".into(),
            name: "Bounded actual preparation costs".into(),
            circular: false,
            method: library::Method::Ligation {
                allow_nicked: false,
            },
            bins: vec![
                library::Bin {
                    id: "large".into(),
                    label: "Large source".into(),
                    choices: vec![library::Choice {
                        id: "large-choice".into(),
                        label: "Explicit oversized OE".into(),
                        source: pin,
                        selector: collagene_core::assembly_intent::Selector::ReferenceRange {
                            start: 0,
                            end: 4000,
                            reverse: false,
                        },
                        preparation_reverse: false,
                        generation: None,
                        treatment: None,
                        preparation: Some(preparation),
                    }],
                },
                library::Bin {
                    id: "small".into(),
                    label: "Other".into(),
                    choices: vec![library::Choice {
                        id: "small-choice".into(),
                        label: "Other".into(),
                        source: event(&w, 2),
                        selector: collagene_core::assembly_intent::Selector::ReferenceRange {
                            start: 0,
                            end: 4,
                            reverse: false,
                        },
                        preparation_reverse: false,
                        generation: None,
                        treatment: None,
                        preparation: Some(Preparation::WholeLinear),
                    }],
                },
            ],
        };
        let rows = (0..20)
            .map(|i| {
                batch::RowRequest::Ready(library::SelectedConstruct {
                    id: format!("row-{i}"),
                    name: format!("Bounded declaration {i}"),
                    choices: vec![
                        library::ChoiceSelection {
                            bin_id: "large".into(),
                            choice_id: "large-choice".into(),
                        },
                        library::ChoiceSelection {
                            bin_id: "small".into(),
                            choice_id: "small-choice".into(),
                        },
                    ],
                    target: library::TargetDeclaration::OriginalRegions {
                        junction_edits: vec![],
                        substitutions: vec![],
                        edits: vec![],
                    },
                    overrides: None,
                })
            })
            .collect();
        let mut session = batch::Session::begin(
            w.clone(),
            batch::Request {
                design: Some(d),
                rows,
            },
        )
        .unwrap();
        let context = session.summary().context_revision;
        let mut first_withheld = None;
        for i in 0..20 {
            let row = session.advance(&context, i).unwrap().row.unwrap();
            assert_eq!(row.state, batch::RowState::Error, "family {family},row {i}");
            assert!(
                row.status.is_none(),
                "preparation is rejected before native PCR"
            );
            let withheld = row
                .error
                .as_ref()
                .unwrap()
                .contains("Not computed: aggregate batch work budget");
            if withheld && first_withheld.is_none() {
                first_withheld = Some(i);
            }
            if first_withheld.is_some() {
                assert!(withheld, "budget remains exhausted");
            }
        }
        // Each native OE request reserves >127M declared search units and is
        // individually rejected by the recipe's 100M preparation limit. The
        // first15 consume <2B aggregate units; the sixteenth must be withheld.
        assert_eq!(
            first_withheld,
            Some(15),
            "wrapper family {family} must not hide child work"
        );
        assert_eq!(serde_json::to_vec(&w).unwrap(), before);
    }
}
