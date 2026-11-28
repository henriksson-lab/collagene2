//! Named saved material feeds the ordinary shared compiler through fresh origin replay.
use collagene_core::{
    chemistry::{EndTreatment, FivePrime, Terminus},
    cloning::{
        assembly::recipe::{FragmentSpec, Preparation, RunPreparation},
        ligation::recipe as lig,
    },
    documents::workspace::Workspace,
    material::{
        capsule::Capsule,
        lineage::{EventBinding, InputEventBinding},
        outputs::{self as cloning_outputs, Options},
        run::Selection,
    },
    planning::plan::{DesignPlan, Operation, OrientedSource, PlanInput, Source},
    sequence::analysis::Primer,
    verification::evidence::Status,
};
use collagene_model::{Molecule, Topology};
fn binding(w: &Workspace, i: usize) -> EventBinding {
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
fn saved(source: Molecule, target: Molecule, steps: Vec<Operation>) -> (Molecule, Capsule) {
    let plan = DesignPlan {
        name: "Independent original preparation".into(),
        target_revision: target.revision(),
        inputs: vec![input("original", source.clone())],
        steps,
    };
    let w = Workspace::from_documents("Retained origin", vec![source, target], 1).unwrap();
    let commit = cloning_outputs::finalize(
        &w,
        &plan,
        binding(&w, 1),
        vec![InputEventBinding {
            input_id: "original".into(),
            event: binding(&w, 0),
        }],
        &Options {
            name: "Named reusable material".into(),
            show_prepared: false,
            save_primers: false,
            allow_incomplete: false,
        },
    )
    .unwrap();
    let event = commit.receipt.construct.clone();
    let m = commit.workspace.revisions[&event.revision].clone();
    let origin = collagene_core::material::capsule::capsule(
        &commit.workspace,
        &Selection {
            run_id: commit.receipt.run.run_id,
            step_index: plan.steps.len() - 1,
            event,
        },
    )
    .unwrap();
    (m, origin)
}
fn pcr_saved() -> (Molecule, Capsule) {
    saved(
        Molecule::dna("Original6bp", "AAACCC", Topology::Linear).unwrap(),
        Molecule::dna("Independent6bp", "AAACCC", Topology::Linear).unwrap(),
        vec![Operation::Pcr {
            template: Source::Input("original".into()),
            forward: Primer::new("F", "AAA", "").unwrap(),
            reverse: Primer::new("R", "GGG", "").unwrap(),
            template_strand: collagene_core::chemistry::Strand::Top,
        }],
    )
}
fn circle_saved() -> (Molecule, Capsule) {
    saved(
        Molecule::dna("Original12bp", "GGAATTCTTACG", Topology::Linear).unwrap(),
        Molecule::dna("Independent12bp circle", "GGAATTCTTACG", Topology::Circular).unwrap(),
        vec![
            Operation::DeclareSourceTermini {
                template: Source::Input("original".into()),
                left: Terminus::CLEAVED,
                right: Terminus::CLEAVED,
            },
            Operation::Ligate {
                fragments: vec![OrientedSource {
                    source: Source::Step(0),
                    reverse: false,
                }],
                circular: true,
                allow_nicked: false,
                options: Default::default(),
            },
        ],
    )
}
fn row(id: &str, preparation: Preparation, treatment: Option<EndTreatment>) -> lig::Fragment {
    lig::Fragment {
        source: FragmentSpec {
            generation: None,
            role: format!("{id}-role"),
            input_id: id.into(),
            preparation,
            reverse: false,
        },
        treatment,
    }
}
fn selected(name: &str) -> collagene_core::cloning::restriction::catalog::CatalogSelection {
    let c = collagene_core::cloning::restriction::catalog::builtin_catalog().unwrap();
    let id = c
        .edition()
        .enzymes
        .iter()
        .find(|e| e.name == name)
        .unwrap()
        .id
        .clone();
    collagene_core::cloning::restriction::catalog::select(&c, &[id]).unwrap()
}
#[test]
fn saved_circle_is_freshly_digested_and_ligated_with_two_real_oligos() {
    let (m, origin) = circle_saved();
    let inputs = vec![
        input("saved", m),
        input(
            "top",
            Molecule::dna("Top", "AATTGCTACG", Topology::Linear).unwrap(),
        ),
        input(
            "bottom",
            Molecule::dna("Bottom", "AATTCGTAGC", Topology::Linear).unwrap(),
        ),
    ];
    let target = Molecule::dna(
        "Independent22bp circle",
        "AATTCTTACGGGAATTGCTACG",
        Topology::Circular,
    )
    .unwrap();
    let recipe = lig::Recipe {
        name: "Serial digest and ligation".into(),
        circular: true,
        allow_nicked: false,
        fragments: vec![
            row(
                "saved",
                Preparation::FromRun {
                    origin,
                    preparation: RunPreparation::CatalogDigest {
                        conditions: None,
                        selection: selected("EcoRI"),
                        fragment_index: 0,
                        repair: None,
                    },
                },
                None,
            ),
            row(
                "top",
                Preparation::AnnealedOligos {
                    bottom_input_id: "bottom".into(),
                    options: collagene_core::chemistry::AnnealOptions {
                        min_paired_bases: 6,
                        ..Default::default()
                    },
                    candidate_index: 0,
                },
                Some(EndTreatment::Phosphorylate),
            ),
        ],
    };
    let preview = lig::preview(
        &inputs,
        &target,
        &recipe,
        collagene_core::cloning::assembly::recipe::CleavagePolicy::AssumedComplete,
    )
    .unwrap();
    assert_eq!(
        preview.verification.status,
        Status::ConsistentUnderAssumptions,
        "{:?}",
        preview.verification.findings
    );
    assert_eq!(
        preview
            .verification
            .final_product()
            .expect("final declared step must reconstruct")
            .sequence
            .as_str(),
        "AATTCTTACGGGAATTGCTACG"
    );
    assert!(matches!(
        preview.compiled.plan.steps[0],
        Operation::ReconstructRunProduct { .. }
    ));
    assert!(matches!(
        preview.compiled.plan.steps[1],
        Operation::CatalogDigest {
            template: Source::Step(0),
            ..
        }
    ));
    assert_eq!(
        preview.compiled.rows[0].preparation.preparation_step,
        Some(0)
    );
    let diagnostics = collagene_core::cloning::assembly::diagnostics::ligation(&preview).unwrap();
    assert_eq!(diagnostics.rows[0].total_steps, 2);
    assert!(preview.product_physical.unwrap().covalently_closed_in_model);
}
#[test]
fn saved_material_pcr_and_whole_material_keep_real_oh_ends_until_treatment() {
    for pcr in [false, true] {
        let (m, origin) = pcr_saved();
        let inputs = vec![
            input("saved", m),
            input(
                "other",
                Molecule::dna("Other6bp", "GGGTTT", Topology::Linear).unwrap(),
            ),
        ];
        let target =
            Molecule::dna("Independent12bp circle", "AAACCCGGGTTT", Topology::Circular).unwrap();
        let child = if pcr {
            RunPreparation::ExactPcr {
                forward: Primer::new("NewF", "AAA", "").unwrap(),
                reverse: Primer::new("NewR", "GGG", "").unwrap(),
                template_strand: collagene_core::chemistry::Strand::Top,
            }
        } else {
            RunPreparation::WholeLinear
        };
        let mut recipe = lig::Recipe {
            name: "Serial material reuse".into(),
            circular: true,
            allow_nicked: false,
            fragments: vec![
                row(
                    "saved",
                    Preparation::FromRun {
                        origin,
                        preparation: child,
                    },
                    None,
                ),
                row(
                    "other",
                    Preparation::DeclaredDuplex {
                        left: Terminus::CLEAVED,
                        right: Terminus::CLEAVED,
                    },
                    None,
                ),
            ],
        };
        let untreated = lig::preview(
            &inputs,
            &target,
            &recipe,
            collagene_core::cloning::assembly::recipe::CleavagePolicy::AssumedComplete,
        )
        .unwrap();
        assert!(untreated.product_physical.is_none());
        assert_eq!(
            untreated.physical[0]
                .as_ref()
                .unwrap()
                .ends
                .as_ref()
                .unwrap()
                .0
                .chemistry
                .five_prime,
            FivePrime::Hydroxyl
        );
        recipe.fragments[0].treatment = Some(EndTreatment::Phosphorylate);
        let preview = lig::preview(
            &inputs,
            &target,
            &recipe,
            collagene_core::cloning::assembly::recipe::CleavagePolicy::AssumedComplete,
        )
        .unwrap();
        assert_eq!(
            preview.verification.status,
            Status::ConsistentUnderAssumptions,
            "{:?}",
            preview.verification.findings
        );
        assert!(
            preview
                .product_physical
                .as_ref()
                .unwrap()
                .covalently_closed_in_model
        );
        preview
            .validate_binding(
                &lig::preview_revision(
                    &inputs,
                    &target,
                    &recipe,
                    collagene_core::cloning::assembly::recipe::CleavagePolicy::AssumedComplete,
                )
                .unwrap(),
            )
            .unwrap();
        collagene_core::cloning::assembly::diagnostics::ligation(&preview).unwrap();
    }
}
#[test]
fn stale_origin_and_cross_row_or_skipped_replay_chain_cannot_pass_readonly_validation() {
    let (m, origin) = pcr_saved();
    let inputs = vec![input("saved", m)];
    let target = Molecule::dna("Independent reference", "AAACCC", Topology::Linear).unwrap();
    let mut recipe = lig::Recipe {
        name: "Serial PCR".into(),
        circular: false,
        allow_nicked: false,
        fragments: vec![row(
            "saved",
            Preparation::FromRun {
                origin,
                preparation: RunPreparation::ExactPcr {
                    forward: Primer::new("F", "AAA", "").unwrap(),
                    reverse: Primer::new("R", "GGG", "").unwrap(),
                    template_strand: collagene_core::chemistry::Strand::Top,
                },
            },
            None,
        )],
    };
    let preview = lig::preview(
        &inputs,
        &target,
        &recipe,
        collagene_core::cloning::assembly::recipe::CleavagePolicy::AssumedComplete,
    )
    .unwrap();
    assert_eq!(
        preview.verification.status,
        Status::ConsistentUnderAssumptions
    );
    let mut forged = preview.clone();
    if let Operation::Pcr { template, .. } = &mut forged.compiled.plan.steps[1] {
        *template = Source::Input("saved".into());
    }
    forged.verification.plan_revision =
        collagene_core::planning::plan::DesignPlan::revision(&forged.compiled.plan);
    assert!(forged.validate().is_err());
    if let Preparation::FromRun { origin, .. } = &mut recipe.fragments[0].source.preparation {
        origin.selection.event.revision = "0".repeat(64);
    }
    assert!(
        lig::compile(
            &inputs,
            &target,
            &recipe,
            collagene_core::cloning::assembly::recipe::CleavagePolicy::AssumedComplete
        )
        .is_err()
    );
}

#[test]
fn replay_pcr_catalog_repair_tail_and_treatment_preserve_the_full_typed_chain() {
    use collagene_core::{
        chemistry::EndRepairMethod,
        cloning::ligation::{tailing::TailingBase, tailing_recipe as tail},
    };
    let (m, origin) = pcr_saved();
    let inputs = vec![
        input("saved", m),
        input(
            "other",
            Molecule::dna("Other6bp", "GGGTTT", Topology::Linear).unwrap(),
        ),
    ];
    // SmaI sites in the two literal tails cut CCC|GGG at boundaries3 and15.
    // The selected original top calls are GGG + AAACCC + CCC = GGGAAACCCCCC.
    // Add T, then the independently specified GGGTTT insert with its A tail.
    let target = Molecule::dna(
        "Independent20bp tailing circle",
        "GGGAAACCCCCCTGGGTTTA",
        Topology::Circular,
    )
    .unwrap();
    let recipe = tail::Recipe {
        profile: collagene_core::cloning::ligation::tailing_recipe::Profile::Ta,

        name: "Serial prepared material with complete downstream chain".into(),
        circular: true,
        allow_nicked: false,
        fragments: vec![
            tail::Fragment {
                source: FragmentSpec {
                    generation: None,
                    role: "reused".into(),
                    input_id: "saved".into(),
                    reverse: false,
                    preparation: Preparation::FromRun {
                        origin,
                        preparation: RunPreparation::PcrCatalogDigest {
                            conditions: None,
                            forward: Primer::new("AddedSiteF", "AAA", "CCCGGG").unwrap(),
                            reverse: Primer::new("AddedSiteR", "GGG", "CCCGGG").unwrap(),
                            selection: selected("SmaI"),
                            fragment_index: 1,
                            repair: Some(EndRepairMethod::Blunt),
                            template_strand: collagene_core::chemistry::Strand::Top,
                        },
                    },
                },
                action: collagene_core::cloning::ligation::tailing_recipe::Action::AddTail {
                    base: TailingBase::T,
                },
                treatment: None,
            },
            tail::Fragment {
                source: FragmentSpec {
                    generation: None,
                    role: "new-insert".into(),
                    input_id: "other".into(),
                    reverse: false,
                    preparation: Preparation::ExactPcr {
                        forward: Primer::new("NewF", "GGG", "").unwrap(),
                        reverse: Primer::new("NewR", "AAA", "").unwrap(),
                        template_strand: collagene_core::chemistry::Strand::Top,
                    },
                },
                action: collagene_core::cloning::ligation::tailing_recipe::Action::AddTail {
                    base: TailingBase::A,
                },
                treatment: Some(EndTreatment::Phosphorylate),
            },
        ],
    };
    let preview = tail::preview(&inputs, &target, &recipe).unwrap();
    assert_eq!(
        preview.verification.status,
        Status::ConsistentUnderAssumptions,
        "{:?}",
        preview.verification.findings
    );
    assert_eq!(
        preview
            .verification
            .final_product()
            .expect("final declared step must reconstruct")
            .sequence
            .as_str(),
        "GGGAAACCCCCCTGGGTTTA"
    );
    assert_eq!(
        preview.compiled.rows[0].preparation.preparation_step,
        Some(0)
    );
    assert_eq!(preview.compiled.rows[0].preparation.repair_step, Some(3));
    assert!(matches!(
        preview.compiled.plan.steps[2],
        Operation::CatalogDigest {
            template: Source::Step(1),
            ..
        }
    ));
    assert!(matches!(
        preview.compiled.plan.steps[3],
        Operation::RepairEnds {
            template: Source::Step(2),
            ..
        }
    ));
    preview
        .validate_binding(&tail::preview_revision(&inputs, &target, &recipe).unwrap())
        .unwrap();
    let mut forged = preview.clone();
    if let Operation::RepairEnds { template, .. } = &mut forged.compiled.plan.steps[3] {
        *template = Source::Step(0);
    }
    forged.verification.plan_revision =
        collagene_core::planning::plan::DesignPlan::revision(&forged.compiled.plan);
    assert!(forged.validate().is_err());
}
