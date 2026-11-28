//! Independently declared scars, saved actual tails, and explicit per-row decisions.
use collagene_core::{
    chemistry::{EndTreatment, Polarity},
    cloning::{
        assembly::recipe::{FragmentSpec, Preparation, RunPreparation},
        ligation::{
            tailing::TailingBase,
            tailing_recipe::{self as t, Action, Profile, ProfileStatus},
        },
    },
    documents::workspace::Workspace,
    material::{
        lineage::{self as run_lineage, EventBinding, InputEventBinding},
        run::Selection,
    },
    planning::plan::{DesignPlan, Operation, PlanInput, Source},
    sequence::analysis::Primer,
    verification::{self, evidence::Status},
};
use collagene_model::{Molecule, Topology};
fn input(id: &str, molecule: Molecule) -> PlanInput {
    PlanInput {
        id: id.into(),
        revision: molecule.revision(),
        molecule,
    }
}
fn pcr(f: &str, r: &str) -> Preparation {
    Preparation::ExactPcr {
        forward: Primer::new("F", f, "").unwrap(),
        reverse: Primer::new("R", r, "").unwrap(),
        template_strand: collagene_core::chemistry::Strand::Top,
    }
}
fn row(
    id: &str,
    preparation: Preparation,
    action: Action,
    treatment: Option<EndTreatment>,
) -> t::Fragment {
    t::Fragment {
        source: FragmentSpec {
            generation: None,
            role: id.into(),
            input_id: id.into(),
            preparation,
            reverse: false,
        },
        action,
        treatment,
    }
}
fn event(w: &Workspace, i: usize) -> EventBinding {
    let d = &w.documents[i];
    EventBinding {
        document_id: d.id.clone(),
        event_id: d.head_event.clone().unwrap(),
        revision: d.head.clone(),
    }
}
fn saved_tail(
    letters: &str,
    base: TailingBase,
    f: &str,
    r: &str,
) -> (Molecule, collagene_core::material::capsule::Capsule) {
    let original = Molecule::dna("Original stock", letters, Topology::Linear).unwrap();
    let target = Molecule::dna(
        "Independently declared one-base tailed stock",
        &format!("{letters}{}", base.as_str()),
        Topology::Linear,
    )
    .unwrap();
    let w = Workspace::from_documents(
        "Native tailed stock",
        vec![original.clone(), target.clone()],
        1,
    )
    .unwrap();
    let plan = DesignPlan {
        standard: None,
        name: "Historical stock preparation".into(),
        target_revision: target.revision(),
        inputs: vec![input("stock", original)],
        steps: vec![
            Operation::Pcr {
                template: Source::Input("stock".into()),
                forward: Primer::new("Historical F", f, "").unwrap(),
                reverse: Primer::new("Historical R", r, "").unwrap(),
                template_strand: collagene_core::chemistry::Strand::Top,
            },
            Operation::Tail3Prime {
                template: Source::Step(0),
                base,
            },
            Operation::TreatEnds {
                template: Source::Step(1),
                treatment: EndTreatment::Phosphorylate,
            },
        ],
    };
    let commit = run_lineage::record_run(
        &w,
        &plan,
        event(&w, 1),
        vec![InputEventBinding {
            input_id: "stock".into(),
            event: event(&w, 0),
        }],
        false,
    )
    .unwrap();
    let product = &commit.receipt.products[2];
    let selected = Selection {
        run_id: commit.receipt.run_id.clone(),
        step_index: 2,
        event: EventBinding {
            document_id: product.document_id.clone(),
            event_id: product.event_id.clone(),
            revision: product.revision.clone(),
        },
    };
    let source = commit.workspace.revisions[&product.revision].clone();
    let origin = collagene_core::material::capsule::capsule(&commit.workspace, &selected).unwrap();
    (source, origin)
}
fn same_circle(a: &str, b: &str) -> bool {
    a.len() == b.len() && format!("{a}{a}").contains(b)
}
#[test]
fn keep_saved_t_vector_and_add_a_insert_preserves_actual_material_and_only_new_oligos() {
    let (vector, origin) = saved_tail("AAACCC", TailingBase::T, "AAA", "GGG");
    let saved_modifications = vector.metadata.get("tracked_modifications").cloned();
    let inputs = vec![
        input("vector", vector),
        input(
            "insert",
            Molecule::dna("Insert stock", "GGGTTT", Topology::Linear).unwrap(),
        ),
    ];
    let target = Molecule::dna(
        "Independent mixed TA circle",
        "AAACCCTGGGTTTA",
        Topology::Circular,
    )
    .unwrap();
    let recipe = t::Recipe {
        standard: None,
        name: "Keep actual vector tails".into(),
        profile: Profile::Ta,
        fragments: vec![
            row(
                "vector",
                Preparation::FromRun {
                    origin,
                    preparation: RunPreparation::WholeLinear,
                },
                Action::KeepExisting,
                None,
            ),
            row(
                "insert",
                pcr("GGG", "AAA"),
                Action::AddTail {
                    base: TailingBase::A,
                },
                Some(EndTreatment::Phosphorylate),
            ),
        ],
        circular: true,
        allow_nicked: false,
    };
    let before = serde_json::to_vec(&inputs).unwrap();
    let p = t::preview(&inputs, &target, &recipe).unwrap();
    assert_eq!(
        p.profile_status,
        ProfileStatus::Matched,
        "{:?}",
        p.profile_findings
    );
    assert_eq!(
        p.verification.status,
        Status::ConsistentUnderAssumptions,
        "{:?}",
        p.verification.findings
    );
    assert!(p.compiled.rows[0].tail_step.is_none());
    assert!(p.compiled.rows[1].tail_step.is_some());
    assert_eq!(
        p.compiled
            .plan
            .steps
            .iter()
            .filter(|o| matches!(o, Operation::Tail3Prime { .. }))
            .count(),
        1
    );
    assert_eq!(
        p.prepared[0]
            .as_ref()
            .expect("prepared row must reconstruct")
            .sequence
            .as_str(),
        "AAACCCT"
    );
    assert_eq!(
        p.prepared[0]
            .as_ref()
            .expect("prepared row must reconstruct")
            .metadata
            .get("tracked_modifications")
            .cloned(),
        saved_modifications
    );
    assert!(same_circle(
        p.verification
            .final_product()
            .expect("final declared step must reconstruct")
            .sequence
            .as_str(),
        "AAACCCTGGGTTTA"
    ));
    assert!(
        p.product_physical
            .as_ref()
            .unwrap()
            .covalently_closed_in_model
    );
    let oligos =
        collagene_core::planning::oligos::build(&p.compiled.plan, &target, &Default::default())
            .unwrap();
    assert_eq!(
        oligos.oligos.records.len(),
        2,
        "historical vector PCR oligos are not new orders"
    );
    for (i, base) in ["T", "A"].into_iter().enumerate() {
        let (left, right) = p.physical[i].as_ref().unwrap().ends.as_ref().unwrap();
        assert_eq!(left.polarity, Polarity::ThreePrime);
        assert_eq!(left.overhang, base);
        assert_eq!(right.overhang, base);
    }
    p.validate_for(
        &inputs,
        &target,
        &recipe,
        &t::preview_revision(&inputs, &target, &recipe).unwrap(),
    )
    .unwrap();
    assert_eq!(before, serde_json::to_vec(&inputs).unwrap());
    let mut changed = p.clone();
    changed.compiled.rows[0].tail_step = Some(0);
    assert!(
        changed
            .validate_for(&inputs, &target, &recipe, &p.context_revision)
            .is_err()
    );
    // The separately declared alternative target never chooses the insert orientation.
    let opposite_target = Molecule::dna(
        "Independent opposite orientation",
        "AAACCCTAAACCCA",
        Topology::Circular,
    )
    .unwrap();
    let wrong = t::preview(&inputs, &opposite_target, &recipe).unwrap();
    assert_eq!(wrong.verification.status, Status::Invalid);
    assert!(same_circle(
        wrong
            .verification
            .final_product()
            .expect("final declared step must reconstruct")
            .sequence
            .as_str(),
        "AAACCCTGGGTTTA"
    ));
    let outcomes = wrong
        .verification
        .physical_processing_reports
        .iter()
        .find_map(|step| match &step.evidence {
            verification::evidence::PhysicalProcessingEvidence::LigationOutcomes(e) => Some(e),
            _ => None,
        })
        .unwrap();
    assert!(
        outcomes
            .report
            .alternatives
            .iter()
            .any(|product| same_circle(product.sequence.as_str(), "AAACCCTAAACCCA"))
    );
    let mut explicitly_reversed = recipe.clone();
    explicitly_reversed.fragments[1].source.reverse = true;
    let reverse = t::preview(&inputs, &opposite_target, &explicitly_reversed).unwrap();
    assert_eq!(reverse.profile_status, ProfileStatus::Matched);
    assert_eq!(
        reverse.verification.status,
        Status::ConsistentUnderAssumptions,
        "{:?}",
        reverse.verification.findings
    );
    assert!(same_circle(
        reverse
            .verification
            .final_product()
            .expect("final declared step must reconstruct")
            .sequence
            .as_str(),
        "AAACCCTAAACCCA"
    ));
    // Correct tail geometry alone does not supply the insert's missing 5′ phosphate.
    let mut unphosphorylated = recipe.clone();
    unphosphorylated.fragments[1].treatment = None;
    let bad = t::preview(&inputs, &target, &unphosphorylated).unwrap();
    assert_eq!(bad.profile_status, ProfileStatus::Matched);
    assert_ne!(bad.verification.status, Status::ConsistentUnderAssumptions);
    assert!(bad.product_physical.is_none());
}
#[test]
fn both_existing_saved_tails_do_not_repeat_tailing_or_order_upstream_primers() {
    let (vector, vo) = saved_tail("AAACCC", TailingBase::T, "AAA", "GGG");
    let (insert, io) = saved_tail("GGGTTT", TailingBase::A, "GGG", "AAA");
    let inputs = vec![input("vector", vector), input("insert", insert)];
    let target = Molecule::dna(
        "Independent both-existing TA circle",
        "AAACCCTGGGTTTA",
        Topology::Circular,
    )
    .unwrap();
    let recipe = t::Recipe {
        standard: None,
        name: "Hold both actual tailed substrates".into(),
        profile: Profile::Ta,
        fragments: vec![
            row(
                "vector",
                Preparation::FromRun {
                    origin: vo,
                    preparation: RunPreparation::WholeLinear,
                },
                Action::KeepExisting,
                None,
            ),
            row(
                "insert",
                Preparation::FromRun {
                    origin: io,
                    preparation: RunPreparation::WholeLinear,
                },
                Action::KeepExisting,
                None,
            ),
        ],
        circular: true,
        allow_nicked: false,
    };
    let p = t::preview(&inputs, &target, &recipe).unwrap();
    assert_eq!(p.profile_status, ProfileStatus::Matched);
    assert_eq!(
        p.verification.status,
        Status::ConsistentUnderAssumptions,
        "{:?}",
        p.verification.findings
    );
    assert!(p.compiled.rows.iter().all(|r| r.tail_step.is_none()));
    assert_eq!(p.compiled.plan.steps.len(), 3);
    assert!(
        !p.compiled
            .plan
            .steps
            .iter()
            .any(|op| matches!(op, Operation::Tail3Prime { .. } | Operation::Pcr { .. }))
    );
    assert!(
        collagene_core::planning::oligos::build(&p.compiled.plan, &target, &Default::default())
            .unwrap()
            .oligos
            .records
            .is_empty()
    );
    assert!(same_circle(
        p.verification
            .final_product()
            .expect("final declared step must reconstruct")
            .sequence
            .as_str(),
        "AAACCCTGGGTTTA"
    ));
}
#[test]
fn all_four_explicit_bases_match_independent_ta_gc_scar_sequences() {
    for (profile, a, b, letters) in [
        (
            Profile::Ta,
            TailingBase::T,
            TailingBase::A,
            "AAACCCTGGGTTTA",
        ),
        (
            Profile::Gc,
            TailingBase::G,
            TailingBase::C,
            "AAACCCGGGGTTTC",
        ),
    ] {
        let inputs = vec![
            input(
                "vector",
                Molecule::dna("Vector", "AAACCC", Topology::Linear).unwrap(),
            ),
            input(
                "insert",
                Molecule::dna("Insert", "GGGTTT", Topology::Linear).unwrap(),
            ),
        ];
        let target = Molecule::dna("Independent scar", letters, Topology::Circular).unwrap();
        let recipe = t::Recipe {
            standard: None,
            name: "Parity".into(),
            profile,
            fragments: vec![
                row(
                    "vector",
                    pcr("AAA", "GGG"),
                    Action::AddTail { base: a },
                    Some(EndTreatment::Phosphorylate),
                ),
                row(
                    "insert",
                    pcr("GGG", "AAA"),
                    Action::AddTail { base: b },
                    Some(EndTreatment::Phosphorylate),
                ),
            ],
            circular: true,
            allow_nicked: false,
        };
        let p = t::preview(&inputs, &target, &recipe).unwrap();
        assert_eq!(p.profile_status, ProfileStatus::Matched);
        assert_eq!(p.verification.status, Status::ConsistentUnderAssumptions);
        assert!(same_circle(
            p.verification
                .final_product()
                .expect("final declared step must reconstruct")
                .sequence
                .as_str(),
            letters
        ));
        p.validate_for(&inputs, &target, &recipe, &p.context_revision)
            .unwrap();
    }
}
#[test]
fn keep_blunt_is_not_ta_even_when_generic_ligation_is_valid_and_wrong_family_refuses() {
    let inputs = vec![
        input(
            "vector",
            Molecule::dna("V", "AAACCC", Topology::Linear).unwrap(),
        ),
        input(
            "insert",
            Molecule::dna("I", "GGGTTT", Topology::Linear).unwrap(),
        ),
    ];
    let target = Molecule::dna(
        "Independent blunt circle",
        "AAACCCGGGTTT",
        Topology::Circular,
    )
    .unwrap();
    let recipe = t::Recipe {
        standard: None,
        name: "No invented tails".into(),
        profile: Profile::Ta,
        fragments: vec![
            row(
                "vector",
                pcr("AAA", "GGG"),
                Action::KeepExisting,
                Some(EndTreatment::Phosphorylate),
            ),
            row(
                "insert",
                pcr("GGG", "AAA"),
                Action::KeepExisting,
                Some(EndTreatment::Phosphorylate),
            ),
        ],
        circular: true,
        allow_nicked: false,
    };
    let p = t::preview(&inputs, &target, &recipe).unwrap();
    assert_eq!(p.verification.status, Status::ConsistentUnderAssumptions);
    assert_eq!(p.profile_status, ProfileStatus::Incompatible);
    assert!(!p.profile_findings.is_empty());
    let mut wrong = recipe.clone();
    wrong.fragments[1].action = Action::AddTail {
        base: TailingBase::G,
    };
    assert!(t::compile(&inputs, &target, &wrong).is_err());
    let mut missing = recipe.clone();
    missing.fragments[0].source.preparation = Preparation::WholeLinear;
    missing.fragments[0].treatment = None;
    let p = t::preview(&inputs, &target, &missing).unwrap();
    assert_ne!(p.profile_status, ProfileStatus::Matched);
    assert_ne!(p.verification.status, Status::ConsistentUnderAssumptions);
}

#[test]
fn keep_existing_replays_sticky_incomplete_origin_without_upgrading_material_assessment() {
    // Selected native homology is established with exactly enough work; terminal
    // alternatives remain unenumerated. The physical linear product still exists.
    let left = Molecule::dna("Left original", "AGTCCGATGCTA", Topology::Linear).unwrap();
    let right = Molecule::dna("Right original", "GCTATTCGAACC", Topology::Linear).unwrap();
    let stock_target = Molecule::dna(
        "Independent tailed linear stock",
        "AGTCCGATGCTATTCGAACCT",
        Topology::Linear,
    )
    .unwrap();
    let w = Workspace::from_documents(
        "Incomplete native homology stock",
        vec![left.clone(), right.clone(), stock_target.clone()],
        2,
    )
    .unwrap();
    let policy = collagene_core::cloning::overlap::assembly::Options {
        overlaps: vec![4],
        topology: Topology::Linear,
        min_paired_bases: 4,
        max_comparisons: 3 * 24 + 2 * 20 + 4,
        ..Default::default()
    };
    let plan = DesignPlan {
        standard: None,
        name: "Explicit finite-search original stock".into(),
        target_revision: stock_target.revision(),
        inputs: vec![input("left", left), input("right", right)],
        steps: vec![
            Operation::Pcr {
                template: Source::Input("left".into()),
                forward: Primer::new("LF", "AGTC", "").unwrap(),
                reverse: Primer::new("LR", "TAGC", "").unwrap(),
                template_strand: collagene_core::chemistry::Strand::Top,
            },
            Operation::Pcr {
                template: Source::Input("right".into()),
                forward: Primer::new("RF", "GCTA", "").unwrap(),
                reverse: Primer::new("RR", "GGTT", "").unwrap(),
                template_strand: collagene_core::chemistry::Strand::Top,
            },
            Operation::HomologyAssembly {
                fragments: vec![
                    collagene_core::planning::plan::OrientedSource {
                        source: Source::Step(0),
                        reverse: false,
                    },
                    collagene_core::planning::plan::OrientedSource {
                        source: Source::Step(1),
                        reverse: false,
                    },
                ],
                options: Default::default(),
                policy,
            },
            Operation::Tail3Prime {
                template: Source::Step(2),
                base: TailingBase::T,
            },
            Operation::TreatEnds {
                template: Source::Step(3),
                treatment: EndTreatment::Phosphorylate,
            },
        ],
    };
    let report = verification::engine::verify(&plan, &stock_target);
    assert_eq!(report.status, Status::Incomplete, "{:?}", report.findings);
    assert_eq!(report.reconstructed_count(), 5);
    assert_eq!(
        report
            .product(4)
            .expect("selected step must reconstruct")
            .sequence
            .as_str(),
        "AGTCCGATGCTATTCGAACCT"
    );
    let commit = run_lineage::record_run(
        &w,
        &plan,
        event(&w, 2),
        vec![
            InputEventBinding {
                input_id: "left".into(),
                event: event(&w, 0),
            },
            InputEventBinding {
                input_id: "right".into(),
                event: event(&w, 1),
            },
        ],
        true,
    )
    .unwrap();
    let product = &commit.receipt.products[4];
    let selected = Selection {
        run_id: commit.receipt.run_id.clone(),
        step_index: 4,
        event: EventBinding {
            document_id: product.document_id.clone(),
            event_id: product.event_id.clone(),
            revision: product.revision.clone(),
        },
    };
    let stock = commit.workspace.revisions[&product.revision].clone();
    let origin = collagene_core::material::capsule::capsule(&commit.workspace, &selected).unwrap();
    let inputs = vec![
        input("vector", stock),
        input(
            "insert",
            Molecule::dna("Original insert", "GGGTTT", Topology::Linear).unwrap(),
        ),
    ];
    let target = Molecule::dna(
        "Independent downstream TA circle",
        "AGTCCGATGCTATTCGAACCTGGGTTTA",
        Topology::Circular,
    )
    .unwrap();
    let recipe = t::Recipe {
        standard: None,
        name: "Hold incomplete actual tail source".into(),
        profile: Profile::Ta,
        fragments: vec![
            row(
                "vector",
                Preparation::FromRun {
                    origin,
                    preparation: RunPreparation::WholeLinear,
                },
                Action::KeepExisting,
                None,
            ),
            row(
                "insert",
                pcr("GGG", "AAA"),
                Action::AddTail {
                    base: TailingBase::A,
                },
                Some(EndTreatment::Phosphorylate),
            ),
        ],
        circular: true,
        allow_nicked: false,
    };
    let p = t::preview(&inputs, &target, &recipe).unwrap();
    assert_eq!(p.profile_status, ProfileStatus::Matched);
    assert_eq!(
        p.verification.status,
        Status::Incomplete,
        "{:?}",
        p.verification.findings
    );
    assert!(
        p.product_physical
            .as_ref()
            .unwrap()
            .covalently_closed_in_model
    );
    assert!(same_circle(
        p.verification
            .final_product()
            .expect("final declared step must reconstruct")
            .sequence
            .as_str(),
        "AGTCCGATGCTATTCGAACCTGGGTTTA"
    ));
    let origin_report = p
        .verification
        .physical_processing_reports
        .iter()
        .find_map(|step| match &step.evidence {
            verification::evidence::PhysicalProcessingEvidence::RunProduct(e) => Some(
                serde_json::from_slice::<verification::evidence::VerificationReport>(
                    &e.current_report,
                )
                .unwrap(),
            ),
            _ => None,
        })
        .unwrap();
    assert_eq!(origin_report.status, Status::Incomplete);
    p.validate_for(&inputs, &target, &recipe, &p.context_revision)
        .unwrap();
}

mod preparation_and_receipt_tests {
    use collagene_core::{
        chemistry::{EndTreatment, Polarity, Terminus, ThreePrime},
        cloning::{
            assembly::recipe::{FragmentSpec, Preparation},
            ligation::{tailing::TailingBase, tailing_recipe as t},
        },
        planning::plan::{Operation, PlanInput, Source},
        sequence::analysis::Primer,
        verification::evidence::Status,
    };
    use collagene_model::{Molecule, Topology};
    fn input(id: &str, letters: &str) -> PlanInput {
        let molecule = Molecule::dna(id, letters, Topology::Linear).unwrap();
        PlanInput {
            id: id.into(),
            revision: molecule.revision(),
            molecule,
        }
    }
    fn fragment(id: &str, base: TailingBase, forward: &str, reverse: &str) -> t::Fragment {
        t::Fragment {
            source: FragmentSpec {
                generation: None,
                role: format!("{id}-role"),
                input_id: id.into(),
                preparation: Preparation::ExactPcr {
                    forward: Primer::new("Forward", forward, "").unwrap(),
                    reverse: Primer::new("Reverse", reverse, "").unwrap(),
                    template_strand: collagene_core::chemistry::Strand::Top,
                },
                reverse: false,
            },
            action: collagene_core::cloning::ligation::tailing_recipe::Action::AddTail { base },
            treatment: Some(EndTreatment::Phosphorylate),
        }
    }
    fn fixture(
        a: TailingBase,
        b: TailingBase,
        target: &str,
    ) -> (Vec<PlanInput>, Molecule, t::Recipe) {
        (
            vec![input("vector", "AAACCC"), input("insert", "GGGTTT")],
            Molecule::dna("Independent scar intent", target, Topology::Circular).unwrap(),
            t::Recipe {
                standard: None,
                profile: if matches!(a, TailingBase::A | TailingBase::T) {
                    collagene_core::cloning::ligation::tailing_recipe::Profile::Ta
                } else {
                    collagene_core::cloning::ligation::tailing_recipe::Profile::Gc
                },
                name: "Original-source tailing recipe".into(),
                fragments: vec![
                    fragment("vector", a, "AAA", "GGG"),
                    fragment("insert", b, "GGG", "AAA"),
                ],
                circular: true,
                allow_nicked: false,
            },
        )
    }
    fn same_circle(actual: &str, expected: &str) -> bool {
        actual.len() == expected.len() && format!("{actual}{actual}").contains(expected)
    }
    #[test]
    fn ta_and_gc_close_against_literal_targets_with_real_tail_ends_and_two_junctions() {
        for (a, b, letters) in [
            (TailingBase::T, TailingBase::A, "AAACCCTGGGTTTA"),
            (TailingBase::G, TailingBase::C, "AAACCCGGGGTTTC"),
        ] {
            let (inputs, target, recipe) = fixture(a, b, letters);
            let before = serde_json::to_vec(&inputs).unwrap();
            let p = t::preview(&inputs, &target, &recipe).unwrap();
            p.validate_binding(&t::preview_revision(&inputs, &target, &recipe).unwrap())
                .unwrap();
            assert_eq!(
                p.verification.status,
                Status::ConsistentUnderAssumptions,
                "{:?}",
                p.verification.findings
            );
            assert!(p.preparation_complete);
            assert_eq!(p.physical.len(), 2);
            assert_eq!(p.junctions.len(), 2);
            assert!(
                p.product_physical
                    .as_ref()
                    .unwrap()
                    .covalently_closed_in_model
            );
            let actual = p
                .verification
                .final_product()
                .expect("final declared step must reconstruct");
            assert!(same_circle(actual.sequence.as_str(), letters));
            assert_eq!(actual.sequence.len(), 14);
            for (i, base) in [a, b].iter().enumerate() {
                let summary = p.physical[i].as_ref().unwrap();
                let (left, right) = summary.ends.as_ref().unwrap();
                assert_eq!(left.polarity, Polarity::ThreePrime);
                assert_eq!(right.polarity, Polarity::ThreePrime);
                assert_eq!(left.overhang, base.as_str());
                assert_eq!(right.overhang, base.as_str());
                let row = &p.compiled.rows[i];
                assert!(
                    matches!(&p.compiled.plan.steps[row.tail_step.unwrap()],Operation::Tail3Prime{base:actual,..}if actual==base)
                );
            }
            assert_eq!(serde_json::to_vec(&inputs).unwrap(), before);
            let imported: t::Preview =
                serde_json::from_slice(&serde_json::to_vec(&p).unwrap()).unwrap();
            imported.validate_binding(&p.context_revision).unwrap();
        }
    }
    #[test]
    fn wrong_target_incompatible_tails_and_unknown_or_blocked_original_ends_do_not_gain_success() {
        let (inputs, mut target, mut recipe) =
            fixture(TailingBase::T, TailingBase::A, "AAACCCTGGGTTTA");
        target.sequence = collagene_model::Sequence::dna("AAACCCTGGGTTTG").unwrap();
        let wrong = t::preview(&inputs, &target, &recipe).unwrap();
        assert_eq!(wrong.verification.status, Status::Invalid);
        assert!(wrong.preparation_complete);
        assert!(wrong.verification.final_product().is_some());
        recipe.fragments[1].action =
            collagene_core::cloning::ligation::tailing_recipe::Action::AddTail {
                base: TailingBase::T,
            };
        let incompatible = t::preview(&inputs, &target, &recipe).unwrap();
        assert_ne!(
            incompatible.verification.status,
            Status::ConsistentUnderAssumptions
        );
        assert!(incompatible.product_physical.is_none());
        for prep in [
            Preparation::WholeLinear,
            Preparation::DeclaredDuplex {
                left: Terminus {
                    five_prime: collagene_core::chemistry::FivePrime::Phosphate,
                    three_prime: ThreePrime::Blocked,
                },
                right: Terminus::CLEAVED,
            },
        ] {
            recipe.fragments[0].source.preparation = prep;
            let unknown = t::preview(&inputs, &target, &recipe).unwrap();
            assert_ne!(
                unknown.verification.status,
                Status::ConsistentUnderAssumptions
            );
            assert!(!unknown.preparation_complete);
            assert!(unknown.product_physical.is_none());
        }
    }
    #[test]
    fn optional_treatment_never_fabricates_phosphate_or_silently_repairs_retained_nicks() {
        let (inputs, target, mut recipe) =
            fixture(TailingBase::T, TailingBase::A, "AAACCCTGGGTTTA");
        recipe.fragments[1].treatment = None;
        let blocked = t::preview(&inputs, &target, &recipe).unwrap();
        assert!(blocked.product_physical.is_none());
        assert_ne!(
            blocked.verification.status,
            Status::ConsistentUnderAssumptions
        );
        recipe.allow_nicked = true;
        let nicked = t::preview(&inputs, &target, &recipe).unwrap();
        let physical = nicked.product_physical.unwrap();
        assert_eq!(physical.nicks.len(), 2);
        assert!(!physical.covalently_closed_in_model);
        assert!(same_circle(
            nicked
                .verification
                .final_product()
                .expect("final declared step must reconstruct")
                .sequence
                .as_str(),
            "AAACCCTGGGTTTA"
        ));
    }
    #[test]
    fn pinned_catalog_and_pcr_catalog_preparations_reconstruct_their_actual_blunt_arms() {
        let catalog = collagene_core::cloning::restriction::catalog::builtin_catalog().unwrap();
        let id = catalog
            .edition()
            .enzymes
            .iter()
            .find(|e| e.name == "SmaI")
            .unwrap()
            .id
            .clone();
        let selection =
            collagene_core::cloning::restriction::catalog::select(&catalog, &[id]).unwrap();
        // Choose internal fragments between two actual SmaI cuts. Outer original
        // linear/PCR boundaries would retain unknown/5′OH chemistry, not fresh cuts.
        let inputs = vec![
            input("vector", "CCCGGGAAACCCGGG"),
            input("insert", "ATATCCCGGGTTTCCCGGGGCGC"),
        ];
        let target = Molecule::dna(
            "Independent pinned arm intent",
            "GGGAAACCCTGGGTTTCCCA",
            Topology::Circular,
        )
        .unwrap();
        let recipe = t::Recipe {
            standard: None,
            profile: collagene_core::cloning::ligation::tailing_recipe::Profile::Ta,

            name: "Pinned source and whole-PCR digest".into(),
            fragments: vec![
                t::Fragment {
                    source: FragmentSpec {
                        generation: None,
                        role: "vector".into(),
                        input_id: "vector".into(),
                        preparation: Preparation::CatalogDigest {
                            conditions: None,
                            selection: selection.clone(),
                            fragment_index: 1,
                            repair: None,
                        },
                        reverse: false,
                    },
                    action: collagene_core::cloning::ligation::tailing_recipe::Action::AddTail {
                        base: TailingBase::T,
                    },
                    treatment: None,
                },
                t::Fragment {
                    source: FragmentSpec {
                        generation: None,
                        role: "insert".into(),
                        input_id: "insert".into(),
                        preparation: Preparation::PcrCatalogDigest {
                            conditions: None,
                            forward: Primer::new("F", "ATATCCC", "").unwrap(),
                            reverse: Primer::new("R", "GCGCCCC", "").unwrap(),
                            selection: selection.clone(),
                            fragment_index: 1,
                            repair: None,
                            template_strand: collagene_core::chemistry::Strand::Top,
                        },
                        reverse: false,
                    },
                    action: collagene_core::cloning::ligation::tailing_recipe::Action::AddTail {
                        base: TailingBase::A,
                    },
                    treatment: None,
                },
            ],
            circular: true,
            allow_nicked: false,
        };
        let p = t::preview(&inputs, &target, &recipe).unwrap();
        assert_eq!(
            p.verification.status,
            Status::ConsistentUnderAssumptions,
            "{:?}",
            p.verification.findings
        );
        assert!(same_circle(
            p.verification
                .final_product()
                .expect("final declared step must reconstruct")
                .sequence
                .as_str(),
            "GGGAAACCCTGGGTTTCCCA"
        ));
        assert!(
            matches!(&p.compiled.plan.steps[0],Operation::CatalogDigest{selection:actual,..}if actual.catalog_content_hash==selection.catalog_content_hash)
        );
        assert!(matches!(&p.compiled.plan.steps[1], Operation::Pcr { .. }));
        assert!(matches!(
            &p.compiled.plan.steps[2],
            Operation::CatalogDigest {
                template: Source::Step(1),
                ..
            }
        ));
        let mut forged = recipe;
        let Preparation::CatalogDigest { selection, .. } =
            &mut forged.fragments[0].source.preparation
        else {
            unreachable!()
        };
        selection.enzymes[0].enzyme.top_cut += 1;
        assert!(t::preflight(&inputs, &target, &forged).is_err());
    }
    #[test]
    fn receipt_validation_rejects_rebound_tail_treatment_preparation_and_source_context() {
        let (inputs, target, recipe) = fixture(TailingBase::T, TailingBase::A, "AAACCCTGGGTTTA");
        let good = t::preview(&inputs, &target, &recipe).unwrap();
        let expected = good.context_revision.clone();
        let rehash = |p: &mut t::Preview| {
            p.verification.plan_revision =
                collagene_core::planning::plan::DesignPlan::revision(&p.compiled.plan)
        };
        let mut bad = good.clone();
        let Operation::Tail3Prime { base, .. } =
            &mut bad.compiled.plan.steps[bad.compiled.rows[0].tail_step.unwrap()]
        else {
            unreachable!()
        };
        *base = TailingBase::G;
        rehash(&mut bad);
        assert!(bad.validate_binding(&expected).is_err());
        let mut bad = good.clone();
        let Operation::TreatEnds { template, .. } =
            &mut bad.compiled.plan.steps[bad.compiled.rows[0].treatment_step.unwrap()]
        else {
            unreachable!()
        };
        *template = Source::Input("insert".into());
        rehash(&mut bad);
        assert!(bad.validate_binding(&expected).is_err());
        let mut bad = good.clone();
        let Operation::Pcr { template, .. } = &mut bad.compiled.plan.steps[0] else {
            unreachable!()
        };
        *template = Source::Input("insert".into());
        rehash(&mut bad);
        assert!(bad.validate_binding(&expected).is_err());
        let mut bad = good.clone();
        bad.compiled.rows[0].tail_step = Some(usize::MAX);
        assert!(bad.validate_binding(&expected).is_err());
        let mut stale = inputs;
        stale[0].molecule.name = "Different exact source".into();
        assert!(t::preview_revision(&stale, &target, &recipe).is_err());
        let mut bad = good.clone();
        bad.prepared[0]
            .as_mut()
            .expect("prepared row must reconstruct")
            .name = "Foreign prepared document".into();
        assert!(bad.validate_binding(&expected).is_err());
        let mut bad = good;
        bad.junctions[0].upstream_end = None;
        assert!(bad.validate_binding(&expected).is_err());
    }
}
