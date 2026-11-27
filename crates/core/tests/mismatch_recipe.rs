//! Independent literal mutation intent through the ordinary shared preparation compiler.
use collagene_core::{
    Molecule, Topology,
    analysis::Primer,
    assembly_recipe::{self, FragmentSpec, Preparation, RunPreparation},
    chemistry::{EndTreatment, FivePrime},
    cloning_outputs::{self, Options},
    ligation_recipe as lig,
    mutagenesis::PcrOptions,
    run_lineage::{EventBinding, InputEventBinding},
    run_material::{self, Selection},
    verification::{self, DesignPlan, Operation, PlanInput, Source, Status},
    workspace::Workspace,
};
const ORIGINAL: &str = "AGTCCGATGCTACGTAGGACTTCAGACCTGTA";
const INTENDED: &str = "AGACCGATGCTACGTAGGACTTCAGACCTGGA";
fn policy() -> PcrOptions {
    let mut o = PcrOptions::default();
    o.binding.max_mismatches = 1;
    o.binding.min_three_prime_exact = 2;
    o
}
fn primers() -> (Primer, Primer) {
    (
        Primer::new("Forward intentional change", "AGACCGAT", "").unwrap(),
        Primer::new("Reverse intentional change", "TCCAGGTC", "").unwrap(),
    )
}
fn input(molecule: Molecule) -> PlanInput {
    PlanInput {
        id: "source".into(),
        revision: molecule.revision(),
        molecule,
    }
}
fn fragment(preparation: Preparation) -> FragmentSpec {
    FragmentSpec {
        role: "mutated-fragment".into(),
        input_id: "source".into(),
        preparation,
        reverse: false,
    }
}
fn ordinary() -> Preparation {
    let (forward, reverse) = primers();
    Preparation::MismatchPcr {
        forward,
        reverse,
        options: policy(),
    }
}
fn recipe(preparation: Preparation) -> lig::Recipe {
    lig::Recipe {
        name: "Independent mismatch-PCR construct".into(),
        circular: true,
        allow_nicked: false,
        fragments: vec![lig::Fragment {
            source: fragment(preparation),
            treatment: None,
        }],
    }
}
fn binding(w: &Workspace, index: usize) -> EventBinding {
    let d = &w.documents[index];
    EventBinding {
        document_id: d.id.clone(),
        event_id: d.head_event.clone().unwrap(),
        revision: d.head.clone(),
    }
}
#[test]
fn primer_encoded_changes_follow_literal_target_and_real_pcr_chemistry() {
    let inputs = vec![input(
        Molecule::dna("Original32bp", ORIGINAL, Topology::Linear).unwrap(),
    )];
    let target = Molecule::dna(
        "Independent intended32bp circle",
        INTENDED,
        Topology::Circular,
    )
    .unwrap();
    let mut r = recipe(ordinary());
    let untreated = lig::preview(&inputs, &target, &r).unwrap();
    assert!(untreated.preparation_complete);
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
    r.fragments[0].treatment = Some(EndTreatment::Phosphorylate);
    let p = lig::preview(&inputs, &target, &r).unwrap();
    assert_eq!(
        p.verification.status,
        Status::ConsistentUnderAssumptions,
        "{:?}",
        p.verification.findings
    );
    assert_eq!(p.verification.products[0].sequence.as_str(), INTENDED);
    assert!(
        p.product_physical
            .as_ref()
            .unwrap()
            .covalently_closed_in_model
    );
    assert!(matches!(
        p.compiled.plan.steps[0],
        Operation::PcrWithMismatches {
            template: Source::Input(_),
            ..
        }
    ));
    let d = collagene_core::assembly_diagnostics::ligation(&p).unwrap();
    assert_eq!(d.rows[0].total_steps, 2);
    assert!(d.failure.is_none());
    let wrong = Molecule::dna(
        "Original letters are not the declared mutant",
        ORIGINAL,
        Topology::Circular,
    )
    .unwrap();
    let wrong = lig::preview(&inputs, &wrong, &r).unwrap();
    assert_eq!(wrong.verification.status, Status::Invalid);
    assert!(wrong.product_physical.is_some());
}
#[test]
fn mismatch_three_prime_and_search_budgets_remain_explicit_and_cannot_inherit_a_preview() {
    let inputs = vec![input(
        Molecule::dna("Original32bp", ORIGINAL, Topology::Linear).unwrap(),
    )];
    let target = Molecule::dna("Independent mutant", INTENDED, Topology::Circular).unwrap();
    let mut r = recipe(ordinary());
    r.fragments[0].treatment = Some(EndTreatment::Phosphorylate);
    let original = lig::preview(&inputs, &target, &r).unwrap();
    for case in 0..3 {
        let mut changed = r.clone();
        if let Preparation::MismatchPcr { options, .. } =
            &mut changed.fragments[0].source.preparation
        {
            match case {
                0 => options.binding.max_mismatches = 0,
                1 => options.binding.min_three_prime_exact = 6,
                _ => options.binding.max_work = 1,
            }
        }
        let context = lig::preview_revision(&inputs, &target, &changed).unwrap();
        assert_ne!(context, original.context_revision);
        assert!(original.validate_binding(&context).is_err());
        let checked = lig::preview(&inputs, &target, &changed).unwrap();
        assert_ne!(
            checked.verification.status,
            Status::ConsistentUnderAssumptions
        );
        assert!(checked.product_physical.is_none());
    }
    if let Preparation::MismatchPcr { options, .. } = &mut r.fragments[0].source.preparation {
        options.max_products = usize::MAX;
    }
    assert!(lig::source_preflight(&[("source", &inputs[0].molecule)], &target, &r).is_err());
}
#[test]
fn mismatch_preparation_reuses_a_fresh_saved_origin_and_preserves_the_replay_dependency() {
    let source = Molecule::dna("Original32bp", ORIGINAL, Topology::Linear).unwrap();
    let origin_target =
        Molecule::dna("Independent original PCR", ORIGINAL, Topology::Linear).unwrap();
    let plan = DesignPlan {
        name: "Origin exactPCR".into(),
        target_revision: origin_target.revision(),
        inputs: vec![input(source.clone())],
        steps: vec![Operation::Pcr {
            template: Source::Input("source".into()),
            forward: Primer::new("OriginF", "AGTCCGAT", "").unwrap(),
            reverse: Primer::new("OriginR", "TACAGGTC", "").unwrap(),
        }],
    };
    let w = Workspace::from_documents("Original inputs", vec![source, origin_target], 1).unwrap();
    let commit = cloning_outputs::finalize(
        &w,
        &plan,
        binding(&w, 1),
        vec![InputEventBinding {
            input_id: "source".into(),
            event: binding(&w, 0),
        }],
        &Options {
            name: "Saved exactPCR template".into(),
            show_prepared: false,
            save_primers: false,
            allow_incomplete: false,
        },
    )
    .unwrap();
    let event = commit.receipt.construct;
    let origin = run_material::capsule(
        &commit.workspace,
        &Selection {
            run_id: commit.receipt.run.run_id,
            step_index: 0,
            event: event.clone(),
        },
    )
    .unwrap();
    let saved = commit.workspace.revisions[&event.revision].clone();
    let (forward, reverse) = primers();
    let mut r = recipe(Preparation::FromRun {
        origin,
        preparation: RunPreparation::MismatchPcr {
            forward,
            reverse,
            options: policy(),
        },
    });
    r.fragments[0].treatment = Some(EndTreatment::Phosphorylate);
    let inputs = vec![input(saved)];
    let target = Molecule::dna("Independent mutant circle", INTENDED, Topology::Circular).unwrap();
    let p = lig::preview(&inputs, &target, &r).unwrap();
    assert_eq!(
        p.verification.status,
        Status::ConsistentUnderAssumptions,
        "{:?}",
        p.verification.findings
    );
    assert!(matches!(
        p.compiled.plan.steps[0],
        Operation::ReconstructRunProduct { .. }
    ));
    assert!(matches!(
        p.compiled.plan.steps[1],
        Operation::PcrWithMismatches {
            template: Source::Step(0),
            ..
        }
    ));
    assert_eq!(p.compiled.rows[0].preparation.preparation_step, Some(0));
    collagene_core::assembly_diagnostics::ligation(&p).unwrap();
    let mut forged = p.clone();
    if let Operation::PcrWithMismatches { template, .. } = &mut forged.compiled.plan.steps[1] {
        *template = Source::Input("source".into());
    }
    forged.verification.plan_revision = verification::plan_revision(&forged.compiled.plan);
    assert!(forged.validate().is_err());
}
#[test]
fn shared_overlap_compilation_keeps_mismatch_policy_without_inferring_target() {
    let inputs = vec![input(
        Molecule::dna("Original32bp", ORIGINAL, Topology::Linear).unwrap(),
    )];
    let target = Molecule::dna("Independent linear mutant", INTENDED, Topology::Linear).unwrap();
    let r = assembly_recipe::Recipe {
        name: "Explicit overlap recipe with mutagenic preparation".into(),
        fragments: vec![fragment(ordinary())],
        overlaps: vec![],
        circular: false,
        options: Default::default(),
    };
    let compiled = assembly_recipe::compile(&inputs, &target, &r).unwrap();
    assert_eq!(compiled.plan.target_revision, target.revision());
    assert!(matches!(
        compiled.plan.steps[0],
        Operation::PcrWithMismatches { .. }
    ));
    let p = assembly_recipe::preview(&inputs, &target, &r).unwrap();
    assert_eq!(p.verification.products[0].sequence.as_str(), INTENDED);
    p.validate_binding(&assembly_recipe::preview_revision(&inputs, &target, &r).unwrap())
        .unwrap();
}

#[test]
fn tailing_recipe_keeps_literal_mismatch_policy_and_checks_it_readonly() {
    use collagene_core::{tailing::TailingBase, tailing_recipe as tail};
    let inputs = vec![
        input(Molecule::dna("Original32bp", ORIGINAL, Topology::Linear).unwrap()),
        PlanInput {
            id: "insert".into(),
            molecule: Molecule::dna("Literal6bp insert", "GGGTTT", Topology::Linear).unwrap(),
            revision: Molecule::dna("Literal6bp insert", "GGGTTT", Topology::Linear)
                .unwrap()
                .revision(),
        },
    ];
    let target = Molecule::dna(
        "Independent40bp mutant circle",
        "AGACCGATGCTACGTAGGACTTCAGACCTGGATGGGTTTA",
        Topology::Circular,
    )
    .unwrap();
    let recipe = tail::Recipe {
        profile: collagene_core::tailing_recipe::Profile::Ta,

        name: "Mismatch preparation and explicit complementary tails".into(),
        circular: true,
        allow_nicked: false,
        fragments: vec![
            tail::Fragment {
                source: fragment(ordinary()),
                action: collagene_core::tailing_recipe::Action::AddTail {
                    base: TailingBase::T,
                },
                treatment: Some(EndTreatment::Phosphorylate),
            },
            tail::Fragment {
                source: FragmentSpec {
                    role: "insert-role".into(),
                    input_id: "insert".into(),
                    reverse: false,
                    preparation: Preparation::ExactPcr {
                        forward: Primer::new("IF", "GGG", "").unwrap(),
                        reverse: Primer::new("IR", "AAA", "").unwrap(),
                    },
                },
                action: collagene_core::tailing_recipe::Action::AddTail {
                    base: TailingBase::A,
                },
                treatment: Some(EndTreatment::Phosphorylate),
            },
        ],
    };
    let mut p = tail::preview(&inputs, &target, &recipe).unwrap();
    assert_eq!(
        p.verification.status,
        Status::ConsistentUnderAssumptions,
        "{:?}",
        p.verification.findings
    );
    assert_eq!(
        p.verification.products.last().unwrap().sequence.as_str(),
        "AGACCGATGCTACGTAGGACTTCAGACCTGGATGGGTTTA"
    );
    p.validate_binding(&tail::preview_revision(&inputs, &target, &recipe).unwrap())
        .unwrap();
    if let Operation::PcrWithMismatches { options, .. } = &mut p.compiled.plan.steps[0] {
        options.binding.max_mismatches = 0;
    }
    p.verification.plan_revision = verification::plan_revision(&p.compiled.plan);
    assert!(p.validate().is_err());
}
