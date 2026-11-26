//! Literal TA/GC circular targets and actual pinned preparation chains; no returned product supplies intent.
use collagene_core::{
    Molecule, Topology,
    analysis::Primer,
    assembly_recipe::{FragmentSpec, Preparation},
    chemistry::{EndTreatment, Polarity, Terminus, ThreePrime},
    tailing::TailingBase,
    tailing_recipe as t,
    verification::{self, Operation, PlanInput, Source, Status},
};
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
            role: format!("{id}-role"),
            input_id: id.into(),
            preparation: Preparation::ExactPcr {
                forward: Primer::new("Forward", forward, "").unwrap(),
                reverse: Primer::new("Reverse", reverse, "").unwrap(),
            },
            reverse: false,
        },
        base,
        treatment: Some(EndTreatment::Phosphorylate),
    }
}
fn fixture(a: TailingBase, b: TailingBase, target: &str) -> (Vec<PlanInput>, Molecule, t::Recipe) {
    (
        vec![input("vector", "AAACCC"), input("insert", "GGGTTT")],
        Molecule::dna("Independent scar intent", target, Topology::Circular).unwrap(),
        t::Recipe {
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
        let actual = p.verification.products.last().unwrap();
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
                matches!(&p.compiled.plan.steps[row.tail_step],Operation::Tail3Prime{base:actual,..}if actual==base)
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
    target.sequence = collagene_core::Sequence::dna("AAACCCTGGGTTTG").unwrap();
    let wrong = t::preview(&inputs, &target, &recipe).unwrap();
    assert_eq!(wrong.verification.status, Status::Invalid);
    assert!(wrong.preparation_complete);
    assert!(wrong.verification.products.last().is_some());
    recipe.fragments[1].base = TailingBase::T;
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
    let (inputs, target, mut recipe) = fixture(TailingBase::T, TailingBase::A, "AAACCCTGGGTTTA");
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
            .products
            .last()
            .unwrap()
            .sequence
            .as_str(),
        "AAACCCTGGGTTTA"
    ));
}
#[test]
fn pinned_catalog_and_pcr_catalog_preparations_reconstruct_their_actual_blunt_arms() {
    let catalog = collagene_core::catalog_digest::builtin_catalog().unwrap();
    let id = catalog
        .edition()
        .enzymes
        .iter()
        .find(|e| e.name == "SmaI")
        .unwrap()
        .id
        .clone();
    let selection = collagene_core::catalog_digest::select(&catalog, &[id]).unwrap();
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
        name: "Pinned source and whole-PCR digest".into(),
        fragments: vec![
            t::Fragment {
                source: FragmentSpec {
                    role: "vector".into(),
                    input_id: "vector".into(),
                    preparation: Preparation::CatalogDigest {
                        selection: selection.clone(),
                        fragment_index: 1,
                        repair: None,
                    },
                    reverse: false,
                },
                base: TailingBase::T,
                treatment: None,
            },
            t::Fragment {
                source: FragmentSpec {
                    role: "insert".into(),
                    input_id: "insert".into(),
                    preparation: Preparation::PcrCatalogDigest {
                        forward: Primer::new("F", "ATATCCC", "").unwrap(),
                        reverse: Primer::new("R", "GCGCCCC", "").unwrap(),
                        selection: selection.clone(),
                        fragment_index: 1,
                        repair: None,
                    },
                    reverse: false,
                },
                base: TailingBase::A,
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
        p.verification.products.last().unwrap().sequence.as_str(),
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
    let Preparation::CatalogDigest { selection, .. } = &mut forged.fragments[0].source.preparation
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
        p.verification.plan_revision = verification::plan_revision(&p.compiled.plan)
    };
    let mut bad = good.clone();
    let Operation::Tail3Prime { base, .. } =
        &mut bad.compiled.plan.steps[bad.compiled.rows[0].tail_step]
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
    bad.compiled.rows[0].tail_step = usize::MAX;
    assert!(bad.validate_binding(&expected).is_err());
    let mut stale = inputs;
    stale[0].molecule.name = "Different exact source".into();
    assert!(t::preview_revision(&stale, &target, &recipe).is_err());
    let mut bad = good.clone();
    bad.prepared[0].name = "Foreign prepared document".into();
    assert!(bad.validate_binding(&expected).is_err());
    let mut bad = good;
    bad.junctions[0].upstream_end = None;
    assert!(bad.validate_binding(&expected).is_err());
}
