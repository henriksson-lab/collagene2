use collagene_core::{
    Molecule, Topology,
    analysis::Primer,
    assembly_recipe::{FragmentSpec, Preparation},
    chemistry::{EndTreatment, FivePrime, ThreePrime},
    ligation_recipe as l,
    verification::{Operation, PlanInput, Source, Status},
};
fn molecule(name: &str, dna: &str) -> Molecule {
    Molecule::dna(name, dna, Topology::Linear).unwrap()
}
fn input(id: &str, m: Molecule) -> PlanInput {
    PlanInput {
        id: id.into(),
        revision: m.revision(),
        molecule: m,
    }
}
fn pcr(role: &str, id: &str, f: &str, r: &str, treatment: Option<EndTreatment>) -> l::Fragment {
    l::Fragment {
        source: FragmentSpec {
            role: role.into(),
            input_id: id.into(),
            preparation: Preparation::ExactPcr {
                forward: Primer::new("actual forward", f, "").unwrap(),
                reverse: Primer::new("actual reverse", r, "").unwrap(),
            },
            reverse: false,
        },
        treatment,
    }
}
fn fixture() -> (Vec<PlanInput>, Molecule, l::Recipe) {
    let inputs = vec![
        input("left", molecule("left original", "TTGCAATCGTCGATGCTAGC")),
        input("right", molecule("right original", "ATCGTCGATGCTAGCGGTAC")),
    ];
    let target = molecule(
        "independently supplied target",
        "TTGCAATCGTCGATGCTAGCATCGTCGATGCTAGCGGTAC",
    );
    let recipe = l::Recipe {
        name: "prepared blunt ligation".into(),
        fragments: vec![
            pcr(
                "first",
                "left",
                "TTGCAATC",
                "GCTAGCAT",
                Some(EndTreatment::Phosphorylate),
            ),
            pcr(
                "second",
                "right",
                "ATCGTCGA",
                "GTACCGCT",
                Some(EndTreatment::Phosphorylate),
            ),
        ],
        circular: false,
        allow_nicked: false,
    };
    (inputs, target, recipe)
}
#[test]
fn actual_pcr_phosphorylation_and_junctions_come_from_fresh_physical_trace() {
    let (inputs, target, recipe) = fixture();
    let p = l::preview(&inputs, &target, &recipe).unwrap();
    assert_eq!(
        p.verification.status,
        Status::ConsistentUnderAssumptions,
        "{:?}",
        p.verification.findings
    );
    assert!(p.preparation_complete);
    assert_eq!(p.prepared.len(), 2);
    assert_eq!(p.junctions.len(), 1);
    assert_eq!(p.junctions[0].upstream, 0);
    assert_eq!(p.junctions[0].downstream, 1);
    let ends = p.physical[0].as_ref().unwrap().ends.as_ref().unwrap();
    assert_eq!(ends.0.chemistry.five_prime, FivePrime::Phosphate);
    assert_eq!(ends.0.chemistry.three_prime, ThreePrime::Hydroxyl);
    assert!(p.product_physical.as_ref().unwrap().nicks.is_empty());
    assert!(matches!(
        p.compiled.plan.steps[p.compiled.rows[0].treatment_step.unwrap()],
        Operation::TreatEnds {
            template: Source::Step(0),
            treatment: EndTreatment::Phosphorylate
        }
    ));
    assert_eq!(
        p.verification.products.last().unwrap().sequence,
        target.sequence
    );
    p.validate_binding(&l::preview_revision(&inputs, &target, &recipe).unwrap())
        .unwrap();
    let mut forged = p.clone();
    forged.junctions[0].upstream_end.as_mut().unwrap().overhang = "A".into();
    assert!(forged.validate().is_err());
}
#[test]
fn nick_acceptance_circular_closure_and_wrong_target_are_explicit() {
    let (inputs, mut target, mut recipe) = fixture();
    recipe.fragments[1].treatment = None;
    let rejected = l::preview(&inputs, &target, &recipe).unwrap();
    assert_eq!(rejected.verification.status, Status::Invalid);
    assert!(rejected.product_physical.is_none());
    recipe.allow_nicked = true;
    let nicked = l::preview(&inputs, &target, &recipe).unwrap();
    assert_eq!(
        nicked.verification.status,
        Status::ConsistentUnderAssumptions,
        "{:?}",
        nicked.verification.findings
    );
    assert_eq!(nicked.product_physical.as_ref().unwrap().nicks.len(), 1);
    recipe.circular = true;
    target.topology = Topology::Circular;
    let circle = l::preview(&inputs, &target, &recipe).unwrap();
    assert_eq!(
        circle.verification.status,
        Status::ConsistentUnderAssumptions,
        "{:?}",
        circle.verification.findings
    );
    assert_eq!(circle.product_physical.as_ref().unwrap().nicks.len(), 2);
    assert!(circle.junctions[1].closing);
    assert_eq!(
        (circle.junctions[1].upstream, circle.junctions[1].downstream),
        (1, 0)
    );
    let wrong = Molecule::dna(
        "wrong independenttarget",
        "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA",
        Topology::Circular,
    )
    .unwrap();
    let mismatch = l::preview(&inputs, &wrong, &recipe).unwrap();
    assert_eq!(mismatch.verification.status, Status::Invalid);
    assert!(mismatch.product_physical.is_some());
    assert!(mismatch.preparation_complete);
}
#[test]
fn orientation_and_unphosphorylated_products_cannot_be_auto_repaired() {
    let (inputs, target, mut recipe) = fixture();
    recipe.fragments[1].source.reverse = true;
    let reverse = l::preview(&inputs, &target, &recipe).unwrap();
    assert_eq!(reverse.verification.status, Status::Invalid);
    assert!(reverse.product_physical.is_some());
    for f in &mut recipe.fragments {
        f.treatment = None;
        f.source.reverse = false;
    }
    recipe.allow_nicked = true;
    let unsealed = l::preview(&inputs, &target, &recipe).unwrap();
    assert_eq!(unsealed.verification.status, Status::Invalid);
    assert!(unsealed.product_physical.is_none());
    assert!(
        unsealed
            .verification
            .findings
            .iter()
            .any(|f| f.contains("neither strand"))
    );
}
#[test]
fn plain_reference_unknown_ends_and_mixed_strands_are_not_physical_authority() {
    let (inputs, target, mut recipe) = fixture();
    for f in &mut recipe.fragments {
        f.source.preparation = Preparation::WholeLinear;
    }
    let p = l::preview(&inputs, &target, &recipe).unwrap();
    assert_eq!(
        p.verification.status,
        Status::Incomplete,
        "{:?}",
        p.verification.findings
    );
    let summary = p.physical[0].as_ref().unwrap();
    assert_eq!(
        summary.ends.as_ref().unwrap().0.chemistry.three_prime,
        ThreePrime::Unknown
    );
    assert!(p.product_physical.is_none());
    let mut mixed = inputs.clone();
    mixed[0]
        .molecule
        .metadata
        .insert("genbank.molecule_type".into(), "ss-DNA".into());
    mixed[0].revision = mixed[0].molecule.revision();
    recipe.fragments[0].treatment = None;
    let p = l::preview(&mixed, &target, &recipe).unwrap();
    assert_eq!(p.verification.status, Status::Incomplete);
    assert!(p.physical[0].is_none());
}

#[test]
fn three_fresh_cohesive_digest_arms_close_but_opposite_polarities_do_not() {
    let mut inputs = vec![];
    let mut fragments = vec![];
    for (index, dna) in ["AAAAGAATTCTTTT", "CCCCGAATTCGGGG", "ACACGAATTCTGTG"]
        .iter()
        .enumerate()
    {
        let id = format!("circle-{index}");
        let mut m = molecule(&id, dna);
        m.topology = Topology::Circular;
        inputs.push(input(&id, m));
        fragments.push(l::Fragment {
            source: FragmentSpec {
                role: format!("arm-{index}"),
                input_id: id,
                preparation: Preparation::Digest {
                    enzymes: vec!["EcoRI".into()],
                    fragment_index: 0,
                    repair: None,
                },
                reverse: false,
            },
            treatment: None,
        });
    }
    let target = Molecule::dna(
        "independentthreewaytarget",
        "AATTCTTTTAAAAGAATTCGGGGCCCCGAATTCTGTGACACG",
        Topology::Circular,
    )
    .unwrap();
    let mut recipe = l::Recipe {
        name: "actualthreewaycohesiveassembly".into(),
        fragments,
        circular: true,
        allow_nicked: false,
    };
    let joined = l::preview(&inputs, &target, &recipe).unwrap();
    assert_eq!(
        joined.verification.status,
        Status::ConsistentUnderAssumptions,
        "{:?}",
        joined.verification.findings
    );
    assert_eq!(joined.junctions.len(), 3);
    assert!(
        joined
            .product_physical
            .as_ref()
            .unwrap()
            .covalently_closed_in_model
    );
    assert!(
        joined
            .junctions
            .iter()
            .all(|j| j.upstream_end.as_ref().unwrap().overhang == "AATT")
    );
    inputs[1].molecule = Molecule::dna(
        "oppositepolaritysource",
        "CCCCCTGCAGGGGG",
        Topology::Circular,
    )
    .unwrap();
    inputs[1].revision = inputs[1].molecule.revision();
    if let Preparation::Digest { enzymes, .. } = &mut recipe.fragments[1].source.preparation {
        *enzymes = vec!["PstI".into()];
    }
    let mismatch = l::preview(&inputs, &target, &recipe).unwrap();
    assert!(mismatch.preparation_complete);
    assert_eq!(mismatch.verification.status, Status::Invalid);
    assert!(mismatch.product_physical.is_none());
    assert!(
        mismatch
            .verification
            .findings
            .iter()
            .any(|f| f.contains("polarity"))
    );
}
