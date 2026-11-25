use collagene_core::{
    Feature, Location, Molecule, Topology, assembly_intent,
    assembly_junction_edit::*,
    assembly_recipe::{self, Preparation},
    overlap_primer_design::Fragment,
    primers::DesignOptions,
    verification::{self, PlanInput, Status},
};
const L: &str = "ATGACCCCGGGGTTTCCAGA";
const R: &str = "ATCCGTTAGCGACCTTGGAA";
fn input(id: &str, dna: &str) -> PlanInput {
    let molecule = Molecule::dna(id, dna, Topology::Linear).unwrap();
    PlanInput {
        id: id.into(),
        revision: molecule.revision(),
        molecule,
    }
}
fn fixture() -> (Vec<PlanInput>, Request) {
    let inputs = vec![input("l", L), input("r", R)];
    let base = assembly_intent::Request {
        name: "Explicit edited assembly".into(),
        fragments: vec![
            Fragment {
                role: "left".into(),
                input_id: "l".into(),
                start: 0,
                end: 20,
                reverse: false,
            },
            Fragment {
                role: "right".into(),
                input_id: "r".into(),
                start: 0,
                end: 20,
                reverse: false,
            },
        ],
        circular: false,
    };
    let request = Request {
        base,
        edits: vec![Edit {
            junction: 0,
            remove_left: 2,
            remove_right: 2,
            inserted: "GGA".into(),
        }],
        overlap: 8,
        primer_options: DesignOptions {
            min_length: 8,
            max_length: 8,
            min_gc_fraction: 0.,
            max_gc_fraction: 1.,
            max_tm_difference_c: 100.,
            max_homopolymer: 200,
            max_three_prime_complementarity: 200,
            ..Default::default()
        },
        profile_options: Default::default(),
    };
    (inputs, request)
}
fn verify(inputs: &[PlanInput], r: &Report) {
    let compiled = assembly_recipe::compile(inputs, &r.intent.target, &r.recipe).unwrap();
    let report = verification::verify(&compiled.plan, &r.intent.target);
    assert_eq!(
        report.status,
        Status::ConsistentUnderAssumptions,
        "{:?}",
        report.findings
    );
    assert_eq!(
        report.products.last().unwrap().sequence,
        r.intent.target.sequence
    );
}
#[test]
fn literal_linear_trim_and_insert_updates_independent_intent_and_full_tail_metrics() {
    let (inputs, r) = fixture();
    let before = serde_json::to_value(&inputs).unwrap();
    let proposed = design(&inputs, &r).unwrap();
    assert_eq!(
        proposed.intent.target.sequence.as_str(),
        "ATGACCCCGGGGTTTCCAGGACCGTTAGCGACCTTGGAA"
    );
    assert_eq!(
        (
            proposed.intent.adjusted_request.fragments[0].start,
            proposed.intent.adjusted_request.fragments[0].end
        ),
        (0, 18)
    );
    assert_eq!(
        (
            proposed.intent.adjusted_request.fragments[1].start,
            proposed.intent.adjusted_request.fragments[1].end
        ),
        (2, 20)
    );
    assert_eq!(
        (
            proposed.intent.fragments[0].output_start,
            proposed.intent.fragments[0].output_end
        ),
        (0, 21)
    );
    assert_eq!(
        (
            proposed.intent.fragments[1].output_start,
            proposed.intent.fragments[1].output_end
        ),
        (21, 39)
    );
    let pair = &proposed.designs[0].pairs[0];
    assert_eq!(pair.reverse.primer.tail, "GCTAACGGTCC");
    assert_eq!(pair.reverse.primer.annealing.as_str(), "TGGAAACC");
    assert_eq!(
        pair.reverse.analysis.ordered_sequence,
        "GCTAACGGTCCTGGAAACC"
    );
    assert_eq!(
        proposed.designs[0].options.reverse_tail,
        pair.reverse.primer.tail
    );
    assert!(proposed.designs[1].pairs[0].reverse.primer.tail.is_empty());
    assert_eq!(serde_json::to_value(&inputs).unwrap(), before);
    proposed
        .validate_binding(&design_revision(&inputs, &r).unwrap())
        .unwrap();
    verify(&inputs, &proposed);
    let wrong = Molecule::dna(
        "independent wrong target",
        "ATGACCCCGGGGTTTCCAGGACCGTTAGCGACCTTGGAT",
        Topology::Linear,
    )
    .unwrap();
    let plan = assembly_recipe::compile(&inputs, &wrong, &proposed.recipe)
        .unwrap()
        .plan;
    assert_eq!(verification::verify(&plan, &wrong).status, Status::Invalid);
}
#[test]
fn reverse_row_uses_adjusted_original_coordinates_and_original_forward_tail() {
    let (inputs, mut r) = fixture();
    r.base.fragments[0].reverse = true;
    let p = design(&inputs, &r).unwrap();
    assert_eq!(
        p.intent.target.sequence.as_str(),
        "TCTGGAAACCCCGGGGTCGGACCGTTAGCGACCTTGGAA"
    );
    assert_eq!(
        (
            p.intent.adjusted_request.fragments[0].start,
            p.intent.adjusted_request.fragments[0].end
        ),
        (2, 20)
    );
    let pair = &p.designs[0].pairs[0];
    assert_eq!(pair.forward.primer.tail, "GCTAACGGTCC");
    assert!(pair.reverse.primer.tail.is_empty());
    assert!(p.recipe.fragments[0].reverse);
    verify(&inputs, &p);
}
#[test]
fn circle_closure_suffix_preserves_the_first_adjusted_core_origin() {
    let (inputs, mut r) = fixture();
    r.base.circular = true;
    r.edits.push(Edit {
        junction: 1,
        remove_left: 1,
        remove_right: 1,
        inserted: "TTA".into(),
    });
    let p = design(&inputs, &r).unwrap();
    assert_eq!(
        p.intent.target.sequence.as_str(),
        "TGACCCCGGGGTTTCCAGGACCGTTAGCGACCTTGGATTA"
    );
    assert_eq!(p.intent.target.topology, Topology::Circular);
    assert_eq!(p.intent.fragments[0].output_start, 0);
    assert_eq!(p.intent.fragments.last().unwrap().output_end, 40);
    assert_eq!(p.designs[1].pairs[0].reverse.primer.tail, "CGGGGTCATAA");
    verify(&inputs, &p);
    let mut original = input("wrap", "AACCGGTT");
    original.molecule.topology = Topology::Circular;
    original.revision = original.molecule.revision();
    let base = assembly_intent::Request {
        name: "Trim beyond original origin".into(),
        fragments: vec![Fragment {
            role: "wrap".into(),
            input_id: "wrap".into(),
            start: 6,
            end: 12,
            reverse: false,
        }],
        circular: true,
    };
    let intent = build_intent(
        &[original],
        &base,
        &[Edit {
            junction: 0,
            remove_left: 1,
            remove_right: 3,
            inserted: "G".into(),
        }],
    )
    .unwrap();
    assert_eq!(intent.target.sequence.as_str(), "ACG");
    assert_eq!(
        (
            intent.adjusted_request.fragments[0].start,
            intent.adjusted_request.fragments[0].end
        ),
        (1, 3)
    );
}
#[test]
fn trimmed_annotations_are_reviewed_and_inserted_calls_are_explicitly_synthetic() {
    let (mut inputs, r) = fixture();
    inputs[0].molecule.features.push(Feature {
        kind: "CDS".into(),
        location: Location::range(4, 20),
        qualifiers: vec![
            ("label".into(), Some("source gene".into())),
            (
                "translation".into(),
                Some("old annotation assertion".into()),
            ),
        ],
        needs_review: false,
    });
    inputs[0].revision = inputs[0].molecule.revision();
    let before = serde_json::to_value(&inputs).unwrap();
    let intent = build_intent(&inputs, &r.base, &r.edits).unwrap();
    assert!(
        intent
            .target
            .features
            .iter()
            .any(|f| f.kind == "CDS" && f.needs_review)
    );
    let synthetic = intent
        .target
        .features
        .iter()
        .find(|f| {
            f.qualifiers
                .iter()
                .any(|(k, v)| k == "label" && v.as_deref() == Some("Declared junction 1 insert"))
        })
        .unwrap();
    assert_eq!(synthetic.location, Location::range(18, 21));
    assert!(synthetic.needs_review);
    assert_eq!(serde_json::to_value(&inputs).unwrap(), before);
    assert!(collagene_core::chemistry::Duplex::from_molecule(&intent.target).is_err());
}
#[test]
fn stale_context_duplicate_overtrim_ambiguity_and_full_tail_work_profiles_fail_honestly() {
    let (inputs, r) = fixture();
    let intent = build_intent(&inputs, &r.base, &r.edits).unwrap();
    let mut changed = r.clone();
    changed.edits[0].inserted = "GGT".into();
    assert_ne!(
        intent.context_revision,
        intent_revision(&inputs, &changed.base, &changed.edits).unwrap()
    );
    assert!(
        intent
            .validate_binding(&intent_revision(&inputs, &changed.base, &changed.edits).unwrap())
            .is_err()
    );
    let mut stale = inputs.clone();
    stale[0].revision = "0".repeat(64);
    assert!(design(&stale, &r).is_err());
    let mut duplicate = r.clone();
    duplicate.edits.push(duplicate.edits[0].clone());
    assert!(
        intent_preflight(
            &inputs
                .iter()
                .map(|i| (i.id.as_str(), &i.molecule))
                .collect::<Vec<_>>(),
            &duplicate.base,
            &duplicate.edits
        )
        .is_err()
    );
    let mut over = r.clone();
    over.edits[0].remove_left = 20;
    assert!(build_intent(&inputs, &over.base, &over.edits).is_err());
    let mut ambiguous = r.clone();
    ambiguous.edits[0].inserted = "N".into();
    assert!(build_intent(&inputs, &ambiguous.base, &ambiguous.edits).is_ok());
    assert!(design(&inputs, &ambiguous).is_err());
    let mut wide = r.clone();
    wide.primer_options.max_length = 71;
    assert!(
        preflight(
            &inputs
                .iter()
                .map(|i| (i.id.as_str(), &i.molecule))
                .collect::<Vec<_>>(),
            &wide
        )
        .is_err()
    );
    let mut long = r.clone();
    long.edits[0].inserted = "A".repeat(201);
    assert!(build_intent(&inputs, &long.base, &long.edits).is_err());
    let mut tails = r;
    tails.primer_options.forward_tail = "AC".into();
    assert!(design(&inputs, &tails).is_err());
    let p = design(&inputs, &fixture().1).unwrap();
    let Preparation::ExactPcr { .. } = &p.recipe.fragments[0].preparation else {
        panic!("actual original-source PCR required")
    };
    let mut forged = p.clone();
    forged.designs[0].pairs[0].reverse.analysis.ordered_sequence = "ACGT".into();
    assert!(forged.validate_binding(&p.context_revision).is_err());
}
