use collagene_core::{
    cloning::{
        assembly::{
            intent::{self as assembly_intent, Contribution, Selector},
            junction_edit::*,
            recipe::{self as assembly_recipe, Preparation},
        },
        pcr::primers::DesignOptions,
    },
    planning::plan::PlanInput,
    verification::{self, evidence::Status},
};
use collagene_model::{Feature, Location, Molecule, Topology};
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
        material_bindings: vec![],
        name: "Explicit edited assembly".into(),
        fragments: vec![
            Contribution {
                role: "left".into(),
                input_id: "l".into(),
                selector: Selector::ReferenceRange {
                    start: 0,
                    end: 20,
                    reverse: false,
                },
            },
            Contribution {
                role: "right".into(),
                input_id: "r".into(),
                selector: Selector::ReferenceRange {
                    start: 0,
                    end: 20,
                    reverse: false,
                },
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

#[test]
fn literal_linear_terminal_intent_maps_last_suffix_without_extending_source_annotations() {
    let mut inputs = vec![input("vector", "AAACCC"), input("insert", "GGGTTT")];
    for source in &mut inputs {
        source.molecule.features.push(Feature {
            kind: "CDS".into(),
            location: Location::range(0, 6),
            qualifiers: vec![("label".into(), Some(format!("{} core", source.id)))],
            needs_review: false,
        });
        source.revision = source.molecule.revision();
    }
    let before = serde_json::to_value(&inputs).unwrap();
    let base = assembly_intent::Request {
        material_bindings: vec![],
        name: "Independently declared linear target".into(),
        circular: false,
        fragments: inputs
            .iter()
            .map(|source| Contribution {
                role: source.id.clone(),
                input_id: source.id.clone(),
                selector: Selector::ReferenceRange {
                    start: 0,
                    end: 6,
                    reverse: false,
                },
            })
            .collect(),
    };
    let edits = vec![
        Edit {
            junction: 0,
            remove_left: 0,
            remove_right: 0,
            inserted: "T".into(),
        },
        Edit {
            junction: 1,
            remove_left: 0,
            remove_right: 0,
            inserted: "a".into(),
        },
    ];
    intent_preflight(
        &inputs
            .iter()
            .map(|input| (input.id.as_str(), &input.molecule))
            .collect::<Vec<_>>(),
        &base,
        &edits,
        None,
    )
    .unwrap();
    let intent = build_intent(&inputs, &base, &edits, None).unwrap();
    assert_eq!(intent.target.sequence.as_str(), "AAACCCTGGGTTTA");
    assert_eq!(intent.target.topology, Topology::Linear);
    assert_eq!(
        (
            intent.fragments[0].output_start,
            intent.fragments[0].output_end
        ),
        (0, 7)
    );
    assert_eq!(
        (
            intent.fragments[1].output_start,
            intent.fragments[1].output_end
        ),
        (7, 14)
    );
    for (label, expected) in [("vector core", (0, 6)), ("insert core", (7, 13))] {
        let feature = intent
            .target
            .features
            .iter()
            .find(|feature| feature.label() == label)
            .unwrap();
        assert_eq!(
            feature.location.spans(),
            vec![(expected.0, expected.1, false)]
        );
    }
    assert!(
        intent
            .target
            .features
            .iter()
            .any(|feature| feature.needs_review && feature.location == Location::range(13, 14))
    );
    intent
        .validate_binding(&intent_revision(&inputs, &base, &edits, None).unwrap())
        .unwrap();
    let mut other = edits.clone();
    other[1].inserted = "G".into();
    assert!(
        intent
            .validate_binding(&intent_revision(&inputs, &base, &other, None).unwrap())
            .is_err()
    );
    assert_eq!(serde_json::to_value(&inputs).unwrap(), before);
    assert!(collagene_core::chemistry::Duplex::from_molecule(&intent.target).is_err());
}

#[test]
fn terminal_trim_respects_reverse_origin_spanning_core_without_wrapping_into_first_row() {
    let mut inputs = vec![input("first", "AACCGGTT"), input("last", "AACCGGTT")];
    for source in &mut inputs {
        source.molecule.topology = Topology::Circular;
        source.revision = source.molecule.revision();
    }
    inputs[1].molecule.features.push(Feature {
        kind: "misc_feature".into(),
        location: Location::range(0, 4),
        qualifiers: vec![("label".into(), Some("Last original core".into()))],
        needs_review: false,
    });
    inputs[1].revision = inputs[1].molecule.revision();
    let base = assembly_intent::Request {
        material_bindings: vec![],
        name: "Original digest origins retained".into(),
        circular: false,
        fragments: vec![
            Contribution {
                role: "first".into(),
                input_id: "first".into(),
                selector: Selector::ReferenceRange {
                    start: 6,
                    end: 12,
                    reverse: false,
                },
            },
            Contribution {
                role: "last".into(),
                input_id: "last".into(),
                selector: Selector::ReferenceRange {
                    start: 6,
                    end: 12,
                    reverse: true,
                },
            },
        ],
    };
    let edits = vec![
        Edit {
            junction: 0,
            remove_left: 0,
            remove_right: 0,
            inserted: "G".into(),
        },
        Edit {
            junction: 1,
            remove_left: 2,
            remove_right: 0,
            inserted: "C".into(),
        },
    ];
    let intent = build_intent(&inputs, &base, &edits, None).unwrap();
    // First unrolled core TTAACC, literal G, reverse trimmed core GGTT, literal C.
    assert_eq!(intent.target.sequence.as_str(), "TTAACCGGGTTC");
    assert_eq!(
        (
            intent.adjusted_fragments[0]
                .selector
                .reference_interval()
                .unwrap()
                .0,
            intent.adjusted_fragments[0]
                .selector
                .reference_interval()
                .unwrap()
                .1
        ),
        (6, 12)
    );
    assert_eq!(
        (
            intent.adjusted_fragments[1]
                .selector
                .reference_interval()
                .unwrap()
                .0,
            intent.adjusted_fragments[1]
                .selector
                .reference_interval()
                .unwrap()
                .1
        ),
        (0, 4)
    );
    let feature = intent
        .target
        .features
        .iter()
        .find(|feature| feature.label() == "Last original core")
        .unwrap();
    assert_eq!(feature.location.spans(), vec![(7, 11, true)]);
    let mut downstream = edits.clone();
    downstream[1].remove_right = 1;
    assert!(build_intent(&inputs, &base, &downstream, None).is_err());
    let mut empty_core = edits;
    empty_core[1].remove_left = 6;
    assert!(build_intent(&inputs, &base, &empty_core, None).is_err());
}

#[test]
fn terminal_literal_is_intent_only_and_not_a_silent_overlap_primer_tail() {
    let (inputs, mut request) = fixture();
    request.edits.push(Edit {
        junction: 1,
        remove_left: 0,
        remove_right: 0,
        inserted: "A".into(),
    });
    assert!(build_intent(&inputs, &request.base, &request.edits, None).is_ok());
    assert!(
        preflight(
            &inputs
                .iter()
                .map(|input| (input.id.as_str(), &input.molecule))
                .collect::<Vec<_>>(),
            &request
        )
        .is_err()
    );
    assert!(design_revision(&inputs, &request).is_err());
    assert!(design(&inputs, &request).is_err());
    request.edits[1].remove_right = 1;
    assert!(intent_revision(&inputs, &request.base, &request.edits, None).is_err());
}
fn verify(inputs: &[PlanInput], r: &Report) {
    let compiled = assembly_recipe::compile(inputs, &r.intent.target, &r.recipe).unwrap();
    let report = verification::engine::verify(&compiled.plan, &r.intent.target);
    assert_eq!(
        report.status,
        Status::ConsistentUnderAssumptions,
        "{:?}",
        report.findings
    );
    assert_eq!(
        report
            .final_product()
            .expect("final declared step must reconstruct")
            .sequence,
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
            proposed.intent.adjusted_fragments[0]
                .selector
                .reference_interval()
                .unwrap()
                .0,
            proposed.intent.adjusted_fragments[0]
                .selector
                .reference_interval()
                .unwrap()
                .1
        ),
        (0, 18)
    );
    assert_eq!(
        (
            proposed.intent.adjusted_fragments[1]
                .selector
                .reference_interval()
                .unwrap()
                .0,
            proposed.intent.adjusted_fragments[1]
                .selector
                .reference_interval()
                .unwrap()
                .1
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
    assert_eq!(
        verification::engine::verify(&plan, &wrong).status,
        Status::Invalid
    );
}
#[test]
fn reverse_row_uses_adjusted_original_coordinates_and_original_forward_tail() {
    let (inputs, mut r) = fixture();
    let Selector::ReferenceRange {
        reverse: selected_reverse,
        ..
    } = &mut r.base.fragments[0].selector
    else {
        panic!("reference selector")
    };
    *selected_reverse = true;
    let p = design(&inputs, &r).unwrap();
    assert_eq!(
        p.intent.target.sequence.as_str(),
        "TCTGGAAACCCCGGGGTCGGACCGTTAGCGACCTTGGAA"
    );
    assert_eq!(
        (
            p.intent.adjusted_fragments[0]
                .selector
                .reference_interval()
                .unwrap()
                .0,
            p.intent.adjusted_fragments[0]
                .selector
                .reference_interval()
                .unwrap()
                .1
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
        material_bindings: vec![],
        name: "Trim beyond original origin".into(),
        fragments: vec![Contribution {
            role: "wrap".into(),
            input_id: "wrap".into(),
            selector: Selector::ReferenceRange {
                start: 6,
                end: 12,
                reverse: false,
            },
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
        None,
    )
    .unwrap();
    assert_eq!(intent.target.sequence.as_str(), "ACG");
    assert_eq!(
        (
            intent.adjusted_fragments[0]
                .selector
                .reference_interval()
                .unwrap()
                .0,
            intent.adjusted_fragments[0]
                .selector
                .reference_interval()
                .unwrap()
                .1
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
    let intent = build_intent(&inputs, &r.base, &r.edits, None).unwrap();
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
    let intent = build_intent(&inputs, &r.base, &r.edits, None).unwrap();
    let mut changed = r.clone();
    changed.edits[0].inserted = "GGT".into();
    assert_ne!(
        intent.context_revision,
        intent_revision(&inputs, &changed.base, &changed.edits, None).unwrap()
    );
    assert!(
        intent
            .validate_binding(
                &intent_revision(&inputs, &changed.base, &changed.edits, None).unwrap()
            )
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
            &duplicate.edits,
            None
        )
        .is_err()
    );
    let mut over = r.clone();
    over.edits[0].remove_left = 20;
    assert!(build_intent(&inputs, &over.base, &over.edits, None).is_err());
    let mut ambiguous = r.clone();
    ambiguous.edits[0].inserted = "N".into();
    assert!(build_intent(&inputs, &ambiguous.base, &ambiguous.edits, None).is_ok());
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
    assert!(build_intent(&inputs, &long.base, &long.edits, None).is_err());
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
