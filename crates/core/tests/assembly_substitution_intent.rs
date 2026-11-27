//! Independent literal source-oriented substitutions, not a target reconstructed from primers.
use collagene_core::{
    Feature, Location, Molecule, Topology,
    analysis::Primer,
    assembly_intent,
    assembly_intent::{Contribution, ContributionPosition, Selector},
    assembly_junction_edit::Edit,
    assembly_substitution_intent as s,
    verification::{self, DesignPlan, Operation, PlanInput, Source, Status},
};
fn input(molecule: Molecule) -> PlanInput {
    PlanInput {
        id: "original".into(),
        revision: molecule.revision(),
        molecule,
    }
}
fn fragment(role: &str, start: usize, end: usize, reverse: bool) -> Contribution {
    Contribution {
        role: role.into(),
        input_id: "original".into(),
        selector: Selector::ReferenceRange {
            start,
            end,
            reverse,
        },
    }
}
fn edit(role: &str, offset: usize, expected: &str, replacement: &str) -> s::Substitution {
    s::Substitution {
        role: role.into(),
        offset,
        expected: expected.into(),
        replacement: replacement.into(),
    }
}
fn request(fragments: Vec<Contribution>, substitutions: Vec<s::Substitution>) -> s::Request {
    s::Request {
        base: assembly_intent::Request {
            material_bindings: vec![],
            name: "Independently declared mutant".into(),
            fragments,
            circular: false,
        },
        junction_edits: vec![],
        substitutions,
    }
}
#[test]
fn repeated_originals_have_distinct_roles_and_reverse_edits_use_construct_letters() {
    let inputs = vec![input(
        Molecule::dna("Original palindrome", "AACCGGTT", Topology::Linear).unwrap(),
    )];
    let before = serde_json::to_vec(&inputs).unwrap();
    let request = request(
        vec![
            fragment("forward", 0, 8, false),
            fragment("reverse", 0, 8, true),
        ],
        vec![edit("forward", 2, "C", "A"), edit("reverse", 2, "C", "T")],
    );
    let result = s::build(&inputs, &request, None).unwrap();
    assert_eq!(result.target.sequence.as_str(), "AAACGGTTAATCGGTT");
    assert_eq!(
        result
            .applied
            .iter()
            .map(|a| {
                let (position, reverse) = reference_position(&a.source_position);
                (position, a.output_position, reverse)
            })
            .collect::<Vec<_>>(),
        vec![(2, 2, false), (5, 10, true)]
    );
    assert_eq!(serde_json::to_vec(&inputs).unwrap(), before);
    result
        .validate_binding(&s::intent_revision(&inputs, &request, None).unwrap())
        .unwrap();
    let mut only_forward = request.clone();
    only_forward.substitutions.pop();
    assert_eq!(
        s::build(&inputs, &only_forward, None)
            .unwrap()
            .target
            .sequence
            .as_str(),
        "AAACGGTTAACCGGTT"
    );
}
#[test]
fn circular_original_contributions_cross_origin_without_moving_declared_output_origin() {
    let inputs = vec![input(
        Molecule::dna("Original6bp circle", "AAGTCC", Topology::Circular).unwrap(),
    )];
    let forward = request(
        vec![fragment("wrap", 4, 8, false)],
        vec![edit("wrap", 2, "A", "G")],
    );
    let result = s::build(&inputs, &forward, None).unwrap();
    assert_eq!(result.target.sequence.as_str(), "CCGA");
    assert_eq!(
        (
            reference_position(&result.applied[0].source_position).0,
            result.applied[0].output_position
        ),
        (0, 2)
    );
    assert_eq!(
        result.fragments[0].selector.reference_interval().unwrap(),
        (4, 8, false)
    );
    assert_eq!(
        (
            result.fragments[0].retained_start,
            result.fragments[0].retained_end
        ),
        (0, 4)
    );
    let reverse = request(
        vec![fragment("wrap", 4, 8, true)],
        vec![edit("wrap", 0, "T", "A")],
    );
    let result = s::build(&inputs, &reverse, None).unwrap();
    assert_eq!(result.target.sequence.as_str(), "ATGG");
    assert_eq!(reference_position(&result.applied[0].source_position).0, 1);
    assert!(reference_position(&result.applied[0].source_position).1);
}
#[test]
fn retained_cores_and_literal_junctions_keep_separate_offsets_and_trimmed_edits_fail() {
    let inputs = vec![input(
        Molecule::dna("Original8bp", "AACCGGTT", Topology::Linear).unwrap(),
    )];
    let mut request = request(
        vec![
            fragment("left", 0, 4, false),
            fragment("right", 4, 8, false),
        ],
        vec![edit("left", 1, "A", "T"), edit("right", 2, "T", "C")],
    );
    request.junction_edits = vec![Edit {
        junction: 0,
        remove_left: 1,
        remove_right: 1,
        inserted: "AAA".into(),
    }];
    let result = s::build(&inputs, &request, None).unwrap();
    assert_eq!(result.target.sequence.as_str(), "ATCAAAGCT");
    assert_eq!(
        result
            .fragments
            .iter()
            .map(|f| (
                f.retained_start,
                f.retained_end,
                f.output_start,
                f.output_end
            ))
            .collect::<Vec<_>>(),
        vec![(0, 3, 0, 6), (1, 4, 6, 9)]
    );
    assert_eq!(result.applied[1].output_position, 7);
    for substitution in [edit("left", 3, "C", "A"), edit("right", 0, "G", "A")] {
        let mut changed = request.clone();
        changed.substitutions = vec![substitution];
        assert!(s::intent_revision(&inputs, &changed, None).is_err());
    }
    // A final linear suffix remains literal synthetic DNA, never a source offset.
    request.junction_edits.push(Edit {
        junction: 1,
        remove_left: 0,
        remove_right: 0,
        inserted: "GG".into(),
    });
    assert_eq!(
        s::build(&inputs, &request, None)
            .unwrap()
            .target
            .sequence
            .as_str(),
        "ATCAAAGCTGG"
    );
}
#[test]
fn affected_annotations_require_review_and_keep_full_original_qualifiers() {
    let mut source = Molecule::dna("Annotated original", "AACCGGTT", Topology::Linear).unwrap();
    source.features = vec![
        Feature {
            kind: "CDS".into(),
            location: Location::range(0, 8),
            qualifiers: vec![
                ("label".into(), Some("Original coding feature".into())),
                ("translation".into(), Some("FAKE".into())),
                ("note".into(), Some("Original source note".into())),
            ],
            needs_review: false,
        },
        Feature {
            kind: "misc_feature".into(),
            location: Location::Join(vec![Location::range(0, 2), Location::range(6, 8)]),
            qualifiers: vec![("label".into(), Some("Untouched compound".into()))],
            needs_review: false,
        },
    ];
    let inputs = vec![input(source.clone())];
    let r = request(
        vec![fragment("whole", 0, 8, false)],
        vec![edit("whole", 2, "C", "A")],
    );
    let result = s::build(&inputs, &r, None).unwrap();
    assert!(result.target.features[0].needs_review);
    assert!(
        !result.target.features[0]
            .qualifiers
            .iter()
            .any(|(key, _)| key.eq_ignore_ascii_case("translation"))
    );
    assert!(
        result.target.features[0]
            .qualifiers
            .iter()
            .any(|(key, value)| key == "note" && value.as_deref() == Some("Original source note"))
    );
    assert_eq!(result.target.features[1], source.features[1]);
    assert!(
        source.features[0]
            .qualifiers
            .iter()
            .any(|(key, value)| key == "translation" && value.as_deref() == Some("FAKE"))
    );
    assert_eq!(inputs[0].molecule, source);
    let mut unsafe_source = source;
    unsafe_source.features[0]
        .qualifiers
        .push(("transl_except".into(), Some("(pos:3..3,aa:Sec)".into())));
    assert!(s::build(&[input(unsafe_source)], &r, None).is_err());
}
#[test]
fn physical_reconstruction_is_checked_against_predeclared_internal_changes() {
    const ORIGINAL: &str = "AGTCCGATGCTACGTAGGACTTCAGACCTGTA";
    const INTENDED: &str = "AGACCGATGCTACGTAGGACTTCAGACCTGGA";
    let inputs = vec![input(
        Molecule::dna("Original32bp", ORIGINAL, Topology::Linear).unwrap(),
    )];
    let r = request(
        vec![fragment("whole", 0, 32, false)],
        vec![edit("whole", 2, "T", "A"), edit("whole", 30, "T", "G")],
    );
    let target = s::build(&inputs, &r, None).unwrap().target;
    assert_eq!(target.sequence.as_str(), INTENDED);
    let mut plan = DesignPlan {
        name: "Check independently declared substitutions".into(),
        target_revision: target.revision(),
        inputs: inputs.clone(),
        steps: vec![Operation::Pcr {
            template: Source::Input("original".into()),
            forward: Primer::new("OriginalF", "AGTCCGAT", "").unwrap(),
            reverse: Primer::new("OriginalR", "TACAGGTC", "").unwrap(),
        }],
    };
    assert_eq!(verification::verify(&plan, &target).status, Status::Invalid);
    let mut options = collagene_core::mutagenesis::PcrOptions::default();
    options.binding.max_mismatches = 1;
    options.binding.min_three_prime_exact = 2;
    plan.steps = vec![Operation::PcrWithMismatches {
        template: Source::Input("original".into()),
        forward: Primer::new("MutantF", "AGACCGAT", "").unwrap(),
        reverse: Primer::new("MutantR", "TCCAGGTC", "").unwrap(),
        options,
    }];
    let checked = verification::verify(&plan, &target);
    assert_eq!(
        checked.status,
        Status::ConsistentUnderAssumptions,
        "{:?}",
        checked.findings
    );
}
#[test]
fn stale_bases_roles_sources_report_tampering_and_hostile_counts_fail_explicitly() {
    let inputs = vec![input(
        Molecule::dna("Original8bp", "AACCGGTT", Topology::Linear).unwrap(),
    )];
    let r = request(
        vec![fragment("whole", 0, 8, false)],
        vec![edit("whole", 2, "C", "A")],
    );
    let result = s::build(&inputs, &r, None).unwrap();
    for edit in [
        edit("whole", 2, "A", "T"),
        edit("missing", 2, "C", "A"),
        edit("whole", usize::MAX, "C", "A"),
        edit("whole", 2, "c", "A"),
        edit("whole", 2, "C", "NN"),
    ] {
        let mut changed = r.clone();
        changed.substitutions = vec![edit];
        assert!(s::intent_revision(&inputs, &changed, None).is_err());
    }
    let mut duplicate = r.clone();
    duplicate
        .substitutions
        .push(duplicate.substitutions[0].clone());
    assert!(s::intent_revision(&inputs, &duplicate, None).is_err());
    let mut huge = r.clone();
    huge.substitutions = vec![edit("whole", 2, "C", "A"); 257];
    assert!(s::source_preflight(&[("original", &inputs[0].molecule)], &huge, None).is_err());
    let mut stale = inputs.clone();
    stale[0].revision = "0".repeat(64);
    assert!(s::intent_revision(&stale, &r, None).is_err());
    let mut tampered = result.clone();
    tampered.applied[0].source_position = ContributionPosition::Reference {
        position: 5,
        reverse: false,
    };
    assert!(tampered.validate_binding(&result.context_revision).is_err());
    let mut tampered = result.clone();
    tampered.target.sequence = collagene_core::Sequence::dna("AATCGGTT").unwrap();
    assert!(tampered.validate_binding(&result.context_revision).is_err());
    let mut tampered = result.clone();
    tampered.target.features[0].needs_review = false;
    assert!(tampered.validate_binding(&result.context_revision).is_err());
    let decoded: s::Report = serde_json::from_slice(&serde_json::to_vec(&result).unwrap()).unwrap();
    decoded.validate_binding(&result.context_revision).unwrap();
}

#[test]
fn normalized_circular_trimmed_cores_keep_original_oriented_offsets_on_both_strands() {
    let inputs = vec![input(
        Molecule::dna("Original8bp circle", "AACCGGTT", Topology::Circular).unwrap(),
    )];
    // Forward TTAA at6..10: incoming trim3 retains original offset3 (source1).
    let mut forward = request(
        vec![
            fragment("sibling", 0, 4, false),
            fragment("wrap", 6, 10, false),
        ],
        vec![edit("sibling", 1, "A", "T"), edit("wrap", 3, "A", "C")],
    );
    forward.junction_edits = vec![Edit {
        junction: 0,
        remove_left: 0,
        remove_right: 3,
        inserted: "GG".into(),
    }];
    let result = s::build(&inputs, &forward, None).unwrap();
    assert_eq!(result.target.sequence.as_str(), "ATCCGGC");
    assert_eq!(
        (
            result.fragments[1].retained_start,
            result.fragments[1].retained_end
        ),
        (3, 4)
    );
    assert!(matches!(
        result.fragments[1].source_position(3).unwrap(),
        ContributionPosition::Reference {
            position: 1,
            reverse: false
        }
    ));
    assert_eq!(
        (
            reference_position(&result.applied[1].source_position).0,
            result.applied[1].output_position
        ),
        (1, 6)
    );
    result
        .validate_binding(&s::intent_revision(&inputs, &forward, None).unwrap())
        .unwrap();
    let mut trimmed = forward.clone();
    trimmed.substitutions = vec![edit("wrap", 2, "A", "G")];
    assert!(s::intent_revision(&inputs, &trimmed, None).is_err());
    // Reverse TTAA: outgoing trim3 retains original reverse offset0 (source1,
    // expected T is the complement of original A). The core again normalizes1..2.
    let mut reverse = request(
        vec![
            fragment("wrap", 6, 10, true),
            fragment("sibling", 0, 4, false),
        ],
        vec![edit("wrap", 0, "T", "C"), edit("sibling", 2, "C", "T")],
    );
    reverse.junction_edits = vec![Edit {
        junction: 0,
        remove_left: 3,
        remove_right: 1,
        inserted: "GG".into(),
    }];
    let result = s::build(&inputs, &reverse, None).unwrap();
    assert_eq!(result.target.sequence.as_str(), "CGGATC");
    assert_eq!(
        (
            result.fragments[0].retained_start,
            result.fragments[0].retained_end
        ),
        (0, 1)
    );
    assert!(matches!(
        result.fragments[0].source_position(0).unwrap(),
        ContributionPosition::Reference {
            position: 1,
            reverse: true
        }
    ));
    assert_eq!(
        (
            reference_position(&result.applied[0].source_position).0,
            result.applied[0].output_position
        ),
        (1, 0)
    );
    assert!(reference_position(&result.applied[0].source_position).1);
    result
        .validate_binding(&s::intent_revision(&inputs, &reverse, None).unwrap())
        .unwrap();
    let mut trimmed = reverse.clone();
    trimmed.substitutions = vec![edit("wrap", 1, "T", "A")];
    assert!(s::intent_revision(&inputs, &trimmed, None).is_err());
}

fn reference_position(position: &ContributionPosition) -> (usize, bool) {
    match position {
        ContributionPosition::Reference { position, reverse } => (*position, *reverse),
        ContributionPosition::Material { .. } => panic!("expected an original reference position"),
    }
}
