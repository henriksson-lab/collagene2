//! Literal reference edits are source-bound declarations, never primer-derived products.
use collagene_core::{
    Feature, Location, Molecule, Topology, assembly_edit_intent as e, assembly_intent,
    assembly_intent::{Contribution, ContributionPosition, Selector},
    assembly_junction_edit, assembly_substitution_intent as s,
    verification::PlanInput,
};
fn source(dna: &str, topology: Topology) -> Vec<PlanInput> {
    let molecule = Molecule::dna("Original", dna, topology).unwrap();
    vec![PlanInput {
        id: "source".into(),
        revision: molecule.revision(),
        molecule,
    }]
}
fn fragment(role: &str, start: usize, end: usize, reverse: bool) -> Contribution {
    Contribution {
        role: role.into(),
        input_id: "source".into(),
        selector: Selector::ReferenceRange {
            start,
            end,
            reverse,
        },
    }
}
fn edit(role: &str, start: usize, end: usize, expected: &str, replacement: &str) -> e::Edit {
    e::Edit {
        role: role.into(),
        start,
        end,
        expected: expected.into(),
        replacement: replacement.into(),
    }
}
fn request(fragments: Vec<Contribution>, edits: Vec<e::Edit>) -> e::Request {
    e::Request {
        base: s::Request {
            base: assembly_intent::Request {
                material_bindings: vec![],
                name: "Literal edited intent".into(),
                fragments,
                circular: false,
            },
            junction_edits: vec![],
            substitutions: vec![],
        },
        edits,
    }
}
#[test]
fn simultaneous_original_offsets_support_insert_delete_replace_and_retained_receipts() {
    let inputs = source("AACCGGTT", Topology::Linear);
    let before = serde_json::to_vec(&inputs).unwrap();
    let req = request(
        vec![fragment("insert", 0, 8, false)],
        vec![
            edit("insert", 6, 8, "TT", "A"),
            edit("insert", 0, 0, "", "GG"),
            edit("insert", 2, 4, "CC", ""),
        ],
    );
    let result = e::build(&inputs, &req, None).unwrap();
    assert_eq!(result.target.sequence.as_str(), "GGAAGGA");
    assert_eq!(
        result
            .applied
            .iter()
            .map(|a| (a.output_start, a.output_end))
            .collect::<Vec<_>>(),
        vec![(6, 7), (0, 2), (4, 4)]
    );
    assert_eq!(
        (
            result.fragments[0].output_start,
            result.fragments[0].output_end
        ),
        (0, 7)
    );
    assert_eq!(serde_json::to_vec(&inputs).unwrap(), before);
    result.validate_binding(&result.context_revision).unwrap();
    result
        .validate_sources(&inputs, &result.context_revision, None)
        .unwrap();
    assert_eq!(
        result.target.metadata["assembly:intent:reference_only"],
        "true"
    );
}
#[test]
fn repeated_sources_reverse_orientation_and_disjoint_substitution_stay_independent() {
    let inputs = source("AACCGGTA", Topology::Linear);
    let mut req = request(
        vec![
            fragment("forward", 0, 8, false),
            fragment("reverse", 0, 8, true),
        ],
        vec![edit("reverse", 1, 3, "AC", "TTT")],
    );
    req.base.substitutions.push(s::Substitution {
        role: "forward".into(),
        offset: 0,
        expected: "A".into(),
        replacement: "G".into(),
    });
    let result = e::build(&inputs, &req, None).unwrap();
    // Reverse-complement original is TACCGGTT; original offsets1..3 are AC.
    assert_eq!(result.target.sequence.as_str(), "GACCGGTATTTTCGGTT");
    assert_eq!(
        (
            result.applied[0].original_output_start,
            result.applied[0].output_start,
            result.applied[0].output_end
        ),
        (9, 9, 12)
    );
    assert!(result.applied[0].selector.reference_interval().unwrap().2);
    assert_eq!(result.fragments[1].output_end, 17);
    result.validate_binding(&result.context_revision).unwrap();
}
#[test]
fn circular_trim_normalization_preserves_original_offsets_and_outgoing_spacer() {
    let inputs = source("AACCGGTT", Topology::Circular);
    let mut req = request(
        vec![
            fragment("head", 2, 4, false),
            fragment("wrap", 6, 10, false),
        ],
        vec![edit("wrap", 3, 4, "A", "CC")],
    );
    req.base.junction_edits.push(assembly_junction_edit::Edit {
        junction: 0,
        remove_left: 0,
        remove_right: 3,
        inserted: "GG".into(),
    });
    let result = e::build(&inputs, &req, None).unwrap();
    // Head CC, literal junction GG, retained wrap source base1 A replaced CC.
    assert_eq!(result.target.sequence.as_str(), "CCGGCC");
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
            result.applied[0].original_output_start,
            result.applied[0].output_start
        ),
        (4, 4)
    );
    let mut trimmed = req.clone();
    trimmed.edits[0] = edit("wrap", 2, 3, "A", "C");
    assert!(e::build(&inputs, &trimmed, None).is_err());
    result.validate_binding(&result.context_revision).unwrap();
}
#[test]
fn spanning_cds_includes_inserted_bases_and_drops_stale_translation() {
    let mut inputs = source("AACCGGTT", Topology::Linear);
    inputs[0].molecule.features.push(Feature {
        kind: "CDS".into(),
        location: Location::range(1, 7),
        qualifiers: vec![("translation".into(), Some("stale".into()))],
        needs_review: false,
    });
    inputs[0].revision = inputs[0].molecule.revision();
    let req = request(
        vec![fragment("gene", 0, 8, false)],
        vec![edit("gene", 4, 4, "", "AAA")],
    );
    let result = e::build(&inputs, &req, None).unwrap();
    assert_eq!(result.target.sequence.as_str(), "AACCAAAGGTT");
    let cds = result
        .target
        .features
        .iter()
        .find(|f| f.kind == "CDS")
        .unwrap();
    assert_eq!(cds.location, Location::range(1, 10));
    assert!(cds.needs_review);
    assert!(!cds.qualifiers.iter().any(|(k, _)| k == "translation"));
    result
        .validate_sources(&inputs, &result.context_revision, None)
        .unwrap();
}
#[test]
fn stale_literals_conflicts_boundaries_and_oversized_edits_are_rejected() {
    let inputs = source("AACCGGTT", Topology::Linear);
    let base = request(
        vec![fragment("one", 0, 8, false)],
        vec![edit("one", 2, 4, "CC", "")],
    );
    let mut bad = base.clone();
    bad.edits[0].expected = "AA".into();
    assert!(e::build(&inputs, &bad, None).is_err());
    for extra in [
        edit("one", 2, 2, "", "G"),
        edit("one", 4, 4, "", "G"),
        edit("one", 3, 5, "CG", "T"),
    ] {
        let mut bad = base.clone();
        bad.edits.push(extra);
        assert!(e::build(&inputs, &bad, None).is_err());
    }
    let mut bad = base.clone();
    bad.base.substitutions.push(s::Substitution {
        role: "one".into(),
        offset: 3,
        expected: "C".into(),
        replacement: "A".into(),
    });
    assert!(e::build(&inputs, &bad, None).is_err());
    let mut bad = base.clone();
    bad.edits[0].replacement = "A".repeat(100001);
    assert!(e::build(&inputs, &bad, None).is_err());
    let mut bad = base.clone();
    bad.edits[0].role = "missing".into();
    assert!(e::build(&inputs, &bad, None).is_err());
    let mut bad = base.clone();
    bad.edits = vec![edit("one", 0, 8, "AACCGGTT", "")];
    assert!(e::build(&inputs, &bad, None).is_err());
    let mut adjacent = base.clone();
    adjacent.edits.push(edit("one", 4, 6, "GG", "A"));
    assert_eq!(
        e::build(&inputs, &adjacent, None)
            .unwrap()
            .target
            .sequence
            .as_str(),
        "AAATT"
    );
}
#[test]
fn changed_source_target_receipts_and_context_cannot_revalidate() {
    let inputs = source("AACCGGTT", Topology::Linear);
    let req = request(
        vec![fragment("one", 0, 8, false)],
        vec![edit("one", 2, 4, "CC", "A")],
    );
    let result = e::build(&inputs, &req, None).unwrap();
    let mut changed = inputs.clone();
    changed[0].molecule.name = "Renamed original".into();
    changed[0].revision = changed[0].molecule.revision();
    assert!(
        result
            .validate_sources(&changed, &result.context_revision, None)
            .is_err()
    );
    let mut changed = result.clone();
    changed.target.sequence = collagene_core::Sequence::dna("AAAGGTA").unwrap();
    assert!(changed.validate_binding(&result.context_revision).is_err());
    let mut changed = result.clone();
    changed.applied[0].output_start += 1;
    assert!(changed.validate_binding(&result.context_revision).is_err());
    let mut changed = result.clone();
    changed.fragments[0].output_end += 1;
    assert!(changed.validate_binding(&result.context_revision).is_err());
    let mut forged = result.clone();
    forged.fragments[0].retained_start += 1;
    let mut descriptors: serde_json::Value =
        serde_json::from_str(&forged.target.metadata["assembly_edit_intent:contributions"])
            .unwrap();
    descriptors[0]["retained_start"] = serde_json::json!(forged.fragments[0].retained_start);
    forged.target.metadata.insert(
        "assembly_edit_intent:contributions".into(),
        serde_json::to_string(&descriptors).unwrap(),
    );
    assert!(forged.validate_binding(&result.context_revision).is_err());
    assert!(result.validate_binding(&"0".repeat(64)).is_err());
}

#[test]
fn protein_alphabet_cannot_masquerade_as_dna_with_identical_letters_and_receipts() {
    let inputs = source("AACCGGTT", Topology::Linear);
    let req = request(
        vec![fragment("one", 0, 8, false)],
        vec![edit("one", 2, 4, "CC", "A")],
    );
    let report = e::build(&inputs, &req, None).unwrap();
    let mut forged = report.clone();
    forged.target.sequence = collagene_core::Sequence::new(
        collagene_core::Alphabet::Protein,
        report.target.sequence.as_str(),
    )
    .unwrap();
    assert_eq!(
        forged.target.sequence.as_str(),
        report.target.sequence.as_str()
    );
    assert_eq!(forged.target.metadata, report.target.metadata);
    assert!(forged.validate_binding(&report.context_revision).is_err());
    report.validate_binding(&report.context_revision).unwrap();
}
