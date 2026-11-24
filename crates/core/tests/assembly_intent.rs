use collagene_core::{
    Feature, Location, Molecule, Topology, assembly_intent as i, overlap_primer_design::Fragment,
    verification::PlanInput,
};
fn input(id: &str, dna: &str) -> PlanInput {
    let molecule = Molecule::dna(id, dna, Topology::Linear).unwrap();
    PlanInput {
        id: id.into(),
        revision: molecule.revision(),
        molecule,
    }
}
fn row(role: &str, id: &str, start: usize, end: usize, reverse: bool) -> Fragment {
    Fragment {
        role: role.into(),
        input_id: id.into(),
        start,
        end,
        reverse,
    }
}
#[test]
fn literal_contributions_preserve_safe_annotations_and_originals_without_simulation() {
    let mut inputs = vec![input("a", "AAACCCGGGTTT"), input("b", "ATGCCCAA")];
    inputs[0].molecule.features.push(Feature {
        kind: "misc_feature".into(),
        location: Location::range(3, 6),
        qualifiers: vec![("label".into(), Some("retained".into()))],
        needs_review: false,
    });
    inputs[0].molecule.original_genbank =
        Some("original raw source remains external evidence".into());
    inputs[0].revision = inputs[0].molecule.revision();
    let request = i::Request {
        name: "frozen independent intended reference".into(),
        fragments: vec![
            row("insert", "a", 3, 9, false),
            row("reverse arm", "b", 0, 4, true),
        ],
        circular: false,
    };
    let before = serde_json::to_value(&inputs).unwrap();
    let report = i::build(&inputs, &request).unwrap();
    assert_eq!(report.target.sequence.as_str(), "CCCGGGGCAT");
    assert_eq!(
        (
            report.fragments[0].output_start,
            report.fragments[0].output_end
        ),
        (0, 6)
    );
    assert_eq!(
        (
            report.fragments[1].output_start,
            report.fragments[1].output_end
        ),
        (6, 10)
    );
    assert_eq!(report.target.features[0].location, Location::range(0, 3));
    assert_eq!(
        report.target.features[0]
            .location
            .extract(&report.target.sequence)
            .unwrap(),
        "CCC"
    );
    assert_eq!(serde_json::to_value(&inputs).unwrap(), before);
    assert!(report.target.original_genbank.is_none());
    assert_eq!(
        report.target.metadata["assembly:intent:reference_only"],
        "true"
    );
    assert!(collagene_core::chemistry::Duplex::from_molecule(&report.target).is_err());
    report
        .validate_binding(&i::intent_revision(&inputs, &request).unwrap())
        .unwrap();
}
#[test]
fn circular_origin_and_reverse_are_explicit_reference_contributions() {
    let mut inputs = vec![input("circle", "AACCGGTT")];
    inputs[0].molecule.topology = Topology::Circular;
    inputs[0].revision = inputs[0].molecule.revision();
    let request = i::Request {
        name: "circular intended origin".into(),
        fragments: vec![row("origin-arm", "circle", 5, 10, true)],
        circular: true,
    };
    let report = i::build(&inputs, &request).unwrap();
    assert_eq!(report.target.sequence.as_str(), "TTAAC");
    assert_eq!(report.target.topology, Topology::Circular);
    assert_eq!(report.fragments[0].output_end, 5);
    assert_eq!(inputs[0].molecule.sequence.as_str(), "AACCGGTT");
    assert!(report.warnings.iter().any(|s| s.contains("not a PCR")));
}
#[test]
fn readable_physical_source_and_iupac_remain_reference_only_and_immutable() {
    let mut inputs = vec![input("reference", "ACGN")];
    inputs[0].molecule.metadata.insert(
        "duplex_model".into(),
        "historical physical producer evidence".into(),
    );
    inputs[0]
        .molecule
        .metadata
        .insert("variant:reference_only".into(), "true".into());
    inputs[0].revision = inputs[0].molecule.revision();
    let before = serde_json::to_value(&inputs).unwrap();
    let request = i::Request {
        name: "reference intent".into(),
        fragments: vec![row("readable", "reference", 0, 4, false)],
        circular: false,
    };
    let report = i::build(&inputs, &request).unwrap();
    assert_eq!(report.target.sequence.as_str(), "ACGN");
    assert!(report.warnings.iter().any(|w| w.contains("IUPAC")));
    assert_eq!(serde_json::to_value(&inputs).unwrap(), before);
    assert_eq!(
        report.target.metadata["duplex_model"],
        "Reference-only assembly design intent; physical material state unavailable"
    );
    assert!(
        !report
            .target
            .metadata
            .contains_key("variant:reference_only")
    );
}
#[test]
fn ranges_revisions_requests_and_result_receipts_cannot_change_silently() {
    let mut inputs = vec![input("a", "AAACCCGGGTTT")];
    let request = i::Request {
        name: "explicit ranges".into(),
        fragments: vec![row("piece", "a", 0, 6, false)],
        circular: false,
    };
    let report = i::build(&inputs, &request).unwrap();
    let mut bad = request.clone();
    bad.fragments[0].end = 13;
    assert!(i::build(&inputs, &bad).is_err());
    bad = request.clone();
    bad.fragments[0].start = 6;
    assert!(i::build(&inputs, &bad).is_err());
    bad = request.clone();
    bad.fragments[0].end = 7;
    assert_ne!(
        i::intent_revision(&inputs, &bad).unwrap(),
        report.context_revision
    );
    let mut changed = report.clone();
    changed.fragments[0].output_end = 5;
    assert!(changed.validate().is_err());
    changed = report;
    changed.target.sequence = collagene_core::Sequence::dna("CCCCCC").unwrap();
    assert!(changed.validate().is_err());
    inputs[0].molecule.name = "advanced actual document".into();
    assert!(i::build(&inputs, &request).is_err());
}
#[test]
fn clipped_or_unresolved_annotations_are_reviewed_or_retained_as_source_evidence() {
    let mut inputs = vec![input("a", "AAACCCGGGTTT")];
    inputs[0].molecule.features = vec![
        Feature {
            kind: "CDS".into(),
            location: Location::range(0, 9),
            qualifiers: vec![
                ("label".into(), Some("clipped CDS".into())),
                ("codon_start".into(), Some("1".into())),
                ("translation".into(), Some("KPG".into())),
            ],
            needs_review: false,
        },
        Feature {
            kind: "misc_feature".into(),
            location: Location::Unresolved("unplaced original".into()),
            qualifiers: vec![("label".into(), Some("unmapped".into()))],
            needs_review: true,
        },
    ];
    inputs[0].revision = inputs[0].molecule.revision();
    let request = i::Request {
        name: "partial contribution".into(),
        fragments: vec![row("partial", "a", 3, 6, false)],
        circular: false,
    };
    let report = i::build(&inputs, &request).unwrap();
    assert_eq!(report.target.sequence.as_str(), "CCC");
    assert!(
        report
            .target
            .features
            .iter()
            .filter(|f| f.kind == "CDS")
            .all(|f| f.needs_review
                && !f
                    .qualifiers
                    .iter()
                    .any(|(k, _)| k == "translation" || k == "codon_start"))
    );
    assert!(
        report
            .warnings
            .iter()
            .any(|w| w.contains("review") || w.contains("withheld"))
    );
    assert_eq!(inputs[0].molecule.features.len(), 2);
}
