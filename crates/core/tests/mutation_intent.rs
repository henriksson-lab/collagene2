use collagene_core::{Feature, Location, Molecule, Topology, mutation_intent::*};
fn request(start: usize, end: usize, inserted: &str) -> Request {
    Request {
        name: "explicit mutation target".into(),
        edit_start: start,
        edit_end: end,
        inserted: inserted.into(),
    }
}
fn source(topology: Topology) -> Molecule {
    Molecule::dna("original", "AAAACCCCGGGGTTTT", topology).unwrap()
}
#[test]
fn linear_substitution_insertion_and_deletion_are_literal_independent_references() {
    let s = source(Topology::Linear);
    let before = serde_json::to_value(&s).unwrap();
    for (r, expected, deleted) in [
        (request(4, 8, "AT"), "AAAAATGGGGTTTT", "CCCC"),
        (request(4, 4, "AT"), "AAAAATCCCCGGGGTTTT", ""),
        (request(4, 8, ""), "AAAAGGGGTTTT", "CCCC"),
        (request(16, 16, "AC"), "AAAACCCCGGGGTTTTAC", ""),
    ] {
        let report = build(&s, &r).unwrap();
        assert_eq!(report.target.sequence.as_str(), expected);
        assert_eq!(report.deleted, deleted);
        assert_eq!(report.origin_shift, 0);
        assert_eq!(report.source_revision, s.revision());
        report
            .validate_binding(&intent_revision(&s, &r).unwrap())
            .unwrap();
        assert_eq!(
            report.target.metadata["mutation:intent:reference_only"],
            "true"
        );
    }
    assert_eq!(serde_json::to_value(&s).unwrap(), before);
}
#[test]
fn circular_wrap_has_explicit_new_origin_and_retained_source_path() {
    let s = source(Topology::Circular);
    let r = build(&s, &request(14, 18, "GC")).unwrap();
    assert_eq!(r.deleted, "TTAA");
    assert_eq!(r.target.sequence.as_str(), "GCAACCCCGGGGTT");
    assert_eq!(r.origin_shift, 14);
    assert_eq!(r.target.topology, Topology::Circular);
    let terminal_caret = build(&s, &request(16, 16, "GC")).unwrap();
    assert_eq!(
        terminal_caret.target.sequence.as_str(),
        "AAAACCCCGGGGTTTTGC"
    );
    assert_eq!(terminal_caret.origin_shift, 0);
    assert!(terminal_caret.deleted.is_empty());
    assert!(build(&s, &request(16, 17, "A")).is_err());
    let source24 = Molecule::dna(
        "24nt source",
        "AGTCCGATGCTACGTAGGACTTCA",
        Topology::Circular,
    )
    .unwrap();
    let target27 = Molecule::dna(
        "independent 27nt target",
        "AGTCCGATGCTACGTAGGACTTCAGGA",
        Topology::Circular,
    )
    .unwrap();
    let options = collagene_core::mutagenesis::InversePcrOptions {
        edit_start: 24,
        edit_end: 24,
        inserted: "GGA".into(),
        forward_length: 8,
        reverse_length: 8,
        reverse_tail_bases: 0,
    };
    let designed =
        collagene_core::mutagenesis::design_inverse_pcr(&source24, &target27, &options).unwrap();
    assert_eq!(
        designed.verification.status,
        collagene_core::verification::Status::ConsistentUnderAssumptions,
        "{:?}",
        designed.verification.findings
    );
    assert_eq!(designed.options.edit_start, 24);
    assert_eq!(designed.options.edit_end, 24);
    let internal = build(&s, &request(4, 8, "AT")).unwrap();
    assert_eq!(internal.target.sequence.as_str(), "AAAAATGGGGTTTT");
    assert_eq!(internal.origin_shift, 0);
    let full = build(&s, &request(3, 19, "AC")).unwrap();
    assert_eq!(full.target.sequence.as_str(), "AC");
    assert_eq!(full.origin_shift, 3);
    assert_eq!(full.deleted, "ACCCCGGGGTTTTAAA");
    assert!(build(&s, &request(3, 19, "")).is_err());
}
#[test]
fn affected_annotations_are_reviewed_and_original_feature_evidence_is_retained() {
    let mut s = source(Topology::Linear);
    s.features = vec![
        Feature {
            kind: "CDS".into(),
            location: Location::range(0, 12),
            qualifiers: vec![
                ("label".into(), Some("gene".into())),
                ("translation".into(), Some("old protein".into())),
            ],
            needs_review: false,
        },
        Feature {
            kind: "misc_feature".into(),
            location: Location::range(4, 8),
            qualifiers: vec![("label".into(), Some("removed region".into()))],
            needs_review: false,
        },
    ];
    s.original_genbank =
        Some("LOCUS       original 16 bp DNA linear\nORIGIN\n 1 aaaaccccggggtttt\n//\n".into());
    let before = s.clone();
    let r = build(&s, &request(4, 8, "")).unwrap();
    assert_eq!(r.target.features.len(), 1);
    assert!(r.target.features[0].needs_review);
    assert!(
        !r.target.features[0]
            .qualifiers
            .iter()
            .any(|(k, _)| k == "translation")
    );
    let evidence: Vec<Feature> = serde_json::from_str(
        &r.target.metadata[&format!("mutation_intent:source_features:{}", s.revision())],
    )
    .unwrap();
    assert_eq!(evidence, before.features);
    assert_eq!(r.target.original_genbank, before.original_genbank);
    assert_eq!(s, before);
    s.features[0]
        .qualifiers
        .push(("transl_except".into(), Some("(pos:1..3,aa:Sec)".into())));
    assert!(build(&s, &request(4, 8, "AT")).is_err());
}
#[test]
fn readable_physical_reference_can_supply_intent_but_never_material_authority() {
    let mut s = source(Topology::Circular);
    s.metadata
        .insert("duplex_model".into(), "original physical evidence".into());
    s.metadata
        .insert("source:snapgene:chemistry_unmapped".into(), "true".into());
    let before = s.clone();
    let r = build(&s, &request(14, 18, "GC")).unwrap();
    assert_eq!(r.target.sequence.as_str(), "GCAACCCCGGGGTT");
    assert_eq!(
        r.target.metadata[&format!("mutation_intent:source_duplex_model:{}", s.revision())],
        "original physical evidence"
    );
    assert_eq!(
        r.target.metadata["source:snapgene:chemistry_unmapped"],
        "true"
    );
    assert!(collagene_core::chemistry::Duplex::from_molecule(&r.target).is_err());
    assert_eq!(s, before);
}
#[test]
fn stale_request_and_tampered_result_fail_readonly_binding_checks() {
    let s = source(Topology::Linear);
    let r = request(4, 8, "GC");
    let report = build(&s, &r).unwrap();
    let context = intent_revision(&s, &r).unwrap();
    assert_ne!(context, intent_revision(&s, &request(4, 8, "GT")).unwrap());
    let mut changed = s.clone();
    changed.name = "another source revision".into();
    assert_ne!(context, intent_revision(&changed, &r).unwrap());
    let mut bad = report.clone();
    bad.request.inserted = "GT".into();
    assert!(bad.validate_binding(&context).is_err());
    let mut bad = report.clone();
    bad.deleted = "AAAA".into();
    assert!(bad.validate_binding(&context).is_err());
    let mut bad = report.clone();
    bad.target.sequence = collagene_core::Sequence::dna("AAAAGTGGGGTTTT").unwrap();
    assert!(bad.validate_binding(&context).is_err());
    let mut bad = report;
    bad.target.topology = Topology::Circular;
    assert!(bad.validate_binding(&context).is_err());
}
#[test]
fn explicit_bounds_and_unsupported_alphabets_fail_before_building_an_intent() {
    let linear = source(Topology::Linear);
    for r in [
        request(5, 4, "A"),
        request(0, 17, "A"),
        request(0, 16, ""),
        request(0, 0, &"A".repeat(100_001)),
        request(0, 0, &" ".repeat(200_001)),
    ] {
        assert!(source_preflight(&linear, &r).is_err());
    }
    let circle = source(Topology::Circular);
    for r in [request(1, 18, "A"), request(usize::MAX, usize::MAX, "A")] {
        assert!(source_preflight(&circle, &r).is_err());
    }
    let mut protein = linear;
    protein.sequence =
        collagene_core::Sequence::new(collagene_core::Alphabet::Protein, "AAAA").unwrap();
    assert!(build(&protein, &request(0, 1, "A")).is_err());
}
