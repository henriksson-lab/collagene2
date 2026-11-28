use collagene_core::cloning::pcr::mutation_intent::*;
use collagene_model::{Feature, Location, Molecule, Topology};
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
    let options = collagene_core::cloning::pcr::mutagenesis::InversePcrOptions {
        edit_start: 24,
        edit_end: 24,
        inserted: "GGA".into(),
        forward_length: 8,
        reverse_length: 8,
        reverse_tail_bases: 0,
    };
    let designed = collagene_core::cloning::pcr::mutagenesis::design_inverse_pcr(
        &source24, &target27, None, &options,
    )
    .unwrap();
    assert_eq!(
        designed.verification.status,
        collagene_core::verification::evidence::Status::ConsistentUnderAssumptions,
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
    let affected = build(&s, &request(4, 8, "AT")).unwrap();
    assert!(affected.target.features[0].needs_review);
    assert!(values(&affected.target.features[0], "transl_except").is_empty());
    assert!(
        values(&affected.target.features[0], "note")
            .iter()
            .any(|note| note.contains("(pos:1..3,aa:Sec)"))
    );
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
    bad.target.sequence = collagene_model::Sequence::dna("AAAAGTGGGGTTTT").unwrap();
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
        collagene_model::Sequence::new(collagene_model::Alphabet::Protein, "AAAA").unwrap();
    assert!(build(&protein, &request(0, 1, "A")).is_err());
}

fn values<'a>(feature: &'a Feature, key: &str) -> Vec<&'a str> {
    feature
        .qualifiers
        .iter()
        .filter(|(k, _)| k == key)
        .filter_map(|(_, v)| v.as_deref())
        .collect()
}
fn coordinate_source(topology: Topology) -> Molecule {
    let mut s = source(topology);
    s.features = vec![
        Feature {
            kind: "CDS".into(),
            location: Location::range(0, 6),
            qualifiers: vec![
                ("transl_except".into(), Some("(pos:1..3,aa:Sec)".into())),
                ("translation".into(), Some("UT".into())),
                ("codon_start".into(), Some("1".into())),
            ],
            needs_review: false,
        },
        Feature {
            kind: "tRNA".into(),
            location: Location::range(8, 11),
            qualifiers: vec![(
                "anticodon".into(),
                Some("(pos:9..11,aa:Gly,seq:ccc)".into()),
            )],
            needs_review: false,
        },
    ];
    s
}
fn archived_features(report: &Report, source: &Molecule) -> Vec<Feature> {
    serde_json::from_str(
        &report.target.metadata[&format!("mutation_intent:source_features:{}", source.revision())],
    )
    .unwrap()
}
#[test]
fn supported_untouched_claims_survive_literal_edits_and_prefix_coordinate_shifts() {
    let s = coordinate_source(Topology::Linear);
    let original = s.clone();
    let outside = build(&s, &request(12, 14, "AC")).unwrap();
    assert_eq!(outside.target.sequence.as_str(), "AAAACCCCGGGGACTT");
    assert_eq!(outside.target.features, s.features);
    let shifted = build(&s, &request(0, 0, "TTT")).unwrap();
    assert_eq!(shifted.target.sequence.as_str(), "TTTAAAACCCCGGGGTTTT");
    assert_eq!(shifted.target.features[0].location, Location::range(3, 9));
    assert_eq!(shifted.target.features[1].location, Location::range(11, 14));
    assert_eq!(
        values(&shifted.target.features[0], "transl_except"),
        ["(pos:4..6,aa:Sec)"]
    );
    assert_eq!(
        values(&shifted.target.features[1], "anticodon"),
        ["(pos:12..14,aa:Gly,seq:ccc)"]
    );
    assert!(
        shifted
            .target
            .features
            .iter()
            .all(|feature| !feature.needs_review)
    );
    assert_eq!(archived_features(&shifted, &s), s.features);
    assert!(collagene_core::chemistry::Duplex::from_molecule(&shifted.target).is_err());
    assert_eq!(s, original);
}
#[test]
fn touched_codon_and_anticodon_claims_are_source_evidence_not_active_destination_assertions() {
    let s = coordinate_source(Topology::Linear);
    let original = s.clone();
    let codon = build(&s, &request(1, 2, "T")).unwrap();
    assert_eq!(codon.target.sequence.as_str(), "ATAACCCCGGGGTTTT");
    assert!(codon.target.features[0].needs_review);
    for key in ["transl_except", "translation", "codon_start"] {
        assert!(values(&codon.target.features[0], key).is_empty());
    }
    assert!(
        values(&codon.target.features[0], "note")
            .iter()
            .any(|value| value.contains("(pos:1..3,aa:Sec)"))
    );
    assert_eq!(codon.target.features[1], s.features[1]);
    let anticodon = build(&s, &request(9, 10, "A")).unwrap();
    assert_eq!(anticodon.target.sequence.as_str(), "AAAACCCCGAGGTTTT");
    assert_eq!(anticodon.target.features[0], s.features[0]);
    assert!(anticodon.target.features[1].needs_review);
    assert!(values(&anticodon.target.features[1], "anticodon").is_empty());
    assert!(
        values(&anticodon.target.features[1], "note")
            .iter()
            .any(|value| value.contains("(pos:9..11,aa:Gly,seq:ccc)"))
    );
    assert_eq!(archived_features(&codon, &s), s.features);
    assert_eq!(archived_features(&anticodon, &s), s.features);
    assert_eq!(s, original);
}
#[test]
fn unsupported_retained_claims_still_refuse_but_wholly_deleted_features_need_no_mapping() {
    for (key, value) in [
        ("rpt_unit_range", "9..11"),
        ("coded_by", "AB123456.1:9..11"),
        ("transl_except", "not a coordinate wrapper"),
    ] {
        let mut s = coordinate_source(Topology::Linear);
        s.features[1].qualifiers = vec![(key.into(), Some(value.into()))];
        let original = s.clone();
        assert!(build(&s, &request(12, 14, "AC")).is_err(), "retained {key}");
        let deleted = build(&s, &request(8, 11, "")).unwrap();
        assert_eq!(deleted.target.sequence.as_str(), "AAAACCCCGTTTT");
        assert_eq!(deleted.target.features, vec![s.features[0].clone()]);
        assert_eq!(archived_features(&deleted, &s), s.features);
        assert_eq!(s, original);
    }
}
#[test]
fn circular_origin_wrap_projects_supported_claims_through_the_canonical_rotation_and_edit() {
    let mut s = coordinate_source(Topology::Circular);
    s.features[0].location = Location::range(4, 10);
    s.features[0].qualifiers[0].1 = Some("(pos:5..7,aa:Sec)".into());
    s.features[0].qualifiers[1].1 = Some("UR".into());
    let original = s.clone();
    let r = build(&s, &request(14, 18, "G")).unwrap();
    assert_eq!(r.target.sequence.as_str(), "GAACCCCGGGGTT");
    assert_eq!(r.deleted, "TTAA");
    assert_eq!(r.origin_shift, 14);
    assert_eq!(r.target.features[0].location, Location::range(3, 9));
    assert_eq!(r.target.features[1].location, Location::range(7, 10));
    assert_eq!(
        values(&r.target.features[0], "transl_except"),
        ["(pos:4..6,aa:Sec)"]
    );
    assert_eq!(
        values(&r.target.features[1], "anticodon"),
        ["(pos:8..10,aa:Gly,seq:ccc)"]
    );
    assert!(
        r.target
            .features
            .iter()
            .all(|feature| !feature.needs_review)
    );
    assert_eq!(archived_features(&r, &s), s.features);
    let mut stale = r.clone();
    stale
        .target
        .metadata
        .insert("mutation_intent:engine".into(), "mutation-intent-3".into());
    assert!(stale.validate_binding(&r.context_revision).is_err());
    assert!(collagene_core::chemistry::Duplex::from_molecule(&r.target).is_err());
    assert_eq!(s, original);
}
