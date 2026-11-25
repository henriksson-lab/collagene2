use collagene_core::{
    Feature, Location, Molecule, Topology,
    annotation_batch::{self as batch, ColorEdit, Request},
};

fn feature(location: Location) -> Feature {
    Feature {
        kind: "CDS".into(),
        location,
        qualifiers: vec![
            ("label".into(), Some("retained feature".into())),
            ("translation".into(), Some("MK".into())),
            ("codon_start".into(), Some("1".into())),
            ("vendor_unknown".into(), Some("first".into())),
            ("vendor_unknown".into(), Some("second".into())),
            ("vendor_flag".into(), None),
        ],
        needs_review: false,
    }
}
fn source() -> Molecule {
    let mut source = Molecule::dna("batch fixture", "ATGAAATAACCC", Topology::Circular).unwrap();
    source.original_genbank = Some("original source evidence".into());
    source
        .metadata
        .insert("unknown evidence".into(), "untouched".into());
    source.features = vec![
        feature(Location::range(0, 9)),
        feature(Location::range(9, 12)),
    ];
    source
}
fn request(source: &Molecule) -> Request {
    Request {
        source_revision: source.revision(),
        indices: vec![0],
        kind: None,
        flip_strand: false,
        color: ColorEdit::Keep,
    }
}

#[test]
fn type_color_batch_is_atomic_and_preserves_unrelated_data_and_multiplicity() {
    let source = source();
    let mut r = request(&source);
    r.kind = Some("misc_feature".into());
    r.color = ColorEdit::Set("#A0b1C2".into());
    let edited = batch::apply(&source, &r).unwrap();
    let mut expected = source.clone();
    expected.features[0].kind = "misc_feature".into();
    expected.features[0]
        .qualifiers
        .push(("color".into(), Some("#a0b1c2".into())));
    assert_eq!(edited, expected);
    assert_eq!(
        batch::feature_rgb(&edited.features[0]),
        Some([160, 177, 194])
    );
    assert_eq!(source.features[0].kind, "CDS");
    assert_eq!(source.features[0].qualifiers.len(), 6);
    assert_ne!(edited.revision(), source.revision());
}

#[test]
fn strand_flip_toggles_outer_complement_without_flattening_mixed_fuzzy_tree() {
    let mut source = source();
    let tree = Location::Join(vec![
        Location::Range {
            start: 9,
            end: 12,
            before: true,
            after: false,
        },
        Location::Complement(Box::new(Location::range(0, 6))),
    ]);
    source.features[0].location = tree.clone();
    source.features[0]
        .qualifiers
        .push(("TRANSLATION".into(), Some("another assertion".into())));
    source.features[0]
        .qualifiers
        .push(("CODON_START".into(), None));
    let mut r = request(&source);
    r.flip_strand = true;
    let edited = batch::apply(&source, &r).unwrap();
    assert_eq!(
        edited.features[0].location,
        Location::Complement(Box::new(tree.clone()))
    );
    assert!(edited.features[0].needs_review);
    assert!(
        !edited.features[0]
            .qualifiers
            .iter()
            .any(|(key, _)| key.eq_ignore_ascii_case("translation")
                || key.eq_ignore_ascii_case("codon_start"))
    );
    assert_eq!(
        edited.features[0]
            .qualifiers
            .iter()
            .filter(|(key, _)| key == "note")
            .count(),
        4
    );
    assert!(
        edited.features[0]
            .qualifiers
            .iter()
            .any(|(_, value)| value.as_ref().is_some_and(|text| text
                .contains("another assertion")
                && text.contains("source evidence")))
    );
    assert_eq!(edited.features[1], source.features[1]);
    assert_eq!(edited.sequence, source.sequence);
    let mut second = request(&edited);
    second.flip_strand = true;
    let flipped_back = batch::apply(&edited, &second).unwrap();
    assert_eq!(flipped_back.features[0].location, tree);
    assert!(flipped_back.features[0].needs_review);
    assert_eq!(
        flipped_back.features[0].qualifiers,
        edited.features[0].qualifiers
    );
}

#[test]
fn exact_strand_reversal_changes_extraction_without_changing_reference() {
    let source = source();
    assert_eq!(
        source.features[0]
            .location
            .extract(&source.sequence)
            .unwrap(),
        "ATGAAATAA"
    );
    let mut r = request(&source);
    r.flip_strand = true;
    let edited = batch::apply(&source, &r).unwrap();
    assert_eq!(
        edited.features[0]
            .location
            .extract(&edited.sequence)
            .unwrap(),
        "TTATTTCAT"
    );
    assert_eq!(edited.sequence, source.sequence);
}

#[test]
fn coordinate_linked_qualifiers_block_only_strand_edits_and_do_not_partially_apply() {
    for key in ["transl_except", "ANTICODON", "rpt_unit_range", "Coded_By"] {
        let mut source = source();
        source.features[1]
            .qualifiers
            .push((key.into(), Some("(pos:10..12,aa:Sec)".into())));
        let original = source.clone();
        let mut r = request(&source);
        r.indices = vec![0, 1];
        r.kind = Some("gene".into());
        r.color = ColorEdit::Set("#112233".into());
        r.flip_strand = true;
        assert!(batch::apply(&source, &r).is_err());
        assert_eq!(source, original);
        r.flip_strand = false;
        let edited = batch::apply(&source, &r).unwrap();
        assert!(
            edited.features[1]
                .qualifiers
                .contains(&(key.into(), Some("(pos:10..12,aa:Sec)".into())))
        );
    }
}

#[test]
fn color_replacement_and_clear_touch_only_declared_color_keys() {
    let mut source = source();
    for (key, value) in [
        ("Color", "#010203"),
        ("COLOR", "unsafe CSS"),
        ("ApEinfo_fwdcolor", "#123456"),
        ("apeinfo_REVCOLOR", "#654321"),
        ("color_extra", "keep exactly"),
    ] {
        source.features[0]
            .qualifiers
            .push((key.into(), Some(value.into())));
    }
    let mut r = request(&source);
    r.color = ColorEdit::Set("#Ff0012".into());
    let edited = batch::apply(&source, &r).unwrap();
    assert_eq!(
        edited.features[0]
            .qualifiers
            .iter()
            .filter(|(key, _)| key.eq_ignore_ascii_case("color"))
            .count(),
        1
    );
    assert_eq!(batch::feature_rgb(&edited.features[0]), Some([255, 0, 18]));
    assert!(
        edited.features[0]
            .qualifiers
            .contains(&("color_extra".into(), Some("keep exactly".into())))
    );
    let mut clear = request(&source);
    clear.color = ColorEdit::Clear;
    let cleared = batch::apply(&source, &clear).unwrap();
    let mut expected = source.features[0].qualifiers[..6].to_vec();
    expected.push(("color_extra".into(), Some("keep exactly".into())));
    assert_eq!(cleared.features[0].qualifiers, expected);
    assert_eq!(batch::feature_rgb(&cleared.features[0]), None);
    assert_eq!(cleared.features[1], source.features[1]);
}

#[test]
fn stale_duplicate_empty_invalid_and_oversized_requests_fail_atomically() {
    let source = source();
    let original = source.clone();
    let mut cases = Vec::new();
    let mut r = request(&source);
    r.source_revision = "0".repeat(64);
    cases.push(r);
    let mut r = request(&source);
    r.indices = vec![0, 0];
    cases.push(r);
    let mut r = request(&source);
    r.indices.clear();
    cases.push(r);
    let mut r = request(&source);
    r.indices = vec![0, usize::MAX];
    cases.push(r);
    let mut r = request(&source);
    r.indices = vec![0; 4097];
    cases.push(r);
    let mut r = request(&source);
    r.kind = Some("bad type".into());
    cases.push(r);
    let mut r = request(&source);
    r.kind = Some(String::new());
    cases.push(r);
    let mut r = request(&source);
    r.kind = Some("a".repeat(65));
    cases.push(r);
    let mut r = request(&source);
    r.color = ColorEdit::Set("url(evil)".into());
    cases.push(r);
    for r in cases {
        assert!(batch::apply(&source, &r).is_err());
        assert_eq!(source, original);
    }
}

#[test]
fn borrowed_source_is_bounded_before_revision_hashing_or_cloning() {
    let mut source = source();
    let r = request(&source);
    source.features = vec![source.features[0].clone(); 4097];
    assert!(batch::apply(&source, &r).is_err());
    let mut source = super_deep_source();
    assert!(batch::apply(&source, &r).is_err());
    // Avoid recursive drop of a deliberately invalid tree in a future model.
    source.features.clear();
}

fn super_deep_source() -> Molecule {
    let mut source = source();
    for _ in 0..33 {
        source.features[0].location =
            Location::Complement(Box::new(source.features[0].location.clone()));
    }
    source
}

#[test]
fn strict_hex_helpers_reject_style_text_and_unsupported_color_syntax() {
    assert_eq!(batch::parse_rgb("#aBcD09"), Some([171, 205, 9]));
    for invalid in [
        "red",
        "#fff",
        "#ffffff00",
        " #010203",
        "#010203 ",
        "rgb(1,2,3)",
        "#GG0000",
        "#aé001",
        "#000000;display:none",
        "url(javascript:evil)",
    ] {
        assert_eq!(batch::parse_rgb(invalid), None);
    }
    let mut feature = feature(Location::range(0, 3));
    feature
        .qualifiers
        .push(("CoLoR".into(), Some("url(evil)".into())));
    assert_eq!(batch::feature_rgb(&feature), None);
    feature
        .qualifiers
        .push(("ApEinfo_FwdColor".into(), Some("#010203".into())));
    assert_eq!(batch::feature_rgb(&feature), Some([1, 2, 3]));
    feature
        .qualifiers
        .push(("color".into(), Some("#040506".into())));
    assert_eq!(batch::feature_rgb(&feature), Some([4, 5, 6]));
    feature
        .qualifiers
        .push(("COLOR".into(), Some("#070809".into())));
    assert_eq!(batch::feature_rgb(&feature), None);
}

#[test]
fn imported_directional_colors_follow_strand_and_preserve_mixed_ambiguity() {
    let mut feature = feature(Location::range(0, 3));
    feature
        .qualifiers
        .push(("ApEinfo_fwdcolor".into(), Some("#112233".into())));
    feature
        .qualifiers
        .push(("ApEinfo_revcolor".into(), Some("#445566".into())));
    assert_eq!(batch::feature_rgb(&feature), Some([17, 34, 51]));
    feature.location = Location::Complement(Box::new(feature.location));
    assert_eq!(batch::feature_rgb(&feature), Some([68, 85, 102]));
    feature.location = Location::Join(vec![
        Location::range(0, 3),
        Location::Complement(Box::new(Location::range(3, 6))),
    ]);
    assert_eq!(batch::feature_rgb(&feature), None);
    assert_eq!(
        batch::feature_rgb_for_strand(&feature, false),
        Some([17, 34, 51])
    );
    assert_eq!(
        batch::feature_rgb_for_strand(&feature, true),
        Some([68, 85, 102])
    );
    feature.location = Location::Unresolved("remote:1..3".into());
    assert_eq!(batch::feature_rgb(&feature), None);
    feature
        .qualifiers
        .push(("color".into(), Some("#aabbcc".into())));
    assert_eq!(batch::feature_rgb(&feature), Some([170, 187, 204]));
}

#[test]
fn no_op_preserves_document_and_request_roundtrips_without_unknown_fields() {
    let source = source();
    let r = request(&source);
    assert_eq!(batch::apply(&source, &r).unwrap(), source);
    let json = serde_json::to_string(&r).unwrap();
    assert_eq!(serde_json::from_str::<Request>(&json).unwrap(), r);
    let mut value = serde_json::to_value(&r).unwrap();
    value["verified"] = serde_json::json!(true);
    assert!(serde_json::from_value::<Request>(value).is_err());
}
