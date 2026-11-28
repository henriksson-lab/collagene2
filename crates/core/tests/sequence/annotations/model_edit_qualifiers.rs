use collagene_model::{Feature, Location, Molecule, Topology};

fn source() -> Molecule {
    let mut source = Molecule::dna("reference", "ACGTACGT", Topology::Linear).unwrap();
    source.features.push(Feature {
        kind: "CDS".into(),
        location: Location::range(0, 4),
        needs_review: false,
        qualifiers: vec![
            ("TRANSLATION".into(), Some("old assertion".into())),
            ("translation".into(), Some("another assertion".into())),
        ],
    });
    source.features.push(Feature {
        kind: "CDS".into(),
        location: Location::range(4, 8),
        needs_review: false,
        qualifiers: vec![("Translation".into(), Some("untouched assertion".into()))],
    });
    source
}

#[test]
fn native_reference_edits_remove_case_variants_only_from_touched_features() {
    let source = source();
    let before = source.revision();
    let edited = source.replace(1, 2, "T").unwrap();
    assert_eq!(edited.sequence.as_str(), "ATGTACGT");
    assert!(edited.features[0].needs_review);
    assert!(edited.features[0].qualifiers.is_empty());
    assert!(!edited.features[1].needs_review);
    assert_eq!(edited.features[1].qualifiers, source.features[1].qualifiers);
    assert_eq!(source.revision(), before);
}

#[test]
fn coordinate_bearing_case_variants_do_not_bypass_native_edit_rejection() {
    for key in ["TRANSL_EXCEPT", "AnTiCoDoN", "RPT_UNIT_RANGE", "Coded_By"] {
        let mut source = source();
        source.features[0]
            .qualifiers
            .push((key.into(), Some("source coordinates".into())));
        let before = source.revision();
        assert!(source.replace(1, 2, "T").is_err(), "{key}");
        assert_eq!(source.revision(), before);
    }
}
