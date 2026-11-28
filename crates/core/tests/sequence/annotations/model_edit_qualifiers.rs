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
    assert!(
        !edited.features[0]
            .qualifiers
            .iter()
            .any(|(key, _)| key.eq_ignore_ascii_case("translation"))
    );
    let notes = edited.features[0]
        .qualifiers
        .iter()
        .filter(|(key, _)| key == "note")
        .filter_map(|(_, value)| value.as_deref())
        .collect::<Vec<_>>();
    assert_eq!(notes.len(), 2);
    assert!(
        notes
            .iter()
            .any(|value| value.contains("/TRANSLATION=old assertion")
                && value.contains("source evidence"))
    );
    assert!(
        notes
            .iter()
            .any(|value| value.contains("/translation=another assertion")
                && value.contains("source evidence"))
    );
    assert!(!edited.features[1].needs_review);
    assert_eq!(edited.features[1].qualifiers, source.features[1].qualifiers);
    assert_eq!(source.revision(), before);
}

#[test]
fn coordinate_bearing_case_variants_obey_the_supported_and_unsupported_edit_profiles() {
    for key in ["TRANSL_EXCEPT", "AnTiCoDoN", "RPT_UNIT_RANGE", "Coded_By"] {
        let mut source = source();
        source.features[0]
            .qualifiers
            .push((key.into(), Some("source coordinates".into())));
        let before = source.revision();
        if key.eq_ignore_ascii_case("transl_except") || key.eq_ignore_ascii_case("anticodon") {
            let edited = source.replace(1, 2, "T").unwrap();
            assert!(edited.features[0].needs_review);
            assert!(
                !edited.features[0]
                    .qualifiers
                    .iter()
                    .any(|(name, _)| name.eq_ignore_ascii_case(key))
            );
            assert!(
                edited.features[0]
                    .qualifiers
                    .iter()
                    .any(|(name, value)| name == "note"
                        && value.as_deref().is_some_and(|text| text.contains(key)
                            && text.contains("source coordinates")
                            && text.contains("source evidence")))
            );
        } else {
            assert!(source.replace(1, 2, "T").is_err(), "{key}");
        }
        assert_eq!(source.revision(), before);
    }
}
