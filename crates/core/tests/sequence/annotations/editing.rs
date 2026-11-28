use collagene_core::sequence::editing::*;
use collagene_model::annotations::insdc::{location_text, parse_location};
use collagene_model::{Error, Feature, Location, Molecule, Sequence, Topology};
fn molecule() -> Molecule {
    Molecule::dna("fixture", "AAAACCCCGGGGTTTT", Topology::Circular).unwrap()
}
#[test]
fn annotation_form_accepts_compound_fuzzy_and_duplicate_qualifiers_without_flattening() {
    let source = molecule();
    let location = parse_location("complement(join(13..16,1..4))", 16).unwrap();
    assert_eq!(location.extract(&source.sequence).unwrap(), "TTTTAAAA");
    assert_eq!(
        location_text(&location).unwrap(),
        "complement(join(13..16,1..4))"
    );
    let feature = Feature {
        kind: "CDS".into(),
        location,
        qualifiers: vec![
            ("label".into(), Some("compound".into())),
            ("note".into(), Some("first\nsecond".into())),
            ("note".into(), Some("third".into())),
            ("pseudo".into(), None),
        ],
        needs_review: false,
    };
    let edited = set_feature(&source, None, feature.clone()).unwrap();
    assert_eq!(edited.features[0], feature);
    let binary = collagene_core::protocol::save_document(&edited).unwrap();
    assert_eq!(
        collagene_core::protocol::load_document(&binary).unwrap(),
        edited
    );
    let gb = collagene_core::formats::sequence::export_genbank(&edited).unwrap();
    let restored = collagene_core::formats::sequence::import_genbank(&gb).unwrap();
    assert_eq!(restored.molecules[0].features[0].location, feature.location);
    let fuzzy = parse_location("<2..>8", 16).unwrap();
    assert!(fuzzy.extract(&source.sequence).is_err());
    assert!(parse_location("0..8", 16).is_err());
    assert!(parse_location("1..17", 16).is_err());
    assert!(parse_location("1..8\nORIGIN", 16).is_err());
    let nested = format!("{}1..4{}", "complement(".repeat(33), ")".repeat(33));
    assert!(parse_location(&nested, 16).is_err());
}
#[test]
fn location_changes_clear_stale_translation_and_edits_are_reversible_snapshots() {
    let mut source = molecule();
    source.features.push(Feature {
        kind: "CDS".into(),
        location: Location::range(0, 9),
        qualifiers: vec![("translation".into(), Some("old".into()))],
        needs_review: false,
    });
    let mut replacement = source.features[0].clone();
    replacement.location = Location::range(3, 12);
    let changed = set_feature(&source, Some(0), replacement).unwrap();
    assert!(changed.features[0].needs_review);
    assert!(
        changed.features[0]
            .qualifiers
            .iter()
            .all(|(key, _)| !key.eq_ignore_ascii_case("translation")
                && !key.eq_ignore_ascii_case("codon_start"))
    );
    assert!(
        changed.features[0]
            .qualifiers
            .iter()
            .any(|(key, value)| key == "note"
                && value
                    .as_deref()
                    .is_some_and(|text| text.contains("/translation=old")
                        && text.contains("source evidence")))
    );
    assert_eq!(source.features[0].qualifiers[0].1.as_deref(), Some("old"));
    assert!(delete_feature(&source, 1).is_err());
    assert!(delete_feature(&source, 0).unwrap().features.is_empty());
    assert_eq!(rename(&source, " new name ").unwrap().name, "new name");
    assert!(rename(&source, "a\nb").is_err());
    let mut bad = source.features[0].clone();
    bad.kind = "CDS\nORIGIN".into();
    assert!(set_feature(&source, None, bad).is_err());
}
#[test]
fn extraction_is_bounded_before_materializing_repeated_valid_ranges() {
    let seq = Sequence::dna(&"A".repeat(1000)).unwrap();
    let repeated = Location::Join(vec![Location::range(0, 1000); 10_001]);
    repeated.validate(seq.len()).unwrap();
    assert!(matches!(repeated.extract(&seq), Err(Error::TooLarge)));
    let nested = Location::Complement(Box::new(Location::Join(vec![
        Location::Complement(Box::new(Location::range(0, 2))),
        Location::range(2, 4),
    ])));
    assert_eq!(
        nested.extract(&Sequence::dna("ACGT").unwrap()).unwrap(),
        "ACAC"
    );
}
