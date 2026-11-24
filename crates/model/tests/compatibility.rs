//! Literal wire/revision fixtures captured from the pre-extraction core model.
//! Expected JSON/Postcard/revisions are committed bytes, never regenerated here.
use collagene_model::*;
fn rich() -> Molecule {
    let mut m = Molecule::dna("golden-rich", "AAAACCCCGGGGTTTT", Topology::Circular).unwrap();
    m.features.push(Feature {
        kind: "CDS".into(),
        location: Location::Complement(Box::new(Location::Join(vec![
            Location::range(12, 16),
            Location::range(0, 4),
        ]))),
        qualifiers: vec![
            ("label".into(), Some("origin".into())),
            ("pseudo".into(), None),
        ],
        needs_review: false,
    });
    m.metadata.insert("zeta".into(), "last".into());
    m.metadata.insert("alpha".into(), "first".into());
    m.original_genbank = Some("LOCUS retained evidence\n//\n".into());
    m
}
#[test]
fn literal_simple_json_postcard_and_revision_stay_identical() {
    let m = Molecule::dna("golden", "a c\nGt", Topology::Linear).unwrap();
    assert_eq!(
        serde_json::to_vec(&m).unwrap(),
        include_bytes!("../fixtures/simple.json")
    );
    assert_eq!(
        postcard::to_allocvec(&m).unwrap(),
        include_bytes!("../fixtures/simple.postcard")
    );
    assert_eq!(
        m.revision(),
        "76dcf1814413551f814f42095931436c0f5397e8392d8843609a9ff9ffb80f5c"
    );
    let json: Molecule = serde_json::from_slice(include_bytes!("../fixtures/simple.json")).unwrap();
    let wire: Molecule =
        postcard::from_bytes(include_bytes!("../fixtures/simple.postcard")).unwrap();
    assert_eq!(json, m);
    assert_eq!(wire, m);
    m.validate().unwrap();
}
#[test]
fn literal_compound_features_sorted_metadata_and_original_evidence_keep_the_same_revision() {
    let m = rich();
    assert_eq!(
        serde_json::to_vec(&m).unwrap(),
        include_bytes!("../fixtures/rich.json")
    );
    assert_eq!(
        postcard::to_allocvec(&m).unwrap(),
        include_bytes!("../fixtures/rich.postcard")
    );
    assert_eq!(
        m.revision(),
        "ecd1e2bc9497600cfee7ccbbbe42a779c4e3dbc70ec6bae9c54f0306bee63ed1"
    );
    let wire: Molecule = postcard::from_bytes(include_bytes!("../fixtures/rich.postcard")).unwrap();
    assert_eq!(wire, m);
    assert_eq!(
        wire.features[0].location.extract(&wire.sequence).unwrap(),
        "TTTTAAAA"
    );
    assert_eq!(wire.features[0].label(), "origin");
    assert_eq!(
        wire.original_genbank.as_deref(),
        Some("LOCUS retained evidence\n//\n")
    );
}
#[test]
fn literal_location_variants_and_fuzzy_semantics_are_unchanged() {
    let l = Location::Order(vec![
        Location::Range {
            start: 1,
            end: 3,
            before: true,
            after: false,
        },
        Location::Unresolved("REMOTE.1:2..4".into()),
    ]);
    assert_eq!(
        serde_json::to_vec(&l).unwrap(),
        include_bytes!("../fixtures/location.json")
    );
    assert_eq!(
        postcard::to_allocvec(&l).unwrap(),
        include_bytes!("../fixtures/location.postcard")
    );
    let decoded: Location =
        postcard::from_bytes(include_bytes!("../fixtures/location.postcard")).unwrap();
    assert_eq!(decoded, l);
    decoded.validate(8).unwrap();
    assert!(
        decoded
            .extract(&Sequence::dna("ACGTACGT").unwrap())
            .is_err()
    );
    let fuzzy = Location::Range {
        start: 1,
        end: 3,
        before: true,
        after: false,
    };
    assert!(fuzzy.extract(&Sequence::dna("ACGT").unwrap()).is_err());
}
#[test]
fn alphabet_normalization_complements_and_invalid_calls_preserve_existing_policy() {
    let dna = Sequence::new(Alphabet::Dna, "a c\nGtRYSWKMBDHVN").unwrap();
    assert_eq!(dna.as_str(), "ACGTRYSWKMBDHVN");
    assert_eq!(
        dna.reverse_complement().unwrap().as_str(),
        "NBDHVKMWSRYACGT"
    );
    assert!(Sequence::dna("ACGU").is_err());
    assert!(Sequence::new(Alphabet::Rna, "ACGT").is_err());
    assert_eq!(
        Sequence::new(Alphabet::Rna, "a cgu")
            .unwrap()
            .reverse_complement()
            .unwrap()
            .as_str(),
        "ACGU"
    );
    let protein = Sequence::new(Alphabet::Protein, "m x z *").unwrap();
    assert_eq!(protein.as_str(), "MXZ*");
    assert!(protein.reverse_complement().is_err());
    assert!(Sequence::new(Alphabet::Protein, "M-*").is_err());
    assert!(matches!(Sequence::dna(" \n\t"), Err(Error::EmptySequence)));
    assert!(matches!(
        Sequence::dna("Aé"),
        Err(Error::InvalidBase {
            alphabet: Alphabet::Dna,
            character: 'é',
            position: 2
        })
    ));
    assert!(matches!(
        Sequence::dna(" a ?"),
        Err(Error::InvalidBase {
            character: '?',
            position: 4,
            ..
        })
    ));
}
#[test]
fn deserialization_still_validates_alphabet_and_caps_recursive_location_depth() {
    for bad in [
        r#"{"alphabet":"Dna","bases":"U"}"#,
        r#"{"alphabet":"Rna","bases":"T"}"#,
        r#"{"alphabet":"Protein","bases":"-"}"#,
        r#"{"alphabet":"Dna","bases":""}"#,
        r#"{"alphabet":"Dna","bases":"ACGT","unexpected":true}"#,
    ] {
        assert!(serde_json::from_str::<Sequence>(bad).is_err());
    }
    let mut deep = r#"{"Range":{"start":0,"end":1,"before":false,"after":false}}"#.to_string();
    for _ in 0..32 {
        deep = format!("{{\"Complement\":{deep}}}");
    }
    assert!(serde_json::from_str::<Location>(&deep).is_err());
    let simple: Location =
        serde_json::from_str(r#"{"Range":{"start":0,"end":1,"before":false,"after":false}}"#)
            .unwrap();
    simple.validate(1).unwrap();
}
#[test]
fn invalid_coordinates_and_document_constraints_fail_without_invented_locations() {
    for l in [
        Location::range(1, 1),
        Location::range(3, 2),
        Location::range(0, 5),
        Location::Join(vec![]),
        Location::Order(vec![]),
    ] {
        assert!(l.validate(4).is_err());
    }
    let mut m = rich();
    m.features[0].location = Location::range(0, 17);
    assert!(m.validate().is_err());
    let mut protein = Molecule::dna("protein", "ACGT", Topology::Circular).unwrap();
    protein.sequence = Sequence::new(Alphabet::Protein, "MKW").unwrap();
    assert!(protein.validate().is_err());
    let mut unnamed = Molecule::dna(" ", "ACGT", Topology::Linear).unwrap();
    assert!(unnamed.validate().is_err());
    unnamed.name = "named".into();
    unnamed.validate().unwrap();
}
#[test]
fn exact_edit_behavior_keeps_the_parent_evidence_and_review_contract() {
    let m = rich();
    let changed = m.replace(1, 3, "GG").unwrap();
    assert_eq!(m.sequence.as_str(), "AAAACCCCGGGGTTTT");
    assert_eq!(changed.sequence.as_str(), "AGGACCCCGGGGTTTT");
    assert!(changed.features[0].needs_review);
    let mut guarded = m.clone();
    guarded.features[0]
        .qualifiers
        .push(("transl_except".into(), Some("(pos:1..3,aa:Sec)".into())));
    assert!(guarded.replace(0, 1, "G").is_err());
}
