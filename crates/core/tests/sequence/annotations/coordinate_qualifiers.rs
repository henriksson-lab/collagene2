//! Literal source/product geometry for coordinate-bearing qualifier projection.
use collagene_core::{
    formats::sequence,
    sequence::{
        annotations::{self, MappingSegment},
        editing::set_feature,
    },
};
use collagene_model::{Feature, Location, Molecule, Topology};

fn feature(kind: &str, location: Location, qualifiers: &[(&str, Option<&str>)]) -> Feature {
    Feature {
        kind: kind.into(),
        location,
        needs_review: false,
        qualifiers: qualifiers
            .iter()
            .map(|(k, v)| ((*k).into(), v.map(str::to_owned)))
            .collect(),
    }
}

fn sec_source() -> Molecule {
    let mut molecule = Molecule::dna(
        "Literal Sec annotation",
        "CCCATGTGATAACCC",
        Topology::Linear,
    )
    .unwrap();
    molecule.features.push(feature(
        "CDS",
        Location::range(3, 12),
        &[
            ("label", Some("Literal CDS")),
            ("transl_except", Some("(pos:7..9,aa:Sec)")),
            ("translation", Some("MU")),
            ("codon_start", Some("1")),
            ("x_opaque", Some("first source value")),
            ("x_opaque", Some("second source value")),
        ],
    ));
    molecule
}

fn values<'a>(feature: &'a Feature, key: &str) -> Vec<&'a str> {
    feature
        .qualifiers
        .iter()
        .filter(|(k, _)| k.eq_ignore_ascii_case(key))
        .filter_map(|(_, v)| v.as_deref())
        .collect()
}

fn mapped(source: &Molecule, offset: usize) -> annotations::Projection {
    annotations::project_features(
        source,
        &[MappingSegment {
            source_start: 0,
            source_end: source.sequence.len(),
            destination_start: offset,
            reverse: false,
        }],
        source.sequence.len() + offset,
    )
    .unwrap()
}

fn roundtrip(product: &Molecule) -> Molecule {
    let text = sequence::export_genbank(product).unwrap();
    let mut parsed = sequence::import_genbank(&text).unwrap();
    assert_eq!(parsed.molecules.len(), 1);
    let parsed = parsed.molecules.remove(0);
    assert_eq!(parsed.sequence, product.sequence);
    let bytes = collagene_core::protocol::save_document(product).unwrap();
    assert_eq!(
        collagene_core::protocol::load_document(&bytes).unwrap(),
        *product
    );
    parsed
}

#[test]
fn cropping_actual_source_shifts_both_exception_coordinates_and_preserves_unknown_multiplicity() {
    let mut source = sec_source();
    source.features[0]
        .qualifiers
        .push(("transl_except".into(), Some("(pos:10..12,aa:TERM)".into())));
    let original = source.clone();
    let product = annotations::extract_region(&source, 3, 12).unwrap();
    assert_eq!(product.sequence.as_str(), "ATGTGATAA");
    assert_eq!(product.features[0].location, Location::range(0, 9));
    assert_eq!(
        values(&product.features[0], "transl_except"),
        vec!["(pos:4..6,aa:Sec)", "(pos:7..9,aa:TERM)"]
    );
    assert_eq!(values(&product.features[0], "translation"), vec!["MU"]);
    assert_eq!(values(&product.features[0], "codon_start"), vec!["1"]);
    assert_eq!(
        values(&product.features[0], "x_opaque"),
        vec!["first source value", "second source value"]
    );
    assert!(!product.features[0].needs_review);
    assert_eq!(roundtrip(&product).features, product.features);
    assert_eq!(source, original);
}

#[test]
fn reverse_complement_changes_qualifier_strand_without_changing_its_biological_codon() {
    let source = sec_source();
    let product = annotations::reverse_complement(&source).unwrap();
    assert_eq!(product.sequence.as_str(), "GGGTTATCACATGGG");
    assert_eq!(
        product.features[0].location,
        Location::Complement(Box::new(Location::range(3, 12)))
    );
    assert_eq!(
        product.features[0]
            .location
            .extract(&product.sequence)
            .unwrap(),
        "ATGTGATAA"
    );
    assert_eq!(
        values(&product.features[0], "transl_except"),
        vec!["(pos:complement(7..9),aa:Sec)"]
    );
    assert_eq!(roundtrip(&product).features, product.features);
    assert_eq!(
        annotations::reverse_complement(&product).unwrap().features,
        source.features
    );
}

#[test]
fn circular_rotation_keeps_a_split_codon_while_the_same_linear_cut_deactivates_its_claim() {
    let mut source = Molecule::dna("Circle", "ATGTGATAACCC", Topology::Circular).unwrap();
    source.features.push(feature(
        "CDS",
        Location::range(0, 9),
        &[
            ("transl_except", Some("(pos:4..6,aa:Sec)")),
            ("translation", Some("MU")),
            ("codon_start", Some("1")),
        ],
    ));
    let rotated = annotations::rotate_origin(&source, 5).unwrap();
    assert_eq!(rotated.sequence.as_str(), "ATAACCCATGTG");
    assert_eq!(
        rotated.features[0]
            .location
            .extract(&rotated.sequence)
            .unwrap(),
        "ATGTGATAA"
    );
    assert_eq!(
        values(&rotated.features[0], "transl_except"),
        vec!["(pos:join(11..12,1),aa:Sec)"]
    );
    assert!(!rotated.features[0].needs_review);
    assert_eq!(roundtrip(&rotated).features, rotated.features);
    let linear = annotations::extract_region(&source, 5, 17).unwrap();
    assert_eq!(linear.sequence, rotated.sequence);
    assert_eq!(linear.topology, Topology::Linear);
    assert!(linear.features[0].needs_review);
    for key in ["transl_except", "translation", "codon_start"] {
        assert!(values(&linear.features[0], key).is_empty());
    }
    assert!(
        values(&linear.features[0], "note")
            .iter()
            .any(|note| note.contains("(pos:4..6,aa:Sec)"))
    );
}

#[test]
fn a_spliced_codon_uses_biological_join_order_and_literal_original_coordinates() {
    let mut source = Molecule::dna("Spliced TGA", "TGCCAATA", Topology::Linear).unwrap();
    source.features.push(feature(
        "CDS",
        Location::Join(vec![Location::range(0, 2), Location::range(4, 8)]),
        &[
            ("transl_except", Some("(pos:join(1..2,5),aa:Sec)")),
            ("codon_start", Some("1")),
        ],
    ));
    assert_eq!(
        source.features[0]
            .location
            .extract(&source.sequence)
            .unwrap(),
        "TGAATA"
    );
    let projected = mapped(&source, 2);
    assert_eq!(
        values(&projected.features[0], "transl_except"),
        vec!["(pos:join(3..4,7),aa:Sec)"]
    );
    let mut product =
        Molecule::dna("Shifted spliced reference", "CCTGCCAATA", Topology::Linear).unwrap();
    product.features = projected.features;
    assert_eq!(
        product.features[0]
            .location
            .extract(&product.sequence)
            .unwrap(),
        "TGAATA"
    );
    assert_eq!(roundtrip(&product).features, product.features);
    for invalid in ["(pos:join(5,1..2),aa:Sec)", "(pos:join(1,3..4),aa:Sec)"] {
        source.features[0].qualifiers[0].1 = Some(invalid.into());
        assert!(
            annotations::project_features(
                &source,
                &[MappingSegment {
                    source_start: 0,
                    source_end: 8,
                    destination_start: 2,
                    reverse: false
                }],
                10
            )
            .is_err(),
            "{invalid}"
        );
    }
}

#[test]
fn codon_start_and_short_terminal_stop_exceptions_preserve_literal_phase_and_aa_payload() {
    let mut phased = Molecule::dna(
        "One incomplete leading base",
        "AATGTGATAA",
        Topology::Linear,
    )
    .unwrap();
    phased.features.push(feature(
        "CDS",
        Location::range(0, 10),
        &[
            ("codon_start", Some("2")),
            ("transl_except", Some("(pos:5..7,aa:Sec)")),
        ],
    ));
    assert_eq!(
        values(&mapped(&phased, 2).features[0], "transl_except"),
        vec!["(pos:7..9,aa:Sec)"]
    );
    assert_eq!(
        values(&mapped(&phased, 2).features[0], "codon_start"),
        vec!["2"]
    );
    for (dna, exception, expected) in [
        ("ATGT", "(pos:4,aa:TERM)", "(pos:6,aa:TERM)"),
        ("ATGTA", "(pos:4..5,aa:TERM)", "(pos:6..7,aa:TERM)"),
    ] {
        let mut source = Molecule::dna(
            "Polyadenylation-completable terminal stop",
            dna,
            Topology::Linear,
        )
        .unwrap();
        source.features.push(feature(
            "CDS",
            Location::range(0, dna.len()),
            &[
                ("transl_except", Some(exception)),
                ("codon_start", Some("1")),
            ],
        ));
        assert_eq!(
            values(&mapped(&source, 2).features[0], "transl_except"),
            vec![expected]
        );
    }
}

#[test]
fn anticodon_keeps_optional_sequence_and_other_payload_bytes_without_synthesizing_fields() {
    let mut source = sec_source();
    source.features[0] = feature(
        "tRNA",
        Location::range(3, 12),
        &[
            ("ANTICODON", Some("(aa:Ser,pos:7..9,seq:tga)")),
            ("anticodon", Some("(pos:7..9,aa:Ser)")),
        ],
    );
    let product = annotations::extract_region(&source, 3, 12).unwrap();
    assert_eq!(
        values(&product.features[0], "anticodon"),
        vec!["(aa:Ser,pos:4..6,seq:tga)", "(pos:4..6,aa:Ser)"]
    );
    assert_eq!(roundtrip(&product).features, product.features);
}

#[test]
fn clipped_or_interrupted_parent_removes_active_coordinate_and_translation_assertions() {
    let mut source = sec_source();
    source.features[0].qualifiers[1].1 = Some("(pos:NC_UNKNOWN:7..9,aa:Sec)".into());
    let before = source.clone();
    let product = annotations::extract_region(&source, 4, 12).unwrap();
    assert_eq!(product.sequence.as_str(), "TGTGATAA");
    assert!(product.features[0].needs_review);
    for key in ["transl_except", "translation", "codon_start"] {
        assert!(values(&product.features[0], key).is_empty());
    }
    assert!(
        values(&product.features[0], "note")
            .iter()
            .any(|note| note.contains("NC_UNKNOWN:7..9"))
    );
    assert_eq!(
        values(&product.features[0], "x_opaque"),
        vec!["first source value", "second source value"]
    );
    assert_eq!(source, before);
}

#[test]
fn malformed_or_unsupported_active_wrappers_do_not_guess_a_mapped_coordinate() {
    for invalid in [
        "(pos:7..9,aa:Sec)junk",
        "(pos:7..9,pos:7..9,aa:Sec)",
        "(pos:7..9,aa:Sec,extra:1)",
        "(pos:7..9)",
        "(aa:Sec)",
        "(pos:NC_000001.1:7..9,aa:Sec)",
        "(pos:<7..9,aa:Sec)",
        "(pos:order(7,8,9),aa:Sec)",
        "(pos:complement(7..9),aa:Sec)",
        "(pos:1..3,aa:Sec)",
        "(pos:5..7,aa:Sec)",
        "(pos:9,aa:Sec)",
    ] {
        let mut source = sec_source();
        source.features[0].qualifiers[1].1 = Some(invalid.into());
        let before = source.clone();
        assert!(
            annotations::extract_region(&source, 3, 12).is_err(),
            "{invalid}"
        );
        assert_eq!(source, before);
    }
    let mut flag = sec_source();
    flag.features[0].qualifiers[1].1 = None;
    assert!(annotations::extract_region(&flag, 3, 12).is_err());
}

#[test]
fn identity_keeps_source_bytes_and_dropped_annotations_remain_unmapped_evidence() {
    let mut source = sec_source();
    source.features[0]
        .qualifiers
        .push(("coded_by".into(), Some("NC_REMOTE:1..3".into())));
    source.features[0].qualifiers[1].1 = Some("(pos:remote:7..9,aa:Sec)unknown".into());
    let identity = annotations::project_features(
        &source,
        &[MappingSegment {
            source_start: 0,
            source_end: 15,
            destination_start: 0,
            reverse: false,
        }],
        15,
    )
    .unwrap();
    assert_eq!(identity.features, source.features);
    let absent = annotations::project_features(
        &source,
        &[MappingSegment {
            source_start: 0,
            source_end: 3,
            destination_start: 0,
            reverse: false,
        }],
        3,
    )
    .unwrap();
    assert!(absent.features.is_empty());
    assert_eq!(absent.dropped_features, 1);
    assert!(absent.unmapped_features.is_empty());
    assert_eq!(source.features, identity.features);
    for key in ["rpt_unit_range", "coded_by"] {
        let mut unsupported = sec_source();
        unsupported.features[0]
            .qualifiers
            .push((key.into(), Some("7..9".into())));
        assert!(annotations::extract_region(&unsupported, 3, 12).is_err());
    }
}

#[test]
fn wrapper_and_shared_projection_work_bounds_refuse_without_mutating_the_sec_source() {
    let mut oversized = sec_source();
    oversized.features[0].qualifiers[1].1 = Some("x".repeat(65_537));
    assert!(matches!(
        annotations::extract_region(&oversized, 3, 12),
        Err(collagene_model::Error::TooLarge)
    ));
    // Even nesting preserves the original strand and is below the wrapper-depth
    // ceiling. The aggregate structural work exceeds the shared budget.
    let mut aggregate = sec_source();
    let position = format!("{}7..9{}", "complement(".repeat(20), ")".repeat(20));
    let value = format!("(pos:{position},aa:Sec)");
    aggregate.features[0].qualifiers = vec![("transl_except".into(), Some(value)); 10_000];
    let before = aggregate.revision();
    let error = annotations::extract_region(&aggregate, 3, 12).unwrap_err();
    assert!(error.to_string().contains("budget"), "{error}");
    assert_eq!(aggregate.revision(), before);
}

#[test]
fn manual_location_change_cannot_hide_coordinate_fields_by_removing_or_adding_them() {
    let source = sec_source();
    let before = source.clone();
    let mut replacement = source.features[0].clone();
    replacement.location = Location::range(4, 12);
    replacement
        .qualifiers
        .retain(|(key, _)| !key.eq_ignore_ascii_case("transl_except"));
    assert!(set_feature(&source, Some(0), replacement).is_err());
    let mut plain = source.clone();
    plain.features[0]
        .qualifiers
        .retain(|(key, _)| !key.eq_ignore_ascii_case("transl_except"));
    for key in ["TRANSL_EXCEPT", "AnTiCoDoN", "RPT_UNIT_RANGE", "Coded_By"] {
        let mut added = plain.features[0].clone();
        added.location = Location::range(4, 12);
        added
            .qualifiers
            .push((key.into(), Some("source-coordinate claim".into())));
        assert!(set_feature(&plain, Some(0), added).is_err(), "{key}");
    }
    assert_eq!(source, before);
}

#[test]
fn other_manual_location_changes_deactivate_case_variants_but_preserve_unrelated_qualifiers() {
    let mut source = sec_source();
    source.features[0]
        .qualifiers
        .retain(|(key, _)| !key.eq_ignore_ascii_case("transl_except"));
    source.features[0]
        .qualifiers
        .push(("TRANSLATION".into(), Some("second stale assertion".into())));
    source.features[0]
        .qualifiers
        .push(("CoDoN_sTaRt".into(), Some("2".into())));
    let before = source.clone();
    let mut replacement = source.features[0].clone();
    replacement.location = Location::range(4, 12);
    let edited = set_feature(&source, Some(0), replacement).unwrap();
    assert!(edited.features[0].needs_review);
    for key in ["translation", "codon_start"] {
        assert!(values(&edited.features[0], key).is_empty());
    }
    assert!(
        values(&edited.features[0], "note")
            .iter()
            .any(|note| note.contains("second stale assertion"))
    );
    assert_eq!(
        values(&edited.features[0], "x_opaque"),
        vec!["first source value", "second source value"]
    );
    let mut cosmetic = source.features[0].clone();
    cosmetic
        .qualifiers
        .push(("note".into(), Some("Only a label/documentary edit".into())));
    let same_location = set_feature(&source, Some(0), cosmetic.clone()).unwrap();
    assert_eq!(same_location.features[0], cosmetic);
    assert_eq!(source, before);
}

#[test]
fn exact_pcr_tails_move_the_exception_into_the_actual_amplicon_and_genbank_record() {
    use collagene_core::{chemistry::Strand, sequence::analysis::Primer, verification::pcr};
    let mut source =
        Molecule::dna("PCR source", "AGTCCGATATGTGATAACGATACGA", Topology::Linear).unwrap();
    source.features.push(feature(
        "CDS",
        Location::range(8, 17),
        &[
            ("transl_except", Some("(pos:12..14,aa:Sec)")),
            ("codon_start", Some("1")),
        ],
    ));
    let before = source.clone();
    let products = pcr::simulate_pcr_documents(
        &source,
        &Primer::new("F", "AGTCCGAT", "GG").unwrap(),
        &Primer::new("R", "TCGTATCG", "CCC").unwrap(),
        Strand::Top,
    )
    .unwrap();
    assert_eq!(products.len(), 1);
    let product = &products[0].molecule;
    assert_eq!(product.sequence.as_str(), "GGAGTCCGATATGTGATAACGATACGAGGG");
    assert_eq!(product.features[0].location, Location::range(10, 19));
    assert_eq!(
        values(&product.features[0], "transl_except"),
        vec!["(pos:14..16,aa:Sec)"]
    );
    assert_eq!(
        product.features[0]
            .location
            .extract(&product.sequence)
            .unwrap(),
        "ATGTGATAA"
    );
    assert_eq!(roundtrip(product).features, product.features);
    assert_eq!(source, before);
}

#[test]
fn imported_genbank_source_stays_exact_while_the_product_retains_documentary_original() {
    let raw = "LOCUS       QUAL_LITERAL              15 bp    DNA     linear   UNK 28-NOV-2026\nFEATURES             Location/Qualifiers\n     CDS             4..12\n                     /transl_except=(pos:7..9,aa:Sec)\n                     /db_xref=\"Unresolved:opaque\"\nORIGIN\n        1 cccatgtgat aaccc\n//\n";
    let mut imported = sequence::import_genbank(raw).unwrap();
    let source = imported.molecules.remove(0);
    assert_eq!(source.original_genbank.as_deref(), Some(raw));
    let before = source.clone();
    let product = annotations::extract_region(&source, 3, 12).unwrap();
    assert_eq!(product.sequence.as_str(), "ATGTGATAA");
    assert_eq!(
        values(&product.features[0], "transl_except"),
        vec!["(pos:4..6,aa:Sec)"]
    );
    assert_eq!(
        values(&product.features[0], "db_xref"),
        vec!["Unresolved:opaque"]
    );
    assert_eq!(
        product
            .metadata
            .get("annotation_mapping:original_genbank")
            .map(String::as_str),
        Some(raw)
    );
    assert!(product.original_genbank.is_none());
    assert_eq!(roundtrip(&product).features, product.features);
    assert_eq!(source, before);
}

#[test]
fn uncertain_parent_deactivates_coordinate_claims_before_interpreting_their_payload() {
    for location in [
        Location::Range {
            start: 3,
            end: 12,
            before: true,
            after: false,
        },
        Location::Order(vec![Location::range(3, 6), Location::range(6, 12)]),
    ] {
        let mut source = sec_source();
        source.features[0].location = location;
        source.features[0].qualifiers[1].1 =
            Some("(pos:REMOTE:7..9,aa:Sec)unparsed source evidence".into());
        source.features[0].qualifiers.push((
            "anticodon".into(),
            Some("unparsed anticodon evidence".into()),
        ));
        let before = source.clone();
        let product = annotations::extract_region(&source, 3, 12).unwrap();
        assert!(product.features[0].needs_review);
        for key in ["transl_except", "anticodon", "translation", "codon_start"] {
            assert!(values(&product.features[0], key).is_empty(), "{key}");
        }
        let notes = values(&product.features[0], "note").join("\n");
        assert!(notes.contains("REMOTE:7..9"));
        assert!(notes.contains("unparsed anticodon evidence"));
        assert_eq!(
            values(&product.features[0], "x_opaque"),
            vec!["first source value", "second source value"]
        );
        assert_eq!(source, before);
    }
}

#[test]
fn replacing_at_the_feature_limit_does_not_count_the_old_annotation_twice() {
    let mut source = Molecule::dna("At the feature limit", "ACGTACGT", Topology::Linear).unwrap();
    source.features = vec![
        feature(
            "misc_feature",
            Location::range(0, 1),
            &[("label", Some("Original"))]
        );
        4096
    ];
    let mut replacement = source.features[0].clone();
    replacement.qualifiers[0].1 = Some("Updated annotation".into());
    let updated = set_feature(&source, Some(0), replacement).unwrap();
    assert_eq!(updated.features.len(), 4096);
    assert_eq!(updated.features[0].label(), "Updated annotation");
    assert_eq!(source.features[0].label(), "Original");
}
