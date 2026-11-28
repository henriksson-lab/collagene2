//! Literal native edit, clipboard and strand-curation expectations. Source bases
//! and biological exon order are authored here, independently of projection code.
use collagene_core::sequence::{
    annotations::batch::{self, ColorEdit},
    editing::set_feature,
    reference_transfer::{self as transfer, Part, RequestCopy, RequestPaste, Selection},
};
use collagene_model::{Error, Feature, Location, Molecule, Topology};

fn feature(location: Location, qualifiers: &[(&str, Option<&str>)]) -> Feature {
    Feature {
        kind: "CDS".into(),
        location,
        qualifiers: qualifiers
            .iter()
            .map(|(k, v)| ((*k).into(), v.map(str::to_owned)))
            .collect(),
        needs_review: false,
    }
}
fn source() -> Molecule {
    let mut source =
        Molecule::dna("Literal native edit", "CCCATGTGATAACCC", Topology::Linear).unwrap();
    source.features.push(feature(
        Location::range(3, 12),
        &[
            ("label", Some("Same biological CDS")),
            ("transl_except", Some("(pos:7..9,aa:Sec)")),
            ("anticodon", Some("(aa:Ser,pos:7..9,seq:tga)")),
            ("translation", Some("MU")),
            ("codon_start", Some("1")),
            ("x_opaque", Some("first")),
            ("x_opaque", Some("second")),
        ],
    ));
    source
}
fn values<'a>(feature: &'a Feature, key: &str) -> Vec<&'a str> {
    feature
        .qualifiers
        .iter()
        .filter(|(k, _)| k.eq_ignore_ascii_case(key))
        .filter_map(|(_, v)| v.as_deref())
        .collect()
}
fn inactive(feature: &Feature, original: &[&str]) {
    assert!(feature.needs_review);
    for key in ["transl_except", "anticodon", "translation", "codon_start"] {
        assert!(values(feature, key).is_empty(), "stale active /{key}");
    }
    for raw in original {
        assert!(
            values(feature, "note").iter().any(|n| n.contains(raw)),
            "missing source evidence {raw}"
        );
    }
}
fn roundtrip(product: &Molecule) {
    let native = collagene_core::protocol::save_document(product).unwrap();
    assert_eq!(
        collagene_core::protocol::load_document(&native).unwrap(),
        *product
    );
    let gb = collagene_core::formats::sequence::export_genbank(product).unwrap();
    let again = collagene_core::formats::sequence::import_genbank(&gb)
        .unwrap()
        .molecules
        .remove(0);
    assert_eq!(again.sequence, product.sequence);
    assert_eq!(again.features, product.features);
}
#[test]
fn native_edits_before_and_after_the_parent_remap_claims_without_changing_biology() {
    let source = source();
    let before = source.clone();
    for (start, end, replacement, dna, location, exception, anticodon) in [
        (
            0,
            0,
            "GG",
            "GGCCCATGTGATAACCC",
            Location::range(5, 14),
            "(pos:9..11,aa:Sec)",
            "(aa:Ser,pos:9..11,seq:tga)",
        ),
        (
            0,
            3,
            "",
            "ATGTGATAACCC",
            Location::range(0, 9),
            "(pos:4..6,aa:Sec)",
            "(aa:Ser,pos:4..6,seq:tga)",
        ),
        (
            15,
            15,
            "GG",
            "CCCATGTGATAACCCGG",
            Location::range(3, 12),
            "(pos:7..9,aa:Sec)",
            "(aa:Ser,pos:7..9,seq:tga)",
        ),
    ] {
        let product = source.replace(start, end, replacement).unwrap();
        let f = &product.features[0];
        assert_eq!(product.sequence.as_str(), dna);
        assert_eq!(f.location, location);
        assert_eq!(f.location.extract(&product.sequence).unwrap(), "ATGTGATAA");
        assert_eq!(values(f, "transl_except"), vec![exception]);
        assert_eq!(values(f, "anticodon"), vec![anticodon]);
        assert_eq!(values(f, "translation"), vec!["MU"]);
        assert_eq!(values(f, "codon_start"), vec!["1"]);
        assert_eq!(values(f, "x_opaque"), vec!["first", "second"]);
        assert!(!f.needs_review);
        roundtrip(&product);
    }
    assert_eq!(source, before);
}
#[test]
fn reverse_parent_edit_shifts_complemented_positions_without_reversing_their_payload() {
    let mut source = Molecule::dna(
        "Reverse biological CDS",
        "GGGTTATCACATGGG",
        Topology::Linear,
    )
    .unwrap();
    source.features.push(feature(
        Location::Complement(Box::new(Location::range(3, 12))),
        &[
            ("transl_except", Some("(pos:complement(7..9),aa:Sec)")),
            ("anticodon", Some("(pos:complement(7..9),aa:Ser,seq:tga)")),
            ("codon_start", Some("1")),
        ],
    ));
    let product = source.replace(0, 0, "CC").unwrap();
    let f = &product.features[0];
    assert_eq!(
        f.location,
        Location::Complement(Box::new(Location::range(5, 14)))
    );
    assert_eq!(f.location.extract(&product.sequence).unwrap(), "ATGTGATAA");
    assert_eq!(
        values(f, "transl_except"),
        vec!["(pos:complement(9..11),aa:Sec)"]
    );
    assert_eq!(
        values(f, "anticodon"),
        vec!["(pos:complement(9..11),aa:Ser,seq:tga)"]
    );
    assert!(!f.needs_review);
    roundtrip(&product);
}
#[test]
fn insertion_in_a_spliced_intronic_gap_preserves_the_same_ordered_exception_codon() {
    let mut source =
        Molecule::dna("Spliced literal codon", "CCTGCCAATA", Topology::Linear).unwrap();
    source.features.push(feature(
        Location::Join(vec![Location::range(2, 4), Location::range(6, 10)]),
        &[
            ("transl_except", Some("(pos:join(3..4,7),aa:Sec)")),
            ("codon_start", Some("1")),
        ],
    ));
    let product = source.replace(4, 4, "GG").unwrap();
    let f = &product.features[0];
    assert_eq!(product.sequence.as_str(), "CCTGGGCCAATA");
    assert_eq!(
        f.location,
        Location::Join(vec![Location::range(2, 4), Location::range(8, 12)])
    );
    assert_eq!(f.location.extract(&product.sequence).unwrap(), "TGAATA");
    assert_eq!(
        values(f, "transl_except"),
        vec!["(pos:join(3..4,9),aa:Sec)"]
    );
    assert_eq!(values(f, "codon_start"), vec!["1"]);
    assert!(!f.needs_review);
}
#[test]
fn phase_and_partial_terminal_stop_assertions_survive_edits_outside_the_parent() {
    for (dna, location, phase, exception, expected) in [
        (
            "CAATGTGATAAC",
            Location::range(1, 11),
            "2",
            "(pos:6..8,aa:Sec)",
            "(pos:8..10,aa:Sec)",
        ),
        (
            "CCATGTCC",
            Location::range(2, 6),
            "1",
            "(pos:6,aa:TERM)",
            "(pos:8,aa:TERM)",
        ),
        (
            "CCATGTACC",
            Location::range(2, 7),
            "1",
            "(pos:6..7,aa:TERM)",
            "(pos:8..9,aa:TERM)",
        ),
    ] {
        let mut source = Molecule::dna("Literal phase and TERM", dna, Topology::Linear).unwrap();
        source.features.push(feature(
            location,
            &[
                ("codon_start", Some(phase)),
                ("transl_except", Some(exception)),
            ],
        ));
        let product = source.replace(0, 0, "GG").unwrap();
        let f = &product.features[0];
        assert_eq!(values(f, "codon_start"), vec![phase]);
        assert_eq!(values(f, "transl_except"), vec![expected]);
        assert!(!f.needs_review);
    }
}
#[test]
fn boundary_insertions_exclude_new_bases_while_interior_edits_expand_the_native_envelope() {
    let source = source();
    let before = source.replace(3, 3, "GG").unwrap();
    assert_eq!(before.features[0].location, Location::range(5, 14));
    assert_eq!(
        values(&before.features[0], "transl_except"),
        vec!["(pos:9..11,aa:Sec)"]
    );
    assert!(!before.features[0].needs_review);
    let after = source.replace(12, 12, "GG").unwrap();
    assert_eq!(after.features[0].location, Location::range(3, 12));
    assert!(!after.features[0].needs_review);
    for (start, end, replacement, new_end) in [(6, 6, "A", 13), (4, 5, "C", 12), (4, 5, "", 11)] {
        let product = source.replace(start, end, replacement).unwrap();
        let f = &product.features[0];
        assert_eq!(f.location, Location::range(3, new_end));
        inactive(f, &["(pos:7..9,aa:Sec)", "(aa:Ser,pos:7..9,seq:tga)", "MU"]);
        assert_eq!(values(f, "x_opaque"), vec!["first", "second"]);
    }
}
#[test]
fn normalized_identity_preserves_raw_claims_and_review_flags_even_for_unparsed_source_text() {
    let mut source = source();
    source.features[0].qualifiers[1].1 =
        Some("(pos:remote:7..9,aa:Sec) trailing source text".into());
    source.features[0].qualifiers.push((
        "coded_by".into(),
        Some("unparsed original source evidence".into()),
    ));
    assert_eq!(source.replace(3, 6, " a t g \n").unwrap(), source);
    assert_eq!(source.replace(0, 0, "").unwrap(), source);
}
#[test]
fn retained_invalid_or_unsupported_claims_refuse_remapping_but_wholly_deleted_features_do_not() {
    for (key, value) in [
        ("transl_except", Some("(pos:7..9,aa:Sec)junk")),
        ("transl_except", Some("(pos:NC_REMOTE:7..9,aa:Sec)")),
        ("transl_except", Some("(pos:5..7,aa:Sec)")),
        ("anticodon", None),
        ("RPT_UNIT_RANGE", Some("7..9")),
        ("CoDeD_By", Some("NC_REMOTE:1..9")),
    ] {
        let mut source = source();
        source.features[0].qualifiers = vec![(key.into(), value.map(str::to_owned))];
        let before = source.clone();
        assert!(source.replace(0, 0, "GG").is_err(), "{key} {value:?}");
        assert_eq!(source, before);
        let deleted = source.replace(3, 12, "").unwrap();
        assert_eq!(deleted.sequence.as_str(), "CCCCCC");
        assert!(deleted.features.is_empty());
    }
}
#[test]
fn interrupted_supported_claims_become_evidence_without_parsing_stale_source_wrappers() {
    let mut source = source();
    source.features[0].qualifiers[1].1 = Some("(pos:NC_REMOTE:7..9,aa:Sec) unparsed source".into());
    source.features[0].qualifiers[2].1 = None;
    let before = source.clone();
    let product = source.replace(4, 5, "C").unwrap();
    inactive(
        &product.features[0],
        &["(pos:NC_REMOTE:7..9,aa:Sec) unparsed source", "MU"],
    );
    assert!(
        values(&product.features[0], "note")
            .iter()
            .any(|note| note.contains("anticodon") && note.contains("(flag)"))
    );
    assert_eq!(source, before);
}

#[test]
fn native_edits_preserve_rich_annotations_beyond_the_common_consumer_feature_limit() {
    let mut source = source();
    source.features = vec![source.features[0].clone(); 4_097];
    let original = source.clone();
    let product = source.replace(0, 0, "GG").unwrap();
    assert_eq!(product.sequence.as_str(), "GGCCCATGTGATAACCC");
    assert_eq!(product.features.len(), 4_097);
    for feature in &product.features {
        assert_eq!(feature.location, Location::range(5, 14));
        assert_eq!(values(feature, "transl_except"), vec!["(pos:9..11,aa:Sec)"]);
        assert_eq!(
            values(feature, "anticodon"),
            vec!["(aa:Ser,pos:9..11,seq:tga)"]
        );
        assert_eq!(values(feature, "translation"), vec!["MU"]);
        assert_eq!(values(feature, "x_opaque"), vec!["first", "second"]);
        assert!(!feature.needs_review);
    }
    assert_eq!(source, original);
}
#[test]
fn native_edits_retain_original_source_attachments_larger_than_four_mebibytes() {
    let mut source = source();
    source.original_genbank = Some(format!("COMMENT     {}\n", "x".repeat(5 * 1024 * 1024)));
    let original_annotations = source.features.clone();
    let original_raw = source.original_genbank.as_deref().unwrap();
    assert!(original_raw.len() > 4 * 1024 * 1024);
    assert!(!original_raw.contains("REFERENCE"));
    let product = source.replace(0, 0, "GG").unwrap();
    assert_eq!(product.sequence.as_str(), "GGCCCATGTGATAACCC");
    assert_eq!(product.original_genbank.as_deref(), Some(original_raw));
    assert_eq!(product.features[0].location, Location::range(5, 14));
    assert_eq!(
        values(&product.features[0], "transl_except"),
        vec!["(pos:9..11,aa:Sec)"]
    );
    assert_eq!(source.sequence.as_str(), "CCCATGTGATAACCC");
    assert_eq!(source.features, original_annotations);
    assert_eq!(source.original_genbank.as_deref(), Some(original_raw));
    assert!(source.metadata.is_empty());
}

#[test]
fn source_bounds_and_expanding_source_evidence_refuse_without_mutating_the_document() {
    let mut too_many = Molecule::dna(
        "Native edit feature ceiling",
        "CCCATGTGATAACCC",
        Topology::Linear,
    )
    .unwrap();
    let original_feature = feature(Location::range(3, 12), &[]);
    too_many.features = vec![original_feature.clone(); 100_001];
    assert!(matches!(too_many.replace(0, 0, "GG"), Err(Error::TooLarge)));
    assert_eq!(too_many.name, "Native edit feature ceiling");
    assert_eq!(too_many.sequence.as_str(), "CCCATGTGATAACCC");
    assert_eq!(too_many.features.len(), 100_001);
    assert!(
        too_many
            .features
            .iter()
            .all(|feature| feature == &original_feature)
    );
    assert!(too_many.metadata.is_empty());
    assert!(too_many.original_genbank.is_none());
    let mut deep = source();
    for _ in 0..33 {
        deep.features[0].location =
            Location::Complement(Box::new(deep.features[0].location.clone()));
    }
    assert!(matches!(deep.replace(0, 0, "GG"), Err(Error::TooLarge)));
    let mut note = source();
    note.features[0].qualifiers = vec![("TrAnSlAtIoN".into(), Some("X".repeat(1_048_576)))];
    let before = note.clone();
    assert!(note.replace(4, 5, "C").is_err());
    assert_eq!(note, before);
}
fn paste_request(destination: Molecule, start: usize, end: usize) -> RequestPaste {
    let source = Molecule::dna("Explicit inserted source", "GG", Topology::Linear).unwrap();
    let selection = Selection {
        parts: vec![Part {
            start: 0,
            end: 2,
            reverse: false,
        }],
        reverse: false,
    };
    let fragment_revision = transfer::copy_revision(&RequestCopy {
        source: source.clone(),
        selection: selection.clone(),
    })
    .unwrap();
    RequestPaste {
        source,
        destination,
        selection,
        fragment_revision,
        start,
        end,
    }
}
#[test]
fn destination_paste_preserves_or_deactivates_assertions_using_the_native_edit_policy() {
    let untouched = paste_request(source(), 0, 0);
    let before = untouched.destination.clone();
    let product = transfer::paste(&untouched).unwrap().product;
    assert_eq!(product.sequence.as_str(), "GGCCCATGTGATAACCC");
    assert_eq!(product.features[0].location, Location::range(5, 14));
    assert_eq!(
        values(&product.features[0], "transl_except"),
        vec!["(pos:9..11,aa:Sec)"]
    );
    assert!(!product.features[0].needs_review);
    assert_eq!(untouched.destination, before);
    let touched = paste_request(source(), 6, 6);
    let product = transfer::paste(&touched).unwrap().product;
    assert_eq!(product.features[0].location, Location::range(3, 14));
    inactive(&product.features[0], &["(pos:7..9,aa:Sec)", "MU"]);
}
#[test]
fn paste_replacement_and_explicit_deletion_keep_later_claims_bound_to_their_actual_feature() {
    let mut destination = source();
    destination.features.insert(
        0,
        feature(Location::range(0, 3), &[("translation", Some("DELETED"))]),
    );
    let request = paste_request(destination, 0, 3);
    let before = request.destination.clone();
    let product = transfer::paste(&request).unwrap().product;
    assert_eq!(product.sequence.as_str(), "GGATGTGATAACCC");
    assert_eq!(product.features.len(), 2);
    assert_eq!(product.features[0].location, Location::range(0, 2));
    inactive(&product.features[0], &["DELETED"]);
    assert_eq!(product.features[1].location, Location::range(2, 11));
    assert_eq!(
        values(&product.features[1], "transl_except"),
        vec!["(pos:6..8,aa:Sec)"]
    );
    assert_eq!(values(&product.features[1], "translation"), vec!["MU"]);
    assert!(!product.features[1].needs_review);
    // Empty replacement really removes the first feature; a subsequent paste
    // must bind the surviving CDS by its current geometry, not its former index.
    let deleted = before.replace(0, 3, "").unwrap();
    assert_eq!(deleted.features.len(), 1);
    assert_eq!(deleted.features[0].location, Location::range(0, 9));
    let inserted = transfer::paste(&paste_request(deleted, 0, 0))
        .unwrap()
        .product;
    assert_eq!(inserted.sequence.as_str(), "GGATGTGATAACCC");
    assert_eq!(inserted.features.len(), 1);
    assert_eq!(inserted.features[0].location, Location::range(2, 11));
    assert_eq!(
        values(&inserted.features[0], "transl_except"),
        vec!["(pos:6..8,aa:Sec)"]
    );
    assert_eq!(values(&inserted.features[0], "translation"), vec!["MU"]);
    assert!(!inserted.features[0].needs_review);
    assert_eq!(request.destination, before);
}
#[test]
fn annotation_only_strand_flips_deactivate_coordinate_claims_instead_of_remapping_them() {
    let mut source = source();
    source.features[0].qualifiers.extend([
        ("RPT_UNIT_RANGE".into(), Some("7..9".into())),
        ("CoDeD_By".into(), Some("NC_REMOTE:1..9".into())),
        ("TRANSLATION".into(), Some("Second translation".into())),
        ("CODON_START".into(), Some("2".into())),
    ]);
    let before = source.clone();
    let product = batch::apply(
        &source,
        &batch::Request {
            source_revision: source.revision(),
            indices: vec![0],
            kind: None,
            flip_strand: true,
            color: ColorEdit::Keep,
        },
    )
    .unwrap();
    assert_eq!(product.sequence, source.sequence);
    assert_eq!(
        product.features[0].location,
        Location::Complement(Box::new(Location::range(3, 12)))
    );
    let f = &product.features[0];
    inactive(
        f,
        &["(pos:7..9,aa:Sec)", "Second translation", "NC_REMOTE:1..9"],
    );
    for key in ["rpt_unit_range", "coded_by"] {
        assert!(values(f, key).is_empty());
    }
    assert_eq!(values(f, "x_opaque"), vec!["first", "second"]);
    assert_eq!(source, before);
}
#[test]
fn arbitrary_annotation_relocation_still_requires_explicit_coordinate_claim_review() {
    let source = source();
    let before = source.clone();
    let mut proposed = source.features[0].clone();
    proposed.location = Location::range(0, 9);
    assert!(set_feature(&source, Some(0), proposed.clone()).is_err());
    proposed.qualifiers.clear();
    assert!(set_feature(&source, Some(0), proposed).is_err());
    assert_eq!(source, before);
}

#[test]
fn non_ascii_whitespace_is_invalid_replacement_data_and_cannot_silently_delete_bases() {
    let source = source();
    let original = source.clone();
    for replacement in ["\u{00a0}", "\u{2003}", "\u{202f}"] {
        assert!(matches!(
            source.replace(4, 5, replacement),
            Err(Error::InvalidBase { .. })
        ));
        assert_eq!(source, original);
    }
    let deleted = source.replace(4, 5, " \n\t").unwrap();
    assert_eq!(deleted.sequence.as_str(), "CCCAGTGATAACCC");
    assert_eq!(source, original);
}

#[test]
fn inserted_bases_between_previously_adjacent_join_parts_invalidate_claims_on_both_strands() {
    for reverse in [false, true] {
        let mut document = source();
        let parts = Location::Join(vec![Location::range(3, 6), Location::range(6, 12)]);
        document.features[0].location = if reverse {
            Location::Complement(Box::new(parts))
        } else {
            parts
        };
        if reverse {
            for (key, value) in &mut document.features[0].qualifiers {
                if key == "transl_except" {
                    *value = Some("(pos:complement(7..9),aa:Sec)".into());
                }
                if key == "anticodon" {
                    *value = Some("(aa:Ser,pos:complement(7..9),seq:tca)".into());
                }
            }
        }
        let original = document.clone();
        let edited = document.replace(6, 6, "GG").unwrap();
        assert_eq!(edited.sequence.as_str(), "CCCATGGGTGATAACCC");
        let parts = Location::Join(vec![Location::range(3, 6), Location::range(8, 14)]);
        assert_eq!(
            edited.features[0].location,
            if reverse {
                Location::Complement(Box::new(parts))
            } else {
                parts
            }
        );
        inactive(&edited.features[0], &["aa:Sec", "MU"]);
        assert_eq!(document, original);
    }
}

#[test]
fn insertion_at_the_circular_origin_interrupts_a_crossing_gene_on_both_strands() {
    for reverse in [false, true] {
        let mut document =
            Molecule::dna("Origin-spanning gene", "GTGATAACCCAT", Topology::Circular).unwrap();
        let parts = Location::Join(vec![Location::range(9, 12), Location::range(0, 6)]);
        document.features.push(feature(
            if reverse {
                Location::Complement(Box::new(parts))
            } else {
                parts
            },
            &[
                (
                    "transl_except",
                    Some(if reverse {
                        "(pos:complement(1..3),aa:Sec)"
                    } else {
                        "(pos:1..3,aa:Sec)"
                    }),
                ),
                ("translation", Some("SOURCE_CDS")),
                ("codon_start", Some("1")),
            ],
        ));
        let original = document.clone();
        let edited = document.replace(0, 0, "GG").unwrap();
        assert_eq!(edited.sequence.as_str(), "GGGTGATAACCCAT");
        let parts = Location::Join(vec![Location::range(11, 14), Location::range(2, 8)]);
        assert_eq!(
            edited.features[0].location,
            if reverse {
                Location::Complement(Box::new(parts))
            } else {
                parts
            }
        );
        inactive(&edited.features[0], &["aa:Sec", "SOURCE_CDS"]);
        assert_eq!(document, original);
    }
}

#[test]
fn full_native_feature_capacity_remains_usable_for_untouched_single_span_annotations() {
    let mut document = Molecule::dna("Annotation capacity", "ACGT", Topology::Linear).unwrap();
    let annotation = Feature {
        kind: "misc_feature".into(),
        location: Location::range(0, 4),
        needs_review: false,
        qualifiers: vec![
            ("label".into(), Some("capacity".into())),
            ("note".into(), Some("retained".into())),
        ],
    };
    document.features = vec![annotation.clone(); 100_000];
    let edited = document.replace(0, 0, "GG").unwrap();
    assert_eq!(edited.sequence.as_str(), "GGACGT");
    assert_eq!(edited.features.len(), 100_000);
    assert!(
        edited
            .features
            .iter()
            .all(|feature| feature.location == Location::range(2, 6)
                && feature.kind == annotation.kind
                && feature.qualifiers == annotation.qualifiers
                && !feature.needs_review)
    );
    assert_eq!(document.sequence.as_str(), "ACGT");
    assert!(
        document
            .features
            .iter()
            .all(|feature| feature == &annotation)
    );
}
