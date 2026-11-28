//! Independent location fixtures. Normative descriptor examples are from INSDC's
//! Feature Table Definition v11.4 (April 2026), sections 3.4.2 and 3.4.3:
//! https://www.insdc.org/submitting-standards/feature-table/
//! Other cases are synthetic adversarial inputs or explicitly retained grammar
//! extensions; they are not claims that every extension is a current INSDC form.
use collagene_core::{formats::sequence, protocol};
use collagene_model::{Feature, Location, Molecule, Topology};

fn range(start: usize, end: usize, before: bool, after: bool) -> Location {
    Location::Range {
        start,
        end,
        before,
        after,
    }
}

fn rejected_without_panic(text: &str, length: usize) {
    let outcome = std::panic::catch_unwind(|| sequence::parse_location(text, length));
    assert!(outcome.is_ok(), "location parser panicked for {text:?}");
    assert!(
        outcome.unwrap().is_err(),
        "accepted malformed/out-of-bounds location {text:?}"
    );
}

fn genbank(features: &str) -> String {
    format!(
        "LOCUS       integrity                 32 bp    DNA     linear   SYN 28-NOV-2026\nDEFINITION  Synthetic location integrity fixture.\nFEATURES             Location/Qualifiers\n{features}ORIGIN\n        1 acgtacgtacgtacgtacgtacgtacgtacgt\n//\n"
    )
}

fn rejected_record_without_panic(features: &str) {
    let record = genbank(features);
    let outcome = std::panic::catch_unwind(|| sequence::import_genbank(&record));
    assert!(outcome.is_ok(), "GenBank import panicked for {features:?}");
    assert!(
        outcome.unwrap().is_err(),
        "GenBank import accepted invalid feature text {features:?}"
    );
}

#[test]
fn location_parser_consumes_every_numeric_compound_and_remote_token() {
    for text in [
        "1junk",
        "1..3junk",
        "1..3..7",
        "1..3,7..9",
        "1..3)",
        "join(1..3,7..9)junk",
        "order(1..3,7..9)junk",
        "complement(1..3)junk",
        "complement(join(1..3,7..9))junk",
        "J00194.1:100..202junk",
        "J00194.1:100..202)",
        "join(1..3,J00194.1:100..202)junk",
        "A(:1..3)junk",
        "gap(100)junk",
        "bond(12,15)junk",
        "123^124junk",
        "102.110junk",
    ] {
        rejected_without_panic(text, 6000);
    }
}

#[test]
fn broken_location_shapes_never_become_prefix_ranges_or_unresolved_successes() {
    for text in [
        "",
        "junk",
        "unknown(1..3)",
        "join()",
        "order()",
        "complement()",
        "complement(1..3,7..9)",
        "join(,1..3)",
        "join(1..3,)",
        "join(1..3,,7..9)",
        "join(1..3,7..9",
        "1..3)",
        "(1..3)",
        "J00194.1:",
        ":1..3",
        "J00194.1:1..3:4..6",
        "3^7",
        "1..\n3",
        "1..\t3",
        "1..3\0junk",
        "1..3é",
    ] {
        rejected_without_panic(text, 6000);
    }
}

#[test]
fn signed_minima_and_unsigned_overflows_return_errors_without_panics() {
    for text in [
        "-9223372036854775808",
        "-9223372036854775808..3",
        "complement(-9223372036854775808..3)",
        "join(1..3,-9223372036854775808..4)",
        "-1",
        "+1",
        "0",
        "0..3",
        "3..0",
        "9..3",
        "18446744073709551615",
        "18446744073709551616",
        "1..18446744073709551616",
        "J00194.1:18446744073709551616",
        "gap(18446744073709551616)",
        "gap(unk18446744073709551616)",
        "bond(1,18446744073709551616)",
    ] {
        rejected_without_panic(text, 32);
    }
    for text in ["33", "1..33", "complement(33)", "join(1..3,32..33)"] {
        rejected_without_panic(text, 32);
    }
}

#[test]
fn location_resource_bounds_cover_depth_width_and_text_size() {
    let deep = format!("{}1..3{}", "complement(".repeat(33), ")".repeat(33));
    let wide = format!("join({})", vec!["1"; 4096].join(","));
    let long_number = "9".repeat(65_537);
    let long_remote = format!("{}:1..3", "A".repeat(65_537));
    for text in [&deep, &wide, &long_number, &long_remote] {
        rejected_without_panic(text, 6000);
    }
}

#[test]
fn insdc_local_examples_preserve_exact_tree_coordinates_and_roundtrip() {
    let cases = [
        ("467", range(466, 467, false, false)),
        ("340..565", range(339, 565, false, false)),
        ("<345..500", range(344, 500, true, false)),
        ("<1..888", range(0, 888, true, false)),
        ("1..>888", range(0, 888, false, true)),
        (
            "join(12..78,134..202)",
            Location::Join(vec![Location::range(11, 78), Location::range(133, 202)]),
        ),
        (
            "order(12..78,134..202)",
            Location::Order(vec![Location::range(11, 78), Location::range(133, 202)]),
        ),
        (
            "complement(34..126)",
            Location::Complement(Box::new(Location::range(33, 126))),
        ),
        (
            "complement(join(2691..4571,4918..5163))",
            Location::Complement(Box::new(Location::Join(vec![
                Location::range(2690, 4571),
                Location::range(4917, 5163),
            ]))),
        ),
        (
            "join(complement(4918..5163),complement(2691..4571))",
            Location::Join(vec![
                Location::Complement(Box::new(Location::range(4917, 5163))),
                Location::Complement(Box::new(Location::range(2690, 4571))),
            ]),
        ),
        (
            "join(544..589,688..>1032)",
            Location::Join(vec![
                Location::range(543, 589),
                range(687, 1032, false, true),
            ]),
        ),
    ];
    for (text, expected) in cases {
        let actual = sequence::parse_location(text, 6000).unwrap();
        assert_eq!(actual, expected, "{text}");
        let rendered = sequence::location_text(&actual).unwrap();
        assert_eq!(
            sequence::parse_location(&rendered, 6000).unwrap(),
            expected,
            "{text} -> {rendered}"
        );
    }
}

#[test]
fn unsupported_valid_descriptors_retain_whole_text_without_local_sequence_authority() {
    // Initial examples exercise INSDC remote/between/uncertain descriptors.
    // Gap/bond/one-of forms below exercise the explicitly retained parser profile.
    for text in [
        "J00194.1:100..202",
        "join(1..100,J00194.1:100..202)",
        "123^124",
        "128^1",
        "102.110",
        "<23",
        ">23",
        "gap()",
        "gap(100)",
        "gap(unk100)",
        "bond(12,15)",
        "one-of(12,15)",
        "complement(J00194.1:100..202)",
        "join(1..3,gap(100),7..9)",
    ] {
        let location = sequence::parse_location(text, 128).unwrap();
        assert_eq!(location, Location::Unresolved(text.into()), "{text}");
        assert_eq!(sequence::location_text(&location).unwrap(), text);
        assert!(
            location.spans().is_empty(),
            "unresolved coordinates must not become local spans"
        );
        let mut molecule = Molecule::dna(
            "Unresolved location",
            &"ACGT".repeat(32),
            Topology::Circular,
        )
        .unwrap();
        assert!(
            location.extract(&molecule.sequence).is_err(),
            "unresolved sequence extraction must be refused"
        );
        molecule.features.push(Feature {
            kind: "misc_feature".into(),
            location,
            qualifiers: vec![],
            needs_review: true,
        });
        let bytes = protocol::save_document(&molecule).unwrap();
        assert_eq!(protocol::load_document(&bytes).unwrap(), molecule);
    }
}

#[test]
fn raw_genbank_feature_locations_are_checked_before_any_prefix_can_be_accepted() {
    for location in [
        "1..3junk",
        "complement(1..3)junk",
        "join(1..3,7..9)junk",
        "J00194.1:1..3junk",
        "A(:1..3)junk",
        "-9223372036854775808..3",
        "complement(-9223372036854775808..3)",
        "1..18446744073709551616",
    ] {
        rejected_record_without_panic(&format!(
            "     misc_feature    {location}\n                     /label=\"invalid\"\n"
        ));
    }
    for features in [
        "     misc_feature    complement(join(1..3,\n                     7..9))junk\n                     /label=\"invalid continuation\"\n",
        "     misc_feature    1..3\n                     junk\n                     /label=\"trailing continuation\"\n",
        "     misc_feature    join(1..3,\n                     7..9\n                     /label=\"unclosed\"\n",
        "     misc_feature    \n                     /label=\"missing location\"\n",
    ] {
        rejected_record_without_panic(features);
    }
}

#[test]
fn continued_genbank_locations_and_qualifier_lookalikes_keep_their_own_boundaries() {
    let features = "     misc_feature    complement(join(1..4,\n                     13..16))\n                     /label=\"continued\"\n                     /note=\"A(:1..3)junk is qualifier text, not a location\"\n     misc_feature    order(17..20,\n                     25..28)\n                     /label=\"ordered\"\n";
    let source = genbank(features);
    let imported = sequence::import_genbank(&source).unwrap();
    assert_eq!(imported.molecules.len(), 1);
    let molecule = &imported.molecules[0];
    assert_eq!(molecule.original_genbank.as_deref(), Some(source.as_str()));
    assert_eq!(molecule.features.len(), 2);
    let first = Location::Complement(Box::new(Location::Join(vec![
        Location::range(0, 4),
        Location::range(12, 16),
    ])));
    let second = Location::Order(vec![Location::range(16, 20), Location::range(24, 28)]);
    assert_eq!(molecule.features[0].location, first);
    assert_eq!(molecule.features[1].location, second);
    assert!(
        molecule.features[0]
            .qualifiers
            .iter()
            .any(|(key, value)| key == "note"
                && value.as_deref() == Some("A(:1..3)junk is qualifier text, not a location"))
    );
    let output = sequence::export_genbank(molecule).unwrap();
    let reopened = sequence::import_genbank(&output).unwrap();
    assert_eq!(reopened.molecules[0].features[0].location, first);
    assert_eq!(reopened.molecules[0].features[1].location, second);
}

#[test]
fn raw_record_location_budgets_and_later_invalid_records_do_not_return_partial_success() {
    let wide = format!("join({})", vec!["1"; 4096].join(","));
    let deep = format!("{}1..3{}", "complement(".repeat(33), ")".repeat(33));
    let oversized = "9".repeat(65_537);
    for text in [&wide, &deep, &oversized] {
        rejected_record_without_panic(&format!("     misc_feature    {text}\n"));
    }
    let valid =
        genbank("     misc_feature    1..3\n                     /label=\"valid first record\"\n");
    let invalid = genbank("     misc_feature    complement(1..3)junk\n");
    let records = format!("{valid}{invalid}");
    let outcome = std::panic::catch_unwind(|| sequence::import_genbank(&records));
    assert!(outcome.is_ok(), "later malformed records must not panic");
    assert!(
        outcome.unwrap().is_err(),
        "a later invalid location must reject the complete import"
    );
}

#[test]
fn raw_genbank_contig_descriptors_cannot_reach_an_unsafe_prefix_parser() {
    let valid = genbank(
        "     misc_feature    1..3\n                     /label=\"ordinary local feature\"\n",
    );
    for descriptor in [
        "-9223372036854775808",
        "-9223372036854775808..3",
        "join(-9223372036854775808..3,4..6)",
        "J00194.1:1..3junk",
        "join(J00194.1:1..3)junk",
        "A(:1..3)junk",
        "J00194.1:1..18446744073709551616",
    ] {
        let record = valid.replace("ORIGIN\n", &format!("CONTIG      {descriptor}\nORIGIN\n"));
        let outcome = std::panic::catch_unwind(|| sequence::import_genbank(&record));
        assert!(outcome.is_ok(), "CONTIG parser panicked for {descriptor:?}");
        assert!(
            outcome.unwrap().is_err(),
            "invalid CONTIG descriptor was accepted: {descriptor:?}"
        );
    }
    let continued = valid.replace(
        "ORIGIN\n",
        "CONTIG      join(J00194.1:1..3,\n            J00194.1:7..9)junk\nORIGIN\n",
    );
    let outcome = std::panic::catch_unwind(|| sequence::import_genbank(&continued));
    assert!(
        outcome.is_ok(),
        "continued CONTIG descriptor must not panic"
    );
    assert!(
        outcome.unwrap().is_err(),
        "trailing CONTIG continuation must be rejected"
    );
}

#[test]
fn malformed_qualifiers_never_recover_as_feature_locations_or_panic() {
    for qualifier in [
        "/note=\"unclosed -9223372036854775808",
        "/note=\"ok\" junk -9223372036854775808",
        "/pseudo junk -9223372036854775808",
        "/note bare -9223372036854775808",
    ] {
        rejected_record_without_panic(&format!(
            "     misc_feature    1..3\n                     {qualifier}\n"
        ));
    }
    rejected_record_without_panic(concat!(
        "     misc_feature    1..3\n",
        "                     /note=\"unterminated\n",
        "     misc_feature    -9223372036854775808\n",
    ));
}

#[test]
fn location_like_and_signed_minimum_quoted_annotations_remain_literal_data() {
    let first = "literal -9223372036854775808 and A(:1..3)junk";
    let second = "literal join(1..3)junk /note=not-a-new-qualifier";
    let source = genbank(&format!(
        "     misc_feature    1..3\n                     /note=\"{first}\"\n                     /note=\"{second}\"\n"
    ));
    let result = std::panic::catch_unwind(|| sequence::import_genbank(&source));
    assert!(
        result.is_ok(),
        "quoted annotation data must not reach any location parser"
    );
    let imported = result.unwrap().unwrap();
    let molecule = &imported.molecules[0];
    assert_eq!(molecule.features.len(), 1);
    assert_eq!(molecule.features[0].location, Location::range(0, 3));
    assert_eq!(
        molecule.features[0].qualifiers,
        vec![
            ("note".into(), Some(first.into())),
            ("note".into(), Some(second.into())),
        ]
    );
    let exported = sequence::export_genbank(molecule).unwrap();
    let reopened = sequence::import_genbank(&exported).unwrap();
    assert_eq!(reopened.molecules[0].features, molecule.features);
}

#[test]
fn tab_separated_features_headers_receive_the_same_location_safety_checks() {
    for header in [
        "FEATURES\tLocation/Qualifiers",
        "FEATURES \t Location/Qualifiers",
    ] {
        let source =
            genbank("     misc_feature    1..3\n                     /label=\"tabbed header\"\n")
                .replace("FEATURES             Location/Qualifiers", header);
        let valid = std::panic::catch_unwind(|| sequence::import_genbank(&source));
        assert!(
            valid.is_ok(),
            "valid tab-separated FEATURES header must not panic"
        );
        let imported = valid.unwrap().unwrap();
        assert_eq!(imported.molecules[0].features.len(), 1);
        assert_eq!(
            imported.molecules[0].features[0].location,
            Location::range(0, 3)
        );
        let invalid = source.replace(
            "misc_feature    1..3",
            "misc_feature    -9223372036854775808",
        );
        let rejected = std::panic::catch_unwind(|| sequence::import_genbank(&invalid));
        assert!(
            rejected.is_ok(),
            "tab-separated FEATURES must not bypass location preflight"
        );
        assert!(
            rejected.unwrap().is_err(),
            "invalid location after tab-separated FEATURES was accepted"
        );
    }
}

#[test]
fn ordered_repeated_qualifiers_flags_quotes_and_continuations_roundtrip() {
    let source = genbank(concat!(
        "     misc_feature    1..9\n",
        "                     /note=\"first\"\n",
        "                     /pseudo\n",
        "                     /note=\"two \"\"quoted\"\" words\"\n",
        "                     /gene=\"example\"\n",
        "                     /note=\"first prose line\n",
        "                     second prose line\"\n",
        "                     /translation=\"TY\n",
        "                     V\"\n",
        "                     /codon_start=1\n",
        "                     /ribosomal_slippage\n",
        "                     /note=\"last\"\n",
    ));
    let imported = sequence::import_genbank(&source).unwrap();
    let molecule = &imported.molecules[0];
    let expected = [
        ("note", Some("first")),
        ("pseudo", None),
        ("note", Some("two \"quoted\" words")),
        ("gene", Some("example")),
        ("note", Some("first prose line\nsecond prose line")),
        ("translation", Some("TYV")),
        ("codon_start", Some("1")),
        ("ribosomal_slippage", None),
        ("note", Some("last")),
    ]
    .into_iter()
    .map(|(key, value)| (key.to_owned(), value.map(str::to_owned)))
    .collect::<Vec<_>>();
    assert_eq!(molecule.features[0].qualifiers, expected);
    assert_eq!(molecule.original_genbank.as_deref(), Some(source.as_str()));
    let native = protocol::save_document(molecule).unwrap();
    assert_eq!(protocol::load_document(&native).unwrap(), *molecule);
    let exported = sequence::export_genbank(molecule).unwrap();
    let reopened = sequence::import_genbank(&exported).unwrap();
    assert_eq!(reopened.molecules[0].features[0].qualifiers, expected);
    assert_eq!(
        reopened.molecules[0].features[0].location,
        molecule.features[0].location
    );
}

#[test]
fn valid_uninterpreted_locations_survive_complete_genbank_and_native_roundtrips() {
    for text in [
        "J00194.1:100..202",
        "join(1..3,J00194.1:100..202)",
        "complement(J00194.1:100..202)",
        "3^4",
        "3.9",
        "<23",
        ">23",
        "gap()",
        "gap(100)",
        "gap(unk100)",
        "bond(12,15)",
        "one-of(12,15)",
        "join(1..3,gap(100),7..9)",
    ] {
        let source = genbank(&format!(
            "     misc_feature    {text}\n                     /label=\"opaque source evidence\"\n"
        ));
        let result = std::panic::catch_unwind(|| sequence::import_genbank(&source));
        assert!(
            result.is_ok(),
            "valid uninterpreted location import panicked: {text}"
        );
        let imported = result.unwrap().unwrap();
        let molecule = &imported.molecules[0];
        let expected = Location::Unresolved(text.into());
        assert_eq!(molecule.features[0].location, expected, "{text}");
        assert!(
            molecule.features[0]
                .location
                .extract(&molecule.sequence)
                .is_err()
        );
        let native = protocol::save_document(molecule).unwrap();
        assert_eq!(protocol::load_document(&native).unwrap(), *molecule);
        let exported = sequence::export_genbank(molecule).unwrap();
        let reopened = sequence::import_genbank(&exported).unwrap();
        assert_eq!(
            reopened.molecules[0].features[0].location, expected,
            "{text}"
        );
        assert_eq!(
            reopened.molecules[0].features[0].qualifiers,
            molecule.features[0].qualifiers
        );
    }
}

#[test]
fn blank_qualifier_line_flood_is_bounded_before_structural_allocation() {
    let mut features =
        "     misc_feature    1..3\n                     /note=\"bounded\"\n".to_owned();
    features.push_str(&"\n".repeat(1_000_001));
    let source = genbank(&features);
    assert!(source.len() < sequence::MAX_IMPORT_BYTES);
    let result = std::panic::catch_unwind(|| sequence::import_genbank(&source));
    assert!(result.is_ok());
    assert!(matches!(
        result.unwrap(),
        Err(collagene_model::Error::TooLarge)
    ));
}

#[test]
fn native_feature_fields_cannot_inject_flat_file_structure() {
    for kind in ["gene\0", "gene\t", "gene\nORIGIN", "key longer than field"] {
        rejected_record_without_panic(&format!("     {kind}    1..3\n"));
        let mut molecule = Molecule::dna("Untrusted fields", "ACGT", Topology::Linear).unwrap();
        molecule.features.push(Feature {
            kind: kind.into(),
            location: Location::range(0, 3),
            qualifiers: vec![],
            needs_review: false,
        });
        assert!(sequence::export_genbank(&molecule).is_err());
    }
    for (key, value) in [("note\nORIGIN", "literal"), ("note", "control\0text")] {
        let mut molecule = Molecule::dna("Untrusted fields", "ACGT", Topology::Linear).unwrap();
        molecule.features.push(Feature {
            kind: "misc_feature".into(),
            location: Location::range(0, 3),
            qualifiers: vec![(key.into(), Some(value.into()))],
            needs_review: false,
        });
        assert!(sequence::export_genbank(&molecule).is_err());
    }
    for (key, value) in [
        ("k".repeat(100), "long-key"),
        ("anticodon".into(), "(pos:1..3,aa:🧬)"),
        ("transl_except".into(), "(pos:1..3,\naa:Sec)"),
    ] {
        let mut molecule = Molecule::dna("Untrusted widths", "ACGT", Topology::Linear).unwrap();
        molecule.features.push(Feature {
            kind: "misc_feature".into(),
            location: Location::range(0, 3),
            qualifiers: vec![(key, Some(value.into()))],
            needs_review: false,
        });
        let result = std::panic::catch_unwind(|| sequence::export_genbank(&molecule));
        assert!(result.is_ok());
        assert!(result.unwrap().is_err());
    }
}

#[test]
fn quoted_value_boundaries_preserve_literal_characters_without_writer_panics() {
    for value in [
        format!("{}\"", "A".repeat(50)),
        format!("{}\"Ω", "λ".repeat(70)),
    ] {
        let mut molecule = Molecule::dna("Boundary quotes", "ACGT", Topology::Linear).unwrap();
        molecule.features.push(Feature {
            kind: "misc_feature".into(),
            location: Location::range(0, 3),
            qualifiers: vec![("note".into(), Some(value.clone()))],
            needs_review: false,
        });
        let result = std::panic::catch_unwind(|| sequence::export_genbank(&molecule));
        assert!(result.is_ok());
        let text = result.unwrap().unwrap();
        let imported = sequence::import_genbank(&text).unwrap();
        let actual = imported.molecules[0].features[0].qualifiers[0]
            .1
            .as_ref()
            .unwrap();
        assert_eq!(actual.replace('\n', ""), value);
    }
}
