//! Hand-specified genetic-code/coordinate fixtures, not outputs of another cloning application.
use collagene_core::sequence::translation::*;
use collagene_model::{Alphabet, Error, Feature, Location, Molecule, Sequence, Topology};
fn dna(s: &str) -> Sequence {
    Sequence::dna(s).unwrap()
}
fn molecule(s: &str, topology: Topology) -> Molecule {
    Molecule::dna("translation fixture", s, topology).unwrap()
}
fn cds(location: Location) -> Feature {
    Feature {
        kind: "CDS".into(),
        location,
        qualifiers: vec![("label".into(), Some("coding feature".into()))],
        needs_review: false,
    }
}
fn orf_options() -> OrfOptions {
    OrfOptions {
        min_amino_acids: 1,
        both_strands: false,
        ..OrfOptions::default()
    }
}

#[test]
fn standard_amino_acids_stops_and_unique_iupac_codons_match_independent_expectations() {
    for (codon, amino) in [
        ("GCT", "A"),
        ("TGT", "C"),
        ("GAT", "D"),
        ("GAA", "E"),
        ("TTT", "F"),
        ("GGT", "G"),
        ("CAT", "H"),
        ("ATT", "I"),
        ("AAA", "K"),
        ("CTT", "L"),
        ("ATG", "M"),
        ("AAT", "N"),
        ("CCT", "P"),
        ("CAA", "Q"),
        ("CGT", "R"),
        ("TCT", "S"),
        ("ACT", "T"),
        ("GTT", "V"),
        ("TGG", "W"),
        ("TAT", "Y"),
        ("TAA", "*"),
        ("TAG", "*"),
        ("TGA", "*"),
    ] {
        assert_eq!(
            translate(&dna(codon), &TranslationOptions::default())
                .unwrap()
                .protein,
            amino,
            "{codon}"
        );
    }
    let report = translate(&dna("GCNATHATNTARTANGAR"), &TranslationOptions::default()).unwrap();
    assert_eq!(report.protein, "AIX*XE");
    assert_eq!(report.ambiguous_codons, 6);
    assert_eq!(report.unresolved_codons, 2);
    assert_eq!(report.possible_stop_codons, 1);
    assert_eq!(
        translate(&dna("TRA"), &TranslationOptions::default())
            .unwrap()
            .protein,
        "*"
    );
}

#[test]
fn frames_reverse_rna_and_incomplete_codons_are_explicit() {
    assert_eq!(
        translate(
            &dna("AATGAAATAA"),
            &TranslationOptions {
                frame: 1,
                ..TranslationOptions::default()
            }
        )
        .unwrap()
        .protein,
        "MK*"
    );
    assert_eq!(
        translate(
            &dna("CCATGAAATAA"),
            &TranslationOptions {
                frame: 2,
                ..TranslationOptions::default()
            }
        )
        .unwrap()
        .protein,
        "MK*"
    );
    assert_eq!(
        translate(
            &dna("TTATTTCAT"),
            &TranslationOptions {
                reverse: true,
                ..TranslationOptions::default()
            }
        )
        .unwrap()
        .protein,
        "MK*"
    );
    let rna = Sequence::new(Alphabet::Rna, "AUGGCNUAAA").unwrap();
    let report = translate(&rna, &TranslationOptions::default()).unwrap();
    assert_eq!(report.protein, "MA*");
    assert_eq!(report.trailing_bases, 1);
    assert!(
        translate(
            &dna("ATG"),
            &TranslationOptions {
                frame: 3,
                ..TranslationOptions::default()
            }
        )
        .is_err()
    );
    assert!(
        translate(
            &Sequence::new(Alphabet::Protein, "MAK").unwrap(),
            &TranslationOptions::default()
        )
        .is_err()
    );
    let empty = translate(&dna("A"), &TranslationOptions::default()).unwrap();
    assert!(empty.protein.is_empty());
    assert_eq!(empty.trailing_bases, 1);
}

#[test]
fn alternative_start_recoding_and_code_two_differences_are_recorded() {
    let ordinary = translate(
        &dna("GTGGTGTAA"),
        &TranslationOptions {
            genetic_code: 11,
            ..TranslationOptions::default()
        },
    )
    .unwrap();
    assert_eq!(ordinary.protein, "VV*");
    let started = translate(
        &dna("GTGGTGTAA"),
        &TranslationOptions {
            genetic_code: 11,
            initiation: InitiationPolicy::RequireRecognizedStart,
            ..TranslationOptions::default()
        },
    )
    .unwrap();
    assert_eq!(started.protein, "MV*");
    assert!(started.initiator_methionine_applied);
    assert!(
        translate(
            &dna("GTGGTGTAA"),
            &TranslationOptions {
                initiation: InitiationPolicy::RequireRecognizedStart,
                ..TranslationOptions::default()
            }
        )
        .is_err()
    );
    let mitochondrial = TranslationOptions {
        genetic_code: 2,
        ..TranslationOptions::default()
    };
    assert_eq!(
        translate(&dna("ATATGAAGAGGT"), &mitochondrial)
            .unwrap()
            .protein,
        "MW*G"
    );
    assert_eq!(
        translate(&dna("AGRTGR"), &mitochondrial).unwrap().protein,
        "*W"
    );
    assert!(
        translate(
            &dna("ATG"),
            &TranslationOptions {
                genetic_code: 99,
                ..TranslationOptions::default()
            }
        )
        .is_err()
    );
    assert_eq!(
        genetic_codes().iter().map(|c| c.id).collect::<Vec<_>>(),
        [
            1, 2, 3, 4, 5, 6, 9, 10, 11, 12, 13, 14, 15, 16, 21, 22, 23, 24, 25, 26, 27, 28, 29,
            30, 31, 32, 33
        ]
    );
}

#[test]
fn cds_join_and_reverse_origin_features_translate_in_biological_order() {
    let m = molecule("AAATAACCCATG", Topology::Circular);
    let feature = cds(Location::Join(vec![
        Location::range(9, 12),
        Location::range(0, 6),
    ]));
    let report = translate_cds(&m, &feature).unwrap();
    assert_eq!(report.protein, "MK");
    assert!(report.complete_under_model);
    assert!(report.has_terminal_stop);
    let reverse = molecule("CATGGGTTATTT", Topology::Circular);
    let feature = cds(Location::Complement(Box::new(Location::Join(vec![
        Location::range(6, 12),
        Location::range(0, 3),
    ]))));
    assert_eq!(translate_cds(&reverse, &feature).unwrap().protein, "MK");
}

#[test]
fn cds_phase_internal_stops_and_stored_translation_do_not_become_silent_passes() {
    let m = molecule("CCGTGAAATAA", Topology::Linear);
    let mut f = cds(Location::range(0, 11));
    f.qualifiers.extend([
        ("codon_start".into(), Some("3".into())),
        ("transl_table".into(), Some("11".into())),
    ]);
    let report = translate_cds(&m, &f).unwrap();
    assert_eq!(report.protein, "VK");
    assert!(!report.translation.initiator_methionine_applied);
    assert!(!report.complete_under_model);
    let m = molecule("ATGTAAGCTTAA", Topology::Linear);
    let mut f = cds(Location::range(0, 12));
    f.qualifiers.push(("translation".into(), Some("MA".into())));
    let report = translate_cds(&m, &f).unwrap();
    assert_eq!(report.protein, "M*A");
    assert_eq!(report.internal_stops, 1);
    assert_eq!(report.stored_translation_matches, Some(false));
    assert!(!report.complete_under_model);
    let m = molecule("GTGAAATAA", Topology::Linear);
    let mut f = cds(Location::range(0, 9));
    f.qualifiers
        .push(("transl_table".into(), Some("11".into())));
    assert_eq!(translate_cds(&m, &f).unwrap().protein, "MK");
}

#[test]
fn unsupported_cds_biology_and_malformed_qualifiers_are_rejected() {
    let m = molecule("ATGAAATAA", Topology::Linear);
    for key in [
        "transl_except",
        "exception",
        "ribosomal_slippage",
        "frameshift",
        "pseudo",
    ] {
        let mut f = cds(Location::range(0, 9));
        f.qualifiers.push((key.into(), Some("data".into())));
        assert!(
            matches!(translate_cds(&m, &f), Err(Error::Unsupported(_))),
            "{key}"
        );
    }
    for (key, value) in [
        ("codon_start", "0"),
        ("codon_start", "4"),
        ("codon_start", "abc"),
        ("transl_table", ""),
        ("transl_table", "-1"),
    ] {
        let mut f = cds(Location::range(0, 9));
        f.qualifiers.push((key.into(), Some(value.into())));
        assert!(translate_cds(&m, &f).is_err());
    }
    let mut f = cds(Location::range(0, 9));
    f.qualifiers.extend([
        ("codon_start".into(), Some("1".into())),
        ("codon_start".into(), Some("1".into())),
    ]);
    assert!(translate_cds(&m, &f).is_err());
    let fuzzy = cds(Location::Range {
        start: 0,
        end: 9,
        before: true,
        after: false,
    });
    assert!(translate_cds(&m, &fuzzy).is_err());
    let ordered = cds(Location::Order(vec![
        Location::range(0, 3),
        Location::range(3, 9),
    ]));
    assert!(translate_cds(&m, &ordered).is_err());
    let mut reviewed = cds(Location::range(0, 9));
    reviewed.needs_review = true;
    assert!(translate_cds(&m, &reviewed).is_err());
}

#[test]
fn linear_and_reverse_orfs_have_exact_reference_coordinates() {
    let report = find_orfs(
        &molecule("CCCATGAAATAACCC", Topology::Linear),
        &orf_options(),
    )
    .unwrap();
    assert_eq!(report.orfs.len(), 1);
    let orf = &report.orfs[0];
    assert_eq!(orf.location, Location::range(3, 12));
    assert_eq!(orf.protein, "MK");
    assert_eq!(orf.nucleotide_length, 9);
    assert!(!orf.partial_5prime && !orf.partial_3prime);
    let m = molecule("GGGTTATTTCATGGG", Topology::Linear);
    let options = OrfOptions {
        both_strands: true,
        ..orf_options()
    };
    let report = find_orfs(&m, &options).unwrap();
    let orf = report.orfs.iter().find(|f| f.reverse).unwrap();
    assert_eq!(
        orf.location,
        Location::Complement(Box::new(Location::range(3, 12)))
    );
    assert_eq!(orf.protein, "MK");
    assert_eq!(orf.location.extract(&m.sequence).unwrap(), "ATGAAATAA");
}

#[test]
fn circular_origin_orfs_respect_single_turn_even_when_length_is_not_divisible_by_three() {
    for (bases, start, n) in [("AAATAACCCATG", 9, 12), ("AAATAACCATG", 8, 11)] {
        let m = molecule(bases, Topology::Circular);
        let report = find_orfs(&m, &orf_options()).unwrap();
        let orf = report
            .orfs
            .iter()
            .find(|o| o.oriented_start == start)
            .unwrap();
        assert!(orf.crosses_origin);
        assert_eq!(orf.protein, "MK");
        assert_eq!(orf.nucleotide_length, 9);
        assert_eq!(
            orf.location,
            Location::Join(vec![Location::range(start, n), Location::range(0, 6)])
        );
        assert_eq!(orf.location.extract(&m.sequence).unwrap(), "ATGAAATAA");
    }
    let m = molecule("CATGGGTTATTT", Topology::Circular);
    let report = find_orfs(
        &m,
        &OrfOptions {
            both_strands: true,
            ..orf_options()
        },
    )
    .unwrap();
    let orf = report
        .orfs
        .iter()
        .find(|o| o.reverse && o.crosses_origin)
        .unwrap();
    assert_eq!(orf.location.extract(&m.sequence).unwrap(), "ATGAAATAA");
    assert_eq!(orf.protein, "MK");
    let m = molecule("ATGAAAAA", Topology::Circular);
    assert!(find_orfs(&m, &orf_options()).unwrap().orfs.is_empty());
    let report = find_orfs(
        &m,
        &OrfOptions {
            incomplete: IncompletePolicy::IncludeAllPartials,
            ..orf_options()
        },
    )
    .unwrap();
    let orf = &report.orfs[0];
    assert_eq!(orf.nucleotide_length, 6);
    assert_eq!(orf.protein, "MK");
    assert_eq!(orf.termination, Termination::OneTurnLimit);
    assert!(orf.partial_3prime);
}

#[test]
fn partial_and_ambiguity_policies_are_explicit_and_do_not_invent_complete_orfs() {
    let m = molecule("CCCAAATAA", Topology::Linear);
    assert!(find_orfs(&m, &orf_options()).unwrap().orfs.is_empty());
    let report = find_orfs(
        &m,
        &OrfOptions {
            incomplete: IncompletePolicy::IncludeLinearPartials,
            ..orf_options()
        },
    )
    .unwrap();
    let prefix = report.orfs.iter().find(|o| o.oriented_start == 0).unwrap();
    assert_eq!(prefix.protein, "PK");
    assert!(prefix.partial_5prime);
    assert!(prefix.has_stop_codon);
    let m = molecule("ATGAAATANTAA", Topology::Linear);
    assert!(find_orfs(&m, &orf_options()).unwrap().orfs.is_empty());
    let report = find_orfs(
        &m,
        &OrfOptions {
            incomplete: IncompletePolicy::IncludeAllPartials,
            ..orf_options()
        },
    )
    .unwrap();
    let orf = &report.orfs[0];
    assert_eq!(orf.protein, "MK");
    assert_eq!(orf.nucleotide_length, 6);
    assert_eq!(orf.termination, Termination::AmbiguousStop);
    assert!(orf.partial_3prime);
    let m = molecule("ATHAAATAA", Topology::Linear);
    let options = OrfOptions {
        genetic_code: 11,
        start_policy: StartPolicy::CodeStarts,
        ..orf_options()
    };
    assert_eq!(find_orfs(&m, &options).unwrap().orfs[0].protein, "MK");
    assert!(find_orfs(&m, &orf_options()).unwrap().orfs.is_empty());
}

#[test]
fn nested_starts_and_result_limits_report_truncation() {
    let m = molecule("ATGATGTAA", Topology::Linear);
    let report = find_orfs(&m, &orf_options()).unwrap();
    assert_eq!(
        report
            .orfs
            .iter()
            .map(|o| o.protein.as_str())
            .collect::<Vec<_>>(),
        ["MM", "M"]
    );
    assert!(!report.truncated);
    let report = find_orfs(
        &m,
        &OrfOptions {
            max_orfs: 1,
            ..orf_options()
        },
    )
    .unwrap();
    assert_eq!(report.orfs.len(), 1);
    assert!(report.truncated);
    assert!(
        find_orfs(
            &m,
            &OrfOptions {
                max_orfs: 0,
                ..orf_options()
            }
        )
        .is_err()
    );
    assert!(
        find_orfs(
            &m,
            &OrfOptions {
                max_orfs: 10_001,
                ..orf_options()
            }
        )
        .is_err()
    );
    assert!(
        find_orfs(&molecule("AA", Topology::Circular), &orf_options())
            .unwrap()
            .orfs
            .is_empty()
    );
}

#[test]
fn circular_complete_orf_boundaries_match_independent_direct_walks() {
    // Exhaustive small canonical strings; reference walks at most one turn per start.
    // It does not use the implementation's cached codons or rolling boundary scan.
    for mut value in 0..6561usize {
        // all 3^8 sequences over A/T/G
        let mut bases = [b'A'; 8];
        for base in &mut bases {
            *base = b"ATG"[value % 3];
            value /= 3;
        }
        let text = std::str::from_utf8(&bases).unwrap();
        let m = molecule(text, Topology::Circular);
        let report = find_orfs(&m, &orf_options()).unwrap();
        let mut expected = Vec::new();
        for start in 0..bases.len() {
            let read = |p: usize| {
                [
                    bases[p % bases.len()],
                    bases[(p + 1) % bases.len()],
                    bases[(p + 2) % bases.len()],
                ]
            };
            if read(start) != *b"ATG" {
                continue;
            }
            for codon in 1..bases.len() / 3 {
                let sequence = read(start + codon * 3);
                if matches!(&sequence, b"TAA" | b"TAG" | b"TGA") {
                    expected.push((start, (codon + 1) * 3));
                    break;
                }
            }
        }
        let actual: Vec<_> = report
            .orfs
            .iter()
            .map(|orf| (orf.oriented_start, orf.nucleotide_length))
            .collect();
        assert_eq!(actual, expected, "{text}");
        for orf in &report.orfs {
            assert_eq!(
                orf.location.extract(&m.sequence).unwrap().len(),
                orf.nucleotide_length
            );
            assert_eq!(orf.protein, "M");
        }
    }
}

#[test]
fn repeated_cds_ranges_are_bounded_before_sequence_materialization() {
    let m = molecule(&"A".repeat(1000), Topology::Linear);
    let feature = cds(Location::Join(vec![Location::range(0, 1000); 10_001]));
    assert!(matches!(translate_cds(&m, &feature), Err(Error::TooLarge)));
}

#[test]
fn ncbi_catalogue_has_all_extant_ids_and_explicit_context_dependent_stop_metadata() {
    let expected = [
        1, 2, 3, 4, 5, 6, 9, 10, 11, 12, 13, 14, 15, 16, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30,
        31, 32, 33,
    ];
    let codes = genetic_codes();
    assert_eq!(codes.iter().map(|c| c.id).collect::<Vec<_>>(), expected);
    for entry in codes {
        assert_eq!(entry, CodeDescriptor::for_id(entry.id).unwrap());
        assert!(!entry.name.is_empty() && entry.name.len() <= 256);
        assert!(entry.start_codons.iter().any(|c| c == "ATG"));
        let dual = match entry.id {
            27 => vec!["TGA"],
            28 => vec!["TAA", "TAG", "TGA"],
            31 => vec!["TAA", "TAG"],
            _ => vec![],
        };
        assert_eq!(entry.context_dependent_stop_codons, dual);
    }
    for absent in [0, 7, 8, 17, 18, 19, 20, 34, 99, 255] {
        assert!(CodeDescriptor::for_id(absent).is_err());
    }
    assert_eq!(
        CodeDescriptor::for_id(3).unwrap().start_codons,
        ["ATA", "ATG", "GTG"]
    );
    assert_eq!(
        CodeDescriptor::for_id(25).unwrap().start_codons,
        ["TTG", "ATG", "GTG"]
    );
}

#[test]
fn added_tables_translate_literal_reassigned_codons_without_silent_standard_fallback() {
    // Literal codon expectations taken from the primary NCBI difference tables;
    // this does not regenerate expected amino acids from the implementation table.
    for (id, text, expected) in [
        (3, "CTTATA", "TM"),
        (4, "TGA", "W"),
        (5, "AGAATA", "SM"),
        (6, "TAATAG", "QQ"),
        (9, "AAAAGA", "NS"),
        (10, "TGA", "C"),
        (12, "CTG", "S"),
        (13, "AGAAGG", "GG"),
        (14, "TAAAAA", "YN"),
        (15, "TAG", "Q"),
        (16, "TAG", "L"),
        (21, "AAAATA", "NM"),
        (22, "TCATAG", "*L"),
        (23, "TTA", "*"),
        (24, "AGAAGG", "SK"),
        (25, "TGA", "G"),
        (26, "CTG", "A"),
        (27, "TAATAGTGA", "QQX"),
        (28, "TAATAGTGA", "XXX"),
        (29, "TAATAG", "YY"),
        (30, "TAATAG", "EE"),
        (31, "TAATAGTGA", "XXW"),
        (32, "TAG", "W"),
        (33, "TAAAGG", "YK"),
    ] {
        let options = TranslationOptions {
            genetic_code: id,
            ..Default::default()
        };
        let report = translate(&dna(text), &options).unwrap();
        assert_eq!(report.protein, expected, "table {id}, {text}");
        let rna = Sequence::new(Alphabet::Rna, &text.replace('T', "U")).unwrap();
        assert_eq!(translate(&rna, &options).unwrap().protein, expected);
        let reverse = dna(text).reverse_complement().unwrap();
        assert_eq!(
            translate(
                &reverse,
                &TranslationOptions {
                    reverse: true,
                    ..options
                }
            )
            .unwrap()
            .protein,
            expected
        );
    }
    let options = TranslationOptions {
        genetic_code: 3,
        initiation: InitiationPolicy::RequireRecognizedStart,
        ..Default::default()
    };
    assert_eq!(
        translate(&dna("GTGCTTTAA"), &options).unwrap().protein,
        "MT*"
    );
    assert!(translate(&dna("ATC"), &options).is_err());
}

#[test]
fn context_dependent_stops_remain_uncertain_in_literal_and_iupac_translation() {
    for (id, text, expected, possible, ambiguous) in [
        (27, "TAATAGTGATGG", "QQXW", 1, 0),
        (28, "TAATAGTGATGG", "XXXW", 3, 0),
        (31, "TAATAGTGATGG", "XXWW", 2, 0),
        (27, "TGR", "X", 1, 1),
        (28, "TAR", "X", 1, 1),
        (31, "TAR", "X", 1, 1),
        (27, "TAR", "Q", 0, 1),
    ] {
        let report = translate(
            &dna(text),
            &TranslationOptions {
                genetic_code: id,
                ..Default::default()
            },
        )
        .unwrap();
        assert_eq!(report.protein, expected);
        assert_eq!(report.possible_stop_codons, possible);
        assert_eq!(report.ambiguous_codons, ambiguous);
        assert_eq!(
            report.unresolved_codons,
            expected.bytes().filter(|b| *b == b'X').count()
        );
        if possible > 0 {
            assert!(
                report
                    .warnings
                    .iter()
                    .any(|s| s.contains("context-dependent termination"))
            );
        }
    }
}

#[test]
fn context_dependent_terminal_codons_do_not_assert_complete_cds_or_complete_orfs() {
    for (id, ending) in [
        (27, "TGA"),
        (28, "TAA"),
        (28, "TAG"),
        (28, "TGA"),
        (31, "TAA"),
        (31, "TAG"),
    ] {
        let m = molecule(&format!("ATGAAA{ending}"), Topology::Linear);
        let mut feature = cds(Location::range(0, 9));
        feature
            .qualifiers
            .push(("transl_table".into(), Some(id.to_string())));
        let report = translate_cds(&m, &feature).unwrap();
        assert_eq!(report.protein, "MKX");
        assert!(!report.has_terminal_stop && !report.complete_under_model);
        assert_eq!(report.translation.possible_stop_codons, 1);
        let mut options = orf_options();
        options.genetic_code = id;
        assert!(find_orfs(&m, &options).unwrap().orfs.is_empty());
        options.incomplete = IncompletePolicy::IncludeAllPartials;
        let partial = find_orfs(&m, &options)
            .unwrap()
            .orfs
            .into_iter()
            .find(|o| o.oriented_start == 0)
            .unwrap();
        assert_eq!(partial.termination, Termination::AmbiguousStop);
        assert_eq!(partial.protein, "MK");
        assert!(!partial.has_stop_codon && partial.partial_3prime);
        assert_eq!(partial.location, Location::range(0, 6));
        let mut circular = m.clone();
        circular.topology = Topology::Circular;
        let circle = find_orfs(&circular, &options).unwrap();
        let partial = circle.orfs.iter().find(|o| o.oriented_start == 0).unwrap();
        assert_eq!(partial.termination, Termination::AmbiguousStop);
        assert!(!partial.has_stop_codon && partial.partial_3prime);
    }
    // New table 6 has definite TGA termination: unlike dual-use tables it can
    // establish the ordinary codon-model completeness predicate.
    let m = molecule("ATGTAATGA", Topology::Linear);
    let mut feature = cds(Location::range(0, 9));
    feature
        .qualifiers
        .push(("transl_table".into(), Some("6".into())));
    let report = translate_cds(&m, &feature).unwrap();
    assert_eq!(report.protein, "MQ");
    assert!(report.has_terminal_stop && report.complete_under_model);
}

#[test]
fn captured_ncbi_catalogue_retains_all_primary_rows_and_provenance() {
    let data = include_str!("../../data/genetic_codes/gc.prt");
    assert!(data.contains("Version 4.6"));
    // Provenance is separate from scientific literal-codon expectations above.
    let manifest: serde_json::Value =
        serde_json::from_str(include_str!("../../data/genetic_codes/provenance.json")).unwrap();
    assert_eq!(manifest["bytes"], data.len());
    let mut ids = Vec::new();
    for block in data.split("\n {\n").skip(1) {
        let id = block
            .lines()
            .find_map(|line| line.trim().strip_prefix("id "))
            .unwrap()
            .split_whitespace()
            .next()
            .unwrap()
            .parse::<u8>()
            .unwrap();
        ids.push(id);
        let row = |key: &str| {
            block
                .lines()
                .find_map(|line| line.trim().strip_prefix(key))
                .unwrap()
                .split('"')
                .nth(1)
                .unwrap()
        };
        let aas = row("ncbieaa").as_bytes();
        let starts = row("sncbieaa").as_bytes();
        assert_eq!(aas.len(), 64);
        assert_eq!(starts.len(), 64);
        let mut expected_starts = Vec::new();
        for (index, &amino) in aas.iter().enumerate() {
            let order = b"TCAG";
            let codon = String::from_utf8(vec![
                order[index / 16],
                order[index / 4 % 4],
                order[index % 4],
            ])
            .unwrap();
            let dual = starts[index] == b'*' && amino != b'*';
            let report = translate(
                &dna(&codon),
                &TranslationOptions {
                    genetic_code: id,
                    ..Default::default()
                },
            )
            .unwrap();
            assert_eq!(
                report.protein.as_bytes(),
                &[if dual { b'X' } else { amino }],
                "table {id}, {codon}"
            );
            assert_eq!(report.possible_stop_codons, usize::from(dual));
            if starts[index] == b'M' {
                expected_starts.push(codon);
            }
        }
        assert_eq!(
            CodeDescriptor::for_id(id).unwrap().start_codons,
            expected_starts
        );
    }
    assert_eq!(
        ids,
        genetic_codes().iter().map(|c| c.id).collect::<Vec<_>>()
    );
}
