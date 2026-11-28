//! Hand-specified genetic-code/coordinate fixtures, not outputs of another cloning application.
use collagene_core::translation::*;
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
                genetic_code: 14,
                ..TranslationOptions::default()
            }
        )
        .is_err()
    );
    assert_eq!(
        genetic_codes().iter().map(|c| c.id).collect::<Vec<_>>(),
        [1, 2, 11]
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
