use collagene_core::{alignment::*, *};
fn dna(s: &str) -> Sequence {
    Sequence::dna(s).unwrap()
}
fn options(mode: AlignmentMode) -> AlignmentOptions {
    AlignmentOptions {
        mode,
        ..AlignmentOptions::default()
    }
}

#[test]
fn independent_global_gap_fixture_uses_one_open_plus_remaining_extensions() {
    let a = align(
        &dna("AAAACCCCGGGG"),
        &dna("AAAAGGGG"),
        &AlignmentOptions::default(),
    )
    .unwrap();
    assert_eq!(a.score, 9);
    assert_eq!(a.cigar, "4=4D4=");
    assert_eq!(a.reference_aligned, "AAAACCCCGGGG");
    assert_eq!(a.query_aligned, "AAAA----GGGG");
    assert_eq!(a.differences.len(), 1);
    let d = &a.differences[0];
    assert_eq!(
        (
            d.kind,
            d.reference_start,
            d.reference_end,
            d.query_start,
            d.query_end
        ),
        (DifferenceKind::Deletion, 4, 8, 4, 4)
    );
    assert_eq!(a.stats.paired_bases, 8);
    assert_eq!(a.stats.reference_coverage, 8. / 12.);
    assert_eq!(a.stats.query_coverage, 1.);
    let b = align(
        &dna("AAAAGGGG"),
        &dna("AAAACCCCGGGG"),
        &AlignmentOptions::default(),
    )
    .unwrap();
    assert_eq!(b.score, 9);
    assert_eq!(b.cigar, "4=4I4=");
    assert_eq!(b.differences[0].kind, DifferenceKind::Insertion);
}
#[test]
fn local_and_query_end_to_end_modes_report_unaligned_flanks_without_false_coverage() {
    let a = align(
        &dna("TTTACGTGGG"),
        &dna("CCACGTAAA"),
        &options(AlignmentMode::Local),
    )
    .unwrap();
    assert_eq!(a.score, 8);
    assert_eq!(a.cigar, "4=");
    assert_eq!(
        (
            a.reference_start,
            a.reference_end,
            a.query_start,
            a.query_end
        ),
        (3, 7, 2, 6)
    );
    assert_eq!(a.stats.reference_coverage, 0.4);
    assert_eq!(a.stats.query_coverage, 4. / 9.);
    let b = align(
        &dna("TTTACGTGGG"),
        &dna("ACGT"),
        &options(AlignmentMode::QueryInReference),
    )
    .unwrap();
    assert_eq!(b.score, 8);
    assert_eq!(
        (
            b.reference_start,
            b.reference_end,
            b.query_start,
            b.query_end
        ),
        (3, 7, 0, 4)
    );
    assert_eq!(b.cigar, "4=");
    assert_eq!(b.stats.query_coverage, 1.);
    let empty = align(&dna("AAAA"), &dna("CCCC"), &options(AlignmentMode::Local)).unwrap();
    assert_eq!(empty.score, 0);
    assert!(empty.cigar.is_empty());
    assert_eq!(empty.stats.reference_coverage, 0.);
    assert_eq!(empty.stats.certain_identity, None);
    assert!(empty.warnings.iter().any(|s| s.contains("No positive")));
}
#[test]
fn compatible_ambiguity_is_distinct_from_known_identity_and_disjoint_mismatches() {
    let a = align(&dna("ACRN"), &dna("ACAT"), &AlignmentOptions::default()).unwrap();
    assert_eq!(a.cigar, "2=2M");
    assert_eq!(a.score, 4);
    assert_eq!(a.stats.certain_matches, 2);
    assert_eq!(a.stats.ambiguous_compatible, 2);
    assert_eq!(a.stats.mismatches, 0);
    assert_eq!(a.stats.certain_identity, Some(0.5));
    assert_eq!(a.stats.compatible_identity, Some(1.));
    assert_eq!(a.differences[0].kind, DifferenceKind::AmbiguousCompatible);
    let b = align(&dna("RY"), &dna("CA"), &AlignmentOptions::default()).unwrap();
    assert_eq!(b.cigar, "2X");
    assert_eq!(b.score, -6);
    assert_eq!(b.stats.mismatches, 2);
    let rna = Sequence::new(Alphabet::Rna, "ACGU").unwrap();
    assert_eq!(
        align(&rna, &rna, &AlignmentOptions::default())
            .unwrap()
            .cigar,
        "4="
    );
    assert!(align(&dna("ACGT"), &rna, &AlignmentOptions::default()).is_err());
}

// Exhaustive path enumeration, not the affine DP recurrence used by the implementation.
fn enumerate(a: &[u8], b: &[u8], previous: u8, s: Scoring) -> i32 {
    if a.is_empty() && b.is_empty() {
        return 0;
    }
    let mut best = i32::MIN;
    if !a.is_empty() && !b.is_empty() {
        best = best.max(
            (if a[0] == b[0] {
                s.match_score
            } else {
                s.mismatch_score
            }) + enumerate(&a[1..], &b[1..], 0, s),
        );
    }
    if !b.is_empty() {
        best = best.max(
            (if previous == 1 {
                s.gap_extend
            } else {
                s.gap_open
            }) + enumerate(a, &b[1..], 1, s),
        );
    }
    if !a.is_empty() {
        best = best.max(
            (if previous == 2 {
                s.gap_extend
            } else {
                s.gap_open
            }) + enumerate(&a[1..], b, 2, s),
        );
    }
    best
}
fn rescore(a: &str, b: &str, s: Scoring) -> i32 {
    let mut total = 0;
    let mut previous = 0;
    for (a, b) in a.bytes().zip(b.bytes()) {
        let state = if a == b'-' {
            1
        } else if b == b'-' {
            2
        } else {
            0
        };
        total += match state {
            1 | 2 => {
                if state == previous {
                    s.gap_extend
                } else {
                    s.gap_open
                }
            }
            _ => {
                if a == b {
                    s.match_score
                } else {
                    s.mismatch_score
                }
            }
        };
        previous = state;
    }
    total
}
fn words() -> Vec<String> {
    (1..=4)
        .flat_map(|n| {
            (0..(1 << n)).map(move |bits| {
                (0..n)
                    .map(|i| if bits & (1 << i) == 0 { 'A' } else { 'C' })
                    .collect()
            })
        })
        .collect()
}
#[test]
fn affine_global_scores_match_exhaustive_alignment_paths_including_opposing_gaps() {
    for scoring in [
        Scoring::default(),
        Scoring {
            match_score: 3,
            ambiguous_score: 0,
            mismatch_score: -9,
            gap_open: -2,
            gap_extend: -1,
        },
    ] {
        let options = AlignmentOptions {
            scoring,
            ..AlignmentOptions::default()
        };
        for a in words() {
            for b in words() {
                let result = align(&dna(&a), &dna(&b), &options).unwrap();
                assert_eq!(
                    result.score,
                    enumerate(a.as_bytes(), b.as_bytes(), 0, scoring),
                    "{a} vs {b}"
                );
                assert_eq!(
                    result.score,
                    rescore(&result.reference_aligned, &result.query_aligned, scoring),
                    "traceback {a}/{b}"
                );
                assert_eq!(result.reference_aligned.replace('-', ""), a);
                assert_eq!(result.query_aligned.replace('-', ""), b);
            }
        }
    }
}
#[test]
fn local_scores_match_exhaustive_substring_global_paths() {
    let scoring = Scoring::default();
    for a in words().into_iter().filter(|s| s.len() <= 3) {
        for b in words().into_iter().filter(|s| s.len() <= 3) {
            let mut expected = 0;
            for rs in 0..a.len() {
                for re in rs + 1..=a.len() {
                    for qs in 0..b.len() {
                        for qe in qs + 1..=b.len() {
                            expected = expected.max(enumerate(
                                &a.as_bytes()[rs..re],
                                &b.as_bytes()[qs..qe],
                                0,
                                scoring,
                            ));
                        }
                    }
                }
            }
            let aligned = align(&dna(&a), &dna(&b), &options(AlignmentMode::Local)).unwrap();
            assert_eq!(aligned.score, expected, "{a}/{b}");
            assert_eq!(
                aligned.score,
                rescore(&aligned.reference_aligned, &aligned.query_aligned, scoring)
            );
            assert_eq!(
                aligned.reference_aligned.replace('-', ""),
                a[aligned.reference_start..aligned.reference_end]
            );
            assert_eq!(
                aligned.query_aligned.replace('-', ""),
                b[aligned.query_start..aligned.query_end]
            );
        }
    }
}
#[test]
fn query_in_reference_scores_match_exhaustive_reference_substring_paths() {
    let scoring = Scoring::default();
    for a in words().into_iter().filter(|s| s.len() <= 3) {
        for b in words().into_iter().filter(|s| s.len() <= 3) {
            let mut expected = i32::MIN;
            for rs in 0..=a.len() {
                for re in rs..=a.len() {
                    expected =
                        expected.max(enumerate(&a.as_bytes()[rs..re], b.as_bytes(), 0, scoring));
                }
            }
            let aligned = align(
                &dna(&a),
                &dna(&b),
                &options(AlignmentMode::QueryInReference),
            )
            .unwrap();
            assert_eq!(aligned.score, expected, "{a}/{b}");
            assert_eq!(
                aligned.score,
                rescore(&aligned.reference_aligned, &aligned.query_aligned, scoring)
            );
            assert_eq!(
                aligned.reference_aligned.replace('-', ""),
                a[aligned.reference_start..aligned.reference_end]
            );
            assert_eq!(
                aligned.query_aligned.replace('-', ""),
                b[aligned.query_start..aligned.query_end]
            );
            assert_eq!(aligned.query_start, 0);
            assert_eq!(aligned.query_end, b.len());
        }
    }
}
#[test]
fn invalid_scoring_and_matrix_budgets_fail_before_large_allocations() {
    let mut o = AlignmentOptions::default();
    o.max_cells = 3;
    assert!(align(&dna("ACGT"), &dna("AGGT"), &o).is_err());
    o.max_cells = HARD_MAX_CELLS + 1;
    assert!(align(&dna("A"), &dna("A"), &o).is_err());
    o = AlignmentOptions::default();
    o.scoring.gap_extend = 0;
    assert!(align(&dna("A"), &dna("A"), &o).is_err());
    o.scoring.gap_extend = -1;
    o.scoring.match_score = i32::MAX;
    assert!(align(&dna("A"), &dna("A"), &o).is_err());
    let skewed = dna(&"A".repeat(4_000_000));
    assert!(
        align(
            &dna("A"),
            &skewed,
            &options(AlignmentMode::QueryInReference)
        )
        .unwrap_err()
        .to_string()
        .contains("allocation budget")
    );
}

#[test]
fn circular_identical_rotation_and_reverse_complement_have_explicit_source_mapping() {
    let reference = Molecule::dna("reference", "ACGATTCG", Topology::Circular).unwrap();
    let query = Molecule::dna("query", "ATTCGACG", Topology::Circular).unwrap();
    let result = compare_molecules(&reference, &query, &ComparisonOptions::default()).unwrap();
    assert_eq!(result.reference_origin, 3);
    assert!(!result.query_reversed);
    assert_eq!(result.alignment.cigar, "8=");
    assert!(result.optimality_proven);
    assert!(!result.exhaustive_reference_origins);
    let reverse = Molecule::dna("reverse", "CGTCGAAT", Topology::Circular).unwrap();
    let result = compare_molecules(&reference, &reverse, &ComparisonOptions::default()).unwrap();
    assert_eq!(result.reference_origin, 3);
    assert!(result.query_reversed);
    assert_eq!(result.alignment.score, 16);
    assert_eq!(result.query_revision, reverse.revision());
    assert_eq!(result.reference_revision, reference.revision());
    let mut options = ComparisonOptions::default();
    options.both_query_strands = false;
    options.reference_rotations = false;
    let fixed = compare_molecules(&reference, &query, &options).unwrap();
    assert_eq!(fixed.reference_origin, 0);
    assert!(fixed.alignment.score < 16);
}
#[test]
fn small_circular_search_matches_independently_enumerated_rotation_scores() {
    let reference = Molecule::dna("ref", "AACG", Topology::Circular).unwrap();
    let query = Molecule::dna("query", "CGAT", Topology::Circular).unwrap();
    let result = compare_molecules(&reference, &query, &ComparisonOptions::default()).unwrap();
    let mut expected = i32::MIN;
    for oriented in [
        query.sequence.clone(),
        query.sequence.reverse_complement().unwrap(),
    ] {
        for origin in 0..reference.sequence.len() {
            let r = format!(
                "{}{}",
                &reference.sequence.as_str()[origin..],
                &reference.sequence.as_str()[..origin]
            );
            expected = expected.max(enumerate(
                r.as_bytes(),
                oriented.as_str().as_bytes(),
                0,
                Scoring::default(),
            ));
        }
    }
    assert_eq!(result.alignment.score, expected);
    assert!(result.exhaustive_reference_origins);
    assert!(result.exhaustive_query_orientations);
    assert!(result.optimality_proven);
    assert_eq!(result.evaluated_candidates, 8);
}
#[test]
fn heuristic_origin_and_total_budget_limits_never_claim_exhaustiveness() {
    let reference = Molecule::dna("ref", "ACGATTCG", Topology::Circular).unwrap();
    let query = Molecule::dna("query", "ATTCGACA", Topology::Circular).unwrap();
    let mut options = ComparisonOptions::default();
    options.max_origins = 1;
    options.max_total_cells = 81;
    let result = compare_molecules(&reference, &query, &options).unwrap();
    assert_eq!(result.evaluated_candidates, 1);
    assert!(!result.exhaustive_reference_origins);
    assert!(!result.exhaustive_query_orientations);
    assert!(!result.optimality_proven);
    assert!(result.warnings.iter().any(|s| s.contains("incomplete")));
    options.max_total_cells = 80;
    assert!(compare_molecules(&reference, &query, &options).is_err());
    let bases = "AAACGTCAGTTCGGACCTAGGCATTCACGTGA";
    let reference = Molecule::dna("larger", bases, Topology::Circular).unwrap();
    let mut edited = format!("{}{}", &bases[7..], &bases[..7]);
    let last = edited.len() - 1;
    let replacement = if &edited[last..] == "A" { "C" } else { "A" };
    edited.replace_range(last.., replacement);
    let query = Molecule::dna("edited", &edited, Topology::Circular).unwrap();
    let result = compare_molecules(
        &reference,
        &query,
        &ComparisonOptions {
            max_origins: 2,
            ..ComparisonOptions::default()
        },
    )
    .unwrap();
    assert_eq!(result.reference_origin, 7);
    assert_eq!(result.alignment.stats.mismatches, 1);
    assert!(!result.optimality_proven);
}

#[test]
fn identical_global_input_uses_a_proven_linear_shortcut_without_hiding_ambiguity() {
    let sequence = dna(&"ACGT".repeat(12_500));
    let options = AlignmentOptions {
        max_cells: 1,
        ..AlignmentOptions::default()
    };
    let result = align(&sequence, &sequence, &options).unwrap();
    assert_eq!(result.score, 100_000);
    assert_eq!(result.cigar, "50000=");
    let uncertain = dna("ACNN");
    let result = align(&uncertain, &uncertain, &options).unwrap();
    assert_eq!(result.cigar, "2=2M");
    assert_eq!(result.stats.ambiguous_compatible, 2);
    let low = Scoring {
        ambiguous_score: -3,
        mismatch_score: -3,
        gap_open: -1,
        gap_extend: -1,
        ..Scoring::default()
    };
    let result = align(
        &dna("N"),
        &dna("N"),
        &AlignmentOptions {
            scoring: low,
            ..AlignmentOptions::default()
        },
    )
    .unwrap();
    assert_eq!(result.score, -2);
    assert_eq!(result.stats.paired_bases, 0);
}
