//! Independent finite-path oracle and storage/resource owner tests.
//! The oracle enumerates legal moves; it does not implement a DP recurrence.
use super::*;
use crate::alignment::{self, AlignmentOptions};
use collagene_model::Sequence;

fn mask(base: u8) -> u8 {
    match base {
        b'A' => 1,
        b'C' => 2,
        b'G' => 4,
        b'T' => 8,
        b'N' => 15,
        b'R' => 5,
        b'Y' => 10,
        _ => panic!("unsupported oracle symbol"),
    }
}
fn substitution(a: u8, b: u8, s: Scoring) -> i32 {
    let (a, b) = (mask(a), mask(b));
    if a & b == 0 {
        s.mismatch_score
    } else if a.count_ones() == 1 && b.count_ones() == 1 {
        s.match_score
    } else {
        s.ambiguous_score
    }
}
struct Oracle<'a> {
    a: &'a [u8],
    b: &'a [u8],
    s: Scoring,
    lower: i64,
    upper: i64,
    best: i32,
    outside: Option<i32>,
}
impl Oracle<'_> {
    fn visit(&mut self, i: usize, j: usize, last: u8, score: i32, escaped: bool) {
        let d = i as i64 - j as i64;
        let escaped = escaped || d < self.lower || d > self.upper;
        if i == self.a.len() && j == self.b.len() {
            self.best = self.best.max(score);
            if escaped {
                self.outside = Some(self.outside.map_or(score, |best| best.max(score)));
            }
            return;
        }
        if i < self.a.len() && j < self.b.len() {
            self.visit(
                i + 1,
                j + 1,
                0,
                score + substitution(self.a[i], self.b[j], self.s),
                escaped,
            );
        }
        if j < self.b.len() {
            let cost = if last == 1 {
                self.s.gap_extend
            } else {
                self.s.gap_open
            };
            self.visit(i, j + 1, 1, score + cost, escaped);
        }
        if i < self.a.len() {
            let cost = if last == 2 {
                self.s.gap_extend
            } else {
                self.s.gap_open
            };
            self.visit(i + 1, j, 2, score + cost, escaped);
        }
    }
}
fn oracle(a: &str, b: &str, s: Scoring, radius: usize) -> (i32, Option<i32>) {
    let delta = a.len() as i64 - b.len() as i64;
    let mut oracle = Oracle {
        a: a.as_bytes(),
        b: b.as_bytes(),
        s,
        lower: delta.min(0) - radius as i64,
        upper: delta.max(0) + radius as i64,
        best: i32::MIN,
        outside: None,
    };
    oracle.visit(0, 0, 3, 0, false);
    (oracle.best, oracle.outside)
}
fn words(length: usize) -> Vec<String> {
    (0..1usize << length)
        .map(|bits| {
            (0..length)
                .map(|i| if bits & (1 << i) == 0 { 'A' } else { 'C' })
                .collect()
        })
        .collect()
}
fn scorings() -> [Scoring; 3] {
    [
        Scoring::default(),
        Scoring {
            match_score: 3,
            ambiguous_score: -2,
            mismatch_score: -9,
            gap_open: -1,
            gap_extend: -7,
        },
        Scoring {
            match_score: 4,
            ambiguous_score: 2,
            mismatch_score: -50,
            gap_open: -6,
            gap_extend: -2,
        },
    ]
}
#[test]
fn exhaustive_tiny_legal_paths_bound_every_actual_omitted_path() {
    for s in scorings() {
        for n in 0..=3 {
            for m in 0..=3 {
                for a in words(n) {
                    for b in words(m) {
                        for radius in 0..=n.min(m) + 1 {
                            let (_, outside) = oracle(&a, &b, s, radius);
                            let bound =
                                Domain::Corridor { radius }.off_band_bound(n, m, s).unwrap();
                            match (bound, outside) {
                                (OffBandBound::Finite(upper), Some(actual)) => assert!(
                                    i64::from(actual) <= upper,
                                    "{a}/{b} radius{radius}: {actual}>{upper}"
                                ),
                                (OffBandBound::NoEscapingPath, Some(_)) => panic!(
                                    "An actual omitted path contradicts no-escaping certificate"
                                ),
                                (_, None) => {}
                            }
                        }
                    }
                }
            }
        }
    }
}
#[test]
fn full_and_adaptive_current_kernel_scores_match_independent_path_enumeration() {
    // Both policies use the same native recurrence; the independent oracle
    // enumerates paths without sharing production storage or transitions.
    for s in scorings() {
        for n in 0..=3 {
            for m in 0..=3 {
                for a in words(n) {
                    for b in words(m) {
                        let expected = oracle(&a, &b, s, 0).0;
                        let reference = Sequence::dna(&a);
                        let query = Sequence::dna(&b);
                        if a.is_empty() || b.is_empty() {
                            // Empty endpoints remain in the independent geometric
                            // proof oracle, but the native sequence model refuses
                            // them before any alignment kernel can be called.
                            assert_eq!(
                                matches!(&reference, Err(Error::EmptySequence)),
                                a.is_empty()
                            );
                            assert_eq!(matches!(&query, Err(Error::EmptySequence)), b.is_empty());
                            continue;
                        }
                        let reference = reference.unwrap();
                        let query = query.unwrap();
                        let mut options = AlignmentOptions {
                            mode: AlignmentMode::Global,
                            scoring: s,
                            max_cells: 100_000,
                            global_storage: GlobalStorage::Full,
                        };
                        let full = alignment::align(&reference, &query, &options).unwrap();
                        assert_eq!(full.score, expected);
                        assert!(full.search.unrestricted_optimum_proven);
                        options.global_storage = GlobalStorage::Adaptive {
                            initial_radius: 0,
                            maximum_radius: n.min(m),
                            max_attempts: 8,
                        };
                        let adaptive = alignment::align(&reference, &query, &options).unwrap();
                        assert_eq!(adaptive.score, expected);
                        assert!(adaptive.search.unrestricted_optimum_proven);
                        assert_eq!(adaptive.cigar, full.cigar);
                        options.global_storage = GlobalStorage::Adaptive {
                            initial_radius: 0,
                            maximum_radius: 0,
                            max_attempts: 1,
                        };
                        let capped = alignment::align(&reference, &query, &options).unwrap();
                        assert!(capped.score <= expected);
                        if capped.search.unrestricted_optimum_proven {
                            assert_eq!(capped.score, expected);
                        }
                    }
                }
            }
        }
    }
}
#[test]
fn independent_ambiguity_and_affine_gap_opening_cases_cover_nondefault_scores() {
    for s in scorings() {
        for (a, b) in [("ANR", "AGG"), ("A", "TT"), ("AA", "C"), ("AC", "CA")] {
            let expected = oracle(a, b, s, 0).0;
            let result = alignment::align(
                &Sequence::dna(a).unwrap(),
                &Sequence::dna(b).unwrap(),
                &AlignmentOptions {
                    mode: AlignmentMode::Global,
                    scoring: s,
                    max_cells: 10_000,
                    global_storage: GlobalStorage::Adaptive {
                        initial_radius: 0,
                        maximum_radius: 3,
                        max_attempts: 8,
                    },
                },
            )
            .unwrap();
            assert_eq!(result.score, expected);
            assert!(result.search.unrestricted_optimum_proven);
        }
    }
}
#[test]
fn equality_retains_an_actual_off_band_tie_as_unproven() {
    let scoring = Scoring {
        match_score: 2,
        ambiguous_score: 0,
        mismatch_score: -2,
        gap_open: -1,
        gap_extend: -1,
    };
    let (best, outside) = oracle("A", "T", scoring, 0);
    assert_eq!((best, outside), (-2, Some(-2)));
    let attempt = Attempt::completed(Domain::Corridor { radius: 0 }, 1, 1, scoring, -2, 2).unwrap();
    assert_eq!(attempt.off_band_bound, OffBandBound::Finite(-2));
    assert!(!attempt.unrestricted_optimum_proven);
    assert!(
        Attempt::completed(Domain::Corridor { radius: 1 }, 1, 1, scoring, -2, 4)
            .unwrap()
            .unrestricted_optimum_proven
    );
    assert!(Attempt::completed(Domain::Full, 1, 1, scoring, NEG, 4).is_err());
}
#[test]
fn unequal_endpoint_corridors_store_exact_original_j_coordinates() {
    for (n, m) in [(8, 3), (3, 8), (0, 4), (4, 0), (4, 4)] {
        for radius in 0..=3 {
            let domain = Domain::Corridor { radius };
            let layout = RowLayout::new(n, m, domain, 10_000).unwrap();
            let delta = n as i64 - m as i64;
            let mut independently_counted = 0;
            for i in 0..=n {
                for j in 0..=m {
                    let d = i as i64 - j as i64;
                    let admitted =
                        delta.min(0) - radius as i64 <= d && d <= delta.max(0) + radius as i64;
                    if admitted {
                        assert_eq!(layout.index(i, j).unwrap(), independently_counted);
                        independently_counted += 1;
                    } else {
                        assert!(layout.index(i, j).is_err());
                    }
                }
            }
            assert_eq!(layout.cells(), independently_counted);
            assert!(layout.index(0, 0).is_ok() && layout.index(n, m).is_ok());
            let mut grid = Grid::new(layout).unwrap();
            grid.begin(0).unwrap();
            grid.rows().1.set(0, [0, NEG, NEG]).unwrap();
            grid.finish_row();
            assert_eq!(grid.final_row().get(0), [0, NEG, NEG]);
            grid.set_predecessor(n, m, 0, 2).unwrap();
            assert_eq!(grid.predecessor(n, m, 0).unwrap(), 2);
            assert!(grid.predecessor(n + 1, m, 0).is_err());
            assert!(grid.layout().bytes() <= MAX_MATRIX_BYTES);
        }
    }
}
#[test]
fn omitted_neighbors_are_unreachable_and_never_restart_or_wrap_a_score_row() {
    let mut grid =
        Grid::new(RowLayout::new(4, 4, Domain::Corridor { radius: 0 }, 5).unwrap()).unwrap();
    grid.begin(0).unwrap();
    grid.rows().1.set(0, [0, NEG, NEG]).unwrap();
    grid.finish_row();
    grid.begin(1).unwrap();
    let (previous, current) = grid.rows();
    assert_eq!(previous.get(1), [NEG; 3]);
    assert_eq!(current.get(0), [NEG; 3]);
    assert!(current.set(0, [0; 3]).is_err());
    assert!(grid.set_predecessor(1, 0, 0, 0).is_err());
}
#[test]
fn retries_and_new_candidates_share_actual_cells_without_refunds_or_reset() {
    let mut shared = CellWork::new(20).unwrap();
    {
        let mut alignment = CellBudget::new(&mut shared, 15).unwrap();
        let first =
            RowLayout::new(4, 4, Domain::Corridor { radius: 0 }, alignment.remaining()).unwrap();
        alignment.admit(first.cells()).unwrap();
        for i in 0..=4 {
            alignment.charge_row(first.row_width(i).unwrap()).unwrap();
        }
        assert_eq!(alignment.used(), 5);
        let second = RowLayout::new(4, 4, Domain::Corridor { radius: 1 }, 10);
        assert!(second.is_err());
        assert_eq!(alignment.used(), 5);
    }
    assert_eq!(shared.used(), 5);
    let mut next = CellBudget::new(&mut shared, 20).unwrap();
    next.charge_row(15).unwrap();
    assert_eq!(next.remaining(), 0);
    assert!(next.charge_row(1).is_err());
    assert_eq!(next.used(), 15);
    assert_eq!(shared.used(), 20);
    assert_eq!(shared.remaining(), 0);
}
#[test]
fn layout_and_proof_reject_arithmetic_and_storage_amplification_before_allocating() {
    assert!(RowLayout::new(usize::MAX, 1, Domain::Full, 100).is_err());
    assert!(RowLayout::new(20_000, 20_000, Domain::Full, 64_000_000).is_err());
    assert!(
        Domain::Corridor { radius: usize::MAX }
            .off_band_bound(1, 1, Scoring::default())
            .is_err()
    );
    let mut invalid = Scoring::default();
    invalid.gap_extend = 0;
    assert!(
        Domain::Corridor { radius: 0 }
            .off_band_bound(1, 1, invalid)
            .is_err()
    );
    assert!(
        GlobalStorage::default()
            .validate(AlignmentMode::Local)
            .is_err()
    );
    assert!(GlobalStorage::Full.validate(AlignmentMode::Local).is_ok());
}

#[test]
fn evidence_validator_rejects_forged_cells_score_certificate_and_attempt_order() {
    let scoring = Scoring {
        match_score: 2,
        ambiguous_score: 0,
        mismatch_score: -2,
        gap_open: -1,
        gap_extend: -1,
    };
    let first = Attempt::completed(Domain::Corridor { radius: 0 }, 1, 1, scoring, -2, 2).unwrap();
    let second = Attempt::completed(Domain::Corridor { radius: 1 }, 1, 1, scoring, -2, 4).unwrap();
    let evidence = SearchEvidence {
        storage: GlobalStorage::Adaptive {
            initial_radius: 0,
            maximum_radius: 1,
            max_attempts: 2,
        },
        attempts: vec![first, second],
        evaluated_cells: 6,
        unrestricted_optimum_proven: true,
        stop_reason: StopReason::Certified,
    };
    evidence
        .validate(1, 1, AlignmentMode::Global, scoring, -2, 6)
        .unwrap();
    let mut changed = evidence.clone();
    changed.attempts[0].evaluated_cells = 1;
    assert!(
        changed
            .validate(1, 1, AlignmentMode::Global, scoring, -2, 6)
            .is_err()
    );
    let mut changed = evidence.clone();
    changed.attempts[0].unrestricted_optimum_proven = true;
    assert!(
        changed
            .validate(1, 1, AlignmentMode::Global, scoring, -2, 6)
            .is_err()
    );
    let mut changed = evidence.clone();
    changed.attempts.swap(0, 1);
    assert!(
        changed
            .validate(1, 1, AlignmentMode::Global, scoring, -2, 6)
            .is_err()
    );
    assert!(
        evidence
            .validate(1, 1, AlignmentMode::Global, scoring, -1, 6)
            .is_err()
    );
    assert!(
        evidence
            .validate(1, 1, AlignmentMode::Global, scoring, -2, 5)
            .is_err()
    );
}
