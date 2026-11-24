use collagene_view_layout::*;
use std::collections::BTreeMap;

fn item(id: u64, spans: &[(u64, u64)]) -> Item {
    Item {
        id,
        spans: spans
            .iter()
            .map(|&(start, end)| Interval { start, end })
            .collect(),
    }
}

fn options(lanes: usize) -> Limits {
    Limits {
        max_lanes: lanes,
        ..Limits::default()
    }
}

#[test]
fn literal_adjacent_nested_joined_and_overflow() {
    let inputs = [
        item(20, &[(4, 8)]),
        item(10, &[(0, 4)]),
        item(30, &[(1, 3)]),
        item(40, &[(7, 9), (0, 1)]),
        item(50, &[(0, 10)]),
    ];
    let result = pack(&inputs, 10, &options(2)).unwrap();
    assert_eq!(
        result
            .assignments
            .iter()
            .map(|a| (a.id, a.placement))
            .collect::<Vec<_>>(),
        vec![
            (20, Placement::Lane(1)),
            (10, Placement::Lane(1)),
            (30, Placement::Lane(0)),
            (40, Placement::Lane(0)),
            (50, Placement::Overflow)
        ]
    );
    assert_eq!(result.lane_count, 2);
    assert_eq!(result.overflow_ids().collect::<Vec<_>>(), vec![50]);
    assert_eq!(
        result.assignments[3].spans,
        vec![Interval { start: 0, end: 1 }, Interval { start: 7, end: 9 }]
    );
}

#[test]
fn own_spans_normalize_without_envelope_filling() {
    let input = [
        item(9, &[(8, 10), (2, 4), (0, 2), (3, 5), (8, 9)]),
        item(2, &[(5, 8)]),
    ];
    let result = pack(&input, 10, &options(1)).unwrap();
    assert_eq!(
        result.assignments[0].spans,
        vec![
            Interval { start: 0, end: 5 },
            Interval { start: 8, end: 10 }
        ]
    );
    assert!(
        result
            .assignments
            .iter()
            .all(|a| a.placement == Placement::Lane(0))
    );
    assert_eq!(result.lane_count, 1);
}

#[test]
fn circular_split_keeps_one_lane_and_gap_is_available() {
    let input = [
        item(1, &[(90, 100), (0, 10)]),
        item(2, &[(10, 90)]),
        item(3, &[(95, 100)]),
    ];
    let result = pack(&input, 100, &options(1)).unwrap();
    assert_eq!(
        result
            .assignments
            .iter()
            .map(|a| a.placement)
            .collect::<Vec<_>>(),
        vec![Placement::Lane(0), Placement::Lane(0), Placement::Overflow]
    );
}

#[test]
fn equal_geometry_id_tie_break_does_not_follow_input_order() {
    let a = [item(99, &[(1, 4)]), item(3, &[(1, 4)])];
    let b = [a[1].clone(), a[0].clone()];
    assert_eq!(
        pack(&a, 4, &options(1))
            .unwrap()
            .overflow_ids()
            .collect::<Vec<_>>(),
        vec![99]
    );
    assert_eq!(
        pack(&b, 4, &options(1))
            .unwrap()
            .overflow_ids()
            .collect::<Vec<_>>(),
        vec![99]
    );
}

#[test]
fn empty_zero_lane_and_u64_boundary() {
    assert_eq!(pack(&[], 0, &options(0)).unwrap().lane_count, 0);
    let input = [item(u64::MAX, &[(u64::MAX - 2, u64::MAX)])];
    assert_eq!(
        pack(&input, u64::MAX, &options(0))
            .unwrap()
            .overflow_ids()
            .collect::<Vec<_>>(),
        vec![u64::MAX]
    );
    assert_eq!(pack(&input, u64::MAX, &options(1)).unwrap().lane_count, 1);
    assert_eq!(
        Interval::from_start_length(u64::MAX, 1),
        Err(Error::CoordinateOverflow)
    );
    assert_eq!(Interval::from_start_length(0, 0), Err(Error::EmptyInterval));
    assert_eq!(
        Interval::from_start_length(u64::MAX - 1, 1).unwrap().end,
        u64::MAX
    );
}

#[test]
fn rejects_duplicate_missing_empty_reversed_outside_source() {
    assert_eq!(
        pack(&[item(1, &[(0, 1)]), item(1, &[(2, 3)])], 3, &options(1)),
        Err(Error::DuplicateId(1))
    );
    assert_eq!(
        pack(&[item(1, &[])], 3, &options(1)),
        Err(Error::MissingIntervals(1))
    );
    for spans in [
        &[(0, 0)][..],
        &[(3, 2)][..],
        &[(1, 4)][..],
        &[(0, u64::MAX)][..],
    ] {
        assert_eq!(
            pack(&[item(1, spans)], 3, &options(1)),
            Err(Error::InvalidInterval { id: 1, index: 0 })
        );
    }
    assert!(pack(&[item(1, &[(0, 1)])], 0, &options(1)).is_err());
}

#[test]
fn limits_and_exact_work_budget_are_checked_before_input_validation() {
    let input = [item(1, &[(0, 1)]), item(2, &[(1, 2)])];
    let mut limits = options(1);
    let bound = work_bound(&input, &limits).unwrap();
    limits.max_work = bound;
    assert_eq!(pack(&input, 2, &limits).unwrap().work_bound, bound);
    limits.max_work = bound - 1;
    assert_eq!(
        pack(&input, 2, &limits),
        Err(Error::WorkLimit {
            required: bound,
            allowed: bound - 1
        })
    );
    limits = options(1);
    limits.max_items = 1;
    assert_eq!(pack(&input, 2, &limits), Err(Error::ItemLimit));
    limits = options(1);
    limits.max_intervals = 1;
    assert_eq!(pack(&input, 2, &limits), Err(Error::IntervalLimit));
    for invalid in [
        Limits {
            max_items: MAX_ITEMS + 1,
            ..options(1)
        },
        Limits {
            max_intervals: MAX_INTERVALS + 1,
            ..options(1)
        },
        Limits {
            max_lanes: MAX_LANES + 1,
            ..options(1)
        },
        Limits {
            max_work: MAX_WORK + 1,
            ..options(1)
        },
    ] {
        assert_eq!(pack(&[], 0, &invalid), Err(Error::InvalidLimits));
    }
    // A request cap must fail before duplicate-ID / span validation work.
    assert_eq!(
        pack(
            &[item(1, &[]), item(1, &[])],
            0,
            &Limits {
                max_items: 0,
                ..options(1)
            }
        ),
        Err(Error::ItemLimit)
    );
}

#[test]
fn many_overlaps_are_retained_under_cap() {
    let input: Vec<_> = (0..1000).rev().map(|id| item(id, &[(0, 1)])).collect();
    let result = pack(&input, 1, &options(3)).unwrap();
    assert_eq!(result.assignments.len(), 1000);
    assert_eq!(result.lane_count, 3);
    assert_eq!(result.overflow_ids().count(), 997);
    assert_eq!(
        result.assignments.last().unwrap().placement,
        Placement::Lane(0)
    );
}

#[test]
fn ten_thousand_sparse_annotations_on_one_megabase_use_one_lane() {
    let inputs: Vec<_> = (0..10_000)
        .rev()
        .map(|id| item(id, &[(id * 100 + 10, id * 100 + 30)]))
        .collect();
    let result = pack(&inputs, 1_000_000, &Limits::default()).unwrap();
    assert_eq!(result.assignments.len(), 10_000);
    assert_eq!(result.lane_count, 1);
    assert_eq!(result.overflow_ids().count(), 0);
    for (assignment, original) in result.assignments.iter().zip(&inputs) {
        assert_eq!(assignment.id, original.id);
        assert_eq!(assignment.spans, original.spans);
        assert_eq!(assignment.placement, Placement::Lane(0));
    }
    assert!(result.work_bound <= Limits::default().max_work);
    assert!(memory_charge(&inputs, &Limits::default()).unwrap() < 4 * 1024 * 1024);
}

#[test]
fn ten_thousand_mixed_annotations_keep_joined_overflow_and_all_ids() {
    let mut inputs = Vec::new();
    for group in 0..2500 {
        let b = group * 20;
        inputs.extend([
            item(group * 4, &[(b, b + 10)]),
            item(group * 4 + 1, &[(b + 2, b + 8)]),
            item(group * 4 + 2, &[(b + 1, b + 3), (b + 7, b + 11)]),
            item(group * 4 + 3, &[(b + 11, b + 15)]),
        ]);
    }
    inputs.reverse();
    let result = pack(&inputs, 50_000, &options(2)).unwrap();
    assert_eq!(result.assignments.len(), 10_000);
    assert_eq!(result.lane_count, 2);
    assert_eq!(result.overflow_ids().count(), 2500);
    for (assignment, original) in result.assignments.iter().zip(&inputs) {
        assert_eq!(assignment.id, original.id);
        assert_eq!(assignment.spans, original.spans);
        assert_eq!(
            assignment.placement,
            match assignment.id % 4 {
                1 => Placement::Overflow,
                2 => Placement::Lane(1),
                _ => Placement::Lane(0),
            }
        );
    }
}

#[test]
fn rejected_joined_item_does_not_partially_occupy_lane() {
    let inputs = [
        item(1, &[(0, 2), (8, 10)]),
        item(2, &[(2, 4), (9, 11)]),
        item(3, &[(3, 4)]),
    ];
    let result = pack(&inputs, 11, &options(1)).unwrap();
    assert_eq!(
        result
            .assignments
            .iter()
            .map(|a| a.placement)
            .collect::<Vec<_>>(),
        vec![Placement::Lane(0), Placement::Overflow, Placement::Lane(0)]
    );
}

// This oracle uses a four-cell bitset, not the production span-conflict code.
// Only tests have a per-position representation.
fn runs(mask: u8) -> Vec<Interval> {
    let mut result = vec![];
    let mut at = 0;
    while at < 4 {
        if mask & (1 << at) == 0 {
            at += 1;
            continue;
        }
        let start = at;
        while at < 4 && mask & (1 << at) != 0 {
            at += 1;
        }
        result.push(Interval { start, end: at });
    }
    result
}

fn mask(spans: &[Interval]) -> u8 {
    let mut result = 0;
    for s in spans {
        for i in s.start..s.end {
            result |= 1 << i;
        }
    }
    result
}

#[test]
fn exhaustive_small_domain_first_fit_collision_coverage_and_permutation_oracle() {
    for a in 1u8..16 {
        for b in 1u8..16 {
            for c in 1u8..16 {
                let specs = [(19, a), (2, b), (80, c)];
                // Deliberately represent adjacent covered cells as separate, reversed
                // spans; normalization must independently recover each union.
                let input: Vec<_> = specs
                    .iter()
                    .map(|&(id, bits)| Item {
                        id,
                        spans: (0..4)
                            .rev()
                            .filter(|&i| bits & (1 << i) != 0)
                            .map(|i| Interval {
                                start: i,
                                end: i + 1,
                            })
                            .collect(),
                    })
                    .collect();
                let mut ordered = specs;
                ordered.sort_by_key(|&(id, bits)| (runs(bits), id));
                for cap in 0..4 {
                    let mut expected = BTreeMap::new();
                    let mut occupied = vec![0u8; cap];
                    for (id, bits) in ordered {
                        let lane = occupied.iter().position(|&used| used & bits == 0);
                        if let Some(lane) = lane {
                            occupied[lane] |= bits;
                        }
                        expected
                            .insert(id, lane.map(Placement::Lane).unwrap_or(Placement::Overflow));
                    }
                    for candidate in [input.clone(), input.iter().rev().cloned().collect()] {
                        let actual = pack(&candidate, 4, &options(cap)).unwrap();
                        assert_eq!(actual.assignments.len(), 3);
                        assert_eq!(
                            actual.lane_count,
                            occupied.iter().filter(|&&bits| bits != 0).count()
                        );
                        let mut checked = vec![0u8; cap];
                        for (assignment, original) in actual.assignments.iter().zip(&candidate) {
                            assert_eq!(assignment.id, original.id);
                            let bits = mask(&original.spans);
                            assert_eq!(assignment.spans, runs(bits));
                            assert_eq!(assignment.placement, expected[&original.id]);
                            if let Placement::Lane(lane) = assignment.placement {
                                assert_eq!(checked[lane] & bits, 0);
                                checked[lane] |= bits;
                            }
                        }
                    }
                }
            }
        }
    }
}
