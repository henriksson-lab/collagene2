use collagene_map_label_layout::*;
use std::collections::BTreeMap;

fn viewport() -> Viewport {
    Viewport {
        width: 300.,
        height: 200.,
        padding: 10.,
    }
}
fn column(side: Side) -> Column {
    Column {
        side,
        rect: match side {
            Side::Left => Rect {
                x: 10.,
                y: 30.,
                width: 70.,
                height: 140.,
            },
            Side::Right => Rect {
                x: 220.,
                y: 30.,
                width: 70.,
                height: 140.,
            },
            Side::Above => Rect {
                x: 90.,
                y: 10.,
                width: 120.,
                height: 40.,
            },
            Side::Below => Rect {
                x: 90.,
                y: 150.,
                width: 120.,
                height: 40.,
            },
        },
    }
}
fn label(id: u64, side: Side) -> Label {
    Label {
        id,
        stable_rank: id,
        anchor: Point { x: 150., y: 100. },
        preferred_center: Point { x: 150., y: 100. },
        side,
        measured: Size {
            width: 40.,
            height: 14.,
        },
        priority: 0,
    }
}
fn placed(result: &Layout) -> Vec<(u64, Rect, Leader)> {
    result
        .labels
        .iter()
        .filter_map(|r| match r.outcome {
            Outcome::Placed { rect, leader } => Some((r.id, rect, leader)),
            _ => None,
        })
        .collect()
}
fn independent_check(inputs: &[Label], columns: &[Column], result: &Layout, gap: f64) {
    assert_eq!(result.labels.len(), inputs.len());
    let mut boxes = vec![];
    for (input, result) in inputs.iter().zip(&result.labels) {
        assert_eq!(input.id, result.id);
        if let Outcome::Placed { rect, leader } = result.outcome {
            let region = columns.iter().find(|c| c.side == input.side).unwrap().rect;
            assert!(rect.x >= region.x && rect.y >= region.y);
            assert!(rect.right() <= region.right() && rect.bottom() <= region.bottom());
            assert_eq!(rect.size(), input.measured);
            assert_eq!(leader.anchor, input.anchor);
            assert!(leader.endpoint.x >= rect.x && leader.endpoint.x <= rect.right());
            assert!(leader.endpoint.y >= rect.y && leader.endpoint.y <= rect.bottom());
            for &(side, other) in &boxes {
                let other: Rect = other;
                assert!(
                    rect.right() <= other.x
                        || other.right() <= rect.x
                        || rect.bottom() <= other.y
                        || other.bottom() <= rect.y
                );
                if side == input.side {
                    if matches!(side, Side::Left | Side::Right) {
                        assert!(rect.bottom() + gap <= other.y || other.bottom() + gap <= rect.y);
                    } else {
                        assert!(rect.right() + gap <= other.x || other.right() + gap <= rect.x);
                    }
                }
            }
            boxes.push((input.side, rect));
        }
    }
    assert_eq!(result.displayed, boxes.len());
}

#[test]
fn six_colliding_linear_labels_fit_without_text_collision() {
    let labels: Vec<_> = (0..6)
        .map(|id| Label {
            measured: Size {
                width: 16.,
                height: 14.,
            },
            ..label(id, Side::Above)
        })
        .collect();
    let columns = [column(Side::Above)];
    let result = layout(
        &labels,
        viewport(),
        &columns,
        &Options {
            gap: 4.,
            ..Options::default()
        },
    )
    .unwrap();
    assert_eq!(result.displayed, 6);
    assert_eq!(placed(&result)[0].1.x, 94.);
    assert_eq!(placed(&result)[5].1.right(), 210.);
    independent_check(&labels, &columns, &result, 4.);
}

#[test]
fn circular_sides_boundary_preferences_and_original_anchors() {
    let columns = [column(Side::Left), column(Side::Right)];
    let labels: Vec<_> = (0..12)
        .map(|id| Label {
            anchor: Point {
                x: if id % 2 == 0 { -500. } else { 500. },
                y: 100. + id as f64,
            },
            preferred_center: Point {
                x: if id % 2 == 0 { -100. } else { 900. },
                y: if id % 3 == 0 { -500. } else { 900. },
            },
            ..label(id, if id % 2 == 0 { Side::Left } else { Side::Right })
        })
        .collect();
    let result = layout(&labels, viewport(), &columns, &Options::default()).unwrap();
    assert_eq!(result.displayed, 12);
    independent_check(&labels, &columns, &result, 4.);
    for (id, rect, leader) in placed(&result) {
        assert_eq!(
            leader.endpoint.x,
            if id % 2 == 0 { rect.right() } else { rect.x }
        );
        assert_eq!(
            leader.anchor.x.to_bits(),
            labels[id as usize].anchor.x.to_bits()
        );
    }
}

#[test]
fn explicit_size_density_missing_and_zero_budget_overflow() {
    let inputs = [
        label(0, Side::Left),
        Label {
            measured: Size {
                width: 71.,
                height: 14.,
            },
            ..label(1, Side::Left)
        },
        Label {
            measured: Size {
                width: 20.,
                height: 141.,
            },
            ..label(2, Side::Left)
        },
        label(3, Side::Below),
    ];
    let result = layout(
        &inputs,
        viewport(),
        &[column(Side::Left)],
        &Options::default(),
    )
    .unwrap();
    assert!(matches!(
        result.labels[1].outcome,
        Outcome::Overflow {
            reason: OverflowReason::TooWide,
            available: Some(Size {
                width: 70.,
                height: 140.
            })
        }
    ));
    assert!(matches!(
        result.labels[2].outcome,
        Outcome::Overflow {
            reason: OverflowReason::TooTall,
            ..
        }
    ));
    assert!(matches!(
        result.labels[3].outcome,
        Outcome::Overflow {
            reason: OverflowReason::MissingColumn,
            ..
        }
    ));
    let dense: Vec<_> = (0..20).map(|id| label(id, Side::Left)).collect();
    let result = layout(
        &dense,
        viewport(),
        &[column(Side::Left)],
        &Options::default(),
    )
    .unwrap();
    assert_eq!(result.displayed, 8);
    assert_eq!(
        result
            .labels
            .iter()
            .filter(|r| matches!(
                r.outcome,
                Outcome::Overflow {
                    reason: OverflowReason::Dense,
                    ..
                }
            ))
            .count(),
        12
    );
    let zero = layout(
        &dense,
        viewport(),
        &[column(Side::Left)],
        &Options {
            max_displayed: 0,
            ..Options::default()
        },
    )
    .unwrap();
    assert_eq!(zero.considered, 0);
    assert_eq!(zero.displayed, 0);
    assert!(zero.labels.iter().all(|r| matches!(
        r.outcome,
        Outcome::Overflow {
            reason: OverflowReason::DisplayBudget,
            ..
        }
    )));
}

#[test]
fn literal_rank_stratification_and_explicit_priority() {
    let mut labels: Vec<_> = (0..10).map(|id| label(id, Side::Left)).collect();
    let options = Options {
        max_displayed: 3,
        ..Options::default()
    };
    let selected = |result: Layout| {
        result
            .labels
            .into_iter()
            .filter(|r| {
                !matches!(
                    r.outcome,
                    Outcome::Overflow {
                        reason: OverflowReason::DisplayBudget,
                        ..
                    }
                )
            })
            .map(|r| r.id)
            .collect::<Vec<_>>()
    };
    assert_eq!(
        selected(layout(&labels, viewport(), &[column(Side::Left)], &options).unwrap()),
        vec![1, 5, 8]
    );
    labels[0].priority = 1;
    assert_eq!(
        selected(layout(&labels, viewport(), &[column(Side::Left)], &options).unwrap()),
        vec![0, 3, 7]
    );
}

#[test]
fn budget_ids_are_invariant_under_camera_projection_side_and_input_permutations() {
    let original: Vec<_> = (0..100)
        .map(|id| Label {
            stable_rank: 1000 - id,
            ..label(id, Side::Left)
        })
        .collect();
    let options = Options {
        max_displayed: 12,
        ..Options::default()
    };
    let columns = [column(Side::Left), column(Side::Right)];
    let considered = |result: Layout| {
        result
            .labels
            .into_iter()
            .filter(|r| {
                !matches!(
                    r.outcome,
                    Outcome::Overflow {
                        reason: OverflowReason::DisplayBudget,
                        ..
                    }
                )
            })
            .map(|r| (r.id, true))
            .collect::<BTreeMap<_, _>>()
    };
    let expected = considered(layout(&original, viewport(), &columns, &options).unwrap());
    for step in 0..24 {
        let mut rotated = original.clone();
        for label in &mut rotated {
            let angle = (label.id as f64 / 100. + step as f64 / 24.) * std::f64::consts::TAU;
            label.anchor = Point {
                x: 150. + 80. * angle.cos(),
                y: 100. + 80. * angle.sin(),
            };
            label.preferred_center = label.anchor;
            label.side = if angle.cos() < 0. {
                Side::Left
            } else {
                Side::Right
            };
        }
        if step % 2 == 0 {
            rotated.reverse();
        }
        let result = layout(&rotated, viewport(), &columns, &options).unwrap();
        independent_check(&rotated, &columns, &result, 4.);
        assert_eq!(considered(result), expected);
    }
}

#[test]
fn finite_validation_duplicate_ids_and_resource_policy() {
    let columns = [column(Side::Left)];
    for value in [
        f64::NAN,
        f64::INFINITY,
        f64::NEG_INFINITY,
        MAX_COORDINATE + 1.,
    ] {
        let bad = Label {
            anchor: Point { x: value, y: 0. },
            ..label(0, Side::Left)
        };
        assert!(matches!(
            layout(&[bad], viewport(), &columns, &Options::default()),
            Err(Error::InvalidLabel(0))
        ));
    }
    for width in [0., -1., f64::NAN, f64::INFINITY, MAX_VIEWPORT + 1.] {
        assert!(
            layout(
                &[Label {
                    measured: Size { width, height: 14. },
                    ..label(0, Side::Left)
                }],
                viewport(),
                &columns,
                &Options::default()
            )
            .is_err()
        );
    }
    assert_eq!(
        layout(
            &[label(1, Side::Left), label(1, Side::Right)],
            viewport(),
            &columns,
            &Options::default()
        ),
        Err(Error::DuplicateId(1))
    );
    assert_eq!(
        layout(
            &[label(1, Side::Left)],
            viewport(),
            &columns,
            &Options {
                max_inputs: 0,
                ..Options::default()
            }
        ),
        Err(Error::InputLimit)
    );
    assert_eq!(
        layout(
            &[],
            viewport(),
            &columns,
            &Options {
                max_displayed: MAX_DISPLAYED + 1,
                ..Options::default()
            }
        ),
        Err(Error::InvalidOptions)
    );
    assert!(
        layout(
            &[],
            Viewport {
                padding: 100.,
                ..viewport()
            },
            &columns,
            &Options::default()
        )
        .is_err()
    );
    assert_eq!(
        layout(
            &[],
            viewport(),
            &[columns[0], columns[0]],
            &Options::default()
        ),
        Err(Error::DuplicateSide)
    );
    assert_eq!(
        layout(
            &[],
            viewport(),
            &[
                columns[0],
                Column {
                    side: Side::Right,
                    ..columns[0]
                }
            ],
            &Options::default()
        ),
        Err(Error::OverlappingColumns)
    );
    assert!(
        layout(
            &[],
            viewport(),
            &[Column {
                rect: Rect {
                    x: 0.,
                    ..columns[0].rect
                },
                ..columns[0]
            }],
            &Options::default()
        )
        .is_err()
    );
}

#[test]
fn small_exhaustive_preferences_sizes_and_all_four_sides() {
    let columns = [
        column(Side::Left),
        column(Side::Right),
        column(Side::Above),
        column(Side::Below),
    ];
    for side in [Side::Left, Side::Right, Side::Above, Side::Below] {
        for a in [-100., 10., 100., 250.] {
            for b in [-100., 10., 100., 250.] {
                for c in [-100., 10., 100., 250.] {
                    for size in [10., 23., 60.] {
                        let labels: Vec<_> = [a, b, c]
                            .into_iter()
                            .enumerate()
                            .map(|(i, p)| Label {
                                preferred_center: Point { x: p, y: p },
                                measured: Size {
                                    width: size,
                                    height: size,
                                },
                                ..label(i as u64, side)
                            })
                            .collect();
                        let options = Options {
                            gap: 3.,
                            ..Options::default()
                        };
                        let result = layout(&labels, viewport(), &columns, &options).unwrap();
                        independent_check(&labels, &columns, &result, 3.);
                        let mut reversed = labels.clone();
                        reversed.reverse();
                        let by_id = |result: Layout| {
                            result
                                .labels
                                .into_iter()
                                .map(|r| (r.id, r.outcome))
                                .collect::<BTreeMap<_, _>>()
                        };
                        assert_eq!(
                            by_id(layout(&reversed, viewport(), &columns, &options).unwrap()),
                            by_id(result)
                        );
                    }
                }
            }
        }
    }
}

#[test]
fn ten_thousand_dense_and_sparse_inputs_are_bounded_and_loss_explicit() {
    let columns = [column(Side::Left), column(Side::Right)];
    for dense in [false, true] {
        let labels: Vec<_> = (0..10_000)
            .rev()
            .map(|id| Label {
                preferred_center: Point {
                    x: 30.,
                    y: if dense { 100. } else { id as f64 / 50. },
                },
                side: if id % 2 == 0 { Side::Left } else { Side::Right },
                ..label(id, Side::Left)
            })
            .collect();
        let result = layout(&labels, viewport(), &columns, &Options::default()).unwrap();
        assert_eq!(result.labels.len(), 10_000);
        assert_eq!(result.considered, 128);
        assert_eq!(
            result
                .labels
                .iter()
                .filter(|r| matches!(
                    r.outcome,
                    Outcome::Overflow {
                        reason: OverflowReason::DisplayBudget,
                        ..
                    }
                ))
                .count(),
            9872
        );
        assert_eq!(result.displayed, 16);
        independent_check(&labels, &columns, &result, 4.);
    }
}

#[test]
fn fractional_exact_fit_and_tiny_boxes_do_not_invert_clamp_bounds() {
    for side in [Side::Left, Side::Right, Side::Above, Side::Below] {
        for (origin, extent) in [
            (0.1, 0.7),
            (0.7, 0.1),
            (10.1, 0.7),
            (10., f64::MIN_POSITIVE),
        ] {
            let region = Column {
                side,
                rect: Rect {
                    x: origin,
                    y: origin,
                    width: extent,
                    height: extent,
                },
            };
            let input = Label {
                measured: Size {
                    width: extent,
                    height: extent,
                },
                ..label(1, side)
            };
            let result = layout(
                &[input],
                Viewport {
                    width: 100.,
                    height: 100.,
                    padding: 0.,
                },
                &[region],
                &Options {
                    gap: 0.,
                    ..Options::default()
                },
            )
            .unwrap();
            if extent == f64::MIN_POSITIVE {
                assert_eq!(result.displayed, 0);
                assert!(matches!(
                    result.labels[0].outcome,
                    Outcome::Overflow {
                        reason: OverflowReason::Precision,
                        ..
                    }
                ));
            } else {
                assert_eq!(result.displayed, 1);
            }
            independent_check(&[input], &[region], &result, 0.);
        }
    }
}
