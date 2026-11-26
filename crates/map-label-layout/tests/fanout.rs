use collagene_map_label_layout::{
    Error, Label, Options, OverflowReason, Point, Rect, Side, Size, Viewport,
    radial::{self, FanLayout, FanOutcome},
};
use std::{
    collections::BTreeSet,
    f64::consts::{PI, TAU},
};
fn view() -> (Viewport, Point, f64) {
    (
        Viewport {
            width: 760.,
            height: 600.,
            padding: 12.,
        },
        Point { x: 380., y: 300. },
        150.,
    )
}
fn labels(count: usize, start: f64, sweep: f64, center: Point, radius: f64) -> Vec<Label> {
    (0..count)
        .map(|i| {
            let angle = start + sweep * i as f64 / count as f64;
            Label {
                id: i as u64,
                stable_rank: i as u64,
                anchor: Point {
                    x: center.x + radius * angle.cos(),
                    y: center.y + radius * angle.sin(),
                },
                // The circular API uses actual anchors, rather than a caller's
                // unrelated or stale side/desired-column coordinates.
                preferred_center: Point { x: 0., y: 0. },
                side: Side::Above,
                measured: Size {
                    width: 72. + (i % 3) as f64 * 4.,
                    height: 12.,
                },
                priority: 1,
            }
        })
        .collect()
}
fn rect_distance(c: Point, r: Rect) -> f64 {
    let dx = if c.x < r.x {
        r.x - c.x
    } else if c.x > r.right() {
        c.x - r.right()
    } else {
        0.
    };
    let dy = if c.y < r.y {
        r.y - c.y
    } else if c.y > r.bottom() {
        c.y - r.bottom()
    } else {
        0.
    };
    dx.hypot(dy)
}
fn segment_distance(c: Point, a: Point, b: Point) -> f64 {
    let dx = b.x - a.x;
    let dy = b.y - a.y;
    let dot = (a.x - c.x) * dx + (a.y - c.y) * dy;
    let square = dx * dx + dy * dy;
    if dot >= 0. || square == 0. {
        return (a.x - c.x).hypot(a.y - c.y);
    }
    if dot + square <= 0. {
        return (b.x - c.x).hypot(b.y - c.y);
    }
    ((a.x - c.x) * dy - (a.y - c.y) * dx).abs() / square.sqrt()
}
fn verify(input: &[Label], layout: &FanLayout, viewport: Viewport, center: Point, radius: f64) {
    assert_eq!(layout.labels.len(), input.len());
    assert_eq!(
        layout.displayed,
        layout
            .labels
            .iter()
            .filter(|r| matches!(r.outcome, FanOutcome::Placed { .. }))
            .count()
    );
    let mut placed = Vec::new();
    for (source, result) in input.iter().zip(&layout.labels) {
        assert_eq!(result.id, source.id);
        if let FanOutcome::Placed { rect, leader } = result.outcome {
            assert_eq!(leader.anchor, source.anchor);
            assert_eq!(rect.width, source.measured.width);
            assert_eq!(rect.height, source.measured.height.max(14.));
            assert!(rect.x >= viewport.padding && rect.y >= viewport.padding);
            assert!(
                rect.right() <= viewport.width - viewport.padding
                    && rect.bottom() <= viewport.height - viewport.padding
            );
            assert!(rect_distance(center, rect) >= radius - 1e-7);
            let dx = leader.knee.x - leader.anchor.x;
            let dy = leader.knee.y - leader.anchor.y;
            let ax = leader.anchor.x - center.x;
            let ay = leader.anchor.y - center.y;
            assert!(
                (dx.hypot(dy) - 8.).abs() < 1e-7,
                "actual radial stem length"
            );
            assert!(
                (ax * dy - ay * dx).abs() < 1e-6,
                "stem must be radial, not a horizontal-first L"
            );
            assert!(ax * dx + ay * dy > 0., "stem points outwards");
            assert!(segment_distance(center, leader.anchor, leader.knee) >= radius - 1e-7);
            assert!(segment_distance(center, leader.knee, leader.endpoint) >= radius - 1e-7);
            assert_eq!(leader.endpoint.y, rect.y + rect.height / 2.);
            if leader.anchor.x >= center.x {
                assert_eq!(leader.endpoint.x, rect.x);
            } else {
                assert_eq!(leader.endpoint.x, rect.right());
            }
            placed.push((source, rect));
        }
    }
    for (i, (a, first)) in placed.iter().enumerate() {
        for (b, second) in &placed[i + 1..] {
            assert!(
                first.right() <= second.x
                    || second.right() <= first.x
                    || first.bottom() <= second.y
                    || second.bottom() <= first.y,
                "measured labels collide"
            );
            if (a.anchor.x >= center.x) == (b.anchor.x >= center.x) {
                let angular_a =
                    (a.anchor.y - center.y) / (a.anchor.x - center.x).hypot(a.anchor.y - center.y);
                let angular_b =
                    (b.anchor.y - center.y) / (b.anchor.x - center.x).hypot(b.anchor.y - center.y);
                if angular_a < angular_b {
                    assert!(first.bottom() + 2. <= second.y + 1e-7, "angular row order");
                }
                if angular_a > angular_b {
                    assert!(second.bottom() + 2. <= first.y + 1e-7, "angular row order");
                }
            }
        }
    }
}
#[test]
fn twelve_sparse_compass_sites_keep_radial_stems_and_local_contours() {
    let (viewport, center, radius) = view();
    let input = labels(12, -PI / 2., TAU, center, radius);
    let layout = radial::fan_out(
        &input,
        viewport,
        center,
        radius,
        &Options {
            gap: 0.,
            ..Default::default()
        },
    )
    .unwrap();
    verify(&input, &layout, viewport, center, radius);
    assert_eq!(layout.displayed, 12);
    assert_eq!(layout.considered, 12);
    let mut edges = Vec::new();
    for result in &layout.labels {
        let FanOutcome::Placed { leader, .. } = result.outcome else {
            unreachable!()
        };
        edges.push((leader.endpoint.x - center.x).abs());
        // Sparse sites remain close to the plasmid, rather than all acquiring a
        // long route to one common outside column.
        assert!((leader.endpoint.x - leader.knee.x).hypot(leader.endpoint.y - leader.knee.y) < 50.);
    }
    edges.sort_unstable_by(f64::total_cmp);
    edges.dedup_by(|a, b| (*a - *b).abs() < 1.);
    assert!(edges.len() >= 4, "contour-following horizontal positions");
}
#[test]
fn crowded_quadrant_fans_both_ways_and_reports_every_omitted_site() {
    let (viewport, center, radius) = view();
    let input = labels(40, -1.35, 1., center, radius);
    let layout = radial::fan_out(
        &input,
        viewport,
        center,
        radius,
        &Options {
            gap: 2.,
            ..Default::default()
        },
    )
    .unwrap();
    verify(&input, &layout, viewport, center, radius);
    assert_eq!(layout.considered, 40);
    assert!(
        layout.displayed >= 16,
        "a usable fan must retain a substantial subset: {}",
        layout.displayed
    );
    assert!(
        layout.displayed < 40,
        "forty14px rows cannot fit in this side's bounded fan"
    );
    let overflows = layout
        .labels
        .iter()
        .filter(|r| matches!(r.outcome, FanOutcome::Overflow { .. }))
        .count();
    assert_eq!(layout.displayed + overflows, 40);
    let rows: Vec<_> = layout
        .labels
        .iter()
        .filter_map(|r| match r.outcome {
            FanOutcome::Placed { rect, .. } => Some(rect.y + rect.height / 2.),
            _ => None,
        })
        .collect();
    assert!(rows.iter().copied().fold(f64::INFINITY, f64::min) < input[0].anchor.y);
    assert!(
        rows.iter().copied().fold(f64::NEG_INFINITY, f64::max) > input.last().unwrap().anchor.y
    );
}
#[test]
fn rotation_rebuilds_actual_anchor_geometry_without_reordering_ids() {
    let (viewport, center, radius) = view();
    for angle in [0., 0.17, 0.63, PI / 2., PI, TAU - 0.03] {
        let input = labels(12, angle, TAU, center, radius);
        let layout = radial::fan_out(
            &input,
            viewport,
            center,
            radius,
            &Options {
                gap: 2.,
                ..Default::default()
            },
        )
        .unwrap();
        verify(&input, &layout, viewport, center, radius);
        assert_eq!(layout.displayed, 12, "sparse rotation {angle}");
        assert_eq!(
            layout,
            radial::fan_out(
                &input,
                viewport,
                center,
                radius,
                &Options {
                    gap: 2.,
                    ..Default::default()
                }
            )
            .unwrap()
        );
    }
}
#[test]
fn small_viewport_and_measured_widths_produce_explicit_bounded_overflow() {
    let viewport = Viewport {
        width: 96.,
        height: 92.,
        padding: 4.,
    };
    let center = Point { x: 48., y: 46. };
    let input = labels(12, 0., TAU, center, 32.);
    let layout = radial::fan_out(&input, viewport, center, 32., &Options::default()).unwrap();
    assert_eq!(layout.considered, 12);
    assert_eq!(layout.displayed, 0);
    for result in &layout.labels {
        let FanOutcome::Overflow {
            reason,
            available: Some(size),
        } = result.outcome
        else {
            panic!("explicit measured overflow");
        };
        assert_eq!(reason, OverflowReason::TooWide);
        assert!(size.width >= 0. && size.width < 72. && size.height == 84.);
    }
    let mut tall = labels(1, 0., 0., center, 32.);
    tall[0].measured.height = 85.;
    let layout = radial::fan_out(&tall, viewport, center, 32., &Options::default()).unwrap();
    assert!(matches!(
        layout.labels[0].outcome,
        FanOutcome::Overflow {
            reason: OverflowReason::TooTall,
            ..
        }
    ));
}
#[test]
fn priority_budget_is_stable_across_input_order_and_camera_rotation() {
    let (viewport, center, radius) = view();
    let options = Options {
        max_displayed: 3,
        gap: 2.,
        ..Default::default()
    };
    let mut input = labels(40, 0., TAU, center, radius);
    input[7].priority = 100;
    let layout = radial::fan_out(&input, viewport, center, radius, &options).unwrap();
    verify(&input, &layout, viewport, center, radius);
    assert_eq!(layout.considered, 3);
    assert_eq!(layout.displayed, 3);
    assert!(matches!(
        layout.labels[7].outcome,
        FanOutcome::Placed { .. }
    ));
    let selected = |layout: &FanLayout| {
        layout
            .labels
            .iter()
            .filter_map(|r| {
                (!matches!(
                    r.outcome,
                    FanOutcome::Overflow {
                        reason: OverflowReason::DisplayBudget,
                        ..
                    }
                ))
                .then_some(r.id)
            })
            .collect::<BTreeSet<_>>()
    };
    let selected_before = selected(&layout);
    input.reverse();
    let reversed = radial::fan_out(&input, viewport, center, radius, &options).unwrap();
    assert_eq!(selected(&reversed), selected_before);
    for result in &layout.labels {
        assert_eq!(
            reversed
                .labels
                .iter()
                .find(|r| r.id == result.id)
                .unwrap()
                .outcome,
            result.outcome
        );
    }
    let mut rotated = labels(40, 0.4, TAU, center, radius);
    rotated[7].priority = 100;
    let rotated_layout = radial::fan_out(&rotated, viewport, center, radius, &options).unwrap();
    verify(&rotated, &rotated_layout, viewport, center, radius);
    assert_eq!(selected(&rotated_layout), selected_before);
}
#[test]
fn invalid_anchors_ids_and_resource_limits_fail_before_placement() {
    let (viewport, center, radius) = view();
    let mut input = labels(2, 0., TAU, center, radius);
    input[0].anchor = center;
    assert_eq!(
        radial::fan_out(&input, viewport, center, radius, &Options::default()),
        Err(Error::InvalidLabel(0))
    );
    let mut input = labels(2, 0., TAU, center, radius);
    input[1].id = 0;
    assert_eq!(
        radial::fan_out(&input, viewport, center, radius, &Options::default()),
        Err(Error::DuplicateId(0))
    );
    assert_eq!(
        radial::fan_out(
            &input,
            viewport,
            center,
            radius,
            &Options {
                max_inputs: 1,
                ..Default::default()
            }
        ),
        Err(Error::InputLimit)
    );
    let input = labels(2, 0., TAU, center, radius);
    let zero = radial::fan_out(
        &input,
        viewport,
        center,
        radius,
        &Options {
            max_displayed: 0,
            ..Default::default()
        },
    )
    .unwrap();
    assert_eq!(zero.considered, 0);
    assert_eq!(zero.displayed, 0);
    assert!(zero.labels.iter().all(|r| matches!(
        r.outcome,
        FanOutcome::Overflow {
            reason: OverflowReason::DisplayBudget,
            ..
        }
    )));
}
