use collagene_map_label_layout::{radial, *};
use std::f64::consts::TAU;
fn label(id: u64, angle: f64) -> Label {
    Label {
        id,
        stable_rank: id,
        anchor: Point {
            x: 400. + 200. * angle.cos(),
            y: 350. + 200. * angle.sin(),
        },
        preferred_center: Point { x: 0., y: 0. },
        side: Side::Above,
        measured: Size {
            width: 76.,
            height: 18.,
        },
        priority: 0,
    }
}
fn view() -> Viewport {
    Viewport {
        width: 800.,
        height: 700.,
        padding: 10.,
    }
}
fn options() -> Options {
    Options {
        max_inputs: 128,
        max_displayed: 64,
        gap: 4.,
    }
}
fn distance_to_segment(p: Point, a: Point, b: Point) -> f64 {
    let (dx, dy) = (b.x - a.x, b.y - a.y);
    let n = dx * dx + dy * dy;
    let t = if n > 0. {
        ((p.x - a.x) * dx + (p.y - a.y) * dy) / n
    } else {
        0.
    }
    .clamp(0., 1.);
    (p.x - a.x - t * dx).hypot(p.y - a.y - t * dy)
}
#[test]
fn circumferential_labels_keep_true_anchors_and_both_elbow_segments_outside_backbone() {
    let labels = (0..16)
        .map(|i| label(i, (i as f64 + 0.15) * TAU / 16.))
        .collect::<Vec<_>>();
    let out = radial::outside_circle(
        &labels,
        view(),
        Point { x: 400., y: 350. },
        200.,
        &options(),
    )
    .unwrap();
    assert_eq!(out.displayed, labels.len());
    for (i, r) in out.labels.iter().enumerate() {
        let Outcome::Placed { rect, leader } = r.outcome else {
            panic!("sparse labels should fit")
        };
        assert_eq!(leader.anchor, labels[i].anchor);
        let elbow = Point {
            x: leader.endpoint.x,
            y: leader.anchor.y,
        };
        assert!(
            distance_to_segment(Point { x: 400., y: 350. }, leader.anchor, elbow) >= 200. - 1e-8
        );
        assert!(
            distance_to_segment(Point { x: 400., y: 350. }, elbow, leader.endpoint) >= 200. - 1e-8
        );
        assert!(rect.x >= 10. && rect.y >= 10. && rect.right() <= 790. && rect.bottom() <= 690.);
        let near = Point {
            x: 400f64.clamp(rect.x, rect.right()),
            y: 350f64.clamp(rect.y, rect.bottom()),
        };
        assert!((near.x - 400.).hypot(near.y - 350.) >= 200.);
        assert_eq!(leader.endpoint.y, rect.y + rect.height / 2.);
    }
    for side in [false, true] {
        let mut ordered = labels
            .iter()
            .enumerate()
            .filter(|(_, l)| (l.anchor.x >= 400.) == side)
            .collect::<Vec<_>>();
        ordered.sort_by(|a, b| a.1.anchor.y.total_cmp(&b.1.anchor.y));
        for pair in ordered.windows(2) {
            let a = out.labels[pair[0].0].outcome;
            let b = out.labels[pair[1].0].outcome;
            if let (Outcome::Placed { rect: a, .. }, Outcome::Placed { rect: b, .. }) = (a, b) {
                assert!(a.bottom() <= b.y);
            }
        }
    }
}
#[test]
fn dense_side_retains_collision_free_hit_rectangles_without_opposite_edge_bands() {
    let labels = (0..20)
        .map(|i| label(i, 0.0002 * i as f64))
        .collect::<Vec<_>>();
    let out = radial::outside_circle(
        &labels,
        view(),
        Point { x: 400., y: 350. },
        200.,
        &options(),
    )
    .unwrap();
    assert_eq!(out.displayed, 20);
    let mut boxes = vec![];
    for r in out.labels {
        let Outcome::Placed { rect, leader } = r.outcome else {
            panic!()
        };
        assert!(rect.x > 600.);
        assert_eq!(leader.anchor.x >= 400., leader.endpoint.x >= 400.);
        boxes.push(rect);
    }
    boxes.sort_by(|a, b| a.y.total_cmp(&b.y));
    for p in boxes.windows(2) {
        assert!(p[0].bottom() + 4. <= p[1].y + 1e-8);
    }
}
#[test]
fn short_and_reverse_arrow_heads_remain_inside_the_literal_feature_and_inner_lane() {
    for (start, end) in [(0., 1.), (2., 2.000001), (-1.5, -1.49)] {
        for reverse in [false, true] {
            let (tip, base) = radial::arrow_head(start, end, reverse, 0.055).unwrap();
            assert_eq!(tip, if reverse { start } else { end });
            assert!(base > start && base < end);
        }
    }
    assert_eq!(radial::arrow_head(1., 1., false, 0.05), None);
    assert_eq!(radial::arrow_head(0., 1., false, f64::NAN), None);
    for lane in 0..3 {
        let (inner, outer) = radial::feature_band(lane);
        assert!(inner > 0. && outer < 0.65 && outer - inner >= 0.0399);
        if lane > 0 {
            assert!(outer < radial::feature_band(lane - 1).0);
        }
    }
}
#[test]
fn radial_shape_guards_and_explicit_overflow_preserve_all_ids() {
    let mut l = label(3, 0.);
    l.anchor = Point { x: 400., y: 350. };
    assert!(
        radial::outside_circle(&[l], view(), Point { x: 400., y: 350. }, 200., &options()).is_err()
    );
    assert!(
        radial::outside_circle(
            &[label(3, 0.)],
            view(),
            Point { x: 400., y: 350. },
            f64::NAN,
            &options()
        )
        .is_err()
    );
    let mut l = label(3, 0.);
    l.measured.width = 350.;
    let out =
        radial::outside_circle(&[l], view(), Point { x: 400., y: 350. }, 200., &options()).unwrap();
    assert_eq!(out.labels[0].id, 3);
    assert!(matches!(
        out.labels[0].outcome,
        Outcome::Overflow {
            reason: OverflowReason::TooWide,
            ..
        }
    ));
}
