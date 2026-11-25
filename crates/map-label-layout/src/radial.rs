//! Circumferential text placement with two-segment leaders that cannot cross the disk.
use crate::{
    Column, Error, Label, Layout, Options, Outcome, OverflowReason, Point, Rect, Side, Size,
    Viewport, layout,
};
/// Place horizontal labels near their actual angular position. Left/right vertical
/// order follows their anchors. The final elbow is at (endpoint.x, anchor.y), so
/// both leader segments remain outside the excluded disk; displaced labels never
/// connect to an arbitrary opposite-edge band.
pub fn outside_circle(
    labels: &[Label],
    viewport: Viewport,
    center: Point,
    radius: f64,
    options: &Options,
) -> Result<Layout, Error> {
    crate::preflight(labels, viewport, &[], options)?;
    if !radius.is_finite()
        || radius <= 0.
        || radius > crate::MAX_VIEWPORT
        || !center.x.is_finite()
        || !center.y.is_finite()
    {
        return Err(Error::InvalidOptions);
    }
    for l in labels {
        if (l.anchor.x - center.x).hypot(l.anchor.y - center.y) + 1e-6 < radius {
            return Err(Error::InvalidLabel(l.id));
        }
    }
    // Existing layout owns finite coordinates, input bounds, priority sampling,
    // deterministic collision spacing and explicit overflow.
    let mut input = labels.to_vec();
    for l in &mut input {
        l.side = if l.anchor.x < center.x {
            Side::Left
        } else {
            Side::Right
        };
        let dy = l.anchor.y - center.y;
        let dx = l.anchor.x - center.x;
        let scale = if radius > 0. {
            (radius + 18.) / radius
        } else {
            1.
        };
        l.preferred_center = Point {
            x: center.x + dx * scale,
            y: center.y + dy * scale,
        };
    }
    let pad = viewport.padding;
    let width = (viewport.width / 2. - pad - 8.).max(0.);
    let columns = [
        Column {
            side: Side::Left,
            rect: Rect {
                x: pad,
                y: pad,
                width,
                height: (viewport.height - 2. * pad).max(0.),
            },
        },
        Column {
            side: Side::Right,
            rect: Rect {
                x: viewport.width / 2. + 8.,
                y: pad,
                width,
                height: (viewport.height - 2. * pad).max(0.),
            },
        },
    ];
    let mut out = layout(&input, viewport, &columns, options)?;
    let clearance = radius + 10.;
    out.displayed = 0;
    for (i, result) in out.labels.iter_mut().enumerate() {
        let Outcome::Placed {
            mut rect,
            mut leader,
        } = result.outcome
        else {
            continue;
        };
        let anchor = input[i].anchor;
        let lo = anchor.y.min(rect.y);
        let hi = anchor.y.max(rect.bottom());
        let nearest = if lo <= center.y && hi >= center.y {
            0.
        } else {
            (lo - center.y).abs().min((hi - center.y).abs())
        };
        let x = (clearance * clearance - nearest * nearest)
            .max(0.)
            .sqrt()
            .max((anchor.x - center.x).abs() + 8.);
        let right = input[i].side == Side::Right;
        let edge = center.x + if right { x } else { -x };
        let available = if right {
            viewport.width - pad - edge
        } else {
            edge - pad
        };
        if rect.width > available || available <= 0. {
            result.outcome = Outcome::Overflow {
                reason: OverflowReason::TooWide,
                available: Some(Size {
                    width: available.max(0.),
                    height: rect.height,
                }),
            };
            continue;
        }
        rect.x = if right { edge } else { edge - rect.width };
        leader.endpoint = Point {
            x: edge,
            y: rect.y + rect.height / 2.,
        };
        result.outcome = Outcome::Placed { rect, leader };
        out.displayed += 1;
    }
    Ok(out)
}
/// A terminal head is contained in the exact span, even for a very short feature.
/// Returns (tip, shoulder); strand orientation changes which true endpoint is tip.
pub fn arrow_head(start: f64, end: f64, reverse: bool, maximum: f64) -> Option<(f64, f64)> {
    if !start.is_finite()
        || !end.is_finite()
        || !maximum.is_finite()
        || end <= start
        || maximum <= 0.
    {
        return None;
    }
    let head = maximum.min((end - start) * 0.6);
    Some(if reverse {
        (start, start + head)
    } else {
        (end, end - head)
    })
}
/// Three inner feature lanes, wholly within the .65-normalized backbone radius.
pub fn feature_band(lane: usize) -> (f64, f64) {
    let outer = 0.615 - lane as f64 * 0.065;
    (outer - 0.04, outer)
}
