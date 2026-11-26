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
/// A short outward radial stem followed by a fan segment to horizontal text.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FanLeader {
    pub anchor: Point,
    pub knee: Point,
    pub endpoint: Point,
}
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum FanOutcome {
    Placed {
        rect: Rect,
        leader: FanLeader,
    },
    Overflow {
        reason: OverflowReason,
        available: Option<Size>,
    },
}
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FanLabelResult {
    pub id: u64,
    pub outcome: FanOutcome,
}
#[derive(Clone, Debug, PartialEq)]
pub struct FanLayout {
    /// Every original ID exactly once, in original input order.
    pub labels: Vec<FanLabelResult>,
    pub displayed: usize,
    pub considered: usize,
}
const FAN_STEM: f64 = 8.;
const FAN_CLEARANCE: f64 = 8.;
const FAN_PREFERRED: f64 = 18.;
const FAN_MIN_HEIGHT: f64 = 14.;
#[derive(Clone, Copy)]
struct FanEntry {
    index: usize,
    right: bool,
    ordinate: f64,
    preferred: f64,
    low: f64,
    high: f64,
    height: f64,
    max_edge: f64,
    ux: f64,
    uy: f64,
    knee: Point,
}
fn fan_entry(
    index: usize,
    label: &Label,
    viewport: Viewport,
    center: Point,
    radius: f64,
) -> Result<FanEntry, (OverflowReason, Option<Size>)> {
    let height = label.measured.height.max(FAN_MIN_HEIGHT);
    let available_height = viewport.height - 2. * viewport.padding;
    let dx = label.anchor.x - center.x;
    let dy = label.anchor.y - center.y;
    let distance = dx.hypot(dy);
    let ux = dx.abs() / distance;
    let uy = dy / distance;
    let right = dx >= 0.;
    let knee = Point {
        x: label.anchor.x + dx / distance * FAN_STEM,
        y: label.anchor.y + uy * FAN_STEM,
    };
    let knee_x = (knee.x - center.x).abs();
    let extent = if right {
        viewport.width - viewport.padding - center.x
    } else {
        center.x - viewport.padding
    };
    // Keep text and the fan segment on their actual anchor's side. A minimum
    // two-pixel horizontal advance prevents an inward folded-back terminal fan.
    let available_width = (extent - knee_x - 2.).max(0.);
    let available = Some(Size {
        width: available_width,
        height: available_height,
    });
    if height > available_height {
        return Err((OverflowReason::TooTall, available));
    }
    if label.measured.width > available_width {
        return Err((OverflowReason::TooWide, available));
    }
    let max_edge = extent - label.measured.width;
    let mut low = viewport.padding + height / 2.;
    let mut high = viewport.height - viewport.padding - height / 2.;
    let clear_radius = radius + FAN_CLEARANCE;
    if max_edge < clear_radius {
        // A wide label can still fit near a pole. Its entire measured box,
        // rather than just its midpoint, must clear the circle's contour.
        let nearest = (clear_radius * clear_radius - max_edge * max_edge)
            .max(0.)
            .sqrt();
        if uy <= 0. {
            high = high.min(center.y - nearest - height / 2.);
        } else {
            low = low.max(center.y + nearest + height / 2.);
        }
    }
    // The knee is exactly radial. Requiring the following segment to start in
    // its outward tangent half-plane analytically excludes the backbone disk.
    // Density displacement is allowed only while that connection remains safe.
    let horizontal_progress = ux * (max_edge - knee_x);
    if uy > 1e-12 {
        low = low.max(knee.y - horizontal_progress / uy);
    } else if uy < -1e-12 {
        high = high.min(knee.y - horizontal_progress / uy);
    }
    if low > high {
        return Err((
            OverflowReason::Dense,
            Some(Size {
                width: available_width,
                height: (high - low + height).max(0.),
            }),
        ));
    }
    Ok(FanEntry {
        index,
        right,
        ordinate: uy,
        preferred: label.anchor.y + uy * FAN_PREFERRED,
        low,
        high,
        height,
        max_edge,
        ux,
        uy,
        knee,
    })
}
fn fan_fits(entries: &[FanEntry], gap: f64) -> bool {
    let mut previous_bottom = f64::NEG_INFINITY;
    for entry in entries {
        let middle = entry.low.max(previous_bottom + gap + entry.height / 2.);
        if middle > entry.high {
            return false;
        }
        previous_bottom = middle + entry.height / 2.;
    }
    true
}
#[derive(Clone, Copy)]
struct FanBlock {
    first: usize,
    end: usize,
    sum: f64,
    low: f64,
    high: f64,
}
impl FanBlock {
    fn value(self) -> f64 {
        (self.sum / (self.end - self.first) as f64)
            .max(self.low)
            .min(self.high)
    }
}
/// Bounded isotonic packing minimizes squared row displacement for the admitted
/// labels. Unlike a top-to-bottom push, it fans a dense group in both directions.
fn fan_rows(entries: &[FanEntry], gap: f64) -> Vec<f64> {
    let mut offsets = Vec::with_capacity(entries.len());
    let mut offset = 0.;
    let mut blocks: Vec<FanBlock> = Vec::with_capacity(entries.len());
    for (i, entry) in entries.iter().enumerate() {
        let middle = offset + entry.height / 2.;
        offsets.push(middle);
        offset += entry.height + gap;
        blocks.push(FanBlock {
            first: i,
            end: i + 1,
            sum: entry.preferred - middle,
            low: entry.low - middle,
            high: entry.high - middle,
        });
        while blocks.len() >= 2 {
            let last = blocks.len() - 1;
            if blocks[last - 1].value() <= blocks[last].value() {
                break;
            }
            let right = blocks.pop().expect("two existing blocks");
            let left = blocks.pop().expect("two existing blocks");
            blocks.push(FanBlock {
                first: left.first,
                end: right.end,
                sum: left.sum + right.sum,
                low: left.low.max(right.low),
                high: left.high.min(right.high),
            });
        }
    }
    let mut rows = vec![0.; entries.len()];
    for block in blocks {
        let value = block.value();
        for i in block.first..block.end {
            rows[i] = value + offsets[i];
        }
    }
    rows
}
fn segment_distance(center: Point, a: Point, b: Point) -> f64 {
    let dx = b.x - a.x;
    let dy = b.y - a.y;
    let length = dx * dx + dy * dy;
    let t = if length == 0. {
        0.
    } else {
        (((center.x - a.x) * dx + (center.y - a.y) * dy) / length).clamp(0., 1.)
    };
    (a.x + t * dx - center.x).hypot(a.y + t * dy - center.y)
}
fn fan_position(
    entry: FanEntry,
    middle: f64,
    label: &Label,
    center: Point,
    radius: f64,
    viewport: Viewport,
) -> Option<(Rect, FanLeader)> {
    if middle < entry.low || middle > entry.high {
        return None;
    }
    let nearest_y = ((middle - center.y).abs() - entry.height / 2.).max(0.);
    let clear_radius = radius + FAN_CLEARANCE;
    let circle_edge = (clear_radius * clear_radius - nearest_y * nearest_y)
        .max(0.)
        .sqrt();
    let knee_x = (entry.knee.x - center.x).abs();
    let tangent = if entry.ux > 1e-12 {
        knee_x - entry.uy * (middle - entry.knee.y) / entry.ux
    } else {
        knee_x
    };
    let preferred = (label.anchor.x - center.x).abs() + entry.ux * FAN_PREFERRED + 4.;
    let edge = preferred
        .max(circle_edge)
        .max(knee_x + 2.)
        .max(tangent)
        .min(entry.max_edge);
    let rect = Rect {
        x: if entry.right {
            center.x + edge
        } else {
            center.x - edge - label.measured.width
        },
        y: middle - entry.height / 2.,
        width: label.measured.width,
        height: entry.height,
    };
    let leader = FanLeader {
        anchor: label.anchor,
        knee: entry.knee,
        endpoint: Point {
            x: if entry.right { rect.x } else { rect.right() },
            y: middle,
        },
    };
    let near = Point {
        x: center.x.clamp(rect.x, rect.right()),
        y: center.y.clamp(rect.y, rect.bottom()),
    };
    let epsilon = radius * 1e-10 + 1e-7;
    if rect.x < viewport.padding
        || rect.y < viewport.padding
        || rect.right() > viewport.width - viewport.padding
        || rect.bottom() > viewport.height - viewport.padding
        || rect.right() <= rect.x
        || rect.bottom() <= rect.y
        || (near.x - center.x).hypot(near.y - center.y) < radius
        || segment_distance(center, leader.anchor, leader.knee) + epsilon < radius
        || segment_distance(center, leader.knee, leader.endpoint) + epsilon < radius
    {
        return None;
    }
    Some((rect, leader))
}
fn segment_hits_rect(a: Point, b: Point, rect: Rect) -> bool {
    let mut low: f64 = 0.;
    let mut high: f64 = 1.;
    for (origin, delta, minimum, maximum) in [
        (a.x, b.x - a.x, rect.x, rect.right()),
        (a.y, b.y - a.y, rect.y, rect.bottom()),
    ] {
        if delta.abs() < 1e-12 {
            if origin < minimum || origin > maximum {
                return false;
            }
        } else {
            let first = (minimum - origin) / delta;
            let second = (maximum - origin) / delta;
            low = low.max(first.min(second));
            high = high.min(first.max(second));
            if low > high {
                return false;
            }
        }
    }
    true
}
fn proper_cross(a: Point, b: Point, c: Point, d: Point) -> bool {
    let cross = |ax: f64, ay: f64, bx: f64, by: f64| ax * by - ay * bx;
    let dx = b.x - a.x;
    let dy = b.y - a.y;
    let ex = d.x - c.x;
    let ey = d.y - c.y;
    let denominator = cross(dx, dy, ex, ey);
    if denominator.abs() < 1e-9 {
        return false;
    }
    let t = cross(c.x - a.x, c.y - a.y, ex, ey) / denominator;
    let u = cross(c.x - a.x, c.y - a.y, dx, dy) / denominator;
    t > 1e-8 && t < 1. - 1e-8 && u > 1e-8 && u < 1. - 1e-8
}
fn fan_conflict(a: Rect, first: FanLeader, b: Rect, second: FanLeader) -> bool {
    if crate::overlaps(a, b) {
        return true;
    }
    let first_segments = [(first.anchor, first.knee), (first.knee, first.endpoint)];
    let second_segments = [(second.anchor, second.knee), (second.knee, second.endpoint)];
    first_segments
        .iter()
        .any(|&(a, c)| segment_hits_rect(a, c, b))
        || second_segments
            .iter()
            .any(|&(c, d)| segment_hits_rect(c, d, a))
        || first_segments.iter().any(|&(a, b)| {
            second_segments
                .iter()
                .any(|&(c, d)| proper_cross(a, b, c, d))
        })
}
/// Horizontal restriction labels fanned around their true circular anchors.
///
/// Priority/stable-rank selection is camera-independent. Within each side,
/// angular order is retained, measured boxes use at least14px rows plus at least
///2px gaps, and an isotonic fit favors each anchor's nearby radial position.
/// X follows the local circle contour rather than one common distant column.
/// Every leader starts with an exact8px outward radial segment; its fan segment
/// is admitted only outside the disk. Unsafe/crowded lower-priority labels remain
/// explicit overflow. Admission and conflict checks are bounded by1024²; this
/// greedy policy does not claim maximum cardinality or globally optimal routing.
pub fn fan_out(
    labels: &[Label],
    viewport: Viewport,
    center: Point,
    radius: f64,
    options: &Options,
) -> Result<FanLayout, Error> {
    crate::preflight(labels, viewport, &[], options)?;
    if !radius.is_finite()
        || radius <= 0.
        || radius > crate::MAX_VIEWPORT
        || !crate::coordinate(center.x)
        || !crate::coordinate(center.y)
    {
        return Err(Error::InvalidOptions);
    }
    let charge = labels
        .len()
        .checked_mul(std::mem::size_of::<FanLabelResult>() + std::mem::size_of::<usize>())
        .and_then(|n| {
            options
                .max_displayed
                .checked_mul(
                    3 * std::mem::size_of::<FanEntry>()
                        + 2 * std::mem::size_of::<FanBlock>()
                        + 4 * std::mem::size_of::<usize>(),
                )
                .and_then(|s| n.checked_add(s))
        })
        .ok_or(Error::MemoryLimit)?;
    if charge > crate::MAX_MEMORY_CHARGE {
        return Err(Error::MemoryLimit);
    }
    for label in labels {
        let distance = (label.anchor.x - center.x).hypot(label.anchor.y - center.y);
        if distance == 0. || distance + radius * 1e-10 + 1e-7 < radius {
            return Err(Error::InvalidLabel(label.id));
        }
    }
    let mut order: Vec<usize> = (0..labels.len()).collect();
    order.sort_unstable_by_key(|&i| labels[i].id);
    for pair in order.windows(2) {
        if labels[pair[0]].id == labels[pair[1]].id {
            return Err(Error::DuplicateId(labels[pair[0]].id));
        }
    }
    order.sort_unstable_by(|&a, &b| {
        labels[b]
            .priority
            .cmp(&labels[a].priority)
            .then_with(|| labels[a].stable_rank.cmp(&labels[b].stable_rank))
            .then_with(|| labels[a].id.cmp(&labels[b].id))
    });
    let mut selected = Vec::with_capacity(options.max_displayed.min(labels.len()));
    let mut start = 0;
    while start < order.len() && selected.len() < options.max_displayed {
        let mut end = start + 1;
        while end < order.len() && labels[order[end]].priority == labels[order[start]].priority {
            end += 1;
        }
        let count = end - start;
        let take = count.min(options.max_displayed - selected.len());
        for k in 0..take {
            let relative = (((2 * k + 1) as u64 * count as u64) / (2 * take) as u64) as usize;
            selected.push(order[start + relative]);
        }
        start = end;
    }
    let mut result = FanLayout {
        labels: labels
            .iter()
            .map(|label| FanLabelResult {
                id: label.id,
                outcome: FanOutcome::Overflow {
                    reason: OverflowReason::DisplayBudget,
                    available: None,
                },
            })
            .collect(),
        displayed: 0,
        considered: selected.len(),
    };
    let gap = options.gap.max(2.);
    let mut sides: [Vec<FanEntry>; 2] = std::array::from_fn(|_| Vec::new());
    for &index in &selected {
        let entry = match fan_entry(index, &labels[index], viewport, center, radius) {
            Ok(entry) => entry,
            Err((reason, available)) => {
                result.labels[index].outcome = FanOutcome::Overflow { reason, available };
                continue;
            }
        };
        let side = &mut sides[usize::from(entry.right)];
        let at = side
            .binary_search_by(|other| {
                other
                    .ordinate
                    .total_cmp(&entry.ordinate)
                    .then_with(|| {
                        labels[other.index]
                            .stable_rank
                            .cmp(&labels[index].stable_rank)
                    })
                    .then_with(|| labels[other.index].id.cmp(&labels[index].id))
            })
            .unwrap_or_else(|at| at);
        side.insert(at, entry);
        result.labels[index].outcome = FanOutcome::Overflow {
            reason: OverflowReason::Dense,
            available: Some(Size {
                width: (entry.max_edge + labels[index].measured.width
                    - (entry.knee.x - center.x).abs()
                    - 2.)
                    .max(0.),
                height: (entry.high - entry.low + entry.height).max(0.),
            }),
        };
        if !fan_fits(side, gap) {
            side.remove(at);
        }
    }
    for side in &sides {
        let rows = fan_rows(side, gap);
        for (&entry, middle) in side.iter().zip(rows) {
            result.labels[entry.index].outcome = match fan_position(
                entry,
                middle,
                &labels[entry.index],
                center,
                radius,
                viewport,
            ) {
                Some((rect, leader)) => FanOutcome::Placed { rect, leader },
                None => FanOutcome::Overflow {
                    reason: OverflowReason::Precision,
                    available: None,
                },
            };
        }
    }
    // Retain higher-priority routing when density makes a lower-priority fan
    // intersect another label or leader. Shared anchors/trunks are intentional;
    // only proper segment crossings count. Positions are not repeatedly rebuilt.
    for (rank, &higher) in selected.iter().enumerate() {
        let FanOutcome::Placed {
            rect: a,
            leader: first,
        } = result.labels[higher].outcome
        else {
            continue;
        };
        for &lower in &selected[rank + 1..] {
            let FanOutcome::Placed {
                rect: b,
                leader: second,
            } = result.labels[lower].outcome
            else {
                continue;
            };
            if fan_conflict(a, first, b, second) {
                result.labels[lower].outcome = FanOutcome::Overflow {
                    reason: OverflowReason::Dense,
                    available: None,
                };
            }
        }
    }
    result.displayed = result
        .labels
        .iter()
        .filter(|label| matches!(label.outcome, FanOutcome::Placed { .. }))
        .count();
    Ok(result)
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
