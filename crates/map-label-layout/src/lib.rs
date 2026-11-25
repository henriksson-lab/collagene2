//! Measured-label layout, independent of rendering and biological coordinates.
//! Original anchors are retained; overflow is evidence, never silent deletion.
#![doc = include_str!("../README.md")]

pub mod radial;

use std::{cmp::Ordering, fmt};

pub const MAX_INPUTS: usize = 131_072;
pub const MAX_DISPLAYED: usize = 1024;
pub const MAX_VIEWPORT: f64 = 1_000_000.;
pub const MAX_COORDINATE: f64 = 1_000_000_000.;
pub const MAX_MEMORY_CHARGE: usize = 64 * 1024 * 1024;

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Point {
    pub x: f64,
    pub y: f64,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Size {
    pub width: f64,
    pub height: f64,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Rect {
    pub x: f64,
    pub y: f64,
    pub width: f64,
    pub height: f64,
}

impl Rect {
    pub fn right(self) -> f64 {
        self.x + self.width
    }
    pub fn bottom(self) -> f64 {
        self.y + self.height
    }
    pub fn size(self) -> Size {
        Size {
            width: self.width,
            height: self.height,
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Viewport {
    pub width: f64,
    pub height: f64,
    pub padding: f64,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub enum Side {
    Left,
    Right,
    Above,
    Below,
}

impl Side {
    fn vertical(self) -> bool {
        matches!(self, Self::Left | Self::Right)
    }
}

/// Disjoint caller-owned label region inside the padded viewport. Left/Right
/// flow vertically; Above/Below flow horizontally. No map exclusion is inferred.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Column {
    pub side: Side,
    pub rect: Rect,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Label {
    pub id: u64,
    /// Camera-independent caller rank used only for tied-priority LOD sampling.
    pub stable_rank: u64,
    pub anchor: Point,
    pub preferred_center: Point,
    pub side: Side,
    /// Actual caller-measured full text box, including desired text padding.
    pub measured: Size,
    /// Explicit presentation priority only. Higher values enter the budget first.
    pub priority: u32,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Options {
    pub max_inputs: usize,
    pub max_displayed: usize,
    pub gap: f64,
}

impl Default for Options {
    fn default() -> Self {
        Self {
            max_inputs: 32_768,
            max_displayed: 128,
            gap: 4.,
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum OverflowReason {
    DisplayBudget,
    MissingColumn,
    TooWide,
    TooTall,
    Dense,
    /// The finite measured box cannot be represented inside its region without
    /// a boundary/collision error at the supplied floating-point coordinates.
    Precision,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Leader {
    /// The supplied projected anchor, bit-for-bit unchanged (even offscreen).
    pub anchor: Point,
    /// A point on the placed text-box edge facing the map side.
    pub endpoint: Point,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Outcome {
    Placed {
        rect: Rect,
        leader: Leader,
    },
    Overflow {
        reason: OverflowReason,
        available: Option<Size>,
    },
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct LabelResult {
    pub id: u64,
    pub outcome: Outcome,
}

#[derive(Clone, Debug, PartialEq)]
pub struct Layout {
    /// Every supplied ID exactly once, in original input order.
    pub labels: Vec<LabelResult>,
    pub displayed: usize,
    /// Labels considered after explicit priority/stable-rank budget selection. This
    /// can exceed displayed when some measured boxes do not fit their columns.
    pub considered: usize,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub enum Error {
    InvalidOptions,
    InputLimit,
    MemoryLimit,
    InvalidViewport,
    InvalidColumn,
    DuplicateSide,
    OverlappingColumns,
    InvalidLabel(u64),
    DuplicateId(u64),
}

impl fmt::Display for Error {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "map label layout: {self:?}")
    }
}
impl std::error::Error for Error {}

fn positive(value: f64, cap: f64) -> bool {
    value.is_finite() && value > 0. && value <= cap
}
fn coordinate(value: f64) -> bool {
    value.is_finite() && value.abs() <= MAX_COORDINATE
}
fn overlaps(a: Rect, b: Rect) -> bool {
    a.x < b.right() && b.x < a.right() && a.y < b.bottom() && b.y < a.bottom()
}

/// Allocation-free count, finite-coordinate and geometry validation. Duplicate
/// ID detection uses a bounded sorted index after this resource preflight.
pub fn preflight(
    labels: &[Label],
    viewport: Viewport,
    columns: &[Column],
    options: &Options,
) -> Result<(), Error> {
    if options.max_inputs > MAX_INPUTS
        || options.max_displayed > MAX_DISPLAYED
        || !options.gap.is_finite()
        || options.gap < 0.
        || options.gap > MAX_VIEWPORT
    {
        return Err(Error::InvalidOptions);
    }
    if labels.len() > options.max_inputs {
        return Err(Error::InputLimit);
    }
    // Output, full ordering index, selected index and per-column selected lists.
    // sort_unstable allocates no scratch. Each selected ID belongs to one column.
    let per_input = std::mem::size_of::<LabelResult>() + std::mem::size_of::<usize>();
    let bytes = labels
        .len()
        .checked_mul(per_input)
        .and_then(|n| {
            options
                .max_displayed
                .checked_mul(2 * std::mem::size_of::<usize>())
                .and_then(|s| n.checked_add(s))
        })
        .ok_or(Error::MemoryLimit)?;
    if bytes > MAX_MEMORY_CHARGE {
        return Err(Error::MemoryLimit);
    }
    if !positive(viewport.width, MAX_VIEWPORT)
        || !positive(viewport.height, MAX_VIEWPORT)
        || !viewport.padding.is_finite()
        || viewport.padding < 0.
        || viewport.padding * 2. >= viewport.width
        || viewport.padding * 2. >= viewport.height
    {
        return Err(Error::InvalidViewport);
    }
    if columns.len() > 4 {
        return Err(Error::InvalidColumn);
    }
    for (i, column) in columns.iter().enumerate() {
        let r = column.rect;
        if !coordinate(r.x)
            || !coordinate(r.y)
            || !positive(r.width, MAX_VIEWPORT)
            || !positive(r.height, MAX_VIEWPORT)
            || r.x < viewport.padding
            || r.y < viewport.padding
            || r.right() > viewport.width - viewport.padding
            || r.bottom() > viewport.height - viewport.padding
        {
            return Err(Error::InvalidColumn);
        }
        for previous in &columns[..i] {
            if previous.side == column.side {
                return Err(Error::DuplicateSide);
            }
            if overlaps(previous.rect, r) {
                return Err(Error::OverlappingColumns);
            }
        }
    }
    for label in labels {
        if !coordinate(label.anchor.x)
            || !coordinate(label.anchor.y)
            || !coordinate(label.preferred_center.x)
            || !coordinate(label.preferred_center.y)
            || !positive(label.measured.width, MAX_VIEWPORT)
            || !positive(label.measured.height, MAX_VIEWPORT)
        {
            return Err(Error::InvalidLabel(label.id));
        }
    }
    Ok(())
}

fn spatial(a: &Label, b: &Label) -> Ordering {
    a.side
        .cmp(&b.side)
        .then_with(|| {
            if a.side.vertical() {
                a.preferred_center
                    .y
                    .total_cmp(&b.preferred_center.y)
                    .then_with(|| a.preferred_center.x.total_cmp(&b.preferred_center.x))
            } else {
                a.preferred_center
                    .x
                    .total_cmp(&b.preferred_center.x)
                    .then_with(|| a.preferred_center.y.total_cmp(&b.preferred_center.y))
            }
        })
        .then_with(|| a.id.cmp(&b.id))
}

fn overflow(reason: OverflowReason, available: Option<Size>) -> Outcome {
    Outcome::Overflow { reason, available }
}

/// Deterministic bounded placement. Selection uses explicit presentation priority
/// and evenly spaced caller-provided stable ranks within each equal-priority tier. No feature
/// order is interpreted as biological importance. Box fitting is greedy and does
/// not claim maximum cardinality, optimal displacement or uncrossed leader lines.
pub fn layout(
    labels: &[Label],
    viewport: Viewport,
    columns: &[Column],
    options: &Options,
) -> Result<Layout, Error> {
    preflight(labels, viewport, columns, options)?;
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
    let mut begin = 0;
    while begin < order.len() && selected.len() < options.max_displayed {
        let priority = labels[order[begin]].priority;
        let mut end = begin + 1;
        while end < order.len() && labels[order[end]].priority == priority {
            end += 1;
        }
        let count = end - begin;
        let take = count.min(options.max_displayed - selected.len());
        for k in 0..take {
            // Midpoints of equal rank strata; all arithmetic fits u64 under caps.
            let relative = (((2 * k + 1) as u64 * count as u64) / (2 * take) as u64) as usize;
            selected.push(order[begin + relative]);
        }
        begin = end;
    }
    let mut output = Layout {
        labels: labels
            .iter()
            .map(|label| LabelResult {
                id: label.id,
                outcome: overflow(OverflowReason::DisplayBudget, None),
            })
            .collect(),
        displayed: 0,
        considered: selected.len(),
    };
    let mut accepted: [Vec<usize>; 4] = std::array::from_fn(|index| {
        let count = columns.get(index).map_or(0, |column| {
            selected
                .iter()
                .filter(|&&i| labels[i].side == column.side)
                .count()
        });
        Vec::with_capacity(count)
    });
    let mut used = [0.; 4];
    for &i in &selected {
        let label = &labels[i];
        let Some(column_index) = columns.iter().position(|column| column.side == label.side) else {
            output.labels[i].outcome = overflow(OverflowReason::MissingColumn, None);
            continue;
        };
        let column = columns[column_index];
        let available = Some(column.rect.size());
        let reason = if label.measured.width > column.rect.width {
            Some(OverflowReason::TooWide)
        } else if label.measured.height > column.rect.height {
            Some(OverflowReason::TooTall)
        } else {
            None
        };
        if let Some(reason) = reason {
            output.labels[i].outcome = overflow(reason, available);
            continue;
        }
        let extent = if column.side.vertical() {
            label.measured.height
        } else {
            label.measured.width
        };
        let limit = if column.side.vertical() {
            column.rect.height
        } else {
            column.rect.width
        };
        let total = used[column_index]
            + extent
            + if accepted[column_index].is_empty() {
                0.
            } else {
                options.gap
            };
        if total > limit {
            output.labels[i].outcome = overflow(OverflowReason::Dense, available);
            continue;
        }
        used[column_index] = total;
        accepted[column_index].push(i);
    }
    for (index, column) in columns.iter().enumerate() {
        let vertical = column.side.vertical();
        let origin = if vertical {
            column.rect.y
        } else {
            column.rect.x
        };
        let limit = if vertical {
            column.rect.bottom()
        } else {
            column.rect.right()
        };
        let members = &mut accepted[index];
        members.sort_unstable_by(|&a, &b| spatial(&labels[a], &labels[b]));
        let mut cursor = origin;
        for &i in members.iter() {
            let label = labels[i];
            let extent = if vertical {
                label.measured.height
            } else {
                label.measured.width
            };
            let preferred = if vertical {
                label.preferred_center.y
            } else {
                label.preferred_center.x
            };
            let axis_size = if vertical {
                column.rect.height
            } else {
                column.rect.width
            };
            let at = (preferred - extent / 2.)
                .clamp(origin, origin + (axis_size - extent))
                .max(cursor);
            let rect = if vertical {
                Rect {
                    x: (label.preferred_center.x - label.measured.width / 2.).clamp(
                        column.rect.x,
                        column.rect.x + (column.rect.width - label.measured.width),
                    ),
                    y: at,
                    width: label.measured.width,
                    height: label.measured.height,
                }
            } else {
                Rect {
                    x: at,
                    y: (label.preferred_center.y - label.measured.height / 2.).clamp(
                        column.rect.y,
                        column.rect.y + (column.rect.height - label.measured.height),
                    ),
                    width: label.measured.width,
                    height: label.measured.height,
                }
            };
            output.labels[i].outcome = Outcome::Placed {
                rect,
                leader: Leader {
                    anchor: label.anchor,
                    endpoint: label.anchor,
                },
            };
            cursor = at + extent + options.gap;
        }
        // Backward sweep fits the tail while retaining order. Total measured
        // extents+gaps were admitted first, so an in-bounds solution exists.
        cursor = limit;
        for &i in members.iter().rev() {
            let Outcome::Placed { mut rect, .. } = output.labels[i].outcome else {
                unreachable!()
            };
            if vertical {
                rect.y = rect.y.min(cursor - rect.height).max(column.rect.y);
                cursor = rect.y - options.gap;
            } else {
                rect.x = rect.x.min(cursor - rect.width).max(column.rect.x);
                cursor = rect.x - options.gap;
            }
            let anchor = labels[i].anchor;
            let endpoint = match column.side {
                Side::Left => Point {
                    x: rect.right(),
                    y: anchor.y.clamp(rect.y, rect.bottom()),
                },
                Side::Right => Point {
                    x: rect.x,
                    y: anchor.y.clamp(rect.y, rect.bottom()),
                },
                Side::Above => Point {
                    x: anchor.x.clamp(rect.x, rect.right()),
                    y: rect.bottom(),
                },
                Side::Below => Point {
                    x: anchor.x.clamp(rect.x, rect.right()),
                    y: rect.y,
                },
            };
            output.labels[i].outcome = Outcome::Placed {
                rect,
                leader: Leader { anchor, endpoint },
            };
        }
        // Floating-point addition/subtraction can lose a tiny box at a large
        // origin, or round an exact-fit boundary outward. Never publish an
        // overlapping/out-of-bounds rectangle as successfully placed. This
        // linear adjacent check retains explicit evidence rather than silently
        // accepting an epsilon violation or running an unbounded repair search.
        let mut previous_end = None;
        for &i in members.iter() {
            let Outcome::Placed { rect, .. } = output.labels[i].outcome else {
                unreachable!()
            };
            let start = if vertical { rect.y } else { rect.x };
            let end = if vertical {
                rect.bottom()
            } else {
                rect.right()
            };
            let valid = rect.x >= column.rect.x
                && rect.y >= column.rect.y
                && rect.right() <= column.rect.right()
                && rect.bottom() <= column.rect.bottom()
                && rect.right() > rect.x
                && rect.bottom() > rect.y
                && previous_end.is_none_or(|previous| start >= previous + options.gap);
            if valid {
                previous_end = Some(end);
                output.displayed += 1;
            } else {
                output.labels[i].outcome =
                    overflow(OverflowReason::Precision, Some(column.rect.size()));
            }
        }
    }
    Ok(output)
}
