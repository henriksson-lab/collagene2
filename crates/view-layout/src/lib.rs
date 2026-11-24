//! Static interval layout, independent of pixels, sequences and biological models.
//!
//! [`pack`] validates borrowed input and budgets worst-case work before allocating.
//! Results preserve input order. Processing order is normalized interval-vector
//! order followed by stable integer ID, so input permutations cannot change lanes.
//! First-fit packing is deterministic, not an optimal graph coloring algorithm.
#![doc = include_str!("../README.md")]

use std::{cmp::Ordering, collections::BTreeMap, fmt};

pub const MAX_ITEMS: usize = 131072;
pub const MAX_INTERVALS: usize = 524288;
pub const MAX_LANES: usize = 256;
pub const MAX_WORK: u64 = 2_000_000_000;
pub const MAX_MEMORY_CHARGE: u64 = 256 * 1024 * 1024;

/// Half-open source interval. Public DTOs are always checked again by `pack`.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub struct Interval {
    pub start: u64,
    pub end: u64,
}

impl Interval {
    /// Construct without wrapping a start + length computation.
    pub fn from_start_length(start: u64, length: u64) -> Result<Self, Error> {
        if length == 0 {
            return Err(Error::EmptyInterval);
        }
        let end = start.checked_add(length).ok_or(Error::CoordinateOverflow)?;
        Ok(Self { start, end })
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Item {
    pub id: u64,
    /// Joined / circular-split spans belong to this one item and one lane.
    pub spans: Vec<Interval>,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Limits {
    pub max_items: usize,
    pub max_intervals: usize,
    pub max_lanes: usize,
    /// Conservative work units, not wall-clock milliseconds. See `work_bound`.
    pub max_work: u64,
}

impl Default for Limits {
    fn default() -> Self {
        Self {
            max_items: 32768,
            max_intervals: 131072,
            max_lanes: 64,
            max_work: 128_000_000,
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Placement {
    /// Zero-based lane; every span of this item occupies this lane.
    Lane(usize),
    /// Retained evidence that this whole item did not fit the requested lane cap.
    Overflow,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Assignment {
    pub id: u64,
    /// Sorted disjoint nonadjacent intervals; no biological meaning is inferred.
    pub spans: Vec<Interval>,
    pub placement: Placement,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Layout {
    /// Every input ID appears exactly once, in input order (including overflow).
    pub assignments: Vec<Assignment>,
    /// Used lanes, not the requested cap. Empty input / zero cap yields zero.
    pub lane_count: usize,
    /// Conservative preflight work bound charged even if actual work is smaller.
    pub work_bound: u64,
}

impl Layout {
    pub fn overflow_ids(&self) -> impl Iterator<Item = u64> + '_ {
        self.assignments
            .iter()
            .filter_map(|a| (a.placement == Placement::Overflow).then_some(a.id))
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub enum Error {
    InvalidLimits,
    ItemLimit,
    IntervalLimit,
    WorkLimit { required: u64, allowed: u64 },
    MemoryLimit { required: u64, allowed: u64 },
    ArithmeticOverflow,
    DuplicateId(u64),
    MissingIntervals(u64),
    EmptyInterval,
    InvalidInterval { id: u64, index: usize },
    CoordinateOverflow,
}

impl fmt::Display for Error {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidLimits => f.write_str("limits exceed hard library caps"),
            Self::ItemLimit => f.write_str("item count exceeds the requested limit"),
            Self::IntervalLimit => f.write_str("interval count exceeds the requested limit"),
            Self::WorkLimit { required, allowed } => {
                write!(f, "layout work bound {required} exceeds {allowed}")
            }
            Self::MemoryLimit { required, allowed } => {
                write!(f, "layout memory charge {required} exceeds {allowed}")
            }
            Self::ArithmeticOverflow => f.write_str("resource-bound arithmetic overflow"),
            Self::DuplicateId(id) => write!(f, "duplicate item ID {id}"),
            Self::MissingIntervals(id) => write!(f, "item {id} has no intervals"),
            Self::EmptyInterval => f.write_str("interval length must be positive"),
            Self::InvalidInterval { id, index } => write!(
                f,
                "item {id} interval {index} is empty, reversed or outside the source"
            ),
            Self::CoordinateOverflow => f.write_str("interval end exceeds u64 coordinates"),
        }
    }
}

impl std::error::Error for Error {}

fn add(a: u64, b: u64) -> Result<u64, Error> {
    a.checked_add(b).ok_or(Error::ArithmeticOverflow)
}

fn mul(a: u64, b: u64) -> Result<u64, Error> {
    a.checked_mul(b).ok_or(Error::ArithmeticOverflow)
}

fn log_bound(n: u64) -> u64 {
    if n <= 1 {
        0
    } else {
        u64::from(64 - (n - 1).leading_zeros())
    }
}

fn counts(items: &[Item], limits: &Limits) -> Result<(u64, u64, u64), Error> {
    if limits.max_items > MAX_ITEMS
        || limits.max_intervals > MAX_INTERVALS
        || limits.max_lanes > MAX_LANES
        || limits.max_work > MAX_WORK
    {
        return Err(Error::InvalidLimits);
    }
    if items.len() > limits.max_items {
        return Err(Error::ItemLimit);
    }
    let mut total = 0usize;
    let mut largest = 0usize;
    for item in items {
        total = total
            .checked_add(item.spans.len())
            .ok_or(Error::ArithmeticOverflow)?;
        if total > limits.max_intervals {
            return Err(Error::IntervalLimit);
        }
        largest = largest.max(item.spans.len());
    }
    Ok((items.len() as u64, total as u64, largest as u64))
}

/// Conservative storage charge including assignments, original-sized normalized
/// span vectors, one largest-vector merge scratch, order + order scratch, lane
/// headers and up to one B-tree entry per original span. No allocation occurs.
/// B-tree charge is 256 bytes per stored interval; allocator bookkeeping and
/// caller-owned inputs are outside this library's managed-storage accounting.
pub fn memory_charge(items: &[Item], limits: &Limits) -> Result<u64, Error> {
    let (n, s, largest) = counts(items, limits)?;
    let mut bytes = mul(n, std::mem::size_of::<Assignment>() as u64)?;
    bytes = add(
        bytes,
        mul(add(s, largest)?, std::mem::size_of::<Interval>() as u64)?,
    )?;
    bytes = add(bytes, mul(n, (2 * std::mem::size_of::<usize>()) as u64)?)?;
    bytes = add(
        bytes,
        mul(
            limits.max_lanes.min(items.len()) as u64,
            std::mem::size_of::<BTreeMap<u64, u64>>() as u64,
        )?,
    )?;
    bytes = add(bytes, mul(s, 256)?)?;
    if bytes > MAX_MEMORY_CHARGE {
        return Err(Error::MemoryLimit {
            required: bytes,
            allowed: MAX_MEMORY_CHARGE,
        });
    }
    Ok(bytes)
}

/// Preallocation budget for validation, explicit mergesorts, and ordered disjoint
/// interval-index operations. Original spans are charged before normalization.
/// An index operation costs 4*(ceil(log2(S+1))+1) work units; an item probes at
/// most two neighbors per span per candidate lane, then inserts each span once.
/// These are conservative logical work units, not instruction/time guarantees.
pub fn work_bound(items: &[Item], limits: &Limits) -> Result<u64, Error> {
    let (n, s, largest) = counts(items, limits)?;
    memory_charge(items, limits)?;
    let mut work = add(n, mul(2, s)?)?;
    work = add(work, mul(mul(4, n)?, log_bound(n))?)?; // duplicate-ID sort
    work = add(work, mul(mul(4, s)?, log_bound(largest))?)?;
    // Each merge comparison consumes one key. Lexicographic comparisons inspect
    // no more spans than the consumed key, so each level charges O(N+S), not N*S.
    work = add(work, mul(mul(4, add(n, s)?)?, log_bound(n))?)?;
    let index_cost = mul(4, add(log_bound(add(s, 1)?), 1)?)?;
    let probes = mul(mul(2, limits.max_lanes as u64)?, s)?;
    work = add(work, mul(add(probes, s)?, index_cost)?)?;
    work = add(work, mul(n, limits.max_lanes as u64)?)?;
    if work > limits.max_work {
        return Err(Error::WorkLimit {
            required: work,
            allowed: limits.max_work,
        });
    }
    Ok(work)
}

fn merge_sort<T: Copy>(values: &mut [T], scratch: &mut [T], compare: impl Fn(&T, &T) -> Ordering) {
    let n = values.len();
    let mut width = 1;
    while width < n {
        let mut start = 0;
        while start < n {
            let middle = (start + width).min(n);
            let end = (middle + width).min(n);
            let (mut left, mut right) = (start, middle);
            for output in &mut scratch[start..end] {
                if left < middle
                    && (right == end || compare(&values[left], &values[right]) != Ordering::Greater)
                {
                    *output = values[left];
                    left += 1;
                } else {
                    *output = values[right];
                    right += 1;
                }
            }
            start = end;
        }
        values.copy_from_slice(&scratch[..n]);
        width *= 2; // n is bounded by the library's checked item/span caps.
    }
}

fn normalize(spans: &[Interval], scratch: &mut [Interval]) -> Vec<Interval> {
    let mut out = spans.to_vec();
    merge_sort(&mut out, scratch, Interval::cmp);
    let mut kept = 0;
    for read in 0..out.len() {
        let next = out[read];
        if kept > 0 && next.start <= out[kept - 1].end {
            out[kept - 1].end = out[kept - 1].end.max(next.end);
        } else {
            out[kept] = next;
            kept += 1;
        }
    }
    out.truncate(kept);
    out
}

fn fits(index: &BTreeMap<u64, u64>, spans: &[Interval]) -> bool {
    spans.iter().all(|span| {
        // A lane index is disjoint: only the predecessor and successor can
        // overlap a new span. Half-open equality at either boundary is allowed.
        !index
            .range(..=span.start)
            .next_back()
            .is_some_and(|(_, &end)| end > span.start)
            && !index
                .range(span.start..)
                .next()
                .is_some_and(|(&start, _)| start < span.end)
    })
}

/// Assign whole items to the first nonconflicting lane, retaining explicit
/// overflow. Adjacent intervals of different items do not conflict.
///
/// All count/work/storage budgets are checked before allocation or clone.
/// Duplicate-ID sorting uses bounded allocated scratch, after borrowed span
/// validation. Errors publish no partial layout. Zero lanes retain all overflow.
pub fn pack(items: &[Item], source_length: u64, limits: &Limits) -> Result<Layout, Error> {
    let work_bound = work_bound(items, limits)?;
    for item in items {
        if item.spans.is_empty() {
            return Err(Error::MissingIntervals(item.id));
        }
        for (index, span) in item.spans.iter().enumerate() {
            if span.start >= span.end || span.end > source_length {
                return Err(Error::InvalidInterval { id: item.id, index });
            }
        }
    }
    let mut order: Vec<usize> = (0..items.len()).collect();
    let mut order_scratch = vec![0; items.len()];
    merge_sort(&mut order, &mut order_scratch, |&a, &b| {
        items[a].id.cmp(&items[b].id)
    });
    for pair in order.windows(2) {
        if items[pair[0]].id == items[pair[1]].id {
            return Err(Error::DuplicateId(items[pair[0]].id));
        }
    }
    let largest = items.iter().map(|item| item.spans.len()).max().unwrap_or(0);
    let mut span_scratch = vec![Interval { start: 0, end: 0 }; largest];
    let mut assignments: Vec<_> = items
        .iter()
        .map(|item| Assignment {
            id: item.id,
            spans: normalize(&item.spans, &mut span_scratch),
            placement: Placement::Overflow,
        })
        .collect();
    merge_sort(&mut order, &mut order_scratch, |&a, &b| {
        (&assignments[a].spans, assignments[a].id).cmp(&(&assignments[b].spans, assignments[b].id))
    });
    let mut lanes: Vec<BTreeMap<u64, u64>> = Vec::with_capacity(limits.max_lanes.min(items.len()));
    for position in order {
        let lane = lanes
            .iter()
            .position(|lane| fits(lane, &assignments[position].spans));
        let lane = match lane {
            Some(lane) => lane,
            None if lanes.len() < limits.max_lanes => {
                lanes.push(BTreeMap::new());
                lanes.len() - 1
            }
            None => continue,
        };
        // Check ALL joined spans before inserting ANY, so a rejected item cannot
        // poison a lane with partial occupancy or affect later items.
        for span in &assignments[position].spans {
            lanes[lane].insert(span.start, span.end);
        }
        assignments[position].placement = Placement::Lane(lane);
    }
    Ok(Layout {
        assignments,
        lane_count: lanes.len(),
        work_bound,
    })
}
