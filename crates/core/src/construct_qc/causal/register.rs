//! Register changes require actual DNA adjacency/introduced-call evidence.
//! A different annotation start/order alone cannot establish a DNA-junction cause.
use super::*;
use super::{
    coding::Coding,
    comparison::{CausePool, Cell, CellKind, OutputBudget, ParentCall},
};
use crate::construct_qc::Verdict;
use std::collections::BTreeSet;

pub(super) struct Anchor {
    pub at: usize,
    pub affected_start: usize,
    pub phase: u8,
    pub cause_sets: Vec<usize>,
}
pub(super) struct RegisterAnalysis {
    pub segments: Vec<RegisterSegment>,
    pub anchors: Vec<Anchor>,
    pub verdict: Verdict,
    pub reason: String,
}
fn phase(index: usize, parent: ParentCall) -> u8 {
    ((index as i64 - parent.ordinal as i64).rem_euclid(3)) as u8
}
fn product_adjacent(product: &Coding, a: usize, b: usize) -> bool {
    product
        .path
        .call_at(product.offset + a)
        .and_then(|call| call.next(product.source_length, product.topology))
        == product.path.call_at(product.offset + b)
}
fn new_block(cell: &Cell) -> bool {
    matches!(cell.kind, CellKind::Introduced | CellKind::Foreign)
}
fn source_interval(before: &Coding, left: ParentCall, right: ParentCall, work: &mut Work) -> bool {
    if left.ordinal >= right.ordinal || left.call.reverse != right.call.reverse {
        return false;
    }
    let distance =
        |call: collagene_model::annotations::biological_path::BiologicalCall| -> Option<usize> {
            if call.reverse != left.call.reverse {
                return None;
            }
            let direct = if call.reverse {
                left.call.position.checked_sub(call.position)
            } else {
                call.position.checked_sub(left.call.position)
            };
            direct.or_else(|| {
                (before.topology == collagene_model::Topology::Circular).then(|| {
                    if call.reverse {
                        left.call.position + before.source_length - call.position
                    } else {
                        before.source_length - left.call.position + call.position
                    }
                })
            })
        };
    let Some(end) = distance(right.call).filter(|end| *end > 0) else {
        return false;
    };
    let mut previous = 0;
    for ordinal in left.ordinal + 1..=right.ordinal {
        if !work.charge(1) {
            return false;
        }
        let Some(call) = before.path.call_at(before.offset + ordinal) else {
            return false;
        };
        let Some(current) = distance(call) else {
            return false;
        };
        if current <= previous || current > end {
            return false;
        }
        previous = current;
    }
    previous == end
}
pub(super) fn inspect(
    before: &Coding,
    after: &Coding,
    cells: &[Cell],
    pool: &CausePool,
    request: &Request,
    output: &mut OutputBudget,
    work: &mut Work,
) -> RegisterAnalysis {
    let mut anchors = Vec::new();
    let mut unknown = cells.len() != after.coding_len() || cells.is_empty();
    for index in 1..cells.len() {
        if !work.charge(1) {
            unknown = true;
            break;
        }
        if let (CellKind::Parent(left), CellKind::Parent(right)) =
            (cells[index - 1].kind, cells[index].kind)
        {
            if product_adjacent(after, index - 1, index)
                && left.call.next(before.source_length, before.topology) != Some(right.call)
            {
                if anchors.len() >= request.limits.max_segments {
                    output.exhausted = true;
                    unknown = true;
                    break;
                }
                let cause_sets = vec![cells[index - 1].causes, cells[index].causes];
                if cause_sets.iter().any(|index| !pool.get(*index).is_empty()) {
                    anchors.push(Anchor {
                        at: index,
                        affected_start: index,
                        phase: phase(index, right),
                        cause_sets,
                    });
                } else {
                    unknown = true;
                }
            }
        }
    }
    let mut cursor = 0;
    while cursor < cells.len() {
        if !work.charge(1) {
            unknown = true;
            break;
        }
        if !new_block(&cells[cursor]) {
            cursor += 1;
            continue;
        }
        let start = cursor;
        while cursor < cells.len() && new_block(&cells[cursor]) {
            if !work.charge(1) {
                unknown = true;
                break;
            }
            cursor += 1;
        }
        if work.remaining() == 0 {
            break;
        }
        let end = cursor;
        let left = start
            .checked_sub(1)
            .and_then(|index| match cells[index].kind {
                CellKind::Parent(parent) => Some(parent),
                _ => None,
            });
        let right = cells.get(end).and_then(|cell| match cell.kind {
            CellKind::Parent(parent) => Some(parent),
            _ => None,
        });
        let first = start.saturating_sub(1);
        let last = if end < cells.len() {
            end
        } else {
            end.saturating_sub(1)
        };
        let mut contiguous = true;
        for index in first..last {
            if !work.charge(1) || !product_adjacent(after, index, index + 1) {
                contiguous = false;
                break;
            }
        }
        let supported = contiguous
            && match (left, right) {
                (Some(left), Some(right)) => source_interval(before, left, right, work),
                (None, Some(right)) => start == 0 && right.ordinal == 0,
                (Some(left), None) => {
                    end == after.coding_len()
                        && left.ordinal + 1 == before.coding_len()
                        && (end - start).is_multiple_of(3)
                        && before.coding_len().is_multiple_of(3)
                }
                _ => false,
            };
        if supported {
            let indices = (first..=last)
                .filter_map(|index| cells.get(index).map(|cell| cell.causes))
                .collect::<BTreeSet<_>>();
            if indices.len() > 64 || anchors.len() >= request.limits.max_segments {
                unknown = true;
                output.exhausted = true;
                break;
            }
            let cause_sets = indices.into_iter().collect::<Vec<_>>();
            if !cause_sets.iter().any(|index| !pool.get(*index).is_empty()) {
                unknown = true;
            } else {
                let offset = right.map(|parent| phase(end, parent)).unwrap_or(0);
                anchors.push(Anchor {
                    at: end,
                    affected_start: start,
                    phase: offset,
                    cause_sets,
                });
            }
        } else {
            unknown = true;
        }
    }
    anchors.sort_by_key(|anchor| anchor.at);
    let mut segments: Vec<RegisterSegment> = Vec::new();
    let mut active: Option<&Anchor> = None;
    let mut next_anchor = 0;
    let mut previous_parent: Option<usize> = None;
    let mut has_parent = false;
    let mut failure = false;
    for (index, cell) in cells.iter().enumerate() {
        if !work.charge(1) {
            unknown = true;
            break;
        }
        let mut starts_anchor = false;
        while next_anchor < anchors.len() && anchors[next_anchor].at <= index {
            active = Some(&anchors[next_anchor]);
            next_anchor += 1;
            starts_anchor = true;
        }
        let CellKind::Parent(parent) = cell.kind else {
            if matches!(cell.kind, CellKind::Unknown | CellKind::Outside) {
                unknown = true;
            }
            previous_parent = None;
            continue;
        };
        has_parent = true;
        let offset = phase(index, parent);
        let cause = active.filter(|anchor| anchor.phase == offset);
        let proven = cause.is_some();
        if offset != 0 && !proven {
            unknown = true;
        }
        let merge = !starts_anchor
            && segments.last().is_some_and(|segment| {
                segment.product_coding_end == index
                    && segment.offset_mod3 == offset
                    && segment.native_junction_proven == proven
                    && previous_parent.and_then(|previous| previous.checked_add(1))
                        == Some(parent.ordinal)
            });
        if merge {
            let last = segments.last_mut().expect("existing register segment");
            last.product_coding_end += 1;
            last.last_parent_coding_ordinal = parent.ordinal;
        } else {
            if output.segments >= request.limits.max_segments {
                output.exhausted = true;
                unknown = true;
                break;
            }
            let causes = cause.map_or_else(Vec::new, |anchor| {
                pool.combine(anchor.cause_sets.iter().copied())
            });
            if proven && causes.is_empty() {
                unknown = true;
            }
            let segment = RegisterSegment {
                product_coding_start: index,
                product_coding_end: index + 1,
                first_parent_coding_ordinal: parent.ordinal,
                last_parent_coding_ordinal: parent.ordinal,
                offset_mod3: offset,
                native_junction_proven: proven && !causes.is_empty(),
                causes,
            };
            if !output.retain(&segment) {
                unknown = true;
                break;
            }
            failure |= segment.offset_mod3 != 0 && segment.native_junction_proven;
            output.segments += 1;
            segments.push(segment);
        }
        previous_parent = Some(parent.ordinal);
    }
    let verdict = if failure {
        Verdict::Fail
    } else if unknown || !has_parent || output.exhausted {
        Verdict::Indeterminate
    } else {
        Verdict::Pass
    };
    let reason = match verdict {
        Verdict::Fail => "Actual native DNA junction/introduced-call anchors establish a changed coding register under the paired declarations. Later compensating junctions may restore phase; the segments retain both states.",
        Verdict::Pass => "All resolved paired-parent coding calls retain register; explicitly bounded inserted blocks are in frame. This does not certify CDS identity, expression or function.",
        _ => "DNA/role correspondence, native causality or bounded evidence is insufficient for a complete register conclusion. Annotation redefinition alone is not a DNA-junction cause.",
    }.into();
    RegisterAnalysis {
        segments,
        anchors,
        verdict,
        reason,
    }
}
