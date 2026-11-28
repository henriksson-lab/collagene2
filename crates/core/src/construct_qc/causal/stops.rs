//! New stop comparisons use native ancestry and the canonical CDS translator.
use super::*;
use super::{
    coding::Coding,
    comparison::{CausePool, Cell, CellKind, OutputBudget},
    register::RegisterAnalysis,
};
use crate::construct_qc::Verdict;

pub(super) struct StopAnalysis {
    pub evidence: Vec<StopEvidence>,
    pub verdict: Verdict,
    pub reason: String,
}
fn new_call(cell: &Cell) -> bool {
    matches!(cell.kind, CellKind::Introduced | CellKind::Foreign)
}
pub(super) fn inspect(
    before: &Coding,
    after: &Coding,
    cells: &[Cell],
    pool: &CausePool,
    register: &RegisterAnalysis,
    native_annotation: bool,
    request: &Request,
    output: &mut OutputBudget,
    work: &mut Work,
) -> StopAnalysis {
    let mut evidence = Vec::new();
    let mut incomplete = cells.len() != after.coding_len()
        || after.translated.translation.possible_stop_codons > 0
        || after.translated.translation.trailing_bases > 0;
    let mut has_parent = native_annotation;
    for cell in cells {
        if !work.charge(1) {
            incomplete = true;
            break;
        }
        has_parent |= matches!(cell.kind, CellKind::Parent(_));
        incomplete |= matches!(cell.kind, CellKind::Unknown | CellKind::Outside);
    }
    incomplete |= !has_parent || register.verdict == Verdict::Indeterminate;
    let protein = after.translated.translation.protein.as_bytes();
    let mut failure = false;
    for (index, &amino) in protein.iter().enumerate() {
        if !work.charge(1) || output.exhausted {
            incomplete = true;
            break;
        }
        if amino != b'*' || index + 1 == protein.len() {
            continue;
        }
        let start = index * 3;
        let end = start + 3;
        let Some(group) = cells.get(start..end) else {
            incomplete = true;
            break;
        };
        let mut result = StopEvidence {
            verdict: Verdict::Indeterminate, product_codon_index: index,
            product_coordinates: after.codon_coordinates(index).expect("translated native codon"),
            product_codon: after.codon(index).expect("translated native codon").into(),
            parent_codon_index: None, parent_coordinates: None, parent_codon: None, parent_amino_acid: None,
            reason: "A definite product internal stop is present, but its creation relative to the declared parent is unresolved.".into(),
            causes: vec![],
        };
        let parents = group
            .iter()
            .filter_map(|cell| match cell.kind {
                CellKind::Parent(parent) => Some(parent),
                _ => None,
            })
            .collect::<Vec<_>>();
        let all_parent = parents.len() == 3;
        let consecutive = all_parent
            && parents
                .windows(2)
                .all(|pair| pair[0].ordinal.checked_add(1) == Some(pair[1].ordinal));
        let known = group
            .iter()
            .all(|cell| !matches!(cell.kind, CellKind::Unknown | CellKind::Outside));
        let changes = group.iter().any(|cell| cell.changed);
        let added = group.iter().any(new_call);
        let candidate = if consecutive {
            Some(parents[0].ordinal)
        } else if known && !parents.is_empty() {
            let mut candidate = None;
            let mut coherent = true;
            for (offset, cell) in group.iter().enumerate() {
                if let CellKind::Parent(parent) = cell.kind {
                    match parent.ordinal.checked_sub(offset) {
                        Some(value) if candidate.is_none_or(|old| old == value) => {
                            candidate = Some(value)
                        }
                        _ => coherent = false,
                    }
                }
            }
            if coherent { candidate } else { None }
        } else {
            None
        };
        let aligned_parent = candidate
            .filter(|ordinal| ordinal.is_multiple_of(3))
            .map(|ordinal| ordinal / 3)
            .filter(|codon| before.codon(*codon).is_some());
        if let Some(parent_index) = aligned_parent {
            result.parent_codon_index = Some(parent_index);
            result.parent_coordinates = before.codon_coordinates(parent_index);
            result.parent_codon = before.codon(parent_index).map(str::to_owned);
            result.parent_amino_acid = before.amino(parent_index);
            match result.parent_amino_acid {
                Some('*') if all_parent && consecutive && !changes => {
                    result.verdict = Verdict::Informational;
                    result.reason = "This exact stop codon is retained from the parent coding group, including a terminal stop moved into the selected interior. It is not asserted to be newly created.".into();
                }
                Some(amino)
                    if amino != '*'
                        && amino != 'X'
                        && known
                        && (changes || added)
                        && has_parent =>
                {
                    result.causes = pool.combine(group.iter().map(|cell| cell.causes));
                    if !result.causes.is_empty() {
                        result.verdict = Verdict::Fail;
                        result.reason = "Actual native changed/introduced calls convert the retained parent codon anchors from a definite sense codon into a definite internal stop, under the paired declarations/code.".into();
                    }
                }
                _ => {
                    result.reason = "The anchored parent codon is already a stop, uncertain, or lacks a complete native changed-call cause; new stop creation is not established.".into();
                }
            }
        } else if known && has_parent {
            let mut shifted_segment = None;
            if consecutive {
                for segment in &register.segments {
                    if !work.charge(1) {
                        incomplete = true;
                        break;
                    }
                    if segment.product_coding_start <= start
                        && segment.product_coding_end >= end
                        && segment.native_junction_proven
                        && segment.offset_mod3 != 0
                    {
                        shifted_segment = Some(segment);
                        break;
                    }
                }
            }
            let mut junction = None;
            if shifted_segment.is_none() {
                for anchor in &register.anchors {
                    if !work.charge(1) {
                        incomplete = true;
                        break;
                    }
                    if anchor.at > start && anchor.at < end
                        || added && anchor.affected_start < end && anchor.at >= start
                    {
                        junction = Some(anchor);
                        break;
                    }
                }
            }
            if let Some(segment) = shifted_segment {
                result.causes = segment.causes.clone();
                result.verdict = Verdict::Fail;
                result.reason = format!(
                    "A native DNA junction established the shifted register before this codon. Parent coding calls {}..{} now form a new definite stop group; distance from the junction does not erase its recorded cause.",
                    parents[0].ordinal,
                    parents[2].ordinal + 1
                );
            } else if let Some(anchor) = junction {
                result.causes = pool.combine(anchor.cause_sets.iter().copied());
                if !result.causes.is_empty() {
                    result.verdict = Verdict::Fail;
                    result.reason = "The definite stop group crosses or lies inside a native inserted/reference-joining block within the paired coding path; it is not a copied intact parent codon.".into();
                }
            }
        }
        if result.verdict == Verdict::Indeterminate {
            incomplete = true;
        }
        if output.stops >= request.limits.max_stops {
            output.exhausted = true;
            incomplete = true;
            break;
        }
        if !output.retain(&result) {
            incomplete = true;
            break;
        }
        failure |= result.verdict == Verdict::Fail;
        output.stops += 1;
        evidence.push(result);
    }
    let verdict = if failure {
        Verdict::Fail
    } else if incomplete || output.exhausted {
        Verdict::Indeterminate
    } else {
        Verdict::Pass
    };
    let reason = match verdict {
        Verdict::Fail => "At least one retained evidence item proves a new internal stop under the declared CDS/code and actual native changed-call/junction history. No expression or experimental outcome is inferred.",
        Verdict::Pass => "No new definite internal stop was established in the complete observed comparison; intact inherited stops remain explicitly distinguished. This does not certify overall CDS integrity.",
        _ => "Native ancestry, codon certainty, boundaries or work/evidence limits leave new-stop comparison incomplete; absence is not asserted.",
    }.into();
    StopAnalysis {
        evidence,
        verdict,
        reason,
    }
}
