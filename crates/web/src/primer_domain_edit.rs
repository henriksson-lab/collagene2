//! Explicit edits to a selected reference interval, without binding searches.
//! Retained literal calls are never silently replaced by template DNA.
use crate::primer_pair_model::{Comparison, OligoDraft, SelectionAnchor, anchor_reference};
use collagene_core::{analysis::base_mask, primers::MAX_OLIGO, run_lineage::EventBinding};
use collagene_model::{Molecule, Topology};
use std::collections::BTreeMap;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Edit {
    /// Source-forward coordinates; orientation stays with the existing anchor.
    Endpoints { start: usize, end: usize },
    /// Length of the synthetic 5′ portion of the complete ordered oligo.
    TailBoundary { tail_length: usize },
}
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RemovedEdit {
    pub source_position: usize,
    pub old_domain_index: usize,
    pub reference: char,
    pub literal: char,
    pub ambiguous: bool,
    pub mismatch: bool,
}
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Preview {
    pub source: EventBinding,
    pub epoch: u64,
    pub before: OligoDraft,
    pub after: OligoDraft,
    pub edit: Edit,
    /// Literal differences leaving the anchored annealing domain. Tail-boundary
    /// edits retain these calls in the tail instead of deleting ordered DNA.
    pub removed_edits: Vec<RemovedEdit>,
    pub warnings: Vec<String>,
    pub comparison: Comparison,
}
fn geometry(m: &Molecule, start: usize, end: usize) -> Result<usize, String> {
    let length = end
        .checked_sub(start)
        .ok_or("Primer endpoints are reversed.")?;
    let n = m.sequence.len();
    if start >= n
        || length == 0
        || length > n
        || length > MAX_OLIGO
        || (m.topology == Topology::Linear && end > n)
    {
        return Err("Choose a nonempty interval of at most 200 bases and one source turn; circular starts must be normalized.".into());
    }
    Ok(length)
}
fn coordinate(anchor: &SelectionAnchor, index: usize, length: usize) -> usize {
    if anchor.reverse {
        (anchor.end - 1 - index) % length
    } else {
        (anchor.start + index) % length
    }
}
fn boundary_geometry(
    m: &Molecule,
    anchor: &SelectionAnchor,
    old_tail: usize,
    new_tail: usize,
    ordered: usize,
) -> Result<(usize, usize), String> {
    let domain = ordered
        .checked_sub(new_tail)
        .filter(|n| *n > 0)
        .ok_or("The annealing domain cannot be empty.")?;
    if domain > m.sequence.len() || domain > MAX_OLIGO {
        return Err("The annealing domain cannot exceed one source turn or 200 bases.".into());
    }
    let n = m.sequence.len();
    let (start, end) = if new_tail >= old_tail {
        let trim = new_tail - old_tail;
        if anchor.reverse {
            (
                anchor.start,
                anchor
                    .end
                    .checked_sub(trim)
                    .ok_or("The annealing domain cannot be empty.")?,
            )
        } else {
            let start = anchor
                .start
                .checked_add(trim)
                .ok_or("Primer endpoint overflow.")?;
            let start = if m.topology == Topology::Circular {
                start % n
            } else {
                start
            };
            (
                start,
                start
                    .checked_add(domain)
                    .ok_or("Primer endpoint overflow.")?,
            )
        }
    } else {
        let extension = old_tail - new_tail;
        if anchor.reverse {
            (
                anchor.start,
                anchor
                    .end
                    .checked_add(extension)
                    .ok_or("Primer endpoint overflow.")?,
            )
        } else {
            let start = if m.topology == Topology::Circular {
                anchor
                    .start
                    .checked_add(n)
                    .and_then(|p| p.checked_sub(extension))
                    .ok_or("Primer endpoint overflow.")?
                    % n
            } else {
                anchor
                    .start
                    .checked_sub(extension)
                    .ok_or("Promoted tail bases extend beyond the linear source's 5′ context.")?
            };
            (
                start,
                start
                    .checked_add(domain)
                    .ok_or("Primer endpoint overflow.")?,
            )
        }
    };
    geometry(m, start, end)?;
    Ok((start, end))
}

/// Build a reviewable, source-bound proposal. Incomplete/invalid editable text
/// stays in the caller's draft; it is never repaired or silently normalized here.
pub fn preview(
    m: &Molecule,
    source: &EventBinding,
    epoch: u64,
    draft: &OligoDraft,
    edit: Edit,
) -> Result<Preview, String> {
    let old = draft.anchor.as_ref().ok_or(
        "Choose an exact reference anchor before editing primer endpoints or tail boundaries.",
    )?;
    if &old.source != source {
        return Err(
            "The primer anchor belongs to another source event; select the intended region again."
                .into(),
        );
    }
    let old_length = geometry(m, old.start, old.end)?;
    if old.reference.len() != old_length {
        return Err("The retained anchor reference length is invalid.".into());
    }
    // primer() bounds raw input before constructing its canonical DNA calls.
    let primer = draft.primer()?;
    if primer.annealing.len() != old_length {
        return Err("The annealing domain length no longer matches its selected reference; select the intended region again.".into());
    }
    let ordered_length = primer.tail.len() + primer.annealing.len();
    let (start, end) = match edit {
        Edit::Endpoints { start, end } => {
            let length = geometry(m, start, end)?;
            if length
                .checked_add(primer.tail.len())
                .is_none_or(|n| n > MAX_OLIGO)
            {
                return Err("The complete ordered oligo would exceed 200 bases.".into());
            }
            (start, end)
        }
        Edit::TailBoundary { tail_length } => {
            boundary_geometry(m, old, primer.tail.len(), tail_length, ordered_length)?
        }
    };
    // The shared helper guards the entire borrowed source before hashing or
    // extracting it, then verifies the actual molecular revision and geometry.
    if anchor_reference(m, source, old.start, old.end, old.reverse)? != old.reference {
        return Err("The retained anchor differs from the actual source reference.".into());
    }
    let reference = anchor_reference(m, source, start, end, old.reverse)?;
    let anchor = SelectionAnchor {
        source: source.clone(),
        start,
        end,
        reverse: old.reverse,
        reference,
    };
    let n = m.sequence.len();
    let old_calls = (0..old_length)
        .map(|i| (coordinate(old, i, n), i))
        .collect::<BTreeMap<_, _>>();
    let new_calls = (0..end - start)
        .map(|i| (coordinate(&anchor, i, n), i))
        .collect::<BTreeMap<_, _>>();
    let (annealing, tail) = match edit {
        Edit::Endpoints { .. } => {
            let mut domain = anchor.reference.as_bytes().to_vec();
            for (position, index) in &new_calls {
                if let Some(old_index) = old_calls.get(position) {
                    domain[*index] = primer.annealing.as_str().as_bytes()[*old_index];
                }
            }
            (
                String::from_utf8(domain).map_err(|_| "Primer calls must be canonical DNA.")?,
                draft.tail.clone(),
            )
        }
        Edit::TailBoundary { tail_length } => {
            let ordered = format!("{}{}", primer.tail, primer.annealing.as_str());
            (
                ordered[tail_length..].to_owned(),
                ordered[..tail_length].to_owned(),
            )
        }
    };
    let after = OligoDraft {
        name: draft.name.clone(),
        annealing,
        tail,
        anchor: Some(anchor),
    };
    let mut removed_edits = Vec::new();
    for i in 0..old_length {
        let position = coordinate(old, i, n);
        let literal = primer.annealing.as_str().as_bytes()[i];
        let reference = old.reference.as_bytes()[i];
        if !new_calls.contains_key(&position) && literal != reference {
            let a = base_mask(reference);
            let b = base_mask(literal);
            removed_edits.push(RemovedEdit {
                source_position: position,
                old_domain_index: i,
                reference: reference as char,
                literal: literal as char,
                ambiguous: a.count_ones() != 1 || b.count_ones() != 1,
                mismatch: a & b == 0,
            });
        }
    }
    let comparison = after
        .comparison()?
        .ok_or("The proposed primer has no reference comparison.")?;
    let mut warnings=vec!["This is an explicit selected-reference comparison, not a binding search or evidence of PCR amplification.".into()];
    if !removed_edits.is_empty() {
        warnings.push(match edit {Edit::Endpoints {..}=>format!("{} literal changed calls leave the ordered oligo because their source positions were removed.",removed_edits.len()),Edit::TailBoundary {..}=>format!("{} literal changed calls leave the annealing domain and remain in the synthetic tail; the complete ordered oligo is unchanged.",removed_edits.len())});
    }
    if let Edit::TailBoundary { tail_length } = edit {
        warnings.push("Moving the tail boundary preserves every ordered DNA call; added tail DNA is not source-derived binding evidence.".into());
        if tail_length < primer.tail.len() {
            warnings.push("Promoted tail calls remain literal DNA. Their newly selected source context can mismatch or contain ambiguous calls.".into());
        }
        let updated = after.primer()?;
        if format!("{}{}", updated.tail, updated.annealing.as_str())
            != format!("{}{}", primer.tail, primer.annealing.as_str())
        {
            return Err("Tail boundary changed the complete ordered oligo.".into());
        }
    }
    Ok(Preview {
        source: source.clone(),
        epoch,
        before: draft.clone(),
        after,
        edit,
        removed_edits,
        warnings,
        comparison,
    })
}

/// Recompute the whole proposal from the current immutable source and draft;
/// preview data and names never authorize edited letters or an inherited anchor.
pub fn accept(
    proposed: &Preview,
    m: &Molecule,
    source: &EventBinding,
    epoch: u64,
    current: &OligoDraft,
) -> Result<OligoDraft, String> {
    if &proposed.source != source || proposed.epoch != epoch || &proposed.before != current {
        return Err(
            "The source event, workspace or primer draft changed; preview the edit again.".into(),
        );
    }
    let recomputed = preview(m, source, epoch, current, proposed.edit)?;
    if &recomputed != proposed {
        return Err(
            "The primer edit preview differs from the freshly reconstructed proposal.".into(),
        );
    }
    Ok(recomputed.after)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::primer_pair_model::MAX_RAW_BASE_BYTES;
    fn source(sequence: &str, topology: Topology) -> (Molecule, EventBinding) {
        let m = Molecule::dna("literal endpoint fixture", sequence, topology).unwrap();
        let b = EventBinding {
            document_id: "literal-source".into(),
            event_id: format!("event-{}", "a".repeat(64)),
            revision: m.revision(),
        };
        (m, b)
    }
    fn draft(
        m: &Molecule,
        b: &EventBinding,
        start: usize,
        end: usize,
        reverse: bool,
        literal: &str,
        tail: &str,
    ) -> OligoDraft {
        OligoDraft {
            name: "Deliberate primer".into(),
            annealing: literal.into(),
            tail: tail.into(),
            anchor: Some(SelectionAnchor {
                source: b.clone(),
                start,
                end,
                reverse,
                reference: anchor_reference(m, b, start, end, reverse).unwrap(),
            }),
        }
    }
    fn ordered(draft: &OligoDraft) -> String {
        let p = draft.primer().unwrap();
        format!("{}{}", p.tail, p.annealing.as_str())
    }
    #[test]
    fn forward_extension_preserves_compatible_ambiguity_and_mismatched_literal_calls() {
        let (m, b) = source("ACGTTGCAAGTC", Topology::Linear);
        let before = draft(&m, &b, 2, 8, false, "GYTGTA", "NN");
        let p = preview(&m, &b, 4, &before, Edit::Endpoints { start: 1, end: 10 }).unwrap();
        assert_eq!(p.after.annealing, "CGYTGTAAG");
        assert_eq!(p.after.tail, "NN");
        assert_eq!(p.after.name, before.name);
        assert_eq!(p.comparison.mismatches, vec![5]);
        assert_eq!(p.comparison.ambiguous, vec![2]);
        assert!(p.removed_edits.is_empty());
        assert_eq!(accept(&p, &m, &b, 4, &before).unwrap(), p.after);
        assert_eq!(before.annealing, "GYTGTA");
        assert_eq!(m.sequence.as_str(), "ACGTTGCAAGTC");
    }
    #[test]
    fn three_prime_trim_reports_removed_changed_calls_in_both_orientations() {
        let (m, b) = source("ACGTTGCAAGTC", Topology::Linear);
        let before = draft(&m, &b, 2, 8, false, "GYTGTA", "");
        let p = preview(&m, &b, 0, &before, Edit::Endpoints { start: 2, end: 6 }).unwrap();
        assert_eq!(p.after.annealing, "GYTG");
        assert_eq!(
            p.removed_edits,
            vec![RemovedEdit {
                source_position: 6,
                old_domain_index: 4,
                reference: 'C',
                literal: 'T',
                ambiguous: false,
                mismatch: true
            }]
        );
        let before = draft(&m, &b, 2, 8, true, "TRCAAT", "GG");
        let extension = preview(&m, &b, 0, &before, Edit::Endpoints { start: 1, end: 9 }).unwrap();
        assert_eq!(extension.after.annealing, "TTRCAATG");
        assert_eq!(
            extension.after.anchor.as_ref().unwrap().reference,
            "TTGCAACG"
        );
        let p = preview(&m, &b, 0, &before, Edit::Endpoints { start: 4, end: 8 }).unwrap();
        assert_eq!(p.after.annealing, "TRCA");
        assert_eq!(p.removed_edits[0].source_position, 2);
        assert_eq!(p.removed_edits[0].old_domain_index, 5);
        assert_eq!(p.removed_edits[0].literal, 'T');
        assert!(p.removed_edits[0].mismatch);
    }
    #[test]
    fn circular_overlap_is_by_actual_coordinate_across_normalized_origin() {
        let (m, b) = source("ACGTTGCA", Topology::Circular);
        let before = draft(&m, &b, 6, 10, false, "CARY", "AT");
        let p = preview(&m, &b, 0, &before, Edit::Endpoints { start: 7, end: 12 }).unwrap();
        assert_eq!(p.after.annealing, "ARYGT");
        assert_eq!(p.after.anchor.as_ref().unwrap().reference, "AACGT");
        let before = draft(&m, &b, 6, 10, true, "RYTG", "");
        let p = preview(&m, &b, 0, &before, Edit::Endpoints { start: 5, end: 10 }).unwrap();
        assert_eq!(p.after.annealing, "RYTGC");
        assert_eq!(p.after.anchor.as_ref().unwrap().reference, "GTTGC");
    }
    #[test]
    fn moving_boundary_preserves_ordered_dna_and_reports_calls_moved_into_tail() {
        let (m, b) = source("ACGTTGCAAGTC", Topology::Linear);
        let before = draft(&m, &b, 2, 8, false, "GYTGTA", "NN");
        let p = preview(&m, &b, 0, &before, Edit::TailBoundary { tail_length: 4 }).unwrap();
        assert_eq!(p.after.tail, "NNGY");
        assert_eq!(p.after.annealing, "TGTA");
        assert_eq!(ordered(&p.after), ordered(&before));
        assert_eq!(
            (
                p.after.anchor.as_ref().unwrap().start,
                p.after.anchor.as_ref().unwrap().end
            ),
            (4, 8)
        );
        assert_eq!(p.removed_edits.len(), 1);
        assert_eq!(p.removed_edits[0].literal, 'Y');
        assert!(p.removed_edits[0].ambiguous);
        assert!(!p.removed_edits[0].mismatch);
        let before = draft(&m, &b, 2, 8, true, "TRCAAT", "NN");
        let p = preview(&m, &b, 0, &before, Edit::TailBoundary { tail_length: 3 }).unwrap();
        assert_eq!(p.after.tail, "NNT");
        assert_eq!(p.after.annealing, "RCAAT");
        assert_eq!(ordered(&p.after), ordered(&before));
        assert_eq!(p.after.anchor.as_ref().unwrap().end, 7);
    }
    #[test]
    fn promoted_tail_calls_remain_literal_with_fresh_forward_and_reverse_circular_context() {
        let (m, b) = source("ACGTTGCAAGTC", Topology::Linear);
        let before = draft(&m, &b, 4, 8, false, "TRCA", "AA");
        let p = preview(&m, &b, 0, &before, Edit::TailBoundary { tail_length: 0 }).unwrap();
        assert_eq!(p.after.annealing, "AATRCA");
        assert_eq!(p.after.anchor.as_ref().unwrap().reference, "GTTGCA");
        assert_eq!(p.comparison.mismatches, vec![0, 1]);
        assert_eq!(p.comparison.ambiguous, vec![3]);
        assert_eq!(ordered(&p.after), ordered(&before));
        let (m, b) = source("ACGTTGCA", Topology::Circular);
        let before = draft(&m, &b, 6, 9, true, "TYG", "GC");
        let p = preview(&m, &b, 0, &before, Edit::TailBoundary { tail_length: 1 }).unwrap();
        assert_eq!(p.after.tail, "G");
        assert_eq!(p.after.annealing, "CTYG");
        assert_eq!(p.after.anchor.as_ref().unwrap().reference, "GTTG");
        assert_eq!(p.comparison.mismatches, vec![0]);
        assert_eq!(p.comparison.ambiguous, vec![2]);
        assert_eq!(ordered(&p.after), ordered(&before));
    }
    #[test]
    fn acceptance_recomputes_and_refuses_stale_or_tampered_previews() {
        let (m, b) = source("ACGTTGCAAGTC", Topology::Linear);
        let before = draft(&m, &b, 2, 8, false, "GYTGTA", "NN");
        let p = preview(&m, &b, 3, &before, Edit::Endpoints { start: 2, end: 6 }).unwrap();
        assert!(accept(&p, &m, &b, 4, &before).is_err());
        let mut current = before.clone();
        current.tail = "GG".into();
        assert!(accept(&p, &m, &b, 3, &current).is_err());
        let mut other = b.clone();
        other.event_id = format!("event-{}", "b".repeat(64));
        assert!(accept(&p, &m, &other, 3, &before).is_err());
        let mut changed = m.clone();
        changed.name = "another actual revision".into();
        assert!(accept(&p, &changed, &b, 3, &before).is_err());
        let mut tampered = p.clone();
        tampered.after.annealing = "AAAA".into();
        assert!(accept(&tampered, &m, &b, 3, &before).is_err());
        let mut tampered = p.clone();
        tampered.removed_edits.clear();
        assert!(accept(&tampered, &m, &b, 3, &before).is_err());
        let mut tampered = p.clone();
        tampered.edit = Edit::Endpoints { start: 1, end: 6 };
        assert!(accept(&tampered, &m, &b, 3, &before).is_err());
    }
    #[test]
    fn missing_stale_length_changed_and_invalid_anchors_remain_unresolved() {
        let (m, b) = source("ACGTTGCA", Topology::Circular);
        let before = draft(&m, &b, 6, 10, false, "CAAC", "");
        let edit = Edit::Endpoints { start: 0, end: 3 };
        let mut invalid = before.clone();
        invalid.anchor = None;
        assert!(preview(&m, &b, 0, &invalid, edit).is_err());
        let mut invalid = before.clone();
        invalid.anchor.as_mut().unwrap().reference = "NNNN".into();
        assert!(preview(&m, &b, 0, &invalid, edit).is_err());
        let mut invalid = before.clone();
        invalid.annealing.push('A');
        assert!(preview(&m, &b, 0, &invalid, edit).is_err());
        let mut invalid = before.clone();
        invalid.annealing = "?".into();
        assert!(preview(&m, &b, 0, &invalid, edit).is_err());
        for edit in [
            Edit::Endpoints { start: 8, end: 10 },
            Edit::Endpoints { start: 6, end: 15 },
            Edit::Endpoints { start: 0, end: 0 },
            Edit::Endpoints {
                start: usize::MAX,
                end: 0,
            },
            Edit::TailBoundary { tail_length: 4 },
            Edit::TailBoundary {
                tail_length: usize::MAX,
            },
        ] {
            assert!(preview(&m, &b, 0, &before, edit).is_err());
        }
        let mut invalid = before.clone();
        invalid.annealing = "A".repeat(MAX_RAW_BASE_BYTES + 1);
        assert!(preview(&m, &b, 0, &invalid, edit).is_err());
    }
    #[test]
    fn complete_oligo_and_linear_or_single_turn_context_limits_are_enforced() {
        let (m, b) = source(&"A".repeat(300), Topology::Linear);
        let before = draft(&m, &b, 0, 190, false, &"A".repeat(190), &"G".repeat(10));
        assert!(preview(&m, &b, 0, &before, Edit::Endpoints { start: 0, end: 191 }).is_err());
        assert!(preview(&m, &b, 0, &before, Edit::TailBoundary { tail_length: 9 }).is_err());
        let (m, b) = source("ACGT", Topology::Circular);
        let before = draft(&m, &b, 3, 7, false, "TACG", "NN");
        assert!(preview(&m, &b, 0, &before, Edit::TailBoundary { tail_length: 1 }).is_err());
        let (m, b) = source("ACGT", Topology::Linear);
        let before = draft(&m, &b, 1, 3, true, "CG", "NN");
        assert!(preview(&m, &b, 0, &before, Edit::TailBoundary { tail_length: 0 }).is_err());
    }
}
