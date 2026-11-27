//! Literal internal contribution intent, independent of primers and preparation.
//!
//! UI entry uses one-based half-open boundaries: `3..5` replaces contribution
//! bases3 and4, `3..3` inserts before base3, and `length+1` appends at the end.
//! Stored offsets are zero-based in the ORIGINAL contribution's construct
//! orientation, before junction trims. Source documents are never modified.
use crate::assembly_draft_model::Draft;
use collagene_core::{Alphabet, Sequence, Topology, assembly_edit_intent::Edit};

#[derive(Clone)]
pub struct SourceReview {
    contributions: Vec<(
        crate::source_identity::SourceInfo,
        Option<crate::run_origin_model::PinnedOrigin>,
        collagene_core::assembly_intent::PlacedFragment,
    )>,
}
impl PartialEq for SourceReview {
    fn eq(&self, other: &Self) -> bool {
        self.contributions.len() == other.contributions.len()
            && self.contributions.iter().zip(&other.contributions).all(
                |((a, origin, fragment), (b, other_origin, other_fragment))| {
                    a == b
                        && origin == other_origin
                        && fragment.selector == other_fragment.selector
                        && fragment.source_event == other_fragment.source_event
                        && fragment.source_state_revision == other_fragment.source_state_revision
                        && fragment.original_sequence == other_fragment.original_sequence
                        && fragment.top_mappings == other_fragment.top_mappings
                        && fragment.source_domain == other_fragment.source_domain
                        && fragment.incomplete == other_fragment.incomplete
                },
            )
    }
}
impl Eq for SourceReview {}
impl SourceReview {
    pub fn capture(
        draft: &Draft,
        fragments: &[collagene_core::assembly_intent::PlacedFragment],
    ) -> Result<Self, String> {
        draft.validate()?;
        if fragments.len() != draft.slots.len() {
            return Err("Source contribution count changed".into());
        }
        let contributions = draft
            .slots
            .iter()
            .zip(&draft.methods.rows)
            .zip(fragments)
            .map(|((slot, row), fragment)| {
                let capture = slot.captured.as_ref().ok_or("Pin the source first")?;
                if fragment.input_revision != capture.source.binding.revision
                    || fragment.selector != slot.selector
                    || fragment.original_sequence.len()
                        != slot.selector.length().map_err(|e| e.to_string())?
                    || fragment.original_sequence.alphabet() != Alphabet::Dna
                    || fragment.top_mappings.iter().any(|mapping| {
                        mapping.source_end <= mapping.source_start
                            || mapping.source_end > capture.source.length
                            || mapping
                                .destination_start
                                .checked_add(mapping.source_end - mapping.source_start)
                                .is_none_or(|end| end > fragment.original_sequence.len())
                    })
                {
                    return Err(
                        "Resolved original contribution differs from the pinned source".into(),
                    );
                }
                if matches!(
                    &slot.selector,
                    collagene_core::assembly_intent::Selector::MaterialStrandRange { .. }
                ) && fragment.source_event.as_ref() != Some(&capture.source.binding)
                {
                    return Err("Actual strand belongs to another exact source event".into());
                }
                Ok((
                    capture.source.clone(),
                    row.run_origin.clone(),
                    fragment.clone(),
                ))
            })
            .collect::<Result<Vec<_>, String>>()?;
        Ok(Self { contributions })
    }
    pub fn get<'a>(
        &'a self,
        draft: &Draft,
        index: usize,
    ) -> Result<&'a collagene_core::assembly_intent::PlacedFragment, String> {
        let (source, origin, fragment) = self
            .contributions
            .get(index)
            .ok_or("Review the original source contribution first")?;
        let slot = draft
            .slots
            .get(index)
            .ok_or("Source contribution is unavailable")?;
        let captured = slot.captured.as_ref().ok_or("Pin the source first")?;
        if captured.source != *source
            || slot.selector != fragment.selector
            || draft
                .methods
                .rows
                .get(index)
                .ok_or("Source preparation row is unavailable")?
                .run_origin
                != *origin
        {
            return Err(
                "Source, strand or material origin changed; review the intended construct again"
                    .into(),
            );
        }
        Ok(fragment)
    }
}
pub fn top_regions(
    draft: &Draft,
    slot: usize,
    start: usize,
    end: usize,
    review: Option<&SourceReview>,
) -> Result<Vec<(usize, usize, bool)>, String> {
    let length = draft
        .slots
        .get(slot)
        .ok_or("Source contribution is unavailable")?
        .selector
        .length()
        .map_err(|e| e.to_string())?;
    if start > end || end > length {
        return Err("Navigation lies outside the original contribution".into());
    }
    if let collagene_core::assembly_intent::Selector::ReferenceRange { reverse, .. } =
        &draft.slots[slot].selector
    {
        let (a, b) = source_interval(draft, slot, start, end)?;
        return Ok(vec![(a, b, *reverse)]);
    }
    let fragment = review
        .ok_or("Review the original actual strand before navigating it")?
        .get(draft, slot)?;
    let mut ranges = Vec::new();
    for mapping in &fragment.top_mappings {
        let length = mapping
            .source_end
            .checked_sub(mapping.source_start)
            .ok_or("Invalid source mapping")?;
        let stop = mapping
            .destination_start
            .checked_add(length)
            .ok_or("Source mapping overflow")?;
        let a = start.max(mapping.destination_start);
        let b = end.min(stop);
        if a >= b {
            continue;
        }
        let (lo, hi) = if mapping.reverse {
            (
                mapping.source_end - (b - mapping.destination_start),
                mapping.source_end - (a - mapping.destination_start),
            )
        } else {
            (
                mapping.source_start + a - mapping.destination_start,
                mapping.source_start + b - mapping.destination_start,
            )
        };
        ranges.push((lo, hi, mapping.reverse));
    }
    Ok(ranges)
}

fn contribution(
    draft: &Draft,
    slot: usize,
) -> Result<(&crate::assembly_draft_model::Slot, usize), String> {
    let s = draft
        .slots
        .get(slot)
        .ok_or("The source contribution is unavailable.")?;
    let capture = s
        .captured
        .as_ref()
        .ok_or("Pin an original source region first.")?;
    let n = capture.molecule.sequence.len();
    let (source_start, source_end, _) =
        s.selector.reference_interval().map_err(|e| e.to_string())?;
    let length = source_end
        .checked_sub(source_start)
        .filter(|n| *n > 0)
        .ok_or("The original contribution interval is invalid.")?;
    if n == 0
        || source_start >= n
        || length > n
        || capture.molecule.sequence.alphabet() != Alphabet::Dna
        || (capture.molecule.topology == Topology::Linear && source_end > n)
    {
        return Err("Choose a valid original DNA contribution of at most one source turn.".into());
    }
    Ok((s, length))
}

/// Exact source interval in construct orientation; never derives letters from
/// already edited targets, proposed primers or reconstructed products.
pub fn source_text(
    draft: &Draft,
    slot: usize,
    start: usize,
    end: usize,
    review: Option<&SourceReview>,
) -> Result<String, String> {
    if matches!(
        &draft
            .slots
            .get(slot)
            .ok_or("Source contribution is unavailable")?
            .selector,
        collagene_core::assembly_intent::Selector::MaterialStrandRange { .. }
    ) {
        let fragment = review
            .ok_or("Review the intended construct once to read its original actual strand")?
            .get(draft, slot)?;
        if start > end || end > fragment.original_sequence.len() || end - start > 100_000 {
            return Err("Choose at most100000 original contribution bases".into());
        }
        return Ok(fragment.original_sequence.as_str()[start..end].to_owned());
    }
    let (s, length) = contribution(draft, slot)?;
    if start > end || end > length || end - start > 100_000 {
        return Err("Choose at most100000 bases within the original contribution.".into());
    }
    let capture = s.captured.as_ref().unwrap();
    let (source_start, _, reverse) = s.selector.reference_interval().map_err(|e| e.to_string())?;
    let n = capture.molecule.sequence.len();
    let (a, b) = if reverse {
        (length - end, length - start)
    } else {
        (start, end)
    };
    let bases = capture.molecule.sequence.as_str().as_bytes();
    let text: String = (a..b)
        .map(|offset| bases[(source_start + offset) % n] as char)
        .collect();
    let text = if reverse && !text.is_empty() {
        Sequence::dna(&text)
            .and_then(|s| s.reverse_complement())
            .map_err(|e| e.to_string())?
            .as_str()
            .to_owned()
    } else {
        text
    };
    if !text.bytes().all(|b| b"ACGTRYSWKMBDHVN".contains(&b)) {
        return Err("The original contribution contains unsupported DNA symbols.".into());
    }
    Ok(text)
}

/// Map an exact source range/caret (zero-based; circular end may be unrolled)
/// into a contiguous original contribution interval, in construct orientation.
/// The caller must separately check the selected document's pinned event identity.
pub fn selection_interval(
    draft: &Draft,
    slot: usize,
    start: usize,
    end: usize,
    review: Option<&SourceReview>,
) -> Result<(usize, usize), String> {
    if matches!(
        &draft
            .slots
            .get(slot)
            .ok_or("Source contribution is unavailable")?
            .selector,
        collagene_core::assembly_intent::Selector::MaterialStrandRange { .. }
    ) {
        let fragment = review
            .ok_or("Review the actual strand before copying source selections")?
            .get(draft, slot)?;
        let mut selected = None;
        for mapping in &fragment.top_mappings {
            if start < mapping.source_start || end > mapping.source_end || start > end {
                continue;
            }
            let interval = if mapping.reverse {
                (
                    mapping.destination_start + mapping.source_end - end,
                    mapping.destination_start + mapping.source_end - start,
                )
            } else {
                (
                    mapping.destination_start + start - mapping.source_start,
                    mapping.destination_start + end - mapping.source_start,
                )
            };
            if selected.replace(interval).is_some() {
                return Err("Source selection maps to multiple actual strand locations; enter contribution positions explicitly".into());
            }
        }
        return selected.ok_or("Selection includes bases outside the displayed top mapping; enter actual contribution positions explicitly".into());
    }
    let (s, length) = contribution(draft, slot)?;
    let capture = s.captured.as_ref().unwrap();
    let (source_start, _, reverse) = s.selector.reference_interval().map_err(|e| e.to_string())?;
    let n = capture.molecule.sequence.len();
    let count = end
        .checked_sub(start)
        .filter(|v| *v <= n)
        .ok_or("Select at most one original source turn.")?;
    if start > n
        || (start == n && count > 0)
        || (capture.molecule.topology == Topology::Linear && end > n)
    {
        return Err("Selection is outside the original source.".into());
    }
    let mut unrolled = start;
    if capture.molecule.topology == Topology::Circular && unrolled < source_start {
        unrolled = unrolled
            .checked_add(n)
            .ok_or("Selection coordinates overflow.")?;
    }
    let offset = unrolled
        .checked_sub(source_start)
        .ok_or("Select inside the original contribution.")?;
    let limit = offset
        .checked_add(count)
        .filter(|v| *v <= length)
        .ok_or("Selection crosses the chosen contribution boundary.")?;
    let interval = if reverse {
        (length - limit, length - offset)
    } else {
        (offset, limit)
    };
    source_text(draft, slot, interval.0, interval.1, review)?;
    Ok(interval)
}

/// Map stored construct-oriented boundaries back to the original source. The
/// returned source end stays unrolled across a circular origin for exact Go to.
pub fn source_interval(
    draft: &Draft,
    slot: usize,
    start: usize,
    end: usize,
) -> Result<(usize, usize), String> {
    let (s, length) = contribution(draft, slot)?;
    if start > end || end > length {
        return Err("Navigation lies outside the original contribution".into());
    }
    let capture = s.captured.as_ref().unwrap();
    let (source_start, _, reverse) = s.selector.reference_interval().map_err(|e| e.to_string())?;
    let n = capture.molecule.sequence.len();
    let offset = if reverse { length - end } else { start };
    let raw = source_start
        .checked_add(offset)
        .ok_or("Source coordinates overflow.")?;
    let source_start = if capture.molecule.topology == Topology::Circular {
        raw % n
    } else {
        raw
    };
    let source_end = source_start
        .checked_add(end - start)
        .ok_or("Source coordinates overflow.")?;
    Ok((source_start, source_end))
}

fn boundary(text: &str) -> Result<usize, String> {
    if text.is_empty() || text.len() > 20 || !text.bytes().all(|b| b.is_ascii_digit()) {
        return Err("Enter one-based contribution boundaries.".into());
    }
    text.parse::<usize>()
        .ok()
        .and_then(|v| v.checked_sub(1))
        .ok_or_else(|| "Contribution boundaries start at1.".into())
}

pub fn declare(
    draft: &Draft,
    slot: usize,
    start: &str,
    end: &str,
    replacement: &str,
    review: Option<&SourceReview>,
) -> Result<Draft, String> {
    crate::assembly_alternatives::require_active_choice_match(draft)?;
    let start = boundary(start)?;
    let end = boundary(end)?;
    let expected = source_text(draft, slot, start, end, review)?;
    if replacement.len() > 200_000 || !replacement.is_ascii() {
        return Err("Provide bounded DNA with optional ASCII whitespace.".into());
    }
    let replacement: String = replacement
        .bytes()
        .filter(|b| !b.is_ascii_whitespace())
        .map(|b| (b as char).to_ascii_uppercase())
        .collect();
    if replacement.len() > 100_000 || !replacement.bytes().all(|b| b"ACGTRYSWKMBDHVN".contains(&b))
    {
        return Err(
            "Use at most100000 IUPAC replacement DNA bases, or leave empty to delete.".into(),
        );
    }
    if replacement == expected {
        return Err("Declare a change to the original contribution.".into());
    }
    let roles = draft.methods.role_ids(draft.slots.len())?;
    let role = roles
        .get(slot)
        .ok_or("The source role is unavailable.")?
        .role
        .clone();
    let mut next = draft.clone();
    next.methods.source_roles = Some(roles);
    next.edits.push(Edit {
        role,
        start,
        end,
        expected,
        replacement,
    });
    next.sync_rows();
    next.validate()?;
    Ok(next)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{assembly_draft_model::CapturedSource, source_identity::SourceInfo};
    use collagene_core::{Molecule, run_lineage::EventBinding, workspace::Workspace};
    use std::sync::Arc;
    fn fixture(topology: Topology) -> Draft {
        let molecule = Molecule::dna("Original", "AACCGGTT", topology).unwrap();
        let workspace = Workspace::from_documents("Original", vec![molecule.clone()], 0).unwrap();
        let doc = &workspace.documents[0];
        let capture = Arc::new(CapturedSource {
            source: SourceInfo {
                binding: EventBinding {
                    document_id: doc.id.clone(),
                    event_id: doc.head_event.clone().unwrap(),
                    revision: doc.head.clone(),
                },
                name: molecule.name.clone(),
                length: 8,
                alphabet: Alphabet::Dna,
                topology,
            },
            molecule,
        });
        let mut draft = Draft::default();
        for slot in &mut draft.slots {
            slot.captured = Some(capture.clone());
            {
                let collagene_core::assembly_intent::Selector::ReferenceRange {
                    start: value, ..
                } = &mut slot.selector
                else {
                    panic!("fixture requires reference selector")
                };
                *value = 0;
            }
            {
                let collagene_core::assembly_intent::Selector::ReferenceRange {
                    end: value, ..
                } = &mut slot.selector
                else {
                    panic!("fixture requires reference selector")
                };
                *value = 8;
            }
        }
        draft.sync_rows();
        draft.validate().unwrap();
        draft
    }
    #[test]
    fn declarations_preserve_originals_and_reference_boundaries() {
        let original = fixture(Topology::Linear);
        let mut d = declare(&original, 0, "3", "5", "tag", None).unwrap();
        assert_eq!(d.version, crate::assembly_draft_model::CURRENT_VERSION);
        assert_eq!(d.edits[0].expected, "CC");
        assert_eq!(d.edits[0].replacement, "TAG");
        assert!(d.slots == original.slots && d.methods.rows == original.methods.rows);
        d = declare(&d, 0, "9", "9", "AC", None).unwrap();
        assert_eq!((d.edits[1].start, d.edits[1].end), (8, 8));
        assert_eq!(d.edits[1].expected, "");
        assert!(declare(&d, 0, "4", "4", "G", None).is_err());
        assert!(declare(&d, 0, "5", "5", "G", None).is_err());
        assert!(declare(&d, 0, "3", "5", "", None).is_err());
        let deletion = declare(&original, 0, "2", "4", "", None).unwrap();
        assert_eq!(deletion.edits[0].expected, "AC");
        assert_eq!(deletion.edits[0].replacement, "");
        assert!(declare(&original, 0, "0", "1", "T", None).is_err());
        assert!(declare(&original, 0, "2", "4", "AC", None).is_err());
    }
    #[test]
    fn precise_circular_reverse_ranges_and_carets_map_to_original_contribution() {
        let mut d = fixture(Topology::Circular);
        {
            let collagene_core::assembly_intent::Selector::ReferenceRange { start: value, .. } =
                &mut d.slots[0].selector
            else {
                panic!("fixture requires reference selector")
            };
            *value = 6;
        }
        {
            let collagene_core::assembly_intent::Selector::ReferenceRange { end: value, .. } =
                &mut d.slots[0].selector
            else {
                panic!("fixture requires reference selector")
            };
            *value = 10;
        } // TT|AA
        assert_eq!(selection_interval(&d, 0, 7, 9, None).unwrap(), (1, 3));
        assert_eq!(source_text(&d, 0, 0, 4, None).unwrap(), "TTAA");
        {
            let collagene_core::assembly_intent::Selector::ReferenceRange {
                reverse: value, ..
            } = &mut d.slots[0].selector
            else {
                panic!("fixture requires reference selector")
            };
            *value = true;
        }
        assert_eq!(selection_interval(&d, 0, 0, 1, None).unwrap(), (1, 2));
        assert_eq!(selection_interval(&d, 0, 0, 0, None).unwrap(), (2, 2));
        assert_eq!(source_text(&d, 0, 1, 2, None).unwrap(), "T");
        assert_eq!(source_interval(&d, 0, 1, 3).unwrap(), (7, 9));
        assert_eq!(source_interval(&d, 0, 2, 2).unwrap(), (0, 0));
        assert!(selection_interval(&d, 0, 1, 4, None).is_err());
        assert!(selection_interval(&d, 0, 8, 9, None).is_err());
        assert_eq!(
            selection_interval(&fixture(Topology::Linear), 0, 8, 8, None).unwrap(),
            (8, 8)
        );
    }
    #[test]
    fn role_bound_edits_roundtrip_and_clear_only_changed_source() {
        let d = fixture(Topology::Linear);
        let mut edited = declare(&d, 0, "3", "5", "T", None).unwrap();
        edited = declare(&edited, 1, "1", "1", "G", None).unwrap();
        let old = edited.edits.clone();
        let restored = crate::assembly_draft_model::decode(
            &crate::assembly_draft_model::encode_raw(&edited).unwrap(),
        )
        .unwrap();
        assert!(restored == edited);
        edited.slots.swap(0, 1);
        edited.methods.swap_rows(0, 1);
        edited.sync_rows();
        assert_eq!(edited.edits, old);
        edited.clear_slot_substitutions(0).unwrap();
        assert_eq!(edited.edits, vec![old[0].clone()]);
        edited.clear_slot_substitutions(1).unwrap();
        assert!(edited.edits.is_empty());
    }
}
