//! Literal internal contribution intent, independent of primers and preparation.
//!
//! UI entry uses one-based half-open boundaries: `3..5` replaces contribution
//! bases3 and4, `3..3` inserts before base3, and `length+1` appends at the end.
//! Stored offsets are zero-based in the ORIGINAL contribution's construct
//! orientation, before junction trims. Source documents are never modified.
use crate::assembly_draft_model::Draft;
use collagene_core::{Alphabet, Sequence, Topology, assembly_edit_intent::Edit};

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
    let length = s
        .end
        .checked_sub(s.start)
        .filter(|n| *n > 0)
        .ok_or("The original contribution interval is invalid.")?;
    if n == 0
        || s.start >= n
        || length > n
        || capture.molecule.sequence.alphabet() != Alphabet::Dna
        || (capture.molecule.topology == Topology::Linear && s.end > n)
    {
        return Err("Choose a valid original DNA contribution of at most one source turn.".into());
    }
    Ok((s, length))
}

/// Exact source interval in construct orientation; never derives letters from
/// already edited targets, proposed primers or reconstructed products.
pub fn source_text(draft: &Draft, slot: usize, start: usize, end: usize) -> Result<String, String> {
    let (s, length) = contribution(draft, slot)?;
    if start > end || end > length || end - start > 100_000 {
        return Err("Choose at most100000 bases within the original contribution.".into());
    }
    let capture = s.captured.as_ref().unwrap();
    let n = capture.molecule.sequence.len();
    let (a, b) = if s.reverse {
        (length - end, length - start)
    } else {
        (start, end)
    };
    let bases = capture.molecule.sequence.as_str().as_bytes();
    let text: String = (a..b)
        .map(|offset| bases[(s.start + offset) % n] as char)
        .collect();
    let text = if s.reverse && !text.is_empty() {
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
) -> Result<(usize, usize), String> {
    let (s, length) = contribution(draft, slot)?;
    let capture = s.captured.as_ref().unwrap();
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
    if capture.molecule.topology == Topology::Circular && unrolled < s.start {
        unrolled = unrolled
            .checked_add(n)
            .ok_or("Selection coordinates overflow.")?;
    }
    let offset = unrolled
        .checked_sub(s.start)
        .ok_or("Select inside the original contribution.")?;
    let limit = offset
        .checked_add(count)
        .filter(|v| *v <= length)
        .ok_or("Selection crosses the chosen contribution boundary.")?;
    let interval = if s.reverse {
        (length - limit, length - offset)
    } else {
        (offset, limit)
    };
    source_text(draft, slot, interval.0, interval.1)?;
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
    source_text(draft, slot, start, end)?;
    let capture = s.captured.as_ref().unwrap();
    let n = capture.molecule.sequence.len();
    let offset = if s.reverse { length - end } else { start };
    let raw = s
        .start
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
) -> Result<Draft, String> {
    if draft.alternatives.is_some() {
        return Err("Internal contribution edits currently require an ordinary construct.".into());
    }
    let start = boundary(start)?;
    let end = boundary(end)?;
    let expected = source_text(draft, slot, start, end)?;
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
            slot.start = 0;
            slot.end = 8;
        }
        draft.sync_rows();
        draft.validate().unwrap();
        draft
    }
    #[test]
    fn declarations_preserve_originals_and_reference_boundaries() {
        let original = fixture(Topology::Linear);
        let mut d = declare(&original, 0, "3", "5", "tag").unwrap();
        assert_eq!(d.version, 7);
        assert_eq!(d.edits[0].expected, "CC");
        assert_eq!(d.edits[0].replacement, "TAG");
        assert!(d.slots == original.slots && d.methods.rows == original.methods.rows);
        d = declare(&d, 0, "9", "9", "AC").unwrap();
        assert_eq!((d.edits[1].start, d.edits[1].end), (8, 8));
        assert_eq!(d.edits[1].expected, "");
        assert!(declare(&d, 0, "4", "4", "G").is_err());
        assert!(declare(&d, 0, "5", "5", "G").is_err());
        assert!(declare(&d, 0, "3", "5", "").is_err());
        let deletion = declare(&original, 0, "2", "4", "").unwrap();
        assert_eq!(deletion.edits[0].expected, "AC");
        assert_eq!(deletion.edits[0].replacement, "");
        assert!(declare(&original, 0, "0", "1", "T").is_err());
        assert!(declare(&original, 0, "2", "4", "AC").is_err());
    }
    #[test]
    fn precise_circular_reverse_ranges_and_carets_map_to_original_contribution() {
        let mut d = fixture(Topology::Circular);
        d.slots[0].start = 6;
        d.slots[0].end = 10; // TT|AA
        assert_eq!(selection_interval(&d, 0, 7, 9).unwrap(), (1, 3));
        assert_eq!(source_text(&d, 0, 0, 4).unwrap(), "TTAA");
        d.slots[0].reverse = true;
        assert_eq!(selection_interval(&d, 0, 0, 1).unwrap(), (1, 2));
        assert_eq!(selection_interval(&d, 0, 0, 0).unwrap(), (2, 2));
        assert_eq!(source_text(&d, 0, 1, 2).unwrap(), "T");
        assert_eq!(source_interval(&d, 0, 1, 3).unwrap(), (7, 9));
        assert_eq!(source_interval(&d, 0, 2, 2).unwrap(), (0, 0));
        assert!(selection_interval(&d, 0, 1, 4).is_err());
        assert!(selection_interval(&d, 0, 8, 9).is_err());
        assert_eq!(
            selection_interval(&fixture(Topology::Linear), 0, 8, 8).unwrap(),
            (8, 8)
        );
    }
    #[test]
    fn role_bound_edits_roundtrip_and_clear_only_changed_source() {
        let d = fixture(Topology::Linear);
        let mut edited = declare(&d, 0, "3", "5", "T").unwrap();
        edited = declare(&edited, 1, "1", "1", "G").unwrap();
        let old = edited.edits.clone();
        let restored = crate::assembly_draft_model::decode_version(
            &crate::assembly_draft_model::encode_raw(&edited).unwrap(),
            7,
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
