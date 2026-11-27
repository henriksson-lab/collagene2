//! Editing declarations only; primers and simulated products never supply intent.
use crate::{
    assembly_draft_model::{Draft, Slot},
    assembly_role_ids::Role,
    cloning_methods_model::{Mode, OverlapInput, OverlapPreparation},
};
use collagene_core::assembly_substitution_intent::Substitution;

/// Coordinate in the captured original, with its letter in construct orientation.
pub fn source_base(draft: &Draft, slot: usize, offset: usize) -> Result<(usize, char), String> {
    let s = draft
        .slots
        .get(slot)
        .ok_or("The source slot is unavailable.")?;
    let capture = s.captured.as_ref().ok_or("Pin a source region first.")?;
    let n = capture.molecule.sequence.len();
    let length = s
        .end
        .checked_sub(s.start)
        .filter(|&v| v > 0 && v <= n)
        .ok_or("The source interval is invalid.")?;
    if n == 0
        || s.start >= n
        || offset >= length
        || (capture.molecule.topology == collagene_core::Topology::Linear && s.end > n)
    {
        return Err("The edit position is outside this contribution.".into());
    }
    let position = (s.start
        + if s.reverse {
            length - 1 - offset
        } else {
            offset
        })
        % n;
    let mut base = capture.molecule.sequence.as_str().as_bytes()[position];
    if s.reverse {
        base = match base {
            b'A' => b'T',
            b'T' => b'A',
            b'C' => b'G',
            b'G' => b'C',
            _ => base,
        };
    }
    if !b"ACGT".contains(&base) {
        return Err("This position has an ambiguous source base; choose a concrete reference before declaring a substitution.".into());
    }
    Ok((position, base as char))
}

/// Convert an exact original-source base selection into a contribution offset.
pub fn selection_offset(draft: &Draft, slot: usize, position: usize) -> Result<usize, String> {
    let s = draft
        .slots
        .get(slot)
        .ok_or("The source slot is unavailable.")?;
    let capture = s.captured.as_ref().ok_or("Pin a source region first.")?;
    let n = capture.molecule.sequence.len();
    if n == 0 || position >= n || s.start >= n {
        return Err("Selection is outside the original source.".into());
    }
    let length = s
        .end
        .checked_sub(s.start)
        .filter(|&v| v > 0 && v <= n)
        .ok_or("The source interval is invalid.")?;
    let forward = if position >= s.start {
        position - s.start
    } else {
        n - s.start + position
    };
    if forward >= length {
        return Err("Select one base inside the chosen source contribution.".into());
    }
    let offset = if s.reverse {
        length - 1 - forward
    } else {
        forward
    };
    source_base(draft, slot, offset)?;
    Ok(offset)
}

pub fn declare(
    draft: &Draft,
    slot: usize,
    position: &str,
    replacement: &str,
) -> Result<Draft, String> {
    if draft.alternatives.is_some() {
        return Err("Internal substitutions currently belong to an ordinary construct.".into());
    }
    if position.len() > 20 || !position.bytes().all(|b| b.is_ascii_digit()) {
        return Err("Enter a one-based contribution position.".into());
    }
    let offset = position
        .parse::<usize>()
        .ok()
        .and_then(|p| p.checked_sub(1))
        .ok_or("Contribution positions start at 1.")?;
    let replacement = replacement.to_ascii_uppercase();
    if replacement.len() != 1 || !b"ACGT".contains(&replacement.as_bytes()[0]) {
        return Err("Choose one replacement DNA base.".into());
    }
    let (_, expected) = source_base(draft, slot, offset)?;
    if replacement.as_bytes()[0] == expected as u8 {
        return Err("Choose a different base from the original reference.".into());
    }
    let roles = draft.methods.role_ids(draft.slots.len())?;
    let role = roles
        .get(slot)
        .ok_or("The fragment role is unavailable.")?
        .role
        .clone();
    if draft
        .substitutions
        .iter()
        .any(|e| e.role == role && e.offset == offset)
    {
        return Err("An edit already exists at this contribution position. Remove it before declaring another.".into());
    }
    if draft.substitutions.len() >= 256 {
        return Err("A construct supports at most 256 internal substitutions.".into());
    }
    let mut next = draft.clone();
    next.methods.source_roles = Some(roles);
    next.substitutions.push(Substitution {
        role,
        offset,
        expected: expected.to_string(),
        replacement,
    });
    next.sync_rows();
    next.validate()?;
    Ok(next)
}

/// Extend ordinary named roles without renumbering existing edit bindings.
pub fn append_slot(draft: &Draft) -> Result<Draft, String> {
    if draft.alternatives.is_some()
        || draft.slots.len() >= 12
        || draft.methods.mode == Mode::Tailing
    {
        return Err("No additional source slot is available for this construct.".into());
    }
    let mut next = draft.clone();
    let roles = next.methods.role_ids(next.slots.len())?;
    if next.methods.source_roles.is_some() {
        let suffix = (1..=1000)
            .find(|n| {
                !roles.iter().any(|r| {
                    r.input_id == format!("source-{n}") || r.role == format!("fragment-{n}")
                })
            })
            .ok_or("No unused fragment identity is available.")?;
        let mut extended = roles;
        extended.push(Role {
            input_id: format!("source-{suffix}"),
            role: format!("fragment-{suffix}"),
        });
        next.methods.rows.push(OverlapInput {
            preparation: OverlapPreparation::ExactPcr,
            ..Default::default()
        });
        next.methods.source_roles = Some(extended);
        if let Some(form) = &mut next.methods.golden_gate_pool {
            form.choices.resize(next.methods.rows.len(), 0);
            form.candidate = None;
        }
    } else {
        next.methods.push_row(OverlapInput {
            preparation: OverlapPreparation::ExactPcr,
            ..Default::default()
        })?;
    }
    next.slots
        .push(Slot::empty(&format!("Fragment {}", next.slots.len() + 1)));
    next.active_slot = next.slots.len() - 1;
    next.junction_edits.clear();
    next.sync_rows();
    next.validate()?;
    Ok(next)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{assembly_draft_model::CapturedSource, source_identity::SourceInfo};
    use collagene_core::{Molecule, Topology, run_lineage::EventBinding, workspace::Workspace};
    use std::sync::Arc;
    fn draft(topology: Topology) -> Draft {
        let molecule = Molecule::dna("Original source", "AACCGGTT", topology).unwrap();
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
                alphabet: molecule.sequence.alphabet(),
                topology,
            },
            molecule,
        });
        let mut d = Draft::default();
        for s in &mut d.slots {
            s.captured = Some(capture.clone());
            s.start = 0;
            s.end = 8;
        }
        d.sync_rows();
        d.validate().unwrap();
        d
    }
    #[test]
    fn exact_base_selection_maps_reverse_and_origin_crossing_contributions() {
        let mut d = draft(Topology::Circular);
        d.slots[0].start = 6;
        d.slots[0].end = 10;
        // Source TT|AA, reverse complement TT|AA. Original coordinate0 maps
        // to offset2 forward, but offset1 in construct reverse orientation.
        assert_eq!(selection_offset(&d, 0, 0).unwrap(), 2);
        assert_eq!(source_base(&d, 0, 2).unwrap(), (0, 'A'));
        d.slots[0].reverse = true;
        d.sync_rows();
        assert_eq!(selection_offset(&d, 0, 0).unwrap(), 1);
        assert_eq!(source_base(&d, 0, 1).unwrap(), (0, 'T'));
        assert!(selection_offset(&d, 0, 3).is_err());
        assert!(source_base(&d, 0, 4).is_err());
    }
    #[test]
    fn declarations_and_new_slots_keep_original_role_pins_through_reordering() {
        let d = draft(Topology::Linear);
        let mut edited = declare(&d, 0, "3", "T").unwrap();
        assert_eq!(edited.version, crate::assembly_draft_model::CURRENT_VERSION);
        assert_eq!(edited.substitutions[0].expected, "C");
        assert!(edited.slots == d.slots);
        assert!(edited.methods.rows == d.methods.rows);
        let role = edited.substitutions[0].role.clone();
        edited = append_slot(&edited).unwrap();
        assert_eq!(edited.slots.len(), 3);
        assert_eq!(edited.substitutions[0].role, role);
        edited.slots.swap(0, 1);
        edited.methods.swap_rows(0, 1);
        edited.sync_rows();
        assert_eq!(edited.methods.role_ids(3).unwrap()[1].role, role);
        edited.clear_slot_substitutions(0).unwrap();
        assert_eq!(edited.substitutions.len(), 1);
        edited.clear_slot_substitutions(1).unwrap();
        assert!(edited.substitutions.is_empty());
    }
    #[test]
    fn invalid_or_duplicate_declarations_do_not_change_the_draft() {
        let d = declare(&draft(Topology::Linear), 0, "3", "T").unwrap();
        let before = d.clone();
        for (position, replacement) in [
            ("3", "G"),
            ("0", "T"),
            ("9", "T"),
            ("2e", "T"),
            ("2", "A"),
            ("2", "NN"),
        ] {
            assert!(declare(&d, 0, position, replacement).is_err());
            assert!(d == before);
        }
    }
}
