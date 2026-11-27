//! Editing declarations only; primers and simulated products never supply intent.
use crate::{
    assembly_draft_model::{Draft, Slot},
    assembly_role_ids::Role,
    cloning_methods_model::{Mode, OverlapInput, OverlapPreparation},
};
use collagene_core::assembly_substitution_intent::Substitution;

/// Coordinate in the captured original, with its letter in construct orientation.
pub fn source_base(
    draft: &Draft,
    slot: usize,
    offset: usize,
    review: Option<&crate::assembly_edit_model::SourceReview>,
) -> Result<(collagene_core::assembly_intent::ContributionPosition, char), String> {
    use collagene_core::assembly_intent::{ContributionPosition, Selector};
    let end = offset
        .checked_add(1)
        .ok_or("Contribution offset overflow")?;
    let text = crate::assembly_edit_model::source_text(draft, slot, offset, end, review)?;
    let base = *text
        .as_bytes()
        .first()
        .ok_or("Choose an original contribution base")?;
    if !b"ACGT".contains(&base) {
        return Err("Choose a concrete original base before declaring a substitution".into());
    }
    let s = &draft.slots[slot];
    let position = match &s.selector {
        Selector::ReferenceRange {
            start,
            end,
            reverse,
        } => {
            let length = end - start;
            let n = s
                .captured
                .as_ref()
                .ok_or("Pin a source first")?
                .source
                .length;
            ContributionPosition::Reference {
                position: (start
                    + if *reverse {
                        length - 1 - offset
                    } else {
                        offset
                    })
                    % n,
                reverse: *reverse,
            }
        }
        Selector::MaterialStrandRange { .. } => review
            .ok_or("Review the original actual strand first")?
            .get(draft, slot)?
            .source_position(offset)
            .map_err(|e| e.to_string())?,
    };
    Ok((position, base as char))
}

pub fn selection_offset(
    draft: &Draft,
    slot: usize,
    position: usize,
    review: Option<&crate::assembly_edit_model::SourceReview>,
) -> Result<usize, String> {
    let end = position.checked_add(1).ok_or("Selection overflow")?;
    let (start, stop) =
        crate::assembly_edit_model::selection_interval(draft, slot, position, end, review)?;
    if stop != start + 1 {
        return Err("Choose one original source base".into());
    }
    source_base(draft, slot, start, review)?;
    Ok(start)
}

pub fn declare(
    draft: &Draft,
    slot: usize,
    position: &str,
    replacement: &str,
    review: Option<&crate::assembly_edit_model::SourceReview>,
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
    let (_, expected) = source_base(draft, slot, offset, review)?;
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
            s.selector = collagene_core::assembly_intent::Selector::ReferenceRange {
                start: 0,
                end: 8,
                reverse: false,
            };
        }
        d.sync_rows();
        d.validate().unwrap();
        d
    }
    #[test]
    fn exact_base_selection_maps_reverse_and_origin_crossing_contributions() {
        let mut d = draft(Topology::Circular);
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
        }
        // Source TT|AA, reverse complement TT|AA. Original coordinate0 maps
        // to offset2 forward, but offset1 in construct reverse orientation.
        assert_eq!(selection_offset(&d, 0, 0, None).unwrap(), 2);
        assert_eq!(
            source_base(&d, 0, 2, None).unwrap(),
            (
                collagene_core::assembly_intent::ContributionPosition::Reference {
                    position: 0,
                    reverse: false
                },
                'A'
            )
        );
        {
            let collagene_core::assembly_intent::Selector::ReferenceRange {
                reverse: value, ..
            } = &mut d.slots[0].selector
            else {
                panic!("fixture requires reference selector")
            };
            *value = true;
        }
        d.sync_rows();
        assert_eq!(selection_offset(&d, 0, 0, None).unwrap(), 1);
        assert_eq!(
            source_base(&d, 0, 1, None).unwrap(),
            (
                collagene_core::assembly_intent::ContributionPosition::Reference {
                    position: 0,
                    reverse: true
                },
                'T'
            )
        );
        assert!(selection_offset(&d, 0, 3, None).is_err());
        assert!(source_base(&d, 0, 4, None).is_err());
    }
    #[test]
    fn declarations_and_new_slots_keep_original_role_pins_through_reordering() {
        let d = draft(Topology::Linear);
        let mut edited = declare(&d, 0, "3", "T", None).unwrap();
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
        let d = declare(&draft(Topology::Linear), 0, "3", "T", None).unwrap();
        let before = d.clone();
        for (position, replacement) in [
            ("3", "G"),
            ("0", "T"),
            ("9", "T"),
            ("2e", "T"),
            ("2", "A"),
            ("2", "NN"),
        ] {
            assert!(declare(&d, 0, position, replacement, None).is_err());
            assert!(d == before);
        }
    }
}
