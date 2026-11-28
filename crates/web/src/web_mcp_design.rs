//! Browser-agent design suggestions reuse normal declaration commands. Only the
//! mounted UI can apply a proposal; this module never supplies verified material.
use crate::{
    assembly_draft_model::Draft,
    assembly_edit_model::{self, SourceReview},
    web_mcp_assembly::{self, Pin},
};
use serde::{Deserialize, Serialize};
#[derive(Clone, Debug, Deserialize, Serialize)]
#[serde(tag = "operation", rename_all = "snake_case", deny_unknown_fields)]
pub(crate) enum Change {
    Replace {
        role: String,
        start: usize,
        end: usize,
        expected: String,
        replacement: String,
    },
    Substitution {
        role: String,
        offset: usize,
        expected: String,
        replacement: String,
    },
}
#[derive(Clone, Debug, Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub(crate) struct Request {
    pub expected: Pin,
    pub changes: Vec<Change>,
}
#[derive(Clone, Debug, PartialEq, Eq, Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub(crate) struct Token {
    pub owner: u64,
    pub generation: u32,
    pub revision: String,
}
#[derive(Clone, Deserialize)]
#[serde(deny_unknown_fields)]
pub(crate) struct DiscardRequest {
    pub expected: Pin,
    pub proposal: Token,
}
#[derive(Clone)]
pub(crate) struct Proposal {
    expected: Pin,
    changes: Vec<Change>,
    next: Draft,
    revision: String,
}
impl Change {
    pub(crate) fn description(&self) -> String {
        match self {
            Self::Replace {
                role,
                start,
                end,
                expected,
                replacement,
            } => format!(
                "{role}: original contribution[{start},{end}) {expected:?} → {replacement:?}"
            ),
            Self::Substitution {
                role,
                offset,
                expected,
                replacement,
            } => format!("{role}: original contribution base{offset} {expected} → {replacement}"),
        }
    }
    fn role(&self) -> &str {
        match self {
            Self::Replace { role, .. } | Self::Substitution { role, .. } => role,
        }
    }
}
impl Proposal {
    pub(crate) fn prepare(
        draft: &Draft,
        epoch: u64,
        request: Request,
        review: Option<&SourceReview>,
    ) -> Result<Self, String> {
        web_mcp_assembly::require_pin(draft, epoch, &request.expected)?;
        draft.validate()?;
        crate::assembly_alternatives::require_active_choice_match(draft)?;
        if request.changes.is_empty() || request.changes.len() > 8 {
            return Err("Propose 1..8 literal contribution changes per review.".into());
        }
        let roles = draft.methods.role_ids(draft.slots.len())?;
        let mut next = draft.clone();
        for change in &request.changes {
            if change.role().is_empty() || change.role().len() > 128 {
                return Err("Use an exact bounded contribution role from assembly context.".into());
            }
            let slot = roles
                .iter()
                .position(|r| r.role == change.role())
                .ok_or("The proposed contribution role is not in this draft.")?;
            let (start, end, expected, replacement) = match change {
                Change::Replace {
                    start,
                    end,
                    expected,
                    replacement,
                    ..
                } => (*start, *end, expected, replacement),
                Change::Substitution {
                    offset,
                    expected,
                    replacement,
                    ..
                } => (
                    *offset,
                    offset
                        .checked_add(1)
                        .ok_or("Contribution offset overflow")?,
                    expected,
                    replacement,
                ),
            };
            if start > end
                || end - start > 512
                || expected.len() != end - start
                || replacement.len() > 512
                || !expected
                    .bytes()
                    .chain(replacement.bytes())
                    .all(|b| b"ACGTRYSWKMBDHVN".contains(&b))
            {
                return Err("Use exact uppercase IUPAC literals of at most 512 bases in bounded original contribution coordinates.".into());
            }
            let actual = assembly_edit_model::source_text(&next, slot, start, end, review)?;
            if actual != *expected {
                return Err("The expected literal differs from the exact original contribution; no draft changes were applied.".into());
            }
            next = match change {
                Change::Replace { .. } => assembly_edit_model::declare(
                    &next,
                    slot,
                    &(start + 1).to_string(),
                    &(end + 1).to_string(),
                    replacement,
                    review,
                )?,
                Change::Substitution { .. } => crate::assembly_substitution_model::declare(
                    &next,
                    slot,
                    &(start + 1).to_string(),
                    replacement,
                    review,
                )?,
            };
        }
        next.validate()?;
        let revision = blake3::hash(&serde_json::to_vec(&request).map_err(|e| e.to_string())?)
            .to_hex()
            .to_string();
        Ok(Self {
            expected: request.expected,
            changes: request.changes,
            next,
            revision,
        })
    }
    pub(crate) fn revision(&self) -> &str {
        &self.revision
    }
    pub(crate) fn changes(&self) -> &[Change] {
        &self.changes
    }
    pub(crate) fn is_current(&self, current: &Draft, epoch: u64) -> bool {
        web_mcp_assembly::require_pin(current, epoch, &self.expected).is_ok()
    }
    pub(crate) fn apply(&self, current: &Draft, epoch: u64) -> Result<Draft, String> {
        web_mcp_assembly::require_pin(current, epoch, &self.expected)?;
        self.next.validate()?;
        let mut next = self.next.clone();
        next.active_slot = current.active_slot;
        Ok(next)
    }
}
#[cfg(test)]
mod tests {
    use super::*;
    use crate::{assembly_draft_model::CapturedSource, source_identity::SourceInfo};
    use collagene_core::history::DocumentHistory;
    use collagene_model::{Molecule, Topology};
    use std::sync::Arc;
    fn draft() -> Draft {
        let source = DocumentHistory::new(
            "original",
            Molecule::dna("Original", "AACCGGTT", Topology::Linear).unwrap(),
        )
        .unwrap();
        let capture = Arc::new(CapturedSource {
            source: SourceInfo {
                binding: collagene_core::run_lineage::EventBinding {
                    document_id: source.id.clone(),
                    event_id: source.head_event().into(),
                    revision: source.head_revision().into(),
                },
                name: source.molecule.name.clone(),
                length: 8,
                alphabet: source.molecule.sequence.alphabet(),
                topology: Topology::Linear,
            },
            molecule: source.molecule.clone(),
        });
        let mut draft = Draft::default();
        for slot in &mut draft.slots {
            slot.captured = Some(capture.clone());
            slot.selector = collagene_core::assembly_intent::Selector::ReferenceRange {
                start: 0,
                end: 8,
                reverse: false,
            };
        }
        draft.sync_rows();
        draft.validate().unwrap();
        draft
    }
    #[test]
    fn proposed_changes_remain_private_until_explicit_apply_and_wrong_literal_is_atomic() {
        let original = draft();
        let role = original.methods.role_ids(original.slots.len()).unwrap()[0]
            .role
            .clone();
        let request = Request {
            expected: web_mcp_assembly::pin(&original, 4).unwrap(),
            changes: vec![
                Change::Substitution {
                    role: role.clone(),
                    offset: 1,
                    expected: "A".into(),
                    replacement: "T".into(),
                },
                Change::Replace {
                    role: role.clone(),
                    start: 4,
                    end: 6,
                    expected: "GG".into(),
                    replacement: "ACA".into(),
                },
            ],
        };
        let before = crate::assembly_draft_model::encode(&original).unwrap();
        let proposal = Proposal::prepare(&original, 4, request.clone(), None).unwrap();
        assert_eq!(
            crate::assembly_draft_model::encode(&original).unwrap(),
            before
        );
        let applied = proposal.apply(&original, 4).unwrap();
        assert_eq!(applied.substitutions.len(), 1);
        assert_eq!(applied.edits.len(), 1);
        assert!(applied.slots == original.slots);
        assert!(applied.methods.rows == original.methods.rows);
        let mut wrong = request;
        let Change::Replace { expected, .. } = &mut wrong.changes[1] else {
            unreachable!()
        };
        *expected = "GT".into();
        assert!(Proposal::prepare(&original, 4, wrong, None).is_err());
        assert_eq!(
            crate::assembly_draft_model::encode(&original).unwrap(),
            before
        );
        assert!(proposal.apply(&original, 5).is_err());
        let mut navigation = original.clone();
        navigation.active_slot = 1;
        assert!(proposal.is_current(&navigation, 4));
        assert_eq!(proposal.apply(&navigation, 4).unwrap().active_slot, 1);
        let mut changed = original.clone();
        changed.methods.overlaps.push('1');
        assert!(proposal.apply(&changed, 4).is_err());
    }
    #[test]
    fn role_orientation_and_edit_conflicts_use_existing_declaration_guards() {
        let mut original = draft();
        original.slots[0].selector = collagene_core::assembly_intent::Selector::ReferenceRange {
            start: 0,
            end: 4,
            reverse: true,
        };
        original.sync_rows();
        let role = original.methods.role_ids(original.slots.len()).unwrap()[0]
            .role
            .clone();
        let change = Change::Replace {
            role: role.clone(),
            start: 0,
            end: 2,
            expected: "GG".into(),
            replacement: "TT".into(),
        };
        let request = Request {
            expected: web_mcp_assembly::pin(&original, 2).unwrap(),
            changes: vec![change.clone()],
        };
        let proposed = Proposal::prepare(&original, 2, request.clone(), None).unwrap();
        assert_eq!(
            proposed.apply(&original, 2).unwrap().edits[0].expected,
            "GG"
        );
        let mut overlap = request.clone();
        overlap.changes.push(change);
        assert!(Proposal::prepare(&original, 2, overlap, None).is_err());
        let mut unknown = request;
        let Change::Replace { role, .. } = &mut unknown.changes[0] else {
            unreachable!()
        };
        *role = "unknown".into();
        assert!(Proposal::prepare(&original, 2, unknown, None).is_err());
    }
}
