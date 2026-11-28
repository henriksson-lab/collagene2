//! Exact local sequence-edit navigation receipts. No view state, scientific result
//! or content-equal alias can manufacture an edit transition.
use collagene_core::documents::history::DocumentHistory;
#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) struct SourcePin {
    pub epoch: u64,
    pub document_id: String,
    pub event_id: String,
    pub revision: String,
    pub length: usize,
}
impl SourcePin {
    pub fn of(epoch: u64, document: &DocumentHistory) -> Self {
        Self {
            epoch,
            document_id: document.id.clone(),
            event_id: document.head_event().into(),
            revision: document.head_revision().into(),
            length: document.molecule.sequence.len(),
        }
    }
}
/// Boundary after a successful replacement, measured from normalized committed
/// sequence lengths rather than a raw string or clipboard byte count.
pub(crate) fn committed_caret(
    before_len: usize,
    after_len: usize,
    start: usize,
    end: usize,
) -> Option<usize> {
    if start > end || end > before_len {
        return None;
    }
    let retained = before_len.checked_sub(end - start)?;
    let inserted = after_len.checked_sub(retained)?;
    let caret = start.checked_add(inserted)?;
    (caret <= after_len).then_some(caret)
}
#[derive(Clone, Debug)]
pub(crate) struct EditTransition {
    before: SourcePin,
    after: SourcePin,
    caret: usize,
}
impl EditTransition {
    /// The caller obtains both pins from its actual history controller around one
    /// successful local replace/commit. Length is the normalized committed delta,
    /// never raw clipboard bytes or an imported producer receipt.
    pub fn attest(
        before: SourcePin,
        after: SourcePin,
        start: usize,
        end: usize,
        replacement_len: usize,
    ) -> Option<Self> {
        if before.epoch != after.epoch
            || before.document_id != after.document_id
            || start > end
            || end > before.length
        {
            return None;
        }
        let expected = before
            .length
            .checked_sub(end - start)?
            .checked_add(replacement_len)?;
        let caret = start.checked_add(replacement_len)?;
        if expected != after.length || caret > after.length {
            return None;
        }
        Some(Self {
            before,
            after,
            caret,
        })
    }
}
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum Decision {
    Unchanged,
    Clear,
    PreserveCaret(usize),
}
#[derive(Debug)]
pub(crate) struct Tracker {
    observed: Option<SourcePin>,
    pending: Option<EditTransition>,
    invalidated: bool,
}
impl Tracker {
    pub fn new(observed: Option<SourcePin>) -> Self {
        Self {
            observed,
            pending: None,
            invalidated: false,
        }
    }
    pub fn invalidate(&mut self) {
        self.pending = None;
        self.invalidated = true;
    }
    /// Coalesced key events may precede one reactive invalidation tick. Retain only
    /// a contiguous exact-head chain rooted at the last observed source.
    pub fn record(&mut self, edit: EditTransition) {
        if self.invalidated {
            return;
        }
        match self.pending.take() {
            Some(mut chain) if chain.after == edit.before => {
                chain.after = edit.after;
                chain.caret = edit.caret;
                self.pending = Some(chain);
            }
            None if self.observed.as_ref() == Some(&edit.before) => self.pending = Some(edit),
            _ => self.invalidate(),
        }
    }
    pub fn observe(&mut self, current: Option<SourcePin>) -> Decision {
        let pending = self.pending.take();
        let result = if self.invalidated {
            Decision::Clear
        } else if let Some(edit) = pending {
            if self.observed.as_ref() == Some(&edit.before) && current.as_ref() == Some(&edit.after)
            {
                Decision::PreserveCaret(edit.caret)
            } else {
                Decision::Clear
            }
        } else if self.observed == current {
            Decision::Unchanged
        } else {
            Decision::Clear
        };
        self.observed = current;
        self.invalidated = false;
        result
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use collagene_model::{Molecule, Topology};
    fn doc(id: &str) -> DocumentHistory {
        DocumentHistory::new(
            id,
            Molecule::dna("Literal editor source", "ACGT", Topology::Linear).unwrap(),
        )
        .unwrap()
    }
    fn replace(
        document: &mut DocumentHistory,
        tracker: &mut Tracker,
        start: usize,
        end: usize,
        text: &str,
    ) {
        let before = SourcePin::of(7, document);
        document
            .commit(
                vec![
                    collagene_core::documents::reference_edit::Operation::ReplaceSequence {
                        start,
                        end,
                        replacement: text.into(),
                    },
                ],
                "Local keyboard/clipboard edit",
            )
            .unwrap();
        let after = SourcePin::of(7, document);
        let n = after.length - (before.length - (end - start));
        tracker.record(EditTransition::attest(before, after, start, end, n).unwrap());
    }
    #[test]
    fn repeated_real_key_commits_preserve_new_caret_after_each_reactive_tick() {
        let mut d = doc("document-a");
        let mut tracker = Tracker::new(Some(SourcePin::of(7, &d)));
        for (n, text) in ["t", "t", "t"].into_iter().enumerate() {
            replace(&mut d, &mut tracker, n, n, text);
            assert_eq!(
                tracker.observe(Some(SourcePin::of(7, &d))),
                Decision::PreserveCaret(n + 1)
            );
        }
        assert_eq!(d.molecule.sequence.as_str(), "TTTACGT");
    }
    #[test]
    fn coalesced_edits_normalized_clipboard_cut_and_noop_use_the_actual_committed_delta() {
        let mut d = doc("document-a");
        let mut tracker = Tracker::new(Some(SourcePin::of(7, &d)));
        replace(&mut d, &mut tracker, 0, 0, "t");
        replace(&mut d, &mut tracker, 1, 1, "g");
        assert_eq!(
            tracker.observe(Some(SourcePin::of(7, &d))),
            Decision::PreserveCaret(2)
        );
        replace(&mut d, &mut tracker, 1, 3, "a c\ng");
        assert_eq!(d.molecule.sequence.as_str(), "TACGCGT");
        assert_eq!(committed_caret(6, d.molecule.sequence.len(), 1, 3), Some(4));
        assert_eq!(committed_caret(6, 0, 1, 3), None);
        assert_eq!(
            tracker.observe(Some(SourcePin::of(7, &d))),
            Decision::PreserveCaret(4)
        );
        replace(&mut d, &mut tracker, 1, 4, "");
        assert_eq!(
            tracker.observe(Some(SourcePin::of(7, &d))),
            Decision::PreserveCaret(1)
        );
        let event = d.head_event().to_string();
        replace(&mut d, &mut tracker, 0, 1, "T");
        assert_eq!(d.head_event(), event);
        assert_eq!(
            tracker.observe(Some(SourcePin::of(7, &d))),
            Decision::PreserveCaret(1)
        );
    }
    #[test]
    fn unrelated_document_writes_and_first_tick_preserve_active_selection_without_a_receipt() {
        let active = doc("document-a");
        let mut other = doc("document-b");
        let mut tracker = Tracker::new(Some(SourcePin::of(7, &active)));
        assert_eq!(
            tracker.observe(Some(SourcePin::of(7, &active))),
            Decision::Unchanged
        );
        other
            .commit(
                vec![
                    collagene_core::documents::reference_edit::Operation::ReplaceSequence {
                        start: 0,
                        end: 0,
                        replacement: "G".into(),
                    },
                ],
                "Other tab changed",
            )
            .unwrap();
        assert_eq!(
            tracker.observe(Some(SourcePin::of(7, &active))),
            Decision::Unchanged
        );
    }
    #[test]
    fn external_edit_undo_redo_and_restore_clear_instead_of_reusing_old_caret() {
        let mut d = doc("document-a");
        let mut tracker = Tracker::new(Some(SourcePin::of(7, &d)));
        replace(&mut d, &mut tracker, 1, 1, "T");
        assert_eq!(
            tracker.observe(Some(SourcePin::of(7, &d))),
            Decision::PreserveCaret(2)
        );
        d.commit(
            vec![
                collagene_core::documents::reference_edit::Operation::ReplaceSequence {
                    start: 0,
                    end: 0,
                    replacement: "G".into(),
                },
            ],
            "External annotation/editor route",
        )
        .unwrap();
        assert_eq!(tracker.observe(Some(SourcePin::of(7, &d))), Decision::Clear);
        assert!(d.undo().unwrap());
        tracker.invalidate();
        assert_eq!(tracker.observe(Some(SourcePin::of(7, &d))), Decision::Clear);
        assert!(d.redo().unwrap());
        tracker.invalidate();
        assert_eq!(tracker.observe(Some(SourcePin::of(7, &d))), Decision::Clear);
        assert_eq!(tracker.observe(Some(SourcePin::of(8, &d))), Decision::Clear);
    }
    #[test]
    fn equal_content_other_document_or_other_event_and_interleaved_external_heads_cannot_consume_receipt()
     {
        let mut d = doc("document-a");
        let mut tracker = Tracker::new(Some(SourcePin::of(7, &d)));
        replace(&mut d, &mut tracker, 0, 0, "T");
        let clone = d.fork_event(d.head_event(), "document-b").unwrap();
        assert_eq!(clone.molecule, d.molecule);
        assert_eq!(
            tracker.observe(Some(SourcePin::of(7, &clone))),
            Decision::Clear
        );
        let mut d = doc("same-content-events");
        d.commit_with_source_events(
            vec![
                collagene_core::documents::reference_edit::Operation::Rename {
                    name: d.molecule.name.clone(),
                },
            ],
            "Explicit new event at equal content",
            &[d.head_event().into()],
        )
        .unwrap();
        let mut tracker = Tracker::new(Some(SourcePin::of(7, &d)));
        let before = SourcePin::of(7, &d);
        d.commit_with_source_events(
            vec![
                collagene_core::documents::reference_edit::Operation::Rename {
                    name: d.molecule.name.clone(),
                },
            ],
            "Another source event at equal content",
            &[d.head_event().into()],
        )
        .unwrap();
        let after = SourcePin::of(7, &d);
        assert_eq!(before.revision, after.revision);
        assert_ne!(before.event_id, after.event_id);
        assert_eq!(tracker.observe(Some(after)), Decision::Clear);
        let mut d = doc("interleaved");
        let mut tracker = Tracker::new(Some(SourcePin::of(7, &d)));
        replace(&mut d, &mut tracker, 0, 0, "T");
        d.commit(
            vec![
                collagene_core::documents::reference_edit::Operation::ReplaceSequence {
                    start: 0,
                    end: 0,
                    replacement: "G".into(),
                },
            ],
            "External intervening commit",
        )
        .unwrap();
        replace(&mut d, &mut tracker, 0, 0, "A");
        assert_eq!(tracker.observe(Some(SourcePin::of(7, &d))), Decision::Clear);
    }
    #[test]
    fn forged_geometry_stale_sources_or_navigation_cannot_attest_valid_local_continuity() {
        let d = doc("a");
        let before = SourcePin::of(7, &d);
        let mut after = before.clone();
        after.length += 1;
        assert!(EditTransition::attest(before.clone(), after.clone(), 0, 0, 2).is_none());
        assert!(EditTransition::attest(before.clone(), after.clone(), 3, 8, 5).is_none());
        after.document_id = "other".into();
        assert!(EditTransition::attest(before.clone(), after, 0, 0, 1).is_none());
        let mut tracker = Tracker::new(Some(before.clone()));
        tracker.invalidate();
        assert_eq!(tracker.observe(Some(before)), Decision::Clear);
    }
}
