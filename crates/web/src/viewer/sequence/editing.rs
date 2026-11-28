//! Source-pinned sequence editing, local caret receipts and undo/redo transitions.
use crate::viewer::sequence::selection::ViewerSelection;
use crate::workspace::session::WorkspaceSession;
use collagene_model::Molecule;
use leptos::prelude::*;
#[derive(Clone, Copy)]
pub(crate) struct SequenceEditor {
    pub(crate) sequence_transitions:
        StoredValue<crate::viewer::sequence::sequence_edit_transition::Tracker>,
    pub(crate) sequence_on_edit: Callback<(String, usize, usize, String), Result<(), String>>,
    pub(crate) sequence_on_history: Callback<bool, Result<(), String>>,
}
impl SequenceEditor {
    pub(crate) fn new(session: WorkspaceSession, viewer: ViewerSelection) -> Self {
        let WorkspaceSession {
            documents,
            active,
            workspace_epoch,
            workspace_restoring,
            ..
        } = session;
        let ViewerSelection {
            molecule,
            selection,
            selected_feature,
            search_regions,
            ..
        } = viewer;
        let sequence_transitions = StoredValue::new(
            crate::viewer::sequence::sequence_edit_transition::Tracker::new(
                documents.with_untracked(|docs| {
                    docs.get(active.get_untracked()).map(|doc| {
                        crate::viewer::sequence::sequence_edit_transition::SourcePin::of(
                            workspace_epoch.get_untracked(),
                            doc,
                        )
                    })
                }),
            ),
        );
        let sequence_on_edit = Callback::new(
            move |(base, a, b, replacement): (String, usize, usize, String)| {
                let current = molecule.get_untracked();
                if current.revision() != base {
                    return Err("The sequence changed before this edit could be applied.".into());
                }
                let before = documents
                    .with_untracked(|docs| {
                        docs.get(active.get_untracked()).map(|doc| {
                            crate::viewer::sequence::sequence_edit_transition::SourcePin::of(
                                workspace_epoch.get_untracked(),
                                doc,
                            )
                        })
                    })
                    .ok_or_else(|| "The source document is no longer open.".to_string())?;
                session.commit_document(
                    &base,
                    vec![
                        collagene_core::documents::reference_edit::Operation::ReplaceSequence {
                            start: a,
                            end: b,
                            replacement,
                        },
                    ],
                )?;
                let after = documents.with_untracked(|docs| {
                    docs.get(active.get_untracked()).map(|doc| {
                        crate::viewer::sequence::sequence_edit_transition::SourcePin::of(
                            workspace_epoch.get_untracked(),
                            doc,
                        )
                    })
                });
                // Normalized committed lengths, rather than raw clipboard bytes,
                // determine the resulting caret after the actual transaction.
                let replacement_len = after.as_ref().and_then(|after| {
                    before
                        .length
                        .checked_sub(b.saturating_sub(a))
                        .and_then(|retained| after.length.checked_sub(retained))
                });
                let transition = after.zip(replacement_len).and_then(|(after, n)| {
                    crate::viewer::sequence::sequence_edit_transition::EditTransition::attest(
                        before, after, a, b, n,
                    )
                });
                sequence_transitions.update_value(|tracker| {
                    if let Some(transition) = transition {
                        tracker.record(transition);
                    } else {
                        tracker.invalidate();
                    }
                });
                Ok(())
            },
        );
        let sequence_on_history = Callback::new(move |redo: bool| {
            if workspace_restoring.get_untracked() {
                return Err("Workspace restoration is in progress.".into());
            }
            let mut result = Ok(false);
            documents.update(|docs| {
                let doc = &mut docs[active.get_untracked()];
                result = if redo { doc.redo() } else { doc.undo() }.map_err(|e| e.to_string());
            });
            let changed = result?;
            if changed {
                sequence_transitions.update_value(|tracker| tracker.invalidate());
                batch(|| {
                    selection.set(None);
                    selected_feature.set(None);
                    search_regions.set(None);
                });
            }
            Ok(())
        });
        Self {
            sequence_transitions,
            sequence_on_edit,
            sequence_on_history,
        }
    }
}
