use collagene_core::documents::transfer;
use collagene_core::material::lineage::EventBinding;

use leptos::prelude::*;
use wasm_bindgen::JsCast;

use crate::workspace::session::WorkspaceSession;

use crate::workspace::session::current_event_binding;

use crate::viewer::sequence::selection::ViewerSelection;
use crate::workspace::shell_state::ShellState;

#[derive(Clone, Copy)]
pub(crate) struct ComparisonController {
    pub(crate) session: WorkspaceSession,
    pub(crate) viewer: ViewerSelection,
    pub(crate) shell: ShellState,
}
impl ComparisonController {
    pub(crate) fn saved_comparisons(
        &self,
        candidate: RwSignal<Option<transfer::Candidate>>,
        stored_candidates: RwSignal<Vec<transfer::Candidate>>,
    ) -> Memo<Vec<collagene_core::alignment::history::Summary>> {
        let WorkspaceSession {
            documents,
            active,
            attachments,
            workspace_name,
            ..
        } = self.session;
        let ShellState { message, .. } = self.shell;
        Memo::new(move |_| {
            let frame = WorkspaceSession::snapshot(
                &documents.get(),
                active.get_untracked(),
                candidate.get_untracked().as_ref(),
                &stored_candidates.get_untracked(),
                &attachments.get(),
                &workspace_name.get_untracked(),
            );
            match frame.and_then(|frame| {
                collagene_core::alignment::history::list(&frame).map_err(|e| e.to_string())
            }) {
                Ok(summaries) => summaries,
                Err(error) => {
                    message.set(format!("Saved comparisons need attention: {error}"));
                    Vec::new()
                }
            }
        })
    }
    pub(crate) fn save_comparison(
        &self,
        contig_snapshot: Callback<
            (),
            Result<collagene_core::documents::workspace::Workspace, String>,
        >,
    ) -> Callback<collagene_core::alignment::history::Record, Result<(), String>> {
        let WorkspaceSession {
            workspace_restoring,
            attachments,
            ..
        } = self.session;
        let ShellState { message, .. } = self.shell;
        Callback::new(
            move |record: collagene_core::alignment::history::Record| -> Result<(), String> {
                if workspace_restoring.get_untracked() {
                    return Err("Wait for the workspace update to finish.".into());
                }
                let frame = contig_snapshot.run(())?;
                let (staged, _) = collagene_core::alignment::history::capture(
                    &frame,
                    &record.reference,
                    &record.query,
                    &record.report,
                )
                .map_err(|e| e.to_string())?;
                staged.validate().map_err(|e| e.to_string())?;
                attachments.set(staged.attachments);
                message.set("Comparison retained with its exact source versions. Reopening shows historical evidence until explicitly recomputed.".into());
                Ok(())
            },
        )
    }
    pub(crate) fn reopen_comparison(
        &self,
        contig_snapshot: Callback<
            (),
            Result<collagene_core::documents::workspace::Workspace, String>,
        >,
    ) -> Callback<String, Result<collagene_core::alignment::history::Record, String>> {
        let WorkspaceSession {
            workspace_restoring,
            ..
        } = self.session;
        Callback::new(
            move |hash: String| -> Result<collagene_core::alignment::history::Record, String> {
                if workspace_restoring.get_untracked() {
                    return Err("Wait for the workspace update to finish.".into());
                }
                collagene_core::alignment::history::load(&contig_snapshot.run(())?, &hash)
                    .map_err(|e| e.to_string())
            },
        )
    }
    pub(crate) fn guard_similarity(
        &self,
    ) -> Callback<crate::viewer::alignment::sequence_similarity_model::Capture, Result<(), String>>
    {
        let WorkspaceSession {
            documents,
            workspace_restoring,
            workspace_epoch,
            ..
        } = self.session;
        Callback::new(
        move |capture: crate::viewer::alignment::sequence_similarity_model::Capture| -> Result<(), String> {
            if workspace_restoring.get() {
                return Err("The workspace is being restored.".into());
            }
            documents.with(|docs| {
                let d = docs
                    .iter()
                    .find(|d| d.id == capture.query.document_id)
                    .ok_or("The pinned query is unavailable.")?;
                let m = d
                    .revision_molecule_ref(&capture.query.event_id)
                    .map_err(|e| e.to_string())?;
                if m.revision() != capture.query.revision {
                    return Err("The pinned query changed.".into());
                }
                let current = capture
                    .subjects
                    .iter()
                    .map(|p| {
                        docs.get(p.index)
                            .map(|d| (p.index, current_event_binding(d)))
                            .ok_or_else(|| "A search subject is unavailable.".to_owned())
                    })
                    .collect::<Result<Vec<_>, String>>()?;
                capture.validate_current(workspace_epoch.get(), &capture.query, &current)
            })
        },
    )
    }
    pub(crate) fn prepare_similarity(
        &self,
    ) -> Callback<
        (
            EventBinding,
            Vec<usize>,
            collagene_core::sequence::similarity::Options,
        ),
        Result<crate::viewer::alignment::sequence_similarity_model::Prepared, String>,
    > {
        let session = self.session;
        let WorkspaceSession {
            documents,
            workspace_restoring,
            workspace_epoch,
            ..
        } = self.session;
        let ViewerSelection { molecule, .. } = self.viewer;
        Callback::new(
            move |(query, indices, options): (
                EventBinding,
                Vec<usize>,
                collagene_core::sequence::similarity::Options,
            )|
                  -> Result<
                crate::viewer::alignment::sequence_similarity_model::Prepared,
                String,
            > {
                if workspace_restoring.get_untracked()
                    || indices.len() > collagene_core::sequence::similarity::MAX_SUBJECTS
                {
                    return Err("Choose 1–128 documents in a ready workspace.".into());
                }
                documents.with_untracked(|docs| {
                    let query_document = docs
                        .iter()
                        .find(|d| d.id == query.document_id)
                        .ok_or("Pinned query document is unavailable.")?;
                    let molecule = query_document
                        .revision_molecule_ref(&query.event_id)
                        .map_err(|e| e.to_string())?;
                    let subjects = indices
                        .into_iter()
                        .map(|index| {
                            docs.get(index)
                                .map(|d| (index, current_event_binding(d), &d.molecule))
                                .ok_or_else(|| "A selected document is unavailable.".to_owned())
                        })
                        .collect::<Result<Vec<_>, String>>()?;
                    crate::viewer::alignment::sequence_similarity_model::prepare(
                        workspace_epoch.get_untracked(),
                        query,
                        molecule,
                        &subjects,
                        options,
                    )
                })
            },
        )
    }
    pub(crate) fn navigate_comparison(
        &self,
        pending_comparison_source: RwSignal<
            Option<crate::viewer::alignment::pairwise_source_view::Snapshot>,
        >,
    ) -> Callback<crate::viewer::alignment::compare::Navigation, Result<(), String>> {
        let session = self.session;
        let WorkspaceSession {
            documents,
            active,
            workspace_restoring,
            workspace_epoch,
            ..
        } = self.session;
        let ViewerSelection {
            selection,
            selected_feature,
            search_regions,
            molecule,
            has_open_document,
            ..
        } = self.viewer;
        let ShellState {
            modal, modal_error, ..
        } = self.shell;
        Callback::new(
            move |request: crate::viewer::alignment::compare::Navigation| -> Result<(), String> {
                if workspace_restoring.get_untracked()
                    || request.epoch != workspace_epoch.get_untracked()
                {
                    return Err("The workspace changed; reopen the comparison.".into());
                }
                let (index, current, molecule) =
                    documents.with_untracked(|docs| -> Result<_, String> {
                        let index = docs
                            .iter()
                            .position(|doc| doc.id == request.event.document_id)
                            .ok_or("The retained source document is unavailable.")?;
                        let doc = &docs[index];
                        let molecule = doc
                            .revision_molecule(&request.event.event_id)
                            .map_err(|e| e.to_string())?;
                        if molecule.revision() != request.event.revision
                            || request.spans.is_empty()
                            || request.spans.len() > 2
                            || request
                                .spans
                                .iter()
                                .any(|&(start, end)| start > end || end > molecule.sequence.len())
                        {
                            return Err(
                                "The comparison position differs from its exact source version."
                                    .into(),
                            );
                        }
                        Ok((index, current_event_binding(doc) == request.event, molecule))
                    })?;
                if !current || request.spans.len() > 1 {
                    pending_comparison_source.set(Some(
                        crate::viewer::alignment::pairwise_source_view::Snapshot {
                            source: request.event,
                            epoch: request.epoch,
                            molecule: std::sync::Arc::new(molecule),
                            spans: request.spans,
                        },
                    ));
                    modal_error.set(String::new());
                    modal.set("comparison_source");
                    return Ok(());
                }
                let apply = move |deferred: bool| {
                    if workspace_epoch.try_get_untracked() != Some(request.epoch)
                        || workspace_restoring.try_get_untracked().unwrap_or(true)
                        || !documents.with_untracked(|docs| {
                            docs.get(active.get_untracked())
                                .is_some_and(|doc| current_event_binding(doc) == request.event)
                        })
                    {
                        return;
                    }
                    if deferred
                        && (selection.get_untracked().is_some()
                            || selected_feature.get_untracked().is_some())
                    {
                        return;
                    }
                    batch(|| {
                        selected_feature.set(None);
                        search_regions.set(None);
                        selection.set(request.spans.first().copied());
                    });
                };
                if active.get_untracked() == index && has_open_document.get_untracked() {
                    apply(false);
                } else {
                    active.set(index);
                    let callback =
                        wasm_bindgen::closure::Closure::once_into_js(move || apply(true));
                    web_sys::window()
                        .ok_or("No browser window")?
                        .request_animation_frame(callback.unchecked_ref())
                        .map_err(|_| "Comparison navigation could not be scheduled.")?;
                }
                Ok(())
            },
        )
    }
}

/// Retained alignment revisions and pending alignment/source dialogs.
#[derive(Clone, Copy)]
pub(crate) struct ComparisonState {
    pub(crate) alignments:
        RwSignal<Vec<crate::viewer::alignment::alignment_view::AlignmentRevision>>,
    pub(crate) pending_alignment:
        RwSignal<Option<crate::viewer::alignment::alignment_view::PendingAlignment>>,
    pub(crate) pending_comparison_source:
        RwSignal<Option<crate::viewer::alignment::pairwise_source_view::Snapshot>>,
}
impl ComparisonState {
    pub(crate) fn new() -> Self {
        Self {
            alignments: RwSignal::new(Vec::<
                crate::viewer::alignment::alignment_view::AlignmentRevision,
            >::new()),
            pending_alignment: RwSignal::new(
                None::<crate::viewer::alignment::alignment_view::PendingAlignment>,
            ),
            pending_comparison_source: RwSignal::new(
                None::<crate::viewer::alignment::pairwise_source_view::Snapshot>,
            ),
        }
    }
    pub(crate) fn restore(
        &self,
        alignments: Vec<crate::viewer::alignment::alignment_view::AlignmentRevision>,
    ) {
        self.alignments.set(alignments);
        self.pending_alignment.set(None);
        self.pending_comparison_source.set(None);
    }
}
