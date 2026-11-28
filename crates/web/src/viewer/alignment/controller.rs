use collagene_core::material::lineage::EventBinding;
use collagene_core::verification;
use collagene_core::{
    alignment, documents::workspace, material::lineage as run_lineage,
    sequence::similarity as sequence_similarity,
};

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
