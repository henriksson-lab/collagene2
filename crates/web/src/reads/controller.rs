use collagene_core::documents::transfer;
use collagene_core::material::lineage::EventBinding;
use collagene_core::planning::plan::DesignPlan;
use collagene_core::verification;
use collagene_core::verification::evidence::VerificationReport;
use collagene_core::{
    documents::{history, workspace},
    material::lineage as run_lineage,
    review::reads::{review as read_review, trace_validation},
};
use collagene_model::Molecule;
use leptos::prelude::*;
use wasm_bindgen::JsCast;

use crate::workspace::session::RestoredWorkspace;
use crate::workspace::session::WorkspaceSession;

use crate::workspace::session::MAX_WORKSPACE_BYTES;

use crate::workspace::session::current_event_binding;

use crate::workspace::session::document_bytes;

use crate::viewer::sequence::selection::ViewerSelection;
use crate::workspace::shell_state::ShellState;

#[derive(Clone)]
pub(crate) struct PendingReadReview {
    pub(crate) request: collagene_core::review::reads::review::Request,
    pub(crate) owner: collagene_core::material::lineage::EventBinding,
    pub(crate) epoch: u64,
}
#[derive(Clone, Copy)]
pub(crate) struct ReadController {
    pub(crate) session: WorkspaceSession,
    pub(crate) viewer: ViewerSelection,
    pub(crate) shell: ShellState,
}
impl ReadController {
    pub(crate) fn fork_run_event(
        &self,
        contig_snapshot: Callback<
            (),
            Result<collagene_core::documents::workspace::Workspace, String>,
        >,
        last_plan: RwSignal<Option<(DesignPlan, Molecule)>>,
        last_plan_events: RwSignal<
            Option<(
                String,
                Vec<collagene_core::material::lineage::InputEventBinding>,
            )>,
        >,
        last_plan_target: RwSignal<
            Option<(String, u64, collagene_core::material::lineage::EventBinding)>,
        >,
        pending_run_event: RwSignal<
            Option<std::sync::Arc<crate::history::run_event_model::Snapshot>>,
        >,
        plan_report: RwSignal<Option<VerificationReport>>,
    ) -> Callback<crate::history::run_event_model::Request, Result<(), String>> {
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
            molecule,
            ..
        } = self.viewer;
        let ShellState {
            panel,
            message,
            modal,
            ..
        } = self.shell;
        Callback::new(
            move |request: crate::history::run_event_model::Request| -> Result<(), String> {
                if workspace_restoring.get_untracked()
                    || request.epoch != workspace_epoch.get_untracked()
                {
                    return Err("The workspace changed; reopen the cloning run version.".into());
                }
                if !pending_run_event.with_untracked(|pending| {
                    pending
                        .as_ref()
                        .is_some_and(|snapshot| snapshot.request == request)
                }) {
                    return Err(
                        "The displayed run version changed; open it again before forking.".into(),
                    );
                }
                let frame = contig_snapshot.run(())?;
                let prefix = format!(
                    "run-version-{}",
                    request
                        .binding
                        .revision
                        .chars()
                        .take(16)
                        .collect::<String>()
                );
                let id = (1..=1001)
                    .map(|suffix| format!("{prefix}-{suffix}"))
                    .find(|id| !frame.documents.iter().any(|doc| doc.id == *id))
                    .ok_or("No document identity is available for this copy.")?;
                let staged = crate::history::run_event_model::fork(&frame, &request, &id)?;
                let state = RestoredWorkspace::from_frame(&staged)?;
                if state.documents.len() > 1000
                    || state
                        .documents
                        .iter()
                        .map(|doc| document_bytes(&doc.molecule))
                        .sum::<usize>()
                        > MAX_WORKSPACE_BYTES
                {
                    return Err(
                        "This version copy exceeds the workspace document budget; nothing changed."
                            .into(),
                    );
                }
                batch(|| {
                    documents.set(state.documents);
                    active.set(state.active);
                    selected_feature.set(None);
                    selection.set(None);
                    plan_report.set(None);
                    last_plan.set(None);
                    last_plan_target.set(None);
                    last_plan_events.set(None);
                    pending_run_event.set(None);
                    modal.set("");
                    panel.set("history");
                    message.set("Opened an editable copy of the exact recorded version. Original document heads and recorded assessments are preserved.".into());
                });
                Ok(())
            },
        )
    }
    pub(crate) fn record_read_check(
        &self,
        pending_read_check: RwSignal<Option<crate::reads::trace_validation::RecordRequest>>,
        read_check_reference: Memo<EventBinding>,
    ) -> Callback<crate::reads::trace_validation::RecordRequest, Result<(), String>> {
        let session = self.session;
        let WorkspaceSession {
            documents,
            workspace_restoring,
            workspace_epoch,
            ..
        } = self.session;
        let ShellState { modal, .. } = self.shell;
        Callback::new(
            move |request: crate::reads::trace_validation::RecordRequest| -> Result<(), String> {
                if workspace_restoring.get_untracked()
                    || request.epoch != workspace_epoch.get_untracked()
                {
                    return Err("The workspace changed. Recheck the reads before saving.".into());
                }
                if request.reference != read_check_reference.get_untracked() {
                    return Err(
                        "The construct reference changed. Recheck the reads before saving.".into(),
                    );
                }
                documents.with_untracked(|docs| {
                    if std::iter::once(&request.reference)
                        .chain(request.reads.iter())
                        .any(|binding| {
                            !docs
                                .iter()
                                .any(|doc| current_event_binding(doc) == *binding)
                        })
                    {
                        return Err(
                            "A reference or read event changed. Recheck before saving.".to_string()
                        );
                    }
                    Ok(())
                })?;
                pending_read_check.set(Some(request));
                modal.set("read_check_capture");
                Ok(())
            },
        )
    }
    pub(crate) fn navigate_check_reference(
        &self,
        candidate: RwSignal<Option<transfer::Candidate>>,
        stored_candidates: RwSignal<Vec<transfer::Candidate>>,
    ) -> Callback<crate::reads::trace_validation::ReferenceNavigation, Result<(), String>> {
        let session = self.session;
        let WorkspaceSession {
            documents,
            active,
            workspace_restoring,
            workspace_epoch,
            workspace_name,
            attachments,
            ..
        } = self.session;
        let ViewerSelection {
            selection,
            selected_feature,
            molecule,
            has_open_document,
            ..
        } = self.viewer;
        Callback::new(
        move |request: crate::reads::trace_validation::ReferenceNavigation| -> Result<(), String> {
            if workspace_restoring.get_untracked()
                || request.epoch != workspace_epoch.get_untracked()
            {
                return Err("The workspace changed; reopen the read check.".into());
            }
            let mut docs = documents.get_untracked();
            let index = docs
                .iter()
                .position(|doc| doc.id == request.reference.document_id)
                .ok_or("The reference document is unavailable.")?;
            let original = &docs[index];
            let molecule = original
                .revision_molecule(&request.reference.event_id)
                .map_err(|e| e.to_string())?;
            if molecule.revision() != request.reference.revision
                || request
                    .selection
                    .is_some_and(|(a, b)| a > b || b > molecule.sequence.len())
            {
                return Err(
                    "The reference event or requested coordinates differ from the saved check."
                        .into(),
                );
            }
            let index = if original.head_event() == request.reference.event_id {
                index
            } else {
                if docs.len() >= 1000 {
                    return Err(
                        "Opening the historical reference exceeds this workspace's document limit."
                            .into(),
                    );
                }
                let prefix = format!("read-check-reference-{}", &request.reference.revision[..16]);
                let id = (1..=1001)
                    .map(|n| format!("{prefix}-{n}"))
                    .find(|id| !docs.iter().any(|doc| doc.id == *id))
                    .ok_or("No historical reference identity is available.")?;
                let fork = original
                    .fork_event(&request.reference.event_id, id)
                    .map_err(|e| e.to_string())?;
                let index = docs.len();
                docs.push(fork);
                index
            };
            let id = docs[index].id.clone();
            let event = docs[index].head_event().to_owned();
            WorkspaceSession::snapshot(
                &docs,
                index,
                candidate.get_untracked().as_ref(),
                &stored_candidates.get_untracked(),
                &attachments.get_untracked(),
                &workspace_name.get_untracked(),
            )?;
            let apply = move |deferred: bool| {
                if workspace_epoch.try_get_untracked() != Some(request.epoch)
                    || workspace_restoring.try_get_untracked().unwrap_or(true)
                {
                    return;
                }
                if !documents.with_untracked(|docs| {
                    docs.get(active.get_untracked())
                        .is_some_and(|doc| doc.id == id && doc.head_event() == event)
                }) {
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
                    selection.set(request.selection);
                });
            };
            documents.set(docs);
            if active.get_untracked() == index && has_open_document.get_untracked() {
                apply(false);
            } else {
                active.set(index);
                let frame = wasm_bindgen::closure::Closure::once_into_js(move || apply(true));
                if let Some(window) = web_sys::window() {
                    window
                        .request_animation_frame(frame.unchecked_ref())
                        .map_err(|_| "Reference navigation could not be scheduled.")?;
                }
            }
            Ok(())
        },
    )
    }
    pub(crate) fn navigate_read_call(
        &self,
        candidate: RwSignal<Option<transfer::Candidate>>,
        stored_candidates: RwSignal<Vec<transfer::Candidate>>,
    ) -> Callback<crate::reads::contigs::ReadNavigation, ()> {
        let session = self.session;
        let WorkspaceSession {
            documents,
            active,
            workspace_restoring,
            workspace_epoch,
            workspace_name,
            attachments,
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
            panel,
            message,
            modal,
            ..
        } = self.shell;
        Callback::new(move |request: crate::reads::contigs::ReadNavigation| {
            if workspace_restoring.get_untracked()
                || workspace_epoch.get_untracked() != request.epoch
            {
                return;
            }
            let prepared = (|| -> Result<_, String> {
                let mut docs = documents.get_untracked();
                let source_index = docs
                    .iter()
                    .position(|d| d.id == request.binding.document_id)
                    .ok_or("The source read document is unavailable.")?;
                let source = &docs[source_index];
                let molecule = source
                    .revision_molecule(&request.binding.event_id)
                    .map_err(|e| e.to_string())?;
                if molecule.revision() != request.binding.revision
                    || request.call_index >= molecule.sequence.len()
                {
                    return Err("The source call no longer matches its retained event.".into());
                }
                if source.head_event() == request.binding.event_id {
                    return Ok((
                        source_index,
                        source.id.clone(),
                        source.head_event().to_string(),
                        molecule.revision(),
                    ));
                }
                if docs.len() >= 1000
                    || docs
                        .iter()
                        .map(|d| document_bytes(&d.molecule))
                        .sum::<usize>()
                        + document_bytes(&molecule)
                        > MAX_WORKSPACE_BYTES
                {
                    return Err(
                        "Opening the retained read would exceed this workspace's capacity.".into(),
                    );
                }
                let base = format!("read-history-{}", &request.binding.revision[..16]);
                let mut suffix = 1usize;
                let id = loop {
                    let id = format!("{base}-{suffix}");
                    if docs.iter().all(|d| d.id != id) {
                        break id;
                    }
                    suffix += 1;
                };
                let fork = source
                    .fork_event(&request.binding.event_id, id.clone())
                    .map_err(|e| e.to_string())?;
                let event = fork.head_event().to_string();
                let revision = fork.molecule.revision();
                let next = docs.len();
                docs.push(fork);
                let frame = WorkspaceSession::snapshot(
                    &docs,
                    next,
                    candidate.get_untracked().as_ref(),
                    &stored_candidates.get_untracked(),
                    &attachments.get_untracked(),
                    &workspace_name.get_untracked(),
                )?;
                frame.validate().map_err(|e| e.to_string())?;
                documents.set(docs);
                message.set("Opened the retained read as a new document; its original source head is unchanged.".into());
                Ok((next, id, event, revision))
            })();
            let (index, id, event, revision) = match prepared {
                Ok(v) => v,
                Err(error) => {
                    message.set(error);
                    return;
                }
            };
            let expected_epoch = request.epoch;
            let at = request.call_index;
            let apply = move |deferred: bool| {
                if workspace_epoch.try_get_untracked() != Some(expected_epoch)
                    || workspace_restoring.try_get_untracked().unwrap_or(true)
                {
                    return;
                }
                let Some(current) = active.try_get_untracked() else {
                    return;
                };
                if !documents.with_untracked(|docs| {
                    docs.get(current).is_some_and(|d| {
                        d.id == id && d.head_event() == event && d.molecule.revision() == revision
                    })
                }) {
                    return;
                }
                if deferred
                    && (selected_feature.get_untracked().is_some()
                        || selection.get_untracked().is_some()
                        || search_regions.get_untracked().is_some())
                {
                    return;
                }
                batch(|| {
                    selected_feature.set(None);
                    selection.set(Some((at, at + 1)));
                    search_regions.set(None);
                });
            };
            if active.get_untracked() == index && has_open_document.get_untracked() {
                apply(false);
            } else {
                active.set(index);
                let frame = wasm_bindgen::closure::Closure::once_into_js(move || apply(true));
                if let Some(window) = web_sys::window() {
                    let _ = window.request_animation_frame(frame.unchecked_ref());
                }
            }
            panel.set("features");
            modal.set("");
        })
    }
}

#[derive(Clone, Copy)]
pub(crate) struct ReadReviewState {
    pub(crate) read_reference: RwSignal<Option<collagene_core::material::lineage::EventBinding>>,
    pub(crate) pending_read_check: RwSignal<Option<crate::reads::trace_validation::RecordRequest>>,
    pub(crate) reviewed_check: RwSignal<Option<String>>,
    pub(crate) review_seed: RwSignal<Option<collagene_core::review::reads::review::Request>>,
    pub(crate) pending_read_review: RwSignal<Option<PendingReadReview>>,
}
impl ReadReviewState {
    pub(crate) fn new() -> Self {
        let read_reference = RwSignal::new(None::<collagene_core::material::lineage::EventBinding>);
        let pending_read_check =
            RwSignal::new(None::<crate::reads::trace_validation::RecordRequest>);
        let reviewed_check = RwSignal::new(None::<String>);
        let review_seed = RwSignal::new(None::<collagene_core::review::reads::review::Request>);
        let pending_read_review = RwSignal::new(None::<PendingReadReview>);
        Self {
            read_reference,
            pending_read_check,
            reviewed_check,
            review_seed,
            pending_read_review,
        }
    }
}
pub(crate) struct ReadReviewActions {
    pub(crate) read_check_snapshot:
        Callback<(), Result<collagene_core::documents::workspace::Workspace, String>>,
    pub(crate) restore_read_check: Callback<collagene_core::documents::workspace::Workspace>,
    pub(crate) open_read_review: Callback<String, Result<(), String>>,
    pub(crate) reopen_read_review: Callback<String, Result<(), String>>,
    pub(crate) save_read_review:
        Callback<collagene_core::review::reads::review::Request, Result<(), String>>,
    pub(crate) read_review_snapshot:
        Callback<(), Result<collagene_core::documents::workspace::Workspace, String>>,
    pub(crate) restore_read_review: Callback<collagene_core::review::reads::review::Commit>,
}
impl ReadController {
    pub(crate) fn review_actions(
        &self,
        state: ReadReviewState,
        contig_snapshot: Callback<
            (),
            Result<collagene_core::documents::workspace::Workspace, String>,
        >,
        read_check_reference: Memo<collagene_core::material::lineage::EventBinding>,
        read_check_attachment_available: Callback<String, bool>,
        apply_restored: Callback<RestoredWorkspace>,
    ) -> ReadReviewActions {
        let WorkspaceSession {
            documents,
            workspace_epoch,
            workspace_restoring,
            ..
        } = self.session;
        let ShellState {
            panel,
            message,
            modal,
            ..
        } = self.shell;
        let ReadReviewState {
            read_reference,
            pending_read_check,
            reviewed_check,
            review_seed,
            pending_read_review,
        } = state;
        let read_check_snapshot = Callback::new(
            move |()| -> Result<collagene_core::documents::workspace::Workspace, String> {
                let pending = pending_read_check
                    .get_untracked()
                    .ok_or("No sequencing check is pending.")?;
                if pending.epoch != workspace_epoch.get_untracked()
                    || pending.reference != read_check_reference.get_untracked()
                    || !documents.with_untracked(|docs| {
                        std::iter::once(&pending.reference)
                            .chain(pending.reads.iter())
                            .all(|binding| {
                                docs.iter()
                                    .any(|doc| current_event_binding(doc) == *binding)
                            })
                    })
                {
                    return Err(
                        "The reference or read events changed. Recheck before saving.".into(),
                    );
                }
                contig_snapshot.run(())
            },
        );
        let restore_read_check = Callback::new(
            move |frame: collagene_core::documents::workspace::Workspace| {
                let pinned = read_reference.get_untracked();
                match RestoredWorkspace::from_frame(&frame) {
                    Ok(state) => batch(|| {
                        apply_restored.run(state);
                        read_reference.set(pinned);
                        panel.set("reads");
                    }),
                    Err(error) => message.set(error),
                }
            },
        );
        let open_read_review = Callback::new(move |hash: String| -> Result<(), String> {
            if workspace_restoring.get_untracked() {
                return Err("Wait for the workspace operation to finish.".into());
            }
            let frame = contig_snapshot.run(())?;
            let saved = collagene_core::review::reads::checks::load(&frame, &hash)
                .map_err(|e| e.to_string())?;
            if saved.reference.document_id != read_check_reference.get_untracked().document_id {
                return Err("Choose the construct that owns this check.".into());
            }
            reviewed_check.set(Some(hash));
            review_seed.set(None);
            modal.set("read_review");
            Ok(())
        });
        let reopen_read_review = Callback::new(move |hash: String| -> Result<(), String> {
            if workspace_restoring.get_untracked() {
                return Err("Wait for the workspace operation to finish.".into());
            }
            let frame = contig_snapshot.run(())?;
            let saved = collagene_core::review::reads::review::load(&frame, &hash)
                .map_err(|e| e.to_string())?;
            if saved.reference.document_id != read_check_reference.get_untracked().document_id {
                return Err(
                    "Choose the original checked construct before reopening its decisions.".into(),
                );
            }
            reviewed_check.set(Some(saved.request.attachment_hash.clone()));
            review_seed.set(Some(saved.request));
            modal.set("read_review");
            Ok(())
        });
        let save_read_review = Callback::new(
            move |request: collagene_core::review::reads::review::Request| -> Result<(), String> {
                if workspace_restoring.get_untracked() {
                    return Err("Wait for the workspace operation to finish.".into());
                }
                if reviewed_check.get_untracked().as_deref()
                    != Some(request.attachment_hash.as_str())
                {
                    return Err(
                        "The selected read check changed; preview the decisions again.".into(),
                    );
                }
                let frame = contig_snapshot.run(())?;
                let saved =
                    collagene_core::review::reads::checks::load(&frame, &request.attachment_hash)
                        .map_err(|e| e.to_string())?;
                let owner = read_check_reference.get_untracked();
                if saved.reference.document_id != owner.document_id {
                    return Err(
                        "The checked construct changed; preview the decisions again.".into(),
                    );
                }
                pending_read_review.set(Some(PendingReadReview {
                    request,
                    owner,
                    epoch: workspace_epoch.get_untracked(),
                }));
                modal.set("read_review_commit");
                Ok(())
            },
        );
        let read_review_snapshot = Callback::new(
            move |()| -> Result<collagene_core::documents::workspace::Workspace, String> {
                let pending = pending_read_review
                    .get_untracked()
                    .ok_or("No read review is pending.")?;
                if pending.epoch != workspace_epoch.get_untracked()
                    || pending.owner != read_check_reference.get_untracked()
                    || !read_check_attachment_available.run(pending.request.attachment_hash.clone())
                {
                    return Err(
                        "The construct/check context changed; preview the decisions again.".into(),
                    );
                }
                contig_snapshot.run(())
            },
        );
        let restore_read_review = Callback::new(
            move |commit: collagene_core::review::reads::review::Commit| {
                let pinned = commit.target.or_else(|| read_reference.get_untracked());
                match RestoredWorkspace::from_frame(&commit.workspace) {
                    Ok(state) => batch(|| {
                        apply_restored.run(state);
                        read_reference.set(pinned);
                        panel.set("reads");
                    }),
                    Err(error) => message.set(error),
                }
            },
        );
        ReadReviewActions {
            read_check_snapshot,
            restore_read_check,
            open_read_review,
            reopen_read_review,
            save_read_review,
            read_review_snapshot,
            restore_read_review,
        }
    }
}
