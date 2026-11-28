use collagene_core::documents::transfer;
use collagene_core::verification;
use collagene_core::{documents::workspace, sequence::reference_transfer};

use leptos::prelude::*;

use crate::workspace::session::WorkspaceSession;

use crate::workspace::session::MAX_WORKSPACE_BYTES;

use crate::workspace::session::current_event_binding;

use crate::workspace::session::document_bytes;

use crate::viewer::sequence::selection::ViewerSelection;
use crate::workspace::shell_state::ShellState;

#[derive(Clone, Copy)]
pub(crate) struct ReferenceController {
    pub(crate) session: WorkspaceSession,
    pub(crate) viewer: ViewerSelection,
    pub(crate) shell: ShellState,
}
impl ReferenceController {
    pub(crate) fn reference_copy_selection(
        &self,
    ) -> Callback<(), Result<collagene_core::sequence::reference_transfer::Selection, String>> {
        let session = self.session;
        let ViewerSelection {
            selection,
            selected_feature,
            search_regions,
            molecule,
            ..
        } = self.viewer;
        Callback::new(
            move |()| -> Result<collagene_core::sequence::reference_transfer::Selection, String> {
                use collagene_core::sequence::reference_transfer::{Part, Selection};
                let source = molecule.get_untracked();
                if let Some(index) = selected_feature.get_untracked() {
                    let feature = source
                        .features
                        .get(index)
                        .ok_or("The selected annotation changed.")?;
                    return Selection::from_location(&feature.location, source.sequence.len())
                        .map_err(|e| e.to_string());
                }
                let parts = if let Some((revision, spans)) = search_regions
                    .get_untracked()
                    .filter(|(revision, spans)| *revision == source.revision() && !spans.is_empty())
                {
                    let _ = revision;
                    spans
                } else {
                    vec![
                        selection
                            .get_untracked()
                            .ok_or("Select reference bases or an exact annotation to copy.")?,
                    ]
                };
                Ok(Selection {
                    parts: parts
                        .into_iter()
                        .map(|(start, end)| Part {
                            start,
                            end,
                            reverse: false,
                        })
                        .collect(),
                    reverse: false,
                })
            },
        )
    }
    pub(crate) fn reference_paste_guard(
        &self,
        internal_clipboard: RwSignal<
            Option<std::sync::Arc<crate::workspace::reference_transfer_view::Copied>>,
        >,
    ) -> Callback<crate::workspace::reference_transfer_view::PasteOpen, Result<(), String>> {
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
            ..
        } = self.viewer;
        Callback::new(
        move |opened: crate::workspace::reference_transfer_view::PasteOpen| -> Result<(), String> {
            if workspace_restoring.get_untracked()
                || opened.epoch != workspace_epoch.get_untracked()
                || opened.copied.epoch != opened.epoch
                || !internal_clipboard.with_untracked(|current| {
                    current
                        .as_ref()
                        .is_some_and(|current| std::sync::Arc::ptr_eq(current, &opened.copied))
                })
            {
                return Err(
                    "The copied reference or workspace changed; copy and preview again.".into(),
                );
            }
            if selected_feature.get_untracked().is_some()
                || selection.get_untracked() != Some((opened.start, opened.end))
                || search_regions.with_untracked(|regions| {
                    regions.as_ref().is_some_and(|(_, spans)| spans.len() > 1)
                })
            {
                return Err("The destination selection changed. Choose one interval or insertion caret and preview again.".into());
            }
            documents.with_untracked(|docs| {
                let destination = docs
                    .get(active.get_untracked())
                    .ok_or("The destination document is unavailable.")?;
                if current_event_binding(destination) != opened.destination
                    || destination.molecule.revision() != opened.destination_molecule.revision()
                {
                    return Err("The exact destination version changed; preview again.".into());
                }
                let source = docs
                    .iter()
                    .find(|doc| doc.id == opened.copied.source.document_id)
                    .ok_or("The copied source document is no longer retained.")?;
                let original = source
                    .revision_molecule(&opened.copied.source.event_id)
                    .map_err(|e| e.to_string())?;
                if original.revision() != opened.copied.source.revision
                    || original.revision() != opened.copied.source_molecule.revision()
                {
                    return Err(
                        "The copied source differs from its retained original event.".into(),
                    );
                }
                Ok(())
            })
        },
    )
    }
    pub(crate) fn apply_reference_paste(
        &self,
        candidate: RwSignal<Option<transfer::Candidate>>,
        pending_reference_paste: RwSignal<
            Option<crate::workspace::reference_transfer_view::PasteOpen>,
        >,
        reference_paste_guard: Callback<
            crate::workspace::reference_transfer_view::PasteOpen,
            Result<(), String>,
        >,
        stored_candidates: RwSignal<Vec<transfer::Candidate>>,
    ) -> Callback<crate::workspace::reference_transfer_view::Apply, Result<(), String>> {
        let session = self.session;
        let WorkspaceSession {
            documents,
            active,
            workspace_name,
            attachments,
            ..
        } = self.session;
        let ViewerSelection {
            selection,
            selected_feature,
            search_regions,
            molecule,
            ..
        } = self.viewer;
        let ShellState { message, modal, .. } = self.shell;
        Callback::new(
            move |value: crate::workspace::reference_transfer_view::Apply| -> Result<(), String> {
                reference_paste_guard.run(value.opened.clone())?;
                let opened = &value.opened;
                let mut staged_docs = documents.get_untracked();
                let index = active.get_untracked();
                let source = staged_docs
                    .iter()
                    .find(|doc| doc.id == opened.copied.source.document_id)
                    .ok_or("The copied source document is unavailable.")?;
                let source_molecule = source
                    .revision_molecule(&opened.copied.source.event_id)
                    .map_err(|e| e.to_string())?;
                let destination = staged_docs
                    .get(index)
                    .ok_or("The destination is unavailable.")?;
                let request = collagene_core::sequence::reference_transfer::RequestPaste {
                    source: source_molecule,
                    destination: destination.molecule.clone(),
                    selection: opened.copied.fragment.selection.clone(),
                    fragment_revision: opened.copied.fragment.context_revision.clone(),
                    start: opened.start,
                    end: opened.end,
                };
                let expected =
                    collagene_core::sequence::reference_transfer::paste_revision(&request)
                        .map_err(|e| e.to_string())?;
                value
                    .report
                    .validate_binding(&expected)
                    .map_err(|e| e.to_string())?;
                if value.report.source_revision != opened.copied.source.revision
                    || value.report.destination_revision != opened.destination.revision
                    || value.report.fragment_revision != opened.copied.fragment.context_revision
                    || value.report.start != opened.start
                    || value.report.end != opened.end
                {
                    return Err(
                        "The paste result differs from the reviewed source/destination selection."
                            .into(),
                    );
                }
                let destination = staged_docs
                    .get_mut(index)
                    .ok_or("The destination is unavailable.")?;
                destination.commit_with_source_events(value.report.product.clone(),
            "Annotated reference paste from an exact source event; no physical cloning reaction is asserted",
            &[opened.copied.source.event_id.clone()]).map_err(|e|e.to_string())?;
                if staged_docs.len() > 1000
                    || staged_docs
                        .iter()
                        .map(|doc| document_bytes(&doc.molecule))
                        .sum::<usize>()
                        > MAX_WORKSPACE_BYTES
                {
                    return Err(
                    "The reference paste exceeds the workspace document budget; nothing changed."
                        .into(),
                );
                }
                let receipt=serde_json::to_vec(&serde_json::json!({
            "version":1,"source":opened.copied.source,"destination_before":opened.destination,
            "destination_after":current_event_binding(&staged_docs[index]),"context_revision":value.report.context_revision,
            "fragment_revision":value.report.fragment_revision,"selection":value.report.selection,
            "start":value.report.start,"end":value.report.end,"inserted_start":value.report.inserted_start,"inserted_end":value.report.inserted_end,
            "warnings":value.report.warnings,"model":"Reference edit with annotation projection; no physical end chemistry or cloning success asserted"
        })).map_err(|e|e.to_string())?;
                let mut assets = attachments.get_untracked();
                assets.insert(
                    blake3::hash(&receipt).to_hex().to_string(),
                    workspace::WorkspaceAttachment {
                        kind: "reference-paste-v1".into(),
                        name: "Annotated reference paste".into(),
                        media_type: "application/json".into(),
                        bytes: receipt,
                    },
                );
                WorkspaceSession::snapshot(
                    &staged_docs,
                    index,
                    candidate.get_untracked().as_ref(),
                    &stored_candidates.get_untracked(),
                    &assets,
                    &workspace_name.get_untracked(),
                )?
                .validate()
                .map_err(|e| e.to_string())?;
                batch(|| {
                    documents.set(staged_docs);
                    attachments.set(assets);
                    selected_feature.set(None);
                    selection.set(None);
                    search_regions.set(None);
                    pending_reference_paste.set(None);
                    modal.set("");
                    message.set("Annotated reference pasted in one reversible edit. Original source and destination events remain retained; this does not verify a physical cloning reaction.".into());
                });
                Ok(())
            },
        )
    }
}
