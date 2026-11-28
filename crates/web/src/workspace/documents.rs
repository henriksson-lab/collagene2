use collagene_core::documents::transfer;
use collagene_core::documents::workspace;
use collagene_core::verification;

use leptos::prelude::*;

use crate::workspace::session::WorkspaceSession;

use crate::workspace::session::MAX_WORKSPACE_BYTES;

use crate::workspace::session::document_bytes;

use crate::viewer::sequence::selection::ViewerSelection;
use crate::workspace::shell_state::ShellState;

#[derive(Clone, Copy)]
pub(crate) struct DocumentController {
    pub(crate) session: WorkspaceSession,
    pub(crate) viewer: ViewerSelection,
    pub(crate) shell: ShellState,
}
impl DocumentController {
    pub(crate) fn close_document_tab(
        &self,
        assembly_visible: RwSignal<bool>,
        open_tabs: RwSignal<crate::workspace::document_tabs_model::Tabs>,
        show_library: RwSignal<bool>,
    ) -> Callback<String, ()> {
        let session = self.session;
        let WorkspaceSession {
            documents,
            active,
            workspace_restoring,
            ..
        } = self.session;
        let ViewerSelection {
            selection,
            selected_feature,
            search_regions,
            ..
        } = self.viewer;
        let ShellState { panel, .. } = self.shell;
        Callback::new(move |id: String| {
            if workspace_restoring.get_untracked() {
                return;
            }
            let current = open_tabs.get_untracked();
            let mut next = current.clone();
            next.close(&id);
            if next == current {
                return;
            }
            let changed_selection = next.active() != current.active();
            let next_index = next.active().and_then(|id| {
                documents.with_untracked(|docs| docs.iter().position(|doc| doc.id == id))
            });
            batch(|| {
                open_tabs.set(next);
                if changed_selection {
                    if let Some(index) = next_index {
                        active.set(index);
                    } else {
                        // Keep the retained archive document intact, but expose no
                        // editor/source-capture controls until a tab is opened.
                        selection.set(None);
                        selected_feature.set(None);
                        search_regions.set(None);
                        panel.set("");
                        assembly_visible.set(false);
                        show_library.set(true);
                    }
                }
            });
        })
    }
    pub(crate) fn fork_revision(
        &self,
        candidate: RwSignal<Option<transfer::Candidate>>,
        stored_candidates: RwSignal<Vec<transfer::Candidate>>,
    ) -> Callback<String, Result<(), String>> {
        let session = self.session;
        let WorkspaceSession {
            documents,
            active,
            workspace_restoring,
            workspace_name,
            attachments,
            ..
        } = self.session;
        let ViewerSelection { molecule, .. } = self.viewer;
        Callback::new(move |event: String| {
            if workspace_restoring.get_untracked() {
                return Err("Wait for workspace restoration to finish.".into());
            }
            let mut docs = documents.get_untracked();
            let first = docs.len();
            if first >= 1000 {
                return Err("The workspace already contains1,000 documents.".into());
            }
            let source = docs
                .get(active.get_untracked())
                .ok_or("No active document")?;
            let molecule = source
                .revision_molecule(&event)
                .map_err(|e| e.to_string())?;
            if docs
                .iter()
                .map(|d| document_bytes(&d.molecule))
                .sum::<usize>()
                + document_bytes(&molecule)
                > MAX_WORKSPACE_BYTES
            {
                return Err("Forking this event exceeds the current workspace size limit.".into());
            }
            let prefix = format!("document-{}", &molecule.revision()[..16]);
            let mut suffix = 1;
            let id = loop {
                let id = format!("{prefix}-{suffix}");
                if docs.iter().all(|d| d.id != id) {
                    break id;
                }
                suffix += 1;
            };
            let fork = source.fork_event(&event, id).map_err(|e| e.to_string())?;
            docs.push(fork);
            let frame = WorkspaceSession::snapshot(
                &docs,
                first,
                candidate.get_untracked().as_ref(),
                &stored_candidates.get_untracked(),
                &attachments.get_untracked(),
                &workspace_name.get_untracked(),
            )?;
            frame.validate().map_err(|e| e.to_string())?;
            batch(|| {
                documents.set(docs);
                active.set(first);
            });
            Ok(())
        })
    }
}
