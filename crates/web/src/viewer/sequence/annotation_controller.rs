use collagene_core::documents::transfer;
use collagene_core::verification;
use collagene_core::{documents::workspace, sequence::annotations::batch as annotation_batch};

use leptos::prelude::*;

use crate::workspace::session::WorkspaceSession;

use crate::workspace::session::MAX_WORKSPACE_BYTES;

use crate::workspace::session::current_event_binding;

use crate::workspace::session::document_bytes;

use crate::viewer::sequence::selection::ViewerSelection;
use crate::workspace::shell_state::ShellState;

#[derive(Clone, Copy)]
pub(crate) struct AnnotationController {
    pub(crate) session: WorkspaceSession,
    pub(crate) viewer: ViewerSelection,
    pub(crate) shell: ShellState,
}
impl AnnotationController {
    pub(crate) fn set_annotation_visibility(
        &self,
        contig_snapshot: Callback<
            (),
            Result<collagene_core::documents::workspace::Workspace, String>,
        >,
    ) -> Callback<(Vec<usize>, bool), Result<(), String>> {
        let session = self.session;
        let WorkspaceSession {
            documents,
            active,
            workspace_restoring,
            attachments,
            ..
        } = self.session;
        let ViewerSelection { molecule, .. } = self.viewer;
        Callback::new(
            move |(indices, hide): (Vec<usize>, bool)| -> Result<(), String> {
                if workspace_restoring.get_untracked() {
                    return Err("Wait for the workspace update to finish.".into());
                }
                let mut frame = contig_snapshot.run(())?;
                let (source, count) = documents.with_untracked(|docs| {
                    let doc = &docs[active.get_untracked()];
                    (current_event_binding(doc), doc.molecule.features.len())
                });
                let mut hidden = crate::viewer::sequence::annotation_visibility::load(
                    &frame.attachments,
                    &source,
                    count,
                )?;
                for index in indices {
                    if index >= count {
                        return Err("Annotation indices changed; select them again.".into());
                    }
                    if hide {
                        hidden.insert(index);
                    } else {
                        hidden.remove(&index);
                    }
                }
                frame.attachments = crate::viewer::sequence::annotation_visibility::store(
                    &frame.attachments,
                    crate::viewer::sequence::annotation_visibility::Record {
                        schema_version: 1,
                        source,
                        feature_count: count,
                        hidden,
                    },
                )?;
                frame.validate().map_err(|e| e.to_string())?;
                crate::viewer::sequence::annotation_visibility::validate_workspace(&frame)?;
                attachments.set(frame.attachments);
                Ok(())
            },
        )
    }
    pub(crate) fn apply_annotation_batch(
        &self,
        candidate: RwSignal<Option<transfer::Candidate>>,
        stored_candidates: RwSignal<Vec<transfer::Candidate>>,
    ) -> Callback<collagene_core::sequence::annotations::batch::Request, Result<(), String>> {
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
        let ShellState { message, .. } = self.shell;
        Callback::new(
        move |request: collagene_core::sequence::annotations::batch::Request| -> Result<(), String> {
            if workspace_restoring.get_untracked() {
                return Err("Wait for the workspace update to finish.".into());
            }
            let index = active.get_untracked();
            let mut staged_docs = documents.get_untracked();
            let doc = staged_docs
                .get_mut(index)
                .ok_or("Active annotation document is unavailable.")?;
            let previous = current_event_binding(doc);
            let hidden = attachments.with_untracked(|assets| {
                crate::viewer::sequence::annotation_visibility::load(assets, &previous, doc.molecule.features.len())
            })?;
            doc.commit(vec![collagene_core::documents::reference_edit::Operation::EditAnnotations {request}], "Batch edited annotations")
                .map_err(|e| e.to_string())?;
            let source = current_event_binding(doc);
            let count = doc.molecule.features.len();
            if staged_docs
                .iter()
                .map(|d| document_bytes(&d.molecule))
                .sum::<usize>()
                > MAX_WORKSPACE_BYTES
            {
                return Err("This edit exceeds the current workspace size limit.".into());
            }
            let mut frame = WorkspaceSession::snapshot(
                &staged_docs,
                index,
                candidate.get_untracked().as_ref(),
                &stored_candidates.get_untracked(),
                &attachments.get_untracked(),
                &workspace_name.get_untracked(),
            )?;
            frame.attachments = crate::viewer::sequence::annotation_visibility::store(
                &frame.attachments,
                crate::viewer::sequence::annotation_visibility::Record {
                    schema_version: 1,
                    source,
                    feature_count: count,
                    hidden,
                },
            )?;
            frame.validate().map_err(|e| e.to_string())?;
            crate::viewer::sequence::annotation_visibility::validate_workspace(&frame)?;
            batch(|| {
                documents.set(staged_docs);
                attachments.set(frame.attachments);
            });
            message.set("Annotation changes saved as one undoable edit.".into());
            Ok(())
        },
    )
    }
}

/// Retained annotation libraries and the pending annotation-import source.
#[derive(Clone, Copy)]
pub(crate) struct AnnotationState {
    pub(crate) libraries: RwSignal<Vec<crate::viewer::sequence::feature_tools::LibraryRevision>>,
    pub(crate) annotation_source: RwSignal<String>,
}
impl AnnotationState {
    pub(crate) fn new() -> Self {
        Self {
            libraries: RwSignal::new(
                Vec::<crate::viewer::sequence::feature_tools::LibraryRevision>::new(),
            ),
            annotation_source: RwSignal::new(String::new()),
        }
    }
    pub(crate) fn restore(
        &self,
        libraries: Vec<crate::viewer::sequence::feature_tools::LibraryRevision>,
    ) {
        self.libraries.set(libraries);
        self.annotation_source.set(String::new());
    }
}
