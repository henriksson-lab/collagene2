use collagene_core::documents::transfer;
use collagene_core::material::lineage::EventBinding;
use collagene_core::verification;
use collagene_core::{
    alignment, cloning::restriction::selection as restriction_selection, documents::workspace,
    material::lineage as run_lineage,
    protein::liability_annotation as protein_liability_annotation, protein_liabilities,
    sequence::editing,
};
use collagene_model::Location;
use leptos::prelude::*;

use crate::workspace::session::WorkspaceSession;

use crate::workspace::session::MAX_WORKSPACE_BYTES;

use crate::workspace::session::current_event_binding;

use crate::workspace::session::document_bytes;

use crate::viewer::sequence::selection::ViewerSelection;
use crate::workspace::shell_state::ShellState;

#[derive(Clone, Copy)]
pub(crate) struct BiologyController {
    pub(crate) session: WorkspaceSession,
    pub(crate) viewer: ViewerSelection,
    pub(crate) shell: ShellState,
}
impl BiologyController {
    pub(crate) fn open_translated_protein(
        &self,
        candidate: RwSignal<Option<transfer::Candidate>>,
        stored_candidates: RwSignal<Vec<transfer::Candidate>>,
    ) -> Callback<crate::biology::biology::ProteinOutput, Result<(), String>> {
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
            ..
        } = self.viewer;
        let ShellState { message, .. } = self.shell;
        Callback::new(
            move |output: crate::biology::biology::ProteinOutput| -> Result<(), String> {
                if workspace_restoring.get_untracked()
                    || workspace_epoch.get_untracked() != output.epoch
                {
                    return Err("The workspace changed; translate again.".into());
                }
                let index = active.get_untracked();
                let mut docs = documents.get_untracked();
                let source = docs
                    .get(index)
                    .ok_or("The translation source is unavailable.")?;
                if current_event_binding(source) != output.source
                    || source.molecule != output.request.source
                {
                    return Err(
                        "The exact translation source event changed; translate again.".into(),
                    );
                }
                output
                    .outputs
                    .validate_sources(&output.request, &output.name)
                    .map_err(|e| e.to_string())?;
                let mut protein = output.outputs.protein.clone();
                protein.metadata.insert(
                    "translation:source_document".into(),
                    output.source.document_id.clone(),
                );
                protein.metadata.insert(
                    "translation:source_event".into(),
                    output.source.event_id.clone(),
                );
                let location = serde_json::to_string(&output.outputs.report.effective_location)
                    .map_err(|e| e.to_string())?;
                if location.len() > 512 * 1024 {
                    return Err("Translation location evidence exceeds its output limit.".into());
                }
                protein
                    .metadata
                    .insert("translation:source_location".into(), location);
                let selection_evidence =
                    serde_json::to_string(&output.request.selection).map_err(|e| e.to_string())?;
                if selection_evidence.len() > 512 * 1024 {
                    return Err("Translation policy evidence exceeds its output limit.".into());
                }
                protein
                    .metadata
                    .insert("translation:selection".into(), selection_evidence);
                protein.metadata.insert(
                    "translation:options".into(),
                    serde_json::to_string(&output.outputs.report.translation.options)
                        .map_err(|e| e.to_string())?,
                );
                protein
                    .metadata
                    .insert("translation:output_name".into(), output.name.clone());
                protein.metadata.insert(
                    "translation:warnings".into(),
                    serde_json::to_string(&output.outputs.report.warnings)
                        .map_err(|e| e.to_string())?,
                );
                let prefix = format!("protein-{}", &protein.revision()[..16]);
                let mut suffix = 1usize;
                let id = loop {
                    let id = format!("{prefix}-{suffix}");
                    if !docs.iter().any(|d| d.id == id) {
                        break id;
                    }
                    suffix += 1;
                };
                let derived=source.derive_event(&output.source.event_id,id,protein,
            "Computed reference translation; source codon policy retained, expression and experimental outcome unverified")
            .map_err(|e|e.to_string())?;
                let new_index = docs.len();
                docs.push(derived);
                if docs.len() > 1000
                    || docs
                        .iter()
                        .map(|d| document_bytes(&d.molecule))
                        .sum::<usize>()
                        > MAX_WORKSPACE_BYTES
                {
                    return Err(
                        "The linked protein exceeds this workspace's document or size limit."
                            .into(),
                    );
                }
                let frame = WorkspaceSession::snapshot(
                    &docs,
                    new_index,
                    candidate.get_untracked().as_ref(),
                    &stored_candidates.get_untracked(),
                    &attachments.get_untracked(),
                    &workspace_name.get_untracked(),
                )?;
                frame.validate().map_err(|e| e.to_string())?;
                batch(|| {
                    documents.set(docs);
                    active.set(new_index);
                    selected_feature.set(None);
                    selection.set(None);
                    search_regions.set(None);
                });
                message.set(
                    "Opened the computed protein with its exact nucleotide source event retained."
                        .into(),
                );
                Ok(())
            },
        )
    }
    pub(crate) fn retain_translation_cds(
        &self,
    ) -> Callback<crate::biology::biology::AnnotationOutput, Result<(), String>> {
        let session = self.session;
        let WorkspaceSession {
            documents,
            active,
            workspace_restoring,
            workspace_epoch,
            ..
        } = self.session;
        let ViewerSelection { molecule, .. } = self.viewer;
        let ShellState { message, .. } = self.shell;
        Callback::new(
            move |output: crate::biology::biology::AnnotationOutput| -> Result<(), String> {
                if workspace_restoring.get_untracked()
                    || workspace_epoch.get_untracked() != output.epoch
                {
                    return Err("The workspace changed; translate again.".into());
                }
                let current = documents.with_untracked(|docs| {
                    docs.get(active.get_untracked()).map(current_event_binding)
                });
                if current.as_ref() != Some(&output.source)
                    || molecule.with_untracked(|m| m != &output.request.source)
                {
                    return Err(
                        "The exact translation source event changed; translate again.".into(),
                    );
                }
                output
                    .outputs
                    .validate_sources(&output.request, &output.name)
                    .map_err(|e| e.to_string())?;
                let feature = output
                    .outputs
                    .annotation
                    .clone()
                    .ok_or("This existing-CDS translation has no new annotation proposal.")?;
                let edited = editing::set_feature(&output.request.source, None, feature)
                    .map_err(|e| e.to_string())?;
                session.commit_document(&output.source.revision, edited)?;
                message.set("Retained a CDS proposal requiring review. Its notes retain the chosen frame policy; biological CDS initiation can differ.".into());
                Ok(())
            },
        )
    }
    pub(crate) fn open_translation_origin(
        &self,
        contig_snapshot: Callback<
            (),
            Result<collagene_core::documents::workspace::Workspace, String>,
        >,
        pending_translation_source: RwSignal<
            Option<crate::viewer::alignment::pairwise_source_view::Snapshot>,
        >,
        restriction_selection: Memo<Option<Vec<(usize, usize)>>>,
    ) -> Callback<collagene_core::material::lineage::EventBinding, Result<(), String>> {
        let session = self.session;
        let WorkspaceSession {
            documents,
            active,
            workspace_restoring,
            workspace_epoch,
            ..
        } = self.session;
        let ViewerSelection { molecule, .. } = self.viewer;
        let ShellState { modal, .. } = self.shell;
        Callback::new(
            move |source: collagene_core::material::lineage::EventBinding| -> Result<(), String> {
                if workspace_restoring.get_untracked() {
                    return Err("Wait for the workspace update to finish.".into());
                }
                let frame = contig_snapshot.run(())?;
                let protein = documents
                    .with_untracked(|docs| current_event_binding(&docs[active.get_untracked()]));
                let original =
                    crate::biology::biology_selection::origin(&frame, &protein, &source)?;
                let spans = molecule
                    .with_untracked(|m| m.metadata.get("translation:source_location").cloned())
                    .filter(|text| text.len() <= 512 * 1024)
                    .map(|text| {
                        serde_json::from_str::<Location>(&text)
                            .map_err(|e| e.to_string())
                            .and_then(|location| {
                                collagene_core::cloning::restriction::selection::footprint(
                                    &location,
                                    original.sequence.len(),
                                )
                                .map_err(|e| e.to_string())
                            })
                    })
                    .transpose()?
                    .unwrap_or_default();
                pending_translation_source.set(Some(
                    crate::viewer::alignment::pairwise_source_view::Snapshot {
                        source,
                        epoch: workspace_epoch.get_untracked(),
                        molecule: std::sync::Arc::new(original.clone()),
                        spans,
                    },
                ));
                modal.set("translation_source");
                Ok(())
            },
        )
    }
    pub(crate) fn select_protein_liability(
        &self,
        biology_directed: RwSignal<Option<crate::biology::biology_selection::DirectedSelection>>,
    ) -> Callback<
        (
            crate::biology::protein_liability_model::Capture,
            collagene_core::protein_liabilities::Hit,
        ),
        Result<(), String>,
    > {
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
            move |(capture, hit): (
                crate::biology::protein_liability_model::Capture,
                collagene_core::protein_liabilities::Hit,
            )|
                  -> Result<(), String> {
                if workspace_restoring.get_untracked() {
                    return Err("Wait for the workspace update to finish.".into());
                }
                documents.with_untracked(|docs| {
                    let doc = &docs[active.get_untracked()];
                    capture.guard(
                        workspace_epoch.get_untracked(),
                        &current_event_binding(doc),
                        &capture.options,
                    )?;
                    collagene_core::protein_liabilities::validate_hit(
                        &doc.molecule,
                        &capture.options,
                        &hit,
                    )
                    .map_err(|e| e.to_string())
                })?;
                batch(|| {
                    selected_feature.set(None);
                    search_regions.set(None);
                    biology_directed.set(None);
                    selection.set(Some((hit.start, hit.end)));
                });
                Ok(())
            },
        )
    }
    pub(crate) fn annotate_protein_liability(
        &self,
        candidate: RwSignal<Option<transfer::Candidate>>,
    ) -> Callback<
        (
            crate::biology::protein_liability_model::Capture,
            collagene_core::protein_liabilities::Hit,
        ),
        Result<(), String>,
    > {
        let session = self.session;
        let WorkspaceSession {
            documents,
            active,
            workspace_restoring,
            workspace_epoch,
            ..
        } = self.session;
        let ViewerSelection { molecule, .. } = self.viewer;
        Callback::new(
            move |(capture, hit): (
                crate::biology::protein_liability_model::Capture,
                collagene_core::protein_liabilities::Hit,
            )|
                  -> Result<(), String> {
                if workspace_restoring.get_untracked() {
                    return Err("Wait for the workspace update to finish.".into());
                }
                let edited = documents.with_untracked(|docs| {
                    let doc = &docs[active.get_untracked()];
                    capture.guard(
                        workspace_epoch.get_untracked(),
                        &current_event_binding(doc),
                        &capture.options,
                    )?;
                    let feature = collagene_core::protein::liability_annotation::candidate(
                        &doc.molecule,
                        &capture.options,
                        &hit,
                    )
                    .map_err(|e| e.to_string())?;
                    editing::set_feature(&doc.molecule, None, feature).map_err(|e| e.to_string())
                })?;
                session.commit_document(&capture.source.revision, edited)
            },
        )
    }
}
