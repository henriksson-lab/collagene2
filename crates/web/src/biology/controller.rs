use collagene_core::documents::history::transition::{AnalysisOrigin, AnalysisProducer};
use collagene_core::documents::transfer;
use collagene_core::material::lineage::EventBinding;
use collagene_model::{Location, Molecule};
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
    pub(crate) fn antibody_publication(
        &self,
    ) -> (
        Memo<Option<std::sync::Arc<crate::biology::antibody::antibody_overlay::Overlay>>>,
        Callback<
            Option<std::sync::Arc<crate::biology::antibody::antibody_overlay::Overlay>>,
            Result<(), String>,
        >,
    ) {
        let WorkspaceSession {
            documents,
            active,
            workspace_restoring,
            workspace_epoch,
            ..
        } = self.session;
        let published_antibody = RwSignal::new(
            None::<std::sync::Arc<crate::biology::antibody::antibody_overlay::Overlay>>,
        );
        let antibody_overlay = Memo::new(move |_| {
            if workspace_restoring.get() {
                return None;
            }
            let p = published_antibody.get()?;
            documents.with(|docs| {
                let doc = docs.get(active.get())?;
                p.capture
                    .guard(
                        workspace_epoch.get(),
                        &current_event_binding(doc),
                        &p.capture.options,
                        &p.capture.limits,
                        &p.capture.selection,
                    )
                    .ok()?;
                Some(p)
            })
        });
        let publish_antibody = Callback::new(
            move |value: Option<
                std::sync::Arc<crate::biology::antibody::antibody_overlay::Overlay>,
            >|
                  -> Result<(), String> {
                if let Some(p) = &value {
                    if workspace_restoring.get_untracked() {
                        return Err("Wait for the workspace update to finish.".into());
                    }
                    documents.with_untracked(|docs| {
                        p.capture.guard(
                            workspace_epoch.get_untracked(),
                            &current_event_binding(&docs[active.get_untracked()]),
                            &p.capture.options,
                            &p.capture.limits,
                            &p.capture.selection,
                        )
                    })?;
                }
                published_antibody.set(value);
                Ok(())
            },
        );
        (antibody_overlay, publish_antibody)
    }
    pub(crate) fn prepare_antibody(&self) -> Callback<EventBinding, Result<Molecule, String>> {
        let WorkspaceSession {
            documents,
            active,
            workspace_restoring,
            ..
        } = self.session;
        Callback::new(move |binding: EventBinding| -> Result<Molecule, String> {
            if workspace_restoring.get_untracked() {
                return Err("Wait for the workspace update to finish.".into());
            }
            documents.with_untracked(|docs| {
                let doc = &docs[active.get_untracked()];
                if current_event_binding(doc) != binding {
                    return Err("The active protein source event changed.".into());
                }
                collagene_core::antibody_numbering::source_revision(&doc.molecule)
                    .map_err(|e| e.to_string())?;
                Ok(doc.molecule.clone())
            })
        })
    }
    pub(crate) fn select_antibody_residue(
        &self,
        biology_directed: RwSignal<Option<crate::biology::biology_selection::DirectedSelection>>,
    ) -> Callback<
        (
            crate::biology::antibody::antibody_model::Capture,
            usize,
            usize,
        ),
        Result<(), String>,
    > {
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
            ..
        } = self.viewer;
        Callback::new(
            move |(capture, start, end): (
                crate::biology::antibody::antibody_model::Capture,
                usize,
                usize,
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
                        &capture.limits,
                        &capture.selection,
                    )?;
                    if start >= end || end > doc.molecule.sequence.len() {
                        return Err("Invalid source residue range.".into());
                    }
                    Ok::<_, String>(())
                })?;
                batch(|| {
                    selected_feature.set(None);
                    search_regions.set(None);
                    biology_directed.set(None);
                    selection.set(Some((start, end)));
                });
                Ok(())
            },
        )
    }
    pub(crate) fn prepare_protein_liabilities(
        &self,
    ) -> Callback<
        (EventBinding, collagene_core::protein_liabilities::Options),
        Result<Molecule, String>,
    > {
        let WorkspaceSession {
            documents,
            active,
            workspace_restoring,
            ..
        } = self.session;
        Callback::new(
            move |(binding, options): (
                EventBinding,
                collagene_core::protein_liabilities::Options,
            )|
                  -> Result<Molecule, String> {
                if workspace_restoring.get_untracked() {
                    return Err("Wait for the workspace update to finish.".into());
                }
                documents.with_untracked(|docs| {
                    let doc = &docs[active.get_untracked()];
                    if current_event_binding(doc) != binding {
                        return Err("The active protein source event changed.".into());
                    }
                    collagene_core::protein_liabilities::request_revision(&doc.molecule, &options)
                        .map_err(|e| e.to_string())?;
                    Ok(doc.molecule.clone())
                })
            },
        )
    }
    pub(crate) fn saved_rna(
        &self,
        candidate: RwSignal<Option<transfer::Candidate>>,
        stored_candidates: RwSignal<Vec<transfer::Candidate>>,
    ) -> Memo<Result<Vec<collagene_core::documents::rna_history::Summary>, String>> {
        let WorkspaceSession {
            documents,
            active,
            attachments,
            workspace_name,
            ..
        } = self.session;
        Memo::new(move |_| {
            WorkspaceSession::snapshot(
                &documents.get(),
                active.get_untracked(),
                candidate.get_untracked().as_ref(),
                &stored_candidates.get_untracked(),
                &attachments.get(),
                &workspace_name.get_untracked(),
            )
            .and_then(|w| {
                collagene_core::documents::rna_history::list(&w).map_err(|e| e.to_string())
            })
        })
    }
    pub(crate) fn resolve_rna(&self) -> Callback<EventBinding, Result<Molecule, String>> {
        let WorkspaceSession {
            documents,
            workspace_restoring,
            ..
        } = self.session;
        Callback::new(move |binding: EventBinding| -> Result<Molecule, String> {
            if workspace_restoring.get_untracked() {
                return Err("Wait for the workspace update to finish.".into());
            }
            documents.with_untracked(|docs| {
                let d = docs
                    .iter()
                    .find(|d| d.id == binding.document_id)
                    .ok_or("The pinned RNA document is unavailable.")?;
                let m = d
                    .revision_molecule_ref(&binding.event_id)
                    .map_err(|e| e.to_string())?;
                if collagene_core::rna_chemistry::source_revision(m).map_err(|e| e.to_string())?
                    != binding.revision
                {
                    return Err("The pinned RNA source revision differs.".into());
                }
                Ok(m.clone())
            })
        })
    }
    pub(crate) fn save_rna(
        &self,
        contig_snapshot: Callback<
            (),
            Result<collagene_core::documents::workspace::Workspace, String>,
        >,
    ) -> Callback<(u64, collagene_core::documents::rna_history::Record), Result<(), String>> {
        let WorkspaceSession {
            workspace_restoring,
            workspace_epoch,
            attachments,
            ..
        } = self.session;
        Callback::new(
        move |(epoch, record): (u64, collagene_core::documents::rna_history::Record)| -> Result<(), String> {
            if workspace_restoring.get_untracked() || workspace_epoch.get_untracked() != epoch {
                return Err("The workspace changed; calculate RNA chemistry again.".into());
            }
            let (next, _) =
                collagene_core::documents::rna_history::capture(&contig_snapshot.run(())?, &record)
                    .map_err(|e| e.to_string())?;
            attachments.set(next.attachments);
            Ok(())
        },
    )
    }
    pub(crate) fn load_rna(
        &self,
        contig_snapshot: Callback<
            (),
            Result<collagene_core::documents::workspace::Workspace, String>,
        >,
    ) -> Callback<String, Result<collagene_core::documents::rna_history::Record, String>> {
        let WorkspaceSession {
            workspace_restoring,
            ..
        } = self.session;
        Callback::new(
            move |hash: String| -> Result<collagene_core::documents::rna_history::Record, String> {
                if workspace_restoring.get_untracked() {
                    return Err("Wait for the workspace update to finish.".into());
                }
                collagene_core::documents::rna_history::load(&contig_snapshot.run(())?, &hash)
                    .map_err(|e| e.to_string())
            },
        )
    }
    pub(crate) fn delete_rna(
        &self,
        contig_snapshot: Callback<
            (),
            Result<collagene_core::documents::workspace::Workspace, String>,
        >,
    ) -> Callback<String, Result<(), String>> {
        let WorkspaceSession {
            workspace_restoring,
            attachments,
            ..
        } = self.session;
        Callback::new(move |hash: String| -> Result<(), String> {
            if workspace_restoring.get_untracked() {
                return Err("Wait for the workspace update to finish.".into());
            }
            let next =
                collagene_core::documents::rna_history::delete(&contig_snapshot.run(())?, &hash)
                    .map_err(|e| e.to_string())?;
            attachments.set(next.attachments);
            Ok(())
        })
    }
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
                let derived={
                    let analysis_origin = AnalysisOrigin::preview(AnalysisProducer::Translation, protein.revision());
                    source.derive_event(&output.source.event_id, id, protein, "Computed reference translation; source codon policy retained, expression and experimental outcome unverified", analysis_origin)
                }
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
                session.commit_document(
                    &output.source.revision,
                    vec![
                        collagene_core::documents::reference_edit::Operation::SetFeature {
                            index: None,
                            feature,
                        },
                    ],
                )?;
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
                let feature = documents.with_untracked(|docs| {
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
                    Ok::<_, String>(feature)
                })?;
                session.commit_document(
                    &capture.source.revision,
                    vec![
                        collagene_core::documents::reference_edit::Operation::SetFeature {
                            index: None,
                            feature,
                        },
                    ],
                )
            },
        )
    }
}
