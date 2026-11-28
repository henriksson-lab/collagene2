pub(crate) mod genome;
pub(crate) mod genome_store;
pub(crate) mod presentation;
pub mod references;
#[cfg(test)]
mod tests;

use crate::server::response::bad;
use crate::server::response::bounded_json;
use crate::server::response::bounded_json_limit;
use crate::server::response::output;
use crate::session::presentation::alignment_summary;
use crate::session::presentation::summary;
use crate::session::references::DocumentRef;
use crate::session::references::EvidenceRef;
use crate::session::references::MaterialDerivationOrigin;
use crate::session::references::MaterialOriginWorkspace;
use crate::session::references::WorkspaceRef;
use collagene_core::documents::transfer;
use collagene_core::documents::workspace;
use collagene_core::formats::alignment as alignment_formats;
use collagene_core::formats::annotations as annotation_formats;
use collagene_core::formats::oligos as oligo_formats;
use collagene_core::planning::plan::DesignPlan;
use collagene_core::planning::plan::Operation;
use collagene_core::planning::plan::PlanInput;
use collagene_core::sequence::annotations::library as feature_library;
use collagene_model::Molecule;
use rmcp::ErrorData;
use rmcp::model::CallToolResult;
use serde::Serialize;
use serde_json::Value;
use serde_json::json;
use std::collections::BTreeMap;
use std::collections::BTreeSet;
use std::sync::Arc;
use std::sync::Mutex;

pub(crate) const MAX_STORE_BYTES: usize = 64 * 1024 * 1024;
pub(crate) const MAX_ITEMS: usize = 128;
pub(crate) const MAX_DOCUMENT_REVISIONS: usize = 10_000;
pub(crate) const MAX_UPLOADS: usize = 4;
pub(crate) const MAX_UPLOAD_CHUNK: usize = 256 * 1024;
#[derive(Clone)]
pub(crate) struct StoredPlan {
    pub(crate) plan: DesignPlan,
    pub(crate) target: Molecule,
}
pub(crate) enum ArtifactData {
    Text(String),
    Binary(Vec<u8>),
}
pub(crate) struct Artifact {
    pub(crate) data: ArtifactData,
    pub(crate) mime_type: &'static str,
    pub(crate) encoding: &'static str,
    pub(crate) filename: &'static str,
}
#[derive(Clone, Serialize)]
pub(crate) struct StoredCandidate {
    pub(crate) plan_revision: String,
    pub(crate) name: String,
    pub(crate) target_revision: String,
    pub(crate) inputs: Vec<(String, String)>,
    pub(crate) steps: Vec<Operation>,
}
#[derive(Clone, Serialize)]
pub(crate) struct StoredWorkspace {
    pub(crate) name: String,
    pub(crate) active_document: String,
    pub(crate) documents: Vec<workspace::WorkspaceDocument>,
    pub(crate) revision_ids: Vec<String>,
    pub(crate) history: Vec<workspace::RevisionEvent>,
    pub(crate) candidates: Vec<StoredCandidate>,
    pub(crate) active_candidate: Option<String>,
    pub(crate) attachments: BTreeMap<String, workspace::WorkspaceAttachment>,
}
impl StoredWorkspace {
    pub(crate) fn split(value: workspace::Workspace) -> (Self, BTreeMap<String, Molecule>) {
        let revisions = value.revisions;
        (
            Self {
                name: value.name,
                active_document: value.active_document,
                documents: value.documents,
                revision_ids: revisions.keys().cloned().collect(),
                history: value.history,
                candidates: value
                    .candidates
                    .into_iter()
                    .map(|c| StoredCandidate {
                        plan_revision: c.plan_revision,
                        name: c.plan.name,
                        target_revision: c.plan.target_revision,
                        inputs: c
                            .plan
                            .inputs
                            .into_iter()
                            .map(|i| (i.id, i.revision))
                            .collect(),
                        steps: c.plan.steps,
                    })
                    .collect(),
                active_candidate: value.active_candidate,
                attachments: value.attachments,
            },
            revisions,
        )
    }
    pub(crate) fn restore(&self, store: &Store) -> Result<workspace::Workspace, ErrorData> {
        let molecule = |revision: &String| {
            store
                .documents
                .get(revision)
                .cloned()
                .ok_or_else(|| bad("workspace revision unavailable"))
        };
        Ok(workspace::Workspace {
            name: self.name.clone(),
            active_document: self.active_document.clone(),
            documents: self.documents.clone(),
            revisions: self
                .revision_ids
                .iter()
                .map(|r| Ok((r.clone(), molecule(r)?)))
                .collect::<Result<_, ErrorData>>()?,
            history: self.history.clone(),
            active_candidate: self.active_candidate.clone(),
            attachments: self.attachments.clone(),
            candidates: self
                .candidates
                .iter()
                .map(|c| {
                    Ok(transfer::Candidate {
                        plan_revision: c.plan_revision.clone(),
                        plan: DesignPlan {
                            name: c.name.clone(),
                            target_revision: c.target_revision.clone(),
                            steps: c.steps.clone(),
                            inputs: c
                                .inputs
                                .iter()
                                .map(|(id, r)| {
                                    Ok(PlanInput {
                                        id: id.clone(),
                                        revision: r.clone(),
                                        molecule: molecule(r)?,
                                    })
                                })
                                .collect::<Result<_, ErrorData>>()?,
                        },
                    })
                })
                .collect::<Result<_, ErrorData>>()?,
        })
    }
    pub(crate) fn summary(&self, id: &str) -> Value {
        json!({"id":id,"revision":id,"name":self.name,"active_document":self.active_document,
            "document_count":self.documents.len(),"revision_count":self.revision_ids.len(),"history_event_count":self.history.len(),
            "attachment_count":self.attachments.len(),"active_candidate":self.active_candidate,
            "candidate_count":self.candidates.len(),"candidate_summaries_truncated":self.candidates.len()>20,"candidates":self.candidates.iter().take(20).map(|c|json!({"plan_revision":c.plan_revision,"name":c.name.chars().take(128).collect::<String>(),"target_revision":c.target_revision,"state":"candidate","verified":false})).collect::<Vec<_>>(),
            "authority":"Archived candidates are preserved data. To compute a plan, explicitly rebind its inputs and independent target with propose_plan; import does not assert verification."})
    }
}
#[derive(Clone, Serialize)]
pub(crate) struct StoredSelection {
    pub(crate) collection_id: String,
    pub(crate) source_rows: Vec<usize>,
}
pub(crate) struct WorkspaceUpload {
    pub(crate) expected: usize,
    pub(crate) hash: String,
    pub(crate) bytes: Vec<u8>,
}
#[derive(Default)]
pub(crate) struct Store {
    pub(crate) documents: BTreeMap<String, Molecule>,
    pub(crate) plans: BTreeMap<String, StoredPlan>,
    pub(crate) artifacts: BTreeMap<String, Artifact>,
    pub(crate) used: usize,
    pub(crate) workspaces: BTreeMap<String, StoredWorkspace>,
    pub(crate) oligos: BTreeMap<String, oligo_formats::OligoTable>,
    pub(crate) uploads: BTreeMap<String, WorkspaceUpload>,
    pub(crate) oligo_selections: BTreeMap<String, StoredSelection>,
    pub(crate) upload_serial: u64,
    pub(crate) feature_libraries: BTreeMap<String, feature_library::FeatureLibrary>,
    pub(crate) feature_scans: BTreeMap<String, feature_library::ScanReport>,
    pub(crate) annotation_datasets: BTreeMap<String, annotation_formats::AnnotationImportReport>,
    pub(crate) alignments: BTreeMap<String, Arc<alignment_formats::AlignmentDocument>>,
    pub(crate) read_alignments:
        BTreeMap<String, Arc<collagene_core::review::reads::store::Collection>>,
    pub(crate) variants:
        BTreeMap<String, Arc<collagene_core::review::variants::reference::exchange::Collection>>,
    pub(crate) contig_reports: BTreeMap<String, String>,
    pub(crate) sbol_sources: BTreeSet<String>,
}
impl Store {
    pub(crate) fn variants(
        &self,
        reference: &EvidenceRef,
    ) -> Result<Arc<collagene_core::review::variants::reference::exchange::Collection>, ErrorData>
    {
        if reference.id != reference.revision {
            return Err(bad("stale variant collection revision"));
        }
        self.variants
            .get(&reference.id)
            .cloned()
            .ok_or_else(|| bad("unknown variant collection"))
    }
    pub(crate) fn read_alignments(
        &self,
        reference: &EvidenceRef,
    ) -> Result<Arc<collagene_core::review::reads::store::Collection>, ErrorData> {
        if reference.id != reference.revision {
            return Err(bad("stale read alignment collection revision"));
        }
        self.read_alignments
            .get(&reference.id)
            .cloned()
            .ok_or_else(|| bad("unknown read alignment collection"))
    }
    pub(crate) fn alignment(
        &self,
        reference: &EvidenceRef,
    ) -> Result<Arc<alignment_formats::AlignmentDocument>, ErrorData> {
        if reference.id != reference.revision {
            return Err(bad("stale alignment revision"));
        }
        self.alignments
            .get(&reference.id)
            .cloned()
            .ok_or_else(|| bad("unknown alignment"))
    }
    pub(crate) fn document(&self, reference: &DocumentRef) -> Result<&Molecule, ErrorData> {
        let mol = self
            .documents
            .get(&reference.id)
            .ok_or_else(|| bad("unknown document ID"))?;
        if reference.id != reference.revision {
            return Err(bad("stale document revision"));
        }
        Ok(mol)
    }
    pub(crate) fn feature_library(
        &self,
        reference: &EvidenceRef,
    ) -> Result<&feature_library::FeatureLibrary, ErrorData> {
        if reference.id != reference.revision {
            return Err(bad("stale feature library revision"));
        }
        self.feature_libraries
            .get(&reference.id)
            .ok_or_else(|| bad("unknown feature library"))
    }
    pub(crate) fn feature_scan(
        &self,
        reference: &EvidenceRef,
    ) -> Result<&feature_library::ScanReport, ErrorData> {
        if reference.id != reference.revision {
            return Err(bad("stale feature scan revision"));
        }
        self.feature_scans
            .get(&reference.id)
            .ok_or_else(|| bad("unknown feature scan"))
    }
    pub(crate) fn annotation_dataset(
        &self,
        reference: &EvidenceRef,
    ) -> Result<&annotation_formats::AnnotationImportReport, ErrorData> {
        if reference.id != reference.revision {
            return Err(bad("stale annotation dataset revision"));
        }
        self.annotation_datasets
            .get(&reference.id)
            .ok_or_else(|| bad("unknown annotation dataset"))
    }
    pub(crate) fn reserve(&mut self, bytes: usize) -> Result<(), ErrorData> {
        if bytes > MAX_STORE_BYTES.saturating_sub(self.used) {
            return Err(bad(
                "session storage limit reached; restart or use fewer documents/plans",
            ));
        }
        self.used += bytes;
        Ok(())
    }
}
#[derive(Clone, Default)]
pub(crate) struct Session {
    pub(crate) store: Arc<Mutex<Store>>,
}

impl Session {
    pub(crate) fn text_artifact(
        &self,
        text: String,
        mime_type: &'static str,
        filename: &'static str,
    ) -> Result<String, ErrorData> {
        let id = blake3::hash(text.as_bytes()).to_hex().to_string();
        let mut store = self.store.lock().unwrap();
        if !store.artifacts.contains_key(&id) {
            if store.artifacts.len() >= MAX_ITEMS {
                return Err(bad("at most 128 artifacts"));
            }
            store.reserve(text.len())?;
            store.artifacts.insert(
                id.clone(),
                Artifact {
                    data: ArtifactData::Text(text),
                    mime_type,
                    encoding: "utf-8",
                    filename,
                },
            );
        }
        Ok(id)
    }
    // Install an immutable codon report and optional candidate document together.
    // Capacity checks precede every mutation, including artifact reservation.
    pub(crate) fn codon_artifact<T: Serialize>(
        &self,
        report: &T,
        candidate: Option<Molecule>,
    ) -> Result<(String, Option<Value>), ErrorData> {
        let text = bounded_json(report)?;
        let artifact_id = blake3::hash(text.as_bytes()).to_hex().to_string();
        let candidate = candidate.map(|m| (m.revision(), m));
        let mut store = self.store.lock().unwrap();
        let new_artifact = !store.artifacts.contains_key(&artifact_id);
        if new_artifact && store.artifacts.len() >= MAX_ITEMS {
            return Err(bad("at most 128 artifacts"));
        }
        let mut bytes = if new_artifact { text.len() } else { 0 };
        let candidate_summary = if let Some((id, m)) = &candidate {
            m.validate().map_err(bad)?;
            if !store.documents.contains_key(id) {
                if store.documents.len() >= MAX_DOCUMENT_REVISIONS {
                    return Err(bad("at most 10000 session document revisions"));
                }
                bytes = bytes
                    .checked_add(bounded_json(m)?.len())
                    .ok_or_else(|| bad("session capacity"))?;
            }
            Some(summary(id, m))
        } else {
            None
        };
        store.reserve(bytes)?;
        if new_artifact {
            store.artifacts.insert(
                artifact_id.clone(),
                Artifact {
                    data: ArtifactData::Text(text),
                    mime_type: "application/json",
                    encoding: "utf-8",
                    filename: "codon-design.json",
                },
            );
        }
        if let Some((id, m)) = candidate {
            store.documents.entry(id).or_insert(m);
        }
        Ok((artifact_id, candidate_summary))
    }
    // Retain a fresh assembly report and every source-derived consensus together.
    // Normal capacity errors cannot leave an orphan document or report behind.
    pub(crate) fn contig_artifacts<T: Serialize>(
        &self,
        report: &T,
        candidates: Vec<Molecule>,
        batch_revision: &str,
    ) -> Result<(String, Vec<Value>), ErrorData> {
        let text = bounded_json(report)?;
        let artifact_id = blake3::hash(text.as_bytes()).to_hex().to_string();
        let mut rows = Vec::with_capacity(candidates.len());
        let mut unique = BTreeMap::new();
        for molecule in candidates {
            molecule.validate().map_err(bad)?;
            let encoded_size = bounded_json(&molecule)?.len();
            let id = molecule.revision();
            rows.push(summary(&id, &molecule));
            unique.entry(id).or_insert((molecule, encoded_size));
        }
        let mut store = self.store.lock().unwrap();
        let new_artifact = match store.artifacts.get(&artifact_id) {
            Some(existing) => {
                if !matches!(existing.data, ArtifactData::Text(_))
                    || existing.mime_type != "application/json"
                    || existing.filename != "contig-assembly.json"
                {
                    return Err(bad(
                        "assembly artifact already has a different representation",
                    ));
                }
                false
            }
            None => true,
        };
        if new_artifact && store.artifacts.len() >= MAX_ITEMS {
            return Err(bad("at most 128 artifacts"));
        }
        let new_documents = unique
            .keys()
            .filter(|id| !store.documents.contains_key(*id))
            .count();
        if store.documents.len().saturating_add(new_documents) > MAX_DOCUMENT_REVISIONS {
            return Err(bad("at most 10000 session document revisions"));
        }
        let new_cache = !store.contig_reports.contains_key(batch_revision);
        let mut bytes = if new_artifact { text.len() } else { 0 };
        if new_cache {
            bytes = bytes
                .checked_add(batch_revision.len() + artifact_id.len() + 64)
                .ok_or_else(|| bad("session capacity"))?;
        }
        for (id, (_, size)) in &unique {
            if !store.documents.contains_key(id) {
                bytes = bytes
                    .checked_add(*size)
                    .ok_or_else(|| bad("session capacity"))?;
            }
        }
        store.reserve(bytes)?;
        if new_artifact {
            store.artifacts.insert(
                artifact_id.clone(),
                Artifact {
                    data: ArtifactData::Text(text),
                    mime_type: "application/json",
                    encoding: "utf-8",
                    filename: "contig-assembly.json",
                },
            );
        }
        for (id, (molecule, _)) in unique {
            store.documents.entry(id).or_insert(molecule);
        }
        store
            .contig_reports
            .insert(batch_revision.to_owned(), artifact_id.clone());
        Ok((artifact_id, rows))
    }
    // Candidate donor installation is atomic with its source-bound computed report.
    pub(crate) fn edit_template_artifact<T: Serialize>(
        &self,
        report: &T,
        donor: Option<Molecule>,
    ) -> Result<(String, Option<Value>), ErrorData> {
        self.candidate_document_artifact(report, donor, "edit-template.json")
    }
    pub(crate) fn candidate_document_artifact<T: Serialize>(
        &self,
        report: &T,
        donor: Option<Molecule>,
        filename: &'static str,
    ) -> Result<(String, Option<Value>), ErrorData> {
        self.candidate_document_artifact_bounded(report, donor, filename, MAX_STORE_BYTES)
    }
    pub(crate) fn candidate_document_artifact_bounded<T: Serialize>(
        &self,
        report: &T,
        donor: Option<Molecule>,
        filename: &'static str,
        max_report_bytes: usize,
    ) -> Result<(String, Option<Value>), ErrorData> {
        let text = bounded_json_limit(report, max_report_bytes)?;
        let artifact_id = blake3::hash(text.as_bytes()).to_hex().to_string();
        let donor = donor
            .map(|m| {
                m.validate().map_err(bad)?;
                let bytes = bounded_json(&m)?.len();
                Ok::<_, ErrorData>((m.revision(), m, bytes))
            })
            .transpose()?;
        let mut store = self.store.lock().unwrap();
        let new_artifact = match store.artifacts.get(&artifact_id) {
            Some(a) => {
                if !matches!(a.data, ArtifactData::Text(_))
                    || a.mime_type != "application/json"
                    || a.filename != filename
                {
                    return Err(bad("candidate artifact already has another representation"));
                }
                false
            }
            None => true,
        };
        if new_artifact && store.artifacts.len() >= MAX_ITEMS {
            return Err(bad("at most 128 artifacts"));
        }
        let mut bytes = if new_artifact { text.len() } else { 0 };
        let donor_summary = if let Some((id, m, size)) = &donor {
            if !store.documents.contains_key(id) {
                if store.documents.len() >= MAX_DOCUMENT_REVISIONS {
                    return Err(bad("at most 10000 session document revisions"));
                }
                bytes = bytes
                    .checked_add(*size)
                    .ok_or_else(|| bad("session capacity"))?;
            }
            Some(summary(id, m))
        } else {
            None
        };
        store.reserve(bytes)?;
        if new_artifact {
            store.artifacts.insert(
                artifact_id.clone(),
                Artifact {
                    data: ArtifactData::Text(text),
                    mime_type: "application/json",
                    encoding: "utf-8",
                    filename,
                },
            );
        }
        if let Some((id, m, _)) = donor {
            store.documents.entry(id).or_insert(m);
        }
        Ok((artifact_id, donor_summary))
    }
    // Full dot-plot JSON and SVG are installed together after all capacity checks.
    pub(crate) fn dotplot_artifacts(
        &self,
        report: &collagene_core::alignment::dotplot::DotPlotReport,
        svg: String,
    ) -> Result<(String, String), ErrorData> {
        let text = bounded_json(report)?;
        let entries = [
            (text, "application/json", "dotplot.json"),
            (svg, "image/svg+xml", "dotplot.svg"),
        ];
        self.paired_text_artifacts(entries)
    }
    // Publish the documentary receipt, binary candidate and unverified session
    // plan atomically. Neither workspace heads nor historical status are adopted.
    pub(crate) fn saved_review_artifacts(
        &self,
        review: collagene_core::material::review::PreparedReview,
        reference: WorkspaceRef,
    ) -> Result<CallToolResult, ErrorData> {
        let text = bounded_json_limit(&review, collagene_core::material::review::MAX_REVIEW_BYTES)?;
        let bytes = transfer::save_candidate(&review.candidate).map_err(bad)?;
        let report_id = blake3::hash(text.as_bytes()).to_hex().to_string();
        let candidate_id = blake3::hash(&bytes).to_hex().to_string();
        let plan_revision = review.candidate.plan_revision.clone();
        let plan_bytes =
            bounded_json_limit(&(&review.candidate.plan, &review.target), MAX_STORE_BYTES)?.len();
        let inputs = review.inputs.iter().take(100).collect::<Vec<_>>();
        let assessment = &review.recorded_assessment;
        let recorded_assessment = json!({"engine_version":assessment.engine_version,"status":assessment.status,"findings":assessment.findings.iter().take(20).map(|s|s.chars().take(1024).collect::<String>()).collect::<Vec<_>>(),"finding_count":assessment.findings.len(),"findings_truncated":assessment.findings.len()>20||assessment.findings.iter().take(20).any(|s|s.chars().count()>1024),"assumptions":assessment.assumptions.iter().take(20).map(|s|s.chars().take(1024).collect::<String>()).collect::<Vec<_>>(),"assumption_count":assessment.assumptions.len(),"assumptions_truncated":assessment.assumptions.len()>20||assessment.assumptions.iter().take(20).any(|s|s.chars().count()>1024),"exhaustive_alternative_search":assessment.exhaustive_alternative_search,"authority":"historical only"});
        let result = output(
            json!({"workspace":reference,"run_id":review.run_id,"context_revision":review.context_revision,"source_workspace_fingerprint":review.source_workspace_fingerprint,"run_attachment_hash":review.run_attachment_hash,"state":"candidate","verified":false,"plan_revision":plan_revision,"target_event":review.target_event,"target":{"id":review.target.revision(),"revision":review.target.revision()},"input_count":review.inputs.len(),"inputs":inputs,"inputs_truncated":review.inputs.len()>100,"recorded_assessment":recorded_assessment,"recorded_engine_matches":review.recorded_engine_matches,"current_verification_engine":review.current_verification_engine,"has_saved_material_origins":review.has_saved_material_origins,"warnings":review.warnings,"artifact_id":report_id,"candidate_artifact_id":candidate_id,"candidate_filename":"saved-run-review.cg2plan","candidate_encoding":"base64","authority":"Documentary preparation of exact retained run inputs and independent target. The session candidate is unverified; call verify_plan for fresh checks. No historical status/material authority was adopted and no workspace head or native catalog was changed."}),
        )?;
        let mut store = self.store.lock().unwrap();
        let new_plan = !store.plans.contains_key(&plan_revision);
        if let Some(existing) = store.plans.get(&plan_revision) {
            if existing.target != review.target
                || bounded_json(&existing.plan)? != bounded_json(&review.candidate.plan)?
            {
                return Err(bad(
                    "saved plan revision already has different exact source/target data",
                ));
            }
        }
        let representations = [
            (
                &report_id,
                false,
                "application/json",
                "saved-run-review.json",
            ),
            (
                &candidate_id,
                true,
                "application/octet-stream",
                "saved-run-review.cg2plan",
            ),
        ];
        let mut added = 0usize;
        let mut charge = if new_plan { plan_bytes } else { 0 };
        for (id, binary, mime, filename) in representations {
            if let Some(existing) = store.artifacts.get(id) {
                if matches!(&existing.data, ArtifactData::Binary(_)) != binary
                    || existing.mime_type != mime
                    || existing.filename != filename
                {
                    return Err(bad("saved review artifact has another representation"));
                }
            } else {
                added += 1;
                charge = charge
                    .checked_add(if binary { bytes.len() } else { text.len() })
                    .ok_or_else(|| bad("session capacity"))?;
            }
        }
        if (new_plan && store.plans.len() >= MAX_ITEMS)
            || store.artifacts.len().saturating_add(added) > MAX_ITEMS
        {
            return Err(bad(
                "session item limit; no saved-review artifacts or plan were published",
            ));
        }
        store.reserve(charge)?;
        store.artifacts.entry(report_id).or_insert(Artifact {
            data: ArtifactData::Text(text),
            mime_type: "application/json",
            encoding: "utf-8",
            filename: "saved-run-review.json",
        });
        store.artifacts.entry(candidate_id).or_insert(Artifact {
            data: ArtifactData::Binary(bytes),
            mime_type: "application/octet-stream",
            encoding: "base64",
            filename: "saved-run-review.cg2plan",
        });
        if new_plan {
            store.plans.insert(
                plan_revision,
                StoredPlan {
                    plan: review.candidate.plan,
                    target: review.target,
                },
            );
        }
        Ok(result)
    }
    pub(crate) fn paired_text_artifacts(
        &self,
        entries: [(String, &'static str, &'static str); 2],
    ) -> Result<(String, String), ErrorData> {
        let ids = entries
            .iter()
            .map(|(text, _, _)| blake3::hash(text.as_bytes()).to_hex().to_string())
            .collect::<Vec<_>>();
        let mut store = self.store.lock().unwrap();
        let mut bytes = 0usize;
        let mut added = 0usize;
        for ((text, mime_type, filename), id) in entries.iter().zip(&ids) {
            if let Some(existing) = store.artifacts.get(id) {
                if !matches!(existing.data, ArtifactData::Text(_))
                    || existing.mime_type != *mime_type
                    || existing.filename != *filename
                {
                    return Err(bad(
                        "artifact content already exists with a different representation",
                    ));
                }
            } else {
                bytes = bytes
                    .checked_add(text.len())
                    .ok_or_else(|| bad("session capacity"))?;
                added += 1;
            }
        }
        if store.artifacts.len().saturating_add(added) > MAX_ITEMS {
            return Err(bad("at most 128 artifacts"));
        }
        store.reserve(bytes)?;
        for ((text, mime_type, filename), id) in entries.into_iter().zip(&ids) {
            store.artifacts.entry(id.clone()).or_insert(Artifact {
                data: ArtifactData::Text(text),
                mime_type,
                encoding: "utf-8",
                filename,
            });
        }
        Ok((ids[0].clone(), ids[1].clone()))
    }
    pub(crate) fn structure_output(
        &self,
        original: String,
        report: &collagene_core::protein::structure::Prepared,
    ) -> Result<CallToolResult, ErrorData> {
        use collagene_core::protein::structure::coordinates::Format;
        if blake3::hash(original.as_bytes()).to_hex().as_str() != report.source_hash {
            return Err(bad("structure report does not bind its original source"));
        }
        let encoded = bounded_json(report)?;
        let (mime, filename) = match report.format {
            Format::Pdb => ("chemical/x-pdb", "structure.pdb"),
            Format::Mmcif => ("chemical/x-mmcif", "structure.cif"),
        };
        let (source, artifact) = self.paired_text_artifacts([
            (original, mime, filename),
            (encoded, "application/json", "structure-view.json"),
        ])?;
        let models = report
            .atoms
            .iter()
            .map(|a| a.model)
            .collect::<BTreeSet<_>>();
        output(
            json!({"source":{"id":source,"revision":source},"artifact_id":artifact,
            "source_hash":report.source_hash,"context_revision":report.revision,"engine":report.engine,"format":report.format,
            "atom_count":report.atoms.len(),"models":models,"selection":report.scene.selection,
            "selected_atom_count":report.scene.points.len(),"declared_connection_count":report.notes.connection_declarations,
            "display_edge_count":report.scene.edges.len(),"unresolved_connections":report.scene.unresolved_connections,
            "withheld_symmetry_or_relation":report.scene.withheld_symmetry_or_relation,
            "atoms_without_element":report.notes.atoms_without_element,"atoms_with_alternates":report.notes.atoms_with_alternates,
            "atom_preview":report.scene.points.iter().take(32).map(|p|json!({"source_atom_index":p.atom_index,"atom":report.atoms[p.atom_index as usize]})).collect::<Vec<_>>(),
            "atom_preview_truncated":report.scene.points.len()>32,
            "interpretation":"Imported coordinates only; no fold prediction, sequence correspondence, experimental validity or cloning/material authority."}),
        )
    }
    pub(crate) fn analysis_artifact(&self, value: &Value) -> Result<String, ErrorData> {
        let text = serde_json::to_string(value).map_err(bad)?;
        let id = blake3::hash(text.as_bytes()).to_hex().to_string();
        let mut store = self.store.lock().unwrap();
        if !store.artifacts.contains_key(&id) {
            if store.artifacts.len() >= MAX_ITEMS {
                return Err(bad("at most 128 artifacts"));
            }
            store.reserve(text.len())?;
            store.artifacts.insert(
                id.clone(),
                Artifact {
                    data: ArtifactData::Text(text),
                    mime_type: "application/json",
                    encoding: "utf-8",
                    filename: "sequence-analysis.json",
                },
            );
        }
        Ok(id)
    }
    pub(crate) fn binary_artifact(
        &self,
        bytes: Vec<u8>,
        mime_type: &'static str,
        filename: &'static str,
    ) -> Result<String, ErrorData> {
        // Existing clients use binary artifact IDs as the decoded-content digest.
        let id = blake3::hash(&bytes).to_hex().to_string();
        let mut store = self.store.lock().unwrap();
        if let Some(existing) = store.artifacts.get(&id) {
            if !matches!(existing.data, ArtifactData::Binary(_))
                || existing.mime_type != mime_type
                || existing.filename != filename
            {
                return Err(bad(
                    "artifact content already exists with a different representation",
                ));
            }
        }
        if !store.artifacts.contains_key(&id) {
            if store.artifacts.len() >= MAX_ITEMS {
                return Err(bad("at most 128 artifacts"));
            }
            store.reserve(bytes.len())?;
            store.artifacts.insert(
                id.clone(),
                Artifact {
                    data: ArtifactData::Binary(bytes),
                    mime_type,
                    encoding: "base64",
                    filename,
                },
            );
        }
        Ok(id)
    }
    pub(crate) fn insert_alignment(
        &self,
        document: alignment_formats::AlignmentDocument,
    ) -> Result<CallToolResult, ErrorData> {
        let bytes = alignment_formats::save_alignment(&document).map_err(bad)?;
        let id = blake3::hash(&bytes).to_hex().to_string();
        let mut result = alignment_summary(&id, &document);
        result["source_is_current"] = json!(document.source_is_current().map_err(bad)?);
        let mut store = self.store.lock().unwrap();
        if !store.alignments.contains_key(&id) {
            if store.alignments.len() >= MAX_ITEMS {
                return Err(bad("at most128 alignments"));
            }
            store.reserve(bytes.len())?;
            store.alignments.insert(id, Arc::new(document));
        }
        output(result)
    }
    pub(crate) fn insert_feature_library(
        &self,
        library: feature_library::FeatureLibrary,
    ) -> Result<CallToolResult, ErrorData> {
        let bytes = feature_library::save_library(&library).map_err(bad)?;
        let id = library.revision();
        let mut store = self.store.lock().unwrap();
        let result = json!({"library":{"id":id,"revision":id},"name":library.name,"version":library.version,"entry_count":library.entries.len(),"entries":library.entries.iter().take(100).map(|entry|json!({"id":entry.id,"revision":entry.revision,"name":entry.name.chars().take(256).collect::<String>(),"length":entry.sequence.len(),"kind":entry.kind,"source":{"document_revision":entry.source.document_revision,"feature_index":entry.source.feature_index,"document_name":entry.source.document_name.chars().take(256).collect::<String>(),"source_needs_review":entry.source.source_needs_review}})).collect::<Vec<_>>(),"entries_truncated":library.entries.len()>100,"persistence":"session only"});
        if !store.feature_libraries.contains_key(&id) {
            if store.feature_libraries.len() >= MAX_ITEMS {
                return Err(bad("at most 128 feature libraries"));
            }
            store.reserve(bytes.len())?;
            store.feature_libraries.insert(id, library);
        }
        output(result)
    }
    /// Attest a documentary saved-source capsule against one exact retained
    /// workspace without synthesizing a plan or running scientific reconstruction.
    pub(crate) fn attest_run_source(
        &self,
        reference: &WorkspaceRef,
        origin: &collagene_core::material::capsule::Capsule,
        source: &Molecule,
    ) -> Result<(), ErrorData> {
        let frame = self.restore_workspace(reference)?;
        let recorded =
            collagene_core::material::lineage::load_run(&frame, &origin.selection.run_id)
                .map_err(bad)?;
        let candidate = transfer::load_candidate(&origin.candidate).map_err(bad)?;
        if candidate.plan_revision != recorded.plan_revision
            || origin.target_event != recorded.target
            || origin.inputs != recorded.inputs
            || origin.products != recorded.products
            || frame.revisions.get(&recorded.target.revision) != Some(&origin.target)
            || frame.revisions.get(&origin.selection.event.revision) != Some(source)
        {
            return Err(bad(
                "saved-run capsule differs from the exact retained run/source/target events",
            ));
        }
        let origins =
            collagene_core::material::run::originating_products(&frame, &origin.selection.event)
                .map_err(bad)?;
        if !origins.contains(&origin.selection) {
            return Err(bad(
                "saved-run origin selection is not retained for this exact source event",
            ));
        }
        Ok(())
    }
    pub(crate) fn resolve_material_derivation_origin(
        &self,
        workspace: &workspace::Workspace,
        event: &collagene_core::material::lineage::EventBinding,
        origin: &MaterialDerivationOrigin,
    ) -> Result<collagene_core::material::derivation::SourceOrigin, ErrorData> {
        use collagene_core::material::derivation::SourceOrigin;
        match origin {
            MaterialDerivationOrigin::OriginalReference => {
                SourceOrigin::OriginalReference
                    .validate_shallow_for(event)
                    .map_err(bad)?;
                Ok(SourceOrigin::OriginalReference)
            }
            MaterialDerivationOrigin::Saved(selection) => self
                .resolve_material_origin(workspace, event, selection)
                .map(SourceOrigin::Saved),
        }
    }
    pub(crate) fn resolve_material_origin(
        &self,
        workspace: &workspace::Workspace,
        event: &collagene_core::material::lineage::EventBinding,
        selection: &collagene_core::material::run::MaterialSelection,
    ) -> Result<collagene_core::material::run::MaterialOrigin, ErrorData> {
        use collagene_core::material::run as run_material;

        selection.validate().map_err(bad)?;
        if selection.event() != event {
            return Err(bad(
                "Selected material handle names a different exact source event.",
            ));
        }
        let requested = Some(selection.clone());
        let expected = run_material::material_inspection_revision(workspace, event, &requested)
            .map_err(bad)?;
        let inspection =
            run_material::inspect_material(workspace, event, &requested).map_err(bad)?;
        inspection.validate_binding(&expected).map_err(bad)?;
        let selected=inspection.selected.ok_or_else(||bad("Selected retained material cannot currently be reconstructed; inspect its origin findings."))?;
        if selected.selection() != *selection {
            return Err(bad("Fresh origin differs from explicit selected handle."));
        }
        selected
            .validate_for_workspace(workspace, event)
            .map_err(bad)?;
        Ok(selected)
    }
    pub(crate) fn attest_material_origins(
        &self,
        inputs: &[PlanInput],
        bindings: &[MaterialOriginWorkspace],
        origins: &[(
            &str,
            &str,
            collagene_core::material::capsule::MaterialOriginRef<'_>,
        )],
    ) -> Result<(), ErrorData> {
        use collagene_core::material::{
            capsule::MaterialOriginRef,
            run::{MaterialOrigin, MaterialSelection},
        };
        if bindings.len() != origins.len() || bindings.len() > 128 {
            return Err(bad(
                "Every saved reaction role requires exactly one retained workspace binding.",
            ));
        }
        let mut seen = BTreeSet::new();
        for binding in bindings {
            if binding.role.is_empty()
                || binding.role.len() > 128
                || !seen.insert(binding.role.as_str())
                || !origins.iter().any(|(role, _, _)| *role == binding.role)
            {
                return Err(bad(
                    "Invalid, duplicate or unrelated saved reaction workspace role.",
                ));
            }
        }
        for (role, input_id, origin) in origins {
            let binding = bindings
                .iter()
                .find(|b| b.role == *role)
                .ok_or_else(|| bad("Saved reaction workspace is unavailable."))?;
            let selection = match origin {
                MaterialOriginRef::Run(c) => MaterialSelection::Run(c.selection.clone()),
                MaterialOriginRef::Derived(c) => MaterialSelection::Derived(c.product.clone()),
            };
            let frame = self.restore_workspace(&binding.workspace)?;
            let actual = self.resolve_material_origin(&frame, selection.event(), &selection)?;
            let source = inputs
                .iter()
                .find(|i| i.id == *input_id)
                .ok_or_else(|| bad("Saved reaction source role is unavailable."))?;
            if frame.revisions.get(&selection.event().revision) != Some(&source.molecule) {
                return Err(bad(
                    "Saved reaction source differs from its exact retained event.",
                ));
            }
            let same = match (&actual, origin) {
                (MaterialOrigin::Run(a), MaterialOriginRef::Run(b)) => {
                    bounded_json_limit(a, 32 * 1024 * 1024)?
                        == bounded_json_limit(b, 32 * 1024 * 1024)?
                }
                (MaterialOrigin::Derived(a), MaterialOriginRef::Derived(b)) => {
                    bounded_json_limit(a, 32 * 1024 * 1024)?
                        == bounded_json_limit(b, 32 * 1024 * 1024)?
                }
                _ => false,
            };
            if !same {
                return Err(bad(
                    "Supplied reaction origin differs from the freshly resolved retained native origin.",
                ));
            }
        }
        Ok(())
    }
    pub(crate) fn material_derivation_source<'a>(
        &self,
        workspace: &'a workspace::Workspace,
        event: &collagene_core::material::lineage::EventBinding,
        origin: &collagene_core::material::derivation::SourceOrigin,
    ) -> Result<&'a Molecule, ErrorData> {
        origin.validate_shallow_for(event).map_err(bad)?;
        workspace.validate().map_err(bad)?;
        let exact = workspace
            .history
            .iter()
            .find(|e| e.id == event.event_id)
            .ok_or_else(|| bad("material source event is absent"))?;
        if exact.document_id != event.document_id || exact.revision != event.revision {
            return Err(bad("material source event/document/revision disagree"));
        }
        let document = workspace
            .documents
            .iter()
            .find(|d| d.id == event.document_id)
            .ok_or_else(|| bad("material source document is absent"))?;
        if document.head != event.revision
            || document.head_event.as_deref() != Some(event.event_id.as_str())
        {
            return Err(bad(
                "material source must be the exact retained current document head",
            ));
        }
        if let collagene_core::material::derivation::SourceOrigin::Saved(selected) = origin {
            selected
                .validate_for_workspace(workspace, event)
                .map_err(bad)?;
        }
        let source = workspace
            .revisions
            .get(&event.revision)
            .ok_or_else(|| bad("material source revision is absent"))?;
        collagene_model::bounds::model_preflight(&[source]).map_err(bad)?;
        Ok(source)
    }
    pub(crate) fn restore_workspace(
        &self,
        reference: &WorkspaceRef,
    ) -> Result<workspace::Workspace, ErrorData> {
        let store = self.store.lock().unwrap();
        if reference.id != reference.revision {
            return Err(bad("stale workspace revision"));
        }
        store
            .workspaces
            .get(&reference.id)
            .ok_or_else(|| bad("unknown workspace"))?
            .restore(&store)
    }
    pub(crate) fn import_archive_bytes(&self, bytes: &[u8]) -> Result<CallToolResult, ErrorData> {
        let value = workspace::load_workspace(bytes).map_err(bad)?;
        self.commit_workspace(value)
    }
    pub(crate) fn commit_workspace(
        &self,
        workspace: workspace::Workspace,
    ) -> Result<CallToolResult, ErrorData> {
        self.commit_workspace_artifact(workspace, None)
    }
    pub(crate) fn commit_workspace_artifact(
        &self,
        workspace: workspace::Workspace,
        run_artifact: Option<(String, Value)>,
    ) -> Result<CallToolResult, ErrorData> {
        self.commit_workspace_outputs(workspace, run_artifact, None)
    }
    // Stage every workspace/revision, run receipt and optional native output
    // before reserving capacity or publishing any session entry.
    pub(crate) fn commit_workspace_outputs(
        &self,
        workspace: workspace::Workspace,
        run_artifact: Option<(String, Value)>,
        native_artifact: Option<Vec<u8>>,
    ) -> Result<CallToolResult, ErrorData> {
        let finalized = native_artifact.is_some();
        let receipt_filename = if finalized {
            "cloning-outputs.json"
        } else {
            "cloning-run.json"
        };
        self.commit_workspace_outputs_named(
            workspace,
            run_artifact,
            native_artifact,
            receipt_filename,
            "historical run data; replay recomputes current evidence",
        )
    }
    pub(crate) fn commit_workspace_outputs_named(
        &self,
        workspace: workspace::Workspace,
        run_artifact: Option<(String, Value)>,
        native_artifact: Option<Vec<u8>>,
        receipt_filename: &'static str,
        receipt_authority: &'static str,
    ) -> Result<CallToolResult, ErrorData> {
        let mut pending_alignments = BTreeMap::new();
        for attachment in workspace
            .attachments
            .values()
            .filter(|a| a.kind == alignment_formats::WORKSPACE_ATTACHMENT_KIND)
        {
            let document = alignment_formats::load_alignment(&attachment.bytes).map_err(bad)?;
            let bytes = alignment_formats::save_alignment(&document).map_err(bad)?;
            let id = blake3::hash(&bytes).to_hex().to_string();
            pending_alignments.insert(id, (Arc::new(document), bytes.len()));
        }
        let (structure, pending) = StoredWorkspace::split(workspace);
        let metadata = serde_json::to_vec(&structure).map_err(bad)?;
        let id = blake3::hash(&metadata).to_hex().to_string();
        let mut store = self.store.lock().unwrap();
        let result = structure.summary(&id);
        let artifact_id = run_artifact
            .as_ref()
            .map(|(text, _)| blake3::hash(text.as_bytes()).to_hex().to_string());
        let native_id = native_artifact
            .as_ref()
            .map(|bytes| blake3::hash(bytes).to_hex().to_string());
        let mut added_artifacts = 0usize;
        let mut artifact_bytes = 0usize;
        if let Some((text, _)) = &run_artifact {
            let key = artifact_id
                .as_ref()
                .ok_or_else(|| bad("missing staged receipt ID"))?;
            if let Some(existing) = store.artifacts.get(key) {
                if !matches!(existing.data, ArtifactData::Text(_))
                    || existing.mime_type != "application/json"
                    || existing.filename != receipt_filename
                {
                    return Err(bad(
                        "receipt artifact already has another representation; no workspace outputs published",
                    ));
                }
            } else {
                added_artifacts += 1;
                artifact_bytes = text.len();
            }
        }
        if let Some(bytes) = &native_artifact {
            let key = native_id
                .as_ref()
                .ok_or_else(|| bad("missing staged native archive ID"))?;
            if let Some(existing) = store.artifacts.get(key) {
                if !matches!(existing.data, ArtifactData::Binary(_))
                    || existing.mime_type != "application/octet-stream"
                    || existing.filename != "workspace.cg2workspace"
                {
                    return Err(bad(
                        "native workspace artifact already has another representation; no outputs published",
                    ));
                }
            } else {
                added_artifacts += 1;
                artifact_bytes = artifact_bytes
                    .checked_add(bytes.len())
                    .ok_or_else(|| bad("artifact accounting overflow"))?;
            }
        }
        if store.artifacts.len().saturating_add(added_artifacts) > MAX_ITEMS {
            return Err(bad(
                "final workspace artifacts exceed item capacity; no outputs published",
            ));
        }
        if store.workspaces.contains_key(&id) && run_artifact.is_none() && native_artifact.is_none()
        {
            return output(
                json!({"workspace":result,"deduplicated":true,"persistence":"session only"}),
            );
        }
        let deduplicated = store.workspaces.contains_key(&id);
        if !deduplicated && store.workspaces.len() >= MAX_ITEMS {
            return Err(bad("at most 128 workspaces"));
        }
        let (count, bytes) = pending
            .iter()
            .filter(|(r, _)| !store.documents.contains_key(*r))
            .try_fold(
                (0usize, if deduplicated { 0 } else { metadata.len() }),
                |(count, bytes), (_, m)| {
                    Ok::<_, ErrorData>((
                        count + 1,
                        bytes
                            .checked_add(serde_json::to_vec(m).map_err(bad)?.len())
                            .ok_or_else(|| bad("session accounting overflow"))?,
                    ))
                },
            )?;
        if store.documents.len() + count > MAX_DOCUMENT_REVISIONS {
            return Err(bad(
                "at most 10000 session document revisions; workspace import is atomic",
            ));
        }
        let (alignment_count, alignment_bytes) = pending_alignments
            .iter()
            .filter(|(id, _)| !store.alignments.contains_key(*id))
            .try_fold((0usize, 0usize), |(count, bytes), (_, (_, size))| {
                Ok::<_, ErrorData>((
                    count + 1,
                    bytes
                        .checked_add(*size)
                        .ok_or_else(|| bad("alignment accounting overflow"))?,
                ))
            })?;
        if store.alignments.len() + alignment_count > MAX_ITEMS {
            return Err(bad(
                "workspace import exceeds128 alignment collections; no state installed",
            ));
        }
        let total = bytes
            .checked_add(alignment_bytes)
            .and_then(|n| n.checked_add(artifact_bytes))
            .ok_or_else(|| bad("session accounting overflow"))?;
        store.reserve(total)?;
        store.alignments.extend(
            pending_alignments
                .into_iter()
                .map(|(id, (document, _))| (id, document)),
        );
        store.documents.extend(pending);
        store.workspaces.insert(id, structure);
        let mut response =
            json!({"workspace":result,"deduplicated":deduplicated,"persistence":"session only"});
        if let Some((text, receipt)) = run_artifact {
            let artifact_id = artifact_id.unwrap();
            store
                .artifacts
                .entry(artifact_id.clone())
                .or_insert(Artifact {
                    data: ArtifactData::Text(text),
                    mime_type: "application/json",
                    encoding: "utf-8",
                    filename: receipt_filename,
                });
            response["receipt"] = receipt;
            response["artifact_id"] = json!(artifact_id);
            response["authority"] = json!(receipt_authority);
        }
        if let Some(bytes) = native_artifact {
            let native_id = native_id.ok_or_else(|| bad("missing native artifact ID"))?;
            let length = bytes.len();
            store
                .artifacts
                .entry(native_id.clone())
                .or_insert(Artifact {
                    data: ArtifactData::Binary(bytes),
                    mime_type: "application/octet-stream",
                    encoding: "base64",
                    filename: "workspace.cg2workspace",
                });
            response["native_artifact_id"] = json!(native_id);
            response["native_bytes"] = json!(length);
            response["encoding"] = json!("base64");
            response["filename"] = json!("workspace.cg2workspace");
            response["authority"] = json!(
                "Finalized declared-model outputs in a new immutable session workspace only. All native preparation history is retained; presentation visibility and saved order data do not establish experimental material or inherited verification authority. Replay reconstructs again."
            );
        }
        output(response)
    }
    pub(crate) fn insert_documents(
        &self,
        molecules: Vec<Molecule>,
        source_warnings: Vec<String>,
    ) -> Result<CallToolResult, ErrorData> {
        let mut pending = BTreeMap::new();
        let warnings_truncated = source_warnings.len() > 100;
        let warnings: Vec<_> = source_warnings
            .into_iter()
            .take(100)
            .map(|w| w.chars().take(1024).collect::<String>())
            .collect();
        for m in molecules {
            m.validate().map_err(bad)?;
            pending.insert(m.revision(), m);
        }
        let mut store = self.store.lock().unwrap();
        let mut bytes = 0;
        let mut count = 0;
        for (id, m) in &pending {
            if !store.documents.contains_key(id) {
                bytes += serde_json::to_vec(m).map_err(bad)?.len();
                count += 1;
            }
        }
        if store.documents.len() + count > MAX_DOCUMENT_REVISIONS {
            return Err(bad("at most 10000 session document revisions"));
        }
        store.reserve(bytes)?;
        let result: Vec<_> = pending.iter().map(|(id, m)| summary(id, m)).collect();
        store.documents.extend(pending);
        output(
            json!({"documents":result,"persistence":"session only","warnings":warnings,"warnings_truncated":warnings_truncated}),
        )
    }
}
