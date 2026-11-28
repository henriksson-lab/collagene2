use crate::platform::worker;
use collagene_core::documents::transfer;
use collagene_core::material::lineage::EventBinding;
use collagene_core::protocol;
use collagene_core::sequence::analysis::Primer;
use collagene_core::verification;
use collagene_core::{
    alignment,
    cloning::{pcr::primers, restriction::setups as digest_setups},
    documents::{history, workspace},
    formats::oligos as oligo_formats,
    material::{lineage as run_lineage, outputs as cloning_outputs},
    review::contigs::history as contig_history,
    sequence::{analysis, annotations::library as feature_library},
};
use collagene_model::Molecule;
use leptos::prelude::*;

pub(crate) type OpenDocument = collagene_core::documents::history::DocumentHistory;
pub(crate) const MAX_WORKSPACE_BYTES: usize = 64 * 1024 * 1024;
pub(crate) fn assembly_draft_kind(kind: &str) -> bool {
    kind == crate::assembly::draft::model::ATTACHMENT_KIND
}
pub(crate) fn current_event_binding(
    doc: &OpenDocument,
) -> collagene_core::material::lineage::EventBinding {
    let revision = doc.head_revision().to_owned();
    collagene_core::material::lineage::EventBinding {
        document_id: doc.id.clone(),
        event_id: doc.head_event().into(),
        revision,
    }
}
pub(crate) fn recombination_source_info(
    doc: &OpenDocument,
) -> crate::sources::identity::SourceInfo {
    crate::sources::identity::SourceInfo {
        binding: current_event_binding(doc),
        name: doc.molecule.name.clone(),
        length: doc.molecule.sequence.len(),
        alphabet: doc.molecule.sequence.alphabet(),
        topology: doc.molecule.topology,
    }
}

pub(crate) fn bound_primer_source<'a>(
    documents: &'a [OpenDocument],
    binding: &collagene_core::material::lineage::EventBinding,
) -> Result<&'a Molecule, String> {
    documents
        .iter()
        .find(|doc| current_event_binding(doc) == *binding)
        .map(|doc| &doc.molecule)
        .ok_or_else(|| {
            "The saved-primer document event changed; select its current revision.".into()
        })
}

pub(crate) type Attachments =
    std::collections::BTreeMap<String, collagene_core::documents::workspace::WorkspaceAttachment>;
pub(crate) fn open_histories(
    molecules: Vec<Molecule>,
    existing: &[OpenDocument],
) -> Result<Vec<OpenDocument>, String> {
    let mut ids = existing
        .iter()
        .map(|d| d.id.clone())
        .collect::<std::collections::BTreeSet<_>>();
    let mut result = Vec::new();
    for molecule in molecules {
        let prefix = format!("document-{}", &molecule.revision()[..16]);
        let mut suffix = 1;
        let id = loop {
            let id = format!("{prefix}-{suffix}");
            if ids.insert(id.clone()) {
                break id;
            }
            suffix += 1;
        };
        result.push(OpenDocument::new(id, molecule).map_err(|e| e.to_string())?);
    }
    Ok(result)
}
pub(crate) struct RestoredWorkspace {
    pub(crate) name: String,
    pub(crate) documents: Vec<OpenDocument>,
    pub(crate) active: usize,
    pub(crate) candidate: Option<transfer::Candidate>,
    pub(crate) candidates: Vec<transfer::Candidate>,
    pub(crate) attachments: Attachments,
    pub(crate) oligos: Vec<std::sync::Arc<collagene_core::formats::oligos::OligoTable>>,
    pub(crate) libraries: Vec<crate::viewer::sequence::feature_tools::LibraryRevision>,
    pub(crate) alignments: Vec<crate::viewer::alignment::alignment_view::AlignmentRevision>,
    pub(crate) runs: Vec<collagene_core::material::lineage::RunSummary>,
    pub(crate) assemblies: Vec<collagene_core::review::contigs::history::AssemblySummary>,
    pub(crate) assembly_draft: crate::assembly::draft::model::Draft,
    pub(crate) pair_draft: crate::primers::primer_pair_model::PairDraft,
    pub(crate) hidden_products: std::collections::BTreeSet<String>,
}
impl RestoredWorkspace {
    pub(crate) fn from_frame(
        frame: &collagene_core::documents::workspace::Workspace,
    ) -> Result<RestoredWorkspace, String> {
        crate::viewer::sequence::annotation_visibility::validate_workspace(frame)?;
        let documents = OpenDocument::restore_all(frame).map_err(|e| e.to_string())?;
        if frame
            .attachments
            .values()
            .any(|asset| asset.kind == collagene_core::cloning::restriction::setups::KIND)
        {
            collagene_core::cloning::restriction::setups::list(frame).map_err(|e| e.to_string())?;
        }
        if documents
            .iter()
            .map(|d| document_bytes(&d.molecule))
            .sum::<usize>()
            > MAX_WORKSPACE_BYTES
        {
            return Err("Restored current documents exceed this web workspace's 64 MiB limit; the archive has not been changed.".into());
        }
        let active = documents
            .iter()
            .position(|d| d.id == frame.active_document)
            .ok_or("Active workspace document is unavailable")?;
        let mut oligos = vec![];
        for attachment in frame
            .attachments
            .values()
            .filter(|a| a.kind == "oligo-table-v1")
        {
            let table: collagene_core::formats::oligos::OligoTable =
                serde_json::from_slice(&attachment.bytes)
                    .map_err(|_| "Invalid oligo collection in workspace")?;
            table.validate().map_err(|e| e.to_string())?;
            oligos.push(std::sync::Arc::new(table));
        }
        let libraries = frame
            .attachments
            .values()
            .filter(|attachment| attachment.kind == "feature-library-v1")
            .map(|attachment| {
                collagene_core::sequence::annotations::library::load_library(&attachment.bytes)
                    .map(crate::viewer::sequence::feature_tools::LibraryRevision::new)
                    .map_err(|error| error.to_string())
            })
            .collect::<Result<Vec<_>, _>>()?;
        let alignments = frame
            .attachments
            .values()
            .filter(|a| a.kind == crate::viewer::alignment::alignment_view::ATTACHMENT_KIND)
            .map(|a| {
                crate::viewer::alignment::alignment_view::AlignmentRevision::from_bytes(&a.bytes)
            })
            .collect::<Result<Vec<_>, _>>()?;
        if frame.attachments.values().any(|asset| {
            (asset.kind.starts_with("assembly-draft-v")
                && asset.kind != crate::assembly::draft::model::ATTACHMENT_KIND)
                || (asset.kind.starts_with("primer-pair-draft-v")
                    && asset.kind != crate::primers::primer_pair_storage::ATTACHMENT_KIND)
        }) {
            return Err("Workspace contains an unsupported internal draft format.".into());
        }
        let draft_assets = frame
            .attachments
            .values()
            .filter(|a| assembly_draft_kind(&a.kind))
            .collect::<Vec<_>>();
        if draft_assets.len() > 1 {
            return Err("Workspace contains multiple active assembly drafts.".into());
        }
        let assembly_draft = match draft_assets.first() {
            Some(asset) if asset.media_type == "application/json" => {
                crate::assembly::draft::model::decode(&asset.bytes)?
            }
            Some(_) => return Err("Assembly draft has the wrong media type.".into()),
            None => crate::assembly::draft::model::Draft::default(),
        };
        for captured in crate::assembly::draft::assembly_alternatives::all_captures(&assembly_draft)
        {
            let binding = &captured.source.binding;
            let source = documents
                .iter()
                .find(|d| d.id == binding.document_id)
                .ok_or("Assembly draft source document is not retained.")?;
            let exact = source
                .revision_molecule(&binding.event_id)
                .map_err(|e| e.to_string())?;
            if exact.revision() != binding.revision
                || exact.revision() != captured.molecule.revision()
            {
                return Err(
                    "Assembly draft source event differs from its retained snapshot.".into(),
                );
            }
        }
        let pair_assets = frame
            .attachments
            .values()
            .filter(|a| a.kind == crate::primers::primer_pair_storage::ATTACHMENT_KIND)
            .collect::<Vec<_>>();
        if pair_assets.len() > 1 {
            return Err("Workspace contains multiple active primer pair drafts.".into());
        }
        let pair_draft = match pair_assets.first() {
            Some(asset) if asset.media_type == "application/json" => {
                crate::primers::primer_pair_storage::decode(&asset.bytes, |binding| {
                    let source = documents
                        .iter()
                        .find(|d| d.id == binding.document_id)
                        .ok_or("Primer pair draft source document is not retained.")?;
                    source
                        .revision_molecule(&binding.event_id)
                        .map_err(|e| e.to_string())
                })?
            }
            Some(_) => return Err("Primer pair draft has the wrong media type.".into()),
            None => crate::primers::primer_pair_model::PairDraft::default(),
        };
        Ok(RestoredWorkspace {
            hidden_products: collagene_core::material::outputs::load_hidden(frame)
                .map_err(|e| e.to_string())?,
            assembly_draft,
            pair_draft,
            runs: collagene_core::material::lineage::list_runs(frame).map_err(|e| e.to_string())?,
            assemblies: collagene_core::review::contigs::history::list(frame)
                .map_err(|e| e.to_string())?,
            name: frame.name.clone(),
            documents,
            active,
            candidate: frame
                .active_candidate
                .as_ref()
                .and_then(|id| frame.candidates.iter().find(|c| &c.plan_revision == id))
                .cloned(),
            candidates: frame.candidates.clone(),
            attachments: frame.attachments.clone(),
            oligos,
            libraries,
            alignments,
        })
    }
}
pub(crate) fn document_bytes(m: &Molecule) -> usize {
    m.sequence.len()
        + m.name.len()
        + m.metadata
            .iter()
            .map(|(key, value)| key.len() + value.len())
            .sum::<usize>()
        + m.original_genbank.as_ref().map_or(0, String::len)
        + m.features
            .iter()
            .map(|f| {
                128 + f
                    .qualifiers
                    .iter()
                    .map(|(k, v)| k.len() + v.as_ref().map_or(0, String::len))
                    .sum::<usize>()
            })
            .sum::<usize>()
}
pub(crate) async fn encode_workspace(
    frame: collagene_core::documents::workspace::Workspace,
    id: u32,
) -> Result<Vec<u8>, String> {
    let revision = frame
        .head(&frame.active_document)
        .map_err(|e| e.to_string())?
        .revision();
    match worker::analyze(protocol::request::Request {
        id,
        revision: revision.clone(),
        job: protocol::schema::AnalysisJob::SaveWorkspace { workspace: frame },
    })
    .await?
    {
        protocol::response::Response {
            id: returned,
            revision: bound,
            result: protocol::schema::AnalysisResult::WorkspaceSaved(bytes),
        } if returned == id && bound == revision => Ok(bytes),
        protocol::response::Response {
            result: protocol::schema::AnalysisResult::Error(error),
            ..
        } => Err(error),
        _ => Err("Unexpected workspace-save response".into()),
    }
}
pub(crate) async fn decode_workspace(
    bytes: Vec<u8>,
    id: u32,
) -> Result<collagene_core::documents::workspace::Workspace, String> {
    let revision = blake3::hash(&bytes).to_hex().to_string();
    match worker::analyze(protocol::request::Request {
        id,
        revision: revision.clone(),
        job: protocol::schema::AnalysisJob::LoadWorkspace { bytes },
    })
    .await?
    {
        protocol::response::Response {
            id: returned,
            revision: bound,
            result: protocol::schema::AnalysisResult::WorkspaceLoaded(frame),
        } if returned == id && bound == revision => Ok(frame),
        protocol::response::Response {
            result: protocol::schema::AnalysisResult::Error(error),
            ..
        } => Err(error),
        _ => Err("Unexpected workspace-load response".into()),
    }
}

#[derive(Clone, Copy)]
pub(crate) struct WorkspaceSession {
    pub(crate) documents: RwSignal<Vec<OpenDocument>>,
    pub(crate) active: RwSignal<usize>,
    pub(crate) workspace_restoring: RwSignal<bool>,
    pub(crate) workspace_epoch: RwSignal<u64>,
    pub(crate) workspace_name: RwSignal<String>,
    pub(crate) attachments: RwSignal<Attachments>,
}
impl WorkspaceSession {
    pub(crate) fn new(molecule: Molecule) -> Self {
        Self {
            documents: RwSignal::new(vec![
                OpenDocument::new("document-demo", molecule).expect("valid demonstration history"),
            ]),
            active: RwSignal::new(0usize),
            workspace_restoring: RwSignal::new(true),
            workspace_epoch: RwSignal::new(0u64),
            workspace_name: RwSignal::new("Collagene workspace".to_string()),
            attachments: RwSignal::new(Attachments::new()),
        }
    }
    /// Advance identity before installing any state from a replacement frame.
    pub(crate) fn begin_replacement(&self) {
        self.workspace_epoch
            .update(|epoch| *epoch = epoch.wrapping_add(1));
    }
    /// The frame has already passed RestoredWorkspace's full bounded validation.
    pub(crate) fn install_documents(
        &self,
        name: String,
        documents: Vec<OpenDocument>,
        active: usize,
        attachments: Attachments,
    ) {
        self.documents.set(documents);
        self.active.set(active);
        self.workspace_name.set(name);
        self.attachments.set(attachments);
    }
    pub(crate) fn snapshot(
        documents: &[OpenDocument],
        active: usize,
        current: Option<&transfer::Candidate>,
        candidates: &[transfer::Candidate],
        attachments: &Attachments,
        name: &str,
    ) -> Result<collagene_core::documents::workspace::Workspace, String> {
        use collagene_core::documents::workspace::Workspace;
        let active_document = documents
            .get(active)
            .ok_or("No active workspace document")?
            .id
            .clone();
        let mut frame = Workspace {
            name: name.into(),
            active_document,
            documents: vec![],
            revisions: std::collections::BTreeMap::new(),
            history: vec![],
            candidates: vec![],
            active_candidate: current.map(|c| c.plan_revision.clone()),
            attachments: attachments.clone(),
        };
        for doc in documents {
            doc.contribute(&mut frame).map_err(|e| e.to_string())?;
        }
        for candidate in candidates.iter().chain(current) {
            if frame
                .candidates
                .iter()
                .any(|c| c.plan_revision == candidate.plan_revision)
            {
                continue;
            }
            for input in &candidate.plan.inputs {
                frame
                    .revisions
                    .entry(input.revision.clone())
                    .or_insert_with(|| input.molecule.clone());
            }
            frame.candidates.push(candidate.clone());
        }
        Ok(frame)
    }
    pub(crate) fn add_documents(&self, molecules: Vec<Molecule>) -> Result<(), String> {
        let Self {
            documents,
            active,
            workspace_restoring,
            ..
        } = *self;
        if workspace_restoring.get_untracked() {
            return Err("Workspace restoration is in progress; retry after it finishes.".into());
        }
        if molecules.is_empty() {
            return Ok(());
        }
        let first = documents.with_untracked(Vec::len);
        let existing = documents.with_untracked(|docs| {
            docs.iter()
                .map(|d| document_bytes(&d.molecule))
                .sum::<usize>()
        });
        if first + molecules.len() > 1000
            || existing + molecules.iter().map(document_bytes).sum::<usize>() > MAX_WORKSPACE_BYTES
        {
            return Err("This prototype supports at most 1,000 open documents and 64 MiB of document data. Export your work and open fewer records.".into());
        }
        let additions = documents.with_untracked(|docs| open_histories(molecules, docs))?;
        documents.update(|docs| docs.extend(additions));
        active.set(first);
        Ok(())
    }

    pub(crate) fn export_context_guard(&self) -> Callback<(), Result<(), String>> {
        let Self {
            documents,
            active,
            workspace_restoring: restoring,
            workspace_epoch: epoch,
            ..
        } = *self;
        let opened_epoch = epoch.get_untracked();
        let pin = documents.with_untracked(|docs| {
            docs.get(active.get_untracked()).map(|d| {
                (
                    d.id.clone(),
                    d.head_event().to_string(),
                    d.molecule.revision(),
                )
            })
        });
        Callback::new(move |()| {
            if epoch.try_get_untracked() != Some(opened_epoch)
                || restoring.try_get_untracked().unwrap_or(true)
            {
                return Err(
                    "The workspace changed; reopen the export against its current source.".into(),
                );
            }
            let index = active
                .try_get_untracked()
                .ok_or("The export view closed.")?;
            let current = documents
                .try_with_untracked(|docs| {
                    docs.get(index).map(|d| {
                        (
                            d.id.clone(),
                            d.head_event().to_string(),
                            d.molecule.revision(),
                        )
                    })
                })
                .ok_or("The export view closed.")?;
            if current != pin {
                return Err(
                    "The source document/version changed; prepare the export again.".into(),
                );
            }
            Ok(())
        })
    }

    pub(crate) fn commit_document(&self, base: &str, edited: Molecule) -> Result<(), String> {
        let Self {
            documents,
            active,
            workspace_restoring,
            ..
        } = *self;
        if workspace_restoring.get_untracked() {
            return Err("Workspace restoration is in progress; retry after it finishes.".into());
        }
        edited.validate().map_err(|e| e.to_string())?;
        let index = active.get_untracked();
        let valid = documents.with_untracked(|docs| {
            docs.get(index)
                .is_some_and(|d| d.molecule.revision() == base)
        });
        if !valid {
            return Err("The active document changed while the editor was open. Reopen the editor against its current revision.".into());
        }
        let bytes = documents.with_untracked(|docs| {
            docs.iter()
                .enumerate()
                .filter(|(i, _)| *i != index)
                .map(|(_, d)| document_bytes(&d.molecule))
                .sum::<usize>()
        }) + document_bytes(&edited);
        if bytes > MAX_WORKSPACE_BYTES {
            return Err("This edit exceeds the current workspace size limit.".into());
        }
        let mut result = Ok(());
        documents.update(|docs| {
            result = docs[index]
                .commit(edited, "Edited document")
                .map_err(|e| e.to_string());
        });
        result
    }
}
