use crate::{renderer::MapRenderer, storage, worker};
fn dismiss_viewer_menus(event: web_sys::MouseEvent) {
    let Some(target) = event
        .target()
        .and_then(|t| t.dyn_into::<web_sys::Element>().ok())
    else {
        return;
    };
    let menus =
        ".selection-action-menu,.toolbar-action-menu,.document-file-menu,.viewer-tools-menu";
    let inside = target.closest(menus).ok().flatten().is_some();
    if inside
        && (target.closest("summary").ok().flatten().is_some()
            || target.closest("button").ok().flatten().is_none())
    {
        return;
    }
    if let Some(document) = web_sys::window().and_then(|w| w.document()) {
        if let Ok(open)=document.query_selector_all(".selection-action-menu[open],.toolbar-action-menu[open],.document-file-menu[open],.viewer-tools-menu[open]") {
            for i in 0..open.length(){if let Some(node)=open.item(i).and_then(|n|n.dyn_into::<web_sys::Element>().ok()){let _=node.remove_attribute("open");}}
        }
    }
}

use collagene_core::{
    analysis::Primer,
    formats, protocol, transfer,
    verification::{self, DesignPlan, Operation, PlanInput, Source, VerificationReport},
    *,
};
use leptos::{html, prelude::*};
use wasm_bindgen::JsCast;
use wasm_bindgen_futures::{JsFuture, spawn_local};
use web_sys::{HtmlCanvasElement, HtmlInputElement};

type OpenDocument = collagene_core::history::DocumentHistory;
const MAX_WORKSPACE_BYTES: usize = 64 * 1024 * 1024;
fn assembly_draft_kind(kind: &str) -> bool {
    kind == crate::assembly_draft_model::ATTACHMENT_KIND
}
fn current_event_binding(doc: &OpenDocument) -> collagene_core::run_lineage::EventBinding {
    let revision = doc.head_revision().to_owned();
    collagene_core::run_lineage::EventBinding {
        document_id: doc.id.clone(),
        event_id: doc.head_event().into(),
        revision,
    }
}
fn recombination_source_info(doc: &OpenDocument) -> crate::digest::SourceInfo {
    crate::digest::SourceInfo {
        binding: current_event_binding(doc),
        name: doc.molecule.name.clone(),
        length: doc.molecule.sequence.len(),
        alphabet: doc.molecule.sequence.alphabet(),
        topology: doc.molecule.topology,
    }
}

fn bound_primer_source<'a>(
    documents: &'a [OpenDocument],
    binding: &collagene_core::run_lineage::EventBinding,
) -> Result<&'a Molecule, String> {
    documents
        .iter()
        .find(|doc| current_event_binding(doc) == *binding)
        .map(|doc| &doc.molecule)
        .ok_or_else(|| {
            "The saved-primer document event changed; select its current revision.".into()
        })
}

type Attachments =
    std::collections::BTreeMap<String, collagene_core::workspace::WorkspaceAttachment>;
fn open_histories(
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
fn workspace_snapshot(
    documents: &[OpenDocument],
    active: usize,
    current: Option<&transfer::Candidate>,
    candidates: &[transfer::Candidate],
    attachments: &Attachments,
    name: &str,
) -> Result<collagene_core::workspace::Workspace, String> {
    use collagene_core::workspace::Workspace;
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
#[derive(Clone)]
struct PendingFinalize {
    workspace: std::sync::Arc<collagene_core::workspace::Workspace>,
    plan: DesignPlan,
    sources: Vec<collagene_core::run_lineage::InputEventBinding>,
    source_fingerprint: String,
}
#[derive(Clone)]
struct PendingReadReview {
    request: collagene_core::read_review::Request,
    owner: collagene_core::run_lineage::EventBinding,
    epoch: u64,
}
struct RestoredWorkspace {
    name: String,
    documents: Vec<OpenDocument>,
    active: usize,
    candidate: Option<transfer::Candidate>,
    candidates: Vec<transfer::Candidate>,
    attachments: Attachments,
    oligos: Vec<std::sync::Arc<collagene_core::oligo_formats::OligoTable>>,
    libraries: Vec<crate::feature_tools::LibraryRevision>,
    alignments: Vec<crate::alignment_view::AlignmentRevision>,
    runs: Vec<collagene_core::run_lineage::RunSummary>,
    assemblies: Vec<collagene_core::contig_history::AssemblySummary>,
    assembly_draft: crate::assembly_workspace::Draft,
    pair_draft: crate::primer_pair_model::PairDraft,
    hidden_products: std::collections::BTreeSet<String>,
}
fn restore_workspace_state(
    frame: &collagene_core::workspace::Workspace,
) -> Result<RestoredWorkspace, String> {
    crate::annotation_visibility::validate_workspace(frame)?;
    let documents = OpenDocument::restore_all(frame).map_err(|e| e.to_string())?;
    if frame
        .attachments
        .values()
        .any(|asset| asset.kind == collagene_core::digest_setups::KIND)
    {
        collagene_core::digest_setups::list(frame).map_err(|e| e.to_string())?;
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
        let table: collagene_core::oligo_formats::OligoTable =
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
            collagene_core::feature_library::load_library(&attachment.bytes)
                .map(crate::feature_tools::LibraryRevision::new)
                .map_err(|error| error.to_string())
        })
        .collect::<Result<Vec<_>, _>>()?;
    let alignments = frame
        .attachments
        .values()
        .filter(|a| a.kind == crate::alignment_view::ATTACHMENT_KIND)
        .map(|a| crate::alignment_view::AlignmentRevision::from_bytes(&a.bytes))
        .collect::<Result<Vec<_>, _>>()?;
    if frame.attachments.values().any(|asset| {
        (asset.kind.starts_with("assembly-draft-v")
            && asset.kind != crate::assembly_draft_model::ATTACHMENT_KIND)
            || (asset.kind.starts_with("primer-pair-draft-v")
                && asset.kind != crate::primer_pair_storage::ATTACHMENT_KIND)
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
            crate::assembly_draft_model::decode(&asset.bytes)?
        }
        Some(_) => return Err("Assembly draft has the wrong media type.".into()),
        None => crate::assembly_workspace::Draft::default(),
    };
    for captured in crate::assembly_alternatives::all_captures(&assembly_draft) {
        let binding = &captured.source.binding;
        let source = documents
            .iter()
            .find(|d| d.id == binding.document_id)
            .ok_or("Assembly draft source document is not retained.")?;
        let exact = source
            .revision_molecule(&binding.event_id)
            .map_err(|e| e.to_string())?;
        if exact.revision() != binding.revision || exact.revision() != captured.molecule.revision()
        {
            return Err("Assembly draft source event differs from its retained snapshot.".into());
        }
    }
    let pair_assets = frame
        .attachments
        .values()
        .filter(|a| a.kind == crate::primer_pair_storage::ATTACHMENT_KIND)
        .collect::<Vec<_>>();
    if pair_assets.len() > 1 {
        return Err("Workspace contains multiple active primer pair drafts.".into());
    }
    let pair_draft = match pair_assets.first() {
        Some(asset) if asset.media_type == "application/json" => {
            crate::primer_pair_storage::decode(&asset.bytes, |binding| {
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
        None => crate::primer_pair_model::PairDraft::default(),
    };
    Ok(RestoredWorkspace {
        hidden_products: collagene_core::cloning_outputs::load_hidden(frame)
            .map_err(|e| e.to_string())?,
        assembly_draft,
        pair_draft,
        runs: collagene_core::run_lineage::list_runs(frame).map_err(|e| e.to_string())?,
        assemblies: collagene_core::contig_history::list(frame).map_err(|e| e.to_string())?,
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
fn document_bytes(m: &Molecule) -> usize {
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
    frame: collagene_core::workspace::Workspace,
    id: u32,
) -> Result<Vec<u8>, String> {
    let revision = frame
        .head(&frame.active_document)
        .map_err(|e| e.to_string())?
        .revision();
    match worker::analyze(protocol::Request {
        id,
        revision: revision.clone(),
        job: protocol::AnalysisJob::SaveWorkspace { workspace: frame },
    })
    .await?
    {
        protocol::Response {
            id: returned,
            revision: bound,
            result: protocol::AnalysisResult::WorkspaceSaved(bytes),
        } if returned == id && bound == revision => Ok(bytes),
        protocol::Response {
            result: protocol::AnalysisResult::Error(error),
            ..
        } => Err(error),
        _ => Err("Unexpected workspace-save response".into()),
    }
}
pub(crate) async fn decode_workspace(
    bytes: Vec<u8>,
    id: u32,
) -> Result<collagene_core::workspace::Workspace, String> {
    let revision = blake3::hash(&bytes).to_hex().to_string();
    match worker::analyze(protocol::Request {
        id,
        revision: revision.clone(),
        job: protocol::AnalysisJob::LoadWorkspace { bytes },
    })
    .await?
    {
        protocol::Response {
            id: returned,
            revision: bound,
            result: protocol::AnalysisResult::WorkspaceLoaded(frame),
        } if returned == id && bound == revision => Ok(frame),
        protocol::Response {
            result: protocol::AnalysisResult::Error(error),
            ..
        } => Err(error),
        _ => Err("Unexpected workspace-load response".into()),
    }
}

fn demonstration() -> Molecule {
    let mut seed = 271828u32;
    let mut bases: String = (0..4800)
        .map(|_| {
            seed = seed.wrapping_mul(1664525).wrapping_add(1013904223);
            b"ACGT"[(seed >> 24) as usize % 4] as char
        })
        .collect();
    for (pos, site) in [
        (250, "GAATTC"),
        (1250, "GGATCC"),
        (2700, "AAGCTT"),
        (4050, "GCGGCCGC"),
    ] {
        bases.replace_range(pos..pos + site.len(), site);
    }
    let mut mol = Molecule::dna("pCG_demo", &bases, Topology::Circular).expect("synthetic DNA");
    for (kind, label, start, end, reverse) in [
        ("promoter", "Promoter example", 120, 320, false),
        ("CDS", "Reporter example", 450, 1250, false),
        ("terminator", "Terminator example", 1300, 1420, false),
        ("CDS", "Marker example", 1900, 2750, true),
        (
            "rep_origin",
            "Replication origin example",
            3200,
            3700,
            false,
        ),
    ] {
        let location = Location::range(start, end);
        mol.features.push(Feature {
            kind: kind.into(),
            location: if reverse {
                Location::Complement(Box::new(location))
            } else {
                location
            },
            qualifiers: vec![("label".into(), Some(label.into()))],
            needs_review: false,
        });
    }
    mol.features.push(Feature {
        kind: "misc_feature".into(),
        location: Location::Join(vec![Location::range(4650, 4800), Location::range(0, 60)]),
        qualifiers: vec![("label".into(), Some("Across-origin example".into()))],
        needs_review: false,
    });
    mol.metadata
        .insert("collagene:example".into(), "interface-demonstration".into());
    mol.metadata.insert("description".into(),"Synthetic interface demonstration. Annotations illustrate the UI and do not describe a lab-ready construct.".into());
    mol
}

thread_local! { static EXPORT_NOTICE: std::cell::RefCell<Option<RwSignal<String>>> = const { std::cell::RefCell::new(None) }; }
/// An owned, bounded export snapshot. Native admission precedes its byte copy.
pub(crate) struct ExportArtifact {
    name: String,
    mime: String,
    bytes: Vec<u8>,
    reservation: Option<crate::desktop_rpc::ExportReservation>,
    notice: Option<RwSignal<String>>,
}
pub(crate) fn prepare_download(
    name: &str,
    mime: &str,
    bytes: &[u8],
) -> Result<ExportArtifact, String> {
    if bytes.len() > 64 * 1024 * 1024 {
        return Err("Export exceeds the 64 MiB file limit.".into());
    }
    let reservation = if crate::desktop_rpc::is_desktop() {
        Some(crate::desktop_rpc::reserve_export()?)
    } else {
        None
    };
    let name = if reservation.is_some() {
        let mut safe: String = name
            .chars()
            .map(|c| {
                if c.is_control() || matches!(c, '/' | '\\' | ':') {
                    '_'
                } else {
                    c
                }
            })
            .collect();
        while safe.len() > 256 {
            safe.pop();
        }
        if safe.trim().is_empty() {
            safe = "export.bin".into();
        }
        safe
    } else {
        name.into()
    };
    Ok(ExportArtifact {
        name,
        mime: mime.into(),
        bytes: bytes.to_vec(),
        reservation,
        notice: EXPORT_NOTICE.with(|notice| *notice.borrow()),
    })
}
impl ExportArtifact {
    pub async fn save(self, guard: impl Fn() -> Result<(), String>) -> Result<(), String> {
        if let Some(reservation) = self.reservation {
            let name =
                crate::desktop_rpc::export_reserved(reservation, self.name, self.bytes, guard)
                    .await?;
            if let Some(notice) = self.notice {
                let _ = notice.try_set(format!("Saved export: exports/{name}"));
            }
            Ok(())
        } else {
            guard()?;
            browser_download(&self.name, &self.mime, &self.bytes)
        }
    }
    pub fn complete(self, done: impl FnOnce(Result<(), String>) + 'static) {
        if self.reservation.is_none() {
            done(browser_download(&self.name, &self.mime, &self.bytes));
        } else {
            wasm_bindgen_futures::spawn_local(async move { done(self.save(|| Ok(())).await) });
        }
    }
    pub fn complete_guarded(
        self,
        guard: impl Fn() -> Result<(), String> + 'static,
        done: impl FnOnce(Result<(), String>) + 'static,
    ) {
        if self.reservation.is_none() {
            done(guard().and_then(|()| browser_download(&self.name, &self.mime, &self.bytes)));
        } else {
            wasm_bindgen_futures::spawn_local(async move { done(self.save(guard).await) });
        }
    }
}
pub(crate) async fn download(name: &str, mime: &str, bytes: &[u8]) -> Result<(), String> {
    prepare_download(name, mime, bytes)?.save(|| Ok(())).await
}
pub(crate) async fn download_guarded(
    name: &str,
    mime: &str,
    bytes: &[u8],
    guard: impl Fn() -> Result<(), String>,
) -> Result<(), String> {
    prepare_download(name, mime, bytes)?.save(guard).await
}
/// Completion runs after the native durable receipt or browser dispatch, never
/// merely after queueing. Capturing owned bytes keeps later edits out of the file.
pub(crate) fn download_with_completion(
    name: &str,
    mime: &str,
    bytes: &[u8],
    done: impl FnOnce(Result<(), String>) + 'static,
) {
    match prepare_download(name, mime, bytes) {
        Ok(artifact) => artifact.complete(done),
        Err(error) => done(Err(error)),
    }
}
pub(crate) fn download_with_completion_guarded(
    name: &str,
    mime: &str,
    bytes: &[u8],
    guard: impl Fn() -> Result<(), String> + 'static,
    done: impl FnOnce(Result<(), String>) + 'static,
) {
    complete_download_guarded(prepare_download(name, mime, bytes), guard, done);
}
pub(crate) fn complete_download(
    result: Result<ExportArtifact, String>,
    done: impl FnOnce(Result<(), String>) + 'static,
) {
    match result {
        Ok(artifact) => artifact.complete(done),
        Err(error) => done(Err(error)),
    }
}
pub(crate) fn complete_download_guarded(
    result: Result<ExportArtifact, String>,
    guard: impl Fn() -> Result<(), String> + 'static,
    done: impl FnOnce(Result<(), String>) + 'static,
) {
    match result {
        Ok(artifact) => artifact.complete_guarded(guard, done),
        Err(error) => done(Err(error)),
    }
}

fn browser_download(name: &str, mime: &str, bytes: &[u8]) -> Result<(), String> {
    let array = js_sys::Uint8Array::from(bytes);
    let parts = js_sys::Array::new();
    parts.push(&array);
    let options = web_sys::BlobPropertyBag::new();
    options.set_type(mime);
    let blob = web_sys::Blob::new_with_u8_array_sequence_and_options(&parts, &options)
        .map_err(|_| "Cannot prepare download")?;
    let url =
        web_sys::Url::create_object_url_with_blob(&blob).map_err(|_| "Cannot create download")?;
    let window = web_sys::window().ok_or("No browser window")?;
    let document = window.document().ok_or("No browser document")?;
    let anchor = document
        .create_element("a")
        .map_err(|_| "Cannot create download link")?;
    anchor
        .set_attribute("href", &url)
        .map_err(|_| "Cannot set download URL")?;
    anchor
        .set_attribute("download", name)
        .map_err(|_| "Cannot set filename")?;
    let anchor = anchor.unchecked_into::<web_sys::HtmlElement>();
    if let Some(body) = document.body() {
        let _ = body.append_child(&anchor);
    }
    anchor.click();
    anchor.remove();
    // The click queues navigation; defer revocation until the browser consumes the URL.
    let cleanup = wasm_bindgen::closure::Closure::once_into_js(move || {
        let _ = web_sys::Url::revoke_object_url(&url);
    });
    let _ =
        window.set_timeout_with_callback_and_timeout_and_arguments_0(cleanup.unchecked_ref(), 1000);
    Ok(())
}

fn add_documents(
    documents: RwSignal<Vec<OpenDocument>>,
    active: RwSignal<usize>,
    workspace_restoring: RwSignal<bool>,
    molecules: Vec<Molecule>,
) -> Result<(), String> {
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

fn export_context_guard(
    documents: RwSignal<Vec<OpenDocument>>,
    active: RwSignal<usize>,
    restoring: RwSignal<bool>,
    epoch: RwSignal<u64>,
) -> Callback<(), Result<(), String>> {
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
            return Err("The source document/version changed; prepare the export again.".into());
        }
        Ok(())
    })
}

fn commit_document(
    documents: RwSignal<Vec<OpenDocument>>,
    active: RwSignal<usize>,
    workspace_restoring: RwSignal<bool>,
    base: &str,
    edited: Molecule,
) -> Result<(), String> {
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

pub(crate) fn describe_step(step: &Operation) -> String {
    let source = |s: &Source| match s {
        Source::Input(id) => format!("input {id}"),
        Source::Step(i) => i
            .checked_add(1)
            .map(|number| format!("product of step {number}"))
            .unwrap_or_else(|| format!("invalid step index {i}")),
    };
    match step {
        Operation::OverlapExtension {
            left,
            right,
            options,
        } => format!(
            "Reciprocal overlap extension of {} ({}) / {} ({})\nSelected terminal homology: {} bp\nRetain actual source strands and extend paired 3′ ends; alternative terminal pairings remain diagnostic.",
            source(&left.source),
            if left.reverse {
                "reverse complement"
            } else {
                "forward"
            },
            source(&right.source),
            if right.reverse {
                "reverse complement"
            } else {
                "forward"
            },
            options.overlap
        ),
        Operation::ReconstructMaterial { template, origin } => format!(
            "Reconstruct reusable PCR/digest material from {} · exact event {}",
            source(template),
            origin.product.event_id
        ),

        Operation::ReconstructRunProduct { template, origin } => format!(
            "Reconstruct saved material from {}\nOrigin run: {} · product {}\nUse the unchanged originating recipe and target; retain nicks, termini, positional modifications and all incomplete findings.",
            source(template),
            origin.selection.run_id,
            origin.selection.step_index + 1
        ),
        Operation::CatalogDigest {
            template,
            selection,
            fragment_index,
        } => format!(
            "Catalog-pinned complete digest of {}\nEnzymes: {}\nCatalog edition: {}\nZero-based fragment: {}\nCleavage assumes compatible conditions and no blocking methylation; no experimental success is asserted.",
            source(template),
            selection
                .enzymes
                .iter()
                .map(|e| e.enzyme.name.as_str())
                .collect::<Vec<_>>()
                .join(", "),
            selection.catalog_edition,
            fragment_index
        ),
        Operation::DeclareSourceTermini {
            template,
            left,
            right,
        } => format!(
            "Original intact linear duplex: {}\nDeclared left terminus: {:?}\nDeclared right terminus: {:?}\nOperator assumptions only; no measurement, synthesis or physical graph repair is inferred.",
            source(template),
            left,
            right
        ),
        Operation::TopoActions { sources, request } => format!(
            "Accounted TOPO chemical actions\n{}\nLiteral declarations: {}\nIdeal hydrolytic profiles are explicit model assumptions; excised material remains evidence, not a claimed isolated product.",
            sources
                .iter()
                .map(|s| format!("{} = {}", s.role, source(&s.source)))
                .collect::<Vec<_>>()
                .join("\n"),
            serde_json::to_string_pretty(request).unwrap_or_default()
        ),
        Operation::ProcessedTopo { sources, request } => format!(
            "Source-backed TOPO transfer, explicit ligase treatment and faithful duplex bridge\n{}\nLiteral processing request: {}\nOther components remain evidence; no experimental isolation or measured success is inferred.",
            sources
                .iter()
                .map(|s| format!("{} = {}", s.role, source(&s.source)))
                .collect::<Vec<_>>()
                .join("\n"),
            serde_json::to_string_pretty(request).unwrap_or_default()
        ),
        Operation::CheckedTopo {
            sources,
            request,
            selection,
            policy,
        } => format!(
            "Final source-backed TOPO graph with explicit independent-target acceptance\n{}\nPreparation/transfers: {}\nSelected component/seed: {}\nAcceptance policy: {}\nAccepted allowances do not repair nicks, eliminate other components or establish experimental success. The sequence view remains chemically restricted.",
            sources
                .iter()
                .map(|entry| format!("{} = {}", entry.role, source(&entry.source)))
                .collect::<Vec<_>>()
                .join("\n"),
            serde_json::to_string_pretty(request).unwrap_or_default(),
            serde_json::to_string(selection).unwrap_or_default(),
            serde_json::to_string_pretty(policy).unwrap_or_default()
        ),
        Operation::Topo {
            sources,
            request,
            selection,
        } => format!(
            "Source-backed TOPO strand transfers\n{}\nPreparation and transfer recipe: {}\nSelected component: {}\nPhysical graph and paired-contour evidence are retained in the reconstruction report. A reference preview does not establish repaired, intact duplex chemistry or experimental success.",
            sources
                .iter()
                .map(|s| format!("{} = {}", s.role, source(&s.source)))
                .collect::<Vec<_>>()
                .join("\n"),
            serde_json::to_string_pretty(request).unwrap_or_default(),
            serde_json::to_string(selection).unwrap_or_default()
        ),
        Operation::MultiSiteRecombineReference {
            sources,
            request,
            selection,
        } => format!(
            "Multi-parent {:?} reference exchange\n{}\nProfile {} / {}\nSpecificities {:?}\nSelected candidate {}, component {} from batch {}. All reciprocal components retained; no chemical strand authority.",
            request.reaction.reaction,
            sources
                .iter()
                .map(|s| format!("{} = {}", s.role, source(&s.source)))
                .collect::<Vec<_>>()
                .join("\n"),
            request.reaction.profile.id,
            request.reaction.profile.edition,
            request.reaction.specificities,
            selection.candidate_id + 1,
            selection.component_id + 1,
            selection.batch_revision
        ),
        Operation::RecombineReference {
            sources,
            request,
            selection,
        } => format!(
            "Reference-only {:?} exchange of {} and {}\nExplicit profile {} / {} (revision {})\nSelected candidate {}, component {} from batch {}.\nAll reciprocal components are computed under the declared full-site profile and bounded limits. Selection does not establish physical strand state, biological repair, efficacy or an experimental winner.",
            request.reaction,
            source(&sources[0]),
            source(&sources[1]),
            request.profile.profile.id,
            request.profile.profile.edition,
            request.profile.revision,
            selection.candidate_id,
            selection.component_id,
            selection.batch_revision
        ),
        Operation::DeclareSourceModifications {
            template,
            assumptions,
        } => format!(
            "Declare source generation for {}\nDam: {}\nDcm: {}\nCpG: {}\nClosed Dam/Dcm/CpG model: modification is assigned only to eligible positions in the original source motifs. Other original positions and newly synthesized bases are unmodified under this scope. Descendant motifs are not rescanned to infer inheritance. This is an assumption, not imported physical evidence.",
            source(template),
            crate::checked_cloning::state_label(assumptions.dam),
            crate::checked_cloning::state_label(assumptions.dcm),
            crate::checked_cloning::state_label(assumptions.cpg)
        ),
        Operation::CheckedDigest {
            template,
            selection,
            options,
            fragment_index,
        } => format!(
            "Checked restriction digest of {}\nEnzymes: {}\nPinned positional supplier profile revision: {}\nSelect fragment {} in ascending reference top-cut order. Unknown or impaired sites withhold undetermined products; blocked sites remain uncut under this model. Reagent conditions and experimental success remain outside the check.",
            source(template),
            selection
                .enzymes
                .iter()
                .map(|e| e.enzyme.name.as_str())
                .collect::<Vec<_>>()
                .join(", "),
            format!(
                "{} · catalog {} · {} comparisons / {} sites",
                modified_duplex::profile_revision(),
                selection.catalog_content_hash,
                options.max_comparisons,
                options.max_sites
            ),
            fragment_index
                .checked_add(1)
                .map(|n| n.to_string())
                .unwrap_or_else(|| "out of range".into())
        ),
        Operation::Tail3Prime { template, base } => format!(
            "Substrate: {}\nAdd one non-templated {} to both actual 3′ ends. Source end chemistry and nicks are retained; enzyme activity and completion are assumptions.",
            source(template),
            base.as_str()
        ),
        Operation::ResectThreePrime {
            source: substrate,
            cuts,
        } => format!(
            "Explicit 3′ resection of {}\nLeft bottom-strand cut: {} bp · right top-strand cut: {} bp\nAccessible 3′ suffixes are removed; original 5′ groups remain. No fill or sealing is inferred.",
            source(substrate),
            cuts.left,
            cuts.right
        ),
        Operation::HomologyAnnealing {
            fragments, policy, ..
        } => format!(
            "3′ resection and homology annealing\n{} oriented fragments · {:?}\nJoining homology: {:?} bp\nExplicit cuts: {:?}\nThe joined reference remains incomplete. Unfilled gaps and unsealed nicks are retained; no reusable assembled material or host processing is inferred.",
            fragments.len(),
            policy.topology,
            policy.overlaps,
            policy.cuts
        ),
        Operation::HomologyAssembly {
            fragments, policy, ..
        } => format!(
            "Physical homology assembly\n{} oriented fragments · {:?}\nJoining homology: {:?} bp\n5′ resection: {:?}\nOuter-end policy: {:?}\nGeneric resection/annealing/gap-fill/ligation assumptions; inspect processing evidence before saving reusable material.",
            fragments.len(),
            policy.topology,
            policy.overlaps,
            policy.resection,
            policy.outer_end_policy
        ),
        Operation::ProfiledOverlap {
            fragments,
            overlaps,
            circular,
            options,
        } => format!(
            "{:?} overlap profile\nOrdered fragments: {}\nJunction homology: {:?}bp\nTopology: {}\nPublished length guidance is advisory; unique reaction outcomes and reagent activity are not established.",
            options.method,
            fragments
                .iter()
                .map(|f| format!(
                    "{} ({})",
                    source(&f.source),
                    if f.reverse {
                        "reverse complement"
                    } else {
                        "forward"
                    }
                ))
                .collect::<Vec<_>>()
                .join(" → "),
            overlaps,
            if *circular { "circular" } else { "linear" }
        ),
        Operation::RepairEnds { template, method } => format!(
            "Substrate: {}\nEnd repair: {:?}\nOnly modeled 3′ fill-in/exonuclease changes occur; 5′ chemistry and internal nicks remain explicit.",
            source(template),
            method
        ),
        Operation::AnnealOligos {
            top,
            bottom,
            options,
            candidate_index,
        } => format!(
            "Ordered 5′→3′ oligos: {} / {}\nExplicit zero-based annealing candidate {}\nDeclared chemistry/search policy: {:?}\nAnnealing creates no new phosphodiester bonds.",
            source(top),
            source(bottom),
            candidate_index,
            options
        ),
        Operation::PcrWithMismatches {
            template,
            forward,
            reverse,
            options,
        } => format!(
            "Template: {}\nOligos 5′→3′: {} / {}\nBinding/search policy: {:?}",
            source(template),
            forward.ordered_sequence(),
            reverse.ordered_sequence(),
            options
        ),
        Operation::SelectedGoldenGate {
            templates,
            options,
            candidate_index,
            directed_sequence_key,
            physical_sequence_key,
        } => format!(
            "Golden Gate inputs: {}\nExplicit zero-based candidate {}\nDirected sequence key {}\nPhysical sequence key {}\nPolicy: {:?}",
            templates.iter().map(source).collect::<Vec<_>>().join(", "),
            candidate_index,
            directed_sequence_key,
            physical_sequence_key,
            options
        ),
        Operation::GoldenGate {
            templates,
            options,
            candidate_index,
        } => format!(
            "Golden Gate inputs: {}\nExplicit zero-based candidate {}\nPolicy: {:?}",
            templates.iter().map(source).collect::<Vec<_>>().join(", "),
            candidate_index,
            options
        ),
        Operation::Pcr {
            template,
            forward,
            reverse,
        } => format!(
            "Template: {}\nForward primer {} (5′→3′): {}\nReverse primer {} (5′→3′): {}",
            source(template),
            forward.name,
            forward.ordered_sequence(),
            reverse.name,
            reverse.ordered_sequence()
        ),
        Operation::Digest {
            template,
            enzymes,
            fragment_index,
        } => format!(
            "Template: {}\nComplete digest with: {}\nSelected zero-based fragment index: {} (ascending top-cut order)",
            source(template),
            enzymes.join(", "),
            fragment_index
        ),
        Operation::TreatEnds {
            template,
            treatment,
        } => format!(
            "Substrate: {}\nTerminal treatment: {:?}\nInternal nick repair is not assumed.",
            source(template),
            treatment
        ),
        Operation::Ligate {
            fragments,
            circular,
            allow_nicked,
            options,
        } => format!(
            "Ordered substrates: {}\nProduct topology: {}\nRetained nicks explicitly allowed: {}\nAlternative-enumeration bounds: {:?}",
            fragments
                .iter()
                .map(|f| format!(
                    "{} ({})",
                    source(&f.source),
                    if f.reverse { "reversed" } else { "forward" }
                ))
                .collect::<Vec<_>>()
                .join(" → "),
            if *circular { "circular" } else { "linear" },
            allow_nicked,
            options
        ),
        Operation::Overlap {
            fragments,
            overlaps,
            circular,
        } => format!(
            "Fragment order: {}\nJunction overlaps (bp): {}\nProduct topology: {}",
            fragments.iter().map(source).collect::<Vec<_>>().join(" → "),
            overlaps
                .iter()
                .map(usize::to_string)
                .collect::<Vec<_>>()
                .join(", "),
            if *circular { "circular" } else { "linear" }
        ),
    }
}

#[component]
fn RestrictionFields(
    side: &'static str,
    enzymes: RwSignal<String>,
    fragment: RwSignal<String>,
    treatment: RwSignal<String>,
    repair: RwSignal<String>,
    reverse: RwSignal<bool>,
) -> impl IntoView {
    view! {
        <h4>{format!("{side} substrate")}</h4>
        <label>{format!("{side} enzymes (comma separated)")}<input aria-label=format!("{side} enzymes") prop:value=move||enzymes.get() on:input=move|e|enzymes.set(event_target_value(&e))/></label>
        <label>{format!("{side} digest fragment")}<input aria-label=format!("{side} fragment number") type="number" min="1" prop:value=move||fragment.get() on:input=move|e|fragment.set(event_target_value(&e))/></label>
        <label>{format!("{side} end repair")}<select aria-label=format!("{side} end repair") prop:value=move||repair.get() on:change=move|e|repair.set(event_target_value(&e))><option value="none">"No end repair"</option><option value="blunt">"Blunt: fill 5′ overhangs and remove 3′ overhangs"</option><option value="fill">"Fill 5′ overhangs only"</option><option value="chew">"Remove 3′ overhangs only"</option></select></label>
        <label>{format!("{side} end treatment")}<select aria-label=format!("{side} end treatment") prop:value=move||treatment.get() on:change=move|e|treatment.set(event_target_value(&e))><option value="none">"No treatment"</option><option value="phosphorylate">"Phosphorylate 5′ termini"</option><option value="dephosphorylate">"Dephosphorylate 5′ termini"</option></select></label>
        <label class="checkbox-label"><input type="checkbox" aria-label=format!("Reverse {side} fragment") prop:checked=move||reverse.get() on:change=move|e|reverse.set(event_target::<HtmlInputElement>(&e).checked())/>{format!("Reverse {side} fragment")}</label>
    }
}

#[component]
pub fn App() -> impl IntoView {
    let native_export_notice = RwSignal::new(String::new());
    EXPORT_NOTICE.with(|notice| *notice.borrow_mut() = Some(native_export_notice));
    on_cleanup(move || {
        EXPORT_NOTICE.with(|notice| {
            if *notice.borrow() == Some(native_export_notice) {
                *notice.borrow_mut() = None;
            }
        })
    });
    crate::transient_menus::install();
    let documents = RwSignal::new(vec![
        OpenDocument::new("document-demo", demonstration()).expect("valid demonstration history"),
    ]);
    let workspace_restoring = RwSignal::new(true);
    let workspace_epoch = RwSignal::new(0u64);
    let workspace_name = RwSignal::new("Collagene workspace".to_string());
    let attachments = RwSignal::new(Attachments::new());
    let genome_links = Memo::new(move |_| {
        attachments.with(|assets| {
            assets
                .values()
                .filter(|a| a.kind == collagene_core::genome::WORKSPACE_ATTACHMENT_KIND)
                .filter_map(|a| collagene_core::genome::load_link(&a.bytes).ok())
                .collect::<Vec<_>>()
        })
    });
    let recorded_runs = RwSignal::new(Vec::<collagene_core::run_lineage::RunSummary>::new());
    let saved_assemblies =
        RwSignal::new(Vec::<collagene_core::contig_history::AssemblySummary>::new());
    let assembly_page = RwSignal::new(0usize);
    let pending_assembly = RwSignal::new(None::<crate::contigs::ArchiveRequest>);
    let pending_donor = RwSignal::new(None::<crate::crispr::DonorRequest>);
    let pending_digest = RwSignal::new(None::<crate::digest_product::DigestOpenRequest>);
    let pending_material = RwSignal::new(None::<crate::digest_product::MaterialOpenRequest>);
    let gel_digest_seeds = RwSignal::new(Vec::<crate::gel::DigestLaneSeed>::new());
    let inspected_assembly = RwSignal::new(String::new());
    let selected_run = RwSignal::new(String::new());
    let internal_clipboard =
        RwSignal::new(None::<std::sync::Arc<crate::reference_transfer_view::Copied>>);
    let pending_reference_paste = RwSignal::new(None::<crate::reference_transfer_view::PasteOpen>);
    let pending_comparison_source = RwSignal::new(None::<crate::pairwise_source_view::Snapshot>);
    let pending_translation_source = RwSignal::new(None::<crate::pairwise_source_view::Snapshot>);
    let pending_pcr_source = RwSignal::new(None::<crate::pcr_source_view::Snapshot>);
    let pending_run_event = RwSignal::new(None::<std::sync::Arc<crate::run_event_model::Snapshot>>);
    let run_page = RwSignal::new(0usize);
    let oligo_tables = RwSignal::new(Vec::<
        std::sync::Arc<collagene_core::oligo_formats::OligoTable>,
    >::new());
    let libraries = RwSignal::new(Vec::<crate::feature_tools::LibraryRevision>::new());
    let alignments = RwSignal::new(Vec::<crate::alignment_view::AlignmentRevision>::new());
    let pending_alignment = RwSignal::new(None::<crate::alignment_view::PendingAlignment>);
    let pending_sbol = RwSignal::new(None::<std::sync::Arc<Vec<u8>>>);
    let method_config = RwSignal::new(crate::cloning_methods::MethodsState::default());
    let assembly_draft = RwSignal::new(crate::assembly_workspace::Draft::default());
    let previous_assembly_draft =
        RwSignal::new(None::<std::sync::Arc<crate::assembly_workspace::Draft>>);
    Effect::new(move |_| {
        let _ = workspace_epoch.get();
        previous_assembly_draft.set(None);
    });
    let assembly_visible = RwSignal::new(false);
    let pending_finalize = RwSignal::new(None::<std::sync::Arc<PendingFinalize>>);
    let pending_release = RwSignal::new(None::<String>);
    let selected_release_manifest = RwSignal::new(None::<String>);
    let read_reference = RwSignal::new(None::<collagene_core::run_lineage::EventBinding>);
    let pending_read_check = RwSignal::new(None::<crate::trace_validation::RecordRequest>);
    let reviewed_check = RwSignal::new(None::<String>);
    let review_seed = RwSignal::new(None::<collagene_core::read_review::Request>);
    let pending_read_review = RwSignal::new(None::<PendingReadReview>);
    let hidden_products = RwSignal::new(std::collections::BTreeSet::<String>::new());
    let show_prepared_history = RwSignal::new(false);
    let digest_boundary_pick = RwSignal::new(None::<crate::digest::BoundaryPick>);
    let digest_run_origin = RwSignal::new(None::<crate::run_origin_model::PinnedOrigin>);
    Effect::new(move || {
        let _ = workspace_epoch.get();
        digest_run_origin.set(None);
    });
    let annotation_source = RwSignal::new(String::new());
    let annotation_format =
        RwSignal::new(collagene_core::annotation_formats::AnnotationFormat::Gff3);
    let existing_oligo = RwSignal::new(None::<collagene_core::oligo_formats::OligoTable>);
    let pending_oligo = RwSignal::new(None::<(String, oligo_formats::Delimiter)>);
    let pending_fastq = RwSignal::new(None::<std::sync::Arc<Vec<u8>>>);
    let pending_variants = RwSignal::new(
        None::<(
            std::sync::Arc<Vec<u8>>,
            collagene_core::variants_session::FileFormat,
            String,
        )>,
    );
    let pending_read_alignment = RwSignal::new(
        None::<(
            std::sync::Arc<Vec<u8>>,
            collagene_core::read_alignment_session::FileFormat,
            String,
        )>,
    );
    let fastq_encoding = RwSignal::new(None::<oligo_formats::QualityEncoding>);
    let stored_candidates = RwSignal::new(Vec::<transfer::Candidate>::new());
    let active = RwSignal::new(0usize);
    // Working tabs are presentation state, separate from retained histories and
    // the archive's last working document index (which must remain valid).
    let open_tabs = RwSignal::new(crate::document_tabs_model::Tabs::default());
    let tab_epoch = StoredValue::new(None::<u64>);
    let tab_document = StoredValue::new(None::<String>);
    let retained_document_ids = Memo::new(move |_| {
        documents.with(|docs| docs.iter().map(|doc| doc.id.clone()).collect::<Vec<_>>())
    });
    let has_open_document = Memo::new(move |_| open_tabs.with(|tabs| tabs.active().is_some()));
    Effect::new(move |_| {
        let epoch = workspace_epoch.get();
        let index = active.get();
        let id = documents.with_untracked(|docs| docs.get(index).map(|doc| doc.id.clone()));
        let current = open_tabs.get_untracked();
        let mut next = current.clone();
        if tab_epoch.get_value() != Some(epoch) {
            // A finalized assembly also replaces the workspace frame. Keep
            // working source tabs whose IDs remain in that frame; a different
            // workspace naturally removes IDs it does not retain.
            let retained = documents
                .with_untracked(|docs| docs.iter().map(|doc| doc.id.clone()).collect::<Vec<_>>());
            next.reconcile(&retained);
            if next.ids().is_empty() {
                next.reset(id.as_deref());
            } else if let Some(id) = &id {
                next.open(id);
            }
        } else if let Some(id) = &id {
            next.open(id);
        }
        tab_epoch.set_value(Some(epoch));
        tab_document.set_value(id);
        if next != current {
            open_tabs.set(next);
        }
    });
    Effect::new(move |_| {
        let retained = retained_document_ids.get();
        let id = documents
            .with_untracked(|docs| docs.get(active.get_untracked()).map(|doc| doc.id.clone()));
        let current = open_tabs.get_untracked();
        let mut next = current.clone();
        next.reconcile(&retained);
        // A replacement can change the active document ID without changing its
        // array index. Ordinary revision edits do not reopen closed tabs.
        if tab_document.get_value() != id {
            if let Some(id) = &id {
                next.open(id);
            }
            tab_document.set_value(id);
        }
        if next != current {
            open_tabs.set(next);
        }
    });
    let molecule = Memo::new(move |_| {
        documents.with(|docs| docs[active.get().min(docs.len() - 1)].molecule.clone())
    });
    let molecule_revision = Memo::new(move |_| molecule.with(Molecule::revision));
    let selection = RwSignal::new(None::<(usize, usize)>);
    let selected_feature = RwSignal::new(None::<usize>);
    let feature_for_editor = RwSignal::new(None::<usize>);
    let feature_initial_location = RwSignal::new(None::<Location>);
    let primer_seed = RwSignal::new(None::<Primer>);
    let pair_draft = RwSignal::new(crate::primer_pair_model::PairDraft::default());
    let mutation_seed = RwSignal::new(None::<crate::mutagenesis::MutationSelection>);
    let inspected_feature = Memo::new(move |_| {
        molecule.with(|m| {
            selected_feature
                .get()
                .and_then(|i| m.features.get(i).cloned())
        })
    });
    let search_regions = RwSignal::new(None::<(String, Vec<(usize, usize)>)>);
    // Navigation receipts are local UI state, never imported scientific evidence.
    // Initialize before effects so the first tick cannot erase a pointer selection.
    let sequence_transitions = StoredValue::new(crate::sequence_edit_transition::Tracker::new(
        documents.with_untracked(|docs| {
            docs.get(active.get_untracked()).map(|doc| {
                crate::sequence_edit_transition::SourcePin::of(workspace_epoch.get_untracked(), doc)
            })
        }),
    ));
    let sequence_on_edit = Callback::new(
        move |(base, a, b, replacement): (String, usize, usize, String)| {
            let current = molecule.get_untracked();
            if current.revision() != base {
                return Err("The sequence changed before this edit could be applied.".into());
            }
            let before = documents
                .with_untracked(|docs| {
                    docs.get(active.get_untracked()).map(|doc| {
                        crate::sequence_edit_transition::SourcePin::of(
                            workspace_epoch.get_untracked(),
                            doc,
                        )
                    })
                })
                .ok_or_else(|| "The source document is no longer open.".to_string())?;
            let edited = current
                .replace(a, b, &replacement)
                .map_err(|e| e.to_string())?;
            // Sequence::new normalizes whitespace/case. Raw clipboard byte length
            // must never determine the committed cursor boundary.
            let replacement_len = before
                .length
                .checked_sub(b.saturating_sub(a))
                .and_then(|retained| edited.sequence.len().checked_sub(retained));
            commit_document(documents, active, workspace_restoring, &base, edited)?;
            let after = documents.with_untracked(|docs| {
                docs.get(active.get_untracked()).map(|doc| {
                    crate::sequence_edit_transition::SourcePin::of(
                        workspace_epoch.get_untracked(),
                        doc,
                    )
                })
            });
            let transition = after.zip(replacement_len).and_then(|(after, n)| {
                crate::sequence_edit_transition::EditTransition::attest(before, after, a, b, n)
            });
            sequence_transitions.update_value(|tracker| {
                if let Some(transition) = transition {
                    tracker.record(transition);
                } else {
                    tracker.invalidate();
                }
            });
            Ok(())
        },
    );
    let sequence_on_history = Callback::new(move |redo: bool| {
        if workspace_restoring.get_untracked() {
            return Err("Workspace restoration is in progress.".into());
        }
        let mut result = Ok(false);
        documents.update(|docs| {
            let doc = &mut docs[active.get_untracked()];
            result = if redo { doc.redo() } else { doc.undo() }.map_err(|e| e.to_string());
        });
        let changed = result?;
        if changed {
            sequence_transitions.update_value(|tracker| tracker.invalidate());
            batch(|| {
                selection.set(None);
                selected_feature.set(None);
                search_regions.set(None);
            });
        }
        Ok(())
    });
    let biology_directed = RwSignal::new(None::<crate::biology_selection::DirectedSelection>);
    let highlighted_spans = Memo::new(move |_| {
        if let Some(feature) = inspected_feature.get() {
            return feature
                .location
                .spans()
                .into_iter()
                .map(|(a, b, _)| (a, b))
                .collect::<Vec<_>>();
        }
        if let Some(span) = selection.get() {
            return vec![span];
        }
        search_regions.with(|regions| {
            regions
                .as_ref()
                .filter(|(revision, _)| *revision == molecule_revision.get())
                .map(|(_, parts)| parts.clone())
                .unwrap_or_default()
        })
    });
    let selection_status = Memo::new(move |_| {
        if let Some((a, b)) = selection.get().filter(|(a, b)| a == b) {
            let _ = b;
            return format!("Insertion boundary {a} · type to insert");
        }
        let parts = highlighted_spans.get();
        if parts.len() > 1 {
            let positions = parts
                .iter()
                .take(4)
                .map(|(a, b)| format!("{}–{}", a + 1, b))
                .collect::<Vec<_>>()
                .join(", ");
            return format!(
                "Selection {positions}{} · {} reference spans",
                if parts.len() > 4 { "…" } else { "" },
                parts.len()
            );
        }
        selection
            .get()
            .map(|(a, b)| {
                format!(
                    "Selection {}–{} · {} {}",
                    a + 1,
                    b,
                    b - a,
                    molecule.with(|m| m.sequence.alphabet().map_unit())
                )
            })
            .unwrap_or_else(|| {
                molecule.with(|m| {
                    if m.sequence.alphabet() == Alphabet::Protein {
                        "Click a residue or feature to select".into()
                    } else {
                        "Click a base or feature to select".into()
                    }
                })
            })
    });
    let panel = RwSignal::new("");
    let hovered_restriction = RwSignal::new(None::<usize>);
    let message = RwSignal::new(String::new());
    let modal = RwSignal::new("");
    let modal_error = RwSignal::new(String::new());
    let open_bulk_release = Callback::new(move |hash: String| {
        if workspace_restoring.get_untracked() {
            return;
        }
        pending_release.set(Some(hash));
        modal_error.set(String::new());
        modal.set("combinatorial_release");
    });
    let released_bulk = Callback::new(move |hash: String| {
        selected_release_manifest.set(Some(hash));
        assembly_visible.set(true);
        pending_release.set(None);
        modal.set("");
    });
    let import_text = RwSignal::new(String::new());
    let pending_compat_import = RwSignal::new(
        None::<(
            std::sync::Arc<Vec<u8>>,
            String,
            crate::sequence_compat_import::Profile,
        )>,
    );
    let edit_text = RwSignal::new(String::new());
    let trace_view = RwSignal::new(false);
    let file_input = NodeRef::<html::Input>::new();
    let gl_canvas = NodeRef::<html::Canvas>::new();
    let label_canvas = NodeRef::<html::Canvas>::new();
    let renderer = StoredValue::new_local(None::<MapRenderer>);
    let render_error = RwSignal::new(String::new());
    let map_labels = RwSignal::new(crate::map_labels::LabelStats::default());
    let editor_view = RwSignal::new("split");
    let show_library = RwSignal::new(true);
    let layout = crate::workspace_layout::use_layout(
        editor_view,
        show_library,
        Memo::new(move |_| workspace_restoring.get()),
    );
    let restriction_identity = Memo::new(move |_| {
        documents.with(|docs| {
            let doc = &docs[active.get()];
            format!("{}:{}:{}", workspace_epoch.get(), doc.id, doc.head_event())
        })
    });
    let annotation_hidden = Memo::new(move |_| {
        documents.with(|docs| {
            let doc = &docs[active.get()];
            attachments.with(|assets| {
                match crate::annotation_visibility::load(
                    assets,
                    &current_event_binding(doc),
                    doc.molecule.features.len(),
                ) {
                    Ok(hidden) => hidden,
                    Err(error) => {
                        message.set(format!(
                            "Annotation display preferences need attention: {error}"
                        ));
                        Default::default()
                    }
                }
            })
        })
    });
    let restriction_selection = Memo::new(move |_| {
        molecule.with(|m| {
            if let Some(index) = selected_feature.get() {
                return m.features.get(index).and_then(|feature| {
                    collagene_core::restriction_selection::footprint(
                        &feature.location,
                        m.sequence.len(),
                    )
                    .ok()
                });
            }
            if let Some((revision, spans)) = search_regions.get() {
                if revision == molecule_revision.get() && !spans.is_empty() {
                    return Some(spans);
                }
            }
            selection.get().and_then(|(start, end)| {
                (start < end && end <= m.sequence.len()).then_some(vec![(start, end)])
            })
        })
    });
    let restrictions = crate::restriction_view::use_restrictions(
        molecule,
        restriction_identity,
        Memo::new(move |_| workspace_restoring.get()),
        restriction_selection,
    );
    let navigate_restriction = Callback::new(
        move |site: collagene_core::restriction_inspection::InspectionSite| {
            let parts = site
                .spans
                .iter()
                .map(|p| (p.start, p.end))
                .collect::<Vec<_>>();
            batch(|| {
                selected_feature.set(None);
                restrictions.selected.set(Some(site.index));
                selection.set(if parts.len() == 1 {
                    parts.first().copied()
                } else {
                    None
                });
                search_regions.set(Some((molecule_revision.get_untracked(), parts)));
                if assembly_visible.get_untracked() {
                    let enzyme = restrictions.catalog.with_value(|catalog| {
                        catalog
                            .iter()
                            .find(|enzyme| enzyme.name == site.enzyme)
                            .map(
                                |enzyme| collagene_core::restriction_inspection::EnzymeSpec {
                                    name: enzyme.name.clone(),
                                    recognition: enzyme.recognition.clone(),
                                    top_cut: enzyme.top_cut,
                                    bottom_cut: enzyme.bottom_cut,
                                },
                            )
                    });
                    if let Some(enzyme) = enzyme {
                        let source = documents.with_untracked(|docs| {
                            current_event_binding(&docs[active.get_untracked()])
                        });
                        digest_boundary_pick.set(Some(crate::digest::BoundaryPick {
                            source,
                            enzyme,
                            recognition_start: site.recognition_start,
                            reverse: site.reverse,
                        }));
                    }
                    panel.set("enzymes");
                } else {
                    panel.set("restriction");
                }
            });
        },
    );
    let document_items = Memo::new(move |_| {
        documents.with(|docs| {
            docs.iter()
                .map(|doc| {
                    let provider = doc.molecule.metadata.get("source:provider");
                    let accession = doc
                        .molecule
                        .metadata
                        .get("source:accession_version")
                        .or_else(|| doc.molecule.metadata.get("source:identifier"));
                    crate::document_explorer::DocumentItem {
                        id: doc.id.clone(),
                        name: doc.molecule.name.clone(),
                        length: doc.molecule.sequence.len(),
                        alphabet: doc.molecule.sequence.alphabet(),
                        topology: doc.molecule.topology,
                        source: if doc
                            .molecule
                            .metadata
                            .get("assembly:intent:reference_only")
                            .or_else(|| doc.molecule.metadata.get("mutation:intent:reference_only"))
                            .or_else(|| doc.molecule.metadata.get("read_review:reference_only"))
                            .is_some_and(|value| value == "true")
                        {
                            Some("Intended reference".into())
                        } else {
                            provider.map(|provider| {
                                accession
                                    .map(|id| format!("{provider} · {id}"))
                                    .unwrap_or_else(|| provider.clone())
                            })
                        },
                    }
                })
                .collect::<Vec<_>>()
        })
    });
    let explorer_active =
        Memo::new(move |_| open_tabs.with(|tabs| tabs.active().unwrap_or_default().to_owned()));
    let explorer_documents = Memo::new(move |_| {
        let show_all = show_prepared_history.get();
        let hidden = hidden_products.get();
        let current = explorer_active.get();
        document_items.with(|items| {
            items
                .iter()
                .filter(|doc| show_all || !hidden.contains(&doc.id) || doc.id == current)
                .cloned()
                .collect::<Vec<_>>()
        })
    });
    let tab_items = Memo::new(move |_| {
        let index = document_items.with(|items| {
            items
                .iter()
                .map(|item| (item.id.clone(), item.clone()))
                .collect::<std::collections::BTreeMap<_, _>>()
        });
        open_tabs.with(|tabs| {
            tabs.ids()
                .iter()
                .filter_map(|id| index.get(id).cloned())
                .collect::<Vec<_>>()
        })
    });
    let explorer_locked = Memo::new(move |_| workspace_restoring.get());
    let explorer_activate = Callback::new(move |id: String| {
        if workspace_restoring.get_untracked() {
            return;
        }
        if let Some(index) =
            documents.with_untracked(|docs| docs.iter().position(|doc| doc.id == id))
        {
            batch(|| {
                open_tabs.update(|tabs| tabs.open(&id));
                active.set(index);
            });
        }
    });
    let close_document_tab = Callback::new(move |id: String| {
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
    });

    let db = StoredValue::new_local(None::<crate::storage::Database>);
    let native_library_context = RwSignal::new(None::<collagene_library_protocol::LabContext>);
    let loaded = RwSignal::new(false);
    let save_status = RwSignal::new("Opening local storage…".to_string());
    let save_generation = RwSignal::new(0u64);
    let pcr_forward_name = RwSignal::new("Forward".to_string());
    let pcr_reverse_name = RwSignal::new("Reverse".to_string());
    let pcr_forward = RwSignal::new(String::new());
    let pcr_reverse = RwSignal::new(String::new());
    let pcr_forward_tail = RwSignal::new(String::new());
    let pcr_reverse_tail = RwSignal::new(String::new());
    let pcr_generation = RwSignal::new(0u32);
    let pcr_review = RwSignal::new(None::<crate::pcr_outputs::Review>);
    let pcr_origin = RwSignal::new(None::<crate::run_origin_model::PinnedOrigin>);
    Effect::new(move |_| {
        let current = documents.with(|docs| docs.get(active.get()).map(current_event_binding));
        if pcr_origin
            .get_untracked()
            .as_ref()
            .is_some_and(|origin| current.as_ref() != Some(origin.event()))
        {
            pcr_origin.set(None);
        }
    });
    let published_pcr_overlay = RwSignal::new(None::<crate::pcr_overlay::Publication>);
    let pcr_output_busy = RwSignal::new(false);
    let pcr_output_worker = StoredValue::new_local(None::<std::rc::Rc<worker::CachedWorker>>);
    let pcr_alive = StoredValue::new(std::sync::Arc::new(std::sync::atomic::AtomicBool::new(
        true,
    )));
    let pcr_cleanup_alive = pcr_alive.get_value();
    on_cleanup(move || {
        pcr_cleanup_alive.store(false, std::sync::atomic::Ordering::Relaxed);
        if let Some(worker) = pcr_output_worker.get_value() {
            worker.close();
        }
    });
    Effect::new(move |_| {
        let _ = (
            pcr_forward_name.get(),
            pcr_reverse_name.get(),
            pcr_forward.get(),
            pcr_reverse.get(),
            pcr_forward_tail.get(),
            pcr_reverse_tail.get(),
            workspace_epoch.get(),
            workspace_restoring.get(),
            pcr_origin.get(),
        );
        documents.with(|docs| {
            let _ = docs.get(active.get()).map(current_event_binding);
        });
        pcr_generation.update(|g| *g = g.wrapping_add(1));
        pcr_review.set(None);
        if published_pcr_overlay.with_untracked(|p| {
            p.as_ref()
                .and_then(|p| p.overlay.as_ref())
                .is_some_and(|overlay| overlay.kind == crate::pcr_overlay::Kind::ExactProduct)
        }) {
            published_pcr_overlay.set(None);
        }
        pcr_output_busy.set(false);
        if let Some(worker) = pcr_output_worker.get_value() {
            worker.close();
        }
        pcr_output_worker.set_value(None);
    });
    Effect::new(move |_| {
        let _ = workspace_epoch.get();
        published_pcr_overlay.set(None);
        pending_pcr_source.set(None);
    });
    let pcr_overlay = Memo::new(move |_| {
        if workspace_restoring.get() {
            return None;
        }
        let epoch = workspace_epoch.get();
        let inputs = (
            pcr_forward_name.get(),
            pcr_reverse_name.get(),
            pcr_forward.get(),
            pcr_reverse.get(),
            pcr_forward_tail.get(),
            pcr_reverse_tail.get(),
        );
        documents.with(|docs| {
            let doc = &docs[active.get()];
            let binding = current_event_binding(doc);
            if let Some(overlay) = published_pcr_overlay.get().and_then(|p| p.overlay) {
                let fresh = overlay.is_current(&binding, epoch)
                    && (overlay.kind != crate::pcr_overlay::Kind::ExactProduct
                        || pcr_review.with(|review| {
                            review.as_ref().is_some_and(|r| {
                                r.source == binding && r.epoch == epoch && r.inputs == inputs
                            })
                        }));
                if fresh {
                    return Some(overlay);
                }
            }
            let forward = Primer::new(&inputs.0, &inputs.2, &inputs.4).ok()?;
            let reverse = Primer::new(&inputs.1, &inputs.3, &inputs.5).ok()?;
            pair_draft.with(|draft| {
                crate::pcr_overlay::captured(
                    &doc.molecule,
                    &binding,
                    epoch,
                    draft,
                    &forward,
                    &reverse,
                )
                .ok()
                .flatten()
                .map(std::sync::Arc::new)
            })
        })
    });
    let publish_pcr_overlay = Callback::new(move |publication: crate::pcr_overlay::Publication| {
        match publication.overlay.as_ref() {
            Some(overlay) => {
                let valid = !workspace_restoring.get_untracked()
                    && workspace_epoch.get_untracked() == overlay.epoch
                    && documents.with_untracked(|docs| {
                        docs.iter().any(|d| {
                            current_event_binding(d) == overlay.source
                                && overlay.spans.iter().all(|span| {
                                    span.start < span.end && span.end <= d.molecule.sequence.len()
                                })
                        })
                    });
                if valid {
                    published_pcr_overlay.set(Some(publication));
                }
            }
            None => {
                if published_pcr_overlay
                    .with_untracked(|p| p.as_ref().is_some_and(|p| p.owner == publication.owner))
                {
                    published_pcr_overlay.set(None);
                }
            }
        }
    });
    let primer_sources = Memo::new(move |_| {
        documents.with(|docs| {
            docs.iter()
                .filter(|doc| doc.molecule.sequence.alphabet() == Alphabet::Dna)
                .map(recombination_source_info)
                .collect::<Vec<_>>()
        })
    });
    let list_saved_primers = Callback::new(move |binding: collagene_core::run_lineage::EventBinding| -> Result<collagene_core::primer_reuse::List,String> {
        if workspace_restoring.get_untracked() {return Err("Wait for workspace restoration to finish.".into());}
        documents.with_untracked(|docs| {
            let source=bound_primer_source(docs,&binding)?;
            let list=collagene_core::primer_reuse::list(source,1000).map_err(|e|e.to_string())?;
            if list.source_revision!=binding.revision {return Err("The saved-primer source revision changed.".into());}
            Ok(list)
        })
    });
    let resolve_saved_primer = Callback::new(move |request: crate::primer_picker::ResolveRequest| -> Result<collagene_core::primer_reuse::Resolved,String> {
        if workspace_restoring.get_untracked() {return Err("Wait for workspace restoration to finish.".into());}
        documents.with_untracked(|docs| {
            let source=bound_primer_source(docs,&request.source)?;
            collagene_core::primer_reuse::resolve_at_revision(source,request.feature_index,&request.source.revision).map_err(|e|e.to_string())
        })
    });

    Effect::new(move |_| {
        if panel.get() == "reads" && read_reference.get().is_none() {
            read_reference.set(Some(documents.with_untracked(|docs| {
                current_event_binding(&docs[active.get_untracked()])
            })));
        }
    });
    let read_check_source = Memo::new(move |_| {
        let pinned = read_reference.get();
        documents.with(|docs| {
            let doc = pinned
                .as_ref()
                .and_then(|pin| docs.iter().find(|doc| doc.id == pin.document_id))
                .unwrap_or_else(|| &docs[active.get()]);
            recombination_source_info(doc)
        })
    });
    let read_check_attachment_available = Callback::new(move |hash: String| {
        attachments.with(|assets| {
            assets
                .get(&hash)
                .is_some_and(|asset| asset.kind == collagene_core::read_checks::ATTACHMENT_KIND)
        })
    });
    let read_check_reference = Memo::new(move |_| read_check_source.get().binding);
    let read_check_molecule = Memo::new(move |_| {
        let id = read_check_source.get().binding.document_id;
        documents.with(|docs| {
            docs.iter()
                .find(|doc| doc.id == id)
                .unwrap_or_else(|| &docs[active.get()])
                .molecule
                .clone()
        })
    });
    let primer_context =
        Memo::new(move |_| documents.with(|docs| current_event_binding(&docs[active.get()])));
    let pair_selection = Callback::new(move |()| {
        if workspace_restoring.get_untracked() {
            return Err("Wait for the workspace operation to finish.".into());
        }
        documents.with_untracked(|docs| {
            let source = &docs[active.get_untracked()];
            let (start, end) = crate::selection_toolbar::primer_range(
                &source.molecule,
                selection.get_untracked(),
                selected_feature.get_untracked(),
                &highlighted_spans.get_untracked(),
            )?;
            crate::primer_pair_model::SelectionSeed::capture(
                &source.molecule,
                &current_event_binding(source),
                start,
                end,
            )
        })
    });
    let navigate_primer_binding = Callback::new(
        move |(source, spans): (
            collagene_core::run_lineage::EventBinding,
            Vec<(usize, usize)>,
        )|
              -> Result<(), String> {
            if workspace_restoring.get_untracked() || primer_context.get_untracked() != source {
                return Err(
                    "The source event changed; rerun primer analysis before navigating.".into(),
                );
            }
            let n = molecule.with_untracked(|m| m.sequence.len());
            if spans.is_empty()
                || spans.len() > 2
                || spans.iter().any(|(start, end)| start >= end || *end > n)
            {
                return Err("The primer placement has invalid source spans.".into());
            }
            batch(|| {
                selected_feature.set(None);
                selection.set(if spans.len() == 1 {
                    spans.first().copied()
                } else {
                    None
                });
                search_regions.set(Some((source.revision, spans)));
            });
            Ok(())
        },
    );
    // Index event identities on workspace changes, not full DNA on each oligo keystroke.
    let primer_event_index = Memo::new(move |_| {
        documents.with(|docs| {
            docs.iter()
                .flat_map(|doc| {
                    doc.summaries()
                        .into_iter()
                        .map(|event| (doc.id.clone(), event.event_id, event.revision))
                        .collect::<Vec<_>>()
                })
                .collect::<std::collections::BTreeSet<_>>()
        })
    });
    let mutation_selection = Callback::new(move |()| {
        if workspace_restoring.get_untracked() {
            return Err("Wait for the workspace operation to finish.".into());
        }
        documents.with_untracked(|docs| {
            let document = active.get_untracked();
            let source = &docs[document];
            let (start, end) = crate::selection_toolbar::mutation_range(
                &source.molecule,
                selection.get_untracked(),
                selected_feature.get_untracked(),
                &highlighted_spans.get_untracked(),
            )?;
            Ok(crate::mutagenesis::MutationSelection {
                source: current_event_binding(source),
                document,
                start,
                end,
            })
        })
    });
    let primer_anchor_available =
        Callback::new(move |binding: collagene_core::run_lineage::EventBinding| {
            primer_event_index.with(|events| {
                events.contains(&(binding.document_id, binding.event_id, binding.revision))
            })
        });

    let second_input = RwSignal::new(0usize);
    let first_input = RwSignal::new(0usize);
    let overlap = RwSignal::new("20".to_string());
    let closure_overlap = RwSignal::new("20".to_string());
    let candidate = RwSignal::new(None::<transfer::Candidate>);
    Effect::new(move |_| {
        if let Some(value) = candidate.get() {
            stored_candidates.update(|values| {
                if !values
                    .iter()
                    .any(|old| old.plan_revision == value.plan_revision)
                {
                    values.push(value);
                }
            });
        }
    });
    let plan_method = RwSignal::new("overlap");
    let first_enzymes = RwSignal::new("EcoRI".to_owned());
    let second_enzymes = RwSignal::new("EcoRI".to_owned());
    let first_fragment = RwSignal::new("1".to_owned());
    let second_fragment = RwSignal::new("1".to_owned());
    let first_treatment = RwSignal::new("none".to_owned());
    let second_treatment = RwSignal::new("none".to_owned());
    let first_reverse = RwSignal::new(false);
    let second_reverse = RwSignal::new(false);
    let restriction_policy = RwSignal::new("legacy");
    let restriction_final = RwSignal::new("ligate");
    let second_restriction = RwSignal::new(true);
    let first_repair = RwSignal::new("none".to_owned());
    let second_repair = RwSignal::new("none".to_owned());
    let first_modifications = crate::checked_cloning::AssumptionSignals::new();
    let second_modifications = crate::checked_cloning::AssumptionSignals::new();
    let last_plan_events =
        RwSignal::new(None::<(String, Vec<collagene_core::run_lineage::InputEventBinding>)>);
    let candidate_events = RwSignal::new(
        None::<(
            String,
            u64,
            collagene_core::run_lineage::EventBinding,
            Vec<collagene_core::run_lineage::InputEventBinding>,
        )>,
    );

    let allow_nicked = RwSignal::new(false);
    let plan_report = RwSignal::new(None::<VerificationReport>);
    let last_plan = RwSignal::new(None::<(DesignPlan, Molecule)>);
    let last_plan_target =
        RwSignal::new(None::<(String, u64, collagene_core::run_lineage::EventBinding)>);
    let capture_junction_plan = Callback::new(
        move |()| -> Result<(DesignPlan, Molecule), String> {
            if workspace_restoring.get_untracked() {
                return Err("Wait for workspace restoration.".into());
            }
            let (bound_plan, epoch, binding) = last_plan_target
                .get_untracked()
                .ok_or("Reconstruct against the exact intended-target event first.")?;
            if epoch != workspace_epoch.get_untracked()
                || !documents.with_untracked(|docs| {
                    docs.get(active.get_untracked())
                        .is_some_and(|doc| current_event_binding(doc) == binding)
                })
            {
                return Err(
                    "The intended-target event changed; select it and reconstruct again.".into(),
                );
            }
            last_plan.with_untracked(|last| {
            let (plan,target)=last.as_ref().ok_or("Reconstruct the current cloning plan before inspecting junctions.")?;
            collagene_core::assembly_junctions::source_preflight(plan,target).map_err(|e|e.to_string())?;
            let revision=target.revision();
            if !documents.with_untracked(|docs|docs.get(active.get_untracked()).is_some_and(|doc|doc.head_revision()==revision)) {
                return Err("The independently selected target changed; reconstruct the plan again.".into());
            }
            let expected=verification::plan_revision(plan);
            if expected!=bound_plan {return Err("The reviewed plan changed; reconstruct again.".into());}
            if !plan_report.with_untracked(|report|report.as_ref().is_some_and(|r|r.plan_revision==expected&&r.target_revision==revision&&r.engine_version==verification::ENGINE_VERSION)) {
                return Err("The current reconstruction changed; verify the plan again before inspecting junctions.".into());
            }
            Ok((plan.clone(),target.clone()))
        })
        },
    );
    let plan_generation = RwSignal::new(0u32);
    Effect::new(move |_| {
        let _ = (
            first_input.get(),
            second_input.get(),
            overlap.get(),
            closure_overlap.get(),
        );
        let _ = (
            plan_method.get(),
            first_enzymes.get(),
            second_enzymes.get(),
            first_fragment.get(),
            second_fragment.get(),
            first_treatment.get(),
            second_treatment.get(),
            first_reverse.get(),
            second_reverse.get(),
            allow_nicked.get(),
            restriction_policy.get(),
            restriction_final.get(),
            second_restriction.get(),
            first_repair.get(),
            second_repair.get(),
        );
        first_modifications.track();
        second_modifications.track();
        plan_generation.set(plan_generation.get_untracked().wrapping_add(1));
        plan_report.set(None);
        last_plan.set(None);
        last_plan_target.set(None);
    });
    let apply_restored = Callback::new(move |state: RestoredWorkspace| {
        workspace_epoch.update(|epoch| *epoch = epoch.wrapping_add(1));
        pcr_generation.update(|g| *g = g.wrapping_add(1));
        existing_oligo.set(None);
        pending_oligo.set(None);
        pending_fastq.set(None);
        pending_read_alignment.set(None);
        pending_sbol.set(None);
        pending_variants.set(None);
        pending_finalize.set(None);
        pending_release.set(None);
        selected_release_manifest.set(None);
        pending_run_event.set(None);
        pending_comparison_source.set(None);
        internal_clipboard.set(None);
        pending_reference_paste.set(None);
        read_reference.set(None);
        pending_read_check.set(None);
        reviewed_check.set(None);
        review_seed.set(None);
        pending_read_review.set(None);
        pair_draft.set(state.pair_draft);
        primer_seed.set(None);
        mutation_seed.set(None);
        hidden_products.set(state.hidden_products);
        show_prepared_history.set(false);
        documents.set(state.documents);
        active.set(state.active);
        workspace_name.set(state.name);
        stored_candidates.set(state.candidates);
        attachments.set(state.attachments);
        oligo_tables.set(state.oligos);
        libraries.set(state.libraries);
        alignments.set(state.alignments);
        recorded_runs.set(state.runs);
        saved_assemblies.set(state.assemblies);
        assembly_page.set(0);
        pending_assembly.set(None);
        pending_donor.set(None);
        pending_digest.set(None);
        pending_material.set(None);
        inspected_assembly.set(String::new());
        run_page.set(0);
        pending_alignment.set(None);
        annotation_source.set(String::new());
        method_config.set(crate::cloning_methods::MethodsState::default());
        assembly_draft.set(state.assembly_draft);
        assembly_visible.set(false);
        digest_boundary_pick.set(None);
        candidate_events.set(None);
        candidate.set(state.candidate);
        plan_generation.update(|g| *g = g.wrapping_add(1));
        plan_report.set(None);
        last_plan.set(None);
        last_plan_target.set(None);
        selection.set(None);
        selected_feature.set(None);
        modal.set("");
    });
    // Persistence callbacks outlive the dock/modal that initiated a transaction.
    // A committed write must reach App even after that child owner is disposed.
    let validate_workspace_frame = Callback::new(|frame: collagene_core::workspace::Workspace| {
        restore_workspace_state(&frame).map(|_| ())
    });
    let restore_workspace_frame =
        Callback::new(move |frame: collagene_core::workspace::Workspace| {
            match restore_workspace_state(&frame) {
                Ok(state) => batch(|| apply_restored.run(state)),
                Err(error) => message.set(error),
            }
        });
    let restore_library_frame =
        Callback::new(move |frame: collagene_core::workspace::Workspace| {
            match restore_workspace_state(&frame) {
                Ok(state) => {
                    let visible = assembly_visible.get_untracked();
                    batch(|| {
                        apply_restored.run(state);
                        assembly_visible.set(visible);
                    });
                }
                Err(error) => message.set(error),
            }
        });
    let library_current = Memo::new(move |_| {
        let epoch = workspace_epoch.get();
        documents.with(|docs| {
            let document = docs.get(active.get())?;
            Some(crate::library_view::LocalDocument {
                id: document.id.clone(),
                name: document.molecule.name.clone(),
                revision: collagene_library_protocol::Digest::new(document.head_revision()).ok()?,
                head_event: Some(document.head_event().into()),
                workspace_epoch: epoch,
            })
        })
    });
    let library_bindings = Memo::new(move |_| {
        attachments.with(|assets| crate::library_documents::bindings(assets).unwrap_or_default())
    });
    let native_catalog_refresh = RwSignal::new(0u64);
    // Native activation status must be checked before enabling saved-library
    // writes or synchronization. Working-editor navigation/autosave stays live.
    let native_catalog_pending = RwSignal::new(crate::desktop_rpc::is_desktop());
    let native_catalog_locked =
        Memo::new(move |_| workspace_restoring.get() || native_catalog_pending.get());
    let native_destination = RwSignal::new(None::<collagene_library_protocol::FolderId>);
    let native_catalog_epoch = RwSignal::new(0u64);
    let catalog_pending_changed = Callback::new(move |pending: bool| {
        let token = native_catalog_epoch.get_untracked().wrapping_add(1);
        native_catalog_epoch.set(token);
        if pending {
            native_catalog_pending.set(true);
            crate::library_client::browser().disconnect();
            native_library_context.set(None);
            native_destination.set(None);
        } else if native_library_context.get_untracked().is_some() {
            native_catalog_pending.set(false);
            native_catalog_refresh.update(|generation| *generation = generation.wrapping_add(1));
        } else {
            // Rebind the actual native catalog after an unknown/failed adoption,
            // never revive a cached browser context just from a UI status flag.
            native_catalog_pending.set(true);
            spawn_local(async move {
                let bootstrap =
                    crate::desktop_rpc::call(collagene_desktop_protocol::Command::Bootstrap).await;
                if native_catalog_epoch.try_get_untracked() != Some(token) {
                    return;
                }
                if let Ok(collagene_desktop_protocol::Success::Bootstrap(value)) = bootstrap {
                    if matches!(
                        value.catalog,
                        collagene_desktop_protocol::CatalogStatus::Ready { .. }
                    ) {
                        if let Some(context) = value.library_context {
                            native_library_context.set(Some(context));
                            native_catalog_pending.set(false);
                            native_catalog_refresh
                                .update(|generation| *generation = generation.wrapping_add(1));
                        }
                    }
                }
            });
        }
    });
    let native_destination_memo = Memo::new(move |_| native_destination.get());
    let native_folder_changed = Callback::new(move |folder| native_destination.set(folder));
    let native_catalog_changed = Callback::new(move |()| {
        native_catalog_refresh.update(|generation| *generation = generation.wrapping_add(1));
    });
    let native_saved = Memo::new(move |_| {
        let _ = loaded.get();
        let _ = native_catalog_refresh.get();
        let context = native_library_context.get()?;
        let current = library_current.get()?;
        let remote = crate::library_bindings::Remote {
            server: "collagene-desktop://local".into(),
            lab: context.lab_id,
        };
        let binding =
            library_bindings.with(|bindings| bindings.get(&remote, &current.id).cloned())?;
        Some(crate::remote_library_view::SavedLocalDocument {
            id: binding.remote_document,
            revision: binding.remote_revision,
            name: current.name,
            dirty: binding.local_revision != current.revision
                || binding.local_head_event != current.head_event,
        })
    });
    let library_snapshot = Callback::new(move |()| {
        documents.with_untracked(|docs| {
            workspace_snapshot(
                docs,
                active.get_untracked(),
                candidate.get_untracked().as_ref(),
                &stored_candidates.get_untracked(),
                &attachments.get_untracked(),
                &workspace_name.get_untracked(),
            )
        })
    });
    let library_capture = Callback::new(
        move |request: crate::library_documents::CaptureRequest| -> crate::library_view::LocalResult<crate::library_documents::CaptureOutcome> {
            let prepared = (|| {
                if workspace_restoring.get_untracked() || !loaded.get_untracked() {
                    return Err("Wait for local workspace storage to finish opening.".to_string());
                }
                Ok((library_snapshot.run(())?, workspace_epoch.get_untracked()))
            })();
            Box::pin(async move {
                let (frame, epoch) = prepared?;
                let captured = if crate::desktop_rpc::is_desktop() {
                    crate::library_documents::capture_private_local(frame, epoch, request).await?
                } else {
                    crate::library_documents::capture(frame, epoch, request).await?
                };
                if workspace_epoch.try_get_untracked() != Some(epoch)
                    || workspace_restoring.try_get_untracked().unwrap_or(true)
                {
                    return Err(
                        "Workspace changed while preparing the shared plasmid; nothing was sent."
                            .into(),
                    );
                }
                Ok(captured)
            })
        },
    );
    let library_record_binding = Callback::new(
        move |(epoch, binding): (u64, crate::library_bindings::Binding)| {
            if workspace_epoch.get_untracked() != epoch || workspace_restoring.get_untracked() {
                return Err("The server saved its revision, but the local workspace changed; no local link was overwritten.".to_string());
            }
            let retained = documents.with_untracked(|docs| {
                docs.iter()
                    .find(|document| document.id == binding.local_document)
                    .is_some_and(|document| {
                        document.summaries().iter().any(|event| {
                            event.revision == binding.local_revision.as_str()
                                && binding.local_head_event.as_deref()
                                    == Some(event.event_id.as_str())
                        })
                    })
            });
            if !retained {
                return Err("The remote revision is saved, but its captured local history is no longer present; reopen it from the lab library.".into());
            }
            let unchanged = attachments.with_untracked(|assets| {
                crate::library_documents::bindings(assets).map(|links| {
                    links.get(&binding.remote, &binding.local_document) == Some(&binding)
                })
            })?;
            if unchanged {
                return Ok(());
            }
            attachments
                .try_update(|assets| crate::library_documents::record_binding(assets, binding))
                .ok_or_else(|| {
                    "Workspace closed before its remote link could be recorded.".to_string()
                })?
        },
    );
    let library_open = Callback::new(move |opened:crate::library_view::Opened| -> crate::library_view::LocalResult<crate::library_view::LocalDocument> {
        let prepared = (|| {
            if workspace_restoring.get_untracked() || !loaded.get_untracked() || !opened.guard.is_current() {
                return Err("Wait for the current workspace operation to finish, then open the library revision again.".to_string());
            }
            let database=db.get_value().ok_or("Local browser storage is unavailable; the library revision was not imported.")?;
            let backup=library_snapshot.run(())?;
            let fingerprint=collagene_core::run_lineage::workspace_fingerprint(&backup).map_err(|error|error.to_string())?;
            let epoch=workspace_epoch.get_untracked();
            workspace_restoring.set(true);
            save_generation.update(|generation|*generation=generation.wrapping_add(1));
            Ok((database,backup,fingerprint,epoch))
        })();
        Box::pin(async move {
            let (database,backup,fingerprint,epoch)=prepared?;
            let result=async {
                let decoded=crate::library_documents::decode(&opened).await?;
                let guard=||->Result<(),String>{
                    if !opened.guard.is_current() || workspace_epoch.get_untracked()!=epoch
                        || collagene_core::run_lineage::workspace_fingerprint(&library_snapshot.run(())?).map_err(|error|error.to_string())?!=fingerprint
                    {return Err("Local workspace or library session changed; the import was not committed.".into());}
                    Ok(())
                };
                guard()?;
                let mut next=crate::library_documents::merge(&backup,&decoded,&opened)?;
                let mut state=restore_workspace_state(&next)?;
                let document=state.documents.get(state.active).ok_or("Imported document is unavailable")?;
                let mut info=crate::library_view::LocalDocument{
                    id:document.id.clone(),name:document.molecule.name.clone(),
                    revision:collagene_library_protocol::Digest::new(document.head_revision()).map_err(|error|error.to_string())?,
                    head_event:Some(document.head_event().into()),workspace_epoch:epoch.wrapping_add(1),
                };
                crate::library_documents::record_binding(&mut next.attachments,crate::library_bindings::Binding{
                    remote:opened.remote.clone(),local_document:info.id.clone(),local_revision:info.revision.clone(),
                    local_head_event:info.head_event.clone(),remote_document:opened.document.id.clone(),remote_revision:opened.revision.id.clone(),
                })?;
                state.attachments=next.attachments.clone();
                let current=encode_workspace(next,u32::MAX-303).await?;
                let previous=encode_workspace(backup,u32::MAX-304).await?;
                guard()?;
                crate::storage::replace_workspace(&database,&current,&previous).await?;
                // A committed storage transaction reaches App even if its sidebar
                // has closed; a late view cancellation must not undo durable work.
                let visible=assembly_visible.get_untracked();
                batch(||{apply_restored.run(state);assembly_visible.set(visible);});
                info.workspace_epoch=workspace_epoch.get_untracked();
                Ok(info)
            }.await;
            workspace_restoring.set(false);
            result
        })
    });
    Effect::new(move |_| {
        let draft = assembly_draft.get();
        if workspace_restoring.get() {
            return;
        }
        match crate::assembly_draft_model::encode_raw(&draft) {
            Ok(bytes) => {
                let draft_kind = crate::assembly_draft_model::ATTACHMENT_KIND;
                let hash = blake3::hash(&bytes).to_hex().to_string();
                if attachments.with_untracked(|assets| {
                    assets.get(&hash).is_some_and(|a| a.kind == draft_kind)
                        && assets
                            .values()
                            .filter(|a| assembly_draft_kind(&a.kind))
                            .count()
                            == 1
                }) {
                    return;
                }
                attachments.update(|assets| {
                    assets.retain(|_, a| !assembly_draft_kind(&a.kind));
                    assets.insert(
                        hash,
                        collagene_core::workspace::WorkspaceAttachment {
                            kind: draft_kind.into(),
                            name: "assembly-draft.json".into(),
                            media_type: "application/json".into(),
                            bytes,
                        },
                    );
                });
            }
            Err(error) => {
                message.set(format!("Assembly draft could not be saved: {error}"));
            }
        }
    });
    Effect::new(move |_| {
        let draft = pair_draft.get();
        if workspace_restoring.get() {
            return;
        }
        match crate::primer_pair_storage::encode(&draft) {
            Ok(bytes) => {
                let kind = crate::primer_pair_storage::ATTACHMENT_KIND;
                let hash = blake3::hash(&bytes).to_hex().to_string();
                if attachments
                    .with_untracked(|assets| assets.get(&hash).is_some_and(|a| a.kind == kind))
                {
                    return;
                }
                let is_default = draft == crate::primer_pair_model::PairDraft::default();
                if is_default
                    && !attachments.with_untracked(|assets| assets.values().any(|a| a.kind == kind))
                {
                    return;
                }
                attachments.update(|assets| {
                    assets.retain(|_, a| a.kind != kind);
                    if !is_default {
                        assets.insert(
                            hash,
                            collagene_core::workspace::WorkspaceAttachment {
                                kind: kind.into(),
                                name: "primer-pair-draft.json".into(),
                                media_type: "application/json".into(),
                                bytes,
                            },
                        );
                    }
                });
            }
            Err(error) => message.set(format!("Primer pair draft could not be saved: {error}")),
        }
    });
    spawn_local(async move {
        match storage::open().await {
            Ok(database) => {
                native_library_context.set(database.native_context());
                db.set_value(Some(database.clone()));
                match storage::load(&database).await {
                    Ok(Some(bytes)) => {
                        let frame = decode_workspace(bytes, 0).await;
                        match frame.and_then(|frame| restore_workspace_state(&frame)) {
                            Ok(state) => apply_restored.run(state),
                            Err(error) => {
                                message.set(format!("Autosave could not be read. The stored record has been left intact; open a saved workspace to recover. {error}"));
                                save_status.set("Recovery needed".into());
                                loaded.set(false);
                                workspace_restoring.set(false);
                                return;
                            }
                        }
                    }
                    Ok(None) => {}
                    Err(error) => {
                        message.set(error);
                        save_status.set("Autosave unavailable".into());
                        workspace_restoring.set(false);
                        return;
                    }
                }
                loaded.set(true);
                workspace_restoring.set(false);
            }
            Err(error) => {
                message.set(format!(
                    "{error}. Export a workspace file to keep your work."
                ));
                save_status.set("Autosave unavailable".into());
                workspace_restoring.set(false);
            }
        }
    });
    Effect::new(move |_| {
        if !loaded.get() || workspace_restoring.get() {
            return;
        }
        let index = active.get();
        let name = workspace_name.get();
        let frame = documents.with(|docs| {
            stored_candidates.with(|plans| {
                attachments.with(|assets| {
                    candidate.with(|current| {
                        workspace_snapshot(docs, index, current.as_ref(), plans, assets, &name)
                    })
                })
            })
        });
        let generation = save_generation.get_untracked().wrapping_add(1);
        save_generation.set(generation);
        save_status.set("Saving locally…".into());
        let frame = match frame {
            Ok(frame) => frame,
            Err(error) => {
                save_status.set("Autosave failed".into());
                message.set(error);
                return;
            }
        };
        if let Some(database) = db.get_value() {
            spawn_local(async move {
                let result = encode_workspace(frame, generation as u32).await;
                // Work can finish out of order. Only the current generation may start a
                // write transaction, including after native admission waits.
                if generation != save_generation.get_untracked()
                    || workspace_restoring.get_untracked()
                {
                    return;
                }
                let result = match result {
                    Ok(bytes) => {
                        storage::save(&database, &bytes, move || {
                            if generation != save_generation.get_untracked()
                                || workspace_restoring.get_untracked()
                            {
                                Err("Autosave superseded".into())
                            } else {
                                Ok(())
                            }
                        })
                        .await
                    }
                    Err(error) => Err(error),
                };
                if generation != save_generation.get_untracked()
                    || workspace_restoring.get_untracked()
                {
                    return;
                }
                match result {
                    Ok(()) => save_status.set(
                        if crate::desktop_rpc::is_desktop() {
                            "Saved locally"
                        } else {
                            "Saved in this browser"
                        }
                        .into(),
                    ),
                    Err(error) => {
                        save_status.set("Autosave failed".into());
                        message.set(error);
                    }
                }
            });
        }
    });
    Effect::new(move |_| {
        if !has_open_document.get() {
            if let Some(map) = renderer.get_value() {
                map.pointer_cancel();
            }
            return;
        }
        molecule.with(|current| {
            let selected = highlighted_spans.get();
            let selected_index = selected_feature.get();
            if let (Some(canvas), Some(overlay)) = (gl_canvas.get(), label_canvas.get()) {
                if let Some(map) = renderer.get_value() {
                    map.set_document_spans(
                        current,
                        &molecule_revision.get(),
                        &selected,
                        selected_index,
                    );
                    map.set_hidden_features(&annotation_hidden.get());
                } else {
                    match MapRenderer::new(
                        canvas.unchecked_ref::<HtmlCanvasElement>().clone(),
                        overlay.unchecked_ref::<HtmlCanvasElement>().clone(),
                        current.clone(),
                        std::rc::Rc::new(move |stats| {
                            if map_labels
                                .try_get_untracked()
                                .is_some_and(|current| current != stats)
                            {
                                let _ = map_labels.try_set(stats);
                            }
                        }),
                    ) {
                        Ok(map) => {
                            map.set_restriction_hover_callback(std::rc::Rc::new(move |site| {
                                if hovered_restriction.try_get_untracked() != Some(site) {
                                    let _ = hovered_restriction.try_set(site);
                                }
                            }));
                            map.set_selection_mode(false);
                            map.set_hidden_features(&annotation_hidden.get());
                            renderer.set_value(Some(map));
                        }
                        Err(error) => render_error.set(error),
                    }
                }
            }
        });
    });
    Effect::new(move |_| {
        let _ = restriction_identity.get();
        hovered_restriction.set(None);
    });
    Effect::new(move |_| {
        let context = restriction_identity.get();
        let sites = restrictions
            .displayed
            .get()
            .unwrap_or_else(|| std::sync::Arc::new(Vec::new()));
        let selected = hovered_restriction.get().or(restrictions.selected.get());
        let revision = molecule_revision.get();
        if let Some(map) = renderer.get_value() {
            map.set_interaction_context(&context);
            map.set_restriction_sites(&revision, sites, selected);
            map.set_pcr_overlay(pcr_overlay.get());
        }
    });
    Effect::new(move |_| {
        // All document writes still invalidate multi-source scientific reports.
        // Selection follows only the exact active document/event and restore epoch.
        let epoch = workspace_epoch.get();
        let index = active.get();
        let pin = documents.with(|docs| {
            docs.get(index)
                .map(|doc| crate::sequence_edit_transition::SourcePin::of(epoch, doc))
        });
        let mut decision = crate::sequence_edit_transition::Decision::Unchanged;
        sequence_transitions.update_value(|tracker| decision = tracker.observe(pin));
        match decision {
            crate::sequence_edit_transition::Decision::Unchanged => {}
            crate::sequence_edit_transition::Decision::Clear => {
                selected_feature.set(None);
                selection.set(None);
                search_regions.set(None);
            }
            crate::sequence_edit_transition::Decision::PreserveCaret(caret) => {
                // Read untracked: caret changes must not retrigger source processing.
                // A later pointer selection takes precedence over the local edit.
                if selection.get_untracked() == Some((caret, caret)) {
                    selected_feature.set(None);
                    search_regions.set(None);
                    selection.set(Some((caret, caret)));
                }
            }
        }
        plan_generation.set(plan_generation.get_untracked().wrapping_add(1));
        plan_report.set(None);
        last_plan.set(None);
        last_plan_target.set(None);
    });

    let open_file = move |event: web_sys::Event| {
        let input = event_target::<HtmlInputElement>(&event);
        let Some(files) = input.files() else {
            return;
        };
        let files: Vec<_> = (0..files.length()).filter_map(|i| files.get(i)).collect();
        if workspace_restoring.get_untracked() {
            input.set_value("");
            return;
        }
        let opened_epoch = workspace_epoch.get_untracked();
        spawn_local(async move {
            let mut imported = vec![];
            let mut warnings = vec![];
            let file_count = files.len();
            let mut requested_target = None::<String>;
            for file in files {
                if file.size() > (formats::MAX_IMPORT_BYTES + 42) as f64 {
                    message.set("File exceeds the 32 MiB payload import limit.".into());
                    return;
                }
                let bytes = match JsFuture::from(file.array_buffer()).await {
                    Ok(buffer) => js_sys::Uint8Array::new(&buffer).to_vec(),
                    Err(_) => {
                        message.set("Could not read file.".into());
                        return;
                    }
                };
                if workspace_restoring.get_untracked()
                    || workspace_epoch.get_untracked() != opened_epoch
                {
                    message.set("The workspace changed while reading files. Reopen the files to import them here.".into());
                    return;
                }
                let extension = file
                    .name()
                    .rsplit('.')
                    .next()
                    .unwrap_or("")
                    .to_ascii_lowercase();
                let lower_name = file.name().to_ascii_lowercase();
                if extension == "sbol"
                    || ((extension == "jsonld" || extension == "json")
                        && crate::sbol_exchange::looks_like_sbol(&bytes))
                {
                    if file_count != 1 {
                        message.set("Open an SBOL graph separately to review its component and sequence choices; no files were imported.".into());
                        return;
                    }
                    if bytes.len() > collagene_design_exchange::MAX_BYTES {
                        message.set("SBOL source exceeds the 16 MiB exchange bound.".into());
                        return;
                    }
                    pending_sbol.set(Some(std::sync::Arc::new(bytes)));
                    modal_error.set(String::new());
                    modal.set("sbol");
                    return;
                }

                let vcf_source = lower_name.ends_with(".vcf")
                    || lower_name.ends_with(".vcf.gz")
                    || bytes.starts_with(b"##fileformat=VCF");
                let bcf_source = lower_name.ends_with(".bcf") || lower_name.ends_with(".bcf.gz");
                if vcf_source || bcf_source {
                    if file_count != 1 {
                        message.set("Open a VCF/BCF source separately to inspect its declared alleles; no files were imported.".into());
                        return;
                    }
                    if bytes.len() > 16 * 1024 * 1024 {
                        message.set(
                            "This browser variant profile accepts source files up to 16 MiB."
                                .into(),
                        );
                        return;
                    }
                    let format = if bcf_source {
                        collagene_core::variants_session::FileFormat::Bcf
                    } else {
                        collagene_core::variants_session::FileFormat::Vcf
                    };
                    pending_variants.set(Some((std::sync::Arc::new(bytes), format, file.name())));
                    modal_error.set(String::new());
                    modal.set("variants");
                    return;
                }
                let sam_header = [
                    b"@HD\t".as_slice(),
                    b"@SQ\t".as_slice(),
                    b"@RG\t".as_slice(),
                    b"@PG\t".as_slice(),
                    b"@CO\t".as_slice(),
                ]
                .iter()
                .any(|prefix| bytes.starts_with(prefix));
                if matches!(extension.as_str(), "sam" | "bam")
                    || (sam_header && !matches!(extension.as_str(), "fastq" | "fq"))
                {
                    if file_count != 1 {
                        message.set("Open a SAM/BAM source separately to inspect its records and explicitly bind DNA references; no files were imported.".into());
                        return;
                    }
                    if bytes.len() > 16 * 1024 * 1024 {
                        message.set("This browser read-alignment profile accepts source files up to 16 MiB.".into());
                        return;
                    }
                    let format = if extension == "bam" {
                        collagene_core::read_alignment_session::FileFormat::Bam
                    } else {
                        collagene_core::read_alignment_session::FileFormat::Sam
                    };
                    pending_read_alignment.set(Some((
                        std::sync::Arc::new(bytes),
                        format,
                        file.name(),
                    )));
                    modal_error.set(String::new());
                    modal.set("read_alignments");
                    return;
                }
                if bytes.first() == Some(&b'@') {
                    if file_count != 1 {
                        message.set("Open FASTQ separately to choose its supplied quality encoding; no files were imported.".into());
                        return;
                    }
                    pending_fastq.set(Some(std::sync::Arc::new(bytes)));
                    fastq_encoding.set(None);
                    modal_error.set(String::new());
                    modal.set("fastq");
                    return;
                }
                if matches!(
                    extension.as_str(),
                    "aln"
                        | "clustal"
                        | "sto"
                        | "stockholm"
                        | "afa"
                        | "afasta"
                        | "cg2aln"
                        | "phy"
                        | "phylip"
                        | "msf"
                        | "nex"
                        | "nexus"
                        | "pir"
                        | "nbrf"
                        | "gde"
                        | "slx"
                        | "selex"
                        | "sel"
                ) || bytes.starts_with(b"CLUSTAL")
                    || bytes.starts_with(b"# STOCKHOLM")
                    || bytes.starts_with(b"!!AA_MULTIPLE_ALIGNMENT")
                    || bytes.starts_with(b"!!NA_MULTIPLE_ALIGNMENT")
                    || [b">P1;", b">F1;", b">DL;", b">RL;", b">DC;", b">RC;"]
                        .iter()
                        .any(|header| bytes.starts_with(*header))
                    || bytes
                        .get(..6)
                        .is_some_and(|header| header.eq_ignore_ascii_case(b"#NEXUS"))
                {
                    if file_count != 1 {
                        message.set("Open an alignment separately so its alphabet can be chosen explicitly; no files were imported.".into());
                        return;
                    }
                    let format = if bytes.starts_with(b"CLUSTAL")
                        || matches!(extension.as_str(), "aln" | "clustal")
                    {
                        Some(collagene_core::alignment_formats::AlignmentFormat::Clustal)
                    } else if bytes.starts_with(b"# STOCKHOLM")
                        || matches!(extension.as_str(), "sto" | "stockholm")
                    {
                        Some(collagene_core::alignment_formats::AlignmentFormat::Stockholm)
                    } else if matches!(extension.as_str(), "afa" | "afasta") {
                        Some(collagene_core::alignment_formats::AlignmentFormat::AlignedFasta)
                    } else if extension == "msf"
                        || bytes.starts_with(b"!!AA_MULTIPLE_ALIGNMENT")
                        || bytes.starts_with(b"!!NA_MULTIPLE_ALIGNMENT")
                    {
                        Some(collagene_core::alignment_formats::AlignmentFormat::Msf)
                    } else if matches!(extension.as_str(), "nex" | "nexus")
                        || bytes
                            .get(..6)
                            .is_some_and(|header| header.eq_ignore_ascii_case(b"#NEXUS"))
                    {
                        Some(collagene_core::alignment_formats::AlignmentFormat::Nexus)
                    } else if matches!(extension.as_str(), "pir" | "nbrf")
                        || [b">P1;", b">F1;", b">DL;", b">RL;", b">DC;", b">RC;"]
                            .iter()
                            .any(|header| bytes.starts_with(*header))
                    {
                        Some(collagene_core::alignment_formats::AlignmentFormat::Pir)
                    } else if extension == "gde" {
                        Some(collagene_core::alignment_formats::AlignmentFormat::Gde)
                    } else if matches!(extension.as_str(), "slx" | "selex" | "sel") {
                        Some(collagene_core::alignment_formats::AlignmentFormat::Selex)
                    } else {
                        None
                    };
                    pending_alignment.set(Some(crate::alignment_view::PendingAlignment {
                        name: file.name(),
                        bytes: std::sync::Arc::new(bytes),
                        format,
                        native: extension == "cg2aln",
                    }));
                    modal_error.set(String::new());
                    modal.set("alignments");
                    return;
                }
                if matches!(extension.as_str(), "gff3" | "gff" | "gtf" | "bed") {
                    if bytes.len() > collagene_core::annotation_formats::MAX_ANNOTATION_BYTES {
                        message.set("Annotation source exceeds the 8 MiB import limit.".into());
                        return;
                    }
                    if file_count != 1 {
                        message.set("Open an annotation file separately so its reference mapping can be reviewed; no files were imported.".into());
                        return;
                    }
                    match String::from_utf8(bytes) {
                        Ok(text) => {
                            annotation_source.set(text);
                            annotation_format.set(match extension.as_str() {
                                "bed" => collagene_core::annotation_formats::AnnotationFormat::Bed,
                                "gtf" => collagene_core::annotation_formats::AnnotationFormat::Gtf,
                                _ => collagene_core::annotation_formats::AnnotationFormat::Gff3,
                            });
                            panel.set("annotation_files");
                            message.set("Choose the exact annotation reference ID and review the proposed coordinates before applying.".into());
                        }
                        Err(_) => message.set("Annotation files must be UTF-8 text.".into()),
                    }
                    return;
                }
                if extension == "csv" || extension == "tsv" {
                    if file_count != 1 {
                        message.set("Open an oligo table separately to map its columns; no files were imported.".into());
                        return;
                    }
                    match String::from_utf8(bytes) {
                        Ok(text) => {
                            pending_oligo.set(Some((
                                text,
                                if extension == "tsv" {
                                    oligo_formats::Delimiter::Tab
                                } else {
                                    oligo_formats::Delimiter::Comma
                                },
                            )));
                            existing_oligo.set(None);
                            modal_error.set(String::new());
                            modal.set("oligos");
                        }
                        Err(_) => message.set("Oligo tables must be UTF-8 text.".into()),
                    }
                    return;
                }
                if bytes.starts_with(b"CGP\0") {
                    if file_count != 1 {
                        message.set("Open a candidate separately from sequence files.".into());
                        return;
                    }
                    match transfer::load_candidate(&bytes) {
                        Ok(imported_candidate) => {
                            candidate_events.set(None);
                            candidate.set(Some(imported_candidate));
                            plan_generation.update(|g| *g = g.wrapping_add(1));
                            plan_report.set(None);
                            last_plan.set(None);
                            last_plan_target.set(None);
                            panel.set("plan");
                            message.set("Candidate opened for review. Select its intended target and verify against your open documents.".into());
                        }
                        Err(e) => message.set(e.to_string()),
                    }
                    return;
                } else if bytes.starts_with(b"CGD\0") {
                    match transfer::load_workspace(&bytes) {
                        Ok(snapshot) => {
                            requested_target = Some(snapshot.target_revision);
                            imported.extend(snapshot.documents);
                        }
                        Err(e) => {
                            message.set(e.to_string());
                            return;
                        }
                    }
                } else if bytes.starts_with(b"CG2\0") {
                    match protocol::load_document(&bytes) {
                        Ok(m) => imported.push(m),
                        Err(e) => {
                            message.set(e.to_string());
                            return;
                        }
                    }
                } else {
                    let compat_profile = if collagene_core::gcg_sequence::is_gcg(&bytes) {
                        Some(crate::sequence_compat_import::Profile::Gcg)
                    } else if std::str::from_utf8(&bytes).ok().is_some_and(|text| {
                        text.trim_start_matches('\u{feff}')
                            .trim_start()
                            .starts_with('>')
                    }) {
                        let hint = match extension.as_str() {
                            "faa" => Some(Alphabet::Protein),
                            "frn" => Some(Alphabet::Rna),
                            _ => None,
                        };
                        Some(crate::sequence_compat_import::Profile::Fasta(hint))
                    } else if collagene_core::legacy_formats::is_collagene_xml(&bytes)
                        && collagene_core::collagene_xml_compat::needs_lower_only_review(&bytes)
                            .unwrap_or(false)
                    {
                        Some(crate::sequence_compat_import::Profile::CollageneXml)
                    } else {
                        None
                    };
                    if let Some(profile) = compat_profile {
                        if file_count != 1 {
                            message.set("Open this sequence file separately to review its import profile; no files were imported.".into());
                            return;
                        }
                        pending_compat_import.set(Some((
                            std::sync::Arc::new(bytes),
                            file.name(),
                            profile,
                        )));
                        modal_error.set(String::new());
                        modal.set("compat_import");
                        return;
                    }
                    let report = if collagene_core::legacy_formats::is_xdna(&bytes) {
                        collagene_core::legacy_formats::import_xdna(&bytes, &file.name())
                    } else {
                        formats::import_bytes(&bytes)
                    };
                    match report {
                        Ok(report) => {
                            imported.extend(report.molecules);
                            warnings.extend(report.warnings);
                        }
                        Err(e) => {
                            message.set(format!("{}: {e}", file.name()));
                            return;
                        }
                    }
                }
                if imported.len() > 1000
                    || imported.iter().map(document_bytes).sum::<usize>() > MAX_WORKSPACE_BYTES
                {
                    message.set("Import exceeds this prototype's 1,000-document / 64 MiB workspace limit. Open a smaller selection of files.".into());
                    return;
                }
            }
            if let Err(error) = add_documents(documents, active, workspace_restoring, imported) {
                message.set(error);
                return;
            }
            if let Some(revision) = requested_target {
                if let Some(index) = documents.with_untracked(|docs| {
                    docs.iter().position(|d| d.molecule.revision() == revision)
                }) {
                    active.set(index);
                }
            }
            message.set(warnings.join(" "));
        });
        input.set_value("");
    };
    let save_file = move |_: web_sys::MouseEvent| {
        let mol = molecule.get_untracked();
        match protocol::save_document(&mol) {
            Ok(bytes) => crate::app::download_with_completion(
                &format!("{}.cg2", mol.name),
                "application/octet-stream",
                &bytes,
                move |result| {
                    if message.try_get_untracked().is_none() {
                        return;
                    }
                    if let Err(e) = result {
                        message.set(e)
                    }
                },
            ),
            Err(e) => message.set(e.to_string()),
        }
    };
    let do_pcr = move |_: web_sys::MouseEvent| {
        if [
            pcr_forward_name.get_untracked(),
            pcr_reverse_name.get_untracked(),
        ]
        .iter()
        .any(|name| {
            name.trim().is_empty() || name.len() > 256 || name.chars().any(char::is_control)
        }) {
            message.set("Primer names must be nonempty, at most 256 bytes, and contain no control characters.".into());
            return;
        }
        let primers = Primer::new(
            pcr_forward_name.get_untracked(),
            &pcr_forward.get_untracked(),
            &pcr_forward_tail.get_untracked(),
        )
        .and_then(|f| {
            Primer::new(
                pcr_reverse_name.get_untracked(),
                &pcr_reverse.get_untracked(),
                &pcr_reverse_tail.get_untracked(),
            )
            .map(|r| (f, r))
        });
        let (forward, reverse) = match primers {
            Ok(p) => p,
            Err(error) => {
                message.set(error.to_string());
                return;
            }
        };
        let mol = molecule.get_untracked();
        let source_revision = mol.revision();
        let origin = pcr_origin.get_untracked();
        let revision = match origin.as_ref().map_or_else(
            || Ok(source_revision.clone()),
            |origin| {
                collagene_core::material_pcr::revision(&mol, origin.as_ref(), &forward, &reverse)
                    .map_err(|e| e.to_string())
            },
        ) {
            Ok(revision) => revision,
            Err(error) => {
                message.set(error);
                return;
            }
        };
        let pcr_source =
            documents.with_untracked(|docs| current_event_binding(&docs[active.get_untracked()]));
        if let Some(origin) = &origin {
            if origin.event() != &pcr_source {
                message.set("PCR source changed; choose its material origin again.".into());
                return;
            }
            let checked = library_snapshot.run(()).and_then(|workspace| {
                origin
                    .as_ref()
                    .validate_for_workspace(&workspace, &pcr_source)
                    .map_err(|e| e.to_string())
            });
            if let Err(error) = checked {
                message.set(error);
                return;
            }
        }
        let pcr_epoch = workspace_epoch.get_untracked();
        let pcr_inputs = (
            pcr_forward_name.get_untracked(),
            pcr_reverse_name.get_untracked(),
            pcr_forward.get_untracked(),
            pcr_reverse.get_untracked(),
            pcr_forward_tail.get_untracked(),
            pcr_reverse_tail.get_untracked(),
        );
        let generation = pcr_generation.get_untracked().wrapping_add(1);
        pcr_generation.set(generation);
        pcr_review.set(None);
        pcr_output_busy.set(false);
        if let Some(worker) = pcr_output_worker.get_value() {
            worker.close();
        }
        pcr_output_worker.set_value(None);
        let reviewed_forward = forward.clone();
        let reviewed_reverse = reverse.clone();
        let request = protocol::Request {
            id: generation,
            revision: revision.clone(),
            job: if let Some(origin) = &origin {
                protocol::AnalysisJob::MaterialPcrDocuments {
                    molecule: mol,
                    origin: origin.as_ref().clone(),
                    forward,
                    reverse,
                }
            } else {
                protocol::AnalysisJob::PcrDocuments {
                    molecule: mol,
                    forward,
                    reverse,
                }
            },
        };
        message.set("Computing PCR products in background…".into());
        let alive = pcr_alive.get_value();
        spawn_local(async move {
            let result = worker::analyze_cancellable(request, alive.clone()).await;
            if !alive.load(std::sync::atomic::Ordering::Relaxed) {
                return;
            }
            if workspace_restoring.get_untracked()
                || workspace_epoch.get_untracked() != pcr_epoch
                || !documents.with_untracked(|docs| {
                    docs.get(active.get_untracked())
                        .is_some_and(|doc| current_event_binding(doc) == pcr_source)
                })
                || molecule.get_untracked().revision() != source_revision
                || pcr_origin.get_untracked() != origin
                || pcr_generation.get_untracked() != generation
                || pcr_inputs
                    != (
                        pcr_forward_name.get_untracked(),
                        pcr_reverse_name.get_untracked(),
                        pcr_forward.get_untracked(),
                        pcr_reverse.get_untracked(),
                        pcr_forward_tail.get_untracked(),
                        pcr_reverse_tail.get_untracked(),
                    )
            {
                return;
            }
            match (origin.is_some(), result) {
                (
                    false,
                    Ok(protocol::Response {
                        result: protocol::AnalysisResult::PcrDocuments(products),
                        id,
                        revision: response_revision,
                    }),
                )
                | (
                    true,
                    Ok(protocol::Response {
                        result: protocol::AnalysisResult::MaterialPcrDocuments(products),
                        id,
                        revision: response_revision,
                    }),
                ) if id == generation && response_revision == revision => {
                    let count = products.len();
                    pcr_review.set(Some(crate::pcr_outputs::Review {
                        source: pcr_source,
                        origin,
                        epoch: pcr_epoch,
                        generation,
                        forward: reviewed_forward,
                        reverse: reviewed_reverse,
                        products: std::sync::Arc::new(products),
                        inputs: pcr_inputs,
                    }));
                    message.set(format!("{count} exact-annealing PCR geometries. Review and choose outputs; no documents have been created."));
                }
                (
                    _,
                    Ok(protocol::Response {
                        result: protocol::AnalysisResult::Error(error),
                        ..
                    }),
                )
                | (_, Err(error)) => message.set(error),
                _ => message.set("Unexpected worker response".into()),
            }
        });
    };
    let pcr_review_current = Callback::new(move |review: crate::pcr_outputs::Review| {
        !workspace_restoring.get_untracked()
            && workspace_epoch.get_untracked() == review.epoch
            && pcr_generation.get_untracked() == review.generation
            && pcr_origin.get_untracked() == review.origin
            && review.inputs
                == (
                    pcr_forward_name.get_untracked(),
                    pcr_reverse_name.get_untracked(),
                    pcr_forward.get_untracked(),
                    pcr_reverse.get_untracked(),
                    pcr_forward_tail.get_untracked(),
                    pcr_reverse_tail.get_untracked(),
                )
            && documents.with_untracked(|docs| {
                docs.get(active.get_untracked())
                    .is_some_and(|doc| current_event_binding(doc) == review.source)
            })
            && pcr_review.with_untracked(|current| {
                current.as_ref().is_some_and(|current| {
                    current.generation == review.generation
                        && std::sync::Arc::ptr_eq(&current.products, &review.products)
                })
            })
    });
    let open_material_derivation = Callback::new(
        move |request: crate::digest_product::MaterialOpenRequest| -> Result<(), String> {
            if workspace_restoring.get_untracked()
                || workspace_epoch.get_untracked() != request.epoch
            {
                return Err("The material source workspace changed; review again.".into());
            }
            request.guard.run(())?;
            let original = request.guard;
            let bound_epoch = request.epoch;
            let guarded = Callback::new(move |()| -> Result<(), String> {
                if workspace_restoring.get_untracked()
                    || workspace_epoch.get_untracked() != bound_epoch
                    || modal.get_untracked() != "material_product"
                {
                    return Err("The material installation context changed.".into());
                }
                original.run(())
            });
            pending_material.set(Some(crate::digest_product::MaterialOpenRequest {
                guard: guarded,
                ..request
            }));
            modal_error.set(String::new());
            modal.set("material_product");
            Ok(())
        },
    );
    let create_pcr_outputs = Callback::new(move |choice: crate::pcr_outputs::Choice| {
        if pcr_output_busy.get_untracked() || !pcr_review_current.run(choice.review.clone()) {
            message.set("The PCR draft or template changed; compute and review again.".into());
            return;
        }
        let Some(selected) = choice.review.products.get(choice.product_index) else {
            message.set("Choose a reported PCR geometry.".into());
            return;
        };
        let Some(source) = documents.with_untracked(|docs| {
            docs.get(active.get_untracked())
                .map(|doc| doc.molecule.clone())
        }) else {
            return;
        };
        if choice.save_product {
            let prepared = (|| -> Result<crate::digest_product::MaterialOpenRequest, String> {
                let pin = crate::material_derivation_model::SourcePin::capture(
                    choice.review.source.clone(),
                    choice.review.epoch,
                    &source,
                    choice
                        .review
                        .origin
                        .as_ref()
                        .map(|origin| std::sync::Arc::new(origin.as_ref().clone())),
                )?;
                let captured = crate::material_derivation_model::PcrChoice::capture(
                    pin,
                    choice.review.forward.clone(),
                    choice.review.reverse.clone(),
                    choice.product_index,
                    std::sync::Arc::new(selected.clone()),
                )?;
                let frame = library_snapshot.run(())?;
                let pending =
                    captured.request_with_primers(&frame, &choice.name, choice.save_primers)?;
                let bound = choice.review.clone();
                let guard = Callback::new(move |()| {
                    if modal.get() != "material_product"
                        || workspace_epoch.get() != bound.epoch
                        || pcr_origin.get() != bound.origin
                        || bound.inputs
                            != (
                                pcr_forward_name.get(),
                                pcr_reverse_name.get(),
                                pcr_forward.get(),
                                pcr_reverse.get(),
                                pcr_forward_tail.get(),
                                pcr_reverse_tail.get(),
                            )
                        || !documents.with(|docs| {
                            docs.get(active.get())
                                .is_some_and(|doc| current_event_binding(doc) == bound.source)
                        })
                    {
                        return Err("PCR source, primers or material origin changed; recompute and select again.".into());
                    }
                    Ok(())
                });
                Ok(crate::digest_product::MaterialOpenRequest {
                    request: pending.request().clone(),
                    epoch: choice.review.epoch,
                    guard,
                })
            })();
            match prepared {
                Ok(request) => {
                    pending_material.set(Some(request));
                    modal_error.set(String::new());
                    modal.set("material_product");
                }
                Err(error) => message.set(error),
            }
            return;
        }
        let request = collagene_core::pcr_outputs::Request {
            source,
            forward: choice.review.forward.clone(),
            reverse: choice.review.reverse.clone(),
            selected_index: choice.product_index,
            expected_product_revision: selected.molecule.revision(),
            options: collagene_core::pcr_outputs::Options {
                name: choice.name.clone(),
                save_product: choice.save_product,
                save_primers: choice.save_primers,
            },
        };
        let preview_revision = request.expected_product_revision.clone();
        let revision = match choice.review.origin.as_ref().map_or_else(
            || collagene_core::pcr_outputs::revision(&request).map_err(|e| e.to_string()),
            |origin| {
                collagene_core::material_pcr::outputs_revision(origin.as_ref(), &request)
                    .map_err(|e| e.to_string())
            },
        ) {
            Ok(revision) => revision,
            Err(error) => {
                message.set(error.to_string());
                return;
            }
        };
        let w = match worker::CachedWorker::new() {
            Ok(w) => std::rc::Rc::new(w),
            Err(error) => {
                message.set(error);
                return;
            }
        };
        pcr_output_worker.set_value(Some(w.clone()));
        pcr_output_busy.set(true);
        message.set("Reconstructing the selected PCR geometry before creating outputs…".into());
        let alive = pcr_alive.get_value();
        spawn_local(async move {
            let result = w
                .ask(protocol::Request {
                    id: choice.review.generation,
                    revision: revision.clone(),
                    job: if let Some(origin) = &choice.review.origin {
                        protocol::AnalysisJob::PrepareMaterialPcrOutputs {
                            origin: origin.as_ref().clone(),
                            request,
                        }
                    } else {
                        protocol::AnalysisJob::PreparePcrOutputs { request }
                    },
                })
                .await;
            w.close();
            if !alive.load(std::sync::atomic::Ordering::Relaxed)
                || !pcr_review_current.run(choice.review.clone())
            {
                return;
            }
            pcr_output_worker.set_value(None);
            pcr_output_busy.set(false);
            let outcome = (|| -> Result<(), String> {
                let response = result?;
                if response.id != choice.review.generation || response.revision != revision {
                    return Err("PCR output response does not match this request.".into());
                }
                let prepared = match response.result {
                    protocol::AnalysisResult::PcrOutputs(prepared)
                        if choice.review.origin.is_none() =>
                    {
                        prepared
                    }
                    protocol::AnalysisResult::MaterialPcrOutputsPrepared(prepared)
                        if choice.review.origin.is_some() =>
                    {
                        prepared
                    }
                    protocol::AnalysisResult::Error(error) => return Err(error),
                    _ => return Err("Unexpected PCR output response.".into()),
                };
                prepared
                    .validate_binding(&revision)
                    .map_err(|e| e.to_string())?;
                if prepared.context_revision != revision
                    || prepared.selected_index != choice.product_index
                    || prepared.product.is_some() != choice.save_product
                    || prepared.primers.is_some() != choice.save_primers
                {
                    return Err("PCR output selection does not match the reviewed request.".into());
                }
                let output_findings = prepared.findings.join(" ");
                let mut staged_docs = documents.get_untracked();
                let mut staged_active = active.get_untracked();
                let mut staged_assets = attachments.get_untracked();
                let mut product_event = None;
                let mut primer_asset = None;
                if let Some(mut product) = prepared.product {
                    product.metadata.insert(
                        "pcr:template_document".into(),
                        choice.review.source.document_id.clone(),
                    );
                    product.metadata.insert(
                        "pcr:template_event".into(),
                        choice.review.source.event_id.clone(),
                    );
                    let source = staged_docs
                        .iter()
                        .find(|doc| current_event_binding(doc) == choice.review.source)
                        .ok_or("The exact PCR template event changed; recompute products.")?;
                    let prefix = format!("pcr-{}", &product.revision()[..16]);
                    let mut suffix = 1usize;
                    let id = loop {
                        let id = format!("{prefix}-{suffix}");
                        if !staged_docs.iter().any(|doc| doc.id == id) {
                            break id;
                        }
                        suffix += 1;
                    };
                    let derived = source.derive_event(&choice.review.source.event_id,id,product,
                        "Selected exact-PCR reference preview; ordered primers retained, physical chemistry and experimental success unverified")
                        .map_err(|e|e.to_string())?;
                    product_event = Some(current_event_binding(&derived));
                    staged_active = staged_docs.len();
                    staged_docs.push(derived);
                }
                if staged_docs.len() > 1000
                    || staged_docs
                        .iter()
                        .map(|d| document_bytes(&d.molecule))
                        .sum::<usize>()
                        > MAX_WORKSPACE_BYTES
                {
                    return Err("Selected PCR output exceeds the workspace document budget; nothing created.".into());
                }
                let retained_table = prepared.primers.map(std::sync::Arc::new);
                if let Some(table) = &retained_table {
                    table.validate().map_err(|e| e.to_string())?;
                    let bytes = serde_json::to_vec(table.as_ref()).map_err(|e| e.to_string())?;
                    let asset_hash = blake3::hash(&bytes).to_hex().to_string();
                    primer_asset = Some(asset_hash.clone());
                    staged_assets.insert(
                        asset_hash,
                        workspace::WorkspaceAttachment {
                            kind: "oligo-table-v1".into(),
                            name: format!("{} primers", choice.name),
                            media_type: "application/json".into(),
                            bytes,
                        },
                    );
                }
                let receipt = serde_json::to_vec(&serde_json::json!({
                    "schema_version":1,"source":choice.review.source,"context_revision":revision,
                    "selected_index":choice.product_index,"original_preview_revision":preview_revision,
                    "product":product_event,"primer_attachment":primer_asset,
                    "geometry":{"forward_start":prepared.geometry.forward_start,
                        "reverse_start":prepared.geometry.reverse_start,"crosses_origin":prepared.geometry.crosses_origin,
                        "forward_primer":prepared.geometry.forward_primer,"reverse_primer":prepared.geometry.reverse_primer},
                    "forward":choice.review.forward,"reverse":choice.review.reverse,
                    "output_name":choice.name,"save_product":choice.save_product,"save_primers":choice.save_primers,
                    "model":"Exact-annealing reference preview; no independent target verification or experimental success asserted"
                })).map_err(|e|e.to_string())?;
                staged_assets.insert(
                    blake3::hash(&receipt).to_hex().to_string(),
                    workspace::WorkspaceAttachment {
                        kind: "pcr-output-reference-v1".into(),
                        name: format!("{} PCR source", choice.name),
                        media_type: "application/json".into(),
                        bytes: receipt,
                    },
                );
                workspace_snapshot(
                    &staged_docs,
                    staged_active,
                    candidate.get_untracked().as_ref(),
                    &stored_candidates.get_untracked(),
                    &staged_assets,
                    &workspace_name.get_untracked(),
                )?
                .validate()
                .map_err(|e| e.to_string())?;
                batch(|| {
                    documents.set(staged_docs);
                    active.set(staged_active);
                    attachments.set(staged_assets);
                    if let Some(table) = retained_table {
                        oligo_tables.update(|tables| {
                            if !tables.iter().any(|old| {
                                old.source_blake3 == table.source_blake3
                                    && old.options == table.options
                            }) {
                                tables.push(table);
                            }
                        });
                    }
                    pcr_review.set(None);
                    message.set(format!("Selected PCR outputs created. Sequence-reference previews; physical chemistry and experimental success are unverified. {output_findings}"));
                });
                Ok(())
            })();
            if let Err(error) = outcome {
                message.set(error);
            }
        });
    });
    let run_plan =
        move |plan: DesignPlan,
              target: Molecule,
              exact_sources: Option<Vec<collagene_core::run_lineage::InputEventBinding>>| {
            let expected_epoch = workspace_epoch.get_untracked();
            let expected_target = documents
                .with_untracked(|docs| docs.get(active.get_untracked()).map(current_event_binding));
            let expected_sources = documents.with_untracked(|docs| {
                docs.iter()
                    .map(|d| {
                        (
                            d.id.clone(),
                            d.head_event().to_string(),
                            d.head_revision().to_string(),
                        )
                    })
                    .collect::<Vec<_>>()
            });
            last_plan_events.set(None);
            plan_report.set(None);
            last_plan.set(None);
            last_plan_target.set(None);
            let revision = target.revision();
            let generation = plan_generation.get_untracked();
            let request = protocol::Request {
                id: 0,
                revision: revision.clone(),
                job: protocol::AnalysisJob::Verify {
                    target: target.clone(),
                    plan: plan.clone(),
                },
            };
            message.set("Reconstructing the assembly in background…".into());
            spawn_local(async move {
                let result = worker::analyze(request).await;
                if molecule.get_untracked().revision() != revision
                    || plan_generation.get_untracked() != generation
                    || workspace_epoch.get_untracked() != expected_epoch
                    || documents.with_untracked(|docs| {
                        docs.get(active.get_untracked()).map(current_event_binding)
                    }) != expected_target
                    || documents.with_untracked(|docs| {
                        docs.iter()
                            .map(|d| {
                                (
                                    d.id.clone(),
                                    d.head_event().to_string(),
                                    d.head_revision().to_string(),
                                )
                            })
                            .collect::<Vec<_>>()
                    }) != expected_sources
                {
                    return;
                }
                match result {
                    Ok(protocol::Response {
                        result: protocol::AnalysisResult::Verification(report),
                        ..
                    }) => {
                        plan_report.set(Some(report));
                        last_plan_events.set(
                            exact_sources
                                .map(|bindings| (verification::plan_revision(&plan), bindings)),
                        );
                        last_plan_target.set(expected_target.map(|event| {
                            (verification::plan_revision(&plan), expected_epoch, event)
                        }));
                        last_plan.set(Some((plan, target)));
                        message.set("Assembly reconstruction finished.".into());
                    }
                    Ok(protocol::Response {
                        result: protocol::AnalysisResult::Error(error),
                        ..
                    })
                    | Err(error) => message.set(error),
                    _ => message.set("Unexpected worker response".into()),
                }
            });
        };
    let verify_candidate = move |_: web_sys::MouseEvent| {
        let Some(candidate) = candidate.get_untracked() else {
            return;
        };
        let preferences = candidate_events
            .get_untracked()
            .filter(|(revision, _, _, _)| *revision == candidate.plan_revision);
        let is_reference = candidate.plan.steps.iter().any(|step| {
            matches!(
                step,
                Operation::RecombineReference { .. }
                    | Operation::MultiSiteRecombineReference { .. }
                    | Operation::Topo { .. }
                    | Operation::CheckedTopo { .. }
                    | Operation::ProcessedTopo { .. }
                    | Operation::TopoActions { .. }
            )
        });
        // Exact native source preferences may name retained events even when a
        // different version is currently open. The native verifier reconstructs
        // their declarations; reading history grants no physical authority.
        let exact_preferences = preferences.is_some();
        if let Some((_, epoch, target, sources)) = &preferences {
            let valid = *epoch == workspace_epoch.get_untracked()
                && documents.with_untracked(|docs| {
                    docs.get(active.get_untracked())
                        .is_some_and(|d| current_event_binding(d) == *target)
                        && sources.iter().all(|s| {
                            docs.iter().any(|d| {
                                if exact_preferences {
                                    d.id == s.event.document_id
                                        && d.revision_molecule(&s.event.event_id)
                                            .is_ok_and(|m| m.revision() == s.event.revision)
                                } else {
                                    current_event_binding(d) == s.event
                                }
                            })
                        })
                });
            if !valid {
                message.set("Pinned source or target events are unavailable or changed. Inspect the sources again before verification.".into());
                return;
            }
        }
        let captured = if exact_preferences {
            documents.with_untracked(|docs| -> Result<(Molecule, Vec<Molecule>), String> {
                let target = docs
                    .get(active.get_untracked())
                    .ok_or("No active intended reference")?;
                let sources = &preferences
                    .as_ref()
                    .ok_or("Missing pinned source events")?
                    .3;
                let originals = crate::assembly_review::retained_inputs(docs, sources, &candidate)?;
                Ok((target.molecule.clone(), originals))
            })
        } else if is_reference {
            documents.with_untracked(|open| -> Result<(Molecule, Vec<Molecule>), String> {
                let target = open.get(active.get_untracked()).ok_or("No active target")?;
                if target.head_revision() != candidate.plan.target_revision {
                    return Err(
                        "The active target differs from this candidate's bound revision.".into(),
                    );
                }
                let bound = candidate
                    .plan
                    .inputs
                    .iter()
                    .map(|input| {
                        if let Some((_, _, _, events)) = &preferences {
                            let event = events
                                .iter()
                                .find(|e| e.input_id == input.id)
                                .ok_or("Missing exact source event preference")?;
                            if event.event.revision != input.revision {
                                return Err("Candidate/source-event revision mismatch".to_string());
                            }
                            open.iter()
                                .find(|d| current_event_binding(d) == event.event)
                                .ok_or_else(|| "The exact source event changed".to_string())
                        } else {
                            open.iter()
                                .find(|d| d.head_revision() == input.revision)
                                .ok_or_else(|| {
                                    format!("Input '{}' is missing or changed", input.id)
                                })
                        }
                    })
                    .collect::<Result<Vec<_>, String>>()?;
                let mut models = bound.iter().map(|d| &d.molecule).collect::<Vec<_>>();
                models.push(&target.molecule);
                models.extend(candidate.plan.inputs.iter().map(|i| &i.molecule));
                collagene_core::recombination::reference::model_preflight(&models)
                    .map_err(|e| e.to_string())?;
                for step in &candidate.plan.steps {
                    if let Operation::Topo{sources,..}|Operation::CheckedTopo{sources,..}|Operation::ProcessedTopo{sources,..}|Operation::TopoActions{sources,..}=step {
                        let mut source_bases=0usize;
                        for source in sources {
                            if let Source::Input(id)=&source.source {
                                let molecule=candidate.plan.inputs.iter().position(|input|input.id==*id).and_then(|index|bound.get(index)).map(|doc|&doc.molecule).ok_or("Missing TOPO source binding")?;
                                if molecule.sequence.alphabet()!=Alphabet::Dna{return Err("TOPO sources must be DNA.".into());}
                                source_bases=source_bases.checked_add(molecule.sequence.len()).filter(|&n|n<=100_000).ok_or("TOPO source DNA exceeds its 100,000-base profile; no reference is copied.")?;
                            }
                        }
                    }
                    if let Operation::RecombineReference {
                        sources, request, ..
                    } = step
                    {
                        let resolve = |source: &Source| match source {
                            Source::Input(id) => candidate
                                .plan
                                .inputs
                                .iter()
                                .position(|i| i.id == *id)
                                .and_then(|i| bound.get(i))
                                .map(|d| &d.molecule),
                            Source::Step(_) => None,
                        };
                        if let (Some(a), Some(b)) = (resolve(&sources[0]), resolve(&sources[1])) {
                            collagene_core::recombination::source_target_preflight(
                                [a, b],
                                request,
                                &target.molecule,
                            )
                            .map_err(|e| e.to_string())?;
                        }
                    }
                }
                Ok((
                    target.molecule.clone(),
                    bound.iter().map(|d| d.molecule.clone()).collect(),
                ))
            })
        } else {
            Ok((
                molecule.get_untracked(),
                documents.with_untracked(|docs| docs.iter().map(|d| d.molecule.clone()).collect()),
            ))
        };
        let (target, docs) = match captured {
            Ok(value) => value,
            Err(error) => {
                message.set(error);
                return;
            }
        };
        match candidate.bind(&target, &docs) {
            Ok(plan) => run_plan(plan, target, preferences.map(|(_, _, _, sources)| sources)),
            Err(e) => {
                plan_report.set(None);
                last_plan.set(None);
                last_plan_target.set(None);
                message.set(e.to_string());
            }
        }
    };
    let export_workspace = move |_: web_sys::MouseEvent| {
        let mut unique = std::collections::BTreeMap::new();
        documents.with_untracked(|docs| {
            for d in docs {
                unique
                    .entry(d.molecule.revision())
                    .or_insert_with(|| d.molecule.clone());
            }
        });
        let snapshot = transfer::PlanningWorkspace {
            target_revision: molecule.get_untracked().revision(),
            documents: unique.into_values().collect(),
        };
        match transfer::save_workspace(&snapshot) {
            Ok(bytes) if bytes.len()<=8*1024*1024=>{crate::app::download_with_completion("planning-workspace.cg2design","application/octet-stream",&bytes,move|result|{if message.try_get_untracked().is_none(){return;}if let Err(e)=result {message.set(e);}})},
            Ok(_)=>message.set("The initial MCP exchange supports snapshots up to 8 MiB. Open only the target and needed inputs, then export again.".into()),
            Err(e)=>message.set(e.to_string()),
        }
    };
    let do_plan = move |_: web_sys::MouseEvent| {
        let docs = documents.get_untracked();
        let a = first_input.get_untracked();
        let b = second_input.get_untracked();
        let Some(first) = docs.get(a) else {
            message.set("Select an open first source document.".into());
            return;
        };
        let second = if plan_method.get_untracked() == "restriction"
            && (!second_restriction.get_untracked()
                || restriction_final.get_untracked() != "ligate")
        {
            first
        } else {
            let Some(second) = docs.get(b) else {
                message.set("Select an open second source document.".into());
                return;
            };
            second
        };
        if plan_method.get_untracked() == "restriction" {
            let parse_index = |s: String| s.parse::<usize>().ok().and_then(|i| i.checked_sub(1));
            let Some(a_index) = parse_index(first_fragment.get_untracked()) else {
                message.set("Choose fragment numbers starting at 1.".into());
                return;
            };
            let ligate = restriction_final.get_untracked() == "ligate";
            let use_second = second_restriction.get_untracked() && ligate;
            let b_index = if use_second {
                let Some(i) = parse_index(second_fragment.get_untracked()) else {
                    message.set("Choose fragment numbers starting at 1.".into());
                    return;
                };
                i
            } else {
                0
            };
            let checked = restriction_policy.get_untracked() == "checked";
            let names = |s: String| {
                s.split(',')
                    .map(str::trim)
                    .filter(|s| !s.is_empty())
                    .map(str::to_owned)
                    .collect::<Vec<_>>()
            };
            let mut steps = Vec::new();
            let mut references = Vec::new();
            let mut inputs = Vec::new();
            let mut bindings = Vec::new();
            for (
                ordinal,
                id,
                doc,
                enzymes,
                fragment_index,
                treatment,
                repair,
                reverse,
                assumptions,
            ) in [
                (
                    0,
                    "first",
                    first,
                    first_enzymes.get_untracked(),
                    a_index,
                    first_treatment.get_untracked(),
                    first_repair.get_untracked(),
                    first_reverse.get_untracked(),
                    first_modifications.current(),
                ),
                (
                    1,
                    "second",
                    second,
                    second_enzymes.get_untracked(),
                    b_index,
                    second_treatment.get_untracked(),
                    second_repair.get_untracked(),
                    second_reverse.get_untracked(),
                    second_modifications.current(),
                ),
            ] {
                if ordinal == 1 && !use_second {
                    continue;
                }
                inputs.push(PlanInput {
                    id: id.into(),
                    revision: doc.molecule.revision(),
                    molecule: doc.molecule.clone(),
                });
                bindings.push(collagene_core::run_lineage::InputEventBinding {
                    input_id: id.into(),
                    event: collagene_core::run_lineage::EventBinding {
                        document_id: doc.id.clone(),
                        event_id: doc.head_event().into(),
                        revision: doc.molecule.revision(),
                    },
                });
                let source = if checked {
                    let declaration = steps.len();
                    steps.push(Operation::DeclareSourceModifications {
                        template: Source::Input(id.into()),
                        assumptions,
                    });
                    Source::Step(declaration)
                } else {
                    Source::Input(id.into())
                };
                let digest_index = steps.len();
                steps.push(if checked {
                    let selection = match (|| -> Result<_, String> {
                            let catalog=collagene_core::catalog_digest::builtin_catalog().map_err(|e|e.to_string())?;
                            let ids=names(enzymes).into_iter().map(|name|{
                                let mut found=catalog.edition().enzymes.iter().filter(|e|e.name==name);
                                let record=found.next().ok_or_else(||format!("Choose an exact catalog enzyme name; {name} is unavailable."))?;
                                if found.next().is_some(){return Err(format!("{name} is ambiguous; choose a unique canonical formulation."));}Ok(record.id.clone())
                            }).collect::<Result<Vec<_>,String>>()?;
                            collagene_core::catalog_digest::select(&catalog,&ids).map_err(|e|e.to_string())
                    })() {
                        Ok(selection) => selection,
                        Err(error) => { message.set(error); return; }
                    };
                    Operation::CheckedDigest {
                        template: source,
                        selection,
                        options:collagene_core::modified_duplex::checked_inspection_options(),
                        fragment_index,
                    }
                } else {
                    Operation::Digest {
                        template: source,
                        enzymes: names(enzymes),
                        fragment_index,
                    }
                });
                let mut product_index = digest_index;
                if repair != "none" {
                    let i = steps.len();
                    steps.push(Operation::RepairEnds {
                        template: Source::Step(product_index),
                        method: match repair.as_str() {
                            "fill" => chemistry::EndRepairMethod::FillFivePrimeOverhangs,
                            "chew" => chemistry::EndRepairMethod::RemoveThreePrimeOverhangs,
                            _ => chemistry::EndRepairMethod::Blunt,
                        },
                    });
                    product_index = i;
                }
                if treatment != "none" {
                    let i = steps.len();
                    steps.push(Operation::TreatEnds {
                        template: Source::Step(product_index),
                        treatment: if treatment == "phosphorylate" {
                            chemistry::EndTreatment::Phosphorylate
                        } else {
                            chemistry::EndTreatment::Dephosphorylate
                        },
                    });
                    product_index = i;
                }
                references.push(verification::OrientedSource {
                    source: Source::Step(product_index),
                    reverse,
                });
            }
            let target = molecule.get_untracked();
            if ligate {
                steps.push(collagene_core::ligation_outcomes::operation(
                    references,
                    target.topology == Topology::Circular,
                    allow_nicked.get_untracked(),
                ));
            }
            let plan = DesignPlan {
                name: format!(
                    "{} restriction {} to {}",
                    if checked {
                        "Checked"
                    } else {
                        "Assumption-based"
                    },
                    if ligate { "ligation" } else { "fragment" },
                    target.name
                ),
                target_revision: target.revision(),
                inputs,
                steps,
            };
            run_plan(plan, target, Some(bindings));
            return;
        }
        let Ok(overlap) = overlap.get_untracked().parse::<usize>() else {
            message.set("Overlap must be a positive integer.".into());
            return;
        };
        let target = molecule.get_untracked();
        let circular = target.topology == Topology::Circular;
        let mut overlaps = vec![overlap];
        if circular {
            let Ok(close) = closure_overlap.get_untracked().parse::<usize>() else {
                message.set("Closure overlap must be a positive integer.".into());
                return;
            };
            overlaps.push(close);
        }
        let plan = DesignPlan {
            name: format!("Assembly to {}", target.name),
            target_revision: target.revision(),
            inputs: vec![
                PlanInput {
                    id: "first".into(),
                    revision: first.molecule.revision(),
                    molecule: first.molecule.clone(),
                },
                PlanInput {
                    id: "second".into(),
                    revision: second.molecule.revision(),
                    molecule: second.molecule.clone(),
                },
            ],
            steps: vec![Operation::Overlap {
                fragments: vec![
                    Source::Input("first".into()),
                    Source::Input("second".into()),
                ],
                overlaps,
                circular,
            }],
        };
        let events = [("first", first), ("second", second)]
            .into_iter()
            .map(|(id, doc)| collagene_core::run_lineage::InputEventBinding {
                input_id: id.into(),
                event: collagene_core::run_lineage::EventBinding {
                    document_id: doc.id.clone(),
                    event_id: doc.head_event().into(),
                    revision: doc.molecule.revision(),
                },
            })
            .collect();
        run_plan(plan, target, Some(events));
    };

    let use_digest_in_assembly = Callback::new(
        move |fragment: crate::digest::AssemblyFragment| -> Result<(), String> {
            if workspace_restoring.get_untracked()
                || workspace_epoch.get_untracked() != fragment.epoch
            {
                return Err("The workspace changed; run the digest again.".into());
            }
            documents.with_untracked(|docs| -> Result<(), String> {
                let doc = docs
                    .iter()
                    .find(|d| current_event_binding(d) == fragment.source.binding)
                    .ok_or("The original digest source event changed; run the digest again.")?;
                if doc.molecule.revision() != fragment.molecule.revision() {
                    return Err("The digest source snapshot changed.".into());
                }
                Ok(())
            })?;
            let mut staged = assembly_draft.get_untracked();
            crate::assembly_workspace::set_digest_fragment(&mut staged, fragment)?;
            assembly_draft.set(staged);
            assembly_visible.set(true);
            message.set("Selected digest fragment placed in the active assembly slot. Its original source and pinned digest preparation are retained.".into());
            Ok(())
        },
    );
    let finalize_assembly = Callback::new(
        move |request: crate::assembly_workspace::ReviewRequest| -> Result<(), String> {
            if workspace_restoring.get_untracked() {
                return Err("Wait for the workspace operation to finish.".into());
            }
            let frame = documents.with_untracked(|docs| {
                workspace_snapshot(
                    docs,
                    active.get_untracked(),
                    candidate.get_untracked().as_ref(),
                    &stored_candidates.get_untracked(),
                    &attachments.get_untracked(),
                    &workspace_name.get_untracked(),
                )
            })?;
            let fingerprint = collagene_core::run_lineage::workspace_fingerprint(&frame)
                .map_err(|e| e.to_string())?;
            let staged = crate::assembly_review::stage_any(
                &frame,
                &request.intent,
                &request.sources,
                &request.candidate,
            )?;
            restore_workspace_state(&staged.workspace)?;
            pending_finalize.set(Some(std::sync::Arc::new(PendingFinalize {
                workspace: std::sync::Arc::new(staged.workspace),
                plan: request.candidate.plan,
                sources: staged.sources,
                source_fingerprint: fingerprint,
            })));
            modal.set("assembly_finalize");
            Ok(())
        },
    );
    let finalize_reviewed_candidate = Callback::new(
        move |value: transfer::Candidate| -> Result<(), String> {
            let (plan, _) = capture_junction_plan.run(())?;
            if value.plan_revision != verification::plan_revision(&plan) {
                return Err("The reviewed reconstruction changed.".into());
            }
            let (_, epoch, target) = last_plan_target
                .get_untracked()
                .ok_or("Missing exact target event")?;
            if epoch != workspace_epoch.get_untracked() {
                return Err("The workspace changed; reconstruct again.".into());
            }
            let sources=last_plan_events.get_untracked().filter(|(revision,_)|*revision==value.plan_revision)
            .map(|(_,sources)|sources).ok_or("Exact source events are unavailable. Review a source-bound design before creating outputs.")?;
            let frame = documents.with_untracked(|docs| {
                workspace_snapshot(
                    docs,
                    active.get_untracked(),
                    candidate.get_untracked().as_ref(),
                    &stored_candidates.get_untracked(),
                    &attachments.get_untracked(),
                    &workspace_name.get_untracked(),
                )
            })?;
            let fingerprint = collagene_core::run_lineage::workspace_fingerprint(&frame)
                .map_err(|e| e.to_string())?;
            let staged = crate::assembly_review::stage_existing(&frame, &value, &target, &sources)?;
            restore_workspace_state(&staged.workspace)?;
            pending_finalize.set(Some(std::sync::Arc::new(PendingFinalize {
                workspace: std::sync::Arc::new(staged.workspace),
                plan: value.plan,
                sources: staged.sources,
                source_fingerprint: fingerprint,
            })));
            modal.set("assembly_finalize");
            Ok(())
        },
    );
    let reviewed_product_current = Memo::new(move |_| {
        let _ = (workspace_epoch.get(), workspace_restoring.get());
        let Some((revision, epoch, binding)) = last_plan_target.get() else {
            return false;
        };
        epoch == workspace_epoch.get()
            && !workspace_restoring.get()
            && documents.with(|docs| {
                docs.get(active.get())
                    .is_some_and(|doc| current_event_binding(doc) == binding)
            })
            && plan_report.with(|report| {
                report.as_ref().is_some_and(|report| {
                    report.plan_revision == revision
                        && report.target_revision == binding.revision
                        && report.engine_version == verification::ENGINE_VERSION
                })
            })
            && last_plan.with(Option::is_some)
    });
    let reviewed_product_guard = Callback::new(move |()| capture_junction_plan.run(()).map(|_| ()));
    let finalize_snapshot = Callback::new(
        move |()| -> Result<collagene_core::workspace::Workspace, String> {
            let pending = pending_finalize
                .get_untracked()
                .ok_or("No assembly output choice is pending.")?;
            let current = documents.with_untracked(|docs| {
                workspace_snapshot(
                    docs,
                    active.get_untracked(),
                    candidate.get_untracked().as_ref(),
                    &stored_candidates.get_untracked(),
                    &attachments.get_untracked(),
                    &workspace_name.get_untracked(),
                )
            })?;
            if collagene_core::run_lineage::workspace_fingerprint(&current)
                .map_err(|e| e.to_string())?
                != pending.source_fingerprint
            {
                return Err(
                    "The workspace changed while choosing outputs. Review this assembly again."
                        .into(),
                );
            }
            Ok(pending.workspace.as_ref().clone())
        },
    );
    let review_assembly = Callback::new(
        move |request: crate::assembly_workspace::ReviewRequest| -> Result<(), String> {
            if workspace_restoring.get_untracked() {
                return Err("Wait for the workspace operation to finish.".into());
            }
            let frame = documents.with_untracked(|docs| {
                workspace_snapshot(
                    docs,
                    active.get_untracked(),
                    candidate.get_untracked().as_ref(),
                    &stored_candidates.get_untracked(),
                    &attachments.get_untracked(),
                    &workspace_name.get_untracked(),
                )
            })?;
            let staged = crate::assembly_review::stage_any(
                &frame,
                &request.intent,
                &request.sources,
                &request.candidate,
            )?;
            let state = restore_workspace_state(&staged.workspace)?;
            let preferred = (
                request.candidate.plan_revision.clone(),
                workspace_epoch.get_untracked(),
                staged.target,
                staged.sources,
            );
            batch(|| {
                documents.set(state.documents);
                active.set(state.active);
                stored_candidates.set(state.candidates);
                candidate_events.set(None);
                candidate.set(Some(request.candidate));
                candidate_events.set(Some(preferred));
                plan_report.set(None);
                last_plan.set(None);
                last_plan_target.set(None);
                last_plan_events.set(None);
                plan_generation.update(|g| *g = g.wrapping_add(1));
                assembly_visible.set(false);
                panel.set("plan");
                message.set("Opened the intended reference and proposed preparation. Verify the assembly to reconstruct it from the pinned source events.".into());
            });
            Ok(())
        },
    );
    let navigate_assembly = Callback::new(
        move |request: crate::assembly_workspace::SelectionNavigation| {
            if workspace_restoring.get_untracked() {
                return;
            }
            let opened =
                documents.with_untracked(|docs| -> Result<(Vec<OpenDocument>, usize), String> {
                    let index = docs
                        .iter()
                        .position(|d| d.id == request.source.document_id)
                        .ok_or("The pinned source document is unavailable.")?;
                    let source = &docs[index];
                    let exact = source
                        .revision_molecule(&request.source.event_id)
                        .map_err(|e| e.to_string())?;
                    if exact.revision() != request.source.revision {
                        return Err(
                            "The source event no longer matches the selected fragment.".into()
                        );
                    }
                    if current_event_binding(source) == request.source {
                        return Ok((docs.clone(), index));
                    }
                    if docs.len() >= 1000 {
                        return Err("The workspace document limit has been reached.".into());
                    }
                    let prefix = format!("assembly-source-{}", &request.source.revision[..16]);
                    let id = (1..=1001)
                        .map(|n| format!("{prefix}-{n}"))
                        .find(|id| !docs.iter().any(|d| d.id == *id))
                        .ok_or("No source document identity is available.")?;
                    let fork = source
                        .fork_event(&request.source.event_id, id)
                        .map_err(|e| e.to_string())?;
                    let mut staged = docs.clone();
                    let index = staged.len();
                    staged.push(fork);
                    workspace_snapshot(
                        &staged,
                        index,
                        candidate.get_untracked().as_ref(),
                        &stored_candidates.get_untracked(),
                        &attachments.get_untracked(),
                        &workspace_name.get_untracked(),
                    )?
                    .validate()
                    .map_err(|e| e.to_string())?;
                    Ok((staged, index))
                });
            let (staged, index) = match opened {
                Ok(v) => v,
                Err(error) => {
                    message.set(error);
                    return;
                }
            };
            let expected = staged[index].head_revision().to_owned();
            batch(|| {
                documents.set(staged);
                active.set(index);
            });
            let apply = wasm_bindgen::closure::Closure::once_into_js(move || {
                if documents.with_untracked(|docs| {
                    docs.get(active.get_untracked())
                        .is_none_or(|d| d.head_revision() != expected)
                }) {
                    return;
                }
                let n = molecule.with_untracked(|m| m.sequence.len());
                if n == 0
                    || request.start >= n
                    || request.end <= request.start
                    || request.end > request.start.saturating_add(n)
                {
                    return;
                }
                selected_feature.set(None);
                let parts = if request.end <= n {
                    vec![(request.start, request.end)]
                } else {
                    vec![(request.start, n), (0, request.end - n)]
                };
                selection.set(if parts.len() == 1 {
                    Some(parts[0])
                } else {
                    None
                });
                search_regions.set(Some((expected, parts)));
                message.set(if request.reverse {
                    "Selected the pinned source region; this fragment is reversed in the assembly."
                        .into()
                } else {
                    "Selected the pinned assembly source region.".into()
                });
            });
            if let Some(window) = web_sys::window() {
                let _ = window.request_animation_frame(apply.unchecked_ref());
            }
        },
    );

    let variant_asset =
        move || -> Result<(String, collagene_core::workspace::WorkspaceAttachment), String> {
            let Some((source, format, name)) = pending_variants.get_untracked() else {
                return Err("No variant source is open.".into());
            };
            collagene_core::variants_archive::make_attachment(&source, format, &name)
                .map_err(|e| e.to_string())
        };
    let keep_variant_source = Callback::new(move |()| -> Result<(), String> {
        if workspace_restoring.get_untracked() {
            return Err("Wait for workspace restoration.".into());
        }
        let (hash, asset) = variant_asset()?;
        let mut assets = attachments.get_untracked();
        if let Some(existing) = assets.get(&hash) {
            if existing.kind != asset.kind {
                return Err(
                    "These original bytes are already retained under a different source format."
                        .into(),
                );
            }
            return Ok(());
        }
        assets.insert(hash, asset);
        let frame = documents.with_untracked(|docs| {
            workspace_snapshot(
                docs,
                active.get_untracked(),
                candidate.get_untracked().as_ref(),
                &stored_candidates.get_untracked(),
                &assets,
                &workspace_name.get_untracked(),
            )
        })?;
        frame.validate().map_err(|e| e.to_string())?;
        attachments.set(assets);
        Ok(())
    });
    let open_variant_preview = Callback::new(
        move |(mut molecule, identity): (Molecule, String)| -> Result<(), String> {
            if workspace_restoring.get_untracked() {
                return Err("Wait for workspace restoration.".into());
            }
            collagene_core::model::bounds::model_preflight(&[&molecule])
                .map_err(|e| e.to_string())?;
            if molecule
                .metadata
                .get("source:legacy:chemistry_unmapped")
                .map(String::as_str)
                != Some("true")
                || molecule
                    .metadata
                    .get("variant:reference_only")
                    .map(String::as_str)
                    != Some("true")
            {
                return Err(
                    "The allele preview lacks its reference-only chemistry restriction.".into(),
                );
            }
            let (hash, asset) = variant_asset()?;
            if molecule.metadata.get("variant:unverified_source") != Some(&hash) {
                return Err("The variant source changed. Check the selected allele again.".into());
            }
            molecule.name = format!(
                "{} · allele design",
                molecule.name.chars().take(200).collect::<String>()
            );
            molecule.validate().map_err(|e| e.to_string())?;
            let mut assets = attachments.get_untracked();
            if assets.get(&hash).is_some_and(|a| a.kind != asset.kind) {
                return Err("Original bytes are retained under another source format.".into());
            }
            assets.entry(hash).or_insert(asset);
            let mut staged=documents.with_untracked(|docs|->Result<Vec<OpenDocument>,String>{
            if docs.len()>=1000||docs.iter().map(|d|document_bytes(&d.molecule)).sum::<usize>().saturating_add(document_bytes(&molecule))>MAX_WORKSPACE_BYTES{return Err("This allele design exceeds the open-document budget.".into());}
            let source=docs.iter().find(|doc|format!("{}:{}:{}",workspace_epoch.get_untracked(),doc.id,doc.head_event())==identity).ok_or("The exact reference document event changed. Check the allele again.")?;
            let mut suffix=1usize;let id=loop{let id=format!("variant-design-{}-{suffix}",&source.head_revision()[..12]);if !docs.iter().any(|d|d.id==id){break id;}suffix+=1;};
            let derived=source.derive_event(source.head_event(),id,molecule,"Selected-allele design preview; file evidence is unverified and confers no prepared chemistry").map_err(|e|e.to_string())?;
            let mut staged=docs.clone();staged.push(derived);Ok(staged)
        })?;
            let index = staged.len() - 1;
            let frame = workspace_snapshot(
                &staged,
                index,
                candidate.get_untracked().as_ref(),
                &stored_candidates.get_untracked(),
                &assets,
                &workspace_name.get_untracked(),
            )?;
            frame.validate().map_err(|e| e.to_string())?;
            batch(|| {
                documents.set(std::mem::take(&mut staged));
                attachments.set(assets);
                active.set(index);
                modal.set("");
                message.set("Allele design opened with exact reference-event ancestry and original variant bytes. It remains a reference-only design target.".into());
            });
            Ok(())
        },
    );
    let keep_read_source = Callback::new(move |()| -> Result<(), String> {
        if workspace_restoring.get_untracked() {
            return Err("Wait for workspace restoration.".into());
        }
        let Some((source, format, name)) = pending_read_alignment.get_untracked() else {
            return Err("No read-alignment source is open.".into());
        };
        let format = match format {
            collagene_core::read_alignment_session::FileFormat::Sam => {
                collagene_core::read_alignment_store::Format::Sam
            }
            collagene_core::read_alignment_session::FileFormat::Bam => {
                collagene_core::read_alignment_store::Format::Bam
            }
        };
        let (hash, asset) =
            collagene_core::read_alignment_archive::make_attachment(&source, format, &name)
                .map_err(|e| e.to_string())?;
        if attachments.with_untracked(|assets| assets.contains_key(&hash)) {
            return Ok(());
        }
        let mut assets = attachments.get_untracked();
        assets.insert(hash, asset);
        let frame = documents.with_untracked(|docs| {
            workspace_snapshot(
                docs,
                active.get_untracked(),
                candidate.get_untracked().as_ref(),
                &stored_candidates.get_untracked(),
                &assets,
                &workspace_name.get_untracked(),
            )
        })?;
        frame.validate().map_err(|e| e.to_string())?;
        attachments.set(assets);
        message.set("Original read-alignment source retained in workspace. Reopening reparses the source; reference checks are not saved as authority.".into());
        Ok(())
    });
    let keep_alignment = Callback::new(move |item: crate::alignment_view::AlignmentRevision| {
        if workspace_restoring.get_untracked() {
            return Err("Wait for workspace restoration to finish.".into());
        }
        if let Some(key) = alignments.with_untracked(|items| {
            items
                .iter()
                .find(|a| a.revision == item.revision)
                .map(|a| a.asset_hash.clone())
        }) {
            return Ok(key);
        }
        if alignments.with_untracked(Vec::len) >= 128 {
            return Err("This web workspace retains at most128 alignment collections. Export a workspace and open fewer collections.".into());
        }
        let mut assets = attachments.get_untracked();
        let key = item.asset_hash.clone();
        assets.insert(
            key.clone(),
            collagene_core::workspace::WorkspaceAttachment {
                kind: crate::alignment_view::ATTACHMENT_KIND.into(),
                name: item.prepared.document().name.clone(),
                media_type: "application/json".into(),
                bytes: item.bytes.as_ref().clone(),
            },
        );
        let frame = documents.with_untracked(|docs| {
            workspace_snapshot(
                docs,
                active.get_untracked(),
                candidate.get_untracked().as_ref(),
                &stored_candidates.get_untracked(),
                &assets,
                &workspace_name.get_untracked(),
            )
        })?;
        frame.validate().map_err(|e| e.to_string())?;
        batch(|| {
            attachments.set(assets);
            alignments.update(|items| items.push(item));
        });
        Ok(key)
    });
    let keep_library = Callback::new(
        move |library: collagene_core::feature_library::FeatureLibrary| {
            if workspace_restoring.get_untracked() {
                return Err("Wait for workspace restoration to finish.".into());
            }
            library.validate().map_err(|error| error.to_string())?;
            let revision = library.revision();
            if libraries
                .with_untracked(|items| items.iter().any(|existing| existing.revision == revision))
            {
                return Ok(());
            }
            if libraries.with_untracked(Vec::len) >= 128 {
                return Err("This web workspace retains at most 128 immutable library revisions. Export libraries before opening a smaller workspace.".into());
            }
            let bytes = collagene_core::feature_library::save_library(&library)
                .map_err(|error| error.to_string())?;
            let hash = blake3::hash(&bytes).to_hex().to_string();
            let mut assets = attachments.get_untracked();
            assets.insert(
                hash,
                collagene_core::workspace::WorkspaceAttachment {
                    kind: "feature-library-v1".into(),
                    name: format!("Feature library {}", &revision[..12]),
                    media_type: "application/json".into(),
                    bytes,
                },
            );
            let frame = documents.with_untracked(|docs| {
                workspace_snapshot(
                    docs,
                    active.get_untracked(),
                    candidate.get_untracked().as_ref(),
                    &stored_candidates.get_untracked(),
                    &assets,
                    &workspace_name.get_untracked(),
                )
            })?;
            frame.validate().map_err(|error| error.to_string())?;
            batch(|| {
                attachments.set(assets);
                libraries.update(|items| {
                    items.push(crate::feature_tools::LibraryRevision {
                        revision,
                        library: std::sync::Arc::new(library),
                    })
                });
            });
            Ok(())
        },
    );
    let reference_copy_selection = Callback::new(
        move |()| -> Result<collagene_core::reference_transfer::Selection, String> {
            use collagene_core::reference_transfer::{Part, Selection};
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
    );
    let biology_selection_location = Callback::new(move |()| {
        molecule.with_untracked(|m| {
            let source = documents
                .with_untracked(|docs| current_event_binding(&docs[active.get_untracked()]));
            crate::biology_selection::capture(
                m,
                &source,
                selected_feature.get_untracked(),
                selection.get_untracked(),
                search_regions.get_untracked().as_ref(),
                biology_directed.get_untracked().as_ref(),
            )
        })
    });
    let select_biology_location = Callback::new(
        move |(revision, location): (String, Location)| -> Result<(), String> {
            if workspace_restoring.get_untracked() || molecule_revision.get_untracked() != revision
            {
                return Err("The ORF source changed; rerun the analysis.".into());
            }
            let source = documents
                .with_untracked(|docs| current_event_binding(&docs[active.get_untracked()]));
            let directed = molecule.with_untracked(|m| {
                crate::biology_selection::DirectedSelection::new(source, m, location)
            })?;
            batch(|| {
                selected_feature.set(None);
                selection.set(directed.contiguous());
                search_regions.set(Some((revision, directed.spans.clone())));
                biology_directed.set(Some(directed));
            });
            Ok(())
        },
    );
    Effect::new(move |_| {
        let _ = workspace_epoch.get();
        biology_directed.set(None);
        pending_translation_source.set(None);
    });
    let reference_paste_guard = Callback::new(
        move |opened: crate::reference_transfer_view::PasteOpen| -> Result<(), String> {
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
    );
    let preview_reference_paste = Callback::new(
        move |opened: crate::reference_transfer_view::PasteOpen| -> Result<(), String> {
            reference_paste_guard.run(opened.clone())?;
            pending_reference_paste.set(Some(opened));
            modal_error.set(String::new());
            modal.set("reference_paste");
            Ok(())
        },
    );
    let apply_reference_paste = Callback::new(
        move |value: crate::reference_transfer_view::Apply| -> Result<(), String> {
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
            let request = collagene_core::reference_transfer::RequestPaste {
                source: source_molecule,
                destination: destination.molecule.clone(),
                selection: opened.copied.fragment.selection.clone(),
                fragment_revision: opened.copied.fragment.context_revision.clone(),
                start: opened.start,
                end: opened.end,
            };
            let expected = collagene_core::reference_transfer::paste_revision(&request)
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
            workspace_snapshot(
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
    );
    let commit_annotations = Callback::new(move |(base, edited): (String, Molecule)| {
        commit_document(documents, active, workspace_restoring, &base, edited)?;
        selected_feature.set(None);
        selection.set(None);
        Ok::<_, String>(())
    });
    let import_oligos = Callback::new(
        move |(table, rows): (oligo_formats::OligoTable, Vec<usize>)| {
            if workspace_restoring.get_untracked() {
                return Err("Wait for workspace restoration to finish.".into());
            }
            table.validate().map_err(|e| e.to_string())?;
            let report = table.selected_molecules(&rows).map_err(|e| e.to_string())?;
            let mut staged_docs = documents.get_untracked();
            let first = staged_docs.len();
            let additions = open_histories(report.molecules, &staged_docs)?;
            staged_docs.extend(additions);
            if staged_docs.len() > 1000
                || staged_docs
                    .iter()
                    .map(|d| document_bytes(&d.molecule))
                    .sum::<usize>()
                    > MAX_WORKSPACE_BYTES
            {
                return Err(
                    "Oligo selection exceeds the workspace document limit. Select fewer rows."
                        .into(),
                );
            }
            let mut staged_attachments = attachments.get_untracked();
            let bytes = serde_json::to_vec(&table).map_err(|e| e.to_string())?;
            let label = format!(
                "Oligo collection {}",
                table.source_blake3.chars().take(12).collect::<String>()
            );
            staged_attachments.insert(
                blake3::hash(&bytes).to_hex().to_string(),
                collagene_core::workspace::WorkspaceAttachment {
                    kind: "oligo-table-v1".into(),
                    name: label,
                    media_type: "application/json".into(),
                    bytes,
                },
            );
            let staged_active = if rows.is_empty() {
                active.get_untracked()
            } else {
                first
            };
            let frame = workspace_snapshot(
                &staged_docs,
                staged_active,
                candidate.get_untracked().as_ref(),
                &stored_candidates.get_untracked(),
                &staged_attachments,
                &workspace_name.get_untracked(),
            )?;
            frame.validate().map_err(|e| e.to_string())?;
            let table = std::sync::Arc::new(table);
            batch(|| {
                documents.set(staged_docs);
                active.set(staged_active);
                attachments.set(staged_attachments);
                oligo_tables.update(|tables| {
                    if let Some(old) = tables.iter_mut().find(|old| {
                        old.source_blake3 == table.source_blake3 && old.options == table.options
                    }) {
                        *old = table;
                    } else {
                        tables.push(table);
                    }
                });
                modal.set("");
                message.set(format!(
                    "Oligo collection retained; {} selected primers opened. {}",
                    rows.len(),
                    report.warnings.join(" ")
                ));
            });
            Ok(())
        },
    );
    let open_mismatch_outputs = Callback::new(
        move |opened: crate::mutagenesis::OutputOpen| -> Result<(), String> {
            if workspace_restoring.get_untracked()
                || opened.epoch != workspace_epoch.get_untracked()
            {
                return Err(
                    "The workspace changed; compute and review the current PCR draft again.".into(),
                );
            }
            let mut staged_docs = documents.get_untracked();
            let mut staged_active = active.get_untracked();
            if staged_docs
                .get(staged_active)
                .is_none_or(|doc| current_event_binding(doc) != opened.target)
                || !staged_docs
                    .iter()
                    .any(|doc| current_event_binding(doc) == opened.source)
            {
                return Err(
                    "The original PCR source or viewed reference changed; compute again.".into(),
                );
            }
            opened.options.validate().map_err(|e| e.to_string())?;
            opened
                .prepared
                .validate_binding(&opened.revision)
                .map_err(|e| e.to_string())?;
            let prepared = opened.prepared;
            if prepared.search.source_revision != opened.source.revision
                || prepared.product.is_some() != opened.options.save_product
                || prepared.primers.is_some() != opened.options.save_primers
                || prepared
                    .product
                    .as_ref()
                    .is_some_and(|product| product.name != opened.options.name)
            {
                return Err("PCR output choices differ from the reviewed request.".into());
            }
            let mut assets = attachments.get_untracked();
            let mut product_event = None;
            if let Some(mut product) = prepared.product {
                product.metadata.insert(
                    "pcr:template_document".into(),
                    opened.source.document_id.clone(),
                );
                product
                    .metadata
                    .insert("pcr:template_event".into(), opened.source.event_id.clone());
                let source = staged_docs
                    .iter()
                    .find(|doc| current_event_binding(doc) == opened.source)
                    .ok_or("Original PCR source event is unavailable")?;
                let prefix = format!("mismatch-pcr-{}", &product.revision()[..16]);
                let mut suffix = 1usize;
                let id = loop {
                    let id = format!("{prefix}-{suffix}");
                    if !staged_docs.iter().any(|doc| doc.id == id) {
                        break id;
                    }
                    suffix += 1;
                };
                let derived=source.derive_event(&opened.source.event_id,id,product,
                "Selected intentional-mismatch PCR reference preview; no independent target or physical chemistry verification asserted")
                .map_err(|e|e.to_string())?;
                product_event = Some(current_event_binding(&derived));
                staged_active = staged_docs.len();
                staged_docs.push(derived);
            }
            if staged_docs.len() > 1000
                || staged_docs
                    .iter()
                    .map(|doc| document_bytes(&doc.molecule))
                    .sum::<usize>()
                    > MAX_WORKSPACE_BYTES
            {
                return Err("The selected PCR output exceeds the workspace document budget; nothing created.".into());
            }
            let table = prepared.primers.map(std::sync::Arc::new);
            let mut primer_asset = None;
            if let Some(table) = &table {
                table.validate().map_err(|e| e.to_string())?;
                let bytes = serde_json::to_vec(table.as_ref()).map_err(|e| e.to_string())?;
                let hash = blake3::hash(&bytes).to_hex().to_string();
                primer_asset = Some(hash.clone());
                assets.insert(
                    hash,
                    workspace::WorkspaceAttachment {
                        kind: "oligo-table-v1".into(),
                        name: format!("{} primers", opened.options.name),
                        media_type: "application/json".into(),
                        bytes,
                    },
                );
            }
            let geometry = prepared.geometry;
            let bytes=serde_json::to_vec(&serde_json::json!({
            "schema_version":1,"source":opened.source,"viewed_reference":opened.target,
            "context_revision":opened.revision,"selected_index":prepared.selected_index,
            "outputs":opened.options,"product":product_event,"primer_attachment":primer_asset,
            "search":prepared.search,"findings":prepared.findings,
            "geometry":{"forward_start":geometry.forward_start,"reverse_start":geometry.reverse_start,
                "template_length":geometry.template_length,"crosses_origin":geometry.crosses_origin,
                "forward_primer_index":geometry.forward_primer_index,"reverse_primer_index":geometry.reverse_primer_index,
                "forward_tail_length":geometry.forward_tail_length,"reverse_tail_length":geometry.reverse_tail_length,
                "forward_binding":geometry.forward_binding,"reverse_binding":geometry.reverse_binding,
                "substitutions":geometry.substitutions,"overlapping_domains":geometry.overlapping_domains,"certain":geometry.certain},
            "model":"Intentional-mismatch sequence-reference preview; search limits retained. No independent target, physical chemistry or experimental success verification asserted."
        })).map_err(|e|e.to_string())?;
            assets.insert(
                blake3::hash(&bytes).to_hex().to_string(),
                workspace::WorkspaceAttachment {
                    kind: "mismatch-pcr-output-v1".into(),
                    name: format!("{} PCR source", opened.options.name),
                    media_type: "application/json".into(),
                    bytes,
                },
            );
            workspace_snapshot(
                &staged_docs,
                staged_active,
                candidate.get_untracked().as_ref(),
                &stored_candidates.get_untracked(),
                &assets,
                &workspace_name.get_untracked(),
            )?
            .validate()
            .map_err(|e| e.to_string())?;
            batch(|| {
                documents.set(staged_docs);
                active.set(staged_active);
                attachments.set(assets);
                if let Some(table) = table {
                    oligo_tables.update(|tables| {
                        if !tables.iter().any(|old| {
                            old.source_blake3 == table.source_blake3 && old.options == table.options
                        }) {
                            tables.push(table);
                        }
                    });
                }
                message.set("Selected mismatch-PCR reference outputs retained with the original source event. Search limitations and conditional substitutions remain evidence, not a successful cloning claim.".into());
            });
            Ok(())
        },
    );
    let restore_revision = Callback::new(move |event: String| {
        if workspace_restoring.get_untracked() {
            return Err("Wait for workspace restoration to finish.".into());
        }
        let index = active.get_untracked();
        let mut docs = documents.get_untracked();
        docs.get_mut(index)
            .ok_or("No active document")?
            .checkout(&event)
            .map_err(|e| e.to_string())?;
        if docs
            .iter()
            .map(|d| document_bytes(&d.molecule))
            .sum::<usize>()
            > MAX_WORKSPACE_BYTES
        {
            return Err("Restoring this revision exceeds the current document size limit; open it in a smaller workspace.".into());
        }
        documents.set(docs);
        Ok(())
    });
    let fork_revision = Callback::new(move |event: String| {
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
        let frame = workspace_snapshot(
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
    });

    let contig_snapshot = Callback::new(move |()| {
        documents.with_untracked(|docs| {
            stored_candidates.with_untracked(|plans| {
                attachments.with_untracked(|assets| {
                    candidate.with_untracked(|current| {
                        workspace_snapshot(
                            docs,
                            active.get_untracked(),
                            current.as_ref(),
                            plans,
                            assets,
                            &workspace_name.get_untracked(),
                        )
                    })
                })
            })
        })
    });
    let navigate_pcr_overlay = Callback::new(
        move |(source, epoch, spans): (
            collagene_core::run_lineage::EventBinding,
            u64,
            Vec<(usize, usize)>,
        )|
              -> Result<(), String> {
            if workspace_restoring.get_untracked() || epoch != workspace_epoch.get_untracked() {
                return Err("The PCR workspace changed; select the product again.".into());
            }
            let publication = published_pcr_overlay
                .get_untracked()
                .ok_or("Select a current PCR geometry first.")?;
            let overlay = publication
                .overlay
                .ok_or("The PCR overlay is unavailable.")?;
            if !overlay.is_current(&source, epoch)
                || spans.is_empty()
                || spans
                    .iter()
                    .any(|span| !overlay.spans.iter().any(|s| (s.start, s.end) == *span))
            {
                return Err("The selected PCR geometry changed; select it again.".into());
            }
            let (index, molecule) = documents
                .with_untracked(|docs| {
                    docs.iter()
                        .enumerate()
                        .find(|(_, d)| current_event_binding(d) == source)
                        .map(|(i, d)| (i, d.molecule.clone()))
                })
                .ok_or("The exact PCR template event changed; compute products again.")?;
            if index == active.get_untracked() {
                batch(|| {
                    selected_feature.set(None);
                    selection.set(if spans.len() == 1 {
                        spans.first().copied()
                    } else {
                        None
                    });
                    search_regions.set(Some((source.revision, spans)));
                });
            } else {
                pending_pcr_source.set(Some(crate::pcr_source_view::Snapshot {
                    owner: publication.owner,
                    overlay,
                    molecule: std::sync::Arc::new(molecule),
                    spans,
                }));
                modal.set("pcr_source");
            }
            Ok(())
        },
    );
    let open_translated_protein = Callback::new(
        move |output: crate::biology::ProteinOutput| -> Result<(), String> {
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
                return Err("The exact translation source event changed; translate again.".into());
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
                    "The linked protein exceeds this workspace's document or size limit.".into(),
                );
            }
            let frame = workspace_snapshot(
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
    );
    let retain_translation_cds = Callback::new(
        move |output: crate::biology::AnnotationOutput| -> Result<(), String> {
            if workspace_restoring.get_untracked()
                || workspace_epoch.get_untracked() != output.epoch
            {
                return Err("The workspace changed; translate again.".into());
            }
            let current = documents
                .with_untracked(|docs| docs.get(active.get_untracked()).map(current_event_binding));
            if current.as_ref() != Some(&output.source)
                || molecule.with_untracked(|m| m != &output.request.source)
            {
                return Err("The exact translation source event changed; translate again.".into());
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
            commit_document(
                documents,
                active,
                workspace_restoring,
                &output.source.revision,
                edited,
            )?;
            message.set("Retained a CDS proposal requiring review. Its notes retain the chosen frame policy; biological CDS initiation can differ.".into());
            Ok(())
        },
    );
    let open_translation_origin = Callback::new(
        move |source: collagene_core::run_lineage::EventBinding| -> Result<(), String> {
            if workspace_restoring.get_untracked() {
                return Err("Wait for the workspace update to finish.".into());
            }
            let frame = contig_snapshot.run(())?;
            let protein = documents
                .with_untracked(|docs| current_event_binding(&docs[active.get_untracked()]));
            let original = crate::biology_selection::origin(&frame, &protein, &source)?;
            let spans = molecule
                .with_untracked(|m| m.metadata.get("translation:source_location").cloned())
                .filter(|text| text.len() <= 512 * 1024)
                .map(|text| {
                    serde_json::from_str::<Location>(&text)
                        .map_err(|e| e.to_string())
                        .and_then(|location| {
                            collagene_core::restriction_selection::footprint(
                                &location,
                                original.sequence.len(),
                            )
                            .map_err(|e| e.to_string())
                        })
                })
                .transpose()?
                .unwrap_or_default();
            pending_translation_source.set(Some(crate::pairwise_source_view::Snapshot {
                source,
                epoch: workspace_epoch.get_untracked(),
                molecule: std::sync::Arc::new(original.clone()),
                spans,
            }));
            modal.set("translation_source");
            Ok(())
        },
    );
    let set_annotation_visibility = Callback::new(
        move |(indices, hide): (Vec<usize>, bool)| -> Result<(), String> {
            if workspace_restoring.get_untracked() {
                return Err("Wait for the workspace update to finish.".into());
            }
            let mut frame = contig_snapshot.run(())?;
            let (source, count) = documents.with_untracked(|docs| {
                let doc = &docs[active.get_untracked()];
                (current_event_binding(doc), doc.molecule.features.len())
            });
            let mut hidden =
                crate::annotation_visibility::load(&frame.attachments, &source, count)?;
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
            frame.attachments = crate::annotation_visibility::store(
                &frame.attachments,
                crate::annotation_visibility::Record {
                    schema_version: 1,
                    source,
                    feature_count: count,
                    hidden,
                },
            )?;
            frame.validate().map_err(|e| e.to_string())?;
            crate::annotation_visibility::validate_workspace(&frame)?;
            attachments.set(frame.attachments);
            Ok(())
        },
    );
    let apply_annotation_batch = Callback::new(
        move |request: collagene_core::annotation_batch::Request| -> Result<(), String> {
            if workspace_restoring.get_untracked() {
                return Err("Wait for the workspace update to finish.".into());
            }
            let index = active.get_untracked();
            let mut staged_docs = documents.get_untracked();
            let doc = staged_docs
                .get_mut(index)
                .ok_or("Active annotation document is unavailable.")?;
            let previous = current_event_binding(doc);
            let edited = collagene_core::annotation_batch::apply(&doc.molecule, &request)
                .map_err(|e| e.to_string())?;
            if edited == doc.molecule {
                return Ok(());
            }
            let hidden = attachments.with_untracked(|assets| {
                crate::annotation_visibility::load(assets, &previous, doc.molecule.features.len())
            })?;
            doc.commit(edited, "Batch edited annotations")
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
            let mut frame = workspace_snapshot(
                &staged_docs,
                index,
                candidate.get_untracked().as_ref(),
                &stored_candidates.get_untracked(),
                &attachments.get_untracked(),
                &workspace_name.get_untracked(),
            )?;
            frame.attachments = crate::annotation_visibility::store(
                &frame.attachments,
                crate::annotation_visibility::Record {
                    schema_version: 1,
                    source,
                    feature_count: count,
                    hidden,
                },
            )?;
            frame.validate().map_err(|e| e.to_string())?;
            crate::annotation_visibility::validate_workspace(&frame)?;
            batch(|| {
                documents.set(staged_docs);
                attachments.set(frame.attachments);
            });
            message.set("Annotation changes saved as one undoable edit.".into());
            Ok(())
        },
    );
    let saved_enzyme_sets = Memo::new(move |_| {
        let frame = workspace_snapshot(
            &documents.get(),
            active.get_untracked(),
            candidate.get_untracked().as_ref(),
            &stored_candidates.get_untracked(),
            &attachments.get(),
            &workspace_name.get_untracked(),
        );
        match frame
            .and_then(|frame| collagene_core::enzyme_sets::list(&frame).map_err(|e| e.to_string()))
        {
            Ok(sets) => sets,
            Err(error) => {
                message.set(format!("Saved enzyme sets need attention: {error}"));
                Vec::new()
            }
        }
    });
    let save_enzyme_set = Callback::new(
        move |set: collagene_core::enzyme_sets::Set| -> Result<(), String> {
            if workspace_restoring.get_untracked() {
                return Err("Wait for the workspace update to finish.".into());
            }
            let (staged, _) = collagene_core::enzyme_sets::store(&contig_snapshot.run(())?, &set)
                .map_err(|e| e.to_string())?;
            staged.validate().map_err(|e| e.to_string())?;
            attachments.set(staged.attachments);
            message.set(format!("Saved enzyme set: {}", set.name));
            Ok(())
        },
    );
    let delete_enzyme_set = Callback::new(move |hash: String| -> Result<(), String> {
        if workspace_restoring.get_untracked() {
            return Err("Wait for the workspace update to finish.".into());
        }
        let staged = collagene_core::enzyme_sets::delete(&contig_snapshot.run(())?, &hash)
            .map_err(|e| e.to_string())?;
        staged.validate().map_err(|e| e.to_string())?;
        attachments.set(staged.attachments);
        message.set("Deleted saved enzyme set.".into());
        Ok(())
    });
    let resolve_retained_source = Callback::new(
        move |event: collagene_core::run_lineage::EventBinding| -> Result<Molecule, String> {
            documents.with_untracked(|docs| {
                let doc = docs
                    .iter()
                    .find(|doc| doc.id == event.document_id)
                    .ok_or("The retained source document is unavailable.")?;
                let molecule = doc
                    .revision_molecule(&event.event_id)
                    .map_err(|e| e.to_string())?;
                if molecule.revision() != event.revision {
                    return Err(
                        "The exact retained source event differs from its saved revision.".into(),
                    );
                }
                Ok(molecule)
            })
        },
    );
    let saved_digests = Memo::new(move |_| {
        let assets = attachments.get();
        if !assets
            .values()
            .any(|asset| asset.kind == collagene_core::digest_setups::KIND)
        {
            return Vec::new();
        }
        let frame = workspace_snapshot(
            &documents.get(),
            active.get_untracked(),
            candidate.get_untracked().as_ref(),
            &stored_candidates.get_untracked(),
            &assets,
            &workspace_name.get_untracked(),
        );
        match frame.and_then(|frame| {
            collagene_core::digest_setups::list(&frame).map_err(|e| e.to_string())
        }) {
            Ok(mut setups) => {
                setups.sort_by(|a, b| {
                    a.name
                        .cmp(&b.name)
                        .then(a.attachment_hash.cmp(&b.attachment_hash))
                });
                setups
            }
            Err(error) => {
                message.set(format!("Saved digest setups need attention: {error}"));
                Vec::new()
            }
        }
    });
    let save_digest_setup = Callback::new(
        move |(name, report, source): (
            String,
            collagene_core::catalog_digest::CatalogRestrictionInspection,
            crate::digest::SourceInfo,
        )|
              -> Result<(), String> {
            if workspace_restoring.get_untracked() {
                return Err("Wait for the workspace update to finish.".into());
            }
            if documents.with_untracked(|docs| current_event_binding(&docs[active.get_untracked()]))
                != source.binding
            {
                return Err(
                    "The active digest source event changed; inspect again before saving.".into(),
                );
            }
            let catalog =
                collagene_core::catalog_digest::builtin_catalog().map_err(|e| e.to_string())?;
            let (staged, _) = collagene_core::digest_setups::capture(
                &catalog,
                &contig_snapshot.run(())?,
                &source.binding,
                &name,
                &report,
            )
            .map_err(|e| e.to_string())?;
            staged.validate().map_err(|e| e.to_string())?;
            attachments.set(staged.attachments);
            message.set(
                "Named digest setup retained with its exact source event and historical evidence."
                    .into(),
            );
            Ok(())
        },
    );
    let load_digest_setup = Callback::new(
        move |hash: String| -> Result<collagene_core::digest_setups::Record, String> {
            if workspace_restoring.get_untracked() {
                return Err("Wait for the workspace update to finish.".into());
            }
            collagene_core::digest_setups::load(&contig_snapshot.run(())?, &hash)
                .map_err(|e| e.to_string())
        },
    );
    let delete_digest_setup = Callback::new(move |hash: String| -> Result<(), String> {
        if workspace_restoring.get_untracked() {
            return Err("Wait for the workspace update to finish.".into());
        }
        let staged = collagene_core::digest_setups::delete(&contig_snapshot.run(())?, &hash)
            .map_err(|e| e.to_string())?;
        attachments.set(staged.attachments);
        Ok(())
    });
    let queue_digest_gel = Callback::new(move |(name,recipe):(String,collagene_core::digest_setups::Recipe)| ->Result<(),String> {
        if workspace_restoring.get_untracked(){return Err("Wait for the workspace update to finish.".into());}
        recipe.validate().map_err(|e|e.to_string())?;
        resolve_retained_source.run(recipe.source.clone())?;
        if gel_digest_seeds.with_untracked(Vec::len)>=32{return Err("The pending gel queue reached its 32-lane limit.".into());}
        gel_digest_seeds.update(|seeds|seeds.push(crate::gel::DigestLaneSeed{name,recipe,epoch:workspace_epoch.get_untracked()}));
        panel.set("gel");Ok(())
    });
    let fork_digest_source = Callback::new(move |source:collagene_core::run_lineage::EventBinding| ->Result<crate::digest::SourceInfo,String> {
        if workspace_restoring.get_untracked(){return Err("Wait for the workspace update to finish.".into());}
        let mut docs=documents.get_untracked();
        let original=docs.iter().find(|doc|doc.id==source.document_id).ok_or("The retained digest source is unavailable.")?;
        let molecule=original.revision_molecule(&source.event_id).map_err(|e|e.to_string())?;
        if molecule.revision()!=source.revision{return Err("The saved digest source revision differs from its retained event.".into());}
        let prefix=format!("digest-source-{}",&source.revision[..16]);let mut suffix=1usize;
        let id=loop{let id=format!("{prefix}-{suffix}");if !docs.iter().any(|doc|doc.id==id){break id;}suffix+=1;};
        let fork=original.fork_event(&source.event_id,id).map_err(|e|e.to_string())?;
        let info=recombination_source_info(&fork);let index=docs.len();docs.push(fork);
        if docs.len()>1000||docs.iter().map(|doc|document_bytes(&doc.molecule)).sum::<usize>()>MAX_WORKSPACE_BYTES {
            return Err("The source fork exceeds this workspace's document or size limit.".into());
        }
        let frame=workspace_snapshot(&docs,index,candidate.get_untracked().as_ref(),&stored_candidates.get_untracked(),
            &attachments.get_untracked(),&workspace_name.get_untracked())?;
        frame.validate().map_err(|e|e.to_string())?;
        batch(||{documents.set(docs);active.set(index);selection.set(None);selected_feature.set(None);search_regions.set(None);});
        message.set("Opened an editable copy of the exact digest source; original source heads and saved evidence are unchanged.".into());
        Ok(info)
    });
    let saved_comparisons = Memo::new(move |_| {
        let frame = workspace_snapshot(
            &documents.get(),
            active.get_untracked(),
            candidate.get_untracked().as_ref(),
            &stored_candidates.get_untracked(),
            &attachments.get(),
            &workspace_name.get_untracked(),
        );
        match frame.and_then(|frame| {
            collagene_core::pairwise_history::list(&frame).map_err(|e| e.to_string())
        }) {
            Ok(summaries) => summaries,
            Err(error) => {
                message.set(format!("Saved comparisons need attention: {error}"));
                Vec::new()
            }
        }
    });
    let save_comparison = Callback::new(
        move |record: collagene_core::pairwise_history::Record| -> Result<(), String> {
            if workspace_restoring.get_untracked() {
                return Err("Wait for the workspace update to finish.".into());
            }
            let frame = contig_snapshot.run(())?;
            let (staged, _) = collagene_core::pairwise_history::capture(
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
    );
    let reopen_comparison = Callback::new(
        move |hash: String| -> Result<collagene_core::pairwise_history::Record, String> {
            if workspace_restoring.get_untracked() {
                return Err("Wait for the workspace update to finish.".into());
            }
            collagene_core::pairwise_history::load(&contig_snapshot.run(())?, &hash)
                .map_err(|e| e.to_string())
        },
    );
    let navigate_comparison = Callback::new(
        move |request: crate::compare::Navigation| -> Result<(), String> {
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
                            "The comparison position differs from its exact source version.".into(),
                        );
                    }
                    Ok((index, current_event_binding(doc) == request.event, molecule))
                })?;
            if !current || request.spans.len() > 1 {
                pending_comparison_source.set(Some(crate::pairwise_source_view::Snapshot {
                    source: request.event,
                    epoch: request.epoch,
                    molecule: std::sync::Arc::new(molecule),
                    spans: request.spans,
                }));
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
                let callback = wasm_bindgen::closure::Closure::once_into_js(move || apply(true));
                web_sys::window()
                    .ok_or("No browser window")?
                    .request_animation_frame(callback.unchecked_ref())
                    .map_err(|_| "Comparison navigation could not be scheduled.")?;
            }
            Ok(())
        },
    );
    let open_run_event = Callback::new(
        move |request: crate::run_event_model::Request| -> Result<(), String> {
            if workspace_restoring.get_untracked()
                || request.epoch != workspace_epoch.get_untracked()
            {
                return Err("The workspace changed; reopen the cloning run.".into());
            }
            let frame = contig_snapshot.run(())?;
            let resolved = crate::run_event_model::resolve(&frame, &request)?;
            pending_run_event.set(Some(std::sync::Arc::new(resolved)));
            modal_error.set(String::new());
            modal.set("run_event");
            Ok(())
        },
    );
    let fork_run_event = Callback::new(
        move |request: crate::run_event_model::Request| -> Result<(), String> {
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
            let staged = crate::run_event_model::fork(&frame, &request, &id)?;
            let state = restore_workspace_state(&staged)?;
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
    );
    let back_to_run = Callback::new(move |run_id: String| {
        if workspace_restoring.get_untracked() {
            return;
        }
        selected_run.set(run_id);
        pending_run_event.set(None);
        modal_error.set(String::new());
        modal.set("run_details");
    });
    let record_read_check = Callback::new(
        move |request: crate::trace_validation::RecordRequest| -> Result<(), String> {
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
    );
    let read_check_snapshot = Callback::new(
        move |()| -> Result<collagene_core::workspace::Workspace, String> {
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
                return Err("The reference or read events changed. Recheck before saving.".into());
            }
            contig_snapshot.run(())
        },
    );
    let restore_read_check = Callback::new(move |frame: collagene_core::workspace::Workspace| {
        let pinned = read_reference.get_untracked();
        match restore_workspace_state(&frame) {
            Ok(state) => batch(|| {
                apply_restored.run(state);
                read_reference.set(pinned);
                panel.set("reads");
            }),
            Err(error) => message.set(error),
        }
    });
    let navigate_check_reference = Callback::new(
        move |request: crate::trace_validation::ReferenceNavigation| -> Result<(), String> {
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
            workspace_snapshot(
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
    );
    let open_read_review = Callback::new(move |hash: String| -> Result<(), String> {
        if workspace_restoring.get_untracked() {
            return Err("Wait for the workspace operation to finish.".into());
        }
        let frame = contig_snapshot.run(())?;
        let saved = collagene_core::read_checks::load(&frame, &hash).map_err(|e| e.to_string())?;
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
        let saved = collagene_core::read_review::load(&frame, &hash).map_err(|e| e.to_string())?;
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
        move |request: collagene_core::read_review::Request| -> Result<(), String> {
            if workspace_restoring.get_untracked() {
                return Err("Wait for the workspace operation to finish.".into());
            }
            if reviewed_check.get_untracked().as_deref() != Some(request.attachment_hash.as_str()) {
                return Err("The selected read check changed; preview the decisions again.".into());
            }
            let frame = contig_snapshot.run(())?;
            let saved = collagene_core::read_checks::load(&frame, &request.attachment_hash)
                .map_err(|e| e.to_string())?;
            let owner = read_check_reference.get_untracked();
            if saved.reference.document_id != owner.document_id {
                return Err("The checked construct changed; preview the decisions again.".into());
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
        move |()| -> Result<collagene_core::workspace::Workspace, String> {
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
    let restore_read_review = Callback::new(move |commit: collagene_core::read_review::Commit| {
        let pinned = commit.target.or_else(|| read_reference.get_untracked());
        match restore_workspace_state(&commit.workspace) {
            Ok(state) => batch(|| {
                apply_restored.run(state);
                read_reference.set(pinned);
                panel.set("reads");
            }),
            Err(error) => message.set(error),
        }
    });
    let open_mutation_intent = Callback::new(
        move |intent: crate::mutagenesis::IntentOpen| -> Result<(), String> {
            if workspace_restoring.get_untracked() {
                return Err("Wait for the workspace operation to finish.".into());
            }
            let frame = contig_snapshot.run(())?;
            if !documents.with_untracked(|docs| {
                docs.iter()
                    .any(|doc| current_event_binding(doc) == intent.source)
            }) {
                return Err(
                    "The mutation source event changed; preview the intended edit again.".into(),
                );
            }
            let staged = crate::mutation_review::stage(&frame, &intent.source, &intent.report)?;
            restore_workspace_state(&staged)?;
            let restored = OpenDocument::restore_all(&staged).map_err(|e| e.to_string())?;
            let index = restored
                .iter()
                .position(|doc| doc.id == staged.active_document)
                .ok_or("The intended reference is unavailable.")?;
            batch(|| {
                documents.set(restored);
                active.set(index);
                attachments.set(staged.attachments);
                selection.set(None);
                selected_feature.set(None);
                plan_report.set(None);
                last_plan.set(None);
                last_plan_target.set(None);
                message.set("Opened the declared intended reference with its exact source history. Design and verify a preparation against this target next.".into());
            });
            Ok(())
        },
    );
    let open_domestication = Callback::new(
        move |opened: crate::domestication_review::Open| -> Result<(), String> {
            if workspace_restoring.get_untracked()
                || opened.epoch != workspace_epoch.get_untracked()
            {
                return Err(
                    "The domestication workspace changed; inspect the current source again.".into(),
                );
            }
            if !documents.with_untracked(|docs| {
                docs.get(active.get_untracked())
                    .is_some_and(|doc| current_event_binding(doc) == opened.source)
            }) {
                return Err("The exact original source event changed.".into());
            }
            let staged = crate::domestication_review::stage(&contig_snapshot.run(())?, &opened)?;
            restore_workspace_state(&staged.workspace)?;
            let restored =
                OpenDocument::restore_all(&staged.workspace).map_err(|e| e.to_string())?;
            let index = restored
                .iter()
                .position(|doc| doc.id == staged.target.document_id)
                .ok_or("The independently checked target was not retained.")?;
            let preference = staged.candidate.as_ref().map(|value| {
                (
                    value.plan_revision.clone(),
                    opened.epoch,
                    staged.target,
                    staged.inputs,
                )
            });
            batch(|| {
                documents.set(restored);
                active.set(index);
                attachments.set(staged.workspace.attachments);
                stored_candidates.set(staged.workspace.candidates);
                candidate_events.set(None);
                candidate.set(staged.candidate);
                candidate_events.set(preference);
                selection.set(None);
                selected_feature.set(None);
                search_regions.set(None);
                plan_report.set(None);
                last_plan.set(None);
                last_plan_target.set(None);
                last_plan_events.set(None);
                plan_generation.update(|g| *g = g.wrapping_add(1));
                if opened.prepared.is_some() {
                    panel.set("plan");
                    message.set("Opened the original-source mutagenic preparation and independently checked intended reference. Verify the plan, inspect actual products, then choose outputs.".into());
                } else {
                    message.set("Opened checked substitution intent with unchanged annotation coordinates and exact source ancestry. This reference is not prepared DNA.".into());
                }
            });
            Ok(())
        },
    );
    let contig_restore = Callback::new(move |frame: collagene_core::workspace::Workspace| {
        match restore_workspace_state(&frame) {
            Ok(state) => batch(|| apply_restored.run(state)),
            Err(error) => message.set(error),
        }
    });
    let record_assembly = Callback::new(
        move |request: crate::contigs::ArchiveRequest| -> Result<(), String> {
            if workspace_restoring.get_untracked() {
                return Err("Wait for the workspace update to finish.".into());
            }
            if workspace_epoch.get_untracked() != request.epoch {
                return Err("The workspace changed; reopen the assembly.".into());
            }
            pending_assembly.set(Some(request));
            modal_error.set(String::new());
            modal.set("contig_archive");
            Ok(())
        },
    );
    let record_donor = Callback::new(
        move |request: crate::crispr::DonorRequest| -> Result<(), String> {
            if workspace_restoring.get_untracked() {
                return Err("Wait for the native workspace update to finish.".into());
            }
            if request.epoch != workspace_epoch.get_untracked()
                || !documents.with_untracked(|docs| {
                    [&request.source, &request.target].into_iter().all(|info| {
                        docs.iter().any(|d| {
                            d.id == info.binding.document_id
                                && d.head_event() == info.binding.event_id
                                && d.molecule.revision() == info.binding.revision
                        })
                    })
                })
            {
                return Err(
                    "Source or independent target event changed. Design and check again.".into(),
                );
            }
            pending_donor.set(Some(request));
            modal_error.set(String::new());
            modal.set("crispr_product");
            Ok(())
        },
    );
    let open_digest = Callback::new(
        move |request: crate::digest_product::DigestOpenRequest| -> Result<(), String> {
            if workspace_restoring.get_untracked() {
                return Err("Wait for the workspace update to finish.".into());
            }
            if request.epoch != workspace_epoch.get_untracked()
                || !documents.with_untracked(|docs| {
                    docs.get(active.get_untracked()).is_some_and(|d| {
                        d.id == request.source.document_id
                            && d.head_event() == request.source.event_id
                            && molecule_revision.get_untracked() == request.source.revision
                    })
                })
            {
                return Err("The source document/event changed. Run the digest again before opening a fragment.".into());
            }
            pending_digest.set(Some(request));
            modal_error.set(String::new());
            modal.set("digest_product");
            Ok(())
        },
    );
    let navigate_digest = Callback::new(move |request: crate::digest::Navigation| {
        if workspace_restoring.get_untracked()
            || request.epoch != workspace_epoch.get_untracked()
            || !documents.with_untracked(|docs| {
                docs.get(active.get_untracked()).is_some_and(|d| {
                    d.id == request.source.document_id
                        && d.head_event() == request.source.event_id
                        && molecule_revision.get_untracked() == request.source.revision
                })
            })
        {
            return;
        }
        let n = molecule.with_untracked(|m| m.sequence.len());
        if request.spans.is_empty() || request.spans.iter().any(|(a, b)| a >= b || *b > n) {
            return;
        }
        batch(|| {
            selected_feature.set(None);
            selection.set(if request.spans.len() == 1 {
                request.spans.first().copied()
            } else {
                None
            });
            search_regions.set(Some((request.source.revision, request.spans)));
        });
    });
    let navigate_read_call = Callback::new(move |request: crate::contigs::ReadNavigation| {
        if workspace_restoring.get_untracked() || workspace_epoch.get_untracked() != request.epoch {
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
            let frame = workspace_snapshot(
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
    });

    view! {
                        <Show when=move||workspace_restoring.get()&&!matches!(modal.get(),"assembly_finalize"|"read_check_capture"|"read_review_commit")><div role="status" aria-live="polite" style="position:fixed;inset:0;z-index:10000;background:#f5f8f4ed;display:grid;place-items:center;pointer-events:all"><p>"Updating workspace…"</p></div></Show>
                        <div class="application workspace-shell" on:click=dismiss_viewer_menus><div class="workspace-contents" inert=move||workspace_restoring.get()>
                            <header class="topbar">
                                <a class="brand" href="#" aria-label="Collagene workspace"><svg viewBox="0 0 32 32" aria-hidden="true"><circle cx="16" cy="16" r="11" fill="none" stroke="currentColor" stroke-width="3"/><path d="M9 7L23 25M7 10L11 6M21 26L25 22" stroke="currentColor" stroke-width="3"/></svg><span>"collagene"</span><span class="version">"2"</span></a>
                                <div class="top-divider"></div><span class="workspace-label">{move||workspace_name.get()}</span>
                                <div class="top-actions"><button class="button secondary" on:click=move|_|{modal_error.set(String::new());modal.set("workspace");}>"Workspace"</button><button class="button secondary" on:click=move|_|{modal_error.set(String::new());modal.set("import");}>"Paste sequence"</button><button class="button primary" on:click=move|_|if let Some(input)=file_input.get(){input.click();}>"Open files"</button></div>
                                <input node_ref=file_input type="file" multiple accept=".gb,.gbk,.gbff,.genbank,.gp,.gpt,.genpept,.gcg,.seq,.pep,.uniprot,.swiss,.swissprot,.sbol,.jsonld,.json,.ape,.fa,.fasta,.fna,.faa,.frn,.txt,.cg2,.cg2design,.cg2plan,.dna,.prot,.embl,.emb,.ab1,.abi,.scf,.fastq,.fq,.csv,.tsv,.gff3,.gff,.gtf,.bed,.xml,.collagene,.seqfile,.xdna,.xdgn,.xprt,.aln,.clustal,.sto,.stockholm,.afa,.afasta,.cg2aln,.phy,.phylip,.msf,.nex,.nexus,.pir,.nbrf,.gde,.slx,.selex,.sel,.sam,.bam,.vcf,.bcf,.gz" style="display:none" on:change=open_file/>
                            </header>
                            <div class="body-layout" class:library-hidden=move||!show_library.get()>
                                <nav class="workspace-rail" aria-label="Workspace navigation">
                                    <button class="rail-button" class:active=move||show_library.get() title="Documents" aria-label="Documents" aria-expanded=move||show_library.get().to_string() disabled=move||!layout.ready.get() on:click=move|_|show_library.update(|visible|*visible=!*visible)>
                                        <svg viewBox="0 0 24 24" aria-hidden="true"><path d="M3 5h7l2 2h9v13H3z"/><path d="M3 10h18"/></svg>
                                    </button>
                                    <button class="rail-button" title="Import sequence files" aria-label="Import sequence files" on:click=move|_|if let Some(input)=file_input.get(){input.click();}>
                                        <svg viewBox="0 0 24 24" aria-hidden="true"><path d="M12 3v12m-4-4 4 4 4-4M4 16v5h16v-5"/></svg>
                                    </button>
                                    <button class="rail-button" title="Import from Addgene or NCBI" aria-label="Import from Addgene or NCBI" on:click=move|_|{modal_error.set(String::new());modal.set("url");}>
                                        <svg viewBox="0 0 24 24" aria-hidden="true"><path d="m10 13 4-4m-6 6-1 1a4 4 0 0 1-6-6l4-4a4 4 0 0 1 6 0m2 3 1-1a4 4 0 0 1 6 6l-4 4a4 4 0 0 1-6 0"/></svg>
                                    </button>
                                    <span class="rail-spacer"></span>
                                    <button class="rail-button" title="Workspace files and settings" aria-label="Workspace files and settings" on:click=move|_|{modal_error.set(String::new());modal.set("workspace");}>
                                        <svg viewBox="0 0 24 24" aria-hidden="true"><path d="M4 6h16M4 12h16M4 18h16"/><circle cx="9" cy="6" r="2"/><circle cx="15" cy="12" r="2"/><circle cx="8" cy="18" r="2"/></svg>
                                    </button>
                                </nav>
                                <div class="explorer-frame" id="document-explorer-pane"><Show when=crate::library_view::enabled><crate::library_view::LibraryPane client=crate::library_client::browser() current=library_current bindings=library_bindings capture=library_capture on_open=library_open on_binding=library_record_binding locked=native_catalog_locked refresh=native_catalog_refresh on_folder=native_folder_changed/></Show><Show when=crate::desktop_rpc::is_desktop><crate::catalog_backup_view::CatalogBackupControls client=Callback::new(|command|->crate::catalog_backup_view::LocalResult<collagene_desktop_protocol::CatalogBackupSuccess>{Box::pin(async move{match crate::desktop_rpc::call(collagene_desktop_protocol::Command::CatalogBackup(command)).await?{collagene_desktop_protocol::Success::CatalogBackup(value)=>Ok(value),_=>Err("Unexpected native catalog backup reply.".into())}})}) context=Memo::new(move|_|native_library_context.get()) locked=explorer_locked on_catalog_pending=catalog_pending_changed/><crate::remote_library_view::RemoteLibraryPane client=crate::desktop_rpc::sync_client() native=true saved=native_saved destination=native_destination_memo locked=native_catalog_locked on_local_change=native_catalog_changed/></Show><crate::document_explorer::DocumentExplorer documents=explorer_documents active_document=explorer_active locked=explorer_locked on_activate=explorer_activate on_open=Callback::new(move|()|if let Some(input)=file_input.get(){input.click();})/><Show when=move||!hidden_products.with(|ids|ids.is_empty())><label class="preparation-visibility"><input type="checkbox" prop:checked=move||show_prepared_history.get() on:change=move|e|show_prepared_history.set(event_target::<HtmlInputElement>(&e).checked())/>{move||format!("Show preparation history ({})",hidden_products.with(|ids|ids.len()))}</label></Show><crate::workspace_layout::PaneResize layout=layout pane=crate::workspace_layout::Pane::Documents/></div>
                                <main>
                                    <crate::document_explorer::DocumentTabs documents=tab_items active_document=explorer_active locked=explorer_locked on_activate=explorer_activate on_close=close_document_tab/>
                                    <Show when=move||!has_open_document.get()><section class="empty-document-view" id="no-open-documents" tabindex="-1" aria-label="No open documents"><h2>"No documents open"</h2><p>"Choose a sequence from Documents, or import a file. Closing a tab keeps its document and history in the workspace."</p><div class="empty-document-actions"><button class="button primary" disabled=move||workspace_restoring.get() on:click=move|_|if let Some(input)=file_input.get(){input.click();}>"Open sequence…"</button><button class="button secondary" disabled=move||workspace_restoring.get() on:click=move|_|{modal_error.set(String::new());modal.set("url");}>"Import from URL…"</button></div></section></Show>
                                    <div class="active-document" hidden=move||!has_open_document.get() id="active-sequence-document" role="tabpanel" aria-label="Active sequence document" tabindex="0">
                                    <div class="document-header"><div><div class="eyebrow">"SEQUENCE"</div><h1>{move||molecule.get().name}</h1><div class="document-meta"><span>{move||format!("{} {}",molecule.get().sequence.len(),molecule.get().sequence.alphabet().short_unit())}</span><span>"·"</span><span>{move||format!("{} {}",if molecule.get().topology==Topology::Circular{"Circular"}else{"Linear"},match molecule.get().sequence.alphabet(){Alphabet::Dna=>"DNA",Alphabet::Rna=>"RNA",Alphabet::Protein=>"protein"})}</span><span>"·"</span><span>{move||molecule.get().sequence.gc_fraction().map(|v|format!("{:.1}% GC",v*100.)).unwrap_or_else(||"GC unavailable".into())}</span></div><Show when=move||molecule.get().metadata.contains_key("source:provider")><p class="source-provenance">{move||{let m=molecule.get();let mut parts=vec![format!("Source: {}",m.metadata.get("source:provider").map(String::as_str).unwrap_or("unknown")),m.metadata.get("source:accession_version").or_else(||m.metadata.get("source:identifier")).cloned().unwrap_or_else(||"record".into())];if let Some(id)=m.metadata.get("source:sequence_id"){parts.push(format!("sequence {id}"));}if let Some(category)=m.metadata.get("source:addgene_completeness"){parts.push(format!("{category} record (provider label)"));}parts.push(format!("retrieved {}",m.metadata.get("source:retrieved_at").map(String::as_str).unwrap_or("unknown")));parts.join(" · ")}}</p></Show></div>
                                        <div class="document-actions"><details class="document-file-menu" data-transient-menu=""><summary class="button secondary">"File ▾"</summary><div class="document-file-menu-content"><button class="button secondary" on:click=save_file>"Save .cg2"</button><button class="button secondary" on:click=move|_|{let mol=molecule.get_untracked();if mol.sequence.alphabet()==Alphabet::Protein {modal_error.set(String::new());modal.set("protein_genpept");return;}match formats::export_with_report(&mol,formats::TextFormat::GenBank){Ok(report)=>crate::app::download_with_completion(&format!("{}.gb",mol.name),"text/plain",report.text.as_bytes(),move|result|{if message.try_get_untracked().is_none(){return;}if let Err(e)=result {message.set(e)} else {message.set(report.warnings.join(" "));}}),Err(e)=>message.set(e.to_string())}}>{move||if molecule.with(|m|m.sequence.alphabet()==Alphabet::Protein){"Export GenPept"}else{"Export GenBank"}}</button><hr/><div class="document-extra-file-actions"><Show when=move||molecule.with(|m|m.sequence.alphabet()==Alphabet::Protein)><button class="button secondary" on:click=move|_|{modal_error.set(String::new());modal.set("protein_uniprot");}>"Export UniProt"</button></Show><button class="button secondary" disabled=move||workspace_restoring.get()||!has_open_document.get() on:click=move|_|{modal_error.set(String::new());modal.set("map_svg");}>"Export map SVG"</button><button class="button secondary" on:click=move|_|{modal_error.set(String::new());modal.set("compat_snapgene");}>"Export SnapGene"</button><button class="button secondary" on:click=move|_|{modal_error.set(String::new());modal.set("compat_xdna");}>"Export XDNA"</button><button class="button secondary" on:click=move|_|{modal_error.set(String::new());modal.set("compat_gcg");}>"Export GCG"</button><button class="button secondary" on:click=move|_|{let retained=molecule.with_untracked(|m|m.metadata.get("sbol:source_hash").and_then(|hash|attachments.with_untracked(|assets|assets.get(hash).filter(|a|a.kind==crate::sbol_exchange::ATTACHMENT_KIND&&a.media_type==crate::sbol_exchange::MEDIA_TYPE&&blake3::hash(&a.bytes).to_hex().as_str()==hash).map(|a|std::sync::Arc::new(a.bytes.clone())))));pending_sbol.set(retained);modal_error.set(String::new());modal.set("sbol");}>"Export SBOL"</button><button class="button secondary" on:click=move|_|{modal_error.set(String::new());modal.set("document");}>"Document settings"</button><button class="button secondary" on:click=move|_|{let mol=molecule.get_untracked();match formats::export_with_report(&mol,formats::TextFormat::Fasta){Ok(report)=>crate::app::download_with_completion(&format!("{}.fasta",mol.name),"text/plain",report.text.as_bytes(),move|result|{if message.try_get_untracked().is_none(){return;}if let Err(error)=result {message.set(error);} else {message.set(report.warnings.join(" "));}}),Err(e)=>message.set(e.to_string())}}>"FASTA"</button><Show when=move||molecule.with(|m|m.sequence.alphabet()!=Alphabet::Protein)><button class="button secondary" on:click=move|_|match formats::export_with_report(&molecule.get_untracked(),formats::TextFormat::Embl){Ok(report)=>{crate::app::download_with_completion("sequence.embl","text/plain",report.text.as_bytes(),move|result|{if message.try_get_untracked().is_none(){return;}if let Err(e)=result {message.set(e);} else {message.set(report.warnings.join(" "));}})},Err(e)=>message.set(e.to_string())}>"EMBL"</button></Show><Show when=move||molecule.with(|m|m.metadata.contains_key("protein:original:chunks"))><button class="button secondary" on:click=move|_|{let source=molecule.get_untracked();let guard=export_context_guard(documents,active,workspace_restoring,workspace_epoch);match collagene_core::protein_formats::original_bytes(&source){Ok(bytes)=>crate::app::download_with_completion_guarded("original-protein-input.txt","text/plain",&bytes,move||guard.run(()),move|result|{if message.try_get_untracked().is_some(){match result{Ok(())=>message.set("Original protein source saved; it does not represent subsequent edits.".into()),Err(e)=>message.set(e)}}}),Err(e)=>message.set(e.to_string())}} >"Original protein input"</button></Show><Show when=move||molecule.with(|m|m.original_genbank.is_some())><button class="button secondary" title="Download the exact original GenBank input, before edits" on:click=move|_|molecule.with_untracked(|m|{if let Some(raw)=&m.original_genbank{crate::app::download_with_completion("original-input.gb","text/plain",raw.as_bytes(),move|result|{if message.try_get_untracked().is_none(){return;}if let Err(e)=result {message.set(e);}})}})>"Original GenBank input"</button></Show><Show when=move||molecule.get().metadata.contains_key(snapgene::ORIGINAL_BASE64_KEY)><button class="button secondary" on:click=move|_|match snapgene::original_snapgene_bytes(&molecule.get_untracked()){Ok(bytes)=>{crate::app::download_with_completion("original.dna","application/octet-stream",&bytes,move|result|{if message.try_get_untracked().is_none(){return;}if let Err(e)=result {message.set(e);}})},Err(e)=>message.set(e.to_string())}>"Original .dna"</button></Show><Show when=move||molecule.get().metadata.contains_key("trace_original_base64")><button class="button secondary" on:click=move|_|match traces::original_trace_bytes(&molecule.get_untracked()){Ok(bytes)=>{let extension=if bytes.starts_with(b"ABIF"){"ab1"}else{"scf"};crate::app::download_with_completion(&format!("original.{extension}"),"application/octet-stream",&bytes,move|result|{if message.try_get_untracked().is_none(){return;}if let Err(e)=result {message.set(e);}})},Err(e)=>message.set(e.to_string())}>"Original trace"</button></Show><Show when=move||molecule.with(|doc|doc.metadata.contains_key(collagene_core::legacy_formats::ORIGINAL_BASE64_KEY))><button class="button secondary" on:click=move|_|match collagene_core::legacy_formats::original_legacy_bytes(&molecule.get_untracked()){Ok(bytes)=>{let extension=if collagene_core::legacy_formats::is_collagene_xml(&bytes){"xml"}else{"xdna"};crate::app::download_with_completion(&format!("original-legacy.{extension}"),"application/octet-stream",&bytes,move|result|{if message.try_get_untracked().is_none(){return;}if let Err(error)=result {message.set(error);}})},Err(error)=>message.set(error.to_string())}>"Original legacy file"</button></Show><button class="button secondary" on:click=move|_|{let current=molecule.get_untracked();let rendered=if let Some(sites)=restrictions.displayed.get_untracked(){svg::render_map_with_restriction_sites(&current,&current.revision(),&sites)}else{svg::render_map(&current)};match rendered{Ok(text)=>{crate::app::download_with_completion("sequence-map.svg","image/svg+xml",text.as_bytes(),move|result|{if message.try_get_untracked().is_none(){return;}if let Err(e)=result {message.set(e);}})},Err(e)=>message.set(e.to_string())}}>"SVG map"</button></div></div></details></div>
                                    </div>
                                    <crate::digest::DigestProvenance molecule=molecule/>
                                    <Show when=move||molecule.with(|m|m.metadata.contains_key("crispr:producer_json"))><div class="crispr-historical"><strong>"CRISPR producer provenance · unverified on opening"</strong><p>"This saved sequence design does not restore a current source/guide/target check or prove cellular editing. Recompute against the exact original source and independent target before reuse."</p><button class="button secondary" on:click=move|_|molecule.with_untracked(|m|{if let Some(json)=m.metadata.get("crispr:producer_json"){crate::app::download_with_completion("crispr-retained-provenance.json","application/json",json.as_bytes(),move|result|{if message.try_get_untracked().is_none(){return;}if let Err(e)=result {message.set(e);}})}})>"Download retained CRISPR provenance"</button></div></Show>
                                    <div class="editor-toolbar"><div class="tool-group"><button class="icon-button" title="Show documents" disabled=move||!layout.ready.get() on:click=move|_|show_library.update(|v|*v=!*v)>"☷"</button><button class="icon-button" title="Undo" disabled=move||documents.with(|docs|!docs[active.get()].can_undo()) on:click=move|_|documents.update(|docs|{let d=&mut docs[active.get_untracked()];if let Err(error)=d.undo(){message.set(error.to_string());}})>"↶"</button><button class="icon-button" title="Redo" disabled=move||documents.with(|docs|!docs[active.get()].can_redo()) on:click=move|_|documents.update(|docs|{let d=&mut docs[active.get_untracked()];if let Err(error)=d.redo(){message.set(error.to_string());}})>"↷"</button><span class="tool-separator"></span><button class="text-tool" disabled=move||selection.get().is_none()||inspected_feature.get().is_some_and(|f|f.location.spans().len()>1) on:click=move|_|{edit_text.set(selection.get_untracked().and_then(|(a,b)|molecule.get_untracked().sequence.as_str().get(a..b).map(str::to_owned)).unwrap_or_default());modal_error.set(String::new());modal.set("edit");}>"Edit selection"</button></div>
                                        <div class="viewer-command-menus"><details class="toolbar-action-menu" data-transient-menu=""><summary>"Create ▾"</summary><div class="toolbar-action-menu-content"><button class="text-tool" on:click=move |_|{feature_for_editor.set(None);modal.set("feature");}>"Annotation…"</button><button class="text-tool" on:click=move |_|panel.set("pcr")>"PCR product…"</button><button class="text-tool" on:click=move |_|assembly_visible.set(true)>"Assembly…"</button></div></details><details class="toolbar-action-menu" data-transient-menu=""><summary>"Analyze ▾"</summary><div class="toolbar-action-menu-content"><button class="button secondary" disabled=move||workspace_restoring.get()||molecule.with(|doc|doc.sequence.alphabet()!=Alphabet::Dna) on:click=move|_|{read_reference.set(Some(primer_context.get_untracked()));panel.set("reads");}>"Check sequencing reads"</button><button class="text-tool" on:click=move |_|panel.set("biology")>"Translation and sequence analysis…"</button><button class="text-tool" on:click=move |_|panel.set("compare")>"Compare sequences…"</button><button class="text-tool" on:click=move |_|panel.set("enzymes")>"Virtual digest…"</button></div></details></div>
                                        <crate::clipboard::ClipboardTools source=Memo::new(move|_|documents.with(|docs|recombination_source_info(&docs[active.get()]))) internal_clipboard source_selection=reference_copy_selection on_preview=preview_reference_paste molecule=molecule selection=selection selected_feature=selected_feature highlighted_spans locked=workspace_restoring epoch=workspace_epoch identity=Memo::new(move|_|documents.with(|docs|docs.get(active.get()).map(|doc|doc.id.clone()).unwrap_or_default())) on_commit=commit_annotations message=message/>
                                        <crate::search_view::SearchTools molecule=molecule identity=Memo::new(move|_|documents.with(|docs|docs[active.get()].id.clone())) epoch=workspace_epoch selection=selection selected_feature=selected_feature regions=search_regions message=message/>
                                    </div>
                                    <Show when=move||!native_export_notice.get().is_empty()><div class="notice" role="status"><span>{move||native_export_notice.get()}</span><button class="icon-button" aria-label="Dismiss saved export receipt" on:click=move|_|native_export_notice.set(String::new())>"×"</button></div></Show>
                                    <Show when=move||!message.get().is_empty()><div class="notice" role="status"><span>{move||message.get()}</span><button class="icon-button" aria-label="Dismiss message" on:click=move|_|message.set(String::new())>"×"</button></div></Show>
                                    <Show when=move||molecule.with(|m|m.metadata.get("collagene:example").is_some_and(|value|value=="interface-demonstration"))><div class="demo-note mvp-example"><div><strong>"Explore the example, or open your sequence"</strong><span>"Synthetic annotations illustrate the viewer. Open a GenBank file to inspect your own sequence and simulate a digest."</span></div><button class="button secondary" on:click=move|_|if let Some(input)=file_input.get(){input.click();}>"Open a GenBank file"</button></div></Show>
                                    <div class="editor-view-bar"><div class="view-switcher" role="group" aria-label="Editor view">
                                        <button aria-label="Split view" disabled=move||!layout.ready.get() aria-pressed=move||(editor_view.get()=="split").to_string() class:active=move||editor_view.get()=="split" on:click=move|_|editor_view.set("split")>"Split"</button>
                                        <button aria-label="Sequence view" disabled=move||!layout.ready.get() aria-pressed=move||(editor_view.get()=="sequence").to_string() class:active=move||editor_view.get()=="sequence" on:click=move|_|editor_view.set("sequence")>"Sequence"</button>
                                        <button aria-label="Map view" disabled=move||!layout.ready.get() aria-pressed=move||(editor_view.get()=="map").to_string() class:active=move||editor_view.get()=="map" on:click=move|_|editor_view.set("map")>"Map"</button>
                                    </div><div class="restriction-quick-controls"><label><input type="checkbox" aria-label="Show restriction enzymes in editor" prop:checked=move||restrictions.visible.get() disabled=move||molecule.with(|m|m.sequence.alphabet()!=Alphabet::Dna) on:change=move|e|restrictions.visible.set(event_target_checked(&e))/><span>"RE sites"</span></label><button class="text-tool" on:click=move|_|panel.set("restriction")>"Choose enzymes"</button><span class="restriction-quick-count">{move||restrictions.displayed.with(|v|v.as_ref().map(|sites|format!("{} shown",sites.len())).unwrap_or_else(||if restrictions.busy.get(){"Scanning…".into()}else{String::new()}))}</span></div></div>
                                    <crate::selection_toolbar::SelectionToolbar molecule=molecule selection=selection selected_feature=selected_feature highlighted_spans=highlighted_spans locked=Memo::new(move|_|workspace_restoring.get())
                                        on_annotate=Callback::new(move|()|{
                                            let location=molecule.with_untracked(|m|selected_feature.get_untracked().and_then(|i|m.features.get(i).map(|f|f.location.clone()))).or_else(||{let parts=highlighted_spans.get_untracked();if parts.is_empty()||parts.len()>4096{return None;}Some(if parts.len()==1{Location::range(parts[0].0,parts[0].1)}else{Location::Join(parts.into_iter().map(|(a,b)|Location::range(a,b)).collect())})});
                                            if let Some(location)=location{feature_for_editor.set(None);feature_initial_location.set(Some(location));modal_error.set(String::new());modal.set("feature");}
                                        })
                                        on_mutation=Callback::new(move|()|match mutation_selection.run(()){Ok(seed)=>{mutation_seed.set(Some(seed));panel.set("mutagenesis");},Err(error)=>message.set(error)})
                                        on_domestication=Callback::new(move|()|panel.set("domestication"))
                                        on_primer=Callback::new(move|reverse:bool|{match pair_selection.run(()).and_then(|selection|selection.oligo(if reverse{"Reverse"}else{"Forward"},reverse,"")){Ok(oligo)=>{match oligo.primer(){Ok(primer)=>{pair_draft.update(|draft|if reverse{draft.reverse=oligo}else{draft.forward=oligo});primer_seed.set(Some(primer));panel.set("primers");},Err(error)=>message.set(error)}},Err(error)=>message.set(error)}})
                                        on_edit=Callback::new(move|()|{if let Some((a,b))=selection.get_untracked(){edit_text.set(molecule.with_untracked(|m|m.sequence.as_str().get(a..b).unwrap_or_default().to_owned()));modal_error.set(String::new());modal.set("edit");}})
                                        on_clear=Callback::new(move|()|batch(||{selection.set(None);selected_feature.set(None);search_regions.set(None);restrictions.selected.set(None);}))
                                    />
                                    <div class="editor-grid" class:sequence-only=move||editor_view.get()=="sequence" class:map-only=move||editor_view.get()=="map">
                                        <section class="sequence-panel"><div class="panel-heading"><div><h2>"Sequence"</h2></div><span class="panel-subtitle">{move||molecule.with(|doc|if doc.sequence.alphabet()==Alphabet::Protein{"Amino-acid reference"}else{"Reference & complement"})}</span></div><Show when=move||molecule.with(|doc|doc.metadata.contains_key("legacy:upper_aligned"))><details class="legacy-source-rows"><summary>"Original aligned source rows"</summary><p class="helper">"Upper is stored 5′→3′; lower is aligned 3′→5′. These original rows may describe an earlier document revision. The editor displays the reference and its calculated complement."</p><pre>{move||molecule.with(|doc|format!("{}\n{}",doc.metadata.get("legacy:upper_aligned").map(|value|value.chars().take(2000).collect::<String>()).unwrap_or_default(),doc.metadata.get("legacy:lower_aligned").map(|value|value.chars().take(2000).collect::<String>()).unwrap_or_default()))}</pre><p class="helper">"First 2,000 columns shown; original legacy-file download retains all source evidence."</p></details></Show><Show when=move||molecule.get().metadata.contains_key("trace_original_base64")><button class="text-tool" on:click=move|_|trace_view.update(|v|*v=!*v)>{move||if trace_view.get(){"Show sequence"}else{"Show chromatogram"}}</button></Show><Show when=move||trace_view.get()&&molecule.get().metadata.contains_key("trace_original_base64") fallback=move||view!{<crate::sequence_view::SequenceView pcr_overlay=Some(pcr_overlay) hidden_features=Some(annotation_hidden) molecule=molecule selection=selection selected_feature=selected_feature highlighted_spans=highlighted_spans search_regions=search_regions on_show_features=Callback::new(move|()|panel.set("features")) identity=Memo::new(move|_|documents.with(|docs|format!("{}:{}:{}",workspace_epoch.get(),docs[active.get()].id,docs[active.get()].head_event()))) pcr_context=Some(restriction_identity) locked=Memo::new(move|_|workspace_restoring.get()) on_message=Callback::new(move|text:String|message.set(text)) on_history=sequence_on_history restriction_sites=restrictions.displayed restriction_visible=restrictions.visible on_show_restrictions=Callback::new(move|()|panel.set("restriction")) restriction_focus=Some(hovered_restriction) on_select_restriction=navigate_restriction on_edit=sequence_on_edit/>}><crate::trace_view::TraceView molecule=molecule selection=selection source_scope=Memo::new(move |_| {
        let event = documents.with(|docs| docs.get(active.get()).map(current_event_binding));
        format!("{}:{}", workspace_epoch.get(), serde_json::to_string(&event).unwrap_or_default())
    })/></Show></section>
                                        <crate::workspace_layout::PaneResize layout=layout pane=crate::workspace_layout::Pane::Sequence/>
                                        <section class="map-panel"><div class="panel-heading"><div><h2>{move||if molecule.get().topology==Topology::Circular{"Plasmid map"}else{"Linear map"}}</h2></div><div class="map-controls"><button class="icon-button" title="Zoom out" on:click=move|_|if let Some(r)=renderer.get_value(){r.zoom(1./1.15);}>"−"</button><button class="icon-button" title="Reset map view" on:click=move|_|if let Some(r)=renderer.get_value(){r.reset();}>"Fit"</button><button class="icon-button" title="Zoom in" on:click=move|_|if let Some(r)=renderer.get_value(){r.zoom(1.15);}>"+"</button></div></div>
                                            <div class="map-stage"><canvas node_ref=gl_canvas aria-label=move||if molecule.with(|m|m.topology==Topology::Circular){"Interactive circular sequence overview; drag to rotate, scroll to zoom, click annotations to inspect"}else{"Interactive linear sequence overview; scroll to zoom, click annotations to inspect"} on:pointerdown=move|e|{if let Some(canvas)=gl_canvas.get(){let _=canvas.set_pointer_capture(e.pointer_id());}if let Some(r)=renderer.get_value(){r.pointer_down(e.client_x()as f64,e.client_y()as f64);}} on:pointermove=move|e|if let Some(r)=renderer.get_value(){r.pointer_move(e.client_x()as f64,e.client_y()as f64);} on:pointerup=move|e|if let Some(r)=renderer.get_value(){if let Some(hit)=r.pointer_up(e.client_x()as f64,e.client_y()as f64){if let Some(index)=hit.restriction_index{if let Some(site)=restrictions.scan.with_untracked(|s|s.as_ref().and_then(|s|s.sites.iter().find(|site|site.index==index).cloned())){navigate_restriction.run(site);}}else{batch(||{restrictions.selected.set(None);search_regions.set(None);selected_feature.set(hit.feature_index);selection.set(if hit.spans.len()==1{Some(hit.span)}else{None});if hit.spans.len()>1{search_regions.set(Some((molecule_revision.get_untracked(),hit.spans.clone())));}if hit.feature_index.is_some(){panel.set("features");}});}}} on:pointercancel=move|_|if let Some(r)=renderer.get_value(){r.pointer_cancel();} on:lostpointercapture=move|_|if let Some(r)=renderer.get_value(){r.pointer_cancel();} on:wheel=move|e|{e.prevent_default();if let Some(r)=renderer.get_value(){r.zoom(if e.delta_y()>0.{0.95}else{1.05});}}></canvas><canvas node_ref=label_canvas class="map-labels" aria-hidden="true"></canvas><p class="map-geometry-note" role="status">"Overview is simplified. Browse the annotation pages or export the full annotation set."</p><p class="restriction-map-note" role="status">"Some restriction labels are hidden for space. Open RE sites to inspect the selected positions."</p>
                                                <Show when=move||!render_error.get().is_empty()><div class="map-error">{move||render_error.get()}</div></Show>
                                            </div><div class="map-caption"><span>{move||if molecule.get().topology==Topology::Circular{"Drag to rotate · scroll to zoom · select exact ranges in Sequence"}else{"Scroll to zoom · click features or sites to inspect"}}</span></div>
                                            <Show when=move||{let stats=map_labels.get();stats.omitted>0||stats.abbreviated>0||stats.limited}><div class="map-label-note"><span>{move||{let stats=map_labels.get();let mut parts=vec![format!("{} of {} annotation labels shown",stats.displayed,stats.total)];if stats.abbreviated>0{parts.push(format!("{} abbreviated",stats.abbreviated));}if stats.unplaced>0{parts.push(format!("{} without local coordinates",stats.unplaced));}if stats.lane_overflow>0{parts.push(format!("{} overlapping features hidden (3 map lanes)",stats.lane_overflow));}if stats.limited{parts.push("display limit reached".into());}parts.join(" · ")}}</span><button class="text-tool" on:click=move|_|panel.set("features")>"Inspect annotations"</button></div></Show>
                                            <div class="legend"><Show when=move||restrictions.visible.get()><span><i style="background:#8b5ca4"></i>"RE recognition"</span><span><i style="background:#2563a6"></i>"Top cut"</span><span><i style="background:#b45c19"></i>"Bottom cut"</span></Show><span><i style="background:#2b8775"></i>"Coding sequence"</span><span><i style="background:#d6933b"></i>"Promoter"</span><span><i style="background:#666bb3"></i>"Origin"</span></div>
                                        </section>
                                        <aside class="inspector viewer-inspector" class:is-open=move||!panel.get().is_empty()><crate::desktop_shell::InspectorNavigation panel=panel on_assemble=Callback::new(move|()|assembly_visible.set(true)) on_auxiliary=Callback::new(move|action|{use crate::desktop_shell::AuxiliaryAction;modal_error.set(String::new());match action{AuxiliaryAction::GoldenGateWorkflow=>{method_config.update(|config|config.mode=crate::cloning_methods::Mode::GoldenGate);panel.set("methods");},AuxiliaryAction::ImportUrl=>modal.set("url"),AuxiliaryAction::Alignments=>{pending_alignment.set(None);modal.set("alignments");},AuxiliaryAction::Oligos=>{existing_oligo.set(None);pending_oligo.set(None);modal.set("oligos");}}})/><div class="inspector-body"><div class="inspector-drawer-heading"><h2>{move||crate::desktop_shell::panel_title(panel.get())}</h2><button class="icon-button" aria-label="Close tool panel" title="Close tool panel" on:click=move |_|panel.set("")>"×"</button></div>
                                            <Show when=move||panel.get()=="restriction"><crate::restriction_view::RestrictionPanel state=restrictions on_navigate=navigate_restriction saved=saved_enzyme_sets on_save=save_enzyme_set on_delete=delete_enzyme_set/></Show>
                                            <Show when=move||panel.get()=="library"><div class="inspector-content"><crate::feature_tools::FeatureLibraryPanel molecule=molecule selected_feature=selected_feature libraries=libraries on_keep=keep_library on_commit=commit_annotations selection=selection/></div></Show>
                                            <Show when=move||panel.get()=="annotation_files"><div class="inspector-content"><crate::feature_tools::AnnotationFilesPanel molecule=molecule source=annotation_source format=annotation_format on_commit=commit_annotations selection=selection selected_feature=selected_feature/></div></Show>
                                            <Show when=move||panel.get()=="history"><div class="inspector-content"><crate::variant_source_library::VariantSourceLibrary entries=Memo::new(move|_|attachments.with(|assets|assets.iter().filter_map(|(hash,asset)|collagene_core::variants_archive::inspect(asset).ok().map(|info|(hash.clone(),info))).collect())) locked=Memo::new(move|_|workspace_restoring.get()) on_open=Callback::new(move|hash:String|{if workspace_restoring.get_untracked(){return;}let value=attachments.with_untracked(|assets|assets.get(&hash).and_then(|asset|collagene_core::variants_archive::inspect(asset).ok().map(|info|(std::sync::Arc::new(asset.bytes.clone()),info.format,info.name))));if let Some(value)=value{pending_variants.set(Some(value));modal_error.set(String::new());modal.set("variants");}}) on_remove=Callback::new(move|hash:String|{if !workspace_restoring.get_untracked(){attachments.update(|assets|{if assets.get(&hash).is_some_and(|asset|collagene_core::variants_archive::inspect(asset).is_ok()){assets.remove(&hash);}});}})/><crate::read_source_library::ReadSourceLibrary entries=Memo::new(move|_|attachments.with(|assets|assets.iter().filter_map(|(hash,asset)|collagene_core::read_alignment_archive::inspect(asset).ok().map(|info|(hash.clone(),info))).collect())) locked=Memo::new(move|_|workspace_restoring.get()) on_open=Callback::new(move|hash:String|{if workspace_restoring.get_untracked(){return;}let source=attachments.with_untracked(|assets|assets.get(&hash).and_then(|asset|collagene_core::read_alignment_archive::inspect(asset).ok().map(|info|(std::sync::Arc::new(asset.bytes.clone()),match info.format{collagene_core::read_alignment_store::Format::Sam=>collagene_core::read_alignment_session::FileFormat::Sam,collagene_core::read_alignment_store::Format::Bam=>collagene_core::read_alignment_session::FileFormat::Bam},info.name))));if let Some(source)=source{pending_read_alignment.set(Some(source));modal_error.set(String::new());modal.set("read_alignments");}}) on_remove=Callback::new(move|hash:String|{if !workspace_restoring.get_untracked(){attachments.update(|assets|{if assets.get(&hash).is_some_and(|a|collagene_core::read_alignment_archive::inspect(a).is_ok()){assets.remove(&hash);}});}})/><crate::history_view::HistoryPanel molecule=molecule events=Memo::new_owning(move|_|(documents.with(|docs|docs[active.get()].summaries()),true)) current_event=Memo::new(move|_|documents.with(|docs|docs[active.get()].head_event().to_string())) navigation=Memo::new(move|_|documents.with(|docs|docs[active.get()].navigation_event_ids())) external_events=Memo::new(move|_|documents.with(|docs|docs.iter().enumerate().filter(|(i,_)|*i!=active.get()).flat_map(|(_,doc)|{let label=format!("[{}] {}",doc.id,doc.molecule.name.chars().take(100).collect::<String>());doc.summaries().into_iter().map(move|e|(e.event_id,label.clone(),e.revision))}).collect())) on_restore_revision=restore_revision on_fork_revision=fork_revision candidates=Memo::new(move|_|{let current=candidate.get().map(|c|c.plan_revision);stored_candidates.with(|values|values.iter().map(|c|(c.plan_revision.clone(),c.plan.name.clone(),c.plan.target_revision.clone(),current.as_ref()==Some(&c.plan_revision))).collect())}) on_open_candidate=Callback::new(move|revision:String|{if workspace_restoring.get_untracked(){return;}if let Some(value)=stored_candidates.with_untracked(|values|values.iter().find(|c|c.plan_revision==revision).cloned()){candidate_events.set(None);candidate.set(Some(value));plan_report.set(None);last_plan.set(None);last_plan_target.set(None);plan_generation.update(|g|*g=g.wrapping_add(1));panel.set("plan");}}) on_remove_candidate=Callback::new(move|revision:String|{if workspace_restoring.get_untracked(){return;}batch(||{if candidate.with_untracked(|c|c.as_ref().is_some_and(|c|c.plan_revision==revision)){candidate_events.set(None);candidate.set(None);plan_report.set(None);last_plan.set(None);last_plan_target.set(None);plan_generation.update(|g|*g=g.wrapping_add(1));}stored_candidates.update(|values|values.retain(|c|c.plan_revision!=revision));});}) collections=Memo::new(move|_|oligo_tables.with(|tables|tables.iter().enumerate().map(|(i,t)|(i,format!("Oligo collection {}",t.source_blake3.chars().take(12).collect::<String>()),t.records.len())).collect())) on_open_oligos=Callback::new(move|index:usize|{if let Some(table)=oligo_tables.with_untracked(|tables|tables.get(index).cloned()){existing_oligo.set(Some(table.as_ref().clone()));pending_oligo.set(None);modal_error.set(String::new());modal.set("oligos");}})/><Show when=move||!recorded_runs.with(|runs|runs.is_empty())><section class="recorded-runs"><h3>"Recorded reconstructions"</h3><p class="helper">"Saved computational history. Each record needs a fresh replay before its verification can be used."</p>{move||recorded_runs.with(|runs|runs.iter().skip(run_page.get()*20).take(20).map(|run|{let id=run.run_id.clone();view!{<article class="run-product"><strong>{run.plan_name.clone()}</strong><p class="helper">{format!("{} products · recorded {:?} · {} · unverified",run.products.len(),run.recorded_status,run.recorded_engine)}</p><button class="button secondary" on:click=move|_|{selected_run.set(id.clone());modal_error.set(String::new());modal.set("run_details");}>"View cloning run"</button></article>}}).collect_view())}<div class="pagination"><button class="button secondary" disabled=move||run_page.get()==0 on:click=move|_|run_page.update(|p|*p=p.saturating_sub(1))>"Previous records"</button><button class="button secondary" disabled=move||{(run_page.get()+1)*20>=recorded_runs.with(|r|r.len())} on:click=move|_|run_page.update(|p|*p+=1)>"Next records"</button></div></section></Show><Show when=move||!saved_assemblies.with(|rows|rows.is_empty())><section class="saved-assemblies"><h3>"Saved read assemblies"</h3><p class="helper">"Historical evidence only; replay exact retained reads before using a consensus."</p>{move||saved_assemblies.with(|rows|rows.iter().skip(assembly_page.get()*20).take(20).map(|row|{let hash=row.attachment_hash.clone();view!{<article><p>{format!("{} reads · {} contigs · recorded complete {} · {}",row.read_count,row.contig_count,row.recorded_complete,row.engine_version)}</p><button class="button secondary" on:click=move|_|{inspected_assembly.set(hash.clone());modal_error.set(String::new());modal.set("saved_contig");}>"View saved read assembly"</button></article>}}).collect_view())}<button class="button secondary" disabled=move||assembly_page.get()==0 on:click=move|_|assembly_page.update(|p|*p=p.saturating_sub(1))>"Previous saved assemblies"</button><button class="button secondary" disabled=move||{(assembly_page.get()+1)*20>=saved_assemblies.with(Vec::len)} on:click=move|_|assembly_page.update(|p|*p+=1)>"Next saved assemblies"</button></section></Show></div></Show>
                                            <Show when=move||panel.get()=="reads"><div class="inspector-content"><h2>"Read validation"</h2><Show when=move||molecule.with(|m|m.metadata.contains_key(oligo_formats::FASTQ_ORIGINAL_KEY))><section><h3>"FASTQ quality evidence"</h3><p class="helper">{move||match oligo_formats::quality_evidence(&molecule.get()){Ok(e)=>format!("{} calls · {:?} · {}. {}",e.qualities.len(),e.encoding,if e.current_valid{"Quality mapping matches current calls"}else{"Current calls changed: quality mapping is stale"},e.warnings.join(" ")),Err(e)=>e.to_string()}}</p><label>"FASTQ output encoding"<select aria-label="FASTQ output encoding" prop:value=move||match fastq_encoding.get(){None=>"",Some(oligo_formats::QualityEncoding::Phred33)=>"33",Some(oligo_formats::QualityEncoding::Phred64)=>"64"} on:change=move|e|fastq_encoding.set(match event_target_value(&e).as_str(){"33"=>Some(oligo_formats::QualityEncoding::Phred33),"64"=>Some(oligo_formats::QualityEncoding::Phred64),_=>None})><option value="">"Choose encoding"</option><option value="33">"Phred+33"</option><option value="64">"Phred+64"</option></select></label><button class="button secondary full" disabled=move||fastq_encoding.get().is_none() on:click=move|_|{let Some(encoding)=fastq_encoding.get_untracked()else{return;};match oligo_formats::export_fastq(&molecule.get_untracked(),encoding){Ok(report)=>crate::app::download_with_completion("current.fastq","text/plain",report.text.as_bytes(),move|result|{if message.try_get_untracked().is_none(){return;}match result{Ok(())=>message.set(report.warnings.join(" ")),Err(e)=>message.set(e)}}),Err(e)=>message.set(e.to_string())}}>"Export current FASTQ"</button><button class="button secondary full" on:click=move|_|match oligo_formats::original_fastq_bytes(&molecule.get_untracked()){Ok(bytes)=>crate::app::download_with_completion("original.fastq","text/plain",&bytes,move|result|{if message.try_get_untracked().is_none(){return;}if let Err(e)=result {message.set(e)}}),Err(e)=>message.set(e.to_string())}>"Original FASTQ record"</button><button class="button secondary full" disabled=move||selection.get().is_none() on:click=move|_|{let Some((a,b))=selection.get_untracked()else{return;};match oligo_formats::trim_fastq(&molecule.get_untracked(),a,b){Ok(doc)=>if let Err(e)=add_documents(documents,active,workspace_restoring,vec![doc]){message.set(e)},Err(e)=>message.set(e.to_string())}}>"Open selected FASTQ calls"</button></section></Show><label>"Reference construct"<select aria-label="Read check reference construct" prop:value=move||read_check_source.get().binding.document_id on:change=move|e|{let id=event_target_value(&e);if let Some(binding)=documents.with_untracked(|docs|docs.iter().find(|doc|doc.id==id).map(current_event_binding)){read_reference.set(Some(binding));}}>{move||documents.with(|docs|docs.iter().filter(|doc|doc.molecule.sequence.alphabet()==Alphabet::Dna).map(|doc|view!{<option value=doc.id.clone()>{doc.molecule.name.clone()}</option>}).collect_view())}</select></label>
                                            <crate::read_checks_view::SavedReadChecksPanel on_review=open_read_review event_available=primer_anchor_available attachment_available=read_check_attachment_available reference=read_check_reference snapshot=contig_snapshot epoch=workspace_epoch on_read=navigate_read_call on_reference=navigate_check_reference/>
                                            <crate::read_review_history::ReadReviewHistory reference=read_check_reference snapshot=contig_snapshot epoch=workspace_epoch on_open=reopen_read_review on_reference=navigate_check_reference/>
                                            <crate::trace_validation::TraceValidationPanel molecule=read_check_molecule reference_info=read_check_source source_info=Callback::new(move|i:usize|documents.with(|docs|docs.get(i).map(recombination_source_info))) epoch=workspace_epoch on_read=navigate_read_call on_record=record_read_check on_reference=navigate_check_reference query_options=Memo::new(move|_|documents.with(|docs|docs.iter().enumerate().map(|(i,d)|(i,d.molecule.name.clone())).collect())) get_query=Callback::new(move|i:usize|documents.with(|docs|docs.get(i).map(|d|d.molecule.clone()))) selection=selection selected_feature=selected_feature on_consensus=Callback::new(move|doc:Molecule|add_documents(documents,active,workspace_restoring,vec![doc]))/></div></Show>
                                            <Show when=move||panel.get()=="features"><div class="inspector-content"><h2>"Annotations"<span class="count">{move||molecule.get().features.len()}</span></h2><p class="helper">"Select an annotation to inspect all of its parts and biological sequence."</p><crate::annotation_table::AnnotationTable molecule context=restriction_identity locked=Memo::new(move|_|workspace_restoring.get()) selected_feature hidden=annotation_hidden on_visibility=set_annotation_visibility on_apply=apply_annotation_batch on_inspect=Callback::new(move|index:usize|{molecule.with_untracked(|m|{if let Some(feature)=m.features.get(index){batch(||{selected_feature.set(Some(index));selection.set(feature.location.spans().first().map(|(a,b,_)|(*a,*b)));search_regions.set(None);});}});}) on_edit=Callback::new(move|index:Option<usize>|{feature_for_editor.set(index);feature_initial_location.set(None);modal.set("feature");})/>{move||inspected_feature.get().map(|f|{let location=editing::location_text(&f.location).unwrap_or_else(|_|"unresolved".into());let size=f.location.spans().iter().try_fold(0usize,|n,(a,b,_)|n.checked_add(b-a)).unwrap_or(usize::MAX);let exact=if size<=10_000 {f.location.extract(&molecule.get().sequence).unwrap_or_else(|e|e.to_string())}else{"Sequence exceeds this inspector's 10,000-symbol display limit; use a region/sequence export.".into()};view!{<div class="annotation-details"><h3>{f.label().to_string()}</h3><button class="text-tool" on:click=move|_|panel.set("library")>"Use selected annotation in library"</button><p>{location}</p><pre class="annotation-sequence">{exact}</pre><p class="helper">"Compound selections highlight every span. Choose contiguous bases to edit the underlying sequence."</p><button class="button secondary" on:click=move|_|{feature_for_editor.set(selected_feature.get_untracked());modal.set("feature");}>"Edit annotation"</button><button class="text-tool" on:click=move|_|if let Some(i)=selected_feature.get_untracked(){let current=molecule.get_untracked();match editing::delete_feature(&current,i){Ok(edited)=>if let Err(e)=commit_document(documents,active,workspace_restoring,&current.revision(),edited){message.set(e);},Err(e)=>message.set(e.to_string())}}>"Delete annotation"</button><details><summary>"Annotation qualifiers"</summary><dl>{f.qualifiers.into_iter().take(200).map(|(key,value)|view!{<dt>{key}</dt><dd>{value.unwrap_or_else(||"flag".into())}</dd>}).collect_view()}</dl></details></div>}})}</div></Show>
                                            <Show when=move||panel.get()=="enzymes"><crate::digest::DigestPanel origin=digest_run_origin get_workspace=contig_snapshot saved=saved_digests on_save=save_digest_setup on_load=load_digest_setup on_delete=delete_digest_setup on_gel=queue_digest_gel on_fork=fork_digest_source source=Memo::new(move|_|documents.with(|docs|{let d=&docs[active.get()];crate::digest::SourceInfo{binding:collagene_core::run_lineage::EventBinding{document_id:d.id.clone(),event_id:d.head_event().to_string(),revision:molecule_revision.get()},name:d.molecule.name.clone(),length:d.molecule.sequence.len(),alphabet:d.molecule.sequence.alphabet(),topology:d.molecule.topology}})) get_molecule=Callback::new(move|()|molecule.get_untracked()) epoch=workspace_epoch on_navigate=navigate_digest on_open=open_digest on_assembly=use_digest_in_assembly boundary_pick=digest_boundary_pick initial_enzymes=restrictions.names.get_untracked() on_methylation=Callback::new(move|()|panel.set("methylation"))/></Show>
                                            <Show when=move||panel.get()=="pcr"><div class="inspector-content"><h2>"Simulate PCR"</h2><p class="helper">"Enter both primers 5′→3′. Annealing must match exactly; tails are introduced into the product."</p><crate::primer_picker::PrimerPicker label="Forward".to_string() sources=primer_sources list_primers=list_saved_primers resolve_primer=resolve_saved_primer epoch=workspace_epoch on_copy=Callback::new(move|primer:Primer|batch(||{pcr_forward_name.set(primer.name);pcr_forward.set(primer.annealing.as_str().into());pcr_forward_tail.set(primer.tail);pcr_generation.update(|g|*g=g.wrapping_add(1));}))/>
            <crate::primer_picker::PrimerPicker label="Reverse".to_string() sources=primer_sources list_primers=list_saved_primers resolve_primer=resolve_saved_primer epoch=workspace_epoch on_copy=Callback::new(move|primer:Primer|batch(||{pcr_reverse_name.set(primer.name);pcr_reverse.set(primer.annealing.as_str().into());pcr_reverse_tail.set(primer.tail);pcr_generation.update(|g|*g=g.wrapping_add(1));}))/>
            <label>"Forward name"<input maxlength="256" prop:value=move||pcr_forward_name.get() on:input=move|e|pcr_forward_name.set(event_target_value(&e))/></label>
            <crate::run_origin_view::SourceOriginPicker source=Memo::new(move |_|documents.with(|docs|docs.get(active.get()).map(|doc|crate::source_identity::SourceInfo{binding:current_event_binding(doc),name:doc.molecule.name.clone(),length:doc.molecule.sequence.len(),alphabet:doc.molecule.sequence.alphabet(),topology:doc.molecule.topology}))) epoch=workspace_epoch get_workspace=contig_snapshot selected=Memo::new(move |_|pcr_origin.get()) on_select=Callback::new(move |origin:Option<crate::run_origin_model::PinnedOrigin>|{if origin.as_ref().is_some_and(|origin|documents.with_untracked(|docs|docs.get(active.get_untracked()).is_none_or(|doc|current_event_binding(doc)!=*origin.event()))){return Err("PCR source changed; select its material origin again.".into());}pcr_origin.set(origin);Ok(())})/>
            <label>"Reverse name"<input maxlength="256" prop:value=move||pcr_reverse_name.get() on:input=move|e|pcr_reverse_name.set(event_target_value(&e))/></label><label>"Forward annealing"<input prop:value=move||pcr_forward.get() on:input=move|e|pcr_forward.set(event_target_value(&e)) placeholder="ACGT…"/></label><label>"Forward 5′ tail"<input prop:value=move||pcr_forward_tail.get() on:input=move|e|pcr_forward_tail.set(event_target_value(&e)) placeholder="Optional"/></label><label>"Reverse annealing"<input prop:value=move||pcr_reverse.get() on:input=move|e|pcr_reverse.set(event_target_value(&e)) placeholder="ACGT…"/></label><label>"Reverse 5′ tail"<input prop:value=move||pcr_reverse_tail.get() on:input=move|e|pcr_reverse_tail.set(event_target_value(&e)) placeholder="Optional"/></label><button class="button primary full" on:click=do_pcr>"Compute products"</button><p class="helper">"Review reported geometries, then choose a named reference preview and optional primers. This checks sequence geometry, not laboratory success."</p>{move||pcr_review.get().map(|review|view!{<crate::pcr_outputs::PcrOutputsPanel source_molecule=molecule on_overlay=publish_pcr_overlay on_navigate=navigate_pcr_overlay review=review on_create=create_pcr_outputs guard=pcr_review_current busy=pcr_output_busy/>})}</div></Show>
                                            <Show when=move||panel.get()=="plan"><div class="inspector-content"><h2>"Cloning plan"</h2><div class="plan-exchange"><button class="button secondary full" disabled=move||workspace_restoring.get()||molecule.with(|m|m.sequence.alphabet()!=Alphabet::Dna) on:click=move|_|{modal_error.set(String::new());modal.set("plan_new");}>"New manual plan"</button><button class="button secondary full" on:click=export_workspace>"Export for AI"</button><button class="button secondary full" on:click=move|_|if let Some(input)=file_input.get(){input.click();}>"Open candidate"</button></div><p class="helper">"Select a separate intended target, then create a manual plan or open an AI candidate. Review the steps and verify before creating products."</p><Show when=move||candidate.get().is_some()><div class="candidate-review"><h3>"Candidate · review before use"</h3><button class="button secondary full" disabled=move||workspace_restoring.get() on:click=move|_|{modal_error.set(String::new());modal.set("plan_edit");}>"Edit steps and sources"</button>{move||candidate.get().map(|c|view!{<p><strong>{c.plan.name.clone()}</strong></p><p class="helper">{format!("Target revision {}…",&c.plan.target_revision[..12])}</p><h4>"Bound inputs"</h4><ul>{c.plan.inputs.iter().map(|i|view!{<li>{format!("{}: {} · {} bases · {}…",i.id,i.molecule.name,i.molecule.sequence.len(),&i.revision[..12])}</li>}).collect_view()}</ul><h4>"Ordered steps"</h4><ol>{c.plan.steps.iter().map(|step|{let description=match step{Operation::ResectThreePrime{cuts,..}=>format!("Shorten actual 3′ suffixes: left {} bp, right {} bp",cuts.left,cuts.right),Operation::HomologyAnnealing{fragments,..}=>format!("Anneal {} fragments; retain unfilled gaps and unsealed nicks",fragments.len()),Operation::HomologyAssembly{fragments,..}=>format!("Reconstruct homology-assembly strands from {} fragments",fragments.len()),Operation::OverlapExtension{options,..}=>format!("Reciprocal overlap extension across {} bp of terminal homology",options.overlap),Operation::ReconstructRunProduct{..}=>"Reconstruct saved material from its original run".into(),Operation::ReconstructMaterial{..}=>"Reconstruct saved material from its exact native derivation".into(),Operation::CatalogDigest{selection,..}=>format!("Catalog digest with {}",selection.enzymes.iter().map(|e|e.enzyme.name.as_str()).collect::<Vec<_>>().join(", ")),Operation::DeclareSourceTermini{..}=>"Declare original duplex termini".into(),Operation::TopoActions{..}=>"Accounted flap processing and physical duplex reconstruction".into(),Operation::ProcessedTopo{..}=>"TOPO treatment and physical duplex reconstruction".into(),Operation::CheckedTopo{..}=>"TOPO graph with explicit target acceptance".into(),Operation::Topo{..}=>"TOPO physical strand reconstruction".into(),Operation::MultiSiteRecombineReference{request,..}=>format!("Multi-parent {:?} reference exchange",request.reaction.reaction),Operation::RecombineReference{request,..}=>format!("Reference-only {:?} exchange",request.reaction),Operation::DeclareSourceModifications{..}=>"Declare original source modifications".into(),Operation::CheckedDigest{selection,..}=>format!("Checked digest with {}",selection.enzymes.iter().map(|e|e.enzyme.name.as_str()).collect::<Vec<_>>().join(", ")),Operation::Tail3Prime{base,..}=>format!("Add {} to both 3′ ends",base.as_str()),Operation::ProfiledOverlap{options,..}=>format!("Profiled {:?} overlap assembly",options.method),Operation::RepairEnds{method,..}=>format!("End repair: {:?}",method),Operation::AnnealOligos{candidate_index,..}=>format!("Anneal oligos: zero-based candidate {}",candidate_index),Operation::PcrWithMismatches{..}=>"Intentional-mismatch PCR".into(),Operation::SelectedGoldenGate{options,candidate_index,..}=>format!("Golden Gate with {:?} · explicit candidate {} and exact keys",options.enzyme,candidate_index+1),Operation::GoldenGate{options,..}=>format!("Golden Gate with {:?}",options.enzyme),Operation::Digest{enzymes,..}=>format!("Digest with {}",enzymes.join(", ")),Operation::TreatEnds{treatment,..}=>format!("Terminal {treatment:?}"),Operation::Ligate{fragments,circular,..}=>format!("Ligation of {} fragments → {}",fragments.len(),if *circular {"circular"}else{"linear"}),Operation::Pcr{forward,reverse,..}=>format!("PCR with {} and {}",forward.name,reverse.name),Operation::Overlap{fragments,circular,..}=>format!("Overlap assembly of {} fragments → {}",fragments.len(),if *circular {"circular"} else {"linear"})};let details=describe_step(step);view!{<li><details><summary>{description}</summary><pre>{details}</pre></details></li>}}).collect_view()}</ol>})}<button class="button primary full" on:click=verify_candidate>"Verify candidate"</button><button class="button secondary full" on:click=move|_|{candidate_events.set(None);candidate.set(None);plan_report.set(None);last_plan.set(None);last_plan_target.set(None);plan_generation.update(|g|*g=g.wrapping_add(1));}>"Discard candidate"</button></div></Show><Show when=move||candidate.get().is_none()><h3>"Manual cloning plan"</h3><p class="helper">"The active document is the intended target. Choose source documents and a cloning method; the core reconstructs the ordered candidate."</p><label>"First fragment"<select aria-label="First source document" prop:value=move||first_input.get().to_string() on:change=move|e|first_input.set(event_target_value(&e).parse().unwrap_or(0))>{move||documents.get().iter().enumerate().map(|(i,d)|view!{<option value=i.to_string()>{d.molecule.name.clone()}</option>}).collect_view()}</select></label><Show when=move||plan_method.get()!="restriction"||(restriction_final.get()=="ligate"&&second_restriction.get())><label>"Second fragment"<select aria-label="Second source document" prop:value=move||second_input.get().to_string() on:change=move|e|second_input.set(event_target_value(&e).parse().unwrap_or(0))>{move||documents.get().iter().enumerate().map(|(i,d)|view!{<option value=i.to_string()>{d.molecule.name.clone()}</option>}).collect_view()}</select></label></Show><label>"Method"<select aria-label="Cloning method" prop:value=move||plan_method.get() on:change=move|e|plan_method.set(if event_target_value(&e)=="restriction"{"restriction"}else{"overlap"})><option value="overlap">"Overlap assembly"</option><option value="restriction">"Restriction and ligation"</option></select></label><Show when=move||plan_method.get()=="restriction"><div class="checked-restriction-options"><label>"Restriction model"<select aria-label="Restriction model" prop:value=move||restriction_policy.get() on:change=move|e|restriction_policy.set(if event_target_value(&e)=="checked"{"checked"}else{"legacy"})><option value="legacy">"Assume complete cleavage (legacy)"</option><option value="checked">"Check inherited modification state"</option></select></label><label>"Final restriction product"<select aria-label="Restriction final product" prop:value=move||restriction_final.get() on:change=move|e|restriction_final.set(if event_target_value(&e)=="fragment"{"fragment"}else{"ligate"})><option value="ligate">"Ligate selected substrates"</option><option value="fragment">"Single processed fragment (no ligation)"</option></select></label><p class="helper">"Single-fragment mode keeps the selected digest product, or its explicitly repaired/treated descendant, as the final product. A blocked circle can remain an uncut circle without claiming a ligation."</p><Show when=move||restriction_final.get()=="ligate"><label class="checkbox-label"><input type="checkbox" aria-label="Use second restriction substrate" prop:checked=move||second_restriction.get() on:change=move|e|second_restriction.set(event_target_checked(&e))/>"Use a second restriction substrate"</label></Show><Show when=move||restriction_policy.get()=="checked"><p class="helper">"Declare the original source generation under a closed Dam/Dcm/CpG model. Original eligible positions inherit the declared state; positions outside these original motifs and newly synthesized bases are unmodified under this scope. Descendant motifs are not rescanned. Other modification systems and experimental conditions are not modeled. Unknown or impaired checked cleavage withholds undetermined products; blocked sites stay uncut."</p><details class="checked-profile-identities"><summary>"Pinned positional profiles and formulation identities"</summary><p>{format!("{} · revision {}",modified_duplex::PROFILE_VERSION,modified_duplex::profile_revision())}</p>{methylation::profiles().into_iter().map(|p|view!{<p>{format!("{}: {}",p.enzyme,p.supplier_identity)}</p>}).collect_view()}<p class="helper">"This positional profile is distinct from the source-context methylation review. Generic and HF/v2 enzyme formulations are not substituted. Imported summaries cannot authorize checked state."</p></details></Show>
                <RestrictionFields side="First" enzymes=first_enzymes fragment=first_fragment treatment=first_treatment repair=first_repair reverse=first_reverse/><Show when=move||restriction_policy.get()=="checked"><p class="checked-source-binding">{move||documents.with(|docs|docs.get(first_input.get()).map(|d|format!("First source: document {} · event {}",d.id,d.head_event())).unwrap_or_else(||"Choose a first source".into()))}</p><crate::checked_cloning::SourceAssumptions side="First" signals=first_modifications/></Show>
                <Show when=move||second_restriction.get()&&restriction_final.get()=="ligate"><RestrictionFields side="Second" enzymes=second_enzymes fragment=second_fragment treatment=second_treatment repair=second_repair reverse=second_reverse/><Show when=move||restriction_policy.get()=="checked"><p class="checked-source-binding">{move||documents.with(|docs|docs.get(second_input.get()).map(|d|format!("Second source: document {} · event {}",d.id,d.head_event())).unwrap_or_else(||"Choose a second source".into()))}</p><crate::checked_cloning::SourceAssumptions side="Second" signals=second_modifications/></Show></Show><Show when=move||restriction_final.get()=="ligate"><label class="checkbox-label"><input type="checkbox" aria-label="Allow restriction ligation nicks" prop:checked=move||allow_nicked.get() on:change=move|e|allow_nicked.set(event_target::<HtmlInputElement>(&e).checked())/>"Allow explicitly reported nicks"</label></Show><p class="helper">"Fragment numbers follow ascending reference-strand cut order, not size. Intact circles cannot ligate. The active document remains the independent target; declared source state is not measured modification evidence."</p></div></Show><Show when=move||plan_method.get()=="overlap"><label>"Junction overlap (bp)"<input type="number" min="1" prop:value=move||overlap.get() on:input=move|e|overlap.set(event_target_value(&e))/></label><Show when=move||molecule.get().topology==Topology::Circular><label>"Circular closure overlap (bp)"<input type="number" min="1" prop:value=move||closure_overlap.get() on:input=move|e|closure_overlap.set(event_target_value(&e))/></label></Show></Show><button class="button primary full" on:click=do_plan>"Reconstruct and verify"</button></Show>{move||plan_report.with(|report|report.as_ref().map(|report|view!{<div class="verification-result"><strong>{format!("{:?}",report.status)}</strong>{report.findings.iter().cloned().map(|f|view!{<p>{f}</p>}).collect_view()}{report.physical_processing_reports.iter().map(|step|{match &step.evidence{collagene_core::verification::PhysicalProcessingEvidence::ThreePrimeResection(evidence)=>crate::assembly_product_view::three_prime_processing_view(step.step_index,evidence),collagene_core::verification::PhysicalProcessingEvidence::HomologyAnnealing(evidence)=>crate::assembly_product_view::annealing_processing_view(step.step_index,evidence),collagene_core::verification::PhysicalProcessingEvidence::HomologyAssembly(evidence)=>crate::assembly_product_view::homology_processing_view(step.step_index,evidence),      collagene_core::verification::PhysicalProcessingEvidence::TopoActions(evidence)=>{let graph=&evidence.projection.processing.product;view!{<details><summary>{format!("Step {} · accounted chemical actions",step.step_index+1)}</summary><p>{format!("All {} source calls retained · {} components · {} remaining nicks · {} unpaired calls",evidence.projection.processing.source_bases,graph.components.len(),graph.nicks.len(),graph.unpaired.iter().map(|r|r.end-r.start).sum::<usize>())}</p><p class="helper">"Ideal hydrolysis is a declared model, not a named nuclease or host-repair guarantee. Complete original, main and excised component graphs remain in run JSON evidence."</p>{evidence.bridge_notes.iter().map(|s|view!{<p>{s.clone()}</p>}).collect_view()}</details>}.into_any()},     collagene_core::verification::PhysicalProcessingEvidence::OverlapExtension(geometry)=>{view!{<details><summary>{format!("Step {} · reciprocal extension",step.step_index+1)}</summary><p>{format!("Homology {} bp · terminal pairings {:?} · complete {} · unique within terminal model {} · {} comparisons",geometry.options.overlap,geometry.terminal_pairings.overlaps,geometry.terminal_pairings.complete,geometry.terminal_pairings.unique,geometry.terminal_pairings.comparisons)}</p></details>}.into_any()}     collagene_core::verification::PhysicalProcessingEvidence::LigationOutcomes(evidence)=>crate::ligation_outcomes_view::ligation_outcomes_view(step.step_index,evidence), collagene_core::verification::PhysicalProcessingEvidence::DerivedOrigin(evidence)=>{view!{<details><summary>"Reusable material derivation replay"</summary><p>{format!("Exact event {} · inherited incomplete {} · {} retained findings",evidence.origin.product.event_id,evidence.material.origin_incomplete,evidence.material.origin_findings_total)}</p><p class="helper">"The retained target-free PCR/digest recipe was replayed; imported reference or end summaries were not stock."</p></details>}.into_any()}     collagene_core::verification::PhysicalProcessingEvidence::RunProduct(evidence)=>{let origin=evidence.decode_report().ok();view!{<details><summary>{format!("Step {} · saved material reconstruction",step.step_index+1)}</summary><p>{format!("Origin run {} · product {}",evidence.origin.selection.run_id,evidence.origin.selection.step_index+1)}</p><p>{origin.as_ref().map(|r|format!("Origin assessment: {:?} · {}",r.status,r.engine_version))}</p><p class="helper">"Already-prepared material is reconstructed against its original target. Nicks, termini, positional modifications and incomplete findings remain. The complete origin report is retained in reconstruction JSON."</p>{origin.into_iter().flat_map(|r|r.findings.into_iter().take(100)).map(|finding|view!{<p>{finding}</p>}).collect_view()}</details>}.into_any()}     }}).collect_view()}{report.processed_step_reports.iter().map(|step|{let sealing=&step.evidence.projection.sealing;view!{<details class="processed-topo-evidence"><summary>{format!("Step {} · TOPO treatment and duplex evidence",step.step_index+1)}</summary><p>{format!("Original nicks {} · remaining nicks {} · components {} · unpaired runs {} · faithful material available {}",sealing.parent.product.nicks.len(),sealing.product.nicks.len(),sealing.product.components.len(),sealing.product.unpaired.len(),step.evidence.tracked.is_some())}</p>{sealing.steps.iter().map(|seal|view!{<p>{format!("Nick {}: {:?} · {} · 3′ {:?} / 5′ {:?}",seal.step+1,seal.outcome,seal.reason,seal.three_prime,seal.five_prime)}</p>}).collect_view()}{step.evidence.bridge_notes.iter().map(|note|view!{<p>{note.clone()}</p>}).collect_view()}<p class="helper">"Full original and treated graphs remain in reconstruction JSON and the retained run evidence. These results do not establish experimental isolation or success."</p></details>}}).collect_view()}{report.annotation_reports.iter().map(|step|view!{<details class="annotation-evidence"><summary>{format!("Step {} · source annotation correspondence",step.step_index+1)}</summary><p class="helper">"Mapped coordinates follow the actual selected top-reference calls. Review PCR-derived and disrupted annotations. Original unmapped features remain evidence in the reconstruction JSON, not invented product coordinates."</p>{step.evidence.sources.iter().map(|source|view!{<p>{format!("{} · {} projected · {} outside contour · {} unmapped originals · projection supported {}",source.role,source.projected,source.dropped,source.unmapped.len(),source.projection_supported)}</p>}).collect_view()}</details>}).collect_view()}{report.step_reports.iter().filter_map(|step|step.topo.as_ref().map(|topo|{let outcome=match &topo.projection.outcome{collagene_core::topo::physical::projection::Outcome::Proven(contour)=>format!("Paired contour proven; {} nicks retained. No repair is inferred.",contour.nicks.len()),collagene_core::topo::physical::projection::Outcome::Withheld(_)=>"No unique paired contour; physical graph retained and reference withheld.".into()};view!{<details class="physical-evidence"><summary>{format!("Step {} · TOPO physical graph",step.step_index+1)}</summary><p>{outcome}</p><p>{format!("{} physical components · batch {}",topo.projection.physical.product.components.len(),topo.batch_revision)}</p><p class="helper">"Download the reconstruction report for the full physical graph, original source maps and projection evidence. Saved evidence needs fresh replay before use."</p></details>}})).collect_view()}{report.products.iter().enumerate().filter_map(|(i,p)|p.metadata.get("duplex_model").and_then(|s|serde_json::from_str::<chemistry::DuplexSummary>(s).ok()).map(|summary|view!{<details class="physical-evidence"><summary>{format!("Step {} · {}",i+1,p.name)}</summary><pre>{summary.describe()}</pre></details>})).collect_view()}{report.products.iter().enumerate().filter_map(|(i,p)|crate::checked_cloning::modification_summary(p).map(|summary|view!{<details class="checked-modification-summary"><summary>{format!("Step {} · positional modification evidence",i+1)}</summary><p>{summary}</p></details>})).collect_view()}<details><summary>"Assumptions and limits"</summary><ul>{report.assumptions.iter().cloned().map(|a|view!{<li>{a}</li>}).collect_view()}</ul></details><p class="helper">"Alternative assemblies and experimental feasibility are not exhaustively checked."</p></div>}))}<Show when=move||reviewed_product_current.get()>{move||last_plan.get().zip(plan_report.get()).map(|((plan,target),report)|view!{<crate::assembly_product_view::AssemblyProductView plan target report epoch=workspace_epoch guard=Some(reviewed_product_guard) on_products=Callback::new(move|products:Vec<Molecule>|{reviewed_product_guard.run(())?;add_documents(documents,active,workspace_restoring,products)}) on_finalize=Some(finalize_reviewed_candidate)/>})}</Show><button class="button secondary full" disabled=move||plan_report.with(Option::is_none) on:click=move|_|{let result=plan_report.with_untracked(|report|report.as_ref().map(crate::checked_cloning::prepare_report_export));if let Some(result)=result{crate::app::complete_download(result,move|result|{if message.try_get_untracked().is_some(){if let Err(error)=result{message.set(error);}}});}}>"Download reconstruction report"</button><button class="button secondary full" disabled=move||candidate.get().is_none()&&last_plan.get().is_none() on:click=move|_|{let value=candidate.get_untracked().or_else(||last_plan.get_untracked().map(|(plan,_)|transfer::Candidate::new(plan)));if let Some(value)=value {match transfer::save_candidate(&value){Ok(bytes)=>{crate::app::download_with_completion("candidate.cg2plan","application/octet-stream",&bytes,move|result|{if message.try_get_untracked().is_none(){return;}if let Err(e)=result {message.set(e);}})},Err(e)=>message.set(e.to_string())}}}>"Save candidate"</button><button class="button secondary full" disabled=move||last_plan.get().is_none() on:click=move|_|{modal_error.set(String::new());modal.set("record_run");}>"Keep reconstructed products"</button></div></Show>
                                            <Show when=move||panel.get()=="primers"><crate::primers::PrimerPanel on_select_binding=navigate_primer_binding resolve_source=resolve_retained_source on_collection=Callback::new(move|table:oligo_formats::OligoTable|import_oligos.run((table,Vec::new()))) molecule=molecule selection=selection seed=primer_seed source_binding=primer_context pair_draft pair_selection anchor_available=primer_anchor_available primer_sources list_primers=list_saved_primers resolve_primer=resolve_saved_primer epoch=workspace_epoch on_use=Callback::new(move|(forward,reverse):(Primer,Primer)|{pcr_forward_name.set(forward.name);pcr_reverse_name.set(reverse.name);pcr_forward.set(forward.annealing.as_str().into());pcr_forward_tail.set(forward.tail);pcr_reverse.set(reverse.annealing.as_str().into());pcr_reverse_tail.set(reverse.tail);panel.set("pcr");message.set("Candidate primers loaded. Compute PCR products to inspect alternative binding geometries.".into());}) on_annotate=Callback::new(move|(binding,feature):(collagene_core::run_lineage::EventBinding,Feature)|{if primer_context.get_untracked()!=binding{return Err("The document event changed; rerun primer analysis before annotating it.".into());}let current=molecule.get_untracked();let edited=editing::set_feature(&current,None,feature).map_err(|e|e.to_string())?;commit_document(documents,active,workspace_restoring,&binding.revision,edited)})/></Show>
                                            <div style:display=move||{if panel.get()=="gel"{"contents"}else{"none"}}><crate::gel::GelPanel seeds=gel_digest_seeds epoch=workspace_epoch resolve_event=resolve_retained_source molecule=molecule query_options=Memo::new(move|_|documents.with(|docs|docs.iter().enumerate().map(|(i,d)|(i,d.molecule.name.clone())).collect())) get_query=Callback::new(move|i:usize|documents.with(|docs|docs.get(i).map(|d|d.molecule.clone())))/></div>
                                            <Show when=move||panel.get()=="methods"><crate::cloning_methods::CloningMethodsPanel molecule=molecule get_workspace=Some(contig_snapshot) selected_region=Callback::new(move|():()|{if inspected_feature.get_untracked().is_some_and(|f|f.location.spans().len()>1){return None;}selection.get_untracked().filter(|(a,b)|a<b).map(|(a,b)|(active.get_untracked(),a,b))}) target_info=Memo::new(move|_|crate::cloning_methods::TargetContext::Document(documents.with(|docs|recombination_source_info(&docs[active.get()])))) source_info=Callback::new(move|i:usize|documents.with(|docs|docs.get(i).map(recombination_source_info))) primer_sources=primer_sources list_primers=list_saved_primers resolve_primer=resolve_saved_primer query_options=Memo::new(move|_|documents.with(|docs|docs.iter().enumerate().map(|(i,d)|(i,d.molecule.name.clone())).collect())) get_query=Callback::new(move|i:usize|documents.with(|docs|docs.get(i).map(|d|d.molecule.clone()))) config=method_config epoch=workspace_epoch on_candidate=Callback::new(move|value:transfer::Candidate|{
                                            if workspace_restoring.get_untracked(){return;}
                                            let state=method_config.get_untracked();

                                            let captured=documents.with_untracked(|docs|->Result<_,String>{
                                                let indices=state.source_indices()?;
                                                let target=docs.get(active.get_untracked()).map(current_event_binding).ok_or("No independent target is active")?;
                                                if target.revision!=value.plan.target_revision||indices.len()!=value.plan.inputs.len(){return Err("Cloning source/target choices changed; reconstruct the method evidence again.".into());}
                                                let bindings=value.plan.inputs.iter().zip(indices).map(|(input,index)|{
                                                    let doc=docs.get(index).ok_or("A selected source document is unavailable")?;
                                                    let event=current_event_binding(doc);
                                                    if event.revision!=input.revision{return Err("A source event changed; reconstruct the method evidence again.".to_string());}
                                                    Ok(collagene_core::run_lineage::InputEventBinding{input_id:input.id.clone(),event})
                                                }).collect::<Result<Vec<_>,String>>()?;
                                                Ok((target,bindings))
                                            });
                                            let (target,bindings)=match captured{Ok(value)=>value,Err(error)=>{message.set(error);return;}};
                                            let preferred=(value.plan_revision.clone(),workspace_epoch.get_untracked(),target,bindings);
                                            candidate_events.set(None);candidate.set(Some(value));candidate_events.set(Some(preferred));plan_report.set(None);last_plan.set(None);last_plan_target.set(None);last_plan_events.set(None);plan_generation.update(|g|*g=g.wrapping_add(1));panel.set("plan");message.set("Method candidate opened with exact source-event bindings for independent target verification. No experimental success or mixture uniqueness is established.".into());
                                        }) on_products=Callback::new(move|products:Vec<Molecule>|add_documents(documents,active,workspace_restoring,products))/></Show>
                                            <Show when=move||panel.get()=="ends"><crate::end_repair::EndRepairPanel epoch=workspace_epoch get_workspace=library_snapshot on_material=open_material_derivation source_info=Callback::new(move|i:usize|documents.with(|docs|docs.get(i).map(recombination_source_info))) molecule=molecule query_options=Memo::new(move|_|documents.with(|docs|docs.iter().enumerate().map(|(i,doc)|(i,doc.molecule.name.clone())).collect())) get_query=Callback::new(move|index:usize|documents.with(|docs|docs.get(index).map(|doc|doc.molecule.clone()))) on_candidate=Callback::new(move|value:transfer::Candidate|{if workspace_restoring.get_untracked(){return;}candidate_events.set(None);candidate.set(Some(value));plan_report.set(None);last_plan.set(None);last_plan_target.set(None);plan_generation.update(|generation|*generation=generation.wrapping_add(1));panel.set("plan");message.set("Physical candidate opened for review. Verify it against the independently selected active target; alternatives and chemistry remain explicit.".into());}) on_products=Callback::new(move|products:Vec<Molecule>|add_documents(documents,active,workspace_restoring,products))/></Show>
                                            <Show when=move||panel.get()=="domestication"><crate::domestication_view::DomesticationView source=Memo::new(move|_|documents.with(|docs|recombination_source_info(&docs[active.get()]))) molecule selected_feature epoch=workspace_epoch on_open=open_domestication on_navigate=Callback::new(move|(start,end):(usize,usize)|->Result<(),String>{
                                                let n=molecule.with_untracked(|m|m.sequence.len());
                                                if workspace_restoring.get_untracked()||start>=n||end<=start||end-start>n{return Err("The selected source range is no longer available.".into());}
                                                let spans=if end>n{vec![(start,n),(0,end-n)]}else{vec![(start,end)]};
                                                batch(||{selected_feature.set(None);selection.set(if spans.len()==1{spans.first().copied()}else{None});search_regions.set(Some((molecule.with_untracked(Molecule::revision),spans)));});
                                                Ok(())
                                            })/></Show>
                                            <Show when=move||panel.get()=="mutagenesis"><crate::mutagenesis::MutagenesisPanel get_workspace=library_snapshot on_material=open_material_derivation resolve_source=resolve_retained_source on_overlay=Some(publish_pcr_overlay) on_navigate=Some(navigate_pcr_overlay) molecule=molecule mutation_selection mutation_seed on_intent=open_mutation_intent pair_draft pair_selection anchor_available=primer_anchor_available primer_sources list_primers=list_saved_primers resolve_primer=resolve_saved_primer epoch=workspace_epoch source_info=Callback::new(move|i:usize|documents.with(|docs|docs.get(i).map(recombination_source_info))) target_binding=primer_context query_options=Memo::new(move|_|documents.with(|docs|docs.iter().enumerate().map(|(i,d)|(i,d.molecule.name.clone())).collect())) get_query=Callback::new(move|i:usize|documents.with(|docs|docs.get(i).map(|d|d.molecule.clone()))) on_candidate=Callback::new(move|opened:crate::mutagenesis::CandidateOpen|{
                                                let result=(||->Result<_,String>{
                                                    if workspace_restoring.get_untracked() || opened.epoch!=workspace_epoch.get_untracked(){return Err("The mutation workspace changed; review the current design again.".into());}
                                                    let frame=documents.with_untracked(|docs|workspace_snapshot(docs,active.get_untracked(),candidate.get_untracked().as_ref(),&stored_candidates.get_untracked(),&attachments.get_untracked(),&workspace_name.get_untracked()))?;
                                                    crate::assembly_review::stage_existing(&frame,&opened.candidate,&opened.target,&opened.sources)?;
                                                    Ok(())
                                                })();
                                                if let Err(error)=result{message.set(error);return;}
                                                let preferred=(opened.candidate.plan_revision.clone(),opened.epoch,opened.target,opened.sources);
                                                batch(||{candidate_events.set(None);candidate.set(Some(opened.candidate));candidate_events.set(Some(preferred));plan_report.set(None);last_plan.set(None);last_plan_target.set(None);last_plan_events.set(None);plan_generation.update(|g|*g=g.wrapping_add(1));panel.set("plan");});
                                                message.set("Inverse-PCR candidate opened with its original source and intended-target events. Verify, inspect the reconstructed product, then choose outputs.".into());
                                            }) on_outputs=open_mismatch_outputs/></Show>

                                            <Show when=move||panel.get()=="topo"><crate::topo_view::TopoPanel
                                                target=Memo::new(move|_|documents.with(|docs|recombination_source_info(&docs[active.get()])))
                                                sources=Memo::new(move|_|documents.with(|docs|docs.iter().map(recombination_source_info).collect()))
                                                capture=Callback::new(move|capture:crate::topo_view::CaptureRequest|->Result<crate::topo_view::Captured,String>{
                                                    if workspace_restoring.get_untracked(){return Err("Wait for the workspace operation to finish.".into());}
                                                    if capture.sources.is_empty()||capture.sources.len()>32{return Err("Choose 1–32 explicit source roles.".into());}
                                                    let template_origin=capture.template_origin.as_ref().map(|selection|{
                                                        let (_,binding)=capture.sources.iter().find(|(role,_)|role=="insert_template").ok_or("The selected PCR template role is unavailable.")?;
                                                        if selection.event()!=binding{return Err("The saved template origin belongs to another source event.".to_string());}
                                                        let frame=library_snapshot.run(())?;
                                                        let inspection=collagene_core::run_material::inspect_material(&frame,binding,&Some(selection.clone())).map_err(|error|error.to_string())?;
                                                        let origin=inspection.selected.ok_or("The exact saved material origin is unavailable.")?;
                                                        origin.validate_for_workspace(&frame,binding).map_err(|error|error.to_string())?;
                                                        Ok::<_,String>(origin)
                                                    }).transpose()?;
                                                    documents.with_untracked(|docs|{
                                                        let resolve=|binding:&collagene_core::run_lineage::EventBinding|docs.iter().find(|d|current_event_binding(d)==*binding).ok_or_else(||"A selected source event changed. Inspect again.".to_string());
                                                        let target=resolve(&capture.target)?;
                                                        if docs.get(active.get_untracked()).map(current_event_binding).as_ref()!=Some(&capture.target){return Err("The independent target changed.".into());}
                                                        if target.molecule.sequence.alphabet()!=Alphabet::Dna||target.molecule.sequence.len()>100_000{return Err("Choose a DNA target within the 100,000-base TOPO profile.".into());}
                                                        let mut molecules=Vec::new();let mut roles=Vec::new();
                                                        for (role,binding) in &capture.sources{
                                                            if binding.document_id==capture.target.document_id{return Err("Choose an independent target document, separate from the sources.".into());}
                                                            let source=resolve(binding)?;molecules.push(&source.molecule);roles.push(verification::TopoSource{role:role.clone(),source:Source::Input(role.clone())});
                                                        }
                                                        let mut models=molecules.clone();models.push(&target.molecule);
                                                        collagene_core::model::bounds::model_preflight(&models).map_err(|e|e.to_string())?;
                                                        let inputs=collagene_core::topo_run_preparation::source_inputs(&roles,&molecules,&capture.request,"insert_template",template_origin.as_ref()).map_err(|e|e.to_string())?;
                                                        target.molecule.validate().map_err(|e|e.to_string())?;
                                                        Ok(crate::topo_view::Captured{inputs,target:target.molecule.clone(),template_origin})
                                                    })
                                                }) get_workspace=library_snapshot epoch=workspace_epoch locked=Memo::new(move|_|workspace_restoring.get())
                                                on_candidate=Callback::new(move|(value,bindings):(transfer::Candidate,Vec<collagene_core::run_lineage::InputEventBinding>)|{
                                                    if workspace_restoring.get_untracked(){return;}
                                                    let target=documents.with_untracked(|docs|docs.get(active.get_untracked()).map(current_event_binding));let Some(target)=target else{return;};
                                                    if value.plan.target_revision!=target.revision||bindings.len()!=value.plan.inputs.len()||!documents.with_untracked(|docs|bindings.iter().all(|binding|value.plan.inputs.iter().any(|input|input.id==binding.input_id&&input.revision==binding.event.revision)&&docs.iter().any(|doc|current_event_binding(doc)==binding.event))){message.set("TOPO source/target event bindings changed. Inspect again.".into());return;}
                                                    let preferred=(value.plan_revision.clone(),workspace_epoch.get_untracked(),target,bindings);
                                                    candidate_events.set(None);candidate.set(Some(value));candidate_events.set(Some(preferred));plan_report.set(None);last_plan.set(None);last_plan_target.set(None);last_plan_events.set(None);plan_generation.update(|g|*g=g.wrapping_add(1));panel.set("plan");message.set("TOPO candidate opened for source-bound review and independent verification. Full physical evidence remains distinct from reference DNA.".into());
                                                })/></Show>
                                            <Show when=move||panel.get()=="recombination"><details><summary>"Multi-parent reference exchange"</summary><crate::multisite::MultiSitePanel
                                                target=Memo::new(move|_|documents.with(|docs|recombination_source_info(&docs[active.get()])))
                                                sources=Memo::new(move|_|documents.with(|docs|docs.iter().map(recombination_source_info).collect()))
                                                capture=Callback::new(move|capture:crate::multisite::CaptureRequest| {
                                                    if workspace_restoring.get_untracked(){return Err("Wait for the workspace operation to finish.".into());}
                                                    documents.with_untracked(|docs| {
                                                        let resolve=|binding:&collagene_core::run_lineage::EventBinding|docs.iter().find(|d|d.id==binding.document_id&&current_event_binding(d)==*binding).ok_or_else(||"The selected document event changed. Inspect again.".to_string());
                                                        if !(2..=5).contains(&capture.sources.len()){return Err("Choose two to five original parents.".into());}
                                                        let parents=capture.sources.iter().map(resolve).collect::<Result<Vec<_>,_>>()?;
                                                        let target=resolve(&capture.target)?;
                                                        if docs.get(active.get_untracked()).map(current_event_binding).as_ref()!=Some(&capture.target){return Err("The active independent target changed.".into());}
                                                        let molecules=parents.iter().map(|d|&d.molecule).collect::<Vec<_>>();
                                                        collagene_core::multisite_recombination::source_preflight(&molecules,&capture.request).map_err(|e|e.to_string())?;
                                                        if target.molecule.sequence.len()>2_000_000{return Err("Independent target exceeds 2,000,000 bases.".into());}
                                                        collagene_core::model::bounds::model_preflight(&[&target.molecule]).map_err(|e|e.to_string())?;
                                                        target.molecule.validate().map_err(|e|e.to_string())?;
                                                        let inputs=parents.iter().enumerate().map(|(i,d)|collagene_core::recombination::reference::Input{role:format!("parent_{}",i+1),molecule:d.molecule.clone()}).collect();
                                                        Ok(crate::multisite::Captured{inputs,target:target.molecule.clone()})
                                                    })
                                                }) epoch=workspace_epoch locked=Memo::new(move|_|workspace_restoring.get())
                                                on_candidate=Callback::new(move|(value,bindings):(transfer::Candidate,Vec<collagene_core::run_lineage::InputEventBinding>)| {
                                                    if workspace_restoring.get_untracked(){return;}
                                                    let Some(target)=documents.with_untracked(|docs|docs.get(active.get_untracked()).map(current_event_binding)) else{return;};
                                                    if value.plan.target_revision!=target.revision||bindings.len()!=value.plan.inputs.len()||!documents.with_untracked(|docs|bindings.iter().all(|binding|value.plan.inputs.iter().any(|input|input.id==binding.input_id&&input.revision==binding.event.revision)&&docs.iter().any(|doc|current_event_binding(doc)==binding.event))){message.set("Multi-parent source or target event bindings changed. Inspect again.".into());return;}
                                                    let preferred=(value.plan_revision.clone(),workspace_epoch.get_untracked(),target,bindings);
                                                    candidate_events.set(None);candidate.set(Some(value));candidate_events.set(Some(preferred));plan_report.set(None);last_plan.set(None);last_plan_target.set(None);last_plan_events.set(None);plan_generation.update(|g|*g=g.wrapping_add(1));panel.set("plan");message.set("Multi-parent reference candidate opened. Verify against the independent target; reference ancestry remains chemically restricted.".into());
                                                })/></details></Show>
                                            <Show when=move||panel.get()=="recombination"><crate::recombination::RecombinationPanel
                                                target=Memo::new(move|_|documents.with(|docs|recombination_source_info(&docs[active.get()])))
                                                sources=Memo::new(move|_|documents.with(|docs|docs.iter().map(recombination_source_info).collect()))
                                                capture=Callback::new(move|capture:crate::recombination::CaptureRequest| {
                                                    if workspace_restoring.get_untracked(){return Err("Wait for the workspace operation to finish.".into());}
                                                    documents.with_untracked(|docs| {
                                                        let resolve=|binding:&collagene_core::run_lineage::EventBinding|docs.iter().find(|d|d.id==binding.document_id && current_event_binding(d)==*binding).ok_or_else(||"The selected document event changed. Inspect again.".to_string());
                                                        let a=resolve(&capture.sources[0])?;
                                                        let b=resolve(&capture.sources[1])?;
                                                        let target=resolve(&capture.target)?;
                                                        if docs.get(active.get_untracked()).map(current_event_binding).as_ref()!=Some(&capture.target){return Err("The active independent target changed.".into());}
                                                        collagene_core::recombination::source_target_preflight([&a.molecule,&b.molecule],&capture.request,&target.molecule).map_err(|e|e.to_string())?;
                                                        let inputs=verification::recombination_inputs(&[Source::Input("source_1".into()),Source::Input("source_2".into())],[&a.molecule,&b.molecule],&capture.request).map_err(|e|e.to_string())?;
                                                        Ok(crate::recombination::Captured{inputs,target:target.molecule.clone()})
                                                    })
                                                }) epoch=workspace_epoch locked=Memo::new(move|_|workspace_restoring.get())
                                                on_candidate=Callback::new(move|(value,bindings):(transfer::Candidate,Vec<collagene_core::run_lineage::InputEventBinding>)| {
                                                    if workspace_restoring.get_untracked(){return;}
                                                    let target=documents.with_untracked(|docs|docs.get(active.get_untracked()).map(current_event_binding));
                                                    let Some(target)=target else{return;};
                                                    if value.plan.target_revision!=target.revision || bindings.len()!=2 || !documents.with_untracked(|docs|bindings.iter().all(|b|docs.iter().any(|d|current_event_binding(d)==b.event))){message.set("Reference-exchange event bindings changed. Inspect again.".into());return;}
                                                    let preferred=(value.plan_revision.clone(),workspace_epoch.get_untracked(),target,bindings);
                                                    candidate_events.set(None);candidate.set(Some(value));candidate_events.set(Some(preferred));plan_report.set(None);last_plan.set(None);last_plan_target.set(None);last_plan_events.set(None);plan_generation.update(|g|*g=g.wrapping_add(1));panel.set("plan");message.set("Reference-only candidate opened for review. Verify against the independent target before recording source-event ancestry.".into());
                                                })/></Show>
                                            <Show when=move||panel.get()=="methylation"><crate::methylation::MethylationPanel source=Memo::new(move|_|documents.with(|docs|{let d=&docs[active.get()];crate::methylation::SourceInfo{binding:crate::contigs::ReadBinding{document_id:d.id.clone(),event_id:d.head_event().to_string(),revision:d.molecule.revision()},name:d.molecule.name.clone(),length:d.molecule.sequence.len(),alphabet:d.molecule.sequence.alphabet()}})) get_molecule=Callback::new(move|()|molecule.get_untracked()) epoch=workspace_epoch on_navigate=navigate_read_call on_digest=Callback::new(move|()|panel.set("enzymes"))/></Show>
                                            <Show when=move||panel.get()=="genomes"><crate::genomes::GenomePanel links=genome_links epoch=workspace_epoch on_link=Callback::new(move|link:collagene_core::genome::ReferenceLink|->Result<(),String>{
                                                if workspace_restoring.get_untracked(){return Err("Wait for workspace restoration.".into());}
                                                let(hash,attachment)=collagene_core::genome::make_attachment(&link).map_err(|e|e.to_string())?;
                                                let mut assets=attachments.get_untracked();assets.retain(|_,a|a.kind!=collagene_core::genome::WORKSPACE_ATTACHMENT_KIND||collagene_core::genome::load_link(&a.bytes).is_ok_and(|old|old.role!=link.role));assets.insert(hash,attachment);
                                                let frame=documents.with_untracked(|docs|workspace_snapshot(docs,active.get_untracked(),candidate.get_untracked().as_ref(),&stored_candidates.get_untracked(),&assets,&workspace_name.get_untracked()))?;frame.validate().map_err(|e|e.to_string())?;attachments.set(assets);Ok(())
                                            }) on_unlink=Callback::new(move|role:String|->Result<(),String>{if workspace_restoring.get_untracked(){return Err("Wait for workspace restoration.".into());}attachments.update(|assets|assets.retain(|_,a|a.kind!=collagene_core::genome::WORKSPACE_ATTACHMENT_KIND||collagene_core::genome::load_link(&a.bytes).is_ok_and(|l|l.role!=role)));Ok(())})/></Show>
                                            <Show when=move||panel.get()=="crispr"><crate::crispr::CrisprPanel source=Memo::new(move|_|documents.with(|docs|{let d=&docs[active.get()];crate::crispr::DocumentInfo{binding:crate::contigs::ReadBinding{document_id:d.id.clone(),event_id:d.head_event().to_string(),revision:d.molecule.revision()},name:d.molecule.name.clone(),length:d.molecule.sequence.len(),alphabet:d.molecule.sequence.alphabet(),topology:d.molecule.topology}})) documents=Memo::new(move|_|documents.with(|docs|docs.iter().map(|d|crate::crispr::DocumentInfo{binding:crate::contigs::ReadBinding{document_id:d.id.clone(),event_id:d.head_event().to_string(),revision:d.molecule.revision()},name:d.molecule.name.clone(),length:d.molecule.sequence.len(),alphabet:d.molecule.sequence.alphabet(),topology:d.molecule.topology}).collect())) get_document=Callback::new(move|id:String|documents.with_untracked(|docs|docs.iter().find(|d|d.id==id).map(|d|d.molecule.clone()))) epoch=workspace_epoch on_navigate=navigate_read_call on_record=record_donor/></Show>
                                            <Show when=move||panel.get()=="contigs"><crate::contigs::ContigPanel read_options=Memo::new(move|_|documents.with(|docs|docs.iter().map(|d|crate::contigs::ReadInfo{document_id:d.id.clone(),event_id:d.head_event().to_string(),revision:d.molecule.revision(),name:d.molecule.name.clone(),length:d.molecule.sequence.len(),eligible:d.molecule.sequence.alphabet()==Alphabet::Dna&&d.molecule.topology==Topology::Linear,quality_label:if d.molecule.metadata.contains_key("trace_original_base64"){"Sanger trace · native caller evidence".into()}else if d.molecule.metadata.contains_key(collagene_core::oligo_formats::FASTQ_ORIGINAL_KEY){format!("FASTQ · {}",d.molecule.metadata.get("fastq:encoding").map(String::as_str).unwrap_or("declared encoding"))}else{"No quality evidence".into()}}).collect())) get_read=Callback::new(move|id:String|documents.with_untracked(|docs|docs.iter().find(|d|d.id==id).map(|d|d.molecule.clone()))) epoch=workspace_epoch on_record=record_assembly on_navigate=navigate_read_call/></Show>
                                            <Show when=move||panel.get()=="codons"><crate::codon_design::CodonDesignPanel molecule=molecule source_identity=Memo::new(move|_|documents.with(|docs|format!("{}:{}",docs[active.get()].id,docs[active.get()].head_event()))) epoch=workspace_epoch on_open=Callback::new(move|product:Molecule|{
                                                if workspace_restoring.get_untracked(){return Err("Wait for the workspace update to finish.".into());}
                                                let mut docs=documents.get_untracked();let first=docs.len();
                                                if first>=1000||docs.iter().map(|d|document_bytes(&d.molecule)).sum::<usize>()+document_bytes(&product)>MAX_WORKSPACE_BYTES{return Err("This derived document exceeds the workspace document or memory limit.".into());}
                                                let source=docs.get(active.get_untracked()).ok_or("No active source document")?;
                                                let prefix=format!("document-{}",&product.revision()[..16]);let mut suffix=1;
                                                let id=loop{let id=format!("{prefix}-{suffix}");if docs.iter().all(|d|d.id!=id){break id;}suffix+=1;};
                                                let derived=source.derive_event(source.head_event(),id,product,"Opened codon candidate after independent translation and constraint checks").map_err(|e|e.to_string())?;
                                                docs.push(derived);
                                                let frame=workspace_snapshot(&docs,first,candidate.get_untracked().as_ref(),&stored_candidates.get_untracked(),&attachments.get_untracked(),&workspace_name.get_untracked())?;
                                                frame.validate().map_err(|e|e.to_string())?;
                                                batch(||{documents.set(docs);active.set(first);});Ok(())
                                            })/></Show>
                                            <Show when=move||panel.get()=="biology"><Show when=move||molecule.with(|m|m.sequence.alphabet()!=Alphabet::Protein)><crate::biology::BiologyPanel source=Memo::new(move|_|documents.with(|docs|recombination_source_info(&docs[active.get()]))) epoch=workspace_epoch selection_location=biology_selection_location on_select_location=select_biology_location on_open_protein=open_translated_protein on_retain_cds=retain_translation_cds molecule=molecule selected_feature=selected_feature selection=selection on_annotate=Callback::new(move|(base,feature):(String,Feature)|{let current=molecule.get_untracked();if current.revision()!=base{return Err("The document changed; rerun ORF analysis before annotating it.".into());}let edited=editing::set_feature(&current,None,feature).map_err(|e|e.to_string())?;commit_document(documents,active,workspace_restoring,&base,edited)})/></Show><crate::biology_source_view::TranslationSourceLink molecule on_open=open_translation_origin/><crate::protein_properties::ProteinPanel molecule=molecule/></Show>
                                            <Show when=move||panel.get()=="dotplot"><crate::dotplot::DotPlotPanel molecule=molecule active_index=Memo::new(move|_|active.get()) source_identity=Memo::new(move|_|documents.with(|docs|format!("{}:{}",docs[active.get()].id,docs[active.get()].head_event()))) query_options=Memo::new(move|_|documents.with(|docs|docs.iter().enumerate().map(|(i,d)|(i,format!("{} [{}]",d.molecule.name,d.id))).collect())) query_info=Callback::new(move|i:usize|documents.with(|docs|docs.get(i).map(|d|crate::dotplot::QueryInfo{identity:format!("{}:{}",d.id,d.head_event()),revision:d.molecule.revision(),name:d.molecule.name.clone(),length:d.molecule.sequence.len(),alphabet:d.molecule.sequence.alphabet()}))) get_query=Callback::new(move|i:usize|documents.with(|docs|docs.get(i).map(|d|d.molecule.clone()))) epoch=workspace_epoch on_navigate=Callback::new(move|(query,index,spans):(bool,usize,Vec<(usize,usize)>)|{
                                                if workspace_restoring.get_untracked(){return;}
                                                let index=if query{index}else{active.get_untracked()};
                                                let Some((document_id,event_id,revision))=documents.with_untracked(|docs|docs.get(index).map(|d|(d.id.clone(),d.head_event().to_string(),d.molecule.revision())))else{return;};
                                                let expected_epoch=workspace_epoch.get_untracked();
                                                let apply=move|deferred:bool|{
                                                    if workspace_epoch.try_get_untracked()!=Some(expected_epoch)||workspace_restoring.try_get_untracked().unwrap_or(true){return;}
                                                    let Some(current_index)=active.try_get_untracked()else{return;};
                                                    if documents.with_untracked(|docs|docs.iter().position(|d|d.id==document_id&&d.head_event()==event_id))!=Some(current_index){return;}
                                                    if deferred&&(selected_feature.get_untracked().is_some()||selection.get_untracked().is_some()||search_regions.get_untracked().is_some()){return;}
                                                    batch(||{selected_feature.set(None);selection.set(if spans.len()==1{spans.first().copied()}else{None});search_regions.set(Some((revision,spans)));});
                                                };
                                                if active.get_untracked()==index&&has_open_document.get_untracked(){apply(false);}else{
                                                    active.set(index);
                                                    // App clears old selections when the active source changes.
                                                    // Install this source-bound selection after those effects settle.
                                                    let frame=wasm_bindgen::closure::Closure::once_into_js(move||apply(true));
                                                    if let Some(window)=web_sys::window(){let _=window.request_animation_frame(frame.unchecked_ref());}
                                                }
                                                if query{panel.set("features");}
                                            })/></Show><Show when=move||panel.get()=="compare"><crate::compare::ComparePanel molecule=molecule query_options=Memo::new(move|_|documents.with(|docs|docs.iter().enumerate().map(|(i,d)|(i,d.molecule.name.clone())).collect())) get_query=Callback::new(move|i:usize|documents.with(|docs|docs.get(i).map(|d|d.molecule.clone()))) active_source=Memo::new(move|_|documents.with(|docs|recombination_source_info(&docs[active.get()]))) source_info=Callback::new(move|index:usize|documents.with(|docs|docs.get(index).map(recombination_source_info))) resolve=resolve_retained_source epoch=workspace_epoch saved=saved_comparisons on_save=save_comparison on_reopen=reopen_comparison on_navigate=navigate_comparison/></Show>
                                        </div></aside>
                                    </div>
                                    <Show when=move||assembly_visible.get()><section class="assembly-dock" aria-label="Assembly workspace"><div class="assembly-dock-bar"><span>"Assembly design"</span><button class="text-tool" on:click=move|_|assembly_visible.set(false)>"Hide assembly"</button></div>
                                        <fieldset aria-label="Assembly source and preparation controls" disabled=move||workspace_restoring.get() style="border:0;padding:0;margin:0;min-width:0">
                                        <crate::assembly_workspace::AssemblyWorkspace draft=assembly_draft previous_draft=previous_assembly_draft epoch=workspace_epoch snapshot=contig_snapshot workspace_busy=workspace_restoring
                                            apply_saved_attachments=Callback::new(move|(before,after):(workspace::Workspace,workspace::Workspace)|{
                                                if workspace_restoring.get_untracked(){return Err("Wait for the workspace update to finish.".to_owned());}
                                                let current=contig_snapshot.run(())?;
                                                crate::assembly_saved_designs_commit::validate(&before,&after,&current)?;
                                                attachments.set(after.attachments);
                                                Ok(())
                                            })
                                            sources=Memo::new(move|_|documents.with(|docs|docs.iter().map(recombination_source_info).collect()))
                                            active_source=Memo::new(move|_|documents.with(|docs|recombination_source_info(&docs[active.get()])))
                                            selected_region=Callback::new(move|():()|{if inspected_feature.get_untracked().is_some_and(|f|f.location.spans().len()>1){return None;}selection.get_untracked()})
                                            get_source=Callback::new(move|binding:collagene_core::run_lineage::EventBinding|documents.with_untracked(|docs|{
                                                let doc=docs.iter().find(|d|d.id==binding.document_id).ok_or("The selected source document is unavailable.")?;
                                                let exact=doc.revision_molecule(&binding.event_id).map_err(|e|e.to_string())?;
                                                if exact.revision()!=binding.revision{return Err("The selected source event changed.".into());}
                                                Ok(exact)
                                            }))
                                            on_navigate=navigate_assembly on_review=review_assembly on_finalize=finalize_assembly
                                            on_products=Callback::new(move|products:Vec<Molecule>|add_documents(documents,active,workspace_restoring,products))
                                            primer_sources=primer_sources list_primers=list_saved_primers resolve_primer=resolve_saved_primer/>
                                        </fieldset>
                                        <crate::assembly_library_view::AssemblyLibraryView draft=assembly_draft epoch=workspace_epoch snapshot=contig_snapshot
                                            sources=Memo::new(move|_|documents.with(|docs|docs.iter().map(recombination_source_info).collect()))
                                            get_source=Callback::new(move|binding:collagene_core::run_lineage::EventBinding|documents.with_untracked(|docs|{
                                                let doc=docs.iter().find(|doc|doc.id==binding.document_id).ok_or("The source document is unavailable.")?;
                                                let source=doc.revision_molecule(&binding.event_id).map_err(|error|error.to_string())?;
                                                if source.revision()!=binding.revision{return Err("The selected source event changed.".into());}
                                                Ok(source)
                                            }))
                                            validate=validate_workspace_frame
                                            on_restore=restore_library_frame
                                            db workspace_busy=workspace_restoring loaded save_generation message/>
                                        <crate::assembly_pool_view::AssemblyPoolView draft=assembly_draft snapshot=contig_snapshot epoch=workspace_epoch busy=workspace_restoring on_review=review_assembly on_finalize=finalize_assembly/>
                                        <crate::assembly_batch_view::AssemblyBatchView draft=assembly_draft snapshot=contig_snapshot
                                            validate=validate_workspace_frame
                                            on_restore=restore_workspace_frame
                                            attachment_available=Callback::new(move|hash:String|attachments.with_untracked(|assets|assets.get(&hash).is_some_and(|asset|asset.kind==collagene_core::combinatorial_batch::KIND)))
                                            db workspace_busy=workspace_restoring loaded save_generation epoch=workspace_epoch message on_release=open_bulk_release/>
                                        <crate::assembly_archive_view::AssemblyArchiveView snapshot=contig_snapshot epoch=workspace_epoch selected=selected_release_manifest busy=workspace_restoring
                                            attachment_available=Callback::new(move|hash:String|attachments.with_untracked(|assets|assets.get(&hash).is_some_and(|asset|asset.kind==collagene_core::combinatorial_release::KIND)))/>
                                    </section></Show>
                                    <footer class="statusbar"><span>{move||save_status.get()}</span><Show when=move||!layout.status.get().is_empty()><span role="status">{move||layout.status.get()}</span></Show><span>{move||selection_status.get()}</span><span>{move||documents.with(|docs|format!("{} {}",docs.len(),if docs.len()==1{"document"}else{"documents"}))}</span></footer>
                                    </div>
                                </main>
                            </div>
                            </div>
                            <Show when=move||modal.get()!="">{view!{<div class="modal-backdrop"><section class="modal" class:alignment-modal=move||modal.get()=="alignments" class:read-review-modal=move||matches!(modal.get(),"plan_edit"|"plan_new")||modal.get()=="combinatorial_release"||modal.get()=="read_review"||modal.get()=="run_event"||modal.get()=="comparison_source"||modal.get()=="translation_source"||modal.get()=="pcr_source" role="dialog" aria-modal="true" aria-labelledby="dialog-title"><div class="modal-heading"><h2 id="dialog-title">{move||if matches!(modal.get(),"plan_edit"|"plan_new"){"Edit cloning plan"}else if modal.get()=="variants"{"Variant source"}else if modal.get()=="read_alignments"{"Read alignments"}else if modal.get()=="material_product"{"Open reusable PCR product"}else if modal.get()=="digest_product"{"Open virtual digest fragment"}else if modal.get()=="crispr_product"{"Open checked sequence design"}else if modal.get()=="contig_archive"{"Record read assembly"}else if modal.get()=="saved_contig"{"Saved read assembly"}else if modal.get()=="read_review"{"Review sequencing differences"}else if modal.get()=="read_review_commit"{"Save read-review decisions"}else if modal.get()=="read_check_capture"{"Save sequencing check"}else if modal.get()=="combinatorial_release"{"Create selected construct outputs"}else if modal.get()=="assembly_finalize"{"Create assembly outputs"}else if modal.get()=="record_run"{"Keep reconstructed products"}else if modal.get()=="run_details"{"Recorded reconstruction"}else if modal.get()=="run_event"{"Recorded sequence version"}else if modal.get()=="comparison_source"{"Comparison source version"}else if modal.get()=="translation_source"{"Original translation source"}else if modal.get()=="pcr_source"{"PCR template placements"}else if modal.get()=="reference_paste"{"Paste annotated reference"}else if modal.get()=="alignments"{"Alignments"}else if modal.get()=="workspace"{"Workspace"}else if modal.get()=="oligos"{"Oligo table"}else if modal.get()=="fastq"{"Import FASTQ reads"}else if modal.get()=="import"{"Open a sequence"}else if modal.get()=="url"{"Import from a sequence source"}else if modal.get()=="feature"{"Edit annotation"}else if modal.get()=="compat_import"{"Review sequence import"}else if modal.get()=="compat_snapgene"{"Export SnapGene sequence"}else if modal.get()=="compat_xdna"{"Export XDNA sequence"}else if modal.get()=="compat_gcg"{"Export GCG sequence"}else if modal.get()=="protein_genpept"{"Export GenPept protein"}else if modal.get()=="protein_uniprot"{"Export UniProt protein"}else if modal.get()=="sbol"{"SBOL design exchange"}else if modal.get()=="document"{"Document settings"}else if molecule.with(|m|m.sequence.alphabet()==Alphabet::Protein){"Edit selected residues"}else{"Edit selected bases"}}</h2><button class="icon-button" aria-label="Close dialog" disabled=move||workspace_restoring.get() on:click=move|_|{if !workspace_restoring.get_untracked(){modal.set("");pending_finalize.set(None);pending_release.set(None);pending_run_event.set(None);pending_comparison_source.set(None);pending_reference_paste.set(None);pending_read_check.set(None);pending_read_review.set(None);reviewed_check.set(None);review_seed.set(None);}}>"×"</button></div><Show when=move||!modal_error.get().is_empty()><p class="dialog-error" role="alert">{move||modal_error.get()}</p></Show><Show when=move||modal.get()=="material_product">{move||pending_material.get_untracked().map(|request|view!{<crate::digest_product::MaterialOpenPanel request snapshot=contig_snapshot validate=validate_workspace_frame on_restore=restore_workspace_frame db=db busy=workspace_restoring loaded=loaded save_generation=save_generation epoch=workspace_epoch message=message/>})}</Show><Show when=move||modal.get()=="digest_product">{move||pending_digest.get_untracked().map(|request|view!{<crate::digest_product::DigestOpenPanel request=request snapshot=contig_snapshot validate=validate_workspace_frame on_restore=contig_restore db=db busy=workspace_restoring loaded=loaded save_generation=save_generation epoch=workspace_epoch message=message/>})}</Show><Show when=move||modal.get()=="crispr_product">{move||pending_donor.get_untracked().map(|request|view!{<crate::crispr::DonorOpenPanel request=request snapshot=contig_snapshot validate=validate_workspace_frame on_restore=contig_restore db=db busy=workspace_restoring loaded=loaded save_generation=save_generation epoch=workspace_epoch message=message/>})}</Show><Show when=move||modal.get()=="contig_archive">{move||pending_assembly.get_untracked().map(|request|view!{<crate::contigs::ContigArchivePanel request=request snapshot=contig_snapshot validate=validate_workspace_frame on_restore=contig_restore db=db busy=workspace_restoring loaded=loaded save_generation=save_generation epoch=workspace_epoch message=message/>})}</Show><Show when=move||modal.get()=="saved_contig"><crate::contigs::SavedAssemblyPanel attachment_hash=inspected_assembly.get_untracked() snapshot=contig_snapshot epoch=workspace_epoch on_record=record_assembly on_navigate=navigate_read_call/></Show><Show when=move||modal.get()=="record_run">{move||last_plan.get_untracked().map(|(plan,_)|view!{<crate::run_capture::RunCapturePanel plan=plan preferred_bindings=last_plan_events.get_untracked().filter(|(revision,_)|last_plan.get_untracked().is_some_and(|(p,_)|verification::plan_revision(&p)==*revision)).map(|(_,bindings)|bindings).unwrap_or_default() snapshot=Callback::new(move|()|documents.with_untracked(|docs|stored_candidates.with_untracked(|plans|attachments.with_untracked(|assets|candidate.with_untracked(|current|workspace_snapshot(docs,active.get_untracked(),current.as_ref(),plans,assets,&workspace_name.get_untracked())))))) validate=validate_workspace_frame on_restore=restore_workspace_frame db=db busy=workspace_restoring loaded=loaded save_generation=save_generation epoch=workspace_epoch message=message/>})}</Show><Show when=move||matches!(modal.get(),"plan_edit"|"plan_new")><crate::plan_editing_host::ManualPlanEditor candidate preferences=candidate_events documents active epoch=workspace_epoch locked=Memo::new(move|_|workspace_restoring.get()) creating=modal.get_untracked()=="plan_new" on_close=Callback::new(move|()|modal.set("")) on_apply=Callback::new(move|(value,target,bindings):(transfer::Candidate,collagene_core::run_lineage::EventBinding,Vec<collagene_core::run_lineage::InputEventBinding>)|->Result<(),String>{
                                    if workspace_restoring.get_untracked(){return Err("Wait for workspace restoration.".into());}
                                    let preferred=(value.plan_revision.clone(),workspace_epoch.get_untracked(),target,bindings);
                                    batch(||{
                                        plan_generation.update(|g|*g=g.wrapping_add(1));
                                        plan_report.set(None);last_plan.set(None);last_plan_target.set(None);last_plan_events.set(None);
                                        pending_finalize.set(None);candidate_events.set(None);candidate.set(Some(value));candidate_events.set(Some(preferred));
                                        modal.set("");message.set("Edited candidate saved for review. Verify it again before creating products or protocols.".into());
                                    });Ok(())
                                })/></Show><Show when=move||modal.get()=="reference_paste">{move||pending_reference_paste.get_untracked().map(|opened|view!{<crate::reference_transfer_view::ReferencePastePanel opened guard=reference_paste_guard on_apply=apply_reference_paste on_cancel=Callback::new(move|()|{pending_reference_paste.set(None);modal.set("");})/>})}</Show><Show when=move||modal.get()=="pcr_source">{move||pending_pcr_source.get().map(|snapshot|{
                                let owner=snapshot.owner.clone();let bound=snapshot.overlay.clone();
                                let current=Memo::new(move|_|workspace_epoch.get()==bound.epoch&&!workspace_restoring.get()
                                    &&published_pcr_overlay.with(|p|p.as_ref().is_some_and(|p|p.owner==owner&&p.overlay.as_ref().is_some_and(|o|std::sync::Arc::ptr_eq(o,&bound))))
                                    &&documents.with(|docs|docs.iter().any(|d|current_event_binding(d)==bound.source)));
                                view!{<crate::pcr_source_view::PcrSourceView snapshot current on_back=Callback::new(move|()|{pending_pcr_source.set(None);modal.set("");})/>}
                            })}</Show><Show when=move||modal.get()=="translation_source">{move||pending_translation_source.get().map(|snapshot|view!{<crate::biology_source_view::TranslationSourceView snapshot epoch=workspace_epoch on_back=Callback::new(move|()|{pending_translation_source.set(None);modal.set("");})/>})}</Show><Show when=move||modal.get()=="comparison_source">{move||pending_comparison_source.get().map(|snapshot|view!{<crate::pairwise_source_view::PairwiseSourceView snapshot epoch=workspace_epoch on_back=Callback::new(move|()|{pending_comparison_source.set(None);modal.set("");})/>})}</Show><Show when=move||modal.get()=="run_event">{move||pending_run_event.get().map(|snapshot|view!{<crate::run_event_view::RunEventPanel snapshot epoch=workspace_epoch on_fork=fork_run_event on_back=back_to_run/>})}</Show><Show when=move||modal.get()=="run_details"><crate::run_capture::RunDetailsPanel snapshot=Callback::new(move|()|documents.with_untracked(|docs|stored_candidates.with_untracked(|plans|attachments.with_untracked(|assets|candidate.with_untracked(|current|workspace_snapshot(docs,active.get_untracked(),current.as_ref(),plans,assets,&workspace_name.get_untracked())))))) run_id=selected_run.get_untracked() epoch=workspace_epoch on_open_event=open_run_event/></Show><Show when=move||modal.get()=="workspace"><crate::workspace_view::WorkspacePanel snapshot=Callback::new(move|()|documents.with_untracked(|docs|stored_candidates.with_untracked(|plans|attachments.with_untracked(|assets|candidate.with_untracked(|current|workspace_snapshot(docs,active.get_untracked(),current.as_ref(),plans,assets,&workspace_name.get_untracked())))))) validate=validate_workspace_frame on_restore=restore_workspace_frame db=db busy=workspace_restoring loaded=loaded save_generation=save_generation message=message/></Show>
                                <Show when=move||modal.get()=="read_review">{move||reviewed_check.get().map(|hash|view!{<crate::read_review_view::ReadReviewPanel attachment_hash=hash reference=read_check_reference snapshot=contig_snapshot epoch=workspace_epoch on_commit=save_read_review on_read=navigate_read_call on_reference=navigate_check_reference event_available=primer_anchor_available attachment_available=read_check_attachment_available initial_request=review_seed.get_untracked()/>})}</Show>
                                <Show when=move||modal.get()=="read_review_commit">{move||pending_read_review.get().map(|pending|view!{<crate::read_review_capture::ReadReviewCapturePanel request=pending.request snapshot=read_review_snapshot validate=validate_workspace_frame on_restore=restore_read_review db busy=workspace_restoring loaded save_generation epoch=workspace_epoch message/>})}</Show>
                                <Show when=move||modal.get()=="read_check_capture">{move||pending_read_check.get().map(|request|view!{<crate::read_checks_view::ReadCheckCapturePanel request snapshot=read_check_snapshot validate=validate_workspace_frame on_restore=restore_read_check db busy=workspace_restoring loaded save_generation epoch=workspace_epoch message/>})}</Show>
                                <Show when=move||modal.get()=="combinatorial_release">{move||pending_release.get().map(|hash|view!{<crate::assembly_release_view::AssemblyReleaseView batch_hash=hash draft=assembly_draft snapshot=contig_snapshot
                                    validate=validate_workspace_frame
                                    on_restore=restore_workspace_frame
                                    db workspace_busy=workspace_restoring loaded save_generation epoch=workspace_epoch message on_success=released_bulk/>})}</Show>
                                <Show when=move||modal.get()=="assembly_finalize">{move||pending_finalize.get().map(|pending|view!{<crate::assembly_finalize::AssemblyFinalizePanel plan=pending.plan.clone() preferred_bindings=pending.sources.clone() snapshot=finalize_snapshot recovery_snapshot=contig_snapshot validate=validate_workspace_frame on_restore=restore_workspace_frame db busy=workspace_restoring loaded save_generation epoch=workspace_epoch message/>})}</Show>
                                <Show when=move||modal.get()=="alignments"><crate::alignment_view::AlignmentPanel alignments=alignments pending=pending_alignment.get_untracked() epoch=workspace_epoch on_keep={let opened_epoch=workspace_epoch.get_untracked();Callback::new(move|item:crate::alignment_view::AlignmentRevision|{if opened_epoch!=workspace_epoch.get_untracked(){return Err("The workspace changed during alignment import. Open the alignment again.".into());}keep_alignment.run(item)})} on_extract={let opened_epoch=workspace_epoch.get_untracked();Callback::new(move|molecules:Vec<Molecule>|{if opened_epoch!=workspace_epoch.get_untracked(){return Err("The workspace changed during row extraction.".into());}add_documents(documents,active,workspace_restoring,molecules)})}/></Show>
                                <Show when=move||modal.get()=="oligos"><crate::oligo_import::OligoImport on_import=import_oligos existing=existing_oligo.get_untracked() initial_text=pending_oligo.get_untracked()/></Show>
                                <Show when=move||modal.get()=="fastq"><p>"Select the quality encoding supplied by the sequencing source. Phred+33 and Phred+64 overlap; it cannot be inferred from printable characters."</p><label>"Supplied FASTQ encoding"<select aria-label="Supplied FASTQ encoding" prop:value=move||match fastq_encoding.get(){None=>"",Some(oligo_formats::QualityEncoding::Phred33)=>"33",Some(oligo_formats::QualityEncoding::Phred64)=>"64"} on:change=move|e|fastq_encoding.set(match event_target_value(&e).as_str(){"33"=>Some(oligo_formats::QualityEncoding::Phred33),"64"=>Some(oligo_formats::QualityEncoding::Phred64),_=>None})><option value="">"Choose supplied encoding"</option><option value="33">"Phred+33"</option><option value="64">"Phred+64"</option></select></label><button class="button primary" disabled=move||fastq_encoding.get().is_none() on:click=move|_|{let(Some(bytes),Some(encoding))=(pending_fastq.get_untracked(),fastq_encoding.get_untracked())else{return;};match oligo_formats::import_fastq(&bytes,encoding){Ok(report)=>match add_documents(documents,active,workspace_restoring,report.molecules){Ok(())=>{message.set(report.warnings.join(" "));pending_fastq.set(None);modal.set("");panel.set("reads");},Err(e)=>modal_error.set(e)},Err(e)=>modal_error.set(e.to_string())}}>"Import FASTQ with selected encoding"</button></Show>
                                <Show when=move||modal.get()=="variants"><p class="helper">"Declared file alleles and genotype columns are source evidence. Choose an actual DNA reference and a one-based alternative explicitly before checking or creating a design."</p><button class="button secondary" disabled=move||workspace_restoring.get() on:click=move|_|if let Err(error)=keep_variant_source.run(()){modal_error.set(error)}>"Keep original source in workspace"</button>{move||pending_variants.get().map(|(source,format,name)|view!{<crate::variants_view::VariantPanel source=source format=format name=name documents=Memo::new(move|_|documents.with(|docs|docs.iter().enumerate().filter(|(_,doc)|doc.molecule.sequence.alphabet()==Alphabet::Dna).map(|(index,doc)|(index,doc.molecule.name.clone(),format!("{}:{}:{}",workspace_epoch.get(),doc.id,doc.head_event()))).collect())) get_document=Callback::new(move|index:usize|documents.with_untracked(|docs|docs.get(index).and_then(|doc|collagene_core::model::bounds::model_preflight(&[&doc.molecule]).ok().map(|_|doc.molecule.clone())))) on_navigate=Callback::new(move|(index,parts):(usize,Vec<(usize,usize)>)|{if workspace_restoring.get_untracked(){return;}let Some((id,event,revision))=documents.with_untracked(|docs|docs.get(index).map(|doc|(doc.id.clone(),doc.head_event().to_string(),doc.molecule.revision())))else{return;};let epoch=workspace_epoch.get_untracked();let apply=move|deferred:bool|{if workspace_epoch.try_get_untracked()!=Some(epoch)||workspace_restoring.try_get_untracked().unwrap_or(true){return;}if !documents.with_untracked(|docs|docs.get(active.get_untracked()).is_some_and(|doc|doc.id==id&&doc.head_event()==event)){return;}if deferred&&(selection.get_untracked().is_some()||selected_feature.get_untracked().is_some()||search_regions.get_untracked().is_some()){return;}batch(||{selected_feature.set(None);selection.set(if parts.len()==1{parts.first().copied()}else{None});search_regions.set(Some((revision,parts)));});};modal.set("");if active.get_untracked()==index&&has_open_document.get_untracked(){apply(false);}else{active.set(index);let frame=wasm_bindgen::closure::Closure::once_into_js(move||apply(true));if let Some(window)=web_sys::window(){let _=window.request_animation_frame(frame.unchecked_ref());}}}) on_preview=open_variant_preview/>})}</Show><Show when=move||modal.get()=="read_alignments"><p class="helper">"Keep the original source in this workspace to reopen it later. Mapping claims remain unverified until checked against an explicitly selected DNA document."</p><button class="button secondary" disabled=move||workspace_restoring.get() on:click=move|_|if let Err(error)=keep_read_source.run(()){modal_error.set(error)}>"Keep original source in workspace"</button>{move||pending_read_alignment.get().map(|(source,format,name)|view!{
                                    <crate::read_alignments_view::ReadAlignmentPanel source=source format=format name=name documents=Memo::new(move|_|documents.with(|docs|docs.iter().enumerate().filter(|(_,doc)|doc.molecule.sequence.alphabet()==Alphabet::Dna).map(|(index,doc)|(index,doc.molecule.name.clone(),format!("{}:{}:{}",workspace_epoch.get(),doc.id,doc.head_event()))).collect())) get_document=Callback::new(move|index:usize|documents.with_untracked(|docs|docs.get(index).and_then(|doc|collagene_core::model::bounds::model_preflight(&[&doc.molecule]).ok().map(|_|doc.molecule.clone())))) on_navigate=Callback::new(move|(index,parts):(usize,Vec<(usize,usize)>)|{if workspace_restoring.get_untracked(){return;}let Some((id,event,revision))=documents.with_untracked(|docs|docs.get(index).map(|doc|(doc.id.clone(),doc.head_event().to_string(),doc.molecule.revision())))else{return;};
                                        let epoch=workspace_epoch.get_untracked();
                                        let apply=move|deferred:bool|{if workspace_epoch.try_get_untracked()!=Some(epoch)||workspace_restoring.try_get_untracked().unwrap_or(true){return;}if !documents.with_untracked(|docs|docs.get(active.get_untracked()).is_some_and(|doc|doc.id==id&&doc.head_event()==event)){return;}if deferred&&(selection.get_untracked().is_some()||selected_feature.get_untracked().is_some()||search_regions.get_untracked().is_some()){return;}batch(||{selected_feature.set(None);selection.set(if parts.len()==1{parts.first().copied()}else{None});search_regions.set(Some((revision,parts)));});};
                                        modal.set("");if active.get_untracked()==index&&has_open_document.get_untracked(){apply(false);}else{active.set(index);let frame=wasm_bindgen::closure::Closure::once_into_js(move||apply(true));if let Some(window)=web_sys::window(){let _=window.request_animation_frame(frame.unchecked_ref());}}})/>
                                })}</Show><Show when=move||modal.get()=="import"><p>"Paste GenBank, GenPept, UniProt, EMBL, FASTA, or plain DNA. No project setup needed."</p><button class="button secondary" disabled=move||workspace_restoring.get() on:click=move|_|{let text=import_text.get_untracked();if !text.trim_start_matches('\u{feff}').trim_start().starts_with('>'){modal_error.set("Paste FASTA records beginning with > before choosing their alphabet.".into());return;}pending_compat_import.set(Some((std::sync::Arc::new(text.into_bytes()),"Pasted FASTA".into(),crate::sequence_compat_import::Profile::Fasta(None))));modal_error.set(String::new());modal.set("compat_import");}>"Choose FASTA alphabet"</button><textarea aria-label="Sequence input" rows="12" prop:value=move||import_text.get() on:input=move|e|import_text.set(event_target_value(&e)) placeholder=">my_sequence\nACGT…"></textarea><button class="button primary" on:click=move|_|{let text=import_text.get_untracked();if text.starts_with("##fileformat=VCF"){if text.len()>16*1024*1024{modal_error.set("Variant source exceeds the 16 MiB browser profile.".into());return;}pending_variants.set(Some((std::sync::Arc::new(text.into_bytes()),collagene_core::variants_session::FileFormat::Vcf,"pasted.vcf".into())));modal_error.set(String::new());modal.set("variants");return;}if [b"@HD\t".as_slice(),b"@SQ\t".as_slice(),b"@RG\t".as_slice(),b"@PG\t".as_slice(),b"@CO\t".as_slice()].iter().any(|prefix|text.as_bytes().starts_with(prefix)){if text.len()>collagene_core::read_alignment_session::MAX_SOURCE_BYTES{modal_error.set("This browser read-alignment profile accepts sources up to 16 MiB.".into());return;}pending_read_alignment.set(Some((std::sync::Arc::new(text.into_bytes()),collagene_core::read_alignment_session::FileFormat::Sam,"pasted.sam".into())));modal_error.set(String::new());modal.set("read_alignments");return;}if text.as_bytes().first()==Some(&b'@'){pending_fastq.set(Some(std::sync::Arc::new(text.into_bytes())));fastq_encoding.set(None);modal.set("fastq");return;}match formats::import_text(&text){Ok(report)=>{if let Err(error)=add_documents(documents,active,workspace_restoring,report.molecules){modal_error.set(error);return;}message.set(report.warnings.join(" "));modal.set("");import_text.set(String::new());},Err(e)=>modal_error.set(e.to_string())}}>"Open sequence"</button></Show><Show when=move||modal.get()=="url"><crate::source_import::SourceImport on_import={let epoch=workspace_epoch.get_untracked();Callback::new(move |report:formats::ImportReport|{if epoch!=workspace_epoch.get_untracked(){return Err("The workspace changed during retrieval. Import the source again in this workspace.".into());}add_documents(documents,active,workspace_restoring,report.molecules)?;message.set(report.warnings.join(" "));modal.set("");Ok(())})}/></Show><Show when=move||modal.get()=="feature"><crate::editors::FeatureEditor molecule=molecule.get_untracked() index=feature_for_editor.get_untracked() selection=selection.get_untracked() initial_location=feature_initial_location.get_untracked() on_commit=Callback::new(move|(base,edited):(String,Molecule)|{commit_document(documents,active,workspace_restoring,&base,edited)?;modal.set("");Ok(())})/></Show><Show when=move||matches!(modal.get(),"protein_genpept"|"protein_uniprot")>{move||{let profile=if modal.get()=="protein_genpept"{collagene_core::protein_formats::Profile::GenPept}else{collagene_core::protein_formats::Profile::UniProt};let guard=export_context_guard(documents,active,workspace_restoring,workspace_epoch);view!{<crate::protein_export::ProteinExportDialog source=molecule.get_untracked() profile guard on_saved=Callback::new(move|notice:String|{if message.try_get_untracked().is_some(){message.set(notice);}})/>}}}</Show>
        <Show when=move||modal.get()=="compat_import">{move||pending_compat_import.get().map(|(source,name,profile)|{
            let selected=source.clone();let context=export_context_guard(documents,active,workspace_restoring,workspace_epoch);
            let guard=Callback::new(move|()|{context.run(())?;if modal.get_untracked()!="compat_import"{return Err("The import review closed.".into());}if !pending_compat_import.with_untracked(|p|p.as_ref().is_some_and(|(current,_,_)|std::sync::Arc::ptr_eq(current,&selected))){return Err("The selected import source changed; reopen its review.".into());}Ok(())});
            view!{<crate::sequence_compat_import::SequenceCompatImport source name profile guard on_import=Callback::new(move|report:formats::ImportReport|{let notice=report.warnings.join(" ");add_documents(documents,active,workspace_restoring,report.molecules)?;message.set(notice);modal.set("");pending_compat_import.set(None);Ok(())})/>}
        })}</Show>
        <Show when=move||modal.get()=="map_svg">{move||{
            let context=export_context_guard(documents,active,workspace_restoring,workspace_epoch);
            let snapshot=Callback::new(move|()|->Result<crate::map_export_model::Snapshot,String>{
                context.run(())?;
                if modal.try_get_untracked()!=Some("map_svg") || !loaded.get_untracked() {
                    return Err("The map export view closed or its workspace is unavailable.".into());
                }
                let map=renderer.get_value().ok_or("Wait for the source map to finish opening.")?.export_view()?;
                let (source,event)=documents.with_untracked(|docs|docs.get(active.get_untracked()).map(|doc|(doc.molecule.clone(),current_event_binding(doc)))).ok_or("The selected map document is unavailable.")?;
                let hidden=annotation_hidden.get_untracked().into_iter().collect::<Vec<_>>();
                if map.source_revision!=event.revision || map.interaction_context!=restriction_identity.get_untracked()
                    || map.hidden_features!=hidden || map.selected_feature!=selected_feature.get_untracked()
                    || map.highlighted_spans!=highlighted_spans.get_untracked() {
                    return Err("Map source or display settings are still updating; prepare the export again.".into());
                }
                let restriction_layer=if restrictions.visible.get_untracked(){
                    if restrictions.busy.get_untracked(){return Err("Wait for the displayed restriction scan before exporting it.".into());}
                    let scan=restrictions.scan.get_untracked().ok_or("No exact restriction scan is ready; scan or hide restriction sites before exporting.")?;
                    if scan.revision!=event.revision || scan.context!=restriction_identity.get_untracked(){return Err("Restriction scan belongs to a different source event.".into());}
                    let sites=restrictions.displayed.get_untracked().ok_or("Displayed restriction sites are still updating.")?;
                    if sites.as_ref()!=map.restriction_sites.as_ref(){return Err("Restriction map filters are still updating; prepare the export again.".into());}
                    Some(crate::map_export_model::RestrictionLayer{
                        source_revision:event.revision.clone(),
                        display_context:{let settings=serde_json::to_vec(&(scan.context.as_str(),restrictions.names.get_untracked(),restrictions.count_filter.get_untracked(),restrictions.region_filter.get_untracked(),format!("{:?}",restrictions.basis.get_untracked()))).map_err(|error|error.to_string())?;blake3::hash(&settings).to_hex().to_string()},
                        sites:sites.as_ref().clone(),selected_site:map.selected_site,
                    })
                }else{
                    if !map.restriction_sites.is_empty(){return Err("Restriction visibility is still updating; prepare the export again.".into());}
                    None
                };
                let view=crate::map_export_model::ViewSettings{
                    selection:selection.get_untracked().map(|(start,end)|crate::map_export_model::Selection{start,end,wrap:start>end}),
                    rotation_radians:map.rotation_radians,zoom:map.zoom,selected_feature:map.selected_feature,
                    highlighted_spans:map.highlighted_spans,hidden_features:map.hidden_features,
                    display_context:serde_json::to_string(&(map.interaction_context,editor_view.get_untracked(),map.overlay_context)).map_err(|error|error.to_string())?,
                };
                Ok(crate::map_export_model::Snapshot{source,event,epoch:workspace_epoch.get_untracked(),view,restrictions:restriction_layer})
            });
            view!{<crate::map_export::MapExportDialog snapshot render=crate::map_export::worker_renderer() on_saved=Callback::new(move|notice:String|{if message.try_get_untracked().is_some(){message.set(notice);}})/>}
        }}</Show>
        <Show when=move||matches!(modal.get(),"compat_snapgene"|"compat_xdna"|"compat_gcg")>{move||{let expected_modal=modal.get();let mode=match expected_modal{"compat_snapgene"=>crate::sequence_export_model::Mode::SnapGene,"compat_xdna"=>crate::sequence_export_model::Mode::Xdna,_=>crate::sequence_export_model::Mode::Gcg};let context=export_context_guard(documents,active,workspace_restoring,workspace_epoch);let guard=Callback::new(move|()|{context.run(())?;if modal.try_get_untracked().as_ref()!=Some(&expected_modal){return Err("The export view closed or changed format.".into());}Ok(())});view!{<crate::sequence_compat_export::SequenceCompatExport source=molecule.get_untracked() mode guard on_saved=Callback::new(move|notice:String|{if message.try_get_untracked().is_some(){message.set(notice);}})/>}}}</Show>
        <Show when=move||modal.get()=="sbol">{move||{
            let source=pending_sbol.get();
            let expected_source=source.clone();
            let context_guard=export_context_guard(documents,active,workspace_restoring,workspace_epoch);
            let guard=Callback::new(move|()|{context_guard.run(())?;let current=pending_sbol.try_get_untracked().ok_or("The SBOL source view closed.")?;let same=match(&expected_source,&current){(None,None)=>true,(Some(a),Some(b))=>std::sync::Arc::ptr_eq(a,b),_=>false};if !same{return Err("The selected SBOL source changed; reopen its review.".into());}Ok(())});
            let document_id=documents.with_untracked(|docs|docs.get(active.get_untracked()).map(|d|d.id.clone()).unwrap_or_else(||"no-active-document".into()));
            view!{<crate::sbol_exchange::SbolExchangeDialog bytes=source current=std::sync::Arc::new(molecule.get_untracked()) document_id epoch=workspace_epoch guard on_close=Callback::new(move|()|{pending_sbol.set(None);modal.set("");}) on_import=Callback::new(move|(product,attachment):(Molecule,Option<collagene_core::workspace::WorkspaceAttachment>)|->Result<(),String>{
                guard.run(())?;
                let mut docs=documents.get_untracked();
                if docs.len()>=1000 || docs.iter().map(|d|document_bytes(&d.molecule)).sum::<usize>()+document_bytes(&product)>MAX_WORKSPACE_BYTES {return Err("This SBOL projection exceeds the workspace document/memory limit.".into());}
                let index=docs.len();let additions=open_histories(vec![product],&docs)?;docs.extend(additions);
                let mut assets=attachments.get_untracked();
                if let Some(asset)=attachment {
                    let hash=blake3::hash(&asset.bytes).to_hex().to_string();
                    if let Some(existing)=assets.get(&hash) {if existing.kind!=asset.kind || existing.media_type!=asset.media_type || existing.bytes!=asset.bytes {return Err("The SBOL source collides with an existing attachment; no workspace data changed.".into());}}
                    else {assets.insert(hash,asset);}
                }
                let frame=workspace_snapshot(&docs,index,candidate.get_untracked().as_ref(),&stored_candidates.get_untracked(),&assets,&workspace_name.get_untracked())?;
                frame.validate().map_err(|e|e.to_string())?;
                batch(||{documents.set(docs);active.set(index);attachments.set(assets);pending_sbol.set(None);});
                message.set("Reviewed SBOL component imported with its original graph when available; source annotations remain documentary.".into());modal.set("");Ok(())
            })/>}
        }}</Show><Show when=move||modal.get()=="document"><crate::editors::DocumentEditor molecule=molecule.get_untracked() selection=selection.get_untracked() on_commit=Callback::new(move|(base,edited):(String,Molecule)|{commit_document(documents,active,workspace_restoring,&base,edited)?;modal.set("");Ok(())})/></Show><Show when=move||modal.get()=="edit"><p>{move||selection.get().map(|(a,b)|if a==b{format!("Insert at boundary {a}. Affected annotations will be marked for review.")}else{format!("Replace {} {}–{}; clear the field to delete. Affected annotations will be marked for review.",molecule.with(|m|m.sequence.alphabet().selected_symbols()),a+1,b)}).unwrap_or_else(||format!("Select {} in the sequence view first.",molecule.with(|m|m.sequence.alphabet().selected_symbols())))}</p><textarea aria-label="Replacement sequence" rows="6" prop:value=move||edit_text.get() on:input=move|e|edit_text.set(event_target_value(&e))></textarea><button class="button primary" disabled=move||selection.get().is_none() on:click=move|_|if let Some((a,b))=selection.get_untracked(){match molecule.get_untracked().replace(a,b,&edit_text.get_untracked()){Ok(edited)=>{let base=molecule.get_untracked().revision();if let Err(error)=commit_document(documents,active,workspace_restoring,&base,edited){modal_error.set(error);return;}selection.set(None);plan_report.set(None);last_plan.set(None);last_plan_target.set(None);modal.set("");},Err(e)=>modal_error.set(e.to_string())}}>"Apply edit"</button></Show></section></div>}.into_any()}</Show>
                        </div>
                    }
}
