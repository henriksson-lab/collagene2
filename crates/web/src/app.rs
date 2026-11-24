use crate::{
    renderer::{MapRenderer, css_color},
    storage, worker,
};
use collagene_core::{
    analysis::Primer,
    formats, protocol, transfer,
    verification::{self, DesignPlan, Operation, PlanInput, Source, VerificationReport},
    *,
};
use leptos::{html, prelude::*};
use serde::{Deserialize, Serialize};
use wasm_bindgen::JsCast;
use wasm_bindgen_futures::{JsFuture, spawn_local};
use web_sys::{HtmlCanvasElement, HtmlInputElement};

type OpenDocument = collagene_core::history::DocumentHistory;
const MAX_WORKSPACE_BYTES: usize = 64 * 1024 * 1024;
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
}
fn restore_workspace_state(
    frame: &collagene_core::workspace::Workspace,
) -> Result<RestoredWorkspace, String> {
    let documents = OpenDocument::restore_all(frame).map_err(|e| e.to_string())?;
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
    Ok(RestoredWorkspace {
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
#[derive(Serialize, Deserialize)]
struct Snapshot {
    version: u32,
    documents: Vec<Molecule>,
    active: usize,
    #[serde(default)]
    candidate: Option<transfer::Candidate>,
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

pub(crate) fn download(name: &str, mime: &str, bytes: &[u8]) -> Result<(), String> {
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
            enzymes,
            profile_revision,
            fragment_index,
        } => format!(
            "Checked restriction digest of {}\nEnzymes: {}\nPinned positional supplier profile revision: {}\nSelect fragment {} in ascending reference top-cut order. Unknown or impaired sites withhold undetermined products; blocked sites remain uncut under this model. Reagent conditions and experimental success remain outside the check.",
            source(template),
            enzymes.join(", "),
            profile_revision,
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
        } => format!(
            "Ordered substrates: {}\nProduct topology: {}\nRetained nicks explicitly allowed: {}",
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
            allow_nicked
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
    let inspected_assembly = RwSignal::new(String::new());
    let selected_run = RwSignal::new(String::new());
    let run_page = RwSignal::new(0usize);
    let oligo_tables = RwSignal::new(Vec::<
        std::sync::Arc<collagene_core::oligo_formats::OligoTable>,
    >::new());
    let libraries = RwSignal::new(Vec::<crate::feature_tools::LibraryRevision>::new());
    let alignments = RwSignal::new(Vec::<crate::alignment_view::AlignmentRevision>::new());
    let pending_alignment = RwSignal::new(None::<crate::alignment_view::PendingAlignment>);
    let method_config = RwSignal::new(crate::cloning_methods::MethodsState::default());
    let annotation_source = RwSignal::new(String::new());
    let annotation_format =
        RwSignal::new(collagene_core::annotation_formats::AnnotationFormat::Gff3);
    let existing_oligo = RwSignal::new(None::<collagene_core::oligo_formats::OligoTable>);
    let pending_oligo = RwSignal::new(None::<(String, oligo_formats::Delimiter)>);
    let pending_fastq = RwSignal::new(None::<std::sync::Arc<Vec<u8>>>);
    let fastq_encoding = RwSignal::new(None::<oligo_formats::QualityEncoding>);
    let stored_candidates = RwSignal::new(Vec::<transfer::Candidate>::new());
    let active = RwSignal::new(0usize);
    let molecule = Memo::new(move |_| {
        documents.with(|docs| docs[active.get().min(docs.len() - 1)].molecule.clone())
    });
    let molecule_revision = Memo::new(move |_| molecule.with(Molecule::revision));
    let selection = RwSignal::new(None::<(usize, usize)>);
    let selected_feature = RwSignal::new(None::<usize>);
    let feature_page = RwSignal::new(0usize);
    let feature_for_editor = RwSignal::new(None::<usize>);
    let inspected_feature = Memo::new(move |_| {
        molecule.with(|m| {
            selected_feature
                .get()
                .and_then(|i| m.features.get(i).cloned())
        })
    });
    Effect::new(move |_| {
        if let Some(index) = selected_feature.get() {
            feature_page.set(index / 500);
        }
    });
    let search_regions = RwSignal::new(None::<(String, Vec<(usize, usize)>)>);
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
    let panel = RwSignal::new("features");
    let message = RwSignal::new(String::new());
    let modal = RwSignal::new("");
    let modal_error = RwSignal::new(String::new());
    let import_text = RwSignal::new(String::new());
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
    let explorer_documents = Memo::new(move |_| {
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
                        source: provider.map(|provider| {
                            accession
                                .map(|id| format!("{provider} · {id}"))
                                .unwrap_or_else(|| provider.clone())
                        }),
                    }
                })
                .collect::<Vec<_>>()
        })
    });
    let explorer_active = Memo::new(move |_| {
        documents.with(|docs| {
            docs.get(active.get())
                .map(|doc| doc.id.clone())
                .unwrap_or_default()
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
            if active.get_untracked() != index {
                active.set(index);
            }
        }
    });

    let db = StoredValue::new_local(None::<web_sys::IdbDatabase>);
    let loaded = RwSignal::new(false);
    let save_status = RwSignal::new("Opening local storage…".to_string());
    let save_generation = RwSignal::new(0u64);
    let pcr_forward = RwSignal::new(String::new());
    let pcr_reverse = RwSignal::new(String::new());
    let pcr_forward_tail = RwSignal::new(String::new());
    let pcr_reverse_tail = RwSignal::new(String::new());
    let pcr_generation = RwSignal::new(0u32);
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
    let pdf_busy = RwSignal::new(false);
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
    });
    let apply_restored = Callback::new(move |state: RestoredWorkspace| {
        workspace_epoch.update(|epoch| *epoch = epoch.wrapping_add(1));
        pcr_generation.update(|g| *g = g.wrapping_add(1));
        existing_oligo.set(None);
        pending_oligo.set(None);
        pending_fastq.set(None);
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
        inspected_assembly.set(String::new());
        run_page.set(0);
        pending_alignment.set(None);
        annotation_source.set(String::new());
        method_config.set(crate::cloning_methods::MethodsState::default());
        candidate_events.set(None);
        candidate.set(state.candidate);
        plan_generation.update(|g| *g = g.wrapping_add(1));
        plan_report.set(None);
        last_plan.set(None);
        selection.set(None);
        selected_feature.set(None);
        modal.set("");
    });
    spawn_local(async move {
        match storage::open().await {
            Ok(database) => {
                db.set_value(Some(database.clone()));
                match storage::load(&database).await {
                    Ok(Some(bytes)) => {
                        let frame = if bytes.starts_with(b"PK") {
                            decode_workspace(bytes, 0).await
                        } else {
                            (|| {
                                let snapshot: Snapshot = serde_json::from_slice(&bytes)
                                    .map_err(|_| "Invalid legacy browser snapshot".to_string())?;
                                if snapshot.version != 1 || snapshot.documents.is_empty() {
                                    return Err("Unsupported legacy browser snapshot".into());
                                }
                                let index = snapshot.active.min(snapshot.documents.len() - 1);
                                let mut frame =
                                    collagene_core::workspace::Workspace::from_documents(
                                        "Recovered browser workspace",
                                        snapshot.documents,
                                        index,
                                    )
                                    .map_err(|e| e.to_string())?;
                                if let Some(c) = snapshot.candidate {
                                    for input in &c.plan.inputs {
                                        frame
                                            .revisions
                                            .entry(input.revision.clone())
                                            .or_insert_with(|| input.molecule.clone());
                                    }
                                    frame.active_candidate = Some(c.plan_revision.clone());
                                    frame.candidates.push(c);
                                }
                                frame.validate().map_err(|e| e.to_string())?;
                                Ok(frame)
                            })()
                        };
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
                // write transaction; IndexedDB then serializes same-store transactions.
                if generation != save_generation.get_untracked()
                    || workspace_restoring.get_untracked()
                {
                    return;
                }
                let result = match result {
                    Ok(bytes) => storage::save(&database, &bytes).await,
                    Err(error) => Err(error),
                };
                if generation != save_generation.get_untracked()
                    || workspace_restoring.get_untracked()
                {
                    return;
                }
                match result {
                    Ok(()) => save_status.set("Saved in this browser".into()),
                    Err(error) => {
                        save_status.set("Autosave failed".into());
                        message.set(error);
                    }
                }
            });
        }
    });
    Effect::new(move |_| {
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
                        Ok(map) => renderer.set_value(Some(map)),
                        Err(error) => render_error.set(error),
                    }
                }
            }
        });
    });
    Effect::new(move |_| {
        let _ = active.get();
        plan_generation.update(|g| *g = g.wrapping_add(1));
        selected_feature.set(None);
        feature_page.set(0);
        selection.set(None);
        search_regions.set(None);
        plan_report.set(None);
        last_plan.set(None);
    });
    Effect::new(move |_| {
        let _ = molecule.get().revision();
        selected_feature.set(None);
        feature_page.set(0);
        let _ = documents.with(|docs| {
            docs.iter()
                .map(|d| {
                    (
                        d.id.clone(),
                        d.head_event().to_string(),
                        d.molecule.revision(),
                    )
                })
                .collect::<Vec<_>>()
        });
        plan_generation.set(plan_generation.get_untracked().wrapping_add(1));
        selection.set(None);
        search_regions.set(None);
        plan_report.set(None);
        last_plan.set(None);
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
                let extension = file
                    .name()
                    .rsplit('.')
                    .next()
                    .unwrap_or("")
                    .to_ascii_lowercase();
                if matches!(
                    extension.as_str(),
                    "aln" | "clustal" | "sto" | "stockholm" | "afa" | "afasta" | "cg2aln"
                ) || bytes.starts_with(b"CLUSTAL")
                    || bytes.starts_with(b"# STOCKHOLM")
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
                if matches!(extension.as_str(), "gff3" | "gff" | "bed") {
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
                            annotation_format.set(if extension == "bed" {
                                collagene_core::annotation_formats::AnnotationFormat::Bed
                            } else {
                                collagene_core::annotation_formats::AnnotationFormat::Gff3
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
            Ok(bytes) => {
                if let Err(e) = download(
                    &format!("{}.cg2", mol.name),
                    "application/octet-stream",
                    &bytes,
                ) {
                    message.set(e)
                }
            }
            Err(e) => message.set(e.to_string()),
        }
    };
    let do_pcr = move |_: web_sys::MouseEvent| {
        let primers = Primer::new(
            "Forward",
            &pcr_forward.get_untracked(),
            &pcr_forward_tail.get_untracked(),
        )
        .and_then(|f| {
            Primer::new(
                "Reverse",
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
        let revision = mol.revision();
        let pcr_inputs = (
            pcr_forward.get_untracked(),
            pcr_reverse.get_untracked(),
            pcr_forward_tail.get_untracked(),
            pcr_reverse_tail.get_untracked(),
        );
        let generation = pcr_generation.get_untracked().wrapping_add(1);
        pcr_generation.set(generation);
        let request = protocol::Request {
            id: generation,
            revision: revision.clone(),
            job: protocol::AnalysisJob::PcrDocuments {
                molecule: mol,
                forward,
                reverse,
            },
        };
        message.set("Computing PCR products in background…".into());
        spawn_local(async move {
            let result = worker::analyze(request).await;
            if molecule.get_untracked().revision() != revision
                || pcr_generation.get_untracked() != generation
                || pcr_inputs
                    != (
                        pcr_forward.get_untracked(),
                        pcr_reverse.get_untracked(),
                        pcr_forward_tail.get_untracked(),
                        pcr_reverse_tail.get_untracked(),
                    )
            {
                return;
            }
            match result {
                Ok(protocol::Response {
                    result: protocol::AnalysisResult::PcrDocuments(products),
                    ..
                }) => {
                    let count = products.len();
                    let findings = products
                        .iter()
                        .flat_map(|p| p.findings.iter())
                        .take(20)
                        .cloned()
                        .collect::<Vec<_>>()
                        .join(" ");
                    let docs = products.into_iter().map(|p| p.molecule).collect();
                    if let Err(error) = add_documents(documents, active, workspace_restoring, docs)
                    {
                        message.set(error);
                        return;
                    }
                    message.set(format!("{count} exact-annealing PCR product geometry/ies with projected annotations. Reaction feasibility is not evaluated. {findings}"));
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
    let run_plan =
        move |plan: DesignPlan,
              target: Molecule,
              exact_sources: Option<Vec<collagene_core::run_lineage::InputEventBinding>>| {
            let expected_epoch = workspace_epoch.get_untracked();
            let expected_target = documents.with_untracked(|docs| {
                docs.get(active.get_untracked())
                    .map(|d| (d.id.clone(), d.head_event().to_string()))
            });
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
                        docs.get(active.get_untracked())
                            .map(|d| (d.id.clone(), d.head_event().to_string()))
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
        if let Some((_, epoch, target, sources)) = &preferences {
            let valid = *epoch == workspace_epoch.get_untracked()
                && documents.with_untracked(|docs| {
                    docs.get(active.get_untracked())
                        .is_some_and(|d| current_event_binding(d) == *target)
                        && sources
                            .iter()
                            .all(|s| docs.iter().any(|d| current_event_binding(d) == s.event))
                });
            if !valid {
                message.set("Reference-exchange source or target events changed. Inspect the current sources again before verification.".into());
                return;
            }
        }
        let is_reference = candidate
            .plan
            .steps
            .iter()
            .any(|step| matches!(step, Operation::RecombineReference { .. }));
        let captured = if is_reference {
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
            Ok(bytes) if bytes.len()<=8*1024*1024=>{if let Err(e)=download("planning-workspace.cg2design","application/octet-stream",&bytes) {message.set(e);}},
            Ok(_)=>message.set("The initial MCP exchange supports snapshots up to 8 MiB. Open only the target and needed inputs, then export again.".into()),
            Err(e)=>message.set(e.to_string()),
        }
    };
    let export_pdf = move |_: web_sys::MouseEvent| {
        if pdf_busy.get_untracked() {
            return;
        }
        let Some((plan, target)) = last_plan.get_untracked() else {
            return;
        };
        let revision = target.revision();
        let plan_revision = verification::plan_revision(&plan);
        let generation = plan_generation.get_untracked();
        pdf_busy.set(true);
        message.set("Generating PDF protocol in background…".into());
        spawn_local(async move {
            let result = worker::analyze(protocol::Request {
                id: generation,
                revision: revision.clone(),
                job: protocol::AnalysisJob::ExportProtocolPdf { target, plan },
            })
            .await;
            pdf_busy.set(false);
            if generation != plan_generation.get_untracked()
                || molecule.get_untracked().revision() != revision
                || !last_plan.with_untracked(|p| {
                    p.as_ref()
                        .is_some_and(|(p, _)| verification::plan_revision(p) == plan_revision)
                })
            {
                message.set(
                    "The design changed during PDF generation; export the current design again."
                        .into(),
                );
                return;
            }
            match result {
                Ok(protocol::Response {
                    result: protocol::AnalysisResult::ProtocolPdf(report),
                    ..
                }) if report.plan_revision == plan_revision
                    && report.target_revision == revision =>
                {
                    match download("cloning-protocol.pdf", "application/pdf", &report.bytes) {
                        Ok(()) => message.set(format!(
                            "Exported {} PDF pages · {:?}. {}",
                            report.pages,
                            report.status,
                            report.warnings.join(" ")
                        )),
                        Err(error) => message.set(error),
                    }
                }
                Ok(protocol::Response {
                    result: protocol::AnalysisResult::Error(error),
                    ..
                })
                | Err(error) => message.set(error),
                _ => message.set("Unexpected PDF worker response".into()),
            }
        });
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
                    Operation::CheckedDigest {
                        template: source,
                        enzymes: names(enzymes),
                        profile_revision: modified_duplex::profile_revision(),
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
                steps.push(Operation::Ligate {
                    fragments: references,
                    circular: target.topology == Topology::Circular,
                    allow_nicked: allow_nicked.get_untracked(),
                });
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
        if active.get_untracked() == index {
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
            <Show when=move||workspace_restoring.get()><div role="status" aria-live="polite" style="position:fixed;inset:0;z-index:10000;background:#f5f8f4ed;display:grid;place-items:center;pointer-events:all"><p>"Updating workspace…"</p></div></Show>
            <div class="application workspace-shell" inert=move||workspace_restoring.get()>
                <header class="topbar">
                    <a class="brand" href="#" aria-label="Collagene workspace"><svg viewBox="0 0 32 32" aria-hidden="true"><circle cx="16" cy="16" r="11" fill="none" stroke="currentColor" stroke-width="3"/><path d="M9 7L23 25M7 10L11 6M21 26L25 22" stroke="currentColor" stroke-width="3"/></svg><span>"collagene"</span><span class="version">"2"</span></a>
                    <div class="top-divider"></div><span class="workspace-label">{move||workspace_name.get()}</span>
                    <div class="top-actions"><button class="button secondary" on:click=move|_|{modal_error.set(String::new());modal.set("workspace");}>"Workspace"</button><button class="button secondary" on:click=move|_|{modal_error.set(String::new());modal.set("import");}>"Paste sequence"</button><button class="button primary" on:click=move|_|if let Some(input)=file_input.get(){input.click();}>"Open GenBank"</button></div>
                    <input node_ref=file_input type="file" multiple accept=".gb,.gbk,.gbff,.genbank,.ape,.fa,.fasta,.fna,.faa,.txt,.cg2,.cg2design,.cg2plan,.dna,.embl,.emb,.ab1,.abi,.scf,.fastq,.fq,.csv,.tsv,.gff3,.gff,.bed,.xml,.collagene,.seqfile,.xdna,.xdgn,.xprt,.aln,.clustal,.sto,.stockholm,.afa,.afasta,.cg2aln" style="display:none" on:change=open_file/>
                </header>
                <div class="body-layout" class:library-hidden=move||!show_library.get()>
                    <nav class="workspace-rail" aria-label="Workspace navigation">
                        <button class="rail-button" class:active=move||show_library.get() title="Documents" aria-label="Documents" aria-expanded=move||show_library.get().to_string() on:click=move|_|show_library.update(|visible|*visible=!*visible)>
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
                    <crate::document_explorer::DocumentExplorer documents=explorer_documents active_document=explorer_active locked=explorer_locked on_activate=explorer_activate on_open=Callback::new(move|()|if let Some(input)=file_input.get(){input.click();})/>
                    <main>
                        <crate::document_explorer::DocumentTabs documents=explorer_documents active_document=explorer_active locked=explorer_locked on_activate=explorer_activate/>
                        <div class="active-document" id="active-sequence-document" role="tabpanel" aria-label="Active sequence document" tabindex="0">
                        <div class="document-header"><div><div class="eyebrow">"SEQUENCE"</div><h1>{move||molecule.get().name}</h1><div class="document-meta"><span>{move||format!("{} {}",molecule.get().sequence.len(),molecule.get().sequence.alphabet().short_unit())}</span><span>"·"</span><span>{move||format!("{} {}",if molecule.get().topology==Topology::Circular{"Circular"}else{"Linear"},match molecule.get().sequence.alphabet(){Alphabet::Dna=>"DNA",Alphabet::Rna=>"RNA",Alphabet::Protein=>"protein"})}</span><span>"·"</span><span>{move||molecule.get().sequence.gc_fraction().map(|v|format!("{:.1}% GC",v*100.)).unwrap_or_else(||"GC unavailable".into())}</span></div><Show when=move||molecule.get().metadata.contains_key("source:provider")><p class="source-provenance">{move||{let m=molecule.get();let mut parts=vec![format!("Source: {}",m.metadata.get("source:provider").map(String::as_str).unwrap_or("unknown")),m.metadata.get("source:accession_version").or_else(||m.metadata.get("source:identifier")).cloned().unwrap_or_else(||"record".into())];if let Some(id)=m.metadata.get("source:sequence_id"){parts.push(format!("sequence {id}"));}if let Some(category)=m.metadata.get("source:addgene_completeness"){parts.push(format!("{category} record (provider label)"));}parts.push(format!("retrieved {}",m.metadata.get("source:retrieved_at").map(String::as_str).unwrap_or("unknown")));parts.join(" · ")}}</p></Show></div>
                            <div class="document-actions"><button class="button secondary" on:click=save_file>"Save .cg2"</button><button class="button secondary" on:click=move|_|{let mol=molecule.get_untracked();match formats::export_with_report(&mol,formats::TextFormat::GenBank){Ok(report)=>if let Err(e)=download(&format!("{}.gb",mol.name),"text/plain",report.text.as_bytes()){message.set(e)}else{message.set(report.warnings.join(" "));},Err(e)=>message.set(e.to_string())}}>"Export GenBank"</button><details class="mvp-file-actions"><summary aria-label="More file actions">"More"</summary><div class="mvp-file-menu"><button class="button secondary" on:click=move|_|{modal_error.set(String::new());modal.set("document");}>"Document settings"</button><button class="button secondary" on:click=move|_|{let mol=molecule.get_untracked();match formats::export_with_report(&mol,formats::TextFormat::Fasta){Ok(report)=>if let Err(error)=download(&format!("{}.fasta",mol.name),"text/plain",report.text.as_bytes()){message.set(error);}else{message.set(report.warnings.join(" "));},Err(e)=>message.set(e.to_string())}}>"FASTA"</button><button class="button secondary" on:click=move|_|match formats::export_with_report(&molecule.get_untracked(),formats::TextFormat::Embl){Ok(report)=>{if let Err(e)=download("sequence.embl","text/plain",report.text.as_bytes()){message.set(e);}else{message.set(report.warnings.join(" "));}},Err(e)=>message.set(e.to_string())}>"EMBL"</button><Show when=move||molecule.get().metadata.contains_key(snapgene::ORIGINAL_BASE64_KEY)><button class="button secondary" on:click=move|_|match snapgene::original_snapgene_bytes(&molecule.get_untracked()){Ok(bytes)=>{if let Err(e)=download("original.dna","application/octet-stream",&bytes){message.set(e);}},Err(e)=>message.set(e.to_string())}>"Original .dna"</button></Show><Show when=move||molecule.get().metadata.contains_key("trace_original_base64")><button class="button secondary" on:click=move|_|match traces::original_trace_bytes(&molecule.get_untracked()){Ok(bytes)=>{let extension=if bytes.starts_with(b"ABIF"){"ab1"}else{"scf"};if let Err(e)=download(&format!("original.{extension}"),"application/octet-stream",&bytes){message.set(e);}},Err(e)=>message.set(e.to_string())}>"Original trace"</button></Show><Show when=move||molecule.with(|doc|doc.metadata.contains_key(collagene_core::legacy_formats::ORIGINAL_BASE64_KEY))><button class="button secondary" on:click=move|_|match collagene_core::legacy_formats::original_legacy_bytes(&molecule.get_untracked()){Ok(bytes)=>{let extension=if collagene_core::legacy_formats::is_collagene_xml(&bytes){"xml"}else{"xdna"};if let Err(error)=download(&format!("original-legacy.{extension}"),"application/octet-stream",&bytes){message.set(error);}},Err(error)=>message.set(error.to_string())}>"Original legacy file"</button></Show><button class="button secondary" on:click=move|_|match svg::render_map(&molecule.get_untracked()){Ok(text)=>{if let Err(e)=download("sequence-map.svg","image/svg+xml",text.as_bytes()){message.set(e);}},Err(e)=>message.set(e.to_string())}>"SVG map"</button></div></details></div>
                        </div>
                        <crate::digest::DigestProvenance molecule=molecule/>
                        <Show when=move||molecule.with(|m|m.metadata.contains_key("crispr:producer_json"))><div class="crispr-historical"><strong>"CRISPR producer provenance · unverified on opening"</strong><p>"This saved sequence design does not restore a current source/guide/target check or prove cellular editing. Recompute against the exact original source and independent target before reuse."</p><button class="button secondary" on:click=move|_|molecule.with_untracked(|m|{if let Some(json)=m.metadata.get("crispr:producer_json"){if let Err(e)=download("crispr-retained-provenance.json","application/json",json.as_bytes()){message.set(e);}}})>"Download retained CRISPR provenance"</button></div></Show>
                        <div class="editor-toolbar"><div class="tool-group"><button class="icon-button" title="Show documents" on:click=move|_|show_library.update(|v|*v=!*v)>"☷"</button><button class="icon-button" title="Undo" disabled=move||documents.with(|docs|!docs[active.get()].can_undo()) on:click=move|_|documents.update(|docs|{let d=&mut docs[active.get_untracked()];if let Err(error)=d.undo(){message.set(error.to_string());}})>"↶"</button><button class="icon-button" title="Redo" disabled=move||documents.with(|docs|!docs[active.get()].can_redo()) on:click=move|_|documents.update(|docs|{let d=&mut docs[active.get_untracked()];if let Err(error)=d.redo(){message.set(error.to_string());}})>"↷"</button><span class="tool-separator"></span><button class="text-tool" disabled=move||selection.get().is_none()||inspected_feature.get().is_some_and(|f|f.location.spans().len()>1) on:click=move|_|{edit_text.set(selection.get_untracked().and_then(|(a,b)|molecule.get_untracked().sequence.as_str().get(a..b).map(str::to_owned)).unwrap_or_default());modal_error.set(String::new());modal.set("edit");}>"Edit selection"</button></div>
                            <crate::clipboard::ClipboardTools molecule=molecule selection=selection selected_feature=selected_feature locked=workspace_restoring epoch=workspace_epoch identity=Memo::new(move|_|documents.with(|docs|docs.get(active.get()).map(|doc|doc.id.clone()).unwrap_or_default())) on_commit=commit_annotations message=message/>
                            <crate::search_view::SearchTools molecule=molecule identity=Memo::new(move|_|documents.with(|docs|docs[active.get()].id.clone())) epoch=workspace_epoch selection=selection selected_feature=selected_feature regions=search_regions message=message/>
                        </div>
                        <Show when=move||!message.get().is_empty()><div class="notice" role="status"><span>{move||message.get()}</span><button class="icon-button" aria-label="Dismiss message" on:click=move|_|message.set(String::new())>"×"</button></div></Show>
                        <Show when=move||molecule.with(|m|m.metadata.get("collagene:example").is_some_and(|value|value=="interface-demonstration"))><div class="demo-note mvp-example"><div><strong>"Explore the example, or open your sequence"</strong><span>"Synthetic annotations illustrate the viewer. Open a GenBank file to inspect your own sequence and simulate a digest."</span></div><button class="button secondary" on:click=move|_|if let Some(input)=file_input.get(){input.click();}>"Open a GenBank file"</button></div></Show>
                        <div class="editor-view-bar"><div class="view-switcher" role="group" aria-label="Editor view">
                            <button aria-label="Split view" aria-pressed=move||(editor_view.get()=="split").to_string() class:active=move||editor_view.get()=="split" on:click=move|_|editor_view.set("split")>"Split"</button>
                            <button aria-label="Sequence view" aria-pressed=move||(editor_view.get()=="sequence").to_string() class:active=move||editor_view.get()=="sequence" on:click=move|_|editor_view.set("sequence")>"Sequence"</button>
                            <button aria-label="Map view" aria-pressed=move||(editor_view.get()=="map").to_string() class:active=move||editor_view.get()=="map" on:click=move|_|editor_view.set("map")>"Map"</button>
                        </div></div>
                        <div class="editor-grid" class:sequence-only=move||editor_view.get()=="sequence" class:map-only=move||editor_view.get()=="map">
                            <section class="sequence-panel"><div class="panel-heading"><div><h2>"Sequence"</h2></div><span class="panel-subtitle">{move||molecule.with(|doc|if doc.sequence.alphabet()==Alphabet::Protein{"Amino-acid reference"}else{"Reference & complement"})}</span></div><Show when=move||molecule.with(|doc|doc.metadata.contains_key("legacy:upper_aligned"))><details class="legacy-source-rows"><summary>"Original aligned source rows"</summary><p class="helper">"Upper is stored 5′→3′; lower is aligned 3′→5′. These original rows may describe an earlier document revision. The editor displays the reference and its calculated complement."</p><pre>{move||molecule.with(|doc|format!("{}\n{}",doc.metadata.get("legacy:upper_aligned").map(|value|value.chars().take(2000).collect::<String>()).unwrap_or_default(),doc.metadata.get("legacy:lower_aligned").map(|value|value.chars().take(2000).collect::<String>()).unwrap_or_default()))}</pre><p class="helper">"First 2,000 columns shown; original legacy-file download retains all source evidence."</p></details></Show><Show when=move||molecule.get().metadata.contains_key("trace_original_base64")><button class="text-tool" on:click=move|_|trace_view.update(|v|*v=!*v)>{move||if trace_view.get(){"Show sequence"}else{"Show chromatogram"}}</button></Show><Show when=move||trace_view.get()&&molecule.get().metadata.contains_key("trace_original_base64") fallback=move||view!{<crate::sequence_view::SequenceView molecule=molecule selection=selection selected_feature=selected_feature highlighted_spans=highlighted_spans search_regions=search_regions on_show_features=Callback::new(move|()|panel.set("features")) identity=Memo::new(move|_|documents.with(|docs|format!("{}:{}",workspace_epoch.get(),docs[active.get()].id)))/>}><crate::trace_view::TraceView molecule=molecule selection=selection/></Show></section>
                            <section class="map-panel"><div class="panel-heading"><div><h2>{move||if molecule.get().topology==Topology::Circular{"Plasmid map"}else{"Linear map"}}</h2></div><div class="map-controls"><button class="icon-button" title="Zoom out" on:click=move|_|if let Some(r)=renderer.get_value(){r.zoom(1./1.15);}>"−"</button><button class="icon-button" title="Reset map view" on:click=move|_|if let Some(r)=renderer.get_value(){r.reset();}>"Fit"</button><button class="icon-button" title="Zoom in" on:click=move|_|if let Some(r)=renderer.get_value(){r.zoom(1.15);}>"+"</button></div></div>
                                <div class="map-stage"><canvas node_ref=gl_canvas aria-label=move||if molecule.with(|m|m.topology==Topology::Circular){"Interactive circular sequence overview; drag to rotate, use zoom buttons, click to select"}else{"Interactive linear sequence overview; use zoom buttons, click to select"} on:pointerdown=move|e|{if let Some(canvas)=gl_canvas.get(){let _=canvas.set_pointer_capture(e.pointer_id());}if let Some(r)=renderer.get_value(){r.pointer_down(e.client_x()as f64,e.client_y()as f64);}} on:pointermove=move|e|if let Some(r)=renderer.get_value(){r.pointer_move(e.client_x()as f64,e.client_y()as f64);} on:pointerup=move|e|if let Some(r)=renderer.get_value(){if let Some(hit)=r.pointer_up(e.client_x()as f64,e.client_y()as f64){batch(||{search_regions.set(None);selected_feature.set(hit.feature_index);selection.set(Some(hit.span));if hit.feature_index.is_some(){panel.set("features");}});}} on:wheel=move|e|{e.prevent_default();if let Some(r)=renderer.get_value(){r.zoom(if e.delta_y()>0.{0.95}else{1.05});}}></canvas><canvas node_ref=label_canvas class="map-labels" aria-hidden="true"></canvas><p class="map-geometry-note" role="status">"Overview is simplified. Browse the annotation pages or export the full annotation set."</p>
                                    <Show when=move||!render_error.get().is_empty()><div class="map-error">{move||render_error.get()}</div></Show>
                                </div><div class="map-caption">{move||if molecule.get().topology==Topology::Circular{"Drag to rotate · scroll to zoom"}else{"Scroll to zoom · click to select"}}</div>
                                <Show when=move||{let stats=map_labels.get();stats.omitted>0||stats.abbreviated>0||stats.limited}><div class="map-label-note"><span>{move||{let stats=map_labels.get();let mut parts=vec![format!("{} of {} annotation labels shown",stats.displayed,stats.total)];if stats.abbreviated>0{parts.push(format!("{} abbreviated",stats.abbreviated));}if stats.unplaced>0{parts.push(format!("{} without local coordinates",stats.unplaced));}if stats.limited{parts.push("display limit reached".into());}parts.join(" · ")}}</span><button class="text-tool" on:click=move|_|panel.set("features")>"Inspect annotations"</button></div></Show>
                                <div class="legend"><span><i style="background:#2b8775"></i>"Coding sequence"</span><span><i style="background:#d6933b"></i>"Promoter"</span><span><i style="background:#666bb3"></i>"Origin"</span></div>
                            </section>
                            <aside class="inspector"><crate::desktop_shell::InspectorNavigation panel=panel on_auxiliary=Callback::new(move|action|{use crate::desktop_shell::AuxiliaryAction;modal_error.set(String::new());match action{AuxiliaryAction::ImportUrl=>modal.set("url"),AuxiliaryAction::Alignments=>{pending_alignment.set(None);modal.set("alignments");},AuxiliaryAction::Oligos=>{existing_oligo.set(None);pending_oligo.set(None);modal.set("oligos");}}})/>
                                <Show when=move||panel.get()=="library"><div class="inspector-content"><crate::feature_tools::FeatureLibraryPanel molecule=molecule selected_feature=selected_feature libraries=libraries on_keep=keep_library on_commit=commit_annotations selection=selection/></div></Show>
                                <Show when=move||panel.get()=="annotation_files"><div class="inspector-content"><crate::feature_tools::AnnotationFilesPanel molecule=molecule source=annotation_source format=annotation_format on_commit=commit_annotations selection=selection selected_feature=selected_feature/></div></Show>
                                <Show when=move||panel.get()=="history"><div class="inspector-content"><crate::history_view::HistoryPanel molecule=molecule events=Memo::new_owning(move|_|(documents.with(|docs|docs[active.get()].summaries()),true)) current_event=Memo::new(move|_|documents.with(|docs|docs[active.get()].head_event().to_string())) navigation=Memo::new(move|_|documents.with(|docs|docs[active.get()].navigation_event_ids())) external_events=Memo::new(move|_|documents.with(|docs|docs.iter().enumerate().filter(|(i,_)|*i!=active.get()).flat_map(|(_,doc)|{let label=format!("[{}] {}",doc.id,doc.molecule.name.chars().take(100).collect::<String>());doc.summaries().into_iter().map(move|e|(e.event_id,label.clone(),e.revision))}).collect())) on_restore_revision=restore_revision on_fork_revision=fork_revision candidates=Memo::new(move|_|{let current=candidate.get().map(|c|c.plan_revision);stored_candidates.with(|values|values.iter().map(|c|(c.plan_revision.clone(),c.plan.name.clone(),c.plan.target_revision.clone(),current.as_ref()==Some(&c.plan_revision))).collect())}) on_open_candidate=Callback::new(move|revision:String|{if workspace_restoring.get_untracked(){return;}if let Some(value)=stored_candidates.with_untracked(|values|values.iter().find(|c|c.plan_revision==revision).cloned()){candidate_events.set(None);candidate.set(Some(value));plan_report.set(None);last_plan.set(None);plan_generation.update(|g|*g=g.wrapping_add(1));panel.set("plan");}}) on_remove_candidate=Callback::new(move|revision:String|{if workspace_restoring.get_untracked(){return;}batch(||{if candidate.with_untracked(|c|c.as_ref().is_some_and(|c|c.plan_revision==revision)){candidate_events.set(None);candidate.set(None);plan_report.set(None);last_plan.set(None);plan_generation.update(|g|*g=g.wrapping_add(1));}stored_candidates.update(|values|values.retain(|c|c.plan_revision!=revision));});}) collections=Memo::new(move|_|oligo_tables.with(|tables|tables.iter().enumerate().map(|(i,t)|(i,format!("Oligo collection {}",t.source_blake3.chars().take(12).collect::<String>()),t.records.len())).collect())) on_open_oligos=Callback::new(move|index:usize|{if let Some(table)=oligo_tables.with_untracked(|tables|tables.get(index).cloned()){existing_oligo.set(Some(table.as_ref().clone()));pending_oligo.set(None);modal_error.set(String::new());modal.set("oligos");}})/><Show when=move||!recorded_runs.with(|runs|runs.is_empty())><section class="recorded-runs"><h3>"Recorded reconstructions"</h3><p class="helper">"Saved computational history. Each record needs a fresh replay before its verification can be used."</p>{move||recorded_runs.with(|runs|runs.iter().skip(run_page.get()*20).take(20).map(|run|{let id=run.run_id.clone();view!{<article class="run-product"><strong>{run.plan_name.clone()}</strong><p class="helper">{format!("{} products · recorded {:?} · {} · unverified",run.products.len(),run.recorded_status,run.recorded_engine)}</p><button class="button secondary" on:click=move|_|{selected_run.set(id.clone());modal_error.set(String::new());modal.set("run_details");}>"View cloning run"</button></article>}}).collect_view())}<div class="pagination"><button class="button secondary" disabled=move||run_page.get()==0 on:click=move|_|run_page.update(|p|*p=p.saturating_sub(1))>"Previous records"</button><button class="button secondary" disabled=move||{(run_page.get()+1)*20>=recorded_runs.with(|r|r.len())} on:click=move|_|run_page.update(|p|*p+=1)>"Next records"</button></div></section></Show><Show when=move||!saved_assemblies.with(|rows|rows.is_empty())><section class="saved-assemblies"><h3>"Saved read assemblies"</h3><p class="helper">"Historical evidence only; replay exact retained reads before using a consensus."</p>{move||saved_assemblies.with(|rows|rows.iter().skip(assembly_page.get()*20).take(20).map(|row|{let hash=row.attachment_hash.clone();view!{<article><p>{format!("{} reads · {} contigs · recorded complete {} · {}",row.read_count,row.contig_count,row.recorded_complete,row.engine_version)}</p><button class="button secondary" on:click=move|_|{inspected_assembly.set(hash.clone());modal_error.set(String::new());modal.set("saved_contig");}>"View saved read assembly"</button></article>}}).collect_view())}<button class="button secondary" disabled=move||assembly_page.get()==0 on:click=move|_|assembly_page.update(|p|*p=p.saturating_sub(1))>"Previous saved assemblies"</button><button class="button secondary" disabled=move||{(assembly_page.get()+1)*20>=saved_assemblies.with(Vec::len)} on:click=move|_|assembly_page.update(|p|*p+=1)>"Next saved assemblies"</button></section></Show></div></Show>
                                <Show when=move||panel.get()=="reads"><div class="inspector-content"><h2>"Read validation"</h2><Show when=move||molecule.with(|m|m.metadata.contains_key(oligo_formats::FASTQ_ORIGINAL_KEY))><section><h3>"FASTQ quality evidence"</h3><p class="helper">{move||match oligo_formats::quality_evidence(&molecule.get()){Ok(e)=>format!("{} calls · {:?} · {}. {}",e.qualities.len(),e.encoding,if e.current_valid{"Quality mapping matches current calls"}else{"Current calls changed: quality mapping is stale"},e.warnings.join(" ")),Err(e)=>e.to_string()}}</p><label>"FASTQ output encoding"<select aria-label="FASTQ output encoding" prop:value=move||match fastq_encoding.get(){None=>"",Some(oligo_formats::QualityEncoding::Phred33)=>"33",Some(oligo_formats::QualityEncoding::Phred64)=>"64"} on:change=move|e|fastq_encoding.set(match event_target_value(&e).as_str(){"33"=>Some(oligo_formats::QualityEncoding::Phred33),"64"=>Some(oligo_formats::QualityEncoding::Phred64),_=>None})><option value="">"Choose encoding"</option><option value="33">"Phred+33"</option><option value="64">"Phred+64"</option></select></label><button class="button secondary full" disabled=move||fastq_encoding.get().is_none() on:click=move|_|{let Some(encoding)=fastq_encoding.get_untracked()else{return;};match oligo_formats::export_fastq(&molecule.get_untracked(),encoding){Ok(report)=>match download("current.fastq","text/plain",report.text.as_bytes()){Ok(())=>message.set(report.warnings.join(" ")),Err(e)=>message.set(e)},Err(e)=>message.set(e.to_string())}}>"Export current FASTQ"</button><button class="button secondary full" on:click=move|_|match oligo_formats::original_fastq_bytes(&molecule.get_untracked()){Ok(bytes)=>if let Err(e)=download("original.fastq","text/plain",&bytes){message.set(e)},Err(e)=>message.set(e.to_string())}>"Original FASTQ record"</button><button class="button secondary full" disabled=move||selection.get().is_none() on:click=move|_|{let Some((a,b))=selection.get_untracked()else{return;};match oligo_formats::trim_fastq(&molecule.get_untracked(),a,b){Ok(doc)=>if let Err(e)=add_documents(documents,active,workspace_restoring,vec![doc]){message.set(e)},Err(e)=>message.set(e.to_string())}}>"Open selected FASTQ calls"</button></section></Show><crate::trace_validation::TraceValidationPanel molecule=molecule query_options=Memo::new(move|_|documents.with(|docs|docs.iter().enumerate().map(|(i,d)|(i,d.molecule.name.clone())).collect())) get_query=Callback::new(move|i:usize|documents.with(|docs|docs.get(i).map(|d|d.molecule.clone()))) selection=selection selected_feature=selected_feature on_consensus=Callback::new(move|doc:Molecule|add_documents(documents,active,workspace_restoring,vec![doc]))/></div></Show>
                                <Show when=move||panel.get()=="features"><div class="inspector-content"><h2>"Annotations"<span class="count">{move||molecule.get().features.len()}</span></h2><p class="helper">"Select an annotation to inspect all of its parts and biological sequence."</p><button class="button secondary full" on:click=move|_|{feature_for_editor.set(None);modal.set("feature");}>"New annotation"</button>{move||inspected_feature.get().map(|f|{let location=editing::location_text(&f.location).unwrap_or_else(|_|"unresolved".into());let size=f.location.spans().iter().try_fold(0usize,|n,(a,b,_)|n.checked_add(b-a)).unwrap_or(usize::MAX);let exact=if size<=10_000 {f.location.extract(&molecule.get().sequence).unwrap_or_else(|e|e.to_string())}else{"Sequence exceeds this inspector's 10,000-symbol display limit; use a region/sequence export.".into()};view!{<div class="annotation-details"><h3>{f.label().to_string()}</h3><button class="text-tool" on:click=move|_|panel.set("library")>"Use selected annotation in library"</button><p>{location}</p><pre class="annotation-sequence">{exact}</pre><p class="helper">"Compound selections highlight every span. Choose contiguous bases to edit the underlying sequence."</p><button class="button secondary" on:click=move|_|{feature_for_editor.set(selected_feature.get_untracked());modal.set("feature");}>"Edit annotation"</button><button class="text-tool" on:click=move|_|if let Some(i)=selected_feature.get_untracked(){let current=molecule.get_untracked();match editing::delete_feature(&current,i){Ok(edited)=>if let Err(e)=commit_document(documents,active,workspace_restoring,&current.revision(),edited){message.set(e);},Err(e)=>message.set(e.to_string())}}>"Delete annotation"</button><details><summary>"Annotation qualifiers"</summary><dl>{f.qualifiers.into_iter().take(200).map(|(key,value)|view!{<dt>{key}</dt><dd>{value.unwrap_or_else(||"flag".into())}</dd>}).collect_view()}</dl></details></div>}})}<div class="feature-list">{move||molecule.get().features.into_iter().enumerate().skip(feature_page.get()*500).take(500).map(|(index,f)|{let spans=f.location.spans();let span=spans.first().map(|(a,b,_)|(*a,*b));let label=f.label().to_owned();let color=css_color(&f);let coordinates=spans.iter().map(|(a,b,r)|format!("{}–{}{}",a+1,b,if *r{" reverse"}else{""})).collect::<Vec<_>>().join(", ");view!{<button class="feature-entry" class:active=move||selected_feature.get()==Some(index) on:click=move|_|{selected_feature.set(Some(index));selection.set(span);}><i style=format!("background:{color}")></i><span><strong>{label}</strong><small>{coordinates}</small><small>{format!("{}{}",f.kind,if f.needs_review{" · review required"}else{""})}</small></span></button>}}).collect_view()}</div><Show when=move||{molecule.get().features.len()>500}><div class="annotation-pages"><p class="helper">{move||molecule.with(|m|format!("Annotations {}–{} of {}",feature_page.get()*500+1,((feature_page.get()+1)*500).min(m.features.len()),m.features.len()))}</p><button class="button secondary" disabled=move||feature_page.get()==0 on:click=move|_|{selected_feature.set(None);feature_page.update(|page|*page=page.saturating_sub(1));}>"Previous annotations"</button><button class="button secondary" disabled=move||molecule.with(|m|(feature_page.get()+1)*500>=m.features.len()) on:click=move|_|{selected_feature.set(None);feature_page.update(|page|*page+=1);}>"Next annotations"</button></div></Show></div></Show>
                                <Show when=move||panel.get()=="enzymes"><crate::digest::DigestPanel source=Memo::new(move|_|documents.with(|docs|{let d=&docs[active.get()];crate::digest::SourceInfo{binding:collagene_core::run_lineage::EventBinding{document_id:d.id.clone(),event_id:d.head_event().to_string(),revision:molecule_revision.get()},name:d.molecule.name.clone(),length:d.molecule.sequence.len(),alphabet:d.molecule.sequence.alphabet(),topology:d.molecule.topology}})) get_molecule=Callback::new(move|()|molecule.get_untracked()) epoch=workspace_epoch on_navigate=navigate_digest on_open=open_digest on_methylation=Callback::new(move|()|panel.set("methylation"))/></Show>
                                <Show when=move||panel.get()=="pcr"><div class="inspector-content"><h2>"Simulate PCR"</h2><p class="helper">"Enter both primers 5′→3′. Annealing must match exactly; tails are introduced into the product."</p><label>"Forward annealing"<input prop:value=move||pcr_forward.get() on:input=move|e|pcr_forward.set(event_target_value(&e)) placeholder="ACGT…"/></label><label>"Forward 5′ tail"<input prop:value=move||pcr_forward_tail.get() on:input=move|e|pcr_forward_tail.set(event_target_value(&e)) placeholder="Optional"/></label><label>"Reverse annealing"<input prop:value=move||pcr_reverse.get() on:input=move|e|pcr_reverse.set(event_target_value(&e)) placeholder="ACGT…"/></label><label>"Reverse 5′ tail"<input prop:value=move||pcr_reverse_tail.get() on:input=move|e|pcr_reverse_tail.set(event_target_value(&e)) placeholder="Optional"/></label><button class="button primary full" on:click=do_pcr>"Compute products"</button><p class="helper">"Products open as new documents. This checks sequence geometry, not laboratory success."</p></div></Show>
                                <Show when=move||panel.get()=="plan"><div class="inspector-content"><h2>"Cloning plan"</h2><div class="plan-exchange"><button class="button secondary full" on:click=export_workspace>"Export for AI"</button><button class="button secondary full" on:click=move|_|if let Some(input)=file_input.get(){input.click();}>"Open candidate"</button></div><p class="helper">"Export the active target and open inputs for your MCP client. Open its returned candidate here for review."</p><Show when=move||candidate.get().is_some()><div class="candidate-review"><h3>"Candidate · review before use"</h3>{move||candidate.get().map(|c|view!{<p><strong>{c.plan.name.clone()}</strong></p><p class="helper">{format!("Target revision {}…",&c.plan.target_revision[..12])}</p><h4>"Bound inputs"</h4><ul>{c.plan.inputs.iter().map(|i|view!{<li>{format!("{}: {} · {} bases · {}…",i.id,i.molecule.name,i.molecule.sequence.len(),&i.revision[..12])}</li>}).collect_view()}</ul><h4>"Ordered steps"</h4><ol>{c.plan.steps.iter().map(|step|{let description=match step{Operation::RecombineReference{request,..}=>format!("Reference-only {:?} exchange",request.reaction),Operation::DeclareSourceModifications{..}=>"Declare original source modifications".into(),Operation::CheckedDigest{enzymes,..}=>format!("Checked digest with {}",enzymes.join(", ")),Operation::Tail3Prime{base,..}=>format!("Add {} to both 3′ ends",base.as_str()),Operation::ProfiledOverlap{options,..}=>format!("Profiled {:?} overlap assembly",options.method),Operation::RepairEnds{method,..}=>format!("End repair: {:?}",method),Operation::AnnealOligos{candidate_index,..}=>format!("Anneal oligos: zero-based candidate {}",candidate_index),Operation::PcrWithMismatches{..}=>"Intentional-mismatch PCR".into(),Operation::GoldenGate{options,..}=>format!("Golden Gate with {:?}",options.enzyme),Operation::Digest{enzymes,..}=>format!("Digest with {}",enzymes.join(", ")),Operation::TreatEnds{treatment,..}=>format!("Terminal {treatment:?}"),Operation::Ligate{fragments,circular,..}=>format!("Ligation of {} fragments → {}",fragments.len(),if *circular {"circular"}else{"linear"}),Operation::Pcr{forward,reverse,..}=>format!("PCR with {} and {}",forward.name,reverse.name),Operation::Overlap{fragments,circular,..}=>format!("Overlap assembly of {} fragments → {}",fragments.len(),if *circular {"circular"} else {"linear"})};let details=describe_step(step);view!{<li><details><summary>{description}</summary><pre>{details}</pre></details></li>}}).collect_view()}</ol>})}<button class="button primary full" on:click=verify_candidate>"Verify candidate"</button><button class="button secondary full" on:click=move|_|{candidate_events.set(None);candidate.set(None);plan_report.set(None);last_plan.set(None);plan_generation.update(|g|*g=g.wrapping_add(1));}>"Discard candidate"</button></div></Show><Show when=move||candidate.get().is_none()><h3>"Manual cloning plan"</h3><p class="helper">"The active document is the intended target. Choose source documents and a cloning method; the core reconstructs the ordered candidate."</p><label>"First fragment"<select aria-label="First source document" prop:value=move||first_input.get().to_string() on:change=move|e|first_input.set(event_target_value(&e).parse().unwrap_or(0))>{move||documents.get().iter().enumerate().map(|(i,d)|view!{<option value=i.to_string()>{d.molecule.name.clone()}</option>}).collect_view()}</select></label><Show when=move||plan_method.get()!="restriction"||(restriction_final.get()=="ligate"&&second_restriction.get())><label>"Second fragment"<select aria-label="Second source document" prop:value=move||second_input.get().to_string() on:change=move|e|second_input.set(event_target_value(&e).parse().unwrap_or(0))>{move||documents.get().iter().enumerate().map(|(i,d)|view!{<option value=i.to_string()>{d.molecule.name.clone()}</option>}).collect_view()}</select></label></Show><label>"Method"<select aria-label="Cloning method" prop:value=move||plan_method.get() on:change=move|e|plan_method.set(if event_target_value(&e)=="restriction"{"restriction"}else{"overlap"})><option value="overlap">"Overlap assembly"</option><option value="restriction">"Restriction and ligation"</option></select></label><Show when=move||plan_method.get()=="restriction"><div class="checked-restriction-options"><label>"Restriction model"<select aria-label="Restriction model" prop:value=move||restriction_policy.get() on:change=move|e|restriction_policy.set(if event_target_value(&e)=="checked"{"checked"}else{"legacy"})><option value="legacy">"Assume complete cleavage (legacy)"</option><option value="checked">"Check inherited modification state"</option></select></label><label>"Final restriction product"<select aria-label="Restriction final product" prop:value=move||restriction_final.get() on:change=move|e|restriction_final.set(if event_target_value(&e)=="fragment"{"fragment"}else{"ligate"})><option value="ligate">"Ligate selected substrates"</option><option value="fragment">"Single processed fragment (no ligation)"</option></select></label><p class="helper">"Single-fragment mode keeps the selected digest product, or its explicitly repaired/treated descendant, as the final product. A blocked circle can remain an uncut circle without claiming a ligation."</p><Show when=move||restriction_final.get()=="ligate"><label class="checkbox-label"><input type="checkbox" aria-label="Use second restriction substrate" prop:checked=move||second_restriction.get() on:change=move|e|second_restriction.set(event_target_checked(&e))/>"Use a second restriction substrate"</label></Show><Show when=move||restriction_policy.get()=="checked"><p class="helper">"Declare the original source generation under a closed Dam/Dcm/CpG model. Original eligible positions inherit the declared state; positions outside these original motifs and newly synthesized bases are unmodified under this scope. Descendant motifs are not rescanned. Other modification systems and experimental conditions are not modeled. Unknown or impaired checked cleavage withholds undetermined products; blocked sites stay uncut."</p><details class="checked-profile-identities"><summary>"Pinned positional profiles and formulation identities"</summary><p>{format!("{} · revision {}",modified_duplex::PROFILE_VERSION,modified_duplex::profile_revision())}</p>{methylation::profiles().into_iter().map(|p|view!{<p>{format!("{}: {}",p.enzyme,p.supplier_identity)}</p>}).collect_view()}<p class="helper">"This positional profile is distinct from the source-context methylation review. Generic and HF/v2 enzyme formulations are not substituted. Imported summaries cannot authorize checked state."</p></details></Show>
    <RestrictionFields side="First" enzymes=first_enzymes fragment=first_fragment treatment=first_treatment repair=first_repair reverse=first_reverse/><Show when=move||restriction_policy.get()=="checked"><p class="checked-source-binding">{move||documents.with(|docs|docs.get(first_input.get()).map(|d|format!("First source: document {} · event {}",d.id,d.head_event())).unwrap_or_else(||"Choose a first source".into()))}</p><crate::checked_cloning::SourceAssumptions side="First" signals=first_modifications/></Show>
    <Show when=move||second_restriction.get()&&restriction_final.get()=="ligate"><RestrictionFields side="Second" enzymes=second_enzymes fragment=second_fragment treatment=second_treatment repair=second_repair reverse=second_reverse/><Show when=move||restriction_policy.get()=="checked"><p class="checked-source-binding">{move||documents.with(|docs|docs.get(second_input.get()).map(|d|format!("Second source: document {} · event {}",d.id,d.head_event())).unwrap_or_else(||"Choose a second source".into()))}</p><crate::checked_cloning::SourceAssumptions side="Second" signals=second_modifications/></Show></Show><Show when=move||restriction_final.get()=="ligate"><label class="checkbox-label"><input type="checkbox" aria-label="Allow restriction ligation nicks" prop:checked=move||allow_nicked.get() on:change=move|e|allow_nicked.set(event_target::<HtmlInputElement>(&e).checked())/>"Allow explicitly reported nicks"</label></Show><p class="helper">"Fragment numbers follow ascending reference-strand cut order, not size. Intact circles cannot ligate. The active document remains the independent target; declared source state is not measured modification evidence."</p></div></Show><Show when=move||plan_method.get()=="overlap"><label>"Junction overlap (bp)"<input type="number" min="1" prop:value=move||overlap.get() on:input=move|e|overlap.set(event_target_value(&e))/></label><Show when=move||molecule.get().topology==Topology::Circular><label>"Circular closure overlap (bp)"<input type="number" min="1" prop:value=move||closure_overlap.get() on:input=move|e|closure_overlap.set(event_target_value(&e))/></label></Show></Show><button class="button primary full" on:click=do_plan>"Reconstruct and verify"</button></Show>{move||plan_report.get().map(|report|view!{<div class="verification-result"><strong>{format!("{:?}",report.status)}</strong>{report.findings.into_iter().map(|f|view!{<p>{f}</p>}).collect_view()}{report.products.iter().enumerate().filter_map(|(i,p)|p.metadata.get("duplex_model").and_then(|s|serde_json::from_str::<chemistry::DuplexSummary>(s).ok()).map(|summary|view!{<details class="physical-evidence"><summary>{format!("Step {} · {}",i+1,p.name)}</summary><pre>{summary.describe()}</pre></details>})).collect_view()}{report.products.iter().enumerate().filter_map(|(i,p)|crate::checked_cloning::modification_summary(p).map(|summary|view!{<details class="checked-modification-summary"><summary>{format!("Step {} · positional modification evidence",i+1)}</summary><p>{summary}</p></details>})).collect_view()}<details><summary>"Assumptions and limits"</summary><ul>{report.assumptions.into_iter().map(|a|view!{<li>{a}</li>}).collect_view()}</ul></details><p class="helper">"Alternative assemblies and experimental feasibility are not exhaustively checked."</p></div>})}<button class="button secondary full" disabled=move||plan_report.get().is_none() on:click=move|_|if let Some(report)=plan_report.get_untracked(){if let Err(error)=crate::checked_cloning::export_report(&report){message.set(error);}}>"Download reconstruction report"</button><button class="button secondary full" disabled=move||candidate.get().is_none()&&last_plan.get().is_none() on:click=move|_|{let value=candidate.get_untracked().or_else(||last_plan.get_untracked().map(|(plan,_)|transfer::Candidate::new(plan)));if let Some(value)=value {match transfer::save_candidate(&value){Ok(bytes)=>{if let Err(e)=download("candidate.cg2plan","application/octet-stream",&bytes){message.set(e);}},Err(e)=>message.set(e.to_string())}}}>"Save candidate"</button><button class="button secondary full" disabled=move||last_plan.get().is_none() on:click=move|_|{modal_error.set(String::new());modal.set("record_run");}>"Keep reconstructed products"</button><button class="button secondary full" disabled=move||last_plan.get().is_none() on:click=move|_|if let Some((plan,target))=last_plan.get_untracked(){let html=verification::export_protocol_html(&plan,&target);if let Err(e)=download("cloning-protocol.html","text/html",html.as_bytes()){message.set(e);}}>"Export HTML protocol"</button><button class="button secondary full" disabled=move||pdf_busy.get()||last_plan.get().is_none() on:click=export_pdf>"Export PDF protocol"</button></div></Show>
                                <Show when=move||panel.get()=="primers"><crate::primers::PrimerPanel molecule=molecule selection=selection selected_feature=selected_feature on_use=Callback::new(move|(forward,reverse):(Primer,Primer)|{pcr_forward.set(forward.annealing.as_str().into());pcr_forward_tail.set(forward.tail);pcr_reverse.set(reverse.annealing.as_str().into());pcr_reverse_tail.set(reverse.tail);panel.set("pcr");message.set("Candidate primers loaded. Compute PCR products to inspect alternative binding geometries.".into());}) on_annotate=Callback::new(move|(base,feature):(String,Feature)|{let current=molecule.get_untracked();if current.revision()!=base{return Err("The document changed; rerun primer analysis before annotating it.".into());}let edited=editing::set_feature(&current,None,feature).map_err(|e|e.to_string())?;commit_document(documents,active,workspace_restoring,&base,edited)})/></Show>
                                <Show when=move||panel.get()=="gel"><crate::gel::GelPanel molecule=molecule query_options=Memo::new(move|_|documents.with(|docs|docs.iter().enumerate().map(|(i,d)|(i,d.molecule.name.clone())).collect())) get_query=Callback::new(move|i:usize|documents.with(|docs|docs.get(i).map(|d|d.molecule.clone())))/></Show>
                                <Show when=move||panel.get()=="methods"><crate::cloning_methods::CloningMethodsPanel molecule=molecule query_options=Memo::new(move|_|documents.with(|docs|docs.iter().enumerate().map(|(i,d)|(i,d.molecule.name.clone())).collect())) get_query=Callback::new(move|i:usize|documents.with(|docs|docs.get(i).map(|d|d.molecule.clone()))) config=method_config epoch=workspace_epoch on_candidate=Callback::new(move|value:transfer::Candidate|{if workspace_restoring.get_untracked(){return;}candidate_events.set(None);candidate.set(Some(value));plan_report.set(None);last_plan.set(None);plan_generation.update(|g|*g=g.wrapping_add(1));panel.set("plan");message.set("Method candidate opened for independent verification against the active target. No experimental success or mixture uniqueness is established.".into());}) on_products=Callback::new(move|products:Vec<Molecule>|add_documents(documents,active,workspace_restoring,products))/></Show>
                                <Show when=move||panel.get()=="ends"><crate::end_repair::EndRepairPanel molecule=molecule query_options=Memo::new(move|_|documents.with(|docs|docs.iter().enumerate().map(|(i,doc)|(i,doc.molecule.name.clone())).collect())) get_query=Callback::new(move|index:usize|documents.with(|docs|docs.get(index).map(|doc|doc.molecule.clone()))) on_candidate=Callback::new(move|value:transfer::Candidate|{if workspace_restoring.get_untracked(){return;}candidate_events.set(None);candidate.set(Some(value));plan_report.set(None);last_plan.set(None);plan_generation.update(|generation|*generation=generation.wrapping_add(1));panel.set("plan");message.set("Physical candidate opened for review. Verify it against the independently selected active target; alternatives and chemistry remain explicit.".into());}) on_products=Callback::new(move|products:Vec<Molecule>|add_documents(documents,active,workspace_restoring,products))/></Show>
                                <Show when=move||panel.get()=="mutagenesis"><crate::mutagenesis::MutagenesisPanel molecule=molecule query_options=Memo::new(move|_|documents.with(|docs|docs.iter().enumerate().map(|(i,d)|(i,d.molecule.name.clone())).collect())) get_query=Callback::new(move|i:usize|documents.with(|docs|docs.get(i).map(|d|d.molecule.clone()))) on_candidate=Callback::new(move|value:transfer::Candidate|{candidate_events.set(None);candidate.set(Some(value));plan_report.set(None);last_plan.set(None);plan_generation.update(|g|*g=g.wrapping_add(1));panel.set("plan");message.set("Candidate opened for review. Verify it against the independently selected target and current input revisions.".into());}) on_products=Callback::new(move|products:Vec<Molecule>|add_documents(documents,active,workspace_restoring,products))/></Show>
                                <Show when=move||panel.get()=="golden_gate"><crate::golden_gate::GoldenGatePanel molecule=molecule query_options=Memo::new(move|_|documents.with(|docs|docs.iter().enumerate().map(|(i,d)|(i,d.molecule.name.clone())).collect())) get_query=Callback::new(move|i:usize|documents.with(|docs|docs.get(i).map(|d|d.molecule.clone()))) on_candidate=Callback::new(move|value:transfer::Candidate|{candidate_events.set(None);candidate.set(Some(value));plan_report.set(None);last_plan.set(None);plan_generation.update(|g|*g=g.wrapping_add(1));panel.set("plan");message.set("Candidate opened for review. Verify it against the independently selected target and current input revisions.".into());})/></Show>
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
                                        candidate_events.set(None);candidate.set(Some(value));candidate_events.set(Some(preferred));plan_report.set(None);last_plan.set(None);last_plan_events.set(None);plan_generation.update(|g|*g=g.wrapping_add(1));panel.set("plan");message.set("Reference-only candidate opened for review. Verify against the independent target before recording source-event ancestry.".into());
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
                                <Show when=move||panel.get()=="biology"><crate::biology::BiologyPanel molecule=molecule selected_feature=selected_feature selection=selection on_annotate=Callback::new(move|(base,feature):(String,Feature)|{let current=molecule.get_untracked();if current.revision()!=base{return Err("The document changed; rerun ORF analysis before annotating it.".into());}let edited=editing::set_feature(&current,None,feature).map_err(|e|e.to_string())?;commit_document(documents,active,workspace_restoring,&base,edited)})/><crate::protein_properties::ProteinPanel molecule=molecule/></Show>
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
                                    if active.get_untracked()==index{apply(false);}else{
                                        active.set(index);
                                        // App clears old selections when the active source changes.
                                        // Install this source-bound selection after those effects settle.
                                        let frame=wasm_bindgen::closure::Closure::once_into_js(move||apply(true));
                                        if let Some(window)=web_sys::window(){let _=window.request_animation_frame(frame.unchecked_ref());}
                                    }
                                    if query{panel.set("features");}
                                })/></Show><Show when=move||panel.get()=="compare"><crate::compare::ComparePanel molecule=molecule query_options=Memo::new(move|_|documents.with(|docs|docs.iter().enumerate().map(|(i,d)|(i,d.molecule.name.clone())).collect())) get_query=Callback::new(move|i:usize|documents.with(|docs|docs.get(i).map(|d|d.molecule.clone()))) selection=selection/></Show>
                            </aside>
                        </div>
                        <footer class="statusbar"><span>{move||save_status.get()}</span><span>{move||selection.get().map(|(a,b)|format!("Selection {}–{} · {} {}",a+1,b,b-a,molecule.with(|m|m.sequence.alphabet().map_unit()))).unwrap_or_else(||molecule.with(|m|if m.sequence.alphabet()==Alphabet::Protein{"Click a residue or feature to select".into()}else{"Click a base or feature to select".into()}))}</span><span>{move||documents.with(|docs|format!("{} {}",docs.len(),if docs.len()==1{"document"}else{"documents"}))}</span></footer>
                        </div>
                    </main>
                </div>
                <Show when=move||modal.get()!=""><div class="modal-backdrop"><section class="modal" class:alignment-modal=move||modal.get()=="alignments" role="dialog" aria-modal="true" aria-labelledby="dialog-title"><div class="modal-heading"><h2 id="dialog-title">{move||if modal.get()=="digest_product"{"Open virtual digest fragment"}else if modal.get()=="crispr_product"{"Open checked sequence design"}else if modal.get()=="contig_archive"{"Record read assembly"}else if modal.get()=="saved_contig"{"Saved read assembly"}else if modal.get()=="record_run"{"Keep reconstructed products"}else if modal.get()=="run_details"{"Recorded reconstruction"}else if modal.get()=="alignments"{"Alignments"}else if modal.get()=="workspace"{"Workspace"}else if modal.get()=="oligos"{"Oligo table"}else if modal.get()=="fastq"{"Import FASTQ reads"}else if modal.get()=="import"{"Open a sequence"}else if modal.get()=="url"{"Import from a sequence source"}else if modal.get()=="feature"{"Edit annotation"}else if modal.get()=="document"{"Document settings"}else if molecule.with(|m|m.sequence.alphabet()==Alphabet::Protein){"Edit selected residues"}else{"Edit selected bases"}}</h2><button class="icon-button" aria-label="Close dialog" on:click=move|_|modal.set("")>"×"</button></div><Show when=move||!modal_error.get().is_empty()><p class="dialog-error" role="alert">{move||modal_error.get()}</p></Show><Show when=move||modal.get()=="digest_product">{move||pending_digest.get_untracked().map(|request|view!{<crate::digest_product::DigestOpenPanel request=request snapshot=contig_snapshot validate=Callback::new(move|frame:collagene_core::workspace::Workspace|restore_workspace_state(&frame).map(|_|())) on_restore=contig_restore db=db busy=workspace_restoring loaded=loaded save_generation=save_generation epoch=workspace_epoch message=message/>})}</Show><Show when=move||modal.get()=="crispr_product">{move||pending_donor.get_untracked().map(|request|view!{<crate::crispr::DonorOpenPanel request=request snapshot=contig_snapshot validate=Callback::new(move|frame:collagene_core::workspace::Workspace|restore_workspace_state(&frame).map(|_|())) on_restore=contig_restore db=db busy=workspace_restoring loaded=loaded save_generation=save_generation epoch=workspace_epoch message=message/>})}</Show><Show when=move||modal.get()=="contig_archive">{move||pending_assembly.get_untracked().map(|request|view!{<crate::contigs::ContigArchivePanel request=request snapshot=contig_snapshot validate=Callback::new(move|frame:collagene_core::workspace::Workspace|restore_workspace_state(&frame).map(|_|())) on_restore=contig_restore db=db busy=workspace_restoring loaded=loaded save_generation=save_generation epoch=workspace_epoch message=message/>})}</Show><Show when=move||modal.get()=="saved_contig"><crate::contigs::SavedAssemblyPanel attachment_hash=inspected_assembly.get_untracked() snapshot=contig_snapshot epoch=workspace_epoch on_record=record_assembly on_navigate=navigate_read_call/></Show><Show when=move||modal.get()=="record_run">{move||last_plan.get_untracked().map(|(plan,_)|view!{<crate::run_capture::RunCapturePanel plan=plan preferred_bindings=last_plan_events.get_untracked().filter(|(revision,_)|last_plan.get_untracked().is_some_and(|(p,_)|verification::plan_revision(&p)==*revision)).map(|(_,bindings)|bindings).unwrap_or_default() snapshot=Callback::new(move|()|documents.with_untracked(|docs|stored_candidates.with_untracked(|plans|attachments.with_untracked(|assets|candidate.with_untracked(|current|workspace_snapshot(docs,active.get_untracked(),current.as_ref(),plans,assets,&workspace_name.get_untracked())))))) validate=Callback::new(move|frame:collagene_core::workspace::Workspace|restore_workspace_state(&frame).map(|_|())) on_restore=Callback::new(move|frame:collagene_core::workspace::Workspace|match restore_workspace_state(&frame){Ok(state)=>batch(||apply_restored.run(state)),Err(e)=>message.set(e)}) db=db busy=workspace_restoring loaded=loaded save_generation=save_generation epoch=workspace_epoch message=message/>})}</Show><Show when=move||modal.get()=="run_details"><crate::run_capture::RunDetailsPanel snapshot=Callback::new(move|()|documents.with_untracked(|docs|stored_candidates.with_untracked(|plans|attachments.with_untracked(|assets|candidate.with_untracked(|current|workspace_snapshot(docs,active.get_untracked(),current.as_ref(),plans,assets,&workspace_name.get_untracked())))))) run_id=selected_run.get_untracked() epoch=workspace_epoch on_open_product=Callback::new(move|id:String|{if workspace_restoring.get_untracked(){return;}if let Some(index)=documents.with_untracked(|docs|docs.iter().position(|d|d.id==id)){active.set(index);panel.set("history");modal.set("");}})/></Show><Show when=move||modal.get()=="workspace"><crate::workspace_view::WorkspacePanel snapshot=Callback::new(move|()|documents.with_untracked(|docs|stored_candidates.with_untracked(|plans|attachments.with_untracked(|assets|candidate.with_untracked(|current|workspace_snapshot(docs,active.get_untracked(),current.as_ref(),plans,assets,&workspace_name.get_untracked())))))) validate=Callback::new(move|frame:collagene_core::workspace::Workspace|restore_workspace_state(&frame).map(|_|())) on_restore=Callback::new(move|frame:collagene_core::workspace::Workspace|match restore_workspace_state(&frame){Ok(state)=>batch(||apply_restored.run(state)),Err(e)=>message.set(e)}) db=db busy=workspace_restoring loaded=loaded save_generation=save_generation message=message/></Show>
                    <Show when=move||modal.get()=="alignments"><crate::alignment_view::AlignmentPanel alignments=alignments pending=pending_alignment.get_untracked() epoch=workspace_epoch on_keep={let opened_epoch=workspace_epoch.get_untracked();Callback::new(move|item:crate::alignment_view::AlignmentRevision|{if opened_epoch!=workspace_epoch.get_untracked(){return Err("The workspace changed during alignment import. Open the alignment again.".into());}keep_alignment.run(item)})} on_extract={let opened_epoch=workspace_epoch.get_untracked();Callback::new(move|molecules:Vec<Molecule>|{if opened_epoch!=workspace_epoch.get_untracked(){return Err("The workspace changed during row extraction.".into());}add_documents(documents,active,workspace_restoring,molecules)})}/></Show>
                    <Show when=move||modal.get()=="oligos"><crate::oligo_import::OligoImport on_import=import_oligos existing=existing_oligo.get_untracked() initial_text=pending_oligo.get_untracked()/></Show>
                    <Show when=move||modal.get()=="fastq"><p>"Select the quality encoding supplied by the sequencing source. Phred+33 and Phred+64 overlap; it cannot be inferred from printable characters."</p><label>"Supplied FASTQ encoding"<select aria-label="Supplied FASTQ encoding" prop:value=move||match fastq_encoding.get(){None=>"",Some(oligo_formats::QualityEncoding::Phred33)=>"33",Some(oligo_formats::QualityEncoding::Phred64)=>"64"} on:change=move|e|fastq_encoding.set(match event_target_value(&e).as_str(){"33"=>Some(oligo_formats::QualityEncoding::Phred33),"64"=>Some(oligo_formats::QualityEncoding::Phred64),_=>None})><option value="">"Choose supplied encoding"</option><option value="33">"Phred+33"</option><option value="64">"Phred+64"</option></select></label><button class="button primary" disabled=move||fastq_encoding.get().is_none() on:click=move|_|{let(Some(bytes),Some(encoding))=(pending_fastq.get_untracked(),fastq_encoding.get_untracked())else{return;};match oligo_formats::import_fastq(&bytes,encoding){Ok(report)=>match add_documents(documents,active,workspace_restoring,report.molecules){Ok(())=>{message.set(report.warnings.join(" "));pending_fastq.set(None);modal.set("");panel.set("reads");},Err(e)=>modal_error.set(e)},Err(e)=>modal_error.set(e.to_string())}}>"Import FASTQ with selected encoding"</button></Show>
                    <Show when=move||modal.get()=="import"><p>"Paste GenBank, EMBL, FASTA, or plain DNA. No project setup needed."</p><textarea aria-label="Sequence input" rows="12" prop:value=move||import_text.get() on:input=move|e|import_text.set(event_target_value(&e)) placeholder=">my_sequence\nACGT…"></textarea><button class="button primary" on:click=move|_|{let text=import_text.get_untracked();if text.as_bytes().first()==Some(&b'@'){pending_fastq.set(Some(std::sync::Arc::new(text.into_bytes())));fastq_encoding.set(None);modal.set("fastq");return;}match formats::import_text(&text){Ok(report)=>{if let Err(error)=add_documents(documents,active,workspace_restoring,report.molecules){modal_error.set(error);return;}message.set(report.warnings.join(" "));modal.set("");import_text.set(String::new());},Err(e)=>modal_error.set(e.to_string())}}>"Open sequence"</button></Show><Show when=move||modal.get()=="url"><crate::source_import::SourceImport on_import={let epoch=workspace_epoch.get_untracked();Callback::new(move |report:formats::ImportReport|{if epoch!=workspace_epoch.get_untracked(){return Err("The workspace changed during retrieval. Import the source again in this workspace.".into());}add_documents(documents,active,workspace_restoring,report.molecules)?;message.set(report.warnings.join(" "));modal.set("");Ok(())})}/></Show><Show when=move||modal.get()=="feature"><crate::editors::FeatureEditor molecule=molecule.get_untracked() index=feature_for_editor.get_untracked() selection=selection.get_untracked() on_commit=Callback::new(move|(base,edited):(String,Molecule)|{commit_document(documents,active,workspace_restoring,&base,edited)?;modal.set("");Ok(())})/></Show><Show when=move||modal.get()=="document"><crate::editors::DocumentEditor molecule=molecule.get_untracked() selection=selection.get_untracked() on_commit=Callback::new(move|(base,edited):(String,Molecule)|{commit_document(documents,active,workspace_restoring,&base,edited)?;modal.set("");Ok(())})/></Show><Show when=move||modal.get()=="edit"><p>{move||selection.get().map(|(a,b)|format!("Replace {} {}–{}; clear the field to delete. Affected annotations will be marked for review.",molecule.with(|m|m.sequence.alphabet().selected_symbols()),a+1,b)).unwrap_or_else(||format!("Select {} in the sequence view first.",molecule.with(|m|m.sequence.alphabet().selected_symbols())))}</p><textarea aria-label="Replacement sequence" rows="6" prop:value=move||edit_text.get() on:input=move|e|edit_text.set(event_target_value(&e))></textarea><button class="button primary" disabled=move||selection.get().is_none() on:click=move|_|if let Some((a,b))=selection.get_untracked(){match molecule.get_untracked().replace(a,b,&edit_text.get_untracked()){Ok(edited)=>{let base=molecule.get_untracked().revision();if let Err(error)=commit_document(documents,active,workspace_restoring,&base,edited){modal_error.set(error);return;}selection.set(None);plan_report.set(None);last_plan.set(None);modal.set("");},Err(e)=>modal_error.set(e.to_string())}}>"Apply edit"</button></Show></section></div></Show>
            </div>
        }
}
