//! Existing alignments as retained evidence; no alignment engine is run.
#[path = "alignment_index.rs"]
mod index;
use crate::{renderer::ResizeListener, worker};
use collagene_core::{
    alignment_formats::{self, AlignmentDocument, AlignmentExport, AlignmentFormat},
    protocol::{
        request::Request,
        schema::{AnalysisJob, AnalysisResult},
    },
};
use collagene_model::{Alphabet, Molecule};
use leptos::{html, prelude::*};
use std::{
    collections::{BTreeMap, BTreeSet},
    rc::Rc,
    sync::{
        Arc,
        atomic::{AtomicBool, Ordering},
    },
};
use wasm_bindgen::JsCast;
use wasm_bindgen_futures::{JsFuture, spawn_local};
use web_sys::{CanvasRenderingContext2d as Ctx, HtmlInputElement};
pub const ATTACHMENT_KIND: &str = "alignment-document-v1";
#[derive(Clone)]
pub struct AlignmentRevision {
    pub revision: String,
    pub asset_hash: String,
    pub prepared: Arc<index::PreparedAlignment>,
    pub bytes: Arc<Vec<u8>>,
}
impl AlignmentRevision {
    pub fn new(document: AlignmentDocument) -> Result<Self, String> {
        let bytes = serde_json::to_vec(&document).map_err(|e| e.to_string())?;
        Self::prepare(document, bytes)
    }
    pub fn from_bytes(bytes: &[u8]) -> Result<Self, String> {
        let doc = alignment_formats::load_alignment(bytes).map_err(|e| e.to_string())?;
        Self::prepare(doc, bytes.to_vec())
    }
    fn prepare(document: AlignmentDocument, bytes: Vec<u8>) -> Result<Self, String> {
        if bytes.len() > 32 * 1024 * 1024 {
            return Err("Alignment attachment exceeds32MiB".into());
        }
        let canonical = serde_json::to_vec(&document).map_err(|e| e.to_string())?;
        Ok(Self {
            revision: blake3::hash(&canonical).to_hex().to_string(),
            asset_hash: blake3::hash(&bytes).to_hex().to_string(),
            prepared: Arc::new(index::PreparedAlignment::new(document).map_err(|e| e.to_string())?),
            bytes: Arc::new(bytes),
        })
    }
    pub(crate) fn doc(&self) -> &AlignmentDocument {
        self.prepared.document()
    }
}
#[derive(Clone)]
pub struct PendingAlignment {
    pub name: String,
    pub bytes: Arc<Vec<u8>>,
    pub format: Option<AlignmentFormat>,
    pub native: bool,
}
fn format_name(format: AlignmentFormat) -> &'static str {
    match format {
        AlignmentFormat::AlignedFasta => "Aligned FASTA",
        AlignmentFormat::Clustal => "CLUSTAL",
        AlignmentFormat::Stockholm => "Stockholm",
        AlignmentFormat::PhylipSequential => "PHYLIP sequential (strict IDs)",
        AlignmentFormat::PhylipInterleaved => "PHYLIP interleaved (strict IDs)",
        AlignmentFormat::PhylipRelaxedSequential => "PHYLIP sequential (relaxed IDs)",
        AlignmentFormat::PhylipRelaxedInterleaved => "PHYLIP interleaved (relaxed IDs)",
        AlignmentFormat::Msf => "MSF",
        AlignmentFormat::Nexus => "NEXUS",
        AlignmentFormat::Pir => "PIR / NBRF",
        AlignmentFormat::Gde => "GDE flat",
        AlignmentFormat::Selex => "SELEX",
    }
}
fn extension(format: AlignmentFormat) -> &'static str {
    match format {
        AlignmentFormat::AlignedFasta => "fasta",
        AlignmentFormat::Clustal => "aln",
        AlignmentFormat::Stockholm => "sto",
        AlignmentFormat::PhylipSequential
        | AlignmentFormat::PhylipInterleaved
        | AlignmentFormat::PhylipRelaxedSequential
        | AlignmentFormat::PhylipRelaxedInterleaved => "phy",
        AlignmentFormat::Msf => "msf",
        AlignmentFormat::Nexus => "nex",
        AlignmentFormat::Pir => "pir",
        AlignmentFormat::Gde => "gde",
        AlignmentFormat::Selex => "slx",
    }
}
fn parse_format(value: &str) -> AlignmentFormat {
    match value {
        "Clustal" => AlignmentFormat::Clustal,
        "Stockholm" => AlignmentFormat::Stockholm,
        "PhylipSequential" => AlignmentFormat::PhylipSequential,
        "PhylipInterleaved" => AlignmentFormat::PhylipInterleaved,
        "PhylipRelaxedSequential" => AlignmentFormat::PhylipRelaxedSequential,
        "PhylipRelaxedInterleaved" => AlignmentFormat::PhylipRelaxedInterleaved,
        "Msf" => AlignmentFormat::Msf,
        "Nexus" => AlignmentFormat::Nexus,
        "Pir" => AlignmentFormat::Pir,
        "Gde" => AlignmentFormat::Gde,
        "Selex" => AlignmentFormat::Selex,
        _ => AlignmentFormat::AlignedFasta,
    }
}
async fn analyze(
    job: AnalysisJob,
    revision: String,
    id: u32,
    worker: Rc<worker::CachedWorker>,
) -> Result<AnalysisResult, String> {
    let reply = worker
        .ask(Request {
            id,
            revision: revision.clone(),
            job,
        })
        .await;
    worker.close();
    let response = reply?;
    if response.id != id || response.revision != revision {
        return Err("Alignment response has a different input binding".into());
    }
    match response.result {
        AnalysisResult::Error(e) => Err(e),
        result => Ok(result),
    }
}
#[component]
pub fn AlignmentPanel(
    alignments: RwSignal<Vec<AlignmentRevision>>,
    #[prop(optional_no_strip)] pending: Option<PendingAlignment>,
    epoch: RwSignal<u64>,
    on_keep: Callback<AlignmentRevision, Result<String, String>>,
    on_extract: Callback<Vec<Molecule>, Result<(), String>>,
) -> impl IntoView {
    let current = RwSignal::new(alignments.with_untracked(|items| {
        items
            .last()
            .map(|a| a.asset_hash.clone())
            .unwrap_or_default()
    }));
    let source_text = RwSignal::new(String::new());
    let name = RwSignal::new(
        pending
            .as_ref()
            .map(|p| p.name.clone())
            .unwrap_or_else(|| "Alignment".into()),
    );
    let import_format = RwSignal::new(
        pending
            .as_ref()
            .and_then(|p| p.format)
            .unwrap_or(AlignmentFormat::AlignedFasta),
    );
    let alphabet = RwSignal::new(None::<Alphabet>);
    let uploaded = RwSignal::new(pending);
    let file = NodeRef::<html::Input>::new();
    let busy = RwSignal::new(false);
    let status = RwSignal::new(String::new());
    let generation = RwSignal::new(0u32);
    let selected = RwSignal::new(BTreeSet::<String>::new());
    let row_page = RwSignal::new(0usize);
    let annotation_page = RwSignal::new(0usize);
    let focused = RwSignal::new((0usize, 0usize));
    let antibody_report =
        RwSignal::new(None::<(u64, Arc<collagene_core::alignment_antibody::Report>)>);
    let antibody_candidate = RwSignal::new(0usize);
    let export_format = RwSignal::new(AlignmentFormat::AlignedFasta);
    let exported = RwSignal::new(None::<(String, AlignmentFormat, AlignmentExport)>);
    let active_worker = StoredValue::new_local(None::<Rc<worker::CachedWorker>>);
    let start_worker = move || {
        if let Some(old) = active_worker.get_value() {
            old.close();
        }
        let worker = Rc::new(worker::CachedWorker::new()?);
        active_worker.set_value(Some(worker.clone()));
        Ok::<_, String>(worker)
    };
    let alive = Arc::new(AtomicBool::new(true));
    let a = alive.clone();
    on_cleanup(move || {
        a.store(false, Ordering::Relaxed);
        if let Some(worker) = active_worker.get_value() {
            worker.close();
        }
        active_worker.set_value(None);
    });
    let alignment = Memo::new_owning(move |_| {
        (
            alignments.with(|items| {
                items
                    .iter()
                    .find(|a| a.asset_hash == current.get())
                    .cloned()
            }),
            true,
        )
    });
    let antibody = Memo::new(move |_| {
        let (e, r) = antibody_report.get()?;
        let a = alignment.get()?;
        let row = focused.get().0;
        (e == epoch.get()
            && r.source.alignment_revision == a.revision
            && a.doc()
                .rows
                .get(row)
                .is_some_and(|row| row.id == r.source.row_id))
        .then_some(r)
    });
    Effect::new(move |_| {
        let _ = (current.get(), epoch.get());
        if let Some(worker) = active_worker.get_value() {
            worker.close();
        }
        active_worker.set_value(None);
        generation.update(|g| *g = g.wrapping_add(1));
        selected.set(BTreeSet::new());
        row_page.set(0);
        annotation_page.set(0);
        focused.set((0, 0));
        exported.set(None);
        busy.set(false);
    });
    Effect::new(move |_| {
        let _ = (
            import_format.get(),
            alphabet.get(),
            name.get(),
            source_text.get(),
            uploaded.get().as_ref().map(|p| p.name.clone()),
            export_format.get(),
            selected.get(),
        );
        if let Some(worker) = active_worker.get_value() {
            worker.close();
        }
        active_worker.set_value(None);
        generation.update(|g| *g = g.wrapping_add(1));
        exported.set(None);
        busy.set(false);
    });
    let import_alive = alive.clone();
    let import = move || {
        let pending = uploaded.get_untracked();
        let native = pending.as_ref().is_some_and(|p| p.native);
        let bytes = pending
            .as_ref()
            .map(|p| p.bytes.as_ref().clone())
            .unwrap_or_else(|| source_text.get_untracked().into_bytes());
        if bytes.len()
            > if native {
                32 * 1024 * 1024
            } else {
                alignment_formats::MAX_ALIGNMENT_BYTES
            }
        {
            status.set("Alignment exceeds the import size limit".into());
            return;
        }
        if !native && alphabet.get_untracked().is_none() {
            status.set("Choose the supplied alignment alphabet explicitly".into());
            return;
        }
        let revision = blake3::hash(&bytes).to_hex().to_string();
        let job = if native {
            AnalysisJob::LoadAlignment { bytes }
        } else {
            AnalysisJob::ImportAlignment {
                bytes,
                format: import_format.get_untracked(),
                alphabet: alphabet.get_untracked().expect("checked alphabet"),
                name: name.get_untracked(),
            }
        };
        let worker = match start_worker() {
            Ok(w) => w,
            Err(e) => {
                status.set(e);
                return;
            }
        };
        let id = generation.get_untracked().wrapping_add(1);
        generation.set(id);
        busy.set(true);
        status.set("Reading alignment in background…".into());
        let expected_epoch = epoch.get_untracked();
        let alive = import_alive.clone();
        spawn_local(async move {
            let result = analyze(job, revision, id, worker).await;
            if !alive.load(Ordering::Relaxed)
                || epoch.get_untracked() != expected_epoch
                || generation.get_untracked() != id
            {
                return;
            }
            busy.set(false);
            match result {
                Ok(AnalysisResult::AlignmentDocument(doc)) => {
                    match AlignmentRevision::new(doc).and_then(|item| on_keep.run(item)) {
                        Ok(key) => {
                            current.set(key);
                            status.set("Alignment retained in this workspace. Select rows to open ungapped sequence documents.".into());
                        }
                        Err(e) => status.set(e),
                    }
                }
                Ok(_) => status.set("Unexpected alignment import response".into()),
                Err(e) => status.set(e),
            }
        });
    };
    let file_alive = alive.clone();
    let open_file = move |event: web_sys::Event| {
        let input = event_target::<HtmlInputElement>(&event);
        let Some(file) = input.files().and_then(|files| files.get(0)) else {
            return;
        };
        input.set_value("");
        if file.size() > 32. * 1024. * 1024. {
            status.set("Alignment file exceeds32MiB".into());
            return;
        }
        let id = generation.get_untracked().wrapping_add(1);
        generation.set(id);
        let expected_epoch = epoch.get_untracked();
        let alive = file_alive.clone();
        spawn_local(async move {
            let result = JsFuture::from(file.array_buffer()).await;
            if !alive.load(Ordering::Relaxed)
                || epoch.get_untracked() != expected_epoch
                || generation.get_untracked() != id
            {
                return;
            }
            match result {
                Ok(buffer) => {
                    let bytes = js_sys::Uint8Array::new(&buffer).to_vec();
                    let filename = file.name();
                    let ext = filename
                        .rsplit('.')
                        .next()
                        .unwrap_or("")
                        .to_ascii_lowercase();
                    let format = if bytes.starts_with(b"CLUSTAL") {
                        Some(AlignmentFormat::Clustal)
                    } else if bytes.starts_with(b"# STOCKHOLM") {
                        Some(AlignmentFormat::Stockholm)
                    } else if bytes.starts_with(b"!!AA_MULTIPLE_ALIGNMENT")
                        || bytes.starts_with(b"!!NA_MULTIPLE_ALIGNMENT")
                    {
                        Some(AlignmentFormat::Msf)
                    } else if bytes
                        .get(..6)
                        .is_some_and(|header| header.eq_ignore_ascii_case(b"#NEXUS"))
                    {
                        Some(AlignmentFormat::Nexus)
                    } else if matches!(ext.as_str(), "pir" | "nbrf")
                        || [b">P1;", b">F1;", b">DL;", b">RL;", b">DC;", b">RC;"]
                            .iter()
                            .any(|header| bytes.starts_with(*header))
                    {
                        Some(AlignmentFormat::Pir)
                    } else if ext == "gde" {
                        Some(AlignmentFormat::Gde)
                    } else if matches!(ext.as_str(), "slx" | "selex" | "sel") {
                        Some(AlignmentFormat::Selex)
                    } else if matches!(
                        ext.as_str(),
                        "afa" | "afasta" | "fa" | "fasta" | "faa" | "fna"
                    ) {
                        Some(AlignmentFormat::AlignedFasta)
                    } else {
                        None
                    };
                    if let Some(format) = format {
                        import_format.set(format);
                    }
                    name.set(filename.clone());
                    uploaded.set(Some(PendingAlignment {
                        name: filename,
                        bytes: Arc::new(bytes),
                        format,
                        native: ext == "cg2aln",
                    }));
                    status.set(
                        "File ready. Choose its format and supplied alphabet, then import.".into(),
                    );
                }
                Err(_) => status.set("Could not read alignment file".into()),
            }
        });
    };
    let export_alive = alive.clone();
    let prepare_export = Callback::new(move |_: ()| {
        let Some(item) = alignment.get_untracked() else {
            return;
        };
        let source = item.asset_hash.clone();
        let format = export_format.get_untracked();
        let expected_epoch = epoch.get_untracked();
        let worker = match start_worker() {
            Ok(w) => w,
            Err(e) => {
                status.set(e);
                return;
            }
        };
        let id = generation.get_untracked().wrapping_add(1);
        generation.set(id);
        busy.set(true);
        exported.set(None);
        status.set("Preparing alignment export…".into());
        let alive = export_alive.clone();
        spawn_local(async move {
            let result = analyze(
                AnalysisJob::ExportAlignment {
                    document: item.doc().clone(),
                    format,
                },
                item.revision.clone(),
                id,
                worker,
            )
            .await;
            if !alive.load(Ordering::Relaxed)
                || epoch.get_untracked() != expected_epoch
                || generation.get_untracked() != id
                || current.get_untracked() != source
                || export_format.get_untracked() != format
            {
                return;
            }
            busy.set(false);
            match result {
                Ok(AnalysisResult::AlignmentExport(report)) => {
                    exported.set(Some((source.clone(), format, report)));
                    status.set(
                        "Review the export losses, then download the canonical alignment.".into(),
                    );
                }
                Ok(_) => status.set("Unexpected alignment export response".into()),
                Err(e) => status.set(e),
            }
        });
    });
    let extract_alive = alive.clone();
    let extract = Callback::new(move |_: ()| {
        let Some(item) = alignment.get_untracked() else {
            return;
        };
        let ids = selected.get_untracked();
        if ids.is_empty() {
            status.set("Select at least one alignment row".into());
            return;
        }
        let source = item.asset_hash.clone();
        let expected_epoch = epoch.get_untracked();
        let worker = match start_worker() {
            Ok(w) => w,
            Err(e) => {
                status.set(e);
                return;
            }
        };
        let id = generation.get_untracked().wrapping_add(1);
        generation.set(id);
        busy.set(true);
        status.set("Extracting ungapped rows…".into());
        let alive = extract_alive.clone();
        spawn_local(async move {
            let result = analyze(
                AnalysisJob::ExtractAlignmentRows {
                    document: item.doc().clone(),
                    row_ids: ids.iter().cloned().collect(),
                },
                item.revision.clone(),
                id,
                worker,
            )
            .await;
            if !alive.load(Ordering::Relaxed)
                || epoch.get_untracked() != expected_epoch
                || generation.get_untracked() != id
                || current.get_untracked() != source
                || selected.get_untracked() != ids
            {
                return;
            }
            busy.set(false);
            match result {
                Ok(AnalysisResult::AlignmentRows(molecules)) => {
                    let count = molecules.len();
                    match on_extract.run(molecules){Ok(())=>status.set(format!("Opened {count} ungapped sequence documents. Close Alignments to edit them; the original alignment remains retained.")),Err(e)=>status.set(e)}
                }
                Ok(_) => status.set("Unexpected alignment extraction response".into()),
                Err(e) => status.set(e),
            }
        });
    });
    view! {
        <div class="alignment-panel">
            <p class="helper">"Open an existing alignment as a separate collection. Alignment columns are not external reference coordinates. No multiple-alignment engine or confidence score is inferred."</p>
            <details class="alignment-import" open=move||alignment.get().is_none()><summary>"Import an alignment"</summary>
                <input node_ref=file type="file" accept=".aln,.clustal,.sto,.stockholm,.afa,.afasta,.fa,.fasta,.fna,.faa,.phy,.phylip,.msf,.nex,.nexus,.pir,.nbrf,.gde,.slx,.selex,.sel,.cg2aln" style="display:none" on:change=open_file/>
                <button class="button secondary" on:click=move|_|if let Some(input)=file.get(){input.click();}>"Choose alignment file"</button>
                {move||uploaded.get().map(|source|view!{<p class="helper">{format!("File: {} · {} bytes{}",source.name,source.bytes.len(),if source.native{" · saved alignment; stored alphabet retained"}else{""})}</p><button class="text-tool" on:click=move|_|uploaded.set(None)>"Use pasted alignment instead"</button>})}
                <Show when=move||uploaded.get().is_none()><textarea aria-label="Aligned sequence input" rows="6" placeholder=">row_a\nAC-GT\n>row_b\nACCGT" prop:value=move||source_text.get() on:input=move|e|source_text.set(event_target_value(&e))></textarea></Show>
                <div class="alignment-import-fields"><label>"Alignment name"<input aria-label="Alignment name" prop:value=move||name.get() on:input=move|e|name.set(event_target_value(&e))/></label>
                <label>"Source format"<select aria-label="Alignment source format" prop:value=move||format!("{:?}",import_format.get()) disabled=move||uploaded.with(|p|p.as_ref().is_some_and(|p|p.native)) on:change=move|e|import_format.set(parse_format(&event_target_value(&e)))><option value="AlignedFasta">"Aligned FASTA"</option><option value="Clustal">"CLUSTAL"</option><option value="Stockholm">"Stockholm1.0"</option><option value="PhylipSequential">"PHYLIP sequential · strict IDs"</option><option value="PhylipInterleaved">"PHYLIP interleaved · strict IDs"</option><option value="PhylipRelaxedSequential">"PHYLIP sequential · relaxed IDs"</option><option value="PhylipRelaxedInterleaved">"PHYLIP interleaved · relaxed IDs"</option><option value="Msf">"MSF"</option><option value="Nexus">"NEXUS"</option><option value="Pir">"PIR / NBRF"</option><option value="Gde">"GDE flat"</option><option value="Selex">"SELEX"</option></select></label>
                <label>"Supplied alphabet"<select aria-label="Alignment supplied alphabet" prop:value=move||alphabet.get().map(|a|format!("{a:?}")).unwrap_or_default() disabled=move||uploaded.with(|p|p.as_ref().is_some_and(|p|p.native)) on:change=move|e|alphabet.set(match event_target_value(&e).as_str(){"Dna"=>Some(Alphabet::Dna),"Rna"=>Some(Alphabet::Rna),"Protein"=>Some(Alphabet::Protein),_=>None})><option value="">"Choose supplied alphabet"</option><option value="Dna">"DNA"</option><option value="Rna">"RNA"</option><option value="Protein">"Protein"</option></select></label></div>
                <button class="button primary" disabled=move||busy.get() on:click=move|_|import()>"Import alignment"</button>
            </details>
            <p class="alignment-status" role="status">{move||status.get()}</p>
            <Show when=move||!alignments.with(Vec::is_empty)><label class="alignment-collection-label">"Retained alignment"<select aria-label="Retained alignment" prop:value=move||current.get() on:change=move|e|current.set(event_target_value(&e))>{move||alignments.with(|items|items.iter().map(|a|{let key=a.asset_hash.clone();let text=format!("{} · {} rows × {} columns · {:?}",a.doc().name,a.doc().rows.len(),a.doc().column_count(),a.doc().alphabet);view!{<option value=key>{text}</option>}}).collect_view())}</select></label></Show>
            {move||alignment.get().map(|item|{let doc=item.doc();let title=format!("{} · {} rows × {} columns · {:?}",doc.name,doc.rows.len(),doc.column_count(),doc.alphabet);let snapshot=item.clone();let original=item.clone();view!{
                <h3>{title}</h3><p class="helper">{format!("Source: {} · digest {}. Original bytes and imported annotations remain separate from derived sequence documents.",format_name(doc.source.format),&doc.source.blake3[..16])}</p>
                <div class="alignment-actions"><button class="button secondary" on:click=move|_|{crate::app::download_with_completion("alignment.cg2aln","application/json",snapshot.bytes.as_ref(),move|result|{if status.try_get_untracked().is_none(){return;}if let Err(e)=result {status.set(e);}})}>"Save alignment file"</button><button class="button secondary" on:click=move|_|match original.doc().original_bytes(){Ok(bytes)=>{crate::app::download_with_completion(&format!("original-alignment.{}",extension(original.doc().source.format)),"text/plain",&bytes,move|result|{if status.try_get_untracked().is_none(){return;}if let Err(e)=result {status.set(e);}})},Err(e)=>status.set(e.to_string())}>"Original alignment source"</button></div>
                <AlignmentCanvas alignment=item.clone() focused=focused antibody=antibody candidate=antibody_candidate/>
                <Show when=move||alignment.with(|a|a.as_ref().is_some_and(|a|a.doc().alphabet==Alphabet::Protein))>
                    <crate::alignment_antibody_view::AlignmentAntibodyPanel alignment=alignment focused=focused epoch=epoch candidate=antibody_candidate on_report=Callback::new(move|r|antibody_report.set(r))/>
                </Show>
                <p class="helper" aria-live="polite">{move||antibody.get().map(|r|crate::alignment_antibody_model::description(&r,antibody_candidate.get(),focused.get().1)).unwrap_or_default()}</p>
                <div class="alignment-readout" aria-live="polite">{move||alignment.get().map(|a|{let(row,col)=focused.get();let doc=a.doc();let Some(r)=doc.rows.get(row)else{return String::new();};let mut counts=BTreeMap::new();for r in &doc.rows{if let Some(c)=r.aligned.as_bytes().get(col){*counts.entry(*c as char).or_insert(0usize)+=1;}}let position=match a.prepared.position(row,col){Ok(Some(n))=>format!("ungapped row position {}",n+1),Ok(None)=>"gap; no residue position".into(),Err(e)=>e.to_string()};let column_annotations=doc.metadata.annotations.iter().filter(|tag|tag.scope=="GC"||(tag.scope=="GR"&&tag.row_identifier.as_deref()==Some(r.identifier.as_str()))).take(100).filter_map(|tag|tag.value.as_bytes().get(col).map(|symbol|format!("{} {}={}",tag.scope,tag.tag,*symbol as char))).collect::<Vec<_>>().join(", ");let consensus=doc.metadata.clustal_consensus.as_ref().and_then(|value|value.as_bytes().get(col)).map(|c|format!("; imported CLUSTAL consensus: {}",*c as char)).unwrap_or_default();let low=col.saturating_sub(30);let high=(col+31).min(r.aligned.len());format!("Row {} · alignment column {} · {} · {}\nVisible row excerpt (columns {}–{}): {}\nObserved symbols in this column: {}. Counts are not alignment confidence.\nImported column annotations (up to100 tags): {}{}",r.identifier,col+1,r.aligned.as_bytes().get(col).copied().map(char::from).unwrap_or('?'),position,low+1,high,r.aligned.get(low..high).unwrap_or(""),counts.iter().map(|(s,n)|format!("{s}: {n}")).collect::<Vec<_>>().join(", "),column_annotations,consensus)}).unwrap_or_default()}</div>
                <details class="alignment-rows"><summary>"Rows and ungapped sequence extraction"</summary><p class="helper">"Selected rows become new linear documents in their stored5′→3′ / N→C order. Gaps are removed. Features, alignment confidence and external coordinates are not invented. All-gap rows cannot become sequence documents."</p><div class="alignment-actions"><button class="text-tool" on:click=move|_|alignment.with_untracked(|a|if let Some(a)=a{selected.set(a.doc().rows.iter().map(|r|r.id.clone()).collect());})>"Select all alignment rows"</button><button class="text-tool" on:click=move|_|selected.set(BTreeSet::new())>"Clear row selection"</button><button class="button primary" disabled=move||busy.get()||selected.with(BTreeSet::is_empty) on:click=move|_|extract.run(())>"Open selected ungapped rows"</button></div>
                    <p>{move||format!("{} selected rows",selected.with(BTreeSet::len))}</p>
                    <div class="alignment-row-list">{move||alignment.get().map(|a|a.doc().rows.iter().enumerate().skip(row_page.get()*50).take(50).map(|(i,row)|{let id=row.id.clone();let checked=id.clone();let change=id.clone();let count=a.prepared.residue_count(i).unwrap_or(0);let label=format!("{} · {} residues",row.identifier,count);view!{<label class="alignment-row"><input type="checkbox" aria-label=format!("Select alignment row {}",row.identifier) prop:checked=move||selected.with(|s|s.contains(&checked)) on:change=move|e|{let chosen=event_target_checked(&e);selected.update(|s|if chosen{s.insert(change.clone());}else{s.remove(&change);});}/><button class="text-tool" on:click=move|_|focused.update(|f|f.0=i)>{label}</button></label>}}).collect_view())}</div>
                    <div class="alignment-pages"><button class="button secondary" disabled=move||row_page.get()==0 on:click=move|_|row_page.update(|p|*p=p.saturating_sub(1))>"Previous alignment rows"</button><button class="button secondary" disabled=move||{alignment.with(|a|a.as_ref().is_none_or(|a|(row_page.get()+1)*50>=a.doc().rows.len()))} on:click=move|_|row_page.update(|p|*p+=1)>"Next alignment rows"</button></div>
                </details>
                <details class="alignment-metadata"><summary>"Imported alignment annotations"</summary><p class="helper">"Symbols and free-text tags are imported evidence. They are not recalculated confidence, conservation or structure predictions."</p>{move||alignment.get().map(|a|{let doc=a.doc();let consensus=doc.metadata.clustal_consensus.as_ref().map(|s|format!("CLUSTAL consensus (first500 columns): {}",s.chars().take(500).collect::<String>()));let annotations=doc.metadata.annotations.iter().skip(annotation_page.get()*100).take(100).map(|a|format!("{} {} {}: {}{}",a.scope,a.row_identifier.as_deref().unwrap_or(""),a.tag,a.value.chars().take(500).collect::<String>(),if a.value.len()>500{"… (full value retained in saved alignment)"}else{""})).collect::<Vec<_>>();view!{<pre>{consensus}</pre>{annotations.into_iter().map(|line|view!{<p>{line}</p>}).collect_view()}<div class="alignment-pages"><button class="button secondary" disabled=move||annotation_page.get()==0 on:click=move|_|annotation_page.update(|p|*p=p.saturating_sub(1))>"Previous alignment annotations"</button><button class="button secondary" disabled=move||{alignment.with(|a|a.as_ref().is_none_or(|a|(annotation_page.get()+1)*100>=a.doc().metadata.annotations.len()))} on:click=move|_|annotation_page.update(|p|*p+=1)>"Next alignment annotations"</button></div>}})}</details>
                <div class="alignment-export"><label>"Export format"<select aria-label="Alignment export format" prop:value=move||format!("{:?}",export_format.get()) on:change=move|e|export_format.set(parse_format(&event_target_value(&e)))><option value="AlignedFasta">"Aligned FASTA"</option><option value="Clustal">"CLUSTAL"</option><option value="Stockholm">"Stockholm1.0"</option><option value="PhylipSequential">"PHYLIP sequential · strict IDs"</option><option value="PhylipInterleaved">"PHYLIP interleaved · strict IDs"</option><option value="PhylipRelaxedSequential">"PHYLIP sequential · relaxed IDs"</option><option value="PhylipRelaxedInterleaved">"PHYLIP interleaved · relaxed IDs"</option><option value="Msf">"MSF"</option><option value="Nexus">"NEXUS"</option><option value="Pir">"PIR / NBRF"</option><option value="Gde">"GDE flat"</option><option value="Selex">"SELEX"</option></select></label><button class="button secondary" disabled=move||busy.get() on:click=move|_|prepare_export.run(())>"Prepare alignment export"</button></div>
                {move||exported.get().filter(|(key,format,_)|*key==current.get()&&*format==export_format.get()).map(|(_,_,report)|view!{<div class="alignment-export-review"><h4>"Export fidelity"</h4>{report.losses.iter().map(|loss|view!{<p class="helper">{loss.clone()}</p>}).collect_view()}<button class="button primary" on:click=move|_|if let Some((key,format,report))=exported.get_untracked(){if key==current.get_untracked()&&format==export_format.get_untracked(){crate::app::download_with_completion(&format!("alignment.{}",extension(format)),"text/plain",report.text.as_bytes(),move|result|{if status.try_get_untracked().is_none(){return;}if let Err(e)=result {status.set(e);}})}}>"Download alignment export"</button></div>})}
            }})}
        </div>
    }
}

#[component]
fn AlignmentCanvas(
    alignment: AlignmentRevision,
    focused: RwSignal<(usize, usize)>,
    antibody: Memo<Option<Arc<collagene_core::alignment_antibody::Report>>>,
    candidate: RwSignal<usize>,
) -> impl IntoView {
    const LABEL: f64 = 180.;
    const CELL: f64 = 10.;
    const ROW: f64 = 26.;
    const HEADER: f64 = 30.;
    let alignment = StoredValue::new(alignment);
    let container = NodeRef::<html::Div>::new();
    let canvas = NodeRef::<html::Canvas>::new();
    let scroll = RwSignal::new((0f64, 0f64));
    let resize = StoredValue::new_local(None::<ResizeListener>);
    Effect::new(move |_| {
        if let Some(c) = container.get() {
            if resize.with_value(|r| r.is_none()) {
                if let Ok(r) = ResizeListener::new(c.unchecked_ref(), move || {
                    if let Some(s) = scroll.try_get_untracked() {
                        let _ = scroll.try_set(s);
                    }
                }) {
                    resize.set_value(Some(r));
                }
            }
        }
    });
    Effect::new(move |_| {
        let (top, left) = scroll.get();
        let focus = focused.get();
        let analysis = antibody.get();
        let candidate = candidate.get();
        let Some(c) = container.get() else {
            return;
        };
        let Some(canvas) = canvas.get() else {
            return;
        };
        let width = c.client_width().max(220) as f64;
        let height = c.client_height().max(180) as f64;
        let dpr = web_sys::window()
            .map(|w| w.device_pixel_ratio())
            .unwrap_or(1.)
            .min(2.);
        let pw = (width * dpr) as u32;
        let ph = (height * dpr) as u32;
        if canvas.width() != pw {
            canvas.set_width(pw);
        }
        if canvas.height() != ph {
            canvas.set_height(ph);
        }
        let style = canvas.unchecked_ref::<web_sys::HtmlElement>().style();
        let _ = style.set_property("width", &format!("{width}px"));
        let _ = style.set_property("height", &format!("{height}px"));
        let Some(ctx) = canvas
            .get_context("2d")
            .ok()
            .flatten()
            .and_then(|c| c.dyn_into::<Ctx>().ok())
        else {
            return;
        };
        let _ = ctx.set_transform(dpr, 0., 0., dpr, 0., 0.);
        ctx.set_fill_style_str("#fff");
        ctx.fill_rect(0., 0., width, height);
        alignment.with_value(|a| {
            let doc = a.doc();
            let left = index::logical_left(left, doc.column_count(), width);
            let start_row = (top / ROW).floor() as usize;
            let end_row = (start_row + (height / ROW).ceil() as usize + 1).min(doc.rows.len());
            let start_col = (left / CELL).floor() as usize;
            let end_col = (start_col + ((width - LABEL) / CELL).ceil().max(0.) as usize + 1)
                .min(doc.column_count());
            let mut painted = 0;
            ctx.set_font("14px monospace");
            for row in start_row..end_row {
                let y = HEADER + row as f64 * ROW - top;
                let data = &doc.rows[row];
                if row == focus.0 {
                    ctx.set_fill_style_str("#eff8f2");
                    ctx.fill_rect(LABEL, y, width - LABEL, ROW);
                }
                for col in start_col..end_col {
                    let x = LABEL + col as f64 * CELL - left;
                    if x < LABEL - 1. {
                        continue;
                    }
                    if (row, col) == focus {
                        ctx.set_fill_style_str("#badfce");
                        ctx.fill_rect(x, y, CELL, ROW);
                    }
                    if let Some(r) = analysis.as_ref().filter(|r| {
                        r.source.alignment_revision == a.revision && r.source.row_id == data.id
                    }) {
                        if let Some(residue) = r.at_column(candidate, col) {
                            if let Some(color) = residue
                                .region
                                .as_ref()
                                .and_then(crate::antibody_overlay::color)
                            {
                                ctx.set_fill_style_str(color);
                                ctx.fill_rect(x, y, CELL, ROW);
                            }
                            if residue.assignment.is_none() || residue.region.is_none() {
                                ctx.set_fill_style_str("#808080");
                                ctx.fill_rect(x, y + 22., CELL - 1., 2.);
                            }
                        }
                    }
                    let b = data.aligned.as_bytes()[col];
                    ctx.set_fill_style_str(match b {
                        b'-' => "#b9c5bf",
                        b'A' => "#337c6d",
                        b'C' => "#4476a3",
                        b'G' => "#b18b38",
                        b'T' | b'U' => "#ad6b6e",
                        _ => "#566e63",
                    });
                    let _ = ctx.fill_text(&(b as char).to_string(), x + 1., y + 18.);
                    painted += 1;
                }
                ctx.set_fill_style_str("#f4f7f4");
                ctx.fill_rect(0., y, LABEL, ROW);
                ctx.set_fill_style_str("#3e5e50");
                ctx.set_font("12px system-ui");
                let _ = ctx.fill_text(
                    &data.identifier.chars().take(24).collect::<String>(),
                    10.,
                    y + 18.,
                );
                ctx.set_font("14px monospace");
            }
            ctx.set_fill_style_str("#eaf0eb");
            ctx.fill_rect(0., 0., width, HEADER);
            ctx.set_fill_style_str("#66806e");
            ctx.set_font("11px system-ui");
            let _ = ctx.fill_text("ALIGNMENT COLUMNS", 10., 20.);
            for col in start_col..end_col {
                if col % 10 == 0 {
                    let x = LABEL + col as f64 * CELL - left;
                    if x >= LABEL {
                        let _ = ctx.fill_text(&(col + 1).to_string(), x, 20.);
                    }
                }
            }
            let _ = canvas.set_attribute("data-painted-cells", &painted.to_string());
            let _ = canvas.set_attribute("data-first-column", &start_col.to_string());
            let _ = canvas.set_attribute("data-first-row", &start_row.to_string());
        });
    });
    Effect::new(move |_| {
        let (row, col) = focused.get();
        if let Some(c) = container.get() {
            let top = c.scroll_top() as f64;
            let cols = alignment.with_value(|a| a.doc().column_count());
            let left = index::logical_left(c.scroll_left() as f64, cols, c.client_width() as f64);
            let y = row as f64 * ROW;
            let x = col as f64 * CELL;
            if y < top || y + ROW > top + c.client_height() as f64 - HEADER {
                c.set_scroll_top(y as i32);
            }
            if x < left || x + CELL > left + c.client_width() as f64 - LABEL {
                c.set_scroll_left(index::native_left(x, cols, c.client_width() as f64) as i32);
            }
            scroll.set((c.scroll_top() as f64, c.scroll_left() as f64));
        }
    });
    view! {<div class="alignment-viewport" node_ref=container on:scroll=move|_|if let Some(c)=container.get(){scroll.set((c.scroll_top()as f64,c.scroll_left()as f64));}><div class="alignment-spacer" style:height=move||alignment.with_value(|a|format!("{}px",HEADER+a.doc().rows.len()as f64*ROW)) style:width=move||alignment.with_value(|a|format!("{}px",index::spacer_width(a.doc().column_count())))><canvas class="alignment-canvas" node_ref=canvas tabindex="0" aria-label="Existing alignment. Arrow keys select a row and alignment column; Home and End select first and last columns. The selected cell is described below." on:keydown=move|e|{let(row,col)=focused.get_untracked();let(rows,cols)=alignment.with_value(|a|(a.doc().rows.len(),a.doc().column_count()));let next=match e.key().as_str(){"ArrowRight"=>Some((row,(col+1).min(cols-1))),"ArrowLeft"=>Some((row,col.saturating_sub(1))),"ArrowDown"=>Some(((row+1).min(rows-1),col)),"ArrowUp"=>Some((row.saturating_sub(1),col)),"Home"=>Some((row,0)),"End"=>Some((row,cols-1)),"PageDown"=>Some(((row+10).min(rows-1),col)),"PageUp"=>Some((row.saturating_sub(10),col)),_=>None};if let Some(next)=next{e.prevent_default();focused.set(next);}} on:click=move|e|{if let Some(c)=canvas.get(){let r=c.get_bounding_client_rect();let x=e.client_x()as f64-r.left();let y=e.client_y()as f64-r.top();if x>=LABEL&&y>=HEADER{let(top,native_left)=scroll.get_untracked();let cols=alignment.with_value(|a|a.doc().column_count());let viewport=container.get().map(|c|c.client_width()as f64).unwrap_or(220.);let left=index::logical_left(native_left,cols,viewport);let row=((y-HEADER+top)/ROW).floor()as usize;let col=((x-LABEL+left)/CELL).floor()as usize;if alignment.with_value(|a|row<a.doc().rows.len()&&col<a.doc().column_count()){focused.set((row,col));}}}}></canvas></div></div>}
}
