//! Large-reference browser adapter: bounded Blob slices and private worker sessions.
use crate::{
    app::prepare_download as download,
    genome_store::{self, DatasetRecord},
    genome_view::{GenomeWindowView, WindowDisplay},
    worker::CachedWorker,
};
use collagene_core::{
    crispr::{self, GuideOptions, GuideQuery},
    genome::{
        self, BuildInfo, EncodedBlock, ExpectedComponent, FinishedReference, GenomicInterval,
        IngestOptions, Provenance, ReferenceLink,
    },
    genome_crispr::{
        self, BlockRequest, GenomeGuideOptions, GenomeGuideReport, GenomeIntended,
        GenomeIntervalReadOptions, GenomeOffTargetHit, GenomeSearchOptions, GenomeSearchProgress,
    },
    genome_session::{self, GenomeRef, ImportStarted, RegisteredGenome, SearchPage},
    protocol::{
        request::Request,
        schema::{AnalysisJob, AnalysisResult},
    },
};
use collagene_model::Topology;
use leptos::{html, prelude::*};
use serde::Serialize;
use std::{
    cell::Cell,
    rc::Rc,
    sync::{
        Arc,
        atomic::{AtomicBool, Ordering},
    },
};
use wasm_bindgen::JsCast;
use wasm_bindgen_futures::spawn_local;
use web_sys::{Blob, File, IdbDatabase};

struct Connection {
    worker: CachedWorker,
    next: Cell<u32>,
}
impl Connection {
    fn new() -> Result<Rc<Self>, String> {
        Ok(Rc::new(Self {
            worker: CachedWorker::new()?,
            next: Cell::new(1),
        }))
    }
    async fn ask(&self, revision: &str, job: AnalysisJob) -> Result<AnalysisResult, String> {
        let id = self.next.get();
        self.next.set(id.wrapping_add(1));
        let response = self
            .worker
            .ask(Request {
                id,
                revision: revision.into(),
                job,
            })
            .await?;
        if response.id != id || response.revision != revision {
            return Err("Reference worker returned a different request/session binding.".into());
        }
        match response.result {
            AnalysisResult::Error(e) => Err(e),
            result => Ok(result),
        }
    }
    fn close(&self) {
        self.worker.close();
    }
}
#[derive(Clone, Copy)]
struct Operation {
    generation: RwSignal<u64>,
    busy: RwSignal<bool>,
    worker: StoredValue<Option<Rc<Connection>>, LocalStorage>,
    alive: StoredValue<Option<Arc<AtomicBool>>, LocalStorage>,
}
impl Operation {
    fn new() -> Self {
        Self {
            generation: RwSignal::new(0),
            busy: RwSignal::new(false),
            worker: StoredValue::new_local(None),
            alive: StoredValue::new_local(None),
        }
    }
    fn stop(self) {
        self.generation.update(|v| *v = v.wrapping_add(1));
        if let Some(alive) = self.alive.get_value() {
            alive.store(false, Ordering::Relaxed);
        }
        if let Some(w) = self.worker.get_value() {
            w.close();
        }
        self.worker.set_value(None);
        self.alive.set_value(None);
        self.busy.set(false);
    }
    fn start(self) -> Result<(Rc<Connection>, Arc<AtomicBool>, u64), String> {
        self.stop();
        let w = Connection::new()?;
        let alive = Arc::new(AtomicBool::new(true));
        self.worker.set_value(Some(w.clone()));
        self.alive.set_value(Some(alive.clone()));
        self.busy.set(true);
        Ok((w, alive, self.generation.get_untracked()))
    }
    fn complete(self) {
        if let Some(w) = self.worker.get_value() {
            w.close();
        }
        self.worker.set_value(None);
        self.alive.set_value(None);
        self.busy.set(false);
    }
    fn current(self, generation: u64) -> bool {
        self.generation.try_get_untracked() == Some(generation)
    }
}
#[derive(Clone, PartialEq, Eq)]
struct DatasetKey {
    link: ReferenceLink,
    cache_id: String,
}
#[derive(Clone)]
struct LoadedDataset {
    key: DatasetKey,
    reference: FinishedReference,
}
#[derive(Clone)]
struct BoundGuides {
    key: DatasetKey,
    report: Arc<GenomeGuideReport>,
    epoch: u64,
}
#[derive(Clone)]
struct BoundSearch {
    source: DatasetKey,
    reference: DatasetKey,
    job_id: String,
    genome: GenomeRef,
    progress: GenomeSearchProgress,
    interrupted: bool,
    epoch: u64,
}
fn whole<T: std::str::FromStr>(value: &str, label: &str) -> Result<T, String> {
    value
        .trim()
        .parse()
        .map_err(|_| format!("Enter a whole number for {label}."))
}
fn optional(value: &str) -> Option<String> {
    let value = value.trim();
    (!value.is_empty()).then(|| value.into())
}
fn coords(values: &[u64]) -> String {
    values
        .iter()
        .take(80)
        .map(|v| (v + 1).to_string())
        .collect::<Vec<_>>()
        .join(", ")
}
fn boundary(value: Option<u64>) -> String {
    value
        .map(|n| format!("{n} bases from contig origin"))
        .unwrap_or_else(|| "outside available contig".into())
}
fn role_id() -> String {
    genome_store::new_cache_id().replacen("genome-", "reference-", 1)
}
fn same_identity(link: &ReferenceLink, reference: &FinishedReference) -> Result<(), String> {
    let expected = genome::reference_link(reference, &link.role, link.cache_hint.clone())
        .map_err(|e| e.to_string())?;
    if &expected != link {
        return Err("The cache/file does not match both immutable hashes, build identity and descriptor for this reference role.".into());
    }
    Ok(())
}
fn find_record(link: &ReferenceLink, catalog: &[DatasetRecord]) -> Option<DatasetRecord> {
    let hint = link.cache_hint.as_ref()?;
    catalog
        .iter()
        .find(|r| {
            r.cache_id == *hint
                && r.manifest_hash == link.manifest_hash
                && r.content_revision == link.content_revision
                && r.build_name == link.name
                && r.assembly_accession == link.assembly_accession
                && r.version == link.version
                && r.total_bases == link.total_bases
                && r.contig_count == link.contig_count
                && r.native_length == link.file_length
        })
        .cloned()
}
async fn begin_native(
    w: &Connection,
    header: Vec<u8>,
    tail: Vec<u8>,
    length: u64,
    manifest: Vec<u8>,
) -> Result<ImportStarted, String> {
    let revision = genome_session::native_revision(&header, &tail, length, &manifest)
        .map_err(|e| e.to_string())?;
    let AnalysisResult::GenomeImportStarted(started) = w
        .ask(
            &revision,
            AnalysisJob::BeginGenomeNative {
                header,
                tail,
                file_length: length,
                manifest,
            },
        )
        .await?
    else {
        return Err("Unexpected native-reference start response.".into());
    };
    if started.binding != revision
        || started
            .reference
            .as_ref()
            .is_none_or(|r| r.manifest_hash != revision)
    {
        return Err("Native reference metadata binding mismatch.".into());
    }
    Ok(started)
}
async fn cached_metadata(
    db: &IdbDatabase,
    w: &Connection,
    record: &DatasetRecord,
    link: &ReferenceLink,
) -> Result<FinishedReference, String> {
    let bytes = genome_store::manifest_bytes(db, &record.cache_id).await?;
    let started = begin_native(
        w,
        genome::HEADER.to_vec(),
        record.tail.clone(),
        record.native_length,
        bytes,
    )
    .await?;
    let reference = started
        .reference
        .ok_or("Native metadata did not return a manifest.")?;
    same_identity(link, &reference)?;
    let _ = w
        .ask(
            &started.import_id,
            AnalysisJob::AbortGenomeImport {
                import_id: started.import_id.clone(),
            },
        )
        .await?;
    Ok(reference)
}
async fn register(w: &Connection, reference: &FinishedReference) -> Result<GenomeRef, String> {
    let AnalysisResult::GenomeRegistered(registered) = w
        .ask(
            &reference.manifest_hash,
            AnalysisJob::RegisterGenome {
                manifest: reference.manifest.clone(),
            },
        )
        .await?
    else {
        return Err("Unexpected reference registration result.".into());
    };
    if registered.binding != reference.manifest_hash
        || registered.genome.manifest_hash != reference.manifest_hash
        || registered.genome.content_revision != reference.manifest.content_revision
        || registered.manifest != reference.manifest
        || registered.semantic_verified
        || registered.receipt.is_some()
    {
        return Err("Reloaded reference registration must be exact metadata only, without restored semantic verification.".into());
    }
    Ok(registered.genome)
}
async fn stored_block(
    db: &IdbDatabase,
    cache: &str,
    reference: &FinishedReference,
    request: &BlockRequest,
) -> Result<EncodedBlock, String> {
    let contig = reference
        .manifest
        .contigs
        .iter()
        .find(|c| c.id == request.contig_id)
        .ok_or("Worker requested an unknown contig role.")?;
    let descriptor = contig
        .blocks
        .get(usize::try_from(request.index).map_err(|_| "Invalid block index.")?)
        .filter(|b| b.index == request.index)
        .ok_or("Worker requested an unknown block index.")?
        .clone();
    let encoded = genome_store::block_bytes(db, cache, &request.contig_id, request.index).await?;
    Ok(EncodedBlock {
        contig_id: request.contig_id.clone(),
        descriptor,
        encoded,
    })
}
async fn stored_blocks(
    db: &IdbDatabase,
    cache: &str,
    reference: &FinishedReference,
    requests: &[BlockRequest],
) -> Result<Vec<EncodedBlock>, String> {
    if requests.len() > 8 {
        return Err("A bounded window/search requested too many simultaneous blocks.".into());
    }
    let mut blocks = Vec::with_capacity(requests.len());
    for r in requests {
        blocks.push(stored_block(db, cache, reference, r).await?);
    }
    Ok(blocks)
}
async fn verify_cached(
    db: &IdbDatabase,
    w: &Connection,
    cache: &str,
    reference: &FinishedReference,
    alive: &Arc<AtomicBool>,
    progress: Callback<(u64, u64)>,
    retained_parts: Option<&js_sys::Array>,
) -> Result<RegisteredGenome, String> {
    let bytes = genome::encode_manifest(&reference.manifest).map_err(|e| e.to_string())?;
    let tail = genome::encode_tail(reference)
        .map_err(|e| e.to_string())?
        .to_vec();
    let mut started = begin_native(
        w,
        genome::HEADER.to_vec(),
        tail,
        reference.file_length,
        bytes,
    )
    .await?;
    let total = reference
        .manifest
        .contigs
        .iter()
        .map(|c| c.blocks.len() as u64)
        .sum();
    let mut processed = 0u64;
    while let Some(next) = started.next_block.clone() {
        genome_store::ensure_alive(alive)?;
        let stored = genome_store::block_blob(db, cache, &next.contig_id, next.index).await?;
        if stored.size() > genome::BLOCK_BASES as f64 {
            return Err("Stored export/readback block exceeds its bound.".into());
        }
        let encoded = genome_store::slice(&stored, 0, stored.size() as usize).await?;
        let descriptor = reference
            .manifest
            .contigs
            .iter()
            .find(|c| c.id == next.contig_id)
            .and_then(|c| c.blocks.get(next.index as usize))
            .ok_or("Semantic verifier requested an unknown block.")?
            .clone();
        let block = EncodedBlock {
            contig_id: next.contig_id.clone(),
            descriptor,
            encoded,
        };
        genome_store::ensure_alive(alive)?;
        let AnalysisResult::GenomeImportProgress(receipt) = w
            .ask(
                &started.import_id,
                AnalysisJob::PushGenomeNativeBlock {
                    import_id: started.import_id.clone(),
                    block,
                },
            )
            .await?
        else {
            return Err("Unexpected semantic block verification response.".into());
        };
        if receipt.import_id != started.import_id || !receipt.blocks.is_empty() {
            return Err("Semantic readback response does not bind this import session.".into());
        }
        if let Some(parts) = retained_parts {
            parts.push(stored.as_ref());
        }
        started.next_block = receipt.next_block;
        processed += 1;
        if processed > total {
            return Err("Semantic verifier requested duplicate or excess blocks.".into());
        }
        progress.run((processed, total));
    }
    genome_store::ensure_alive(alive)?;
    let AnalysisResult::GenomeRegistered(verified) = w
        .ask(
            &started.import_id,
            AnalysisJob::FinishGenomeNative {
                import_id: started.import_id.clone(),
            },
        )
        .await?
    else {
        return Err("Unexpected semantic reference finish response.".into());
    };
    if verified.binding != started.import_id
        || verified.genome.manifest_hash != reference.manifest_hash
        || verified.genome.content_revision != reference.manifest.content_revision
        || verified.manifest != reference.manifest
        || !verified.semantic_verified
        || verified.receipt.as_ref().is_none_or(|r| {
            r.manifest_hash != reference.manifest_hash
                || r.content_revision != reference.manifest.content_revision
                || r.total_bases != reference.manifest.total_bases
                || r.contigs_verified != reference.manifest.contigs.len()
        })
    {
        return Err("The stored reference did not pass full current semantic readback.".into());
    }
    Ok(verified)
}
async fn file_manifest(file: &File, w: &Connection) -> Result<ImportStarted, String> {
    let blob: &Blob = file.unchecked_ref();
    let length = file.size() as u64;
    if length < genome::HEADER.len() as u64 + genome::TAIL_BYTES as u64 {
        return Err("The .cgref file is too short.".into());
    }
    let header = genome_store::slice(blob, 0, genome::HEADER.len()).await?;
    let tail =
        genome_store::slice(blob, length - genome::TAIL_BYTES as u64, genome::TAIL_BYTES).await?;
    let info = genome::parse_tail(&header, &tail, length).map_err(|e| e.to_string())?;
    let mut manifest = Vec::with_capacity(info.manifest_length as usize);
    let mut at = 0u64;
    while at < info.manifest_length {
        let count = genome_store::READ_CHUNK.min((info.manifest_length - at) as usize);
        manifest.extend(genome_store::slice(blob, info.manifest_offset + at, count).await?);
        at += count as u64;
    }
    begin_native(w, header, tail, length, manifest).await
}
struct Json(Vec<u8>);
impl std::io::Write for Json {
    fn write(&mut self, b: &[u8]) -> std::io::Result<usize> {
        if b.len() > 32 * 1024 * 1024 - self.0.len() {
            return Err(std::io::Error::other(
                "Genome review JSON exceeds the 32 MiB browser bound",
            ));
        }
        self.0.extend_from_slice(b);
        Ok(b.len())
    }
    fn flush(&mut self) -> std::io::Result<()> {
        Ok(())
    }
}
fn json_download(name: &str, value: &impl Serialize) -> Result<crate::app::ExportArtifact, String> {
    let mut bytes = Json(vec![]);
    serde_json::to_writer_pretty(&mut bytes, value).map_err(|e| e.to_string())?;
    download(name, "application/json", &bytes.0)
}

#[component]
pub fn GenomePanel(
    links: Memo<Vec<ReferenceLink>>,
    on_link: Callback<ReferenceLink, Result<(), String>>,
    on_unlink: Callback<String, Result<(), String>>,
    epoch: RwSignal<u64>,
) -> impl IntoView {
    let db = StoredValue::new_local(None::<IdbDatabase>);
    let catalog = RwSignal::new(Vec::<DatasetRecord>::new());
    let storage_status = RwSignal::new("Opening independent large-reference storage…".to_string());
    let reload = RwSignal::new(0u64);
    let view_alive = Arc::new(AtomicBool::new(true));
    let import_op = Operation::new();
    let load_op = Operation::new();
    let window_op = Operation::new();
    let guide_op = Operation::new();
    let search_op = Operation::new();
    let export_op = Operation::new();
    let cleanup = view_alive.clone();
    on_cleanup(move || {
        cleanup.store(false, Ordering::Relaxed);
        for op in [
            import_op, load_op, window_op, guide_op, search_op, export_op,
        ] {
            op.stop();
        }
    });
    let source_role = RwSignal::new(String::new());
    let reference_role = RwSignal::new(String::new());
    let loaded = RwSignal::new(None::<Arc<LoadedDataset>>);
    let loading_status = RwSignal::new(String::new());
    let expanded = RwSignal::new(false);
    let file_input = NodeRef::<html::Input>::new();
    let reattach = RwSignal::new(None::<ReferenceLink>);
    let build_name = RwSignal::new("Imported reference collection".to_string());
    let assembly = RwSignal::new(String::new());
    let build_version = RwSignal::new(String::new());
    let provenance = RwSignal::new(String::new());
    let expected = RwSignal::new("[]".to_string());
    let total_limit = RwSignal::new("16000000000".to_string());
    let source_limit = RwSignal::new("34359738368".to_string());
    let ingest_status = RwSignal::new(String::new());
    let ingest_bytes = RwSignal::new(0u64);
    let ingest_total = RwSignal::new(0u64);
    let contig_id = RwSignal::new(String::new());
    let interval_start = RwSignal::new("1".to_string());
    let interval_end = RwSignal::new("1000".to_string());
    let contig_page = RwSignal::new(0usize);
    let window = RwSignal::new(None::<Arc<WindowDisplay>>);
    let window_status = RwSignal::new(String::new());
    let preset = RwSignal::new("spcas9".to_string());
    let guide_limit = RwSignal::new("4096".to_string());
    let guide_budget = RwSignal::new("20000000".to_string());
    let min_gc = RwSignal::new("0".to_string());
    let max_gc = RwSignal::new("100".to_string());
    let homopolymer = RwSignal::new("64".to_string());
    let guides = RwSignal::new(None::<Arc<BoundGuides>>);
    let selected = RwSignal::new(None::<String>);
    let guide_page = RwSignal::new(0usize);
    let include_rejected = RwSignal::new(false);
    let guide_status = RwSignal::new(String::new());
    let mismatch = RwSignal::new("3".to_string());
    let search_budget = RwSignal::new("1000000000000".to_string());
    let index_limit = RwSignal::new("1000000".to_string());
    let hit_limit = RwSignal::new("4096".to_string());
    let intended = RwSignal::new(false);
    let search = RwSignal::new(None::<Arc<BoundSearch>>);
    let search_status = RwSignal::new(String::new());
    let hit_page = RwSignal::new(None::<Arc<SearchPage>>);
    let fetch_busy = RwSignal::new(false);
    let source_link = Memo::new(move |_| {
        links.with(|links| links.iter().find(|l| l.role == source_role.get()).cloned())
    });
    let reference_link = Memo::new(move |_| {
        links.with(|links| {
            links
                .iter()
                .find(|l| l.role == reference_role.get())
                .cloned()
        })
    });
    let selected_guide = Memo::new(move |_| {
        guides.with(|g| {
            g.as_ref().and_then(|g| {
                selected
                    .get()
                    .and_then(|id| g.report.candidates.iter().find(|c| c.id == id).cloned())
            })
        })
    });
    let interval = Memo::new(move |_| -> Result<GenomicInterval, String> {
        let start = whole::<u64>(&interval_start.get(), "1-based global start")?
            .checked_sub(1)
            .ok_or("Start must be at least 1.")?;
        let end = whole(
            &interval_end.get(),
            "inclusive global end (unwrapped for a circular interval)",
        )?;
        Ok(GenomicInterval {
            contig_id: contig_id.get(),
            start,
            end,
        })
    });
    let guide_options = Memo::new(move |_| -> Result<GenomeGuideOptions, String> {
        let role = source_link
            .get()
            .ok_or("Select a linked genome source.")?
            .role;
        let gc = |s: String| -> Result<u16, String> {
            let n = s
                .trim()
                .parse::<u16>()
                .map_err(|_| "GC filters require whole percentages0–100.".to_string())?;
            n.checked_mul(100)
                .filter(|v| *v <= 10000)
                .ok_or("GC filter outside0–100%.".into())
        };
        Ok(GenomeGuideOptions {
            source_role: role,
            interval: interval.get()?,
            guide: GuideOptions {
                profile: if preset.get() == "cas12a" {
                    crispr::cas12a()
                } else {
                    crispr::spcas9()
                },
                min_gc_bp: gc(min_gc.get())?,
                max_gc_bp: gc(max_gc.get())?,
                max_homopolymer: whole(&homopolymer.get(), "homopolymer limit")?,
                max_candidates: whole(&guide_limit.get(), "candidate limit")?,
                max_comparisons: whole(&guide_budget.get(), "interval guide budget")?,
            },
        })
    });
    let source_record = Memo::new(move |_| {
        source_link
            .get()
            .and_then(|link| catalog.with(|records| find_record(&link, records)))
    });
    let reference_record = Memo::new(move |_| {
        reference_link
            .get()
            .and_then(|link| catalog.with(|records| find_record(&link, records)))
    });
    Effect::new({
        let alive = view_alive.clone();
        move |_| {
            let _ = reload.get();
            let alive = alive.clone();
            spawn_local(async move {
                let result = async {
                    let database = match db.get_value() {
                        Some(v) => v,
                        None => {
                            let v = genome_store::open().await?;
                            db.set_value(Some(v.clone()));
                            v
                        }
                    };
                    genome_store::list(&database).await
                }
                .await;
                if !alive.load(Ordering::Relaxed) {
                    return;
                }
                match result {
                    Ok(records) => {
                        catalog.set(records);
                        storage_status.set("Bulk references use a separate browser cache; workspace archives retain links only.".into());
                    }
                    Err(e) => storage_status.set(e),
                }
            });
        }
    });
    Effect::new(move |_| {
        let values = links.get();
        if source_role.get_untracked().is_empty() {
            if let Some(first) = values.first() {
                source_role.set(first.role.clone());
            }
        }
        if reference_role.get_untracked().is_empty() {
            if let Some(first) = values.first() {
                reference_role.set(first.role.clone());
            }
        }
    });
    Effect::new({
        let alive = view_alive.clone();
        move |_| {
            let link = source_link.get();
            let record = source_record.get();
            let expected_epoch = epoch.get();
            load_op.stop();
            loaded.set(None);
            let Some(link) = link else {
                loading_status.set("Choose or stream a genome-reference link.".into());
                return;
            };
            let Some(record) = record else {
                loading_status.set("Bulk cache missing for this linked reference. Reattach a .cgref matching both hashes; no verification state is restored by the link.".into());
                return;
            };
            let Some(database) = db.get_value() else {
                return;
            };
            let (w, op_alive, generation) = match load_op.start() {
                Ok(v) => v,
                Err(e) => {
                    loading_status.set(e);
                    return;
                }
            };
            loading_status.set("Loading bounded manifest metadata; canonical DNA is not yet verified in this session…".into());
            let alive = alive.clone();
            spawn_local(async move {
                let result = cached_metadata(&database, &w, &record, &link).await;
                if !alive.load(Ordering::Relaxed)
                    || !op_alive.load(Ordering::Relaxed)
                    || !load_op.current(generation)
                    || source_link.get_untracked() != Some(link.clone())
                    || epoch.get_untracked() != expected_epoch
                {
                    return;
                }
                load_op.complete();
                match result {
                    Ok(reference) => {
                        if let Some(c) = reference.manifest.contigs.first() {
                            contig_id.set(c.id.clone());
                            interval_start.set("1".into());
                            interval_end.set(c.length.min(1000).to_string());
                        }
                        loaded.set(Some(Arc::new(LoadedDataset {
                            key: DatasetKey {
                                link,
                                cache_id: record.cache_id,
                            },
                            reference,
                        })));
                        loading_status.set("Manifest identity loaded. Local windows verify requested blocks; whole-reference searches must verify every owner block and canonical contig hash again.".into());
                    }
                    Err(e) => loading_status.set(e),
                }
            });
        }
    });
    Effect::new(move |_| {
        let _ = (
            source_link.get(),
            reference_link.get(),
            loaded.with(|v| v.as_ref().map(|d| d.key.clone())),
            interval.get(),
            epoch.get(),
        );
        window_op.stop();
        window.set(None);
        window_status.set(String::new());
    });
    Effect::new(move |_| {
        let _ = (
            source_link.get(),
            loaded.with(|v| v.as_ref().map(|d| d.key.clone())),
            guide_options.get(),
            epoch.get(),
        );
        guide_op.stop();
        guides.set(None);
        selected.set(None);
        guide_page.set(0);
        guide_status.set(String::new());
    });
    Effect::new(move |_| {
        let _ = (
            source_link.get(),
            reference_link.get(),
            selected_guide.get(),
            guide_options.get(),
            mismatch.get(),
            search_budget.get(),
            index_limit.get(),
            hit_limit.get(),
            intended.get(),
            source_record.get(),
            reference_record.get(),
            epoch.get(),
        );
        search_op.stop();
        fetch_busy.set(false);
        search.set(None);
        hit_page.set(None);
        search_status.set(String::new());
    });
    Effect::new(move |_| {
        let _ = epoch.get();
        import_op.stop();
        reattach.set(None);
    });
    Effect::new(move |_| {
        let _ = (source_link.get(), source_record.get(), epoch.get());
        export_op.stop();
    });
    let ingest = {
        let view_alive = view_alive.clone();
        move |event: web_sys::Event| {
            let input = event_target::<web_sys::HtmlInputElement>(&event);
            let Some(file) = input.files().and_then(|files| files.get(0)) else {
                return;
            };
            input.set_value("");
            let Some(database) = db.get_value() else {
                ingest_status.set("Reference cache is not available yet.".into());
                return;
            };
            if file.name().len() > 4096
                || file.size() == 0.
                || catalog.with_untracked(Vec::len) >= 1000
            {
                ingest_status.set("Choose a nonempty file with a bounded name; at most 1,000 cached datasets are supported.".into());
                return;
            }
            let desired = reattach.get_untracked();
            let prepared = (|| -> Result<_, String> {
                let expected_components: Vec<ExpectedComponent> =
                    serde_json::from_str(&expected.get_untracked())
                        .map_err(|e| format!("Expected-component JSON: {e}"))?;
                let build = BuildInfo {
                    name: build_name.get_untracked(),
                    assembly_accession: optional(&assembly.get_untracked()),
                    version: optional(&build_version.get_untracked()),
                    provenance: optional(&provenance.get_untracked())
                        .map(|source| {
                            vec![Provenance {
                                source,
                                identifier: Some(file.name()),
                                retrieved_at: None,
                                content_hash: None,
                            }]
                        })
                        .unwrap_or_default(),
                    expected_components,
                };
                let opts = IngestOptions {
                    max_source_bytes: whole(&source_limit.get_untracked(), "source byte limit")?,
                    max_total_bases: whole(&total_limit.get_untracked(), "total base limit")?,
                    ..Default::default()
                };
                let (w, alive, generation) = import_op.start()?;
                Ok((build, opts, w, alive, generation, epoch.get_untracked()))
            })();
            let (build, opts, w, alive, generation, expected_epoch) = match prepared {
                Ok(v) => v,
                Err(e) => {
                    ingest_status.set(e);
                    return;
                }
            };
            let cache = genome_store::new_cache_id();
            let logical_role = desired
                .as_ref()
                .map(|l| l.role.clone())
                .unwrap_or_else(role_id);
            let filename = file.name();
            ingest_bytes.set(0);
            ingest_total.set(file.size() as u64);
            ingest_status.set("Staging bounded file slices; no dataset is published until stored blocks pass full semantic readback.".into());
            let view_alive = view_alive.clone();
            spawn_local(async move {
                let result=async{genome_store::stage(&database,&cache,&filename).await?;genome_store::ensure_alive(&alive)?;let source:&Blob=file.unchecked_ref();let prefix=genome_store::slice(source,0,(file.size()as u64).min(genome::HEADER.len()as u64)as usize).await?;let reference=if prefix==genome::HEADER{
            let started=file_manifest(&file,&w).await?;let reference=started.reference.ok_or("Native file has no reference manifest.")?;if let Some(link)=&desired{same_identity(link,&reference)?;}let _=w.ask(&started.import_id,AnalysisJob::AbortGenomeImport{import_id:started.import_id.clone()}).await?;for c in &reference.manifest.contigs{for d in &c.blocks{genome_store::ensure_alive(&alive)?;let bytes=genome_store::slice(source,d.offset,u64::from(d.count).div_ceil(2)as usize).await?;let block=EncodedBlock{contig_id:c.id.clone(),descriptor:d.clone(),encoded:bytes};genome_store::put_block(&database,&cache,&block).await?;if view_alive.load(Ordering::Relaxed)&&import_op.current(generation){ingest_bytes.set(d.offset+block.encoded.len()as u64);}}}reference
        }else{
            if desired.is_some(){return Err("Missing-store reattachment requires the exact hash-matching .cgref file, not a newly labeled FASTA collection.".into());}
            let revision=genome_session::fasta_revision(&build,&opts,file.size()as u64).map_err(|e|e.to_string())?;let AnalysisResult::GenomeImportStarted(started)=w.ask(&revision,AnalysisJob::BeginGenomeFasta{build,options:opts.clone(),expected_source_bytes:file.size()as u64}).await?else{return Err("Unexpected FASTA import start response.".into());};if started.binding!=revision{return Err("FASTA build/options binding mismatch.".into());}let mut offset=0u64;while offset<file.size()as u64{genome_store::ensure_alive(&alive)?;let count=(file.size()as u64-offset).min(genome_store::READ_CHUNK.min(opts.max_push_bytes)as u64)as usize;let bytes=genome_store::slice(source,offset,count).await?;let AnalysisResult::GenomeImportProgress(progress)=w.ask(&started.import_id,AnalysisJob::PushGenomeFasta{import_id:started.import_id.clone(),offset,bytes}).await?else{return Err("Unexpected FASTA block response.".into());};offset+=count as u64;if progress.import_id!=started.import_id||progress.received!=offset{return Err("FASTA source byte offset/session mismatch.".into());}for block in progress.blocks{genome_store::ensure_alive(&alive)?;genome_store::put_block(&database,&cache,&block).await?;}if view_alive.load(Ordering::Relaxed)&&import_op.current(generation){ingest_bytes.set(offset);}}
            let AnalysisResult::GenomeFastaFinished(finished)=w.ask(&started.import_id,AnalysisJob::FinishGenomeFasta{import_id:started.import_id.clone()}).await?else{return Err("Unexpected FASTA finish result.".into());};if finished.import_id!=started.import_id{return Err("FASTA finish belongs to another session.".into());}for block in finished.blocks{genome_store::ensure_alive(&alive)?;genome_store::put_block(&database,&cache,&block).await?;}finished.reference
        };genome_store::ensure_alive(&alive)?;let progress=Callback::new({let view_alive=view_alive.clone();move|(done,total):(u64,u64)|if view_alive.load(Ordering::Relaxed)&&import_op.current(generation){ingest_status.set(format!("Readback verification: {done}/{total} immutable blocks. Nothing is published yet."));}});let _verified=verify_cached(&database,&w,&cache,&reference,&alive,progress,None).await?;genome_store::ensure_alive(&alive)?;if epoch.get_untracked()!=expected_epoch{return Err("Workspace changed during import; staging was not published.".into());}let manifest=genome::encode_manifest(&reference.manifest).map_err(|e|e.to_string())?;let record=genome_store::publish(&database,&cache,&reference,&filename,&manifest).await?;if !alive.load(Ordering::Relaxed)||epoch.get_untracked()!=expected_epoch{genome_store::remove(&database,&cache).await?;return Err("Import cancelled or workspace changed before publication was accepted; staged data was removed.".into());}let link=genome::reference_link(&reference,&logical_role,Some(cache.clone())).map_err(|e|e.to_string())?;Ok::<_,String>((record,link))}.await;
                if result.is_err() {
                    let cleanup = genome_store::remove(&database, &cache).await;
                    if let Err(e) = cleanup {
                        if view_alive.load(Ordering::Relaxed) && import_op.current(generation) {
                            ingest_status.set(format!(
                                "Import failed and staging cleanup needs retry. {e}"
                            ));
                        }
                    }
                }
                if !view_alive.load(Ordering::Relaxed)
                    || !import_op.current(generation)
                    || epoch.get_untracked() != expected_epoch
                {
                    return;
                }
                import_op.complete();
                match result {
                    Ok((_record, link)) => {
                        match on_link.run(link.clone()){Ok(())=>{source_role.set(link.role.clone());reference_role.set(link.role);reattach.set(None);ingest_bytes.set(file.size()as u64);ingest_status.set("Stored blocks passed fresh canonical readback; manifest published atomically. Workspace contains a descriptor only, and no search coverage has been restored.".into());},Err(e)=>ingest_status.set(format!("Cache is verified and published, but workspace link was not retained: {e}. Link it from the cache catalog or export .cgref."))}
                        reload.update(|n| *n = n.wrapping_add(1));
                    }
                    Err(e) => ingest_status.set(format!("Reference was not published. {e}")),
                }
            });
        }
    };
    let read_window = {
        let view_alive = view_alive.clone();
        move |_| {
            let prepared = (|| -> Result<_, String> {
                let data = loaded
                    .get_untracked()
                    .ok_or("Load a linked reference first.")?;
                let opts = GenomeIntervalReadOptions {
                    source_role: data.key.link.role.clone(),
                    interval: interval.get_untracked()?,
                };
                let requests =
                    genome_crispr::read_interval_required_blocks(&data.reference.manifest, &opts)
                        .map_err(|e| e.to_string())?;
                let database = db.get_value().ok_or("Reference cache unavailable.")?;
                let (w, alive, generation) = window_op.start()?;
                Ok((
                    data,
                    opts,
                    requests,
                    database,
                    w,
                    alive,
                    generation,
                    epoch.get_untracked(),
                ))
            })();
            let (data, opts, requests, database, w, alive, generation, expected_epoch) =
                match prepared {
                    Ok(v) => v,
                    Err(e) => {
                        window_status.set(e);
                        return;
                    }
                };
            window.set(None);
            window_status
                .set("Reading only the bounded owner blocks for this global interval…".into());
            let view_alive = view_alive.clone();
            spawn_local(async move {
                let result = async {
                    let genome = register(&w, &data.reference).await?;
                    let blocks =
                        stored_blocks(&database, &data.key.cache_id, &data.reference, &requests)
                            .await?;
                    genome_store::ensure_alive(&alive)?;
                    let AnalysisResult::GenomeIntervalRead(result) = w
                        .ask(
                            &genome.manifest_hash,
                            AnalysisJob::ReadGenomeInterval {
                                genome: genome.clone(),
                                options: opts.clone(),
                                blocks,
                            },
                        )
                        .await?
                    else {
                        return Err("Unexpected global window read result.".into());
                    };
                    if result.genome != genome {
                        return Err("Window read returned a different dataset session.".into());
                    }
                    result
                        .report
                        .validate_bindings_for(&data.reference.manifest, &opts)
                        .map_err(|e| e.to_string())?;
                    let length = data
                        .reference
                        .manifest
                        .contigs
                        .iter()
                        .find(|c| c.id == opts.interval.contig_id)
                        .ok_or("Window contig missing.")?
                        .length;
                    Ok::<_, String>(Arc::new(WindowDisplay {
                        epoch: expected_epoch,
                        link: data.key.link.clone(),
                        cache_id: data.key.cache_id.clone(),
                        interval: opts.interval.clone(),
                        contig_length: length,
                        sequence: result.report.sequence,
                    }))
                }
                .await;
                if !view_alive.load(Ordering::Relaxed)
                    || !alive.load(Ordering::Relaxed)
                    || !window_op.current(generation)
                    || epoch.get_untracked() != expected_epoch
                    || loaded.with_untracked(|d| d.as_ref().map(|d| d.key.clone()))
                        != Some(data.key.clone())
                    || interval.get_untracked() != Ok(opts.interval)
                {
                    return;
                }
                window_op.complete();
                match result {
                    Ok(value) => {
                        window.set(Some(value));
                        window_status.set("Requested encoded blocks verified. This local window is not proof of full-contig canonical integrity or search coverage.".into());
                    }
                    Err(e) => window_status.set(e),
                }
            });
        }
    };
    let design_guides = {
        let view_alive = view_alive.clone();
        move |_| {
            let prepared = (|| -> Result<_, String> {
                let data = loaded
                    .get_untracked()
                    .ok_or("Load a linked reference first.")?;
                let opts = guide_options.get_untracked()?;
                let requests =
                    genome_crispr::interval_required_blocks(&data.reference.manifest, &opts)
                        .map_err(|e| e.to_string())?;
                let database = db.get_value().ok_or("Reference cache unavailable.")?;
                let (w, alive, generation) = guide_op.start()?;
                Ok((
                    data,
                    opts,
                    requests,
                    database,
                    w,
                    alive,
                    generation,
                    epoch.get_untracked(),
                ))
            })();
            let (data, opts, requests, database, w, alive, generation, expected_epoch) =
                match prepared {
                    Ok(v) => v,
                    Err(e) => {
                        guide_status.set(e);
                        return;
                    }
                };
            guides.set(None);
            selected.set(None);
            guide_status.set("Enumerating both orientations with exact global coordinates across verified owner/halo blocks…".into());
            let view_alive = view_alive.clone();
            spawn_local(async move {
                let result = async {
                    let genome = register(&w, &data.reference).await?;
                    let blocks =
                        stored_blocks(&database, &data.key.cache_id, &data.reference, &requests)
                            .await?;
                    genome_store::ensure_alive(&alive)?;
                    let AnalysisResult::GenomeIntervalGuides(result) = w
                        .ask(
                            &genome.manifest_hash,
                            AnalysisJob::GenomeIntervalGuides {
                                genome: genome.clone(),
                                options: opts.clone(),
                                blocks,
                            },
                        )
                        .await?
                    else {
                        return Err("Unexpected genome guide result.".into());
                    };
                    if result.genome != genome {
                        return Err("Guide result returned a different dataset session.".into());
                    }
                    result
                        .report
                        .validate_bindings_for(&data.reference.manifest, &opts)
                        .map_err(|e| e.to_string())?;
                    Ok::<_, String>(result.report)
                }
                .await;
                if !view_alive.load(Ordering::Relaxed)
                    || !alive.load(Ordering::Relaxed)
                    || !guide_op.current(generation)
                    || epoch.get_untracked() != expected_epoch
                    || loaded.with_untracked(|d| d.as_ref().map(|d| d.key.clone()))
                        != Some(data.key.clone())
                    || guide_options.get_untracked() != Ok(opts)
                {
                    return;
                }
                guide_op.complete();
                match result {
                    Ok(report) => {
                        selected.set(
                            report
                                .candidates
                                .iter()
                                .find(|g| g.passes_filters)
                                .map(|g| g.id.clone()),
                        );
                        guide_status.set(format!("{} global candidates · {} interval enumeration; no activity/physical-cut prediction",report.candidates.len(),if report.complete{"complete"}else{"incomplete"}));
                        guides.set(Some(Arc::new(BoundGuides {
                            key: data.key.clone(),
                            report: Arc::new(report),
                            epoch: expected_epoch,
                        })));
                    }
                    Err(e) => guide_status.set(e),
                }
            });
        }
    };
    let start_search = {
        let view_alive = view_alive.clone();
        move |_| {
            let prepared = (|| -> Result<_, String> {
                let source = guides
                    .get_untracked()
                    .ok_or("Enumerate a current source guide first.")?;
                if source.epoch != epoch.get_untracked()
                    || guide_options.get_untracked() != Ok(source.report.options.clone())
                    || loaded.with_untracked(|d| d.as_ref().map(|d| d.key.clone()))
                        != Some(source.key.clone())
                {
                    return Err(
                        "Source role, interval or guide policy changed. Enumerate again.".into(),
                    );
                }
                let guide = selected_guide
                    .get_untracked()
                    .filter(|g| g.passes_filters && g.concrete && g.pam_certain)
                    .ok_or("Choose a concrete DNA spacer with a certain PAM.")?;
                let link = reference_link
                    .get_untracked()
                    .ok_or("Select an explicit search reference role.")?;
                let record = reference_record.get_untracked().ok_or(
                    "Search reference cache is missing. Reattach the matching .cgref file.",
                )?;
                let reference = DatasetKey {
                    link: link.clone(),
                    cache_id: record.cache_id.clone(),
                };
                let locus = if intended.get_untracked() {
                    if source.key != reference {
                        return Err("Intended-site tagging requires the exact source reference role and cache binding, not merely identical sequence content.".into());
                    }
                    Some(GenomeIntended {
                        reference_role: link.role.clone(),
                        manifest_hash: link.manifest_hash.clone(),
                        content_revision: link.content_revision.clone(),
                        contig_id: guide.contig_id.clone(),
                        guide_positions: guide.guide_positions.clone(),
                        pam_positions: guide.pam_positions.clone(),
                        reverse: guide.reverse,
                    })
                } else {
                    None
                };
                let query = GuideQuery {
                    id: guide.id.clone(),
                    sequence: guide.guide_sequence,
                };
                let opts = GenomeSearchOptions {
                    reference_role: link.role.clone(),
                    profile: source.report.options.guide.profile.clone(),
                    max_mismatches: whole(&mismatch.get_untracked(), "maximum Hamming mismatches")?,
                    intended: locus,
                    max_total_comparisons: whole(
                        &search_budget.get_untracked(),
                        "total comparison limit",
                    )?,
                    max_index_entries_per_block: whole(
                        &index_limit.get_untracked(),
                        "per-block index limit",
                    )?,
                    max_hits: whole(&hit_limit.get_untracked(), "retained hit limit")?,
                };
                let database = db.get_value().ok_or("Reference database unavailable.")?;
                let (w, alive, generation) = search_op.start()?;
                Ok((
                    source,
                    reference,
                    record,
                    guide.id,
                    query,
                    opts,
                    database,
                    w,
                    alive,
                    generation,
                    epoch.get_untracked(),
                ))
            })();
            let (
                source,
                reference,
                record,
                guide_id,
                query,
                opts,
                database,
                w,
                alive,
                generation,
                expected_epoch,
            ) = match prepared {
                Ok(v) => v,
                Err(e) => {
                    search_status.set(e);
                    return;
                }
            };
            search.set(None);
            hit_page.set(None);
            search_status.set("Starting a fresh private owner-block cursor. Serialized progress cannot resume or establish coverage.".into());
            let view_alive = view_alive.clone();
            spawn_local(async move {
                let result=async{let metadata=cached_metadata(&database,&w,&record,&reference.link).await?;let genome=register(&w,&metadata).await?;let AnalysisResult::GenomeSearchState(mut state)=w.ask(&genome.manifest_hash,AnalysisJob::StartGenomeSearch{genome:genome.clone(),query:query.clone(),options:opts.clone()}).await?else{return Err("Unexpected whole-reference search start.".into());};if state.binding!=genome.manifest_hash||state.genome!=genome{return Err("Search start is not bound to this private dataset session.".into());}state.progress.validate_bindings_for(&metadata.manifest,&query,&opts).map_err(|e|e.to_string())?;
            loop{genome_store::ensure_alive(&alive)?;if epoch.get_untracked()!=expected_epoch||source_link.get_untracked()!=Some(source.key.link.clone())||reference_link.get_untracked()!=Some(reference.link.clone())||source_record.get_untracked().map(|r|r.cache_id)!=Some(source.key.cache_id.clone())||reference_record.get_untracked().map(|r|r.cache_id)!=Some(reference.cache_id.clone())||guide_options.get_untracked()!=Ok(source.report.options.clone())||selected_guide.get_untracked().map(|g|g.id)!=Some(guide_id.clone()){return Err("Reference role, source guide, options or workspace changed during the streamed search.".into());}
                if view_alive.load(Ordering::Relaxed)&&search_op.current(generation){search.set(Some(Arc::new(BoundSearch{source:source.key.clone(),reference:reference.clone(),job_id:state.job_id.clone(),genome:genome.clone(),progress:state.progress.clone(),interrupted:false,epoch:expected_epoch})));search_status.set(format!("Verified owner blocks {}/{} · {} bases · {} comparisons · {} retained matches. No completed coverage claim while work remains.",state.progress.blocks_processed,state.progress.blocks_total,state.progress.bases_verified,state.progress.comparisons,state.progress.hits_retained));}
                if state.progress.finished{break;}if state.progress.next_blocks.is_empty(){return Err("Private search has unfinished work but did not request owner/halo blocks.".into());}let blocks=stored_blocks(&database,&reference.cache_id,&metadata,&state.progress.next_blocks).await?;genome_store::ensure_alive(&alive)?;let AnalysisResult::GenomeSearchState(next)=w.ask(&state.job_id,AnalysisJob::AdvanceGenomeSearch{job_id:state.job_id.clone(),blocks}).await?else{return Err("Unexpected owner-block search response.".into());};if next.binding!=state.job_id||next.job_id!=state.job_id||next.genome!=genome{return Err("Owner-block progress does not bind the same live search.".into());}next.progress.validate_bindings_for(&metadata.manifest,&query,&opts).map_err(|e|e.to_string())?;if next.progress.blocks_processed<=state.progress.blocks_processed&&!next.progress.finished{return Err("Private search failed to advance its owner cursor.".into());}state=next;
            }
            let AnalysisResult::GenomeSearchPage(page)=w.ask(&state.job_id,AnalysisJob::PageGenomeSearch{job_id:state.job_id.clone(),offset:0,limit:25}).await?else{return Err("Unexpected final hit page.".into());};if page.job_id!=state.job_id||page.genome!=genome||page.progress!=state.progress||page.offset!=0||page.hits.len()>25{return Err("Final page does not match this completed search snapshot.".into());}Ok::<_,String>(page)}.await;
                if !view_alive.load(Ordering::Relaxed) || !search_op.current(generation) {
                    return;
                }
                search_op.busy.set(false);
                match result {
                    Ok(page) => {
                        search_status.set(if page.progress.exhaustive_for_manifest{if page.progress.complete{"Every declared expected component was present and every manifest owner was searched with canonical hashes verified. This does not authenticate provider release labels or population genomes.".into()}else{"Every imported contig owner was searched and canonical hashes verified, but expected-build conditions are absent/incomplete. No whole-build absence claim is permitted.".into()}}else{"Search ended without complete verified coverage. Review cancellation, hit/index/work limits and missing blocks.".into()});
                        hit_page.set(Some(Arc::new(page)));
                    }
                    Err(e) => {
                        search.update(|v| {
                            if let Some(old) = v.as_ref() {
                                let mut next = (**old).clone();
                                next.interrupted = true;
                                *v = Some(Arc::new(next));
                            }
                        });
                        search_status
                            .set(format!("Search cannot establish complete coverage. {e}"));
                        if let Some(w) = search_op.worker.get_value() {
                            w.close();
                        }
                        search_op.worker.set_value(None);
                    }
                }
            });
        }
    };
    let cancel_search = move |_| {
        search_op.stop();
        fetch_busy.set(false);
        search.update(|v| {
            if let Some(old) = v.as_ref() {
                let mut next = (**old).clone();
                next.interrupted = true;
                *v = Some(Arc::new(next));
            }
        });
        search_status.set("Search cancelled and private cursor disposed. Unsearched/missing data cannot support absence claims.".into());
    };
    let page_search = {
        let view_alive = view_alive.clone();
        Callback::new(move |offset: usize| {
            if search_op.busy.get_untracked() || fetch_busy.get_untracked() {
                return;
            }
            let Some(state) = search.get_untracked() else {
                return;
            };
            if state.interrupted
                || state.epoch != epoch.get_untracked()
                || source_link.get_untracked() != Some(state.source.link.clone())
                || reference_link.get_untracked() != Some(state.reference.link.clone())
            {
                search_status.set("This search scope is no longer current.".into());
                return;
            }
            let Some(w) = search_op.worker.get_value() else {
                search_status
                    .set("The private search session is unavailable. Start a fresh search.".into());
                return;
            };
            let generation = search_op.generation.get_untracked();
            fetch_busy.set(true);
            let view_alive = view_alive.clone();
            spawn_local(async move {
                let result = w
                    .ask(
                        &state.job_id,
                        AnalysisJob::PageGenomeSearch {
                            job_id: state.job_id.clone(),
                            offset,
                            limit: 25,
                        },
                    )
                    .await;
                if !view_alive.load(Ordering::Relaxed) || !search_op.current(generation) {
                    return;
                }
                fetch_busy.set(false);
                match result {
                    Ok(AnalysisResult::GenomeSearchPage(page))
                        if page.job_id == state.job_id
                            && page.genome == state.genome
                            && page.progress == state.progress
                            && page.offset == offset
                            && page.hits.len() <= 25 =>
                    {
                        hit_page.set(Some(Arc::new(page)))
                    }
                    Ok(AnalysisResult::GenomeSearchPage(_)) => {
                        search_status.set("Search page binding mismatch.".into())
                    }
                    Ok(_) => search_status.set("Unexpected search page response.".into()),
                    Err(e) => search_status.set(e),
                }
            });
        })
    };
    let export_native = {
        let view_alive = view_alive.clone();
        move |_| {
            let Some(data) = loaded.get_untracked() else {
                storage_status.set("Load a linked reference before exporting.".into());
                return;
            };
            let Some(database) = db.get_value() else {
                return;
            };
            let (w, alive, generation) = match export_op.start() {
                Ok(v) => v,
                Err(e) => {
                    storage_status.set(e);
                    return;
                }
            };
            let expected_epoch = epoch.get_untracked();
            storage_status.set(
                "Verifying stored Blob snapshots before constructing a streaming .cgref download…"
                    .into(),
            );
            let view_alive = view_alive.clone();
            spawn_local(async move {
                let result=async{let parts=js_sys::Array::new();parts.push(genome_store::blob(&genome::HEADER)?.as_ref());let progress=Callback::new({let view_alive=view_alive.clone();move|(done,total):(u64,u64)|if view_alive.load(Ordering::Relaxed)&&export_op.current(generation){storage_status.set(format!("Native export semantic readback {done}/{total} blocks; byte buffers stay bounded."));}});let _=verify_cached(&database,&w,&data.key.cache_id,&data.reference,&alive,progress,Some(&parts)).await?;genome_store::ensure_alive(&alive)?;parts.push(genome_store::blob(&genome::encode_manifest(&data.reference.manifest).map_err(|e|e.to_string())?)?.as_ref());parts.push(genome_store::blob(&genome::encode_tail(&data.reference).map_err(|e|e.to_string())?)?.as_ref());let output=Blob::new_with_blob_sequence(&parts).map_err(|_|"Cannot construct native reference Blob from immutable backing parts.")?;if output.size()as u64!=data.reference.file_length{return Err("Exported native file length differs from verified manifest framing.".into());}Ok::<_,String>(output)}.await;
                if !view_alive.load(Ordering::Relaxed)
                    || !alive.load(Ordering::Relaxed)
                    || !export_op.current(generation)
                    || epoch.get_untracked() != expected_epoch
                    || loaded.with_untracked(|d| d.as_ref().map(|d| d.key.clone()))
                        != Some(data.key.clone())
                {
                    return;
                }
                export_op.complete();
                match result{Ok(blob)=>match genome_store::download_blob("reference.cgref",&blob){Ok(())=>storage_status.set("Exported hash-checked .cgref from Blob-backed immutable parts; full DNA was never allocated as one buffer.".into()),Err(e)=>storage_status.set(e)},Err(e)=>storage_status.set(format!("Reference export was withheld. {e}"))}
            });
        }
    };
    let browse_position = {
        let view_alive = view_alive.clone();
        Callback::new(
            move |(key, contig, position, navigation_epoch): (DatasetKey, String, u64, u64)| {
                if navigation_epoch != epoch.get_untracked()
                    || !links.with_untracked(|links| links.contains(&key.link))
                {
                    window_status.set("This result reference role or workspace epoch changed. Recompute before navigation.".into());
                    return;
                }
                let Some(database) = db.get_value() else {
                    return;
                };
                let record = catalog.with_untracked(|records| find_record(&key.link, records));
                let Some(record) = record.filter(|r| r.cache_id == key.cache_id) else {
                    window_status
                        .set("Requested reference role/cache is no longer available.".into());
                    return;
                };
                let (w, alive, generation) = match window_op.start() {
                    Ok(v) => v,
                    Err(e) => {
                        window_status.set(e);
                        return;
                    }
                };
                let expected_epoch = epoch.get_untracked();
                let view_alive = view_alive.clone();
                window.set(None);
                window_status.set(
                    "Loading a bounded global window for the explicitly named result role…".into(),
                );
                spawn_local(async move {
                    let result = async {
                        let metadata = cached_metadata(&database, &w, &record, &key.link).await?;
                        let c = metadata
                            .manifest
                            .contigs
                            .iter()
                            .find(|c| c.id == contig)
                            .ok_or("Result contig is absent from its manifest.")?;
                        if position >= c.length {
                            return Err("Result position lies outside its original contig.".into());
                        }
                        let width = c.length.min(1000);
                        let start = if c.topology == Topology::Circular {
                            (position + c.length - 50.min(c.length)) % c.length
                        } else {
                            position.saturating_sub(50).min(c.length - width)
                        };
                        let opts = GenomeIntervalReadOptions {
                            source_role: key.link.role.clone(),
                            interval: GenomicInterval {
                                contig_id: contig,
                                start,
                                end: start + width,
                            },
                        };
                        let required =
                            genome_crispr::read_interval_required_blocks(&metadata.manifest, &opts)
                                .map_err(|e| e.to_string())?;
                        let blocks =
                            stored_blocks(&database, &key.cache_id, &metadata, &required).await?;
                        genome_store::ensure_alive(&alive)?;
                        let genome = register(&w, &metadata).await?;
                        let AnalysisResult::GenomeIntervalRead(result) = w
                            .ask(
                                &genome.manifest_hash,
                                AnalysisJob::ReadGenomeInterval {
                                    genome: genome.clone(),
                                    options: opts.clone(),
                                    blocks,
                                },
                            )
                            .await?
                        else {
                            return Err("Unexpected result-position window response.".into());
                        };
                        if result.genome != genome {
                            return Err("Result window dataset binding mismatch.".into());
                        }
                        result
                            .report
                            .validate_bindings_for(&metadata.manifest, &opts)
                            .map_err(|e| e.to_string())?;
                        Ok::<_, String>(Arc::new(WindowDisplay {
                            epoch: expected_epoch,
                            link: key.link.clone(),
                            cache_id: key.cache_id.clone(),
                            interval: opts.interval,
                            contig_length: c.length,
                            sequence: result.report.sequence,
                        }))
                    }
                    .await;
                    if !view_alive.load(Ordering::Relaxed)
                        || !alive.load(Ordering::Relaxed)
                        || !window_op.current(generation)
                        || epoch.get_untracked() != expected_epoch
                        || !links.with_untracked(|links| links.contains(&key.link))
                    {
                        return;
                    }
                    window_op.complete();
                    match result {
                        Ok(value) => {
                            window.set(Some(value));
                            window_status.set("Result navigation uses original global coordinates in its explicit reference role; source guide/search scope is unchanged.".into());
                        }
                        Err(e) => window_status.set(e),
                    }
                });
            },
        )
    };
    let link_cached = {
        let view_alive = view_alive.clone();
        Callback::new(move |record: DatasetRecord| {
            let Some(database) = db.get_value() else {
                return;
            };
            let (w, alive, generation) = match load_op.start() {
                Ok(v) => v,
                Err(e) => {
                    loading_status.set(e);
                    return;
                }
            };
            let view_alive = view_alive.clone();
            spawn_local(async move {
                let result = async {
                    let bytes = genome_store::manifest_bytes(&database, &record.cache_id).await?;
                    let started = begin_native(
                        &w,
                        genome::HEADER.to_vec(),
                        record.tail.clone(),
                        record.native_length,
                        bytes,
                    )
                    .await?;
                    let reference = started
                        .reference
                        .ok_or("Cached reference has no manifest.")?;
                    genome_store::ensure_alive(&alive)?;
                    let link = genome::reference_link(
                        &reference,
                        &role_id(),
                        Some(record.cache_id.clone()),
                    )
                    .map_err(|e| e.to_string())?;
                    let _ = w
                        .ask(
                            &started.import_id,
                            AnalysisJob::AbortGenomeImport {
                                import_id: started.import_id.clone(),
                            },
                        )
                        .await?;
                    Ok::<_, String>(link)
                }
                .await;
                if !view_alive.load(Ordering::Relaxed) || !load_op.current(generation) {
                    return;
                }
                load_op.complete();
                match result {
                    Ok(link) => match on_link.run(link.clone()) {
                        Ok(()) => {
                            source_role.set(link.role);
                            loading_status.set("Linked exact manifest metadata under a distinct logical role. Current canonical/search verification was not restored.".into());
                        }
                        Err(e) => loading_status.set(e),
                    },
                    Err(e) => loading_status.set(e),
                }
            });
        })
    };
    let remove_cached = {
        let view_alive = view_alive.clone();
        Callback::new(move |cache: String| {
            if import_op.busy.get_untracked()
                || search_op.busy.get_untracked()
                || export_op.busy.get_untracked()
            {
                storage_status.set(
                    "Finish or cancel active reference work before removing cached data.".into(),
                );
                return;
            }
            let Some(database) = db.get_value() else {
                return;
            };
            let view_alive = view_alive.clone();
            spawn_local(async move {
                let result = genome_store::remove(&database, &cache).await;
                if !view_alive.load(Ordering::Relaxed) {
                    return;
                }
                match result {
                    Ok(()) => {
                        reload.update(|v| *v = v.wrapping_add(1));
                        storage_status.set("Removed local bulk data. Workspace links remain explicit and require matching reattachment.".into());
                    }
                    Err(e) => storage_status.set(e),
                }
            });
        })
    };
    let export_guides = move |_| {
        if let Some(data) = guides.get_untracked() {
            if data.epoch != epoch.get_untracked()
                || source_link.get_untracked() != Some(data.key.link.clone())
                || source_record.get_untracked().map(|r| r.cache_id)
                    != Some(data.key.cache_id.clone())
                || guide_options.get_untracked() != Ok(data.report.options.clone())
            {
                guide_status
                    .set("Guide reference role, cache or options changed. Enumerate again.".into());
                return;
            }
            #[derive(Serialize)]
            struct Artifact<'a> {
                kind: &'static str,
                reference: &'a ReferenceLink,
                report: &'a GenomeGuideReport,
            }
            crate::app::complete_download(
                json_download(
                    "genome-guides.json",
                    &Artifact {
                        kind: "genome-guide-review-v1",
                        reference: &data.key.link,
                        report: &data.report,
                    },
                ),
                move |result| {
                    if guide_status.try_get_untracked().is_none() {
                        return;
                    }
                    if let Err(e) = result {
                        guide_status.set(e);
                    }
                },
            );
        }
    };
    let export_search = {
        let view_alive = view_alive.clone();
        move |_| {
            let Some(state) = search.get_untracked() else {
                return;
            };
            if state.interrupted
                || search_op.busy.get_untracked()
                || fetch_busy.get_untracked()
                || state.epoch != epoch.get_untracked()
                || source_link.get_untracked() != Some(state.source.link.clone())
                || reference_link.get_untracked() != Some(state.reference.link.clone())
            {
                search_status.set("Finish a current live search before exporting full retained evidence; cancelled state cannot resume from JSON.".into());
                return;
            }
            let Some(w) = search_op.worker.get_value() else {
                return;
            };
            let generation = search_op.generation.get_untracked();
            fetch_busy.set(true);
            let view_alive = view_alive.clone();
            spawn_local(async move {
                let result = async {
                    let mut hits = Vec::with_capacity(state.progress.hits_retained);
                    let mut offset = 0;
                    loop {
                        let AnalysisResult::GenomeSearchPage(page) = w
                            .ask(
                                &state.job_id,
                                AnalysisJob::PageGenomeSearch {
                                    job_id: state.job_id.clone(),
                                    offset,
                                    limit: 256,
                                },
                            )
                            .await?
                        else {
                            return Err("Unexpected full search evidence page.".into());
                        };
                        if page.job_id != state.job_id
                            || page.genome != state.genome
                            || page.progress != state.progress
                            || page.offset != offset
                            || page.hits.len() > 256
                        {
                            return Err("Search evidence page binding mismatch.".into());
                        }
                        hits.extend(page.hits);
                        if page.done {
                            break;
                        }
                        if page.next_offset <= offset {
                            return Err("Search evidence pagination did not advance.".into());
                        }
                        offset = page.next_offset;
                    }
                    if hits.len() != state.progress.hits_retained {
                        return Err("Retained search evidence count mismatch.".into());
                    }
                    Ok::<_, String>(hits)
                }
                .await;
                if !view_alive.load(Ordering::Relaxed) || !search_op.current(generation) {
                    return;
                }
                fetch_busy.set(false);
                match result {
                    Ok(hits) => {
                        #[derive(Serialize)]
                        struct Artifact<'a> {
                            kind: &'static str,
                            authority: &'static str,
                            source: &'a ReferenceLink,
                            reference: &'a ReferenceLink,
                            progress: &'a GenomeSearchProgress,
                            hits: &'a [GenomeOffTargetHit],
                        }
                        crate::app::complete_download(
                            json_download(
                                "genome-search.json",
                                &Artifact {
                                    kind: "genome-search-review-v1",
                                    authority: "Historical evidence only; serialized progress cannot restore a private coverage cursor",
                                    source: &state.source.link,
                                    reference: &state.reference.link,
                                    progress: &state.progress,
                                    hits: &hits,
                                },
                            ),
                            move |result| {
                                if search_status.try_get_untracked().is_none() {
                                    return;
                                }
                                if let Err(e) = result {
                                    search_status.set(e);
                                }
                            },
                        );
                    }
                    Err(e) => search_status.set(e),
                }
            });
        }
    };
    let position_notice = Callback::new(move |(role, cache, position): (String, String, u64)| {
        if window.with_untracked(|v| {
            v.as_ref().is_some_and(|v| {
                v.link.role == role && v.cache_id == cache && v.epoch == epoch.get_untracked()
            })
        }) {
            window_status.set(format!(
                "Selected original global base {} in reference role {role}.",
                position + 1
            ));
        }
    });
    view! {
        <section class="inspector-content genome-panel" class:genome-expanded=move||expanded.get()>
        <div class="genome-heading"><h2>"Genome references"</h2><button class="button secondary" on:click=move|_|expanded.update(|v|*v=!*v)>{move||if expanded.get(){"Compact panel"}else{"Expand workspace"}}</button></div>
        <p class="helper">"Large DNA references live in a separate browser cache. Windows and guide coordinates refer to original contigs; a chromosome is never opened as a small sequence document. Build labels are declarations, not verified provider releases."</p>
        <details class="genome-import" open><summary>"Stream FASTA or native .cgref"</summary>
        <label>"Build / collection name"<input aria-label="Genome build name" prop:value=move||build_name.get() on:input=move|e|build_name.set(event_target_value(&e))/></label>
        <div class="genome-fields"><label>"Assembly accession (optional)"<input aria-label="Genome assembly accession" prop:value=move||assembly.get() on:input=move|e|assembly.set(event_target_value(&e))/></label><label>"Version (optional)"<input aria-label="Genome build version" prop:value=move||build_version.get() on:input=move|e|build_version.set(event_target_value(&e))/></label></div>
        <label>"Source / provenance (optional)"<input aria-label="Genome provenance" prop:value=move||provenance.get() on:input=move|e|provenance.set(event_target_value(&e))/></label>
        <details><summary>"Declared expected components and import limits"</summary><p class="helper">"Without an expected-component list, exhaustive searches cover this imported collection only. Expected components require id, source_id, accession (or null), length, sequence_hash (or null), and topology (Linear or Circular). A native file retains its own immutable declaration."</p>
        <label>"Expected components (JSON)"<textarea aria-label="Genome expected components" prop:value=move||expected.get() on:input=move|e|expected.set(event_target_value(&e))/></label>
        <label>"Maximum source bytes"<input aria-label="Genome source byte limit" prop:value=move||source_limit.get() on:input=move|e|source_limit.set(event_target_value(&e))/></label><label>"Maximum canonical bases"<input aria-label="Genome total base limit" prop:value=move||total_limit.get() on:input=move|e|total_limit.set(event_target_value(&e))/></label><p class="helper">"One MiB source slices; 65,536-base immutable packed blocks. Up to 10,000 contigs and 250,000 blocks. Quota, malformed input or readback failure prevents publication."</p></details>
        <Show when=move||reattach.get().is_some()><p class="genome-warning">{move||reattach.get().map(|l|format!("Reattach exact .cgref for retained role {}. Both content and manifest hashes must match.",l.role))}</p><button class="button secondary" on:click=move|_|reattach.set(None)>"Cancel reattachment choice"</button></Show>
        <label class="genome-file">"Choose reference file"<input type="file" aria-label="Genome reference file" node_ref=file_input accept=".fa,.fasta,.fna,.ffn,.cgref,text/plain" disabled=move||import_op.busy.get() on:change=ingest/></label>
        <progress max=move||ingest_total.get().max(1).to_string() value=move||ingest_bytes.get().to_string()/><p class="helper genome-ingest-status" role="status">{move||ingest_status.get()}</p>
        <button class="button secondary" disabled=move||!import_op.busy.get() on:click=move|_|{import_op.stop();ingest_status.set("Import cancelled. Its unpublished staging namespace is being removed; no search or semantic authority is retained.".into());}>"Cancel reference import"</button>
        </details>
        <h3>"Workspace reference roles"</h3><p class="helper genome-storage-status" role="status">{move||storage_status.get()}</p>
        {move||links.get().into_iter().map(|link|{let unlink=link.role.clone();let missing=find_record(&link,&catalog.get()).is_none();let attach=StoredValue::new(link.clone());view!{<details class="genome-link"><summary>{format!("{} · {} · {} bases{}",link.name,link.role,link.total_bases,if missing{" · cache missing"}else{""})}</summary><p>{format!("{} contigs · accession {} · version {}",link.contig_count,link.assembly_accession.as_deref().unwrap_or("undeclared"),link.version.as_deref().unwrap_or("undeclared"))}</p><p class="genome-hash">{format!("Manifest: {}\nCanonical content: {}",link.manifest_hash,link.content_revision)}</p><Show when=move||missing><button class="button secondary" on:click=move|_|{reattach.set(Some(attach.get_value()));if let Some(input)=file_input.get(){input.click();}}>"Reattach matching .cgref"</button></Show><button class="button secondary" on:click=move|_|if let Err(e)=on_unlink.run(unlink.clone()){storage_status.set(e);}>"Remove workspace link"</button></details>}}).collect_view()}
        <details><summary>"Local bulk cache catalog"</summary><p class="helper">"Linking cached data creates a distinct reference role. Removing cached data leaves saved workspace links explicitly missing."</p>{move||catalog.get().into_iter().take(1000).map(|r|{let link=r.clone();let cache=r.cache_id.clone();view!{<div class="genome-cache-entry"><p>{format!("{} · {} bases · {} contigs",r.build_name,r.total_bases,r.contig_count)}</p><code>{r.cache_id}</code><button class="button secondary" on:click=move|_|link_cached.run(link.clone())>"Link under new role"</button><button class="button secondary" on:click=move|_|remove_cached.run(cache.clone())>"Remove cached data"</button></div>}}).collect_view()}<p class="helper">"Failed imports remove their own staging. Interrupted browser sessions may leave unpublished storage; clearing site data removes all local caches."</p></details>
        <label>"Source reference role"<select aria-label="Genome source reference" prop:value=move||source_role.get() on:change=move|e|source_role.set(event_target_value(&e))><option value="">"Choose source reference"</option>{move||links.get().into_iter().map(|l|view!{<option value=l.role.clone()>{format!("{} · {}",l.name,l.role)}</option>}).collect_view()}</select></label>
        <p class="helper genome-loading-status" role="status">{move||loading_status.get()}</p>
        {move||loaded.get().map(|d|{let n=d.reference.manifest.contigs.len();view!{<div class="genome-manifest-summary"><p>{format!("{} canonical bases · {} contigs · {} immutable blocks",d.reference.manifest.total_bases,n,d.reference.manifest.contigs.iter().map(|c|c.blocks.len()).sum::<usize>())}</p><p class="genome-warning">{if d.reference.manifest.expected_complete{"Declared expected components are present. Search still requires fresh canonical owner coverage.".to_string()}else{format!("Expected build completeness is not established ({} missing declared components). Exhaustive results cover imported contigs only.",d.reference.manifest.missing_expected.len())}}</p><div class="genome-contigs">{d.reference.manifest.contigs.iter().skip(contig_page.get()*25).take(25).map(|c|{let id=c.id.clone();let length=c.length;view!{<button class="result-entry" on:click=move|_|{contig_id.set(id.clone());interval_start.set("1".into());interval_end.set(length.min(1000).to_string());}>{format!("{} · {} bases · {}",c.id,c.length,if c.topology==Topology::Circular{"circular"}else{"linear"})}</button>}}).collect_view()}</div><div class="genome-actions"><button class="button secondary" disabled=move||contig_page.get()==0 on:click=move|_|contig_page.update(|p|*p=p.saturating_sub(1))>"Previous contigs"</button><button class="button secondary" disabled={move||(contig_page.get()+1)*25>=n} on:click=move|_|contig_page.update(|p|*p+=1)>"Next contigs"</button></div></div>}})}
        <h3>"Global interval"</h3><label>"Contig ID"<input aria-label="Genome contig ID" prop:value=move||contig_id.get() on:input=move|e|contig_id.set(event_target_value(&e))/></label><div class="genome-fields"><label>"First global base (1-based)"<input aria-label="Genome interval start" prop:value=move||interval_start.get() on:input=move|e|interval_start.set(event_target_value(&e))/></label><label>"Last global base (inclusive)"<input aria-label="Genome interval end" prop:value=move||interval_end.get() on:input=move|e|interval_end.set(event_target_value(&e))/></label></div>
        <p class="helper">"At most 200,000 owner starts per interval. Circular end coordinates may continue past the origin for one traversal. Guide ownership follows the forward full-footprint start, including reverse-strand guides."</p>
        <button class="button secondary" disabled=move||window_op.busy.get()||loaded.get().is_none() on:click=read_window>"Read global window"</button><p class="helper genome-window-status" role="status">{move||window_status.get()}</p>{move||window.get().map(|data|view!{<GenomeWindowView data=data on_position=position_notice/>})}
        <h3>"Interval guide design"</h3><label>"Nominal nuclease profile"<select aria-label="Genome nuclease profile" prop:value=move||preset.get() on:change=move|e|preset.set(event_target_value(&e))><option value="spcas9">"SpCas9 · NGG"</option><option value="cas12a">"Cas12a · TTTV"</option></select></label>
        <div class="genome-fields"><label>"Minimum GC %"<input aria-label="Genome minimum GC" prop:value=move||min_gc.get() on:input=move|e|min_gc.set(event_target_value(&e))/></label><label>"Maximum GC %"<input aria-label="Genome maximum GC" prop:value=move||max_gc.get() on:input=move|e|max_gc.set(event_target_value(&e))/></label></div><label>"Longest homopolymer"<input aria-label="Genome homopolymer limit" prop:value=move||homopolymer.get() on:input=move|e|homopolymer.set(event_target_value(&e))/></label>
        <details><summary>"Guide resource limits"</summary><label>"Maximum retained candidates"<input aria-label="Genome guide candidate limit" prop:value=move||guide_limit.get() on:input=move|e|guide_limit.set(event_target_value(&e))/></label><label>"Maximum guide comparisons"<input aria-label="Genome guide comparison limit" prop:value=move||guide_budget.get() on:input=move|e|guide_budget.set(event_target_value(&e))/></label></details>
        <button class="button primary" disabled=move||guide_op.busy.get()||loaded.get().is_none() on:click=design_guides>"Enumerate interval guides"</button><button class="button secondary" disabled=move||!guide_op.busy.get() on:click=move|_|{guide_op.stop();guides.set(None);guide_status.set("Guide enumeration cancelled.".into());}>"Cancel guide enumeration"</button><p class="helper genome-guide-status" role="status">{move||guide_status.get()}</p>
        <label class="checkbox-label"><input type="checkbox" aria-label="Genome show rejected guides" prop:checked=move||include_rejected.get() on:change=move|e|{include_rejected.set(event_target_checked(&e));guide_page.set(0);}/>"Show rejected or ambiguous candidates"</label>
        {move||guides.get().map(|data|{let visible=data.report.candidates.iter().filter(|g|include_rejected.get()||g.passes_filters).cloned().collect::<Vec<_>>();let count=visible.len();view!{<div class="genome-guide-results"><p>{format!("{} visible candidates · enumeration {}",count,if data.report.complete{"finished under declared scope"}else{"incomplete"})}</p>{data.report.warnings.iter().map(|w|view!{<p class="genome-warning">{w.clone()}</p>}).collect_view()}{data.report.limits.iter().map(|w|view!{<p class="genome-warning">{w.clone()}</p>}).collect_view()}{visible.into_iter().skip(guide_page.get()*25).take(25).map(|g|{let choose=g.id.clone();let key=data.key.clone();let navigation_epoch=data.epoch;let navcontig=g.contig_id.clone();let position=g.guide_positions.first().copied().unwrap_or(0);view!{<div class="genome-guide-entry" class:selected=move||selected.get().as_deref()==Some(choose.as_str())><button class="result-entry" disabled=!g.passes_filters on:click={let id=g.id.clone();move|_|selected.set(Some(id.clone()))}><strong>{g.guide_sequence.clone()}</strong>{format!(" · PAM {} · {} · {}",g.pam_sequence,if g.reverse{"reverse"}else{"forward"},if g.passes_filters{"eligible"}else{"rejected"})}</button><p class="helper">"DNA spacer (5′→3′), not a complete guide RNA/scaffold"</p><p>{format!("{} · original global guide bases {} · PAM bases {}",g.contig_id,coords(&g.guide_positions),coords(&g.pam_positions))}</p><p>{format!("Nominal top / bottom boundaries: {} / {}",boundary(g.nominal_top_cut),boundary(g.nominal_bottom_cut))}</p><p class="genome-warning">{g.rejections.join("; ")}</p><button class="button secondary" on:click=move|_|browse_position.run((key.clone(),navcontig.clone(),position,navigation_epoch))>"Browse original guide locus"</button></div>}}).collect_view()}<div class="genome-actions"><button class="button secondary" disabled=move||guide_page.get()==0 on:click=move|_|guide_page.update(|p|*p=p.saturating_sub(1))>"Previous guides"</button><button class="button secondary" disabled={move||(guide_page.get()+1)*25>=count} on:click=move|_|guide_page.update(|p|*p+=1)>"Next guides"</button><button class="button secondary" on:click=export_guides>"Export genome guides JSON"</button></div></div>}})}
        <h3>"Whole imported-reference search"</h3><p class="helper">"Hamming mismatches only; no bulges or PAMs outside the selected profile. Ambiguous reference bases retain a mismatch range and uncertainty. No efficacy or cellular editing predictions. Every owner block and canonical contig hash must be checked before exhaustive imported-collection coverage can be reported."</p>
        <label>"Explicit search reference role"<select aria-label="Genome search reference" prop:value=move||reference_role.get() on:change=move|e|reference_role.set(event_target_value(&e))><option value="">"Choose search reference"</option>{move||links.get().into_iter().map(|l|view!{<option value=l.role.clone()>{format!("{} · {}",l.name,l.role)}</option>}).collect_view()}</select></label><label>"Maximum Hamming mismatches"<input aria-label="Genome maximum mismatches" prop:value=move||mismatch.get() on:input=move|e|mismatch.set(event_target_value(&e))/></label><label class="checkbox-label"><input type="checkbox" aria-label="Genome tag exact intended locus" prop:checked=move||intended.get() on:change=move|e|intended.set(event_target_checked(&e))/>"Tag the exact source-role locus as intended (retain the hit)"</label>
        <details><summary>"Search resource limits"</summary><label>"Total comparisons"<input aria-label="Genome search comparison limit" prop:value=move||search_budget.get() on:input=move|e|search_budget.set(event_target_value(&e))/></label><label>"Index entries per owner block"<input aria-label="Genome search index limit" prop:value=move||index_limit.get() on:input=move|e|index_limit.set(event_target_value(&e))/></label><label>"Retained hit limit"<input aria-label="Genome search hit limit" prop:value=move||hit_limit.get() on:input=move|e|hit_limit.set(event_target_value(&e))/></label></details>
        <div class="genome-actions"><button class="button primary" disabled=move||search_op.busy.get()||selected_guide.get().is_none() on:click=start_search>"Search all imported contigs"</button><button class="button secondary" disabled=move||!search_op.busy.get() on:click=cancel_search>"Cancel genome search"</button></div><p class="helper genome-search-status" role="status">{move||search_status.get()}</p>
        {move||search.get().map(|s|view!{<div class="genome-search-progress" data-interrupted=s.interrupted.to_string()><progress max=s.progress.blocks_total.max(1).to_string() value=s.progress.blocks_processed.to_string()/><p>{format!("Owner blocks {}/{} · canonical bases verified {} · comparisons {} · retained hits {}",s.progress.blocks_processed,s.progress.blocks_total,s.progress.bases_verified,s.progress.comparisons,s.progress.hits_retained)}</p><p>{format!("Exhaustive for imported manifest: {} · canonical source integrity verified: {} · declared expected components complete: {} · complete under all declared conditions: {}",s.progress.exhaustive_for_manifest&&!s.interrupted,s.progress.source_integrity_verified&&!s.interrupted,s.progress.expected_components_complete,s.progress.complete&&!s.interrupted)}</p>{s.progress.warnings.iter().chain(s.progress.limits.iter()).map(|w|view!{<p class="genome-warning">{w.clone()}</p>}).collect_view()}<Show when={let interrupted=s.interrupted;move||interrupted}><p class="genome-warning">"Interrupted evidence cannot establish absence or resume a private search cursor."</p></Show></div>})}
        {move||hit_page.get().map(|page|{let offset=page.offset;let next=page.next_offset;let done=page.done;view!{<div class="genome-hit-results">{page.hits.iter().cloned().map(|hit|{let key=search.get_untracked().map(|s|(s.reference.clone(),s.epoch));let contig=hit.contig_id.clone();let position=hit.guide_positions.first().copied().unwrap_or(0);view!{<div class="genome-hit-entry"><strong>{format!("{} · {} · mismatches {}–{}{}",hit.contig_id,if hit.reverse{"reverse"}else{"forward"},hit.min_mismatches,hit.max_mismatches,if hit.is_intended{" · exact intended locus"}else{""})}</strong><p>{format!("DNA spacer {} · PAM {} · {}",hit.reference_guide,hit.pam_sequence,if hit.definite{"certain match"}else{"possible match / ambiguous evidence"})}</p><p>{format!("Original global guide bases {} · PAM bases {}",coords(&hit.guide_positions),coords(&hit.pam_positions))}</p><button class="button secondary" on:click=move|_|if let Some((key,navigation_epoch))=&key{browse_position.run((key.clone(),contig.clone(),position,*navigation_epoch));}>"Browse original hit locus"</button></div>}}).collect_view()}<div class="genome-actions"><button class="button secondary" disabled=move||offset==0||fetch_busy.get() on:click=move|_|page_search.run(offset.saturating_sub(25))>"Previous hits"</button><button class="button secondary" disabled=move||done||fetch_busy.get() on:click=move|_|page_search.run(next)>"Next hits"</button></div></div>}})}
        <button class="button secondary" disabled=move||search_op.busy.get()||fetch_busy.get()||search.with(|s|s.as_ref().is_none_or(|s|s.interrupted)) on:click=export_search>"Export genome search JSON"</button>
        <h3>"Portable reference cache"</h3><p class="helper">"Workspace/native archives contain only small reference links. Export .cgref separately to retain bulk blocks. Reopening a link or a report restores no current verified search authority."</p><button class="button secondary" disabled=move||export_op.busy.get()||loaded.get().is_none() on:click=export_native>"Export reference .cgref"</button><button class="button secondary" disabled=move||!export_op.busy.get() on:click=move|_|{export_op.stop();storage_status.set("Reference export cancelled; no verified download was committed.".into());}>"Cancel reference export"</button>
        </section>
    }
}
