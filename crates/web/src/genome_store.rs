//! Blob-backed large-reference cache, deliberately separate from workspace storage.
use collagene_core::genome::{EncodedBlock, FinishedReference};
use serde::{Deserialize, Serialize};
use std::{
    cell::RefCell,
    rc::Rc,
    sync::{
        Arc,
        atomic::{AtomicBool, Ordering},
    },
};
use wasm_bindgen::{JsCast, JsValue, closure::Closure};
use wasm_bindgen_futures::JsFuture;
use web_sys::{
    Blob, Event, IdbDatabase, IdbKeyRange, IdbRequest, IdbTransaction, IdbTransactionMode,
};
const DATABASE: &str = "collagene-genomes";
const CATALOG: &str = "catalog";
const MANIFESTS: &str = "manifests";
const BLOCKS: &str = "blocks";
const STAGING: &str = "staging";
pub const READ_CHUNK: usize = 1024 * 1024;
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct DatasetRecord {
    pub cache_id: String,
    pub manifest_hash: String,
    pub content_revision: String,
    pub build_name: String,
    pub assembly_accession: Option<String>,
    pub version: Option<String>,
    pub contig_count: usize,
    pub total_bases: u64,
    pub native_length: u64,
    pub manifest_offset: u64,
    pub tail: Vec<u8>,
    pub expected_complete: bool,
    pub missing_expected_count: usize,
    pub source_filename: String,
}
#[derive(Debug, Clone, Serialize, Deserialize)]
struct StagingRecord {
    cache_id: String,
    source_filename: String,
}
async fn request(request: &IdbRequest) -> Result<JsValue, String> {
    let (tx, rx) = futures_channel::oneshot::channel();
    let tx = Rc::new(RefCell::new(Some(tx)));
    let success = {
        let tx = tx.clone();
        let request = request.clone();
        Closure::<dyn FnMut(Event)>::new(move |_| {
            if let Some(tx) = tx.borrow_mut().take() {
                let _ = tx.send(
                    request
                        .result()
                        .map_err(|_| "Genome cache read failed.".into()),
                );
            }
        })
    };
    let error = {
        let tx = tx.clone();
        Closure::<dyn FnMut(Event)>::new(move |_| {
            if let Some(tx) = tx.borrow_mut().take() {
                let _ = tx.send(Err(
                    "Genome cache request failed (possibly quota or storage access).".into(),
                ));
            }
        })
    };
    request.set_onsuccess(Some(success.as_ref().unchecked_ref()));
    request.set_onerror(Some(error.as_ref().unchecked_ref()));
    let result = rx
        .await
        .map_err(|_| "Genome cache request was interrupted.".to_string());
    request.set_onsuccess(None);
    request.set_onerror(None);
    result?
}
async fn committed(tx: &IdbTransaction) -> Result<(), String> {
    let (sender, receiver) = futures_channel::oneshot::channel();
    let sender = Rc::new(RefCell::new(Some(sender)));
    let complete = {
        let sender = sender.clone();
        Closure::<dyn FnMut(Event)>::new(move |_| {
            if let Some(s) = sender.borrow_mut().take() {
                let _ = s.send(Ok(()));
            }
        })
    };
    let aborted = {
        let sender = sender.clone();
        Closure::<dyn FnMut(Event)>::new(move |_| {
            if let Some(s) = sender.borrow_mut().take() {
                let _ = s.send(Err(
                    "Genome cache transaction failed; no completed dataset was published.".into(),
                ));
            }
        })
    };
    tx.set_oncomplete(Some(complete.as_ref().unchecked_ref()));
    tx.set_onabort(Some(aborted.as_ref().unchecked_ref()));
    tx.set_onerror(Some(aborted.as_ref().unchecked_ref()));
    let result = receiver
        .await
        .map_err(|_| "Genome cache transaction was interrupted.".to_string());
    tx.set_oncomplete(None);
    tx.set_onabort(None);
    tx.set_onerror(None);
    result?
}
fn stores(db: &IdbDatabase, names: &[&str]) -> Result<IdbTransaction, String> {
    let list = js_sys::Array::new();
    for name in names {
        list.push(&JsValue::from_str(name));
    }
    db.transaction_with_str_sequence_and_mode(&list, IdbTransactionMode::Readwrite)
        .map_err(|_| "Cannot start genome cache transaction.".into())
}
fn queue_error(tx: &IdbTransaction, error: &str) -> String {
    let _ = tx.abort();
    error.into()
}
pub async fn open() -> Result<IdbDatabase, String> {
    let factory = web_sys::window()
        .ok_or("No browser window.")?
        .indexed_db()
        .map_err(|_| "IndexedDB is unavailable.")?
        .ok_or("IndexedDB is unavailable.")?;
    let pending = factory
        .open_with_u32(DATABASE, 1)
        .map_err(|_| "Cannot open the large-reference cache.")?;
    let upgrade = {
        let pending = pending.clone();
        Closure::<dyn FnMut(Event)>::new(move |_| {
            if let Ok(db) = pending.result().and_then(|v| v.dyn_into::<IdbDatabase>()) {
                for store in [CATALOG, MANIFESTS, BLOCKS, STAGING] {
                    let _ = db.create_object_store(store);
                }
            }
        })
    };
    pending.set_onupgradeneeded(Some(upgrade.as_ref().unchecked_ref()));
    let result = request(pending.unchecked_ref()).await;
    pending.set_onupgradeneeded(None);
    result?
        .dyn_into::<IdbDatabase>()
        .map_err(|_| "Invalid genome cache connection.".into())
}
pub fn new_cache_id() -> String {
    format!(
        "genome-{:x}-{:08x}-{:08x}",
        js_sys::Date::now() as u64,
        (js_sys::Math::random() * u32::MAX as f64) as u32,
        (js_sys::Math::random() * u32::MAX as f64) as u32
    )
}
fn block_key(cache: &str, contig: &str, index: u64) -> String {
    format!("{cache}/{contig}/{index:016x}")
}
fn valid_cache(cache: &str) -> bool {
    cache.starts_with("genome-")
        && cache.len() <= 100
        && cache
            .bytes()
            .all(|b| b.is_ascii_alphanumeric() || b == b'-')
}
pub fn blob(bytes: &[u8]) -> Result<Blob, String> {
    let parts = js_sys::Array::new();
    parts.push(&js_sys::Uint8Array::from(bytes));
    Blob::new_with_u8_array_sequence(&parts)
        .map_err(|_| "Cannot create a bounded reference block Blob.".into())
}
pub async fn slice(file: &Blob, start: u64, length: usize) -> Result<Vec<u8>, String> {
    if length > READ_CHUNK {
        return Err("Reference read exceeds the 1 MiB chunk bound.".into());
    }
    let end = start
        .checked_add(length as u64)
        .ok_or("Invalid reference byte range.")?;
    if end > file.size() as u64 {
        return Err("Reference byte range is outside the file.".into());
    }
    let part = file
        .slice_with_f64_and_f64(start as f64, end as f64)
        .map_err(|_| "Cannot slice the reference file.".to_string())?;
    let buffer = JsFuture::from(part.array_buffer())
        .await
        .map_err(|_| "Cannot read the bounded reference file slice.".to_string())?;
    Ok(js_sys::Uint8Array::new(&buffer).to_vec())
}
pub async fn stage(db: &IdbDatabase, cache: &str, source_filename: &str) -> Result<(), String> {
    if !valid_cache(cache) {
        return Err("Invalid reference cache identity.".into());
    }
    let value = serde_json::to_string(&StagingRecord {
        cache_id: cache.into(),
        source_filename: source_filename.into(),
    })
    .map_err(|e| e.to_string())?;
    let tx = stores(db, &[STAGING])?;
    if tx
        .object_store(STAGING)
        .map_err(|_| queue_error(&tx, "Cannot access reference staging store."))?
        .put_with_key(&JsValue::from_str(&value), &JsValue::from_str(cache))
        .is_err()
    {
        return Err(queue_error(&tx, "Cannot stage the reference import."));
    }
    committed(&tx).await
}
pub async fn put_block(db: &IdbDatabase, cache: &str, block: &EncodedBlock) -> Result<(), String> {
    if !valid_cache(cache) || block.encoded.len() > collagene_core::genome::BLOCK_BASES {
        return Err("Invalid bounded encoded genome block.".into());
    }
    let value = blob(&block.encoded)?;
    let tx = stores(db, &[BLOCKS])?;
    if tx
        .object_store(BLOCKS)
        .map_err(|_| queue_error(&tx, "Cannot access reference blocks."))?
        .put_with_key(
            value.as_ref(),
            &JsValue::from_str(&block_key(cache, &block.contig_id, block.descriptor.index)),
        )
        .is_err()
    {
        return Err(queue_error(
            &tx,
            "Cannot write reference block; quota or storage access may be exhausted.",
        ));
    }
    committed(&tx).await
}
pub async fn block_blob(
    db: &IdbDatabase,
    cache: &str,
    contig: &str,
    index: u64,
) -> Result<Blob, String> {
    let tx = db
        .transaction_with_str(BLOCKS)
        .map_err(|_| "Cannot read reference block store.")?;
    let pending = tx
        .object_store(BLOCKS)
        .map_err(|_| "Cannot access reference block store.")?
        .get(&JsValue::from_str(&block_key(cache, contig, index)))
        .map_err(|_| "Cannot request reference block.")?;
    let value = request(&pending).await?;
    if value.is_undefined() {
        return Err(format!(
            "Cached reference block is missing ({contig}, block {index}). Reattach the exact .cgref file; no complete search can be claimed."
        ));
    }
    value
        .dyn_into::<Blob>()
        .map_err(|_| "Reference block is not a valid stored Blob.".into())
}
pub async fn block_bytes(
    db: &IdbDatabase,
    cache: &str,
    contig: &str,
    index: u64,
) -> Result<Vec<u8>, String> {
    let value = block_blob(db, cache, contig, index).await?;
    if value.size() > collagene_core::genome::BLOCK_BASES as f64 {
        return Err("Stored reference block exceeds its bound.".into());
    }
    slice(&value, 0, value.size() as usize).await
}
pub async fn publish(
    db: &IdbDatabase,
    cache: &str,
    finished: &FinishedReference,
    source_filename: &str,
    manifest_bytes: &[u8],
) -> Result<DatasetRecord, String> {
    if !valid_cache(cache) || manifest_bytes.len() > collagene_core::genome::MAX_MANIFEST_BYTES {
        return Err("Invalid bounded reference manifest.".into());
    }
    let record = DatasetRecord {
        cache_id: cache.into(),
        manifest_hash: finished.manifest_hash.clone(),
        content_revision: finished.manifest.content_revision.clone(),
        build_name: finished.manifest.build.name.clone(),
        assembly_accession: finished.manifest.build.assembly_accession.clone(),
        version: finished.manifest.build.version.clone(),
        contig_count: finished.manifest.contigs.len(),
        total_bases: finished.manifest.total_bases,
        native_length: finished.file_length,
        manifest_offset: finished.manifest_offset,
        tail: collagene_core::genome::encode_tail(finished)
            .map_err(|e| e.to_string())?
            .to_vec(),
        expected_complete: finished.manifest.expected_complete,
        missing_expected_count: finished.manifest.missing_expected.len(),
        source_filename: source_filename.into(),
    };
    let descriptor = serde_json::to_string(&record).map_err(|e| e.to_string())?;
    let manifest = blob(manifest_bytes)?;
    let tx = stores(db, &[CATALOG, MANIFESTS, STAGING])?;
    let capacity = tx
        .object_store(CATALOG)
        .map_err(|_| queue_error(&tx, "Cannot check catalog capacity."))?
        .count()
        .map_err(|_| queue_error(&tx, "Cannot request catalog capacity."))?;
    let pending = tx
        .object_store(STAGING)
        .map_err(|_| queue_error(&tx, "Cannot check staging owner."))?
        .get(&JsValue::from_str(cache))
        .map_err(|_| queue_error(&tx, "Cannot request staging owner."))?;
    let publication_error = Rc::new(RefCell::new(None::<String>));
    let queue = {
        let publication_error = publication_error.clone();
        let tx = tx.clone();
        let pending = pending.clone();
        let cache = cache.to_string();
        Closure::<dyn FnMut(Event)>::new(move |_| {
            let result = (|| -> Result<(), String> {
                if capacity
                    .result()
                    .ok()
                    .and_then(|v| v.as_f64())
                    .is_none_or(|n| n >= 1000.)
                {
                    return Err("Reference cache capacity reached; remove a cached dataset before publishing another.".into());
                }
                let marker = pending.result().map_err(|_| "Cannot read staging owner.")?;
                if marker
                    .as_string()
                    .and_then(|s| serde_json::from_str::<StagingRecord>(&s).ok())
                    .is_none_or(|owner| owner.cache_id != cache)
                {
                    return Err("Reference staging was removed before publication.".into());
                }
                tx.object_store(MANIFESTS)
                    .map_err(|_| "Cannot access manifest store.")?
                    .put_with_key(manifest.as_ref(), &JsValue::from_str(&cache))
                    .map_err(|_| "Cannot write completed manifest.")?;
                tx.object_store(CATALOG)
                    .map_err(|_| "Cannot access catalog store.")?
                    .put_with_key(&JsValue::from_str(&descriptor), &JsValue::from_str(&cache))
                    .map_err(|_| "Cannot publish reference descriptor.")?;
                tx.object_store(STAGING)
                    .map_err(|_| "Cannot access staging store.")?
                    .delete(&JsValue::from_str(&cache))
                    .map_err(|_| "Cannot finish staging.")?;
                Ok(())
            })();
            if let Err(error) = result {
                *publication_error.borrow_mut() = Some(error);
                let _ = tx.abort();
            }
        })
    };
    pending.set_onsuccess(Some(queue.as_ref().unchecked_ref()));
    let outcome = committed(&tx).await;
    pending.set_onsuccess(None);
    outcome.map_err(|fallback| publication_error.borrow_mut().take().unwrap_or(fallback))?;
    Ok(record)
}
async fn all(db: &IdbDatabase, store: &str) -> Result<js_sys::Array, String> {
    let tx = db
        .transaction_with_str(store)
        .map_err(|_| "Cannot read reference cache catalog.")?;
    let pending = tx
        .object_store(store)
        .map_err(|_| "Cannot access reference cache catalog.")?
        .get_all_with_key_and_limit(&JsValue::UNDEFINED, 1001)
        .map_err(|_| "Cannot request reference cache catalog.")?;
    Ok(js_sys::Array::from(&request(&pending).await?))
}
pub async fn list(db: &IdbDatabase) -> Result<Vec<DatasetRecord>, String> {
    let values = all(db, CATALOG).await?;
    if values.length() > 1000 {
        return Err("Large-reference catalog exceeds the 1,000-dataset browser bound.".into());
    }
    values
        .iter()
        .map(|v| {
            let s = v
                .as_string()
                .ok_or("Invalid cached reference descriptor.")?;
            let record: DatasetRecord = serde_json::from_str(&s).map_err(|e| e.to_string())?;
            if !valid_cache(&record.cache_id) {
                return Err("Invalid cached dataset identity.".into());
            }
            Ok(record)
        })
        .collect()
}
pub async fn manifest_blob(db: &IdbDatabase, cache: &str) -> Result<Blob, String> {
    let tx = db
        .transaction_with_str(MANIFESTS)
        .map_err(|_| "Cannot access reference manifest store.")?;
    let pending = tx
        .object_store(MANIFESTS)
        .map_err(|_| "Cannot access reference manifest store.")?
        .get(&JsValue::from_str(cache))
        .map_err(|_| "Cannot request reference manifest.")?;
    let value = request(&pending).await?;
    if value.is_undefined() {
        return Err("The linked genome cache is unavailable. Reattach the hash-matching .cgref file; workspace links do not include bulk data.".into());
    }
    value
        .dyn_into::<Blob>()
        .map_err(|_| "Invalid cached manifest Blob.".into())
}
pub async fn manifest_bytes(db: &IdbDatabase, cache: &str) -> Result<Vec<u8>, String> {
    let value = manifest_blob(db, cache).await?;
    if value.size() > collagene_core::genome::MAX_MANIFEST_BYTES as f64 {
        return Err("Cached manifest exceeds its bound.".into());
    }
    let mut bytes = Vec::with_capacity(value.size() as usize);
    let mut at = 0;
    while at < value.size() as u64 {
        let n = READ_CHUNK.min((value.size() as u64 - at) as usize);
        bytes.extend(slice(&value, at, n).await?);
        at += n as u64;
    }
    Ok(bytes)
}
pub async fn remove(db: &IdbDatabase, cache: &str) -> Result<(), String> {
    if !valid_cache(cache) {
        return Err("Invalid reference cache identity.".into());
    }
    let low = JsValue::from_str(&format!("{cache}/"));
    let high = JsValue::from_str(&format!("{cache}0"));
    let range = IdbKeyRange::bound(&low, &high)
        .map_err(|_| "Cannot define reference block cleanup range.")?;
    let tx = stores(db, &[CATALOG, MANIFESTS, STAGING, BLOCKS])?;
    let queued = (|| -> Result<(), String> {
        for store in [CATALOG, MANIFESTS, STAGING] {
            tx.object_store(store)
                .map_err(|_| "Cannot access reference cleanup store.")?
                .delete(&JsValue::from_str(cache))
                .map_err(|_| "Cannot remove reference descriptor.")?;
        }
        tx.object_store(BLOCKS)
            .map_err(|_| "Cannot access reference block cleanup.")?
            .delete(range.as_ref())
            .map_err(|_| "Cannot remove staged reference blocks.")?;
        Ok(())
    })();
    if let Err(e) = queued {
        return Err(queue_error(&tx, &e));
    }
    committed(&tx).await
}
pub fn ensure_alive(alive: &Arc<AtomicBool>) -> Result<(), String> {
    if alive.load(Ordering::Relaxed) {
        Ok(())
    } else {
        Err("Reference operation cancelled; staging data will be removed.".into())
    }
}
pub fn download_blob(name: &str, value: &Blob) -> Result<(), String> {
    let url = web_sys::Url::create_object_url_with_blob(value)
        .map_err(|_| "Cannot create reference download URL.")?;
    let document = web_sys::window()
        .and_then(|w| w.document())
        .ok_or("No browser document.")?;
    let a = document
        .create_element("a")
        .map_err(|_| "Cannot create reference download.")?;
    a.set_attribute("href", &url)
        .map_err(|_| "Cannot set reference download URL.")?;
    a.set_attribute("download", name)
        .map_err(|_| "Cannot name reference download.")?;
    a.unchecked_ref::<web_sys::HtmlElement>().click();
    let callback = Closure::once_into_js(move || {
        let _ = web_sys::Url::revoke_object_url(&url);
    });
    if let Some(window) = web_sys::window() {
        let _ = window
            .set_timeout_with_callback_and_timeout_and_arguments_0(callback.unchecked_ref(), 1000);
    }
    Ok(())
}
