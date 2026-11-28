//! File-menu documentary images: originals, private worker previews, exact-byte downloads.
use super::{
    download::ExportArtifact,
    image_evidence_model::{Context, FileCapture, Install, Snapshot},
};
use crate::platform::worker;
use collagene_core::{
    documents::image_evidence::{self, Format, ImageEvidence, Prepared, Preview, Summary, Warning},
    protocol::{
        request::Request,
        schema::{AnalysisJob, AnalysisResult},
    },
};
use leptos::prelude::*;
use std::sync::{
    Arc,
    atomic::{AtomicBool, Ordering},
};
use std::{future::Future, pin::Pin};
use wasm_bindgen::JsCast;
use wasm_bindgen_futures::{JsFuture, spawn_local};
use web_sys::{Blob, BlobPropertyBag, HtmlInputElement, Url};

pub type AttachResult = Pin<Box<dyn Future<Output = Result<Summary, String>>>>;

const STYLE: &str = ".image-evidence-dialog{display:flex;flex-direction:column;gap:12px}.image-evidence-dialog label{display:flex;flex-direction:column;gap:5px}.image-evidence-dialog .image-evidence-actions{display:flex;gap:8px;flex-wrap:wrap}.image-evidence-records{max-height:22rem;overflow:auto;display:flex;flex-direction:column;gap:8px}.image-evidence-record{border:1px solid #dce4df;border-radius:5px;padding:10px;overflow-wrap:anywhere}.image-evidence-record p{margin:4px 0}.image-evidence-preview{margin:0;padding:10px;border:1px solid #dce4df;border-radius:5px;background:#f7faf7}.image-evidence-preview img{display:block;max-width:100%;max-height:55vh;object-fit:contain;margin:auto}.image-evidence-preview figcaption{padding-top:8px;overflow-wrap:anywhere}.image-evidence-identity{font-size:11px;overflow-wrap:anywhere}.image-evidence-warnings{padding-left:18px;font-size:12px;line-height:1.5}";
#[derive(Clone)]
struct SelectedFile {
    filename: String,
    hash: String,
    bytes: Arc<Vec<u8>>,
}
struct DisplayPreview {
    summary: Summary,
    url: String,
}
impl DisplayPreview {
    fn new(preview: Preview, expected: &Summary) -> Result<Self, String> {
        preview
            .validate_binding(expected)
            .map_err(|e| e.to_string())?;
        let parts = js_sys::Array::new();
        let bytes = js_sys::Uint8Array::from(preview.png.as_slice());
        parts.push(&bytes.buffer());
        let options = BlobPropertyBag::new();
        options.set_type("image/png");
        let blob = Blob::new_with_u8_array_sequence_and_options(&parts, &options)
            .map_err(|_| "Could not create the bounded image preview.")?;
        let url = Url::create_object_url_with_blob(&blob)
            .map_err(|_| "Could not open the bounded image preview.")?;
        Ok(Self {
            summary: preview.summary,
            url,
        })
    }
}
impl Drop for DisplayPreview {
    fn drop(&mut self) {
        let _ = Url::revoke_object_url(&self.url);
    }
}
fn format_label(format: Format) -> &'static str {
    match format {
        Format::Jpeg => "JPEG",
        Format::Tiff => "TIFF",
        Format::Png => "PNG",
    }
}
fn warning_label(warning: Warning) -> &'static str {
    match warning {
        Warning::DocumentaryPreview => {
            "This image is documentary evidence; it does not verify a cloning result or establish physical material."
        }
        Warning::MetadataNotTransferred => {
            "The preview omits source metadata. Exact original bytes remain available for download."
        }
        Warning::UncalibratedDisplay => {
            "Preview color and intensity are uncalibrated display values, not quantitative measurements."
        }
        Warning::AdditionalImagesNotPreviewed => {
            "Only the first image/page is previewed. Additional images remain in the original file."
        }
        Warning::PrecisionReducedForPreview => {
            "The preview has reduced sample precision. Original sample values remain in the source file."
        }
        Warning::DownsampledForPreview => {
            "The preview was downsampled to the display budget; the original dimensions are retained."
        }
    }
}
#[component]
fn Warnings(summary: Summary) -> impl IntoView {
    view! {<ul class="image-evidence-warnings">{summary.record.warnings.into_iter().map(|warning|view!{<li>{warning_label(warning)}</li>}).collect_view()}</ul>}
}

#[component]
pub fn ImageEvidenceDialog(
    snapshot: Callback<(), Result<Snapshot, String>>,
    on_attach: Callback<Install, AttachResult>,
    on_saved: Arc<dyn Fn(String) + Send + Sync>,
) -> impl IntoView {
    let initial = (|| {
        let source = snapshot.run(())?;
        let records = ImageEvidence::new(&source.workspace)
            .and_then(|images| images.list())
            .map_err(|e| e.to_string())?;
        let records = records
            .into_iter()
            .filter(|record| {
                record.record.association.document_id == source.context.association.document_id
            })
            .collect::<Vec<_>>();
        Ok::<_, String>((source.context, records))
    })();
    let opening = StoredValue::new(
        initial
            .as_ref()
            .map(|(context, _)| context.clone())
            .map_err(Clone::clone),
    );
    let records = StoredValue::new(initial.map(|(_, records)| records).unwrap_or_default());
    let selected = RwSignal::new(None::<Arc<SelectedFile>>);
    let name = RwSignal::new(String::new());
    let prepared = RwSignal::new(None::<Arc<Prepared>>);
    let preview = RwSignal::new(None::<Arc<DisplayPreview>>);
    let busy = RwSignal::new(false);
    let generation = RwSignal::new(0u32);
    let status = RwSignal::new(opening.get_value().err().unwrap_or_default());
    let alive = Arc::new(AtomicBool::new(true));
    let disposed = alive.clone();
    on_cleanup(move || disposed.store(false, Ordering::SeqCst));
    let current = Callback::new(move |()| -> Result<Snapshot, String> {
        let current = snapshot.run(())?;
        opening.get_value()?.matches(&current.context)?;
        Ok(current)
    });
    let choose = {
        let alive = alive.clone();
        move |event: web_sys::Event| {
            if busy.get_untracked() {
                return;
            }
            let input = event
                .target()
                .and_then(|node| node.dyn_into::<HtmlInputElement>().ok());
            let Some(input) = input else {
                return;
            };
            let files = input.files();
            generation.update(|value| *value = value.wrapping_add(1));
            selected.set(None);
            prepared.set(None);
            preview.set(None);
            name.set(String::new());
            status.set(String::new());
            let Some(files) = files else {
                return;
            };
            if files.length() != 1 {
                status.set("Choose exactly one JPEG, PNG or TIFF image.".into());
                return;
            }
            let Some(file) = files.get(0) else {
                return;
            };
            // Retain the File before clearing its browser-owned selection list.
            input.set_value("");
            let context = match current.run(()) {
                Ok(source) => source.context,
                Err(error) => {
                    status.set(error);
                    return;
                }
            };
            if !file.size().is_finite()
                || file.size() <= 0.
                || file.size() > image_evidence::MAX_SOURCE_BYTES as f64
            {
                status.set("Choose a nonempty image of at most 16 MiB.".into());
                return;
            }
            let expected_size = file.size() as usize;
            let filename = file.name();
            let token = generation.get_untracked().wrapping_add(1);
            generation.set(token);
            busy.set(true);
            selected.set(None);
            prepared.set(None);
            preview.set(None);
            status.set("Reading the explicitly selected image…".into());
            let alive = alive.clone();
            spawn_local(async move {
                let result = JsFuture::from(file.array_buffer()).await;
                if !alive.load(Ordering::SeqCst) || generation.try_get_untracked() != Some(token) {
                    return;
                }
                busy.set(false);
                let result = (|| {
                    context.matches(&current.run(())?.context)?;
                    let result =
                        result.map_err(|_| "Could not read the selected image.".to_string())?;
                    let array = js_sys::Uint8Array::new(&result);
                    if array.length() as usize != expected_size
                        || array.length() as usize > image_evidence::MAX_SOURCE_BYTES
                    {
                        return Err("Image bytes changed or exceeded their size bound.".into());
                    }
                    let bytes = array.to_vec();
                    let hash = blake3::hash(&bytes).to_hex().to_string();
                    Ok::<_, String>(SelectedFile {
                        filename: filename.clone(),
                        hash,
                        bytes: Arc::new(bytes),
                    })
                })();
                match result {
                    Ok(source) => {
                        name.set(filename);
                        selected.set(Some(Arc::new(source)));
                        status.set("Review the image name, then prepare its preview.".into());
                    }
                    Err(error) => status.set(error),
                }
            });
        }
    };
    let prepare = {
        let alive = alive.clone();
        move |_| {
            if busy.get_untracked() {
                return;
            }
            let Some(file) = selected.get_untracked() else {
                status.set("Choose an image file first.".into());
                return;
            };
            let capture = match current.run(()).and_then(|source| {
                FileCapture::new(
                    source.context,
                    name.get_untracked(),
                    file.hash.clone(),
                    file.bytes.len(),
                )
            }) {
                Ok(capture) => capture,
                Err(error) => {
                    status.set(error);
                    return;
                }
            };
            let token = generation.get_untracked().wrapping_add(1);
            generation.set(token);
            busy.set(true);
            prepared.set(None);
            preview.set(None);
            status.set("Preparing a bounded preview in the private worker…".into());
            let alive = alive.clone();
            spawn_local(async move {
                let result = worker::analyze_cancellable(
                    Request {
                        id: token,
                        revision: capture.revision.clone(),
                        job: AnalysisJob::PrepareImageEvidence {
                            name: capture.name.clone(),
                            association: capture.context.association.clone(),
                            original: file.bytes.as_ref().clone(),
                        },
                    },
                    alive.clone(),
                )
                .await;
                if !alive.load(Ordering::SeqCst) || generation.try_get_untracked() != Some(token) {
                    return;
                }
                busy.set(false);
                let result = (|| {
                    capture.context.matches(&current.run(())?.context)?;
                    if name.get_untracked() != capture.name
                        || !selected.with_untracked(|source| {
                            source.as_ref().is_some_and(|source| {
                                source.hash == capture.original_hash
                                    && source.bytes.len() == capture.source_bytes
                            })
                        })
                    {
                        return Err("The selected image or name changed; prepare it again.".into());
                    }
                    let response = result?;
                    if response.id != token || response.revision != capture.revision {
                        return Err("Image worker response belongs to another capture.".into());
                    }
                    let packet = match response.result {
                        AnalysisResult::ImageEvidencePrepared(packet) => packet,
                        AnalysisResult::Error(error) => return Err(error),
                        _ => return Err("Unexpected image preparation response.".into()),
                    };
                    let summary = capture.admit(&packet)?;
                    let display = DisplayPreview::new(
                        packet.preview().map_err(|e| e.to_string())?,
                        &summary,
                    )?;
                    Ok::<_, String>((Arc::new(packet), Arc::new(display)))
                })();
                match result {
                    Ok((packet, display)) => {
                        prepared.set(Some(packet));
                        preview.set(Some(display));
                        status.set("Preview prepared. Attach only if this is the intended image and document.".into());
                    }
                    Err(error) => status.set(error),
                }
            });
        }
    };
    let attach = {
        let alive = alive.clone();
        let on_saved = on_saved.clone();
        move |_| {
            if busy.get_untracked() {
                return;
            }
            let Some(packet) = prepared.get_untracked() else {
                return;
            };
            let context = match current.run(()) {
                Ok(source) => source.context,
                Err(error) => {
                    status.set(error);
                    return;
                }
            };
            let document_name = context.document_name.clone();
            let token = generation.get_untracked().wrapping_add(1);
            generation.set(token);
            busy.set(true);
            status.set("Checking that the complete native workspace can retain this image…".into());
            let operation = on_attach.run(Install {
                context,
                prepared: packet,
                alive: alive.clone(),
            });
            let alive = alive.clone();
            let on_saved = on_saved.clone();
            spawn_local(async move {
                let result = operation.await;
                if let Ok(summary) = &result {
                    // Adoption closes this dialog. The owned notice and captured label remain
                    // valid after disposal; do not read the child arena on durable success.
                    on_saved(format!(
                        "Attached image '{}' to '{}'. Native workspace saveability was checked; local autosave will retain its original and preview.",
                        summary.record.name, document_name
                    ));
                }
                if !alive.load(Ordering::SeqCst) || generation.try_get_untracked() != Some(token) {
                    return;
                }
                busy.set(false);
                if let Err(error) = result {
                    status.set(error);
                }
            });
        }
    };
    let reopen = {
        let alive = alive.clone();
        Callback::new(move |summary: Summary| {
            if busy.get_untracked() {
                return;
            }
            prepared.set(None);
            preview.set(None);
            let source = match current.run(()) {
                Ok(source) => source,
                Err(error) => {
                    status.set(error);
                    return;
                }
            };
            let result = (|| {
                let packet = ImageEvidence::new(&source.workspace)
                    .and_then(|images| images.prepared(&summary.manifest_hash))
                    .map_err(|e| e.to_string())?;
                let retained = packet.summary().map_err(|e| e.to_string())?;
                if retained.manifest_hash != summary.manifest_hash
                    || retained.record != summary.record
                {
                    return Err("The selected image record changed.".into());
                }
                let revision = summary.reopen_revision().map_err(|e| e.to_string())?;
                Ok::<_, String>((packet, revision))
            })();
            let (packet, revision) = match result {
                Ok(value) => value,
                Err(error) => {
                    status.set(error);
                    return;
                }
            };
            let context = source.context;
            let token = generation.get_untracked().wrapping_add(1);
            generation.set(token);
            busy.set(true);
            prepared.set(None);
            preview.set(None);
            status.set("Reopening the exact retained original in the private worker…".into());
            let alive = alive.clone();
            spawn_local(async move {
                let result = worker::analyze_cancellable(
                    Request {
                        id: token,
                        revision: revision.clone(),
                        job: AnalysisJob::ReopenImageEvidence { prepared: packet },
                    },
                    alive.clone(),
                )
                .await;
                if !alive.load(Ordering::SeqCst) || generation.try_get_untracked() != Some(token) {
                    return;
                }
                busy.set(false);
                let result = (|| {
                    context.matches(&current.run(())?.context)?;
                    let response = result?;
                    if response.id != token || response.revision != revision {
                        return Err("Image preview response belongs to another record.".into());
                    }
                    match response.result {
                        AnalysisResult::ImageEvidencePreview(preview) => {
                            DisplayPreview::new(preview, &summary)
                        }
                        AnalysisResult::Error(error) => Err(error),
                        _ => Err("Unexpected image preview response.".into()),
                    }
                })();
                match result {
                    Ok(display) => {
                        preview.set(Some(Arc::new(display)));
                        status.set("Retained original decoded and preview checked. Original bytes remain unchanged.".into());
                    }
                    Err(error) => status.set(error),
                }
            });
        })
    };
    let original = {
        let alive = alive.clone();
        let on_saved = on_saved.clone();
        Callback::new(move |summary: Summary| {
            if busy.get_untracked() {
                return;
            }
            let source = match current.run(()) {
                Ok(source) => source,
                Err(error) => {
                    status.set(error);
                    return;
                }
            };
            let result = (|| {
                let images = ImageEvidence::new(&source.workspace).map_err(|e| e.to_string())?;
                let original = images
                    .original(&summary.manifest_hash)
                    .map_err(|e| e.to_string())?;
                if original.content_hash != summary.record.original.hash
                    || original.bytes.len() as u64 != summary.record.original.bytes
                {
                    return Err("Retained original differs from the selected record.".into());
                }
                ExportArtifact::prepare(&original.name, original.media_type, original.bytes)
            })();
            let artifact = match result {
                Ok(artifact) => artifact,
                Err(error) => {
                    status.set(error);
                    return;
                }
            };
            let context = source.context;
            let token = generation.get_untracked().wrapping_add(1);
            generation.set(token);
            busy.set(true);
            let saving = alive.clone();
            let complete = alive.clone();
            let on_saved = on_saved.clone();
            let label = summary.record.name;
            artifact.complete_guarded(move||{if !saving.load(Ordering::SeqCst)||generation.try_get_untracked()!=Some(token){return Err("The image evidence view closed before save admission.".into());}context.matches(&current.run(())?.context)},move|result|{
            if result.is_ok(){on_saved(format!("Downloaded exact original image '{label}'."));}
            if !complete.load(Ordering::SeqCst)||busy.try_get_untracked().is_none(){return;}busy.set(false);status.set(match result{Ok(())=>"Exact original bytes downloaded; preview transformations were not applied.".into(),Err(error)=>error});
        });
        })
    };
    view! {<section class="image-evidence-dialog"><style>{STYLE}</style>
        <p>{move||opening.get_value().map(|context|format!("Images for {} · revision {}",context.document_name,&context.association.revision[..12])).unwrap_or_else(|error|error)}</p>
        <p class="helper">"Documentary images are saved with the native workspace (.cg2design), not standalone .cg2, GenBank or FASTA sequence exports. Originals remain exact; previews are uncalibrated display images. This does not compare observed and predicted gels."</p>
        <label>"Choose JPEG, PNG or TIFF (maximum 16 MiB)"<input type="file" accept=".jpg,.jpeg,.png,.tif,.tiff,image/jpeg,image/png,image/tiff" disabled=move||busy.get() on:change=choose/></label>
        {move||selected.get().map(|source|view!{<p class="helper">{format!("{} · {} bytes",source.filename,source.bytes.len())}</p>})}
        <label>"Image name"<input type="text" maxlength="256" prop:value=move||name.get() disabled=move||busy.get() on:input=move|event|{name.set(event_target_value(&event));prepared.set(None);preview.set(None);}/></label>
        <div class="image-evidence-actions"><button class="button secondary" disabled=move||busy.get()||selected.with(Option::is_none) on:click=prepare>"Prepare preview"</button><button class="button primary" disabled=move||busy.get()||prepared.with(Option::is_none) on:click=attach>"Attach image"</button></div>
        {move||preview.get().map(|display|{let summary=display.summary.clone();view!{<figure class="image-evidence-preview"><img src=display.url.clone() alt=format!("Documentary preview of {}",summary.record.name)/><figcaption>{format!("{} · {} · original {}×{} · preview {}×{}",summary.record.name,format_label(summary.record.format),summary.record.original.width,summary.record.original.height,summary.record.preview.width,summary.record.preview.height)}</figcaption><Warnings summary=summary.clone()/><details class="image-evidence-identity"><summary>"Source identity"</summary><p>{format!("Original {} · associated document {} / revision {}",summary.record.original.hash,summary.record.association.document_id,summary.record.association.revision)}</p></details></figure>}})}
        <h3>"Retained images"</h3>
        <p class="helper">"These records keep their explicitly associated revision, including older retained revisions of this document."</p>
        <div class="image-evidence-records">{records.get_value().into_iter().map(|summary|{let open=summary.clone();let download=summary.clone();view!{<article class="image-evidence-record"><strong>{summary.record.name.clone()}</strong><p>{format!("{} · {}×{} · {} bytes · associated revision {}",format_label(summary.record.format),summary.record.original.width,summary.record.original.height,summary.record.original.bytes,&summary.record.association.revision[..12])}</p><div class="image-evidence-actions"><button class="button secondary" disabled=move||busy.get() on:click=move|_|reopen.run(open.clone())>"Preview"</button><button class="button secondary" disabled=move||busy.get() on:click=move|_|original.run(download.clone())>"Download original"</button></div></article>}}).collect_view()}</div>
        <Show when=move||records.with_value(Vec::is_empty)><p class="helper">"No documentary images are retained for this document."</p></Show>
        <p role="status" aria-live="polite">{move||status.get()}</p>
    </section>}
}
