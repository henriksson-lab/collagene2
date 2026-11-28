//! Worker-backed, explicitly reviewed regeneration of the displayed actual trace.
use crate::platform::worker::CachedWorker;
use crate::reads::scf_export_model::Review;
use crate::reads::scf_export_model::Scope;
use crate::reads::scf_export_model::safe_filename;
use collagene_core::{
    formats::scf::{self as scf_export, LossPolicy, Options, SamplePrecision},
    protocol::{
        request::Request,
        schema::{AnalysisJob, AnalysisResult},
    },
    review::reads::{
        traces::QualityEvidence,
        traces::display::{DisplayTrace, PlotData},
    },
};
use collagene_model::Molecule;
use leptos::prelude::*;
use std::{
    rc::Rc,
    sync::{
        Arc,
        atomic::{AtomicBool, Ordering},
    },
};
use wasm_bindgen_futures::spawn_local;
struct Prepared {
    review: Review,
    report: scf_export::DocumentExport,
    name: String,
    fallback_name: bool,
}
fn scope(revision: &str, trace: &DisplayTrace, precision: SamplePrecision) -> Scope {
    Scope {
        document_revision: revision.into(),
        source_hash: trace.source_hash.clone(),
        orientation: trace.orientation,
        precision,
    }
}
#[component]
pub fn ScfExportControls(
    document_revision: Memo<String>,
    molecule: Memo<Molecule>,
    display: Memo<Option<Arc<PlotData>>>,
) -> impl IntoView {
    let precision = RwSignal::new(SamplePrecision::Smallest);
    let busy = RwSignal::new(false);
    let status = RwSignal::new(String::new());
    let acknowledged = RwSignal::new(false);
    let prepared = RwSignal::new(None::<Arc<Prepared>>);
    let generation = RwSignal::new(0u32);
    let worker = StoredValue::new_local(None::<Rc<CachedWorker>>);
    let alive = Arc::new(AtomicBool::new(true));
    let clean = alive.clone();
    on_cleanup(move || {
        clean.store(false, Ordering::Relaxed);
        if let Some(w) = worker.get_value() {
            w.close();
        }
    });
    let current = Memo::new(move |_| {
        display
            .get()
            .map(|t| scope(&document_revision.get(), &t.trace, precision.get()))
    });
    Effect::new(move |_| {
        let _ = current.get();
        generation.update(|n| *n = n.wrapping_add(1));
        if let Some(w) = worker.get_value() {
            w.close();
        }
        worker.set_value(None);
        prepared.set(None);
        acknowledged.set(false);
        busy.set(false);
        status.set(String::new());
    });
    let prepare = {
        let alive = alive.clone();
        move |_| {
            if busy.get_untracked() {
                return;
            }
            let Some(trace) = display.get_untracked() else {
                return;
            };
            let Some(capture) = current.get_untracked() else {
                return;
            };
            let source = molecule.get_untracked();
            if source.revision() != trace.source_revision {
                status.set("The displayed trace no longer matches the document.".into());
                return;
            }
            let source_revision = source.revision();
            let options = Options {
                precision: precision.get_untracked(),
                loss_policy: LossPolicy::AllowDocumentaryOmissions,
            };
            let context = match scf_export::document_revision(&source, &options) {
                Ok(c) => c,
                Err(e) => {
                    status.set(e.to_string());
                    return;
                }
            };
            let active = match CachedWorker::new() {
                Ok(w) => Rc::new(w),
                Err(e) => {
                    status.set(e);
                    return;
                }
            };
            let id = generation.get_untracked().wrapping_add(1);
            generation.set(id);
            worker.set_value(Some(active.clone()));
            prepared.set(None);
            acknowledged.set(false);
            busy.set(true);
            status.set("Preparing SCF from the retained trace evidence…".into());
            let alive = alive.clone();
            spawn_local(async move {
                let response = active
                    .ask(Request {
                        id,
                        revision: context.clone(),
                        job: AnalysisJob::ExportScfTrace {
                            molecule: source,
                            options,
                        },
                    })
                    .await;
                active.close();
                if !alive.load(Ordering::Relaxed)
                    || generation.get_untracked() != id
                    || current.get_untracked().as_ref() != Some(&capture)
                {
                    return;
                }
                worker.set_value(None);
                busy.set(false);
                let result = (|| {
                    let response = response?;
                    if response.id != id || response.revision != context {
                        return Err("SCF worker response identity changed.".into());
                    }
                    let report = match response.result {
                        AnalysisResult::ScfTraceExported(r) => r,
                        AnalysisResult::Error(e) => return Err(e),
                        _ => return Err("Unexpected SCF export response.".into()),
                    };
                    report
                        .validate_binding(&context)
                        .map_err(|e| e.to_string())?;
                    if report.source_revision != source_revision
                        || report.report.source_hash != capture.source_hash
                        || report.report.orientation != capture.orientation
                    {
                        return Err("SCF response names another trace/orientation.".into());
                    }
                    if !matches!(
                        (capture.precision, report.report.sample_size),
                        (SamplePrecision::Smallest, 1 | 2)
                            | (SamplePrecision::EightBit, 1)
                            | (SamplePrecision::SixteenBit, 2)
                    ) {
                        return Err("SCF response sample precision changed.".into());
                    }
                    let review = Review::new(capture, context, !report.report.losses.is_empty());
                    let name = safe_filename(
                        &trace.trace.name,
                        &trace.trace.source_hash,
                        trace.trace.orientation,
                    );
                    let fallback_name = name != format!("{}.scf", trace.trace.name);
                    Ok::<_, String>(Prepared {
                        review,
                        report,
                        name,
                        fallback_name,
                    })
                })();
                match result {
                    Ok(value) => {
                        status.set(if value.report.report.losses.is_empty(){"SCF prepared. Review the summary, then save."}else{"Review the documentary omissions and acknowledge them before saving."}.into());
                        prepared.set(Some(Arc::new(value)));
                    }
                    Err(e) => status.set(e),
                }
            });
        }
    };
    let save = {
        let alive = alive.clone();
        move |_| {
            if busy.get_untracked() {
                return;
            }
            let Some(value) = prepared.get_untracked() else {
                return;
            };
            let Some(now) = current.get_untracked() else {
                return;
            };
            if let Err(e) = value.review.permit(&now, acknowledged.get_untracked()) {
                status.set(e);
                return;
            }
            if let Err(e) = value.report.validate_binding(value.review.context()) {
                status.set(e.to_string());
                return;
            }
            let artifact = match crate::io::download::prepare_download(
                &value.name,
                "application/octet-stream",
                &value.report.report.bytes,
            ) {
                Ok(v) => v,
                Err(e) => {
                    status.set(e);
                    return;
                }
            };
            let token = generation.get_untracked();
            busy.set(true);
            status.set("Saving the reviewed SCF export…".into());
            let valid_value = value.clone();
            let valid_alive = alive.clone();
            let done_alive = alive.clone();
            artifact.complete_guarded(
                move || {
                    if !valid_alive.load(Ordering::Relaxed)
                        || generation.try_get_untracked() != Some(token)
                    {
                        return Err("SCF export panel changed or closed.".into());
                    }
                    let current = current
                        .try_get_untracked()
                        .flatten()
                        .ok_or("Trace is no longer available.")?;
                    valid_value
                        .review
                        .permit(&current, acknowledged.try_get_untracked().unwrap_or(false))
                },
                move |result| {
                    if !done_alive.load(Ordering::Relaxed)
                        || generation.try_get_untracked() != Some(token)
                    {
                        return;
                    }
                    let Some(now) = current.try_get_untracked().flatten() else {
                        return;
                    };
                    if value
                        .review
                        .permit(&now, acknowledged.try_get_untracked().unwrap_or(false))
                        .is_err()
                    {
                        return;
                    }
                    busy.set(false);
                    status.set(match result {
                        Ok(()) => "Reviewed SCF export saved.".into(),
                        Err(e) => e,
                    });
                },
            );
        }
    };
    let cancel = move |_| {
        generation.update(|n| *n = n.wrapping_add(1));
        if let Some(w) = worker.get_value() {
            w.close();
        }
        worker.set_value(None);
        busy.set(false);
        prepared.set(None);
        acknowledged.set(false);
        status.set("SCF preparation cancelled.".into());
    };
    view! {<details class="trace-export"><summary>"Export trace…"</summary>
        <div class="trace-export-content"><p>"Regenerate SCF 3.00 from the displayed actual signal, calls, peaks and four stored accuracy vectors. Original attachment export remains separate. No synthetic trace or quality conversion."</p>
        <label>"Sample precision "<select disabled=move||busy.get() on:change=move|e|precision.set(match event_target_value(&e).as_str(){"8"=>SamplePrecision::EightBit,"16"=>SamplePrecision::SixteenBit,_=>SamplePrecision::Smallest})><option value="auto">"Smallest lossless (8 or 16 bits)"</option><option value="8">"8 bits (refuse larger signal)"</option><option value="16">"16 bits"</option></select></label>
        <Show when=move||display.get().is_some_and(|t|!matches!(&t.trace.qualities,Some(QualityEvidence::ScfAccuracy(_))))><p>"This trace has no four-channel SCF accuracy evidence. ABIF scores or missing qualities cannot be converted by this strict profile."</p></Show>
        <div class="trace-export-actions"><button disabled=move||busy.get()||display.get().is_none() on:click=prepare>"Prepare export review"</button><Show when=move||busy.get()&&prepared.get().is_none()><button on:click=cancel>"Cancel"</button></Show></div>
        {move||prepared.get().map(|value|{let report=&value.report.report;let filename = format!("Filename: {}{}",value.name,if value.fallback_name{". A portable basename was selected; the SCF sample NAME and original attachment are unchanged."}else{""});let summary=format!("{} bytes · {}-bit samples · {:?} orientation · source {}",report.bytes.len(),report.sample_size*8,report.orientation,report.source_hash);let losses=report.losses.iter().map(|l|format!("{} ({} bytes): {}",l.field,l.bytes,l.detail)).collect::<Vec<_>>();let has_losses=!losses.is_empty();let warnings=report.warnings.clone();view!{<div class="trace-export-review"><p>{filename}</p><p>{summary}</p><ul>{warnings.into_iter().map(|w|view!{<li>{w}</li>}).collect_view()}</ul><Show when=move||has_losses><ul>{losses.iter().cloned().map(|l|view!{<li>{l}</li>}).collect_view()}</ul><label><input type="checkbox" prop:checked=move||acknowledged.get() disabled=move||busy.get() on:change=move|e|acknowledged.set(event_target_checked(&e))/>"I accept these documentary omissions; retain the original attachment separately."</label></Show><button disabled=move||busy.get()||(has_losses&&!acknowledged.get()) on:click=save.clone()>"Save reviewed SCF…"</button></div>}})}
        <p role="status" aria-live="polite">{move||status.get()}</p></div>
    </details>}
}
