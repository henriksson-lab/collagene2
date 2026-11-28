//! Cancellable binary-worker exports from the exact retained comparison.
use crate::{pairwise_result_view::Review, worker::CachedWorker};
use collagene_core::{
    pairwise_export::{self, Export, Format},
    protocol::{
        request::Request,
        schema::{AnalysisJob, AnalysisResult},
    },
};
use leptos::prelude::*;
use std::{
    rc::Rc,
    sync::{
        Arc,
        atomic::{AtomicBool, Ordering},
    },
};
use wasm_bindgen_futures::spawn_local;

#[component]
pub fn PairwiseExportView(
    review: Arc<Review>,
    guard: Callback<Arc<Review>, Result<(), String>>,
) -> impl IntoView {
    let bound = StoredValue::new(review);
    let format = RwSignal::new(Format::AlignedFasta);
    let busy = RwSignal::new(false);
    let status = RwSignal::new(String::new());
    let prepared = RwSignal::new(None::<Arc<Export>>);
    let generation = RwSignal::new(0u32);
    let worker = StoredValue::new_local(None::<Rc<CachedWorker>>);
    let alive = Arc::new(AtomicBool::new(true));
    let cleanup_alive = alive.clone();
    on_cleanup(move || {
        cleanup_alive.store(false, Ordering::Relaxed);
        if let Some(w) = worker.get_value() {
            w.close();
        }
    });
    Effect::new(move || {
        let _ = format.get();
        generation.update(|v| *v = v.wrapping_add(1));
        if let Some(w) = worker.get_value() {
            w.close();
        }
        worker.set_value(None);
        busy.set(false);
        prepared.set(None);
        status.set(String::new());
    });
    let create = {
        let alive = alive.clone();
        move |_| {
            if busy.get_untracked() {
                return;
            }
            let input = (|| -> Result<_, String> {
                let review = bound.get_value();
                guard.run(review.clone())?;
                let chosen = format.get_untracked();
                let revision = pairwise_export::revision(
                    &review.record,
                    &review.reference_name,
                    &review.query_name,
                    review.fresh,
                    chosen,
                )
                .map_err(|e| e.to_string())?;
                let w = Rc::new(CachedWorker::new()?);
                Ok((review, chosen, revision, w))
            })();
            let (review, chosen, revision, w) = match input {
                Ok(value) => value,
                Err(error) => {
                    status.set(error);
                    return;
                }
            };
            let token = generation.get_untracked().wrapping_add(1);
            generation.set(token);
            worker.set_value(Some(w.clone()));
            busy.set(true);
            prepared.set(None);
            status.set("Preparing the complete comparison export in the Rust worker…".into());
            let alive = alive.clone();
            spawn_local(async move {
                let response = w
                    .ask(Request {
                        id: token,
                        revision: revision.clone(),
                        job: AnalysisJob::ExportPairwiseComparison {
                            record: review.record.clone(),
                            reference_name: review.reference_name.clone(),
                            query_name: review.query_name.clone(),
                            fresh: review.fresh,
                            format: chosen,
                        },
                    })
                    .await;
                w.close();
                if !alive.load(Ordering::Relaxed) || generation.get_untracked() != token {
                    return;
                }
                worker.set_value(None);
                busy.set(false);
                let outcome = (|| -> Result<Export, String> {
                    guard.run(review.clone())?;
                    let response = response?;
                    if response.id != token
                        || response.revision != revision
                        || format.get_untracked() != chosen
                    {
                        return Err("Export inputs changed; prepare the export again.".into());
                    }
                    let output = match response.result {
                        AnalysisResult::ExportedPairwiseComparison(output) => output,
                        AnalysisResult::Error(error) => return Err(error),
                        _ => return Err("Unexpected comparison-export response.".into()),
                    };
                    output.validate().map_err(|e| e.to_string())?;
                    let expected = match chosen {
                        Format::AlignedFasta => "fasta",
                        Format::Clustal => "aln",
                        Format::Html => "html",
                        Format::Pdf => "pdf",
                    };
                    if output.extension != expected {
                        return Err("Export format differs from the requested format.".into());
                    }
                    guard.run(review)?;
                    Ok(output)
                })();
                match outcome {
                    Ok(output) => {
                        status.set(
                            "Export prepared. Review the format limitations before downloading."
                                .into(),
                        );
                        prepared.set(Some(Arc::new(output)));
                    }
                    Err(error) => status.set(error),
                }
            });
        }
    };
    let save = move |_| {
        let result = (|| {
            guard.run(bound.get_value())?;
            let output = prepared.get_untracked().ok_or("Prepare an export first.")?;
            crate::app::prepare_download(
                &format!("pairwise-comparison.{}", output.extension),
                &output.mime,
                &output.bytes,
            )
        })();
        let captured_generation = generation.get_untracked();
        let captured_format = format.get_untracked();
        match result {
            Ok(artifact) => {
                let check = move || {
                    if status.try_get_untracked().is_none() {
                        return Err("Export review closed.".into());
                    }
                    if generation.get_untracked() != captured_generation
                        || format.get_untracked() != captured_format
                    {
                        return Err("Comparison export settings changed.".into());
                    }
                    guard.run(bound.get_value())
                };
                let completion = check.clone();
                artifact.complete_guarded(check, move |result| {
                    if completion().is_err() { return; }
                    status.set(match result { Ok(()) => "Comparison export saved. Historical/fresh labels and sequence differences do not establish experimental verification.".into(), Err(error) => error });
                });
            }
            Err(error) => status.set(error),
        }
    };
    view!{<details class="pairwise-export"><summary>"Export aligned sequences or printable review"</summary>
        <p class="helper">"Standard alignment files retain the displayed orientation, gaps and aligned segment. HTML/PDF retain the source-event and comparison context. Original source coordinates can differ from displayed coordinates. Complete native evidence remains available as JSON."</p>
        <label>"Export format"<select prop:value=move||match format.get(){Format::AlignedFasta=>"fasta",Format::Clustal=>"clustal",Format::Html=>"html",Format::Pdf=>"pdf"} on:change=move|event|format.set(match event_target_value(&event).as_str(){"clustal"=>Format::Clustal,"html"=>Format::Html,"pdf"=>Format::Pdf,_=>Format::AlignedFasta})><option value="fasta">"Aligned FASTA"</option><option value="clustal">"Clustal"</option><option value="html">"Printable HTML"</option><option value="pdf">"PDF review"</option></select></label>
        <button class="button secondary" disabled=move||busy.get() on:click=create>"Prepare comparison export"</button>
        <button class="text-tool" disabled=move||!busy.get() on:click=move|_|{generation.update(|v|*v=v.wrapping_add(1));if let Some(w)=worker.get_value(){w.close();}worker.set_value(None);busy.set(false);status.set("Export preparation cancelled. No download was created.".into());}>"Cancel export"</button>
        {move||prepared.get().map(|output|view!{<section><h4>"Export limitations"</h4><ul>{output.warnings.iter().map(|warning|view!{<li>{warning.clone()}</li>}).collect_view()}</ul><p>{format!("{} bytes · {}",output.bytes.len(),output.mime)}</p><button class="button primary" on:click=save>"Download prepared comparison"</button></section>})}
        <p class="helper" role="status">{move||status.get()}</p>
    </details>}.into_any()
}
