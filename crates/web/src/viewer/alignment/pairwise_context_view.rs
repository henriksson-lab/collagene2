//! Optional, cached source annotation/CDS context for computed pairwise results.
use crate::platform::worker::CachedWorker;
use crate::viewer::alignment::pairwise_result_view::Review;
use collagene_core::{
    alignment::context::{self as pairwise_context, Context, Page},
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

const TRACK_PAGE: usize = 20;

#[component]
pub fn PairwiseContextView(
    review: Arc<Review>,
    page: RwSignal<usize>,
    guard: Callback<Arc<Review>, Result<(), String>>,
) -> impl IntoView {
    let review = StoredValue::new(review);
    let annotations = RwSignal::new(true);
    let translations = RwSignal::new(false);
    let prepared = RwSignal::new(None::<Result<Arc<Context>, String>>);
    let tracks_page = RwSignal::new(0usize);
    let warnings_page = RwSignal::new(0usize);
    let busy = RwSignal::new(false);
    let status = RwSignal::new(String::new());
    let generation = RwSignal::new(0u32);
    let active = StoredValue::new_local(None::<Rc<CachedWorker>>);
    let alive = Arc::new(AtomicBool::new(true));
    let cleanup_alive = alive.clone();
    on_cleanup(move || {
        cleanup_alive.store(false, Ordering::Relaxed);
        if let Some(worker) = active.get_value() {
            worker.close();
        }
    });
    Effect::new(move |_| {
        let _ = translations.get();
        generation.update(|n| *n = n.wrapping_add(1));
        if let Some(worker) = active.get_value() {
            worker.close();
        }
        active.set_value(None);
        busy.set(false);
        prepared.set(None);
        status.set(String::new());
        tracks_page.set(0);
        warnings_page.set(0);
    });
    Effect::new(move |_| {
        let _ = page.get();
        tracks_page.set(0);
    });
    let load = {
        let alive = alive.clone();
        move |_| {
            if busy.get_untracked() {
                return;
            }
            let input = (|| -> Result<_, String> {
                let captured = review.get_value();
                guard.run(captured.clone())?;
                let translate = translations.get_untracked();
                let revision = pairwise_context::revision(
                    &captured.record.report,
                    &captured.reference,
                    &captured.query,
                    translate,
                )
                .map_err(|e| e.to_string())?;
                let worker = Rc::new(CachedWorker::new()?);
                Ok((captured, translate, revision, worker))
            })();
            let (captured, translate, revision, worker) = match input {
                Ok(value) => value,
                Err(error) => {
                    status.set(error);
                    return;
                }
            };
            let token = generation.get_untracked().wrapping_add(1);
            generation.set(token);
            active.set_value(Some(worker.clone()));
            busy.set(true);
            prepared.set(None);
            status.set("Preparing annotation/CDS context in the Rust worker…".into());
            let alive = alive.clone();
            spawn_local(async move {
                let response = worker
                    .ask(Request {
                        id: token,
                        revision: revision.clone(),
                        job: AnalysisJob::PreparePairwiseContext {
                            report: captured.record.report.clone(),
                            reference: (*captured.reference).clone(),
                            query: (*captured.query).clone(),
                            translate,
                        },
                    })
                    .await;
                worker.close();
                if !alive.load(Ordering::Relaxed) || generation.get_untracked() != token {
                    return;
                }
                active.set_value(None);
                busy.set(false);
                let outcome = (|| -> Result<Arc<Context>, String> {
                    guard.run(captured.clone())?;
                    let response = response?;
                    if response.id != token
                        || response.revision != revision
                        || translations.get_untracked() != translate
                    {
                        return Err("Context inputs changed; prepare the context again.".into());
                    }
                    let context = match response.result {
                        AnalysisResult::PreparedPairwiseContext(context) => context,
                        AnalysisResult::Error(error) => return Err(error),
                        _ => return Err("Unexpected pairwise context response.".into()),
                    };
                    context.validate().map_err(|e| e.to_string())?;
                    if context.context_revision != revision || context.translated != translate {
                        return Err(
                            "Context receipt differs from the requested comparison/options.".into(),
                        );
                    }
                    guard.run(captured)?;
                    Ok(Arc::new(context))
                })();
                prepared.set(Some(outcome));
                tracks_page.set(0);
                warnings_page.set(0);
                status.set("Context request finished. These source annotations and translations are reference context, not experimental verification.".into());
            });
        }
    };
    let projected = Memo::new(move |_| -> Option<Result<Arc<Page>, String>> {
        prepared.get().map(|result| {
            result.and_then(|context| {
                review.with_value(|review| {
                    context
                        .page(&review.record.report, page.get())
                        .map(Arc::new)
                })
            })
        })
    });
    view!{<details class="pairwise-context"><summary>"Annotation and CDS context"</summary>
        <style>".pairwise-context .pairwise-context-scroll{overflow-x:auto;max-width:100%}.pairwise-context .pairwise-context-scroll pre{white-space:pre;font:12px/1.6 monospace;max-height:none;overflow:visible;margin:6px 0}.pairwise-context .pairwise-context-track{border-top:1px solid #e1e7e0;padding:8px 0}.pairwise-context label{display:block;margin:6px 0}"</style>
        <p class="helper">"Context follows this alignment's current column page and exact retained source versions. It does not edit annotations or infer gene function. Gaps have no original source base. Source feature indices stay stable."</p>
        <label><input type="checkbox" prop:checked=move||annotations.get() on:change=move|e|annotations.set(event_target_checked(&e))/>"Show source annotation directions"</label>
        <label><input type="checkbox" prop:checked=move||translations.get() on:change=move|e|translations.set(event_target_checked(&e))/>"Include native CDS translations"</label>
        <button class="button secondary" disabled=move||busy.get() on:click=load>"Load aligned source context"</button>
        <button class="text-tool" disabled=move||!busy.get() on:click=move|_|{generation.update(|n|*n=n.wrapping_add(1));if let Some(worker)=active.get_value(){worker.close();}active.set_value(None);busy.set(false);prepared.set(None);status.set("Context preparation cancelled.".into());}>"Cancel context preparation"</button>
        <p class="helper" role="status">{move||status.get()}</p>
        <p class="helper">"Translation uses the complete original CDS, its genetic code and codon_start. It never restarts at a page, gap or circular origin. Residues occupy their middle codon base; dots mark the other codon calls. Read reverse CDS in biological direction. '*' includes terminal stop context; '!' marks conflicting repeated mappings."</p>
        <p class="helper">"Translation preparation is bounded to 1,500,000 aggregate source positions and 64,000,000 bases of source/CDS preparation work. A withheld CDS retains its annotation and an explicit reason. Context is prepared once; column navigation reuses cached mappings."</p>
        {move||match projected.get() {
            None=>view!{<p class="helper">"Load context after choosing whether to include CDS translations."</p>}.into_any(),
            Some(Err(error))=>view!{<p role="alert">{error}</p>}.into_any(),
            Some(Ok(current))=>{
                let pages=current.tracks.len().div_ceil(TRACK_PAGE).max(1);
                let current_page=tracks_page.get().min(pages-1);
                let labels=format!("Columns {}–{} · {} overlapping source features · feature page {} of {pages}",if current.end==0{0}else{current.start+1},current.end,current.tracks.len(),current_page+1);
                let base_rows=format!("R {}\nQ {}",current.reference,current.query);
                let position_label=|positions:&[Option<usize>]|{
                    let first=positions.iter().flatten().next().map(|p|(p+1).to_string()).unwrap_or_else(||"none".into());
                    let last=positions.iter().rev().flatten().next().map(|p|(p+1).to_string()).unwrap_or_else(||"none".into());
                    format!("first {first}, last {last}")
                };
                let coordinates=format!("Original source bases (1-based, column order): R {}; Q {}. Circular/reverse rows may be nonascending.",position_label(&current.reference_positions),position_label(&current.query_positions));
                let rows=current.tracks.iter().skip(current_page*TRACK_PAGE).take(TRACK_PAGE).map(|track|{
                    let label=format!("{} feature {} · {} · {}",track.side.label(),track.feature_index+1,track.kind,track.label);
                    let mut rows=String::new();
                    if annotations.get(){rows.push_str(&format!("A {}\n",track.annotation));}
                    if translations.get(){if let Some(protein)=&track.translation {rows.push_str(&format!("T {protein}\n"));}}
                    let notes=track.notes.clone();
                    view!{<article class="pairwise-context-track"><p>{label}</p><pre>{rows}</pre>{notes.into_iter().map(|note|view!{<p class="helper">{note}</p>}).collect_view()}</article>}
                }).collect_view();
                view!{<section><p>{labels}</p><p class="helper">{coordinates}</p>
                    <div class="method-row-actions"><button class="text-tool" disabled=current_page==0 on:click=move|_|tracks_page.set(current_page.saturating_sub(1))>"Previous features"</button><button class="text-tool" disabled={current_page+1>=pages} on:click=move|_|tracks_page.set((current_page+1).min(pages-1))>"Next features"</button></div>
                    <div class="pairwise-context-scroll"><pre>{base_rows}</pre>{rows}</div>
                </section>}.into_any()
            }
        }}
        {move||prepared.get().and_then(Result::ok).map(|context|{
            let count=context.warnings.len();let pages=count.div_ceil(TRACK_PAGE).max(1);let current=warnings_page.get().min(pages-1);
            let warnings=context.warnings.iter().skip(current*TRACK_PAGE).take(TRACK_PAGE).cloned().map(|warning|view!{<li>{warning}</li>}).collect_view();
            view!{<details><summary>{format!("Context limitations and withheld features · {count}")}</summary><ul>{warnings}</ul><p>{format!("Note page {} of {pages}",current+1)}</p><button class="text-tool" disabled=current==0 on:click=move|_|warnings_page.set(current.saturating_sub(1))>"Previous notes"</button><button class="text-tool" disabled={current+1>=pages} on:click=move|_|warnings_page.set((current+1).min(pages-1))>"Next notes"</button></details>}
        })}
    </details>}.into_any()
}
