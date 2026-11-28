//! Search retained alignment rows in the worker and navigate original coordinates.
use crate::{compare::Navigation, pairwise_result_view::Review, worker::CachedWorker};
use collagene_core::{
    pairwise_search::{self, Matches, Mode, Options, Row, Rows},
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
pub fn PairwiseSearchView(
    review: Arc<Review>,
    guard: Callback<Arc<Review>, Result<(), String>>,
    page: RwSignal<usize>,
    on_navigate: Callback<Navigation, Result<(), String>>,
) -> impl IntoView {
    let bound = StoredValue::new(review);
    let pattern = RwSignal::new(String::new());
    let rows = RwSignal::new("both".to_string());
    let mode = RwSignal::new("literal".to_string());
    let work = RwSignal::new(Options::default().max_work.to_string());
    let hit_limit = RwSignal::new(pairwise_search::MAX_HITS.to_string());
    let busy = RwSignal::new(false);
    let status = RwSignal::new(String::new());
    let results = RwSignal::new(None::<Arc<Matches>>);
    let results_page = RwSignal::new(0usize);
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
        let _ = (
            pattern.get(),
            rows.get(),
            mode.get(),
            work.get(),
            hit_limit.get(),
        );
        generation.update(|value| *value = value.wrapping_add(1));
        if let Some(w) = worker.get_value() {
            w.close();
        }
        worker.set_value(None);
        busy.set(false);
        results.set(None);
        results_page.set(0);
        status.set(String::new());
    });
    let run = {
        let alive = alive.clone();
        move |_| {
            if busy.get_untracked() {
                return;
            }
            let prepared = (|| -> Result<_, String> {
                let review = bound.get_value();
                guard.run(review.clone())?;
                let options = Options {
                    rows: match rows.get_untracked().as_str() {
                        "reference" => Rows::Reference,
                        "query" => Rows::Query,
                        _ => Rows::Both,
                    },
                    mode: if mode.get_untracked() == "compatible" {
                        Mode::Compatible
                    } else {
                        Mode::Literal
                    },
                    max_work: work
                        .get_untracked()
                        .parse()
                        .map_err(|_| "Enter a positive search-work bound.")?,
                    max_hits: hit_limit
                        .get_untracked()
                        .parse()
                        .map_err(|_| "Enter a positive returned-hit bound.")?,
                };
                let text = pattern.get_untracked();
                let revision = pairwise_search::revision(&review.record.report, &text, &options)
                    .map_err(|e| e.to_string())?;
                let w = Rc::new(CachedWorker::new()?);
                Ok((review, text, options, revision, w))
            })();
            let (review, text, options, revision, w) = match prepared {
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
            results.set(None);
            results_page.set(0);
            status.set("Searching the displayed rows in the Rust worker…".into());
            let alive = alive.clone();
            spawn_local(async move {
                let response = w
                    .ask(Request {
                        id: token,
                        revision: revision.clone(),
                        job: AnalysisJob::SearchPairwiseComparison {
                            report: review.record.report.clone(),
                            pattern: text.clone(),
                            options,
                        },
                    })
                    .await;
                w.close();
                if !alive.load(Ordering::Relaxed) || generation.get_untracked() != token {
                    return;
                }
                worker.set_value(None);
                busy.set(false);
                let outcome = (|| -> Result<Matches, String> {
                    guard.run(review.clone())?;
                    let response = response?;
                    if response.id != token || response.revision != revision {
                        return Err(
                            "Search response differs from the retained comparison request.".into(),
                        );
                    }
                    let found = match response.result {
                        AnalysisResult::SearchedPairwiseComparison(found) => found,
                        AnalysisResult::Error(error) => return Err(error),
                        _ => return Err("Unexpected comparison-search response.".into()),
                    };
                    found.validate().map_err(|e| e.to_string())?;
                    if found.pattern != text.to_ascii_uppercase() || found.options != options {
                        return Err(
                            "Search results differ from the selected pattern or policy.".into()
                        );
                    }
                    let report = &review.record.report;
                    for hit in &found.hits {
                        let spans = match hit.row {
                            Row::Reference => report
                                .map_reference(hit.oriented.0, hit.oriented.1)
                                .map_err(|e| e.to_string())?,
                            Row::Query => vec![
                                report
                                    .map_query(hit.oriented.0, hit.oriented.1)
                                    .map_err(|e| e.to_string())?,
                            ],
                        };
                        if hit.columns.1 > report.comparison.alignment.reference_aligned.len()
                            || spans != hit.source_spans
                        {
                            return Err(
                                "Search-hit coordinates differ from the retained comparison."
                                    .into(),
                            );
                        }
                    }
                    guard.run(review)?;
                    Ok(found)
                })();
                match outcome {
                    Ok(found) => {
                        status.set(format!(
                            "{} returned row matches · {} · {} work units",
                            found.hits.len(),
                            if found.complete {
                                "search complete over the selected displayed rows"
                            } else {
                                "search incomplete; additional matches remain unknown"
                            },
                            found.work
                        ));
                        results.set(Some(Arc::new(found)));
                    }
                    Err(error) => status.set(error),
                }
            });
        }
    };
    let select = Callback::new(move |(index, navigate): (usize, bool)| {
        let outcome = (|| -> Result<(), String> {
            let review = bound.get_value();
            guard.run(review.clone())?;
            let result = results
                .get_untracked()
                .ok_or("Run the alignment search first.")?;
            let hit = result
                .hits
                .get(index)
                .ok_or("Search match is unavailable.")?;
            page.set(hit.columns.0 / 600);
            if navigate {
                on_navigate.run(Navigation {
                    event: if hit.row == Row::Reference {
                        review.record.reference.clone()
                    } else {
                        review.record.query.clone()
                    },
                    epoch: review.epoch,
                    spans: hit.source_spans.clone(),
                })?;
            }
            Ok(())
        })();
        if let Err(error) = outcome {
            status.set(error);
        }
    });
    view!{<details class="pairwise-search"><summary>"Find sequence in aligned rows"</summary>
        <p class="helper">"Search the displayed row orientation, skipping alignment gaps. Hits may span intervening gap columns. Unaligned source flanks and a new circular wrap beyond the displayed row are not searched. Reverse-complement or origin changes require an explicit new comparison."</p>
        <label>"Sequence (1–4096 IUPAC calls)"<input maxlength="4096" spellcheck="false" prop:value=move||pattern.get() on:input=move|event|pattern.set(event_target_value(&event))/></label>
        <label>"Rows"<select prop:value=move||rows.get() on:change=move|event|rows.set(event_target_value(&event))><option value="both">"Reference and query"</option><option value="reference">"Reference"</option><option value="query">"Query"</option></select></label>
        <label>"Call comparison"<select prop:value=move||mode.get() on:change=move|event|mode.set(event_target_value(&event))><option value="literal">"Literal IUPAC symbols"</option><option value="compatible">"Compatible IUPAC call sets"</option></select></label>
        <details><summary>"Search bounds"</summary><label>"Work units"<input type="number" min="1" max=pairwise_search::MAX_WORK prop:value=move||work.get() on:input=move|event|work.set(event_target_value(&event))/></label><label>"Returned hits"<input type="number" min="1" max=pairwise_search::MAX_HITS prop:value=move||hit_limit.get() on:input=move|event|hit_limit.set(event_target_value(&event))/></label><p class="helper">"Limits are explicit. An interrupted search cannot establish absence or uniqueness; reduce the rows/pattern ambiguity or revise the bounded settings and rerun."</p></details>
        <button class="button secondary" disabled=move||busy.get() on:click=run>"Search aligned rows"</button>
        <button class="text-tool" disabled=move||!busy.get() on:click=move|_|{generation.update(|v|*v=v.wrapping_add(1));if let Some(w)=worker.get_value(){w.close();}worker.set_value(None);busy.set(false);status.set("Search cancelled.".into());}>"Cancel search"</button>
        {move||results.get().map(|found|{
            let pages=found.hits.len().div_ceil(25).max(1);let current=results_page.get().min(pages-1);
            let matches=found.hits.iter().enumerate().skip(current*25).take(25).map(|(index,hit)|{
                view!{<article><p>{format!("{} · {} · columns {}–{} · original source {}",index+1,if hit.row==Row::Reference{"Reference"}else{"Query"},hit.columns.0+1,hit.columns.1,hit.source_spans.iter().map(|(a,b)|format!("{}–{}",a+1,b)).collect::<Vec<_>>().join(", "))}</p><button class="text-tool" on:click=move|_|select.run((index,false))>"Show alignment columns"</button><button class="text-tool" on:click=move|_|select.run((index,true))>"Show original source bases"</button></article>}
            }).collect_view();
            view!{<div>{matches}<p>{format!("Search-hit page {} of {pages}",current+1)}</p><button class="text-tool" disabled=current==0 on:click=move|_|results_page.set(current.saturating_sub(1))>"Previous hits"</button><button class="text-tool" disabled={current+1>=pages} on:click=move|_|results_page.set(current+1)>"Next hits"</button></div>}
        })}
        <p role="status" class="helper">{move||status.get()}</p>
    </details>}.into_any()
}
