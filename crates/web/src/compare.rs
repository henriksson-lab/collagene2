//! Exact-version pairwise review. Source navigation does not retarget a pinned reference.
use crate::{
    digest::SourceInfo,
    pairwise_result_view::{PairwiseResultView, Review},
    worker::CachedWorker,
};
use collagene_core::{
    Molecule,
    alignment::{AlignmentMode, AlignmentOptions, ComparisonOptions, Scoring},
    pairwise_history::{Record, Summary},
    pairwise_review::{self, Options},
    protocol::{AnalysisJob, AnalysisResult, Request},
    run_lineage::EventBinding,
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

#[derive(Clone)]
pub struct Navigation {
    pub event: EventBinding,
    pub epoch: u64,
    pub spans: Vec<(usize, usize)>,
}
#[derive(Clone)]
struct Pinned {
    binding: EventBinding,
    molecule: Arc<Molecule>,
}
#[derive(Clone, PartialEq, Eq)]
struct Context {
    reference: EventBinding,
    query: EventBinding,
    options: Options,
    epoch: u64,
}

#[component]
pub fn ComparePanel(
    molecule: Memo<Molecule>,
    query_options: Memo<Vec<(usize, String)>>,
    get_query: Callback<usize, Option<Molecule>>,
    active_source: Memo<SourceInfo>,
    source_info: Callback<usize, Option<SourceInfo>>,
    resolve: Callback<EventBinding, Result<Molecule, String>>,
    epoch: RwSignal<u64>,
    prepare_similarity: Callback<
        (
            EventBinding,
            Vec<usize>,
            collagene_core::sequence_similarity::Options,
        ),
        Result<crate::sequence_similarity_model::Prepared, String>,
    >,
    guard_similarity: Callback<crate::sequence_similarity_model::Capture, Result<(), String>>,
    saved: Memo<Vec<Summary>>,
    on_save: Callback<Record, Result<(), String>>,
    on_reopen: Callback<String, Result<Record, String>>,
    on_navigate: Callback<Navigation, Result<(), String>>,
) -> impl IntoView {
    let reference = RwSignal::new(Pinned {
        binding: active_source.get_untracked().binding,
        molecule: Arc::new(molecule.get_untracked()),
    });
    let saved_query = RwSignal::new(None::<Pinned>);
    let selected = RwSignal::new(
        query_options
            .get_untracked()
            .first()
            .map(|(index, _)| *index)
            .unwrap_or(0),
    );
    let policy = RwSignal::new("stored".to_string());
    let mode = RwSignal::new("global".to_string());
    let match_score = RwSignal::new("2".to_string());
    let ambiguity = RwSignal::new("0".to_string());
    let mismatch = RwSignal::new("-3".to_string());
    let gap_open = RwSignal::new("-4".to_string());
    let gap_extend = RwSignal::new("-1".to_string());
    let max_cells = RwSignal::new("16000000".to_string());
    let both_strands = RwSignal::new(true);
    let rotations = RwSignal::new(true);
    let origins = RwSignal::new("16".to_string());
    let total_cells = RwSignal::new("64000000".to_string());
    let report = RwSignal::new(None::<Arc<Review>>);
    let status = RwSignal::new(String::new());
    let busy = RwSignal::new(false);
    let generation = RwSignal::new(0u32);
    let saved_page = RwSignal::new(0usize);
    let worker = StoredValue::new_local(None::<Rc<CachedWorker>>);
    let alive = Arc::new(AtomicBool::new(true));
    let cleanup_alive = alive.clone();
    on_cleanup(move || {
        cleanup_alive.store(false, Ordering::Relaxed);
        if let Some(w) = worker.get_value() {
            w.close();
        }
    });
    let options = Memo::new(move |_| {
        let number = |value: String, label: &str| {
            value
                .parse::<i32>()
                .map_err(|_| format!("{label} must be an integer"))
        };
        let alignment = AlignmentOptions {
            mode: match mode.get().as_str() {
                "local" => AlignmentMode::Local,
                "query" => AlignmentMode::QueryInReference,
                _ => AlignmentMode::Global,
            },
            scoring: Scoring {
                match_score: number(match_score.get(), "Match score")?,
                ambiguous_score: number(ambiguity.get(), "Ambiguity score")?,
                mismatch_score: number(mismatch.get(), "Mismatch score")?,
                gap_open: number(gap_open.get(), "Gap opening")?,
                gap_extend: number(gap_extend.get(), "Gap extension")?,
            },
            max_cells: max_cells
                .get()
                .parse()
                .map_err(|_| "Enter a positive matrix cell budget".to_string())?,
        };
        if policy.get() == "whole" {
            if alignment.mode != AlignmentMode::Global {
                return Err("Whole-construct comparison requires Global mode. Select Global explicitly or use stored orientation/origin.".to_string());
            }
            Ok(Options::Whole(ComparisonOptions {
                alignment,
                both_query_strands: both_strands.get(),
                reference_rotations: rotations.get(),
                max_origins: origins
                    .get()
                    .parse()
                    .map_err(|_| "Enter an origin-candidate budget".to_string())?,
                max_total_cells: total_cells
                    .get()
                    .parse()
                    .map_err(|_| "Enter a total cell budget".to_string())?,
            }))
        } else {
            Ok(Options::Stored(alignment))
        }
    });
    // Historical reopening explicitly pins both versions. A newly chosen query
    // returns to current-head selection, while reference navigation stays pinned.
    let context = Memo::new(move |_| -> Result<Context, String> {
        let pinned = reference.get();
        let historical = saved_query.get();
        let choices = query_options.get();
        let query = if let Some(query) = historical {
            query.binding
        } else {
            let current_ref = choices
                .iter()
                .filter_map(|(index, _)| source_info.run(*index))
                .find(|info| info.binding.document_id == pinned.binding.document_id)
                .ok_or("The pinned reference document is unavailable")?;
            if current_ref.binding != pinned.binding {
                return Err("The reference head changed. Open the intended reference and choose Use active document as reference, or reopen an explicit saved version.".into());
            }
            source_info
                .run(selected.get())
                .ok_or("Choose an open query document")?
                .binding
        };
        Ok(Context {
            reference: pinned.binding,
            query,
            options: options.get()?,
            epoch: epoch.get(),
        })
    });
    Effect::new(move |_| {
        let current = context.get();
        generation.update(|g| *g = g.wrapping_add(1));
        if let Some(w) = worker.get_value() {
            w.close();
        }
        worker.set_value(None);
        busy.set(false);
        let keep = report.with_untracked(|shown| {
            shown.as_ref().is_some_and(|shown| {
                current.as_ref().is_ok_and(|key| {
                    key.reference == shown.record.reference
                        && key.query == shown.record.query
                        && key.options == shown.record.report.options
                        && key.epoch == shown.epoch
                })
            })
        });
        if !keep {
            report.set(None);
        }
        status.set(current.err().unwrap_or_default());
    });
    let pin_active = move |_| {
        let info = active_source.get_untracked();
        let source = molecule.get_untracked();
        if info.binding.revision != source.revision() {
            status.set("The active reference changed; select it again.".into());
            return;
        }
        batch(|| {
            reference.set(Pinned {
                binding: info.binding,
                molecule: Arc::new(source),
            });
            saved_query.set(None);
            report.set(None);
        });
    };
    let guard = Callback::new(move |shown: Arc<Review>| -> Result<(), String> {
        let key = context.get_untracked()?;
        if !report.with_untracked(|current| {
            current
                .as_ref()
                .is_some_and(|current| Arc::ptr_eq(current, &shown))
        }) || key.epoch != shown.epoch
            || key.reference != shown.record.reference
            || key.query != shown.record.query
            || key.options != shown.record.report.options
        {
            return Err(
                "The source versions, options or displayed comparison changed; compare again."
                    .into(),
            );
        }
        let reference = resolve.run(key.reference)?;
        let query = resolve.run(key.query)?;
        if reference.revision() != shown.record.report.reference_revision
            || query.revision() != shown.record.report.query_revision
        {
            return Err("An exact comparison input event is unavailable or changed.".into());
        }
        Ok(())
    });
    let cancel = Callback::new(move |()| {
        generation.update(|g| *g = g.wrapping_add(1));
        if let Some(w) = worker.get_value() {
            w.close();
        }
        worker.set_value(None);
        busy.set(false);
        status.set("Comparison cancelled. Saved evidence is unchanged.".into());
    });
    let compute = {
        let alive = alive.clone();
        move |_| {
            if busy.get_untracked() {
                return;
            }
            let prepared = (|| -> Result<_, String> {
                let key = context.get_untracked()?;
                let r = resolve.run(key.reference.clone())?;
                let q = if saved_query.get_untracked().is_some() {
                    resolve.run(key.query.clone())?
                } else {
                    get_query
                        .run(selected.get_untracked())
                        .ok_or("Choose an open query document")?
                };
                if r.revision() != key.reference.revision || q.revision() != key.query.revision {
                    return Err("Comparison sources changed; select them again.".into());
                }
                let revision =
                    pairwise_review::revision(&r, &q, &key.options).map_err(|e| e.to_string())?;
                let w = Rc::new(CachedWorker::new()?);
                Ok((key, Arc::new(r), Arc::new(q), revision, w))
            })();
            let (key, r, q, revision, w) = match prepared {
                Ok(v) => v,
                Err(error) => {
                    status.set(error);
                    return;
                }
            };
            let token = generation.get_untracked().wrapping_add(1);
            generation.set(token);
            worker.set_value(Some(w.clone()));
            busy.set(true);
            report.set(None);
            status.set("Comparing the pinned versions in the Rust worker…".into());
            let alive = alive.clone();
            spawn_local(async move {
                let response = w
                    .ask(Request {
                        id: token,
                        revision: revision.clone(),
                        job: AnalysisJob::CompareSequences {
                            reference: r.as_ref().clone(),
                            query: q.as_ref().clone(),
                            options: key.options,
                        },
                    })
                    .await;
                w.close();
                if !alive.load(Ordering::Relaxed) || generation.get_untracked() != token {
                    return;
                }
                worker.set_value(None);
                busy.set(false);
                if context.get_untracked().as_ref() != Ok(&key) {
                    status.set("Sources or options changed; compare again.".into());
                    return;
                }
                let outcome = (|| -> Result<_, String> {
                    let response = response?;
                    if response.id != token || response.revision != revision {
                        return Err("Comparison response differs from its request.".into());
                    }
                    let result = match response.result {
                        AnalysisResult::ReviewedComparison(result) => result,
                        AnalysisResult::Error(error) => return Err(error),
                        _ => return Err("Unexpected comparison response".into()),
                    };
                    result
                        .validate_binding(&revision)
                        .map_err(|e| e.to_string())?;
                    let current_r = resolve.run(key.reference.clone())?;
                    let current_q = resolve.run(key.query.clone())?;
                    if current_r.revision() != result.reference_revision
                        || current_q.revision() != result.query_revision
                    {
                        return Err("Exact input versions changed during comparison.".into());
                    }
                    Ok(Arc::new(Review {
                        record: Record {
                            schema_version: 1,
                            reference: key.reference,
                            query: key.query,
                            report: result,
                        },
                        reference_name: r.name.clone(),
                        query_name: q.name.clone(),
                        reference: r.clone(),
                        query: q.clone(),
                        fresh: true,
                        epoch: key.epoch,
                    }))
                })();
                match outcome {
                    Ok(shown) => {
                        report.set(Some(shown));
                        status.set("Fresh comparison of the pinned versions. Review coverage, ambiguity and unaligned flanks.".into());
                    }
                    Err(error) => status.set(error),
                }
            });
        }
    };
    let reopen = Callback::new(move |hash: String| {
        let result = (|| -> Result<_, String> {
            let record = on_reopen.run(hash)?;
            let r = resolve.run(record.reference.clone())?;
            let q = resolve.run(record.query.clone())?;
            Ok((record, r, q))
        })();
        let (record, r, q) = match result {
            Ok(v) => v,
            Err(error) => {
                status.set(error);
                return;
            }
        };
        let alignment = match record.report.options {
            Options::Stored(a) => a,
            Options::Whole(c) => c.alignment,
        };
        batch(|| {
            reference.set(Pinned {
                binding: record.reference.clone(),
                molecule: Arc::new(r.clone()),
            });
            saved_query.set(Some(Pinned {
                binding: record.query.clone(),
                molecule: Arc::new(q.clone()),
            }));
            match record.report.options {
                Options::Stored(_) => policy.set("stored".into()),
                Options::Whole(c) => {
                    policy.set("whole".into());
                    both_strands.set(c.both_query_strands);
                    rotations.set(c.reference_rotations);
                    origins.set(c.max_origins.to_string());
                    total_cells.set(c.max_total_cells.to_string());
                }
            }
            mode.set(
                match alignment.mode {
                    AlignmentMode::Global => "global",
                    AlignmentMode::Local => "local",
                    AlignmentMode::QueryInReference => "query",
                }
                .into(),
            );
            match_score.set(alignment.scoring.match_score.to_string());
            ambiguity.set(alignment.scoring.ambiguous_score.to_string());
            mismatch.set(alignment.scoring.mismatch_score.to_string());
            gap_open.set(alignment.scoring.gap_open.to_string());
            gap_extend.set(alignment.scoring.gap_extend.to_string());
            max_cells.set(alignment.max_cells.to_string());
            report.set(Some(Arc::new(Review {
                reference_name: r.name.clone(),
                query_name: q.name.clone(),
                reference: Arc::new(r),
                query: Arc::new(q),
                record,
                fresh: false,
                epoch: epoch.get_untracked(),
            })));
            status.set("Historical saved comparison opened. Compute explicitly to replay these pinned versions with the current engine.".into());
        });
    });
    view!{
        <style>".compare-panel .compare-scores{display:grid;grid-template-columns:1fr 1fr;gap:8px}.compare-panel pre{font:11px/1.6 monospace;white-space:pre;overflow:auto;max-height:420px}.compare-panel .compare-result{font-size:13px;line-height:1.6;overflow-wrap:anywhere}.compare-panel .compare-difference{padding:8px 0;border-bottom:1px solid #edf1ed}.compare-panel .compare-difference code{display:block;white-space:pre-wrap;overflow-wrap:anywhere}.compare-panel .compare-warning{color:#796a3e}.compare-panel details{margin:10px 0}"</style>
        <div class="inspector-content compare-panel">
            <crate::sequence_similarity_view::SimilaritySearchPanel active_source=active_source subjects=query_options source_info=source_info epoch=epoch prepare=prepare_similarity guard=guard_similarity on_navigate=on_navigate/>
            <h2>"Pairwise comparison"</h2>
            <p class="helper">{move||reference.with(|r|format!("Pinned reference: {}. Source navigation keeps this reference fixed.",r.molecule.name))}</p>
            <button class="button secondary full" on:click=pin_active>"Use active document as reference"</button>
            <label>"Query document"<select prop:value=move||selected.get().to_string() on:change=move|event|{selected.set(event_target_value(&event).parse().unwrap_or(0));saved_query.set(None);}>
                {move||query_options.get().into_iter().map(|(index,name)|view!{<option value=index.to_string()>{name}</option>}).collect_view()}
            </select></label>
            {move||saved_query.get().map(|query|view!{<p class="helper">{format!("Pinned saved query: {}. Choosing a query document switches to its current version.",query.molecule.name)}</p>})}
            <label>"Comparison policy"<select prop:value=move||policy.get() on:change=move|event|policy.set(event_target_value(&event))><option value="stored">"Stored origin and orientation"</option><option value="whole">"Whole construct · explicit origin/strand search"</option></select></label>
            <label>"Alignment mode"<select prop:value=move||mode.get() on:change=move|event|mode.set(event_target_value(&event))><option value="global">"Global · end-to-end"</option><option value="local">"Local · matching regions"</option><option value="query">"Query end-to-end in reference"</option></select></label>
            <Show when=move||policy.get()=="whole"><p class="helper">"Whole-construct search requires Global mode. Circular-origin candidates can be seed-limited; best score does not establish a unique alignment or physical topology."</p><label><input type="checkbox" prop:checked=move||both_strands.get() on:change=move|event|both_strands.set(event_target_checked(&event))/>"Compare both query strands"</label><label><input type="checkbox" prop:checked=move||rotations.get() on:change=move|event|rotations.set(event_target_checked(&event))/>"Search circular reference origins"</label><label>"Maximum origin candidates"<input type="number" min="1" max="64" prop:value=move||origins.get() on:input=move|event|origins.set(event_target_value(&event))/></label><label>"Total matrix cells"<input type="number" min="1" prop:value=move||total_cells.get() on:input=move|event|total_cells.set(event_target_value(&event))/></label></Show>
            <details><summary>"Scoring and compute budget"</summary><p class="helper">"Gap length k costs opening + (k−1) × extension. Compatible ambiguity is separate from a certain match."</p><div class="compare-scores">{[("Certain match",match_score),("Compatible ambiguity",ambiguity),("Mismatch",mismatch),("Gap opening",gap_open),("Gap extension",gap_extend)].into_iter().map(|(label,field)|view!{<label>{label}<input type="number" prop:value=move||field.get() on:input=move|event|field.set(event_target_value(&event))/></label>}).collect_view()}</div><label>"Matrix cells per alignment"<input type="number" min="1" prop:value=move||max_cells.get() on:input=move|event|max_cells.set(event_target_value(&event))/></label></details>
            <button class="button primary full" disabled=move||busy.get() on:click=compute>"Compare pinned versions"</button><button class="text-tool" disabled=move||!busy.get() on:click=move|_|cancel.run(())>"Cancel comparison"</button>
            <p class="helper" role="status">{move||status.get()}</p>
            {move||report.get().map(|review|view!{<PairwiseResultView review guard on_save on_navigate/>})}
            <details><summary>"Saved comparisons"</summary><p class="helper">"Saved results reopen as historical evidence. They never restore a current verification result."</p>{move||{
                let all=saved.get();let pages=all.len().div_ceil(10).max(1);let current=saved_page.get().min(pages-1);
                let rows=all.into_iter().skip(current*10).take(10).map(|summary|{let hash=summary.attachment_hash;view!{<article><p>{format!("{} → {}",summary.reference_name,summary.query_name)}</p><button class="text-tool" on:click=move|_|reopen.run(hash.clone())>"Open historical comparison"</button></article>}}).collect_view();
                view!{<div>{rows}</div><div class="method-row-actions"><button class="text-tool" disabled=current==0 on:click=move|_|saved_page.set(current.saturating_sub(1))>"Previous comparisons"</button><span>{format!("Page {} of {pages}",current+1)}</span><button class="text-tool" disabled={current+1>=pages} on:click=move|_|saved_page.set((current+1).min(pages-1))>"Next comparisons"</button></div>}
            }}</details>
        </div>
    }.into_any()
}
