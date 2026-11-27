use crate::{
    sequence_similarity_model::{Capture, Prepared},
    source_identity::SourceInfo,
    worker::CachedWorker,
};
use collagene_core::{
    Alphabet,
    protocol::{AnalysisJob, AnalysisResult, Request},
    run_lineage::EventBinding,
    sequence_similarity as native,
};
use leptos::prelude::*;
use std::{
    collections::BTreeSet,
    rc::Rc,
    sync::{
        Arc,
        atomic::{AtomicBool, Ordering},
    },
};
#[derive(Clone)]
struct Reviewed {
    capture: Capture,
    report: native::Report,
}
#[component]
pub fn SimilaritySearchPanel(
    active_source: Memo<SourceInfo>,
    subjects: Memo<Vec<(usize, String)>>,
    source_info: Callback<usize, Option<SourceInfo>>,
    epoch: RwSignal<u64>,
    prepare: Callback<(EventBinding, Vec<usize>, native::Options), Result<Prepared, String>>,
    guard: Callback<Capture, Result<(), String>>,
    on_navigate: Callback<crate::compare::Navigation, Result<(), String>>,
) -> impl IntoView {
    let pinned = RwSignal::new(active_source.get_untracked());
    let selected = RwSignal::new(None::<BTreeSet<String>>);
    let choices = Memo::new(move |_| {
        subjects
            .get()
            .into_iter()
            .filter_map(|(i, _)| source_info.run(i).map(|s| (i, s)))
            .collect::<Vec<_>>()
    });
    let chosen = Memo::new(move |_| {
        let q = pinned.get();
        let choices = choices.get();
        let ids = selected.get().unwrap_or_else(|| {
            choices
                .iter()
                .filter(|(_, s)| {
                    s.alphabet == q.alphabet && s.binding.document_id != q.binding.document_id
                })
                .map(|(_, s)| s.binding.document_id.clone())
                .collect()
        });
        if ids
            .iter()
            .any(|id| !choices.iter().any(|(_, s)| s.binding.document_id == *id))
        {
            return Err("A selected document is unavailable. Review the scope again.".to_owned());
        }
        let indices = choices
            .iter()
            .filter(|(_, s)| ids.contains(&s.binding.document_id))
            .map(|(i, _)| *i)
            .collect::<Vec<_>>();
        if indices.is_empty() || indices.len() > native::MAX_SUBJECTS {
            return Err("Choose 1–128 documents to search.".to_owned());
        }
        Ok(indices)
    });
    let both = RwSignal::new(true);
    let origins = RwSignal::new("8".to_owned());
    let cells = RwSignal::new("8000000".to_owned());
    let work = RwSignal::new("64000000".to_owned());
    let hit_limit = RwSignal::new("100".to_owned());
    let minimum = RwSignal::new("1".to_owned());
    let matches = RwSignal::new("2".to_owned());
    let ambiguous = RwSignal::new("0".to_owned());
    let mismatch = RwSignal::new("-3".to_owned());
    let opening = RwSignal::new("-4".to_owned());
    let extension = RwSignal::new("-1".to_owned());
    let options = Memo::new(move |_| -> Result<native::Options, String> {
        let signed = |v: String| {
            v.parse::<i32>()
                .map_err(|_| "Scores must be integers.".to_owned())
        };
        let unsigned = |v: String| {
            v.parse::<usize>()
                .map_err(|_| "Search limits must be positive integers.".to_owned())
        };
        let o = native::Options {
            scoring: collagene_core::alignment::Scoring {
                match_score: signed(matches.get())?,
                ambiguous_score: signed(ambiguous.get())?,
                mismatch_score: signed(mismatch.get())?,
                gap_open: signed(opening.get())?,
                gap_extend: signed(extension.get())?,
            },
            both_query_strands: both.get(),
            max_origins: unsigned(origins.get())?,
            max_cells_per_alignment: unsigned(cells.get())?,
            max_work: unsigned(work.get())?,
            max_hits: unsigned(hit_limit.get())?,
            minimum_score: signed(minimum.get())?,
        };
        o.validate().map_err(|e| e.to_string())?;
        Ok(o)
    });
    let report = RwSignal::new(None::<Arc<Reviewed>>);
    let page = RwSignal::new(0usize);
    let busy = RwSignal::new(false);
    let status = RwSignal::new(String::new());
    let generation = RwSignal::new(0u32);
    let worker = StoredValue::new_local(None::<Rc<CachedWorker>>);
    let alive = Arc::new(AtomicBool::new(true));
    let cleanup = alive.clone();
    on_cleanup(move || {
        cleanup.store(false, Ordering::Relaxed);
        if let Some(w) = worker.get_value() {
            w.close();
        }
    });
    Effect::new(move |_| {
        let _ = (
            pinned.get(),
            choices.get(),
            selected.get(),
            options.get(),
            epoch.get(),
        );
        generation.update(|n| *n = n.wrapping_add(1).max(1));
        if let Some(w) = worker.get_value() {
            w.close();
        }
        worker.set_value(None);
        busy.set(false);
        report.set(None);
        page.set(0);
    });
    let run = {
        let alive = alive.clone();
        move |_| {
            if busy.get_untracked() {
                return;
            }
            let request = (|| -> Result<_, String> {
                let o = options.get_untracked()?;
                let indices = chosen.get_untracked()?;
                let q = pinned.get_untracked();
                if q.alphabet == Alphabet::Protein {
                    return Err("Choose a DNA or RNA query.".into());
                }
                let prepared = prepare.run((q.binding, indices, o))?;
                guard.run(prepared.capture.clone())?;
                Ok(prepared)
            })();
            let prepared = match request {
                Ok(p) => p,
                Err(e) => {
                    status.set(e);
                    return;
                }
            };
            let session = match CachedWorker::new() {
                Ok(w) => Rc::new(w),
                Err(e) => {
                    status.set(e);
                    return;
                }
            };
            let token = generation.get_untracked().wrapping_add(1).max(1);
            generation.set(token);
            busy.set(true);
            report.set(None);
            status.set("Searching selected documents…".into());
            worker.set_value(Some(session.clone()));
            let alive = alive.clone();
            leptos::task::spawn_local(async move {
                let capture = prepared.capture;
                let reply = session
                    .ask(Request {
                        id: token,
                        revision: capture.context_revision.clone(),
                        job: AnalysisJob::FindSimilarSequences {
                            query: prepared.query,
                            subjects: prepared.subjects,
                            options: capture.options.clone(),
                        },
                    })
                    .await;
                session.close();
                if !alive.load(Ordering::Relaxed) || generation.get_untracked() != token {
                    return;
                }
                worker.set_value(None);
                busy.set(false);
                let checked = (|| -> Result<native::Report, String> {
                    guard.run(capture.clone())?;
                    let reply = reply?;
                    if reply.id != token || reply.revision != capture.context_revision {
                        return Err("Search response context changed.".into());
                    }
                    match reply.result {
                        AnalysisResult::SequenceSimilarity(r) => {
                            capture.validate_report(&r)?;
                            Ok(r)
                        }
                        AnalysisResult::Error(e) => Err(e),
                        _ => Err("Unexpected sequence search response.".into()),
                    }
                })();
                match checked {
                    Ok(r) => {
                        status.set(format!(
                            "{} matches · {}{}",
                            r.hits.len(),
                            if r.complete {
                                "best scores proven for the selected scope"
                            } else {
                                "incomplete search; review limits"
                            },
                            if r.hits_truncated {
                                " · result list limited"
                            } else {
                                ""
                            }
                        ));
                        report.set(Some(Arc::new(Reviewed { capture, report: r })));
                        page.set(0);
                    }
                    Err(e) => status.set(e),
                }
            });
        }
    };
    let cancel = move |_| {
        generation.update(|n| *n = n.wrapping_add(1).max(1));
        if let Some(w) = worker.get_value() {
            w.close();
        }
        worker.set_value(None);
        busy.set(false);
        report.set(None);
        status.set("Search cancelled. Documents are unchanged.".into());
    };
    let download = {
        let alive = alive.clone();
        move |_| {
            let Some(r) = report.get_untracked() else {
                return;
            };
            if let Err(e) = guard.run(r.capture.clone()) {
                status.set(e);
                return;
            }
            let bytes = match serde_json::to_vec_pretty(
                &serde_json::json!({"capture":r.capture,"report":r.report}),
            ) {
                Ok(b) => b,
                Err(e) => {
                    status.set(e.to_string());
                    return;
                }
            };
            let token = generation.get_untracked();
            let alive = alive.clone();
            let capture = r.capture.clone();
            crate::app::download_with_completion_guarded(
                "sequence-similarity.json",
                "application/json",
                &bytes,
                move || {
                    if !alive.load(Ordering::Relaxed) || generation.get_untracked() != token {
                        return Err("Search context changed; export cancelled.".into());
                    }
                    guard.run(capture.clone())
                },
                move |result| {
                    if let Err(e) = result {
                        status.set(e);
                    }
                },
            );
        }
    };
    view! {
        <details class="similarity-search"><summary>"Search similar sequences in this workspace"</summary>
            <p class="helper">"DNA/RNA: one best local match per document. This compares reference sequences; it does not verify a cloning plan or sample."</p>
            <p>{move||format!("Pinned query: {} · {} bases",pinned.get().name,pinned.get().length)}</p>
            <Show when=move||pinned.get().alphabet==Alphabet::Protein><p class="helper">"Choose a DNA or RNA document as the query."</p></Show>
            <button class="text-tool" on:click=move|_|{pinned.set(active_source.get_untracked());selected.set(None);}>"Use active document as query"</button>
            <details><summary>"Choose documents"</summary>
                <button class="text-tool" on:click=move|_|selected.set(None)>"Select compatible documents"</button>
                <button class="text-tool" on:click=move|_|selected.set(Some(BTreeSet::new()))>"Clear selection"</button>
                <div class="similarity-scope">{move||choices.get().into_iter().map(|(_,info)|{
                    let id=info.binding.document_id.clone();let check=id.clone();let toggle=id.clone();let alpha=info.alphabet;
                    view!{<label><input type="checkbox" disabled=move||alpha!=pinned.get().alphabet prop:checked=move||{
                        selected.with(|s|s.as_ref().map_or(check!=pinned.get().binding.document_id&&alpha==pinned.get().alphabet,|s|s.contains(&check)))
                    } on:change=move|event|{
                        let mut ids=selected.get_untracked().unwrap_or_else(||choices.get_untracked().iter().filter(|(_,s)|s.alphabet==pinned.get_untracked().alphabet&&s.binding.document_id!=pinned.get_untracked().binding.document_id).map(|(_,s)|s.binding.document_id.clone()).collect());
                        if event_target_checked(&event){ids.insert(toggle.clone());}else{ids.remove(&toggle);}selected.set(Some(ids));
                    }/>{info.name}</label>}
                }).collect_view()}</div>
            </details>
            <label class="check"><input type="checkbox" prop:checked=move||both.get() on:change=move|e|both.set(event_target_checked(&e))/>"Both query strands"</label>
            <details><summary>"Scoring and search limits"</summary>
                <label>"Match score"<input prop:value=move||matches.get() on:input=move|e|matches.set(event_target_value(&e))/></label>
                <label>"Compatible ambiguity score"<input prop:value=move||ambiguous.get() on:input=move|e|ambiguous.set(event_target_value(&e))/></label>
                <label>"Mismatch score"<input prop:value=move||mismatch.get() on:input=move|e|mismatch.set(event_target_value(&e))/></label>
                <label>"Gap opening"<input prop:value=move||opening.get() on:input=move|e|opening.set(event_target_value(&e))/></label>
                <label>"Gap extension"<input prop:value=move||extension.get() on:input=move|e|extension.set(event_target_value(&e))/></label>
                <label>"Circular origin limit"<input prop:value=move||origins.get() on:input=move|e|origins.set(event_target_value(&e))/></label>
                <label>"Cells per alignment"<input prop:value=move||cells.get() on:input=move|e|cells.set(event_target_value(&e))/></label>
                <label>"Total work limit"<input prop:value=move||work.get() on:input=move|e|work.set(event_target_value(&e))/></label>
                <label>"Result limit"<input prop:value=move||hit_limit.get() on:input=move|e|hit_limit.set(event_target_value(&e))/></label>
                <label>"Minimum score"<input prop:value=move||minimum.get() on:input=move|e|minimum.set(event_target_value(&e))/></label>
                <p class="helper">"Query bases use their stored order. Subject circular origins are bounded; repeated matches and tied locations are not enumerated."</p>
            </details>
            <button class="button secondary" disabled=move||busy.get()||chosen.get().is_err()||options.get().is_err()||pinned.get().alphabet==Alphabet::Protein on:click=run>"Search selected documents"</button>
            <button class="text-tool" disabled=move||!busy.get() on:click=cancel>"Cancel"</button>
            <p class="helper">{move||chosen.get().map(|v|format!("{} documents selected",v.len())).unwrap_or_else(|e|e)}</p>
            <p class="helper">{move||options.get().err().unwrap_or_default()}</p>
            <p role="status" aria-live="polite">{move||status.get()}</p>
            {move||report.get().map(|review|{
                let result=review.report.clone();let count=result.hits.len();let pins=review.capture.subjects.clone();let captured_epoch=review.capture.epoch;
                view!{
                    <ol class="similarity-hits">{result.hits.into_iter().skip(page.get()*20).take(20).map(|hit|{
                        let pin=pins.iter().find(|p|p.binding.document_id==hit.subject_key).expect("validated result subject").binding.clone();
                        let name=choices.get_untracked().iter().find(|(_,s)|s.binding.document_id==hit.subject_key).map(|(_,s)|s.name.clone()).unwrap_or_else(||hit.subject_key.clone());
                        let spans=hit.reference_spans.iter().map(|s|(s.start,s.end)).collect::<Vec<_>>();let capture=review.capture.clone();
                        let label=format!("{} · score {} · {:.1}% query paired · {}{}",name,hit.score,hit.statistics.query_coverage*100.0,if hit.query_reversed{"reverse"}else{"forward"},if hit.reference_spans.len()>1{" · crosses origin"}else{""});
                        view!{<li><button class="text-tool" on:click=move|_|{
                            let result=guard.run(capture.clone()).and_then(|_|on_navigate.run(crate::compare::Navigation{event:pin.clone(),epoch:captured_epoch,spans:spans.clone()}));if let Err(e)=result{status.set(e);}
                        }>{label}</button><small>{format!("{} certain matches; {} ambiguous-compatible columns",hit.statistics.certain_matches,hit.statistics.ambiguous_compatible)}</small></li>}
                    }).collect_view()}</ol>
                    <div><button class="text-tool" disabled=move||page.get()==0 on:click=move|_|page.update(|p|*p=p.saturating_sub(1))>"Previous matches"</button><button class="text-tool" disabled=move||(page.get()+1)*20>=count on:click=move|_|page.update(|p|*p+=1)>"Next matches"</button></div>
                    <details><summary>"Search evidence and limits"</summary><ul>{result.subjects.into_iter().filter(|s|!s.limits.is_empty()).map(|s|view!{<li>{format!("{}: {}",s.key,s.limits.iter().map(limit_label).collect::<Vec<_>>().join(", "))}</li>}).collect_view()}</ul><p class="helper">{result.warnings.join(" ")}</p></details>
                }
            })}
            <button class="text-tool" disabled=move||report.with(Option::is_none) on:click=download>"Download search evidence"</button>
        </details>
    }
}
fn limit_label(limit: &native::Limit) -> String {
    match limit {
        native::Limit::Cancelled => "cancelled".into(),
        native::Limit::Work => "work limit".into(),
        native::Limit::AlignmentCells => "alignment size limit".into(),
        native::Limit::CircularOrigins => "circular origins not exhausted".into(),
        native::Limit::AlignmentUnavailable(e) => format!("alignment unavailable: {e}"),
        native::Limit::IncompatibleAlphabet => "different sequence alphabet".into(),
    }
}
