//! Published-profile ranking stays advisory and uses the existing explicit materializer.
use crate::platform::worker::CachedWorker;
use collagene_core::{
    cloning::{
        golden_gate::{fusion_design as fusion, fusion_rank as rank},
        ligation::profiles,
    },
    planning::plan::PlanInput,
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
pub struct Capture {
    pub inputs: Vec<PlanInput>,
    pub request: fusion::Request,
    pub search_revision: String,
}
struct Ranked {
    choice: usize,
    revision: String,
    report: rank::Report,
}
struct Detail {
    candidate_id: String,
    report: profiles::Report,
}
enum Work {
    Rank,
    Detail(String),
}
fn estimate(value: Option<f64>) -> String {
    value
        .map(|v| format!("{:.4}%", v * 100.0))
        .unwrap_or_else(|| "Unknown".into())
}
#[component]
pub fn FusionProfileView(
    capture: Callback<(), Result<Capture, String>>,
    guard: Callback<(), Result<(), String>>,
    on_materialize: Callback<String>,
) -> impl IntoView {
    // Cached fixed metadata only. XML/matrix preparation happens in the worker.
    let (catalog, error) = match profiles::list() {
        Ok(v) => (v, String::new()),
        Err(e) => (vec![], e.to_string()),
    };
    let catalog = StoredValue::new(catalog);
    let choice = RwSignal::new(None::<usize>);
    let ranked = RwSignal::new(None::<Arc<Ranked>>);
    let detail = RwSignal::new(None::<Arc<Detail>>);
    let pending = RwSignal::new(false);
    let status = RwSignal::new(error);
    let page = RwSignal::new(0usize);
    let pair_page = RwSignal::new(0usize);
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
    let cancel = Callback::new(move |clear: bool| {
        generation.update(|v| *v = v.wrapping_add(1));
        if let Some(w) = worker.get_value() {
            w.close();
        }
        worker.set_value(None);
        pending.set(false);
        detail.set(None);
        pair_page.set(0);
        if clear {
            ranked.set(None);
            page.set(0);
        }
    });
    Effect::new(move |_| {
        let _ = choice.get();
        cancel.run(true);
    });
    let run_alive = alive.clone();
    let run = Callback::new(move |work: Work| {
        if pending.get_untracked() {
            return;
        }
        let result = (|| -> Result<_, String> {
            guard.run(())?;
            let selected = choice
                .get_untracked()
                .ok_or("Choose an explicitly measured experiment.")?;
            let profile = catalog
                .with_value(|all| all.get(selected).cloned())
                .ok_or("The profile is unavailable.")?;
            let (revision, job, search_context, parent, selected_id) = match work {
                Work::Rank => {
                    let captured = capture.run(())?;
                    let request = rank::Request {
                        design: captured.request,
                        profile: profile.profile.clone(),
                        enzyme_formulation: profile.enzyme_formulation.clone(),
                    };
                    let revision = rank::inspect_revision(&captured.inputs, &request)
                        .map_err(|e| e.to_string())?;
                    (
                        revision,
                        AnalysisJob::RankGoldenGateFusions {
                            inputs: captured.inputs,
                            request,
                        },
                        Some(captured.search_revision),
                        None,
                        None,
                    )
                }
                Work::Detail(id) => {
                    let saved = ranked
                        .get_untracked()
                        .ok_or("Rank current placements first.")?;
                    if saved.choice != selected {
                        return Err("The experiment changed; rank again.".into());
                    }
                    rank::select(&saved.report, &saved.revision, &id).map_err(|e| e.to_string())?;
                    let scope = collagene_core::cloning::golden_gate::fusion_design::end_scope(
                        &saved.report.search,
                        &id,
                    )
                    .map_err(|e| e.to_string())?;
                    let ends = scope.intended_junctions;
                    let physical_competitors = scope.physical_competitors;
                    let revision = profiles::score_revision_with_competitors(
                        &profile.profile,
                        &ends,
                        &physical_competitors,
                    )
                    .map_err(|e| e.to_string())?;
                    (
                        revision,
                        AnalysisJob::ScorePublishedLigation {
                            profile: profile.profile.clone(),
                            ends,
                            physical_competitors,
                        },
                        None,
                        Some(saved),
                        Some(id),
                    )
                }
            };
            Ok((
                selected,
                revision,
                job,
                search_context,
                parent,
                selected_id,
                Rc::new(CachedWorker::new()?),
            ))
        })();
        let (selected, revision, job, search_context, parent, selected_id, w) = match result {
            Ok(v) => v,
            Err(e) => {
                status.set(e);
                return;
            }
        };
        cancel.run(search_context.is_some());
        let token = generation.get_untracked().wrapping_add(1);
        generation.set(token);
        pending.set(true);
        worker.set_value(Some(w.clone()));
        status.set("Computing the explicitly selected published-frequency model in Rust…".into());
        let alive = run_alive.clone();
        spawn_local(async move {
            let response = w
                .ask(Request {
                    id: token,
                    revision: revision.clone(),
                    job,
                })
                .await;
            w.close();
            if !alive.load(Ordering::Relaxed) || generation.get_untracked() != token {
                return;
            }
            pending.set(false);
            worker.set_value(None);
            let result = (|| -> Result<(), String> {
                guard.run(())?;
                if choice.get_untracked() != Some(selected) {
                    return Err("The selected experiment changed.".into());
                }
                let response = response?;
                if response.id != token || response.revision != revision {
                    return Err("Published-profile response context changed.".into());
                }
                match response.result {
                    AnalysisResult::RankedGoldenGateFusions(report) => {
                        report
                            .validate_binding(&revision)
                            .map_err(|e| e.to_string())?;
                        if Some(&report.search.context_revision) != search_context.as_ref() {
                            return Err(
                                "Ranked sets differ from the original source search.".into()
                            );
                        }
                        ranked.set(Some(Arc::new(Ranked {
                            choice: selected,
                            revision,
                            report,
                        })));
                        page.set(0);
                        status.set("Review ranked sets, then explicitly materialize one through the original primer/digest checks. No preparation was applied.".into());
                    }
                    AnalysisResult::PublishedLigationScore(report) => {
                        report
                            .validate_binding(&revision)
                            .map_err(|e| e.to_string())?;
                        let saved = parent.ok_or("No originating ranked set.")?;
                        if !ranked.with_untracked(|current| {
                            current
                                .as_ref()
                                .is_some_and(|current| Arc::ptr_eq(current, &saved))
                        }) {
                            return Err("The originating ranking changed.".into());
                        }
                        detail.set(Some(Arc::new(Detail {
                            candidate_id: selected_id.ok_or("No selected end set.")?,
                            report,
                        })));
                        pair_page.set(0);
                        status.set("Full conditional pair-weight evidence is ready. Missing or palindromic estimates remain unknown.".into());
                    }
                    AnalysisResult::Error(error) => return Err(error),
                    _ => return Err("Unexpected published-profile response.".into()),
                }
                Ok(())
            })();
            if let Err(error) = result {
                status.set(error);
            }
        });
    });
    let materialize = Callback::new(move |(saved, id): (Arc<Ranked>, String)| {
        let result = (|| -> Result<(), String> {
            guard.run(())?;
            if choice.get_untracked() != Some(saved.choice)
                || !ranked.with_untracked(|current| {
                    current
                        .as_ref()
                        .is_some_and(|current| Arc::ptr_eq(current, &saved))
                })
            {
                return Err("The ranking is stale.".into());
            }
            let selected =
                rank::select(&saved.report, &saved.revision, &id).map_err(|e| e.to_string())?;
            on_materialize.run(selected.candidate_id);
            Ok(())
        })();
        if let Err(error) = result {
            status.set(error);
        }
    });
    view!{<details class="assembly-profile-ranking"><summary>"Rank fusion sets with published experimental data"</summary>
        <p class="helper">"Choose the measured experiment to assume for this estimate. Ranking does not change your recipe's enzyme identity, conditions, intended DNA or primers. A good score does not establish reaction yield."</p>
        <label>"Measured experiment"<select aria-label="Published ligation experiment" prop:value=move||choice.get().map(|i|i.to_string()).unwrap_or_default() on:change=move|e|choice.set(event_target_value(&e).parse::<usize>().ok())>
            <option value="">"Choose a published profile"</option>
            {catalog.with_value(|all|all.iter().enumerate().map(|(i,p)|view!{<option value=i.to_string()>{format!("{} · {}-base ends",p.name,p.overhang_length)}</option>}).collect_view())}
        </select></label>
        {move||choice.get().and_then(|i|catalog.with_value(|all|all.get(i).cloned())).map(|p|view!{
            <p>{format!("Formulation assumed for this estimate: {}",p.enzyme_formulation)}</p>
            <p>{format!("Ligase {} · buffer {} · temperature protocol {} · incubation {}",p.conditions.ligase,p.conditions.buffer,p.conditions.temperature_protocol,p.conditions.incubation)}</p>
            <details><summary>"Publication and data attribution"</summary><p>{p.provenance.citation}</p><p>{p.provenance.license}</p><p>{p.provenance.attribution}</p><p>{p.provenance.normalization}</p><a href=p.provenance.source_url target="_blank" rel="noopener noreferrer">"Original supplementary data"</a></details>
        })}
        <button class="button secondary" disabled=move||pending.get()||choice.get().is_none() on:click=move|_|run.run(Work::Rank)>"Rank current fusion sets"</button>
        <button class="text-tool" disabled=move||!pending.get() on:click=move|_|{cancel.run(false);status.set("Published-profile computation cancelled.".into());}>"Cancel profile computation"</button>
        {move||ranked.get().map(|saved|{
            let total=saved.report.ranked.len();let pages=total.max(1).div_ceil(20);let current=page.get().min(pages-1);
            let rows=saved.report.ranked.iter().skip(current*20).take(20).map(|row|{
                let candidate=saved.report.search.candidates.iter().find(|c|c.id==row.candidate_id);
                let open=saved.clone();let id=row.candidate_id.clone();let details_id=id.clone();
                view!{<article><h4>{format!("Rank {} · conditional estimate {}",row.rank,estimate(row.score.overall_fidelity))}</h4>
                    <p>{candidate.map(|c|format!("Signed shifts {:?} · fusion words {}",c.shifts,c.fusions.join(" · "))).unwrap_or_else(||"Candidate unavailable".into())}</p>
                    <p>{format!("Score class {:?} · log estimate {:?}",row.score.class,row.score.log_overall_fidelity)}</p>
                    <button class="button secondary" disabled=move||pending.get() on:click=move|_|materialize.run((open.clone(),id.clone()))>"Materialize this set"</button>
                    <button class="text-tool" disabled=move||pending.get() on:click=move|_|run.run(Work::Detail(details_id.clone()))>"Inspect pair weights"</button>
                </article>}
            }).collect_view();
            view!{<section><p>{format!("{} retained sets ranked · placement search complete under its constraints: {}",total,saved.report.search.search_complete)}</p>{rows}
                <p>{format!("Ranking page {} of {pages}",current+1)}</p><button class="text-tool" disabled=current==0 on:click=move|_|page.set(current.saturating_sub(1))>"Previous ranked sets"</button><button class="text-tool" disabled={current+1>=pages} on:click=move|_|page.set(current+1)>"Next ranked sets"</button>
                <details><summary>"Assumptions, limits and exact profile identity"</summary>{saved.report.limits.iter().chain(&saved.report.warnings).map(|v|view!{<p>{v.clone()}</p>}).collect_view()}<code>{saved.report.context_revision.clone()}</code><p>{format!("Data revision {} · source BLAKE3 {} · source SHA256 {}",saved.report.metadata.data_revision,saved.report.metadata.provenance.source_blake3,saved.report.metadata.source_sha256)}</p></details>
            </section>}
        })}
        {move||detail.get().map(|saved|{
            let total=saved.report.pairs.len();let pages=total.max(1).div_ceil(20);let current=pair_page.get().min(pages-1);
            view!{<section><h4>"Selected set: published pair-frequency weights"</h4><p>{format!("Conditional estimate {} · {} requested ends",estimate(saved.report.overall_fidelity),saved.report.requested_ends.len())}</p>
                {saved.report.pairs.iter().skip(current*20).take(20).map(|pair|view!{<article><h5>{format!("{} ↔ {}",pair.canonical_end,pair.reverse_complement)}</h5>
                    <p>{match &pair.estimate{profiles::Estimate::Known{correct_weight,total_weight,fidelity}=>format!("Correct-pair weight {} / total competing weight {} · {:.4}%",correct_weight,total_weight,fidelity*100.0),profiles::Estimate::Unknown{reason}=>format!("Estimate unknown: {reason:?}")}}</p>
                    {pair.orientations.iter().map(|o|view!{<p>{format!("{} → {} · correct weight {:?} · total weight {:?} · missing partners {}",o.end_5to3,o.correct_partner_5to3,o.correct_weight,o.total_weight,o.missing_partners.join(", "))}</p>}).collect_view()}
                </article>}).collect_view()}
                <p>{format!("Pair page {} of {pages}",current+1)}</p><button class="text-tool" disabled=current==0 on:click=move|_|pair_page.set(current.saturating_sub(1))>"Previous pairs"</button><button class="text-tool" disabled={current+1>=pages} on:click=move|_|pair_page.set(current+1)>"Next pairs"</button>
                <details><summary>"Score warnings and selected identity"</summary>{saved.report.warnings.iter().map(|w|view!{<p>{w.clone()}</p>}).collect_view()}<code>{saved.candidate_id.clone()}</code></details>
            </section>}
        })}
        <p role="status">{move||status.get()}</p>
    </details>}.into_any()
}
