//! Native fusion placement discovery and explicit preparation materialization.
use crate::{
    assembly_draft_model::{CapturedSource, Draft},
    assembly_fusion_design::{SearchForm, WindowForm},
    worker::CachedWorker,
};
use collagene_core::{
    golden_gate_fusion_design as fusion, golden_gate_junction_edit,
    protocol::{
        request::Request,
        schema::{AnalysisJob, AnalysisResult},
    },
    verification::PlanInput,
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
    pub before: Draft,
    pub epoch: u64,
    pub inputs: Vec<PlanInput>,
    pub base: golden_gate_junction_edit::Request,
    pub sources: Vec<Arc<CapturedSource>>,
    pub regenerate: Vec<bool>,
}
struct SearchEvidence {
    capture: Arc<Capture>,
    form: SearchForm,
    request: fusion::Request,
    report: fusion::Report,
}
pub struct MaterializedEvidence {
    pub capture: Arc<Capture>,
    pub report: fusion::Materialized,
}
#[component]
pub fn AssemblyFusionView(
    epoch: RwSignal<u64>,
    invalidation: RwSignal<u32>,
    junction_labels: Memo<Vec<String>>,
    capture: Callback<(), Result<Capture, String>>,
    guard: Callback<Arc<Capture>, Result<(), String>>,
    on_apply: Callback<Arc<MaterializedEvidence>, Result<(), String>>,
) -> impl IntoView {
    let form = RwSignal::new(SearchForm::new(junction_labels.with_untracked(Vec::len)));
    let evidence = RwSignal::new(None::<Arc<SearchEvidence>>);
    let materialized = RwSignal::new(None::<Arc<MaterializedEvidence>>);
    let dataset_text = RwSignal::new(String::new());
    let fidelity_candidate = RwSignal::new(None::<String>);
    let pending = RwSignal::new(false);
    let page = RwSignal::new(0usize);
    let status = RwSignal::new(String::new());
    let generation = RwSignal::new(0u32);
    let worker = StoredValue::new_local(None::<Rc<CachedWorker>>);
    let alive = Arc::new(AtomicBool::new(true));
    let cleanup = alive.clone();
    on_cleanup(move || {
        cleanup.store(false, Ordering::Relaxed);
        if let Some(worker) = worker.get_value() {
            worker.close();
        }
    });
    let cancel = Callback::new(move |()| {
        generation.update(|value| *value = value.wrapping_add(1));
        if let Some(worker) = worker.get_value() {
            worker.close();
        }
        worker.set_value(None);
        pending.set(false);
        evidence.set(None);
        materialized.set(None);
        fidelity_candidate.set(None);
        page.set(0);
    });
    Effect::new(move |_| {
        let _ = (epoch.get(), invalidation.get(), form.get());
        let count = junction_labels.get().len();
        cancel.run(());
        if form.with_untracked(|value| value.windows.len() != count) {
            form.update(|value| {
                value.windows = vec![
                    WindowForm::default();
                    count.min(crate::assembly_fusion_design::MAX_JUNCTIONS)
                ]
            });
        }
    });
    let search_guard = Callback::new(move |saved: Arc<SearchEvidence>| -> Result<(), String> {
        guard.run(saved.capture.clone())?;
        if saved.form != form.get_untracked()
            || !evidence.with_untracked(|current| {
                current
                    .as_ref()
                    .is_some_and(|current| Arc::ptr_eq(current, &saved))
            })
        {
            return Err("Fusion settings or source context changed; search again.".into());
        }
        Ok(())
    });
    let search_alive = alive.clone();
    let search = move |_| {
        if pending.get_untracked() {
            return;
        }
        let result = (|| -> Result<_, String> {
            let bound = Arc::new(capture.run(())?);
            guard.run(bound.clone())?;
            let settings = form.get_untracked();
            let request = settings.request(bound.base.clone())?;
            let revision = fusion::inspect_revision(&bound.inputs, &request)
                .map_err(|error| error.to_string())?;
            let worker = Rc::new(CachedWorker::new()?);
            Ok((bound, settings, request, revision, worker))
        })();
        let (bound, settings, request, revision, w) = match result {
            Ok(value) => value,
            Err(error) => {
                status.set(error);
                return;
            }
        };
        cancel.run(());
        let token = generation.get_untracked().wrapping_add(1);
        generation.set(token);
        pending.set(true);
        worker.set_value(Some(w.clone()));
        status.set("Searching declared fusion-position windows in Rust; no preparation or target has been changed.".into());
        let alive = search_alive.clone();
        spawn_local(async move {
            let response = w
                .ask(Request {
                    id: token,
                    revision: revision.clone(),
                    job: AnalysisJob::SearchGoldenGateFusions {
                        inputs: bound.inputs.clone(),
                        request: request.clone(),
                    },
                })
                .await;
            w.close();
            if !alive.load(Ordering::Relaxed) || generation.get_untracked() != token {
                return;
            }
            pending.set(false);
            worker.set_value(None);
            let outcome = (|| -> Result<Arc<SearchEvidence>, String> {
                guard.run(bound.clone())?;
                if form.get_untracked() != settings {
                    return Err("Fusion search settings changed.".into());
                }
                let response = response?;
                if response.id != token || response.revision != revision {
                    return Err("Fusion search response differs from its request.".into());
                }
                let report = match response.result {
                    AnalysisResult::GoldenGateFusions(report) => report,
                    AnalysisResult::Error(error) => return Err(error),
                    _ => return Err("Unexpected fusion search response.".into()),
                };
                report
                    .validate_binding(&revision)
                    .map_err(|error| error.to_string())?;
                Ok(Arc::new(SearchEvidence {
                    capture: bound,
                    form: settings,
                    request,
                    report,
                }))
            })();
            match outcome {
                Ok(saved) => {
                    evidence.set(Some(saved));
                    status.set("Review candidate placements, then explicitly materialize one to check its original-source PCR and digest preparation.".into());
                }
                Err(error) => status.set(error),
            }
        });
    };
    let materialize_alive = alive.clone();
    let materialize = Callback::new(move |candidate_id: String| {
        if pending.get_untracked() {
            return;
        }
        let prepared = (|| -> Result<_, String> {
            let saved = evidence
                .get_untracked()
                .ok_or("Search current fusion windows first.")?;
            search_guard.run(saved.clone())?;
            if !saved
                .report
                .candidates
                .iter()
                .any(|candidate| candidate.id == candidate_id)
            {
                return Err("Choose a candidate from this exact search.".into());
            }
            let selection = fusion::Selection {
                context_revision: saved.report.context_revision.clone(),
                candidate_id,
            };
            let revision =
                fusion::materialize_revision(&saved.capture.inputs, &saved.request, &selection)
                    .map_err(|error| error.to_string())?;
            Ok((saved, selection, revision, Rc::new(CachedWorker::new()?)))
        })();
        let (saved, selection, revision, w) = match prepared {
            Ok(value) => value,
            Err(error) => {
                status.set(error);
                return;
            }
        };
        let token = generation.get_untracked().wrapping_add(1);
        generation.set(token);
        pending.set(true);
        materialized.set(None);
        worker.set_value(Some(w.clone()));
        status.set("Materializing the selected placement through native original-source primers and actual Type IIS digestion…".into());
        let alive = materialize_alive.clone();
        spawn_local(async move {
            let response = w
                .ask(Request {
                    id: token,
                    revision: revision.clone(),
                    job: AnalysisJob::MaterializeGoldenGateFusions {
                        inputs: saved.capture.inputs.clone(),
                        request: saved.request.clone(),
                        selection: selection.clone(),
                    },
                })
                .await;
            w.close();
            if !alive.load(Ordering::Relaxed) || generation.get_untracked() != token {
                return;
            }
            pending.set(false);
            worker.set_value(None);
            let outcome = (|| -> Result<Arc<MaterializedEvidence>, String> {
                search_guard.run(saved.clone())?;
                let response = response?;
                if response.id != token || response.revision != revision {
                    return Err("Fusion materialization response differs from its request.".into());
                }
                let report = match response.result {
                    AnalysisResult::GoldenGateFusionMaterialized(report) => report,
                    AnalysisResult::Error(error) => return Err(error),
                    _ => return Err("Unexpected fusion materialization response.".into()),
                };
                report
                    .validate_binding(&revision)
                    .map_err(|error| error.to_string())?;
                if report.parent_context != saved.report.context_revision
                    || report.selection.context_revision != selection.context_revision
                    || report.selection.candidate_id != selection.candidate_id
                    || report.candidate.id != selection.candidate_id
                {
                    return Err(
                        "Materialized fusion differs from the explicitly selected candidate."
                            .into(),
                    );
                }
                Ok(Arc::new(MaterializedEvidence {
                    capture: saved.capture.clone(),
                    report,
                }))
            })();
            match outcome {
                Ok(value) => {
                    materialized.set(Some(value));
                    status.set("Materialized preparation ready for explicit Apply. Original source contributions and independent target remain unchanged.".into());
                }
                Err(error) => status.set(error),
            }
        });
    });
    let apply = move |_| {
        let result = (|| -> Result<(), String> {
            let selected = materialized
                .get_untracked()
                .ok_or("Materialize a fusion candidate first.")?;
            let saved = evidence
                .get_untracked()
                .ok_or("The parent fusion search is unavailable.")?;
            search_guard.run(saved.clone())?;
            guard.run(selected.capture.clone())?;
            if selected.report.parent_context != saved.report.context_revision {
                return Err("The materialized placement is stale.".into());
            }
            on_apply.run(selected)
        })();
        if let Err(error) = result {
            status.set(error);
        }
    };
    view!{<details class="assembly-fusion-design"><summary>"Explore Golden Gate fusion positions"</summary>
        <p class="helper">"Shift fusion positions within signed windows around the original junctions, or fix an exact shift. This changes preparation design while keeping the original intended DNA independent. Retained manual preparations stay literal. Optional empirical scoring uses supplied count data and remains separate from these sequence constraints."</p>
        {move||junction_labels.get().into_iter().enumerate().map(|(index,label)|view!{
            <fieldset><legend>{label}</legend>
                <label>"Minimum shift"<input type="number" aria-label=format!("Fusion minimum shift {}",index+1) prop:value=move||form.with(|value|value.windows.get(index).map(|window|window.minimum.clone()).unwrap_or_default()) on:input=move|event|form.update(|value|if let Some(window)=value.windows.get_mut(index){window.minimum=event_target_value(&event)})/></label>
                <label>"Maximum shift"<input type="number" aria-label=format!("Fusion maximum shift {}",index+1) prop:value=move||form.with(|value|value.windows.get(index).map(|window|window.maximum.clone()).unwrap_or_default()) on:input=move|event|form.update(|value|if let Some(window)=value.windows.get_mut(index){window.maximum=event_target_value(&event)})/></label>
                <label>"Fixed shift (optional)"<input type="number" aria-label=format!("Fusion fixed shift {}",index+1) prop:value=move||form.with(|value|value.windows.get(index).map(|window|window.fixed.clone()).unwrap_or_default()) on:input=move|event|form.update(|value|if let Some(window)=value.windows.get_mut(index){window.fixed=event_target_value(&event)})/></label>
            </fieldset>
        }).collect_view()}
        <p class="helper">"Negative shifts move toward upstream sequence; positive shifts move downstream. Coordinates are relative to the declared junction, not original source numbering."</p>
        <label><input type="checkbox" prop:checked=move||form.with(|value|value.reject_palindromes) on:change=move|event|form.update(|value|value.reject_palindromes=event_target_checked(&event))/>"Exclude palindromic fusion words"</label>
        <label><input type="checkbox" prop:checked=move||form.with(|value|value.distinct_reverse_complements) on:change=move|event|form.update(|value|value.distinct_reverse_complements=event_target_checked(&event))/>"Require distinct reverse-complement classes"</label>
        <label>"Minimum Hamming distance"<input type="number" min="0" prop:value=move||form.with(|value|value.minimum_hamming.clone()) on:input=move|event|form.update(|value|value.minimum_hamming=event_target_value(&event))/></label>
        <label>"Maximum candidate placements"<input type="number" min="1" prop:value=move||form.with(|value|value.max_candidates.clone()) on:input=move|event|form.update(|value|value.max_candidates=event_target_value(&event))/></label>
        <label>"Search-state budget"<input type="number" min="1" prop:value=move||form.with(|value|value.max_states.clone()) on:input=move|event|form.update(|value|value.max_states=event_target_value(&event))/></label>
        <button class="button secondary" disabled=move||pending.get() on:click=search>"Search fusion placements"</button>
        <button class="text-tool" disabled=move||!pending.get() on:click=move|_|{cancel.run(());status.set("Fusion design cancelled; no draft changed.".into());}>"Cancel fusion design"</button>
        {move||evidence.get().map(|saved|{
            let total=saved.report.candidates.len();let pages=total.max(1).div_ceil(20);let current=page.get().min(pages-1);
            let ranking_source=saved.clone();
            let ranking_capture=Callback::new(move|()|->Result<crate::assembly_profile_view::Capture,String>{
                search_guard.run(ranking_source.clone())?;
                Ok(crate::assembly_profile_view::Capture{inputs:ranking_source.capture.inputs.clone(),request:ranking_source.request.clone(),search_revision:ranking_source.report.context_revision.clone()})
            });
            let ranking_source=saved.clone();
            let ranking_guard=Callback::new(move|()|search_guard.run(ranking_source.clone()));
            let candidates=saved.report.candidates.iter().skip(current*20).take(20).map(|candidate|{
                let id=candidate.id.clone();let score_id=id.clone();view!{<article><h4>"Candidate placement"</h4>
                    <p>{format!("Signed shifts {:?} · fusion words {}",candidate.shifts,candidate.fusions.join(" · "))}</p>
                    <details><summary>"Exact placement identity"</summary><code>{candidate.id.clone()}</code></details>
                    <button class="button secondary" disabled=move||pending.get() on:click=move|_|materialize.run(id.clone())>"Materialize this placement"</button>
                    <button class="text-tool" disabled=move||pending.get() on:click=move|_|fidelity_candidate.set(Some(score_id.clone()))>"Review empirical score for this set"</button>
                </article>}
            }).collect_view();
            view!{<section><p>{format!("{} candidate placements · {} states · finite search complete {}",total,saved.report.states_examined,saved.report.search_complete)}</p><crate::assembly_profile_view::FusionProfileView capture=ranking_capture guard=ranking_guard on_materialize=materialize/>
                {saved.report.limits.iter().chain(&saved.report.warnings).map(|warning|view!{<p>{warning.clone()}</p>}).collect_view()}{candidates}
                <p>{format!("Candidate page {} of {pages}",current+1)}</p>
                <button class="text-tool" disabled=current==0 on:click=move|_|page.set(current.saturating_sub(1))>"Previous placements"</button>
                <button class="text-tool" disabled={current+1>=pages} on:click=move|_|page.set(current+1)>"Next placements"</button>
            </section>}
        })}
        {move||fidelity_candidate.get().and_then(|id|evidence.get().and_then(|saved|saved.report.candidates.iter().find(|candidate|candidate.id==id).and_then(|candidate|fusion::end_scope(&saved.report,&candidate.id).ok()).map(|scope|{
            view!{<FidelityView saved=saved.clone() ends=scope.intended_junctions physical_competitors=scope.physical_competitors dataset_text guard=search_guard/>}
        })))}
        {move||materialized.get().map(|saved|view!{<section><h4>"Materialized original-source preparations"</h4>
            <p>{format!("Independent target: {} · {} bases",saved.report.original_intent.target.name,saved.report.original_intent.target.sequence.len())}</p>
            {saved.report.placements.iter().enumerate().map(|(index,placement)|view!{<p>{format!("Source row {}: {}–{}{} · prefix {} · suffix {}",index+1,placement.start+1,placement.end,if placement.reverse{" reverse"}else{" forward"},placement.prefix,placement.suffix)}</p>}).collect_view()}
            {saved.report.designs.iter().enumerate().filter_map(|(index,row)|row.as_ref().map(|row|(index,row))).map(|(index,row)|view!{<article><h5>{format!("Regenerated row {}",index+1)}</h5>
                <pre>{format!("{}: {}{}\n{}: {}{}",row.forward.name,row.forward.tail,row.forward.annealing.as_str(),row.reverse.name,row.reverse.tail,row.reverse.annealing.as_str())}</pre>
                <p>{format!("Pinned digest fragment {} · left {:?} {} · right {:?} {}",row.fragment_index+1,row.left.polarity,row.left.overhang,row.right.polarity,row.right.overhang)}</p>
            </article>}).collect_view()}
            {saved.report.warnings.iter().map(|warning|view!{<p>{warning.clone()}</p>}).collect_view()}
            <button class="button primary" disabled=move||pending.get() on:click=apply>"Apply fusion placement and selected primers"</button>
        </section>})}
        <p role="status">{move||status.get()}</p>
    </details>}.into_any()
}

#[component]
fn FidelityView(
    saved: Arc<SearchEvidence>,
    ends: Vec<String>,
    physical_competitors: Vec<String>,
    dataset_text: RwSignal<String>,
    guard: Callback<Arc<SearchEvidence>, Result<(), String>>,
) -> impl IntoView {
    use collagene_core::ligation_fidelity as fidelity;
    let saved = StoredValue::new(saved);
    let ends = StoredValue::new(ends);
    let physical_competitors = StoredValue::new(physical_competitors);
    let report = RwSignal::new(None::<Arc<fidelity::CompetitionReport>>);
    let pending = RwSignal::new(false);
    let status = RwSignal::new(String::new());
    let generation = RwSignal::new(0u32);
    let page = RwSignal::new(0usize);
    let worker = StoredValue::new_local(None::<Rc<CachedWorker>>);
    let alive = Arc::new(AtomicBool::new(true));
    let cleanup = alive.clone();
    on_cleanup(move || {
        cleanup.store(false, Ordering::Relaxed);
        if let Some(worker) = worker.get_value() {
            worker.close();
        }
    });
    let cancel = Callback::new(move |()| {
        generation.update(|value| *value = value.wrapping_add(1));
        if let Some(worker) = worker.get_value() {
            worker.close();
        }
        worker.set_value(None);
        pending.set(false);
        report.set(None);
        page.set(0);
    });
    Effect::new(move |_| {
        let _ = dataset_text.get();
        cancel.run(());
    });
    let upload_alive = alive.clone();
    let upload = move |event| {
        let input = event_target::<web_sys::HtmlInputElement>(&event);
        let Some(file) = input.files().and_then(|files| files.get(0)) else {
            return;
        };
        if file.size() > fidelity::MAX_BYTES as f64 {
            status.set("Provided matrix exceeds16MiB.".into());
            return;
        }
        if let Err(error) = guard.run(saved.get_value()) {
            status.set(error);
            return;
        }
        cancel.run(());
        let token = generation.get_untracked().wrapping_add(1);
        generation.set(token);
        pending.set(true);
        status.set("Reading the selected count matrix; scoring waits for this file.".into());
        let alive = upload_alive.clone();
        spawn_local(async move {
            let value = wasm_bindgen_futures::JsFuture::from(file.text()).await;
            if !alive.load(Ordering::Relaxed) || generation.get_untracked() != token {
                return;
            }
            pending.set(false);
            let result = (|| -> Result<String, String> {
                guard.run(saved.get_value())?;
                let text = value
                    .map_err(|_| "Could not read the selected matrix file.")?
                    .as_string()
                    .ok_or("Matrix file must contain text JSON.")?;
                if text.len() > fidelity::MAX_BYTES {
                    return Err("Provided matrix exceeds16MiB.".into());
                }
                Ok(text)
            })();
            match result {
                Ok(text) => {
                    dataset_text.set(text);
                    status.set("Provided matrix loaded as editable text. Explicitly score this selected end set.".into());
                }
                Err(error) => status.set(error),
            }
        });
    };
    let score_alive = alive.clone();
    let score = move |_| {
        if pending.get_untracked() {
            return;
        }
        let prepared = (|| -> Result<_, String> {
            guard.run(saved.get_value())?;
            let text = dataset_text.get_untracked();
            let dataset = fidelity::decode(text.as_bytes()).map_err(|error| error.to_string())?;
            let requested = ends.get_value();
            let revision = fidelity::competition_revision(
                &dataset,
                &requested,
                &physical_competitors.get_value(),
            )
            .map_err(|error| error.to_string())?;
            Ok((
                text,
                dataset,
                requested,
                revision,
                Rc::new(CachedWorker::new()?),
            ))
        })();
        let (text, dataset, requested, revision, w) = match prepared {
            Ok(value) => value,
            Err(error) => {
                status.set(error);
                return;
            }
        };
        cancel.run(());
        let token = generation.get_untracked().wrapping_add(1);
        generation.set(token);
        worker.set_value(Some(w.clone()));
        pending.set(true);
        status.set("Scoring the exact selected fusion set using the supplied directed counts and conditions…".into());
        let alive = score_alive.clone();
        spawn_local(async move {
            let response = w
                .ask(Request {
                    id: token,
                    revision: revision.clone(),
                    job: AnalysisJob::ScoreLigationCompetition {
                        dataset: dataset.clone(),
                        ends: requested.clone(),
                        physical_competitors: physical_competitors.get_value(),
                    },
                })
                .await;
            w.close();
            if !alive.load(Ordering::Relaxed) || generation.get_untracked() != token {
                return;
            }
            pending.set(false);
            worker.set_value(None);
            let result = (|| -> Result<Arc<fidelity::CompetitionReport>, String> {
                guard.run(saved.get_value())?;
                if dataset_text.get_untracked() != text {
                    return Err("Provided matrix text changed; score again.".into());
                }
                let response = response?;
                if response.id != token || response.revision != revision {
                    return Err("Empirical score response differs from its request.".into());
                }
                let scored = match response.result {
                    AnalysisResult::LigationCompetition(value) => value,
                    AnalysisResult::Error(error) => return Err(error),
                    _ => return Err("Unexpected empirical count score response.".into()),
                };
                scored
                    .validate_binding(&revision)
                    .map_err(|error| error.to_string())?;
                if scored.dataset_revision != dataset.content_revision
                    || scored.requested_ends != requested
                    || scored.physical_competitors != physical_competitors.get_value()
                {
                    return Err("Score is bound to different counts or fusion ends.".into());
                }
                Ok(Arc::new(scored))
            })();
            match result {
                Ok(value) => {
                    report.set(Some(value));
                    status.set("Condition-specific count score ready. It does not select a winner or replace physical preparation checks.".into());
                }
                Err(error) => status.set(error),
            }
        });
    };
    view!{<details open class="ligation-fidelity"><summary>"Provided empirical count matrix"</summary>
        <p class="helper">"Supply a licensed matrix in Collagene directed-count JSON format, including conditions and provenance. Missing counts stay unknown. No bundled dataset is assumed, and no winner is selected automatically."</p>
        <p>{ends.with_value(|ends|format!("Intended junctions: {}",ends.join(" · ")))}</p>
        <p>{physical_competitors.with_value(|ends|format!("Additional physical competing ends (5′→3′): {}",ends.join(" · ")))}</p>
        <label>"Matrix JSON"<textarea maxlength=fidelity::MAX_BYTES.to_string() prop:value=move||dataset_text.get() on:input=move|event|{
            let input=event_target::<web_sys::HtmlTextAreaElement>(&event);let value=input.value();
            if value.len()>fidelity::MAX_BYTES{input.set_value(&dataset_text.get_untracked());status.set("Provided matrix exceeds16MiB; previous text retained.".into());}else{dataset_text.set(value);}
        }/></label>
        <label>"Load a provided JSON matrix"<input type="file" accept=".json,application/json" on:change=upload/></label>
        <button class="button secondary" disabled=move||pending.get() on:click=score>"Score this selected set"</button>
        <button class="text-tool" disabled=move||!pending.get() on:click=move|_|{cancel.run(());status.set("Empirical scoring cancelled.".into());}>"Cancel empirical scoring"</button>
        {move||report.get().map(|value|{
            let count=value.pairs.len();let pages=count.max(1).div_ceil(20);let current=page.get().min(pages-1);
            view!{<section><h4>{value.overall_fidelity.map(|score|format!("Count-model combined estimate {:.4}%",score*100.0)).unwrap_or_else(||"Combined estimate unknown".into())}</h4>
                <p>{format!("Model {} · natural log of combined estimate {:?}",value.model,value.log_overall_fidelity)}</p>
                <p>{format!("Ligase {} · buffer {} · temperature protocol {} · incubation {} · restriction enzymes {} · substrate context {}",value.conditions.ligase,value.conditions.buffer,value.conditions.temperature_protocol,value.conditions.incubation,value.conditions.restriction_enzymes.join(", "),value.conditions.substrate_context)}</p>
                <details><summary>"Dataset provenance and declared normalization"</summary><p>{value.provenance.source_url.clone()}</p><p>{value.provenance.citation.clone()}</p><p>{format!("License {} · attribution {}",value.provenance.license,value.provenance.attribution)}</p><p>{value.provenance.normalization.clone()}</p><p>{format!("Source data hash {} · dataset revision {}",value.provenance.source_blake3,value.dataset_revision)}</p></details>
                {value.pairs.iter().skip(current*20).take(20).map(|pair|view!{<article><h5>{format!("{} ↔ {}",pair.canonical_end,pair.reverse_complement)}</h5>
                    <p>{match &pair.estimate{fidelity::Estimate::Known{correct_count,total_count,fidelity}=>format!("Correct directed observations {} / total {} · count-model estimate {:.4}%",correct_count,total_count,fidelity*100.0),fidelity::Estimate::Unknown{reason}=>format!("Estimate unknown: {reason:?}")}}</p>
                    {pair.orientations.iter().map(|orientation|view!{<p>{format!("{} → {} · correct count {:?} · total count {:?} · missing partners {}",orientation.end_5to3,orientation.correct_partner_5to3,orientation.correct_count,orientation.total_count,orientation.missing_partners.join(", "))}</p>}).collect_view()}
                </article>}).collect_view()}
                {value.warnings.iter().map(|warning|view!{<p>{warning.clone()}</p>}).collect_view()}
                <p>{format!("Pair score page {} of {pages}",current+1)}</p><button class="text-tool" disabled=current==0 on:click=move|_|page.set(current.saturating_sub(1))>"Previous pair scores"</button><button class="text-tool" disabled={current+1>=pages} on:click=move|_|page.set(current+1)>"Next pair scores"</button>
            </section>}
        })}
        <p role="status">{move||status.get()}</p>
    </details>}.into_any()
}
