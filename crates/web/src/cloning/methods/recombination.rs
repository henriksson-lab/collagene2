//! Exact-event-bound reference exchange. Positive chemistry remains outside this model.
use crate::platform::worker;
use crate::sources::identity::SourceInfo;
use collagene_core::{
    cloning::recombination::{self, ProfileSnapshot, RecombinationReport, Selection, reference},
    documents::transfer::Candidate,
    material::lineage::{EventBinding, InputEventBinding},
    planning::plan::{DesignPlan, Operation, PlanInput, Source},
    protocol::{
        request::Request as WireRequest,
        schema::{AnalysisJob, AnalysisResult},
    },
};
use collagene_model::Molecule;
use leptos::prelude::*;
use std::sync::{
    Arc,
    atomic::{AtomicBool, Ordering},
};
use wasm_bindgen_futures::spawn_local;

/// The application resolves exact current events and applies borrowed guards before
/// copying any rich source/target document into a worker request.
pub struct CaptureRequest {
    pub sources: [EventBinding; 2],
    pub target: EventBinding,
    pub request: recombination::Request,
}
pub struct Captured {
    pub inputs: Vec<reference::Input>,
    pub target: Molecule,
}
#[derive(Clone, PartialEq, Eq)]
struct Context {
    sources: [SourceInfo; 2],
    target: SourceInfo,
    request: recombination::Request,
    epoch: u64,
}
struct Bound {
    context: Context,
    captured: Arc<Captured>,
    report: RecombinationReport,
}
fn whole<T: std::str::FromStr>(text: &str, name: &str) -> Result<T, String> {
    text.trim()
        .parse()
        .map_err(|_| format!("Enter a whole number for {name}."))
}
fn input_sources() -> [Source; 2] {
    [
        Source::Input("source_1".into()),
        Source::Input("source_2".into()),
    ]
}

#[component]
pub fn RecombinationPanel(
    target: Memo<SourceInfo>,
    sources: Memo<Vec<SourceInfo>>,
    capture: Callback<CaptureRequest, Result<Captured, String>>,
    epoch: RwSignal<u64>,
    locked: Memo<bool>,
    on_candidate: Callback<(Candidate, Vec<InputEventBinding>)>,
) -> impl IntoView {
    let first = RwSignal::new(String::new());
    let second = RwSignal::new(String::new());
    let profile_enabled = RwSignal::new(false);
    let reaction = RwSignal::new(reference::Reaction::BP);
    let comparisons = RwSignal::new("20000000".to_string());
    let candidates = RwSignal::new("64".to_string());
    let product_bases = RwSignal::new("4000000".to_string());
    let profile =
        StoredValue::new(ProfileSnapshot::manufacturer_context().map_err(|e| e.to_string()));
    let context = Memo::new(move |_| -> Result<Context, String> {
        if locked.get() {
            return Err("Wait for the workspace operation to finish.".into());
        }
        if !profile_enabled.get() {
            return Err(
                "Choose the explicit manufacturer-context profile before searching.".into(),
            );
        }
        let available = sources.get();
        let find = |id: &str| {
            available
                .iter()
                .find(|s| s.binding.document_id == id)
                .cloned()
                .ok_or_else(|| "Choose both source documents.".to_string())
        };
        let selected = [find(&first.get())?, find(&second.get())?];
        let target = target.get();
        if selected
            .iter()
            .any(|s| s.binding.document_id == target.binding.document_id)
        {
            return Err("Choose an independent target document as the active document.".into());
        }
        let mut request = recombination::Request::new(profile.get_value()?, reaction.get());
        request.limits.max_comparisons = whole(&comparisons.get(), "comparison budget")?;
        request.limits.max_candidates = whole(&candidates.get(), "candidate limit")?;
        request.limits.max_product_bases = whole(&product_bases.get(), "copied-base budget")?;
        request.validate().map_err(|e| e.to_string())?;
        Ok(Context {
            sources: selected,
            target,
            request,
            epoch: epoch.get(),
        })
    });
    let report = RwSignal::new(None::<Arc<Bound>>);
    let assessment = RwSignal::new(None::<Arc<recombination::TargetAssessment>>);
    let selection = RwSignal::new(None::<Selection>);
    let page = RwSignal::new(0usize);
    let status = RwSignal::new(String::new());
    let busy = RwSignal::new(false);
    let generation = RwSignal::new(0u32);
    let active_job = RwSignal::new(None::<Arc<AtomicBool>>);
    let alive = Arc::new(AtomicBool::new(true));
    {
        let alive = alive.clone();
        on_cleanup(move || {
            alive.store(false, Ordering::Relaxed);
            if let Some(job) = active_job.get_untracked() {
                job.store(false, Ordering::Relaxed);
            }
        });
    }
    Effect::new(move |_| {
        let _ = context.get();
        if let Some(job) = active_job.get_untracked() {
            job.store(false, Ordering::Relaxed);
        }
        generation.update(|g| *g = g.wrapping_add(1));
        report.set(None);
        assessment.set(None);
        selection.set(None);
        busy.set(false);
        status.set(String::new());
        page.set(0);
    });
    let cancel = move |_| {
        if let Some(job) = active_job.get_untracked() {
            job.store(false, Ordering::Relaxed);
        }
        generation.update(|g| *g = g.wrapping_add(1));
        busy.set(false);
        status.set("Computation cancelled.".into());
    };
    let run_alive = alive.clone();
    let run = move |_| {
        let ctx = match context.get_untracked() {
            Ok(c) => c,
            Err(e) => {
                status.set(e);
                return;
            }
        };
        let captured = match capture.run(CaptureRequest {
            sources: ctx.sources.clone().map(|s| s.binding),
            target: ctx.target.binding.clone(),
            request: ctx.request.clone(),
        }) {
            Ok(value) => Arc::new(value),
            Err(e) => {
                status.set(e);
                return;
            }
        };
        let revision = match recombination::batch_revision(&captured.inputs, &ctx.request) {
            Ok(v) => v,
            Err(e) => {
                status.set(e.to_string());
                return;
            }
        };
        if let Some(job) = active_job.get_untracked() {
            job.store(false, Ordering::Relaxed);
        }
        let job_alive = Arc::new(AtomicBool::new(true));
        active_job.set(Some(job_alive.clone()));
        let token = generation.get_untracked().wrapping_add(1);
        generation.set(token);
        busy.set(true);
        report.set(None);
        assessment.set(None);
        selection.set(None);
        page.set(0);
        status.set("Searching full site contexts and reconstructing reciprocal components…".into());
        let alive = run_alive.clone();
        spawn_local(async move {
            let result = worker::analyze_cancellable(
                WireRequest {
                    id: token,
                    revision: revision.clone(),
                    job: AnalysisJob::InspectRecombination {
                        inputs: captured.inputs.clone(),
                        request: ctx.request.clone(),
                    },
                },
                job_alive,
            )
            .await;
            if !alive.load(Ordering::Relaxed) || generation.get_untracked() != token {
                return;
            }
            if context.get_untracked().as_ref().ok() != Some(&ctx) {
                return;
            }
            busy.set(false);
            match result {
                Ok(response) if response.id == token && response.revision == revision => {
                    match response.result {
                        AnalysisResult::RecombinationInspection(value) => {
                            if let Err(e) =
                                value.validate_bindings_for(&captured.inputs, &ctx.request)
                            {
                                status.set(e.to_string());
                                return;
                            }
                            status.set(format!("{} candidates · {} reciprocal components · search complete within this profile: {}. Select a component explicitly.", value.reference.candidates.len(), value.documents.len(), value.reference.complete));
                            report.set(Some(Arc::new(Bound {
                                context: ctx,
                                captured,
                                report: value,
                            })));
                        }
                        AnalysisResult::Error(e) => status.set(e),
                        _ => status.set("Unexpected recombination response.".into()),
                    }
                }
                Ok(_) => status.set("Recombination response binding changed.".into()),
                Err(e) => status.set(e),
            }
        });
    };
    let check_alive = alive.clone();
    let check_target = Callback::new(move |(): ()| {
        let Some(bound) = report.get_untracked() else {
            return;
        };
        if context.get_untracked().as_ref().ok() != Some(&bound.context) {
            return;
        }
        let policy = reference::TargetPolicy {
            allow_reverse_complement: false,
        };
        let revision = match recombination::check_revision(
            &bound.captured.inputs,
            &bound.context.request,
            &bound.captured.target,
            &policy,
        ) {
            Ok(v) => v,
            Err(e) => {
                status.set(e.to_string());
                return;
            }
        };
        if let Some(job) = active_job.get_untracked() {
            job.store(false, Ordering::Relaxed);
        }
        let job_alive = Arc::new(AtomicBool::new(true));
        active_job.set(Some(job_alive.clone()));
        let token = generation.get_untracked().wrapping_add(1);
        generation.set(token);
        busy.set(true);
        assessment.set(None);
        status.set("Reconstructing sources afresh and comparing the independent target…".into());
        let alive = check_alive.clone();
        spawn_local(async move {
            let result = worker::analyze_cancellable(
                WireRequest {
                    id: token,
                    revision: revision.clone(),
                    job: AnalysisJob::CheckRecombinationTarget {
                        inputs: bound.captured.inputs.clone(),
                        request: bound.context.request.clone(),
                        target: bound.captured.target.clone(),
                        policy,
                    },
                },
                job_alive,
            )
            .await;
            if !alive.load(Ordering::Relaxed) || generation.get_untracked() != token {
                return;
            }
            if context.get_untracked().as_ref().ok() != Some(&bound.context) {
                return;
            }
            busy.set(false);
            match result {
                Ok(response) if response.id == token && response.revision == revision => {
                    match response.result {
                        AnalysisResult::RecombinationTargetCheck(value) => {
                            if let Err(e) = value.validate_bindings_for(
                                &bound.captured.inputs,
                                &bound.context.request,
                                &bound.captured.target,
                                &policy,
                            ) {
                                status.set(e.to_string());
                                return;
                            }
                            status.set(format!("{} target matches within the declared model. Target comparison does not select a component or establish chemistry.", value.reference.matches.len()));
                            assessment.set(Some(Arc::new(value)));
                        }
                        AnalysisResult::Error(e) => status.set(e),
                        _ => status.set("Unexpected target-check response.".into()),
                    }
                }
                Ok(_) => status.set("Target-check response binding changed.".into()),
                Err(e) => status.set(e),
            }
        });
    });
    let plan = move |_| {
        let (Some(bound), Some(selected)) = (report.get_untracked(), selection.get_untracked())
        else {
            return;
        };
        if busy.get_untracked() || context.get_untracked().as_ref().ok() != Some(&bound.context) {
            return;
        }
        if !bound
            .report
            .documents
            .iter()
            .any(|d| d.selection == selected)
        {
            return;
        }
        let inputs = bound
            .captured
            .inputs
            .iter()
            .enumerate()
            .map(|(i, input)| PlanInput {
                id: format!("source_{}", i + 1),
                revision: input.molecule.revision(),
                molecule: input.molecule.clone(),
            })
            .collect();
        let preferred = bound
            .context
            .sources
            .iter()
            .enumerate()
            .map(|(i, source)| InputEventBinding {
                input_id: format!("source_{}", i + 1),
                event: source.binding.clone(),
            })
            .collect();
        on_candidate.run((
            Candidate::new(DesignPlan {
                name: format!(
                    "{:?} reference exchange to {}",
                    bound.context.request.reaction, bound.context.target.name
                ),
                target_revision: bound.context.target.binding.revision.clone(),
                inputs,
                steps: vec![Operation::RecombineReference {
                    sources: input_sources(),
                    request: bound.context.request.clone(),
                    selection: selected,
                }],
            }),
            preferred,
        ));
    };
    let export = move |_| {
        let Some(bound) = report.get_untracked() else {
            return;
        };
        if context.get_untracked().as_ref().ok() != Some(&bound.context) {
            return;
        }
        match serde_json::to_vec_pretty(&bound.report) {
            Ok(bytes) => crate::io::download::download_with_completion(
                "recombination-reference-evidence.json",
                "application/json",
                &bytes,
                move |result| {
                    if status.try_get_untracked().is_none() {
                        return;
                    }
                    if let Err(e) = result {
                        status.set(e);
                    }
                },
            ),
            Err(e) => status.set(e.to_string()),
        }
    };
    view! {
        <div class="inspector-content recombination-panel">
            <h2>"Recombination"</h2>
            <p class="helper">"BP/LR reference exchange under an explicit full-site profile. All reciprocal components remain visible. End chemistry, host selection and experimental success require further modeling."</p>
            <p><strong>"Independent target: "</strong>{move || target.get().name}</p>
            <label>"First source"<select aria-label="Recombination first source" prop:value=move || first.get() on:change=move |e| first.set(event_target_value(&e))>
                <option value="">"Choose source"</option>
                {move || sources.get().into_iter().map(|s| view! {<option value=s.binding.document_id>{format!("{} · {} bases",s.name,s.length)}</option>}).collect_view()}
            </select></label>
            <label>"Second source"<select aria-label="Recombination second source" prop:value=move || second.get() on:change=move |e| second.set(event_target_value(&e))>
                <option value="">"Choose source"</option>
                {move || sources.get().into_iter().map(|s| view! {<option value=s.binding.document_id>{format!("{} · {} bases",s.name,s.length)}</option>}).collect_view()}
            </select></label>
            <label>"Reaction"<select aria-label="Recombination reaction" prop:value=move || if reaction.get()==reference::Reaction::BP {"BP"} else {"LR"} on:change=move |e| reaction.set(if event_target_value(&e)=="LR" {reference::Reaction::LR} else {reference::Reaction::BP})>
                <option value="BP">"BP (B/P reference exchange)"</option><option value="LR">"LR (L/R reference exchange)"</option>
            </select></label>
            <label class="checkbox-label"><input type="checkbox" aria-label="Use manufacturer recombination context profile" prop:checked=move || profile_enabled.get() on:change=move |e| profile_enabled.set(event_target_checked(&e))/>"Use manufacturer vector contexts · 2026-11-24.1"</label>
            <details><summary>"Profile and search limits"</summary>
                <p class="helper">"Eight full-site variants from the cited pDONR221, pDEST14, pENTR-D-TOPO and attB primer references. A missing context does not establish biochemical inactivity. MultiSite and other unregistered contexts remain unsupported."</p>
                {profile.get_value().ok().map(|p| view! {<p class="helper">{format!("{} · edition {} · content {}",p.profile.id,p.profile.edition,p.revision)}</p>})}
                <label>"Comparison budget"<input aria-label="Recombination comparison budget" type="text" inputmode="numeric" prop:value=move || comparisons.get() on:input=move |e| comparisons.set(event_target_value(&e))/></label>
                <label>"Candidate limit"<input aria-label="Recombination candidate limit" type="text" inputmode="numeric" prop:value=move || candidates.get() on:input=move |e| candidates.set(event_target_value(&e))/></label>
                <label>"Cumulative copied-base budget"<input aria-label="Recombination copied-base budget" type="text" inputmode="numeric" prop:value=move || product_bases.get() on:input=move |e| product_bases.set(event_target_value(&e))/></label>
                <p class="helper">"Other bounds: 1,000,000 source bases; 256 sites; 4,096 site combinations; four components per candidate; 16 MiB reference report; bounded annotation projection. Limits can leave the search incomplete."</p>
            </details>
            <button class="button full" disabled=move || busy.get() || context.get().is_err() on:click=run>"Reconstruct reciprocal components"</button>
            <Show when=move || busy.get()><button class="button secondary full" on:click=cancel>"Cancel recombination"</button></Show>
            <p class="helper" role="status">{move || if status.get().is_empty() {context.get().err().unwrap_or_default()} else {status.get()}}</p>
            {move || report.get().map(|bound| {
                let reference=&bound.report.reference;
                view! {<section>
                    <h3>"Reciprocal components"</h3>
                    <p class="helper">{format!("{} exact/possible site hits · {} attempted combinations · {} materialized bases. Annotation work: {} units.",reference.sites.len(),reference.combinations,reference.materialized_bases,bound.report.annotation_work)}</p>
                    <ul>{reference.diagnostics.iter().chain(bound.report.warnings.iter()).map(|w| view! {<li>{w.clone()}</li>}).collect_view()}</ul>
                    {move || report.get().map(|b| b.report.documents.iter().skip(page.get()*8).take(8).map(|doc| {
                        let choice=doc.selection.clone(); let id=choice.clone();
                        let label=format!("Candidate {}, component {} · {} bases · {:?}",choice.candidate_id+1,choice.component_id+1,doc.molecule.sequence.len(),doc.molecule.topology);
                        let warnings=doc.warnings.clone();
                        view! {<div class="candidate-review"><label class="checkbox-label"><input type="radio" name="recombination-component" aria-label=label.clone() prop:checked=move || selection.get().as_ref()==Some(&id) on:change=move |_| selection.set(Some(choice.clone()))/>{label}</label>
                            <p class="helper">{format!("{} projected annotations · {} reciprocal components retained",doc.molecule.features.len(),doc.reciprocal_components.len())}</p>
                            <details><summary>"Annotation and model findings"</summary><ul>{warnings.into_iter().map(|w| view! {<li>{w.clone()}</li>}).collect_view()}</ul>
                            <ul>{doc.annotation_evidence.iter().map(|e| view! {<li>{format!("{}: {} mapped, {} reviewed, {} dropped, {} retained as unmapped original evidence. {}",e.source.role,e.mapped_features,e.reviewed_features,e.dropped_features,e.unmapped_features.len(),e.warnings.join(" "))}</li>}).collect_view()}</ul></details>
                        </div>}
                    }).collect_view())}
                    <div class="plan-exchange"><button class="button secondary" disabled=move || page.get()==0 on:click=move |_| page.update(|p| *p=p.saturating_sub(1))>"Previous components"</button><button class="button secondary" disabled=move || report.get().is_none_or(|b| (page.get()+1)*8>=b.report.documents.len()) on:click=move |_| page.update(|p| *p+=1)>"Next components"</button></div>
                    <button class="button secondary full" disabled=move || busy.get() on:click=move |_| check_target.run(())>"Check independent target"</button>
                    <button class="button secondary full" disabled=move || busy.get() || selection.get().is_none() on:click=plan>"Review selected component in Plan"</button>
                    <button class="button secondary full" on:click=export>"Download reference evidence JSON"</button>
                </section>}
            })}
            {move || assessment.get().map(|a| view! {<section class="recombination-target-check"><h3>"Independent target check"</h3>
                <p>{format!("{} matches · search complete: {} · applicability: {:?}",a.reference.matches.len(),a.reference.search_complete,a.reference.applicability)}</p>
                <p class="helper">{format!("Possible sites: {} · unsupported geometry: {} · ambiguous source DNA: {} · ambiguous target DNA: {} · absence proven within the declared model: {}",a.reference.possible_sites,a.reference.unsupported_geometry,a.reference.ambiguous_source_dna,a.reference.ambiguous_target_dna,a.reference.no_match_proven_under_model)}</p>
                <ul>{a.reference.matches.iter().map(|m| view! {<li>{format!("Candidate {}, component {} · reverse complement: {} · fully specified DNA: {}",m.candidate+1,m.product+1,m.reverse,m.fully_specified_dna)}</li>}).collect_view()}</ul>
                <ul>{a.reference.diagnostics.iter().chain(a.reference.warnings.iter()).map(|w| view! {<li>{w.clone()}</li>}).collect_view()}</ul>
            </section>})}
        </div>
    }
}
