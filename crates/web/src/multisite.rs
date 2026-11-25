//! Ordered parent contributions, explicit full-context profiles and reference-only results.
use crate::{app::download, digest::SourceInfo, worker};
use collagene_core::{
    Molecule,
    multisite_recombination::{self as science, Request, Selection},
    protocol::{AnalysisJob, AnalysisResult, Request as WireRequest},
    recombination::{ProfileSnapshot, reference},
    run_lineage::{EventBinding, InputEventBinding},
    transfer::Candidate,
    verification::{DesignPlan, Operation, PlanInput, Source, TopoSource},
};
use leptos::prelude::*;
use std::sync::{
    Arc,
    atomic::{AtomicBool, Ordering},
};
use wasm_bindgen_futures::{JsFuture, spawn_local};
pub struct CaptureRequest {
    pub sources: Vec<EventBinding>,
    pub target: EventBinding,
    pub request: Request,
}
pub struct Captured {
    pub inputs: Vec<reference::Input>,
    pub target: Molecule,
}
struct Bound {
    key: String,
    sources: Vec<SourceInfo>,
    target: SourceInfo,
    captured: Captured,
    request: Request,
    batch: science::Batch,
}
fn initial_profile() -> Result<reference::multisite::Profile, String> {
    let old = ProfileSnapshot::manufacturer_context()
        .map_err(|e| e.to_string())?
        .profile;
    Ok(reference::multisite::Profile { id: old.id, edition: old.edition, sources: old.sources, variants: old.variants, limitations: vec!["Registered manufacturer contexts for specificities 1 and 2 only. This is not a complete MultiSite kit profile. Import a reviewed full-context profile for other sites.".into()] })
}
#[component]
pub fn MultiSitePanel(
    target: Memo<SourceInfo>,
    sources: Memo<Vec<SourceInfo>>,
    capture: Callback<CaptureRequest, Result<Captured, String>>,
    epoch: RwSignal<u64>,
    locked: Memo<bool>,
    on_candidate: Callback<(Candidate, Vec<InputEventBinding>)>,
) -> impl IntoView {
    let slots = RwSignal::new(vec![String::new(), String::new()]);
    let profile = RwSignal::new(initial_profile());
    let reaction = RwSignal::new(reference::Reaction::LR);
    let specificity = RwSignal::new("1,2".to_string());
    let acknowledge = RwSignal::new(false);
    let import_generation = RwSignal::new(0u32);
    let status = RwSignal::new(String::new());
    let busy = RwSignal::new(false);
    let generation = RwSignal::new(0u32);
    let cancellation = RwSignal::new(None::<Arc<AtomicBool>>);
    let alive = Arc::new(AtomicBool::new(true));
    {
        let alive = alive.clone();
        on_cleanup(move || {
            alive.store(false, Ordering::Relaxed);
            if let Some(job) = cancellation.get_untracked() {
                job.store(false, Ordering::Relaxed);
            }
        });
    }
    let draft = Memo::new_with_compare(
        move |_| -> Result<(String, Vec<SourceInfo>, SourceInfo, Request), String> {
            if locked.get() {
                return Err("Wait for the workspace operation to finish.".into());
            }
            if !acknowledge.get() {
                return Err("Review and accept the declared site-profile scope.".into());
            }
            let available = sources.get();
            let mut selected = Vec::new();
            for id in slots.get() {
                selected.push(
                    available
                        .iter()
                        .find(|s| s.binding.document_id == id)
                        .cloned()
                        .ok_or("Choose each original parent document.")?,
                );
            }
            let target = target.get();
            if selected
                .iter()
                .any(|s| s.binding.document_id == target.binding.document_id)
            {
                return Err("Use a separate active document for the intended target.".into());
            }
            let mut unique = std::collections::BTreeSet::new();
            if selected
                .iter()
                .any(|s| !unique.insert(s.binding.document_id.clone()))
            {
                return Err("Choose distinct parent documents.".into());
            }
            let ids = specificity
                .get()
                .split(',')
                .map(|s| {
                    s.trim()
                        .parse::<u8>()
                        .map_err(|_| "Enter comma-separated specificity numbers.".to_string())
                })
                .collect::<Result<Vec<_>, _>>()?;
            let mut limits = reference::Limits::default();
            limits.max_components_per_candidate = 32;
            let request = Request {
                reaction: reference::multisite::Request {
                    profile: profile.get()?,
                    reaction: reaction.get(),
                    specificities: ids,
                    limits,
                },
                annotation_limits: science::AnnotationLimits::default(),
            };
            science::request_preflight(&request).map_err(|e| e.to_string())?;
            let key = serde_json::to_string(&(
                epoch.get(),
                selected.iter().map(|s| &s.binding).collect::<Vec<_>>(),
                &target.binding,
                &request,
            ))
            .map_err(|e| e.to_string())?;
            Ok((key, selected, target, request))
        },
        |old, new| match (old, new) {
            (Some(Ok(a)), Some(Ok(b))) => a.0 != b.0,
            (Some(Err(a)), Some(Err(b))) => a != b,
            _ => true,
        },
    );
    let current =
        Callback::new(move |key: String| draft.get_untracked().as_ref().is_ok_and(|v| v.0 == key));
    let result = RwSignal::new(None::<Arc<Bound>>);
    let choice = RwSignal::new(None::<Selection>);
    let page = RwSignal::new(0usize);
    let checked = RwSignal::new(None::<Arc<science::CheckedProduct>>);
    Effect::new(move |_| {
        let _ = draft.get();
        if let Some(job) = cancellation.get_untracked() {
            job.store(false, Ordering::Relaxed);
        }
        generation.update(|g| *g = g.wrapping_add(1));
        result.set(None);
        choice.set(None);
        checked.set(None);
        busy.set(false);
        page.set(0);
    });
    let run_alive = alive.clone();
    let run = move |_| {
        let (key, selected, target, request) = match draft.get_untracked() {
            Ok(v) => v,
            Err(e) => {
                status.set(e);
                return;
            }
        };
        let captured = match capture.run(CaptureRequest {
            sources: selected.iter().map(|s| s.binding.clone()).collect(),
            target: target.binding.clone(),
            request: request.clone(),
        }) {
            Ok(v) => v,
            Err(e) => {
                status.set(e);
                return;
            }
        };
        let revision = match science::revision(&captured.inputs, &request) {
            Ok(v) => v,
            Err(e) => {
                status.set(e.to_string());
                return;
            }
        };
        let token = generation.get_untracked().wrapping_add(1);
        generation.set(token);
        busy.set(true);
        result.set(None);
        choice.set(None);
        checked.set(None);
        let job = Arc::new(AtomicBool::new(true));
        cancellation.set(Some(job.clone()));
        let alive = run_alive.clone();
        status.set(
            "Classifying all registered contexts and reconstructing whole reciprocal networks…"
                .into(),
        );
        spawn_local(async move {
            let response = worker::analyze_cancellable(
                WireRequest {
                    id: token,
                    revision: revision.clone(),
                    job: AnalysisJob::SimulateMultiSiteRecombination {
                        inputs: captured.inputs.clone(),
                        request: request.clone(),
                    },
                },
                job,
            )
            .await;
            if !alive.load(Ordering::Relaxed)
                || generation.get_untracked() != token
                || !current.run(key.clone())
            {
                return;
            }
            busy.set(false);
            match response {
                Ok(r) if r.id == token && r.revision == revision => match r.result {
                    AnalysisResult::MultiSiteRecombination(batch) => {
                        if let Err(e) = batch.validate_binding(&revision) {
                            status.set(e.to_string());
                            return;
                        }
                        status.set(format!("{} candidate networks; {} components. Search complete: {}. Select a component explicitly.",batch.reference.candidates.len(),batch.documents.len(),batch.reference.complete));
                        result.set(Some(Arc::new(Bound {
                            key,
                            sources: selected,
                            target,
                            captured,
                            request,
                            batch,
                        })));
                    }
                    AnalysisResult::Error(e) => status.set(e),
                    _ => status.set("Unexpected multi-parent response.".into()),
                },
                Ok(_) => status.set("Response binding changed.".into()),
                Err(e) => status.set(e),
            }
        });
    };
    let check_alive = alive.clone();
    let check = Callback::new(move |(): ()| {
        let (Some(bound), Some(selection)) = (result.get_untracked(), choice.get_untracked())
        else {
            return;
        };
        if !current.run(bound.key.clone()) {
            return;
        }
        let revision = match science::check_revision(
            &bound.captured.inputs,
            &bound.request,
            &selection,
            &bound.captured.target,
        ) {
            Ok(v) => v,
            Err(e) => {
                status.set(e.to_string());
                return;
            }
        };
        let token = generation.get_untracked().wrapping_add(1);
        generation.set(token);
        busy.set(true);
        checked.set(None);
        let job = Arc::new(AtomicBool::new(true));
        cancellation.set(Some(job.clone()));
        let alive = check_alive.clone();
        spawn_local(async move {
            let response = worker::analyze_cancellable(
                WireRequest {
                    id: token,
                    revision: revision.clone(),
                    job: AnalysisJob::CheckMultiSiteRecombination {
                        inputs: bound.captured.inputs.clone(),
                        request: bound.request.clone(),
                        selection: selection.clone(),
                        target: bound.captured.target.clone(),
                    },
                },
                job,
            )
            .await;
            if !alive.load(Ordering::Relaxed)
                || generation.get_untracked() != token
                || !current.run(bound.key.clone())
                || choice.get_untracked() != Some(selection.clone())
            {
                return;
            }
            busy.set(false);
            match response {
                Ok(r) if r.id == token && r.revision == revision => match r.result {
                    AnalysisResult::MultiSiteRecombinationChecked(value) => {
                        if value.context_revision().ok().as_ref() != Some(&revision) {
                            status.set("Target-check binding changed.".into());
                            return;
                        }
                        status.set("Fresh independent-target comparison complete.".into());
                        checked.set(Some(Arc::new(value)));
                    }
                    AnalysisResult::Error(e) => status.set(e),
                    _ => status.set("Unexpected target response.".into()),
                },
                Ok(_) => status.set("Target response binding changed.".into()),
                Err(e) => status.set(e),
            }
        });
    });
    let plan = move |_| {
        let (Some(bound), Some(selection)) = (result.get_untracked(), choice.get_untracked())
        else {
            return;
        };
        if busy.get_untracked()
            || !current.run(bound.key.clone())
            || !bound
                .batch
                .documents
                .iter()
                .any(|d| d.selection == selection)
        {
            return;
        }
        let inputs = bound
            .captured
            .inputs
            .iter()
            .enumerate()
            .map(|(i, s)| PlanInput {
                id: format!("parent_{}", i + 1),
                revision: s.molecule.revision(),
                molecule: s.molecule.clone(),
            })
            .collect();
        let bindings = bound
            .sources
            .iter()
            .enumerate()
            .map(|(i, s)| InputEventBinding {
                input_id: format!("parent_{}", i + 1),
                event: s.binding.clone(),
            })
            .collect();
        let sources = bound
            .captured
            .inputs
            .iter()
            .enumerate()
            .map(|(i, s)| TopoSource {
                role: s.role.clone(),
                source: Source::Input(format!("parent_{}", i + 1)),
            })
            .collect();
        on_candidate.run((
            Candidate::new(DesignPlan {
                name: format!(
                    "Multi-parent {:?} reference exchange to {}",
                    bound.request.reaction.reaction, bound.target.name
                ),
                target_revision: bound.target.binding.revision.clone(),
                inputs,
                steps: vec![Operation::MultiSiteRecombineReference {
                    sources,
                    request: bound.request.clone(),
                    selection,
                }],
            }),
            bindings,
        ));
    };
    let import_profile = move |event| {
        let input = event_target::<web_sys::HtmlInputElement>(&event);
        let Some(file) = input.files().and_then(|files| files.get(0)) else {
            return;
        };
        if file.size() > 2.0 * 1024.0 * 1024.0 {
            status.set("Site profiles must be at most 2 MiB.".into());
            return;
        }
        acknowledge.set(false);
        import_generation.update(|g| *g = g.wrapping_add(1));
        let import_token = import_generation.get_untracked();
        let alive = alive.clone();
        spawn_local(async move {
            let value = JsFuture::from(file.text())
                .await
                .map_err(|_| "Unable to read profile.".to_string())
                .and_then(|v| v.as_string().ok_or("Profile text is unavailable.".into()))
                .and_then(|s| {
                    serde_json::from_str::<reference::multisite::Profile>(&s)
                        .map_err(|e| e.to_string())
                })
                .and_then(|p| {
                    p.validate().map_err(|e| e.to_string())?;
                    Ok(p)
                });
            if !alive.load(Ordering::Relaxed) || import_generation.get_untracked() != import_token {
                return;
            }
            match value {
                Ok(p) => {
                    acknowledge.set(false);
                    profile.set(Ok(p));
                    status.set(
                        "Imported profile. Review its sources and limitations before accepting."
                            .into(),
                    );
                }
                Err(e) => status.set(e),
            }
        });
    };
    view! {<section class="inspector-content">
        <h2>"Multi-parent reference exchange"</h2><p class="helper">"Choose original parents, a declared full-context site profile and an independent intended target. This predicts reference DNA and all reciprocal components; chemical strand state and experimental selection remain unestablished."</p>
        <p><strong>"Independent target: "</strong>{move||target.get().name}</p>
        {move||slots.get().into_iter().enumerate().map(|(index,id)|view!{<label>{format!("Parent {}",index+1)}<select prop:value=id on:change=move|e|slots.update(|v|v[index]=event_target_value(&e))><option value="">"Choose original parent"</option>{move||sources.get().into_iter().map(|s|view!{<option value=s.binding.document_id>{format!("{} · {} bases",s.name,s.length)}</option>}).collect_view()}</select></label>}).collect_view()}
        <div class="plan-exchange"><button class="button secondary" disabled=move||slots.with(|v|v.len()>=5)||busy.get() on:click=move |_|slots.update(|v|v.push(String::new()))>"Add parent"</button><button class="button secondary" disabled=move||slots.with(|v|v.len()<=2)||busy.get() on:click=move |_|slots.update(|v|{v.pop();})>"Remove last parent"</button></div>
        <label>"Reaction"<select on:change=move|e|reaction.set(if event_target_value(&e)=="BP"{reference::Reaction::BP}else{reference::Reaction::LR})><option value="LR">"LR · L/R exchange"</option><option value="BP">"BP · B/P exchange"</option></select></label>
        <label>"Specificity numbers"<input type="text" prop:value=move||specificity.get() on:input=move|e|specificity.set(event_target_value(&e))/></label>
        <details><summary>"Site-profile sources and limitations"</summary>{move||profile.get().ok().map(|p|view!{<p>{format!("{} · {} · {} registered contexts",p.id,p.edition,p.variants.len())}</p><ul>{p.limitations.into_iter().map(|v|view!{<li>{v}</li>}).collect_view()}</ul><ul>{p.sources.into_iter().map(|s|view!{<li><a href=s.url target="_blank" rel="noopener noreferrer">{s.vector_name}</a></li>}).collect_view()}</ul>})}<label>"Import reviewed full-context site profile"<input type="file" accept=".json,application/json" on:change=import_profile/></label></details>
        <label class="checkbox-label"><input type="checkbox" prop:checked=move||acknowledge.get() on:change=move|e|acknowledge.set(event_target_checked(&e))/>"Use this declared site-profile scope"</label>
        <button class="button full" disabled=move||busy.get()||draft.get().is_err() on:click=run>"Reconstruct reciprocal networks"</button>
        <Show when=move||busy.get()><button class="button secondary full" on:click=move |_|{if let Some(job)=cancellation.get_untracked(){job.store(false,Ordering::Relaxed);}generation.update(|g|*g=g.wrapping_add(1));busy.set(false);status.set("Computation cancelled.".into());}>"Cancel"</button></Show>
        <p class="helper" role="status">{move||if status.get().is_empty(){draft.get().err().unwrap_or_default()}else{status.get()}}</p>
        {move||result.get().map(|bound|view!{<section><h3>"Candidate networks and reciprocal components"</h3><p>{format!("{} registered site hits · {} explored assignments · {} comparisons",bound.batch.reference.sites.len(),bound.batch.reference.combinations,bound.batch.reference.comparisons)}</p><ul>{bound.batch.reference.warnings.iter().chain(&bound.batch.reference.limits).chain(&bound.batch.reference.unsupported).chain(&bound.batch.warnings).map(|s|view!{<li>{s.clone()}</li>}).collect_view()}</ul>
            {move||result.get().map(|b|b.batch.documents.iter().skip(page.get()*8).take(8).map(|d|{let selection=d.selection.clone();let id=selection.clone();view!{<div class="candidate-review"><label class="checkbox-label"><input type="radio" name="multisite-component" disabled=move||busy.get() prop:checked=move||choice.get().as_ref()==Some(&id) on:change=move |_|{choice.set(Some(selection.clone()));checked.set(None);}/>{format!("Candidate {}, component {} · {} bases · {:?}",d.selection.candidate_id+1,d.selection.component_id+1,d.molecule.sequence.len(),d.molecule.topology)}</label><p class="helper">{format!("{} reciprocal components retained · {} mapped annotations",d.reciprocal_components.len(),d.molecule.features.len())}</p><ul>{d.warnings.iter().map(|w|view!{<li>{w.clone()}</li>}).collect_view()}</ul></div>}}).collect_view())}
            <div class="plan-exchange"><button class="button secondary" disabled=move||page.get()==0 on:click=move |_|page.update(|p|*p=p.saturating_sub(1))>"Previous"</button><button class="button secondary" disabled=move||result.get().is_none_or(|b|(page.get()+1)*8>=b.batch.documents.len()) on:click=move |_|page.update(|p|*p+=1)>"Next"</button></div>
            <button class="button secondary full" disabled=move||busy.get()||choice.get().is_none() on:click=move |_|check.run(())>"Check selected component against target"</button><button class="button secondary full" disabled=move||busy.get()||choice.get().is_none() on:click=plan>"Review selected component in Plan"</button>
            <button class="button secondary full" on:click=move |_|{if let Some(b)=result.get_untracked(){if current.run(b.key.clone()){match serde_json::to_vec_pretty(&b.batch){Ok(bytes)=>if let Err(e)=download("multisite-reference-evidence.json","application/json",&bytes){status.set(e);},Err(e)=>status.set(e.to_string())}}}}>"Download full reference evidence"</button>
        </section>})}
        {move||checked.get().map(|r|view!{<section><h3>"Independent target comparison"</h3><p>{format!("Reference match: {} · complete within declared model: {}",r.assessment.target_matches,r.assessment.complete)}</p><ul>{r.assessment.findings.iter().map(|s|view!{<li>{s.clone()}</li>}).collect_view()}</ul></section>})}
    </section>}
}
