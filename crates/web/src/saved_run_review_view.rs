//! Exact saved-plan review and explicit retained-plan editing handoff.
use crate::{
    saved_run_review_model::{self as model, Scope},
    worker,
};
use collagene_core::{
    protocol::{AnalysisJob, AnalysisResult, Request},
    run_lineage,
    run_review::PreparedReview,
    verification::{PhysicalProcessingEvidence, VerificationReport},
    workspace::Workspace,
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

fn scope(
    snapshot: Callback<(), Result<Workspace, String>>,
    run_id: &str,
    epoch: u64,
) -> Result<(Scope, Workspace), String> {
    let frame = snapshot.run(())?;
    let fingerprint = run_lineage::workspace_fingerprint(&frame).map_err(|e| e.to_string())?;
    let summary = run_lineage::list_runs(&frame)
        .map_err(|e| e.to_string())?
        .into_iter()
        .find(|r| r.run_id == run_id)
        .ok_or("The exact retained run is unavailable.")?;
    Ok((
        Scope {
            epoch,
            workspace_fingerprint: fingerprint,
            run_id: run_id.into(),
            run_attachment_hash: summary.attachment_hash,
        },
        frame,
    ))
}
#[component]
pub fn SavedRunReviewPanel(
    snapshot: Callback<(), Result<Workspace, String>>,
    run_id: String,
    epoch: RwSignal<u64>,
    on_edit_plan: Callback<Arc<PreparedReview>, Result<(), String>>,
) -> impl IntoView {
    let run_id = StoredValue::new(run_id);
    let opened_epoch = epoch.get_untracked();
    let step_page = RwSignal::new(0usize);
    let product_page = RwSignal::new(0usize);
    let pending = RwSignal::new(false);
    let status = RwSignal::new(String::new());
    let generation = RwSignal::new(0u64);
    let prepared = RwSignal::new(None::<Arc<PreparedReview>>);
    let report = RwSignal::new(None::<Arc<VerificationReport>>);
    let bound = StoredValue::new(None::<Scope>);
    let executor = StoredValue::new_local(None::<Rc<worker::CachedWorker>>);
    let alive = Arc::new(AtomicBool::new(true));
    let cleanup = alive.clone();
    on_cleanup(move || {
        cleanup.store(false, Ordering::Relaxed);
        if let Some(w) = executor.get_value() {
            w.close();
        }
        executor.set_value(None);
    });
    Effect::new(move || {
        let _ = epoch.get();
        generation.update(|g| *g = g.wrapping_add(1));
        if let Some(w) = executor.get_value() {
            w.close();
        }
        executor.set_value(None);
        prepared.set(None);
        report.set(None);
        pending.set(false);
        bound.set_value(None);
    });
    let review_alive = alive.clone();
    let review = move |_| {
        if pending.get_untracked() || epoch.get_untracked() != opened_epoch {
            return;
        }
        let (captured, frame) = match scope(snapshot, &run_id.get_value(), opened_epoch) {
            Ok(v) => v,
            Err(e) => {
                status.set(e);
                return;
            }
        };
        let request_context =
            match collagene_core::run_review::request_revision(&frame, &captured.run_id) {
                Ok(value) => value,
                Err(error) => {
                    status.set(error.to_string());
                    return;
                }
            };
        let w = match worker::CachedWorker::new() {
            Ok(w) => Rc::new(w),
            Err(e) => {
                status.set(e);
                return;
            }
        };
        if let Some(old) = executor.get_value() {
            old.close();
        }
        executor.set_value(Some(w.clone()));
        generation.update(|g| *g = g.wrapping_add(1));
        let ticket = generation.get_untracked();
        pending.set(true);
        prepared.set(None);
        report.set(None);
        bound.set_value(None);
        status.set("Preparing the exact retained plan and independent target, then checking them with the current verifier…".into());
        let alive = review_alive.clone();
        spawn_local(async move {
            let result: Result<(PreparedReview, VerificationReport), String> = async {
                let response = w
                    .ask(Request {
                        id: 1,
                        revision: request_context.clone(),
                        job: AnalysisJob::PrepareSavedRunReview {
                            workspace: frame,
                            run_id: captured.run_id.clone(),
                        },
                    })
                    .await?;
                if !alive.load(Ordering::Relaxed)
                    || generation.get_untracked() != ticket
                    || epoch.get_untracked() != opened_epoch
                {
                    return Err("Review was cancelled or its workspace changed.".into());
                }
                let (current, frame) = scope(snapshot, &captured.run_id, epoch.get_untracked())?;
                captured.permit(&current)?;
                if response.id != 1 || response.revision != request_context {
                    return Err("Preparation response belongs to another request.".into());
                }
                let p = match response.result {
                    AnalysisResult::SavedRunReviewPrepared(p) => p,
                    AnalysisResult::Error(e) => return Err(e),
                    _ => return Err("Unexpected saved-run preparation response.".into()),
                };
                captured.validate_prepared(&p)?;
                p.validate_for_workspace(&frame)
                    .map_err(|e| e.to_string())?;
                let response = w
                    .ask(Request {
                        id: 2,
                        revision: p.context_revision.clone(),
                        job: AnalysisJob::Verify {
                            plan: p.candidate.plan.clone(),
                            target: p.target.clone(),
                        },
                    })
                    .await?;
                if response.id != 2 || response.revision != p.context_revision {
                    return Err("Verification response belongs to another review.".into());
                }
                let r = match response.result {
                    AnalysisResult::Verification(r) => r,
                    AnalysisResult::Error(e) => return Err(e),
                    _ => return Err("Unexpected saved-plan verification response.".into()),
                };
                model::validate_report(&p, &r)?;
                Ok((p, r))
            }
            .await;
            w.close();
            if !alive.load(Ordering::Relaxed)
                || generation.get_untracked() != ticket
                || epoch.get_untracked() != opened_epoch
            {
                return;
            }
            executor.set_value(None);
            pending.set(false);
            let current = scope(snapshot, &captured.run_id, epoch.get_untracked())
                .and_then(|(s, _)| captured.permit(&s));
            if let Err(e) = current {
                status.set(e);
                return;
            }
            match result {
                Ok((p, r)) => {
                    bound.set_value(Some(captured));
                    prepared.set(Some(Arc::new(p)));
                    report.set(Some(Arc::new(r)));
                    status.set("Fresh review completed. Editor heads and recorded historical status were not changed.".into());
                }
                Err(e) => status.set(e),
            }
        });
    };
    let cancel = move |_| {
        generation.update(|g| *g = g.wrapping_add(1));
        if let Some(w) = executor.get_value() {
            w.close();
        }
        executor.set_value(None);
        pending.set(false);
        prepared.set(None);
        report.set(None);
        bound.set_value(None);
        status.set("Review cancelled.".into());
    };
    let capture = Callback::new(move |()| {
        if epoch.get_untracked() != opened_epoch {
            return Err("Workspace changed; reopen this saved run.".into());
        }
        let p = prepared
            .get_untracked()
            .ok_or("Prepare and verify the saved plan first.")?;
        let r = report
            .get_untracked()
            .ok_or("No current verification receipt.")?;
        let expected = bound.get_value().ok_or("Review scope is unavailable.")?;
        let (current, frame) = scope(snapshot, &expected.run_id, epoch.get_untracked())?;
        expected.permit(&current)?;
        p.validate_for_workspace(&frame)
            .map_err(|e| e.to_string())?;
        model::validate_report(&p, &r)?;
        Ok((p.candidate.plan.clone(), p.target.clone()))
    });
    let edit = move |_| {
        let result = (|| -> Result<(), String> {
            if epoch.get_untracked() != opened_epoch || pending.get_untracked() {
                return Err("Workspace changed or review is still running.".into());
            }
            let p = prepared
                .get_untracked()
                .ok_or("Prepare the exact saved plan first.")?;
            let expected = bound.get_value().ok_or("Review scope is unavailable.")?;
            let (current, frame) = scope(snapshot, &expected.run_id, epoch.get_untracked())?;
            model::editing_candidate(&expected, &current, &p, &frame)?;
            on_edit_plan.run(p)
        })();
        if let Err(error) = result {
            status.set(error);
        }
    };
    let identity = Memo::new(move |_| {
        prepared
            .get()
            .map(|p| {
                (
                    p.candidate.plan_revision.clone(),
                    p.target_event.revision.clone(),
                )
            })
            .unwrap_or_default()
    });
    view! {<details class="saved-run-review"><summary>"Review saved plan and export instructions"</summary>
        <p class="helper">"Uses this run’s retained independent target and exact original source events, even if current documents have changed. No document is rewound or overwritten. Historical status is not adopted."</p>
        <button class="button secondary" disabled=move||pending.get()||epoch.get()!=opened_epoch on:click=review>"Prepare and verify saved plan"</button>
        <button class="button secondary" disabled={move||pending.get()||prepared.get().is_none()||epoch.get()!=opened_epoch} on:click=edit>"Edit retained plan"</button>
        <Show when=move||pending.get()><button class="text-tool" on:click=cancel>"Cancel review"</button></Show>
        <p role="status">{move||status.get()}</p>
        {move||prepared.get().zip(report.get()).map(|(p,r)|{
            let status_text=format!("Current verification: {:?} · recorded engine matches: {}",r.status,p.recorded_engine_matches);
            let target=format!("Intended target: {} · document {} · event {} · revision {}",p.target.name,p.target_event.document_id,p.target_event.event_id,p.target_event.revision);
            let findings=r.findings.clone();let assumptions=r.assumptions.clone();let warnings=p.warnings.clone();
            let steps=p.candidate.plan.steps.iter().enumerate().skip(step_page.get()*25).take(25).map(|(i,op)|view!{<details><summary>{format!("Step {}",i+1)}</summary><pre style="white-space:pre-wrap;overflow-wrap:anywhere">{crate::app::describe_step(op)}</pre></details>}).collect_view();
            let products=r.products.iter().enumerate().skip(product_page.get()*25).take(25).map(|(i,m)|view!{<details><summary>{format!("Reconstructed step {} · {} bases",i+1,m.sequence.len())}</summary><pre style="white-space:pre-wrap;overflow-wrap:anywhere">{m.sequence.as_str().to_owned()}</pre></details>}).collect_view();
            let outcomes=r.physical_processing_reports.iter().filter_map(|step|match &step.evidence{PhysicalProcessingEvidence::LigationOutcomes(e)=>Some(crate::ligation_outcomes_view::ligation_outcomes_view(step.step_index,e)),_=>None}).collect_view();
            view!{<section><strong>{status_text}</strong><p>{target}</p>{warnings.into_iter().map(|w|view!{<p class="helper">{w}</p>}).collect_view()}{findings.into_iter().map(|f|view!{<p>{f}</p>}).collect_view()}
                <details><summary>"Current reconstruction assumptions"</summary>{assumptions.into_iter().map(|a|view!{<p>{a}</p>}).collect_view()}</details>
                <p class="helper">"Nested saved-material models may be unsupported by the current verifier. An incomplete or failed review remains diagnostic; exporting instructions does not establish physical reuse or experimental success."</p>
                <h5>{format!("{} saved operation steps",p.candidate.plan.steps.len())}</h5>{steps}
                <button class="text-tool" disabled=move||step_page.get()==0 on:click=move|_|step_page.update(|v|*v=v.saturating_sub(1))>"Previous steps"</button>
                <button class="text-tool" disabled=move||(step_page.get()+1)*25>=step_count on:click=move|_|step_page.update(|v|*v+=1)>"Next steps"</button>
                <h5>{format!("{} freshly reconstructed products",r.products.len())}</h5>{products}
                <button class="text-tool" disabled=move||product_page.get()==0 on:click=move|_|product_page.update(|v|*v=v.saturating_sub(1))>"Previous products"</button>
                <button class="text-tool" disabled=move||(product_page.get()+1)*25>=product_count on:click=move|_|product_page.update(|v|*v+=1)>"Next products"</button>
                {outcomes}<crate::assembly_outputs::AssemblyOutputs capture identity epoch/>
            </section>}
        })}
    </details>}
}
