//! Native whole-graph discovery; serialized reports are displayed, never stock.
use crate::{
    topo_preparation::{Draft, SourceMode},
    worker::CachedWorker,
};
use collagene_core::{
    protocol::{AnalysisJob, AnalysisResult, Request},
    topo_graph_history as graph,
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
#[component]
pub fn GraphSourcePicker(
    draft: RwSignal<Draft>,
    epoch: RwSignal<u64>,
    locked: Memo<bool>,
    snapshot: Callback<(), Result<Workspace, String>>,
) -> impl IntoView {
    let records = RwSignal::new(Vec::<graph::Summary>::new());
    let inspected = RwSignal::new(None::<Arc<graph::Inspection>>);
    let pending = RwSignal::new(false);
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
    let key = Memo::new(move |_| {
        (
            draft.with(|d| (d.source_mode, d.graph.clone())),
            epoch.get(),
        )
    });
    Effect::new(move || {
        let _ = key.get();
        generation.update(|g| *g = g.wrapping_add(1));
        if let Some(w) = worker.get_value() {
            w.close();
        }
        worker.set_value(None);
        pending.set(false);
        inspected.set(None);
    });
    let list_alive = alive.clone();
    let list = move |_| {
        if pending.get_untracked() || locked.get_untracked() {
            return;
        }
        let start = (|| -> Result<_, String> {
            let frame = snapshot.run(())?;
            let revision = graph::list_revision(&frame).map_err(|e| e.to_string())?;
            Ok((frame, revision, Rc::new(CachedWorker::new()?)))
        })();
        let (frame, revision, w) = match start {
            Ok(v) => v,
            Err(e) => {
                status.set(e);
                return;
            }
        };
        let captured = key.get_untracked();
        let token = generation.get_untracked().wrapping_add(1);
        generation.set(token);
        pending.set(true);
        worker.set_value(Some(w.clone()));
        let alive = list_alive.clone();
        spawn_local(async move {
            let result = async {
                let response = w
                    .ask(Request {
                        id: token,
                        revision: revision.clone(),
                        job: AnalysisJob::ListTopoGraphs { workspace: frame },
                    })
                    .await?;
                if !alive.load(Ordering::Relaxed)
                    || generation.get_untracked() != token
                    || key.get_untracked() != captured
                {
                    return Err("Graph source list context changed".to_string());
                }
                if response.id != token || response.revision != revision {
                    return Err("Graph list response context mismatch".into());
                }
                let result = match response.result {
                    AnalysisResult::TopoGraphsListed(value) => value,
                    AnalysisResult::Error(e) => return Err(e),
                    _ => return Err("Unexpected graph source list response".into()),
                };
                result
                    .validate_binding(&revision)
                    .map_err(|e| e.to_string())?;
                records.set(result.entries);
                Ok(())
            }
            .await;
            w.close();
            if alive.load(Ordering::Relaxed) && generation.get_untracked() == token {
                pending.set(false);
                worker.set_value(None);
                status.set(match result {
                    Ok(()) => "Choose the exact saved graph, then inspect its strands.".into(),
                    Err(e) => e,
                });
            }
        });
    };
    let inspect_alive = alive.clone();
    let inspect = move |_| {
        if pending.get_untracked() || locked.get_untracked() {
            return;
        }
        let start = (|| -> Result<_, String> {
            let selection = draft
                .with_untracked(|d| d.graph.clone())
                .ok_or("Choose a saved graph")?;
            let frame = snapshot.run(())?;
            let revision =
                graph::inspect_revision(&frame, &selection).map_err(|e| e.to_string())?;
            Ok((frame, selection, revision, Rc::new(CachedWorker::new()?)))
        })();
        let (frame, selection, revision, w) = match start {
            Ok(v) => v,
            Err(e) => {
                status.set(e);
                return;
            }
        };
        let captured = key.get_untracked();
        let token = generation.get_untracked().wrapping_add(1);
        generation.set(token);
        pending.set(true);
        worker.set_value(Some(w.clone()));
        let alive = inspect_alive.clone();
        spawn_local(async move {
            let result = async {
                let response = w
                    .ask(Request {
                        id: token,
                        revision: revision.clone(),
                        job: AnalysisJob::InspectTopoGraph {
                            workspace: frame,
                            selection: selection.clone(),
                        },
                    })
                    .await?;
                if !alive.load(Ordering::Relaxed)
                    || generation.get_untracked() != token
                    || key.get_untracked() != captured
                {
                    return Err("Graph source inspection context changed".to_string());
                }
                if response.id != token || response.revision != revision {
                    return Err("Graph source response context mismatch".into());
                }
                let result = match response.result {
                    AnalysisResult::TopoGraphInspected(value) => value,
                    AnalysisResult::Error(e) => return Err(e),
                    _ => return Err("Unexpected graph source inspection response".into()),
                };
                result
                    .validate_binding(&revision)
                    .map_err(|e| e.to_string())?;
                if result.selection != selection {
                    return Err("Another saved graph was returned".into());
                }
                let result = Arc::new(result);
                inspected.set(Some(result.clone()));
                Ok(())
            }
            .await;
            w.close();
            if alive.load(Ordering::Relaxed) && generation.get_untracked() == token {
                pending.set(false);
                worker.set_value(None);
                status.set(match result {
                    Ok(()) => "Saved graph replayed. Choose its existing strand handles.".into(),
                    Err(e) => e,
                });
            }
        });
    };
    view! {<details open><summary>"Saved graph source"</summary><button class="text-tool" disabled=move||locked.get()||pending.get() on:click=list>"Load saved graphs"</button><label>"Graph"<select disabled=move||locked.get()||pending.get() prop:value=move||draft.with(|d|d.graph.as_ref().map(|s|s.record_hash.clone()).unwrap_or_default()) on:change=move|e|{let hash=event_target_value(&e);let selection=records.with_untracked(|entries|entries.iter().find(|entry|entry.selection.record_hash==hash).map(|entry|entry.selection.clone()));draft.update(|d|{if d.graph!=selection{d.graph=selection;d.vector_strands=Default::default();d.insert_strands=Default::default();}});}>
    <option value="">"Choose a saved graph"</option>{move||records.get().into_iter().map(|entry|view!{<option value=entry.selection.record_hash.clone()>{format!("{} · {} calls · {}…",entry.name,entry.source_bases,entry.selection.record_hash.chars().take(10).collect::<String>())}</option>}).collect_view()}</select></label><button class="button secondary" disabled=move||locked.get()||pending.get()||draft.with(|d|d.graph.is_none()) on:click=inspect>"Inspect selected graph"</button>
    {move||inspected.get().map(|value|{let summary=value.summary.clone();view!{<p>{format!("{} source calls · {} chains · {} components · {} remaining activations",summary.source_bases,summary.chains,summary.components,summary.remaining_activations.len())}</p><p class="helper">"The entire graph is retained, including byproducts. A component view does not isolate material."</p><details><summary>"Activation and source evidence"</summary><p>{format!("Remaining: {}",summary.remaining_activations.join(", "))}</p><p>{format!("Released: {}",summary.released_activations.join(", "))}</p><p>{format!("Inherited limits present: {}",summary.incomplete)}</p>{summary.assumptions.into_iter().chain(summary.findings).map(|text|view!{<p>{text}</p>}).collect_view()}<p class="helper">"Recorded source and reference-check events:"</p>{summary.source_events.into_iter().map(|event|view!{<p style="overflow-wrap:anywhere;">{format!("{} · {} · {}",event.document_id,event.event_id,event.revision)}</p>}).collect_view()}</details>
    <GraphHandles draft inspected=value locked/>
    }})}<p role="status">{move||status.get()}</p></details>}
}
#[component]
fn GraphHandles(
    draft: RwSignal<Draft>,
    inspected: Arc<graph::Inspection>,
    locked: Memo<bool>,
) -> impl IntoView {
    let saved = StoredValue::new(inspected);
    view! {<div class="primer-grid">{[false,true].into_iter().map(|insert|view!{<Show when=move||!insert||draft.with(|d|d.source_mode==SourceMode::ResumeGraph)>{[0usize,1].into_iter().map(|index|view!{<label>{format!("Existing {} {} strand",if insert{"insert"}else{"vector"},if index==0{"top"}else{"bottom"})}<select disabled=move||locked.get() prop:value=move||draft.with(|d|if insert{d.insert_strands[index].clone()}else{d.vector_strands[index].clone()}) on:change=move|e|draft.update(|d|if insert{d.insert_strands[index]=event_target_value(&e)}else{d.vector_strands[index]=event_target_value(&e)})><option value="">"Choose an existing source strand"</option>{saved.with_value(|inspection|inspection.report.preparation.substrates.strands.iter().enumerate().map(|(i,strand)|view!{<option value=i.to_string()>{format!("{} · {} · {} calls",strand.role,strand.id,strand.sequence.len())}</option>}).collect_view())}</select></label>}).collect_view()}</Show>}).collect_view()}</div>}
}

pub struct CaptureCandidate {
    pub request: graph::CaptureRequest,
    pub guard: Callback<(), Result<(), String>>,
}
#[component]
pub fn GraphCaptureControls(
    capture: Callback<
        (crate::topo_graph_model::CaptureKind, String),
        Result<CaptureCandidate, String>,
    >,
    initial: Memo<bool>,
    pending_transfers: Memo<Option<usize>>,
    locked: Memo<bool>,
    snapshot: Callback<(), Result<Workspace, String>>,
    validate: Callback<Workspace, Result<(), String>>,
    on_restore: Callback<Workspace>,
    epoch: RwSignal<u64>,
    db: StoredValue<Option<crate::storage::Database>, LocalStorage>,
    workspace_busy: RwSignal<bool>,
    save_generation: RwSignal<u64>,
    message: RwSignal<String>,
) -> impl IntoView {
    use crate::topo_graph_model::{CaptureKind, Scope};
    use std::cell::{Cell, RefCell};
    let name = RwSignal::new(String::new());
    let prefix = RwSignal::new("0".to_string());
    let status = RwSignal::new(String::new());
    let pending = RwSignal::new(false);
    let nonce = StoredValue::new_local(Rc::new(Cell::new(0u64)));
    let worker = StoredValue::new_local(None::<Rc<CachedWorker>>);
    let tx = StoredValue::new_local(Rc::new(RefCell::new(None::<crate::storage::Transaction>)));
    let lease = StoredValue::new_local(Rc::new(crate::commit_lease::CommitLease::default()));
    let alive = Arc::new(AtomicBool::new(true));
    let cleanup = alive.clone();
    on_cleanup(move || {
        cleanup.store(false, Ordering::Relaxed);
        if let Some(w) = worker.get_value() {
            w.close();
        }
        if let Some(transaction) = tx.with_value(|t| t.borrow().clone()) {
            let _ = transaction.abort();
        } else if lease.with_value(|l| l.cancel_unstarted()) {
            workspace_busy.set(false);
        }
    });
    let cancel = move |_| {
        if let Some(transaction) = tx.with_value(|t| t.borrow().clone()) {
            status.set(match transaction.abort() {
                Ok(()) => "Save cancellation requested; waiting for storage.".into(),
                Err(e) => e,
            });
            return;
        }
        nonce.with_value(|n| n.set(n.get().wrapping_add(1)));
        if let Some(w) = worker.get_value() {
            w.close();
        }
        worker.set_value(None);
        if lease.with_value(|l| l.cancel_unstarted()) {
            workspace_busy.set(false);
        }
        pending.set(false);
        status.set("Graph save cancelled before storage.".into());
    };
    let save_alive = alive.clone();
    let save = Callback::new(move |kind: CaptureKind| {
        if pending.get_untracked() || locked.get_untracked() || workspace_busy.get_untracked() {
            return;
        }
        let prepared = (|| -> Result<_, String> {
            let candidate = capture.run((kind, name.get_untracked()))?;
            candidate.guard.run(())?;
            let frame = snapshot.run(())?;
            let scope = Scope::capture(&frame, &candidate.request, epoch.get_untracked())?;
            let database = db
                .get_value()
                .ok_or("Workspace recovery storage is unavailable")?;
            Ok((
                candidate,
                frame,
                scope,
                database,
                Rc::new(CachedWorker::new()?),
            ))
        })();
        let (candidate, frame, scope, database, w) = match prepared {
            Ok(v) => v,
            Err(e) => {
                status.set(e);
                return;
            }
        };
        let n = nonce.get_value();
        let token = n.get().wrapping_add(1);
        n.set(token);
        let lock = lease.get_value();
        if let Err(e) = lock.acquire(token) {
            status.set(e);
            return;
        }
        workspace_busy.set(true);
        save_generation.update(|g| *g = g.wrapping_add(1));
        pending.set(true);
        worker.set_value(Some(w.clone()));
        tx.with_value(|t| *t.borrow_mut() = None);
        let pending_tx = tx.get_value();
        let alive = save_alive.clone();
        spawn_local(async move {
            let scope = Arc::new(scope);
            let guard_scope = scope.clone();
            let source_guard = candidate.guard;
            let guard_alive = alive.clone();
            let guard_n = n.clone();
            let guard = move || -> Result<(), String> {
                if !guard_alive.load(Ordering::Relaxed) || guard_n.get() != token {
                    return Err("Graph save cancelled before commit".into());
                }
                source_guard.run(())?;
                guard_scope.permit(&snapshot.run(())?, epoch.get_untracked())
            };
            let result=async{
                guard()?;
                let response=w.ask(Request{id:token as u32,revision:scope.context_revision.clone(),job:AnalysisJob::CaptureTopoGraph{workspace:frame.clone(),request:candidate.request}}).await?;
                guard()?;
                if response.id!=token as u32||response.revision!=scope.context_revision{return Err("Graph capture response context mismatch".into());}
                let result=match response.result{AnalysisResult::TopoGraphCaptured(value)=>value,AnalysisResult::Error(e)=>return Err(e),_=>return Err("Unexpected native graph capture response".into())};
                result.receipt.validate_binding(&scope.context_revision).map_err(|e|e.to_string())?;
                if result.receipt.source_workspace_fingerprint!=scope.fingerprint{return Err("Graph capture source workspace changed".into());}
                let saved_name=result.receipt.summary.name.clone();
                let next=crate::workspace_commit::persist(&database,result.workspace,frame,validate,w.clone(),guard,pending_tx).await?;
                if message.try_get_untracked().is_some(){on_restore.run(next);message.set(format!("Graph saved: {saved_name}. Select it under saved graph source to continue. No sequence document was replaced."));}
                Ok::<_,String>(())
            }.await;
            w.close();
            if lock.release(token) && workspace_busy.try_get_untracked().is_some() {
                workspace_busy.set(false);
            }
            if alive.load(Ordering::Relaxed) && n.get() == token {
                pending.set(false);
                worker.set_value(None);
                status.set(match result {
                    Ok(()) => "Graph saved.".into(),
                    Err(e) => e,
                });
            }
        });
    });
    view! {<details><summary>"Save elementary graph"</summary><p class="helper">"Retain the entire declared graph and its source history. Component selection does not isolate material. Flap cleavage and ligase actions are separate and are not included in this elementary record."</p><label>"Graph name"<input maxlength="256" disabled=move||locked.get()||pending.get() prop:value=move||name.get() on:input=move|e|name.set(event_target_value(&e))/></label>
    <Show when=move||initial.get()><button class="button secondary" disabled=move||locked.get()||pending.get() on:click=move|_|save.run(CaptureKind::InitialVector)>"Save initial vector graph"</button></Show>
    <Show when=move||pending_transfers.get().is_some()><label>"New elementary transfers to retain"<select disabled=move||locked.get()||pending.get() prop:value=move||prefix.get() on:change=move|e|prefix.set(event_target_value(&e))>{move||(0..=pending_transfers.get().unwrap_or(0)).map(|n|view!{<option value=n.to_string()>{format!("First {n} pending transfers")}</option>}).collect_view()}</select></label><button class="button secondary" disabled=move||locked.get()||pending.get() on:click={move|_|{match prefix.get_untracked().parse::<usize>(){Ok(count)=>save.run(CaptureKind::Reaction{pending_transfers:count}),Err(_)=>status.set("Choose an explicit transfer prefix".into())}}}>"Save selected elementary state"</button></Show>
    <Show when=move||pending.get()><button class="text-tool" on:click=cancel>"Cancel graph save"</button></Show><p role="status">{move||status.get()}</p></details>}
}
