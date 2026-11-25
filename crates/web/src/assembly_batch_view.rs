//! Cancellable selected-row evaluation. Only a private live worker can capture
//! fresh evidence; historical records replay explicitly and never restore a badge.
use crate::{assembly_batch_projection, assembly_draft_model::Draft, worker::CachedWorker};
use collagene_core::{
    combinatorial_batch::{self as batch, Record, RowEvidence, RowRequest},
    protocol::{AnalysisJob, AnalysisResult, Request},
    run_lineage,
    workspace::Workspace,
};
use leptos::prelude::*;
use std::{
    cell::{Cell, RefCell},
    rc::Rc,
    sync::{
        Arc,
        atomic::{AtomicBool, Ordering},
    },
};
use wasm_bindgen_futures::spawn_local;

#[derive(Clone)]
enum Origin {
    Current(String),
    Saved(String),
}
#[derive(Clone)]
struct Binding {
    context: String,
    epoch: u64,
    origin: Origin,
}
fn row_id(row: &RowRequest) -> &str {
    match row {
        RowRequest::Ready(row) => &row.id,
        RowRequest::DraftError { id, .. } => id,
    }
}

#[component]
pub fn AssemblyBatchView(
    draft: RwSignal<Draft>,
    epoch: RwSignal<u64>,
    snapshot: Callback<(), Result<Workspace, String>>,
    validate: Callback<Workspace, Result<(), String>>,
    on_restore: Callback<Workspace>,
    attachment_available: Callback<String, bool>,
    db: StoredValue<Option<web_sys::IdbDatabase>, LocalStorage>,
    workspace_busy: RwSignal<bool>,
    loaded: RwSignal<bool>,
    save_generation: RwSignal<u64>,
    message: RwSignal<String>,
) -> impl IntoView {
    let key = Memo::new(move |_| {
        draft.with(|value| assembly_batch_projection::project(value).map(|p| p.revision))
    });
    let scope = RwSignal::new(None::<Arc<Binding>>);
    let rows = RwSignal::new(Vec::<RowEvidence>::new());
    let total = RwSignal::new(0usize);
    let evaluating = RwSignal::new(false);
    let capturing = RwSignal::new(false);
    let finished = RwSignal::new(false);
    let historical = RwSignal::new(false);
    let status = RwSignal::new(String::new());
    let history = RwSignal::new(Vec::<batch::Summary>::new());
    let historical_record = RwSignal::new(None::<Arc<Record>>);
    let worker = StoredValue::new_local(None::<Rc<CachedWorker>>);
    let generation = StoredValue::new_local(Rc::new(Cell::new(0u64)));
    let transaction =
        StoredValue::new_local(Rc::new(RefCell::new(None::<web_sys::IdbTransaction>)));
    let owns_lock = StoredValue::new_local(Rc::new(crate::commit_lease::CommitLease::default()));
    let alive = Arc::new(AtomicBool::new(true));
    let cancel = Callback::new(move |()| {
        generation.with_value(|g| g.set(g.get().wrapping_add(1)));
        if let Some(w) = worker.get_value() {
            w.close();
        }
        worker.set_value(None);
        if let Some(tx) = transaction.with_value(|value| value.borrow().clone()) {
            let _ = tx.abort();
        } else if owns_lock.with_value(|value| value.cancel_unstarted()) {
            workspace_busy.set(false);
        }
        scope.set(None);
        evaluating.set(false);
        capturing.set(false);
        finished.set(false);
        status.set("Batch operation cancelled. Completed rows remain display-only; no partial record was installed.".into());
    });
    let cleanup_alive = alive.clone();
    on_cleanup(move || {
        cleanup_alive.store(false, Ordering::Relaxed);
        generation.with_value(|g| g.set(g.get().wrapping_add(1)));
        if let Some(w) = worker.get_value() {
            w.close();
        }
        if let Some(tx) = transaction.with_value(|value| value.borrow().clone()) {
            let _ = tx.abort();
        } else if owns_lock.with_value(|value| value.cancel_unstarted()) {
            workspace_busy.set(false);
        }
    });
    let current = Callback::new(move |binding: Arc<Binding>| -> Result<(), String> {
        if epoch.get_untracked() != binding.epoch {
            return Err("The workspace was restored; rerun the batch.".into());
        }
        match &binding.origin {
            Origin::Current(expected) if key.get_untracked().as_ref() != Ok(expected) => {
                return Err(
                    "Selected constructs or editable settings changed; rerun the batch.".into(),
                );
            }
            Origin::Saved(hash) if !attachment_available.run(hash.clone()) => {
                return Err("The historical batch attachment is no longer retained.".into());
            }
            _ => {}
        }
        Ok(())
    });
    Effect::new(move || {
        let now = epoch.get();
        let declaration = key.get();
        if let Some(binding) = scope.get_untracked() {
            if now != binding.epoch
                || matches!(&binding.origin,Origin::Current(expected) if declaration.as_ref()!=Ok(expected))
            {
                cancel.run(());
                status.set(
                    "The batch context changed. Run again; previous rows are display-only.".into(),
                );
            }
        }
    });
    let start = {
        let alive = alive.clone();
        Callback::new(move |(matrix, origin): (batch::Request, Origin)| {
            if evaluating.get_untracked() || workspace_busy.get_untracked() {
                return;
            }
            let prepared = (|| -> Result<_, String> {
                let frame = snapshot.run(())?;
                let revision = batch::begin_revision(&frame, &matrix).map_err(|e| e.to_string())?;
                let binding = Arc::new(Binding {
                    context: revision.clone(),
                    epoch: epoch.get_untracked(),
                    origin,
                });
                current.run(binding.clone())?;
                let w = Rc::new(CachedWorker::new()?);
                Ok((frame, revision, binding, w))
            })();
            let (frame, revision, binding, w) = match prepared {
                Ok(v) => v,
                Err(e) => {
                    status.set(e);
                    return;
                }
            };
            if let Some(old) = worker.get_value() {
                old.close();
            }
            let nonce = generation.get_value();
            let token = nonce.get().wrapping_add(1);
            nonce.set(token);
            worker.set_value(Some(w.clone()));
            scope.set(Some(binding.clone()));
            evaluating.set(true);
            finished.set(false);
            historical.set(false);
            historical_record.set(None);
            rows.set(vec![]);
            total.set(matrix.rows.len());
            status.set("Starting selected-construct evaluation in the Rust worker…".into());
            let alive = alive.clone();
            spawn_local(async move {
                let guard = || -> Result<(), String> {
                    if !alive.load(Ordering::Relaxed) || nonce.get() != token {
                        return Err("Batch cancelled.".into());
                    }
                    current.run(binding.clone())
                };
                let result = async {
                    let response = w
                        .ask(Request {
                            id: token as u32,
                            revision: revision.clone(),
                            job: AnalysisJob::BeginCombinatorialBatch {
                                workspace: frame,
                                matrix: matrix.clone(),
                            },
                        })
                        .await?;
                    guard()?;
                    if response.id != token as u32 || response.revision != revision {
                        return Err("Batch start response differs from its request.".into());
                    }
                    let began = match response.result {
                        AnalysisResult::CombinatorialBatchStarted(value) => value,
                        AnalysisResult::Error(error) => return Err(error),
                        _ => return Err("Unexpected batch-start response.".into()),
                    };
                    if began.context_revision != revision
                        || began.next_index != 0
                        || began.total_rows != matrix.rows.len()
                        || began.finished
                    {
                        return Err(
                            "Batch start context/count differs from the selected matrix.".into(),
                        );
                    }
                    for index in 0..matrix.rows.len() {
                        guard()?;
                        status.set(format!(
                            "Evaluating construct {} of {}…",
                            index + 1,
                            matrix.rows.len()
                        ));
                        let id = (token as u32).wrapping_add(index as u32 + 1);
                        let response = w
                            .ask(Request {
                                id,
                                revision: revision.clone(),
                                job: AnalysisJob::AdvanceCombinatorialBatch {
                                    context: revision.clone(),
                                    index,
                                },
                            })
                            .await?;
                        guard()?;
                        if response.id != id || response.revision != revision {
                            return Err("Batch row response differs from its request.".into());
                        }
                        let progress = match response.result {
                            AnalysisResult::CombinatorialBatchAdvanced(value) => value,
                            AnalysisResult::Error(error) => return Err(error),
                            _ => return Err("Unexpected batch-row response.".into()),
                        };
                        if progress.context_revision != revision
                            || progress.next_index != index + 1
                            || progress.total_rows != matrix.rows.len()
                            || progress.finished != (index + 1 == matrix.rows.len())
                        {
                            return Err(
                                "Batch progress differs from the requested construct order.".into(),
                            );
                        }
                        let row = progress.row.ok_or("Batch row evidence is missing.")?;
                        if row.index != index || row.construct_id != row_id(&matrix.rows[index]) {
                            return Err(
                                "Row evidence belongs to another selected construct.".into()
                            );
                        }
                        rows.update(|all| all.push(row));
                    }
                    guard()?;
                    Ok::<(), String>(())
                }
                .await;
                if alive.load(Ordering::Relaxed) && nonce.get() == token {
                    evaluating.set(false);
                    match result {
                        Ok(()) => {
                            finished.set(true);
                            status.set("All selected rows have a result. Inspect evaluated, failed and uncomputed rows before saving evidence. No construct outputs were created.".into());
                        }
                        Err(error) => {
                            w.close();
                            worker.set_value(None);
                            scope.set(None);
                            status.set(error);
                        }
                    }
                } else {
                    w.close();
                }
            });
        })
    };
    let run_current = move |_| match draft.with_untracked(assembly_batch_projection::project) {
        Ok(projected) => start.run((projected.request, Origin::Current(projected.revision))),
        Err(error) => status.set(error),
    };
    let capture = {
        let alive = alive.clone();
        move |_| {
            if !finished.get_untracked()
                || evaluating.get_untracked()
                || workspace_busy.get_untracked()
            {
                return;
            }
            let prepared = (|| -> Result<_, String> {
                let binding = scope
                    .get_untracked()
                    .ok_or("Rerun the matrix before saving its evidence.")?;
                current.run(binding.clone())?;
                let w = worker
                    .get_value()
                    .ok_or("The private batch worker is unavailable; rerun.")?;
                let database = db
                    .get_value()
                    .ok_or("Browser recovery storage is unavailable; no record was saved.")?;
                let frame = snapshot.run(())?;
                let fingerprint =
                    run_lineage::workspace_fingerprint(&frame).map_err(|e| e.to_string())?;
                let revision =
                    batch::capture_revision(&frame, &binding.context).map_err(|e| e.to_string())?;
                Ok((binding, w, database, frame, fingerprint, revision))
            })();
            let (binding, w, database, frame, fingerprint, revision) = match prepared {
                Ok(v) => v,
                Err(e) => {
                    status.set(e);
                    return;
                }
            };
            let nonce = generation.get_value();
            let token = nonce.get().wrapping_add(1);
            nonce.set(token);
            let pending = transaction.get_value();
            let lock = owns_lock.get_value();
            if let Err(error) = lock.acquire(token) {
                status.set(error);
                return;
            }
            workspace_busy.set(true);
            capturing.set(true);
            save_generation.update(|n| *n = n.wrapping_add(1));
            status.set(
                "Saving immutable batch evidence and browser recovery in one transaction…".into(),
            );
            let alive = alive.clone();
            spawn_local(async move {
                let guard = || -> Result<(), String> {
                    if !alive.load(Ordering::Relaxed) || nonce.get() != token {
                        return Err("Batch capture cancelled.".into());
                    }
                    current.run(binding.clone())?;
                    if run_lineage::workspace_fingerprint(&snapshot.run(())?)
                        .map_err(|e| e.to_string())?
                        != fingerprint
                    {
                        return Err(
                            "The workspace changed while saving; no batch record was installed."
                                .into(),
                        );
                    }
                    Ok(())
                };
                let result = async {
                    let response = w
                        .ask(Request {
                            id: u32::MAX - 550,
                            revision: revision.clone(),
                            job: AnalysisJob::CaptureCombinatorialBatch {
                                workspace: frame.clone(),
                                current_context: binding.context.clone(),
                            },
                        })
                        .await?;
                    guard()?;
                    if response.id != u32::MAX - 550 || response.revision != revision {
                        return Err("Batch capture response differs from its request.".into());
                    }
                    let commit = match response.result {
                        AnalysisResult::CombinatorialBatchCaptured(commit) => commit,
                        AnalysisResult::Error(error) => return Err(error),
                        _ => return Err("Unexpected batch capture response.".into()),
                    };
                    commit
                        .validate_binding(&revision)
                        .map_err(|e| e.to_string())?;
                    if commit.context_revision != binding.context
                        || commit.source_workspace_fingerprint != fingerprint
                    {
                        return Err(
                            "Batch receipt differs from the current workspace/session.".into()
                        );
                    }
                    crate::workspace_commit::persist(
                        &database,
                        commit.workspace,
                        frame,
                        validate,
                        w.clone(),
                        guard,
                        pending,
                    )
                    .await
                }
                .await;
                match result {
                    Ok(frame) => {
                        loaded.set(true);
                        message.set("Batch evidence saved as a historical native record. Construct release still requires fresh reconstruction.".into());
                        on_restore.run(frame);
                    }
                    Err(error) => {
                        if alive.load(Ordering::Relaxed) && nonce.get() == token {
                            status.set(error);
                        }
                    }
                }
                if lock.release(token) {
                    workspace_busy.set(false);
                }
                if alive.load(Ordering::Relaxed) && nonce.get() == token {
                    capturing.set(false);
                }
            });
        }
    };
    let load_history = move |_| match snapshot
        .run(())
        .and_then(|frame| batch::list(&frame).map_err(|e| e.to_string()))
    {
        Ok(all) => {
            history.set(all);
            status.set("Saved records are historical; opening one does not recompute it.".into());
        }
        Err(error) => status.set(error),
    };
    let open_record = Callback::new(move |hash: String| {
        let result = snapshot
            .run(())
            .and_then(|frame| batch::load(&frame, &hash).map_err(|e| e.to_string()));
        match result {
            Ok(record) => {
                cancel.run(());
                rows.set(record.rows.clone());
                total.set(record.rows.len());
                historical.set(true);
                historical_record.set(Some(Arc::new(record)));
                status.set("Historical batch opened without recomputation. Replay explicitly to obtain current row assessments.".into());
            }
            Err(error) => status.set(error),
        }
    });
    let replay = Callback::new(move |hash: String| {
        match snapshot
            .run(())
            .and_then(|frame| batch::load(&frame, &hash).map_err(|e| e.to_string()))
        {
            Ok(record) => start.run((record.request, Origin::Saved(hash))),
            Err(error) => status.set(error),
        }
    });
    let edit_construct = Callback::new(move |id: String| {
        let result = (|| -> Result<(), String> {
            if workspace_busy.get_untracked() {
                return Err("Wait for the workspace operation to finish.".into());
            }
            let binding = scope
                .get_untracked()
                .ok_or("Rerun the current matrix before editing from its results.")?;
            if !matches!(&binding.origin, Origin::Current(_)) {
                return Err("Historical batch results retain their original declaration. Use the current matrix to edit constructs.".into());
            }
            current.run(binding.clone())?;
            let before = draft.get_untracked();
            let next = crate::assembly_alternatives::open_row(&before, &id)?;
            current.run(binding)?;
            if draft.get_untracked() != before {
                return Err("The current assembly changed; choose the result again.".into());
            }
            draft.set(next);
            Ok(())
        })();
        status.set(match result {
            Ok(()) => "Opened the construct in the existing dock. Review its intended reference, then edit its preparation and rerun explicitly.".into(),
            Err(error) => error,
        });
    });
    let editable = Memo::new(move |_| {
        !workspace_busy.get()
            && scope.with(|value| {
                value
                    .as_ref()
                    .is_some_and(|binding| matches!(binding.origin, Origin::Current(_)))
            })
    });
    view!{<details class="assembly-batch"><summary>"Evaluate selected constructs"</summary>
        <p class="helper">"Evaluate the selected matrix using each construct's declared preparations and independent intended reference. This checks prepared-fragment assembly, not a kinetic reaction mixture. Invalid editable rows remain explicitly uncomputed; choosing a row does not imply experimental success."</p>
        <button class="button secondary" disabled=move||evaluating.get()||workspace_busy.get()||draft.with(|d|d.alternatives.is_none()) on:click=run_current>"Evaluate selected matrix"</button>
        <button class="text-tool" disabled=move||!evaluating.get()&&!capturing.get() on:click=move|_|cancel.run(())>"Cancel batch operation"</button>
        <p>{move||format!("{} of {} row results · {}",rows.with(Vec::len),total.get(),if historical.get(){"historical; not replayed"}else if finished.get(){"completed computation on pinned inputs"}else{"transient; not saved"})}</p>
        <BatchRows rows editable on_edit=edit_construct/>
        <button class="button secondary" disabled=move||!finished.get()||scope.get().is_none()||workspace_busy.get() on:click=capture>"Save batch evidence"</button>
        <p class="helper" role="status">{move||status.get()}</p>
        <button class="text-tool" disabled=move||workspace_busy.get() on:click=load_history>"Load saved batch records"</button>
        <BatchHistory history open_record replay/>
        {move||historical_record.get().map(|record|view!{<button class="button secondary" on:click=move|_|match serde_json::to_vec(record.as_ref()){Ok(bytes)=>if let Err(error)=crate::app::download("combinatorial-batch.json","application/json",&bytes){status.set(error)},Err(error)=>status.set(error.to_string())}>"Download historical batch evidence JSON"</button>})}
    </details>}.into_any()
}

#[component]
fn BatchRows(
    rows: RwSignal<Vec<RowEvidence>>,
    editable: Memo<bool>,
    on_edit: Callback<String>,
) -> impl IntoView {
    let page = RwSignal::new(0usize);
    let selected = RwSignal::new(None::<usize>);
    let note_page = RwSignal::new(0usize);
    Effect::new(move || {
        let count = rows.with(Vec::len);
        if selected.get_untracked().is_some_and(|i| i >= count) {
            selected.set(None);
        }
        if count == 0 {
            page.set(0);
        }
    });
    view!{<section class="batch-row-results">
        {move||rows.with(|all|{
            let pages=all.len().div_ceil(20).max(1);let current=page.get().min(pages-1);
            let rendered=all.iter().enumerate().skip(current*20).take(20).map(|(index,row)|{
                let id=row.construct_id.clone();
                view!{<article><button class="text-tool" on:click=move|_|{selected.set(Some(index));note_page.set(0);}>{format!("{} · {} · {:?} · {}",index+1,row.name,row.state,row.status.map(|s|format!("{s:?}")).unwrap_or_else(||"not computed".into()))}</button><p class="helper">{format!("Construct {} · completed steps {}/{} · preparation complete {}",row.construct_id,row.completed_steps,row.planned_steps,row.preparation_complete)}</p>{row.error.clone().map(|error|view!{<p role="alert">{error}</p>})}<button class="text-tool" disabled=move||!editable.get() on:click=move|_|on_edit.run(id.clone())>"Edit construct in assembly dock"</button></article>}
            }).collect_view();
            view!{<div>{rendered}<p>{format!("Result page {} of {pages}",current+1)}</p><button class="text-tool" disabled=current==0 on:click=move|_|page.set(current.saturating_sub(1))>"Previous results"</button><button class="text-tool" disabled={current+1>=pages} on:click=move|_|page.set(current+1)>"Next results"</button></div>}
        })}
        {move||selected.get().and_then(|index|rows.with(|all|all.get(index).cloned())).map(|row|{
            let count=row.findings.len()+row.assumptions.len();let pages=count.div_ceil(50).max(1);let current=note_page.get().min(pages-1);
            let notes=row.findings.iter().map(|s|format!("Finding: {s}")).chain(row.assumptions.iter().map(|s|format!("Assumption: {s}"))).skip(current*50).take(50).map(|s|view!{<li>{s}</li>}).collect_view();
            view!{<section><h4>{row.name.clone()}</h4><p>{format!("Construct {} · combination {} · declaration {}",row.construct_id,row.combination_id.as_deref().unwrap_or("unavailable"),row.declaration_revision.as_deref().unwrap_or("not computed"))}</p><p>{format!("Plan {} · intended target {} · final reconstructed product {}",row.plan_revision.as_deref().unwrap_or("unavailable"),row.target_revision.as_deref().unwrap_or("unavailable"),row.final_product_revision.as_deref().unwrap_or("unavailable"))}</p>
                <details><summary>"Pinned source events"</summary>{row.input_bindings.iter().map(|binding|view!{<p>{format!("{}: document {} · event {} · revision {}",binding.input_id,binding.event.document_id,binding.event.event_id,binding.event.revision)}</p>}).collect_view()}</details>
                <p class="helper">{format!("Exhaustive alternative-product search: {:?}. A consistent sequence model does not establish amplification, sample identity or laboratory success.",row.exhaustive_alternative_search)}</p><ul>{notes}</ul><p>{format!("Evidence note page {} of {pages}",current+1)}</p><button class="text-tool" disabled=current==0 on:click=move|_|note_page.set(current.saturating_sub(1))>"Previous notes"</button><button class="text-tool" disabled={current+1>=pages} on:click=move|_|note_page.set(current+1)>"Next notes"</button>
            </section>}
        })}
    </section>}.into_any()
}

#[component]
fn BatchHistory(
    history: RwSignal<Vec<batch::Summary>>,
    open_record: Callback<String>,
    replay: Callback<String>,
) -> impl IntoView {
    let page = RwSignal::new(0usize);
    view!{<section>{move||history.with(|all|{
        let pages=all.len().div_ceil(10).max(1);let current=page.get().min(pages-1);
        let rows=all.iter().skip(current*10).take(10).map(|summary|{
            let open=summary.attachment_hash.clone();let rerun=open.clone();
            view!{<article><p>{format!("{} · {} rows · {} evaluated · {} errors · recorded {}",summary.design_name,summary.total_rows,summary.evaluated_rows,summary.error_rows,summary.recorded_engine)}</p><button class="text-tool" on:click=move|_|open_record.run(open.clone())>"Open historical batch"</button><button class="text-tool" on:click=move|_|replay.run(rerun.clone())>"Replay pinned matrix"</button></article>}
        }).collect_view();
        view!{<div>{rows}<p>{format!("Saved batch page {} of {pages}",current+1)}</p><button class="text-tool" disabled=current==0 on:click=move|_|page.set(current.saturating_sub(1))>"Previous saved batches"</button><button class="text-tool" disabled={current+1>=pages} on:click=move|_|page.set(current+1)>"Next saved batches"</button></div>}
    })}</section>}.into_any()
}
