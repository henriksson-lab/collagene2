//! Explicit bulk output creation through a fresh private finalizer. Staged rows
//! never enter the live/native workspace until one durable all-or-nothing commit.
use crate::{assembly_draft_model::Draft, commit_lease::CommitLease, worker::CachedWorker};
use collagene_core::{
    cloning_outputs,
    combinatorial_batch::{self as batch, RowRequest},
    combinatorial_release::{self as release, Inspection, ReleasedRow},
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

#[derive(Clone, PartialEq, Eq)]
struct RowForm {
    id: String,
    name: String,
    selected: bool,
    typed: bool,
    show_prepared: bool,
    save_primers: bool,
    allow_incomplete: bool,
}
impl RowForm {
    fn options(&self) -> cloning_outputs::Options {
        cloning_outputs::Options {
            name: self.name.clone(),
            show_prepared: self.show_prepared,
            save_primers: self.save_primers,
            allow_incomplete: self.allow_incomplete,
        }
    }
}
struct Bound {
    context: String,
    fingerprint: String,
    epoch: u64,
    nonce: u64,
    draft: Arc<Draft>,
    forms: Vec<RowForm>,
    request: release::Request,
}

#[component]
pub fn AssemblyReleaseView(
    batch_hash: String,
    draft: RwSignal<Draft>,
    epoch: RwSignal<u64>,
    snapshot: Callback<(), Result<Workspace, String>>,
    validate: Callback<Workspace, Result<(), String>>,
    on_restore: Callback<Workspace>,
    db: StoredValue<Option<crate::storage::Database>, LocalStorage>,
    workspace_busy: RwSignal<bool>,
    loaded: RwSignal<bool>,
    save_generation: RwSignal<u64>,
    message: RwSignal<String>,
    on_success: Callback<String>,
) -> impl IntoView {
    let batch_hash = StoredValue::new(batch_hash);
    let record = RwSignal::new(None::<Arc<batch::Record>>);
    let forms = RwSignal::new(Vec::<RowForm>::new());
    let scope = RwSignal::new(None::<Arc<Bound>>);
    let staged = RwSignal::new(Vec::<ReleasedRow>::new());
    let inspection = RwSignal::new(None::<Arc<Inspection>>);
    let working = RwSignal::new(false);
    let inspecting = RwSignal::new(false);
    let publishing = RwSignal::new(false);
    let finished = RwSignal::new(false);
    let status = RwSignal::new(String::new());
    let worker = StoredValue::new_local(None::<Rc<CachedWorker>>);
    let generation = StoredValue::new_local(Rc::new(Cell::new(0u64)));
    let transaction =
        StoredValue::new_local(Rc::new(RefCell::new(None::<crate::storage::Transaction>)));
    let lease = StoredValue::new_local(Rc::new(CommitLease::default()));
    let alive = Arc::new(AtomicBool::new(true));
    let cancel = Callback::new(move |()| {
        generation.with_value(|value| value.set(value.get().wrapping_add(1)));
        if let Some(w) = worker.get_value() {
            w.close();
        }
        worker.set_value(None);
        if let Some(tx) = transaction.with_value(|value| value.borrow().clone()) {
            let _ = tx.abort();
        } else if lease.with_value(|value| value.cancel_unstarted()) {
            workspace_busy.set(false);
        }
        scope.set(None);
        working.set(false);
        inspecting.set(false);
        publishing.set(false);
        finished.set(false);
        inspection.set(None);
        status.set("Private output creation cancelled. No staged prefix was published.".into());
    });
    let cleanup_alive = alive.clone();
    on_cleanup(move || {
        cleanup_alive.store(false, Ordering::Relaxed);
        generation.with_value(|value| value.set(value.get().wrapping_add(1)));
        if let Some(w) = worker.get_value() {
            w.close();
        }
        if let Some(tx) = transaction.with_value(|value| value.borrow().clone()) {
            let _ = tx.abort();
        } else if lease.with_value(|value| value.cancel_unstarted()) {
            workspace_busy.set(false);
        }
    });
    let current = Callback::new(move |bound: Arc<Bound>| -> Result<(), String> {
        if epoch.get_untracked() != bound.epoch
            || generation.with_value(|value| value.get()) != bound.nonce
            || !scope.with_untracked(|value| {
                value
                    .as_ref()
                    .is_some_and(|current| Arc::ptr_eq(current, &bound))
            })
            || forms.get_untracked() != bound.forms
            || !draft.with_untracked(|value| value == bound.draft.as_ref())
        {
            return Err("Output selection, options, assembly or workspace changed; start a new private reconstruction.".into());
        }
        Ok(())
    });
    let current_frame = Callback::new(move |bound: Arc<Bound>| -> Result<(), String> {
        current.run(bound.clone())?;
        if run_lineage::workspace_fingerprint(&snapshot.run(())?).map_err(|e| e.to_string())?
            != bound.fingerprint
        {
            return Err(
                "The source workspace changed; staged outputs cannot replace it. Start again."
                    .into(),
            );
        }
        Ok(())
    });
    Effect::new(move || {
        let _ = epoch.get();
        cancel.run(());
        match snapshot.run(()).and_then(|frame| {
            batch::load(&frame, &batch_hash.get_value()).map_err(|e| e.to_string())
        }) {
            Ok(value) => {
                let rows = value
                    .request
                    .rows
                    .iter()
                    .map(|row| {
                        let (id, name, typed) = match row {
                            RowRequest::Ready(row) => (row.id.clone(), row.name.clone(), true),
                            RowRequest::DraftError { id, name, .. } => {
                                (id.clone(), name.clone(), false)
                            }
                        };
                        RowForm {
                            id,
                            name,
                            selected: false,
                            typed,
                            show_prepared: false,
                            save_primers: false,
                            allow_incomplete: false,
                        }
                    })
                    .collect();
                record.set(Some(Arc::new(value)));
                forms.set(rows);
                staged.set(vec![]);
                status.set("Choose output rows and names. Recorded batch statuses are historical; every selected construct will be reconstructed and finalized afresh.".into());
            }
            Err(error) => {
                record.set(None);
                forms.set(vec![]);
                staged.set(vec![]);
                status.set(error);
            }
        }
    });
    Effect::new(move || {
        let now = epoch.get();
        let form = forms.get();
        let same = draft.with(|value| {
            scope.with_untracked(|bound| {
                bound
                    .as_ref()
                    .is_none_or(|bound| value == bound.draft.as_ref())
            })
        });
        if scope.with_untracked(|bound| {
            bound
                .as_ref()
                .is_some_and(|bound| bound.epoch != now || bound.forms != form || !same)
        }) {
            cancel.run(());
            status.set("Selection or editable context changed. Private outputs were discarded; reconstruct again.".into());
        }
    });
    let start = {
        let alive = alive.clone();
        move |_| {
            if workspace_busy.get_untracked()
                || working.get_untracked()
                || publishing.get_untracked()
            {
                return;
            }
            let prepared = (|| -> Result<_, String> {
                let form = forms.get_untracked();
                let rows=form.iter().filter(|row|row.selected).map(|row|{
                    if !row.typed{return Err("Uncomputed editable rows need a corrected design and a new batch record before output creation.".into());}
                    let options=row.options();options.validate().map_err(|e|format!("{}: {e}",row.id))?;
                    Ok(release::Selection{construct_id:row.id.clone(),options})
                }).collect::<Result<Vec<_>,String>>()?;
                let request = release::Request {
                    batch_hash: batch_hash.get_value(),
                    rows,
                };
                request.validate().map_err(|e| e.to_string())?;
                let frame = snapshot.run(())?;
                let fingerprint =
                    run_lineage::workspace_fingerprint(&frame).map_err(|e| e.to_string())?;
                let context = release::revision(&frame, &request).map_err(|e| e.to_string())?;
                let nonce = generation.get_value();
                let token = nonce.get().wrapping_add(1);
                let lock = lease.get_value();
                lock.acquire(token)?;
                nonce.set(token);
                let bound = Arc::new(Bound {
                    context,
                    fingerprint,
                    epoch: epoch.get_untracked(),
                    nonce: token,
                    draft: Arc::new(draft.get_untracked()),
                    forms: form,
                    request,
                });
                let w = match CachedWorker::new() {
                    Ok(w) => Rc::new(w),
                    Err(error) => {
                        lock.release(token);
                        return Err(error);
                    }
                };
                Ok((frame, bound, w, nonce, lock))
            })();
            let (frame, bound, w, nonce, lock) = match prepared {
                Ok(value) => value,
                Err(error) => {
                    status.set(error);
                    return;
                }
            };
            workspace_busy.set(true);
            save_generation.update(|value| *value = value.wrapping_add(1));
            scope.set(Some(bound.clone()));
            worker.set_value(Some(w.clone()));
            working.set(true);
            finished.set(false);
            inspection.set(None);
            staged.set(vec![]);
            status.set("Reconstructing selected outputs privately. Nothing will be published until every row succeeds and you explicitly publish.".into());
            let alive = alive.clone();
            spawn_local(async move {
                let token = bound.nonce;
                let guard = || -> Result<(), String> {
                    if !alive.load(Ordering::Relaxed) || nonce.get() != token {
                        return Err("Private reconstruction cancelled.".into());
                    }
                    current.run(bound.clone())
                };
                let result = async {
                    let response = w
                        .ask(Request {
                            id: token as u32,
                            revision: bound.context.clone(),
                            job: AnalysisJob::BeginCombinatorialRelease {
                                workspace: frame,
                                request: bound.request.clone(),
                            },
                        })
                        .await?;
                    guard()?;
                    if response.id != token as u32 || response.revision != bound.context {
                        return Err("Release start response differs from its request.".into());
                    }
                    let began = match response.result {
                        AnalysisResult::CombinatorialReleaseStarted(value) => value,
                        AnalysisResult::Error(error) => return Err(error),
                        _ => return Err("Unexpected release-start response.".into()),
                    };
                    if began.context_revision != bound.context
                        || began.next_index != 0
                        || began.total_rows != bound.request.rows.len()
                        || began.finished
                        || began.failed.is_some()
                    {
                        return Err("Private release start context/count is inconsistent.".into());
                    }
                    for index in 0..bound.request.rows.len() {
                        guard()?;
                        status.set(format!(
                            "Freshly finalizing output {} of {}…",
                            index + 1,
                            bound.request.rows.len()
                        ));
                        let id = (token as u32).wrapping_add(index as u32 + 1);
                        let response = w
                            .ask(Request {
                                id,
                                revision: bound.context.clone(),
                                job: AnalysisJob::AdvanceCombinatorialRelease {
                                    context: bound.context.clone(),
                                    index,
                                },
                            })
                            .await?;
                        guard()?;
                        if response.id != id || response.revision != bound.context {
                            return Err("Release row response differs from its request.".into());
                        }
                        let progress = match response.result {
                            AnalysisResult::CombinatorialReleaseAdvanced(value) => value,
                            AnalysisResult::Error(error) => return Err(error),
                            _ => return Err("Unexpected private release-row response.".into()),
                        };
                        if progress.context_revision != bound.context
                            || progress.total_rows != bound.request.rows.len()
                        {
                            return Err(
                                "Private release progress belongs to another selection.".into()
                            );
                        }
                        if let Some(error) = progress.failed {
                            return Err(format!(
                                "Output {} failed; the entire private stage was discarded. {error}",
                                index + 1
                            ));
                        }
                        if progress.next_index != index + 1
                            || progress.finished != (index + 1 == bound.request.rows.len())
                        {
                            return Err("Private release row order is inconsistent.".into());
                        }
                        let row = progress.row.ok_or("Private output receipt is missing.")?;
                        if row.ordinal != index
                            || row.construct_id != bound.request.rows[index].construct_id
                        {
                            return Err(
                                "Private output receipt belongs to another construct.".into()
                            );
                        }
                        row.validate().map_err(|e| e.to_string())?;
                        if serde_json::to_vec(&row.options).map_err(|e| e.to_string())?
                            != serde_json::to_vec(&bound.request.rows[index].options)
                                .map_err(|e| e.to_string())?
                        {
                            return Err(
                                "Private output options differ from the explicit selection.".into(),
                            );
                        }
                        staged.update(|rows| rows.push(row));
                    }
                    current_frame.run(bound.clone())?;
                    guard()?;
                    Ok::<(), String>(())
                }
                .await;
                if alive.load(Ordering::Relaxed) && nonce.get() == token {
                    working.set(false);
                    match result {
                        Ok(()) => {
                            finished.set(true);
                            status.set("All selected outputs were freshly finalized in private. Inspect their actual products and intended references, then publish all together or cancel.".into());
                        }
                        Err(error) => {
                            w.close();
                            worker.set_value(None);
                            scope.set(None);
                            inspection.set(None);
                            finished.set(false);
                            if lock.release(token) {
                                workspace_busy.set(false);
                            }
                            status.set(error);
                        }
                    }
                } else {
                    w.close();
                    if lock.release(token) {
                        workspace_busy.set(false);
                    }
                }
            });
        }
    };
    let inspect = {
        let alive = alive.clone();
        Callback::new(move |index: usize| {
            if working.get_untracked() || publishing.get_untracked() || inspecting.get_untracked() {
                return;
            }
            let prepared = (|| -> Result<_, String> {
                let bound = scope
                    .get_untracked()
                    .ok_or("Reconstruct selected outputs first.")?;
                current_frame.run(bound.clone())?;
                let w = worker
                    .get_value()
                    .ok_or("The private release worker is unavailable.")?;
                Ok((bound, w))
            })();
            let (bound, w) = match prepared {
                Ok(value) => value,
                Err(error) => {
                    status.set(error);
                    return;
                }
            };
            inspecting.set(true);
            inspection.set(None);
            status.set(
                "Loading the actual private product and its independent intended reference…".into(),
            );
            let nonce = generation.get_value();
            let lock = lease.get_value();
            let alive = alive.clone();
            spawn_local(async move {
                let token = bound.nonce;
                let result = async {
                    let response = w
                        .ask(Request {
                            id: u32::MAX - 575,
                            revision: bound.context.clone(),
                            job: AnalysisJob::InspectCombinatorialRelease {
                                context: bound.context.clone(),
                                index,
                            },
                        })
                        .await?;
                    if !alive.load(Ordering::Relaxed) || nonce.get() != token {
                        return Err("Inspection cancelled.".into());
                    }
                    current_frame.run(bound.clone())?;
                    if response.id != u32::MAX - 575 || response.revision != bound.context {
                        return Err("Inspection response differs from its private request.".into());
                    }
                    let value = match response.result {
                        AnalysisResult::CombinatorialReleaseInspected(value) => value,
                        AnalysisResult::Error(error) => return Err(error),
                        _ => return Err("Unexpected private inspection response.".into()),
                    };
                    value.validate().map_err(|e| e.to_string())?;
                    if value.context_revision != bound.context
                        || value.index != index
                        || bound
                            .request
                            .rows
                            .get(index)
                            .is_none_or(|row| row.construct_id != value.row.construct_id)
                    {
                        return Err(
                            "Private product inspection belongs to another selected row.".into(),
                        );
                    }
                    Ok::<_, String>(Arc::new(value))
                }
                .await;
                if alive.load(Ordering::Relaxed) && nonce.get() == token {
                    inspecting.set(false);
                    match result {
                        Ok(value) => {
                            inspection.set(Some(value));
                            status.set("Actual staged product and intended reference loaded. These read-only models are not yet published documents.".into());
                        }
                        Err(error) => {
                            w.close();
                            worker.set_value(None);
                            scope.set(None);
                            inspection.set(None);
                            finished.set(false);
                            if lock.release(token) {
                                workspace_busy.set(false);
                            }
                            status.set(error);
                        }
                    }
                } else {
                    w.close();
                    if lock.release(token) {
                        workspace_busy.set(false);
                    }
                }
            });
        })
    };
    let publish = {
        let alive = alive.clone();
        move |_| {
            if !finished.get_untracked()
                || working.get_untracked()
                || inspecting.get_untracked()
                || publishing.get_untracked()
            {
                return;
            }
            let prepared = (|| -> Result<_, String> {
                let bound = scope
                    .get_untracked()
                    .ok_or("Reconstruct the selected outputs first.")?;
                current_frame.run(bound.clone())?;
                let w = worker
                    .get_value()
                    .ok_or("The private release worker is unavailable.")?;
                let database = db
                    .get_value()
                    .ok_or("Browser recovery storage is unavailable; no outputs were published.")?;
                let frame = snapshot.run(())?;
                Ok((bound, w, database, frame))
            })();
            let (bound, w, database, frame) = match prepared {
                Ok(value) => value,
                Err(error) => {
                    status.set(error);
                    return;
                }
            };
            publishing.set(true);
            status.set("Publishing all selected outputs with current/previous browser recovery in one transaction…".into());
            let nonce = generation.get_value();
            let lock = lease.get_value();
            let pending = transaction.get_value();
            let alive = alive.clone();
            spawn_local(async move {
                let token = bound.nonce;
                let guard = || -> Result<(), String> {
                    if !alive.load(Ordering::Relaxed) || nonce.get() != token {
                        return Err("Output publication cancelled.".into());
                    }
                    current_frame.run(bound.clone())
                };
                let result = async {
                    let response = w
                        .ask(Request {
                            id: u32::MAX - 570,
                            revision: bound.context.clone(),
                            job: AnalysisJob::FinishCombinatorialRelease {
                                workspace: frame.clone(),
                                current_context: bound.context.clone(),
                            },
                        })
                        .await?;
                    guard()?;
                    if response.id != u32::MAX - 570 || response.revision != bound.context {
                        return Err("Private publication response differs from its request.".into());
                    }
                    let commit = match response.result {
                        AnalysisResult::CombinatorialReleaseFinished(value) => value,
                        AnalysisResult::Error(error) => return Err(error),
                        _ => return Err("Unexpected release publication response.".into()),
                    };
                    commit
                        .validate_binding(&bound.context, &bound.fingerprint)
                        .map_err(|e| e.to_string())?;
                    let manifest = commit.manifest_hash.clone();
                    let frame = crate::workspace_commit::persist(
                        &database,
                        commit.workspace,
                        frame,
                        validate,
                        w.clone(),
                        guard,
                        pending,
                    )
                    .await?;
                    Ok::<_, String>((frame, manifest))
                }
                .await;
                match result {
                    Ok((frame, manifest)) => {
                        loaded.set(true);
                        message.set("Selected constructs, retained history and requested oligos were published atomically. Laboratory success is not asserted.".into());
                        on_restore.run(frame);
                        on_success.run(manifest);
                    }
                    Err(error) => {
                        if alive.load(Ordering::Relaxed) && nonce.get() == token {
                            status.set(error);
                            finished.set(false);
                            scope.set(None);
                            inspection.set(None);
                            worker.set_value(None);
                        }
                    }
                }
                w.close();
                if lock.release(token) {
                    workspace_busy.set(false);
                }
                if alive.load(Ordering::Relaxed) && nonce.get() == token {
                    publishing.set(false);
                }
            });
        }
    };
    let inspection_current = Memo::new(move |_| {
        let now = epoch.get();
        let form = forms.get();
        inspection.with(|inspection| {
            scope.with(|bound| {
                inspection
                    .as_ref()
                    .zip(bound.as_ref())
                    .is_some_and(|(inspection, bound)| {
                        now == bound.epoch
                            && inspection.context_revision == bound.context
                            && form == bound.forms
                            && draft.with(|value| value == bound.draft.as_ref())
                    })
            })
        })
    });
    let rows_page = RwSignal::new(0usize);
    view!{<section class="assembly-release"><h3>"Choose and freshly reconstruct bulk outputs"</h3>
        <p class="helper">"Select constructs and output names explicitly. Recorded batch results are historical; each selected construct is reconstructed and passed through the current finalizer again. All outputs stay private until every row succeeds and you publish them together. A failed row discards the entire stage."</p>
        <p class="helper">"Source controls are locked while this private stage is open. Inspect the products here, publish all, or cancel to return to editing."</p>
        <div class="method-row-actions"><button class="text-tool" disabled=move||publishing.get() on:click=move|_|forms.update(|rows|for row in rows{row.selected=row.typed;})>"Select all typed rows"</button><button class="text-tool" disabled=move||publishing.get() on:click=move|_|forms.update(|rows|for row in rows{row.selected=false;})>"Clear selection"</button></div>
        {move||forms.with(|all|{
            let pages=all.len().div_ceil(20).max(1);let page=rows_page.get().min(pages-1);
            let rendered=all.iter().enumerate().skip(page*20).take(20).map(|(index,row)|{
                let recorded=record.with(|entry|{
                    entry.as_ref().and_then(|entry|entry.rows.get(index).map(|value|format!("Recorded {:?} · {}",value.state,value.status.map(|value|format!("{value:?}")).unwrap_or_else(||"not computed".into())))).unwrap_or_else(||"No historical row assessment".into())
                });
                view!{<fieldset><legend>{format!("{} · {}",index+1,row.id)}</legend>
                    <label><input type="checkbox" disabled=move||publishing.get()||!forms.with(|rows|rows.get(index).is_some_and(|row|row.typed)) prop:checked=move||forms.with(|rows|rows.get(index).is_some_and(|row|row.selected)) on:change=move|event|forms.update(|rows|if let Some(row)=rows.get_mut(index){row.selected=event_target_checked(&event);})/>"Create this construct"</label>
                    <p class="helper">{recorded}</p>{(!row.typed).then(||view!{<p class="helper">"This row has unfinished editable inputs. Correct its design and retain a new batch before creating outputs."</p>})}
                    <label>"Output name"<input maxlength="256" disabled=move||publishing.get() prop:value=move||forms.with(|rows|rows.get(index).map(|row|row.name.clone()).unwrap_or_default()) on:input=move|event|forms.update(|rows|if let Some(row)=rows.get_mut(index){row.name=event_target_value(&event);})/></label>
                    <label><input type="checkbox" disabled=move||publishing.get() prop:checked=move||forms.with(|rows|rows.get(index).is_some_and(|row|row.show_prepared)) on:change=move|event|forms.update(|rows|if let Some(row)=rows.get_mut(index){row.show_prepared=event_target_checked(&event);})/>"Show preparation documents"</label>
                    <label><input type="checkbox" disabled=move||publishing.get() prop:checked=move||forms.with(|rows|rows.get(index).is_some_and(|row|row.save_primers)) on:change=move|event|forms.update(|rows|if let Some(row)=rows.get_mut(index){row.save_primers=event_target_checked(&event);})/>"Create primer collection"</label>
                    <label><input type="checkbox" disabled=move||publishing.get() prop:checked=move||forms.with(|rows|rows.get(index).is_some_and(|row|row.allow_incomplete)) on:change=move|event|forms.update(|rows|if let Some(row)=rows.get_mut(index){row.allow_incomplete=event_target_checked(&event);})/>"Permit the native incomplete-assessment policy when a final construct is actually reconstructed"</label>
                </fieldset>}
            }).collect_view();
            view!{<section>{rendered}<p>{format!("Output option page {} of {pages}",page+1)}</p><button class="text-tool" disabled=page==0 on:click=move|_|rows_page.set(page.saturating_sub(1))>"Previous output rows"</button><button class="text-tool" disabled={page+1>=pages} on:click=move|_|rows_page.set(page+1)>"Next output rows"</button></section>}
        })}
        <button class="button primary" disabled=move||workspace_busy.get()||working.get()||publishing.get()||forms.with(|rows|!rows.iter().any(|row|row.selected)) on:click=start>"Freshly reconstruct selected outputs"</button>
        <button class="button secondary" disabled=move||scope.get().is_none() on:click=move|_|cancel.run(())>"Cancel and discard private outputs"</button>
        <p>{move||{let count=staged.with(Vec::len);if scope.get().is_some(){format!("{count} private outputs staged; none published yet.")}else{format!("{count} discarded staging receipts shown; no documents published.")}}}</p>
        {move||staged.with(|rows|rows.iter().enumerate().map(|(index,row)|view!{<article><p>{format!("{} · {} · staged construct {}",index+1,row.options.name,row.output.construct.document_id)}</p><button class="text-tool" disabled=move||working.get()||inspecting.get()||publishing.get()||scope.get().is_none() on:click=move|_|inspect.run(index)>"Inspect actual staged product and intended reference"</button></article>}).collect_view())}
        {move||inspection.get().map(|value|{
            let product=value.product.clone();let target=value.target.clone();let context=value.context_revision.clone();let index=value.index;
            view!{<details open><summary>{format!("Private product review · {}",value.row.options.name)}</summary><p class="helper">"Compare the reconstructed product with its independent intended reference. Both views are read-only; the staged documents remain unpublished until you confirm."</p>
                <h4>"Actual reconstructed product · not yet published"</h4><crate::assembly_product_view::ReferenceCanvas product identity=format!("release-product-{context}-{index}") current=inspection_current on_message=Callback::new(move|value|status.set(value))/>
                <h4>"Independent intended reference"</h4><crate::assembly_product_view::ReferenceCanvas product=target identity=format!("release-target-{context}-{index}") current=inspection_current on_message=Callback::new(move|value|status.set(value))/>
            </details>}
        })}
        <button class="button primary" disabled=move||!finished.get()||working.get()||inspecting.get()||publishing.get()||scope.get().is_none() on:click=publish>"Publish all selected outputs atomically"</button>
        <p class="helper" role="status">{move||status.get()}</p>
    </section>}.into_any()
}
