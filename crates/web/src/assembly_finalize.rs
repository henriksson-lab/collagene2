//! Explicit output presentation choices followed by fresh reconstruction and one
//! atomic recovery commit. Imported receipts never authorize installation.
use crate::worker::CachedWorker;
use collagene_core::{
    cloning_outputs::{self, Options},
    protocol::{AnalysisJob, AnalysisResult, Request},
    run_lineage::{self, EventBinding, InputEventBinding},
    verification::DesignPlan,
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
struct Binding {
    fingerprint: String,
    target: EventBinding,
    inputs: Vec<InputEventBinding>,
}
fn capture_binding(
    workspace: &Workspace,
    plan: &DesignPlan,
    inputs: Vec<InputEventBinding>,
) -> Result<Binding, String> {
    let fingerprint = run_lineage::workspace_fingerprint(workspace).map_err(|e| e.to_string())?;
    let document = workspace
        .documents
        .iter()
        .find(|d| d.id == workspace.active_document)
        .ok_or("Choose the independent intended-target document before creating outputs.")?;
    if document.head != plan.target_revision {
        return Err(
            "The active target differs from this plan. Review the intended construct again.".into(),
        );
    }
    let target = EventBinding {
        document_id: document.id.clone(),
        event_id: document
            .head_event
            .clone()
            .ok_or("The target has no retained history event.")?,
        revision: document.head.clone(),
    };
    let exists = |binding: &EventBinding| {
        workspace.history.iter().any(|event| {
            event.id == binding.event_id
                && event.document_id == binding.document_id
                && event.revision == binding.revision
        }) && workspace.revisions.contains_key(&binding.revision)
    };
    if !exists(&target)
        || inputs.len() != plan.inputs.len()
        || plan.inputs.iter().zip(&inputs).any(|(input, binding)| {
            input.id != binding.input_id
                || input.revision != binding.event.revision
                || !exists(&binding.event)
        })
    {
        return Err("Exact source-event preferences are missing or no longer match the ordered plan inputs. Review the assembly again; identical DNA does not identify a source document.".into());
    }
    Ok(Binding {
        fingerprint,
        target,
        inputs,
    })
}

#[component]
pub fn AssemblyFinalizePanel(
    plan: DesignPlan,
    preferred_bindings: Vec<InputEventBinding>,
    snapshot: Callback<(), Result<Workspace, String>>,
    recovery_snapshot: Callback<(), Result<Workspace, String>>,
    validate: Callback<Workspace, Result<(), String>>,
    on_restore: Callback<Workspace>,
    db: StoredValue<Option<web_sys::IdbDatabase>, LocalStorage>,
    busy: RwSignal<bool>,
    loaded: RwSignal<bool>,
    save_generation: RwSignal<u64>,
    epoch: RwSignal<u64>,
    message: RwSignal<String>,
) -> impl IntoView {
    let frame = snapshot.run(());
    let default_name = frame
        .as_ref()
        .ok()
        .and_then(|workspace| workspace.head(&workspace.active_document).ok())
        .map(|m| m.name.clone())
        .filter(|name| !name.trim().is_empty() && name.len() <= 256)
        .unwrap_or_else(|| "Assembly".into());
    let name = RwSignal::new(default_name);
    let initial = frame.and_then(|frame| capture_binding(&frame, &plan, preferred_bindings));
    let status = RwSignal::new(initial.as_ref().err().cloned().unwrap_or_default());
    let binding = StoredValue::new(initial.ok());
    let plan = StoredValue::new(plan);
    let captured_epoch = epoch.get_untracked();
    let show_prepared = RwSignal::new(false);
    let save_primers = RwSignal::new(false);
    let allow_incomplete = RwSignal::new(false);
    let completed = RwSignal::new(false);
    let alive = Arc::new(AtomicBool::new(true));
    let worker = StoredValue::new_local(None::<Rc<CachedWorker>>);
    let transaction =
        StoredValue::new_local(Rc::new(RefCell::new(None::<web_sys::IdbTransaction>)));
    let owns_lock = StoredValue::new_local(Rc::new(Cell::new(false)));
    let generation = StoredValue::new_local(Rc::new(Cell::new(0u64)));
    let cancel = Callback::new(move |()| {
        generation.with_value(|g| g.set(g.get().wrapping_add(1)));
        if let Some(worker) = worker.get_value() {
            worker.close();
        }
        let pending = transaction.with_value(|pending| pending.borrow().clone());
        if let Some(tx) = pending {
            let _ = tx.abort();
        }
        status.set("Creating outputs was cancelled. No partial outputs were installed.".into());
    });
    let cleanup_alive = alive.clone();
    on_cleanup(move || {
        cleanup_alive.store(false, Ordering::Relaxed);
        generation.with_value(|g| g.set(g.get().wrapping_add(1)));
        if let Some(worker) = worker.get_value() {
            worker.close();
        }
        let pending = transaction.with_value(|pending| pending.borrow().clone());
        if let Some(tx) = pending {
            let _ = tx.abort();
        } else if owns_lock.with_value(|lock| lock.replace(false)) {
            busy.set(false);
        }
    });
    Effect::new(move |_| {
        if epoch.get() != captured_epoch && owns_lock.with_value(|lock| lock.get()) {
            cancel.run(());
        }
    });
    let run = {
        let alive = alive.clone();
        move || {
            if busy.get_untracked() || completed.get_untracked() {
                return;
            }
            let prepared = (|| -> Result<_, String> {
                if epoch.get_untracked() != captured_epoch {
                    return Err(
                        "The workspace changed. Review the assembly again before creating outputs."
                            .into(),
                    );
                }
                let database = db
                    .get_value()
                    .ok_or("Browser recovery storage is unavailable; no outputs were created.")?;
                let binding = binding
                    .get_value()
                    .ok_or("Exact source events are unavailable; review the assembly again.")?;
                let name = name.get_untracked();
                if name.trim().is_empty() || name.len() > 256 {
                    return Err("Enter a construct name containing 1–256 bytes of text.".into());
                }
                let options = Options {
                    name,
                    show_prepared: show_prepared.get_untracked(),
                    save_primers: save_primers.get_untracked(),
                    allow_incomplete: allow_incomplete.get_untracked(),
                };
                options.validate().map_err(|e| e.to_string())?;
                let frame = snapshot.run(())?;
                if run_lineage::workspace_fingerprint(&frame).map_err(|e| e.to_string())?
                    != binding.fingerprint
                {
                    return Err(
                        "The workspace changed while choosing outputs. Review the assembly again."
                            .into(),
                    );
                }
                let previous = recovery_snapshot.run(())?;
                let plan = plan.get_value();
                let revision = cloning_outputs::finalize_revision(
                    &frame,
                    &plan,
                    binding.target.clone(),
                    binding.inputs.clone(),
                    &options,
                )
                .map_err(|e| e.to_string())?;
                let worker = Rc::new(CachedWorker::new()?);
                Ok((
                    database, binding, options, frame, previous, plan, revision, worker,
                ))
            })();
            let (database, binding, options, frame, previous, plan, revision, w) = match prepared {
                Ok(value) => value,
                Err(error) => {
                    status.set(error);
                    return;
                }
            };
            worker.set_value(Some(w.clone()));
            let generation = generation.get_value();
            let owns_lock = owns_lock.get_value();
            let transaction = transaction.get_value();
            generation.set(generation.get().wrapping_add(1));
            let request_generation = generation.get();
            owns_lock.set(true);
            busy.set(true);
            save_generation.update(|g| *g = g.wrapping_add(1));
            status.set("Reconstructing the plan and retaining every product with its exact source history…".into());
            let alive = alive.clone();
            spawn_local(async move {
                let result = async {
                    let response = w.ask(Request {
                        id: u32::MAX - 210,
                        revision: revision.clone(),
                        job: AnalysisJob::FinalizeCloningOutputs {workspace:frame.clone(),plan,target:binding.target.clone(),inputs:binding.inputs.clone(),options},
                    }).await?;
                    if response.id != u32::MAX - 210 || response.revision != revision {
                        return Err("Output response belongs to another request.".into());
                    }
                    let commit = match response.result {
                        AnalysisResult::CloningOutputsFinalized(commit) => commit,
                        AnalysisResult::Error(error) => return Err(error),
                        _ => return Err("Unexpected assembly-output result.".into()),
                    };
                    commit.validate_binding(&revision).map_err(|e|e.to_string())?;
                    if commit.receipt.source_workspace_fingerprint != binding.fingerprint
                        || commit.receipt.run.recorded_engine != collagene_core::verification::ENGINE_VERSION {
                        return Err("Output receipt does not match the source workspace and current reconstruction engine.".into());
                    }
                    let guard = || -> Result<(), String> {
                        if !alive.load(Ordering::Relaxed) || epoch.get_untracked() != captured_epoch
                            || generation.get() != request_generation {
                            return Err("The operation was cancelled or its workspace changed.".into());
                        }
                        if run_lineage::workspace_fingerprint(&snapshot.run(())?).map_err(|e|e.to_string())? != binding.fingerprint {
                            return Err("The workspace changed during reconstruction; no outputs were installed.".into());
                        }
                        Ok(())
                    };
                    guard()?;
                    validate.run(commit.workspace.clone())?;
                    let current = crate::app::encode_workspace(commit.workspace.clone(),u32::MAX - 211).await?;
                    let previous = crate::app::encode_workspace(previous,u32::MAX - 212).await?;
                    guard()?;
                    crate::storage::replace_workspace_cancellable(&database,&current,&previous,transaction.clone()).await?;
                    Ok::<_,String>(commit.workspace)
                }.await;
                match result {
                    Ok(workspace) => {
                        // The durable transaction is the commit point. Restore even
                        // if the modal was closed just after that transaction committed.
                        loaded.set(true);
                        message.set("Assembly outputs created with full reconstruction history. Document visibility is a presentation choice; experimental success is not established.".into());
                        if alive.load(Ordering::Relaxed) {
                            completed.set(true);
                            status.set("Assembly outputs created.".into());
                        }
                        on_restore.run(workspace);
                    }
                    Err(error) => {
                        if alive.load(Ordering::Relaxed) {
                            status.set(format!("Outputs were not created. {error}"));
                        }
                    }
                }
                if owns_lock.replace(false) {
                    busy.set(false);
                }
                w.close();
            });
        }
    };
    let run = StoredValue::new_local(run);
    view! {<section class="assembly-finalize-panel">
        <h2>"Create assembly outputs"</h2>
        <p>"Reconstruct the plan again and save the final construct with its exact source history. Every reconstructed preparation remains retained, even when hidden from Documents."</p>
        <label>"Construct name"<input aria-label="Assembly output name" maxlength="256" prop:value=move||name.get() disabled=move||busy.get()||completed.get() on:input=move|e|name.set(event_target_value(&e))/></label>
        <label><input type="checkbox" prop:checked=move||show_prepared.get() disabled=move||busy.get()||completed.get() on:change=move|e|show_prepared.set(event_target_checked(&e))/>"Show prepared fragments in Documents"</label>
        <label><input type="checkbox" prop:checked=move||save_primers.get() disabled=move||busy.get()||completed.get() on:change=move|e|save_primers.set(event_target_checked(&e))/>"Save primer collection"</label>
        <p class="helper">"Primer collections contain the declared oligos, not proof of synthesized or experimentally tested primers."</p>
        <details><summary>"Incomplete reconstruction"</summary>
            <p class="helper">"Normally, incomplete reconstructions are refused. Enable this only to retain a reconstruction with explicitly unresolved evidence. Invalid plans and missing final products remain refused; this does not turn an incomplete assessment into a verified outcome."</p>
            <label><input type="checkbox" prop:checked=move||allow_incomplete.get() disabled=move||busy.get()||completed.get() on:change=move|e|allow_incomplete.set(event_target_checked(&e))/>"Keep incomplete reconstruction"</label>
        </details>
        <p role="status">{move||status.get()}</p>
        <button class="button primary" disabled=move||busy.get()||completed.get()||binding.with_value(Option::is_none)||epoch.get()!=captured_epoch on:click=move|_|run.with_value(|run|run())>"Create outputs"</button>
        <button class="button secondary" disabled=move||!busy.get()||!owns_lock.with_value(|lock|lock.get()) on:click=move|_|cancel.run(())>"Cancel creating outputs"</button>
    </section>}
}
