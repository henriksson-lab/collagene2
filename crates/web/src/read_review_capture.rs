//! Durable read-review decisions and optional revised reference; no restored verdict.
use crate::worker::CachedWorker;
use collagene_core::{
    protocol::{AnalysisJob, AnalysisResult, Request},
    read_review, run_lineage,
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
fn fingerprint(w: &Workspace) -> Result<String, String> {
    run_lineage::workspace_fingerprint(w).map_err(|e| e.to_string())
}
#[component]
pub fn ReadReviewCapturePanel(
    request: collagene_core::read_review::Request,
    snapshot: Callback<(), Result<Workspace, String>>,
    validate: Callback<Workspace, Result<(), String>>,
    on_restore: Callback<read_review::Commit>,
    db: StoredValue<Option<web_sys::IdbDatabase>, LocalStorage>,
    busy: RwSignal<bool>,
    loaded: RwSignal<bool>,
    save_generation: RwSignal<u64>,
    epoch: RwSignal<u64>,
    message: RwSignal<String>,
) -> impl IntoView {
    let captured_epoch = epoch.get_untracked();
    let request = StoredValue::new(request);
    let status = RwSignal::new(String::new());
    let completed = RwSignal::new(false);
    let alive = Arc::new(AtomicBool::new(true));
    let worker = StoredValue::new_local(None::<Rc<CachedWorker>>);
    let transaction =
        StoredValue::new_local(Rc::new(RefCell::new(None::<web_sys::IdbTransaction>)));
    let generation = StoredValue::new_local(Rc::new(Cell::new(0u64)));
    let owns_lock = StoredValue::new_local(Rc::new(Cell::new(false)));
    let cancel = Callback::new(move |()| {
        generation.with_value(|g| g.set(g.get().wrapping_add(1)));
        if let Some(w) = worker.get_value() {
            w.close();
        }
        if let Some(tx) = transaction.with_value(|t| t.borrow().clone()) {
            let _ = tx.abort();
        }
        status.set(
            "Review saving cancelled. No partial review or revised reference was installed.".into(),
        );
    });
    let cleanup = alive.clone();
    on_cleanup(move || {
        cleanup.store(false, Ordering::Relaxed);
        generation.with_value(|g| g.set(g.get().wrapping_add(1)));
        if let Some(w) = worker.get_value() {
            w.close();
        }
        if let Some(tx) = transaction.with_value(|t| t.borrow().clone()) {
            let _ = tx.abort();
        } else if owns_lock.with_value(|o| o.replace(false)) {
            busy.set(false);
        }
    });
    Effect::new(move |_| {
        if epoch.get() != captured_epoch && owns_lock.with_value(|o| o.get()) {
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
                        "Workspace changed; preview the review decisions again before saving."
                            .into(),
                    );
                }
                let database = db
                    .get_value()
                    .ok_or("Recovery storage is unavailable. No review was saved.")?;
                let frame = snapshot.run(())?;
                let original = fingerprint(&frame)?;
                let r = request.get_value();
                let revision =
                    read_review::commit_revision(&frame, &r).map_err(|e| e.to_string())?;
                let w = Rc::new(CachedWorker::new()?);
                Ok((database, frame, original, r, revision, w))
            })();
            let (database, frame, original, r, revision, w) = match prepared {
                Ok(p) => p,
                Err(e) => {
                    status.set(e);
                    return;
                }
            };
            worker.set_value(Some(w.clone()));
            let generation = generation.get_value();
            let owns_lock = owns_lock.get_value();
            let transaction = transaction.get_value();
            generation.set(generation.get().wrapping_add(1));
            let token = generation.get();
            owns_lock.set(true);
            busy.set(true);
            save_generation.update(|g| *g = g.wrapping_add(1));
            status.set(
                "Replaying source evidence and reconstructing the explicit review decisions…"
                    .into(),
            );
            let alive = alive.clone();
            spawn_local(async move {
                let result = async {
                    let response = w
                        .ask(Request {
                            id: u32::MAX - 240,
                            revision: revision.clone(),
                            job: AnalysisJob::CommitReadReview {
                                workspace: frame.clone(),
                                request: r.clone(),
                            },
                        })
                        .await?;
                    if response.id != u32::MAX - 240 || response.revision != revision {
                        return Err("Read-review response belongs to another request.".into());
                    }
                    let commit = match response.result {
                        AnalysisResult::ReadReviewCommitted(c) => c,
                        AnalysisResult::Error(e) => return Err(e),
                        _ => return Err("Unexpected read-review commit result.".into()),
                    };
                    commit
                        .validate_binding(&revision)
                        .map_err(|e| e.to_string())?;
                    let guard = || -> Result<(), String> {
                        if !alive.load(Ordering::Relaxed)
                            || epoch.get_untracked() != captured_epoch
                            || generation.get() != token
                        {
                            return Err("Capture cancelled or context changed.".into());
                        }
                        if fingerprint(&snapshot.run(())?)? != original {
                            return Err(
                                "Workspace changed while reconstructing the review. No review was installed."
                                    .into(),
                            );
                        }
                        Ok(())
                    };
                    guard()?;
                    validate.run(commit.workspace.clone())?;
                    let current =
                        crate::app::encode_workspace(commit.workspace.clone(), u32::MAX - 241)
                            .await?;
                    let previous = crate::app::encode_workspace(frame, u32::MAX - 242).await?;
                    guard()?;
                    crate::storage::replace_workspace_cancellable(
                        &database,
                        &current,
                        &previous,
                        transaction.clone(),
                    )
                    .await?;
                    Ok::<_, String>(commit)
                }
                .await;
                match result {
                    Ok(commit) => {
                        // Durable completion is the commit point; restore even if the modal closed immediately afterwards.
                        loaded.set(true);
                        message.set("Read review saved with exact source evidence. Any revised reference requires a new read check; review decisions do not establish experimental identity or chemical state.".into());
                        if alive.load(Ordering::Relaxed) {
                            completed.set(true);
                            status.set("Read review saved.".into());
                        }
                        on_restore.run(commit);
                    }
                    Err(e) => {
                        if alive.load(Ordering::Relaxed) {
                            status.set(format!("Read review was not saved. {e}"));
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
    view! {<section class="read-review-capture"><h2>"Save read review"</h2><p>"Reconstruct the decisions from fresh source evidence and save them in one recovery transaction. Original read/reference heads remain unchanged. A revised reference is a separate derived design and must be checked again."</p><p>{request.with_value(|r|format!("{} decisions · {} · source check {}",r.decisions.len(),r.name,r.attachment_hash))}</p><p role="status">{move||status.get()}</p><button class="button primary" disabled=move||busy.get()||completed.get()||epoch.get()!=captured_epoch on:click=move|_|run.with_value(|r|r())>"Save read review"</button><button class="button secondary" disabled=move||!busy.get()||!owns_lock.with_value(|o|o.get()) on:click=move|_|cancel.run(())>"Cancel saving review"</button></section>}
}
