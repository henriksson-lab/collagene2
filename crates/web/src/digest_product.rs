//! Fresh digest derivation and atomic native installation, separate from panel rendering.
use crate::worker::CachedWorker;
use collagene_core::{
    catalog_digest::{self, CatalogSelection},
    protocol::{AnalysisJob, AnalysisResult, Request},
    restriction_history::SourceBinding,
    restriction_inspection::InspectionOptions,
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
pub struct DigestOpenRequest {
    pub source: SourceBinding,
    pub selection: CatalogSelection,
    pub options: InspectionOptions,
    pub fragment_index: usize,
    pub epoch: u64,
}
#[component]
pub fn DigestOpenPanel(
    request: DigestOpenRequest,
    snapshot: Callback<(), Result<Workspace, String>>,
    validate: Callback<Workspace, Result<(), String>>,
    on_restore: Callback<Workspace>,
    db: StoredValue<Option<web_sys::IdbDatabase>, LocalStorage>,
    busy: RwSignal<bool>,
    loaded: RwSignal<bool>,
    save_generation: RwSignal<u64>,
    epoch: RwSignal<u64>,
    message: RwSignal<String>,
) -> impl IntoView {
    let request = StoredValue::new(request);
    let status = RwSignal::new("Preparing the selected fragment…".to_owned());
    let failed = RwSignal::new(false);
    let started = RwSignal::new(false);
    let alive = Arc::new(AtomicBool::new(true));
    let worker = StoredValue::new_local(None::<Rc<CachedWorker>>);
    let transaction = Rc::new(RefCell::new(None::<web_sys::IdbTransaction>));
    let owns_lock = Rc::new(Cell::new(false));
    let cleanup = alive.clone();
    let pending = StoredValue::new_local(transaction.clone());
    let lock = StoredValue::new_local(owns_lock.clone());
    on_cleanup(move || {
        cleanup.store(false, Ordering::Relaxed);
        if let Some(w) = worker.get_value() {
            w.close();
        }
        let pending = pending.get_value();
        let lock = lock.get_value();
        let active_transaction = pending.borrow().clone();
        if let Some(tx) = active_transaction {
            let _ = tx.abort();
        } else if lock.replace(false) {
            busy.set(false);
        }
    });
    let watch_alive = alive.clone();
    let watch_transaction = transaction.clone();
    let watch_lock = owns_lock.clone();
    Effect::new(move |_| {
        let current = epoch.get();
        if watch_lock.get() && current != request.with_value(|r| r.epoch) {
            watch_alive.store(false, Ordering::Relaxed);
            if let Some(w) = worker.get_value() {
                w.close();
            }
            if let Some(tx) = watch_transaction.borrow().as_ref() {
                let _ = tx.abort();
            }
            status.set("The workspace changed; fragment opening was cancelled.".into());
        }
    });
    let run = {
        let alive = alive.clone();
        let transaction = transaction.clone();
        let owns_lock = owns_lock.clone();
        move || {
            if busy.get_untracked() {
                return;
            }
            let Some(database) = db.get_value() else {
                status.set(
                    "Browser recovery storage is unavailable. The workspace was not changed."
                        .into(),
                );
                failed.set(true);
                return;
            };
            let input = request.get_value();
            if input.epoch != epoch.get_untracked() {
                status.set(
                    "The workspace changed. Run the digest again before opening a fragment.".into(),
                );
                failed.set(true);
                return;
            }
            busy.set(true);
            owns_lock.set(true);
            save_generation.update(|v| *v = v.wrapping_add(1));
            failed.set(false);
            let prepared = (|| -> Result<_, String> {
                let frame = snapshot.run(())?;
                let fingerprint =
                    run_lineage::workspace_fingerprint(&frame).map_err(|e| e.to_string())?;
                let catalog = catalog_digest::builtin_catalog().map_err(|e| e.to_string())?;
                let revision = catalog_digest::installation_revision(
                    &catalog,
                    &input.selection,
                    &frame,
                    &input.source,
                    &input.options,
                    input.fragment_index,
                )
                .map_err(|e| e.to_string())?;
                let w = Rc::new(CachedWorker::new()?);
                worker.set_value(Some(w.clone()));
                Ok((frame, fingerprint, revision, w, catalog))
            })();
            let (frame, fingerprint, revision, w, catalog) = match prepared {
                Ok(v) => v,
                Err(e) => {
                    owns_lock.set(false);
                    busy.set(false);
                    status.set(e);
                    failed.set(true);
                    return;
                }
            };
            status.set("Reconstructing the fragment and preserving its source history…".into());
            let alive = alive.clone();
            let transaction = transaction.clone();
            let owns_lock = owns_lock.clone();
            spawn_local(async move {
                let result = async {
                    let response = w
                        .ask(Request {
                            id: u32::MAX - 60,
                            revision: revision.clone(),
                            job: AnalysisJob::InstallCatalogDigestFragment {
                                workspace: frame.clone(),
                                source: input.source.clone(),
                                selection: input.selection.clone(),
                                options: input.options.clone(),
                                fragment_index: input.fragment_index,
                            },
                        })
                        .await?;
                    if response.id != u32::MAX - 60 || response.revision != revision {
                        return Err(
                            "Digest installation response belongs to another request.".into()
                        );
                    }
                    let installed = match response.result {
                        AnalysisResult::CatalogDigestInstalled(p) => p,
                        AnalysisResult::Error(e) => return Err(e),
                        _ => return Err("Unexpected digest installation result.".into()),
                    };
                    installed
                        .validate_bindings_for(
                            &catalog,
                            &input.selection,
                            &frame,
                            &input.source,
                            &input.options,
                            input.fragment_index,
                        )
                        .map_err(|e| e.to_string())?;
                    let guard = || -> Result<(), String> {
                        if !alive.load(Ordering::Relaxed) || epoch.get_untracked() != input.epoch {
                            return Err(
                                "The operation was cancelled or its workspace changed.".into()
                            );
                        }
                        let current = snapshot.run(())?;
                        if run_lineage::workspace_fingerprint(&current)
                            .map_err(|e| e.to_string())?
                            != fingerprint
                        {
                            return Err(
                                "The workspace changed while the fragment was being prepared."
                                    .into(),
                            );
                        }
                        Ok(())
                    };
                    guard()?;
                    let next = crate::app::decode_workspace(
                        installed.installed.archive.clone(),
                        u32::MAX - 61,
                    )
                    .await?;
                    installed
                        .validate_archive_bindings_for(
                            &catalog,
                            &input.selection,
                            &frame,
                            &input.source,
                            &input.options,
                            input.fragment_index,
                            &next,
                        )
                        .map_err(|e| e.to_string())?;
                    validate.run(next.clone())?;
                    let previous = crate::app::encode_workspace(frame, u32::MAX - 62).await?;
                    guard()?;
                    crate::storage::replace_workspace_cancellable(
                        &database,
                        &installed.installed.archive,
                        &previous,
                        transaction.clone(),
                    )
                    .await?;
                    Ok::<_, String>(next)
                }
                .await;
                match result {
                    Ok(next) => {
                        loaded.set(true);
                        message.set("Opened the digest fragment with its exact source-event parent. Cleavage conditions remain assumptions; saved end summaries are producer evidence.".into());
                        on_restore.run(next);
                    }
                    Err(e) => {
                        if alive.load(Ordering::Relaxed) {
                            status.set(format!("The fragment was not opened. {e}"));
                            failed.set(true);
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
    Effect::new(move |_| {
        if !started.get_untracked() {
            started.set(true);
            run.with_value(|start| start());
        }
    });
    view! {<section class="digest-open-panel"><p>"Opening a freshly reconstructed fragment. Its source sequence stays in the workspace."</p><p role="status">{move||status.get()}</p><Show when=move||failed.get()><button class="button primary" disabled=move||busy.get() on:click=move|_|run.with_value(|start|start())>"Retry fragment opening"</button></Show></section>}
}
