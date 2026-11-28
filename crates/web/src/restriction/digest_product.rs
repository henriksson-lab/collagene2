//! Fresh material derivation and atomic native installation, separate from panel rendering.
use crate::platform::worker::CachedWorker;
use collagene_core::{
    cloning::restriction::{catalog::CatalogSelection, inspection::InspectionOptions},
    documents::workspace::Workspace,
    material::{
        derivation::{self as material_derivation, Derivation, SourceOrigin},
        lineage::{self as run_lineage, EventBinding as SourceBinding},
    },
    protocol::{
        request::Request,
        schema::{AnalysisJob, AnalysisResult},
    },
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
    pub origin: Option<crate::history::run_origin_model::PinnedOrigin>,
    pub selection: CatalogSelection,
    pub options: InspectionOptions,
    pub fragment_index: usize,
    pub epoch: u64,
    pub guard: Callback<(), Result<(), String>>,
}
#[derive(Clone)]
pub struct MaterialOpenRequest {
    pub request: material_derivation::Request,
    pub epoch: u64,
    pub guard: Callback<(), Result<(), String>>,
}
#[component]
pub fn DigestOpenPanel(
    request: DigestOpenRequest,
    snapshot: Callback<(), Result<Workspace, String>>,
    validate: Callback<Workspace, Result<(), String>>,
    on_restore: Callback<Workspace>,
    db: StoredValue<Option<crate::platform::storage::Database>, LocalStorage>,
    busy: RwSignal<bool>,
    loaded: RwSignal<bool>,
    save_generation: RwSignal<u64>,
    epoch: RwSignal<u64>,
    message: RwSignal<String>,
) -> impl IntoView {
    let request = MaterialOpenRequest {
        request: material_derivation::Request {
            sources: vec![material_derivation::BoundSource {
                role: material_derivation::SourceRole::Template,
                event: request.source,
                origin: request
                    .origin
                    .as_ref()
                    .map_or(SourceOrigin::OriginalReference, |origin| {
                        SourceOrigin::Saved(origin.as_ref().clone())
                    }),
            }],
            derivation: Derivation::CatalogDigest {
                selection: request.selection,
                options: request.options,
                fragment_index: request.fragment_index,
            },
            name: format!("Digest fragment {}", request.fragment_index + 1),
            save_primers: false,
        },
        epoch: request.epoch,
        guard: request.guard,
    };
    view! {<MaterialOpenPanel request snapshot validate on_restore db busy loaded save_generation epoch message/>}
}
#[component]
pub fn MaterialOpenPanel(
    request: MaterialOpenRequest,
    snapshot: Callback<(), Result<Workspace, String>>,
    validate: Callback<Workspace, Result<(), String>>,
    on_restore: Callback<Workspace>,
    db: StoredValue<Option<crate::platform::storage::Database>, LocalStorage>,
    busy: RwSignal<bool>,
    loaded: RwSignal<bool>,
    save_generation: RwSignal<u64>,
    epoch: RwSignal<u64>,
    message: RwSignal<String>,
) -> impl IntoView {
    let request = StoredValue::new(request);
    let status = RwSignal::new("Preparing the selected product…".to_owned());
    let failed = RwSignal::new(false);
    let started = RwSignal::new(false);
    let alive = Arc::new(AtomicBool::new(true));
    let worker = StoredValue::new_local(None::<Rc<CachedWorker>>);
    let transaction = Rc::new(RefCell::new(None::<crate::platform::storage::Transaction>));
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
        let valid = request.with_value(|r| current == r.epoch && r.guard.run(()).is_ok());
        if watch_lock.get() && !valid {
            watch_alive.store(false, Ordering::Relaxed);
            if let Some(w) = worker.get_value() {
                w.close();
            }
            if let Some(tx) = watch_transaction.borrow().as_ref() {
                let _ = tx.abort();
            }
            status.set("The workspace changed; product opening was cancelled.".into());
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
            if input.epoch != epoch.get_untracked() || input.guard.run(()).is_err() {
                status.set(
                    "The workspace changed. Refresh the preview before opening a product.".into(),
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
                input.guard.run(())?;
                let revision = material_derivation::revision(&frame, &input.request)
                    .map_err(|e| e.to_string())?;
                let w = Rc::new(CachedWorker::new()?);
                worker.set_value(Some(w.clone()));
                Ok((frame, fingerprint, revision, w))
            })();
            let (frame, fingerprint, revision, w) = match prepared {
                Ok(v) => v,
                Err(e) => {
                    owns_lock.set(false);
                    busy.set(false);
                    status.set(e);
                    failed.set(true);
                    return;
                }
            };
            status.set("Reconstructing the product and preserving its source history…".into());
            let alive = alive.clone();
            let transaction = transaction.clone();
            let owns_lock = owns_lock.clone();
            spawn_local(async move {
                let result = async {
                    let response = w
                        .ask(Request {
                            id: u32::MAX - 60,
                            revision: revision.clone(),
                            job: AnalysisJob::InstallMaterialDerivation {
                                workspace: frame.clone(),
                                request: input.request.clone(),
                            },
                        })
                        .await?;
                    if response.id != u32::MAX - 60 || response.revision != revision {
                        return Err(
                            "Material installation response belongs to another request.".into()
                        );
                    }
                    let installed = match response.result {
                        AnalysisResult::MaterialDerivationInstalled(p) => p,
                        AnalysisResult::Error(e) => return Err(e),
                        _ => return Err("Unexpected material derivation result.".into()),
                    };
                    installed
                        .validate_binding(&revision)
                        .map_err(|e| e.to_string())?;
                    if installed.receipt.sources != input.request.source_bindings() {
                        return Err(
                            "Material derivation receipt source roles or exact events changed."
                                .into(),
                        );
                    }
                    let archive = installed.archive.clone();
                    let guard = || -> Result<(), String> {
                        if !alive.load(Ordering::Relaxed) || epoch.get_untracked() != input.epoch {
                            return Err(
                                "The operation was cancelled or its workspace changed.".into()
                            );
                        }
                        input.guard.run(())?;
                        let current = snapshot.run(())?;
                        if run_lineage::workspace_fingerprint(&current)
                            .map_err(|e| e.to_string())?
                            != fingerprint
                        {
                            return Err(
                                "The workspace changed while the product was being prepared."
                                    .into(),
                            );
                        }
                        Ok(())
                    };
                    guard()?;
                    let next =
                        crate::workspace::session::decode_workspace(archive.clone(), u32::MAX - 61)
                            .await?;
                    installed
                        .receipt
                        .validate_archive(&next)
                        .map_err(|e| e.to_string())?;
                    validate.run(next.clone())?;
                    let previous =
                        crate::workspace::session::encode_workspace(frame, u32::MAX - 62).await?;
                    guard()?;
                    crate::platform::storage::replace_workspace_cancellable(
                        &database,
                        &archive,
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
                        message.set("Opened reusable material with its exact source-event derivation record. No independent target verification or experimental success is asserted.".into());
                        on_restore.run(next);
                    }
                    Err(e) => {
                        if alive.load(Ordering::Relaxed) {
                            status.set(format!("The product was not opened. {e}"));
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
    view! {<section class="digest-open-panel"><p>"Opening freshly reconstructed material with its source and replay record. Nothing is opened until the recovery transaction succeeds."</p><p role="status">{move||status.get()}</p><Show when=move||failed.get()><button class="button primary" disabled=move||busy.get() on:click=move|_|run.with_value(|start|start())>"Retry product opening"</button></Show></section>}
}
