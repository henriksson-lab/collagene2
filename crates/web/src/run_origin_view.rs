//! Explicit source-bound saved-material origin choice, beneath preparation controls.
use crate::{
    cloning_methods_model::MethodsState, run_origin_model, source_identity::SourceInfo,
    worker::CachedWorker,
};
use collagene_core::{
    protocol::{AnalysisJob, AnalysisResult, Request},
    run_material::{self, Inspection, Selection},
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

#[derive(Clone, PartialEq, Eq)]
struct Context {
    source: SourceInfo,
    document: usize,
    epoch: u64,
}

#[component]
pub fn RunOriginPicker(
    config: RwSignal<MethodsState>,
    index: usize,
    source_info: Callback<usize, Option<SourceInfo>>,
    epoch: RwSignal<u64>,
    get_workspace: Callback<(), Result<Workspace, String>>,
) -> impl IntoView {
    let status = RwSignal::new(String::new());
    let choices = RwSignal::new(None::<Arc<Inspection>>);
    let busy = RwSignal::new(false);
    let generation = RwSignal::new(0u32);
    let active = StoredValue::new_local(None::<Rc<CachedWorker>>);
    let alive = Arc::new(AtomicBool::new(true));
    let cleanup = alive.clone();
    on_cleanup(move || {
        cleanup.store(false, Ordering::Relaxed);
        if let Some(worker) = active.get_value() {
            worker.close();
        }
    });
    let context = Memo::new(move |_| {
        let document = config.with(|s| s.rows.get(index).map(|r| r.document))?;
        Some(Context {
            source: source_info.run(document)?,
            document,
            epoch: epoch.get(),
        })
    });
    Effect::new(move |_| {
        let current = context.get();
        generation.update(|g| *g = g.wrapping_add(1));
        if let Some(worker) = active.get_value() {
            worker.close();
        }
        active.set_value(None);
        busy.set(false);
        choices.set(None);
        let stale = config.with_untracked(|s| {
            s.rows
                .get(index)
                .and_then(|r| r.run_origin.as_ref())
                .is_some_and(|pin| {
                    current
                        .as_ref()
                        .is_none_or(|c| pin.selection.event != c.source.binding)
                })
        });
        if stale {
            config.update(|s| {
                if let Some(row) = s.rows.get_mut(index) {
                    row.run_origin = None;
                }
            });
            status.set("The source event changed. Its saved material origin was cleared.".into());
        }
    });
    let run = {
        let alive = alive.clone();
        Callback::new(move |selection: Option<Selection>| {
            let Some(bound) = context.get_untracked() else {
                status.set("Choose an available exact source event.".into());
                return;
            };
            if selection
                .as_ref()
                .is_some_and(|s| s.event != bound.source.binding)
            {
                status.set("This origin belongs to a different source event.".into());
                return;
            }
            let workspace = match get_workspace.run(()) {
                Ok(w) => w,
                Err(e) => {
                    status.set(e);
                    return;
                }
            };
            let revision = match run_material::inspection_revision(
                &workspace,
                &bound.source.binding,
                &selection,
            ) {
                Ok(r) => r,
                Err(e) => {
                    status.set(e.to_string());
                    return;
                }
            };
            let worker = match CachedWorker::new() {
                Ok(w) => Rc::new(w),
                Err(e) => {
                    status.set(e);
                    return;
                }
            };
            if let Some(old) = active.get_value() {
                old.close();
            }
            active.set_value(Some(worker.clone()));
            let token = generation.get_untracked().wrapping_add(1);
            generation.set(token);
            busy.set(true);
            let requested = selection.clone();
            let alive = alive.clone();
            status.set(
                if selection.is_some() {
                    "Reconstructing the originating run against its original target…"
                } else {
                    "Finding retained origins for this exact source event…"
                }
                .into(),
            );
            spawn_local(async move {
                let response = worker
                    .ask(Request {
                        id: token,
                        revision: revision.clone(),
                        job: AnalysisJob::InspectRunMaterial {
                            workspace,
                            event: bound.source.binding.clone(),
                            selection,
                        },
                    })
                    .await;
                worker.close();
                if !alive.load(Ordering::Relaxed)
                    || generation.get_untracked() != token
                    || context.get_untracked() != Some(bound.clone())
                {
                    return;
                }
                active.set_value(None);
                busy.set(false);
                let result = match response {
                    Ok(response) if response.id == token && response.revision == revision => {
                        match response.result {
                            AnalysisResult::RunMaterialInspection(value) => Ok(value),
                            AnalysisResult::Error(e) => Err(e),
                            _ => Err("Unexpected saved-material response.".into()),
                        }
                    }
                    Ok(_) => Err("Stale saved-material response withheld.".into()),
                    Err(e) => Err(e),
                };
                match result {
                    Err(e) => status.set(e),
                    Ok(inspection) => {
                        if let Err(e) = inspection.validate_binding(&revision) {
                            status.set(e.to_string());
                            return;
                        }
                        if inspection.event != bound.source.binding {
                            status.set("The inspection belongs to another saved event.".into());
                            return;
                        }
                        if let Some(requested) = requested {
                            if inspection
                                .selected
                                .as_ref()
                                .is_some_and(|p| p.selection != requested)
                            {
                                status.set(
                                    "The reconstructed origin differs from the requested choice."
                                        .into(),
                                );
                                return;
                            }
                            if inspection.selected.is_some() {
                                let result = config.with_untracked(|s| {
                                    run_origin_model::apply(
                                        s,
                                        index,
                                        bound.document,
                                        &bound.source.binding,
                                        &inspection,
                                        &revision,
                                    )
                                });
                                match result {
                                    Ok(next) => {
                                        config.set(next);
                                        status.set("Saved material selected. Each method preview reconstructs the origin again before applying this row’s preparation.".into());
                                    }
                                    Err(e) => {
                                        status.set(e);
                                        return;
                                    }
                                }
                            } else {
                                status.set(inspection.findings.join(" "));
                            }
                        } else {
                            status.set(inspection.findings.join(" "));
                        }
                        choices.set(Some(Arc::new(inspection)));
                    }
                }
            });
        })
    };
    view!{<details class="saved-run-origin"><summary>{move||if config.with(|s|s.rows.get(index).is_some_and(|r|r.run_origin.is_some())){"Saved material · origin selected"}else{"Reuse a saved cloning product…"}}</summary>
        <p class="helper">"Select the actual saved product’s originating run. Its recipe is reconstructed before this row’s PCR or digest; saved end-state metadata is not used as material."</p>
        {move||config.with(|s|s.rows.get(index).and_then(|r|r.run_origin.clone())).map(|pin|view!{<p title=pin.selection.run_id.clone()>{format!("Run {}… · product {}",&pin.selection.run_id[..12],pin.selection.step_index+1)}</p><button class="text-tool" disabled=move||busy.get() on:click=move|_|{generation.update(|g|*g=g.wrapping_add(1));if let Some(w)=active.get_value(){w.close();}active.set_value(None);busy.set(false);config.update(|s|if let Some(row)=s.rows.get_mut(index){row.run_origin=None;});status.set("Origin removed. The source is again only a reference input.".into());}>"Remove saved-material origin"</button>})}
        <button class="button secondary" disabled=move||busy.get()||context.get().is_none() on:click=move|_|run.run(None)>"Find originating runs"</button>
        {move||choices.get().map(|inspection|inspection.origins.iter().take(100).cloned().map(|selection|{let label=format!("Reconstruct run {}… · product {}",&selection.run_id[..12],selection.step_index+1);let title=selection.run_id.clone();view!{<button class="text-tool" title=title disabled=move||busy.get() on:click=move|_|run.run(Some(selection.clone()))>{label}</button>}}).collect_view())}
        <p class="helper" role="status" aria-live="polite">{move||status.get()}</p>
    </details>}.into_any()
}
