//! Explicit source-bound saved-material origin choice, beneath preparation controls.
use crate::{
    cloning_methods_model::{MethodsState, OverlapPreparation},
    run_origin_model::PinnedOrigin,
    source_identity::SourceInfo,
    worker::CachedWorker,
};
use collagene_core::{
    protocol::{AnalysisJob, AnalysisResult, Request},
    run_material::{self, MaterialInspection, MaterialSelection},
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
    let row_source = Memo::new(move |_| {
        config.with(|s| s.rows.get(index).and_then(|r| source_info.run(r.document)))
    });
    let selected =
        Memo::new(move |_| config.with(|s| s.rows.get(index).and_then(|r| r.run_origin.clone())));
    Effect::new(move |_| {
        let current = row_source.get();
        if selected.get_untracked().is_some_and(|pin| {
            current
                .as_ref()
                .is_none_or(|source| *pin.event() != source.binding)
        }) {
            config.update(|state| {
                if let Some(row) = state.rows.get_mut(index) {
                    row.run_origin = None;
                }
            });
        }
    });
    let on_select = Callback::new(move |pin: Option<PinnedOrigin>| {
        let current = row_source
            .get_untracked()
            .ok_or("The source row is unavailable.")?;
        let mut next = config.get_untracked();
        let row = next
            .rows
            .get_mut(index)
            .ok_or("The preparation row is unavailable.")?;
        if let Some(pin) = &pin {
            if *pin.event() != current.binding
                || matches!(
                    row.preparation,
                    OverlapPreparation::DeclaredDuplex | OverlapPreparation::AnnealedOligos
                )
            {
                return Err(
                    "The source or preparation changed; select saved material again.".into(),
                );
            }
        }
        row.run_origin = pin;
        config.set(next);
        Ok(())
    });
    view! {<SourceOriginPicker source=row_source epoch get_workspace selected on_select/>}
}

#[component]
pub fn SourceOriginPicker(
    source: Memo<Option<SourceInfo>>,
    epoch: RwSignal<u64>,
    get_workspace: Callback<(), Result<Workspace, String>>,
    selected: Memo<Option<PinnedOrigin>>,
    on_select: Callback<Option<PinnedOrigin>, Result<(), String>>,
) -> impl IntoView {
    let status = RwSignal::new(String::new());
    let choices = RwSignal::new(None::<Arc<MaterialInspection>>);
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
        Some(Context {
            source: source.get()?,
            epoch: epoch.get(),
        })
    });
    Effect::new(move |_| {
        let current = context.get();
        let selected_pin = selected.get();
        generation.update(|g| *g = g.wrapping_add(1));
        if let Some(worker) = active.get_value() {
            worker.close();
        }
        active.set_value(None);
        busy.set(false);
        choices.set(None);
        if selected_pin.is_some_and(|pin| {
            current
                .as_ref()
                .is_none_or(|c| *pin.event() != c.source.binding)
        }) {
            if let Err(error) = on_select.run(None) {
                status.set(error);
            } else {
                status
                    .set("The source event changed. Its saved material origin was cleared.".into());
            }
        }
    });
    let run = {
        let alive = alive.clone();
        Callback::new(move |selection: Option<MaterialSelection>| {
            let Some(bound) = context.get_untracked() else {
                status.set("Choose an available exact source event.".into());
                return;
            };
            if selection
                .as_ref()
                .is_some_and(|s| *s.event() != bound.source.binding)
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
            let revision = match run_material::material_inspection_revision(
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
                        job: AnalysisJob::InspectMaterialOrigin {
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
                            AnalysisResult::MaterialOriginInspection(value) => Ok(value),
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
                                .is_some_and(|p| p.selection() != requested)
                            {
                                status.set(
                                    "The reconstructed origin differs from the requested choice."
                                        .into(),
                                );
                                return;
                            }
                            if inspection.selected.is_some() {
                                let result = inspection
                                    .selected
                                    .as_ref()
                                    .cloned()
                                    .ok_or("No fresh material was selected.".to_string())
                                    .and_then(PinnedOrigin::from_material)
                                    .and_then(|pin| on_select.run(Some(pin)));
                                match result {
                                    Ok(()) => {
                                        status.set("Saved material selected. Each new inspection or preparation reconstructs this origin again.".into());
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
    view!{<details class="saved-run-origin"><summary>{move||if selected.with(Option::is_some){"Saved material · origin selected"}else{"Reuse a saved cloning product…"}}</summary>
        <p class="helper">"Select the actual saved product’s retained material origin. Its recorded cloning steps are checked again before further preparation, retaining its ends and modifications."</p>
        {move||selected.get().map(|pin|view!{<p>{pin.label()}</p><button class="text-tool" disabled=move||busy.get() on:click=move|_|{generation.update(|g|*g=g.wrapping_add(1));if let Some(w)=active.get_value(){w.close();}active.set_value(None);busy.set(false);match on_select.run(None){Ok(())=>status.set("Origin removed. The source is again only a reference input.".into()),Err(error)=>status.set(error)};}>"Remove saved-material origin"</button>})}
        <button class="button secondary" disabled=move||busy.get()||context.get().is_none() on:click=move|_|run.run(None)>"Find saved material origins"</button>
        {move||choices.get().map(|inspection|inspection.origins.iter().take(100).cloned().map(|selection|{let (label,title)=match &selection{MaterialSelection::Run(s)=>(format!("Reconstruct run {}… · product {}",s.run_id.chars().take(12).collect::<String>(),s.step_index+1),s.run_id.clone()),MaterialSelection::Digest(e)=>("Reconstruct saved digest fragment".to_string(),e.event_id.clone()),MaterialSelection::Derived(e)=>("Reconstruct reusable PCR/digest material".to_string(),e.event_id.clone())};view!{<button class="text-tool" title=title disabled=move||busy.get() on:click=move|_|run.run(Some(selection.clone()))>{label}</button>}}).collect_view())}
        <p class="helper" role="status" aria-live="polite">{move||status.get()}</p>
    </details>}.into_any()
}
