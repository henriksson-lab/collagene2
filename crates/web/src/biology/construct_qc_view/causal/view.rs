//! Lazy history selection and guarded native-worker comparison; no main-thread replay.
use super::{
    model::Capture,
    selection::{self, EventChoice, Form},
};
use crate::{
    io::download::ExportArtifact, platform::worker::CachedWorker, sources::identity::SourceInfo,
};
use collagene_core::{
    construct_qc::causal,
    documents::workspace::Workspace,
    protocol::{
        request::Request,
        schema::{AnalysisJob, AnalysisResult},
    },
    sequence::translation,
};
use leptos::prelude::*;
use std::{
    rc::Rc,
    sync::{
        Arc,
        atomic::{AtomicBool, Ordering},
    },
};

#[component]
pub(crate) fn CausalQcPanel(
    source: Memo<SourceInfo>,
    epoch: RwSignal<u64>,
    workspace_busy: RwSignal<bool>,
    get_workspace: Callback<(), Result<Workspace, String>>,
    on_select: Callback<Vec<(usize, usize)>>,
) -> impl IntoView {
    let form = RwSignal::new(Form::default());
    let loaded = RwSignal::new(None::<Arc<Workspace>>);
    let choices = RwSignal::new(Vec::<EventChoice>::new());
    let filter = RwSignal::new(String::new());
    let source_filter = RwSignal::new(String::new());
    let product_filter = RwSignal::new(String::new());
    let status = RwSignal::new(
        "Load retained history to choose a source CDS and compare it with the active document."
            .to_owned(),
    );
    let busy = RwSignal::new(false);
    let generation = RwSignal::new(0u32);
    let report = RwSignal::new(None::<(Capture, Arc<causal::Report>)>);
    let page = RwSignal::new(0usize);
    let worker = StoredValue::new_local(None::<Rc<CachedWorker>>);
    let alive = Arc::new(AtomicBool::new(true));
    Effect::new(move |_| {
        let _ = (source.get(), epoch.get());
        loaded.set(None);
        choices.set(vec![]);
        form.set(Form::default());
        filter.set(String::new());
        source_filter.set(String::new());
        product_filter.set(String::new());
    });
    Effect::new(move |_| {
        let _ = (source.get(), epoch.get(), form.get(), workspace_busy.get());
        generation.update(|value| *value = value.wrapping_add(1));
        if let Some(worker) = worker.get_value() {
            worker.close();
        }
        worker.set_value(None);
        busy.set(false);
        report.set(None);
        page.set(0);
    });
    {
        let alive = alive.clone();
        on_cleanup(move || {
            alive.store(false, Ordering::Relaxed);
            if let Some(worker) = worker.get_value() {
                worker.close();
            }
        });
    }
    let guard = Callback::new(move |capture: Capture| -> Result<(), String> {
        if workspace_busy.try_get_untracked() != Some(false) {
            return Err("The workspace is changing; try again afterward.".into());
        }
        let current = get_workspace.run(())?;
        capture.guard(
            &current,
            &source
                .try_get_untracked()
                .ok_or("The comparison panel was closed.")?
                .binding,
            epoch
                .try_get_untracked()
                .ok_or("The comparison panel was closed.")?,
            &form
                .try_get_untracked()
                .ok_or("The comparison panel was closed.")?,
        )
    });
    let load = move |_| {
        if workspace_busy.get_untracked() || busy.get_untracked() {
            return;
        }
        match get_workspace.run(()) {
            Ok(workspace) => {
                if let Err(error) = selection::molecule(&workspace, &source.get_untracked().binding)
                {
                    status.set(error);
                    return;
                }
                choices.set(selection::events(&workspace));
                loaded.set(Some(Arc::new(workspace)));
                form.set(Form::default());
                status.set("Choose a source version and both CDS annotations. State why their roles should be compared.".into());
            }
            Err(error) => status.set(error),
        }
    };
    let selected_source = Memo::new(move |_| form.with(|form| form.source_event.clone()));
    let source_cds = Memo::new(move |_| {
        let selected = selected_source.get();
        loaded
            .with(|workspace| {
                workspace.as_ref().and_then(|workspace| {
                    choices.with(|choices| {
                        choices
                            .iter()
                            .find(|choice| choice.binding.event_id == selected)
                            .map(|choice| &choice.binding)
                            .and_then(|binding| selection::molecule(workspace, binding).ok())
                            .map(selection::cds)
                    })
                })
            })
            .unwrap_or_default()
    });
    let product_cds = Memo::new(move |_| {
        loaded
            .with(|workspace| {
                workspace
                    .as_ref()
                    .and_then(|workspace| {
                        selection::molecule(workspace, &source.get().binding).ok()
                    })
                    .map(selection::cds)
            })
            .unwrap_or_default()
    });
    let run = {
        let alive = alive.clone();
        Callback::new(move |(): ()| {
            if workspace_busy.get_untracked() || busy.get_untracked() {
                return;
            }
            let prepared = loaded.with_untracked(|workspace| {
                let workspace = workspace.as_ref().ok_or("Load retained history first.")?;
                let capture = choices.with_untracked(|choices| {
                    Capture::prepare(
                        workspace,
                        choices,
                        source.get_untracked().binding,
                        epoch.get_untracked(),
                        form.get_untracked(),
                    )
                })?;
                guard.run(capture.clone())?;
                Ok::<_, String>((capture, workspace.clone()))
            });
            let (capture, workspace) = match prepared {
                Ok(value) => value,
                Err(error) => {
                    status.set(error);
                    return;
                }
            };
            let current_worker = match CachedWorker::new() {
                Ok(worker) => Rc::new(worker),
                Err(error) => {
                    status.set(error);
                    return;
                }
            };
            let token = generation.get_untracked().wrapping_add(1);
            generation.set(token);
            busy.set(true);
            report.set(None);
            page.set(0);
            worker.set_value(Some(current_worker.clone()));
            status.set("Reconstructing retained history and comparing CDS consequences…".into());
            let alive = alive.clone();
            wasm_bindgen_futures::spawn_local(async move {
                let response = current_worker
                    .ask(Request {
                        id: token,
                        revision: capture.revision.clone(),
                        job: AnalysisJob::AssessConstructCausalQc {
                            workspace: workspace.as_ref().clone(),
                            request: capture.request.clone(),
                        },
                    })
                    .await;
                current_worker.close();
                if !alive.load(Ordering::Relaxed) || generation.get_untracked() != token {
                    return;
                }
                worker.set_value(None);
                busy.set(false);
                let checked = (|| {
                    guard.run(capture.clone())?;
                    match response?.result {
                        AnalysisResult::ConstructCausalQc(value) => {
                            capture.validate_report(workspace.as_ref(), &value)?;
                            Ok(value)
                        }
                        AnalysisResult::Error(error) => Err(error),
                        _ => Err("Unexpected causal comparison response.".into()),
                    }
                })();
                match checked {
                    Ok(value) => {
                        report.set(Some((capture, Arc::new(value))));
                        status.set("History comparison ready. Review each finding and its assumptions; missing evidence remains indeterminate.".into());
                    }
                    Err(error) => status.set(error),
                }
            });
        })
    };
    let export = {
        let alive = alive.clone();
        Callback::new(move |(): ()| {
            let Some((capture, value)) = report.get_untracked() else {
                return;
            };
            let prepared = (|| {
                guard.run(capture.clone())?;
                let bytes =
                    serde_json::to_vec(value.as_ref()).map_err(|error| error.to_string())?;
                if bytes.len() > causal::MAX_REPORT_BYTES {
                    return Err("The comparison exceeds its bounded export size.".into());
                }
                ExportArtifact::prepare("causal-cds-qc.json", "application/json", &bytes)
            })();
            match prepared {
                Ok(artifact) => {
                    let alive = alive.clone();
                    artifact.complete_guarded(
                        move || {
                            if !alive.load(Ordering::Relaxed) {
                                return Err("The comparison panel was closed.".into());
                            }
                            guard.run(capture.clone())
                        },
                        move |result| {
                            if let Err(error) = result {
                                let _ = status.try_set(error);
                            }
                        },
                    );
                }
                Err(error) => status.set(error),
            }
        })
    };
    let select_alive = alive.clone();
    let select = Callback::new(move |(capture, ranges): (Capture, Vec<(usize, usize)>)| {
        if !select_alive.load(Ordering::Relaxed) {
            return;
        }
        match guard.run(capture) {
            Ok(()) => {
                let Some(current) = source.try_get_untracked() else {
                    return;
                };
                if selection::navigation_ranges_valid(&ranges, current.length) {
                    on_select.run(ranges);
                } else {
                    status.set(
                        "These coordinates are outside the current product; compare again.".into(),
                    );
                }
            }
            Err(error) => status.set(error),
        }
    });
    view! {
        <details class="qc-causal">
            <summary>"Compare CDS changes through history"</summary>
            <p class="helper">"Compare one explicitly paired source and product CDS. Actual retained edits and cloning operations establish correspondence; labels and sequence similarity do not. Annotation roles remain your assumptions."</p>
            <button class="button" disabled=move ||busy.get()||workspace_busy.get() on:click=load>
                "Load retained history"
            </button>
            <Show when=move ||loaded.with(Option::is_some)>
                <p class="helper">{move ||format!("Product: {} · {}",source.get().name,source.get().binding.event_id)}</p>
                <label>
                    "Find source version"
                    <input maxlength="128" prop:value=move ||filter.get()
                        on:input=move |event|filter.set(event_target_value(&event))/>
                </label>
                <label>
                    "Saved source version"
                    <select disabled=move ||busy.get()||workspace_busy.get()
                        prop:value=move ||form.with(|form|form.source_event.clone())
                        on:change=move |event|form.update(|form| {
                            form.source_event=event_target_value(&event);
                            form.source_feature.clear();
                            source_filter.set(String::new());
                        })>
                        <option value="">"Choose saved version"</option>
                        {move || {
                            let selected=form.with(|form|form.source_event.clone());
                            choices.with(|choices|selection::event_menu(choices,&filter.get(),&selected)
                                .into_iter().map(|choice|view! {
                                    <option value=choice.binding.event_id.clone()>{choice.label.clone()}</option>
                                }).collect_view())
                        }}
                    </select>
                </label>
                <p class="helper">{move ||format!(
                    "{} saved DNA versions. The menu shows at most 100 choices, including your selection; search by document name or version ID.",
                    choices.with(Vec::len))}</p>
                <Show when=move ||source_cds.with(|choices|choices.len()>100)>
                    <label>"Find source CDS"
                        <input maxlength="128" prop:value=move ||source_filter.get()
                            on:input=move |event|source_filter.set(event_target_value(&event))/>
                    </label>
                </Show>
                <label>
                    "Source CDS"
                    <select prop:value=move ||form.with(|form|form.source_feature.clone())
                        on:change=move |event|form.update(|form|form.source_feature=event_target_value(&event))>
                        <option value="">"Choose source annotation"</option>
                        {move ||source_cds.with(|choices|selection::cds_menu(choices,&source_filter.get(),&form.with(|form|form.source_feature.clone()))
                            .into_iter().map(|choice|view! {
                                <option value=choice.index.to_string()>{choice.label.clone()}</option>
                            }).collect_view())}
                    </select>
                </label>
                <Show when=move ||product_cds.with(|choices|choices.len()>100)>
                    <label>"Find product CDS"
                        <input maxlength="128" prop:value=move ||product_filter.get()
                            on:input=move |event|product_filter.set(event_target_value(&event))/>
                    </label>
                </Show>
                <label>
                    "Active product CDS"
                    <select prop:value=move ||form.with(|form|form.product_feature.clone())
                        on:change=move |event|form.update(|form|form.product_feature=event_target_value(&event))>
                        <option value="">"Choose product annotation"</option>
                        {move ||product_cds.with(|choices|selection::cds_menu(choices,&product_filter.get(),&form.with(|form|form.product_feature.clone()))
                            .into_iter().map(|choice|view! {
                                <option value=choice.index.to_string()>{choice.label.clone()}</option>
                            }).collect_view())}
                    </select>
                </label>
                <label>
                    "Genetic code"
                    <select prop:value=move ||form.with(|form|form.genetic_code.clone())
                        on:change=move |event|form.update(|form|form.genetic_code=event_target_value(&event))>
                        {translation::genetic_codes().into_iter().map(|code|view! {
                            <option value=code.id.to_string()>{format!("{} · {}",code.id,code.name)}</option>
                        }).collect_view()}
                    </select>
                </label>
                <label>
                    "Declared CDS role correspondence"
                    <input maxlength="1024"
                        placeholder="Why should these annotations represent the CDS roles being compared?"
                        prop:value=move ||form.with(|form|form.rationale.clone())
                        on:input=move |event|form.update(|form|form.rationale=event_target_value(&event))/>
                </label>
                <button class="button primary" disabled=move ||busy.get()||workspace_busy.get() on:click=move |_|run.run(())>
                    {move ||if busy.get(){"Comparing…"}else{"Compare through history"}}
                </button>
            </Show>
            <p class="helper" role="status">{move ||status.get()}</p>
            {move ||report.with(|stored|stored.as_ref().map(|(capture,report)| {
                report.pairs.iter().map(|pair|super::rendering::pair(pair,capture,page.get(),select)).collect_view()
            }))}
            <Show when=move ||report.with(Option::is_some)>
                <p class="helper">{move ||format!(
                    "Evidence page {} (20 frame segments and stop comparisons per page). Full details and recorded causes are included in JSON.",
                    page.get()+1)}</p>
                <button class="text-tool" disabled=move ||page.get()==0
                    on:click=move |_|page.update(|page|*page=page.saturating_sub(1))>
                    "Previous evidence"
                </button>
                <button class="text-tool" disabled=move ||report.with(|stored| {
                    stored.as_ref().is_none_or(|(_,report)|report.pairs.iter().all(|pair| {
                        (page.get()+1)*20>=pair.register_segments.len().max(pair.new_stops.len())
                    }))
                }) on:click=move |_|page.update(|page|*page+=1)>
                    "Next evidence"
                </button>
                <button class="button" on:click=move |_|export.run(())>"Export full history comparison JSON"</button>
                {move ||report.with(|stored|stored.as_ref().map(|(_,report)|view! {
                    <p class="helper">{report.statement.clone()}</p>
                }))}
            </Show>
        </details>
    }
}
