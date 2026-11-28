use super::{
    editors::{CandidateEditor, ComponentEditor, SubjectEditor},
    model::{self, Candidate, Capture, Draft, PcrDraft},
    rendering,
};
use crate::{platform::worker::CachedWorker, sources::identity::SourceInfo};
use collagene_core::{
    cloning::screening as native,
    documents::workspace::Workspace,
    protocol::{
        request::Request,
        schema::{AnalysisJob, AnalysisResult},
    },
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
pub(crate) fn ScreeningPanel(
    source: Memo<SourceInfo>,
    epoch: RwSignal<u64>,
    workspace_busy: RwSignal<bool>,
    get_workspace: Callback<(), Result<Workspace, String>>,
    on_select: Callback<Vec<(usize, usize)>>,
) -> impl IntoView {
    let draft = RwSignal::new(Draft::default());
    let workspace = RwSignal::new(None::<Arc<Workspace>>);
    let query = RwSignal::new(String::new());
    let event = RwSignal::new(String::new());
    let busy = RwSignal::new(false);
    let status = RwSignal::new(String::new());
    let generation = RwSignal::new(0u32);
    let report = RwSignal::new(None::<(Arc<Capture>, Arc<native::Report>)>);
    let selected_test = RwSignal::new(0usize);
    let worker = StoredValue::new_local(None::<Rc<CachedWorker>>);
    let alive = Arc::new(AtomicBool::new(true));
    let stop = Callback::new(move |()| {
        generation.update(|value| *value = value.wrapping_add(1));
        busy.set(false);
        report.set(None);
        if let Some(worker) = worker.get_value() {
            worker.close();
        }
        worker.set_value(None);
    });
    Effect::new(move |_| {
        let _ = (source.get(), epoch.get());
        stop.run(());
        workspace.set(None);
        draft.set(Draft::default());
        event.set(String::new());
    });
    Effect::new(move |_| {
        let _ = draft.get();
        stop.run(());
    });
    Effect::new(move |_| {
        if workspace_busy.get() {
            stop.run(());
        }
    });
    let cleanup = alive.clone();
    on_cleanup(move || {
        cleanup.store(false, Ordering::Relaxed);
        if let Some(worker) = worker.get_value() {
            worker.close();
        }
    });
    let load = Callback::new(move |()| {
        if workspace_busy.get_untracked() {
            return;
        }
        stop.run(());
        match get_workspace.run(()) {
            Ok(value) => {
                let active = source.get_untracked();
                if let Err(error) = model::molecule(&value, &active.binding) {
                    status.set(error);
                    return;
                }
                let mut next = Draft::default();
                if let Err(error) = next.add_component(&active) {
                    status.set(error);
                    return;
                }
                workspace.set(Some(Arc::new(value)));
                draft.set(next);
                event.set(String::new());
                status.set("Retained sources loaded. Add actual alternative components, define named outcomes, then choose candidate tests.".into());
            }
            Err(error) => status.set(error),
        }
    });
    let add_component = Callback::new(move |()| {
        let value = workspace.get_untracked();
        let Some(value) = value else {
            status.set("Load the retained sources first.".into());
            return;
        };
        let Some(info) = model::info(&value, &event.get_untracked()) else {
            status.set("Choose an exact retained reference event.".into());
            return;
        };
        if let Err(error) = draft
            .try_update(|d| d.add_component(&info))
            .unwrap_or_else(|| Err("The screening panel was closed.".into()))
        {
            status.set(error);
        }
    });
    let run_alive = alive.clone();
    let run = Callback::new(move |()| {
        if busy.get_untracked() || workspace_busy.get_untracked() {
            return;
        }
        let fresh = match get_workspace.run(()) {
            Ok(value) => value,
            Err(error) => {
                status.set(error);
                return;
            }
        };
        let capture = match Capture::prepare(
            &fresh,
            &source.get_untracked(),
            epoch.get_untracked(),
            draft.get_untracked(),
        ) {
            Ok(value) => Arc::new(value),
            Err(error) => {
                status.set(error);
                return;
            }
        };
        let active = match CachedWorker::new() {
            Ok(value) => Rc::new(value),
            Err(error) => {
                status.set(error);
                return;
            }
        };
        let token = generation.get_untracked().wrapping_add(1);
        generation.set(token);
        busy.set(true);
        report.set(None);
        worker.set_value(Some(active.clone()));
        status.set("Comparing explicit named outcomes under the ideal size model…".into());
        let alive = run_alive.clone();
        wasm_bindgen_futures::spawn_local(async move {
            let response = active
                .ask(Request {
                    id: token,
                    revision: capture.revision.clone(),
                    job: AnalysisJob::AssessScreeningStrategy {
                        request: capture.request.clone(),
                    },
                })
                .await;
            active.close();
            if !alive.load(Ordering::Relaxed) || generation.get_untracked() != token {
                return;
            }
            busy.set(false);
            worker.set_value(None);
            let received = (|| {
                if workspace_busy.get_untracked() {
                    return Err("The workspace is changing; run the analysis afterward.".into());
                }
                let current = get_workspace.run(())?;
                capture.guard(
                    &current,
                    &source.get_untracked(),
                    epoch.get_untracked(),
                    &draft.get_untracked(),
                )?;
                match response?.result {
                    AnalysisResult::ScreeningStrategy(value) => {
                        capture.validate_report(&value)?;
                        Ok(value)
                    }
                    AnalysisResult::Error(error) => Err(error),
                    _ => Err("Unexpected screening analysis response.".into()),
                }
            })();
            match received {
                Ok(value) => {
                    selected_test.set(value.ranked_test_indices.first().copied().unwrap_or(0));
                    report.set(Some((capture, Arc::new(value))));
                    status.set("Ideal-size comparisons ready. Ranking compares the declared size model; it is not an experimental recommendation.".into());
                }
                Err(error) => status.set(error),
            }
        });
    });
    let export_alive = alive.clone();
    let export = Callback::new(move |()| {
        if workspace_busy.get_untracked() {
            return;
        }
        let Some((capture, value)) = report.get_untracked() else {
            return;
        };
        let current = match get_workspace.run(()) {
            Ok(value) => value,
            Err(error) => {
                status.set(error);
                return;
            }
        };
        if let Err(error) = capture.guard(
            &current,
            &source.get_untracked(),
            epoch.get_untracked(),
            &draft.get_untracked(),
        ) {
            status.set(error);
            return;
        }
        let bytes = match capture.export_json(value.as_ref()) {
            Ok(bytes) => bytes,
            Err(error) => {
                status.set(error);
                return;
            }
        };
        let alive = export_alive.clone();
        let token = generation.get_untracked();
        wasm_bindgen_futures::spawn_local(async move {
            let result = crate::io::download::download_guarded(
                "screening-strategy.json",
                "application/json",
                &bytes,
                move || {
                    if !alive.load(Ordering::Relaxed)
                        || generation.try_get_untracked() != Some(token)
                        || workspace_busy.try_get_untracked() != Some(false)
                    {
                        return Err(
                            "The screening panel or workspace changed; export the current result."
                                .into(),
                        );
                    }
                    capture.guard(
                        &get_workspace.run(())?,
                        &source
                            .try_get_untracked()
                            .ok_or("The selected source is unavailable.")?,
                        epoch
                            .try_get_untracked()
                            .ok_or("The workspace is unavailable.")?,
                        &draft
                            .try_get_untracked()
                            .ok_or("The screening setup is unavailable.")?,
                    )
                },
            )
            .await;
            if let Err(error) = result {
                let _ = status.try_set(error);
            }
        });
    });
    let navigate = Callback::new(move |(cell, site): (usize, usize)| {
        if workspace_busy.get_untracked() {
            return;
        }
        let Some((capture, value)) = report.get_untracked() else {
            return;
        };
        let result = (|| {
            let current = get_workspace.run(())?;
            let active = source.get_untracked();
            capture.guard(
                &current,
                &active,
                epoch.get_untracked(),
                &draft.get_untracked(),
            )?;
            capture.recognition_ranges(&value, cell, site, &active)
        })();
        match result {
            Ok(ranges) => on_select.run(ranges),
            Err(error) => status.set(error),
        }
    });
    view! {
        <details class="biology-section screening-strategy">
            <summary>"Screening strategy"</summary>
            <p class="helper">
                "Compare explicitly named intended and alternative products, including mixtures, using ideal modeled sizes. "
                "This does not predict gel resolution, concentrations or experimental success."
            </p>
            <button disabled=move || workspace_busy.get() || busy.get() on:click=move |_| load.run(())>"Load retained sources / start with this reference"</button>
            <Show when=move || workspace.get().is_some()>
                <details class="biology-section" open>
                    <summary>"Retained components"</summary>
                    <label>"Find a source by name, document or event"
                        <input maxlength="160" prop:value=move || query.get() on:input=move |event| query.set(event_target_value(&event))/>
                    </label>
                    <label>"Exact retained reference event"
                        <select prop:value=move || event.get() on:change=move |value| event.set(event_target_value(&value))>
                            <option value="">"Choose a reference…"</option>
                            {move || workspace.with(|w| w.as_ref().map(|w| model::event_choices(w, &query.get(), &event.get()).into_iter().map(|s| view! {
                                <option value=s.binding.event_id.clone()>{format!("{} · {} bp · {}", s.name, s.length, s.binding.event_id)}</option>
                            }).collect_view()))}
                        </select>
                    </label>
                    <p class="helper">"At most 100 matching references are shown. Search by exact event to find an older reference. Separate named outcomes can contain identical DNA."</p>
                    <button disabled=move || draft.with(|d| d.components.len() >= native::MAX_COMPONENTS) on:click=move |_| add_component.run(())>"Add retained component"</button>
                    <For each=move || draft.with(|d| d.components.iter().map(|c| c.id.clone()).collect::<Vec<_>>()) key=|id| id.clone()
                        children=move |id| view! { <ComponentEditor id draft workspace epoch get_workspace/> }/>
                </details>
                <details class="biology-section" open>
                    <summary>"Named intended and alternative outcomes"</summary>
                    <For each=move || draft.with(|d| d.subjects.iter().map(|s| s.id.clone()).collect::<Vec<_>>()) key=|id| id.clone()
                        children=move |id| view! { <SubjectEditor id draft/> }/>
                    <button disabled=move || draft.with(|d| d.subjects.len() >= native::MAX_SUBJECTS)
                        on:click=move |_| {
                            if let Some(Err(error)) = draft.try_update(Draft::add_decoy) {
                                status.set(error);
                            }
                        }>"Add named alternative outcome"</button>
                </details>
                <details class="biology-section">
                    <summary>"Candidate digest and PCR tests"</summary>
                    <button disabled=move || busy.get() || workspace_busy.get()
                        on:click=move |_| {
                            match draft.try_update(Draft::add_catalogue_single_digests) {
                                Some(Ok(count)) => status.set(format!("Added {count} single-enzyme candidates. Review the tests and run the size comparison.")),
                                Some(Err(error)) => status.set(error),
                                None => (),
                            }
                        }>"Add catalogue single-enzyme digests"</button>
                    <For each=move || draft.with(|d| (0..d.candidates.len()).collect::<Vec<_>>()) key=|index| *index
                        children=move |index| view! { <CandidateEditor index draft/> }/>
                    <button disabled=move || draft.with(|d| d.candidates.len() >= native::MAX_TESTS)
                        on:click=move |_| draft.update(|d| {
                            if d.candidates.len() < native::MAX_TESTS {
                                d.candidates.push(Candidate::Digest {
                                    label: format!("Digest {}", d.candidates.len() + 1),
                                    enzymes: vec![],
                                });
                            }
                        })>"Add digest candidate"</button>
                    <button disabled=move || draft.with(|d| d.candidates.len() >= native::MAX_TESTS)
                        on:click=move |_| draft.update(|d| {
                            if d.candidates.len() < native::MAX_TESTS {
                                d.candidates.push(Candidate::Pcr(PcrDraft {
                                    name: format!("PCR {}", d.candidates.len() + 1),
                                    forward: String::new(),
                                    forward_tail: String::new(),
                                    reverse: String::new(),
                                    reverse_tail: String::new(),
                                    bottom: false,
                                }));
                            }
                        })>"Add PCR primer pair"</button>
                </details>
                <details class="biology-section">
                    <summary>"Ideal size comparison and search limits"</summary>
                    <p class="helper">
                        "Sizes are mean strand lengths. Exact sizes are the default; a nonzero tolerance is your declared comparison rule, not calibrated gel resolution. "
                        "Equal-size products and mixture components do not produce intensity evidence."
                    </p>
                    {[
                        ("Counterpart tolerance (half-bp units; 2 = 1 bp)", 0),
                        ("Shared restriction search work", 1),
                        ("Shared PCR search work", 2),
                        ("Panel comparison work", 3),
                    ].into_iter().map(move |(label, field)| view! {
                        <label>{label}
                            <input inputmode="numeric" maxlength="20"
                                prop:value=move || draft.with(|d| match field {
                                    0 => d.tolerance.clone(),
                                    1 => d.cleavage_work.clone(),
                                    2 => d.pcr_work.clone(),
                                    _ => d.panel_work.clone(),
                                })
                                on:input=move |event| draft.update(|d| {
                                    let value = event_target_value(&event);
                                    match field {
                                        0 => d.tolerance = value,
                                        1 => d.cleavage_work = value,
                                        2 => d.pcr_work = value,
                                        _ => d.panel_work = value,
                                    }
                                })/>
                        </label>
                    }).collect_view()}
                </details>
                <button disabled=move || busy.get() || workspace_busy.get() on:click=move |_| run.run(())>"Compare named outcomes"</button>
                <button disabled=move || !busy.get() on:click=move |_| { stop.run(()); status.set("Screening request cancelled.".into()); }>"Cancel"</button>
            </Show>
            <p role="status">{move || status.get()}</p>
            <Show when=move || report.get().is_some()>
                <button on:click=move |_| export.run(())>"Export complete JSON evidence"</button>
                <label>"Inspect candidate"
                    <select prop:value=move || selected_test.get().to_string()
                        on:change=move |event| { if let Ok(index) = event_target_value(&event).parse::<usize>() { selected_test.set(index); } }>
                        {move || report.get().map(|(_, value)| {
                            value.tests.iter().map(|test| view! {
                                <option value=test.test_index.to_string()>{test.definition.label.clone()}</option>
                            }).collect_view()
                        })}
                    </select>
                </label>
                {move || report.get().map(|(_, value)| rendering::results(value, selected_test.get(), navigate))}
            </Show>
        </details>
    }
}
