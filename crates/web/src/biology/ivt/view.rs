//! Explicit IVT template controls in the Biology inspector.
use super::model::{Capture, Draft};
use crate::{platform::worker::CachedWorker, sources::identity::SourceInfo};
use collagene_core::{
    protocol::{
        request::Request,
        schema::{AnalysisJob, AnalysisResult},
    },
    rna_chemistry::ivt,
};
use collagene_model::{Location, Molecule};
use leptos::prelude::*;
use std::{
    rc::Rc,
    sync::{
        Arc,
        atomic::{AtomicBool, Ordering},
    },
};

fn coordinate(
    draft: RwSignal<Draft>,
    label: &'static str,
    get: fn(&Draft) -> &str,
    set: fn(&mut Draft, String),
) -> impl IntoView {
    view! {
        <label>{label}
            <input inputmode="numeric" maxlength="20"
                prop:value=move ||draft.with(|d|get(d).to_string())
                on:input=move |e|draft.update(|d|set(d,event_target_value(&e)))/>
        </label>
    }
}
fn source_location(spans: &[ivt::MappedSpan]) -> Option<Location> {
    let mut locations = spans
        .iter()
        .map(|span| {
            let range = Location::range(span.start, span.end);
            if span.reverse {
                Location::Complement(Box::new(range))
            } else {
                range
            }
        })
        .collect::<Vec<_>>();
    match locations.len() {
        0 => None,
        1 => locations.pop(),
        _ => Some(Location::Join(locations)),
    }
}

#[component]
pub(crate) fn IvtPanel(
    molecule: Memo<Molecule>,
    source: Memo<SourceInfo>,
    epoch: RwSignal<u64>,
    workspace_busy: RwSignal<bool>,
    on_select: Callback<(String, Location), Result<(), String>>,
) -> impl IntoView {
    let draft = RwSignal::new(Draft::default());
    let status = RwSignal::new(String::new());
    let busy = RwSignal::new(false);
    let generation = RwSignal::new(0u32);
    let report = RwSignal::new(None::<(Capture, Arc<ivt::Report>)>);
    let worker = StoredValue::new_local(None::<Rc<CachedWorker>>);
    let alive = Arc::new(AtomicBool::new(true));
    Effect::new(move |_| {
        source.get();
        epoch.get();
        draft.set(Draft::default());
    });
    Effect::new(move |_| {
        let _ = (source.get(), epoch.get(), draft.get());
        generation.update(|n| *n = n.wrapping_add(1));
        if let Some(w) = worker.get_value() {
            w.close();
        }
        worker.set_value(None);
        busy.set(false);
        report.set(None);
        status.set(
            "Declare the transcript interval and requirements, then check this template.".into(),
        );
    });
    {
        let alive = alive.clone();
        on_cleanup(move || {
            alive.store(false, Ordering::Relaxed);
            if let Some(w) = worker.get_value() {
                w.close();
            }
        });
    }
    let run = {
        let alive = alive.clone();
        move |_| {
            if busy.get_untracked() || workspace_busy.get_untracked() {
                return;
            }
            let captured = molecule.with_untracked(|m| {
                ivt::source_preflight(m).map_err(|error| error.to_string())?;
                let captured = Capture::prepare(
                    source.get_untracked().binding,
                    epoch.get_untracked(),
                    draft.get_untracked(),
                    m,
                )?;
                Ok::<_, String>((m.clone(), captured))
            });
            let (m, captured) = match captured {
                Ok(value) => value,
                Err(error) => {
                    status.set(error);
                    return;
                }
            };
            let w = match CachedWorker::new() {
                Ok(value) => Rc::new(value),
                Err(error) => {
                    status.set(error);
                    return;
                }
            };
            worker.set_value(Some(w.clone()));
            let token = generation.get_untracked().wrapping_add(1);
            generation.set(token);
            busy.set(true);
            report.set(None);
            status.set("Checking the declared IVT template…".into());
            let alive = alive.clone();
            wasm_bindgen_futures::spawn_local(async move {
                let response = w
                    .ask(Request {
                        id: token,
                        revision: captured.revision.clone(),
                        job: AnalysisJob::AssessIvtTemplate {
                            molecule: m,
                            request: captured.request.clone(),
                        },
                    })
                    .await;
                w.close();
                if !alive.load(Ordering::Relaxed) || generation.get_untracked() != token {
                    return;
                }
                busy.set(false);
                worker.set_value(None);
                let received = (|| {
                    captured.guard(
                        &source.get_untracked().binding,
                        epoch.get_untracked(),
                        &draft.get_untracked(),
                    )?;
                    match response?.result {
                        AnalysisResult::IvtTemplate(value) => {
                            captured.validate_report(&value)?;
                            Ok(value)
                        }
                        AnalysisResult::Error(error) => Err(error),
                        _ => Err("Unexpected IVT template response.".into()),
                    }
                })();
                match received {
                    Ok(value) => {
                        report.set(Some((captured, Arc::new(value))));
                        status.set("Template checks ready. Review each conflict, assumption and source below.".into());
                    }
                    Err(error) => status.set(error),
                }
            });
        }
    };
    let export = move |_| {
        let Some((capture, value)) = report.get_untracked() else {
            return;
        };
        if workspace_busy.get_untracked() {
            return;
        }
        if let Err(error) = capture.guard(
            &source.get_untracked().binding,
            epoch.get_untracked(),
            &draft.get_untracked(),
        ) {
            status.set(error);
            return;
        }
        let bytes = match serde_json::to_vec_pretty(value.as_ref()) {
            Ok(value) => value,
            Err(error) => {
                status.set(error.to_string());
                return;
            }
        };
        wasm_bindgen_futures::spawn_local(async move {
            let saved = crate::io::download::download_guarded(
                "ivt-template.json",
                "application/json",
                &bytes,
                move || {
                    if workspace_busy.try_get_untracked() != Some(false) {
                        return Err("The workspace is changing; export again afterward.".into());
                    }
                    let current = source
                        .try_get_untracked()
                        .ok_or("The IVT panel was closed.")?;
                    let current_epoch = epoch
                        .try_get_untracked()
                        .ok_or("The IVT panel was closed.")?;
                    let current_draft = draft
                        .try_get_untracked()
                        .ok_or("The IVT panel was closed.")?;
                    capture.guard(&current.binding, current_epoch, &current_draft)
                },
            )
            .await;
            if let Err(error) = saved {
                let _ = status.try_set(error);
            }
        });
    };
    view! {
        <details class="biology-section">
            <summary>"IVT template"</summary>
            <p class="helper">"Check a declared T7 coding strand and transcript interval. Coordinates are 1-based; reverse coding proceeds toward decreasing source positions."</p>
            {coordinate(draft,"First transcribed source base (+1)",|d|&d.first_base,|d,v|d.first_base=v)}
            {coordinate(draft,"Transcript length (bases)",|d|&d.length,|d,v|d.length=v)}
            <label>"Coding strand"
                <select prop:value=move ||if draft.get().reverse {"reverse"} else {"forward"}
                    on:change=move |e|draft.update(|d|d.reverse=event_target_value(&e)=="reverse")>
                    <option value="forward">"Forward"</option><option value="reverse">"Reverse"</option>
                </select>
            </label>
            <label>"Cap chemistry"
                <select prop:value=move ||match draft.get().cap {
                    ivt::CapChoice::NotSpecified=>"unknown", ivt::CapChoice::Uncapped=>"uncapped",
                    ivt::CapChoice::EnzymaticPostTranscriptional=>"enzymatic",ivt::CapChoice::CleanCapAgN7113=>"ag",
                } on:change=move |e|draft.update(|d|d.cap=match event_target_value(&e).as_str() {
                    "uncapped"=>ivt::CapChoice::Uncapped,"enzymatic"=>ivt::CapChoice::EnzymaticPostTranscriptional,
                    "ag"=>ivt::CapChoice::CleanCapAgN7113,_=>ivt::CapChoice::NotSpecified,
                })>
                    <option value="unknown">"Unspecified"</option><option value="uncapped">"Uncapped"</option>
                    <option value="enzymatic">"Enzymatic, after transcription"</option><option value="ag">"CleanCap AG · original N-7113"</option>
                </select>
            </label>
            <label>"Additional declared initiator requirement"
                <select prop:value=move ||match draft.get().initiator {
                    ivt::InitiatorRequirement::NoAdditionalRequirement=>"none",ivt::InitiatorRequirement::G=>"g",
                    ivt::InitiatorRequirement::Gg=>"gg",ivt::InitiatorRequirement::Ggg=>"ggg",
                } on:change=move |e|draft.update(|d|d.initiator=match event_target_value(&e).as_str() {
                    "g"=>ivt::InitiatorRequirement::G,"gg"=>ivt::InitiatorRequirement::Gg,
                    "ggg"=>ivt::InitiatorRequirement::Ggg,_=>ivt::InitiatorRequirement::NoAdditionalRequirement,
                })>
                    <option value="none">"None"</option><option value="g">"G"</option>
                    <option value="gg">"GG"</option><option value="ggg">"GGG"</option>
                </select>
            </label>
            <details><summary>"Runoff geometry and encoded poly(A)"</summary>
                <p class="helper">"Optional cut offsets count bonds from +1 along the coding direction: 0 is before +1. Enter both strand cuts; this declares geometry without simulating digestion."</p>
                {coordinate(draft,"Coding-strand cut offset",|d|&d.coding_cut,|d,v|d.coding_cut=v)}
                {coordinate(draft,"Template-strand cut offset",|d|&d.template_cut,|d,v|d.template_cut=v)}
                <p class="helper">"Optional poly(A) interval uses transcript positions, 1-based and inclusive. The check reports exact encoding and interruptions."</p>
                {coordinate(draft,"First poly(A) transcript base",|d|&d.poly_a_first,|d,v|d.poly_a_first=v)}
                {coordinate(draft,"Last poly(A) transcript base",|d|&d.poly_a_last,|d,v|d.poly_a_last=v)}
            </details>
            <button class="text-tool" disabled=move ||busy.get()||workspace_busy.get() on:click=run>"Check template"</button>
            <p class="helper" role="status">{move ||status.get()}</p>
            {move ||report.get().map(|(capture,r)| {
                let revision = capture.binding.revision.clone();
                let checks = r.findings.iter().map(|finding| {
                    let label = format!("{:?} · {:?}",finding.rule,finding.verdict);
                    let message = finding.message.clone();
                    let location = source_location(&finding.locations);
                    let revision = revision.clone();
                    let saved_capture = capture.clone();
                    let references = finding.source_ids.iter().filter_map(|id|r.profile.sources.iter().find(|s|s.id==*id))
                        .map(|s|view! {<p class="helper"><a href=s.url.clone() target="_blank" rel="noopener noreferrer">{s.id.clone()}</a>{format!(" · {}",s.locator)}</p>}).collect_view();
                    view! {<article class="orf-entry"><strong>{label}</strong><p>{message}</p>{references}
                        {location.map(|location|view! {<button class="text-tool" on:click=move |_| {
                            let selected = (|| {
                                if workspace_busy.try_get_untracked() != Some(false) {
                                    return Err("The workspace is changing; select again afterward.".into());
                                }
                                let current = source.try_get_untracked().ok_or("The IVT panel was closed.")?;
                                let current_epoch = epoch.try_get_untracked().ok_or("The IVT panel was closed.")?;
                                let current_draft = draft.try_get_untracked().ok_or("The IVT panel was closed.")?;
                                saved_capture.guard(&current.binding,current_epoch,&current_draft)?;
                                on_select.run((revision.clone(),location.clone()))
                            })();
                            if let Err(error) = selected {let _ = status.try_set(error);}
                        }>"Show region"</button>})}
                    </article>}
                }).collect_view();
                let evidence = serde_json::to_string_pretty(r.as_ref()).unwrap_or_else(|e|e.to_string());
                view! {<div>{checks}
                    <details><summary>"Assumptions and complete evidence"</summary>
                        {r.assumptions.iter().map(|text|view! {<p class="helper">{text.clone()}</p>}).collect_view()}
                        <pre>{evidence}</pre>
                    </details>
                    <button class="text-tool" on:click=export>"Export IVT checks (JSON)"</button>
                </div>}
            })}
        </details>
    }
}
