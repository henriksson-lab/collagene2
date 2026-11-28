//! Worker-reviewed sequence exchange; saving awaits the real guarded file receipt.
use crate::{
    sequence_export_model::{Capture, Mode, Review},
    worker::CachedWorker,
};
use collagene_core::protocol::{
    request::Request,
    schema::{AnalysisJob, AnalysisResult},
};
use collagene_model::Molecule;
use leptos::prelude::*;
use std::{
    rc::Rc,
    sync::{
        Arc,
        atomic::{AtomicBool, Ordering},
    },
};
use wasm_bindgen_futures::spawn_local;
#[component]
pub fn SequenceCompatExport(
    source: Molecule,
    mode: Mode,
    guard: Callback<(), Result<(), String>>,
    on_saved: Callback<String>,
) -> impl IntoView {
    let id = RwSignal::new(source.name.clone());
    let date = RwSignal::new(
        js_sys::Date::new_0()
            .to_iso_string()
            .as_string()
            .unwrap_or_default()
            .chars()
            .take(10)
            .collect::<String>(),
    );
    let allow = RwSignal::new(false);
    let error = RwSignal::new(String::new());
    let review = RwSignal::new(None::<Arc<Review>>);
    let busy = RwSignal::new(false);
    let generation = RwSignal::new(0u32);
    let worker = StoredValue::new_local(None::<Rc<CachedWorker>>);
    let source = Arc::new(source);
    let alive = Arc::new(AtomicBool::new(true));
    let cleanup = alive.clone();
    on_cleanup(move || {
        cleanup.store(false, Ordering::Relaxed);
        if let Some(w) = worker.get_value() {
            w.close();
        }
    });
    let invalidate = Callback::new(move |()| {
        generation.update(|n| *n = n.wrapping_add(1));
        if let Some(w) = worker.get_value() {
            w.close();
        }
        worker.set_value(None);
        busy.set(false);
        review.set(None);
        allow.set(false);
        error.set(String::new());
    });
    let inspect = {
        let source = source.clone();
        let alive = alive.clone();
        move |_| {
            if busy.get_untracked() {
                return;
            }
            invalidate.run(());
            if let Err(e) = guard.run(()) {
                error.set(e);
                return;
            }
            let capture = match Capture::new(
                &source,
                mode.format(&id.get_untracked(), &date.get_untracked()),
            ) {
                Ok(c) => c,
                Err(e) => {
                    error.set(e);
                    return;
                }
            };
            let w = match CachedWorker::new() {
                Ok(w) => Rc::new(w),
                Err(e) => {
                    error.set(e);
                    return;
                }
            };
            worker.set_value(Some(w.clone()));
            busy.set(true);
            let token = generation.get_untracked();
            let source = source.clone();
            let alive = alive.clone();
            spawn_local(async move {
                let response = w
                    .ask(Request {
                        id: token,
                        revision: capture.context_revision.clone(),
                        job: AnalysisJob::ReviewSequenceExport {
                            molecule: source.as_ref().clone(),
                            format: capture.format.clone(),
                        },
                    })
                    .await;
                w.close();
                if !alive.load(Ordering::Relaxed) || generation.try_get_untracked() != Some(token) {
                    return;
                }
                worker.set_value(None);
                busy.set(false);
                if let Err(e) = guard.run(()) {
                    error.set(e);
                    return;
                }
                if !capture.matches(
                    &source,
                    &mode.format(&id.get_untracked(), &date.get_untracked()),
                ) {
                    error.set("Export options changed; review again.".into());
                    return;
                }
                let result = (|| -> Result<Review, String> {
                    let response = response?;
                    if response.id != token || response.revision != capture.context_revision {
                        return Err("Export assessment response belongs to another request.".into());
                    }
                    let assessment = match response.result {
                        AnalysisResult::SequenceExportAssessment(value) => value,
                        AnalysisResult::Error(e) => return Err(e),
                        _ => return Err("Unexpected export assessment response.".into()),
                    };
                    Review::new(capture, assessment)
                })();
                match result {
                    Ok(value) => review.set(Some(Arc::new(value))),
                    Err(e) => error.set(e),
                }
            });
        }
    };
    let save = {
        let source = source.clone();
        let alive = alive.clone();
        move |_| {
            if busy.get_untracked() {
                return;
            }
            let Some(bound) = review.get_untracked() else {
                return;
            };
            let accepted = allow.get_untracked();
            if let Err(e) = guard.run(()).and_then(|()| bound.allow_prepare(accepted)) {
                error.set(e);
                return;
            }
            if !bound.capture.matches(
                &source,
                &mode.format(&id.get_untracked(), &date.get_untracked()),
            ) {
                invalidate.run(());
                error.set("Export options changed; review again.".into());
                return;
            }
            let w = match CachedWorker::new() {
                Ok(w) => Rc::new(w),
                Err(e) => {
                    error.set(e);
                    return;
                }
            };
            worker.set_value(Some(w.clone()));
            busy.set(true);
            error.set(String::new());
            let token = generation.get_untracked();
            let source = source.clone();
            let alive = alive.clone();
            spawn_local(async move {
                let response = w
                    .ask(Request {
                        id: token,
                        revision: bound.capture.context_revision.clone(),
                        job: AnalysisJob::PrepareSequenceExport {
                            molecule: source.as_ref().clone(),
                            format: bound.capture.format.clone(),
                            expected_context_revision: bound.capture.context_revision.clone(),
                            allow_omissions: accepted,
                        },
                    })
                    .await;
                w.close();
                if !alive.load(Ordering::Relaxed) || generation.try_get_untracked() != Some(token) {
                    return;
                }
                worker.set_value(None);
                let current = {
                    let bound = bound.clone();
                    let source = source.clone();
                    let alive = alive.clone();
                    move || -> Result<(), String> {
                        if !alive.load(Ordering::Relaxed)
                            || generation.try_get_untracked() != Some(token)
                        {
                            return Err("Export review closed or was replaced.".into());
                        }
                        guard.run(())?;
                        if !bound.capture.matches(
                            &source,
                            &mode.format(&id.get_untracked(), &date.get_untracked()),
                        ) || allow.get_untracked() != accepted
                        {
                            return Err(
                                "Export options or omission approval changed; review again.".into(),
                            );
                        }
                        if !review
                            .with_untracked(|r| r.as_ref().is_some_and(|r| Arc::ptr_eq(r, &bound)))
                        {
                            return Err("Export assessment was replaced; review again.".into());
                        }
                        Ok(())
                    }
                };
                let result = async {
                    current()?;
                    let response = response?;
                    if response.id != token || response.revision != bound.capture.context_revision {
                        return Err("Prepared export response belongs to another request.".into());
                    }
                    let prepared = match response.result {
                        AnalysisResult::SequenceExportPrepared(value) => value,
                        AnalysisResult::Error(e) => return Err(e),
                        _ => return Err("Unexpected prepared export response.".into()),
                    };
                    bound.validate_prepared(&prepared, accepted)?;
                    current()?;
                    crate::app::download_guarded(
                        &prepared.assessment.filename,
                        &prepared.assessment.mime,
                        &prepared.bytes,
                        current,
                    )
                    .await?;
                    Ok::<_, String>(prepared.assessment.filename)
                }
                .await;
                if !alive.load(Ordering::Relaxed) || generation.try_get_untracked() != Some(token) {
                    return;
                }
                busy.set(false);
                match result {
                    Ok(name) => {
                        if guard.run(()).is_ok() {
                            on_saved.run(format!("Saved {name}. Native history and original evidence remain in the native document."));
                        }
                    }
                    Err(e) => error.set(e),
                }
            });
        }
    };
    view! {<div class="sequence-compat-export">
        <p class="helper">{format!("{} reference exchange. Review representability and omissions before preparing the file. Native history, original attachments and physical reconstruction remain in the native document.",mode.label())}</p>
        <Show when=move||mode==Mode::Gcg><label>"Record ID"<input prop:value=move||id.get() disabled=move||busy.get() on:input=move|e|{id.set(event_target_value(&e));invalidate.run(());}/></label><label>"Export date (UTC)"<input prop:value=move||date.get() disabled=move||busy.get() on:input=move|e|{date.set(event_target_value(&e));invalidate.run(());}/></label></Show>
        <button class="button secondary" disabled=move||busy.get() on:click=inspect>"Review export"</button>
        {move||review.get().map(|bound|{let assessment=&bound.assessment;view!{<p>{format!("{} · {}",assessment.profile,assessment.filename)}</p><ul>{assessment.losses.iter().map(|loss|view!{<li>{format!("{}: {}",loss.path,loss.message)}</li>}).collect_view()}{assessment.warnings.iter().map(|notice|view!{<li>{notice.clone()}</li>}).collect_view()}</ul>}})}
        <label><input type="checkbox" prop:checked=move||allow.get() disabled=move||busy.get()||review.get().is_none() on:change=move|e|allow.set(event_target_checked(&e))/>"I accept the listed omissions for this exchange file"</label>
        <p role="alert">{move||error.get()}</p>
        <button class="button" disabled=move||busy.get()||review.get().is_none_or(|r|r.allow_prepare(allow.get()).is_err()) on:click=save>{move||if busy.get(){"Preparing…"}else{"Save reviewed export"}}</button>
    </div>}
}
