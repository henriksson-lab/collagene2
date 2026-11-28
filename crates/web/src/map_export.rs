//! Source-pinned portable map export through the native Rust worker.
use crate::map_export_model::{Capture, Prepared, RenderRequest, Rendered, Snapshot};
use leptos::prelude::*;
use std::{
    future::Future,
    pin::Pin,
    sync::{
        Arc,
        atomic::{AtomicBool, Ordering},
    },
};
use wasm_bindgen_futures::spawn_local;
pub type LocalResult<T> = Pin<Box<dyn Future<Output = Result<T, String>>>>;
/// One isolated, timeout-bounded native Rust worker per explicit preparation.
pub fn worker_renderer() -> Callback<RenderRequest, LocalResult<Rendered>> {
    let counter = Arc::new(std::sync::atomic::AtomicU32::new(1));
    Callback::new(move |request: RenderRequest| -> LocalResult<Rendered> {
        let id = counter.fetch_add(1, Ordering::Relaxed);
        Box::pin(async move {
            use collagene_core::protocol::{
                request::Request,
                schema::{AnalysisJob, AnalysisResult},
            };
            let context = request.context_revision.clone();
            let source_revision = request.source_revision.clone();
            let response = crate::worker::analyze(Request {
                id,
                revision: context.clone(),
                job: AnalysisJob::RenderMapSvg {
                    request: collagene_core::svg_export::Request {
                        source: request.source,
                        source_revision: source_revision.clone(),
                        restriction_sites: request.restriction_sites,
                    },
                },
            })
            .await?;
            if response.id != id || response.revision != context {
                return Err("Map worker response belongs to a different captured view.".into());
            }
            match response.result {
                AnalysisResult::MapSvgExported(value)
                    if value.source_revision == source_revision =>
                {
                    value
                        .validate_binding(&source_revision)
                        .map_err(|error| error.to_string())?;
                    Ok(Rendered {
                        source_revision: value.source_revision,
                        context_revision: context,
                        svg: value.svg,
                        warnings: value.warnings,
                    })
                }
                AnalysisResult::Error(error) => Err(error),
                _ => Err("Unexpected map export response; no file was prepared.".into()),
            }
        })
    })
}

/// The renderer callback must use a worker/native adapter, not synchronous SVG
/// drawing in an event handler. Snapshot reads the exact live document AND view.
#[component]
pub fn MapExportDialog(
    snapshot: Callback<(), Result<Snapshot, String>>,
    render: Callback<RenderRequest, LocalResult<Rendered>>,
    on_saved: Callback<String>,
) -> impl IntoView {
    let opening = StoredValue::new(snapshot.run(()).and_then(Capture::new));
    let prepared = RwSignal::new(None::<Arc<Prepared>>);
    let busy = RwSignal::new(false);
    let generation = RwSignal::new(0u64);
    let status = RwSignal::new(String::new());
    let alive = Arc::new(AtomicBool::new(true));
    let cleanup = alive.clone();
    on_cleanup(move || {
        cleanup.store(false, Ordering::Relaxed);
    });
    let current = Callback::new(move |()| -> Result<Capture, String> {
        let captured = Capture::new(snapshot.run(())?)?;
        opening.get_value()?.matches(captured.snapshot.clone())?;
        Ok(captured)
    });
    let prepare = {
        let alive = alive.clone();
        move |_| {
            if busy.get_untracked() {
                return;
            }
            prepared.set(None);
            status.set(String::new());
            let capture = match current.run(()) {
                Ok(capture) => capture,
                Err(error) => {
                    status.set(error);
                    return;
                }
            };
            let token = generation.get_untracked().wrapping_add(1);
            generation.set(token);
            busy.set(true);
            let alive = alive.clone();
            spawn_local(async move {
                let reply = render.run(capture.request()).await;
                if !alive.load(Ordering::Relaxed) || generation.try_get_untracked() != Some(token) {
                    return;
                }
                busy.set(false);
                let result = (|| {
                    capture.matches(current.run(())?.snapshot)?;
                    capture.accept(reply?)
                })();
                match result {
                    Ok(value) => prepared.set(Some(Arc::new(value))),
                    Err(error) => status.set(error),
                }
            });
        }
    };
    let save = {
        let alive = alive.clone();
        move |_| {
            if busy.get_untracked() {
                return;
            }
            let Some(value) = prepared.get_untracked() else {
                return;
            };
            if let Err(error) = current
                .run(())
                .and_then(|now| value.capture.matches(now.snapshot))
            {
                prepared.set(None);
                status.set(error);
                return;
            }
            let artifact = match crate::app::prepare_download(
                &value.filename(),
                "image/svg+xml",
                value.rendered.svg.as_bytes(),
            ) {
                Ok(artifact) => artifact,
                Err(error) => {
                    status.set(error);
                    return;
                }
            };
            busy.set(true);
            let token = generation.get_untracked();
            let guard_capture = value.capture.clone();
            let saving = alive.clone();
            let completion = alive.clone();
            let description = format!(
                "Saved SVG map of {} · document {} · event {}.",
                value.capture.snapshot.source.name,
                value.capture.snapshot.event.document_id,
                value.capture.snapshot.event.event_id
            );
            artifact.complete_guarded(
                move || {
                    if !saving.load(Ordering::Relaxed)
                        || generation.try_get_untracked() != Some(token)
                    {
                        return Err("Map export dialog changed before save admission.".into());
                    }
                    guard_capture.matches(current.run(())?.snapshot)
                },
                move |result| {
                    // Acknowledged bytes belong to the captured source. A later
                    // source change cannot turn real native success into failure.
                    if result.is_ok() {
                        on_saved.run(description.clone());
                    }
                    if !completion.load(Ordering::Relaxed) || busy.try_get_untracked().is_none() {
                        return;
                    }
                    busy.set(false);
                    match result {
                        Ok(()) => {
                            status.set(description);
                            prepared.set(None);
                        }
                        Err(error) => status.set(error),
                    }
                },
            );
        }
    };
    view! {<section class="map-export-dialog">
        <h3>"Export map as SVG"</h3>
        <p>"Export a standalone vector reference map of the whole document. Annotations and supplied restriction coordinates describe the source; they do not verify experimental cleavage or cloning outcomes."</p>
        <p class="helper">"The portable renderer uses its own canonical orientation and feature/label layout. Viewport rotation, zoom, selected-range shading, hidden-feature filters, search/PCR overlays and hovered-site highlights are not reproduced. Displayed restriction-site filters are preserved when supplied. The source, event and current view must remain unchanged until save admission."</p>
        <p class="helper">"Rendering is bounded to 512 features, 4096 feature spans, 32 feature labels and 64 printed restriction labels. Additional annotations may be omitted or abbreviated; review renderer warnings below. The native workspace retains full annotations."</p>
        <button class="button secondary" disabled=move||busy.get() on:click=prepare>"Prepare SVG"</button>
        {move||prepared.get().map(|value|view!{<div>
            <p>{format!("{} · {} bytes · source {}",value.filename(),value.rendered.svg.len(),value.capture.snapshot.event.revision)}</p>
            <ul>{value.rendered.warnings.iter().cloned().map(|warning|view!{<li>{warning}</li>}).collect_view()}</ul>
            <button class="button primary" disabled=move||busy.get() on:click=save.clone()>"Save SVG"</button>
        </div>})}
        <p role="status" aria-live="polite">{move||status.get()}</p>
    </section>}
}
