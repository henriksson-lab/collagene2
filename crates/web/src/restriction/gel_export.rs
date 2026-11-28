//! One immutable displayed gel report is prepared in a worker and saved under a live-pane guard.
use crate::{io::download::ExportArtifact, platform::worker};
use collagene_core::{
    output::{
        gel::GelReport,
        gel_export,
        raster::{Format, PngOptions},
    },
    protocol::{
        request::Request,
        schema::{AnalysisJob, AnalysisResult},
    },
};
use leptos::prelude::*;
use std::sync::{
    Arc,
    atomic::{AtomicBool, Ordering},
};
use wasm_bindgen_futures::spawn_local;

#[component]
pub fn GelPngExport(
    report: GelReport,
    is_current: Callback<(), bool>,
    on_notice: Callback<String>,
) -> impl IntoView {
    let source = StoredValue::new(report);
    let width = RwSignal::new(1600u32);
    let busy = RwSignal::new(false);
    let generation = RwSignal::new(0u32);
    let alive = Arc::new(AtomicBool::new(true));
    let disposed = alive.clone();
    on_cleanup(move || disposed.store(false, Ordering::Relaxed));
    let prepare = {
        let alive = alive.clone();
        move |_| {
            if busy.get_untracked() || !is_current.run(()) {
                return;
            }
            let request = match gel_export::Request::new(
                source.get_value(),
                Format::Png(PngOptions {
                    width: width.get_untracked(),
                }),
            ) {
                Ok(request) => request,
                Err(error) => {
                    on_notice.run(error.to_string());
                    return;
                }
            };
            let revision = match request.revision() {
                Ok(revision) => revision,
                Err(error) => {
                    on_notice.run(error.to_string());
                    return;
                }
            };
            let token = generation.get_untracked().wrapping_add(1);
            generation.set(token);
            busy.set(true);
            on_notice.run("Preparing the displayed predicted gel as PNG in the background…".into());
            let alive = alive.clone();
            spawn_local(async move {
                let response = worker::analyze(Request {
                    id: token,
                    revision: revision.clone(),
                    job: AnalysisJob::RenderGelImage { request },
                })
                .await;
                if !alive.load(Ordering::Relaxed) || generation.try_get_untracked() != Some(token) {
                    return;
                }
                busy.set(false);
                if !is_current.run(()) {
                    on_notice
                        .run("The displayed gel changed; prepare its current report again.".into());
                    return;
                }
                let rendered = match response {
                    Ok(response) if response.id == token && response.revision == revision => {
                        match response.result {
                            AnalysisResult::GelImageExported(rendered)
                                if matches!(
                                    &rendered.image,
                                    collagene_core::output::raster::Image::Png(_)
                                ) =>
                            {
                                rendered
                            }
                            AnalysisResult::Error(error) => {
                                on_notice.run(error);
                                return;
                            }
                            _ => {
                                on_notice.run(
                                    "Unexpected gel image response; no file was saved.".into(),
                                );
                                return;
                            }
                        }
                    }
                    Ok(_) => {
                        on_notice.run("Gel image response belongs to another request.".into());
                        return;
                    }
                    Err(error) => {
                        on_notice.run(error);
                        return;
                    }
                };
                if let Err(error) = rendered.validate_binding(&revision) {
                    on_notice.run(error.to_string());
                    return;
                }
                let artifact = match ExportArtifact::prepare(
                    "predicted-gel.png",
                    "image/png",
                    rendered.image.bytes(),
                ) {
                    Ok(artifact) => artifact,
                    Err(error) => {
                        on_notice.run(error);
                        return;
                    }
                };
                let notice = format!("Saved predicted gel PNG. {}", rendered.warnings.join(" "));
                busy.set(true);
                let saving = alive.clone();
                let complete = alive.clone();
                artifact.complete_guarded(
                    move || {
                        if !saving.load(Ordering::Relaxed)
                            || generation.try_get_untracked() != Some(token)
                            || !is_current.run(())
                        {
                            return Err("Gel view changed before save admission.".into());
                        }
                        Ok(())
                    },
                    move |result| {
                        if !complete.load(Ordering::Relaxed) || busy.try_get_untracked().is_none() {
                            return;
                        }
                        busy.set(false);
                        on_notice.run(match result {
                            Ok(()) => notice,
                            Err(error) => error,
                        });
                    },
                );
            });
        }
    };
    view! {<div class="gel-png-export">
        <label>"PNG width (pixels) "<input type="number" min="16" max="4096" prop:value=move||width.get() disabled=move||busy.get() on:input=move|event|width.set(event_target_value(&event).parse().unwrap_or(0))/></label>
        <button class="button secondary full" disabled=move||busy.get() on:click=prepare>"Download gel PNG"</button>
        <p class="helper">"Static predicted image with the bundled font. Height preserves aspect ratio; at most 4 million pixels. Keep the JSON report for complete source and calibration evidence."</p>
    </div>}
}
