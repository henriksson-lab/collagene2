//! Export the exact currently reviewed graph; scientific work remains in a worker.
use crate::platform::worker::CachedWorker;
use collagene_core::{
    formats::oligos::Delimiter,
    planning::plan::DesignPlan,
    protocol::{
        request::Request,
        schema::{AnalysisJob, AnalysisResult},
    },
    verification,
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
#[derive(Clone, Copy)]
enum Format {
    Csv,
    Tsv,
    Html,
    Pdf,
}
#[component]
pub fn AssemblyOutputs(
    capture: Callback<(), Result<(DesignPlan, Molecule), String>>,
    identity: Memo<(String, String)>,
    epoch: RwSignal<u64>,
) -> impl IntoView {
    let busy = RwSignal::new(false);
    let status = RwSignal::new(String::new());
    let generation = RwSignal::new(0u32);
    let worker = StoredValue::new_local(None::<Rc<CachedWorker>>);
    let alive = Arc::new(AtomicBool::new(true));
    let cleanup = alive.clone();
    on_cleanup(move || {
        cleanup.store(false, Ordering::Relaxed);
        if let Some(w) = worker.get_value() {
            w.close();
        }
        worker.set_value(None);
    });
    Effect::new(move |_| {
        let _ = (identity.get(), epoch.get());
        generation.update(|g| *g = g.wrapping_add(1));
        if let Some(w) = worker.get_value() {
            w.close();
        }
        worker.set_value(None);
        busy.set(false);
        status.set(String::new());
    });
    let export = Callback::new(move |format: Format| {
        if busy.get_untracked() {
            return;
        }
        let bound = identity.get_untracked();
        let bound_epoch = epoch.get_untracked();
        let (plan, target) = match capture.run(()) {
            Ok(v) => v,
            Err(error) => {
                status.set(error);
                return;
            }
        };
        if collagene_core::planning::plan::DesignPlan::revision(&plan) != bound.0
            || target.revision() != bound.1
        {
            status.set("The reviewed plan or target changed; reconstruct again.".into());
            return;
        }
        let options = collagene_core::planning::oligos::Options::default();
        let (job, revision) = match format {
            Format::Csv | Format::Tsv => {
                let delimiter = if matches!(format, Format::Csv) {
                    Delimiter::Comma
                } else {
                    Delimiter::Tab
                };
                let revision = match collagene_core::planning::oligos::export_revision(
                    &plan, &target, &options, delimiter,
                ) {
                    Ok(v) => v,
                    Err(error) => {
                        status.set(error.to_string());
                        return;
                    }
                };
                (
                    AnalysisJob::ExportPlanOligos {
                        plan,
                        target,
                        options,
                        delimiter,
                    },
                    revision,
                )
            }
            Format::Html => (
                AnalysisJob::ExportProtocolHtml { plan, target },
                bound.1.clone(),
            ),
            Format::Pdf => (
                AnalysisJob::ExportProtocolPdf { plan, target },
                bound.1.clone(),
            ),
        };
        let w = match CachedWorker::new() {
            Ok(v) => Rc::new(v),
            Err(error) => {
                status.set(error);
                return;
            }
        };
        if let Some(previous) = worker.get_value() {
            previous.close();
        }
        worker.set_value(Some(w.clone()));
        let token = generation.get_untracked().wrapping_add(1);
        generation.set(token);
        busy.set(true);
        status.set("Preparing the reviewed design export…".into());
        let alive = alive.clone();
        spawn_local(async move {
            let result = w
                .ask(Request {
                    id: token,
                    revision: revision.clone(),
                    job,
                })
                .await;
            w.close();
            if !alive.load(Ordering::Relaxed)
                || generation.get_untracked() != token
                || epoch.get_untracked() != bound_epoch
                || identity.get_untracked() != bound
            {
                return;
            }
            worker.set_value(None);
            let export_guard = || -> Result<(), String> {
                if !alive.load(Ordering::Relaxed)
                    || generation.get_untracked() != token
                    || epoch.get_untracked() != bound_epoch
                    || identity.get_untracked() != bound
                {
                    return Err(
                        "The reviewed design changed before native export acceptance.".into(),
                    );
                }
                let (plan, target) = capture.run(())?;
                if collagene_core::planning::plan::DesignPlan::revision(&plan) != bound.0
                    || target.revision() != bound.1
                {
                    return Err("The source or target changed before export acceptance.".into());
                }
                Ok(())
            };
            let outcome: Result<String, String> = async {
                let response = result?;
                // Recheck host-owned source/target event guards after the worker,
                // including switches between equal-content documents.
                let (current_plan, current_target) = capture.run(())?;
                if collagene_core::planning::plan::DesignPlan::revision(&current_plan) != bound.0
                    || current_target.revision() != bound.1
                {
                    return Err("The reviewed plan or target changed during export.".into());
                }
                if response.id != token || response.revision != revision {
                    return Err("Export response does not match the reviewed design.".into());
                }
                match (format, response.result) {
                    (Format::Csv | Format::Tsv, AnalysisResult::PlanOligosExported(result)) => {
                        if result.report.plan_revision != bound.0
                            || result.report.target_revision != bound.1
                        {
                            return Err("Oligo list belongs to another design.".into());
                        }
                        if !result.report.complete {
                            return Err(format!(
                                "Order sheet needs review: {}",
                                result
                                    .report
                                    .unavailable
                                    .iter()
                                    .take(12)
                                    .cloned()
                                    .collect::<Vec<_>>()
                                    .join(" ")
                            ));
                        }
                        let (name, mime) = if matches!(format, Format::Csv) {
                            ("cloning-primers.csv", "text/csv")
                        } else {
                            ("cloning-primers.tsv", "text/tab-separated-values")
                        };
                        crate::io::download::download_guarded(
                            name,
                            mime,
                            result.export.text.as_bytes(),
                            &export_guard,
                        )
                        .await?;
                        Ok(format!(
                            "Exported {} oligos. {}",
                            result.report.oligos.records.len(),
                            result
                                .report
                                .warnings
                                .iter()
                                .chain(&result.export.warnings)
                                .take(8)
                                .cloned()
                                .collect::<Vec<_>>()
                                .join(" ")
                        ))
                    }
                    (Format::Html, AnalysisResult::ProtocolHtml(html)) => {
                        crate::io::download::download_guarded(
                            "cloning-protocol.html",
                            "text/html",
                            html.as_bytes(),
                            &export_guard,
                        )
                        .await?;
                        Ok("Exported HTML steps from this plan and target.".into())
                    }
                    (Format::Pdf, AnalysisResult::ProtocolPdf(pdf)) => {
                        if pdf.plan_revision != bound.0 || pdf.target_revision != bound.1 {
                            return Err("PDF belongs to another design.".into());
                        }
                        crate::io::download::download_guarded(
                            "cloning-protocol.pdf",
                            "application/pdf",
                            &pdf.bytes,
                            &export_guard,
                        )
                        .await?;
                        Ok(format!(
                            "Exported {} PDF pages · {:?}.",
                            pdf.pages, pdf.status
                        ))
                    }
                    (_, AnalysisResult::Error(error)) => Err(error),
                    _ => Err("Unexpected design export response.".into()),
                }
            }
            .await;
            if !alive.load(Ordering::Relaxed)
                || generation.get_untracked() != token
                || epoch.get_untracked() != bound_epoch
            {
                return;
            }
            busy.set(false);
            status.set(outcome.unwrap_or_else(|e| e));
        });
    });
    view! {<section class="assembly-outputs"><h4>"Design outputs"</h4><div class="assembly-output-actions"><button class="button secondary" disabled=move||busy.get() on:click=move|_|export.run(Format::Csv)>"Primer CSV"</button><button class="button secondary" disabled=move||busy.get() on:click=move|_|export.run(Format::Tsv)>"Primer TSV"</button><button class="button secondary" disabled=move||busy.get() on:click=move|_|export.run(Format::Html)>"Steps HTML"</button><button class="button secondary" disabled=move||busy.get() on:click=move|_|export.run(Format::Pdf)>"Steps PDF"</button><button class="text-tool" disabled=move||!busy.get() on:click=move|_|{generation.update(|g|*g=g.wrapping_add(1));if let Some(w)=worker.get_value(){w.close();}worker.set_value(None);busy.set(false);status.set("Export view cancelled; a native write already accepted may still finish.".into());}>"Cancel export"</button></div><p class="helper">"Order sheets list declared oligos, including 5′ tails. Reaction treatments do not become supplier modification requests."</p><p role="status">{move||status.get()}</p></section>}
}
