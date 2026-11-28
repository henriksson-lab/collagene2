//! Session-local annotated reference clipboard and explicit paste review.
//! Prepared material and physical cloning success are never inferred here.
use crate::platform::worker::CachedWorker;
use collagene_core::{
    material::lineage::EventBinding,
    protocol::{
        request::Request,
        schema::{AnalysisJob, AnalysisResult},
    },
    sequence::reference_transfer::{self, Fragment, Report, RequestPaste},
};
use collagene_model::Molecule;
use leptos::prelude::*;
use std::{
    cell::Cell,
    rc::Rc,
    sync::{
        Arc,
        atomic::{AtomicBool, Ordering},
    },
};
use wasm_bindgen_futures::spawn_local;

#[derive(Clone)]
pub struct Copied {
    pub source: EventBinding,
    pub source_molecule: Arc<Molecule>,
    pub fragment: Arc<Fragment>,
    pub epoch: u64,
}

#[derive(Clone)]
pub struct PasteOpen {
    pub copied: Arc<Copied>,
    pub destination: EventBinding,
    pub destination_molecule: Arc<Molecule>,
    pub start: usize,
    pub end: usize,
    pub epoch: u64,
}

#[derive(Clone)]
pub struct Apply {
    pub opened: PasteOpen,
    pub report: Report,
}

impl PasteOpen {
    fn request(&self) -> RequestPaste {
        RequestPaste {
            source: self.copied.source_molecule.as_ref().clone(),
            destination: self.destination_molecule.as_ref().clone(),
            selection: self.copied.fragment.selection.clone(),
            fragment_revision: self.copied.fragment.context_revision.clone(),
            start: self.start,
            end: self.end,
        }
    }
}

fn location_label(location: &collagene_model::Location) -> String {
    use collagene_model::Location;
    match location {
        Location::Range {
            start,
            end,
            before,
            after,
        } => format!(
            "{}{}..{}{}",
            if *before { "<" } else { "" },
            start,
            if *after { ">" } else { "" },
            end
        ),
        Location::Complement(inner) => match inner.as_ref() {
            Location::Range {
                start,
                end,
                before,
                after,
            } => format!(
                "complement({}{}..{}{})",
                if *before { "<" } else { "" },
                start,
                if *after { ">" } else { "" },
                end
            ),
            _ => "complemented compound location; full structure in JSON".into(),
        },
        Location::Join(parts) => format!("join · {} parts; full structure in JSON", parts.len()),
        Location::Order(parts) => format!("order · {} parts; full structure in JSON", parts.len()),
        Location::Unresolved(_) => "unresolved location; original evidence retained".into(),
    }
}
fn feature_label(feature: &collagene_model::Feature) -> String {
    format!(
        "{} · {} · {}{}",
        feature.label().chars().take(200).collect::<String>(),
        feature.kind.chars().take(64).collect::<String>(),
        location_label(&feature.location),
        if feature.needs_review {
            " · review required"
        } else {
            ""
        }
    )
}

#[component]
pub fn ReferencePastePanel(
    opened: PasteOpen,
    guard: Callback<PasteOpen, Result<(), String>>,
    on_apply: Callback<Apply, Result<(), String>>,
    on_cancel: Callback<()>,
) -> impl IntoView {
    let opened = StoredValue::new(opened);
    let report = RwSignal::new(None::<Arc<Report>>);
    let status = RwSignal::new(String::new());
    let busy = RwSignal::new(false);
    let generation = RwSignal::new(0u32);
    let sequence_page = RwSignal::new(0usize);
    let features_page = RwSignal::new(0usize);
    let warnings_page = RwSignal::new(0usize);
    let unmapped_page = RwSignal::new(0usize);
    let worker = StoredValue::new_local(None::<Rc<CachedWorker>>);
    let alive = Arc::new(AtomicBool::new(true));
    let cleanup_alive = alive.clone();
    on_cleanup(move || {
        cleanup_alive.store(false, Ordering::Relaxed);
        if let Some(worker) = worker.get_value() {
            worker.close();
        }
    });
    let cancel = Callback::new(move |()| {
        generation.update(|value| *value = value.wrapping_add(1));
        if let Some(worker) = worker.get_value() {
            worker.close();
        }
        worker.set_value(None);
        busy.set(false);
        on_cancel.run(());
    });
    let run = {
        let alive = alive.clone();
        Callback::new(move |apply: bool| {
            if busy.get_untracked() {
                return;
            }
            let captured = opened.get_value();
            let prepared = (|| -> Result<_, String> {
                guard.run(captured.clone())?;
                let request = captured.request();
                let revision =
                    reference_transfer::paste_revision(&request).map_err(|e| e.to_string())?;
                let worker = Rc::new(CachedWorker::new()?);
                Ok((request, revision, worker))
            })();
            let (request, revision, w) = match prepared {
                Ok(value) => value,
                Err(error) => {
                    status.set(error);
                    return;
                }
            };
            let token = generation.get_untracked().wrapping_add(1);
            generation.set(token);
            busy.set(true);
            worker.set_value(Some(w.clone()));
            status.set(
                if apply {
                    "Rebuilding the reference edit before applying…"
                } else {
                    "Preparing annotation-aware reference preview…"
                }
                .into(),
            );
            let alive = alive.clone();
            spawn_local(async move {
                let result = w
                    .ask(Request {
                        id: token,
                        revision: revision.clone(),
                        job: AnalysisJob::PasteReferenceFragment { request },
                    })
                    .await;
                w.close();
                if !alive.load(Ordering::Relaxed) || generation.get_untracked() != token {
                    return;
                }
                worker.set_value(None);
                busy.set(false);
                let outcome = (|| -> Result<Report, String> {
                    guard.run(captured.clone())?;
                    let response = result?;
                    if response.id != token || response.revision != revision {
                        return Err("The reference-paste response belongs to another edit.".into());
                    }
                    let returned = match response.result {
                        AnalysisResult::ReferencePaste(returned) => returned,
                        AnalysisResult::Error(error) => return Err(error),
                        _ => return Err("Unexpected reference-paste response.".into()),
                    };
                    returned
                        .validate_binding(&revision)
                        .map_err(|e| e.to_string())?;
                    if returned.source_revision != captured.copied.source.revision
                        || returned.destination_revision != captured.destination.revision
                        || returned.fragment_revision != captured.copied.fragment.context_revision
                        || returned.selection != captured.copied.fragment.selection
                        || returned.start != captured.start
                        || returned.end != captured.end
                    {
                        return Err("The reference-paste evidence differs from the selected versions or interval.".into());
                    }
                    Ok(returned)
                })();
                match outcome {
                    Ok(returned) if apply => {
                        let outcome = on_apply.run(Apply {
                            opened: captured,
                            report: returned,
                        });
                        if !alive.load(Ordering::Relaxed) {
                            return;
                        }
                        if let Err(error) = outcome {
                            status.set(error);
                        }
                    }
                    Ok(returned) => {
                        report.set(Some(Arc::new(returned)));
                        status.set("Reference preview ready. Review retained, changed and unmapped annotations before applying.".into());
                    }
                    Err(error) => status.set(error),
                }
            });
        })
    };
    let started = StoredValue::new_local(Cell::new(false));
    Effect::new(move |_| {
        if !started.with_value(|started| started.replace(true)) {
            run.run(false);
        }
    });
    let download_preview = move |_| {
        let Some(shown) = report.get_untracked() else {
            return;
        };
        let result = (|| {
            guard.run(opened.get_value())?;
            let bytes = serde_json::to_vec(shown.as_ref()).map_err(|e| e.to_string())?;
            if bytes.len() > collagene_core::formats::sequence::MAX_IMPORT_BYTES {
                return Err("Reference preview exceeds the JSON export limit.".into());
            }
            crate::io::download::ExportArtifact::prepare(
                "annotated-reference-paste.json",
                "application/json",
                &bytes,
            )
        })();
        match result {
            Ok(artifact) => {
                let check = move || {
                    if status.try_get_untracked().is_none() {
                        return Err("Export review closed.".into());
                    }
                    if !report.with_untracked(|current| {
                        current.as_ref().is_some_and(|r| Arc::ptr_eq(r, &shown))
                    }) {
                        return Err("Reference preview changed.".into());
                    }
                    guard.run(opened.get_value())
                };
                let completion = check.clone();
                artifact.complete_guarded(check, move |result| {
                    if completion().is_err() {
                        return;
                    }
                    status.set(match result {
                        Ok(()) => "Reference preview saved.".into(),
                        Err(error) => error,
                    });
                });
            }
            Err(error) => status.set(error),
        }
    };
    view! {
        <section class="reference-paste-review">
            <h3>"Paste reference with annotations"</h3>
            <p class="helper">"This edits a sequence reference. It does not reconstruct digestion, ligation, synthesis or physical material. Source documents remain unchanged; Undo restores the destination's prior version."</p>
            <p style="overflow-wrap:anywhere">{opened.with_value(|opened|format!("{} → {} · replace zero-based interval [{}, {})",opened.copied.source_molecule.name,opened.destination_molecule.name,opened.start,opened.end))}</p>
            <details style="overflow-wrap:anywhere"><summary>"Exact source and destination versions"</summary>
                <p>{opened.with_value(|opened|format!("Source {} · event {} · revision {}",opened.copied.source.document_id,opened.copied.source.event_id,opened.copied.source.revision))}</p>
                <p>{opened.with_value(|opened|format!("Destination {} · event {} · revision {}",opened.destination.document_id,opened.destination.event_id,opened.destination.revision))}</p>
            </details>
            <p class="helper" role="status">{move||status.get()}</p>
            {move||report.get().map(|shown|{
                let total=shown.inserted_end-shown.inserted_start;
                let sequence_pages=total.div_ceil(600).max(1);
                let page=sequence_page.get().min(sequence_pages-1);
                let start=(page*600).min(total);let end=(start+600).min(total);
                let letters=shown.product.sequence.as_str()[shown.inserted_start+start..shown.inserted_start+end].to_string();
                let features_count=shown.product.features.len();let feature_pages=features_count.div_ceil(25).max(1);let feature_page=features_page.get().min(feature_pages-1);
                let warning_pages=shown.warnings.len().div_ceil(20).max(1);let warning_page=warnings_page.get().min(warning_pages-1);
                let evidence_pages=shown.unmapped_features.len().div_ceil(20).max(1);let evidence_page=unmapped_page.get().min(evidence_pages-1);
                view!{
                    <p>{format!("Result {} bases · {} total annotations · {} copied annotations · {} unmapped source annotations · {} fully omitted source annotations",shown.product.sequence.len(),features_count,shown.copied_features,shown.unmapped_features.len(),shown.dropped_features)}</p>
                    <p class="helper">{format!("Inserted destination interval [{}, {}). Cropped, interrupted or unsupported annotations remain review evidence; source IDs/names do not prove function.",shown.inserted_start,shown.inserted_end)}</p>
                    <details open><summary>"Inserted reference bases · all pages"</summary><pre style="white-space:pre-wrap;overflow-wrap:anywhere">{letters}</pre><p>{format!("Fragment bases {}–{} of {total} · page {} of {sequence_pages}",start+1,end,page+1)}</p>
                        <div class="method-row-actions"><button class="text-tool" disabled=move||busy.get()||page==0 on:click=move|_|sequence_page.set(page.saturating_sub(1))>"Previous bases"</button><button class="text-tool" disabled=move||{busy.get()||page+1>=sequence_pages} on:click=move|_|sequence_page.set((page+1).min(sequence_pages-1))>"Next bases"</button></div>
                    </details>
                    <details><summary>"Result annotations"</summary><ul>{shown.product.features.iter().skip(feature_page*25).take(25).map(|feature|view!{<li>{feature_label(feature)}</li>}).collect_view()}</ul><p>{format!("Annotation page {} of {feature_pages}",feature_page+1)}</p>
                        <div class="method-row-actions"><button class="text-tool" disabled=move||busy.get()||feature_page==0 on:click=move|_|features_page.set(feature_page.saturating_sub(1))>"Previous annotations"</button><button class="text-tool" disabled=move||{busy.get()||feature_page+1>=feature_pages} on:click=move|_|features_page.set((feature_page+1).min(feature_pages-1))>"Next annotations"</button></div>
                    </details>
                    <div>{shown.warnings.iter().skip(warning_page*20).take(20).cloned().map(|warning|view!{<p class="helper">{warning}</p>}).collect_view()}</div><p>{format!("Warning page {} of {warning_pages}",warning_page+1)}</p>
                    <div class="method-row-actions"><button class="text-tool" disabled=move||busy.get()||warning_page==0 on:click=move|_|warnings_page.set(warning_page.saturating_sub(1))>"Previous warnings"</button><button class="text-tool" disabled=move||{busy.get()||warning_page+1>=warning_pages} on:click=move|_|warnings_page.set((warning_page+1).min(warning_pages-1))>"Next warnings"</button></div>
                    <details><summary>"Unmapped annotation evidence"</summary><ul>{shown.unmapped_features.iter().skip(evidence_page*20).take(20).map(|feature|view!{<li>{feature_label(feature)}</li>}).collect_view()}</ul><p>{format!("Evidence page {} of {evidence_pages}",evidence_page+1)}</p>
                        <div class="method-row-actions"><button class="text-tool" disabled=move||busy.get()||evidence_page==0 on:click=move|_|unmapped_page.set(evidence_page.saturating_sub(1))>"Previous evidence"</button><button class="text-tool" disabled=move||{busy.get()||evidence_page+1>=evidence_pages} on:click=move|_|unmapped_page.set((evidence_page+1).min(evidence_pages-1))>"Next evidence"</button></div>
                    </details>
                }
            })}
            <p class="helper">"Location labels summarize compound annotations; full locations, qualifiers and source evidence remain in the preview JSON and retained native reference."</p>
            <button class="button secondary" disabled=move||busy.get()||report.get().is_none() on:click=download_preview>"Download complete reference preview"</button>
            <div class="method-row-actions"><button class="button primary" disabled=move||busy.get()||report.get().is_none() on:click=move|_|run.run(true)>"Apply reference paste"</button><button class="button secondary" on:click=move|_|cancel.run(())>"Cancel"</button></div>
        </section>
    }.into_any()
}
