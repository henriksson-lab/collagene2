//! Reviewed sequence exchange. Native history and source evidence remain separate.
use collagene_core::{Molecule, gcg_sequence, snapgene_export};
use leptos::prelude::*;
use std::sync::Arc;
use wasm_bindgen_futures::spawn_local;
#[derive(Clone)]
struct Prepared {
    bytes: Vec<u8>,
    name: String,
    mime: String,
    warnings: Vec<String>,
}
#[component]
pub fn SequenceCompatExport(
    source: Molecule,
    snapgene: bool,
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
    let warnings = RwSignal::new(Vec::<String>::new());
    let prepared = RwSignal::new(None::<Arc<Prepared>>);
    let busy = RwSignal::new(false);
    let source = Arc::new(source);
    let prepare = {
        let source = source.clone();
        move |_| {
            prepared.set(None);
            error.set(String::new());
            warnings.set(vec![]);
            if let Err(e) = guard.run(()) {
                error.set(e);
                return;
            }
            let result: Result<Prepared, String> = (|| {
                if snapgene {
                    let review = snapgene_export::review(&source).map_err(|e| e.to_string())?;
                    let mut summary = review.warnings;
                    summary.extend(
                        review
                            .losses
                            .iter()
                            .map(|l| format!("{}: {}", l.path, l.message)),
                    );
                    warnings.set(summary);
                    if !review.losses.is_empty() && !allow.get_untracked() {
                        return Err("Review the listed omissions and acknowledge them before preparing this export.".into());
                    }
                    let report = snapgene_export::export_snapgene(
                        &source,
                        if allow.get_untracked() {
                            snapgene_export::Policy::AllowOmissions
                        } else {
                            snapgene_export::Policy::Strict
                        },
                    )
                    .map_err(|e| e.to_string())?;
                    let mut notices = report.warnings;
                    notices.extend(
                        report
                            .losses
                            .into_iter()
                            .map(|l| format!("{}: {}", l.path, l.message)),
                    );
                    Ok(Prepared {
                        bytes: report.bytes,
                        name: format!("{}{}", source.name, report.extension),
                        mime: report.mime,
                        warnings: notices,
                    })
                } else {
                    let report = gcg_sequence::export_gcg(
                        &source,
                        &gcg_sequence::ExportOptions {
                            record_id: Some(id.get_untracked()),
                            date_label: Some(date.get_untracked()),
                        },
                    )
                    .map_err(|e| e.to_string())?;
                    Ok(Prepared {
                        bytes: report.text.into_bytes(),
                        name: format!("{}.gcg", id.get_untracked()),
                        mime: "text/plain".into(),
                        warnings: report.warnings,
                    })
                }
            })();
            match result {
                Ok(value) => {
                    warnings.set(value.warnings.clone());
                    prepared.set(Some(Arc::new(value)));
                }
                Err(e) => error.set(e),
            }
        }
    };
    let save = move |_| {
        if busy.get_untracked() {
            return;
        }
        let Some(value) = prepared.get_untracked() else {
            return;
        };
        if !value.warnings.is_empty() && !allow.get_untracked() {
            error.set("Acknowledge the listed omissions before saving.".into());
            return;
        }
        if let Err(e) = guard.run(()) {
            error.set(e);
            return;
        }
        busy.set(true);
        error.set(String::new());
        spawn_local(async move {
            let result =
                crate::app::download_guarded(&value.name, &value.mime, &value.bytes, move || {
                    guard.run(())
                })
                .await;
            if busy.try_get_untracked().is_some() {
                busy.set(false);
            }
            match result {
                Ok(()) => {
                    if busy.try_get_untracked().is_some() && guard.run(()).is_ok() {
                        on_saved.run(format!(
                            "Saved {}. {}",
                            value.name,
                            value.warnings.join(" ")
                        ));
                    }
                }
                Err(e) => {
                    if error.try_get_untracked().is_some() {
                        error.set(e);
                    }
                }
            }
        });
    };
    view! {<div class="sequence-compat-export"><p class="helper">{if snapgene{"Selected SnapGene reference-sequence profile. Features are written only when their locations and qualifiers can be represented. The original file, cloning history and physical chemistry remain in the native document."}else{"GCG single-sequence exchange retains letters and a checksum. Features, circular topology, cloning history and original attachments require the native document."}}</p>
    <Show when=move||!snapgene><label>"Record ID"<input prop:value=move||id.get() disabled=move||busy.get() on:input=move|e|{id.set(event_target_value(&e));prepared.set(None);}/></label><label>"Export date (UTC)"<input prop:value=move||date.get() disabled=move||busy.get() on:input=move|e|{date.set(event_target_value(&e));prepared.set(None);}/></label></Show>
    <label><input type="checkbox" prop:checked=move||allow.get() disabled=move||busy.get() on:change=move|e|{allow.set(event_target_checked(&e));if snapgene{prepared.set(None);}}/>"I accept the listed omissions for this exchange file"</label>
    <button class="button secondary" disabled=move||busy.get() on:click=prepare>"Review export"</button>
    <ul>{move||warnings.get().into_iter().map(|w|view!{<li>{w}</li>}).collect_view()}</ul><p role="alert">{move||error.get()}</p>
    <button class="button" disabled=move||busy.get()||prepared.get().is_none()||(!warnings.get().is_empty()&&!allow.get()) on:click=save>"Save reviewed export"</button>
    </div>}
}
