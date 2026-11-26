//! Explicit annotated-protein exchange, isolated from cloning/material authority.
use collagene_core::{
    Molecule,
    formats::ExportReport,
    protein_formats::{self, ExportOptions, InsdcDate, Profile},
};
use leptos::prelude::*;
use std::sync::Arc;
use wasm_bindgen_futures::spawn_local;
fn utc_date() -> String {
    let now = js_sys::Date::new_0();
    let months = [
        "JAN", "FEB", "MAR", "APR", "MAY", "JUN", "JUL", "AUG", "SEP", "OCT", "NOV", "DEC",
    ];
    format!(
        "{:02}-{}-{:04}",
        now.get_utc_date(),
        months.get(now.get_utc_month() as usize).unwrap_or(&"JAN"),
        now.get_utc_full_year()
    )
}
#[component]
pub fn ProteinExportDialog(
    source: Molecule,
    profile: Profile,
    guard: Callback<(), Result<(), String>>,
    on_saved: Callback<String>,
) -> impl IntoView {
    let original_accession = source
        .metadata
        .get("protein:accessions")
        .and_then(|v| v.split(';').next())
        .unwrap_or("")
        .to_string();
    let accession = RwSignal::new(original_accession.clone());
    let original_entry = source
        .metadata
        .get("protein:entry_name")
        .cloned()
        .unwrap_or_else(|| source.name.clone());
    let entry = RwSignal::new(original_entry.clone());
    let date = RwSignal::new(utc_date());
    let prepared = RwSignal::new(None::<Arc<ExportReport>>);
    let error = RwSignal::new(String::new());
    let busy = RwSignal::new(false);
    let saved = RwSignal::new(false);
    let source = Arc::new(source);
    let filename = move || {
        format!(
            "{}.{}",
            entry.get_untracked(),
            if profile == Profile::GenPept {
                "gp"
            } else {
                "uniprot"
            }
        )
    };
    let prepare = {
        let source = source.clone();
        move |_| {
            prepared.set(None);
            saved.set(false);
            error.set(String::new());
            if let Err(e) = guard.run(()) {
                error.set(e);
                return;
            }
            let options = ExportOptions {
                locus_date: if profile == Profile::GenPept {
                    match InsdcDate::new(&date.get_untracked()) {
                        Ok(date) => Some(date),
                        Err(e) => {
                            error.set(e.to_string());
                            return;
                        }
                    }
                } else {
                    None
                },
                entry_name: (entry.get_untracked() != original_entry)
                    .then(|| entry.get_untracked()),
                uniprot_accession: if profile == Profile::UniProt
                    && accession.get_untracked() != original_accession
                {
                    Some(accession.get_untracked())
                } else {
                    None
                },
            };
            match protein_formats::export_with_options(
                std::slice::from_ref(source.as_ref()),
                profile,
                &options,
            ) {
                Ok(report) => prepared.set(Some(Arc::new(report))),
                Err(e) => error.set(e.to_string()),
            }
        }
    };
    let save = move |_| {
        let Some(report) = prepared.get_untracked() else {
            return;
        };
        if busy.get_untracked() {
            return;
        }
        busy.set(true);
        error.set(String::new());
        let name = filename();
        spawn_local(async move {
            let result = crate::app::download_guarded(
                &name,
                "text/plain",
                report.text.as_bytes(),
                move || guard.run(()),
            )
            .await;
            if busy.try_get_untracked().is_some() {
                busy.set(false);
            }
            match result {
                Ok(()) => {
                    if saved.try_get_untracked().is_some() {
                        saved.set(true);
                    }
                    on_saved.run(format!("Saved {name}. {}", report.warnings.join(" ")));
                }
                Err(e) => {
                    if error.try_get_untracked().is_some() {
                        error.set(e);
                    }
                }
            }
        });
    };
    view! {
        <p class="helper">"Exports the captured protein sequence and annotations. Source identifiers are documentary; this creates no verified material or database review."</p>
        <label>"Entry identifier"<input prop:value=move||entry.get() disabled=move||busy.get() on:input=move|e|{entry.set(event_target_value(&e));prepared.set(None);saved.set(false);}/></label>
        <Show when=move||profile==Profile::GenPept><label>"Export date (UTC; DD-MMM-YYYY)"<input prop:value=move||date.get() disabled=move||busy.get() on:input=move|e|{date.set(event_target_value(&e));prepared.set(None);saved.set(false);}/></label></Show>
        <Show when=move||profile==Profile::UniProt><label>"UniProt accession"<input prop:value=move||accession.get() disabled=move||busy.get() on:input=move|e|{accession.set(event_target_value(&e));prepared.set(None);saved.set(false);}/></label><p class="helper">"Enter an explicit valid accession when absent. Collagene does not allocate a UniProt accession; generated exchange records are Unreviewed. FASTA/native export remains available."</p></Show>
        <button class="button secondary" disabled=move||busy.get() on:click=prepare>"Prepare export"</button>
        <Show when=move||!error.get().is_empty()><p class="dialog-error" role="alert">{move||error.get()}</p></Show>
        {move||prepared.get().map(|report|view!{<p>{format!("{} bytes of normalized annotated protein text.",report.text.len())}</p><details open><summary>"Export notes"</summary><ul>{report.warnings.iter().cloned().map(|warning|view!{<li>{warning}</li>}).collect_view()}</ul></details>})}
        <button class="button primary" disabled=move||busy.get()||prepared.get().is_none()||saved.get() on:click=save>{move||if busy.get(){"Saving…"}else if saved.get(){"Saved"}else{"Save file"}}</button>
    }
}
