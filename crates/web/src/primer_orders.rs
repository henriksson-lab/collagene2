//! Explicit selected-oligo exports using the standard table and collection format.
use collagene_core::{
    analysis::Primer,
    oligo_formats::{self, Delimiter, OligoTable, TableExportMode},
};
use leptos::prelude::*;

#[derive(Clone, Copy)]
enum Action {
    Csv,
    Tsv,
    Save,
}
#[component]
pub fn PrimerOrdersPanel(
    primers: Vec<Primer>,
    guard: Callback<(), Result<(), String>>,
    on_collection: Callback<OligoTable, Result<(), String>>,
) -> impl IntoView {
    let primers = StoredValue::new(primers);
    let status = RwSignal::new(String::new());
    let action = Callback::new(move |action: Action| {
        let result = (|| -> Result<String, String> {
            guard.run(())?;
            let table = primers
                .with_value(|primers| collagene_core::oligo_orders::from_primers(primers))
                .map_err(|e| e.to_string())?;
            guard.run(())?;
            match action {
                Action::Save => {
                    let evidence_only = table
                        .records
                        .iter()
                        .filter(|row| row.primer.is_none())
                        .count();
                    on_collection.run(table)?;
                    Ok(if evidence_only > 0 {
                        format!(
                            "Oligo collection retained. {evidence_only} rows remain order evidence outside the reusable-primer profile."
                        )
                    } else {
                        "Named oligo collection retained.".into()
                    })
                }
                Action::Csv | Action::Tsv => {
                    let (delimiter, filename, mime) = match action {
                        Action::Csv => (
                            Delimiter::Comma,
                            "selected-primers.csv",
                            "text/csv;charset=utf-8",
                        ),
                        _ => (
                            Delimiter::Tab,
                            "selected-primers.tsv",
                            "text/tab-separated-values;charset=utf-8",
                        ),
                    };
                    let report = oligo_formats::export_oligos(
                        &table,
                        delimiter,
                        TableExportMode::SpreadsheetSafe,
                    )
                    .map_err(|e| e.to_string())?;
                    guard.run(())?;
                    crate::app::download(filename, mime, report.text.as_bytes())?;
                    Ok(format!(
                        "Named oligo order reference exported. {}",
                        report.warnings.join(" ")
                    ))
                }
            }
        })();
        status.set(result.unwrap_or_else(|error| error));
    });
    view! {
        <section class="primer-order-outputs" aria-label="Selected oligo outputs">
            <h4>"Selected oligo outputs"</h4>
            {primers.with_value(|primers|primers.iter().map(|primer|view!{
                <p class="helper">{format!("{} · ordered {} bases · tail {} · annealing {}",primer.name,primer.ordered_sequence().len(),if primer.tail.is_empty(){"(none)"}else{&primer.tail},primer.annealing.as_str())}</p>
            }).collect_view())}
            <div class="primer-actions">
                <button class="button secondary" on:click=move|_|action.run(Action::Csv)>"Primer CSV"</button>
                <button class="button secondary" on:click=move|_|action.run(Action::Tsv)>"Primer TSV"</button>
                <button class="button secondary" on:click=move|_|action.run(Action::Save)>"Save oligo collection"</button>
            </div>
            <p class="helper">"Full ordered 5′→3′ bases and domain/tail splits are retained. Orders do not establish amplification, synthesis or material availability."</p>
            <p class="helper" role="status">{move||status.get()}</p>
        </section>
    }
}
