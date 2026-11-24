//! Explicit oligo column mapping and retained collection/row selection.
use collagene_core::oligo_formats::{
    self, ColumnMapping, Delimiter, OligoImportOptions, OligoTable, TableExportMode,
};
use leptos::prelude::*;
use std::{collections::BTreeSet, sync::Arc};

#[component]
fn ColumnChoice(
    label: &'static str,
    value: RwSignal<String>,
    headers: RwSignal<Vec<String>>,
    on_change: Callback<()>,
) -> impl IntoView {
    view! {<label>{label}<select aria-label=label prop:value=move||value.get() on:change=move|e|{value.set(event_target_value(&e));on_change.run(());}><option value="">"Not mapped"</option>{move||headers.get().into_iter().enumerate().map(|(i,title)|view!{<option value=i.to_string()>{format!("{} · {}",i+1,title.chars().take(100).collect::<String>())}</option>}).collect_view()}</select></label>}
}

#[component]
pub fn OligoImport(
    on_import: Callback<(OligoTable, Vec<usize>), Result<(), String>>,
    #[prop(optional_no_strip)] existing: Option<OligoTable>,
    #[prop(optional_no_strip)] initial_text: Option<(String, Delimiter)>,
) -> impl IntoView {
    let mut initial_options = existing
        .as_ref()
        .map(|t| t.options.clone())
        .unwrap_or_default();
    let initial_source = existing.as_ref().map(|t| {
        t.original_bytes().and_then(|b| {
            String::from_utf8(b).map_err(|e| collagene_core::Error::Format(e.to_string()))
        })
    });
    let initial_error = initial_source
        .as_ref()
        .and_then(|r| r.as_ref().err())
        .map(ToString::to_string)
        .unwrap_or_default();
    if let Some((_, delimiter)) = &initial_text {
        initial_options.delimiter = *delimiter;
    }
    let source = RwSignal::new(
        initial_text
            .map(|(text, _)| text)
            .or_else(|| initial_source.and_then(Result::ok))
            .unwrap_or_default(),
    );
    let delimiter = RwSignal::new(initial_options.delimiter);
    let has_header = RwSignal::new(initial_options.has_header);
    let name = RwSignal::new(initial_options.columns.name.to_string());
    let val = |index: Option<usize>| index.map(|i| i.to_string()).unwrap_or_default();
    let ordered = RwSignal::new(val(initial_options.columns.ordered_sequence));
    let domain = RwSignal::new(val(initial_options.columns.annealing_sequence));
    let tail = RwSignal::new(val(initial_options.columns.tail));
    let modification = RwSignal::new(val(initial_options.columns.modification));
    let amount = RwSignal::new(val(initial_options.columns.amount));
    let unit = RwSignal::new(val(initial_options.columns.unit));
    let headers = RwSignal::new(
        existing
            .as_ref()
            .map(|t| t.headers.clone())
            .unwrap_or_default(),
    );
    let table = RwSignal::new(if initial_error.is_empty() {
        existing.map(Arc::new)
    } else {
        None
    });
    let table_state = table;
    let selected = RwSignal::new(BTreeSet::<usize>::new());
    let status = RwSignal::new(initial_error);
    let page = RwSignal::new(0usize);
    let reset = Callback::new(move |()| {
        table.set(None);
        selected.set(BTreeSet::new());
        page.set(0);
        status.set(String::new());
    });
    let read = move |_| {
        reset.run(());
        match oligo_formats::inspect_table(
            &source.get_untracked(),
            delimiter.get_untracked(),
            has_header.get_untracked(),
        ) {
            Ok(preview) => {
                let count = preview
                    .rows
                    .iter()
                    .map(Vec::len)
                    .max()
                    .unwrap_or(0)
                    .max(preview.headers.len());
                let mut labels = preview.headers;
                for i in labels.len()..count {
                    labels.push(format!("Column {} (no header)", i + 1));
                }
                headers.set(labels);
                status.set(format!(
                    "{} source rows. Review every mapped role, then interpret the table.",
                    preview.rows.len()
                ));
            }
            Err(e) => {
                headers.set(Vec::new());
                status.set(e.to_string());
            }
        }
    };
    let interpret = move |_| {
        let optional = |v: String| {
            if v.is_empty() {
                Ok(None)
            } else {
                v.parse::<usize>()
                    .map(Some)
                    .map_err(|_| "Invalid column index".to_string())
            }
        };
        let options = (|| {
            Ok::<_, String>(OligoImportOptions {
                delimiter: delimiter.get_untracked(),
                has_header: has_header.get_untracked(),
                columns: ColumnMapping {
                    name: name
                        .get_untracked()
                        .parse()
                        .map_err(|_| "Choose a name column".to_string())?,
                    ordered_sequence: optional(ordered.get_untracked())?,
                    annealing_sequence: optional(domain.get_untracked())?,
                    tail: optional(tail.get_untracked())?,
                    modification: optional(modification.get_untracked())?,
                    amount: optional(amount.get_untracked())?,
                    unit: optional(unit.get_untracked())?,
                },
            })
        })();
        match options.and_then(|o| {
            oligo_formats::import_oligos(&source.get_untracked(), &o).map_err(|e| e.to_string())
        }) {
            Ok(result) => {
                selected.set(BTreeSet::new());
                page.set(0);
                status.set(format!("{} source rows · {} convertible unmodified DNA oligos. Choose rows to open, or keep the collection alone.",result.records.len(),result.records.iter().filter(|r|r.primer.is_some()).count()));
                table.set(Some(Arc::new(result)));
            }
            Err(e) => status.set(e),
        }
    };
    let save = move |rows: Vec<usize>| {
        let Some(table) = table.get_untracked() else {
            return;
        };
        if let Err(e) = on_import.run((table.as_ref().clone(), rows)) {
            status.set(e);
        }
    };
    view! {
        <style>".oligo-import .oligo-grid{display:grid;grid-template-columns:1fr 1fr;gap:10px}.oligo-import .oligo-grid label{margin:0}.oligo-import .oligo-table-wrap{max-height:420px;overflow:auto}.oligo-import table{border-collapse:collapse;font-size:12px;min-width:620px;width:100%}.oligo-import th,.oligo-import td{text-align:left;vertical-align:top;padding:7px;border-bottom:1px solid #e1e8e0;overflow-wrap:anywhere;max-width:280px}.oligo-import table input[type=checkbox]{width:auto}.oligo-import .oligo-actions{display:flex;flex-wrap:wrap;gap:8px;margin:12px 0}.oligo-import .oligo-check{display:flex;align-items:center;gap:8px}.oligo-import .oligo-check input{width:auto}.oligo-import textarea{font:12px/1.5 monospace}.oligo-import .oligo-status{font-size:13px;line-height:1.6;overflow-wrap:anywhere}"</style>
        <div class="oligo-import"><p>"Map the source columns explicitly. Modified or unsupported oligos remain in the collection; they are never stripped into ordinary DNA. Quantities are source fields, not stock management."</p>
        <label>"CSV/TSV source"<textarea aria-label="Oligo table source" rows="7" prop:value=move||source.get() on:input=move|e|{source.set(event_target_value(&e));headers.set(Vec::new());reset.run(());}/></label>
        <div class="oligo-grid"><label>"Delimiter"<select aria-label="Oligo table delimiter" prop:value=move||if delimiter.get()==Delimiter::Tab{"tab"}else{"comma"} on:change=move|e|{delimiter.set(if event_target_value(&e)=="tab"{Delimiter::Tab}else{Delimiter::Comma});headers.set(Vec::new());reset.run(());}><option value="comma">"Comma (CSV)"</option><option value="tab">"Tab (TSV)"</option></select></label><label class="oligo-check"><input type="checkbox" aria-label="Oligo table has header" prop:checked=move||has_header.get() on:change=move|e|{has_header.set(event_target::<web_sys::HtmlInputElement>(&e).checked());headers.set(Vec::new());reset.run(());}/>"First record is a header"</label></div>
        <button class="button secondary" on:click=read>"Read oligo columns"</button>
        <Show when=move||!headers.with(|h|h.is_empty())><div class="oligo-grid"><ColumnChoice label="Oligo name column" value=name headers=headers on_change=reset/><ColumnChoice label="Oligo ordered sequence column" value=ordered headers=headers on_change=reset/><ColumnChoice label="Oligo annealing column" value=domain headers=headers on_change=reset/><ColumnChoice label="Oligo tail column" value=tail headers=headers on_change=reset/><ColumnChoice label="Oligo modification column" value=modification headers=headers on_change=reset/><ColumnChoice label="Oligo amount column" value=amount headers=headers on_change=reset/><ColumnChoice label="Oligo unit column" value=unit headers=headers on_change=reset/></div><button class="button primary" on:click=interpret>"Interpret oligo table"</button></Show>
        <p class="oligo-status" role="status">{move||status.get()}</p>
        {move||table.get().map(|table|{let collection=table.clone();let original=table.clone();let csv=table.clone();let tsv=table.clone();let count=table.records.len();let start=page.get()*100;let rows=table.records.iter().skip(start).take(100).cloned().collect::<Vec<_>>();view!{
            <section class="oligo-result"><h3>{format!("{} retained source rows",count)}</h3><ul>{table.warnings.iter().cloned().map(|w|view!{<li>{w}</li>}).collect_view()}</ul>
            <div class="oligo-actions"><button class="button secondary" on:click=move|_|{let Some(t)=table_state.get_untracked()else{return;};let start=page.get_untracked()*100;selected.update(|s|for r in t.records.iter().skip(start).take(100).filter(|r|r.primer.is_some()){s.insert(r.row_index);});}>"Select convertible rows on page"</button><button class="button secondary" on:click=move|_|selected.set(BTreeSet::new())>"Clear oligo selection"</button></div>
            <div class="oligo-table-wrap"><table><thead><tr><th>"Open"</th><th>"Source row / name"</th><th>"Ordered source / status"</th></tr></thead><tbody>{rows.into_iter().map(|row|{let index=row.row_index;let convertible=row.primer.is_some();let name=row.name;let sequence=row.raw_ordered_sequence.or(row.raw_annealing_sequence).unwrap_or_default();view!{<tr><td><input type="checkbox" aria-label=format!("Select oligo row {}",index+1) disabled=!convertible prop:checked=move||selected.with(|s|s.contains(&index)) on:change=move|e|{let checked=event_target::<web_sys::HtmlInputElement>(&e).checked();selected.update(|s|{if checked{s.insert(index);}else{s.remove(&index);}});}/></td><td>{format!("{} · {}",index+1,name.chars().take(160).collect::<String>())}</td><td><code>{sequence.chars().take(200).collect::<String>()}</code><p>{if convertible{"Unmodified DNA fields are consistent"}else{"Evidence only · unsupported notation/modification or inconsistent fields"}}</p><details><summary>"Row diagnostics"</summary><ul>{row.diagnostics.into_iter().map(|w|view!{<li>{w}</li>}).collect_view()}</ul></details></td></tr>}}).collect_view()}</tbody></table></div>
            <div class="oligo-actions"><button class="button secondary" disabled=move||page.get()==0 on:click=move|_|page.update(|p|*p=p.saturating_sub(1))>"Previous oligo rows"</button><button class="button secondary" disabled=move||{(page.get()+1)*100>=count} on:click=move|_|page.update(|p|*p+=1)>"Next oligo rows"</button><span>{move||format!("{} selected",selected.with(BTreeSet::len))}</span></div>
            <div class="oligo-actions"><button class="button primary" disabled=move||selected.with(BTreeSet::is_empty) on:click=move|_|save(selected.get_untracked().into_iter().collect())>"Keep collection and open selected primers"</button><button class="button secondary" on:click=move|_|save(Vec::new())>"Keep oligo collection only"</button></div>
            <div class="oligo-actions"><button class="button secondary" on:click=move|_|match serde_json::to_vec(collection.as_ref()){Ok(bytes)=>if let Err(e)=crate::app::download("oligo-collection.json","application/json",&bytes){status.set(e)},Err(e)=>status.set(e.to_string())}>"Download oligo collection"</button><button class="button secondary" on:click=move|_|match original.original_bytes(){Ok(bytes)=>if let Err(e)=crate::app::download("original-oligos.txt","text/plain",&bytes){status.set(e)},Err(e)=>status.set(e.to_string())}>"Download original oligo table"</button><button class="button secondary" on:click=move|_|match oligo_formats::export_oligos(&csv,Delimiter::Comma,TableExportMode::SpreadsheetSafe){Ok(report)=>{if let Err(e)=crate::app::download("oligo-order-sheet.csv","text/csv",report.text.as_bytes()){status.set(e)}else{status.set(report.warnings.join(" "))}},Err(e)=>status.set(e.to_string())}>"Spreadsheet-safe oligo CSV"</button><button class="button secondary" on:click=move|_|match oligo_formats::export_oligos(&tsv,Delimiter::Tab,TableExportMode::SpreadsheetSafe){Ok(report)=>{if let Err(e)=crate::app::download("oligo-order-sheet.tsv","text/tab-separated-values",report.text.as_bytes()){status.set(e)}else{status.set(report.warnings.join(" "))}},Err(e)=>status.set(e.to_string())}>"Spreadsheet-safe oligo TSV"</button></div>
            <p class="helper">"Spreadsheet-safe order sheets prefix formula-like cells with an apostrophe and report that change. Original fields and table JSON retain literal values. Whole collections are retained once in workspace archives; individual DNA documents retain their own row evidence."</p>
            </section>
        }})}
        </div>
    }
}
