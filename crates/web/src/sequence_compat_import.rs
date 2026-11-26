//! Explicit source/alphabet review for formats which cannot safely be guessed.
use collagene_core::{Alphabet, formats::ImportReport};
use leptos::prelude::*;
use std::sync::Arc;
#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Profile {
    Gcg,
    CollageneXml,
}
#[component]
pub fn SequenceCompatImport(
    source: Arc<Vec<u8>>,
    name: String,
    profile: Profile,
    guard: Callback<(), Result<(), String>>,
    on_import: Callback<ImportReport, Result<(), String>>,
) -> impl IntoView {
    let alphabet = RwSignal::new("dna".to_string());
    let lower = RwSignal::new(false);
    let prepared = RwSignal::new(None::<Arc<ImportReport>>);
    let error = RwSignal::new(String::new());
    let source = StoredValue::new(source);
    let filename = StoredValue::new(name.clone());
    let prepare = move |_| {
        prepared.set(None);
        error.set(String::new());
        if let Err(e) = guard.run(()) {
            error.set(e);
            return;
        }
        let bytes = source.get_value();
        let result = match profile {
            Profile::Gcg => collagene_core::gcg_sequence::import_gcg(
                &bytes,
                &collagene_core::gcg_sequence::ImportOptions {
                    alphabet: match alphabet.get_untracked().as_str() {
                        "rna" => Alphabet::Rna,
                        "protein" => Alphabet::Protein,
                        _ => Alphabet::Dna,
                    },
                    fallback_name: filename.get_value(),
                },
            ),
            Profile::CollageneXml => collagene_core::collagene_xml_compat::import_with_policy(
                &bytes,
                if lower.get_untracked() {
                    collagene_core::collagene_xml_compat::ReferencePolicy::AllowLowerOnly
                } else {
                    collagene_core::collagene_xml_compat::ReferencePolicy::UpperOnly
                },
            ),
        };
        match result {
            Ok(r) => prepared.set(Some(Arc::new(r))),
            Err(e) => error.set(e.to_string()),
        }
    };
    let commit = move |_| {
        if let Err(e) = guard.run(()) {
            error.set(e);
            return;
        }
        if let Some(report) = prepared.get_untracked() {
            if let Err(e) = on_import.run((*report).clone()) {
                error.set(e);
            }
        }
    };
    view! {<div class="sequence-compat-import"><p>{format!("Source: {name}")}</p>
    <Show when=move||profile==Profile::Gcg><label>"Sequence alphabet"<select prop:value=move||alphabet.get() on:change=move|e|{alphabet.set(event_target_value(&e));prepared.set(None);}><option value="dna">"DNA"</option><option value="rna">"RNA"</option><option value="protein">"Protein"</option></select></label><p class="helper">"GCG nucleotide records do not distinguish DNA from RNA. Choose the alphabet explicitly; the declared type, length and checksum will be checked."</p></Show>
    <Show when=move||profile==Profile::CollageneXml><p class="helper">"This project contains a record with only a lower strand. Importing its literal 5′→3′ sequence requires review. Original aligned rows and source bytes remain retained; this migration does not reconstruct cloning chemistry."</p><label><input type="checkbox" prop:checked=move||lower.get() on:change=move|e|{lower.set(event_target_checked(&e));prepared.set(None);}/>"Allow lower-only sequence references; retain uncertain annotations for review"</label></Show>
    <button class="button secondary" on:click=prepare>"Review import"</button>
    <p role="alert">{move||error.get()}</p>
    {move||prepared.get().map(|report|{let count=report.molecules.len();let warnings=report.warnings.clone();view!{<div><p>{format!("{count} sequence record(s) ready. No records have been added yet.")}</p><ul>{warnings.into_iter().map(|w|view!{<li>{w}</li>}).collect_view()}</ul><button class="button" on:click=commit>"Import reviewed records"</button></div>}})}
    </div>}
}
