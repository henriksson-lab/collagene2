//! Revision-bound unified search. Results are navigation evidence, never edits.
use crate::worker;
use collagene_core::{
    Molecule,
    protocol::{AnalysisJob, AnalysisResult, Request},
    search::{SearchHit, SearchMode, SearchOptions, SearchReport},
};
use leptos::prelude::*;
use wasm_bindgen_futures::spawn_local;

#[component]
pub fn SearchTools(
    molecule: Memo<Molecule>,
    identity: Memo<String>,
    epoch: RwSignal<u64>,
    selection: RwSignal<Option<(usize, usize)>>,
    selected_feature: RwSignal<Option<usize>>,
    regions: RwSignal<Option<(String, Vec<(usize, usize)>)>>,
    message: RwSignal<String>,
) -> impl IntoView {
    let query = RwSignal::new(String::new());
    let mode = RwSignal::new(SearchMode::All);
    let both = RwSignal::new(true);
    let case = RwSignal::new(false);
    let words = RwSignal::new(false);
    let open = RwSignal::new(false);
    let busy = RwSignal::new(false);
    let report = RwSignal::new(None::<SearchReport>);
    let page = RwSignal::new(0usize);
    let generation = RwSignal::new(0u32);
    let revision = Memo::new(move |_| molecule.with(|m| m.revision()));
    let options = Memo::new(move |_| SearchOptions {
        query: query.get(),
        mode: mode.get(),
        both_strands: both.get(),
        case_sensitive: case.get(),
        whole_word: words.get(),
        ..Default::default()
    });
    Effect::new(move |_| {
        let _ = (revision.get(), identity.get(), epoch.get());
        query.set(String::new());
        open.set(false);
    });
    Effect::new(move |_| {
        let _ = (revision.get(), identity.get(), epoch.get(), options.get());
        generation.update(|n| *n = n.wrapping_add(1));
        report.set(None);
        busy.set(false);
        page.set(0);
    });
    let select = Callback::new(move |hit: SearchHit| {
        let parts = hit
            .parts
            .iter()
            .map(|p| (p.start, p.end))
            .collect::<Vec<_>>();
        batch(|| {
            selected_feature.set(hit.feature_index);
            // Compound search hits retain every highlight without making an editable
            // span that silently includes unrelated bases or crosses the origin.
            selection.set(if parts.len() == 1 || hit.feature_index.is_some() {
                parts.first().copied()
            } else {
                None
            });
            regions.set(Some((revision.get_untracked(), parts)));
        });
    });
    let run = move || {
        let opts = options.get_untracked();
        if opts.query.trim().is_empty() {
            report.set(None);
            return;
        }
        let source = revision.get_untracked();
        let doc_id = identity.get_untracked();
        let workspace_epoch = epoch.get_untracked();
        let id = generation.get_untracked().wrapping_add(1);
        generation.set(id);
        busy.set(true);
        open.set(true);
        report.set(None);
        page.set(0);
        message.set("Searching in background…".into());
        spawn_local(async move {
            let result = worker::analyze(Request {
                id,
                revision: source.clone(),
                job: AnalysisJob::SearchDocument {
                    molecule: molecule.get_untracked(),
                    options: opts.clone(),
                },
            })
            .await;
            if generation.get_untracked() != id
                || revision.get_untracked() != source
                || identity.get_untracked() != doc_id
                || epoch.get_untracked() != workspace_epoch
                || options.get_untracked() != opts
            {
                return;
            }
            busy.set(false);
            match result {
                Ok(response) => match response.result {
                    AnalysisResult::DocumentSearch(value)
                        if value.source_revision == source && value.options == opts =>
                    {
                        if let Some(first) = value.hits.first() {
                            select.run(first.clone());
                        }
                        message.set(format!("{} match(es) · {} search. Ambiguous sites remain possible rather than certain.",value.hits.len(),if value.complete{"complete"}else{"incomplete; limits reached"}));
                        report.set(Some(value));
                    }
                    AnalysisResult::Error(error) => message.set(error),
                    _ => message.set("Unexpected search response".into()),
                },
                Err(error) => message.set(error),
            }
        });
    };
    view! {
        <div class="search-tools">
            <div class="search"><input aria-label="Find nucleotide sequence" placeholder="Find sequence or annotation…" prop:value=move||query.get() on:input=move|e|query.set(event_target_value(&e)) on:keydown=move|e|if e.key()=="Enter"{run();}/><button title="Find on both strands" disabled=move||busy.get() on:click=move|_|run()>"Find"</button><button title="Search options and results" aria-expanded=move||open.get() on:click=move|_|open.update(|v|*v=!*v)>"⌄"</button></div>
            <Show when=move||open.get()><section class="search-results" aria-label="Search options and results">
                <div class="search-heading"><strong>"Find in this document"</strong><button class="icon-button" aria-label="Close search results" on:click=move|_|open.set(false)>"×"</button></div>
                <label>"Search source"<select aria-label="Search source" prop:value=move||format!("{:?}",mode.get()) on:change=move|e|mode.set(match event_target_value(&e).as_str(){"Sequence"=>SearchMode::Sequence,"Annotations"=>SearchMode::Annotations,"Primers"=>SearchMode::Primers,"Enzymes"=>SearchMode::Enzymes,_=>SearchMode::All})><option value="All">"All"</option><option value="Sequence">"Sequence"</option><option value="Annotations">"Annotations and qualifiers"</option><option value="Primers">"Saved primers"</option><option value="Enzymes">"Enzyme recognition sites"</option></select></label>
                <div class="search-options"><label><input type="checkbox" prop:checked=move||both.get() on:change=move|e|both.set(event_target_checked(&e))/><span>"Both strands"</span></label><label><input type="checkbox" prop:checked=move||case.get() on:change=move|e|case.set(event_target_checked(&e))/><span>"Case-sensitive text"</span></label><label><input type="checkbox" prop:checked=move||words.get() on:change=move|e|words.set(event_target_checked(&e))/><span>"Whole text words"</span></label></div>
                <p class="helper">"Sequence queries use this document’s alphabet; nucleotide ambiguity symbols are supported. Enzymes use the starter catalog. Saved primers are primer_bind annotations. Changing options clears the previous results."</p>
                <Show when=move||busy.get()><p role="status">"Searching…"</p></Show>
                {move||report.get().map(|r|{let count=r.hits.len();let complete=r.complete;let has_warnings=!r.warnings.is_empty();view!{
                    <p class="search-summary" role="status">{format!("{count} results · {} · page {} of {}",if complete{"complete"}else{"incomplete"},page.get()+1,count.div_ceil(50).max(1))}</p>
                    <div class="search-hit-list">{r.hits.into_iter().skip(page.get()*50).take(50).map(|hit|{
                        let coords=hit.parts.iter().map(|p|format!("{}–{}{}",p.start+1,p.end,if p.reverse{" reverse"}else{""})).collect::<Vec<_>>().join(", ");
                        let cut=hit.enzyme.as_ref().map(|_|format!(" · cuts: top {}, bottom {}",boundary_label(hit.top_cut),boundary_label(hit.bottom_cut))).unwrap_or_default();
                        let label=hit.label.clone();let kind=hit.kind;let id=hit.id.clone();let certain=hit.certain;
                        view!{<button class="result-entry search-hit" data-hit-id=id on:click=move|_|select.run(hit.clone())><strong>{label}</strong><small>{format!("{kind:?} · {}{} · {}{cut}",if coords.is_empty(){"unresolved location"}else{&coords},if certain{""}else{" · possible / review"},if hit_parts_compound(&coords){"all parts highlighted"}else{"reference coordinates"})}</small></button>}
                    }).collect_view()}</div>
                    <Show when=move||{count>50}><div class="search-pages"><button disabled=move||page.get()==0 on:click=move|_|page.update(|p|*p=p.saturating_sub(1))>"Previous results"</button><button disabled=move||{(page.get()+1)*50>=count} on:click=move|_|page.update(|p|*p+=1)>"Next results"</button></div></Show>
                    <Show when=move||has_warnings><details><summary>"Search warnings"</summary>{r.warnings.clone().into_iter().map(|s|view!{<p class="helper">{s}</p>}).collect_view()}</details></Show>
                }})}
            </section></Show>
        </div>
    }
}
fn hit_parts_compound(coordinates: &str) -> bool {
    coordinates.contains(", ")
}

pub(crate) fn boundary_label(position: Option<usize>) -> String {
    match position {
        Some(0) => "at origin (before base 1)".into(),
        Some(p) => format!("after base {p}"),
        None => "outside reference".into(),
    }
}
