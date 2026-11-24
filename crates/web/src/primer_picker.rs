//! Explicit source/event-bound reuse of saved primer annotations. No automatic
//! oligo selection or imported-qualifier authority; copies remain editable.
use crate::digest::SourceInfo;
use collagene_core::{
    analysis::Primer,
    primer_reuse::{List, Resolved},
    run_lineage::EventBinding,
};
use leptos::prelude::*;
use std::sync::Arc;
#[derive(Clone, Debug)]
pub struct ResolveRequest {
    pub source: EventBinding,
    pub feature_index: usize,
}
#[component]
pub fn PrimerPicker(
    label: String,
    sources: Memo<Vec<SourceInfo>>,
    list_primers: Callback<EventBinding, Result<List, String>>,
    resolve_primer: Callback<ResolveRequest, Result<Resolved, String>>,
    epoch: RwSignal<u64>,
    on_copy: Callback<Primer>,
) -> impl IntoView {
    let selected = RwSignal::new(None::<EventBinding>);
    let listed = RwSignal::new(None::<Arc<List>>);
    let copied = RwSignal::new(None::<Arc<Resolved>>);
    let status = RwSignal::new(String::new());
    let query = RwSignal::new(String::new());
    let page = RwSignal::new(0usize);
    Effect::new(move |_| {
        let current = sources.get();
        let _ = epoch.get();
        if selected
            .get_untracked()
            .is_some_and(|binding| !current.iter().any(|s| s.binding == binding))
        {
            selected.set(None);
            listed.set(None);
            copied.set(None);
            page.set(0);
            status.set(
                "Primer source changed; select its current revision before copying again.".into(),
            );
        }
    });
    Effect::new(move |_| {
        let _ = query.get();
        page.set(0);
    });
    let choose = move |event: web_sys::Event| {
        let id = event_target_value(&event);
        selected.set(None);
        listed.set(None);
        copied.set(None);
        page.set(0);
        status.set(String::new());
        if id.is_empty() {
            return;
        }
        let Some(source) = sources
            .with_untracked(|items| items.iter().find(|s| s.binding.document_id == id).cloned())
        else {
            status.set("Primer source is unavailable.".into());
            return;
        };
        let binding = source.binding;
        match list_primers.run(binding.clone()) {
            Ok(list)
                if list.source_revision == binding.revision
                    && sources
                        .with_untracked(|items| items.iter().any(|s| s.binding == binding)) =>
            {
                selected.set(Some(binding));
                listed.set(Some(Arc::new(list)));
            }
            Ok(_) => status.set("Primer source changed while listing; select it again.".into()),
            Err(error) => status.set(error),
        }
    };
    let matches = Memo::new(move |_| {
        let q = query.get().to_lowercase();
        listed.with(|list| {
            list.as_ref()
                .map(|list| {
                    list.primers
                        .iter()
                        .filter(|p| {
                            q.is_empty()
                                || p.name.to_lowercase().contains(&q)
                                || p.feature_index.to_string().contains(&q)
                        })
                        .cloned()
                        .collect::<Vec<_>>()
                })
                .unwrap_or_default()
        })
    });
    let select_primer = Callback::new(move |feature_index: usize| {
        let Some(binding) = selected.get_untracked() else {
            return;
        };
        let captured_epoch = epoch.get_untracked();
        if !sources.with_untracked(|items| items.iter().any(|s| s.binding == binding)) {
            status.set("Primer source changed; choose it again.".into());
            return;
        }
        match resolve_primer.run(ResolveRequest {
            source: binding.clone(),
            feature_index,
        }) {
            Ok(value)
                if value.source_revision == binding.revision
                    && value.feature_index == feature_index
                    && epoch.get_untracked() == captured_epoch
                    && sources
                        .with_untracked(|items| items.iter().any(|s| s.binding == binding)) =>
            {
                on_copy.run(value.primer.clone());
                status.set("Copied explicit primer name, annealing domain and 5′ tail. Review the editable fields.".into());
                copied.set(Some(Arc::new(value)));
            }
            Ok(_) => status.set("Primer identity changed before copying; choose it again.".into()),
            Err(error) => status.set(error),
        }
    });
    let source_label = format!("{label} saved primer source");
    let search_label = format!("Search {label} saved primers");
    view! {<details class="saved-primer-picker"><summary>{format!("Reuse saved {label} primer")}</summary>
        <p class="helper">"Choose an open source explicitly, then copy a saved primer annotation. Imported oligo documents are supported when their annotations retain consistent sequence domains. No primer is selected automatically."</p>
        <label>"Saved primer document"<select aria-label=source_label prop:value=move||selected.get().map(|b|b.document_id).unwrap_or_default() on:change=choose><option value="">"Choose a document…"</option>{move||sources.get().into_iter().map(|s|view!{<option value=s.binding.document_id>{s.name}</option>}).collect_view()}</select></label>
        <Show when=move||listed.with(Option::is_some)><label>"Find saved primer"<input aria-label=search_label.clone() prop:value=move||query.get() on:input=move|e|query.set(event_target_value(&e))/></label>
            {move||listed.get().map(|list|view!{<p class="helper">{format!("{} primer annotations; {} listed; {} omitted by the bounded list limit. Source revision {}.",list.total,list.primers.len(),list.omitted,list.source_revision)}</p>})}
            {move||{let choices=matches.get();let pages=choices.len().max(1).div_ceil(50);let current=page.get().min(pages-1);view!{<p>{format!("{} matching primers · page {} of {}",choices.len(),current+1,pages)}</p>{choices.into_iter().skip(current*50).take(50).map(|p|{let index=p.feature_index;view!{<article class="saved-primer-choice"><strong>{p.name}</strong><p>{format!("Annotation {} · annealing {} bp · tail {} bp · ordered {} bp{}",index+1,p.annealing_bases.map(|n|n.to_string()).unwrap_or_else(||"?".into()),p.tail_bases.map(|n|n.to_string()).unwrap_or_else(||"?".into()),p.ordered_bases.map(|n|n.to_string()).unwrap_or_else(||"?".into()),if p.needs_review{" · source review flag"}else{""})}</p><p class="helper">{p.diagnostic.unwrap_or_default()}</p><button class="button secondary" disabled=!p.reusable on:click=move|_|select_primer.run(index)>"Copy this primer"</button></article>}}).collect_view()}<div class="method-row-actions"><button class="text-tool" disabled=current==0 on:click=move|_|page.set(current.saturating_sub(1))>"Previous saved primers"</button><button class="text-tool" disabled={current+1>=pages} on:click=move|_|page.set((current+1).min(pages-1))>"Next saved primers"</button></div>}}}
        </Show>
        <p role="status">{move||status.get()}</p>
        {move||copied.get().map(|resolved|view!{<details open><summary>"Last explicit copy — review full oligo"</summary><p>{resolved.primer.name.clone()}</p><pre>{format!("Ordered 5′→3′: {}{}\nAnnealing 5′→3′: {}\n5′ tail: {}",resolved.primer.tail,resolved.primer.annealing.as_str(),resolved.primer.annealing.as_str(),if resolved.primer.tail.is_empty(){"(none)"}else{&resolved.primer.tail})}</pre>{resolved.warnings.iter().map(|w|view!{<p class="helper">{w.clone()}</p>}).collect_view()}<p class="helper">"Copied sequence domains are editable design inputs; annotation identity and source review flags do not establish experimental binding or success."</p></details>})}
    </details>}
}
