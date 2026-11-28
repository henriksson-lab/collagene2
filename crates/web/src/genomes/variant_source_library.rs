//! Workspace sources are reparsed on open, never treated as checked alignments.
use collagene_core::review::variants::archive::SourceInfo;
use leptos::prelude::*;
#[component]
pub fn VariantSourceLibrary(
    entries: Memo<Vec<(String, SourceInfo)>>,
    locked: Memo<bool>,
    on_open: Callback<String>,
    on_remove: Callback<String>,
) -> impl IntoView {
    let page = RwSignal::new(0usize);
    let visible = Memo::new(move |_| {
        let rows = entries.get();
        let last = rows.len().saturating_sub(1) / 10;
        let current = page.get().min(last);
        rows.into_iter()
            .skip(current * 10)
            .take(10)
            .collect::<Vec<_>>()
    });
    view! {<Show when=move||!entries.with(Vec::is_empty)><section class="saved-read-sources"><h3>"Variant sources"</h3><p class="helper">"Original VCF/BCF files. Open to parse records and check selected alleles against explicitly mapped DNA. Saved files do not establish allele presence or genotype."</p>
        {move||visible.get().into_iter().map(|(hash,info)|{let remove=hash.clone();view!{<article><strong>{info.name}</strong><p class="helper">{format!("{:?} · {} bytes · {}",info.format,info.source_bytes,info.source_blake3.chars().take(12).collect::<String>())}</p><button class="button secondary" disabled=move||locked.get() on:click=move|_|on_open.run(hash.clone())>"Open variant source"</button><button class="button secondary" disabled=move||locked.get() on:click=move|_|on_remove.run(remove.clone())>"Remove variant source"</button></article>}}).collect_view()}
        <div class="pagination"><button class="button secondary" disabled=move||page.get()==0 on:click=move|_|page.update(|p|*p=p.saturating_sub(1))>"Previous variant files"</button><button class="button secondary" disabled=move||(page.get()+1)*10>=entries.with(Vec::len) on:click=move|_|page.update(|p|*p+=1)>"Next variant files"</button></div>
    </section></Show>}
}
