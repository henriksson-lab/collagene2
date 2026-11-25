//! Read-only projections of native finite-pool evidence; no verification is invented here.
use collagene_core::{
    Topology,
    chemistry::End,
    golden_gate::{GoldenGateFragment, GoldenGateReport},
};
use leptos::prelude::*;
use std::sync::Arc;

const PAGE: usize = 40;
fn pages(count: usize) -> usize {
    count.max(1).div_ceil(PAGE)
}
fn end_text(end: &End) -> String {
    let preview = end.overhang.chars().take(80).collect::<String>();
    format!(
        "{:?}; overhang 5′→3′ {}{} ({} bases); 5′ {:?}, 3′ {:?}",
        end.polarity,
        if preview.is_empty() { "none" } else { &preview },
        if end.overhang.len() > 80 {
            "… (preview)"
        } else {
            ""
        },
        end.overhang.len(),
        end.chemistry.five_prime,
        end.chemistry.three_prime
    )
}
fn fragment_text(pool: &GoldenGateReport, index: usize, fragment: &GoldenGateFragment) -> String {
    let source = pool
        .sources
        .iter()
        .find(|source| source.source_index == fragment.source_index);
    let name = source
        .map(|source| source.name.chars().take(160).collect::<String>())
        .unwrap_or_else(|| "Unavailable source".into());
    format!(
        "Pool fragment {} · source role {} ({}) · digest fragment {}\nSource revision {}\nOriginal source boundaries: top {} (zero-based), bottom {} (signed/unwrapped). Top {} bases, bottom {} bases, bottom offset {}, {:?}.\nLeft end: {}\nRight end: {}",
        index + 1,
        fragment.source_index + 1,
        name,
        fragment.fragment_index + 1,
        fragment.source_revision,
        fragment.top_start,
        fragment.bottom_start,
        fragment.duplex.top.len(),
        fragment.duplex.bottom.len(),
        fragment.duplex.bottom_start,
        fragment.duplex.topology,
        fragment
            .left_end
            .as_ref()
            .map(end_text)
            .unwrap_or_else(|| "No free end".into()),
        fragment
            .right_end
            .as_ref()
            .map(end_text)
            .unwrap_or_else(|| "No free end".into())
    )
}
#[component]
fn EvidencePager(page: RwSignal<usize>, total: usize, label: &'static str) -> impl IntoView {
    let count = pages(total);
    view! {<div class="evidence-pager"><p>{move||format!("{} {} · page {} of {}",total,label,page.get().min(count-1)+1,count)}</p>
        <button class="text-tool" disabled=move||page.get()==0 on:click=move|_|page.update(|value|*value=value.saturating_sub(1))>{format!("Previous {label}")}</button>
        <button class="text-tool" disabled=move||page.get()+1>=count on:click=move|_|page.update(|value|*value=(*value+1).min(count-1))>{format!("Next {label}")}</button>
    </div>}
}

#[component]
pub fn GoldenGateCandidateView(pool: Arc<GoldenGateReport>, index: usize) -> impl IntoView {
    let pool = StoredValue::new(pool);
    let path_page = RwSignal::new(0usize);
    let recut_page = RwSignal::new(0usize);
    let nick_page = RwSignal::new(0usize);
    let Some((directed,physical,geometry,ends,path_count,recut_count,nick_count))=pool.with_value(|pool|pool.candidates.get(index).map(|candidate|{
        let duplex=&candidate.duplex;
        let ends=if duplex.topology==Topology::Circular { "Circular duplex: no free ends.".into() } else {
            match duplex.ends() { Ok((left,right))=>format!("Left end: {}\nRight end: {}",end_text(&left),end_text(&right)),Err(error)=>format!("Free-end geometry unavailable: {error}") }
        };
        (candidate.directed_sequence_key.clone(),candidate.physical_sequence_key.clone(),
            format!("{:?} duplex · top {} bases · bottom {} bases · bottom start {} in zero-based top-reference coordinates. Recognition-site-free under the selected profile: {}.",duplex.topology,duplex.top.len(),duplex.bottom.len(),duplex.bottom_start,candidate.restriction_site_free),
            ends,candidate.path.len(),candidate.recut_sites.len(),duplex.nicks.len())
    })) else { return view!{<p>"Candidate unavailable in this evidence snapshot."</p>}.into_any(); };
    view!{<details class="golden-gate-candidate-evidence"><summary>"Candidate identity, physical geometry and source path"</summary>
        <p>"Exact directed sequence key"</p><code>{directed}</code><p>"Exact physical sequence key"</p><code>{physical}</code>
        <p class="helper">"Physical grouping does not identify a unique source or annotation history. Candidate evidence is a computational result under the declared finite-pool assumptions."</p>
        <p>{geometry}</p><pre>{ends}</pre>
        <details><summary>"Ordered raw source-fragment path"</summary>
            <p class="helper">"Path orientation applies to the raw fragment geometry shown below; the source boundaries and free ends remain in that fragment's original orientation."</p>
            {move||pool.with_value(|pool|pool.candidates.get(index).map(|candidate|candidate.path.iter().enumerate().skip(path_page.get()*PAGE).take(PAGE).map(|(step,part)|{
                let detail=pool.fragments.get(part.fragment).map(|fragment|fragment_text(pool,part.fragment,fragment)).unwrap_or_else(||"Unavailable pool fragment".into());
                view!{<article><h5>{format!("Path position {} · {} orientation",step+1,if part.reverse{"reverse"}else{"forward"})}</h5><pre>{detail}</pre></article>}
            }).collect_view()))}
            <EvidencePager page=path_page total=path_count label="path fragments"/>
        </details>
        <details><summary>"Recognition and recutting observations"</summary>
            <p class="helper">"Recognition starts below are one-based. Cut positions are zero-based boundaries; circular boundaries wrap modulo product length, while linear boundaries remain signed and unclipped."</p>
            {move||pool.with_value(|pool|pool.candidates.get(index).map(|candidate|candidate.recut_sites.iter().enumerate().skip(recut_page.get()*PAGE).take(PAGE).map(|(number,site)|view!{
                <p>{format!("Observation {} · recognition start {} · {} strand · top cut {} · bottom cut {} · fully double-stranded {} · cleavable under model {}",number+1,site.recognition_start+1,if site.reverse{"reverse"}else{"forward"},site.top_cut,site.bottom_cut,site.fully_double_stranded,site.cleavable_under_model)}</p>
            }).collect_view()))}
            <EvidencePager page=recut_page total=recut_count label="recut observations"/>
        </details>
        <details><summary>"Retained nicks"</summary>
            {move||pool.with_value(|pool|pool.candidates.get(index).map(|candidate|candidate.duplex.nicks.iter().skip(nick_page.get()*PAGE).take(PAGE).map(|nick|view!{<p>{format!("{:?} strand · boundary {} (zero-based top reference)",nick.strand,nick.position)}</p>}).collect_view()))}
            <EvidencePager page=nick_page total=nick_count label="nicks"/>
        </details>
    </details>}.into_any()
}

#[component]
pub fn PoolSearchNotes(pool: Arc<GoldenGateReport>) -> impl IntoView {
    let pool = StoredValue::new(pool);
    let note_page = RwSignal::new(0usize);
    let fragment_page = RwSignal::new(0usize);
    let carryover_page = RwSignal::new(0usize);
    let (notes,fragments,carryovers,summary)=pool.with_value(|pool|(
        pool.truncation_reasons.len()+pool.warnings.len(),pool.fragments.len(),pool.carryover_sources.len(),
        format!("{} directed candidates · {} physical sequence groups · {} states examined. Declared finite search complete: {}. Exhaustive real-mixture search: {}. Model: {}.",pool.candidates.len(),pool.stable_physical_product_count,pool.states_examined,pool.enumeration_complete,pool.exhaustive_mixture_search,pool.model)));
    view!{<section class="golden-gate-pool-notes"><p>{summary}</p>
        <details><summary>"Search warnings and truncation reasons"</summary>
            {move||pool.with_value(|pool|pool.truncation_reasons.iter().map(|text|("Truncation",text)).chain(pool.warnings.iter().map(|text|("Warning",text))).skip(note_page.get()*PAGE).take(PAGE).map(|(kind,text)|view!{<p>{format!("{kind}: {text}")}</p>}).collect_view())}
            <EvidencePager page=note_page total=notes label="search notes"/>
        </details>
        <details><summary>"All raw pool fragments"</summary>
            <p class="helper">"These include competing pieces, before path orientation. Source names are display labels (up to 160 characters); exact source revisions are retained below."</p>
            {move||pool.with_value(|pool|pool.fragments.iter().enumerate().skip(fragment_page.get()*PAGE).take(PAGE).map(|(index,fragment)|view!{<pre>{fragment_text(pool,index,fragment)}</pre>}).collect_view())}
            <EvidencePager page=fragment_page total=fragments label="raw fragments"/>
        </details>
        <details><summary>"Uncut circular carryover sources"</summary>
            <p class="helper">"Intact substrate carryover is separate from newly assembled candidates."</p>
            {move||pool.with_value(|pool|pool.carryover_sources.iter().skip(carryover_page.get()*PAGE).take(PAGE).map(|index|{
                let detail=pool.sources.iter().find(|source|source.source_index==*index).map(|source|format!("Role {} · {} · revision {}",index+1,source.name.chars().take(160).collect::<String>(),source.revision)).unwrap_or_else(||format!("Unavailable source role {}",index+1));
                view!{<p>{detail}</p>}
            }).collect_view())}
            <EvidencePager page=carryover_page total=carryovers label="carryover sources"/>
        </details>
    </section>}.into_any()
}
