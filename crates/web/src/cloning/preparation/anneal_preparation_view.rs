//! Inspect an actual two-original oligo preparation before pinning its geometry.
use crate::cloning::preparation::cloning_methods_model::MethodsState;
use crate::cloning::preparation::cloning_methods_model::OverlapPreparation;
use crate::platform::worker::CachedWorker;
use crate::sources::identity::SourceInfo;
use collagene_core::{
    chemistry::{self, AnnealOptions, FivePrime},
    protocol::{
        request::Request,
        schema::{AnalysisJob, AnalysisResult},
    },
    verification::terminal::AnnealPreview,
};
use collagene_model::Molecule;
use leptos::prelude::*;
use std::{
    rc::Rc,
    sync::{
        Arc,
        atomic::{AtomicBool, Ordering},
    },
};
use wasm_bindgen_futures::spawn_local;

#[derive(Clone, PartialEq, Eq)]
struct Context {
    top: SourceInfo,
    bottom: SourceInfo,
    top_document: usize,
    bottom_document: usize,
    options: AnnealOptions,
    epoch: u64,
}
#[component]
pub fn AnnealPreparationFields(
    config: RwSignal<MethodsState>,
    index: usize,
    query_options: Memo<Vec<(usize, String)>>,
    get_query: Callback<usize, Option<Molecule>>,
    source_info: Callback<usize, Option<SourceInfo>>,
    epoch: RwSignal<u64>,
    fixed_sources: bool,
) -> impl IntoView {
    let row = move || config.with(|c| c.rows.get(index).cloned().unwrap_or_default());
    let context = Memo::new(move |_| {
        let row = config.with(|c| c.rows.get(index).cloned())?;
        if row.preparation != OverlapPreparation::AnnealedOligos {
            return None;
        }
        let bottom_document = row.anneal_bottom_document?;
        Some(Context {
            top: source_info.run(row.document)?,
            bottom: source_info.run(bottom_document)?,
            top_document: row.document,
            bottom_document,
            options: row.anneal_options,
            epoch: epoch.get(),
        })
    });
    let status = RwSignal::new(String::new());
    let proof = RwSignal::new(None::<Arc<(Context, AnnealPreview)>>);
    let busy = RwSignal::new(false);
    let generation = RwSignal::new(0u32);
    let page = RwSignal::new(0usize);
    let active = StoredValue::new_local(None::<Rc<CachedWorker>>);
    let alive = Arc::new(AtomicBool::new(true));
    let cleanup = alive.clone();
    on_cleanup(move || {
        cleanup.store(false, Ordering::Relaxed);
        if let Some(worker) = active.get_value() {
            worker.close();
        }
    });
    let previous = StoredValue::new(None::<Context>);
    Effect::new(move |_| {
        let current = context.get();
        if previous.with_value(|old| old.is_some() && old != &current) {
            config.update(|c| {
                if let Some(row) = c.rows.get_mut(index) {
                    row.anneal_candidate = None;
                }
            });
        }
        previous.set_value(current);
        generation.update(|g| *g = g.wrapping_add(1));
        if let Some(worker) = active.get_value() {
            worker.close();
        }
        active.set_value(None);
        proof.set(None);
        busy.set(false);
        page.set(0);
    });
    let compute = Callback::new(move |()| {
        let Some(bound) = context.get_untracked() else {
            status.set("Pin both original oligos before inspecting annealing.".into());
            return;
        };
        if [bound.top.length, bound.bottom.length]
            .iter()
            .any(|n| *n == 0 || *n > chemistry::MAX_ANNEAL_OLIGO_BASES)
        {
            status.set("Each oligo must contain1–10000 DNA bases.".into());
            return;
        }
        let (Some(top), Some(bottom)) = (
            get_query.run(bound.top_document),
            get_query.run(bound.bottom_document),
        ) else {
            status.set("One pinned oligo is unavailable.".into());
            return;
        };
        if top.revision() != bound.top.binding.revision
            || bottom.revision() != bound.bottom.binding.revision
            || context.get_untracked() != Some(bound.clone())
        {
            status.set("An oligo source changed. Capture it again.".into());
            return;
        }
        let worker = match CachedWorker::new() {
            Ok(w) => Rc::new(w),
            Err(error) => {
                status.set(error);
                return;
            }
        };
        if let Some(old) = active.get_value() {
            old.close();
        }
        active.set_value(Some(worker.clone()));
        let token = generation.get_untracked().wrapping_add(1);
        generation.set(token);
        busy.set(true);
        proof.set(None);
        status.set("Reconstructing antiparallel oligo offsets…".into());
        let top_hash = blake3::hash(top.sequence.as_str().as_bytes())
            .to_hex()
            .to_string();
        let bottom_hash = blake3::hash(bottom.sequence.as_str().as_bytes())
            .to_hex()
            .to_string();
        let alive = alive.clone();
        spawn_local(async move {
            let response = worker
                .ask(Request {
                    id: token,
                    revision: bound.top.binding.revision.clone(),
                    job: AnalysisJob::AnnealPreview {
                        molecule: top,
                        bottom,
                        options: bound.options.clone(),
                    },
                })
                .await;
            worker.close();
            if !alive.load(Ordering::Relaxed)
                || generation.get_untracked() != token
                || context.get_untracked() != Some(bound.clone())
            {
                return;
            }
            active.set_value(None);
            busy.set(false);
            match response {
                Ok(response) if response.id==token&&response.revision==bound.top.binding.revision=>match response.result {
                    AnalysisResult::AnnealPreview(preview) if preview.top_revision==bound.top.binding.revision&&preview.bottom_revision==bound.bottom.binding.revision&&preview.annealing.options==bound.options&&preview.annealing.engine_version==chemistry::ANNEAL_ENGINE_VERSION&&preview.annealing.top_sequence_blake3==top_hash&&preview.annealing.bottom_sequence_blake3==bottom_hash=>{
                        status.set(format!("{} geometries · search complete {} · unique within this model {}. Choose an explicit geometry; stability and reaction success are not predicted.",preview.annealing.candidates.len(),preview.annealing.complete,preview.annealing.unique));
                        proof.set(Some(Arc::new((bound,preview))));
                    },
                    AnalysisResult::Error(error)=>status.set(error),
                    _=>status.set("The annealing response differs from the current source/settings binding.".into()),
                },
                Ok(_)=>status.set("Stale annealing response withheld.".into()),
                Err(error)=>status.set(error),
            }
        });
    });
    let choose = Callback::new(move |candidate: usize| {
        let Some(saved) = proof.get_untracked() else {
            return;
        };
        if context.get_untracked() != Some(saved.0.clone())
            || saved.1.annealing.candidates.get(candidate).is_none()
        {
            status.set(
                "Sources or settings changed; inspect the annealing alternatives again.".into(),
            );
            return;
        }
        config.update(|c| {
            if let Some(row) = c.rows.get_mut(index) {
                row.anneal_candidate = Some(candidate);
            }
        });
        status.set(format!("Geometry {} selected. Compute the complete source preparation and ligation against the independent target.",candidate+1));
    });
    view!{<section class="anneal-preparation-fields"><h4>"Two ordered oligos"</h4>
        <p class="helper">"Both original strands are written 5′→3′. Antiparallel offset and terminal chemistry are reconstructed, without automatic best-placement selection. Intended DNA remains the separately reviewed source-region declaration."</p>
        <Show when=move||!fixed_sources><label>"Bottom oligo source"<select aria-label=format!("Bottom oligo source {}",index+1) prop:value=move||row().anneal_bottom_document.map(|n|n.to_string()).unwrap_or_default() on:change=move|e|config.update(|c|if let Some(row)=c.rows.get_mut(index){row.anneal_bottom_document=event_target_value(&e).parse().ok();row.anneal_candidate=None;})><option value="">"Choose bottom oligo…"</option>{move||query_options.get().into_iter().map(|(i,name)|view!{<option value=i.to_string()>{name}</option>}).collect_view()}</select></label></Show>
        <p>{move||context.get().map(|c|format!("Top: {} · bottom: {}",c.top.name,c.bottom.name)).unwrap_or_else(||"Capture both original strands before reviewing this preparation.".into())}</p>
        <details><summary>"Annealing search and declared termini"</summary>
        <label>"Minimum paired bases"<input type="number" min="1" max="10000" aria-label=format!("Oligo minimum paired bases {}",index+1) prop:value=move||row().anneal_options.min_paired_bases on:input=move|e|config.update(|c|if let Some(row)=c.rows.get_mut(index){row.anneal_options.min_paired_bases=event_target_value(&e).parse().unwrap_or(0);row.anneal_candidate=None;})/></label>
        <label>"Comparison budget"<input type="number" min="1" max="200000000" aria-label=format!("Oligo comparison budget {}",index+1) prop:value=move||row().anneal_options.max_comparisons on:input=move|e|config.update(|c|if let Some(row)=c.rows.get_mut(index){row.anneal_options.max_comparisons=event_target_value(&e).parse().unwrap_or(0);row.anneal_candidate=None;})/></label>
        <label>"Maximum geometries"<input type="number" min="1" max="1000" aria-label=format!("Oligo maximum geometries {}",index+1) prop:value=move||row().anneal_options.max_candidates on:input=move|e|config.update(|c|if let Some(row)=c.rows.get_mut(index){row.anneal_options.max_candidates=event_target_value(&e).parse().unwrap_or(0);row.anneal_candidate=None;})/></label>
        <label class="checkbox-label"><input type="checkbox" prop:checked=move||row().anneal_options.top_chemistry.five_prime==FivePrime::Phosphate on:change=move|e|config.update(|c|if let Some(row)=c.rows.get_mut(index){row.anneal_options.top_chemistry.five_prime=if event_target_checked(&e){FivePrime::Phosphate}else{FivePrime::Hydroxyl};row.anneal_candidate=None;})/>"Top oligo ordered with 5′ phosphate"</label>
        <label class="checkbox-label"><input type="checkbox" prop:checked=move||row().anneal_options.bottom_chemistry.five_prime==FivePrime::Phosphate on:change=move|e|config.update(|c|if let Some(row)=c.rows.get_mut(index){row.anneal_options.bottom_chemistry.five_prime=if event_target_checked(&e){FivePrime::Phosphate}else{FivePrime::Hydroxyl};row.anneal_candidate=None;})/>"Bottom oligo ordered with 5′ phosphate"</label>
        <p class="helper">"The oligos are declared linear single strands with accessible 3′ OH. Phosphorylation is explicit. Terminal treatment after annealing is a separate preparation choice below."</p></details>
        <button class="button secondary" disabled=move||busy.get() on:click=move|_|compute.run(())>"Inspect annealing geometries"</button>
        <p class="helper" role="status">{move||status.get()}</p><p>{move||row().anneal_candidate.map(|i|format!("Selected geometry: {}",i+1)).unwrap_or_else(||"No geometry selected.".into())}</p>
        {move||proof.get().map(|saved|{let start=page.get().min(saved.1.annealing.candidates.len().saturating_sub(1)/20)*20;view!{<div class="anneal-candidate-review">{saved.1.annealing.warnings.iter().map(|warning|view!{<p class="helper">{warning.clone()}</p>}).collect_view()}{saved.1.annealing.candidates.iter().enumerate().skip(start).take(20).map(|(i,candidate)|view!{<details><summary>{format!("Geometry {} · {} paired bases · top-reference offset {}",i+1,candidate.paired_bases,candidate.duplex.bottom_start)}</summary>{candidate.duplex.summary().ok().as_ref().map(crate::cloning::preparation::end_view::summary)}<pre>{format!("Top 5′→3′: {}\nBottom 5′→3′: {}",candidate.duplex.top.as_str(),candidate.duplex.bottom.as_str())}</pre><button class="button secondary" on:click=move|_|choose.run(i)>"Use this geometry"</button></details>}).collect_view()}<div class="method-row-actions"><button class="text-tool" disabled=move||page.get()==0 on:click=move|_|page.update(|p|*p=p.saturating_sub(1))>"Previous geometries"</button><button class="text-tool" disabled=move||(page.get()+1)*20>=saved.1.annealing.candidates.len() on:click=move|_|page.update(|p|*p=p.saturating_add(1))>"Next geometries"</button></div></div>}})}
    </section>}.into_any()
}
