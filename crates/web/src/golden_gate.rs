//! Revision-bound Type IIS pool exploration and explicit candidate-plan review.
use crate::worker;
use collagene_core::{
    Molecule,
    golden_gate::{GoldenGateEnzyme, GoldenGateOptions, GoldenGateReport},
    protocol::{AnalysisJob, AnalysisResult, Request, Response},
    transfer::Candidate,
    verification::{DesignPlan, Operation, PlanInput, Source},
};
use leptos::prelude::*;
use std::sync::{
    Arc,
    atomic::{AtomicBool, Ordering},
};
use wasm_bindgen_futures::spawn_local;

#[derive(Clone)]
struct BoundReport {
    target_revision: String,
    indices: Vec<usize>,
    templates: Vec<Molecule>,
    report: GoldenGateReport,
}
fn number(value: String, label: &str) -> Result<usize, String> {
    value
        .parse()
        .map_err(|_| format!("{label} must be a nonnegative integer"))
}
fn source_revisions(
    get_query: Callback<usize, Option<Molecule>>,
    indices: &[usize],
) -> Option<Vec<String>> {
    indices
        .iter()
        .map(|i| get_query.run(*i).map(|m| m.revision()))
        .collect()
}

#[component]
pub fn GoldenGatePanel(
    molecule: Memo<Molecule>,
    query_options: Memo<Vec<(usize, String)>>,
    get_query: Callback<usize, Option<Molecule>>,
    on_candidate: Callback<Candidate>,
) -> impl IntoView {
    let selected = RwSignal::new(Vec::<usize>::new());
    let enzyme = RwSignal::new(GoldenGateEnzyme::BsaI);
    let circular = RwSignal::new(true);
    let both = RwSignal::new(true);
    let minimum = RwSignal::new("1".to_string());
    let maximum = RwSignal::new("16".to_string());
    let max_products = RwSignal::new("128".to_string());
    let max_states = RwSignal::new("50000".to_string());
    let report = RwSignal::new(None::<BoundReport>);
    let page = RwSignal::new(0usize);
    let status = RwSignal::new(String::new());
    let busy = RwSignal::new(false);
    let generation = RwSignal::new(0u32);
    let alive = Arc::new(AtomicBool::new(true));
    {
        let alive = alive.clone();
        on_cleanup(move || alive.store(false, Ordering::Relaxed));
    }
    let options = Memo::new(move |_| {
        Ok::<_, String>(GoldenGateOptions {
            enzyme: enzyme.get(),
            circular: circular.get(),
            both_orientations: both.get(),
            min_fragments: number(minimum.get(), "Minimum fragments")?,
            max_fragments: number(maximum.get(), "Maximum fragments")?,
            max_products: number(max_products.get(), "Product limit")?,
            max_states: number(max_states.get(), "State limit")?,
        })
    });
    Effect::new(move |_| {
        let _ = molecule.get().revision();
        let _ = query_options.get();
        let _ = source_revisions(get_query, &selected.get());
        let _ = options.get();
        generation.update(|g| *g = g.wrapping_add(1));
        report.set(None);
        page.set(0);
        status.set(String::new());
        busy.set(false);
    });
    let compute = move |_| {
        let target = molecule.get_untracked();
        let target_revision = target.revision();
        let indices = selected.get_untracked();
        let templates = untrack(|| {
            indices
                .iter()
                .map(|i| get_query.run(*i))
                .collect::<Option<Vec<_>>>()
        });
        let Some(templates) = templates else {
            status.set("An input document disappeared; choose the current templates.".into());
            return;
        };
        if templates.is_empty() {
            status.set("Choose at least one open template to digest.".into());
            return;
        }
        let parameters = match options.get_untracked() {
            Ok(o) => o,
            Err(e) => {
                status.set(e);
                return;
            }
        };
        let expected_revisions = templates.iter().map(Molecule::revision).collect::<Vec<_>>();
        let token = generation.get_untracked().wrapping_add(1);
        generation.set(token);
        busy.set(true);
        status.set("Digesting and enumerating finite pool products in background…".into());
        let alive = alive.clone();
        spawn_local(async move {
            let result = worker::analyze(Request {
                id: token,
                revision: target_revision.clone(),
                job: AnalysisJob::GoldenGate {
                    molecule: target,
                    templates: templates.clone(),
                    options: parameters,
                },
            })
            .await;
            if !alive.load(Ordering::Relaxed) {
                return;
            }
            if generation.get_untracked() != token
                || molecule.get_untracked().revision() != target_revision
                || selected.get_untracked() != indices
                || untrack(|| source_revisions(get_query, &indices)) != Some(expected_revisions)
            {
                return;
            }
            busy.set(false);
            match result {
                Ok(Response {
                    id,
                    revision,
                    result: AnalysisResult::GoldenGate(computed),
                }) if id == token && revision == target_revision => {
                    status.set(format!("{} directed candidates returned; finite search {}. Real mixture enumeration is not exhaustive.",computed.candidates.len(),if computed.enumeration_complete{"complete under its one-copy/orientation scope"}else{"incomplete"}));
                    report.set(Some(BoundReport {
                        target_revision,
                        indices,
                        templates,
                        report: computed,
                    }));
                }
                Ok(Response {
                    result: AnalysisResult::Error(error),
                    ..
                })
                | Err(error) => status.set(error),
                _ => status.set("Unexpected Golden Gate worker response".into()),
            }
        });
    };
    view! {
        <style>".gate-panel .gate-template{display:flex;align-items:flex-start;gap:7px;font-size:12px;padding:6px 0}.gate-panel .gate-template input{width:auto;flex-shrink:0;margin-top:3px}.gate-panel .gate-grid{display:grid;grid-template-columns:1fr 1fr;gap:8px}.gate-panel .gate-grid label{margin:0}.gate-panel .gate-candidate,.gate-panel .gate-fragment{border-top:1px solid #e1e7e0;padding:10px 0;font-size:12px;line-height:1.6;overflow-wrap:anywhere}.gate-panel pre{font:11px/1.6 monospace;white-space:pre-wrap;overflow-wrap:anywhere;max-height:160px;overflow:auto}.gate-panel .gate-warning{color:#8b4d24;font-size:12px;line-height:1.6}.gate-panel ul{padding-left:18px}.gate-panel .gate-status{font-size:12px;line-height:1.6;overflow-wrap:anywhere}.gate-panel details{margin:12px 0}"</style>
        <div class="inspector-content gate-panel"><h2>"Golden Gate"</h2><p class="helper">"Choose templates for complete Type IIS digestion. The active document is the independently defined target, not an inferred assembly result."</p>
        <details open><summary>"Input templates"</summary>{move||query_options.get().into_iter().take(1000).map(|(index,name)|view!{<label class="gate-template"><input type="checkbox" aria-label=format!("Golden Gate template {}",index+1) prop:checked=move||selected.with(|indices|indices.contains(&index)) on:change=move|e|{let checked=event_target::<web_sys::HtmlInputElement>(&e).checked();selected.update(|indices|{if checked{if !indices.contains(&index){indices.push(index);indices.sort_unstable();}}else{indices.retain(|i|*i!=index);}});}/><span>{name.chars().take(160).collect::<String>()}</span></label>}).collect_view()}</details>
        <label>"Type IIS enzyme"<select aria-label="Golden Gate enzyme" prop:value=move||format!("{:?}",enzyme.get()) on:change=move|e|enzyme.set(match event_target_value(&e).as_str(){"BsmBI"=>GoldenGateEnzyme::BsmBI,"BbsI"=>GoldenGateEnzyme::BbsI,"SapI"=>GoldenGateEnzyme::SapI,_=>GoldenGateEnzyme::BsaI})><option value="BsaI">"BsaI · GGTCTC(1/5)"</option><option value="BsmBI">"BsmBI · CGTCTC(1/5)"</option><option value="BbsI">"BbsI · GAAGAC(2/6)"</option><option value="SapI">"SapI · GCTCTTC(1/4)"</option></select></label>
        <label class="gate-template"><input type="checkbox" aria-label="Golden Gate circular product" prop:checked=move||circular.get() on:change=move|e|circular.set(event_target::<web_sys::HtmlInputElement>(&e).checked())/>"Require a circular product"</label>
        <label class="gate-template"><input type="checkbox" aria-label="Golden Gate both orientations" prop:checked=move||both.get() on:change=move|e|both.set(event_target::<web_sys::HtmlInputElement>(&e).checked())/>"Explore both fragment orientations"</label>
        <details><summary>"Finite search bounds"</summary><div class="gate-grid">{[("Minimum fragments","Golden Gate minimum fragments",minimum),("Maximum fragments","Golden Gate maximum fragments",maximum),("Product limit","Golden Gate product limit",max_products),("State limit","Golden Gate state limit",max_states)].into_iter().map(|(label,aria,value)|view!{<label>{label}<input aria-label=aria type="number" prop:value=move||value.get() on:input=move|e|value.set(event_target_value(&e))/></label>}).collect_view()}</div></details>
        <p class="helper">"Includes self-ligation by default. Each complete-digest fragment is available at most once; concatemers, repeated copies, kinetics and partial digestion remain outside this search."</p>
        <button class="button primary full" disabled=move||busy.get() on:click=compute>"Enumerate Golden Gate products"</button>
        <p class="gate-status" role="status">{move||status.get()}</p>
        {move||report.get().map(|bound|{
            let bound=Arc::new(bound);let download=bound.clone();
            let count=bound.report.candidates.len();let physical=bound.report.stable_physical_product_count;let start=page.get()*50;
            let complete=bound.report.enumeration_complete;let carryover=bound.report.carryover_sources.len();
            let warnings=bound.report.warnings.clone();let truncation=bound.report.truncation_reasons.clone();
            let fragments=bound.report.fragments.clone();
            view!{<section class="gate-report"><h3>{format!("{count} directed assembly candidates")}</h3><p class="helper">{format!("{physical} physical sequence groups among returned site-free candidates; {carryover} uncut circular carryover sources. Sequence grouping does not establish a unique source/annotation history.")}</p>
                <p class="gate-completeness">{if complete{"Finite one-copy search complete under the requested orientation scope."}else{"Search incomplete: additional finite-pool candidates may exist."}}</p><ul>{truncation.into_iter().chain(warnings).map(|warning|view!{<li class="gate-warning">{warning}</li>}).collect_view()}</ul>
                <details><summary>{format!("{} complete-digest fragments",fragments.len())}</summary>{fragments.into_iter().enumerate().map(|(index,fragment)|{let left=fragment.left_end.map(|e|format!("{:?} {}",e.polarity,e.overhang)).unwrap_or_else(||"no free end".into());let right=fragment.right_end.map(|e|format!("{:?} {}",e.polarity,e.overhang)).unwrap_or_else(||"no free end".into());view!{<article class="gate-fragment"><strong>{format!("Pool {} · template {} / digest {}",index+1,fragment.source_index+1,fragment.fragment_index+1)}</strong><p>{format!("Top {} nt · bottom {} nt · source top boundary {} · bottom boundary {}",fragment.duplex.top.len(),fragment.duplex.bottom.len(),fragment.top_start,fragment.bottom_start)}</p><p>{format!("Left: {left}; right: {right}")}</p></article>}}).collect_view()}</details>
                {bound.report.candidates.iter().enumerate().skip(start).take(50).map(|(index,candidate)|{
                    let source=bound.clone();let source_for_path=bound.clone();
                    let path=candidate.path.iter().map(|p|{let f=&source_for_path.report.fragments[p.fragment];format!("T{} / F{}{}",f.source_index+1,f.fragment_index+1,if p.reverse{" reverse"}else{" forward"})}).collect::<Vec<_>>().join(" → ");
                    let recognition=candidate.recut_sites.iter().take(20).map(|s|format!("Recognition at {}{}; paired {}; cleavable {}",s.recognition_start+1,if s.reverse{" reverse"}else{" forward"},s.fully_double_stranded,s.cleavable_under_model)).collect::<Vec<_>>().join("\n");
                    let has_recognition=!recognition.is_empty();
                    view!{<article class="gate-candidate"><strong>{format!("Candidate {} · {} bp · {} fragments",index+1,candidate.duplex.top.len(),candidate.path.len())}</strong><p class="gate-site-status">{if candidate.restriction_site_free{"Restriction-site-free reference"}else{"Recognition retained/recreated: review as a recuttable or unresolved candidate"}}</p><p>{path}</p><Show when=move||has_recognition><pre>{recognition.clone()}</pre></Show><button class="button secondary full" on:click=move|_|{
                        if molecule.get_untracked().revision()!=source.target_revision||untrack(||source_revisions(get_query,&source.indices))!=Some(source.templates.iter().map(Molecule::revision).collect()){
                            status.set("Target or templates changed; rerun Golden Gate enumeration.".into());return;
                        }
                        let inputs=source.templates.iter().enumerate().map(|(i,m)|PlanInput{id:format!("template_{}",i+1),revision:m.revision(),molecule:m.clone()}).collect::<Vec<_>>();
                        let templates=inputs.iter().map(|i|Source::Input(i.id.clone())).collect();
                        on_candidate.run(Candidate::new(DesignPlan{name:format!("Golden Gate candidate {}",index+1),target_revision:source.target_revision.clone(),inputs,steps:vec![Operation::GoldenGate{templates,options:source.report.options.clone(),candidate_index:index}]}));
                    }>{format!("Review Golden Gate candidate {}",index+1)}</button></article>}
                }).collect_view()}
                <Show when=move||{count>50}><p class="helper">{format!("Showing candidates {}–{} of {count}; download retains the full returned report. Capped enumeration remains explicitly incomplete.",start+1,(start+50).min(count))}</p><div class="gate-grid"><button class="button secondary" disabled=move||page.get()==0 on:click=move|_|page.update(|p|*p=p.saturating_sub(1))>"Previous candidates"</button><button class="button secondary" disabled=move||{(page.get()+1)*50>=count} on:click=move|_|page.update(|p|*p+=1)>"Next candidates"</button></div></Show>
                <button class="button secondary full" on:click=move|_|match serde_json::to_vec_pretty(&download.report){Ok(bytes)=>crate::app::download_with_completion("golden-gate-report.json","application/json",&bytes,move|result|{if status.try_get_untracked().is_none(){return;}if let Err(e)=result {status.set(e)}}),Err(e)=>status.set(e.to_string())}>"Download Golden Gate report"</button>
                <p class="helper">"Review opens a typed proposal. Independently verify it in the Plan panel against the current target; importing a candidate grants no scientific status."</p>
            </section>}
        })}
        </div>
    }
}
