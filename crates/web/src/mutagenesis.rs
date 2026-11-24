//! Revision-bound mutation design and conditional mismatch-PCR reference previews.
use crate::worker;
use collagene_core::{
    Molecule,
    analysis::Primer,
    mutagenesis::{InversePcrDesign, InversePcrOptions, PcrOptions, PcrReport},
    primers::BindingOptions,
    protocol::{AnalysisJob, AnalysisResult, Request, Response},
    transfer::Candidate,
};
use leptos::prelude::*;
use std::sync::{
    Arc,
    atomic::{AtomicBool, Ordering},
};
use wasm_bindgen_futures::spawn_local;
fn number<T: std::str::FromStr>(value: String, label: &str) -> Result<T, String> {
    value
        .trim()
        .parse()
        .map_err(|_| format!("{label} is not a valid number"))
}
#[derive(Clone)]
struct PcrEvidence {
    report: PcrReport,
    documents: Vec<Molecule>,
    findings: Vec<String>,
}
fn preview_document(bound: &PcrEvidence, index: usize) -> Result<Molecule, String> {
    bound
        .documents
        .get(index)
        .cloned()
        .ok_or_else(|| "Product preview no longer exists".into())
}
#[component]
pub fn MutagenesisPanel(
    molecule: Memo<Molecule>,
    query_options: Memo<Vec<(usize, String)>>,
    get_query: Callback<usize, Option<Molecule>>,
    on_candidate: Callback<Candidate>,
    on_products: Callback<Vec<Molecule>, Result<(), String>>,
) -> impl IntoView {
    let template_index = RwSignal::new(0usize);
    let source = Memo::new(move |_| get_query.run(template_index.get()));
    let f_domain = RwSignal::new(String::new());
    let r_domain = RwSignal::new(String::new());
    let f_tail = RwSignal::new(String::new());
    let r_tail = RwSignal::new(String::new());
    let mismatches = RwSignal::new("1".to_string());
    let exact = RwSignal::new("1".to_string());
    let edit_start = RwSignal::new("1".to_string());
    let edit_end = RwSignal::new("1".to_string());
    let inserted = RwSignal::new(String::new());
    let f_length = RwSignal::new("20".to_string());
    let r_length = RwSignal::new("20".to_string());
    let split = RwSignal::new("0".to_string());
    let pcr = RwSignal::new(None::<Arc<PcrEvidence>>);
    let inverse = RwSignal::new(None::<Arc<InversePcrDesign>>);
    let status = RwSignal::new(String::new());
    let busy = RwSignal::new(false);
    let generation = RwSignal::new(0u32);
    let alive = Arc::new(AtomicBool::new(true));
    {
        let alive = alive.clone();
        on_cleanup(move || alive.store(false, Ordering::Relaxed));
    }
    Effect::new(move |_| {
        let _ = (
            molecule.get().revision(),
            source.get().map(|m| m.revision()),
            template_index.get(),
            f_domain.get(),
            r_domain.get(),
            f_tail.get(),
            r_tail.get(),
            mismatches.get(),
            exact.get(),
            edit_start.get(),
            edit_end.get(),
            inserted.get(),
            f_length.get(),
            r_length.get(),
            split.get(),
        );
        generation.update(|n| *n = n.wrapping_add(1));
        pcr.set(None);
        inverse.set(None);
        busy.set(false);
        status.set(String::new());
    });
    let run = {
        let alive = alive.clone();
        move |job: AnalysisJob, expected_source: Molecule, target_revision: Option<String>| {
            let token = generation.get_untracked().wrapping_add(1);
            generation.set(token);
            busy.set(true);
            status.set("Computing mutation evidence in background…".into());
            let revision = expected_source.revision();
            let alive = alive.clone();
            spawn_local(async move {
                let response = worker::analyze(Request {
                    id: token,
                    revision: revision.clone(),
                    job,
                })
                .await;
                if !alive.load(Ordering::Relaxed) || generation.get_untracked() != token {
                    return;
                }
                if source
                    .get_untracked()
                    .is_none_or(|m| m.revision() != revision)
                    || target_revision.is_some_and(|r| molecule.get_untracked().revision() != r)
                {
                    return;
                }
                busy.set(false);
                match response {
                    Ok(Response {
                        id,
                        revision: returned,
                        result,
                    }) if id == token && returned == revision => match result {
                        AnalysisResult::MismatchPcrPreview(preview) => {
                            let report = preview.report;
                            status.set(format!("{} candidate geometries · search complete: {} · complete under declared model: {}. No empirical amplification claim.",report.products.len(),report.search_complete,report.enumeration_complete));
                            pcr.set(Some(Arc::new(PcrEvidence {
                                report,
                                documents: preview.documents,
                                findings: preview.findings,
                            })));
                        }
                        AnalysisResult::InversePcr(report) => {
                            status.set(format!("Independent inverse-PCR plan check: {:?}. Review the bound steps and target before use.",report.verification.status));
                            inverse.set(Some(Arc::new(report)));
                        }
                        AnalysisResult::Error(error) => status.set(error),
                        _ => status.set("Unexpected mutation worker result".into()),
                    },
                    Err(error) => status.set(error),
                    _ => status.set("Mutation worker identity changed".into()),
                }
            });
        }
    };
    let run_pcr = run.clone();
    let simulate = move |_| {
        let result = (|| -> Result<_, String> {
            let source = source.get_untracked().ok_or("Choose a template document")?;
            let forward = Primer::new(
                "Mutagenic forward",
                &f_domain.get_untracked(),
                &f_tail.get_untracked(),
            )
            .map_err(|e| e.to_string())?;
            let reverse = Primer::new(
                "Mutagenic reverse",
                &r_domain.get_untracked(),
                &r_tail.get_untracked(),
            )
            .map_err(|e| e.to_string())?;
            let options = PcrOptions {
                binding: BindingOptions {
                    max_mismatches: number(mismatches.get_untracked(), "Allowed mismatches")?,
                    min_three_prime_exact: number(exact.get_untracked(), "Exact 3′ bases")?,
                    ..Default::default()
                },
                ..Default::default()
            };
            Ok((source, forward, reverse, options))
        })();
        match result {
            Ok((source, forward, reverse, options)) => run_pcr(
                AnalysisJob::MismatchPcrPreview {
                    molecule: source.clone(),
                    forward: forward.clone(),
                    reverse: reverse.clone(),
                    options,
                },
                source,
                None,
            ),
            Err(e) => status.set(e),
        }
    };
    let design = move |_| {
        let result = (|| -> Result<_, String> {
            let source = source.get_untracked().ok_or("Choose a source template")?;
            let target = molecule.get_untracked();
            let start: usize = number(edit_start.get_untracked(), "Edit start")?;
            let end: usize = number(edit_end.get_untracked(), "Edit end")?;
            let options = InversePcrOptions {
                edit_start: start.checked_sub(1).ok_or("Edit start is one-based")?,
                edit_end: end.checked_sub(1).ok_or("Edit end is one-based")?,
                inserted: inserted.get_untracked(),
                forward_length: number(f_length.get_untracked(), "Forward domain length")?,
                reverse_length: number(r_length.get_untracked(), "Reverse domain length")?,
                reverse_tail_bases: number(split.get_untracked(), "Reverse-tail split")?,
            };
            Ok((source, target, options))
        })();
        match result {
            Ok((source, target, options)) => {
                let revision = target.revision();
                run(
                    AnalysisJob::InversePcr {
                        molecule: source.clone(),
                        target,
                        options,
                    },
                    source,
                    Some(revision),
                );
            }
            Err(e) => status.set(e),
        }
    };
    view! {
        <style>".mutation-panel .mutation-grid{display:grid;grid-template-columns:1fr 1fr;gap:8px}.mutation-panel pre{white-space:pre-wrap;overflow-wrap:anywhere;font:11px/1.6 monospace;background:#f3f6f3;padding:9px;max-height:250px;overflow:auto}.mutation-panel .mutation-product{padding:10px 0;border-bottom:1px solid #e2e9e2;font-size:12px;line-height:1.6}.mutation-panel .mutation-status{font-size:12px;line-height:1.6;overflow-wrap:anywhere}.mutation-panel ul{padding-left:18px}.mutation-panel details{margin:10px 0}"</style>
        <div class="inspector-content mutation-panel"><h2>"Mutagenesis"</h2><label>"Source template"<select aria-label="Mutagenesis template" prop:value=move||template_index.get().to_string() on:change=move|e|template_index.set(event_target_value(&e).parse().unwrap_or(0))>{move||query_options.get().into_iter().map(|(i,name)|view!{<option value=i.to_string()>{name}</option>}).collect_view()}</select></label><p class="mutation-status" role="status">{move||status.get()}</p>
        <h3>"PCR with declared mismatches"</h3><p class="helper">"Enter actual ordered oligo domains 5′→3′, including intended substitutions. Binding mismatches are conditional rules, not predicted polymerase success."</p>
        {[("Forward domain","Mutagenesis forward domain",f_domain),("Forward 5′ tail","Mutagenesis forward tail",f_tail),("Reverse domain","Mutagenesis reverse domain",r_domain),("Reverse 5′ tail","Mutagenesis reverse tail",r_tail)].into_iter().map(|(label,aria,value)|view!{<label>{label}<input aria-label=aria prop:value=move||value.get() on:input=move|e|value.set(event_target_value(&e))/></label>}).collect_view()}
        <div class="mutation-grid"><label>"Allowed disjoint mismatches"<input aria-label="Mutagenesis allowed mismatches" type="number" min="0" prop:value=move||mismatches.get() on:input=move|e|mismatches.set(event_target_value(&e))/></label><label>"Required exact 3′ bases"<input aria-label="Mutagenesis exact 3prime bases" type="number" min="0" prop:value=move||exact.get() on:input=move|e|exact.set(event_target_value(&e))/></label></div><button class="button primary full" disabled=move||busy.get()||source.get().is_none() on:click=simulate>"Simulate mismatch PCR"</button>
        {move||pcr.get().map(|bound|{let report=&bound.report;let download=bound.clone();view!{<section class="mutation-pcr-result"><h3>{format!("{} conditional PCR products",report.products.len())}</h3><p class="helper">{format!("Search complete: {} · complete model: {} · unsupported geometries {} · conflicting overlaps {}",report.search_complete,report.enumeration_complete,report.unsupported_geometries,report.conflicting_overlaps)}</p><ul>{report.warnings.iter().chain(bound.findings.iter()).cloned().map(|w|view!{<li>{w}</li>}).collect_view()}</ul>{report.products.iter().take(50).enumerate().map(|(index,p)|{let action=bound.clone();view!{<article class="mutation-product"><strong>{format!("Product {} · {} bp · {} primer substitutions",index+1,p.sequence.len(),p.substitutions.len())}</strong><p>{format!("Source starts F {} / R {} · domains overlap: {} · concrete: {}",p.forward_start+1,p.reverse_start+1,p.overlapping_domains,p.certain)}</p><pre>{p.sequence.as_str().chars().take(160).collect::<String>()}</pre><button class="button secondary" on:click=move|_|match preview_document(&action,index).and_then(|m|on_products.run(vec![m])){Ok(())=>status.set("Reference preview opened with source annotation mapping. Changed or uncertain features require review.".into()),Err(e)=>status.set(e)}>"Open product reference preview"</button></article>}}).collect_view()}<p class="helper">"First 50 geometries shown. Previews retain source annotation mappings and binding evidence; changed or uncertain features require review. They assert neither target verification nor physical chemistry. JSON retains all returned products and limits."</p><button class="button secondary full" on:click=move|_|match serde_json::to_vec_pretty(&download.report){Ok(bytes)=>if let Err(e)=crate::app::download("mismatch-pcr.json","application/json",&bytes){status.set(e)},Err(e)=>status.set(e.to_string())}>"Download mismatch PCR evidence"</button></section>}})}
        <h3>"Inverse-PCR replacement"</h3><p class="helper">{move||format!("Active document is the independently specified intended circular target: {}. Choose the original circular source above.",molecule.get().name)}</p><p class="helper">"Edit coordinates are one-based boundaries: start inclusive, end exclusive. Equal values insert without deletion. A circular end can exceed the source length to cross the origin. Example: start 4, end 6 replaces source bases 4–5."</p>
        <div class="mutation-grid">{[("Edit start boundary","Mutation edit start",edit_start),("Edit end boundary","Mutation edit end",edit_end),("Forward annealing length","Mutation forward length",f_length),("Reverse annealing length","Mutation reverse length",r_length),("Inserted bases on reverse tail","Mutation tail split",split)].into_iter().map(|(label,aria,value)|view!{<label>{label}<input aria-label=aria type="number" prop:value=move||value.get() on:input=move|e|value.set(event_target_value(&e))/></label>}).collect_view()}</div><label>"Replacement/inserted sequence"<textarea aria-label="Mutation inserted sequence" rows="3" prop:value=move||inserted.get() on:input=move|e|inserted.set(event_target_value(&e))></textarea></label><button class="button primary full" disabled=move||busy.get()||source.get().is_none() on:click=design>"Design inverse-PCR candidate"</button>
        {move||inverse.get().map(|report|{let candidate=Candidate::new(report.plan.clone());let download=report.clone();view!{<section class="inverse-pcr-result"><h3>{format!("Inverse-PCR check: {:?}",report.verification.status)}</h3><pre>{format!("F  {}\nR  {}",report.forward.ordered_sequence(),report.reverse.ordered_sequence())}</pre><p class="helper">{format!("Retained template {} bp · intended target {}…",report.retained_template_length,report.target_revision.chars().take(12).collect::<String>())}</p><ol><li>"Exact back-to-back PCR with explicit tails"</li><li>"Explicit terminal phosphorylation"</li><li>"Sealed blunt circular ligation"</li></ol><ul>{report.verification.findings.iter().chain(report.warnings.iter()).cloned().map(|w|view!{<li>{w}</li>}).collect_view()}</ul><button class="button primary full" on:click=move|_|on_candidate.run(candidate.clone())>"Review inverse-PCR candidate"</button><button class="button secondary full" on:click=move|_|match serde_json::to_vec_pretty(download.as_ref()){Ok(bytes)=>if let Err(e)=crate::app::download("inverse-pcr-design.json","application/json",&bytes){status.set(e)},Err(e)=>status.set(e.to_string())}>"Download inverse-PCR evidence"</button></section>}})}
        </div>
    }
}
