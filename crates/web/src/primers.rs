//! Worker-backed oligo evidence and anchored candidate selection.
use crate::worker;
use collagene_core::{
    Alphabet, Feature, Molecule,
    analysis::Primer,
    primers::{BindingOptions, DesignOptions, DesignReport, PrimerAnalysis, ThermoConditions},
    protocol::{AnalysisJob, AnalysisResult, PrimerEvidence, Request, Response},
};
use leptos::prelude::*;
use std::sync::{
    Arc,
    atomic::{AtomicBool, Ordering},
};
use wasm_bindgen_futures::spawn_local;

fn parse<T: std::str::FromStr>(value: String, label: &str) -> Result<T, String> {
    value
        .trim()
        .parse()
        .map_err(|_| format!("{label} is not a valid number"))
}
fn annotation(
    primer: &Primer,
    location: collagene_core::Location,
    metrics: &PrimerAnalysis,
    review: bool,
    source_revision: &str,
) -> Feature {
    let mut qualifiers=vec![("label".into(),Some(primer.name.clone())),("oligo_sequence".into(),Some(primer.ordered_sequence())),("annealing_sequence".into(),Some(primer.annealing.as_str().into())),("5prime_tail".into(),Some(primer.tail.clone())),("note".into(),Some("Primer annotation from reference-sequence evidence; amplification and specificity are not established.".into()))];
    qualifiers.push(("source_revision".into(), Some(source_revision.into())));
    if let Some(tm) = &metrics.tm {
        qualifiers.push(("note".into(),Some(format!("{}; Tm {:.3} C; sodium {} M; total strands {} M. Perfect annealing-domain duplex model.",tm.model,tm.tm_c,tm.conditions.sodium_molar,tm.conditions.total_strand_molar))));
    }
    Feature {
        kind: "primer_bind".into(),
        location,
        qualifiers,
        needs_review: review,
    }
}
#[component]
pub fn PrimerPanel(
    molecule: Memo<Molecule>,
    selection: RwSignal<Option<(usize, usize)>>,
    selected_feature: RwSignal<Option<usize>>,
    on_use: Callback<(Primer, Primer)>,
    on_annotate: Callback<(String, Feature), Result<(), String>>,
    #[prop(optional_no_strip)] seed: Option<String>,
) -> impl IntoView {
    let name = RwSignal::new("primer".to_string());
    let annealing = RwSignal::new(seed.unwrap_or_default());
    let tail = RwSignal::new(String::new());
    let sodium = RwSignal::new("50".to_string());
    let concentration = RwSignal::new("500".to_string());
    let mismatches = RwSignal::new("0".to_string());
    let three_prime = RwSignal::new("0".to_string());
    let start = RwSignal::new("1".to_string());
    let end = RwSignal::new(molecule.get_untracked().sequence.len().to_string());
    let min_length = RwSignal::new("18".to_string());
    let max_length = RwSignal::new("30".to_string());
    let target_tm = RwSignal::new("60".to_string());
    let max_delta = RwSignal::new("5".to_string());
    let min_gc = RwSignal::new("20".to_string());
    let max_gc = RwSignal::new("80".to_string());
    let max_homopolymer = RwSignal::new("5".to_string());
    let max_complementarity = RwSignal::new("3".to_string());
    let forward_tail = RwSignal::new(String::new());
    let reverse_tail = RwSignal::new(String::new());
    let reuse_index = RwSignal::new(String::new());
    let reuse_notice = RwSignal::new(String::new());
    let saved_primers =
        RwSignal::new(None::<Result<Arc<collagene_core::primer_reuse::List>, String>>);
    Effect::new(move |_| {
        let result = molecule
            .with(|source| collagene_core::primer_reuse::list(source, 1000))
            .map(Arc::new)
            .map_err(|error| error.to_string());
        saved_primers.set(Some(result));
        reuse_index.set(String::new());
        reuse_notice.set(String::new());
    });
    let kept_forward = RwSignal::new(None::<Primer>);
    let kept_reverse = RwSignal::new(None::<Primer>);
    let evidence = RwSignal::new(None::<PrimerEvidence>);
    let design = RwSignal::new(None::<DesignReport>);
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
            name.get(),
            annealing.get(),
            tail.get(),
            sodium.get(),
            concentration.get(),
            mismatches.get(),
            three_prime.get(),
            start.get(),
            end.get(),
            min_length.get(),
            max_length.get(),
            target_tm.get(),
            max_delta.get(),
            min_gc.get(),
            max_gc.get(),
            max_homopolymer.get(),
            max_complementarity.get(),
            forward_tail.get(),
            reverse_tail.get(),
        );
        generation.update(|g| *g = g.wrapping_add(1));
        evidence.set(None);
        design.set(None);
        status.set(String::new());
        busy.set(false);
    });
    let conditions = move || -> Result<ThermoConditions, String> {
        Ok(ThermoConditions {
            sodium_molar: parse::<f64>(sodium.get_untracked(), "Sodium")? / 1000.,
            total_strand_molar: parse::<f64>(concentration.get_untracked(), "Total strands")?
                * 1e-9,
        })
    };
    let run = {
        let alive = alive.clone();
        move |job: AnalysisJob, revision: String| {
            let token = generation.get_untracked().wrapping_add(1);
            generation.set(token);
            busy.set(true);
            evidence.set(None);
            design.set(None);
            status.set("Calculating in background…".into());
            let alive = alive.clone();
            spawn_local(async move {
                let response = worker::analyze(Request {
                    id: token,
                    revision: revision.clone(),
                    job,
                })
                .await;
                if !alive.load(Ordering::Relaxed)
                    || generation.get_untracked() != token
                    || molecule.get_untracked().revision() != revision
                {
                    return;
                }
                busy.set(false);
                match response {
                    Ok(Response {
                        id,
                        revision: returned,
                        result,
                    }) if id == token && returned == revision => {
                        match result {
                            AnalysisResult::PrimerEvidence(report) => {
                                status.set(format!("{} reference binding hit(s){}; inspect alternatives before PCR.",report.bindings.hits.len(),if report.bindings.truncated{" · incomplete enumeration"}else{""}));
                                evidence.set(Some(report));
                            }
                            AnalysisResult::PrimerDesign(report) => {
                                status.set(format!("{} ranked pair(s); binding uniqueness has not been established.",report.pairs.len()));
                                design.set(Some(report));
                            }
                            AnalysisResult::Error(error) => status.set(error),
                            _ => status.set("Unexpected primer worker response".into()),
                        }
                    }
                    Ok(_) => {
                        status.set("Worker response did not match this source revision".into())
                    }
                    Err(error) => status.set(error),
                }
            });
        }
    };
    let analyze = {
        let run = run.clone();
        move |_| {
            let request = (|| -> Result<AnalysisJob, String> {
                let primer = Primer::new(
                    name.get_untracked(),
                    &annealing.get_untracked(),
                    &tail.get_untracked(),
                )
                .map_err(|e| e.to_string())?;
                Ok(AnalysisJob::AnalyzePrimer {
                    molecule: molecule.get_untracked(),
                    primer,
                    conditions: conditions()?,
                    binding_options: BindingOptions {
                        max_mismatches: parse(mismatches.get_untracked(), "Allowed mismatches")?,
                        min_three_prime_exact: parse(
                            three_prime.get_untracked(),
                            "Required 3′ matches",
                        )?,
                        ..BindingOptions::default()
                    },
                })
            })();
            match request {
                Ok(job) => run(job, molecule.get_untracked().revision()),
                Err(error) => status.set(error),
            }
        }
    };
    let reuse = move |_| {
        let result = (|| -> Result<collagene_core::primer_reuse::Resolved, String> {
            let index = reuse_index
                .get_untracked()
                .parse::<usize>()
                .map_err(|_| "Select a reusable saved primer annotation".to_string())?;
            let list = saved_primers
                .get_untracked()
                .ok_or("Saved primer list is not ready")??;
            molecule
                .with_untracked(|source| {
                    collagene_core::primer_reuse::resolve_at_revision(
                        source,
                        index,
                        &list.source_revision,
                    )
                })
                .map_err(|error| error.to_string())
        })();
        match result {
            Ok(resolved) => {
                let primer = resolved.primer;
                name.set(primer.name);
                annealing.set(primer.annealing.as_str().into());
                tail.set(primer.tail);
                reuse_notice.set(format!("Saved oligo loaded from feature {} at exact source revision {}. Analyze binding against the current template before use. {}",
                    resolved.feature_index,resolved.source_revision,resolved.warnings.iter().take(20).cloned().collect::<Vec<_>>().join(" ")));
            }
            Err(error) => status.set(error),
        }
    };
    let select_design = move |_| {
        let request = (|| -> Result<AnalysisJob, String> {
            Ok(AnalysisJob::DesignPrimers {
                molecule: molecule.get_untracked(),
                options: DesignOptions {
                    target_start: parse::<usize>(start.get_untracked(), "Target start")?
                        .checked_sub(1)
                        .ok_or("Target starts at base 1 or greater")?,
                    target_end: parse(end.get_untracked(), "Target end")?,
                    min_length: parse(min_length.get_untracked(), "Minimum length")?,
                    max_length: parse(max_length.get_untracked(), "Maximum length")?,
                    target_tm_c: parse(target_tm.get_untracked(), "Target Tm")?,
                    max_tm_difference_c: parse(max_delta.get_untracked(), "Pair Tm difference")?,
                    min_gc_fraction: parse::<f64>(min_gc.get_untracked(), "Minimum GC")? / 100.,
                    max_gc_fraction: parse::<f64>(max_gc.get_untracked(), "Maximum GC")? / 100.,
                    max_homopolymer: parse(max_homopolymer.get_untracked(), "Homopolymer limit")?,
                    max_three_prime_complementarity: parse(
                        max_complementarity.get_untracked(),
                        "3′ complementarity limit",
                    )?,
                    forward_tail: forward_tail.get_untracked(),
                    reverse_tail: reverse_tail.get_untracked(),
                    conditions: conditions()?,
                    ..DesignOptions::default()
                },
            })
        })();
        match request {
            Ok(job) => run(job, molecule.get_untracked().revision()),
            Err(error) => status.set(error),
        }
    };
    view! {
        <style>".primer-panel .primer-grid{display:grid;grid-template-columns:1fr 1fr;gap:8px}.primer-panel .primer-grid label{margin:0}.primer-panel pre{font:11px/1.6 monospace;white-space:pre-wrap;overflow-wrap:anywhere;background:#f2f6f2;padding:8px;max-height:220px;overflow:auto}.primer-panel .primer-hit,.primer-panel .primer-pair{padding:10px 0;border-bottom:1px solid #e1e7e0;font-size:12px;line-height:1.6}.primer-panel .primer-actions{display:flex;gap:6px;flex-wrap:wrap}.primer-panel .primer-actions button{font-size:11px;padding:6px 8px}.primer-panel ul{padding-left:18px}.primer-status{font-size:12px;line-height:1.6;overflow-wrap:anywhere}"</style>
        <div class="inspector-content primer-panel"><h2>"Primers & oligos"</h2><p class="helper">"All oligos are written 5′→3′. Initial binding and Tm use the annealing domain; structure screens include tails."</p>
        <details open><summary>"Calculation conditions"</summary><div class="primer-grid"><label>"Sodium (mM)"<input aria-label="Primer sodium mM" type="number" prop:value=move||sodium.get() on:input=move|e|sodium.set(event_target_value(&e))/></label><label>"Total strands (nM)"<input aria-label="Primer total strands nM" type="number" prop:value=move||concentration.get() on:input=move|e|concentration.set(event_target_value(&e))/></label></div><p class="helper">"Perfect DNA-duplex model; magnesium, dNTPs and additives are not inferred. Tm is not a recommended PCR temperature."</p></details>
        <h3>"Analyze one primer"</h3><label>"Saved primer annotations"<select aria-label="Saved primer annotations" prop:value=move||reuse_index.get() on:change=move|e|reuse_index.set(event_target_value(&e))><option value="">"Choose a saved primer…"</option>{move||saved_primers.get().and_then(Result::ok).map(|list|list.primers.iter().map(|entry|{let title=entry.diagnostic.clone().unwrap_or_default();view!{<option value=entry.feature_index.to_string() disabled=!entry.reusable title=title>{format!("{}{}",entry.name,if entry.reusable{""}else{" · cannot reuse"})}</option>}}).collect_view())}</select></label><button class="button secondary full" on:click=reuse>"Load saved primer"</button><p class="helper">{move||match saved_primers.get(){Some(Ok(list)) if list.omitted>0=>format!("Showing {} of {} saved primers; {} omitted by the bounded list.",list.primers.len(),list.total,list.omitted),Some(Err(error))=>error,_=>String::new()}}</p><p class="helper primer-reuse-notice" role="status">{move||reuse_notice.get()}</p><label>"Name"<input aria-label="Primer name" prop:value=move||name.get() on:input=move|e|name.set(event_target_value(&e))/></label><label>"Annealing domain"<input aria-label="Primer annealing domain" prop:value=move||annealing.get() on:input=move|e|annealing.set(event_target_value(&e))/></label><label>"5′ tail"<input aria-label="Primer tail" prop:value=move||tail.get() on:input=move|e|tail.set(event_target_value(&e))/></label>
        <button class="button secondary full" disabled=move||selection.get().is_none() on:click=move|_|{if let Some((a,b))=selection.get_untracked(){if let Some(text)=molecule.get_untracked().sequence.as_str().get(a..b){annealing.set(text.to_string());}}}>"Use selected reference bases"</button>
        <div class="primer-grid"><label>"Allowed mismatches"<input aria-label="Primer allowed mismatches" type="number" min="0" prop:value=move||mismatches.get() on:input=move|e|mismatches.set(event_target_value(&e))/></label><label>"Required exact 3′ bases"<input aria-label="Primer exact 3prime bases" type="number" min="0" prop:value=move||three_prime.get() on:input=move|e|three_prime.set(event_target_value(&e))/></label></div>
        <button class="button primary full" disabled=move||busy.get()||molecule.get().sequence.alphabet()!=Alphabet::Dna on:click=analyze>"Analyze primer & binding"</button>
        <p class="primer-status" role="status">{move||status.get()}</p>
        {move||evidence.get().map(|report|{let revision=report.bindings.source_revision.clone();let metrics=report.analysis.clone();let primer=Primer::new(&metrics.name,&metrics.annealing,&metrics.tail).expect("worker validated primer");let for_download=report.clone();let keep_forward=primer.clone();let keep_reverse=primer.clone();view!{<section class="primer-evidence"><h3>"Oligo evidence"</h3><pre>{report.analysis.ordered_sequence.clone()}</pre><div class="primer-actions"><button class="button secondary" on:click=move|_|kept_forward.set(Some(keep_forward.clone()))>"Keep as forward primer"</button><button class="button secondary" on:click=move|_|kept_reverse.set(Some(keep_reverse.clone()))>"Keep as reverse primer"</button></div><p>{format!("Annealing {} nt · tail {} nt · GC {:.1}–{:.1}%",report.analysis.annealing.len(),report.analysis.tail.len(),report.analysis.annealing_gc.minimum_fraction*100.,report.analysis.annealing_gc.maximum_fraction*100.)}</p><p class="primer-tm">{report.analysis.tm.as_ref().map(|t|format!("Tm {:.2} °C · {}",t.tm_c,t.model)).unwrap_or_else(||report.analysis.tm_unavailable.clone().unwrap_or_else(||"Tm unavailable".into()))}</p><p>{format!("Contiguous possible runs: hairpin {} · homodimer {} · homodimer at a 3′ end {}",metrics.hairpin.max_run,metrics.homodimer.max_run,metrics.homodimer.max_three_prime_run)}</p><ul>{report.analysis.warnings.into_iter().chain(report.bindings.warnings.clone()).map(|w|view!{<li>{w}</li>}).collect_view()}</ul><h4>{format!("{} binding hits{}",report.bindings.hits.len(),if report.bindings.truncated{" · truncated"}else{""})}</h4>{report.bindings.hits.into_iter().take(100).map(|hit|{let spans=hit.location.spans();let first=spans.first().map(|(a,b,_)|(*a,*b));let saved=annotation(&primer,hit.location.clone(),&metrics,!hit.certain,&revision);let source=revision.clone();view!{<article class="primer-hit"><strong>{format!("{}–{}{}{}",hit.start+1,hit.start+hit.length,if hit.reverse{" · reverse"}else{" · forward"},if hit.crosses_origin{" · crosses origin"}else{""})}</strong><p>{format!("{} mismatches · {} ambiguous positions · {} exact 3′ bases",hit.mismatches.len(),hit.ambiguous_positions.len(),hit.three_prime_exact)}</p><div class="primer-actions"><button class="button secondary" on:click=move|_|{selected_feature.set(None);selection.set(first);status.set(if spans.len()>1{"Origin-spanning binding: first reference span selected; annotation preserves every span.".into()}else{"Binding interval selected.".into()});}>"Select binding"</button><button class="button secondary" on:click=move|_|match on_annotate.run((source.clone(),saved.clone())){Ok(())=>status.set("Primer annotation saved; binding source revision is retained in its qualifiers.".into()),Err(error)=>status.set(error)}>"Save primer annotation"</button></div></article>}}).collect_view()}<p class="helper">"Showing up to 100 hits; the JSON report retains all returned hits. Compatibility is not amplification or genome-wide specificity."</p><button class="button secondary full" on:click=move|_|match serde_json::to_vec_pretty(&for_download){Ok(bytes)=>if let Err(error)=crate::app::download("primer-evidence.json","application/json",&bytes){status.set(error)},Err(error)=>status.set(error.to_string())}>"Download primer evidence"</button></section>}})}
        <h3>"Reuse a primer pair"</h3><p class="helper">{move||format!("Forward: {} · reverse: {}",kept_forward.get().map(|p|p.name).unwrap_or_else(||"not chosen".into()),kept_reverse.get().map(|p|p.name).unwrap_or_else(||"not chosen".into()))}</p><p class="helper">"Keep analyzed or loaded oligos in either role, then independently simulate their products against the active template."</p><button class="button secondary full" disabled=move||kept_forward.get().is_none()||kept_reverse.get().is_none() on:click=move|_|if let (Some(f),Some(r))=(kept_forward.get_untracked(),kept_reverse.get_untracked()){on_use.run((f,r));}>"Use kept pair for PCR"</button>
        <h3>"Select an anchored primer pair"</h3><p class="helper">"One-based inclusive target coordinates. For a circle, end may exceed the document length to cross the origin, up to one turn."</p><button class="button secondary full" disabled=move||selection.get().is_none() on:click=move|_|if let Some((a,b))=selection.get_untracked(){start.set((a+1).to_string());end.set(b.to_string());}>"Use selection as target"</button>
        <div class="primer-grid">{[("Target start","Primer target start",start),("Target end","Primer target end",end),("Min length","Primer minimum length",min_length),("Max length","Primer maximum length",max_length),("Target Tm (°C)","Primer target Tm",target_tm),("Pair ΔTm (°C)","Primer pair Tm difference",max_delta),("Min GC (%)","Primer minimum GC",min_gc),("Max GC (%)","Primer maximum GC",max_gc),("Homopolymer limit","Primer homopolymer limit",max_homopolymer),("3′ run limit","Primer complementarity limit",max_complementarity)].into_iter().map(|(label,aria,value)|view!{<label>{label}<input aria-label=aria type="number" prop:value=move||value.get() on:input=move|e|value.set(event_target_value(&e))/></label>}).collect_view()}</div>
        <label>"Forward 5′ tail"<input aria-label="Design forward tail" prop:value=move||forward_tail.get() on:input=move|e|forward_tail.set(event_target_value(&e))/></label><label>"Reverse 5′ tail"<input aria-label="Design reverse tail" prop:value=move||reverse_tail.get() on:input=move|e|reverse_tail.set(event_target_value(&e))/></label>
        <button class="button primary full" disabled=move||busy.get()||molecule.get().sequence.alphabet()!=Alphabet::Dna on:click=select_design>"Select primer pairs"</button>
        {move||design.get().map(|report|{let source=report.source_revision.clone();let for_download=report.clone();view!{<section class="primer-design"><h3>{format!("{} ranked primer pairs",report.pairs.len())}</h3><ul>{report.warnings.into_iter().map(|w|view!{<li>{w}</li>}).collect_view()}</ul>{report.pairs.into_iter().enumerate().map(|(i,pair)|{let f=pair.forward.primer.clone();let r=pair.reverse.primer.clone();let save_f=annotation(&f,pair.forward.location.clone(),&pair.forward.analysis,false,&source);let save_r=annotation(&r,pair.reverse.location.clone(),&pair.reverse.analysis,false,&source);let source_f=source.clone();let source_r=source.clone();view!{<article class="primer-pair"><strong>{format!("Pair {} · {} bp including tails",i+1,pair.product_length)}</strong><pre>{format!("F  {}\nR  {}",f.ordered_sequence(),r.ordered_sequence())}</pre><p>{format!("Tm {:.2} / {:.2} °C · Δ {:.2} °C · score {:.3}",pair.forward.analysis.tm.as_ref().unwrap().tm_c,pair.reverse.analysis.tm.as_ref().unwrap().tm_c,pair.tm_difference_c,pair.score)}</p><div class="primer-actions"><button class="button primary" on:click=move|_|on_use.run((f.clone(),r.clone()))>"Use pair for PCR"</button><button class="button secondary" on:click=move|_|match on_annotate.run((source_f.clone(),save_f.clone())){Ok(())=>status.set("Forward primer annotation saved. Re-run selection against the new document revision to save another.".into()),Err(error)=>status.set(error)}>"Save forward primer"</button><button class="button secondary" on:click=move|_|match on_annotate.run((source_r.clone(),save_r.clone())){Ok(())=>status.set("Reverse primer annotation saved. Re-run selection against the new document revision to save another.".into()),Err(error)=>status.set(error)}>"Save reverse primer"</button></div></article>}}).collect_view()}<button class="button secondary full" on:click=move|_|match serde_json::to_vec_pretty(&for_download){Ok(bytes)=>if let Err(error)=crate::app::download("primer-pairs.json","application/json",&bytes){status.set(error)},Err(error)=>status.set(error.to_string())}>"Download primer pairs"</button></section>}})}
        </div>
    }
}
