//! Read-only sequence screening and explicit, review-required candidate annotations.
use crate::{
    protein_liability_model::{Capture, ProvidedContext, interval},
    source_identity::SourceInfo,
    worker::CachedWorker,
};
use collagene_core::{
    Molecule, protein_liabilities as native,
    protocol::{AnalysisJob, AnalysisResult, Request},
    run_lineage::EventBinding,
};
use leptos::prelude::*;
use std::{
    rc::Rc,
    sync::{
        Arc,
        atomic::{AtomicBool, Ordering},
    },
};
fn certainty(c: &native::Certainty) -> &'static str {
    match c {
        native::Certainty::CertainPattern => "Certain pattern",
        native::Certainty::PossiblePattern => "Possible pattern",
    }
}
fn cdr_label(c: &native::CdrMembership, o: &native::Options) -> String {
    match c {
        native::CdrMembership::Unassessed => "CDR context unassessed".into(),
        native::CdrMembership::OutsideProvided => "Outside supplied regions".into(),
        native::CdrMembership::Inside { region_index } => o
            .provided_cdrs
            .as_ref()
            .and_then(|c| c.regions.get(*region_index))
            .map(|r| {
                format!(
                    "{} · {}",
                    r.domain_id,
                    match r.kind {
                        native::CdrKind::Cdr1 => "CDR1",
                        native::CdrKind::Cdr2 => "CDR2",
                        native::CdrKind::Cdr3 => "CDR3",
                    }
                )
            })
            .unwrap_or_else(|| "Invalid supplied region".into()),
    }
}
#[component]
pub fn ProteinLiabilityPanel(
    source: Memo<SourceInfo>,
    epoch: RwSignal<u64>,
    workspace_busy: RwSignal<bool>,
    prepare: Callback<(EventBinding, native::Options), Result<Molecule, String>>,
    on_select: Callback<(Capture, native::Hit), Result<(), String>>,
    on_annotate: Callback<(Capture, native::Hit), Result<(), String>>,
    selection: Memo<Option<(usize, usize)>>,
) -> impl IntoView {
    let maximum = RwSignal::new("1000".to_string());
    let context = RwSignal::new(None::<ProvidedContext>);
    let options = Memo::new(move |_| -> Result<native::Options, String> {
        let max_hits = maximum
            .get()
            .trim()
            .parse::<usize>()
            .map_err(|_| "Enter a whole retained-hit limit.".to_owned())?;
        if max_hits > native::MAX_HITS {
            return Err(format!(
                "Retain at most {} hits; zero requests counts only.",
                native::MAX_HITS
            ));
        }
        if let Some(c) = context.get() {
            if c.context.definition_label.trim().is_empty()
                || c.context.provenance.trim().is_empty()
            {
                return Err("Describe the supplied CDR definition and provenance.".into());
            }
            c.options(&source.get().binding, max_hits)
        } else {
            Ok(native::Options {
                max_hits,
                provided_cdrs: None,
            })
        }
    });
    let status = RwSignal::new(String::new());
    let busy = RwSignal::new(false);
    let generation = RwSignal::new(0u32);
    let report = RwSignal::new(None::<(Capture, Arc<native::Report>)>);
    let page = RwSignal::new(0usize);
    let rule_filter = RwSignal::new(None::<native::RuleId>);
    let possible_only = RwSignal::new(false);
    let indices = Memo::new(move |_| {
        let selected = rule_filter.get();
        let possible = possible_only.get();
        report.with(|r| {
            r.as_ref()
                .map(|(_, r)| {
                    r.hits
                        .iter()
                        .enumerate()
                        .filter(|(_, h)| {
                            selected.as_ref().is_none_or(|id| *id == h.rule)
                                && (!possible || h.certainty == native::Certainty::PossiblePattern)
                        })
                        .map(|(i, _)| i)
                        .collect::<Vec<_>>()
                })
                .unwrap_or_default()
        })
    });
    Effect::new(move |_| {
        indices.with(|_| ());
        page.set(0);
    });
    let worker = StoredValue::new_local(None::<Rc<CachedWorker>>);
    let alive = Arc::new(AtomicBool::new(true));
    Effect::new(move |_| {
        let _ = (
            source.get(),
            options.get(),
            epoch.get(),
            workspace_busy.get(),
        );
        generation.update(|n| *n = n.wrapping_add(1));
        if let Some(w) = worker.get_value() {
            w.close();
        }
        worker.set_value(None);
        busy.set(false);
        report.set(None);
    });
    {
        let alive = alive.clone();
        on_cleanup(move || {
            alive.store(false, Ordering::Relaxed);
            if let Some(w) = worker.get_value() {
                w.close();
            }
        });
    }
    let run = {
        let alive = alive.clone();
        move |_| {
            if busy.get_untracked() || workspace_busy.get_untracked() {
                return;
            }
            let prepared = (|| {
                let o = options.get_untracked()?;
                let s = source.get_untracked();
                let m = prepare.run((s.binding.clone(), o.clone()))?;
                let capture = Capture::prepare(epoch.get_untracked(), s.binding, &m, o)?;
                Ok::<_, String>((m, capture))
            })();
            let (m, c) = match prepared {
                Ok(v) => v,
                Err(e) => {
                    status.set(e);
                    return;
                }
            };
            let w = match CachedWorker::new() {
                Ok(w) => Rc::new(w),
                Err(e) => {
                    status.set(e);
                    return;
                }
            };
            if let Some(old) = worker.get_value() {
                old.close();
            }
            worker.set_value(Some(w.clone()));
            let token = generation.get_untracked().wrapping_add(1);
            generation.set(token);
            report.set(None);
            busy.set(true);
            status.set("Screening sequence patterns…".into());
            let alive = alive.clone();
            wasm_bindgen_futures::spawn_local(async move {
                let result = w
                    .ask(Request {
                        id: token,
                        revision: c.context_revision.clone(),
                        job: AnalysisJob::ScreenProteinLiabilities {
                            molecule: m,
                            options: c.options.clone(),
                        },
                    })
                    .await;
                w.close();
                if !alive.load(Ordering::Relaxed) || generation.get_untracked() != token {
                    return;
                }
                busy.set(false);
                worker.set_value(None);
                let result = (|| {
                    c.guard(
                        epoch.get_untracked(),
                        &source.get_untracked().binding,
                        &options.get_untracked()?,
                    )?;
                    let r = result?;
                    match r.result {
                        AnalysisResult::ProteinLiabilities(value) => {
                            c.validate_report(&value)?;
                            Ok(value)
                        }
                        AnalysisResult::Error(e) => Err(e),
                        _ => Err("Unexpected protein liability response.".into()),
                    }
                })();
                match result {
                    Ok(r) => {
                        status.set("Full-source counts ready. Matches are candidates, not measured modifications.".into());
                        report.set(Some((c, Arc::new(r))));
                    }
                    Err(e) => status.set(e),
                }
            });
        }
    };
    let domain = RwSignal::new(String::new());
    let kind = RwSignal::new(native::CdrKind::Cdr1);
    let start = RwSignal::new("1".to_string());
    let end = RwSignal::new("1".to_string());
    let add_region = move |_| {
        let result = (|| {
            let s = source.get_untracked();
            let mut c = context
                .get_untracked()
                .ok_or("Start an explicit supplied-CDR context first.")?;
            let (a, b) = interval(&start.get_untracked(), &end.get_untracked(), s.length)?;
            c.context.regions.push(native::CdrRegion {
                domain_id: domain.get_untracked(),
                kind: kind.get_untracked(),
                start: a,
                end: b,
            });
            let o = c.options(&s.binding, 0)?;
            prepare.run((s.binding, o))?;
            context.set(Some(c));
            Ok::<_, String>(())
        })();
        status.set(
            result
                .map(|_| {
                    "Supplied region added. Its definition remains documentary context.".into()
                })
                .unwrap_or_else(|e| e),
        );
    };
    view! {<section class="inspector-content protein-liabilities" aria-label="Protein sequence liabilities" style="overflow-wrap:anywhere"><details><summary>"Sequence liability screening"</summary>
        <p class="helper">"Literal sequence candidates for review. This screen does not predict measured PTMs, risk probabilities or activity. Automatic antibody domain detection and numbering are not yet available."</p>
        <fieldset disabled=move||workspace_busy.get() style="border:0;padding:0;min-width:0">
        <label>"Retained hit limit"<input aria-label="Protein liability retained hit limit" inputmode="numeric" prop:value=move||maximum.get() on:input=move|e|maximum.set(event_target_value(&e))/></label><p class="helper">"Zero requests counts only. Every source position is scanned even after the retained-hit cap."</p>
        <details><summary>"Optional supplied CDR context"</summary><p class="helper">"Use one-based inclusive coordinates under your explicitly stated definition. These regions are not validated Kabat, Chothia or IMGT numbering and do not prove complete domain coverage."</p>
            <button class="button secondary full" on:click=move|_|{let s=source.get_untracked();context.set(Some(ProvidedContext{source:s.binding.clone(),context:native::ProvidedCdrs{source_revision:s.binding.revision,definition_label:String::new(),provenance:String::new(),regions:vec![]}}));status.set("New supplied-CDR context: enter the definition and provenance before adding regions.".into());}>"Start or reset CDR context for this source"</button>
            <button class="text-tool" disabled=move||context.with(Option::is_none) on:click=move|_|context.set(None)>"Clear supplied CDR context"</button>
            <Show when=move||context.with(Option::is_some)>
            <label>"CDR definition description"<input aria-label="Supplied CDR definition" maxlength="128" prop:value=move||context.with(|c|c.as_ref().map(|c|c.context.definition_label.clone()).unwrap_or_default()) on:input=move|e|context.update(|c|if let Some(c)=c{c.context.definition_label=event_target_value(&e);})/></label>
            <label>"CDR provenance"<input aria-label="Supplied CDR provenance" maxlength="1024" prop:value=move||context.with(|c|c.as_ref().map(|c|c.context.provenance.clone()).unwrap_or_default()) on:input=move|e|context.update(|c|if let Some(c)=c{c.context.provenance=event_target_value(&e);})/></label>
            <label>"Domain label"<input aria-label="Supplied CDR domain label" maxlength="64" prop:value=move||domain.get() on:input=move|e|domain.set(event_target_value(&e))/></label>
            <label>"CDR"<select aria-label="Supplied CDR kind" prop:value=move||match kind.get(){native::CdrKind::Cdr1=>"1",native::CdrKind::Cdr2=>"2",native::CdrKind::Cdr3=>"3"} on:change=move|e|match event_target_value(&e).as_str(){"1"=>kind.set(native::CdrKind::Cdr1),"2"=>kind.set(native::CdrKind::Cdr2),"3"=>kind.set(native::CdrKind::Cdr3),_=>status.set("Unsupported CDR label.".into())}><option value="1">"CDR1"</option><option value="2">"CDR2"</option><option value="3">"CDR3"</option></select></label>
            <label>"Start residue"<input aria-label="Supplied CDR start" inputmode="numeric" prop:value=move||start.get() on:input=move|e|start.set(event_target_value(&e))/></label><label>"End residue (inclusive)"<input aria-label="Supplied CDR end" inputmode="numeric" prop:value=move||end.get() on:input=move|e|end.set(event_target_value(&e))/></label>
            <button class="text-tool" disabled=move||selection.get().is_none() on:click=move|_|if let Some((a,b))=selection.get_untracked(){start.set((a+1).to_string());end.set(b.to_string());}>"Use current selection"</button><button class="button secondary full" on:click=add_region>"Add supplied CDR"</button>
            {move||context.with(|c|c.as_ref().map(|c|c.context.regions.clone()).unwrap_or_default()).into_iter().enumerate().map(|(i,r)|view!{<p class="helper">{format!("{} · {} · {}–{}",r.domain_id,match r.kind{native::CdrKind::Cdr1=>"CDR1",native::CdrKind::Cdr2=>"CDR2",native::CdrKind::Cdr3=>"CDR3"},r.start+1,r.end)}<button class="text-tool" on:click=move|_|context.update(|c|if let Some(c)=c{if i<c.context.regions.len(){c.context.regions.remove(i);}})>"Remove region"</button></p>}).collect_view()}
            </Show>
        </details>
        {move||options.get().err().map(|e|view!{<p class="helper" role="alert">{e}</p>})}
        <button class="button primary full" disabled=move||busy.get()||options.get().is_err() on:click=run.clone()>"Screen protein liabilities"</button><Show when=move||busy.get()><button class="text-tool" on:click=move|_|{generation.update(|n|*n=n.wrapping_add(1));if let Some(w)=worker.get_value(){w.close();}worker.set_value(None);busy.set(false);status.set("Protein screening cancelled.".into());}>"Cancel screening"</button></Show>
        </fieldset>
        {move||report.with(|r|r.as_ref().map(|(_,r)|view!{<div><h3>"Full-source pattern counts"</h3><table style="width:100%;font-size:12px"><thead><tr><th>"Rule"</th><th>"Certain"</th><th>"Possible"</th></tr></thead><tbody>{r.counts.iter().map(|c|view!{<tr><th scope="row">{c.rule.label()}</th><td>{c.certain}</td><td>{c.possible}</td></tr>}).collect_view()}</tbody></table><p class="helper">{format!("{} source positions · {} ambiguous residues · {} stop symbols · {} retained hits{}",r.source_length,r.ambiguous_residues,r.stop_symbols,r.hits.len(),if r.hits_truncated{" (hit list capped; counts remain complete)"}else{""})}</p><p class="helper">{if r.cdr_rule_unassessed{"The supplied-CDR rule is unassessed: no regions were supplied."}else{"CDR membership refers only to supplied regions; other CDRs/domains remain unassessed."}}</p></div>}))}
        <Show when=move||report.with(Option::is_some)>
            <label>"Visible rule"<select aria-label="Protein liability rule filter" prop:value=move||rule_filter.with(|id|id.as_ref().map(|id|native::rules().iter().position(|r|r.id==*id).expect("current rule").to_string()).unwrap_or_else(||"all".into())) on:change=move|e|{let v=event_target_value(&e);if v=="all"{rule_filter.set(None);}else if let Some(r)=v.parse::<usize>().ok().and_then(|i|native::rules().get(i).cloned()){rule_filter.set(Some(r.id));}else{status.set("Unsupported rule filter.".into());}}><option value="all">"All rules"</option>{native::rules().into_iter().enumerate().map(|(i,r)|view!{<option value=i.to_string()>{r.id.label()}</option>}).collect_view()}</select></label>
            <label class="check-row"><input type="checkbox" prop:checked=move||possible_only.get() on:change=move|e|possible_only.set(event_target_checked(&e))/>"Possible matches only"</label>
            {move||{let rows=indices.with(|rows|rows.iter().skip(page.get().saturating_mul(50)).take(50).copied().collect::<Vec<_>>());report.with(|r|r.as_ref().map(|(capture,r)|rows.into_iter().map(|i|{let h=r.hits[i].clone();let for_select=h.clone();let for_annotation=h.clone();let select_capture=capture.clone();let annotate_capture=capture.clone();view!{<div class="protein-liability-hit"><p class="helper">{format!("{}–{} · {} · {} · {}",h.start+1,h.end,h.residues,h.rule.label(),certainty(&h.certainty))}</p><p class="helper">{cdr_label(&h.cdr,&r.options)}</p><button class="text-tool" disabled=move||workspace_busy.get() on:click=move|_|{let result=options.get_untracked().and_then(|o|{select_capture.guard(epoch.get_untracked(),&source.get_untracked().binding,&o)?;on_select.run((select_capture.clone(),for_select.clone()))});if let Err(e)=result{status.set(e);}}>"Select match"</button><button class="text-tool" disabled=move||workspace_busy.get() on:click=move|_|{let result=options.get_untracked().and_then(|o|{annotate_capture.guard(epoch.get_untracked(),&source.get_untracked().binding,&o)?;on_annotate.run((annotate_capture.clone(),for_annotation.clone()))});status.set(result.map(|_|"Added a review-required sequence-pattern candidate, not a measured modification.".into()).unwrap_or_else(|e|e));}>"Annotate candidate"</button></div>}}).collect_view()))}}
            <p class="helper">{move||format!("{} retained hits match this filter · page {}",indices.with(Vec::len),page.get()+1)}</p><button class="text-tool" disabled=move||page.get()==0 on:click=move|_|page.update(|p|*p=p.saturating_sub(1))>"Previous matches"</button><button class="text-tool" disabled=move||page.get().saturating_add(1).saturating_mul(50)>=indices.with(Vec::len) on:click=move|_|page.update(|p|*p+=1)>"Next matches"</button>
            <button class="button secondary full" disabled=move||workspace_busy.get() on:click=move|_|{let result=(||{let (c,r)=report.get_untracked().ok_or("Screen the source first.")?;c.guard(epoch.get_untracked(),&source.get_untracked().binding,&options.get_untracked()?)?;serde_json::to_vec_pretty(r.as_ref()).map_err(|e|e.to_string())})();match result{Ok(bytes)=>crate::app::download_with_completion("protein-liabilities.json","application/json",&bytes,move|result|if status.try_get_untracked().is_some(){if let Err(e)=result{status.set(e);}}),Err(e)=>status.set(e)}}>"Download liability report"</button>
        </Show><p class="helper" role="status">{move||status.get()}</p>
    </details></section>}
}
