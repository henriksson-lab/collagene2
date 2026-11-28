//! Analysis of the focused row, distinct from an extracted sequence document.
use crate::biology::antibody::alignment_antibody_model::Capture;
use crate::biology::antibody::antibody_model;
use crate::platform::worker::CachedWorker;
use crate::viewer::alignment::alignment_view::AlignmentRevision;
use collagene_core::antibody_numbering::domain_tracks as t;
use collagene_core::{
    protein::antibody_alignment as a,
    protocol::{
        request::Request,
        schema::{AnalysisJob, AnalysisResult},
    },
};
use leptos::prelude::*;
use std::{
    rc::Rc,
    sync::{
        Arc,
        atomic::{AtomicBool, Ordering},
    },
};
#[component]
pub fn AlignmentAntibodyPanel(
    alignment: Memo<Option<AlignmentRevision>>,
    focused: RwSignal<(usize, usize)>,
    epoch: RwSignal<u64>,
    candidate: RwSignal<usize>,
    on_report: Callback<Option<(u64, Arc<a::Report>)>>,
) -> impl IntoView {
    let scheme = RwSignal::new(t::NumberingScheme::Imgt);
    let definition = RwSignal::new(t::RegionDefinition::Imgt);
    let choice = Memo::new(move |_| t::Selection {
        numbering_scheme: scheme.get(),
        region_definition: definition.get(),
    });
    let minimum = RwSignal::new("300".to_string());
    let substitutions = RwSignal::new("1".to_string());
    let loop_max = RwSignal::new("60".to_string());
    let retained = RwSignal::new("32".to_string());
    let settings = Memo::new(move |_| -> Result<_, String> {
        Ok(antibody_model::preset(
            minimum
                .get()
                .parse()
                .map_err(|_| "Enter an integer V score.")?,
            substitutions
                .get()
                .parse()
                .map_err(|_| "Enter a whole FR4 substitution count.")?,
            loop_max
                .get()
                .parse()
                .map_err(|_| "Enter a whole CDR3 length.")?,
            retained
                .get()
                .parse()
                .map_err(|_| "Enter a whole candidate limit.")?,
        ))
    });
    let pin = Memo::new(move |_| {
        alignment.with(|a| {
            a.as_ref().and_then(|a| {
                a.doc().rows.get(focused.get().0).map(|r| a::Source {
                    alignment_revision: a.revision.clone(),
                    row_id: r.id.clone(),
                })
            })
        })
    });
    let busy = RwSignal::new(false);
    let status = RwSignal::new(String::new());
    let report = RwSignal::new(None::<(Capture, Arc<a::Report>)>);
    let generation = RwSignal::new(0u32);
    let worker = StoredValue::new_local(None::<Rc<CachedWorker>>);
    let alive = Arc::new(AtomicBool::new(true));
    Effect::new(move |_| {
        let _ = (pin.get(), epoch.get(), settings.get(), choice.get());
        generation.update(|n| *n = n.wrapping_add(1));
        if let Some(w) = worker.get_value() {
            w.close();
        }
        worker.set_value(None);
        busy.set(false);
        report.set(None);
        candidate.set(0);
        on_report.run(None);
        status.set(String::new());
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
            if busy.get_untracked() {
                return;
            }
            let prepared = (|| {
                let source = pin
                    .get_untracked()
                    .ok_or("Focus a protein alignment row first.")?;
                let snap = alignment
                    .get_untracked()
                    .ok_or("Open an alignment first.")?;
                let (o, l) = settings.get_untracked()?;
                let c = Capture::prepare(
                    epoch.get_untracked(),
                    snap.doc(),
                    source,
                    o,
                    l,
                    choice.get_untracked(),
                )?;
                Ok::<_, String>((snap, c))
            })();
            let (snap, c) = match prepared {
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
            busy.set(true);
            report.set(None);
            on_report.run(None);
            status.set("Analyzing the exact focused protein row…".into());
            let alive = alive.clone();
            wasm_bindgen_futures::spawn_local(async move {
                let result = w
                    .ask(Request {
                        id: token,
                        revision: c.revision.clone(),
                        job: AnalysisJob::AnalyzeAlignmentAntibody {
                            alignment: snap.doc().clone(),
                            source: c.source.clone(),
                            options: c.options.clone(),
                            limits: c.limits.clone(),
                            selection: c.selection.clone(),
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
                    let (o, l) = settings.get_untracked()?;
                    c.guard(
                        epoch.get_untracked(),
                        &pin.get_untracked().ok_or("The focused row changed.")?,
                        &o,
                        &l,
                        &choice.get_untracked(),
                    )?;
                    match result?.result {
                        AnalysisResult::AlignmentAntibody(r) => {
                            c.validate_report(snap.doc(), &r)?;
                            Ok(r)
                        }
                        AnalysisResult::Error(e) => Err(e),
                        _ => Err("Unexpected alignment antibody response.".into()),
                    }
                })();
                match result {
                    Ok(r) => {
                        let r = Arc::new(r);
                        status.set(format!(
                            "{} retained conditional candidates; choose one to color this row.",
                            r.analysis.tracks.len()
                        ));
                        candidate.set(0);
                        on_report.run(Some((c.epoch, r.clone())));
                        report.set(Some((c, r)));
                    }
                    Err(e) => status.set(e),
                }
            });
        }
    };
    view! {<details class="alignment-rows"><summary>"Antibody analysis of focused row"</summary>
        <p class="helper">"Human V/J references · conditional numbering and independently selected regions · up to 2,000 ungapped residues. Alignment gaps never receive labels. Focus another row to change source; moving between columns keeps the current row analysis."</p>
        <label>"Numbering scheme"<select prop:value=move||format!("{:?}",scheme.get()) on:change=move|e|scheme.set(match event_target_value(&e).as_str(){"Kabat1991"=>t::NumberingScheme::Kabat1991,"Chothia1997"=>t::NumberingScheme::Chothia1997,_=>t::NumberingScheme::Imgt})><option value="Imgt">"IMGT"</option><option value="Kabat1991">"Kabat (1991 profiles)"</option><option value="Chothia1997">"Chothia (1997 profiles)"</option></select></label>
        <label>"Region definition"<select prop:value=move||format!("{:?}",definition.get()) on:change=move|e|definition.set(match event_target_value(&e).as_str(){"Kabat1991"=>t::RegionDefinition::Kabat1991,"ChothiaCanonical"=>t::RegionDefinition::ChothiaCanonical,_=>t::RegionDefinition::Imgt})><option value="Imgt">"IMGT"</option><option value="Kabat1991">"Kabat (1991)"</option><option value="ChothiaCanonical">"Chothia (canonical charts)"</option></select></label>
        <p class="helper">"Number labels and CDR definitions are independent. Unsupported windows remain unresolved; the full report retains diagnostics."</p>
        <details><summary>"Analysis settings"</summary><p class="helper">"Heuristic V weights 5/8/−6, affine gaps −12/−2; minimum CDR3 length 2."</p>
            <label>"Minimum V score"<input type="number" prop:value=move||minimum.get() on:input=move|e|minimum.set(event_target_value(&e))/></label>
            <label>"Allowed FR4 substitutions (0–3)"<input type="number" min="0" max="3" prop:value=move||substitutions.get() on:input=move|e|substitutions.set(event_target_value(&e))/></label>
            <label>"Maximum CDR3 residues"<input type="number" min="2" max="2000" prop:value=move||loop_max.get() on:input=move|e|loop_max.set(event_target_value(&e))/></label>
            <label>"Retained candidates (0–64)"<input type="number" min="0" max="64" prop:value=move||retained.get() on:input=move|e|retained.set(event_target_value(&e))/></label>
        </details>
        <button class="button secondary" disabled=move||busy.get() on:click=run>"Analyze focused row"</button>
        <Show when=move||busy.get()><button class="text-tool" on:click=move|_|{generation.update(|g|*g=g.wrapping_add(1));if let Some(w)=worker.get_value(){w.close();}worker.set_value(None);busy.set(false);status.set("Analysis cancelled.".into());}>"Cancel"</button></Show>
        <p class="helper" role="status">{move||status.get()}</p>
        <Show when=move||report.with(|r|r.is_some())>
            <label>"Conditional candidate"<select prop:value=move||candidate.get().to_string() on:change=move|e|candidate.set(event_target_value(&e).parse().unwrap_or(0))>{move||report.with(|r|r.as_ref().map(|(_,r)|r.analysis.tracks.iter().enumerate().map(|(i,t)|view!{<option value=i.to_string()>{format!("{} · {:?} · {} unassigned residues",i+1,t.locus,t.unassigned_residues)}</option>}).collect_view()))}</select></label>
            <p class="helper">{move||report.with(|r|r.as_ref().map(|(_,r)|format!("{} pairs in retained V scope; V endpoints truncated: {}; pairs truncated: {}. CDR colors are conditional; grey underlines remain unresolved. Focus cells for exact labels under the selected scheme and their origins.",r.analysis.evidence.qualified_pairs_in_retained_v,r.analysis.evidence.v_retention_truncated,r.analysis.evidence.pair_retention_truncated)).unwrap_or_default())}</p>
            {move||report.with(|r|r.as_ref().map(|(_,r)|crate::biology::antibody::antibody_overlay::cdr_summaries(&r.analysis,candidate.get()).into_iter().map(|s|view!{<p class="helper">{s}</p>}).collect_view()))}
            <details><summary>"Mapping diagnostics"</summary>{move||report.with(|r|r.as_ref().map(|(_,r)|crate::biology::antibody::antibody_overlay::diagnostics(&r.analysis,candidate.get()).into_iter().map(|s|view!{<p class="helper">{s}</p>}).collect_view()))}</details>
            <button class="button secondary" on:click=move|_|{let result=(||{let(c,r)=report.get_untracked().ok_or("Analyze the row first.")?;let(o,l)=settings.get_untracked()?;c.guard(epoch.get_untracked(),&pin.get_untracked().ok_or("The row changed.")?,&o,&l, &choice.get_untracked())?;serde_json::to_vec_pretty(r.as_ref()).map_err(|e|e.to_string())})();match result {Ok(bytes)=>crate::io::download::download_with_completion("alignment-antibody.json","application/json",&bytes,move|r|if status.try_get_untracked().is_some(){if let Err(e)=r{status.set(e);}}),Err(e)=>status.set(e)}}>"Download full row analysis"</button>
        </Show>
    </details>}
}
