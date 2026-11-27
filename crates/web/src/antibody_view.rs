//! Conditional antibody reference analysis, kept off the UI thread.
use crate::{
    antibody_model::{self, Capture},
    source_identity::SourceInfo,
    worker::CachedWorker,
};
use collagene_core::{
    Molecule,
    antibody_numbering::{domain_candidates as d, domain_tracks as n},
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
fn loop_text(name: &str, outcome: &d::LoopAssignment) -> String {
    match outcome {
        d::LoopAssignment::Numbered(n) => format!(
            "{name}: {}{}",
            n.residues.iter().map(|r| r.residue).collect::<String>(),
            if n.residues.is_empty() { "(empty)" } else { "" }
        ),
        d::LoopAssignment::Unresolved {
            issue: d::LoopIssue::MissingInvariantBoundaries { positions },
            ..
        } => format!(
            "{name}: unresolved boundary {}",
            positions
                .iter()
                .map(u16::to_string)
                .collect::<Vec<_>>()
                .join(", ")
        ),
        d::LoopAssignment::Unresolved {
            issue: d::LoopIssue::UnsupportedLength { residues },
            ..
        } => format!("{name}: unsupported length {residues}"),
    }
}
fn region(r: &n::Region) -> &'static str {
    match r {
        n::Region::Fr1 => "FR1",
        n::Region::Cdr1 => "CDR1",
        n::Region::Fr2 => "FR2",
        n::Region::Cdr2 => "CDR2",
        n::Region::Fr3 => "FR3",
        n::Region::Cdr3 => "CDR3",
        n::Region::Fr4 => "FR4",
    }
}
fn origin(e: &n::AssignmentEvidence) -> &'static str {
    match e {
        n::AssignmentEvidence::InvariantVReference => {
            "Reference correspondence shared by all optimal V paths"
        }
        n::AssignmentEvidence::ConditionalLoopRule => "Conditional published loop-numbering rule",
        n::AssignmentEvidence::JFr4Reference => "Matched complete J/FR4 reference",
    }
}
#[component]
pub fn AntibodyPanel(
    source: Memo<SourceInfo>,
    epoch: RwSignal<u64>,
    workspace_busy: RwSignal<bool>,
    prepare: Callback<EventBinding, Result<Molecule, String>>,
    on_select: Callback<(Capture, usize, usize), Result<(), String>>,
) -> impl IntoView {
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
                .map_err(|_| "Enter a whole maximum CDR3 length.")?,
            retained
                .get()
                .parse()
                .map_err(|_| "Enter a whole candidate limit.")?,
        ))
    });
    let status = RwSignal::new(String::new());
    let busy = RwSignal::new(false);
    let generation = RwSignal::new(0u32);
    let report = RwSignal::new(None::<(Capture, Arc<n::Report>)>);
    let candidate = RwSignal::new(0usize);
    let page = RwSignal::new(0usize);
    let worker = StoredValue::new_local(None::<Rc<CachedWorker>>);
    let alive = Arc::new(AtomicBool::new(true));
    Effect::new(move |_| {
        let _ = (
            source.get(),
            epoch.get(),
            settings.get(),
            workspace_busy.get(),
        );
        generation.update(|g| *g = g.wrapping_add(1));
        if let Some(w) = worker.get_value() {
            w.close();
        }
        worker.set_value(None);
        busy.set(false);
        report.set(None);
        candidate.set(0);
        page.set(0);
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
                let (o, l) = settings.get_untracked()?;
                let s = source.get_untracked();
                let m = prepare.run(s.binding.clone())?;
                let c = Capture::prepare(epoch.get_untracked(), s.binding, &m, o, l)?;
                Ok::<_, String>((m, c))
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
            status.set("Analyzing human V/J reference candidates…".into());
            let alive = alive.clone();
            wasm_bindgen_futures::spawn_local(async move {
                let result = w
                    .ask(Request {
                        id: token,
                        revision: c.context_revision.clone(),
                        job: AnalysisJob::AnalyzeAntibody {
                            molecule: m,
                            options: c.options.clone(),
                            limits: c.limits.clone(),
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
                        &source.get_untracked().binding,
                        &o,
                        &l,
                    )?;
                    match result?.result {
                        AnalysisResult::Antibody(r) => {
                            c.validate_report(&r)?;
                            Ok(r)
                        }
                        AnalysisResult::Error(e) => Err(e),
                        _ => Err("Unexpected antibody analysis response.".into()),
                    }
                })();
                match result {
                    Ok(r) => {
                        status.set(format!("{} retained candidates. Conditional reference evidence; review unresolved assignments.",r.tracks.len()));
                        report.set(Some((c, Arc::new(r))));
                    }
                    Err(e) => status.set(e),
                }
            });
        }
    };
    view! {<section class="inspector-content" aria-label="Antibody reference analysis" style="overflow-wrap:anywhere"><details><summary>"Antibody numbering and CDRs"</summary>
        <p class="helper">"Human V/J references · up to 2,000 residues · conditional IMGT numbering and IMGT regions. Unresolved residues remain unassigned. This is not a calibrated classifier or validated Kabat/Chothia conversion."</p>
        <details><summary>"Analysis settings"</summary><p class="helper">"V weights: observed 5, frequency 8, unseen −6; affine gaps −12/−2. Minimum CDR3 length 2. Score thresholds are heuristic."</p>
            <label>"Minimum V score"<input type="number" prop:value=move||minimum.get() on:input=move|e|minimum.set(event_target_value(&e))/></label>
            <label>"Allowed FR4 substitutions (0–3)"<input type="number" min="0" max="3" prop:value=move||substitutions.get() on:input=move|e|substitutions.set(event_target_value(&e))/></label>
            <label>"Maximum CDR3 residues"<input type="number" min="2" max="2000" prop:value=move||loop_max.get() on:input=move|e|loop_max.set(event_target_value(&e))/></label>
            <label>"Retained candidates (0–64)"<input type="number" min="0" max="64" prop:value=move||retained.get() on:input=move|e|retained.set(event_target_value(&e))/></label>
            <p class="helper">"Zero retained candidates requests counts only. At most 128 optimal V endpoints per locus are retained; omitted endpoints and suboptimal paths are not searched for products."</p>
        </details>
        <button class="button secondary full" disabled=move||busy.get()||workspace_busy.get() on:click=run>"Analyze antibody candidates"</button>
        <Show when=move||busy.get()><button class="text-tool" on:click=move|_|{generation.update(|g|*g=g.wrapping_add(1));if let Some(w)=worker.get_value(){w.close();}worker.set_value(None);busy.set(false);status.set("Analysis cancelled.".into());}>"Cancel"</button></Show>
        <p class="helper" role="status">{move||status.get()}</p>
        <Show when=move||report.with(|r|r.is_some())>
            <p class="helper">{move||report.with(|r|r.as_ref().map(|(_,r)|format!("{} qualifying pairs within retained V ends. V endpoints truncated: {} · pairs truncated: {}.",r.evidence.qualified_pairs_in_retained_v,r.evidence.v_retention_truncated,r.evidence.pair_retention_truncated)).unwrap_or_default())}</p>
            <label>"Candidate"<select prop:value=move||candidate.get().to_string() on:change=move|e|{candidate.set(event_target_value(&e).parse().unwrap_or(0));page.set(0);}>
                {move||report.with(|r|r.as_ref().map(|(_,r)|r.tracks.iter().enumerate().map(|(i,t)|view!{<option value=i.to_string()>{format!("{} · {:?} · supported {}–{}",i+1,t.locus,t.first_supported_source_index+1,t.source_end)}</option>}).collect_view()))}
            </select></label>
            {move||report.with(|r|r.as_ref().and_then(|(_,r)|r.evidence.candidates.get(candidate.get())).map(|c|view!{<p class="helper">{loop_text("CDR1",&c.cdr1)}</p><p class="helper">{loop_text("CDR2",&c.cdr2)}</p><p class="helper">{format!("CDR3: {}",c.cdr3.residues.iter().map(|r|r.residue).collect::<String>())}</p>}))}
            <p class="helper">{move||report.with(|r|r.as_ref().and_then(|(_,r)|r.tracks.get(candidate.get())).map(|t|format!("{} unassigned residues · position 1: {}. Coordinates below are one-based source residues.",t.unassigned_residues,t.invariant_position_one.map(|i|(i+1).to_string()).unwrap_or_else(||"unresolved N terminus".into()))).unwrap_or_default())}</p>
            <table class="data-table"><thead><tr><th>"Source"</th><th>"AA"</th><th>"IMGT"</th><th>"Region"</th></tr></thead><tbody>
                {move||report.with(|r|r.as_ref().and_then(|(c,r)|r.tracks.get(candidate.get()).map(|t|t.residues.iter().skip(page.get()*60).take(60).map(|row|{let index=row.source_index;let capture=c.clone();let assignment=row.assignment.as_ref();view!{<tr><td><button class="text-tool" disabled=move||workspace_busy.get() on:click=move|_|{let result=(||{let(o,l)=settings.get_untracked()?;capture.guard(epoch.get_untracked(),&source.get_untracked().binding,&o,&l)?;on_select.run((capture.clone(),index,index+1))})();if let Err(e)=result{status.set(e);}}>{index+1}</button></td><td>{row.residue.to_string()}</td><td title=assignment.map(|a|origin(&a.evidence)).unwrap_or("No supported assignment")>{assignment.map(|a|a.position.to_string()).unwrap_or_else(||"—".into())}</td><td>{assignment.map(|a|region(&a.region)).unwrap_or("unresolved")}</td></tr>}}).collect_view())))}
            </tbody></table>
            <button class="text-tool" disabled=move||page.get()==0 on:click=move|_|page.update(|p|*p=p.saturating_sub(1))>"Previous residues"</button>
            <button class="text-tool" disabled=move||report.with(|r|r.as_ref().and_then(|(_,r)|r.tracks.get(candidate.get())).is_none_or(|t|(page.get()+1)*60>=t.residues.len())) on:click=move|_|page.update(|p|*p+=1)>"Next residues"</button>
            <button class="button secondary full" disabled=move||workspace_busy.get() on:click=move|_|{let result=(||{let(c,r)=report.get_untracked().ok_or("Analyze the source first.")?;let(o,l)=settings.get_untracked()?;c.guard(epoch.get_untracked(),&source.get_untracked().binding,&o,&l)?;serde_json::to_vec_pretty(r.as_ref()).map_err(|e|e.to_string())})();match result {Ok(bytes)=>crate::app::download_with_completion("antibody-analysis.json","application/json",&bytes,move|r|if status.try_get_untracked().is_some(){if let Err(e)=r{status.set(e);}}),Err(e)=>status.set(e)}}>"Download full analysis"</button>
        </Show>
    </details></section>}
}
