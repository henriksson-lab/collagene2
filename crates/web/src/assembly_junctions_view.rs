//! Explicit fresh assembly-junction review. Scientific reconstruction/translation
//! stays in a private worker; the host capture callback guards before cloning.
use crate::{app::download, worker::CachedWorker};
use collagene_core::{
    Molecule,
    assembly_junctions::{self, Frame, Options, Report},
    protocol::{AnalysisJob, AnalysisResult, Request},
    verification::DesignPlan,
};
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
struct Form {
    flank: String,
    code: u8,
    frames: [bool; 6],
    junctions: String,
    cds: bool,
}
impl Default for Form {
    fn default() -> Self {
        Self {
            flank: "30".into(),
            code: 1,
            frames: [true; 6],
            junctions: "128".into(),
            cds: true,
        }
    }
}
impl Form {
    fn options(&self) -> Result<Options, String> {
        let options = Options {
            flank_bases: self.flank.parse().map_err(|_| "Invalid DNA flank size")?,
            genetic_code: self.code,
            frames: self
                .frames
                .iter()
                .enumerate()
                .filter(|(_, on)| **on)
                .map(|(i, _)| Frame {
                    reverse: i >= 3,
                    frame: (i % 3) as u8,
                })
                .collect(),
            max_junctions: self
                .junctions
                .parse()
                .map_err(|_| "Invalid junction limit")?,
            include_cds: self.cds,
        };
        options.validate().map_err(|e| e.to_string())?;
        Ok(options)
    }
}
struct Evidence {
    report: Report,
    identity: (String, String),
    epoch: u64,
    form: Form,
}
#[component]
pub fn AssemblyJunctionsPanel(
    capture: Callback<(), Result<(DesignPlan, Molecule), String>>,
    identity: Memo<(String, String)>,
    epoch: RwSignal<u64>,
) -> impl IntoView {
    let form = RwSignal::new(Form::default());
    let status = RwSignal::new(String::new());
    let busy = RwSignal::new(false);
    let generation = RwSignal::new(0u32);
    let evidence = RwSignal::new(None::<Arc<Evidence>>);
    let page = RwSignal::new(0usize);
    let worker = StoredValue::new_local(None::<Rc<CachedWorker>>);
    let alive = Arc::new(AtomicBool::new(true));
    let cleanup = alive.clone();
    on_cleanup(move || {
        cleanup.store(false, Ordering::Relaxed);
        if let Some(w) = worker.get_value() {
            w.close();
        }
        worker.set_value(None);
    });
    Effect::new(move |_| {
        let _ = (identity.get(), epoch.get(), form.get());
        generation.update(|g| *g = g.wrapping_add(1));
        if let Some(w) = worker.get_value() {
            w.close();
        }
        worker.set_value(None);
        evidence.set(None);
        busy.set(false);
        page.set(0);
        status.set(String::new());
    });
    let inspect = Callback::new(move |()| {
        let selected = form.get_untracked();
        let options = match selected.options() {
            Ok(v) => v,
            Err(e) => {
                status.set(e);
                return;
            }
        };
        let bound = identity.get_untracked();
        let bound_epoch = epoch.get_untracked();
        if bound.0.is_empty() || bound.1.is_empty() {
            status.set(
                "Verify a current plan against its selected target before inspecting junctions."
                    .into(),
            );
            return;
        }
        let (plan, target) = match capture.run(()) {
            Ok(v) => v,
            Err(e) => {
                status.set(e);
                return;
            }
        };
        let revision = match assembly_junctions::inspection_revision(&plan, &target, &options) {
            Ok(v) => v,
            Err(e) => {
                status.set(e.to_string());
                return;
            }
        };
        let target_revision = target.revision();
        if identity.get_untracked() != bound || epoch.get_untracked() != bound_epoch {
            status.set("Plan or target changed before inspection.".into());
            return;
        }
        let fresh = match CachedWorker::new() {
            Ok(v) => Rc::new(v),
            Err(e) => {
                status.set(e);
                return;
            }
        };
        if let Some(old) = worker.get_value() {
            old.close();
        }
        worker.set_value(Some(fresh.clone()));
        let token = generation.get_untracked().wrapping_add(1);
        generation.set(token);
        evidence.set(None);
        busy.set(true);
        status.set("Reconstructing original sources and inspecting actual junction DNA…".into());
        let alive = alive.clone();
        spawn_local(async move {
            let response = fresh
                .ask(Request {
                    id: token,
                    revision: revision.clone(),
                    job: AnalysisJob::InspectAssemblyJunctions {
                        plan,
                        target,
                        options,
                    },
                })
                .await;
            fresh.close();
            if !alive.load(Ordering::Relaxed)
                || generation.get_untracked() != token
                || identity.get_untracked() != bound
                || epoch.get_untracked() != bound_epoch
                || form.get_untracked() != selected
            {
                return;
            }
            worker.set_value(None);
            busy.set(false);
            match response {
                Ok(response) => match response.result {
                    AnalysisResult::AssemblyJunctionInspection(report)
                        if report.context_revision == revision
                            && report.plan_revision == bound.0
                            && report.target_revision == target_revision =>
                    {
                        page.set(0);
                        status.set("Fresh junction review ready. Frame translations describe sequence only.".into());
                        evidence.set(Some(Arc::new(Evidence {
                            report,
                            identity: bound,
                            epoch: bound_epoch,
                            form: selected,
                        })));
                    }
                    AnalysisResult::Error(error) => status.set(error),
                    _ => status.set(
                        "Junction response did not match the requested plan, target and options."
                            .into(),
                    ),
                },
                Err(error) => status.set(error),
            }
        });
    });
    let cancel = move |_| {
        generation.update(|g| *g = g.wrapping_add(1));
        if let Some(w) = worker.get_value() {
            w.close();
        }
        worker.set_value(None);
        busy.set(false);
        evidence.set(None);
        status.set("Junction inspection cancelled.".into());
    };
    view! {<details class="assembly-junction-review"><summary>"Review junction DNA and reading frames"</summary>
        <p class="helper">"Explicitly reconstruct the original plan and inspect junction flanks, homology and reading frames. This does not establish expression, protein function, a unique mixture outcome or experimental success."</p>
        <label>"DNA flank on each side (bases)"<input type="number" min="1" max="150" aria-label="Junction DNA flank bases" prop:value=move||form.with(|f|f.flank.clone()) on:input=move|e|form.update(|f|f.flank=event_target_value(&e))/></label>
        <label>"Genetic code"<select aria-label="Junction genetic code" prop:value=move||form.with(|f|f.code.to_string()) on:change=move|e|form.update(|f|f.code=event_target_value(&e).parse().unwrap_or(1))><option value="1">"1 — Standard"</option><option value="2">"2 — Vertebrate mitochondrial"</option><option value="11">"11 — Bacterial, archaeal and plant plastid"</option></select></label>
        <fieldset><legend>"Window-relative reading frames"</legend>{(0..6).map(|i|{let label=format!("{} frame {}",if i>=3{"Reverse-complement"}else{"Forward"},i%3+1);view!{<label class="checkbox-label"><input type="checkbox" aria-label=label.clone() prop:checked=move||form.with(|f|f.frames[i]) on:change=move|e|form.update(|f|f.frames[i]=event_target_checked(&e))/>{label}</label>}}).collect_view()}</fieldset>
        <p class="helper">"Frames are relative to the returned local DNA window (or its reverse complement), not inferred gene frames. * marks a stop; X marks unresolved translation."</p>
        <label>"Maximum junctions (1–128)"<input type="number" min="1" max="128" aria-label="Maximum assembly junctions" prop:value=move||form.with(|f|f.junctions.clone()) on:input=move|e|form.update(|f|f.junctions=event_target_value(&e))/></label>
        <label class="checkbox-label"><input type="checkbox" aria-label="Review CDS annotations at junctions" prop:checked=move||form.with(|f|f.cds) on:change=move|e|form.update(|f|f.cds=event_target_checked(&e))/>"Review available CDS annotations"</label>
        <button class="button secondary full" disabled=move||busy.get()||identity.with(|i|i.0.is_empty()||i.1.is_empty()) on:click=move|_|inspect.run(())>"Inspect assembly junctions"</button><button class="text-tool" disabled=move||!busy.get() on:click=cancel>"Cancel junction inspection"</button><p class="assembly-junction-status" role="status">{move||status.get()}</p>
        {move||evidence.get().map(|shown|{let saved=shown.clone();let count=shown.report.junctions.len();let pages=count.max(1).div_ceil(8);let current=page.get().min(pages-1);view!{<div class="assembly-junction-report"><p>{format!("Inspection complete within model: {} · {} junctions retained · {} omitted by limits",shown.report.complete,count,shown.report.omitted_junctions)}</p><p>{format!("Product reference: {}",shown.report.product_revision.as_deref().unwrap_or("unavailable"))}</p><p>{format!("Independent-target assessment: {}",match shown.report.verification.status{collagene_core::verification::Status::Invalid=>"inconsistent",collagene_core::verification::Status::Incomplete=>"incomplete",collagene_core::verification::Status::ConsistentUnderAssumptions=>"consistent under stated assumptions"})}</p><details><summary>"Source and calculation details"</summary><p class="helper">{format!("Plan {} · target {} · engine {}",shown.report.plan_revision,shown.report.target_revision,shown.report.engine_version)}</p></details><details><summary>"Findings and limitations"</summary>{shown.report.warnings.iter().chain(&shown.report.verification.findings).map(|w|view!{<p>{w.clone()}</p>}).collect_view()}</details>{shown.report.junctions.iter().skip(current*8).take(8).map(|junction|{
            let dna=junction.window.dna.clone();let boundary=junction.window.boundary_offset;let marked=dna.split_at_checked(boundary).map(|(left,right)|format!("{left}|{right}")).unwrap_or_else(||dna.clone());
            view!{<article class="assembly-junction"><h4>{format!("Junction {} · fragment {} → {}{}",junction.index+1,junction.upstream+1,junction.downstream+1,if junction.closing{" · circle closure"}else{""})}</h4><p>{format!("Product boundary {} (zero-based, between bases)",junction.boundary)}</p><p>{format!("DNA window product ranges: {}",junction.window.parts.iter().map(|r|format!("[{}, {})",r.start,r.end)).collect::<Vec<_>>().join(", "))}</p><p class="helper">{format!("Zero-based half-open product coordinates; {} left / {} right bases; window clipped: {}",junction.window.left_bases,junction.window.right_bases,junction.window.clipped)}</p><pre>{format!("Window 5′→3′: {marked}\n| marks the junction at window offset {boundary}")}</pre><p>{format!("Homology product ranges: {}",if junction.homology.is_empty(){"none".into()}else{junction.homology.iter().map(|r|format!("[{}, {})",r.start,r.end)).collect::<Vec<_>>().join(", ")})}</p>
                {junction.frames.iter().map(|frame|view!{<details><summary>{format!("{} frame {} · {} amino acids",if frame.spec.reverse{"Reverse-complement window"}else{"Forward window"},frame.spec.frame+1,frame.translation.protein.len())}</summary><pre>{frame.translation.protein.clone()}</pre><p>{format!("Junction in oriented window: offset {} · {} trailing bases · {} unresolved codons · {} possible stop codons",frame.boundary_offset,frame.translation.trailing_bases,frame.translation.unresolved_codons,frame.translation.possible_stop_codons)}</p>{frame.crossing_codon.as_ref().map(|codon|view!{<p>{format!("Boundary-crossing codon {} → {} at oriented-window offset {}",codon.dna,codon.amino_acid,codon.window_start)}</p>})}{frame.translation.warnings.iter().map(|w|view!{<p>{w.clone()}</p>}).collect_view()}</details>}).collect_view()}
                {junction.cds.iter().map(|cds|view!{<details><summary>{format!("CDS annotation {}: {}{}",cds.feature_index+1,cds.label,if cds.needs_review{" · review required"}else{""})}</summary><p>{format!("Annotation spans boundary: {}",cds.spans_boundary)}</p>{cds.warning.as_ref().map(|w|view!{<p>{w.clone()}</p>})}{cds.translation.as_ref().map(|translation|view!{<pre>{format!("Coding DNA: {}\nPeptide: {}",translation.dna,translation.protein)}</pre><p>{format!("CDS codon_start {} · coding offset {} · first codon {} · {} codons · {} · code {} ({}) · complete under translation model: {}",translation.codon_start,translation.junction_coding_offset,translation.first_codon,translation.codon_count,if translation.reverse{"reverse strand"}else{"forward strand"},translation.genetic_code,translation.code_name,translation.complete_under_model)}</p>})}<p class="helper">"A local codon window and annotation correspondence are not evidence of a functional fusion or expression."</p></details>}).collect_view()}<p class="helper">{format!("Additional CDS annotations omitted: {}",junction.cds_omitted)}</p>
            </article>}
        }).collect_view()}<div class="method-row-actions"><button class="text-tool" disabled=current==0 on:click=move|_|page.set(current.saturating_sub(1))>"Previous junctions"</button><span>{format!("Page {} of {}",current+1,pages)}</span><button class="text-tool" disabled={current+1>=pages} on:click=move|_|page.set((current+1).min(pages-1))>"Next junctions"</button></div><button class="button secondary full" on:click=move|_|{
            if identity.get_untracked()!=saved.identity||epoch.get_untracked()!=saved.epoch||form.get_untracked()!=saved.form{status.set("Context changed; inspect junctions again before exporting current evidence.".into());return;}
            match serde_json::to_vec_pretty(&saved.report){Ok(bytes)=>if let Err(error)=download("assembly-junction-review.json","application/json",&bytes){status.set(error)},Err(error)=>status.set(error.to_string())}
        }>"Download junction evidence JSON"</button></div>}})}
    </details>}
}
