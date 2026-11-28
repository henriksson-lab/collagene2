//! One explicit oligo proposal for the selected overlap construct.
use crate::assembly::draft::model::Draft;
use crate::assembly::primer_design::assembly_edited_primer_model as model;
use crate::assembly::review::assembly_review::IntendedReference;
use crate::platform::worker::CachedWorker;
use collagene_core::{
    cloning::assembly::recipe::{Preparation, RunPreparation},
    documents::workspace::Workspace,
    protocol::{
        request::Request,
        schema::{AnalysisJob, AnalysisResult},
    },
    sequence::analysis::Primer,
    verification::evidence::Status,
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
struct Settings {
    fields: [String; 10],
    epoch: u64,
    current: bool,
    draft: Draft,
}
fn status_label(status: &Status) -> &'static str {
    match status {
        Status::ConsistentUnderAssumptions => "consistent under stated assumptions",
        Status::Incomplete => "incomplete",
        Status::Invalid => "inconsistent",
    }
}
fn integer(text: &str, name: &str) -> Result<usize, String> {
    if text.is_empty() || text.len() > 20 || !text.bytes().all(|b| b.is_ascii_digit()) {
        return Err(format!("Enter a whole number for {name}."));
    }
    text.parse().map_err(|_| format!("{name} is too large."))
}
fn temperature(text: &str) -> Result<f64, String> {
    if text.len() > 32 {
        return Err("The target temperature field is too long.".into());
    }
    let value = text
        .parse::<f64>()
        .map_err(|_| "Enter a target annealing Tm.")?;
    if !value.is_finite() {
        return Err("The target annealing Tm must be finite.".into());
    }
    Ok(value)
}
fn oe_oligos(
    prep: &collagene_core::cloning::overlap::extension::pcr::Request,
) -> Vec<(&'static str, &Primer)> {
    vec![
        ("Left PCR forward", &prep.left.forward),
        ("Left PCR reverse", &prep.left.reverse),
        ("Right PCR forward", &prep.right.forward),
        ("Right PCR reverse", &prep.right.reverse),
        ("Outer PCR forward", &prep.outer.forward),
        ("Outer PCR reverse", &prep.outer.reverse),
    ]
}
fn oligos(prep: &Preparation) -> Vec<(&'static str, &Primer)> {
    match prep {
        Preparation::ExactPcr {
            template_strand: _,
            forward,
            reverse,
        } => {
            vec![("Forward", forward), ("Reverse", reverse)]
        }
        Preparation::OverlapExtensionPcr { request } => oe_oligos(request),
        Preparation::FromRun { preparation, .. } | Preparation::FromDerived { preparation, .. } => {
            match preparation {
                RunPreparation::ExactPcr {
                    template_strand: _,
                    forward,
                    reverse,
                } => {
                    vec![("Forward", forward), ("Reverse", reverse)]
                }
                RunPreparation::OverlapExtensionPcr { request } => oe_oligos(request),
                _ => vec![],
            }
        }
        _ => vec![],
    }
}
#[component]
pub fn EditedAssemblyPrimerView(
    draft: RwSignal<Draft>,
    intent: Arc<IntendedReference>,
    epoch: RwSignal<u64>,
    snapshot: Callback<(), Result<Workspace, String>>,
    current: Memo<bool>,
    guard: Callback<(), Result<(), String>>,
    on_use: Callback<(Draft, Draft), Result<(), String>>,
) -> impl IntoView {
    let defaults = model::Options::default();
    let overlap = RwSignal::new(defaults.overlap.to_string());
    let minimum = RwSignal::new(defaults.preparation.primer_options.min_length.to_string());
    let maximum = RwSignal::new(defaults.preparation.primer_options.max_length.to_string());
    let tm = RwSignal::new(defaults.preparation.primer_options.target_tm_c.to_string());
    let flank = RwSignal::new(defaults.preparation.homology_flank.to_string());
    let max_oligo = RwSignal::new(defaults.preparation.max_oligo_length.to_string());
    let min_paired = RwSignal::new(defaults.preparation.min_paired_bases.to_string());
    let comparisons = RwSignal::new(defaults.preparation.max_comparisons.to_string());
    let candidates = RwSignal::new(defaults.preparation.max_candidates.to_string());
    let pairs = RwSignal::new(defaults.preparation.primer_options.max_pairs.to_string());
    let intent = StoredValue::new(intent);
    let pending = RwSignal::new(false);
    let status = RwSignal::new(String::new());
    let proposal = RwSignal::new(None::<Arc<model::Proposal>>);
    let generation = RwSignal::new(0u32);
    let worker = StoredValue::new_local(None::<Rc<CachedWorker>>);
    let alive = Arc::new(AtomicBool::new(true));
    let cleanup = alive.clone();
    on_cleanup(move || {
        cleanup.store(false, Ordering::Relaxed);
        if let Some(worker) = worker.get_value() {
            worker.close();
        }
    });
    let settings = Memo::new(move |_| Settings {
        fields: [
            overlap.get(),
            minimum.get(),
            maximum.get(),
            tm.get(),
            flank.get(),
            max_oligo.get(),
            min_paired.get(),
            comparisons.get(),
            candidates.get(),
            pairs.get(),
        ],
        epoch: epoch.get(),
        current: current.get(),
        draft: model::context_draft(draft.get()),
    });
    Effect::new(move |_| {
        let _ = settings.get();
        generation.update(|g| *g = g.wrapping_add(1));
        if let Some(worker) = worker.get_value() {
            worker.close();
        }
        worker.set_value(None);
        pending.set(false);
        proposal.set(None);
        status.set(String::new());
    });
    let eligible = Memo::new(move |_| draft.with(model::eligible));
    let read_options = move || {
        let mut options = model::Options::default();
        options.overlap = integer(&overlap.get_untracked(), "uniform junction overlap")?;
        let p = &mut options.preparation;
        p.primer_options.min_length =
            integer(&minimum.get_untracked(), "minimum annealing length")?;
        p.primer_options.max_length =
            integer(&maximum.get_untracked(), "maximum annealing length")?;
        p.primer_options.target_tm_c = temperature(&tm.get_untracked())?;
        p.primer_options.max_pairs = integer(&pairs.get_untracked(), "ranked primer pairs")?;
        p.homology_flank = integer(&flank.get_untracked(), "unchanged OE flank")?;
        p.max_oligo_length = integer(&max_oligo.get_untracked(), "maximum total oligo length")?;
        p.min_paired_bases = integer(&min_paired.get_untracked(), "minimum terminal pairing")?;
        p.max_comparisons = integer(&comparisons.get_untracked(), "terminal comparison budget")?;
        p.max_candidates = integer(&candidates.get_untracked(), "retained terminal pairings")?;
        Ok::<_, String>(options)
    };
    let prepare = Callback::new(move |()| {
        let captured = (|| {
            guard.run(())?;
            if !current.get_untracked() {
                return Err("Review the unchanged intended construct before proposing all assembly primers.".into());
            }
            let baseline = draft.get_untracked();
            let workspace = snapshot.run(())?;
            let events = model::attest(&baseline, &workspace)?;
            let options = read_options()?;
            let (capture, inputs, request) = intent.with_value(|i| {
                model::request(&baseline, i.as_ref(), &events, &options, Some(&workspace))
            })?;
            Ok::<_, String>((
                baseline, events, options, capture, inputs, request, workspace,
            ))
        })();
        let (baseline, events, options, capture, inputs, request, workspace) = match captured {
            Ok(v) => v,
            Err(error) => {
                status.set(error);
                return;
            }
        };
        let live = match CachedWorker::new() {
            Ok(w) => Rc::new(w),
            Err(error) => {
                status.set(error);
                return;
            }
        };
        if let Some(old) = worker.get_value() {
            old.close();
        }
        worker.set_value(Some(live.clone()));
        generation.update(|g| *g = g.wrapping_add(1));
        let token = generation.get_untracked();
        let saved_settings = settings.get_untracked();
        pending.set(true);
        proposal.set(None);
        status.set("Proposing all original-source oligos, including edited insert preparations, and checking the complete assembly…".into());
        let alive = alive.clone();
        spawn_local(async move {
            let revision = capture.design_revision.clone();
            let response = live
                .ask(Request {
                    id: token,
                    revision: revision.clone(),
                    job: AnalysisJob::DesignEditedOverlapPrimers {
                        inputs,
                        request,
                        workspace: Some(workspace),
                    },
                })
                .await;
            live.close();
            if !alive.load(Ordering::Relaxed)
                || generation.get_untracked() != token
                || settings.get_untracked() != saved_settings
                || !current.get_untracked()
            {
                return;
            }
            worker.set_value(None);
            pending.set(false);
            let checked = (|| {
                guard.run(())?;
                let now = draft.get_untracked();
                if model::context_draft(now.clone()) != model::context_draft(baseline.clone()) {
                    return Err("The assembly draft changed while proposing primers.".into());
                }
                let fresh_workspace = snapshot.run(())?;
                let sources = model::attest(&now, &fresh_workspace)?;
                if sources != events {
                    return Err(
                        "An exact original-source event changed while proposing primers.".into(),
                    );
                }
                let response = response?;
                if response.id != token || response.revision != revision {
                    return Err("The all-primer reply belongs to another request.".into());
                }
                let report = match response.result {
                    AnalysisResult::EditedOverlapPrimersDesigned(report) => report,
                    AnalysisResult::Error(error) => return Err(error),
                    _ => return Err("Unexpected whole-assembly primer response.".into()),
                };
                let proposed = intent.with_value(|i| {
                    model::from_design(
                        &now,
                        i.as_ref(),
                        &sources,
                        &options,
                        report,
                        Some(&fresh_workspace),
                    )
                })?;
                if proposed.capture != capture {
                    return Err("The original whole-assembly declaration changed.".into());
                }
                Ok::<_, String>(proposed)
            })();
            match checked {
                Ok(p) => {
                    status.set("Review every ordered oligo and the complete assembly check. Applying changes preparations and the requested uniform overlap; sources and intended DNA stay fixed.".into());
                    proposal.set(Some(Arc::new(p)));
                }
                Err(error) => status.set(error),
            }
        });
    });
    let apply = Callback::new(move |saved: Arc<model::Proposal>| {
        if !proposal.with_untracked(|p| p.as_ref().is_some_and(|p| Arc::ptr_eq(p, &saved))) {
            return Err("This whole-assembly proposal is no longer current.".into());
        }
        guard.run(())?;
        if !current.get_untracked() {
            return Err("The independent intended construct is stale.".into());
        }
        let before = draft.get_untracked();
        let workspace = snapshot.run(())?;
        let events = model::attest(&before, &workspace)?;
        let options = read_options()?;
        let next = intent.with_value(|i| {
            saved.apply(&before, i.as_ref(), &events, &options, Some(&workspace))
        })?;
        on_use.run((before, next))?;
        proposal.set(None);
        Ok::<_, String>(())
    });
    view! {<Show when=move||draft.with(|d|d.methods.mode==crate::cloning::preparation::cloning_methods_model::Mode::Overlap)><details class="assembly-fragment-options"><summary>"Propose primers for this construct"</summary>
        <p class="helper">"Design every contribution together. Edited contributions use two source PCRs, reciprocal extension and outer PCR; other contributions use exact PCR. Junction tails use your literal inserts and the next edited contribution’s homology. Original sources and intended DNA remain unchanged."</p>
        <Show when=move||!eligible.get()><p class="helper">"Use two to twelve overlap contributions. Keep active digest, paired-oligo and declared-chemistry preparations in their existing workflow, or explicitly choose plain PCR before using this proposal."</p></Show>
        <div class="primer-grid">
            <label>"Uniform junction homology (bp)"<input type="number" min="1" max="80" prop:value=move||overlap.get() on:input=move|e|overlap.set(event_target_value(&e))/></label>
            <label>"Minimum annealing length"<input type="number" min="8" max="80" prop:value=move||minimum.get() on:input=move|e|minimum.set(event_target_value(&e))/></label>
            <label>"Maximum annealing length"<input type="number" min="8" max="80" prop:value=move||maximum.get() on:input=move|e|maximum.set(event_target_value(&e))/></label>
            <label>"Target annealing Tm (°C)"<input type="number" min="-50" max="100" step="0.1" prop:value=move||tm.get() on:input=move|e|tm.set(event_target_value(&e))/></label>
            <label>"Unchanged OE flank on each side (bp)"<input type="number" min="1" max="200" prop:value=move||flank.get() on:input=move|e|flank.set(event_target_value(&e))/></label>
            <label>"Maximum total oligo length"<input type="number" min="8" max="200" prop:value=move||max_oligo.get() on:input=move|e|max_oligo.set(event_target_value(&e))/></label>
        </div>
        <details><summary>"Search limits and existing tails"</summary>
            <p class="helper">"Nonempty active outer-PCR tails must already match the assembly-derived tails. Conflicts require an explicit change in that contribution’s editor; this proposal does not silently replace them."</p>
            <div class="primer-grid">
                <label>"Ranked pairs per primer call"<input type="number" min="1" max="20" prop:value=move||pairs.get() on:input=move|e|pairs.set(event_target_value(&e))/></label>
                <label>"Minimum OE terminal pairing (bp)"<input type="number" min="1" prop:value=move||min_paired.get() on:input=move|e|min_paired.set(event_target_value(&e))/></label>
                <label>"Terminal comparison budget"<input type="number" min="1" max="50000000" prop:value=move||comparisons.get() on:input=move|e|comparisons.set(event_target_value(&e))/></label>
                <label>"Retained terminal pairings"<input type="number" min="1" max="1000" prop:value=move||candidates.get() on:input=move|e|candidates.set(event_target_value(&e))/></label>
            </div>
        </details>
        <button class="button secondary" disabled=move||pending.get()||!current.get()||!eligible.get() on:click=move|_|prepare.run(())>"Propose primers for all contributions"</button>
        <Show when=move||pending.get()><button class="text-tool" on:click=move|_|{generation.update(|g|*g=g.wrapping_add(1));if let Some(w)=worker.get_value(){w.close();}worker.set_value(None);pending.set(false);proposal.set(None);status.set("Whole-assembly primer proposal canceled.".into());}>"Cancel"</button></Show>
        <p role="status">{move||status.get()}</p>
        {move||proposal.get().map(|saved|{
            let invalid=saved.verification.status==Status::Invalid;
            let diagnostic=format!("Complete assembly check: {}. This check uses the canonical declared reference. After applying, rerun the selected assembly against your unchanged intended construct.",status_label(&saved.verification.status));
            let saved=StoredValue::new(saved);
            view!{<section class="assembly-edited-primer-proposal"><p>{diagnostic}</p>
                {saved.with_value(|p|p.recipe.fragments.iter().enumerate().map(|(index,fragment)|{
                    let primers=oligos(&fragment.preparation);let count=primers.len();
                    view!{<details><summary>{format!("Fragment {} · {} · {} ordered oligos",index+1,fragment.role,count)}</summary>
                        <table><thead><tr><th>"Stage"</th><th>"Primer"</th><th>"Ordered sequence · 5′→3′"</th></tr></thead><tbody>{primers.into_iter().map(|(stage,primer)|view!{<tr><td>{stage}</td><td>{primer.name.clone()}</td><td><pre>{primer.ordered_sequence()}</pre></td></tr>}).collect_view()}</tbody></table>
                        {p.rows().get(index).map(|row|view!{<p>{format!("PCR source frame: {:?}",row.template_strand)}</p>})}{p.tails().get(index).map(|tail|view!{<p class="helper">{format!("Declared junction insert {} bp · outgoing homology {} bp{}",tail.outgoing_insert.len(),tail.homology.len(),tail.downstream_role.as_ref().map(|role|format!(" toward {role}")).unwrap_or_default())}</p>})}
                        {p.rows().get(index).map(|row|view!{<details><summary>"Primer search notes"</summary>{row.search.iter().map(|search|view!{<section><p>{format!("{} · examined {} · candidates {} · rejected {} · capped {} · unique geometry {}",search.call,search.examined,search.candidates,search.rejected,search.capped,search.unique_geometry)}</p><p>{format!("Showing {} of {} notes",search.notes.len().min(20),search.notes.len())}</p>{search.notes.iter().take(20).map(|note|view!{<p class="helper">{note.clone()}</p>}).collect_view()}</section>}).collect_view()}</details>})}
                        {p.source_frame(&fragment.input_id).map(|frame|view!{<details><summary>{format!("Actual {:?} template · {} calls",frame.template_strand,frame.sequence.len())}</summary><pre>{frame.sequence.as_str().chars().take(200).collect::<String>()}</pre><p class="helper">{format!("{} of {} native template calls shown. Original top-reference positions exist only where mapped by the source evidence.",frame.sequence.len().min(200),frame.sequence.len())}</p><code>{frame.state_revision.clone()}</code></details>})}
                    </details>}
                }).collect_view())}
                {saved.with_value(|p|p.warnings.iter().take(24).map(|warning|view!{<p class="helper">{warning.clone()}</p>}).collect_view())}
                <details><summary>"Whole-recipe findings"</summary>{saved.with_value(|p| crate::protocol::verification_view::findings(&p.verification.findings))}</details>
                <button class="button primary" disabled=move||invalid||!current.get() on:click=move|_|match apply.run(saved.get_value()){Ok(())=>status.set("All proposed oligos entered. Reconstruct and verify the complete assembly against the unchanged intended construct.".into()),Err(error)=>status.set(error)}>"Use all proposed oligos"</button>
            </section>}
        })}
    </details></Show>}
}
