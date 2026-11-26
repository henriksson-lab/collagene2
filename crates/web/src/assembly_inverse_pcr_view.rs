//! Dock preparation for a literal internal edit, using the existing native worker.
use crate::{
    assembly_draft_model::Draft, assembly_inverse_pcr_model as model,
    assembly_review::IntendedReference, worker::CachedWorker,
};
use collagene_core::{
    history::DocumentHistory,
    protocol::{AnalysisJob, AnalysisResult, Request},
    run_lineage::EventBinding,
    transfer::Candidate,
    verification::Status,
    workspace::Workspace,
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

fn options(forward: &str, reverse: &str, split: &str) -> Result<model::Options, String> {
    fn number(text: &str) -> Result<usize, String> {
        if text.len() > 20 || text.is_empty() || !text.bytes().all(|b| b.is_ascii_digit()) {
            return Err("Enter whole-number primer lengths and insertion-tail split.".into());
        }
        text.parse()
            .map_err(|_| "Primer setting is too large.".into())
    }
    Ok(model::Options {
        forward_length: number(forward)?,
        reverse_length: number(reverse)?,
        reverse_tail_bases: number(split)?,
    })
}

/// Match the actual retained historical event, not an equal sequence in another document.
fn attest(draft: &Draft, frame: &Workspace) -> Result<EventBinding, String> {
    let source = draft
        .slots
        .first()
        .and_then(|s| s.captured.as_ref())
        .ok_or("Pin the original circular source first.")?;
    let docs = DocumentHistory::restore_all(frame).map_err(|e| e.to_string())?;
    let doc = docs
        .iter()
        .find(|d| d.id == source.source.binding.document_id)
        .ok_or("The exact source document is unavailable.")?;
    let molecule = doc
        .revision_molecule(&source.source.binding.event_id)
        .map_err(|e| e.to_string())?;
    if molecule != source.molecule || molecule.revision() != source.source.binding.revision {
        return Err("The pinned original event no longer matches this preparation.".into());
    }
    Ok(source.source.binding.clone())
}

#[component]
pub fn InversePcrPreparationView(
    draft: RwSignal<Draft>,
    intent: Arc<IntendedReference>,
    epoch: RwSignal<u64>,
    snapshot: Callback<(), Result<Workspace, String>>,
    current: Memo<bool>,
    guard: Callback<(), Result<(), String>>,
    on_review: Callback<Candidate>,
    on_finalize: Callback<Candidate, Result<(), String>>,
) -> impl IntoView {
    // Presentation eligibility only; native request validation remains authoritative.
    let single_edit = Memo::new(move |_| {
        draft.with(|d| {
            d.alternatives.is_none()
                && d.circular
                && d.slots.len() == 1
                && ((d.edits.len() == 1 && d.substitutions.is_empty())
                    || (d.edits.is_empty() && d.substitutions.len() == 1))
                && d.junction_edits.is_empty()
                && d.slots[0].captured.as_ref().is_some_and(|c| {
                    c.molecule.topology == collagene_core::Topology::Circular
                        && d.slots[0].start == 0
                        && d.slots[0].end == c.molecule.sequence.len()
                        && !d.slots[0].reverse
                        && d.slots[0].digest_geometry.is_none()
                        && d.slots[0].anneal_bottom.is_none()
                })
                && d.methods
                    .rows
                    .first()
                    .is_some_and(|r| r.run_origin.is_none())
        })
    });
    let intent = StoredValue::new(intent);
    let forward = RwSignal::new("20".to_owned());
    let reverse = RwSignal::new("20".to_owned());
    let split = RwSignal::new("0".to_owned());
    let status = RwSignal::new(String::new());
    let pending = RwSignal::new(false);
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
    let settings = Memo::new(move |_| {
        (
            forward.get(),
            reverse.get(),
            split.get(),
            epoch.get(),
            current.get(),
        )
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
    let read_options = move || {
        options(
            &forward.get_untracked(),
            &reverse.get_untracked(),
            &split.get_untracked(),
        )
    };
    let prepare = move |_| {
        let result = (|| {
            guard.run(())?;
            if !current.get_untracked() {
                return Err("Review the current independent intent first.".to_owned());
            }
            let declared = draft.get_untracked();
            let source = attest(&declared, &snapshot.run(())?)?;
            let options = read_options()?;
            let (capture, molecule, native_options) = intent.with_value(|i| {
                model::request(&declared, i.as_ref(), &source, i.target(), &options)
            })?;
            Ok::<_, String>((source, options, capture, molecule, native_options))
        })();
        let (source, options, capture, molecule, native_options) = match result {
            Ok(value) => value,
            Err(error) => {
                status.set(error);
                return;
            }
        };
        let live = match CachedWorker::new() {
            Ok(value) => Rc::new(value),
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
        status.set("Proposing original-source inverse PCR and checking its physical steps…".into());
        let alive = alive.clone();
        spawn_local(async move {
            let revision = capture.source.revision.clone();
            let response = live
                .ask(Request {
                    id: token,
                    revision: revision.clone(),
                    job: AnalysisJob::InversePcr {
                        molecule,
                        target: intent.with_value(|i| i.target().clone()),
                        options: native_options,
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
                let declared = draft.get_untracked();
                let now = attest(&declared, &snapshot.run(())?)?;
                if now != source {
                    return Err("The source event changed while preparing inverse PCR.".into());
                }
                let response = response?;
                if response.id != token || response.revision != revision {
                    return Err("The inverse-PCR reply belongs to another request.".into());
                }
                let native = match response.result {
                    AnalysisResult::InversePcr(report) => report,
                    AnalysisResult::Error(error) => return Err(error),
                    _ => return Err("Unexpected inverse-PCR response.".into()),
                };
                let result = intent.with_value(|i| {
                    model::from_design(&declared, i.as_ref(), &now, i.target(), &options, native)
                })?;
                if result.capture != capture {
                    return Err("The inverse-PCR declaration changed.".into());
                }
                Ok::<_, String>(result)
            })();
            match checked {
                Ok(result) => {
                    status.set("Review both ordered primers and the proposed preparation, then open fresh candidate verification.".into());
                    proposal.set(Some(Arc::new(result)));
                }
                Err(error) => status.set(error),
            }
        });
    };
    let checked_candidate = Callback::new(move |saved: Arc<model::Proposal>| {
        guard.run(())?;
        if !current.get_untracked() {
            return Err("The independent intent is stale.".into());
        }
        let declared = draft.get_untracked();
        let source = attest(&declared, &snapshot.run(())?)?;
        let options = read_options()?;
        intent.with_value(|i| {
            saved.validate_binding(&declared, i.as_ref(), &source, i.target(), &options)
        })?;
        Ok::<_, String>(saved.candidate.clone())
    });
    view! {<details class="assembly-fragment-options"><summary>"Prepare an internal edit by inverse PCR"</summary>
        <p class="helper">"For one whole forward circular source and one declared edit. The proposed steps amplify the retained original DNA, phosphorylate its ends and close the circle. Your intended sequence stays fixed. More complex designs need explicitly staged preparations."</p>
        <Show when=move||!single_edit.get()><p class="helper">"This draft needs staged preparation. Direct inverse PCR requires a single whole forward circular contribution and a single internal edit."</p></Show>
        <div class="primer-grid">
            <label>"Forward annealing length"<input type="number" min="8" max="80" prop:value=move||forward.get() on:input=move|e|forward.set(event_target_value(&e))/></label>
            <label>"Reverse annealing length"<input type="number" min="8" max="80" prop:value=move||reverse.get() on:input=move|e|reverse.set(event_target_value(&e))/></label>
            <label>"Inserted bases on reverse primer"<input type="number" min="0" prop:value=move||split.get() on:input=move|e|split.set(event_target_value(&e))/></label>
        </div>
        <button class="button secondary" disabled=move||pending.get()||!current.get()||!single_edit.get() on:click=prepare>"Propose inverse-PCR preparation"</button>
        <Show when=move||pending.get()><button class="text-tool" on:click=move|_|{
            generation.update(|g|*g=g.wrapping_add(1));
            if let Some(worker)=worker.get_value(){worker.close();}
            worker.set_value(None); pending.set(false); proposal.set(None);
            status.set("Inverse-PCR preparation canceled.".into());
        }>"Cancel"</button></Show>
        <p role="status">{move||status.get()}</p>
        {move||proposal.get().map(|saved|{
            let invalid=saved.verification.status==Status::Invalid;
            let result = match saved.verification.status {
                Status::ConsistentUnderAssumptions => "Sequence and chemistry are consistent under the stated assumptions.",
                Status::Invalid => "The preparation does not satisfy the intended design.",
                Status::Incomplete => "The preparation could not be fully checked.",
            };
            let text=format!("{result} Review the assembly candidate for its current source-bound verification.");
            let saved=StoredValue::new(saved);
            view!{<section class="assembly-inverse-pcr-proposal"><p>{text}</p>
                <ol><li>"PCR from the pinned original circular template"</li><li>"Phosphorylate both amplicon 5′ ends"</li><li>"Ligate the amplicon into a circle"</li></ol>
                <table><thead><tr><th>"Primer"</th><th>"Ordered sequence · 5′→3′"</th></tr></thead><tbody>
                    <tr><td>"Forward"</td><td><pre>{saved.with_value(|p|p.preparation.forward.ordered_sequence())}</pre></td></tr>
                    <tr><td>"Reverse"</td><td><pre>{saved.with_value(|p|p.preparation.reverse.ordered_sequence())}</pre></td></tr>
                </tbody></table>
                {saved.with_value(|p|p.preparation.warnings.iter().take(20).map(|s|view!{<p class="helper">{s.clone()}</p>}).collect_view())}
                {saved.with_value(|p|(p.preparation.warnings.len()>20).then(||view!{<p class="helper">"Additional preparation warnings are retained for candidate review."</p>}))}
                <button class="button primary" disabled=move||!current.get() on:click=move|_|{
                    match checked_candidate.run(saved.get_value()) {
                        Ok(candidate)=>on_review.run(candidate),Err(error)=>status.set(error),
                    }
                }>"Review cloning candidate"</button>
                <button class="button secondary" disabled=move||invalid||!current.get() on:click=move|_|{
                    match checked_candidate.run(saved.get_value()).and_then(|c|on_finalize.run(c)) {
                        Ok(())=>status.set("Preparation submitted to the reviewed output workflow.".into()),Err(error)=>status.set(error),
                    }
                }>"Review outputs"</button>
                {invalid.then(||view!{<p class="dialog-error">"The proposed preparation is invalid against the intended sequence. Inspect its candidate findings before changing the design."</p>})}
            </section>}
        })}
    </details>}
}
