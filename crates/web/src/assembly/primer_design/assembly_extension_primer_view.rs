//! Compact ordinary-dock proposal UI. Worker results are diagnostic declarations;
//! installing six oligos never installs a product or current verification status.
use crate::assembly::draft::model::Draft;
use crate::assembly::primer_design::assembly_extension_primer_model as model;
use crate::assembly::review::assembly_review::IntendedReference;
use crate::platform::worker::CachedWorker;
use collagene_core::{
    cloning::overlap::extension::primers as native,
    documents::workspace::Workspace,
    protocol::{
        request::Request,
        schema::{AnalysisJob, AnalysisResult},
    },
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
fn integer(text: &str, name: &str) -> Result<usize, String> {
    if text.is_empty() || text.len() > 20 || !text.bytes().all(|b| b.is_ascii_digit()) {
        return Err(format!("Enter a whole number for {name}."));
    }
    text.parse()
        .map_err(|_| format!("{name} exceeds the supported range."))
}
fn temperature(text: &str) -> Result<f64, String> {
    if text.len() > 32 {
        return Err("The target annealing temperature is too long.".into());
    }
    let tm = text
        .parse::<f64>()
        .map_err(|_| "Enter a target annealing temperature.")?;
    if !tm.is_finite() {
        return Err("The target annealing temperature must be finite.".into());
    }
    Ok(tm)
}
#[component]
pub fn ExtensionPrimerProposalView(
    draft: RwSignal<Draft>,
    intent: Arc<IntendedReference>,
    epoch: RwSignal<u64>,
    snapshot: Callback<(), Result<Workspace, String>>,
    current: Memo<bool>,
    guard: Callback<(), Result<(), String>>,
    on_use: Callback<(Draft, Draft), Result<(), String>>,
) -> impl IntoView {
    let defaults = native::Options::default();
    let minimum = RwSignal::new(defaults.primer_options.min_length.to_string());
    let maximum = RwSignal::new(defaults.primer_options.max_length.to_string());
    let tm = RwSignal::new(defaults.primer_options.target_tm_c.to_string());
    let flank = RwSignal::new(defaults.homology_flank.to_string());
    let max_oligo = RwSignal::new(defaults.max_oligo_length.to_string());
    let min_paired = RwSignal::new(defaults.min_paired_bases.to_string());
    let comparisons = RwSignal::new(defaults.max_comparisons.to_string());
    let candidates = RwSignal::new(defaults.max_candidates.to_string());
    let pairs = RwSignal::new(defaults.primer_options.max_pairs.to_string());
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
    let intent = StoredValue::new(intent);
    let settings = Memo::new(move |_| {
        (
            minimum.get(),
            maximum.get(),
            tm.get(),
            flank.get(),
            max_oligo.get(),
            min_paired.get(),
            comparisons.get(),
            candidates.get(),
            pairs.get(),
            epoch.get(),
            current.get(),
            draft.get(),
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
    let eligible = Memo::new(move |_| {
        draft.with(|d| {
            if d.alternatives.is_some()
                || d.methods.mode == crate::cloning::preparation::cloning_methods_model::Mode::MultiSiteReference
            {
                return false;
            }
            let Some(slot) = d.slots.get(d.active_slot) else {
                return false;
            };
            let Some(row) = d.methods.rows.get(d.active_slot) else {
                return false;
            };
            if slot.captured.is_none()
                || slot.anneal_bottom.is_some()
                || row.preparation
                    == crate::cloning::preparation::cloning_methods_model::OverlapPreparation::AnnealedOligos
            {
                return false;
            }
            let Ok(roles) = d.methods.role_ids(d.slots.len()) else {
                return false;
            };
            let Some(role) = roles.get(d.active_slot) else {
                return false;
            };
            d.edits.iter().any(|edit| edit.role == role.role)
                || d.substitutions.iter().any(|edit| edit.role == role.role)
        })
    });
    let read_options = move |declared: &Draft| {
        let mut options = native::Options::default();
        options.primer_options.min_length =
            integer(&minimum.get_untracked(), "minimum annealing length")?;
        options.primer_options.max_length =
            integer(&maximum.get_untracked(), "maximum annealing length")?;
        options.primer_options.target_tm_c = temperature(&tm.get_untracked())?;
        options.primer_options.max_pairs =
            integer(&pairs.get_untracked(), "ranked pairs per native call")?;
        options.homology_flank = integer(&flank.get_untracked(), "unchanged flank length")?;
        options.max_oligo_length =
            integer(&max_oligo.get_untracked(), "maximum total oligo length")?;
        options.min_paired_bases =
            integer(&min_paired.get_untracked(), "minimum terminal pairing")?;
        options.max_comparisons =
            integer(&comparisons.get_untracked(), "terminal comparison budget")?;
        options.max_candidates =
            integer(&candidates.get_untracked(), "retained terminal-pair budget")?;
        let tails = model::current_outer_tails(declared, declared.active_slot)?;
        options.outer_forward_tail = tails.0;
        options.outer_reverse_tail = tails.1;
        Ok::<_, String>(options)
    };
    let prepare = move |_| {
        let captured = (|| {
            guard.run(())?;
            if !current.get_untracked() {
                return Err("Review the current independent edit intent first.".to_owned());
            }
            let baseline = draft.get_untracked();
            let workspace = snapshot.run(())?;
            let source = model::attest(&baseline, baseline.active_slot, &workspace)?;
            let options = read_options(&baseline)?;
            let (capture, molecule, request) = intent.with_value(|i| {
                model::request(
                    &baseline,
                    i.as_ref(),
                    baseline.active_slot,
                    &source,
                    &options,
                    Some(&workspace),
                )
            })?;
            Ok::<_, String>((
                baseline, source, options, capture, molecule, request, workspace,
            ))
        })();
        let (baseline, source, options, capture, molecule, request, workspace) = match captured {
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
        status.set(
            "Proposing native exact-source strand oligos and reconstructing the two PCR branches…"
                .into(),
        );
        let alive = alive.clone();
        spawn_local(async move {
            let revision = capture.design_revision.clone();
            let response = live
                .ask(Request {
                    id: token,
                    revision: revision.clone(),
                    job: AnalysisJob::DesignOverlapExtensionPrimers {
                        source: molecule,
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
                if now != baseline {
                    return Err("The assembly draft changed while proposing primers.".into());
                }
                let workspace = snapshot.run(())?;
                let event = model::attest(&now, now.active_slot, &workspace)?;
                if event != source {
                    return Err("The exact source event changed while proposing primers.".into());
                }
                let response = response?;
                if response.id != token || response.revision != revision {
                    return Err("The primer reply belongs to another request.".into());
                }
                let report = match response.result {
                    AnalysisResult::OverlapExtensionPrimersDesigned(report) => report,
                    AnalysisResult::Error(error) => return Err(error),
                    _ => return Err("Unexpected overlap-extension primer response.".into()),
                };
                let result = intent.with_value(|i| {
                    model::from_design(
                        &now,
                        i.as_ref(),
                        now.active_slot,
                        &event,
                        &options,
                        report,
                        Some(&workspace),
                    )
                })?;
                if result.capture != capture {
                    return Err("The original source/edit/options declaration changed.".into());
                }
                Ok::<_, String>(result)
            })();
            match checked {
                Ok(result) => {
                    status.set("Review all six ordered oligos. Using them changes this contribution’s preparation only; verify the complete assembly afterward.".into());
                    proposal.set(Some(Arc::new(result)));
                }
                Err(error) => status.set(error),
            }
        });
    };
    let use_proposal = Callback::new(move |saved: Arc<model::Proposal>| {
        if !proposal.with_untracked(|active| {
            active
                .as_ref()
                .is_some_and(|active| Arc::ptr_eq(active, &saved))
        }) {
            return Err("This primer proposal is no longer current.".into());
        }
        guard.run(())?;
        if !current.get_untracked() {
            return Err("The independent edit intent is stale.".into());
        }
        let before = draft.get_untracked();
        let workspace = snapshot.run(())?;
        let source = model::attest(&before, before.active_slot, &workspace)?;
        let options = read_options(&before)?;
        let next = intent.with_value(|i| {
            saved.apply(&before, i.as_ref(), &source, &options, Some(&workspace))
        })?;
        on_use.run((before, next))?;
        proposal.set(None);
        Ok::<_, String>(())
    });
    view! {<Show when=move||draft.with(|d|d.methods.mode!=crate::cloning::preparation::cloning_methods_model::Mode::MultiSiteReference)><details class="assembly-fragment-options"><summary>"Propose overlap-extension primers"</summary>
        <p class="helper">"Prepare the active contribution’s internal edits using two source PCRs, reciprocal extension and an outer PCR. Partial, reversed and origin-crossing contributions retain their original source coordinates. Your intended construct stays fixed."</p>
        <Show when=move||!eligible.get()><p class="helper">"Select a pinned contribution with an internal edit or substitution. Paired-oligo contributions retain their separate preparation."</p></Show>
        <div class="primer-grid">
            <label>"Minimum annealing length"<input type="number" min="8" max="80" prop:value=move||minimum.get() on:input=move|e|minimum.set(event_target_value(&e))/></label>
            <label>"Maximum annealing length"<input type="number" min="8" max="80" prop:value=move||maximum.get() on:input=move|e|maximum.set(event_target_value(&e))/></label>
            <label>"Target annealing Tm (°C)"<input type="number" min="-50" max="100" step="0.1" prop:value=move||tm.get() on:input=move|e|tm.set(event_target_value(&e))/></label>
            <label>"Unchanged flank on each side (bp)"<input type="number" min="1" max="200" prop:value=move||flank.get() on:input=move|e|flank.set(event_target_value(&e))/></label>
            <label>"Maximum total oligo length"<input type="number" min="8" max="200" prop:value=move||max_oligo.get() on:input=move|e|max_oligo.set(event_target_value(&e))/></label>
        </div>
        <details><summary>"Retained outer tails and search bounds"</summary>
            <p class="helper">"Current outer-PCR tail sequences are preserved using standard oligo case and whitespace normalization. Edit the preparation fields to change their bases explicitly."</p>
            <pre>{move||draft.with(|d|match model::current_outer_tails(d,d.active_slot){Ok((forward,reverse))=>format!("Forward 5′ tail: {}\nReverse 5′ tail: {}",if forward.is_empty(){"(none)"}else{&forward},if reverse.is_empty(){"(none)"}else{&reverse}),Err(error)=>error})}</pre>
            <div class="primer-grid">
                <label>"Ranked pairs per native call"<input type="number" min="1" max="20" prop:value=move||pairs.get() on:input=move|e|pairs.set(event_target_value(&e))/></label>
                <label>"Minimum terminal pairing (bp)"<input type="number" min="1" prop:value=move||min_paired.get() on:input=move|e|min_paired.set(event_target_value(&e))/></label>
                <label>"Terminal comparison budget"<input type="number" min="1" max="50000000" prop:value=move||comparisons.get() on:input=move|e|comparisons.set(event_target_value(&e))/></label>
                <label>"Retained terminal pairings"<input type="number" min="1" max="1000" prop:value=move||candidates.get() on:input=move|e|candidates.set(event_target_value(&e))/></label>
            </div>
        </details>
        <button class="button secondary" disabled=move||pending.get()||!current.get()||!eligible.get() on:click=prepare.clone()>"Propose primers for active contribution"</button>
        <Show when=move||pending.get()><button class="text-tool" on:click=move|_|{generation.update(|g|*g=g.wrapping_add(1));if let Some(worker)=worker.get_value(){worker.close();}worker.set_value(None);pending.set(false);proposal.set(None);status.set("Primer proposal canceled.".into());}>"Cancel"</button></Show>
        <p role="status">{move||status.get()}</p>
        {move||proposal.get().map(|saved|{
            let invalid=saved.verification.status==Status::Invalid;
            let outcome=format!("Native preparation diagnostic: {:?}. The complete assembly still requires fresh ordinary verification, including any saved material origin.",saved.verification.status);
            let saved=StoredValue::new(saved);
            view!{<section class="assembly-extension-primer-proposal"><p>{outcome}</p>
                <ol><li>"Left PCR from the pinned original source"</li><li>"Right PCR from the same original source"</li><li>"Reciprocal extension of the actual PCR strands"</li><li>"Outer-primer PCR"</li></ol>
                <table><thead><tr><th>"Pair"</th><th>"Primer"</th><th>"Ordered sequence · 5′→3′"</th></tr></thead><tbody>
                    {saved.with_value(|p|[(&p.preparation.left,"Left source PCR"),(&p.preparation.right,"Right source PCR"),(&p.preparation.outer,"Outer PCR")].into_iter().flat_map(|(pair,label)|[(&pair.forward,"Forward"),(&pair.reverse,"Reverse")].into_iter().map(move|(primer,direction)|view!{<tr><td>{label}</td><td>{format!("{direction} · {}",primer.name)}</td><td><pre>{primer.ordered_sequence()}</pre></td></tr>})).collect_view())}
                </tbody></table>
                {saved.with_value(|p|p.warnings.iter().take(20).map(|warning|view!{<p class="helper">{warning.clone()}</p>}).collect_view())}
                <details><summary>"Native diagnostic findings"</summary>{saved.with_value(|p| crate::protocol::verification_view::report(&p.verification))}</details>
                <button class="button primary" disabled=move||invalid||!current.get() on:click=move|_|{match use_proposal.run(saved.get_value()){Ok(())=>status.set("Oligos entered for the selected preparation. Reconstruct and verify the complete cloning plan.".into()),Err(error)=>status.set(error)}}>"Use oligos in this preparation"</button>
            </section>}
        })}
    </details></Show>}
}
