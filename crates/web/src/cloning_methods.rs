//! Explicit preparation graphs for TA/GC ligation and oriented overlap profiles.
use crate::{app::download, golden_gate_pool_form::PoolForm, worker::CachedWorker};
use collagene_core::{
    Molecule,
    analysis::Primer,
    chemistry::{EndRepairMethod, EndTreatment, FivePrime, Terminus, ThreePrime},
    overlap_profiles::{OverlapMethod, ProfileOptions, ProfileReport},
    protocol::{AnalysisJob, AnalysisResult, Request},
    tailing::TailingBase,
    transfer::Candidate,
    verification::{
        DesignPlan, Operation, OrientedSource, PlanInput, Source, TailPreview, TailSubstrate,
    },
};
use leptos::{html, prelude::*};
use serde::{Deserialize, Serialize};
use std::{
    rc::Rc,
    sync::{
        Arc,
        atomic::{AtomicBool, Ordering},
    },
};
use wasm_bindgen::JsCast;
use wasm_bindgen_futures::spawn_local;

pub use crate::cloning_methods_model::*;

#[derive(Clone, PartialEq, Eq, Serialize, Deserialize)]
struct PlannerForm {
    overlap: String,
    minimum: String,
    maximum: String,
    tm: String,
}
impl Default for PlannerForm {
    fn default() -> Self {
        Self {
            overlap: "20".into(),
            minimum: "18".into(),
            maximum: "30".into(),
            tm: "60".into(),
        }
    }
}
/// Actual document identity or an explicitly reviewed reference intention.
/// An intent has no invented document/event ancestry.
#[derive(Clone, PartialEq, Eq, Serialize, Deserialize)]
pub enum TargetContext {
    Document(crate::digest::SourceInfo),
    Intent {
        context_revision: String,
        revision: String,
    },
}
impl TargetContext {
    pub fn revision(&self) -> &str {
        match self {
            Self::Document(info) => &info.binding.revision,
            Self::Intent { revision, .. } => revision,
        }
    }
}
struct Suggestions {
    report: collagene_core::overlap_primer_design::Report,
    config: MethodsState,
    form: PlannerForm,
    context: Vec<crate::digest::SourceInfo>,
    target: String,
    target_context: TargetContext,
    epoch: u64,
}
struct GoldenGateSuggestions {
    report: collagene_core::golden_gate_primer_design::Report,
    config: MethodsState,
    context: Vec<crate::digest::SourceInfo>,
    target_context: TargetContext,
    epoch: u64,
}
#[derive(Clone)]
struct Evidence {
    config: MethodsState,
    revisions: Vec<String>,
    context: Vec<crate::digest::SourceInfo>,
    target_context: TargetContext,
    target: String,
    epoch: u64,
    tails: Vec<TailPreview>,
    profile: Option<ProfileReport>,
    overlap: Option<collagene_core::assembly_recipe::Preview>,
    ligation: Option<collagene_core::ligation_recipe::Preview>,
    tailing: Option<collagene_core::tailing_recipe::Preview>,
}
fn preparation_label(mode: Mode) -> &'static str {
    match mode {
        Mode::Tailing => "Tailing",
        Mode::Ligation => "Ligation",
        Mode::GoldenGate => "Golden Gate",
        Mode::Overlap => "Overlap",
        Mode::MultiSiteReference => "Reference",
    }
}
fn guidance(minimum: Option<usize>, maximum: Option<usize>) -> String {
    match (minimum, maximum) {
        (Some(a), Some(b)) => format!("{a}–{b} bp"),
        (Some(a), None) => format!("at least {a} bp; no upper bound specified"),
        _ => "not established".into(),
    }
}
fn short(sequence: &str) -> String {
    let s = sequence.chars().take(300).collect::<String>();
    if s.len() < sequence.len() {
        format!("{s}…")
    } else {
        s
    }
}
fn append_preparation(
    steps: &mut Vec<Operation>,
    input_id: &str,
    input: &TailInput,
) -> Result<Source, String> {
    let source = Source::Input(input_id.into());
    match input.recipe()? {
        TailSubstrate::Digest {
            enzymes,
            fragment_index,
            blunt_before,
        } => {
            steps.push(Operation::Digest {
                template: source,
                enzymes,
                fragment_index,
            });
            if blunt_before {
                let template = Source::Step(steps.len() - 1);
                steps.push(Operation::RepairEnds {
                    template,
                    method: EndRepairMethod::Blunt,
                });
            }
        }
        TailSubstrate::ExactPcr { forward, reverse } => steps.push(Operation::Pcr {
            template: source,
            forward,
            reverse,
        }),
    };
    let template = Source::Step(steps.len() - 1);
    steps.push(Operation::Tail3Prime {
        template,
        base: input.base,
    });
    let mut source = Source::Step(steps.len() - 1);
    let treatment = match input.treatment {
        Treatment::Keep => None,
        Treatment::Phosphorylate => Some(EndTreatment::Phosphorylate),
        Treatment::Dephosphorylate => Some(EndTreatment::Dephosphorylate),
    };
    if let Some(treatment) = treatment {
        steps.push(Operation::TreatEnds {
            template: source,
            treatment,
        });
        source = Source::Step(steps.len() - 1);
    }
    Ok(source)
}
async fn ask(
    worker: &CachedWorker,
    job: AnalysisJob,
    revision: String,
    id: u32,
) -> Result<AnalysisResult, String> {
    let response = worker.ask(Request { id, revision, job }).await?;
    match response.result {
        AnalysisResult::Error(e) => Err(e),
        result => Ok(result),
    }
}
/// Full own-5′→3′ domains are deliberately inside an explicit review disclosure.
fn proposed_oligo(direction: &'static str, primer: &Primer) -> impl IntoView + use<> {
    let name = format!("{} · {}", direction, primer.name);
    let ordered = format!("{}{}", primer.tail, primer.annealing.as_str());
    let annealing = primer.annealing.as_str().to_owned();
    let tail = if primer.tail.is_empty() {
        "None".to_owned()
    } else {
        primer.tail.clone()
    };
    view! {<section class="primer-review-oligo"><strong>{name}</strong><dl>
        <dt>"Ordered oligo · 5′→3′"</dt><dd><code>{ordered}</code></dd>
        <dt>"Annealing domain · 5′→3′"</dt><dd><code>{annealing}</code></dd>
        <dt>"Added 5′ tail · 5′→3′"</dt><dd><code>{tail}</code></dd>
    </dl></section>}
}

#[component]
pub fn CloningMethodsPanel(
    molecule: Memo<Molecule>,
    query_options: Memo<Vec<(usize, String)>>,
    get_query: Callback<usize, Option<Molecule>>,
    source_info: Callback<usize, Option<crate::digest::SourceInfo>>,
    target_info: Memo<TargetContext>,
    selected_region: Callback<(), Option<(usize, usize, usize)>>,
    config: RwSignal<MethodsState>,
    epoch: RwSignal<u64>,
    on_candidate: Callback<Candidate>,
    on_products: Callback<Vec<Molecule>, Result<(), String>>,
    primer_sources: Memo<Vec<crate::digest::SourceInfo>>,
    list_primers: Callback<
        collagene_core::run_lineage::EventBinding,
        Result<collagene_core::primer_reuse::List, String>,
    >,
    resolve_primer: Callback<
        crate::primer_picker::ResolveRequest,
        Result<collagene_core::primer_reuse::Resolved, String>,
    >,
    #[prop(optional_no_strip)] get_workspace: Option<
        Callback<(), Result<collagene_core::workspace::Workspace, String>>,
    >,
    #[prop(optional_no_strip)] on_finalize: Option<Callback<Candidate, Result<(), String>>>,
    #[prop(optional_no_strip)] on_source: Option<Callback<(usize, usize), Result<(), String>>>,
    // View state is separate from preparation/verification identity.
    #[prop(optional_no_strip)] active_row: Option<Memo<usize>>,
    #[prop(optional_no_strip)] on_edit_row: Option<Callback<usize>>,
    #[prop(default = false)] fixed_sources: bool,
    #[prop(default = false)] disable_plain_primer_design: bool,
) -> impl IntoView {
    let row_nodes = StoredValue::new_local(
        (0..12)
            .map(|_| NodeRef::<html::Div>::new())
            .collect::<Vec<_>>(),
    );
    let focus_row = Callback::new(move |index: usize| {
        if index >= config.with_untracked(|c| c.rows.len()).min(12) {
            return;
        }
        if let Some(select) = on_edit_row {
            select.run(index);
        }
        // The selected editor mounts in a reactive update. Focus after that
        // update; an old frame must not focus a different source/workspace.
        let expected_epoch = epoch.get_untracked();
        let expected_config = config.with_untracked(MethodsState::preparation_context);
        let frame = wasm_bindgen::closure::Closure::once_into_js(move || {
            if epoch.try_get_untracked() != Some(expected_epoch)
                || config.try_with_untracked(MethodsState::preparation_context)
                    != Some(expected_config)
                || active_row.is_some_and(|row| row.try_get_untracked() != Some(index))
            {
                return;
            }
            if let Some(node) = row_nodes
                .try_with_value(|nodes| nodes.get(index).and_then(|n| n.get()))
                .flatten()
            {
                node.scroll_into_view();
                let _ = node.focus();
            }
        });
        if let Some(window) = web_sys::window() {
            let _ = window.request_animation_frame(frame.unchecked_ref());
        }
    });
    let source_indices = Memo::new(move |_| {
        config.with(|c| {
            if c.mode == Mode::Tailing && !fixed_sources {
                Some(vec![c.vector.document, c.insert.document])
            } else {
                c.source_indices().ok()
            }
        })
    });
    let source_context = Memo::new(move |_| {
        source_indices.get().and_then(|indices| {
            indices
                .into_iter()
                .map(|i| source_info.run(i))
                .collect::<Option<Vec<_>>>()
        })
    });
    let sources = Memo::new(move |_| {
        source_indices.get().and_then(|indices| {
            indices
                .into_iter()
                .map(|i| get_query.run(i))
                .collect::<Option<Vec<_>>>()
        })
    });
    let source_revisions = Memo::new(move |_| {
        sources.with(|s| {
            s.as_ref()
                .map(|m| m.iter().map(Molecule::revision).collect::<Vec<_>>())
        })
    });
    let target = Memo::new(move |_| molecule.with(Molecule::revision));
    let busy = RwSignal::new(false);
    let status = RwSignal::new(String::new());
    let generation = RwSignal::new(0u32);
    let evidence = RwSignal::new(None::<Arc<Evidence>>);
    let planner_form = RwSignal::new(PlannerForm::default());
    let suggestions = RwSignal::new(None::<Arc<Suggestions>>);
    let gate_suggestions = RwSignal::new(None::<Arc<GoldenGateSuggestions>>);

    let active = StoredValue::new_local(None::<Rc<CachedWorker>>);
    let alive = Arc::new(AtomicBool::new(true));
    let cleanup = alive.clone();
    on_cleanup(move || {
        cleanup.store(false, Ordering::Relaxed);
        if let Some(w) = active.get_value() {
            w.close();
        }
        active.set_value(None);
    });
    Effect::new(move |_| {
        let _ = (
            config.get(),
            source_revisions.get(),
            source_context.get(),
            target_info.get(),
            planner_form.get(),
            target.get(),
            epoch.get(),
        );
        generation.update(|n| *n = n.wrapping_add(1));
        evidence.set(None);
        suggestions.set(None);
        gate_suggestions.set(None);
        if let Some(w) = active.get_value() {
            w.close();
        }
        active.set_value(None);
        busy.set(false);
        status.set(String::new());
    });
    let gate_alive = alive.clone();
    let suggest_gate = Callback::new(move |()| {
        if busy.get_untracked() {
            return;
        }
        let state = config.get_untracked();
        if state.mode != Mode::GoldenGate || disable_plain_primer_design {
            status.set("Use the junction editor to preserve the active edited cores when proposing primers.".into());
            return;
        }
        let (Some(originals), Some(context)) =
            (sources.get_untracked(), source_context.get_untracked())
        else {
            status.set("Choose the original source documents.".into());
            return;
        };
        let target_molecule = molecule.get_untracked();
        let target_context = target_info.get_untracked();
        let captured_epoch = epoch.get_untracked();
        let (inputs, request) = match golden_gate_request(&state, originals, &target_molecule) {
            Ok(value) => value,
            Err(error) => {
                status.set(error);
                return;
            }
        };
        let revisions = inputs
            .iter()
            .map(|p| p.revision.clone())
            .collect::<Vec<_>>();
        if revisions
            != context
                .iter()
                .map(|s| s.binding.revision.clone())
                .collect::<Vec<_>>()
            || target_context.revision() != target_molecule.revision()
        {
            status.set("Source or target identity changed before primer planning.".into());
            return;
        }
        let revision = match collagene_core::golden_gate_primer_design::preview_revision(
            &inputs,
            &target_molecule,
            &request,
        ) {
            Ok(value) => value,
            Err(error) => {
                status.set(error.to_string());
                return;
            }
        };
        let worker = match CachedWorker::new() {
            Ok(worker) => Rc::new(worker),
            Err(error) => {
                status.set(error);
                return;
            }
        };
        if let Some(old) = active.get_value() {
            old.close();
        }
        active.set_value(Some(worker.clone()));
        let token = generation.get_untracked().wrapping_add(1);
        generation.set(token);
        busy.set(true);
        evidence.set(None);
        suggestions.set(None);
        gate_suggestions.set(None);
        status.set(
            "Proposing explicit Type IIS PCR tails and freshly checking their digested ends…"
                .into(),
        );
        let alive = gate_alive.clone();
        spawn_local(async move {
            let result = ask(
                &worker,
                AnalysisJob::DesignGoldenGatePrimers {
                    inputs,
                    target: target_molecule,
                    request,
                },
                revision.clone(),
                token,
            )
            .await;
            worker.close();
            if !alive.load(Ordering::Relaxed)
                || generation.get_untracked() != token
                || config.get_untracked() != state
                || source_context.get_untracked() != Some(context.clone())
                || source_revisions.get_untracked() != Some(revisions.clone())
                || target_info.get_untracked() != target_context
                || target.get_untracked() != target_context.revision()
                || epoch.get_untracked() != captured_epoch
            {
                return;
            }
            active.set_value(None);
            busy.set(false);
            match result {
                Ok(AnalysisResult::GoldenGatePrimersDesigned(report)) => {
                    match report.validate_binding(&revision) {
                        Ok(()) if report.source_revisions == revisions => {
                            gate_suggestions.set(Some(Arc::new(GoldenGateSuggestions {
                                report,
                                config: state,
                                context,
                                target_context,
                                epoch: captured_epoch,
                            })));
                            status.set("Type IIS primer proposals ready. Review tails and actual digested ends, then explicitly apply.".into());
                        }
                        Ok(()) => status
                            .set("Primer proposals belong to different original sources.".into()),
                        Err(error) => status.set(error.to_string()),
                    }
                }
                Ok(AnalysisResult::Error(error)) | Err(error) => status.set(error),
                _ => status.set("Unexpected Golden Gate primer response.".into()),
            }
        });
    });
    let apply_gate = Callback::new(move |()| {
        let Some(saved) = gate_suggestions.get_untracked() else {
            return;
        };
        if config.get_untracked() != saved.config
            || source_context.get_untracked() != Some(saved.context.clone())
            || target_info.get_untracked() != saved.target_context
            || epoch.get_untracked() != saved.epoch
        {
            status.set("Sources or settings changed; request fresh Type IIS primers.".into());
            return;
        }
        let mut next = saved.config.clone();
        let roles = match next.role_ids(next.rows.len()) {
            Ok(roles) => roles,
            Err(error) => {
                status.set(error);
                return;
            }
        };
        for fragment in &saved.report.fragments {
            let Some(index) = roles
                .iter()
                .position(|role| fragment.role == role.role && fragment.input_id == role.input_id)
            else {
                status.set("Proposal has an unknown source role; nothing applied.".into());
                return;
            };
            let collagene_core::assembly_recipe::Preparation::PcrCatalogDigest {
                forward,
                reverse,
                selection,
                fragment_index,
                repair,
            } = &fragment.preparation
            else {
                status.set("Proposal is not an explicit PCR/digest chain; nothing applied.".into());
                return;
            };
            let row = &mut next.rows[index];
            if row.preparation != OverlapPreparation::PcrCatalogDigest
                || row.reverse != fragment.reverse
            {
                status.set("Preparation or orientation changed; nothing applied.".into());
                return;
            }
            row.forward_name = forward.name.clone();
            row.forward = forward.annealing.as_str().into();
            row.forward_tail = forward.tail.clone();
            row.reverse_name = reverse.name.clone();
            row.reverse_primer = reverse.annealing.as_str().into();
            row.reverse_tail = reverse.tail.clone();
            row.catalog_selection = Some(selection.clone());
            row.fragment = (fragment_index + 1).to_string();
            row.repair = *repair;
        }
        config.set(next);
        gate_suggestions.set(None);
        evidence.set(None);
        status.set("Proposed PCR/digest preparations applied explicitly. Existing-site rows were preserved. Compute method evidence to check actual junctions and the intended target.".into());
    });
    let suggestion_alive = alive.clone();
    let suggest = Callback::new(move |()| {
        let state = config.get_untracked();
        let selected = planner_form.get_untracked();
        if state.mode != Mode::Overlap || disable_plain_primer_design {
            status.set("Choose overlap assembly before suggesting primers.".into());
            return;
        }
        let (Some(inputs), Some(context)) =
            (sources.get_untracked(), source_context.get_untracked())
        else {
            status.set("Choose available original source documents.".into());
            return;
        };
        let bound_target = molecule.get_untracked();
        let target_context = target_info.get_untracked();
        let bound_epoch = epoch.get_untracked();
        let (inputs, request) =
            match primer_design_request(&state, inputs, &bound_target, &selected) {
                Ok(v) => v,
                Err(e) => {
                    status.set(e);
                    return;
                }
            };
        let revisions = inputs
            .iter()
            .map(|p| p.revision.clone())
            .collect::<Vec<_>>();
        if revisions
            != context
                .iter()
                .map(|s| s.binding.revision.clone())
                .collect::<Vec<_>>()
        {
            status.set("Original sources changed before primer planning.".into());
            return;
        }
        let revision = match collagene_core::overlap_primer_design::design_revision(
            &inputs,
            &bound_target,
            &request,
        ) {
            Ok(v) => v,
            Err(e) => {
                status.set(e.to_string());
                return;
            }
        };
        let target_revision = bound_target.revision();
        if target_revision != target_context.revision() {
            status.set("Target changed before primer planning.".into());
            return;
        }
        if source_context.get_untracked() != Some(context.clone())
            || target_info.get_untracked() != target_context
            || config.get_untracked() != state
            || planner_form.get_untracked() != selected
            || epoch.get_untracked() != bound_epoch
        {
            status.set("Sources or settings changed before primer planning.".into());
            return;
        }
        let worker = match CachedWorker::new() {
            Ok(v) => Rc::new(v),
            Err(e) => {
                status.set(e);
                return;
            }
        };
        if let Some(old) = active.get_value() {
            old.close();
        }
        active.set_value(Some(worker.clone()));
        let token = generation.get_untracked().wrapping_add(1);
        generation.set(token);
        busy.set(true);
        evidence.set(None);
        suggestions.set(None);
        gate_suggestions.set(None);
        status.set("Suggesting deterministic annealing domains and homology tails from original source regions…".into());
        let alive = suggestion_alive.clone();
        spawn_local(async move {
            let response = ask(
                &worker,
                AnalysisJob::DesignOverlapPrimers {
                    inputs,
                    target: bound_target,
                    request,
                },
                revision.clone(),
                token,
            )
            .await;
            worker.close();
            if !alive.load(Ordering::Relaxed)
                || generation.get_untracked() != token
                || config.get_untracked() != state
                || planner_form.get_untracked() != selected
                || source_context.get_untracked() != Some(context.clone())
                || target_info.get_untracked() != target_context
                || target.get_untracked() != target_revision
                || epoch.get_untracked() != bound_epoch
            {
                return;
            }
            active.set_value(None);
            busy.set(false);
            match response {
                Ok(AnalysisResult::OverlapPrimersDesigned(report)) => match report
                    .validate_binding(&revision)
                {
                    Ok(()) if report.source_revisions == revisions => {
                        suggestions.set(Some(Arc::new(Suggestions {
                            report,
                            config: state,
                            form: selected,
                            context,
                            target: target_revision,
                            target_context,
                            epoch: bound_epoch,
                        })));
                        status.set("Primer proposals ready. Review ordered sequences and explicitly apply before computing assembly evidence.".into());
                    }
                    Ok(()) => status.set(
                        "Primer suggestions did not match the original source revisions.".into(),
                    ),
                    Err(error) => status.set(error.to_string()),
                },
                Ok(_) => status.set("Unexpected overlap primer design response.".into()),
                Err(error) => status.set(error),
            }
        });
    });
    let apply_suggestions = Callback::new(move |()| {
        let Some(proposal) = suggestions.get_untracked() else {
            return;
        };
        if busy.get_untracked()
            || config.get_untracked() != proposal.config
            || planner_form.get_untracked() != proposal.form
            || source_context.get_untracked() != Some(proposal.context.clone())
            || target_info.get_untracked() != proposal.target_context
            || target.get_untracked() != proposal.target
            || epoch.get_untracked() != proposal.epoch
        {
            status.set("Sources, target or settings changed; request fresh primer suggestions before applying.".into());
            return;
        }
        if proposal.report.recipe.fragments.len() != proposal.config.rows.len() {
            status.set("Suggested primer rows do not match the original fragment order.".into());
            return;
        }
        let mut updated = proposal.config.clone();
        let roles = match updated.role_ids(updated.rows.len()) {
            Ok(roles) => roles,
            Err(error) => {
                status.set(error);
                return;
            }
        };
        for (i, fragment) in proposal.report.recipe.fragments.iter().enumerate() {
            if fragment.input_id != roles[i].input_id
                || fragment.role != roles[i].role
                || fragment.reverse != updated.rows[i].reverse
            {
                status
                    .set("Suggested source order or orientation changed; nothing applied.".into());
                return;
            }
            let collagene_core::assembly_recipe::Preparation::ExactPcr { forward, reverse } =
                &fragment.preparation
            else {
                status
                    .set("Suggestion is not an explicit PCR preparation; nothing applied.".into());
                return;
            };
            let row = &mut updated.rows[i];
            row.preparation = OverlapPreparation::ExactPcr;
            row.forward_name = forward.name.clone();
            row.forward = forward.annealing.as_str().into();
            row.forward_tail = forward.tail.clone();
            row.reverse_name = reverse.name.clone();
            row.reverse_primer = reverse.annealing.as_str().into();
            row.reverse_tail = reverse.tail.clone();
        }
        updated.overlaps = proposal
            .report
            .recipe
            .overlaps
            .iter()
            .map(usize::to_string)
            .collect::<Vec<_>>()
            .join(",");
        config.set(updated);
        suggestions.set(None);
        gate_suggestions.set(None);
        evidence.set(None);
        status.set("Suggested primers copied explicitly. Compute method evidence and verify the prepared assembly against the independent target.".into());
    });
    let compute = Callback::new(move |()| {
        let state = config.get_untracked();
        if state.mode == Mode::MultiSiteReference {
            status.set(
                "Use the ordinary assembly reference-exchange review against its frozen intent."
                    .into(),
            );
            return;
        }
        if state.mode == Mode::GoldenGate
            && state.source_roles.is_some()
            && state.golden_gate_pool.is_some()
        {
            status.set("This construct declares a finite reaction pool. Use Review selected reaction and junctions in its Construct reaction candidate panel, or explicitly return to ordered prepared assembly.".into());
            return;
        }
        let Some(inputs) = sources.get_untracked() else {
            status.set("Choose available source documents".into());
            return;
        };
        let Some(bound_context) = source_context.get_untracked() else {
            status.set("Choose original source events before computing method evidence.".into());
            return;
        };
        if state.rows.iter().enumerate().any(|(index, row)| {
            row.run_origin.as_ref().is_some_and(|pin| {
                bound_context
                    .get(index)
                    .is_none_or(|source| source.binding != pin.selection.event)
            })
        }) {
            status.set("A saved-material origin no longer matches its exact source event. Select it again before preparing fragments.".into());
            return;
        }
        let bound_target_context = target_info.get_untracked();
        if source_revisions.get_untracked()
            != Some(
                bound_context
                    .iter()
                    .map(|info| info.binding.revision.clone())
                    .collect(),
            )
            || target.get_untracked() != bound_target_context.revision()
        {
            status.set("Source or target event revisions differ from their current models. Select the current documents again.".into());
            return;
        }

        let recipes = if state.mode == Mode::Tailing && !fixed_sources {
            match state
                .vector
                .recipe()
                .and_then(|v| state.insert.recipe().map(|i| vec![v, i]))
            {
                Ok(r) => r,
                Err(e) => {
                    status.set(e);
                    return;
                }
            }
        } else {
            vec![]
        };
        let profile = if state.mode == Mode::Overlap {
            match state.profile() {
                Ok(p) => Some(p),
                Err(e) => {
                    status.set(e);
                    return;
                }
            }
        } else {
            None
        };
        let worker = match CachedWorker::new() {
            Ok(w) => Rc::new(w),
            Err(e) => {
                status.set(e);
                return;
            }
        };
        if let Some(old) = active.get_value() {
            old.close();
        }
        active.set_value(Some(worker.clone()));
        let token = generation.get_untracked().wrapping_add(1);
        generation.set(token);
        busy.set(true);
        evidence.set(None);
        status.set("Reconstructing method evidence…".into());
        let bound_target = target.get_untracked();
        let bound_target_molecule = molecule.get_untracked();
        let bound_epoch = epoch.get_untracked();
        let alive = alive.clone();
        spawn_local(async move {
            let revs = inputs.iter().map(Molecule::revision).collect::<Vec<_>>();
            let result = async {
                if state.mode == Mode::Tailing && fixed_sources {
                    let (bound_inputs, recipe) =
                        tailing_request(&state, inputs, &bound_target_molecule)?;
                    let revision = collagene_core::tailing_recipe::preview_revision(
                        &bound_inputs,
                        &bound_target_molecule,
                        &recipe,
                    )
                    .map_err(|e| e.to_string())?;
                    let response = ask(
                        &worker,
                        AnalysisJob::PreviewTailingRecipe {
                            inputs: bound_inputs,
                            target: bound_target_molecule,
                            recipe,
                        },
                        revision.clone(),
                        token,
                    )
                    .await?;
                    if let AnalysisResult::TailingRecipePreview(preview) = response {
                        preview
                            .validate_binding(&revision)
                            .map_err(|e| e.to_string())?;
                        Ok((vec![], None, None, None, Some(preview)))
                    } else {
                        Err("Unexpected tailing/ligation response".into())
                    }
                } else if matches!(state.mode, Mode::Ligation | Mode::GoldenGate) {
                    let (bound_inputs, recipe) =
                        ligation_request(&state, inputs, &bound_target_molecule)?;
                    let revision = collagene_core::ligation_recipe::preview_revision(
                        &bound_inputs,
                        &bound_target_molecule,
                        &recipe,
                    )
                    .map_err(|e| e.to_string())?;
                    let response = ask(
                        &worker,
                        AnalysisJob::PreviewLigationRecipe {
                            inputs: bound_inputs,
                            target: bound_target_molecule,
                            recipe,
                        },
                        revision.clone(),
                        token,
                    )
                    .await?;
                    if let AnalysisResult::LigationRecipePreview(preview) = response {
                        preview
                            .validate_binding(&revision)
                            .map_err(|e| e.to_string())?;
                        Ok((vec![], None, None, Some(preview), None))
                    } else {
                        Err("Unexpected ordered ligation response".into())
                    }
                } else if let Some((overlaps, options)) = profile {
                    let _ = (overlaps, options); // shared recipe reuses the exact parsed profile.
                    let (bound_inputs, recipe) =
                        overlap_request(&state, inputs, &bound_target_molecule)?;
                    let revision = collagene_core::assembly_recipe::preview_revision(
                        &bound_inputs,
                        &bound_target_molecule,
                        &recipe,
                    )
                    .map_err(|e| e.to_string())?;
                    let response = ask(
                        &worker,
                        AnalysisJob::PreviewOverlapRecipe {
                            inputs: bound_inputs,
                            target: bound_target_molecule,
                            recipe,
                        },
                        revision.clone(),
                        token,
                    )
                    .await?;
                    if let AnalysisResult::OverlapRecipePreview(preview) = response {
                        preview
                            .validate_binding(&revision)
                            .map_err(|e| e.to_string())?;
                        Ok((vec![], preview.profile.clone(), Some(preview), None, None))
                    } else {
                        Err("Unexpected overlap preparation response".into())
                    }
                } else {
                    let mut previews = vec![];
                    for (i, (input, recipe)) in inputs.into_iter().zip(recipes).enumerate() {
                        let base = if i == 0 {
                            state.vector.base
                        } else {
                            state.insert.base
                        };
                        let result = ask(
                            &worker,
                            AnalysisJob::TailPreview {
                                molecule: input,
                                recipe,
                                base,
                            },
                            revs[i].clone(),
                            token.wrapping_add(i as u32),
                        )
                        .await
                        .map_err(|e| {
                            format!("{}: {e}", if i == 0 { "Vector" } else { "Insert" })
                        })?;
                        if let AnalysisResult::TailPreview(preview) = result {
                            previews.push(preview);
                        } else {
                            return Err("Unexpected tailing preview response".into());
                        }
                    }
                    Ok((previews, None, None, None, None))
                }
            }
            .await;
            worker.close();
            if !alive.load(Ordering::Relaxed)
                || generation.get_untracked() != token
                || config.get_untracked() != state
                || source_revisions.get_untracked() != Some(revs.clone())
                || source_context.get_untracked() != Some(bound_context.clone())
                || target_info.get_untracked() != bound_target_context
                || target.get_untracked() != bound_target
                || epoch.get_untracked() != bound_epoch
            {
                return;
            }
            active.set_value(None);
            busy.set(false);
            match result {
                Ok((tails, profile, overlap, ligation, tailing)) => {
                    evidence.set(Some(Arc::new(Evidence {
                        config: state,
                        revisions: revs,
                        context: bound_context,
                        target_context: bound_target_context,
                        target: bound_target,
                        epoch: bound_epoch,
                        tails,
                        profile,
                        overlap,
                        ligation,
                        tailing,
                    })));
                    status.set("Evidence ready. Review assumptions and build a candidate against the independently selected target.".into());
                }
                Err(e) => status.set(e),
            }
        });
    });
    let diagnostic_guard = Callback::new(move |saved: Arc<Evidence>| -> Result<(), String> {
        if !evidence.with_untracked(|current| {
            current
                .as_ref()
                .is_some_and(|value| Arc::ptr_eq(value, &saved))
        }) || config.get_untracked() != saved.config
            || source_context.get_untracked() != Some(saved.context.clone())
            || target_info.get_untracked() != saved.target_context
            || source_revisions.get_untracked() != Some(saved.revisions.clone())
            || target.get_untracked() != saved.target
            || epoch.get_untracked() != saved.epoch
        {
            return Err("Sources or preparation settings changed. Compute fresh evidence before navigating its diagnostics.".into());
        }
        Ok(())
    });
    let review = Callback::new(move |()| {
        if config.with_untracked(|state| {
            state.mode == Mode::GoldenGate
                && state.source_roles.is_some()
                && state.golden_gate_pool.is_some()
        }) {
            status.set("Review this construct's explicit pool choice in the Construct reaction candidate panel; ordered evidence cannot replace that choice.".into());
            return;
        }
        let (Some(proof), Some(inputs)) = (evidence.get_untracked(), sources.get_untracked())
        else {
            return;
        };
        if let Err(error) = diagnostic_guard.run(proof.clone()) {
            status.set(error);
            return;
        }
        let state = config.get_untracked();
        let current_target = molecule.get_untracked();
        if proof.config != state
            || proof.revisions != inputs.iter().map(Molecule::revision).collect::<Vec<_>>()
            || proof.target != current_target.revision()
            || proof.epoch != epoch.get_untracked()
        {
            status.set(
                "Inputs changed; recompute method evidence before reviewing a candidate".into(),
            );
            return;
        }
        if state.mode == Mode::Tailing && fixed_sources {
            match tailing_request(&state, inputs, &current_target).and_then(|(bound, recipe)| {
                collagene_core::tailing_recipe::compile(&bound, &current_target, &recipe)
                    .map_err(|e| e.to_string())
            }) {
                Ok(compiled) => on_candidate.run(Candidate::new(compiled.plan)),
                Err(error) => status.set(error),
            }
            return;
        }
        let mut steps = vec![];
        let mut bound_inputs = vec![];
        if state.mode == Mode::Tailing {
            for (id, mol) in ["vector", "insert"].into_iter().zip(inputs) {
                bound_inputs.push(PlanInput {
                    id: id.into(),
                    revision: mol.revision(),
                    molecule: mol,
                });
            }
            let vector = match append_preparation(&mut steps, "vector", &state.vector) {
                Ok(s) => s,
                Err(e) => {
                    status.set(e);
                    return;
                }
            };
            let insert = match append_preparation(&mut steps, "insert", &state.insert) {
                Ok(s) => s,
                Err(e) => {
                    status.set(e);
                    return;
                }
            };
            steps.push(Operation::Ligate {
                fragments: vec![
                    OrientedSource {
                        source: vector,
                        reverse: state.vector.reversed,
                    },
                    OrientedSource {
                        source: insert,
                        reverse: state.insert.reversed,
                    },
                ],
                circular: state.circular,
                allow_nicked: state.allow_nicked,
            });
        } else if matches!(state.mode, Mode::Ligation | Mode::GoldenGate) {
            let (bound, recipe) = match ligation_request(&state, inputs, &current_target) {
                Ok(v) => v,
                Err(e) => {
                    status.set(e);
                    return;
                }
            };
            match collagene_core::ligation_recipe::compile(&bound, &current_target, &recipe) {
                Ok(compiled) => {
                    on_candidate.run(Candidate::new(compiled.plan));
                    return;
                }
                Err(error) => {
                    status.set(error.to_string());
                    return;
                }
            }
        } else {
            let (bound, recipe) = match overlap_request(&state, inputs, &current_target) {
                Ok(v) => v,
                Err(e) => {
                    status.set(e);
                    return;
                }
            };
            match collagene_core::assembly_recipe::compile(&bound, &current_target, &recipe) {
                Ok(compiled) => {
                    on_candidate.run(Candidate::new(compiled.plan));
                    return;
                }
                Err(error) => {
                    status.set(error.to_string());
                    return;
                }
            }
        }
        on_candidate.run(Candidate::new(DesignPlan {
            name: format!(
                "{} to {}",
                if state.mode == Mode::Tailing {
                    "Declared tailing and ligation"
                } else {
                    "Profiled overlap assembly"
                },
                current_target.name
            ),
            target_revision: current_target.revision(),
            inputs: bound_inputs,
            steps,
        }));
    });
    let reaction_context = Memo::new(move |_| {
        let state = config.with(MethodsState::preparation_context);
        if state.mode != Mode::GoldenGate {
            return None;
        }
        let context = source_context.get()?;
        let revisions = source_revisions.get()?;
        if state.rows.is_empty() || state.rows.len() > 12 || context.len() != state.rows.len() {
            return None;
        }
        state.role_ids(state.rows.len()).ok()?;
        Some(Arc::new(PoolContext {
            config: state,
            context,
            revisions,
            target_context: target_info.get(),
            target: target.get(),
            epoch: epoch.get(),
        }))
    });
    let reaction_guard = Callback::new(move |saved: Arc<PoolContext>| -> Result<(), String> {
        if config.with_untracked(MethodsState::preparation_context) != saved.config
            || source_context.get_untracked() != Some(saved.context.clone())
            || source_revisions.get_untracked() != Some(saved.revisions.clone())
            || target_info.get_untracked() != saved.target_context
            || target.get_untracked() != saved.target
            || epoch.get_untracked() != saved.epoch
        {
            return Err("Reaction source events, target or settings changed.".into());
        }
        Ok(())
    });
    let reaction_capture = Callback::new(move |saved: Arc<PoolContext>| {
        reaction_guard.run(saved.clone())?;
        let originals = sources
            .get_untracked()
            .ok_or("Choose available original reaction sources.")?;
        let roles = saved.config.role_ids(originals.len())?;
        let inputs = roles
            .into_iter()
            .zip(originals)
            .map(|(role, molecule)| PlanInput {
                id: role.input_id,
                revision: molecule.revision(),
                molecule,
            })
            .collect();
        Ok::<_, String>((inputs, molecule.get_untracked()))
    });
    view! {<div class="inspector-content cloning-methods-panel" class:active-fragment-preparation=active_row.is_some()><Show when=move||!fixed_sources><h2>"Cloning methods"</h2></Show><p class="helper">{if fixed_sources {"Verify the prepared fragments against the separately reviewed intended construct."} else {"The active document is the independently chosen target. Method evidence never selects a target or establishes experimental success."}}</p><Show when=move||!fixed_sources>
        <label>"Method workflow"<select aria-label="Cloning workflow" disabled=disable_plain_primer_design prop:value=move||match config.get().mode{Mode::Tailing=>"tailing",Mode::Overlap=>"overlap",Mode::Ligation=>"ligation",Mode::GoldenGate=>"golden_gate",Mode::MultiSiteReference=>"multisite"} on:change=move|e|config.update(|c|c.mode=match event_target_value(&e).as_str(){"overlap"=>Mode::Overlap,"ligation"=>Mode::Ligation,"golden_gate"=>Mode::GoldenGate,_=>Mode::Tailing})><option value="tailing">"TA / GC tailing and ligation"</option><option value="overlap">"Gibson / HiFi / In-Fusion overlap"</option><option value="ligation">"Ordered multi-fragment ligation"</option><option value="golden_gate">"Golden Gate — reaction and prepared assembly"</option></select></label></Show>
        <Show when=move||config.get().mode==Mode::Tailing><p class="helper">"Append one declared nucleotide to both actual 3′ ends of an established blunt duplex. This is an ideal endpoint operation; polymerase activity, mixtures and proprietary kit chemistry are not predicted. TA/GC products include two added junction bases."</p><div class="method-shortcuts"><button class="text-tool" on:click=move|_|config.update(|c|{c.vector.base=TailingBase::T;c.insert.base=TailingBase::A;})>"Use T/A pair"</button><button class="text-tool" on:click=move|_|config.update(|c|{c.vector.base=TailingBase::C;c.insert.base=TailingBase::G;})>"Use C/G pair"</button></div><Show when=move||active_row.is_none_or(|row|row.get()==0)><TailFields vector=true config=config query_options=query_options primer_sources list_primers resolve_primer epoch fixed_sources/></Show><Show when=move||active_row.is_none_or(|row|row.get()==1)><TailFields vector=false config=config query_options=query_options primer_sources list_primers resolve_primer epoch fixed_sources/></Show><label class="checkbox-label"><input type="checkbox" aria-label="Allow retained tailing ligation nicks" prop:checked=move||config.get().allow_nicked on:change=move|e|config.update(|c|c.allow_nicked=event_target_checked(&e))/>"Allow explicitly retained nicks"</label><p class="helper">"Exact PCR has 5′ OH termini. A sealed ligation may need an explicit phosphorylation step. Factory-vector chemistry and host nick repair are not inferred."</p></Show>
        <Show when=move||config.get().mode!=Mode::Tailing||fixed_sources><Show when=move||config.get().mode==Mode::Overlap><label>"Overlap profile"<select aria-label="Overlap profile method" prop:value=move||format!("{:?}",config.get().method) on:change=move|e|config.update(|c|c.method=match event_target_value(&e).as_str(){"Gibson"=>OverlapMethod::Gibson,"InFusion"=>OverlapMethod::InFusion,_=>OverlapMethod::HiFi})><option value="Gibson">"Gibson"</option><option value="HiFi">"HiFi"</option><option value="InFusion">"In-Fusion"</option></select></label><p class="helper">"Fragments are listed in product order with explicit orientation. Homology geometry is checked exactly; published length guidance is advisory. No reaction completion or unique mixture product is established."</p></Show>
            <Show when=move||matches!(config.get().mode,Mode::Ligation|Mode::GoldenGate)><p class="helper">"List fragments in intended ligation order. Preparation, end repair and terminal treatment are explicit source-bound operations. Whole linear DNA retains unknown terminal chemistry; its sequence alone does not imply cleaved or PCR ends. Unknown ends can withhold ligation even after phosphorylation."</p><label class="checkbox-label"><input type="checkbox" aria-label="Allow retained ordered ligation nicks" prop:checked=move||config.get().allow_nicked on:change=move|e|config.update(|c|c.allow_nicked=event_target_checked(&e))/>"Allow explicitly retained nicks"</label></Show>
            <h3>{if fixed_sources{"Prepare active fragment"}else{"1. Choose source fragments and regions"}}</h3>
            <Show when=move||config.with(|c|c.rows.len()>12)><p role="alert">"This preparation exceeds the supported 12 rows. No partial row list is shown."</p></Show><For each=move||config.with(|c|if c.rows.len()<=12{(0..c.rows.len()).collect::<Vec<_>>()}else{Vec::new()}) key=|i|*i children=move|i|view!{<Show when=move||active_row.is_none_or(|row|row.get()==i)><div node_ref=row_nodes.with_value(|nodes|nodes[i]) tabindex="-1" class="method-preparation-focus"><OverlapFields config index=i query_options primer_sources list_primers resolve_primer epoch selected_region fixed_sources edited_junctions=disable_plain_primer_design get_query source_info get_workspace/></div></Show>}/>
            <Show when=move||!fixed_sources><button class="button secondary full" disabled=move||{fixed_sources||config.with(|c|c.rows.len()>=12)} on:click=move|_|{let mut result=Ok(());config.update(|c|result=c.push_row(OverlapInput::default()));if let Err(error)=result{status.set(error);}}>{move||if matches!(config.get().mode,Mode::Ligation|Mode::GoldenGate){"Add ligation fragment"}else{"Add overlap fragment"}}</button><p class="helper">"This panel supports up to 12 ordered fragments; the typed core has its separate documented work bounds."</p></Show><Show when=move||config.get().mode==Mode::Overlap><label>"Overlap lengths in junction order (bp)"<input aria-label="Profile overlap lengths" prop:value=move||config.get().overlaps on:input=move|e|config.update(|c|c.overlaps=event_target_value(&e))/></label><label class="checkbox-label"><input type="checkbox" aria-label="Scan repeated overlap words" prop:checked=move||config.get().scan_repeats on:change=move|e|config.update(|c|c.scan_repeats=event_target_checked(&e))/>"Search supplied fragments for repeated homology words"</label><label>"Repeat comparison budget"<input type="number" min="1" aria-label="Overlap repeat budget" prop:value=move||config.get().comparisons on:input=move|e|config.update(|c|c.comparisons=event_target_value(&e))/></label><label>"Maximum reported homology occurrences"<input type="number" min="1" aria-label="Overlap occurrence limit" prop:value=move||config.get().occurrences on:input=move|e|config.update(|c|c.occurrences=event_target_value(&e))/></label></Show>
        </Show>
        <Show when=move||config.get().mode==Mode::GoldenGate&&!disable_plain_primer_design><section class="golden-gate-primer-planner"><h3>"2. Prepare Type IIS fragments"</h3>
            <p class="helper">"Use existing pinned digest fragments, or choose Exact PCR → pinned catalog digest for rows needing new sites. Proposals use each chosen original source region and explicit right fusion. They do not replace existing-site rows or choose a target. This ordered reconstruction does not enumerate every mixture outcome."</p>
            <details class="assembly-primer-settings"><summary>"Type IIS primer design settings"</summary><div class="primer-grid">
                <label>"Type IIS enzyme"<select aria-label="Golden Gate preparation enzyme" prop:value=move||format!("{:?}",config.with(|c|c.golden_gate.enzyme)) on:change=move|e|config.update(|c|c.golden_gate.enzyme=match event_target_value(&e).as_str(){"BsmBI"=>collagene_core::golden_gate::GoldenGateEnzyme::BsmBI,"BbsI"=>collagene_core::golden_gate::GoldenGateEnzyme::BbsI,"SapI"=>collagene_core::golden_gate::GoldenGateEnzyme::SapI,_=>collagene_core::golden_gate::GoldenGateEnzyme::BsaI})><option value="BsaI">"BsaI — 4-base fusion"</option><option value="BsmBI">"BsmBI — 4-base fusion"</option><option value="BbsI">"BbsI — 4-base fusion"</option><option value="SapI">"SapI — 3-base fusion"</option></select></label>
                <label>"Pre-recognition clamp"<input maxlength="12" aria-label="Golden Gate clamp" prop:value=move||config.with(|c|c.golden_gate.clamp.clone()) on:input=move|e|config.update(|c|c.golden_gate.clamp=event_target_value(&e))/></label>
                <label>"Cleavage spacer"<input maxlength="2" aria-label="Golden Gate spacer" prop:value=move||config.with(|c|c.golden_gate.spacer.clone()) on:input=move|e|config.update(|c|c.golden_gate.spacer=event_target_value(&e))/></label>
                <label>"Forward annealing length"<input type="number" min="1" max="200" aria-label="Golden Gate forward annealing length" prop:value=move||config.with(|c|c.golden_gate.forward_length.clone()) on:input=move|e|config.update(|c|c.golden_gate.forward_length=event_target_value(&e))/></label>
                <label>"Reverse annealing length"<input type="number" min="1" max="200" aria-label="Golden Gate reverse annealing length" prop:value=move||config.with(|c|c.golden_gate.reverse_length.clone()) on:input=move|e|config.update(|c|c.golden_gate.reverse_length=event_target_value(&e))/></label>
            </div>
            <p class="helper">"The clamp and spacer are explicit synthesized DNA. BbsI needs a 2-base spacer; the other listed enzymes need 1. These settings do not establish reagent compatibility, primer thermodynamics or amplification success."</p>
            </details>
            <button class="button secondary full" disabled=move||busy.get() on:click=move|_|suggest_gate.run(())>"Suggest Type IIS PCR primers"</button>
            {move||gate_suggestions.get().map(|saved|{
                let guarded=StoredValue::new(saved.clone());
                let edit=Callback::new(move|index:usize|guarded.with_value(|saved|{
                    if config.get_untracked()!=saved.config
                        ||source_context.get_untracked()!=Some(saved.context.clone())
                        ||target_info.get_untracked()!=saved.target_context||epoch.get_untracked()!=saved.epoch {
                        status.set("Sources or settings changed; request fresh Type IIS primers.".into());return;
                    }
                    focus_row.run(index);
                }));
                let roles=saved.config.role_ids(saved.config.rows.len()).unwrap_or_default();
                let bounded=saved.report.rows.len()<=12;
                let row_saved=saved.clone();
                view!{<section class="golden-gate-primer-proposals compact-primer-review"><h4>"Proposed PCR/digest chains — not applied"</h4>
                    <Show when=move||bounded fallback=||view!{<p role="alert">"Proposal exceeds the supported 12 preparation rows. No partial primer list is shown."</p>}>
                        <div class="primer-review-table-wrap"><table class="primer-review-table"><caption>"Type IIS primer proposals · expand a row to inspect complete oligos"</caption><thead><tr><th scope="col">"Fragment"</th><th scope="col">"Oligos"</th><th scope="col">"Actual digested ends"</th><th scope="col">"Review"</th></tr></thead><tbody>
                        {row_saved.report.rows.iter().map(|row|{
                            let index=roles.iter().position(|role|role.role==row.role&&role.input_id==row.input_id);
                            let source=index.and_then(|i|row_saved.context.get(i)).map(|source|source.name.clone()).unwrap_or_else(||row.role.clone());
                            view!{<tr><th scope="row"><button class="text-tool primer-review-fragment" disabled=index.is_none() on:click=move|_|if let Some(index)=index{edit.run(index)}>{format!("{} · {}",index.map(|i|format!("Fragment {}",i+1)).unwrap_or_else(||row.role.clone()),source)}</button><small>{format!("Digest fragment {}",row.fragment_index+1)}</small></th>
                                <td><span>{format!("F {} nt · R {} nt",row.forward.tail.len()+row.forward.annealing.len(),row.reverse.tail.len()+row.reverse.annealing.len())}</span><small>{format!("Annealing {} / {} nt",row.forward.annealing.len(),row.reverse.annealing.len())}</small></td>
                                <td><small>{format!("Left: {}",method_end(Some(&row.left)))}</small><small>{format!("Right: {}",method_end(Some(&row.right)))}</small></td>
                                <td><details class="primer-review-details"><summary>{format!("Inspect primers · {}",row.role)}</summary><div class="primer-review-domains">{proposed_oligo("Forward",&row.forward)}{proposed_oligo("Reverse",&row.reverse)}</div>
                                    <dl class="primer-domain-segments"><dt>"Forward tail · own 5′→3′"</dt><dd>{format!("clamp [{}] · motif [{}] · spacer [{}] · fusion [{}]",row.forward_tail.clamp,row.forward_tail.motif,row.forward_tail.spacer,row.forward_tail.fusion)}</dd><dt>"Reverse tail · own 5′→3′"</dt><dd>{format!("clamp [{}] · motif [{}] · spacer [{}] · fusion [{}]",row.reverse_tail.clamp,row.reverse_tail.motif,row.reverse_tail.spacer,row.reverse_tail.fusion)}</dd></dl>
                                    {crate::end_view::pair(Some(&row.left),Some(&row.right))}<p>{format!("Assembly-oriented top DNA: {}",short(&row.oriented_top))}</p>
                                </details></td></tr>}
                        }).collect_view()}
                        </tbody></table></div>
                    </Show>
                    {saved.report.warnings.iter().map(|warning|view!{<p class="primer-review-warning">{warning.clone()}</p>}).collect_view()}
                    <p class="helper">"Inspect complete oligos and cleaved ends before applying. Applying copies these named primers, the exact catalog pin and PCR-product fragment choice; preparation evidence must then be recomputed."</p>
                    <button class="button primary full" disabled=move||busy.get() on:click=move|_|apply_gate.run(())>"Apply Type IIS primer proposals"</button>
                </section>}
            })}
        </section></Show>
        {move||{let snapshot=reaction_context.get()?;
            let captured=StoredValue::new(snapshot.clone());
            let edit=Callback::new(move|row:usize|{match reaction_guard.run(captured.get_value()){Ok(())=>focus_row.run(row),Err(error)=>status.set(error)}});
            let view_source=on_source.map(|navigate|Callback::new(move|row:usize|{let result=captured.with_value(|ctx|{reaction_guard.run(ctx.clone())?;let input=ctx.config.rows.get(row).ok_or("Source role unavailable.")?;navigate.run((input.document,row))});if let Err(error)=result{status.set(error);}}));
            Some(view!{<PreparedPool proof=snapshot config guard=reaction_guard capture=reaction_capture epoch on_candidate on_finalize on_edit=edit on_source=view_source/>})
        }}
        <Show when=move||config.get().mode==Mode::Overlap&&!disable_plain_primer_design><section class="overlap-primer-planner"><h3>"2. Suggest and review PCR primers"</h3><p class="helper">"Use the selected source regions and order to suggest one primer pair per fragment, including added homology tails. Your current primers stay unchanged until you explicitly apply the proposal. The independent target is never generated or selected by this step."</p><details class="assembly-primer-settings"><summary>"Overlap primer design settings"</summary><div class="primer-grid"><label>"Junction homology (bp)"<input type="number" min="1" aria-label="Suggested overlap homology bases" prop:value=move||planner_form.with(|f|f.overlap.clone()) on:input=move|e|planner_form.update(|f|f.overlap=event_target_value(&e))/></label><label>"Minimum annealing length"<input type="number" min="1" max="200" aria-label="Suggested primer minimum length" prop:value=move||planner_form.with(|f|f.minimum.clone()) on:input=move|e|planner_form.update(|f|f.minimum=event_target_value(&e))/></label><label>"Maximum annealing length"<input type="number" min="1" max="200" aria-label="Suggested primer maximum length" prop:value=move||planner_form.with(|f|f.maximum.clone()) on:input=move|e|planner_form.update(|f|f.maximum=event_target_value(&e))/></label><label>"Target annealing Tm (°C)"<input type="number" step="0.1" aria-label="Suggested primer target Tm" prop:value=move||planner_form.with(|f|f.tm.clone()) on:input=move|e|planner_form.update(|f|f.tm=event_target_value(&e))/></label></div><p class="helper">"Default filters: GC 20–80%, pair ΔTm ≤5°C, homopolymer run ≤5 and 3′ complementary run ≤3. Tm assumes 50 mM sodium and 500 nM total strands; these are calculations, not reagent or experimental performance predictions."</p></details><button class="button secondary full" disabled=move||busy.get() on:click=move|_|suggest.run(())>"Suggest PCR primers"</button>
        {move||suggestions.get().map(|proposal|{
            let saved=proposal.clone();
            let guarded=StoredValue::new(proposal.clone());
            let edit=Callback::new(move|index:usize|guarded.with_value(|saved|{
                if config.get_untracked()!=saved.config||planner_form.get_untracked()!=saved.form
                    ||source_context.get_untracked()!=Some(saved.context.clone())
                    ||target_info.get_untracked()!=saved.target_context||target.get_untracked()!=saved.target
                    ||epoch.get_untracked()!=saved.epoch {
                    status.set("Sources or settings changed; request fresh suggestions.".into());return;
                }
                focus_row.run(index);
            }));
            let bounded=proposal.report.designs.len()<=12;
            let row_proposal=proposal.clone();
            view!{<div class="overlap-primer-suggestions compact-primer-review"><h4>"Proposed primer pairs — not yet applied"</h4>
                <Show when=move||bounded fallback=||view!{<p role="alert">"Proposal exceeds the supported 12 preparation rows. No partial primer list is shown."</p>}>
                    <div class="primer-review-table-wrap"><table class="primer-review-table"><caption>"Overlap primer proposals · expand a row to inspect complete oligos"</caption><thead><tr><th scope="col">"Fragment"</th><th scope="col">"Annealing Tm"</th><th scope="col">"PCR product"</th><th scope="col">"Review"</th></tr></thead><tbody>
                    {row_proposal.report.designs.iter().enumerate().map(|(i,design)|{
                        let interval=row_proposal.report.request.fragments.get(i);
                        let name=row_proposal.context.get(i).map(|source|source.name.as_str()).unwrap_or("source unavailable");
                        view!{<tr><th scope="row"><button class="text-tool primer-review-fragment" on:click=move|_|edit.run(i)>{format!("Fragment {} · {}",i+1,name)}</button>{interval.map(|region|view!{<small>{format!("{}–{} · {}",region.start+1,region.end,if region.reverse{"reverse"}else{"forward"})}</small>})}</th>
                            <td>{design.pairs.first().map(|pair|view!{<span>{format!("{} / {}",pair.forward.analysis.tm.as_ref().map(|tm|format!("{:.1}°C",tm.tm_c)).unwrap_or_else(||"unavailable".into()),pair.reverse.analysis.tm.as_ref().map(|tm|format!("{:.1}°C",tm.tm_c)).unwrap_or_else(||"unavailable".into()))}</span><small>{format!("Pair ΔTm {:.1}°C",pair.tm_difference_c)}</small>})}</td>
                            <td>{design.pairs.first().map(|pair|view!{<span>{format!("{} bp",pair.product_length)}</span><small>{format!("F {} nt · R {} nt",pair.forward.primer.tail.len()+pair.forward.primer.annealing.len(),pair.reverse.primer.tail.len()+pair.reverse.primer.annealing.len())}</small>})}</td>
                            <td>{design.pairs.is_empty().then(||view!{<span>"No proposed primer pair"</span>})}{design.pairs.first().map(|pair|view!{<details class="primer-review-details"><summary>{format!("Inspect primers · fragment {}",i+1)}</summary><div class="primer-review-domains">{proposed_oligo("Forward",&pair.forward.primer)}{proposed_oligo("Reverse",&pair.reverse.primer)}</div>{interval.map(|region|view!{<p>{format!("Source interval {}–{} (1-based inclusive, unrolled) · {} in assembly",region.start+1,region.end,if region.reverse{"reverse"}else{"forward"})}</p>})}</details>})}</td>
                        </tr>{(!design.warnings.is_empty()).then(||view!{<tr class="primer-review-warning-row"><td colspan="4">{design.warnings.iter().map(|warning|view!{<p class="primer-review-warning">{format!("Fragment {}: {}",i+1,warning)}</p>}).collect_view()}</td></tr>})}
                    }}).collect_view()}
                    </tbody></table></div>
                </Show>
                {proposal.report.warnings.iter().map(|warning|view!{<p class="primer-review-warning">{warning.clone()}</p>}).collect_view()}
                <p class="helper">"A deterministic proposal under these filters is not a unique or experimentally validated primer design. Applying changes the preparation to exact PCR; verify the resulting assembly next."</p>
                <button class="button primary full" disabled=move||busy.get() on:click=move|_|apply_suggestions.run(())>"Apply suggested primers"</button>
                <button class="button secondary full" on:click=move|_|{if config.get_untracked()!=saved.config||planner_form.get_untracked()!=saved.form||source_context.get_untracked()!=Some(saved.context.clone())||target_info.get_untracked()!=saved.target_context||target.get_untracked()!=saved.target||epoch.get_untracked()!=saved.epoch{status.set("Sources or settings changed; request fresh suggestions.".into());return;}match serde_json::to_vec_pretty(&serde_json::json!({"sources":saved.context,"target":saved.target_context,"report":saved.report})){Ok(bytes)=>if let Err(e)=download("overlap-primer-suggestions.json","application/json",&bytes){status.set(e)},Err(e)=>status.set(e.to_string())}}>"Download suggested primer evidence"</button>
            </div>}
        })}
        </section></Show><Show when=move||config.get().mode!=Mode::Tailing><h3>"3. Verify the prepared assembly"</h3><p class="helper">"Review the applied preparation and orientation, then compute method evidence against your independent target. Suggested primers alone do not verify a cloning design."</p></Show>
        <label class="checkbox-label"><input disabled=fixed_sources type="checkbox" aria-label="Method product circular" prop:checked=move||config.get().circular on:change=move|e|config.update(|c|c.circular=event_target_checked(&e))/>"Circular candidate product"</label><button class="button primary full" disabled=move||busy.get() on:click=move|_|compute.run(())>"Compute method evidence"</button><button class="text-tool" disabled=move||!busy.get() on:click=move|_|{generation.update(|n|*n=n.wrapping_add(1));if let Some(w)=active.get_value(){w.close();}active.set_value(None);busy.set(false);suggestions.set(None);
        gate_suggestions.set(None);status.set("Method analysis cancelled".into());}>"Cancel method analysis"</button><p class="method-status" role="status">{move||status.get()}</p>
        {move||evidence.get().map(|proof|{
            let report=proof.clone();
            let current=StoredValue::new(proof.clone());
            let edit=Callback::new(move|index:usize|{
                let outcome=current.with_value(|saved|{diagnostic_guard.run(saved.clone())?;
                    if index>=saved.config.rows.len(){return Err("The preparation row is unavailable.".into());}
                    focus_row.run(index);
                    Ok::<_,String>(())
                });if let Err(error)=outcome{status.set(error);}
            });
            let view_source=on_source.map(|navigate|Callback::new(move|row:usize|{
                let outcome=current.with_value(|saved|{diagnostic_guard.run(saved.clone())?;
                    let input=saved.config.rows.get(row).ok_or("The original source row is unavailable.")?;
                    navigate.run((input.document,row))
                });if let Err(error)=outcome{status.set(error);}
            }));
            let diagnostics=proof.overlap.as_ref().map(collagene_core::assembly_diagnostics::overlap).or_else(||proof.ligation.as_ref().map(collagene_core::assembly_diagnostics::ligation));
            let names=proof.context.iter().map(|source|source.name.clone()).collect::<Vec<_>>();
            let preparations=proof.config.rows.iter().map(OverlapInput::summary).collect::<Vec<_>>();
            view!{<div class="method-evidence">
            {move||current.with_value(|proof|{
                let physical=proof.ligation.as_ref().map(|p|&p.physical).or_else(||proof.tailing.as_ref().map(|p|&p.physical))?;
                let index=active_row.map(|row|row.get()).unwrap_or(0);
                let label=proof.context.get(index).map(|source|source.name.clone()).unwrap_or_else(||format!("Fragment {}",index+1));
                Some(view!{<section class="active-fragment-end-preview"><h3>{format!("Fragment {} · {} · reconstructed ends",index+1,label)}</h3>
                    {physical.get(index).and_then(Option::as_ref).map(crate::end_view::summary).unwrap_or_else(||view!{<p class="physical-end-unavailable">"Physical ends unavailable. Review the blocking preparation result below."</p>}.into_any())}
                    <p class="helper">"Shown in assembly orientation after declared preparation and terminal treatment."</p>
                </section>})
            })}
            {diagnostics.map(|result|match result{Ok(diagnostics)=>view!{<MethodDiagnostics diagnostics names preparations on_edit=edit on_source=view_source/>}.into_any(),Err(error)=>view!{<p class="helper">{format!("Structured preparation diagnostics unavailable: {error}. The full fresh findings remain below.")}</p>}.into_any()})}
            {proof.ligation.as_ref().map(|p|(&p.compiled.plan,&p.verification)).or_else(||proof.overlap.as_ref().map(|p|(&p.compiled.plan,&p.verification))).or_else(||proof.tailing.as_ref().map(|p|(&p.compiled.plan,&p.verification))).map(|(plan,report)|{
                let captured = StoredValue::new(proof.clone());
                let captured_sources = StoredValue::new(Some(proof.context.clone()));
                let captured_target = StoredValue::new(proof.target_context.clone());
                let guarded_finalize = on_finalize.map(|finalize|Callback::new(move|candidate:Candidate| {
                    let fresh = captured.with_value(|saved| {
                        evidence.with_untracked(|current|current.as_ref().is_some_and(|current|Arc::ptr_eq(current,saved)))
                            && config.get_untracked()==saved.config
                            && source_revisions.get_untracked()==Some(saved.revisions.clone())
                            && target.get_untracked()==saved.target
                            && epoch.get_untracked()==saved.epoch
                    });
                    if !fresh || source_context.get_untracked()!=captured_sources.get_value() || target_info.get_untracked()!=captured_target.get_value() {
                        return Err("Sources or preparation settings changed. Compute fresh method evidence before creating outputs.".into());
                    }
                    finalize.run(candidate)
                }));
                view!{<crate::assembly_product_view::AssemblyProductView plan=plan.clone() target=molecule.get_untracked() report=report.clone() epoch on_products on_finalize=guarded_finalize/>}
            })}

            {proof.tailing.as_ref().map(|preview|view!{<article class="tail-evidence"><h3>"Reconstructed tailing and ligation"</h3><p>{format!("Preparation reconstructed: {} · independent-target assessment: {}",preview.preparation_complete,method_status(&preview.verification.status))}</p><p class="helper">"Ends are reconstructed from the pinned originals, then tailed and treated as declared. Junction bases are compared with your separately reviewed intended DNA."</p>{preview.junctions.iter().map(|junction|view!{<details><summary>{format!("Junction {} · fragment {} → {}{}",junction.index+1,junction.upstream+1,junction.downstream+1,if junction.closing{" · circle closure"}else{""})}</summary><div class="physical-junction-ends"><section><h4>"Upstream right end"</h4>{crate::end_view::right(junction.upstream_end.as_ref())}</section><section><h4>"Downstream left end"</h4>{crate::end_view::left(junction.downstream_end.as_ref())}</section></div></details>}).collect_view()}{preview.verification.findings.iter().map(|finding|view!{<p>{finding.clone()}</p>}).collect_view()}</article>})}
            {proof.tails.iter().enumerate().map(|(i,preview)|{let opened=preview.molecule.clone();view!{<article class="tail-evidence"><h3>{format!("{} · {}-tailed duplex",if i==0{"Vector"}else{"Insert"},preview.base.as_str())}</h3><pre>{format!("Top 5′→3′: {}\nBottom 5′→3′: {}\n{}",short(preview.tailing.duplex.top.as_str()),short(preview.tailing.duplex.bottom.as_str()),preview.tailing.duplex.summary().map(|s|s.describe()).unwrap_or_else(|e|e.to_string()))}</pre><details><summary>"Source, additions and assumptions"</summary><pre>{format!("Source top 5′→3′: {}\nSource bottom 5′→3′: {}",short(preview.tailing.source.top.as_str()),short(preview.tailing.source.bottom.as_str()))}</pre>{preview.tailing.additions.iter().map(|a|view!{<p>{format!("New {:?} residue {}–{}: {} (no template ancestry)",a.strand,a.product_start+1,a.product_end,a.sequence)}</p>}).collect_view()}{preview.findings.iter().chain(&preview.tailing.warnings).map(|w|view!{<p>{w.clone()}</p>}).collect_view()}</details><button class="button secondary full" on:click=move|_|match on_products.run(vec![opened.clone()]){Ok(())=>status.set("Reference preview opened; reconstruct its original graph for later physical operations.".into()),Err(e)=>status.set(e)}>{format!("Open {} reference preview",if i==0{"vector"}else{"insert"})}</button></article>}}).collect_view()}
            {proof.ligation.as_ref().map(|preview|view!{<article class="ordered-ligation-evidence"><h3>"Reconstructed ligation evidence"</h3><p>{format!("Preparation reconstructed: {} · independent-target assessment: {}",preview.preparation_complete,method_status(&preview.verification.status))}</p><p class="helper">"These ends and junctions were freshly reconstructed from the original source operations. An end-compatible ordered candidate is not a unique mixture product or demonstrated reaction outcome."</p>{preview.prepared.iter().enumerate().map(|(i,fragment)|view!{<details><summary>{format!("Fragment {} · {} bp · {}",i+1,fragment.sequence.len(),proof.config.rows.get(i).map(OverlapInput::summary).unwrap_or_default())}</summary><pre>{format!("Prepared source before orientation: {}",short(fragment.sequence.as_str()))}</pre>{preview.physical.get(i).and_then(Option::as_ref).map(crate::end_view::summary)}<p>{preview.physical.get(i).and_then(Option::as_ref).map(method_duplex).unwrap_or_else(||"Physical end geometry unavailable; no chemistry inferred".into())}</p><p class="helper">"End summary is in assembly orientation, after any explicit terminal treatment."</p></details>}).collect_view()}{preview.junctions.iter().map(|junction|view!{<details><summary>{format!("Junction {} · fragment {} → {}{}",junction.index+1,junction.upstream+1,junction.downstream+1,if junction.closing{" · circle closure"}else{""})}</summary><div class="physical-junction-ends"><section><h4>"Upstream fragment: right end"</h4>{crate::end_view::right(junction.upstream_end.as_ref())}</section><section><h4>"Downstream fragment: left end"</h4>{crate::end_view::left(junction.downstream_end.as_ref())}</section></div><p>{format!("Upstream right: {}",method_end(junction.upstream_end.as_ref()))}</p><p>{format!("Downstream left: {}",method_end(junction.downstream_end.as_ref()))}</p></details>}).collect_view()}{preview.product_physical.as_ref().map(|physical|view!{<section><h4>"Reconstructed product"</h4><p>{method_duplex(physical)}</p>{physical.nicks.iter().map(|nick|view!{<p>{format!("Retained {} strand nick at top-reference boundary {}",match nick.strand{collagene_core::chemistry::Strand::Top=>"top",collagene_core::chemistry::Strand::Bottom=>"bottom"},nick.position)}</p>}).collect_view()}</section>})}{preview.verification.findings.iter().map(|finding|view!{<p>{finding.clone()}</p>}).collect_view()}</article>})}
            {proof.overlap.as_ref().map(|preview|view!{<article class="overlap-preparation-evidence"><h3>"Prepared fragments from original sources"</h3><p>{format!("Preparation reconstructed: {} · independent-target assessment: {}",preview.preparation_complete,method_status(&preview.verification.status))}</p><p class="helper">"Prepared reference sequences below are reconstructed from the declared originals. They are not substituted as plan inputs, confirmed amplicons or physical reaction outcomes."</p>{preview.prepared.iter().enumerate().map(|(i,fragment)|view!{<details><summary>{format!("Fragment {} · {} bp · {}",i+1,fragment.sequence.len(),proof.config.rows.get(i).map(OverlapInput::summary).unwrap_or_default())}</summary><pre>{format!("Prepared source 5′→3′ (before orientation): {}",short(fragment.sequence.as_str()))}</pre></details>}).collect_view()}{preview.verification.findings.iter().map(|finding|view!{<p>{finding.clone()}</p>}).collect_view()}</article>})}
            {proof.profile.as_ref().map(|p|view!{<article class="overlap-profile-evidence"><h3>"Overlap evidence"</h3><p>{format!("Geometry consistent: {} · analysis complete: {} · published length guidance: {}",p.geometry_valid,p.complete,match p.recommendation_compliant{Some(true)=>"met",Some(false)=>"outside published guidance",None=>"not established"})}</p><p class="helper">{format!("{} comparisons; no exhaustive alternative-product search",p.comparisons)}</p>{p.junctions.iter().map(|j|view!{<details class="profile-junction"><summary>{format!("Junction {} · {}bp · {}",j.index+1,j.overlap_bases,if j.exact_match{"exact end homology"}else{"incompatible homology"})}</summary><pre>{format!("Upstream {}\nDownstream {}",j.homology,j.downstream_homology)}</pre><p>{format!("GC: {} · advisory length: {} · repeat search complete: {}",j.gc_fraction.map(|f|format!("{:.1}%",f*100.)).unwrap_or_else(||"unavailable".into()),guidance(j.recommended_min,j.recommended_max),j.repeat_search_complete)}</p>{j.occurrences.iter().take(40).map(|o|view!{<p>{format!("Source {} · {}–{}{}{}",o.fragment+1,o.start+1,o.end,if o.reverse{" reverse"}else{" forward"},if o.assigned_junction{" · assigned end"}else{" · additional occurrence"})}</p>}).collect_view()}<p class="helper">"Up to 40 occurrence positions shown; full evidence download retains every reported position."</p></details>}).collect_view()}<details><summary>"Profile findings and assumptions"</summary>{p.findings.iter().chain(&p.warnings).map(|w|view!{<p>{w.clone()}</p>}).collect_view()}</details></article>})}
            <button class="button primary full" on:click=move|_|review.run(())>"Review method candidate"</button><button class="button secondary full" on:click=move|_|{let result=if let Some(preview)=&report.tailing{serde_json::to_vec_pretty(preview)}else if let Some(preview)=&report.ligation{serde_json::to_vec_pretty(preview)}else if let Some(preview)=&report.overlap{serde_json::to_vec_pretty(preview)}else if let Some(profile)=&report.profile{serde_json::to_vec_pretty(profile)}else{serde_json::to_vec_pretty(&report.tails)};match result{Ok(bytes)=>if let Err(e)=download("cloning-method-evidence.json","application/json",&bytes){status.set(e)},Err(e)=>status.set(e.to_string())}}>"Download method evidence"</button>
        </div>}})}
    </div>}
}

#[component]
fn TailFields(
    vector: bool,
    config: RwSignal<MethodsState>,
    query_options: Memo<Vec<(usize, String)>>,
    primer_sources: Memo<Vec<crate::digest::SourceInfo>>,
    list_primers: Callback<
        collagene_core::run_lineage::EventBinding,
        Result<collagene_core::primer_reuse::List, String>,
    >,
    resolve_primer: Callback<
        crate::primer_picker::ResolveRequest,
        Result<collagene_core::primer_reuse::Resolved, String>,
    >,
    epoch: RwSignal<u64>,
    fixed_sources: bool,
) -> impl IntoView {
    let label = if fixed_sources {
        if vector { "Fragment 1" } else { "Fragment 2" }
    } else if vector {
        "Vector"
    } else {
        "Insert"
    };
    view! {<fieldset class="tail-source-fields"><legend>{label}</legend><Show when=move||!fixed_sources><label>"Source document"<select aria-label=format!("{label} tailing source") prop:value=move||config.with(|c|c.input(vector).document.to_string()) on:change=move|e|config.update(|c|c.input_mut(vector).document=event_target_value(&e).parse().unwrap_or(0))>{move||query_options.get().into_iter().map(|(i,name)|view!{<option value=i.to_string()>{name}</option>}).collect_view()}</select></label><label>"Preparation"<select aria-label=format!("{label} tailing preparation") prop:value=move||if config.with(|c|c.input(vector).pcr){"pcr"}else{"digest"} on:change=move|e|config.update(|c|c.input_mut(vector).pcr=event_target_value(&e)=="pcr")><option value="digest">"Restriction fragment"</option><option value="pcr">"Unique exact PCR product"</option></select></label>
        <Show when=move||!config.with(|c|c.input(vector).pcr)><label>"Restriction enzymes"<input aria-label=format!("{label} tailing enzymes") prop:value=move||config.with(|c|c.input(vector).enzymes.clone()) on:input=move|e|config.update(|c|c.input_mut(vector).enzymes=event_target_value(&e))/></label><label>"Digest fragment (cut order, 1-based)"<input type="number" min="1" aria-label=format!("{label} tailing fragment") prop:value=move||config.with(|c|c.input(vector).fragment.clone()) on:input=move|e|config.update(|c|c.input_mut(vector).fragment=event_target_value(&e))/></label><label class="checkbox-label"><input type="checkbox" aria-label=format!("{label} blunt before tailing") prop:checked=move||config.with(|c|c.input(vector).blunt) on:change=move|e|config.update(|c|c.input_mut(vector).blunt=event_target_checked(&e))/>"Blunt staggered ends before tailing"</label></Show>
        <Show when=move||config.with(|c|c.input(vector).pcr)>
        <crate::primer_picker::PrimerPicker label=format!("{label} tailing forward") sources=primer_sources list_primers resolve_primer epoch on_copy=Callback::new(move|primer:Primer|config.update(|c|{let r=c.input_mut(vector);r.forward_name=primer.name;r.forward=primer.annealing.as_str().into();r.forward_tail=primer.tail;}))/>
        <label>"Forward primer name"<input maxlength="256" aria-label=format!("{label} tailing forward primer name") prop:value=move||config.with(|c|c.input(vector).forward_name.clone()) on:input=move|e|config.update(|c|c.input_mut(vector).forward_name=event_target_value(&e))/></label>
        <crate::primer_picker::PrimerPicker label=format!("{label} tailing reverse") sources=primer_sources list_primers resolve_primer epoch on_copy=Callback::new(move|primer:Primer|config.update(|c|{let r=c.input_mut(vector);r.reverse_name=primer.name;r.reverse=primer.annealing.as_str().into();r.reverse_tail=primer.tail;}))/>
        <label>"Reverse primer name"<input maxlength="256" aria-label=format!("{label} tailing reverse primer name") prop:value=move||config.with(|c|c.input(vector).reverse_name.clone()) on:input=move|e|config.update(|c|c.input_mut(vector).reverse_name=event_target_value(&e))/></label>
    <label>"Forward annealing sequence 5′→3′"<input aria-label=format!("{label} tailing forward primer") prop:value=move||config.with(|c|c.input(vector).forward.clone()) on:input=move|e|config.update(|c|c.input_mut(vector).forward=event_target_value(&e))/></label><label>"Reverse annealing sequence 5′→3′"<input aria-label=format!("{label} tailing reverse primer") prop:value=move||config.with(|c|c.input(vector).reverse.clone()) on:input=move|e|config.update(|c|c.input_mut(vector).reverse=event_target_value(&e))/></label><label>"Forward primer 5′ extension"<input aria-label=format!("{label} tailing forward extension") prop:value=move||config.with(|c|c.input(vector).forward_tail.clone()) on:input=move|e|config.update(|c|c.input_mut(vector).forward_tail=event_target_value(&e))/></label><label>"Reverse primer 5′ extension"<input aria-label=format!("{label} tailing reverse extension") prop:value=move||config.with(|c|c.input(vector).reverse_tail.clone()) on:input=move|e|config.update(|c|c.input_mut(vector).reverse_tail=event_target_value(&e))/></label><details><summary>"Review current primer domains"</summary><pre>{move||config.with(|c|{let r=c.input(vector);format!("{}\nOrdered 5′→3′: {}{}\nAnnealing: {}\n5′ tail: {}\n\n{}\nOrdered 5′→3′: {}{}\nAnnealing: {}\n5′ tail: {}",r.forward_name,r.forward_tail,r.forward,r.forward,r.forward_tail,r.reverse_name,r.reverse_tail,r.reverse,r.reverse,r.reverse_tail)})}</pre></details></Show>
        </Show><label>"One nucleotide added to BOTH 3′ ends"<select aria-label=format!("{label} tailing base") prop:value=move||config.with(|c|c.input(vector).base.as_str()) on:change=move|e|config.update(|c|c.input_mut(vector).base=match event_target_value(&e).as_str(){"C"=>TailingBase::C,"G"=>TailingBase::G,"T"=>TailingBase::T,_=>TailingBase::A})><option value="A">"A"</option><option value="C">"C"</option><option value="G">"G"</option><option value="T">"T"</option></select></label><label>"Explicit terminal treatment after tailing"<select aria-label=format!("{label} tailing terminal treatment") prop:value=move||config.with(|c|match c.input(vector).treatment{Treatment::Keep=>"keep",Treatment::Phosphorylate=>"phosphorylate",Treatment::Dephosphorylate=>"dephosphorylate"}) on:change=move|e|config.update(|c|c.input_mut(vector).treatment=match event_target_value(&e).as_str(){"phosphorylate"=>Treatment::Phosphorylate,"dephosphorylate"=>Treatment::Dephosphorylate,_=>Treatment::Keep})><option value="keep">"Keep reconstructed 5′ chemistry"</option><option value="phosphorylate">"Add 5′ phosphorylation"</option><option value="dephosphorylate">"Add 5′ dephosphorylation"</option></select></label><Show when=move||!fixed_sources><label class="checkbox-label"><input type="checkbox" aria-label=format!("Reverse {label} tailing product") prop:checked=move||config.with(|c|c.input(vector).reversed) on:change=move|e|config.update(|c|c.input_mut(vector).reversed=event_target_checked(&e))/>"Reverse fragment for ligation"</label></Show></fieldset>}
}

#[component]
fn OverlapFields(
    config: RwSignal<MethodsState>,
    index: usize,
    query_options: Memo<Vec<(usize, String)>>,
    primer_sources: Memo<Vec<crate::digest::SourceInfo>>,
    list_primers: Callback<
        collagene_core::run_lineage::EventBinding,
        Result<collagene_core::primer_reuse::List, String>,
    >,
    resolve_primer: Callback<
        crate::primer_picker::ResolveRequest,
        Result<collagene_core::primer_reuse::Resolved, String>,
    >,
    epoch: RwSignal<u64>,
    selected_region: Callback<(), Option<(usize, usize, usize)>>,
    fixed_sources: bool,
    edited_junctions: bool,
    get_query: Callback<usize, Option<Molecule>>,
    source_info: Callback<usize, Option<crate::digest::SourceInfo>>,
    get_workspace: Option<Callback<(), Result<collagene_core::workspace::Workspace, String>>>,
) -> impl IntoView {
    let row = move || config.with(|c| c.rows.get(index).cloned().unwrap_or_default());
    let region_status = RwSignal::new(String::new());

    view! {<fieldset class="overlap-source-row"><legend>{format!("Fragment {}",index+1)}</legend>
        <Show when=move||fixed_sources><p class="preparation-source-summary">{move||{
            let r=row();
            let name=query_options.with(|options|options.iter().find(|(i,_)|*i==r.document).map(|(_,name)|name.clone()).unwrap_or_else(||"Source unavailable".into()));
            format!("{} · source {}–{} (1-based, unrolled) · {}",name,r.interval_start,if r.interval_end.is_empty(){"end"}else{&r.interval_end},if r.reverse{"reverse"}else{"forward"})
        }}</p></Show>
        <Show when=move||!fixed_sources>
        <label>"Original source"<select disabled=fixed_sources aria-label=move||format!("{} source {}",preparation_label(config.get().mode),index+1) prop:value=move||row().document.to_string() on:change=move|e|config.update(|c|if let Some(r)=c.rows.get_mut(index){r.document=event_target_value(&e).parse().unwrap_or(0)})>{move||query_options.get().into_iter().map(|(n,label)|view!{<option value=n.to_string()>{label}</option>}).collect_view()}</select></label>
        <Show when=move||matches!(config.get().mode,Mode::Overlap|Mode::GoldenGate)><button class="button secondary" disabled=fixed_sources aria-label=format!("Set fragment {} from sequence selection",index+1) on:click=move|_|{match selected_region.run(()){Some((document,start,end)) if start<end&&start.checked_add(1).is_some()=>{config.update(|c|if let Some(row)=c.rows.get_mut(index){row.document=document;row.interval_start=(start+1).to_string();row.interval_end=end.to_string();});region_status.set("Selected source region copied. Review preparation and orientation; primers have not been changed.".into());},_=>region_status.set("Select one contiguous region in the sequence view first. Compound or origin-spanning feature selections need an explicit interval.".into())}} >"Set from sequence selection"</button><label>"Primer design interval start (1-based inclusive)"<input disabled=fixed_sources type="number" min="1" aria-label=format!("Overlap interval start {}",index+1) prop:value=move||row().interval_start on:input=move|e|config.update(|c|if let Some(r)=c.rows.get_mut(index){r.interval_start=event_target_value(&e)})/></label><label>"Primer design interval end (1-based inclusive; blank = source end)"<input disabled=fixed_sources type="number" min="1" aria-label=format!("Overlap interval end {}",index+1) prop:value=move||row().interval_end on:input=move|e|config.update(|c|if let Some(r)=c.rows.get_mut(index){r.interval_end=event_target_value(&e)})/></label><p class="helper">"Intervals declare intended regions for primer proposals. For a circular source, an end beyond its length describes one unrolled origin crossing. Existing whole-source, digest and manual PCR preparations are not silently sliced."</p></Show>
        </Show>
        <p class="helper preparation-region-status" role="status">{move||region_status.get()}</p>
        <Show when=move||config.get().mode==Mode::GoldenGate&&!edited_junctions><label>"Right fusion end (assembly direction, 5′→3′)"<input maxlength="4" aria-label=format!("Golden Gate right fusion {}",index+1) prop:value=move||row().right_fusion on:input=move|e|config.update(|c|if let Some(row)=c.rows.get_mut(index){row.right_fusion=event_target_value(&e)})/></label><p class="helper">{move||format!("Enter {} concrete bases explicitly. The left fusion comes from this oriented source region; enter the next fragment’s first fusion bases explicitly for a matching junction. The right fusion is added design DNA. Existing-site compatibility is checked from its actual ends.",if config.with(|c|c.golden_gate.enzyme)==collagene_core::golden_gate::GoldenGateEnzyme::SapI{3}else{4})}</p></Show>
        <label>"Prepare fragment"<select aria-label=move||format!("{} preparation {}",preparation_label(config.get().mode),index+1) prop:value=move||match row().preparation{OverlapPreparation::WholeLinear=>"whole",OverlapPreparation::DeclaredDuplex=>"declared",OverlapPreparation::Digest=>"digest",OverlapPreparation::CatalogDigest=>"catalog",OverlapPreparation::ExactPcr=>"pcr",OverlapPreparation::MismatchPcr=>"mismatch_pcr",OverlapPreparation::PcrCatalogDigest=>"pcr_catalog",OverlapPreparation::AnnealedOligos=>"anneal"} on:change=move|e|config.update(|c|if let Some(r)=c.rows.get_mut(index){r.preparation=match event_target_value(&e).as_str(){"declared"=>OverlapPreparation::DeclaredDuplex,"digest"=>OverlapPreparation::Digest,"catalog"=>OverlapPreparation::CatalogDigest,"pcr"=>OverlapPreparation::ExactPcr,"mismatch_pcr"=>OverlapPreparation::MismatchPcr,"pcr_catalog"=>OverlapPreparation::PcrCatalogDigest,"anneal"=>OverlapPreparation::AnnealedOligos,_=>OverlapPreparation::WholeLinear}})><option value="anneal" disabled=move||row().run_origin.is_some()||config.get().mode!=Mode::Ligation||(fixed_sources&&row().anneal_bottom_document.is_none())>"Annealed oligo pair"</option><option value="whole">"Whole linear DNA"</option><option value="declared" disabled=move||row().run_origin.is_some()>"Declare intact linear duplex"</option><option value="digest">"Manual starter-enzyme digest"</option><option value="catalog">"Pinned catalog digest fragment"</option><option value="pcr">"Exact PCR product"</option><option value="mismatch_pcr">"Intentional-mismatch PCR product"</option><option value="pcr_catalog">"Exact PCR → pinned catalog digest"</option></select></label>
        {get_workspace.map(|snapshot|view!{<Show when=move||row().run_origin.is_some()||!matches!(row().preparation,OverlapPreparation::DeclaredDuplex|OverlapPreparation::AnnealedOligos)><crate::run_origin_view::RunOriginPicker config index source_info epoch get_workspace=snapshot/></Show>})}
        <Show when=move||row().preparation==OverlapPreparation::WholeLinear><p class="helper">"Use an already linear source. Circular DNA must first be prepared by restriction digest or exact PCR; its origin is not an assumed cut."</p></Show>
        <Show when=move||row().preparation==OverlapPreparation::AnnealedOligos><crate::anneal_preparation_view::AnnealPreparationFields config index query_options get_query source_info epoch fixed_sources/></Show>
        <Show when=move||row().preparation==OverlapPreparation::DeclaredDuplex><p class="helper">"Operator assumption for an original, intact linear duplex: declare its two terminal groups explicitly. This does not measure chemistry, cut or repair DNA, change its reference sequence, or reuse producer metadata. Circular, single-stranded and imported physical/reference-only intermediate sources are rejected by the core."</p><DeclaredEndFields config index left=true/><DeclaredEndFields config index left=false/><p class="helper">"All groups initially remain unknown. Orientation is applied after this declaration; any ligation terminal treatment remains a separate explicit step."</p></Show>
        <Show when=move||matches!(row().preparation,OverlapPreparation::Digest|OverlapPreparation::CatalogDigest|OverlapPreparation::PcrCatalogDigest)><Show when=move||row().preparation==OverlapPreparation::Digest>
            <label>"Restriction enzymes (comma-separated)"<input aria-label=move||format!("{} digest enzymes {}",preparation_label(config.get().mode),index+1) prop:value=move||row().enzymes on:input=move|e|config.update(|c|if let Some(r)=c.rows.get_mut(index){r.enzymes=event_target_value(&e)})/></label>
            <label>"Digest fragment number (1-based)"<input type="number" min="1" aria-label=move||format!("{} digest fragment {}",preparation_label(config.get().mode),index+1) prop:value=move||row().fragment on:input=move|e|config.update(|c|if let Some(r)=c.rows.get_mut(index){r.fragment=event_target_value(&e)})/></label>
            </Show><Show when=move||row().preparation==OverlapPreparation::PcrCatalogDigest><crate::pcr_digest_view::PcrDigestPicker config index get_query source_info epoch/></Show><Show when=move||matches!(row().preparation,OverlapPreparation::CatalogDigest|OverlapPreparation::PcrCatalogDigest)>{move||row().catalog_selection.map(|selection|view!{<p>{format!("Selected fragment {} · {}",row().fragment,selection.enzymes.iter().map(|e|e.enzyme.name.as_str()).collect::<Vec<_>>().join(", "))}</p><details><summary>"Pinned catalog identity"</summary><p>{format!("Edition {} · content {}",selection.catalog_edition,selection.catalog_content_hash)}</p>{selection.enzymes.iter().map(|entry|view!{<p>{format!("{} · {} · {} · cuts {}/{}",entry.id,entry.enzyme.name,entry.enzyme.recognition,entry.enzyme.top_cut,entry.enzyme.bottom_cut)}</p>}).collect_view()}</details>})}<p class="helper">"The exact catalog selection and fragment index are retained without starter-enzyme name substitution. For an original-source fragment, choose it in Virtual digest. PCR-product fragments use the explicit PCR/digest settings here. Saved pins describe a request; they do not restore current digest status or prove cleavage."</p></Show>
            <label>"Optional explicit end repair"<select aria-label=move||format!("{} end repair {}",preparation_label(config.get().mode),index+1) prop:value=move||match row().repair{None=>"none",Some(EndRepairMethod::FillFivePrimeOverhangs)=>"fill",Some(EndRepairMethod::RemoveThreePrimeOverhangs)=>"remove",Some(EndRepairMethod::Blunt)=>"blunt"} on:change=move|e|config.update(|c|if let Some(r)=c.rows.get_mut(index){r.repair=match event_target_value(&e).as_str(){"fill"=>Some(EndRepairMethod::FillFivePrimeOverhangs),"remove"=>Some(EndRepairMethod::RemoveThreePrimeOverhangs),"blunt"=>Some(EndRepairMethod::Blunt),_=>None}})><option value="none">"None — retain cohesive ends"</option><option value="fill">"Fill 5′ overhangs"</option><option value="remove">"Remove 3′ overhangs"</option><option value="blunt">"Blunt ends"</option></select></label>
            <p class="helper">"Fragments are reconstructed in core digest order using the stated complete-cleavage assumption. This is not evidence of enzyme activity, methylation permission or reagent compatibility."</p>
        </Show>
        <Show when=move||row().is_pcr()>
            <crate::primer_picker::PrimerPicker label=format!("Fragment {} forward",index+1) sources=primer_sources list_primers resolve_primer epoch on_copy=Callback::new(move|primer:Primer|config.update(|c|if let Some(r)=c.rows.get_mut(index){r.forward_name=primer.name;r.forward=primer.annealing.as_str().into();r.forward_tail=primer.tail;}))/>
            <label>"Forward primer name"<input maxlength="256" aria-label=move||format!("{} forward primer name {}",preparation_label(config.get().mode),index+1) prop:value=move||row().forward_name on:input=move|e|config.update(|c|if let Some(r)=c.rows.get_mut(index){r.forward_name=event_target_value(&e)})/></label>
            <crate::primer_picker::PrimerPicker label=format!("Fragment {} reverse",index+1) sources=primer_sources list_primers resolve_primer epoch on_copy=Callback::new(move|primer:Primer|config.update(|c|if let Some(r)=c.rows.get_mut(index){r.reverse_name=primer.name;r.reverse_primer=primer.annealing.as_str().into();r.reverse_tail=primer.tail;}))/>
            <label>"Reverse primer name"<input maxlength="256" aria-label=move||format!("{} reverse primer name {}",preparation_label(config.get().mode),index+1) prop:value=move||row().reverse_name on:input=move|e|config.update(|c|if let Some(r)=c.rows.get_mut(index){r.reverse_name=event_target_value(&e)})/></label>

            <label>"Forward annealing sequence 5′→3′"<input aria-label=move||format!("{} forward primer {}",preparation_label(config.get().mode),index+1) prop:value=move||row().forward on:input=move|e|config.update(|c|if let Some(r)=c.rows.get_mut(index){r.forward=event_target_value(&e)})/></label>
            <label>"Forward primer 5′ extension"<input aria-label=move||format!("{} forward extension {}",preparation_label(config.get().mode),index+1) prop:value=move||row().forward_tail on:input=move|e|config.update(|c|if let Some(r)=c.rows.get_mut(index){r.forward_tail=event_target_value(&e)})/></label>
            <label>"Reverse annealing sequence 5′→3′"<input aria-label=move||format!("{} reverse primer {}",preparation_label(config.get().mode),index+1) prop:value=move||row().reverse_primer on:input=move|e|config.update(|c|if let Some(r)=c.rows.get_mut(index){r.reverse_primer=event_target_value(&e)})/></label>
            <label>"Reverse primer 5′ extension"<input aria-label=move||format!("{} reverse extension {}",preparation_label(config.get().mode),index+1) prop:value=move||row().reverse_tail on:input=move|e|config.update(|c|if let Some(r)=c.rows.get_mut(index){r.reverse_tail=event_target_value(&e)})/></label>
            <Show when=move||row().preparation==OverlapPreparation::MismatchPcr>
                <div class="primer-settings-grid">
                    <label>"Maximum mismatches per annealing domain"<input type="number" min="0" max="200" aria-label=format!("PCR maximum mismatches {}",index+1) prop:value=move||row().mismatch_max_mismatches on:input=move|e|config.update(|c|if let Some(r)=c.rows.get_mut(index){r.mismatch_max_mismatches=event_target_value(&e)})/></label>
                    <label>"Required exact 3′ bases"<input type="number" min="0" max="200" aria-label=format!("PCR exact three prime bases {}",index+1) prop:value=move||row().mismatch_three_prime_exact on:input=move|e|config.update(|c|if let Some(r)=c.rows.get_mut(index){r.mismatch_three_prime_exact=event_target_value(&e)})/></label>
                </div>
                <p class="helper">"Primer letters specify the introduced changes. Verification compares the reconstructed product with your independently declared target; primer settings never change that target. Internal changes must already be present in the intended reference. Binding compatibility does not establish amplification."</p>
            </Show>
            <details><summary>"Review current primer domains"</summary><pre>{move||{let r=row();format!("{}\nOrdered 5′→3′: {}{}\nAnnealing: {}\n5′ tail: {}\n\n{}\nOrdered 5′→3′: {}{}\nAnnealing: {}\n5′ tail: {}",r.forward_name,r.forward_tail,r.forward,r.forward,r.forward_tail,r.reverse_name,r.reverse_tail,r.reverse_primer,r.reverse_primer,r.reverse_tail)}}</pre></details>
            <p class="helper">"Enter annealing domains and added tails separately. Both oligos are written 5′→3′; no primer design or unique amplicon is assumed."</p>
        </Show>
        <Show when=move||!fixed_sources><label class="checkbox-label"><input disabled=fixed_sources type="checkbox" aria-label=format!("Reverse overlap source {}",index+1) prop:checked=move||row().reverse on:change=move|e|config.update(|c|if let Some(r)=c.rows.get_mut(index){r.reverse=event_target_checked(&e)})/>"Reverse complement prepared fragment"</label></Show>
        <Show when=move||matches!(config.get().mode,Mode::Ligation|Mode::GoldenGate)><label>"Explicit 5′ terminal treatment"<select aria-label=format!("Ligation terminal treatment {}",index+1) prop:value=move||match row().treatment{Treatment::Keep=>"none",Treatment::Phosphorylate=>"phosphorylate",Treatment::Dephosphorylate=>"dephosphorylate"} on:change=move|e|config.update(|c|if let Some(r)=c.rows.get_mut(index){r.treatment=match event_target_value(&e).as_str(){"phosphorylate"=>Treatment::Phosphorylate,"dephosphorylate"=>Treatment::Dephosphorylate,_=>Treatment::Keep}})><option value="none">"None — keep reconstructed chemistry"</option><option value="phosphorylate">"Add 5′ phosphorylation"</option><option value="dephosphorylate">"Add 5′ dephosphorylation"</option></select></label><p class="helper">"Treatment does not resolve unknown 3′ chemistry, repair ends or imply host nick repair."</p></Show>
        <p class="helper">{move||row().summary()}</p>
        <Show when=move||!fixed_sources><div class="method-row-actions"><button class="text-tool" disabled=move||fixed_sources||index==0 on:click=move|_|config.update(|c|{if index>0&&index<c.rows.len(){c.swap_rows(index,index-1);}})>"Move up"</button><button class="text-tool" disabled=move||fixed_sources||config.with(|c|index+1>=c.rows.len()) on:click=move|_|config.update(|c|{if index+1<c.rows.len(){c.swap_rows(index,index+1);}})>"Move down"</button><button class="text-tool" disabled=move||fixed_sources||config.with(|c|c.rows.len()<=1) on:click=move|_|config.update(|c|{if c.rows.len()>1&&index<c.rows.len(){c.remove_row(index);}})>"Remove fragment"</button></div></Show>
    </fieldset>}
}

pub(crate) use crate::cloning_preparations::{prepared_rows, tailing_request};

fn overlap_request(
    state: &MethodsState,
    inputs: Vec<Molecule>,
    target: &Molecule,
) -> Result<(Vec<PlanInput>, collagene_core::assembly_recipe::Recipe), String> {
    use collagene_core::assembly_recipe::Recipe;
    let (overlaps, options) = state.profile()?;
    let (ids, fragments) = prepared_rows(state, inputs.len())?;
    let recipe = Recipe {
        name: format!("Profiled overlap assembly to {}", target.name),
        fragments,
        overlaps,
        circular: state.circular,
        options,
    };
    let borrowed = ids
        .iter()
        .zip(&inputs)
        .map(|(id, m)| (id.as_str(), m))
        .collect::<Vec<_>>();
    collagene_core::assembly_recipe::source_preflight(&borrowed, target, &recipe)
        .map_err(|e| e.to_string())?;
    let bound = ids
        .into_iter()
        .zip(inputs)
        .map(|(id, molecule)| PlanInput {
            revision: molecule.revision(),
            id,
            molecule,
        })
        .collect();
    Ok((bound, recipe))
}

fn ligation_request(
    state: &MethodsState,
    inputs: Vec<Molecule>,
    target: &Molecule,
) -> Result<(Vec<PlanInput>, collagene_core::ligation_recipe::Recipe), String> {
    use collagene_core::ligation_recipe::{Fragment, Recipe};
    let (ids, fragments) = prepared_rows(state, inputs.len())?;
    if state.mode == Mode::GoldenGate {
        let geometry = state.golden_gate.enzyme.profile();
        for fragment in &fragments {
            let selection=match &fragment.preparation {
                collagene_core::assembly_recipe::Preparation::CatalogDigest{selection,..} | collagene_core::assembly_recipe::Preparation::PcrCatalogDigest{selection,..}=>selection,
                _=>return Err("Golden Gate rows require an existing pinned Type IIS digest fragment or an explicit PCR → catalog digest preparation. Use ordered ligation for other preparation chains.".into()),
            };
            if selection.enzymes.len() != 1
                || !selection.enzymes.iter().all(|entry| {
                    entry.enzyme.recognition == geometry.recognition
                        && entry.enzyme.top_cut == geometry.top_cut
                        && entry.enzyme.bottom_cut == geometry.bottom_cut
                })
            {
                return Err("A row's pinned enzyme geometry differs from the chosen Type IIS enzyme. Choose the matching existing digest fragment or apply fresh PCR-site proposals.".into());
            }
        }
    }
    let recipe = Recipe {
        name: format!(
            "{} to {}",
            if state.mode == Mode::GoldenGate {
                "Ordered Golden Gate preparation/ligation"
            } else {
                "Ordered ligation"
            },
            target.name
        ),
        fragments: fragments
            .into_iter()
            .zip(&state.rows)
            .map(|(source, row)| Fragment {
                source,
                treatment: row.treatment.operation(),
            })
            .collect(),
        circular: state.circular,
        allow_nicked: state.allow_nicked,
    };
    let borrowed = ids
        .iter()
        .zip(&inputs)
        .map(|(id, m)| (id.as_str(), m))
        .collect::<Vec<_>>();
    collagene_core::ligation_recipe::source_preflight(&borrowed, target, &recipe)
        .map_err(|e| e.to_string())?;
    let bound = ids
        .into_iter()
        .zip(inputs)
        .map(|(id, molecule)| PlanInput {
            revision: molecule.revision(),
            id,
            molecule,
        })
        .collect();
    Ok((bound, recipe))
}
fn method_status(status: &collagene_core::verification::Status) -> &'static str {
    match status {
        collagene_core::verification::Status::ConsistentUnderAssumptions => {
            "consistent under stated assumptions"
        }
        collagene_core::verification::Status::Incomplete => "incomplete",
        collagene_core::verification::Status::Invalid => "inconsistent",
    }
}
pub(crate) fn method_end(end: Option<&collagene_core::chemistry::End>) -> String {
    use collagene_core::chemistry::{FivePrime, Polarity, ThreePrime};
    let Some(end) = end else {
        return "End unavailable; no terminal chemistry inferred".into();
    };
    let shape = match end.polarity {
        Polarity::Blunt => "blunt".into(),
        Polarity::FivePrime => format!("5′ overhang {} (5′→3′)", end.overhang),
        Polarity::ThreePrime => format!("3′ overhang {} (5′→3′)", end.overhang),
    };
    format!(
        "{shape}; 5′ {}; 3′ {}",
        match end.chemistry.five_prime {
            FivePrime::Phosphate => "phosphate",
            FivePrime::Hydroxyl => "hydroxyl",
            FivePrime::Unknown => "chemistry unknown",
        },
        match end.chemistry.three_prime {
            ThreePrime::Hydroxyl => "hydroxyl",
            ThreePrime::Blocked => "blocked",
            ThreePrime::Unknown => "chemistry unknown",
        }
    )
}
fn method_duplex(summary: &collagene_core::chemistry::DuplexSummary) -> String {
    format!(
        "Top {} bases; bottom {} bases; bottom reference offset {}. {}. Retained nicks: {}. Covalently closed in model: {}.",
        summary.top_length,
        summary.bottom_length,
        summary.bottom_start,
        match summary.ends.as_ref() {
            Some((left, right)) => format!(
                "Left: {}. Right: {}",
                method_end(Some(left)),
                method_end(Some(right))
            ),
            None => "No free ends on the reconstructed circular duplex".into(),
        },
        summary.nicks.len(),
        summary.covalently_closed_in_model
    )
}

#[component]
fn DeclaredEndFields(config: RwSignal<MethodsState>, index: usize, left: bool) -> impl IntoView {
    let end = move || {
        config.with(|c| {
            c.rows
                .get(index)
                .map(|row| {
                    if left {
                        row.left_terminus
                    } else {
                        row.right_terminus
                    }
                })
                .unwrap_or(Terminus::UNKNOWN)
        })
    };
    let side = if left { "Left" } else { "Right" };
    view! {<fieldset><legend>{format!("{side} source end")}</legend><p class="helper">{if left{"Top-strand 5′ terminus and bottom-strand 3′ terminus, before orientation."}else{"Bottom-strand 5′ terminus and top-strand 3′ terminus, before orientation."}}</p>
        <label>"Declared 5′ group"<select aria-label=format!("Fragment {} {side} declared 5 prime",index+1) prop:value=move||match end().five_prime{FivePrime::Phosphate=>"phosphate",FivePrime::Hydroxyl=>"hydroxyl",FivePrime::Unknown=>"unknown"} on:change=move|e|{let value=match event_target_value(&e).as_str(){"phosphate"=>FivePrime::Phosphate,"hydroxyl"=>FivePrime::Hydroxyl,_=>FivePrime::Unknown};config.update(|c|if let Some(row)=c.rows.get_mut(index){if left{row.left_terminus.five_prime=value}else{row.right_terminus.five_prime=value}})}><option value="unknown">"Unknown"</option><option value="phosphate">"Phosphate (P)"</option><option value="hydroxyl">"Hydroxyl (OH)"</option></select></label>
        <label>"Declared 3′ group"<select aria-label=format!("Fragment {} {side} declared 3 prime",index+1) prop:value=move||match end().three_prime{ThreePrime::Hydroxyl=>"hydroxyl",ThreePrime::Blocked=>"blocked",ThreePrime::Unknown=>"unknown"} on:change=move|e|{let value=match event_target_value(&e).as_str(){"hydroxyl"=>ThreePrime::Hydroxyl,"blocked"=>ThreePrime::Blocked,_=>ThreePrime::Unknown};config.update(|c|if let Some(row)=c.rows.get_mut(index){if left{row.left_terminus.three_prime=value}else{row.right_terminus.three_prime=value}})}><option value="unknown">"Unknown"</option><option value="hydroxyl">"Hydroxyl (OH)"</option><option value="blocked">"Blocked"</option></select></label>
    </fieldset>}
}

fn primer_design_request(
    state: &MethodsState,
    inputs: Vec<Molecule>,
    target: &Molecule,
    form: &PlannerForm,
) -> Result<
    (
        Vec<PlanInput>,
        collagene_core::overlap_primer_design::Request,
    ),
    String,
> {
    use collagene_core::overlap_primer_design::{Fragment, Request};
    if state.mode != Mode::Overlap
        || inputs.len() != state.rows.len()
        || inputs.is_empty()
        || inputs.len() > 12
    {
        return Err("Choose 1–12 ordered source regions for overlap assembly.".into());
    }
    if [&form.overlap, &form.minimum, &form.maximum, &form.tm]
        .iter()
        .any(|s| s.len() > 32)
    {
        return Err("Primer settings exceed their field limits.".into());
    }
    let roles = state.role_ids(inputs.len())?;
    let ids = roles
        .iter()
        .map(|role| role.input_id.clone())
        .collect::<Vec<_>>();
    let mut fragments = Vec::with_capacity(inputs.len());
    for (i, (row, molecule)) in state.rows.iter().zip(&inputs).enumerate() {
        if row.interval_start.len() > 20 || row.interval_end.len() > 20 {
            return Err("Invalid primer design interval.".into());
        }
        let start = row
            .interval_start
            .trim()
            .parse::<usize>()
            .ok()
            .and_then(|n| n.checked_sub(1))
            .ok_or_else(|| {
                format!(
                    "Fragment {} start must be a 1-based source position.",
                    i + 1
                )
            })?;
        let end = if row.interval_end.trim().is_empty() {
            molecule.sequence.len()
        } else {
            row.interval_end
                .trim()
                .parse::<usize>()
                .map_err(|_| format!("Invalid end position for fragment {}.", i + 1))?
        };
        fragments.push(Fragment {
            role: roles[i].role.clone(),
            input_id: ids[i].clone(),
            start,
            end,
            reverse: row.reverse,
        });
    }
    let profile_options = ProfileOptions {
        method: state.method,
        scan_repeated_overlaps: state.scan_repeats,
        max_comparisons: state
            .comparisons
            .parse()
            .map_err(|_| "Invalid repeat comparison budget")?,
        max_occurrences: state
            .occurrences
            .parse()
            .map_err(|_| "Invalid repeat occurrence limit")?,
    };
    profile_options.validate().map_err(|e| e.to_string())?;
    let request = Request {
        name: format!("Suggested overlap PCR primers for {}", target.name),
        fragments,
        overlap: form
            .overlap
            .parse()
            .map_err(|_| "Invalid homology length")?,
        circular: state.circular,
        primer_options: collagene_core::primers::DesignOptions {
            min_length: form
                .minimum
                .parse()
                .map_err(|_| "Invalid minimum annealing length")?,
            max_length: form
                .maximum
                .parse()
                .map_err(|_| "Invalid maximum annealing length")?,
            target_tm_c: form.tm.parse().map_err(|_| "Invalid target annealing Tm")?,
            max_pairs: 1,
            ..Default::default()
        },
        profile_options,
    };
    let borrowed = ids
        .iter()
        .zip(&inputs)
        .map(|(id, m)| (id.as_str(), m))
        .collect::<Vec<_>>();
    collagene_core::overlap_primer_design::preflight(&borrowed, target, &request)
        .map_err(|e| e.to_string())?;
    let bound = ids
        .into_iter()
        .zip(inputs)
        .map(|(id, molecule)| PlanInput {
            revision: molecule.revision(),
            id,
            molecule,
        })
        .collect();
    Ok((bound, request))
}

fn golden_gate_request(
    state: &MethodsState,
    inputs: Vec<Molecule>,
    target: &Molecule,
) -> Result<
    (
        Vec<PlanInput>,
        collagene_core::golden_gate_primer_design::Request,
    ),
    String,
> {
    use collagene_core::golden_gate_primer_design::{Fragment, Request};
    if inputs.len() != state.rows.len() || inputs.is_empty() || inputs.len() > 12 {
        return Err("Choose 1–12 ordered original sources.".into());
    }
    let form = &state.golden_gate;
    if form.clamp.len() > 12
        || form.spacer.len() > 2
        || form.forward_length.len() > 4
        || form.reverse_length.len() > 4
    {
        return Err("Type IIS primer settings exceed their field bounds.".into());
    }
    let forward_length = form
        .forward_length
        .parse()
        .map_err(|_| "Enter a forward annealing length.")?;
    let reverse_length = form
        .reverse_length
        .parse()
        .map_err(|_| "Enter a reverse annealing length.")?;
    let roles = state.role_ids(inputs.len())?;
    let ids = roles
        .iter()
        .map(|role| role.input_id.clone())
        .collect::<Vec<_>>();
    let mut fragments = Vec::new();
    for (i, row) in state.rows.iter().enumerate() {
        if row.preparation != OverlapPreparation::PcrCatalogDigest {
            continue;
        }
        if row.interval_start.len() > 20
            || row.interval_end.len() > 20
            || row.right_fusion.len() > 16
        {
            return Err("Primer source interval or fusion field exceeds its bound.".into());
        }
        let start = row
            .interval_start
            .trim()
            .parse::<usize>()
            .ok()
            .and_then(|n| n.checked_sub(1))
            .ok_or("Source interval starts at a 1-based position.")?;
        let end = if row.interval_end.trim().is_empty() {
            inputs[i].sequence.len()
        } else {
            row.interval_end
                .trim()
                .parse()
                .map_err(|_| "Invalid source interval end.")?
        };
        fragments.push(Fragment {
            role: roles[i].role.clone(),
            input_id: ids[i].clone(),
            start,
            end,
            reverse: row.reverse,
            forward_length,
            reverse_length,
            right_overhang: row.right_fusion.clone(),
        });
    }
    if fragments.is_empty() {
        return Err("Choose Exact PCR → pinned catalog digest for at least one row before proposing primers. Existing-site rows remain unchanged.".into());
    }
    let request = Request {
        name: format!("Type IIS primer preparation for {}", target.name),
        fragments,
        enzyme: form.enzyme,
        clamp: form.clamp.clone(),
        spacer: form.spacer.clone(),
    };
    let borrowed = ids
        .iter()
        .zip(&inputs)
        .map(|(id, m)| (id.as_str(), m))
        .collect::<Vec<_>>();
    collagene_core::golden_gate_primer_design::source_preflight(&borrowed, target, &request)
        .map_err(|e| e.to_string())?;
    let inputs = ids
        .into_iter()
        .zip(inputs)
        .map(|(id, molecule)| PlanInput {
            revision: molecule.revision(),
            id,
            molecule,
        })
        .collect();
    Ok((inputs, request))
}

#[component]
fn MethodDiagnostics(
    diagnostics: collagene_core::assembly_diagnostics::Report,
    names: Vec<String>,
    preparations: Vec<String>,
    on_edit: Callback<usize>,
    #[prop(optional_no_strip)] on_source: Option<Callback<usize>>,
) -> impl IntoView {
    use collagene_core::assembly_diagnostics::{JunctionState, RowState, Stage};
    let status = method_status(&diagnostics.status);
    let stage_name = |stage: Stage| match stage {
        Stage::Pcr => "PCR",
        Stage::Digest => "restriction digest",
        Stage::Repair => "end repair",
        Stage::EndTreatment => "terminal treatment",
        Stage::SourceDeclaration => "source assumption",
        Stage::Annealing => "oligo annealing",
        Stage::MaterialReconstruction => "saved material reconstruction",
        Stage::OverlapAssembly => "overlap assembly",
        Stage::LigationAssembly => "ligation",
        Stage::TargetComparison => "independent-target comparison",
    };
    view! {<section class="method-diagnostics" aria-label="Preparation and junction diagnostics">

        <h3>"Preparation and junction review"</h3><p>{status}</p>
        <p class="helper">"These locations come from the fresh typed step graph and retained product prefix. Later blocked rows have not been reconstructed; they are not additional failed experiments. Editing a field clears this evidence; viewing its editor does not change the draft."</p>
        {diagnostics.failure.map(|failure|view!{<article class="diagnostic-failure"><h4>{format!("First blocking result: step {} · {}",failure.step+1,stage_name(failure.stage))}</h4><p>{failure.message}</p><p class="helper">{failure.recovery}</p>
            <div class="method-row-actions">{failure.row.map(|row|view!{<button class="button secondary" on:click=move|_|on_edit.run(row)>"Edit preparation"</button>})}{failure.upstream.map(|row|view!{<button class="text-tool" on:click=move|_|on_edit.run(row)>"Edit upstream preparation"</button>})}{failure.downstream.map(|row|view!{<button class="text-tool" on:click=move|_|on_edit.run(row)>"Edit downstream preparation"</button>})}</div>
        </article>})}
        <div class="diagnostic-table"><table><thead><tr><th>"Original fragment / direction"</th><th>"Declared preparation"</th><th>"Reconstruction"</th><th>"Review"</th></tr></thead><tbody>
            {diagnostics.rows.into_iter().map(|row|{let index=row.index;let state=match row.state{RowState::InputReference=>"Original reference; no preparation step",RowState::Reconstructed=>"Reconstructed",RowState::Failed=>"Failed at its preparation",RowState::Blocked=>"Blocked by an earlier result"};view!{<tr><td>{format!("{} · {} · {}",index+1,names.get(index).map(String::as_str).unwrap_or(&row.input_id),if row.reverse{"reverse after preparation"}else{"forward"})}</td><td>{preparations.get(index).cloned().unwrap_or_else(||"Preparation unavailable".into())}</td><td>{state}<p>{format!("{}/{} declared steps reconstructed",row.completed_steps,row.total_steps)}</p>{row.preparation_step.map(|step|view!{<p>{format!("Preparation step {}",step+1)}</p>})}{row.repair_step.map(|step|view!{<p>{format!("Repair step {}",step+1)}</p>})}</td><td><div class="method-row-actions"><button class="text-tool" aria-label=format!("Edit fragment {} preparation",index+1) on:click=move|_|on_edit.run(index)>"Edit preparation"</button>{on_source.map(|source|view!{<button class="text-tool" on:click=move|_|source.run(index)>"View original source"</button>})}</div></td></tr>}}).collect_view()}
        </tbody></table></div>
        {(!diagnostics.junctions.is_empty()).then(||view!{<h4>"Fresh junction evidence"</h4>})}
        {diagnostics.junctions.into_iter().map(|junction|{let upstream=junction.upstream;let downstream=junction.downstream;let state=match junction.state{JunctionState::ExactHomology=>"Exact homology",JunctionState::HomologyMismatch=>"Incompatible homology",JunctionState::HomologyUnavailable=>"Homology unavailable",JunctionState::EndpointsAvailable=>"Ends available; compatibility is not inferred here",JunctionState::EndpointsUnavailable=>"Physical ends unavailable"};view!{<details><summary>{format!("Junction {} · fragment {} → {}{} · {}",junction.index+1,upstream+1,downstream+1,if junction.closing{" (circle closure)"}else{""},state)}</summary><p>{junction.message}</p><p class="helper">{junction.recovery}</p><div class="method-row-actions"><button class="text-tool" on:click=move|_|on_edit.run(upstream)>"Edit upstream preparation"</button><button class="text-tool" on:click=move|_|on_edit.run(downstream)>"Edit downstream preparation"</button>{on_source.map(|source|view!{<button class="text-tool" on:click=move|_|source.run(upstream)>"View upstream source"</button><button class="text-tool" on:click=move|_|source.run(downstream)>"View downstream source"</button>})}</div></details>}}).collect_view()}
        <ul>{diagnostics.warnings.into_iter().map(|warning|view!{<li class="helper">{warning}</li>}).collect_view()}</ul>
        <details><summary>"Full fresh findings"</summary><ul>{diagnostics.findings.into_iter().map(|finding|view!{<li>{finding}</li>}).collect_view()}</ul></details>
    </section>}
}

#[derive(Clone, PartialEq, Eq)]
struct PoolContext {
    config: MethodsState,
    context: Vec<crate::digest::SourceInfo>,
    target_context: TargetContext,
    target: String,
    revisions: Vec<String>,
    epoch: u64,
}
#[derive(Clone)]
enum PoolReceipt {
    Prepared(Box<collagene_core::prepared_golden_gate::Report>),
    Reaction(Box<collagene_core::golden_gate_reaction::Report>),
}
impl PoolReceipt {
    fn pool(&self) -> Option<&collagene_core::golden_gate::GoldenGateReport> {
        match self {
            Self::Prepared(r) => r.pool.as_ref(),
            Self::Reaction(r) => r.pool.as_ref(),
        }
    }
    fn assessments(&self) -> &[collagene_core::prepared_golden_gate::Assessment] {
        match self {
            Self::Prepared(r) => &r.assessments,
            Self::Reaction(r) => &r.assessments,
        }
    }
    fn context(&self) -> &str {
        match self {
            Self::Prepared(r) => &r.context_revision,
            Self::Reaction(r) => &r.context_revision,
        }
    }
    fn warnings(&self) -> &[String] {
        match self {
            Self::Prepared(r) => &r.warnings,
            Self::Reaction(r) => &r.warnings,
        }
    }
    fn verification(&self) -> &collagene_core::verification::VerificationReport {
        match self {
            Self::Prepared(r) => &r.preview.verification,
            Self::Reaction(r) => &r.verification,
        }
    }
    fn missing(&self) -> Vec<usize> {
        match self {
            Self::Prepared(r) => r.missing_physical.clone(),
            Self::Reaction(r) => r
                .materials
                .iter()
                .filter_map(|m| match m.source {
                    Source::Step(i) if i >= r.verification.products.len() => Some(m.index),
                    _ => None,
                })
                .collect(),
        }
    }
    fn diagnostics(
        &self,
    ) -> Option<Result<collagene_core::assembly_diagnostics::Report, collagene_core::Error>> {
        match self {
            Self::Prepared(r) => Some(collagene_core::assembly_diagnostics::ligation(&r.preview)),
            Self::Reaction(_) => None,
        }
    }
    fn validate(&self, revision: &str) -> Result<(), String> {
        match self {
            Self::Prepared(r) => r.validate_binding(revision),
            Self::Reaction(r) => r.validate_binding(revision),
        }
        .map_err(|e| e.to_string())
    }
    fn selection_revision(
        &self,
        s: &collagene_core::prepared_golden_gate::Selection,
    ) -> Result<String, String> {
        match self {
            Self::Prepared(r) => collagene_core::prepared_golden_gate::selection_revision(r, s),
            Self::Reaction(r) => collagene_core::golden_gate_reaction::selection_revision(r, s),
        }
        .map_err(|e| e.to_string())
    }
    fn select_job(&self, s: collagene_core::prepared_golden_gate::Selection) -> AnalysisJob {
        match self {
            Self::Prepared(r) => AnalysisJob::SelectPreparedGoldenGateCandidate {
                report: r.clone(),
                selection: s,
            },
            Self::Reaction(r) => AnalysisJob::SelectGoldenGateReaction {
                report: r.clone(),
                selection: s,
            },
        }
    }
    fn json(&self) -> Result<Vec<u8>, String> {
        match self {
            Self::Prepared(r) => serde_json::to_vec_pretty(r),
            Self::Reaction(r) => serde_json::to_vec_pretty(r),
        }
        .map_err(|e| e.to_string())
    }
}
/// A projection into the saved draft. Reading or editing options never stores
/// scientific evidence, and restoring a draft cannot restore a successful receipt.
#[derive(Clone, Copy)]
struct PoolSettings(RwSignal<MethodsState>);
impl PoolSettings {
    fn get(&self) -> PoolForm {
        self.0.with(MethodsState::pool_form)
    }
    fn get_untracked(&self) -> PoolForm {
        self.0.with_untracked(MethodsState::pool_form)
    }
    fn with<T>(&self, f: impl FnOnce(&PoolForm) -> T) -> T {
        self.0.with(|c| f(&c.pool_form()))
    }
    fn update(&self, f: impl FnOnce(&mut PoolForm)) {
        self.0.update(|c| {
            let mut next = c.pool_form();
            f(&mut next);
            c.golden_gate_pool = Some(next);
        });
    }
}
struct PoolEvidence {
    report: PoolReceipt,
    form: PoolForm,
    epoch: u64,
}
#[component]
fn PreparedPool(
    proof: Arc<PoolContext>,
    config: RwSignal<MethodsState>,
    guard: Callback<Arc<PoolContext>, Result<(), String>>,
    capture: Callback<Arc<PoolContext>, Result<(Vec<PlanInput>, Molecule), String>>,
    epoch: RwSignal<u64>,
    on_candidate: Callback<Candidate>,
    #[prop(optional_no_strip)] on_finalize: Option<Callback<Candidate, Result<(), String>>>,
    on_edit: Callback<usize>,
    #[prop(optional_no_strip)] on_source: Option<Callback<usize>>,
) -> impl IntoView {
    use collagene_core::prepared_golden_gate as pool;
    let form = PoolSettings(config);
    let proof = StoredValue::new(proof);
    let construct_matrix = proof.with_value(|proof| proof.config.source_roles.is_some());
    let evidence = RwSignal::new(None::<Arc<PoolEvidence>>);
    let page = RwSignal::new(0usize);
    let status = RwSignal::new(String::new());
    let pending = RwSignal::new(false);
    let generation = RwSignal::new(0u32);
    let worker = StoredValue::new_local(None::<Rc<CachedWorker>>);
    let alive = Arc::new(AtomicBool::new(true));
    let cleanup = alive.clone();
    on_cleanup(move || {
        cleanup.store(false, Ordering::Relaxed);
        if let Some(w) = worker.get_value() {
            w.close();
        }
    });
    Effect::new(move |_| {
        let _ = (form.get(), epoch.get());
        generation.update(|g| *g = g.wrapping_add(1));
        if let Some(w) = worker.get_value() {
            w.close();
        }
        worker.set_value(None);
        pending.set(false);
        evidence.set(None);
        page.set(0);
    });
    let pool_guard = Callback::new(move |saved: Arc<PoolEvidence>| -> Result<(), String> {
        guard.run(proof.get_value())?;
        if epoch.get_untracked() != saved.epoch
            || form.get_untracked() != saved.form
            || !evidence.with_untracked(|current| {
                current
                    .as_ref()
                    .is_some_and(|current| Arc::ptr_eq(current, &saved))
            })
        {
            return Err(
                "Prepared pool options or source context changed. Inspect the pool again.".into(),
            );
        }
        Ok(())
    });
    let inspect_alive = alive.clone();
    let inspect = move |_| {
        if construct_matrix {
            status.set("These options belong to a source-matrix construct. Inspect and choose its candidate in the Construct reaction candidate panel below.".into());
            return;
        }
        let options = form.get_untracked();
        let captured_epoch = epoch.get_untracked();
        let prepared = (|| -> Result<_, String> {
            guard.run(proof.get_value())?;
            options.validate()?;
            if options.choices.len() != proof.with_value(|p| p.config.rows.len()) {
                return Err("Reaction substrate choices differ from the source roles.".into());
            }
            let (inputs, target) = capture.run(proof.get_value())?;
            let pool_options =
                options.options(proof.with_value(|p| p.config.golden_gate.enzyme))?;
            let (job, revision) = if options.reaction {
                let config = proof.with_value(|p| p.config.clone());
                let mut mapped = config.clone();
                for (i, row) in mapped.rows.iter_mut().enumerate() {
                    match options.choices.get(i).copied().unwrap_or(0) {
                        1 => {
                            *row = OverlapInput {
                                document: row.document,
                                ..Default::default()
                            }
                        }
                        2 => {
                            row.use_complete_pcr_amplicon()?;
                        }
                        _ => {}
                    }
                }
                let (_, fragments) = prepared_rows(&mapped, inputs.len())?;
                let rows = fragments
                    .into_iter()
                    .enumerate()
                    .map(|(i, f)| {
                        let substrate = match options.choices.get(i).copied().unwrap_or(0) {
                            1 => collagene_core::golden_gate_reaction::Substrate::WholeSource,
                            2 => collagene_core::golden_gate_reaction::Substrate::Prepared {
                                preparation: f.preparation,
                                treatment: None,
                            },
                            _ => collagene_core::golden_gate_reaction::Substrate::Prepared {
                                preparation: f.preparation,
                                treatment: match config.rows[i].treatment {
                                    Treatment::Keep => None,
                                    Treatment::Phosphorylate => Some(EndTreatment::Phosphorylate),
                                    Treatment::Dephosphorylate => {
                                        Some(EndTreatment::Dephosphorylate)
                                    }
                                },
                            },
                        };
                        collagene_core::golden_gate_reaction::Row {
                            role: f.role,
                            input_id: f.input_id,
                            substrate,
                        }
                    })
                    .collect();
                let request = collagene_core::golden_gate_reaction::Request {
                    name: "Declared Golden Gate reaction".into(),
                    rows,
                    options: pool_options,
                };
                let revision = collagene_core::golden_gate_reaction::inspect_revision(
                    &inputs, &target, &request,
                )
                .map_err(|e| e.to_string())?;
                (
                    AnalysisJob::InspectGoldenGateReaction {
                        inputs,
                        target,
                        request,
                    },
                    revision,
                )
            } else {
                let originals = inputs.into_iter().map(|i| i.molecule).collect();
                let (inputs, recipe) =
                    proof.with_value(|p| ligation_request(&p.config, originals, &target))?;
                let request = pool::Request {
                    recipe,
                    options: pool_options,
                };
                let revision = pool::inspect_revision(&inputs, &target, &request)
                    .map_err(|e| e.to_string())?;
                (
                    AnalysisJob::InspectPreparedGoldenGatePool {
                        inputs,
                        target,
                        request,
                    },
                    revision,
                )
            };
            let w = Rc::new(CachedWorker::new()?);
            Ok((job, revision, w))
        })();
        let (job, revision, w) = match prepared {
            Ok(value) => value,
            Err(error) => {
                status.set(error);
                return;
            }
        };
        if let Some(old) = worker.get_value() {
            old.close();
        }
        worker.set_value(Some(w.clone()));
        let token = generation.get_untracked().wrapping_add(1);
        generation.set(token);
        pending.set(true);
        evidence.set(None);
        page.set(0);
        status.set(
            "Reconstructing every preparation and exploring the finite raw-substrate pool…".into(),
        );
        let alive = inspect_alive.clone();
        spawn_local(async move {
            let response = w
                .ask(Request {
                    id: token,
                    revision: revision.clone(),
                    job,
                })
                .await;
            w.close();
            if !alive.load(Ordering::Relaxed)
                || generation.get_untracked() != token
                || epoch.get_untracked() != captured_epoch
                || form.get_untracked() != options
            {
                return;
            }
            pending.set(false);
            worker.set_value(None);
            let outcome = (|| -> Result<_, String> {
                guard.run(proof.get_value())?;
                let response = response?;
                if response.id != token || response.revision != revision {
                    return Err("Prepared pool response identity changed.".into());
                }
                let report = match response.result {
                    AnalysisResult::PreparedGoldenGatePool(report) if !options.reaction => {
                        PoolReceipt::Prepared(Box::new(report))
                    }
                    AnalysisResult::GoldenGateReaction(report) if options.reaction => {
                        PoolReceipt::Reaction(Box::new(report))
                    }
                    AnalysisResult::Error(error) => return Err(error),
                    _ => return Err("Unexpected prepared pool result.".into()),
                };
                report.validate(&revision)?;
                Ok(Arc::new(PoolEvidence {
                    report,
                    form: options,
                    epoch: captured_epoch,
                }))
            })();
            match outcome {
                Ok(bound) => {
                    let withheld = bound.report.pool().is_none();
                    evidence.set(Some(bound));
                    status.set(if withheld{"Prepared physical state is missing; pool search is withheld. Review the fresh preparation diagnostics."}else{"Finite pool evidence ready. No candidate is automatically preferred; review alternatives and the independent target assessment."}.into());
                }
                Err(error) => status.set(error),
            }
        });
    };
    let select_alive = alive.clone();
    let select = Callback::new(
        move |(saved, index, finalize): (Arc<PoolEvidence>, usize, bool)| {
            let prepared = (|| -> Result<_, String> {
                if construct_matrix {
                    return Err("Use the construct's exact source-bound pool discovery and candidate choice below.".into());
                }
                pool_guard.run(saved.clone())?;
                let candidate = saved
                    .report
                    .pool()
                    .and_then(|report| report.candidates.get(index))
                    .ok_or("The selected pool candidate is unavailable.")?;
                if !candidate.restriction_site_free {
                    return Err("Only site-free candidates can be selected. Retained/recreated recognition remains visible for diagnosis.".into());
                }
                if finalize && on_finalize.is_none() {
                    return Err("Output creation is unavailable in this standalone context; review the candidate instead.".into());
                }
                let selection = pool::Selection {
                    context_revision: saved.report.context().to_string(),
                    candidate_index: index,
                    directed_sequence_key: candidate.directed_sequence_key.clone(),
                    physical_sequence_key: candidate.physical_sequence_key.clone(),
                };
                let revision = saved
                    .report
                    .selection_revision(&selection)
                    .map_err(|e| e.to_string())?;
                let w = Rc::new(CachedWorker::new()?);
                Ok((selection, revision, w))
            })();
            let (selection, revision, w) = match prepared {
                Ok(value) => value,
                Err(error) => {
                    status.set(error);
                    return;
                }
            };
            if let Some(old) = worker.get_value() {
                old.close();
            }
            worker.set_value(Some(w.clone()));
            let token = generation.get_untracked().wrapping_add(1);
            generation.set(token);
            pending.set(true);
            status.set(
                "Binding the exact candidate keys into an unverified source-backed plan…".into(),
            );
            let alive = select_alive.clone();
            spawn_local(async move {
                let response = w
                    .ask(Request {
                        id: token,
                        revision: revision.clone(),
                        job: saved.report.select_job(selection.clone()),
                    })
                    .await;
                w.close();
                if !alive.load(Ordering::Relaxed) || generation.get_untracked() != token {
                    return;
                }
                pending.set(false);
                worker.set_value(None);
                let outcome = (|| -> Result<(), String> {
                    pool_guard.run(saved.clone())?;
                    let response = response?;
                    if response.id != token || response.revision != revision {
                        return Err("Pool selection response identity changed.".into());
                    }
                    let result = match response.result {
                        AnalysisResult::PreparedGoldenGateCandidate(result)
                        | AnalysisResult::GoldenGateReactionCandidate(result) => result,
                        AnalysisResult::Error(error) => return Err(error),
                        _ => return Err("Unexpected prepared pool selection result.".into()),
                    };
                    result
                        .validate_binding(&revision)
                        .map_err(|e| e.to_string())?;
                    if result.parent_context != saved.report.context()
                        || result.selection.context_revision != selection.context_revision
                        || result.selection.candidate_index != selection.candidate_index
                        || result.selection.directed_sequence_key != selection.directed_sequence_key
                        || result.selection.physical_sequence_key != selection.physical_sequence_key
                    {
                        return Err(
                            "Selected candidate keys or parent pool context changed.".into()
                        );
                    }
                    let candidate = Candidate::new(result.plan);
                    if finalize {
                        on_finalize
                            .ok_or("Output creation context is unavailable.")?
                            .run(candidate)
                    } else {
                        on_candidate.run(candidate);
                        Ok(())
                    }
                })();
                match outcome{Ok(())=>status.set("Unverified candidate submitted for fresh plan reconstruction. A target mismatch can remain Invalid; no mixture uniqueness or experimental success was established.".into()),Err(error)=>status.set(error)}
            });
        },
    );
    view! {<section class="prepared-gate-pool" aria-label="Prepared Golden Gate pool alternatives"><h3>"Explore Golden Gate alternatives"</h3>
        <p class="helper">"Inspect selected prepared fragments, or explicitly choose what enters the reaction for every role. The finite search reports competing products, recutting and uncut circular carryover without selecting a winner or predicting experimental yield."</p>
        <p class="helper">"Reaction choices and search limits are kept with the assembly draft. Reopen the draft and inspect again to reconstruct current evidence."</p>
        {construct_matrix.then(||view!{<p class="helper">"These are this construct's pool settings. Use the Construct reaction candidate panel below for exact source-bound discovery, selection, product/junction review and output creation. Ordered method evidence is unavailable until you explicitly return to ordered prepared assembly."</p>})}
        <label class="checkbox-label"><input type="checkbox" aria-label="Inspect whole Golden Gate reaction" prop:checked=move||form.with(|f|f.reaction) on:change=move|e|form.update(|f|f.reaction=event_target_checked(&e))/>"Declare reaction substrates (whole DNA / whole amplicons / prepared fragments)"</label>
        <Show when=move||form.with(|f|f.reaction)><p class="helper">"Choose what enters the tube for every role. Whole source includes its circular backbone and all cut pieces. Whole amplicon uses the current primers before digestion and omits post-digest repair/treatment. Retaining a prepared fragment keeps its current explicit operations. Primers and the intended target are unchanged."</p>
        {proof.with_value(|p|p.config.rows.iter().enumerate().map(|(i,row)|{let pcr=row.is_pcr();view!{<label>{format!("{}: {}",i+1,p.context.get(i).map(|c|c.name.as_str()).unwrap_or("source"))}<select aria-label=format!("Reaction substrate {}",i+1) prop:value=move||form.with(|f|f.choices.get(i).copied().unwrap_or(0).to_string()) on:change=move|e|form.update(|f|if let Some(choice)=f.choices.get_mut(i){*choice=event_target_value(&e).parse().unwrap_or(0)})><option value="0">"Current selected prepared fragment"</option><option value="1">"Whole original source"</option><option value="2" disabled=!pcr>"Complete PCR amplicon before digestion"</option></select></label>}}).collect_view())}</Show>
        <details><summary>"Finite pool scope and limits"</summary><p>{proof.with_value(|p|format!("Recut profile: {:?}. Bounds do not model kinetics, partial digestion, repeated copies or all real mixture outcomes.",p.config.golden_gate.enzyme))}</p>
            <label class="checkbox-label"><input type="checkbox" aria-label="Prepared pool circular products" prop:checked=move||form.with(|f|f.circular) on:change=move|e|form.update(|f|f.circular=event_target_checked(&e))/>"Search circular products (otherwise linear)"</label>
            <label class="checkbox-label"><input type="checkbox" aria-label="Prepared pool both raw orientations" prop:checked=move||form.with(|f|f.both) on:change=move|e|form.update(|f|f.both=event_target_checked(&e))/>"Explore both raw substrate orientations"</label><p class="helper">"If unchecked, each prepared substrate uses its original, unreversed preparation orientation. Ordered-row reverse checkboxes do not constrain this pool search."</p>
            <div class="primer-grid"><label>"Minimum fragments (1 includes background)"<input type="number" min="1" aria-label="Prepared pool minimum fragments" prop:value=move||form.with(|f|f.minimum.clone()) on:input=move|e|form.update(|f|f.minimum=event_target_value(&e))/></label><label>"Maximum fragments"<input type="number" min="1" aria-label="Prepared pool maximum fragments" prop:value=move||form.with(|f|f.maximum.clone()) on:input=move|e|form.update(|f|f.maximum=event_target_value(&e))/></label><label>"Product limit"<input type="number" min="1" aria-label="Prepared pool product limit" prop:value=move||form.with(|f|f.products.clone()) on:input=move|e|form.update(|f|f.products=event_target_value(&e))/></label><label>"Search-state limit"<input type="number" min="1" aria-label="Prepared pool state limit" prop:value=move||form.with(|f|f.states.clone()) on:input=move|e|form.update(|f|f.states=event_target_value(&e))/></label></div>
        </details>
        <button class="button secondary" disabled=move||pending.get()||construct_matrix on:click=inspect>{move||if construct_matrix{"Use construct reaction panel below"}else if form.with(|f|f.reaction){"Inspect declared reaction"}else{"Inspect prepared pool"}}</button><button class="text-tool" disabled=move||!pending.get() on:click=move|_|{generation.update(|g|*g=g.wrapping_add(1));if let Some(w)=worker.get_value(){w.close();}worker.set_value(None);pending.set(false);status.set("Prepared-pool computation cancelled.".into());}>"Cancel prepared pool"</button><p role="status">{move||status.get()}</p>
        {move||evidence.get().map(|bound|{let request=StoredValue::new(bound.clone());let diagnostics=bound.report.diagnostics();let names=proof.with_value(|p|p.context.iter().map(|s|s.name.clone()).collect::<Vec<_>>());let preparations=proof.with_value(|p|p.config.rows.iter().map(OverlapInput::summary).collect::<Vec<_>>());view!{<section class="prepared-gate-report">
            <h4>"Fresh pool reconstruction"</h4><p>{format!("Native diagnostic assessment: {}. For a reaction this is the index0 probe, not a selected winner.",method_status(&bound.report.verification().status))}</p>
            {(!bound.report.missing().is_empty()).then(||view!{<p role="alert">{format!("Pool withheld: physical substrate state is unavailable for rows {}. No partial search is presented as complete.",bound.report.missing().iter().map(|i|(i+1).to_string()).collect::<Vec<_>>().join(", "))}</p>})}
            <ul>{bound.report.warnings().iter().map(|warning|view!{<li>{warning.clone()}</li>}).collect_view()}</ul>
            {match diagnostics{Some(Ok(diagnostics))=>view!{<details><summary>"Recomputed row preparation diagnostics"</summary><MethodDiagnostics diagnostics names preparations on_edit on_source/></details>}.into_any(),Some(Err(error))=>view!{<p class="helper">{format!("Structured diagnostics unavailable: {error}")}</p>}.into_any(),None=>view!{<details><summary>"Native preparation and probe findings"</summary>{bound.report.verification().findings.iter().map(|finding|view!{<p>{finding.clone()}</p>}).collect_view()}</details>}.into_any()}}
            {bound.report.pool().map(|pool|{let count=pool.candidates.len();let pages=count.max(1).div_ceil(50);let start=page.get().min(pages-1)*50;view!{<h4>{format!("{} directed candidates · {} physical sequence groups",count,pool.stable_physical_product_count)}</h4><p>{if pool.enumeration_complete{"Declared finite substrate-pool search completed within the selected options."}else{"Search incomplete: additional declared substrate-pool alternatives may exist."}}</p><p class="helper">{format!("States examined {}. Repeated-copy/real-mixture exhaustiveness is outside this model. Uncut circular substrate carryover roles: {} (separate from assembled candidates). Sequence grouping is not a unique source/annotation history.",pool.states_examined,if pool.carryover_sources.is_empty(){"none".into()}else{pool.carryover_sources.iter().map(|row|(row+1).to_string()).collect::<Vec<_>>().join(", ")})}</p><ul>{pool.truncation_reasons.iter().chain(pool.warnings.iter()).map(|warning|view!{<li>{warning.clone()}</li>}).collect_view()}</ul>
                {pool.candidates.iter().enumerate().skip(start).take(50).map(|(index,candidate)|{let assessment=bound.report.assessments().iter().find(|a|a.candidate_index==index);let allowed=candidate.restriction_site_free;view!{<article class="prepared-gate-candidate"><h4>{format!("Candidate {} · {} top-strand bases · {} pool fragments",index+1,candidate.duplex.top.len(),candidate.path.len())}</h4><p>{if allowed{"Recognition-site-free under the chosen profile"}else{"Recognition retained/recreated — diagnosis only"}}</p><p>{format!("Independent target match: {} · used recipe rows (unique): {}",assessment.map(|a|a.target_matches.to_string()).unwrap_or_else(||"unavailable".into()),assessment.map(|a|a.source_roles.iter().map(|&row|proof.with_value(|p|format!("{} ({})",row+1,p.context.get(row).map(|s|s.name.as_str()).unwrap_or("unknown source")))).collect::<Vec<_>>().join(", ")).unwrap_or_else(||"unavailable".into()))}</p><p class="helper">{format!("Raw path: {}",candidate.path.iter().map(|step|pool.fragments.get(step.fragment).map(|fragment|format!("pool fragment {} · row {} / digest {}{}",step.fragment+1,fragment.source_index+1,fragment.fragment_index+1,if step.reverse{" reverse"}else{" forward"})).unwrap_or_else(||"unavailable pool fragment".into())).collect::<Vec<_>>().join(" → "))}</p><pre>{short(candidate.duplex.top.as_str())}</pre><p>{format!("{} recut/recognition observations",candidate.recut_sites.len())}</p><PoolCandidateDetails bound=bound.clone() index/>
                    <button class="button secondary" disabled=move||pending.get()||!allowed on:click=move|_|select.run((request.get_value(),index,false))>"Review pool candidate"</button>{on_finalize.map(|_|view!{<button class="button primary" disabled=move||pending.get()||!allowed on:click=move|_|select.run((request.get_value(),index,true))>"Create outputs from this candidate"</button>})}
                </article>}}).collect_view()}
                <p>{format!("Page {} of {} · candidates {}–{} of {}",page.get().min(pages-1)+1,pages,if count==0{0}else{start+1},(start+50).min(count),count)}</p><button class="text-tool" disabled=move||page.get()==0 on:click=move|_|page.update(|p|*p=p.saturating_sub(1))>"Previous pool candidates"</button><button class="text-tool" disabled=move||{page.get()+1>=pages} on:click=move|_|page.update(|p|*p=(*p+1).min(pages-1))>"Next pool candidates"</button>
            }})}
            <button class="button secondary" on:click=move|_|{let saved=request.get_value();let result=pool_guard.run(saved.clone()).and_then(|()|saved.report.json()).and_then(|bytes|download("prepared-golden-gate-pool.json","application/json",&bytes));if let Err(error)=result{status.set(error);}}>"Download full prepared-pool evidence"</button>
        </section>}})}
    </section>}
}

#[component]
fn PoolCandidateDetails(bound: Arc<PoolEvidence>, index: usize) -> impl IntoView {
    let bound = StoredValue::new(bound);
    let page = RwSignal::new(0usize);
    let total = bound.with_value(|saved| {
        saved
            .report
            .pool()
            .and_then(|pool| pool.candidates.get(index))
            .map(|candidate| candidate.recut_sites.len())
            .unwrap_or(0)
    });
    let pages = total.max(1).div_ceil(25);
    view! {<details><summary>"Candidate identity and recut evidence"</summary>
        {bound.with_value(|saved|saved.report.pool().and_then(|pool|pool.candidates.get(index)).map(|candidate|view!{<p>{format!("Directed key {} · physical key {}",candidate.directed_sequence_key,candidate.physical_sequence_key)}</p>}))}
        {move||bound.with_value(|saved|saved.report.pool().and_then(|pool|pool.candidates.get(index)).map(|candidate|candidate.recut_sites.iter().skip(page.get()*25).take(25).map(|site|view!{<p>{format!("Recognition at {}{} · paired {} · cleavable under model {}",site.recognition_start+1,if site.reverse{" reverse"}else{" forward"},site.fully_double_stranded,site.cleavable_under_model)}</p>}).collect_view()))}
        <p>{move||format!("{} recut observations · page {} of {}",total,page.get()+1,pages)}</p>
        <button class="text-tool" disabled=move||page.get()==0 on:click=move|_|page.update(|p|*p=p.saturating_sub(1))>"Previous recut observations"</button><button class="text-tool" disabled=move||{page.get()+1>=pages} on:click=move|_|page.update(|p|*p=(*p+1).min(pages-1))>"Next recut observations"</button>
    </details>}
}
