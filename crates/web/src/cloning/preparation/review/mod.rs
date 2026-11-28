//! Source-bound primer proposals and prepared-method evidence review.
//! One worker lifecycle invalidates proposals and evidence together; receipts own
//! their admission and application invariants, independently of presentation.
use crate::cloning::preparation::{
    cloning_methods_model::*,
    cloning_preparations::{ligation_request, overlap_request, tailing_request},
    diagnostics::{MethodDiagnostics, guidance, method_duplex, method_end, method_status, short},
    pool::PreparedPoolSection,
};
use crate::platform::worker::CachedWorker;
use collagene_core::{
    cloning::{
        ligation::tailing_recipe::{Action, ProfileStatus},
        overlap::profiles::{ProfileOptions, ProfileReport},
    },
    documents::transfer::Candidate,
    planning::plan::PlanInput,
    protocol::{
        request::Request,
        schema::{AnalysisJob, AnalysisResult},
    },
    sequence::analysis::Primer,
};
use collagene_model::Molecule;
use leptos::prelude::*;
use serde::{Deserialize, Serialize};
use std::{
    rc::Rc,
    sync::{
        Arc,
        atomic::{AtomicBool, Ordering},
    },
};
use wasm_bindgen_futures::spawn_local;
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
    Document(crate::sources::identity::SourceInfo),
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
    report: collagene_core::cloning::overlap::primer_design::Report,
    config: MethodsState,
    form: PlannerForm,
    context: Vec<crate::sources::identity::SourceInfo>,
    target: String,
    target_context: TargetContext,
    epoch: u64,
}
struct GoldenGateSuggestions {
    report: collagene_core::cloning::golden_gate::primers::Report,
    config: MethodsState,
    context: Vec<crate::sources::identity::SourceInfo>,
    target_context: TargetContext,
    epoch: u64,
}
#[derive(Clone)]
struct Evidence {
    config: MethodsState,
    revisions: Vec<String>,
    context: Vec<crate::sources::identity::SourceInfo>,
    target_context: TargetContext,
    target: String,
    epoch: u64,
    profile: Option<ProfileReport>,
    overlap: Option<collagene_core::cloning::assembly::recipe::Preview>,
    ligation: Option<collagene_core::cloning::ligation::recipe::Preview>,
    tailing: Option<collagene_core::cloning::ligation::tailing_recipe::Preview>,
}

impl GoldenGateSuggestions {
    fn is_current(
        &self,
        config: &MethodsState,
        context: Option<Vec<crate::sources::identity::SourceInfo>>,
        target_context: &TargetContext,
        epoch: u64,
    ) -> bool {
        config == &self.config
            && context == Some(self.context.clone())
            && target_context == &self.target_context
            && epoch == self.epoch
    }
    fn apply(&self) -> Result<MethodsState, String> {
        let mut next = self.config.clone();
        let roles = match next.role_ids(next.rows.len()) {
            Ok(roles) => roles,
            Err(error) => {
                return Err(error);
            }
        };
        for fragment in &self.report.fragments {
            let Some(index) = roles
                .iter()
                .position(|role| fragment.role == role.role && fragment.input_id == role.input_id)
            else {
                return Err("Proposal has an unknown source role; nothing applied.".into());
            };
            let collagene_core::cloning::assembly::recipe::Preparation::PcrCatalogDigest {
                conditions,
                template_strand,
                forward,
                reverse,
                selection,
                fragment_index,
                repair,
            } = &fragment.preparation
            else {
                return Err(
                    "Proposal is not an explicit PCR/digest chain; nothing applied.".into(),
                );
            };
            let row = &mut next.rows[index];
            if row.preparation != OverlapPreparation::PcrCatalogDigest
                || row.reverse != fragment.reverse
            {
                return Err("Preparation or orientation changed; nothing applied.".into());
            }
            row.template_strand = *template_strand;
            row.forward_name = forward.name.clone();
            row.forward = forward.annealing.as_str().into();
            row.forward_tail = forward.tail.clone();
            row.reverse_name = reverse.name.clone();
            row.reverse_primer = reverse.annealing.as_str().into();
            row.reverse_tail = reverse.tail.clone();
            row.catalog_selection = Some(selection.clone());
            if let Some(conditions) = conditions {
                row.reaction_conditions =
                    crate::restriction::reaction_conditions_model::Draft::from_saved(Some(
                        conditions.clone(),
                    ));
            }
            row.fragment = (fragment_index + 1).to_string();
            row.repair = *repair;
        }
        Ok(next)
    }
}
impl Suggestions {
    fn is_current(
        &self,
        config: &MethodsState,
        form: &PlannerForm,
        context: Option<Vec<crate::sources::identity::SourceInfo>>,
        target_context: &TargetContext,
        target: &str,
        epoch: u64,
    ) -> bool {
        config == &self.config
            && form == &self.form
            && context == Some(self.context.clone())
            && target_context == &self.target_context
            && target == self.target
            && epoch == self.epoch
    }
    fn apply(&self) -> Result<MethodsState, String> {
        if self.report.recipe.fragments.len() != self.config.rows.len() {
            return Err("Suggested primer rows do not match the original fragment order.".into());
        }
        let mut updated = self.config.clone();
        let roles = match updated.role_ids(updated.rows.len()) {
            Ok(roles) => roles,
            Err(error) => {
                return Err(error);
            }
        };
        for (i, fragment) in self.report.recipe.fragments.iter().enumerate() {
            if fragment.input_id != roles[i].input_id
                || fragment.role != roles[i].role
                || fragment.reverse != updated.rows[i].reverse
            {
                return Err(
                    "Suggested source order or orientation changed; nothing applied.".into(),
                );
            }
            let collagene_core::cloning::assembly::recipe::Preparation::ExactPcr {
                template_strand,
                forward,
                reverse,
            } = &fragment.preparation
            else {
                return Err(
                    "Suggestion is not an explicit PCR preparation; nothing applied.".into(),
                );
            };
            let row = &mut updated.rows[i];
            row.preparation = OverlapPreparation::ExactPcr;
            row.template_strand = *template_strand;
            row.forward_name = forward.name.clone();
            row.forward = forward.annealing.as_str().into();
            row.forward_tail = forward.tail.clone();
            row.reverse_name = reverse.name.clone();
            row.reverse_primer = reverse.annealing.as_str().into();
            row.reverse_tail = reverse.tail.clone();
        }
        updated.overlaps = self
            .report
            .recipe
            .overlaps
            .iter()
            .map(usize::to_string)
            .collect::<Vec<_>>()
            .join(",");
        Ok(updated)
    }
}
impl Evidence {
    fn is_current(
        &self,
        config: &MethodsState,
        context: Option<Vec<crate::sources::identity::SourceInfo>>,
        revisions: Option<Vec<String>>,
        target_context: &TargetContext,
        target: &str,
        epoch: u64,
    ) -> bool {
        config == &self.config
            && context == Some(self.context.clone())
            && revisions == Some(self.revisions.clone())
            && target_context == &self.target_context
            && target == self.target
            && epoch == self.epoch
    }
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
pub(super) fn MethodsReview(
    molecule: Memo<Molecule>,
    get_query: Callback<usize, Option<Molecule>>,
    source_info: Callback<usize, Option<crate::sources::identity::SourceInfo>>,
    target_info: Memo<TargetContext>,
    config: RwSignal<MethodsState>,
    epoch: RwSignal<u64>,
    on_candidate: Callback<Candidate>,
    on_products: Callback<Vec<Molecule>, Result<(), String>>,
    #[prop(optional_no_strip)] on_finalize: Option<Callback<Candidate, Result<(), String>>>,
    #[prop(optional_no_strip)] on_source: Option<Callback<(usize, usize), Result<(), String>>>,
    focus_row: Callback<usize>,
    #[prop(optional_no_strip)] active_row: Option<Memo<usize>>,
    fixed_sources: bool,
    status: RwSignal<String>,
    disable_plain_primer_design: bool,
) -> impl IntoView {
    let source_indices = Memo::new(move |_| config.with(|state| state.source_indices().ok()));
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
        let revision = match collagene_core::cloning::golden_gate::primers::preview_revision(
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
        if !saved.is_current(
            &config.get_untracked(),
            source_context.get_untracked(),
            &target_info.get_untracked(),
            epoch.get_untracked(),
        ) {
            status.set("Sources or settings changed; request fresh Type IIS primers.".into());
            return;
        }
        let next = match saved.apply() {
            Ok(next) => next,
            Err(error) => {
                status.set(error);
                return;
            }
        };
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
        let (inputs, request) = match selected.request(&state, inputs, &bound_target) {
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
        let revision = match collagene_core::cloning::overlap::primer_design::design_revision(
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
            || !proposal.is_current(
                &config.get_untracked(),
                &planner_form.get_untracked(),
                source_context.get_untracked(),
                &target_info.get_untracked(),
                &target.get_untracked(),
                epoch.get_untracked(),
            )
        {
            status.set("Sources, target or settings changed; request fresh primer suggestions before applying.".into());
            return;
        }
        let updated = match proposal.apply() {
            Ok(updated) => updated,
            Err(error) => {
                status.set(error);
                return;
            }
        };
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
                    .is_none_or(|source| source.binding != *pin.event())
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
                if state.mode == Mode::Tailing {
                    let (bound_inputs, recipe) =
                        tailing_request(&state, inputs, &bound_target_molecule)?;
                    let revision =
                        collagene_core::cloning::ligation::tailing_recipe::preview_revision(
                            &bound_inputs,
                            &bound_target_molecule,
                            &recipe,
                        )
                        .map_err(|e| e.to_string())?;
                    let response = ask(
                        &worker,
                        AnalysisJob::PreviewTailingRecipe {
                            inputs: bound_inputs.clone(),
                            target: bound_target_molecule.clone(),
                            recipe: recipe.clone(),
                        },
                        revision.clone(),
                        token,
                    )
                    .await?;
                    if let AnalysisResult::TailingRecipePreview(preview) = response {
                        preview
                            .validate_for(&bound_inputs, &bound_target_molecule, &recipe, &revision)
                            .map_err(|e| e.to_string())?;
                        Ok((None, None, None, Some(preview)))
                    } else {
                        Err("Unexpected tailing/ligation response".into())
                    }
                } else if matches!(state.mode, Mode::Ligation | Mode::GoldenGate) {
                    let (bound_inputs, recipe) =
                        ligation_request(&state, inputs, &bound_target_molecule)?;
                    let revision = collagene_core::cloning::ligation::recipe::preview_revision(
                        &bound_inputs,
                        &bound_target_molecule,
                        &recipe,
                        state.cleavage_policy(),
                    )
                    .map_err(|e| e.to_string())?;
                    let response = ask(
                        &worker,
                        AnalysisJob::PreviewLigationRecipe {
                            policy: state.cleavage_policy(),
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
                        Ok((None, None, Some(preview), None))
                    } else {
                        Err("Unexpected ordered ligation response".into())
                    }
                } else if let Some((overlaps, options)) = profile {
                    let _ = (overlaps, options); // shared recipe reuses the exact parsed profile.
                    let (bound_inputs, recipe) =
                        overlap_request(&state, inputs, &bound_target_molecule)?;
                    let (revision, job) = if let Some(form) =
                        state.annealing.as_ref().filter(|f| f.enabled)
                    {
                        let policy = form.policy(
                            recipe.overlaps.clone(),
                            recipe.circular,
                            recipe.fragments.len(),
                        )?;
                        let revision =
                            collagene_core::cloning::assembly::recipe::annealing_preview_revision(
                                &bound_inputs,
                                &bound_target_molecule,
                                &recipe,
                                &policy,
                            )
                            .map_err(|e| e.to_string())?;
                        (
                            revision,
                            AnalysisJob::PreviewAnnealingRecipe {
                                inputs: bound_inputs,
                                target: bound_target_molecule,
                                recipe,
                                policy,
                            },
                        )
                    } else if let Some(form) = state.homology.as_ref().filter(|f| f.enabled) {
                        let policy = form.policy(
                            recipe.overlaps.clone(),
                            recipe.circular,
                            recipe.fragments.len(),
                        )?;
                        let revision =
                            collagene_core::cloning::assembly::recipe::homology_preview_revision(
                                &bound_inputs,
                                &bound_target_molecule,
                                &recipe,
                                &policy,
                            )
                            .map_err(|e| e.to_string())?;
                        (
                            revision,
                            AnalysisJob::PreviewHomologyRecipe {
                                inputs: bound_inputs,
                                target: bound_target_molecule,
                                recipe,
                                policy,
                            },
                        )
                    } else {
                        let revision = collagene_core::cloning::assembly::recipe::preview_revision(
                            &bound_inputs,
                            &bound_target_molecule,
                            &recipe,
                        )
                        .map_err(|e| e.to_string())?;
                        (
                            revision,
                            AnalysisJob::PreviewOverlapRecipe {
                                inputs: bound_inputs,
                                target: bound_target_molecule,
                                recipe,
                            },
                        )
                    };
                    let response = ask(&worker, job, revision.clone(), token).await?;
                    if let AnalysisResult::OverlapRecipePreview(preview) = response {
                        preview
                            .validate_binding(&revision)
                            .map_err(|e| e.to_string())?;
                        Ok((preview.profile.clone(), Some(preview), None, None))
                    } else {
                        Err("Unexpected overlap preparation response".into())
                    }
                } else {
                    Err("Choose a supported explicit cloning recipe.".into())
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
                Ok((profile, overlap, ligation, tailing)) => {
                    evidence.set(Some(Arc::new(Evidence {
                        config: state,
                        revisions: revs,
                        context: bound_context,
                        target_context: bound_target_context,
                        target: bound_target,
                        epoch: bound_epoch,
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
        }) || !saved.is_current(
            &config.get_untracked(),
            source_context.get_untracked(),
            source_revisions.get_untracked(),
            &target_info.get_untracked(),
            &target.get_untracked(),
            epoch.get_untracked(),
        ) {
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
        if state.mode == Mode::Tailing {
            let reviewed = (|| -> Result<Candidate, String> {
                let (bound, recipe) = tailing_request(&state, inputs, &current_target)?;
                let context = collagene_core::cloning::ligation::tailing_recipe::preview_revision(
                    &bound,
                    &current_target,
                    &recipe,
                )
                .map_err(|error| error.to_string())?;
                let preview = proof
                    .tailing
                    .as_ref()
                    .ok_or("Compute the complete tailing recipe before reviewing a candidate.")?;
                preview
                    .validate_for(&bound, &current_target, &recipe, &context)
                    .map_err(|error| error.to_string())?;
                Ok(Candidate::new(preview.compiled.plan.clone()))
            })();
            match reviewed {
                Ok(candidate) => on_candidate.run(candidate),
                Err(error) => status.set(error),
            }
            return;
        } else if matches!(state.mode, Mode::Ligation | Mode::GoldenGate) {
            let (bound, recipe) = match ligation_request(&state, inputs, &current_target) {
                Ok(v) => v,
                Err(e) => {
                    status.set(e);
                    return;
                }
            };
            match collagene_core::cloning::ligation::recipe::compile(
                &bound,
                &current_target,
                &recipe,
                state.cleavage_policy(),
            ) {
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
            let compiled = if let Some(form) = state.annealing.as_ref().filter(|form| form.enabled)
            {
                form.policy(
                    recipe.overlaps.clone(),
                    recipe.circular,
                    recipe.fragments.len(),
                )
                .and_then(|policy| {
                    collagene_core::cloning::assembly::recipe::compile_annealing(
                        &bound,
                        &current_target,
                        &recipe,
                        &policy,
                    )
                    .map_err(|e| e.to_string())
                })
            } else if let Some(form) = state.homology.as_ref().filter(|form| form.enabled) {
                form.policy(
                    recipe.overlaps.clone(),
                    recipe.circular,
                    recipe.fragments.len(),
                )
                .and_then(|policy| {
                    collagene_core::cloning::assembly::recipe::compile_homology(
                        &bound,
                        &current_target,
                        &recipe,
                        &policy,
                    )
                    .map_err(|e| e.to_string())
                })
            } else {
                collagene_core::cloning::assembly::recipe::compile(&bound, &current_target, &recipe)
                    .map_err(|e| e.to_string())
            };
            match compiled {
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
    });
    view! {<>
        <Show when=move||config.get().mode==Mode::GoldenGate&&!disable_plain_primer_design><section class="golden-gate-primer-planner"><h3>"2. Prepare Type IIS fragments"</h3>
            <p class="helper">"Use existing pinned digest fragments, or choose Exact PCR → pinned catalog digest for rows needing new sites. Proposals use each chosen original source region and explicit right fusion. They do not replace existing-site rows or choose a target. This ordered reconstruction does not enumerate every mixture outcome."</p>
            <details class="assembly-primer-settings"><summary>"Type IIS primer design settings"</summary><div class="primer-grid">
                <label>"Type IIS enzyme"<select aria-label="Golden Gate preparation enzyme" prop:value=move||format!("{:?}",config.with(|c|c.golden_gate.enzyme)) on:change=move|e|config.update(|c|c.golden_gate.enzyme=match event_target_value(&e).as_str(){"BsmBI"=>collagene_core::cloning::golden_gate::GoldenGateEnzyme::BsmBI,"BbsI"=>collagene_core::cloning::golden_gate::GoldenGateEnzyme::BbsI,"SapI"=>collagene_core::cloning::golden_gate::GoldenGateEnzyme::SapI,_=>collagene_core::cloning::golden_gate::GoldenGateEnzyme::BsaI})><option value="BsaI">"BsaI — 4-base fusion"</option><option value="BsmBI">"BsmBI — 4-base fusion"</option><option value="BbsI">"BbsI — 4-base fusion"</option><option value="SapI">"SapI — 3-base fusion"</option></select></label>
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
                    if !saved.is_current(&config.get_untracked(), source_context.get_untracked(),
                        &target_info.get_untracked(), epoch.get_untracked()) {
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
                                    {crate::cloning::preparation::end_view::pair(Some(&row.left),Some(&row.right))}<p>{format!("Assembly-oriented top DNA: {}",short(&row.oriented_top))}</p>
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
        <PreparedPoolSection molecule sources source_context source_revisions target target_info config epoch on_candidate on_finalize on_source focus_row status/>
        <Show when=move||config.get().mode==Mode::Overlap&&!disable_plain_primer_design><section class="overlap-primer-planner"><h3>"2. Suggest and review PCR primers"</h3><p class="helper">"Use the selected source regions and order to suggest one primer pair per fragment, including added homology tails. Your current primers stay unchanged until you explicitly apply the proposal. The independent target is never generated or selected by this step."</p><details class="assembly-primer-settings"><summary>"Overlap primer design settings"</summary><div class="primer-grid"><label>"Junction homology (bp)"<input type="number" min="1" aria-label="Suggested overlap homology bases" prop:value=move||planner_form.with(|f|f.overlap.clone()) on:input=move|e|planner_form.update(|f|f.overlap=event_target_value(&e))/></label><label>"Minimum annealing length"<input type="number" min="1" max="200" aria-label="Suggested primer minimum length" prop:value=move||planner_form.with(|f|f.minimum.clone()) on:input=move|e|planner_form.update(|f|f.minimum=event_target_value(&e))/></label><label>"Maximum annealing length"<input type="number" min="1" max="200" aria-label="Suggested primer maximum length" prop:value=move||planner_form.with(|f|f.maximum.clone()) on:input=move|e|planner_form.update(|f|f.maximum=event_target_value(&e))/></label><label>"Target annealing Tm (°C)"<input type="number" step="0.1" aria-label="Suggested primer target Tm" prop:value=move||planner_form.with(|f|f.tm.clone()) on:input=move|e|planner_form.update(|f|f.tm=event_target_value(&e))/></label></div><p class="helper">"Default filters: GC 20–80%, pair ΔTm ≤5°C, homopolymer run ≤5 and 3′ complementary run ≤3. Tm assumes 50 mM sodium and 500 nM total strands; these are calculations, not reagent or experimental performance predictions."</p></details><button class="button secondary full" disabled=move||busy.get() on:click=move|_|suggest.run(())>"Suggest PCR primers"</button>
        {move||suggestions.get().map(|proposal|{
            let saved=proposal.clone();
            let guarded=StoredValue::new(proposal.clone());
            let edit=Callback::new(move|index:usize|guarded.with_value(|saved|{
                if !saved.is_current(&config.get_untracked(), &planner_form.get_untracked(),
                    source_context.get_untracked(), &target_info.get_untracked(),
                    &target.get_untracked(), epoch.get_untracked()) {
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
                <button class="button secondary full" on:click=move|_|{if !saved.is_current(&config.get_untracked(), &planner_form.get_untracked(), source_context.get_untracked(), &target_info.get_untracked(), &target.get_untracked(), epoch.get_untracked()){status.set("Sources or settings changed; request fresh suggestions.".into());return;}match serde_json::to_vec_pretty(&serde_json::json!({"sources":saved.context,"target":saved.target_context,"report":saved.report})){Ok(bytes)=>crate::io::download::download_with_completion("overlap-primer-suggestions.json","application/json",&bytes,move|result|{if status.try_get_untracked().is_none(){return;}if let Err(e)=result {status.set(e)}}),Err(e)=>status.set(e.to_string())}}>"Download suggested primer evidence"</button>
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
            let on_download=Callback::new(move|()|{let result=if let Some(preview)=&report.tailing{serde_json::to_vec_pretty(preview)}else if let Some(preview)=&report.ligation{serde_json::to_vec_pretty(preview)}else if let Some(preview)=&report.overlap{serde_json::to_vec_pretty(preview)}else if let Some(profile)=&report.profile{serde_json::to_vec_pretty(profile)}else{status.set("No current native method evidence is available.".into());return;};match result{Ok(bytes)=>crate::io::download::download_with_completion("cloning-method-evidence.json","application/json",&bytes,move|result|{if status.try_get_untracked().is_none(){return;}if let Err(e)=result {status.set(e)}}),Err(e)=>status.set(e.to_string())}});
            let product=proof.ligation.as_ref().map(|p|(&p.compiled.plan,&p.verification)).or_else(||proof.overlap.as_ref().map(|p|(&p.compiled.plan,&p.verification))).or_else(||proof.tailing.as_ref().map(|p|(&p.compiled.plan,&p.verification))).map(|(plan,report)|{
                let captured = StoredValue::new(proof.clone());
                let captured_sources = StoredValue::new(Some(proof.context.clone()));
                let captured_target = StoredValue::new(proof.target_context.clone());
                let guarded_finalize = on_finalize.map(|finalize|Callback::new(move|candidate:Candidate| {
                    let fresh = captured.with_value(|saved| {
                        evidence.with_untracked(|current|current.as_ref().is_some_and(|current|Arc::ptr_eq(current,saved)))
                            && saved.is_current(&config.get_untracked(), source_context.get_untracked(),
                                source_revisions.get_untracked(), &target_info.get_untracked(),
                                &target.get_untracked(), epoch.get_untracked())
                    });
                    if !fresh || source_context.get_untracked()!=captured_sources.get_value() || target_info.get_untracked()!=captured_target.get_value() {
                        return Err("Sources or preparation settings changed. Compute fresh method evidence before creating outputs.".into());
                    }
                    finalize.run(candidate)
                }));
                view!{<crate::assembly::review::assembly_product_view::AssemblyProductView plan=plan.clone() target=molecule.get_untracked() report=report.clone() epoch on_products on_finalize=guarded_finalize/>}.into_any()
            });
            view!{<MethodEvidence proof current active_row on_edit=edit on_source=view_source product on_review=review on_download/>}})}
    </>}
}
// Erase the large evidence display type at the component boundary. The panel retains
// scientific acceptance guards, worker state and all actions that can change outputs.
#[component]
fn MethodEvidence(
    proof: Arc<Evidence>,
    current: StoredValue<Arc<Evidence>>,
    #[prop(optional_no_strip)] active_row: Option<Memo<usize>>,
    on_edit: Callback<usize>,
    #[prop(optional_no_strip)] on_source: Option<Callback<usize>>,
    product: Option<AnyView>,
    on_review: Callback<()>,
    on_download: Callback<()>,
) -> AnyView {
    let diagnostics = proof
        .overlap
        .as_ref()
        .map(collagene_core::cloning::assembly::diagnostics::overlap)
        .or_else(|| {
            proof
                .ligation
                .as_ref()
                .map(collagene_core::cloning::assembly::diagnostics::ligation)
        });
    let names = proof
        .context
        .iter()
        .map(|source| source.name.clone())
        .collect::<Vec<_>>();
    let preparations = proof
        .config
        .rows
        .iter()
        .map(OverlapInput::summary)
        .collect::<Vec<_>>();
    view!{<div class="method-evidence">
            {move||current.with_value(|proof|{
                let physical=proof.ligation.as_ref().map(|p|&p.physical).or_else(||proof.tailing.as_ref().map(|p|&p.physical))?;
                let index=active_row.map(|row|row.get()).unwrap_or(0);
                let label=proof.context.get(index).map(|source|source.name.clone()).unwrap_or_else(||format!("Fragment {}",index+1));
                Some(view!{<section class="active-fragment-end-preview"><h3>{format!("Fragment {} · {} · reconstructed ends",index+1,label)}</h3>
                    {physical.get(index).and_then(Option::as_ref).map(crate::cloning::preparation::end_view::summary).unwrap_or_else(||view!{<p class="physical-end-unavailable">"Physical ends unavailable. Review the blocking preparation result below."</p>}.into_any())}
                    <p class="helper">"Shown in assembly orientation after declared preparation and terminal treatment."</p>
                </section>})
            })}
            {diagnostics.map(|result|match result{Ok(diagnostics)=>view!{<MethodDiagnostics diagnostics names preparations on_edit on_source/>}.into_any(),Err(error)=>view!{<p class="helper">{format!("Structured preparation diagnostics unavailable: {error}. The full fresh findings remain below.")}</p>}.into_any()})}
            {product}

            {proof.tailing.as_ref().map(|preview|view!{<article class="tail-evidence"><h3>"Reconstructed tailing and ligation"</h3><p>{format!("Preparation reconstructed: {} · independent-target assessment: {}",preview.preparation_complete,method_status(&preview.verification.status))}</p><p>{format!("Declared {:?} end-family assessment: {}",preview.recipe.profile,match preview.profile_status{ProfileStatus::Matched=>"matched actual ends",ProfileStatus::Incompatible=>"incompatible actual ends",ProfileStatus::Undetermined=>"undetermined"})}</p>{preview.profile_findings.iter().map(|finding|view!{<p>{finding.clone()}</p>}).collect_view()}<ul>{preview.compiled.rows.iter().enumerate().map(|(index,row)|view!{<li>{format!("Fragment {}: {} · {}",index+1,match row.action{Action::KeepExisting=>"keep existing prepared ends".into(),Action::AddTail{base}=>format!("add {} to both 3′ ends",base.as_str())},row.tail_step.map(|step|format!("tailing step {}",step+1)).unwrap_or_else(||"no tailing operation".into()))}</li>}).collect_view()}</ul><p class="helper">"Source preparation, optional addition, treatment and orientation are reconstructed exactly. Profile matching is separate from ligation chemistry and the independently reviewed intended DNA."</p>{preview.junctions.iter().map(|junction|view!{<details><summary>{format!("Junction {} · fragment {} → {}{}",junction.index+1,junction.upstream+1,junction.downstream+1,if junction.closing{" · circle closure"}else{""})}</summary><div class="physical-junction-ends"><section><h4>"Upstream right end"</h4>{crate::cloning::preparation::end_view::right(junction.upstream_end.as_ref())}</section><section><h4>"Downstream left end"</h4>{crate::cloning::preparation::end_view::left(junction.downstream_end.as_ref())}</section></div></details>}).collect_view()}{preview.verification.findings.iter().map(|finding|view!{<p>{finding.clone()}</p>}).collect_view()}</article>})}

            {proof.ligation.as_ref().map(|preview|view!{<article class="ordered-ligation-evidence"><h3>{if proof.config.mode==Mode::GoldenGate{"Ordered Type IIS restriction/ligation evidence"}else{"Reconstructed ligation evidence"}}</h3><p>{format!("Preparation reconstructed: {} · independent-target assessment: {}",preview.preparation_complete,method_status(&preview.verification.status))}</p><p class="helper">"These ends and junctions were freshly reconstructed from the original source operations. An end-compatible ordered candidate is not a unique mixture product or demonstrated reaction outcome."</p>{preview.prepared.iter().enumerate().map(|(i,fragment)|view!{<details><summary>{format!("Fragment {} · {} bp · {}",i+1,fragment.sequence.len(),proof.config.rows.get(i).map(OverlapInput::summary).unwrap_or_default())}</summary><pre>{format!("Prepared source before orientation: {}",short(fragment.sequence.as_str()))}</pre>{preview.physical.get(i).and_then(Option::as_ref).map(crate::cloning::preparation::end_view::summary)}<p>{preview.physical.get(i).and_then(Option::as_ref).map(method_duplex).unwrap_or_else(||"Physical end geometry unavailable; no chemistry inferred".into())}</p><p class="helper">"End summary is in assembly orientation, after any explicit terminal treatment."</p></details>}).collect_view()}{preview.junctions.iter().map(|junction|view!{<details><summary>{format!("Junction {} · fragment {} → {}{}",junction.index+1,junction.upstream+1,junction.downstream+1,if junction.closing{" · circle closure"}else{""})}</summary><div class="physical-junction-ends"><section><h4>"Upstream fragment: right end"</h4>{crate::cloning::preparation::end_view::right(junction.upstream_end.as_ref())}</section><section><h4>"Downstream fragment: left end"</h4>{crate::cloning::preparation::end_view::left(junction.downstream_end.as_ref())}</section></div><p>{format!("Upstream right: {}",method_end(junction.upstream_end.as_ref()))}</p><p>{format!("Downstream left: {}",method_end(junction.downstream_end.as_ref()))}</p></details>}).collect_view()}{preview.product_physical.as_ref().map(|physical|view!{<section><h4>"Reconstructed product"</h4><p>{method_duplex(physical)}</p>{physical.nicks.iter().map(|nick|view!{<p>{format!("Retained {} strand nick at top-reference boundary {}",match nick.strand{collagene_core::chemistry::Strand::Top=>"top",collagene_core::chemistry::Strand::Bottom=>"bottom"},nick.position)}</p>}).collect_view()}</section>})}{preview.verification.findings.iter().map(|finding|view!{<p>{finding.clone()}</p>}).collect_view()}</article>})}
            {proof.overlap.as_ref().map(|preview|view!{<article class="overlap-preparation-evidence"><h3>"Prepared fragments from original sources"</h3>{crate::assembly::review::assembly_product_view::processing_views(&preview.verification)}<p>{format!("Preparation reconstructed: {} · independent-target assessment: {}",preview.preparation_complete,method_status(&preview.verification.status))}</p><p class="helper">"Prepared reference sequences below are reconstructed from the declared originals. They are not substituted as plan inputs, confirmed amplicons or physical reaction outcomes."</p>{preview.prepared.iter().enumerate().map(|(i,fragment)|view!{<details><summary>{format!("Fragment {} · {} bp · {}",i+1,fragment.sequence.len(),proof.config.rows.get(i).map(OverlapInput::summary).unwrap_or_default())}</summary><pre>{format!("Prepared source 5′→3′ (before orientation): {}",short(fragment.sequence.as_str()))}</pre></details>}).collect_view()}{preview.verification.findings.iter().map(|finding|view!{<p>{finding.clone()}</p>}).collect_view()}</article>})}
            {proof.profile.as_ref().map(|p|view!{<article class="overlap-profile-evidence"><h3>"Overlap evidence"</h3><p>{format!("Geometry consistent: {} · analysis complete: {} · published length guidance: {}",p.geometry_valid,p.complete,match p.recommendation_compliant{Some(true)=>"met",Some(false)=>"outside published guidance",None=>"not established"})}</p><p class="helper">{format!("{} comparisons; no exhaustive alternative-product search",p.comparisons)}</p>{p.junctions.iter().map(|j|view!{<details class="profile-junction"><summary>{format!("Junction {} · {}bp · {}",j.index+1,j.overlap_bases,if j.exact_match{"exact end homology"}else{"incompatible homology"})}</summary><pre>{format!("Upstream {}\nDownstream {}",j.homology,j.downstream_homology)}</pre><p>{format!("GC: {} · advisory length: {} · repeat search complete: {}",j.gc_fraction.map(|f|format!("{:.1}%",f*100.)).unwrap_or_else(||"unavailable".into()),guidance(j.recommended_min,j.recommended_max),j.repeat_search_complete)}</p>{j.occurrences.iter().take(40).map(|o|view!{<p>{format!("Source {} · {}–{}{}{}",o.fragment+1,o.start+1,o.end,if o.reverse{" reverse"}else{" forward"},if o.assigned_junction{" · assigned end"}else{" · additional occurrence"})}</p>}).collect_view()}<p class="helper">"Up to 40 occurrence positions shown; full evidence download retains every reported position."</p></details>}).collect_view()}<details><summary>"Profile findings and assumptions"</summary>{p.findings.iter().chain(&p.warnings).map(|w|view!{<p>{w.clone()}</p>}).collect_view()}</details></article>})}
            <button class="button primary full" on:click=move|_|on_review.run(())>"Review method candidate"</button><button class="button secondary full" on:click=move|_|on_download.run(())>"Download method evidence"</button>
        </div>}.into_any()
}

impl PlannerForm {
    fn request(
        &self,
        state: &MethodsState,
        inputs: Vec<Molecule>,
        target: &Molecule,
    ) -> Result<
        (
            Vec<PlanInput>,
            collagene_core::cloning::overlap::primer_design::Request,
        ),
        String,
    > {
        let form = self;
        use collagene_core::cloning::overlap::primer_design::{Fragment, Request};
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
            primer_options: collagene_core::cloning::pcr::primers::DesignOptions {
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
        collagene_core::cloning::overlap::primer_design::preflight(&borrowed, target, &request)
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
}

fn golden_gate_request(
    state: &MethodsState,
    inputs: Vec<Molecule>,
    target: &Molecule,
) -> Result<
    (
        Vec<PlanInput>,
        collagene_core::cloning::golden_gate::primers::Request,
    ),
    String,
> {
    use collagene_core::cloning::golden_gate::primers::{Fragment, Request};
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
    collagene_core::cloning::golden_gate::primers::source_preflight(&borrowed, target, &request)
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
