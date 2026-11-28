//! Finite Golden Gate pool capture, guarded evidence and explicit candidate selection.
use crate::cloning::preparation::{
    cloning_methods_model::*,
    cloning_preparations::{ligation_request, prepared_rows},
    diagnostics::{MethodDiagnostics, method_status, short},
    golden_gate_pool_form::PoolForm,
    review::TargetContext,
};
use crate::{io::download::prepare_download, platform::worker::CachedWorker};
use collagene_core::{
    chemistry::EndTreatment,
    documents::transfer::Candidate,
    planning::plan::PlanInput,
    protocol::{
        request::Request,
        schema::{AnalysisJob, AnalysisResult},
    },
};
use collagene_model::Molecule;
use leptos::prelude::*;
use std::{
    rc::Rc,
    sync::{
        Arc,
        atomic::{AtomicBool, Ordering},
    },
};
use wasm_bindgen_futures::spawn_local;

#[component]
pub(super) fn PreparedPoolSection(
    molecule: Memo<Molecule>,
    sources: Memo<Option<Vec<Molecule>>>,
    source_context: Memo<Option<Vec<crate::sources::identity::SourceInfo>>>,
    source_revisions: Memo<Option<Vec<String>>>,
    target: Memo<String>,
    target_info: Memo<TargetContext>,
    config: RwSignal<MethodsState>,
    epoch: RwSignal<u64>,
    on_candidate: Callback<Candidate>,
    #[prop(optional_no_strip)] on_finalize: Option<Callback<Candidate, Result<(), String>>>,
    #[prop(optional_no_strip)] on_source: Option<Callback<(usize, usize), Result<(), String>>>,
    focus_row: Callback<usize>,
    status: RwSignal<String>,
) -> impl IntoView {
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
        saved.admit(
            &config.with_untracked(MethodsState::preparation_context),
            source_context.get_untracked(),
            source_revisions.get_untracked(),
            &target_info.get_untracked(),
            &target.get_untracked(),
            epoch.get_untracked(),
        )?;
        Ok(())
    });
    let reaction_capture = Callback::new(move |saved: Arc<PoolContext>| {
        reaction_guard.run(saved.clone())?;
        let originals = sources
            .get_untracked()
            .ok_or("Choose available original reaction sources.")?;
        let inputs = saved.bind_inputs(originals)?;
        Ok::<_, String>((inputs, molecule.get_untracked()))
    });
    move || {
        let snapshot = reaction_context.get()?;
        let captured = StoredValue::new(snapshot.clone());
        let edit =
            Callback::new(
                move |row: usize| match reaction_guard.run(captured.get_value()) {
                    Ok(()) => focus_row.run(row),
                    Err(error) => status.set(error),
                },
            );
        let view_source = on_source.map(|navigate| {
            Callback::new(move |row: usize| {
                let result = captured.with_value(|ctx| {
                    reaction_guard.run(ctx.clone())?;
                    let input = ctx.config.rows.get(row).ok_or("Source role unavailable.")?;
                    navigate.run((input.document, row))
                });
                if let Err(error) = result {
                    status.set(error);
                }
            })
        });
        Some(
            view! {<PreparedPool proof=snapshot config guard=reaction_guard capture=reaction_capture epoch on_candidate on_finalize on_edit=edit on_source=view_source/>},
        )
    }
}

#[derive(Clone, PartialEq, Eq)]
struct PoolContext {
    config: MethodsState,
    context: Vec<crate::sources::identity::SourceInfo>,
    target_context: TargetContext,
    target: String,
    revisions: Vec<String>,
    epoch: u64,
}
impl PoolContext {
    fn admit(
        &self,
        config: &MethodsState,
        context: Option<Vec<crate::sources::identity::SourceInfo>>,
        revisions: Option<Vec<String>>,
        target_context: &TargetContext,
        target: &str,
        epoch: u64,
    ) -> Result<(), String> {
        if config != &self.config
            || context != Some(self.context.clone())
            || revisions != Some(self.revisions.clone())
            || target_context != &self.target_context
            || target != self.target
            || epoch != self.epoch
        {
            return Err("Reaction source events, target or settings changed.".into());
        }
        Ok(())
    }
    fn bind_inputs(&self, originals: Vec<Molecule>) -> Result<Vec<PlanInput>, String> {
        let roles = self.config.role_ids(originals.len())?;
        let inputs = roles
            .into_iter()
            .zip(originals)
            .map(|(role, molecule)| PlanInput {
                id: role.input_id,
                revision: molecule.revision(),
                molecule,
            })
            .collect();
        Ok(inputs)
    }
}
#[derive(Clone)]
enum PoolReceipt {
    Prepared(Box<collagene_core::cloning::golden_gate::prepared::Report>),
    Reaction(Box<collagene_core::cloning::golden_gate::reaction::Report>),
}
impl PoolReceipt {
    fn pool(&self) -> Option<&collagene_core::cloning::golden_gate::GoldenGateReport> {
        match self {
            Self::Prepared(r) => r.pool.as_ref(),
            Self::Reaction(r) => r.pool.as_ref(),
        }
    }
    fn assessments(&self) -> &[collagene_core::cloning::golden_gate::prepared::Assessment] {
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
    fn verification(&self) -> &collagene_core::verification::evidence::VerificationReport {
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
                    collagene_core::planning::plan::Source::Step(i)
                        if i >= r.verification.products.len() =>
                    {
                        Some(m.index)
                    }
                    _ => None,
                })
                .collect(),
        }
    }
    fn diagnostics(
        &self,
    ) -> Option<
        Result<collagene_core::cloning::assembly::diagnostics::Report, collagene_model::Error>,
    > {
        match self {
            Self::Prepared(r) => Some(collagene_core::cloning::assembly::diagnostics::ligation(
                &r.preview,
            )),
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
        s: &collagene_core::cloning::golden_gate::prepared::Selection,
    ) -> Result<String, String> {
        match self {
            Self::Prepared(r) => {
                collagene_core::cloning::golden_gate::prepared::selection_revision(r, s)
            }
            Self::Reaction(r) => {
                collagene_core::cloning::golden_gate::reaction::selection_revision(r, s)
            }
        }
        .map_err(|e| e.to_string())
    }
    fn select_job(
        &self,
        s: collagene_core::cloning::golden_gate::prepared::Selection,
    ) -> AnalysisJob {
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
    use collagene_core::cloning::golden_gate::prepared as pool;
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
                                run_origin: row.run_origin.clone(),
                                source_generation: row.source_generation.clone(),
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
                            1 => collagene_core::cloning::golden_gate::reaction::Substrate::WholeSource {
                                origin: config.rows[i]
                                    .run_origin
                                    .as_ref()
                                    .map(|p| p.as_ref().clone()),
                            },
                            2 => collagene_core::cloning::golden_gate::reaction::Substrate::Prepared {
                                preparation: f.preparation,
                                treatment: None,
                            },
                            _ => collagene_core::cloning::golden_gate::reaction::Substrate::Prepared {
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
                        collagene_core::cloning::golden_gate::reaction::Row {
                            generation: f.generation,
                            role: f.role,
                            input_id: f.input_id,
                            substrate,
                        }
                    })
                    .collect();
                let request = collagene_core::cloning::golden_gate::reaction::Request {
                    name: "Declared Golden Gate reaction".into(),
                    rows,
                    options: pool_options,
                };
                let revision = collagene_core::cloning::golden_gate::reaction::inspect_revision(
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
    view! {<section class="prepared-gate-pool" aria-label="Prepared Golden Gate pool alternatives"><h3>"Finite Golden Gate reaction and recutting"</h3>
        <p class="helper">"Inspect selected prepared fragments, or explicitly choose what enters the reaction for every role. The finite search reports competing products, recutting and uncut circular carryover without selecting a winner or predicting experimental yield."</p>
        <p class="helper">"Reaction choices and search limits are kept with the assembly draft. Reopen the draft and inspect again to reconstruct current evidence."</p>
        {construct_matrix.then(||view!{<p class="helper">"These are this construct's pool settings. Use the Construct reaction candidate panel below for exact source-bound discovery, selection, product/junction review and output creation. Ordered method evidence is unavailable until you explicitly return to ordered prepared assembly."</p>})}
        <label class="checkbox-label"><input type="checkbox" aria-label="Inspect whole Golden Gate reaction" prop:checked=move||form.with(|f|f.reaction) on:change=move|e|form.update(|f|f.reaction=event_target_checked(&e))/>"Declare reaction substrates (whole DNA / whole amplicons / prepared fragments)"</label>
        <Show when=move||form.with(|f|f.reaction)><p class="helper">"Choose what enters the tube for every role. Whole source includes its circular backbone and all cut pieces. Whole amplicon uses the current primers before digestion and omits post-digest repair/treatment. Retaining a prepared fragment keeps its current explicit operations. Primers and the intended target are unchanged."</p>
        {proof.with_value(|p|p.config.rows.iter().enumerate().map(|(i,row)|{let pcr=row.has_pcr_amplicon();view!{<label>{format!("{}: {}",i+1,p.context.get(i).map(|c|c.name.as_str()).unwrap_or("source"))}<select aria-label=format!("Reaction substrate {}",i+1) prop:value=move||form.with(|f|f.choices.get(i).copied().unwrap_or(0).to_string()) on:change=move|e|form.update(|f|if let Some(choice)=f.choices.get_mut(i){*choice=event_target_value(&e).parse().unwrap_or(0)})><option value="0">"Current selected prepared fragment"</option><option value="1">"Whole original source"</option><option value="2" disabled=!pcr>"Complete PCR amplicon before digestion"</option></select></label>}}).collect_view())}</Show>
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
            <p class="helper">"Positional cleavage uses only fresh inherited or explicitly declared source state. Unknown or impaired sites withhold an undetermined input pool; blocked sites remain uncut. Candidate recut uncertainty is not a site-free or isolation claim. Standard newly synthesized PCR DNA is unmodified under this model."</p>
            {bound.report.pool().map(|pool|{let determined=pool.products_determined;view!{<p>{format!("Initial pool products determined {} · {} shared cleavage comparisons · pinned profile {}",pool.products_determined,pool.comparisons,pool.profile_revision)}</p><Show when=move||!determined><p class="helper">"The requested pool is withheld. No partial subset of available sources is enumerated as the whole reaction."</p></Show>{pool.source_cleavage.iter().enumerate().map(|(i,e)|view!{<CheckedCleavageDetails evidence=Arc::new(e.clone()) label=format!("Source {} positional cleavage",i+1)/>}).collect_view()}}})}

            {match diagnostics{Some(Ok(diagnostics))=>view!{<details><summary>"Recomputed row preparation diagnostics"</summary><MethodDiagnostics diagnostics names preparations on_edit on_source/></details>}.into_any(),Some(Err(error))=>view!{<p class="helper">{format!("Structured diagnostics unavailable: {error}")}</p>}.into_any(),None=>view!{<details><summary>"Native preparation and probe findings"</summary>{bound.report.verification().findings.iter().map(|finding|view!{<p>{finding.clone()}</p>}).collect_view()}</details>}.into_any()}}
            {bound.report.pool().map(|pool|{let count=pool.candidates.len();let pages=count.max(1).div_ceil(50);let start=page.get().min(pages-1)*50;view!{<h4>{format!("{} directed candidates · {} physical sequence groups",count,pool.stable_physical_product_count)}</h4><p>{if pool.enumeration_complete{"Declared finite substrate-pool search completed within the selected options."}else{"Search incomplete: additional declared substrate-pool alternatives may exist."}}</p><p class="helper">{format!("States examined {}. Repeated-copy/real-mixture exhaustiveness is outside this model. Uncut circular substrate carryover roles: {} (separate from assembled candidates). Sequence grouping is not a unique source/annotation history.",pool.states_examined,if pool.carryover_sources.is_empty(){"none".into()}else{pool.carryover_sources.iter().map(|row|(row+1).to_string()).collect::<Vec<_>>().join(", ")})}</p><ul>{pool.truncation_reasons.iter().chain(pool.warnings.iter()).map(|warning|view!{<li>{warning.clone()}</li>}).collect_view()}</ul>
                {pool.candidates.iter().enumerate().skip(start).take(50).map(|(index,candidate)|{let assessment=bound.report.assessments().iter().find(|a|a.candidate_index==index);let allowed=candidate.restriction_site_free;view!{<article class="prepared-gate-candidate"><h4>{format!("Candidate {} · {} top-strand bases · {} pool fragments",index+1,candidate.duplex.top.len(),candidate.path.len())}</h4><p>{if allowed{"Recognition-site-free under the chosen profile"}else{"Recognition retained/recreated — diagnosis only"}}</p><p>{format!("Independent target match: {} · used recipe rows (unique): {}",assessment.map(|a|a.target_matches.to_string()).unwrap_or_else(||"unavailable".into()),assessment.map(|a|a.source_roles.iter().map(|&row|proof.with_value(|p|format!("{} ({})",row+1,p.context.get(row).map(|s|s.name.as_str()).unwrap_or("unknown source")))).collect::<Vec<_>>().join(", ")).unwrap_or_else(||"unavailable".into()))}</p><p class="helper">{format!("Raw path: {}",candidate.path.iter().map(|step|pool.fragments.get(step.fragment).map(|fragment|format!("pool fragment {} · row {} / digest {}{}",step.fragment+1,fragment.source_index+1,fragment.fragment_index+1,if step.reverse{" reverse"}else{" forward"})).unwrap_or_else(||"unavailable pool fragment".into())).collect::<Vec<_>>().join(" → "))}</p><pre>{short(candidate.duplex.top.as_str())}</pre><p>{format!("{} recut/recognition observations",candidate.recut_sites.len())}</p><PoolCandidateDetails bound=bound.clone() index/>
                    <button class="button secondary" disabled=move||pending.get()||!allowed on:click=move|_|select.run((request.get_value(),index,false))>"Review pool candidate"</button>{on_finalize.map(|_|view!{<button class="button primary" disabled=move||pending.get()||!allowed on:click=move|_|select.run((request.get_value(),index,true))>"Create outputs from this candidate"</button>})}
                </article>}}).collect_view()}
                <p>{format!("Page {} of {} · candidates {}–{} of {}",page.get().min(pages-1)+1,pages,if count==0{0}else{start+1},(start+50).min(count),count)}</p><button class="text-tool" disabled=move||page.get()==0 on:click=move|_|page.update(|p|*p=p.saturating_sub(1))>"Previous pool candidates"</button><button class="text-tool" disabled=move||{page.get()+1>=pages} on:click=move|_|page.update(|p|*p=(*p+1).min(pages-1))>"Next pool candidates"</button>
            }})}
            <button class="button secondary" on:click=move|_|{let saved=request.get_value();let result=pool_guard.run(saved.clone()).and_then(|()|saved.report.json()).and_then(|bytes|prepare_download("prepared-golden-gate-pool.json","application/json",&bytes));crate::io::download::complete_download_guarded(result,move||pool_guard.run(saved.clone()),move|result|{if status.try_get_untracked().is_some(){if let Err(error)=result{status.set(error);}}});}>"Download full prepared-pool evidence"</button>
        </section>}})}
    </section>}
}

fn cleavage_label(outcome: collagene_core::chemistry::methylation::Outcome) -> &'static str {
    use collagene_core::chemistry::methylation::Outcome;
    match outcome {
        Outcome::PermittedUnderProfile => "permitted under the pinned profile",
        Outcome::Blocked => "blocked",
        Outcome::Impaired => "impaired — cleavage unresolved",
        Outcome::Unknown => "unknown — cleavage unresolved",
    }
}
#[component]
fn CheckedCleavageDetails(
    evidence: Arc<collagene_core::chemistry::tracked_duplex::CheckedCleavageEvidence>,
    label: String,
) -> impl IntoView {
    let page = RwSignal::new(0usize);
    let total = evidence.sites.len();
    let pages = total.max(1).div_ceil(25);
    let saved = StoredValue::new(evidence);
    view! {<details><summary>{label}</summary>
        {saved.with_value(|e|view!{<p>{format!("Source {} · scan complete {} · products determined {} · {} comparisons",e.source_revision,e.complete,e.products_determined,e.comparisons)}</p><p>{format!("Pinned positional profile {} · {} · catalog {}",e.profile_version,e.profile_revision,e.options.selection.catalog_content_hash)}</p><ul>{e.warnings.iter().chain(&e.limits).map(|s|view!{<li>{s.clone()}</li>}).collect_view()}</ul>})}
        {move||saved.with_value(|e|e.sites.iter().skip(page.get()*25).take(25).map(|site|view!{<p>{format!("{} at source base {}{} · {} · top/bottom cuts {:?}/{:?} · {}",site.enzyme,site.recognition_start+1,if site.reverse{" reverse"}else{" forward"},cleavage_label(site.outcome),site.top_cut,site.bottom_cut,site.reason)}</p>}).collect_view())}
        <p>{move||format!("{} recognition observations · page {} of {}. Full positional evidence is retained in JSON.",total,page.get()+1,pages)}</p>
        <button class="text-tool" disabled={move||page.get()==0} on:click=move|_|page.update(|p|*p=p.saturating_sub(1))>"Previous cleavage observations"</button><button class="text-tool" disabled={move||page.get()+1>=pages} on:click=move|_|page.update(|p|*p=(*p+1).min(pages-1))>"Next cleavage observations"</button>
    </details>}
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
        {bound.with_value(|saved|saved.report.pool().and_then(|pool|pool.candidates.get(index)).map(|candidate|view!{<p>{format!("Directed key {} · physical key {}",candidate.directed_sequence_key,candidate.physical_sequence_key)}</p><CheckedCleavageDetails evidence=Arc::new(candidate.cleavage.clone()) label="Candidate positional recut assessment".to_string()/>}))}
        {move||bound.with_value(|saved|saved.report.pool().and_then(|pool|pool.candidates.get(index)).map(|candidate|candidate.recut_sites.iter().skip(page.get()*25).take(25).map(|site|view!{<p>{format!("Recognition at {}{} · paired {} · {} · definite cleavage under model {}",site.recognition_start+1,if site.reverse{" reverse"}else{" forward"},site.fully_double_stranded,cleavage_label(site.outcome),site.cleavable_under_model)}</p>}).collect_view()))}
        <p>{move||format!("{} recut observations · page {} of {}",total,page.get()+1,pages)}</p>
        <button class="text-tool" disabled=move||page.get()==0 on:click=move|_|page.update(|p|*p=p.saturating_sub(1))>"Previous recut observations"</button><button class="text-tool" disabled=move||{page.get()+1>=pages} on:click=move|_|page.update(|p|*p=(*p+1).min(pages-1))>"Next recut observations"</button>
    </details>}
}
