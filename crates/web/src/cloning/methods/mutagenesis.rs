//! Revision-bound mutation design and conditional mismatch-PCR reference previews.
use crate::platform::worker;
use crate::platform::worker::CachedWorker;
use collagene_core::{
    cloning::pcr::{
        mutagenesis::{InversePcrDesign, InversePcrOptions, PcrOptions},
        primers::BindingOptions,
    },
    documents::transfer::Candidate,
    material::lineage::{EventBinding, InputEventBinding},
    protocol::{
        request::Request,
        response::Response,
        schema::{AnalysisJob, AnalysisResult},
    },
    sequence::analysis::Primer,
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
fn number<T: std::str::FromStr>(value: String, label: &str) -> Result<T, String> {
    value
        .trim()
        .parse()
        .map_err(|_| format!("{label} is not a valid number"))
}
#[derive(Clone, PartialEq, Eq)]
pub struct MutationSelection {
    pub source: EventBinding,
    pub document: usize,
    pub start: usize,
    pub end: usize,
}
#[derive(Clone)]
pub struct IntentOpen {
    pub source: EventBinding,
    pub report: collagene_core::cloning::pcr::mutation_intent::Report,
}
#[derive(Clone)]
pub struct CandidateOpen {
    pub candidate: Candidate,
    pub sources: Vec<InputEventBinding>,
    pub target: EventBinding,
    pub epoch: u64,
}
type InverseForm = (String, String, String, String, String, String);
struct BoundInverse {
    report: InversePcrDesign,
    origin: Option<crate::history::run_origin_model::PinnedOrigin>,
    source: EventBinding,
    target: EventBinding,
    epoch: u64,
    generation: u32,
    form: InverseForm,
    draft: crate::primers::primer_pair_model::PairDraft,
}
struct IntendedEdit {
    report: collagene_core::cloning::pcr::mutation_intent::Report,
    source: EventBinding,
    document: usize,
    form: (String, String, String, String),
    epoch: u64,
}
fn dna_preview(dna: &str) -> String {
    let shown = dna.chars().take(200).collect::<String>();
    if dna.len() > shown.len() {
        format!("{shown}… ({} bases total)", dna.len())
    } else {
        shown
    }
}
#[derive(Clone)]
struct AppliedPair {
    primers: (Primer, Primer),
    draft: crate::primers::primer_pair_model::PairDraft,
    source: EventBinding,
    target: EventBinding,
    epoch: u64,
}
pub struct OutputOpen {
    pub source: EventBinding,
    pub target: EventBinding,
    pub epoch: u64,
    pub revision: String,
    pub options: collagene_core::cloning::pcr::outputs::Options,
    pub prepared: collagene_core::cloning::pcr::mismatch_outputs::Prepared,
}
struct PcrEvidence {
    review: Arc<crate::primers::pcr::mismatch_pcr_outputs::Review>,
    generation: u32,
    form: (String, String),
    draft: crate::primers::primer_pair_model::PairDraft,
    inverse_form: InverseForm,
}
#[component]
pub fn MutagenesisPanel(
    molecule: Memo<Molecule>,
    query_options: Memo<Vec<(usize, String)>>,
    get_query: Callback<usize, Option<Molecule>>,
    on_candidate: Callback<CandidateOpen>,
    on_outputs: Callback<OutputOpen, Result<(), String>>,
    get_workspace: Callback<(), Result<collagene_core::documents::workspace::Workspace, String>>,
    on_material: Callback<
        crate::restriction::digest_product::MaterialOpenRequest,
        Result<(), String>,
    >,
    pair_draft: RwSignal<crate::primers::primer_pair_model::PairDraft>,
    pair_selection: Callback<(), Result<crate::primers::primer_pair_model::SelectionSeed, String>>,
    anchor_available: Callback<EventBinding, bool>,
    resolve_source: Callback<EventBinding, Result<Molecule, String>>,
    primer_sources: Memo<Vec<crate::sources::identity::SourceInfo>>,
    list_primers: Callback<EventBinding, Result<collagene_core::cloning::pcr::reuse::List, String>>,
    resolve_primer: Callback<
        crate::primers::primer_picker::ResolveRequest,
        Result<collagene_core::cloning::pcr::reuse::Resolved, String>,
    >,
    epoch: RwSignal<u64>,
    source_info: Callback<usize, Option<crate::sources::identity::SourceInfo>>,
    target_binding: Memo<EventBinding>,
    mutation_selection: Callback<(), Result<MutationSelection, String>>,
    mutation_seed: RwSignal<Option<MutationSelection>>,
    on_intent: Callback<IntentOpen, Result<(), String>>,
    #[prop(optional_no_strip)] on_overlay: Option<
        Callback<crate::primers::pcr::pcr_overlay::Publication>,
    >,
    #[prop(optional_no_strip)] on_navigate: Option<
        Callback<(EventBinding, u64, Vec<(usize, usize)>), Result<(), String>>,
    >,
) -> impl IntoView {
    let on_overlay = on_overlay.unwrap_or_else(|| Callback::new(|_| {}));
    let on_navigate = on_navigate.unwrap_or_else(|| {
        Callback::new(|_| Err("PCR source navigation is unavailable in this context.".into()))
    });
    let template_index = RwSignal::new(0usize);
    let source = Memo::new(move |_| get_query.run(template_index.get()));
    let inverse_origin = RwSignal::new(None::<crate::history::run_origin_model::PinnedOrigin>);
    let applied_pair = RwSignal::new(None::<AppliedPair>);
    let source_binding = Memo::new(move |_| {
        let _ = epoch.get();
        source_info
            .run(template_index.get())
            .map(|info| info.binding)
    });
    let attest_inverse = Callback::new(
        move |(event, origin): (
            EventBinding,
            Option<crate::history::run_origin_model::PinnedOrigin>,
        )|
              -> Result<(), String> {
            let frame = get_workspace.run(())?;
            let retained = crate::protocol::plan_editing_sources::workspace_target(
                &frame,
                &event,
                &event.revision,
            )?;
            let visible = resolve_source.run(event.clone())?;
            if retained != &visible {
                return Err("The exact inverse-PCR source snapshot changed.".into());
            }
            if let Some(origin) = origin {
                origin
                    .as_ref()
                    .validate_for_workspace(&frame, &event)
                    .map_err(|e| e.to_string())?;
            }
            Ok(())
        },
    );
    let mismatches = RwSignal::new("1".to_string());
    let exact = RwSignal::new("1".to_string());
    let edit_start = RwSignal::new("1".to_string());
    let edit_end = RwSignal::new("1".to_string());
    let inserted = RwSignal::new(String::new());
    let intended_name = RwSignal::new("Intended mutation".to_string());
    let intended = RwSignal::new(None::<Arc<IntendedEdit>>);
    let intent_busy = RwSignal::new(false);
    let intent_status = RwSignal::new(String::new());
    let intent_generation = RwSignal::new(0u32);
    let intent_worker = StoredValue::new_local(None::<Rc<CachedWorker>>);
    let intent_form = Callback::new(move |()| {
        (
            intended_name.get_untracked(),
            edit_start.get_untracked(),
            edit_end.get_untracked(),
            inserted.get_untracked(),
        )
    });
    let f_length = RwSignal::new("20".to_string());
    let r_length = RwSignal::new("20".to_string());
    let split = RwSignal::new("0".to_string());
    let pcr = RwSignal::new(None::<Arc<PcrEvidence>>);
    let inverse = RwSignal::new(None::<Arc<BoundInverse>>);
    let inverse_form = Callback::new(move |()| {
        (
            edit_start.get_untracked(),
            edit_end.get_untracked(),
            inserted.get_untracked(),
            f_length.get_untracked(),
            r_length.get_untracked(),
            split.get_untracked(),
        )
    });
    let status = RwSignal::new(String::new());
    let busy = RwSignal::new(false);
    let generation = RwSignal::new(0u32);
    let alive = Arc::new(AtomicBool::new(true));
    {
        let alive = alive.clone();
        on_cleanup(move || alive.store(false, Ordering::Relaxed));
    }
    on_cleanup(move || {
        if let Some(worker) = intent_worker.get_value() {
            worker.close();
        }
    });
    Effect::new(move |_| {
        let _ = (
            intended_name.get(),
            edit_start.get(),
            edit_end.get(),
            inserted.get(),
            source_binding.get(),
            inverse_origin.get(),
            template_index.get(),
            epoch.get(),
        );
        intent_generation.update(|value| *value = value.wrapping_add(1));
        if let Some(worker) = intent_worker.get_value() {
            worker.close();
        }
        intent_worker.set_value(None);
        intended.set(None);
        intent_busy.set(false);
        intent_status.set(String::new());
    });
    let use_selected = Callback::new(move |selected: MutationSelection| {
        let captured_epoch = epoch.get_untracked();
        let result = (|| -> Result<_, String> {
            if selected.start > selected.end
                || epoch.get_untracked() != captured_epoch
                || source_info
                    .run(selected.document)
                    .is_none_or(|info| info.binding != selected.source)
                || !anchor_available.run(selected.source.clone())
            {
                return Err("The selection source changed or its range is unavailable. Select the original source again.".into());
            }
            let start = selected
                .start
                .checked_add(1)
                .ok_or("Selection coordinate overflow.")?;
            let end = selected
                .end
                .checked_add(1)
                .ok_or("Selection coordinate overflow.")?;
            Ok((selected.document, start, end))
        })();
        match result {
            Ok((document, start, end)) => {
                template_index.set(document);
                edit_start.set(start.to_string());
                edit_end.set(end.to_string());
                intent_status.set("Copied the selected source and edit boundaries. Replacement DNA was preserved; preview the intended edit next.".into());
            }
            Err(error) => intent_status.set(error),
        }
    });
    let use_mutation_selection = move |_| match mutation_selection.run(()) {
        Ok(selected) => use_selected.run(selected),
        Err(error) => intent_status.set(error),
    };
    Effect::new(move |_| {
        if let Some(selected) = mutation_seed.get() {
            mutation_seed.set(None);
            use_selected.run(selected);
        }
    });
    let intent_alive = alive.clone();
    let preview_intent = move |_| {
        let captured_epoch = epoch.get_untracked();
        let form = intent_form.run(());
        let prepared = (|| -> Result<_, String> {
            let source_binding = source_binding
                .get_untracked()
                .ok_or("Choose an original source document.")?;
            if !anchor_available.run(source_binding.clone()) {
                return Err("The original source event is unavailable.".into());
            }
            let source = source
                .get_untracked()
                .ok_or("The original source is unavailable.")?;
            let request = collagene_core::cloning::pcr::mutation_intent::Request {
                name: form.0.clone(),
                edit_start: number::<usize>(form.1.clone(), "Edit start")?
                    .checked_sub(1)
                    .ok_or("Edit start is one-based.")?,
                edit_end: number::<usize>(form.2.clone(), "Edit end")?
                    .checked_sub(1)
                    .ok_or("Edit end is one-based.")?,
                inserted: form.3.clone(),
            };
            collagene_core::cloning::pcr::mutation_intent::source_preflight(&source, &request)
                .map_err(|e| e.to_string())?;
            let revision =
                collagene_core::cloning::pcr::mutation_intent::intent_revision(&source, &request)
                    .map_err(|e| e.to_string())?;
            if source.revision() != source_binding.revision {
                return Err("Original source revision changed before preview.".into());
            }
            let worker = Rc::new(CachedWorker::new()?);
            Ok((source, source_binding, request, revision, worker))
        })();
        let (source_molecule, binding, request, revision, w) = match prepared {
            Ok(value) => value,
            Err(error) => {
                intent_status.set(error);
                return;
            }
        };
        if let Some(old) = intent_worker.get_value() {
            old.close();
        }
        intent_worker.set_value(Some(w.clone()));
        let token = intent_generation.get_untracked().wrapping_add(1);
        intent_generation.set(token);
        let document = template_index.get_untracked();
        intended.set(None);
        intent_busy.set(true);
        intent_status
            .set("Building the declared intended reference from the original source…".into());
        let alive = intent_alive.clone();
        spawn_local(async move {
            let response = w
                .ask(Request {
                    id: token,
                    revision: revision.clone(),
                    job: AnalysisJob::BuildMutationIntent {
                        molecule: source_molecule,
                        request,
                    },
                })
                .await;
            w.close();
            if !alive.load(Ordering::Relaxed)
                || intent_generation.get_untracked() != token
                || epoch.get_untracked() != captured_epoch
                || template_index.get_untracked() != document
                || source_binding.get_untracked().as_ref() != Some(&binding)
                || intent_form.run(()) != form
                || !anchor_available.run(binding.clone())
            {
                return;
            }
            intent_worker.set_value(None);
            intent_busy.set(false);
            match response {
                Ok(response) if response.id == token && response.revision == revision => {
                    match response.result {
                        AnalysisResult::MutationIntent(report) => {
                            match report.validate_binding(&revision) {
                                Ok(()) if report.source_revision == binding.revision => {
                                    intended.set(Some(Arc::new(IntendedEdit {
                                        report,
                                        source: binding,
                                        document,
                                        form,
                                        epoch: captured_epoch,
                                    })));
                                    intent_status.set("Intended edit ready for review. Open it explicitly as the independent target before inverse-PCR design.".into());
                                }
                                Ok(()) => intent_status.set(
                                    "Intended edit belongs to a different original source.".into(),
                                ),
                                Err(error) => intent_status.set(error.to_string()),
                            }
                        }
                        AnalysisResult::Error(error) => intent_status.set(error),
                        _ => intent_status.set("Unexpected intended-edit result.".into()),
                    }
                }
                Err(error) => intent_status.set(error),
                _ => intent_status.set("Intended-edit response identity changed.".into()),
            }
        });
    };
    let open_intent = move |_| {
        let Some(bound) = intended.get_untracked() else {
            return;
        };
        if epoch.get_untracked() != bound.epoch
            || template_index.get_untracked() != bound.document
            || source_binding.get_untracked().as_ref() != Some(&bound.source)
            || intent_form.run(()) != bound.form
            || !anchor_available.run(bound.source.clone())
        {
            intent_status.set(
                "Original source or edit settings changed. Preview the intended edit again.".into(),
            );
            return;
        }
        match on_intent.run(IntentOpen{source:bound.source.clone(),report:bound.report.clone()}) {
            Ok(())=>intent_status.set("Intended reference opened. It is a design requirement, not a reconstructed experimental product; design and verify the inverse-PCR plan next.".into()),Err(error)=>intent_status.set(error),
        }
    };
    Effect::new(move |_| {
        let _ = (
            pair_draft.get(),
            source_binding.get(),
            inverse_origin.get(),
            target_binding.get(),
            epoch.get(),
        );
        applied_pair.set(None);
    });
    let apply_pair = Callback::new(move |pair: (Primer, Primer)| {
        let binding = source_binding
            .get_untracked()
            .ok_or("Choose a current source template before applying primers.")?;
        if !anchor_available.run(binding.clone()) {
            return Err(
                "The source event is no longer available. Select the template again.".into(),
            );
        }
        applied_pair.set(Some(AppliedPair {
            primers: pair,
            draft: pair_draft.get_untracked(),
            source: binding,
            target: target_binding.get_untracked(),
            epoch: epoch.get_untracked(),
        }));
        status.set("Primer pair applied. Simulate with the declared mismatch rules.".into());
        Ok(())
    });
    Effect::new(move |_| {
        let _ = (
            molecule.get().revision(),
            source.get().map(|m| m.revision()),
            template_index.get(),
            applied_pair.get(),
            pair_draft.get(),
            source_binding.get(),
            inverse_origin.get(),
            target_binding.get(),
            epoch.get(),
            mismatches.get(),
            exact.get(),
            edit_start.get(),
            edit_end.get(),
            inserted.get(),
            f_length.get(),
            r_length.get(),
            split.get(),
        );
        generation.update(|n| *n = n.wrapping_add(1));
        pcr.set(None);
        inverse.set(None);
        busy.set(false);
        status.set(String::new());
    });
    let run = {
        let alive = alive.clone();
        move |job: AnalysisJob, expected_source: Molecule, target_revision: Option<String>| {
            let Some(binding) = source_binding.get_untracked() else {
                status.set("Choose a current source template.".into());
                return;
            };
            if binding.revision != expected_source.revision()
                || !anchor_available.run(binding.clone())
                || !anchor_available.run(target_binding.get_untracked())
            {
                status.set(
                    "Source or target event changed. Select the current documents again.".into(),
                );
                return;
            }
            let token = generation.get_untracked().wrapping_add(1);
            generation.set(token);
            busy.set(true);
            status.set("Computing mutation evidence in background…".into());
            let source_revision = expected_source.revision();
            let inverse_job = matches!(&job, AnalysisJob::InversePcr { .. });
            let captured_origin = inverse_origin.get_untracked();
            let revision = match &job {
                AnalysisJob::InversePcr {
                    molecule,
                    target,
                    origin,
                    options,
                } => match collagene_core::cloning::pcr::mutagenesis::inverse_pcr_revision(
                    molecule,
                    target,
                    origin.as_ref(),
                    options,
                ) {
                    Ok(v) => v,
                    Err(e) => {
                        busy.set(false);
                        status.set(e.to_string());
                        return;
                    }
                },
                _ => source_revision.clone(),
            };
            let captured_source = source_binding.get_untracked();
            let captured_target = target_binding.get_untracked();
            let captured_epoch = epoch.get_untracked();
            let captured_draft = pair_draft.get_untracked();
            let captured_form = inverse_form.run(());
            let captured_pcr_form = (mismatches.get_untracked(), exact.get_untracked());
            let captured_oligos = match &job {
                AnalysisJob::MismatchPcrPreview {
                    forward, reverse, ..
                } => Some((forward.clone(), reverse.clone())),
                _ => None,
            };
            let captured_molecule = Arc::new(expected_source);
            let alive = alive.clone();
            spawn_local(async move {
                let response = worker::analyze_cancellable(
                    Request {
                        id: token,
                        revision: revision.clone(),
                        job,
                    },
                    alive.clone(),
                )
                .await;
                if !alive.load(Ordering::Relaxed)
                    || generation.get_untracked() != token
                    || epoch.get_untracked() != captured_epoch
                    || pair_draft.get_untracked() != captured_draft
                    || inverse_form.run(()) != captured_form
                    || (mismatches.get_untracked(), exact.get_untracked()) != captured_pcr_form
                    || (inverse_job && inverse_origin.get_untracked() != captured_origin)
                    || source_binding.get_untracked() != captured_source
                    || target_binding.get_untracked() != captured_target
                {
                    return;
                }
                if source
                    .get_untracked()
                    .is_none_or(|m| m.revision() != source_revision)
                    || target_revision.is_some_and(|r| molecule.get_untracked().revision() != r)
                {
                    return;
                }
                busy.set(false);
                match response {
                    Ok(Response {
                        id,
                        revision: returned,
                        result,
                    }) if id == token && returned == revision => match result {
                        AnalysisResult::MismatchPcrPreview(preview) => {
                            let report = preview.report;
                            status.set(format!("{} candidate geometries · search complete: {} · complete under declared model: {}. No empirical amplification claim.",report.products.len(),report.search_complete,report.enumeration_complete));
                            let Some(source) = captured_source else {
                                status.set(
                                    "The original mismatch-PCR source event is missing.".into(),
                                );
                                return;
                            };
                            let Some((forward, reverse)) = captured_oligos else {
                                status.set("Mismatch-PCR oligo bindings are unavailable.".into());
                                return;
                            };
                            if report.source_revision != source.revision
                                || report.products.len() != preview.documents.len()
                            {
                                status.set("Mismatch-PCR preview differs from the requested source or geometry count.".into());
                                return;
                            }
                            pcr.set(Some(Arc::new(PcrEvidence {
                                review: Arc::new(
                                    crate::primers::pcr::mismatch_pcr_outputs::Review {
                                        source,
                                        target: captured_target,
                                        epoch: captured_epoch,
                                        source_molecule: captured_molecule,
                                        forward,
                                        reverse,
                                        report,
                                        documents: preview.documents,
                                        findings: preview.findings,
                                    },
                                ),
                                generation: token,
                                form: captured_pcr_form,
                                draft: captured_draft,
                                inverse_form: captured_form,
                            })));
                        }
                        AnalysisResult::InversePcr(report) => {
                            if let Err(e) = report.validate_binding(&revision) {
                                status.set(e.to_string());
                                return;
                            }
                            if let Some(event) = &captured_source {
                                if let Err(e) =
                                    attest_inverse.run((event.clone(), captured_origin.clone()))
                                {
                                    status.set(e);
                                    return;
                                }
                            }
                            status.set(format!("Independent inverse-PCR plan check: {:?}. Review the bound steps and target before use.",report.verification.status));
                            let Some(source) = captured_source else {
                                status.set(
                                    "The original inverse-PCR source event is missing.".into(),
                                );
                                return;
                            };
                            if report.source_revision != source.revision
                                || report.target_revision != captured_target.revision
                                || report.plan.target_revision != captured_target.revision
                                || report.plan.inputs.len() != 1
                                || report.plan.inputs[0].revision != source.revision
                            {
                                status.set(
                                    "Inverse-PCR source or intended-target binding changed.".into(),
                                );
                                return;
                            }
                            inverse.set(Some(Arc::new(BoundInverse {
                                report,
                                origin: captured_origin,
                                source,
                                target: captured_target,
                                epoch: captured_epoch,
                                generation: token,
                                form: captured_form,
                                draft: captured_draft,
                            })));
                        }
                        AnalysisResult::Error(error) => status.set(error),
                        _ => status.set("Unexpected mutation worker result".into()),
                    },
                    Err(error) => status.set(error),
                    _ => status.set("Mutation worker identity changed".into()),
                }
            });
        }
    };
    let run_pcr = run.clone();
    let simulate = move |_| {
        let result = (|| -> Result<_, String> {
            let source = source.get_untracked().ok_or("Choose a template document")?;
            let applied = applied_pair
                .get_untracked()
                .ok_or("Apply a primer pair before simulating mismatch PCR.")?;
            if applied.draft != pair_draft.get_untracked()
                || source_binding.get_untracked().as_ref() != Some(&applied.source)
                || applied.target != target_binding.get_untracked()
                || applied.epoch != epoch.get_untracked()
            {
                return Err("Primer draft or source/target event changed. Apply the pair again before simulation.".into());
            }
            let (forward, reverse) = applied.primers;
            let options = PcrOptions {
                binding: BindingOptions {
                    max_mismatches: number(mismatches.get_untracked(), "Allowed mismatches")?,
                    min_three_prime_exact: number(exact.get_untracked(), "Exact 3′ bases")?,
                    ..Default::default()
                },
                ..Default::default()
            };
            Ok((source, forward, reverse, options))
        })();
        match result {
            Ok((source, forward, reverse, options)) => run_pcr(
                AnalysisJob::MismatchPcrPreview {
                    molecule: source.clone(),
                    forward: forward.clone(),
                    reverse: reverse.clone(),
                    options,
                },
                source,
                None,
            ),
            Err(e) => status.set(e),
        }
    };
    let design = move |_| {
        let result = (|| -> Result<_, String> {
            let source = source.get_untracked().ok_or("Choose a source template")?;
            let target = molecule.get_untracked();
            let start: usize = number(edit_start.get_untracked(), "Edit start")?;
            let end: usize = number(edit_end.get_untracked(), "Edit end")?;
            let options = InversePcrOptions {
                edit_start: start.checked_sub(1).ok_or("Edit start is one-based")?,
                edit_end: end.checked_sub(1).ok_or("Edit end is one-based")?,
                inserted: inserted.get_untracked(),
                forward_length: number(f_length.get_untracked(), "Forward domain length")?,
                reverse_length: number(r_length.get_untracked(), "Reverse domain length")?,
                reverse_tail_bases: number(split.get_untracked(), "Reverse-tail split")?,
            };
            let event = source_binding
                .get_untracked()
                .ok_or("Choose an exact source event")?;
            let origin = inverse_origin.get_untracked();
            attest_inverse.run((event, origin.clone()))?;
            Ok((source, target, options, origin))
        })();
        match result {
            Ok((source, target, options, origin)) => {
                let revision = target.revision();
                run(
                    AnalysisJob::InversePcr {
                        origin: origin.as_ref().map(|o| o.as_ref().clone()),
                        molecule: source.clone(),
                        target,
                        options,
                    },
                    source,
                    Some(revision),
                );
            }
            Err(e) => status.set(e),
        }
    };
    let current_pcr = Callback::new(move |bound: Arc<PcrEvidence>| -> Result<(), String> {
        let r = &bound.review;
        if generation.get_untracked() != bound.generation
            || epoch.get_untracked() != r.epoch
            || pair_draft.get_untracked() != bound.draft
            || inverse_form.run(()) != bound.inverse_form
            || (mismatches.get_untracked(), exact.get_untracked()) != bound.form
            || source_binding.get_untracked().as_ref() != Some(&r.source)
            || target_binding.get_untracked() != r.target
            || !anchor_available.run(r.source.clone())
            || !anchor_available.run(r.target.clone())
            || source
                .get_untracked()
                .is_none_or(|source| source.revision() != r.source.revision)
            || !pcr.with_untracked(|current| {
                current
                    .as_ref()
                    .is_some_and(|current| Arc::ptr_eq(current, &bound))
            })
        {
            return Err("The primer pair, template or PCR settings changed; apply and simulate again before using this preview.".into());
        }
        let applied = applied_pair
            .get_untracked()
            .ok_or("Apply the named primer pair again before using PCR outputs.")?;
        if applied.draft != bound.draft
            || applied.source != r.source
            || applied.target != r.target
            || applied.epoch != r.epoch
        {
            return Err(
                "The applied primer pair no longer belongs to this source and reference context."
                    .into(),
            );
        }
        Ok(())
    });
    let review_inverse = Callback::new(move |bound: Arc<BoundInverse>| {
        if epoch.get_untracked() != bound.epoch
            || inverse_origin.get_untracked() != bound.origin
            || generation.get_untracked() != bound.generation
            || inverse_form.run(()) != bound.form
            || pair_draft.get_untracked() != bound.draft
            || source_binding.get_untracked().as_ref() != Some(&bound.source)
            || target_binding.get_untracked() != bound.target
            || !anchor_available.run(bound.source.clone())
            || !anchor_available.run(bound.target.clone())
            || source
                .get_untracked()
                .is_none_or(|source| source.revision() != bound.source.revision)
            || molecule.get_untracked().revision() != bound.target.revision
            || !inverse.with_untracked(|current| {
                current
                    .as_ref()
                    .is_some_and(|current| Arc::ptr_eq(current, &bound))
            })
        {
            status.set(
                "The inverse-PCR draft, source or intended target changed; design it again.".into(),
            );
            return;
        }
        if let Err(e) = attest_inverse.run((bound.source.clone(), bound.origin.clone())) {
            status.set(e);
            return;
        }
        if let Err(e) = bound
            .report
            .validate_binding(&bound.report.context_revision)
        {
            status.set(e.to_string());
            return;
        }
        let candidate = Candidate::new(bound.report.plan.clone());
        if let Err(error) = candidate.validate() {
            status.set(error.to_string());
            return;
        }
        if candidate.plan.inputs.len() != 1
            || candidate.plan.inputs[0].revision != bound.source.revision
            || candidate.plan.inputs[0].molecule.revision() != bound.source.revision
            || candidate.plan.target_revision != bound.target.revision
        {
            status.set(
                "The inverse-PCR candidate differs from its original source or intended target."
                    .into(),
            );
            return;
        }
        let input_id = candidate.plan.inputs[0].id.clone();
        on_candidate.run(CandidateOpen {
            candidate,
            sources: vec![InputEventBinding {
                input_id,
                event: bound.source.clone(),
            }],
            target: bound.target.clone(),
            epoch: bound.epoch,
        });
    });
    let reusable_settings = Callback::new(
        move |()| -> Result<crate::restriction::reusable_material_view::Settings, String> {
            let applied = applied_pair
                .get_untracked()
                .ok_or("Apply the actual ordered primer pair first.")?;
            if applied.draft != pair_draft.get_untracked()
                || source_binding.get_untracked().as_ref() != Some(&applied.source)
                || applied.epoch != epoch.get_untracked()
                || applied.target != target_binding.get_untracked()
            {
                return Err(
                    "The primer draft or source context changed; apply the primer pair again."
                        .into(),
                );
            }
            Ok(
                crate::restriction::reusable_material_view::Settings::Mismatch {
                    forward: applied.primers.0,
                    reverse: applied.primers.1,
                    options: PcrOptions {
                        binding: BindingOptions {
                            max_mismatches: number(
                                mismatches.get_untracked(),
                                "Allowed mismatches",
                            )?,
                            min_three_prime_exact: number(exact.get_untracked(), "Exact 3′ bases")?,
                            ..Default::default()
                        },
                        ..Default::default()
                    },
                },
            )
        },
    );
    let reusable_settings_revision = Memo::new(move |_| {
        serde_json::to_string(&(
            pair_draft.get(),
            applied_pair.get().map(|p| p.primers),
            mismatches.get(),
            exact.get(),
            source_binding.get(),
            inverse_origin.get(),
            target_binding.get(),
            epoch.get(),
        ))
        .unwrap_or_default()
    });
    view! {
        <style>".mutation-panel .mutation-grid{display:grid;grid-template-columns:1fr 1fr;gap:8px}.mutation-panel pre{white-space:pre-wrap;overflow-wrap:anywhere;font:11px/1.6 monospace;background:#f3f6f3;padding:9px;max-height:250px;overflow:auto}.mutation-panel .mutation-product{padding:10px 0;border-bottom:1px solid #e2e9e2;font-size:12px;line-height:1.6}.mutation-panel .mutation-status{font-size:12px;line-height:1.6;overflow-wrap:anywhere}.mutation-panel ul{padding-left:18px}.mutation-panel details{margin:10px 0}"</style>
        <div class="inspector-content mutation-panel"><h2>"Mutagenesis"</h2><label>"Source template"<select aria-label="Mutagenesis template" prop:value=move||template_index.get().to_string() on:change=move|e|template_index.set(event_target_value(&e).parse().unwrap_or(0))>{move||query_options.get().into_iter().map(|(i,name)|view!{<option value=i.to_string()>{name}</option>}).collect_view()}</select></label><p class="mutation-status" role="status">{move||status.get()}</p>
        <section class="mutation-intent"><h3>"1. Set the intended edit"</h3>
            <p class="helper">"Select the original DNA range to replace, or a caret for insertion. Preview the declared intended reference, open it as your target, then design inverse PCR below. An independently imported target can also be used; creating one here is optional."</p>
            <button class="button secondary full" on:click=use_mutation_selection>"Use selected source range/caret"</button>
            <label>"Intended reference name"<input maxlength="256" aria-label="Mutation intended reference name" prop:value=move||intended_name.get() on:input=move|e|intended_name.set(event_target_value(&e))/></label>
            <p class="helper">"Edit coordinates are one-based boundaries: start inclusive, end exclusive. Equal values insert without deletion. A circular end can exceed the source length to cross the origin. Example: start 4, end 6 replaces source bases 4–5. Selection size is independent of primer length."</p>
            <div class="mutation-grid">{[("Edit start boundary","Mutation edit start",edit_start),("Edit end boundary","Mutation edit end",edit_end)].into_iter().map(|(label,aria,value)|view!{<label>{label}<input aria-label=aria type="number" min="1" prop:value=move||value.get() on:input=move|e|value.set(event_target_value(&e))/></label>}).collect_view()}</div>
            <label>"Replacement/inserted sequence"<textarea maxlength="200000" aria-label="Mutation inserted sequence" rows="3" prop:value=move||inserted.get() on:input=move|e|inserted.set(event_target_value(&e))></textarea></label>
            <button class="button primary full" disabled=move||intent_busy.get()||busy.get()||source.get().is_none() on:click=preview_intent>"Preview intended edit"</button>
            <button class="text-tool" disabled=move||!intent_busy.get() on:click=move|_|{intent_generation.update(|g|*g=g.wrapping_add(1));if let Some(w)=intent_worker.get_value(){w.close();}intent_worker.set_value(None);intended.set(None);intent_busy.set(false);intent_status.set("Intended-edit preview cancelled.".into());}>"Cancel intended edit"</button>
            <p class="mutation-status" role="status">{move||intent_status.get()}</p>
            {move||intended.get().map(|bound|{let replacement=bound.report.request.inserted.chars().filter(|call|!call.is_ascii_whitespace()).map(|call|call.to_ascii_uppercase()).collect::<String>();view!{<article class="mutation-intent-preview"><h4>"2. Review the intended reference"</h4>
                <p>{format!("Source edit boundaries {}–{} (1-based, end exclusive) · {} deleted bases · {} inserted bases · {} desired bases",bound.report.request.edit_start+1,bound.report.request.edit_end+1,bound.report.deleted.len(),replacement.len(),bound.report.target.sequence.len())}</p>
                <pre>{format!("Deleted source DNA: {}\nReplacement DNA: {}\nDesired target DNA: {}",dna_preview(&bound.report.deleted),dna_preview(&replacement),dna_preview(bound.report.target.sequence.as_str()))}</pre>
                <p class="helper">{format!("DNA previews show at most 200 bases each. Origin shift: {} source bases (zero-based); an origin-crossing edit can rotate the intended reference. This is source plus your literal edit, not a PCR or ligation result.",bound.report.origin_shift)}</p>
                <details><summary>"Original source event"</summary><p>{format!("Document {} · event {} · revision {}",bound.source.document_id,bound.source.event_id,bound.source.revision)}</p></details>
                {bound.report.warnings.iter().map(|warning|view!{<p class="helper">{warning.clone()}</p>}).collect_view()}
                <button class="button primary full" disabled=move||intent_busy.get()||busy.get() on:click=open_intent>"Open intended reference"</button>
            </article>}})}
        </section>
        <crate::cloning::preparation::overlap_extension_preparation_view::ReusableOverlapExtensionMaterialView source source_info=Memo::new(move|_|source_info.run(template_index.get())) epoch get_workspace on_open=on_material/>
        <h3>"PCR with declared mismatches"</h3><p class="helper">"Enter actual ordered oligo domains 5′→3′, including intended substitutions. Binding mismatches are conditional rules, not predicted polymerase success."</p>
        <crate::primers::primer_pair_editor::PrimerPairEditor draft=pair_draft selection=pair_selection anchor_available=anchor_available resolve_source primer_sources=primer_sources list_primers=list_primers resolve_primer=resolve_primer epoch=epoch on_apply=apply_pair apply_label="Apply pair to mismatch PCR".to_string()/>
        <p class="helper" role="status">{move || applied_pair.get().map(|pair| format!("Applied: {} / {}. Draft edits require Apply again.",pair.primers.0.name,pair.primers.1.name)).unwrap_or_else(|| "No applied primer pair. Edit or reuse primers, then Apply before simulation.".into())}</p>
        <div class="mutation-grid"><label>"Allowed disjoint mismatches"<input aria-label="Mutagenesis allowed mismatches" type="number" min="0" prop:value=move||mismatches.get() on:input=move|e|mismatches.set(event_target_value(&e))/></label><label>"Required exact 3′ bases"<input aria-label="Mutagenesis exact 3prime bases" type="number" min="0" prop:value=move||exact.get() on:input=move|e|exact.set(event_target_value(&e))/></label></div><button class="button primary full" disabled=move||busy.get()||intent_busy.get()||source.get().is_none()||applied_pair.get().is_none() on:click=simulate>"Simulate mismatch PCR"</button>
        <details class="assembly-fragment-options"><summary>"Inverse-PCR template material"</summary><crate::history::run_origin_view::SourceOriginPicker source=Memo::new(move|_|source_info.run(template_index.get())) epoch get_workspace selected=Memo::new(move|_|inverse_origin.get()) on_select=Callback::new(move|value:Option<crate::history::run_origin_model::PinnedOrigin>|{if let Some(pin)=&value{if source_binding.get_untracked().as_ref()!=Some(pin.event()){return Err("Selected material belongs to another exact template event.".into());}}inverse_origin.set(value);Ok(())})/><p class="helper">{move||inverse_origin.get().map(|o|format!("{} · source event {}…",o.label(),o.event().event_id.chars().take(12).collect::<String>())).unwrap_or_else(||"Explicit original-reference template; no saved material selected.".into())}</p><p class="helper">"Inverse PCR requires an intact, fully paired circular template under the native model. Saved material is freshly replayed. Unsupported nicks or staggered geometry are refused; no repair or linearization is inserted automatically."</p></details>
        <crate::restriction::reusable_material_view::ReusableMaterialView source source_info=Memo::new(move|_|source_info.run(template_index.get())) settings=reusable_settings settings_revision=reusable_settings_revision kind=crate::restriction::reusable_material_view::Kind::MismatchPcr origin=RwSignal::new(None) epoch get_workspace on_open=on_material/>
        {move||pcr.get().map(|bound|{
            let review=bound.review.clone();
            let guard=Callback::new(move|()|current_pcr.run(bound.clone()));
            view!{<crate::primers::pcr::mismatch_pcr_outputs::MismatchPcrOutputsPanel review guard on_outputs on_overlay on_navigate/>}
        })}

        <h3>"3. Design inverse-PCR preparation"</h3><p class="helper">{move||format!("Active document is the independently specified intended circular target: {}. Choose the original circular source above.",molecule.get().name)}</p>
        <div class="mutation-grid">{[("Forward annealing length","Mutation forward length",f_length),("Reverse annealing length","Mutation reverse length",r_length),("Inserted bases on reverse tail","Mutation tail split",split)].into_iter().map(|(label,aria,value)|view!{<label>{label}<input aria-label=aria type="number" prop:value=move||value.get() on:input=move|e|value.set(event_target_value(&e))/></label>}).collect_view()}</div><button class="button primary full" disabled=move||busy.get()||intent_busy.get()||source.get().is_none() on:click=design>"Design inverse-PCR candidate"</button>
        {move||inverse.get().map(|bound|{let report=&bound.report;let action=bound.clone();let download=bound.clone();view!{<section class="inverse-pcr-result"><h3>{format!("Inverse-PCR check: {:?}",report.verification.status)}</h3><pre>{format!("F  {}\nR  {}",report.forward.ordered_sequence(),report.reverse.ordered_sequence())}</pre><p class="helper">{format!("Retained template {} bp · intended target {}…",report.retained_template_length,report.target_revision.chars().take(12).collect::<String>())}</p><ol><li>"Exact back-to-back PCR with explicit tails"</li><li>"Explicit terminal phosphorylation"</li><li>"Sealed blunt circular ligation"</li></ol><ul>{report.verification.findings.iter().chain(report.warnings.iter()).cloned().map(|w|view!{<li>{w}</li>}).collect_view()}</ul><button class="button primary full" on:click=move|_|review_inverse.run(action.clone())>"Review inverse-PCR candidate"</button><button class="button secondary full" on:click=move|_|match serde_json::to_vec_pretty(&download.report){Ok(bytes)=>crate::io::download::download_with_completion("inverse-pcr-design.json","application/json",&bytes,move|result|{if status.try_get_untracked().is_none(){return;}if let Err(e)=result {status.set(e)}}),Err(e)=>status.set(e.to_string())}>"Download inverse-PCR evidence"</button></section>}})}
        </div>
    }
}
