//! Explicit decisions over one native relationship preview, committed as an edit.
use super::{
    attribution::contribution_text,
    model::{self, Context},
};
use crate::{platform::worker, sources::identity::SourceInfo};
use collagene_core::{
    documents::reference_edit::Operation,
    material::lineage::EventBinding,
    protocol::{
        request::Request,
        schema::{AnalysisJob, AnalysisResult},
    },
    sequence::annotations::merge::{self, Candidate, Decision, Preview, RelationKind, Source},
};
use collagene_model::Molecule;
use leptos::prelude::*;
use std::{
    collections::BTreeMap,
    sync::{
        Arc,
        atomic::{AtomicBool, Ordering},
    },
};
use wasm_bindgen_futures::spawn_local;

#[derive(Clone)]
struct Bound {
    request: merge::Request,
    preview: Preview,
    context: Context,
}

#[component]
pub(crate) fn MergeReview(
    molecule: Memo<Molecule>,
    info: Memo<SourceInfo>,
    epoch: RwSignal<u64>,
    locked: Memo<bool>,
    input: Memo<(u32, Vec<String>)>,
    prepare: Callback<(), Result<(Source, Vec<String>), String>>,
    on_commit: Callback<(EventBinding, Vec<Operation>), Result<(), String>>,
) -> impl IntoView {
    let shown = RwSignal::new(None::<Arc<Bound>>);
    let choices = RwSignal::new(BTreeMap::<String, Decision>::new());
    let generation = RwSignal::new(0u32);
    let page = RwSignal::new(0usize);
    let busy = RwSignal::new(false);
    let status = RwSignal::new(String::new());
    let alive = Arc::new(AtomicBool::new(true));
    {
        let alive = alive.clone();
        on_cleanup(move || alive.store(false, Ordering::Relaxed));
    }
    Effect::new(move |_| {
        let _ = (info.get(), epoch.get(), input.get(), locked.get());
        generation.update(|value| *value = value.wrapping_add(1));
        shown.set(None);
        choices.set(BTreeMap::new());
        busy.set(false);
        page.set(0);
    });
    let guard = Callback::new(move |bound: Arc<Bound>| -> Result<(), String> {
        let info = info.try_get_untracked().ok_or("The merge panel closed.")?;
        let epoch = epoch.try_get_untracked().ok_or("The merge panel closed.")?;
        let input = input.try_get_untracked().ok_or("The merge panel closed.")?;
        bound.context.validate_current(
            &info.binding,
            epoch,
            &input,
            locked.try_get_untracked().unwrap_or(true),
        )
    });
    let review = {
        let alive = alive.clone();
        Callback::new(move |()| {
            if locked.get_untracked() || busy.get_untracked() {
                return;
            }
            let prepared = (|| -> Result<_, String> {
                let (source, ids) = prepare.run(())?;
                source.preflight().map_err(|error| error.to_string())?;
                if ids.is_empty() {
                    return Err("Select at least one proposal.".into());
                }
                let request = merge::Request {
                    schema_version: 1,
                    target_revision: info.get_untracked().binding.revision,
                    source,
                    selections: ids
                        .into_iter()
                        .map(|proposal_id| merge::Selection {
                            proposal_id,
                            decision: Decision::Skip,
                        })
                        .collect(),
                    // The library scanner already admits up to this shared work
                    // ceiling. Keep room for both its fresh replay and target
                    // relationships rather than refusing ordinary larger scans
                    // through the smaller generic API default.
                    limits: merge::Limits {
                        max_comparisons:
                            collagene_core::sequence::annotations::library::MAX_SCAN_COMPARISONS,
                        ..merge::Limits::default()
                    },
                };
                let (target, revision) = molecule.with_untracked(|target| {
                    request
                        .validate(target)
                        .map_err(|error| error.to_string())?;
                    let revision =
                        merge::revision(target, &request).map_err(|error| error.to_string())?;
                    Ok::<_, String>((target.clone(), revision))
                })?;
                Ok((
                    target,
                    request,
                    revision,
                    info.get_untracked().binding,
                    epoch.get_untracked(),
                    input.get_untracked(),
                ))
            })();
            let (target, request, revision, event, captured_epoch, captured_input) = match prepared
            {
                Ok(value) => value,
                Err(error) => {
                    status.set(error);
                    return;
                }
            };
            let id = generation.get_untracked().wrapping_add(1);
            generation.set(id);
            busy.set(true);
            shown.set(None);
            status.set("Inspecting annotation relationships in the Rust worker…".into());
            let alive = alive.clone();
            spawn_local(async move {
                let response = worker::analyze(Request {
                    id,
                    revision: revision.clone(),
                    job: AnalysisJob::PreviewAnnotationMerge {
                        molecule: target,
                        request: request.clone(),
                    },
                })
                .await;
                if !alive.load(Ordering::Relaxed) || generation.get_untracked() != id {
                    return;
                }
                busy.set(false);
                let accepted = (|| -> Result<Arc<Bound>, String> {
                    let response = response?;
                    if response.id != id || response.revision != revision {
                        return Err("Annotation preview binding changed.".into());
                    }
                    let preview = match response.result {
                        AnalysisResult::AnnotationMergePreview(preview) => preview,
                        AnalysisResult::Error(error) => return Err(error),
                        _ => return Err("Unexpected annotation preview response.".into()),
                    };
                    preview
                        .validate_binding(&revision)
                        .map_err(|error| error.to_string())?;
                    molecule
                        .with_untracked(|target| preview.validate_request(target, &request))
                        .map_err(|error| error.to_string())?;
                    let bound = Arc::new(Bound {
                        request,
                        preview,
                        context: Context {
                            event,
                            epoch: captured_epoch,
                            input: captured_input,
                        },
                    });
                    guard.run(bound.clone())?;
                    Ok(bound)
                })();
                match accepted {
                    Ok(bound) => {
                        choices.set(
                            bound
                                .request
                                .selections
                                .iter()
                                .map(|selection| (selection.proposal_id.clone(), Decision::Skip))
                                .collect(),
                        );
                        shown.set(Some(bound));
                        status.set("Choose append, skip or an eligible existing annotation for every proposal. Applying rechecks the exact choices.".into());
                    }
                    Err(error) => status.set(error),
                }
            });
        })
    };
    let apply = {
        let alive = alive.clone();
        Callback::new(move |()| {
            if busy.get_untracked() {
                return;
            }
            let prepared = (|| -> Result<_, String> {
                let bound = shown
                    .get_untracked()
                    .ok_or("Preview merge relationships first.")?;
                guard.run(bound.clone())?;
                let mut request = bound.request.clone();
                let decisions = choices.get_untracked();
                request.selections = model::selections(&request.selections, &decisions)?;
                if request
                    .selections
                    .iter()
                    .all(|selection| matches!(selection.decision, Decision::Skip))
                {
                    return Err("Every proposal is skipped; no edit is needed.".into());
                }
                let (target, revision) = molecule.with_untracked(|target| {
                    request
                        .validate(target)
                        .map_err(|error| error.to_string())?;
                    let revision =
                        merge::revision(target, &request).map_err(|error| error.to_string())?;
                    Ok::<_, String>((target.clone(), revision))
                })?;
                Ok((bound, request, target, revision, decisions))
            })();
            let (bound, request, target, revision, decisions) = match prepared {
                Ok(value) => value,
                Err(error) => {
                    status.set(error);
                    return;
                }
            };
            let id = generation.get_untracked().wrapping_add(1);
            generation.set(id);
            busy.set(true);
            let alive = alive.clone();
            spawn_local(async move {
                let response = worker::analyze(Request {
                    id,
                    revision: revision.clone(),
                    job: AnalysisJob::ApplyAnnotationMerge {
                        molecule: target,
                        request: request.clone(),
                    },
                })
                .await;
                if !alive.load(Ordering::Relaxed) || generation.get_untracked() != id {
                    return;
                }
                busy.set(false);
                let outcome = (|| -> Result<(), String> {
                    guard.run(bound.clone())?;
                    if choices.get_untracked() != decisions {
                        return Err("Merge choices changed while checking. Apply again.".into());
                    }
                    let response = response?;
                    if response.id != id || response.revision != revision {
                        return Err("Annotation apply binding changed.".into());
                    }
                    let applied = match response.result {
                        AnalysisResult::AnnotationMergeApplied(applied) => applied,
                        AnalysisResult::Error(error) => return Err(error),
                        _ => return Err("Unexpected annotation apply response.".into()),
                    };
                    molecule
                        .with_untracked(|target| applied.validate_request(target, &request))
                        .map_err(|error| error.to_string())?;
                    if applied.review.context_revision != revision {
                        return Err("Annotation apply evidence belongs to another request.".into());
                    }
                    // The canonical native edit executes the user's actual request;
                    // a returned molecule is never installed as a snapshot edit.
                    on_commit.run((
                        bound.context.event.clone(),
                        vec![Operation::MergeAnnotations { request }],
                    ))
                })();
                match outcome {
                    Ok(()) => {
                        shown.set(None);
                        status.set(
                            "Merge saved as one undoable annotation edit. DNA is unchanged.".into(),
                        );
                    }
                    Err(error) => status.set(error),
                }
            });
        })
    };
    view! {
        <section class="annotation-merge">
            <button class="button secondary" disabled=move || busy.get() || locked.get()
                on:click=move |_| review.run(())>"Preview selected annotation relationships"</button>
            {move || shown.get().map(|bound| {
                let count = bound.preview.candidates.len();
                let pages = count.div_ceil(20).max(1);
                let current = page.get().min(pages - 1);
                let rows = bound.preview.candidates.iter().skip(current * 20).take(20)
                    .map(|candidate| view! { <MergeCandidate candidate=candidate.clone() choices/> }).collect_view();
                view! {
                    <p>{format!("Relationship search complete: {} · source proposals complete: {} · {} comparisons",
                        bound.preview.complete, bound.preview.source_complete, bound.preview.comparisons)}</p>
                    {bound.preview.warnings.iter().cloned().map(|warning| view! { <p>{warning}</p> }).collect_view()}
                    {rows}
                    <button class="text-tool" disabled=current == 0 on:click=move |_| page.set(current.saturating_sub(1))>"Previous proposals"</button>
                    <span>{format!("Page {} of {pages}", current + 1)}</span>
                    <button class="text-tool" disabled=current + 1 >= pages on:click=move |_| page.set((current + 1).min(pages - 1))>"Next proposals"</button>
                    <p class="helper">"Choose a decision for each proposal, then apply the merge."</p>
                    <details><summary>"What these choices mean"</summary>
                        <p>"Append creates an annotation requiring review. Skip changes nothing.
                        Adding provenance preserves the existing annotation and its review state;
                        it does not reconcile conflicting claims or verify biology.
                        An incomplete relationship search cannot authorize a change."</p>
                    </details>
                    <button class="button primary" disabled=move || busy.get() || locked.get()
                        on:click=move |_| apply.run(())>"Apply explicit merge choices"</button>
                }
            })}
            <p class="helper" role="status">{move || status.get()}</p>
        </section>
    }
}

#[component]
fn MergeCandidate(
    candidate: Candidate,
    choices: RwSignal<BTreeMap<String, Decision>>,
) -> impl IntoView {
    let candidate = StoredValue::new(candidate);
    let page = RwSignal::new(0usize);
    let decision = Memo::new(move |_| {
        candidate.with_value(|candidate| {
            choices.with(|choices| choices.get(&candidate.proposal_id).cloned())
        })
    });
    view! {
        <article>
            <h4>{candidate.with_value(|candidate| format!("{} · {}", candidate.feature.label(), candidate.feature.kind))}</h4>
            <details><summary>"Source contribution"</summary>
                <p>{candidate.with_value(|candidate| contribution_text(&candidate.contribution))}</p>
            </details>
            <p>{move || match decision.get() {
                Some(Decision::Append) => "Decision: append a new annotation".into(),
                Some(Decision::Skip) => "Decision: skip".into(),
                Some(Decision::AddProvenance { target_index, .. }) => format!("Decision: add provenance to annotation {}", target_index + 1),
                None => "No decision selected".into(),
            }}</p>
            <button class="text-tool" on:click=move |_| candidate.with_value(|candidate| choices.update(|choices| {
                choices.insert(candidate.proposal_id.clone(), Decision::Skip);
            }))>"Skip"</button>
            <button class="text-tool" disabled=candidate.with_value(|candidate| !candidate.relations_complete)
                on:click=move |_| candidate.with_value(|candidate| choices.update(|choices| {
                choices.insert(candidate.proposal_id.clone(), Decision::Append);
            }))>"Append separately"</button>
            <details>
                <summary>"Existing annotation relationships"</summary>
                {move || candidate.with_value(|candidate| {
                    let pages = candidate.relations.len().div_ceil(20).max(1);
                    let current = page.get().min(pages - 1);
                    let rows = candidate.relations.iter().skip(current * 20).take(20).map(|relation| {
                        let eligible = candidate.relations_complete && matches!(relation.kind, RelationKind::ExactDuplicate | RelationKind::EquivalentAnnotation);
                        let proposal = candidate.proposal_id.clone();
                        let choice = Decision::AddProvenance { target_index: relation.target_index, annotation_revision: relation.annotation_revision.clone() };
                        let choose = Callback::new(move |()| choices.update(|choices| { choices.insert(proposal.clone(), choice.clone()); }));
                        view! {
                            <p>{format!("Annotation {}: {:?}", relation.target_index + 1, relation.kind)}</p>
                            <button class="text-tool" disabled=!eligible on:click=move |_| choose.run(())>"Add provenance to this annotation"</button>
                        }
                    }).collect_view();
                    view! {
                        <p>{format!("Relationships complete: {}", candidate.relations_complete)}</p>
                        {rows}
                        <button class="text-tool" disabled=current == 0 on:click=move |_| page.set(current.saturating_sub(1))>"Previous relationships"</button>
                        <span>{format!("Page {} of {pages}", current + 1)}</span>
                        <button class="text-tool" disabled=current + 1 >= pages on:click=move |_| page.set((current + 1).min(pages - 1))>"Next relationships"</button>
                    }
                })}
            </details>
        </article>
    }
}
