//! Lazy workspace capture and guarded worker response; no alignment on the UI thread.
use super::model::{Capture, Context, Draft};
use crate::viewer::alignment::pairwise_result_view::PairwiseResultView;
use crate::viewer::alignment::pairwise_result_view::Review;
use crate::{platform::worker::CachedWorker, sources::identity::SourceInfo};
use collagene_core::alignment::history::Record;
use collagene_core::{
    alignment::ComparisonOptions,
    documents::workspace::Workspace,
    material::lineage::EventBinding,
    protocol::{
        request::Request,
        schema::{AnalysisJob, AnalysisResult},
    },
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

#[derive(Clone)]
pub(super) struct Shown {
    pub capture: Capture,
    pub review: Arc<Review>,
}
#[derive(Clone)]
struct FeatureChoice {
    index: usize,
    label: String,
}

#[component]
pub(crate) fn DesignConsensusPanel(
    design: Memo<EventBinding>,
    consensus_options: Memo<Vec<(usize, String)>>,
    source_info: Callback<usize, Option<SourceInfo>>,
    get_workspace: Callback<(), Result<Workspace, String>>,
    workspace_busy: Memo<bool>,
    epoch: RwSignal<u64>,
    comparison: Memo<Result<ComparisonOptions, String>>,
    on_save: Callback<Record, Result<(), String>>,
    on_navigate: Callback<crate::viewer::alignment::compare::Navigation, Result<(), String>>,
) -> impl IntoView {
    let selected = RwSignal::new(None::<usize>);
    let source_page = RwSignal::new(0usize);
    let draft = RwSignal::new(Draft::default());
    let features = RwSignal::new(Vec::<FeatureChoice>::new());
    let feature_source = RwSignal::new(None::<EventBinding>);
    let report = RwSignal::new(None::<Arc<Shown>>);
    let busy = RwSignal::new(false);
    let status = RwSignal::new(String::new());
    let generation = RwSignal::new(0u32);
    let worker = StoredValue::new_local(None::<Rc<CachedWorker>>);
    let alive = Arc::new(AtomicBool::new(true));
    {
        let alive = alive.clone();
        on_cleanup(move || {
            alive.store(false, Ordering::Relaxed);
            if let Some(worker) = worker.get_value() {
                worker.close();
            }
        });
    }
    let context = Memo::new(move |_| -> Result<Context, String> {
        if workspace_busy.get() {
            return Err("Wait for the workspace update to finish.".into());
        }
        // Reading the source menu tracks head changes without cloning a workspace.
        let choices = consensus_options.get();
        let index = selected
            .get()
            .ok_or("Choose the supplied consensus document.")?;
        if !choices.iter().any(|(candidate, _)| *candidate == index) {
            return Err("The supplied consensus document is unavailable.".into());
        }
        let consensus = source_info
            .run(index)
            .ok_or("The supplied consensus source is unavailable.")?
            .binding;
        Ok(Context {
            design: design.get(),
            consensus,
            epoch: epoch.get(),
            comparison: comparison.get()?,
            draft: draft.get(),
        })
    });
    Effect::new(move |_| {
        let _ = context.get();
        generation.update(|value| *value = value.wrapping_add(1));
        if let Some(worker) = worker.get_value() {
            worker.close();
        }
        worker.set_value(None);
        report.set(None);
        busy.set(false);
    });
    Effect::new(move |_| {
        let current = design.get();
        if feature_source.get_untracked().as_ref() != Some(&current) {
            features.set(Vec::new());
            feature_source.set(None);
        }
    });
    let guard = Callback::new(move |shown: Arc<Shown>| -> Result<(), String> {
        let current = context
            .try_get_untracked()
            .ok_or("The comparison panel closed.")??;
        if !report
            .try_with_untracked(|current| {
                current
                    .as_ref()
                    .is_some_and(|current| Arc::ptr_eq(current, &shown))
            })
            .unwrap_or(false)
        {
            return Err("The displayed comparison changed or closed.".into());
        }
        let workspace = get_workspace.run(())?;
        shown.capture.guard(&workspace, &current)
    });
    let load_features = Callback::new(move |()| {
        let outcome = (|| -> Result<(), String> {
            let mut current = context.get_untracked()?;
            // Loading annotation choices does not require a CDS declaration yet.
            current.draft.coding.enabled = false;
            let workspace = get_workspace.run(())?;
            let capture = Capture::prepare(&workspace, current)?;
            let (design_molecule, _) = capture
                .request
                .sources(&workspace)
                .map_err(|error| error.to_string())?;
            let choices = design_molecule
                .features
                .iter()
                .enumerate()
                .filter(|(_, feature)| feature.kind.eq_ignore_ascii_case("CDS"))
                .take(100)
                .map(|(index, feature)| FeatureChoice {
                    index,
                    label: format!(
                        "{} · {}{}",
                        index + 1,
                        feature.label().chars().take(96).collect::<String>(),
                        if feature.needs_review {
                            " · needs review"
                        } else {
                            ""
                        },
                    ),
                })
                .collect();
            features.set(choices);
            feature_source.set(Some(capture.context.design));
            Ok(())
        })();
        status.set(match outcome {
            Ok(()) => "Exact design annotations loaded. Choose the intended CDS explicitly.".into(),
            Err(error) => error,
        });
    });
    let cancel = Callback::new(move |()| {
        generation.update(|value| *value = value.wrapping_add(1));
        if let Some(worker) = worker.get_value() {
            worker.close();
        }
        worker.set_value(None);
        busy.set(false);
        status.set("Comparison cancelled.".into());
    });
    let compute = {
        let alive = alive.clone();
        Callback::new(move |()| {
            if busy.get_untracked() {
                return;
            }
            let prepared = (|| -> Result<_, String> {
                let current = context.get_untracked()?;
                let workspace = get_workspace.run(())?;
                let capture = Capture::prepare(&workspace, current)?;
                let worker = Rc::new(CachedWorker::new()?);
                Ok((workspace, capture, worker))
            })();
            let (workspace, capture, pending) = match prepared {
                Ok(prepared) => prepared,
                Err(error) => {
                    status.set(error);
                    return;
                }
            };
            let token = generation.get_untracked().wrapping_add(1);
            generation.set(token);
            worker.set_value(Some(pending.clone()));
            busy.set(true);
            report.set(None);
            status.set(
                "Comparing the exact design and supplied consensus in the Rust worker…".into(),
            );
            let alive = alive.clone();
            spawn_local(async move {
                let response = pending
                    .ask(Request {
                        id: token,
                        revision: capture.revision.clone(),
                        job: AnalysisJob::AssessDesignConsensus {
                            workspace,
                            request: capture.request.clone(),
                        },
                    })
                    .await;
                pending.close();
                if !alive.load(Ordering::Relaxed) || generation.get_untracked() != token {
                    return;
                }
                worker.set_value(None);
                busy.set(false);
                let accepted = (|| -> Result<Arc<Review>, String> {
                    let current = context.get_untracked()?;
                    let workspace = get_workspace.run(())?;
                    capture.guard(&workspace, &current)?;
                    let response = response?;
                    if response.id != token || response.revision != capture.revision {
                        return Err("The comparison response belongs to another request.".into());
                    }
                    let result = match response.result {
                        AnalysisResult::DesignConsensus(result) => result,
                        AnalysisResult::Error(error) => return Err(error),
                        _ => return Err("Unexpected design/consensus response.".into()),
                    };
                    capture.validate_report(&workspace, &result)?;
                    let (design, consensus) = capture
                        .request
                        .sources(&workspace)
                        .map_err(|error| error.to_string())?;
                    let record =
                        Record::from_design_consensus(result).map_err(|error| error.to_string())?;
                    Ok(Arc::new(Review {
                        record,
                        reference: Arc::new(design.clone()),
                        query: Arc::new(consensus.clone()),
                        reference_name: capture.design_name.clone(),
                        query_name: capture.consensus_name.clone(),
                        fresh: true,
                        epoch: capture.context.epoch,
                    }))
                })();
                match accepted {
                    Ok(review) => {
                        report.set(Some(Arc::new(Shown { capture, review })));
                        status.set("Comparison ready. Coding summaries are conditional on the selected alignment path.".into());
                    }
                    Err(error) => status.set(error),
                }
            });
        })
    };
    view! {
        <details class="design-consensus">
            <summary>"Design versus supplied consensus"</summary>
            <p class="helper">
                "Use the pinned comparison reference as the design. Choose an explicit supplied consensus;
                this does not load reads, quality scores or experimental verification.
                Whole-construct search settings above are reused."
            </p>
            <p>{move || format!("Design event: {}", design.get().event_id)}</p>
            <label>
                "Supplied consensus"
                <select
                    prop:value=move || selected.get().map(|value| value.to_string()).unwrap_or_default()
                    on:change=move |event| selected.set(event_target_value(&event).parse().ok())
                >
                    <option value="">"Choose a document"</option>
                    {move || {
                        let all = consensus_options.get();
                        let current = source_page.get().min(all.len().div_ceil(100).max(1) - 1);
                        all.into_iter().skip(current * 100).take(100).map(|(index, name)| {
                            view! { <option value=index.to_string()>{name}</option> }
                        }).collect_view()
                    }}
                </select>
            </label>
            <div class="method-row-actions">
                <button class="text-tool" disabled=move || source_page.get() == 0
                    on:click=move |_| source_page.update(|page| *page = page.saturating_sub(1))>
                    "Previous documents"
                </button>
                <span>{move || format!("Document page {}", source_page.get() + 1)}</span>
                <button class="text-tool"
                    disabled=move || source_page.get() + 1 >= consensus_options.get().len().div_ceil(100)
                    on:click=move |_| source_page.update(|page| *page += 1)>
                    "Next documents"
                </button>
            </div>
            <details>
                <summary>"Consensus axis and conditional coding summary"</summary>
                <label>
                    <input type="checkbox" prop:checked=move || draft.get().single_circular_turn
                        on:change=move |event| draft.update(|draft| {
                            draft.single_circular_turn = event_target_checked(&event);
                        })/>
                    "Declare that this consensus contains one circular turn"
                </label>
                <p class="helper">"Stored topology is used unless you explicitly supply this assumption. It does not establish a circular sample."</p>
                <label>
                    "Circular-turn rationale"
                    <input maxlength="1024" prop:value=move || draft.get().axis_rationale
                        on:input=move |event| draft.update(|draft| {
                            draft.axis_rationale = event_target_value(&event);
                        })/>
                </label>
                <label>
                    <input type="checkbox" prop:checked=move || draft.get().coding.enabled
                        on:change=move |event| draft.update(|draft| {
                            draft.coding.enabled = event_target_checked(&event);
                        })/>
                    "Project an explicitly selected design CDS onto this alignment"
                </label>
                <Show when=move || draft.get().coding.enabled>
                    <button class="text-tool" on:click=move |_| load_features.run(())>
                        "Load exact design CDS choices"
                    </button>
                    <label>
                        "Design CDS"
                        <select prop:value=move || draft.get().coding.feature
                            on:change=move |event| draft.update(|draft| {
                                draft.coding.feature = event_target_value(&event);
                            })>
                            <option value="">"Choose a CDS"</option>
                            {move || features.get().into_iter()
                                .map(|feature| view! {
                                    <option value=feature.index.to_string()>{feature.label}</option>
                                }).collect_view()}
                        </select>
                    </label>
                    <p class="helper">
                        "The menu shows at most the first 100 CDS annotations. Enter any exact annotation
                        number below; it is validated against the native source before comparison."
                    </p>
                    <label>
                        "Design annotation number (1-based)"
                        <input type="number" min="1"
                            prop:value={move || {
                                draft.get().coding.feature.parse::<usize>()
                                    .ok().and_then(|value| value.checked_add(1))
                                    .map(|value| value.to_string()).unwrap_or_default()
                            }}
                            on:input=move |event| draft.update(|draft| {
                                draft.coding.feature = event_target_value(&event).parse::<usize>()
                                    .ok().and_then(|value| value.checked_sub(1))
                                    .map(|value| value.to_string()).unwrap_or_default();
                            })/>
                    </label>
                    <label>
                        "Genetic code number"
                        <input type="number" min="1" max="255"
                            prop:value=move || draft.get().coding.genetic_code
                            on:input=move |event| draft.update(|draft| {
                                draft.coding.genetic_code = event_target_value(&event);
                            })/>
                    </label>
                    <label>
                        "Why this annotation is the intended design CDS"
                        <input maxlength="1024" prop:value=move || draft.get().coding.rationale
                            on:input=move |event| draft.update(|draft| {
                                draft.coding.rationale = event_target_value(&event);
                            })/>
                    </label>
                    <p class="helper">"The supplied sequence is interpreted conditionally through the chosen alignment. This does not create or verify a consensus annotation."</p>
                </Show>
            </details>
            <button class="button primary" disabled=move || busy.get() || context.get().is_err()
                on:click=move |_| compute.run(())>
                "Compare design and supplied consensus"
            </button>
            <button class="text-tool" disabled=move || !busy.get()
                on:click=move |_| cancel.run(())>
                "Cancel"
            </button>
            <p class="helper" role="status">{move || context.get().err().unwrap_or_else(|| status.get())}</p>
            {move || report.get().map(|shown| {
                let bound = StoredValue::new(shown);
                let pairwise_guard = Callback::new(move |review: Arc<Review>| {
                    let shown = bound.get_value();
                    if !Arc::ptr_eq(&review, &shown.review) {
                        return Err("The displayed comparison changed.".into());
                    }
                    guard.run(shown)
                });
                view! {
                    <PairwiseResultView review=bound.with_value(|shown| shown.review.clone())
                        guard=pairwise_guard on_save on_navigate/>
                }
            })}
        </details>
    }
}
