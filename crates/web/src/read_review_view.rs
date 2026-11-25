//! Human review of freshly replayed saved read evidence. A decision is a design
//! choice; it never changes the reference/read originals or proves sample identity.
use crate::{
    contigs::{ReadBinding, ReadNavigation},
    trace_validation::ReferenceNavigation,
    worker::CachedWorker,
};
use collagene_core::{
    protocol::{AnalysisJob, AnalysisResult, Request as WorkerRequest},
    read_checks,
    read_review::{self, Choice, Decision, EvidenceKey, Request},
    run_lineage::EventBinding,
    trace_validation::{ColumnKind, SiteState, TraceValidationReport},
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

struct Replayed {
    reference: EventBinding,
    reads: Vec<EventBinding>,
    report: TraceValidationReport,
    selected_context: EventBinding,
    epoch: u64,
}
#[derive(Clone, PartialEq, Eq)]
struct Nomination {
    value: String,
    label: String,
    choice: Choice,
}
struct Proposed {
    request: Request,
    preview: read_review::Preview,
}
fn label(key: &EvidenceKey) -> String {
    match key {
        EvidenceKey::Site { position } => format!("Reference base {}", position.saturating_add(1)),
        EvidenceKey::Insertion { index } => {
            format!("Insertion observation {}", index.saturating_add(1))
        }
    }
}
fn choice_value(choice: &Choice) -> String {
    match choice {
        Choice::RetainReference => "retain".into(),
        Choice::Unresolved => "unresolved".into(),
        Choice::UseBase {
            read_index,
            call_index,
            ..
        } => format!("base:{read_index}:{call_index}"),
        Choice::DeleteBase { read_index } => format!("delete:{read_index}"),
        Choice::UseInsertion { .. } => "insertion".into(),
    }
}
fn choice_label(choice: &Choice) -> String {
    match choice {
        Choice::RetainReference => "Retain reference".into(),
        Choice::Unresolved => "Leave unresolved".into(),
        Choice::UseBase {
            read_index,
            call_index,
            base,
        } => format!(
            "Use observed {} from read {} call {}",
            base,
            read_index + 1,
            call_index + 1
        ),
        Choice::DeleteBase { read_index } => {
            format!("Delete base, nominated read {}", read_index + 1)
        }
        Choice::UseInsertion {
            read_index,
            boundary,
            sequence,
        } => format!(
            "Use observed insertion {} at caret {} from read {}",
            preview_dna(sequence),
            boundary,
            read_index.saturating_add(1)
        ),
    }
}
fn default_reason(choice: &Choice) -> String {
    match choice {
        Choice::RetainReference => "Retained reference by explicit review".into(),
        Choice::Unresolved => "Left unresolved by explicit review".into(),
        Choice::UseBase {
            read_index,
            call_index,
            base,
        } => format!(
            "Selected observed {base} from read {}, call {}",
            read_index.saturating_add(1),
            call_index.saturating_add(1)
        ),
        Choice::DeleteBase { read_index } => format!(
            "Selected observed deletion from read {}",
            read_index.saturating_add(1)
        ),
        Choice::UseInsertion {
            read_index,
            boundary,
            ..
        } => format!(
            "Selected observed insertion from read {} at boundary {boundary}",
            read_index.saturating_add(1)
        ),
    }
}
fn preview_dna(dna: &str) -> String {
    let shown = dna.chars().take(200).collect::<String>();
    if shown.len() < dna.len() {
        format!("{shown}… ({} bases total)", dna.len())
    } else {
        shown
    }
}
fn nominations(report: &TraceValidationReport, key: &EvidenceKey) -> Vec<Nomination> {
    let mut result = vec![
        Nomination {
            value: "retain".into(),
            label: "Retain the original reference".into(),
            choice: Choice::RetainReference,
        },
        Nomination {
            value: "unresolved".into(),
            label: "Leave unresolved — no revised target".into(),
            choice: Choice::Unresolved,
        },
    ];
    match key {
        EvidenceKey::Site { position } => {
            for (read_index, read) in report.reads.iter().enumerate() {
                let mut deletion = false;
                for column in &read.columns {
                    if !column.eligible || column.reference_position != Some(*position) {
                        continue;
                    }
                    let choice = if column.kind == ColumnKind::Deletion {
                        if deletion {
                            continue;
                        }
                        deletion = true;
                        Some(Choice::DeleteBase { read_index })
                    } else if let (Some(call_index), Some(base)) =
                        (column.call_index, column.observed_base)
                    {
                        matches!(base, 'A' | 'C' | 'G' | 'T').then_some(Choice::UseBase {
                            read_index,
                            call_index,
                            base,
                        })
                    } else {
                        None
                    };
                    if let Some(choice) = choice {
                        let text = match choice {
                            Choice::UseBase { call_index, .. } => format!(
                                "Read {} · call {} · observed {}",
                                read_index + 1,
                                call_index + 1,
                                column.observed_base.unwrap_or('?')
                            ),
                            _ => format!("Read {} · eligible deletion", read_index + 1),
                        };
                        result.push(Nomination {
                            value: choice_value(&choice),
                            label: text,
                            choice,
                        });
                    }
                }
            }
        }
        EvidenceKey::Insertion { index } => {
            if let Some(insertion) = report.insertions.get(*index).filter(|i| i.eligible) {
                result.push(Nomination {
                    value: "insertion".into(),
                    label: "Use this eligible observed insertion".into(),
                    choice: Choice::UseInsertion {
                        read_index: insertion.read_index,
                        boundary: insertion.reference_boundary,
                        sequence: insertion.sequence.clone(),
                    },
                });
            }
        }
    }
    result
}
#[component]
pub fn ReadReviewPanel(
    attachment_hash: String,
    reference: Memo<EventBinding>,
    snapshot: Callback<(), Result<Workspace, String>>,
    epoch: RwSignal<u64>,
    on_commit: Callback<Request, Result<(), String>>,
    on_read: Callback<ReadNavigation>,
    on_reference: Callback<ReferenceNavigation, Result<(), String>>,
    event_available: Callback<EventBinding, bool>,
    attachment_available: Callback<String, bool>,
    #[prop(optional_no_strip)] initial_request: Option<Request>,
) -> impl IntoView {
    let initial_request =
        initial_request.filter(|request| request.attachment_hash == attachment_hash);
    let hash = StoredValue::new(attachment_hash);
    let loaded = RwSignal::new(None::<Arc<Replayed>>);
    let name = RwSignal::new(
        initial_request
            .as_ref()
            .map(|request| request.name.clone())
            .unwrap_or_else(|| "Reviewed reference".into()),
    );
    let decisions = RwSignal::new(
        initial_request
            .map(|request| request.decisions)
            .unwrap_or_default(),
    );
    let initialized = RwSignal::new(false);
    let selected = RwSignal::new(None::<EvidenceKey>);
    let chosen = RwSignal::new(String::new());
    let reason = RwSignal::new(String::new());
    let page = RwSignal::new(0usize);
    let edit_page = RwSignal::new(0usize);
    let preview = RwSignal::new(None::<Arc<Proposed>>);
    let status = RwSignal::new(
        "Load and freshly replay this saved check before choosing decisions.".to_string(),
    );
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
        let _ = (reference.get(), epoch.get());
        generation.update(|g| *g = g.wrapping_add(1));
        if let Some(w) = worker.get_value() {
            w.close();
        }
        worker.set_value(None);
        loaded.set(None);
        preview.set(None);
        pending.set(false);
        selected.set(None);
        chosen.set(String::new());
        reason.set(String::new());
        if initialized.get_untracked() {
            decisions.set(Vec::new());
        } else {
            initialized.set(true);
        }
        page.set(0);
    });
    Effect::new(move |_| {
        let _ = (name.get(), decisions.get());
        generation.update(|g| *g = g.wrapping_add(1));
        if let Some(active) = worker.get_value() {
            active.close();
        }
        worker.set_value(None);
        pending.set(false);
        preview.set(None);
        edit_page.set(0);
    });
    let current_request = Callback::new(move |()| Request {
        attachment_hash: hash.get_value(),
        name: name.get_untracked(),
        decisions: decisions.get_untracked(),
    });
    let guard = Callback::new(move |bound: Arc<Replayed>| -> Result<(), String> {
        if bound.epoch != epoch.get_untracked()
            || bound.selected_context != reference.get_untracked()
        {
            return Err(
                "The pinned construct changed. Replay this check again before reviewing it.".into(),
            );
        }
        if !attachment_available.run(hash.get_value())
            || !event_available.run(bound.reference.clone())
            || bound
                .reads
                .iter()
                .any(|read| !event_available.run(read.clone()))
        {
            return Err("The saved check or an exact original source event is unavailable. No decision was applied.".into());
        }
        Ok(())
    });
    let replay_alive = alive.clone();
    let replay = move |_| {
        let prepared = (|| -> Result<_, String> {
            let frame = snapshot.run(())?;
            let selected_context = reference.get_untracked();
            let attachment = hash.get_value();
            let saved = read_checks::load(&frame, &attachment).map_err(|e| e.to_string())?;
            if saved.reference.document_id != selected_context.document_id {
                return Err("This saved check belongs to a different construct document.".into());
            }
            let revision =
                read_checks::replay_revision(&frame, &attachment).map_err(|e| e.to_string())?;
            let w = Rc::new(CachedWorker::new()?);
            Ok((frame, attachment, selected_context, revision, w))
        })();
        let (frame, attachment, context, revision, w) = match prepared {
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
        let captured_epoch = epoch.get_untracked();
        pending.set(true);
        loaded.set(None);
        preview.set(None);
        status.set(
            "Replaying the exact retained reference and read events with the current engine…"
                .into(),
        );
        let alive = replay_alive.clone();
        spawn_local(async move {
            let response = w
                .ask(WorkerRequest {
                    id: token,
                    revision: revision.clone(),
                    job: AnalysisJob::ReplayReadCheck {
                        workspace: frame,
                        attachment_hash: attachment,
                    },
                })
                .await;
            w.close();
            if !alive.load(Ordering::Relaxed)
                || generation.get_untracked() != token
                || epoch.get_untracked() != captured_epoch
                || reference.get_untracked() != context
            {
                return;
            }
            pending.set(false);
            worker.set_value(None);
            let outcome = (|| -> Result<_, String> {
                let response = response?;
                if response.id != token || response.revision != revision {
                    return Err("Read replay response identity changed.".into());
                }
                let checked = match response.result {
                    AnalysisResult::ReadCheckReplayed(checked) => checked,
                    AnalysisResult::Error(error) => return Err(error),
                    _ => return Err("Unexpected read replay response.".into()),
                };
                checked
                    .validate_binding(&revision)
                    .map_err(|e| e.to_string())?;
                if checked.report.engine_version != collagene_core::verification::ENGINE_VERSION
                    || checked.report.reference_revision != checked.saved.reference.revision
                    || checked.report.read_revisions
                        != checked
                            .saved
                            .reads
                            .iter()
                            .map(|read| read.revision.clone())
                            .collect::<Vec<_>>()
                {
                    return Err("Fresh replay has different source identities or engine.".into());
                }
                let bound = Arc::new(Replayed {
                    reference: checked.saved.reference,
                    reads: checked.saved.reads,
                    report: checked.report,
                    selected_context: context,
                    epoch: captured_epoch,
                });
                guard.run(bound.clone())?;
                Ok(bound)
            })();
            match outcome {
                Ok(bound) => {
                    loaded.set(Some(bound));
                    status.set("Fresh evidence ready. Choose each decision explicitly and optionally add a review note. Conflicting observations remain retained.".into());
                }
                Err(error) => status.set(error),
            }
        });
    };
    let choices = Memo::new(move |_| match (loaded.get(), selected.get()) {
        (Some(bound), Some(key)) => nominations(&bound.report, &key),
        _ => Vec::new(),
    });
    let choose_evidence = Callback::new(move |key: EvidenceKey| {
        let Some(bound) = loaded.get_untracked() else {
            return;
        };
        if let Err(error) = guard.run(bound) {
            status.set(error);
            return;
        }
        let existing =
            decisions.with_untracked(|items| items.iter().find(|d| d.key == key).cloned());
        selected.set(Some(key));
        chosen.set(
            existing
                .as_ref()
                .map(|d| choice_value(&d.choice))
                .unwrap_or_default(),
        );
        reason.set(
            existing
                .map(|d| {
                    if d.reason == default_reason(&d.choice) {
                        String::new()
                    } else {
                        d.reason
                    }
                })
                .unwrap_or_default(),
        );
    });
    let add = move |_| {
        let outcome = (|| -> Result<(), String> {
            let bound = loaded
                .get_untracked()
                .ok_or("Freshly replay the check first.")?;
            guard.run(bound)?;
            let key = selected
                .get_untracked()
                .ok_or("Choose a site or insertion observation.")?;
            let value = chosen.get_untracked();
            let choice = choices
                .with_untracked(|items| {
                    items
                        .iter()
                        .find(|n| n.value == value)
                        .map(|n| n.choice.clone())
                })
                .ok_or("Choose an available observed alternative explicitly.")?;
            let note = reason.get_untracked();
            if note.len() > 1024 {
                return Err("Keep the optional review note within 1024 bytes.".into());
            }
            let why = if note.trim().is_empty() {
                default_reason(&choice)
            } else {
                note.trim().to_string()
            };
            if decisions
                .with_untracked(|items| items.len() >= 256 && !items.iter().any(|d| d.key == key))
            {
                return Err("A review supports at most 256 explicit decisions.".into());
            }
            decisions.update(|items| {
                let decision = Decision {
                    key: key.clone(),
                    choice,
                    reason: why,
                };
                if let Some(old) = items.iter_mut().find(|d| d.key == key) {
                    *old = decision;
                } else {
                    items.push(decision);
                }
            });
            Ok(())
        })();
        match outcome {
            Ok(()) => status.set(
                "Decision recorded in this draft. Preview the complete decision set before saving."
                    .into(),
            ),
            Err(error) => status.set(error),
        }
    };
    let preview_alive = alive.clone();
    let run_preview = Callback::new(move |()| {
        let prepared = (|| -> Result<_, String> {
            let bound = loaded
                .get_untracked()
                .ok_or("Freshly replay the saved check first.")?;
            guard.run(bound.clone())?;
            let request = current_request.run(());
            let frame = snapshot.run(())?;
            read_review::source_preflight(&frame, &request).map_err(|e| e.to_string())?;
            let revision =
                read_review::preview_revision(&frame, &request).map_err(|e| e.to_string())?;
            let w = Rc::new(CachedWorker::new()?);
            Ok((bound, request, frame, revision, w))
        })();
        let (bound, request, frame, revision, w) = match prepared {
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
        preview.set(None);
        status.set("Replaying sources and checking every nominated observation before constructing a review preview…".into());
        let alive = preview_alive.clone();
        spawn_local(async move {
            let response = w
                .ask(WorkerRequest {
                    id: token,
                    revision: revision.clone(),
                    job: AnalysisJob::PreviewReadReview {
                        workspace: frame,
                        request: request.clone(),
                    },
                })
                .await;
            w.close();
            if !alive.load(Ordering::Relaxed) || generation.get_untracked() != token {
                return;
            }
            pending.set(false);
            worker.set_value(None);
            let outcome = (|| -> Result<_, String> {
                guard.run(bound.clone())?;
                if current_request.run(()) != request {
                    return Err("Review decisions changed. Preview them again.".into());
                }
                let response = response?;
                if response.id != token || response.revision != revision {
                    return Err("Review preview response identity changed.".into());
                }
                let computed = match response.result {
                    AnalysisResult::ReadReviewPreview(preview) => preview,
                    AnalysisResult::Error(error) => return Err(error),
                    _ => return Err("Unexpected review preview response.".into()),
                };
                computed
                    .validate_binding(&revision)
                    .map_err(|e| e.to_string())?;
                if computed.request != request
                    || computed.reference != bound.reference
                    || computed.reads != bound.reads
                {
                    return Err("Review preview has different source events or choices.".into());
                }
                Ok(Arc::new(Proposed {
                    request,
                    preview: computed,
                }))
            })();
            match outcome {
                Ok(value) => {
                    preview.set(Some(value));
                    status.set("Review preview ready. This is a human interpretation, not a new measurement or a confirmed genotype.".into());
                }
                Err(error) => status.set(error),
            }
        });
    });
    let commit = move |_| {
        let outcome = (|| -> Result<(), String> {
            let bound = loaded.get_untracked().ok_or("Replay the check first.")?;
            guard.run(bound)?;
            let proposed = preview
                .get_untracked()
                .ok_or("Preview the current decisions first.")?;
            if proposed.request != current_request.run(()) {
                return Err("Decisions changed; preview them again before saving.".into());
            }
            on_commit.run(proposed.request.clone())
        })();
        if let Err(error) = outcome {
            status.set(error);
        }
    };
    let navigate_reference = Callback::new(move |range: (usize, usize)| {
        let outcome = (|| -> Result<(), String> {
            let bound = loaded.get_untracked().ok_or("Replay the check first.")?;
            guard.run(bound.clone())?;
            on_reference.run(ReferenceNavigation {
                reference: bound.reference.clone(),
                selection: Some(range),
                epoch: bound.epoch,
            })
        })();
        if let Err(error) = outcome {
            status.set(error);
        }
    });
    let navigate_call = Callback::new(move |(index, call_index): (usize, usize)| {
        let outcome = (|| -> Result<(), String> {
            let bound = loaded.get_untracked().ok_or("Replay the check first.")?;
            guard.run(bound.clone())?;
            let source = bound
                .reads
                .get(index)
                .ok_or("The nominated read event is unavailable.")?;
            on_read.run(ReadNavigation {
                binding: ReadBinding {
                    document_id: source.document_id.clone(),
                    event_id: source.event_id.clone(),
                    revision: source.revision.clone(),
                },
                call_index,
                epoch: bound.epoch,
            });
            Ok(())
        })();
        if let Err(error) = outcome {
            status.set(error);
        }
    });
    view! {<section class="read-review-panel"><h2>"Review saved read evidence"</h2>
        <p class="helper">"Fresh replay is required before choosing observations. Reference and read originals remain immutable. Your decisions express an interpretation; conflicting alleles and unknown sample identity are not resolved by choosing one call."</p>
        <button class="button primary" disabled=move||pending.get() on:click=replay>"Load and replay saved check"</button>
        <button class="text-tool" disabled=move||!pending.get() on:click=move|_|{generation.update(|g|*g=g.wrapping_add(1));if let Some(w)=worker.get_value(){w.close();}worker.set_value(None);pending.set(false);preview.set(None);status.set("Review computation cancelled; source evidence and saved records are unchanged.".into());}>"Cancel review computation"</button>
        <p role="status">{move||status.get()}</p>
        <Show when=move||loaded.get().is_some()>
            <label>"Revised reference name"<input maxlength="256" aria-label="Read review target name" prop:value=move||name.get() on:input=move|e|name.set(event_target_value(&e))/></label>
            {move||loaded.get().map(|bound|{let site_indices=bound.report.sites.iter().enumerate().filter(|(_,s)|s.state!=SiteState::SupportedMatch).map(|(i,_)|i).collect::<Vec<_>>();let total=site_indices.len()+bound.report.insertions.len();let pages=total.max(1).div_ceil(50);let start=(page.get().min(pages-1))*50;view!{<h3>"Observed differences and unresolved evidence"</h3><p>{format!("{} observations · page {} of {}",total,page.get().min(pages-1)+1,pages)}</p>
                {(start..(start+50).min(total)).map(|ordinal|{if ordinal<site_indices.len(){let site=&bound.report.sites[site_indices[ordinal]];let position=site.position;view!{<article class="read-observation"><button class="text-tool" on:click=move|_|choose_evidence.run(EvidenceKey::Site{position})>{format!("Base {} · {:?} · reference {} · consensus {}",position+1,site.state,site.reference_base,site.consensus_base)}</button><button class="text-tool" on:click=move|_|navigate_reference.run((position,position+1))>"Show reference base"</button><p class="helper">{format!("Eligible files {:?}; low-quality {}; ambiguous {}; eligible deletion observations {}.",site.eligible_read_indices,site.low_quality_calls,site.ambiguous_calls,site.eligible_deletions)}</p>{(site.state==SiteState::ConflictingCalls).then(||view!{<p class="helper">"Conflicting alleles remain in the source evidence even if you nominate one observed call."</p>})}</article>}.into_any()}else{let index=ordinal-site_indices.len();let insertion=&bound.report.insertions[index];let boundary=insertion.reference_boundary;view!{<article class="read-observation"><button class="text-tool" on:click=move|_|choose_evidence.run(EvidenceKey::Insertion{index})>{format!("Insertion {} · caret {} · read {} · eligible {}",index+1,boundary,insertion.read_index+1,insertion.eligible)}</button><pre>{preview_dna(&insertion.sequence)}</pre><button class="text-tool" on:click=move|_|navigate_reference.run((boundary,boundary))>"Show reference insertion caret"</button></article>}.into_any()}}).collect_view()}
                <button class="text-tool" disabled=move||page.get()==0 on:click=move|_|page.update(|p|*p=p.saturating_sub(1))>"Previous observations"</button><button class="text-tool" disabled=move||{page.get()+1>=pages} on:click=move|_|page.update(|p|*p=(*p+1).min(pages-1))>"Next observations"</button>
            }})}
            {move||selected.get().map(|key|view!{<section class="read-review-decision"><h3>{label(&key)}</h3><label>"Decision"<select aria-label="Read review decision" prop:value=move||chosen.get() on:change=move|e|chosen.set(event_target_value(&e))><option value="">"Choose an interpretation…"</option>{move||choices.get().into_iter().map(|choice|view!{<option value=choice.value>{choice.label}</option>}).collect_view()}</select></label><p class="helper">"Only exact eligible observed calls/deletions can be nominated. Retain reference and unresolved are separate human choices. A deletion has no invented read call."</p>
                <button class="text-tool" disabled=move||!choices.with(|items|items.iter().any(|n|n.value==chosen.get()&&matches!(n.choice,Choice::UseBase{..}))) on:click=move|_|{if let Some(Choice::UseBase{read_index,call_index,..})=choices.with_untracked(|items|items.iter().find(|n|n.value==chosen.get_untracked()).map(|n|n.choice.clone())){navigate_call.run((read_index,call_index));}}>"Show nominated original read call"</button>
                <label>"Optional review note"<textarea maxlength="1024" aria-label="Optional read review note" prop:value=move||reason.get() on:input=move|e|reason.set(event_target_value(&e))></textarea></label><button class="button secondary" on:click=add>"Add / update decision"</button>
            </section>})}
            <h3>{move||format!("{} explicit decisions",decisions.with(Vec::len))}</h3>
            {move||decisions.get().into_iter().map(|decision|{let key=StoredValue::new(decision.key.clone());view!{<article class="read-observation"><p>{format!("{} · {} · {}",label(&decision.key),choice_label(&decision.choice),decision.reason)}</p><button class="text-tool" on:click=move|_|choose_evidence.run(key.get_value())>"Edit decision"</button><button class="text-tool" on:click=move|_|decisions.update(|items|items.retain(|d|d.key!=key.get_value()))>"Remove decision"</button></article>}}).collect_view()}
            <button class="button primary" disabled=move||pending.get() on:click=move|_|run_preview.run(())>"Preview review decisions"</button>
            {move||preview.get().map(|proposed|{let pages=proposed.preview.edits.len().max(1).div_ceil(50);view!{<section class="read-review-preview"><h3>"Literal edit preview"</h3><p class="helper">"Coordinates below are zero-based original-reference intervals [start,end). Insertions have equal boundaries. A full revised reference is withheld when decisions remain unresolved."</p>{proposed.preview.edits.iter().skip(edit_page.get()*50).take(50).map(|edit|view!{<article><p>{format!("[{}, {}) → {} · {}",edit.start,edit.end,preview_dna(&edit.replacement),edit.reason)}</p></article>}).collect_view()}<button class="text-tool" disabled=move||edit_page.get()==0 on:click=move|_|edit_page.update(|p|*p=p.saturating_sub(1))>"Previous edits"</button><button class="text-tool" disabled=move||{edit_page.get()+1>=pages} on:click=move|_|edit_page.update(|p|*p=(*p+1).min(pages-1))>"Next edits"</button>
                {proposed.preview.target.as_ref().map(|target|view!{<p>{format!("{} · {} bases · {:?}",target.name,target.sequence.len(),target.topology)}</p><pre>{preview_dna(target.sequence.as_str())}</pre>})}
                {proposed.preview.warnings.iter().map(|warning|view!{<p class="helper">{warning.clone()}</p>}).collect_view()}
                <p class="helper">"DNA displays are limited to 200 bases; source originals and the saved full decision/evidence record are retained. A revised reference is an interpretation, not a confirmed genotype or physical construct."</p>
                <button class="button primary" disabled=move||pending.get() on:click=commit>{if proposed.preview.target.is_some(){"Save decisions and create revised reference"}else{"Save unresolved decisions"}}</button>
            </section>}})}
        </Show>
    </section>}
}
