//! Source-event-bound read checks: durable capture and explicitly fresh replay.
use crate::{
    contigs::{ReadBinding, ReadNavigation},
    worker::CachedWorker,
};
use collagene_core::{
    protocol::{AnalysisJob, AnalysisResult, Request},
    read_checks,
    run_lineage::{self, EventBinding},
    trace_validation::{SiteState, TraceValidationReport},
    workspace::Workspace,
};
use leptos::prelude::*;
use std::{
    cell::{Cell, RefCell},
    rc::Rc,
    sync::{
        Arc,
        atomic::{AtomicBool, Ordering},
    },
};
use wasm_bindgen_futures::spawn_local;
fn fingerprint(w: &Workspace) -> Result<String, String> {
    run_lineage::workspace_fingerprint(w).map_err(|e| e.to_string())
}
#[component]
pub fn ReadCheckCapturePanel(
    request: crate::trace_validation::RecordRequest,
    snapshot: Callback<(), Result<Workspace, String>>,
    validate: Callback<Workspace, Result<(), String>>,
    on_restore: Callback<Workspace>,
    db: StoredValue<Option<crate::storage::Database>, LocalStorage>,
    busy: RwSignal<bool>,
    loaded: RwSignal<bool>,
    save_generation: RwSignal<u64>,
    epoch: RwSignal<u64>,
    message: RwSignal<String>,
) -> impl IntoView {
    let captured_epoch = request.epoch;
    let request = StoredValue::new(request);
    let status = RwSignal::new(String::new());
    let completed = RwSignal::new(false);
    let alive = Arc::new(AtomicBool::new(true));
    let worker = StoredValue::new_local(None::<Rc<CachedWorker>>);
    let transaction =
        StoredValue::new_local(Rc::new(RefCell::new(None::<crate::storage::Transaction>)));
    let generation = StoredValue::new_local(Rc::new(Cell::new(0u64)));
    let owns_lock = StoredValue::new_local(Rc::new(Cell::new(false)));
    let cancel = Callback::new(move |()| {
        generation.with_value(|g| g.set(g.get().wrapping_add(1)));
        if let Some(w) = worker.get_value() {
            w.close();
        }
        if let Some(tx) = transaction.with_value(|t| t.borrow().clone()) {
            let _ = tx.abort();
        }
        status.set("Read-check capture cancelled. No partial check was installed.".into());
    });
    let cleanup = alive.clone();
    on_cleanup(move || {
        cleanup.store(false, Ordering::Relaxed);
        generation.with_value(|g| g.set(g.get().wrapping_add(1)));
        if let Some(w) = worker.get_value() {
            w.close();
        }
        if let Some(tx) = transaction.with_value(|t| t.borrow().clone()) {
            let _ = tx.abort();
        } else if owns_lock.with_value(|o| o.replace(false)) {
            busy.set(false);
        }
    });
    Effect::new(move |_| {
        if epoch.get() != captured_epoch && owns_lock.with_value(|o| o.get()) {
            cancel.run(());
        }
    });
    let run = {
        let alive = alive.clone();
        move || {
            if busy.get_untracked() || completed.get_untracked() {
                return;
            }
            let prepared = (|| -> Result<_, String> {
                if epoch.get_untracked() != captured_epoch {
                    return Err(
                        "Workspace changed; validate the selected reference and reads again."
                            .into(),
                    );
                }
                let database = db
                    .get_value()
                    .ok_or("Recovery storage is unavailable. No check was saved.")?;
                let frame = snapshot.run(())?;
                let original = fingerprint(&frame)?;
                let r = request.get_value();
                let revision =
                    read_checks::capture_revision(&frame, &r.reference, &r.reads, &r.options)
                        .map_err(|e| e.to_string())?;
                let w = Rc::new(CachedWorker::new()?);
                Ok((database, frame, original, r, revision, w))
            })();
            let (database, frame, original, r, revision, w) = match prepared {
                Ok(p) => p,
                Err(e) => {
                    status.set(e);
                    return;
                }
            };
            worker.set_value(Some(w.clone()));
            let generation = generation.get_value();
            let owns_lock = owns_lock.get_value();
            let transaction = transaction.get_value();
            generation.set(generation.get().wrapping_add(1));
            let token = generation.get();
            owns_lock.set(true);
            busy.set(true);
            save_generation.update(|g| *g = g.wrapping_add(1));
            status.set(
                "Reconstructing read evidence against exact retained reference/read events…".into(),
            );
            let alive = alive.clone();
            spawn_local(async move {
                let result = async {
                    let response = w
                        .ask(Request {
                            id: u32::MAX - 230,
                            revision: revision.clone(),
                            job: AnalysisJob::CaptureReadCheck {
                                workspace: frame.clone(),
                                reference: r.reference.clone(),
                                reads: r.reads.clone(),
                                options: r.options.clone(),
                            },
                        })
                        .await?;
                    if response.id != u32::MAX - 230 || response.revision != revision {
                        return Err("Read-check response belongs to another request.".into());
                    }
                    let commit = match response.result {
                        AnalysisResult::ReadCheckCaptured(c) => c,
                        AnalysisResult::Error(e) => return Err(e),
                        _ => return Err("Unexpected read-check capture result.".into()),
                    };
                    commit
                        .validate_binding(&revision)
                        .map_err(|e| e.to_string())?;
                    let guard = || -> Result<(), String> {
                        if !alive.load(Ordering::Relaxed)
                            || epoch.get_untracked() != captured_epoch
                            || generation.get() != token
                        {
                            return Err("Capture cancelled or context changed.".into());
                        }
                        if fingerprint(&snapshot.run(())?)? != original {
                            return Err(
                                "Workspace changed during read validation. No check was installed."
                                    .into(),
                            );
                        }
                        Ok(())
                    };
                    guard()?;
                    validate.run(commit.workspace.clone())?;
                    let current =
                        crate::app::encode_workspace(commit.workspace.clone(), u32::MAX - 231)
                            .await?;
                    let previous = crate::app::encode_workspace(frame, u32::MAX - 232).await?;
                    guard()?;
                    crate::storage::replace_workspace_cancellable(
                        &database,
                        &current,
                        &previous,
                        transaction.clone(),
                    )
                    .await?;
                    Ok::<_, String>(commit.workspace)
                }
                .await;
                match result {
                    Ok(frame) => {
                        // Durable completion is the commit point; restore even if the modal closed immediately afterwards.
                        loaded.set(true);
                        message.set("Read check saved with exact reference/read source events. Stored checks are historical until explicitly replayed; no experimental sample confirmation is implied.".into());
                        if alive.load(Ordering::Relaxed) {
                            completed.set(true);
                            status.set("Read check saved.".into());
                        }
                        on_restore.run(frame);
                    }
                    Err(e) => {
                        if alive.load(Ordering::Relaxed) {
                            status.set(format!("Read check was not saved. {e}"));
                        }
                    }
                }
                if owns_lock.replace(false) {
                    busy.set(false);
                }
                w.close();
            });
        }
    };
    let run = StoredValue::new_local(run);
    view! {<section class="read-check-capture"><h2>"Save read check"</h2><p>"Recompute against exact retained reference and read events, then save the check in one recovery transaction. Original reads and document heads remain unchanged."</p><p>{request.with_value(|r|format!("{} reads · reference event {}",r.reads.len(),r.reference.event_id))}</p><p role="status">{move||status.get()}</p><button class="button primary" disabled=move||busy.get()||completed.get()||epoch.get()!=captured_epoch on:click=move|_|run.with_value(|r|r())>"Save read check"</button><button class="button secondary" disabled=move||!busy.get()||!owns_lock.with_value(|o|o.get()) on:click=move|_|cancel.run(())>"Cancel saving check"</button></section>}
}
#[derive(Clone)]
struct Review {
    hash: String,
    reference: EventBinding,
    reads: Vec<EventBinding>,
    report: Arc<TraceValidationReport>,
    fresh: bool,
    epoch: u64,
    selected_context: EventBinding,
}
fn historical(
    frame: &Workspace,
    hash: String,
    epoch: u64,
    selected_context: EventBinding,
) -> Result<Review, String> {
    let saved = read_checks::load(frame, &hash).map_err(|e| e.to_string())?;
    Ok(Review {
        hash,
        reference: saved.reference,
        reads: saved.reads,
        report: Arc::new(saved.report),
        fresh: false,
        epoch,
        selected_context,
    })
}
fn read_navigation(binding: &EventBinding, call_index: usize, epoch: u64) -> ReadNavigation {
    ReadNavigation {
        binding: ReadBinding {
            document_id: binding.document_id.clone(),
            event_id: binding.event_id.clone(),
            revision: binding.revision.clone(),
        },
        call_index,
        epoch,
    }
}
#[component]
pub fn SavedReadChecksPanel(
    reference: Memo<EventBinding>,
    snapshot: Callback<(), Result<Workspace, String>>,
    epoch: RwSignal<u64>,
    on_read: Callback<ReadNavigation>,
    on_reference: Callback<crate::trace_validation::ReferenceNavigation, Result<(), String>>,
    event_available: Callback<EventBinding, bool>,
    attachment_available: Callback<String, bool>,
    on_review: Callback<String, Result<(), String>>,
) -> impl IntoView {
    let action_guard = Callback::new(move |r: Arc<Review>| -> Result<(), String> {
        if r.epoch != epoch.get_untracked()
            || r.selected_context != reference.get_untracked()
            || !attachment_available.run(r.hash.clone())
            || !event_available.run(r.reference.clone())
            || r.reads
                .iter()
                .any(|binding| !event_available.run(binding.clone()))
        {
            return Err("The pinned reference, saved attachment or exact source events are no longer available. Reopen the check.".into());
        }
        Ok(())
    });
    let listed = RwSignal::new(Vec::<read_checks::Summary>::new());
    let review = RwSignal::new(None::<Arc<Review>>);
    let status = RwSignal::new(String::new());
    let pending = RwSignal::new(false);
    let page = RwSignal::new(0usize);
    let site_page = RwSignal::new(0usize);
    let column_page = RwSignal::new(0usize);
    let read_index = RwSignal::new(0usize);
    let worker = StoredValue::new_local(None::<Rc<CachedWorker>>);
    let generation = StoredValue::new_local(Rc::new(Cell::new(0u64)));
    let alive = Arc::new(AtomicBool::new(true));
    let cleanup = alive.clone();
    on_cleanup(move || {
        cleanup.store(false, Ordering::Relaxed);
        generation.with_value(|g| g.set(g.get().wrapping_add(1)));
        if let Some(w) = worker.get_value() {
            w.close();
        }
    });
    Effect::new(move |_| {
        let binding = reference.get();
        let _ = epoch.get();
        generation.with_value(|g| g.set(g.get().wrapping_add(1)));
        if let Some(w) = worker.get_value() {
            w.close();
        }
        review.set(None);
        pending.set(false);
        page.set(0);
        site_page.set(0);
        column_page.set(0);
        read_index.set(0);
        match snapshot
            .run(())
            .and_then(|frame| read_checks::list(&frame).map_err(|e| e.to_string()))
        {
            Ok(items) => {
                listed.set(
                    items
                        .into_iter()
                        .filter(|s| s.reference.document_id == binding.document_id)
                        .collect(),
                );
                status.set(String::new());
            }
            Err(e) => {
                listed.set(Vec::new());
                status.set(e);
            }
        }
    });
    let open = Callback::new(move |hash: String| {
        generation.with_value(|g| g.set(g.get().wrapping_add(1)));
        if let Some(w) = worker.get_value() {
            w.close();
        }
        pending.set(false);
        site_page.set(0);
        column_page.set(0);
        read_index.set(0);
        match snapshot.run(()).and_then(|frame| {
            historical(
                &frame,
                hash,
                epoch.get_untracked(),
                reference.get_untracked(),
            )
        }) {
            Ok(bound) => {
                review.set(Some(Arc::new(bound)));
                status.set(
                    "Historical record opened. Replay explicitly for current-engine evidence."
                        .into(),
                );
            }
            Err(e) => status.set(e),
        }
    });
    let cancel = Callback::new(move |()| {
        generation.with_value(|g| g.set(g.get().wrapping_add(1)));
        if let Some(w) = worker.get_value() {
            w.close();
        }
        pending.set(false);
        status.set("Replay cancelled. The stored historical record is unchanged.".into());
    });
    let replay = {
        let alive = alive.clone();
        move |_| {
            let prepared = (|| -> Result<_, String> {
                let old = review.get_untracked().ok_or("Open a saved check first.")?;
                let frame = snapshot.run(())?;
                let fp = fingerprint(&frame)?;
                action_guard.run(old.clone())?;
                let revision =
                    read_checks::replay_revision(&frame, &old.hash).map_err(|e| e.to_string())?;
                let w = Rc::new(CachedWorker::new()?);
                Ok((old, frame, fp, revision, w))
            })();
            let (old, frame, fp, revision, w) = match prepared {
                Ok(p) => p,
                Err(e) => {
                    status.set(e);
                    return;
                }
            };
            if let Some(previous) = worker.get_value() {
                previous.close();
            }
            worker.set_value(Some(w.clone()));
            let generation = generation.get_value();
            generation.set(generation.get().wrapping_add(1));
            let token = generation.get();
            let bound_reference = reference.get_untracked();
            pending.set(true);
            let alive = alive.clone();
            status.set("Replaying exact historical source events with the current engine…".into());
            spawn_local(async move {
                let result=async{
 let response=w.ask(Request{id:u32::MAX-233,revision:revision.clone(),job:AnalysisJob::ReplayReadCheck{workspace:frame,attachment_hash:old.hash.clone()}}).await?;
 if response.id!=u32::MAX-233||response.revision!=revision{return Err("Replay response belongs to another request.".into())}
 let checked=match response.result{AnalysisResult::ReadCheckReplayed(c)=>c,AnalysisResult::Error(e)=>return Err(e),_=>return Err("Unexpected saved read-check result.".into())};checked.validate_binding(&revision).map_err(|e|e.to_string())?;
 if checked.saved.reference!=old.reference||checked.saved.reads!=old.reads||checked.report.reference_revision!=old.reference.revision||checked.report.read_revisions!=old.reads.iter().map(|r|r.revision.clone()).collect::<Vec<_>>()||checked.report.engine_version!=collagene_core::verification::ENGINE_VERSION{return Err("Replay source bindings or current engine differ from the selected check.".into())}
 Ok::<_,String>(checked)
 }.await;
                if !alive.load(Ordering::Relaxed) || generation.get() != token {
                    return;
                }
                pending.set(false);
                if epoch.get_untracked() != old.epoch
                    || reference.get_untracked() != bound_reference
                    || snapshot
                        .run(())
                        .and_then(|f| fingerprint(&f))
                        .ok()
                        .as_deref()
                        != Some(fp.as_str())
                {
                    status.set(
                        "Workspace or selected reference changed; replay result discarded.".into(),
                    );
                    w.close();
                    return;
                }
                match result {
                    Ok(checked) => {
                        review.set(Some(Arc::new(Review {
                            hash: old.hash.clone(),
                            reference: checked.saved.reference,
                            reads: checked.saved.reads,
                            report: Arc::new(checked.report),
                            fresh: true,
                            epoch: old.epoch,
                            selected_context: bound_reference,
                        })));
                        status.set("Fresh current-engine evidence under the recorded options. This does not establish experimental sample identity or physical chemistry.".into());
                    }
                    Err(e) => status.set(e),
                }
                w.close();
            });
        }
    };
    let current = move || {
        review.get().filter(|r| {
            r.epoch == epoch.get()
                && r.selected_context == reference.get()
                && action_guard.run(r.clone()).is_ok()
        })
    };
    let download = Callback::new(move |()| {
        let Some(r) = review.get_untracked() else {
            return;
        };
        let result = (|| {
            action_guard.run(r.clone())?;
            let bytes=serde_json::to_vec_pretty(&serde_json::json!({"attachment_hash":r.hash,"reference":r.reference,"reads":r.reads,"fresh_replay":r.fresh,"scope":"read evidence under recorded thresholds; not experimental confirmation","report":r.report.as_ref()})).map_err(|e|e.to_string())?;
            if bytes.len() > 32 << 20 {
                return Err("JSON exceeds the32MiB download limit.".into());
            }
            crate::app::prepare_download("read-check.json", "application/json", &bytes)
        })();
        match result {
            Ok(artifact) => {
                let check = move || {
                    if status.try_get_untracked().is_none() {
                        return Err("Export review closed.".into());
                    }
                    if !review.with_untracked(|current| {
                        current.as_ref().is_some_and(|v| Arc::ptr_eq(v, &r))
                    }) {
                        return Err("Selected read-check evidence changed.".into());
                    }
                    action_guard.run(r.clone())
                };
                let completion = check.clone();
                artifact.complete_guarded(check, move |result| {
                    if completion().is_err() {
                        return;
                    }
                    status.set(match result {
                        Ok(()) => "Read-check evidence saved.".into(),
                        Err(error) => error,
                    });
                });
            }
            Err(error) => status.set(error),
        }
    });
    view! {<section class="saved-read-checks"><h2>"Saved read checks"</h2><p class="helper">"Checks retain exact historical reference and read events. Opening a record does not restore a current verification result."</p><p role="status">{move||status.get()}</p>
    <div class="read-check-list">{move||listed.get().into_iter().skip(page.get()*20).take(20).map(|s|{let hash=s.attachment_hash.clone();view!{<article><p>{format!("{} reads · recorded {:?} · engine {} · reference event {}",s.read_count,s.recorded_status,s.recorded_engine,s.reference.event_id)}</p><button class="button secondary" on:click=move|_|open.run(hash.clone())>"Open historical check"</button></article>}}).collect_view()}</div>
    <button disabled=move||page.get()==0 on:click=move|_|page.update(|p|*p=p.saturating_sub(1))>"Previous checks"</button><button disabled={move||(page.get()+1)*20>=listed.with(Vec::len)} on:click=move|_|page.update(|p|*p+=1)>"Next checks"</button>
    {move||current().map(|r|{let target=r.reference.clone();let report=r.report.clone();let context=r.clone();view!{<section class="read-check-review"><h3>{if r.fresh{"Fresh replay"}else{"Historical evidence"}}</h3><p>{format!("{:?} · eligible coverage {}/{} · source engine {}",report.status,report.eligible_coverage,report.consensus.len(),report.engine_version)}</p><p class="helper">"Uncovered bases, conflicts, low-quality calls and unresolved placements remain separate evidence. Native caller scores are not independently calibrated probabilities."</p><button class="button secondary" on:click=move|_|{let result=action_guard.run(context.clone()).and_then(|()|on_reference.run(crate::trace_validation::ReferenceNavigation{reference:target.clone(),selection:None,epoch:context.epoch}));if let Err(e)=result{status.set(e)}}>"Open exact reference event"</button><ul>{report.warnings.iter().take(32).map(|w|view!{<li>{w.clone()}</li>}).collect_view()}</ul></section>}})}
    <button class="button primary" disabled=move||pending.get()||current().is_none() on:click=replay>"Replay saved read check"</button><button class="button secondary" disabled=move||!pending.get() on:click=move|_|cancel.run(())>"Cancel replay"</button><button class="button secondary" disabled=move||current().is_none() on:click=move|_|download.run(())>"Download read-check evidence"</button>
    <button class="button primary" disabled=move||pending.get()||current().is_none() on:click=move|_|{if let Some(record)=current(){if let Err(error)=action_guard.run(record.clone()).and_then(|()|on_review.run(record.hash.clone())){status.set(error)}}}>"Review differences and propose a revision"</button>
    <h3>"Discrepancies and unresolved sites"</h3><p class="helper">"25 sites per page; exact supporting call columns are paged below. Supported matches remain available in the full report."</p>
    {move||current().map(|r|r.report.sites.iter().filter(|s|s.state!=SiteState::SupportedMatch).skip(site_page.get()*25).take(25).map(|s|{
        let context=r.clone();let position=s.position;let ref_binding=r.reference.clone();let ref_epoch=r.epoch;let reference_context=r.clone();
        view!{<article><p>{format!("Reference {}: {} → {} · {:?} · eligible reads {} · deletions {}",s.position+1,s.reference_base,s.consensus_base,s.state,s.eligible_read_indices.len(),s.eligible_deletions)}</p>
        <button class="button secondary" on:click=move|_|{
            let result=action_guard.run(reference_context.clone()).and_then(|()|on_reference.run(crate::trace_validation::ReferenceNavigation{reference:ref_binding.clone(),selection:Some((position,position+1)),epoch:ref_epoch}));if let Err(e)=result{status.set(e);}
        }>"Open reference position"</button>
        <button class="button secondary" on:click=move|_|{
            if let Err(e)=action_guard.run(context.clone()){status.set(e);return}
            if let Some((ri,ci))=context.report.reads.iter().enumerate().find_map(|(ri,read)|read.columns.iter().position(|c|c.reference_position==Some(position)).map(|ci|(ri,ci))){read_index.set(ri);column_page.set(ci/50);status.set(format!("Call column {} for reference base {} is shown below; open its exact source call.",ci+1,position+1));}
            else{status.set("No mapped source call at this site; absence of coverage is not a measured base.".into());}
        }>"Show supporting call columns"</button></article>}
    }).collect_view())}
    <button disabled=move||site_page.get()==0 on:click=move|_|site_page.update(|p|*p=p.saturating_sub(1))>"Previous sites"</button><button disabled={move||current().is_none_or(|r|(site_page.get()+1)*25>=r.report.sites.iter().filter(|s|s.state!=SiteState::SupportedMatch).count())} on:click=move|_|site_page.update(|p|*p+=1)>"Next sites"</button>
    <h3>"Original read call mappings"</h3><select aria-label="Saved check read" prop:value=move||read_index.get().to_string() on:change=move|e|{read_index.set(event_target_value(&e).parse().unwrap_or(0));column_page.set(0);}>{move||current().map(|r|r.report.reads.iter().enumerate().map(|(i,s)|view!{<option value=i.to_string()>{format!("{} · {}",i+1,s.name.chars().take(160).collect::<String>())}</option>}).collect_view())}</select>
    {move||current().and_then(|r|r.report.reads.get(read_index.get()).map(|read|{let binding=r.reads.get(read_index.get()).cloned();let expected_epoch=r.epoch;let warnings=read.warnings.iter().take(32).map(|w|view!{<p>{w.clone()}</p>}).collect_view();let rows=read.columns.iter().skip(column_page.get()*50).take(50).map(|c|{let binding=binding.clone();let call=c.call_index;let context=r.clone();view!{<article class="read-check-call"><p>{format!("Reference base {:?} / boundary {} · {:?}: {:?}→{:?} · call {:?}, original {:?}, sample {:?}/original {:?} · raw {:?} · eligible {}",c.reference_position.map(|p|p+1),c.reference_boundary,c.kind,c.reference_base,c.observed_base,call.map(|v|v+1),c.stored_call_index.map(|v|v+1),c.sample_index,c.stored_sample_index,c.raw_confidence,c.eligible)}</p><button disabled=call.is_none()||binding.is_none() on:click=move|_|{if let Err(e)=action_guard.run(context.clone()){status.set(e);return}if let(Some(b),Some(at))=(&binding,call){on_read.run(read_navigation(b,at,expected_epoch));}}>"Open exact source call"</button></article>}}).collect_view();view!{<p>{format!("Placement unique: {} · search complete: {} · current read reversed: {} · original read reversed: {}",read.placement_unique_proven,read.search_complete,read.query_reversed,read.original_read_reverse_to_reference)}</p>{warnings}{read.error.clone().map(|e|view!{<p role="alert">{e}</p>})}{rows}}}))}
    <button disabled=move||column_page.get()==0 on:click=move|_|column_page.update(|p|*p=p.saturating_sub(1))>"Previous call columns"</button><button disabled={move||current().is_none_or(|r|r.report.reads.get(read_index.get()).is_none_or(|read|(column_page.get()+1)*50>=read.columns.len()))} on:click=move|_|column_page.update(|p|*p+=1)>"Next call columns"</button>
    </section>}
}
