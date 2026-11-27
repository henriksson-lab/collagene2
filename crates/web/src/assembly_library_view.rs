//! Reviewed library edits publish exact recovery and replacement together.
use crate::{
    assembly_draft_model::Draft,
    assembly_library_commit as bridge,
    assembly_library_topology::{self as topology, Change, Proposal},
    commit_lease::CommitLease,
    digest::SourceInfo,
    worker::CachedWorker,
};
use collagene_core::{
    Molecule,
    run_lineage::{self, EventBinding},
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

struct Review {
    proposal: Arc<Proposal>,
    epoch: u64,
    fingerprint: String,
    title: String,
}

#[component]
pub fn AssemblyLibraryView(
    draft: RwSignal<Draft>,
    epoch: RwSignal<u64>,
    sources: Memo<Vec<SourceInfo>>,
    get_source: Callback<EventBinding, Result<Molecule, String>>,
    snapshot: Callback<(), Result<Workspace, String>>,
    validate: Callback<Workspace, Result<(), String>>,
    on_restore: Callback<Workspace>,
    db: StoredValue<Option<crate::storage::Database>, LocalStorage>,
    workspace_busy: RwSignal<bool>,
    loaded: RwSignal<bool>,
    save_generation: RwSignal<u64>,
    message: RwSignal<String>,
) -> impl IntoView {
    let review = RwSignal::new(None::<Arc<Review>>);
    let acknowledged = RwSignal::new(false);
    let pending = RwSignal::new(false);
    let status = RwSignal::new(String::new());
    let edit_locked = RwSignal::new(false);
    let generation = StoredValue::new_local(Rc::new(Cell::new(0u64)));
    let lease = StoredValue::new_local(Rc::new(CommitLease::default()));
    let transaction =
        StoredValue::new_local(Rc::new(RefCell::new(None::<crate::storage::Transaction>)));
    let worker = StoredValue::new_local(None::<Rc<CachedWorker>>);
    let alive = Arc::new(AtomicBool::new(true));
    let cancel = Callback::new(move |()| {
        generation.with_value(|value| value.set(value.get().wrapping_add(1)));
        if let Some(worker) = worker.get_value() {
            worker.close();
        }
        worker.set_value(None);
        if let Some(tx) = transaction.with_value(|value| value.borrow().clone()) {
            let _ = tx.abort();
        } else if lease.with_value(|value| value.cancel_unstarted()) {
            workspace_busy.set(false);
        }
        pending.set(false);
        review.set(None);
        acknowledged.set(false);
    });
    let cleanup_alive = alive.clone();
    on_cleanup(move || {
        cleanup_alive.store(false, Ordering::Relaxed);
        generation.with_value(|value| value.set(value.get().wrapping_add(1)));
        if let Some(worker) = worker.get_value() {
            worker.close();
        }
        if let Some(tx) = transaction.with_value(|value| value.borrow().clone()) {
            let _ = tx.abort();
        } else if lease.with_value(|value| value.cancel_unstarted()) {
            workspace_busy.set(false);
        }
    });
    Effect::new(move |_| edit_locked.set(workspace_busy.get() || review.get().is_some()));
    Effect::new(move |_| {
        let now = epoch.get();
        let same = draft.with(|current| {
            review.with_untracked(|saved| {
                saved
                    .as_ref()
                    .is_none_or(|saved| saved.epoch == now && saved.proposal.before() == current)
            })
        });
        if !same {
            cancel.run(());
            status.set(
                "The draft or workspace changed. Review a new proposal before applying it.".into(),
            );
        }
    });
    let current = Callback::new(move |saved: Arc<Review>| -> Result<Workspace, String> {
        if epoch.get_untracked() != saved.epoch
            || !draft.with_untracked(|current| current == saved.proposal.before())
            || !review.with_untracked(|current| {
                current
                    .as_ref()
                    .is_some_and(|current| Arc::ptr_eq(current, &saved))
            })
        {
            return Err("The library edit proposal is stale; nothing was installed.".into());
        }
        let frame = snapshot.run(())?;
        if run_lineage::workspace_fingerprint(&frame).map_err(|e| e.to_string())?
            != saved.fingerprint
        {
            return Err("The workspace changed after this proposal. Review it again.".into());
        }
        Ok(frame)
    });
    let propose = Callback::new(
        move |(before, saved_epoch, change): (Arc<Draft>, u64, Change)| {
            let result = (|| -> Result<Arc<Review>, String> {
                if workspace_busy.get_untracked() || review.get_untracked().is_some() {
                    return Err("Apply or discard the existing proposal first.".into());
                }
                if epoch.get_untracked() != saved_epoch
                    || !draft.with_untracked(|current| current == before.as_ref())
                {
                    return Err(
                        "Source or draft context changed; configure this edit again.".into(),
                    );
                }
                let frame = snapshot.run(())?;
                let fingerprint =
                    run_lineage::workspace_fingerprint(&frame).map_err(|e| e.to_string())?;
                let title = match &change {
                    Change::Add { label, .. } => format!("Add source bin {label}"),
                    Change::Remove { .. } => "Remove source bin".into(),
                    Change::Reorder { .. } => "Reorder source bins".into(),
                };
                let proposal = Arc::new(topology::propose(&before, change)?);
                Ok(Arc::new(Review {
                    proposal,
                    epoch: saved_epoch,
                    fingerprint,
                    title,
                }))
            })();
            match result {
                Ok(saved) => {
                    review.set(Some(saved));
                    acknowledged.set(false);
                    status.set("Review every affected construct. Applying saves the exact current draft as recovery together with the replacement.".into());
                }
                Err(error) => status.set(error),
            }
        },
    );
    let restore_proposal = Callback::new(move |hash: String| {
        let result = (|| -> Result<Arc<Review>, String> {
            if workspace_busy.get_untracked() || review.get_untracked().is_some() {
                return Err("Apply or discard the existing proposal first.".into());
            }
            let before = draft.get_untracked();
            let frame = snapshot.run(())?;
            let fingerprint =
                run_lineage::workspace_fingerprint(&frame).map_err(|e| e.to_string())?;
            let proposal = Arc::new(bridge::proposal_restore(&frame, &before, &hash)?);
            let title = format!("Restore library draft: {}", proposal.after().name);
            Ok(Arc::new(Review {
                proposal,
                epoch: epoch.get_untracked(),
                fingerprint,
                title,
            }))
        })();
        match result {
            Ok(saved) => {
                review.set(Some(saved));
                acknowledged.set(false);
                status.set("Restoration is a new reviewed edit. Current raw work will be saved first, and identifier allocation cannot move backwards.".into());
            }
            Err(error) => status.set(error),
        }
    });
    let apply_alive = alive.clone();
    let apply = move |_| {
        if pending.get_untracked() || workspace_busy.get_untracked() {
            return;
        }
        let prepared = (|| -> Result<_, String> {
            if !acknowledged.get_untracked() {
                return Err("Review and acknowledge this proposal before applying it.".into());
            }
            let saved = review
                .get_untracked()
                .ok_or("Create an edit proposal first.")?;
            let frame = current.run(saved.clone())?;
            let staged = bridge::stage(
                &frame,
                &draft.get_untracked(),
                &saved.proposal,
                saved.proposal.revision(),
                &saved.fingerprint,
            )?;
            if staged.source_workspace_fingerprint != saved.fingerprint
                || staged.proposal_revision != saved.proposal.revision()
                || &staged.draft != saved.proposal.after()
                || staged.after_hash
                    != blake3::hash(&bridge::raw_encode(saved.proposal.after())?)
                        .to_hex()
                        .to_string()
            {
                return Err("Staged library edit differs from the reviewed proposal.".into());
            }
            let database = db
                .get_value()
                .ok_or("Browser recovery storage is unavailable; no edit was installed.")?;
            let encoder = Rc::new(CachedWorker::new()?);
            let nonce = generation.get_value();
            let token = nonce.get().wrapping_add(1);
            let lock = lease.get_value();
            lock.acquire(token)?;
            nonce.set(token);
            Ok((saved, frame, staged, database, encoder, nonce, token, lock))
        })();
        let (saved, previous, staged, database, encoder, nonce, token, lock) = match prepared {
            Ok(value) => value,
            Err(error) => {
                status.set(error);
                return;
            }
        };
        workspace_busy.set(true);
        save_generation.update(|value| *value = value.wrapping_add(1));
        worker.set_value(Some(encoder.clone()));
        pending.set(true);
        status.set("Saving the reviewed edit, exact recovery snapshot and previous workspace in one transaction…".into());
        let pending_transaction = transaction.get_value();
        let alive = apply_alive.clone();
        spawn_local(async move {
            let guard = || -> Result<(), String> {
                if !alive.load(Ordering::Relaxed) || nonce.get() != token {
                    return Err("Library edit cancelled.".into());
                }
                current.run(saved.clone()).map(|_| ())
            };
            let record_hash = staged.record_hash.clone();
            let result = crate::workspace_commit::persist(
                &database,
                staged.workspace,
                previous,
                validate,
                encoder.clone(),
                guard,
                pending_transaction,
            )
            .await;
            match result {
                Ok(frame) => {
                    loaded.set(true);
                    message.set(format!("Library edit and exact recovery saved together. Recovery {} is available in library history. Scientific evidence must be recomputed.",record_hash.get(..12).unwrap_or(&record_hash)));
                    // Durable success must reach App even if this panel closes
                    // immediately after the transaction's commit event.
                    on_restore.run(frame);
                }
                Err(error) => {
                    if alive.load(Ordering::Relaxed) && nonce.get() == token {
                        status.set(error);
                    }
                }
            }
            encoder.close();
            if lock.release(token) {
                workspace_busy.set(false);
            }
            if alive.load(Ordering::Relaxed) && nonce.get() == token {
                pending.set(false);
                worker.set_value(None);
            }
        });
    };
    view!{<details class="assembly-library-edit"><summary>"Edit source-bin layout and recover drafts"</summary>
        <p class="helper">"Add, remove or reorder bins across the selected construct library. Source identities and literal preparations are retained. Changed neighbors require primer/junction review; exact previous raw work is saved before publication."</p>
        <fieldset disabled=move||edit_locked.get() style="border:0;margin:0;padding:0;min-width:0">
            {move||{
                let saved=Arc::new(draft.get());let saved_epoch=epoch.get();
                saved.alternatives.as_ref().map(|library|{
                    let ids=library.bins.iter().map(|bin|bin.id.clone()).collect::<Vec<_>>();
                    let count=ids.len();
                    library.bins.iter().enumerate().map(|(index,bin)|{
                        let up=saved.clone();let down=saved.clone();let remove=saved.clone();
                        let up_ids=ids.clone();let down_ids=ids.clone();let bin_id=bin.id.clone();
                        view!{<article><p>{format!("{} · {} · {} source alternatives",index+1,bin.label,bin.choices.len())}</p>
                            <button class="text-tool" disabled=index==0 on:click=move|_|{let mut ids=up_ids.clone();if index==0||index>=ids.len(){return;}ids.swap(index,index-1);propose.run((up.clone(),saved_epoch,Change::Reorder{bin_ids:ids}));}>"Review moving up"</button>
                            <button class="text-tool" disabled={index+1>=count} on:click=move|_|{let mut ids=down_ids.clone();if index+1>=ids.len(){return;}ids.swap(index,index+1);propose.run((down.clone(),saved_epoch,Change::Reorder{bin_ids:ids}));}>"Review moving down"</button>
                            <button class="text-tool" disabled=count<=1 on:click=move|_|propose.run((remove.clone(),saved_epoch,Change::Remove{bin_id:bin_id.clone()}))>"Review removing bin"</button>
                        </article>}
                    }).collect_view()
                })
            }}
            <crate::assembly_library_add_view::AddLibraryBinView draft epoch busy=edit_locked sources get_source on_propose=propose/>
        </fieldset>
        {move||review.get().map(|saved|view!{<section><h3>{saved.title.clone()}</h3>
            <p>{format!("{} → {} bins · {} affected construct records",saved.proposal.before().slots.len(),saved.proposal.after().slots.len(),saved.proposal.impacts().len())}</p>
            <p>{format!("Before: {}",saved.proposal.before().slots.iter().map(|slot|slot.label.clone()).collect::<Vec<_>>().join(" → "))}</p>
            <p>{format!("After: {}",saved.proposal.after().slots.iter().map(|slot|slot.label.clone()).collect::<Vec<_>>().join(" → "))}</p>
            <details><summary>"Exact proposal identity"</summary><code>{saved.proposal.revision().to_owned()}</code><p>{saved.fingerprint.clone()}</p></details>
            <p class="helper">"Removed settings and exact raw formatting remain in the recovery snapshot. Restored candidate choices remain declarations only. Review every impact page before acknowledgement; changing the draft invalidates this proposal."</p>
            <LibraryImpacts proposal=saved.proposal.clone()/>
            <label><input type="checkbox" disabled=move||pending.get() prop:checked=move||acknowledged.get() on:change=move|event|acknowledged.set(event_target_checked(&event))/>"I reviewed this exact library change and its affected constructs."</label>
        </section>})}
        <button class="button primary" disabled=move||review.get().is_none()||!acknowledged.get()||pending.get()||workspace_busy.get() on:click=apply>"Apply reviewed change with recovery"</button>
        <button class="text-tool" disabled=move||review.get().is_none()&&!pending.get() on:click=move|_|{cancel.run(());status.set("Library change cancelled. No unpublished draft or staged recovery was installed.".into());}>"Cancel / discard proposal"</button>
        <p role="status">{move||status.get()}</p>
        <crate::assembly_library_history_view::LibraryHistoryView snapshot epoch busy=edit_locked on_restore_proposal=restore_proposal/>
    </details>}.into_any()
}

fn bin_label(proposal: &Proposal, id: &str) -> String {
    for draft in [proposal.after(), proposal.before()] {
        if let Some(library) = &draft.alternatives {
            if let Some(index) = library.bins.iter().position(|bin| bin.id == id) {
                let label = draft
                    .slots
                    .get(index)
                    .map(|slot| slot.label.as_str())
                    .unwrap_or(&library.bins[index].label);
                return format!("{label} ({id})");
            }
        }
    }
    id.to_owned()
}

#[component]
fn LibraryImpacts(proposal: Arc<Proposal>) -> impl IntoView {
    let page = RwSignal::new(0usize);
    let proposal = StoredValue::new(proposal);
    let count = proposal.with_value(|proposal| proposal.impacts().len());
    let pages = count.div_ceil(20).max(1);
    view!{<section>{move||proposal.with_value(|proposal|proposal.impacts().iter().skip(page.get().min(pages-1)*20).take(20).map(|row|view!{<article><h4>{format!("{} · {}{}{}",row.name,row.construct_id,if row.construct_removed{" · removed construct"}else{""},if row.construct_added{" · restored construct"}else{""})}</h4>
        <p>{format!("Primer neighbor review: {} · removed preparation bins: {}",row.primer_review_bins.iter().map(|id|bin_label(proposal,id)).collect::<Vec<_>>().join(", "),row.removed_preparation_bins.iter().map(|id|bin_label(proposal,id)).collect::<Vec<_>>().join(", "))}</p>
        <p>{format!("Internal changes retained {} · previous changes kept in recovery {}",row.retained_internal_edits,row.recovered_internal_edits)}</p>
<p>{format!("Retained literal junctions {} · previous junction declarations kept in recovery {} · new junctions (one-based) {:?}",row.retained_junctions,row.recovered_junctions,row.new_junctions.iter().map(|index|index+1).collect::<Vec<_>>())}</p>
        <p>{format!("Overlap text changed {} · source-role order changed {} · previous candidate declaration cleared {}",row.overlap_text_changed,row.role_order_changed,row.candidate_cleared)}</p>
    </article>}).collect_view())}<p>{move||format!("Impact page {} of {pages} · {count} construct records",page.get().min(pages-1)+1)}</p>
    <button class="text-tool" disabled=move||page.get()==0 on:click=move|_|page.update(|value|*value=value.saturating_sub(1))>"Previous impacts"</button>
    <button class="text-tool" disabled=move||page.get()+1>=pages on:click=move|_|page.update(|value|*value+=1)>"Next impacts"</button>
    </section>}.into_any()
}
