//! Explicit native-only local unlink; the backend owns vault/SQLite ordering.
use crate::library::remote_removal_model as model;
use crate::library::remote_removal_model::PinnedRemoval;
use crate::platform::desktop_sync_client::Client;
use crate::platform::desktop_sync_client::failure_message;
use collagene_desktop_sync_protocol as wire;
use leptos::prelude::*;
use std::{cell::Cell, rc::Rc};
use wasm_bindgen_futures::spawn_local;
#[component]
pub fn RemoteRemovalView(
    client: Client,
    selected: RwSignal<Option<wire::Remote>>,
    busy: RwSignal<bool>,
    locked: Memo<bool>,
    on_removed: Callback<wire::RemoteId>,
) -> impl IntoView {
    let client = StoredValue::new_local(client);
    let life = StoredValue::new_local(Rc::new(Cell::new(0u64)));
    let alive = StoredValue::new_local(Rc::new(Cell::new(true)));
    on_cleanup(move || alive.with_value(|flag| flag.set(false)));
    let pending = RwSignal::new(None::<PinnedRemoval>);
    let confirm = RwSignal::new(false);
    let status = RwSignal::new(String::new());
    let uncertain = RwSignal::new(None::<wire::Remote>);
    Effect::new(move |_| {
        let _ = selected.get();
        pending.set(None);
        confirm.set(false);
    });
    let prepare_alive = alive.clone();
    let prepare = move |_| {
        if busy.get_untracked() || locked.get_untracked() {
            return;
        }
        let Some(remote) = selected.get_untracked() else {
            return;
        };
        let client = client.get_value();
        let generation = life.get_value();
        let ticket = generation.get().wrapping_add(1);
        generation.set(ticket);
        busy.set(true);
        confirm.set(false);
        pending.set(None);
        let alive = prepare_alive.get_value();
        spawn_local(async move {
            let result = client
                .request(
                    wire::Command::PrepareRemoteRemoval {
                        remote: remote.id.clone(),
                    },
                    Some(remote.clone()),
                )
                .await;
            busy.try_set(false);
            if !alive.get()
                || generation.get() != ticket
                || selected.try_get_untracked().flatten().as_ref() != Some(&remote)
            {
                return;
            }
            match result {
                Ok(wire::Success::RemoteRemovalPrepared { review }) => {
                    match model::pin(&remote, review) {
                        Ok(p) => {
                            pending.set(Some(p));
                            status.set(
                                "Review the exact local unlink and key-forget scope below.".into(),
                            );
                        }
                        Err(e) => status.set(e),
                    }
                }
                Ok(_) => status.set("Unexpected removal review response.".into()),
                Err(e) => status.set(e),
            }
        });
    };
    let confirm_alive = alive.clone();
    let remove = move |_| {
        if busy.get_untracked() || locked.get_untracked() || !confirm.get_untracked() {
            return;
        }
        let Some(p) = pending.get_untracked() else {
            return;
        };
        if !model::is_current(&p, selected.get_untracked().as_ref()) {
            status.set("Database selection changed; prepare a fresh removal review.".into());
            return;
        }
        let client = client.get_value();
        let generation = life.get_value();
        let ticket = generation.get().wrapping_add(1);
        generation.set(ticket);
        busy.set(true);
        confirm.set(false);
        uncertain.set(Some(p.remote.clone()));
        let alive = confirm_alive.get_value();
        spawn_local(async move {
            let result = client
                .request_outcome(
                    wire::Command::ConfirmRemoteRemoval {
                        remote: p.remote.id.clone(),
                        review: p.review.id.clone(),
                    },
                    Some(p.remote.clone()),
                )
                .await;
            busy.try_set(false);
            if !alive.get() || generation.get() != ticket {
                return;
            }
            match result{Ok(wire::Outcome::Success(wire::Success::RemoteRemoved(receipt)))=>match model::validate_receipt(&p.remote,Some(&p),&receipt){Ok(())=>{pending.set(None);uncertain.set(None);status.set("Database unlinked locally and its key forgotten. All local designs remain.".into());on_removed.run(p.remote.id.clone());},Err(e)=>status.set(e)},Ok(wire::Outcome::Failure(failure))=>{status.set(failure_message(&failure));if matches!(failure,wire::Failure::Error{code:wire::ErrorCode::ReviewChanged|wire::ErrorCode::CredentialForgetFailed,..}){pending.set(None);}},Ok(_)=>status.set("Unlink outcome could not be confirmed. Check the durable outcome before retrying.".into()),Err(e)=>status.set(format!("{e} Check the durable unlink outcome before preparing another action."))}
        });
    };
    let check_alive = alive.clone();
    let check = move |_| {
        if busy.get_untracked() || locked.get_untracked() {
            return;
        }
        let Some(remote) = uncertain
            .get_untracked()
            .or_else(|| selected.get_untracked())
        else {
            return;
        };
        let p = pending.get_untracked().filter(|p| p.remote == remote);
        let client = client.get_value();
        let generation = life.get_value();
        let ticket = generation.get().wrapping_add(1);
        generation.set(ticket);
        busy.set(true);
        let alive = check_alive.get_value();
        spawn_local(async move {
            let result = client
                .request(
                    wire::Command::CheckRemoteRemovalOutcome {
                        remote: remote.id.clone(),
                    },
                    Some(remote.clone()),
                )
                .await;
            busy.try_set(false);
            if !alive.get() || generation.get() != ticket {
                return;
            }
            match result {
                Ok(wire::Success::RemoteRemovalOutcome {
                    remote: actual,
                    receipt: Some(receipt),
                }) if actual == remote.id => {
                    match model::validate_receipt(&remote, p.as_ref(), &receipt) {
                        Ok(()) => {
                            pending.set(None);
                            uncertain.set(None);
                            status.set("Durable outcome confirms this database was unlinked locally; local designs remain.".into());
                            on_removed.run(remote.id.clone());
                        }
                        Err(e) => status.set(e),
                    }
                }
                Ok(wire::Success::RemoteRemovalOutcome {
                    remote: actual,
                    receipt: None,
                }) if actual == remote.id => {
                    uncertain.set(None);
                    status.set("No completed local unlink receipt was found. Configuration may remain after its key was forgotten; refresh and prepare a new review.".into());
                }
                Ok(_) => status.set("Unexpected removal outcome response.".into()),
                Err(e) => status.set(e),
            }
        });
    };
    let discard_alive = alive;
    let dismiss = move |_| {
        if busy.get_untracked() {
            return;
        }
        let p = pending.get_untracked();
        pending.set(None);
        confirm.set(false);
        if let Some(p) = p {
            let client = client.get_value();
            let alive = discard_alive.get_value();
            spawn_local(async move {
                let result = client
                    .request(
                        wire::Command::DiscardRemoteRemovalReview {
                            review: p.review.id,
                        },
                        None,
                    )
                    .await;
                if alive.get() {
                    if let Err(e) = result {
                        status.set(e);
                    }
                }
            });
        }
    };
    view! {<details><summary>"Remove database from this device"</summary>
       <p class="helper">"This forgets its exact scoped access key and unlinks local remote bookkeeping. Local plasmids, folders and full history remain; nothing is deleted from the remote server. Pending uploads and unfinished folder clones must be resolved first."</p>
       <button class="text-tool" disabled=move||busy.get()||locked.get()||selected.get().is_none() on:click=prepare>"Prepare removal review"</button>
       {move||pending.get().map(|p|{let d=p.review.description;view!{<section class="library-revisions"><strong>{d.remote.name.into_string()}</strong><p class="library-status">{format!("{} · remote lab {} · local repository {}",d.remote.endpoint.as_str(),d.remote.lab.as_str(),d.local_lab.as_str())}</p><p>{format!("Unlink {} plasmid bindings, {} folder bindings and {} completed folder-clone records. Preserve every local design and archive.",d.document_bindings,d.folder_bindings,d.completed_tree_jobs)}</p><p class="helper">"Key-forget and local unlink are separate operations. If the native vault refuses, configuration stays. If a concurrent change blocks unlink after key-forget, configuration may stay without its key; the result reports this explicitly."</p><label class="checkbox-label"><input type="checkbox" prop:checked=move||confirm.get() disabled=move||busy.get()||locked.get() on:change=move|e|confirm.set(event_target_checked(&e))/>"I reviewed this exact database and want to forget its key and unlink it locally."</label><button class="button secondary" disabled=move||busy.get()||locked.get()||!confirm.get() on:click=remove.clone()>"Forget key and unlink database"</button><button class="text-tool" disabled=move||busy.get() on:click=dismiss.clone()>"Dismiss review"</button></section>}})}
       <p class="helper">{move||uncertain.get().map(|r|format!("Last unlink attempt: {} · {}",r.name.as_str(),r.endpoint.as_str())).unwrap_or_default()}</p>
    <button class="text-tool" disabled=move||busy.get()||locked.get()||(selected.get().is_none()&&uncertain.get().is_none()) on:click=check>"Check last unlink outcome"</button><p class="library-status" role="status">{move||status.get()}</p></details>}
}
