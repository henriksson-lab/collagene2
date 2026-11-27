//! Native-only remote database accordion. Clone/pull save into the local catalog;
//! App refreshes that catalog but deliberately does not replace its working editor.
use crate::desktop_sync_client::Client;
use collagene_desktop_sync_protocol::{self as wire, library as base};
use leptos::prelude::*;
use std::{cell::Cell, rc::Rc};
use wasm_bindgen_futures::spawn_local;

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct SavedLocalDocument {
    pub id: base::DocumentId,
    pub revision: base::RevisionId,
    pub name: String,
    /// Working editor differs from this saved local record. Save locally first.
    pub dirty: bool,
}
#[derive(Clone, Default)]
struct Navigation {
    path: Vec<base::Folder>,
    query: String,
    offset: u32,
    generation: Option<u64>,
}
impl Navigation {
    fn parent(&self) -> Option<base::FolderId> {
        self.path.last().map(|f| f.id.clone())
    }
    fn command(&self, remote: wire::RemoteId) -> Result<wire::Command, String> {
        let page = base::PageRequest {
            offset: self.offset,
            limit: 100,
            expected_catalog_generation: self.generation,
        };
        if self.query.trim().is_empty() {
            Ok(wire::Command::Browse {
                remote,
                parent: self.parent(),
                page,
            })
        } else {
            Ok(wire::Command::Search {
                remote,
                query: base::Text::new(self.query.trim()).map_err(|_| "Search is too long.")?,
                scope: base::SearchScope::AllFolders,
                page,
            })
        }
    }
}
#[derive(Clone)]
struct History {
    document: base::Document,
    offset: u32,
    generation: Option<u64>,
    next: Option<u32>,
    rows: Vec<base::RevisionSummary>,
}
#[derive(Clone)]
struct PullReview {
    local: SavedLocalDocument,
    remote: wire::Remote,
    failure: wire::Failure,
}
fn review_is_current(
    review: &PullReview,
    local: Option<&SavedLocalDocument>,
    remote: Option<&wire::Remote>,
) -> bool {
    local.is_some_and(|local| {
        !local.dirty && local.id == review.local.id && local.revision == review.local.revision
    }) && remote == Some(&review.remote)
}
fn history_conflict(failure: &wire::Failure) -> bool {
    matches!(
        failure,
        wire::Failure::Error {
            code: wire::ErrorCode::HistoryConflict,
            ..
        }
    )
}
struct Action {
    command: wire::Command,
    navigation: Option<Navigation>,
}
impl Action {
    fn command(command: wire::Command) -> Self {
        Self {
            command,
            navigation: None,
        }
    }
}
fn metadata_changed(success: &wire::Success) -> bool {
    matches!(
        success,
        wire::Success::Cloned { .. } | wire::Success::Pulled { .. } | wire::Success::TreeUpdated(_)
    )
}
fn note(success: &wire::Success) -> String {
    match success {
        wire::Success::Cloned{..}|wire::Success::Pulled{..}=>"Saved to the local library. Reopen it there to review; the working editor has not changed.".into(),
        wire::Success::Queued(_)=>"Push queued durably. Resume it from Pending pushes.".into(),
        wire::Success::Resumed(r)=>match &r.outcome {
            wire::PushOutcome::Applied=>"Remote accepted the saved mutation.".into(),
            wire::PushOutcome::Rejected(f)=>format!("Remote rejected the saved mutation ({:?}): {}",f.code,f.message.as_str()),
        },
        wire::Success::CredentialSet{retention,..}=>match retention {wire::CredentialRetention::Session=>"Access key kept in native session memory.".into(),wire::CredentialRetention::Remember=>"Access key stored in the native credential vault.".into()},
        wire::Success::CredentialDeleted{..}=>"Remote credential forgotten.".into(),
        wire::Success::CredentialStatus{state,..}=>format!("Credential: {state:?}"),
        wire::Success::SessionsCleared{..}=>"Native session credential cleared; remembered vault entries are retained.".into(),
        wire::Success::CancelledUnattempted{..}=>"Unattempted push cancelled.".into(),
        wire::Success::DiscardedRejected{..}=>"Rejected push discarded. Local data is retained.".into(),
        wire::Success::TreeStarted(_)=>"Folder clone saved. Resume in Folder clones; each step is bounded.".into(),
        wire::Success::TreeUpdated(t)=>format!("Folder clone {:?}: {} folders, {} plasmids saved locally.{}",t.phase,t.folders_created,t.documents_created,if t.issue.is_some(){" Review its issue before continuing."}else{""}),
        _=>String::new(),
    }
}

#[component]
pub fn RemoteLibraryPane(
    client: Client,
    /// Host-derived native context. Ordinary web must never mount this pane.
    native: bool,
    saved: Memo<Option<SavedLocalDocument>>,
    /// Selected local destination. None means local catalog root.
    destination: Memo<Option<base::FolderId>>,
    locked: Memo<bool>,
    on_local_change: Callback<()>,
) -> impl IntoView {
    let client = StoredValue::new_local(client);
    let life = StoredValue::new_local(Rc::new(Cell::new(0u64)));
    let busy = RwSignal::new(false);
    let status = RwSignal::new(String::new());
    let remotes = RwSignal::new(Vec::<wire::Remote>::new());
    let selected = RwSignal::new(None::<wire::Remote>);
    let name = RwSignal::new(String::new());
    let endpoint = RwSignal::new(String::new());
    let expected_lab = RwSignal::new(String::new());
    let key = NodeRef::<leptos::html::Input>::new();
    let remember = RwSignal::new(false);
    let navigation = RwSignal::new(Navigation::default());
    let query = RwSignal::new(String::new());
    let entries = RwSignal::new(Vec::<base::Entry>::new());
    let next = RwSignal::new(None::<u32>);
    let history = RwSignal::new(None::<History>);
    let pending = RwSignal::new(Vec::<wire::PendingPush>::new());
    let pending_offset = RwSignal::new(0u16);
    let pending_next = RwSignal::new(None::<u16>);
    let trees = RwSignal::new(Vec::<wire::TreeStatus>::new());
    let tree_offset = RwSignal::new(0u16);
    let tree_next = RwSignal::new(None::<u16>);
    let checkout = RwSignal::new(false);
    let review = RwSignal::new(None::<PullReview>);
    let bookmarks = RwSignal::new(None::<wire::BookmarkPage>);
    let bookmark_offset = RwSignal::new(0u16);
    // A review never follows a different editor selection, saved head, dirty
    // state or remote. Its acceptance always names the original exact pins.
    Effect::new(move |_| {
        let _ = (saved.get(), selected.get());
        checkout.set(false);
        review.set(None);
    });
    Effect::new(move |_| {
        let _ = selected.get();
        history.set(None);
    });
    on_cleanup(move || {
        life.with_value(|v| v.set(v.get().wrapping_add(1)));
        if let Some(input) = key.get() {
            input.set_value("");
        }
    });
    let send = Callback::new(move |action: Action| {
        if !native || busy.get_untracked() || locked.get_untracked() {
            return;
        }
        if matches!(
            &action.command,
            wire::Command::SetCredential { .. }
                | wire::Command::DeleteCredential { .. }
                | wire::Command::ClearSession { .. }
                | wire::Command::ConfirmRemoteRemoval { .. }
        ) {
            history.set(None);
        }
        let client = client.get_value();
        let life = life.get_value();
        let epoch = life.get().wrapping_add(1);
        life.set(epoch);
        let remote_id = action
            .command
            .remote_id()
            .cloned()
            .or_else(|| match &action.command {
                wire::Command::ResumePush { job }
                | wire::Command::CancelUnattempted { job }
                | wire::Command::DiscardRejected { job } => pending.with(|rows| {
                    rows.iter()
                        .find(|p| &p.job == job)
                        .map(|p| p.remote.clone())
                }),
                wire::Command::ResumeTreeClone { job, .. }
                | wire::Command::RetryTreeClone { job }
                | wire::Command::QueueTreeFolderMetadata { job, .. }
                | wire::Command::QueueTreeDocumentMetadata { job, .. } => trees.with(|rows| {
                    rows.iter()
                        .find(|t| &t.job == job)
                        .map(|t| t.remote.clone())
                }),
                _ => None,
            });
        let expected =
            remote_id.and_then(|id| remotes.with(|rows| rows.iter().find(|r| r.id == id).cloned()));
        let review_context = match &action.command {
            wire::Command::Pull {
                local_document,
                expected_local_revision,
                ..
            }
            | wire::Command::QueuePush {
                local_document,
                expected_local_revision,
                ..
            } => saved
                .get_untracked()
                .filter(|local| {
                    !local.dirty
                        && &local.id == local_document
                        && &local.revision == expected_local_revision
                })
                .zip(expected.clone()),
            _ => None,
        };
        if let Some(nav) = &action.navigation {
            navigation.set(nav.clone());
        }
        // Credentials are handed directly to native execution, never retained in a
        // signal, status, debug formatter, retry closure or local storage.
        busy.set(true);
        status.set("Working…".into());
        spawn_local(async move {
            let result = client.request_outcome(action.command, expected).await;
            // A real native mutation may finish after navigation changed. Refresh
            // the local catalog even then, without loading anything into the editor.
            if matches!(&result,Ok(wire::Outcome::Success(success)) if metadata_changed(success)) {
                on_local_change.run(());
            }
            if life.get() != epoch {
                return;
            }
            busy.set(false);
            match result {
                Err(error) => status.set(error),
                Ok(wire::Outcome::Failure(failure)) => {
                    status.set(crate::desktop_sync_client::failure_message(&failure));
                    if matches!(
                        &failure,
                        wire::Failure::Error {
                            code: wire::ErrorCode::BookmarkChanged,
                            ..
                        }
                    ) {
                        bookmarks.set(None);
                        bookmark_offset.set(0);
                    }
                    if history_conflict(&failure)
                        || matches!(&failure, wire::Failure::LocalConflict { .. })
                    {
                        if let Some((local, remote)) = review_context {
                            checkout.set(false);
                            let pinned = PullReview {
                                local,
                                remote,
                                failure,
                            };
                            if review_is_current(
                                &pinned,
                                saved.get_untracked().as_ref(),
                                selected.get_untracked().as_ref(),
                            ) {
                                review.set(Some(pinned));
                            }
                        }
                    }
                }
                Ok(wire::Outcome::Success(success)) => {
                    status.set(note(&success));
                    match success {
                        wire::Success::Remotes(rows) => remotes.set(rows.into_vec()),
                        wire::Success::Bookmarks(page) => bookmarks.set(Some(page)),
                        wire::Success::AddedRemote(remote) => {
                            remotes.update(|rs| rs.push(remote.clone()));
                            selected.set(Some(remote));
                            navigation.set(Navigation::default());
                            query.set(String::new());
                            entries.set(vec![]);
                            next.set(None);
                            history.set(None);
                            checkout.set(false);
                            name.set(String::new());
                            endpoint.set(String::new());
                            expected_lab.set(String::new());
                        }
                        wire::Success::Browsed(page) | wire::Success::Searched(page) => {
                            navigation.update(|n| n.generation = Some(page.catalog_generation));
                            entries.set(page.entries.into_vec());
                            next.set(page.next_offset);
                        }
                        wire::Success::Revisions(page) => history.update(|current| {
                            if let Some(h) = current {
                                h.generation = Some(page.catalog_generation);
                                h.next = page.next_offset;
                                h.rows = page.revisions.into_vec();
                            }
                        }),
                        wire::Success::Queued(push) => pending.update(|rows| {
                            rows.retain(|p| p.job != push.job);
                            rows.push(push);
                        }),
                        wire::Success::Pending(page) => {
                            pending.set(page.items.into_vec());
                            pending_next.set(page.next_offset);
                        }
                        wire::Success::Resumed(receipt) => {
                            pending.update(|rows| match receipt.outcome {
                                wire::PushOutcome::Applied => rows.retain(|p| p.job != receipt.job),
                                wire::PushOutcome::Rejected(rejection) => {
                                    if let Some(p) = rows.iter_mut().find(|p| p.job == receipt.job)
                                    {
                                        p.attempted = true;
                                        p.rejection = Some(rejection);
                                    }
                                }
                            })
                        }
                        wire::Success::CancelledUnattempted { job }
                        | wire::Success::DiscardedRejected { job } => {
                            pending.update(|rows| rows.retain(|p| p.job != job))
                        }
                        wire::Success::TreeStarted(tree) | wire::Success::TreeUpdated(tree) => {
                            trees.update(|rows| {
                                rows.retain(|p| p.job != tree.job);
                                rows.push(tree);
                            })
                        }
                        wire::Success::TreeJobs(page) => {
                            trees.set(page.items.into_vec());
                            tree_next.set(page.next_offset);
                        }
                        _ => {}
                    }
                }
            }
        });
    });
    Effect::new(move |_| {
        if native {
            send.run(Action::command(wire::Command::Remotes));
        }
    });
    let browse = Callback::new(move |nav: Navigation| {
        let Some(remote) = selected.get_untracked() else {
            return;
        };
        match nav.command(remote.id) {
            Ok(command) => {
                history.set(None);
                send.run(Action {
                    command,
                    navigation: Some(nav),
                });
            }
            Err(error) => status.set(error),
        }
    });
    let read_history = Callback::new(move |h: History| {
        let Some(remote) = selected.get_untracked() else {
            return;
        };
        let command = wire::Command::ListRevisions {
            remote: remote.id,
            document: h.document.id.clone(),
            page: base::PageRequest {
                offset: h.offset,
                limit: 100,
                expected_catalog_generation: h.generation,
            },
        };
        history.set(Some(h));
        send.run(Action::command(command));
    });
    let choose_remote = move |event: leptos::ev::Event| {
        life.with_value(|life| life.set(life.get().wrapping_add(1)));
        busy.set(false);
        status.set(String::new());
        let id = event_target_value(&event);
        selected.set(remotes.with(|rows| rows.iter().find(|r| r.id.as_str() == id).cloned()));
        query.set(String::new());
        navigation.set(Navigation::default());
        entries.set(vec![]);
        next.set(None);
        history.set(None);
        checkout.set(false);
        if let Some(input) = key.get() {
            input.set_value("");
        }
    };
    let submit_key = move |event: leptos::ev::SubmitEvent| {
        event.prevent_default();
        // Clear the actual DOM value before validation or asynchronous dispatch.
        let value = key
            .get()
            .map(|input| {
                let value = input.value();
                input.set_value("");
                value
            })
            .unwrap_or_default();
        let Some(remote) = selected.get_untracked() else {
            status.set("Select a database first.".into());
            return;
        };
        match wire::Secret::new(&value) {
            Ok(secret) => send.run(Action::command(wire::Command::SetCredential {
                remote: remote.id,
                secret,
                retention: if remember.get_untracked() {
                    wire::CredentialRetention::Remember
                } else {
                    wire::CredentialRetention::Session
                },
            })),
            Err(_) => status
                .set("Enter a valid access key (at most 2048 nonspace ASCII characters).".into()),
        }
    };
    let add_remote = move |event: leptos::ev::SubmitEvent| {
        event.prevent_default();
        let command = (|| {
            Ok::<_, String>(wire::Command::AddRemote {
                name: base::Text::new(name.get_untracked())
                    .map_err(|_| "Database name is too long.")?,
                endpoint: base::Text::new(endpoint.get_untracked())
                    .map_err(|_| "Origin is too long.")?,
                expected_lab: base::LabId::new(expected_lab.get_untracked().trim())
                    .map_err(|_| "Enter the expected lab ID supplied by its administrator.")?,
            })
        })();
        match command {
            Ok(command) => send.run(Action::command(command)),
            Err(error) => status.set(error),
        }
    };
    let push = move |_| {
        let Some(local) = saved.get_untracked() else {
            status.set("Save this plasmid in the local library first.".into());
            return;
        };
        if local.dirty {
            status.set("Save current changes to the local library before pushing.".into());
            return;
        }
        let Some(remote) = selected.get_untracked() else {
            return;
        };
        send.run(Action::command(wire::Command::QueuePush {
            local_document: local.id,
            expected_local_revision: local.revision,
            remote: remote.id,
        }));
    };
    let pull = move |_| {
        let Some(local) = saved.get_untracked() else {
            return;
        };
        if local.dirty {
            status.set("Save current changes to the local library before pulling.".into());
            return;
        }
        let Some(remote) = selected.get_untracked() else {
            return;
        };
        review.set(None);
        send.run(Action::command(wire::Command::Pull {
            local_document: local.id,
            expected_local_revision: local.revision,
            remote: remote.id,
            mode: wire::PullMode::CleanFastForward,
        }));
    };
    let accept_checkout = Callback::new(move |()| {
        let Some(pinned) = review.get_untracked() else {
            return;
        };
        if !checkout.get_untracked()
            || !history_conflict(&pinned.failure)
            || !review_is_current(
                &pinned,
                saved.get_untracked().as_ref(),
                selected.get_untracked().as_ref(),
            )
        {
            status.set(
                "The reviewed selection changed. Review a fresh clean pull before checkout.".into(),
            );
            checkout.set(false);
            return;
        }
        checkout.set(false);
        review.set(None);
        send.run(Action::command(wire::Command::Pull {
            local_document: pinned.local.id,
            expected_local_revision: pinned.local.revision,
            remote: pinned.remote.id,
            mode: wire::PullMode::ExplicitCheckout,
        }));
    });
    view! {
        <Show when=move||native>
        <details class="lab-library-pane remote-library-pane"><summary>"Databases"</summary>
            <label>"Configured database"<select aria-label="Configured database" on:change=choose_remote prop:value=move||selected.get().map(|r|r.id.as_str().to_owned()).unwrap_or_default()>
                <option value="">"Select a database"</option>
                {move||remotes.get().into_iter().map(|r|view!{<option value=r.id.as_str().to_owned()>{r.name.into_string()}</option>}).collect_view()}
            </select></label>
            <details><summary>"Add database"</summary><form class="library-connect" on:submit=add_remote>
                <label>"Name"<input maxlength="256" prop:value=move||name.get() on:input=move|e|name.set(event_target_value(&e))/></label>
                <label>"Server origin"<input type="url" placeholder="https://lab.example" maxlength="2048" prop:value=move||endpoint.get() on:input=move|e|endpoint.set(event_target_value(&e))/></label>
                <label>"Expected lab ID"<input placeholder="lab-…" maxlength="36" prop:value=move||expected_lab.get() on:input=move|e|expected_lab.set(event_target_value(&e))/></label>
                <button class="button secondary" type="submit" disabled=move||busy.get()||locked.get()>"Add database"</button>
            </form></details>
            <crate::remote_removal_view::RemoteRemovalView
                client=client.get_value() selected busy locked
                on_removed=Callback::new(move |id: wire::RemoteId| {
                    remotes.update(|items|items.retain(|remote|remote.id!=id));
                    if selected.with_untracked(|remote|remote.as_ref().is_some_and(|remote|remote.id==id)) {
                        selected.set(None); navigation.set(Navigation::default()); entries.set(vec![]); history.set(None);
                    }
                    status.set("Database unlinked from this device. All local designs and history remain.".into());
                })
            />
            <Show when=move||selected.with(Option::is_some)>
                <p class="helper">{move||selected.get().map(|r|format!("{} · {}",r.endpoint.as_str(),r.expected_lab.as_str())).unwrap_or_default()}</p>
                <details><summary>"Access key"</summary>
                    <form class="library-connect" on:submit=submit_key>
                        <label>"Access key"<input node_ref=key type="password" autocomplete="off" maxlength="2048" aria-label="Database access key" disabled=move||busy.get()||locked.get()/></label>
                        <label><input type="checkbox" prop:checked=move||remember.get() on:change=move|e|remember.set(event_target_checked(&e))/>"Remember in native credential vault"</label>
                        <p class="helper">"Unchecked: native session memory only. The field clears on submit."</p>
                        <button class="button secondary" type="submit" disabled=move||busy.get()||locked.get()>"Store key"</button>
                    </form>
                    <button class="text-tool" disabled=move||busy.get()||locked.get() on:click=move|_|if let Some(r)=selected.get_untracked(){send.run(Action::command(wire::Command::CredentialStatus{remote:r.id}));}>"Credential status"</button>
                    <button class="text-tool" disabled=move||busy.get()||locked.get() on:click=move|_|if let Some(r)=selected.get_untracked(){send.run(Action::command(wire::Command::ClearSession{remote:Some(r.id)}));}>"Clear session key"</button>
                    <button class="text-tool" disabled=move||busy.get()||locked.get() on:click=move|_|if let Some(r)=selected.get_untracked(){send.run(Action::command(wire::Command::DeleteCredential{remote:r.id}));}>"Forget key"</button>
                </details>
                <details open><summary>"Remote plasmids"</summary>
                    <div class="library-location"><button class="text-tool" disabled=move||busy.get()||locked.get() on:click=move|_|{query.set(String::new());browse.run(Navigation::default());}>"Root"</button>
                        {move||navigation.with(|n|n.path.iter().enumerate().map(|(i,f)|{let path=n.path[..=i].to_vec();let label=f.name.to_string();view!{<span>" / "</span><button class="text-tool" disabled=move||busy.get()||locked.get() on:click=move|_|{query.set(String::new());browse.run(Navigation{path:path.clone(),..Default::default()});}>{label}</button>}}).collect_view())}
                    </div>
                    <form class="library-search" on:submit=move|e|{e.prevent_default();browse.run(Navigation{query:query.get_untracked(),..Default::default()});}>
                        <input type="search" aria-label="Search remote plasmids" maxlength="256" placeholder="Search database" prop:value=move||query.get() on:input=move|e|query.set(event_target_value(&e))/>
                        <button type="submit" class="text-tool" disabled=move||busy.get()||locked.get()>"Find"</button>
                    </form>
                    <ul class="lab-library-entries">{move||entries.get().into_iter().map(|entry|match entry{
                        base::Entry::Folder(folder)=>{let open=folder.clone();let tree=folder.clone();view!{<li><button class="text-tool" disabled=move||busy.get()||locked.get() on:click=move|_|{let mut path=navigation.get_untracked().path;path.push(open.clone());query.set(String::new());browse.run(Navigation{path,..Default::default()});}>{folder.name.into_string()}</button><details><summary>"Folder actions"</summary><button class="text-tool" disabled=move||busy.get()||locked.get() on:click=move|_|if let Some(r)=selected.get_untracked(){send.run(Action::command(wire::Command::BeginTreeClone{remote:r.id,root:wire::TreeRoot::Folder{id:tree.id.clone(),parent:tree.parent.clone()},destination:destination.get_untracked(),limits:wire::TreeLimits::default()}));}>"Clone folder to local library"</button></details></li>}.into_any()},
                        base::Entry::Document(document)=>{let clone=document.clone();let history_doc=document.clone();view!{<li><span>{document.name.into_string()}</span><details><summary>"Plasmid actions"</summary><button class="text-tool" disabled=move||busy.get()||locked.get() on:click=move|_|if let Some(r)=selected.get_untracked(){send.run(Action::command(wire::Command::CloneDocument{remote:r.id,document:clone.id.clone(),folder:destination.get_untracked()}));}>"Clone to local library"</button><button class="text-tool" disabled=move||busy.get()||locked.get() on:click=move|_|read_history.run(History{document:history_doc.clone(),offset:0,generation:None,next:None,rows:vec![]})>"Revision history"</button></details></li>}.into_any()},
                    }).collect_view()}</ul>
                    <div class="library-pages"><button class="text-tool" disabled=move||busy.get()||locked.get()||navigation.with(|n|n.offset==0) on:click=move|_|{let mut n=navigation.get_untracked();n.offset=n.offset.saturating_sub(100);browse.run(n);}>"Previous"</button><button class="text-tool" disabled=move||busy.get()||locked.get()||next.get().is_none() on:click=move|_|if let Some(offset)=next.get_untracked(){let mut n=navigation.get_untracked();n.offset=offset;browse.run(n);}>"Next"</button></div>
                    <details><summary>"Clone database folders"</summary><p class="helper">"Creates a separate local tree at the selected local destination. Saved local heads are retained."</p><button class="text-tool" disabled=move||busy.get()||locked.get() on:click=move|_|if let Some(r)=selected.get_untracked(){send.run(Action::command(wire::Command::BeginTreeClone{remote:r.id,root:wire::TreeRoot::LabRoot,destination:destination.get_untracked(),limits:wire::TreeLimits::default()}));}>"Begin root tree clone"</button></details>
                </details>
                <Show when=move||history.with(Option::is_some)><details open><summary>"Remote revision history"</summary>
                    <p>{move||history.get().map(|h|h.document.name.into_string()).unwrap_or_default()}</p>
                    <ul>{move||history.get().and_then(|h|selected.get().map(|remote|(h,remote))).map(|(h,remote)|h.rows.into_iter().map(|r|{
                        let target=base::discussion::Target{document:h.document.id.clone(),revision:r.id.clone()};
                        let connection=crate::library_discussion_client::Client::Remote{client:client.get_value(),remote:remote.clone()};
                        view!{<li>{format!("Revision {} · {} · {} bytes",r.number,r.id.as_str(),r.archive_size)}
                            <crate::library_discussion_view::LibraryDiscussion target=target client=connection/>
                        </li>}
                    }).collect_view())}</ul>
                    <p class="helper">"Sequence history is immutable. Comments and human review apply to each exact remote revision; cloning fetches its current head."</p>
                    <button class="text-tool" disabled=move||busy.get()||locked.get()||history.with(|h|h.as_ref().is_none_or(|h|h.offset==0)) on:click=move|_|if let Some(mut h)=history.get_untracked(){h.offset=h.offset.saturating_sub(100);read_history.run(h);}>"Previous"</button>
                    <button class="text-tool" disabled=move||busy.get()||locked.get()||history.with(|h|h.as_ref().is_none_or(|h|h.next.is_none())) on:click=move|_|if let Some(mut h)=history.get_untracked(){if let Some(offset)=h.next{h.offset=offset;read_history.run(h);}}>"Next"</button>
                </details></Show>
                <details><summary>"Saved local plasmid"</summary>
                    <p>{move||saved.get().map(|s|s.name).unwrap_or_else(||"Save a plasmid to the local library first.".into())}</p>
                    <Show when=move||saved.with(|s|s.as_ref().is_some_and(|s|s.dirty))><p class="helper">"Unsaved editor changes: save to the local library before push or pull."</p></Show>
                    <button class="text-tool" disabled=move||busy.get()||locked.get()||saved.with(|s|s.as_ref().is_none_or(|s|s.dirty)) on:click=push>"Queue saved revision push"</button>
                    <button class="text-tool" disabled=move||busy.get()||locked.get()||saved.with(|s|s.as_ref().is_none_or(|s|s.dirty)) on:click=pull>"Pull clean fast-forward"</button>
                    {move||review.get().map(|pinned|{
                        let can_checkout=history_conflict(&pinned.failure);
                        let conflict=match &pinned.failure {
                            wire::Failure::LocalConflict{current}=>format!("The saved local head is now {}@{}. Refresh the catalog and review before another action.",current.id.as_str(),current.head.id.as_str()),
                            _=>"The host could not reconcile these saved native histories as a clean fast-forward. No checkout was performed.".into(),
                        };
                        view!{<section class="library-revisions"><div class="library-heading"><strong>"Pull history review"</strong><button class="text-tool" on:click=move|_|{review.set(None);checkout.set(false);}>"Dismiss"</button></div>
                            <p class="library-status">{conflict}</p>
                            <p class="library-status">{format!("Reviewed local record: {}@{} · database {}",pinned.local.id.as_str(),pinned.local.revision.as_str(),pinned.remote.name.as_str())}</p>
                            <p class="library-status">"Open this plasmid’s Revision history in the Local library to inspect the complete pre-pull archive. The working editor is retained."</p>
                            <button class="text-tool" disabled=move||busy.get()||locked.get() on:click=move|_|on_local_change.run(())>"Refresh local catalog metadata"</button>
                            {can_checkout.then(||view!{
                                <p class="library-status">"Explicit checkout fetches the current remote branch and preserves the entire displaced local branch. The earlier failed pull does not pin a future remote head."</p>
                                <label><input type="checkbox" prop:checked=move||checkout.get() on:change=move|e|checkout.set(event_target_checked(&e))/>"Preserve my local branch and explicitly select the fetched remote branch"</label>
                                <button class="text-tool" disabled=move||busy.get()||locked.get()||!checkout.get()||review.with(|r|r.as_ref().is_none_or(|r|!review_is_current(r,saved.get().as_ref(),selected.get().as_ref()))) on:click=move|_|accept_checkout.run(())>"Pull with explicit checkout"</button>
                            })}
                        </section>}
                    })}
                    <p class="helper">"Pull saves into the local library. Reopen it there to review; it does not replace the working editor."</p>
                </details>
            </Show>
            <details><summary>"Retained local branches"</summary>
                <p class="helper">"Offline repository bookmarks, across all local plasmids. Reviewing these records makes no network request and does not change the working editor."</p>
                <button class="text-tool" disabled=move||busy.get()||locked.get() on:click=move|_|{bookmark_offset.set(0);send.run(Action::command(wire::Command::ListBookmarks{page:wire::LocalPageRequest::default(),expected_checkpoint:None}));}>"Refresh retained branches"</button>
                <ul class="lab-library-entries">{move||bookmarks.get().map(|page|page.items.into_vec().into_iter().map(|bookmark|match bookmark{
                    wire::Bookmark::Snapshot{local_document,revision}=>view!{<li><div><strong>"Saved catalog snapshot"</strong><p class="library-status">{format!("{}@{}",local_document.as_str(),revision.as_str())}</p><p class="helper">"Open this local plasmid’s Revision history to inspect this exact saved revision."</p></div></li>}.into_any(),
                    wire::Bookmark::NativeBranch{document,revision,event,undo_count,redo_count,incoming_revision,incoming_event}=>view!{<li><div><strong>"Retained native navigation branch"</strong><p class="library-status">{format!("Native document {} · revision {} · event {}",document.as_str(),revision.as_str(),event.as_ref().map(|e|e.as_str()).unwrap_or("none"))}</p><p class="library-status">{format!("{} undo / {} redo entries · incoming revision {} · event {}",undo_count,redo_count,incoming_revision.as_str(),incoming_event.as_ref().map(|e|e.as_str()).unwrap_or("none"))}</p><p class="helper">"This bookmark names native history rather than a catalog record. Review the relevant pre-pull archive through Local library Revision history."</p></div></li>}.into_any(),
                }).collect_view())}</ul>
                <Show when=move||bookmarks.with(|p|p.as_ref().is_some_and(|p|p.items.is_empty()))><p class="helper">"No retained branch bookmarks in this repository."</p></Show>
                <div class="library-pages"><button class="text-tool" disabled=move||busy.get()||locked.get()||bookmark_offset.get()==0||bookmarks.with(Option::is_none) on:click=move|_|if let Some(page)=bookmarks.get_untracked(){let offset=bookmark_offset.get_untracked().saturating_sub(100);bookmark_offset.set(offset);send.run(Action::command(wire::Command::ListBookmarks{page:wire::LocalPageRequest{offset,limit:100},expected_checkpoint:page.checkpoint}));}>"Previous"</button>
                    <button class="text-tool" disabled=move||busy.get()||locked.get()||bookmarks.with(|p|p.as_ref().is_none_or(|p|p.next_offset.is_none())) on:click=move|_|if let Some(page)=bookmarks.get_untracked(){if let Some(offset)=page.next_offset{bookmark_offset.set(offset);send.run(Action::command(wire::Command::ListBookmarks{page:wire::LocalPageRequest{offset,limit:100},expected_checkpoint:page.checkpoint}));}}>"Next"</button></div>
            </details>
            <crate::native_publication_view::NativePublicationControls client=client.get_value() saved selected destination locked busy on_queued=Callback::new(move|()|{pending_offset.set(0);send.run(Action::command(wire::Command::Pending{page:wire::LocalPageRequest::default()}));})/>
            <details><summary>"Pending pushes"</summary>
                <button class="text-tool" disabled=move||busy.get()||locked.get() on:click=move|_|{pending_offset.set(0);send.run(Action::command(wire::Command::Pending{page:wire::LocalPageRequest::default()}));}>"Refresh pending"</button>
                <ul class="lab-library-entries">{move||pending.get().into_iter().map(|p|{
                    let resume=p.job.clone();let cancel=p.job.clone();let discard=p.job.clone();
                    let can_cancel=!p.attempted;let rejected=p.rejection.is_some();
                    let label=format!("{} · {}{}",p.job.as_str(),if p.attempted{"attempted"}else{"queued"},if rejected{" · rejected"}else{""});
                    let reason=p.rejection.map(|r|format!("{:?}: {}",r.code,r.message.as_str())).unwrap_or_default();
                    view!{<li><span>{label}</span><p class="helper">{reason}</p><details><summary>"Push actions"</summary>
                        <button class="text-tool" disabled=move||busy.get()||locked.get()||rejected on:click=move|_|send.run(Action::command(wire::Command::ResumePush{job:resume.clone()}))>"Resume / check receipt"</button>
                        <button class="text-tool" disabled=move||busy.get()||locked.get()||!can_cancel on:click=move|_|send.run(Action::command(wire::Command::CancelUnattempted{job:cancel.clone()}))>"Cancel unattempted"</button>
                        <button class="text-tool" disabled=move||busy.get()||locked.get()||!rejected on:click=move|_|send.run(Action::command(wire::Command::DiscardRejected{job:discard.clone()}))>"Discard rejected"</button>
                    </details></li>}
                }).collect_view()}</ul>
                <p class="helper">"An interrupted request may already have committed remotely. Resume checks its durable receipt; attempted pushes cannot be cancelled."</p>
                <button class="text-tool" disabled=move||busy.get()||locked.get()||pending_offset.get()==0 on:click=move|_|{let offset=pending_offset.get_untracked().saturating_sub(100);pending_offset.set(offset);send.run(Action::command(wire::Command::Pending{page:wire::LocalPageRequest{offset,limit:100}}));}>"Previous"</button>
                <button class="text-tool" disabled=move||busy.get()||locked.get()||pending_next.get().is_none() on:click=move|_|if let Some(offset)=pending_next.get_untracked(){pending_offset.set(offset);send.run(Action::command(wire::Command::Pending{page:wire::LocalPageRequest{offset,limit:100}}));}>"Next"</button>
            </details>
            <details><summary>"Folder clones"</summary>
                <button class="text-tool" disabled=move||busy.get()||locked.get() on:click=move|_|{tree_offset.set(0);send.run(Action::command(wire::Command::TreeJobs{page:wire::LocalPageRequest::default()}));}>"Refresh folder clones"</button>
                <ul class="lab-library-entries">{move||trees.get().into_iter().map(|t|{
                    let resume=t.job.clone();let retry=t.job.clone();let doc_metadata=t.job.clone();let folder_metadata=t.job.clone();
                    let stopped=matches!(t.phase,wire::TreePhase::Complete|wire::TreePhase::Blocked|wire::TreePhase::Invalidated);
                    let blocked=t.phase==wire::TreePhase::Blocked;
                    let label=format!("{} · {:?}: {} / {} nodes saved",t.job.as_str(),t.phase,u32::from(t.folders_created)+u32::from(t.documents_created),t.scanned_nodes);
                    let issue=t.issue.map(|i|format!("{i:?}")).unwrap_or_default();
                    view!{<li><span>{label}</span><p class="helper">{issue}</p><details><summary>"Clone actions"</summary>
                        <button class="text-tool" disabled=move||busy.get()||locked.get()||stopped on:click=move|_|send.run(Action::command(wire::Command::ResumeTreeClone{job:resume.clone(),max_steps:1}))>"Resume one step"</button>
                        <button class="text-tool" disabled=move||busy.get()||locked.get()||!blocked on:click=move|_|send.run(Action::command(wire::Command::RetryTreeClone{job:retry.clone()}))>"Retry blocked clone"</button>
                        <button class="text-tool" disabled=move||busy.get()||locked.get()||saved.with(|s|s.as_ref().is_none_or(|s|s.dirty)) on:click=move|_|if let Some(s)=saved.get_untracked(){if !s.dirty{send.run(Action::command(wire::Command::QueueTreeDocumentMetadata{job:doc_metadata.clone(),local_document:s.id}));}}>"Queue saved local plasmid metadata"</button>
                        <button class="text-tool" disabled=move||busy.get()||locked.get()||destination.get().is_none() on:click=move|_|if let Some(local_folder)=destination.get_untracked(){send.run(Action::command(wire::Command::QueueTreeFolderMetadata{job:folder_metadata.clone(),local_folder}));}>"Queue selected local folder metadata"</button>
                    </details></li>}
                }).collect_view()}</ul>
                <p class="helper">"Clones preserve their scanned revisions and completed local edits. Metadata actions require membership in this clone and never turn its local mount into a remote move."</p>
                <button class="text-tool" disabled=move||busy.get()||locked.get()||tree_offset.get()==0 on:click=move|_|{let offset=tree_offset.get_untracked().saturating_sub(100);tree_offset.set(offset);send.run(Action::command(wire::Command::TreeJobs{page:wire::LocalPageRequest{offset,limit:100}}));}>"Previous"</button>
                <button class="text-tool" disabled=move||busy.get()||locked.get()||tree_next.get().is_none() on:click=move|_|if let Some(offset)=tree_next.get_untracked(){tree_offset.set(offset);send.run(Action::command(wire::Command::TreeJobs{page:wire::LocalPageRequest{offset,limit:100}}));}>"Next"</button>
            </details>
            <p class="library-status" role="status" aria-live="polite">{move||status.get()}</p>
        </details>
        </Show>
    }
}
