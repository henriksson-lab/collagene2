//! Shared folders live beside local documents. Capturing/importing native history
//! is delegated to App; this view never reconstructs DNA from display metadata.
use crate::{
    library_bindings::{Binding, Bindings, Remote},
    library_client::Client,
};
use collagene_library_protocol as wire;
use leptos::prelude::*;
use std::{cell::Cell, future::Future, pin::Pin, rc::Rc};
use wasm_bindgen_futures::spawn_local;

pub type LocalResult<T> = Pin<Box<dyn Future<Output = Result<T, String>>>>;

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct LocalDocument {
    pub id: String,
    pub name: String,
    pub revision: wire::Digest,
    pub head_event: Option<String>,
    pub workspace_epoch: u64,
}
pub struct Captured {
    pub local: LocalDocument,
    pub payload: wire::Payload,
}
pub struct Opened {
    pub document: wire::Document,
    pub revision: wire::RevisionSummary,
    pub payload: wire::Payload,
    pub remote: Remote,
    /// App must check immediately before committing an asynchronous import.
    pub guard: OpenGuard,
}

#[derive(Clone)]
pub struct OpenGuard {
    life: Rc<Cell<u64>>,
    sequence: u64,
}
impl OpenGuard {
    pub fn is_current(&self) -> bool {
        self.life.get() == self.sequence
    }
}

#[derive(Clone, Default)]
struct Navigation {
    path: Vec<wire::Folder>,
    query: String,
    offset: u32,
    generation: Option<u64>,
    completion: Option<String>,
}
impl Navigation {
    fn parent(&self) -> Option<wire::FolderId> {
        self.path.last().map(|folder| folder.id.clone())
    }
    fn request(&self) -> Result<wire::Request, String> {
        let page = wire::PageRequest {
            offset: self.offset,
            limit: 100,
            expected_catalog_generation: self.generation,
        };
        if self.query.trim().is_empty() {
            Ok(wire::Request::Browse {
                parent: self.parent(),
                page,
            })
        } else {
            Ok(wire::Request::Search {
                query: wire::Text::new(self.query.trim()).map_err(|error| error.to_string())?,
                scope: wire::SearchScope::AllFolders,
                page,
            })
        }
    }
    fn refresh(&self, message: String) -> Self {
        Self {
            offset: 0,
            generation: None,
            completion: Some(message),
            ..self.clone()
        }
    }
}

#[derive(Clone)]
enum Command {
    Browse(Navigation),
    Request(wire::Request, String),
    Open(wire::Document, Option<wire::RevisionId>),
    Publish,
}
#[derive(Clone)]
struct HistoryPage {
    document: wire::Document,
    revisions: Vec<wire::RevisionSummary>,
    offset: u32,
    generation: u64,
    next_offset: Option<u32>,
}

#[derive(Clone)]
enum NameEditor {
    Folder,
    Rename(wire::Entry),
}

pub fn enabled() -> bool {
    if crate::desktop_rpc::is_desktop() {
        return true;
    }
    web_sys::window()
        .and_then(|window| window.document())
        .and_then(|document| {
            document
                .query_selector("meta[name='collagene-library'][content='enabled']")
                .ok()
                .flatten()
        })
        .is_some()
}

#[component]
pub fn LibraryPane(
    client: Client,
    current: Memo<Option<LocalDocument>>,
    bindings: Memo<Bindings>,
    capture: Callback<(), LocalResult<Captured>>,
    on_open: Callback<Opened, LocalResult<LocalDocument>>,
    on_binding: Callback<(u64, Binding), Result<(), String>>,
    locked: Memo<bool>,
) -> impl IntoView {
    let local = client.is_local();
    let initial_context = client.context();
    let server_key = client.remote_key();
    let session = StoredValue::new_local(client);
    let life = StoredValue::new_local(Rc::new(Cell::new(0u64)));
    let busy = RwSignal::new(false);
    let connected = RwSignal::new(initial_context);
    let key = RwSignal::new(String::new());
    let status = RwSignal::new(String::new());
    let navigation = RwSignal::new(Navigation::default());
    let entries = RwSignal::new(Vec::<wire::Entry>::new());
    let next_offset = RwSignal::new(None::<u32>);
    let search = RwSignal::new(String::new());
    let editor = RwSignal::new(None::<NameEditor>);
    let name = RwSignal::new(String::new());
    let moving = RwSignal::new(None::<wire::Entry>);
    let history = RwSignal::new(None::<HistoryPage>);
    let runner = StoredValue::new_local(None::<Callback<Command>>);
    let origin = StoredValue::new(server_key);
    on_cleanup(move || {
        life.with_value(|life| life.set(life.get().wrapping_add(1)));
        session.with_value(|session| session.disconnect());
    });
    let remote = Memo::new(move |_| {
        connected.get().map(|context| Remote {
            server: origin.get_value(),
            lab: context.lab_id,
        })
    });
    let current_binding = Memo::new(move |_| {
        let remote = remote.get()?;
        let local = current.get()?;
        bindings.with(|bindings| bindings.get(&remote, &local.id).cloned())
    });
    let send = Callback::new(move |command: Command| {
        if busy.get_untracked() || locked.get_untracked() {
            return;
        }
        let client = session.get_value();
        let life = life.get_value();
        let sequence = life.get().wrapping_add(1);
        life.set(sequence);
        busy.set(true);
        status.set(
            match &command {
                Command::Publish => "Preparing plasmid history…",
                Command::Open(..) => "Opening library revision…",
                _ => "Loading library…",
            }
            .into(),
        );
        let nav = navigation.get_untracked();
        spawn_local(async move {
            let result = async {
                let mut captured = None;
                let request = match &command {
                    Command::Browse(target) => target.request()?,
                    Command::Request(request, _) => request.clone(),
                    Command::Open(document, revision) => wire::Request::OpenDocument {
                        id: document.id.clone(),
                        revision: revision.clone(),
                    },
                    Command::Publish => {
                        let value = capture.run(()).await?;
                        if life.get() != sequence {
                            return Err("Library session changed.".into());
                        }
                        let selected = value
                            .payload
                            .selection
                            .as_ref()
                            .ok_or("A shared plasmid must retain its native history.")?;
                        if value.payload.kind != wire::PayloadKind::WorkspaceArchive
                            || value.payload.sharing
                                != Some(wire::WorkspaceSharing::DocumentClosure)
                            || selected.document_id.as_str() != value.local.id
                            || selected.revision != value.local.revision
                            || selected.head_event.as_ref().map(|id| id.as_str())
                                != value.local.head_event.as_deref()
                        {
                            return Err(
                                "Captured plasmid does not match its published history pin.".into(),
                            );
                        }
                        let remote = remote
                            .get_untracked()
                            .ok_or("Connect to the lab library first.")?;
                        let binding = bindings.with_untracked(|bindings| {
                            bindings.get(&remote, &value.local.id).cloned()
                        });
                        let request = if let Some(binding) = binding {
                            wire::Request::SaveRevision {
                                id: binding.remote_document,
                                base_revision: binding.remote_revision,
                                payload: value.payload,
                            }
                        } else {
                            wire::Request::ImportDocument {
                                folder: nav.parent(),
                                name: wire::Text::new(&value.local.name)
                                    .map_err(|error| error.to_string())?,
                                payload: value.payload,
                            }
                        };
                        captured = Some(value.local);
                        request
                    }
                };
                let reply = client.request(request).await?;
                if life.get() != sequence {
                    return Err("Library session changed.".into());
                }
                connected.set(client.context());
                let success = match reply.outcome {
                    wire::Outcome::Success(success) => success,
                    wire::Outcome::Failure(failure) => {
                        return Err(format!(
                            "{}{}",
                            if failure.code == wire::ErrorCode::Conflict {
                                "Conflict: "
                            } else {
                                ""
                            },
                            failure.message.as_str()
                        ));
                    }
                };
                match success {
                    wire::Success::Browsed(page) | wire::Success::Searched(page) => {
                        let Command::Browse(mut target) = command else {
                            return Err("Unexpected library listing.".into());
                        };
                        target.generation = Some(page.catalog_generation);
                        let completion = target.completion.take();
                        navigation.set(target);
                        next_offset.set(page.next_offset);
                        entries.set(page.entries.iter().cloned().collect());
                        status.set(completion.unwrap_or_else(|| "Library is up to date.".into()));
                    }
                    wire::Success::Opened {
                        document,
                        revision,
                        payload,
                    } => {
                        let remote_revision = revision.id.clone();
                        let remote_document = document.id.clone();
                        let opened_remote =
                            remote.get_untracked().ok_or("Library session changed.")?;
                        let local = on_open
                            .run(Opened {
                                document,
                                revision,
                                payload,
                                remote: opened_remote,
                                guard: OpenGuard {
                                    life: life.clone(),
                                    sequence,
                                },
                            })
                            .await?;
                        if life.get() != sequence {
                            return Err("Library session changed.".into());
                        }
                        let remote = remote.get_untracked().ok_or("Library session changed.")?;
                        on_binding.run((
                            local.workspace_epoch,
                            Binding {
                                remote,
                                local_document: local.id,
                                local_revision: local.revision,
                                local_head_event: local.head_event,
                                remote_document,
                                remote_revision,
                            },
                        ))?;
                        status.set(format!("Opened {}.", local.name));
                    }
                    wire::Success::Imported(document) | wire::Success::Saved { document, .. } => {
                        let local =
                            captured.ok_or("Library save has no captured document receipt.")?;
                        let remote = remote.get_untracked().ok_or("Library session changed.")?;
                        on_binding.run((
                            local.workspace_epoch,
                            Binding {
                                remote,
                                local_document: local.id,
                                local_revision: local.revision,
                                local_head_event: local.head_event,
                                remote_document: document.id,
                                remote_revision: document.head.id,
                            },
                        ))?;
                        busy.set(false);
                        runner.get_value().expect("installed library handler").run(
                            Command::Browse(nav.refresh(format!(
                                "Saved {}. Local edits made during the save remain local.",
                                local.name
                            ))),
                        );
                        return Ok(());
                    }
                    wire::Success::Revisions(page) => {
                        let Command::Request(
                            wire::Request::ListRevisions {
                                id,
                                page: requested,
                            },
                            _,
                        ) = &command
                        else {
                            return Err("Unexpected library history.".into());
                        };
                        let document = history
                            .with_untracked(|history| {
                                history
                                    .as_ref()
                                    .filter(|history| &history.document.id == id)
                                    .map(|history| history.document.clone())
                            })
                            .or_else(|| {
                                entries.with_untracked(|entries| {
                                    entries.iter().find_map(|entry| match entry {
                                        wire::Entry::Document(document) if &document.id == id => {
                                            Some(document.clone())
                                        }
                                        _ => None,
                                    })
                                })
                            })
                            .ok_or(
                                "The document is no longer in this listing; refresh the folder.",
                            )?;
                        history.set(Some(HistoryPage {
                            document,
                            revisions: page.revisions.iter().cloned().collect(),
                            offset: requested.offset,
                            generation: page.catalog_generation,
                            next_offset: page.next_offset,
                        }));
                        status.set(
                            "Choose an immutable revision to open. Local changes remain protected."
                                .into(),
                        );
                    }
                    _ => {
                        editor.set(None);
                        moving.set(None);
                        let message = if let Command::Request(_, message) = command {
                            message
                        } else {
                            "Library updated.".into()
                        };
                        busy.set(false);
                        runner
                            .get_value()
                            .expect("installed library handler")
                            .run(Command::Browse(nav.refresh(message)));
                        return Ok(());
                    }
                }
                Ok(())
            }
            .await;
            if life.get() != sequence {
                return;
            }
            busy.set(false);
            if let Err(error) = result {
                connected.set(client.context());
                status.set(error);
            }
        });
    });
    runner.set_value(Some(send));
    if local || connected.with_untracked(Option::is_some) {
        send.run(Command::Browse(Navigation::default()));
    }
    let connect = Callback::new(move |()| {
        let result = session.with_value(|session| session.connect(key.get_untracked()));
        key.set(String::new());
        match result {
            Ok(()) => {
                connected.set(None);
                entries.set(vec![]);
                editor.set(None);
                name.set(String::new());
                search.set(String::new());
                moving.set(None);
                history.set(None);
                send.run(Command::Browse(Navigation::default()));
            }
            Err(error) => status.set(error),
        }
    });
    let disconnect = move |_| {
        life.with_value(|life| life.set(life.get().wrapping_add(1)));
        session.with_value(|session| session.disconnect());
        connected.set(None);
        editor.set(None);
        name.set(String::new());
        search.set(String::new());
        entries.set(vec![]);
        history.set(None);
        moving.set(None);
        busy.set(false);
        key.set(String::new());
        status.set("Disconnected. Local documents and unpublished edits are retained.".into());
    };
    let submit_name = move || {
        let result = wire::Text::new(name.get_untracked())
            .map_err(|error| error.to_string())
            .and_then(|name| {
                wire::validate_name(&name).map_err(|error| error.to_string())?;
                let request = match editor.get_untracked().ok_or("No name action is active")? {
                    NameEditor::Folder => wire::Request::CreateFolder {
                        parent: navigation.get_untracked().parent(),
                        name,
                    },
                    NameEditor::Rename(wire::Entry::Folder(folder)) => {
                        wire::Request::RenameFolder {
                            id: folder.id,
                            expected_generation: folder.generation,
                            name,
                        }
                    }
                    NameEditor::Rename(wire::Entry::Document(document)) => {
                        wire::Request::RenameDocument {
                            id: document.id,
                            expected_generation: document.generation,
                            name,
                        }
                    }
                };
                Ok(request)
            });
        match result {
            Ok(request) => send.run(Command::Request(request, "Library name saved.".into())),
            Err(error) => status.set(error),
        }
    };
    view! {
        <details class="lab-library-pane"><summary>{if local { "Local library" } else { "Lab library" }}</summary>
            <Show when=move||!local&&connected.with(Option::is_none)>
                <form class="library-connect" on:submit=move|event|{event.prevent_default();connect.run(());}>
                    <label>"Lab access key"<input type="password" autocomplete="off" aria-label="Lab library access key" prop:value=move||key.get() disabled=move||busy.get() on:input=move|event|key.set(event_target_value(&event))/></label>
                    <button class="button secondary" type="submit" disabled=move||busy.get()||locked.get()>"Connect"</button>
                    <p class="helper">"The key is kept in memory for this session. Local documents stay on this device."</p>
                </form>
            </Show>
            <Show when=move||connected.with(Option::is_some)>
                <div class="library-heading"><strong>{move||connected.get().map(|context|context.name.into_string()).unwrap_or_default()}</strong><Show when=move||!local><button class="text-tool" on:click=disconnect>"Disconnect"</button></Show></div>
                <div class="library-location"><button class="text-tool" disabled=move||busy.get() on:click=move|_|{search.set(String::new());send.run(Command::Browse(Navigation::default()));}>"Library"</button>
                    {move||navigation.with(|nav|nav.path.iter().enumerate().map(|(index,folder)|{let path=nav.path[..=index].to_vec();let label=folder.name.to_string();view!{<span>" / "</span><button class="text-tool" disabled=move||busy.get() on:click=move|_|{search.set(String::new());send.run(Command::Browse(Navigation{path:path.clone(),..Default::default()}));}>{label}</button>}}).collect_view())}
                </div>
                <form class="library-search" on:submit=move|event|{event.prevent_default();send.run(Command::Browse(Navigation{query:search.get_untracked(),..Default::default()}));}>
                    <input type="search" aria-label="Search lab plasmids and folders" placeholder={if local { "Search local library" } else { "Search lab library" }} prop:value=move||search.get() on:input=move|event|search.set(event_target_value(&event))/><button type="submit" class="text-tool" disabled=move||busy.get()>"Find"</button>
                </form>
                <div class="library-actions"><button class="text-tool" disabled=move||busy.get()||locked.get() on:click=move|_|{name.set(String::new());editor.set(Some(NameEditor::Folder));}>"New folder"</button>
                    <button class="text-tool" disabled=move||busy.get()||locked.get()||current.with(Option::is_none) on:click=move|_|send.run(Command::Publish)>{move||if local { if current_binding.get().is_some(){"Save plasmid revision"}else{"Save current plasmid"} } else if current_binding.get().is_some(){"Push plasmid revision"}else{"Share current plasmid"}}</button>
                    <button class="text-tool" disabled=move||busy.get() on:click=move|_|send.run(Command::Browse(navigation.get_untracked().refresh("Library refreshed.".into())))>"Refresh"</button>
                </div>
                {move||editor.get().map(|_|view!{<form class="library-name-editor" on:submit=move|event|{event.prevent_default();submit_name();}><label>"Name"<input autofocus=true prop:value=move||name.get() on:input=move|event|name.set(event_target_value(&event))/></label><button class="text-tool" disabled=move||busy.get()||locked.get() type="submit">"Save"</button><button class="text-tool" type="button" on:click=move|_|editor.set(None)>"Cancel"</button></form>})}
                {move||moving.get().map(|entry|view!{<div class="library-move"><p>{format!("Choose a destination for {}.",entry.name())}</p><button class="text-tool" disabled=move||busy.get()||locked.get() on:click=move|_|{let parent=navigation.get_untracked().parent();let request=match entry.clone(){wire::Entry::Folder(folder)=>wire::Request::MoveFolder{id:folder.id,expected_generation:folder.generation,parent},wire::Entry::Document(document)=>wire::Request::MoveDocument{id:document.id,expected_generation:document.generation,folder:parent}};send.run(Command::Request(request,"Item moved.".into()));}>"Move here"</button><button class="text-tool" on:click=move|_|moving.set(None)>"Cancel move"</button></div>})}
                <ul class="lab-library-entries" aria-label={if local { "Local plasmids and folders" } else { "Shared plasmids and folders" }}>{move||entries.get().into_iter().map(|entry|{
                    let open=entry.clone();let rename=entry.clone();let move_entry=entry.clone();let history_entry=entry.clone();let delete=entry.clone();let is_folder=matches!(entry,wire::Entry::Folder(_));let label=entry.name().to_owned();
                    view!{<li class:library-folder=is_folder><button class="library-entry-name" disabled=move||busy.get()||locked.get() on:click=move|_|match &open{wire::Entry::Folder(folder)=>{let mut path=if navigation.with_untracked(|nav|nav.query.is_empty()){navigation.get_untracked().path}else{vec![]};path.push(folder.clone());search.set(String::new());send.run(Command::Browse(Navigation{path,..Default::default()}));},wire::Entry::Document(document)=>send.run(Command::Open(document.clone(),None))}><span aria-hidden="true">{if is_folder{"▸"}else{"○"}}</span>{label}</button>
                        <details class="library-entry-menu" data-transient-menu=""><summary aria-label="Item actions">"⋯"</summary><div>
                            <button class="text-tool" disabled=move||busy.get()||locked.get() on:click=move|_|{name.set(rename.name().to_owned());editor.set(Some(NameEditor::Rename(rename.clone())));}>"Rename"</button>
                            <button class="text-tool" disabled=move||busy.get()||locked.get() on:click=move|_|moving.set(Some(move_entry.clone()))>"Move…"</button>
                            {(!is_folder).then(||view!{<button class="text-tool" disabled=move||busy.get() on:click=move|_|if let wire::Entry::Document(document)=&history_entry{send.run(Command::Request(wire::Request::ListRevisions{id:document.id.clone(),page:Default::default()},String::new()));}>"Revision history"</button>})}
                            {is_folder.then(||view!{<button class="text-tool" disabled=move||busy.get()||locked.get() on:click=move|_|if let wire::Entry::Folder(folder)=&delete{send.run(Command::Request(wire::Request::DeleteEmptyFolder{id:folder.id.clone(),expected_generation:folder.generation},"Empty folder removed.".into()));}>"Delete empty folder"</button>})}
                        </div></details>
                    </li>}
                }).collect_view()}</ul>
                <Show when=move||!busy.get()&&entries.with(Vec::is_empty)><p class="helper">"No plasmids or folders in this listing."</p></Show>
                <div class="library-pages"><button class="text-tool" disabled=move||busy.get()||navigation.with(|nav|nav.offset==0) on:click=move|_|{let mut nav=navigation.get_untracked();nav.offset=nav.offset.saturating_sub(100);send.run(Command::Browse(nav));}>"Previous"</button><button class="text-tool" disabled=move||busy.get()||next_offset.with(Option::is_none) on:click=move|_|if let Some(offset)=next_offset.get_untracked(){let mut nav=navigation.get_untracked();nav.offset=offset;send.run(Command::Browse(nav));}>"Next"</button></div>
                {move||history.get().map(|listing|{
                    let previous=listing.clone(); let next=listing.clone();
                    view!{<section class="library-revisions">
                        <div class="library-heading"><strong>{format!("{} revisions",listing.document.name.as_str())}</strong><button class="text-tool" on:click=move|_|history.set(None)>"Close"</button></div>
                        <ul>{listing.revisions.into_iter().map(|revision|{let document=listing.document.clone();view!{<li><button class="text-tool" disabled=move||busy.get()||locked.get() on:click=move|_|send.run(Command::Open(document.clone(),Some(revision.id.clone())))>{format!("Revision {} · {} bytes",revision.number,revision.archive_size)}</button></li>}}).collect_view()}</ul>
                        <div class="library-pages">
                            <button class="text-tool" disabled=move||busy.get()||previous.offset==0 on:click=move|_|send.run(Command::Request(wire::Request::ListRevisions{id:previous.document.id.clone(),page:wire::PageRequest{offset:previous.offset.saturating_sub(100),limit:100,expected_catalog_generation:Some(previous.generation)}},String::new()))>"Previous revisions"</button>
                            <button class="text-tool" disabled=move||busy.get()||next.next_offset.is_none() on:click=move|_|if let Some(offset)=next.next_offset{send.run(Command::Request(wire::Request::ListRevisions{id:next.document.id.clone(),page:wire::PageRequest{offset,limit:100,expected_catalog_generation:Some(next.generation)}},String::new()));}>"Older revisions"</button>
                        </div>
                    </section>}
                })}
            </Show>
            <p class="library-status" role="status" aria-live="polite">{move||status.get()}</p>
        </details>
    }
}
