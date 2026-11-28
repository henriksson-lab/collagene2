//! Shared folders live beside local documents. Capturing/importing native history
//! is delegated to App; this view never reconstructs DNA from display metadata.
use crate::library::library_bindings::Binding;
use crate::library::library_bindings::Bindings;
use crate::library::library_bindings::Remote;
use crate::library::library_client::Client;
use collagene_library_protocol as wire;
use leptos::prelude::*;
use std::{cell::Cell, future::Future, pin::Pin, rc::Rc};
use wasm_bindgen_futures::{JsFuture, spawn_local};

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
    private_local: bool,
}
impl OpenGuard {
    pub fn is_private_local(&self) -> bool {
        self.private_local
    }
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
    Publish(crate::library::library_documents::CaptureRequest),
}
#[derive(Clone)]
struct PendingPublication {
    capture: crate::library::library_documents::CaptureReview,
    context: wire::LabContext,
    parent: Option<wire::FolderId>,
    sequence: u64,
    session_epoch: u64,
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
    if crate::platform::desktop_rpc::is_desktop() {
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
    capture: Callback<
        crate::library::library_documents::CaptureRequest,
        LocalResult<crate::library::library_documents::CaptureOutcome>,
    >,
    on_open: Callback<Opened, LocalResult<LocalDocument>>,
    on_binding: Callback<(u64, Binding), Result<(), String>>,
    locked: Memo<bool>,
    #[prop(optional)] refresh: Option<RwSignal<u64>>,
    #[prop(optional)] on_folder: Option<Callback<Option<wire::FolderId>>>,
) -> impl IntoView {
    let local = client.is_local();
    let native_local = local && crate::platform::desktop_rpc::is_desktop();
    let initial_context = client.context();
    let server_key = client.remote_key();
    let session = StoredValue::new_local(client);
    let life = StoredValue::new_local(Rc::new(Cell::new(0u64)));
    let busy = RwSignal::new(false);
    let connected = RwSignal::new(initial_context);
    let key = RwSignal::new(String::new());
    let status = RwSignal::new(String::new());
    let copy_mcp = Callback::new(move |document: wire::Document| {
        if !native_local {
            return;
        }
        let Some(window) = web_sys::window() else {
            status.set("System clipboard is unavailable.".into());
            return;
        };
        if !window.is_secure_context() {
            status.set("System clipboard access requires a secure application window.".into());
            return;
        }
        let clipboard = window.navigator().clipboard();
        if clipboard.is_undefined() || clipboard.is_null() {
            status.set("System clipboard is unavailable in this application window.".into());
            return;
        }
        let reference = format!("{}@{}", document.id.as_str(), document.head.id.as_str());
        // As in ClipboardTools, start clipboard access synchronously in the
        // actual click gesture, then await its result without retaining DNA.
        let promise = clipboard.write_text(&reference);
        let life = life.get_value();
        let epoch = life.get();
        spawn_local(async move {
            let result = JsFuture::from(promise).await;
            if life.get() != epoch {
                return;
            }
            status.set(if result.is_ok() {
                "MCP reference copied.".into()
            } else {
                "System clipboard access failed; the MCP reference was not copied.".into()
            });
        });
    });
    let navigation = RwSignal::new(Navigation::default());
    let entries = RwSignal::new(Vec::<wire::Entry>::new());
    let next_offset = RwSignal::new(None::<u32>);
    let search = RwSignal::new(String::new());
    Effect::new(move |_| {
        if local {
            if let Some(callback) = on_folder {
                callback.run(navigation.with(Navigation::parent));
            }
        }
    });
    let editor = RwSignal::new(None::<NameEditor>);
    let name = RwSignal::new(String::new());
    let moving = RwSignal::new(None::<wire::Entry>);
    let history = RwSignal::new(None::<HistoryPage>);
    let publication = RwSignal::new(None::<PendingPublication>);
    let publication_choices = RwSignal::new(String::new());
    let sequence_requirements = RwSignal::new(
        None::<(
            LocalDocument,
            Vec<collagene_library_sharing::sbol_subset::ChoiceRequirement>,
        )>,
    );
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
        if let Command::Publish(request) = &command {
            if request.approved.is_some() {
                let valid = publication.with_untracked(|pending| {
                    pending.as_ref().is_some_and(|pending| {
                        current.get_untracked().as_ref() == Some(&pending.capture.local)
                            && client.context().as_ref() == Some(&pending.context)
                            && client.publication_epoch() == pending.session_epoch
                            && life.with_value(|life| life.get()) == pending.sequence
                            && navigation.with_untracked(Navigation::parent) == pending.parent
                    })
                });
                if !valid {
                    publication.set(None);
                    status.set("Publication source or library destination changed; review it again. Nothing was sent.".into());
                    return;
                }
            }
        }
        publication.set(None);
        sequence_requirements.set(None);
        let context_before = client.context();
        let session_epoch_before = client.publication_epoch();
        let life = life.get_value();
        let sequence = life.get().wrapping_add(1);
        life.set(sequence);
        busy.set(true);
        status.set(
            match &command {
                Command::Publish(_) => "Preparing plasmid history…",
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
                    Command::Publish(request) => {
                        let outcome = capture.run(request.clone()).await?;
                        if client.context() != context_before || client.publication_epoch() != session_epoch_before {
                            return Err("Library login changed while preparing publication. Nothing was sent.".into());
                        }
                        if life.get() != sequence {
                            return Err("Library session changed.".into());
                        }
                        let value = match outcome {
                            crate::library::library_documents::CaptureOutcome::Ready(value) => value,
                            crate::library::library_documents::CaptureOutcome::ChoicesRequired { local, requirements } => {
                                if current.get_untracked().as_ref() != Some(&local) {
                                    return Err("Plasmid changed while inspecting SBOL choices. Nothing was sent.".into());
                                }
                                sequence_requirements.set(Some((local, requirements)));
                                status.set("Reachable SBOL components need explicit sequence choices. Fill the declarations below, then prepare publication again. Nothing was sent.".into());
                                return Ok(());
                            }
                            crate::library::library_documents::CaptureOutcome::ReviewRequired(review) => {
                                if current.get_untracked().as_ref() != Some(&review.local) {
                                    return Err("Plasmid changed while preparing documentary review. Nothing was sent.".into());
                                }
                                let context = client.context().ok_or("Connect to the library before reviewing publication.")?;
                                publication.set(Some(PendingPublication { capture: review, context, parent: nav.parent(), sequence, session_epoch: session_epoch_before }));
                                status.set("Review the documentary subset below. No archive has been encoded or sent.".into());
                                return Ok(());
                            }
                        };
                        let selected = value
                            .payload
                            .selection
                            .as_ref()
                            .ok_or("A shared plasmid must retain its native history.")?;
                        if value.payload.kind != wire::PayloadKind::WorkspaceArchive
                            || value.payload.sharing
                                != Some(if client.is_local() {
                                    wire::WorkspaceSharing::ExplicitWholeWorkspace
                                } else {
                                    wire::WorkspaceSharing::DocumentClosure
                                })
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
                                    private_local: client.is_local(),
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
    let observed_refresh = StoredValue::new_local(0u64);
    Effect::new(move |_| {
        let Some(refresh) = refresh else {
            return;
        };
        let generation = refresh.get();
        if generation == observed_refresh.get_value()
            || busy.get()
            || locked.get()
            || connected.get().is_none()
        {
            return;
        }
        observed_refresh.set_value(generation);
        let mut nav = navigation.get_untracked();
        nav.offset = 0;
        nav.generation = None;
        send.run(Command::Browse(nav));
    });
    Effect::new(move |_| {
        let selected = current.get();
        let context = connected.get();
        let blocked = locked.get();
        if publication.with(|pending| {
            pending.as_ref().is_some_and(|pending| {
                blocked
                    || selected.as_ref() != Some(&pending.capture.local)
                    || context.as_ref() != Some(&pending.context)
            })
        }) {
            publication.set(None);
        }
        if sequence_requirements.with(|requirements| {
            requirements
                .as_ref()
                .is_some_and(|(local, _)| blocked || selected.as_ref() != Some(local))
        }) {
            sequence_requirements.set(None);
        }
    });
    let prepare_publication = Callback::new(move |()| {
        if local {
            send.run(Command::Publish(
                crate::library::library_documents::CaptureRequest {
                    expected: current.get_untracked(),
                    ..Default::default()
                },
            ));
            return;
        }
        let text = publication_choices.get_untracked();
        if text.len() > 16 * 1024 {
            status.set("Alternative sequence declarations exceed 16 KiB.".into());
            return;
        }
        let choices = if text.trim().is_empty() {
            Ok(Default::default())
        } else {
            serde_json::from_str::<collagene_library_sharing::sbol_subset::PublicationChoices>(&text)
                .map_err(|_| "Use a JSON object mapping original file hash → component IRI → array of explicitly chosen sequence IRIs.".to_string())
        };
        match choices {
            Ok(choices) => send.run(Command::Publish(
                crate::library::library_documents::CaptureRequest {
                    choices,
                    expected: current.get_untracked(),
                    ..Default::default()
                },
            )),
            Err(error) => status.set(error),
        }
    });
    let approve_publication = Callback::new(move |()| {
        let Some(pending) = publication.get_untracked() else {
            return;
        };
        let approved = pending.capture.review.approve();
        send.run(Command::Publish(
            crate::library::library_documents::CaptureRequest {
                approved: Some(approved),
                expected: Some(pending.capture.local),
                ..Default::default()
            },
        ));
    });
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
                    <label>"Member access key"<input type="password" autocomplete="off" aria-label="Lab member access key" prop:value=move||key.get() disabled=move||busy.get() on:input=move|event|key.set(event_target_value(&event))/></label>
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
                    <button class="text-tool" disabled=move||busy.get()||locked.get()||current.with(Option::is_none) on:click=move|_|prepare_publication.run(())>{move||if local { if current_binding.get().is_some(){"Save plasmid revision"}else{"Save current plasmid"} } else if current_binding.get().is_some(){"Push plasmid revision"}else{"Share current plasmid"}}</button>
                    <Show when=move||!local>
                    <details><summary>"SBOL sequence alternatives"</summary><p class="helper">"Only needed when a reachable SBOL component has alternative sequences. Declare choices explicitly as original file hash → component IRI → sequence IRI array. No alternative is automatically selected."</p><textarea aria-label="Explicit SBOL sequence alternative declarations" maxlength="16384" prop:value=move||publication_choices.get() disabled=move||busy.get() on:input=move|event|{publication_choices.set(event_target_value(&event));publication.set(None);}/></details>
                    {move||sequence_requirements.get().map(|(_,requirements)|view!{<details open><summary>"Reachable SBOL sequence choices"</summary><p class="helper">"These identifiers belong to selected structural dependencies. For each component with no required sequence, explicitly declare one or more available sequences using the JSON field above. Required native source sequences must remain included. No first alternative is chosen."</p>{requirements.into_iter().map(|r|view!{<div><p>"Original file: "<code>{r.original_hash}</code></p><p>"Component: "<code>{r.component_iri}</code></p><p>{format!("Required by native source: {:?}",r.required_sequences)}</p><ul>{r.available_sequences.into_iter().map(|id|view!{<li><code>{id}</code></li>}).collect_view()}</ul></div>}).collect_view()}</details>})}
                    {move||publication.get().map(|pending|view!{<section class="library-publication-review" aria-label="Documentary subset publication review">
                        <h4>"Review SBOL documentary publication"</h4><p>{format!("Plasmid {} · native revision {}",pending.capture.local.name,pending.capture.local.revision.as_str())}</p>
                        <p class="helper">"The exact original SBOL file stays private in your local workspace. Only the reviewed structural subset and a counts-only omission receipt will accompany this plasmid. Native sequence identities, source history and physical verification are unchanged; this documentary projection adds no scientific authority. Retained annotations, descriptions and vocabulary remain included; review their content before sharing."</p>
                        {pending.capture.review.entries().iter().map(|entry|view!{<div><p>"Original file omitted: "<code>{entry.original_hash.clone()}</code></p><p>"Published subset: "<code>{entry.subset_hash.clone()}</code></p><p>{format!("{} retained objects · {} unrelated objects omitted · {} alternative sequence edges omitted",entry.retained_objects,entry.omissions.unrelated_objects,entry.omissions.alternative_sequence_edges)}</p></div>}).collect_view()}
                        <button class="button secondary" disabled=move||busy.get()||locked.get() on:click=move|_|approve_publication.run(())>"Approve subset and save plasmid"</button><button class="text-tool" disabled=move||busy.get() on:click=move|_|{publication.set(None);status.set("Documentary publication cancelled. Nothing was sent.".into());}>"Cancel"</button>
                    </section>})}
                    </Show>
                    <button class="text-tool" disabled=move||busy.get() on:click=move|_|send.run(Command::Browse(navigation.get_untracked().refresh("Library refreshed.".into())))>"Refresh"</button>
                </div>
                {move||editor.get().map(|_|view!{<form class="library-name-editor" on:submit=move|event|{event.prevent_default();submit_name();}><label>"Name"<input autofocus=true prop:value=move||name.get() on:input=move|event|name.set(event_target_value(&event))/></label><button class="text-tool" disabled=move||busy.get()||locked.get() type="submit">"Save"</button><button class="text-tool" type="button" on:click=move|_|editor.set(None)>"Cancel"</button></form>})}
                {move||moving.get().map(|entry|view!{<div class="library-move"><p>{format!("Choose a destination for {}.",entry.name())}</p><button class="text-tool" disabled=move||busy.get()||locked.get() on:click=move|_|{let parent=navigation.get_untracked().parent();let request=match entry.clone(){wire::Entry::Folder(folder)=>wire::Request::MoveFolder{id:folder.id,expected_generation:folder.generation,parent},wire::Entry::Document(document)=>wire::Request::MoveDocument{id:document.id,expected_generation:document.generation,folder:parent}};send.run(Command::Request(request,"Item moved.".into()));}>"Move here"</button><button class="text-tool" on:click=move|_|moving.set(None)>"Cancel move"</button></div>})}
                <ul class="lab-library-entries" aria-label={if local { "Local plasmids and folders" } else { "Shared plasmids and folders" }}>{move||entries.get().into_iter().map(|entry|{
                    let open=entry.clone();let rename=entry.clone();let move_entry=entry.clone();let history_entry=entry.clone();let mcp_entry=entry.clone();let delete=entry.clone();let is_folder=matches!(entry,wire::Entry::Folder(_));let label=entry.name().to_owned();
                    view!{<li class:library-folder=is_folder><button class="library-entry-name" disabled=move||busy.get()||locked.get() on:click=move|_|match &open{wire::Entry::Folder(folder)=>{let mut path=if navigation.with_untracked(|nav|nav.query.is_empty()){navigation.get_untracked().path}else{vec![]};path.push(folder.clone());search.set(String::new());send.run(Command::Browse(Navigation{path,..Default::default()}));},wire::Entry::Document(document)=>send.run(Command::Open(document.clone(),None))}><span aria-hidden="true">{if is_folder{"▸"}else{"○"}}</span>{label}</button>
                        <details class="library-entry-menu" data-transient-menu=""><summary aria-label="Item actions">"⋯"</summary><div>
                            <button class="text-tool" disabled=move||busy.get()||locked.get() on:click=move|_|{name.set(rename.name().to_owned());editor.set(Some(NameEditor::Rename(rename.clone())));}>"Rename"</button>
                            <button class="text-tool" disabled=move||busy.get()||locked.get() on:click=move|_|moving.set(Some(move_entry.clone()))>"Move…"</button>
                            {(!is_folder).then(||view!{<button class="text-tool" disabled=move||busy.get() on:click=move|_|if let wire::Entry::Document(document)=&history_entry{send.run(Command::Request(wire::Request::ListRevisions{id:document.id.clone(),page:Default::default()},String::new()));}>"Revision history"</button>})}
                            {(native_local&&!is_folder).then(||view!{<button class="text-tool" disabled=move||busy.get()||locked.get() on:click=move|_|if let wire::Entry::Document(document)=&mcp_entry{copy_mcp.run(document.clone());}>"Copy MCP reference"</button>})}
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
                        <ul>{listing.revisions.into_iter().map(|revision|{
                            let document=listing.document.clone();
                            let target=wire::discussion::Target{document:document.id.clone(),revision:revision.id.clone()};
                            view!{<li><button class="text-tool" disabled=move||busy.get()||locked.get() on:click=move|_|send.run(Command::Open(document.clone(),Some(revision.id.clone())))>{format!("Revision {} · {} bytes",revision.number,revision.archive_size)}</button>
                                {(!local).then(||view!{<crate::library::library_discussion_view::LibraryDiscussion target=target client=crate::library::library_discussion_client::Client::browser()/>})}
                            </li>}
                        }).collect_view()}</ul>
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
