//! Explicit native publication reviews. Queueing persists a job; it never
//! contacts a remote or replaces the working editor.
use crate::library::remote_library_view::SavedLocalDocument;
use crate::platform::desktop_sync_client::Client;
use crate::platform::native_publication_model as model;
use crate::platform::native_publication_model::Chosen;
use crate::platform::native_publication_model::Context;
use crate::platform::native_publication_model::SavedPin;
use collagene_desktop_sync_protocol::{self as wire, library as base};
use leptos::prelude::*;
use std::{cell::Cell, rc::Rc};
use wasm_bindgen_futures::spawn_local;
#[derive(Clone)]
struct MetadataState {
    context: Context,
    review: wire::MetadataReview,
    deleted_name: Option<String>,
}
#[derive(Clone)]
struct ChoiceState {
    context: Context,
    requirements: wire::PushChoices,
}
#[derive(Clone)]
struct PushState {
    context: Context,
    review: wire::PushPublicationReview,
}
#[derive(Clone)]
struct DeletedState {
    context: Context,
    page: wire::DeletedFolderPage,
    offset: u16,
}
#[derive(Clone, Copy)]
enum Action {
    Folder,
    Document,
    Inspect,
    Prepare,
    QueueMetadata,
    QueuePush,
    DeletedList { offset: u16, refresh: bool },
    DeletedReview { index: usize },
}
fn discard(client: Client, id: wire::ReviewId) {
    spawn_local(async move {
        let _ = client
            .request(wire::Command::DiscardPushReview { review: id }, None)
            .await;
    });
}
fn description(action: &wire::MetadataAction) -> String {
    fn parent(id: &Option<base::FolderId>) -> String {
        id.as_ref()
            .map(|id| format!("folder {}", id.as_str()))
            .unwrap_or_else(|| "the remote root".into())
    }
    match action {
        wire::MetadataAction::CreateFolder { parent: p, name } => {
            format!("Create remote folder “{}” in {}.", name.as_str(), parent(p))
        }
        wire::MetadataAction::RenameFolder {
            id,
            expected_generation,
            name,
        } => format!(
            "Rename remote folder {} to “{}”, only at generation {}.",
            id.as_str(),
            name.as_str(),
            expected_generation
        ),
        wire::MetadataAction::MoveFolder {
            id,
            expected_generation,
            parent: p,
        } => format!(
            "Move remote folder {} to {}, only at generation {}.",
            id.as_str(),
            parent(p),
            expected_generation
        ),
        wire::MetadataAction::RenameDocument {
            id,
            expected_generation,
            name,
        } => format!(
            "Rename remote plasmid {} to “{}”, only at generation {}.",
            id.as_str(),
            name.as_str(),
            expected_generation
        ),
        wire::MetadataAction::MoveDocument {
            id,
            expected_generation,
            folder,
        } => format!(
            "Move remote plasmid {} to {}, only at generation {}.",
            id.as_str(),
            parent(folder),
            expected_generation
        ),
        wire::MetadataAction::DeleteEmptyFolder {
            id,
            expected_generation,
        } => format!(
            "Delete empty remote folder {} only at generation {}.",
            id.as_str(),
            expected_generation
        ),
    }
}
#[component]
pub fn NativePublicationControls(
    client: Client,
    saved: Memo<Option<SavedLocalDocument>>,
    selected: RwSignal<Option<wire::Remote>>,
    destination: Memo<Option<base::FolderId>>,
    locked: Memo<bool>,
    busy: RwSignal<bool>,
    on_queued: Callback<()>,
) -> impl IntoView {
    let client = StoredValue::new_local(client);
    let life = StoredValue::new_local(Rc::new(Cell::new(0u64)));
    let active = StoredValue::new_local(Rc::new(Cell::new(false)));
    let metadata = RwSignal::new(None::<MetadataState>);
    let deleted = RwSignal::new(None::<DeletedState>);
    let alternatives = RwSignal::new(None::<ChoiceState>);
    let chosen = RwSignal::new(Chosen::new());
    let publication = RwSignal::new(None::<PushState>);
    let status = RwSignal::new(String::new());
    let context = Memo::new(move |_| {
        selected.get().map(|remote| Context {
            source: saved.get().map(|s| SavedPin {
                document: s.id,
                revision: s.revision,
                dirty: s.dirty,
                name: s.name,
            }),
            remote,
            destination: destination.get(),
        })
    });
    let observed = StoredValue::new(None::<Context>);
    Effect::new(move |_| {
        let now = context.get();
        let blocked = locked.get();
        if observed.get_value() != now || blocked {
            observed.set_value(now);
            life.with_value(|life| life.set(life.get().wrapping_add(1)));
            metadata.set(None);
            deleted.set(None);
            alternatives.set(None);
            chosen.set(Chosen::new());
            if let Some(old) = publication.get_untracked() {
                discard(client.get_value(), old.review.id);
                publication.set(None);
            }
        }
    });
    on_cleanup(move || {
        life.with_value(|life| life.set(life.get().wrapping_add(1)));
        if let Some(old) = publication.get_untracked() {
            discard(client.get_value(), old.review.id);
        }
        if active.with_value(|active| active.replace(false)) {
            let _ = busy.try_set(false);
        }
    });
    let send = Callback::new(move |action: Action| {
        if !crate::platform::desktop_rpc::is_desktop()
            || busy.get_untracked()
            || locked.get_untracked()
        {
            return;
        }
        let setup = (|| -> Result<(Context, wire::Command), String> {
            let pin = context
                .get_untracked()
                .ok_or("Choose a configured remote first.")?;
            let command=match action{
                Action::DeletedList{offset,refresh}=>{
                    let expected_checkpoint=if refresh{None}else{deleted.get_untracked().filter(|state|state.context==pin).ok_or("Refresh deleted folder choices first.")?.page.checkpoint};
                    if offset>0&&expected_checkpoint.is_none(){return Err("The deleted-folder snapshot changed. Refresh its first page.".into());}
                    wire::Command::ListDeletedFolderPublications{remote:pin.remote.id.clone(),page:wire::LocalPageRequest{offset,limit:25},expected_checkpoint}
                },
                Action::DeletedReview{index}=>{
                    let state=deleted.get_untracked().filter(|state|state.context==pin).ok_or("Refresh deleted folder choices first.")?;
                    let item=state.page.items.get(index).ok_or("Choose a known deleted local folder from this page.")?;
                    wire::Command::PrepareMetadataPublication{target:wire::PublicationTarget::DeletedFolder(item.local_folder.clone()),remote:pin.remote.id.clone()}
                },
                Action::Folder=>wire::Command::PrepareMetadataPublication{target:wire::PublicationTarget::Folder(pin.destination.clone().ok_or("Select a local folder in the local library first. The virtual root cannot be published as a folder.")?),remote:pin.remote.id.clone()},
                Action::Document=>wire::Command::PrepareMetadataPublication{target:wire::PublicationTarget::DocumentMetadata(pin.clean_source()?.document.clone()),remote:pin.remote.id.clone()},
                Action::Inspect=>{let source=pin.clean_source()?;wire::Command::InspectPushChoices{local_document:source.document.clone(),expected_local_revision:source.revision.clone(),remote:pin.remote.id.clone()}},
                Action::Prepare=>{
                    let source=pin.clean_source()?;let state=alternatives.get_untracked().filter(|state|state.context==pin).ok_or("Inspect the saved source choices again.")?;
                    wire::Command::PrepareReviewedPush{local_document:source.document.clone(),expected_local_revision:source.revision.clone(),remote:pin.remote.id.clone(),choices:model::choices(&state.requirements,&chosen.get_untracked())?}
                },
                Action::QueueMetadata=>{let state=metadata.get_untracked().filter(|state|state.context==pin).ok_or("Metadata source or destination changed; prepare a new review.")?;wire::Command::QueueMetadataPublication{review:state.review}},
                Action::QueuePush=>{let source=pin.clean_source()?;let state=publication.get_untracked().filter(|state|state.context==pin).ok_or("Saved source or remote changed; prepare a new review.")?;wire::Command::QueueReviewedPush{local_document:source.document.clone(),expected_local_revision:source.revision.clone(),remote:pin.remote.id.clone(),review:state.review.id}},
            };
            Ok((pin, command))
        })();
        let (pin, command) = match setup {
            Ok(value) => value,
            Err(error) => {
                status.set(error);
                return;
            }
        };
        if !matches!(action, Action::QueuePush) {
            if let Some(old) = publication.get_untracked() {
                discard(client.get_value(), old.review.id);
                publication.set(None);
            }
        }
        if !matches!(action, Action::QueueMetadata) {
            metadata.set(None);
        }
        let deleted_choice = if let Action::DeletedReview { index } = action {
            deleted
                .get_untracked()
                .and_then(|state| state.page.items.get(index).cloned())
        } else {
            None
        };
        let rpc = client.get_value();
        let life = life.get_value();
        let epoch = life.get().wrapping_add(1);
        life.set(epoch);
        let active = active.get_value();
        active.set(true);
        busy.set(true);
        status.set("Preparing native publication review…".into());
        spawn_local(async move {
            let result = rpc.request(command, Some(pin.remote.clone())).await;
            if active.replace(false) {
                let _ = busy.try_set(false);
            }
            let current = life.get() == epoch
                && context.try_get_untracked().flatten().as_ref() == Some(&pin)
                && !locked.try_get_untracked().unwrap_or(true)
                && status.try_get_untracked().is_some();
            if !current {
                match result {
                    Ok(wire::Success::PushPrepared(review)) => discard(rpc, review.id),
                    // A captured, explicitly approved job may have committed
                    // while selection changed. Refresh durable job discovery
                    // without writing a stale review into the new source panel.
                    Ok(wire::Success::Queued(_)) if status.try_get_untracked().is_some() => {
                        on_queued.run(())
                    }
                    _ => (),
                }
                return;
            }
            match result {
                Ok(wire::Success::MetadataPrepared(review))
                    if model::configuration_matches(&review.description.remote, &pin.remote) =>
                {
                    if deleted_choice.as_ref().is_some_and(|item| {
                        !model::deleted_review_matches(item, &review.description)
                    }) {
                        deleted.set(None);
                        status.set("The selected deleted-folder binding changed. Refresh its named choices and review again.".into());
                        return;
                    }
                    status.set("Review the exact remote metadata change, then confirm to queue one operation.".into());
                    metadata.set(Some(MetadataState {
                        context: pin,
                        review,
                        deleted_name: deleted_choice
                            .as_ref()
                            .map(|item| item.remote_folder.name.as_str().to_owned()),
                    }));
                }
                Ok(wire::Success::DeletedFolderPublications(page)) => {
                    let Action::DeletedList { offset, .. } = action else {
                        status.set("Unexpected deleted-folder listing.".into());
                        return;
                    };
                    deleted.set(Some(DeletedState {
                        context: pin,
                        page,
                        offset,
                    }));
                    status.set("Choose a known deleted local folder to review one remote empty-folder deletion. Nothing is queued yet.".into());
                }
                Ok(wire::Success::PushChoices(requirements)) => {
                    chosen.set(Chosen::new());
                    publication.set(None);
                    alternatives.set(Some(ChoiceState {
                        context: pin,
                        requirements,
                    }));
                    status.set("Choose only the reachable sequence alternatives you intend to publish, then prepare the subset review.".into());
                }
                Ok(wire::Success::PushPrepared(review)) => {
                    if !model::configuration_matches(&review.configuration, &pin.remote) {
                        discard(rpc, review.id);
                        status.set("The native review names a different remote configuration. Refresh remotes and review again.".into());
                        return;
                    }
                    publication.set(Some(PushState {
                        context: pin,
                        review,
                    }));
                    status.set("Review the original/subset hashes and omission counts before explicitly queueing this publication.".into());
                }
                Ok(wire::Success::Queued(_)) => {
                    metadata.set(None);
                    deleted.set(None);
                    publication.set(None);
                    alternatives.set(None);
                    chosen.set(Chosen::new());
                    status.set("Publication queued durably. Resume it from Pending pushes when ready. The working editor and private originals are unchanged.".into());
                    on_queued.run(());
                }
                Ok(_) => status.set(
                    "Unexpected native publication response; refresh and review again.".into(),
                ),
                Err(error) => {
                    if matches!(
                        action,
                        Action::DeletedList { .. } | Action::DeletedReview { .. }
                    ) {
                        deleted.set(None);
                    }
                    status.set(error);
                }
            }
        });
    });
    view! {<Show when=crate::platform::desktop_rpc::is_desktop><details class="native-publication-controls"><summary>"Publish folders, names and SBOL sources"</summary>
        <p class="helper">"Each action prepares a review for the selected remote. Confirming queues one durable job; it does not send it. Resume jobs explicitly from Pending pushes."</p>
        <div class="method-row-actions"><button class="text-tool" disabled=move||busy.get()||locked.get()||selected.with(Option::is_none)||destination.with(Option::is_none) on:click=move|_|send.run(Action::Folder)>"Review selected local folder"</button><button class="text-tool" disabled=move||busy.get()||locked.get()||selected.with(Option::is_none)||saved.with(|s|s.as_ref().is_none_or(|s|s.dirty)) on:click=move|_|send.run(Action::Document)>"Review saved plasmid name or folder"</button></div>
        {move||metadata.get().map(|state|{let d=state.review.description;view!{<section class="native-metadata-publication"><h4>"Confirm one remote metadata change"</h4>{state.deleted_name.map(|name|view!{<p>{format!("Delete the empty remote folder last synchronized as “{name}”.")}</p>})}<p>{format!("{} · {} · lab {}",d.remote.name.as_str(),d.remote.endpoint.as_str(),d.remote.lab.as_str())}</p><p>{description(&d.action)}</p><p class="helper">"Local folder placement does not imply a remote move for an unchanged mounted root. The native host will recheck the exact binding and generations before queueing. Additional changes need a separate review."</p><button class="button secondary" disabled=move||busy.get()||locked.get() on:click=move|_|send.run(Action::QueueMetadata)>"Confirm and queue this change"</button><button class="text-tool" disabled=move||busy.get() on:click=move|_|metadata.set(None)>"Cancel review"</button></section>}})}
        <details><summary>"Publish a local folder deletion"</summary><p class="helper">"Only known folder bindings deleted locally are listed. Choose by the last synchronized remote name, then review the exact empty-folder deletion. A nonempty or changed remote folder is protected by native checks and generation conflict handling."</p>
            <button class="text-tool" disabled=move||busy.get()||locked.get()||selected.with(Option::is_none) on:click=move|_|send.run(Action::DeletedList{offset:0,refresh:true})>"Refresh deleted folder choices"</button>
            {move||deleted.get().map(|state|{let previous=state.offset.saturating_sub(25);let next=state.page.next_offset;view!{<section>
                {state.page.items.is_empty().then(||view!{<p class="helper">"No known deleted folder bindings on this page."</p>})}
                <ul>{state.page.items.iter().enumerate().map(|(index,item)|view!{<li><button class="text-tool" disabled=move||busy.get()||locked.get() on:click=move|_|send.run(Action::DeletedReview{index})>{format!("Review deletion of “{}”",item.remote_folder.name.as_str())}</button><small>{format!("Remote folder {} · last synchronized generation {}",item.remote_folder.id.as_str(),item.remote_folder.generation)}</small></li>}).collect_view()}</ul>
                <div class="method-row-actions"><button class="text-tool" disabled=move||busy.get()||locked.get()||state.offset==0 on:click=move|_|send.run(Action::DeletedList{offset:previous,refresh:false})>"Previous choices"</button><button class="text-tool" disabled=move||busy.get()||locked.get()||next.is_none() on:click=move|_|if let Some(offset)=next{send.run(Action::DeletedList{offset,refresh:false});}>"Next choices"</button></div>
            </section>}})}
        </details>
        <details><summary>"Review SBOL documentary publication"</summary><p class="helper">"Use this when the saved native plasmid carries an SBOL source graph that ordinary Push cannot publish. The exact original stays in your private local library. Save dirty editor changes locally before preparing a review."</p>
            <button class="button secondary" disabled=move||busy.get()||locked.get()||context.with(|c|c.as_ref().is_none_or(|c|c.clean_source().is_err())) on:click=move|_|send.run(Action::Inspect)>"Inspect saved source choices"</button>
            {move||alternatives.get().map(|state|view!{<section><p>{format!("Saved library revision {}",state.requirements.local_revision.as_str())}</p>
                {state.requirements.requirements.is_empty().then(||view!{<p>"No reachable component requires an alternative choice. Prepare the documentary subset review to inspect its exact scope."</p>})}
                {state.requirements.requirements.iter().map(|requirement|{let key=(requirement.original_hash.as_str().to_owned(),requirement.component_iri.as_str().to_owned());view!{<fieldset disabled=move||busy.get()||locked.get()><legend>{requirement.component_iri.as_str().to_owned()}</legend><small>{format!("Original file hash {}",requirement.original_hash.as_str())}</small>{requirement.available_sequences.iter().map(|sequence|{let id=sequence.as_str().to_owned();let value=id.clone();let required=requirement.required_sequences.iter().any(|s|s==sequence);let check_key=key.clone();let update_key=key.clone();view!{<label class="checkbox-label"><input type="checkbox" disabled=required prop:checked=move||required||chosen.with(|all|all.get(&check_key).is_some_and(|all|all.contains(&value))) on:change=move|event|{let checked=event_target_checked(&event);chosen.update(|all|{let set=all.entry(update_key.clone()).or_default();if checked{set.insert(id.clone());}else{set.remove(&id);}});if let Some(old)=publication.get_untracked(){discard(client.get_value(),old.review.id);publication.set(None);}}/>{sequence.as_str().to_owned()}{required.then_some(" · required by retained native source")}</label>}}).collect_view()}</fieldset>}}).collect_view()}
                <button class="button secondary" disabled=move||busy.get()||locked.get() on:click=move|_|send.run(Action::Prepare)>"Prepare subset publication review"</button>
            </section>})}
            {move||publication.get().map(|state|view!{<section class="native-sbol-publication"><h4>"Approve this documentary subset"</h4><p>{format!("{} · {} · saved library revision {}",state.context.remote.name.as_str(),state.context.remote.endpoint.as_str(),state.review.local_revision.as_str())}</p>
                {state.review.entries.iter().map(|entry|view!{<div><p>"Original omitted: "<code>{entry.original_hash.as_str().to_owned()}</code></p><p>"Published subset: "<code>{entry.subset_hash.as_str().to_owned()}</code></p><p>{format!("{} retained objects · {} unrelated objects omitted · {} alternative sequence edges omitted",entry.retained_objects,entry.omitted_objects,entry.omitted_sequence_edges)}</p></div>}).collect_view()}
                <p class="helper">"Complete originals remain retained locally. The reviewed selected native history and structural source evidence will be queued; the entries above identify source-file omissions and their receipts. Retained annotations and descriptions remain included. Native source/event identities and physical history are unchanged; this documentary projection adds no scientific authority. The native host rechecks the saved head, source bytes and remote binding."</p>
                <button class="button secondary" disabled=move||busy.get()||locked.get() on:click=move|_|send.run(Action::QueuePush)>"Approve subset and queue push"</button><button class="text-tool" disabled=move||busy.get() on:click=move|_|{if let Some(old)=publication.get_untracked(){discard(client.get_value(),old.review.id);}publication.set(None);}>"Cancel review"</button>
            </section>})}
        </details><p role="status" aria-live="polite">{move||status.get()}</p>
    </details></Show>}
}
