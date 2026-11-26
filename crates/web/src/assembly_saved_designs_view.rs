//! Compact named editable declarations. Host callbacks install attachments or
//! replace the working declaration only; this view never restores a verdict.
use crate::{assembly_draft_model::Draft, assembly_saved_designs_model as model};
use collagene_core::workspace::Workspace;
use leptos::prelude::*;
use std::sync::Arc;

const PAGE_SIZE: usize = 20;

#[derive(Clone, Copy, PartialEq, Eq)]
enum Action {
    Load,
    Delete,
}
#[derive(Clone)]
struct Pending {
    serial: u64,
    action: Action,
    entry: model::Entry,
    epoch: u64,
    before: Arc<Draft>,
    /// A load preview attests exact source events before asking for replacement.
    /// Delete deliberately has no hydration requirement, allowing orphan cleanup.
    loaded: Option<Arc<Draft>>,
}
#[derive(Clone)]
struct Rename {
    serial: u64,
    entry: model::Entry,
    epoch: u64,
}
#[derive(Clone)]
struct Recovery {
    serial: u64,
    entry: model::MalformedEntry,
    epoch: u64,
    before: Arc<Draft>,
}

fn current_entry(frame: &Workspace, expected: &model::Entry) -> Result<model::Entry, String> {
    let entry = model::list(frame)?
        .into_iter()
        .find(|entry| entry.id == expected.id)
        .ok_or("The selected saved design is unavailable. Refresh the list.")?;
    if entry.attachment_hash != expected.attachment_hash {
        return Err(
            "The selected saved design changed. Refresh the list and choose it again.".into(),
        );
    }
    Ok(entry)
}

#[component]
pub(crate) fn SavedDesigns(
    draft: RwSignal<Draft>,
    epoch: RwSignal<u64>,
    busy: Memo<bool>,
    snapshot: Callback<(), Result<Workspace, String>>,
    apply_attachments: Callback<(Workspace, Workspace), Result<(), String>>,
    replace_draft: Callback<(Draft, Draft), Result<(), String>>,
) -> impl IntoView {
    let expanded = RwSignal::new(false);
    let action_serial = RwSignal::new(0u64);
    let next_serial = move || {
        action_serial.update(|serial| *serial = serial.wrapping_add(1));
        action_serial.get_untracked()
    };
    let refresh = RwSignal::new(0u64);
    let entries = RwSignal::new(Vec::<model::Entry>::new());
    let malformed = RwSignal::new(Vec::<model::MalformedEntry>::new());
    let recovery = RwSignal::new(None::<Recovery>);
    let recovery_page = RwSignal::new(0usize);
    let page = RwSignal::new(0usize);
    let name = RwSignal::new(String::new());
    let rename_name = RwSignal::new(String::new());
    let renaming = RwSignal::new(None::<Rename>);
    let pending = RwSignal::new(None::<Pending>);
    let message = RwSignal::new(String::new());
    let list_error = RwSignal::new(String::new());
    let refresh_list =
        Callback::new(move |()| refresh.update(|value| *value = value.wrapping_add(1)));

    Effect::new(move |_| {
        let current_epoch = epoch.get();
        let locked = busy.get();
        let open = expanded.get();
        draft.with(|current| {
            if pending.with_untracked(|saved| {
                saved.as_ref().is_some_and(|saved| {
                    locked
                        || !open
                        || saved.epoch != current_epoch
                        || current != saved.before.as_ref()
                })
            }) {
                pending.set(None);
                message.set(
                    "The draft or workspace changed. Choose the saved-design action again.".into(),
                );
            }
        });
        if recovery.with_untracked(|saved| {
            saved.as_ref().is_some_and(|saved| {
                locked
                    || !open
                    || saved.epoch != current_epoch
                    || !draft.with_untracked(|current| current == saved.before.as_ref())
            })
        }) {
            recovery.set(None);
        }
        if renaming.with_untracked(|saved| {
            saved
                .as_ref()
                .is_some_and(|saved| locked || !open || saved.epoch != current_epoch)
        }) {
            renaming.set(None);
        }
    });
    Effect::new(move |_| {
        let _ = (epoch.get(), refresh.get());
        if !expanded.get() || busy.get() {
            return;
        }
        match snapshot.run(()).and_then(|frame| {
            let listed = model::list(&frame);
            let broken = if listed.is_err() {
                model::malformed(&frame)?
            } else {
                Vec::new()
            };
            Ok((listed, broken))
        }) {
            Ok((Ok(current), _)) => {
                malformed.set(Vec::new());
                recovery.set(None);
                let unchanged = |entry: &model::Entry| {
                    current.iter().any(|now| {
                        now.id == entry.id && now.attachment_hash == entry.attachment_hash
                    })
                };
                if pending.with_untracked(|saved| {
                    saved.as_ref().is_some_and(|saved| !unchanged(&saved.entry))
                }) {
                    pending.set(None);
                    message
                        .set("The selected saved design changed. Choose the action again.".into());
                }
                if renaming.with_untracked(|saved| {
                    saved.as_ref().is_some_and(|saved| !unchanged(&saved.entry))
                }) {
                    renaming.set(None);
                }
                page.update(|page| {
                    *page = (*page).min(current.len().saturating_sub(1) / PAGE_SIZE)
                });
                entries.set(current);
                list_error.set(String::new());
            }
            Ok((Err(error), broken)) => {
                entries.set(Vec::new());
                pending.set(None);
                renaming.set(None);
                if recovery.with_untracked(|saved| {
                    saved
                        .as_ref()
                        .is_some_and(|saved| !broken.iter().any(|entry| entry == &saved.entry))
                }) {
                    recovery.set(None);
                }
                recovery_page
                    .update(|page| *page = (*page).min(broken.len().saturating_sub(1) / PAGE_SIZE));
                malformed.set(broken);
                list_error.set(error);
            }
            Err(error) => {
                entries.set(Vec::new());
                malformed.set(Vec::new());
                pending.set(None);
                renaming.set(None);
                recovery.set(None);
                list_error.set(error);
            }
        }
    });

    let save = move |_| {
        let outcome = (|| {
            if busy.get_untracked() {
                return Err("Wait for the workspace update to finish.".into());
            }
            let saved_epoch = epoch.get_untracked();
            let declared = draft.get_untracked();
            let before = snapshot.run(())?;
            let (after, entry) = model::save(&before, &name.get_untracked(), &declared)?;
            if busy.get_untracked()
                || epoch.get_untracked() != saved_epoch
                || !draft.with_untracked(|current| current == &declared)
            {
                return Err("The draft or workspace changed before saving the design.".into());
            }
            apply_attachments.run((before, after))?;
            Ok::<_, String>(entry)
        })();
        match outcome {
            Ok(entry) => {
                pending.set(None);
                renaming.set(None);
                recovery.set(None);
                name.set(String::new());
                page.set(0);
                message.set(format!(
                    "Saved ‘{}’. Load it later to review intent and recompute preparation.",
                    entry.name
                ));
                refresh_list.run(());
            }
            Err(error) => message.set(error),
        }
    };
    let propose = Callback::new(move |(entry, action): (model::Entry, Action)| {
        let outcome = (|| {
            if busy.get_untracked() {
                return Err("Wait for the workspace update to finish.".into());
            }
            let saved_epoch = epoch.get_untracked();
            let before = Arc::new(draft.get_untracked());
            let frame = snapshot.run(())?;
            let entry = current_entry(&frame, &entry)?;
            let loaded = match action {
                Action::Load => Some(Arc::new(model::load(&frame, &entry.id)?)),
                Action::Delete => None,
            };
            if busy.get_untracked()
                || epoch.get_untracked() != saved_epoch
                || !draft.with_untracked(|current| current == before.as_ref())
            {
                return Err("The draft or workspace changed. Choose the action again.".into());
            }
            Ok::<_, String>(Pending {
                serial: next_serial(),
                action,
                entry,
                epoch: saved_epoch,
                before,
                loaded,
            })
        })();
        match outcome {
            Ok(value) => {
                renaming.set(None);
                recovery.set(None);
                pending.set(Some(value));
                message.set(String::new());
            }
            Err(error) => {
                pending.set(None);
                message.set(error);
                refresh_list.run(());
            }
        }
    });
    let confirm = Callback::new(move |saved: Pending| {
        // An old DOM handler cannot apply or dismiss a newer confirmation.
        if pending.with_untracked(|now| {
            now.as_ref().is_none_or(|now| {
                now.serial != saved.serial
                    || now.action != saved.action
                    || now.entry.id != saved.entry.id
                    || now.entry.attachment_hash != saved.entry.attachment_hash
                    || now.epoch != saved.epoch
                    || now.before.as_ref() != saved.before.as_ref()
            })
        }) {
            message.set(
                "This confirmation is no longer current. Choose the saved-design action again."
                    .into(),
            );
            return;
        }
        let outcome = (|| {
            if busy.get_untracked()
                || epoch.get_untracked() != saved.epoch
                || !draft.with_untracked(|current| current == saved.before.as_ref())
            {
                return Err("The draft or workspace changed. Review this action again.".into());
            }
            let before = snapshot.run(())?;
            current_entry(&before, &saved.entry)?;
            match saved.action {
                Action::Load => {
                    let loaded = model::load(&before, &saved.entry.id)?;
                    if saved
                        .loaded
                        .as_ref()
                        .is_none_or(|preview| preview.as_ref() != &loaded)
                    {
                        return Err("The saved source events or declaration changed. Choose the design again.".into());
                    }
                    if busy.get_untracked()
                        || epoch.get_untracked() != saved.epoch
                        || !draft.with_untracked(|current| current == saved.before.as_ref())
                    {
                        return Err("The draft or workspace changed before loading.".into());
                    }
                    replace_draft.run((saved.before.as_ref().clone(), loaded))?;
                    Ok(format!(
                        "Loaded ‘{}’. Review its intended construct and recompute preparation.",
                        saved.entry.name
                    ))
                }
                Action::Delete => {
                    let after = model::delete(&before, &saved.entry.id)?;
                    if busy.get_untracked()
                        || epoch.get_untracked() != saved.epoch
                        || !draft.with_untracked(|current| current == saved.before.as_ref())
                    {
                        return Err("The draft or workspace changed before deleting.".into());
                    }
                    apply_attachments.run((before, after))?;
                    Ok(format!("Deleted saved design ‘{}’.", saved.entry.name))
                }
            }
        })();
        pending.set(None);
        match outcome {
            Ok(text) => {
                message.set(text);
                refresh_list.run(());
            }
            Err(error) => {
                message.set(error);
                refresh_list.run(());
            }
        }
    });
    let begin_rename = Callback::new(move |entry: model::Entry| {
        if busy.get_untracked() {
            return;
        }
        pending.set(None);
        recovery.set(None);
        rename_name.set(entry.name.clone());
        renaming.set(Some(Rename {
            serial: next_serial(),
            entry,
            epoch: epoch.get_untracked(),
        }));
        message.set(String::new());
    });
    let rename = Callback::new(move |saved: Rename| {
        if renaming.with_untracked(|now| now.as_ref().is_none_or(|now| now.serial != saved.serial))
        {
            return;
        }
        let outcome = (|| {
            if busy.get_untracked() || epoch.get_untracked() != saved.epoch {
                return Err("The workspace changed. Choose the saved design again.".into());
            }
            let before = snapshot.run(())?;
            current_entry(&before, &saved.entry)?;
            let after = model::rename(&before, &saved.entry.id, &rename_name.get_untracked())?;
            if busy.get_untracked() || epoch.get_untracked() != saved.epoch {
                return Err("The workspace changed before renaming.".into());
            }
            apply_attachments.run((before, after))?;
            Ok::<_, String>(())
        })();
        match outcome {
            Ok(()) => {
                renaming.set(None);
                message.set("Saved design renamed.".into());
                refresh_list.run(());
            }
            Err(error) => {
                message.set(error);
                refresh_list.run(());
            }
        }
    });

    let propose_recovery = Callback::new(move |entry: model::MalformedEntry| {
        let outcome = (|| {
            if busy.get_untracked() {
                return Err("Wait for the workspace update to finish.".to_owned());
            }
            let saved_epoch = epoch.get_untracked();
            let before = Arc::new(draft.get_untracked());
            let frame = snapshot.run(())?;
            if !model::malformed(&frame)?.iter().any(|now| now == &entry) {
                return Err("The unreadable saved attachment changed. Refresh the list.".into());
            }
            if busy.get_untracked()
                || epoch.get_untracked() != saved_epoch
                || !draft.with_untracked(|current| current == before.as_ref())
            {
                return Err("The draft or workspace changed. Choose the attachment again.".into());
            }
            Ok::<_, String>(Recovery {
                serial: next_serial(),
                entry,
                epoch: saved_epoch,
                before,
            })
        })();
        match outcome {
            Ok(saved) => {
                pending.set(None);
                renaming.set(None);
                recovery.set(Some(saved));
                message.set(String::new());
            }
            Err(error) => {
                message.set(error);
                refresh_list.run(());
            }
        }
    });
    let discard = Callback::new(move |saved: Recovery| {
        if recovery.with_untracked(|now| now.as_ref().is_none_or(|now| now.serial != saved.serial))
        {
            return;
        }
        let outcome = (|| {
            if busy.get_untracked()
                || epoch.get_untracked() != saved.epoch
                || !draft.with_untracked(|current| current == saved.before.as_ref())
            {
                return Err("The draft or workspace changed. Review removal again.".to_owned());
            }
            let before = snapshot.run(())?;
            if !model::malformed(&before)?
                .iter()
                .any(|now| now == &saved.entry)
            {
                return Err("The unreadable saved attachment changed. Refresh the list.".into());
            }
            let after = model::discard_malformed(&before, &saved.entry.attachment_hash)?;
            if busy.get_untracked()
                || epoch.get_untracked() != saved.epoch
                || !draft.with_untracked(|current| current == saved.before.as_ref())
            {
                return Err("The workspace changed before removing the attachment.".into());
            }
            apply_attachments.run((before, after))?;
            Ok::<_, String>(())
        })();
        recovery.set(None);
        match outcome {
            Ok(()) => message.set("Unreadable named-design attachment removed.".into()),
            Err(error) => message.set(error),
        }
        refresh_list.run(());
    });

    view! {<details class="assembly-saved-designs" on:toggle=move|event:web_sys::Event|
        expanded.set(event_target::<web_sys::HtmlElement>(&event).has_attribute("open"))>
        <summary>"Saved designs"</summary>
        <div class="assembly-draft-menu-content">
            <p class="helper">"Keep named editable snapshots in this workspace. Loading replaces the working draft and retains one temporary previous draft for this session. Save a named snapshot first to retain your current work after restarting; review intent and preparation again after loading."</p>
            <label>"Saved design name"<input maxlength="160" aria-label="Saved assembly design name" disabled=move||busy.get() prop:value=move||name.get() on:input=move|event|name.set(event_target_value(&event))/></label>
            <div class="method-row-actions"><button class="button secondary" disabled=move||busy.get()||!list_error.get().is_empty()||entries.with(|items|items.len()>=model::MAX_DESIGNS) on:click=save>"Save current draft"</button><button class="text-tool" disabled=move||busy.get() on:click=move|_|refresh_list.run(())>"Refresh list"</button></div>
            <p role="status">{move||message.get()}</p>
            <Show when=move||!list_error.get().is_empty()><p class="dialog-error" role="alert">{move||list_error.get()}</p></Show>
            <Show when=move||!list_error.get().is_empty()&&!malformed.with(Vec::is_empty)>
                <section aria-label="Unreadable saved designs"><h4>"Unreadable saved designs"</h4>
                    <p class="helper">"These named-design attachments cannot be read. Remove an unwanted attachment to recover the list; original sequences and history are retained."</p>
                    <ul>{move||malformed.with(|items|items.iter().skip(recovery_page.get().saturating_mul(PAGE_SIZE)).take(PAGE_SIZE).cloned().collect::<Vec<_>>()).into_iter().map(|entry|{
                        let remove=entry.clone();
                        view!{<li><strong>{entry.name}</strong><p>{entry.error}</p><button class="text-tool" disabled=move||busy.get() on:click=move|_|propose_recovery.run(remove.clone())>"Remove unreadable attachment…"</button></li>}
                    }).collect_view()}</ul>
                    <Show when=move||malformed.with(|items|items.len()>PAGE_SIZE)><div class="method-row-actions">
                        <button class="text-tool" disabled=move||recovery_page.get()==0 on:click=move|_|recovery_page.update(|page|*page=page.saturating_sub(1))>"Previous unreadable designs"</button>
                        <button class="text-tool" disabled=move||malformed.with(|items|recovery_page.get().saturating_add(1).saturating_mul(PAGE_SIZE)>=items.len()) on:click=move|_|recovery_page.update(|page|*page=page.saturating_add(1))>"Next unreadable designs"</button>
                    </div></Show>
                </section>
            </Show>
            <Show when=move||list_error.get().is_empty()><p class="helper">{move||entries.with(|items|format!("{} of {} saved designs",items.len(),model::MAX_DESIGNS))}</p></Show>
            <Show when=move||entries.with(Vec::is_empty)&&list_error.get().is_empty()><p>"No named designs saved yet."</p></Show>
            <ul class="assembly-saved-design-list" style="max-height:320px;overflow:auto;list-style:none;padding:0">
                {move||entries.with(|items|items.iter().skip(page.get().saturating_mul(PAGE_SIZE)).take(PAGE_SIZE).cloned().collect::<Vec<_>>()).into_iter().map(|entry|{
                    let load=entry.clone();let delete=entry.clone();let rename_entry=entry.clone();
                    view!{<li><strong>{entry.name}</strong><span class="helper">{format!(" · {} source references",entry.source_count)}</span>
                        <div class="method-row-actions"><button class="text-tool" disabled=move||busy.get() on:click=move|_|propose.run((load.clone(),Action::Load))>"Load…"</button><button class="text-tool" disabled=move||busy.get() on:click=move|_|begin_rename.run(rename_entry.clone())>"Rename"</button><button class="text-tool" disabled=move||busy.get() on:click=move|_|propose.run((delete.clone(),Action::Delete))>"Delete…"</button></div>
                    </li>}
                }).collect_view()}
            </ul>
            <Show when=move||entries.with(|items|items.len()>PAGE_SIZE)><div class="method-row-actions">
                <button class="text-tool" disabled=move||page.get()==0 on:click=move|_|page.update(|page|*page=page.saturating_sub(1))>"Previous designs"</button>
                <span>{move||entries.with(|items|format!("Page {} of {}",page.get()+1,items.len().max(1).div_ceil(PAGE_SIZE)))}</span>
                <button class="text-tool" disabled=move||entries.with(|items|page.get().saturating_add(1).saturating_mul(PAGE_SIZE)>=items.len()) on:click=move|_|page.update(|page|*page=page.saturating_add(1))>"Next designs"</button>
            </div></Show>
            {move||recovery.get().map(|saved|{
                let cancel_serial=saved.serial;
                view!{<section class="assembly-draft-confirm" aria-label="Remove unreadable saved design"><h4>"Remove unreadable saved attachment?"</h4>
                    <p>{saved.entry.name.clone()}</p><p>{saved.entry.error.clone()}</p>
                    <p class="helper">"Only this unreadable named-design attachment will be removed. Source documents, history and the working draft are retained."</p>
                    <div class="method-row-actions"><button class="button primary" disabled=move||busy.get() on:click=move|_|discard.run(saved.clone())>"Remove attachment"</button><button class="text-tool" on:click=move|_|{
                        if recovery.with_untracked(|now|now.as_ref().is_some_and(|now|now.serial==cancel_serial)){recovery.set(None);}
                    }>"Cancel"</button></div>
                </section>}
            })}
            {move||renaming.get().map(|saved|{
                let cancel_serial=saved.serial;
                view!{<section class="assembly-draft-confirm" aria-label="Rename saved design"><h4>{format!("Rename ‘{}’",saved.entry.name)}</h4>
                    <label>"New name"<input maxlength="160" aria-label="New saved assembly design name" disabled=move||busy.get() prop:value=move||rename_name.get() on:input=move|event|rename_name.set(event_target_value(&event))/></label>
                    <div class="method-row-actions"><button class="button primary" disabled=move||busy.get() on:click=move|_|rename.run(saved.clone())>"Rename saved design"</button><button class="text-tool" on:click=move|_|{
                        if renaming.with_untracked(|now|now.as_ref().is_some_and(|now|now.serial==cancel_serial)){renaming.set(None);}
                    }>"Cancel"</button></div>
                </section>}
            })}
            {move||pending.get().map(|saved|{
                let loading=matches!(saved.action,Action::Load);
                let cancel_serial=saved.serial;
                let title=if loading{"Load saved assembly design?"}else{"Delete saved assembly design?"};
                let current_name=saved.before.name.clone();
                let selected_name=saved.entry.name.clone();
                view!{<section class="assembly-draft-confirm" aria-label=title><h4>{title}</h4><p>{format!("Saved design: {selected_name}")}</p>
                    {loading.then(||view!{<p>{format!("Replace current draft: {current_name}")}</p><p class="helper">"The current draft replaces the temporary previous draft for this session. Save a named snapshot first to retain it after restarting. Review intent and preparation after loading."</p>})}
                    {(!loading).then(||view!{<p class="helper">"Only this named snapshot is removed from the workspace."</p>})}
                    <div class="method-row-actions"><button class="button primary" disabled=move||busy.get() on:click=move|_|confirm.run(saved.clone())>{if loading{"Load design"}else{"Delete saved design"}}</button><button class="text-tool" on:click=move|_|{if pending.with_untracked(|now|now.as_ref().is_some_and(|now|now.serial==cancel_serial)){pending.set(None);}}>"Cancel"</button></div>
                </section>}
            })}
        </div>
    </details>}
}
