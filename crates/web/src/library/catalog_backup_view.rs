//! Native-only saved-catalog backup controls. Commands carry known IDs and pins,
//! never a filesystem path or biological backup archive.
use crate::library::catalog_backup_model::AdoptionReview;
use crate::library::catalog_backup_model::CancellationPin;
use crate::library::catalog_backup_model::RecoveryChoice;
use crate::library::catalog_backup_model::RestoreReview;
use crate::library::catalog_backup_model::VerificationPin;
use collagene_desktop_protocol::{self as wire, library};
use leptos::prelude::*;
use std::{cell::Cell, future::Future, pin::Pin, rc::Rc};
use wasm_bindgen::JsValue;
use wasm_bindgen_futures::spawn_local;
pub type LocalResult<T> = Pin<Box<dyn Future<Output = Result<T, String>>>>;
#[derive(Clone)]
struct Listing {
    page: wire::CatalogBackupPage,
    offset: u32,
}
#[derive(Clone)]
struct RecoveryListing {
    page: wire::CatalogRecoveryPage,
    offset: u32,
}
#[derive(Clone)]
enum Action {
    Create,
    List { offset: u32, refresh: bool },
    ListRecoveries { offset: u32, refresh: bool },
    Verify(wire::CatalogBackupEntry),
    Restore,
    ReviewRecovery(wire::CatalogRecoveryReceipt),
    AdoptRecovery,
    ReviewOriginal,
}
fn date(seconds: u64) -> String {
    let date = js_sys::Date::new(&JsValue::from_f64(seconds as f64 * 1000.));
    if date.get_time().is_finite() {
        date.to_iso_string()
            .as_string()
            .unwrap_or_else(|| "Date unavailable".into())
    } else {
        "Date unavailable".into()
    }
}
fn bytes(size: u64) -> String {
    if size >= 1 << 30 {
        format!("{:.2} GiB", size as f64 / (1u64 << 30) as f64)
    } else if size >= 1 << 20 {
        format!("{:.2} MiB", size as f64 / (1u64 << 20) as f64)
    } else if size >= 1 << 10 {
        format!("{:.1} KiB", size as f64 / (1u64 << 10) as f64)
    } else {
        format!("{size} bytes")
    }
}
fn receipt_details(receipt: &wire::CatalogBackupReceipt) -> AnyView {
    view!{<details><summary>"Backup receipt and hashes"</summary><dl>
        <dt>"Managed backup"</dt><dd>{receipt.basename.as_str().to_owned()}</dd>
        <dt>"Created (UTC)"</dt><dd>{date(receipt.created_unix_seconds)}</dd>
        <dt>"Archive / database size"</dt><dd>{format!("{} / {}",bytes(receipt.archive_bytes),bytes(receipt.database_bytes))}</dd>
        <dt>"Saved records"</dt><dd>{format!("{} plasmids · {} folders · {} revisions · {} durable mutation receipts",receipt.documents,receipt.folders,receipt.revisions,receipt.mutation_receipts)}</dd>
        <dt>"Archive hash"</dt><dd><code>{receipt.archive_hash.as_str().to_owned()}</code></dd>
        <dt>"Database hash"</dt><dd><code>{receipt.database_hash.as_str().to_owned()}</code></dd>
        <dt>"Schema"</dt><dd>{format!("Version {}",receipt.schema_version)}<code>{receipt.schema_hash.as_str().to_owned()}</code></dd>
        <dt>"Catalog lab"</dt><dd><code>{receipt.lab.as_str().to_owned()}</code></dd>
    </dl></details>}.into_any()
}
fn location_name(location: &wire::CatalogLocation) -> String {
    match location {
        wire::CatalogLocation::Original => "Original saved library".into(),
        wire::CatalogLocation::Recovery { id } => format!("Managed recovery {}", id.as_str()),
    }
}
fn snapshot_details(snapshot: &wire::CatalogSnapshot) -> AnyView {
    view!{<dl><dt>"Checked target contents"</dt><dd>{format!("{} plasmids · {} folders · {} revisions · {} durable mutation receipts",snapshot.documents,snapshot.folders,snapshot.revisions,snapshot.mutation_receipts)}</dd><dt>"Catalog generation"</dt><dd>{snapshot.catalog_generation}</dd><dt>"Checked database size"</dt><dd>{bytes(snapshot.database_bytes)}</dd><dt>"Checked database hash"</dt><dd><code>{snapshot.database_hash.as_str().to_string()}</code></dd><dt>"Schema hash"</dt><dd><code>{snapshot.schema_hash.as_str().to_string()}</code></dd></dl>}.into_any()
}
#[component]
pub fn CatalogBackupControls(
    client: Callback<wire::CatalogBackupCommand, LocalResult<wire::CatalogBackupSuccess>>,
    context: Memo<Option<library::LabContext>>,
    locked: Memo<bool>,
    on_catalog_pending: Callback<bool>,
) -> impl IntoView {
    let life = StoredValue::new_local(Rc::new(Cell::new(0u64)));
    let observed = StoredValue::new(None::<library::LabContext>);
    let busy = RwSignal::new(false);
    let status = RwSignal::new(String::new());
    let listing = RwSignal::new(None::<Listing>);
    let recoveries = RwSignal::new(None::<RecoveryListing>);
    let latest = RwSignal::new(None::<wire::CatalogBackupEntry>);
    let review = RwSignal::new(None::<RestoreReview>);
    let restored = RwSignal::new(None::<wire::CatalogRecoveryReceipt>);
    let adoption_review = RwSignal::new(None::<AdoptionReview>);
    let activation = RwSignal::new(None::<wire::CatalogActivationStatus>);
    let activation_checked = RwSignal::new(false);
    let catalog_status = RwSignal::new(None::<wire::CatalogStatus>);
    let current_location = Memo::new(move |_| {
        activation.with(|state| match state {
            Some(wire::CatalogActivationStatus::ActiveAfterRestart(receipt)) => {
                receipt.review.target.clone()
            }
            Some(wire::CatalogActivationStatus::PendingRestart(receipt)) => {
                receipt.review.current_catalog.clone()
            }
            Some(wire::CatalogActivationStatus::Cancelled(receipt)) => {
                receipt.accepted.review.current_catalog.clone()
            }
            _ => wire::CatalogLocation::Original,
        })
    });
    let catalog_paused = Memo::new(move |_| {
        if !activation_checked.get()
            || !matches!(
                catalog_status.get(),
                Some(wire::CatalogStatus::Ready { .. })
            )
        {
            return true;
        }
        activation.with(|state| {
            matches!(
                state,
                Some(wire::CatalogActivationStatus::PendingRestart(_))
            )
        })
    });
    Effect::new(move |_| {
        let now = context.get();
        let blocked = locked.get();
        if observed.get_value() != now || blocked {
            observed.set_value(now);
            life.with_value(|life| life.set(life.get().wrapping_add(1)));
            listing.set(None);
            recoveries.set(None);
            latest.set(None);
            review.set(None);
            restored.set(None);
            adoption_review.set(None);
            activation_checked.set(false);
        }
    });
    on_cleanup(move || life.with_value(|life| life.set(life.get().wrapping_add(1))));
    let run = Callback::new(move |action: Action| {
        if !crate::platform::desktop_rpc::is_desktop()
            || busy.get_untracked()
            || locked.get_untracked()
            || activation.with_untracked(|state| {
                matches!(
                    state,
                    Some(wire::CatalogActivationStatus::PendingRestart(_))
                )
            })
        {
            return;
        }
        let setup = (|| -> Result<_, String> {
            let local = context
                .get_untracked()
                .ok_or("Open the local saved library before managing catalog backups.")?;
            let verification = if let Action::Verify(entry) = &action {
                Some(VerificationPin::new(entry, &local)?)
            } else {
                None
            };
            let restoring = if matches!(&action, Action::Restore) {
                Some(
                    review
                        .get_untracked()
                        .ok_or("Verify a known backup before reviewing a restore copy.")?,
                )
            } else {
                None
            };
            let recovery_choice = if let Action::ReviewRecovery(receipt) = &action {
                Some(RecoveryChoice::new(receipt.clone(), &local)?)
            } else {
                None
            };
            let adopting = if matches!(&action, Action::AdoptRecovery) {
                Some(
                    adoption_review
                        .get_untracked()
                        .ok_or("Review the managed recovery copy before choosing it.")?,
                )
            } else {
                None
            };
            let command = match &action {
                Action::ReviewOriginal => wire::CatalogBackupCommand::ReviewOriginal,
                Action::ReviewRecovery(_) => recovery_choice
                    .as_ref()
                    .expect("prepared recovery choice")
                    .command(&local)?,
                Action::AdoptRecovery => adopting
                    .as_ref()
                    .expect("prepared adoption review")
                    .command(&local)?,
                Action::Create => wire::CatalogBackupCommand::Create,
                Action::ListRecoveries { offset, refresh } => {
                    let expected_catalog = if *refresh {
                        None
                    } else {
                        Some(
                            recoveries
                                .get_untracked()
                                .ok_or("Refresh known recovery copies first.")?
                                .page
                                .catalog_revision,
                        )
                    };
                    if *offset > 0 && expected_catalog.is_none() {
                        return Err("Refresh the first recovery page before continuing.".into());
                    }
                    wire::CatalogBackupCommand::ListRecoveries {
                        offset: *offset,
                        limit: 25,
                        expected_catalog,
                    }
                }
                Action::List { offset, refresh } => {
                    let expected_catalog = if *refresh {
                        None
                    } else {
                        Some(
                            listing
                                .get_untracked()
                                .ok_or("Refresh the known backup list first.")?
                                .page
                                .catalog_revision,
                        )
                    };
                    if *offset > 0 && expected_catalog.is_none() {
                        return Err("Refresh the first backup page before continuing.".into());
                    }
                    wire::CatalogBackupCommand::List {
                        offset: *offset,
                        limit: 25,
                        expected_catalog,
                    }
                }
                Action::Verify(_) => verification
                    .as_ref()
                    .expect("Prepared verification pin")
                    .command(),
                Action::Restore => restoring
                    .as_ref()
                    .expect("Prepared restore review")
                    .command(&local)?,
            };
            Ok((
                local,
                command,
                verification,
                restoring,
                recovery_choice,
                adopting,
            ))
        })();
        let (local, command, verification, restoring, recovery_choice, adopting) = match setup {
            Ok(value) => value,
            Err(error) => {
                status.set(error);
                return;
            }
        };
        // Any new verification/list/create requires a new explicit restore
        // confirmation. A list receipt alone never enables restoration.
        review.set(None);
        restored.set(None);
        adoption_review.set(None);
        let life = life.get_value();
        let epoch = life.get().wrapping_add(1);
        life.set(epoch);
        busy.set(true);
        status.set(
            match &action {
                Action::ReviewOriginal => {
                    "Checking the preserved original saved library for rollback…"
                }
                Action::ReviewRecovery(_) => {
                    "Checking the managed recovery copy and current catalog before adoption…"
                }
                Action::AdoptRecovery => {
                    "Selecting the checked recovery catalog for the next desktop start…"
                }
                Action::Create => "Making a streaming backup of the saved local catalog…",
                Action::List { .. } => "Loading known managed backups…",
                Action::ListRecoveries { .. } => "Loading known managed recovery copies…",
                Action::Verify(_) => "Verifying the selected saved backup…",
                Action::Restore => {
                    "Restoring to a new managed copy; current catalog remains in use…"
                }
            }
            .into(),
        );
        spawn_local(async move {
            let result = client.run(command).await;
            if busy.try_get_untracked().is_none() {
                return;
            }
            busy.set(false);
            if life.get() != epoch
                || context.try_get_untracked().flatten().as_ref() != Some(&local)
                || locked.try_get_untracked().unwrap_or(true)
            {
                return;
            }
            let handled = (|| -> Result<(), String> {
                match(result,&action){
                (Ok(wire::CatalogBackupSuccess::Created(receipt)),Action::Create)=>{
                    let entry=wire::CatalogBackupEntry{id:receipt.id.clone(),basename:receipt.basename.clone(),receipt:Some(receipt)};
                    VerificationPin::new(&entry,&local)?;listing.set(None);latest.set(Some(entry));
                    status.set("Saved catalog backup created. Verify it before reviewing a restore copy. Current autosaved editor files are outside this backup.".into());
                },
                (Ok(wire::CatalogBackupSuccess::Listed(page)),Action::List{offset,..})=>{
                    for entry in page.entries.iter(){VerificationPin::new(entry,&local)?;}
                    listing.set(Some(Listing{page,offset:*offset}));status.set("Choose a known backup by name and creation date. Unfinished entries can be inspected with Verify.".into());
                },
                (Ok(wire::CatalogBackupSuccess::Verified(receipt)),Action::Verify(_))=>{
                    let verified=verification.as_ref().expect("Prepared verification pin").accept(receipt.clone(),&local)?;
                    latest.set(Some(wire::CatalogBackupEntry{id:receipt.id.clone(),basename:receipt.basename.clone(),receipt:Some(receipt)}));listing.set(None);review.set(Some(verified));
                    status.set("Backup integrity verified. Review its scope and explicitly confirm a new managed restore copy below.".into());
                },
                (Ok(wire::CatalogBackupSuccess::Restored(receipt)),Action::Restore)=>{
                    restoring.as_ref().expect("Prepared restore review").accept(&receipt,&local)?;
                    restored.set(Some(receipt));status.set("New managed recovery copy created. The current saved catalog and working editor remain unchanged; no automatic switch was made.".into());
                },
                (Ok(wire::CatalogBackupSuccess::RecoveriesListed(page)),Action::ListRecoveries{offset,..})=>{
                    for receipt in page.entries.iter(){RecoveryChoice::new(receipt.clone(),&local)?;}
                    recoveries.set(Some(RecoveryListing{page,offset:*offset}));status.set("Choose a known checked recovery copy to review before adoption.".into());
                },
                (Ok(wire::CatalogBackupSuccess::Reviewed(value)),Action::ReviewOriginal)=>{
                    let checked=AdoptionReview::original(value,&local)?;adoption_review.set(Some(checked));status.set("Review the checked original saved library and explicitly confirm rollback below. Nothing is active yet.".into());
                },
                (Ok(wire::CatalogBackupSuccess::Reviewed(value)),Action::ReviewRecovery(_))=>{
                    let checked=recovery_choice.as_ref().expect("prepared recovery choice").accept(value,&local)?;
                    adoption_review.set(Some(checked));status.set("Review the exact recovery and current catalog below. The recovery library is not active.".into());
                },
                (Ok(wire::CatalogBackupSuccess::AdoptionPending(value)),Action::AdoptRecovery)=>{
                    adopting.as_ref().expect("prepared adoption review").accept_pending(&value,&local)?;
                    activation.set(Some(wire::CatalogActivationStatus::PendingRestart(value)));
                    status.set("Recovery library selected for the next desktop start. Close and reopen Collagene to use it. Library saving and synchronization are paused; the current editor and its autosave remain available.".into());
                    on_catalog_pending.run(true);
                },
                (Err(error),_)=>return Err(error),
                _=>return Err("Native backup response did not match the selected action. Refresh known backups before retrying.".into()),
            }
                Ok(())
            })();
            if let Err(error) = handled {
                if matches!(&action, Action::AdoptRecovery) {
                    // Acceptance may be durable even when reply delivery failed.
                    // Pause old-context actions until explicit native status lookup.
                    on_catalog_pending.run(true);
                    activation_checked.set(true);
                    status.set(format!("Could not confirm library selection. Saved-library controls are paused; check activation status before continuing. {error}"));
                    return;
                }
                if matches!(action, Action::List { .. }) {
                    listing.set(None);
                }
                status.set(error);
            }
        });
    });
    let probe = Callback::new(move |()| {
        if !crate::platform::desktop_rpc::is_desktop()
            || busy.get_untracked()
            || locked.get_untracked()
        {
            return;
        }
        let life = life.get_value();
        let token = life.get().wrapping_add(1);
        life.set(token);
        busy.set(true);
        catalog_status.set(None);
        spawn_local(async move {
            let response = crate::platform::desktop_rpc::call(wire::Command::Bootstrap).await;
            if busy.try_get_untracked().is_none() {
                return;
            }
            busy.set(false);
            if life.get() != token {
                return;
            }
            match response {
                Ok(wire::Success::Bootstrap(value)) => {
                    let paused = matches!(&value.catalog, wire::CatalogStatus::Paused { .. });
                    match &value.catalog {
                        wire::CatalogStatus::Ready { activation: value } => {
                            activation.set(Some(value.clone()))
                        }
                        wire::CatalogStatus::Paused {
                            message,
                            activation: value,
                        } => {
                            activation.set(value.clone());
                            status.set(message.as_str().to_owned());
                        }
                    }
                    catalog_status.set(Some(value.catalog));
                    activation_checked.set(true);
                    on_catalog_pending.run(paused);
                }
                Err(e) => {
                    status.set(e);
                    activation_checked.set(true);
                    on_catalog_pending.run(true);
                }
                _ => {
                    on_catalog_pending.run(true);
                    status.set(
                        "Unexpected catalog activation status; library controls remain paused."
                            .into(),
                    );
                    activation_checked.set(true);
                }
            }
        });
    });
    Effect::new(move |_| {
        if crate::platform::desktop_rpc::is_desktop()
            && !locked.get()
            && !activation_checked.get_untracked()
            && !busy.get()
        {
            activation_checked.set(true);
            probe.run(());
        }
    });
    let cancel_pending = Callback::new(move |()| {
        if busy.get_untracked() || locked.get_untracked() {
            return;
        }
        let pending = match activation.get_untracked() {
            Some(wire::CatalogActivationStatus::PendingRestart(value)) => value,
            _ => return,
        };
        let prepared = CancellationPin::new(pending.clone()).and_then(|pin| {
            let command = pin.command(&pending)?;
            Ok((pin, command))
        });
        let (pin, command) = match prepared {
            Ok(value) => value,
            Err(error) => {
                status.set(error);
                return;
            }
        };
        let life = life.get_value();
        let token = life.get().wrapping_add(1);
        life.set(token);
        busy.set(true);
        spawn_local(async move {
            let response = client.run(command).await;
            if busy.try_get_untracked().is_none() {
                return;
            }
            busy.set(false);
            // A native commit can outlive a changed pane. Preserve the paused
            // library state; only a fresh native bootstrap can rebind context.
            if life.get() != token {
                return;
            }
            match response {
                Ok(wire::CatalogBackupSuccess::ActivationCancelled(receipt)) => {
                    match pin.accept(&receipt) {
                        Ok(()) => {
                            activation.set(Some(wire::CatalogActivationStatus::Cancelled(receipt)));
                            status.set("Pending switch cancelled. Both saved libraries and working-editor files were preserved. Restart Collagene before reopening the previous saved library.".into());
                        }
                        Err(error) => status.set(error),
                    }
                }
                Err(error) => status.set(format!("Could not confirm cancellation. Saved-library controls remain paused; checking native status. {error}")),
                _ => status.set("Unexpected cancellation reply; checking native status. Library controls remain paused.".into()),
            }
            // Handles both a lost accepted reply and a definitive rejection.
            probe.run(());
        });
    });
    view! {<Show when=crate::platform::desktop_rpc::is_desktop><details class="catalog-backup-controls"><summary>"Saved catalog backups"</summary>
        <p class="helper">"Back up the saved local plasmid/folder catalog, its revisions and durable repository state. Save your work in the local library first. This does not back up the current autosaved working-editor files."</p>
        <div class="method-row-actions"><button class="button secondary" disabled=move||busy.get()||locked.get()||catalog_paused.get()||context.with(Option::is_none) on:click=move|_|run.run(Action::Create)>"Make saved catalog backup"</button><button class="text-tool" disabled=move||busy.get()||locked.get()||catalog_paused.get()||context.with(Option::is_none) on:click=move|_|run.run(Action::List{offset:0,refresh:true})>"Refresh known backups"</button></div>
        <button class="text-tool" disabled=move||busy.get()||locked.get()||catalog_paused.get()||context.with(Option::is_none) on:click=move|_|run.run(Action::ListRecoveries{offset:0,refresh:true})>"Browse managed recovery copies"</button>
        {move||recoveries.get().map(|state|{let next=state.page.next_offset;let previous=state.offset.saturating_sub(25);view!{<section><h4>"Managed recovery copies"</h4><ul>{state.page.entries.iter().map(|receipt|{let selected=receipt.clone();let already_current=current_location.get()==(wire::CatalogLocation::Recovery{id:receipt.operation_id.clone()});view!{<li><p>{format!("{} · {}",receipt.directory_basename.as_str(),date(receipt.backup.created_unix_seconds))}</p>{already_current.then(||view!{<strong>"Current saved library"</strong>})}{receipt_details(&receipt.backup)}<button class="text-tool" disabled=move||busy.get()||locked.get()||catalog_paused.get()||already_current on:click=move|_|run.run(Action::ReviewRecovery(selected.clone()))>"Review this recovered library"</button></li>}}).collect_view()}</ul><div><button class="text-tool" disabled=move||busy.get()||locked.get()||catalog_paused.get()||state.offset==0 on:click=move|_|run.run(Action::ListRecoveries{offset:previous,refresh:false})>"Previous recovery copies"</button><button class="text-tool" disabled=move||busy.get()||locked.get()||catalog_paused.get()||next.is_none() on:click=move|_|if let Some(offset)=next{run.run(Action::ListRecoveries{offset,refresh:false});}>"Next recovery copies"</button></div></section>}})}
        {move||latest.get().map(|entry|{let verify=entry.clone();view!{<section><h4>"Latest backup receipt"</h4><p>{entry.basename.as_str().to_owned()}</p>{entry.receipt.as_ref().map(|receipt|view!{<p>{format!("{} · {}",date(receipt.created_unix_seconds),bytes(receipt.archive_bytes))}</p>})}{entry.receipt.as_ref().map(receipt_details)}<button class="text-tool" disabled=move||busy.get()||locked.get()||catalog_paused.get() on:click=move|_|run.run(Action::Verify(verify.clone()))>"Verify this backup"</button></section>}})}
        {move||listing.get().map(|state|{let previous=state.offset.saturating_sub(25);let next=state.page.next_offset;view!{<section><ul>{state.page.entries.iter().map(|entry|{let verify=entry.clone();view!{<li><strong>{entry.basename.as_str().to_owned()}</strong>{entry.receipt.as_ref().map(|receipt|view!{<div><p>{format!("{} · {}",date(receipt.created_unix_seconds),bytes(receipt.archive_bytes))}</p>{receipt_details(receipt)}</div>}.into_any()).unwrap_or_else(||view!{<p class="helper">"Known unfinished creation entry. No completed receipt is claimed; verify to inspect its native state."</p>}.into_any())}<button class="text-tool" disabled=move||busy.get()||locked.get()||catalog_paused.get() on:click=move|_|run.run(Action::Verify(verify.clone()))>"Verify known entry"</button></li>}}).collect_view()}</ul>
            {state.page.entries.is_empty().then(||view!{<p class="helper">"No known backups on this page."</p>})}
            <div class="method-row-actions"><button class="text-tool" disabled=move||busy.get()||locked.get()||catalog_paused.get()||state.offset==0 on:click=move|_|run.run(Action::List{offset:previous,refresh:false})>"Previous backups"</button><button class="text-tool" disabled=move||busy.get()||locked.get()||catalog_paused.get()||next.is_none() on:click=move|_|if let Some(offset)=next{run.run(Action::List{offset,refresh:false});}>"Next backups"</button></div>
        </section>}})}
        {move||review.get().map(|verified|{let receipt=verified.receipt();view!{<section class="catalog-backup-restore-review"><h4>"Confirm restore to a new managed copy"</h4><p>{format!("{} · {} · {}",receipt.basename.as_str(),date(receipt.created_unix_seconds),bytes(receipt.archive_bytes))}</p>{receipt_details(receipt)}
            <p>"Create a separate managed recovery copy from this verified backup. This does not overwrite, switch or merge the current catalog, and does not replace your working editor. Backup verification checks storage integrity; it does not verify cloning designs or experimental outcomes."</p>
            <button class="button secondary" disabled=move||busy.get()||locked.get()||catalog_paused.get() on:click=move|_|run.run(Action::Restore)>"Confirm restore to new copy"</button><button class="text-tool" disabled=move||busy.get() on:click=move|_|{review.set(None);status.set("Restore review cancelled. Current catalog unchanged.".into());}>"Cancel restore review"</button>
        </section>}})}
        {move||restored.get().map(|receipt|{let selected=receipt.clone();view!{<section><h4>"Recovery copy receipt"</h4><p>{format!("Managed recovery copy {} · database {}",receipt.directory_basename.as_str(),receipt.database_basename.as_str())}</p><p class="helper">"The copy is not active. Review it before choosing it for the next desktop start."</p>{receipt_details(&receipt.backup)}<button class="button secondary" disabled=move||busy.get()||locked.get()||catalog_paused.get() on:click=move|_|run.run(Action::ReviewRecovery(selected.clone()))>"Review recovered library"</button></section>}})}
        {move||adoption_review.get().map(|checked|{let selected=checked.review();view!{<section><h4>"Use this saved library after restart?"</h4><p>{format!("Target: {} · current: {} · activation generation {}",location_name(&selected.target),location_name(&selected.current_catalog),selected.activation_generation)}</p>{snapshot_details(&selected.target_snapshot)}{selected.target_recovery.as_ref().map(|receipt|view!{<details><summary>"Original backup source for this managed copy"</summary>{receipt_details(&receipt.backup)}</details>})}<p>"Choose this exact checked saved library for the next desktop start. Your current working editor, recovery files and credential storage stay in their existing private root. The current saved library is preserved for a later explicit rollback. Catalogs are not merged, and this does not validate cloning designs."</p><button class="button primary" disabled=move||busy.get()||locked.get()||catalog_paused.get() on:click=move|_|run.run(Action::AdoptRecovery)>"Use reviewed library after restart"</button><button class="text-tool" disabled=move||busy.get() on:click=move|_|adoption_review.set(None)>"Cancel library selection"</button></section>}})}
        {move||activation.get().map(|state|match state{
            wire::CatalogActivationStatus::None=>view!{<p class="helper">"The original saved library is active."</p>}.into_any(),
            wire::CatalogActivationStatus::PendingRestart(receipt)=>view!{<section role="status"><h4>"Recovery library pending restart"</h4><p>{format!("Selected library: {}",location_name(&receipt.review.target))}</p><p>"Close and reopen Collagene to use this library. It is not active in this window. Library saving and synchronization are paused; working-editor autosave remains available."</p><button class="button secondary" disabled=move||busy.get()||locked.get() on:click=move|_|cancel_pending.run(())>"Cancel pending switch; keep previous library"</button></section>}.into_any(),
            wire::CatalogActivationStatus::Cancelled(_)=>view!{<section role="status"><h4>"Library switch cancelled"</h4><p>"Both saved libraries were preserved. The previous library remains selected. If library controls are paused, close and reopen Collagene; working-editor autosave remains available."</p></section>}.into_any(),
            wire::CatalogActivationStatus::ActiveAfterRestart(receipt)=>view!{<section role="status"><h4>"Recovered library active"</h4><p>{format!("Native startup opened {}",location_name(&receipt.review.target))}</p><p class="helper">"The original saved catalog remains preserved. Your working editor has not been replaced."</p></section>}.into_any(),
        })}
        {move||catalog_status.get().and_then(|state|match state{wire::CatalogStatus::Paused{message,..}=>Some(view!{<p role="alert">{message.as_str().to_owned()}</p>}),_=>None})}
        <button class="text-tool" disabled=move||busy.get()||locked.get()||catalog_paused.get()||context.with(Option::is_none)||current_location.get()==wire::CatalogLocation::Original on:click=move|_|run.run(Action::ReviewOriginal)>"Review rollback to original saved library"</button>
        <button class="text-tool" disabled=move||busy.get()||locked.get() on:click=move|_|probe.run(())>"Check library activation status"</button>
        <p role="status" aria-live="polite">{move||status.get()}</p>
    </details></Show>}
}
