//! Native-only saved-catalog backup controls. Commands carry known IDs and pins,
//! never a filesystem path or biological backup archive.
use crate::catalog_backup_model::{RestoreReview, VerificationPin};
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
enum Action {
    Create,
    List { offset: u32, refresh: bool },
    Verify(wire::CatalogBackupEntry),
    Restore,
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
#[component]
pub fn CatalogBackupControls(
    client: Callback<wire::CatalogBackupCommand, LocalResult<wire::CatalogBackupSuccess>>,
    context: Memo<Option<library::LabContext>>,
    locked: Memo<bool>,
) -> impl IntoView {
    let life = StoredValue::new_local(Rc::new(Cell::new(0u64)));
    let observed = StoredValue::new(None::<library::LabContext>);
    let busy = RwSignal::new(false);
    let status = RwSignal::new(String::new());
    let listing = RwSignal::new(None::<Listing>);
    let latest = RwSignal::new(None::<wire::CatalogBackupEntry>);
    let review = RwSignal::new(None::<RestoreReview>);
    let restored = RwSignal::new(None::<wire::CatalogRecoveryReceipt>);
    Effect::new(move |_| {
        let now = context.get();
        let blocked = locked.get();
        if observed.get_value() != now || blocked {
            observed.set_value(now);
            life.with_value(|life| life.set(life.get().wrapping_add(1)));
            listing.set(None);
            latest.set(None);
            review.set(None);
            restored.set(None);
        }
    });
    on_cleanup(move || life.with_value(|life| life.set(life.get().wrapping_add(1))));
    let run = Callback::new(move |action: Action| {
        if !crate::desktop_rpc::is_desktop() || busy.get_untracked() || locked.get_untracked() {
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
            let command = match &action {
                Action::Create => wire::CatalogBackupCommand::Create,
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
            Ok((local, command, verification, restoring))
        })();
        let (local, command, verification, restoring) = match setup {
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
        let life = life.get_value();
        let epoch = life.get().wrapping_add(1);
        life.set(epoch);
        busy.set(true);
        status.set(
            match &action {
                Action::Create => "Making a streaming backup of the saved local catalog…",
                Action::List { .. } => "Loading known managed backups…",
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
                (Err(error),_)=>return Err(error),
                _=>return Err("Native backup response did not match the selected action. Refresh known backups before retrying.".into()),
            }
                Ok(())
            })();
            if let Err(error) = handled {
                if matches!(action, Action::List { .. }) {
                    listing.set(None);
                }
                status.set(error);
            }
        });
    });
    view! {<Show when=crate::desktop_rpc::is_desktop><details class="catalog-backup-controls"><summary>"Saved catalog backups"</summary>
        <p class="helper">"Back up the saved local plasmid/folder catalog, its revisions and durable repository state. Save your work in the local library first. This does not back up the current autosaved working-editor files."</p>
        <div class="method-row-actions"><button class="button secondary" disabled=move||busy.get()||locked.get()||context.with(Option::is_none) on:click=move|_|run.run(Action::Create)>"Make saved catalog backup"</button><button class="text-tool" disabled=move||busy.get()||locked.get()||context.with(Option::is_none) on:click=move|_|run.run(Action::List{offset:0,refresh:true})>"Refresh known backups"</button></div>
        {move||latest.get().map(|entry|{let verify=entry.clone();view!{<section><h4>"Latest backup receipt"</h4><p>{entry.basename.as_str().to_owned()}</p>{entry.receipt.as_ref().map(|receipt|view!{<p>{format!("{} · {}",date(receipt.created_unix_seconds),bytes(receipt.archive_bytes))}</p>})}{entry.receipt.as_ref().map(receipt_details)}<button class="text-tool" disabled=move||busy.get()||locked.get() on:click=move|_|run.run(Action::Verify(verify.clone()))>"Verify this backup"</button></section>}})}
        {move||listing.get().map(|state|{let previous=state.offset.saturating_sub(25);let next=state.page.next_offset;view!{<section><ul>{state.page.entries.iter().map(|entry|{let verify=entry.clone();view!{<li><strong>{entry.basename.as_str().to_owned()}</strong>{entry.receipt.as_ref().map(|receipt|view!{<div><p>{format!("{} · {}",date(receipt.created_unix_seconds),bytes(receipt.archive_bytes))}</p>{receipt_details(receipt)}</div>}.into_any()).unwrap_or_else(||view!{<p class="helper">"Known unfinished creation entry. No completed receipt is claimed; verify to inspect its native state."</p>}.into_any())}<button class="text-tool" disabled=move||busy.get()||locked.get() on:click=move|_|run.run(Action::Verify(verify.clone()))>"Verify known entry"</button></li>}}).collect_view()}</ul>
            {state.page.entries.is_empty().then(||view!{<p class="helper">"No known backups on this page."</p>})}
            <div class="method-row-actions"><button class="text-tool" disabled=move||busy.get()||locked.get()||state.offset==0 on:click=move|_|run.run(Action::List{offset:previous,refresh:false})>"Previous backups"</button><button class="text-tool" disabled=move||busy.get()||locked.get()||next.is_none() on:click=move|_|if let Some(offset)=next{run.run(Action::List{offset,refresh:false});}>"Next backups"</button></div>
        </section>}})}
        {move||review.get().map(|verified|{let receipt=verified.receipt();view!{<section class="catalog-backup-restore-review"><h4>"Confirm restore to a new managed copy"</h4><p>{format!("{} · {} · {}",receipt.basename.as_str(),date(receipt.created_unix_seconds),bytes(receipt.archive_bytes))}</p>{receipt_details(receipt)}
            <p>"Create a separate managed recovery copy from this verified backup. This does not overwrite, switch or merge the current catalog, and does not replace your working editor. Backup verification checks storage integrity; it does not verify cloning designs or experimental outcomes."</p>
            <button class="button secondary" disabled=move||busy.get()||locked.get() on:click=move|_|run.run(Action::Restore)>"Confirm restore to new copy"</button><button class="text-tool" disabled=move||busy.get() on:click=move|_|{review.set(None);status.set("Restore review cancelled. Current catalog unchanged.".into());}>"Cancel restore review"</button>
        </section>}})}
        {move||restored.get().map(|receipt|view!{<section><h4>"Recovery copy receipt"</h4><p>{format!("Managed recovery copy {} · database {}",receipt.directory_basename.as_str(),receipt.database_basename.as_str())}</p><p class="helper">"Current catalog remains active. No automatic cutover was performed."</p>{receipt_details(&receipt.backup)}</section>})}
        <p role="status" aria-live="polite">{move||status.get()}</p>
    </details></Show>}
}
