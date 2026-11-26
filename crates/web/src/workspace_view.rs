//! Native archive export and atomic browser replacement with a recovery slot.
use crate::{
    app::{decode_workspace, download, encode_workspace},
    storage,
};
use collagene_core::workspace::Workspace;
use leptos::{html, prelude::*};
use std::sync::{
    Arc,
    atomic::{AtomicBool, Ordering},
};
use wasm_bindgen_futures::{JsFuture, spawn_local};

#[component]
pub fn WorkspacePanel(
    snapshot: Callback<(), Result<Workspace, String>>,
    validate: Callback<Workspace, Result<(), String>>,
    on_restore: Callback<Workspace>,
    db: StoredValue<Option<crate::storage::Database>, LocalStorage>,
    busy: RwSignal<bool>,
    loaded: RwSignal<bool>,
    save_generation: RwSignal<u64>,
    message: RwSignal<String>,
) -> impl IntoView {
    let alive = Arc::new(AtomicBool::new(true));
    let cleanup = alive.clone();
    on_cleanup(move || cleanup.store(false, Ordering::Relaxed));
    let input = NodeRef::<html::Input>::new();
    let status = RwSignal::new(String::new());
    let exporting = RwSignal::new(false);
    // This callback locks editing before any await. Validate the new UI state
    // before changing either storage slot; commit both slots in one transaction.
    let replace = Callback::new(move |bytes: Vec<u8>| {
        if busy.get_untracked() {
            return;
        }
        let Some(database) = db.get_value() else {
            status.set("Browser storage is unavailable; export your current work before reopening a workspace in a browser with storage enabled.".into());
            return;
        };
        let backup = match snapshot.run(()) {
            Ok(frame) => frame,
            Err(error) => {
                status.set(error);
                return;
            }
        };
        busy.set(true);
        save_generation.update(|g| *g = g.wrapping_add(1));
        status.set("Checking workspace and retaining your previous workspace…".into());
        spawn_local(async move {
            let result = async {
                let frame = decode_workspace(bytes.clone(), u32::MAX - 1).await?;
                validate.run(frame.clone())?;
                let previous = encode_workspace(backup, u32::MAX - 2).await?;
                storage::replace_workspace(&database, &bytes, &previous).await?;
                Ok::<_, String>(frame)
            }
            .await;
            match result {
                Ok(frame) => {
                    status.set(String::new());
                    loaded.set(true);
                    message.set("Workspace restored. The previous workspace is available in Workspace recovery. Candidate plans require fresh verification.".into());
                    on_restore.run(frame);
                }
                Err(error) => status.set(format!("Workspace was not replaced. {error}")),
            }
            busy.set(false);
        });
    });
    let open_alive = alive.clone();
    let open = move |event: web_sys::Event| {
        let target = event_target::<web_sys::HtmlInputElement>(&event);
        let file = target.files().and_then(|files| files.get(0));
        target.set_value("");
        let Some(file) = file else {
            return;
        };
        if file.size() > 64.0 * 1024.0 * 1024.0 {
            status.set("Workspace archive exceeds 64 MiB.".into());
            return;
        }
        let alive = open_alive.clone();
        spawn_local(async move {
            let result = JsFuture::from(file.array_buffer()).await;
            if !alive.load(Ordering::Relaxed) {
                return;
            }
            match result {
                Ok(buffer) => replace.run(js_sys::Uint8Array::new(&buffer).to_vec()),
                Err(_) => status.set("Workspace file could not be read.".into()),
            }
        });
    };
    let export_alive = alive.clone();
    let export = move |_| {
        if exporting.get_untracked() || busy.get_untracked() {
            return;
        }
        let frame = match snapshot.run(()) {
            Ok(frame) => frame,
            Err(error) => {
                status.set(error);
                return;
            }
        };
        exporting.set(true);
        status.set("Preparing workspace file…".into());
        let alive = export_alive.clone();
        spawn_local(async move {
            let encoded = encode_workspace(frame, u32::MAX - 3).await;
            if !alive.load(Ordering::Relaxed) {
                return;
            }
            let result = match encoded {
                Ok(bytes) => download(
                    "collagene-workspace.cg2workspace",
                    "application/zip",
                    &bytes,
                ),
                Err(error) => Err(error),
            };
            status.set(match result{Ok(())=>"Workspace file prepared with document history, candidates and original evidence.".into(),Err(error)=>error});
            exporting.set(false);
        });
    };
    let recover = move |_| {
        let Some(database) = db.get_value() else {
            status.set("Browser storage is unavailable.".into());
            return;
        };
        let alive = alive.clone();
        spawn_local(async move {
            let result = storage::load_previous(&database).await;
            if !alive.load(Ordering::Relaxed) {
                return;
            }
            match result {
                Ok(Some(bytes)) => replace.run(bytes),
                Ok(None) => status.set("No previous workspace has been retained yet.".into()),
                Err(error) => status.set(error),
            }
        });
    };
    view! {
        <p>"Save all open documents, revision branches, undo/redo, candidate plans and original evidence in one portable file."</p>
        <p class="helper">"Opening a workspace replaces the open documents after validation. Your current workspace is retained as the previous workspace in this browser. Imported candidates must be verified again."</p>
        <button class="button primary" disabled=move||busy.get()||exporting.get() on:click=export>"Save workspace file"</button>
        <button class="button secondary" disabled=move||busy.get() on:click=move|_|if let Some(element)=input.get(){element.click();}>"Open workspace file"</button>
        <input node_ref=input type="file" accept=".cg2workspace,application/zip" style="display:none" on:change=open/>
        <button class="button secondary" disabled=move||busy.get() on:click=recover>"Restore previous workspace"</button>
        <p role="status">{move||status.get()}</p>
    }
}
