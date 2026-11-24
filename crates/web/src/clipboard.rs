//! Clipboard actions apply only to the document and selection captured by the user gesture.
use collagene_core::{Alphabet, Molecule, Sequence};
use leptos::prelude::*;
use std::sync::{
    Arc,
    atomic::{AtomicBool, Ordering},
};
use wasm_bindgen::{JsCast, closure::Closure};
use wasm_bindgen_futures::{JsFuture, spawn_local};

#[derive(Clone, Copy)]
enum Action {
    Copy,
    CopyReverse,
    Cut,
    Paste,
    Insert,
}
fn clipboard() -> Result<web_sys::Clipboard, String> {
    let window = web_sys::window().ok_or("No browser window")?;
    if !window.is_secure_context() {
        return Err("Clipboard access requires HTTPS or localhost. Use the sequence/edit dialogs to paste text here.".into());
    }
    let value = window.navigator().clipboard();
    if value.is_undefined() || value.is_null() {
        return Err("This browser does not expose clipboard access. Use the sequence/edit dialogs to paste text here.".into());
    }
    Ok(value)
}
fn copy_text(
    molecule: &Molecule,
    feature: Option<usize>,
    span: Option<(usize, usize)>,
) -> Result<String, String> {
    if let Some(index) = feature {
        return molecule
            .features
            .get(index)
            .ok_or("Selected annotation changed")?
            .location
            .extract(&molecule.sequence)
            .map_err(|error| error.to_string());
    }
    let (a, b) = span.ok_or("Select bases or an annotation to copy")?;
    molecule
        .sequence
        .as_str()
        .get(a..b)
        .filter(|s| !s.is_empty())
        .map(str::to_string)
        .ok_or("Selected bases are no longer available".into())
}
struct KeyListener {
    window: web_sys::Window,
    callback: Closure<dyn FnMut(web_sys::KeyboardEvent)>,
}
impl Drop for KeyListener {
    fn drop(&mut self) {
        let _ = self
            .window
            .remove_event_listener_with_callback("keydown", self.callback.as_ref().unchecked_ref());
    }
}

#[component]
pub fn ClipboardTools(
    molecule: Memo<Molecule>,
    selection: RwSignal<Option<(usize, usize)>>,
    selected_feature: RwSignal<Option<usize>>,
    locked: RwSignal<bool>,
    epoch: RwSignal<u64>,
    identity: Memo<String>,
    on_commit: Callback<(String, Molecule), Result<(), String>>,
    message: RwSignal<String>,
) -> impl IntoView {
    let busy = RwSignal::new(false);
    let generation = RwSignal::new(0u64);
    let alive = Arc::new(AtomicBool::new(true));
    {
        let alive = alive.clone();
        on_cleanup(move || alive.store(false, Ordering::Relaxed));
    }
    let run = Callback::new(move |action: Action| {
        if locked.get_untracked() || busy.get_untracked() {
            return;
        }
        let source = molecule.get_untracked();
        let revision = source.revision();
        let span = selection.get_untracked();
        let feature = selected_feature.get_untracked();
        let source_epoch = epoch.get_untracked();
        let source_identity = identity.get_untracked();
        let writing = matches!(action, Action::Copy | Action::CopyReverse | Action::Cut);
        if !writing && feature.is_some() {
            message.set("Choose contiguous reference bases before pasting; a compound biological annotation is not a replacement span.".into());
            return;
        }
        if matches!(action, Action::Cut) && feature.is_some() {
            message.set("Choose contiguous reference bases to cut. Annotation copying preserves biological order; cutting needs an explicit reference interval.".into());
            return;
        }
        let clipboard = match clipboard() {
            Ok(value) => value,
            Err(error) => {
                message.set(error);
                return;
            }
        };
        let text = if writing {
            let text = copy_text(&source, feature, span).and_then(|text| {
                if matches!(action, Action::CopyReverse) {
                    let sequence = match source.sequence.alphabet() {
                        Alphabet::Dna => Sequence::dna(&text),
                        Alphabet::Rna => Sequence::new(Alphabet::Rna, &text),
                        Alphabet::Protein => {
                            return Err("Proteins do not have a reverse complement.".into());
                        }
                    }
                    .map_err(|e| e.to_string())?;
                    sequence
                        .reverse_complement()
                        .map(|s| s.as_str().to_string())
                        .map_err(|e| e.to_string())
                } else {
                    Ok(text)
                }
            });
            match text {
                Ok(value) => Some(value),
                Err(error) => {
                    message.set(error);
                    return;
                }
            }
        } else {
            if span.is_none() {
                message.set("Select reference bases to replace, or select an insertion boundary using the sequence view.".into());
                return;
            }
            None
        };
        // Clipboard access starts synchronously in the actual user gesture.
        let promise = if let Some(text) = text.as_ref() {
            clipboard.write_text(text)
        } else {
            clipboard.read_text()
        };
        let id = generation.get_untracked().wrapping_add(1);
        generation.set(id);
        busy.set(true);
        let alive = alive.clone();
        spawn_local(async move {
            let result = JsFuture::from(promise).await;
            if !alive.load(Ordering::Relaxed) || generation.get_untracked() != id {
                return;
            }
            busy.set(false);
            let value = match result {
                Ok(value) => value,
                Err(_) => {
                    message.set("Clipboard access failed. No sequence edit was made; use the paste/edit dialogs if clipboard permission is unavailable.".into());
                    return;
                }
            };
            if matches!(action, Action::Copy | Action::CopyReverse) {
                message.set("Selected sequence copied as plain 5′→3′ text.".into());
                return;
            }
            if locked.get_untracked()
                || epoch.get_untracked() != source_epoch
                || identity.get_untracked() != source_identity
                || molecule.with_untracked(Molecule::revision) != revision
                || selection.get_untracked() != span
                || selected_feature.get_untracked() != feature
            {
                message.set("The document or selection changed while waiting for the clipboard. No sequence edit was made.".into());
                return;
            }
            let (a, b) = span.expect("edit span checked");
            let replacement = if matches!(action, Action::Cut) {
                String::new()
            } else {
                value.as_string().unwrap_or_default()
            };
            if !matches!(action, Action::Cut)
                && replacement.len() > collagene_core::formats::MAX_IMPORT_BYTES
            {
                message.set(
                    "Clipboard sequence exceeds the 32 MiB input limit; no sequence edit was made."
                        .into(),
                );
                return;
            }
            if !matches!(action, Action::Cut) && replacement.trim().is_empty() {
                message.set("Clipboard contains no sequence; no sequence edit was made.".into());
                return;
            }
            let end = if matches!(action, Action::Insert) {
                a
            } else {
                b
            };
            match source
                .replace(a, end, &replacement)
                .map_err(|error| error.to_string())
                .and_then(|edited| on_commit.run((revision, edited)))
            {
                Ok(()) => {
                    selected_feature.set(None);
                    selection.set(None);
                    message.set(if matches!(action,Action::Cut){"Selected reference bases copied and cut. Undo restores the previous document.".into()}else{"Clipboard sequence applied. Affected annotations require review; Undo restores the previous document.".into()});
                }
                Err(error) => message.set(format!("No sequence edit was made. {error}")),
            }
        });
    });
    let listener = StoredValue::new_local(None::<KeyListener>);
    if let Some(window) = web_sys::window() {
        let callback = Closure::<dyn FnMut(web_sys::KeyboardEvent)>::new(
            move |event: web_sys::KeyboardEvent| {
                if !(event.ctrl_key() || event.meta_key()) || event.alt_key() {
                    return;
                }
                let Some(target) = event
                    .target()
                    .and_then(|target| target.dyn_into::<web_sys::Element>().ok())
                else {
                    return;
                };
                // Form controls, dialogs and normal selectable DOM text keep browser behavior.
                if target
                    .closest("input,textarea,select,[contenteditable],[role=dialog]")
                    .ok()
                    .flatten()
                    .is_some()
                {
                    return;
                }
                if web_sys::window()
                    .and_then(|window| window.get_selection().ok().flatten())
                    .is_some_and(|selection| selection.to_string().length() != 0)
                {
                    return;
                }
                if !target
                    .get_attribute("class")
                    .unwrap_or_default()
                    .split_ascii_whitespace()
                    .any(|class| class == "sequence-canvas")
                {
                    return;
                }
                let action = match event.key().to_ascii_lowercase().as_str() {
                    "c" => Action::Copy,
                    "x" => Action::Cut,
                    "v" => Action::Paste,
                    _ => return,
                };
                event.prevent_default();
                run.run(action);
            },
        );
        if window
            .add_event_listener_with_callback("keydown", callback.as_ref().unchecked_ref())
            .is_ok()
        {
            listener.set_value(Some(KeyListener { window, callback }));
        }
    }
    on_cleanup(move || listener.set_value(None));
    view! {<details class="clipboard-tools"><summary>"Clipboard"</summary><div class="clipboard-actions"><p class="helper">"Annotation copying follows its biological orientation and compound parts. Cut and paste edit an explicitly selected reference span."</p><button class="text-tool" disabled=move||busy.get()||locked.get()||(selection.get().is_none()&&selected_feature.get().is_none()) on:click=move|_|run.run(Action::Copy)>"Copy selection"</button><button class="text-tool" disabled=move||busy.get()||locked.get()||molecule.with(|m|m.sequence.alphabet()==Alphabet::Protein)||(selection.get().is_none()&&selected_feature.get().is_none()) on:click=move|_|run.run(Action::CopyReverse)>"Copy reverse complement"</button><button class="text-tool" disabled=move||busy.get()||locked.get()||selection.get().is_none()||selected_feature.get().is_some() on:click=move|_|run.run(Action::Cut)>"Cut reference selection"</button><button class="text-tool" disabled=move||busy.get()||locked.get()||selection.get().is_none()||selected_feature.get().is_some() on:click=move|_|run.run(Action::Paste)>"Paste replacing selection"</button><button class="text-tool" disabled=move||busy.get()||locked.get()||selection.get().is_none()||selected_feature.get().is_some() on:click=move|_|run.run(Action::Insert)>"Paste before selection"</button><p class="helper">"Ctrl/Cmd+C, X and V work when the sequence canvas has keyboard focus. Clipboard access requires browser permission and HTTPS or localhost."</p></div></details>}
}
