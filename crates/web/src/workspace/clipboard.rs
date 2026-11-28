//! Clipboard actions apply only to the document and selection captured by the user gesture.
use crate::platform::worker::CachedWorker;
use crate::sources::identity::SourceInfo;
use crate::workspace::reference_transfer_view::Copied;
use crate::workspace::reference_transfer_view::PasteOpen;
use collagene_core::{
    protocol::{
        request::Request,
        schema::{AnalysisJob, AnalysisResult},
    },
    sequence::reference_transfer::{self, RequestCopy, Selection},
};
use collagene_model::{Alphabet, Molecule, Sequence};
use leptos::prelude::*;
use std::rc::Rc;
use std::sync::{
    Arc,
    atomic::{AtomicBool, Ordering},
};
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
        return Err("Menu clipboard access requires HTTPS or localhost. Focus the sequence and use its Copy/Cut/Paste shortcuts, or use the internal annotated-reference controls.".into());
    }
    let value = window.navigator().clipboard();
    if value.is_undefined() || value.is_null() {
        return Err("This browser does not expose menu clipboard access. Use the focused sequence viewer shortcuts, internal annotated-reference controls, or the edit dialog.".into());
    }
    Ok(value)
}

#[component]
pub fn ClipboardTools(
    molecule: Memo<Molecule>,
    selection: RwSignal<Option<(usize, usize)>>,
    selected_feature: RwSignal<Option<usize>>,
    highlighted_spans: Memo<Vec<(usize, usize)>>,
    locked: RwSignal<bool>,
    epoch: RwSignal<u64>,
    identity: Memo<String>,
    on_commit: Callback<(String, Molecule), Result<(), String>>,
    message: RwSignal<String>,
    source: Memo<SourceInfo>,
    source_selection: Callback<(), Result<Selection, String>>,
    internal_clipboard: RwSignal<Option<Arc<Copied>>>,
    on_preview: Callback<PasteOpen, Result<(), String>>,
) -> impl IntoView {
    let busy = RwSignal::new(false);
    let generation = RwSignal::new(0u64);
    let rich_worker = StoredValue::new_local(None::<Rc<CachedWorker>>);
    let alive = Arc::new(AtomicBool::new(true));
    {
        let alive = alive.clone();
        on_cleanup(move || {
            alive.store(false, Ordering::Relaxed);
            if let Some(worker) = rich_worker.get_value() {
                worker.close();
            }
        });
    }
    Effect::new(move |_| {
        let _ = (
            source.get().binding,
            selection.get(),
            selected_feature.get(),
            highlighted_spans.get(),
            epoch.get(),
            identity.get(),
        );
        generation.update(|value| *value = value.wrapping_add(1));
        if let Some(worker) = rich_worker.get_value() {
            worker.close();
        }
        rich_worker.set_value(None);
        busy.set(false);
    });
    let copy_reference = {
        let alive = alive.clone();
        Callback::new(move |reverse: bool| {
            if locked.get_untracked() || busy.get_untracked() {
                return;
            }
            let model = molecule.get_untracked();
            let info = source.get_untracked();
            let captured_epoch = epoch.get_untracked();
            let declaration = match source_selection.run(()) {
                Ok(mut selection) => {
                    selection.reverse ^= reverse;
                    selection
                }
                Err(error) => {
                    message.set(error);
                    return;
                }
            };
            if model.revision() != info.binding.revision {
                message.set("The source version changed; select its bases again.".into());
                return;
            }
            let source_molecule = Arc::new(model);
            let request = RequestCopy {
                source: source_molecule.as_ref().clone(),
                selection: declaration.clone(),
            };
            let revision = match reference_transfer::copy_revision(&request) {
                Ok(revision) => revision,
                Err(error) => {
                    message.set(error.to_string());
                    return;
                }
            };
            let w = match CachedWorker::new() {
                Ok(worker) => Rc::new(worker),
                Err(error) => {
                    message.set(error);
                    return;
                }
            };
            if let Some(previous) = rich_worker.get_value() {
                previous.close();
            }
            rich_worker.set_value(Some(w.clone()));
            let token = generation.get_untracked().wrapping_add(1);
            generation.set(token);
            busy.set(true);
            message.set("Copying the selected reference and projecting annotations…".into());
            let alive = alive.clone();
            spawn_local(async move {
                let response = w
                    .ask(Request {
                        id: token as u32,
                        revision: revision.clone(),
                        job: AnalysisJob::CopyReferenceFragment { request },
                    })
                    .await;
                w.close();
                if !alive.load(Ordering::Relaxed) || generation.get_untracked() != token {
                    return;
                }
                rich_worker.set_value(None);
                busy.set(false);
                let outcome = (|| -> Result<_, String> {
                    let mut current = source_selection.run(())?;
                    current.reverse ^= reverse;
                    if locked.get_untracked()
                        || epoch.get_untracked() != captured_epoch
                        || source.get_untracked().binding != info.binding
                        || molecule.with_untracked(Molecule::revision) != info.binding.revision
                        || current != declaration
                    {
                        return Err(
                            "The source or selection changed; no annotated reference was copied."
                                .into(),
                        );
                    }
                    let response = response?;
                    if response.id != token as u32 || response.revision != revision {
                        return Err("Reference-copy response belongs to another selection.".into());
                    }
                    let fragment = match response.result {
                        AnalysisResult::ReferenceFragment(fragment) => fragment,
                        AnalysisResult::Error(error) => return Err(error),
                        _ => return Err("Unexpected reference-copy response.".into()),
                    };
                    fragment
                        .validate_binding(&revision)
                        .map_err(|e| e.to_string())?;
                    if fragment.source_revision != info.binding.revision
                        || fragment.selection != declaration
                    {
                        return Err("The copied reference differs from the exact source or biological selection.".into());
                    }
                    Ok(fragment)
                })();
                match outcome {
                    Ok(fragment) => {
                        let length = fragment.molecule.sequence.len();
                        let count = fragment.molecule.features.len();
                        internal_clipboard.set(Some(Arc::new(Copied {
                            source: info.binding,
                            source_molecule,
                            fragment: Arc::new(fragment),
                            epoch: captured_epoch,
                        })));
                        message.set(format!("Copied {length} reference bases and {count} projected annotations. Choose a destination range or caret, then preview annotated paste."));
                    }
                    Err(error) => message.set(error),
                }
            });
        })
    };
    let preview_reference = move |_| {
        if locked.get_untracked() || busy.get_untracked() {
            return;
        }
        let result = (|| -> Result<(), String> {
            let copied = internal_clipboard
                .get_untracked()
                .ok_or("Copy a reference with annotations first.")?;
            let captured_epoch = epoch.get_untracked();
            if copied.epoch != captured_epoch {
                return Err("The workspace changed; copy the source reference again.".into());
            }
            if selected_feature.get_untracked().is_some() {
                return Err("Choose contiguous destination bases or an insertion caret; an annotation is not a replacement interval.".into());
            }
            let (start, end) = selection
                .get_untracked()
                .ok_or("Select destination bases to replace or an insertion caret.")?;
            let destination = source.get_untracked().binding;
            let model = molecule.get_untracked();
            if destination.revision != model.revision() || start > end || end > model.sequence.len()
            {
                return Err("The destination version or selected interval changed.".into());
            }
            on_preview.run(PasteOpen {
                copied,
                destination,
                destination_molecule: Arc::new(model),
                start,
                end,
                epoch: captured_epoch,
            })
        })();
        if let Err(error) = result {
            message.set(error);
        }
    };
    let clear_reference = move |_| {
        generation.update(|value| *value = value.wrapping_add(1));
        if let Some(worker) = rich_worker.get_value() {
            worker.close();
        }
        rich_worker.set_value(None);
        busy.set(false);
        internal_clipboard.set(None);
        message.set("Internal annotated reference clipboard cleared.".into());
    };
    let run = Callback::new(move |action: Action| {
        if locked.get_untracked() || busy.get_untracked() {
            return;
        }
        let source = molecule.get_untracked();
        let revision = source.revision();
        let span = selection.get_untracked();
        let feature = selected_feature.get_untracked();
        let source_parts = highlighted_spans.get_untracked();
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
            let text =
                crate::workspace::clipboard_model::copy_text(&source, feature, span, &source_parts)
                    .and_then(|text| {
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
                    message.set("System clipboard access failed. No sequence edit was made; use the internal annotated-reference controls or the edit dialog for external text.".into());
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
                || highlighted_spans.get_untracked() != source_parts
            {
                message.set("The document or selection changed while waiting for the clipboard. No sequence edit was made.".into());
                return;
            }
            let replacement = if matches!(action, Action::Cut) {
                String::new()
            } else {
                value.as_string().unwrap_or_default()
            };
            let edit = match action {
                Action::Cut => crate::workspace::clipboard_model::cut(
                    &source,
                    feature,
                    span,
                    &source_parts,
                    false,
                    false,
                ),
                Action::Insert => crate::workspace::clipboard_model::insert(
                    &source,
                    feature,
                    span,
                    &source_parts,
                    false,
                    false,
                    &replacement,
                ),
                _ => crate::workspace::clipboard_model::paste(
                    &source,
                    feature,
                    span,
                    &source_parts,
                    false,
                    false,
                    &replacement,
                ),
            };
            let edited = edit.and_then(|edit| {
                source
                    .replace(edit.start, edit.end, &edit.replacement)
                    .map_err(|error| error.to_string())
            });
            match edited.and_then(|edited| on_commit.run((revision, edited))) {
                Ok(()) => {
                    selected_feature.set(None);
                    selection.set(None);
                    message.set(if matches!(action,Action::Cut){"Selected reference bases copied and cut. Undo restores the previous document.".into()}else{"Clipboard sequence applied. Affected annotations require review; Undo restores the previous document.".into()});
                }
                Err(error) => message.set(format!("No sequence edit was made. {error}")),
            }
        });
    });
    view! {<details class="clipboard-tools"><summary>"Clipboard"</summary><div class="clipboard-actions"><p class="helper">"Annotation copying follows its biological orientation and compound parts. Cut and paste edit an explicitly selected reference span."</p><button class="text-tool" disabled=move||busy.get()||locked.get()||(selection.get().is_none()&&selected_feature.get().is_none()) on:click=move|_|run.run(Action::Copy)>"Copy selection"</button><button class="text-tool" disabled=move||busy.get()||locked.get()||molecule.with(|m|m.sequence.alphabet()==Alphabet::Protein)||(selection.get().is_none()&&selected_feature.get().is_none()) on:click=move|_|run.run(Action::CopyReverse)>"Copy reverse complement"</button><button class="text-tool" disabled=move||busy.get()||locked.get()||selection.get().is_none()||selected_feature.get().is_some() on:click=move|_|run.run(Action::Cut)>"Cut reference selection"</button><button class="text-tool" disabled=move||busy.get()||locked.get()||selection.get().is_none()||selected_feature.get().is_some() on:click=move|_|run.run(Action::Paste)>"Paste replacing selection"</button><button class="text-tool" disabled=move||busy.get()||locked.get()||selection.get().is_none()||selected_feature.get().is_some() on:click=move|_|run.run(Action::Insert)>"Paste before selection"</button><p class="helper">"Ctrl/Cmd+C, X and V use the focused sequence viewer. These menu buttons use browser clipboard access, which requires permission and HTTPS or localhost."</p>
        <h4>"Annotated reference clipboard"</h4><p class="helper">"Copy reference bases and annotations within this workspace, without browser clipboard permission. Pasting is an explicit reference edit, not a cloning reaction."</p>
        <button class="text-tool" disabled=move||busy.get()||locked.get() on:click=move|_|copy_reference.run(false)>"Copy with annotations"</button>
        <button class="text-tool" disabled=move||busy.get()||locked.get()||molecule.with(|model|model.sequence.alphabet()==Alphabet::Protein) on:click=move|_|copy_reference.run(true)>"Copy reverse with annotations"</button>
        <button class="text-tool" disabled=move||busy.get()||locked.get()||internal_clipboard.get().is_none()||selection.get().is_none()||selected_feature.get().is_some() on:click=preview_reference>"Paste with annotations…"</button>
        <button class="text-tool" disabled=move||locked.get()||internal_clipboard.get().is_none() on:click=clear_reference>"Clear copied reference"</button>
        <p class="helper" style="overflow-wrap:anywhere">{move||internal_clipboard.with(|copied|copied.as_ref().map(|copied|format!("Copied from {} · {} bases · {} annotations",copied.source_molecule.name,copied.fragment.molecule.sequence.len(),copied.fragment.molecule.features.len())).unwrap_or_else(||"No annotated reference copied.".into()))}</p>
        </div></details>}.into_any()
}
