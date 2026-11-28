//! Owns local storage opening, draft persistence and generation-guarded autosave.
use crate::platform::storage;
use crate::workspace::session::{
    RestoredWorkspace, WorkspaceSession, assembly_draft_kind, decode_workspace, encode_workspace,
};
use crate::workspace::shell_state::ShellState;
use leptos::prelude::*;
use wasm_bindgen_futures::spawn_local;
#[derive(Clone, Copy)]
pub(crate) struct WorkspaceAutosave {
    session: WorkspaceSession,
    shell: ShellState,
    pub(crate) db: StoredValue<Option<crate::platform::storage::Database>, LocalStorage>,
    pub(crate) native_library_context: RwSignal<Option<collagene_library_protocol::LabContext>>,
    pub(crate) loaded: RwSignal<bool>,
    pub(crate) save_status: RwSignal<String>,
    pub(crate) save_generation: RwSignal<u64>,
}
impl WorkspaceAutosave {
    pub(crate) fn new(session: WorkspaceSession, shell: ShellState) -> Self {
        let db = StoredValue::new_local(None::<crate::platform::storage::Database>);
        let native_library_context = RwSignal::new(None::<collagene_library_protocol::LabContext>);
        let loaded = RwSignal::new(false);
        let save_status = RwSignal::new("Opening local storage…".to_string());
        let save_generation = RwSignal::new(0u64);
        Self {
            session,
            shell,
            db,
            native_library_context,
            loaded,
            save_status,
            save_generation,
        }
    }
    pub(crate) fn observe(
        &self,
        assembly_draft: RwSignal<crate::assembly::draft::model::Draft>,
        pair_draft: RwSignal<crate::primers::primer_pair_model::PairDraft>,
        stored_candidates: RwSignal<Vec<collagene_core::documents::transfer::Candidate>>,
        candidate: RwSignal<Option<collagene_core::documents::transfer::Candidate>>,
        apply_restored: Callback<RestoredWorkspace>,
    ) {
        let WorkspaceSession {
            documents,
            active,
            workspace_restoring,
            workspace_name,
            attachments,
            ..
        } = self.session;
        let ShellState { message, .. } = self.shell;
        let Self {
            db,
            native_library_context,
            loaded,
            save_status,
            save_generation,
            ..
        } = *self;
        Effect::new(move |_| {
            let draft = assembly_draft.get();
            if workspace_restoring.get() {
                return;
            }
            match crate::assembly::draft::model::encode_raw(&draft) {
                Ok(bytes) => {
                    let draft_kind = crate::assembly::draft::model::ATTACHMENT_KIND;
                    let hash = blake3::hash(&bytes).to_hex().to_string();
                    if attachments.with_untracked(|assets| {
                        assets.get(&hash).is_some_and(|a| a.kind == draft_kind)
                            && assets
                                .values()
                                .filter(|a| assembly_draft_kind(&a.kind))
                                .count()
                                == 1
                    }) {
                        return;
                    }
                    attachments.update(|assets| {
                        assets.retain(|_, a| !assembly_draft_kind(&a.kind));
                        assets.insert(
                            hash,
                            collagene_core::documents::workspace::WorkspaceAttachment {
                                kind: draft_kind.into(),
                                name: "assembly-draft.json".into(),
                                media_type: "application/json".into(),
                                bytes,
                            },
                        );
                    });
                }
                Err(error) => {
                    message.set(format!("Assembly draft could not be saved: {error}"));
                }
            }
        });
        Effect::new(move |_| {
            let draft = pair_draft.get();
            if workspace_restoring.get() {
                return;
            }
            match crate::primers::primer_pair_storage::encode(&draft) {
                Ok(bytes) => {
                    let kind = crate::primers::primer_pair_storage::ATTACHMENT_KIND;
                    let hash = blake3::hash(&bytes).to_hex().to_string();
                    if attachments
                        .with_untracked(|assets| assets.get(&hash).is_some_and(|a| a.kind == kind))
                    {
                        return;
                    }
                    let is_default =
                        draft == crate::primers::primer_pair_model::PairDraft::default();
                    if is_default
                        && !attachments
                            .with_untracked(|assets| assets.values().any(|a| a.kind == kind))
                    {
                        return;
                    }
                    attachments.update(|assets| {
                        assets.retain(|_, a| a.kind != kind);
                        if !is_default {
                            assets.insert(
                                hash,
                                collagene_core::documents::workspace::WorkspaceAttachment {
                                    kind: kind.into(),
                                    name: "primer-pair-draft.json".into(),
                                    media_type: "application/json".into(),
                                    bytes,
                                },
                            );
                        }
                    });
                }
                Err(error) => message.set(format!("Primer pair draft could not be saved: {error}")),
            }
        });
        spawn_local(async move {
            match storage::open().await {
                Ok(database) => {
                    native_library_context.set(database.native_context());
                    db.set_value(Some(database.clone()));
                    match storage::load(&database).await {
                        Ok(Some(bytes)) => {
                            let frame = decode_workspace(bytes, 0).await;
                            match frame.and_then(|frame| RestoredWorkspace::from_frame(&frame)) {
                                Ok(state) => apply_restored.run(state),
                                Err(error) => {
                                    message.set(format!("Autosave could not be read. The stored record has been left intact; open a saved workspace to recover. {error}"));
                                    save_status.set("Recovery needed".into());
                                    loaded.set(false);
                                    workspace_restoring.set(false);
                                    return;
                                }
                            }
                        }
                        Ok(None) => {}
                        Err(error) => {
                            message.set(error);
                            save_status.set("Autosave unavailable".into());
                            workspace_restoring.set(false);
                            return;
                        }
                    }
                    loaded.set(true);
                    workspace_restoring.set(false);
                }
                Err(error) => {
                    message.set(format!(
                        "{error}. Export a workspace file to keep your work."
                    ));
                    save_status.set("Autosave unavailable".into());
                    workspace_restoring.set(false);
                }
            }
        });
        Effect::new(move |_| {
            if !loaded.get() || workspace_restoring.get() {
                return;
            }
            let index = active.get();
            let name = workspace_name.get();
            let frame = documents.with(|docs| {
                stored_candidates.with(|plans| {
                    attachments.with(|assets| {
                        candidate.with(|current| {
                            WorkspaceSession::snapshot(
                                docs,
                                index,
                                current.as_ref(),
                                plans,
                                assets,
                                &name,
                            )
                        })
                    })
                })
            });
            let generation = save_generation.get_untracked().wrapping_add(1);
            save_generation.set(generation);
            save_status.set("Saving locally…".into());
            let frame = match frame {
                Ok(frame) => frame,
                Err(error) => {
                    save_status.set("Autosave failed".into());
                    message.set(error);
                    return;
                }
            };
            if let Some(database) = db.get_value() {
                spawn_local(async move {
                    let result = encode_workspace(frame, generation as u32).await;
                    // Work can finish out of order. Only the current generation may start a
                    // write transaction, including after native admission waits.
                    if generation != save_generation.get_untracked()
                        || workspace_restoring.get_untracked()
                    {
                        return;
                    }
                    let result = match result {
                        Ok(bytes) => {
                            storage::save(&database, &bytes, move || {
                                if generation != save_generation.get_untracked()
                                    || workspace_restoring.get_untracked()
                                {
                                    Err("Autosave superseded".into())
                                } else {
                                    Ok(())
                                }
                            })
                            .await
                        }
                        Err(error) => Err(error),
                    };
                    if generation != save_generation.get_untracked()
                        || workspace_restoring.get_untracked()
                    {
                        return;
                    }
                    match result {
                        Ok(()) => save_status.set(
                            if crate::platform::desktop_rpc::is_desktop() {
                                "Saved locally"
                            } else {
                                "Saved in this browser"
                            }
                            .into(),
                        ),
                        Err(error) => {
                            save_status.set("Autosave failed".into());
                            message.set(error);
                        }
                    }
                });
            }
        });
    }
}
