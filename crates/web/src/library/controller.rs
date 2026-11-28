use crate::platform::storage;

use collagene_core::documents::transfer;
use collagene_core::verification;
use collagene_core::{
    documents::workspace, material::lineage as run_lineage,
    sequence::annotations::library as feature_library,
};

use collagene_model::Feature;
use leptos::prelude::*;

use wasm_bindgen_futures::spawn_local;

use crate::workspace::session::RestoredWorkspace;
use crate::workspace::session::WorkspaceSession;

use crate::workspace::session::encode_workspace;

use crate::viewer::sequence::selection::ViewerSelection;
use crate::workspace::shell_state::ShellState;

#[derive(Clone, Copy)]
pub(crate) struct LibraryController {
    pub(crate) session: WorkspaceSession,
    pub(crate) viewer: ViewerSelection,
    pub(crate) shell: ShellState,
}
impl LibraryController {
    pub(crate) fn catalog_pending_changed(
        &self,
        native_catalog_epoch: RwSignal<u64>,
        native_catalog_pending: RwSignal<bool>,
        native_catalog_refresh: RwSignal<u64>,
        native_destination: RwSignal<Option<collagene_library_protocol::FolderId>>,
        native_library_context: RwSignal<Option<collagene_library_protocol::LabContext>>,
    ) -> Callback<bool, ()> {
        let session = self.session;
        Callback::new(move |pending: bool| {
            let token = native_catalog_epoch.get_untracked().wrapping_add(1);
            native_catalog_epoch.set(token);
            if pending {
                native_catalog_pending.set(true);
                crate::library::library_client::browser().disconnect();
                native_library_context.set(None);
                native_destination.set(None);
            } else if native_library_context.get_untracked().is_some() {
                native_catalog_pending.set(false);
                native_catalog_refresh
                    .update(|generation| *generation = generation.wrapping_add(1));
            } else {
                // Rebind the actual native catalog after an unknown/failed adoption,
                // never revive a cached browser context just from a UI status flag.
                native_catalog_pending.set(true);
                spawn_local(async move {
                    let bootstrap = crate::platform::desktop_rpc::call(
                        collagene_desktop_protocol::Command::Bootstrap,
                    )
                    .await;
                    if native_catalog_epoch.try_get_untracked() != Some(token) {
                        return;
                    }
                    if let Ok(collagene_desktop_protocol::Success::Bootstrap(value)) = bootstrap {
                        if matches!(
                            value.catalog,
                            collagene_desktop_protocol::CatalogStatus::Ready { .. }
                        ) {
                            if let Some(context) = value.library_context {
                                native_library_context.set(Some(context));
                                native_catalog_pending.set(false);
                                native_catalog_refresh
                                    .update(|generation| *generation = generation.wrapping_add(1));
                            }
                        }
                    }
                });
            }
        })
    }
    pub(crate) fn library_record_binding(
        &self,
        library_bindings: Memo<crate::library::library_bindings::Bindings>,
    ) -> Callback<(u64, crate::library::library_bindings::Binding), Result<(), String>> {
        let session = self.session;
        let WorkspaceSession {
            documents,
            workspace_restoring,
            workspace_epoch,
            attachments,
            ..
        } = self.session;
        Callback::new(
            move |(epoch, binding): (u64, crate::library::library_bindings::Binding)| {
                if workspace_epoch.get_untracked() != epoch || workspace_restoring.get_untracked() {
                    return Err("The server saved its revision, but the local workspace changed; no local link was overwritten.".to_string());
                }
                let retained = documents.with_untracked(|docs| {
                    docs.iter()
                        .find(|document| document.id == binding.local_document)
                        .is_some_and(|document| {
                            document.summaries().iter().any(|event| {
                                event.revision == binding.local_revision.as_str()
                                    && binding.local_head_event.as_deref()
                                        == Some(event.event_id.as_str())
                            })
                        })
                });
                if !retained {
                    return Err("The remote revision is saved, but its captured local history is no longer present; reopen it from the lab library.".into());
                }
                let unchanged = attachments.with_untracked(|assets| {
                    crate::library::library_documents::bindings(assets).map(|links| {
                        links.get(&binding.remote, &binding.local_document) == Some(&binding)
                    })
                })?;
                if unchanged {
                    return Ok(());
                }
                attachments
                    .try_update(|assets| {
                        crate::library::library_documents::record_binding(assets, binding)
                    })
                    .ok_or_else(|| {
                        "Workspace closed before its remote link could be recorded.".to_string()
                    })?
            },
        )
    }
    pub(crate) fn library_open(
        &self,
        apply_restored: Callback<RestoredWorkspace>,
        assembly_visible: RwSignal<bool>,
        db: StoredValue<Option<crate::platform::storage::Database>, leptos::prelude::LocalStorage>,
        library_bindings: Memo<crate::library::library_bindings::Bindings>,
        library_snapshot: Callback<
            (),
            Result<collagene_core::documents::workspace::Workspace, String>,
        >,
        loaded: RwSignal<bool>,
        save_generation: RwSignal<u64>,
    ) -> Callback<
        crate::library::library_view::Opened,
        crate::library::library_view::LocalResult<crate::library::library_view::LocalDocument>,
    > {
        let session = self.session;
        let WorkspaceSession {
            documents,
            active,
            workspace_restoring,
            workspace_epoch,
            attachments,
            ..
        } = self.session;
        let ViewerSelection { molecule, .. } = self.viewer;
        Callback::new(
            move |opened: crate::library::library_view::Opened| -> crate::library::library_view::LocalResult<
                crate::library::library_view::LocalDocument,
            > {
                let prepared = (|| {
                    if workspace_restoring.get_untracked()
                        || !loaded.get_untracked()
                        || !opened.guard.is_current()
                    {
                        return Err("Wait for the current workspace operation to finish, then open the library revision again.".to_string());
                    }
                    let database=db.get_value().ok_or("Local browser storage is unavailable; the library revision was not imported.")?;
                    let backup = library_snapshot.run(())?;
                    let fingerprint =
                        collagene_core::material::lineage::workspace_fingerprint(&backup)
                            .map_err(|error| error.to_string())?;
                    let epoch = workspace_epoch.get_untracked();
                    workspace_restoring.set(true);
                    save_generation.update(|generation| *generation = generation.wrapping_add(1));
                    Ok((database, backup, fingerprint, epoch))
                })();
                Box::pin(async move {
                    let (database, backup, fingerprint, epoch) = prepared?;
                    let result=async {
                let decoded=crate::library::library_documents::decode(&opened).await?;
                let guard=||->Result<(),String>{
                    if !opened.guard.is_current() || workspace_epoch.get_untracked()!=epoch
                        || collagene_core::material::lineage::workspace_fingerprint(&library_snapshot.run(())?).map_err(|error|error.to_string())?!=fingerprint
                    {return Err("Local workspace or library session changed; the import was not committed.".into());}
                    Ok(())
                };
                guard()?;
                let mut next=crate::library::library_documents::merge(&backup,&decoded,&opened)?;
                let mut state=RestoredWorkspace::from_frame(&next)?;
                let document=state.documents.get(state.active).ok_or("Imported document is unavailable")?;
                let mut info=crate::library::library_view::LocalDocument{
                    id:document.id.clone(),name:document.molecule.name.clone(),
                    revision:collagene_library_protocol::Digest::new(document.head_revision()).map_err(|error|error.to_string())?,
                    head_event:Some(document.head_event().into()),workspace_epoch:epoch.wrapping_add(1),
                };
                crate::library::library_documents::record_binding(&mut next.attachments,crate::library::library_bindings::Binding{
                    remote:opened.remote.clone(),local_document:info.id.clone(),local_revision:info.revision.clone(),
                    local_head_event:info.head_event.clone(),remote_document:opened.document.id.clone(),remote_revision:opened.revision.id.clone(),
                })?;
                state.attachments=next.attachments.clone();
                let current=encode_workspace(next,u32::MAX-303).await?;
                let previous=encode_workspace(backup,u32::MAX-304).await?;
                guard()?;
                crate::platform::storage::replace_workspace(&database,&current,&previous).await?;
                // A committed storage transaction reaches App even if its sidebar
                // has closed; a late view cancellation must not undo durable work.
                let visible=assembly_visible.get_untracked();
                batch(||{apply_restored.run(state);assembly_visible.set(visible);});
                info.workspace_epoch=workspace_epoch.get_untracked();
                Ok(info)
            }.await;
                    workspace_restoring.set(false);
                    result
                })
            },
        )
    }
    pub(crate) fn keep_library(
        &self,
        candidate: RwSignal<Option<transfer::Candidate>>,
        libraries: RwSignal<Vec<crate::viewer::sequence::feature_tools::LibraryRevision>>,
        stored_candidates: RwSignal<Vec<transfer::Candidate>>,
    ) -> Callback<collagene_core::sequence::annotations::library::FeatureLibrary, Result<(), String>>
    {
        let session = self.session;
        let WorkspaceSession {
            documents,
            active,
            workspace_restoring,
            workspace_name,
            attachments,
            ..
        } = self.session;
        Callback::new(
            move |library: collagene_core::sequence::annotations::library::FeatureLibrary| {
                if workspace_restoring.get_untracked() {
                    return Err("Wait for workspace restoration to finish.".into());
                }
                library.validate().map_err(|error| error.to_string())?;
                let revision = library.revision();
                if libraries.with_untracked(|items| {
                    items.iter().any(|existing| existing.revision == revision)
                }) {
                    return Ok(());
                }
                if libraries.with_untracked(Vec::len) >= 128 {
                    return Err("This web workspace retains at most 128 immutable library revisions. Export libraries before opening a smaller workspace.".into());
                }
                let bytes = collagene_core::sequence::annotations::library::save_library(&library)
                    .map_err(|error| error.to_string())?;
                let hash = blake3::hash(&bytes).to_hex().to_string();
                let mut assets = attachments.get_untracked();
                assets.insert(
                    hash,
                    collagene_core::documents::workspace::WorkspaceAttachment {
                        kind: "feature-library-v1".into(),
                        name: format!("Feature library {}", &revision[..12]),
                        media_type: "application/json".into(),
                        bytes,
                    },
                );
                let frame = documents.with_untracked(|docs| {
                    WorkspaceSession::snapshot(
                        docs,
                        active.get_untracked(),
                        candidate.get_untracked().as_ref(),
                        &stored_candidates.get_untracked(),
                        &assets,
                        &workspace_name.get_untracked(),
                    )
                })?;
                frame.validate().map_err(|error| error.to_string())?;
                batch(|| {
                    attachments.set(assets);
                    libraries.update(|items| {
                        items.push(crate::viewer::sequence::feature_tools::LibraryRevision {
                            revision,
                            library: std::sync::Arc::new(library),
                        })
                    });
                });
                Ok(())
            },
        )
    }
}
