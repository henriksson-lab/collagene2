//! Explicit reviewed metadata publication. No automatic parents or reconciliation.
use crate::{Error, Pending, RemoteConfig, Repository, Result, Target, durable, store, sync, wire};
use serde::{Deserialize, Serialize};
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub enum PublicationTarget {
    Folder(store::FolderId),
    DocumentMetadata(store::DocumentId),
    DeletedFolder(store::FolderId),
}
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub enum MetadataAction {
    CreateFolder {
        parent: Option<wire::FolderId>,
        name: String,
    },
    RenameFolder {
        id: wire::FolderId,
        expected_generation: u64,
        name: String,
    },
    MoveFolder {
        id: wire::FolderId,
        expected_generation: u64,
        parent: Option<wire::FolderId>,
    },
    RenameDocument {
        id: wire::DocumentId,
        expected_generation: u64,
        name: String,
    },
    MoveDocument {
        id: wire::DocumentId,
        expected_generation: u64,
        folder: Option<wire::FolderId>,
    },
    DeleteEmptyFolder {
        id: wire::FolderId,
        expected_generation: u64,
    },
}
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct PublicationDescription {
    pub local_lab: wire::LabId,
    pub target: PublicationTarget,
    pub remote: RemoteConfig,
    pub expected_checkpoint: Option<store::Digest>,
    pub local_generation: Option<u64>,
    pub local_name: Option<String>,
    pub local_parent: Option<store::FolderId>,
    pub remote_parent: Option<wire::FolderId>,
    pub parent_binding: Option<wire::Folder>,
    pub action: MetadataAction,
}
/// A retained remote binding whose local folder has actually been deleted.
/// Display names come from its acknowledged remote snapshot, not raw ID input.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct DeletedFolderPublication {
    pub local_folder: store::FolderId,
    pub remote_folder: wire::Folder,
}
#[derive(Debug, Clone)]
pub struct PublicationReview {
    description: PublicationDescription,
    digest: wire::Digest,
}
impl PublicationReview {
    pub fn description(&self) -> &PublicationDescription {
        &self.description
    }
    pub fn digest(&self) -> wire::Digest {
        self.digest.clone()
    }
}
impl MetadataAction {
    fn operation(&self) -> Result<sync::MutationOperation> {
        use sync::MutationOperation as O;
        let name = |s: &String| -> Result<wire::Name> {
            if s.len() > 1024 {
                return Err(Error::Invalid("metadata name bounds"));
            }
            let n = wire::Name::new(s).map_err(|_| Error::Invalid("metadata name"))?;
            wire::validate_name(&n).map_err(|_| Error::Invalid("metadata name"))?;
            Ok(n)
        };
        let op = match self {
            Self::CreateFolder { parent, name: n } => O::CreateFolder {
                parent: parent.clone(),
                name: name(n)?,
            },
            Self::RenameFolder {
                id,
                expected_generation,
                name: n,
            } => O::RenameFolder {
                id: id.clone(),
                expected_generation: *expected_generation,
                name: name(n)?,
            },
            Self::MoveFolder {
                id,
                expected_generation,
                parent,
            } => O::MoveFolder {
                id: id.clone(),
                expected_generation: *expected_generation,
                parent: parent.clone(),
            },
            Self::RenameDocument {
                id,
                expected_generation,
                name: n,
            } => O::RenameDocument {
                id: id.clone(),
                expected_generation: *expected_generation,
                name: name(n)?,
            },
            Self::MoveDocument {
                id,
                expected_generation,
                folder,
            } => O::MoveDocument {
                id: id.clone(),
                expected_generation: *expected_generation,
                folder: folder.clone(),
            },
            Self::DeleteEmptyFolder {
                id,
                expected_generation,
            } => O::DeleteEmptyFolder {
                id: id.clone(),
                expected_generation: *expected_generation,
            },
        };
        op.validate()
            .map_err(|_| Error::Invalid("metadata action"))?;
        Ok(op)
    }
}
impl PublicationDescription {
    /// Digest binds the complete reviewed context, not authority supplied by a client.
    pub fn digest(&self) -> Result<wire::Digest> {
        self.remote.validate()?;
        if let Some(n) = &self.local_name {
            let n = wire::Name::new(n).map_err(|_| Error::Invalid("local name"))?;
            wire::validate_name(&n).map_err(|_| Error::Invalid("local name"))?;
        }
        if self.local_generation == Some(0) {
            return Err(Error::Invalid("local generation"));
        }
        if let Some(b) = &self.parent_binding {
            b.validate().map_err(|_| Error::Invalid("parent binding"))?;
        }
        self.action.operation()?;
        let bytes = serde_json::to_vec(&("collagene-reviewed-metadata-publication-1", self))
            .map_err(|_| Error::Invalid("metadata review"))?;
        if bytes.len() > 16 * 1024 {
            return Err(Error::Invalid("metadata review bounds"));
        }
        Ok(wire::Digest::of(&bytes))
    }
    fn target(&self) -> Target {
        match &self.target {
            PublicationTarget::Folder(id) => Target::Folder(id.as_str().into()),
            PublicationTarget::DocumentMetadata(id) => Target::Document(id.as_str().into()),
            PublicationTarget::DeletedFolder(id) => Target::DeleteFolder(id.as_str().into()),
        }
    }
}
impl Repository {
    /// Read-only complete local binding discovery, in stable local-ID order.
    /// Hosts page it with a repository checkpoint pin. Above 10,000 deleted
    /// bindings fails explicitly instead of silently hiding later entries.
    pub fn deleted_folder_publications(
        &mut self,
        remote: &str,
    ) -> Result<Vec<DeletedFolderPublication>> {
        self.reload()?;
        if !self.state.remotes.contains_key(remote) {
            return Err(Error::Invalid("unknown remote"));
        }
        let mut result = Vec::new();
        for (id, bindings) in &self.state.folders {
            let Some(binding) = bindings.get(remote) else {
                continue;
            };
            let id = store::FolderId::parse(id.clone())?;
            match self.db.folder(&self.access, &id) {
                Ok(_) => continue,
                Err(store::Error::NotFound) => {}
                Err(error) => return Err(error.into()),
            }
            result.push(DeletedFolderPublication {
                local_folder: id,
                remote_folder: binding.remote.clone(),
            });
            if result.len() > 10_000 {
                return Err(Error::Invalid(
                    "deleted-folder discovery exceeds 10000 bindings",
                ));
            }
        }
        Ok(result)
    }
    fn reviewed_parent(
        &mut self,
        target: &PublicationTarget,
        parent: Option<&store::FolderId>,
        remote: &str,
    ) -> Result<(Option<wire::FolderId>, Option<wire::Folder>)> {
        // A clone destination is a local mount, not the remote parent. Match
        // actual retained job membership before ordinary parent bindings, even
        // if that mount itself happens to be bound to another remote folder.
        let mut mounted: Option<Option<wire::FolderId>> = None;
        for job in self
            .tree_clone_jobs()
            .map_err(|_| Error::Invalid("clone mount metadata"))?
        {
            if job.remote.id != remote || parent != job.destination.as_ref() {
                continue;
            }
            let interpreted = match target {
                PublicationTarget::Folder(id) if job.folders.values().any(|f| f == id) => {
                    match &job.root {
                        crate::folders::TreeRoot::LabRoot => Some(None),
                        crate::folders::TreeRoot::Folder { id: root, .. }
                            if job.folders.get(root.as_str()) == Some(id) =>
                        {
                            Some(
                                self.folder_binding(id, remote)
                                    .ok_or(Error::Invalid("unbound cloned root"))?
                                    .remote
                                    .parent
                                    .clone(),
                            )
                        }
                        _ => None,
                    }
                }
                PublicationTarget::DocumentMetadata(id)
                    if matches!(job.root, crate::folders::TreeRoot::LabRoot)
                        && job.documents.values().any(|d| d == id) =>
                {
                    Some(None)
                }
                _ => None,
            };
            if let Some(interpretation) = interpreted {
                if mounted.as_ref().is_some_and(|old| old != &interpretation) {
                    return Err(Error::Invalid("ambiguous cloned mount parent"));
                }
                mounted = Some(interpretation);
            }
        }
        if let Some(remote_parent) = mounted {
            return Ok((remote_parent, None));
        }
        match parent {
            None => Ok((None, None)),
            Some(id) => {
                let b = self
                    .folder_binding(id, remote)
                    .ok_or(Error::UnboundFolder)?;
                Ok((Some(b.remote.id.clone()), Some(b.remote.clone())))
            }
        }
    }
    pub fn prepare_metadata_publication(
        &mut self,
        target: &PublicationTarget,
        remote_id: &str,
    ) -> Result<PublicationReview> {
        self.reload()?;
        let remote = self
            .state
            .remotes
            .get(remote_id)
            .cloned()
            .ok_or(Error::Invalid("unknown remote"))?;
        let (generation, name, parent, remote_parent, parent_binding, action) = match target {
            PublicationTarget::Folder(id) => {
                let local = self.db.folder(&self.access, id)?;
                let (parent, binding) =
                    self.reviewed_parent(target, local.parent.as_ref(), remote_id)?;
                let action = match self.folder_binding(id, remote_id) {
                    None => MetadataAction::CreateFolder {
                        parent: parent.clone(),
                        name: local.name.clone(),
                    },
                    Some(b) if local.name != b.remote.name.as_str() => {
                        MetadataAction::RenameFolder {
                            id: b.remote.id.clone(),
                            expected_generation: b.remote.generation,
                            name: local.name.clone(),
                        }
                    }
                    Some(b) if parent != b.remote.parent => MetadataAction::MoveFolder {
                        id: b.remote.id.clone(),
                        expected_generation: b.remote.generation,
                        parent: parent.clone(),
                    },
                    Some(_) => return Err(Error::NoChanges),
                };
                (
                    Some(local.generation),
                    Some(local.name),
                    local.parent,
                    parent,
                    binding,
                    action,
                )
            }
            PublicationTarget::DocumentMetadata(id) => {
                let local = self.db.document(&self.access, id)?;
                let b = self
                    .document_binding(id, remote_id)
                    .cloned()
                    .ok_or(Error::Invalid("unbound document"))?;
                let (parent, binding) =
                    self.reviewed_parent(target, local.folder.as_ref(), remote_id)?;
                let action = if local.name != b.remote.name.as_str() {
                    MetadataAction::RenameDocument {
                        id: b.remote.id.clone(),
                        expected_generation: b.remote.generation,
                        name: local.name.clone(),
                    }
                } else if parent != b.remote.folder {
                    MetadataAction::MoveDocument {
                        id: b.remote.id.clone(),
                        expected_generation: b.remote.generation,
                        folder: parent.clone(),
                    }
                } else {
                    return Err(Error::NoChanges);
                };
                (
                    Some(local.generation),
                    Some(local.name),
                    local.folder,
                    parent,
                    binding,
                    action,
                )
            }
            PublicationTarget::DeletedFolder(id) => {
                match self.db.folder(&self.access, id) {
                    Err(store::Error::NotFound) => {}
                    Ok(_) => {
                        return Err(Error::Invalid("delete local empty folder explicitly first"));
                    }
                    Err(e) => return Err(e.into()),
                }
                let b = self
                    .folder_binding(id, remote_id)
                    .ok_or(Error::Invalid("unbound deleted folder"))?;
                (
                    None,
                    Some(b.remote.name.as_str().into()),
                    None,
                    b.remote.parent.clone(),
                    None,
                    MetadataAction::DeleteEmptyFolder {
                        id: b.remote.id.clone(),
                        expected_generation: b.remote.generation,
                    },
                )
            }
        };
        let description = PublicationDescription {
            local_lab: self.local_context()?.lab_id,
            target: target.clone(),
            remote,
            expected_checkpoint: self.checkpoint.clone(),
            local_generation: generation,
            local_name: name,
            local_parent: parent,
            remote_parent,
            parent_binding,
            action,
        };
        let digest = description.digest()?;
        Ok(PublicationReview {
            description,
            digest,
        })
    }
    /// Validate a bounded existing outbox frame before comparing its immutable
    /// mutation intent. This is read-only replay; no fresh mutation ID is minted.
    fn reviewed_pending_operation(&mut self, p: &Pending) -> Result<sync::MutationOperation> {
        if p.chunks == 0 || p.chunks > 1 {
            return Err(Error::Invalid("outbox chunks"));
        }
        let mut bytes = Vec::new();
        for i in 0..p.chunks {
            let c = self
                .db
                .read_checkpoint(
                    &self.access,
                    &durable::CheckpointKey::repository(format!("outbox-{}-{i}", p.id.as_str())),
                )?
                .ok_or(Error::Invalid("missing outbox chunk"))?;
            if c.value.version != 1
                || bytes
                    .len()
                    .checked_add(c.value.bytes.len())
                    .is_none_or(|n| n > 16 * 1024)
            {
                return Err(Error::Invalid("outbox frame bounds"));
            }
            bytes.extend_from_slice(&c.value.bytes);
        }
        if wire::Digest::of(&bytes) != p.frame_hash {
            return Err(Error::Invalid("outbox integrity"));
        }
        let r = sync::decode_request(&bytes).map_err(|_| Error::Invalid("outbox framing"))?;
        let sync::Request::Apply { mutation } = r.request else {
            return Err(Error::Invalid("outbox request"));
        };
        if r.expected_lab != p.remote.lab || mutation.id != p.id {
            return Err(Error::Invalid("outbox scope"));
        }
        Ok(mutation.operation)
    }
    pub fn queue_metadata_publication(&mut self, review: &PublicationReview) -> Result<String> {
        self.queue_metadata_publication_checked(review.description(), &review.digest())
    }
    pub fn queue_metadata_publication_checked(
        &mut self,
        description: &PublicationDescription,
        expected: &wire::Digest,
    ) -> Result<String> {
        if description.digest()? != *expected {
            return Err(Error::Stale);
        }
        self.reload()?;
        if self.local_context()?.lab_id != description.local_lab {
            return Err(Error::Stale);
        }
        let operation = description.action.operation()?;
        let target = description.target();
        let remote = self
            .state
            .remotes
            .get(&description.remote.id)
            .ok_or(Error::Invalid("unknown remote"))?;
        if serde_json::to_value(remote).map_err(|_| Error::Invalid("remote"))?
            != serde_json::to_value(&description.remote).map_err(|_| Error::Invalid("remote"))?
        {
            return Err(Error::Stale);
        }
        // Resume exactly the existing queued intent, even if its creation advanced the
        // repository checkpoint or the caller lost the reply and restarted.
        let existing = self
            .state
            .pending
            .values()
            .find(|p| {
                let same_target = match (&p.target, &target) {
                    (Target::Document(a), Target::Document(b)) => a == b,
                    (
                        Target::Folder(a) | Target::DeleteFolder(a),
                        Target::Folder(b) | Target::DeleteFolder(b),
                    ) => a == b,
                    _ => false,
                };
                p.remote.id == description.remote.id && same_target
            })
            .cloned();
        if let Some(p) = existing {
            if p.target == target
                && p.native.is_none()
                && p.local_revision.is_none()
                && self
                    .reviewed_pending_operation(&p)?
                    .intent_digest()
                    .map_err(|_| Error::Invalid("pending intent"))?
                    == operation
                        .intent_digest()
                        .map_err(|_| Error::Invalid("metadata intent"))?
            {
                return Ok(p.id.as_str().into());
            }
            return Err(Error::Pending);
        }
        let current =
            self.prepare_metadata_publication(&description.target, &description.remote.id)?;
        if current.digest() != *expected || self.checkpoint != description.expected_checkpoint {
            return Err(Error::Stale);
        }
        self.enqueue(description.remote.clone(), target, operation, None, None)
    }
}

#[cfg(test)]
mod discovery_tests {
    use super::*;
    use std::sync::atomic::{AtomicU64, Ordering};
    static NEXT: AtomicU64 = AtomicU64::new(1);
    #[test]
    fn discovery_does_not_hide_the_201st_binding_or_mutate_checkpoint() {
        let path = std::env::temp_dir().join(format!(
            "collagene-discovery-bound-{}-{}.sqlite",
            std::process::id(),
            NEXT.fetch_add(1, Ordering::Relaxed)
        ));
        let mut repository = Repository::open(&path).unwrap();
        let remote = repository
            .add_remote(
                "Lab",
                "https://lab.example",
                wire::LabId::new(format!("lab-{:032x}", 1)).unwrap(),
            )
            .unwrap();
        let mut next = repository.state.clone();
        for n in 1..=201 {
            let local = store::FolderId::parse(format!("fld-{n:032x}")).unwrap();
            let binding = crate::FolderBinding {
                remote: wire::Folder {
                    id: wire::FolderId::new(format!("fld-{n:032x}")).unwrap(),
                    parent: None,
                    name: wire::Name::new(format!("Deleted folder {n}")).unwrap(),
                    generation: 1,
                },
            };
            next.folders
                .entry(local.as_str().into())
                .or_default()
                .insert(remote.id.clone(), binding);
        }
        repository.commit(next, vec![]).unwrap();
        let checkpoint = repository.checkpoint_version();
        let rows = repository.deleted_folder_publications(&remote.id).unwrap();
        assert_eq!(rows.len(), 201);
        assert_eq!(rows[200].local_folder.as_str(), format!("fld-{:032x}", 201));
        assert_eq!(repository.checkpoint_version(), checkpoint);
        drop(repository);
        for p in [
            path.clone(),
            std::path::PathBuf::from(format!("{}-wal", path.display())),
            std::path::PathBuf::from(format!("{}-shm", path.display())),
        ] {
            let _ = std::fs::remove_file(p);
        }
    }
}
