use collagene_core::{Molecule, Topology, protocol};
use collagene_lab_library::{self as store, mutations as durable};
use collagene_library_protocol as wire;
use collagene_library_sync as sync;
use collagene_repository::{
    self as repo, Channel, Error, RemoteConfig, RemoteTransport, Repository, SyncResult,
    TransportError,
    folder_publication::{MetadataAction, PublicationTarget as Target},
};
use std::{
    path::PathBuf,
    sync::atomic::{AtomicU64, Ordering},
};
static NEXT: AtomicU64 = AtomicU64::new(1);
struct Disk(PathBuf);
impl Disk {
    fn new() -> Self {
        Self(std::env::temp_dir().join(format!(
            "collagene-reviewed-metadata-{}-{}.sqlite",
            std::process::id(),
            NEXT.fetch_add(1, Ordering::Relaxed)
        )))
    }
}
impl Drop for Disk {
    fn drop(&mut self) {
        for path in [
            self.0.clone(),
            PathBuf::from(format!("{}-wal", self.0.display())),
            PathBuf::from(format!("{}-shm", self.0.display())),
        ] {
            let _ = std::fs::remove_file(path);
        }
    }
}
struct Remote {
    db: store::Library,
    access: store::LabAccess,
    lose: bool,
    applies: usize,
}
impl Remote {
    fn new() -> Self {
        let mut db = store::Library::open_in_memory().unwrap();
        let lab = db.create_lab("Shared lab").unwrap();
        let access = db.trusted_host_access(&lab.id).unwrap();
        Self {
            db,
            access,
            lose: false,
            applies: 0,
        }
    }
    fn add(&self, r: &mut Repository) -> RemoteConfig {
        r.add_remote(
            "Lab",
            "https://lab.example",
            wire::LabId::new(self.access.lab_id().as_str()).unwrap(),
        )
        .unwrap()
    }
}
impl RemoteTransport for Remote {
    fn exchange(
        &mut self,
        _: &RemoteConfig,
        c: Channel,
        b: &[u8],
    ) -> Result<Vec<u8>, TransportError> {
        match c {
            Channel::Library => {
                collagene_library_host::dispatch_frame(&mut self.db, &self.access, b, None)
                    .map_err(|_| TransportError::Unavailable)
            }
            Channel::Sync => {
                let q = sync::decode_request(b).map_err(|_| TransportError::Unavailable)?;
                let apply = matches!(q.request, sync::Request::Apply { .. });
                if apply {
                    self.applies += 1;
                }
                let reply = sync::native::dispatch(&mut self.db, &self.access, &q)
                    .map_err(|_| TransportError::Unavailable)?;
                if apply && self.lose {
                    self.lose = false;
                    return Err(TransportError::Unavailable);
                }
                sync::encode_response(&reply).map_err(|_| TransportError::Unavailable)
            }
        }
    }
}
fn folder(r: &mut Repository, parent: Option<store::FolderId>, name: &str) -> store::Folder {
    match r
        .mutate_local(durable::MutationOperation::CreateFolder {
            parent,
            name: name.into(),
        })
        .unwrap()
    {
        durable::Applied::CreatedFolder(f) => f,
        other => panic!("{other:?}"),
    }
}
fn publish_folder(
    r: &mut Repository,
    remote: &mut Remote,
    c: &RemoteConfig,
    f: &store::Folder,
) -> String {
    let review = r
        .prepare_metadata_publication(&Target::Folder(f.id.clone()), &c.id)
        .unwrap();
    let id = r.queue_metadata_publication(&review).unwrap();
    assert!(matches!(
        r.resume_push(&id, remote).unwrap(),
        SyncResult::Applied(_)
    ));
    id
}
#[test]
fn explicitly_reviewed_root_folder_creation_retries_same_nonce_across_lost_reply_and_restart() {
    let disk = Disk::new();
    let mut r = Repository::open(&disk.0).unwrap();
    let mut remote = Remote::new();
    let config = remote.add(&mut r);
    let local = folder(&mut r, None, "New offline folder");
    let review = r
        .prepare_metadata_publication(&Target::Folder(local.id.clone()), &config.id)
        .unwrap();
    let description = review.description().clone();
    let digest = review.digest();
    assert_eq!(description.remote_parent, None);
    assert!(matches!(
        &description.action,
        MetadataAction::CreateFolder { parent: None, .. }
    ));
    let job = r.queue_metadata_publication(&review).unwrap();
    assert_eq!(r.queue_metadata_publication(&review).unwrap(), job);
    remote.lose = true;
    assert!(r.resume_push(&job, &mut remote).is_err());
    drop(r);
    let mut r = Repository::open(&disk.0).unwrap();
    assert_eq!(
        r.queue_metadata_publication_checked(&description, &digest)
            .unwrap(),
        job
    );
    assert!(matches!(
        r.resume_push(&job, &mut remote).unwrap(),
        SyncResult::Applied(_)
    ));
    assert_eq!(remote.applies, 1);
    let page = remote
        .db
        .browse(&remote.access, None, 0, 100, None)
        .unwrap();
    assert_eq!(page.entries.len(), 1);
    assert!(r.folder_binding(&local.id, &config.id).is_some());
}
#[test]
fn nested_parent_requires_explicit_publication_and_changed_local_review_never_queues() {
    let disk = Disk::new();
    let mut r = Repository::open(&disk.0).unwrap();
    let mut remote = Remote::new();
    let config = remote.add(&mut r);
    let parent = folder(&mut r, None, "Parent");
    let child = folder(&mut r, Some(parent.id.clone()), "Child");
    assert!(matches!(
        r.prepare_metadata_publication(&Target::Folder(child.id.clone()), &config.id),
        Err(Error::UnboundFolder)
    ));
    assert!(r.pending_ids().is_empty());
    publish_folder(&mut r, &mut remote, &config, &parent);
    let review = r
        .prepare_metadata_publication(&Target::Folder(child.id.clone()), &config.id)
        .unwrap();
    assert_eq!(
        review.description().remote_parent,
        r.folder_binding(&parent.id, &config.id)
            .map(|b| b.remote.id.clone())
    );
    let mut second = Repository::open(&disk.0).unwrap();
    second
        .mutate_local(durable::MutationOperation::RenameFolder {
            id: child.id.clone(),
            expected_generation: child.generation,
            name: "Changed elsewhere".into(),
        })
        .unwrap();
    assert!(matches!(
        r.queue_metadata_publication(&review),
        Err(Error::Stale)
    ));
    assert!(r.pending_ids().is_empty());
    let fresh = r
        .prepare_metadata_publication(&Target::Folder(child.id.clone()), &config.id)
        .unwrap();
    let mut forged = fresh.description().clone();
    forged.remote_parent = None;
    let hash = forged.digest().unwrap();
    assert!(matches!(
        r.queue_metadata_publication_checked(&forged, &hash),
        Err(Error::Stale)
    ));
    assert!(r.pending_ids().is_empty());
}
#[test]
fn document_rename_then_move_are_separately_reviewed_generation_guarded_operations() {
    let disk = Disk::new();
    let mut r = Repository::open(&disk.0).unwrap();
    let mut remote = Remote::new();
    let config = remote.add(&mut r);
    let first = folder(&mut r, None, "First");
    let second = folder(&mut r, None, "Second");
    publish_folder(&mut r, &mut remote, &config, &first);
    publish_folder(&mut r, &mut remote, &config, &second);
    let bytes = protocol::save_document(
        &Molecule::dna("Original plasmid", "ACGTACGT", Topology::Linear).unwrap(),
    )
    .unwrap();
    let payload = store::Payload {
        kind: store::PayloadKind::NativeDocument,
        archive_hash: store::Digest::parse(wire::Digest::of(&bytes).as_str()).unwrap(),
        bytes,
        selection: None,
        sharing: None,
    };
    let doc = match r
        .mutate_local(durable::MutationOperation::ImportDocument {
            folder: Some(first.id),
            name: "Original plasmid".into(),
            payload,
        })
        .unwrap()
    {
        durable::Applied::Imported(d) => d,
        other => panic!("{other:?}"),
    };
    let job = r.queue_push(&doc.id, &config.id).unwrap();
    assert!(matches!(
        r.resume_push(&job, &mut remote).unwrap(),
        SyncResult::Applied(_)
    ));
    let renamed = match r
        .mutate_local(durable::MutationOperation::RenameDocument {
            id: doc.id.clone(),
            expected_generation: doc.generation,
            name: "Renamed plasmid".into(),
        })
        .unwrap()
    {
        durable::Applied::RenamedDocument(d) => d,
        other => panic!("{other:?}"),
    };
    r.mutate_local(durable::MutationOperation::MoveDocument {
        id: doc.id.clone(),
        expected_generation: renamed.generation,
        folder: Some(second.id.clone()),
    })
    .unwrap();
    let review = r
        .prepare_metadata_publication(&Target::DocumentMetadata(doc.id.clone()), &config.id)
        .unwrap();
    assert!(matches!(
        review.description().action,
        MetadataAction::RenameDocument { .. }
    ));
    let job = r.queue_metadata_publication(&review).unwrap();
    assert!(matches!(
        r.resume_push(&job, &mut remote).unwrap(),
        SyncResult::Applied(_)
    ));
    let review = r
        .prepare_metadata_publication(&Target::DocumentMetadata(doc.id.clone()), &config.id)
        .unwrap();
    assert!(matches!(
        review.description().action,
        MetadataAction::MoveDocument { .. }
    ));
    let job = r.queue_metadata_publication(&review).unwrap();
    assert!(matches!(
        r.resume_push(&job, &mut remote).unwrap(),
        SyncResult::Applied(_)
    ));
    let binding = r.document_binding(&doc.id, &config.id).unwrap();
    assert_eq!(binding.remote.name.as_str(), "Renamed plasmid");
    assert_eq!(
        binding.remote.folder,
        r.folder_binding(&second.id, &config.id)
            .map(|b| b.remote.id.clone())
    );
    assert!(matches!(
        r.prepare_metadata_publication(&Target::DocumentMetadata(doc.id), &config.id),
        Err(Error::NoChanges)
    ));
}
#[test]
fn deleted_bound_folder_review_requires_local_absence_and_remote_nonempty_rejection_preserves_data()
{
    let disk = Disk::new();
    let mut r = Repository::open(&disk.0).unwrap();
    let mut remote = Remote::new();
    let config = remote.add(&mut r);
    let local = folder(&mut r, None, "Empty locally");
    publish_folder(&mut r, &mut remote, &config, &local);
    assert!(
        r.prepare_metadata_publication(&Target::DeletedFolder(local.id.clone()), &config.id)
            .is_err()
    );
    let bound = r
        .folder_binding(&local.id, &config.id)
        .unwrap()
        .remote
        .clone();
    r.mutate_local(durable::MutationOperation::DeleteEmptyFolder {
        id: local.id.clone(),
        expected_generation: local.generation,
    })
    .unwrap();
    let review = r
        .prepare_metadata_publication(&Target::DeletedFolder(local.id.clone()), &config.id)
        .unwrap();
    assert!(matches!(
        review.description().action,
        MetadataAction::DeleteEmptyFolder { .. }
    ));
    let parent = store::FolderId::parse(bound.id.as_str()).unwrap();
    let child = remote
        .db
        .create_folder(&remote.access, Some(&parent), "Other member's child")
        .unwrap();
    let job = r.queue_metadata_publication(&review).unwrap();
    assert!(matches!(
        r.resume_push(&job, &mut remote).unwrap(),
        SyncResult::Rejected(_)
    ));
    assert!(remote.db.folder(&remote.access, &child.id).is_ok());
    assert!(remote.db.folder(&remote.access, &parent).is_ok());
}

#[test]
fn deleted_folder_discovery_is_read_only_remote_scoped_and_requires_actual_local_absence() {
    let disk = Disk::new();
    let mut r = Repository::open(&disk.0).unwrap();
    let mut remote = Remote::new();
    let config = remote.add(&mut r);
    let other = r
        .add_remote(
            "Other independent database",
            "https://other.example",
            wire::LabId::new(format!("lab-{:032x}", 987)).unwrap(),
        )
        .unwrap();
    let present = folder(&mut r, None, "Still present locally");
    let deleted = folder(&mut r, None, "Explicitly deleted locally");
    let unbound = folder(&mut r, None, "Deleted but never published");
    publish_folder(&mut r, &mut remote, &config, &present);
    publish_folder(&mut r, &mut remote, &config, &deleted);
    assert!(
        r.deleted_folder_publications(&config.id)
            .unwrap()
            .is_empty()
    );
    for f in [&deleted, &unbound] {
        r.mutate_local(durable::MutationOperation::DeleteEmptyFolder {
            id: f.id.clone(),
            expected_generation: f.generation,
        })
        .unwrap();
    }
    let before = r.checkpoint_version();
    let pending = r.pending_ids();
    let remote_applies = remote.applies;
    let found = r.deleted_folder_publications(&config.id).unwrap();
    assert_eq!(found.len(), 1);
    assert_eq!(found[0].local_folder, deleted.id);
    assert_eq!(
        found[0].remote_folder.name.as_str(),
        "Explicitly deleted locally"
    );
    assert!(r.deleted_folder_publications(&other.id).unwrap().is_empty());
    assert!(r.deleted_folder_publications("unknown").is_err());
    assert_eq!(r.checkpoint_version(), before);
    assert_eq!(r.pending_ids(), pending);
    assert_eq!(remote.applies, remote_applies);
    let review = r
        .prepare_metadata_publication(
            &Target::DeletedFolder(found[0].local_folder.clone()),
            &config.id,
        )
        .unwrap();
    let job = r.queue_metadata_publication(&review).unwrap();
    assert!(matches!(
        r.resume_push(&job, &mut remote).unwrap(),
        SyncResult::Applied(_)
    ));
    assert!(
        r.deleted_folder_publications(&config.id)
            .unwrap()
            .is_empty()
    );
    assert!(r.folder(&present.id).is_ok());
}

#[test]
fn reviewed_lab_root_mount_preserves_root_folder_and_direct_document_under_unbound_and_bound_destinations()
 {
    use repo::folders::{TreeLimits, TreePhase, TreeRoot};
    let disk = Disk::new();
    let mut r = Repository::open(&disk.0).unwrap();
    let mut remote = Remote::new();
    let config = remote.add(&mut r);
    let source = remote
        .db
        .create_folder(&remote.access, None, "Remote root folder")
        .unwrap();
    let bytes = protocol::save_document(
        &Molecule::dna("Root plasmid", "ACGTACGT", Topology::Linear).unwrap(),
    )
    .unwrap();
    let payload = store::Payload {
        kind: store::PayloadKind::NativeDocument,
        archive_hash: store::Digest::parse(wire::Digest::of(&bytes).as_str()).unwrap(),
        bytes,
        selection: None,
        sharing: None,
    };
    let doc = remote
        .db
        .create_document(&remote.access, None, "Root plasmid", &payload)
        .unwrap();
    let mount = folder(&mut r, None, "Private local mount");
    let job = r
        .begin_tree_clone(
            &config.id,
            TreeRoot::LabRoot,
            Some(mount.id.clone()),
            TreeLimits::default(),
        )
        .unwrap();
    assert_eq!(
        r.resume_tree_clone(&job, 1000, &mut remote).unwrap().phase,
        TreePhase::Complete
    );
    let tree = r.tree_clone_job(&job).unwrap();
    let cloned_folder = tree.folders[source.id.as_str()].clone();
    let cloned_doc = tree.documents[doc.id.as_str()].clone();
    for bound in [false, true] {
        if bound {
            publish_folder(&mut r, &mut remote, &config, &mount);
        }
        assert!(matches!(
            r.prepare_metadata_publication(&Target::Folder(cloned_folder.clone()), &config.id),
            Err(Error::NoChanges)
        ));
        assert!(matches!(
            r.prepare_metadata_publication(
                &Target::DocumentMetadata(cloned_doc.clone()),
                &config.id
            ),
            Err(Error::NoChanges)
        ));
        let current = r.folder(&cloned_folder).unwrap();
        r.mutate_local(durable::MutationOperation::RenameFolder {
            id: cloned_folder.clone(),
            expected_generation: current.generation,
            name: format!("Root folder renamed {bound}"),
        })
        .unwrap();
        let current = r.document(&cloned_doc).unwrap();
        r.mutate_local(durable::MutationOperation::RenameDocument {
            id: cloned_doc.clone(),
            expected_generation: current.generation,
            name: format!("Root plasmid renamed {bound}"),
        })
        .unwrap();
        for target in [
            Target::Folder(cloned_folder.clone()),
            Target::DocumentMetadata(cloned_doc.clone()),
        ] {
            let review = r.prepare_metadata_publication(&target, &config.id).unwrap();
            assert_eq!(review.description().remote_parent, None);
            assert_eq!(review.description().parent_binding, None);
            assert!(matches!(
                review.description().action,
                MetadataAction::RenameFolder { .. } | MetadataAction::RenameDocument { .. }
            ));
            let id = r.queue_metadata_publication(&review).unwrap();
            assert!(matches!(
                r.resume_push(&id, &mut remote).unwrap(),
                SyncResult::Applied(_)
            ));
        }
        assert_eq!(
            r.folder_binding(&cloned_folder, &config.id)
                .unwrap()
                .remote
                .parent,
            None
        );
        assert_eq!(
            r.document_binding(&cloned_doc, &config.id)
                .unwrap()
                .remote
                .folder,
            None
        );
    }
    let destination = folder(&mut r, None, "Explicit remote move destination");
    publish_folder(&mut r, &mut remote, &config, &destination);
    let expected = Some(
        r.folder_binding(&destination.id, &config.id)
            .unwrap()
            .remote
            .id
            .clone(),
    );
    let f = r.folder(&cloned_folder).unwrap();
    r.mutate_local(durable::MutationOperation::MoveFolder {
        id: cloned_folder.clone(),
        expected_generation: f.generation,
        parent: Some(destination.id.clone()),
    })
    .unwrap();
    let d = r.document(&cloned_doc).unwrap();
    r.mutate_local(durable::MutationOperation::MoveDocument {
        id: cloned_doc.clone(),
        expected_generation: d.generation,
        folder: Some(destination.id.clone()),
    })
    .unwrap();
    for target in [
        Target::Folder(cloned_folder.clone()),
        Target::DocumentMetadata(cloned_doc.clone()),
    ] {
        let review = r.prepare_metadata_publication(&target, &config.id).unwrap();
        assert_eq!(review.description().remote_parent, expected);
        assert!(matches!(
            review.description().action,
            MetadataAction::MoveFolder { .. } | MetadataAction::MoveDocument { .. }
        ));
        let id = r.queue_metadata_publication(&review).unwrap();
        assert!(matches!(
            r.resume_push(&id, &mut remote).unwrap(),
            SyncResult::Applied(_)
        ));
    }
    assert_eq!(
        r.folder_binding(&cloned_folder, &config.id)
            .unwrap()
            .remote
            .parent,
        expected
    );
    assert_eq!(
        r.document_binding(&cloned_doc, &config.id)
            .unwrap()
            .remote
            .folder,
        expected
    );
}
#[test]
fn reviewed_selected_folder_mount_retains_external_remote_parent_until_explicit_relocation() {
    use repo::folders::{TreeLimits, TreePhase, TreeRoot};
    let disk = Disk::new();
    let mut r = Repository::open(&disk.0).unwrap();
    let mut remote = Remote::new();
    let config = remote.add(&mut r);
    let outside = remote
        .db
        .create_folder(&remote.access, None, "Not cloned remote parent")
        .unwrap();
    let source = remote
        .db
        .create_folder(&remote.access, Some(&outside.id), "Selected subtree")
        .unwrap();
    let mount = folder(&mut r, None, "Local subtree mount");
    let expected = Some(wire::FolderId::new(outside.id.as_str()).unwrap());
    let job = r
        .begin_tree_clone(
            &config.id,
            TreeRoot::Folder {
                id: wire::FolderId::new(source.id.as_str()).unwrap(),
                parent: expected.clone(),
            },
            Some(mount.id.clone()),
            TreeLimits::default(),
        )
        .unwrap();
    assert_eq!(
        r.resume_tree_clone(&job, 1000, &mut remote).unwrap().phase,
        TreePhase::Complete
    );
    let local = r.tree_clone_job(&job).unwrap().folders[source.id.as_str()].clone();
    for bound in [false, true] {
        if bound {
            publish_folder(&mut r, &mut remote, &config, &mount);
        }
        assert!(matches!(
            r.prepare_metadata_publication(&Target::Folder(local.clone()), &config.id),
            Err(Error::NoChanges)
        ));
        let f = r.folder(&local).unwrap();
        r.mutate_local(durable::MutationOperation::RenameFolder {
            id: local.clone(),
            expected_generation: f.generation,
            name: format!("Subtree rename {bound}"),
        })
        .unwrap();
        let review = r
            .prepare_metadata_publication(&Target::Folder(local.clone()), &config.id)
            .unwrap();
        assert_eq!(review.description().remote_parent, expected);
        assert!(matches!(
            review.description().action,
            MetadataAction::RenameFolder { .. }
        ));
        let id = r.queue_metadata_publication(&review).unwrap();
        assert!(matches!(
            r.resume_push(&id, &mut remote).unwrap(),
            SyncResult::Applied(_)
        ));
        assert_eq!(
            r.folder_binding(&local, &config.id).unwrap().remote.parent,
            expected
        );
    }
    let destination = folder(&mut r, None, "Actual new parent");
    publish_folder(&mut r, &mut remote, &config, &destination);
    let next = Some(
        r.folder_binding(&destination.id, &config.id)
            .unwrap()
            .remote
            .id
            .clone(),
    );
    let f = r.folder(&local).unwrap();
    r.mutate_local(durable::MutationOperation::MoveFolder {
        id: local.clone(),
        expected_generation: f.generation,
        parent: Some(destination.id),
    })
    .unwrap();
    let review = r
        .prepare_metadata_publication(&Target::Folder(local.clone()), &config.id)
        .unwrap();
    assert_eq!(review.description().remote_parent, next);
    assert!(matches!(
        review.description().action,
        MetadataAction::MoveFolder { .. }
    ));
    let id = r.queue_metadata_publication(&review).unwrap();
    assert!(matches!(
        r.resume_push(&id, &mut remote).unwrap(),
        SyncResult::Applied(_)
    ));
    assert_eq!(
        r.folder_binding(&local, &config.id).unwrap().remote.parent,
        next
    );
}
