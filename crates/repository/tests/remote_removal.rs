use collagene_core::workspace::Workspace;
use collagene_lab_library::{self as store, mutations as durable};
use collagene_library_protocol as wire;
use collagene_library_sync as sync;
use collagene_model::{Molecule, Topology};
use collagene_repository::{
    self as repository, Channel, Error, RemoteConfig, RemoteTransport, Repository, SyncResult,
    TransportError,
    folders::{TreeLimits, TreePhase, TreeRoot},
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
            "collagene-remote-removal-{}-{}.sqlite",
            std::process::id(),
            NEXT.fetch_add(1, Ordering::Relaxed)
        )))
    }
}
impl Drop for Disk {
    fn drop(&mut self) {
        for p in [
            self.0.clone(),
            PathBuf::from(format!("{}-wal", self.0.display())),
            PathBuf::from(format!("{}-shm", self.0.display())),
        ] {
            let _ = std::fs::remove_file(p);
        }
    }
}
struct Remote {
    db: store::Library,
    access: store::LabAccess,
    calls: usize,
}
impl Remote {
    fn new() -> Self {
        let mut db = store::Library::open_in_memory().unwrap();
        let lab = db.create_lab("Shared lab").unwrap();
        let access = db.trusted_host_access(&lab.id).unwrap();
        Self {
            db,
            access,
            calls: 0,
        }
    }
    fn add(&self, r: &mut Repository, name: &str) -> RemoteConfig {
        r.add_remote(
            name,
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
        self.calls += 1;
        match c {
            Channel::Discussion => panic!("This fixture does not exercise discussion transport"),
            Channel::Library => {
                collagene_library_host::dispatch_frame(&mut self.db, &self.access, b, None)
                    .map_err(|_| TransportError::Unavailable)
            }
            Channel::Sync => {
                let q = sync::decode_request(b).map_err(|_| TransportError::Unavailable)?;
                let r = sync::native::dispatch(&mut self.db, &self.access, &q)
                    .map_err(|_| TransportError::Unavailable)?;
                sync::encode_response(&r).map_err(|_| TransportError::Unavailable)
            }
        }
    }
}
fn folder(r: &mut Repository, name: &str) -> store::Folder {
    match r
        .mutate_local(durable::MutationOperation::CreateFolder {
            parent: None,
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
    config: &RemoteConfig,
    f: &store::Folder,
) {
    let id = r.queue_folder(&f.id, &config.id).unwrap();
    assert!(matches!(
        r.resume_push(&id, remote).unwrap(),
        SyncResult::Applied(_)
    ));
}
fn native_doc(r: &mut Repository) -> store::Document {
    let w = Workspace::from_documents(
        "Private biological history",
        vec![
            Molecule::dna(
                "Saved plasmid",
                "AGTCCGATGCTACGTAGGACTTCAGACCTGTAGATCGTCA",
                Topology::Linear,
            )
            .unwrap(),
        ],
        0,
    )
    .unwrap();
    let d = &w.documents[0];
    let pin = collagene_library_sharing::Selection {
        document_id: d.id.clone(),
        revision: d.head.clone(),
        head_event: d.head_event.clone(),
    };
    let payload = repository::private_workspace(&w, &pin).unwrap();
    match r
        .mutate_local(durable::MutationOperation::ImportDocument {
            folder: None,
            name: "Saved plasmid".into(),
            payload,
        })
        .unwrap()
    {
        durable::Applied::Imported(d) => d,
        other => panic!("{other:?}"),
    }
}
fn clone_empty(
    r: &mut Repository,
    remote: &mut Remote,
    config: &RemoteConfig,
    name: &str,
) -> (String, store::FolderId) {
    let source = remote.db.create_folder(&remote.access, None, name).unwrap();
    let id = r
        .begin_tree_clone(
            &config.id,
            TreeRoot::Folder {
                id: wire::FolderId::new(source.id.as_str()).unwrap(),
                parent: None,
            },
            None,
            TreeLimits::default(),
        )
        .unwrap();
    assert_eq!(
        r.resume_tree_clone(&id, 1000, remote).unwrap().phase,
        TreePhase::Complete
    );
    let local = r.tree_clone_job(&id).unwrap().folders[source.id.as_str()].clone();
    (id, local)
}
fn same_payload(actual: store::Payload, expected: &store::Payload) {
    assert_eq!(actual.kind, expected.kind);
    assert_eq!(actual.bytes, expected.bytes);
    assert_eq!(actual.archive_hash, expected.archive_hash);
    assert_eq!(actual.selection, expected.selection);
    assert_eq!(actual.sharing, expected.sharing);
}
#[test]
fn reviewed_removal_retires_only_selected_metadata_preserves_native_data_and_replays_after_restart()
{
    let disk = Disk::new();
    let mut r = Repository::open(&disk.0).unwrap();
    let mut a = Remote::new();
    let mut b = Remote::new();
    let ca = a.add(&mut r, "A");
    let cb = b.add(&mut r, "B");
    let local = folder(&mut r, "Local folder");
    publish_folder(&mut r, &mut a, &ca, &local);
    publish_folder(&mut r, &mut b, &cb, &local);
    let doc = native_doc(&mut r);
    for (c, remote) in [(&ca, &mut a), (&cb, &mut b)] {
        let id = r.queue_push(&doc.id, &c.id).unwrap();
        assert!(matches!(
            r.resume_push(&id, remote).unwrap(),
            SyncResult::Applied(_)
        ));
    }
    let (a_job, a_local) = clone_empty(&mut r, &mut a, &ca, "A cloned subtree");
    let (b_job, b_local) = clone_empty(&mut r, &mut b, &cb, "B cloned subtree");
    let source = r.open_document(&doc.id, None).unwrap();
    let folder_before = r.folder(&local.id).unwrap();
    let a_local_before = r.folder(&a_local).unwrap();
    let b_local_before = r.folder(&b_local).unwrap();
    let doc_before = r.document(&doc.id).unwrap();
    let bookmarks = serde_json::to_value(r.bookmarks()).unwrap();
    let snapshots = serde_json::to_value(r.snapshot_bookmarks()).unwrap();
    let a_calls = a.calls;
    let b_calls = b.calls;
    let review = r.prepare_remote_removal(&ca.id).unwrap();
    assert_eq!(review.description().document_bindings, 1);
    assert_eq!(review.description().folder_bindings, 2);
    assert_eq!(review.description().completed_tree_jobs, 1);
    let approved = review.approve();
    let receipt = r.remove_remote_reviewed(&approved).unwrap();
    assert_eq!(receipt.remote, ca);
    assert_eq!(a.calls, a_calls);
    assert_eq!(b.calls, b_calls);
    assert_eq!(
        r.remotes().map(|c| c.id.clone()).collect::<Vec<_>>(),
        vec![cb.id.clone()]
    );
    assert!(r.folder_binding(&local.id, &ca.id).is_none());
    assert!(r.document_binding(&doc.id, &ca.id).is_none());
    assert!(r.folder_binding(&local.id, &cb.id).is_some());
    assert!(r.document_binding(&doc.id, &cb.id).is_some());
    assert!(r.tree_clone_job(&a_job).is_err());
    assert_eq!(r.tree_clone_job(&b_job).unwrap().remote, cb);
    same_payload(r.open_document(&doc.id, None).unwrap(), &source);
    assert_eq!(r.folder(&local.id).unwrap(), folder_before);
    assert_eq!(r.folder(&a_local).unwrap(), a_local_before);
    assert_eq!(r.folder(&b_local).unwrap(), b_local_before);
    assert_eq!(r.document(&doc.id).unwrap(), doc_before);
    assert_eq!(serde_json::to_value(r.bookmarks()).unwrap(), bookmarks);
    assert_eq!(
        serde_json::to_value(r.snapshot_bookmarks()).unwrap(),
        snapshots
    );
    drop(r);
    let mut r = Repository::open(&disk.0).unwrap();
    let checkpoint = r.checkpoint_version();
    assert_eq!(
        r.remote_removal_receipt(&ca.id).unwrap(),
        Some(receipt.clone())
    );
    assert_eq!(r.remove_remote_reviewed(&approved).unwrap(), receipt);
    assert_eq!(r.checkpoint_version(), checkpoint);
    same_payload(r.open_document(&doc.id, None).unwrap(), &source);
    assert_eq!(r.tree_clone_jobs().unwrap().len(), 1);
}
#[test]
fn pending_or_unfinished_jobs_prevent_removal_and_concurrent_metadata_invalidates_approval() {
    let disk = Disk::new();
    let mut r = Repository::open(&disk.0).unwrap();
    let mut remote = Remote::new();
    let config = remote.add(&mut r, "Lab");
    let local = folder(&mut r, "Pending local");
    let pending = r.queue_folder(&local.id, &config.id).unwrap();
    assert!(matches!(
        r.prepare_remote_removal(&config.id),
        Err(Error::Pending)
    ));
    r.cancel_unattempted(&pending).unwrap();
    let source = remote
        .db
        .create_folder(&remote.access, None, "Unfinished subtree")
        .unwrap();
    let job = r
        .begin_tree_clone(
            &config.id,
            TreeRoot::Folder {
                id: wire::FolderId::new(source.id.as_str()).unwrap(),
                parent: None,
            },
            None,
            TreeLimits::default(),
        )
        .unwrap();
    assert!(r.prepare_remote_removal(&config.id).is_err());
    assert!(r.remote_removal_receipt(&config.id).unwrap().is_none());
    assert_eq!(
        r.resume_tree_clone(&job, 1000, &mut remote).unwrap().phase,
        TreePhase::Complete
    );
    let review = r.prepare_remote_removal(&config.id).unwrap();
    let mut second = Repository::open(&disk.0).unwrap();
    folder(&mut second, "Concurrent local metadata");
    assert!(matches!(
        r.remove_remote_reviewed(&review.approve()),
        Err(Error::Stale)
    ));
    assert!(r.remotes().any(|c| c.id == config.id));
    assert!(r.tree_clone_job(&job).is_ok());
    assert!(r.remote_removal_receipt(&config.id).unwrap().is_none());
}
#[test]
fn newly_started_clone_after_review_and_foreign_local_context_cannot_be_removed_by_old_approval() {
    let disk = Disk::new();
    let mut r = Repository::open(&disk.0).unwrap();
    let mut remote = Remote::new();
    let config = remote.add(&mut r, "Lab");
    let review = r.prepare_remote_removal(&config.id).unwrap();
    let approved = review.approve();
    let mut second = Repository::open(&disk.0).unwrap();
    let source = remote
        .db
        .create_folder(&remote.access, None, "Late clone")
        .unwrap();
    let job = second
        .begin_tree_clone(
            &config.id,
            TreeRoot::Folder {
                id: wire::FolderId::new(source.id.as_str()).unwrap(),
                parent: None,
            },
            None,
            TreeLimits::default(),
        )
        .unwrap();
    assert!(r.remove_remote_reviewed(&approved).is_err());
    assert!(r.remotes().any(|c| c.id == config.id));
    assert!(r.tree_clone_job(&job).is_ok());
    let other_disk = Disk::new();
    let mut other = Repository::open(&other_disk.0).unwrap();
    let c = remote.add(&mut other, "Independent repository remote");
    assert!(other.remove_remote_reviewed(&approved).is_err());
    assert!(other.remotes().any(|r| r.id == c.id));
    assert!(other.remote_removal_receipt(&config.id).unwrap().is_none());
}
