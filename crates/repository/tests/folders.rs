use collagene_core::{
    documents::workspace::{self, Workspace},
    material::{
        lineage::{self as run_lineage, EventBinding, InputEventBinding},
        run as run_material,
    },
    planning::plan::{DesignPlan, Operation, PlanInput, Source},
    protocol,
    sequence::analysis::Primer,
};
use collagene_lab_library::{self as store, mutations as durable};
use collagene_library_protocol as wire;
use collagene_library_sharing::{self as sharing, Selection};
use collagene_library_sync as sync;
use collagene_model::{Molecule, Topology};
use collagene_repository::{
    self as repository, Channel, RemoteConfig, RemoteTransport, Repository, SyncResult,
    TransportError,
    folders::{TreeIssue, TreeLimits, TreePhase, TreeRoot},
};
use std::{
    path::PathBuf,
    sync::atomic::{AtomicU64, Ordering},
};
static SERIAL: AtomicU64 = AtomicU64::new(0);
const DNA: &str = "AGTCCGATGCTACGTAGGACTTCAGACCTGTAGATCGTCA";
struct Disk(PathBuf);
impl Disk {
    fn new() -> Self {
        Self(std::env::temp_dir().join(format!(
            "collagene-tree-{}-{}.sqlite",
            std::process::id(),
            SERIAL.fetch_add(1, Ordering::Relaxed)
        )))
    }
}
impl Drop for Disk {
    fn drop(&mut self) {
        let _ = std::fs::remove_file(&self.0);
        let _ = std::fs::remove_file(format!("{}-wal", self.0.display()));
        let _ = std::fs::remove_file(format!("{}-shm", self.0.display()));
    }
}
struct Remote {
    db: store::Library,
    access: store::LabAccess,
    browses: usize,
    opens: usize,
    lose_open: Option<usize>,
    change_after_browse: Option<usize>,
    lose_apply: bool,
    cycle: Option<(wire::FolderId, wire::Folder)>,
}
impl Remote {
    fn new() -> Self {
        let mut db = store::Library::open_in_memory().unwrap();
        let lab = db.create_lab("Real native shared lab").unwrap();
        let access = db.trusted_host_access(&lab.id).unwrap();
        Self {
            db,
            access,
            browses: 0,
            opens: 0,
            lose_open: None,
            change_after_browse: None,
            lose_apply: false,
            cycle: None,
        }
    }
    fn config(&self, repo: &mut Repository) -> RemoteConfig {
        repo.add_remote(
            "Lab",
            "https://lab.example",
            wire::LabId::new(self.access.lab_id().as_str()).unwrap(),
        )
        .unwrap()
    }
    fn folder(&mut self, parent: Option<&store::FolderId>, name: &str) -> store::Folder {
        self.db.create_folder(&self.access, parent, name).unwrap()
    }
    fn document(
        &mut self,
        parent: Option<&store::FolderId>,
        name: &str,
        payload: &store::Payload,
    ) -> store::Document {
        self.db
            .create_document(&self.access, parent, name, payload)
            .unwrap()
    }
}
impl RemoteTransport for Remote {
    fn exchange(
        &mut self,
        _: &RemoteConfig,
        channel: Channel,
        bytes: &[u8],
    ) -> Result<Vec<u8>, TransportError> {
        match channel {
            Channel::Discussion => panic!("This fixture does not exercise discussion transport"),
            Channel::Sync => {
                let request =
                    sync::decode_request(bytes).map_err(|_| TransportError::Unavailable)?;
                let applying = matches!(request.request, sync::Request::Apply { .. });
                let response = sync::native::dispatch(&mut self.db, &self.access, &request)
                    .map_err(|_| TransportError::Unavailable)?;
                if applying && self.lose_apply {
                    self.lose_apply = false;
                    return Err(TransportError::Unavailable);
                }
                sync::encode_response(&response).map_err(|_| TransportError::Unavailable)
            }
            Channel::Library => {
                let request =
                    wire::decode_request(bytes).map_err(|_| TransportError::Unavailable)?;
                let mut response =
                    collagene_library_host::dispatch(&mut self.db, &self.access, &request, None)
                        .map_err(|_| TransportError::Unavailable)?;
                match &request.request {
                    wire::Request::Browse { parent, .. } => {
                        self.browses += 1;
                        if let Some((at, root)) = &self.cycle {
                            if parent.as_ref() == Some(at) {
                                let mut root = root.clone();
                                root.parent = Some(at.clone());
                                if let wire::Outcome::Success(wire::Success::Browsed(page)) =
                                    &mut response.outcome
                                {
                                    page.entries =
                                        wire::Rows::new(vec![wire::Entry::Folder(root)]).unwrap();
                                    page.next_offset = None;
                                }
                            }
                        }
                        if self.change_after_browse == Some(self.browses) {
                            self.change_after_browse = None;
                            self.db
                                .create_folder(&self.access, None, "Concurrent remote folder")
                                .unwrap();
                        }
                    }
                    wire::Request::OpenDocument { revision, .. } => {
                        assert!(
                            revision.is_some(),
                            "recursive clone must pin an immutable scanned revision"
                        );
                        self.opens += 1;
                        if self.lose_open == Some(self.opens) {
                            self.lose_open = None;
                            return Err(TransportError::Unavailable);
                        }
                    }
                    _ => {}
                }
                wire::encode_response(&response).map_err(|_| TransportError::Unavailable)
            }
        }
    }
}
fn native(name: &str, dna: &str) -> store::Payload {
    let m = Molecule::dna(name, dna, Topology::Circular).unwrap();
    let bytes = protocol::save_document(&m).unwrap();
    store::Payload {
        kind: store::PayloadKind::NativeDocument,
        archive_hash: store::Digest::parse(wire::Digest::of(&bytes).as_str()).unwrap(),
        bytes,
        selection: None,
        sharing: None,
    }
}
fn pcr_archive() -> (store::Payload, run_material::Selection) {
    let stock = Molecule::dna("stock", DNA, Topology::Linear).unwrap();
    let target = Molecule::dna("Independent intended PCR", DNA, Topology::Linear).unwrap();
    let w = Workspace::from_documents(
        "Original source history",
        vec![stock.clone(), target.clone()],
        0,
    )
    .unwrap();
    let binding = |index: usize| {
        let d = &w.documents[index];
        EventBinding {
            document_id: d.id.clone(),
            event_id: d.head_event.clone().unwrap(),
            revision: d.head.clone(),
        }
    };
    let plan = DesignPlan {
        name: "Actual native PCR".into(),
        target_revision: target.revision(),
        inputs: vec![PlanInput {
            id: "stock".into(),
            revision: stock.revision(),
            molecule: stock,
        }],
        steps: vec![Operation::Pcr {
            template_strand: collagene_core::chemistry::Strand::Top,
            template: Source::Input("stock".into()),
            forward: Primer::new("F", "AGTCCGAT", "").unwrap(),
            reverse: Primer::new("R", "TGACGATC", "").unwrap(),
        }],
    };
    let saved = run_lineage::record_run(
        &w,
        &plan,
        binding(1),
        vec![InputEventBinding {
            input_id: "stock".into(),
            event: binding(0),
        }],
        false,
    )
    .unwrap();
    let product = &saved.receipt.products[0];
    let origin = run_material::Selection {
        run_id: saved.receipt.run_id.clone(),
        step_index: 0,
        event: EventBinding {
            document_id: product.document_id.clone(),
            event_id: product.event_id.clone(),
            revision: product.revision.clone(),
        },
    };
    let selection = Selection {
        document_id: product.document_id.clone(),
        revision: product.revision.clone(),
        head_event: Some(product.event_id.clone()),
    };
    let export = sharing::document_closure(&saved.workspace, &selection).unwrap();
    (
        store::Payload {
            kind: store::PayloadKind::WorkspaceArchive,
            archive_hash: store::Digest::parse(wire::Digest::of(&export.archive).as_str()).unwrap(),
            bytes: export.archive,
            selection: Some(store::NativeSelection {
                document_id: selection.document_id,
                revision: store::Digest::parse(selection.revision).unwrap(),
                head_event: selection.head_event,
            }),
            sharing: Some(store::WorkspaceSharing::DocumentClosure),
        },
        origin,
    )
}
fn local_folder(
    repo: &mut Repository,
    parent: Option<store::FolderId>,
    name: &str,
) -> store::Folder {
    let durable::Applied::CreatedFolder(f) = repo
        .mutate_local(durable::MutationOperation::CreateFolder {
            parent,
            name: name.into(),
        })
        .unwrap()
    else {
        panic!()
    };
    f
}

#[test]
fn interrupted_subtree_clone_resumes_exact_snapshot_preserves_local_edits_and_real_pcr_ancestry() {
    let disk = Disk::new();
    let mut repo = Repository::open(&disk.0).unwrap();
    let mut remote = Remote::new();
    let config = remote.config(&mut repo);
    let outside = remote.folder(None, "External remote parent");
    let root = remote.folder(Some(&outside.id), "Constructs");
    let child = remote.folder(Some(&root.id), "Vectors");
    let a = remote.document(Some(&root.id), "A plasmid", &native("A", "AGTCCGAT"));
    let (original_b, origin) = pcr_archive();
    let b = remote.document(Some(&child.id), "B saved PCR", &original_b);
    let destination = local_folder(&mut repo, None, "Private local mount");
    let job = repo
        .begin_tree_clone(
            &config.id,
            TreeRoot::Folder {
                id: wire::FolderId::new(root.id.as_str()).unwrap(),
                parent: Some(wire::FolderId::new(outside.id.as_str()).unwrap()),
            },
            Some(destination.id.clone()),
            TreeLimits::default(),
        )
        .unwrap();
    remote.lose_open = Some(2);
    let progress = repo.resume_tree_clone(&job, 1000, &mut remote).unwrap();
    assert_eq!(progress.phase, TreePhase::Importing);
    assert_eq!(progress.folders_created, 2);
    assert_eq!(progress.documents_created, 1);
    assert_eq!(progress.issue, Some(TreeIssue::RemoteUnavailable));
    let captured = repo.tree_clone_job(&job).unwrap();
    let local_a = captured.documents[a.id.as_str()].clone();
    let original_a = repo.document(&local_a).unwrap();
    repo.mutate_local(durable::MutationOperation::SaveRevision {
        id: local_a.clone(),
        base_revision: original_a.head.id,
        payload: native("Unpublished private edit", "TTTTGGGG"),
    })
    .unwrap();
    // The scan is complete. The remote may now move ahead; clone must still open
    // B's immutable scanned workspace revision, not its newly changed head.
    let changed = remote
        .db
        .save(
            &remote.access,
            &b.id,
            &b.head.id,
            &native("Remote newer B", "CCCCAAAA"),
        )
        .unwrap();
    assert_ne!(changed.head.id, b.head.id);
    drop(repo);
    let mut repo = Repository::open(&disk.0).unwrap();
    assert_eq!(repo.tree_clone_jobs().unwrap().len(), 1);
    let progress = repo.resume_tree_clone(&job, 1000, &mut remote).unwrap();
    assert_eq!(progress.phase, TreePhase::Complete);
    assert_eq!(progress.folders_created, 2);
    assert_eq!(progress.documents_created, 2);
    let captured = repo.tree_clone_job(&job).unwrap();
    let local_root = &captured.folders[root.id.as_str()];
    let local_child = &captured.folders[child.id.as_str()];
    let local_b = &captured.documents[b.id.as_str()];
    assert_eq!(
        repo.folder(local_root).unwrap().parent,
        Some(destination.id)
    );
    assert_eq!(
        repo.folder(local_child).unwrap().parent,
        Some(local_root.clone())
    );
    assert_eq!(
        repo.document(local_b).unwrap().folder,
        Some(local_child.clone())
    );
    let private = repo.open_document(&local_a, None).unwrap();
    assert_eq!(
        protocol::load_document(&private.bytes)
            .unwrap()
            .sequence
            .as_str(),
        "TTTTGGGG"
    );
    let cloned = repo.open_document(local_b, None).unwrap();
    assert_eq!(cloned.archive_hash, original_b.archive_hash);
    let reopened = workspace::load_workspace(&cloned.bytes).unwrap();
    let fresh = run_material::reconstruct(&reopened, &origin).unwrap();
    assert_eq!(fresh.reference().sequence.as_str(), DNA);
    let binding = repo.document_binding(local_b, &config.id).unwrap();
    assert_eq!(binding.remote.head.id.as_str(), b.head.id.as_str());
    assert_ne!(binding.remote.head.id.as_str(), changed.head.id.as_str());
    assert!(matches!(
        repo.queue_tree_folder_metadata(&job, local_root),
        Err(collagene_repository::folders::FolderError::Repository(
            repository::Error::NoChanges
        ))
    ));
    assert_eq!(
        repo.resume_tree_clone(&job, 1000, &mut remote)
            .unwrap()
            .documents_created,
        2
    );
}

#[test]
fn catalog_changes_invalidate_scan_before_import_and_explicit_restart_uses_one_new_snapshot() {
    let disk = Disk::new();
    let mut repo = Repository::open(&disk.0).unwrap();
    let mut remote = Remote::new();
    let config = remote.config(&mut repo);
    remote.folder(None, "A unrelated");
    let root = remote.folder(None, "Selected root");
    remote.document(
        Some(&root.id),
        "Selected plasmid",
        &native("Selected", "AGTCCGAT"),
    );
    let job = repo
        .begin_tree_clone(
            &config.id,
            TreeRoot::Folder {
                id: wire::FolderId::new(root.id.as_str()).unwrap(),
                parent: None,
            },
            None,
            TreeLimits {
                page_size: 1,
                ..Default::default()
            },
        )
        .unwrap();
    remote.change_after_browse = Some(1);
    let progress = repo.resume_tree_clone(&job, 1000, &mut remote).unwrap();
    assert_eq!(progress.phase, TreePhase::Invalidated);
    assert_eq!(progress.folders_created, 0);
    assert_eq!(progress.documents_created, 0);
    assert!(repo.browse(None, 0, 200, None).unwrap().entries.is_empty());
    repo.restart_tree_scan(&job).unwrap();
    let done = repo.resume_tree_clone(&job, 1000, &mut remote).unwrap();
    assert_eq!(done.phase, TreePhase::Complete);
    assert_eq!(done.folders_created, 1);
    assert_eq!(done.documents_created, 1);
    assert_eq!(repo.browse(None, 0, 200, None).unwrap().entries.len(), 1);
}

#[test]
fn mount_placement_is_not_remote_move_and_offline_rename_uses_snapshot_generation_conflict() {
    let disk = Disk::new();
    let mut repo = Repository::open(&disk.0).unwrap();
    let mut remote = Remote::new();
    let config = remote.config(&mut repo);
    let outside = remote.folder(None, "Outside");
    let root = remote.folder(Some(&outside.id), "Root");
    let destination = local_folder(&mut repo, None, "Unbound local destination");
    let job = repo
        .begin_tree_clone(
            &config.id,
            TreeRoot::Folder {
                id: wire::FolderId::new(root.id.as_str()).unwrap(),
                parent: Some(wire::FolderId::new(outside.id.as_str()).unwrap()),
            },
            Some(destination.id),
            TreeLimits::default(),
        )
        .unwrap();
    assert_eq!(
        repo.resume_tree_clone(&job, 1000, &mut remote)
            .unwrap()
            .phase,
        TreePhase::Complete
    );
    let local = repo.tree_clone_job(&job).unwrap().folders[root.id.as_str()].clone();
    assert!(matches!(
        repo.queue_tree_folder_metadata(&job, &local),
        Err(collagene_repository::folders::FolderError::Repository(
            repository::Error::NoChanges
        ))
    ));
    let local_folder = repo.folder(&local).unwrap();
    repo.mutate_local(durable::MutationOperation::RenameFolder {
        id: local.clone(),
        expected_generation: local_folder.generation,
        name: "Private offline root name".into(),
    })
    .unwrap();
    let queued = repo.queue_tree_folder_metadata(&job, &local).unwrap();
    remote
        .db
        .rename_folder(
            &remote.access,
            &root.id,
            root.generation,
            "Other member name",
        )
        .unwrap();
    let result = repo.resume_push(&queued, &mut remote).unwrap();
    assert!(matches!(
        result,
        SyncResult::Rejected(sync::Receipt {
            outcome: sync::MutationOutcome::Rejected(sync::Rejection {
                code: sync::RejectionCode::Conflict,
                ..
            }),
            ..
        })
    ));
    assert_eq!(
        repo.folder(&local).unwrap().name,
        "Private offline root name"
    );
    let remote_root = remote.db.folder(&remote.access, &root.id).unwrap();
    assert_eq!(remote_root.name, "Other member name");
    assert_eq!(remote_root.parent, Some(outside.id));
}

#[test]
fn mounted_root_document_rename_replays_without_moving_and_explicit_moves_use_bindings() {
    let disk = Disk::new();
    let mut repo = Repository::open(&disk.0).unwrap();
    let mut remote = Remote::new();
    let config = remote.config(&mut repo);
    let child = remote.folder(None, "Shared destination");
    let document = remote.document(None, "Root plasmid", &native("Original", DNA));
    let destination = local_folder(&mut repo, None, "Unbound local mount");
    let job = repo
        .begin_tree_clone(
            &config.id,
            TreeRoot::LabRoot,
            Some(destination.id.clone()),
            TreeLimits::default(),
        )
        .unwrap();
    assert_eq!(
        repo.resume_tree_clone(&job, 1000, &mut remote)
            .unwrap()
            .phase,
        TreePhase::Complete
    );
    let cloned = repo.tree_clone_job(&job).unwrap();
    let local = cloned.documents[document.id.as_str()].clone();
    let local_child = cloned.folders[child.id.as_str()].clone();
    assert!(repo.folder_binding(&destination.id, &config.id).is_none());
    assert!(matches!(
        repo.queue_tree_document_metadata(&job, &local),
        Err(repository::folders::FolderError::Repository(
            repository::Error::NoChanges
        ))
    ));
    let original = repo.document(&local).unwrap();
    repo.mutate_local(durable::MutationOperation::RenameDocument {
        id: local.clone(),
        expected_generation: original.generation,
        name: "Offline renamed plasmid".into(),
    })
    .unwrap();
    let queued = repo.queue_tree_document_metadata(&job, &local).unwrap();
    remote.lose_apply = true;
    assert!(repo.resume_push(&queued, &mut remote).is_err());
    let committed = remote.db.document(&remote.access, &document.id).unwrap();
    assert_eq!(committed.name, "Offline renamed plasmid");
    assert_eq!(committed.folder, None);
    assert_eq!(committed.head.id, document.head.id);
    drop(repo);
    let mut repo = Repository::open(&disk.0).unwrap();
    assert!(matches!(
        repo.resume_push(&queued, &mut remote).unwrap(),
        SyncResult::Applied(_)
    ));
    assert_eq!(
        remote.db.document(&remote.access, &document.id).unwrap(),
        committed
    );
    assert_eq!(repo.document(&local).unwrap().head.id, original.head.id);
    assert!(matches!(
        repo.queue_tree_document_metadata(&job, &local),
        Err(repository::folders::FolderError::Repository(
            repository::Error::NoChanges
        ))
    ));
    let current = repo.document(&local).unwrap();
    repo.mutate_local(durable::MutationOperation::MoveDocument {
        id: local.clone(),
        expected_generation: current.generation,
        folder: Some(local_child),
    })
    .unwrap();
    let queued = repo.queue_tree_document_metadata(&job, &local).unwrap();
    assert!(matches!(
        repo.resume_push(&queued, &mut remote).unwrap(),
        SyncResult::Applied(_)
    ));
    assert_eq!(
        remote
            .db
            .document(&remote.access, &document.id)
            .unwrap()
            .folder,
        Some(child.id)
    );

    // The same path uses the acknowledged remote generation, not a fresh head
    // that would silently overwrite another member's concurrent metadata edit.
    let current = repo.document(&local).unwrap();
    repo.mutate_local(durable::MutationOperation::RenameDocument {
        id: local.clone(),
        expected_generation: current.generation,
        name: "Private second rename".into(),
    })
    .unwrap();
    let queued = repo.queue_tree_document_metadata(&job, &local).unwrap();
    let current = remote.db.document(&remote.access, &document.id).unwrap();
    remote
        .db
        .rename_document(
            &remote.access,
            &document.id,
            current.generation,
            "Other member rename",
        )
        .unwrap();
    assert!(matches!(
        repo.resume_push(&queued, &mut remote).unwrap(),
        SyncResult::Rejected(sync::Receipt {
            outcome: sync::MutationOutcome::Rejected(sync::Rejection {
                code: sync::RejectionCode::Conflict,
                ..
            }),
            ..
        })
    ));
    assert_eq!(repo.document(&local).unwrap().name, "Private second rename");
    assert_eq!(
        remote
            .db
            .document(&remote.access, &document.id)
            .unwrap()
            .name,
        "Other member rename"
    );
}

#[test]
fn duplicate_cycle_reply_is_blocked_before_any_local_nodes_are_created() {
    let disk = Disk::new();
    let mut repo = Repository::open(&disk.0).unwrap();
    let mut remote = Remote::new();
    let config = remote.config(&mut repo);
    let root = remote.folder(None, "Root");
    let child = remote.folder(Some(&root.id), "Child");
    remote.cycle = Some((
        wire::FolderId::new(child.id.as_str()).unwrap(),
        wire::Folder {
            id: wire::FolderId::new(root.id.as_str()).unwrap(),
            parent: None,
            name: wire::Text::new("Root").unwrap(),
            generation: root.generation,
        },
    ));
    let job = repo
        .begin_tree_clone(
            &config.id,
            TreeRoot::Folder {
                id: wire::FolderId::new(root.id.as_str()).unwrap(),
                parent: None,
            },
            None,
            TreeLimits::default(),
        )
        .unwrap();
    let blocked = repo.resume_tree_clone(&job, 1000, &mut remote).unwrap();
    assert_eq!(blocked.phase, TreePhase::Blocked);
    assert_eq!(blocked.issue, Some(TreeIssue::InvalidGraph));
    assert_eq!(blocked.folders_created, 0);
    assert!(repo.browse(None, 0, 200, None).unwrap().entries.is_empty());
    remote.cycle = None;
    repo.restart_tree_scan(&job).unwrap();
    assert_eq!(
        repo.resume_tree_clone(&job, 1000, &mut remote)
            .unwrap()
            .phase,
        TreePhase::Complete
    );
}

#[test]
fn virtual_root_maps_to_explicit_local_destination_and_changed_parent_needs_review() {
    let disk = Disk::new();
    let mut repo = Repository::open(&disk.0).unwrap();
    let mut remote = Remote::new();
    let config = remote.config(&mut repo);
    let top = remote.folder(None, "Remote top folder");
    remote.document(None, "Root plasmid", &native("Root", "AGTCCGAT"));
    let destination = local_folder(&mut repo, None, "Explicit local mount");
    let job = repo
        .begin_tree_clone(
            &config.id,
            TreeRoot::LabRoot,
            Some(destination.id.clone()),
            TreeLimits::default(),
        )
        .unwrap();
    let scanned = repo.resume_tree_clone(&job, 2, &mut remote).unwrap();
    assert_eq!(scanned.phase, TreePhase::Importing);
    assert_eq!(scanned.folders_created, 0);
    let durable::Applied::RenamedFolder(renamed) = repo
        .mutate_local(durable::MutationOperation::RenameFolder {
            id: destination.id.clone(),
            expected_generation: destination.generation,
            name: "Keep my private destination name".into(),
        })
        .unwrap()
    else {
        panic!()
    };
    let blocked = repo.resume_tree_clone(&job, 1000, &mut remote).unwrap();
    assert_eq!(blocked.phase, TreePhase::Blocked);
    assert_eq!(
        blocked.issue,
        Some(TreeIssue::LocalParentChanged {
            folder: destination.id.clone()
        })
    );
    assert_eq!(blocked.folders_created, 0);
    repo.accept_tree_parent_change(&job, &destination.id, renamed.generation)
        .unwrap();
    repo.retry_tree_clone(&job).unwrap();
    let complete = repo.resume_tree_clone(&job, 1000, &mut remote).unwrap();
    assert_eq!(complete.phase, TreePhase::Complete);
    assert_eq!(complete.folders_created, 1);
    assert_eq!(complete.documents_created, 1);
    let saved = repo.tree_clone_job(&job).unwrap();
    let local_top = &saved.folders[top.id.as_str()];
    assert_eq!(
        repo.folder(local_top).unwrap().parent,
        Some(destination.id.clone())
    );
    assert_eq!(
        repo.folder(&destination.id).unwrap().name,
        "Keep my private destination name"
    );
    assert!(matches!(
        repo.queue_tree_folder_metadata(&job, local_top),
        Err(collagene_repository::folders::FolderError::Repository(
            repository::Error::NoChanges
        ))
    ));
    assert_eq!(
        repo.browse(Some(&destination.id), 0, 200, None)
            .unwrap()
            .entries
            .len(),
        2
    );
}
