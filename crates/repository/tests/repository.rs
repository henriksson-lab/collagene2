use collagene_core::documents::{
    history::transition::Transition, reference_edit::Operation as Edit,
};
use collagene_core::{
    documents::workspace::{self, RevisionEvent, Workspace},
    material::{
        lineage::{self as run_lineage, EventBinding, InputEventBinding},
        run as run_material,
    },
    planning::plan::{DesignPlan, Operation, PlanInput, Source},
    sequence::analysis::Primer,
    verification::evidence::Status,
};
use collagene_lab_library::{self as store, mutations as durable};
use collagene_library_protocol as wire;
use collagene_library_sharing::{
    Selection,
    import::{ImportOptions, SelectedPolicy},
};
use collagene_library_sync as sync;
use collagene_model::{Molecule, Topology};
use collagene_repository::{
    self as repository, Channel, RemoteConfig, RemoteTransport, Repository, SyncResult,
    TransportError,
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
            "collagene-repository-{}-{}.sqlite",
            std::process::id(),
            NEXT.fetch_add(1, Ordering::Relaxed)
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
struct NativeRemote {
    db: store::Library,
    access: store::LabAccess,
    fail_after_apply: bool,
    applies: usize,
}
impl NativeRemote {
    fn new() -> Self {
        let mut db = store::Library::open_in_memory().unwrap();
        let lab = db.create_lab("Actual shared lab").unwrap();
        let access = db.trusted_host_access(&lab.id).unwrap();
        Self {
            db,
            access,
            fail_after_apply: false,
            applies: 0,
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
}
impl RemoteTransport for NativeRemote {
    fn exchange(
        &mut self,
        _: &RemoteConfig,
        channel: Channel,
        bytes: &[u8],
    ) -> Result<Vec<u8>, TransportError> {
        match channel {
            Channel::Discussion => panic!("This fixture does not exercise discussion transport"),
            Channel::Library => {
                collagene_library_host::dispatch_frame(&mut self.db, &self.access, bytes, None)
                    .map_err(|_| TransportError::Unavailable)
            }
            Channel::Sync => {
                let request =
                    sync::decode_request(bytes).map_err(|_| TransportError::Unavailable)?;
                let applying = matches!(request.request, sync::Request::Apply { .. });
                if applying {
                    self.applies += 1;
                }
                let response = sync::native::dispatch(&mut self.db, &self.access, &request)
                    .map_err(|_| TransportError::Unavailable)?;
                if applying && self.fail_after_apply {
                    self.fail_after_apply = false;
                    return Err(TransportError::Unavailable);
                }
                sync::encode_response(&response).map_err(|_| TransportError::Unavailable)
            }
        }
    }
}
fn pin(frame: &Workspace, id: &str) -> Selection {
    let d = frame.documents.iter().find(|d| d.id == id).unwrap();
    Selection {
        document_id: d.id.clone(),
        revision: d.head.clone(),
        head_event: d.head_event.clone(),
    }
}
fn base() -> Workspace {
    Workspace::from_documents(
        "Private offline editing",
        vec![
            Molecule::dna(
                "Plasmid",
                "AGTCCGATGCTACGTAGGACTTCAGACCTGTAGATCGTCA",
                Topology::Linear,
            )
            .unwrap(),
        ],
        0,
    )
    .unwrap()
}
fn import(
    repo: &mut Repository,
    frame: &Workspace,
    name: &str,
    folder: Option<store::FolderId>,
) -> store::Document {
    let payload = repository::private_workspace(frame, &pin(frame, "document-1")).unwrap();
    match repo
        .mutate_local(durable::MutationOperation::ImportDocument {
            folder,
            name: name.into(),
            payload,
        })
        .unwrap()
    {
        durable::Applied::Imported(d) => d,
        other => panic!("{other:?}"),
    }
}
fn edit(frame: &mut Workspace, text: &str) {
    let old = pin(frame, "document-1");
    let previous = frame.head("document-1").unwrap();
    let actions = vec![
        Edit::ReplaceSequence {
            start: 0,
            end: previous.sequence.len(),
            replacement: text.into(),
        },
        Edit::Rename {
            name: "Private edited plasmid".into(),
        },
        Edit::SetTopology {
            topology: Topology::Linear,
        },
    ];
    let m = Edit::apply_all(previous, &actions).unwrap();
    let revision = m.revision();
    let e = RevisionEvent::new(
        "document-1",
        &revision,
        old.head_event.iter().cloned().collect(),
        "Offline private edit",
        Transition::ReferenceEdit {
            source: collagene_core::material::lineage::EventBinding {
                document_id: "document-1".into(),
                event_id: old.head_event.clone().unwrap(),
                revision: old.revision.clone(),
            },
            actions,
            supplemental_sources: vec![],
        },
    )
    .unwrap();
    let d = &mut frame.documents[0];
    d.undo.push(d.head.clone());
    d.undo_events.push(d.head_event.clone().unwrap());
    d.head = revision.clone();
    d.head_event = Some(e.id.clone());
    frame.revisions.insert(revision, m);
    frame.history.push(e);
    frame.validate().unwrap();
}
fn save(repo: &mut Repository, doc: &store::Document, frame: &Workspace) -> store::Document {
    match repo
        .mutate_local(durable::MutationOperation::SaveRevision {
            id: doc.id.clone(),
            base_revision: doc.head.id.clone(),
            payload: repository::private_workspace(frame, &pin(frame, "document-1")).unwrap(),
        })
        .unwrap()
    {
        durable::Applied::Saved { document, .. } => document,
        other => panic!("{other:?}"),
    }
}
#[test]
fn offline_crud_and_private_archives_survive_restart_without_transport_or_credentials() {
    let disk = Disk::new();
    let mut repo = Repository::open(&disk.0).unwrap();
    let folder = match repo
        .mutate_local(durable::MutationOperation::CreateFolder {
            parent: None,
            name: "Local plasmids".into(),
        })
        .unwrap()
    {
        durable::Applied::CreatedFolder(f) => f,
        other => panic!("{other:?}"),
    };
    let mut frame = base();
    let bytes = b"PRIVATE UNPUBLISHED UI STATE".to_vec();
    frame.attachments.insert(
        wire::Digest::of(&bytes).as_str().into(),
        workspace::WorkspaceAttachment {
            kind: "lab-library-bindings-v1".into(),
            name: "Private state".into(),
            media_type: "application/json".into(),
            bytes,
        },
    );
    let doc = import(
        &mut repo,
        &frame,
        "Offline plasmid",
        Some(folder.id.clone()),
    );
    repo.mutate_local(durable::MutationOperation::RenameFolder {
        id: folder.id.clone(),
        expected_generation: folder.generation,
        name: "Renamed offline".into(),
    })
    .unwrap();
    drop(repo);
    let mut reopened = Repository::open(&disk.0).unwrap();
    assert_eq!(reopened.folder(&folder.id).unwrap().name, "Renamed offline");
    assert_eq!(
        reopened
            .browse(Some(&folder.id), 0, 20, None)
            .unwrap()
            .entries
            .len(),
        1
    );
    let stored = reopened.open_document(&doc.id, None).unwrap();
    let loaded = workspace::load_workspace(&stored.bytes).unwrap();
    assert_eq!(
        serde_json::to_value(&loaded).unwrap(),
        serde_json::to_value(&frame).unwrap()
    );
    assert_eq!(reopened.remotes().count(), 0);
    assert!(
        reopened
            .add_remote(
                "No secrets",
                "https://user:password@lab.example",
                wire::LabId::new("lab-00000000000000000000000000000000").unwrap()
            )
            .is_err()
    );
}
#[test]
fn interrupted_creation_restarts_with_same_nonce_no_duplicate_and_keeps_later_private_edits() {
    let disk = Disk::new();
    let mut repo = Repository::open(&disk.0).unwrap();
    let mut remote = NativeRemote::new();
    let config = remote.config(&mut repo);
    let mut frame = base();
    let bytes = b"PRIVATE LOCAL UI".to_vec();
    frame.attachments.insert(
        wire::Digest::of(&bytes).as_str().into(),
        workspace::WorkspaceAttachment {
            kind: "lab-library-bindings-v1".into(),
            name: "Private UI".into(),
            media_type: "application/json".into(),
            bytes,
        },
    );
    let doc = import(&mut repo, &frame, "Shared", None);
    let queued = repo.queue_push(&doc.id, &config.id).unwrap();
    remote.fail_after_apply = true;
    assert!(repo.resume_push(&queued, &mut remote).is_err());
    assert_eq!(remote.applies, 1);
    edit(&mut frame, "AGTCCGATGCTACGTAGGACTTCAGACCTGTAGATCGTCT");
    let changed = save(&mut repo, &doc, &frame);
    drop(repo);
    let mut repo = Repository::open(&disk.0).unwrap();
    assert_eq!(repo.pending_ids(), vec![queued.clone()]);
    assert!(matches!(
        repo.resume_push(&queued, &mut remote).unwrap(),
        SyncResult::Applied(_)
    ));
    assert_eq!(remote.applies, 1);
    assert!(repo.pending_ids().is_empty());
    let binding = repo.document_binding(&doc.id, &config.id).unwrap();
    let remote_id = store::DocumentId::parse(binding.remote.id.as_str()).unwrap();
    let remote_revision = store::RevisionId::parse(binding.remote.head.id.as_str()).unwrap();
    let remote_body = remote
        .db
        .payload(&remote.access, &remote_id, &remote_revision)
        .unwrap();
    assert!(
        workspace::load_workspace(&remote_body.bytes)
            .unwrap()
            .attachments
            .is_empty()
    );
    assert!(
        !workspace::load_workspace(&repo.open_document(&doc.id, None).unwrap().bytes)
            .unwrap()
            .attachments
            .is_empty()
    );
    assert_ne!(binding.common_local_revision, changed.head.id);
    assert_eq!(repo.document(&doc.id).unwrap().head.id, changed.head.id);
    assert_eq!(
        remote
            .db
            .browse(&remote.access, None, 0, 200, None)
            .unwrap()
            .entries
            .len(),
        1
    );
    let next = repo.queue_push(&doc.id, &config.id).unwrap();
    assert_ne!(next, queued);
    assert!(matches!(
        repo.resume_push(&next, &mut remote).unwrap(),
        SyncResult::Applied(_)
    ));
}
#[test]
fn two_independent_remotes_and_folder_metadata_bindings_survive_offline_rename_move() {
    let disk = Disk::new();
    let mut repo = Repository::open(&disk.0).unwrap();
    let mut first = NativeRemote::new();
    let mut second = NativeRemote::new();
    let a = first.config(&mut repo);
    let b = second.config(&mut repo);
    let parent = match repo
        .mutate_local(durable::MutationOperation::CreateFolder {
            parent: None,
            name: "Parent".into(),
        })
        .unwrap()
    {
        durable::Applied::CreatedFolder(f) => f,
        other => panic!("{other:?}"),
    };
    let child = match repo
        .mutate_local(durable::MutationOperation::CreateFolder {
            parent: Some(parent.id.clone()),
            name: "Children".into(),
        })
        .unwrap()
    {
        durable::Applied::CreatedFolder(f) => f,
        other => panic!("{other:?}"),
    };
    let q = repo.queue_folder(&parent.id, &a.id).unwrap();
    repo.resume_push(&q, &mut first).unwrap();
    let q = repo.queue_folder(&child.id, &a.id).unwrap();
    repo.resume_push(&q, &mut first).unwrap();
    let prior_generation = repo
        .folder_binding(&child.id, &a.id)
        .unwrap()
        .remote
        .generation;
    repo.mutate_local(durable::MutationOperation::RenameFolder {
        id: child.id.clone(),
        expected_generation: child.generation,
        name: "Offline rename".into(),
    })
    .unwrap();
    let current = repo.folder(&child.id).unwrap();
    repo.mutate_local(durable::MutationOperation::MoveFolder {
        id: child.id.clone(),
        expected_generation: current.generation,
        parent: None,
    })
    .unwrap();
    let q = repo.queue_folder(&child.id, &a.id).unwrap();
    repo.resume_push(&q, &mut first).unwrap();
    let q = repo.queue_folder(&child.id, &a.id).unwrap();
    repo.resume_push(&q, &mut first).unwrap();
    let bound = repo.folder_binding(&child.id, &a.id).unwrap();
    assert!(bound.remote.generation > prior_generation);
    assert!(bound.remote.parent.is_none());
    assert_eq!(bound.remote.name.as_str(), "Offline rename");
    let doc = import(&mut repo, &base(), "Multiple remotes", None);
    let q = repo.queue_push(&doc.id, &a.id).unwrap();
    repo.resume_push(&q, &mut first).unwrap();
    let q = repo.queue_push(&doc.id, &b.id).unwrap();
    repo.resume_push(&q, &mut second).unwrap();
    assert_ne!(
        repo.document_binding(&doc.id, &a.id).unwrap().remote.id,
        repo.document_binding(&doc.id, &b.id).unwrap().remote.id
    );
    drop(repo);
    let repo = Repository::open(&disk.0).unwrap();
    assert_eq!(repo.remotes().count(), 2);
    assert!(repo.document_binding(&doc.id, &a.id).is_some());
    assert!(repo.document_binding(&doc.id, &b.id).is_some());
}
#[test]
fn competing_push_is_rejected_and_explicit_pull_durably_bookmarks_private_branch() {
    let da = Disk::new();
    let db = Disk::new();
    let mut a = Repository::open(&da.0).unwrap();
    let mut b = Repository::open(&db.0).unwrap();
    let mut remote = NativeRemote::new();
    let ca = remote.config(&mut a);
    let cb = remote.config(&mut b);
    let mut fa = base();
    let local = import(&mut a, &fa, "Shared", None);
    let q = a.queue_push(&local.id, &ca.id).unwrap();
    a.resume_push(&q, &mut remote).unwrap();
    let rid = a
        .document_binding(&local.id, &ca.id)
        .unwrap()
        .remote
        .id
        .clone();
    let other = b.clone_from(&cb.id, rid, None, &mut remote).unwrap();
    let mut fb = base();
    edit(&mut fa, "AGTCCGATGCTACGTAGGACTTCAGACCTGTAGATCGTCT");
    let local = save(&mut a, &local, &fa);
    let qa = a.queue_push(&local.id, &ca.id).unwrap();
    edit(&mut fb, "AGTCCGATGCTACGTAGGACTTCAGACCTGTAGATCGTCC");
    let other = save(&mut b, &other, &fb);
    let qb = b.queue_push(&other.id, &cb.id).unwrap();
    b.resume_push(&qb, &mut remote).unwrap();
    assert!(matches!(
        a.resume_push(&qa, &mut remote).unwrap(),
        SyncResult::Rejected(_)
    ));
    assert_eq!(a.document(&local.id).unwrap().head.id, local.head.id);
    assert_eq!(a.pending_ids(), vec![qa.clone()]);
    a.discard_rejected(&qa).unwrap();
    assert!(matches!(
        a.pull(&local.id, &ca.id, &ImportOptions::default(), &mut remote),
        Err(repository::Error::Merge(_))
    ));
    let options = ImportOptions {
        selected_policy: SelectedPolicy::Checkout {
            expected_local: pin(&fa, "document-1"),
        },
        activate_selected: true,
    };
    let updated = a.pull(&local.id, &ca.id, &options, &mut remote).unwrap();
    assert_ne!(updated.head.id, local.head.id);
    assert_eq!(a.bookmarks().len(), 1);
    drop(a);
    let a = Repository::open(&da.0).unwrap();
    assert_eq!(a.bookmarks().len(), 1);
    assert_eq!(
        a.bookmarks()[0].document.head,
        pin(&fa, "document-1").revision
    );
    let old = a.open_document(&local.id, Some(&local.head.id)).unwrap();
    assert_eq!(
        workspace::load_workspace(&old.bytes).unwrap().documents[0].head,
        pin(&fa, "document-1").revision
    );
}
fn event(w: &Workspace, index: usize) -> EventBinding {
    let d = &w.documents[index];
    EventBinding {
        document_id: d.id.clone(),
        event_id: d.head_event.clone().unwrap(),
        revision: d.head.clone(),
    }
}
#[test]
fn saved_product_evidence_survives_queue_remote_clone_restart_and_fresh_replay() {
    let source = Molecule::dna(
        "Source",
        "AGTCCGATGCTACGTAGGACTTCAGACCTGTAGATCGTCA",
        Topology::Linear,
    )
    .unwrap();
    let target = Molecule::dna(
        "Independent requirement",
        source.sequence.as_str(),
        Topology::Linear,
    )
    .unwrap();
    let frame =
        Workspace::from_documents("Scientific origin", vec![source.clone(), target.clone()], 0)
            .unwrap();
    let plan = DesignPlan {
        name: "Original native PCR".into(),
        target_revision: target.revision(),
        inputs: vec![PlanInput {
            id: "stock".into(),
            revision: source.revision(),
            molecule: source,
        }],
        steps: vec![Operation::Pcr {
            template_strand: collagene_core::chemistry::Strand::Top,
            template: Source::Input("stock".into()),
            forward: Primer::new("F", "AGTCCGAT", "").unwrap(),
            reverse: Primer::new("R", "TGACGATC", "").unwrap(),
        }],
    };
    let saved = run_lineage::record_run(
        &frame,
        &plan,
        event(&frame, 1),
        vec![InputEventBinding {
            input_id: "stock".into(),
            event: event(&frame, 0),
        }],
        false,
    )
    .unwrap();
    let p = &saved.receipt.products[0];
    let origin = run_material::Selection {
        run_id: saved.receipt.run_id.clone(),
        step_index: 0,
        event: EventBinding {
            document_id: p.document_id.clone(),
            event_id: p.event_id.clone(),
            revision: p.revision.clone(),
        },
    };
    let selection = Selection {
        document_id: p.document_id.clone(),
        revision: p.revision.clone(),
        head_event: Some(p.event_id.clone()),
    };
    let da = Disk::new();
    let db = Disk::new();
    let mut a = Repository::open(&da.0).unwrap();
    let mut b = Repository::open(&db.0).unwrap();
    let mut remote = NativeRemote::new();
    let ca = remote.config(&mut a);
    let cb = remote.config(&mut b);
    let payload = repository::private_workspace(&saved.workspace, &selection).unwrap();
    let local = match a
        .mutate_local(durable::MutationOperation::ImportDocument {
            folder: None,
            name: "Reusable PCR stock".into(),
            payload,
        })
        .unwrap()
    {
        durable::Applied::Imported(d) => d,
        other => panic!("{other:?}"),
    };
    let q = a.queue_push(&local.id, &ca.id).unwrap();
    a.resume_push(&q, &mut remote).unwrap();
    let remote_id = a
        .document_binding(&local.id, &ca.id)
        .unwrap()
        .remote
        .id
        .clone();
    let clone = b.clone_from(&cb.id, remote_id, None, &mut remote).unwrap();
    drop(b);
    let b = Repository::open(&db.0).unwrap();
    let payload = b.open_document(&clone.id, None).unwrap();
    let loaded = workspace::load_workspace(&payload.bytes).unwrap();
    let fresh = run_material::reconstruct(&loaded, &origin).unwrap();
    assert_eq!(fresh.report().status, Status::ConsistentUnderAssumptions);
}
#[test]
fn concurrent_checkpoint_change_rolls_back_local_mutation_and_never_sends_network() {
    let disk = Disk::new();
    let mut first = Repository::open(&disk.0).unwrap();
    let mut stale = Repository::open(&disk.0).unwrap();
    first
        .mutate_local(durable::MutationOperation::CreateFolder {
            parent: None,
            name: "Committed".into(),
        })
        .unwrap();
    assert!(
        stale
            .mutate_local(durable::MutationOperation::CreateFolder {
                parent: None,
                name: "Must rollback".into()
            })
            .is_err()
    );
    stale.reload().unwrap();
    let page = stale.browse(None, 0, 20, None).unwrap();
    assert_eq!(page.entries.len(), 1);
    assert!(matches!(&page.entries[0],store::Entry::Folder(f) if f.name=="Committed"));
}

#[test]
fn selective_folder_clone_is_generation_bound_and_does_not_import_unselected_children() {
    let a_disk = Disk::new();
    let b_disk = Disk::new();
    let mut a = Repository::open(&a_disk.0).unwrap();
    let mut b = Repository::open(&b_disk.0).unwrap();
    let mut remote = NativeRemote::new();
    let ca = remote.config(&mut a);
    let cb = remote.config(&mut b);
    let folder = match a
        .mutate_local(durable::MutationOperation::CreateFolder {
            parent: None,
            name: "Selected folder".into(),
        })
        .unwrap()
    {
        durable::Applied::CreatedFolder(f) => f,
        other => panic!("{other:?}"),
    };
    let q = a.queue_folder(&folder.id, &ca.id).unwrap();
    a.resume_push(&q, &mut remote).unwrap();
    let bound = a.folder_binding(&folder.id, &ca.id).unwrap().remote.clone();
    let doc = import(&mut a, &base(), "Unselected child", Some(folder.id));
    let q = a.queue_push(&doc.id, &ca.id).unwrap();
    a.resume_push(&q, &mut remote).unwrap();
    let cloned = b
        .clone_folder(&cb.id, bound.id.clone(), None, None, &mut remote)
        .unwrap();
    assert_eq!(
        b.folder_binding(&cloned.id, &cb.id)
            .unwrap()
            .remote
            .generation,
        bound.generation
    );
    assert!(
        b.browse(Some(&cloned.id), 0, 20, None)
            .unwrap()
            .entries
            .is_empty()
    );
    drop(b);
    let mut b = Repository::open(&b_disk.0).unwrap();
    assert_eq!(b.browse(None, 0, 20, None).unwrap().entries.len(), 1);
    assert_eq!(
        b.folder_binding(&cloned.id, &cb.id).unwrap().remote.id,
        bound.id
    );
}
fn bare(text: &str) -> store::Payload {
    let bytes = collagene_core::protocol::save_document(
        &Molecule::dna("Bare native plasmid", text, Topology::Linear).unwrap(),
    )
    .unwrap();
    store::Payload {
        kind: store::PayloadKind::NativeDocument,
        archive_hash: store::Digest::parse(wire::Digest::of(&bytes).as_str()).unwrap(),
        bytes,
        selection: None,
        sharing: None,
    }
}
#[test]
fn explicit_snapshot_pull_bookmarks_eventless_private_body_without_inventing_ancestry() {
    let a_disk = Disk::new();
    let b_disk = Disk::new();
    let mut a = Repository::open(&a_disk.0).unwrap();
    let mut b = Repository::open(&b_disk.0).unwrap();
    let mut remote = NativeRemote::new();
    let ca = remote.config(&mut a);
    let cb = remote.config(&mut b);
    let a_doc = match a
        .mutate_local(durable::MutationOperation::ImportDocument {
            folder: None,
            name: "Bare plasmid".into(),
            payload: bare("ACGTACGT"),
        })
        .unwrap()
    {
        durable::Applied::Imported(d) => d,
        other => panic!("{other:?}"),
    };
    let q = a.queue_push(&a_doc.id, &ca.id).unwrap();
    a.resume_push(&q, &mut remote).unwrap();
    let b_doc = b
        .clone_from(
            &cb.id,
            a.document_binding(&a_doc.id, &ca.id)
                .unwrap()
                .remote
                .id
                .clone(),
            None,
            &mut remote,
        )
        .unwrap();
    let b_doc = match b
        .mutate_local(durable::MutationOperation::SaveRevision {
            id: b_doc.id,
            base_revision: b_doc.head.id,
            payload: bare("ACGTACGA"),
        })
        .unwrap()
    {
        durable::Applied::Saved { document, .. } => document,
        other => panic!("{other:?}"),
    };
    let q = b.queue_push(&b_doc.id, &cb.id).unwrap();
    b.resume_push(&q, &mut remote).unwrap();
    let private = match a
        .mutate_local(durable::MutationOperation::SaveRevision {
            id: a_doc.id.clone(),
            base_revision: a_doc.head.id.clone(),
            payload: bare("ACGTACGC"),
        })
        .unwrap()
    {
        durable::Applied::Saved { document, .. } => document,
        other => panic!("{other:?}"),
    };
    assert!(
        a.pull_snapshot(&a_doc.id, &ca.id, &a_doc.head.id, &mut remote)
            .is_err()
    );
    a.pull_snapshot(&a_doc.id, &ca.id, &private.head.id, &mut remote)
        .unwrap();
    drop(a);
    let a = Repository::open(&a_disk.0).unwrap();
    assert_eq!(a.snapshot_bookmarks().len(), 1);
    assert_eq!(a.snapshot_bookmarks()[0].revision, private.head.id);
    let old = a.open_document(&a_doc.id, Some(&private.head.id)).unwrap();
    assert_eq!(
        collagene_core::protocol::load_document(&old.bytes)
            .unwrap()
            .sequence
            .as_str(),
        "ACGTACGC"
    );
    let current = a.open_document(&a_doc.id, None).unwrap();
    assert_eq!(
        collagene_core::protocol::load_document(&current.bytes)
            .unwrap()
            .sequence
            .as_str(),
        "ACGTACGA"
    );
    assert!(current.selection.is_none());
}

struct NoNetwork {
    calls: usize,
}
impl RemoteTransport for NoNetwork {
    fn exchange(
        &mut self,
        _: &RemoteConfig,
        _: Channel,
        _: &[u8],
    ) -> Result<Vec<u8>, TransportError> {
        self.calls += 1;
        Err(TransportError::Unavailable)
    }
}

#[test]
fn native_saved_head_guards_precede_network_and_outbox_views_survive_restart() {
    let disk = Disk::new();
    let mut repo = Repository::open(&disk.0).unwrap();
    let context = repo.local_context().unwrap();
    let mut remote = NativeRemote::new();
    let config = remote.config(&mut repo);
    let mut frame = base();
    let original = import(&mut repo, &frame, "Guarded plasmid", None);
    let checkpoint = repo.checkpoint_version();
    let first = repo
        .queue_push_pinned(&original.id, &original.head.id, &config.id, &checkpoint)
        .unwrap();
    assert!(matches!(
        repo.resume_push(&first, &mut remote).unwrap(),
        SyncResult::Applied(_)
    ));
    edit(&mut frame, "AGTCCGATGCTACGTAGGACTTCAGACCTGTAGATCGTCAGG");
    let current = save(&mut repo, &original, &frame);
    let checkpoint = repo.checkpoint_version();
    assert!(matches!(
        repo.queue_push_pinned(&current.id, &original.head.id, &config.id, &checkpoint),
        Err(repository::Error::Stale)
    ));
    let mut network = NoNetwork { calls: 0 };
    assert!(matches!(
        repo.pull_pinned(
            &current.id,
            &original.head.id,
            &config.id,
            &ImportOptions::default(),
            &mut network,
            &checkpoint,
        ),
        Err(repository::Error::Stale)
    ));
    assert_eq!(network.calls, 0);
    assert!(repo.outbox_metadata().unwrap().is_empty());
    let job = repo
        .queue_push_pinned(&current.id, &current.head.id, &config.id, &checkpoint)
        .unwrap();
    let pending = repo.outbox_metadata().unwrap();
    assert_eq!(pending.len(), 1);
    assert_eq!(pending[0].id, job);
    assert_eq!(pending[0].remote_id, config.id);
    assert_eq!(pending[0].local_revision.as_ref(), Some(&current.head.id));
    assert!(
        matches!(&pending[0].target, repository::OutboxTarget::Document(id) if id == &current.id)
    );
    assert!(!pending[0].attempted);
    assert!(pending[0].rejection.is_none());
    assert!(repo.resume_push(&job, &mut network).is_err());
    assert!(network.calls > 0);
    drop(repo);
    let reopened = Repository::open(&disk.0).unwrap();
    assert_eq!(reopened.local_context().unwrap(), context);
    assert_eq!(
        reopened.document(&current.id).unwrap().head.id,
        current.head.id
    );
    let pending = reopened.outbox_metadata().unwrap();
    assert_eq!(pending.len(), 1);
    assert_eq!(pending[0].id, job);
    assert!(pending[0].attempted);
    assert_eq!(pending[0].local_revision.as_ref(), Some(&current.head.id));
}

#[test]
fn native_internal_checkpoint_rejects_intervening_writer_without_changing_displayed_head() {
    let disk = Disk::new();
    let mut repo = Repository::open(&disk.0).unwrap();
    let remote = NativeRemote::new();
    let config = remote.config(&mut repo);
    let document = import(&mut repo, &base(), "Original label", None);
    let checkpoint = repo.checkpoint_version();
    let mut writer = Repository::open(&disk.0).unwrap();
    writer
        .mutate_local(durable::MutationOperation::RenameDocument {
            id: document.id.clone(),
            expected_generation: document.generation,
            name: "Renamed in another window".into(),
        })
        .unwrap();
    assert!(matches!(
        repo.queue_push_pinned(&document.id, &document.head.id, &config.id, &checkpoint),
        Err(repository::Error::Stale)
    ));
    assert!(repo.outbox_metadata().unwrap().is_empty());
    assert_eq!(
        repo.document(&document.id).unwrap().head.id,
        document.head.id
    );
    let refreshed_checkpoint = repo.checkpoint_version();
    assert_ne!(checkpoint, refreshed_checkpoint);
    repo.queue_push_pinned(
        &document.id,
        &document.head.id,
        &config.id,
        &refreshed_checkpoint,
    )
    .unwrap();
    assert_eq!(repo.outbox_metadata().unwrap().len(), 1);
}
