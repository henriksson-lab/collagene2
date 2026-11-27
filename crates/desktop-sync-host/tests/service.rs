use collagene_core::{
    Molecule, Topology,
    analysis::Primer,
    run_lineage::{self, EventBinding, InputEventBinding},
    run_material,
    verification::{DesignPlan, Operation, PlanInput, Source, Status},
    workspace::{self, RevisionEvent, Workspace},
};
use collagene_desktop_credentials as creds;
use collagene_desktop_sync_host::{SyncAdapter, SyncHost};
use collagene_desktop_sync_protocol as wire;
use collagene_lab_library::{self as store, mutations as durable};
use collagene_library_protocol as library;
use collagene_library_sharing::Selection;
use collagene_library_sync as sync;
use collagene_repository::{
    self as repository, Channel, RemoteConfig, RemoteTransport, Repository, TransportError,
};
use std::{
    collections::BTreeMap,
    path::PathBuf,
    sync::{
        Arc,
        atomic::{AtomicBool, AtomicU64, AtomicUsize, Ordering},
    },
};
static NEXT: AtomicU64 = AtomicU64::new(1);
struct Disk(PathBuf);
impl Disk {
    fn new() -> Self {
        Self(std::env::temp_dir().join(format!(
            "collagene-sync-service-{}-{}.sqlite",
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
#[derive(Default, Clone)]
struct Vault {
    values: BTreeMap<creds::Scope, Vec<u8>>,
    locked: Arc<AtomicBool>,
    stores: Arc<AtomicUsize>,
}
impl creds::Backend for Vault {
    fn read(&mut self, scope: &creds::Scope) -> creds::Result<creds::Secret> {
        if self.locked.load(Ordering::Relaxed) {
            return Err(creds::Error::Locked);
        }
        creds::Secret::new(
            self.values
                .get(scope)
                .cloned()
                .ok_or(creds::Error::Missing)?,
        )
    }
    fn store(&mut self, scope: &creds::Scope, secret: &creds::Secret) -> creds::Result<()> {
        self.stores.fetch_add(1, Ordering::Relaxed);
        if self.locked.load(Ordering::Relaxed) {
            return Err(creds::Error::Locked);
        }
        self.values
            .insert(scope.clone(), secret.as_bytes().to_vec());
        Ok(())
    }
    fn delete(&mut self, scope: &creds::Scope) -> creds::Result<()> {
        if self.locked.load(Ordering::Relaxed) {
            return Err(creds::Error::Locked);
        }
        self.values.remove(scope);
        Ok(())
    }
}
struct Server {
    db: store::Library,
    access: store::LabAccess,
    applies: usize,
    lose_apply: bool,
    wrong_reply: bool,
    wrong_lab_reply: bool,
}
impl Server {
    fn new() -> Self {
        let mut db = store::Library::open_in_memory().unwrap();
        let lab = db.create_lab("Actual native server lab").unwrap();
        let access = db.trusted_host_access(&lab.id).unwrap();
        Self {
            db,
            access,
            applies: 0,
            lose_apply: false,
            wrong_reply: false,
            wrong_lab_reply: false,
        }
    }
}
struct Adapter {
    credentials: creds::Credentials<Vault>,
    server: Server,
    calls: usize,
    vault: Vault,
}
impl Adapter {
    fn new(local: library::LabId) -> Self {
        let vault = Vault::default();
        Self {
            credentials: creds::Credentials::with_backend(local, vault.clone()),
            server: Server::new(),
            calls: 0,
            vault,
        }
    }
}
impl SyncAdapter for Adapter {
    fn local_repository_lab(&mut self) -> library::LabId {
        self.credentials.local_repository_lab().clone()
    }
    fn credential_status(&mut self, r: &RemoteConfig) -> creds::Result<creds::CredentialStatus> {
        self.credentials.status(r)
    }
    fn store_credential(
        &mut self,
        r: &RemoteConfig,
        s: creds::Secret,
        p: creds::Persistence,
    ) -> creds::Result<()> {
        self.credentials.store(r, s, p)
    }
    fn delete_credential(&mut self, r: &RemoteConfig) -> creds::Result<()> {
        self.credentials.delete(r)
    }
    fn clear_session(&mut self, r: Option<&RemoteConfig>) -> creds::Result<()> {
        match r {
            Some(r) => self.credentials.forget_session(r),
            None => {
                self.credentials.clear_session();
                Ok(())
            }
        }
    }
}
impl RemoteTransport for Adapter {
    fn exchange(
        &mut self,
        r: &RemoteConfig,
        c: Channel,
        bytes: &[u8],
    ) -> Result<Vec<u8>, TransportError> {
        self.calls += 1;
        self.credentials
            .bearer(r)
            .map_err(|_| TransportError::Unavailable)?;
        match c {
            Channel::Library => {
                let req =
                    library::decode_request(bytes).map_err(|_| TransportError::Unavailable)?;
                let mut response = collagene_library_host::dispatch(
                    &mut self.server.db,
                    &self.server.access,
                    &req,
                    None,
                )
                .map_err(|_| TransportError::Unavailable)?;
                if self.server.wrong_reply {
                    response.request_id += 1;
                }
                if self.server.wrong_lab_reply {
                    response.context.as_mut().unwrap().lab_id =
                        library::LabId::new("lab-ffffffffffffffffffffffffffffffff").unwrap();
                }
                library::encode_response(&response).map_err(|_| TransportError::Unavailable)
            }
            Channel::Sync => {
                let req = sync::decode_request(bytes).map_err(|_| TransportError::Unavailable)?;
                let applying = matches!(req.request, sync::Request::Apply { .. });
                if applying {
                    self.server.applies += 1;
                }
                let reply = sync::native::dispatch(&mut self.server.db, &self.server.access, &req)
                    .map_err(|_| TransportError::Unavailable)?;
                if applying && self.server.lose_apply {
                    self.server.lose_apply = false;
                    return Err(TransportError::Unavailable);
                }
                sync::encode_response(&reply).map_err(|_| TransportError::Unavailable)
            }
        }
    }
}
fn base() -> Workspace {
    Workspace::from_documents(
        "Full private design",
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
    .unwrap()
}
fn selection(w: &Workspace, id: &str) -> Selection {
    let d = w.documents.iter().find(|d| d.id == id).unwrap();
    Selection {
        document_id: d.id.clone(),
        revision: d.head.clone(),
        head_event: d.head_event.clone(),
    }
}
fn local_import(
    repo: &mut Repository,
    w: &Workspace,
    pin: &Selection,
    name: &str,
) -> store::Document {
    match repo
        .mutate_local(durable::MutationOperation::ImportDocument {
            folder: None,
            name: name.into(),
            payload: repository::private_workspace(w, pin).unwrap(),
        })
        .unwrap()
    {
        durable::Applied::Imported(d) => d,
        other => panic!("{other:?}"),
    }
}
fn setup() -> (Disk, SyncHost<Adapter>, store::Document) {
    let disk = Disk::new();
    let mut repo = Repository::open(&disk.0).unwrap();
    let frame = base();
    let doc = local_import(
        &mut repo,
        &frame,
        &selection(&frame, "document-1"),
        "Saved plasmid",
    );
    let adapter = Adapter::new(repo.local_context().unwrap().lab_id);
    (disk, SyncHost::new(repo, adapter).unwrap(), doc)
}
fn command(n: u64, command: wire::Command) -> wire::RequestEnvelope {
    wire::RequestEnvelope {
        request_id: n,
        operation_id: wire::OperationId::new(format!("op-{n:032x}")).unwrap(),
        command,
    }
}
fn call(service: &mut SyncHost<Adapter>, n: u64, c: wire::Command) -> wire::ResponseEnvelope {
    let request = command(n, c);
    let bytes = wire::encode_request(&request).unwrap();
    let reply = service.dispatch_raw(&bytes, &mut || true).unwrap();
    assert!(reply.len() <= wire::MAX_RESPONSE_BYTES);
    wire::decode_response_for(&reply, &request).unwrap()
}
fn success(response: wire::ResponseEnvelope) -> wire::Success {
    match response.outcome {
        wire::Outcome::Success(s) => s,
        other => panic!("{other:?}"),
    }
}
fn add(service: &mut SyncHost<Adapter>) -> wire::Remote {
    let lab = library::LabId::new(service.adapter_mut().server.access.lab_id().as_str()).unwrap();
    match success(call(
        service,
        1,
        wire::Command::AddRemote {
            name: library::Text::new("Native lab").unwrap(),
            endpoint: library::Text::new("https://configured.example").unwrap(),
            expected_lab: lab,
        },
    )) {
        wire::Success::AddedRemote(r) => r,
        other => panic!("{other:?}"),
    }
}
fn token(service: &mut SyncHost<Adapter>, remote: &wire::Remote) {
    let req = command(
        2,
        wire::Command::SetCredential {
            remote: remote.id.clone(),
            secret: wire::Secret::new("fixture-secret-no-response-copy").unwrap(),
            retention: wire::CredentialRetention::Session,
        },
    );
    let encoded = wire::encode_request(&req).unwrap();
    let reply = service.dispatch_raw(&encoded, &mut || true).unwrap();
    assert!(
        !reply
            .windows(b"fixture-secret-no-response-copy".len())
            .any(|w| w == b"fixture-secret-no-response-copy")
    );
    wire::decode_response_for(&reply, &req).unwrap();
}
fn queue(
    service: &mut SyncHost<Adapter>,
    doc: &store::Document,
    r: &wire::Remote,
    n: u64,
) -> wire::PendingPush {
    match success(call(
        service,
        n,
        wire::Command::QueuePush {
            local_document: library::DocumentId::new(doc.id.as_str()).unwrap(),
            expected_local_revision: library::RevisionId::new(doc.head.id.as_str()).unwrap(),
            remote: r.id.clone(),
        },
    )) {
        wire::Success::Queued(p) => p,
        other => panic!("{other:?}"),
    }
}
#[test]
fn local_metadata_is_offline_and_native_credentials_are_scoped_secret_free_and_fail_closed() {
    let (_disk, mut service, _) = setup();
    let req = command(1, wire::Command::Remotes);
    service
        .dispatch_raw(&wire::encode_request(&req).unwrap(), &mut || {
            panic!("read-only metadata accepted a mutation")
        })
        .unwrap();
    assert_eq!(service.adapter_mut().calls, 0);
    let remote = add(&mut service);
    token(&mut service, &remote);
    assert_eq!(
        service.adapter_mut().vault.stores.load(Ordering::Relaxed),
        0
    );
    assert!(matches!(
        success(call(
            &mut service,
            3,
            wire::Command::CredentialStatus {
                remote: remote.id.clone()
            }
        )),
        wire::Success::CredentialStatus {
            state: wire::CredentialState::Session,
            ..
        }
    ));
    service
        .adapter_mut()
        .vault
        .locked
        .store(true, Ordering::Relaxed);
    let response = call(
        &mut service,
        4,
        wire::Command::SetCredential {
            remote: remote.id.clone(),
            secret: wire::Secret::new("never-in-metadata").unwrap(),
            retention: wire::CredentialRetention::Remember,
        },
    );
    assert!(matches!(response.outcome, wire::Outcome::Failure(_)));
    service
        .adapter_mut()
        .vault
        .locked
        .store(false, Ordering::Relaxed);
    assert!(matches!(
        success(call(
            &mut service,
            5,
            wire::Command::CredentialStatus {
                remote: remote.id.clone()
            }
        )),
        wire::Success::CredentialStatus {
            state: wire::CredentialState::Missing,
            ..
        }
    ));
    token(&mut service, &remote);
    success(call(
        &mut service,
        6,
        wire::Command::ClearSession {
            remote: Some(remote.id.clone()),
        },
    ));
    assert!(matches!(
        success(call(
            &mut service,
            7,
            wire::Command::CredentialStatus {
                remote: remote.id.clone()
            }
        )),
        wire::Success::CredentialStatus {
            state: wire::CredentialState::Missing,
            ..
        }
    ));
    success(call(
        &mut service,
        8,
        wire::Command::DeleteCredential { remote: remote.id },
    ));
    let disk = Disk::new();
    let repo = Repository::open(&disk.0).unwrap();
    let adapter =
        Adapter::new(library::LabId::new("lab-ffffffffffffffffffffffffffffffff").unwrap());
    assert!(SyncHost::new(repo, adapter).is_err());
}
#[test]
fn saved_head_guard_cancel_and_lost_upload_response_preserve_the_immutable_outbox() {
    let (disk, mut service, doc) = setup();
    let remote = add(&mut service);
    token(&mut service, &remote);
    let mut other = Repository::open(&disk.0).unwrap();
    let payload = other.open_document(&doc.id, None).unwrap();
    let current = match other
        .mutate_local(durable::MutationOperation::SaveRevision {
            id: doc.id.clone(),
            base_revision: doc.head.id.clone(),
            payload,
        })
        .unwrap()
    {
        durable::Applied::Saved { document, .. } => document,
        other => panic!("{other:?}"),
    };
    let stale = command(
        3,
        wire::Command::QueuePush {
            local_document: library::DocumentId::new(doc.id.as_str()).unwrap(),
            expected_local_revision: library::RevisionId::new(doc.head.id.as_str()).unwrap(),
            remote: remote.id.clone(),
        },
    );
    let reply = service
        .dispatch_raw(&wire::encode_request(&stale).unwrap(), &mut || {
            panic!("stale saved head reached mutation acceptance")
        })
        .unwrap();
    assert!(matches!(
        wire::decode_response_for(&reply, &stale).unwrap().outcome,
        wire::Outcome::Failure(wire::Failure::LocalConflict { .. })
    ));
    assert_eq!(service.adapter_mut().calls, 0);
    let denied = command(
        4,
        wire::Command::QueuePush {
            local_document: library::DocumentId::new(doc.id.as_str()).unwrap(),
            expected_local_revision: library::RevisionId::new(current.head.id.as_str()).unwrap(),
            remote: remote.id.clone(),
        },
    );
    let reply = service
        .dispatch_raw(&wire::encode_request(&denied).unwrap(), &mut || false)
        .unwrap();
    assert!(matches!(
        wire::decode_response_for(&reply, &denied).unwrap().outcome,
        wire::Outcome::Failure(wire::Failure::Error {
            code: wire::ErrorCode::Cancelled,
            ..
        })
    ));
    assert!(service.repository().pending_ids().is_empty());
    let queued = queue(&mut service, &current, &remote, 5);
    service.adapter_mut().server.lose_apply = true;
    assert!(matches!(
        call(
            &mut service,
            6,
            wire::Command::ResumePush {
                job: queued.job.clone()
            }
        )
        .outcome,
        wire::Outcome::Failure(_)
    ));
    assert_eq!(service.adapter_mut().server.applies, 1);
    let pending = match success(call(
        &mut service,
        7,
        wire::Command::Pending {
            page: Default::default(),
        },
    )) {
        wire::Success::Pending(p) => p,
        other => panic!("{other:?}"),
    };
    assert_eq!(pending.items.len(), 1);
    assert!(pending.items[0].attempted);
    assert!(matches!(
        call(
            &mut service,
            8,
            wire::Command::CancelUnattempted {
                job: queued.job.clone()
            }
        )
        .outcome,
        wire::Outcome::Failure(_)
    ));
    let (repo, adapter) = service.into_parts();
    drop(repo);
    let repo = Repository::open(&disk.0).unwrap();
    let mut service = SyncHost::new(repo, adapter).unwrap();
    assert!(matches!(
        success(call(
            &mut service,
            9,
            wire::Command::ResumePush { job: queued.job }
        )),
        wire::Success::Resumed(wire::PushReceipt {
            outcome: wire::PushOutcome::Applied,
            ..
        })
    ));
    assert_eq!(service.adapter_mut().server.applies, 1);
    assert!(service.repository().pending_ids().is_empty());
}
#[test]
fn remote_reads_are_fully_bound_and_pending_pagination_never_drops_a_cursor() {
    let (_disk, mut service, doc) = setup();
    let remote = add(&mut service);
    token(&mut service, &remote);
    let first = queue(&mut service, &doc, &remote, 3);
    let (repo, adapter) = service.into_parts();
    let mut repo = repo;
    let mut frame = base();
    frame.name = "Another private frame".into();
    let second = local_import(
        &mut repo,
        &frame,
        &selection(&frame, "document-1"),
        "Second local catalog item",
    );
    let mut service = SyncHost::new(repo, adapter).unwrap();
    queue(&mut service, &second, &remote, 4);
    let first_page = match success(call(
        &mut service,
        5,
        wire::Command::Pending {
            page: wire::LocalPageRequest {
                offset: 0,
                limit: 1,
            },
        },
    )) {
        wire::Success::Pending(p) => p,
        other => panic!("{other:?}"),
    };
    assert_eq!(first_page.items.len(), 1);
    let next = first_page.next_offset.unwrap();
    let second_page = match success(call(
        &mut service,
        6,
        wire::Command::Pending {
            page: wire::LocalPageRequest {
                offset: next,
                limit: 1,
            },
        },
    )) {
        wire::Success::Pending(p) => p,
        other => panic!("{other:?}"),
    };
    assert_eq!(second_page.items.len(), 1);
    assert!(second_page.next_offset.is_none());
    assert_ne!(first_page.items[0].job, second_page.items[0].job);
    success(call(
        &mut service,
        7,
        wire::Command::CancelUnattempted { job: first.job },
    ));
    let response = call(
        &mut service,
        8,
        wire::Command::Browse {
            remote: remote.id.clone(),
            parent: None,
            page: Default::default(),
        },
    );
    assert_eq!(response.context.as_ref().unwrap().remote, remote);
    assert!(matches!(
        response.outcome,
        wire::Outcome::Success(wire::Success::Browsed(_))
    ));
    let remaining = wire::PushJobId::new(service.repository().pending_ids()[0].clone()).unwrap();
    success(call(
        &mut service,
        20,
        wire::Command::ResumePush { job: remaining },
    ));
    assert!(
        matches!(success(call(&mut service,21,wire::Command::Search{remote:remote.id.clone(),query:library::Text::new("Second").unwrap(),scope:library::SearchScope::AllFolders,page:Default::default()})),wire::Success::Searched(p) if p.entries.len()==1)
    );
    let bound = service
        .repository()
        .document_binding(&second.id, remote.id.as_str())
        .unwrap()
        .remote
        .id
        .clone();
    assert!(
        matches!(success(call(&mut service,22,wire::Command::ListRevisions{remote:remote.id.clone(),document:bound,page:Default::default()})),wire::Success::Revisions(p) if p.revisions.len()==1)
    );
    service.adapter_mut().server.wrong_reply = true;
    assert!(matches!(
        call(
            &mut service,
            9,
            wire::Command::Browse {
                remote: remote.id.clone(),
                parent: None,
                page: Default::default()
            }
        )
        .outcome,
        wire::Outcome::Failure(wire::Failure::Error {
            code: wire::ErrorCode::Unavailable,
            ..
        })
    ));
    service.adapter_mut().server.wrong_reply = false;
    service.adapter_mut().server.wrong_lab_reply = true;
    assert!(matches!(
        call(
            &mut service,
            23,
            wire::Command::Browse {
                remote: remote.id,
                parent: None,
                page: Default::default()
            }
        )
        .outcome,
        wire::Outcome::Failure(wire::Failure::Error {
            code: wire::ErrorCode::Unavailable,
            ..
        })
    ));
}
fn edit(frame: &mut Workspace, id: &str, text: &str) {
    let old = selection(frame, id);
    let molecule = Molecule::dna("Edited native branch", text, Topology::Linear).unwrap();
    let revision = molecule.revision();
    let event = RevisionEvent::new(
        id,
        &revision,
        old.head_event.iter().cloned().collect(),
        "Private branch edit",
    )
    .unwrap();
    let doc = frame.documents.iter_mut().find(|d| d.id == id).unwrap();
    doc.undo.push(doc.head.clone());
    doc.undo_events.push(doc.head_event.clone().unwrap());
    doc.head = revision.clone();
    doc.head_event = Some(event.id.clone());
    frame.history.push(event);
    frame.revisions.insert(revision, molecule);
    frame.validate().unwrap();
}
#[test]
fn clone_and_explicit_pull_preserve_private_branches_and_durable_native_bookmarks() {
    let (disk, mut service, doc) = setup();
    let remote = add(&mut service);
    token(&mut service, &remote);
    let queued = queue(&mut service, &doc, &remote, 3);
    success(call(
        &mut service,
        4,
        wire::Command::ResumePush { job: queued.job },
    ));
    let remote_id = service
        .repository()
        .document_binding(&doc.id, remote.id.as_str())
        .unwrap()
        .remote
        .id
        .clone();
    let clone_disk = Disk::new();
    let clone_repo = Repository::open(&clone_disk.0).unwrap();
    let (_, mut adapter) = service.into_parts();
    adapter.credentials = creds::Credentials::with_backend(
        clone_repo.local_context().unwrap().lab_id,
        adapter.vault.clone(),
    );
    let mut clone_service = SyncHost::new(clone_repo, adapter).unwrap();
    let remote = add(&mut clone_service);
    token(&mut clone_service, &remote);
    let cloned = match success(call(
        &mut clone_service,
        5,
        wire::Command::CloneDocument {
            remote: remote.id.clone(),
            document: remote_id.clone(),
            folder: None,
        },
    )) {
        wire::Success::Cloned { local, .. } => local,
        other => panic!("{other:?}"),
    };
    let mut other = Repository::open(&clone_disk.0).unwrap();
    let id = store::DocumentId::parse(cloned.id.as_str()).unwrap();
    let body = other.open_document(&id, None).unwrap();
    let mut private = workspace::load_workspace(&body.bytes).unwrap();
    edit(
        &mut private,
        "document-1",
        "AGTCCGATGCTACGTAGGACTTCAGACCTGTAGATCGTCT",
    );
    let current = match other
        .mutate_local(durable::MutationOperation::SaveRevision {
            id: id.clone(),
            base_revision: store::RevisionId::parse(cloned.head.id.as_str()).unwrap(),
            payload: repository::private_workspace(&private, &selection(&private, "document-1"))
                .unwrap(),
        })
        .unwrap()
    {
        durable::Applied::Saved { document, .. } => document,
        other => panic!("{other:?}"),
    };
    let mut original = Repository::open(&disk.0).unwrap();
    let mut changed = base();
    edit(
        &mut changed,
        "document-1",
        "AGTCCGATGCTACGTAGGACTTCAGACCTGTAGATCGTCC",
    );
    let old = original.document(&doc.id).unwrap();
    original
        .mutate_local(durable::MutationOperation::SaveRevision {
            id: doc.id.clone(),
            base_revision: old.head.id,
            payload: repository::private_workspace(&changed, &selection(&changed, "document-1"))
                .unwrap(),
        })
        .unwrap();
    let private_queue = queue(&mut clone_service, &current, &remote, 40);
    // Advance the actual remote through another reliable native mutation.
    let exported =
        collagene_library_sharing::document_closure(&changed, &selection(&changed, "document-1"))
            .unwrap();
    let bytes = exported.archive;
    let selected = selection(&exported.workspace, "document-1");
    let payload = library::Payload {
        kind: library::PayloadKind::WorkspaceArchive,
        archive_hash: library::Digest::of(&bytes),
        bytes: library::Blob::new(bytes).unwrap(),
        selection: Some(library::NativeSelection {
            document_id: library::NativeId::new(&selected.document_id).unwrap(),
            revision: library::Digest::new(&selected.revision).unwrap(),
            head_event: selected
                .head_event
                .as_ref()
                .map(|id| library::NativeId::new(id).unwrap()),
        }),
        sharing: Some(library::WorkspaceSharing::DocumentClosure),
    };
    let adapter = clone_service.adapter_mut();
    let current_remote = adapter
        .server
        .db
        .document(
            &adapter.server.access,
            &store::DocumentId::parse(remote_id.as_str()).unwrap(),
        )
        .unwrap();
    let mutation = sync::Mutation {
        id: sync::MutationId::new("mut-ffffffffffffffffffffffffffffffff").unwrap(),
        operation: sync::MutationOperation::SaveRevision {
            id: remote_id,
            base_revision: library::RevisionId::new(current_remote.head.id.as_str()).unwrap(),
            payload,
        },
    };
    let request = sync::RequestEnvelope {
        request_id: 1,
        expected_lab: remote.expected_lab.clone(),
        request: sync::Request::Apply { mutation },
    };
    let response =
        sync::native::dispatch(&mut adapter.server.db, &adapter.server.access, &request).unwrap();
    assert!(matches!(
        response.outcome,
        sync::ResponseOutcome::Receipt(_)
    ));
    assert!(matches!(
        success(call(
            &mut clone_service,
            41,
            wire::Command::ResumePush {
                job: private_queue.job.clone()
            }
        )),
        wire::Success::Resumed(wire::PushReceipt {
            outcome: wire::PushOutcome::Rejected(_),
            ..
        })
    ));
    assert!(
        matches!(success(call(&mut clone_service,42,wire::Command::Pending{page:Default::default()})),wire::Success::Pending(p) if p.items.len()==1&&p.items[0].rejection.is_some())
    );
    success(call(
        &mut clone_service,
        43,
        wire::Command::DiscardRejected {
            job: private_queue.job,
        },
    ));
    let expected = library::RevisionId::new(current.head.id.as_str()).unwrap();
    assert!(matches!(
        call(
            &mut clone_service,
            6,
            wire::Command::Pull {
                local_document: cloned.id.clone(),
                expected_local_revision: expected.clone(),
                remote: remote.id.clone(),
                mode: wire::PullMode::CleanFastForward
            }
        )
        .outcome,
        wire::Outcome::Failure(wire::Failure::Error {
            code: wire::ErrorCode::HistoryConflict,
            ..
        })
    ));
    assert!(matches!(
        success(call(
            &mut clone_service,
            7,
            wire::Command::Pull {
                local_document: cloned.id,
                expected_local_revision: expected,
                remote: remote.id,
                mode: wire::PullMode::ExplicitCheckout
            }
        )),
        wire::Success::Pulled { .. }
    ));
    assert_eq!(clone_service.repository().bookmarks().len(), 1);
    let calls = clone_service.adapter_mut().calls;
    let checkpoint = clone_service.repository().checkpoint_version();
    let req = command(
        80,
        wire::Command::ListBookmarks {
            page: Default::default(),
            expected_checkpoint: None,
        },
    );
    let bytes = clone_service
        .dispatch_raw(&wire::encode_request(&req).unwrap(), &mut || {
            panic!("local bookmark read must not accept a mutation")
        })
        .unwrap();
    let response = wire::decode_response_for(&bytes, &req).unwrap();
    assert!(response.context.is_none());
    let wire::Success::Bookmarks(bookmarks) = success(response) else {
        panic!("expected local bookmarks")
    };
    assert_eq!(bookmarks.items.len(), 1);
    assert_eq!(
        clone_service.adapter_mut().calls,
        calls,
        "bookmark listing must be fully offline"
    );
    assert_eq!(
        clone_service.repository().checkpoint_version(),
        checkpoint,
        "bookmark read must not modify the repository"
    );
    let wire::Bookmark::NativeBranch {
        document,
        revision,
        event,
        undo_count,
        incoming_revision,
        ..
    } = &bookmarks.items[0]
    else {
        panic!("expected preserved graph branch")
    };
    let pin = selection(&private, "document-1");
    assert_eq!(document.as_str(), pin.document_id);
    assert_eq!(revision.as_str(), pin.revision);
    assert_eq!(
        event.as_ref().map(|e| e.as_str()),
        pin.head_event.as_deref()
    );
    assert_eq!(*undo_count, private.documents[0].undo.len() as u32);
    assert_eq!(
        incoming_revision.as_str(),
        selection(&changed, "document-1").revision
    );
    // Any concurrent checkpoint change invalidates the next page, even when
    // it changes unrelated repository metadata rather than a bookmark itself.
    let mut concurrent = Repository::open(&clone_disk.0).unwrap();
    concurrent
        .add_remote(
            "Other database",
            "https://other.example",
            library::LabId::new("lab-eeeeeeeeeeeeeeeeeeeeeeeeeeeeeeee").unwrap(),
        )
        .unwrap();
    assert!(matches!(
        call(
            &mut clone_service,
            81,
            wire::Command::ListBookmarks {
                page: Default::default(),
                expected_checkpoint: bookmarks.checkpoint
            }
        )
        .outcome,
        wire::Outcome::Failure(wire::Failure::Error {
            code: wire::ErrorCode::BookmarkChanged,
            ..
        })
    ));
    let (repo, adapter) = clone_service.into_parts();
    drop(repo);
    let reopened = Repository::open(&clone_disk.0).unwrap();
    assert_eq!(reopened.bookmarks().len(), 1);
    assert_eq!(
        reopened.bookmarks()[0].document.head,
        selection(&private, "document-1").revision
    );
    drop(adapter);
}
fn event(frame: &Workspace, index: usize) -> EventBinding {
    let d = &frame.documents[index];
    EventBinding {
        document_id: d.id.clone(),
        revision: d.head.clone(),
        event_id: d.head_event.clone().unwrap(),
    }
}
#[test]
fn saved_product_native_history_crosses_dsq_queue_and_clone_without_becoming_metadata_authority() {
    let disk = Disk::new();
    let mut repo = Repository::open(&disk.0).unwrap();
    let source = Molecule::dna(
        "Source",
        "AGTCCGATGCTACGTAGGACTTCAGACCTGTAGATCGTCA",
        Topology::Linear,
    )
    .unwrap();
    let target = Molecule::dna(
        "Separate target",
        source.sequence.as_str(),
        Topology::Linear,
    )
    .unwrap();
    let frame = Workspace::from_documents(
        "Actual saved stock",
        vec![source.clone(), target.clone()],
        0,
    )
    .unwrap();
    let plan = DesignPlan {
        name: "Native source PCR".into(),
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
    let product = &saved.receipt.products[0];
    let selected = Selection {
        document_id: product.document_id.clone(),
        revision: product.revision.clone(),
        head_event: Some(product.event_id.clone()),
    };
    let origin = run_material::Selection {
        run_id: saved.receipt.run_id.clone(),
        step_index: 0,
        event: EventBinding {
            document_id: product.document_id.clone(),
            event_id: product.event_id.clone(),
            revision: product.revision.clone(),
        },
    };
    let doc = local_import(&mut repo, &saved.workspace, &selected, "Reusable stock");
    let adapter = Adapter::new(repo.local_context().unwrap().lab_id);
    let mut service = SyncHost::new(repo, adapter).unwrap();
    let remote = add(&mut service);
    token(&mut service, &remote);
    let q = queue(&mut service, &doc, &remote, 3);
    success(call(
        &mut service,
        4,
        wire::Command::ResumePush { job: q.job },
    ));
    let remote_id = service
        .repository()
        .document_binding(&doc.id, remote.id.as_str())
        .unwrap()
        .remote
        .id
        .clone();
    let other_disk = Disk::new();
    let other_repo = Repository::open(&other_disk.0).unwrap();
    let (_, mut adapter) = service.into_parts();
    adapter.credentials = creds::Credentials::with_backend(
        other_repo.local_context().unwrap().lab_id,
        adapter.vault.clone(),
    );
    let mut other = SyncHost::new(other_repo, adapter).unwrap();
    let remote = add(&mut other);
    token(&mut other, &remote);
    let local = match success(call(
        &mut other,
        5,
        wire::Command::CloneDocument {
            remote: remote.id,
            document: remote_id,
            folder: None,
        },
    )) {
        wire::Success::Cloned { local, .. } => local,
        other => panic!("{other:?}"),
    };
    let payload = other
        .repository()
        .open_document(&store::DocumentId::parse(local.id.as_str()).unwrap(), None)
        .unwrap();
    let loaded = workspace::load_workspace(&payload.bytes).unwrap();
    assert_eq!(
        run_material::reconstruct(&loaded, &origin)
            .unwrap()
            .report()
            .status,
        Status::ConsistentUnderAssumptions
    );
}

#[test]
fn tree_jobs_keep_mount_placement_and_queue_only_selected_metadata_updates() {
    let (disk, mut service, _) = setup();
    let remote = add(&mut service);
    token(&mut service, &remote);
    let adapter = service.adapter_mut();
    let root = adapter
        .server
        .db
        .create_folder(&adapter.server.access, None, "Remote selected root")
        .unwrap();
    let body = collagene_core::protocol::save_document(
        &Molecule::dna("Tree DNA", "ACGTACGT", Topology::Linear).unwrap(),
    )
    .unwrap();
    let payload = library::Payload {
        kind: library::PayloadKind::NativeDocument,
        archive_hash: library::Digest::of(&body),
        bytes: library::Blob::new(body).unwrap(),
        selection: None,
        sharing: None,
    };
    let request = library::RequestEnvelope {
        request_id: 1,
        request: library::Request::ImportDocument {
            folder: Some(library::FolderId::new(root.id.as_str()).unwrap()),
            name: library::Text::new("Child").unwrap(),
            payload,
        },
    };
    let reply = collagene_library_host::dispatch(
        &mut adapter.server.db,
        &adapter.server.access,
        &request,
        None,
    )
    .unwrap();
    assert!(matches!(reply.outcome, library::Outcome::Success(_)));
    let mut other = Repository::open(&disk.0).unwrap();
    let mount = match other
        .mutate_local(durable::MutationOperation::CreateFolder {
            parent: None,
            name: "Private local mount".into(),
        })
        .unwrap()
    {
        durable::Applied::CreatedFolder(f) => f,
        other => panic!("{other:?}"),
    };
    let start = match success(call(
        &mut service,
        3,
        wire::Command::BeginTreeClone {
            remote: remote.id.clone(),
            root: wire::TreeRoot::Folder {
                id: library::FolderId::new(root.id.as_str()).unwrap(),
                parent: None,
            },
            destination: Some(library::FolderId::new(mount.id.as_str()).unwrap()),
            limits: Default::default(),
        },
    )) {
        wire::Success::TreeStarted(t) => t,
        other => panic!("{other:?}"),
    };
    let mut latest = start.clone();
    for n in 4..20 {
        latest = match success(call(
            &mut service,
            n,
            wire::Command::ResumeTreeClone {
                job: start.job.clone(),
                max_steps: 8,
            },
        )) {
            wire::Success::TreeUpdated(t) => t,
            other => panic!("{other:?}"),
        };
        if latest.phase == wire::TreePhase::Complete {
            break;
        }
    }
    assert_eq!(latest.phase, wire::TreePhase::Complete);
    assert_eq!(latest.folders_created, 1);
    assert_eq!(latest.documents_created, 1);
    let (mut repo, adapter) = service.into_parts();
    let job = repo.tree_clone_job(start.job.as_str()).unwrap();
    let local_root = job.folders[root.id.as_str()].clone();
    let current = repo.folder(&local_root).unwrap();
    assert_eq!(current.parent, Some(mount.id));
    repo.mutate_local(durable::MutationOperation::RenameFolder {
        id: local_root.clone(),
        expected_generation: current.generation,
        name: "Only requested rename".into(),
    })
    .unwrap();
    let mut service = SyncHost::new(repo, adapter).unwrap();
    let queued = match success(call(
        &mut service,
        30,
        wire::Command::QueueTreeFolderMetadata {
            job: start.job.clone(),
            local_folder: library::FolderId::new(local_root.as_str()).unwrap(),
        },
    )) {
        wire::Success::Queued(p) => p,
        other => panic!("{other:?}"),
    };
    assert!(queued.local_revision.is_none());
    success(call(
        &mut service,
        31,
        wire::Command::ResumePush { job: queued.job },
    ));
    let adapter = service.adapter_mut();
    let actual = adapter
        .server
        .db
        .folder(&adapter.server.access, &root.id)
        .unwrap();
    assert_eq!(actual.name, "Only requested rename");
    assert!(actual.parent.is_none());
    let page = match success(call(
        &mut service,
        32,
        wire::Command::TreeJobs {
            page: wire::LocalPageRequest {
                offset: 0,
                limit: 1,
            },
        },
    )) {
        wire::Success::TreeJobs(p) => p,
        other => panic!("{other:?}"),
    };
    assert_eq!(page.items.len(), 1);
    assert!(page.next_offset.is_none());
}
