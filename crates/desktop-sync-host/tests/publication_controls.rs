//! Actual bounded DSQ/DSR binary controls over SQLite and an injected native server.
use collagene_core::{
    Molecule, Topology,
    workspace::{Workspace, WorkspaceAttachment},
};
use collagene_desktop_credentials as credentials;
use collagene_desktop_sync_host::{SyncAdapter, SyncHost};
use collagene_desktop_sync_protocol as wire;
use collagene_lab_library::{self as store, mutations as durable};
use collagene_library_protocol as library;
use collagene_library_sync as sync;
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
            "collagene-publication-controls-{}-{}.sqlite",
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
struct Adapter {
    local: library::LabId,
    db: store::Library,
    access: store::LabAccess,
    calls: usize,
    applies: usize,
}
impl Adapter {
    fn new(local: library::LabId) -> Self {
        let mut db = store::Library::open_in_memory().unwrap();
        let lab = db.create_lab("Native fixture lab").unwrap();
        let access = db.trusted_host_access(&lab.id).unwrap();
        Self {
            local,
            db,
            access,
            calls: 0,
            applies: 0,
        }
    }
}
impl SyncAdapter for Adapter {
    fn local_repository_lab(&mut self) -> library::LabId {
        self.local.clone()
    }
    fn credential_status(
        &mut self,
        _: &RemoteConfig,
    ) -> credentials::Result<credentials::CredentialStatus> {
        Ok(credentials::CredentialStatus::Session)
    }
    fn store_credential(
        &mut self,
        _: &RemoteConfig,
        _: credentials::Secret,
        _: credentials::Persistence,
    ) -> credentials::Result<()> {
        Ok(())
    }
    fn delete_credential(&mut self, _: &RemoteConfig) -> credentials::Result<()> {
        Ok(())
    }
    fn clear_session(&mut self, _: Option<&RemoteConfig>) -> credentials::Result<()> {
        Ok(())
    }
}
impl RemoteTransport for Adapter {
    fn exchange(
        &mut self,
        _: &RemoteConfig,
        channel: Channel,
        bytes: &[u8],
    ) -> Result<Vec<u8>, TransportError> {
        self.calls += 1;
        match channel {
            Channel::Library => {
                collagene_library_host::dispatch_frame(&mut self.db, &self.access, bytes, None)
                    .map_err(|_| TransportError::Unavailable)
            }
            Channel::Sync => {
                let q = sync::decode_request(bytes).map_err(|_| TransportError::Unavailable)?;
                if matches!(q.request, sync::Request::Apply { .. }) {
                    self.applies += 1;
                }
                let r = sync::native::dispatch(&mut self.db, &self.access, &q)
                    .map_err(|_| TransportError::Unavailable)?;
                sync::encode_response(&r).map_err(|_| TransportError::Unavailable)
            }
        }
    }
}
fn setup() -> (Disk, Repository, Adapter, wire::Remote) {
    let disk = Disk::new();
    let mut repo = Repository::open(&disk.0).unwrap();
    let adapter = Adapter::new(repo.local_context().unwrap().lab_id);
    let r = repo
        .add_remote(
            "Lab",
            "https://lab.example",
            library::LabId::new(adapter.access.lab_id().as_str()).unwrap(),
        )
        .unwrap();
    let remote = wire::Remote {
        id: wire::RemoteId::new(&r.id).unwrap(),
        name: library::Text::new(&r.name).unwrap(),
        endpoint: library::Text::new(&r.endpoint).unwrap(),
        expected_lab: r.lab,
    };
    (disk, repo, adapter, remote)
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
fn request(command: wire::Command) -> wire::RequestEnvelope {
    let n = NEXT.fetch_add(1, Ordering::Relaxed);
    wire::RequestEnvelope {
        request_id: n,
        operation_id: wire::OperationId::new(format!("op-{n:032x}")).unwrap(),
        command,
    }
}
fn call(host: &mut SyncHost<Adapter>, command: wire::Command) -> wire::ResponseEnvelope {
    let q = request(command);
    let bytes = wire::encode_request(&q).unwrap();
    let reply = host.dispatch_raw(&bytes, &mut || true).unwrap();
    wire::decode_response_for(&reply, &q).unwrap()
}
fn success(r: wire::ResponseEnvelope) -> wire::Success {
    match r.outcome {
        wire::Outcome::Success(s) => s,
        other => panic!("{other:?}"),
    }
}
fn changed(r: wire::ResponseEnvelope) {
    assert!(
        matches!(
            r.outcome,
            wire::Outcome::Failure(wire::Failure::Error {
                code: wire::ErrorCode::ReviewChanged,
                ..
            })
        ),
        "{:?}",
        r.outcome
    );
}
fn metadata(
    host: &mut SyncHost<Adapter>,
    remote: &wire::Remote,
    target: wire::PublicationTarget,
) -> wire::MetadataReview {
    match success(call(
        host,
        wire::Command::PrepareMetadataPublication {
            target,
            remote: remote.id.clone(),
        },
    )) {
        wire::Success::MetadataPrepared(r) => r,
        other => panic!("{other:?}"),
    }
}
fn queued(r: wire::ResponseEnvelope) -> wire::PendingPush {
    match success(r) {
        wire::Success::Queued(p) => p,
        other => panic!("{other:?}"),
    }
}
#[test]
fn metadata_new_folder_review_current_binary_version_and_lost_queue_reply_reuse_same_outbox_nonce()
{
    let (disk, mut repo, adapter, remote) = setup();
    let local = folder(&mut repo, "Offline folder");
    let mut host = SyncHost::new(repo, adapter).unwrap();
    let q = request(wire::Command::PrepareMetadataPublication {
        target: wire::PublicationTarget::Folder(library::FolderId::new(local.id.as_str()).unwrap()),
        remote: remote.id.clone(),
    });
    let frame = wire::encode_request(&q).unwrap();
    assert_eq!(
        u16::from_le_bytes([frame[4], frame[5]]),
        wire::SCHEMA_VERSION
    );
    let mut unsupported = frame.clone();
    unsupported[4..6].copy_from_slice(&1u16.to_le_bytes());
    assert!(wire::decode_request(&unsupported).is_err());
    let reply = host.dispatch_raw(&frame, &mut || true).unwrap();
    let review = match success(wire::decode_response_for(&reply, &q).unwrap()) {
        wire::Success::MetadataPrepared(r) => r,
        other => panic!("{other:?}"),
    };
    assert!(matches!(
        review.description.action,
        wire::MetadataAction::CreateFolder { parent: None, .. }
    ));
    assert_eq!(host.adapter_mut().calls, 0);
    let first = queued(call(
        &mut host,
        wire::Command::QueueMetadataPublication {
            review: review.clone(),
        },
    ));
    let (repo, adapter) = host.into_parts();
    drop(repo);
    let mut host = SyncHost::new(Repository::open(&disk.0).unwrap(), adapter).unwrap();
    let retried = queued(call(
        &mut host,
        wire::Command::QueueMetadataPublication { review },
    ));
    assert_eq!(first.job, retried.job);
    assert_eq!(
        host.repository().pending_ids(),
        vec![first.job.as_str().to_string()]
    );
    assert_eq!(host.adapter_mut().calls, 0);
    success(call(
        &mut host,
        wire::Command::ResumePush { job: first.job },
    ));
    assert_eq!(host.adapter_mut().applies, 1);
    assert!(
        host.repository()
            .folder_binding(&local.id, remote.id.as_str())
            .is_some()
    );
    let current = request(wire::Command::Remotes);
    let frame = wire::encode_request(&current).unwrap();
    assert_eq!(
        u16::from_le_bytes([frame[4], frame[5]]),
        wire::SCHEMA_VERSION
    );
    let reply = host.dispatch_raw(&frame, &mut || true).unwrap();
    assert_eq!(
        u16::from_le_bytes([reply[4], reply[5]]),
        wire::SCHEMA_VERSION
    );
    wire::decode_response_for(&reply, &current).unwrap();
}
#[test]
fn changed_local_metadata_rejects_previous_review_without_new_outbox_or_network() {
    let (disk, mut repo, adapter, remote) = setup();
    let local = folder(&mut repo, "Before");
    let mut host = SyncHost::new(repo, adapter).unwrap();
    let review = metadata(
        &mut host,
        &remote,
        wire::PublicationTarget::Folder(library::FolderId::new(local.id.as_str()).unwrap()),
    );
    let mut second = Repository::open(&disk.0).unwrap();
    second
        .mutate_local(durable::MutationOperation::RenameFolder {
            id: local.id,
            expected_generation: local.generation,
            name: "After".into(),
        })
        .unwrap();
    changed(call(
        &mut host,
        wire::Command::QueueMetadataPublication { review },
    ));
    assert!(host.repository().pending_ids().is_empty());
    assert_eq!(host.adapter_mut().calls, 0);
}
#[test]
fn deleted_folder_choices_reach_201st_and_reject_mixed_checkpoint_continuation() {
    let (disk, mut repo, mut adapter, remote) = setup();
    // Complete real bindings/deletes; a page is not populated by forged client IDs.
    for n in 0..201 {
        let f = folder(&mut repo, &format!("Deleted folder {n:03}"));
        let job = repo.queue_folder(&f.id, remote.id.as_str()).unwrap();
        assert!(matches!(
            repo.resume_push(&job, &mut adapter).unwrap(),
            SyncResult::Applied(_)
        ));
        repo.mutate_local(durable::MutationOperation::DeleteEmptyFolder {
            id: f.id,
            expected_generation: f.generation,
        })
        .unwrap();
    }
    let mut host = SyncHost::new(repo, adapter).unwrap();
    let calls = host.adapter_mut().calls;
    let page = match success(call(
        &mut host,
        wire::Command::ListDeletedFolderPublications {
            remote: remote.id.clone(),
            page: wire::LocalPageRequest {
                offset: 0,
                limit: 200,
            },
            expected_checkpoint: None,
        },
    )) {
        wire::Success::DeletedFolderPublications(p) => p,
        other => panic!("{other:?}"),
    };
    assert_eq!(page.items.len(), 200);
    assert_eq!(page.next_offset, Some(200));
    let next = wire::Command::ListDeletedFolderPublications {
        remote: remote.id.clone(),
        page: wire::LocalPageRequest {
            offset: 200,
            limit: 200,
        },
        expected_checkpoint: page.checkpoint.clone(),
    };
    let second = match success(call(&mut host, next.clone())) {
        wire::Success::DeletedFolderPublications(p) => p,
        other => panic!("{other:?}"),
    };
    assert_eq!(second.items.len(), 1);
    assert_eq!(second.next_offset, None);
    assert!(
        !page
            .items
            .iter()
            .any(|p| p.local_folder == second.items[0].local_folder)
    );
    assert_eq!(host.adapter_mut().calls, calls);
    let mut other = Repository::open(&disk.0).unwrap();
    folder(&mut other, "Concurrent local change");
    changed(call(&mut host, next));
    assert_eq!(host.adapter_mut().calls, calls);
}
fn saved_sbol(repo: &mut Repository) -> (store::Document, String, Vec<u8>) {
    let dna = "AGTCCGATGCTACGTAGGACTTCAGACCTGTAGATCGTCA";
    let raw=serde_json::to_vec(&serde_json::json!({"@context":{"sbol":"http://sbols.org/v3#"},"@graph":[
 {"@id":"https://example.org/c","@type":"sbol:Component","sbol:hasNamespace":{"@id":"https://example.org/"},"sbol:type":{"@id":"https://identifiers.org/SBO:0000251"},"sbol:hasSequence":{"@id":"https://example.org/s"}},
 {"@id":"https://example.org/s","@type":"sbol:Sequence","sbol:hasNamespace":{"@id":"https://example.org/"},"sbol:encoding":{"@id":"https://identifiers.org/edam:format_1207"},"sbol:elements":dna},
 {"@id":"https://private.invalid/note","@type":"https://private.invalid/Activity","https://private.invalid/value":"PRIVATE-ORIGINAL-RETAINED-LOCALLY"}]})).unwrap();
    let hash = library::Digest::of(&raw).as_str().to_string();
    let mut m = Molecule::dna("SBOL plasmid", dna, Topology::Linear).unwrap();
    for (k, v) in [
        ("sbol:source_hash", hash.as_str()),
        ("sbol:component_iri", "https://example.org/c"),
        ("sbol:sequence_iri", "https://example.org/s"),
        (
            "sbol:exchange_profile",
            "collagene-sbol-3.1-structural-jsonld-1",
        ),
    ] {
        m.metadata.insert(k.into(), v.into());
    }
    let mut w = Workspace::from_documents("Private full snapshot", vec![m], 0).unwrap();
    w.attachments.insert(
        hash.clone(),
        WorkspaceAttachment {
            kind: "sbol-design-exchange-v1".into(),
            name: "Full private original.jsonld".into(),
            media_type: "application/ld+json".into(),
            bytes: raw.clone(),
        },
    );
    let d = &w.documents[0];
    let pin = collagene_library_sharing::Selection {
        document_id: d.id.clone(),
        revision: d.head.clone(),
        head_event: d.head_event.clone(),
    };
    let payload = repository::private_workspace(&w, &pin).unwrap();
    let doc = match repo
        .mutate_local(durable::MutationOperation::ImportDocument {
            folder: None,
            name: "SBOL plasmid".into(),
            payload,
        })
        .unwrap()
    {
        durable::Applied::Imported(d) => d,
        other => panic!("{other:?}"),
    };
    (doc, hash, raw)
}
fn prepare(
    host: &mut SyncHost<Adapter>,
    remote: &wire::Remote,
    doc: &store::Document,
) -> wire::PushPublicationReview {
    match success(call(
        host,
        wire::Command::PrepareReviewedPush {
            local_document: library::DocumentId::new(doc.id.as_str()).unwrap(),
            expected_local_revision: library::RevisionId::new(doc.head.id.as_str()).unwrap(),
            remote: remote.id.clone(),
            choices: library::Rows::new(vec![]).unwrap(),
        },
    )) {
        wire::Success::PushPrepared(r) => r,
        other => panic!("{other:?}"),
    }
}
fn approve(remote: &wire::Remote, doc: &store::Document, id: wire::ReviewId) -> wire::Command {
    wire::Command::QueueReviewedPush {
        local_document: library::DocumentId::new(doc.id.as_str()).unwrap(),
        expected_local_revision: library::RevisionId::new(doc.head.id.as_str()).unwrap(),
        remote: remote.id.clone(),
        review: id,
    }
}
#[test]
fn native_sbol_review_lifetime_and_queue_preserve_private_original_and_publish_only_subset() {
    let (disk, mut repo, adapter, remote) = setup();
    let (doc, original, raw) = saved_sbol(&mut repo);
    let mut host = SyncHost::new(repo, adapter).unwrap();
    let review = prepare(&mut host, &remote, &doc);
    assert_eq!(review.entries.len(), 1);
    assert_eq!(review.entries[0].original_hash.as_str(), original);
    // A token from another native service instance is not an approval capability.
    let foreign_repo = Repository::open(&disk.0).unwrap();
    let foreign_adapter = Adapter::new(foreign_repo.local_context().unwrap().lab_id);
    let mut foreign = SyncHost::new(foreign_repo, foreign_adapter).unwrap();
    changed(call(
        &mut foreign,
        approve(&remote, &doc, review.id.clone()),
    ));
    success(call(
        &mut host,
        wire::Command::DiscardPushReview {
            review: review.id.clone(),
        },
    ));
    changed(call(&mut host, approve(&remote, &doc, review.id)));
    let review = prepare(&mut host, &remote, &doc);
    success(call(
        &mut host,
        wire::Command::ClearSession {
            remote: Some(remote.id.clone()),
        },
    ));
    changed(call(&mut host, approve(&remote, &doc, review.id)));
    let review = prepare(&mut host, &remote, &doc);
    let mut changed_repo = Repository::open(&disk.0).unwrap();
    changed_repo
        .mutate_local(durable::MutationOperation::RenameDocument {
            id: doc.id.clone(),
            expected_generation: doc.generation,
            name: "Renamed locally".into(),
        })
        .unwrap();
    changed(call(&mut host, approve(&remote, &doc, review.id)));
    let current = changed_repo.document(&doc.id).unwrap();
    let review = prepare(&mut host, &remote, &current);
    let token = review.id.clone();
    let job = queued(call(&mut host, approve(&remote, &current, token.clone())));
    changed(call(&mut host, approve(&remote, &current, token)));
    let private = host.repository().open_document(&doc.id, None).unwrap();
    let frame = collagene_core::workspace::load_workspace(&private.bytes).unwrap();
    assert_eq!(frame.attachments[&original].bytes, raw);
    assert_eq!(host.adapter_mut().calls, 0);
    success(call(&mut host, wire::Command::ResumePush { job: job.job }));
    assert_eq!(host.adapter_mut().applies, 1);
    let binding = host
        .repository()
        .document_binding(&doc.id, remote.id.as_str())
        .unwrap()
        .remote
        .clone();
    let adapter = host.adapter_mut();
    let rid = store::DocumentId::parse(binding.id.as_str()).unwrap();
    let rev = store::RevisionId::parse(binding.head.id.as_str()).unwrap();
    let public = adapter.db.payload(&adapter.access, &rid, &rev).unwrap();
    let frame = collagene_core::workspace::load_workspace(&public.bytes).unwrap();
    assert!(!frame.attachments.contains_key(&original));
    assert!(
        frame
            .attachments
            .values()
            .all(|a| !String::from_utf8_lossy(&a.bytes).contains("private.invalid"))
    );
    assert!(
        frame
            .revisions
            .values()
            .any(|m| m.metadata.get("sbol:source_hash") == Some(&original))
    );
}

#[test]
fn cancelled_preparation_cannot_hide_tokens_and_oldest_reviews_expire_without_exhaustion() {
    let (_disk, mut repo, adapter, remote) = setup();
    let (doc, _, _) = saved_sbol(&mut repo);
    let mut host = SyncHost::new(repo, adapter).unwrap();
    for _ in 0..6 {
        let q = request(wire::Command::PrepareReviewedPush {
            local_document: library::DocumentId::new(doc.id.as_str()).unwrap(),
            expected_local_revision: library::RevisionId::new(doc.head.id.as_str()).unwrap(),
            remote: remote.id.clone(),
            choices: library::Rows::new(vec![]).unwrap(),
        });
        let reply = host
            .dispatch_raw(&wire::encode_request(&q).unwrap(), &mut || false)
            .unwrap();
        assert!(matches!(
            wire::decode_response_for(&reply, &q).unwrap().outcome,
            wire::Outcome::Failure(wire::Failure::Error {
                code: wire::ErrorCode::Cancelled,
                ..
            })
        ));
    }
    let first = prepare(&mut host, &remote, &doc);
    for _ in 0..4 {
        prepare(&mut host, &remote, &doc);
    }
    changed(call(&mut host, approve(&remote, &doc, first.id)));
    assert_eq!(host.adapter_mut().calls, 0);
    assert!(host.repository().outbox_metadata().unwrap().is_empty());
}
