use collagene_core::documents::workspace::Workspace;
use collagene_desktop_credentials as credentials;
use collagene_desktop_sync_host::{SyncAdapter, SyncHost};
use collagene_desktop_sync_protocol as wire;
use collagene_lab_library::{self as store, mutations as durable};
use collagene_model::{Molecule, Topology};
use collagene_repository::{
    self as repository, Channel, RemoteConfig, RemoteTransport, Repository, TransportError,
};
use std::{
    path::PathBuf,
    sync::atomic::{AtomicU64, Ordering},
};
use wire::library;
static NEXT: AtomicU64 = AtomicU64::new(1);
struct Disk(PathBuf);
impl Disk {
    fn new() -> Self {
        Self(std::env::temp_dir().join(format!(
            "collagene-unlink-service-{}-{}.sqlite",
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
#[derive(Clone)]
enum Mode {
    Success,
    Missing,
    Locked,
    Concurrent(PathBuf),
}
struct Adapter {
    local: library::LabId,
    mode: Mode,
    forgot: Vec<RemoteConfig>,
    calls: usize,
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
    fn delete_credential(&mut self, r: &RemoteConfig) -> credentials::Result<()> {
        self.forgot.push(r.clone());
        match self.mode.clone() {
            Mode::Success => Ok(()),
            Mode::Missing => Err(credentials::Error::Missing),
            Mode::Locked => Err(credentials::Error::Locked),
            Mode::Concurrent(path) => {
                let mut r = Repository::open(path).unwrap();
                r.mutate_local(durable::MutationOperation::CreateFolder {
                    parent: None,
                    name: "Concurrent edit after key forget".into(),
                })
                .unwrap();
                Ok(())
            }
        }
    }
    fn clear_session(&mut self, _: Option<&RemoteConfig>) -> credentials::Result<()> {
        Ok(())
    }
}
impl RemoteTransport for Adapter {
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
fn setup() -> (Disk, SyncHost<Adapter>, RemoteConfig, store::Document) {
    let disk = Disk::new();
    let mut r = Repository::open(&disk.0).unwrap();
    let remote = r
        .add_remote(
            "Chosen lab",
            "https://lab.example",
            library::LabId::new(format!("lab-{:032x}", 7)).unwrap(),
        )
        .unwrap();
    let frame = Workspace::from_documents(
        "Private native design",
        vec![Molecule::dna("Saved plasmid", "ACGTACGTACGT", Topology::Linear).unwrap()],
        0,
    )
    .unwrap();
    let d = &frame.documents[0];
    let pin = collagene_library_sharing::Selection {
        document_id: d.id.clone(),
        revision: d.head.clone(),
        head_event: d.head_event.clone(),
    };
    let doc = match r
        .mutate_local(durable::MutationOperation::ImportDocument {
            folder: None,
            name: "Saved plasmid".into(),
            payload: repository::private_workspace(&frame, &pin).unwrap(),
        })
        .unwrap()
    {
        durable::Applied::Imported(d) => d,
        other => panic!("{other:?}"),
    };
    let adapter = Adapter {
        local: r.local_context().unwrap().lab_id,
        mode: Mode::Success,
        forgot: vec![],
        calls: 0,
    };
    (disk, SyncHost::new(r, adapter).unwrap(), remote, doc)
}
fn remote(r: &RemoteConfig) -> wire::RemoteId {
    wire::RemoteId::new(&r.id).unwrap()
}
fn request(command: wire::Command) -> wire::RequestEnvelope {
    let n = NEXT.fetch_add(1, Ordering::Relaxed);
    wire::RequestEnvelope {
        request_id: n,
        operation_id: wire::OperationId::new(format!("op-{n:032x}")).unwrap(),
        command,
    }
}
fn gated(host: &mut SyncHost<Adapter>, c: wire::Command, accept: bool) -> wire::ResponseEnvelope {
    let q = request(c);
    let bytes = wire::encode_request(&q).unwrap();
    let response = host.dispatch_raw(&bytes, &mut || accept).unwrap();
    wire::decode_response_for(&response, &q).unwrap()
}
fn call(h: &mut SyncHost<Adapter>, c: wire::Command) -> wire::ResponseEnvelope {
    gated(h, c, true)
}
fn success(r: wire::ResponseEnvelope) -> wire::Success {
    match r.outcome {
        wire::Outcome::Success(s) => s,
        other => panic!("{other:?}"),
    }
}
fn error(r: wire::ResponseEnvelope, code: wire::ErrorCode) {
    match r.outcome {
        wire::Outcome::Failure(wire::Failure::Error { code: actual, .. }) => {
            assert_eq!(actual, code)
        }
        other => panic!("{other:?}"),
    }
}
fn prepare(h: &mut SyncHost<Adapter>, r: &RemoteConfig) -> wire::RemoteRemovalReview {
    match success(call(
        h,
        wire::Command::PrepareRemoteRemoval { remote: remote(r) },
    )) {
        wire::Success::RemoteRemovalPrepared { review } => review,
        other => panic!("{other:?}"),
    }
}
fn confirm(r: &RemoteConfig, id: wire::ReviewId) -> wire::Command {
    wire::Command::ConfirmRemoteRemoval {
        remote: remote(r),
        review: id,
    }
}
fn same_payload(actual: store::Payload, expected: &store::Payload) {
    assert_eq!(actual.kind, expected.kind);
    assert_eq!(actual.bytes, expected.bytes);
    assert_eq!(actual.archive_hash, expected.archive_hash);
    assert_eq!(actual.selection, expected.selection);
    assert_eq!(actual.sharing, expected.sharing);
}
#[test]
fn approved_unlink_forgets_exact_key_before_metadata_commit_retains_native_data_and_has_durable_outcome()
 {
    let (disk, mut host, r, doc) = setup();
    let original = host.repository().open_document(&doc.id, None).unwrap();
    let review = prepare(&mut host, &r);
    // Refused owner acceptance must not insert/evict review state or forget a key.
    for _ in 0..5 {
        error(
            gated(
                &mut host,
                wire::Command::PrepareRemoteRemoval { remote: remote(&r) },
                false,
            ),
            wire::ErrorCode::Cancelled,
        );
    }
    error(
        gated(&mut host, confirm(&r, review.id.clone()), false),
        wire::ErrorCode::Cancelled,
    );
    assert!(host.adapter_mut().forgot.is_empty());
    let receipt = match success(call(&mut host, confirm(&r, review.id.clone()))) {
        wire::Success::RemoteRemoved(r) => r,
        other => panic!("{other:?}"),
    };
    assert_eq!(host.adapter_mut().forgot, vec![r.clone()]);
    assert_eq!(host.adapter_mut().calls, 0);
    assert!(host.repository().remotes().next().is_none());
    same_payload(
        host.repository().open_document(&doc.id, None).unwrap(),
        &original,
    );
    error(
        call(&mut host, confirm(&r, review.id)),
        wire::ErrorCode::ReviewChanged,
    );
    let (repo, adapter) = host.into_parts();
    drop(repo);
    let mut host = SyncHost::new(Repository::open(&disk.0).unwrap(), adapter).unwrap();
    match success(call(
        &mut host,
        wire::Command::CheckRemoteRemovalOutcome { remote: remote(&r) },
    )) {
        wire::Success::RemoteRemovalOutcome {
            receipt: Some(actual),
            ..
        } => assert_eq!(actual, receipt),
        other => panic!("{other:?}"),
    };
    assert_eq!(host.adapter_mut().forgot.len(), 1);
    same_payload(
        host.repository().open_document(&doc.id, None).unwrap(),
        &original,
    );
}
#[test]
fn locked_vault_or_source_change_keeps_configuration_and_cas_after_forget_is_reported_truthfully() {
    let (disk, mut host, r, doc) = setup();
    let original = host.repository().open_document(&doc.id, None).unwrap();
    let review = prepare(&mut host, &r);
    host.adapter_mut().mode = Mode::Locked;
    error(
        call(&mut host, confirm(&r, review.id.clone())),
        wire::ErrorCode::CredentialForgetFailed,
    );
    assert!(host.repository().remotes().any(|c| c == &r));
    same_payload(
        host.repository().open_document(&doc.id, None).unwrap(),
        &original,
    );
    // Local metadata changed before acceptance: vault must not be touched at all.
    let mut second = Repository::open(&disk.0).unwrap();
    second
        .mutate_local(durable::MutationOperation::CreateFolder {
            parent: None,
            name: "Changed before vault".into(),
        })
        .unwrap();
    let before = host.adapter_mut().forgot.len();
    error(
        call(&mut host, confirm(&r, review.id)),
        wire::ErrorCode::ReviewChanged,
    );
    assert_eq!(host.adapter_mut().forgot.len(), before);
    let review = prepare(&mut host, &r);
    host.adapter_mut().mode = Mode::Concurrent(disk.0.clone());
    error(
        call(&mut host, confirm(&r, review.id)),
        wire::ErrorCode::RemovalChangedAfterForget,
    );
    assert_eq!(host.adapter_mut().forgot.len(), before + 1);
    assert!(host.repository().remotes().any(|c| c == &r));
    same_payload(
        host.repository().open_document(&doc.id, None).unwrap(),
        &original,
    );
    match success(call(
        &mut host,
        wire::Command::CheckRemoteRemovalOutcome { remote: remote(&r) },
    )) {
        wire::Success::RemoteRemovalOutcome { receipt: None, .. } => {}
        other => panic!("{other:?}"),
    };
    assert_eq!(host.adapter_mut().calls, 0);
}
#[test]
fn held_review_ids_expire_oldest_and_cannot_cross_native_instance_or_session_reset() {
    let (disk, mut host, r, _) = setup();
    let first = prepare(&mut host, &r);
    for _ in 0..4 {
        prepare(&mut host, &r);
    }
    error(
        call(&mut host, confirm(&r, first.id)),
        wire::ErrorCode::ReviewChanged,
    );
    assert!(host.adapter_mut().forgot.is_empty());
    let held = prepare(&mut host, &r);
    let repo = Repository::open(&disk.0).unwrap();
    let adapter = Adapter {
        local: repo.local_context().unwrap().lab_id,
        mode: Mode::Success,
        forgot: vec![],
        calls: 0,
    };
    let mut foreign = SyncHost::new(repo, adapter).unwrap();
    error(
        call(&mut foreign, confirm(&r, held.id.clone())),
        wire::ErrorCode::ReviewChanged,
    );
    assert!(foreign.adapter_mut().forgot.is_empty());
    success(call(
        &mut host,
        wire::Command::ClearSession {
            remote: Some(remote(&r)),
        },
    ));
    error(
        call(&mut host, confirm(&r, held.id)),
        wire::ErrorCode::ReviewChanged,
    );
    host.adapter_mut().mode = Mode::Missing;
    let review = prepare(&mut host, &r);
    assert!(matches!(
        success(call(&mut host, confirm(&r, review.id))),
        wire::Success::RemoteRemoved(_)
    ));
}

#[test]
fn failed_direct_key_forget_expires_existing_publication_review_for_that_scope() {
    let (_disk, host, r, _) = setup();
    let (mut repo, adapter) = host.into_parts();
    let raw=br#"{"@context":{"sbol":"http://sbols.org/v3#"},"@graph":[{"@id":"https://example.org/c","@type":"sbol:Component","sbol:hasNamespace":{"@id":"https://example.org/"},"sbol:type":{"@id":"https://identifiers.org/SBO:0000251"},"sbol:hasSequence":{"@id":"https://example.org/s"}},{"@id":"https://example.org/s","@type":"sbol:Sequence","sbol:hasNamespace":{"@id":"https://example.org/"},"sbol:encoding":{"@id":"https://identifiers.org/edam:format_1207"},"sbol:elements":"ACGT"},{"@id":"https://private.invalid/note","@type":"https://private.invalid/Activity"}]}"#.to_vec();
    let hash = library::Digest::of(&raw).as_str().to_owned();
    let mut molecule = Molecule::dna("Documentary review", "ACGT", Topology::Linear).unwrap();
    for (k, v) in [
        ("sbol:source_hash", hash.as_str()),
        ("sbol:component_iri", "https://example.org/c"),
        ("sbol:sequence_iri", "https://example.org/s"),
        (
            "sbol:exchange_profile",
            "collagene-sbol-3.1-structural-jsonld-1",
        ),
    ] {
        molecule.metadata.insert(k.into(), v.into());
    }
    let mut w = Workspace::from_documents("Private original", vec![molecule], 0).unwrap();
    w.attachments.insert(
        hash,
        collagene_core::documents::workspace::WorkspaceAttachment {
            kind: "sbol-design-exchange-v1".into(),
            name: "Original.jsonld".into(),
            media_type: "application/ld+json".into(),
            bytes: raw,
        },
    );
    let d = &w.documents[0];
    let pin = collagene_library_sharing::Selection {
        document_id: d.id.clone(),
        revision: d.head.clone(),
        head_event: d.head_event.clone(),
    };
    let doc = match repo
        .mutate_local(durable::MutationOperation::ImportDocument {
            folder: None,
            name: "Documentary review".into(),
            payload: repository::private_workspace(&w, &pin).unwrap(),
        })
        .unwrap()
    {
        durable::Applied::Imported(d) => d,
        other => panic!("{other:?}"),
    };
    let mut host = SyncHost::new(repo, adapter).unwrap();
    let push = match success(call(
        &mut host,
        wire::Command::PrepareReviewedPush {
            local_document: library::DocumentId::new(doc.id.as_str()).unwrap(),
            expected_local_revision: library::RevisionId::new(doc.head.id.as_str()).unwrap(),
            remote: remote(&r),
            choices: library::Rows::new(vec![]).unwrap(),
        },
    )) {
        wire::Success::PushPrepared(p) => p,
        other => panic!("{other:?}"),
    };
    let removal = prepare(&mut host, &r);
    host.adapter_mut().mode = Mode::Locked;
    error(
        call(&mut host, confirm(&r, removal.id)),
        wire::ErrorCode::CredentialForgetFailed,
    );
    error(
        call(
            &mut host,
            wire::Command::QueueReviewedPush {
                local_document: library::DocumentId::new(doc.id.as_str()).unwrap(),
                expected_local_revision: library::RevisionId::new(doc.head.id.as_str()).unwrap(),
                remote: remote(&r),
                review: push.id,
            },
        ),
        wire::ErrorCode::ReviewChanged,
    );
    assert!(host.repository().pending_ids().is_empty());
    assert!(host.repository().remotes().any(|c| c == &r));
}
