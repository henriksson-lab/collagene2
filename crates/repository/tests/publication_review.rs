//! Actual saved-native PCR publication, opaque acknowledgement and durable retry.
use collagene_core::{
    Molecule, Topology,
    analysis::Primer,
    cloning_outputs,
    run_lineage::{EventBinding, InputEventBinding},
    run_material,
    verification::{DesignPlan, Operation, PlanInput, Source},
    workspace::{self, Workspace, WorkspaceAttachment},
};
use collagene_design_exchange as sbol;
use collagene_lab_library::{self as store, mutations as durable};
use collagene_library_protocol as wire;
use collagene_library_sharing::{
    self as sharing, Selection,
    sbol_subset::{PublicationChoices, RECEIPT_KIND, SUBSET_KIND},
};
use collagene_library_sync as sync;
use collagene_repository::{
    self as repo, Channel, RemoteConfig, RemoteTransport, Repository, SyncResult, TransportError,
};
use serde_json::json;
use std::{
    collections::{BTreeMap, BTreeSet},
    path::PathBuf,
    sync::atomic::{AtomicU64, Ordering},
};
const DNA: &str = "AGTCCGATGCTACGTAGGACTTCAGACCTGTAGATCGTCA";
const COMPONENT: &str = "https://example.test/c";
const SEQUENCE: &str = "https://example.test/s";
const ALTERNATIVE: &str = "https://example.test/alternative";
static NEXT: AtomicU64 = AtomicU64::new(1);
struct Disk(PathBuf);
impl Disk {
    fn new() -> Self {
        Self(std::env::temp_dir().join(format!(
            "collagene-reviewed-push-{}-{}.sqlite",
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
    lose_reply: bool,
    applies: usize,
}
impl NativeRemote {
    fn new() -> Self {
        let mut db = store::Library::open_in_memory().unwrap();
        let lab = db.create_lab("Shared actual lab").unwrap();
        let access = db.trusted_host_access(&lab.id).unwrap();
        Self {
            db,
            access,
            lose_reply: false,
            applies: 0,
        }
    }
    fn config(&self, repository: &mut Repository) -> RemoteConfig {
        repository
            .add_remote(
                "Shared lab",
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
                let apply = matches!(request.request, sync::Request::Apply { .. });
                if apply {
                    self.applies += 1;
                }
                let reply = sync::native::dispatch(&mut self.db, &self.access, &request)
                    .map_err(|_| TransportError::Unavailable)?;
                if apply && self.lose_reply {
                    self.lose_reply = false;
                    return Err(TransportError::Unavailable);
                }
                sync::encode_response(&reply).map_err(|_| TransportError::Unavailable)
            }
        }
    }
}
fn event(w: &Workspace, index: usize) -> EventBinding {
    let d = &w.documents[index];
    EventBinding {
        document_id: d.id.clone(),
        event_id: d.head_event.clone().unwrap(),
        revision: d.head.clone(),
    }
}
fn fixture(alternatives: bool) -> (Workspace, Selection, run_material::Selection, String) {
    let mut graph = json!({"@context":{"sbol":sbol::SBOL},"@graph":[
        {"@id":COMPONENT,"@type":"sbol:Component","sbol:hasNamespace":{"@id":"https://example.test"},"sbol:type":[{"@id":sbol::DNA},{"@id":sbol::LINEAR}],"sbol:hasSequence":{"@id":SEQUENCE}},
        {"@id":SEQUENCE,"@type":"sbol:Sequence","sbol:hasNamespace":{"@id":"https://example.test"},"sbol:encoding":{"@id":sbol::NUCLEIC_ENCODING},"sbol:elements":DNA},
        {"@id":"https://private.invalid/activity","@type":"https://private.invalid/Activity","https://private.invalid/note":"DO-NOT-PUBLISH-PRIVATE-SENTINEL"}
    ]});
    if alternatives {
        graph["@graph"][0]["sbol:hasSequence"] = json!([{"@id":SEQUENCE},{"@id":ALTERNATIVE}]);
        graph["@graph"].as_array_mut().unwrap().push(json!({"@id":ALTERNATIVE,"@type":"sbol:Sequence","sbol:hasNamespace":{"@id":"https://example.test"},"sbol:encoding":{"@id":sbol::NUCLEIC_ENCODING},"sbol:elements":"TTTTCCCCAAAAGGGG"}));
    }
    let raw = serde_json::to_vec(&graph).unwrap();
    let hash = wire::Digest::of(&raw).as_str().to_string();
    let stock = sbol::import(&raw)
        .unwrap()
        .project(
            &sbol::Selection {
                component: COMPONENT.into(),
                sequence: SEQUENCE.into(),
                topology: Some(Topology::Linear),
            },
            sbol::ProjectionPolicy::AllowOmissions,
        )
        .unwrap()
        .molecule;
    let target = Molecule::dna("Independent PCR target", DNA, Topology::Linear).unwrap();
    let private = Molecule::dna(
        "UNRELATED PRIVATE PLASMID",
        "CCCGGGTTTAAA",
        Topology::Circular,
    )
    .unwrap();
    let mut w = Workspace::from_documents(
        "Exact private sources",
        vec![stock.clone(), target.clone(), private],
        1,
    )
    .unwrap();
    w.attachments.insert(
        hash.clone(),
        WorkspaceAttachment {
            kind: sharing::SBOL_SOURCE_ATTACHMENT_KIND.into(),
            name: "original-private.sbol.jsonld".into(),
            media_type: "application/ld+json".into(),
            bytes: raw,
        },
    );
    let plan = DesignPlan {
        name: "Saved actual PCR from SBOL projected stock".into(),
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
    let commit = cloning_outputs::finalize(
        &w,
        &plan,
        event(&w, 1),
        vec![InputEventBinding {
            input_id: "stock".into(),
            event: event(&w, 0),
        }],
        &cloning_outputs::Options {
            name: "Saved native PCR".into(),
            show_prepared: false,
            save_primers: false,
            allow_incomplete: false,
        },
    )
    .unwrap();
    let native = run_material::Selection {
        run_id: commit.receipt.run.run_id.clone(),
        step_index: 0,
        event: commit.receipt.construct.clone(),
    };
    let selection = Selection {
        document_id: native.event.document_id.clone(),
        revision: native.event.revision.clone(),
        head_event: Some(native.event.event_id.clone()),
    };
    (commit.workspace, selection, native, hash)
}
fn import(repository: &mut Repository, w: &Workspace, pin: &Selection) -> store::Document {
    match repository
        .mutate_local(durable::MutationOperation::ImportDocument {
            folder: None,
            name: "Offline plasmid".into(),
            payload: repo::private_workspace(w, pin).unwrap(),
        })
        .unwrap()
    {
        durable::Applied::Imported(document) => document,
        other => panic!("{other:?}"),
    }
}
#[test]
fn approved_subset_survives_restart_lost_reply_and_fresh_native_material_replay() {
    let disk = Disk::new();
    let mut repository = Repository::open(&disk.0).unwrap();
    let mut remote = NativeRemote::new();
    let config = remote.config(&mut repository);
    let (w, pin, origin, original) = fixture(false);
    let document = import(&mut repository, &w, &pin);
    let before = repository.open_document(&document.id, None).unwrap();
    assert!(matches!(
        repository.queue_push(&document.id, &config.id),
        Err(repo::Error::SbolSource)
    ));
    let checkpoint = repository.checkpoint_version();
    let review = repository
        .prepare_push_review(
            &document.id,
            &document.head.id,
            &config.id,
            &PublicationChoices::new(),
        )
        .unwrap();
    assert_eq!(review.local_revision(), &document.head.id);
    assert_eq!(review.remote(), &config);
    assert_eq!(review.checkpoint(), &checkpoint);
    assert_eq!(review.review().entries().len(), 1);
    assert_eq!(review.review().entries()[0].omissions.unrelated_objects, 1);
    assert_eq!(repository.checkpoint_version(), checkpoint);
    assert_eq!(repository.pending_ids().len(), 0);
    let approved = review.approve();
    let job = repository.queue_reviewed_push(&approved).unwrap();
    assert!(matches!(
        repository.queue_reviewed_push(&approved),
        Err(repo::Error::Stale)
    ));
    assert_eq!(
        repository.open_document(&document.id, None).unwrap().bytes,
        before.bytes
    );
    drop(repository);
    let mut repository = Repository::open(&disk.0).unwrap();
    remote.lose_reply = true;
    assert!(repository.resume_push(&job, &mut remote).is_err());
    assert_eq!(remote.applies, 1);
    drop(repository);
    let mut repository = Repository::open(&disk.0).unwrap();
    assert!(matches!(
        repository.resume_push(&job, &mut remote).unwrap(),
        SyncResult::Applied(_)
    ));
    assert_eq!(remote.applies, 1);
    assert_eq!(repository.pending_ids().len(), 0);
    let binding = repository
        .document_binding(&document.id, &config.id)
        .unwrap();
    let remote_id = store::DocumentId::parse(binding.remote.id.as_str()).unwrap();
    let remote_head = remote.db.document(&remote.access, &remote_id).unwrap();
    let published = remote
        .db
        .payload(&remote.access, &remote_id, &remote_head.head.id)
        .unwrap();
    let frame = workspace::load_workspace(&published.bytes).unwrap();
    assert!(!frame.attachments.contains_key(&original));
    assert!(frame.attachments.values().any(|a| a.kind == RECEIPT_KIND));
    assert!(frame.attachments.values().any(|a| a.kind == SUBSET_KIND));
    for attachment in frame.attachments.values() {
        assert!(
            !String::from_utf8_lossy(&attachment.bytes).contains("DO-NOT-PUBLISH-PRIVATE-SENTINEL")
        );
    }
    assert_eq!(frame.revisions[&pin.revision], w.revisions[&pin.revision]);
    assert_eq!(
        run_material::reconstruct(&frame, &origin)
            .unwrap()
            .reference()
            .sequence
            .as_str(),
        DNA
    );
    let local = repository.open_document(&document.id, None).unwrap();
    assert_eq!(local.bytes, before.bytes);
    let local = workspace::load_workspace(&local.bytes).unwrap();
    assert!(local.attachments.contains_key(&original));
    assert!(
        local
            .documents
            .iter()
            .any(|d| local.revisions[&d.head].name == "UNRELATED PRIVATE PLASMID")
    );
}
#[test]
fn review_guards_saved_revision_metadata_checkpoint_and_remote_binding() {
    let disk = Disk::new();
    let mut repository = Repository::open(&disk.0).unwrap();
    let remote = NativeRemote::new();
    let config = remote.config(&mut repository);
    let (w, pin, _, _) = fixture(false);
    let document = import(&mut repository, &w, &pin);
    let review = repository
        .prepare_push_review(
            &document.id,
            &document.head.id,
            &config.id,
            &PublicationChoices::new(),
        )
        .unwrap();
    let mut concurrent = Repository::open(&disk.0).unwrap();
    concurrent
        .mutate_local(durable::MutationOperation::RenameDocument {
            id: document.id.clone(),
            expected_generation: document.generation,
            name: "Changed offline name".into(),
        })
        .unwrap();
    assert!(matches!(
        repository.queue_reviewed_push(&review.approve()),
        Err(repo::Error::Stale)
    ));
    assert_eq!(repository.pending_ids().len(), 0);
    let renamed = repository.document(&document.id).unwrap();
    let review = repository
        .prepare_push_review(
            &renamed.id,
            &renamed.head.id,
            &config.id,
            &PublicationChoices::new(),
        )
        .unwrap();
    let next = match repository
        .mutate_local(durable::MutationOperation::SaveRevision {
            id: renamed.id.clone(),
            base_revision: renamed.head.id.clone(),
            payload: repo::private_workspace(&w, &pin).unwrap(),
        })
        .unwrap()
    {
        durable::Applied::Saved { document, .. } => document,
        other => panic!("{other:?}"),
    };
    assert_ne!(next.head.id, renamed.head.id);
    assert!(matches!(
        repository.queue_reviewed_push(&review.approve()),
        Err(repo::Error::Stale)
    ));
    assert!(matches!(
        repository.prepare_push_review(
            &document.id,
            &document.head.id,
            &config.id,
            &PublicationChoices::new()
        ),
        Err(repo::Error::Stale)
    ));
}
#[test]
fn explicit_alternative_documentary_update_can_publish_without_changing_native_head() {
    let disk = Disk::new();
    let mut repository = Repository::open(&disk.0).unwrap();
    let mut remote = NativeRemote::new();
    let config = remote.config(&mut repository);
    let (w, pin, _, original) = fixture(true);
    let document = import(&mut repository, &w, &pin);
    let requirements = repository
        .push_review_choices(&document.id, &document.head.id, &config.id)
        .unwrap();
    assert_eq!(requirements.len(), 1);
    assert_eq!(
        requirements[0].required_sequences,
        BTreeSet::from([SEQUENCE.to_string()])
    );
    let first = repository
        .prepare_push_review(
            &document.id,
            &document.head.id,
            &config.id,
            &PublicationChoices::new(),
        )
        .unwrap();
    let job = repository.queue_reviewed_push(&first.approve()).unwrap();
    repository.resume_push(&job, &mut remote).unwrap();
    let before = repository
        .document_binding(&document.id, &config.id)
        .unwrap()
        .remote
        .clone();
    let choices = BTreeMap::from([(
        original,
        BTreeMap::from([(
            COMPONENT.into(),
            BTreeSet::from([SEQUENCE.into(), ALTERNATIVE.into()]),
        )]),
    )]);
    let review = repository
        .prepare_push_review(&document.id, &document.head.id, &config.id, &choices)
        .unwrap();
    assert_eq!(
        review.review().entries()[0]
            .omissions
            .alternative_sequence_edges,
        0
    );
    let job = repository.queue_reviewed_push(&review.approve()).unwrap();
    assert!(matches!(
        repository.resume_push(&job, &mut remote).unwrap(),
        SyncResult::Applied(_)
    ));
    let after = repository
        .document_binding(&document.id, &config.id)
        .unwrap()
        .remote
        .clone();
    assert_eq!(before.id, after.id);
    assert_eq!(after.head.number, before.head.number + 1);
    assert_ne!(after.head.archive_hash, before.head.archive_hash);
    assert_eq!(
        repository.document(&document.id).unwrap().head.id,
        document.head.id
    );
    let duplicate = repository
        .prepare_push_review(&document.id, &document.head.id, &config.id, &choices)
        .unwrap();
    assert!(matches!(
        repository.queue_reviewed_push(&duplicate.approve()),
        Err(repo::Error::NoChanges)
    ));
}
