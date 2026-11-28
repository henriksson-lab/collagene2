use collagene_core::documents::reference_edit::Operation as Edit;
use collagene_core::{
    alignment::{
        AlignmentOptions,
        history::{self as pairwise_history, Record},
        review::{self as pairwise_review, Options, Report},
    },
    documents::{
        history::DocumentHistory,
        workspace::{Workspace, WorkspaceAttachment},
    },
    material::lineage::EventBinding,
};
use collagene_model::{Molecule, Topology};

fn binding(frame: &Workspace, index: usize) -> EventBinding {
    let doc = &frame.documents[index];
    EventBinding {
        document_id: doc.id.clone(),
        event_id: doc.head_event.clone().unwrap(),
        revision: doc.head.clone(),
    }
}
fn fixture() -> (Workspace, EventBinding, EventBinding, Report) {
    let reference = Molecule::dna("reference", "ACGTTAGC", Topology::Linear).unwrap();
    let query = Molecule::dna("query", "ACGTTCGC", Topology::Linear).unwrap();
    let report = pairwise_review::compare(
        &reference,
        &query,
        &Options::Stored(AlignmentOptions::default()),
    )
    .unwrap();
    let frame = Workspace::from_documents(
        "retained pairwise comparison",
        vec![reference.clone(), reference, query],
        0,
    )
    .unwrap();
    let reference = binding(&frame, 0);
    let query = binding(&frame, 2);
    (frame, reference, query, report)
}
fn advance_reference(frame: &Workspace) -> Workspace {
    let mut docs = DocumentHistory::restore_all(frame).unwrap();
    let edit_end = docs[0].molecule.sequence.len();
    docs[0]
        .commit(
            vec![
                Edit::ReplaceSequence {
                    start: 0,
                    end: edit_end,
                    replacement: "TCGTTAGC".into(),
                },
                Edit::Rename {
                    name: "later reference".into(),
                },
                Edit::SetTopology {
                    topology: Topology::Linear,
                },
            ],
            "Later reference edit",
        )
        .unwrap();
    let mut next = frame.clone();
    next.documents.clear();
    next.history.clear();
    for doc in docs {
        doc.contribute(&mut next).unwrap();
    }
    next.validate().unwrap();
    next
}
fn insert_record(frame: &Workspace, value: &serde_json::Value) -> (Workspace, String) {
    let bytes = serde_json::to_vec(value).unwrap();
    let hash = blake3::hash(&bytes).to_hex().to_string();
    let mut next = frame.clone();
    next.attachments.insert(
        hash.clone(),
        WorkspaceAttachment {
            kind: pairwise_history::KIND.into(),
            name: "modified recorded evidence".into(),
            media_type: "application/json".into(),
            bytes,
        },
    );
    (next, hash)
}

#[test]
fn exact_event_capture_reopens_historically_and_deduplicates_without_mutating_documents() {
    let (frame, reference, query, report) = fixture();
    let record = Record {
        schema_version: pairwise_history::VERSION,
        reference: reference.clone(),
        query: query.clone(),
        report: report.clone(),
        design_consensus: None,
    };
    let before = serde_json::to_value(&frame).unwrap();
    let (saved, hash) = pairwise_history::capture(&frame, &record).unwrap();
    assert_eq!(serde_json::to_value(&frame).unwrap(), before);
    assert_eq!(saved.attachments.len(), frame.attachments.len() + 1);
    let mut without_record = saved.clone();
    without_record.attachments.remove(&hash);
    assert_eq!(serde_json::to_value(without_record).unwrap(), before);
    let reopened = pairwise_history::load(&saved, &hash).unwrap();
    assert_eq!(reopened.schema_version, pairwise_history::VERSION);
    assert_eq!(reopened.reference, reference);
    assert_eq!(reopened.query, query);
    assert_eq!(reopened.report.comparison.alignment.differences.len(), 1);
    let changed = &reopened.report.comparison.alignment.differences[0];
    assert_eq!((changed.reference_start, changed.reference_end), (5, 6));
    assert_eq!(changed.reference_bases, "A");
    assert_eq!(changed.query_bases, "C");
    // Record identity carries no persisted session freshness or verification bit.
    let fields = serde_json::to_value(&reopened).unwrap();
    assert_eq!(fields.as_object().unwrap().len(), 5);
    assert!(fields.get("design_consensus").unwrap().is_null());
    assert!(fields.get("fresh").is_none());
    assert!(fields.get("verified").is_none());
    let (retry, retry_hash) = pairwise_history::capture(&saved, &record).unwrap();
    assert_eq!(retry_hash, hash);
    assert_eq!(
        serde_json::to_value(retry).unwrap(),
        serde_json::to_value(&saved).unwrap()
    );
    let listed = pairwise_history::list(&saved).unwrap();
    assert_eq!(listed.len(), 1);
    assert_eq!(listed[0].attachment_hash, hash);
    assert_eq!(listed[0].reference_name, "reference");
    assert_eq!(listed[0].query_name, "query");
}

#[test]
fn historical_bindings_survive_head_changes_but_equal_content_does_not_guess_an_event() {
    let (frame, reference, query, report) = fixture();
    let record = Record {
        schema_version: pairwise_history::VERSION,
        reference: reference.clone(),
        query: query.clone(),
        report: report.clone(),
        design_consensus: None,
    };
    let (saved, hash) = pairwise_history::capture(&frame, &record).unwrap();
    let later = advance_reference(&saved);
    assert_ne!(later.documents[0].head, reference.revision);
    assert_eq!(
        pairwise_history::load(&later, &hash).unwrap().reference,
        reference
    );
    let (retained, retained_hash) = pairwise_history::capture(&later, &record).unwrap();
    assert_eq!(retained_hash, hash);
    assert_eq!(retained.documents[0].head, later.documents[0].head);
    let before = serde_json::to_value(&later).unwrap();
    assert!(
        pairwise_history::capture(
            &later,
            &Record {
                reference: binding(&later, 0),
                ..record.clone()
            }
        )
        .is_err()
    );
    let mut spoof = reference.clone();
    spoof.document_id = later.documents[1].id.clone();
    assert_eq!(later.documents[1].head, reference.revision);
    assert!(
        pairwise_history::capture(
            &later,
            &Record {
                reference: spoof,
                ..record.clone()
            }
        )
        .is_err()
    );
    assert_eq!(serde_json::to_value(&later).unwrap(), before);

    // An explicitly chosen equal-content document is a distinct, valid binding;
    // it receives a distinct evidence identity instead of replacing the original.
    let alias = binding(&later, 1);
    let (both, alias_hash) = pairwise_history::capture(
        &later,
        &Record {
            reference: alias.clone(),
            ..record.clone()
        },
    )
    .unwrap();
    assert_ne!(alias_hash, hash);
    assert_eq!(
        pairwise_history::load(&both, &hash).unwrap().reference,
        reference
    );
    assert_eq!(
        pairwise_history::load(&both, &alias_hash)
            .unwrap()
            .reference,
        alias
    );
    assert_eq!(pairwise_history::list(&both).unwrap().len(), 2);
}

#[test]
fn tampered_bindings_scores_versions_or_attachment_roles_fail_without_publishing_partial_state() {
    let (frame, reference, query, report) = fixture();
    let record = Record {
        schema_version: pairwise_history::VERSION,
        reference: reference.clone(),
        query: query.clone(),
        report: report.clone(),
        design_consensus: None,
    };
    let before = serde_json::to_value(&frame).unwrap();
    assert!(
        pairwise_history::capture(
            &frame,
            &Record {
                reference: query.clone(),
                query: reference.clone(),
                ..record.clone()
            }
        )
        .is_err()
    );
    let mut wrong_score = report.clone();
    wrong_score.comparison.alignment.score += 1;
    assert!(
        pairwise_history::capture(
            &frame,
            &Record {
                report: wrong_score,
                ..record.clone()
            }
        )
        .is_err()
    );
    assert_eq!(serde_json::to_value(&frame).unwrap(), before);
    let (saved, hash) = pairwise_history::capture(&frame, &record).unwrap();
    let original: serde_json::Value =
        serde_json::from_slice(&saved.attachments[&hash].bytes).unwrap();
    let mut changed = original.clone();
    changed["reference"]["document_id"] = saved.documents[1].id.clone().into();
    let (spoofed, spoofed_hash) = insert_record(&saved, &changed);
    assert!(pairwise_history::load(&spoofed, &spoofed_hash).is_err());
    assert!(pairwise_history::list(&spoofed).is_err());
    changed = original.clone();
    changed["schema_version"] = (pairwise_history::VERSION + 1).into();
    let (version, version_hash) = insert_record(&saved, &changed);
    assert!(pairwise_history::load(&version, &version_hash).is_err());
    let mut missing_semantics = original.clone();
    missing_semantics
        .as_object_mut()
        .unwrap()
        .remove("design_consensus");
    let (missing_semantics, missing_hash) = insert_record(&saved, &missing_semantics);
    assert!(pairwise_history::load(&missing_semantics, &missing_hash).is_err());
    changed = original;
    changed["verified"] = true.into();
    let (extra, extra_hash) = insert_record(&saved, &changed);
    assert!(pairwise_history::load(&extra, &extra_hash).is_err());
    let mut wrong_kind = saved.clone();
    wrong_kind.attachments.get_mut(&hash).unwrap().kind = "other-source".into();
    assert!(pairwise_history::load(&wrong_kind, &hash).is_err());
    let mut wrong_media = saved.clone();
    wrong_media.attachments.get_mut(&hash).unwrap().media_type = "application/octet-stream".into();
    assert!(pairwise_history::load(&wrong_media, &hash).is_err());
    let mut wrong_hash = saved.clone();
    wrong_hash
        .attachments
        .get_mut(&hash)
        .unwrap()
        .bytes
        .push(b' ');
    assert!(pairwise_history::load(&wrong_hash, &hash).is_err());

    let record: Record = pairwise_history::load(&saved, &hash).unwrap();
    assert_eq!(record.reference.event_id, reference.event_id);
    let mut missing = saved.clone();
    missing
        .history
        .retain(|event| event.id != reference.event_id);
    assert!(pairwise_history::load(&missing, &hash).is_err());
}
