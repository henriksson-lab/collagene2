use collagene_core::{
    documents::{
        history::DocumentHistory,
        rna_history as native,
        workspace::{self, Workspace, WorkspaceAttachment},
    },
    material::lineage::EventBinding,
    rna_chemistry as rna,
};
use collagene_model::{Alphabet, Molecule, Sequence, Topology};
fn fixture() -> (Workspace, native::Record) {
    let m = Molecule {
        name: "RNA".into(),
        sequence: Sequence::new(Alphabet::Rna, "AC").unwrap(),
        topology: Topology::Linear,
        features: vec![],
        metadata: Default::default(),
        original_genbank: None,
    };
    let w = Workspace::from_documents("RNA records", vec![m.clone(), m], 0).unwrap();
    let d = &w.documents[0];
    let r = native::Record {
        version: native::VERSION,
        name: "Intended RNA model".into(),
        source: EventBinding {
            document_id: d.id.clone(),
            event_id: d.head_event.clone().unwrap(),
            revision: d.head.clone(),
        },
        declaration: rna::Declaration {
            version: rna::DECLARATION_VERSION,
            source_revision: d.head.clone(),
            coverage: rna::Coverage::ClosedModel,
            ends: rna::Ends::Linear {
                five: rna::FivePrime::Hydroxyl,
                three: rna::ThreePrime::Hydroxyl,
            },
            sites: vec![],
            links: vec![],
            annotation_reviews: vec![],
            provenance: rna::Provenance {
                statement: "Declared fixture hypothesis".into(),
                reference: None,
            },
        },
    };
    (w, r)
}
#[test]
fn declarations_roundtrip_exact_historical_events_without_cached_quantities_or_source_changes() {
    let (w, r) = fixture();
    let before = serde_json::to_value(&w).unwrap();
    let (saved, hash) = native::capture(&w, &r).unwrap();
    assert_eq!(serde_json::to_value(&w).unwrap(), before);
    let body: serde_json::Value = serde_json::from_slice(&saved.attachments[&hash].bytes).unwrap();
    assert_eq!(body.as_object().unwrap().len(), 4);
    assert!(body.get("report").is_none());
    assert_eq!(native::capture(&saved, &r).unwrap().0.attachments.len(), 1);
    let loaded = workspace::load_workspace(&workspace::save_workspace(&saved).unwrap()).unwrap();
    let mut histories = DocumentHistory::restore_all(&loaded).unwrap();
    let mut changed = histories[0].molecule.clone();
    changed.name = "Later edit".into();
    histories[0].commit(changed, "Change display name").unwrap();
    let mut advanced = loaded.clone();
    advanced.documents.clear();
    advanced.history.clear();
    for d in histories {
        d.contribute(&mut advanced).unwrap();
    }
    advanced.validate().unwrap();
    let reopened = native::load(&advanced, &hash).unwrap();
    assert_eq!(reopened.source, r.source);
    assert_ne!(advanced.documents[0].head, reopened.source.revision);
    assert_eq!(native::list(&advanced).unwrap()[0].source_name, "RNA");
    let source = native::source(&advanced, &reopened.source).unwrap();
    assert_eq!(
        rna::analyze(source, &reopened.declaration)
            .unwrap()
            .formula
            .unwrap()
            .hill(),
        "C19H25N8O11P"
    );
    let deleted = native::delete(&advanced, &hash).unwrap();
    assert!(deleted.attachments.is_empty());
    assert_eq!(deleted.documents[0].head, advanced.documents[0].head);
}
fn inject(w: &Workspace, value: &serde_json::Value) -> (Workspace, String) {
    let bytes = serde_json::to_vec(value).unwrap();
    let hash = blake3::hash(&bytes).to_hex().to_string();
    let mut next = w.clone();
    next.attachments.insert(
        hash.clone(),
        WorkspaceAttachment {
            kind: native::KIND.into(),
            name: "Injected fixture".into(),
            media_type: "application/json".into(),
            bytes,
        },
    );
    (next, hash)
}
#[test]
fn current_archive_refuses_foreign_equal_content_events_corrupt_roles_and_noncurrent_records() {
    let (w, r) = fixture();
    let mut wrong = r.clone();
    wrong.source.document_id = w.documents[1].id.clone();
    assert!(native::capture(&w, &wrong).is_err());
    wrong = r.clone();
    wrong.declaration.source_revision = "f".repeat(64);
    assert!(native::capture(&w, &wrong).is_err());
    let value = serde_json::to_value(&r).unwrap();
    for version in [0, 2, 65535] {
        let mut invalid = value.clone();
        invalid["version"] = version.into();
        let (bad, hash) = inject(&w, &invalid);
        assert!(bad.validate().is_err());
        assert!(native::load(&bad, &hash).is_err());
    }
    let mut invalid = value.clone();
    invalid["declaration"]["version"] = 1.into();
    assert!(inject(&w, &invalid).0.validate().is_err());
    let mut invalid = value;
    invalid["report"] = serde_json::json!({"mass":123});
    assert!(inject(&w, &invalid).0.validate().is_err());
    let (mut saved, hash) = native::capture(&w, &r).unwrap();
    saved.attachments.get_mut(&hash).unwrap().kind = "other-documentary-role".into();
    assert!(native::capture(&saved, &r).is_err());
    assert!(native::load(&saved, &hash).is_err());
    let mut wrong = r;
    wrong.name = "x".repeat(257);
    assert!(native::capture(&w, &wrong).is_err());
}
