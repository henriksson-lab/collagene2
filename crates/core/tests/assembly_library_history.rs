//! Recovery is raw editable data, atomically retained with one active draft.
use collagene_core::{
    Feature, Location, Molecule, Topology,
    assembly_library_history::{self as history, Impact},
    run_lineage::{self, EventBinding},
    workspace::{RevisionEvent, Workspace, WorkspaceAttachment},
};
fn hash(bytes: &[u8]) -> String {
    blake3::hash(bytes).to_hex().to_string()
}
// Deliberately only a probe-valid header: core does not own WebDraft semantics.
fn draft(counter: u64, text: &str) -> Vec<u8> {
    serde_json::to_vec(&serde_json::json!({"version":3,"alternatives":{"id":"literal-library","next_id":counter},"working":text})).unwrap()
}
fn fixture() -> Workspace {
    let mut w = Workspace::from_documents(
        "library recovery",
        vec![
            Molecule::dna("original", "AACC", Topology::Linear).unwrap(),
            Molecule::dna("original", "AACC", Topology::Linear).unwrap(),
        ],
        0,
    )
    .unwrap();
    let bytes = draft(8, "persisted saved-row text");
    w.attachments.insert(
        hash(&bytes),
        WorkspaceAttachment {
            kind: history::ACTIVE_DRAFT_KIND.into(),
            name: "assembly-draft.json".into(),
            media_type: "application/json".into(),
            bytes,
        },
    );
    w
}
fn binding(w: &Workspace, index: usize) -> EventBinding {
    let d = &w.documents[index];
    EventBinding {
        document_id: d.id.clone(),
        event_id: d.head_event.clone().unwrap(),
        revision: d.head.clone(),
    }
}
fn impact() -> Impact {
    Impact {
        construct_id: "construct-7".into(),
        name: "Literal row".into(),
        primer_review_bins: vec!["bin-2".into()],
        removed_preparation_bins: vec!["bin-3".into()],
        retained_junctions: 1,
        recovered_junctions: 2,
        new_junctions: vec![0],
        overlap_text_changed: true,
        role_order_changed: true,
        candidate_cleared: true,
        construct_removed: false,
        construct_added: false,
    }
}
fn proposal() -> String {
    hash(b"explicit reviewed edit")
}
#[test]
fn exact_unsaved_before_bytes_survive_atomic_replacement_and_native_archive() {
    let w = fixture();
    let snapshot = serde_json::to_vec(&w).unwrap();
    let before=b"{\n  \"working\": \"invalid overlap text: 12,,7; intentional primer edit\",\n  \"version\":3,\"alternatives\":{\"id\":\"literal-library\",\"next_id\":9}\n}".to_vec();
    let after = draft(12, "new reviewed order");
    let pin = binding(&w, 0);
    let record = history::create(
        &w,
        before.clone(),
        &after,
        &proposal(),
        vec![pin.clone(), pin],
        vec![impact()],
    )
    .unwrap();
    let mut extended = history::encode_record(&record).unwrap();
    extended.push(0);
    assert!(history::load_record(&extended).is_err());
    assert_ne!(
        record.before_hash, record.before_active_hash,
        "persisted canonical row is not the unsaved working recovery"
    );
    assert_eq!(record.before_next_id, 9);
    assert_eq!(record.after_next_id, 12);
    assert_eq!(
        record.sources.len(),
        1,
        "duplicate exact pins are normalized, not content-selected"
    );
    let commit = history::install(w.clone(), record, after.clone()).unwrap();
    assert_eq!(serde_json::to_vec(&w).unwrap(), snapshot);
    let active: Vec<_> = commit
        .workspace
        .attachments
        .values()
        .filter(|a| a.kind == history::ACTIVE_DRAFT_KIND)
        .collect();
    assert_eq!(active.len(), 1);
    assert_eq!(active[0].bytes, after);
    let saved = history::load(&commit.workspace, &commit.record_hash).unwrap();
    assert_eq!(saved.before_draft, before);
    assert!(saved.impacts[0].candidate_cleared);
    assert_eq!(saved.impacts[0].removed_preparation_bins, vec!["bin-3"]);
    assert_eq!(
        history::list(&commit.workspace).unwrap()[0].allocation_highwater,
        12
    );
    let archive = collagene_core::workspace::save_workspace(&commit.workspace).unwrap();
    let imported = collagene_core::workspace::load_workspace(&archive).unwrap();
    assert_eq!(
        history::load(&imported, &commit.record_hash)
            .unwrap()
            .before_draft,
        before
    );
    assert!(imported.candidates.is_empty());
    assert!(
        imported
            .attachments
            .values()
            .all(|a| a.kind != run_lineage::WORKSPACE_ATTACHMENT_KIND)
    );
}
#[test]
fn stale_install_replacement_tampering_or_ambiguous_active_drafts_publish_nothing() {
    let w = fixture();
    let after = draft(9, "after");
    let record = history::create(
        &w,
        draft(8, "before"),
        &after,
        &proposal(),
        vec![binding(&w, 0)],
        vec![impact()],
    )
    .unwrap();
    let mut changed = w.clone();
    changed.name = "unrelated concurrent edit".into();
    assert!(history::install(changed, record.clone(), after.clone()).is_err());
    assert!(history::install(w.clone(), record.clone(), draft(9, "tampered after")).is_err());
    let mut forged = record;
    forged.before_draft = draft(8, "rewritten raw before");
    assert!(history::install(w.clone(), forged, after.clone()).is_err());
    let mut ambiguous = w.clone();
    let legacy = b"legacy".to_vec();
    ambiguous.attachments.insert(
        hash(&legacy),
        WorkspaceAttachment {
            kind: "assembly-draft-v2".into(),
            name: "second active".into(),
            media_type: "application/json".into(),
            bytes: legacy,
        },
    );
    assert!(
        history::create(
            &ambiguous,
            draft(8, "before"),
            &after,
            &proposal(),
            vec![],
            vec![]
        )
        .is_err()
    );
    assert_eq!(w.attachments.len(), 1);
}
#[test]
fn historical_exact_sources_survive_later_heads_but_equal_dna_never_proves_membership() {
    let w = fixture();
    let first = binding(&w, 0);
    let other = binding(&w, 1);
    assert_eq!(first.revision, other.revision);
    let record = history::create(
        &w,
        draft(8, "before"),
        &draft(9, "after"),
        &proposal(),
        vec![first.clone(), other.clone()],
        vec![],
    )
    .unwrap();
    assert_eq!(
        record.sources.len(),
        2,
        "equal DNA distinct documents remain separate pins"
    );
    let commit = history::install(w, record, draft(9, "after")).unwrap();
    let mut later = commit.workspace;
    let molecule = Molecule::dna("edited later", "GGCC", Topology::Linear).unwrap();
    let rev = molecule.revision();
    let event = RevisionEvent::new(
        first.document_id.clone(),
        rev.clone(),
        vec![first.event_id.clone()],
        "later source edit",
    )
    .unwrap();
    later.revisions.insert(rev.clone(), molecule);
    let d = &mut later.documents[0];
    d.undo.push(d.head.clone());
    d.undo_events.push(d.head_event.clone().unwrap());
    d.head = rev;
    d.head_event = Some(event.id.clone());
    later.history.push(event);
    later.validate().unwrap();
    let historical = history::load(&later, &commit.record_hash).unwrap();
    assert!(
        historical
            .sources
            .iter()
            .any(|b| b.event_id == first.event_id)
    );
    assert_ne!(
        historical.source_workspace_fingerprint,
        run_lineage::workspace_fingerprint(&later).unwrap()
    );
    let mut mismatched = other;
    mismatched.document_id = "not-retained".into();
    assert!(
        history::create(
            &later,
            draft(9, "old"),
            &draft(10, "new"),
            &proposal(),
            vec![mismatched],
            vec![]
        )
        .is_err()
    );
    let mut tampered = historical;
    tampered.sources[0].document_id = "foreign-document".into();
    let bytes = history::encode_record(&tampered).unwrap();
    let key = hash(&bytes);
    later.attachments.insert(
        key.clone(),
        WorkspaceAttachment {
            kind: history::KIND.into(),
            name: "forged declared pins".into(),
            media_type: history::MEDIA_TYPE.into(),
            bytes,
        },
    );
    assert!(history::load(&later, &key).is_err());
}
#[test]
fn recovery_preserves_allocation_highwater_and_rejects_header_regression() {
    let w = fixture();
    let after = draft(40, "new allocated choices");
    let record =
        history::create(&w, draft(8, "old raw"), &after, &proposal(), vec![], vec![]).unwrap();
    let commit = history::install(w, record, after).unwrap();
    assert_eq!(
        history::allocation_highwater(&commit.workspace, "literal-library").unwrap(),
        40
    );
    assert!(
        history::create(
            &commit.workspace,
            draft(8, "recovered old working draft"),
            &draft(9, "unsafe old counter"),
            &proposal(),
            vec![],
            vec![]
        )
        .is_err()
    );
    let restored = draft(40, "recovered exact working text with counter raised");
    let record = history::create(
        &commit.workspace,
        draft(40, "current"),
        &restored,
        &hash(b"restore proposal"),
        vec![],
        vec![],
    )
    .unwrap();
    let restored = history::install(commit.workspace, record, restored).unwrap();
    assert_eq!(history::list(&restored.workspace).unwrap().len(), 2);
    assert_eq!(
        history::allocation_highwater(&restored.workspace, "literal-library").unwrap(),
        40
    );
    assert!(
        history::create(
            &restored.workspace,
            draft(40, "current"),
            br#"{"version":3,"alternatives":{"id":"other-library","next_id":41}}"#,
            &proposal(),
            vec![],
            vec![]
        )
        .is_err()
    );
    assert!(
        history::create(
            &restored.workspace,
            draft(40, "current"),
            br#"{"version":2,"alternatives":{"id":"literal-library","next_id":41}}"#,
            &proposal(),
            vec![],
            vec![]
        )
        .is_err()
    );
}
#[test]
fn rich_native_sources_and_restoration_union_impacts_are_not_scientific_profiles() {
    let mut w = fixture();
    let old = w.documents[0].head.clone();
    let mut rich = w.revisions[&old].clone();
    rich.features = (0..10_000)
        .map(|_| Feature {
            kind: "misc_feature".into(),
            location: Location::range(0, 1),
            qualifiers: vec![],
            needs_review: false,
        })
        .collect();
    let rich = Workspace::from_documents("rich untouched", vec![rich], 0).unwrap();
    // Use a separately valid rich frame; source pinning is documentary native work.
    w = rich;
    let bytes = draft(8, "raw");
    w.attachments.insert(
        hash(&bytes),
        WorkspaceAttachment {
            kind: history::ACTIVE_DRAFT_KIND.into(),
            name: "assembly-draft.json".into(),
            media_type: "application/json".into(),
            bytes,
        },
    );
    let impacts = (0..history::MAX_IMPACTS)
        .map(|i| {
            let mut impact = impact();
            impact.construct_id = format!("row-{i}");
            impact.construct_removed = i < 256;
            impact.construct_added = i >= 256;
            impact
        })
        .collect();
    let record = history::create(
        &w,
        draft(8, "raw"),
        &draft(9, "after"),
        &proposal(),
        vec![binding(&w, 0)],
        impacts,
    )
    .unwrap();
    assert_eq!(record.impacts.len(), 512);
    let mut value = serde_json::to_value(&record).unwrap();
    let impacts = value["impacts"].as_array_mut().unwrap();
    impacts.push(impacts[0].clone());
    assert!(serde_json::from_value::<history::Record>(value).is_err());
    let committed = history::install(w, record, draft(9, "after")).unwrap();
    assert_eq!(
        history::load(&committed.workspace, &committed.record_hash)
            .unwrap()
            .impacts
            .len(),
        512
    );
}
#[test]
fn exhausted_allocation_counter_still_supports_nonallocating_recovery_edits() {
    let w = fixture();
    let maximum = draft(u64::MAX, "last allocation already consumed");
    let record = history::create(
        &w,
        draft(u64::MAX - 1, "before last allocation"),
        &maximum,
        &proposal(),
        vec![],
        vec![impact()],
    )
    .unwrap();
    assert_eq!(record.allocation_highwater, u64::MAX);
    let committed = history::install(w, record, maximum).unwrap();
    let reordered = draft(u64::MAX, "reordered existing bins without allocating");
    let record = history::create(
        &committed.workspace,
        draft(u64::MAX, "exact exhausted working draft"),
        &reordered,
        &hash(b"reorder at exhausted counter"),
        vec![],
        vec![impact()],
    )
    .unwrap();
    let committed = history::install(committed.workspace, record, reordered).unwrap();
    let summaries = history::list(&committed.workspace).unwrap();
    assert_eq!(summaries.len(), 2);
    assert!(
        summaries
            .iter()
            .all(|s| s.after_next_id == u64::MAX && s.allocation_highwater == u64::MAX)
    );
    assert_eq!(
        history::allocation_highwater(&committed.workspace, "literal-library").unwrap(),
        u64::MAX
    );
    let recovered = history::load(&committed.workspace, &committed.record_hash).unwrap();
    assert_eq!(recovered.before_next_id, u64::MAX);
    assert_eq!(recovered.after_next_id, u64::MAX);
    assert!(
        history::create(
            &committed.workspace,
            draft(u64::MAX, "current"),
            &draft(u64::MAX - 1, "counter rewind"),
            &proposal(),
            vec![],
            vec![]
        )
        .is_err()
    );
}
#[test]
fn record_retention_count_exhaustion_refuses_the_whole_edit_without_pruning() {
    let mut w = fixture();
    let after = draft(9, "after");
    let template =
        history::create(&w, draft(8, "before"), &after, &proposal(), vec![], vec![]).unwrap();
    for i in 0..history::MAX_RECORDS {
        let mut record = template.clone();
        record.proposal_revision = hash(format!("past-{i}").as_bytes());
        let bytes = history::encode_record(&record).unwrap();
        w.attachments.insert(
            hash(&bytes),
            WorkspaceAttachment {
                kind: history::KIND.into(),
                name: format!("past-{i}"),
                media_type: history::MEDIA_TYPE.into(),
                bytes,
            },
        );
    }
    let snapshot = serde_json::to_vec(&w).unwrap();
    let proposed =
        history::create(&w, draft(8, "before"), &after, &proposal(), vec![], vec![]).unwrap();
    assert!(history::install(w.clone(), proposed, after).is_err());
    assert_eq!(history::list(&w).unwrap().len(), 128);
    assert_eq!(serde_json::to_vec(&w).unwrap(), snapshot);
}

#[test]
fn newer_ordinary_drafts_cannot_be_ignored_as_second_active_library_drafts() {
    let original = fixture();
    for version in 5..=10 {
        let mut w = original.clone();
        let raw = format!("ordinary editable draft {version}").into_bytes();
        w.attachments.insert(
            hash(&raw),
            WorkspaceAttachment {
                kind: format!("assembly-draft-v{version}"),
                name: "second active draft".into(),
                media_type: "application/json".into(),
                bytes: raw,
            },
        );
        let before = serde_json::to_vec(&w).unwrap();
        assert!(
            history::create(
                &w,
                draft(8, "before"),
                &draft(9, "after"),
                &proposal(),
                vec![],
                vec![]
            )
            .is_err()
        );
        assert!(history::allocation_highwater(&w, "literal-library").is_err());
        assert_eq!(serde_json::to_vec(&w).unwrap(), before);
    }
}
