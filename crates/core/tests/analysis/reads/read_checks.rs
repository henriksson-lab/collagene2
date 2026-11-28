use collagene_core::{
    history::DocumentHistory,
    oligo_formats::{self, QualityEncoding},
    read_checks::*,
    run_lineage::EventBinding,
    trace_validation::{TraceValidationOptions, TrimMode, ValidationStatus},
    workspace::{Workspace, WorkspaceAttachment},
};
use collagene_model::{Molecule, Topology};
fn binding(w: &Workspace, n: usize) -> EventBinding {
    let d = &w.documents[n];
    EventBinding {
        document_id: d.id.clone(),
        event_id: d.head_event.clone().unwrap(),
        revision: d.head.clone(),
    }
}
fn fixture() -> Workspace {
    let reference = Molecule::dna("construct", "TTTCGATACGGG", Topology::Linear).unwrap();
    let read = oligo_formats::import_fastq(b"@read\nCGATAC\n+\nIIIIII\n", QualityEncoding::Phred33)
        .unwrap()
        .molecules
        .remove(0);
    Workspace::from_documents(
        "checks",
        vec![reference.clone(), reference, read.clone(), read],
        0,
    )
    .unwrap()
}
fn options() -> TraceValidationOptions {
    let mut o = TraceValidationOptions::default();
    o.trim.mode = TrimMode::AllCalls;
    o.trim.min_retained_bases = 1;
    o
}
#[test]
fn capture_is_construct_owned_partial_coverage_not_reference_filled_and_retry_dedupes() {
    let w = fixture();
    let before = serde_json::to_value(&w).unwrap();
    let r = binding(&w, 0);
    let reads = vec![binding(&w, 2)];
    let token = capture_revision(&w, &r, &reads, &options()).unwrap();
    let c = capture(&w, r.clone(), reads.clone(), &options()).unwrap();
    c.validate_binding(&token).unwrap();
    assert_eq!(serde_json::to_value(&w).unwrap(), before);
    assert_eq!(c.workspace.documents.len(), 4);
    assert_eq!(c.workspace.history.len(), w.history.len());
    assert_eq!(c.workspace.active_document, w.active_document);
    let s = load(&c.workspace, &c.attachment_hash).unwrap();
    assert_eq!(s.reference, r);
    assert_eq!(s.reads, reads);
    assert_eq!(s.report.status, ValidationStatus::Incomplete);
    assert_eq!(s.report.consensus, "NNNCGATACNNN");
    assert_eq!(s.report.eligible_coverage, 6);
    assert_eq!(s.report.uncovered_ranges, vec![(0, 3), (9, 12)]);
    let retry = capture(&c.workspace, r, reads, &options()).unwrap();
    assert_eq!(retry.attachment_hash, c.attachment_hash);
    assert_eq!(
        retry.workspace.attachments.len(),
        c.workspace.attachments.len()
    );
    assert_ne!(retry.context_revision, c.context_revision);
    let summary = list(&retry.workspace).unwrap();
    assert_eq!(summary.len(), 1);
    assert_eq!(summary[0].read_count, 1);
    assert_eq!(summary[0].recorded_status, ValidationStatus::Incomplete);
    assert!(c.validate_binding(&"0".repeat(64)).is_err());
}
#[test]
fn equal_content_does_not_replace_explicit_reference_or_read_events() {
    let w = fixture();
    assert_eq!(binding(&w, 0).revision, binding(&w, 1).revision);
    assert_ne!(binding(&w, 0).event_id, binding(&w, 1).event_id);
    let first = capture(&w, binding(&w, 0), vec![binding(&w, 2)], &options()).unwrap();
    let alias = capture(
        &first.workspace,
        binding(&w, 1),
        vec![binding(&w, 2)],
        &options(),
    )
    .unwrap();
    assert_ne!(first.attachment_hash, alias.attachment_hash);
    let read_alias = capture(
        &alias.workspace,
        binding(&w, 0),
        vec![binding(&w, 3)],
        &options(),
    )
    .unwrap();
    assert_ne!(first.attachment_hash, read_alias.attachment_hash);
    assert_eq!(list(&read_alias.workspace).unwrap().len(), 3);
    let duplicates = capture(
        &w,
        binding(&w, 0),
        vec![binding(&w, 2), binding(&w, 3)],
        &options(),
    )
    .unwrap();
    let saved = load(&duplicates.workspace, &duplicates.attachment_hash).unwrap();
    assert!(
        saved.report.reads[1]
            .error
            .as_ref()
            .unwrap()
            .contains("duplicate")
    );
    assert_eq!(saved.report.eligible_coverage, 6);
}
#[test]
fn replay_uses_exact_historical_inputs_after_user_head_and_navigation_changes() {
    let w = fixture();
    let reference = binding(&w, 0);
    let read = binding(&w, 2);
    let c = capture(&w, reference.clone(), vec![read.clone()], &options()).unwrap();
    let mut changed = c.workspace.clone();
    let mut histories = DocumentHistory::restore_all(&changed).unwrap();
    histories[0]
        .commit(
            Molecule::dna("user edited construct", "AAAA", Topology::Linear).unwrap(),
            "User edit",
        )
        .unwrap();
    histories[2]
        .commit(
            Molecule::dna("user edited read", "CCCC", Topology::Linear).unwrap(),
            "User edit",
        )
        .unwrap();
    changed.documents.clear();
    changed.history.clear();
    for h in &histories {
        h.contribute(&mut changed).unwrap();
    }
    changed.active_document = read.document_id.clone();
    let checked = replay(&changed, &c.attachment_hash).unwrap();
    checked
        .validate_binding(&replay_revision(&changed, &c.attachment_hash).unwrap())
        .unwrap();
    assert_eq!(checked.saved.reference, reference);
    assert_eq!(checked.saved.reads, vec![read]);
    assert_eq!(checked.report.reference_revision, binding(&w, 0).revision);
    assert_eq!(checked.report.consensus, "NNNCGATACNNN");
    assert_eq!(
        changed
            .head(&reference.document_id)
            .unwrap()
            .sequence
            .as_str(),
        "AAAA"
    );
}
#[test]
fn archived_engine_is_historical_only_and_replay_is_fresh_current_engine() {
    let w = fixture();
    let c = capture(&w, binding(&w, 0), vec![binding(&w, 2)], &options()).unwrap();
    let mut saved = load(&c.workspace, &c.attachment_hash).unwrap();
    saved.report.engine_version = "archived-untrusted-engine".into();
    saved.context_revision = blake3::hash(
        &serde_json::to_vec(&(
            "read-check-evidence-1",
            &saved.reference,
            &saved.reads,
            &saved.report.options,
            &saved.report.engine_version,
            &saved.report.model,
        ))
        .unwrap(),
    )
    .to_hex()
    .to_string();
    let bytes = serde_json::to_vec(&saved).unwrap();
    let hash = blake3::hash(&bytes).to_hex().to_string();
    let mut archive = c.workspace.clone();
    archive.attachments.insert(
        hash.clone(),
        WorkspaceAttachment {
            kind: ATTACHMENT_KIND.into(),
            name: "Imported historical record".into(),
            media_type: "application/json".into(),
            bytes,
        },
    );
    assert_eq!(
        load(&archive, &hash).unwrap().report.engine_version,
        "archived-untrusted-engine"
    );
    let current = replay(&archive, &hash).unwrap();
    assert_eq!(
        current.report.engine_version,
        collagene_core::verification::ENGINE_VERSION
    );
    assert_eq!(
        current.saved.report.engine_version,
        "archived-untrusted-engine"
    );
    let mut forged = current.clone();
    forged.report.status = ValidationStatus::ConsistentUnderRecordedThresholds;
    assert!(forged.validate_binding(&current.context_revision).is_err());
    let mut position = current;
    position.report.sites[0].position = usize::MAX;
    let token = position.context_revision.clone();
    assert!(position.validate_binding(&token).is_err());
}
#[test]
fn wrong_event_bindings_options_and_attachment_bytes_are_rejected_atomically() {
    let w = fixture();
    let before = serde_json::to_value(&w).unwrap();
    let mut wrong = binding(&w, 0);
    wrong.document_id = binding(&w, 1).document_id;
    assert!(capture(&w, wrong, vec![binding(&w, 2)], &options()).is_err());
    let mut stale = binding(&w, 2);
    stale.revision = "0".repeat(64);
    assert!(capture_revision(&w, &binding(&w, 0), &[stale], &options()).is_err());
    let mut invalid = options();
    invalid.max_total_cells = usize::MAX;
    assert!(capture_revision(&w, &binding(&w, 0), &[binding(&w, 2)], &invalid).is_err());
    assert!(capture_revision(&w, &binding(&w, 0), &vec![binding(&w, 2); 33], &options()).is_err());
    assert_eq!(serde_json::to_value(&w).unwrap(), before);
    let c = capture(&w, binding(&w, 0), vec![binding(&w, 2)], &options()).unwrap();
    let mut broken = c.workspace;
    broken
        .attachments
        .get_mut(&c.attachment_hash)
        .unwrap()
        .bytes
        .push(b' ');
    assert!(load(&broken, &c.attachment_hash).is_err());
}
