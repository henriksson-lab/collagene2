use collagene_core::{Molecule, Topology, history::DocumentHistory, workspace::Workspace};
#[test]
fn transformed_document_retains_exact_source_event_and_survives_archive_roundtrip() {
    let original = Molecule::dna("original", "ACGT", Topology::Linear).unwrap();
    let mut source = DocumentHistory::new("source", original.clone()).unwrap();
    let equal = DocumentHistory::new("equal-content-other-document", original).unwrap();
    let event = source.head_event().to_string();
    source
        .commit(
            Molecule::dna("edited", "AGGT", Topology::Linear).unwrap(),
            "substitution",
        )
        .unwrap();
    let navigation = source.navigation_event_ids();
    let head = source.head_event().to_string();
    let product = Molecule::dna("derived", "AAATTT", Topology::Linear).unwrap();
    let derived = source
        .derive_event(
            &event,
            "derived",
            product.clone(),
            "Source-bound design; no physical authority",
        )
        .unwrap();
    assert_eq!(derived.summaries()[0].parents, vec![event.clone()]);
    assert!(!derived.can_undo());
    assert!(!derived.can_redo());
    assert_eq!(source.head_event(), head);
    assert_eq!(source.navigation_event_ids(), navigation);
    let other = equal
        .derive_event(
            equal.head_event(),
            "derived-other",
            product,
            "same sequence, distinct source",
        )
        .unwrap();
    assert_ne!(derived.summaries()[0].parents, other.summaries()[0].parents);
    let mut frame = Workspace::from_documents("test", vec![source.molecule.clone()], 0).unwrap();
    frame.documents.clear();
    frame.history.clear();
    frame.revisions.clear();
    source.contribute(&mut frame).unwrap();
    equal.contribute(&mut frame).unwrap();
    derived.contribute(&mut frame).unwrap();
    other.contribute(&mut frame).unwrap();
    frame.active_document = derived.id.clone();
    frame.validate().unwrap();
    let bytes = collagene_core::workspace::save_workspace(&frame).unwrap();
    let loaded = collagene_core::workspace::load_workspace(&bytes).unwrap();
    let restored = DocumentHistory::restore_all(&loaded).unwrap();
    let result = restored.iter().find(|d| d.id == "derived").unwrap();
    assert_eq!(result.summaries()[0].parents, vec![event]);
    assert_eq!(
        restored
            .iter()
            .find(|d| d.id == "source")
            .unwrap()
            .navigation_event_ids(),
        navigation
    );
}
#[test]
fn missing_source_and_reused_identity_are_rejected_without_mutation() {
    let source = DocumentHistory::new(
        "source",
        Molecule::dna("source", "ACGT", Topology::Linear).unwrap(),
    )
    .unwrap();
    let head = source.head_event().to_string();
    let product = Molecule::dna("derived", "GATTACA", Topology::Linear).unwrap();
    assert!(
        source
            .derive_event("event-missing", "derived", product.clone(), "design")
            .is_err()
    );
    assert!(
        source
            .derive_event(&head, "source", product.clone(), "design")
            .is_err()
    );
    assert!(
        source
            .derive_event(&head, "derived", product, "x".repeat(4097))
            .is_err()
    );
    assert_eq!(source.head_event(), head);
    assert_eq!(source.event_count(), 1);
}
