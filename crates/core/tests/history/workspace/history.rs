use collagene_core::documents::{
    history::DocumentHistory,
    workspace::{self, RevisionEvent, Workspace},
};
use collagene_model::{Feature, Location, Molecule, Topology};
use std::collections::BTreeMap;
fn m(name: &str) -> Molecule {
    Molecule::dna(name, "ACGTACGT", Topology::Linear).unwrap()
}
fn archive(docs: &[DocumentHistory]) -> Workspace {
    let mut w = Workspace {
        name: "History fixture".into(),
        active_document: docs[0].id.clone(),
        documents: vec![],
        revisions: BTreeMap::new(),
        history: vec![],
        candidates: vec![],
        active_candidate: None,
        attachments: BTreeMap::new(),
    };
    for d in docs {
        d.contribute(&mut w).unwrap();
    }
    w.validate().unwrap();
    w
}
#[test]
fn duplicate_content_commits_retain_exact_navigation_events_after_native_reopen() {
    let mut d = DocumentHistory::new("doc", m("A")).unwrap();
    let root = d.head_event().to_owned();
    d.commit(m("B"), "rename").unwrap();
    let middle = d.head_event().to_owned();
    d.commit(m("A"), "return to A").unwrap();
    let returned = d.head_event().to_owned();
    assert_ne!(root, returned);
    assert_eq!(d.event_count(), 3);
    let bytes = workspace::save_workspace(&archive(&[d])).unwrap();
    let w = workspace::load_workspace(&bytes).unwrap();
    assert_eq!(w.revisions.len(), 2);
    let mut d = DocumentHistory::restore_all(&w).unwrap().remove(0);
    assert_eq!(d.head_event(), returned);
    assert!(d.undo().unwrap());
    assert_eq!(d.head_event(), middle);
    assert!(d.undo().unwrap());
    assert_eq!(d.head_event(), root);
    assert!(!d.undo().unwrap());
    assert!(d.redo().unwrap());
    assert_eq!(d.head_event(), middle);
    assert!(d.redo().unwrap());
    assert_eq!(d.head_event(), returned);
}
#[test]
fn editing_after_undo_retains_the_old_branch_and_checkout_is_reversible() {
    let mut d = DocumentHistory::new("doc", m("A")).unwrap();
    d.commit(m("B"), "B").unwrap();
    let b = d.head_event().to_owned();
    d.commit(m("C"), "C").unwrap();
    let c = d.head_event().to_owned();
    d.undo().unwrap();
    d.commit(m("D"), "D").unwrap();
    let current = d.head_event().to_owned();
    assert!(!d.can_redo());
    assert_eq!(d.event_count(), 4);
    let w = archive(&[d]);
    assert_eq!(
        w.history
            .iter()
            .filter(|e| e.parents == vec![b.clone()])
            .count(),
        2
    );
    let mut restored = DocumentHistory::restore_all(&w).unwrap().remove(0);
    restored.checkout(&c).unwrap();
    assert_eq!(restored.molecule.name, "C");
    restored.undo().unwrap();
    assert_eq!(restored.head_event(), current);
    assert_eq!(restored.molecule.name, "D");
    assert!(restored.checkout("unavailable").is_err());
    assert_eq!(restored.molecule.name, "D");
}
#[test]
fn navigation_window_does_not_destroy_older_revision_evidence() {
    let mut d = DocumentHistory::new("doc", m("initial")).unwrap();
    let first = d.head_event().to_owned();
    for i in 0..105 {
        d.commit(m(&format!("state_{i}")), "rename").unwrap();
    }
    assert_eq!(d.event_count(), 106);
    let w = archive(&[d.clone()]);
    assert_eq!(w.documents[0].undo.len(), 100);
    assert_eq!(w.history.len(), 106);
    d.checkout(&first).unwrap();
    assert_eq!(d.molecule.name, "initial");
}
#[test]
fn invalid_commits_are_atomic_and_external_unrecorded_changes_are_detected() {
    let mut d = DocumentHistory::new("doc", m("initial")).unwrap();
    let head = d.head_event().to_owned();
    let mut invalid = m("invalid");
    invalid.features.push(Feature {
        kind: "misc_feature".into(),
        location: Location::range(0, 99),
        qualifiers: vec![],
        needs_review: false,
    });
    assert!(d.commit(invalid, "bad").is_err());
    assert_eq!(d.head_event(), head);
    assert_eq!(d.event_count(), 1);
    assert!(!d.can_undo());
    d.molecule.name = "unrecorded".into();
    assert!(d.commit(m("next"), "edit").is_err());
}
#[test]
fn native_navigation_rejects_missing_event_ids_instead_of_inventing_recovery_roots() {
    let mut document = DocumentHistory::new("exact", m("root")).unwrap();
    document.commit(m("middle"), "first edit").unwrap();
    document.commit(m("last"), "second edit").unwrap();
    document.undo().unwrap();
    let current = archive(&[document]);
    assert_eq!(current.documents[0].undo_events.len(), 1);
    assert_eq!(current.documents[0].redo_events.len(), 1);
    for missing in ["head", "undo", "redo"] {
        let mut incomplete = current.clone();
        match missing {
            "head" => incomplete.documents[0].head_event = None,
            "undo" => incomplete.documents[0].undo_events.clear(),
            "redo" => incomplete.documents[0].redo_events.clear(),
            _ => unreachable!(),
        };
        assert!(
            incomplete.validate().is_err(),
            "Missing {missing} event identity was accepted"
        );
        assert!(DocumentHistory::restore_all(&incomplete).is_err());
        assert!(workspace::save_workspace(&incomplete).is_err());
    }
    let mut no_history = current.clone();
    no_history.history.clear();
    assert!(DocumentHistory::restore_all(&no_history).is_err());
    let fresh = archive(&[DocumentHistory::new("imported", m("External import")).unwrap()]);
    let restored = DocumentHistory::restore_all(&fresh).unwrap();
    assert_eq!(
        restored[0].head_event(),
        fresh.documents[0].head_event.as_deref().unwrap()
    );
    assert_eq!(restored[0].event_count(), 1);
    assert!(!restored[0].can_undo());
    assert!(!restored[0].can_redo());
    assert_eq!(
        serde_json::to_value(&archive(&restored).history).unwrap(),
        serde_json::to_value(&fresh.history).unwrap()
    );
}

#[test]
fn shared_revisions_and_cross_document_parents_survive_restore_and_new_edits() {
    let mut w = Workspace::from_documents("shared", vec![m("same"), m("same")], 0).unwrap();
    let origin = w.documents[0].head_event.clone().unwrap();
    let derived = RevisionEvent::new(
        w.documents[1].id.clone(),
        w.documents[1].head.clone(),
        vec![origin.clone()],
        "Derived reference snapshot",
    )
    .unwrap();
    w.documents[1].head_event = Some(derived.id.clone());
    w.history.push(derived);
    w.validate().unwrap();
    let mut docs = DocumentHistory::restore_all(&w).unwrap();
    assert_eq!(archive(&docs).revisions.len(), 1);
    docs[0].commit(m("changed"), "edit").unwrap();
    let restored = archive(&docs);
    assert!(
        restored
            .history
            .iter()
            .any(|e| e.parents == vec![origin.clone()])
    );
    assert_eq!(restored.documents.len(), 2);
}

#[test]
fn exact_navigation_and_cross_document_fork_do_not_conflate_repeated_content() {
    let mut source = DocumentHistory::new("source", m("A")).unwrap();
    let root = source.head_event().to_owned();
    source.commit(m("B"), "rename").unwrap();
    let middle = source.head_event().to_owned();
    source.commit(m("A"), "return").unwrap();
    let returned = source.head_event().to_owned();
    assert_ne!(root, returned);
    assert_eq!(
        source.navigation_event_ids(),
        (vec![root.clone(), middle.clone()], vec![])
    );
    source.undo().unwrap();
    assert_eq!(
        source.navigation_event_ids(),
        (vec![root.clone()], vec![returned.clone()])
    );
    source.checkout(&root).unwrap();
    assert_eq!(
        source.navigation_event_ids(),
        (vec![root.clone(), middle], vec![])
    );
    let fork = source.fork_event(&returned, "fork").unwrap();
    assert_eq!(
        fork.molecule.revision(),
        source.revision_molecule(&returned).unwrap().revision()
    );
    assert_eq!(fork.navigation_event_ids(), (vec![], vec![]));
    assert_eq!(fork.summaries()[0].parents, vec![returned.clone()]);
    assert_ne!(fork.head_event(), returned);
    assert!(source.fork_event(&returned, "source").is_err());
    assert!(source.fork_event("missing", "other").is_err());
    let archive = archive(&[source, fork]);
    let bytes = workspace::save_workspace(&archive).unwrap();
    let reopened =
        DocumentHistory::restore_all(&workspace::load_workspace(&bytes).unwrap()).unwrap();
    assert_eq!(reopened[1].summaries()[0].parents, vec![returned]);
    assert_eq!(
        reopened[0].navigation_event_ids(),
        (
            vec![root, archive.documents[0].undo_events[1].clone()],
            vec![]
        )
    );
}
