use collagene_core::{
    Alphabet, Molecule, Sequence, Topology, rna_chemistry as rna, rna_history,
    run_lineage::EventBinding, workspace::Workspace,
};
use collagene_library_sharing::{self as sharing, Selection};
#[test]
fn only_declarations_on_retained_exact_events_follow_selected_document_sharing() {
    let m = Molecule {
        name: "RNA".into(),
        sequence: Sequence::new(Alphabet::Rna, "AC").unwrap(),
        topology: Topology::Linear,
        features: vec![],
        metadata: Default::default(),
        original_genbank: None,
    };
    let mut w = Workspace::from_documents("RNA publication", vec![m.clone(), m], 0).unwrap();
    let mut hashes = vec![];
    for i in 0..2 {
        let d = &w.documents[i];
        let record = rna_history::Record {
            version: rna_history::VERSION,
            name: format!("Model {}", i + 1),
            source: EventBinding {
                document_id: d.id.clone(),
                event_id: d.head_event.clone().unwrap(),
                revision: d.head.clone(),
            },
            declaration: rna::Declaration {
                version: rna::DECLARATION_VERSION,
                source_revision: d.head.clone(),
                coverage: rna::Coverage::Partial,
                ends: rna::Ends::Linear {
                    five: rna::FivePrime::Hydroxyl,
                    three: rna::ThreePrime::Hydroxyl,
                },
                sites: vec![],
                links: vec![],
                annotation_reviews: vec![],
                provenance: rna::Provenance {
                    statement: "Independent hypothesis for one document event".into(),
                    reference: None,
                },
            },
        };
        let (next, hash) = rna_history::capture(&w, &record).unwrap();
        w = next;
        hashes.push(hash);
    }
    let d = &w.documents[0];
    let selection = Selection {
        document_id: d.id.clone(),
        revision: d.head.clone(),
        head_event: d.head_event.clone(),
    };
    let export = sharing::document_closure(&w, &selection).unwrap();
    assert!(export.workspace.attachments.contains_key(&hashes[0]));
    assert!(!export.workspace.attachments.contains_key(&hashes[1]));
    assert_eq!(export.workspace.documents.len(), 1);
    assert_eq!(rna_history::list(&export.workspace).unwrap().len(), 1);
    let loaded = collagene_core::workspace::load_workspace(&export.archive).unwrap();
    assert_eq!(
        rna_history::load(&loaded, &hashes[0])
            .unwrap()
            .source
            .document_id,
        selection.document_id
    );
    assert_eq!(w.attachments.len(), 2);
}

#[test]
fn ancestor_declarations_follow_a_new_head_without_becoming_a_model_for_that_head() {
    let m = Molecule {
        name: "Original RNA".into(),
        sequence: Sequence::new(Alphabet::Rna, "AC").unwrap(),
        topology: Topology::Linear,
        features: vec![],
        metadata: Default::default(),
        original_genbank: None,
    };
    let w = Workspace::from_documents("RNA ancestry", vec![m], 0).unwrap();
    let d = &w.documents[0];
    let r = rna_history::Record {
        version: rna_history::VERSION,
        name: "Original model".into(),
        source: EventBinding {
            document_id: d.id.clone(),
            event_id: d.head_event.clone().unwrap(),
            revision: d.head.clone(),
        },
        declaration: rna::Declaration {
            version: rna::DECLARATION_VERSION,
            source_revision: d.head.clone(),
            coverage: rna::Coverage::Partial,
            ends: rna::Ends::Linear {
                five: rna::FivePrime::Hydroxyl,
                three: rna::ThreePrime::Hydroxyl,
            },
            sites: vec![],
            links: vec![],
            annotation_reviews: vec![],
            provenance: rna::Provenance {
                statement: "Model for the original version only".into(),
                reference: None,
            },
        },
    };
    let (mut w, hash) = rna_history::capture(&w, &r).unwrap();
    let mut histories = collagene_core::history::DocumentHistory::restore_all(&w).unwrap();
    let mut later = histories[0].molecule.clone();
    later.sequence = Sequence::new(Alphabet::Rna, "AG").unwrap();
    histories[0].commit(later, "Edit RNA sequence").unwrap();
    w.documents.clear();
    w.history.clear();
    for d in histories {
        d.contribute(&mut w).unwrap();
    }
    w.validate().unwrap();
    let d = &w.documents[0];
    let selected = Selection {
        document_id: d.id.clone(),
        revision: d.head.clone(),
        head_event: d.head_event.clone(),
    };
    let exported = sharing::document_closure(&w, &selected).unwrap();
    let saved = rna_history::load(&exported.workspace, &hash).unwrap();
    assert_eq!(saved.source, r.source);
    assert_ne!(saved.source.revision, selected.revision);
    assert!(
        exported
            .workspace
            .history
            .iter()
            .any(|e| e.id == r.source.event_id)
    );
    assert!(
        rna::analyze(
            exported.workspace.head(&selected.document_id).unwrap(),
            &saved.declaration
        )
        .is_err()
    );
}
