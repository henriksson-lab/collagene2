use collagene_core::documents::reference_edit::Operation as Edit;
use collagene_core::{
    documents::{history::DocumentHistory, workspace::Workspace},
    material::{
        lineage::{self as run_lineage, EventBinding, InputEventBinding},
        outputs::*,
    },
    planning::plan::{DesignPlan, Operation, PlanInput, Source},
    sequence::analysis::Primer,
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
fn fixture() -> (Workspace, DesignPlan) {
    let s = Molecule::dna("original", "ACGACCCCGGGGTAAC", Topology::Linear).unwrap();
    let t = Molecule::dna("independent target", "ACGACCCCGGGGTAAC", Topology::Linear).unwrap();
    let p = DesignPlan {
        standard: None,
        name: "literal PCR".into(),
        target_revision: t.revision(),
        inputs: vec![PlanInput {
            id: "source".into(),
            revision: s.revision(),
            molecule: s.clone(),
        }],
        steps: vec![Operation::Pcr {
            template: Source::Input("source".into()),
            forward: Primer::new("forward", "ACGACCCC", "").unwrap(),
            reverse: Primer::new("reverse", "GTTACCCC", "").unwrap(),
            template_strand: collagene_core::chemistry::Strand::Top,
        }],
    };
    (
        Workspace::from_documents("outputs", vec![s, t], 1).unwrap(),
        p,
    )
}
fn capture(w: &Workspace, p: &DesignPlan, o: &Options) -> Commit {
    finalize(
        w,
        p,
        binding(w, 1),
        vec![InputEventBinding {
            input_id: "source".into(),
            event: binding(w, 0),
        }],
        o,
    )
    .unwrap()
}
#[test]
fn named_output_retains_original_run_events_and_is_idempotent() {
    let (w, p) = fixture();
    let before = run_lineage::workspace_fingerprint(&w).unwrap();
    let o = Options {
        name: "My construct".into(),
        save_primers: true,
        ..Options::default()
    };
    let c = capture(&w, &p, &o);
    assert_eq!(run_lineage::workspace_fingerprint(&w).unwrap(), before);
    assert_eq!(
        c.workspace
            .head(&c.receipt.construct.document_id)
            .unwrap()
            .name,
        "My construct"
    );
    assert_eq!(c.workspace.documents.len(), 3);
    assert_ne!(
        c.receipt.construct.event_id,
        c.receipt.run.products[0].event_id
    );
    let replay = run_lineage::replay_run(&c.workspace, &c.receipt.run.run_id).unwrap();
    assert!(replay.products_match);
    let a = c
        .workspace
        .attachments
        .get(c.receipt.primer_attachment.as_ref().unwrap())
        .unwrap();
    let table: collagene_core::formats::oligos::OligoTable =
        serde_json::from_slice(&a.bytes).unwrap();
    assert_eq!(table.records.len(), 2);
    assert_eq!(
        table.records[0].raw_ordered_sequence.as_deref(),
        Some("ACGACCCC")
    );
    assert_eq!(
        table.records[1].raw_ordered_sequence.as_deref(),
        Some("GTTACCCC")
    );
    let again = capture(&c.workspace, &p, &o);
    assert!(again.receipt.run.already_recorded);
    assert_eq!(again.workspace.history.len(), c.workspace.history.len());
    assert_eq!(
        again.workspace.attachments.len(),
        c.workspace.attachments.len()
    );
    assert_eq!(again.receipt.construct, c.receipt.construct);
    assert!(c.validate_binding(&"0".repeat(64)).is_err());
}
#[test]
fn prepared_documents_are_hidden_only_by_presentation_attachment() {
    let (w, mut p) = fixture();
    p.steps.push(Operation::TreatEnds {
        template: Source::Step(0),
        treatment: collagene_core::chemistry::EndTreatment::Phosphorylate,
    });
    let c = capture(&w, &p, &Options::default());
    assert_eq!(c.receipt.hidden_products.len(), 1);
    assert_eq!(c.workspace.documents.len(), 4);
    assert_eq!(load_hidden(&c.workspace).unwrap().len(), 1);
    assert!(
        run_lineage::replay_run(&c.workspace, &c.receipt.run.run_id)
            .unwrap()
            .products_match
    );
    let show = capture(
        &c.workspace,
        &p,
        &Options {
            show_prepared: true,
            ..Options::default()
        },
    );
    assert!(load_hidden(&show.workspace).unwrap().is_empty());
    assert_eq!(show.workspace.history.len(), c.workspace.history.len());
}
#[test]
fn edited_output_is_preserved_and_original_product_fork_is_reused() {
    let (w, p) = fixture();
    let c = capture(&w, &p, &Options::default());
    let old = c.receipt.construct.document_id.clone();
    let mut changed = c.workspace.clone();
    let mut histories = DocumentHistory::restore_all(&changed).unwrap();
    let h = histories.iter_mut().find(|h| h.id == old).unwrap();
    let edit_end = h.molecule.sequence.len();
    h.commit(
        vec![
            Edit::ReplaceSequence {
                start: 0,
                end: edit_end,
                replacement: "TTTT".into(),
            },
            Edit::Rename {
                name: "user changed construct".into(),
            },
            Edit::SetTopology {
                topology: Topology::Linear,
            },
        ],
        "User changed DNA",
    )
    .unwrap();
    let edited = h.molecule.clone();
    assert_eq!(edited.sequence.as_str(), "TTTT");
    changed.documents.clear();
    changed.history.clear();
    for h in &histories {
        h.contribute(&mut changed).unwrap();
    }
    let newer = capture(&changed, &p, &Options::default());
    assert_ne!(newer.receipt.construct.document_id, old);
    assert_eq!(newer.workspace.head(&old).unwrap(), &edited);
    assert!(newer.receipt.visible_products.contains(&old));
    assert!(
        run_lineage::replay_run(&newer.workspace, &newer.receipt.run.run_id)
            .unwrap()
            .products_match
    );
    let again = capture(&newer.workspace, &p, &Options::default());
    assert_eq!(again.receipt.construct, newer.receipt.construct);
    assert_eq!(
        again.workspace.documents.len(),
        newer.workspace.documents.len()
    );
    assert_eq!(again.workspace.history.len(), newer.workspace.history.len());
}
#[test]
fn stale_bindings_wrong_target_and_empty_oligo_collection_are_atomic_errors() {
    let (w, mut p) = fixture();
    let before = run_lineage::workspace_fingerprint(&w).unwrap();
    let mut stale = binding(&w, 0);
    stale.revision = "0".repeat(64);
    assert!(
        finalize(
            &w,
            &p,
            binding(&w, 1),
            vec![InputEventBinding {
                input_id: "source".into(),
                event: stale
            }],
            &Options::default()
        )
        .is_err()
    );
    p.target_revision = "0".repeat(64);
    assert!(finalize_revision(&w, &p, binding(&w, 1), vec![], &Options::default()).is_err());
    let (w, mut no_primers) = fixture();
    no_primers.steps = vec![Operation::DeclareSourceTermini {
        template: Source::Input("source".into()),
        left: collagene_core::chemistry::Terminus::CLEAVED,
        right: collagene_core::chemistry::Terminus::CLEAVED,
    }];
    assert!(
        finalize(
            &w,
            &no_primers,
            binding(&w, 1),
            vec![InputEventBinding {
                input_id: "source".into(),
                event: binding(&w, 0)
            }],
            &Options {
                save_primers: true,
                ..Options::default()
            }
        )
        .is_err()
    );
    let (w, p) = fixture();
    let first = capture(&w, &p, &Options::default());
    let mut malformed = first.clone();
    malformed.receipt.run.products[0].revision = "0".repeat(64);
    assert!(
        malformed
            .validate_binding(&first.receipt.context_revision)
            .is_err()
    );
    assert_eq!(run_lineage::workspace_fingerprint(&w).unwrap(), before);
}

#[test]
fn historical_source_binding_survives_source_head_edits_and_context_changes() {
    let (w, p) = fixture();
    let source = binding(&w, 0);
    let target = binding(&w, 1);
    let inputs = vec![InputEventBinding {
        input_id: "source".into(),
        event: source.clone(),
    }];
    let old_context =
        finalize_revision(&w, &p, target.clone(), inputs.clone(), &Options::default()).unwrap();
    let mut changed = w.clone();
    let mut histories = DocumentHistory::restore_all(&changed).unwrap();
    let edit_end = histories[0].molecule.sequence.len();
    histories[0]
        .commit(
            vec![
                Edit::ReplaceSequence {
                    start: 0,
                    end: edit_end,
                    replacement: "CCCC".into(),
                },
                Edit::Rename {
                    name: "advanced source".into(),
                },
                Edit::SetTopology {
                    topology: Topology::Linear,
                },
            ],
            "User edit after selecting source",
        )
        .unwrap();
    let edited = histories[0].molecule.clone();
    assert_eq!(edited.sequence.as_str(), "CCCC");
    changed.documents.clear();
    changed.history.clear();
    for h in &histories {
        h.contribute(&mut changed).unwrap();
    }
    let context = finalize_revision(
        &changed,
        &p,
        target.clone(),
        inputs.clone(),
        &Options::default(),
    )
    .unwrap();
    assert_ne!(old_context, context);
    let c = finalize(&changed, &p, target, inputs, &Options::default()).unwrap();
    assert_eq!(c.workspace.head(&source.document_id).unwrap(), &edited);
    let run = run_lineage::load_run(&c.workspace, &c.receipt.run.run_id).unwrap();
    assert_eq!(run.inputs[0].event, source);
    assert!(
        run_lineage::replay_run(&c.workspace, &run.run_id)
            .unwrap()
            .products_match
    );
    let mut broken = c.workspace.clone();
    let old = broken
        .attachments
        .iter()
        .find(|(_, a)| a.kind == VISIBILITY_KIND)
        .unwrap()
        .0
        .clone();
    broken.attachments.remove(&old);
    let bytes = br#"{"version":1,"hidden":["missing-document"]}"#.to_vec();
    broken.attachments.insert(
        blake3::hash(&bytes).to_hex().to_string(),
        collagene_core::documents::workspace::WorkspaceAttachment {
            kind: VISIBILITY_KIND.into(),
            name: "bad references".into(),
            media_type: "application/json".into(),
            bytes,
        },
    );
    assert!(load_hidden(&broken).is_err());
}

#[test]
fn ordinary_capture_does_not_apply_the_optional_oligo_order_length_profile() {
    let dna = format!("{}{}", "A".repeat(201), "C".repeat(201));
    let source = Molecule::dna("long domains", &dna, Topology::Linear).unwrap();
    let target = Molecule::dna("independent 402nt target", &dna, Topology::Linear).unwrap();
    let plan = DesignPlan {
        standard: None,
        name: "long exact primers".into(),
        target_revision: target.revision(),
        inputs: vec![PlanInput {
            id: "source".into(),
            revision: source.revision(),
            molecule: source.clone(),
        }],
        steps: vec![Operation::Pcr {
            template: Source::Input("source".into()),
            forward: Primer::new("forward", &"A".repeat(201), "").unwrap(),
            reverse: Primer::new("reverse", &"G".repeat(201), "").unwrap(),
            template_strand: collagene_core::chemistry::Strand::Top,
        }],
    };
    let w = Workspace::from_documents("long primer capture", vec![source, target], 1).unwrap();
    let before = run_lineage::workspace_fingerprint(&w).unwrap();
    let c = capture(&w, &plan, &Options::default());
    assert_eq!(
        c.workspace
            .head(&c.receipt.construct.document_id)
            .unwrap()
            .sequence
            .len(),
        402
    );
    let result = finalize(
        &w,
        &plan,
        binding(&w, 1),
        vec![InputEventBinding {
            input_id: "source".into(),
            event: binding(&w, 0),
        }],
        &Options {
            save_primers: true,
            ..Options::default()
        },
    );
    assert!(result.is_err());
    assert_eq!(run_lineage::workspace_fingerprint(&w).unwrap(), before);
}
