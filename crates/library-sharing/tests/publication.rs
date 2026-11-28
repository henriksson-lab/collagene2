use collagene_core::documents::{
    history::transition::Transition, reference_edit::Operation as Edit,
};
use collagene_core::{
    documents::{
        transfer::Candidate,
        workspace::{self, RevisionEvent, Workspace, WorkspaceAttachment, WorkspaceDocument},
    },
    material::{
        lineage::{self as run_lineage, EventBinding, InputEventBinding},
        run as run_material,
    },
    planning::plan::{DesignPlan, Operation, PlanInput, Source},
    sequence::analysis::Primer,
    verification::{self, evidence::Status},
};
use collagene_library_sharing::{self as sharing, Selection};
use collagene_model::{Molecule, Topology};
const DNA: &str = "AGTCCGATGCTACGTAGGACTTCAGACCTGTAGATCGTCA";
fn binding(w: &Workspace, id: &str) -> EventBinding {
    let d = w.documents.iter().find(|d| d.id == id).unwrap();
    EventBinding {
        document_id: d.id.clone(),
        event_id: d.head_event.clone().unwrap(),
        revision: d.head.clone(),
    }
}
fn pcr(template: Source) -> Operation {
    Operation::Pcr {
        template,
        forward: Primer::new("F", "AGTCCGAT", "").unwrap(),
        reverse: Primer::new("R", "TGACGATC", "").unwrap(),
        template_strand: collagene_core::chemistry::Strand::Top,
    }
}
fn add(w: &mut Workspace, id: &str, m: Molecule) {
    let revision = m.revision();
    let e = RevisionEvent::new(
        id,
        revision.clone(),
        vec![],
        "Imported document",
        Transition::IndependentReference,
    )
    .unwrap();
    w.documents.push(WorkspaceDocument {
        id: id.into(),
        head: revision.clone(),
        head_event: Some(e.id.clone()),
        undo: vec![],
        redo: vec![],
        undo_events: vec![],
        redo_events: vec![],
    });
    w.history.push(e);
    w.revisions.insert(revision, m);
}
fn fixture() -> (Workspace, Selection, run_material::Selection, String) {
    let source = Molecule::dna("original stock", DNA, Topology::Linear).unwrap();
    let target = Molecule::dna("independent first PCR reference", DNA, Topology::Linear).unwrap();
    let private = Molecule::dna(
        "UNRELATED PRIVATE PLASMID",
        "CCCCAAAAGGGGTTTT",
        Topology::Circular,
    )
    .unwrap();
    let mut w = Workspace::from_documents(
        "Unpublished lab notebook",
        vec![source.clone(), target.clone(), private.clone()],
        0,
    )
    .unwrap();
    let plan = DesignPlan {
        standard: None,
        name: "first native PCR".into(),
        target_revision: target.revision(),
        inputs: vec![PlanInput {
            id: "stock".into(),
            revision: source.revision(),
            molecule: source,
        }],
        steps: vec![pcr(Source::Input("stock".into()))],
    };
    let first = run_lineage::record_run(
        &w,
        &plan,
        binding(&w, "document-2"),
        vec![InputEventBinding {
            input_id: "stock".into(),
            event: binding(&w, "document-1"),
        }],
        false,
    )
    .unwrap();
    w = first.workspace;
    let product = first.receipt.products[0].clone();
    let origin = run_material::Selection {
        run_id: first.receipt.run_id,
        step_index: 0,
        event: EventBinding {
            document_id: product.document_id,
            event_id: product.event_id,
            revision: product.revision,
        },
    };
    let capsule = collagene_core::material::capsule::capsule(&w, &origin).unwrap();
    let prepared = w.revisions[&origin.event.revision].clone();
    let target2 = Molecule::dna(
        "independent downstream PCR reference",
        DNA,
        Topology::Linear,
    )
    .unwrap();
    add(&mut w, "second-target", target2.clone());
    let plan2 = DesignPlan {
        standard: None,
        name: "saved material then native PCR".into(),
        target_revision: target2.revision(),
        inputs: vec![PlanInput {
            id: "prepared".into(),
            revision: prepared.revision(),
            molecule: prepared,
        }],
        steps: vec![
            Operation::ReconstructRunProduct {
                template: Source::Input("prepared".into()),
                origin: capsule,
            },
            pcr(Source::Step(0)),
        ],
    };
    let second = run_lineage::record_run(
        &w,
        &plan2,
        binding(&w, "second-target"),
        vec![InputEventBinding {
            input_id: "prepared".into(),
            event: origin.event,
        }],
        false,
    )
    .unwrap();
    w = second.workspace;
    let final_product = second.receipt.products[1].clone();
    let material = run_material::Selection {
        run_id: second.receipt.run_id,
        step_index: 1,
        event: EventBinding {
            document_id: final_product.document_id,
            event_id: final_product.event_id,
            revision: final_product.revision,
        },
    };
    // Unrelated proposals and later source edits are local work, not publication
    // dependencies. The source run still binds its exact original import event.
    w.candidates.push(Candidate::new(DesignPlan {
        standard: None,
        name: "private future design".into(),
        target_revision: private.revision(),
        inputs: vec![PlanInput {
            id: "private".into(),
            revision: private.revision(),
            molecule: private,
        }],
        steps: vec![pcr(Source::Input("private".into()))],
    }));
    let old = binding(&w, "document-1");
    let previous_source = w.revisions.get(&old.revision).unwrap();
    let actions = vec![
        Edit::ReplaceSequence {
            start: 0,
            end: previous_source.sequence.len(),
            replacement: "TTTTGGGGCCCCAAAA".into(),
        },
        Edit::Rename {
            name: "PRIVATE UNPUBLISHED SOURCE EDIT".into(),
        },
        Edit::SetTopology {
            topology: Topology::Linear,
        },
    ];
    let edit = Edit::apply_all(previous_source, &actions).unwrap();
    let revision = edit.revision();
    let event = RevisionEvent::new(
        "document-1",
        revision.clone(),
        vec![old.event_id.clone()],
        "Unpublished private source edit",
        Transition::ReferenceEdit {
            source: old.clone(),
            actions,
            supplemental_sources: vec![],
        },
    )
    .unwrap();
    let d = w
        .documents
        .iter_mut()
        .find(|d| d.id == "document-1")
        .unwrap();
    d.head = revision.clone();
    d.head_event = Some(event.id.clone());
    d.undo = vec![old.revision];
    d.undo_events = vec![old.event_id];
    w.revisions.insert(revision.clone(), edit);
    w.history.push(event);
    let ui = WorkspaceAttachment {
        // App's current homology-capable autosave kind (app.rs). Closure treats
        // working bytes opaquely, including unfinished/recovery-only drafts.
        kind: collagene_core::cloning::assembly::library_history::ACTIVE_DRAFT_KIND.into(),
        name: "Private draft".into(),
        media_type: "application/json".into(),
        bytes: br#"{"version":14,"name":"PRIVATE WORKING DRAFT","unpublished_source":"DO NOT PUBLISH"}"#.to_vec(),
    };
    // Attachment identity uses the exact bytes, not a private draft's semantics.
    let ui_hash = hash_bytes(&ui.bytes);
    w.attachments.insert(ui_hash, ui);
    w.validate().unwrap();
    let selected = Selection {
        document_id: material.event.document_id.clone(),
        revision: material.event.revision.clone(),
        head_event: Some(material.event.event_id.clone()),
    };
    (w, selected, material, revision)
}

#[test]
fn published_pcr_history_reopens_and_replays_nested_saved_material_without_private_documents() {
    let (original, selection, material, private_revision) = fixture();
    let fingerprint = run_lineage::workspace_fingerprint(&original).unwrap();
    let export = sharing::document_closure(&original, &selection).unwrap();
    assert!(!export.inventory.documents.contains("document-3"));
    assert!(!export.inventory.revisions.contains(&private_revision));
    assert!(export.inventory.candidates.is_empty());
    assert!(export.workspace.candidates.is_empty());
    assert!(
        !export.workspace.attachments.values().any(
            |a| a.kind == collagene_core::cloning::assembly::library_history::ACTIVE_DRAFT_KIND
        )
    );
    assert_eq!(export.workspace.name, "Plasmid publication");
    assert_eq!(
        export
            .workspace
            .documents
            .iter()
            .find(|d| d.id == "document-1")
            .unwrap()
            .head,
        original
            .documents
            .iter()
            .find(|d| d.id == "document-1")
            .unwrap()
            .undo[0]
    );
    assert_eq!(
        run_lineage::workspace_fingerprint(&original).unwrap(),
        fingerprint
    );
    let loaded = workspace::load_workspace(&export.archive).unwrap();
    let before = run_material::reconstruct(&original, &material).unwrap();
    let after = run_material::reconstruct(&loaded, &material).unwrap();
    assert_eq!(before.reference(), after.reference());
    assert_eq!(after.reference().sequence.as_str(), DNA);
    let fresh = collagene_core::material::capsule::capsule(&loaded, &material).unwrap();
    let actual = after.reference().clone();
    let intended =
        Molecule::dna("third independently declared target", DNA, Topology::Linear).unwrap();
    let downstream = DesignPlan {
        standard: None,
        name: "after library download".into(),
        target_revision: intended.revision(),
        inputs: vec![PlanInput {
            id: "saved".into(),
            revision: actual.revision(),
            molecule: actual,
        }],
        steps: vec![
            Operation::ReconstructRunProduct {
                template: Source::Input("saved".into()),
                origin: fresh,
            },
            pcr(Source::Step(0)),
        ],
    };
    let report = verification::engine::verify(&downstream, &intended);
    assert_eq!(
        report.status,
        Status::ConsistentUnderAssumptions,
        "{:?}",
        report.findings
    );
    assert_eq!(
        report
            .final_product()
            .expect("final declared step must reconstruct")
            .sequence
            .as_str(),
        DNA
    );
    let derived = sharing::inspect_document_closure(&loaded, &selection).unwrap();
    assert_eq!(derived.documents, export.inventory.documents);
    assert_eq!(derived.revisions, export.inventory.revisions);
    assert_eq!(derived.events, export.inventory.events);
    assert_eq!(derived.attachments, export.inventory.attachments);
    assert_eq!(derived.document_heads, export.inventory.document_heads);
}

#[test]
fn plain_plasmid_projection_is_exact_and_does_not_export_other_local_branches() {
    let mut w = Workspace::from_documents(
        "Local workspace",
        vec![
            Molecule::dna("p1", "AGTCCGAT", Topology::Circular).unwrap(),
            Molecule::dna("private", "CCCCAAAA", Topology::Circular).unwrap(),
        ],
        0,
    )
    .unwrap();
    let d = &w.documents[0];
    let selection = Selection {
        document_id: d.id.clone(),
        revision: d.head.clone(),
        head_event: d.head_event.clone(),
    };
    let inventory = sharing::inspect_document_closure(&w, &selection).unwrap();
    assert_eq!(inventory.documents.len(), 1);
    assert_eq!(inventory.events.len(), 1);
    assert_eq!(inventory.revisions.len(), 1);
    assert!(inventory.attachments.is_empty());
    let p = sharing::project_document_closure(&w, &selection).unwrap();
    assert_eq!(p.workspace.documents.len(), 1);
    assert_eq!(
        p.workspace
            .head(&selection.document_id)
            .unwrap()
            .sequence
            .as_str(),
        "AGTCCGAT"
    );
    let mut stale = selection.clone();
    stale.head_event = None;
    assert!(sharing::inspect_document_closure(&w, &stale).is_err());
    // Unknown extension schemas cannot silently hide or publish foreign refs.
    let a = WorkspaceAttachment {
        kind: "unknown-graph-v1".into(),
        name: "unknown".into(),
        media_type: "application/octet-stream".into(),
        bytes: vec![1, 2, 3],
    };
    // Existing workspace serialization performs its own byte-identity checking;
    // construct this fixture using a known checked hash via the protocol crate.
    let hash = hash_bytes(&a.bytes);
    w.attachments.insert(hash, a);
    assert!(matches!(
        sharing::inspect_document_closure(&w, &selection),
        Err(sharing::Error::UnsupportedAttachment(_))
    ));
}
fn hash_bytes(bytes: &[u8]) -> String {
    blake3::hash(bytes).to_hex().to_string()
}

#[test]
fn current_and_historical_app_state_is_local_but_future_or_scientific_schemas_are_not_guessed() {
    let mut frame = Workspace::from_documents(
        "Local library",
        vec![
            Molecule::dna("Shared plasmid", "AGTCCGAT", Topology::Circular).unwrap(),
            Molecule::dna(
                "UNRELATED PRIVATE PLASMID",
                "CCCCAAAAGGGGTTTT",
                Topology::Circular,
            )
            .unwrap(),
        ],
        0,
    )
    .unwrap();
    let selected = &frame.documents[0];
    let selection = Selection {
        document_id: selected.id.clone(),
        revision: selected.head.clone(),
        head_event: selected.head_event.clone(),
    };
    let mut kinds =
        vec![collagene_core::cloning::assembly::library_history::ACTIVE_DRAFT_KIND.to_owned()];
    kinds.extend(
        [
            "assembly-design-v1",
            collagene_core::cloning::assembly::library_history::KIND,
            "primer-pair-draft-v2",
            "annotation-display-v1",
            "lab-library-bindings-v1",
            "lab-library-branches-v1",
        ]
        .map(str::to_owned),
    );
    for kind in &kinds {
        // These are recognized App-side local state families, not native product
        // evidence. No decoding or interpretation of their private working DNA,
        // source pins or future design fields is necessary to omit them safely.
        let bytes=serde_json::to_vec(&serde_json::json!({"kind":kind,"unpublished_source":frame.revisions[&frame.documents[1].head],"private_working_text":"DO NOT PUBLISH","unfinished":true})).unwrap();
        frame.attachments.insert(
            hash_bytes(&bytes),
            WorkspaceAttachment {
                kind: kind.clone(),
                name: format!("Local {kind}"),
                media_type: "application/json".into(),
                bytes,
            },
        );
    }
    frame.validate().unwrap();
    let export = sharing::document_closure(&frame, &selection).unwrap();
    assert_eq!(frame.attachments.len(), kinds.len());
    assert!(export.inventory.attachments.is_empty());
    assert!(export.workspace.attachments.is_empty());
    assert_eq!(export.inventory.omitted.attachments, kinds.len());
    let reopened = workspace::load_workspace(&export.archive).unwrap();
    assert!(reopened.attachments.is_empty());
    assert_eq!(reopened.documents.len(), 1);
    assert!(
        reopened
            .revisions
            .values()
            .all(|m| m.sequence.as_str() != "CCCCAAAAGGGGTTTT")
    );
    assert_eq!(
        reopened
            .head(&selection.document_id)
            .unwrap()
            .sequence
            .as_str(),
        "AGTCCGAT"
    );

    let future_draft = format!(
        "assembly-draft-v{}",
        collagene_core::cloning::assembly::library_history::ACTIVE_DRAFT_VERSION + 1
    );
    for kind in [
        future_draft.as_str(),
        "assembly-draft-v13",
        "primer-pair-draft-v1",
        "primer-pair-draft-v3",
        "assembly-design-v2",
        "pcr-output-reference-v1",
        "mismatch-pcr-output-v1",
        "feature-library-v1",
        "future-scientific-graph-v1",
    ] {
        let mut unsupported = frame.clone();
        let bytes = format!("unsupported {kind}").into_bytes();
        unsupported.attachments.insert(
            hash_bytes(&bytes),
            WorkspaceAttachment {
                kind: kind.into(),
                name: "Not a known private UI contract".into(),
                media_type: "application/octet-stream".into(),
                bytes,
            },
        );
        let result = sharing::inspect_document_closure(&unsupported, &selection);
        if kind.starts_with("assembly-draft-") {
            assert!(
                matches!(result, Err(sharing::Error::Native(collagene_model::Error::Format(message)))
                if message.contains("unsupported internal assembly/library attachment version")),
                "{kind} must fail native admission before closure traversal"
            );
        } else {
            assert!(
                matches!(result,Err(sharing::Error::UnsupportedAttachment(actual)) if actual==kind),
                "{kind} was guessed to be presentation state"
            );
        }
    }
}
