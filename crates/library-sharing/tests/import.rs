use collagene_core::documents::{
    history::transition::Transition, reference_edit::Operation as Edit,
};
use collagene_core::{
    documents::workspace::{RevisionEvent, Workspace, WorkspaceDocument},
    material::{
        lineage::{self as run_lineage, EventBinding, InputEventBinding},
        run as run_material,
    },
    planning::plan::{DesignPlan, Operation, PlanInput, Source},
    sequence::analysis::Primer,
    verification::evidence::Status,
};
use collagene_library_sharing::{self as sharing, Selection};
use collagene_model::{Molecule, Topology};
// Module registration is parent-owned until the current coordinated native batch
// ends. This test target then imports the normal public module.
use collagene_library_sharing::import::{
    self as importer, IdentityKind, ImportError, ImportOptions, SelectedPolicy,
};
const DNA: &str = "AGTCCGATGCTACGTAGGACTTCAGACCTGTAGATCGTCA";
fn pin(w: &Workspace, id: &str) -> Selection {
    let d = w.documents.iter().find(|d| d.id == id).unwrap();
    Selection {
        document_id: d.id.clone(),
        revision: d.head.clone(),
        head_event: d.head_event.clone(),
    }
}
fn event(w: &Workspace, index: usize) -> EventBinding {
    let d = &w.documents[index];
    EventBinding {
        document_id: d.id.clone(),
        revision: d.head.clone(),
        event_id: d.head_event.clone().unwrap(),
    }
}
fn add(w: &mut Workspace, id: &str, name: &str, text: &str) {
    let m = Molecule::dna(name, text, Topology::Linear).unwrap();
    let revision = m.revision();
    let e = RevisionEvent::new(
        id,
        revision.clone(),
        vec![],
        "Local import",
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
    w.revisions.insert(revision, m);
    w.history.push(e);
}
fn advance(w: &mut Workspace, id: &str, text: &str) -> Selection {
    let before = pin(w, id);
    let previous = w.head(id).unwrap();
    let actions = vec![
        Edit::ReplaceSequence {
            start: 0,
            end: previous.sequence.len(),
            replacement: text.into(),
        },
        Edit::Rename {
            name: "Local branch".into(),
        },
        Edit::SetTopology {
            topology: Topology::Linear,
        },
    ];
    let m = Edit::apply_all(previous, &actions).unwrap();
    let revision = m.revision();
    let e = RevisionEvent::new(
        id,
        revision.clone(),
        before.head_event.iter().cloned().collect(),
        "Local edit",
        Transition::ReferenceEdit {
            source: collagene_core::material::lineage::EventBinding {
                document_id: id.into(),
                event_id: before.head_event.clone().unwrap(),
                revision: before.revision.clone(),
            },
            actions,
            supplemental_sources: vec![],
        },
    )
    .unwrap();
    let d = w.documents.iter_mut().find(|d| d.id == id).unwrap();
    d.undo.push(d.head.clone());
    if let Some(event) = d.head_event.clone() {
        d.undo_events.push(event);
    }
    d.head = revision.clone();
    d.head_event = Some(e.id.clone());
    d.redo.clear();
    d.redo_events.clear();
    w.revisions.insert(revision, m);
    w.history.push(e);
    w.validate().unwrap();
    before
}
fn pcr(template: Source) -> Operation {
    Operation::Pcr {
        template,
        forward: Primer::new("F", "AGTCCGAT", "").unwrap(),
        reverse: Primer::new("R", "TGACGATC", "").unwrap(),
        template_strand: collagene_core::chemistry::Strand::Top,
    }
}
fn run_fixture() -> (Workspace, Workspace, Selection, run_material::Selection) {
    let source = Molecule::dna("Original stock", DNA, Topology::Linear).unwrap();
    let target = Molecule::dna("Independent product", DNA, Topology::Linear).unwrap();
    let base = Workspace::from_documents("Local library", vec![source.clone(), target.clone()], 0)
        .unwrap();
    let plan = DesignPlan {
        name: "PCR of original stock".into(),
        target_revision: target.revision(),
        inputs: vec![PlanInput {
            id: "stock".into(),
            revision: source.revision(),
            molecule: source,
        }],
        steps: vec![pcr(Source::Input("stock".into()))],
    };
    let first = run_lineage::record_run(
        &base,
        &plan,
        event(&base, 1),
        vec![InputEventBinding {
            input_id: "stock".into(),
            event: event(&base, 0),
        }],
        false,
    )
    .unwrap();
    let product = &first.receipt.products[0];
    let origin = run_material::Selection {
        run_id: first.receipt.run_id.clone(),
        step_index: 0,
        event: EventBinding {
            document_id: product.document_id.clone(),
            event_id: product.event_id.clone(),
            revision: product.revision.clone(),
        },
    };
    let mut remote = first.workspace;
    let stock = remote.revisions[&origin.event.revision].clone();
    let target = Molecule::dna("Independent descendant", DNA, Topology::Linear).unwrap();
    add(&mut remote, "downstream-target", &target.name, DNA);
    let required = pin(&remote, "downstream-target");
    let capsule = collagene_core::material::capsule::capsule(&remote, &origin).unwrap();
    let plan = DesignPlan {
        name: "Fresh use of saved PCR stock".into(),
        target_revision: required.revision.clone(),
        inputs: vec![PlanInput {
            id: "saved".into(),
            revision: stock.revision(),
            molecule: stock,
        }],
        steps: vec![
            Operation::ReconstructRunProduct {
                template: Source::Input("saved".into()),
                origin: capsule,
            },
            pcr(Source::Step(0)),
        ],
    };
    let next = run_lineage::record_run(
        &remote,
        &plan,
        EventBinding {
            document_id: required.document_id,
            event_id: required.head_event.unwrap(),
            revision: required.revision,
        },
        vec![InputEventBinding {
            input_id: "saved".into(),
            event: origin.event,
        }],
        false,
    )
    .unwrap();
    let product = &next.receipt.products[1];
    let material = run_material::Selection {
        run_id: next.receipt.run_id.clone(),
        step_index: 1,
        event: EventBinding {
            document_id: product.document_id.clone(),
            event_id: product.event_id.clone(),
            revision: product.revision.clone(),
        },
    };
    let selected = Selection {
        document_id: material.event.document_id.clone(),
        revision: material.event.revision.clone(),
        head_event: Some(material.event.event_id.clone()),
    };
    let fetched = sharing::project_document_closure(&next.workspace, &selected)
        .unwrap()
        .workspace;
    (base, fetched, selected, material)
}
fn json(w: &Workspace) -> Vec<u8> {
    serde_json::to_vec(w).unwrap()
}
#[test]
fn imported_fromrun_closure_keeps_private_local_heads_and_exact_material_ancestry() {
    let (mut local, fetched, selected, material) = run_fixture();
    let original_source = pin(&local, "document-1");
    advance(&mut local, "document-1", "TTTTGGGGCCCCAAAA");
    add(
        &mut local,
        "unrelated-private",
        "Unpublished private plasmid",
        "CCGGAATTCCGG",
    );
    local.active_document = "unrelated-private".into();
    let before = json(&local);
    let fetched_before = json(&fetched);
    let private_head = local
        .documents
        .iter()
        .find(|d| d.id == "document-1")
        .unwrap()
        .clone();
    let unrelated = local
        .documents
        .iter()
        .find(|d| d.id == "unrelated-private")
        .unwrap()
        .clone();
    let report =
        importer::merge_document_closure(&local, &fetched, &selected, &Default::default()).unwrap();
    assert_eq!(json(&local), before);
    assert_eq!(json(&fetched), fetched_before);
    assert_eq!(
        serde_json::to_value(
            report
                .workspace
                .documents
                .iter()
                .find(|d| d.id == "document-1")
                .unwrap()
        )
        .unwrap(),
        serde_json::to_value(&private_head).unwrap()
    );
    assert_eq!(
        serde_json::to_value(
            report
                .workspace
                .documents
                .iter()
                .find(|d| d.id == "unrelated-private")
                .unwrap()
        )
        .unwrap(),
        serde_json::to_value(&unrelated).unwrap()
    );
    assert_eq!(report.workspace.active_document, selected.document_id);
    assert!(report.displaced_selected.is_none());
    assert!(
        report
            .preserved_dependency_heads
            .contains(&pin(&local, "document-1"))
    );
    assert!(
        report
            .workspace
            .history
            .iter()
            .any(|e| Some(&e.id) == original_source.head_event.as_ref())
    );
    let fresh = run_material::reconstruct(&report.workspace, &material).unwrap();
    assert_eq!(fresh.report().status, Status::ConsistentUnderAssumptions);
    let repeated = importer::merge_document_closure(
        &report.workspace,
        &fetched,
        &selected,
        &Default::default(),
    )
    .unwrap();
    assert!(repeated.inserted_documents.is_empty());
    assert_eq!(json(&repeated.workspace), json(&report.workspace));
    let preserved_active = importer::merge_document_closure(
        &local,
        &fetched,
        &selected,
        &ImportOptions {
            activate_selected: false,
            ..Default::default()
        },
    )
    .unwrap();
    assert_eq!(
        preserved_active.workspace.active_document,
        "unrelated-private"
    );
}
#[test]
fn selected_conflict_is_explicit_and_checkout_returns_entire_displaced_branch() {
    let source = Workspace::from_documents(
        "Original",
        vec![Molecule::dna("Plasmid", "ACGTACGT", Topology::Linear).unwrap()],
        0,
    )
    .unwrap();
    let mut local = source.clone();
    let mut remote = source.clone();
    advance(&mut local, "document-1", "ACGTACGA");
    advance(&mut remote, "document-1", "ACGTACGC");
    let selected = pin(&remote, "document-1");
    let fetched = sharing::project_document_closure(&remote, &selected)
        .unwrap()
        .workspace;
    let before = json(&local);
    let expected = pin(&local, "document-1");
    let error = importer::merge_document_closure(&local, &fetched, &selected, &Default::default())
        .unwrap_err();
    assert!(
        matches!(error,ImportError::SelectedConflict{local:actual,incoming} if actual==expected&&incoming==selected)
    );
    assert_eq!(json(&local), before);
    let ff = ImportOptions {
        selected_policy: SelectedPolicy::FastForward {
            expected_local: expected.clone(),
        },
        activate_selected: true,
    };
    assert!(matches!(
        importer::merge_document_closure(&local, &fetched, &selected, &ff),
        Err(ImportError::NotFastForward { .. })
    ));
    let checkout = ImportOptions {
        selected_policy: SelectedPolicy::Checkout {
            expected_local: expected.clone(),
        },
        activate_selected: true,
    };
    let report = importer::merge_document_closure(&local, &fetched, &selected, &checkout).unwrap();
    let bookmark = report.displaced_selected.as_ref().unwrap();
    let serialized = serde_json::to_vec(bookmark).unwrap();
    let restored: importer::HeadBookmark = serde_json::from_slice(&serialized).unwrap();
    assert_eq!(
        serde_json::to_value(&restored.document).unwrap(),
        serde_json::to_value(&local.documents[0]).unwrap()
    );
    assert_eq!(restored.imported_selection, selected);
    assert_eq!(pin(&report.workspace, "document-1"), selected);
    assert!(report.workspace.documents[0].undo.is_empty());
    assert!(report.workspace.revisions.contains_key(&expected.revision));
    assert!(
        report
            .workspace
            .history
            .iter()
            .any(|e| Some(&e.id) == expected.head_event.as_ref())
    );
    assert_eq!(json(&local), before);
    let mut changed = local.clone();
    advance(&mut changed, "document-1", "ACGTACGG");
    assert!(matches!(
        importer::merge_document_closure(&changed, &fetched, &selected, &checkout),
        Err(ImportError::StaleLocal { .. })
    ));
}
#[test]
fn fastforward_requires_exact_expected_event_ancestry_even_when_sequence_matches() {
    let local = Workspace::from_documents(
        "Original",
        vec![Molecule::dna("Plasmid", "ACGTACGT", Topology::Linear).unwrap()],
        0,
    )
    .unwrap();
    let mut remote = local.clone();
    let expected = pin(&local, "document-1");
    advance(&mut remote, "document-1", "ACGTACGA");
    let selected = pin(&remote, "document-1");
    let fetched = sharing::project_document_closure(&remote, &selected)
        .unwrap()
        .workspace;
    let options = ImportOptions {
        selected_policy: SelectedPolicy::FastForward {
            expected_local: expected.clone(),
        },
        activate_selected: true,
    };
    let merged = importer::merge_document_closure(&local, &fetched, &selected, &options).unwrap();
    assert_eq!(pin(&merged.workspace, "document-1"), selected);
    assert!(merged.displaced_selected.is_some());
    let mut same_dna = local.clone();
    let d = &mut same_dna.documents[0];
    let alternate = RevisionEvent::new(
        &d.id,
        &d.head,
        vec![],
        "Different root event for same exact molecule",
        Transition::IndependentReference,
    )
    .unwrap();
    d.head_event = Some(alternate.id.clone());
    same_dna.history = vec![alternate];
    same_dna.validate().unwrap();
    let incoming = pin(&same_dna, "document-1");
    let fetched = sharing::project_document_closure(&same_dna, &incoming)
        .unwrap()
        .workspace;
    assert!(matches!(
        importer::merge_document_closure(&local, &fetched, &incoming, &Default::default()),
        Err(ImportError::SelectedConflict { .. })
    ));
    assert!(matches!(
        importer::merge_document_closure(&local, &fetched, &incoming, &options),
        Err(ImportError::NotFastForward { .. })
    ));
}
#[test]
fn immutable_attachment_collision_and_nonclosure_archive_leave_inputs_unchanged() {
    let (local, fetched, selected, _) = run_fixture();
    let mut collision = fetched.clone();
    let key = collision.attachments.keys().next().unwrap().clone();
    collision.attachments.get_mut(&key).unwrap().name =
        "A different retained local attachment name".into();
    collision.validate().unwrap();
    let before = json(&collision);
    let remote_before = json(&fetched);
    assert!(
        matches!(importer::merge_document_closure(&collision,&fetched,&selected,&Default::default()),Err(ImportError::Collision{kind:IdentityKind::Attachment,id}) if id==key)
    );
    assert_eq!(json(&collision), before);
    assert_eq!(json(&fetched), remote_before);
    let mut broader = fetched.clone();
    add(
        &mut broader,
        "unrelated-remote",
        "Unshared data",
        "TTGGAACC",
    );
    broader.validate().unwrap();
    assert!(matches!(
        importer::merge_document_closure(&local, &broader, &selected, &Default::default()),
        Err(ImportError::InvalidClosure)
    ));
    let mut invalid = fetched.clone();
    invalid
        .revisions
        .get_mut(&selected.revision)
        .unwrap()
        .name
        .push_str(" tampered");
    assert!(matches!(
        importer::merge_document_closure(&local, &invalid, &selected, &Default::default()),
        Err(ImportError::InvalidFetched)
    ));
}
