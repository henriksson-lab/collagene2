//! Native-local reopen retains exact documentary bytes without public authority.
use collagene_core::workspace::{self, RevisionEvent, Workspace, WorkspaceAttachment};
use collagene_design_exchange as sbol;
use collagene_library_sharing::{
    self as sharing, Selection,
    import::{self as importer, ImportError, ImportOptions, SelectedPolicy},
};
use collagene_model::{Molecule, Topology};
use serde_json::json;
fn fixture() -> (Workspace, Selection, String) {
    let raw=serde_json::to_vec(&json!({"@context":{"sbol":sbol::SBOL},"@graph":[
 {"@id":"https://example.org/c","@type":"sbol:Component","sbol:hasNamespace":{"@id":"https://example.org/"},"sbol:type":{"@id":sbol::DNA},"sbol:hasSequence":{"@id":"https://example.org/s"}},
 {"@id":"https://example.org/s","@type":"sbol:Sequence","sbol:hasNamespace":{"@id":"https://example.org/"},"sbol:encoding":{"@id":sbol::NUCLEIC_ENCODING},"sbol:elements":"AGTCCGATGCTACGTAGGACTTCAGACCTGTAGATCGTCA"},
 {"@id":"https://private.invalid/note","@type":"https://private.invalid/Activity","https://private.invalid/note":"EXACT-LOCAL-ORIGINAL-MUST-SURVIVE"}
 ]})).unwrap();
    let hash = blake3::hash(&raw).to_hex().to_string();
    let molecule = sbol::import(&raw)
        .unwrap()
        .project(
            &sbol::Selection {
                component: "https://example.org/c".into(),
                sequence: "https://example.org/s".into(),
                topology: Some(Topology::Linear),
            },
            sbol::ProjectionPolicy::Strict,
        )
        .unwrap()
        .molecule;
    let mut w = Workspace::from_documents(
        "Full private workspace",
        vec![
            molecule,
            Molecule::dna(
                "Unrelated catalog snapshot document",
                "CCCGGGTTTAAA",
                Topology::Circular,
            )
            .unwrap(),
        ],
        0,
    )
    .unwrap();
    w.attachments.insert(
        hash.clone(),
        WorkspaceAttachment {
            kind: sharing::SBOL_SOURCE_ATTACHMENT_KIND.into(),
            name: "Exact local source.jsonld".into(),
            media_type: "application/ld+json".into(),
            bytes: raw,
        },
    );
    let d = &w.documents[0];
    let pin = Selection {
        document_id: d.id.clone(),
        revision: d.head.clone(),
        head_event: d.head_event.clone(),
    };
    (w, pin, hash)
}
#[test]
fn private_reopen_keeps_exact_original_but_restores_only_selected_native_dependency_graph() {
    let (w, pin, hash) = fixture();
    let before = serde_json::to_value(&w).unwrap();
    assert!(sharing::project_document_closure(&w, &pin).is_err());
    let out = sharing::project_private_document_closure(&w, &pin).unwrap();
    assert_eq!(out.workspace.documents.len(), 1);
    assert_eq!(out.workspace.revisions.len(), 1);
    assert_eq!(
        out.workspace.attachments[&hash].bytes,
        w.attachments[&hash].bytes
    );
    assert_eq!(
        out.workspace.revisions[&pin.revision],
        w.revisions[&pin.revision]
    );
    assert_eq!(out.workspace.documents[0].head_event, pin.head_event);
    let mut unrelated = Workspace::from_documents(
        "Live current editor",
        vec![Molecule::dna("Keep editor branch", "GGGGAAAACCCC", Topology::Linear).unwrap()],
        0,
    )
    .unwrap();
    // Native document IDs are chosen by Workspace deterministically; use a stable
    // independent local document identity to exercise a real union, not a collision.
    let local_id = "independent-current-editor".to_string();
    let local_m = unrelated.revisions[&unrelated.documents[0].head].clone();
    // Rebrand through native history constructor instead of rewriting immutable events.
    let history = collagene_core::history::DocumentHistory::new(&local_id, local_m).unwrap();
    unrelated.documents.clear();
    unrelated.history.clear();
    unrelated.revisions.clear();
    unrelated.active_document = local_id.clone();
    history.contribute(&mut unrelated).unwrap();
    let local_before = serde_json::to_value(&unrelated).unwrap();
    let merged = importer::merge_private_document_closure(
        &unrelated,
        &out.workspace,
        &pin,
        &ImportOptions::default(),
    )
    .unwrap();
    assert_eq!(merged.workspace.documents.len(), 2);
    assert!(merged.workspace.documents.iter().any(|d| d.id == local_id));
    assert!(merged.workspace.attachments.contains_key(&hash));
    assert!(
        importer::merge_document_closure(
            &unrelated,
            &out.workspace,
            &pin,
            &ImportOptions::default()
        )
        .is_err()
    );
    assert_eq!(serde_json::to_value(&unrelated).unwrap(), local_before);
    assert_eq!(serde_json::to_value(&w).unwrap(), before);
    let archive = workspace::save_workspace(&merged.workspace).unwrap();
    let reopened = workspace::load_workspace(&archive).unwrap();
    assert_eq!(
        reopened.attachments[&hash].bytes,
        w.attachments[&hash].bytes
    );
}
#[test]
fn private_import_retains_public_collision_head_and_stale_checkout_guards() {
    let (w, pin, _) = fixture();
    let incoming = sharing::project_private_document_closure(&w, &pin)
        .unwrap()
        .workspace;
    let mut local = incoming.clone();
    let mut edited = local.revisions[&pin.revision].clone();
    edited.name = "Unpublished local edit".into();
    let revision = edited.revision();
    let event = RevisionEvent::new(
        &pin.document_id,
        &revision,
        vec![pin.head_event.clone().unwrap()],
        "Local edit",
    )
    .unwrap();
    local.revisions.insert(revision.clone(), edited);
    local.history.push(event.clone());
    local.documents[0].head = revision;
    local.documents[0].head_event = Some(event.id);
    let before = serde_json::to_value(&local).unwrap();
    assert!(matches!(
        importer::merge_private_document_closure(
            &local,
            &incoming,
            &pin,
            &ImportOptions::default()
        ),
        Err(ImportError::SelectedConflict { .. })
    ));
    let stale = ImportOptions {
        selected_policy: SelectedPolicy::Checkout {
            expected_local: pin.clone(),
        },
        activate_selected: true,
    };
    assert!(matches!(
        importer::merge_private_document_closure(&local, &incoming, &pin, &stale),
        Err(ImportError::StaleLocal { .. })
    ));
    assert_eq!(serde_json::to_value(&local).unwrap(), before);
}

#[test]
fn private_open_preserves_real_saved_run_and_fresh_capsule_pcr_reuse() {
    use collagene_core::{
        analysis::Primer,
        run_lineage::{self, EventBinding, InputEventBinding},
        run_material,
        verification::{DesignPlan, Operation, PlanInput, Source, Status},
    };
    use collagene_model::Molecule;
    let (mut w, _, original) = fixture();
    let stock = w.revisions[&w.documents[0].head].clone();
    let target = Molecule::dna(
        "Independent PCR target",
        stock.sequence.as_str(),
        Topology::Linear,
    )
    .unwrap();
    // Rebuild only the native document list; the SBOL source hash and original
    // molecule are unchanged and the unrelated graph remains local.
    let attachments = w.attachments.clone();
    w = Workspace::from_documents("Local run", vec![stock.clone(), target.clone()], 0).unwrap();
    w.attachments = attachments;
    let bind = |i: usize| {
        let d = &w.documents[i];
        EventBinding {
            document_id: d.id.clone(),
            revision: d.head.clone(),
            event_id: d.head_event.clone().unwrap(),
        }
    };
    let pcr = |template| Operation::Pcr {
        template,
        forward: Primer::new("F", "AGTCCGAT", "").unwrap(),
        reverse: Primer::new("R", "TGACGATC", "").unwrap(),
        template_strand: collagene_core::chemistry::Strand::Top,
    };
    let plan = DesignPlan {
        name: "Actual PCR".into(),
        target_revision: target.revision(),
        inputs: vec![PlanInput {
            id: "stock".into(),
            revision: stock.revision(),
            molecule: stock.clone(),
        }],
        steps: vec![pcr(Source::Input("stock".into()))],
    };
    let recorded = run_lineage::record_run(
        &w,
        &plan,
        bind(1),
        vec![InputEventBinding {
            input_id: "stock".into(),
            event: bind(0),
        }],
        false,
    )
    .unwrap();
    let product = &recorded.receipt.products[0];
    let material = run_material::Selection {
        run_id: recorded.receipt.run_id.clone(),
        step_index: 0,
        event: EventBinding {
            document_id: product.document_id.clone(),
            event_id: product.event_id.clone(),
            revision: product.revision.clone(),
        },
    };
    let s = Selection {
        document_id: product.document_id.clone(),
        revision: product.revision.clone(),
        head_event: Some(product.event_id.clone()),
    };
    let before = run_material::reconstruct(&recorded.workspace, &material).unwrap();
    let out = sharing::project_private_document_closure(&recorded.workspace, &s).unwrap();
    assert!(out.workspace.attachments.contains_key(&original));
    assert_eq!(out.workspace.revisions[&stock.revision()], stock);
    let reopened =
        workspace::load_workspace(&workspace::save_workspace(&out.workspace).unwrap()).unwrap();
    let after = run_material::reconstruct(&reopened, &material).unwrap();
    assert_eq!(before.selection(), after.selection());
    assert_eq!(before.reference(), after.reference());
    assert_eq!(before.tracked(), after.tracked());
    let prepared = reopened.revisions[&material.event.revision].clone();
    let capsule = run_material::capsule(&reopened, &material).unwrap();
    let downstream = DesignPlan {
        name: "Fresh physical reuse".into(),
        target_revision: target.revision(),
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
    let verified = collagene_core::verification::verify(&downstream, &target);
    assert_eq!(
        verified.status,
        Status::ConsistentUnderAssumptions,
        "{:?}",
        verified.findings
    );
}
