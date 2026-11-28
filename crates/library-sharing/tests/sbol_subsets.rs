use collagene_core::workspace::{self, Workspace, WorkspaceAttachment};
use collagene_design_exchange::{self as sbol, ProjectionPolicy};
use collagene_library_sharing::{
    self as sharing, Selection,
    sbol_subset::{self as review, PublicationChoices, Receipt},
};
use collagene_model::Topology;
use serde_json::{Value, json};
fn original() -> Vec<u8> {
    serde_json::to_vec(&json!({"@context":{"sbol":sbol::SBOL},"@graph":[
 {"@id":"https://example.org/c","@type":"sbol:Component","sbol:hasNamespace":{"@id":"https://example.org/"},"sbol:type":{"@id":sbol::DNA},"sbol:hasSequence":{"@id":"https://example.org/s"}},
 {"@id":"https://example.org/s","@type":"sbol:Sequence","sbol:hasNamespace":{"@id":"https://example.org/"},"sbol:encoding":{"@id":sbol::NUCLEIC_ENCODING},"sbol:elements":"AGTCCGATGCTACGTAGGACTTCAGACCTGTAGATCGTCA"},
 {"@id":"https://private.invalid/component","@type":"sbol:Component","sbol:hasNamespace":{"@id":"https://private.invalid/"},"sbol:name":"DO_NOT_SHARE_PRIVATE_SENTINEL","sbol:type":{"@id":sbol::DNA},"sbol:hasSequence":{"@id":"https://private.invalid/sequence"}},
 {"@id":"https://private.invalid/sequence","@type":"sbol:Sequence","sbol:hasNamespace":{"@id":"https://private.invalid/"},"sbol:encoding":{"@id":sbol::NUCLEIC_ENCODING},"sbol:elements":"TTTTCCCCAAAAGGGG"}
 ]})).unwrap()
}
fn setup() -> (Workspace, Selection, String) {
    let raw = original();
    let hash = blake3::hash(&raw).to_hex().to_string();
    let molecule = sbol::import(&raw)
        .unwrap()
        .project(
            &sbol::Selection {
                component: "https://example.org/c".into(),
                sequence: "https://example.org/s".into(),
                topology: Some(Topology::Linear),
            },
            ProjectionPolicy::Strict,
        )
        .unwrap()
        .molecule;
    let mut w = Workspace::from_documents("Private local", vec![molecule], 0).unwrap();
    w.attachments.insert(
        hash.clone(),
        WorkspaceAttachment {
            kind: sharing::SBOL_SOURCE_ATTACHMENT_KIND.into(),
            name: "private-original.jsonld".into(),
            media_type: "application/ld+json".into(),
            bytes: raw,
        },
    );
    let d = &w.documents[0];
    let selected = Selection {
        document_id: d.id.clone(),
        revision: d.head.clone(),
        head_event: d.head_event.clone(),
    };
    (w, selected, hash)
}
#[test]
fn reviewed_publish_preserves_native_identity_and_original_locally_while_excluding_private_graph() {
    let (w, s, original) = setup();
    let before = workspace::save_workspace(&w).unwrap();
    assert!(sharing::project_document_closure(&w, &s).is_err());
    let r = review::prepare_documentary_subset(&w, &s, &PublicationChoices::new()).unwrap();
    assert_eq!(r.entries().len(), 1);
    assert_eq!(r.entries()[0].omissions.unrelated_objects, 2);
    let out = review::document_closure_reviewed(&w, &s, &r.approve()).unwrap();
    assert!(!out.workspace.attachments.contains_key(&original));
    assert_eq!(out.workspace.revisions, w.revisions);
    assert_eq!(
        serde_json::to_value(&out.workspace.history).unwrap(),
        serde_json::to_value(&w.history).unwrap()
    );
    assert_eq!(out.workspace.documents[0].id, w.documents[0].id);
    let reopened = workspace::load_workspace(&out.archive).unwrap();
    assert!(sharing::inspect_document_closure(&reopened, &s).is_ok());
    for a in reopened.attachments.values() {
        let t = std::str::from_utf8(&a.bytes).unwrap();
        assert!(!t.contains("DO_NOT_SHARE_PRIVATE_SENTINEL"));
        assert!(!t.contains("private.invalid"));
    }
    assert_eq!(workspace::save_workspace(&w).unwrap(), before);
    let receipt: Receipt = serde_json::from_slice(
        &reopened
            .attachments
            .values()
            .find(|a| a.kind == review::RECEIPT_KIND)
            .unwrap()
            .bytes,
    )
    .unwrap();
    assert_eq!(receipt.original_hash, original);
    assert!(!receipt.original_bytes_included);
    assert_eq!(receipt.derivation, "unverified-documentary");
    assert_ne!(receipt.subset_hash, receipt.original_hash);
    assert_eq!(
        reopened.revisions[&s.revision].metadata["sbol:source_hash"],
        original
    );
    assert!(sharing::document_closure(&reopened, &s).is_ok());
}
#[test]
fn altered_review_context_and_receipt_claims_cannot_be_accepted() {
    let (w, s, _) = setup();
    let r = review::prepare_documentary_subset(&w, &s, &PublicationChoices::new()).unwrap();
    let mut changed = w.clone();
    changed.attachments.insert(
        "00".repeat(32),
        WorkspaceAttachment {
            kind: "sbol-design-exchange-v1".into(),
            name: "extra".into(),
            media_type: "application/ld+json".into(),
            bytes: vec![],
        },
    );
    assert!(review::project_document_closure_reviewed(&changed, &s, &r.approve()).is_err());
    let out = review::project_document_closure_reviewed(&w, &s, &r.approve()).unwrap();
    for field in [
        "original_bytes_included",
        "subset_policy",
        "original_hash",
        "request",
    ] {
        let mut bad = out.workspace.clone();
        let (key, mut a) = bad
            .attachments
            .iter()
            .find(|(_, a)| a.kind == review::RECEIPT_KIND)
            .map(|(k, a)| (k.clone(), a.clone()))
            .unwrap();
        bad.attachments.remove(&key);
        let mut v: Value = serde_json::from_slice(&a.bytes).unwrap();
        match field {
            "original_bytes_included" => v[field] = json!(true),
            "subset_policy" => v[field] = json!("unknown"),
            "original_hash" => v[field] = json!("ab".repeat(32)),
            _ => v["request"]["roots"] = json!({}),
        };
        // Encode via the typed schema, so failures test policy rather than map key ordering.
        let r: Receipt = serde_json::from_value(v).unwrap();
        a.bytes = serde_json::to_vec(&r).unwrap();
        let hash = blake3::hash(&a.bytes).to_hex().to_string();
        bad.attachments.insert(hash, a);
        assert!(
            sharing::inspect_document_closure(&bad, &s).is_err(),
            "{field}"
        );
    }
}
#[test]
fn extra_subset_object_and_overlapping_receipts_fail_closed() {
    let (w, s, _) = setup();
    let r = review::prepare_documentary_subset(&w, &s, &PublicationChoices::new()).unwrap();
    let mut out = review::project_document_closure_reviewed(&w, &s, &r.approve())
        .unwrap()
        .workspace;
    let (key, mut a) = out
        .attachments
        .iter()
        .find(|(_, a)| a.kind == review::RECEIPT_KIND)
        .map(|(k, a)| (k.clone(), a.clone()))
        .unwrap();
    let mut r: Receipt = serde_json::from_slice(&a.bytes).unwrap();
    r.omissions.unrelated_objects += 1;
    a.bytes = serde_json::to_vec(&r).unwrap();
    let key2 = blake3::hash(&a.bytes).to_hex().to_string();
    assert_ne!(key, key2);
    out.attachments.insert(key2, a);
    assert!(sharing::inspect_document_closure(&out, &s).is_err());
}

#[test]
fn documentary_projection_preserves_real_native_run_and_fresh_capsule_pcr_reuse() {
    use collagene_core::{
        analysis::Primer,
        run_lineage::{self, EventBinding, InputEventBinding},
        run_material,
        verification::{DesignPlan, Operation, PlanInput, Source, Status},
    };
    use collagene_model::Molecule;
    let (mut w, _, original) = setup();
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
    let review =
        review::prepare_documentary_subset(&recorded.workspace, &s, &PublicationChoices::new())
            .unwrap();
    let out =
        review::document_closure_reviewed(&recorded.workspace, &s, &review.approve()).unwrap();
    assert!(!out.workspace.attachments.contains_key(&original));
    assert_eq!(out.workspace.revisions[&stock.revision()], stock);
    let reopened = workspace::load_workspace(&out.archive).unwrap();
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
