use collagene_core::{
    documents::workspace::{self, Workspace, WorkspaceAttachment},
    material::{
        lineage::{self as run_lineage, EventBinding, InputEventBinding},
        run as run_material,
    },
    planning::plan::{DesignPlan, Operation, PlanInput, Source},
    sequence::analysis::Primer,
    verification::evidence::Status,
};
use collagene_design_exchange::{self as sbol, ProjectionPolicy};
use collagene_library_sharing::{self as sharing, Selection};
use collagene_model::{Molecule, Topology};
use serde_json::{Value, json};
const DNA: &str = "AGTCCGATGCTACGTAGGACTTCAGACCTGTAGATCGTCA";
const C: &str = "https://example.test/sbol/component";
const S: &str = "https://example.test/sbol/sequence";
const ALT: &str = "https://example.test/sbol/alternative";
fn graph() -> Value {
    json!({"@context":{"sbol":sbol::SBOL},"@graph":[
        {"@id":C,"@type":"sbol:Component","sbol:hasNamespace":{"@id":"https://example.test/sbol"},"sbol:name":"Imported literal stock","sbol:type":[{"@id":sbol::DNA},{"@id":sbol::LINEAR}],"sbol:hasSequence":{"@id":S}},
        {"@id":S,"@type":"sbol:Sequence","sbol:hasNamespace":{"@id":"https://example.test/sbol"},"sbol:encoding":{"@id":sbol::NUCLEIC_ENCODING},"sbol:elements":DNA}
    ]})
}
fn bytes(value: &Value) -> Vec<u8> {
    serde_json::to_vec(value).unwrap()
}
fn project(raw: &[u8], sequence: &str, policy: ProjectionPolicy) -> Molecule {
    sbol::import(raw)
        .unwrap()
        .project(
            &sbol::Selection {
                component: C.into(),
                sequence: sequence.into(),
                topology: Some(Topology::Linear),
            },
            policy,
        )
        .unwrap()
        .molecule
}
fn attach(w: &mut Workspace, raw: &[u8]) -> String {
    let hash = blake3::hash(raw).to_hex().to_string();
    w.attachments.insert(
        hash.clone(),
        WorkspaceAttachment {
            kind: sharing::SBOL_SOURCE_ATTACHMENT_KIND.into(),
            name: "original-design.sbol.jsonld".into(),
            media_type: "application/ld+json".into(),
            bytes: raw.to_vec(),
        },
    );
    hash
}
fn selected(w: &Workspace, index: usize) -> Selection {
    let d = &w.documents[index];
    Selection {
        document_id: d.id.clone(),
        revision: d.head.clone(),
        head_event: d.head_event.clone(),
    }
}
fn binding(w: &Workspace, index: usize) -> EventBinding {
    let d = &w.documents[index];
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

#[test]
fn exact_sbol_original_survives_real_saved_pcr_and_native_capsule_replay() {
    let raw = bytes(&graph());
    let stock = project(&raw, S, ProjectionPolicy::Strict);
    let independent = Molecule::dna("Independent intended PCR", DNA, Topology::Linear).unwrap();
    let private = Molecule::dna(
        "UNRELATED PRIVATE STOCK",
        "CCCCAAAAGGGGTTTT",
        Topology::Circular,
    )
    .unwrap();
    let mut w = Workspace::from_documents(
        "Local unpublished work",
        vec![stock.clone(), independent.clone(), private.clone()],
        0,
    )
    .unwrap();
    let source_hash = attach(&mut w, &raw);
    let unneeded = bytes(
        &json!({"@context":{"sbol":sbol::SBOL},"@graph":[{"@id":"https://example.test/sbol/private_activity","@type":"https://example.test/private/Activity","https://example.test/private/note":"UNRELATED PRIVATE SBOL"}]}),
    );
    let unneeded_hash = attach(&mut w, &unneeded);
    let plan = DesignPlan {
        name: "Native PCR from projected stock".into(),
        target_revision: independent.revision(),
        inputs: vec![PlanInput {
            id: "stock".into(),
            revision: stock.revision(),
            molecule: stock.clone(),
        }],
        steps: vec![pcr(Source::Input("stock".into()))],
    };
    let first = run_lineage::record_run(
        &w,
        &plan,
        binding(&w, 1),
        vec![InputEventBinding {
            input_id: "stock".into(),
            event: binding(&w, 0),
        }],
        false,
    )
    .unwrap();
    let product = first.receipt.products[0].clone();
    let material = run_material::Selection {
        run_id: first.receipt.run_id.clone(),
        step_index: 0,
        event: EventBinding {
            document_id: product.document_id.clone(),
            event_id: product.event_id.clone(),
            revision: product.revision.clone(),
        },
    };
    let before = run_material::reconstruct(&first.workspace, &material).unwrap();
    assert_eq!(
        before.report().status,
        Status::ConsistentUnderAssumptions,
        "{:?}",
        before.report().findings
    );
    let pin = Selection {
        document_id: product.document_id,
        revision: product.revision,
        head_event: Some(product.event_id),
    };
    let source_fingerprint = run_lineage::workspace_fingerprint(&first.workspace).unwrap();
    let export = sharing::document_closure(&first.workspace, &pin).unwrap();
    assert!(export.inventory.attachments.contains(&source_hash));
    assert!(!export.inventory.attachments.contains(&unneeded_hash));
    assert!(!export.inventory.revisions.contains(&private.revision()));
    assert_eq!(export.workspace.attachments[&source_hash].bytes, raw);
    assert_eq!(
        run_lineage::workspace_fingerprint(&first.workspace).unwrap(),
        source_fingerprint
    );
    let reopened = workspace::load_workspace(&export.archive).unwrap();
    assert_eq!(reopened.revisions[&stock.revision()], stock);
    let after = run_material::reconstruct(&reopened, &material).unwrap();
    assert_eq!(before.report().status, after.report().status);
    assert_eq!(before.selection(), after.selection());
    assert_eq!(before.reference(), after.reference());
    assert_eq!(before.tracked(), after.tracked());
    assert_eq!(after.reference().sequence.as_str(), DNA);
    assert_eq!(after.tracked().duplex.top.as_str(), DNA);
    assert_eq!(
        after.tracked().duplex.left,
        collagene_core::chemistry::Terminus::PCR
    );
    assert_eq!(
        after.tracked().duplex.right,
        collagene_core::chemistry::Terminus::PCR
    );
    assert!(after.tracked().duplex.nicks.is_empty());
    let capsule = collagene_core::material::capsule::capsule(&reopened, &material).unwrap();
    let prepared = reopened.revisions[&material.event.revision].clone();
    let downstream = DesignPlan {
        name: "Fresh use of retained PCR stock".into(),
        target_revision: independent.revision(),
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
    let verified = collagene_core::verification::engine::verify(&downstream, &independent);
    assert_eq!(
        verified.status,
        Status::ConsistentUnderAssumptions,
        "{:?}",
        verified.findings
    );
    assert_eq!(verified.products.last().unwrap().sequence.as_str(), DNA);
}

#[test]
fn unrelated_original_graph_objects_refuse_without_whole_file_fallback() {
    let mut value = graph();
    value["@graph"].as_array_mut().unwrap().push(json!({"@id":"https://example.test/sbol/private","@type":"sbol:Component","sbol:hasNamespace":{"@id":"https://example.test/sbol"},"sbol:name":"PRIVATE OTHER COMPONENT","sbol:type":{"@id":sbol::DNA},"sbol:hasSequence":{"@id":"https://example.test/sbol/private_sequence"}}));
    value["@graph"].as_array_mut().unwrap().push(json!({"@id":"https://example.test/sbol/private_sequence","@type":"sbol:Sequence","sbol:hasNamespace":{"@id":"https://example.test/sbol"},"sbol:encoding":{"@id":sbol::NUCLEIC_ENCODING},"sbol:elements":"CCCCAAAAGGGGTTTT"}));
    let raw = bytes(&value);
    let molecule = project(&raw, S, ProjectionPolicy::Strict);
    let mut w = Workspace::from_documents("Private graph", vec![molecule], 0).unwrap();
    let hash = attach(&mut w, &raw);
    let original = run_lineage::workspace_fingerprint(&w).unwrap();
    assert!(
        matches!(sharing::document_closure(&w,&selected(&w,0)),Err(sharing::Error::SbolSource(reason))if reason.contains("unrelated graph objects"))
    );
    assert_eq!(run_lineage::workspace_fingerprint(&w).unwrap(), original);
    assert_eq!(w.attachments[&hash].bytes, raw);
}

#[test]
fn alternative_sequences_require_explicit_retained_native_input_selections() {
    let mut value = graph();
    value["@graph"][0]["sbol:hasSequence"] = json!([{"@id":S},{"@id":ALT}]);
    value["@graph"].as_array_mut().unwrap().push(json!({"@id":ALT,"@type":"sbol:Sequence","sbol:hasNamespace":{"@id":"https://example.test/sbol"},"sbol:encoding":{"@id":sbol::NUCLEIC_ENCODING},"sbol:elements":DNA}));
    let raw = bytes(&value);
    let a = project(&raw, S, ProjectionPolicy::AllowOmissions);
    let b = project(&raw, ALT, ProjectionPolicy::AllowOmissions);
    let target = Molecule::dna("Independent PCR target", DNA, Topology::Linear).unwrap();
    let mut w = Workspace::from_documents(
        "Two local representations",
        vec![a.clone(), b.clone(), target.clone()],
        0,
    )
    .unwrap();
    let hash = attach(&mut w, &raw);
    // Merely retaining a second local document is not selected publication scope.
    assert!(
        matches!(sharing::document_closure(&w,&selected(&w,0)),Err(sharing::Error::SbolSource(reason))if reason.contains("alternative"))
    );
    let plan = DesignPlan {
        name: "Both declared representations used as original PCR inputs".into(),
        target_revision: target.revision(),
        inputs: vec![
            PlanInput {
                id: "a".into(),
                revision: a.revision(),
                molecule: a,
            },
            PlanInput {
                id: "b".into(),
                revision: b.revision(),
                molecule: b,
            },
        ],
        steps: vec![
            pcr(Source::Input("a".into())),
            pcr(Source::Input("b".into())),
        ],
    };
    let kept = run_lineage::record_run(
        &w,
        &plan,
        binding(&w, 2),
        vec![
            InputEventBinding {
                input_id: "a".into(),
                event: binding(&w, 0),
            },
            InputEventBinding {
                input_id: "b".into(),
                event: binding(&w, 1),
            },
        ],
        false,
    )
    .unwrap();
    let product = &kept.receipt.products[1];
    let pin = Selection {
        document_id: product.document_id.clone(),
        revision: product.revision.clone(),
        head_event: Some(product.event_id.clone()),
    };
    let export = sharing::document_closure(&kept.workspace, &pin).unwrap();
    assert_eq!(export.workspace.attachments[&hash].bytes, raw);
    assert!(export.inventory.revisions.contains(&w.documents[0].head));
    assert!(export.inventory.revisions.contains(&w.documents[1].head));
}

#[test]
fn selected_unknown_properties_classes_and_missing_documentary_pins_are_not_guessed() {
    let mut property = graph();
    property["@graph"][0]["https://example.test/private/relationship"] =
        json!({"@id":"https://example.test/private/unknown"});
    let mut classes = graph();
    classes["@graph"][0]["@type"] = json!([
        "sbol:Component",
        "https://example.test/private/AdditionalType"
    ]);
    for value in [property, classes] {
        let raw = bytes(&value);
        let m = project(&raw, S, ProjectionPolicy::AllowOmissions);
        let mut w =
            Workspace::from_documents("Unsupported reference contract", vec![m], 0).unwrap();
        attach(&mut w, &raw);
        assert!(matches!(
            sharing::document_closure(&w, &selected(&w, 0)),
            Err(sharing::Error::SbolSource(_))
        ));
    }
    let raw = bytes(&graph());
    let m = project(&raw, S, ProjectionPolicy::Strict);
    let missing = Workspace::from_documents("No original attached", vec![m.clone()], 0).unwrap();
    assert!(matches!(
        sharing::document_closure(&missing, &selected(&missing, 0)),
        Err(sharing::Error::SbolSource(_))
    ));
    let mut wrong = m;
    wrong.metadata.insert(
        "sbol:sequence_iri".into(),
        "https://example.test/sbol/not_declared".into(),
    );
    let mut w = Workspace::from_documents("Wrong documentary selection", vec![wrong], 0).unwrap();
    attach(&mut w, &raw);
    assert!(matches!(
        sharing::document_closure(&w, &selected(&w, 0)),
        Err(sharing::Error::SbolSource(_))
    ));
}

#[test]
fn unrelated_sbol_attachment_does_not_block_legacy_publication_and_other_unknown_kinds_still_do() {
    let normal = Molecule::dna("Ordinary local plasmid", DNA, Topology::Linear).unwrap();
    let mut w = Workspace::from_documents("Local private documentary collections", vec![normal], 0)
        .unwrap();
    let raw = bytes(&graph());
    let hash = attach(&mut w, &raw);
    let export = sharing::document_closure(&w, &selected(&w, 0)).unwrap();
    assert!(!export.inventory.attachments.contains(&hash));
    assert!(export.workspace.attachments.is_empty());
    let unknown = WorkspaceAttachment {
        kind: "unknown-foreign-graph-v1".into(),
        name: "unknown".into(),
        media_type: "application/json".into(),
        bytes: b"{}".to_vec(),
    };
    let hash = blake3::hash(&unknown.bytes).to_hex().to_string();
    w.attachments.insert(hash, unknown);
    assert!(matches!(
        sharing::document_closure(&w, &selected(&w, 0)),
        Err(sharing::Error::UnsupportedAttachment(_))
    ));
}

#[test]
fn documentary_source_hash_is_not_a_fresh_sequence_match_or_material_bootstrap() {
    let raw = bytes(&graph());
    let mut edited = project(&raw, S, ProjectionPolicy::Strict);
    edited.sequence = collagene_model::Sequence::dna("TTTTGGGGCCCCAAAA").unwrap();
    edited.name = "Edited after SBOL projection".into();
    let hash = edited.metadata["sbol:source_hash"].clone();
    let revision = edited.revision();
    let mut w = Workspace::from_documents("Edited native stock", vec![edited.clone()], 0).unwrap();
    attach(&mut w, &raw);
    let export = sharing::document_closure(&w, &selected(&w, 0)).unwrap();
    assert_eq!(export.workspace.revisions[&revision], edited);
    assert_eq!(
        export.workspace.revisions[&revision].metadata["sbol:source_hash"],
        hash
    );
    assert_eq!(export.workspace.attachments[&hash].bytes, raw);
    assert!(
        export
            .workspace
            .attachments
            .values()
            .all(|a| a.kind != collagene_core::material::lineage::WORKSPACE_ATTACHMENT_KIND)
    );
}
