//! In-process actual SDK flow, with no HTTP, filesystem or host-catalog writes.
use base64::{Engine as _, engine::general_purpose::STANDARD};
use collagene_core::{
    cloning::restriction::catalog as catalog_digest,
    documents::workspace::{self, Workspace},
    material::{
        derivation as native,
        lineage::EventBinding,
        run::{self as run_material, MaterialOrigin},
    },
    sequence::analysis::Primer,
};
use collagene_mcp::server::CollageneServer;
use collagene_mcp::session::references::MaterialDerivationOrigin as Origin;
use collagene_mcp::tools::workspace::MaterialDerivationInstallRequest as InstallRequest;
use collagene_mcp::tools::workspace::MaterialDerivationOperation as Derivation;
use collagene_mcp::tools::workspace::MaterialPcrChoice;
use collagene_model::{Molecule, Topology};
use rmcp::{RoleClient, ServiceExt, model::CallToolRequestParams, service::RunningService};
use serde_json::{Value, json};
type Client = RunningService<RoleClient, ()>;
const STOCK: &str = "AGTCCGATGAATTCCGTAGGACTTCAGACCTGTAGATCGTCA";
async fn setup() -> (Client, tokio::task::JoinHandle<()>) {
    let (a, b) = tokio::io::duplex(65536);
    let server = tokio::spawn(async move {
        CollageneServer::default()
            .serve(a)
            .await
            .unwrap()
            .waiting()
            .await
            .unwrap();
    });
    (().serve(b).await.unwrap(), server)
}
async fn call(c: &Client, name: &str, args: Value) -> Value {
    let reply = c
        .call_tool(
            CallToolRequestParams::new(name.to_owned())
                .with_arguments(args.as_object().unwrap().clone()),
        )
        .await
        .unwrap();
    assert_ne!(reply.is_error, Some(true), "{name}: {:?}", reply.content);
    reply.structured_content.unwrap()
}
async fn reject(c: &Client, args: Value) {
    match c
        .call_tool(
            CallToolRequestParams::new("material_derivation")
                .with_arguments(args.as_object().unwrap().clone()),
        )
        .await
    {
        Err(_) => {}
        Ok(reply) => assert_eq!(
            reply.is_error,
            Some(true),
            "Unexpected success: {:?}",
            reply
        ),
    }
}
async fn artifact(c: &Client, id: &Value) -> (String, String) {
    let (mut at, mut text) = (0, String::new());
    let encoding = loop {
        let page = call(
            c,
            "read_artifact",
            json!({"artifact_id":id,"start":at,"max_bytes":16384}),
        )
        .await;
        text.push_str(page["text"].as_str().unwrap());
        if page["done"] == true {
            break page["encoding"].as_str().unwrap().to_owned();
        }
        at = page["next_start"].as_u64().unwrap();
    };
    (text, encoding)
}
async fn json_artifact(c: &Client, id: &Value) -> Value {
    let (text, encoding) = artifact(c, id).await;
    assert_eq!(encoding, "utf-8");
    serde_json::from_str(&text).unwrap()
}
async fn native_archive(c: &Client, id: &Value) -> Workspace {
    let (text, encoding) = artifact(c, id).await;
    assert_eq!(encoding, "base64");
    workspace::load_workspace(&STANDARD.decode(text).unwrap()).unwrap()
}
fn reference(reply: &Value) -> Value {
    json!({"id":reply["workspace"]["id"],"revision":reply["workspace"]["revision"]})
}
fn event(w: &Workspace, index: usize) -> EventBinding {
    let d = &w.documents[index];
    EventBinding {
        document_id: d.id.clone(),
        event_id: d.head_event.clone().unwrap(),
        revision: d.head.clone(),
    }
}
fn fixture() -> Workspace {
    Workspace::from_documents(
        "Immutable selected source plus another active tab",
        vec![
            Molecule::dna("Original source", STOCK, Topology::Linear).unwrap(),
            Molecule::dna("Unrelated private document", "TTGACCGATC", Topology::Linear).unwrap(),
        ],
        1,
    )
    .unwrap()
}
// Descriptions may discuss geometry; only schema properties admit caller payloads.
fn schema_has_property(value: &serde_json::Value, name: &str) -> bool {
    match value {
        serde_json::Value::Object(fields) => {
            fields
                .get("properties")
                .and_then(serde_json::Value::as_object)
                .is_some_and(|properties| properties.contains_key(name))
                || fields
                    .values()
                    .any(|child| schema_has_property(child, name))
        }
        serde_json::Value::Array(items) => {
            items.iter().any(|child| schema_has_property(child, name))
        }
        _ => false,
    }
}
#[tokio::test]
async fn ordinary_pcr_native_install_redigest_and_derived_discovery_preserve_immutable_sources() {
    let (c, server) = setup().await;
    let tools = c.list_tools(Default::default()).await.unwrap().tools;
    assert_eq!(tools.len(), collagene_mcp::CURRENT_TOOL_COUNT);
    assert!(tools.iter().any(|t| t.name == "material_derivation"));
    let schema = serde_json::to_value(
        &tools
            .iter()
            .find(|t| t.name == "material_derivation")
            .unwrap()
            .input_schema,
    )
    .unwrap();
    assert!(
        !schema_has_property(&schema, "record_bytes"),
        "Closed capsules are not MCP arguments"
    );
    assert!(
        schema_has_property(&schema, "expected_product_revision"),
        "OE selection binds a native product revision without accepting geometry"
    );
    assert!(
        !schema_has_property(&schema, "geometry"),
        "Amplicon geometry is resolved in-process"
    );
    let original = fixture();
    let input_bytes = workspace::save_workspace(&original).unwrap();
    let pin = event(&original, 0);
    let imported = call(
        &c,
        "import_workspace_archive",
        json!({"data_base64":STANDARD.encode(&input_bytes)}),
    )
    .await;
    let input = reference(&imported);
    let before_count = call(&c, "list_documents", json!({})).await["total"].clone();
    let forward = Primer::new("Literal F", "AGTCCGAT", "").unwrap();
    let reverse = Primer::new("Literal R", "TGACGATC", "").unwrap();
    let preview=call(&c,"material_derivation",json!({"workspace":input,"task":{"PreviewPcr":{"template_strand":"Top","source":pin,"origin":"OriginalReference","forward":forward,"reverse":reverse}}})).await;
    assert_eq!(preview["count"], 1);
    assert_eq!(
        call(&c, "list_documents", json!({})).await["total"],
        before_count
    );
    assert_eq!(preview["previews"][0]["reusable"], true);
    let request = InstallRequest {
        sources: vec![collagene_mcp::tools::workspace::MaterialDerivationSource {
            role: collagene_core::material::derivation::SourceRole::Template,
            event: pin.clone(),
            origin: Origin::OriginalReference,
        }],
        derivation: Derivation::ExactPcr {
            forward: forward.clone(),
            reverse: reverse.clone(),
            selection: MaterialPcrChoice {
                selected_index: preview["previews"][0]["index"].as_u64().unwrap() as usize,
                expected_preview_revision: preview["previews"][0]["preview_revision"]
                    .as_str()
                    .unwrap()
                    .into(),
            },
            template_strand: collagene_core::chemistry::Strand::Top,
        },
        name: "Standalone native PCR".into(),
        save_primers: true,
    };
    let small = serde_json::to_string(&request).unwrap();
    assert!(small.len() < 4096);
    assert!(!small.contains("record_bytes"));
    assert!(!small.contains("geometry"));
    assert!(!small.contains(STOCK));
    let mut bad_request = request.clone();
    let Derivation::ExactPcr { selection, .. } = &mut bad_request.derivation else {
        unreachable!()
    };
    selection.expected_preview_revision = "0".repeat(64);
    reject(
        &c,
        json!({"workspace":input,"task":{"Install":{"request":bad_request}}}),
    )
    .await;
    assert_eq!(
        call(&c, "list_documents", json!({})).await["total"],
        before_count
    );
    let mut missing_index = request.clone();
    let Derivation::ExactPcr { selection, .. } = &mut missing_index.derivation else {
        unreachable!()
    };
    selection.selected_index = 99;
    reject(
        &c,
        json!({"workspace":input,"task":{"Install":{"request":missing_index}}}),
    )
    .await;
    let mut stale = input.clone();
    stale["revision"] = json!("0".repeat(64));
    reject(&c,json!({"workspace":stale,"task":{"PreviewPcr":{"template_strand":"Top","source":pin,"origin":"OriginalReference","forward":forward,"reverse":reverse}}})).await;
    let mut injected = serde_json::to_value(&request).unwrap();
    injected["derivation"]["ExactPcr"]["selection"]["geometry"] = json!({"sequence":STOCK});
    reject(
        &c,
        json!({"workspace":input,"task":{"Install":{"request":injected}}}),
    )
    .await;
    let installed = call(
        &c,
        "material_derivation",
        json!({"workspace":input,"task":{"Install":{"request":request}}}),
    )
    .await;
    assert_ne!(reference(&installed), input);
    assert_eq!(installed["source_workspace"], input);
    assert!(
        installed["authority"]
            .as_str()
            .unwrap()
            .contains("No independent target")
    );
    let receipt: native::Receipt = serde_json::from_value(installed["receipt"].clone()).unwrap();
    assert!(receipt.primer_attachment.is_some());
    let w = native_archive(&c, &installed["native_artifact_id"]).await;
    receipt.validate_archive(&w).unwrap();
    assert_eq!(w.documents.len(), 3);
    assert!(
        collagene_core::material::lineage::list_runs(&w)
            .unwrap()
            .is_empty()
    );
    assert_eq!(w.documents[0].head, original.documents[0].head);
    assert_eq!(w.documents[1].head_event, original.documents[1].head_event);
    let selected = run_material::MaterialSelection::Derived(receipt.product.clone());
    let inspection = call(
        &c,
        "inspect_material_origin",
        json!({"workspace":reference(&installed),"event":receipt.product,"selection":selected}),
    )
    .await;
    assert_eq!(inspection["origin_available"], true);
    let full = json_artifact(&c, &inspection["artifact_id"]).await;
    let inspection: run_material::MaterialInspection = serde_json::from_value(full).unwrap();
    let closed = inspection.selected.unwrap();
    assert!(matches!(closed, MaterialOrigin::Derived(_)));
    closed.validate_for_workspace(&w, &receipt.product).unwrap();
    let origin = run_material::MaterialSelection::Derived(receipt.product.clone());
    // Even a handle with the right event must name its actual origin family.
    reject(&c,json!({"workspace":reference(&installed),"task":{"PreviewPcr":{"template_strand":"Top","source":receipt.product,"origin":{"Saved":run_material::MaterialSelection::Run(run_material::Selection { run_id: "0".repeat(64), step_index: 0, event: receipt.product.clone() })},"forward":forward,"reverse":reverse}}})).await;
    let catalog = catalog_digest::builtin_catalog().unwrap();
    let eco = catalog
        .edition()
        .enzymes
        .iter()
        .find(|e| e.name == "EcoRI")
        .unwrap();
    let selection = catalog_digest::select(&catalog, &[eco.id.clone()]).unwrap();
    let digest=call(&c,"material_derivation",json!({"workspace":reference(&installed),"task":{"PreviewDigest":{"source":receipt.product,"origin":{"Saved":origin},"selection":selection,"options":collagene_core::cloning::restriction::inspection::InspectionOptions::default()}}})).await;
    assert_eq!(digest["products_determined"], true);
    assert_eq!(digest["fragment_count"], 2);
    let full = json_artifact(&c, &digest["artifact_id"]).await;
    let report: native::digest::DigestPreview = serde_json::from_value(full).unwrap();
    report
        .validate_binding(digest["context_revision"].as_str().unwrap())
        .unwrap();
    assert_eq!(report.fragment_modifications.len(), 2);
    let request = InstallRequest {
        sources: vec![collagene_mcp::tools::workspace::MaterialDerivationSource {
            role: collagene_core::material::derivation::SourceRole::Template,
            event: receipt.product.clone(),
            origin: Origin::Saved(origin.clone()),
        }],
        derivation: Derivation::CatalogDigest {
            conditions: None,
            selection,
            options: Default::default(),
            fragment_index: 1,
        },
        name: "Native right cut fragment".into(),
        save_primers: false,
    };
    let cut = call(
        &c,
        "material_derivation",
        json!({"workspace":reference(&installed),"task":{"Install":{"request":request}}}),
    )
    .await;
    let cut_receipt: native::Receipt = serde_json::from_value(cut["receipt"].clone()).unwrap();
    assert!(cut_receipt.primer_attachment.is_none());
    let cut_w = native_archive(&c, &cut["native_artifact_id"]).await;
    cut_receipt.validate_archive(&cut_w).unwrap();
    let fresh = native::reconstruct(&cut_w, &cut_receipt.product).unwrap();
    assert_eq!(
        fresh.reference().sequence.as_str(),
        "AATTCCGTAGGACTTCAGACCTGTAGATCGTCA"
    );
    assert_eq!(fresh.tracked().duplex.bottom_start, 4);
    let inspection=call(&c,"inspect_material_origin",json!({"workspace":reference(&cut),"event":cut_receipt.product,"selection":run_material::MaterialSelection::Derived(cut_receipt.product.clone())})).await;
    assert_eq!(inspection["origin_available"], true);
    reject(&c,json!({"workspace":input,"task":{"PreviewPcr":{"template_strand":"Top","source":pin,"origin":{"Saved":origin},"forward":forward,"reverse":reverse}}})).await;
    let mut wrong = pin.clone();
    wrong.event_id = event(&original, 1).event_id;
    reject(&c,json!({"workspace":input,"task":{"PreviewPcr":{"template_strand":"Top","source":wrong,"origin":"OriginalReference","forward":forward,"reverse":reverse}}})).await;
    reject(&c,json!({"workspace":input,"task":{"PreviewPcr":{"template_strand":"Top","source":pin,"origin":"OriginalReference","forward":forward,"reverse":reverse,"claimed_duplex":{}}}})).await;
    let exported = call(&c, "export_workspace_archive", input.clone()).await;
    let kept = native_archive(&c, &exported["artifact_id"]).await;
    assert_eq!(workspace::save_workspace(&kept).unwrap(), input_bytes);
    c.cancel().await.unwrap();
    server.await.unwrap();
}

#[tokio::test]
async fn mismatch_and_saved_fragment_repair_use_same_target_free_tool_and_exact_context() {
    use collagene_core::{chemistry::EndRepairMethod, cloning::pcr::mutagenesis::PcrOptions};
    let (c, server) = setup().await;
    let original = fixture();
    let bytes = workspace::save_workspace(&original).unwrap();
    let imported = call(
        &c,
        "import_workspace_archive",
        json!({"data_base64":STANDARD.encode(&bytes)}),
    )
    .await;
    let input = reference(&imported);
    let pin = event(&original, 0);
    let forward = Primer::new("Mismatch F", "AGTTCGAT", "CC").unwrap();
    let reverse = Primer::new("Exact R", "TGACGATC", "GG").unwrap();
    let mut options = PcrOptions::default();
    options.binding.min_three_prime_exact = 4;
    let args = json!({"workspace":input,"task":{"PreviewMismatchPcr":{"source":pin,"origin":"OriginalReference","forward":forward,"reverse":reverse,"options":options}}});
    let preview = call(&c, "material_derivation", args.clone()).await;
    assert_eq!(preview["count"], 1);
    assert_eq!(preview["reusable"], true);
    assert_eq!(preview["previews"][0]["substitution_count"], 1);
    let full = json_artifact(&c, &preview["artifact_id"]).await;
    assert_eq!(full["workspace"], input);
    assert_eq!(full["source"], json!(pin));
    let native_preview: native::pcr::MismatchPreview =
        serde_json::from_value(full["preview"].clone()).unwrap();
    native_preview
        .validate_binding(preview["context_revision"].as_str().unwrap())
        .unwrap();
    assert_eq!(
        native_preview.preview.documents[0].sequence.as_str(),
        "CCAGTTCGATGAATTCCGTAGGACTTCAGACCTGTAGATCGTCACC"
    );
    let request = InstallRequest {
        sources: vec![collagene_mcp::tools::workspace::MaterialDerivationSource {
            role: collagene_core::material::derivation::SourceRole::Template,
            event: pin.clone(),
            origin: Origin::OriginalReference,
        }],
        derivation: Derivation::MismatchPcr {
            forward: forward.clone(),
            reverse: reverse.clone(),
            options: options.clone(),
            selection: MaterialPcrChoice {
                selected_index: 0,
                expected_preview_revision: preview["previews"][0]["preview_revision"]
                    .as_str()
                    .unwrap()
                    .into(),
            },
        },
        name: "Reusable tailed substitution".into(),
        save_primers: true,
    };
    let mut changed = request.clone();
    let Derivation::MismatchPcr { options, .. } = &mut changed.derivation else {
        unreachable!()
    };
    options.binding.min_three_prime_exact = 5;
    reject(
        &c,
        json!({"workspace":input,"task":{"Install":{"request":changed}}}),
    )
    .await;
    let mut injected = serde_json::to_value(&request).unwrap();
    injected["derivation"]["MismatchPcr"]["selection"]["geometry"] = json!({"sequence":STOCK});
    reject(
        &c,
        json!({"workspace":input,"task":{"Install":{"request":injected}}}),
    )
    .await;
    let mut wrong = args;
    wrong["task"]["PreviewMismatchPcr"]["source"]["event_id"] = json!(event(&original, 1).event_id);
    reject(&c, wrong).await;
    let installed = call(
        &c,
        "material_derivation",
        json!({"workspace":input,"task":{"Install":{"request":request}}}),
    )
    .await;
    let receipt: native::Receipt = serde_json::from_value(installed["receipt"].clone()).unwrap();
    assert!(receipt.primer_attachment.is_some());
    let w = native_archive(&c, &installed["native_artifact_id"]).await;
    receipt.validate_archive(&w).unwrap();
    assert!(
        collagene_core::material::lineage::list_runs(&w)
            .unwrap()
            .is_empty()
    );
    let capsule = native::capsule(&w, &receipt.product).unwrap();
    let record = native::load_record(&capsule.record_bytes).unwrap();
    let orders = native::primer_table(&record).unwrap().unwrap();
    assert_eq!(
        orders.records[0]
            .primer
            .as_ref()
            .unwrap()
            .ordered_sequence(),
        "CCAGTTCGAT"
    );
    assert_eq!(
        orders.records[1]
            .primer
            .as_ref()
            .unwrap()
            .ordered_sequence(),
        "GGTGACGATC"
    );
    let native::Derivation::MismatchPcr {
        options, selection, ..
    } = record.request.derivation
    else {
        panic!("full mismatch recipe")
    };
    assert_eq!(options.binding.min_three_prime_exact, 4);
    assert_eq!(selection.geometry.substitutions.len(), 1);
    let catalog = catalog_digest::builtin_catalog().unwrap();
    let eco = catalog
        .edition()
        .enzymes
        .iter()
        .find(|e| e.name == "EcoRI")
        .unwrap();
    let request = InstallRequest {
        sources: vec![collagene_mcp::tools::workspace::MaterialDerivationSource {
            role: collagene_core::material::derivation::SourceRole::Template,
            event: receipt.product.clone(),
            origin: Origin::Saved(run_material::MaterialSelection::Derived(receipt.product)),
        }],
        derivation: Derivation::CatalogDigest {
            conditions: None,
            selection: catalog_digest::select(&catalog, &[eco.id.clone()]).unwrap(),
            options: Default::default(),
            fragment_index: 1,
        },
        name: "Cut mismatch-PCR product".into(),
        save_primers: false,
    };
    let cut = call(
        &c,
        "material_derivation",
        json!({"workspace":reference(&installed),"task":{"Install":{"request":request}}}),
    )
    .await;
    let receipt: native::Receipt = serde_json::from_value(cut["receipt"].clone()).unwrap();
    let origin = Origin::Saved(run_material::MaterialSelection::Derived(
        receipt.product.clone(),
    ));
    let preview=call(&c,"material_derivation",json!({"workspace":reference(&cut),"task":{"PreviewRepair":{"source":receipt.product,"origin":origin,"method":EndRepairMethod::FillFivePrimeOverhangs}}})).await;
    assert_eq!(preview["repair"]["source"]["bottom_start"], 4);
    assert_eq!(preview["repair"]["duplex"]["bottom_start"], 0);
    assert_eq!(preview["repair"]["bottom_added_3prime"], "AATT");
    let full = json_artifact(&c, &preview["artifact_id"]).await;
    assert_eq!(full["source"], json!(receipt.product));
    let full: native::end_repair::RepairPreview =
        serde_json::from_value(full["preview"].clone()).unwrap();
    full.validate_binding(preview["context_revision"].as_str().unwrap())
        .unwrap();
    reject(&c,json!({"workspace":reference(&cut),"task":{"PreviewRepair":{"source":receipt.product,"origin":"OriginalReference","method":"Blunt"}}})).await;
    let request = InstallRequest {
        sources: vec![collagene_mcp::tools::workspace::MaterialDerivationSource {
            role: collagene_core::material::derivation::SourceRole::Template,
            event: receipt.product,
            origin: origin,
        }],
        derivation: Derivation::EndRepair {
            method: EndRepairMethod::FillFivePrimeOverhangs,
        },
        name: "Repaired reusable source".into(),
        save_primers: false,
    };
    let mut wrong = request.clone();
    wrong.save_primers = true;
    reject(
        &c,
        json!({"workspace":reference(&cut),"task":{"Install":{"request":wrong}}}),
    )
    .await;
    let repaired = call(
        &c,
        "material_derivation",
        json!({"workspace":reference(&cut),"task":{"Install":{"request":request}}}),
    )
    .await;
    let receipt: native::Receipt = serde_json::from_value(repaired["receipt"].clone()).unwrap();
    let w = native_archive(&c, &repaired["native_artifact_id"]).await;
    receipt.validate_archive(&w).unwrap();
    assert!(receipt.primer_attachment.is_none());
    let fresh = native::reconstruct(&w, &receipt.product).unwrap();
    assert_eq!(
        fresh.reference().sequence.as_str(),
        "AATTCCGTAGGACTTCAGACCTGTAGATCGTCACC"
    );
    assert_eq!(fresh.tracked().duplex.bottom_start, 0);
    let kept = call(&c, "export_workspace_archive", input).await;
    assert_eq!(
        workspace::save_workspace(&native_archive(&c, &kept["artifact_id"]).await).unwrap(),
        bytes
    );
    c.cancel().await.unwrap();
    server.await.unwrap();
}
