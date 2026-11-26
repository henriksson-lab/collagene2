//! Real MCP SDK routing with literal checksums, packet framing and source pins.
use base64::{Engine as _, engine::general_purpose::STANDARD};
use collagene_mcp::CollageneServer;
use rmcp::{RoleClient, ServiceExt, model::CallToolRequestParams, service::RunningService};
use serde_json::{Value, json};
type Client = RunningService<RoleClient, ()>;
async fn call(c: &Client, name: &str, args: Value) -> Value {
    let r = c
        .call_tool(
            CallToolRequestParams::new(name.to_string())
                .with_arguments(args.as_object().unwrap().clone()),
        )
        .await
        .unwrap();
    assert_ne!(r.is_error, Some(true), "{name}: {:?}", r.content);
    r.structured_content.unwrap()
}
async fn reject(c: &Client, name: &str, args: Value) {
    let r = c
        .call_tool(
            CallToolRequestParams::new(name.to_string())
                .with_arguments(args.as_object().unwrap().clone()),
        )
        .await;
    assert!(r.is_err() || r.unwrap().is_error == Some(true));
}
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
fn document(result: &Value) -> Value {
    let d = &result["documents"][0];
    json!({"id":d["id"],"revision":d["revision"]})
}
async fn bytes(c: &Client, result: &Value) -> Vec<u8> {
    let chunk = call(
        c,
        "read_artifact",
        json!({"artifact_id":result["artifact_id"],"start":0,"max_bytes":16384}),
    )
    .await;
    assert_eq!(chunk["done"], true);
    STANDARD.decode(chunk["text"].as_str().unwrap()).unwrap()
}
#[tokio::test]
async fn gcg_and_lower_xml_use_explicit_profiles_and_preserve_session_sources() {
    let (c, server) = setup().await;
    let tools = c.list_all_tools().await.unwrap();
    assert_eq!(tools.len(), collagene_mcp::CURRENT_TOOL_COUNT);
    for name in ["import_compat_sequence", "export_gcg", "export_snapgene"] {
        assert_eq!(
            tools.iter().find(|t| t.name == name).unwrap().input_schema["additionalProperties"],
            false
        );
    }
    let text = "!!NA_SEQUENCE 1.0\n\nsample  Length: 11  November 26, 2026  Type: N  Check: 4920  ..\n\n 1 GATTCCGATTG\n";
    let imported = call(
        &c,
        "import_compat_sequence",
        json!({"text":text,"profile":{"Gcg":{"alphabet":"Dna","fallback_name":"literal"}}}),
    )
    .await;
    let d = document(&imported);
    let before = call(&c, "export_document", d.clone()).await;
    let review = call(
        &c,
        "export_gcg",
        json!({"document":d,"options":{"record_id":"new_record","date_label":"2026-11-26"}}),
    )
    .await;
    assert_eq!(review["state"], "review_required");
    assert!(review["artifact_id"].is_null());
    let export=call(&c,"export_gcg",json!({"document":d,"options":{"record_id":"new_record","date_label":"2026-11-26"},"accept_omissions":true})).await;
    let out = String::from_utf8(bytes(&c, &export).await).unwrap();
    assert!(out.contains("Check: 4920"));
    assert!(out.contains("Type: N"));
    let body = out.split_once("..\n").unwrap().1;
    let calls: String = body
        .split_whitespace()
        .filter(|word| word.bytes().all(|b| b.is_ascii_alphabetic()))
        .collect();
    assert_eq!(calls, "GATTCCGATTG");
    assert_eq!(call(&c, "export_document", d.clone()).await, before);
    reject(
        &c,
        "import_compat_sequence",
        json!({"text":text,"profile":{"Gcg":{"alphabet":"Protein","fallback_name":"bad"}}}),
    )
    .await;
    reject(&c,"import_compat_sequence",json!({"text":text.replace("4920","4921"),"profile":{"Gcg":{"alphabet":"Dna","fallback_name":"bad"}}})).await;
    let xml = "<f><sequence name=\"lower\" su=\"    \" sl=\"ACGA\" circular=\"false\" notes=\"literal lower 3 to 5\"/></f>";
    reject(
        &c,
        "import_compat_sequence",
        json!({"text":xml,"profile":{"CollageneXml":{"policy":"UpperOnly"}}}),
    )
    .await;
    let lower = call(
        &c,
        "import_compat_sequence",
        json!({"text":xml,"profile":{"CollageneXml":{"policy":"AllowLowerOnly"}}}),
    )
    .await;
    let native = call(&c, "export_document", document(&lower)).await;
    let m = collagene_core::protocol::load_document(&bytes(&c, &native).await).unwrap();
    assert_eq!(m.sequence.as_str(), "AGCA");
    assert_eq!(m.metadata["source:legacy:chemistry_unmapped"], "true");
    assert_eq!(
        collagene_core::legacy_formats::original_legacy_bytes(&m).unwrap(),
        xml.as_bytes()
    );
    reject(
        &c,
        "export_snapgene",
        json!({"document":document(&lower),"accept_omissions":true}),
    )
    .await;
    c.cancel().await.unwrap();
    server.await.unwrap();
}
#[tokio::test]
async fn snapgene_export_requires_loss_review_and_emits_fresh_reference_packets() {
    let (c, server) = setup().await;
    let molecule =
        collagene_core::Molecule::dna("circle", "ACGTTGCA", collagene_core::Topology::Circular)
            .unwrap();
    let gb = collagene_core::formats::export_genbank(&molecule).unwrap();
    let imported = call(&c, "import_document", json!({"text":gb})).await;
    let d = document(&imported);
    let before = call(&c, "export_document", d.clone()).await;
    let review = call(&c, "export_snapgene", json!({"document":d})).await;
    assert_eq!(review["state"], "review_required");
    assert!(review["artifact_id"].is_null());
    assert!(!review["review"]["losses"].as_array().unwrap().is_empty());
    let result = call(
        &c,
        "export_snapgene",
        json!({"document":d,"accept_omissions":true}),
    )
    .await;
    assert_eq!(result["filename"], "sequence.dna");
    let raw = bytes(&c, &result).await;
    assert_eq!(raw[0], 9);
    assert_eq!(u32::from_be_bytes(raw[1..5].try_into().unwrap()), 14);
    assert_eq!(&raw[5..13], b"SnapGene");
    let read = collagene_core::snapgene::import_snapgene(&raw).unwrap();
    assert_eq!(read.molecules[0].sequence.as_str(), "ACGTTGCA");
    assert_eq!(
        read.molecules[0].topology,
        collagene_core::Topology::Circular
    );
    assert_eq!(call(&c, "export_document", d).await, before);
    c.cancel().await.unwrap();
    server.await.unwrap();
}
