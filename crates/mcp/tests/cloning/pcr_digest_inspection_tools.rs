//! Inspect actual PCR-added catalog sites, with no target/product publication.
use collagene_core::{
    chemistry::FivePrime,
    cloning::{pcr::digest::Preview, restriction::inspection::InspectionOptions},
    sequence::analysis::Primer,
};
use collagene_mcp::server::CollageneServer;
use rmcp::{RoleClient, ServiceExt, model::CallToolRequestParams, service::RunningService};
use serde_json::{Value, json};
type Client = RunningService<RoleClient, ()>;
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
    let result = c
        .call_tool(
            CallToolRequestParams::new(name.to_owned())
                .with_arguments(args.as_object().unwrap().clone()),
        )
        .await
        .unwrap();
    assert_ne!(result.is_error, Some(true), "{name}: {:?}", result.content);
    result.structured_content.unwrap()
}
async fn reject(c: &Client, name: &str, args: Value) {
    let result = c
        .call_tool(
            CallToolRequestParams::new(name.to_owned())
                .with_arguments(args.as_object().unwrap().clone()),
        )
        .await;
    assert!(result.is_err() || result.unwrap().is_error == Some(true));
}
async fn evidence(c: &Client, id: &Value) -> Value {
    let mut text = String::new();
    let mut start = 0;
    loop {
        let page = call(
            c,
            "read_artifact",
            json!({"artifact_id":id,"start":start,"max_bytes":16384}),
        )
        .await;
        text.push_str(page["text"].as_str().unwrap());
        if page["done"] == true {
            break;
        }
        start = page["next_start"].as_u64().unwrap();
    }
    serde_json::from_str(&text).unwrap()
}
fn named(imported: &Value, name: &str) -> Value {
    let row = imported["documents"]
        .as_array()
        .unwrap()
        .iter()
        .find(|d| d["name"] == name)
        .unwrap();
    json!({"id":row["id"],"revision":row["revision"]})
}
#[tokio::test]
async fn actual_added_sites_unique_versus_ambiguous_and_stale_originals() {
    let (client, server) = setup().await;
    let tools = client.list_tools(Default::default()).await.unwrap().tools;
    assert_eq!(tools.len(), collagene_mcp::CURRENT_TOOL_COUNT);
    assert!(tools.iter().any(|t| t.name == "inspect_pcr_digest"));
    let imported = call(
        &client,
        "import_document",
        json!({"text":">unique\nAAACCCGGGTTT\n>ambiguous\nAAACCCGGGTTTAAACCCGGGTTT\n"}),
    )
    .await;
    let catalog = collagene_core::cloning::restriction::catalog::builtin_catalog().unwrap();
    let id = catalog
        .edition()
        .enzymes
        .iter()
        .find(|e| e.name == "EcoRI")
        .unwrap()
        .id
        .clone();
    let selection = collagene_core::cloning::restriction::catalog::select(&catalog, &[id]).unwrap();
    let request = json!({"template_strand":"Top","source_role":"original","forward":Primer::new("Added forward site","AAA","ATGCGAATTC").unwrap(),"reverse":Primer::new("Added reverse site","AAA","ATGCGAATTC").unwrap(),"selection":selection,"options":InspectionOptions::default()});
    let args = json!({"document":named(&imported,"unique"),"request":request,"fragment_start":0,"max_fragments":2});
    let result = call(&client, "inspect_pcr_digest", args.clone()).await;
    assert_eq!(result["pcr_geometry_count"], 1);
    assert_eq!(result["amplicon_length"], 32);
    assert_eq!(result["products_determined"], true);
    assert_eq!(result["fragment_count"], 3);
    assert_eq!(result["fragments"].as_array().unwrap().len(), 2);
    assert_eq!(result["fragments_done"], false);
    let full = evidence(&client, &result["artifact_id"]).await;
    assert_eq!(
        full["amplicon"]["sequence"]["bases"],
        "ATGCGAATTCAAACCCGGGTTTGAATTCGCAT"
    );
    let decoded: Preview = serde_json::from_value(full.clone()).unwrap();
    decoded
        .validate_binding(result["context_revision"].as_str().unwrap())
        .unwrap();
    let fragments = &decoded.inspection.as_ref().unwrap().fragments;
    assert_eq!(
        fragments[0]
            .duplex
            .ends
            .as_ref()
            .unwrap()
            .0
            .chemistry
            .five_prime,
        FivePrime::Hydroxyl
    );
    assert_eq!(
        fragments[1]
            .duplex
            .ends
            .as_ref()
            .unwrap()
            .0
            .chemistry
            .five_prime,
        FivePrime::Phosphate
    );
    let mut next = args.clone();
    next["fragment_start"] = result["next_fragment_start"].clone();
    let page = call(&client, "inspect_pcr_digest", next).await;
    assert_eq!(page["fragments"].as_array().unwrap().len(), 1);
    assert_eq!(page["artifact_id"], result["artifact_id"]);
    let mut ambiguous = args.clone();
    ambiguous["document"] = named(&imported, "ambiguous");
    let ambiguous = call(&client, "inspect_pcr_digest", ambiguous).await;
    assert!(ambiguous["pcr_geometry_count"].as_u64().unwrap() > 1);
    assert_eq!(ambiguous["products_determined"], false);
    assert_eq!(ambiguous["fragment_count"], 0);
    assert_eq!(ambiguous["amplicon_length"], Value::Null);
    let mut stale = args.clone();
    stale["document"]["revision"] = json!("stale");
    reject(&client, "inspect_pcr_digest", stale).await;
    let mut bad_pin = args.clone();
    bad_pin["request"]["selection"]["catalog_content_hash"] = json!("0".repeat(64));
    reject(&client, "inspect_pcr_digest", bad_pin).await;
    let mut oversized = args;
    oversized["max_fragments"] = json!(101);
    reject(&client, "inspect_pcr_digest", oversized).await;
    assert_eq!(call(&client, "list_documents", json!({})).await["total"], 2);
    drop(client);
    server.abort();
}
