//! Actual SDK dispatch for bounded source-bound inspection, without imported report authority.
use collagene_core::{formats, restriction_inspection::*};
use collagene_mcp::CollageneServer;
use rmcp::{RoleClient, ServiceExt, model::CallToolRequestParams, service::RunningService};
use serde_json::{Value, json};
type Client = RunningService<RoleClient, ()>;
async fn call(c: &Client, name: &str, args: Value) -> Value {
    let r = c
        .call_tool(
            CallToolRequestParams::new(name.to_owned())
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
            CallToolRequestParams::new(name.to_owned())
                .with_arguments(args.as_object().unwrap().clone()),
        )
        .await;
    assert!(r.is_err() || r.unwrap().is_error == Some(true));
}
async fn setup() -> (Client, tokio::task::JoinHandle<()>) {
    let (a, b) = tokio::io::duplex(65536);
    let s = tokio::spawn(async move {
        CollageneServer::default()
            .serve(a)
            .await
            .unwrap()
            .waiting()
            .await
            .unwrap();
    });
    (().serve(b).await.unwrap(), s)
}
async fn import(c: &Client, dna: &str) -> Value {
    let r = call(
        c,
        "import_document",
        json!({"text":format!(">source\n{dna}\n")}),
    )
    .await;
    json!({"id":r["documents"][0]["id"],"revision":r["documents"][0]["revision"]})
}
async fn artifact(c: &Client, id: &Value) -> Value {
    let mut text = String::new();
    let mut start = 0;
    loop {
        let p = call(
            c,
            "read_artifact",
            json!({"artifact_id":id,"start":start,"max_bytes":1024}),
        )
        .await;
        text.push_str(p["text"].as_str().unwrap());
        if p["done"] == true {
            break;
        }
        start = p["next_start"].as_u64().unwrap();
    }
    serde_json::from_str(&text).unwrap()
}
#[tokio::test]
async fn inspection_sdk_pages_exact_geometry_and_source_without_creating_products() {
    let (c, s) = setup().await;
    let tools = c.list_tools(Default::default()).await.unwrap().tools;
    assert_eq!(tools.len(), 93);
    let tool = tools
        .iter()
        .find(|t| t.name == "inspect_restriction")
        .unwrap();
    assert_eq!(tool.input_schema["additionalProperties"], false);
    let dna = "AAAAGAATTCCCCCGGATCCAAAA";
    let d = import(&c, dna).await;
    let before = call(&c, "list_documents", json!({})).await;
    let args = json!({"document":d,"source_role":"source-event","enzymes":["ecori","BamHI"],"max_sites":1,"max_fragments":1});
    let first = call(&c, "inspect_restriction", args.clone()).await;
    assert_eq!(first["scan_complete"], true);
    assert_eq!(first["products_determined"], true);
    assert_eq!(first["site_count"], 2);
    assert_eq!(first["fragment_count"], 3);
    assert_eq!(first["sites_done"], false);
    assert_eq!(first["fragments_done"], false);
    assert_eq!(first["sites"][0]["recognition_length"], 6);
    assert_eq!(first["fragments"][0]["duplex"]["top_length"], 5);
    let report: RestrictionInspection =
        serde_json::from_value(artifact(&c, &first["artifact_id"]).await).unwrap();
    let molecule = formats::import_text(&format!(">source\n{dna}\n"))
        .unwrap()
        .molecules
        .remove(0);
    let enzymes = collagene_core::analysis::starter_enzymes()
        .into_iter()
        .filter(|e| matches!(e.name.as_str(), "EcoRI" | "BamHI"))
        .collect::<Vec<_>>();
    report
        .validate_bindings_for(
            "source-event",
            &molecule,
            &enzymes,
            &InspectionOptions::default(),
        )
        .unwrap();
    assert_eq!(
        report
            .fragments
            .iter()
            .map(|f| f.duplex.top_length)
            .collect::<Vec<_>>(),
        vec![5, 10, 9]
    );
    let mut next = args.clone();
    next["site_start"] = json!(1);
    next["fragment_start"] = json!(1);
    let next = call(&c, "inspect_restriction", next).await;
    assert_eq!(next["artifact_id"], first["artifact_id"]);
    assert_eq!(next["sites"][0]["recognition_start"], 14);
    assert_eq!(next["fragments"][0]["top_start"], 5);
    let mut other = args.clone();
    other["source_role"] = json!("other-event");
    let other = call(&c, "inspect_restriction", other).await;
    assert_ne!(other["batch_revision"], first["batch_revision"]);
    for (key, val) in [
        ("site_start", json!(999)),
        ("fragment_start", json!(999)),
        ("max_sites", json!(0)),
        ("max_fragments", json!(101)),
        ("enzymes", json!(["unknown"])),
        ("trusted_report", json!(true)),
    ] {
        let mut invalid = args.clone();
        invalid[key] = val;
        reject(&c, "inspect_restriction", invalid).await;
    }
    let mut stale = args;
    stale["document"]["revision"] = json!("0".repeat(64));
    reject(&c, "inspect_restriction", stale).await;
    assert_eq!(call(&c, "list_documents", json!({})).await, before);
    c.cancel().await.unwrap();
    s.await.unwrap();
}
#[tokio::test]
async fn possible_edge_and_incomplete_evidence_survive_while_legacy_tool_stays_strict() {
    let (c, s) = setup().await;
    for (dna, name) in [("AAAAGAANTCAAAA", "EcoRI"), ("GGTCTCA", "BsaI")] {
        let d = import(&c, dna).await;
        let r = call(
            &c,
            "inspect_restriction",
            json!({"document":d,"source_role":"source","enzymes":[name]}),
        )
        .await;
        assert_eq!(r["scan_complete"], true);
        assert_eq!(r["products_determined"], false);
        assert_eq!(r["fragment_count"], 0);
        assert!(!r["sites"][0]["diagnostics"].as_array().unwrap().is_empty());
        let report: RestrictionInspection =
            serde_json::from_value(artifact(&c, &r["artifact_id"]).await).unwrap();
        assert!(!report.products_determined);
        reject(
            &c,
            "analyze_restriction",
            json!({"document":d,"enzyme":name}),
        )
        .await;
    }
    let d = import(&c, "AAAAGAATTCCCCCGGATCCAAAA").await;
    let mut o = InspectionOptions::default();
    o.max_comparisons = 1;
    let r = call(
        &c,
        "inspect_restriction",
        json!({"document":d,"source_role":"source","enzymes":["EcoRI"],"options":o}),
    )
    .await;
    assert_eq!(r["scan_complete"], false);
    assert_eq!(r["products_determined"], false);
    assert_eq!(r["limits"][0], "Comparisons");
    c.cancel().await.unwrap();
    s.await.unwrap();
}
