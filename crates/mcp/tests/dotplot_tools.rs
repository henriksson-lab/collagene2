//! Actual SDK routing: independently specified coordinates and immutable artifacts.
use collagene_core::dotplot::*;
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
    assert_ne!(r.is_error, Some(true), "{name}:{:?}", r.content);
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
async fn text(c: &Client, id: &Value) -> String {
    let mut result = String::new();
    let mut at = 0;
    loop {
        let page = call(
            c,
            "read_artifact",
            json!({"artifact_id":id,"start":at,"max_bytes":16384}),
        )
        .await;
        result.push_str(page["text"].as_str().unwrap());
        if page["done"] == true {
            break;
        }
        at = page["next_start"].as_u64().unwrap();
    }
    result
}
#[tokio::test]
async fn dotplot_sdk_retains_bound_full_evidence_svg_and_exact_paginated_hits() {
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
    let c = ().serve(b).await.unwrap();
    let imported = call(
        &c,
        "import_document",
        json!({"text":">reference\nACGTTT\n>query\nACG\n"}),
    )
    .await;
    let doc = |name: &str| {
        let d = imported["documents"]
            .as_array()
            .unwrap()
            .iter()
            .find(|d| d["name"] == name)
            .unwrap();
        json!({"id":d["id"],"revision":d["revision"]})
    };
    let r = doc("reference");
    let q = doc("query");
    let options = DotPlotOptions {
        window: 3,
        min_matches: 3,
        ..Default::default()
    };
    let first = call(
        &c,
        "analyze_dotplot",
        json!({"reference":r,"query":q,"options":options,"max_hits":1}),
    )
    .await;
    assert_eq!(first["reference_revision"], r["revision"]);
    assert_eq!(first["query_revision"], q["revision"]);
    assert_eq!(first["complete"], true);
    assert_eq!(first["hit_count"], 2);
    assert_eq!(first["hits"].as_array().unwrap().len(), 1);
    assert_eq!(first["done"], false);
    let full: DotPlotReport = serde_json::from_str(&text(&c, &first["artifact_id"]).await).unwrap();
    assert_eq!(full.hits.len(), 2);
    assert!(
        full.hits
            .iter()
            .any(|h| h.reference_start == 1 && h.query_start == 2 && h.reverse)
    );
    let svg = text(&c, &first["svg_artifact_id"]).await;
    assert!(svg.contains("<svg") && svg.contains("dot"));
    let next = call(
        &c,
        "analyze_dotplot",
        json!({"reference":r,"query":q,"options":options,"hit_start":1,"max_hits":1}),
    )
    .await;
    assert_eq!(first["artifact_id"], next["artifact_id"]);
    assert_eq!(first["svg_artifact_id"], next["svg_artifact_id"]);
    assert_eq!(next["done"], true);
    assert_ne!(first["hits"][0]["id"], next["hits"][0]["id"]);
    let mut sampled = options.clone();
    sampled.stride = 2;
    let report = call(
        &c,
        "analyze_dotplot",
        json!({"reference":r,"query":q,"options":sampled}),
    )
    .await;
    assert_eq!(report["sampled"], true);
    assert_eq!(report["complete"], false);
    let mut limited = options.clone();
    limited.max_work = 1;
    let report = call(
        &c,
        "analyze_dotplot",
        json!({"reference":r,"query":q,"options":limited}),
    )
    .await;
    assert_eq!(report["complete"], false);
    let mut stale = r.clone();
    stale["revision"] = json!("0".repeat(64));
    reject(
        &c,
        "analyze_dotplot",
        json!({"reference":stale,"query":q,"options":options}),
    )
    .await;
    let mut unknown = serde_json::to_value(&options).unwrap();
    unknown["invented_assumption"] = json!(true);
    reject(
        &c,
        "analyze_dotplot",
        json!({"reference":r,"query":q,"options":unknown}),
    )
    .await;
    reject(
        &c,
        "analyze_dotplot",
        json!({"reference":r,"query":q,"options":options,"max_hits":201}),
    )
    .await;
    let documents = call(&c, "list_documents", json!({})).await;
    assert_eq!(documents["documents"].as_array().unwrap().len(), 2);
    c.cancel().await.unwrap();
    server.await.unwrap();
}
#[tokio::test]
async fn dotplot_json_and_svg_installation_is_atomic_at_artifact_capacity() {
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
    let c = ().serve(b).await.unwrap();
    let imported = call(
        &c,
        "import_document",
        json!({"text":">reference\nACGTTT\n>query\nACG\n>protein [molecule=protein]\nACW\n"}),
    )
    .await;
    let rows = imported["documents"].as_array().unwrap();
    let doc = |name: &str| {
        let row = rows
            .iter()
            .find(|d| {
                d["name"]
                    .as_str()
                    .is_some_and(|n| n.split_whitespace().next() == Some(name))
            })
            .unwrap();
        json!({"id":row["id"],"revision":row["revision"]})
    };
    for n in 0..127 {
        let options = collagene_core::protein_properties::ProteinOptions {
            charge_ph: n as f64 / 10.0,
            ..Default::default()
        };
        call(
            &c,
            "analyze_protein",
            json!({"document":doc("protein"),"options":options}),
        )
        .await;
    }
    let options = DotPlotOptions {
        window: 3,
        min_matches: 3,
        ..Default::default()
    };
    reject(
        &c,
        "analyze_dotplot",
        json!({"reference":doc("reference"),"query":doc("query"),"options":options}),
    )
    .await;
    // If either dot-plot artifact had leaked before the second failed, this
    // unrelated new report could no longer occupy the final available slot.
    let options = collagene_core::protein_properties::ProteinOptions {
        charge_ph: 13.0,
        ..Default::default()
    };
    call(
        &c,
        "analyze_protein",
        json!({"document":doc("protein"),"options":options}),
    )
    .await;
    let options = collagene_core::protein_properties::ProteinOptions {
        charge_ph: 13.1,
        ..Default::default()
    };
    reject(
        &c,
        "analyze_protein",
        json!({"document":doc("protein"),"options":options}),
    )
    .await;
    c.cancel().await.unwrap();
    server.await.unwrap();
}
