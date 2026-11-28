//! Actual SDK over memory pipes; no browser, HTTP, user files or scientific process.
use collagene_mcp::server::CollageneServer;
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
    match c
        .call_tool(
            CallToolRequestParams::new(name.to_owned())
                .with_arguments(args.as_object().unwrap().clone()),
        )
        .await
    {
        Ok(r) => assert_eq!(r.is_error, Some(true)),
        Err(rmcp::service::ServiceError::McpError(e)) => assert_eq!(e.code.0, -32602),
        Err(e) => panic!("unexpected transport failure: {e}"),
    }
}

#[tokio::test]
async fn typed_discovery_source_pins_full_artifacts_and_counts_only_requests_remain_readonly() {
    use base64::{Engine as _, engine::general_purpose::STANDARD};
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
    let tools = c.list_all_tools().await.unwrap();
    assert_eq!(tools.len(), collagene_mcp::CURRENT_TOOL_COUNT);
    for name in [
        "list_protein_liability_profiles",
        "screen_protein_liabilities",
    ] {
        assert!(tools.iter().any(|t| t.name == name));
    }
    let schema = serde_json::to_string(
        &tools
            .iter()
            .find(|t| t.name == "screen_protein_liabilities")
            .unwrap()
            .input_schema,
    )
    .unwrap();
    assert!(schema.contains("provided_cdrs") && schema.contains("source_revision"));
    let profile = call(&c, "list_protein_liability_profiles", json!({})).await;
    assert_eq!(profile["rules"].as_array().unwrap().len(), 6);
    let imported = call(
        &c,
        "import_compat_sequence",
        json!({"text":">source\nNATMWX\n","profile":{"Fasta":{"alphabet":"Protein"}}}),
    )
    .await;
    let d = &imported["documents"][0];
    let doc = json!({"id":d["id"],"revision":d["revision"]});
    let before = call(&c, "list_documents", json!({})).await["total"].clone();
    let options = json!({"max_hits":1,"provided_cdrs":{"source_revision":doc["revision"],"definition_label":"Supplied coordinate fixture","provenance":"Manually authored context, no numbering claim","regions":[{"domain_id":"variable-a","kind":"Cdr1","start":0,"end":1}]}});
    let result = call(
        &c,
        "screen_protein_liabilities",
        json!({"document":doc,"options":options}),
    )
    .await;
    assert_eq!(result["retained_hit_count"], 1);
    assert_eq!(result["hits_truncated"], true);
    assert_eq!(result["cdr_rule_unassessed"], false);
    let full = call(
        &c,
        "read_artifact",
        json!({"artifact_id":result["artifact_id"],"start":0,"max_bytes":16384}),
    )
    .await;
    assert_eq!(full["done"], true);
    let report: collagene_core::protein_liabilities::Report =
        serde_json::from_str(full["text"].as_str().unwrap()).unwrap();
    report
        .validate_binding(result["context_revision"].as_str().unwrap())
        .unwrap();
    let exported = call(&c, "export_document", doc.clone()).await;
    let binary = call(
        &c,
        "read_artifact",
        json!({"artifact_id":exported["artifact_id"],"start":0,"max_bytes":16384}),
    )
    .await;
    assert_eq!(binary["done"], true);
    let m = collagene_core::protocol::load_document(
        &STANDARD.decode(binary["text"].as_str().unwrap()).unwrap(),
    )
    .unwrap();
    let native_options = serde_json::from_value(options.clone()).unwrap();
    report.validate_for(&m, &native_options).unwrap();
    let mut no_hits = options.clone();
    no_hits["max_hits"] = 0.into();
    let counted = call(
        &c,
        "screen_protein_liabilities",
        json!({"document":doc,"options":no_hits}),
    )
    .await;
    assert_eq!(counted["retained_hit_count"], 0);
    assert_eq!(counted["counts"], result["counts"]);
    assert_eq!(counted["hits_truncated"], true);
    let mut stale = doc.clone();
    stale["revision"] = json!("f".repeat(64));
    reject(
        &c,
        "screen_protein_liabilities",
        json!({"document":stale,"options":options}),
    )
    .await;
    let mut foreign = options.clone();
    foreign["provided_cdrs"]["source_revision"] = json!("f".repeat(64));
    reject(
        &c,
        "screen_protein_liabilities",
        json!({"document":doc,"options":foreign}),
    )
    .await;
    let artifact_doc = json!({"id":result["artifact_id"],"revision":result["artifact_id"]});
    reject(
        &c,
        "screen_protein_liabilities",
        json!({"document":artifact_doc,"options":options}),
    )
    .await;
    assert_eq!(call(&c, "list_documents", json!({})).await["total"], before);
    c.cancel().await.unwrap();
    server.await.unwrap();
}
