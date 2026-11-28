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
async fn rna_profiles_annotation_review_and_native_reanalysis_bind_immutable_documents() {
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
        "list_rna_chemistry_profiles",
        "inspect_rna_annotations",
        "analyze_rna_chemistry",
    ] {
        assert!(tools.iter().any(|t| t.name == name));
    }
    let schema = serde_json::to_string(
        &tools
            .iter()
            .find(|t| t.name == "analyze_rna_chemistry")
            .unwrap()
            .input_schema,
    )
    .unwrap();
    assert!(
        schema.contains("annotation_reviews")
            && schema.contains("coverage")
            && schema.contains("Phosphorothioate")
    );
    let profiles = call(&c, "list_rna_chemistry_profiles", json!({})).await;
    assert_eq!(profiles["declaration_version"], 2);
    assert_eq!(
        profiles["annotation_profiles"].as_array().unwrap().len(),
        10
    );
    let gb = "LOCUS       reviewed_rna               2 bp    mRNA    linear   SYN 27-NOV-2026\nFEATURES             Location/Qualifiers\n     modified_base   1\n                     /mod_base=\"m6a\"\nORIGIN\n        1 ac\n//\n";
    let imported = call(&c, "import_document", json!({"text":gb})).await;
    let doc = json!({"id":imported["documents"][0]["id"],"revision":imported["documents"][0]["revision"]});
    let before = call(&c, "list_documents", json!({})).await;
    let found = call(&c, "inspect_rna_annotations", json!({"document":doc})).await;
    assert_eq!(found["count"], 1);
    let candidate = &found["candidates"][0];
    assert_eq!(candidate["position"], 0);
    let mut d = json!({"version":2,"source_revision":doc["revision"],"coverage":"ClosedModel","ends":{"Linear":{"five":"Hydroxyl","three":"Hydroxyl"}},"sites":[candidate["suggested_site"]],"links":[],"annotation_reviews":[],"provenance":{"statement":"Explicit covalent model fixture","reference":null}});
    let unresolved = call(
        &c,
        "analyze_rna_chemistry",
        json!({"document":doc,"declaration":d}),
    )
    .await;
    assert_eq!(unresolved["formula_complete"], false);
    assert!(unresolved["mass"].is_null());
    d["annotation_reviews"] = json!([{"feature_index":0,"expected_feature_revision":candidate["feature_revision"],"site_entry":0}]);
    let analyzed = call(
        &c,
        "analyze_rna_chemistry",
        json!({"document":doc,"declaration":d}),
    )
    .await;
    assert_eq!(analyzed["formula_complete"], true);
    assert_eq!(analyzed["formula"]["c"], 20);
    assert_eq!(analyzed["unreviewed_evidence_count"], 0);
    assert_eq!(analyzed["retained_evidence_count"], 1);
    let artifact = call(
        &c,
        "read_artifact",
        json!({"artifact_id":analyzed["artifact_id"],"start":0,"max_bytes":16384}),
    )
    .await;
    assert_eq!(artifact["done"], true);
    let report: collagene_core::rna_chemistry::Report =
        serde_json::from_str(artifact["text"].as_str().unwrap()).unwrap();
    let source = collagene_core::formats::sequence::import_text(gb)
        .unwrap()
        .molecules
        .remove(0);
    let declaration = serde_json::from_value(d.clone()).unwrap();
    report.validate_for(&source, &declaration).unwrap();
    let repeated = call(
        &c,
        "analyze_rna_chemistry",
        json!({"document":doc,"declaration":d}),
    )
    .await;
    assert_eq!(repeated["artifact_id"], analyzed["artifact_id"]);
    let mut stale = doc.clone();
    stale["revision"] = json!("f".repeat(64));
    reject(&c, "inspect_rna_annotations", json!({"document":stale})).await;
    let mut wrong = d.clone();
    wrong["sites"][0]["nucleoside"] = json!("Canonical");
    reject(
        &c,
        "analyze_rna_chemistry",
        json!({"document":doc,"declaration":wrong}),
    )
    .await;
    let mut wrong = d.clone();
    wrong["annotation_reviews"][0]["expected_feature_revision"] = json!("f".repeat(64));
    reject(
        &c,
        "analyze_rna_chemistry",
        json!({"document":doc,"declaration":wrong}),
    )
    .await;
    let mut old = d.clone();
    old["version"] = json!(1);
    reject(
        &c,
        "analyze_rna_chemistry",
        json!({"document":doc,"declaration":old}),
    )
    .await;
    let dna = call(&c, "import_document", json!({"text":">dna\nAC\n"})).await;
    let dna = json!({"id":dna["documents"][0]["id"],"revision":dna["documents"][0]["revision"]});
    reject(&c, "inspect_rna_annotations", json!({"document":dna})).await;
    let after = call(&c, "list_documents", json!({})).await;
    assert_eq!(
        after["total"].as_u64().unwrap(),
        before["total"].as_u64().unwrap() + 1
    );
    c.cancel().await.unwrap();
    server.await.unwrap();
}
