//! Literal insertion/deletion requirements remain source-bound, reference-only MCP documents.
use collagene_core::assembly_edit_intent::Report;
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
    match c
        .call_tool(
            CallToolRequestParams::new(name.to_owned())
                .with_arguments(args.as_object().unwrap().clone()),
        )
        .await
    {
        Err(_) => {}
        Ok(r) => assert_eq!(
            r.is_error,
            Some(true),
            "{name} unexpectedly succeeded: {:?}",
            r
        ),
    }
}
async fn artifact(c: &Client, id: &Value) -> Value {
    let (mut at, mut out) = (0, String::new());
    loop {
        let p = call(
            c,
            "read_artifact",
            json!({"artifact_id":id,"start":at,"max_bytes":16384}),
        )
        .await;
        out.push_str(p["text"].as_str().unwrap());
        if p["done"] == true {
            break;
        }
        at = p["next_start"].as_u64().unwrap();
    }
    serde_json::from_str(&out).unwrap()
}
#[tokio::test]
async fn literal_contribution_edits_keep_original_sources_and_publish_only_unverified_intent() {
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
    let tools = c.list_tools(Default::default()).await.unwrap().tools;
    assert_eq!(tools.len(), collagene_mcp::CURRENT_TOOL_COUNT);
    assert!(
        tools
            .iter()
            .any(|t| t.name == "build_contribution_edit_intent")
    );
    let imported = call(
        &c,
        "import_document",
        json!({"text":">original\nAACCGGTT\n"}),
    )
    .await;
    let d = &imported["documents"][0];
    let source = json!({"id":d["id"],"revision":d["revision"]});
    // Original AACC plus its independently chosen reverse-complement GGTT.
    // Insert TT into left and replace GT by A in right: AATTCC + GAT.
    let args = json!({"origin_workspace":null,"inputs":[{"id":"source","document":source}],"request":{
        "base":{"base":{"name":"Independent insertion and deletion requirement","fragments":[
            {"role":"left","input_id":"source","selector":{"ReferenceRange":{"start":0,"end":4,"reverse":false}}},
            {"role":"right","input_id":"source","selector":{"ReferenceRange":{"start":0,"end":4,"reverse":true}}}],"circular":false,"material_bindings":[]},
            "junction_edits":[],"substitutions":[]},
        "edits":[{"role":"left","start":2,"end":2,"expected":"","replacement":"TT"},
            {"role":"right","start":1,"end":3,"expected":"GT","replacement":"A"}]}});
    let result = call(&c, "build_contribution_edit_intent", args.clone()).await;
    assert_eq!(result["state"], "intended-reference");
    assert_eq!(result["verified"], false);
    assert_eq!(result["preparation_simulated"], false);
    assert_eq!(result["intended_length"], 9);
    let full = artifact(&c, &result["artifact_id"]).await;
    let report: Report = serde_json::from_value(full["report"].clone()).unwrap();
    report
        .validate_binding(result["context_revision"].as_str().unwrap())
        .unwrap();
    assert_eq!(report.target.sequence.as_str(), "AATTCCGAT");
    assert_eq!(report.fragments.len(), 2);
    assert_eq!(report.applied.len(), 2);
    assert!(collagene_core::chemistry::has_unmapped_physical_evidence(
        &report.target
    ));
    assert_eq!(
        call(
            &c,
            "read_sequence_region",
            json!({"document":result["target"],"start":0,"end":9})
        )
        .await["sequence"],
        "AATTCCGAT"
    );
    assert!(full["origin_workspace"].is_null());
    assert!(
        full["report"]["request"]["base"]["base"]["material_bindings"]
            .as_array()
            .unwrap()
            .is_empty()
    );
    let retry = call(&c, "build_contribution_edit_intent", args.clone()).await;
    assert_eq!(retry["artifact_id"], result["artifact_id"]);
    assert_eq!(retry["target"], result["target"]);
    let mut wrong = args.clone();
    wrong["request"]["base"]["base"]["name"] = json!("Must not publish wrong expected bases");
    wrong["request"]["edits"][1]["expected"] = json!("AC");
    reject(&c, "build_contribution_edit_intent", wrong).await;
    let mut stale = args.clone();
    stale["inputs"][0]["document"]["revision"] = json!("0".repeat(64));
    reject(&c, "build_contribution_edit_intent", stale).await;
    let mut conflicting = args;
    conflicting["request"]["edits"]
        .as_array_mut()
        .unwrap()
        .push(json!({"role":"right","start":2,"end":3,"expected":"T","replacement":"C"}));
    reject(&c, "build_contribution_edit_intent", conflicting).await;
    assert_eq!(
        call(
            &c,
            "read_sequence_region",
            json!({"document":source,"start":0,"end":8})
        )
        .await["sequence"],
        "AACCGGTT"
    );
    assert_eq!(call(&c, "list_documents", json!({})).await["total"], 2);
    drop(c);
    server.abort();
}
