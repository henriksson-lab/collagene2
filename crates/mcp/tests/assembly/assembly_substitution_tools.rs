//! Independent expected-base declarations can use repeated immutable source roles.
use collagene_core::cloning::assembly::substitution_intent::Report;
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
async fn declared_substitutions_keep_distinct_repeated_roles_and_reject_stale_sources_before_publication()
 {
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
            .any(|t| t.name == "build_substituted_assembly_intent")
    );
    let imported = call(
        &c,
        "import_document",
        json!({"text":">original\nAACCGGTT\n"}),
    )
    .await;
    let d = &imported["documents"][0];
    let source = json!({"id":d["id"],"revision":d["revision"]});
    // Two independent copies of the same selected source bases have different
    // roles. Their position1 requirements differ without mutating the original.
    let args = json!({"origin_workspace":null,"inputs":[{"id":"source","document":source}],"request":{"base":{"name":"Independent two-copy mutant requirement","fragments":[{"role":"left","input_id":"source","selector":{"ReferenceRange":{"start":0,"end":4,"reverse":false}}},{"role":"right","input_id":"source","selector":{"ReferenceRange":{"start":0,"end":4,"reverse":false}}}],"circular":false,"material_bindings":[]},"junction_edits":[],"substitutions":[{"role":"left","offset":1,"expected":"A","replacement":"T"},{"role":"right","offset":1,"expected":"A","replacement":"G"}]}});
    let result = call(&c, "build_substituted_assembly_intent", args.clone()).await;
    assert_eq!(result["state"], "intended-reference");
    assert_eq!(result["verified"], false);
    assert_eq!(result["preparation_simulated"], false);
    assert_eq!(result["intended_length"], 8);
    let full = artifact(&c, &result["artifact_id"]).await;
    let report: Report = serde_json::from_value(full["report"].clone()).unwrap();
    report
        .validate_binding(result["context_revision"].as_str().unwrap())
        .unwrap();
    assert_eq!(report.target.sequence.as_str(), "ATCCAGCC");
    assert_eq!(report.fragments.len(), 2);
    assert_eq!(report.applied.len(), 2);
    assert_eq!(report.applied[0].role, "left");
    assert_eq!(report.applied[1].role, "right");
    assert_eq!(
        report.applied[0].source_position,
        collagene_core::cloning::assembly::intent::ContributionPosition::Reference {
            position: 1,
            reverse: false
        }
    );
    assert_eq!(
        report.applied[1].source_position,
        collagene_core::cloning::assembly::intent::ContributionPosition::Reference {
            position: 1,
            reverse: false
        }
    );
    assert_eq!(report.applied[0].output_position, 1);
    assert_eq!(report.applied[1].output_position, 5);
    for change in &report.applied {
        assert_eq!(change.input_id, "source");
        assert_eq!(change.input_revision, source["revision"].as_str().unwrap());
    }
    assert!(collagene_core::chemistry::has_unmapped_physical_evidence(
        &report.target
    ));
    assert_eq!(
        call(
            &c,
            "read_sequence_region",
            json!({"document":result["target"],"start":0,"end":8})
        )
        .await["sequence"],
        "ATCCAGCC"
    );
    let retry = call(&c, "build_substituted_assembly_intent", args.clone()).await;
    assert_eq!(retry["artifact_id"], result["artifact_id"]);
    assert_eq!(retry["target"], result["target"]);
    let mut stale_expected = args.clone();
    stale_expected["request"]["base"]["name"] = json!("Must not be retained");
    stale_expected["request"]["substitutions"][0]["expected"] = json!("C");
    reject(&c, "build_substituted_assembly_intent", stale_expected).await;
    let mut stale_source = args;
    stale_source["inputs"][0]["document"]["revision"] = json!("0".repeat(64));
    reject(&c, "build_substituted_assembly_intent", stale_source).await;
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
