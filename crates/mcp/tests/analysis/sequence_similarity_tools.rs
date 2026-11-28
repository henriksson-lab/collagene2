//! Real SDK over memory pipes; no network, browser, external aligner or filesystem.
use collagene_mcp::server::CollageneServer;
use rmcp::{RoleClient, ServiceExt, model::CallToolRequestParams, service::RunningService};
use serde_json::{Value, json};
type Client = RunningService<RoleClient, ()>;
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
#[tokio::test]
async fn immutable_subject_discovery_retains_bounds_and_never_installs_documents() {
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
    let tools = c.list_tools(Default::default()).await.unwrap();
    assert_eq!(tools.tools.len(), collagene_mcp::CURRENT_TOOL_COUNT);
    assert!(
        tools
            .tools
            .iter()
            .any(|t| t.name == "search_sequence_library")
    );
    let imported = call(
        &c,
        "import_document",
        json!({"text":">query\nCGAA\n>subject\nTTCGAATT\n>other\nACGTTGAC\n"}),
    )
    .await;
    let reference = |name: &str| {
        let row = imported["documents"]
            .as_array()
            .unwrap()
            .iter()
            .find(|d| d["name"] == name)
            .unwrap();
        json!({"id":row["id"],"revision":row["revision"]})
    };
    let q = reference("query");
    let s = reference("subject");
    let other = reference("other");
    let before = call(&c, "list_documents", json!({})).await["total"].clone();
    let report = call(
        &c,
        "search_sequence_library",
        json!({"query":q,"subjects":[s,other]}),
    )
    .await;
    assert_eq!(report["complete"], true);
    assert_eq!(report["hits"][0]["score"], 8);
    assert_eq!(report["subjects"].as_array().unwrap().len(), 2);
    let artifact = call(
        &c,
        "read_artifact",
        json!({"artifact_id":report["artifact_id"],"start":0,"max_bytes":16384}),
    )
    .await;
    let evidence: Value = serde_json::from_str(artifact["text"].as_str().unwrap()).unwrap();
    assert_eq!(evidence["query"], q);
    assert_eq!(
        evidence["report"]["engine_version"],
        "sequence-similarity-2"
    );
    let options = collagene_core::sequence::similarity::Options {
        max_work: 1,
        ..Default::default()
    };
    let partial = call(
        &c,
        "search_sequence_library",
        json!({"query":q,"subjects":[s],"options":options}),
    )
    .await;
    assert_eq!(partial["complete"], false);
    assert!(
        !partial["subjects"][0]["limits"]
            .as_array()
            .unwrap()
            .is_empty()
    );
    let mut stale = q.clone();
    stale["revision"] = json!("a".repeat(64));
    let reject = c
        .call_tool(
            CallToolRequestParams::new("search_sequence_library").with_arguments(
                json!({"query":stale,"subjects":[s]})
                    .as_object()
                    .unwrap()
                    .clone(),
            ),
        )
        .await;
    match reject {
        Ok(reply) => assert_eq!(reply.is_error, Some(true)),
        Err(rmcp::service::ServiceError::McpError(error)) => {
            assert_eq!(error.code.0, -32602);
            assert!(error.message.contains("stale"));
        }
        Err(error) => panic!("unexpected SDK transport failure: {error}"),
    }
    assert_eq!(call(&c, "list_documents", json!({})).await["total"], before);
    c.cancel().await.unwrap();
    server.await.unwrap();
}
