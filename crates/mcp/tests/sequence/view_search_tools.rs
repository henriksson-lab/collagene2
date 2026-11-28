//! Unified search and browser tracks use the real MCP SDK and immutable refs.
use collagene_core::{
    formats::sequence as formats,
    sequence::search::{SearchMode, SearchOptions},
    tracks::view::{TrackOptions, TrackWindow},
};
use collagene_mcp::server::CollageneServer;
use collagene_model::{Molecule, Topology};
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
async fn call(client: &Client, name: &str, args: Value) -> Value {
    let result = client
        .call_tool(
            CallToolRequestParams::new(name.to_owned())
                .with_arguments(args.as_object().unwrap().clone()),
        )
        .await
        .unwrap();
    assert_ne!(result.is_error, Some(true), "{name}: {:?}", result.content);
    result.structured_content.unwrap()
}
async fn reject(client: &Client, name: &str, args: Value) {
    let result = client
        .call_tool(
            CallToolRequestParams::new(name.to_owned())
                .with_arguments(args.as_object().unwrap().clone()),
        )
        .await;
    assert!(
        result.is_err() || result.unwrap().is_error == Some(true),
        "accepted invalid {name}"
    );
}
fn reference(value: &Value) -> Value {
    json!({"id":value["id"],"revision":value["revision"]})
}
async fn import(client: &Client, name: &str, sequence: &str, topology: Topology) -> Value {
    let m = Molecule::dna(name, sequence, topology).unwrap();
    let text = formats::export_genbank(&m).unwrap();
    reference(&call(client, "import_document", json!({"text":text})).await["documents"][0])
}
async fn artifact(client: &Client, id: &Value) -> Value {
    let mut start = 0;
    let mut text = String::new();
    loop {
        let result = call(
            client,
            "read_artifact",
            json!({"artifact_id":id,"start":start,"max_bytes":16384}),
        )
        .await;
        text.push_str(result["text"].as_str().unwrap());
        if result["done"] == true {
            break;
        }
        start = result["next_start"].as_u64().unwrap();
    }
    serde_json::from_str(&text).unwrap()
}
#[tokio::test]
async fn search_has_bounded_preview_full_artifact_and_exact_revision() {
    let (client, server) = setup().await;
    let source = import(&client, "source", &"A".repeat(240), Topology::Linear).await;
    let options = SearchOptions {
        query: "A".into(),
        mode: SearchMode::Sequence,
        both_strands: false,
        ..Default::default()
    };
    let result = call(
        &client,
        "search_document",
        json!({"document":source,"options":options}),
    )
    .await;
    assert_eq!(result["hit_count"], 240);
    assert_eq!(result["hits"].as_array().unwrap().len(), 100);
    assert_eq!(result["preview_truncated"], true);
    assert_eq!(result["complete"], true);
    assert_eq!(result["source_revision"], source["revision"]);
    let full = artifact(&client, &result["artifact_id"]).await;
    assert_eq!(full["hits"].as_array().unwrap().len(), 240);
    assert_eq!(full["options"], json!(options));
    let mut stale = source.clone();
    stale["revision"] = json!("stale");
    reject(
        &client,
        "search_document",
        json!({"document":stale,"options":options}),
    )
    .await;
    let limited = SearchOptions {
        max_hits: 2,
        ..options.clone()
    };
    let result = call(
        &client,
        "search_document",
        json!({"document":source,"options":limited}),
    )
    .await;
    assert_eq!(result["hit_count"], 2);
    assert_eq!(result["complete"], false);
    let invalid = SearchOptions {
        max_work: 0,
        ..options
    };
    reject(
        &client,
        "search_document",
        json!({"document":source,"options":invalid}),
    )
    .await;
    client.cancel().await.unwrap();
    server.await.unwrap();
}
#[tokio::test]
async fn track_tool_exposes_the_same_real_cut_coordinates_and_rejects_invalid_window() {
    let (client, server) = setup().await;
    let source = import(&client, "source", "AAAAGAATTCCCC", Topology::Circular).await;
    let options = TrackOptions {
        enzyme_recognition: true,
        enzyme_cuts: true,
        enzyme_names: vec!["EcoRI".into()],
        ..Default::default()
    };
    let result = call(
        &client,
        "view_tracks",
        json!({"document":source,"window":TrackWindow{start:0,end:13},"options":options}),
    )
    .await;
    assert_eq!(result["source_revision"], source["revision"]);
    assert_eq!(result["complete"], true);
    assert_eq!(result["item_count"], 1);
    let full = artifact(&client, &result["artifact_id"]).await;
    assert_eq!(full["items"][0]["cuts"]["top"], 5);
    assert_eq!(full["items"][0]["cuts"]["bottom"], 9);
    assert_eq!(full["items"][0]["parts"][0]["start"], 4);
    reject(
        &client,
        "view_tracks",
        json!({"document":source,"window":{"start":14,"end":15},"options":options}),
    )
    .await;
    let listed = call(&client, "list_documents", json!({})).await;
    assert_eq!(listed["documents"].as_array().unwrap().len(), 1);
    client.cancel().await.unwrap();
    server.await.unwrap();
}
