//! Real SDK calls use small retained selections, never supplied closed graph authority.
#[path = "../../../../core/tests/common/topo_graph.rs"]
mod literal;
use base64::{Engine, engine::general_purpose::STANDARD};
use collagene_core::{cloning::topo::graph_history as g, documents::workspace};
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
async fn reject(c: &Client, args: Value) {
    let r = c
        .call_tool(
            CallToolRequestParams::new("topo_graph")
                .with_arguments(args.as_object().unwrap().clone()),
        )
        .await;
    assert!(r.is_err() || r.unwrap().is_error == Some(true));
}
async fn import(c: &Client, w: &workspace::Workspace) -> Value {
    let r = call(
        c,
        "import_workspace_archive",
        json!({"data_base64":STANDARD.encode(workspace::save_workspace(w).unwrap())}),
    )
    .await;
    json!({"id":r["workspace"]["id"],"revision":r["workspace"]["revision"]})
}
#[tokio::test]
async fn native_sdk_graph_capture_is_immutable_paged_and_requires_exact_retained_selector() {
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
    let client = ().serve(b).await.unwrap();
    let (w, request) = literal::fixture();
    let original = import(&client, &w).await;
    let result = call(
        &client,
        "topo_graph",
        json!({"workspace":original,"task":{"Capture":{"request":request}}}),
    )
    .await;
    assert_eq!(result["receipt"]["verified"], false);
    assert!(result["native_artifact_id"].is_string());
    let saved = json!({"id":result["workspace"]["id"],"revision":result["workspace"]["revision"]});
    let selection = result["receipt"]["selection"].clone();
    let list = call(
        &client,
        "topo_graph",
        json!({"workspace":saved,"task":{"List":{"offset":0,"limit":1}}}),
    )
    .await;
    assert_eq!(list["total"], 1);
    assert_eq!(list["entries"][0]["source_bases"], 26);
    assert_eq!(
        list["entries"][0]["remaining_activations"]
            .as_array()
            .unwrap()
            .len(),
        2
    );
    let old = call(
        &client,
        "topo_graph",
        json!({"workspace":original,"task":{"List":{"offset":0,"limit":1}}}),
    )
    .await;
    assert_eq!(old["total"], 0);
    reject(
        &client,
        json!({"workspace":original,"task":{"Inspect":{"selection":selection}}}),
    )
    .await;
    let inspected = call(
        &client,
        "topo_graph",
        json!({"workspace":saved,"task":{"Inspect":{"selection":selection}}}),
    )
    .await;
    assert_eq!(inspected["verified"], false);
    assert_eq!(inspected["summary"]["chains"], 2);
    let mut wrong = selection.clone();
    wrong["state_revision"] = json!("0".repeat(64));
    reject(
        &client,
        json!({"workspace":saved,"task":{"Inspect":{"selection":wrong}}}),
    )
    .await;
    // A JSON capsule/report is not an accepted graph selection shape.
    let first = literal::captured();
    let capsule = g::documentary_capsule(&first.workspace, &first.receipt.selection).unwrap();
    reject(
        &client,
        json!({"workspace":saved,"task":{"Inspect":{"selection":capsule}}}),
    )
    .await;
    // Resume capture carries only the retained selection; native hydration restores
    // the whole graph under this exact immutable session workspace.
    let partial = literal::partial(&first.workspace, &first.receipt.selection);
    let projected = json!({"recipe":{"schema_version":1,"materials":[{"SavedGraph":{"material_id":"retained vector","selection":selection}},partial.request.recipe.materials[1]]},"steps":partial.request.steps,"engine_limits":partial.request.engine_limits,"projection_options":partial.request.projection_options});
    let resumed=call(&client,"topo_graph",json!({"workspace":saved,"task":{"Capture":{"request":{"name":"One transfer retained whole graph","bindings":partial.bindings,"request":projected}}}})).await;
    assert_eq!(resumed["receipt"]["summary"]["source_bases"], 40);
    assert_eq!(
        resumed["receipt"]["summary"]["remaining_activations"]
            .as_array()
            .unwrap()
            .len(),
        1
    );
    assert_eq!(
        resumed["receipt"]["summary"]["released_activations"]
            .as_array()
            .unwrap()
            .len(),
        1
    );
    let unchanged = call(
        &client,
        "topo_graph",
        json!({"workspace":saved,"task":{"List":{"offset":0,"limit":1}}}),
    )
    .await;
    assert_eq!(unchanged["total"], 1);
    client.cancel().await.unwrap();
    server.await.unwrap();
}
