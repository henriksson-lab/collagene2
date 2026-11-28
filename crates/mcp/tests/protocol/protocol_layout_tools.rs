#[path = "../../../core/tests/common/protocol_layout.rs"]
mod fixture;
use base64::{Engine as _, engine::general_purpose::STANDARD};
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
            CallToolRequestParams::new("protocol_layout")
                .with_arguments(args.as_object().unwrap().clone()),
        )
        .await;
    assert!(r.is_err() || r.unwrap().is_error == Some(true));
}
#[tokio::test]
async fn actual_binary_workspace_preview_capture_exports_and_explicit_projection_keep_identity() {
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
    assert_eq!(
        c.list_tools(Default::default()).await.unwrap().tools.len(),
        collagene_mcp::CURRENT_TOOL_COUNT
    );
    let (w, r) = fixture::fixture();
    let original = collagene_core::documents::workspace::save_workspace(&w).unwrap();
    let imported = call(
        &c,
        "import_workspace_archive",
        json!({"data_base64":STANDARD.encode(&original)}),
    )
    .await;
    let reference =
        json!({"id":imported["workspace"]["id"],"revision":imported["workspace"]["revision"]});
    let preview = call(
        &c,
        "protocol_layout",
        json!({"workspace":reference,"task":{"Preview":{"request":r}}}),
    )
    .await;
    assert_eq!(preview["verified"], false);
    assert_eq!(preview["assignment_count"], 1);
    let mut changed = r.clone();
    changed.assignments[0].well.column = 1;
    reject(&c,json!({"workspace":reference,"task":{"Capture":{"request":changed,"expected_context_revision":preview["context_revision"]}}})).await;
    let saved=call(&c,"protocol_layout",json!({"workspace":reference,"task":{"Capture":{"request":r,"expected_context_revision":preview["context_revision"]}}})).await;
    let saved_ref =
        json!({"id":saved["workspace"]["id"],"revision":saved["workspace"]["revision"]});
    let hash = saved["receipt"]["attachment_hash"].clone();
    let list = call(
        &c,
        "protocol_layout",
        json!({"workspace":saved_ref,"task":{"List":{"offset":0,"limit":20}}}),
    )
    .await;
    assert_eq!(list["total"], 1);
    for format in ["Csv", "Html", "Pdf"] {
        let export = call(
            &c,
            "protocol_layout",
            json!({"workspace":saved_ref,"task":{"Export":{"layout_hash":hash,"format":format}}}),
        )
        .await;
        let data = call(
            &c,
            "read_artifact",
            json!({"artifact_id":export["artifact_id"],"start":0,"max_bytes":16384}),
        )
        .await;
        let bytes = STANDARD.decode(data["text"].as_str().unwrap()).unwrap();
        assert!(!bytes.is_empty());
        if format == "Pdf" {
            assert!(bytes.starts_with(b"%PDF"));
        } else {
            assert!(String::from_utf8(bytes).unwrap().contains("clone-A"));
        }
    }
    let review = call(
        &c,
        "protocol_layout",
        json!({"workspace":saved_ref,"task":{"ReviewPublication":{"layout_hash":hash}}}),
    )
    .await;
    assert_eq!(review["documents"], 1);
    reject(&c,json!({"workspace":saved_ref,"task":{"ProjectPublication":{"layout_hash":hash,"expected_scope_revision":"0".repeat(64)}}})).await;
    let projected=call(&c,"protocol_layout",json!({"workspace":saved_ref,"task":{"ProjectPublication":{"layout_hash":hash,"expected_scope_revision":review["scope_revision"]}}})).await;
    assert!(projected["workspace"].is_object());
    assert_eq!(
        collagene_core::documents::workspace::save_workspace(&w).unwrap(),
        original
    );
    c.cancel().await.unwrap();
    server.await.unwrap();
}
