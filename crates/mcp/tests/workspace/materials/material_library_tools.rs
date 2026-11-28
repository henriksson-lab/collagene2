#[path = "../../../../core/tests/common/material_library.rs"]
mod material_library;
use base64::{Engine, engine::general_purpose::STANDARD};
use collagene_core::workspace;
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
async fn reject(c: &Client, args: Value) {
    let r = c
        .call_tool(
            CallToolRequestParams::new("inspect_combinatorial_construct")
                .with_arguments(args.as_object().unwrap().clone()),
        )
        .await;
    assert!(r.is_err() || r.unwrap().is_error == Some(true));
}
#[tokio::test]
async fn compact_material_selection_hydrates_exact_saved_source_and_preserves_independent_bottom_target()
 {
    let (a, b) = tokio::io::duplex(65_536);
    let server = tokio::spawn(async move {
        CollageneServer::default()
            .serve(a)
            .await
            .unwrap()
            .waiting()
            .await
            .unwrap();
    });
    let c: Client = ().serve(b).await.unwrap();
    let (w, d, row) = material_library::fixture();
    let imported = call(
        &c,
        "import_workspace_archive",
        json!({"data_base64":STANDARD.encode(workspace::save_workspace(&w).unwrap())}),
    )
    .await;
    let wr = json!({"id":imported["workspace"]["id"],"revision":imported["workspace"]["revision"]});
    let mut declared = serde_json::to_value(&d).unwrap();
    declared["material_bindings"] = json!([{"event":d.material_bindings[0].event,"selection":d.material_bindings[0].origin.selection()}]);
    let inspected = call(
        &c,
        "inspect_combinatorial_construct",
        json!({"workspace":wr,"design":declared,"construct":row}),
    )
    .await;
    assert_eq!(inspected["target"]["length"], 10);
    assert_eq!(inspected["verified"], false);
    let mut text = String::new();
    let mut start = 0;
    loop {
        let page = call(
            &c,
            "read_artifact",
            json!({"artifact_id":inspected["artifact_id"],"start":start,"max_bytes":16384}),
        )
        .await;
        text.push_str(page["text"].as_str().unwrap());
        if page["done"] == true {
            break;
        }
        start = page["next_start"].as_u64().unwrap();
    }
    let evidence: Value = serde_json::from_str(&text).unwrap();
    assert_eq!(
        evidence["independent_target"]["sequence"],
        serde_json::to_value(collagene_model::Sequence::dna(material_library::EXPECTED).unwrap())
            .unwrap()
    );
    let mut forged = declared.clone();
    forged["material_bindings"][0]["event"] =
        serde_json::to_value(material_library::event(&w, 2)).unwrap();
    reject(&c, json!({"workspace":wr,"design":forged,"construct":row})).await;
    reject(&c, json!({"workspace":wr,"design":d,"construct":row})).await; // closed client capsule is not a source selector
    c.cancel().await.unwrap();
    server.await.unwrap();
}
