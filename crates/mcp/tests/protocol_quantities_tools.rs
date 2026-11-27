#[path = "../../core/tests/common/protocol_quantities.rs"]
mod fixture;
use base64::{Engine as _, engine::general_purpose::STANDARD};
use collagene_core::{protocol_quantities as q, run_lineage, workspace};
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
            CallToolRequestParams::new("protocol_quantities")
                .with_arguments(args.as_object().unwrap().clone()),
        )
        .await;
    assert!(r.is_err() || r.unwrap().is_error == Some(true));
}
async fn artifact(c: &Client, id: &Value) -> Vec<u8> {
    let mut text = String::new();
    let mut start = 0;
    let encoding;
    loop {
        let part = call(
            c,
            "read_artifact",
            json!({"artifact_id":id,"start":start,"max_bytes":16384}),
        )
        .await;
        text.push_str(part["text"].as_str().unwrap());
        if part["done"] == true {
            encoding = part["encoding"].as_str().unwrap().to_owned();
            break;
        }
        start = part["next_start"].as_u64().unwrap();
    }
    let bytes = if encoding == "base64" {
        STANDARD.decode(text).unwrap()
    } else {
        text.into_bytes()
    };
    assert_eq!(blake3::hash(&bytes).to_hex().as_str(), id.as_str().unwrap());
    bytes
}
#[tokio::test]
async fn real_session_saved_run_quantities_are_atomic_documentary_artifacts_not_new_science() {
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
    let original = workspace::save_workspace(&w).unwrap();
    let imported = call(
        &c,
        "import_workspace_archive",
        json!({"data_base64":STANDARD.encode(&original)}),
    )
    .await;
    let reference =
        json!({"id":imported["workspace"]["id"],"revision":imported["workspace"]["revision"]});
    // Exercise the discoverable tool example over the actual SDK transport,
    // replacing only exact native identity pins rather than writing another schema.
    let tools = c.list_tools(Default::default()).await.unwrap();
    let tool = tools
        .tools
        .iter()
        .find(|t| t.name == "protocol_quantities")
        .unwrap();
    let description = tool.description.as_ref().unwrap();
    let sample = description
        .split("MINIMAL_PREVIEW=")
        .nth(1)
        .unwrap()
        .split(" END_PREVIEW")
        .next()
        .unwrap();
    let mut example: Value = serde_json::from_str(sample).unwrap();
    example["workspace"] = reference.clone();
    let request = &mut example["task"]["Preview"]["request"];
    request["plan_source"] = serde_json::to_value(&r.plan_source).unwrap();
    request["plan_revision"] = json!(r.plan_revision);
    request["target"] = serde_json::to_value(&r.target).unwrap();
    request["inputs"] = serde_json::to_value(&r.inputs).unwrap();
    let example_result = call(&c, "protocol_quantities", example).await;
    assert_eq!(example_result["arithmetic_complete"], true);
    let example_preview: q::Preview =
        serde_json::from_slice(&artifact(&c, &example_result["artifact_id"]).await).unwrap();
    assert_eq!(
        example_preview.reactions[0].fill.as_ref().unwrap().volume,
        fixture::volume(18)
    );
    let preview = call(
        &c,
        "protocol_quantities",
        json!({"workspace":reference,"task":{"Preview":{"request":r}}}),
    )
    .await;
    assert_eq!(preview["verified"], false);
    assert_eq!(preview["arithmetic_complete"], true);
    let full: q::Preview =
        serde_json::from_slice(&artifact(&c, &preview["artifact_id"]).await).unwrap();
    assert_eq!(
        full.reactions[0].fill.as_ref().unwrap().volume,
        fixture::volume(15)
    );
    let mut changed = r.clone();
    changed.reactions[0].final_volume = Some(fixture::volume(21));
    reject(&c,json!({"workspace":reference,"task":{"Capture":{"request":changed,"expected_context_revision":preview["context_revision"]}}})).await;
    for target in [false, true] {
        let mut alias = r.clone();
        if target {
            alias.target = fixture::pin(&w, 3);
        } else {
            alias.inputs[0].event = fixture::pin(&w, 2);
        }
        reject(
            &c,
            json!({"workspace":reference,"task":{"Preview":{"request":alias}}}),
        )
        .await;
    }
    let captured=call(&c,"protocol_quantities",json!({"workspace":reference,"task":{"Capture":{"request":r,"expected_context_revision":preview["context_revision"]}}})).await;
    assert_eq!(captured["receipt"]["verified"], false);
    assert_eq!(captured["receipt"]["arithmetic_complete"], true);
    let saved_ref =
        json!({"id":captured["workspace"]["id"],"revision":captured["workspace"]["revision"]});
    let hash = captured["receipt"]["attachment_hash"].clone();
    let native_id = captured["native_artifact_id"].clone();
    let native = artifact(&c, &native_id).await;
    let reopened = workspace::load_workspace(&native).unwrap();
    assert_eq!(
        serde_json::to_value(&reopened.documents).unwrap(),
        serde_json::to_value(&w.documents).unwrap()
    );
    assert_eq!(
        serde_json::to_value(&reopened.history).unwrap(),
        serde_json::to_value(&w.history).unwrap()
    );
    assert_eq!(
        serde_json::to_value(run_lineage::list_runs(&reopened).unwrap()).unwrap(),
        serde_json::to_value(run_lineage::list_runs(&w).unwrap()).unwrap()
    );
    let list = call(
        &c,
        "protocol_quantities",
        json!({"workspace":saved_ref,"task":{"List":{"offset":0,"limit":20}}}),
    )
    .await;
    assert_eq!(list["total"], 1);
    let opened = call(
        &c,
        "protocol_quantities",
        json!({"workspace":saved_ref,"task":{"Open":{"quantity_hash":hash}}}),
    )
    .await;
    assert_ne!(
        opened["context_revision"],
        opened["historical_context_revision"]
    );
    assert_eq!(opened["arithmetic_complete"], true);
    let stored: q::Preview =
        serde_json::from_slice(&artifact(&c, &opened["artifact_id"]).await).unwrap();
    assert_eq!(stored.request.target, r.target);
    assert_eq!(stored.request.inputs, r.inputs);
    for format in ["Csv", "Html", "Pdf"] {
        let exported = call(
            &c,
            "protocol_quantities",
            json!({"workspace":saved_ref,"task":{"Export":{"quantity_hash":hash,"format":format}}}),
        )
        .await;
        assert_eq!(exported["verified"], false);
        let bytes = artifact(&c, &exported["artifact_id"]).await;
        if format == "Pdf" {
            assert!(bytes.starts_with(b"%PDF-"));
        } else {
            assert!(
                String::from_utf8(bytes)
                    .unwrap()
                    .contains("PRIVATE EXPLICIT QUANTITY PLAN")
            );
        }
    }
    let mut draft = r.clone();
    draft.reactions[0].final_volume = None;
    let draft = call(
        &c,
        "protocol_quantities",
        json!({"workspace":reference,"task":{"Preview":{"request":draft}}}),
    )
    .await;
    assert_eq!(draft["arithmetic_complete"], false);
    assert!(draft["finding_count"].as_u64().unwrap() > 0);
    assert_eq!(workspace::save_workspace(&w).unwrap(), original);
    c.cancel().await.unwrap();
    server.await.unwrap();
}
