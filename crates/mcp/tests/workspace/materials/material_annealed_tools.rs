//! Actual MCP SDK: both supplied oligos survive atomic native installation and restart.
use base64::{Engine as _, engine::general_purpose::STANDARD};
use collagene_core::{
    chemistry::AnnealOptions,
    material_derivation as native,
    run_lineage::EventBinding,
    workspace::{self, Workspace},
};
use collagene_mcp::CollageneServer;
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
    if let Ok(r) = c
        .call_tool(
            CallToolRequestParams::new("material_derivation")
                .with_arguments(args.as_object().unwrap().clone()),
        )
        .await
    {
        assert_eq!(r.is_error, Some(true), "{:?}", r.content)
    }
}
async fn artifact(c: &Client, id: &Value) -> (String, String) {
    let (mut start, mut text) = (0, String::new());
    loop {
        let p = call(
            c,
            "read_artifact",
            json!({"artifact_id":id,"start":start,"max_bytes":16384}),
        )
        .await;
        text.push_str(p["text"].as_str().unwrap());
        if p["done"] == true {
            return (text, p["encoding"].as_str().unwrap().into());
        }
        start = p["next_start"].as_u64().unwrap();
    }
}
fn event(w: &Workspace, i: usize) -> EventBinding {
    let d = &w.documents[i];
    EventBinding {
        document_id: d.id.clone(),
        event_id: d.head_event.clone().unwrap(),
        revision: d.head.clone(),
    }
}
fn reference(v: &Value) -> Value {
    json!({"id":v["workspace"]["id"],"revision":v["workspace"]["revision"]})
}
#[tokio::test]
async fn two_role_adapter_requires_exact_small_source_handles_and_replays_after_new_server_import()
{
    let original = Workspace::from_documents(
        "Supplied ssDNA",
        vec![
            Molecule::dna("Top", "AATTCGTACT", Topology::Linear).unwrap(),
            Molecule::dna("Bottom", "AGTACG", Topology::Linear).unwrap(),
        ],
        0,
    )
    .unwrap();
    let bytes = workspace::save_workspace(&original).unwrap();
    let top = event(&original, 0);
    let bottom = event(&original, 1);
    let options = AnnealOptions {
        min_paired_bases: 6,
        ..Default::default()
    };
    let (c, server) = setup().await;
    let imported = call(
        &c,
        "import_workspace_archive",
        json!({"data_base64":STANDARD.encode(&bytes)}),
    )
    .await;
    let input = reference(&imported);
    let before = call(&c, "list_documents", json!({})).await["total"].clone();
    let preview=call(&c,"material_derivation",json!({"workspace":input,"task":{"PreviewAnnealedOligos":{"top":top,"bottom":bottom,"options":options}}})).await;
    assert_eq!(preview["count"], 1);
    assert_eq!(preview["unique"], true);
    let (text, encoding) = artifact(&c, &preview["artifact_id"]).await;
    assert_eq!(encoding, "utf-8");
    let evidence: Value = serde_json::from_str(&text).unwrap();
    let report: native::AnnealedOligosPreview =
        serde_json::from_value(evidence["preview"].clone()).unwrap();
    report
        .validate_binding(preview["context_revision"].as_str().unwrap())
        .unwrap();
    assert_eq!(report.products[0].sequence.as_str(), "AATTCGTACT");
    let request = json!({"sources":[{"role":"top_oligo","event":top,"origin":"OriginalReference"},{"role":"bottom_oligo","event":bottom,"origin":"OriginalReference"}],"derivation":{"AnnealedOligos":{"options":options,"selection":preview["candidates"][0]["selection"]}},"name":"Retained sticky adapter","save_primers":false});
    let mut changed = request.clone();
    changed["sources"][1]["event"] = json!(top);
    reject(
        &c,
        json!({"workspace":input,"task":{"Install":{"request":changed}}}),
    )
    .await;
    let mut changed = request.clone();
    changed["save_primers"] = json!(true);
    reject(
        &c,
        json!({"workspace":input,"task":{"Install":{"request":changed}}}),
    )
    .await;
    assert_eq!(call(&c, "list_documents", json!({})).await["total"], before);
    let installed = call(
        &c,
        "material_derivation",
        json!({"workspace":input,"task":{"Install":{"request":request}}}),
    )
    .await;
    let receipt: native::Receipt = serde_json::from_value(installed["receipt"].clone()).unwrap();
    assert_eq!(receipt.sources.len(), 2);
    assert!(receipt.primer_attachment.is_none());
    let (text, encoding) = artifact(&c, &installed["native_artifact_id"]).await;
    assert_eq!(encoding, "base64");
    let archive = STANDARD.decode(text).unwrap();
    let frame = workspace::load_workspace(&archive).unwrap();
    receipt.validate_archive(&frame).unwrap();
    let fresh = native::reconstruct(&frame, &receipt.product).unwrap();
    assert_eq!(fresh.tracked().duplex.bottom.as_str(), "AGTACG");
    assert_eq!(fresh.tracked().duplex.bottom_start, 4);
    assert!(!fresh.incomplete());
    assert_eq!(workspace::save_workspace(&original).unwrap(), bytes);
    c.cancel().await.unwrap();
    server.await.unwrap();
    let (c, server) = setup().await;
    let imported = call(
        &c,
        "import_workspace_archive",
        json!({"data_base64":STANDARD.encode(&archive)}),
    )
    .await;
    let restored = reference(&imported);
    let inspected=call(&c,"inspect_material_origin",json!({"workspace":restored,"event":receipt.product,"selection":{"Derived":receipt.product}})).await;
    assert_eq!(inspected["origin_available"], true);
    c.cancel().await.unwrap();
    server.await.unwrap();
}
