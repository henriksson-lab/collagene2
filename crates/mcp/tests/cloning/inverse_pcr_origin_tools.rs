#[path = "../../../core/tests/common/inverse_pcr_origin.rs"]
mod fixture;
use base64::{Engine as _, engine::general_purpose::STANDARD};
use collagene_core::documents::workspace::{self, Workspace};
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
            CallToolRequestParams::new("design_inverse_pcr")
                .with_arguments(args.as_object().unwrap().clone()),
        )
        .await;
    assert!(r.is_err() || r.unwrap().is_error == Some(true));
}
async fn artifact(c: &Client, id: &Value) -> Value {
    let mut text = String::new();
    let mut start = 0;
    loop {
        let r = call(
            c,
            "read_artifact",
            json!({"artifact_id":id,"start":start,"max_bytes":16384}),
        )
        .await;
        text.push_str(r["text"].as_str().unwrap());
        if r["done"] == true {
            break;
        }
        start = r["next_start"].as_u64().unwrap();
    }
    assert_eq!(
        blake3::hash(text.as_bytes()).to_hex().as_str(),
        id.as_str().unwrap()
    );
    serde_json::from_str(&text).unwrap()
}
#[tokio::test]
async fn actual_saved_run_and_derived_inverse_templates_are_attested_without_status_or_head_adoption()
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
    assert_eq!(
        c.list_tools(Default::default()).await.unwrap().tools.len(),
        collagene_mcp::CURRENT_TOOL_COUNT
    );
    for derived in [false, true] {
        let (w, origin, event, source, target) = fixture::fixture(derived);
        let before = workspace::save_workspace(&w).unwrap();
        let imported = call(
            &c,
            "import_workspace_archive",
            json!({"data_base64":STANDARD.encode(&before)}),
        )
        .await;
        let wr =
            json!({"id":imported["workspace"]["id"],"revision":imported["workspace"]["revision"]});
        let args = json!({"document":{"id":source.revision(),"revision":source.revision()},"target":{"id":target.revision(),"revision":target.revision()},"origin":origin.selection(),"origin_workspace":wr,"options":fixture::options()});
        let response = call(&c, "design_inverse_pcr", args.clone()).await;
        assert_eq!(response["computed_status"], "ConsistentUnderAssumptions");
        assert_eq!(response["verified"], false);
        assert_eq!(response["source_event"], json!(event));
        let full = artifact(&c, &response["artifact_id"]).await;
        assert_eq!(full["origin"], json!(origin.selection()));
        assert_eq!(full["origin_workspace"], wr);
        assert_eq!(
            full["inverse_pcr"]["plan"]["steps"]
                .as_array()
                .unwrap()
                .len(),
            4
        );
        assert!(full["inverse_pcr"]["plan"]["steps"][1].get("Pcr").is_some());
        let verified = call(
            &c,
            "verify_plan",
            json!({"plan_revision":response["plan_revision"]}),
        )
        .await;
        assert_eq!(verified["status"], "ConsistentUnderAssumptions");
        assert_eq!(verified["steps"].as_array().unwrap().len(), 4);
        assert_eq!(
            call(&c, "design_inverse_pcr", args.clone()).await["artifact_id"],
            response["artifact_id"]
        );
        let mut missing = args.clone();
        missing.as_object_mut().unwrap().remove("origin");
        reject(&c, missing).await;
        let mut unpaired = args.clone();
        unpaired["origin_workspace"] = Value::Null;
        reject(&c, unpaired).await;
        let detached = Workspace::from_documents(
            "Equal molecule has no retained origin",
            vec![source.clone(), target.clone()],
            0,
        )
        .unwrap();
        let detached = call(
            &c,
            "import_workspace_archive",
            json!({"data_base64":STANDARD.encode(workspace::save_workspace(&detached).unwrap())}),
        )
        .await;
        let mut foreign = args.clone();
        foreign["origin_workspace"] =
            json!({"id":detached["workspace"]["id"],"revision":detached["workspace"]["revision"]});
        reject(&c, foreign).await;
        let mut wrong_source = args.clone();
        wrong_source["document"] = json!({"id":w.documents[0].head,"revision":w.documents[0].head});
        reject(&c, wrong_source).await;
        let mut wrong = args.clone();
        wrong["target"] = json!({"id":w.documents[1].head,"revision":w.documents[1].head});
        assert_eq!(
            call(&c, "design_inverse_pcr", wrong).await["computed_status"],
            "Invalid"
        );
        let html = call(
            &c,
            "export_protocol",
            json!({"plan_revision":response["plan_revision"]}),
        )
        .await;
        assert!(html["artifact_id"].is_string());
        let pdf = call(
            &c,
            "export_protocol_pdf",
            json!({"plan_revision":response["plan_revision"]}),
        )
        .await;
        assert_eq!(pdf["status"], "ConsistentUnderAssumptions");
        assert_eq!(workspace::save_workspace(&w).unwrap(), before);
    }
    c.cancel().await.unwrap();
    server.await.unwrap();
}
