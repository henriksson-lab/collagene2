#[path = "../../../../core/tests/common/domestication_oe.rs"]
mod fixture;
use base64::{Engine as _, engine::general_purpose::STANDARD};
use collagene_core::{
    verification,
    workspace::{self, Workspace},
};
use collagene_mcp::CollageneServer;
use collagene_model::Topology;
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
    let r = c
        .call_tool(
            CallToolRequestParams::new(name.to_owned())
                .with_arguments(args.as_object().unwrap().clone()),
        )
        .await;
    assert!(r.is_err() || r.unwrap().is_error == Some(true));
}
async fn artifact(c: &Client, id: &Value) -> (Vec<u8>, String) {
    let mut text = String::new();
    let mut start = 0;
    let encoding;
    loop {
        let r = call(
            c,
            "read_artifact",
            json!({"artifact_id":id,"start":start,"max_bytes":16384}),
        )
        .await;
        text.push_str(r["text"].as_str().unwrap());
        if r["done"] == true {
            encoding = r["encoding"].as_str().unwrap().to_owned();
            break;
        }
        start = r["next_start"].as_u64().unwrap();
    }
    let bytes = if encoding == "base64" {
        STANDARD.decode(text).unwrap()
    } else {
        text.into_bytes()
    };
    assert_eq!(blake3::hash(&bytes).to_hex().as_str(), id.as_str().unwrap());
    (bytes, encoding)
}
#[tokio::test]
async fn real_original_linear_and_circle_oe_keep_checked_target_and_actual_primers_in_one_candidate()
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
    for topology in [Topology::Linear, Topology::Circular] {
        let (input, request) = fixture::fixture(topology);
        let w =
            Workspace::from_documents("Original explicit source", vec![input.molecule.clone()], 0)
                .unwrap();
        let original = workspace::save_workspace(&w).unwrap();
        call(
            &c,
            "import_workspace_archive",
            json!({"data_base64":STANDARD.encode(&original)}),
        )
        .await;
        let source = json!({"id":input.revision,"revision":input.revision});
        let before = call(&c, "export_document", source.clone()).await;
        // Freeze and assert the literal independent patch target before physical planning.
        let checked = call(
            &c,
            "check_domestication",
            json!({"source":source,"request":request.design,"selection":request.selection}),
        )
        .await;
        let target = json!({"id":checked["target"]["id"],"revision":checked["target"]["revision"]});
        assert_eq!(
            call(
                &c,
                "read_sequence_region",
                json!({"document":target,"start":0,"end":40})
            )
            .await["sequence"],
            fixture::INTENDED
        );
        let prepared = call(
            &c,
            "prepare_domestication",
            json!({"source":source,"request":request,"origin_workspace":null}),
        )
        .await;
        assert_eq!(prepared["verified"], false);
        assert_eq!(prepared["proposal_ready"], true);
        assert_eq!(prepared["state"], "candidate");
        assert_eq!(prepared["method"], "OverlapExtension");
        assert_eq!(prepared["target"]["revision"], target["revision"]);
        assert_eq!(
            prepared["verification"]["status"],
            "ConsistentUnderAssumptions"
        );
        assert_eq!(
            prepared["operation_count"],
            if topology == Topology::Circular { 6 } else { 4 }
        );
        let (data, encoding) = artifact(&c, &prepared["artifact_id"]).await;
        assert_eq!(encoding, "utf-8");
        let full: Value = serde_json::from_slice(&data).unwrap();
        let native: collagene_core::domestication_preparation::Prepared =
            serde_json::from_value(full["prepared"].clone()).unwrap();
        assert_eq!(native.target.sequence.as_str(), fixture::INTENDED);
        assert!(full["prepared"].get("preparation").is_none());
        assert_eq!(
            full["prepared"]["method"]["OverlapExtension"]["designs"]
                .as_array()
                .unwrap()
                .len(),
            3
        );
        let plan = json!({"plan_revision":prepared["plan_revision"]});
        assert_eq!(
            call(&c, "verify_plan", plan.clone()).await["status"],
            "ConsistentUnderAssumptions"
        );
        let html = call(&c, "export_protocol", plan.clone()).await;
        let (html, _) = artifact(&c, &html["artifact_id"]).await;
        assert!(String::from_utf8(html).unwrap().contains("PCR"));
        let pdf = call(&c, "export_protocol_pdf", plan.clone()).await;
        let (pdf, _) = artifact(&c, &pdf["artifact_id"]).await;
        assert!(pdf.starts_with(b"%PDF"));
        let oligos = call(
            &c,
            "export_plan_oligos",
            json!({"plan_revision":prepared["plan_revision"],"format":"Csv"}),
        )
        .await;
        assert_eq!(oligos["complete"], true);
        assert_eq!(oligos["use_count"], 6);
        let candidate = call(&c, "export_candidate", plan).await;
        let (bytes, _) = artifact(&c, &candidate["artifact_id"]).await;
        let candidate = collagene_core::transfer::load_candidate(&bytes).unwrap();
        assert_eq!(
            verification::verify(&candidate.plan, &native.target).status,
            verification::Status::ConsistentUnderAssumptions
        );
        let wrong=call(&c,"propose_plan",json!({"name":"Wrong independent original target","target":source,"inputs":[{"id":"template","document":source}],"steps":native.plan.steps})).await;
        assert_eq!(
            call(
                &c,
                "verify_plan",
                json!({"plan_revision":wrong["plan_revision"]})
            )
            .await["status"],
            "Invalid"
        );
        let mut stale = source.clone();
        stale["revision"] = json!("0".repeat(64));
        reject(
            &c,
            "prepare_domestication",
            json!({"source":stale,"request":request,"origin_workspace":null}),
        )
        .await;
        let mut bad = request.clone();
        bad.selection.patches[0].replacement = "G".into();
        reject(
            &c,
            "prepare_domestication",
            json!({"source":source,"request":bad,"origin_workspace":null}),
        )
        .await;
        assert_eq!(call(&c, "export_document", source).await, before);
        assert_eq!(workspace::save_workspace(&w).unwrap(), original);
    }
    c.cancel().await.unwrap();
    server.await.unwrap();
}
