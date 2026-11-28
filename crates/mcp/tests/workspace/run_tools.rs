//! Real SDK routing for immutable retained-run session workflows.
use base64::{Engine, engine::general_purpose::STANDARD};
use collagene_core::workspace::{self, Workspace};
use collagene_mcp::CollageneServer;
use collagene_model::{Molecule, Topology};
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
    assert_ne!(r.is_error, Some(true), "{name}:{:?}", r.content);
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
async fn artifact(c: &Client, id: &Value) -> Value {
    let (mut at, mut text) = (0, String::new());
    loop {
        let r = call(
            c,
            "read_artifact",
            json!({"artifact_id":id,"start":at,"max_bytes":16384}),
        )
        .await;
        text.push_str(r["text"].as_str().unwrap());
        if r["done"] == true {
            break;
        }
        at = r["next_start"].as_u64().unwrap();
    }
    serde_json::from_str(&text).unwrap()
}
#[tokio::test]
async fn immutable_run_record_list_replay_artifacts_and_explicit_event_bindings() {
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
    let source = Molecule::dna("source", "AACCGG", Topology::Linear).unwrap();
    let target = Molecule::dna("target", "AACCGG", Topology::Linear).unwrap();
    let w = Workspace::from_documents("MCP run", vec![source.clone(), target.clone()], 1).unwrap();
    // Optional development fixture used by the real executable stdio smoke.
    if let Ok(path) = std::env::var("COLLAGENE_RUN_SMOKE_FIXTURE") {
        std::fs::write(path, workspace::save_workspace(&w).unwrap()).unwrap();
    }
    let r = call(
        &c,
        "import_workspace_archive",
        json!({"data_base64":STANDARD.encode(workspace::save_workspace(&w).unwrap())}),
    )
    .await;
    let id = r["workspace"]["id"].clone();
    let wr = json!({"id":id,"revision":id});
    let listed = call(&c, "list_workspace_runs", wr.clone()).await;
    assert_eq!(listed["run_count"], 0);
    assert_eq!(listed["events"].as_array().unwrap().len(), 2);
    let events = artifact(&c, &listed["artifact_id"]).await;
    assert_eq!(events["events"].as_array().unwrap().len(), 2);
    let doc = |m: &Molecule| json!({"id":m.revision(),"revision":m.revision()});
    let p=call(&c,"propose_plan",json!({"name":"identity assembly","target":doc(&target),"inputs":[{"id":"s","document":doc(&source)}],"steps":[{"Overlap":{"fragments":[{"Input":"s"}],"overlaps":[],"circular":false}}]})).await;
    let bind = |i: usize| json!({"document_id":w.documents[i].id,"event_id":w.documents[i].head_event,"revision":w.documents[i].head});
    let args = json!({"workspace":wr,"plan_revision":p["plan_revision"],"target_binding":bind(1),"input_bindings":[{"input_id":"s","event":bind(0)}]});
    let mut wrong = args.clone();
    wrong["input_bindings"][0]["event"]["event_id"] =
        json!("event-0000000000000000000000000000000000000000000000000000000000000000");
    reject(&c, "record_workspace_run", wrong).await;
    let r = call(&c, "record_workspace_run", args.clone()).await;
    assert_ne!(r["workspace"]["id"], id);
    assert_eq!(r["receipt"]["products"].as_array().unwrap().len(), 1);
    let full = artifact(&c, &r["artifact_id"]).await;
    let retained: collagene_core::run_lineage::RunCommit = serde_json::from_value(full).unwrap();
    retained.validate().unwrap();
    assert_eq!(retained.workspace.documents.len(), 3);
    assert_eq!(retained.workspace.documents[0].head, w.documents[0].head);
    assert_eq!(
        call(&c, "list_workspace_runs", wr.clone()).await["run_count"],
        0
    );
    let next = json!({"id":r["workspace"]["id"],"revision":r["workspace"]["revision"]});
    let listed = call(&c, "list_workspace_runs", next.clone()).await;
    assert_eq!(listed["run_count"], 1);
    let historical = call(
        &c,
        "inspect_workspace_run_evidence",
        json!({"workspace":next,"run_id":r["receipt"]["run_id"]}),
    )
    .await;
    assert_eq!(historical["state"], "historical");
    assert_eq!(historical["verified"], false);
    assert_eq!(historical["scientific_step_count"], 0);
    assert!(
        historical["artifact_id"].is_null(),
        "ordinary run has no fabricated scientific graph"
    );
    let replay = call(
        &c,
        "replay_workspace_run",
        json!({"workspace":next,"run_id":r["receipt"]["run_id"]}),
    )
    .await;
    assert_eq!(replay["products_match"], true);
    assert_eq!(replay["status_matches"], true);
    assert_eq!(
        artifact(&c, &replay["artifact_id"]).await["current_report"]["status"],
        "ConsistentUnderAssumptions"
    );
    let again = call(&c, "record_workspace_run", args).await;
    assert_eq!(again["workspace"]["id"], r["workspace"]["id"]);
    reject(
        &c,
        "replay_workspace_run",
        json!({"workspace":wr,"run_id":r["receipt"]["run_id"]}),
    )
    .await;
    let _ = c.cancel().await;
    server.await.unwrap();
}
