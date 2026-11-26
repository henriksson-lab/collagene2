//! Direct saved-source catalog inspection keeps exact workspace proof and typed chemistry.
use base64::{Engine, engine::general_purpose::STANDARD};
use collagene_core::{
    Molecule, Topology,
    analysis::Primer,
    cloning_outputs::{self, Commit},
    run_lineage::{EventBinding, InputEventBinding},
    run_material::{self, Selection},
    verification::{DesignPlan, Operation, PlanInput, Source},
    workspace::{self, Workspace},
};
use collagene_mcp::CollageneServer;
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
    let result = c
        .call_tool(
            CallToolRequestParams::new(name.to_owned())
                .with_arguments(args.as_object().unwrap().clone()),
        )
        .await
        .unwrap();
    assert_ne!(result.is_error, Some(true), "{name}: {:?}", result.content);
    result.structured_content.unwrap()
}
async fn reject(c: &Client, name: &str, args: Value) {
    let result = c
        .call_tool(
            CallToolRequestParams::new(name.to_owned())
                .with_arguments(args.as_object().unwrap().clone()),
        )
        .await;
    match result {
        Err(_) => {}
        Ok(reply) => assert_eq!(
            reply.is_error,
            Some(true),
            "{name} unexpectedly succeeded: {:?}",
            reply
        ),
    }
}
async fn evidence(c: &Client, id: &Value) -> Value {
    let mut text = String::new();
    let mut start = 0;
    loop {
        let page = call(
            c,
            "read_artifact",
            json!({"artifact_id":id,"start":start,"max_bytes":16384}),
        )
        .await;
        text.push_str(page["text"].as_str().unwrap());
        if page["done"] == true {
            break;
        }
        start = page["next_start"].as_u64().unwrap();
    }
    serde_json::from_str(&text).unwrap()
}
fn event(w: &Workspace, i: usize) -> EventBinding {
    let d = &w.documents[i];
    EventBinding {
        document_id: d.id.clone(),
        event_id: d.head_event.clone().unwrap(),
        revision: d.head.clone(),
    }
}
fn original() -> Commit {
    let source = Molecule::dna("Original", "ACGAGAATTCCGTAAC", Topology::Linear).unwrap();
    let target = Molecule::dna(
        "Independent requirement",
        "ACGAGAATTCCGTAAC",
        Topology::Linear,
    )
    .unwrap();
    let plan = DesignPlan {
        name: "Actual original PCR".into(),
        target_revision: target.revision(),
        inputs: vec![PlanInput {
            id: "original".into(),
            revision: source.revision(),
            molecule: source.clone(),
        }],
        steps: vec![Operation::Pcr {
            template: Source::Input("original".into()),
            forward: Primer::new("F", "ACGAGAAT", "").unwrap(),
            reverse: Primer::new("R", "GTTACGGA", "").unwrap(),
        }],
    };
    let w = Workspace::from_documents("Recorded PCR", vec![source, target], 1).unwrap();
    cloning_outputs::finalize(
        &w,
        &plan,
        event(&w, 1),
        vec![InputEventBinding {
            input_id: "original".into(),
            event: event(&w, 0),
        }],
        &cloning_outputs::Options {
            name: "Saved PCR".into(),
            show_prepared: false,
            save_primers: false,
            allow_incomplete: false,
        },
    )
    .unwrap()
}
#[tokio::test]
async fn explicit_workspace_saved_catalog_inspection_retains_outer_oh_and_cut_p_without_publishing()
{
    let (client, server) = setup().await;
    assert_eq!(
        client
            .list_tools(Default::default())
            .await
            .unwrap()
            .tools
            .len(),
        collagene_mcp::CURRENT_TOOL_COUNT
    );
    let commit = original();
    let selected = Selection {
        run_id: commit.receipt.run.run_id.clone(),
        step_index: 0,
        event: commit.receipt.construct.clone(),
    };
    let origin = run_material::capsule(&commit.workspace, &selected).unwrap();
    let imported=call(&client,"import_workspace_archive",json!({"data_base64":STANDARD.encode(workspace::save_workspace(&commit.workspace).unwrap())})).await;
    let wr = json!({"id":imported["workspace"]["id"],"revision":imported["workspace"]["revision"]});
    let before = call(&client, "list_documents", json!({})).await["total"].clone();
    let mut args = json!({"document":{"id":selected.event.revision,"revision":selected.event.revision},"source_role":"saved","enzymes":["EcoRI"],"max_sites":10,"max_fragments":1});
    let plain = call(&client, "inspect_restriction", args.clone()).await;
    let recognition: collagene_core::restriction_inspection::RestrictionInspection =
        serde_json::from_value(evidence(&client, &plain["artifact_id"]).await).unwrap();
    assert!(recognition.scan_complete);
    assert_eq!(recognition.sites.len(), 1);
    assert_eq!(recognition.sites[0].recognition_start, 4);
    assert_eq!(recognition.sites[0].top_cut, Some(5));
    assert_eq!(recognition.sites[0].bottom_cut, Some(9));
    assert!(!recognition.products_determined);
    assert!(recognition.fragments.is_empty());
    assert!(
        recognition
            .diagnostics
            .iter()
            .any(|message| message.contains("Physical products withheld"))
    );

    args["origin"] = serde_json::to_value(&origin).unwrap();
    reject(&client, "inspect_restriction", args.clone()).await;
    args["workspace"] = wr.clone();
    let result = call(&client, "inspect_restriction", args.clone()).await;
    assert_eq!(result["fragment_count"], 2);
    assert_eq!(result["products_determined"], true);
    assert_eq!(result["fragments"].as_array().unwrap().len(), 1);
    assert_eq!(result["fragments_done"], false);
    let report: collagene_core::run_catalog_digest::Inspection =
        serde_json::from_value(evidence(&client, &result["artifact_id"]).await).unwrap();
    report
        .validate_binding(result["context_revision"].as_str().unwrap())
        .unwrap();
    assert_eq!(
        report.origin_revision,
        result["origin_revision"].as_str().unwrap()
    );
    assert_eq!(report.fragment_modifications.len(), 2);
    let fragments = &report.inspection.inspection.fragments;
    assert_eq!(fragments[0].duplex.top_length, 5);
    assert_eq!(fragments[1].duplex.top_length, 11);
    assert_eq!(
        fragments[0]
            .duplex
            .ends
            .as_ref()
            .unwrap()
            .0
            .chemistry
            .five_prime,
        collagene_core::chemistry::FivePrime::Hydroxyl
    );
    assert_eq!(
        fragments[1]
            .duplex
            .ends
            .as_ref()
            .unwrap()
            .0
            .chemistry
            .five_prime,
        collagene_core::chemistry::FivePrime::Phosphate
    );
    assert_eq!(result["artifact_id"], result["catalog_artifact_id"]);
    let mut next = args.clone();
    next["fragment_start"] = result["next_fragment_start"].clone();
    let page = call(&client, "inspect_restriction", next).await;
    assert_eq!(page["artifact_id"], result["artifact_id"]);
    assert_eq!(page["fragments_done"], true);
    let mut forged = args.clone();
    forged["origin"]["selection"]["event"]["document_id"] = json!("forged-owner");
    reject(&client, "inspect_restriction", forged).await;
    let mut stale = args.clone();
    stale["workspace"]["revision"] = json!("0".repeat(64));
    reject(&client, "inspect_restriction", stale).await;
    let mut changed = args;
    changed["enzymes"] = json!(["MissingCatalogEnzyme"]);
    reject(&client, "inspect_restriction", changed).await;
    assert_eq!(
        call(&client, "list_documents", json!({})).await["total"],
        before
    );
    drop(client);
    server.abort();
}
