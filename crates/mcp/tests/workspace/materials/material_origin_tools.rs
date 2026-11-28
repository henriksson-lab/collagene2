//! Explicit exact-workspace origin discovery and documentary capsule preparation.
use base64::{Engine, engine::general_purpose::STANDARD};
use collagene_core::{
    documents::workspace::{self, Workspace},
    material::{
        lineage::{EventBinding, InputEventBinding},
        outputs::{self as cloning_outputs, Commit},
        run::{self as run_material, Selection},
    },
    planning::plan::{DesignPlan, Operation, PlanInput, Source},
    sequence::analysis::Primer,
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
    let source = Molecule::dna("Original", "ACGACCCCGGGGTAAC", Topology::Linear).unwrap();
    let target = Molecule::dna(
        "Independent requirement",
        "ACGACCCCGGGGTAAC",
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
            forward: Primer::new("F", "ACGACCCC", "").unwrap(),
            reverse: Primer::new("R", "GTTACCCC", "").unwrap(),
            template_strand: collagene_core::chemistry::Strand::Top,
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
async fn actual_saved_pcr_digest_discovery_and_explicit_replay_are_exact_and_read_only() {
    let (client, server) = setup().await;
    let tools = client.list_tools(Default::default()).await.unwrap().tools;
    assert_eq!(tools.len(), collagene_mcp::CURRENT_TOOL_COUNT);
    let tool = tools
        .iter()
        .find(|tool| tool.name == "inspect_material_origin")
        .unwrap();
    assert_eq!(
        tool.input_schema.get("additionalProperties"),
        Some(&json!(false))
    );
    let mut commit = original();
    let source = commit.receipt.construct.clone();
    commit.workspace.active_document = source.document_id.clone();
    let origin = collagene_core::material::capsule::capsule(
        &commit.workspace,
        &Selection {
            run_id: commit.receipt.run.run_id,
            step_index: 0,
            event: source.clone(),
        },
    )
    .unwrap();
    let catalog = collagene_core::cloning::restriction::catalog::builtin_catalog().unwrap();
    let enzyme = catalog
        .edition()
        .enzymes
        .iter()
        .find(|enzyme| enzyme.name == "SmaI")
        .unwrap();
    let selection =
        collagene_core::cloning::restriction::catalog::select(&catalog, &[enzyme.id.clone()])
            .unwrap();
    let installed = collagene_core::material::derivation::prepare_install(
        &commit.workspace,
        &collagene_core::material::derivation::Request {
            sources: vec![collagene_core::material::derivation::BoundSource {
                role: collagene_core::material::derivation::SourceRole::Template,
                event: source,
                origin: collagene_core::material::derivation::SourceOrigin::Saved(
                    run_material::MaterialOrigin::Run(origin),
                ),
            }],
            derivation: collagene_core::material::derivation::Derivation::CatalogDigest {
                selection,
                options: Default::default(),
                fragment_index: 0,
            },
            name: "Selected saved-PCR SmaI fragment".into(),
            save_primers: false,
        },
    )
    .unwrap();
    let product = installed.receipt.product;
    assert_eq!(
        installed.workspace.revisions[&product.revision]
            .sequence
            .as_str(),
        "ACGACCCC"
    );
    let bytes = workspace::save_workspace(&installed.workspace).unwrap();
    let before = serde_json::to_vec(&installed.workspace).unwrap();
    let imported = call(
        &client,
        "import_workspace_archive",
        json!({"data_base64":STANDARD.encode(bytes)}),
    )
    .await;
    let wr = json!({"id":imported["workspace"]["id"],"revision":imported["workspace"]["revision"]});
    let count = call(&client, "list_documents", json!({})).await["total"].clone();
    let selected = run_material::MaterialSelection::Derived(product.clone());
    let discovery = call(
        &client,
        "inspect_material_origin",
        json!({"workspace":wr,"event":product}),
    )
    .await;
    assert_eq!(discovery["origin_available"], false);
    assert_eq!(discovery["origins"], json!([selected.clone()]));
    let reply = call(
        &client,
        "inspect_material_origin",
        json!({"workspace":wr,"event":product,"selection":selected}),
    )
    .await;
    assert_eq!(reply["origin_available"], true);
    let full = evidence(&client, &reply["artifact_id"]).await;
    let inspection: run_material::MaterialInspection = serde_json::from_value(full).unwrap();
    inspection
        .validate_binding(reply["context_revision"].as_str().unwrap())
        .unwrap();
    let origin = inspection.selected.unwrap();
    assert!(matches!(&origin, run_material::MaterialOrigin::Derived(_)));
    origin
        .validate_for_workspace(&installed.workspace, &product)
        .unwrap();
    let run_material::MaterialOrigin::Derived(capsule) = origin else {
        unreachable!()
    };
    let reference = &installed.workspace.revisions[&product.revision];
    let replay = collagene_core::material::derivation::replay(&capsule, reference).unwrap();
    assert_eq!(replay.tracked().duplex.top.as_str(), "ACGACCCC");
    let wrong = event(&installed.workspace, 0);
    reject(&client,"inspect_material_origin",json!({"workspace":wr,"event":product,"selection":run_material::MaterialSelection::Derived(wrong)})).await;
    reject(
        &client,
        "inspect_material_origin",
        json!({"workspace":wr,"event":product,"selection":selected,"claimed_duplex":{}}),
    )
    .await;
    let mut stale = wr.clone();
    stale["revision"] = json!("0".repeat(64));
    reject(
        &client,
        "inspect_material_origin",
        json!({"workspace":stale,"event":product,"selection":selected}),
    )
    .await;
    assert_eq!(
        call(&client, "list_documents", json!({})).await["total"],
        count
    );
    assert_eq!(serde_json::to_vec(&installed.workspace).unwrap(), before);
    client.cancel().await.unwrap();
    server.await.unwrap();
}
