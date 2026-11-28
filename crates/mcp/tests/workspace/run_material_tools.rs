//! Explicit exact-workspace origin discovery and documentary capsule preparation.
#[path = "../../../core/tests/common/material_library.rs"]
mod material_library;
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
async fn exact_native_origin_discovery_preparation_and_immutable_context_refusal() {
    let (client, server) = setup().await;
    let tools = client.list_tools(Default::default()).await.unwrap().tools;
    assert_eq!(tools.len(), collagene_mcp::CURRENT_TOOL_COUNT);
    assert!(tools.iter().any(|t| t.name == "list_run_material_origins"));
    assert!(tools.iter().any(|t| t.name == "prepare_run_material"));
    let commit = original();
    let selected = Selection {
        run_id: commit.receipt.run.run_id.clone(),
        step_index: 0,
        event: commit.receipt.construct.clone(),
    };
    let imported=call(&client,"import_workspace_archive",json!({"data_base64":STANDARD.encode(workspace::save_workspace(&commit.workspace).unwrap())})).await;
    let wr = json!({"id":imported["workspace"]["id"],"revision":imported["workspace"]["revision"]});
    let count = call(&client, "list_documents", json!({})).await["total"].clone();
    let args = json!({"workspace":wr,"event":selected.event});
    let origins = call(&client, "list_run_material_origins", args.clone()).await;
    assert_eq!(origins["selected"], false);
    assert_eq!(origins["origins"], json!([selected.clone()]));
    let discovery = evidence(&client, &origins["artifact_id"]).await;
    assert_eq!(discovery["selected"], Value::Null);
    let prepared = call(
        &client,
        "prepare_run_material",
        json!({"workspace":wr,"selection":selected}),
    )
    .await;
    assert_eq!(prepared["capsule_available"], true);
    let full = evidence(&client, &prepared["artifact_id"]).await;
    let inspection: run_material::Inspection = serde_json::from_value(full).unwrap();
    inspection
        .validate_binding(prepared["context_revision"].as_str().unwrap())
        .unwrap();
    let capsule = inspection.selected.unwrap();
    capsule.validate_documentary().unwrap();
    assert_eq!(capsule.selection, selected);
    assert_eq!(capsule.inputs.len(), 1);
    assert_eq!(capsule.inputs[0].event, event(&commit.workspace, 0));
    let repeat = call(
        &client,
        "prepare_run_material",
        json!({"workspace":wr,"selection":selected}),
    )
    .await;
    assert_eq!(repeat["artifact_id"], prepared["artifact_id"]);
    // Same original reference letters do not give the original-input event a run origin.
    let unproduced = call(
        &client,
        "list_run_material_origins",
        json!({"workspace":wr,"event":event(&commit.workspace,0)}),
    )
    .await;
    assert_eq!(unproduced["origins"], json!([]));
    // The saved document's physical marker must stay rejected by the plain
    // inspector; only this explicit origin branch reconstructs actual material.
    let catalog = collagene_core::cloning::restriction::catalog::builtin_catalog().unwrap();
    let id = catalog
        .edition()
        .enzymes
        .iter()
        .find(|e| e.name == "EcoRI")
        .unwrap()
        .id
        .clone();
    let request = json!({"template_strand":"Top","source_role":"savedPCR","forward":Primer::new("F","ACGACCCC","ATGCGAATTC").unwrap(),"reverse":Primer::new("R","GTTACCCC","ATGCGAATTC").unwrap(),"selection":collagene_core::cloning::restriction::catalog::select(&catalog,&[id]).unwrap(),"options":collagene_core::cloning::restriction::inspection::InspectionOptions::default()});
    let mut inspect_args = json!({"document":{"id":selected.event.revision,"revision":selected.event.revision},"request":request});
    reject(&client, "inspect_pcr_digest", inspect_args.clone()).await;
    inspect_args["origin"] = serde_json::to_value(
        collagene_core::material::run::MaterialOrigin::Run(capsule.clone()),
    )
    .unwrap();
    reject(&client, "inspect_pcr_digest", inspect_args.clone()).await;
    inspect_args["workspace"] = wr.clone();
    let result = call(&client, "inspect_pcr_digest", inspect_args.clone()).await;
    assert_eq!(result["fragment_count"], 3);
    assert_eq!(result["amplicon_length"], 36);
    assert_eq!(result["products_determined"], true);
    assert_eq!(result["origin_revision"].as_str().unwrap().len(), 64);
    let wrapper: collagene_core::cloning::pcr::digest::MaterialPreview =
        serde_json::from_value(evidence(&client, &result["artifact_id"]).await).unwrap();
    wrapper
        .validate_binding(result["context_revision"].as_str().unwrap())
        .unwrap();
    assert_eq!(
        wrapper.context_revision,
        result["context_revision"].as_str().unwrap()
    );
    assert_eq!(
        wrapper.origin_revision,
        result["origin_revision"].as_str().unwrap()
    );
    assert_eq!(
        wrapper.preview.amplicon.as_ref().unwrap().sequence.as_str(),
        "ATGCGAATTCACGACCCCGGGGTAACGAATTCGCAT"
    );
    let fragments = &wrapper.preview.inspection.as_ref().unwrap().fragments;
    assert_eq!(fragments.len(), 3);
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
    let mut changed = inspect_args;
    changed["origin"]["Run"]["selection"]["event"]["document_id"] = json!("forged-owner");
    reject(&client, "inspect_pcr_digest", changed).await;
    let mut stale = args;
    stale["workspace"]["revision"] = json!("0".repeat(64));
    reject(&client, "list_run_material_origins", stale).await;
    let mut forged = selected.clone();
    forged.event.document_id = "another-owner".into();
    reject(
        &client,
        "prepare_run_material",
        json!({"workspace":wr,"selection":forged}),
    )
    .await;
    let mut bad_step = selected;
    bad_step.step_index = 999;
    reject(
        &client,
        "prepare_run_material",
        json!({"workspace":wr,"selection":bad_step}),
    )
    .await;
    assert_eq!(
        call(&client, "list_documents", json!({})).await["total"],
        count
    );
    drop(client);
    server.abort();
}

#[tokio::test]
async fn derived_bottom_pcr_digest_keeps_native_extra_calls_through_sdk() {
    use collagene_core::{
        chemistry::Strand,
        cloning::restriction::catalog as catalog_digest,
        material::{derivation as derived, run::MaterialOrigin},
    };
    let (workspace, event, origin) =
        material_library::native_source("AGTCCGATGCTAC", "GTAGGTAGCATCGGACT", 13);
    let catalog = catalog_digest::builtin_catalog().unwrap();
    let selection = |name: &str| {
        catalog_digest::select(
            &catalog,
            &[catalog
                .edition()
                .enzymes
                .iter()
                .find(|e| e.name == name)
                .unwrap()
                .id
                .clone()],
        )
        .unwrap()
    };
    let prepared = derived::prepare_install(
        &workspace,
        &derived::Request {
            sources: vec![collagene_core::material::derivation::BoundSource {
                role: collagene_core::material::derivation::SourceRole::Template,
                event: event,
                origin: derived::SourceOrigin::Saved(origin),
            }],
            derivation: derived::Derivation::CatalogDigest {
                selection: selection("SmaI"),
                options: Default::default(),
                fragment_index: 0,
            },
            name: "Saved derived unequal template".into(),
            save_primers: false,
        },
    )
    .unwrap();
    let workspace =
        workspace::load_workspace(&workspace::save_workspace(&prepared.workspace).unwrap())
            .unwrap();
    let origin =
        MaterialOrigin::Derived(derived::capsule(&workspace, &prepared.receipt.product).unwrap());
    let (client, server) = setup().await;
    let uploaded = call(
        &client,
        "import_workspace_archive",
        json!({"data_base64":STANDARD.encode(workspace::save_workspace(&workspace).unwrap())}),
    )
    .await;
    let reference =
        json!({"id":uploaded["workspace"]["id"],"revision":uploaded["workspace"]["revision"]});
    let request = collagene_core::cloning::pcr::digest::Request {
        source_role: "native-bottom".into(),
        template_strand: Strand::Bottom,
        forward: Primer::new("Bottom F", "GTAGGTAG", "CCGAATTC").unwrap(),
        reverse: Primer::new("Bottom R", "AGTCCGAT", "GGGAATTC").unwrap(),
        selection: selection("EcoRI"),
        options: Default::default(),
    };
    let args = json!({"document":{"id":origin.event().revision,"revision":origin.event().revision},"origin":origin,"workspace":reference,"request":request});
    let result = call(&client, "inspect_pcr_digest", args.clone()).await;
    let report: collagene_core::cloning::pcr::digest::MaterialPreview =
        serde_json::from_value(evidence(&client, &result["artifact_id"]).await).unwrap();
    assert_eq!(
        report.preview.amplicon.as_ref().unwrap().sequence.as_str(),
        "CCGAATTCGTAGGTAGCATCGGACTGAATTCCC"
    );
    assert_eq!(
        report.preview.pcr_geometries[0].template_strand,
        Strand::Bottom
    );
    assert_eq!(result["origin_incomplete"], report.origin_incomplete);
    assert_eq!(
        result["origin_findings_total"],
        report.origin_findings_total
    );
    assert_eq!(result["fragment_count"], 3);
    let mut stale = args;
    stale["origin"]["Derived"]["product"]["event_id"] = json!("foreign-event");
    reject(&client, "inspect_pcr_digest", stale).await;
    client.cancel().await.unwrap();
    server.await.unwrap();
}
