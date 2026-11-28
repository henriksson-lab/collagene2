//! Direct saved-source catalog inspection keeps exact workspace proof and typed chemistry.
use base64::{Engine, engine::general_purpose::STANDARD};
use collagene_core::{
    documents::workspace::{self, Workspace},
    material::{
        lineage::{EventBinding, InputEventBinding},
        outputs::{self as cloning_outputs, Commit},
        run::Selection,
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
async fn predicted_catalog_gel_replays_run_and_derived_sources_with_exact_workspace_ancestry() {
    use collagene_core::{
        cloning::restriction::{catalog, inspection::InspectionOptions, setups},
        material::{
            capsule, derivation,
            run::{MaterialOrigin, MaterialSelection},
        },
        output::gel::{Calibration, Preparation},
    };
    let commit = original();
    let selected = Selection {
        run_id: commit.receipt.run.run_id.clone(),
        step_index: 0,
        event: commit.receipt.construct.clone(),
    };
    let origin = capsule::capsule(&commit.workspace, &selected).unwrap();
    let embedded = catalog::builtin_catalog().unwrap();
    let enzymes = catalog::select(&embedded, &["ecori".into()]).unwrap();
    let mut frame = commit.workspace;
    frame.active_document = selected.event.document_id.clone();
    let installed = derivation::install(
        &frame,
        &derivation::Request {
            sources: vec![derivation::BoundSource {
                role: derivation::SourceRole::Template,
                event: selected.event.clone(),
                origin: derivation::SourceOrigin::Saved(MaterialOrigin::Run(origin.clone())),
            }],
            derivation: derivation::Derivation::CatalogDigest {
                selection: enzymes.clone(),
                options: InspectionOptions::default(),
                fragment_index: 0,
            },
            name: "Saved selected digest fragment".into(),
            save_primers: false,
        },
    )
    .unwrap();
    let frame = workspace::load_workspace(&installed.archive).unwrap();
    let selections = [
        MaterialSelection::Run(selected),
        MaterialSelection::Derived(installed.receipt.product),
    ];
    let (client, server) = setup().await;
    let imported = call(
        &client,
        "import_workspace_archive",
        json!({"data_base64":STANDARD.encode(workspace::save_workspace(&frame).unwrap())}),
    )
    .await;
    let workspace_ref =
        json!({"id":imported["workspace"]["id"],"revision":imported["workspace"]["revision"]});
    let before = call(&client, "list_documents", json!({})).await;
    for selection in selections {
        let mut origins = capsule::documentary_origins(&frame, &[selection.clone()], &[]).unwrap();
        let recipe = setups::Recipe {
            source: selection.event().clone(),
            origin: derivation::SourceOrigin::Saved(origins.materials.remove(0)),
            selection: enzymes.clone(),
            options: InspectionOptions::default(),
            cleavage_mode: setups::CleavageMode::AssumedCompleteMethylationUnassessed,
        };
        let query = json!({"lanes":[{"label":"Fresh saved-material prediction",
            "document":{"id":recipe.source.revision,"revision":recipe.source.revision},"workspace":workspace_ref,
            "preparation":Preparation::CatalogDigest{recipe:recipe.clone()},"mass_ng":null}],"calibration":Calibration::default()});
        let predicted = call(&client, "predict_gel", query.clone()).await;
        assert!(predicted["lanes"][0]["error"].is_null(), "{predicted}");
        assert!(predicted["lanes"][0]["fragment_count"].as_u64().unwrap() > 0);
        let source = frame.revisions[&recipe.source.revision].clone();
        let unrelated =
            Workspace::from_documents("No native material owner", vec![source], 0).unwrap();
        let unrelated = call(
            &client,
            "import_workspace_archive",
            json!({"data_base64":STANDARD.encode(workspace::save_workspace(&unrelated).unwrap())}),
        )
        .await;
        let mut missing_owner = query.clone();
        missing_owner["lanes"][0]["workspace"] = json!({"id":unrelated["workspace"]["id"],"revision":unrelated["workspace"]["revision"]});
        reject(&client, "predict_gel", missing_owner).await;
        let mut plain = query;
        plain["lanes"][0]["preparation"]["CatalogDigest"]["recipe"]["origin"] =
            json!("OriginalReference");
        let withheld = call(&client, "predict_gel", plain).await;
        assert_eq!(withheld["lanes"][0]["fragment_count"], 0);
        assert_eq!(withheld["lanes"][0]["complete_under_model"], false);
    }
    assert_eq!(call(&client, "list_documents", json!({})).await, before);
    client.cancel().await.unwrap();
    server.await.unwrap();
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
    let origin = collagene_core::material::capsule::capsule(&commit.workspace, &selected).unwrap();
    let imported=call(&client,"import_workspace_archive",json!({"data_base64":STANDARD.encode(workspace::save_workspace(&commit.workspace).unwrap())})).await;
    let wr = json!({"id":imported["workspace"]["id"],"revision":imported["workspace"]["revision"]});
    let before = call(&client, "list_documents", json!({})).await["total"].clone();
    let mut args = json!({"document":{"id":selected.event.revision,"revision":selected.event.revision},"source_role":"saved","enzymes":["EcoRI"],"max_sites":10,"max_fragments":1});
    let plain = call(&client, "inspect_restriction", args.clone()).await;
    let full: collagene_mcp::tools::cloning::RestrictionInspectionEvidence =
        serde_json::from_value(evidence(&client, &plain["artifact_id"]).await).unwrap();
    let collagene_mcp::tools::cloning::RestrictionInspectionEvidence::Original { inspection } =
        full
    else {
        panic!("unattested source inspection must remain explicitly original-reference evidence");
    };
    let recognition = inspection.inspection;
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
    let full: collagene_mcp::tools::cloning::RestrictionInspectionEvidence =
        serde_json::from_value(evidence(&client, &result["artifact_id"]).await).unwrap();
    let collagene_mcp::tools::cloning::RestrictionInspectionEvidence::SavedRun {
        inspection: report,
    } = full
    else {
        panic!("saved source inspection must retain full replay evidence");
    };
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
    assert!(result.get("catalog_artifact_id").is_none());
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
