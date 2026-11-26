//! Session-scoped saved-template TOPO inspection never restores private material
//! from detached metadata or changes the independently selected target.
use base64::{Engine, engine::general_purpose::STANDARD};
use collagene_core::{
    Molecule, Topology,
    analysis::Primer,
    cloning_outputs::{self, Commit},
    run_lineage::{EventBinding, InputEventBinding},
    run_material,
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
async fn call(client: &Client, name: &str, args: Value) -> Value {
    let result = client
        .call_tool(
            CallToolRequestParams::new(name.to_string())
                .with_arguments(args.as_object().unwrap().clone()),
        )
        .await
        .unwrap();
    assert_ne!(result.is_error, Some(true), "{name}: {:?}", result.content);
    result.structured_content.unwrap()
}
async fn reject(client: &Client, args: Value) {
    let result = client
        .call_tool(
            CallToolRequestParams::new("inspect_saved_topo")
                .with_arguments(args.as_object().unwrap().clone()),
        )
        .await;
    assert!(result.is_err() || result.unwrap().is_error == Some(true));
}
async fn artifact(client: &Client, id: &Value) -> Value {
    let mut text = String::new();
    let mut start = 0;
    loop {
        let page = call(
            client,
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
fn event(workspace: &Workspace, index: usize) -> EventBinding {
    let doc = &workspace.documents[index];
    EventBinding {
        document_id: doc.id.clone(),
        event_id: doc.head_event.clone().unwrap(),
        revision: doc.head.clone(),
    }
}
fn saved() -> Commit {
    let source =
        Molecule::dna("Original literal template", "CACCATGACG", Topology::Linear).unwrap();
    let target = Molecule::dna("Independent first PCR", "CACCATGACG", Topology::Linear).unwrap();
    let plan = DesignPlan {
        name: "Recorded actual PCR".into(),
        target_revision: target.revision(),
        inputs: vec![PlanInput {
            id: "original".into(),
            revision: source.revision(),
            molecule: source.clone(),
        }],
        steps: vec![Operation::Pcr {
            template: Source::Input("original".into()),
            forward: Primer::new("F", "CACCA", "").unwrap(),
            reverse: Primer::new("R", "CGTCA", "").unwrap(),
        }],
    };
    let workspace = Workspace::from_documents("Saved PCR source", vec![source, target], 1).unwrap();
    cloning_outputs::finalize(
        &workspace,
        &plan,
        event(&workspace, 1),
        vec![InputEventBinding {
            input_id: "original".into(),
            event: event(&workspace, 0),
        }],
        &cloning_outputs::Options {
            name: "Saved insert PCR".into(),
            show_prepared: false,
            save_primers: false,
            allow_incomplete: false,
        },
    )
    .unwrap()
}
fn declared_strand(role: &str, length: usize) -> collagene_core::topo::PreparedStrand {
    use collagene_core::topo::{
        self as t,
        physical::{Enzyme, FivePrime, ThreePrime},
    };
    t::PreparedStrand {
        id: role.into(),
        segments: vec![t::SourceInterval {
            source_role: role.into(),
            start: 0,
            length,
            reverse: false,
        }],
        five_prime: FivePrime::Phosphate,
        three_prime: ThreePrime::Phosphotyrosyl {
            enzyme: Enzyme::VacciniaTypeIBTyr274,
            activation_id: role.into(),
        },
    }
}
async fn import(client: &Client, name: &str, dna: &str) -> Value {
    let molecule = Molecule::dna(name, dna, Topology::Linear).unwrap();
    let reply = call(
        client,
        "import_document",
        json!({"text":collagene_core::formats::export_genbank(&molecule).unwrap()}),
    )
    .await;
    json!({"id":reply["documents"][0]["id"],"revision":reply["documents"][0]["revision"]})
}
#[tokio::test]
async fn actual_saved_template_inspection_attests_origin_and_preserves_independent_check_target() {
    use collagene_core::{topo as t, topo_preparation as prep, topo_run_preparation as saved};
    let (client, server) = setup().await;
    assert!(
        client
            .list_tools(Default::default())
            .await
            .unwrap()
            .tools
            .iter()
            .any(|tool| tool.name == "inspect_saved_topo")
    );
    let commit = saved();
    let selection = run_material::Selection {
        run_id: commit.receipt.run.run_id.clone(),
        step_index: 0,
        event: commit.receipt.construct.clone(),
    };
    let imported=call(&client,"import_workspace_archive",json!({"data_base64":STANDARD.encode(workspace::save_workspace(&commit.workspace).unwrap())})).await;
    let workspace =
        json!({"id":imported["workspace"]["id"],"revision":imported["workspace"]["revision"]});
    let prepared = call(
        &client,
        "prepare_run_material",
        json!({"workspace":workspace,"selection":selection}),
    )
    .await;
    let inspection: run_material::Inspection =
        serde_json::from_value(artifact(&client, &prepared["artifact_id"]).await).unwrap();
    let origin = inspection.selected.unwrap();
    let vt = import(&client, "Declared vector top", "AAGGGGCATCCCTT").await;
    let vb = import(&client, "Declared vector bottom", "GGTGAAGGGATGCCCCTT").await;
    let template = json!({"id":selection.event.revision,"revision":selection.event.revision});
    let declaration = prep::PreparationRequest {
        mode: prep::PreparationMode::DirectionalCacc,
        vector: prep::PreparedVectorDeclaration {
            declaration_id: "Explicit vector geometry".into(),
            assumption: "Synthetic activation hypothesis, not an observed stock or kit claim"
                .into(),
            strands: vec![declared_strand("vt", 14), declared_strand("vb", 18)],
            pairing: vec![t::physical::PairSpan {
                strand_a: 0,
                start_a: 0,
                strand_b: 1,
                start_b: 4,
                length: 14,
            }],
        },
        insert: prep::ExactPcrDeclaration {
            material_id: "New PCR from saved template".into(),
            template_role: "template".into(),
            forward: Primer::new("F", "CACCA", "").unwrap(),
            reverse: Primer::new("R", "CGTCA", "").unwrap(),
            tail: None,
        },
        engine_limits: Default::default(),
        projection_options: Default::default(),
    };
    let compiled = prep::compile(&declaration).unwrap();
    let mut args = json!({"inputs":[{"id":"vt","document":vt},{"id":"vb","document":vb},{"id":"template","document":template}],"template_role":"template","origin":origin,"origin_workspace":workspace,
        "task":{"Orientations":{"request":compiled.request,"vector":compiled.vector_strands,"insert":compiled.insert_strands,"directional":compiled.directional}}});
    let before = call(&client, "list_documents", json!({})).await["total"].clone();
    let reply = call(&client, "inspect_saved_topo", args.clone()).await;
    assert_eq!(reply["candidate_saved"], false);
    assert_eq!(reply["document_saved"], false);
    assert_eq!(reply["verified"], false);
    let full = artifact(&client, &reply["artifact_id"]).await;
    let inspected: saved::Inspection = serde_json::from_value(full["inspection"].clone()).unwrap();
    inspected
        .validate_binding(reply["context_revision"].as_str().unwrap())
        .unwrap();
    assert_eq!(
        reply["origin_findings_total"],
        json!(inspected.origin_findings_total)
    );
    assert_eq!(
        reply["origin_findings_truncated"],
        json!(inspected.origin_findings_truncated)
    );
    assert!(!inspected.origin_incomplete);
    assert_eq!(reply["origin_incomplete"], false);
    assert_eq!(
        reply["origin_findings"],
        serde_json::to_value(&inspected.origin_findings).unwrap()
    );
    let saved::Reply::Orientations(oriented) = inspected.reply else {
        panic!("Expected actual saved-template orientations");
    };
    assert_eq!(oriented.reference.candidates.len(), 2);
    assert!(!oriented.reference.biochemical_mixture_exhaustive);
    assert_eq!(
        call(&client, "list_documents", json!({})).await["total"],
        before
    );
    // A retained archive with the same input letters but no selected product/run
    // cannot attest this origin, even though the capsule remains fully valid.
    let wrong = Workspace::from_documents(
        "Unrelated source-only archive",
        vec![commit.workspace.revisions[&selection.event.revision].clone()],
        0,
    )
    .unwrap();
    let wrong = call(
        &client,
        "import_workspace_archive",
        json!({"data_base64":STANDARD.encode(workspace::save_workspace(&wrong).unwrap())}),
    )
    .await;
    let mut bad = args.clone();
    bad["origin_workspace"] =
        json!({"id":wrong["workspace"]["id"],"revision":wrong["workspace"]["revision"]});
    reject(&client, bad).await;
    let mut bad = args.clone();
    bad["inputs"][2]["document"]["revision"] = json!("0".repeat(64));
    reject(&client, bad).await;
    let mut bad = args.clone();
    bad["template_role"] = json!("vt");
    reject(&client, bad).await;
    let mut bad = args.clone();
    bad.as_object_mut().unwrap().remove("origin_workspace");
    reject(&client, bad).await;
    // The selected component is an actual declared graph component, never a
    // component chosen because its sequence fits this deliberately wrong target.
    let target = import(
        &client,
        "Independent deliberately different requirement",
        "TTTTCCCCAAAAGGGG",
    )
    .await;
    let batch = blake3::hash(
        &serde_json::to_vec(&(
            t::ENGINE_VERSION,
            t::physical::ENGINE,
            t::physical::PRIMARY_MODEL,
            t::physical::projection::ENGINE,
            &oriented.sources,
            &compiled.request,
        ))
        .unwrap(),
    )
    .to_hex()
    .to_string();
    args["task"] = json!({"Check":{"request":compiled.request,"selection":{"batch_revision":batch,"component":0,"reference_seed":{"strand":0,"index":0}},"target":target,"policy":t::AcceptancePolicy::default()}});
    let checked = call(&client, "inspect_saved_topo", args.clone()).await;
    assert_eq!(checked["summary"]["kind"], "Checked");
    assert_eq!(checked["summary"]["target_revision"], target["revision"]);
    assert_ne!(
        checked["summary"]["outcome"],
        json!("AcceptedUnderDeclaredGraph")
    );
    let checked_full = artifact(&client, &checked["artifact_id"]).await;
    assert_eq!(checked_full["task"]["Check"]["target"], target);
    let saved::Reply::Checked(assessment) =
        serde_json::from_value::<saved::Inspection>(checked_full["inspection"].clone())
            .unwrap()
            .reply
    else {
        panic!("Expected independent target assessment");
    };
    assert_eq!(
        assessment.target_revision,
        target["revision"].as_str().unwrap()
    );
    args["task"]["Check"]["target"]["revision"] = json!("0".repeat(64));
    reject(&client, args.clone()).await;
    args["task"]["Check"]["target"] = json!({"name":"Injected target","sequence":"TTTT"});
    reject(&client, args).await;
    client.cancel().await.unwrap();
    server.await.unwrap();
}
