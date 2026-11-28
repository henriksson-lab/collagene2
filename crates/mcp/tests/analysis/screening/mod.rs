//! Current F20 SDK fixtures for the registered screening adapter.
use base64::{Engine, engine::general_purpose::STANDARD};
use collagene_core::{
    chemistry::{
        Strand,
        methylation::MethylationState,
        tracked_duplex::{GenerationAssumptions, checked_selection},
    },
    cloning::screening,
    documents::workspace::{self, Workspace},
    material::{
        capsule,
        derivation::SourceOrigin,
        lineage::{EventBinding, InputEventBinding},
        outputs,
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
async fn call(client: &Client, name: &str, args: Value) -> Value {
    let result = client
        .call_tool(
            CallToolRequestParams::new(name.to_owned())
                .with_arguments(args.as_object().unwrap().clone()),
        )
        .await
        .unwrap();
    assert_ne!(result.is_error, Some(true), "{name}: {:?}", result.content);
    result.structured_content.unwrap()
}
async fn reject(client: &Client, args: Value) {
    match client
        .call_tool(
            CallToolRequestParams::new("assess_screening_strategy")
                .with_arguments(args.as_object().unwrap().clone()),
        )
        .await
    {
        Ok(reply) => assert_eq!(reply.is_error, Some(true)),
        Err(rmcp::service::ServiceError::McpError(error)) => assert_eq!(error.code.0, -32602),
        Err(error) => panic!("unexpected transport error: {error}"),
    }
}
async fn artifact(client: &Client, id: &Value) -> screening::Report {
    let (mut start, mut text) = (0, String::new());
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
fn generation() -> GenerationAssumptions {
    GenerationAssumptions {
        host: None,
        dam: MethylationState::Unmodified,
        dcm: MethylationState::Unmodified,
        cpg: MethylationState::Unmodified,
    }
}
fn subjects() -> Vec<screening::Subject> {
    vec![
        screening::Subject {
            id: "intended".into(),
            label: "Declared intended product".into(),
            role: screening::Role::Intended,
            components: vec!["intended".into()],
        },
        screening::Subject {
            id: "carryover".into(),
            label: "Explicit intended plus parent mixture".into(),
            role: screening::Role::Decoy(screening::DecoyKind::ParentCarryover),
            components: vec!["intended".into(), "parent".into()],
        },
    ]
}
fn fixture() -> (Workspace, screening::Request, Value) {
    let workspace = Workspace::from_documents(
        "literal screening SDK",
        vec![
            Molecule::dna("intended", "AAAAGAATTCTTTT", Topology::Linear).unwrap(),
            Molecule::dna("parent", "CCCCCCCCCCCCCCCCCCCC", Topology::Linear).unwrap(),
        ],
        0,
    )
    .unwrap();
    let components = ["intended", "parent"]
        .into_iter()
        .enumerate()
        .map(|(index, id)| {
            let molecule = workspace
                .head(&workspace.documents[index].id)
                .unwrap()
                .clone();
            screening::Component {
                id: id.into(),
                label: id.into(),
                revision: molecule.revision(),
                molecule,
                origin: SourceOrigin::OriginalReference,
                original_generation: Some(generation()),
            }
        })
        .collect();
    let request = screening::Request {
        components,
        subjects: subjects(),
        tests: vec![screening::Test {
            id: "EcoRI".into(),
            label: "EcoRI".into(),
            method: screening::Method::Digest {
                selection: checked_selection("EcoRI").unwrap(),
            },
        }],
        readout: screening::IdealSizeReadout {
            counterpart_tolerance_twice_bp: 0,
        },
        limits: screening::Limits::default(),
    };
    let query = json!({
        "components":[
            {"id":"intended","label":"intended","event":event(&workspace,0),"origin":"OriginalReference","original_generation":generation()},
            {"id":"parent","label":"parent","event":event(&workspace,1),"origin":"OriginalReference","original_generation":generation()}
        ],"subjects":request.subjects,"tests":request.tests,"readout":request.readout,"limits":request.limits,"max_tests":1
    });
    (workspace, request, query)
}
async fn import(client: &Client, workspace: &Workspace) -> Value {
    let imported = call(
        client,
        "import_workspace_archive",
        json!({"data_base64":STANDARD.encode(workspace::save_workspace(workspace).unwrap())}),
    )
    .await;
    let workspace = &imported["workspace"];
    json!({"id": workspace["id"], "revision": workspace["revision"]})
}
#[tokio::test]
async fn literal_mixture_witnesses_full_artifact_paging_and_no_workspace_mutation() {
    let (workspace, mut native, mut args) = fixture();
    let mut second = native.tests[0].clone();
    second.id = "same enzyme second test".into();
    native.tests.push(second);
    args["tests"] = json!(native.tests);
    let (client, server) = setup().await;
    let tools = client.list_all_tools().await.unwrap();
    let schema = tools
        .iter()
        .find(|tool| tool.name == "assess_screening_strategy")
        .unwrap();
    assert_eq!(tools.len(), collagene_mcp::CURRENT_TOOL_COUNT);
    let schema: Value = serde_json::to_value(&schema.input_schema).unwrap();
    let component = schema["$defs"]
        .as_object()
        .unwrap()
        .values()
        .find(|entry| entry["properties"]["original_generation"].is_object())
        .unwrap();
    assert!(
        component["required"]
            .as_array()
            .unwrap()
            .contains(&json!("original_generation"))
    );
    args["workspace"] = import(&client, &workspace).await;
    let before = call(&client, "list_workspaces", json!({})).await;
    let first = call(&client, "assess_screening_strategy", args.clone()).await;
    assert_eq!(first["page_truncated"], true);
    assert_eq!(first["tests"].as_array().unwrap().len(), 1);
    assert_eq!(first["cells"].as_array().unwrap().len(), 2);
    assert_eq!(first["cells"][0]["sizes_twice_bp"], json!([14]));
    assert_eq!(
        first["tests"][0]["signatures"][1]["sizes_twice_bp"],
        json!([14, 40])
    );
    let report = artifact(&client, &first["artifact_id"]).await;
    report.validate_for(&native).unwrap();
    assert_eq!(report.cells[0].native_product_count, 2);
    assert_eq!(
        report.tests[0].comparisons[0].distinction,
        screening::Distinction::Separated(screening::Witness {
            side: screening::WitnessSide::Decoy,
            size_twice_bp: 40
        })
    );
    assert_eq!(report.ranked_test_indices, [0, 1]);
    assert!(report.panel.complete_through_size_three);
    args["test_start"] = json!(1);
    let second = call(&client, "assess_screening_strategy", args).await;
    assert_eq!(second["cells"][0]["test_index"], 1);
    assert_eq!(second["artifact_id"], first["artifact_id"]);
    assert_eq!(call(&client, "list_workspaces", json!({})).await, before);
    client.cancel().await.unwrap();
    server.await.unwrap();
}
#[tokio::test]
async fn stale_component_and_workspace_pins_and_invalid_bounds_refuse_without_writes() {
    let (workspace, _, mut args) = fixture();
    let (client, server) = setup().await;
    args["workspace"] = import(&client, &workspace).await;
    let before = call(&client, "list_workspaces", json!({})).await;
    let mut stale = args.clone();
    stale["workspace"]["revision"] = json!("0".repeat(64));
    reject(&client, stale).await;
    let mut stale = args.clone();
    stale["components"][0]["event"]["revision"] = json!("0".repeat(64));
    reject(&client, stale).await;
    let mut missing = args.clone();
    missing["components"][0]
        .as_object_mut()
        .unwrap()
        .remove("original_generation");
    reject(&client, missing).await;
    let mut huge = args.clone();
    huge["max_tests"] = json!(9);
    reject(&client, huge).await;
    let mut huge = args.clone();
    huge["limits"]["max_pcr_work"] = json!(200_000_001);
    reject(&client, huge).await;
    let mut duplicate = args.clone();
    duplicate["components"][1]["id"] = json!("intended");
    reject(&client, duplicate).await;
    let mut wrong_pin = args;
    wrong_pin["tests"][0]["method"]["Digest"]["selection"]["catalog_content_hash"] =
        json!("0".repeat(64));
    reject(&client, wrong_pin).await;
    assert_eq!(call(&client, "list_workspaces", json!({})).await, before);
    client.cancel().await.unwrap();
    server.await.unwrap();
}
#[tokio::test]
async fn limited_native_search_never_turns_unsearched_cells_into_absence() {
    let (workspace, mut native, mut args) = fixture();
    native.limits.max_cleavage_work = 1;
    args["limits"] = json!(native.limits);
    let (client, server) = setup().await;
    args["workspace"] = import(&client, &workspace).await;
    let result = call(&client, "assess_screening_strategy", args).await;
    let report = artifact(&client, &result["artifact_id"]).await;
    report.validate_for(&native).unwrap();
    assert!(
        report
            .cells
            .iter()
            .any(|cell| cell.state == screening::CellState::Limited)
    );
    assert!(!report.cells[0].enumeration_complete);
    assert_eq!(
        report.tests[0].comparisons[0].distinction,
        screening::Distinction::Unresolved
    );
    assert!(report.panel.separated_decoy_indices.is_empty());
    client.cancel().await.unwrap();
    server.await.unwrap();
}
#[tokio::test]
async fn exact_saved_material_selector_replays_and_marked_reference_cannot_be_plain_original() {
    let source = Molecule::dna("original template", "ACGACCCCGGGGTAAC", Topology::Linear).unwrap();
    let target = Molecule::dna("independent target", "ACGACCCCGGGGTAAC", Topology::Linear).unwrap();
    let workspace = Workspace::from_documents(
        "saved PCR screening",
        vec![
            source.clone(),
            target,
            Molecule::dna(
                "empty vector declaration",
                "CCCCCCCCCCCCCCCCCCCC",
                Topology::Linear,
            )
            .unwrap(),
        ],
        0,
    )
    .unwrap();
    let forward = Primer::new("forward", "ACGACCCC", "").unwrap();
    let reverse = Primer::new("reverse", "GTTACCCC", "").unwrap();
    let plan = DesignPlan {
        standard: None,
        name: "actual saved PCR".into(),
        target_revision: event(&workspace, 1).revision,
        inputs: vec![PlanInput {
            id: "original".into(),
            revision: source.revision(),
            molecule: source,
        }],
        steps: vec![Operation::Pcr {
            template: Source::Input("original".into()),
            forward: forward.clone(),
            reverse: reverse.clone(),
            template_strand: Strand::Top,
        }],
    };
    let commit = outputs::finalize(
        &workspace,
        &plan,
        event(&workspace, 1),
        vec![InputEventBinding {
            input_id: "original".into(),
            event: event(&workspace, 0),
        }],
        &outputs::Options {
            name: "saved product".into(),
            save_primers: false,
            show_prepared: false,
            allow_incomplete: true,
        },
    )
    .unwrap();
    let selection = Selection {
        run_id: commit.receipt.run.run_id,
        step_index: 0,
        event: commit.receipt.construct,
    };
    let saved = commit.workspace.revisions[&selection.event.revision].clone();
    let parent = event(&workspace, 2);
    let native = screening::Request {
        components: vec![
            screening::Component {
                id: "intended".into(),
                label: "saved PCR".into(),
                revision: saved.revision(),
                molecule: saved,
                origin: SourceOrigin::Saved(collagene_core::material::run::MaterialOrigin::Run(
                    capsule::capsule(&commit.workspace, &selection).unwrap(),
                )),
                original_generation: None,
            },
            screening::Component {
                id: "parent".into(),
                label: "explicit empty vector".into(),
                revision: parent.revision.clone(),
                molecule: workspace.revisions[&parent.revision].clone(),
                origin: SourceOrigin::OriginalReference,
                original_generation: None,
            },
        ],
        subjects: vec![
            subjects()[0].clone(),
            screening::Subject {
                id: "empty".into(),
                label: "Named supplied empty-vector DNA".into(),
                role: screening::Role::Decoy(screening::DecoyKind::EmptyVector),
                components: vec!["parent".into()],
            },
        ],
        tests: vec![screening::Test {
            id: "PCR".into(),
            label: "exact screening PCR".into(),
            method: screening::Method::ExactPcr {
                forward,
                reverse,
                template_strand: Strand::Top,
            },
        }],
        readout: screening::IdealSizeReadout {
            counterpart_tolerance_twice_bp: 0,
        },
        limits: screening::Limits::default(),
    };
    let (client, server) = setup().await;
    let mut args = json!({"workspace":import(&client,&commit.workspace).await,"components":[
        {"id":"intended","label":"saved PCR","event":selection.event,"origin":{"Saved":{"Run":selection}},"original_generation":null},
        {"id":"parent","label":"explicit empty vector","event":parent,"origin":"OriginalReference","original_generation":null}],
        "subjects":native.subjects,"tests":native.tests,"readout":native.readout,"limits":native.limits});
    let result = call(&client, "assess_screening_strategy", args.clone()).await;
    let report = artifact(&client, &result["artifact_id"]).await;
    report.validate_for(&native).unwrap();
    assert_eq!(
        report.components[0].origin,
        screening::ComponentOrigin::SavedNative
    );
    assert_eq!(report.cells[0].sizes_twice_bp, [32]);
    assert_eq!(report.cells[1].native_product_count, 0);
    assert!(report.cells[1].enumeration_complete);
    assert_eq!(
        report.subjects[1].role,
        screening::Role::Decoy(screening::DecoyKind::EmptyVector)
    );
    args["components"][0]["origin"] = json!("OriginalReference");
    // The original-reference declaration does not strip saved-product markers.
    // The native cell retains the resulting unsupported evidence, not fake zero signal.
    let refused = call(&client, "assess_screening_strategy", args).await;
    let refused = artifact(&client, &refused["artifact_id"]).await;
    assert!(!refused.cells[0].enumeration_complete);
    assert_eq!(
        refused.tests[0].comparisons[0].distinction,
        screening::Distinction::Unresolved
    );
    client.cancel().await.unwrap();
    server.await.unwrap();
}
