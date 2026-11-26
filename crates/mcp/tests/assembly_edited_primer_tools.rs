//! A whole edited-assembly proposal remains source-bound and separate from target retention.
use base64::{Engine, engine::general_purpose::STANDARD};
use collagene_core::{
    Molecule, Topology,
    analysis::Primer,
    assembly_edit_intent::{self, Edit},
    assembly_edited_primer_design::{Report, Request, RoleOrigin},
    assembly_intent, assembly_junction_edit, assembly_substitution_intent, cloning_outputs,
    overlap_extension_primer_design::Options,
    overlap_primer_design::Fragment,
    primers::DesignOptions,
    run_lineage::{EventBinding, InputEventBinding},
    run_material::{self, Selection},
    verification::{DesignPlan, Operation, PlanInput, Source},
    workspace::{self, Workspace},
};
use collagene_mcp::CollageneServer;
use rmcp::{RoleClient, ServiceExt, model::CallToolRequestParams, service::RunningService};
use serde_json::{Value, json};
type Client = RunningService<RoleClient, ()>;
const ORIGINAL: &str = "AGTCCGATGCTACGTAGGACTTCAGACCTGTAGATCGTCA";
const VECTOR: &str = "GCTTACCGATGACCTAGTACGGCATTCGAAGC";
const TARGET: &str = "AGTCCGATGCTACGTATTTTCAGACCTGTAGATCGTCAATATGCTTACCGATGACCTAGTACGGCATTCGAAGC";
fn request() -> Request {
    Request {
        intent: assembly_edit_intent::Request {
            base: assembly_substitution_intent::Request {
                base: assembly_intent::Request {
                    name: "Independent whole edited assembly".into(),
                    fragments: vec![
                        Fragment {
                            role: "insert".into(),
                            input_id: "original".into(),
                            start: 0,
                            end: 40,
                            reverse: false,
                        },
                        Fragment {
                            role: "vector".into(),
                            input_id: "vector".into(),
                            start: 0,
                            end: 32,
                            reverse: false,
                        },
                    ],
                    circular: false,
                },
                junction_edits: vec![assembly_junction_edit::Edit {
                    junction: 0,
                    remove_left: 0,
                    remove_right: 0,
                    inserted: "ATAT".into(),
                }],
                substitutions: vec![],
            },
            edits: vec![Edit {
                role: "insert".into(),
                start: 16,
                end: 20,
                expected: "GGAC".into(),
                replacement: "TT".into(),
            }],
        },
        overlap: 8,
        options: Options {
            homology_flank: 8,
            min_paired_bases: 8,
            primer_options: DesignOptions {
                min_length: 8,
                max_length: 8,
                min_gc_fraction: 0.,
                max_gc_fraction: 1.,
                max_tm_difference_c: 100.,
                max_homopolymer: 200,
                max_three_prime_complementarity: 200,
                ..Default::default()
            },
            ..Default::default()
        },
        profile_options: Default::default(),
        origins: vec![],
    }
}
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
            return serde_json::from_str(&text).unwrap();
        }
        start = p["next_start"].as_u64().unwrap();
    }
}
fn named(imported: &Value, name: &str) -> Value {
    let d = imported["documents"]
        .as_array()
        .unwrap()
        .iter()
        .find(|d| d["name"] == name)
        .unwrap();
    json!({"id":d["id"],"revision":d["revision"]})
}
#[tokio::test]
async fn whole_proposal_uses_literal_junction_and_native_oe_then_separate_intent_retention_and_review()
 {
    let (client, server) = setup().await;
    let tools = client.list_all_tools().await.unwrap();
    assert_eq!(tools.len(), collagene_mcp::CURRENT_TOOL_COUNT);
    assert!(
        tools
            .iter()
            .any(|t| t.name == "design_edited_overlap_primers")
    );
    let imported = call(
        &client,
        "import_document",
        json!({"text":format!(">original\n{ORIGINAL}\n>vector\n{VECTOR}\n")}),
    )
    .await;
    let original = named(&imported, "original");
    let vector = named(&imported, "vector");
    let inputs = json!([{"id":"original","document":original},{"id":"vector","document":vector}]);
    let before = call(&client, "list_documents", json!({})).await;
    let args = json!({"inputs":inputs,"request":request()});
    let proposal = call(&client, "design_edited_overlap_primers", args.clone()).await;
    assert_eq!(proposal["state"], "advisory");
    assert_eq!(proposal["verified"], false);
    assert_eq!(
        proposal["verification"]["status"], "ConsistentUnderAssumptions",
        "{proposal}"
    );
    assert_eq!(proposal["target_summary"]["length"], 74);
    assert_eq!(proposal["ordered_oligos"].as_array().unwrap().len(), 8);
    assert!(
        proposal["ordered_oligos"]
            .as_array()
            .unwrap()
            .iter()
            .any(|p| p["role"] == "insert"
                && p["stage"] == "outer_pcr"
                && p["direction"] == "reverse"
                && p["ordered_sequence"] == "CGGTAAGCATATTGACGATC")
    );
    let evidence = artifact(&client, &proposal["artifact_id"]).await;
    let report: Report = serde_json::from_value(evidence["report"].clone()).unwrap();
    report
        .validate_binding(proposal["context_revision"].as_str().unwrap())
        .unwrap();
    assert_eq!(report.intent.target.sequence.as_str(), TARGET);
    assert_eq!(
        report
            .preview
            .verification
            .products
            .last()
            .unwrap()
            .sequence
            .as_str(),
        TARGET
    );
    assert_eq!(
        call(&client, "list_documents", json!({})).await,
        before,
        "advisory planning must not install sources, target or products"
    );
    let intended = call(
        &client,
        "build_contribution_edit_intent",
        json!({"inputs":inputs,"request":request().intent}),
    )
    .await;
    assert_eq!(
        intended["context_revision"],
        proposal["intent_context_revision"]
    );
    assert_eq!(
        intended["target"]["revision"],
        proposal["target_summary"]["revision"]
    );
    let prepared = call(
        &client,
        "propose_overlap_assembly",
        json!({"inputs":inputs,"target":intended["target"],"recipe":proposal["recipe"]}),
    )
    .await;
    assert_eq!(prepared["verified"], false);
    assert_eq!(
        prepared["verification"]["status"],
        "ConsistentUnderAssumptions"
    );
    assert_eq!(
        call(
            &client,
            "verify_plan",
            json!({"plan_revision":prepared["plan_revision"]})
        )
        .await["status"],
        "ConsistentUnderAssumptions"
    );
    let wrong = call(
        &client,
        "propose_overlap_assembly",
        json!({"inputs":inputs,"target":original,"recipe":proposal["recipe"]}),
    )
    .await;
    assert_eq!(wrong["verification"]["status"], "Invalid");
    let mut changed = args;
    changed["request"]["options"]["max_candidates"] = json!(99);
    let changed = call(&client, "design_edited_overlap_primers", changed).await;
    assert_ne!(changed["context_revision"], proposal["context_revision"]);
    assert_eq!(
        call(
            &client,
            "read_sequence_region",
            json!({"document":original,"start":0,"end":40})
        )
        .await["sequence"],
        ORIGINAL
    );
    assert_eq!(
        call(
            &client,
            "read_sequence_region",
            json!({"document":vector,"start":0,"end":32})
        )
        .await["sequence"],
        VECTOR
    );
    drop(client);
    server.await.unwrap();
}
fn binding(w: &Workspace, index: usize) -> EventBinding {
    let d = &w.documents[index];
    EventBinding {
        document_id: d.id.clone(),
        event_id: d.head_event.clone().unwrap(),
        revision: d.head.clone(),
    }
}
#[tokio::test]
async fn saved_role_requires_exact_workspace_attestation_before_any_full_design() {
    let (client, server) = setup().await;
    let source = Molecule::dna("Original source", ORIGINAL, Topology::Linear).unwrap();
    let target = Molecule::dna(
        "Independent original source target",
        ORIGINAL,
        Topology::Linear,
    )
    .unwrap();
    let vector = Molecule::dna("Vector source", VECTOR, Topology::Linear).unwrap();
    let plan = DesignPlan {
        name: "Original saved PCR".into(),
        target_revision: target.revision(),
        inputs: vec![PlanInput {
            id: "original".into(),
            revision: source.revision(),
            molecule: source.clone(),
        }],
        steps: vec![Operation::Pcr {
            template: Source::Input("original".into()),
            forward: Primer::new("F", "AGTCCGAT", "").unwrap(),
            reverse: Primer::new("R", "TGACGATC", "").unwrap(),
        }],
    };
    let w = Workspace::from_documents(
        "Saved full-design source",
        vec![source, target, vector.clone()],
        1,
    )
    .unwrap();
    let commit = cloning_outputs::finalize(
        &w,
        &plan,
        binding(&w, 1),
        vec![InputEventBinding {
            input_id: "original".into(),
            event: binding(&w, 0),
        }],
        &cloning_outputs::Options {
            name: "Saved original PCR".into(),
            show_prepared: false,
            save_primers: false,
            allow_incomplete: false,
        },
    )
    .unwrap();
    let selected = Selection {
        run_id: commit.receipt.run.run_id,
        step_index: 0,
        event: commit.receipt.construct,
    };
    let capsule = run_material::capsule(&commit.workspace, &selected).unwrap();
    let imported=call(&client,"import_workspace_archive",json!({"data_base64":STANDARD.encode(workspace::save_workspace(&commit.workspace).unwrap())})).await;
    let wr = json!({"id":imported["workspace"]["id"],"revision":imported["workspace"]["revision"]});
    let mut r = request();
    r.origins.push(RoleOrigin {
        role: "insert".into(),
        origin: capsule,
    });
    let args = json!({"inputs":[{"id":"original","document":{"id":selected.event.revision,"revision":selected.event.revision}},{"id":"vector","document":{"id":vector.revision(),"revision":vector.revision()}}],"request":r});
    reject(&client, "design_edited_overlap_primers", args.clone()).await;
    let mut bound = args;
    bound["origin_workspaces"] = json!([{"role":"insert","workspace":wr}]);
    let before = call(&client, "list_documents", json!({})).await;
    let proposal = call(&client, "design_edited_overlap_primers", bound.clone()).await;
    assert_eq!(
        proposal["verification"]["status"], "ConsistentUnderAssumptions",
        "{proposal}"
    );
    assert!(
        proposal["recipe"]["fragments"][0]["preparation"]
            .get("FromRun")
            .is_some()
    );
    assert_eq!(call(&client, "list_documents", json!({})).await, before);
    let mut wrong = bound.clone();
    wrong["origin_workspaces"][0]["role"] = json!("vector");
    reject(&client, "design_edited_overlap_primers", wrong).await;
    let mut wrong = bound;
    wrong["origin_workspaces"][0]["workspace"]["revision"] = json!("0".repeat(64));
    reject(&client, "design_edited_overlap_primers", wrong).await;
    drop(client);
    server.await.unwrap();
}
