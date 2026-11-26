//! The ordinary MCP recipe tool mutates a real insert before multipart assembly.
use collagene_core::{
    analysis::Primer,
    assembly_recipe::{self, FragmentSpec, Preparation, Recipe},
    overlap_extension,
    overlap_extension_pcr::{PrimerPair, Request as ExtensionRequest},
    verification::{Operation, Source, Status},
};
use collagene_mcp::CollageneServer;
use rmcp::{RoleClient, ServiceExt, model::CallToolRequestParams, service::RunningService};
use serde_json::{Value, json};
type Client = RunningService<RoleClient, ()>;
const ORIGINAL: &str = "AGTCCGATGCTACGTAGGACTTCAGACCTGTAGATCGTCA";
const MUTANT: &str = "AGTCCGATGCTACGTATTTTCAGACCTGTAGATCGTCA";
const VECTOR: &str = "GATCGTCACAGTACGGTACCTAGC";
const TARGET: &str = "AGTCCGATGCTACGTATTTTCAGACCTGTAGATCGTCACAGTACGGTACCTAGC";
const WRONG: &str = "AGTCCGATGCTACGTAGGACTTCAGACCTGTAGATCGTCACAGTACGGTACCTAGC";
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
async fn text(c: &Client, id: &Value) -> String {
    let (mut start, mut out) = (0, String::new());
    loop {
        let chunk = call(
            c,
            "read_artifact",
            json!({"artifact_id":id,"start":start,"max_bytes":16384}),
        )
        .await;
        out.push_str(chunk["text"].as_str().unwrap());
        if chunk["done"] == true {
            return out;
        }
        start = chunk["next_start"].as_u64().unwrap();
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
fn pair(f: &str, ft: &str, r: &str, rt: &str) -> PrimerPair {
    PrimerPair {
        forward: Primer::new("Actual forward", f, ft).unwrap(),
        reverse: Primer::new("Actual reverse", r, rt).unwrap(),
    }
}
fn recipe() -> Recipe {
    Recipe {
        name: "Replace GGAC with TT before inserting into vector".into(),
        fragments: vec![
            FragmentSpec {
                role: "edited_insert".into(),
                input_id: "insert".into(),
                preparation: Preparation::OverlapExtensionPcr {
                    request: ExtensionRequest {
                        left: pair("AGTCCGAT", "", "TACGTAGC", "GGTCTGAAAA"),
                        right: pair("TTCAGACC", "GCTACGTATT", "TGACGATC", ""),
                        outer: pair("AGTCCGAT", "", "TGACGATC", ""),
                        extension: overlap_extension::Options {
                            overlap: 18,
                            min_paired_bases: 8,
                            ..Default::default()
                        },
                    },
                },
                reverse: false,
            },
            FragmentSpec {
                role: "vector".into(),
                input_id: "vector".into(),
                preparation: Preparation::ExactPcr {
                    forward: Primer::new("Vector F", "GATCGTCA", "").unwrap(),
                    reverse: Primer::new("Vector R", "GCTAGGTA", "").unwrap(),
                },
                reverse: false,
            },
        ],
        overlaps: vec![8],
        circular: false,
        options: Default::default(),
    }
}
#[tokio::test]
async fn ordinary_recipe_sdk_retains_both_branches_independent_target_diagnostics_and_ordered_oligos()
 {
    let (server_io, client_io) = tokio::io::duplex(65536);
    let server = tokio::spawn(async move {
        CollageneServer::default()
            .serve(server_io)
            .await
            .unwrap()
            .waiting()
            .await
            .unwrap();
    });
    let client = ().serve(client_io).await.unwrap();
    let tools = client.list_all_tools().await.unwrap();
    let tool = tools
        .iter()
        .find(|tool| tool.name == "propose_overlap_assembly")
        .unwrap();
    assert!(
        serde_json::to_string(&tool.input_schema)
            .unwrap()
            .contains("OverlapExtensionPcr")
    );
    let imported = call(&client, "import_document", json!({"text":format!(">original\n{ORIGINAL}\n>vector\n{VECTOR}\n>independent_target\n{TARGET}\n>wrong_target\n{WRONG}\n")})).await;
    let original = named(&imported, "original");
    let vector = named(&imported, "vector");
    let args = json!({"inputs":[{"id":"insert","document":original},{"id":"vector","document":vector}],"target":named(&imported,"independent_target"),"recipe":recipe()});
    let proposal = call(&client, "propose_overlap_assembly", args.clone()).await;
    assert_eq!(proposal["verified"], false);
    assert_eq!(
        proposal["verification"]["status"], "ConsistentUnderAssumptions",
        "{proposal}"
    );
    let mut evidence: Value =
        serde_json::from_str(&text(&client, &proposal["artifact_id"]).await).unwrap();
    let diagnostics: collagene_core::assembly_diagnostics::Report = serde_json::from_value(
        evidence
            .as_object_mut()
            .unwrap()
            .remove("diagnostics")
            .unwrap(),
    )
    .unwrap();
    let preview: assembly_recipe::Preview = serde_json::from_value(evidence).unwrap();
    preview.validate().unwrap();
    assert_eq!(
        preview.verification.status,
        Status::ConsistentUnderAssumptions
    );
    assert_eq!(diagnostics.context_revision, preview.context_revision);
    assert_eq!(diagnostics.rows[0].total_steps, 4);
    assert_eq!(diagnostics.rows[0].completed_steps, 4);
    assert_eq!(
        preview.compiled.plan.inputs[0].molecule.sequence.as_str(),
        ORIGINAL
    );
    assert_eq!(preview.prepared[0].sequence.as_str(), MUTANT);
    assert_eq!(
        preview
            .verification
            .products
            .last()
            .unwrap()
            .sequence
            .as_str(),
        TARGET
    );
    assert!(
        matches!(&preview.compiled.plan.steps[0], Operation::Pcr { template: Source::Input(id), .. } if id == "insert")
    );
    assert!(
        matches!(&preview.compiled.plan.steps[1], Operation::Pcr { template: Source::Input(id), .. } if id == "insert")
    );
    assert!(
        matches!(&preview.compiled.plan.steps[2], Operation::OverlapExtension { left, right, .. } if matches!(left.source,Source::Step(0)) && matches!(right.source,Source::Step(1)))
    );
    assert!(matches!(
        &preview.compiled.plan.steps[3],
        Operation::Pcr {
            template: Source::Step(2),
            ..
        }
    ));
    let reference = json!({"plan_revision":proposal["plan_revision"]});
    assert_eq!(
        call(&client, "verify_plan", reference.clone()).await["status"],
        "ConsistentUnderAssumptions"
    );
    let html = call(&client, "export_protocol", reference.clone()).await;
    assert!(
        text(&client, &html["artifact_id"])
            .await
            .contains("reanneal")
    );
    let pdf = call(&client, "export_protocol_pdf", reference).await;
    assert_eq!(pdf["mime_type"], "application/pdf");
    assert_eq!(pdf["status"], "ConsistentUnderAssumptions");
    let oligos = call(
        &client,
        "export_plan_oligos",
        json!({"plan_revision":proposal["plan_revision"],"format":"Csv"}),
    )
    .await;
    let table = text(&client, &oligos["table_artifact_id"]).await;
    assert!(table.contains("GGTCTGAAAATACGTAGC"));
    assert!(table.contains("GCTACGTATTTTCAGACC"));
    let mut wrong = args.clone();
    wrong["target"] = named(&imported, "wrong_target");
    let rejected = call(&client, "propose_overlap_assembly", wrong).await;
    assert_eq!(rejected["verification"]["status"], "Invalid");
    let mut broken = args;
    broken["recipe"]["name"] = json!("Bad selected terminal homology is not repaired from target");
    broken["recipe"]["fragments"][0]["preparation"]["OverlapExtensionPcr"]["request"]["extension"]
        ["overlap"] = json!(17);
    let failed = call(&client, "propose_overlap_assembly", broken).await;
    assert_eq!(failed["verification"]["status"], "Incomplete");
    let broken_evidence: Value =
        serde_json::from_str(&text(&client, &failed["artifact_id"]).await).unwrap();
    assert_eq!(
        broken_evidence["diagnostics"]["failure"]["stage"],
        "OverlapExtension"
    );
    assert_eq!(broken_evidence["diagnostics"]["failure"]["row"], 0);
    assert_eq!(
        broken_evidence["verification"]["products"]
            .as_array()
            .unwrap()
            .len(),
        2
    );
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
            json!({"document":vector,"start":0,"end":24})
        )
        .await["sequence"],
        VECTOR
    );
    drop(client);
    server.await.unwrap();
}
