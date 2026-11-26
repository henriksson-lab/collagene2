//! AI proposals use the same freshly reconstructed reciprocal-extension physics.
use collagene_core::{
    analysis::Primer,
    overlap_extension::Options,
    verification::{Operation, OrientedSource, Source},
};
use collagene_mcp::CollageneServer;
use rmcp::{RoleClient, ServiceExt, model::CallToolRequestParams, service::RunningService};
use serde_json::{Value, json};
type Client = RunningService<RoleClient, ()>;
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
async fn document(client: &Client, name: &str, dna: &str) -> Value {
    let imported = call(
        client,
        "import_document",
        json!({"text":format!(">{name}\n{dna}\n")}),
    )
    .await;
    let d = &imported["documents"][0];
    json!({"id":d["id"],"revision":d["revision"]})
}
async fn artifact(client: &Client, id: &Value) -> Value {
    let mut at = 0;
    let mut text = String::new();
    loop {
        let page = call(
            client,
            "read_artifact",
            json!({"artifact_id":id,"start":at,"max_bytes":16384}),
        )
        .await;
        text.push_str(page["text"].as_str().unwrap());
        if page["done"] == true {
            break;
        }
        at = page["next_start"].as_u64().unwrap();
    }
    serde_json::from_str(&text).unwrap()
}
fn operations() -> Vec<Operation> {
    vec![
        Operation::Pcr {
            template: Source::Input("original".into()),
            forward: Primer::new("Left outer", "AGTCCGAT", "").unwrap(),
            reverse: Primer::new("Left mutagenic", "TACGTAGC", "GGTCTGAAAA").unwrap(),
        },
        Operation::Pcr {
            template: Source::Input("original".into()),
            forward: Primer::new("Right mutagenic", "TTCAGACC", "GCTACGTATT").unwrap(),
            reverse: Primer::new("Right outer", "TGACGATC", "").unwrap(),
        },
        Operation::OverlapExtension {
            left: OrientedSource {
                source: Source::Step(0),
                reverse: false,
            },
            right: OrientedSource {
                source: Source::Step(1),
                reverse: false,
            },
            options: Options {
                overlap: 18,
                min_paired_bases: 8,
                ..Default::default()
            },
        },
        Operation::Pcr {
            template: Source::Step(2),
            forward: Primer::new("Outer forward", "AGTCCGAT", "").unwrap(),
            reverse: Primer::new("Outer reverse", "TGACGATC", "").unwrap(),
        },
    ]
}
#[tokio::test]
async fn sdk_proposal_verification_and_exports_use_actual_extension_not_target_stock() {
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
    let client = ().serve(b).await.unwrap();
    let tools = client.list_tools(Default::default()).await.unwrap().tools;
    assert_eq!(tools.len(), collagene_mcp::CURRENT_TOOL_COUNT);
    let propose = tools
        .iter()
        .find(|tool| tool.name == "propose_plan")
        .unwrap();
    assert!(
        serde_json::to_string(&propose.input_schema)
            .unwrap()
            .contains("OverlapExtension")
    );
    assert!(
        propose
            .description
            .as_ref()
            .unwrap()
            .contains("OverlapExtension")
    );
    let original = document(
        &client,
        "Original unchanged stock",
        "AGTCCGATGCTACGTAGGACTTCAGACCTGTAGATCGTCA",
    )
    .await;
    let target = document(
        &client,
        "Independent GGAC to TT requirement",
        "AGTCCGATGCTACGTATTTTCAGACCTGTAGATCGTCA",
    )
    .await;
    let args = json!({"name":"Actual reciprocal extension","target":target,"inputs":[{"id":"original","document":original}],"steps":operations()});
    let candidate = call(&client, "propose_plan", args.clone()).await;
    assert_eq!(candidate["state"], "candidate");
    assert_eq!(candidate["verified"], false);
    let reference = json!({"plan_revision":candidate["plan_revision"]});
    let verified = call(&client, "verify_plan", reference.clone()).await;
    assert_eq!(
        verified["status"], "ConsistentUnderAssumptions",
        "{verified}"
    );
    assert!(
        verified["engine_version"]
            .as_str()
            .unwrap()
            .ends_with("/science-23")
    );
    let products = verified["products"].as_array().unwrap();
    assert_eq!(products.len(), 4);
    assert_eq!(products[0]["length"], 26);
    assert_eq!(products[1]["length"], 30);
    assert_eq!(products[2]["length"], 38);
    assert_eq!(products[3]["length"], 38);
    let evidence = artifact(&client, &verified["artifact_id"]).await;
    let extension = evidence["physical_processing_reports"]
        .as_array()
        .unwrap()
        .iter()
        .find(|step| step["step_index"] == 2)
        .unwrap();
    let geometry = &extension["evidence"]["OverlapExtension"];
    assert_eq!(geometry["options"]["overlap"], 18);
    assert_eq!(geometry["terminal_pairings"]["complete"], true);
    assert_eq!(geometry["terminal_pairings"]["unique"], true);
    let html = call(&client, "export_protocol", reference.clone()).await;
    assert_eq!(html["mime_type"], "text/html; charset=utf-8");
    let page = call(
        &client,
        "read_artifact",
        json!({"artifact_id":html["artifact_id"],"start":0,"max_bytes":16384}),
    )
    .await;
    assert!(
        page["text"]
            .as_str()
            .unwrap()
            .to_lowercase()
            .contains("extension")
    );
    let pdf = call(&client, "export_protocol_pdf", reference).await;
    assert_eq!(pdf["mime_type"], "application/pdf");
    assert_eq!(pdf["status"], "ConsistentUnderAssumptions");
    assert!(pdf["bytes"].as_u64().unwrap() > 100);
    let mut wrong_target = args.clone();
    wrong_target["name"] = json!("Wrong independent target");
    wrong_target["target"] = original.clone();
    let wrong = call(&client, "propose_plan", wrong_target).await;
    let wrong_report = call(
        &client,
        "verify_plan",
        json!({"plan_revision":wrong["plan_revision"]}),
    )
    .await;
    assert_eq!(wrong_report["status"], "Invalid");
    // A sequence-only overlap cannot manufacture the physical template needed
    // by the final PCR, even when its sequence matches the independent target.
    let mut reference_only = operations();
    reference_only[2] = Operation::Overlap {
        fragments: vec![Source::Step(0), Source::Step(1)],
        overlaps: vec![18],
        circular: false,
    };
    let mut documentary = args;
    documentary["name"] = json!("Reference overlap cannot feed physical PCR");
    documentary["steps"] = serde_json::to_value(reference_only).unwrap();
    let refused = call(&client, "propose_plan", documentary).await;
    let withheld = call(
        &client,
        "verify_plan",
        json!({"plan_revision":refused["plan_revision"]}),
    )
    .await;
    assert_eq!(withheld["status"], "Incomplete");
    assert_eq!(withheld["products"].as_array().unwrap().len(), 3);
    assert!(
        withheld["findings"]
            .as_array()
            .unwrap()
            .iter()
            .any(|finding| {
                finding
                    .as_str()
                    .is_some_and(|text| text.contains("established physical template"))
            }),
        "{withheld}"
    );
    assert_eq!(
        call(
            &client,
            "read_sequence_region",
            json!({"document":original,"start":0,"end":40})
        )
        .await["sequence"],
        "AGTCCGATGCTACGTAGGACTTCAGACCTGTAGATCGTCA"
    );
    drop(client);
    server.abort();
}
