//! Native auto-primers remain a read-only proposal until ordinary recipe verification.
use collagene_core::{
    assembly_edit_intent::Edit,
    overlap_extension_primer_design::{Options, Request},
    overlap_primer_design::Fragment,
    overlap_profiles::ProfileOptions,
    primers::DesignOptions,
};
use collagene_mcp::CollageneServer;
use rmcp::{RoleClient, ServiceExt, model::CallToolRequestParams, service::RunningService};
use serde_json::{Value, json};
type Client = RunningService<RoleClient, ()>;
const ORIGINAL: &str = "AGTCCGATGCTACGTAGGACTTCAGACCTGTAGATCGTCA";
const EXPECTED: &str = "AGTCCGATGCTACGTATTTTCAGACCTGTAGATCGTCA";
fn request() -> Request {
    Request {
        fragment: Fragment {
            role: "insert".into(),
            input_id: "original".into(),
            start: 0,
            end: 40,
            reverse: false,
        },
        edits: vec![Edit {
            role: "insert".into(),
            start: 16,
            end: 20,
            expected: "GGAC".into(),
            replacement: "TT".into(),
        }],
        options: Options {
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
            homology_flank: 8,
            min_paired_bases: 8,
            ..Default::default()
        },
    }
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
    let result = c
        .call_tool(
            CallToolRequestParams::new(name.to_owned())
                .with_arguments(args.as_object().unwrap().clone()),
        )
        .await;
    assert!(result.is_err() || result.unwrap().is_error == Some(true));
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
async fn six_auto_oligos_keep_original_stock_independent_intent_and_full_recipe_authority_separate()
{
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
    let tools = client.list_all_tools().await.unwrap();
    assert_eq!(tools.len(), collagene_mcp::CURRENT_TOOL_COUNT);
    let tool = tools
        .iter()
        .find(|t| t.name == "design_overlap_extension_primers")
        .unwrap();
    let schema = serde_json::to_string(&tool.input_schema).unwrap();
    assert!(schema.contains("edits") && schema.contains("homology_flank"));
    let imported = call(
        &client,
        "import_document",
        json!({"text":format!(">original\n{ORIGINAL}\n>independent_target\n{EXPECTED}\n")}),
    )
    .await;
    let original = named(&imported, "original");
    let target = named(&imported, "independent_target");
    let documents_before = call(&client, "list_documents", json!({})).await;
    let args = json!({"document":original,"request":request()});
    let proposal = call(&client, "design_overlap_extension_primers", args.clone()).await;
    assert_eq!(proposal["state"], "advisory");
    assert_eq!(proposal["verified"], false);
    assert_eq!(
        proposal["verification"]["status"], "ConsistentUnderAssumptions",
        "{proposal}"
    );
    assert_eq!(proposal["expected"]["length"], 38);
    assert_eq!(proposal["design_count"], 3);
    assert_eq!(proposal["source_revision"], original["revision"]);
    let oligos = proposal["ordered_oligos"].as_array().unwrap();
    assert_eq!(oligos.len(), 6);
    assert!(oligos.iter().any(|p| p["stage"] == "left_pcr"
        && p["direction"] == "reverse"
        && p["ordered_sequence"] == "GGTCTGAAAATACGTAGC"));
    assert!(oligos.iter().any(|p| p["stage"] == "right_pcr"
        && p["direction"] == "forward"
        && p["ordered_sequence"] == "GCTACGTATTTTCAGACC"));
    let evidence = artifact(&client, &proposal["artifact_id"]).await;
    assert_eq!(evidence["document"], original);
    let report: collagene_core::overlap_extension_primer_design::Report =
        serde_json::from_value(evidence["report"].clone()).unwrap();
    report
        .validate_binding(proposal["context_revision"].as_str().unwrap())
        .unwrap();
    assert_eq!(report.expected.sequence.as_str(), EXPECTED);
    assert_eq!(report.plan.inputs[0].molecule.sequence.as_str(), ORIGINAL);
    assert_eq!(
        report
            .verification
            .products
            .iter()
            .map(|m| m.sequence.len())
            .collect::<Vec<_>>(),
        vec![26, 30, 38, 38]
    );
    assert_eq!(
        call(&client, "list_documents", json!({})).await,
        documents_before,
        "the proposal must not install stock, a target or a product document"
    );
    let mut changed = args.clone();
    changed["request"]["options"]["max_candidates"] = json!(99);
    let changed = call(&client, "design_overlap_extension_primers", changed).await;
    assert_ne!(changed["context_revision"], proposal["context_revision"]);
    let sequences = |value: &Value| {
        value["ordered_oligos"]
            .as_array()
            .unwrap()
            .iter()
            .map(|p| p["ordered_sequence"].clone())
            .collect::<Vec<_>>()
    };
    assert_eq!(sequences(&changed), sequences(&proposal));
    let mut stale = args.clone();
    stale["document"]["revision"] = json!("0".repeat(64));
    reject(&client, "design_overlap_extension_primers", stale).await;
    let mut bad = args;
    bad["request"]["edits"][0]["expected"] = json!("AAAA");
    reject(&client, "design_overlap_extension_primers", bad).await;
    let recipe = json!({"name":"Use advisory auto-primers in an independently targeted assembly","fragments":[{"role":"insert","input_id":"original","preparation":{"OverlapExtensionPcr":{"request":proposal["preparation"]}},"reverse":false}],"overlaps":[],"circular":false,"options":ProfileOptions::default()});
    let assembled = call(
        &client,
        "propose_overlap_assembly",
        json!({"inputs":[{"id":"original","document":original}],"target":target,"recipe":recipe}),
    )
    .await;
    assert_eq!(assembled["verified"], false);
    assert_eq!(
        assembled["verification"]["status"], "ConsistentUnderAssumptions",
        "{assembled}"
    );
    let verified = call(
        &client,
        "verify_plan",
        json!({"plan_revision":assembled["plan_revision"]}),
    )
    .await;
    assert_eq!(verified["status"], "ConsistentUnderAssumptions");
    let wrong = call(
        &client,
        "propose_overlap_assembly",
        json!({"inputs":[{"id":"original","document":original}],"target":original,"recipe":recipe}),
    )
    .await;
    assert_eq!(wrong["verification"]["status"], "Invalid");
    assert_eq!(
        call(
            &client,
            "read_sequence_region",
            json!({"document":original,"start":0,"end":40})
        )
        .await["sequence"],
        ORIGINAL
    );
    drop(client);
    server.await.unwrap();
}
