//! Actual SDK dispatch for bounded source-bound inspection, without imported report authority.
use collagene_core::{cloning::restriction::inspection::*, formats::sequence as formats};
use collagene_mcp::server::CollageneServer;
use collagene_mcp::tools::cloning::RestrictionInspectionEvidence;
use rmcp::{RoleClient, ServiceExt, model::CallToolRequestParams, service::RunningService};
use serde_json::{Value, json};
type Client = RunningService<RoleClient, ()>;
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
async fn setup() -> (Client, tokio::task::JoinHandle<()>) {
    let (a, b) = tokio::io::duplex(65536);
    let s = tokio::spawn(async move {
        CollageneServer::default()
            .serve(a)
            .await
            .unwrap()
            .waiting()
            .await
            .unwrap();
    });
    (().serve(b).await.unwrap(), s)
}
async fn import(c: &Client, dna: &str) -> Value {
    let r = call(
        c,
        "import_document",
        json!({"text":format!(">source\n{dna}\n")}),
    )
    .await;
    json!({"id":r["documents"][0]["id"],"revision":r["documents"][0]["revision"]})
}
async fn artifact(c: &Client, id: &Value) -> Value {
    let mut text = String::new();
    let mut start = 0;
    loop {
        let p = call(
            c,
            "read_artifact",
            json!({"artifact_id":id,"start":start,"max_bytes":1024}),
        )
        .await;
        text.push_str(p["text"].as_str().unwrap());
        if p["done"] == true {
            break;
        }
        start = p["next_start"].as_u64().unwrap();
    }
    serde_json::from_str(&text).unwrap()
}

#[tokio::test]
async fn scan_sdk_preserves_recognition_uncertainty_and_never_constructs_products() {
    let (c, s) = setup().await;
    let tools = c.list_tools(Default::default()).await.unwrap().tools;
    assert_eq!(tools.len(), collagene_mcp::CURRENT_TOOL_COUNT);
    let tool = tools.iter().find(|t| t.name == "scan_restriction").unwrap();
    assert_eq!(tool.input_schema["additionalProperties"], false);
    assert!(
        tool.input_schema["properties"]
            .get("fragment_start")
            .is_none()
    );
    let exact = import(&c, "AAAAGAATTCCCCCGGATCCAAAA").await;
    let ambiguous = import(&c, "AAAAGAANTCAAAA").await;
    let edge = import(&c, "GGTCTCA").await;
    let before = call(&c, "list_documents", json!({})).await;
    let args = json!({"document":exact,"source_role":"scan-source-event","enzymes":["EcoRI","BamHI","NotI"],"max_sites":1});
    let first = call(&c, "scan_restriction", args.clone()).await;
    assert_eq!(first["scan_complete"], true);
    assert_eq!(first["products_determined"], false);
    assert_eq!(first["fragment_count"], 0);
    assert_eq!(first["ref_top_lengths"], Value::Null);
    assert_eq!(first["site_count"], 2);
    assert_eq!(first["sites_done"], false);
    assert_eq!(first["sites"][0]["recognition_start"], 4);
    assert_eq!(first["sites"][0]["top_cut"], 5);
    assert_eq!(first["sites"][0]["bottom_cut"], 9);
    assert_eq!(first["enzyme_counts"][0]["certain_observed"], 1);
    assert_eq!(first["enzyme_counts"][0]["possible_observed"], 0);
    assert_eq!(first["enzyme_counts"][0]["single_cutter"], true);
    assert_eq!(first["enzyme_counts"][2]["exact_site_count"], 0);
    assert_eq!(first["enzyme_counts"][2]["single_cutter"], false);
    let evidence = artifact(&c, &first["artifact_id"]).await;
    let report: RestrictionInspection = serde_json::from_value(evidence["report"].clone()).unwrap();
    report.validate().unwrap();
    assert!(
        report.fragments.is_empty()
            && report.ref_top_lengths.is_none()
            && !report.products_determined
    );
    assert_eq!(evidence["catalog_pin"], first["catalog_pin"]);
    assert_eq!(evidence["enzyme_counts"], first["enzyme_counts"]);
    let mut second = args.clone();
    second["site_start"] = json!(1);
    second["catalog_pin"] = first["catalog_pin"].clone();
    let second = call(&c, "scan_restriction", second).await;
    assert_eq!(second["artifact_id"], first["artifact_id"]);
    assert_eq!(second["sites"][0]["recognition_start"], 14);
    let possible = call(
        &c,
        "scan_restriction",
        json!({"document":ambiguous,"source_role":"scan-source-event","enzymes":["EcoRI"]}),
    )
    .await;
    assert_eq!(possible["scan_complete"], true);
    assert_eq!(possible["enzyme_counts"][0]["possible_observed"], 1);
    assert_eq!(
        possible["enzyme_counts"][0]["exact_site_count"],
        Value::Null
    );
    assert_eq!(possible["enzyme_counts"][0]["single_cutter"], Value::Null);
    let edge = call(
        &c,
        "scan_restriction",
        json!({"document":edge,"source_role":"scan-source-event","enzymes":["BsaI"]}),
    )
    .await;
    assert_eq!(edge["scan_complete"], true);
    assert_eq!(edge["sites"][0]["top_cut"], 7);
    assert_eq!(edge["sites"][0]["bottom_cut"], Value::Null);
    assert_eq!(edge["sites"][0]["bottom_boundary"], 11);
    let mut partial = args.clone();
    let mut options = InspectionOptions::default();
    options.max_comparisons = 1;
    partial["options"] = serde_json::to_value(options).unwrap();
    let partial = call(&c, "scan_restriction", partial).await;
    assert_eq!(partial["scan_complete"], false);
    assert_eq!(partial["site_count"], 0);
    for count in partial["enzyme_counts"].as_array().unwrap() {
        assert_eq!(count["counts_complete"], false);
        assert_eq!(count["exact_site_count"], Value::Null);
        assert_eq!(count["single_cutter"], Value::Null);
    }
    for (key, value) in [
        ("max_sites", json!(0)),
        ("site_start", json!(999)),
        ("fragment_start", json!(0)),
        ("enzymes", json!(["not-an-enzyme"])),
    ] {
        let mut invalid = args.clone();
        invalid[key] = value;
        reject(&c, "scan_restriction", invalid).await;
    }
    let mut stale = args.clone();
    stale["document"]["revision"] = json!("0".repeat(64));
    reject(&c, "scan_restriction", stale).await;
    let mut forged = args;
    forged["catalog_pin"] = first["catalog_pin"].clone();
    forged["catalog_pin"]["enzymes"][0]["enzyme"]["top_cut"] = json!(2);
    reject(&c, "scan_restriction", forged).await;
    assert_eq!(call(&c, "list_documents", json!({})).await, before);
    c.cancel().await.unwrap();
    s.await.unwrap();
}

#[tokio::test]
async fn inspection_sdk_pages_exact_geometry_and_source_without_creating_products() {
    let (c, s) = setup().await;
    let tools = c.list_tools(Default::default()).await.unwrap().tools;
    assert_eq!(tools.len(), collagene_mcp::CURRENT_TOOL_COUNT);
    let tool = tools
        .iter()
        .find(|t| t.name == "inspect_restriction")
        .unwrap();
    assert_eq!(tool.input_schema["additionalProperties"], false);
    let dna = "AAAAGAATTCCCCCGGATCCAAAA";
    let d = import(&c, dna).await;
    let before = call(&c, "list_documents", json!({})).await;
    let args = json!({"document":d,"source_role":"source-event","enzymes":["ecori","BamHI"],"max_sites":1,"max_fragments":1});
    let first = call(&c, "inspect_restriction", args.clone()).await;
    assert_eq!(first["scan_complete"], true);
    assert_eq!(first["products_determined"], true);
    assert_eq!(first["site_count"], 2);
    assert_eq!(first["fragment_count"], 3);
    assert_eq!(first["sites_done"], false);
    assert_eq!(first["fragments_done"], false);
    assert_eq!(first["sites"][0]["recognition_length"], 6);
    assert_eq!(first["fragments"][0]["duplex"]["top_length"], 5);
    assert!(first.get("catalog_artifact_id").is_none());
    let evidence: RestrictionInspectionEvidence =
        serde_json::from_value(artifact(&c, &first["artifact_id"]).await).unwrap();
    let RestrictionInspectionEvidence::Original { inspection: pinned } = evidence else {
        panic!("Original reference inspection must retain its explicit tagged catalog evidence");
    };
    assert_eq!(pinned.batch_revision, first["catalog_batch_revision"]);
    assert_eq!(
        serde_json::to_value(&pinned.selection).unwrap(),
        first["catalog_pin"]
    );
    let molecule = formats::import_text(&format!(">source\n{dna}\n"))
        .unwrap()
        .molecules
        .remove(0);
    let catalog = collagene_core::cloning::restriction::catalog::builtin_catalog().unwrap();
    let selection: collagene_core::cloning::restriction::catalog::CatalogSelection =
        serde_json::from_value(first["catalog_pin"].clone()).unwrap();
    let mut selected_names = selection
        .enzymes
        .iter()
        .map(|e| e.enzyme.name.as_str())
        .collect::<Vec<_>>();
    selected_names.sort();
    assert_eq!(selected_names, vec!["BamHI", "EcoRI"]);
    pinned
        .validate_bindings_for(
            &catalog,
            &selection,
            "source-event",
            &molecule,
            &InspectionOptions::default(),
        )
        .unwrap();
    let report = &pinned.inspection;
    assert_eq!(report.batch_revision, first["batch_revision"]);
    assert_eq!(report.source_revision, d["revision"]);
    assert_eq!(
        report
            .fragments
            .iter()
            .map(|f| f.duplex.top_length)
            .collect::<Vec<_>>(),
        vec![5, 10, 9]
    );
    let mut next = args.clone();
    next["site_start"] = json!(1);
    next["fragment_start"] = json!(1);
    let next = call(&c, "inspect_restriction", next).await;
    assert_eq!(next["artifact_id"], first["artifact_id"]);
    assert_eq!(next["sites"][0]["recognition_start"], 14);
    assert_eq!(next["fragments"][0]["top_start"], 5);
    let mut other = args.clone();
    other["source_role"] = json!("other-event");
    let other = call(&c, "inspect_restriction", other).await;
    assert_ne!(other["batch_revision"], first["batch_revision"]);
    for (key, val) in [
        ("site_start", json!(999)),
        ("fragment_start", json!(999)),
        ("max_sites", json!(0)),
        ("max_fragments", json!(101)),
        ("enzymes", json!(["unknown"])),
        ("trusted_report", json!(true)),
    ] {
        let mut invalid = args.clone();
        invalid[key] = val;
        reject(&c, "inspect_restriction", invalid).await;
    }
    let mut stale = args;
    stale["document"]["revision"] = json!("0".repeat(64));
    reject(&c, "inspect_restriction", stale).await;
    assert_eq!(call(&c, "list_documents", json!({})).await, before);
    c.cancel().await.unwrap();
    s.await.unwrap();
}
#[tokio::test]
async fn possible_edge_and_incomplete_evidence_survive_while_legacy_tool_stays_strict() {
    let (c, s) = setup().await;
    for (dna, name) in [("AAAAGAANTCAAAA", "EcoRI"), ("GGTCTCA", "BsaI")] {
        let d = import(&c, dna).await;
        let r = call(
            &c,
            "inspect_restriction",
            json!({"document":d,"source_role":"source","enzymes":[name]}),
        )
        .await;
        assert_eq!(r["scan_complete"], true);
        assert_eq!(r["products_determined"], false);
        assert_eq!(r["fragment_count"], 0);
        assert!(!r["sites"][0]["diagnostics"].as_array().unwrap().is_empty());
        assert!(r.get("catalog_artifact_id").is_none());
        let evidence: RestrictionInspectionEvidence =
            serde_json::from_value(artifact(&c, &r["artifact_id"]).await).unwrap();
        let RestrictionInspectionEvidence::Original { inspection: pinned } = evidence else {
            panic!("Ordinary input must retain Original evidence");
        };
        pinned.validate().unwrap();
        assert_eq!(pinned.batch_revision, r["catalog_batch_revision"]);
        assert_eq!(
            serde_json::to_value(&pinned.selection).unwrap(),
            r["catalog_pin"]
        );
        assert_eq!(pinned.inspection.source_revision, d["revision"]);
        assert!(!pinned.inspection.products_determined);
        assert!(pinned.inspection.fragments.is_empty());
        reject(
            &c,
            "analyze_restriction",
            json!({"document":d,"enzyme":name}),
        )
        .await;
    }
    let d = import(&c, "AAAAGAATTCCCCCGGATCCAAAA").await;
    let mut o = InspectionOptions::default();
    o.max_comparisons = 1;
    let r = call(
        &c,
        "inspect_restriction",
        json!({"document":d,"source_role":"source","enzymes":["EcoRI"],"options":o}),
    )
    .await;
    assert_eq!(r["scan_complete"], false);
    assert_eq!(r["products_determined"], false);
    assert_eq!(r["limits"][0], "Comparisons");
    c.cancel().await.unwrap();
    s.await.unwrap();
}
