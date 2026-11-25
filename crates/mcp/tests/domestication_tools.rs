//! Real SDK route from source features through native checking and portable steps.
use base64::{Engine, engine::general_purpose::STANDARD};
use collagene_core::{
    Feature, Location, Molecule, Topology, catalog_digest, codon_design, domestication,
};
use collagene_mcp::CollageneServer;
use rmcp::{RoleClient, ServiceExt, model::CallToolRequestParams, service::RunningService};
use serde_json::{Value, json};
type Client = RunningService<RoleClient, ()>;
async fn call(c: &Client, name: &str, args: Value) -> Value {
    let r = c
        .call_tool(
            CallToolRequestParams::new(name.to_string())
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
            CallToolRequestParams::new(name.to_string())
                .with_arguments(args.as_object().unwrap().clone()),
        )
        .await;
    assert!(r.is_err() || r.unwrap().is_error == Some(true));
}
async fn artifact(c: &Client, id: &Value) -> (String, String) {
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
            return (text, p["encoding"].as_str().unwrap().to_string());
        }
        start = p["next_start"].as_u64().unwrap();
    }
}
async fn setup() -> (Client, tokio::task::JoinHandle<()>, Value, Value) {
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
    let c = ().serve(b).await.unwrap();
    let mut source = Molecule::dna(
        "Domestication source",
        "ATGGGTCTCTAACCAGTATCGACGATATCACGTAACCTGCACTAGTTCAGGACTGACATCGAACGTTAGTC",
        Topology::Circular,
    )
    .unwrap();
    source.features.push(Feature {
        kind: "CDS".into(),
        location: Location::range(0, 12),
        qualifiers: vec![("gene".into(), Some("source coding feature".into()))],
        needs_review: false,
    });
    let gb = collagene_core::formats::export_with_report(
        &source,
        collagene_core::formats::TextFormat::GenBank,
    )
    .unwrap()
    .text;
    let imported = call(&c, "import_document", json!({"text":gb})).await;
    let row = &imported["documents"][0];
    let binding = json!({"id":row["id"],"revision":row["revision"]});
    let request = json!({"source_revision":row["revision"],"selection":catalog_digest::select(&catalog_digest::builtin_catalog().unwrap(),&["bsai".into()]).unwrap(),
        "scope":[{"start":0,"end":source.sequence.len()}],"retained_sites":[],"coding_features":[],"protected":[],"options":codon_design::DesignOptions::default(),"limits":domestication::Limits::default()});
    (c, server, binding, request)
}
#[tokio::test]
async fn source_pins_synonymous_proposal_real_preparation_and_html_pdf_oligos_share_one_plan() {
    let (c, server, source, mut request) = setup().await;
    assert_eq!(
        c.list_all_tools().await.unwrap().len(),
        collagene_mcp::CURRENT_TOOL_COUNT
    );
    let original = call(&c, "export_document", source.clone()).await;
    let inspected = call(
        &c,
        "inspect_domestication",
        json!({"source":source,"request":request}),
    )
    .await;
    assert_eq!(inspected["site_count"], 1);
    assert_eq!(inspected["verified"], false);
    request["coding_features"] = json!([inspected["features"][0]["pin"]]);
    let proposed = call(
        &c,
        "propose_domestication",
        json!({"source":source,"request":request}),
    )
    .await;
    assert_eq!(proposed["status"], "CandidateFound");
    assert_eq!(proposed["verified"], false);
    let selection = proposed["candidate"].clone();
    assert!(!selection["patches"].as_array().unwrap().is_empty());
    let checked = call(
        &c,
        "check_domestication",
        json!({"source":source,"request":request,"selection":selection}),
    )
    .await;
    assert_eq!(checked["state"], "checked-reference-intent");
    assert_eq!(checked["verified"], false);
    let positions = selection["patches"]
        .as_array()
        .unwrap()
        .iter()
        .map(|p| p["position"].as_u64().unwrap())
        .collect::<Vec<_>>();
    let route = json!({"InversePcr":{"edit_start":positions.iter().min().unwrap(),"edit_end":positions.iter().max().unwrap()+1,"forward_length":8,"reverse_length":8,"reverse_tail_bases":0}});
    let prepared = call(
        &c,
        "prepare_domestication",
        json!({"source":source,"request":{"design":request,"selection":selection,"route":route}}),
    )
    .await;
    assert_eq!(prepared["verified"], false);
    assert_eq!(
        prepared["target"]["revision"],
        checked["target"]["revision"]
    );
    assert_eq!(
        prepared["verification"]["status"],
        "ConsistentUnderAssumptions"
    );
    let plan = json!({"plan_revision":prepared["plan_revision"]});
    let fresh = call(&c, "verify_plan", plan.clone()).await;
    assert_eq!(fresh["status"], "ConsistentUnderAssumptions");
    let html = call(&c, "export_protocol", plan.clone()).await;
    let (text, encoding) = artifact(&c, &html["artifact_id"]).await;
    assert_eq!(encoding, "utf-8");
    assert!(text.contains("PCR"));
    let pdf = call(&c, "export_protocol_pdf", plan.clone()).await;
    let (text, encoding) = artifact(&c, &pdf["artifact_id"]).await;
    assert_eq!(encoding, "base64");
    assert!(STANDARD.decode(text).unwrap().starts_with(b"%PDF-"));
    let oligos = call(
        &c,
        "export_plan_oligos",
        json!({"plan_revision":prepared["plan_revision"],"format":"Csv"}),
    )
    .await;
    assert_eq!(oligos["complete"], true);
    assert_eq!(oligos["use_count"], 2);
    assert_eq!(call(&c, "export_document", source).await, original);
    c.cancel().await.unwrap();
    server.abort();
}
#[tokio::test]
async fn stale_feature_source_or_candidate_and_uncovered_physical_route_are_refused() {
    let (c, server, source, mut request) = setup().await;
    let inspected = call(
        &c,
        "inspect_domestication",
        json!({"source":source,"request":request}),
    )
    .await;
    request["coding_features"] = json!([inspected["features"][0]["pin"]]);
    let proposed = call(
        &c,
        "propose_domestication",
        json!({"source":source,"request":request}),
    )
    .await;
    let selection = proposed["candidate"].clone();
    let mut stale = request.clone();
    stale["coding_features"][0]["revision"] = json!("0".repeat(64));
    reject(
        &c,
        "propose_domestication",
        json!({"source":source,"request":stale}),
    )
    .await;
    let mut wrong = selection.clone();
    wrong["candidate_id"] = json!("0".repeat(64));
    reject(
        &c,
        "check_domestication",
        json!({"source":source,"request":request,"selection":wrong}),
    )
    .await;
    let mut wrong_source = source.clone();
    wrong_source["revision"] = json!("0".repeat(64));
    reject(
        &c,
        "inspect_domestication",
        json!({"source":wrong_source,"request":request}),
    )
    .await;
    reject(&c,"prepare_domestication",json!({"source":source,"request":{"design":request,"selection":selection,"route":{"InversePcr":{"edit_start":20,"edit_end":21,"forward_length":8,"reverse_length":8,"reverse_tail_bases":0}}}})).await;
    c.cancel().await.unwrap();
    server.abort();
}
