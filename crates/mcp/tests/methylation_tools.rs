//! Actual SDK routing for immutable source-bound methylation assessment evidence.
use collagene_core::{Molecule, methylation::*};
use collagene_mcp::CollageneServer;
use rmcp::{RoleClient, ServiceExt, model::CallToolRequestParams, service::RunningService};
use serde_json::{Value, json};
type Client = RunningService<RoleClient, ()>;
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
    assert!(result.is_err() || result.unwrap().is_error == Some(true));
}
async fn artifact(c: &Client, id: &Value) -> Value {
    let mut text = String::new();
    let mut at = 0;
    loop {
        let page = call(
            c,
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
fn doc(imported: &Value, name: &str) -> Value {
    let d = imported["documents"]
        .as_array()
        .unwrap()
        .iter()
        .find(|d| {
            d["name"]
                .as_str()
                .is_some_and(|s| s.split_whitespace().next() == Some(name))
        })
        .unwrap();
    json!({"id":d["id"],"revision":d["revision"]})
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
    let client = ().serve(b).await.unwrap();
    (client, server)
}
#[tokio::test]
async fn assessment_sdk_binds_source_profiles_states_and_full_evidence_pages() {
    let (c, s) = setup().await;
    assert_eq!(
        c.list_tools(Default::default()).await.unwrap().tools.len(),
        91
    );
    let imported = call(
        &c,
        "import_document",
        json!({"text":">source\nTCTAGATCGAATTCGGATCCCTCGAGGGTCTCCGTCTC\n"}),
    )
    .await;
    let source = doc(&imported, "source");
    let options = MethylationOptions::default();
    let before = call(&c, "list_documents", json!({})).await;
    let first = call(
        &c,
        "assess_methylation",
        json!({"document":source,"options":options,"max_sites":1}),
    )
    .await;
    assert_eq!(first["profiles"].as_array().unwrap().len(), 12);
    assert!(first["site_count"].as_u64().unwrap() > 1);
    assert_eq!(first["sites"].as_array().unwrap().len(), 1);
    let report: MethylationReport =
        serde_json::from_value(artifact(&c, &first["artifact_id"]).await).unwrap();
    assert_eq!(report.source_revision, source["revision"]);
    assert_eq!(report.options, options);
    assert!(
        report
            .sites
            .iter()
            .any(|site| site.outcome == Outcome::Unknown)
    );
    let molecule: Molecule =
        collagene_core::formats::import_text(">source\nTCTAGATCGAATTCGGATCCCTCGAGGGTCTCCGTCTC\n")
            .unwrap()
            .molecules
            .remove(0);
    report.validate_bindings_for(&molecule, &options).unwrap();
    replay(&molecule, &options, &report).unwrap();
    let next = call(
        &c,
        "assess_methylation",
        json!({"document":source,"options":options,"site_start":1,"max_sites":1}),
    )
    .await;
    assert_eq!(next["artifact_id"], first["artifact_id"]);
    assert_eq!(next["sites"][0]["index"], report.sites[1].index);
    let mut modified = options.clone();
    modified.dam = MethylationState::FullyModified;
    modified.dcm = MethylationState::FullyModified;
    modified.cpg = MethylationState::FullyModified;
    let modified = call(
        &c,
        "assess_methylation",
        json!({"document":source,"options":modified}),
    )
    .await;
    assert_ne!(modified["batch_revision"], first["batch_revision"]);
    assert!(
        modified["sites"]
            .as_array()
            .unwrap()
            .iter()
            .flat_map(|site| site["effects"].as_array().unwrap())
            .all(|effect| effect["state"] == "FullyModified")
    );
    let mut hemi = options.clone();
    hemi.dam = MethylationState::HemiOrMixed;
    hemi.dcm = MethylationState::HemiOrMixed;
    hemi.cpg = MethylationState::HemiOrMixed;
    let hemi = call(
        &c,
        "assess_methylation",
        json!({"document":source,"options":hemi}),
    )
    .await;
    assert!(
        hemi["sites"]
            .as_array()
            .unwrap()
            .iter()
            .any(|site| site["outcome"] == "Unknown")
    );
    let mut limited = options.clone();
    limited.max_comparisons = 1;
    let limited = call(
        &c,
        "assess_methylation",
        json!({"document":source,"options":limited}),
    )
    .await;
    assert_eq!(limited["complete"], false);
    assert!(!limited["limits"].as_array().unwrap().is_empty());
    let mut stale = source.clone();
    stale["revision"] = json!("0".repeat(64));
    reject(
        &c,
        "assess_methylation",
        json!({"document":stale,"options":options}),
    )
    .await;
    let mut unknown = serde_json::to_value(&options).unwrap();
    unknown["trust_saved_assessment"] = json!(true);
    reject(
        &c,
        "assess_methylation",
        json!({"document":source,"options":unknown}),
    )
    .await;
    let mut unsupported = options.clone();
    unsupported.enzymes = vec!["EcoRI-HF".into()];
    reject(
        &c,
        "assess_methylation",
        json!({"document":source,"options":unsupported}),
    )
    .await;
    reject(
        &c,
        "assess_methylation",
        json!({"document":source,"options":options,"site_start":usize::MAX}),
    )
    .await;
    assert_eq!(call(&c, "list_documents", json!({})).await, before);
    c.cancel().await.unwrap();
    s.await.unwrap();
}
#[tokio::test]
async fn report_capacity_rejection_does_not_mutate_documents_or_existing_artifacts() {
    let (c, s) = setup().await;
    let imported = call(
        &c,
        "import_document",
        json!({"text":">capacity\nTCTAGATC\n"}),
    )
    .await;
    let source = doc(&imported, "capacity");
    let before = call(&c, "list_documents", json!({})).await;
    let mut first = Value::Null;
    for work in 1..=128 {
        let mut options = MethylationOptions::default();
        options.max_comparisons = work;
        let result = call(
            &c,
            "assess_methylation",
            json!({"document":source,"options":options}),
        )
        .await;
        if work == 1 {
            first = result["artifact_id"].clone();
        }
    }
    let mut options = MethylationOptions::default();
    options.max_comparisons = 129;
    reject(
        &c,
        "assess_methylation",
        json!({"document":source,"options":options}),
    )
    .await;
    assert_eq!(call(&c, "list_documents", json!({})).await, before);
    let retained: MethylationReport = serde_json::from_value(artifact(&c, &first).await).unwrap();
    assert_eq!(retained.options.max_comparisons, 1);
    c.cancel().await.unwrap();
    s.await.unwrap();
}
