//! Actual SDK reconstruction of explicitly declared generation and checked digestion plans.
use base64::{Engine, engine::general_purpose::STANDARD};
use collagene_mcp::server::CollageneServer;
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
async fn checked_digest_plans_pin_catalogs_rebuild_physical_state_and_export_without_trusting_metadata()
 {
    let (c, s) = setup().await;
    let imported = call(
        &c,
        "import_document",
        json!({"text":">source\nAAAAGAATTCCCC\n>target\nAAAAG\n"}),
    )
    .await;
    let source = doc(&imported, "source");
    let target = doc(&imported, "target");
    let enzymes = call(&c, "list_enzymes", json!({})).await;
    let profile = enzymes["checked_digest_profile_revision"].clone();
    assert_eq!(profile.as_str().unwrap().len(), 64);
    let catalog = collagene_core::cloning::restriction::catalog::builtin_catalog().unwrap();
    let enzyme = catalog
        .edition()
        .enzymes
        .iter()
        .find(|e| e.name == "EcoRI")
        .unwrap()
        .id
        .clone();
    let selection =
        collagene_core::cloning::restriction::catalog::select(&catalog, &[enzyme]).unwrap();
    let options = collagene_core::chemistry::tracked_duplex::checked_inspection_options();
    let steps = json!([{"DeclareSourceModifications":{"template":{"Input":"s"},"assumptions":{"host":null,"dam":"Unmodified","dcm":"Unmodified","cpg":"Unmodified"}}},{"CheckedDigest":{"conditions":null,"template":{"Step":0},"selection":selection,"options":options,"fragment_index":0}}]);
    let propose = |steps: Value| json!({"name":"checked source generation","target":target,"inputs":[{"id":"s","document":source}],"steps":steps});
    let candidate = call(&c, "propose_plan", propose(steps.clone())).await;
    assert_eq!(candidate["verified"], false);
    let plan = json!({"plan_revision":candidate["plan_revision"]});
    let report = call(&c, "verify_plan", plan.clone()).await;
    assert_eq!(report["status"], "ConsistentUnderAssumptions", "{report}");
    assert_eq!(
        report["engine_version"].as_str().unwrap(),
        collagene_core::verification::ENGINE_VERSION
    );
    assert_eq!(report["steps"].as_array().unwrap().len(), 2);
    assert!(report["steps"][1]["product"]["tracked_modifications"].is_object());
    let html = call(&c, "export_protocol", plan.clone()).await;
    let html = call(
        &c,
        "read_artifact",
        json!({"artifact_id":html["artifact_id"],"start":0,"max_bytes":16384}),
    )
    .await;
    assert!(
        html["text"]
            .as_str()
            .unwrap()
            .contains(collagene_core::verification::ENGINE_VERSION)
    );
    let pdf = call(&c, "export_protocol_pdf", plan.clone()).await;
    let pdf = call(
        &c,
        "read_artifact",
        json!({"artifact_id":pdf["artifact_id"],"start":0,"max_bytes":16384}),
    )
    .await;
    assert!(
        STANDARD
            .decode(pdf["text"].as_str().unwrap())
            .unwrap()
            .starts_with(b"%PDF")
    );
    let exported = call(&c, "export_candidate", plan).await;
    let bytes = call(
        &c,
        "read_artifact",
        json!({"artifact_id":exported["artifact_id"],"start":0,"max_bytes":16384}),
    )
    .await;
    let bytes = STANDARD.decode(bytes["text"].as_str().unwrap()).unwrap();
    assert_eq!(
        &bytes[4..6],
        &collagene_core::protocol::schema::CANDIDATE_SCHEMA_VERSION.to_le_bytes()
    );
    collagene_core::documents::transfer::load_candidate(&bytes).unwrap();
    let mut wrong_catalog = steps.clone();
    wrong_catalog[1]["CheckedDigest"]["selection"]["catalog_content_hash"] = json!("0".repeat(64));
    let candidate = call(&c, "propose_plan", propose(wrong_catalog)).await;
    let report = call(
        &c,
        "verify_plan",
        json!({"plan_revision":candidate["plan_revision"]}),
    )
    .await;
    assert_ne!(report["status"], "ConsistentUnderAssumptions");
    let invalid_pdf = call(
        &c,
        "export_protocol_pdf",
        json!({"plan_revision":candidate["plan_revision"]}),
    )
    .await;
    assert_eq!(invalid_pdf["status"], "Invalid");
    let direct = json!([{"CheckedDigest":{"conditions":null,"template":{"Input":"s"},"selection":selection,"options":options,"fragment_index":0}}]);
    let candidate = call(&c, "propose_plan", propose(direct)).await;
    let report = call(
        &c,
        "verify_plan",
        json!({"plan_revision":candidate["plan_revision"]}),
    )
    .await;
    assert_ne!(report["status"], "ConsistentUnderAssumptions");
    let mut unknown = steps;
    unknown[0]["DeclareSourceModifications"]["assumptions"]["unrelated_modifications_absent"] =
        json!(true);
    reject(&c, "propose_plan", propose(unknown)).await;
    c.cancel().await.unwrap();
    s.await.unwrap();
}
