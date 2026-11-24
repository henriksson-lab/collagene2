//! Actual SDK lookup/evidence and compatible flat plus fully pinned inspection artifacts.
use collagene_core::{
    catalog_digest::*,
    formats,
    restriction_inspection::{InspectionOptions, RestrictionInspection},
};
use collagene_mcp::CollageneServer;
use rmcp::{
    RoleClient, ServiceExt,
    model::{CallToolRequestParams, ReadResourceRequestParams},
    service::RunningService,
};
use serde_json::{Value, json};
type Client = RunningService<RoleClient, ()>;
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
async fn import(c: &Client, name: &str, dna: &str) -> Value {
    let r = call(
        c,
        "import_document",
        json!({"text":format!(">{name}\n{dna}\n")}),
    )
    .await;
    json!({"id":r["documents"][0]["id"],"revision":r["documents"][0]["revision"]})
}
#[tokio::test]
async fn catalog_discovery_is_versioned_and_assay_qualifiers_unknown_and_citations_are_preserved() {
    let (c, s) = setup().await;
    assert_eq!(
        c.list_tools(Default::default()).await.unwrap().tools.len(),
        97
    );
    let legacy = call(&c, "list_enzymes", json!({})).await;
    assert_eq!(legacy["enzymes"].as_array().unwrap().len(), 12);
    assert!(legacy["scope"].as_str().unwrap().contains("geometry only"));
    let first = call(&c, "lookup_enzyme_catalog", json!({"query":"","limit":1})).await;
    assert_eq!(first["enzymes"].as_array().unwrap().len(), 1);
    assert_eq!(first["done"], false);
    let next = call(
        &c,
        "lookup_enzyme_catalog",
        json!({"query":"","limit":1,"offset":1}),
    )
    .await;
    assert_ne!(first["enzymes"][0]["id"], next["enzymes"][0]["id"]);
    assert_eq!(first["catalog_content_hash"], next["catalog_content_hash"]);
    let enzyme = call(
        &c,
        "lookup_enzyme_catalog",
        json!({"query":"HaeIII","mode":"Enzyme"}),
    )
    .await;
    assert_eq!(enzyme["status"], "Unique");
    assert_eq!(enzyme["enzymes"][0]["id"], "haeiii");
    assert_eq!(
        enzyme["enzymes"][0]["geometry"],
        json!({"recognition":"GGCC","top_cut":2,"bottom_cut":2})
    );
    assert!(
        !enzyme["enzymes"][0]["citations"]
            .as_array()
            .unwrap()
            .is_empty()
    );
    let resources = c.list_all_resources().await.unwrap();
    let uri = first["catalog_resource_uri"].as_str().unwrap();
    assert!(uri.contains("/v1/"));
    assert!(uri.ends_with(first["catalog_content_hash"].as_str().unwrap()));
    assert!(resources.iter().any(|r| r.uri == uri));
    let resource = c
        .read_resource(ReadResourceRequestParams::new(uri))
        .await
        .unwrap();
    let resource = serde_json::to_value(resource).unwrap();
    let payload: Value =
        serde_json::from_str(resource["contents"][0]["text"].as_str().unwrap()).unwrap();
    assert_eq!(
        payload["identity"]["catalog_content_hash"],
        first["catalog_content_hash"]
    );
    assert_eq!(
        payload["catalog"]["enzymes"].as_array().unwrap().len(),
        first["count"].as_u64().unwrap() as usize
    );
    assert!(
        c.read_resource(ReadResourceRequestParams::new(
            "collagene://enzyme-catalog/v1/latest"
        ))
        .await
        .is_err()
    );
    let formulation = call(
        &c,
        "lookup_enzyme_catalog",
        json!({"query":"EcoRI-HF","mode":"Formulation"}),
    )
    .await;
    assert_eq!(formulation["formulations"][0]["enzyme"], "ecori");
    assert_eq!(formulation["formulations"][0]["id"], "neb-r3101");
    let known=call(&c,"resolve_enzyme_evidence",json!({"formulation":"neb-r3101","condition":{"kind":"buffer_activity","buffer":"NEBuffer r3.1","temperature_milli_c":37000}})).await;
    assert_eq!(known["formulation_id"], "neb-r3101");
    assert_eq!(known["enzyme_id"], "ecori");
    assert_eq!(known["resolution"]["status"], "Known");
    assert_eq!(
        known["resolution"]["value"],
        json!({"kind":"less_than_percent","value":10})
    );
    let evidence = &known["resolution"]["evidence"][0];
    assert_eq!(evidence["context"], "NEB unit-assay relative activity");
    assert!(
        evidence["qualifiers"][0]
            .as_str()
            .unwrap()
            .contains("not predicted digestion efficiency")
    );
    assert!(!evidence["citations"].as_array().unwrap().is_empty());
    assert!(
        known["sources"]
            .as_array()
            .unwrap()
            .iter()
            .any(|v| v["url"].as_str().unwrap().starts_with("https://"))
    );
    let unknown=call(&c,"resolve_enzyme_evidence",json!({"formulation":"neb-r3101","condition":{"kind":"buffer_activity","buffer":"not-curated","temperature_milli_c":37000}})).await;
    assert_eq!(unknown["resolution"]["status"], "Unknown");
    assert!(unknown["resolution"].get("value").is_none());
    reject(&c,"lookup_enzyme_catalog",json!({"query":"EcoRI","catalog_pin":{"catalog_edition":first["catalog_edition"],"catalog_content_hash":"0".repeat(64)}})).await;
    reject(&c,"resolve_enzyme_evidence",json!({"formulation":"neb-r3101","condition":{"kind":"buffer_activity","buffer":"rCutSmart","temperature_milli_c":-1}})).await;
    c.cancel().await.unwrap();
    s.await.unwrap();
}
#[tokio::test]
async fn expanded_digest_keeps_strict_legacy_artifact_and_adds_exact_catalog_pin() {
    let (c, s) = setup().await;
    let d = import(&c, "source", "AAGGCCTT").await;
    let before = call(&c, "list_documents", json!({})).await;
    let args =
        json!({"document":d,"source_role":"source-event","enzymes":["HaeIII"],"max_fragments":1});
    let first = call(&c, "inspect_restriction", args.clone()).await;
    assert_eq!(first["scan_complete"], true);
    assert_eq!(first["products_determined"], true);
    assert_eq!(first["fragments"][0]["duplex"]["top_length"], 4);
    assert_eq!(first["site_count"], 1);
    assert_eq!(first["fragment_count"], 2);
    assert_eq!(first["catalog_pin"]["enzymes"][0]["id"], "haeiii");
    assert_ne!(first["catalog_batch_revision"], first["batch_revision"]);
    let raw = artifact(&c, &first["artifact_id"]).await;
    assert!(raw.get("catalog_pin").is_none());
    let legacy: RestrictionInspection = serde_json::from_value(raw.clone()).unwrap();
    assert_eq!(legacy.sites[0].top_cut, Some(4));
    assert_eq!(legacy.sites[0].bottom_cut, Some(4));
    let pinned: CatalogRestrictionInspection =
        serde_json::from_value(artifact(&c, &first["catalog_artifact_id"]).await).unwrap();
    assert_eq!(serde_json::to_value(&pinned.inspection).unwrap(), raw);
    assert_eq!(pinned.batch_revision, first["catalog_batch_revision"]);
    let catalog = builtin_catalog().unwrap();
    let selection = select(&catalog, &["haeiii".into()]).unwrap();
    let molecule = formats::import_text(">source\nAAGGCCTT\n")
        .unwrap()
        .molecules
        .remove(0);
    pinned
        .validate_bindings_for(
            &catalog,
            &selection,
            "source-event",
            &molecule,
            &InspectionOptions::default(),
        )
        .unwrap();
    let mut next = args.clone();
    next["fragment_start"] = json!(1);
    next["catalog_pin"] = first["catalog_pin"].clone();
    let next = call(&c, "inspect_restriction", next).await;
    assert_eq!(next["artifact_id"], first["artifact_id"]);
    assert_eq!(next["catalog_artifact_id"], first["catalog_artifact_id"]);
    assert_eq!(next["fragments"][0]["top_start"], 4);
    let mut stale = args.clone();
    stale["catalog_pin"] = first["catalog_pin"].clone();
    stale["catalog_pin"]["catalog_content_hash"] = json!("0".repeat(64));
    reject(&c, "inspect_restriction", stale).await;
    let mut forged = args.clone();
    forged["catalog_pin"] = first["catalog_pin"].clone();
    forged["catalog_pin"]["enzymes"][0]["enzyme"]["top_cut"] = json!(1);
    reject(&c, "inspect_restriction", forged).await;
    reject(
        &c,
        "analyze_restriction",
        json!({"document":d,"enzyme":"HaeIII"}),
    )
    .await;
    assert_eq!(call(&c, "list_documents", json!({})).await, before);
    c.cancel().await.unwrap();
    s.await.unwrap();
}
#[tokio::test]
async fn paired_catalog_and_legacy_artifacts_publish_atomically_at_session_capacity() {
    let (c, s) = setup().await;
    let d = import(&c, "capacity", "AAGGCCTT").await;
    let before = call(&c, "list_documents", json!({})).await;
    for work in 1..=127 {
        let mut o = collagene_core::methylation::MethylationOptions::default();
        o.max_comparisons = work;
        call(&c, "assess_methylation", json!({"document":d,"options":o})).await;
    }
    let catalog = builtin_catalog().unwrap();
    let selection = select(&catalog, &["haeiii".into()]).unwrap();
    let molecule = formats::import_text(">capacity\nAAGGCCTT\n")
        .unwrap()
        .molecules
        .remove(0);
    let report = inspect(
        &catalog,
        &selection,
        "capacity-role",
        &molecule,
        &InspectionOptions::default(),
    )
    .unwrap();
    let raw_id = blake3::hash(&serde_json::to_vec(&report.inspection).unwrap())
        .to_hex()
        .to_string();
    let pinned_id = blake3::hash(&serde_json::to_vec(&report).unwrap())
        .to_hex()
        .to_string();
    reject(
        &c,
        "inspect_restriction",
        json!({"document":d,"source_role":"capacity-role","enzymes":["haeiii"]}),
    )
    .await;
    for id in [raw_id, pinned_id] {
        reject(
            &c,
            "read_artifact",
            json!({"artifact_id":id,"start":0,"max_bytes":1024}),
        )
        .await;
    }
    let mut o = collagene_core::methylation::MethylationOptions::default();
    o.max_comparisons = 128;
    call(&c, "assess_methylation", json!({"document":d,"options":o})).await;
    assert_eq!(call(&c, "list_documents", json!({})).await, before);
    c.cancel().await.unwrap();
    s.await.unwrap();
}
