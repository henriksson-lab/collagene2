//! Actual SDK lookup/evidence and complete catalog-bound inspection artifacts.
use collagene_core::{
    cloning::restriction::{catalog::*, inspection::InspectionOptions},
    formats::sequence as formats,
};
use collagene_mcp::server::CollageneServer;
use collagene_mcp::tools::cloning::RestrictionInspectionEvidence;
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
async fn reaction_conditions_preserve_formulation_differences_unknowns_and_exact_bindings() {
    use collagene_core::cloning::restriction::conditions::{ConditionsReport, ReactionConditions};
    let (client, service) = setup().await;
    let source = import(&client, "unchanged source", "AAGAATTCTT").await;
    let before = call(&client, "list_documents", json!({})).await;
    let catalog = builtin_catalog().unwrap();
    let selection = select(&catalog, &["ecori".into()]).unwrap();
    let options = json!({
        "schema_version":1,
        "formulations":[{"enzyme_id":"ecori","formulation_id":"neb-r0101"}],
        "buffer":"rCutSmart","temperature_milli_c":37000,
        "incubation_minutes":60,"minimum_activity_percent":75,
        "activity_context":"NEB unit-assay relative activity",
        "temperature_context":"vendor recommended incubation"
    });
    let arguments = json!({"selection":selection,"options":options});
    let ordinary = call(&client, "assess_restriction_conditions", arguments.clone()).await;
    assert_eq!(ordinary["sequence_assessed"], false);
    assert_eq!(ordinary["experimental_success_established"], false);
    let full = artifact(&client, &ordinary["artifact_id"]).await;
    assert_eq!(full["documentary_status"], "BelowThreshold");
    let report: ConditionsReport = serde_json::from_value(full.clone()).unwrap();
    let declared: ReactionConditions = serde_json::from_value(options.clone()).unwrap();
    report
        .validate_for(&catalog, &selection, &declared)
        .unwrap();
    report
        .validate_binding(ordinary["revision"].as_str().unwrap())
        .unwrap();
    assert!(
        full["sources"]
            .as_array()
            .unwrap()
            .iter()
            .any(|s| { s["url"].as_str().unwrap().starts_with("https://") })
    );
    assert!(full.to_string().contains("star activity"));
    let repeated = call(&client, "assess_restriction_conditions", arguments.clone()).await;
    assert_eq!(repeated["artifact_id"], ordinary["artifact_id"]);

    let mut high_fidelity = arguments.clone();
    high_fidelity["options"]["formulations"][0]["formulation_id"] = json!("neb-r3101");
    let hf = call(
        &client,
        "assess_restriction_conditions",
        high_fidelity.clone(),
    )
    .await;
    let hf_full = artifact(&client, &hf["artifact_id"]).await;
    assert_eq!(hf_full["documentary_status"], "Known");
    assert_ne!(hf["revision"], ordinary["revision"]);
    let hf_options: ReactionConditions =
        serde_json::from_value(high_fidelity["options"].clone()).unwrap();
    assert!(
        report
            .validate_for(&catalog, &selection, &hf_options)
            .is_err()
    );
    assert!(
        report
            .validate_binding(hf["revision"].as_str().unwrap())
            .is_err()
    );

    let mixed_selection = select(&catalog, &["ecori".into(), "bamhi".into()]).unwrap();
    let mut missing = high_fidelity.clone();
    missing["selection"] = json!(mixed_selection);
    missing["options"]["formulations"] = json!([
        {"enzyme_id":"ecori","formulation_id":"neb-r3101"},
        {"enzyme_id":"bamhi","formulation_id":null}
    ]);
    let mixed = call(&client, "assess_restriction_conditions", missing.clone()).await;
    let mixed_full = artifact(&client, &mixed["artifact_id"]).await;
    assert_eq!(mixed_full["documentary_status"], "Unknown");
    assert_eq!(mixed_full["formulations"][1]["formulation"], Value::Null);
    assert!(
        report
            .validate_for(&catalog, &mixed_selection, &declared)
            .is_err()
    );

    let mut context = arguments.clone();
    context["options"]["activity_context"] = json!("unrepresented assay");
    let changed_context = call(&client, "assess_restriction_conditions", context.clone()).await;
    assert_ne!(changed_context["revision"], ordinary["revision"]);
    let context_options: ReactionConditions =
        serde_json::from_value(context["options"].clone()).unwrap();
    assert!(
        report
            .validate_for(&catalog, &selection, &context_options)
            .is_err()
    );

    let mut stale = arguments.clone();
    stale["selection"]["catalog_content_hash"] = json!("0".repeat(64));
    reject(&client, "assess_restriction_conditions", stale).await;
    let mut forged_geometry = arguments.clone();
    forged_geometry["selection"]["enzymes"][0]["enzyme"]["top_cut"] = json!(2);
    reject(&client, "assess_restriction_conditions", forged_geometry).await;
    let mut wrong_enzyme = arguments.clone();
    wrong_enzyme["options"]["formulations"][0]["formulation_id"] = json!("neb-r0141");
    reject(&client, "assess_restriction_conditions", wrong_enzyme).await;
    let mut unsupported = arguments.clone();
    unsupported["options"]["schema_version"] = json!(2);
    reject(&client, "assess_restriction_conditions", unsupported).await;
    let mut invented_field = arguments;
    invented_field["document"] = source;
    reject(&client, "assess_restriction_conditions", invented_field).await;

    let mut forged_report = full;
    forged_report["sources"][0]["url"] = json!("https://example.invalid/unreviewed");
    let forged: ConditionsReport = serde_json::from_value(forged_report).unwrap();
    assert!(
        forged
            .validate_for(&catalog, &selection, &declared)
            .is_err()
    );
    assert_eq!(call(&client, "list_documents", json!({})).await, before);
    client.cancel().await.unwrap();
    service.await.unwrap();
}

#[tokio::test]
async fn catalog_discovery_is_versioned_and_assay_qualifiers_unknown_and_citations_are_preserved() {
    let (c, s) = setup().await;
    assert_eq!(
        c.list_tools(Default::default()).await.unwrap().tools.len(),
        collagene_mcp::CURRENT_TOOL_COUNT
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
    let current = builtin_catalog().unwrap();
    assert_eq!(
        first["catalog_schema_version"],
        current.edition().schema_version
    );
    assert_eq!(
        uri,
        format!(
            "collagene://enzyme-catalog/v{}/{}/{}",
            current.edition().schema_version,
            current.edition().edition,
            current.content_hash()
        )
    );
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
        c.read_resource(ReadResourceRequestParams::new(format!(
            "collagene://enzyme-catalog/v{}/latest",
            current.edition().schema_version
        )))
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
        json!({"less_than_percent":10})
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
async fn expanded_digest_publishes_one_exact_catalog_bound_artifact() {
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
    assert!(first.get("catalog_artifact_id").is_none());
    let full: RestrictionInspectionEvidence =
        serde_json::from_value(artifact(&c, &first["artifact_id"]).await).unwrap();
    let RestrictionInspectionEvidence::Original { inspection: pinned } = full else {
        panic!("ordinary inspection must retain its explicit source kind");
    };
    assert_eq!(pinned.inspection.sites[0].top_cut, Some(4));
    assert_eq!(pinned.inspection.sites[0].bottom_cut, Some(4));
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
async fn single_catalog_artifact_refuses_full_session_without_partial_publication() {
    let (c, s) = setup().await;
    let d = import(&c, "capacity", "AAGGCCTT").await;
    let before = call(&c, "list_documents", json!({})).await;
    for work in 1..=128 {
        let mut o = collagene_core::chemistry::methylation::MethylationOptions::default();
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
    let evidence = RestrictionInspectionEvidence::Original { inspection: report };
    let artifact_id = blake3::hash(&serde_json::to_vec(&evidence).unwrap())
        .to_hex()
        .to_string();
    reject(
        &c,
        "inspect_restriction",
        json!({"document":d,"source_role":"capacity-role","enzymes":["haeiii"]}),
    )
    .await;
    reject(
        &c,
        "read_artifact",
        json!({"artifact_id":artifact_id,"start":0,"max_bytes":1024}),
    )
    .await;
    // Refusing one new item must not invalidate any already published item.
    let mut o = collagene_core::chemistry::methylation::MethylationOptions::default();
    o.max_comparisons = 128;
    call(&c, "assess_methylation", json!({"document":d,"options":o})).await;
    assert_eq!(call(&c, "list_documents", json!({})).await, before);
    c.cancel().await.unwrap();
    s.await.unwrap();
}
