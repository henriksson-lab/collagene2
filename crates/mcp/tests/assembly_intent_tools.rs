//! Declared selection intent is independent of later primer/PCR reconstruction.
use collagene_core::primers::DesignOptions;
use collagene_mcp::CollageneServer;
use rmcp::{RoleClient, ServiceExt, model::CallToolRequestParams, service::RunningService};
use serde_json::{Value, json};
type Client = RunningService<RoleClient, ()>;
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
async fn evidence(c: &Client, id: &Value) -> Value {
    let mut text = String::new();
    let mut start = 0;
    loop {
        let page = call(
            c,
            "read_artifact",
            json!({"artifact_id":id,"start":start,"max_bytes":16384}),
        )
        .await;
        text.push_str(page["text"].as_str().unwrap());
        if page["done"] == true {
            break;
        }
        start = page["next_start"].as_u64().unwrap();
    }
    serde_json::from_str(&text).unwrap()
}
fn named(imported: &Value, name: &str) -> Value {
    let row = imported["documents"]
        .as_array()
        .unwrap()
        .iter()
        .find(|d| d["name"] == name)
        .unwrap();
    json!({"id":row["id"],"revision":row["revision"]})
}
#[tokio::test]
async fn declared_regions_build_frozen_unverified_target_before_actual_primer_design_and_pcr() {
    let (client, server) = setup().await;
    assert_eq!(
        client
            .list_tools(Default::default())
            .await
            .unwrap()
            .tools
            .len(),
        collagene_mcp::CURRENT_TOOL_COUNT
    );
    let imported = call(
        &client,
        "import_document",
        json!({"text":">left\nATGACCCCGGGGTTTCCAGA\n>right\nATCCGTTAGCGACCTTGGAA\n"}),
    )
    .await;
    let left = named(&imported, "left");
    let right = named(&imported, "right");
    let original = call(&client, "export_document", left.clone()).await;
    let inputs = json!([{"id":"l","document":left},{"id":"r","document":right}]);
    let fragments = json!([{"role":"left","input_id":"l","start":0,"end":20,"reverse":false},{"role":"right","input_id":"r","start":0,"end":20,"reverse":false}]);
    let intent_fragments = json!([{"role":"left","input_id":"l","selector":{"ReferenceRange":{"start":0,"end":20,"reverse":false}}},{"role":"right","input_id":"r","selector":{"ReferenceRange":{"start":0,"end":20,"reverse":false}}}]);
    let request = json!({"name":"Declared intended construct","fragments":intent_fragments,"circular":false,"material_bindings":[]});
    let args = json!({"origin_workspace":null,"inputs":inputs,"request":request});
    let intent = call(&client, "build_assembly_intent", args.clone()).await;
    assert_eq!(intent["state"], "intended-reference");
    assert_eq!(intent["verified"], false);
    assert_eq!(intent["preparation_simulated"], false);
    assert!(intent.get("plan_revision").is_none());
    assert_eq!(intent["target"].as_object().unwrap().len(), 2);
    assert_eq!(intent["intended_length"], 40);
    assert_eq!(intent["fragment_count"], 2);
    let sequence = call(
        &client,
        "read_sequence_region",
        json!({"document":intent["target"],"start":0,"end":40}),
    )
    .await;
    assert_eq!(
        sequence["sequence"],
        "ATGACCCCGGGGTTTCCAGAATCCGTTAGCGACCTTGGAA"
    );
    assert_eq!(intent["fragments"][0]["input_revision"], left["revision"]);
    assert_eq!(intent["fragments"][0]["output_start"], 0);
    assert_eq!(intent["fragments"][0]["output_end"], 20);
    assert_eq!(intent["fragments"][1]["input_revision"], right["revision"]);
    assert_eq!(intent["fragments"][1]["output_start"], 20);
    assert_eq!(intent["fragments"][1]["output_end"], 40);
    let full = evidence(&client, &intent["artifact_id"]).await;
    assert_eq!(full["inputs"], inputs);
    assert_eq!(
        full["report"]["target"]["metadata"]["assembly:intent:reference_only"],
        "true"
    );
    assert!(
        full["report"]["target"]["metadata"]["duplex_model"]
            .as_str()
            .unwrap()
            .contains("Reference-only")
    );
    assert_eq!(
        call(&client, "build_assembly_intent", args.clone()).await["target"],
        intent["target"]
    );
    // No target was copied from a PCR/overlap answer: its40 literal bases above
    // were already fixed solely by the declared source selections.
    let mut options = DesignOptions::default();
    options.target_start = 0;
    options.target_end = 20;
    options.min_length = 8;
    options.max_length = 8;
    options.max_tm_difference_c = 100.;
    options.min_gc_fraction = 0.;
    options.max_gc_fraction = 1.;
    options.max_homopolymer = 200;
    options.max_three_prime_complementarity = 200;
    options.max_pairs = 1;
    let design_request = json!({"name":"Prepare frozen selection intent","fragments":fragments,"overlap":15,"circular":false,"primer_options":options,"profile_options":{"method":"HiFi","scan_repeated_overlaps":false,"max_comparisons":100000,"max_occurrences":100}});
    let designed = call(
        &client,
        "design_overlap_primers",
        json!({"inputs":inputs,"target":intent["target"],"request":design_request}),
    )
    .await;
    assert_eq!(designed["verified"], false);
    let checked = call(
        &client,
        "propose_overlap_assembly",
        json!({"inputs":inputs,"target":intent["target"],"recipe":designed["recipe"]}),
    )
    .await;
    assert_eq!(
        checked["verification"]["status"],
        "ConsistentUnderAssumptions"
    );
    assert_eq!(checked["verified"], false); // stored candidate remains a candidate
    let mut wrong_source = args.clone();
    wrong_source["inputs"][0]["document"] = right.clone();
    let wrong_intent = call(&client, "build_assembly_intent", wrong_source).await;
    assert_eq!(wrong_intent["verified"], false);
    assert_ne!(wrong_intent["context_revision"], intent["context_revision"]);
    assert_eq!(
        wrong_intent["fragments"][0]["input_revision"],
        right["revision"]
    );
    // Rebinding selects different intent, never makes the original preparations
    // satisfy that different target through an identity/name shortcut.
    let wrong = call(
        &client,
        "propose_overlap_assembly",
        json!({"inputs":inputs,"target":wrong_intent["target"],"recipe":designed["recipe"]}),
    )
    .await;
    assert_eq!(wrong["verification"]["status"], "Invalid");
    let mut reversed = args.clone();
    reversed["request"]["name"] = json!("Reverse selected region intent");
    reversed["request"]["fragments"][0]["selector"]["ReferenceRange"]["end"] = json!(4);
    reversed["request"]["fragments"][1]["selector"]["ReferenceRange"]["end"] = json!(4);
    reversed["request"]["fragments"][1]["selector"]["ReferenceRange"]["reverse"] = json!(true);
    let reversed_intent = call(&client, "build_assembly_intent", reversed.clone()).await;
    assert_eq!(
        call(
            &client,
            "read_sequence_region",
            json!({"document":reversed_intent["target"],"start":0,"end":8})
        )
        .await["sequence"],
        "ATGAGGAT"
    );
    reversed["request"]["circular"] = json!(true);
    assert_eq!(
        call(&client, "build_assembly_intent", reversed).await["topology"],
        "Circular"
    );
    let mut stale = args.clone();
    stale["inputs"][0]["document"]["revision"] = json!("stale");
    reject(&client, "build_assembly_intent", stale).await;
    let mut invalid = args.clone();
    invalid["request"]["fragments"][0]["selector"]["ReferenceRange"]["end"] = json!(21);
    reject(&client, "build_assembly_intent", invalid).await;
    let mut missing = args.clone();
    missing["request"]["fragments"][0]["input_id"] = json!("missing");
    reject(&client, "build_assembly_intent", missing).await;
    let mut duplicate = args;
    duplicate["request"]["fragments"][1]["role"] = json!("left");
    reject(&client, "build_assembly_intent", duplicate).await;
    assert_eq!(
        call(&client, "export_document", left.clone()).await["artifact_id"],
        original["artifact_id"]
    );
    assert_eq!(call(&client, "list_documents", json!({})).await["total"], 6);
    drop(client);
    server.abort();
}
