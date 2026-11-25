//! Explicit junction DNA intent precedes original-template PCR proposals and verification.
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
fn options() -> DesignOptions {
    DesignOptions {
        min_length: 8,
        max_length: 8,
        min_gc_fraction: 0.,
        max_gc_fraction: 1.,
        max_tm_difference_c: 100.,
        max_homopolymer: 200,
        max_three_prime_complementarity: 200,
        max_pairs: 1,
        ..Default::default()
    }
}
fn base(reverse: bool, circular: bool) -> Value {
    json!({"name":"Explicit edited assembly","fragments":[{"role":"left","input_id":"l","start":0,"end":20,"reverse":false},{"role":"right","input_id":"r","start":0,"end":20,"reverse":reverse}],"circular":circular})
}
fn edits() -> Value {
    json!([{"junction":0,"remove_left":2,"remove_right":2,"inserted":"GGA"}])
}
fn planner(base: &Value) -> Value {
    json!({"base":base,"edits":edits(),"overlap":15,"primer_options":options(),"profile_options":{"method":"HiFi","scan_repeated_overlaps":false,"max_comparisons":100000,"max_occurrences":100}})
}
#[tokio::test]
async fn declared_junction_trim_insert_targets_precede_primer_proposals_and_actual_assembly() {
    let (c, server) = setup().await;
    assert_eq!(
        c.list_tools(Default::default()).await.unwrap().tools.len(),
        collagene_mcp::CURRENT_TOOL_COUNT
    );
    let imported = call(
        &c,
        "import_document",
        json!({"text":">left\nATGACCCCGGGGTTTCCAGA\n>right\nATCCGTTAGCGACCTTGGAA\n"}),
    )
    .await;
    let left = named(&imported, "left");
    let right = named(&imported, "right");
    let inputs = json!([{"id":"l","document":left},{"id":"r","document":right}]);
    let left_original = call(&c, "export_document", left.clone()).await;
    let right_original = call(&c, "export_document", right.clone()).await;
    for (reverse, circular, literal, tail) in [
        (
            false,
            false,
            "ATGACCCCGGGGTTTCCAGGACCGTTAGCGACCTTGGAA",
            "CAAGGTCGCTAACGGTCC",
        ),
        (
            false,
            true,
            "ATGACCCCGGGGTTTCCAGGACCGTTAGCGACCTTGGAA",
            "CAAGGTCGCTAACGGTCC",
        ),
        (
            true,
            false,
            "ATGACCCCGGGGTTTCCAGGACCAAGGTCGCTAACGGAT",
            "CGTTAGCGACCTTGGTCC",
        ),
    ] {
        let base = base(reverse, circular);
        let args = json!({"inputs":inputs,"base":base,"edits":edits()});
        let intent = call(&c, "build_edited_assembly_intent", args.clone()).await;
        assert_eq!(intent["state"], "intended-reference");
        assert_eq!(intent["verified"], false);
        assert_eq!(intent["preparation_simulated"], false);
        assert!(intent.get("plan_revision").is_none());
        let sequence = call(
            &c,
            "read_sequence_region",
            json!({"document":intent["target"],"start":0,"end":literal.len()}),
        )
        .await;
        assert_eq!(sequence["sequence"], literal);
        assert_eq!(literal.len(), 39);
        let built = evidence(&c, &intent["artifact_id"]).await;
        assert_eq!(built["inputs"], inputs);
        assert_eq!(built["intent"]["base_request"], base);
        assert_eq!(built["intent"]["edits"], edits());
        assert_eq!(built["intent"]["fragments"][0]["start"], 0);
        assert_eq!(built["intent"]["fragments"][0]["end"], 18);
        assert_eq!(built["intent"]["fragments"][0]["output_start"], 0);
        assert_eq!(built["intent"]["fragments"][0]["output_end"], 21);
        assert_eq!(built["intent"]["fragments"][1]["output_start"], 21);
        assert_eq!(built["intent"]["fragments"][1]["output_end"], 39);
        let retry = call(&c, "build_edited_assembly_intent", args).await;
        assert_eq!(retry["target"], intent["target"]);
        assert_eq!(retry["artifact_id"], intent["artifact_id"]);
        // Expected DNA is already fixed above; the planner receives originals, not this target.
        let designed = call(
            &c,
            "design_edited_assembly_primers",
            json!({"inputs":inputs,"request":planner(&base)}),
        )
        .await;
        assert_eq!(designed["state"], "advisory");
        assert_eq!(designed["verified"], false);
        assert_eq!(
            designed["intent_context_revision"],
            intent["context_revision"]
        );
        assert_eq!(designed["target_revision"], intent["target"]["revision"]);
        let report = evidence(&c, &designed["artifact_id"]).await;
        assert_eq!(report["inputs"], inputs);
        assert_eq!(
            report["report"]["intent"]["target"]["sequence"]["bases"],
            literal
        );
        assert_eq!(designed["ordered_oligos"][0]["reverse"]["tail"], tail);
        assert_eq!(designed["recipe"]["fragments"][1]["reverse"], reverse);
        let prepared = call(
            &c,
            "propose_overlap_assembly",
            json!({"inputs":inputs,"target":intent["target"],"recipe":designed["recipe"]}),
        )
        .await;
        assert_eq!(
            prepared["verification"]["status"],
            "ConsistentUnderAssumptions"
        );
        let checked = call(
            &c,
            "verify_plan",
            json!({"plan_revision":prepared["plan_revision"]}),
        )
        .await;
        assert_eq!(checked["status"], "ConsistentUnderAssumptions");
        // An unedited reference is a different independently declared target, not fixed by the simulator.
        let plain = call(
            &c,
            "build_assembly_intent",
            json!({"inputs":inputs,"request":base}),
        )
        .await;
        let wrong = call(
            &c,
            "propose_overlap_assembly",
            json!({"inputs":inputs,"target":plain["target"],"recipe":designed["recipe"]}),
        )
        .await;
        assert_eq!(wrong["verification"]["status"], "Invalid");
    }
    assert_eq!(call(&c, "export_document", left).await, left_original);
    assert_eq!(call(&c, "export_document", right).await, right_original);
    c.cancel().await.unwrap();
    server.await.unwrap();
}
#[tokio::test]
async fn stale_source_bindings_invalid_junctions_and_oversized_trim_never_publish_intent() {
    let (c, server) = setup().await;
    let imported = call(
        &c,
        "import_document",
        json!({"text":">left\nATGACCCCGGGGTTTCCAGA\n>right\nATCCGTTAGCGACCTTGGAA\n"}),
    )
    .await;
    let left = named(&imported, "left");
    let right = named(&imported, "right");
    let inputs = json!([{"id":"l","document":left},{"id":"r","document":right}]);
    let before = call(&c, "list_documents", json!({})).await["documents"]
        .as_array()
        .unwrap()
        .len();
    let mut stale = inputs.clone();
    stale[0]["document"]["revision"] = json!("0".repeat(64));
    reject(
        &c,
        "build_edited_assembly_intent",
        json!({"inputs":stale,"base":base(false,false),"edits":edits()}),
    )
    .await;
    reject(
        &c,
        "design_edited_assembly_primers",
        json!({"inputs":stale,"request":planner(&base(false,false))}),
    )
    .await;
    for invalid in [
        json!([{"junction":2,"remove_left":0,"remove_right":0,"inserted":"GGA"}]),
        json!([{"junction":0,"remove_left":20,"remove_right":0,"inserted":"GGA"}]),
        json!([{"junction":0,"remove_left":0,"remove_right":0,"inserted":"!"}]),
    ] {
        reject(
            &c,
            "build_edited_assembly_intent",
            json!({"inputs":inputs,"base":base(false,false),"edits":invalid}),
        )
        .await;
    }
    let mut wrong_source = base(false, false);
    wrong_source["fragments"][1]["input_id"] = json!("not-bound");
    reject(
        &c,
        "build_edited_assembly_intent",
        json!({"inputs":inputs,"base":wrong_source,"edits":edits()}),
    )
    .await;
    assert_eq!(
        call(&c, "list_documents", json!({})).await["documents"]
            .as_array()
            .unwrap()
            .len(),
        before
    );
    c.cancel().await.unwrap();
    server.await.unwrap();
}
