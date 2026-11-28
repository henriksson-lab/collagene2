//! Independent intended mutation targets precede all PCR/ligation reconstruction.
use collagene_mcp::server::CollageneServer;
use rmcp::{RoleClient, ServiceExt, model::CallToolRequestParams, service::RunningService};
use serde_json::{Value, json};
type Client = RunningService<RoleClient, ()>;
const SOURCE: &str = "AGTCCGATGCTACGTAGGACTTCAGACCTGTACGGATACCTTAGCGCA";
const REPLACED: &str = "AGTGGAGATGCTACGTAGGACTTCAGACCTGTACGGATACCTTAGCGCA";
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
    assert!(
        r.is_err() || r.unwrap().is_error == Some(true),
        "{name} unexpectedly accepted invalid input"
    );
}
async fn unverified(c: &Client, name: &str, args: Value) {
    if let Ok(r) = c
        .call_tool(
            CallToolRequestParams::new(name.to_owned())
                .with_arguments(args.as_object().unwrap().clone()),
        )
        .await
    {
        if r.is_error != Some(true) {
            let v = r.structured_content.unwrap();
            assert_ne!(
                v["computed_status"], "ConsistentUnderAssumptions",
                "reference-only intent became an original physical PCR substrate"
            );
        }
    }
}
async fn artifact(c: &Client, id: &Value) -> Value {
    let mut all = String::new();
    let mut start = 0;
    loop {
        let p = call(
            c,
            "read_artifact",
            json!({"artifact_id":id,"start":start,"max_bytes":16384}),
        )
        .await;
        all.push_str(p["text"].as_str().unwrap());
        if p["done"] == true {
            break;
        }
        start = p["next_start"].as_u64().unwrap();
    }
    serde_json::from_str(&all).unwrap()
}
async fn source(c: &Client) -> Value {
    let gb = format!(
        "LOCUS       mutation_source           48 bp    DNA     circular\nDEFINITION  Immutable physical source.\nFEATURES             Location/Qualifiers\n     misc_feature    10..18\n                     /label=\"retained feature\"\nORIGIN\n        1 {SOURCE}\n//\n"
    );
    let imported = call(c, "import_document", json!({"text":gb})).await;
    let rows = imported["documents"].as_array().unwrap();
    assert_eq!(rows.len(), 1);
    json!({"id":rows[0]["id"],"revision":rows[0]["revision"]})
}
async fn sequence(c: &Client, document: &Value, end: usize) -> String {
    call(
        c,
        "read_sequence_region",
        json!({"document":document,"start":0,"end":end}),
    )
    .await["sequence"]
        .as_str()
        .unwrap()
        .into()
}
fn stale(r: &Value) -> Value {
    json!({"id":r["id"],"revision":"0".repeat(64)})
}
#[tokio::test]
async fn literal_replace_insert_delete_targets_are_source_bound_and_not_preparation_answers() {
    let (c, s) = setup().await;
    assert_eq!(
        c.list_tools(Default::default()).await.unwrap().tools.len(),
        collagene_mcp::CURRENT_TOOL_COUNT
    );
    let original = source(&c).await;
    let original_export = call(&c, "export_document", original.clone()).await;
    let cases = [
        (3, 5, "GGA", REPLACED, 2),
        (
            3,
            3,
            "GGA",
            "AGTGGACCGATGCTACGTAGGACTTCAGACCTGTACGGATACCTTAGCGCA",
            0,
        ),
        (
            3,
            5,
            "",
            "AGTGATGCTACGTAGGACTTCAGACCTGTACGGATACCTTAGCGCA",
            2,
        ),
    ];
    let mut first = None;
    for (start, end, inserted, literal, deleted) in cases {
        let request = json!({"name":"Explicit mutation intent","edit_start":start,"edit_end":end,"inserted":inserted});
        let args = json!({"document":original,"request":request});
        let built = call(&c, "build_mutation_intent", args.clone()).await;
        assert_eq!(built["state"], "intended-reference");
        assert_eq!(built["verified"], false);
        assert_eq!(built["preparation_simulated"], false);
        assert_eq!(built["document"], original);
        assert_eq!(built["request"], request);
        assert!(built.get("plan_revision").is_none());
        assert_eq!(built["deleted_length"], deleted);
        assert_eq!(built["source_origin_shift"], 0);
        assert_eq!(built["target"].as_object().unwrap().len(), 2);
        assert_eq!(sequence(&c, &built["target"], literal.len()).await, literal);
        let stored = artifact(&c, &built["artifact_id"]).await;
        assert_eq!(stored["document"], original);
        assert_eq!(stored["report"]["source_revision"], original["revision"]);
        assert_eq!(stored["report"]["request"], request);
        let again = call(&c, "build_mutation_intent", args).await;
        assert_eq!(again["target"], built["target"]);
        assert_eq!(again["artifact_id"], built["artifact_id"]);
        if first.is_none() {
            first = Some(built);
        }
    }
    assert_eq!(sequence(&c, &original, SOURCE.len()).await, SOURCE);
    assert_eq!(
        call(&c, "export_document", original.clone()).await,
        original_export
    );
    let first = first.unwrap();
    let renamed=call(&c,"build_mutation_intent",json!({"document":original,"request":{"name":"Another intent name","edit_start":3,"edit_end":5,"inserted":"GGA"}})).await;
    assert_ne!(first["context_revision"], renamed["context_revision"]);
    assert_eq!(
        sequence(&c, &renamed["target"], REPLACED.len()).await,
        REPLACED
    );
    for request in [
        json!({"name":"bad","edit_start":48,"edit_end":49,"inserted":"A"}),
        json!({"name":"bad","edit_start":5,"edit_end":3,"inserted":"A"}),
        json!({"name":"bad","edit_start":3,"edit_end":5,"inserted":"!"}),
    ] {
        reject(
            &c,
            "build_mutation_intent",
            json!({"document":original,"request":request}),
        )
        .await;
    }
    reject(&c,"build_mutation_intent",json!({"document":stale(&original),"request":{"name":"stale","edit_start":3,"edit_end":5,"inserted":"GGA"}})).await;
    c.cancel().await.unwrap();
    s.await.unwrap();
}
#[tokio::test]
async fn frozen_literal_intent_then_independent_inverse_pcr_and_wrong_target_check() {
    let (c, s) = setup().await;
    let original = source(&c).await;
    let request =
        json!({"name":"Reviewed substitution intent","edit_start":3,"edit_end":5,"inserted":"GGA"});
    let intent = call(
        &c,
        "build_mutation_intent",
        json!({"document":original,"request":request}),
    )
    .await;
    // Literal sequence asserted BEFORE invoking any preparation simulator.
    assert_eq!(
        sequence(&c, &intent["target"], REPLACED.len()).await,
        REPLACED
    );
    let options = json!({"edit_start":3,"edit_end":5,"inserted":"GGA","forward_length":8,"reverse_length":8,"reverse_tail_bases":0});
    let designed = call(
        &c,
        "design_inverse_pcr",
        json!({"document":original,"target":intent["target"],"origin":null,"origin_workspace":null,"options":options}),
    )
    .await;
    assert_eq!(designed["computed_status"], "ConsistentUnderAssumptions");
    assert_eq!(designed["source_revision"], original["revision"]);
    assert_eq!(designed["target_revision"], intent["target"]["revision"]);
    let checked = call(
        &c,
        "verify_plan",
        json!({"plan_revision":designed["plan_revision"]}),
    )
    .await;
    assert_eq!(checked["status"], "ConsistentUnderAssumptions");
    assert_eq!(checked["products"].as_array().unwrap().len(), 3);
    let full = artifact(&c, &designed["artifact_id"]).await;
    assert_eq!(
        full["inverse_pcr"]["plan"]["steps"][1]["TreatEnds"]["treatment"],
        "Phosphorylate"
    );
    let changed=call(&c,"build_mutation_intent",json!({"document":original,"request":{"name":"Different edit","edit_start":3,"edit_end":5,"inserted":"TTA"}})).await;
    assert_ne!(changed["context_revision"], intent["context_revision"]);
    let wrong = call(
        &c,
        "design_inverse_pcr",
        json!({"document":original,"target":changed["target"],"origin":null,"origin_workspace":null,"options":options}),
    )
    .await;
    assert_eq!(wrong["computed_status"], "Invalid");
    reject(
        &c,
        "design_inverse_pcr",
        json!({"document":original,"target":stale(&intent["target"]),"origin":null,"origin_workspace":null,"options":options}),
    )
    .await;
    unverified(
        &c,
        "design_inverse_pcr",
        json!({"document":intent["target"],"target":original,"origin":null,"origin_workspace":null,"options":options}),
    )
    .await;
    assert_eq!(sequence(&c, &original, SOURCE.len()).await, SOURCE);
    c.cancel().await.unwrap();
    s.await.unwrap();
}
