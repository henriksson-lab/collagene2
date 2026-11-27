//! Selective junction proposals preserve retained manual PCR and catalog digest rows exactly.
use collagene_core::{analysis::Primer, catalog_digest, primers::DesignOptions};
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
const TARGET: &str = "ATGACCCCGGGGTTTCCAGAGCGCATAGTCGACTAGCTTGCAATCATCCGTTAGCGACCTTGGAAGAT";
fn base() -> Value {
    json!({"name":"Selective three-source assembly","fragments":[{"role":"left","input_id":"a","start":0,"end":20,"reverse":false},{"role":"manual","input_id":"b","start":0,"end":20,"reverse":false},{"role":"existing-cut","input_id":"c","start":3,"end":29,"reverse":false}],"circular":false})
}
fn edits() -> Value {
    json!([{"junction":0,"remove_left":0,"remove_right":0,"inserted":"GC"}])
}
fn request() -> Value {
    let catalog = catalog_digest::builtin_catalog().unwrap();
    let pin = catalog_digest::select(&catalog, &["ecorv".into()]).unwrap();
    let manual = json!({"generation":null,"role":"manual","input_id":"b","reverse":false,"preparation":{"ExactPcr":{"forward":Primer::new("Already reviewed named BF","GCATAGTC","").unwrap(),"reverse":Primer::new("Already reviewed named BR","TGCAAGCT","AGGTCGCTAACGGATGAT").unwrap()}}});
    let digest = json!({"generation":null,"role":"existing-cut","input_id":"c","reverse":false,"preparation":{"CatalogDigest":{"selection":pin,"fragment_index":1,"repair":null}}});
    json!({"base":base(),"edits":edits(),"rows":[{"DesignPcr":{"role":"left","input_id":"a","reverse":false}},{"Retain":manual},{"Retain":digest}],"overlaps":[15,18],"primer_options":options(),"profile_options":{"method":"HiFi","scan_repeated_overlaps":false,"max_comparisons":100000,"max_occurrences":100}})
}
async fn inputs(c: &Client) -> (Value, Vec<Value>) {
    let imported = call(
        c,
        "import_document",
        json!({"text":">A
ATGACCCCGGGGTTTCCAGA
>B
GCATAGTCGACTAGCTTGCA
>C
GATATCATCCGTTAGCGACCTTGGAAGATATC
"}),
    )
    .await;
    let docs = vec![
        named(&imported, "A"),
        named(&imported, "B"),
        named(&imported, "C"),
    ];
    let bound = json!([{"id":"a","document":docs[0]},{"id":"b","document":docs[1]},{"id":"c","document":docs[2]}]);
    (bound, docs)
}
#[tokio::test]
async fn replan_only_selected_pcr_preserves_named_oligos_catalog_pin_and_variable_overlaps() {
    let (c, server) = setup().await;
    assert_eq!(
        c.list_tools(Default::default()).await.unwrap().tools.len(),
        collagene_mcp::CURRENT_TOOL_COUNT
    );
    let (bound, docs) = inputs(&c).await;
    let mut originals = Vec::new();
    for doc in &docs {
        originals.push(call(&c, "export_document", doc.clone()).await);
    }
    let intent = call(
        &c,
        "build_edited_assembly_intent",
        json!({"inputs":bound,"base":base(),"edits":edits()}),
    )
    .await;
    assert_eq!(intent["state"], "intended-reference");
    assert_eq!(intent["verified"], false);
    assert_eq!(TARGET.len(), 68);
    // Fix literal expected DNA before any selected-row primer design or preparation.
    assert_eq!(
        call(
            &c,
            "read_sequence_region",
            json!({"document":intent["target"],"start":0,"end":68})
        )
        .await["sequence"],
        TARGET
    );
    let request = request();
    let proposed = call(
        &c,
        "design_selective_assembly_primers",
        json!({"inputs":bound,"request":request}),
    )
    .await;
    assert_eq!(proposed["state"], "advisory");
    assert_eq!(proposed["verified"], false);
    assert_eq!(
        proposed["intent_context_revision"],
        intent["context_revision"]
    );
    assert_eq!(proposed["target_revision"], intent["target"]["revision"]);
    assert_eq!(proposed["recipe"]["overlaps"], json!([15, 18]));
    assert_eq!(proposed["retained_rows"], json!([1, 2]));
    assert_eq!(
        proposed["recipe"]["fragments"][1],
        request["rows"][1]["Retain"]
    );
    assert_eq!(
        proposed["recipe"]["fragments"][2],
        request["rows"][2]["Retain"]
    );
    assert_eq!(
        proposed["recipe"]["fragments"][0]["preparation"]["ExactPcr"]["reverse"]["tail"],
        "GCTAGTCGACTATGCGC"
    );
    let ordered = proposed["ordered_oligos"].as_array().unwrap();
    assert_eq!(ordered.len(), 1);
    assert_eq!(ordered[0]["role"], "left");
    let full = evidence(&c, &proposed["artifact_id"]).await;
    assert_eq!(full["inputs"], bound);
    assert_eq!(full["report"]["request"], request);
    assert_eq!(
        full["report"]["intent"]["target"]["sequence"]["bases"],
        TARGET
    );
    assert!(full["report"]["designs"][0].is_object());
    assert!(full["report"]["designs"][1].is_null());
    assert!(full["report"]["designs"][2].is_null());
    let repeated = call(
        &c,
        "design_selective_assembly_primers",
        json!({"inputs":bound,"request":request}),
    )
    .await;
    assert_eq!(repeated["artifact_id"], proposed["artifact_id"]);
    let assembled = call(
        &c,
        "propose_overlap_assembly",
        json!({"inputs":bound,"target":intent["target"],"recipe":proposed["recipe"]}),
    )
    .await;
    assert_eq!(
        assembled["verification"]["status"],
        "ConsistentUnderAssumptions"
    );
    let fresh = call(
        &c,
        "verify_plan",
        json!({"plan_revision":assembled["plan_revision"]}),
    )
    .await;
    assert_eq!(fresh["status"], "ConsistentUnderAssumptions");
    let unedited = call(
        &c,
        "build_assembly_intent",
        json!({"inputs":bound,"request":base()}),
    )
    .await;
    let mismatch = call(
        &c,
        "propose_overlap_assembly",
        json!({"inputs":bound,"target":unedited["target"],"recipe":proposed["recipe"]}),
    )
    .await;
    assert_eq!(mismatch["verification"]["status"], "Invalid");
    for (doc, original) in docs.iter().zip(originals) {
        assert_eq!(call(&c, "export_document", doc.clone()).await, original);
    }
    c.cancel().await.unwrap();
    server.await.unwrap();
}
#[tokio::test]
async fn role_order_stale_bindings_and_invalid_retained_catalog_preparations_are_rejected() {
    let (c, server) = setup().await;
    let (bound, _) = inputs(&c).await;
    let request = request();
    let count = call(&c, "list_documents", json!({})).await["documents"]
        .as_array()
        .unwrap()
        .len();
    let mut stale = bound.clone();
    stale[1]["document"]["revision"] = json!("0".repeat(64));
    reject(
        &c,
        "design_selective_assembly_primers",
        json!({"inputs":stale,"request":request}),
    )
    .await;
    let mut wrong_role = request.clone();
    wrong_role["rows"][0]["DesignPcr"]["role"] = json!("not-left");
    let mut wrong_order = request.clone();
    wrong_order["rows"].as_array_mut().unwrap().swap(1, 2);
    let mut wrong_pin = request.clone();
    wrong_pin["rows"][2]["Retain"]["preparation"]["CatalogDigest"]["selection"]["catalog_content_hash"] =
        json!("0".repeat(64));
    for invalid in [wrong_role, wrong_order, wrong_pin] {
        reject(
            &c,
            "design_selective_assembly_primers",
            json!({"inputs":bound,"request":invalid}),
        )
        .await;
    }
    assert_eq!(
        call(&c, "list_documents", json!({})).await["documents"]
            .as_array()
            .unwrap()
            .len(),
        count
    );
    c.cancel().await.unwrap();
    server.await.unwrap();
}
