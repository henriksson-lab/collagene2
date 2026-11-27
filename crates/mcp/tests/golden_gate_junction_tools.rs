//! Actual Type IIS geometry with explicit edited intent and selective primer replacement.
use collagene_core::{analysis::Primer, catalog_digest};
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
const A: &str = "ATGACCCCGGGGTTTCCAGA";
const B: &str = "ATCCGTTAGCGACCTTGGAA";
const TARGET: &str = "ATGACCCCGGGGTTTCCAGAGCATCCGTTAGCGACCTTGGAA";
fn base() -> Value {
    json!({"name":"Selective circular Type IIS junction","fragments":[{"role":"left","input_id":"a","start":0,"end":20,"reverse":false},{"role":"manual","input_id":"b","start":0,"end":20,"reverse":false}],"circular":true})
}
fn edits() -> Value {
    json!([{"junction":0,"remove_left":0,"remove_right":0,"inserted":"GC"}])
}
fn request() -> Value {
    let catalog = catalog_digest::builtin_catalog().unwrap();
    let selection = catalog_digest::select(&catalog, &["bsai".into()]).unwrap();
    let manual = json!({"source":{"generation":null,"role":"manual","input_id":"b","reverse":false,"preparation":{"PcrCatalogDigest":{"forward":Primer::new("Retained manual BF","ATCCGTTA","ATGCGGTCTCA").unwrap(),"reverse":Primer::new("Retained manual BR","TTCCAAGG","ATGCGGTCTCATCAT").unwrap(),"selection":selection,"fragment_index":1,"repair":null}}},"treatment":null});
    json!({"base":base(),"edits":edits(),"rows":[{"DesignPcr":{"role":"left","input_id":"a","reverse":false,"forward_length":8,"reverse_length":8}},{"Retain":manual}],"enzyme":"BsaI","clamp":"ATGC","spacer":"A","terminal_right_overhang":null,"allow_nicked":false})
}
async fn sources(c: &Client) -> (Value, Vec<Value>) {
    let imported = call(
        c,
        "import_document",
        json!({"text":format!(">A\n{A}\n>B\n{B}\n")}),
    )
    .await;
    let a = named(&imported, "A");
    let b = named(&imported, "B");
    (
        json!([{"id":"a","document":a},{"id":"b","document":b}]),
        vec![a, b],
    )
}
#[tokio::test]
async fn inserted_junction_replans_only_selected_type_iis_pcr_and_preserves_manual_closure() {
    let (c, server) = setup().await;
    assert_eq!(
        c.list_tools(Default::default()).await.unwrap().tools.len(),
        collagene_mcp::CURRENT_TOOL_COUNT
    );
    let (inputs, docs) = sources(&c).await;
    let mut originals = Vec::new();
    for doc in &docs {
        originals.push(call(&c, "export_document", doc.clone()).await);
    }
    let expected = call(
        &c,
        "build_edited_assembly_intent",
        json!({"inputs":inputs,"base":base(),"edits":edits()}),
    )
    .await;
    assert_eq!(expected["state"], "intended-reference");
    assert_eq!(expected["verified"], false);
    assert_eq!(expected["preparation_simulated"], false);
    assert_eq!(TARGET.len(), 42);
    assert_eq!(
        call(
            &c,
            "read_sequence_region",
            json!({"document":expected["target"],"start":0,"end":42})
        )
        .await["sequence"],
        TARGET
    );
    let count = call(&c, "list_documents", json!({})).await["documents"]
        .as_array()
        .unwrap()
        .len();
    let request = request();
    let proposed = call(
        &c,
        "design_golden_gate_junction_primers",
        json!({"inputs":inputs,"request":request}),
    )
    .await;
    assert_eq!(proposed["state"], "advisory");
    assert_eq!(proposed["verified"], false);
    assert_eq!(
        proposed["intent_context_revision"],
        expected["context_revision"]
    );
    assert_eq!(proposed["target_revision"], expected["target"]["revision"]);
    assert_eq!(proposed["retained_rows"], json!([1]));
    assert_eq!(
        proposed["recipe"]["fragments"][1],
        request["rows"][1]["Retain"]
    );
    assert_eq!(proposed["recipe"]["circular"], true);
    assert_eq!(proposed["recipe"]["allow_nicked"], false);
    let selected = &proposed["recipe"]["fragments"][0]["source"]["preparation"]["PcrCatalogDigest"];
    assert_eq!(selected["forward"]["tail"], "ATGCGGTCTCA");
    assert_eq!(selected["reverse"]["tail"], "ATGCGGTCTCAGGATGC");
    assert_eq!(selected["fragment_index"], 1);
    assert!(selected["repair"].is_null());
    assert!(proposed["recipe"]["fragments"][0]["treatment"].is_null());
    assert_eq!(proposed["ordered_oligos"].as_array().unwrap().len(), 1);
    let full = evidence(&c, &proposed["artifact_id"]).await;
    assert_eq!(full["inputs"], inputs);
    assert_eq!(full["report"]["request"], request);
    assert_eq!(
        full["report"]["intent"]["target"]["sequence"]["bases"],
        TARGET
    );
    let design = &full["report"]["designs"][0];
    assert!(full["report"]["designs"][1].is_null());
    assert_eq!(design["reverse_tail"]["clamp"], "ATGC");
    assert_eq!(design["reverse_tail"]["motif"], "GGTCTC");
    assert_eq!(design["reverse_tail"]["spacer"], "A");
    assert_eq!(design["reverse_tail"]["fusion"], "GGAT");
    assert_eq!(design["reverse_tail"]["inserted"], "GC");
    assert_eq!(design["oriented_top"], "ATGACCCCGGGGTTTCCAGAGC");
    assert_eq!(design["left"]["overhang"], "ATGA");
    assert_eq!(design["right"]["overhang"], "GGAT");
    assert_eq!(
        call(&c, "list_documents", json!({})).await["documents"]
            .as_array()
            .unwrap()
            .len(),
        count
    );
    let repeated = call(
        &c,
        "design_golden_gate_junction_primers",
        json!({"inputs":inputs,"request":request}),
    )
    .await;
    assert_eq!(repeated["artifact_id"], proposed["artifact_id"]);
    // Independent intent is already fixed. Only this separate step reconstructs actual chemistry.
    let assembled = call(
        &c,
        "propose_ligation_assembly",
        json!({"inputs":inputs,"target":expected["target"],"recipe":proposed["recipe"]}),
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
    let wrong_target = call(
        &c,
        "build_assembly_intent",
        json!({"inputs":inputs,"request":base()}),
    )
    .await;
    let wrong = call(
        &c,
        "propose_ligation_assembly",
        json!({"inputs":inputs,"target":wrong_target["target"],"recipe":proposed["recipe"]}),
    )
    .await;
    assert_eq!(wrong["verification"]["status"], "Invalid");
    for (doc, original) in docs.iter().zip(originals) {
        assert_eq!(call(&c, "export_document", doc.clone()).await, original);
    }
    c.cancel().await.unwrap();
    server.await.unwrap();
}
#[tokio::test]
async fn stale_sources_row_role_aliases_order_and_invalid_catalog_pins_are_not_reinterpreted() {
    let (c, server) = setup().await;
    let (inputs, _) = sources(&c).await;
    let request = request();
    let before = call(&c, "list_documents", json!({})).await["documents"]
        .as_array()
        .unwrap()
        .len();
    let mut stale = inputs.clone();
    stale[1]["document"]["revision"] = json!("0".repeat(64));
    reject(
        &c,
        "design_golden_gate_junction_primers",
        json!({"inputs":stale,"request":request}),
    )
    .await;
    let mut wrong_role = request.clone();
    wrong_role["rows"][1]["Retain"]["source"]["role"] = json!("alias-for-manual");
    let mut wrong_order = request.clone();
    wrong_order["rows"].as_array_mut().unwrap().swap(0, 1);
    let mut wrong_pin = request.clone();
    wrong_pin["rows"][1]["Retain"]["source"]["preparation"]["PcrCatalogDigest"]["selection"]["catalog_content_hash"] =
        json!("0".repeat(64));
    let mut alias = request.clone();
    alias["rows"][1]["Retain"]["source"]["preparation"]["PcrCatalogDigest"]["selection"]["enzymes"]
        [0]["id"] = json!("BsaI");
    for invalid in [wrong_role, wrong_order, wrong_pin, alias] {
        reject(
            &c,
            "design_golden_gate_junction_primers",
            json!({"inputs":inputs,"request":invalid}),
        )
        .await;
    }
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
