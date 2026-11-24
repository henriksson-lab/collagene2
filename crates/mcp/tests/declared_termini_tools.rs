//! Explicit original duplex declarations do not restore imported physical authority.
use collagene_core::verification::{DesignPlan, ENGINE_VERSION};
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
async fn text_artifact(client: &Client, id: &Value) -> String {
    let mut text = String::new();
    let mut start = 0;
    loop {
        let page = call(
            client,
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
    text
}
fn terminus() -> Value {
    json!({"five_prime":"Phosphate","three_prime":"Hydroxyl"})
}
fn declaration(id: &str) -> Value {
    json!({"DeclareSourceTermini":{"template":{"Input":id},"left":terminus(),"right":terminus()}})
}
#[tokio::test]
async fn explicit_original_source_termini_enable_declared_ligation_without_mutating_sources() {
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
        json!({"text":">left\nAAAA\n>right\nCCCC\n>independent_target\nAAAACCCC\n"}),
    )
    .await;
    let l = named(&imported, "left");
    let r = named(&imported, "right");
    let target = named(&imported, "independent_target");
    let before = call(&client, "export_document", l.clone()).await;
    let inputs = json!([{"id":"l","document":l},{"id":"r","document":r}]);
    let recipe = json!({"name":"Explicit original duplex assumptions","fragments":[
        {"source":{"role":"left","input_id":"l","preparation":{"DeclaredDuplex":{"left":terminus(),"right":terminus()}},"reverse":false},"treatment":null},
        {"source":{"role":"right","input_id":"r","preparation":{"DeclaredDuplex":{"left":terminus(),"right":terminus()}},"reverse":false},"treatment":null}
    ],"circular":false,"allow_nicked":false});
    let args = json!({"inputs":inputs,"target":target,"recipe":recipe});
    let result = call(&client, "propose_ligation_assembly", args.clone()).await;
    assert_eq!(result["verified"], false);
    assert_eq!(
        result["verification"]["status"],
        "ConsistentUnderAssumptions"
    );
    assert_eq!(result["verification"]["engine_version"], ENGINE_VERSION);
    assert_eq!(result["product_physical"]["top_length"], 8);
    assert_eq!(result["product_physical"]["nick_count"], 0);
    let full = evidence(&client, &result["artifact_id"]).await;
    let plan: DesignPlan = serde_json::from_value(full["compiled"]["plan"].clone()).unwrap();
    assert_eq!(plan.steps.len(), 3);
    assert_eq!(plan.inputs[0].molecule.sequence.as_str(), "AAAA");
    assert!(
        full["compiled"]["plan"]["steps"][0]
            .get("DeclareSourceTermini")
            .is_some()
    );
    let plan_ref = json!({"plan_revision":result["plan_revision"]});
    let verified = call(&client, "verify_plan", plan_ref.clone()).await;
    assert_eq!(verified["status"], "ConsistentUnderAssumptions");
    assert_eq!(verified["engine_version"], ENGINE_VERSION);
    let html = call(&client, "export_protocol", plan_ref.clone()).await;
    assert_eq!(html["plan_revision"], result["plan_revision"]);
    let html = text_artifact(&client, &html["artifact_id"]).await;
    assert!(html.contains(ENGINE_VERSION));
    assert!(html.contains("Phosphate"));
    let pdf = call(&client, "export_protocol_pdf", plan_ref).await;
    assert_eq!(pdf["plan_revision"], result["plan_revision"]);
    assert_eq!(pdf["engine_version"], ENGINE_VERSION);
    assert_eq!(pdf["status"], "ConsistentUnderAssumptions");
    assert_eq!(
        call(&client, "export_document", l.clone()).await["artifact_id"],
        before["artifact_id"]
    );
    let mut omitted = args.clone();
    for fragment in omitted["recipe"]["fragments"].as_array_mut().unwrap() {
        fragment["source"]["preparation"] = json!("WholeLinear");
    }
    assert_eq!(
        call(&client, "propose_ligation_assembly", omitted).await["verification"]["status"],
        "Incomplete"
    );
    let mut unknown = args.clone();
    unknown["recipe"]["fragments"][0]["source"]["preparation"]["DeclaredDuplex"]["right"]["three_prime"] =
        json!("Unknown");
    assert_eq!(
        call(&client, "propose_ligation_assembly", unknown).await["verification"]["status"],
        "Incomplete"
    );
    let mut blocked = args;
    blocked["recipe"]["fragments"][0]["source"]["preparation"]["DeclaredDuplex"]["right"]["three_prime"] =
        json!("Blocked");
    assert_ne!(
        call(&client, "propose_ligation_assembly", blocked).await["verification"]["status"],
        "ConsistentUnderAssumptions"
    );
    // The direct typed tool exposes the same explicit operation; no recipe-only
    // bypass and no ability to redeclare a prior Step as an intact original.
    let direct=call(&client,"propose_plan",json!({"name":"Direct terminal declarations","target":target,"inputs":inputs,"steps":[declaration("l"),declaration("r"),{"Ligate":{"fragments":[{"source":{"Step":0},"reverse":false},{"source":{"Step":1},"reverse":false}],"circular":false,"allow_nicked":false}}]})).await;
    assert_eq!(
        call(
            &client,
            "verify_plan",
            json!({"plan_revision":direct["plan_revision"]})
        )
        .await["status"],
        "ConsistentUnderAssumptions"
    );
    let step=call(&client,"propose_plan",json!({"name":"Step cannot become original","target":l,"inputs":inputs,"steps":[declaration("l"),{"DeclareSourceTermini":{"template":{"Step":0},"left":terminus(),"right":terminus()}}]})).await;
    assert_eq!(
        call(
            &client,
            "verify_plan",
            json!({"plan_revision":step["plan_revision"]})
        )
        .await["status"],
        "Incomplete"
    );
    assert_eq!(call(&client, "list_documents", json!({})).await["total"], 3);
    drop(client);
    server.abort();
}
#[tokio::test]
async fn single_stranded_circular_and_imported_physical_sources_cannot_be_redeclared_as_intact_duplexes()
 {
    let (client, server) = setup().await;
    let text = "LOCUS       ss_source                    8 bp    ss-DNA  linear   SYN 01-JAN-2026\nFEATURES             Location/Qualifiers\nORIGIN\n        1 aaaacccc\n//\nLOCUS       circle                       8 bp    DNA     circular SYN 01-JAN-2026\nFEATURES             Location/Qualifiers\nORIGIN\n        1 aaaacccc\n//\nLOCUS       imported_physical            8 bp    DNA     linear   SYN 01-JAN-2026\nCOMMENT     Collagene physical product v1: unverified sequence-only evidence\nFEATURES             Location/Qualifiers\nORIGIN\n        1 aaaacccc\n//\nLOCUS       reference_only               8 bp    DNA     linear   SYN 01-JAN-2026\nCOMMENT     Collagene reference-only v1: unverified sequence-only evidence\nFEATURES             Location/Qualifiers\nORIGIN\n        1 aaaacccc\n//\n";
    let imported = call(&client, "import_document", json!({"text":text})).await;
    let target = call(
        &client,
        "import_document",
        json!({"text":">independent_target\nAAAACCCC\n"}),
    )
    .await;
    let target = named(&target, "independent_target");
    for name in ["ss_source", "circle", "imported_physical", "reference_only"] {
        let proposal=call(&client,"propose_plan",json!({"name":format!("Reject inherited {name}"),"target":target,"inputs":[{"id":"s","document":named(&imported,name)}],"steps":[declaration("s")]})).await;
        let report = call(
            &client,
            "verify_plan",
            json!({"plan_revision":proposal["plan_revision"]}),
        )
        .await;
        assert_eq!(report["status"], "Incomplete", "{name}: {report}");
        assert!(
            report["products"].as_array().unwrap().is_empty(),
            "{name} manufactured an intermediate"
        );
    }
    drop(client);
    server.abort();
}
