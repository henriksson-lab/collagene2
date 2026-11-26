//! Actual source-bound tailing recipes with independently declared literal targets.
use base64::{Engine, engine::general_purpose::STANDARD};
use collagene_core::{analysis::Primer, verification::DesignPlan};
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
async fn artifact(c: &Client, id: &Value) -> String {
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
    text
}
async fn evidence(c: &Client, id: &Value) -> Value {
    serde_json::from_str(&artifact(c, id).await).unwrap()
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
async fn ta_gc_roundtrip_preserves_originals_target_and_explicit_nicks() {
    let (client, server) = setup().await;
    let tools = client.list_tools(Default::default()).await.unwrap().tools;
    assert_eq!(tools.len(), collagene_mcp::CURRENT_TOOL_COUNT);
    assert!(
        tools
            .iter()
            .any(|tool| tool.name == "propose_tailing_assembly")
    );
    let text = "LOCUS       vector 6 bp DNA linear SYN 01-JAN-2026\nFEATURES             Location/Qualifiers\nORIGIN\n        1 aaaccc\n//\nLOCUS       insert 6 bp DNA linear SYN 01-JAN-2026\nFEATURES             Location/Qualifiers\nORIGIN\n        1 gggttt\n//\nLOCUS       ta_target 14 bp DNA circular SYN 01-JAN-2026\nFEATURES             Location/Qualifiers\nORIGIN\n        1 aaaccctgggttta\n//\nLOCUS       gc_target 14 bp DNA circular SYN 01-JAN-2026\nFEATURES             Location/Qualifiers\nORIGIN\n        1 aaaccccgggtttg\n//\nLOCUS       wrong_target 14 bp DNA circular SYN 01-JAN-2026\nFEATURES             Location/Qualifiers\nORIGIN\n        1 aaaccctgggtttc\n//\n";
    let imported = call(&client, "import_document", json!({"text":text})).await;
    let recipe = json!({"name":"Original PCR products and explicit tails","fragments":[
        {"source":{"role":"vector","input_id":"vector","preparation":{"ExactPcr":{"forward":Primer::new("VF","AAA","").unwrap(),"reverse":Primer::new("VR","GGG","").unwrap()}},"reverse":false},"base":"T","treatment":"Phosphorylate"},
        {"source":{"role":"insert","input_id":"insert","preparation":{"ExactPcr":{"forward":Primer::new("IF","GGG","").unwrap(),"reverse":Primer::new("IR","AAA","").unwrap()}},"reverse":false},"base":"A","treatment":"Phosphorylate"}
    ],"circular":true,"allow_nicked":false});
    let args = json!({"inputs":[{"id":"vector","document":named(&imported,"vector")},{"id":"insert","document":named(&imported,"insert")}],"target":named(&imported,"ta_target"),"recipe":recipe});
    let result = call(&client, "propose_tailing_assembly", args.clone()).await;
    assert_eq!(result["verified"], false);
    assert_eq!(result["preparation_complete"], true);
    assert_eq!(
        result["verification"]["status"],
        "ConsistentUnderAssumptions"
    );
    assert_eq!(result["product_physical"]["top_length"], 14);
    assert_eq!(result["product_physical"]["nick_count"], 0);
    let full = evidence(&client, &result["artifact_id"]).await;
    let plan: DesignPlan = serde_json::from_value(full["compiled"]["plan"].clone()).unwrap();
    assert_eq!(plan.inputs[0].molecule.sequence.as_str(), "AAACCC");
    assert_eq!(plan.inputs[1].molecule.sequence.as_str(), "GGGTTT");
    assert_eq!(
        plan.steps
            .iter()
            .filter(|s| matches!(
                s,
                collagene_core::verification::Operation::Tail3Prime { .. }
            ))
            .count(),
        2
    );
    assert_eq!(full["recipe"], args["recipe"]);
    assert_eq!(full["junctions"].as_array().unwrap().len(), 2);
    let reference = json!({"plan_revision":result["plan_revision"]});
    assert_eq!(
        call(&client, "verify_plan", reference.clone()).await["status"],
        "ConsistentUnderAssumptions"
    );
    let html = call(&client, "export_protocol", reference.clone()).await;
    let text = artifact(&client, &html["artifact_id"]).await;
    assert!(text.contains("Append one T nucleotide"));
    assert!(text.contains("Append one A nucleotide"));
    assert!(text.contains("PCR"));
    let pdf = call(&client, "export_protocol_pdf", reference).await;
    assert_eq!(pdf["mime_type"], "application/pdf");
    let bytes = STANDARD
        .decode(artifact(&client, &pdf["artifact_id"]).await)
        .unwrap();
    assert!(bytes.starts_with(b"%PDF-"));
    assert!(bytes.ends_with(b"%%EOF") || bytes.ends_with(b"%%EOF\n"));
    assert_eq!(bytes.len() as u64, pdf["bytes"].as_u64().unwrap());
    let retry = call(&client, "propose_tailing_assembly", args.clone()).await;
    assert_eq!(retry["plan_revision"], result["plan_revision"]);
    assert_eq!(retry["artifact_id"], result["artifact_id"]);
    let mut gc = args.clone();
    gc["target"] = named(&imported, "gc_target");
    gc["recipe"]["fragments"][0]["base"] = json!("C");
    gc["recipe"]["fragments"][1]["base"] = json!("G");
    assert_eq!(
        call(&client, "propose_tailing_assembly", gc).await["verification"]["status"],
        "ConsistentUnderAssumptions"
    );
    let mut nicked = args.clone();
    nicked["recipe"]["fragments"][1]["treatment"] = Value::Null;
    assert_ne!(
        call(&client, "propose_tailing_assembly", nicked.clone()).await["verification"]["status"],
        "ConsistentUnderAssumptions"
    );
    nicked["recipe"]["allow_nicked"] = json!(true);
    let retained = call(&client, "propose_tailing_assembly", nicked).await;
    assert_eq!(
        retained["verification"]["status"],
        "ConsistentUnderAssumptions"
    );
    assert_eq!(retained["product_physical"]["nick_count"], 2);
    let mut wrong = args.clone();
    wrong["target"] = named(&imported, "wrong_target");
    assert_eq!(
        call(&client, "propose_tailing_assembly", wrong).await["verification"]["status"],
        "Invalid"
    );
    let mut stale = args;
    stale["inputs"][0]["document"]["revision"] = json!("stale");
    reject(&client, "propose_tailing_assembly", stale).await;
    assert_eq!(call(&client, "list_documents", json!({})).await["total"], 5);
    drop(client);
    server.abort();
}
