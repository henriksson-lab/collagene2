//! Literal PCR/digest homology fixture, with an independently supplied target.
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
#[tokio::test]
async fn actual_pcr_and_digest_ordered_recipe_binds_independent_target_and_saved_protocol() {
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
    // EcoRV cuts GAT|ATC. Fragment1 is ATCGTCGATGCTAGCGGTAC;
    // the exact PCR left product is TTGCAATCGTCGATGCTAGC. Shared15nt
    // ATCGTCGATGCTAGC is represented once in the supplied target25nt.
    let imported=call(&client,"import_document",json!({"text":">left_template\nTTGCAATCGTCGATGCTAGC\n>right_digest_source\nGATATCGTCGATGCTAGCGGTAC\n>independent_target\nTTGCAATCGTCGATGCTAGCGGTAC\n>wrong_target\nTTGCAATCGTCGATGCTAGCGGTAT\n"})).await;
    let refs = |name: &str| {
        let row = imported["documents"]
            .as_array()
            .unwrap()
            .iter()
            .find(|d| d["name"] == name)
            .unwrap();
        json!({"id":row["id"],"revision":row["revision"]})
    };
    let inputs = json!([{"id":"template","document":refs("left_template")},{"id":"vector","document":refs("right_digest_source")}]);
    let recipe = json!({"name":"Actual PCR + EcoRV ordered assembly","fragments":[
        {"role":"left","input_id":"template","preparation":{"ExactPcr":{"forward":Primer::new("actual forward","TTGCAATCGT","").unwrap(),"reverse":Primer::new("actual reverse","GCTAGCATCG","").unwrap()}},"reverse":false},
        {"role":"right","input_id":"vector","preparation":{"Digest":{"enzymes":["EcoRV"],"fragment_index":1,"repair":null}},"reverse":false}
    ],"overlaps":[15],"circular":false,"options":{"method":"HiFi","scan_repeated_overlaps":false,"max_comparisons":100000,"max_occurrences":100}});
    let args = json!({"inputs":inputs,"target":refs("independent_target"),"recipe":recipe});
    let result = call(&client, "propose_overlap_assembly", args.clone()).await;
    assert_eq!(result["state"], "candidate");
    assert_eq!(result["verified"], false);
    assert_eq!(result["preparation_complete"], true);
    assert_eq!(
        result["verification"]["status"],
        "ConsistentUnderAssumptions"
    );
    let full = evidence(&client, &result["artifact_id"]).await;
    assert_eq!(full["compiled"]["rows"].as_array().unwrap().len(), 2);
    let plan: DesignPlan = serde_json::from_value(full["compiled"]["plan"].clone()).unwrap();
    assert_eq!(plan.inputs.len(), 2);
    assert_eq!(plan.steps.len(), 3);
    assert_eq!(
        plan.inputs[1].molecule.sequence.as_str(),
        "GATATCGTCGATGCTAGCGGTAC"
    );
    let plan_ref = json!({"plan_revision":result["plan_revision"]});
    assert_eq!(
        call(&client, "verify_plan", plan_ref.clone()).await["status"],
        "ConsistentUnderAssumptions"
    );
    assert!(call(&client, "export_protocol", plan_ref.clone()).await["artifact_id"].is_string());
    assert!(call(&client, "export_protocol_pdf", plan_ref).await["artifact_id"].is_string());
    assert_eq!(
        call(&client, "propose_overlap_assembly", args.clone()).await["plan_revision"],
        result["plan_revision"]
    );
    let mut wrong = args.clone();
    wrong["target"] = refs("wrong_target");
    assert_eq!(
        call(&client, "propose_overlap_assembly", wrong).await["verification"]["status"],
        "Invalid"
    );
    let mut reversed = args.clone();
    reversed["recipe"]["fragments"][1]["reverse"] = json!(true);
    assert_ne!(
        call(&client, "propose_overlap_assembly", reversed).await["verification"]["status"],
        "ConsistentUnderAssumptions"
    );
    let mut duplicate_role = args.clone();
    duplicate_role["recipe"]["fragments"][1]["role"] = json!("left");
    reject(&client, "propose_overlap_assembly", duplicate_role).await;
    let mut stale = args.clone();
    stale["inputs"][0]["document"]["revision"] = json!("stale");
    reject(&client, "propose_overlap_assembly", stale).await;
    let mut swapped = args.clone();
    swapped["inputs"][0]["document"] = refs("right_digest_source");
    assert_ne!(
        call(&client, "propose_overlap_assembly", swapped).await["verification"]["status"],
        "ConsistentUnderAssumptions"
    );
    assert_eq!(call(&client, "list_documents", json!({})).await["total"], 4);
    drop(client);
    server.abort();
}
