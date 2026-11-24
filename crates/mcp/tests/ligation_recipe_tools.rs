//! Explicit physical preparations and ligation targets, specified independently.
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
async fn actual_two_pcr_fragments_require_explicit_phosphorylation_and_independent_target() {
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
    let imported=call(&client,"import_document",json!({"text":">left_template\nACGACCCCGGGGTAAC\n>right_template\nTATGAGCTCCAGGTCA\n>literal_target\nACGACCCCGGGGTAACTATGAGCTCCAGGTCA\n>wrong_target\nACGACCCCGGGGTAACTATGAGCTCCAGGTCT\n"})).await;
    let left = named(&imported, "left_template");
    let right = named(&imported, "right_template");
    let recipe = json!({"name":"Two actual PCR products, explicit phosphorylation","fragments":[
        {"source":{"role":"left","input_id":"left_template","preparation":{"ExactPcr":{"forward":Primer::new("ordered LF","ACGACCCC","").unwrap(),"reverse":Primer::new("ordered LR","GTTACCCC","").unwrap()}},"reverse":false},"treatment":"Phosphorylate"},
        {"source":{"role":"right","input_id":"right_template","preparation":{"ExactPcr":{"forward":Primer::new("ordered RF","TATGAGCT","").unwrap(),"reverse":Primer::new("ordered RR","TGACCTGG","").unwrap()}},"reverse":false},"treatment":"Phosphorylate"}
    ],"circular":false,"allow_nicked":false});
    let args = json!({"inputs":[{"id":"left_template","document":left},{"id":"right_template","document":right}],"target":named(&imported,"literal_target"),"recipe":recipe});
    let result = call(&client, "propose_ligation_assembly", args.clone()).await;
    assert_eq!(result["verified"], false);
    assert_eq!(result["preparation_complete"], true);
    assert_eq!(
        result["verification"]["status"],
        "ConsistentUnderAssumptions"
    );
    assert_eq!(result["product_physical"]["top_length"], 32);
    assert_eq!(result["product_physical"]["nick_count"], 0);
    let full = evidence(&client, &result["artifact_id"]).await;
    let plan: DesignPlan = serde_json::from_value(full["compiled"]["plan"].clone()).unwrap();
    assert_eq!(plan.inputs.len(), 2);
    assert_eq!(plan.steps.len(), 5);
    assert_eq!(full["junctions"].as_array().unwrap().len(), 1);
    assert_eq!(full["junctions"][0]["upstream_end"]["polarity"], "Blunt");
    assert_eq!(
        full["junctions"][0]["downstream_end"]["chemistry"]["five_prime"],
        "Phosphate"
    );
    let plan_ref = json!({"plan_revision":result["plan_revision"]});
    assert_eq!(
        call(&client, "verify_plan", plan_ref.clone()).await["status"],
        "ConsistentUnderAssumptions"
    );
    assert!(call(&client, "export_protocol", plan_ref.clone()).await["artifact_id"].is_string());
    assert!(call(&client, "export_protocol_pdf", plan_ref).await["artifact_id"].is_string());
    assert_eq!(
        call(&client, "propose_ligation_assembly", args.clone()).await["plan_revision"],
        result["plan_revision"]
    );
    let mut untreated = args.clone();
    untreated["recipe"]["fragments"][1]["treatment"] = Value::Null;
    let blocked = call(&client, "propose_ligation_assembly", untreated.clone()).await;
    assert_ne!(
        blocked["verification"]["status"],
        "ConsistentUnderAssumptions"
    );
    // Explicit permission retains the unsealed bond; it does not repair it.
    untreated["recipe"]["allow_nicked"] = json!(true);
    let nicked = call(&client, "propose_ligation_assembly", untreated).await;
    assert_eq!(
        nicked["verification"]["status"],
        "ConsistentUnderAssumptions"
    );
    assert_eq!(nicked["product_physical"]["nick_count"], 1);
    let mut wrong = args.clone();
    wrong["target"] = named(&imported, "wrong_target");
    assert_eq!(
        call(&client, "propose_ligation_assembly", wrong).await["verification"]["status"],
        "Invalid"
    );
    let mut reverse = args.clone();
    reverse["recipe"]["fragments"][1]["source"]["reverse"] = json!(true);
    assert_ne!(
        call(&client, "propose_ligation_assembly", reverse).await["verification"]["status"],
        "ConsistentUnderAssumptions"
    );
    let mut stale = args;
    stale["inputs"][0]["document"]["revision"] = json!("stale");
    reject(&client, "propose_ligation_assembly", stale).await;
    assert_eq!(call(&client, "list_documents", json!({})).await["total"], 4);
    drop(client);
    server.abort();
}
#[tokio::test]
async fn three_actual_circular_ecori_digests_keep_compatible_cohesive_ends_and_all_junctions() {
    let (client, server) = setup().await;
    // Each EcoRI top cut G|AATTC defines the indicated origin. Literal top
    // segments are concatenated, with all three GAATTC junctions reconstructed.
    let sequences = [
        ("vector", "AAAAGAATTCTTTT"),
        ("insert1", "CCCCGAATTCGGGG"),
        ("insert2", "ACACGAATTCTGTG"),
        ("target", "AATTCTTTTAAAAGAATTCGGGGCCCCGAATTCTGTGACACG"),
    ];
    let text=sequences.iter().map(|(name,sequence)|format!("LOCUS       {name:<20} {} bp    DNA     circular SYN 01-JAN-2026\nFEATURES             Location/Qualifiers\nORIGIN\n        1 {}\n//\n",sequence.len(),sequence.to_ascii_lowercase())).collect::<String>();
    let imported = call(&client, "import_document", json!({"text":text})).await;
    let inputs = sequences
        .iter()
        .take(3)
        .map(|(name, _)| json!({"id":name,"document":named(&imported,name)}))
        .collect::<Vec<_>>();
    let fragments=sequences.iter().take(3).map(|(name,_)|json!({"source":{"role":name,"input_id":name,"preparation":{"Digest":{"enzymes":["EcoRI"],"fragment_index":0,"repair":null}},"reverse":false},"treatment":null})).collect::<Vec<_>>();
    let result=call(&client,"propose_ligation_assembly",json!({"inputs":inputs,"target":named(&imported,"target"),"recipe":{"name":"Three compatible EcoRI fragments","fragments":fragments,"circular":true,"allow_nicked":false}})).await;
    assert_eq!(
        result["verification"]["status"],
        "ConsistentUnderAssumptions"
    );
    assert_eq!(result["junction_count"], 3);
    assert_eq!(result["product_physical"]["top_length"], 42);
    assert_eq!(
        result["product_physical"]["covalently_closed_in_model"],
        true
    );
    let full = evidence(&client, &result["artifact_id"]).await;
    assert_eq!(full["junctions"].as_array().unwrap().len(), 3);
    assert_eq!(full["junctions"][2]["closing"], true);
    for junction in full["junctions"].as_array().unwrap() {
        assert_eq!(junction["upstream_end"]["polarity"], "FivePrime");
        assert_eq!(junction["downstream_end"]["polarity"], "FivePrime");
        assert_eq!(junction["upstream_end"]["overhang"], "AATT");
        assert_eq!(junction["downstream_end"]["overhang"], "AATT");
    }
    assert_eq!(call(&client, "list_documents", json!({})).await["total"], 4);
    drop(client);
    server.abort();
}
