//! Independently specified junction DNA; no returned product becomes the target.
use collagene_core::analysis::Primer;
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
fn options() -> Value {
    json!({"flank_bases":5,"genetic_code":1,"frames":[{"reverse":false,"frame":0},{"reverse":false,"frame":1},{"reverse":false,"frame":2},{"reverse":true,"frame":0},{"reverse":true,"frame":1},{"reverse":true,"frame":2}],"max_junctions":128,"include_cds":true})
}
#[tokio::test]
async fn actual_ligation_crossing_codon_keeps_source_cds_conditional_and_target_failure_visible() {
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
    let left = "ATGACCCCGGGGTTTCCAGA";
    let right = "ATCCGTTAGCGACCTTGGAA";
    let text = format!(
        "LOCUS       left_template               20 bp    DNA     linear   SYN 01-JAN-2026\nFEATURES             Location/Qualifiers\n     CDS             1..20\n                     /label=\"source coding assertion\"\n                     /codon_start=1\nORIGIN\n        1 {}\n//\nLOCUS       right_template              20 bp    DNA     linear   SYN 01-JAN-2026\nFEATURES             Location/Qualifiers\nORIGIN\n        1 {}\n//\nLOCUS       target                      40 bp    DNA     linear   SYN 01-JAN-2026\nFEATURES             Location/Qualifiers\nORIGIN\n        1 atgaccccggggtttccagaatccgttagcgaccttggaa\n//\nLOCUS       wrong_target                40 bp    DNA     linear   SYN 01-JAN-2026\nFEATURES             Location/Qualifiers\nORIGIN\n        1 atgaccccggggtttccagaatccgttagcgaccttggat\n//\n",
        left.to_ascii_lowercase(),
        right.to_ascii_lowercase()
    );
    let imported = call(&client, "import_document", json!({"text":text})).await;
    let recipe = json!({"name":"Source-bound junction literal","fragments":[
        {"source":{"generation":null,"role":"left","input_id":"l","preparation":{"ExactPcr":{"forward":Primer::new("LF","ATGACCCC","").unwrap(),"reverse":Primer::new("LR","TCTGGAAA","").unwrap()}},"reverse":false},"treatment":"Phosphorylate"},
        {"source":{"generation":null,"role":"right","input_id":"r","preparation":{"ExactPcr":{"forward":Primer::new("RF","ATCCGTTA","").unwrap(),"reverse":Primer::new("RR","TTCCAAGG","").unwrap()}},"reverse":false},"treatment":"Phosphorylate"}
    ],"circular":false,"allow_nicked":false});
    let args = json!({"inputs":[{"id":"l","document":named(&imported,"left_template")},{"id":"r","document":named(&imported,"right_template")}],"target":named(&imported,"target"),"recipe":recipe});
    let plan = call(&client, "propose_ligation_assembly", args.clone()).await;
    assert_eq!(plan["verification"]["status"], "ConsistentUnderAssumptions");
    let query = json!({"plan_revision":plan["plan_revision"],"options":options()});
    let inspected = call(&client, "inspect_assembly_junctions", query.clone()).await;
    assert_eq!(
        inspected["verification_status"],
        "ConsistentUnderAssumptions"
    );
    assert_eq!(inspected["junction_count"], 1);
    assert_eq!(inspected["product_length"], 40);
    let full = evidence(&client, &inspected["artifact_id"]).await;
    let junction = &full["junctions"][0];
    assert_eq!(junction["boundary"], 20);
    assert_eq!(junction["closing"], false);
    // Five actual bases each side, not target-derived padding. Frame0 is
    // relative to this window: CCA GAA TCC G; GAA crosses the boundary at5.
    assert_eq!(junction["window"]["dna"], "CCAGAATCCG");
    assert_eq!(junction["window"]["boundary_offset"], 5);
    assert_eq!(junction["frames"][0]["crossing_codon"]["dna"], "GAA");
    assert_eq!(junction["frames"][0]["crossing_codon"]["window_start"], 3);
    assert_eq!(junction["physical"]["upstream_end"]["polarity"], "Blunt");
    assert!(
        full["verification"]["products"]
            .as_array()
            .unwrap()
            .last()
            .unwrap()["features"]
            .as_array()
            .unwrap()
            .iter()
            .any(|f| f["kind"] == "CDS")
    );
    assert_eq!(
        call(&client, "inspect_assembly_junctions", query).await["context_revision"],
        inspected["context_revision"]
    );
    let mut wrong = args.clone();
    wrong["target"] = named(&imported, "wrong_target");
    let wrong_plan = call(&client, "propose_ligation_assembly", wrong).await;
    let wrong_report = call(
        &client,
        "inspect_assembly_junctions",
        json!({"plan_revision":wrong_plan["plan_revision"],"options":options()}),
    )
    .await;
    assert_eq!(wrong_report["verification_status"], "Invalid");
    assert_eq!(wrong_report["junctions"][0]["window"]["dna"], "CCAGAATCCG");
    let mut reversed = args;
    reversed["recipe"]["fragments"][1]["source"]["reverse"] = json!(true);
    let reverse_plan = call(&client, "propose_ligation_assembly", reversed).await;
    let reverse_report = call(
        &client,
        "inspect_assembly_junctions",
        json!({"plan_revision":reverse_plan["plan_revision"],"options":options()}),
    )
    .await;
    assert_eq!(reverse_report["verification_status"], "Invalid");
    assert_eq!(
        reverse_report["junctions"][0]["window"]["dna"],
        "CCAGATTCCA"
    );
    reject(
        &client,
        "inspect_assembly_junctions",
        json!({"plan_revision":"unknown","options":options()}),
    )
    .await;
    let mut excessive = options();
    excessive["flank_bases"] = json!(151);
    reject(
        &client,
        "inspect_assembly_junctions",
        json!({"plan_revision":plan["plan_revision"],"options":excessive}),
    )
    .await;
    assert_eq!(call(&client, "list_documents", json!({})).await["total"], 4);
    drop(client);
    server.abort();
}
#[tokio::test]
async fn overlap_homology_and_circular_closing_window_are_actual_not_synthetic_fusions() {
    let (client, server) = setup().await;
    let imported = call(
        &client,
        "import_document",
        json!({"text":">left\nAAAACCCC\n>right\nCCCCGGGG\n>target\nAAAACCCCGGGG\n"}),
    )
    .await;
    let plan=call(&client,"propose_overlap_assembly",json!({"inputs":[{"id":"l","document":named(&imported,"left")},{"id":"r","document":named(&imported,"right")}],"target":named(&imported,"target"),"recipe":{"name":"Literal homology","fragments":[{"generation":null,"role":"l","input_id":"l","preparation":"WholeLinear","reverse":false},{"generation":null,"role":"r","input_id":"r","preparation":"WholeLinear","reverse":false}],"overlaps":[4],"circular":false,"options":{"method":"HiFi","scan_repeated_overlaps":false,"max_comparisons":10000,"max_occurrences":100}}})).await;
    let inspected = call(
        &client,
        "inspect_assembly_junctions",
        json!({"plan_revision":plan["plan_revision"],"options":options()}),
    )
    .await;
    let full = evidence(&client, &inspected["artifact_id"]).await;
    assert_eq!(
        full["junctions"][0]["homology"],
        json!([{"start":4,"end":8}])
    );
    let circles = "LOCUS       vector                      14 bp    DNA     circular SYN 01-JAN-2026\nFEATURES             Location/Qualifiers\nORIGIN\n        1 aaaagaattctttt\n//\nLOCUS       independent_circle          14 bp    DNA     circular SYN 01-JAN-2026\nFEATURES             Location/Qualifiers\nORIGIN\n        1 aaaagaattctttt\n//\n";
    let imported = call(&client, "import_document", json!({"text":circles})).await;
    let circle=call(&client,"propose_ligation_assembly",json!({"inputs":[{"id":"v","document":named(&imported,"vector")}],"target":named(&imported,"independent_circle"),"recipe":{"name":"Actual EcoRI reclosure","fragments":[{"source":{"generation":null,"role":"cut vector","input_id":"v","preparation":{"Digest":{"enzymes":["EcoRI"],"fragment_index":0,"repair":null}},"reverse":false},"treatment":null}],"circular":true,"allow_nicked":false}})).await;
    let inspected = call(
        &client,
        "inspect_assembly_junctions",
        json!({"plan_revision":circle["plan_revision"],"options":options()}),
    )
    .await;
    assert_eq!(
        inspected["verification_status"],
        "ConsistentUnderAssumptions"
    );
    assert_eq!(inspected["product_topology"], "Circular");
    let junction = &inspected["junctions"][0];
    assert_eq!(junction["boundary"], 0);
    assert_eq!(junction["closing"], true);
    assert_eq!(junction["window"]["dna"], "AAAAGAATTC");
    assert_eq!(junction["window"]["parts"].as_array().unwrap().len(), 2);
    assert_eq!(junction["physical"]["upstream_end"]["overhang"], "AATT");
    drop(client);
    server.abort();
}
