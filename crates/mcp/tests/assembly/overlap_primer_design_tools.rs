//! Literal original templates and independent intended product; planner is advisory.
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
async fn source_derived_ordered_oligos_are_advisory_until_actual_pcr_overlap_and_target_verification()
 {
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
    let imported=call(&client,"import_document",json!({"text":">left_template\nATGACCCCGGGGTTTCCAGA\n>right_template\nATCCGTTAGCGACCTTGGAA\n>independent_target\nATGACCCCGGGGTTTCCAGAATCCGTTAGCGACCTTGGAA\n>wrong_target\nATGACCCCGGGGTTTCCAGAATCCGTTAGCGACCTTGGAT\n"})).await;
    let inputs = json!([{"id":"l","document":named(&imported,"left_template")},{"id":"r","document":named(&imported,"right_template")}]);
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
    let request = json!({"name":"Two source-derived homology primers","fragments":[
        {"role":"left","input_id":"l","start":0,"end":20,"reverse":false},
        {"role":"right","input_id":"r","start":0,"end":20,"reverse":false}],"overlap":15,"circular":false,
        "primer_options":options,"profile_options":{"method":"HiFi","scan_repeated_overlaps":false,"max_comparisons":100000,"max_occurrences":100}});
    let args =
        json!({"inputs":inputs,"target":named(&imported,"independent_target"),"request":request});
    let designed = call(&client, "design_overlap_primers", args.clone()).await;
    assert_eq!(designed["state"], "advisory");
    assert_eq!(designed["verified"], false);
    assert!(designed.get("plan_revision").is_none());
    assert_eq!(
        designed["source_revisions"],
        json!([
            named(&imported, "left_template")["revision"],
            named(&imported, "right_template")["revision"]
        ])
    );
    let ordered = &designed["ordered_oligos"];
    assert_eq!(ordered.as_array().unwrap().len(), 2);
    assert_eq!(ordered[0]["forward"]["ordered_sequence"], "ATGACCCC");
    // Actual next-source prefix ATCCGTTAGCGACCT gives a 5' reverse-oligo
    // homology tail AGGTCGCTAACGGAT; no base comes from the target.
    assert_eq!(ordered[0]["reverse"]["primer"]["tail"], "AGGTCGCTAACGGAT");
    assert_eq!(
        ordered[0]["reverse"]["ordered_sequence"],
        "AGGTCGCTAACGGATTCTGGAAA"
    );
    assert_eq!(ordered[1]["forward"]["ordered_sequence"], "ATCCGTTA");
    assert_eq!(ordered[1]["reverse"]["ordered_sequence"], "TTCCAAGG");
    let saved = evidence(&client, &designed["artifact_id"]).await;
    assert_eq!(saved["target"], named(&imported, "independent_target"));
    assert_eq!(saved["report"]["recipe"], designed["recipe"]);
    let proposal=call(&client,"propose_overlap_assembly",json!({"inputs":inputs,"target":named(&imported,"independent_target"),"recipe":designed["recipe"]})).await;
    assert_eq!(proposal["verified"], false);
    assert_eq!(
        proposal["verification"]["status"],
        "ConsistentUnderAssumptions"
    );
    let mut wrong = args.clone();
    wrong["target"] = named(&imported, "wrong_target");
    let advisory_wrong = call(&client, "design_overlap_primers", wrong).await;
    assert_eq!(advisory_wrong["verified"], false);
    assert_eq!(advisory_wrong["recipe"], designed["recipe"]);
    assert_eq!(advisory_wrong["ordered_oligos"], designed["ordered_oligos"]);
    assert_ne!(
        advisory_wrong["context_revision"],
        designed["context_revision"]
    );
    let failed=call(&client,"propose_overlap_assembly",json!({"inputs":inputs,"target":named(&imported,"wrong_target"),"recipe":advisory_wrong["recipe"]})).await;
    assert_eq!(failed["verification"]["status"], "Invalid");
    let mut stale = args.clone();
    stale["inputs"][0]["document"]["revision"] = json!("stale");
    reject(&client, "design_overlap_primers", stale).await;
    let mut bad_interval = args.clone();
    bad_interval["request"]["fragments"][0]["end"] = json!(21);
    reject(&client, "design_overlap_primers", bad_interval).await;
    let mut existing_tail = args.clone();
    existing_tail["request"]["primer_options"]["forward_tail"] = json!("GG");
    reject(&client, "design_overlap_primers", existing_tail).await;
    let mut excessive = args.clone();
    excessive["request"]["overlap"] = json!(81);
    reject(&client, "design_overlap_primers", excessive).await;
    let mut duplicate = args;
    duplicate["request"]["fragments"][1]["role"] = json!("left");
    reject(&client, "design_overlap_primers", duplicate).await;
    assert_eq!(call(&client, "list_documents", json!({})).await["total"], 4);
    drop(client);
    server.abort();
}
