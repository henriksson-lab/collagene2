//! Real SDK contracts for physical previews and independently target-bound plans.
use collagene_core::{
    Molecule, Topology,
    chemistry::{AnnealOptions, EndRepairMethod},
    formats,
};
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
async fn call(client: &Client, name: &str, args: Value) -> Value {
    let result = client
        .call_tool(
            CallToolRequestParams::new(name.to_owned())
                .with_arguments(args.as_object().unwrap().clone()),
        )
        .await
        .unwrap();
    assert_ne!(result.is_error, Some(true), "{name}: {:?}", result.content);
    result.structured_content.unwrap()
}
async fn reject(client: &Client, name: &str, args: Value) {
    let result = client
        .call_tool(
            CallToolRequestParams::new(name.to_owned())
                .with_arguments(args.as_object().unwrap().clone()),
        )
        .await;
    assert!(
        result.is_err() || result.unwrap().is_error == Some(true),
        "accepted invalid {name}"
    );
}
fn reference(value: &Value) -> Value {
    json!({"id":value["id"],"revision":value["revision"]})
}
async fn import(client: &Client, name: &str, sequence: &str, topology: Topology) -> Value {
    let m = Molecule::dna(name, sequence, topology).unwrap();
    let text = formats::export_genbank(&m).unwrap();
    reference(&call(client, "import_document", json!({"text":text})).await["documents"][0])
}
async fn artifact(client: &Client, id: &Value) -> Value {
    let mut start = 0;
    let mut text = String::new();
    loop {
        let result = call(
            client,
            "read_artifact",
            json!({"artifact_id":id,"start":start,"max_bytes":16384}),
        )
        .await;
        text.push_str(result["text"].as_str().unwrap());
        if result["done"] == true {
            break;
        }
        start = result["next_start"].as_u64().unwrap();
    }
    serde_json::from_str(&text).unwrap()
}
#[tokio::test]
async fn digest_repair_preview_preserves_sources_and_a_separate_target_is_reconstructed() {
    let (client, server) = setup().await;
    let source = import(&client, "source", "GAATTCAAAA", Topology::Circular).await;
    let target = import(&client, "target", "AATTCAAAAGAATT", Topology::Circular).await;
    let result = call(
        &client,
        "simulate_end_repair",
        json!({"document":source,"enzymes":["EcoRI"],"fragment_index":0,"method":"Blunt"}),
    )
    .await;
    assert_eq!(result["length"], 14);
    assert_eq!(result["top_added_count"], 4);
    assert_eq!(result["bottom_added_count"], 4);
    let report = artifact(&client, &result["artifact_id"]).await;
    assert_eq!(report["repair"]["source_revision"], source["revision"]);
    assert_eq!(
        report["repair"]["repair"]["duplex"]["top"]["bases"],
        "AATTCAAAAGAATT"
    );
    let plan=call(&client,"propose_plan",json!({"name":"repair and close","target":target,"inputs":[{"id":"s","document":source}],"steps":[{"Digest":{"template":{"Input":"s"},"enzymes":["EcoRI"],"fragment_index":0}},{"RepairEnds":{"template":{"Step":0},"method":"Blunt"}},{"Ligate":{"fragments":[{"source":{"Step":1},"reverse":false}],"circular":true,"allow_nicked":false}}]})).await;
    let checked = call(
        &client,
        "verify_plan",
        json!({"plan_revision":plan["plan_revision"]}),
    )
    .await;
    assert_eq!(checked["status"], "ConsistentUnderAssumptions");
    let mut stale = source.clone();
    stale["revision"] = json!("0".repeat(64));
    reject(
        &client,
        "simulate_end_repair",
        json!({"document":stale,"enzymes":["EcoRI"],"fragment_index":0,"method":"Blunt"}),
    )
    .await;
    reject(
        &client,
        "simulate_end_repair",
        json!({"document":source,"enzymes":[],"fragment_index":0,"method":"Blunt"}),
    )
    .await;
    let documents = call(&client, "list_documents", json!({})).await;
    assert_eq!(documents["total"], 2);
    client.cancel().await.unwrap();
    server.await.unwrap();
}
#[tokio::test]
async fn annealing_is_explicitly_unphosphorylated_until_the_plan_performs_terminal_treatment() {
    let (client, server) = setup().await;
    let top = import(&client, "top", "ACGTCG", Topology::Linear).await;
    let bottom = import(&client, "bottom", "CGACGT", Topology::Linear).await;
    let target = import(&client, "target", "ACGTCG", Topology::Circular).await;
    let options = AnnealOptions {
        min_paired_bases: 4,
        ..Default::default()
    };
    let value = serde_json::to_value(&options).unwrap();
    let result = call(
        &client,
        "simulate_oligo_annealing",
        json!({"top":top,"bottom":bottom,"options":value}),
    )
    .await;
    assert_eq!(result["candidate_count"], 1);
    assert_eq!(result["unique"], true);
    let report = artifact(&client, &result["artifact_id"]).await;
    assert_eq!(report["annealing"]["top_revision"], top["revision"]);
    assert_eq!(
        report["annealing"]["annealing"]["candidates"][0]["duplex"]["left"]["five_prime"],
        "Hydroxyl"
    );
    let steps = json!([{"AnnealOligos":{"top":{"Input":"top"},"bottom":{"Input":"bottom"},"options":value,"candidate_index":0}},{"TreatEnds":{"template":{"Step":0},"treatment":"Phosphorylate"}},{"Ligate":{"fragments":[{"source":{"Step":1},"reverse":false}],"circular":true,"allow_nicked":false}}]);
    let plan=call(&client,"propose_plan",json!({"name":"declared oligos","target":target,"inputs":[{"id":"top","document":top},{"id":"bottom","document":bottom}],"steps":steps})).await;
    let checked = call(
        &client,
        "verify_plan",
        json!({"plan_revision":plan["plan_revision"]}),
    )
    .await;
    assert_eq!(checked["status"], "ConsistentUnderAssumptions");
    let bad_steps = json!([{"AnnealOligos":{"top":{"Input":"top"},"bottom":{"Input":"bottom"},"options":value,"candidate_index":0}},{"Ligate":{"fragments":[{"source":{"Step":0},"reverse":false}],"circular":true,"allow_nicked":false}}]);
    let plan=call(&client,"propose_plan",json!({"name":"missing phosphorylation","target":target,"inputs":[{"id":"top","document":top},{"id":"bottom","document":bottom}],"steps":bad_steps})).await;
    assert_eq!(
        call(
            &client,
            "verify_plan",
            json!({"plan_revision":plan["plan_revision"]})
        )
        .await["status"],
        "Invalid"
    );
    let mut invalid = value;
    invalid["top_chemistry"]["unmodeled_protecting_group"] = json!(true);
    reject(
        &client,
        "simulate_oligo_annealing",
        json!({"top":top,"bottom":bottom,"options":invalid}),
    )
    .await;
    client.cancel().await.unwrap();
    server.await.unwrap();
}
#[tokio::test]
async fn anneal_limits_alternatives_and_known_physical_inputs_remain_explicit() {
    let (client, server) = setup().await;
    let top = import(&client, "repeat top", "AAAAAA", Topology::Linear).await;
    let bottom = import(&client, "repeat bottom", "TTTTTT", Topology::Linear).await;
    let result=call(&client,"simulate_oligo_annealing",json!({"top":top,"bottom":bottom,"options":AnnealOptions{min_paired_bases:2,max_candidates:1,..Default::default()}})).await;
    assert_eq!(result["complete"], false);
    assert_eq!(result["unique"], false);
    assert_eq!(result["candidate_count"], 1);
    let circ = import(
        &client,
        "closed strand reference",
        "ACGTCG",
        Topology::Circular,
    )
    .await;
    reject(&client,"simulate_oligo_annealing",json!({"top":circ,"bottom":bottom,"options":AnnealOptions{min_paired_bases:2,..Default::default()}})).await;
    reject(&client,"simulate_end_repair",json!({"document":top,"enzymes":["NotInCatalog"],"fragment_index":0,"method":EndRepairMethod::Blunt})).await;
    client.cancel().await.unwrap();
    server.await.unwrap();
}
