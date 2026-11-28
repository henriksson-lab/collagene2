//! Saved oligo identity is independently checked before its use in actual PCR.
use collagene_core::sequence::analysis::Primer;
use collagene_mcp::server::CollageneServer;
use rmcp::{RoleClient, ServiceExt, model::CallToolRequestParams, service::RunningService};
use serde_json::{Value, json};
type Client = RunningService<RoleClient, ()>;
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
async fn reject(client: &Client, args: Value) {
    let result = client
        .call_tool(
            CallToolRequestParams::new("resolve_primer")
                .with_arguments(args.as_object().unwrap().clone()),
        )
        .await;
    assert!(result.is_err() || result.unwrap().is_error == Some(true));
}
const SAVED: &str = "LOCUS       saved_primers             16 bp    DNA     linear   SYN 01-JAN-2026\nFEATURES             Location/Qualifiers\n     primer_bind     1..4\n                     /label=\"saved forward\"\n                     /annealing_sequence=\"ACGA\"\n                     /5prime_tail=\"GG\"\n                     /oligo_sequence=\"GGACGA\"\n     primer_bind     complement(13..16)\n                     /label=\"saved reverse\"\n     primer_bind     1..4\n                     /label=\"intentional mismatch\"\n                     /annealing_sequence=\"TTTT\"\n                     /5prime_tail=\"A\"\n                     /oligo_sequence=\"ATTTT\"\n     primer_bind     1..4\n                     /annealing_sequence=\"ACGA\"\n                     /5prime_tail=\"GG\"\n                     /oligo_sequence=\"CCACGA\"\n     primer_bind     <1..4\n                     /label=\"uncertain inferred oligo\"\n     misc_feature    1..4\nORIGIN\n        1 acgaccccggggtaac\n//\n";
#[tokio::test]
async fn exact_saved_domains_and_reverse_location_reuse_keep_sequence_identity_not_binding_authority()
 {
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
    let client = ().serve(b).await.unwrap();
    assert_eq!(
        client
            .list_tools(Default::default())
            .await
            .unwrap()
            .tools
            .len(),
        collagene_mcp::CURRENT_TOOL_COUNT
    );
    let imported = call(&client, "import_document", json!({"text":SAVED})).await;
    let row = &imported["documents"][0];
    let document = json!({"id":row["id"],"revision":row["revision"]});
    let forward = call(
        &client,
        "resolve_primer",
        json!({"document":document,"feature_index":0}),
    )
    .await;
    assert_eq!(forward["resolved"]["source_revision"], document["revision"]);
    assert_eq!(forward["resolved"]["feature_index"], 0);
    assert_eq!(
        forward["resolved"]["primer"],
        serde_json::to_value(Primer::new("saved forward", "ACGA", "GG").unwrap()).unwrap()
    );
    assert_eq!(forward["resolved"]["sequence_origin"], "StoredDomains");
    let reverse = call(
        &client,
        "resolve_primer",
        json!({"document":document,"feature_index":1}),
    )
    .await;
    assert_eq!(
        reverse["resolved"]["primer"],
        serde_json::to_value(Primer::new("saved reverse", "GTTA", "").unwrap()).unwrap()
    );
    assert_eq!(reverse["resolved"]["sequence_origin"], "CurrentLocation");
    let mismatch = call(
        &client,
        "resolve_primer",
        json!({"document":document,"feature_index":2}),
    )
    .await;
    assert_eq!(
        mismatch["resolved"]["primer"],
        serde_json::to_value(Primer::new("intentional mismatch", "TTTT", "A").unwrap()).unwrap()
    );
    assert!(
        !mismatch["resolved"]["warnings"]
            .as_array()
            .unwrap()
            .is_empty()
    );
    for index in [3, 4, 5, 6] {
        reject(&client, json!({"document":document,"feature_index":index})).await;
    }
    let mut stale = document.clone();
    stale["revision"] = json!("stale");
    reject(&client, json!({"document":stale,"feature_index":0})).await;
    // Expected product is independently literal: GG + actual full template,
    // not a sequence taken from resolver output or an old binding annotation.
    let pcr=call(&client,"simulate_pcr",json!({"document":document,"forward":forward["resolved"]["primer"],"reverse":reverse["resolved"]["primer"]})).await;
    assert_eq!(pcr["product_count"], 1);
    assert_eq!(pcr["products"][0]["length"], 18);
    assert_eq!(
        pcr["products"][0]["sequence_hash"],
        blake3::hash(b"GGACGACCCCGGGGTAAC").to_hex().as_str()
    );
    assert_eq!(call(&client, "list_documents", json!({})).await["total"], 1);
    drop(client);
    server.abort();
}
