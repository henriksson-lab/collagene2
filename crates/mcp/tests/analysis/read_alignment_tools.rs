//! Real SDK dispatch with independently authored lexical SAM. No aligner,
//! samtools, reference fetch or imported positive reference-check authority.
use base64::{Engine as _, engine::general_purpose::STANDARD};
use collagene_mcp::server::CollageneServer;
use rmcp::{RoleClient, ServiceExt, model::CallToolRequestParams, service::RunningService};
use serde_json::{Value, json};
type Client = RunningService<RoleClient, ()>;
const SAM:&[u8]=b"@HD\tVN:1.6\tSO:unknown\r\n@SQ\tSN:chr1\tLN:16\tM5:00000000000000000000000000000000\tUR:https://invalid.test/never-fetch\r\nobserved\t0\tchr1\t1\t60\t4M\t*\t0\t0\tACGT\t*\tNM:i:0\tMD:Z:4\r\nmissing\t256\tchr1\t5\t255\t4M\t*\t0\t0\t*\t*\r\nunmapped\t2052\t*\t0\t0\t*\t*\t0\t0\t*\t*";
async fn call(c: &Client, name: &str, args: Value) -> Value {
    let response = c
        .call_tool(
            CallToolRequestParams::new(name.to_owned())
                .with_arguments(args.as_object().unwrap().clone()),
        )
        .await
        .unwrap();
    assert_ne!(
        response.is_error,
        Some(true),
        "{name}: {:?}",
        response.content
    );
    response.structured_content.unwrap()
}
async fn reject(c: &Client, name: &str, args: Value) {
    let response = c
        .call_tool(
            CallToolRequestParams::new(name.to_owned())
                .with_arguments(args.as_object().unwrap().clone()),
        )
        .await;
    assert!(
        response.is_err() || response.unwrap().is_error == Some(true),
        "{name} unexpectedly accepted invalid input"
    );
}
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
async fn artifact(c: &Client, id: &Value) -> (Vec<u8>, String) {
    let mut text = String::new();
    let mut start = 0;
    loop {
        let page = call(
            c,
            "read_artifact",
            json!({"artifact_id":id,"start":start,"max_bytes":1024}),
        )
        .await;
        text.push_str(page["text"].as_str().unwrap());
        if page["done"] == true {
            let encoding = page["encoding"].as_str().unwrap().to_owned();
            return (
                if encoding == "base64" {
                    STANDARD.decode(text).unwrap()
                } else {
                    text.into_bytes()
                },
                encoding,
            );
        }
        start = page["next_start"].as_u64().unwrap();
    }
}
async fn document(c: &Client, dna: &str) -> Value {
    let result = call(
        c,
        "import_document",
        json!({"text":format!(">explicit_actual_DNA\n{dna}\n")}),
    )
    .await;
    json!({"id":result["documents"][0]["id"],"revision":result["documents"][0]["revision"]})
}

#[tokio::test]
async fn sam_bam_sdk_preserves_original_flags_missing_calls_and_actual_reference_evidence() {
    let (client, server) = setup().await;
    let tools = client.list_tools(Default::default()).await.unwrap().tools;
    assert_eq!(tools.len(), collagene_mcp::CURRENT_TOOL_COUNT);
    for name in [
        "import_read_alignments",
        "inspect_read_alignments",
        "check_read_alignment",
        "export_read_alignments",
    ] {
        let tool = tools.iter().find(|t| t.name == name).unwrap();
        assert_eq!(tool.input_schema["additionalProperties"], false);
    }
    let reference = document(&client, "ACGTACGTACGTACGT").await;
    let other = document(&client, "TTTTTTTTTTTTTTTT").await;
    let before = call(&client, "list_documents", json!({})).await;
    let imported = call(
        &client,
        "import_read_alignments",
        json!({"data_base64":STANDARD.encode(SAM),"format":"sam"}),
    )
    .await;
    assert_eq!(imported["file_complete"], true);
    assert_eq!(imported["record_count"], 3);
    assert_eq!(
        imported["source_blake3"],
        blake3::hash(SAM).to_hex().to_string()
    );
    let collection = imported["collection"].clone();
    assert_eq!(collection["id"], collection["revision"]);
    let page = call(
        &client,
        "inspect_read_alignments",
        json!({"collection":collection,"limit":2}),
    )
    .await;
    assert_eq!(page["done"], false);
    assert_eq!(page["next_start"], 2);
    assert_eq!(page["records"][0]["flags"], 0);
    assert_eq!(page["records"][0]["reported_reference_start"], 0);
    assert_eq!(page["records"][0]["sequence_present"], true);
    assert_eq!(page["records"][0]["qualities_present"], false);
    assert_eq!(page["records"][1]["sequence_present"], false);
    assert_eq!(page["records"][1]["sequence_length"], Value::Null);
    assert_eq!(page["records"][1]["qualities_length"], Value::Null);
    assert_eq!(page["records"][1]["mapping_quality"], Value::Null);
    assert_eq!(page["records"][1]["secondary"], true);
    let last = call(
        &client,
        "inspect_read_alignments",
        json!({"collection":collection,"start":2}),
    )
    .await;
    assert_eq!(last["records"][0]["flags"], 2052);
    assert_eq!(last["records"][0]["unmapped"], true);
    assert_eq!(last["records"][0]["supplementary"], true);
    let args = json!({"collection":collection,"record_index":0,"reference":reference,"reference_name":"chr1"});
    let checked = call(&client, "check_read_alignment", args.clone()).await;
    assert_eq!(checked["reference_bound"], true);
    assert_eq!(checked["complete"], true);
    assert_eq!(checked["summary"]["certain_agreements"], 4);
    assert_eq!(checked["summary"]["certain_disagreements"], 0);
    assert_eq!(checked["qualities_present"], false);
    let evidence: Value =
        serde_json::from_slice(&artifact(&client, &checked["artifact_id"]).await.0).unwrap();
    assert_eq!(evidence["source_blake3"], collection["id"]);
    assert_eq!(evidence["record_index"], 0);
    assert_eq!(
        evidence["report"]["binding"]["reference_revision"],
        reference["revision"]
    );
    let mut changed = args.clone();
    changed["reference"] = other;
    let changed = call(&client, "check_read_alignment", changed).await;
    assert_eq!(
        changed["summary"]["certain_disagreements"], 3,
        "NM:i:0/MD:4/M5 cannot override actual DNA"
    );
    let mut missing = args.clone();
    missing["record_index"] = json!(1);
    let missing = call(&client, "check_read_alignment", missing).await;
    assert_eq!(missing["sequence_present"], false);
    assert_eq!(missing["complete"], false);
    assert_eq!(missing["column_count"], 0);
    assert!(
        missing["incomplete"]
            .as_array()
            .unwrap()
            .iter()
            .any(|v| v == "MissingSequence")
    );
    let original = call(
        &client,
        "export_read_alignments",
        json!({"collection":collection,"format":"original"}),
    )
    .await;
    assert_eq!(artifact(&client, &original["artifact_id"]).await.0, SAM);
    let lexical = call(
        &client,
        "export_read_alignments",
        json!({"collection":collection,"format":"sam"}),
    )
    .await;
    assert_eq!(artifact(&client, &lexical["artifact_id"]).await.0, SAM);
    let bam = call(
        &client,
        "export_read_alignments",
        json!({"collection":collection,"format":"bam"}),
    )
    .await;
    assert!(bam["notices"].is_array());
    let binary = artifact(&client, &bam["artifact_id"]).await.0;
    assert_eq!(&binary[..2], &[31, 139]);
    let imported_bam = call(
        &client,
        "import_read_alignments",
        json!({"data_base64":STANDARD.encode(&binary),"format":"bam"}),
    )
    .await;
    let bam_original = call(
        &client,
        "export_read_alignments",
        json!({"collection":imported_bam["collection"],"format":"original"}),
    )
    .await;
    assert_eq!(
        artifact(&client, &bam_original["artifact_id"]).await.0,
        binary
    );
    let bam_page = call(
        &client,
        "inspect_read_alignments",
        json!({"collection":imported_bam["collection"]}),
    )
    .await;
    assert_eq!(bam_page["records"][1]["sequence_present"], false);
    assert_eq!(bam_page["records"][1]["qualities_present"], false);
    assert!(bam_page["records"][0]["virtual_offset"]["compressed"].is_string());
    assert_eq!(call(&client, "list_documents", json!({})).await, before);
    client.cancel().await.unwrap();
    server.await.unwrap();
}

#[tokio::test]
async fn read_alignment_sdk_limits_and_stale_mapping_fail_without_partial_publication() {
    let (client, server) = setup().await;
    let reference = document(&client, "ACGTACGTACGTACGT").await;
    let short = document(&client, "ACGT").await;
    let imported = call(
        &client,
        "import_read_alignments",
        json!({"data_base64":STANDARD.encode(SAM),"format":"sam"}),
    )
    .await;
    let collection = imported["collection"].clone();
    let args = json!({"collection":collection,"record_index":0,"reference":reference,"reference_name":"chr1"});
    let mut stale = args.clone();
    stale["reference"]["revision"] = json!("0".repeat(64));
    reject(&client, "check_read_alignment", stale).await;
    let mut wrong = args.clone();
    wrong["reference_name"] = json!("wrong");
    reject(&client, "check_read_alignment", wrong).await;
    let mut wrong = args;
    wrong["reference"] = short;
    reject(&client, "check_read_alignment", wrong).await;
    reject(
        &client,
        "inspect_read_alignments",
        json!({"collection":collection,"limit":101}),
    )
    .await;
    reject(
        &client,
        "inspect_read_alignments",
        json!({"collection":{"id":collection["id"],"revision":"0".repeat(64)}}),
    )
    .await;
    reject(
        &client,
        "export_read_alignments",
        json!({"collection":collection,"format":"sam","max_bytes":1}),
    )
    .await;
    let mut unique = SAM.to_vec();
    unique.extend_from_slice(b"\r\nbroken\t");
    let hash = blake3::hash(&unique).to_hex().to_string();
    reject(
        &client,
        "import_read_alignments",
        json!({"data_base64":STANDARD.encode(&unique),"format":"sam"}),
    )
    .await;
    reject(
        &client,
        "inspect_read_alignments",
        json!({"collection":{"id":hash,"revision":hash}}),
    )
    .await;
    let mut changed = SAM.to_vec();
    let insertion = changed
        .windows(b"observed\t".len())
        .position(|bytes| bytes == b"observed\t")
        .unwrap();
    changed.splice(
        insertion..insertion,
        b"@CO\tbounded atomic fixture\r\n".iter().copied(),
    );
    let hash = blake3::hash(&changed).to_hex().to_string();
    reject(
        &client,
        "import_read_alignments",
        json!({"data_base64":STANDARD.encode(&changed),"format":"sam","max_records":1}),
    )
    .await;
    reject(
        &client,
        "inspect_read_alignments",
        json!({"collection":{"id":hash,"revision":hash}}),
    )
    .await;
    reject(
        &client,
        "import_read_alignments",
        json!({"data_base64":STANDARD.encode(SAM),"format":"cram"}),
    )
    .await;
    reject(
        &client,
        "import_read_alignments",
        json!({"data_base64":"Q1JBTQ==","format":"bam"}),
    )
    .await;
    assert_eq!(
        call(
            &client,
            "inspect_read_alignments",
            json!({"collection":collection})
        )
        .await["record_count"],
        3
    );
    client.cancel().await.unwrap();
    server.await.unwrap();
}
