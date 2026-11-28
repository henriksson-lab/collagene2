//! Unified search and browser tracks use the real MCP SDK and immutable refs.
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
use base64::{Engine, engine::general_purpose::STANDARD};
use collagene_core::alignment_formats::{self, AlignmentFormat};
use collagene_model::Alphabet;
async fn import_alignment(c: &Client, data: &[u8]) -> Value {
    call(c,"import_alignment",json!({"data_base64":STANDARD.encode(data),"format":"AlignedFasta","alphabet":"Dna","name":"literal"})).await
}
async fn artifact_text(c: &Client, id: &Value) -> String {
    let mut start = 0;
    let mut text = String::new();
    loop {
        let r = call(
            c,
            "read_artifact",
            json!({"artifact_id":id,"start":start,"max_bytes":7}),
        )
        .await;
        text.push_str(r["text"].as_str().unwrap());
        if r["done"] == true {
            return text;
        }
        start = r["next_start"].as_u64().unwrap();
    }
}

#[tokio::test]
async fn nexus_and_phylip_exchange_preserve_gaps_and_require_explicit_name_loss() {
    let (client, server) = setup().await;
    let source=b"#NEXUS\nBEGIN DATA;\nDIMENSIONS NTAX=2 NCHAR=6;\nFORMAT DATATYPE=DNA GAP=-;\nMATRIX\n'alpha_long_identifier' AC-GT-\n'beta_long_identifier' A-TGT-\n;\nEND;\n";
    let imported=call(&client,"import_alignment",json!({"data_base64":STANDARD.encode(source),"format":"Nexus","alphabet":"Dna","name":"exchange"})).await;
    let reference = &imported["alignment"];
    let original = call(
        &client,
        "export_alignment",
        json!({"alignment":reference,"format":"Original"}),
    )
    .await;
    assert_eq!(
        artifact_text(&client, &original["artifact_id"])
            .await
            .as_bytes(),
        source
    );
    reject(
        &client,
        "export_alignment",
        json!({"alignment":reference,"format":"PhylipSequential"}),
    )
    .await;
    let exported=call(&client,"export_alignment",json!({"alignment":reference,"format":"PhylipSequential","options":{"truncate_phylip_identifiers":true,"block_width":4}})).await;
    assert_eq!(exported["renamings"]["alpha_long_identifier"], "alpha_long");
    assert_eq!(exported["renamings"]["beta_long_identifier"], "beta_long_");
    assert!(!exported["losses"].as_array().unwrap().is_empty());
    let text = artifact_text(&client, &exported["artifact_id"]).await;
    let reimported=call(&client,"import_alignment",json!({"data_base64":STANDARD.encode(text.as_bytes()),"format":"PhylipSequential","alphabet":"Dna","name":"reimported"})).await;
    let page=call(&client,"read_alignment",json!({"alignment":reimported["alignment"],"row_start":0,"max_rows":2,"column_start":0,"max_columns":6})).await;
    assert_eq!(page["rows"][0]["aligned"], "AC-GT-");
    assert_eq!(page["rows"][1]["aligned"], "A-TGT-");
    reject(&client,"export_alignment",json!({"alignment":reference,"format":"Native","options":{"truncate_phylip_identifiers":true,"block_width":4}})).await;
    let native = call(
        &client,
        "export_alignment",
        json!({"alignment":reference,"format":"Native"}),
    )
    .await;
    let native = artifact_text(&client, &native["artifact_id"]).await;
    let restored = call(
        &client,
        "import_alignment_document",
        json!({"data_base64":STANDARD.encode(native.as_bytes())}),
    )
    .await;
    assert_eq!(restored["column_count"], 6);
    client.cancel().await.unwrap();
    server.await.unwrap();
}
#[tokio::test]
async fn imported_alignment_window_preserves_gaps_original_bytes_and_native_identity() {
    let (c, server) = setup().await;
    let bytes = b">alpha first\r\nA-CG\r\n>beta\r\nAT-G\r\n";
    let imported = import_alignment(&c, bytes).await;
    let reference = &imported["alignment"];
    assert_eq!(imported["row_count"], 2);
    assert_eq!(imported["column_count"], 4);
    let page = call(
        &c,
        "read_alignment",
        json!({"alignment":reference,"row_start":0,"max_rows":1,"column_start":1,"max_columns":2}),
    )
    .await;
    assert_eq!(page["rows"][0]["aligned"], "-C");
    assert_eq!(page["rows"][0]["ungapped_positions"], json!([null, 1]));
    assert_eq!(page["next_row_start"], 1);
    assert_eq!(page["next_column_start"], 3);
    let original = call(
        &c,
        "export_alignment",
        json!({"alignment":reference,"format":"Original"}),
    )
    .await;
    assert_eq!(
        artifact_text(&c, &original["artifact_id"]).await.as_bytes(),
        bytes
    );
    let native = call(
        &c,
        "export_alignment",
        json!({"alignment":reference,"format":"Native"}),
    )
    .await;
    let native_text = artifact_text(&c, &native["artifact_id"]).await;
    let d = alignment_formats::load_alignment(native_text.as_bytes()).unwrap();
    assert_eq!(d.original_bytes().unwrap(), bytes);
    let reimported = call(
        &c,
        "import_alignment_document",
        json!({"data_base64":STANDARD.encode(native_text)}),
    )
    .await;
    assert_eq!(reimported["alignment"], *reference);
    let list = call(&c, "list_alignments", json!({})).await;
    assert_eq!(list["alignments"].as_array().unwrap().len(), 1);
    let mut stale = reference.clone();
    stale["revision"] = json!("stale");
    reject(
        &c,
        "read_alignment",
        json!({"alignment":stale,"row_start":0,"max_rows":1,"column_start":0,"max_columns":4}),
    )
    .await;
    reject(
        &c,
        "read_alignment",
        json!({"alignment":reference,"row_start":0,"max_rows":51,"column_start":0,"max_columns":4}),
    )
    .await;
    c.cancel().await.unwrap();
    server.await.unwrap();
}
#[tokio::test]
async fn explicit_row_extraction_is_atomic_and_typed_exports_disclose_loss() {
    let (c, server) = setup().await;
    let imported = import_alignment(&c, b">alpha\nA-CG\n>gap\n----\n").await;
    let reference = &imported["alignment"];
    let ids = imported["rows"]
        .as_array()
        .unwrap()
        .iter()
        .map(|r| r["id"].clone())
        .collect::<Vec<_>>();
    let before = call(&c, "list_documents", json!({})).await;
    reject(
        &c,
        "extract_alignment_rows",
        json!({"alignment":reference,"row_ids":ids}),
    )
    .await;
    assert_eq!(call(&c, "list_documents", json!({})).await, before);
    reject(
        &c,
        "extract_alignment_rows",
        json!({"alignment":reference,"row_ids":[ids[0],ids[0]]}),
    )
    .await;
    let rows = call(
        &c,
        "extract_alignment_rows",
        json!({"alignment":reference,"row_ids":[ids[0]]}),
    )
    .await;
    assert_eq!(rows["documents"][0]["length"], 3);
    assert_eq!(rows["documents"][0]["feature_count"], 0);
    let export = call(
        &c,
        "export_alignment",
        json!({"alignment":reference,"format":"Stockholm"}),
    )
    .await;
    assert!(!export["losses"].as_array().unwrap().is_empty());
    let data = artifact_text(&c, &export["artifact_id"]).await;
    let d = alignment_formats::import_alignment(
        data.as_bytes(),
        AlignmentFormat::Stockholm,
        Alphabet::Dna,
        "output",
    )
    .unwrap();
    assert_eq!(d.rows[0].aligned, "A-CG");
    reject(
        &c,
        "import_alignment_document",
        json!({"data_base64":STANDARD.encode(b"{}") }),
    )
    .await;
    reject(&c,"import_alignment",json!({"data_base64":STANDARD.encode(b">bad\nA.CG\n"),"format":"AlignedFasta","alphabet":"Dna","name":"bad"})).await;
    c.cancel().await.unwrap();
    server.await.unwrap();
}
#[tokio::test]
async fn row_previews_are_bounded_and_all_rows_remain_pageable() {
    let (c, server) = setup().await;
    let input = (0..103)
        .map(|i| format!(">row_{i}\nAC-G\n"))
        .collect::<String>();
    let imported = import_alignment(&c, input.as_bytes()).await;
    assert_eq!(imported["rows"].as_array().unwrap().len(), 100);
    assert_eq!(imported["rows_preview_truncated"], true);
    let final_page=call(&c,"read_alignment",json!({"alignment":imported["alignment"],"row_start":100,"max_rows":50,"column_start":0,"max_columns":1000})).await;
    assert_eq!(final_page["rows"].as_array().unwrap().len(), 3);
    assert_eq!(final_page["rows"][2]["identifier"], "row_102");
    assert_eq!(final_page["next_row_start"], Value::Null);
    c.cancel().await.unwrap();
    server.await.unwrap();
}
#[tokio::test]
async fn protein_properties_are_revision_bound_typed_estimates() {
    let (c, server) = setup().await;
    let imported = call(
        &c,
        "import_document",
        json!({"text":">protein [molecule=protein]\nACW"}),
    )
    .await;
    let d = &imported["documents"][0];
    let r = reference(d);
    let options = collagene_core::protein_properties::ProteinOptions::default();
    let result = call(
        &c,
        "analyze_protein",
        json!({"document":r,"options":options}),
    )
    .await;
    assert_eq!(result["source_revision"], d["revision"]);
    assert!(result["artifact_id"].is_string());
    let mut stale = r.clone();
    stale["revision"] = json!("stale");
    reject(
        &c,
        "analyze_protein",
        json!({"document":stale,"options":options}),
    )
    .await;
    let full = artifact(&c, &result["artifact_id"]).await;
    assert_eq!(full["source_revision"], d["revision"]);
    c.cancel().await.unwrap();
    server.await.unwrap();
}
#[tokio::test]
async fn native_workspace_import_exposes_retained_alignment_collections_and_deduplicates() {
    use collagene_core::workspace::{self, Workspace, WorkspaceAttachment};
    use collagene_model::{Molecule, Topology};
    let (c, server) = setup().await;
    let doc = alignment_formats::import_alignment(
        b">alpha\nAC-G\n",
        AlignmentFormat::AlignedFasta,
        Alphabet::Dna,
        "from workspace",
    )
    .unwrap();
    let bytes = alignment_formats::save_alignment(&doc).unwrap();
    let hash = blake3::hash(&bytes).to_hex().to_string();
    let mut frame = Workspace::from_documents(
        "alignment workspace",
        vec![Molecule::dna("reference", "ACGT", Topology::Linear).unwrap()],
        0,
    )
    .unwrap();
    frame.attachments.insert(
        hash.clone(),
        WorkspaceAttachment {
            kind: alignment_formats::WORKSPACE_ATTACHMENT_KIND.into(),
            name: "alignment".into(),
            media_type: "application/json".into(),
            bytes,
        },
    );
    let archive = workspace::save_workspace(&frame).unwrap();
    let imported = call(
        &c,
        "import_workspace_archive",
        json!({"data_base64":STANDARD.encode(&archive)}),
    )
    .await;
    assert_eq!(imported["deduplicated"], false);
    let list = call(&c, "list_alignments", json!({})).await;
    assert_eq!(list["alignments"].as_array().unwrap().len(), 1);
    assert_eq!(list["alignments"][0]["id"], hash);
    let imported = call(
        &c,
        "import_workspace_archive",
        json!({"data_base64":STANDARD.encode(archive)}),
    )
    .await;
    assert_eq!(imported["deduplicated"], true);
    assert_eq!(call(&c, "list_alignments", json!({})).await, list);
    c.cancel().await.unwrap();
    server.await.unwrap();
}
