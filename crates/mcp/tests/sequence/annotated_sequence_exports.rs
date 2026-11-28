//! Real SDK annotated-text writers: explicit alphabets, loss review and immutable sources.
use base64::{Engine as _, engine::general_purpose::STANDARD};
use collagene_mcp::server::CollageneServer;
use rmcp::{RoleClient, ServiceExt, model::CallToolRequestParams, service::RunningService};
use serde_json::{Value, json};
pub(super) type Client = RunningService<RoleClient, ()>;
pub(super) async fn call(c: &Client, name: &str, args: Value) -> Value {
    let r = c
        .call_tool(
            CallToolRequestParams::new(name.to_string())
                .with_arguments(args.as_object().unwrap().clone()),
        )
        .await
        .unwrap();
    assert_ne!(r.is_error, Some(true), "{name}: {:?}", r.content);
    r.structured_content.unwrap()
}
pub(super) async fn reject(c: &Client, name: &str, args: Value) {
    let r = c
        .call_tool(
            CallToolRequestParams::new(name.to_string())
                .with_arguments(args.as_object().unwrap().clone()),
        )
        .await;
    assert!(r.is_err() || r.unwrap().is_error == Some(true));
}
pub(super) async fn setup() -> (Client, tokio::task::JoinHandle<()>) {
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
pub(super) fn document(result: &Value) -> Value {
    let d = &result["documents"][0];
    json!({"id":d["id"],"revision":d["revision"]})
}
pub(super) async fn bytes(c: &Client, result: &Value) -> Vec<u8> {
    let chunk = call(
        c,
        "read_artifact",
        json!({"artifact_id":result["artifact_id"],"start":0,"max_bytes":16384}),
    )
    .await;
    assert_eq!(chunk["done"], true);
    STANDARD.decode(chunk["text"].as_str().unwrap()).unwrap()
}
pub(super) async fn text(c: &Client, result: &Value) -> String {
    let mut out = String::new();
    let mut start = 0;
    loop {
        let page = call(
            c,
            "read_artifact",
            json!({"artifact_id":result["artifact_id"],"start":start,"max_bytes":16384}),
        )
        .await;
        out.push_str(page["text"].as_str().unwrap());
        if page["done"] == true {
            break;
        }
        start = page["next_start"].as_u64().unwrap();
    }
    out
}
#[tokio::test]
async fn nucleotide_annotations_and_rna_use_real_annotated_writers_without_changing_source() {
    let (c, server) = setup().await;
    let tools = c.list_all_tools().await.unwrap();
    assert_eq!(tools.len(), collagene_mcp::CURRENT_TOOL_COUNT);
    assert_eq!(
        tools
            .iter()
            .find(|t| t.name == "export_annotated_sequence")
            .unwrap()
            .input_schema["additionalProperties"],
        false
    );
    let gb = "LOCUS       literal                   16 bp    DNA     circular UNK 26-NOV-2026\nDEFINITION  Literal independent annotation.\nACCESSION   ZZ000001\nFEATURES             Location/Qualifiers\n     misc_feature    join(13..16,1..4)\n                     /label=\"origin spanning\"\n                     /note=\"literal quoted evidence\"\n     primer_bind     complement(5..12)\n                     /label=\"Reverse primer\"\nORIGIN\n        1 acgtacgtac gtacgt\n//\n";
    let imported = call(&c, "import_document", json!({"text":gb})).await;
    let d = document(&imported);
    let before = bytes(&c, &call(&c, "export_document", d.clone()).await).await;
    for format in ["GenBank", "Embl"] {
        let reply = call(
            &c,
            "export_annotated_sequence",
            json!({"document":d,"format":format,"accept_omissions":true}),
        )
        .await;
        assert_eq!(reply["state"], "prepared");
        let output = text(&c, &reply).await;
        assert!(output.contains("join(13..16,1..4)"));
        assert!(output.contains("complement(5..12)"));
        let parsed = collagene_core::formats::sequence::import_text(&output).unwrap();
        let m = &parsed.molecules[0];
        assert_eq!(m.sequence.as_str(), "ACGTACGTACGTACGT");
        assert_eq!(m.topology, collagene_model::Topology::Circular);
        assert_eq!(m.features.len(), 2);
        assert!(m.features.iter().any(|f| f.kind == "primer_bind"));
    }
    assert_eq!(
        before,
        bytes(&c, &call(&c, "export_document", d.clone()).await).await
    );
    let rna = call(
        &c,
        "import_compat_sequence",
        json!({"text":">rna\nAUGGCUUAANRY\n","profile":{"Fasta":{"alphabet":"Rna"}}}),
    )
    .await;
    let r = document(&rna);
    let review = call(
        &c,
        "export_annotated_sequence",
        json!({"document":r,"format":"GenBank"}),
    )
    .await;
    assert_eq!(review["state"], "review_required");
    assert!(review["artifact_id"].is_null());
    assert!(
        review["warnings"]
            .as_array()
            .unwrap()
            .iter()
            .any(|w| w.as_str().unwrap().contains("INSDC T encoding"))
    );
    let reply = call(
        &c,
        "export_annotated_sequence",
        json!({"document":r,"format":"GenBank","accept_omissions":true}),
    )
    .await;
    let output = text(&c, &reply).await;
    let m = collagene_core::formats::sequence::import_text(&output)
        .unwrap()
        .molecules
        .remove(0);
    assert_eq!(m.sequence.alphabet(), collagene_model::Alphabet::Rna);
    assert_eq!(m.sequence.as_str(), "AUGGCUUAANRY");
    reject(&c,"export_annotated_sequence",json!({"document":d,"format":{"GenPept":{"locus_date":"26-NOV-2026"}},"accept_omissions":true})).await;
    drop(c);
    server.abort();
}
#[tokio::test]
async fn protein_profiles_require_real_identifiers_dates_and_explicit_loss_acknowledgement() {
    let (c, server) = setup().await;
    let imported = call(
        &c,
        "import_compat_sequence",
        json!({"text":">literal\nACDEFGHIK\n","profile":{"Fasta":{"alphabet":"Protein"}}}),
    )
    .await;
    let d = document(&imported);
    let before = bytes(&c, &call(&c, "export_document", d.clone()).await).await;
    reject(
        &c,
        "export_annotated_sequence",
        json!({"document":d,"format":"GenBank"}),
    )
    .await;
    reject(
        &c,
        "export_annotated_sequence",
        json!({"document":d,"format":"Embl"}),
    )
    .await;
    reject(
        &c,
        "export_annotated_sequence",
        json!({"document":d,"format":{"GenPept":{}}}),
    )
    .await;
    reject(&c,"export_annotated_sequence",json!({"document":d,"format":{"GenPept":{"locus_date":"31-FEB-2026","entry_name":"LITERAL"}},"accept_omissions":true})).await;
    reject(
        &c,
        "export_annotated_sequence",
        json!({"document":d,"format":{"UniProt":{"entry_name":"LITERAL"}},"accept_omissions":true}),
    )
    .await;
    for (format, profile) in [
        (
            json!({"GenPept":{"locus_date":"26-NOV-2026","entry_name":"LITERAL"}}),
            collagene_core::formats::protein::Profile::GenPept,
        ),
        (
            json!({"UniProt":{"entry_name":"LITERAL","accession":"P69905"}}),
            collagene_core::formats::protein::Profile::UniProt,
        ),
    ] {
        let review = call(
            &c,
            "export_annotated_sequence",
            json!({"document":d,"format":format}),
        )
        .await;
        assert_eq!(review["state"], "review_required");
        assert!(review["artifact_id"].is_null());
        let reply = call(
            &c,
            "export_annotated_sequence",
            json!({"document":d,"format":format,"accept_omissions":true}),
        )
        .await;
        let output = text(&c, &reply).await;
        assert!(!output.contains("Reviewed;"));
        let m = collagene_core::formats::protein::import_bytes(output.as_bytes(), profile)
            .unwrap()
            .molecules
            .remove(0);
        assert_eq!(m.sequence.alphabet(), collagene_model::Alphabet::Protein);
        assert_eq!(m.sequence.as_str(), "ACDEFGHIK");
        if profile == collagene_core::formats::protein::Profile::UniProt {
            assert!(output.contains("Unreviewed;"));
            assert!(output.contains("P69905;"));
        } else {
            assert!(output.contains("26-NOV-2026"));
        }
    }
    assert_eq!(
        before,
        bytes(&c, &call(&c, "export_document", d.clone()).await).await
    );
    reject(&c,"export_annotated_sequence",json!({"document":d,"format":{"GenPept":{"locus_date":"26-NOV-2026","accession":"FAKE"}},"accept_omissions":true})).await;
    reject(&c,"export_annotated_sequence",json!({"document":d,"format":{"UniProt":{"entry_name":"LITERAL","accession":"P69905","reviewed":true}},"accept_omissions":true})).await;
    let mut wrong = d;
    wrong["revision"] = json!("0".repeat(64));
    reject(
        &c,
        "export_annotated_sequence",
        json!({"document":wrong,"format":{"GenPept":{"locus_date":"26-NOV-2026"}}}),
    )
    .await;
    drop(c);
    server.abort();
}
