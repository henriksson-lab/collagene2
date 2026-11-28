use base64::{Engine, engine::general_purpose::STANDARD};
use collagene_core::{
    cloning::pcr::primers::{BindingOptions, DesignOptions, ThermoConditions},
    sequence::analysis::Primer,
};
use collagene_mcp::server::CollageneServer;
use rmcp::{ServiceExt, model::CallToolRequestParams};
use serde_json::{Value, json};

#[tokio::test]
async fn primer_evidence_design_pdf_artifacts_and_provider_selection_contracts() {
    let (a, b) = tokio::io::duplex(65536);
    let server = tokio::spawn(async {
        CollageneServer::default()
            .serve(a)
            .await
            .unwrap()
            .waiting()
            .await
            .unwrap();
    });
    let client = ().serve(b).await.unwrap();
    macro_rules! call {
        ($name:expr,$args:expr) => {{
            let args: Value = $args;
            let result = client
                .call_tool(
                    CallToolRequestParams::new($name)
                        .with_arguments(args.as_object().unwrap().clone()),
                )
                .await
                .unwrap();
            assert_ne!(
                result.is_error,
                Some(true),
                "{}: {:?}",
                $name,
                result.content
            );
            result.structured_content.unwrap()
        }};
    }
    macro_rules! reject {
        ($name:expr,$args:expr) => {{
            let args: Value = $args;
            let result = client
                .call_tool(
                    CallToolRequestParams::new($name)
                        .with_arguments(args.as_object().unwrap().clone()),
                )
                .await;
            assert!(
                result.is_err() || result.unwrap().is_error == Some(true),
                "accepted invalid {}",
                $name
            );
        }};
    }
    let tools = client.list_all_tools().await.unwrap();
    assert_eq!(tools.len(), collagene_mcp::CURRENT_TOOL_COUNT);
    for name in [
        "analyze_primer",
        "design_primers",
        "list_source_sequences",
        "export_protocol_pdf",
    ] {
        assert!(tools.iter().any(|t| t.name == name));
    }
    let schema = &tools
        .iter()
        .find(|t| t.name == "import_url")
        .unwrap()
        .input_schema;
    assert!(schema["properties"].get("sequence_id").is_some());
    assert!(
        !schema["required"]
            .as_array()
            .unwrap()
            .contains(&json!("sequence_id"))
    );
    let imported = call!(
        "import_document",
        json!({"text":">template\nACGACCCCGGGGTAAC\n>left\nAAAACCCC\n>right\nCCCCGGGG\n>target\nAAAACCCCGGGG\n"})
    );
    let reference = |name: &str| {
        let d = imported["documents"]
            .as_array()
            .unwrap()
            .iter()
            .find(|d| d["name"] == name)
            .unwrap();
        json!({"id":d["id"],"revision":d["revision"]})
    };
    let source = reference("template");
    let primer = Primer::new("forward", "ACGA", "GG").unwrap();
    let query = json!({"document":source,"primer":primer,"conditions":ThermoConditions::default(),"binding_options":BindingOptions::default(),"partner":Primer::new("reverse","GTTA","CC").unwrap()});
    let report = call!("analyze_primer", query.clone());
    assert_eq!(report["analysis"]["ordered_sequence"], "GGACGA");
    assert_eq!(report["analysis"]["annealing_gc"]["minimum_fraction"], 0.5);
    assert_eq!(report["binding_count"], 1);
    assert_eq!(report["bindings"][0]["start"], 0);
    assert_eq!(report["bindings"][0]["three_prime_exact"], 4);
    assert_eq!(report["search_truncated"], false);
    let evidence = call!(
        "read_artifact",
        json!({"artifact_id":report["artifact_id"],"start":0,"max_bytes":16384})
    );
    let evidence: Value = serde_json::from_str(evidence["text"].as_str().unwrap()).unwrap();
    assert_eq!(evidence["document"], source);
    assert!(evidence["engine_version"].is_string());
    let mut stale = query;
    stale["document"]["revision"] = json!("stale");
    reject!("analyze_primer", stale);
    let options = DesignOptions {
        target_start: 0,
        target_end: 16,
        min_length: 8,
        max_length: 8,
        target_tm_c: 20.,
        max_tm_difference_c: 100.,
        min_gc_fraction: 0.,
        max_gc_fraction: 1.,
        max_homopolymer: 20,
        max_three_prime_complementarity: 20,
        forward_tail: "AA".into(),
        reverse_tail: "TT".into(),
        ..Default::default()
    };
    let design = call!(
        "design_primers",
        json!({"document":source,"options":options})
    );
    assert_eq!(design["count"], 1);
    assert_eq!(
        design["pairs"][0]["forward"]["primer"]["annealing"]["bases"],
        "ACGACCCC"
    );
    assert_eq!(
        design["pairs"][0]["reverse"]["primer"]["annealing"]["bases"],
        "GTTACCCC"
    );
    assert_eq!(design["pairs"][0]["product_length"], 20);
    reject!(
        "design_primers",
        json!({"document":{"id":source["id"],"revision":"stale"},"options":options})
    );
    // Invalid provider/selection combinations fail before any network access.
    reject!("list_source_sequences", json!({"url":"not-a-provider"}));
    reject!("import_url", json!({"url":"NC_001422.1","sequence_id":9}));
    let plan = call!(
        "propose_plan",
        json!({"name":"PCR α protocol","target":reference("target"),"inputs":[{"id":"left","document":reference("left")},{"id":"right","document":reference("right")}],"steps":[{"Overlap":{"fragments":[{"Input":"left"},{"Input":"right"}],"overlaps":[4],"circular":false}}]})
    );
    let exported = call!(
        "export_protocol_pdf",
        json!({"plan_revision":plan["plan_revision"]})
    );
    assert_eq!(exported["status"], "ConsistentUnderAssumptions");
    assert_eq!(exported["target_revision"], reference("target")["revision"]);
    assert_eq!(exported["mime_type"], "application/pdf");
    let mut encoded = String::new();
    let mut offset = 0u64;
    loop {
        let chunk = call!(
            "read_artifact",
            json!({"artifact_id":exported["artifact_id"],"start":offset,"max_bytes":16384})
        );
        assert_eq!(chunk["encoding"], "base64");
        assert_eq!(chunk["mime_type"], "application/pdf");
        encoded.push_str(chunk["text"].as_str().unwrap());
        offset = chunk["next_start"].as_u64().unwrap();
        if chunk["done"] == true {
            break;
        }
    }
    let pdf = STANDARD.decode(encoded).unwrap();
    assert!(pdf.starts_with(b"%PDF-"));
    assert!(pdf.ends_with(b"%%EOF") || pdf.ends_with(b"%%EOF\n"));
    assert_eq!(pdf.len() as u64, exported["bytes"].as_u64().unwrap());
    assert_eq!(
        blake3::hash(&pdf).to_hex().as_str(),
        exported["artifact_id"].as_str().unwrap()
    );
    assert_eq!(
        call!(
            "export_protocol_pdf",
            json!({"plan_revision":plan["plan_revision"]})
        )["artifact_id"],
        exported["artifact_id"]
    );
    reject!("export_protocol_pdf", json!({"plan_revision":"unknown"}));
    // Read-only scientific/export tools have not created product documents.
    assert_eq!(
        call!("list_documents", json!({}))["documents"]
            .as_array()
            .unwrap()
            .len(),
        4
    );
    client.cancel().await.unwrap();
    server.await.unwrap();
}
