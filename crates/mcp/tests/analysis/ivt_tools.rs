//! Actual in-memory SDK checks of immutable-source IVT declaration evidence.
use collagene_core::rna_chemistry::ivt;
use collagene_mcp::server::CollageneServer;
use rmcp::{RoleClient, ServiceExt, model::CallToolRequestParams, service::RunningService};
use serde_json::{Value, json};

type Client = RunningService<RoleClient, ()>;
const DNA: &str = "TAATACGACTCACTATAAGTTAAAA";
const FASTA: &str = ">IVT declared source\nTAATACGACTCACTATAAGTTAAAA\n";

async fn call(client: &Client, name: &str, args: Value) -> Value {
    let reply = client
        .call_tool(
            CallToolRequestParams::new(name.to_owned())
                .with_arguments(args.as_object().unwrap().clone()),
        )
        .await
        .unwrap();
    assert_ne!(reply.is_error, Some(true), "{name}: {:?}", reply.content);
    reply.structured_content.unwrap()
}
async fn reject(client: &Client, args: Value) {
    match client
        .call_tool(
            CallToolRequestParams::new("assess_ivt_template")
                .with_arguments(args.as_object().unwrap().clone()),
        )
        .await
    {
        Ok(reply) => assert_eq!(reply.is_error, Some(true)),
        Err(rmcp::service::ServiceError::McpError(error)) => assert_eq!(error.code.0, -32602),
        Err(error) => panic!("unexpected transport error: {error}"),
    }
}
async fn artifact(client: &Client, id: &Value) -> Value {
    let (mut start, mut text) = (0, String::new());
    loop {
        let page = call(
            client,
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
fn request(document: &Value, profile: &Value) -> Value {
    json!({
        "source_revision":document["revision"],
        "profile":profile["reference"],
        "transcript":{"start":17,"length":8,"coding_strand":"Forward"},
        "cap":"CleanCapAgN7113", "initiator":"NoAdditionalRequirement",
        "linearization":{"coding_cut_offset":8,"template_cut_offset":8,
            "provenance":{"statement":"Supplied blunt geometry, no cleavage asserted","reference":null}},
        "poly_a":{"start_offset":4,"end_offset":8},
        "provenance":{"statement":"Literal declared coding-strand interval","reference":null}
    })
}

#[tokio::test]
async fn ivt_discovery_assessment_and_artifact_preserve_source_and_typed_scope() {
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
    let tools = client.list_all_tools().await.unwrap();
    assert_eq!(tools.len(), collagene_mcp::CURRENT_TOOL_COUNT);
    for name in ["list_ivt_profiles", "assess_ivt_template"] {
        assert!(tools.iter().any(|tool| tool.name == name));
    }
    let schema = serde_json::to_string(
        &tools
            .iter()
            .find(|tool| tool.name == "assess_ivt_template")
            .unwrap()
            .input_schema,
    )
    .unwrap();
    assert!(
        schema.contains("CleanCapAgN7113")
            && schema.contains("coding_cut_offset")
            && schema.contains("poly_a")
    );
    let schema: Value = serde_json::from_str(&schema).unwrap();
    let declaration_schema = schema["$defs"]
        .as_object()
        .unwrap()
        .values()
        .find(|definition| {
            definition["properties"].get("transcript").is_some()
                && definition["properties"].get("cap").is_some()
        })
        .expect("typed IVT declaration schema");
    let required = declaration_schema["required"].as_array().unwrap();
    for field in ["linearization", "poly_a"] {
        assert!(
            required.iter().any(|value| value == field),
            "{field} is required but nullable"
        );
    }
    let profiles = call(&client, "list_ivt_profiles", json!({})).await;
    assert_eq!(profiles["engine_version"], ivt::ENGINE);
    assert_eq!(
        profiles["catalogue_revision"],
        ivt::catalogue_revision().unwrap()
    );
    assert_eq!(
        profiles["profiles"],
        serde_json::to_value(ivt::profiles().unwrap()).unwrap()
    );
    let imported = call(&client, "import_document", json!({"text":FASTA})).await;
    let document = json!({"id":imported["documents"][0]["id"],"revision":imported["documents"][0]["revision"]});
    let before = call(&client, "list_documents", json!({})).await;
    let declaration = request(&document, &profiles["profiles"][0]);
    let args = json!({"document":document,"request":declaration});
    let result = call(&client, "assess_ivt_template", args.clone()).await;
    assert_eq!(result["document"], document);
    assert_eq!(result["request"], declaration);
    assert_eq!(result["source_revision"], document["revision"]);
    assert_eq!(result["source_length"], DNA.len());
    assert_eq!(result["promoter"]["dna"], "TAATACGACTCACTATA");
    assert_eq!(
        result["initiator_context"]["dna"]
            .as_str()
            .unwrap()
            .get(..2),
        Some("AG")
    );
    assert_eq!(result["poly_a"]["adenosines"], 4);
    assert_eq!(result["poly_a"]["ambiguous_bases"], 0);
    assert_eq!(result["linearization"]["kind"], "Blunt");
    assert!(
        result.get("status").is_none(),
        "individual checks must not invent biological acceptance"
    );
    let source = collagene_core::formats::sequence::import_text(FASTA)
        .unwrap()
        .molecules
        .remove(0);
    let typed: ivt::Request = serde_json::from_value(declaration.clone()).unwrap();
    let report: ivt::Report =
        serde_json::from_value(artifact(&client, &result["artifact_id"]).await).unwrap();
    report.validate_for(&source, &typed).unwrap();
    assert_eq!(result["context_revision"], report.context_revision);
    assert_eq!(result["report_revision"], report.report_revision);
    assert_eq!(
        result["findings"],
        serde_json::to_value(&report.findings).unwrap()
    );
    assert!(
        report
            .findings
            .iter()
            .any(|finding| finding.rule == ivt::Rule::ConsensusPromoter
                && finding.verdict == ivt::Verdict::Pass)
    );
    assert!(
        report
            .findings
            .iter()
            .any(|finding| finding.rule == ivt::Rule::CapInitiator
                && finding.verdict == ivt::Verdict::Pass)
    );
    assert!(
        report
            .profile
            .sources
            .iter()
            .all(|source| source.url.starts_with("https://"))
    );
    let repeated = call(&client, "assess_ivt_template", args.clone()).await;
    assert_eq!(repeated["artifact_id"], result["artifact_id"]);
    assert_eq!(call(&client, "list_documents", json!({})).await, before);
    for path in ["document", "request"] {
        let mut changed = args.clone();
        changed[path][if path == "document" {
            "revision"
        } else {
            "source_revision"
        }] = json!("0".repeat(64));
        reject(&client, changed).await;
    }
    let mut changed = args.clone();
    changed["request"]["profile"]["identity"] = json!("0".repeat(64));
    reject(&client, changed).await;
    let mut changed = args.clone();
    changed["request"]["provenance"]["statement"] = json!("x".repeat(2049));
    reject(&client, changed).await;
    for field in ["linearization", "poly_a"] {
        let mut missing = args.clone();
        missing["request"].as_object_mut().unwrap().remove(field);
        reject(&client, missing).await;
    }
    // A genuine conflicting declaration returns typed evidence, not a tool
    // error or a silently replaced cap/initiator requirement.
    let mut conflict = args;
    conflict["request"]["initiator"] = json!("Gg");
    let conflict = call(&client, "assess_ivt_template", conflict).await;
    assert_ne!(conflict["context_revision"], result["context_revision"]);
    assert!(
        conflict["findings"]
            .as_array()
            .unwrap()
            .iter()
            .any(|finding| finding["rule"] == "ConflictingRequirements"
                && finding["verdict"] == "Fail")
    );
    assert_eq!(call(&client, "list_documents", json!({})).await, before);
    drop(client);
    server.abort();
}
