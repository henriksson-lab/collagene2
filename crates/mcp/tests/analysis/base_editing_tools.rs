//! Actual in-memory SDK coverage of conditional base-editor reference states.
use collagene_core::{
    construct_qc::{AnnotationAssumption, annotation_revision},
    sequence::crispr::{self, base_editing},
};
use collagene_mcp::server::CollageneServer;
use rmcp::{RoleClient, ServiceExt, model::CallToolRequestParams, service::RunningService};
use serde_json::{Value, json};

type Client = RunningService<RoleClient, ()>;
const GUIDE: &str = "ATGCCCCAAAAAAAAAAAAA";

async fn call(client: &Client, name: &str, arguments: Value) -> Value {
    let reply = client
        .call_tool(
            CallToolRequestParams::new(name.to_owned())
                .with_arguments(arguments.as_object().unwrap().clone()),
        )
        .await
        .unwrap();
    assert_ne!(reply.is_error, Some(true), "{name}: {:?}", reply.content);
    reply.structured_content.unwrap()
}
async fn reject(client: &Client, arguments: Value) {
    match client
        .call_tool(
            CallToolRequestParams::new("assess_base_editing")
                .with_arguments(arguments.as_object().unwrap().clone()),
        )
        .await
    {
        Ok(reply) => assert_eq!(reply.is_error, Some(true)),
        Err(rmcp::service::ServiceError::McpError(error)) => assert_eq!(error.code.0, -32602),
        Err(error) => panic!("unexpected transport error: {error}"),
    }
}
async fn artifact(client: &Client, id: &Value) -> Value {
    let mut text = String::new();
    let mut start = 0;
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
fn genbank() -> String {
    format!(
        "LOCUS       EDITOR                    23 bp    DNA     linear   SYN 28-NOV-2026\nDEFINITION  Independently declared coding fixture.\nFEATURES             Location/Qualifiers\n     CDS             1..21\n                     /codon_start=1\n                     /transl_table=1\nORIGIN\n        1 {GUIDE}agg\n//\n"
    )
}
fn cds_is_required(schema: &Value) -> bool {
    match schema {
        Value::Object(object) => {
            (object.get("properties").is_some_and(|properties| {
                properties.get("guide_options").is_some() && properties.get("cds").is_some()
            }) && object
                .get("required")
                .and_then(Value::as_array)
                .is_some_and(|required| required.iter().any(|field| field == "cds")))
                || object.values().any(cds_is_required)
        }
        Value::Array(values) => values.iter().any(cds_is_required),
        _ => false,
    }
}

#[tokio::test]
async fn base_editor_discovery_conditional_pages_codon_calls_and_artifact_are_bound_without_mutation()
 {
    let (client, server) = setup().await;
    let tools = client.list_all_tools().await.unwrap();
    assert_eq!(tools.len(), collagene_mcp::CURRENT_TOOL_COUNT);
    let tool = tools
        .iter()
        .find(|tool| tool.name == "assess_base_editing")
        .unwrap();
    assert!(
        tools
            .iter()
            .any(|tool| tool.name == "list_base_editor_profiles")
    );
    assert!(cds_is_required(
        &serde_json::to_value(&tool.input_schema).unwrap()
    ));
    let catalogue = call(&client, "list_base_editor_profiles", json!({})).await;
    assert_eq!(catalogue["engine_version"], base_editing::ENGINE_VERSION);
    assert_eq!(
        catalogue["catalogue_revision"],
        base_editing::catalogue_revision().unwrap()
    );
    assert_eq!(
        catalogue["profiles"],
        serde_json::to_value(base_editing::catalogue().unwrap()).unwrap()
    );
    assert_eq!(catalogue["profile_count"], 1);
    assert_eq!(catalogue["profiles"][0]["id"], "komor-2016-be1-c4-c8");
    assert_eq!(catalogue["profiles"][0]["window_first"], 4);
    assert_eq!(catalogue["profiles"][0]["window_last"], 8);

    let text = genbank();
    let source = collagene_core::formats::sequence::import_text(&text)
        .unwrap()
        .molecules
        .remove(0);
    let imported = call(&client, "import_document", json!({"text":text})).await;
    let document = json!({"id":imported["documents"][0]["id"],"revision":imported["documents"][0]["revision"]});
    let before = call(&client, "list_documents", json!({})).await;
    let guide_options = crispr::GuideOptions::default();
    let designed = call(
        &client,
        "design_crispr_guides",
        json!({"document":document,"options":guide_options}),
    )
    .await;
    let guide = designed["candidates"]
        .as_array()
        .unwrap()
        .iter()
        .find(|guide| guide["guide_sequence"] == GUIDE && guide["reverse"] == false)
        .unwrap();
    let request = base_editing::Request {
        guide_options,
        guide: serde_json::from_value(guide.clone()).unwrap(),
        profile: serde_json::from_value(catalogue["profiles"][0].clone()).unwrap(),
        intended_positions: vec![4],
        max_outcomes: 4096,
        cds: Some(base_editing::CdsRequest {
            annotation: AnnotationAssumption {
                feature_index: 0,
                source_revision: source.revision(),
                annotation_revision: annotation_revision(&source.features[0]).unwrap(),
                rationale: "Explicit single-interval CDS role assumption, not annotation truth"
                    .into(),
            },
            genetic_code: 1,
        }),
    };
    let args = json!({"document":document,"request":request,"max_states":2});
    let result = call(&client, "assess_base_editing", args.clone()).await;
    assert_eq!(result["document"], document);
    assert_eq!(result["source_revision"], document["revision"]);
    assert_eq!(result["request"], serde_json::to_value(&request).unwrap());
    assert_eq!(result["finite_state_count"], 16);
    assert_eq!(result["retained_state_count"], 16);
    assert_eq!(result["finite_enumeration_complete"], true);
    assert_eq!(result["biological_outcomes_exhaustive"], false);
    assert_eq!(
        result["sites"]
            .as_array()
            .unwrap()
            .iter()
            .map(|site| site["reference_position"].as_u64().unwrap())
            .collect::<Vec<_>>(),
        [3, 4, 5, 6]
    );
    assert_eq!(result["states"].as_array().unwrap().len(), 2);
    assert_eq!(result["states"][0]["edited_positions"], json!([]));
    assert_eq!(result["states"][0]["satisfies_intent"], false);
    assert_eq!(result["states"][1]["edited_positions"], json!([4]));
    assert_eq!(result["states"][1]["bystander_positions"], json!([]));
    assert_eq!(result["states"][1]["satisfies_intent"], true);
    assert_eq!(
        result["states"][1]["codon_consequences"][0]["before"],
        "CCC"
    );
    assert_eq!(result["states"][1]["codon_consequences"][0]["after"], "TCC");
    assert_eq!(result["display_truncated"], true);
    assert!(result.get("status").is_none());
    assert!(result.get("probability").is_none());
    let report: base_editing::Report =
        serde_json::from_value(artifact(&client, &result["artifact_id"]).await).unwrap();
    report.validate_for(&source, &request).unwrap();
    assert_eq!(result["context_revision"], report.context_revision);
    assert_eq!(result["report_revision"], report.report_revision);
    assert_eq!(report.states[15].bystander_positions, [5, 6, 7]);
    assert!(
        report.states[8]
            .codon_consequences
            .iter()
            .any(|consequence| consequence.before == "CAA"
                && consequence.after == "TAA"
                && consequence.kind == base_editing::ConsequenceKind::StopGained)
    );

    let page = call(
        &client,
        "assess_base_editing",
        json!({"document":document,"request":request,"state_start":14,"max_states":2}),
    )
    .await;
    assert_eq!(page["artifact_id"], result["artifact_id"]);
    assert_eq!(page["next_start"], 16);
    assert_eq!(page["done"], true);
    assert_eq!(page["states"][1]["bystander_positions"], json!([5, 6, 7]));
    let mut limited = request.clone();
    limited.max_outcomes = 1;
    let partial = call(
        &client,
        "assess_base_editing",
        json!({"document":document,"request":limited}),
    )
    .await;
    assert_eq!(partial["finite_state_count"], 16);
    assert_eq!(partial["retained_state_count"], 1);
    assert_eq!(partial["finite_enumeration_complete"], false);
    assert_eq!(partial["biological_outcomes_exhaustive"], false);
    assert_eq!(partial["states"][0]["satisfies_intent"], false);
    assert_ne!(partial["context_revision"], result["context_revision"]);
    let partial_report: base_editing::Report =
        serde_json::from_value(artifact(&client, &partial["artifact_id"]).await).unwrap();
    partial_report.validate_for(&source, &limited).unwrap();
    assert_eq!(call(&client, "list_documents", json!({})).await, before);

    for (path, value) in [
        (vec!["document", "revision"], json!("0".repeat(64))),
        (
            vec!["request", "profile", "revision"],
            json!("0".repeat(64)),
        ),
        (vec!["request", "max_outcomes"], json!(4097)),
        (vec!["max_states"], json!(101)),
        (vec!["state_start"], json!(17)),
    ] {
        let mut changed = args.clone();
        let mut field = &mut changed;
        for segment in path {
            field = &mut field[segment];
        }
        *field = value;
        reject(&client, changed).await;
    }
    let mut changed = args.clone();
    changed["request"]["guide"]["guide_positions"][0] = json!(1);
    reject(&client, changed).await;
    let mut changed = args.clone();
    changed["request"].as_object_mut().unwrap().remove("cds");
    reject(&client, changed).await;
    // Even a recomputed hash cannot authenticate an altered builtin window.
    let mut tampered = request.clone();
    tampered.profile.window_first = 5;
    tampered.profile.revision = tampered.profile.compute_revision().unwrap();
    reject(&client, json!({"document":document,"request":tampered})).await;
    assert_eq!(call(&client, "list_documents", json!({})).await, before);
    drop(client);
    server.abort();
}

#[tokio::test]
async fn reverse_guide_sites_preserve_original_guanine_to_adenine_and_undeclared_intent() {
    let (client, server) = setup().await;
    // Independent reverse complement of ATG CCCC A13 AGG.
    let text = format!(">reverse\nCCT{}GGGGCAT\n", "T".repeat(13));
    let source = collagene_core::formats::sequence::import_text(&text)
        .unwrap()
        .molecules
        .remove(0);
    let imported = call(&client, "import_document", json!({"text":text})).await;
    let document = json!({"id":imported["documents"][0]["id"],"revision":imported["documents"][0]["revision"]});
    let before = call(&client, "list_documents", json!({})).await;
    let guide_options = crispr::GuideOptions::default();
    let guides = call(
        &client,
        "design_crispr_guides",
        json!({"document":document,"options":guide_options}),
    )
    .await;
    let guide = guides["candidates"]
        .as_array()
        .unwrap()
        .iter()
        .find(|guide| guide["guide_sequence"] == GUIDE && guide["reverse"] == true)
        .unwrap();
    let request = base_editing::Request {
        guide_options,
        guide: serde_json::from_value(guide.clone()).unwrap(),
        profile: base_editing::be1_2016().unwrap(),
        intended_positions: vec![],
        max_outcomes: 16,
        cds: None,
    };
    let result = call(
        &client,
        "assess_base_editing",
        json!({"document":document,"request":request}),
    )
    .await;
    assert_eq!(
        result["sites"]
            .as_array()
            .unwrap()
            .iter()
            .map(|site| site["reference_position"].as_u64().unwrap())
            .collect::<Vec<_>>(),
        [19, 18, 17, 16]
    );
    assert!(
        result["sites"]
            .as_array()
            .unwrap()
            .iter()
            .all(|site| site["reference_from"] == "G" && site["reference_to"] == "A")
    );
    assert_eq!(result["intent_declared"], false);
    assert!(
        result["states"]
            .as_array()
            .unwrap()
            .iter()
            .all(|state| state["satisfies_intent"].is_null()
                && state["bystander_positions"] == json!([]))
    );
    let report: base_editing::Report =
        serde_json::from_value(artifact(&client, &result["artifact_id"]).await).unwrap();
    report.validate_for(&source, &request).unwrap();
    assert_eq!(call(&client, "list_documents", json!({})).await, before);
    drop(client);
    server.abort();
}
