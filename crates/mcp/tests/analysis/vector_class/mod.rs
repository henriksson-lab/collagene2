//! In-memory SDK coverage using authored synthetic elements, never viral DNA.
use collagene_core::construct_qc::{
    self, AnnotationAssumption, Provenance, ReferenceElement,
    required_parts::{RequiredOrientation, RequiredPart},
    vector_class,
};
use collagene_mcp::server::CollageneServer;
use collagene_model::{Molecule, Sequence};
use rmcp::{RoleClient, ServiceExt, model::CallToolRequestParams, service::RunningService};
use serde_json::{Value, json};

type Client = RunningService<RoleClient, ()>;
const GENBANK: &str = "LOCUS       SYNVEC                    28 bp    DNA     linear   SYN 28-NOV-2026\nDEFINITION  Authored synthetic element graph, not biological vector DNA.\nFEATURES             Location/Qualifiers\n     misc_feature    3..6\n                     /label=first\n     misc_feature    9..12\n                     /label=regulatory\n     misc_feature    15..18\n                     /label=payload\n     misc_feature    complement(21..24)\n                     /label=last\nORIGIN\n        1 ggaaccggcgtattatgcccggttaaaa\n//\n";

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
            CallToolRequestParams::new("assess_vector_class")
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
fn reference(id: &str, dna: &str) -> ReferenceElement {
    let sequence = Sequence::dna(dna).unwrap();
    ReferenceElement {
        id: id.into(),
        sequence_revision: construct_qc::reference_revision(&sequence).unwrap(),
        sequence,
        provenance: Provenance {
            source_url: "https://example.org/authored-synthetic-elements".into(),
            accession: format!("synthetic-{id}-1"),
            retrieved_on: "2026-11-28".into(),
            license: "CC0 authored non-biological fixture".into(),
        },
    }
}
fn request(source: &Molecule) -> vector_class::Request {
    vector_class::Request {
        source_revision: source.revision(),
        profile: vector_class::catalogue()
            .unwrap()
            .into_iter()
            .find(|profile| profile.reference.id == vector_class::profiles::RAAV_ID)
            .unwrap()
            .reference,
        annotations: source
            .features
            .iter()
            .enumerate()
            .map(|(feature_index, feature)| AnnotationAssumption {
                feature_index,
                source_revision: source.revision(),
                annotation_revision: construct_qc::annotation_revision(feature).unwrap(),
                rationale:
                    "Explicit synthetic role/coordinate assumption, not biological annotation truth"
                        .into(),
            })
            .collect(),
        references: vec![
            reference("repeat", "AACC"),
            reference("regulatory", "CGTA"),
            reference("payload", "ATGC"),
        ],
        bindings: [
            ("left", "five-prime-itr", 0, "repeat"),
            ("regulatory", "promoter", 1, "regulatory"),
            ("payload", "payload", 2, "payload"),
            ("right", "three-prime-itr", 3, "repeat"),
        ]
        .into_iter()
        .map(
            |(id, role_id, feature_index, reference_id)| vector_class::RoleBinding {
                id: id.into(),
                role_id: role_id.into(),
                feature_index,
                reference_id: reference_id.into(),
            },
        )
        .collect(),
        boundaries: Some(vector_class::BoundarySelection {
            start_binding: "left".into(),
            end_binding: "right".into(),
            reverse: false,
        }),
        spacing: vec![],
        literal_parts: vec![RequiredPart {
            id: "excluded-synthetic-literal".into(),
            role: "User-excluded literal, not biological role".into(),
            reference_id: "repeat".into(),
            minimum_copies: 0,
            maximum_copies: Some(0),
            orientation: RequiredOrientation::Either,
        }],
        literal_part_limits: Default::default(),
    }
}
fn boundaries_required(schema: &Value) -> bool {
    match schema {
        Value::Object(object) => {
            (object.get("properties").is_some_and(|properties| {
                properties.get("bindings").is_some() && properties.get("boundaries").is_some()
            }) && object
                .get("required")
                .and_then(Value::as_array)
                .is_some_and(|required| required.iter().any(|name| name == "boundaries")))
                || object.values().any(boundaries_required)
        }
        Value::Array(values) => values.iter().any(boundaries_required),
        _ => false,
    }
}

#[tokio::test]
async fn vector_profile_discovery_actual_typed_pages_and_full_artifact_preserve_immutable_premises()
{
    let (client, server) = setup().await;
    let tools = client.list_all_tools().await.unwrap();
    assert_eq!(tools.len(), collagene_mcp::CURRENT_TOOL_COUNT);
    assert!(
        tools
            .iter()
            .any(|tool| tool.name == "list_vector_class_profiles")
    );
    let tool = tools
        .iter()
        .find(|tool| tool.name == "assess_vector_class")
        .unwrap();
    assert!(boundaries_required(
        &serde_json::to_value(&tool.input_schema).unwrap()
    ));
    let catalogue = call(&client, "list_vector_class_profiles", json!({})).await;
    assert_eq!(catalogue["engine_version"], vector_class::ENGINE_VERSION);
    assert_eq!(
        catalogue["catalogue_revision"],
        vector_class::catalogue_revision().unwrap()
    );
    assert_eq!(
        catalogue["profiles"],
        serde_json::to_value(vector_class::catalogue().unwrap()).unwrap()
    );
    assert_eq!(catalogue["profile_count"], 2);
    assert_eq!(catalogue["statement"], construct_qc::STATEMENT);

    let source = collagene_core::formats::sequence::import_text(GENBANK)
        .unwrap()
        .molecules
        .remove(0);
    let imported = call(&client, "import_document", json!({"text":GENBANK})).await;
    let document = json!({"id":imported["documents"][0]["id"],"revision":imported["documents"][0]["revision"]});
    let before = call(&client, "list_documents", json!({})).await;
    let request = request(&source);
    let args = json!({"document":document,"request":request,"max_findings":2});
    let result = call(&client, "assess_vector_class", args.clone()).await;
    assert_eq!(result["document"], document);
    assert_eq!(result["source_revision"], document["revision"]);
    assert_eq!(result["source_feature_count"], 4);
    assert_eq!(
        result["profile"],
        serde_json::to_value(&request.profile).unwrap()
    );
    assert_eq!(result["boundaries"]["outer_bases"], 22);
    assert_eq!(result["boundaries"]["interior_bases"], 14);
    assert_eq!(result["source_length"], 28);
    assert_eq!(
        result["boundaries"]["source_spans"],
        json!([{"start":2,"end":24,"reverse":false}])
    );
    assert_eq!(result["bindings"][0]["reference_verdict"], "Pass");
    assert_eq!(result["bindings"][0]["premise"], "Pass");
    assert_eq!(
        result["bindings"][3]["interval"],
        json!({"start":20,"end":24,"reverse":true})
    );
    assert_eq!(result["literal_parts"][0]["verdict"], "Fail");
    assert_eq!(result["literal_parts"][0]["certain_copies"], 2);
    assert_eq!(result["literal_searches"][0]["complete"], true);
    assert_eq!(result["numeric_context"][0]["minimum"], 4700);
    assert_eq!(result["numeric_context"][0]["maximum"], 8700);
    assert!(result.get("status").is_none());
    assert!(result.get("biological_pass").is_none());
    assert_eq!(result["findings"].as_array().unwrap().len(), 2);
    assert_eq!(result["display_truncated"], true);
    let report: vector_class::Report =
        serde_json::from_value(artifact(&client, &result["artifact_id"]).await).unwrap();
    report.validate_for(&source, &request).unwrap();
    assert_eq!(result["context_revision"], report.context_revision);
    assert_eq!(result["report_revision"], report.report_revision);
    assert_eq!(result["annotations_revision"], report.annotations_revision);
    assert_eq!(result["assumptions_revision"], report.assumptions_revision);
    assert_eq!(result["references_revision"], report.references_revision);
    assert_eq!(
        result["findings"],
        serde_json::to_value(&report.findings[..2]).unwrap()
    );
    assert!(
        report
            .findings
            .iter()
            .filter(|finding| finding.rule == vector_class::Rule::EmpiricalContext)
            .all(|finding| finding.verdict == construct_qc::Verdict::Informational)
    );
    assert_eq!(
        result["literal_findings"],
        serde_json::to_value(&report.literal_findings).unwrap()
    );
    let page = call(
        &client,
        "assess_vector_class",
        json!({"document":document,"request":request,"finding_start":2,"max_findings":100}),
    )
    .await;
    assert_eq!(page["artifact_id"], result["artifact_id"]);
    assert_eq!(page["done"], true);
    assert_eq!(page["finding_start"], 2);
    assert_eq!(
        page["findings"],
        serde_json::to_value(&report.findings[2..]).unwrap()
    );

    let mut limited = request.clone();
    limited.literal_part_limits.max_work = 1;
    let partial = call(
        &client,
        "assess_vector_class",
        json!({"document":document,"request":limited}),
    )
    .await;
    assert_eq!(partial["literal_parts"][0]["verdict"], "Indeterminate");
    assert_eq!(partial["literal_searches"][0]["complete"], false);
    assert!(partial["literal_search_work"].as_u64().unwrap() <= 1);
    assert_ne!(partial["context_revision"], result["context_revision"]);
    let partial_report: vector_class::Report =
        serde_json::from_value(artifact(&client, &partial["artifact_id"]).await).unwrap();
    partial_report.validate_for(&source, &limited).unwrap();
    let mut mismatch = request.clone();
    mismatch.references[0] = reference("repeat", "AAAA");
    let checked = call(
        &client,
        "assess_vector_class",
        json!({"document":document,"request":mismatch}),
    )
    .await;
    assert_eq!(checked["bindings"][0]["reference_verdict"], "Fail");
    assert_eq!(checked["bindings"][0]["premise"], "Indeterminate");
    assert!(checked["boundaries"].is_null());
    let checked_report: vector_class::Report =
        serde_json::from_value(artifact(&client, &checked["artifact_id"]).await).unwrap();
    checked_report.validate_for(&source, &mismatch).unwrap();
    assert!(checked_report.findings.iter().any(|finding| finding.rule
        == vector_class::Rule::ReferenceIntegrity
        && finding.verdict == construct_qc::Verdict::Fail
        && finding.binding_ids == ["left"]));

    for (path, value) in [
        (vec!["document", "revision"], json!("0".repeat(64))),
        (vec!["request", "source_revision"], json!("0".repeat(64))),
        (
            vec!["request", "profile", "identity"],
            json!("0".repeat(64)),
        ),
        (vec!["max_findings"], json!(0)),
        (vec!["max_findings"], json!(101)),
        (vec!["finding_start"], json!(report.findings.len() + 1)),
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
    changed["request"]
        .as_object_mut()
        .unwrap()
        .remove("boundaries");
    reject(&client, changed).await;
    let mut changed = args;
    changed["request"]["bindings"] = json!(vec![request.bindings[0].clone(); 65]);
    reject(&client, changed).await;
    assert_eq!(call(&client, "list_documents", json!({})).await, before);
    drop(client);
    server.abort();
}

#[tokio::test]
async fn missing_annotation_reference_and_boundary_declarations_are_indeterminate_not_absent() {
    let (client, server) = setup().await;
    let source = collagene_core::formats::sequence::import_text(GENBANK)
        .unwrap()
        .molecules
        .remove(0);
    let imported = call(&client, "import_document", json!({"text":GENBANK})).await;
    let document = json!({"id":imported["documents"][0]["id"],"revision":imported["documents"][0]["revision"]});
    let before = call(&client, "list_documents", json!({})).await;
    let mut declaration = request(&source);
    declaration.bindings.clear();
    declaration.annotations.clear();
    declaration.references.clear();
    declaration.boundaries = None;
    declaration.literal_parts.clear();
    let result = call(
        &client,
        "assess_vector_class",
        json!({"document":document,"request":declaration,"max_findings":100}),
    )
    .await;
    assert!(result["boundaries"].is_null());
    assert_eq!(result["bindings"], json!([]));
    assert!(
        result["findings"]
            .as_array()
            .unwrap()
            .iter()
            .filter(|finding| finding["rule"] == "DeclaredMultiplicity")
            .all(|finding| finding["verdict"] == "Indeterminate")
    );
    assert!(result["literal_parts"].as_array().unwrap().is_empty());
    let report: vector_class::Report =
        serde_json::from_value(artifact(&client, &result["artifact_id"]).await).unwrap();
    report.validate_for(&source, &declaration).unwrap();
    let mut stale = request(&source);
    stale.annotations[0].annotation_revision = "0".repeat(64);
    let result = call(
        &client,
        "assess_vector_class",
        json!({"document":document,"request":stale}),
    )
    .await;
    assert_eq!(result["bindings"][0]["premise"], "Indeterminate");
    assert_eq!(result["bindings"][0]["reference_verdict"], "Indeterminate");
    assert!(result["boundaries"].is_null());
    assert_eq!(call(&client, "list_documents", json!({})).await, before);
    drop(client);
    server.abort();
}
