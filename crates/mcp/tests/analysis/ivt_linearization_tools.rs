//! Real SDK checks: read-only source-bound alternatives and complete native evidence.
use base64::{Engine, engine::general_purpose::STANDARD};
use collagene_core::{
    chemistry::{methylation::MethylationState, tracked_duplex::GenerationAssumptions},
    cloning::restriction::{catalog, conditions, ivt_linearization as native},
    documents::transfer::{PlanningWorkspace, save_workspace},
    rna_chemistry::{Provenance, ivt},
};
use collagene_mcp::server::CollageneServer;
use collagene_model::{Molecule, Topology};
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
            CallToolRequestParams::new("assess_ivt_linearization")
                .with_arguments(args.as_object().unwrap().clone()),
        )
        .await
    {
        Ok(reply) => assert_eq!(reply.is_error, Some(true)),
        Err(rmcp::service::ServiceError::McpError(error)) => assert_eq!(error.code.0, -32602),
        Err(error) => panic!("unexpected transport error: {error}"),
    }
}
async fn import(client: &Client, molecule: &Molecule) -> Value {
    let bytes = save_workspace(&PlanningWorkspace {
        documents: vec![molecule.clone()],
        target_revision: molecule.revision(),
    })
    .unwrap();
    let imported = call(
        client,
        "import_workspace",
        json!({"data_base64":STANDARD.encode(bytes)}),
    )
    .await;
    json!({"id":imported["documents"][0]["id"],"revision":molecule.revision()})
}
async fn artifact(client: &Client, id: &Value) -> native::Report {
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
fn request(source: &Molecule, enzyme: &str, length: usize) -> native::Request {
    native::Request {
        ivt: ivt::Request {
            source_revision: source.revision(),
            profile: ivt::profiles().unwrap().remove(0).reference,
            transcript: ivt::Transcript {
                start: 17,
                length,
                coding_strand: ivt::CodingStrand::Forward,
            },
            cap: ivt::CapChoice::CleanCapAgN7113,
            initiator: ivt::InitiatorRequirement::NoAdditionalRequirement,
            linearization: None,
            poly_a: None,
            provenance: Provenance {
                statement: "Explicit conditional IVT source interval".into(),
                reference: None,
            },
        },
        source: native::SourceInput::DeclaredOriginal {
            generation: GenerationAssumptions {
                host: None,
                dam: MethylationState::Unmodified,
                dcm: MethylationState::Unmodified,
                cpg: MethylationState::Unmodified,
            },
        },
        selection: catalog::select(&catalog::builtin_catalog().unwrap(), &[enzyme.into()]).unwrap(),
        conditions: None,
        policy: native::Policy {
            allowed_end_kinds: vec![ivt::EndKind::Blunt, ivt::EndKind::FivePrimeOverhang],
            require_known_conditions: false,
            require_resolved_ivt_checks: false,
            require_complete_origin: true,
        },
        limits: native::Limits {
            max_comparisons: 1_000_000,
            max_sites_per_enzyme: 64,
            max_report_bytes: native::MAX_REPORT_BYTES,
        },
    }
}

#[tokio::test]
async fn linearization_sdk_preserves_alternative_completeness_bounds_and_source_identity() {
    let (client, server) = setup().await;
    let tools = client.list_all_tools().await.unwrap();
    assert_eq!(tools.len(), collagene_mcp::CURRENT_TOOL_COUNT);
    let tool = tools
        .iter()
        .find(|tool| tool.name == "assess_ivt_linearization")
        .unwrap();
    let schema = serde_json::to_value(&tool.input_schema).unwrap();
    let schema_text = schema.to_string();
    assert!(schema_text.contains("DeclaredOriginal") && schema_text.contains("Saved"));
    let request_schema = schema["$defs"]
        .as_object()
        .unwrap()
        .values()
        .find(|definition| {
            definition["properties"].get("selection").is_some()
                && definition["properties"].get("ivt").is_some()
        })
        .unwrap();
    assert!(
        request_schema["required"]
            .as_array()
            .unwrap()
            .iter()
            .any(|f| f == "conditions")
    );
    let source = Molecule::dna(
        "Literal conditional IVT source",
        "TAATACGACTCACTATAAGTTCCCGGGTT",
        Topology::Linear,
    )
    .unwrap();
    let document = import(&client, &source).await;
    let before = call(&client, "list_documents", json!({})).await;
    let mut query = request(&source, "smai", 7);
    query.selection = catalog::select(
        &catalog::builtin_catalog().unwrap(),
        &["smai".into(), "ecori".into()],
    )
    .unwrap();
    let args = json!({"document":document,"request":query});
    let result = call(&client, "assess_ivt_linearization", args.clone()).await;
    assert_eq!(result["alternative_count"], 2);
    assert_eq!(result["document"], document);
    assert_eq!(result["source_revision"], source.revision());
    assert_eq!(
        result["eligibility_policy"],
        serde_json::to_value(&query.policy).unwrap()
    );
    let alternatives = result["alternatives"].as_array().unwrap();
    let smai = alternatives
        .iter()
        .find(|a| a["enzyme_id"] == "smai")
        .unwrap();
    assert_eq!(smai["state"], "EligibleUnderDeclarations");
    assert_eq!(smai["scan_complete"], true);
    assert_eq!(smai["products_determined"], true);
    assert_eq!(smai["mapping"]["top_bond"], 24);
    assert_eq!(smai["mapping"]["bottom_bond"], 24);
    assert_eq!(smai["mapping"]["template_cut_offset"], 7);
    assert_eq!(smai["fragment_count"], 2);
    let ecori = alternatives
        .iter()
        .find(|a| a["enzyme_id"] == "ecori")
        .unwrap();
    assert_eq!(ecori["state"], "Rejected");
    assert_eq!(ecori["scan_complete"], true);
    assert_eq!(ecori["site_count"], 0);
    assert!(
        ecori["issues"]
            .as_array()
            .unwrap()
            .iter()
            .any(|i| i["kind"] == "NoPermittedCut")
    );
    let full = artifact(&client, &result["artifact_id"]).await;
    full.validate_for(&source, &query).unwrap();
    assert_eq!(full.report_revision, result["report_revision"]);
    assert_eq!(full.candidates.len(), 2);
    assert!(full.candidates[0].fragments.len() >= 1);
    assert!(result.get("winner").is_none() && result.get("success").is_none());
    let repeated = call(&client, "assess_ivt_linearization", args.clone()).await;
    assert_eq!(repeated["artifact_id"], result["artifact_id"]);
    let mut limited = query.clone();
    limited.limits.max_comparisons = 1;
    let limited = call(
        &client,
        "assess_ivt_linearization",
        json!({"document":document,"request":limited}),
    )
    .await;
    for candidate in limited["alternatives"].as_array().unwrap() {
        assert_eq!(candidate["state"], "Unresolved");
        assert_eq!(candidate["scan_complete"], false);
        assert!(
            !candidate["issues"]
                .as_array()
                .unwrap()
                .iter()
                .any(|i| i["kind"] == "NoPermittedCut")
        );
    }
    for path in ["document", "request"] {
        let mut stale = args.clone();
        if path == "document" {
            stale["document"]["revision"] = json!("0".repeat(64));
        } else {
            stale["request"]["ivt"]["source_revision"] = json!("0".repeat(64));
        }
        reject(&client, stale).await;
    }
    let mut stale = args.clone();
    stale["request"]["selection"]["catalog_content_hash"] = json!("0".repeat(64));
    reject(&client, stale).await;
    let mut missing = args.clone();
    missing["request"]
        .as_object_mut()
        .unwrap()
        .remove("conditions");
    reject(&client, missing).await;
    let mut too_small = args;
    too_small["request"]["limits"]["max_report_bytes"] = json!(1024);
    reject(&client, too_small).await;
    assert_eq!(call(&client, "list_documents", json!({})).await, before);
    drop(client);
    server.abort();
}

#[tokio::test]
async fn linearization_sdk_keeps_ordinary_hf_and_missing_formulation_documentary_choices_distinct()
{
    let (client, server) = setup().await;
    let source = Molecule::dna(
        "Literal EcoRI IVT source",
        "TAATACGACTCACTATAAGTTGAATTCTT",
        Topology::Linear,
    )
    .unwrap();
    let document = import(&client, &source).await;
    let before = call(&client, "list_documents", json!({})).await;
    let mut query = request(&source, "ecori", 9);
    query.policy.require_known_conditions = true;
    let mut identities = Vec::new();
    for (formulation, state, documentary) in [
        (None, "Unresolved", "Unknown"),
        (Some("neb-r0101"), "Rejected", "BelowThreshold"),
        (Some("neb-r3101"), "EligibleUnderDeclarations", "Known"),
    ] {
        query.conditions = Some(conditions::ReactionConditions {
            schema_version: 1,
            formulations: vec![conditions::ChosenFormulation {
                enzyme_id: "ecori".into(),
                formulation_id: formulation.map(str::to_owned),
            }],
            buffer: "rCutSmart".into(),
            temperature_milli_c: 37_000,
            incubation_minutes: 60,
            minimum_activity_percent: 75,
            activity_context: Some("NEB unit-assay relative activity".into()),
            temperature_context: Some("vendor recommended incubation".into()),
        });
        let result = call(
            &client,
            "assess_ivt_linearization",
            json!({"document":document,"request":query}),
        )
        .await;
        let candidate = &result["alternatives"][0];
        assert_eq!(candidate["state"], state);
        assert_eq!(candidate["conditions"]["status"], documentary);
        assert_eq!(candidate["conditions"]["star_activity"], "unassessed");
        assert_eq!(
            candidate["conditions"]["incubation_duration_compatibility"],
            "unassessed"
        );
        assert_eq!(candidate["mapping"]["coding_cut_offset"], 5);
        assert_eq!(candidate["mapping"]["template_cut_offset"], 9);
        let report = artifact(&client, &result["artifact_id"]).await;
        report.validate_for(&source, &query).unwrap();
        assert_eq!(
            serde_json::to_value(&report.request.conditions).unwrap(),
            result["reaction_choices"]
        );
        identities.push(report.context_revision);
    }
    assert!(
        identities[0] != identities[1]
            && identities[1] != identities[2]
            && identities[0] != identities[2]
    );
    assert_eq!(call(&client, "list_documents", json!({})).await, before);
    drop(client);
    server.abort();
}

#[tokio::test]
async fn linearization_sdk_replays_saved_pcr_and_refuses_original_relabeling_or_foreign_event() {
    use collagene_core::{
        documents::workspace::Workspace,
        material::{
            capsule,
            lineage::{EventBinding, InputEventBinding},
            outputs, run,
        },
        planning::plan::{DesignPlan, Operation, PlanInput, Source},
        sequence::analysis::Primer,
    };
    let original = Molecule::dna(
        "Original PCR template",
        "TAATACGACTCACTATAAGTTCCCGGGTT",
        Topology::Linear,
    )
    .unwrap();
    let mut target = original.clone();
    target.name = "Independent PCR target".into();
    let plan = DesignPlan {
        standard: None,
        name: "Saved literal PCR".into(),
        target_revision: target.revision(),
        inputs: vec![PlanInput {
            id: "original".into(),
            revision: original.revision(),
            molecule: original.clone(),
        }],
        steps: vec![Operation::Pcr {
            template: Source::Input("original".into()),
            forward: Primer::new("literal first8", "TAATACGA", "").unwrap(),
            reverse: Primer::new("literal reverse last8", "AACCCGGG", "").unwrap(),
            template_strand: collagene_core::chemistry::Strand::Top,
        }],
    };
    let workspace =
        Workspace::from_documents("Actual saved source", vec![original, target], 1).unwrap();
    let event = |index: usize| {
        let document = &workspace.documents[index];
        EventBinding {
            document_id: document.id.clone(),
            event_id: document.head_event.clone().unwrap(),
            revision: document.head.clone(),
        }
    };
    let committed = outputs::finalize(
        &workspace,
        &plan,
        event(1),
        vec![InputEventBinding {
            input_id: "original".into(),
            event: event(0),
        }],
        &outputs::Options {
            name: "Saved actual PCR".into(),
            ..Default::default()
        },
    )
    .unwrap();
    let source = committed
        .workspace
        .head(&committed.receipt.construct.document_id)
        .unwrap();
    let origin = capsule::capsule(
        &committed.workspace,
        &run::Selection {
            run_id: committed.receipt.run.run_id.clone(),
            step_index: 0,
            event: committed.receipt.construct.clone(),
        },
    )
    .unwrap();
    let (client, server) = setup().await;
    let document = import(&client, source).await;
    let before = call(&client, "list_documents", json!({})).await;
    let mut query = request(source, "smai", 7);
    let original_route = query.source.clone();
    query.source = native::SourceInput::Saved {
        origin: run::MaterialOrigin::Run(origin.clone()),
    };
    let result = call(
        &client,
        "assess_ivt_linearization",
        json!({"document":document,"request":query}),
    )
    .await;
    assert_eq!(result["source_route"], "fresh_saved_origin_replay");
    assert_eq!(
        result["alternatives"][0]["state"],
        "EligibleUnderDeclarations"
    );
    let report = artifact(&client, &result["artifact_id"]).await;
    report.validate_for(source, &query).unwrap();
    assert!(report.source_modifications.generation.is_empty());
    assert!(!report.source_history_incomplete);
    query.source = original_route;
    reject(&client, json!({"document":document,"request":query})).await;
    let mut foreign = origin;
    foreign.selection.event.revision = "0".repeat(64);
    query.source = native::SourceInput::Saved {
        origin: run::MaterialOrigin::Run(foreign),
    };
    reject(&client, json!({"document":document,"request":query})).await;
    assert_eq!(call(&client, "list_documents", json!({})).await, before);
    drop(client);
    server.abort();
}
