use collagene_mcp::CollageneServer;
use rmcp::{ServiceExt, model::CallToolRequestParams};
use serde_json::json;

// Real SDK client/server negotiation and JSON-RPC routing over an async byte stream.
#[tokio::test]
async fn candidate_verification_and_artifact_workflow() {
    let (server_io, client_io) = tokio::io::duplex(64 * 1024);
    let server = tokio::spawn(async {
        CollageneServer::default()
            .serve(server_io)
            .await
            .unwrap()
            .waiting()
            .await
            .unwrap()
    });
    let client = ().serve(client_io).await.unwrap();
    let tools = client.list_all_tools().await.unwrap();
    assert_eq!(tools.len(), collagene_mcp::CURRENT_TOOL_COUNT);
    let plan_tool = tools.iter().find(|t| t.name == "propose_plan").unwrap();
    assert_eq!(plan_tool.input_schema["additionalProperties"], false);
    let schema = serde_json::to_string(&plan_tool.input_schema).unwrap();
    assert!(schema.contains("Overlap") && schema.contains("Pcr") && schema.contains("target"));
    macro_rules! call {
        ($name:expr,$args:expr) => {{
            client
                .call_tool(
                    CallToolRequestParams::new($name)
                        .with_arguments($args.as_object().unwrap().clone()),
                )
                .await
                .unwrap()
                .structured_content
                .unwrap()
        }};
    }
    let imported = call!(
        "import_document",
        json!({"text":">left\nAAAACCCC\n>right\nCCCCGGGG\n>target\nAAAACCCCGGGG\n"})
    );
    let docs = imported["documents"].as_array().unwrap();
    let reference = |name: &str| {
        let m = docs.iter().find(|m| m["name"] == name).unwrap();
        json!({"id":m["id"],"revision":m["revision"]})
    };
    let region = call!(
        "read_sequence_region",
        json!({"document":reference("target"),"start":4,"end":8})
    );
    assert_eq!(region["sequence"], "CCCC");
    let document_artifact = call!("export_document", reference("target"));
    let chunk = call!(
        "read_artifact",
        json!({"artifact_id":document_artifact["artifact_id"],"start":0,"max_bytes":16384})
    );
    assert_eq!(chunk["done"], true);
    assert_eq!(chunk["encoding"], "base64");
    let native = base64::Engine::decode(
        &base64::engine::general_purpose::STANDARD,
        chunk["text"].as_str().unwrap(),
    )
    .unwrap();
    let restored = collagene_core::protocol::load_document(&native).unwrap();
    assert_eq!(restored.revision(), reference("target")["revision"]);
    assert_eq!(restored.sequence.as_str(), "AAAACCCCGGGG");
    let mut stale = reference("target");
    stale["revision"] = json!("stale");
    assert!(
        client
            .call_tool(
                CallToolRequestParams::new("read_sequence_region").with_arguments(
                    json!({"document":stale,"start":0,"end":1})
                        .as_object()
                        .unwrap()
                        .clone()
                )
            )
            .await
            .is_err()
    );
    let proposal = json!({"name":"Independent expected assembly","target":reference("target"),"inputs":[{"id":"left","document":reference("left")},{"id":"right","document":reference("right")}],"steps":[{"Overlap":{"fragments":[{"Input":"left"},{"Input":"right"}],"overlaps":[4],"circular":false}}]});
    let mut forged = proposal.clone();
    forged["verified"] = json!(true);
    let forged_result = client
        .call_tool(
            CallToolRequestParams::new("propose_plan")
                .with_arguments(forged.as_object().unwrap().clone()),
        )
        .await;
    assert!(forged_result.is_err() || forged_result.unwrap().is_error == Some(true));
    let candidate = call!("propose_plan", proposal.clone());
    assert_eq!(candidate["state"], "candidate");
    assert_eq!(candidate["verified"], false);
    let duplicate = call!("propose_plan", proposal.clone());
    assert_eq!(candidate, duplicate);
    let resources = client.list_all_resources().await.unwrap();
    assert!(resources.iter().any(|r| r.uri
        == format!(
            "collagene://plan/{}",
            candidate["plan_revision"].as_str().unwrap()
        )));
    let resource = client
        .read_resource(rmcp::model::ReadResourceRequestParams::new(format!(
            "collagene://document/{}",
            reference("target")["id"].as_str().unwrap()
        )))
        .await
        .unwrap();
    let encoded = serde_json::to_value(resource).unwrap();
    let summary: serde_json::Value =
        serde_json::from_str(encoded["contents"][0]["text"].as_str().unwrap()).unwrap();
    assert_eq!(summary["length"], 12);
    assert!(summary.get("sequence").is_none());
    assert!(
        client
            .read_resource(rmcp::model::ReadResourceRequestParams::new(
                "file:///etc/passwd"
            ))
            .await
            .is_err()
    );
    let bound = json!({"plan_revision":candidate["plan_revision"]});
    let report = call!("verify_plan", bound.clone());
    assert_eq!(report["status"], "ConsistentUnderAssumptions");
    assert_eq!(report["exhaustive_alternative_search"], false);
    assert_eq!(report["products"][0]["length"], 12);
    let artifact = call!("export_protocol", bound.clone());
    let html = call!(
        "read_artifact",
        json!({"artifact_id":artifact["artifact_id"],"start":0,"max_bytes":16384})
    );
    assert_eq!(html["done"], true);
    let text = html["text"].as_str().unwrap();
    assert!(
        text.contains("ConsistentUnderAssumptions")
            && text.contains("laboratory conditions")
            && text.contains(candidate["plan_revision"].as_str().unwrap())
    );
    // Byte pagination cannot split a Unicode scalar or loop without making progress.
    let unicode_start = text.find('·').unwrap();
    assert!(client.call_tool(CallToolRequestParams::new("read_artifact").with_arguments(json!({"artifact_id":artifact["artifact_id"],"start":unicode_start+1,"max_bytes":4}).as_object().unwrap().clone())).await.is_err());
    assert!(client.call_tool(CallToolRequestParams::new("read_artifact").with_arguments(json!({"artifact_id":artifact["artifact_id"],"start":unicode_start,"max_bytes":1}).as_object().unwrap().clone())).await.is_err());
    let unicode = call!(
        "read_artifact",
        json!({"artifact_id":artifact["artifact_id"],"start":unicode_start,"max_bytes":2})
    );
    assert_eq!(unicode["text"], "·");
    assert_eq!(unicode["next_start"], unicode_start + 2);
    let mut wrong = proposal.clone();
    wrong["steps"][0]["Overlap"]["overlaps"] = json!([3]);
    let wrong = call!("propose_plan", wrong);
    let wrong_report = call!(
        "verify_plan",
        json!({"plan_revision":wrong["plan_revision"]})
    );
    assert_eq!(wrong_report["status"], "Invalid");
    // Export recomputes the wrong candidate too and keeps its invalid status visible.
    let wrong_artifact = call!(
        "export_protocol",
        json!({"plan_revision":wrong["plan_revision"]})
    );
    let wrong_html = call!(
        "read_artifact",
        json!({"artifact_id":wrong_artifact["artifact_id"],"start":0,"max_bytes":16384})
    );
    assert!(
        wrong_html["text"]
            .as_str()
            .unwrap()
            .contains("<strong>Invalid</strong>")
    );
    client.cancel().await.unwrap();
    server.await.unwrap();
}

#[tokio::test]
async fn isolated_sessions_and_atomic_imports() {
    let (a, b) = tokio::io::duplex(65536);
    let task = tokio::spawn(async {
        CollageneServer::default()
            .serve(a)
            .await
            .unwrap()
            .waiting()
            .await
            .unwrap()
    });
    let client = ().serve(b).await.unwrap();
    let args = json!({"text":">valid\nACGT\n>invalid\nXX!\n"});
    assert!(
        client
            .call_tool(
                CallToolRequestParams::new("import_document")
                    .with_arguments(args.as_object().unwrap().clone())
            )
            .await
            .is_err()
    );
    let listed = client
        .call_tool(CallToolRequestParams::new("list_documents"))
        .await
        .unwrap()
        .structured_content
        .unwrap();
    assert_eq!(listed["documents"], json!([]));
    client.cancel().await.unwrap();
    task.await.unwrap();
}

#[tokio::test]
async fn scientific_tools_use_core_models_and_bound_output() {
    let (a, b) = tokio::io::duplex(65536);
    let task = tokio::spawn(async {
        CollageneServer::default()
            .serve(a)
            .await
            .unwrap()
            .waiting()
            .await
            .unwrap()
    });
    let client = ().serve(b).await.unwrap();
    macro_rules! call {
        ($name:expr,$args:expr) => {{
            client
                .call_tool(
                    CallToolRequestParams::new($name)
                        .with_arguments($args.as_object().unwrap().clone()),
                )
                .await
                .unwrap()
                .structured_content
                .unwrap()
        }};
    }
    let imported = call!(
        "import_document",
        json!({"text":">template\nACGACCCCGGGGTAAC\n>restriction\nAAAAGAATTCTTTT\n"})
    );
    let docs = imported["documents"].as_array().unwrap();
    let reference = |name: &str| {
        let m = docs.iter().find(|m| m["name"] == name).unwrap();
        json!({"id":m["id"],"revision":m["revision"]})
    };
    let digest = call!(
        "analyze_restriction",
        json!({"document":reference("restriction"),"enzyme":"EcoRI"})
    );
    assert_eq!(digest["site_count"], 1);
    assert_eq!(digest["fragment_lengths"], json!([9, 5]));
    let primer = |name: &str, annealing: &str, tail: &str| json!({"name":name,"annealing":{"alphabet":"Dna","bases":annealing},"tail":tail});
    let pcr = call!(
        "simulate_pcr",
        json!({"document":reference("template"),"forward":primer("F","ACGA","CGT"),"reverse":primer("R","GTTA","GAC")})
    );
    assert_eq!(pcr["product_count"], 1);
    assert_eq!(pcr["products"][0]["length"], 22);
    assert_eq!(
        pcr["products"][0]["sequence_hash"],
        blake3::hash(b"CGTACGACCCCGGGGTAACGTC").to_hex().as_str()
    );
    let invalid = client
        .call_tool(
            CallToolRequestParams::new("read_sequence_region").with_arguments(
                json!({"document":reference("template"),"start":0,"end":10001})
                    .as_object()
                    .unwrap()
                    .clone(),
            ),
        )
        .await;
    assert!(invalid.is_err());
    // A new server session cannot read a document imported into this one.
    let (c, d) = tokio::io::duplex(65536);
    let second_task = tokio::spawn(async {
        CollageneServer::default()
            .serve(c)
            .await
            .unwrap()
            .waiting()
            .await
            .unwrap()
    });
    let second = ().serve(d).await.unwrap();
    assert!(
        second
            .call_tool(
                CallToolRequestParams::new("read_sequence_region").with_arguments(
                    json!({"document":reference("template"),"start":0,"end":1})
                        .as_object()
                        .unwrap()
                        .clone()
                )
            )
            .await
            .is_err()
    );
    second.cancel().await.unwrap();
    second_task.await.unwrap();
    client.cancel().await.unwrap();
    task.await.unwrap();
}

#[tokio::test]
async fn binary_snapshot_to_candidate_preserves_revisions_and_has_no_trusted_verdict() {
    use base64::{Engine as _, engine::general_purpose::STANDARD};
    use collagene_core::{Molecule, Topology, transfer, verification};
    let mut left = Molecule::dna("left", "AAAACCCC", Topology::Linear).unwrap();
    left.metadata
        .insert("provenance".into(), "retained browser metadata".into());
    let right = Molecule::dna("right", "CCCCGGGG", Topology::Linear).unwrap();
    let target = Molecule::dna("target", "AAAACCCCGGGG", Topology::Linear).unwrap();
    let workspace = transfer::PlanningWorkspace {
        target_revision: target.revision(),
        documents: vec![left.clone(), right.clone(), target.clone()],
    };
    let snapshot = transfer::save_workspace(&workspace).unwrap();
    let (a, b) = tokio::io::duplex(65536);
    let task = tokio::spawn(async {
        CollageneServer::default()
            .serve(a)
            .await
            .unwrap()
            .waiting()
            .await
            .unwrap()
    });
    let client = ().serve(b).await.unwrap();
    macro_rules! call {
        ($name:expr,$args:expr) => {{
            client
                .call_tool(
                    CallToolRequestParams::new($name)
                        .with_arguments($args.as_object().unwrap().clone()),
                )
                .await
                .unwrap()
                .structured_content
                .unwrap()
        }};
    }
    let imported = call!(
        "import_workspace",
        json!({"data_base64":STANDARD.encode(snapshot)})
    );
    assert_eq!(imported["target"]["revision"], target.revision());
    let reference = |m: &Molecule| json!({"id":m.revision(),"revision":m.revision()});
    let candidate = call!(
        "propose_plan",
        json!({"name":"Binary exchange","target":reference(&target),"inputs":[{"id":"left","document":reference(&left)},{"id":"right","document":reference(&right)}],"steps":[{"Overlap":{"fragments":[{"Input":"left"},{"Input":"right"}],"overlaps":[4],"circular":false}}]})
    );
    let exported = call!(
        "export_candidate",
        json!({"plan_revision":candidate["plan_revision"]})
    );
    assert_eq!(exported["encoding"], "base64");
    assert_eq!(exported["state"], "candidate");
    let mut encoded = String::new();
    let mut start = 0;
    loop {
        let chunk = call!(
            "read_artifact",
            json!({"artifact_id":exported["artifact_id"],"start":start,"max_bytes":100})
        );
        assert_eq!(chunk["encoding"], "base64");
        encoded.push_str(chunk["text"].as_str().unwrap());
        if chunk["done"] == true {
            break;
        }
        start = chunk["next_start"].as_u64().unwrap();
    }
    let package = transfer::load_candidate(&STANDARD.decode(encoded).unwrap()).unwrap();
    assert_eq!(package.plan_revision, candidate["plan_revision"]);
    assert_eq!(
        package.plan.inputs[0].molecule.metadata["provenance"],
        "retained browser metadata"
    );
    let plan = package.bind(&target, &workspace.documents).unwrap();
    assert_eq!(
        verification::verify(&plan, &target).status,
        verification::Status::ConsistentUnderAssumptions
    );
    let mut edited = workspace.documents.clone();
    edited[0] = edited[0].replace(0, 1, "T").unwrap();
    assert!(package.bind(&target, &edited).is_err());
    client.cancel().await.unwrap();
    task.await.unwrap();
}

#[tokio::test]
async fn restriction_ligation_tools_keep_nicks_and_strand_geometry_explicit() {
    let (a, b) = tokio::io::duplex(65536);
    let task = tokio::spawn(async {
        CollageneServer::default()
            .serve(a)
            .await
            .unwrap()
            .waiting()
            .await
            .unwrap()
    });
    let client = ().serve(b).await.unwrap();
    macro_rules! call {
        ($name:expr,$args:expr) => {{
            client
                .call_tool(
                    CallToolRequestParams::new($name)
                        .with_arguments($args.as_object().unwrap().clone()),
                )
                .await
                .unwrap()
                .structured_content
                .unwrap()
        }};
    }
    let imported = call!(
        "import_document",
        json!({"text":include_str!("../../../fixtures/restriction-demo.gb")})
    );
    let docs = imported["documents"].as_array().unwrap();
    let reference = |name: &str| {
        let m = docs.iter().find(|m| m["name"] == name).unwrap();
        json!({"id":m["id"],"revision":m["revision"]})
    };
    let digest = call!(
        "digest_fragments",
        json!({"document":reference("vector_fixture"),"enzymes":["EcoRI"]})
    );
    assert_eq!(digest["fragment_count"], 1);
    assert_eq!(digest["fragments"][0]["duplex"]["top_length"], 14);
    assert_eq!(
        digest["fragments"][0]["duplex"]["ends"][0]["overhang"],
        "AATT"
    );
    assert_eq!(
        digest["fragments"][0]["duplex"]["ends"][0]["chemistry"]["five_prime"],
        "Phosphate"
    );
    let proposal = json!({"name":"Nicked restriction ligation","target":reference("ligation_target"),"inputs":[{"id":"vector","document":reference("vector_fixture")},{"id":"insert","document":reference("insert_fixture")}],"steps":[{"Digest":{"template":{"Input":"vector"},"enzymes":["EcoRI"],"fragment_index":0}},{"Digest":{"template":{"Input":"insert"},"enzymes":["EcoRI"],"fragment_index":0}},{"TreatEnds":{"template":{"Step":0},"treatment":"Dephosphorylate"}},{"Ligate":{"fragments":[{"source":{"Step":2},"reverse":false},{"source":{"Step":1},"reverse":false}],"circular":true,"allow_nicked":true}}]});
    let candidate = call!("propose_plan", proposal.clone());
    let bound = json!({"plan_revision":candidate["plan_revision"]});
    let report = call!("verify_plan", bound.clone());
    assert_eq!(report["status"], "ConsistentUnderAssumptions", "{report}");
    assert_eq!(
        report["products"][3]["duplex"]["nicks"]
            .as_array()
            .unwrap()
            .len(),
        2
    );
    assert_eq!(
        report["products"][3]["duplex"]["covalently_closed_in_model"],
        false
    );
    let artifact = call!("export_protocol", bound);
    let html = call!(
        "read_artifact",
        json!({"artifact_id":artifact["artifact_id"],"start":0,"max_bytes":16384})
    );
    assert!(html["text"].as_str().unwrap().contains("2 retained nicks"));
    let mut forbidden = proposal.clone();
    forbidden["steps"][3]["Ligate"]["allow_nicked"] = json!(false);
    let forbidden = call!("propose_plan", forbidden);
    let failure = call!(
        "verify_plan",
        json!({"plan_revision":forbidden["plan_revision"]})
    );
    assert_eq!(failure["status"], "Invalid");
    assert_eq!(failure["products"].as_array().unwrap().len(), 3);
    client.cancel().await.unwrap();
    task.await.unwrap();
}
