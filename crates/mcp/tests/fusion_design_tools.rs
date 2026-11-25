//! Public SDK route: original sources → explicit fusion → real preparation.
use collagene_core::ligation_fidelity as fidelity;
use collagene_mcp::CollageneServer;
use rmcp::{RoleClient, ServiceExt, model::CallToolRequestParams, service::RunningService};
use serde_json::{Value, json};
type Client = RunningService<RoleClient, ()>;
async fn call(c: &Client, name: &str, args: Value) -> Value {
    let result = c
        .call_tool(
            CallToolRequestParams::new(name.to_string())
                .with_arguments(args.as_object().unwrap().clone()),
        )
        .await
        .unwrap();
    assert_ne!(result.is_error, Some(true), "{name}: {:?}", result.content);
    result.structured_content.unwrap()
}
async fn reject(c: &Client, name: &str, args: Value) {
    let result = c
        .call_tool(
            CallToolRequestParams::new(name.to_string())
                .with_arguments(args.as_object().unwrap().clone()),
        )
        .await;
    assert!(result.is_err() || result.unwrap().is_error == Some(true));
}
fn dataset() -> fidelity::Dataset {
    let mut dataset = fidelity::Dataset {
        schema_version: 1,
        overhang_length: 4,
        content_revision: String::new(),
        conditions: fidelity::Conditions {
            ligase: "synthetic test only".into(),
            buffer: "synthetic".into(),
            temperature_protocol: "not measured".into(),
            incubation: "not measured".into(),
            restriction_enzymes: vec![],
            substrate_context: "synthetic directed observations; no empirical authority".into(),
        },
        provenance: fidelity::Provenance {
            source_url: "https://example.org/synthetic-fixture".into(),
            citation: "Independently constructed fixture".into(),
            license: "MIT".into(),
            attribution: "Collagene 2 test fixture".into(),
            source_blake3: blake3::hash(b"synthetic directed counts")
                .to_hex()
                .to_string(),
            normalization: "literal directed counts; no symmetry inferred".into(),
        },
        cells: vec![
            fidelity::Cell {
                end_5to3: "AACC".into(),
                partner_5to3: "AACC".into(),
                count: 1,
            },
            fidelity::Cell {
                end_5to3: "AACC".into(),
                partner_5to3: "GGTT".into(),
                count: 9,
            },
            fidelity::Cell {
                end_5to3: "GGTT".into(),
                partner_5to3: "AACC".into(),
                count: 18,
            },
            fidelity::Cell {
                end_5to3: "GGTT".into(),
                partner_5to3: "GGTT".into(),
                count: 2,
            },
        ],
    };
    dataset.content_revision = fidelity::content_revision(&dataset).unwrap();
    dataset
}
#[tokio::test]
async fn sdk_search_materialization_and_condition_bound_scores_remain_advisory() {
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
    let c = ().serve(b).await.unwrap();
    assert_eq!(
        c.list_all_tools().await.unwrap().len(),
        collagene_mcp::CURRENT_TOOL_COUNT
    );
    let imported = call(
        &c,
        "import_document",
        json!({"text":">A\nATGACCCCGGGGTTTCCAGA\n>B\nATCCGTTAGCGACCTTGGAA\n"}),
    )
    .await;
    let bindings = imported["documents"]
        .as_array()
        .unwrap()
        .iter()
        .zip(["a", "b"])
        .map(|(row, id)| json!({"id":id,"document":{"id":row["id"],"revision":row["revision"]}}))
        .collect::<Vec<_>>();
    let request = json!({"base":{"base":{"name":"Independent original AB","circular":true,"fragments":[
        {"role":"left","input_id":"a","start":0,"end":20,"reverse":false},{"role":"right","input_id":"b","start":0,"end":20,"reverse":false}]},
        "edits":[],"rows":[{"DesignPcr":{"role":"left","input_id":"a","reverse":false,"forward_length":8,"reverse_length":8}},{"DesignPcr":{"role":"right","input_id":"b","reverse":false,"forward_length":8,"reverse_length":8}}],
        "enzyme":"BsaI","clamp":"ATGC","spacer":"A","terminal_right_overhang":null,"allow_nicked":false},
        "windows":[{"junction":0,"min_shift":-1,"max_shift":1}],
        "constraints":{"reject_palindromes":false,"distinct_reverse_complements":false,"minimum_hamming":0},
        "limits":{"max_states":1000,"max_candidates":32}});
    let search = call(
        &c,
        "search_golden_gate_fusions",
        json!({"inputs":bindings,"request":request}),
    )
    .await;
    assert_eq!(search["verified"], false);
    assert_eq!(search["search_complete"], true);
    let candidate = search["candidates"]
        .as_array()
        .unwrap()
        .iter()
        .find(|c| c["shifts"][0] == -1)
        .unwrap();
    let mut selection =
        json!({"context_revision":search["context_revision"],"candidate_id":candidate["id"]});
    let prepared = call(
        &c,
        "materialize_golden_gate_fusions",
        json!({"inputs":bindings,"request":request,"selection":selection}),
    )
    .await;
    assert_eq!(prepared["verified"], false);
    assert_eq!(prepared["target_revision"], search["target_revision"]);
    assert_eq!(prepared["recipe"]["fragments"].as_array().unwrap().len(), 2);
    selection["context_revision"] = json!("0".repeat(64));
    reject(
        &c,
        "materialize_golden_gate_fusions",
        json!({"inputs":bindings,"request":request,"selection":selection}),
    )
    .await;
    let catalog = call(&c, "list_ligation_profiles", json!({})).await;
    assert_eq!(catalog["profiles"].as_array().unwrap().len(), 5);
    let profile = catalog["profiles"]
        .as_array()
        .unwrap()
        .iter()
        .find(|p| p["enzyme_formulation"] == "BsaI-HFv2")
        .unwrap()["profile"]
        .clone();
    let empirical = call(
        &c,
        "score_published_ligation_profile",
        json!({"profile":profile,"ends":["AACC"],"physical_competitors":["ATGA","ATCC"]}),
    )
    .await;
    assert_eq!(empirical["verified"], false);
    assert_eq!(empirical["report"]["pairs"].as_array().unwrap().len(), 1);
    assert!(
        !empirical["report"]["physical_ends"]
            .as_array()
            .unwrap()
            .contains(&json!("GGAT"))
    );
    let ranked=call(&c,"rank_golden_gate_fusions",json!({"inputs":bindings,"request":{"design":request,"profile":profile,"enzyme_formulation":"BsaI-HFv2"}})).await;
    assert_eq!(ranked["verified"], false);
    assert_eq!(
        ranked["selection_context_revision"],
        search["context_revision"]
    );
    let ranked_preparation=call(&c,"materialize_golden_gate_fusions",json!({"inputs":bindings,"request":request,"selection":{"context_revision":ranked["selection_context_revision"],"candidate_id":ranked["ranked"][0]["candidate_id"]}})).await;
    assert_eq!(
        ranked_preparation["target_revision"],
        search["target_revision"]
    );
    let data = dataset();
    let scored = call(
        &c,
        "score_ligation_fidelity",
        json!({"dataset":data,"ends":["AACC","GGTT","AACC"]}),
    )
    .await;
    assert_eq!(scored["verified"], false);
    assert_eq!(scored["report"]["overall_fidelity"], 0.9);
    assert_eq!(scored["report"]["pairs"].as_array().unwrap().len(), 1);
    let competition = call(
        &c,
        "score_ligation_fidelity",
        json!({"dataset":data,"ends":["AACC"],"physical_competitors":[]}),
    )
    .await;
    assert_eq!(competition["report"]["overall_fidelity"], 0.9);
    assert_ne!(competition["report"]["model"], scored["report"]["model"]);
    // Legal sums of observed u64 cells can exceed JSON's exact integer type.
    let mut wide = data.clone();
    for cell in &mut wide.cells {
        cell.count = u64::MAX;
    }
    wide.content_revision = fidelity::content_revision(&wide).unwrap();
    let wide_score = call(
        &c,
        "score_ligation_fidelity",
        json!({"dataset":wide,"ends":["AACC"]}),
    )
    .await;
    assert_eq!(wide_score["report"]["overall_fidelity"], 0.5);
    assert_eq!(
        wide_score["report"]["pairs"][0]["estimate"]["Known"]["total_count"],
        (u128::from(u64::MAX) * 4).to_string()
    );
    let mut missing = data;
    missing.cells.pop();
    missing.content_revision = fidelity::content_revision(&missing).unwrap();
    let unknown = call(
        &c,
        "score_ligation_fidelity",
        json!({"dataset":missing,"ends":["AACC"]}),
    )
    .await;
    assert!(unknown["report"]["overall_fidelity"].is_null());
    c.cancel().await.unwrap();
    server.abort();
}
