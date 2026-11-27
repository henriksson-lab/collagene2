//! Actual SDK over memory pipes; no browser, HTTP, user files or scientific process.
use collagene_mcp::CollageneServer;
use rmcp::{RoleClient, ServiceExt, model::CallToolRequestParams, service::RunningService};
use serde_json::{Value, json};
type Client = RunningService<RoleClient, ()>;

async fn call(c: &Client, name: &str, args: Value) -> Value {
    let r = c
        .call_tool(
            CallToolRequestParams::new(name.to_owned())
                .with_arguments(args.as_object().unwrap().clone()),
        )
        .await
        .unwrap();
    assert_ne!(r.is_error, Some(true), "{name}: {:?}", r.content);
    r.structured_content.unwrap()
}
async fn reject(c: &Client, name: &str, args: Value) {
    match c
        .call_tool(
            CallToolRequestParams::new(name.to_owned())
                .with_arguments(args.as_object().unwrap().clone()),
        )
        .await
    {
        Ok(r) => assert_eq!(r.is_error, Some(true)),
        Err(rmcp::service::ServiceError::McpError(e)) => assert_eq!(e.code.0, -32602),
        Err(e) => panic!("unexpected transport failure: {e}"),
    }
}

#[tokio::test]
async fn antibody_tools_keep_explicit_limits_exact_sources_full_artifacts_and_readonly_state() {
    use collagene_core::antibody_numbering::{
        self as n, domain_candidates as d, domain_tracks as t,
    };
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
    let tools = c.list_all_tools().await.unwrap();
    assert_eq!(tools.len(), collagene_mcp::CURRENT_TOOL_COUNT);
    for name in [
        "list_antibody_profiles",
        "analyze_antibody",
        "analyze_alignment_antibody",
    ] {
        assert!(tools.iter().any(|t| t.name == name));
    }
    let discovery = call(&c, "list_antibody_profiles", json!({})).await;
    assert_eq!(discovery["v_references"], 354);
    assert_eq!(discovery["j_references"], 23);
    let seq = "EVQLLESGGGLVQPGGSLRLSCAASGFTFSSYAMSWVRQAPGKGLEWVSAISGSGGSTYYADSVKGRFTISRDNSKNTLYLQMNSLRAEDTAVYYCARDYGMDVWGQGTLVTVSS";
    let imported = call(
        &c,
        "import_compat_sequence",
        json!({"text":format!(">source\n{seq}\n"),"profile":{"Fasta":{"alphabet":"Protein"}}}),
    )
    .await;
    let doc = json!({"id":imported["documents"][0]["id"],"revision":imported["documents"][0]["revision"]});
    let before = call(&c, "list_documents", json!({})).await;
    let selection = t::Selection {
        numbering_scheme: t::NumberingScheme::Imgt,
        region_definition: t::RegionDefinition::Imgt,
    };
    let options = d::Options {
        v_scoring: n::v_mapping::Scoring {
            observed_reward: 5,
            frequency_bonus: 8,
            unseen_penalty: -6,
            gap_open: -12,
            gap_extend: -2,
            minimum_score: 300,
        },
        min_cdr3_residues: 2,
        max_cdr3_residues: 60,
        max_fr4_mismatches: 1,
    };
    let mut limits = d::Limits {
        v_mapping: n::v_mapping::Limits {
            max_cells: n::profile::MAX_CELLS,
            max_endpoints_per_locus: n::profile::MAX_ENDPOINTS,
            max_trace_states: n::profile::MAX_TRACE_STATES,
        },
        max_scan_windows: d::MAX_SCAN_WINDOWS,
        max_candidates: 4,
    };
    let r = call(
        &c,
        "analyze_antibody",
        json!({"document":doc,"options":options,"limits":limits,"selection":selection}),
    )
    .await;
    assert_eq!(r["retained_candidates"], 1);
    assert_eq!(r["source_revision"], doc["revision"]);
    let mut text = String::new();
    let mut start = 0;
    loop {
        let chunk = call(
            &c,
            "read_artifact",
            json!({"artifact_id":r["artifact_id"],"start":start,"max_bytes":16384}),
        )
        .await;
        text.push_str(chunk["text"].as_str().unwrap());
        if chunk["done"] == true {
            break;
        }
        start = chunk["next_start"].as_u64().unwrap();
        assert!(text.len() < 4 * 1024 * 1024);
    }
    let report: t::Report = serde_json::from_str(&text).unwrap();
    report
        .validate_binding(r["context_revision"].as_str().unwrap())
        .unwrap();
    assert_eq!(
        report.tracks[0]
            .residues
            .iter()
            .map(|r| r.residue)
            .collect::<String>(),
        seq
    );
    limits.max_candidates = 0;
    let counts = call(
        &c,
        "analyze_antibody",
        json!({"document":doc,"options":options,"limits":limits,"selection":selection}),
    )
    .await;
    assert_eq!(counts["retained_candidates"], 0);
    assert_eq!(counts["pair_retention_truncated"], true);
    assert_eq!(
        counts["qualified_pairs_in_retained_v"],
        r["qualified_pairs_in_retained_v"]
    );
    let mut stale = doc.clone();
    stale["revision"] = json!("f".repeat(64));
    reject(
        &c,
        "analyze_antibody",
        json!({"document":stale,"options":options,"limits":limits,"selection":selection}),
    )
    .await;
    limits.max_candidates = 65;
    reject(
        &c,
        "analyze_antibody",
        json!({"document":doc,"options":options,"limits":limits,"selection":selection}),
    )
    .await;
    reject(
        &c,
        "analyze_antibody",
        json!({"document":doc,"options":options}),
    )
    .await;
    assert_eq!(call(&c, "list_documents", json!({})).await, before);
    limits.max_candidates = 4;
    let mut alternate = selection.clone();
    alternate.numbering_scheme = t::NumberingScheme::Kabat1991;
    let mixed = call(
        &c,
        "analyze_antibody",
        json!({"document":doc,"options":options,"limits":limits,"selection":alternate}),
    )
    .await;
    assert_eq!(mixed["selection"]["numbering_scheme"], "Kabat1991");
    assert_eq!(mixed["candidates"][0]["unassigned_residues"], 8);
    assert_eq!(mixed["candidates"][0]["unassigned_regions"], 0);
    reject(
        &c,
        "analyze_antibody",
        json!({"document":doc,"options":options,"limits":limits}),
    )
    .await;
    // Analyze alignment evidence directly; do not forge an extracted document ref.
    use base64::{Engine as _, engine::general_purpose::STANDARD};
    let input = format!(">row-a\n--{seq}--\n>row-b\n--{seq}--\n");
    let imported=call(&c,"import_alignment",json!({"data_base64":STANDARD.encode(input.as_bytes()),"format":"AlignedFasta","alphabet":"Protein","name":"Antibody row fixture"})).await;
    limits.max_candidates = 4;
    let alignment_before=call(&c,"read_alignment",json!({"alignment":imported["alignment"],"row_start":0,"max_rows":2,"column_start":0,"max_columns":200})).await;
    let row = imported["rows"][0]["id"].clone();
    let ar = call(
        &c,
        "analyze_alignment_antibody",
        json!({"alignment":imported["alignment"],"row_id":row,"options":options,"limits":limits,"selection":selection}),
    )
    .await;
    assert_eq!(ar["retained_candidates"], 1);
    assert_eq!(
        ar["source"]["alignment_revision"],
        imported["alignment"]["revision"]
    );
    let mut text = String::new();
    let mut start = 0;
    loop {
        let chunk = call(
            &c,
            "read_artifact",
            json!({"artifact_id":ar["artifact_id"],"start":start,"max_bytes":16384}),
        )
        .await;
        text.push_str(chunk["text"].as_str().unwrap());
        if chunk["done"] == true {
            break;
        }
        start = chunk["next_start"].as_u64().unwrap();
        assert!(text.len() < 4 * 1024 * 1024);
    }
    let full: collagene_core::alignment_antibody::Report = serde_json::from_str(&text).unwrap();
    full.validate_binding(ar["context_revision"].as_str().unwrap())
        .unwrap();
    assert!(full.at_column(0, 0).is_none());
    assert_eq!(full.at_column(0, 2).unwrap().source_index, 0);
    let other=call(&c,"analyze_alignment_antibody",json!({"alignment":imported["alignment"],"row_id":imported["rows"][1]["id"],"options":options,"limits":limits,"selection":selection})).await;
    assert_ne!(ar["context_revision"], other["context_revision"]);
    reject(&c,"analyze_alignment_antibody",json!({"alignment":imported["alignment"],"row_id":"unknown","options":options,"limits":limits,"selection":selection})).await;
    let mut stale = imported["alignment"].clone();
    stale["revision"] = json!("f".repeat(64));
    reject(
        &c,
        "analyze_alignment_antibody",
        json!({"alignment":stale,"row_id":row,"options":options,"limits":limits,"selection":selection}),
    )
    .await;
    assert_eq!(call(&c,"read_alignment",json!({"alignment":imported["alignment"],"row_start":0,"max_rows":2,"column_start":0,"max_columns":200})).await,alignment_before);
    assert_eq!(call(&c, "list_documents", json!({})).await, before);
    c.cancel().await.unwrap();
    server.await.unwrap();
}
