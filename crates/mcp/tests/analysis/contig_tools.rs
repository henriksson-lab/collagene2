//! Actual SDK routing and immutable source/report/consensus installation.
use collagene_core::review::contigs::*;
use collagene_mcp::server::CollageneServer;
use rmcp::{RoleClient, ServiceExt, model::CallToolRequestParams, service::RunningService};
use serde_json::{Value, json};
type Client = RunningService<RoleClient, ()>;
async fn call(c: &Client, name: &str, args: Value) -> Value {
    let result = c
        .call_tool(
            CallToolRequestParams::new(name.to_owned())
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
            CallToolRequestParams::new(name.to_owned())
                .with_arguments(args.as_object().unwrap().clone()),
        )
        .await;
    assert!(result.is_err() || result.unwrap().is_error == Some(true));
}
async fn artifact(c: &Client, id: &Value) -> Value {
    let mut text = String::new();
    let mut at = 0;
    loop {
        let page = call(
            c,
            "read_artifact",
            json!({"artifact_id":id,"start":at,"max_bytes":16384}),
        )
        .await;
        text.push_str(page["text"].as_str().unwrap());
        if page["done"] == true {
            break;
        }
        at = page["next_start"].as_u64().unwrap();
    }
    serde_json::from_str(&text).unwrap()
}
fn doc(imported: &Value, name: &str) -> Value {
    let d = imported["documents"]
        .as_array()
        .unwrap()
        .iter()
        .find(|d| {
            d["name"]
                .as_str()
                .is_some_and(|s| s.split_whitespace().next() == Some(name))
        })
        .unwrap();
    json!({"id":d["id"],"revision":d["revision"]})
}
const READS: &str = ">left\nAGTCGACTA\n>middle\nACTAGGCTT\n>right\nGCTTACG\n";
fn options() -> ContigOptions {
    let mut o = ContigOptions::default();
    o.trim.min_retained_bases = 1;
    o.min_overlap = 4;
    o.max_error_bp = 0;
    o.both_strands = false;
    o
}
fn refs(imported: &Value) -> Vec<Value> {
    ["left", "middle", "right"]
        .into_iter()
        .map(|id| json!({"id":id,"document":doc(imported,id)}))
        .collect()
}
#[tokio::test]
async fn assembly_sdk_preserves_full_evidence_source_binding_pages_and_uncertain_consensus() {
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
    let imported = call(&c, "import_document", json!({"text":READS})).await;
    let source_refs = refs(&imported);
    let opts = options();
    let first = call(
        &c,
        "assemble_contigs",
        json!({"reads":source_refs,"options":opts,"max_contigs":1,"max_overlaps":1}),
    )
    .await;
    assert_eq!(first["reused_result"], false);
    assert_eq!(first["read_count"], 3);
    assert_eq!(first["contig_count"], 1);
    assert_eq!(first["contigs"][0]["sequence_preview"], "AGTCGACTAGGCTTACG");
    assert_eq!(first["contigs"][0]["placement_count"], 3);
    let full: ContigReport =
        serde_json::from_value(artifact(&c, &first["artifact_id"]).await).unwrap();
    let inputs = collagene_core::formats::sequence::import_text(READS)
        .unwrap()
        .molecules
        .into_iter()
        .map(|molecule| ReadInput {
            id: molecule.name.clone(),
            molecule,
        })
        .collect::<Vec<_>>();
    full.validate_bindings_for(&inputs, &opts).unwrap();
    assert_eq!(
        first["batch_revision"],
        batch_revision(&inputs, &opts).unwrap()
    );
    let expected = consensus_documents(&inputs, &opts, &full).unwrap();
    assert!(
        expected
            .iter()
            .any(|m| m.sequence.as_str() == "AGTCGACTAGGCTTACG")
    );
    let consensus = first["consensus_documents"]
        .as_array()
        .unwrap()
        .iter()
        .find(|d| d["length"] == 17)
        .unwrap();
    let consensus_ref = json!({"id":consensus["id"],"revision":consensus["revision"]});
    let sequence = call(
        &c,
        "read_sequence_region",
        json!({"document":consensus_ref,"start":0,"end":17}),
    )
    .await;
    assert_eq!(sequence["sequence"], "AGTCGACTAGGCTTACG");
    let page = call(
        &c,
        "assemble_contigs",
        json!({"reads":source_refs,"options":opts,"overlap_start":1,"max_overlaps":1}),
    )
    .await;
    assert_eq!(page["reused_result"], true);
    assert_eq!(page["artifact_id"], first["artifact_id"]);
    assert_eq!(page["consensus_documents"], first["consensus_documents"]);
    let mut stale = source_refs.clone();
    stale[0]["document"]["revision"] = json!("0".repeat(64));
    reject(
        &c,
        "assemble_contigs",
        json!({"reads":stale,"options":opts}),
    )
    .await;
    let mut duplicate_role = source_refs.clone();
    duplicate_role[1]["id"] = duplicate_role[0]["id"].clone();
    reject(
        &c,
        "assemble_contigs",
        json!({"reads":duplicate_role,"options":opts}),
    )
    .await;
    reject(
        &c,
        "assemble_contigs",
        json!({"reads":vec![source_refs[0].clone();33],"options":opts}),
    )
    .await;
    let mut unknown = serde_json::to_value(&opts).unwrap();
    unknown["trust_ai_consensus"] = json!(true);
    reject(
        &c,
        "assemble_contigs",
        json!({"reads":source_refs,"options":unknown}),
    )
    .await;
    let mut partial = opts.clone();
    partial.max_total_cells = 1;
    let bounded = call(
        &c,
        "assemble_contigs",
        json!({"reads":source_refs,"options":partial}),
    )
    .await;
    assert_eq!(bounded["complete"], false);
    assert!(!bounded["limits"].as_array().unwrap().is_empty());
    let mut support = opts.clone();
    support.min_support = 2;
    let uncertain = call(
        &c,
        "assemble_contigs",
        json!({"reads":source_refs,"options":support}),
    )
    .await;
    let report: ContigReport =
        serde_json::from_value(artifact(&c, &uncertain["artifact_id"]).await).unwrap();
    assert!(report.contigs.iter().any(|x| x.sequence.contains('N')));
    assert!(
        report
            .contigs
            .iter()
            .flat_map(|x| &x.columns)
            .any(|x| x.state == ConsensusState::InsufficientSupport)
    );
    let conflicting = call(
        &c,
        "import_document",
        json!({"text":">conflict_left\nGGACGATTCG\n>conflict_right\nGACTCGTTA\n"}),
    )
    .await;
    let conflict_inputs = vec![
        json!({"id":"left","document":doc(&conflicting,"conflict_left")}),
        json!({"id":"right","document":doc(&conflicting,"conflict_right")}),
    ];
    let mut conflict_options = opts.clone();
    conflict_options.max_error_bp = 3000;
    let conflict = call(
        &c,
        "assemble_contigs",
        json!({"reads":conflict_inputs,"options":conflict_options}),
    )
    .await;
    let conflict: ContigReport =
        serde_json::from_value(artifact(&c, &conflict["artifact_id"]).await).unwrap();
    assert!(
        conflict
            .contigs
            .iter()
            .flat_map(|c| &c.columns)
            .any(|c| c.state == ConsensusState::ConflictingCalls)
    );
    assert!(conflict.contigs.iter().any(|c| c.sequence.contains('N')));
    let forked = call(
        &c,
        "import_document",
        json!({"text":">fork\nAGTCCGGA\n>arm_a\nCGGATTTG\n>arm_b\nCGGACCCA\n"}),
    )
    .await;
    let fork_inputs = ["fork", "arm_a", "arm_b"]
        .into_iter()
        .map(|id| json!({"id":id,"document":doc(&forked,id)}))
        .collect::<Vec<_>>();
    let fork = call(
        &c,
        "assemble_contigs",
        json!({"reads":fork_inputs,"options":opts}),
    )
    .await;
    let fork: ContigReport =
        serde_json::from_value(artifact(&c, &fork["artifact_id"]).await).unwrap();
    assert_eq!(
        fork.overlaps
            .iter()
            .filter(|e| e.from == 0 && e.accepted && e.kind == OverlapKind::Extension)
            .count(),
        2
    );
    assert_eq!(
        fork.overlaps
            .iter()
            .filter(|e| e.from == 0 && e.selected)
            .count(),
        0
    );
    assert!(fork.contigs.iter().any(|c| !c.layout_unique));
    let repeated = call(
        &c,
        "import_document",
        json!({"text":">repeat_a\nACGTTGCA\n>repeat_b\nACGTTGCA\n"}),
    )
    .await;
    let duplicate_inputs = vec![
        json!({"id":"a","document":doc(&repeated,"repeat_a")}),
        json!({"id":"b","document":doc(&repeated,"repeat_b")}),
    ];
    let duplicate = call(
        &c,
        "assemble_contigs",
        json!({"reads":duplicate_inputs,"options":opts}),
    )
    .await;
    let report: ContigReport =
        serde_json::from_value(artifact(&c, &duplicate["artifact_id"]).await).unwrap();
    assert!(report.reads.iter().any(|r| r.duplicate_group.is_some()));
    c.cancel().await.unwrap();
    server.await.unwrap();
}
#[tokio::test]
async fn full_artifact_capacity_never_installs_orphan_consensus_documents_or_cache() {
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
    let imported = call(
        &c,
        "import_document",
        json!({"text":format!("{READS}>protein [molecule=protein]\nMKW\n")}),
    )
    .await;
    let protein = doc(&imported, "protein");
    for i in 0..128 {
        let options = collagene_core::protein::properties::ProteinOptions {
            charge_ph: i as f64 / 10.,
            ..Default::default()
        };
        call(
            &c,
            "analyze_protein",
            json!({"document":protein,"options":options}),
        )
        .await;
    }
    let before = call(&c, "list_documents", json!({"start":0,"limit":128})).await;
    reject(
        &c,
        "assemble_contigs",
        json!({"reads":refs(&imported),"options":options()}),
    )
    .await;
    let after = call(&c, "list_documents", json!({"start":0,"limit":128})).await;
    assert_eq!(before, after);
    c.cancel().await.unwrap();
    server.await.unwrap();
}
