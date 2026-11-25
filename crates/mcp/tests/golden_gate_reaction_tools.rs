use collagene_core::{Molecule, Topology, analysis::Primer, golden_gate::GoldenGateOptions};
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
fn doc(v: &Value) -> Value {
    json!({"id":v["documents"][0]["id"],"revision":v["documents"][0]["revision"]})
}
async fn import(c: &Client, m: Molecule) -> Value {
    doc(&call(
        c,
        "import_document",
        json!({"text":collagene_core::formats::export_genbank(&m).unwrap()}),
    )
    .await)
}
async fn full(c: &Client, id: &Value) -> Value {
    let mut text = String::new();
    let mut start = 0;
    loop {
        let r = call(
            c,
            "read_artifact",
            json!({"artifact_id":id,"start":start,"max_bytes":16384}),
        )
        .await;
        text.push_str(r["text"].as_str().unwrap());
        if r["done"] == true {
            break;
        }
        start = r["next_start"].as_u64().unwrap();
    }
    serde_json::from_str(&text).unwrap()
}
#[tokio::test]
async fn whole_reaction_tools_retain_vector_stuffer_amplicon_flanks_and_uncut_carryover() {
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
        c.list_tools(Default::default()).await.unwrap().tools.len(),
        collagene_mcp::CURRENT_TOOL_COUNT
    );
    // Independent construct precedes any preparation or reaction inspection.
    let target = import(
        &c,
        Molecule::dna(
            "desired",
            "ATGACCCCGGGGTTTCCAGAATGAGTTAGCGACCTTGGAA",
            Topology::Circular,
        )
        .unwrap(),
    )
    .await;
    let vector = import(
        &c,
        Molecule::dna(
            "whole vector",
            "ATGCGGTCTCAATGACCCCGGGGTTTCCAGAATGATGAGACCGCAT",
            Topology::Circular,
        )
        .unwrap(),
    )
    .await;
    let insert = import(
        &c,
        Molecule::dna("insert template", "ATGAGTTAGCGACCTTGGAA", Topology::Linear).unwrap(),
    )
    .await;
    let carry = import(
        &c,
        Molecule::dna("uncut carryover", "ACGTACGTACGT", Topology::Circular).unwrap(),
    )
    .await;
    let before = call(&c, "export_document", vector.clone()).await;
    let inputs = json!([{"id":"vector","document":vector},{"id":"insert","document":insert},{"id":"carry","document":carry}]);
    let request = json!({"name":"declared tube","rows":[{"role":"vector","input_id":"vector","substrate":"WholeSource"},{"role":"amplicon","input_id":"insert","substrate":{"Prepared":{"preparation":{"ExactPcr":{"forward":Primer::new("manualF","ATGAGTTA","ATGCGGTCTCA").unwrap(),"reverse":Primer::new("manualR","TTCCAAGG","ATGCGGTCTCATCAT").unwrap()}},"treatment":null}}},{"role":"uncut","input_id":"carry","substrate":"WholeSource"}],"options":GoldenGateOptions{max_fragments:3,..Default::default()}});
    let args = json!({"inputs":inputs,"target":target,"request":request});
    let inspected = call(&c, "inspect_golden_gate_reaction", args.clone()).await;
    assert_eq!(inspected["state"], "advisory");
    assert_eq!(inspected["verified"], false);
    assert_eq!(inspected["carryover_sources"], json!([2]));
    assert_eq!(inspected["pool_available"], true);
    let evidence = full(&c, &inspected["artifact_id"]).await;
    assert_eq!(evidence["kind"], "golden-gate-reaction-v1");
    assert_eq!(
        evidence["report"]["pool"]["fragments"]
            .as_array()
            .unwrap()
            .len(),
        6
    );
    assert_eq!(
        evidence["report"]["probe_plan"]["steps"]
            .as_array()
            .unwrap()
            .len(),
        2
    );
    let a = inspected["candidates"]
        .as_array()
        .unwrap()
        .iter()
        .find(|a| a["target_matches"] == true && a["restriction_site_free"] == true)
        .unwrap();
    let selection = json!({"context_revision":inspected["context_revision"],"candidate_index":a["candidate_index"],"directed_sequence_key":a["directed_sequence_key"],"physical_sequence_key":a["physical_sequence_key"]});
    let selected = call(
        &c,
        "select_golden_gate_reaction_candidate",
        json!({"artifact_id":inspected["artifact_id"],"selection":selection}),
    )
    .await;
    assert_eq!(selected["verified"], false);
    assert_eq!(selected["artifact_id"], selected["candidate_artifact_id"]);
    assert_eq!(
        call(
            &c,
            "verify_plan",
            json!({"plan_revision":selected["plan_revision"]})
        )
        .await["status"],
        "ConsistentUnderAssumptions"
    );
    assert_eq!(
        call(&c, "inspect_golden_gate_reaction", args.clone()).await["artifact_id"],
        inspected["artifact_id"]
    );
    let mut capped = args;
    capped["request"]["options"]["max_states"] = 1.into();
    capped["request"]["options"]["min_fragments"] = 2.into();
    let partial = call(&c, "inspect_golden_gate_reaction", capped).await;
    assert_eq!(partial["enumeration_complete"], false);
    assert_eq!(partial["pool_available"], true);
    assert_eq!(partial["candidate_count"], 0);
    let mut bad = selection;
    bad["physical_sequence_key"] = "0".repeat(64).into();
    let result = c
        .call_tool(
            CallToolRequestParams::new("select_golden_gate_reaction_candidate").with_arguments(
                json!({"artifact_id":inspected["artifact_id"],"selection":bad})
                    .as_object()
                    .unwrap()
                    .clone(),
            ),
        )
        .await;
    assert!(result.is_err() || result.unwrap().is_error == Some(true));
    assert_eq!(call(&c, "export_document", vector).await, before);
    c.cancel().await.unwrap();
    server.await.unwrap();
}
