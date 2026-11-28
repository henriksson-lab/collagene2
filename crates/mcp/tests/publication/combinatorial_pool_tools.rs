//! Real SDK handoff from source-bound reaction discovery to explicit native release.
use base64::{Engine, engine::general_purpose::STANDARD};
use collagene_core::{
    cloning::golden_gate::GoldenGateOptions,
    documents::workspace::{self, Workspace},
};
use collagene_mcp::server::CollageneServer;
use collagene_model::{Molecule, Topology};
use rmcp::{RoleClient, ServiceExt, model::CallToolRequestParams, service::RunningService};
use serde_json::{Value, json};
use std::io::{Cursor, Read};
type Client = RunningService<RoleClient, ()>;
const A: &str = "ATGACCCCGGGGTTTCCAGA";
const B: &str = "ATGAGTTAGCGACCTTGGAA";
const AB: &str = "ATGACCCCGGGGTTTCCAGAATGAGTTAGCGACCTTGGAA";
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
async fn artifact(c: &Client, id: &Value) -> (String, String) {
    let (mut start, mut text) = (0, String::new());
    loop {
        let page = call(
            c,
            "read_artifact",
            json!({"artifact_id":id,"start":start,"max_bytes":16384}),
        )
        .await;
        text.push_str(page["text"].as_str().unwrap());
        let encoding = page["encoding"].as_str().unwrap_or("utf8").to_owned();
        if page["done"] == true {
            return (text, encoding);
        }
        start = page["next_start"].as_u64().unwrap();
    }
}
fn reference(v: &Value) -> Value {
    json!({"id":v["workspace"]["id"],"revision":v["workspace"]["revision"]})
}
fn fixture() -> (Workspace, Value, Value) {
    let w = Workspace::from_documents(
        "reaction SDK",
        vec![
            Molecule::dna(
                "whole vector",
                &format!("ATGCGGTCTCA{A}ATGATGAGACCGCAT"),
                Topology::Circular,
            )
            .unwrap(),
            Molecule::dna(
                "whole amplicon",
                &format!("ATGCGGTCTCA{B}ATGATGAGACCGCAT"),
                Topology::Linear,
            )
            .unwrap(),
            Molecule::dna("independent AB target", AB, Topology::Circular).unwrap(),
        ],
        0,
    )
    .unwrap();
    let bindings = w.documents.iter().map(|d| json!({"document_id":d.id,"event_id":d.head_event.as_ref().unwrap(),"revision":d.head})).collect::<Vec<_>>();
    let bins = (0..2).map(|i| json!({"id":format!("bin-{i}"),"label":format!("Bin {i}"),"choices":[{"generation":{"dam":"Unmodified","dcm":"Unmodified","cpg":"Unmodified"},
        "id":format!("choice-{i}"),"label":"Whole reaction material","source":bindings[i],
        "selector":{"ReferenceRange":{"start":0,"end":w.revisions[&w.documents[i].head].sequence.len(),"reverse":false}},"preparation_reverse":false,"preparation":null,"treatment":null
    }]})).collect::<Vec<_>>();
    let design = json!({"id":"reaction-sdk","name":"Whole reaction alternatives","circular":true,"bins":bins,"material_bindings":[],
        "method":{"GoldenGatePool":{"options":GoldenGateOptions {max_fragments:3,..Default::default()},
        "substrates":[{"bin_id":"bin-0","kind":"WholeSource"},{"bin_id":"bin-1","kind":"WholeSource"}],"selection":null}}});
    let row = json!({"id":"construct-one","name":"Explicit Ω path","choices":[{"bin_id":"bin-0","choice_id":"choice-0"},{"bin_id":"bin-1","choice_id":"choice-1"}],"target":{"Existing":bindings[2]},"overrides":null});
    (w, design, row)
}
#[tokio::test]
async fn discovery_explicit_keys_batch_release_and_portable_protocol_keep_one_path() {
    let (a, b) = tokio::io::duplex(65_536);
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
    let tool = tools
        .iter()
        .find(|tool| tool.name == "inspect_combinatorial_pool")
        .unwrap();
    assert_eq!(tool.input_schema["additionalProperties"], false);
    let (w, mut design, row) = fixture();
    let imported = call(
        &c,
        "import_workspace_archive",
        json!({"data_base64":STANDARD.encode(workspace::save_workspace(&w).unwrap())}),
    )
    .await;
    let source = reference(&imported);
    let discovered = call(
        &c,
        "inspect_combinatorial_pool",
        json!({"workspace":source,"design":design,"construct":row}),
    )
    .await;
    assert_eq!(discovered["state"], "fresh-pool-discovery");
    assert_eq!(discovered["verified"], false);
    assert_eq!(discovered["input_bindings"].as_array().unwrap().len(), 2);
    let (text, _) = artifact(&c, &discovered["artifact_id"]).await;
    let evidence: Value = serde_json::from_str(&text).unwrap();
    assert_eq!(evidence["source_revision"], discovered["source_revision"]);
    assert_eq!(
        evidence["report"]["context_revision"],
        discovered["context_revision"]
    );
    assert_eq!(
        evidence["report"]["pool"]["fragments"]
            .as_array()
            .unwrap()
            .len(),
        5,
        "whole circular and linear sources retain every cut fragment"
    );
    let candidate = evidence["report"]["pool"]["candidates"]
        .as_array()
        .unwrap()
        .iter()
        .enumerate()
        .find(|(_, candidate)| {
            candidate["duplex"]["top"]["bases"] == AB && candidate["restriction_site_free"] == true
        })
        .unwrap();
    let selection = json!({"selection":{"context_revision":discovered["context_revision"],"candidate_index":candidate.0,
        "directed_sequence_key":candidate.1["directed_sequence_key"],"physical_sequence_key":candidate.1["physical_sequence_key"]},"source_revision":discovered["source_revision"]});
    design["method"]["GoldenGatePool"]["selection"] = selection.clone();
    let evaluated = call(
        &c,
        "evaluate_combinatorial_batch",
        json!({"workspace":source,"matrix":{"design":design,"rows":[{"Ready":row}]}}),
    )
    .await;
    assert_eq!(
        evaluated["receipt"]["rows"][0]["status"],
        "ConsistentUnderAssumptions"
    );
    let released = call(&c,"release_combinatorial_outputs",json!({"workspace":reference(&evaluated),"request":{"batch_hash":evaluated["receipt"]["attachment_hash"],"rows":[{"construct_id":"construct-one","options":{"name":"Released Ω","show_prepared":false,"save_primers":false,"allow_incomplete":false}}]}})).await;
    let (native, encoding) = artifact(&c, &released["native_artifact_id"]).await;
    assert_eq!(encoding, "base64");
    let retained = workspace::load_workspace(&STANDARD.decode(native).unwrap()).unwrap();
    let runs = collagene_core::material::lineage::list_runs(&retained).unwrap();
    assert_eq!(runs.len(), 1);
    let run = collagene_core::material::lineage::load_run(&retained, &runs[0].run_id).unwrap();
    let op = serde_json::to_value(run.plan.steps.last().unwrap()).unwrap();
    assert_eq!(
        op["SelectedGoldenGate"]["directed_sequence_key"],
        selection["selection"]["directed_sequence_key"]
    );
    assert_eq!(
        op["SelectedGoldenGate"]["physical_sequence_key"],
        selection["selection"]["physical_sequence_key"]
    );
    let exported = call(&c,"export_combinatorial_outputs",json!({"workspace":reference(&released),"manifest_hash":released["receipt"]["manifest_hash"],"options":{"genbank":false,"fasta":true,"primer_csv":false,"primer_tsv":false,"protocol_html":true,"protocol_pdf":false,"construct_ids":null}})).await;
    let (zip, _) = artifact(&c, &exported["artifact_id"]).await;
    let mut zip = zip::ZipArchive::new(Cursor::new(STANDARD.decode(zip).unwrap())).unwrap();
    let mut html_seen = false;
    for i in 0..zip.len() {
        let mut file = zip.by_index(i).unwrap();
        if file.name().ends_with(".html") {
            let mut html = String::new();
            file.read_to_string(&mut html).unwrap();
            assert!(
                html.contains(
                    selection["selection"]["directed_sequence_key"]
                        .as_str()
                        .unwrap()
                )
            );
            assert!(
                html.contains(
                    selection["selection"]["physical_sequence_key"]
                        .as_str()
                        .unwrap()
                )
            );
            html_seen = true;
        }
    }
    assert!(html_seen);
    let mut stale = design.clone();
    stale["method"]["GoldenGatePool"]["selection"]["source_revision"] = json!("0".repeat(64));
    let invalid = call(
        &c,
        "evaluate_combinatorial_batch",
        json!({"workspace":source,"matrix":{"design":stale,"rows":[{"Ready":row}]}}),
    )
    .await;
    assert_ne!(
        invalid["receipt"]["rows"][0]["status"],
        "ConsistentUnderAssumptions"
    );
    reject(&c,"release_combinatorial_outputs",json!({"workspace":reference(&invalid),"request":{"batch_hash":invalid["receipt"]["attachment_hash"],"rows":[{"construct_id":"construct-one","options":{"name":"Never released","show_prepared":false,"save_primers":false,"allow_incomplete":false}}]}})).await;
    drop(c);
    server.abort();
}
