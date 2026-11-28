//! Actual prepared strand pools retain backgrounds and require fresh selected-plan verification.
use base64::{Engine, engine::general_purpose::STANDARD};
use collagene_core::{documents::transfer, formats::sequence as formats};
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
async fn evidence(c: &Client, id: &Value) -> Value {
    let mut text = String::new();
    let mut start = 0;
    loop {
        let page = call(
            c,
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
fn named(imported: &Value, name: &str) -> Value {
    let row = imported["documents"]
        .as_array()
        .unwrap()
        .iter()
        .find(|d| d["name"] == name)
        .unwrap();
    json!({"id":row["id"],"revision":row["revision"]})
}
const A: &str = "ATGACCCCGGGGTTTCCAGA";
const B: &str = "ATGAGTTAGCGACCTTGGAA";
const TARGET: &str = "ATGACCCCGGGGTTTCCAGAATGAGTTAGCGACCTTGGAA";
fn options() -> Value {
    json!({"enzyme":"BsaI","circular":true,"both_orientations":true,"min_fragments":1,"max_fragments":2,"max_products":128,"max_states":50000})
}
fn wr_candidate(pool: &Value, candidate: &Value) -> Value {
    json!({"context_revision":pool["context_revision"],"candidate_index":candidate["candidate_index"],"directed_sequence_key":candidate["directed_sequence_key"],"physical_sequence_key":candidate["physical_sequence_key"]})
}
async fn bytes(c: &Client, id: &Value) -> Vec<u8> {
    let mut encoded = String::new();
    let mut start = 0;
    loop {
        let page = call(
            c,
            "read_artifact",
            json!({"artifact_id":id,"start":start,"max_bytes":16384}),
        )
        .await;
        assert_eq!(page["encoding"], "base64");
        encoded.push_str(page["text"].as_str().unwrap());
        if page["done"] == true {
            break;
        }
        start = page["next_start"].as_u64().unwrap();
    }
    STANDARD.decode(encoded).unwrap()
}
async fn fixture(c: &Client) -> (Value, Value, Vec<Value>, Value) {
    let originals = call(
        c,
        "import_document",
        json!({"text":format!(">A\n{A}\n>B\n{B}\n")}),
    )
    .await;
    let a = named(&originals, "A");
    let b = named(&originals, "B");
    // Independent circular40nt target is imported before any primer design or enumeration.
    let gb = formats::export_genbank(
        &Molecule::dna("Independent_target", TARGET, Topology::Circular).unwrap(),
    )
    .unwrap();
    let imported = call(c, "import_document", json!({"text":gb})).await;
    let target = named(&imported, "Independent_target");
    let inputs = json!([{"id":"a","document":a},{"id":"b","document":b}]);
    let request = json!({"name":"Two explicit same-fusion substrates","fragments":[{"role":"A-role","input_id":"a","start":0,"end":20,"reverse":false,"forward_length":8,"reverse_length":8,"right_overhang":"ATGA"},{"role":"B-role","input_id":"b","start":0,"end":20,"reverse":false,"forward_length":8,"reverse_length":8,"right_overhang":"ATGA"}],"enzyme":"BsaI","clamp":"ATGC","spacer":"A"});
    let proposal = call(
        c,
        "design_golden_gate_primers",
        json!({"inputs":inputs,"target":target,"request":request}),
    )
    .await;
    let recipe = json!({"name":"Ambiguous same-fusion physical pool","fragments":proposal["fragments"].as_array().unwrap().iter().map(|source|json!({"source":source,"treatment":null})).collect::<Vec<_>>(),"circular":true,"allow_nicked":false});
    (
        inputs,
        target,
        vec![a, b],
        json!({"recipe":recipe,"options":options()}),
    )
}
#[tokio::test]
async fn prepared_pool_keeps_self_backgrounds_and_selected_plans_reconstruct_against_independent_target()
 {
    let (c, server) = setup().await;
    assert_eq!(
        c.list_tools(Default::default()).await.unwrap().tools.len(),
        collagene_mcp::CURRENT_TOOL_COUNT
    );
    let (inputs, target, docs, request) = fixture(&c).await;
    assert_eq!(
        call(
            &c,
            "read_sequence_region",
            json!({"document":target,"start":0,"end":40})
        )
        .await["sequence"],
        TARGET
    );
    let mut originals = Vec::new();
    for doc in &docs {
        originals.push(call(&c, "export_document", doc.clone()).await);
    }
    let args = json!({"origin_workspaces":[],"inputs":inputs,"target":target,"request":request});
    let pool = call(&c, "inspect_prepared_golden_gate_pool", args.clone()).await;
    assert_eq!(pool["state"], "advisory");
    assert_eq!(pool["verified"], false);
    assert_eq!(pool["preparation_complete"], true);
    assert_eq!(pool["missing_physical"], json!([]));
    assert_eq!(pool["enumeration_complete"], true);
    assert_eq!(pool["exhaustive_mixture_search"], false);
    let candidates = pool["candidates"].as_array().unwrap();
    let groups = candidates
        .iter()
        .map(|p| p["physical_sequence_key"].as_str().unwrap())
        .collect::<std::collections::BTreeSet<_>>();
    assert!(groups.len() >= 3);
    let actual = candidates
        .iter()
        .find(|p| p["target_matches"] == true && p["length"] == 40)
        .unwrap();
    let self_a = candidates
        .iter()
        .find(|p| p["length"] == 20 && p["source_roles"] == json!([0]))
        .unwrap();
    assert_eq!(self_a["target_matches"], false);
    assert!(
        candidates
            .iter()
            .any(|p| p["length"] == 20 && p["source_roles"] == json!([1]))
    );
    let full = evidence(&c, &pool["artifact_id"]).await;
    assert_eq!(full["kind"], "prepared-golden-gate-pool-v1");
    assert_eq!(full["inputs"], inputs);
    assert_eq!(full["target"], target);
    assert_eq!(full["report"]["context_revision"], pool["context_revision"]);
    assert_eq!(full["report"]["request"], request);
    assert_eq!(
        full["report"]["pool"]["enumeration_complete"],
        pool["enumeration_complete"]
    );
    assert_eq!(
        full["report"]["pool"]["states_examined"],
        pool["states_examined"]
    );
    let original_plan = &full["report"]["preview"]["compiled"]["plan"];
    let original_steps = original_plan["steps"].as_array().unwrap();
    assert!(original_steps.last().unwrap().get("Ligate").is_some());
    let repeat = call(&c, "inspect_prepared_golden_gate_pool", args).await;
    assert_eq!(repeat["artifact_id"], pool["artifact_id"]);
    for (item, status, length) in [
        (actual, "ConsistentUnderAssumptions", 40),
        (self_a, "Invalid", 20),
    ] {
        let selection = wr_candidate(&pool, item);
        let selected = call(
            &c,
            "select_prepared_golden_gate_candidate",
            json!({"artifact_id":pool["artifact_id"],"selection":selection}),
        )
        .await;
        assert_eq!(selected["state"], "candidate");
        assert_eq!(selected["verified"], false);
        assert_eq!(selected["parent_context"], pool["context_revision"]);
        assert_eq!(selected["candidate_artifact_id"], selected["artifact_id"]);
        let candidate =
            transfer::load_candidate(&bytes(&c, &selected["candidate_artifact_id"]).await).unwrap();
        let plan = serde_json::to_value(&candidate.plan).unwrap();
        let steps = plan["steps"].as_array().unwrap();
        assert_eq!(plan["inputs"], original_plan["inputs"]);
        assert_eq!(plan["target_revision"], target["revision"]);
        assert_eq!(
            &steps[..steps.len() - 1],
            &original_steps[..original_steps.len() - 1]
        );
        assert!(steps.last().unwrap().get("GoldenGate").is_some());
        let checked = call(
            &c,
            "verify_plan",
            json!({"plan_revision":selected["plan_revision"]}),
        )
        .await;
        assert_eq!(checked["status"], status);
        assert_eq!(checked["target_revision"], target["revision"]);
        assert_eq!(
            checked["steps"].as_array().unwrap().last().unwrap()["product"]["length"],
            length
        );
        assert_eq!(checked["exhaustive_alternative_search"], false);
        let again = call(
            &c,
            "select_prepared_golden_gate_candidate",
            json!({"artifact_id":pool["artifact_id"],"selection":selection}),
        )
        .await;
        assert_eq!(again["plan_revision"], selected["plan_revision"]);
        assert_eq!(again["artifact_id"], selected["artifact_id"]);
    }
    for (doc, original) in docs.iter().zip(originals) {
        assert_eq!(call(&c, "export_document", doc.clone()).await, original);
    }
    c.cancel().await.unwrap();
    server.await.unwrap();
}
#[tokio::test]
async fn finite_caps_missing_physical_and_stale_candidate_keys_remain_explicit() {
    let (c, server) = setup().await;
    let (inputs, target, _, request) = fixture(&c).await;
    let baseline = call(
        &c,
        "inspect_prepared_golden_gate_pool",
        json!({"origin_workspaces":[],"inputs":inputs,"target":target,"request":request}),
    )
    .await;
    let candidate = &baseline["candidates"][0];
    let good = wr_candidate(&baseline, candidate);
    let mut stale = good.clone();
    stale["context_revision"] = json!("0".repeat(64));
    let mut wrong_directed = good.clone();
    wrong_directed["directed_sequence_key"] = json!("0".repeat(64));
    let mut wrong_physical = good.clone();
    wrong_physical["physical_sequence_key"] = json!("0".repeat(64));
    for selection in [stale, wrong_directed, wrong_physical] {
        reject(
            &c,
            "select_prepared_golden_gate_candidate",
            json!({"artifact_id":baseline["artifact_id"],"selection":selection}),
        )
        .await;
    }
    reject(
        &c,
        "select_prepared_golden_gate_candidate",
        json!({"artifact_id":"0".repeat(64),"selection":good}),
    )
    .await;
    let mut stale_inputs = inputs.clone();
    stale_inputs[1]["document"]["revision"] = json!("0".repeat(64));
    reject(
        &c,
        "inspect_prepared_golden_gate_pool",
        json!({"origin_workspaces":[],"inputs":stale_inputs,"target":target,"request":request}),
    )
    .await;
    let mut capped = request.clone();
    capped["options"]["max_states"] = json!(1);
    let partial = call(
        &c,
        "inspect_prepared_golden_gate_pool",
        json!({"origin_workspaces":[],"inputs":inputs,"target":target,"request":capped}),
    )
    .await;
    assert_eq!(partial["preparation_complete"], true);
    assert_eq!(partial["enumeration_complete"], false);
    assert_eq!(partial["exhaustive_mixture_search"], false);
    assert!(!partial["truncation_reasons"].as_array().unwrap().is_empty());
    assert!(partial["states_examined"].as_u64().unwrap() <= 1);
    let mut failed = request.clone();
    failed["recipe"]["fragments"][1]["source"]["preparation"]["PcrCatalogDigest"]["forward"]["annealing"]
        ["bases"] = json!("AAAAAAAA");
    let withheld = call(
        &c,
        "inspect_prepared_golden_gate_pool",
        json!({"origin_workspaces":[],"inputs":inputs,"target":target,"request":failed}),
    )
    .await;
    assert_eq!(withheld["preparation_complete"], false);
    assert!(!withheld["missing_physical"].as_array().unwrap().is_empty());
    assert_eq!(withheld["candidates"], json!([]));
    let full = evidence(&c, &withheld["artifact_id"]).await;
    assert!(full["report"]["pool"].is_null());
    c.cancel().await.unwrap();
    server.await.unwrap();
}
