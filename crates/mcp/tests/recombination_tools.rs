//! SDK boundary tests use independently specified GenBank expected products, never a returned product as target.
use collagene_mcp::CollageneServer;
use rmcp::{
    RoleClient, ServiceExt,
    model::{CallToolRequestParams, ReadResourceRequestParams},
    service::RunningService,
};
use serde_json::{Value, json};
type Client = RunningService<RoleClient, ()>;
async fn setup() -> (Client, tokio::task::JoinHandle<()>) {
    let (a, b) = tokio::io::duplex(65536);
    let s = tokio::spawn(async move {
        CollageneServer::default()
            .serve(a)
            .await
            .unwrap()
            .waiting()
            .await
            .unwrap();
    });
    (().serve(b).await.unwrap(), s)
}
async fn call(c: &Client, n: &str, a: Value) -> Value {
    let r = c
        .call_tool(
            CallToolRequestParams::new(n.to_owned()).with_arguments(a.as_object().unwrap().clone()),
        )
        .await
        .unwrap();
    assert_ne!(r.is_error, Some(true), "{n}: {:?}", r.content);
    r.structured_content.unwrap()
}
async fn reject(c: &Client, n: &str, a: Value) {
    let r = c
        .call_tool(
            CallToolRequestParams::new(n.to_owned()).with_arguments(a.as_object().unwrap().clone()),
        )
        .await;
    assert!(r.is_err() || r.unwrap().is_error == Some(true));
}
fn doc(v: &Value) -> Value {
    json!({"id":v["id"],"revision":v["revision"]})
}
async fn prepare(c: &Client) -> (Value, Value, Value) {
    let imported = call(
        c,
        "import_document",
        json!({"text":include_str!("../../../fixtures/recombination-bp-demo.gb")}),
    )
    .await;
    let docs = imported["documents"].as_array().unwrap();
    assert_eq!(docs.len(), 4);
    let resource = c
        .read_resource(ReadResourceRequestParams::new(
            "collagene://recombination/profile/manufacturer-context-v1",
        ))
        .await
        .unwrap();
    let value = serde_json::to_value(resource).unwrap();
    let profile: Value =
        serde_json::from_str(value["contents"][0]["text"].as_str().unwrap()).unwrap();
    let request = serde_json::to_value(collagene_core::recombination::Request::new(
        serde_json::from_value(profile["snapshot"].clone()).unwrap(),
        collagene_core::recombination::reference::Reaction::BP,
    ))
    .unwrap();
    let named = |name: &str| {
        docs.iter()
            .find(|d| {
                d["name"]
                    .as_str()
                    .is_some_and(|n| n.split_whitespace().next() == Some(name))
            })
            .unwrap()
    };
    assert_eq!(named("BP_insert")["length"], 59);
    assert_eq!(named("BP_donor")["length"], 4761);
    assert_eq!(named("BP_target")["length"], 2553);
    (
        json!([{"role":"insert","document":doc(named("BP_insert"))},{"role":"donor","document":doc(named("BP_donor"))}]),
        request,
        doc(named("BP_target")),
    )
}
#[tokio::test]
async fn explicit_profile_discovery_all_reciprocal_components_and_fresh_checked_selection() {
    let (c, s) = setup().await;
    assert_eq!(
        c.list_tools(Default::default()).await.unwrap().tools.len(),
        97
    );
    let (sources, request, target) = prepare(&c).await;
    let r = call(
        &c,
        "inspect_recombination",
        json!({"sources":sources,"request":request,"limit":1}),
    )
    .await;
    assert_eq!(r["complete"], true);
    assert_eq!(r["candidate_count"], 1);
    let components = r["candidates"][0]["components"].as_array().unwrap();
    assert_eq!(components.len(), 2);
    assert_eq!(
        components
            .iter()
            .map(|p| p["length"].as_u64().unwrap())
            .sum::<u64>(),
        4820
    );
    let check=call(&c,"check_recombination_target",json!({"sources":sources,"request":request,"target":target,"policy":{"allow_reverse_complement":false}})).await;
    assert_eq!(
        check["assessment"]["reference"]["matches"]
            .as_array()
            .unwrap()
            .len(),
        1
    );
    let m = &check["assessment"]["reference"]["matches"][0];
    assert_eq!(m["fully_specified_dna"], true);
    let selection = json!({"batch_revision":r["batch_revision"],"candidate_id":m["candidate"],"component_id":m["product"]});
    let selected=call(&c,"select_recombination_component",json!({"sources":sources,"request":request,"selection":selection,"target":target,"policy":{"allow_reverse_complement":false}})).await;
    assert_eq!(selected["document"]["length"], 2553);
    assert_eq!(
        selected["reciprocal_components"].as_array().unwrap().len(),
        2
    );
    assert!(
        selected["authority"]
            .as_str()
            .unwrap()
            .contains("not physical")
    );
    let rejected_component = json!({"batch_revision":r["batch_revision"],"candidate_id":m["candidate"],"component_id":1-m["product"].as_u64().unwrap()});
    let before = call(&c, "list_documents", json!({})).await;
    reject(&c,"select_recombination_component",json!({"sources":sources,"request":request,"selection":rejected_component,"target":target,"policy":{"allow_reverse_complement":false}})).await;
    assert_eq!(call(&c, "list_documents", json!({})).await, before);
    c.cancel().await.unwrap();
    s.await.unwrap();
}
#[tokio::test]
async fn missing_explicit_profile_role_forgery_and_stale_batch_do_not_install_products() {
    let (c, s) = setup().await;
    let (sources, request, target) = prepare(&c).await;
    reject(&c, "inspect_recombination", json!({"sources":sources})).await;
    let mut forged = request.clone();
    forged["profile"]["revision"] = json!("0".repeat(64));
    reject(
        &c,
        "inspect_recombination",
        json!({"sources":sources,"request":forged}),
    )
    .await;
    let mut roles = sources.clone();
    roles[1]["role"] = roles[0]["role"].clone();
    reject(
        &c,
        "inspect_recombination",
        json!({"sources":roles,"request":request}),
    )
    .await;
    reject(&c,"select_recombination_component",json!({"sources":sources,"request":request,"target":target,"policy":{"allow_reverse_complement":false},"selection":{"batch_revision":"0".repeat(64),"candidate_id":0,"component_id":0}})).await;
    c.cancel().await.unwrap();
    s.await.unwrap();
}

#[tokio::test]
async fn capped_search_cannot_prove_absence_and_schema_requires_independent_target() {
    let (c, s) = setup().await;
    let (sources, mut request, target) = prepare(&c).await;
    request["limits"]["max_comparisons"] = json!(1);
    let checked = call(&c, "check_recombination_target", json!({"sources":sources,"request":request,"target":target,"policy":{"allow_reverse_complement":false}})).await;
    assert_eq!(checked["assessment"]["reference"]["search_complete"], false);
    assert_eq!(
        checked["assessment"]["reference"]["no_match_proven_under_model"],
        false
    );
    reject(&c,"select_recombination_component",json!({"sources":sources,"request":request,"selection":{"batch_revision":"0".repeat(64),"candidate_id":0,"component_id":0},"policy":{"allow_reverse_complement":false}})).await;
    let tools = c.list_tools(Default::default()).await.unwrap().tools;
    for name in [
        "inspect_recombination",
        "check_recombination_target",
        "select_recombination_component",
    ] {
        let tool = tools.iter().find(|t| t.name == name).unwrap();
        let schema = serde_json::to_value(&tool.input_schema).unwrap();
        let required = schema["required"].as_array().unwrap();
        assert!(required.contains(&json!("request")));
        assert!(required.contains(&json!("sources")));
        if name != "inspect_recombination" {
            assert!(required.contains(&json!("target")));
            assert!(required.contains(&json!("policy")));
        }
        assert!(!schema.to_string().contains("filesystem_path"));
    }
    c.cancel().await.unwrap();
    s.await.unwrap();
}

#[tokio::test]
async fn typed_plan_role_bridge_records_only_exact_causal_source_events() {
    use base64::{Engine as _, engine::general_purpose::STANDARD};
    use collagene_core::{
        history::DocumentHistory,
        run_lineage::EventBinding,
        workspace::{Workspace, save_workspace},
    };
    let (c, s) = setup().await;
    let (sources, request, target) = prepare(&c).await;
    let slots = json!([{"Input":"insert"},{"Input":"donor"}]);
    let inspected = call(
        &c,
        "inspect_recombination",
        json!({"sources":sources,"request":request,"plan_sources":slots}),
    )
    .await;
    assert_eq!(inspected["plan_sources"], slots);
    let checked=call(&c,"check_recombination_target",json!({"sources":sources,"request":request,"plan_sources":slots,"target":target,"policy":{"allow_reverse_complement":false}})).await;
    let m = &checked["assessment"]["reference"]["matches"][0];
    let selection = json!({"batch_revision":inspected["batch_revision"],"candidate_id":m["candidate"],"component_id":m["product"]});
    let proposed=call(&c,"propose_plan",json!({"name":"Source-event-bound reference exchange","target":target,"inputs":[{"id":"insert","document":sources[0]["document"]},{"id":"donor","document":sources[1]["document"]}],"steps":[{"RecombineReference":{"sources":slots,"request":request,"selection":selection}}]})).await;
    let verified = call(
        &c,
        "verify_plan",
        json!({"plan_revision":proposed["plan_revision"]}),
    )
    .await;
    assert_eq!(verified["status"], "Incomplete");
    let molecules = collagene_core::formats::import_genbank(include_str!(
        "../../../fixtures/recombination-bp-demo.gb"
    ))
    .unwrap()
    .molecules;
    let mut workspace = Workspace {
        name: "Independent target is noncausal".into(),
        active_document: "target".into(),
        documents: vec![],
        revisions: Default::default(),
        history: vec![],
        candidates: vec![],
        active_candidate: None,
        attachments: Default::default(),
    };
    let histories = [
        DocumentHistory::new("insert", molecules[0].clone()).unwrap(),
        DocumentHistory::new("donor", molecules[1].clone()).unwrap(),
        DocumentHistory::new("target", molecules[2].clone()).unwrap(),
    ];
    let binding = |i: usize| EventBinding {
        document_id: histories[i].id.clone(),
        event_id: histories[i].head_event().into(),
        revision: histories[i].head_revision().into(),
    };
    for history in &histories {
        history.contribute(&mut workspace).unwrap();
    }
    let imported = call(
        &c,
        "import_workspace_archive",
        json!({"data_base64":STANDARD.encode(save_workspace(&workspace).unwrap())}),
    )
    .await;
    let workspace_ref =
        json!({"id":imported["workspace"]["id"],"revision":imported["workspace"]["revision"]});
    let recorded=call(&c,"record_workspace_run",json!({"workspace":workspace_ref,"plan_revision":proposed["plan_revision"],"target_binding":binding(2),"input_bindings":[{"input_id":"insert","event":binding(0)},{"input_id":"donor","event":binding(1)}],"allow_incomplete":true})).await;
    assert_eq!(recorded["receipt"]["products"].as_array().unwrap().len(), 1);
    let product = &recorded["receipt"]["products"][0];
    let new_ref =
        json!({"id":recorded["workspace"]["id"],"revision":recorded["workspace"]["revision"]});
    let events = call(&c, "list_workspace_runs", new_ref.clone()).await;
    let mut text = String::new();
    let mut offset = 0;
    loop {
        let page = call(
            &c,
            "read_artifact",
            json!({"artifact_id":events["artifact_id"],"start":offset,"max_bytes":16384}),
        )
        .await;
        text.push_str(page["text"].as_str().unwrap());
        if page["done"] == true {
            break;
        }
        offset = page["next_start"].as_u64().unwrap();
    }
    let full: Value = serde_json::from_str(&text).unwrap();
    let product_event = full["events"]
        .as_array()
        .unwrap()
        .iter()
        .find(|e| e["id"] == product["event_id"])
        .unwrap();
    let parents = product_event["parents"].as_array().unwrap();
    assert_eq!(parents.len(), 2);
    assert!(parents.contains(&json!(binding(0).event_id)));
    assert!(parents.contains(&json!(binding(1).event_id)));
    assert!(!parents.contains(&json!(binding(2).event_id)));
    let replay = call(
        &c,
        "replay_workspace_run",
        json!({"workspace":new_ref,"run_id":recorded["receipt"]["run_id"]}),
    )
    .await;
    assert_eq!(replay["current_status"], "Incomplete");
    assert_eq!(replay["products_match"], true);
    c.cancel().await.unwrap();
    s.await.unwrap();
}
