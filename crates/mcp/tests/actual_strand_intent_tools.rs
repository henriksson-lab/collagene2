//! Actual SDK intent requests hydrate only exact small native origin selectors.
#[path = "../../core/tests/common/actual_strand.rs"]
mod fixture;
use base64::{Engine as _, engine::general_purpose::STANDARD};
use collagene_core::{
    Molecule, Topology,
    history::DocumentHistory,
    run_material::MaterialSelection,
    workspace::{self, Workspace},
};
use collagene_mcp::CollageneServer;
use rmcp::{RoleClient, ServiceExt, model::CallToolRequestParams, service::RunningService};
use serde_json::{Value, json};
type Client = RunningService<RoleClient, ()>;
async fn call(c: &Client, name: &str, args: Value) -> Value {
    let reply = c
        .call_tool(
            CallToolRequestParams::new(name.to_owned())
                .with_arguments(args.as_object().unwrap().clone()),
        )
        .await
        .unwrap();
    assert_ne!(reply.is_error, Some(true), "{name}: {:?}", reply.content);
    reply.structured_content.unwrap()
}
async fn reject(c: &Client, name: &str, args: Value) {
    let reply = c
        .call_tool(
            CallToolRequestParams::new(name.to_owned())
                .with_arguments(args.as_object().unwrap().clone()),
        )
        .await;
    assert!(
        reply.is_err() || reply.unwrap().is_error == Some(true),
        "{name} accepted wrong native binding"
    );
}
async fn artifact(c: &Client, id: &Value) -> Value {
    let mut text = String::new();
    let mut offset = 0;
    loop {
        let p = call(
            c,
            "read_artifact",
            json!({"artifact_id":id,"start":offset,"max_bytes":16384}),
        )
        .await;
        text.push_str(p["text"].as_str().unwrap());
        if p["done"] == true {
            break;
        }
        offset = p["next_start"].as_u64().unwrap();
    }
    serde_json::from_str(&text).unwrap()
}
#[tokio::test]
async fn all_four_intents_keep_actual_bottom_before_new_preparation_and_reject_aliases() {
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
    let (mut w, selection, source, native) = fixture::fixture();
    let mut histories = DocumentHistory::restore_all(&w).unwrap();
    for h in &mut histories {
        if h.id == selection.event().document_id || h.id == w.documents[2].id {
            h.commit(
                Molecule::dna("Later unrelated viewer head", "TTTTGGGG", Topology::Linear).unwrap(),
                "Later edit after native event",
            )
            .unwrap();
        }
    }
    w.documents.clear();
    w.history.clear();
    for h in histories {
        h.contribute(&mut w).unwrap();
    }
    let alias = DocumentHistory::new("unrelated-equal-reference", source.molecule.clone()).unwrap();
    let alias_event = collagene_core::run_lineage::EventBinding {
        document_id: alias.id.clone(),
        event_id: alias.head_event().into(),
        revision: alias.head_revision().into(),
    };
    alias.contribute(&mut w).unwrap();
    let bytes = workspace::save_workspace(&w).unwrap();
    let imported = call(
        &c,
        "import_workspace_archive",
        json!({"data_base64":STANDARD.encode(&bytes)}),
    )
    .await;
    let wr = json!({"id":imported["workspace"]["id"],"revision":imported["workspace"]["revision"]});
    let inputs =
        json!([{"id":"saved","document":{"id":source.revision,"revision":source.revision}}]);
    let base = json!({"name":native.name,"fragments":native.fragments,"circular":native.circular,"material_bindings":[{"input_id":"saved","event":selection.event(),"selection":selection}]});
    let substituted = json!({"base":base,"junction_edits":[],"substitutions":[]});
    let edited = json!({"base":substituted,"edits":[]});
    let queries = [
        (
            "build_assembly_intent",
            json!({"inputs":inputs,"request":base,"origin_workspace":wr}),
        ),
        (
            "build_edited_assembly_intent",
            json!({"inputs":inputs,"base":base,"edits":[],"origin_workspace":wr}),
        ),
        (
            "build_substituted_assembly_intent",
            json!({"inputs":inputs,"request":substituted,"origin_workspace":wr}),
        ),
        (
            "build_contribution_edit_intent",
            json!({"inputs":inputs,"request":edited,"origin_workspace":wr}),
        ),
    ];
    for (name, args) in queries {
        let out = call(&c, name, args).await;
        assert_eq!(out["verified"], false);
        let sequence = call(
            &c,
            "read_sequence_region",
            json!({"document":out["target"],"start":0,"end":10}),
        )
        .await;
        assert_eq!(sequence["sequence"], "AATTCGTAAG");
        assert_ne!(sequence["sequence"], "CGTAAGAATT");
        let evidence = artifact(&c, &out["artifact_id"]).await;
        assert!(evidence.to_string().contains(&selection.event().event_id));
        assert_eq!(out["preparation_simulated"], false);
    }
    let valid = json!({"inputs":inputs,"request":base,"origin_workspace":wr});
    let count = call(&c, "list_documents", json!({})).await["total"].clone();
    let mut no_workspace = valid.clone();
    no_workspace["origin_workspace"] = Value::Null;
    reject(&c, "build_assembly_intent", no_workspace).await;
    let mut wrong_state = valid.clone();
    wrong_state["request"]["fragments"][0]["selector"]["MaterialStrandRange"]["expected_state_revision"] =
        json!("0".repeat(64));
    reject(&c, "build_assembly_intent", wrong_state).await;
    let mut alias = valid.clone();
    let mut alias_selection = selection.clone();
    match &mut alias_selection {
        MaterialSelection::Run(s) => s.event = alias_event.clone(),
        MaterialSelection::Derived(e) => *e = alias_event.clone(),
    };
    alias["request"]["material_bindings"][0]["event"] = json!(alias_event);
    alias["request"]["material_bindings"][0]["selection"] = json!(alias_selection);
    reject(&c, "build_assembly_intent", alias).await;
    let mut full_capsule = valid.clone();
    full_capsule["request"]["material_bindings"][0] = json!(native.material_bindings[0]);
    reject(&c, "build_assembly_intent", full_capsule).await;
    let mut missing = valid.clone();
    missing.as_object_mut().unwrap().remove("origin_workspace");
    reject(&c, "build_assembly_intent", missing).await;
    let foreign =
        Workspace::from_documents("No native producer", vec![source.molecule.clone()], 0).unwrap();
    let imported = call(
        &c,
        "import_workspace_archive",
        json!({"data_base64":STANDARD.encode(workspace::save_workspace(&foreign).unwrap())}),
    )
    .await;
    let mut no_owner = valid;
    no_owner["origin_workspace"] =
        json!({"id":imported["workspace"]["id"],"revision":imported["workspace"]["revision"]});
    reject(&c, "build_assembly_intent", no_owner).await;
    assert_eq!(call(&c, "list_documents", json!({})).await["total"], count);
    assert_eq!(workspace::save_workspace(&w).unwrap(), bytes);
    let _ = c.cancel().await;
    server.abort();
}
