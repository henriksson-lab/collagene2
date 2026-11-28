//! Atomic bulk release, native lineage and bounded portable ZIP over the SDK.
use base64::{Engine, engine::general_purpose::STANDARD};
use collagene_core::{
    chemistry::Terminus,
    cloning::assembly::recipe::Preparation,
    documents::workspace::{self, Workspace},
    material::lineage::EventBinding,
    planning::combinatorial::{
        Bin, Choice, ChoiceSelection, Design, Method, SelectedConstruct, TargetDeclaration,
    },
};
use collagene_mcp::server::CollageneServer;
use collagene_model::{Molecule, Topology};
use rmcp::{RoleClient, ServiceExt, model::CallToolRequestParams, service::RunningService};
use serde_json::{Value, json};
use std::io::{Cursor, Read};
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
async fn setup() -> (Client, tokio::task::JoinHandle<()>) {
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
    (().serve(b).await.unwrap(), server)
}
async fn text_artifact(c: &Client, id: &Value) -> Value {
    let (mut start, mut text) = (0, String::new());
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
async fn binary_artifact(c: &Client, id: &Value) -> Vec<u8> {
    let (mut start, mut text) = (0, String::new());
    loop {
        let page = call(
            c,
            "read_artifact",
            json!({"artifact_id":id,"start":start,"max_bytes":16384}),
        )
        .await;
        assert_eq!(page["encoding"], "base64");
        text.push_str(page["text"].as_str().unwrap());
        if page["done"] == true {
            break;
        }
        start = page["next_start"].as_u64().unwrap();
    }
    STANDARD.decode(text).unwrap()
}
fn binding(w: &Workspace, index: usize) -> EventBinding {
    let doc = &w.documents[index];
    EventBinding {
        document_id: doc.id.clone(),
        event_id: doc.head_event.clone().unwrap(),
        revision: doc.head.clone(),
    }
}
fn workspace_ref(value: &Value) -> Value {
    json!({"id":value["workspace"]["id"],"revision":value["workspace"]["revision"]})
}
async fn upload(c: &Client, w: &Workspace) -> Value {
    call(
        c,
        "import_workspace_archive",
        json!({"data_base64":STANDARD.encode(workspace::save_workspace(w).unwrap())}),
    )
    .await
}

fn fixture() -> (Workspace, Value) {
    let w = Workspace::from_documents(
        "bulk SDK",
        vec![
            Molecule::dna("left", "AACC", Topology::Linear).unwrap(),
            Molecule::dna("right", "TTAA", Topology::Linear).unwrap(),
            Molecule::dna("independent wrong target", "GGTT", Topology::Linear).unwrap(),
        ],
        0,
    )
    .unwrap();
    let d = Design {
        material_bindings: vec![],
        id: "bulk-sdk".into(),
        name: "Bulk SDK Ω".into(),
        circular: false,
        method: Method::Ligation {
            allow_nicked: false,
        },
        bins: (0..2)
            .map(|i| Bin {
                id: format!("bin-{i}"),
                label: format!("Bin {i}"),
                choices: vec![Choice {
                    anneal_bottom: None,
                    generation: None,
                    id: format!("choice-{i}"),
                    label: format!("Choice {i}"),
                    source: binding(&w, i),
                    selector: collagene_core::cloning::assembly::intent::Selector::ReferenceRange {
                        start: 0,
                        end: 4,
                        reverse: false,
                    },
                    preparation_reverse: false,
                    preparation: Some(Preparation::DeclaredDuplex {
                        left: Terminus::CLEAVED,
                        right: Terminus::CLEAVED,
                    }),
                    treatment: None,
                }],
            })
            .collect(),
    };
    let choices = (0..2)
        .map(|i| ChoiceSelection {
            bin_id: format!("bin-{i}"),
            choice_id: format!("choice-{i}"),
        })
        .collect::<Vec<_>>();
    let good = SelectedConstruct {
        standard: None,
        id: "..".into(),
        name: "Cloning Ω good".into(),
        choices: choices.clone(),
        target: TargetDeclaration::OriginalRegions {
            substitutions: vec![],
            edits: vec![],
            junction_edits: vec![],
        },
        overrides: None,
    };
    let bad = SelectedConstruct {
        standard: None,
        id: "wrong".into(),
        name: "Independent wrong target".into(),
        choices,
        target: TargetDeclaration::Existing(binding(&w, 2)),
        overrides: None,
    };
    (w, json!({"design":d,"rows":[{"Ready":good},{"Ready":bad}]}))
}
fn options(name: &str) -> Value {
    json!({"name":name,"show_prepared":false,"save_primers":false,"allow_incomplete":false})
}
#[tokio::test]
async fn native_release_and_zip_keep_exact_events_utf8_historical_labels_and_safe_paths() {
    let (c, server) = setup().await;
    let tools = c.list_all_tools().await.unwrap();
    assert_eq!(tools.len(), collagene_mcp::CURRENT_TOOL_COUNT);
    for name in [
        "release_combinatorial_outputs",
        "export_combinatorial_outputs",
        "list_combinatorial_releases",
        "open_combinatorial_release",
    ] {
        let tool = tools.iter().find(|tool| tool.name == name).unwrap();
        assert_eq!(tool.input_schema["additionalProperties"], false);
    }
    let (w, matrix) = fixture();
    let imported = upload(&c, &w).await;
    let source = workspace_ref(&imported);
    let batch = call(
        &c,
        "evaluate_combinatorial_batch",
        json!({"workspace":source,"matrix":matrix}),
    )
    .await;
    let batch_ref = workspace_ref(&batch);
    let requested = json!({"batch_hash":batch["receipt"]["attachment_hash"],"rows":[{"construct_id":"..","options":options("../Released β")} ]});
    let released = call(
        &c,
        "release_combinatorial_outputs",
        json!({"workspace":batch_ref,"request":requested}),
    )
    .await;
    assert_eq!(released["receipt"]["state"], "historical-output-lineage");
    assert_eq!(released["receipt"]["verified"], false);
    let release_ref = workspace_ref(&released);
    let hash = released["receipt"]["manifest_hash"].clone();
    let native =
        workspace::load_workspace(&binary_artifact(&c, &released["native_artifact_id"]).await)
            .unwrap();
    for doc in &w.documents {
        assert_eq!(native.head(&doc.id).unwrap(), w.head(&doc.id).unwrap());
    }
    let manifest = text_artifact(&c, &released["artifact_id"]).await;
    assert_eq!(manifest["rows"][0]["construct_id"], "..");
    let named = &manifest["rows"][0]["output"]["construct"];
    let molecule = &native.revisions[named["revision"].as_str().unwrap()];
    assert_eq!(molecule.sequence.as_str(), "AACCTTAA");
    assert_eq!(molecule.name, "../Released β");
    let listed = call(
        &c,
        "list_combinatorial_releases",
        json!({"workspace":release_ref,"limit":1}),
    )
    .await;
    assert_eq!(listed["total"], 1);
    assert_eq!(listed["verified"], false);
    let opened = call(
        &c,
        "open_combinatorial_release",
        json!({"workspace":release_ref,"manifest_hash":hash}),
    )
    .await;
    assert_eq!(opened["state"], "historical-output-lineage");
    assert_eq!(text_artifact(&c, &opened["artifact_id"]).await, manifest);
    let exported=call(&c,"export_combinatorial_outputs",json!({"workspace":release_ref,"manifest_hash":hash,"options":{"genbank":true,"fasta":true,"primer_csv":false,"primer_tsv":false,"protocol_html":true,"protocol_pdf":true,"construct_ids":[".."]}})).await;
    assert_eq!(exported["mime_type"], "application/zip");
    assert_eq!(exported["verified"], false);
    let bytes = binary_artifact(&c, &exported["artifact_id"]).await;
    assert_eq!(exported["blake3"], blake3::hash(&bytes).to_hex().as_str());
    let mut zip = zip::ZipArchive::new(Cursor::new(bytes)).unwrap();
    assert_eq!(zip.len(), 5);
    for i in 0..zip.len() {
        let mut file = zip.by_index(i).unwrap();
        assert!(file.enclosed_name().is_some());
        assert!(!file.name().contains("../") && !file.name().contains('\\'));
        let name = file.name().to_string();
        let mut data = vec![];
        file.read_to_end(&mut data).unwrap();
        if name.ends_with(".fasta") {
            assert!(std::str::from_utf8(&data).unwrap().contains("AACCTTAA"));
        }
        if name.ends_with(".html") {
            let html = std::str::from_utf8(&data).unwrap();
            assert!(html.contains("Cloning Ω good"));
            assert!(html.contains("ConsistentUnderAssumptions"));
        }
        if name == "manifest.json" {
            let value: Value = serde_json::from_slice(&data).unwrap();
            assert_eq!(value["release_manifest_hash"], hash);
            assert_eq!(
                value["rows"][0]["released"]["inputs"],
                manifest["rows"][0]["inputs"]
            );
        }
    }
    let stale = json!({"id":release_ref["id"],"revision":"0".repeat(64)});
    reject(
        &c,
        "export_combinatorial_outputs",
        json!({"workspace":stale,"manifest_hash":hash}),
    )
    .await;
    reject(
        &c,
        "release_combinatorial_outputs",
        json!({"workspace":source,"request":requested}),
    )
    .await;
    c.cancel().await.unwrap();
    server.await.unwrap();
}
#[tokio::test]
async fn a_late_wrong_target_discards_every_private_prefix_and_the_same_source_can_retry() {
    let (c, server) = setup().await;
    let (w, matrix) = fixture();
    let imported = upload(&c, &w).await;
    let batch = call(
        &c,
        "evaluate_combinatorial_batch",
        json!({"workspace":workspace_ref(&imported),"matrix":matrix}),
    )
    .await;
    let wr = workspace_ref(&batch);
    let before = call(&c, "list_documents", json!({})).await;
    let request = json!({"batch_hash":batch["receipt"]["attachment_hash"],"rows":[{"construct_id":"..","options":options("First private good")},{"construct_id":"wrong","options":options("Late wrong independent target")} ]});
    reject(
        &c,
        "release_combinatorial_outputs",
        json!({"workspace":wr,"request":request}),
    )
    .await;
    assert_eq!(call(&c, "list_documents", json!({})).await, before);
    let list = call(&c, "list_combinatorial_releases", json!({"workspace":wr})).await;
    assert_eq!(list["total"], 0);
    let request = json!({"batch_hash":batch["receipt"]["attachment_hash"],"rows":[{"construct_id":"..","options":options("Successful retry")} ]});
    let released = call(
        &c,
        "release_combinatorial_outputs",
        json!({"workspace":wr,"request":request}),
    )
    .await;
    assert_eq!(
        released["receipt"]["released_rows"]
            .as_array()
            .unwrap()
            .len(),
        1
    );
    let unchanged = call(&c, "list_combinatorial_releases", json!({"workspace":wr})).await;
    assert_eq!(unchanged["total"], 0);
    c.cancel().await.unwrap();
    server.await.unwrap();
}
