//! Real SDK routing, immutable workspace commits and exact source identities.
use base64::{Engine, engine::general_purpose::STANDARD};
use collagene_core::{
    Molecule, Topology,
    assembly_recipe::Preparation,
    chemistry::Terminus,
    combinatorial_design::{
        self as design, Bin, Choice, ChoiceSelection, Design, Method, SelectedConstruct,
        TargetDeclaration,
    },
    history::DocumentHistory,
    run_lineage::EventBinding,
    workspace::{self, Workspace},
};
use collagene_mcp::CollageneServer;
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
async fn native_artifact(c: &Client, id: &Value) -> Workspace {
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
    workspace::load_workspace(&STANDARD.decode(text).unwrap()).unwrap()
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
fn fixture() -> (Workspace, Design, SelectedConstruct) {
    let first = Molecule::dna("left", "AACC", Topology::Linear).unwrap();
    let original = Workspace::from_documents(
        "MCP alternatives",
        vec![
            first.clone(),
            Molecule::dna("right", "TTAA", Topology::Linear).unwrap(),
            Molecule::dna("independent wrong target", "GGTT", Topology::Linear).unwrap(),
            first,
        ],
        0,
    )
    .unwrap();
    let choice = |id: &str, index: usize| Choice {
        generation: None,
        id: id.into(),
        label: id.into(),
        source: binding(&original, index),
        selector: collagene_core::assembly_intent::Selector::ReferenceRange {
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
    };
    let d = Design {
        material_bindings: vec![],
        id: "mcp-library".into(),
        name: "MCP alternative sources".into(),
        circular: false,
        bins: vec![
            Bin {
                id: "left-bin".into(),
                label: "Backbone".into(),
                choices: vec![choice("left-original", 0), choice("left-alias", 3)],
            },
            Bin {
                id: "right-bin".into(),
                label: "Insert".into(),
                choices: vec![choice("right-choice", 1)],
            },
        ],
        method: Method::Ligation {
            allow_nicked: false,
        },
    };
    let row = SelectedConstruct {
        id: "literal-row".into(),
        name: "Literal AACC TTAA".into(),
        choices: vec![
            ChoiceSelection {
                bin_id: "left-bin".into(),
                choice_id: "left-original".into(),
            },
            ChoiceSelection {
                bin_id: "right-bin".into(),
                choice_id: "right-choice".into(),
            },
        ],
        target: TargetDeclaration::OriginalRegions {
            substitutions: vec![],
            edits: vec![],
            junction_edits: vec![],
        },
        overrides: None,
    };
    // The chosen original event remains retained after a later source-head edit.
    let mut histories = DocumentHistory::restore_all(&original).unwrap();
    histories[0]
        .commit(
            Molecule::dna("later left", "CCCC", Topology::Linear).unwrap(),
            "Later edit",
        )
        .unwrap();
    let mut historical = original;
    historical.documents.clear();
    historical.history.clear();
    for doc in histories {
        doc.contribute(&mut historical).unwrap();
    }
    historical.validate().unwrap();
    (historical, d, row)
}

#[tokio::test]
async fn matrix_inspection_keeps_historical_sources_independent_targets_and_aliases_distinct() {
    let (c, server) = setup().await;
    let tools = c.list_all_tools().await.unwrap();
    assert_eq!(tools.len(), collagene_mcp::CURRENT_TOOL_COUNT);
    for name in [
        "page_combinatorial_design",
        "inspect_combinatorial_construct",
        "evaluate_combinatorial_batch",
        "list_combinatorial_batches",
        "open_combinatorial_batch",
        "replay_combinatorial_batch",
    ] {
        let tool = tools.iter().find(|t| t.name == name).unwrap();
        assert_eq!(tool.input_schema["additionalProperties"], false);
    }
    let (w, d, row) = fixture();
    let original_source = d.bins[0].choices[0].source.clone();
    let imported = upload(&c, &w).await;
    let wr = workspace_ref(&imported);
    let page = call(
        &c,
        "page_combinatorial_design",
        json!({"workspace":wr,"design":d,"first":0,"limit":1}),
    )
    .await;
    assert_eq!(page["page"]["total"], 2);
    assert_eq!(page["page"]["combinations"].as_array().unwrap().len(), 1);
    assert_eq!(page["verified"], false);
    let inspected = call(
        &c,
        "inspect_combinatorial_construct",
        json!({"workspace":wr,"design":d,"construct":row}),
    )
    .await;
    let evidence = text_artifact(&c, &inspected["artifact_id"]).await;
    assert_eq!(
        evidence["input_bindings"][0]["event"],
        json!(original_source)
    );
    assert_eq!(
        evidence["inputs"][0]["molecule"]["sequence"]["bases"],
        "AACC"
    );
    assert_eq!(
        evidence["independent_target"]["sequence"]["bases"],
        "AACCTTAA"
    );
    assert_eq!(evidence["verified"], false);
    let mut alias = row.clone();
    alias.choices[0].choice_id = "left-alias".into();
    let second = call(
        &c,
        "inspect_combinatorial_construct",
        json!({"workspace":wr,"design":d,"construct":alias}),
    )
    .await;
    assert_ne!(
        second["inputs"][0]["input_id"],
        inspected["inputs"][0]["input_id"]
    );
    assert_ne!(
        second["input_bindings"][0]["event"],
        inspected["input_bindings"][0]["event"]
    );
    let mut wrong = row.clone();
    wrong.id = "wrong-target".into();
    wrong.target = TargetDeclaration::Existing(binding(&w, 2));
    let wrong = call(
        &c,
        "inspect_combinatorial_construct",
        json!({"workspace":wr,"design":d,"construct":wrong}),
    )
    .await;
    assert_eq!(wrong["target"]["binding"], json!(binding(&w, 2)));
    assert_eq!(wrong["target"]["length"], 4);
    assert_eq!(wrong["verified"], false);
    let mut owner = d.clone();
    owner.bins[0].choices[0].source.document_id = w.documents[3].id.clone();
    reject(
        &c,
        "page_combinatorial_design",
        json!({"workspace":wr,"design":owner,"limit":1}),
    )
    .await;
    reject(
        &c,
        "inspect_combinatorial_construct",
        json!({"workspace":wr,"design":owner,"construct":row}),
    )
    .await;
    reject(
        &c,
        "page_combinatorial_design",
        json!({"workspace":{"id":wr["id"],"revision":"0".repeat(64)},"design":d,"limit":1}),
    )
    .await;
    reject(
        &c,
        "page_combinatorial_design",
        json!({"workspace":wr,"design":d,"limit":65}),
    )
    .await;
    c.cancel().await.unwrap();
    server.await.unwrap();
}

#[tokio::test]
async fn fresh_selected_batches_save_historical_records_and_replay_exact_rows_without_rewinding_heads()
 {
    let (c, server) = setup().await;
    let (w, d, row) = fixture();
    let imported = upload(&c, &w).await;
    let wr = workspace_ref(&imported);
    let mut wrong = row.clone();
    wrong.id = "wrong-row".into();
    wrong.name = "Preserve independent wrong target".into();
    wrong.target = TargetDeclaration::Existing(binding(&w, 2));
    let mut missing = row.clone();
    missing.id = "missing-prep".into();
    missing.name = "Explicit unresolved preparation override".into();
    let mut missing = serde_json::to_value(missing).unwrap();
    missing["overrides"] = json!({"preparations":[{"generation":null,"bin_id":"right-bin","preparation":null,"treatment":null}],"method":null,"circular":null});
    let matrix = json!({"design":d,"rows":[{"Ready":row},{"Ready":wrong},{"Ready":missing},{"DraftError":{"id":"unfinished-row","name":"Unfinished literal primer","choices":row.choices,"message":"Literal primer input has an unfinished ? call; no replacement or threshold relaxation."}}]});
    let args = json!({"workspace":wr,"matrix":matrix});
    let saved = call(&c, "evaluate_combinatorial_batch", args.clone()).await;
    assert_eq!(saved["receipt"]["state"], "historical");
    assert_eq!(saved["receipt"]["verified"], false);
    assert_eq!(saved["receipt"]["recorded_row_count"], 4);
    let full = text_artifact(&c, &saved["artifact_id"]).await;
    assert_eq!(full["schema_version"], 1);
    assert_eq!(full["rows"].as_array().unwrap().len(), 4);
    assert_eq!(full["request"], matrix);
    assert_eq!(full["rows"][0]["state"], "Evaluated");
    assert_eq!(full["rows"][0]["status"], "ConsistentUnderAssumptions");
    assert_eq!(full["rows"][0]["preparation_complete"], true);
    assert_eq!(
        full["rows"][0]["input_bindings"][0]["event"],
        json!(d.bins[0].choices[0].source)
    );
    assert_eq!(full["rows"][1]["state"], "Evaluated");
    assert_eq!(full["rows"][1]["status"], "Invalid");
    assert_eq!(full["rows"][1]["target_binding"], json!(binding(&w, 2)));
    assert_eq!(full["rows"][2]["state"], "Error");
    assert!(full["rows"][2]["status"].is_null());
    assert!(
        full["rows"][2]["error"]
            .as_str()
            .unwrap()
            .contains("preparation")
    );
    assert_eq!(full["rows"][3]["state"], "DraftError");
    assert!(full["rows"][3]["status"].is_null());
    assert!(
        full["rows"][3]["error"]
            .as_str()
            .unwrap()
            .contains("unfinished ? call")
    );
    let native = native_artifact(&c, &saved["native_artifact_id"]).await;
    assert_eq!(native.documents.len(), w.documents.len());
    assert_eq!(native.history.len(), w.history.len());
    assert_eq!(native.active_document, w.active_document);
    for doc in &w.documents {
        assert_eq!(native.head(&doc.id).unwrap(), w.head(&doc.id).unwrap());
    }
    let attachment = saved["receipt"]["attachment_hash"].clone();
    let saved_wr = workspace_ref(&saved);
    assert!(
        native
            .attachments
            .contains_key(attachment.as_str().unwrap())
    );
    let list = call(
        &c,
        "list_combinatorial_batches",
        json!({"workspace":saved_wr,"offset":0,"limit":1}),
    )
    .await;
    assert_eq!(list["state"], "historical");
    assert_eq!(list["total"], 1);
    assert_eq!(list["verified"], false);
    assert_eq!(list["batches"].as_array().unwrap().len(), 1);
    let opened = call(
        &c,
        "open_combinatorial_batch",
        json!({"workspace":saved_wr,"attachment_hash":attachment}),
    )
    .await;
    assert_eq!(opened["state"], "historical");
    assert_eq!(opened["verified"], false);
    assert_eq!(text_artifact(&c, &opened["artifact_id"]).await, full);
    let replay = call(
        &c,
        "replay_combinatorial_batch",
        json!({"workspace":saved_wr,"attachment_hash":attachment,"row_index":0}),
    )
    .await;
    assert_eq!(replay["state"], "fresh-computational-evidence");
    assert_eq!(replay["verified"], false);
    let checked = text_artifact(&c, &replay["artifact_id"]).await;
    assert_eq!(checked["record_hash"], attachment);
    assert_eq!(checked["row_index"], 0);
    assert_eq!(checked["recorded"], full["rows"][0]);
    assert_eq!(checked["current"]["state"], "Evaluated");
    assert_eq!(checked["current"]["status"], "ConsistentUnderAssumptions");
    let imported_again = upload(&c, &native).await;
    let reopened = call(
        &c,
        "open_combinatorial_batch",
        json!({"workspace":workspace_ref(&imported_again),"attachment_hash":attachment}),
    )
    .await;
    assert_eq!(reopened["state"], "historical");
    assert_eq!(reopened["verified"], false);
    let repeated = call(&c, "evaluate_combinatorial_batch", args).await;
    assert_eq!(repeated["receipt"]["attachment_hash"], attachment);
    assert_eq!(workspace_ref(&repeated), saved_wr);
    reject(
        &c,
        "replay_combinatorial_batch",
        json!({"workspace":saved_wr,"attachment_hash":attachment,"row_index":4}),
    )
    .await;
    reject(
        &c,
        "open_combinatorial_batch",
        json!({"workspace":wr,"attachment_hash":attachment}),
    )
    .await;
    let mut forged = matrix.clone();
    forged["verified"] = json!(true);
    reject(
        &c,
        "evaluate_combinatorial_batch",
        json!({"workspace":wr,"matrix":forged}),
    )
    .await;
    reject(
        &c,
        "evaluate_combinatorial_batch",
        json!({"workspace":{"id":wr["id"],"revision":"0".repeat(64)},"matrix":matrix}),
    )
    .await;
    c.cancel().await.unwrap();
    server.await.unwrap();
}

#[tokio::test]
async fn internal_library_declarations_are_row_local_and_never_rewrite_preparations() {
    use collagene_core::{assembly_edit_intent::Edit, assembly_substitution_intent::Substitution};
    let (c, server) = setup().await;
    let (w, d, ordinary) = fixture();
    let imported = upload(&c, &w).await;
    let wr = workspace_ref(&imported);
    let (_, role) =
        design::source_role_ids(&d.id, "snp-row", &d.bins[0].id, &d.bins[0].choices[0].id).unwrap();
    let mut snp = ordinary.clone();
    snp.id = "snp-row".into();
    snp.name = "Only this row requests A to G".into();
    snp.target = TargetDeclaration::OriginalRegions {
        junction_edits: vec![],
        substitutions: vec![Substitution {
            role: role.clone(),
            offset: 1,
            expected: "A".into(),
            replacement: "G".into(),
        }],
        edits: vec![],
    };
    let mut indel = ordinary.clone();
    indel.id = "indel-row".into();
    indel.name = "Only this row replaces AC with T".into();
    indel.target = TargetDeclaration::OriginalRegions {
        junction_edits: vec![],
        substitutions: vec![],
        edits: vec![Edit {
            role: role.clone(),
            start: 1,
            end: 3,
            expected: "AC".into(),
            replacement: "T".into(),
        }],
    };
    let mut stale = snp.clone();
    stale.id = "stale-row".into();
    stale.name = "Stale expected base is a row-local error".into();
    stale.target = TargetDeclaration::OriginalRegions {
        junction_edits: vec![],
        substitutions: vec![Substitution {
            role,
            offset: 1,
            expected: "C".into(),
            replacement: "G".into(),
        }],
        edits: vec![],
    };
    for (row, literal) in [
        (&ordinary, "AACCTTAA"),
        (&snp, "AGCCTTAA"),
        (&indel, "ATCTTAA"),
    ] {
        let inspected = call(
            &c,
            "inspect_combinatorial_construct",
            json!({"workspace":wr,"design":d,"construct":row}),
        )
        .await;
        let evidence = text_artifact(&c, &inspected["artifact_id"]).await;
        assert_eq!(evidence["independent_target"]["sequence"]["bases"], literal);
        assert_eq!(
            evidence["input_bindings"][0]["event"],
            json!(d.bins[0].choices[0].source)
        );
    }
    let before = workspace::save_workspace(&w).unwrap();
    let matrix = json!({"design":d,"rows":[{"Ready":ordinary},{"Ready":snp},{"Ready":indel},{"Ready":stale}]});
    let saved = call(
        &c,
        "evaluate_combinatorial_batch",
        json!({"workspace":wr,"matrix":matrix}),
    )
    .await;
    let record = text_artifact(&c, &saved["artifact_id"]).await;
    assert_eq!(record["request"], matrix);
    assert_eq!(record["rows"][0]["status"], "ConsistentUnderAssumptions");
    for row in &record["rows"].as_array().unwrap()[1..3] {
        assert_eq!(row["state"], "Evaluated");
        assert_eq!(
            row["status"], "Invalid",
            "unedited material cannot satisfy a new intended mutation: {row}"
        );
    }
    assert_eq!(record["rows"][3]["state"], "Error");
    assert!(record["rows"][3]["status"].is_null());
    let restored = native_artifact(&c, &saved["native_artifact_id"]).await;
    assert_eq!(restored.history.len(), w.history.len());
    assert_eq!(restored.documents.len(), w.documents.len());
    for doc in &w.documents {
        assert_eq!(restored.head(&doc.id).unwrap(), w.head(&doc.id).unwrap());
    }
    assert_eq!(workspace::save_workspace(&w).unwrap(), before);
    let mut obsolete_shape = matrix;
    obsolete_shape["rows"][0]["Ready"]["target"]["OriginalRegions"]
        .as_object_mut()
        .unwrap()
        .remove("edits");
    reject(
        &c,
        "evaluate_combinatorial_batch",
        json!({"workspace":wr,"matrix":obsolete_shape}),
    )
    .await;
    c.cancel().await.unwrap();
    server.await.unwrap();
}
