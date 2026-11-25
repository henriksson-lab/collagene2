//! Immutable recovery bytes exposed through the real SDK, without restore authority.
use base64::{Engine, engine::general_purpose::STANDARD};
use collagene_core::{
    Molecule, Topology, assembly_library_history as history,
    run_lineage::EventBinding,
    workspace::{self, Workspace, WorkspaceAttachment},
};
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
    let r = c
        .call_tool(
            CallToolRequestParams::new(name.to_owned())
                .with_arguments(args.as_object().unwrap().clone()),
        )
        .await;
    assert!(r.is_err() || r.unwrap().is_error == Some(true));
}
async fn text(c: &Client, id: &Value) -> String {
    let (mut start, mut value) = (0, String::new());
    loop {
        let p = call(
            c,
            "read_artifact",
            json!({"artifact_id":id,"start":start,"max_bytes":4096}),
        )
        .await;
        assert_eq!(p["encoding"], "utf-8");
        value.push_str(p["text"].as_str().unwrap());
        if p["done"] == true {
            return value;
        }
        start = p["next_start"].as_u64().unwrap();
    }
}
fn binding(w: &Workspace, index: usize) -> EventBinding {
    let d = &w.documents[index];
    EventBinding {
        document_id: d.id.clone(),
        event_id: d.head_event.clone().unwrap(),
        revision: d.head.clone(),
    }
}
#[tokio::test]
async fn sdk_lists_and_opens_exact_raw_recovery_without_restoring_any_state() {
    let mut w = Workspace::from_documents(
        "historical recovery SDK",
        vec![
            Molecule::dna("source one", "AACCGGTT", Topology::Linear).unwrap(),
            Molecule::dna("source two", "TTGGCCAA", Topology::Linear).unwrap(),
        ],
        0,
    )
    .unwrap();
    // Deliberately opaque host data. Native history attests the version/counter
    // header and declared event pins; these bytes are not claimed as a valid WebDraft.
    let persisted =
        br#"{"version":3,"alternatives":{"id":"library-sdk","next_id":7},"pending":"persisted"}"#
            .to_vec();
    let before="{\n  \"version\":3,\n  \"alternatives\":{\"id\":\"library-sdk\",\"next_id\":7},\n  \"pending\":\"unsaved Ω primer: ACGT??  \",\n  \"removed_sources\":[\"source two\"]\n}\n".as_bytes().to_vec();
    let after =
        br#"{"version":3,"alternatives":{"id":"library-sdk","next_id":9},"pending":"after"}"#
            .to_vec();
    let oldhash = blake3::hash(&persisted).to_hex().to_string();
    w.attachments.insert(
        oldhash,
        WorkspaceAttachment {
            kind: history::ACTIVE_DRAFT_KIND.into(),
            name: "assembly-draft.json".into(),
            media_type: "application/json".into(),
            bytes: persisted,
        },
    );
    let original_heads = w
        .documents
        .iter()
        .map(|d| d.head.clone())
        .collect::<Vec<_>>();
    let source = binding(&w, 1);
    let mut wrong = source.clone();
    wrong.document_id = w.documents[0].id.clone();
    assert!(
        history::create(
            &w,
            before.clone(),
            &after,
            &"a".repeat(64),
            vec![wrong],
            vec![]
        )
        .is_err(),
        "event ownership cannot be relabelled to another retained document"
    );
    let record = history::create(
        &w,
        before.clone(),
        &after,
        &"a".repeat(64),
        vec![source.clone()],
        vec![],
    )
    .unwrap();
    assert_ne!(
        record.before_active_hash, record.before_hash,
        "persisted active and raw working draft are separate"
    );
    let commit = history::install(w, record, after.clone()).unwrap();
    let hash = commit.record_hash.clone();
    let native = workspace::save_workspace(&commit.workspace).unwrap();
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
        "list_assembly_library_history",
        "open_assembly_library_history",
    ] {
        let tool = tools.iter().find(|t| t.name == name).unwrap();
        assert_eq!(tool.input_schema["additionalProperties"], false);
    }
    assert!(
        !tools
            .iter()
            .any(|t| t.name == "restore_assembly_library_history")
    );
    let imported = call(
        &c,
        "import_workspace_archive",
        json!({"data_base64":STANDARD.encode(native)}),
    )
    .await;
    let reference =
        json!({"id":imported["workspace"]["id"],"revision":imported["workspace"]["revision"]});
    let listed = call(
        &c,
        "list_assembly_library_history",
        json!({"workspace":reference,"offset":0,"limit":1}),
    )
    .await;
    assert_eq!(listed["total"], 1);
    assert_eq!(listed["verified"], false);
    assert_eq!(listed["records"][0]["attachment_hash"], hash);
    assert_eq!(listed["records"][0]["allocation_highwater"], 9);
    let opened = call(
        &c,
        "open_assembly_library_history",
        json!({"workspace":reference,"attachment_hash":hash}),
    )
    .await;
    assert_eq!(opened["state"], "historical-edit-recovery");
    assert_eq!(opened["verified"], false);
    let raw = text(&c, &opened["before_draft_artifact_id"]).await;
    assert_eq!(raw.as_bytes(), before);
    assert_eq!(
        blake3::hash(raw.as_bytes()).to_hex().as_str(),
        opened["before_hash"]
    );
    let metadata: Value =
        serde_json::from_str(&text(&c, &opened["metadata_artifact_id"]).await).unwrap();
    assert!(
        metadata.get("before_draft").is_none(),
        "raw bytes are never inflated into integer JSON"
    );
    assert_eq!(metadata["sources"], json!([source]));
    assert_eq!(
        metadata["after_hash"],
        blake3::hash(&after).to_hex().as_str()
    );
    assert_eq!(metadata["before_next_id"], 7);
    assert_eq!(metadata["after_next_id"], 9);
    let empty = call(
        &c,
        "list_assembly_library_history",
        json!({"workspace":reference,"offset":1,"limit":1}),
    )
    .await;
    assert_eq!(empty["records"], json!([]));
    reject(
        &c,
        "list_assembly_library_history",
        json!({"workspace":reference,"limit":65}),
    )
    .await;
    reject(
        &c,
        "list_assembly_library_history",
        json!({"workspace":reference,"offset":2}),
    )
    .await;
    reject(
        &c,
        "open_assembly_library_history",
        json!({"workspace":reference,"attachment_hash":"f".repeat(64)}),
    )
    .await;
    // Open/list are immutable: native source heads and the current active bytes
    // remain untouched; the returned raw historical bytes are never installed.
    assert_eq!(
        commit
            .workspace
            .documents
            .iter()
            .map(|d| d.head.clone())
            .collect::<Vec<_>>(),
        original_heads
    );
    assert!(
        commit
            .workspace
            .attachments
            .values()
            .any(|a| a.kind == history::ACTIVE_DRAFT_KIND && a.bytes == after)
    );
    let repeated = call(
        &c,
        "open_assembly_library_history",
        json!({"workspace":reference,"attachment_hash":hash}),
    )
    .await;
    assert_eq!(
        repeated["before_draft_artifact_id"],
        opened["before_draft_artifact_id"]
    );
    drop(c);
    server.abort();
}
