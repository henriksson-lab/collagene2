//! Real MCP transport for exact construct/read event capture and fresh replay.
use base64::{Engine, engine::general_purpose::STANDARD};
use collagene_core::{
    documents::{
        history::DocumentHistory,
        workspace::{self, Workspace, WorkspaceAttachment},
    },
    formats::oligos::{self as oligo_formats, QualityEncoding},
    review::reads::{
        checks::ATTACHMENT_KIND,
        trace_validation::{TraceValidationOptions, TrimMode},
    },
};
use collagene_mcp::server::CollageneServer;
use collagene_model::{Molecule, Topology};
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
    assert_ne!(r.is_error, Some(true), "{name}:{:?}", r.content);
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
async fn artifact(c: &Client, id: &Value) -> Value {
    let (mut at, mut text) = (0, String::new());
    loop {
        let r = call(
            c,
            "read_artifact",
            json!({"artifact_id":id,"start":at,"max_bytes":16384}),
        )
        .await;
        text.push_str(r["text"].as_str().unwrap());
        if r["done"] == true {
            break;
        }
        at = r["next_start"].as_u64().unwrap();
    }
    serde_json::from_str(&text).unwrap()
}
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
fn binding(w: &Workspace, i: usize) -> Value {
    json!({"document_id":w.documents[i].id,"event_id":w.documents[i].head_event,"revision":w.documents[i].head})
}
fn workspace_ref(result: &Value) -> Value {
    json!({"id":result["workspace"]["id"],"revision":result["workspace"]["revision"]})
}
async fn native(c: &Client, id: &Value) -> Vec<u8> {
    let (mut at, mut text) = (0, String::new());
    loop {
        let page = call(
            c,
            "read_artifact",
            json!({"artifact_id":id,"start":at,"max_bytes":16384}),
        )
        .await;
        assert_eq!(page["encoding"], "base64");
        text.push_str(page["text"].as_str().unwrap());
        if page["done"] == true {
            break;
        }
        at = page["next_start"].as_u64().unwrap();
    }
    STANDARD.decode(text).unwrap()
}
fn fixture() -> Workspace {
    let reference = Molecule::dna("construct", "TTTCGATACGGG", Topology::Linear).unwrap();
    let read = oligo_formats::import_fastq(b"@read\nCGATAC\n+\nIIIIII\n", QualityEncoding::Phred33)
        .unwrap()
        .molecules
        .remove(0);
    Workspace::from_documents(
        "native read checks",
        vec![reference.clone(), reference, read.clone(), read],
        0,
    )
    .unwrap()
}
fn options() -> TraceValidationOptions {
    let mut options = TraceValidationOptions::default();
    options.trim.mode = TrimMode::AllCalls;
    options.trim.min_retained_bases = 1;
    options
}
async fn upload(c: &Client, w: &Workspace) -> Value {
    call(
        c,
        "import_workspace_archive",
        json!({"data_base64":STANDARD.encode(workspace::save_workspace(w).unwrap())}),
    )
    .await
}
#[tokio::test]
async fn exact_native_events_capture_partial_evidence_deduplicate_and_replay_without_changing_originals()
 {
    let (c, server) = setup().await;
    assert_eq!(
        c.list_tools(Default::default()).await.unwrap().tools.len(),
        collagene_mcp::CURRENT_TOOL_COUNT
    );
    let w = fixture();
    let bytes = workspace::save_workspace(&w).unwrap();
    let imported = upload(&c, &w).await;
    let wr = workspace_ref(&imported);
    let args = json!({"workspace":wr,"reference":binding(&w,0),"reads":[binding(&w,2)],"options":options()});
    let saved = call(&c, "record_read_check", args.clone()).await;
    assert_eq!(saved["receipt"]["state"], "historical");
    assert_eq!(saved["receipt"]["verified"], false);
    assert_eq!(saved["receipt"]["recorded_status"], "Incomplete");
    assert_eq!(saved["receipt"]["reference"], binding(&w, 0));
    assert_eq!(saved["receipt"]["read_count"], 1);
    let full = artifact(&c, &saved["artifact_id"]).await;
    assert_eq!(full["reference"], binding(&w, 0));
    assert_eq!(full["reads"], json!([binding(&w, 2)]));
    assert_eq!(full["report"]["consensus"], "NNNCGATACNNN");
    assert_eq!(full["report"]["eligible_coverage"], 6);
    assert_eq!(full["report"]["reference_coverage_fraction"], 0.5);
    assert_eq!(full["report"]["uncovered_ranges"], json!([[0, 3], [9, 12]]));
    let kept = workspace::load_workspace(&native(&c, &saved["native_artifact_id"]).await).unwrap();
    assert_eq!(kept.documents.len(), w.documents.len());
    assert_eq!(kept.history.len(), w.history.len());
    assert_eq!(kept.active_document, w.active_document);
    for d in &w.documents {
        assert_eq!(kept.head(&d.id).unwrap(), w.head(&d.id).unwrap());
    }
    let hash = saved["receipt"]["attachment_hash"].clone();
    let saved_wr = workspace_ref(&saved);
    let list = call(&c, "list_read_checks", json!({"workspace":saved_wr})).await;
    assert_eq!(list["state"], "historical");
    assert_eq!(list["verified"], false);
    assert_eq!(list["checks"].as_array().unwrap().len(), 1);
    assert_eq!(list["checks"][0]["attachment_hash"], hash);
    let replay = call(
        &c,
        "replay_read_check",
        json!({"workspace":saved_wr,"attachment_hash":hash}),
    )
    .await;
    assert_eq!(replay["state"], "fresh-computational-evidence");
    assert_eq!(replay["verified"], false);
    assert_eq!(replay["status"], "Incomplete");
    assert_eq!(replay["reference"], binding(&w, 0));
    assert_eq!(replay["eligible_coverage"], 6);
    assert_eq!(replay["reference_coverage_fraction"], 0.5);
    let fresh = artifact(&c, &replay["artifact_id"]).await;
    assert_eq!(fresh["saved"], full);
    assert_eq!(fresh["report"]["consensus"], "NNNCGATACNNN");
    assert_eq!(
        fresh["report"]["engine_version"],
        collagene_core::verification::ENGINE_VERSION
    );
    let mut retry = args.clone();
    retry["workspace"] = saved_wr.clone();
    let again = call(&c, "record_read_check", retry).await;
    assert_eq!(again["receipt"]["attachment_hash"], hash);
    assert_eq!(again["artifact_id"], saved["artifact_id"]);
    assert_eq!(workspace_ref(&again), saved_wr);
    // Equal DNA does not authorize substituting a distinct native reference event.
    let mut alias = args.clone();
    alias["workspace"] = saved_wr.clone();
    alias["reference"] = binding(&w, 1);
    let alias_check = call(&c, "record_read_check", alias).await;
    assert_ne!(alias_check["receipt"]["attachment_hash"], hash);
    assert_eq!(alias_check["receipt"]["reference"], binding(&w, 1));
    assert_eq!(
        call(
            &c,
            "list_read_checks",
            json!({"workspace":workspace_ref(&alias_check)})
        )
        .await["checks"]
            .as_array()
            .unwrap()
            .len(),
        2
    );
    let original_export = call(&c, "export_workspace_archive", wr.clone()).await;
    assert_eq!(native(&c, &original_export["artifact_id"]).await, bytes);
    assert!(
        call(&c, "list_read_checks", json!({"workspace":wr})).await["checks"]
            .as_array()
            .unwrap()
            .is_empty()
    );
    let before = call(&c, "list_workspaces", json!({})).await["workspaces"]
        .as_array()
        .unwrap()
        .len();
    let mut stale = args.clone();
    stale["workspace"]["revision"] = json!("0".repeat(64));
    reject(&c, "record_read_check", stale).await;
    let mut bad = args.clone();
    bad["reference"]["document_id"] = binding(&w, 1)["document_id"].clone();
    reject(&c, "record_read_check", bad).await;
    let mut stale_read = args;
    stale_read["reads"][0]["revision"] = json!("0".repeat(64));
    reject(&c, "record_read_check", stale_read).await;
    reject(
        &c,
        "replay_read_check",
        json!({"workspace":saved_wr,"attachment_hash":"0".repeat(64)}),
    )
    .await;
    assert_eq!(
        call(&c, "list_workspaces", json!({})).await["workspaces"]
            .as_array()
            .unwrap()
            .len(),
        before
    );
    drop(c);
    server.abort();
}
#[tokio::test]
async fn archived_engine_and_user_edited_heads_never_restore_a_current_verdict_or_replace_the_saved_reference()
 {
    let (c, server) = setup().await;
    let w = fixture();
    let imported = upload(&c, &w).await;
    let captured=call(&c,"record_read_check",json!({"workspace":workspace_ref(&imported),"reference":binding(&w,0),"reads":[binding(&w,2)],"options":options()})).await;
    let mut archive =
        workspace::load_workspace(&native(&c, &captured["native_artifact_id"]).await).unwrap();
    let hash = captured["receipt"]["attachment_hash"].as_str().unwrap();
    let mut saved = collagene_core::review::reads::checks::load(&archive, hash).unwrap();
    saved.report.engine_version = "archived-untrusted-engine".into();
    saved.context_revision = blake3::hash(
        &serde_json::to_vec(&(
            "read-check-evidence-1",
            &saved.reference,
            &saved.reads,
            &saved.report.options,
            &saved.report.engine_version,
            &saved.report.model,
        ))
        .unwrap(),
    )
    .to_hex()
    .to_string();
    let bytes = serde_json::to_vec(&saved).unwrap();
    let new_hash = blake3::hash(&bytes).to_hex().to_string();
    archive.attachments.remove(hash);
    archive.attachments.insert(
        new_hash.clone(),
        WorkspaceAttachment {
            kind: ATTACHMENT_KIND.into(),
            name: "Archived read check assertion".into(),
            media_type: "application/json".into(),
            bytes,
        },
    );
    let mut histories = DocumentHistory::restore_all(&archive).unwrap();
    histories[0]
        .commit(
            Molecule::dna("user edited construct", "AAAA", Topology::Linear).unwrap(),
            "User edit after read check",
        )
        .unwrap();
    histories[2]
        .commit(
            Molecule::dna("user edited read", "CCCC", Topology::Linear).unwrap(),
            "User edit after read check",
        )
        .unwrap();
    archive.documents.clear();
    archive.history.clear();
    for h in &histories {
        h.contribute(&mut archive).unwrap();
    }
    archive.active_document = w.documents[2].id.clone();
    let imported = upload(&c, &archive).await;
    let wr = workspace_ref(&imported);
    let listed = call(&c, "list_read_checks", json!({"workspace":wr})).await;
    assert_eq!(listed["state"], "historical");
    assert_eq!(
        listed["checks"][0]["recorded_engine"],
        "archived-untrusted-engine"
    );
    assert_eq!(listed["checks"][0]["reference"], binding(&w, 0));
    let checked = call(
        &c,
        "replay_read_check",
        json!({"workspace":wr,"attachment_hash":new_hash}),
    )
    .await;
    assert_eq!(checked["reference"], binding(&w, 0));
    assert_eq!(checked["reads"], json!([binding(&w, 2)]));
    assert_eq!(
        checked["engine_version"],
        collagene_core::verification::ENGINE_VERSION
    );
    assert_eq!(checked["status"], "Incomplete");
    let full = artifact(&c, &checked["artifact_id"]).await;
    assert_eq!(
        full["saved"]["report"]["engine_version"],
        "archived-untrusted-engine"
    );
    assert_eq!(
        full["report"]["engine_version"],
        collagene_core::verification::ENGINE_VERSION
    );
    assert_eq!(full["report"]["consensus"], "NNNCGATACNNN");
    let exported = call(&c, "export_workspace_archive", wr).await;
    let unchanged = workspace::load_workspace(&native(&c, &exported["artifact_id"]).await).unwrap();
    assert_eq!(unchanged.active_document, w.documents[2].id);
    assert_eq!(
        unchanged
            .head(&w.documents[0].id)
            .unwrap()
            .sequence
            .as_str(),
        "AAAA"
    );
    assert_eq!(
        unchanged
            .head(&w.documents[2].id)
            .unwrap()
            .sequence
            .as_str(),
        "CCCC"
    );
    drop(c);
    server.abort();
}
