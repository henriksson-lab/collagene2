//! Actual native-source read decisions and derived reference history through MCP SDK.
use base64::{Engine, engine::general_purpose::STANDARD};
use collagene_core::{
    documents::workspace::{self, Workspace},
    formats::oligos::{self as oligo_formats, QualityEncoding},
    material::lineage::EventBinding,
    review::reads::{
        checks as read_checks,
        trace_validation::{SiteState, TraceValidationOptions, TrimMode},
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
async fn upload(c: &Client, w: &Workspace) -> Value {
    call(
        c,
        "import_workspace_archive",
        json!({"data_base64":STANDARD.encode(workspace::save_workspace(w).unwrap())}),
    )
    .await
}
fn event(w: &Workspace, i: usize) -> EventBinding {
    serde_json::from_value(binding(w, i)).unwrap()
}
fn fixture(scored: bool) -> (Workspace, String) {
    let reference = Molecule::dna("original construct", "TTTCGATACGGG", Topology::Linear).unwrap();
    let raw = if scored {
        b"@observed\nTTTCGTTACGGG\n+\nIIIIIIIIIIII\n".as_slice()
    } else {
        b"@observed\nTTTCGTTACGGG\n+\n!!!!!!!!!!!!\n".as_slice()
    };
    let read = oligo_formats::import_fastq(raw, QualityEncoding::Phred33)
        .unwrap()
        .molecules
        .remove(0);
    let w = Workspace::from_documents("read review fixture", vec![reference, read], 0).unwrap();
    let mut options = TraceValidationOptions::default();
    options.trim.mode = TrimMode::AllCalls;
    options.trim.min_retained_bases = 1;
    let captured = read_checks::capture(&w, event(&w, 0), vec![event(&w, 1)], &options).unwrap();
    let saved = read_checks::load(&captured.workspace, &captured.attachment_hash).unwrap();
    if scored {
        let s = saved.report.sites.iter().find(|s| s.position == 5).unwrap();
        assert_eq!(s.reference_base, 'A');
        assert_eq!(s.consensus_base, 'T');
        assert_eq!(s.state, SiteState::SupportedDifference);
        let col = saved.report.reads[0]
            .columns
            .iter()
            .find(|c| c.reference_position == Some(5))
            .unwrap();
        assert_eq!(col.call_index, Some(5));
        assert_eq!(col.observed_base, Some('T'));
        assert!(col.eligible);
    }
    (captured.workspace, captured.attachment_hash)
}
fn request(hash: &str, name: &str, choice: Value) -> Value {
    json!({"attachment_hash":hash,"name":name,"decisions":[{"key":{"Site":{"position":5}},"choice":choice,"reason":"Reviewed the exact source call and retained native score evidence."}]})
}
fn assert_originals(before: &Workspace, after: &Workspace) {
    for original in &before.documents {
        let actual = after
            .documents
            .iter()
            .find(|d| d.id == original.id)
            .unwrap();
        assert_eq!(actual.head, original.head);
        assert_eq!(actual.head_event, original.head_event);
        assert_eq!(
            after.revisions[&actual.head],
            before.revisions[&original.head]
        );
    }
}
#[tokio::test]
async fn record_only_decisions_and_observed_base_revision_preserve_sources_and_exact_lineage() {
    let (c, server) = setup().await;
    assert_eq!(
        c.list_tools(Default::default()).await.unwrap().tools.len(),
        collagene_mcp::CURRENT_TOOL_COUNT
    );
    let (w, hash) = fixture(true);
    let original_bytes = workspace::save_workspace(&w).unwrap();
    let imported = upload(&c, &w).await;
    let original_wr = workspace_ref(&imported);
    // Retaining or leaving a call unresolved records a review but synthesizes no target.
    let keep = request(&hash, "Keep original reference", json!("RetainReference"));
    let preview = call(
        &c,
        "preview_read_review",
        json!({"workspace":original_wr,"request":keep}),
    )
    .await;
    assert_eq!(preview["state"], "proposal");
    assert_eq!(preview["verified"], false);
    assert!(preview["preview"]["target"].is_null());
    let recorded = call(
        &c,
        "commit_read_review",
        json!({"workspace":original_wr,"request":keep}),
    )
    .await;
    assert_eq!(recorded["receipt"]["state"], "historical");
    assert_eq!(recorded["receipt"]["verified"], false);
    assert!(recorded["receipt"]["target"].is_null());
    let kept =
        workspace::load_workspace(&native(&c, &recorded["native_artifact_id"]).await).unwrap();
    assert_originals(&w, &kept);
    assert_eq!(kept.documents.len(), w.documents.len());
    assert_eq!(kept.history.len(), w.history.len());
    let kept_wr = workspace_ref(&recorded);
    let unresolved = request(&hash, "Unresolved source call", json!("Unresolved"));
    let uncertain = call(
        &c,
        "commit_read_review",
        json!({"workspace":kept_wr,"request":unresolved}),
    )
    .await;
    assert!(uncertain["receipt"]["target"].is_null());
    let use_base = request(
        &hash,
        "Observed revised reference",
        json!({"UseBase":{"read_index":0,"call_index":5,"base":"T"}}),
    );
    let proposed = call(
        &c,
        "preview_read_review",
        json!({"workspace":workspace_ref(&uncertain),"request":use_base}),
    )
    .await;
    assert_eq!(proposed["verified"], false);
    let full = artifact(&c, &proposed["artifact_id"]).await;
    assert_eq!(full["target"]["sequence"]["bases"], "TTTCGTTACGGG");
    assert_eq!(full["edits"][0]["start"], 5);
    assert_eq!(full["edits"][0]["end"], 6);
    assert_eq!(full["edits"][0]["replacement"], "T");
    let committed = call(
        &c,
        "commit_read_review",
        json!({"workspace":workspace_ref(&uncertain),"request":use_base}),
    )
    .await;
    let target: EventBinding =
        serde_json::from_value(committed["receipt"]["target"].clone()).unwrap();
    let revised =
        workspace::load_workspace(&native(&c, &committed["native_artifact_id"]).await).unwrap();
    assert_originals(&w, &revised);
    assert_eq!(revised.documents.len(), w.documents.len() + 1);
    assert_eq!(
        revised.revisions[&target.revision].sequence.as_str(),
        "TTTCGTTACGGG"
    );
    let created = revised
        .history
        .iter()
        .find(|e| e.id == target.event_id)
        .unwrap();
    assert!(created.parents.contains(&event(&w, 0).event_id));
    assert!(created.parents.contains(&event(&w, 1).event_id));
    let saved = artifact(&c, &committed["artifact_id"]).await;
    assert_eq!(saved["reference"], binding(&w, 0));
    assert_eq!(saved["reads"], json!([binding(&w, 1)]));
    assert_eq!(saved["target"], committed["receipt"]["target"]);
    assert_eq!(saved["request"], use_base);
    let wr = workspace_ref(&committed);
    let retry = call(
        &c,
        "commit_read_review",
        json!({"workspace":wr,"request":use_base}),
    )
    .await;
    assert_eq!(workspace_ref(&retry), wr);
    assert_eq!(
        retry["receipt"]["attachment_hash"],
        committed["receipt"]["attachment_hash"]
    );
    assert_eq!(retry["artifact_id"], committed["artifact_id"]);
    assert_eq!(retry["native_artifact_id"], committed["native_artifact_id"]);
    let listed = call(&c, "list_read_reviews", json!({"workspace":wr})).await;
    assert_eq!(listed["state"], "historical");
    assert_eq!(listed["reviews"].as_array().unwrap().len(), 3);
    let original_export = call(&c, "export_workspace_archive", original_wr.clone()).await;
    assert_eq!(
        native(&c, &original_export["artifact_id"]).await,
        original_bytes
    );
    assert!(
        call(&c, "list_read_reviews", json!({"workspace":original_wr})).await["reviews"]
            .as_array()
            .unwrap()
            .is_empty()
    );
    // A revised design is not silently attached to an old reference read check.
    let still_saved = read_checks::load(&revised, &hash).unwrap();
    assert_eq!(still_saved.reference, event(&w, 0));
    assert_ne!(still_saved.reference, target);
    drop(c);
    server.abort();
}
#[tokio::test]
async fn mismatched_evidence_keys_calls_ineligible_quality_and_stale_workspaces_publish_nothing() {
    let (c, server) = setup().await;
    let (w, hash) = fixture(true);
    let imported = upload(&c, &w).await;
    let wr = workspace_ref(&imported);
    let before = call(&c, "list_workspaces", json!({})).await["workspaces"]
        .as_array()
        .unwrap()
        .len();
    let valid = request(
        &hash,
        "guarded",
        json!({"UseBase":{"read_index":0,"call_index":5,"base":"T"}}),
    );
    let mut wrong_key = valid.clone();
    wrong_key["decisions"][0]["key"] = json!({"Site":{"position":usize::MAX}});
    let mut wrong_call = valid.clone();
    wrong_call["decisions"][0]["choice"] =
        json!({"UseBase":{"read_index":0,"call_index":0,"base":"T"}});
    let mut wrong_allele = valid.clone();
    wrong_allele["decisions"][0]["choice"]["UseBase"]["base"] = json!("A");
    let mut missing = valid.clone();
    missing["attachment_hash"] = json!("0".repeat(64));
    for request in [wrong_key, wrong_call, wrong_allele, missing] {
        reject(
            &c,
            "preview_read_review",
            json!({"workspace":wr,"request":request}),
        )
        .await;
        reject(
            &c,
            "commit_read_review",
            json!({"workspace":wr,"request":request}),
        )
        .await;
    }
    let mut stale = wr.clone();
    stale["revision"] = json!("0".repeat(64));
    reject(
        &c,
        "commit_read_review",
        json!({"workspace":stale,"request":valid}),
    )
    .await;
    assert_eq!(
        call(&c, "list_workspaces", json!({})).await["workspaces"]
            .as_array()
            .unwrap()
            .len(),
        before
    );
    let (low, low_hash) = fixture(false);
    let low_import = upload(&c, &low).await;
    let low_wr = workspace_ref(&low_import);
    let low_request = request(
        &low_hash,
        "cannot promote unscored evidence",
        json!({"UseBase":{"read_index":0,"call_index":5,"base":"T"}}),
    );
    let count = call(&c, "list_workspaces", json!({})).await["workspaces"]
        .as_array()
        .unwrap()
        .len();
    reject(
        &c,
        "commit_read_review",
        json!({"workspace":low_wr,"request":low_request}),
    )
    .await;
    assert_eq!(
        call(&c, "list_workspaces", json!({})).await["workspaces"]
            .as_array()
            .unwrap()
            .len(),
        count
    );
    drop(c);
    server.abort();
}
