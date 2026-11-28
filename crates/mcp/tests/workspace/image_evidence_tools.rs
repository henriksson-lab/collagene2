//! Actual in-memory MCP SDK calls; no HTTP listener, browser or filesystem tool.
use base64::{Engine, engine::general_purpose::STANDARD};
use collagene_core::{
    documents::{
        history::DocumentHistory,
        image_evidence::{self, Binding},
        workspace::{self, Workspace, WorkspaceAttachment},
    },
    material::lineage::workspace_fingerprint,
};
use collagene_mcp::server::CollageneServer;
use collagene_model::{Molecule, Topology};
use rmcp::{RoleClient, ServiceExt, model::CallToolRequestParams, service::RunningService};
use serde_json::{Value, json};
type Client = RunningService<RoleClient, ()>;
const PNG: &[u8] = include_bytes!("../data/image-evidence/known_2x3.png");
const OTHER_PNG: &[u8] = include_bytes!("../data/image-evidence/other_2x3.png");
const TIFF: &[u8] = include_bytes!("../data/image-evidence/known_2x3_orientation6.tif");
const JPEG: &[u8] = include_bytes!("../data/image-evidence/known_8x8.jpg");

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
        "{name} unexpectedly succeeded"
    );
}
async fn binary(c: &Client, id: &Value) -> Vec<u8> {
    let mut text = String::new();
    let mut at = 0;
    loop {
        let page = call(
            c,
            "read_artifact",
            json!({"artifact_id":id,"start":at,"max_bytes":4096}),
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
fn reference(reply: &Value) -> Value {
    json!({"id":reply["workspace"]["id"],"revision":reply["workspace"]["revision"]})
}
fn selection(frame: &Workspace, index: usize) -> Binding {
    Binding {
        document_id: frame.documents[index].id.clone(),
        revision: frame.documents[index].head.clone(),
    }
}
async fn import(c: &Client, frame: &Workspace) -> Value {
    let bytes = workspace::save_workspace(frame).unwrap();
    let reply = call(
        c,
        "import_workspace_archive",
        json!({"data_base64":STANDARD.encode(bytes)}),
    )
    .await;
    reference(&reply)
}
async fn restore(c: &Client, reference: &Value) -> Workspace {
    let reply = call(c, "export_workspace_archive", reference.clone()).await;
    workspace::load_workspace(&binary(c, &reply["artifact_id"]).await).unwrap()
}
fn fixture() -> Workspace {
    // Equal molecules do not make the two document identities interchangeable.
    let source = Molecule::dna("reference", "ACGTTGCA", Topology::Linear).unwrap();
    let mut frame =
        Workspace::from_documents("Image fixture", vec![source.clone(), source], 0).unwrap();
    let note = b"Preserve this unrelated documentary source".to_vec();
    frame.attachments.insert(
        blake3::hash(&note).to_hex().to_string(),
        WorkspaceAttachment {
            kind: "fixture-note".into(),
            name: "retained-note.txt".into(),
            media_type: "text/plain".into(),
            bytes: note,
        },
    );
    frame.validate().unwrap();
    frame
}
fn chemistry_and_unrelated_attachments(frame: &Workspace) -> Value {
    let mut frame = frame.clone();
    frame.attachments.retain(|_, a| {
        a.kind != image_evidence::MANIFEST_KIND && a.kind != image_evidence::IMAGE_KIND
    });
    serde_json::to_value(frame).unwrap()
}

#[tokio::test]
async fn images_roundtrip_original_and_fresh_preview_without_mutating_source_or_chemistry() {
    let (c, server) = setup().await;
    let tools = c.list_tools(Default::default()).await.unwrap().tools;
    assert_eq!(tools.len(), collagene_mcp::CURRENT_TOOL_COUNT);
    for name in [
        "attach_image_evidence",
        "list_image_evidence",
        "export_image_evidence",
    ] {
        assert_eq!(tools.iter().filter(|t| t.name == name).count(), 1);
    }
    let original = fixture();
    let selected = selection(&original, 0);
    let source = import(&c, &original).await;
    let mut current = source.clone();
    let mut last_preview = vec![];
    for (name, bytes, format, width, height, orientation, pw, ph) in [
        ("Observed PNG", PNG, "Png", 2, 3, 1, 2, 3),
        ("Observed TIFF", TIFF, "Tiff", 2, 3, 6, 3, 2),
        ("Observed JPEG", JPEG, "Jpeg", 8, 8, 1, 8, 8),
    ] {
        let listed = call(
            &c,
            "list_image_evidence",
            json!({"workspace":current,"document":selected,"limit":1}),
        )
        .await;
        let before = restore(&c, &current).await;
        assert_eq!(
            listed["source_workspace_fingerprint"],
            workspace_fingerprint(&before).unwrap()
        );
        let args = json!({"workspace":current,"document":selected,"expected_source_fingerprint":listed["source_workspace_fingerprint"],"name":name,"data_base64":STANDARD.encode(bytes)});
        let attached = call(&c, "attach_image_evidence", args.clone()).await;
        assert_eq!(attached["receipt"]["state"], "documentary");
        assert_eq!(attached["receipt"]["verified"], false);
        assert_eq!(attached["receipt"]["source_workspace"], current);
        let record = &attached["receipt"]["image"]["record"];
        assert_eq!(record["format"], format);
        assert_eq!(record["original"]["width"], width);
        assert_eq!(record["original"]["height"], height);
        assert_eq!(record["source_orientation"], orientation);
        assert_eq!(record["preview"]["width"], pw);
        assert_eq!(record["preview"]["height"], ph);
        assert_eq!(
            record["original"]["hash"],
            blake3::hash(bytes).to_hex().to_string()
        );
        let repeated = call(&c, "attach_image_evidence", args).await;
        assert_eq!(repeated["workspace"], attached["workspace"]);
        assert_eq!(repeated["deduplicated"], true);
        current = reference(&attached);
        let manifest = &attached["receipt"]["image"]["manifest_hash"];
        let original_export=call(&c,"export_image_evidence",json!({"workspace":current,"document":selected,"manifest_hash":manifest,"representation":"Original"})).await;
        assert_eq!(original_export["freshly_redecoded"], false);
        assert_eq!(binary(&c, &original_export["artifact_id"]).await, bytes);
        let preview=call(&c,"export_image_evidence",json!({"workspace":current,"document":selected,"manifest_hash":manifest,"representation":"Preview"})).await;
        assert_eq!(preview["freshly_redecoded"], true);
        assert_eq!(preview["verified"], false);
        assert_eq!(preview["calibrated"], false);
        last_preview = binary(&c, &preview["artifact_id"]).await;
        assert!(last_preview.starts_with(b"\x89PNG\r\n\x1a\n"));
        assert_eq!(
            u32::from_be_bytes(last_preview[16..20].try_into().unwrap()),
            pw as u32
        );
        assert_eq!(
            u32::from_be_bytes(last_preview[20..24].try_into().unwrap()),
            ph as u32
        );
        let next = restore(&c, &current).await;
        assert_eq!(
            chemistry_and_unrelated_attachments(&next),
            serde_json::to_value(&original).unwrap()
        );
        assert_eq!(
            serde_json::to_value(restore(&c, &source).await).unwrap(),
            serde_json::to_value(&original).unwrap()
        );
    }
    // A canonical PNG may be both original and preview, with one byte artifact.
    let listed = call(
        &c,
        "list_image_evidence",
        json!({"workspace":current,"document":selected}),
    )
    .await;
    let added=call(&c,"attach_image_evidence",json!({"workspace":current,"document":selected,"expected_source_fingerprint":listed["source_workspace_fingerprint"],"name":"Canonical source PNG","data_base64":STANDARD.encode(&last_preview)})).await;
    let next = reference(&added);
    let manifest = &added["receipt"]["image"]["manifest_hash"];
    let a=call(&c,"export_image_evidence",json!({"workspace":next,"document":selected,"manifest_hash":manifest,"representation":"Original"})).await;
    let b=call(&c,"export_image_evidence",json!({"workspace":next,"document":selected,"manifest_hash":manifest,"representation":"Preview"})).await;
    assert_eq!(a["artifact_id"], b["artifact_id"]);
    assert_eq!(binary(&c, &a["artifact_id"]).await, last_preview);
    c.cancel().await.unwrap();
    server.await.unwrap();
}

#[tokio::test]
async fn image_cas_and_association_refusals_publish_nothing_and_old_revisions_remain_discoverable()
{
    let (c, server) = setup().await;
    let frame = fixture();
    let a = selection(&frame, 0);
    let b = selection(&frame, 1);
    let source = import(&c, &frame).await;
    let listed = call(
        &c,
        "list_image_evidence",
        json!({"workspace":source,"document":a}),
    )
    .await;
    let base = json!({"workspace":source,"document":a,"expected_source_fingerprint":listed["source_workspace_fingerprint"],"name":"Documentary","data_base64":STANDARD.encode(PNG)});
    let before = call(&c, "list_workspaces", json!({})).await;
    let mut wrong = base.clone();
    wrong["expected_source_fingerprint"] = json!("0".repeat(64));
    reject(&c, "attach_image_evidence", wrong).await;
    let mut malformed = base.clone();
    malformed["data_base64"] = json!("not-base64!");
    reject(&c, "attach_image_evidence", malformed).await;
    let mut unknown = base.clone();
    unknown["document"]["document_id"] = json!("document-missing");
    reject(&c, "attach_image_evidence", unknown).await;
    let mut stale = base.clone();
    stale["workspace"]["revision"] = json!("f".repeat(64));
    reject(&c, "attach_image_evidence", stale).await;
    assert_eq!(call(&c, "list_workspaces", json!({})).await, before);
    let attached = call(&c, "attach_image_evidence", base).await;
    let next = reference(&attached);
    let manifest = &attached["receipt"]["image"]["manifest_hash"];
    assert_eq!(
        call(
            &c,
            "list_image_evidence",
            json!({"workspace":next,"document":b})
        )
        .await["total"],
        0
    );
    reject(
        &c,
        "export_image_evidence",
        json!({"workspace":next,"document":b,"manifest_hash":manifest,"representation":"Original"}),
    )
    .await;
    reject(
        &c,
        "list_image_evidence",
        json!({"workspace":next,"document":a,"limit":0}),
    )
    .await;
    reject(
        &c,
        "list_image_evidence",
        json!({"workspace":next,"document":a,"start":2}),
    )
    .await;
    let mut changed = restore(&c, &next).await;
    let mut histories = DocumentHistory::restore_all(&changed).unwrap();
    histories[0]
        .commit(
            Molecule::dna("edited", "ACGTTGCC", Topology::Linear).unwrap(),
            "Sequence edit",
        )
        .unwrap();
    changed.documents.clear();
    changed.revisions.clear();
    changed.history.clear();
    for history in histories {
        history.contribute(&mut changed).unwrap();
    }
    changed.validate().unwrap();
    let current = selection(&changed, 0);
    let changed_ref = import(&c, &changed).await;
    let history_list = call(
        &c,
        "list_image_evidence",
        json!({"workspace":changed_ref,"document":current}),
    )
    .await;
    assert_eq!(history_list["total"], 1);
    assert_eq!(
        history_list["images"][0]["record"]["association"],
        serde_json::to_value(&a).unwrap()
    );
    reject(&c,"export_image_evidence",json!({"workspace":changed_ref,"document":current,"manifest_hash":manifest,"representation":"Preview"})).await;
    let exported=call(&c,"export_image_evidence",json!({"workspace":changed_ref,"document":a,"manifest_hash":manifest,"representation":"Original"})).await;
    assert_eq!(binary(&c, &exported["artifact_id"]).await, PNG);
    reject(&c,"attach_image_evidence",json!({"workspace":changed_ref,"document":a,"expected_source_fingerprint":history_list["source_workspace_fingerprint"],"name":"Stale head","data_base64":STANDARD.encode(PNG)})).await;
    c.cancel().await.unwrap();
    server.await.unwrap();
}

#[tokio::test]
async fn archived_image_preview_claims_need_fresh_decode_and_another_session_has_no_access() {
    let (c, server) = setup().await;
    let frame = fixture();
    let document = selection(&frame, 0);
    let source = import(&c, &frame).await;
    let listed = call(
        &c,
        "list_image_evidence",
        json!({"workspace":source,"document":document}),
    )
    .await;
    let attached=call(&c,"attach_image_evidence",json!({"workspace":source,"document":document,"expected_source_fingerprint":listed["source_workspace_fingerprint"],"name":"Original pixels","data_base64":STANDARD.encode(PNG)})).await;
    let next = reference(&attached);
    let original_manifest = attached["receipt"]["image"]["manifest_hash"]
        .as_str()
        .unwrap();
    let mut forged = restore(&c, &next).await;
    let mut record: image_evidence::Record =
        serde_json::from_slice(&forged.attachments[original_manifest].bytes).unwrap();
    // A valid but unrelated PNG with the same geometry is a metadata/hash-valid
    // archive preview. Import does not decode it or infer an image check.
    let unrelated=call(&c,"attach_image_evidence",json!({"workspace":source,"document":document,"expected_source_fingerprint":listed["source_workspace_fingerprint"],"name":"Other pixels","data_base64":STANDARD.encode(OTHER_PNG)})).await;
    let other=call(&c,"export_image_evidence",json!({"workspace":reference(&unrelated),"document":document,"manifest_hash":unrelated["receipt"]["image"]["manifest_hash"],"representation":"Preview"})).await;
    let replacement = binary(&c, &other["artifact_id"]).await;
    let old_preview = record.preview.hash.clone();
    let hash = blake3::hash(&replacement).to_hex().to_string();
    record.preview.hash = hash.clone();
    record.preview.bytes = replacement.len() as u64;
    record.preview.width = 2;
    record.preview.height = 3;
    if old_preview != record.original.hash {
        forged.attachments.remove(&old_preview);
    }
    forged.attachments.insert(
        hash.clone(),
        WorkspaceAttachment {
            kind: image_evidence::IMAGE_KIND.into(),
            name: format!("image-{hash}.png"),
            media_type: "image/png".into(),
            bytes: replacement,
        },
    );
    forged.attachments.remove(original_manifest);
    let manifest = serde_json::to_vec(&record).unwrap();
    let manifest_hash = blake3::hash(&manifest).to_hex().to_string();
    forged.attachments.insert(
        manifest_hash.clone(),
        WorkspaceAttachment {
            kind: image_evidence::MANIFEST_KIND.into(),
            name: format!("image-evidence-{manifest_hash}.json"),
            media_type: "application/json".into(),
            bytes: manifest,
        },
    );
    forged.validate().unwrap();
    let forged_ref = import(&c, &forged).await;
    reject(&c,"export_image_evidence",json!({"workspace":forged_ref,"document":document,"manifest_hash":manifest_hash,"representation":"Preview"})).await;
    let exact=call(&c,"export_image_evidence",json!({"workspace":forged_ref,"document":document,"manifest_hash":manifest_hash,"representation":"Original"})).await;
    assert_eq!(binary(&c, &exact["artifact_id"]).await, PNG);
    assert_eq!(exact["freshly_redecoded"], false);
    let (other_client, other_server) = setup().await;
    reject(
        &other_client,
        "list_image_evidence",
        json!({"workspace":next,"document":document}),
    )
    .await;
    other_client.cancel().await.unwrap();
    other_server.await.unwrap();
    c.cancel().await.unwrap();
    server.await.unwrap();
}
