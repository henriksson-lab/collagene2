//! Actual SDK routing; no transport mocks, files, network, browser or scientific status adoption.
use base64::{Engine, engine::general_purpose::STANDARD};
use collagene_core::documents::reference_edit::Operation as Edit;
use collagene_core::{
    documents::{
        history::DocumentHistory,
        transfer,
        workspace::{self, Workspace},
    },
    material::lineage::{self as run_lineage, EventBinding, InputEventBinding},
    planning::plan::{DesignPlan, Operation, PlanInput, Source},
    sequence::analysis::Primer,
};
use collagene_mcp::server::CollageneServer;
use collagene_model::{Feature, Location, Molecule, Topology};
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
    assert!(
        r.is_err() || r.unwrap().is_error == Some(true),
        "must reject {name}"
    )
}
async fn artifact(c: &Client, id: &Value) -> (Vec<u8>, Value) {
    let mut text = String::new();
    let mut at = 0;
    loop {
        let r = call(
            c,
            "read_artifact",
            json!({"artifact_id":id,"start":at,"max_bytes":16384}),
        )
        .await;
        text.push_str(r["text"].as_str().unwrap());
        if r["done"] == true {
            return (
                if r["encoding"] == "base64" {
                    STANDARD.decode(text).unwrap()
                } else {
                    text.into_bytes()
                },
                r,
            );
        }
        at = r["next_start"].as_u64().unwrap();
    }
}
fn event(w: &Workspace, i: usize) -> EventBinding {
    let d = &w.documents[i];
    EventBinding {
        document_id: d.id.clone(),
        event_id: d.head_event.clone().unwrap(),
        revision: d.head.clone(),
    }
}
fn fixture() -> (Workspace, String, EventBinding, EventBinding, Molecule) {
    let mut source = Molecule::dna(
        "Original <template> & source",
        "ACGACCCCGGGGTAAC",
        Topology::Linear,
    )
    .unwrap();
    source.features.push(Feature {
        kind: "misc_feature".into(),
        location: Location::range(0, 8),
        qualifiers: vec![("label".into(), Some("Reporter <&> quoted".into()))],
        needs_review: false,
    });
    let target = Molecule::dna(
        "Original intended independent target",
        "ACGACCCCGGGGTAAC",
        Topology::Linear,
    )
    .unwrap();
    let original = Workspace::from_documents(
        "Historical exact event review",
        vec![source.clone(), target.clone()],
        1,
    )
    .unwrap();
    let input_event = event(&original, 0);
    let target_event = event(&original, 1);
    let plan = DesignPlan {
        name: "Retained literal PCR".into(),
        target_revision: target.revision(),
        inputs: vec![PlanInput {
            id: "template".into(),
            revision: source.revision(),
            molecule: source,
        }],
        steps: vec![Operation::Pcr {
            template: Source::Input("template".into()),
            forward: Primer::new("F", "ACGA", "").unwrap(),
            reverse: Primer::new("R", "GTTA", "").unwrap(),
            template_strand: collagene_core::chemistry::Strand::Top,
        }],
    };
    let saved = run_lineage::record_run(
        &original,
        &plan,
        target_event.clone(),
        vec![InputEventBinding {
            input_id: "template".into(),
            event: input_event.clone(),
        }],
        false,
    )
    .unwrap();
    let run = saved.receipt.run_id.clone();
    let mut frame = saved.workspace;
    let mut histories = DocumentHistory::restore_all(&frame).unwrap();
    frame.documents.clear();
    frame.history.clear();
    let mut current_source = None;
    for h in &mut histories {
        if h.id == input_event.document_id {
            let edit_end = h.molecule.sequence.len();
            h.commit(
                vec![
                    Edit::ReplaceSequence {
                        start: 0,
                        end: edit_end,
                        replacement: "TTTTCCCCGGGGAAAA".into(),
                    },
                    Edit::Rename {
                        name: "Current <source> & display".into(),
                    },
                ],
                "Private source edit after run",
            )
            .unwrap();
            current_source = Some(h.molecule.clone());
        } else if h.id == target_event.document_id {
            let edit_end = h.molecule.sequence.len();
            h.commit(
                vec![
                    Edit::ReplaceSequence {
                        start: 0,
                        end: edit_end,
                        replacement: "TTTTTTTTTTTTTTTT".into(),
                    },
                    Edit::Rename {
                        name: "Current independent target changed".into(),
                    },
                ],
                "Target edit after historical run",
            )
            .unwrap();
        }
        h.contribute(&mut frame).unwrap();
    }
    frame.validate().unwrap();
    (
        frame,
        run,
        input_event,
        target_event,
        current_source.unwrap(),
    )
}
#[tokio::test]
async fn retained_run_review_is_exact_unverified_and_map_export_is_bound_escaped_and_side_effect_free()
 {
    let (frame, run, input_event, target_event, source) = fixture();
    assert_ne!(
        frame.head(&target_event.document_id).unwrap().revision(),
        target_event.revision
    );
    assert_ne!(
        frame.head(&input_event.document_id).unwrap().revision(),
        input_event.revision
    );
    let original_bytes = workspace::save_workspace(&frame).unwrap();
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
    assert!(tools.iter().any(|t| t.name == "prepare_saved_run_review"));
    assert!(tools.iter().any(|t| t.name == "export_map"));
    let imported = call(
        &c,
        "import_workspace_archive",
        json!({"data_base64":STANDARD.encode(&original_bytes)}),
    )
    .await;
    let wr = json!({"id":imported["workspace"]["id"],"revision":imported["workspace"]["revision"]});
    let documents = call(&c, "list_documents", json!({"start":0,"limit":128})).await;
    let args = json!({"workspace":wr,"run_id":run});
    let review = call(&c, "prepare_saved_run_review", args.clone()).await;
    assert_eq!(review["state"], "candidate");
    assert_eq!(review["verified"], false);
    assert_eq!(review["target_event"], json!(target_event));
    assert_eq!(review["inputs"][0]["event"], json!(input_event));
    assert_eq!(review["target"]["revision"], target_event.revision);
    assert_ne!(
        review["target"]["revision"],
        frame.head(&target_event.document_id).unwrap().revision()
    );
    let (raw, meta) = artifact(&c, &review["artifact_id"]).await;
    assert_eq!(meta["filename"], "saved-run-review.json");
    let prepared: collagene_core::material::review::PreparedReview =
        serde_json::from_slice(&raw).unwrap();
    prepared.validate_for_workspace(&frame).unwrap();
    assert_eq!(prepared.inputs[0].event, input_event);
    assert_eq!(prepared.target_event, target_event);
    let (bytes, meta) = artifact(&c, &review["candidate_artifact_id"]).await;
    assert_eq!(meta["encoding"], "base64");
    let candidate = transfer::load_candidate(&bytes).unwrap();
    assert_eq!(candidate.plan_revision, prepared.candidate.plan_revision);
    assert_eq!(candidate.plan.inputs[0].revision, input_event.revision);
    assert_eq!(candidate.plan.target_revision, target_event.revision);
    // Prepared session plan uses the original independent target; preparation itself never verified it.
    let checked = call(
        &c,
        "verify_plan",
        json!({"plan_revision":review["plan_revision"]}),
    )
    .await;
    assert_eq!(checked["target_revision"], target_event.revision);
    let svg = call(
        &c,
        "export_map",
        json!({"document":{"id":source.revision(),"revision":source.revision()},"format":"Svg"}),
    )
    .await;
    assert_eq!(svg["source_revision"], source.revision());
    assert_eq!(svg["restriction_overlay"], "none_requested");
    assert_eq!(svg["warnings_truncated"], false);
    let (svg_bytes, meta) = artifact(&c, &svg["artifact_id"]).await;
    assert_eq!(meta["mime_type"], "image/svg+xml");
    let expected = collagene_core::output::map_export::Request {
        source: source.clone(),
        source_revision: source.revision(),
        restriction_sites: None,
        format: collagene_core::output::raster::Format::Svg,
    }
    .render()
    .unwrap();
    assert_eq!(svg_bytes, expected.image.bytes());
    assert_eq!(svg["warnings"], json!(expected.warnings));
    let text = String::from_utf8(svg_bytes).unwrap();
    assert!(text.contains("&lt;"));
    assert!(text.contains("&amp;"));
    assert!(!text.contains("<source>"));
    assert!(!text.contains("<script"));
    let png = call(&c, "export_map", json!({"document":{"id":source.revision(),"revision":source.revision()},"format":{"Png":{"width":800}}})).await;
    let (png_bytes, png_meta) = artifact(&c, &png["artifact_id"]).await;
    assert_eq!(png_meta["mime_type"], "image/png");
    assert_eq!(png_meta["encoding"], "base64");
    assert_eq!(&png_bytes[..8], b"\x89PNG\r\n\x1a\n");
    assert_eq!(
        u32::from_be_bytes(png_bytes[16..20].try_into().unwrap()),
        800
    );
    assert!(
        png_bytes
            .windows(source.revision().len())
            .any(|s| s == source.revision().as_bytes())
    );
    reject(&c, "export_map", json!({"document":{"id":source.revision(),"revision":source.revision()},"format":{"Png":{"width":4097}}})).await;
    let again = call(&c, "prepare_saved_run_review", args.clone()).await;
    assert_eq!(again["artifact_id"], review["artifact_id"]);
    assert_eq!(
        again["candidate_artifact_id"],
        review["candidate_artifact_id"]
    );
    assert_eq!(
        call(&c, "list_documents", json!({"start":0,"limit":128})).await,
        documents,
        "only review artifacts/unverified plan added; no document installed"
    );
    let exported = call(&c, "export_workspace_archive", wr.clone()).await;
    assert_eq!(
        artifact(&c, &exported["artifact_id"]).await.0,
        original_bytes,
        "historical heads/private branches untouched"
    );
    let mut stale = args.clone();
    stale["workspace"]["revision"] = json!("0".repeat(64));
    reject(&c, "prepare_saved_run_review", stale).await;
    let mut missing = args.clone();
    missing["run_id"] = json!("f".repeat(64));
    reject(&c, "prepare_saved_run_review", missing).await;
    let mut extra = args;
    extra["adopt_status"] = json!(true);
    reject(&c, "prepare_saved_run_review", extra).await;
    reject(
        &c,
        "export_map",
        json!({"document":{"id":source.revision(),"revision":"0".repeat(64)},"format":"Svg"}),
    )
    .await;
    reject(
        &c,
        "export_map",
        json!({"document":{"id":source.revision(),"revision":source.revision()},"format":"Svg","restriction_sites":[]}),
    )
    .await;
    c.cancel().await.unwrap();
    server.await.unwrap();
}
