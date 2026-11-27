//! Optional mismatch policy selects a distinct native inspector without changing exact-PCR artifacts.
use base64::{Engine, engine::general_purpose::STANDARD};
use collagene_core::{
    Molecule, Topology,
    analysis::Primer,
    catalog_digest, cloning_outputs,
    mismatch_pcr_digest_inspection::Preview,
    mutagenesis::PcrOptions,
    restriction_inspection::InspectionOptions,
    run_lineage::{EventBinding, InputEventBinding},
    run_material::{self, Selection},
    verification::{DesignPlan, Operation, PlanInput, Source},
    workspace::{self, Workspace},
};
use collagene_mcp::CollageneServer;
use rmcp::{RoleClient, ServiceExt, model::CallToolRequestParams, service::RunningService};
use serde_json::{Value, json};
type Client = RunningService<RoleClient, ()>;
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
    match c
        .call_tool(
            CallToolRequestParams::new(n.to_owned()).with_arguments(a.as_object().unwrap().clone()),
        )
        .await
    {
        Err(_) => {}
        Ok(r) => assert_eq!(r.is_error, Some(true), "{n}: {:?}", r),
    }
}
async fn evidence(c: &Client, id: &Value) -> Value {
    let (mut at, mut text) = (0, String::new());
    loop {
        let p = call(
            c,
            "read_artifact",
            json!({"artifact_id":id,"start":at,"max_bytes":16384}),
        )
        .await;
        text.push_str(p["text"].as_str().unwrap());
        if p["done"] == true {
            break;
        }
        at = p["next_start"].as_u64().unwrap();
    }
    serde_json::from_str(&text).unwrap()
}
fn reference(import: &Value, index: usize) -> Value {
    let d = &import["documents"][index];
    json!({"id":d["id"],"revision":d["revision"]})
}
fn event(w: &Workspace, i: usize) -> EventBinding {
    let d = &w.documents[i];
    EventBinding {
        document_id: d.id.clone(),
        event_id: d.head_event.clone().unwrap(),
        revision: d.head.clone(),
    }
}
fn saved() -> cloning_outputs::Commit {
    let source = Molecule::dna("Original", "ACGACCCCGGGGTAAC", Topology::Linear).unwrap();
    let target = Molecule::dna(
        "Independent exact requirement",
        "ACGACCCCGGGGTAAC",
        Topology::Linear,
    )
    .unwrap();
    let plan = DesignPlan {
        name: "Original exact PCR".into(),
        target_revision: target.revision(),
        inputs: vec![PlanInput {
            id: "source".into(),
            revision: source.revision(),
            molecule: source.clone(),
        }],
        steps: vec![Operation::Pcr {
            template: Source::Input("source".into()),
            forward: Primer::new("F", "ACGACCCC", "").unwrap(),
            reverse: Primer::new("R", "GTTACCCC", "").unwrap(),
            template_strand: collagene_core::chemistry::Strand::Top,
        }],
    };
    let w = Workspace::from_documents("Actual saved PCR", vec![source, target], 1).unwrap();
    cloning_outputs::finalize(
        &w,
        &plan,
        event(&w, 1),
        vec![InputEventBinding {
            input_id: "source".into(),
            event: event(&w, 0),
        }],
        &cloning_outputs::Options {
            name: "Saved PCR".into(),
            show_prepared: false,
            save_primers: false,
            allow_incomplete: false,
        },
    )
    .unwrap()
}
#[tokio::test]
async fn explicit_mismatch_digest_policy_actual_fragments_incomplete_search_and_saved_origin_route()
{
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
    assert_eq!(
        c.list_tools(Default::default()).await.unwrap().tools.len(),
        collagene_mcp::CURRENT_TOOL_COUNT
    );
    let imported=call(&c,"import_document",json!({"text":">original\nACGACCCCGGGGTAAC\n>repeated\nACGACCCCGGGGTAACACGACCCCGGGGTAAC\n"})).await;
    let catalog = catalog_digest::builtin_catalog().unwrap();
    let id = catalog
        .edition()
        .enzymes
        .iter()
        .find(|e| e.name == "EcoRI")
        .unwrap()
        .id
        .clone();
    let pin = catalog_digest::select(&catalog, &[id]).unwrap();
    let mut policy = PcrOptions::default();
    policy.binding.min_three_prime_exact = 2;
    let request = json!({"template_strand":"Top","source_role":"original","forward":Primer::new("Mismatch F","ACTACCCC","CCGAATTC").unwrap(),"reverse":Primer::new("R","GTTACCCC","GGGAATTC").unwrap(),"selection":pin,"options":InspectionOptions::default()});
    let args = json!({"document":reference(&imported,0),"request":request,"pcr_options":policy});
    let result = call(&c, "inspect_pcr_digest", args.clone()).await;
    assert_eq!(result["mode"], "mismatch-pcr");
    assert_eq!(result["incomplete"], false);
    assert_eq!(result["fragment_count"], 3);
    let preview: Preview =
        serde_json::from_value(evidence(&c, &result["artifact_id"]).await).unwrap();
    preview
        .validate_binding(result["context_revision"].as_str().unwrap())
        .unwrap();
    assert_eq!(
        preview.amplicon.as_ref().unwrap().sequence.as_str(),
        "CCGAATTCACTACCCCGGGGTAACGAATTCCC"
    );
    assert_eq!(preview.request.pcr_options.binding.min_three_prime_exact, 2);
    assert_eq!(
        preview
            .inspection
            .as_ref()
            .unwrap()
            .fragments
            .iter()
            .map(|f| f.duplex.top_length)
            .collect::<Vec<_>>(),
        vec![3, 22, 7]
    );
    assert!(preview.origin_revision.is_none());
    let mut bounded = args.clone();
    bounded["document"] = reference(&imported, 1);
    bounded["pcr_options"]["max_products"] = json!(1);
    let partial = call(&c, "inspect_pcr_digest", bounded).await;
    assert_eq!(partial["incomplete"], true);
    let partial: Preview =
        serde_json::from_value(evidence(&c, &partial["artifact_id"]).await).unwrap();
    assert!(!partial.pcr.report.enumeration_complete);
    assert!(partial.incomplete);
    let mut stale = args.clone();
    stale["document"]["revision"] = json!("0".repeat(64));
    reject(&c, "inspect_pcr_digest", stale).await;
    let mut invalid = args.clone();
    invalid["pcr_options"]["max_pair_checks"] = json!(0);
    reject(&c, "inspect_pcr_digest", invalid).await;
    let commit = saved();
    let selection = Selection {
        run_id: commit.receipt.run.run_id.clone(),
        step_index: 0,
        event: commit.receipt.construct.clone(),
    };
    let origin = run_material::capsule(&commit.workspace, &selection).unwrap();
    let uploaded=call(&c,"import_workspace_archive",json!({"data_base64":STANDARD.encode(workspace::save_workspace(&commit.workspace).unwrap())})).await;
    let wr = json!({"id":uploaded["workspace"]["id"],"revision":uploaded["workspace"]["revision"]});
    let before = call(&c, "list_documents", json!({})).await["total"].clone();
    let mut from_run = args;
    from_run["document"] =
        json!({"id":selection.event.revision,"revision":selection.event.revision});
    from_run["origin"] = serde_json::to_value(collagene_core::run_material::MaterialOrigin::Run(
        origin.clone(),
    ))
    .unwrap();
    reject(&c, "inspect_pcr_digest", from_run.clone()).await;
    from_run["workspace"] = wr;
    let inspected = call(&c, "inspect_pcr_digest", from_run.clone()).await;
    let p: Preview = serde_json::from_value(evidence(&c, &inspected["artifact_id"]).await).unwrap();
    p.validate_binding(inspected["context_revision"].as_str().unwrap())
        .unwrap();
    assert!(p.origin_revision.is_some());
    assert_eq!(
        p.amplicon.as_ref().unwrap().sequence.as_str(),
        "CCGAATTCACTACCCCGGGGTAACGAATTCCC"
    );
    let mut forged = from_run;
    forged["origin"]["Run"]["selection"]["event"]["document_id"] = json!("forged-owner");
    reject(&c, "inspect_pcr_digest", forged).await;
    assert_eq!(call(&c, "list_documents", json!({})).await["total"], before);
    assert_eq!(
        call(
            &c,
            "read_sequence_region",
            json!({"document":reference(&imported,0),"start":0,"end":16})
        )
        .await["sequence"],
        "ACGACCCCGGGGTAAC"
    );
    drop(c);
    server.abort();
}
