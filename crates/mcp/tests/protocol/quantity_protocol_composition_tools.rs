#[path = "../../../core/tests/common/protocol_quantities.rs"]
mod fixture;
use base64::{Engine as _, engine::general_purpose::STANDARD};
use collagene_core::{
    documents::{history::DocumentHistory, transfer::Candidate, workspace},
    material::lineage as run_lineage,
    output::quantities as q,
    verification,
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
    assert_ne!(r.is_error, Some(true), "{name}: {:?}", r.content);
    r.structured_content.unwrap()
}
async fn artifact(c: &Client, id: &Value) -> Vec<u8> {
    let mut text = String::new();
    let mut start = 0;
    let encoding;
    loop {
        let p = call(
            c,
            "read_artifact",
            json!({"artifact_id":id,"start":start,"max_bytes":16384}),
        )
        .await;
        text.push_str(p["text"].as_str().unwrap());
        if p["done"] == true {
            encoding = p["encoding"].as_str().unwrap().to_owned();
            break;
        }
        start = p["next_start"].as_u64().unwrap();
    }
    let bytes = if encoding == "base64" {
        STANDARD.decode(text).unwrap()
    } else {
        text.into_bytes()
    };
    assert_eq!(blake3::hash(&bytes).to_hex().as_str(), id.as_str().unwrap());
    bytes
}
#[tokio::test]
async fn combined_quantity_exports_have_fresh_molecular_receipts_not_arithmetic_acceptance() {
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
    let tools = c.list_tools(Default::default()).await.unwrap();
    assert_eq!(tools.tools.len(), collagene_mcp::CURRENT_TOOL_COUNT);
    let description = tools
        .tools
        .iter()
        .find(|t| t.name == "protocol_quantities")
        .unwrap()
        .description
        .as_ref()
        .unwrap();
    assert!(description.contains("ProtocolHtml") && description.contains("ProtocolPdf"));
    for (wrong_target, missing_quantity) in [(false, false), (true, false), (false, true)] {
        let (mut w, mut r) = fixture::fixture();
        if wrong_target {
            let run = match &r.plan_source {
                q::PlanSource::SavedRun { run_id } => run_lineage::load_run(&w, run_id).unwrap(),
                _ => unreachable!(),
            };
            let wrong = Molecule::dna(
                "Independent wrong protocol target",
                "AGTCCGATGCTACGTT",
                Topology::Linear,
            )
            .unwrap();
            DocumentHistory::new("wrong-protocol-target", wrong.clone())
                .unwrap()
                .contribute(&mut w)
                .unwrap();
            let mut plan = run.plan;
            plan.target_revision = wrong.revision();
            r.plan_source = q::PlanSource::Candidate;
            r.plan_revision = collagene_core::planning::plan::DesignPlan::revision(&plan);
            r.target = fixture::pin(&w, w.documents.len() - 1);
            w.candidates.push(Candidate::new(plan));
        }
        if missing_quantity {
            r.reactions[0].final_volume = None;
        }
        let captured = q::capture(&w, &r).unwrap();
        let mut w = captured.workspace;
        fixture::advance(&mut w, 0, "Later unrelated stock head");
        let before = workspace::save_workspace(&w).unwrap();
        let imported = call(
            &c,
            "import_workspace_archive",
            json!({"data_base64":STANDARD.encode(&before)}),
        )
        .await;
        let reference =
            json!({"id":imported["workspace"]["id"],"revision":imported["workspace"]["revision"]});
        for format in ["ProtocolHtml", "ProtocolPdf"] {
            let e=call(&c,"protocol_quantities",json!({"workspace":reference,"task":{"Export":{"quantity_hash":captured.attachment_hash,"format":format}}})).await;
            assert_eq!(e["verified"], false);
            assert_eq!(e["arithmetic_complete"], !missing_quantity);
            assert_eq!(
                e["molecular"]["status"],
                if wrong_target {
                    "Invalid"
                } else {
                    "ConsistentUnderAssumptions"
                }
            );
            assert_eq!(
                e["molecular"]["engine_version"],
                verification::ENGINE_VERSION
            );
            assert_eq!(e["plan_revision"], r.plan_revision);
            assert_eq!(e["target"], serde_json::to_value(&r.target).unwrap());
            assert_eq!(e["inputs"], serde_json::to_value(&r.inputs).unwrap());
            assert!(e["source_workspace_fingerprint"].is_string());
            let bytes = artifact(&c, &e["artifact_id"]).await;
            if format == "ProtocolPdf" {
                assert!(bytes.starts_with(b"%PDF-"));
                assert_eq!(e["filename"], "cloning-protocol-quantities.pdf");
            } else {
                let text = String::from_utf8(bytes).unwrap();
                assert!(
                    text.contains("PCR")
                        && text.contains("template-stock")
                        && text.contains("2000/1 nL")
                        && text.contains(&r.target.event_id)
                );
                assert_eq!(e["filename"], "cloning-protocol-quantities.html");
            }
        }
        let companion=call(&c,"protocol_quantities",json!({"workspace":reference,"task":{"Export":{"quantity_hash":captured.attachment_hash,"format":"Html"}}})).await;
        assert!(companion["molecular"].is_null());
        assert!(companion["source_workspace_fingerprint"].is_null());
        assert_eq!(workspace::save_workspace(&w).unwrap(), before);
    }
    c.cancel().await.unwrap();
    server.await.unwrap();
}
