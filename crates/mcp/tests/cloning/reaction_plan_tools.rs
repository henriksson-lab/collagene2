//! Actual SDK retention/replay/export of source-bound documentary reaction setup.
#[path = "../../../core/tests/common/pdf.rs"]
mod pdf_readback;
use base64::{Engine as _, engine::general_purpose::STANDARD};
use collagene_core::{
    cloning::{
        ligation::recipe::Preview,
        restriction::{
            catalog::{self, CatalogSelection},
            conditions::{EvidenceStatus, ReactionConditions},
        },
    },
    planning::plan::Operation,
};
use collagene_mcp::server::CollageneServer;
use rmcp::{RoleClient, ServiceExt, model::CallToolRequestParams, service::RunningService};
use serde_json::{Value, json};

type Client = RunningService<RoleClient, ()>;
const ORIGINAL: &str = "AAGAATTCTT";
// EcoRI top/bottom boundaries 3/7; filling the selected right fragment's
// recessed strand leaves its original plus-reference DNA unchanged.
const TARGET: &str = "AATTCTT";

async fn setup() -> (Client, tokio::task::JoinHandle<()>) {
    let (server_io, client_io) = tokio::io::duplex(65536);
    let server = tokio::spawn(async move {
        CollageneServer::default()
            .serve(server_io)
            .await
            .unwrap()
            .waiting()
            .await
            .unwrap();
    });
    (().serve(client_io).await.unwrap(), server)
}
async fn call(client: &Client, name: &str, arguments: Value) -> Value {
    let response = client
        .call_tool(
            CallToolRequestParams::new(name.to_owned())
                .with_arguments(arguments.as_object().unwrap().clone()),
        )
        .await
        .unwrap();
    assert_ne!(
        response.is_error,
        Some(true),
        "{name}: {:?}",
        response.content
    );
    response.structured_content.unwrap()
}
async fn reject(client: &Client, name: &str, arguments: Value) {
    let response = client
        .call_tool(
            CallToolRequestParams::new(name.to_owned())
                .with_arguments(arguments.as_object().unwrap().clone()),
        )
        .await;
    assert!(response.is_err() || response.unwrap().is_error == Some(true));
}
async fn text(client: &Client, artifact: &Value) -> String {
    let mut output = String::new();
    let mut start = 0;
    loop {
        let chunk = call(
            client,
            "read_artifact",
            json!({"artifact_id":artifact,"start":start,"max_bytes":16384}),
        )
        .await;
        output.push_str(chunk["text"].as_str().unwrap());
        if chunk["done"] == true {
            return output;
        }
        start = chunk["next_start"].as_u64().unwrap();
    }
}
fn document(imported: &Value, name: &str) -> Value {
    let document = imported["documents"]
        .as_array()
        .unwrap()
        .iter()
        .find(|document| document["name"] == name)
        .unwrap();
    json!({"id":document["id"],"revision":document["revision"]})
}
fn selection() -> CatalogSelection {
    catalog::select(&catalog::builtin_catalog().unwrap(), &["ecori".into()]).unwrap()
}
fn conditions(formulation: Option<&str>) -> ReactionConditions {
    serde_json::from_value(json!({
        "schema_version":1,
        "formulations":[{"enzyme_id":"ecori","formulation_id":formulation}],
        "buffer":"rCutSmart",
        "temperature_milli_c":37000,
        "incubation_minutes":60,
        "minimum_activity_percent":75,
        "activity_context":"NEB unit-assay relative activity",
        "temperature_context":"vendor recommended incubation"
    }))
    .unwrap()
}
fn arguments(source: &Value, target: &Value, setup: Option<&ReactionConditions>) -> Value {
    json!({
        "inputs":[{"id":"source","document":source}],
        "target":target,
        "recipe":{
            "name":"Selected EcoRI product with separately retained documentary setup",
            "fragments":[{
                "source":{
                    "generation":null,
                    "role":"selected_fragment",
                    "input_id":"source",
                    "preparation":{"CatalogDigest":{
                        "selection":selection(),
                        "conditions":setup,
                        "fragment_index":1,
                        "repair":"FillFivePrimeOverhangs"
                    }},
                    "reverse":false
                },
                "treatment":null
            }],
            "circular":false,
            "allow_nicked":false
        }
    })
}
async fn preview(client: &Client, proposal: &Value) -> Preview {
    let mut evidence: Value =
        serde_json::from_str(&text(client, &proposal["artifact_id"]).await).unwrap();
    evidence.as_object_mut().unwrap().remove("diagnostics");
    serde_json::from_value(evidence).unwrap()
}
fn normalized(text: &str) -> String {
    text.split_whitespace().collect::<Vec<_>>().join(" ")
}

#[tokio::test]
async fn retained_recipe_conditions_change_evidence_without_changing_products_or_independent_target()
 {
    let (client, server) = setup().await;
    let imported = call(
        &client,
        "import_document",
        json!({"text":format!(">original\n{ORIGINAL}\n>independent_target\n{TARGET}\n")}),
    )
    .await;
    let source = document(&imported, "original");
    let target = document(&imported, "independent_target");
    let before = call(&client, "export_document", source.clone()).await;
    let inventory = call(&client, "list_documents", json!({})).await;
    let choices = [
        (None, None, "ConsistentUnderAssumptions"),
        (
            Some(conditions(None)),
            Some(EvidenceStatus::Unknown),
            "Incomplete",
        ),
        (
            Some(conditions(Some("neb-r0101"))),
            Some(EvidenceStatus::BelowThreshold),
            "Incomplete",
        ),
        (
            Some(conditions(Some("neb-r3101"))),
            Some(EvidenceStatus::Known),
            "ConsistentUnderAssumptions",
        ),
    ];
    let mut baseline = None;
    let mut prior_evidence =
        None::<collagene_core::cloning::restriction::conditions::source::SourceReport>;
    let mut plan_revisions = std::collections::BTreeSet::new();
    for (setup, expected_documentary_status, expected_status) in choices {
        let proposed = call(
            &client,
            "propose_ligation_assembly",
            arguments(&source, &target, setup.as_ref()),
        )
        .await;
        assert_eq!(proposed["verified"], false);
        assert_eq!(proposed["verification"]["status"], expected_status);
        assert!(plan_revisions.insert(proposed["plan_revision"].as_str().unwrap().to_owned()));
        let retained = preview(&client, &proposed).await;
        retained.validate().unwrap();
        assert_eq!(
            retained.compiled.plan.inputs[0].revision,
            source["revision"].as_str().unwrap()
        );
        assert_eq!(
            retained.compiled.plan.inputs[0].molecule.sequence.as_str(),
            ORIGINAL
        );
        assert_eq!(
            retained.compiled.plan.target_revision,
            target["revision"].as_str().unwrap()
        );
        assert_eq!(
            retained
                .verification
                .final_product()
                .expect("final declared step must reconstruct")
                .sequence
                .as_str(),
            TARGET
        );
        // Plan/source/evidence identity metadata must change when a setup changes.
        // Compare the explicit scientific projection, including both fresh physical
        // summaries; absence of either summary must fail this fixture.
        let reconstructed = retained
            .verification
            .steps
            .iter()
            .map(|step| {
                step.product
                    .as_ref()
                    .expect("each declared step must reconstruct")
            })
            .map(|product| {
                let duplex: collagene_core::chemistry::DuplexSummary = serde_json::from_str(
                    product
                        .metadata
                        .get("duplex_model")
                        .expect("fresh physical duplex summary"),
                )
                .unwrap();
                let modifications: collagene_core::chemistry::tracked_duplex::ModificationSummary =
                    serde_json::from_str(
                        product
                            .metadata
                            .get("tracked_modifications")
                            .expect("fresh tracked-material summary"),
                    )
                    .unwrap();
                json!({
                    "sequence":product.sequence,
                    "topology":product.topology,
                    "template_strand":product.metadata.get("pcr:template_strand"),
                    "features":product.features,
                    "duplex_model":duplex,
                    "tracked_modifications":modifications
                })
            })
            .collect::<Vec<_>>();
        if let Some(baseline) = &baseline {
            assert_eq!(
                &reconstructed, baseline,
                "Documentary uncertainty must not erase or alter physical reconstruction"
            );
        } else {
            baseline = Some(reconstructed);
        }
        let Operation::CatalogDigest {
            conditions: actual,
            selection: chosen,
            ..
        } = &retained.compiled.plan.steps[0]
        else {
            panic!("the recipe retains the catalog-digest operation")
        };
        assert_eq!(actual, &setup);
        assert_eq!(chosen, &selection());
        let documentary = retained
            .verification
            .step_reports
            .iter()
            .filter_map(|step| step.reaction_conditions.as_ref())
            .collect::<Vec<_>>();
        if let Some(expected) = expected_documentary_status {
            assert_eq!(documentary.len(), 1);
            let evidence = documentary[0];
            evidence.validate_historical().unwrap();
            if let Some(prior) = &prior_evidence {
                assert_ne!(prior.context_revision, evidence.context_revision);
                assert!(
                    prior
                        .report
                        .validate_binding(&evidence.report.request_revision)
                        .is_err()
                );
            }
            prior_evidence = Some(evidence.clone());
            let mut forged = evidence.clone();
            forged.source_revision = "0".repeat(64);
            assert!(forged.validate_historical().is_err());
            assert_eq!(
                evidence.source_revision,
                source["revision"].as_str().unwrap()
            );
            assert_eq!(evidence.material_revision.len(), 64);
            assert_eq!(evidence.report.documentary_status, expected);
            assert_eq!(&evidence.report.options, setup.as_ref().unwrap());
            evidence
                .report
                .validate_for(
                    &catalog::builtin_catalog().unwrap(),
                    &selection(),
                    setup.as_ref().unwrap(),
                )
                .unwrap();
            assert!(
                evidence
                    .report
                    .sources
                    .iter()
                    .any(|s| s.url.starts_with("https://"))
            );
        } else {
            assert!(documentary.is_empty());
        }
        let reference = json!({"plan_revision":proposed["plan_revision"]});
        let fresh = call(&client, "verify_plan", reference.clone()).await;
        assert_eq!(fresh["status"], expected_status);
        assert_eq!(
            fresh["steps"].as_array().unwrap().len(),
            retained.verification.reconstructed_count()
        );
        assert_eq!(fresh["target_revision"], target["revision"]);
        if expected_documentary_status.is_some() {
            let full: Value =
                serde_json::from_str(&text(&client, &fresh["artifact_id"]).await).unwrap();
            assert_eq!(
                full["step_reports"],
                serde_json::to_value(&retained.verification.step_reports).unwrap()
            );
        }
        let exported = call(&client, "export_protocol", reference.clone()).await;
        let html = text(&client, &exported["artifact_id"]).await;
        assert!(html.contains("original"));
        assert!(html.contains(target["revision"].as_str().unwrap()));
        let pdf = call(&client, "export_protocol_pdf", reference).await;
        assert_eq!(pdf["mime_type"], "application/pdf");
        assert_eq!(pdf["status"], expected_status);
        assert_eq!(pdf["target_revision"], target["revision"]);
        let bytes = STANDARD
            .decode(text(&client, &pdf["artifact_id"]).await)
            .unwrap();
        assert!(bytes.starts_with(b"%PDF-"));
        pdf_readback::check_xref(&bytes);
        let printed = normalized(&pdf_readback::extracted(&bytes));
        if let Some(expected) = expected_documentary_status {
            for content in [&html, &printed] {
                assert!(content.contains("Documentary reaction-condition evidence"));
                assert!(content.contains("rCutSmart"));
                assert!(content.contains(&format!("{expected:?}")));
                assert!(normalized(content).contains(
                    "Star activity and incubation-duration compatibility remain unassessed"
                ));
                if let Some(id) = &setup.as_ref().unwrap().formulations[0].formulation_id {
                    assert!(
                        content.contains(id),
                        "Exact selected formulation missing from protocol"
                    );
                }
            }
        }
    }
    assert_eq!(call(&client, "list_documents", json!({})).await, inventory);
    assert_eq!(
        call(&client, "export_document", source.clone()).await["artifact_id"],
        before["artifact_id"]
    );
    assert_eq!(
        call(
            &client,
            "read_sequence_region",
            json!({"document":source,"start":0,"end":10})
        )
        .await["sequence"],
        ORIGINAL
    );
    client.cancel().await.unwrap();
    server.await.unwrap();
}

#[tokio::test]
async fn recipe_conditions_require_current_explicit_null_and_reject_mismatched_selection_or_catalog()
 {
    let (client, server) = setup().await;
    let imported = call(
        &client,
        "import_document",
        json!({"text":format!(">original\n{ORIGINAL}\n>independent_target\n{TARGET}\n")}),
    )
    .await;
    let source = document(&imported, "original");
    let target = document(&imported, "independent_target");
    let initial = call(&client, "list_documents", json!({})).await;
    let valid = arguments(&source, &target, Some(&conditions(Some("neb-r3101"))));
    let mut missing = valid.clone();
    missing["recipe"]["fragments"][0]["source"]["preparation"]["CatalogDigest"]
        .as_object_mut()
        .unwrap()
        .remove("conditions");
    reject(&client, "propose_ligation_assembly", missing).await;
    let mut stale = valid.clone();
    stale["recipe"]["fragments"][0]["source"]["preparation"]["CatalogDigest"]["selection"]["catalog_content_hash"] =
        json!("0".repeat(64));
    reject(&client, "propose_ligation_assembly", stale).await;
    let mut unrelated = valid.clone();
    unrelated["recipe"]["fragments"][0]["source"]["preparation"]["CatalogDigest"]["conditions"]["formulations"]
        [0]["enzyme_id"] = json!("bamhi");
    reject(&client, "propose_ligation_assembly", unrelated).await;
    let mut foreign_formulation = valid.clone();
    foreign_formulation["recipe"]["fragments"][0]["source"]["preparation"]["CatalogDigest"]["conditions"]
        ["formulations"][0]["formulation_id"] = json!("neb-r0141");
    reject(&client, "propose_ligation_assembly", foreign_formulation).await;
    let mut wrong_target = valid;
    wrong_target["target"] = source;
    let wrong = call(&client, "propose_ligation_assembly", wrong_target).await;
    assert_eq!(wrong["verification"]["status"], "Invalid");
    assert_eq!(
        call(
            &client,
            "verify_plan",
            json!({"plan_revision":wrong["plan_revision"]})
        )
        .await["status"],
        "Invalid"
    );
    assert_eq!(call(&client, "list_documents", json!({})).await, initial);
    client.cancel().await.unwrap();
    server.await.unwrap();
}
