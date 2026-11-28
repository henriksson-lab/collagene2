//! Real in-memory SDK and binary-worker contracts over typed native history.
use base64::{Engine, engine::general_purpose::STANDARD};
use collagene_core::{
    construct_qc::{self, causal},
    documents::{
        history::{
            DocumentHistory,
            transition::{AnalysisOrigin, AnalysisProducer},
        },
        reference_edit::Operation as Edit,
        workspace::{self, Workspace},
    },
    material::lineage::EventBinding,
    protocol::{
        self,
        request::Request as WorkerRequest,
        response::Response as WorkerResponse,
        schema::{AnalysisJob, AnalysisResult},
    },
};
use collagene_mcp::server::CollageneServer;
use collagene_model::{Feature, Location, Molecule, Topology};
use rmcp::{RoleClient, ServiceExt, model::CallToolRequestParams, service::RunningService};
use serde_json::{Value, json};
type Client = RunningService<RoleClient, ()>;
async fn call(client: &Client, name: &str, args: Value) -> Value {
    let result = client
        .call_tool(
            CallToolRequestParams::new(name.to_owned())
                .with_arguments(args.as_object().unwrap().clone()),
        )
        .await
        .unwrap();
    assert_ne!(result.is_error, Some(true), "{name}: {:?}", result.content);
    result.structured_content.unwrap()
}
async fn reject(client: &Client, args: Value) {
    match client
        .call_tool(
            CallToolRequestParams::new("assess_construct_causal_qc")
                .with_arguments(args.as_object().unwrap().clone()),
        )
        .await
    {
        Ok(result) => assert_eq!(result.is_error, Some(true)),
        Err(rmcp::service::ServiceError::McpError(error)) => assert_eq!(error.code.0, -32602),
        Err(error) => panic!("unexpected transport error: {error}"),
    }
}
async fn artifact(client: &Client, id: &Value) -> Value {
    let mut text = String::new();
    let mut start = 0;
    loop {
        let page = call(
            client,
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
fn binding(history: &DocumentHistory) -> EventBinding {
    EventBinding {
        document_id: history.id.clone(),
        event_id: history.head_event().into(),
        revision: history.head_revision().into(),
    }
}
fn assumption(source: &Molecule) -> construct_qc::AnnotationAssumption {
    construct_qc::AnnotationAssumption{feature_index:0,source_revision:source.revision(),annotation_revision:construct_qc::annotation_revision(&source.features[0]).unwrap(),rationale:"Declared synthetic CDS for exact native coordinate comparison; no biological identity claim".into()}
}
fn fixture() -> (Workspace, causal::Request) {
    let mut source =
        Molecule::dna("original synthetic CDS", "ATGCAACCCTAA", Topology::Linear).unwrap();
    source.features.push(Feature {
        kind: "CDS".into(),
        location: Location::range(0, 12),
        qualifiers: vec![("codon_start".into(), Some("1".into()))],
        needs_review: false,
    });
    let mut history = DocumentHistory::new("cds-document", source.clone()).unwrap();
    let original = binding(&history);
    history
        .commit(
            vec![Edit::Rename {
                name: "renamed synthetic CDS".into(),
            }],
            "Actual presentation edit",
        )
        .unwrap();
    let renamed = binding(&history);
    let renamed_assumption = assumption(&history.molecule);
    let mut preview = history.molecule.clone();
    preview.name = "Unretained computed reference".into();
    let origin = AnalysisOrigin::preview(AnalysisProducer::AssemblyIntent, preview.revision());
    history
        .commit_analysis(preview, "Explicit unretained design preview", origin, &[])
        .unwrap();
    let computed = binding(&history);
    let request = causal::Request {
        schema_version: 1,
        pairs: vec![
            causal::CdsPair {
                id: "native-rename".into(),
                source: original.clone(),
                product: renamed,
                source_cds: assumption(&source),
                product_cds: renamed_assumption,
                genetic_code: 1,
            },
            causal::CdsPair {
                id: "unretained-preview".into(),
                source: original,
                product: computed,
                source_cds: assumption(&source),
                product_cds: assumption(&history.molecule),
                genetic_code: 1,
            },
        ],
        limits: causal::Limits::default(),
    };
    let mut workspace = Workspace::from_documents("causal SDK fixture", vec![source], 0).unwrap();
    workspace.documents.clear();
    workspace.history.clear();
    workspace.revisions.clear();
    history.contribute(&mut workspace).unwrap();
    workspace.active_document = history.id.clone();
    workspace.validate().unwrap();
    (workspace, request)
}
#[tokio::test]
async fn native_history_pairs_roundtrip_full_artifacts_and_keep_preview_gaps_explicit() {
    let (workspace, request) = fixture();
    let context = causal::revision(&workspace, &request).unwrap();
    let wire = WorkerRequest {
        id: 707,
        revision: context.clone(),
        job: AnalysisJob::AssessConstructCausalQc {
            workspace: workspace.clone(),
            request: request.clone(),
        },
    };
    let bytes = wire.encode().unwrap();
    assert_eq!(u16::from_le_bytes(bytes[4..6].try_into().unwrap()), 137);
    let response =
        WorkerResponse::decode(&protocol::execution::process_request(&bytes).unwrap()).unwrap();
    let AnalysisResult::ConstructCausalQc(report) = response.result else {
        panic!("typed causal worker response required");
    };
    report.validate_for(&workspace, &request).unwrap();
    assert!(report.pairs[0].mapping_complete);
    assert!(!report.pairs[1].mapping_complete);
    assert!(
        report.pairs[1]
            .findings
            .iter()
            .any(|finding| finding.rule == causal::Rule::NativeHistory
                && finding.verdict == construct_qc::Verdict::Indeterminate)
    );
    let mut stale = wire.clone();
    stale.revision = "0".repeat(64);
    assert!(stale.encode().is_err());
    let mut obsolete = bytes.clone();
    obsolete[4..6].copy_from_slice(&136u16.to_le_bytes());
    assert!(protocol::execution::process_request(&obsolete).is_err());
    let (client, server) = setup().await;
    let imported = call(
        &client,
        "import_workspace_archive",
        json!({"data_base64":STANDARD.encode(workspace::save_workspace(&workspace).unwrap())}),
    )
    .await;
    let reference = json!({
        "id": imported["workspace"]["id"],
        "revision": imported["workspace"]["revision"],
    });
    let before = call(&client, "list_workspaces", json!({})).await;
    let first = call(
        &client,
        "assess_construct_causal_qc",
        json!({"workspace":reference,"request":request,"max_pairs":1}),
    )
    .await;
    assert_eq!(first["pair_count"], 2);
    assert_eq!(first["next_start"], 1);
    assert_eq!(first["done"], false);
    assert_eq!(first["pairs"][0]["id"], "native-rename");
    assert_eq!(first["pairs"][0]["mapping_complete"], true);
    let full: causal::Report =
        serde_json::from_value(artifact(&client, &first["artifact_id"]).await).unwrap();
    full.validate_for(&workspace, &request).unwrap();
    assert_eq!(full, report);
    let second = call(
        &client,
        "assess_construct_causal_qc",
        json!({"workspace":reference,"request":request,"pair_start":1,"max_pairs":1}),
    )
    .await;
    assert_eq!(second["pairs"][0]["mapping_complete"], false);
    assert_eq!(second["done"], true);
    assert_eq!(first["artifact_id"], second["artifact_id"]);
    assert_eq!(call(&client, "list_workspaces", json!({})).await, before);
    client.cancel().await.unwrap();
    server.abort();
}
#[tokio::test]
async fn stale_sources_oversized_requests_and_invalid_pages_never_change_the_workspace() {
    let (workspace, request) = fixture();
    let (client, server) = setup().await;
    let imported = call(
        &client,
        "import_workspace_archive",
        json!({"data_base64":STANDARD.encode(workspace::save_workspace(&workspace).unwrap())}),
    )
    .await;
    let reference = json!({
        "id": imported["workspace"]["id"],
        "revision": imported["workspace"]["revision"],
    });
    let before = call(&client, "list_workspaces", json!({})).await;
    let mut stale = reference.clone();
    stale["revision"] = json!("0".repeat(64));
    reject(&client, json!({"workspace":stale,"request":request})).await;
    let mut altered = request.clone();
    altered.pairs[0].source.revision = "0".repeat(64);
    let stale_pair = call(
        &client,
        "assess_construct_causal_qc",
        json!({"workspace":reference,"request":altered}),
    )
    .await;
    assert_eq!(stale_pair["pairs"][0]["mapping_complete"], false);
    let mut oversized = request.clone();
    oversized.pairs = vec![request.pairs[0].clone(); 33];
    reject(&client, json!({"workspace":reference,"request":oversized})).await;
    reject(
        &client,
        json!({"workspace":reference,"request":request,"max_pairs":0}),
    )
    .await;
    reject(
        &client,
        json!({"workspace":reference,"request":request,"pair_start":3}),
    )
    .await;
    let mut unsupported = request.clone();
    unsupported.pairs[0].genetic_code = 255;
    reject(
        &client,
        json!({"workspace":reference,"request":unsupported}),
    )
    .await;
    assert_eq!(call(&client, "list_workspaces", json!({})).await, before);
    client.cancel().await.unwrap();
    server.abort();
}
