//! Actual SDK routing of independently retained design/consensus observations.
use base64::{Engine, engine::general_purpose::STANDARD};
use collagene_core::{
    alignment::{
        AlignmentOptions, ComparisonOptions, design_consensus as native, grid::GlobalStorage,
    },
    construct_qc::{AnnotationAssumption, annotation_revision},
    documents::workspace::{self, Workspace},
    material::lineage::EventBinding,
};
use collagene_mcp::server::CollageneServer;
use collagene_model::{Feature, Location, Molecule, Topology};
use rmcp::{RoleClient, ServiceExt, model::CallToolRequestParams, service::RunningService};
use serde_json::{Value, json};
type Client = RunningService<RoleClient, ()>;
async fn call(client: &Client, name: &str, args: Value) -> Value {
    let response = client
        .call_tool(
            CallToolRequestParams::new(name.to_owned())
                .with_arguments(args.as_object().unwrap().clone()),
        )
        .await;
    let response = response.unwrap_or_else(|error| panic!("{name}: {error:?}"));
    assert_ne!(
        response.is_error,
        Some(true),
        "{name}: {:?}",
        response.content
    );
    response.structured_content.unwrap()
}
async fn reject(client: &Client, args: Value) {
    let response = client
        .call_tool(
            CallToolRequestParams::new("assess_design_consensus")
                .with_arguments(args.as_object().unwrap().clone()),
        )
        .await;
    assert!(response.is_err() || response.unwrap().is_error == Some(true));
}
async fn artifact(client: &Client, id: &Value) -> Vec<u8> {
    let (mut offset, mut text, mut encoding) = (0, String::new(), String::new());
    loop {
        let page = call(
            client,
            "read_artifact",
            json!({"artifact_id":id,"start":offset,"max_bytes":16384}),
        )
        .await;
        text.push_str(page["text"].as_str().unwrap());
        encoding.clear();
        encoding.push_str(page["encoding"].as_str().unwrap());
        if page["done"] == true {
            break;
        }
        let next = page["next_start"].as_u64().unwrap();
        assert!(next > offset);
        offset = next;
    }
    if encoding == "base64" {
        STANDARD.decode(text).unwrap()
    } else {
        text.into_bytes()
    }
}
fn event(workspace: &Workspace, index: usize) -> EventBinding {
    let document = &workspace.documents[index];
    EventBinding {
        document_id: document.id.clone(),
        event_id: document.head_event.clone().unwrap(),
        revision: document.head.clone(),
    }
}
fn fixture() -> (Workspace, native::Request) {
    let mut design = Molecule::dna(
        "Independent two-CDS design",
        "ATGCAACCCTAACGTCGAATGGAACTGTAA",
        Topology::Linear,
    )
    .unwrap();
    for location in [Location::range(0, 12), Location::range(18, 30)] {
        design.features.push(Feature {
            kind: "CDS".into(),
            location,
            qualifiers: vec![],
            needs_review: false,
        });
    }
    let consensus = Molecule::dna(
        "Supplied unannotated observation",
        "ATGTAACCCTAACGTCGAATGGAGCTGTAA",
        Topology::Linear,
    )
    .unwrap();
    let cds = design.features.iter().enumerate().map(|(index,feature)| native::CdsSelection {
        id:format!("declared-CDS-{index}"), genetic_code:1,
        design_cds:AnnotationAssumption {feature_index:index,source_revision:design.revision(),annotation_revision:annotation_revision(feature).unwrap(),rationale:"Explicit design CDS coordinates; observed consensus has no inferred gene annotation".into()},
    }).collect();
    let workspace = Workspace::from_documents(
        "Immutable supplied-consensus review",
        vec![design, consensus],
        1,
    )
    .unwrap();
    let request = native::Request {
        schema_version: 1,
        design: event(&workspace, 0),
        consensus: event(&workspace, 1),
        comparison: ComparisonOptions {
            alignment: AlignmentOptions {
                global_storage: GlobalStorage::Full,
                ..Default::default()
            },
            both_query_strands: false,
            reference_rotations: false,
            max_origins: 1,
            max_total_cells: 100_000,
        },
        consensus_axis: native::ConsensusAxis::StoredTopology,
        cds,
        limits: Default::default(),
    };
    (workspace, request)
}

#[tokio::test]
async fn supplied_consensus_sdk_keeps_roles_conditional_coding_pages_and_complete_native_artifact()
{
    let (a, b) = tokio::io::duplex(65536);
    let service = tokio::spawn(async move {
        CollageneServer::default()
            .serve(a)
            .await
            .unwrap()
            .waiting()
            .await
            .unwrap();
    });
    let client = ().serve(b).await.unwrap();
    let (workspace, request) = fixture();
    let imported = call(
        &client,
        "import_workspace_archive",
        json!({"data_base64":STANDARD.encode(workspace::save_workspace(&workspace).unwrap())}),
    )
    .await;
    let reference =
        json!({"id":imported["workspace"]["id"],"revision":imported["workspace"]["revision"]});
    let before = call(&client, "list_documents", json!({})).await;
    let arguments = json!({"workspace":reference,"request":request,"cds_start":0,"max_cds":1});
    let summary = call(&client, "assess_design_consensus", arguments.clone()).await;
    assert_eq!(summary["cds_count"], 2);
    assert_eq!(summary["cds"].as_array().unwrap().len(), 1);
    assert_eq!(summary["next_start"], 1);
    assert_eq!(summary["done"], false);
    assert_eq!(summary["selected_path_conditional"], true);
    assert!(summary.get("verified").is_none());
    let full: native::Report =
        serde_json::from_slice(&artifact(&client, &summary["artifact_id"]).await).unwrap();
    full.validate_request(&workspace, &request).unwrap();
    assert_eq!(full.evidence.request.design, request.design);
    assert_eq!(full.evidence.request.consensus, request.consensus);
    assert_eq!(
        full.evidence.cds[0]
            .projected_translation
            .as_ref()
            .unwrap()
            .protein,
        "M*P*"
    );
    assert_eq!(
        full.evidence.cds[1]
            .projected_translation
            .as_ref()
            .unwrap()
            .protein,
        "MEL*"
    );
    assert!(
        full.evidence.cds[0]
            .consequences
            .iter()
            .any(
                |consequence| consequence.kind == native::ConsequenceKind::DefiniteInternalStop
                    && consequence.projected_codon_index == 1
            )
    );
    assert!(
        full.evidence.cds[1]
            .consequences
            .iter()
            .any(
                |consequence| consequence.kind == native::ConsequenceKind::Synonymous
                    && consequence.projected_codon_index == 1
            )
    );
    assert!(
        full.evidence
            .provenance
            .iter()
            .all(|source| source.original.state == native::OriginalState::Unavailable)
    );
    let mut next = arguments.clone();
    next["cds_start"] = json!(1);
    let next = call(&client, "assess_design_consensus", next).await;
    assert_eq!(next["cds"][0]["id"], full.evidence.cds[1].id);
    assert_eq!(next["done"], true);
    assert_eq!(next["evidence_revision"], summary["evidence_revision"]);
    let mut changed = request.clone();
    changed.cds[0].genetic_code = 6;
    assert!(
        full.validate_binding(&native::revision(&workspace, &changed).unwrap())
            .is_err()
    );
    let mut invalid = arguments.clone();
    invalid["request"]["design"]["document_id"] = json!(request.consensus.document_id);
    reject(&client, invalid).await;
    let mut invalid = arguments.clone();
    invalid["max_cds"] = json!(5);
    reject(&client, invalid).await;
    let mut invalid = arguments.clone();
    invalid["workspace"]["revision"] = json!("0".repeat(64));
    reject(&client, invalid).await;
    let mut limited = arguments;
    limited["request"]["limits"]["max_mapping_work"] = json!(1);
    let limited = call(&client, "assess_design_consensus", limited).await;
    assert_eq!(limited["coding_complete"], false);
    assert_eq!(call(&client, "list_documents", json!({})).await, before);
    let exported = call(&client, "export_workspace_archive", reference).await;
    let restored =
        workspace::load_workspace(&artifact(&client, &exported["artifact_id"]).await).unwrap();
    assert_eq!(
        serde_json::to_value(&restored).unwrap(),
        serde_json::to_value(&workspace).unwrap()
    );
    client.cancel().await.unwrap();
    service.await.unwrap();
}
