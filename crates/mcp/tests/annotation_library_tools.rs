//! Real SDK contracts for immutable annotation/library evidence and explicit new-document application.
use base64::{Engine, engine::general_purpose::STANDARD};
use collagene_core::{
    Molecule,
    annotation_formats::{AnnotationFormat, AnnotationImportReport},
    feature_library::{FeatureLibrary, ScanOptions, ScanReport},
    protocol,
};
use collagene_mcp::CollageneServer;
use rmcp::{RoleClient, ServiceExt, model::CallToolRequestParams, service::RunningService};
use serde_json::{Value, json};
type Client = RunningService<RoleClient, ()>;
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
async fn call(client: &Client, name: &str, args: Value) -> Value {
    let r = client
        .call_tool(
            CallToolRequestParams::new(name.to_owned())
                .with_arguments(args.as_object().unwrap().clone()),
        )
        .await
        .unwrap();
    assert_ne!(r.is_error, Some(true), "{name}: {:?}", r.content);
    r.structured_content.unwrap()
}
async fn reject(client: &Client, name: &str, args: Value) {
    let r = client
        .call_tool(
            CallToolRequestParams::new(name.to_owned())
                .with_arguments(args.as_object().unwrap().clone()),
        )
        .await;
    assert!(
        r.is_err() || r.unwrap().is_error == Some(true),
        "accepted invalid {name}"
    );
}
fn reference(value: &Value) -> Value {
    json!({"id":value["id"],"revision":value["revision"]})
}
async fn import(client: &Client, text: &str) -> Value {
    reference(&call(client, "import_document", json!({"text":text})).await["documents"][0])
}
async fn artifact(client: &Client, id: &Value) -> Vec<u8> {
    let mut offset = 0;
    let mut result = String::new();
    let mut encoding = String::new();
    loop {
        let chunk = call(
            client,
            "read_artifact",
            json!({"artifact_id":id,"start":offset,"max_bytes":16384}),
        )
        .await;
        encoding.clone_from(&chunk["encoding"].as_str().unwrap().to_owned());
        result.push_str(chunk["text"].as_str().unwrap());
        if chunk["done"] == true {
            break;
        }
        let next = chunk["next_start"].as_u64().unwrap();
        assert!(next > offset);
        offset = next;
    }
    if encoding == "base64" {
        STANDARD.decode(result).unwrap()
    } else {
        result.into_bytes()
    }
}
async fn molecule(client: &Client, document: &Value) -> Molecule {
    let output = call(client, "export_document", document.clone()).await;
    protocol::load_document(&artifact(client, &output["artifact_id"]).await).unwrap()
}
const SOURCE: &str = "LOCUS       own_source                 6 bp    DNA     linear   SYN 01-JAN-2026\nFEATURES             Location/Qualifiers\n     misc_feature    1..6\n                     /label=\"Own source part\"\nORIGIN\n        1 cgatac\n//\n";
#[tokio::test]
async fn library_creation_export_scan_and_explicit_apply_keep_original_documents_immutable() {
    let (client, server) = setup().await;
    let source = import(&client, SOURCE).await;
    let target = import(&client, ">target\nTTTCGATACGGG\n").await;
    let created=call(&client,"create_feature_library",json!({"name":"Own parts","version":"1.0","features":[{"document":source,"feature_index":0}]})).await;
    assert_eq!(created["entry_count"], 1);
    let library = created["library"].clone();
    assert_eq!(
        call(&client, "list_feature_libraries", json!({})).await["libraries"]
            .as_array()
            .unwrap()
            .len(),
        1
    );
    let exported = call(
        &client,
        "export_feature_library",
        json!({"library":library}),
    )
    .await;
    let bytes = artifact(&client, &exported["artifact_id"]).await;
    let restored: FeatureLibrary = collagene_core::feature_library::load_library(&bytes).unwrap();
    assert_eq!(restored.entries[0].sequence.as_str(), "CGATAC");
    assert_eq!(
        restored.entries[0].source.document_revision,
        source["revision"].as_str().unwrap()
    );
    let again = call(
        &client,
        "import_feature_library",
        json!({"data_base64":STANDARD.encode(&bytes)}),
    )
    .await;
    assert_eq!(again["library"], library);
    let result = call(
        &client,
        "scan_feature_library",
        json!({"document":target,"library":library,"options":ScanOptions::default()}),
    )
    .await;
    assert_eq!(result["applied"], false);
    assert_eq!(result["complete"], true);
    assert_eq!(result["hit_count"], 1);
    assert_eq!(
        (
            result["hits"][0]["start"].as_u64(),
            result["hits"][0]["length"].as_u64()
        ),
        (Some(3), Some(6))
    );
    let full: ScanReport =
        serde_json::from_slice(&artifact(&client, &result["artifact_id"]).await).unwrap();
    assert_eq!(full.hits[0].concrete_matches, 6);
    assert_eq!(full.hits[0].compatible_ambiguous, 0);
    assert_eq!(full.hits[0].mismatches, 0);
    let annotated=call(&client,"apply_feature_library",json!({"document":target,"library":library,"scan":result["scan"],"selected":[full.hits[0].id]})).await;
    let new_ref = reference(&annotated["documents"][0]);
    assert_ne!(new_ref, target);
    let changed = molecule(&client, &new_ref).await;
    let original = molecule(&client, &target).await;
    assert_eq!(changed.sequence.as_str(), "TTTCGATACGGG");
    assert_eq!(changed.features[0].location.spans(), vec![(3, 9, false)]);
    assert!(changed.features[0].needs_review);
    assert!(original.features.is_empty());
    assert_eq!(call(&client, "list_documents", json!({})).await["total"], 3);
    drop(client);
    server.abort();
}
#[tokio::test]
async fn library_stale_tamper_limits_and_duplicates_reject_without_partial_document_or_library_creation()
 {
    let (client, server) = setup().await;
    let source = import(&client, SOURCE).await;
    let target = import(&client, ">target\nTTTCGATACGGG\n").await;
    reject(&client,"create_feature_library",json!({"name":"bad","version":"1","features":[{"document":source,"feature_index":0},{"document":source,"feature_index":0}]})).await;
    assert!(
        call(&client, "list_feature_libraries", json!({})).await["libraries"]
            .as_array()
            .unwrap()
            .is_empty()
    );
    let created = call(
        &client,
        "create_feature_library",
        json!({"name":"Own","version":"1","features":[{"document":source,"feature_index":0}]}),
    )
    .await;
    let library = created["library"].clone();
    let exported = call(
        &client,
        "export_feature_library",
        json!({"library":library}),
    )
    .await;
    let mut forged: Value =
        serde_json::from_slice(&artifact(&client, &exported["artifact_id"]).await).unwrap();
    forged["entries"][0]["name"] = json!("Unbound changed name");
    reject(
        &client,
        "import_feature_library",
        json!({"data_base64":STANDARD.encode(serde_json::to_vec(&forged).unwrap())}),
    )
    .await;
    let mut stale_library = library.clone();
    stale_library["revision"] = json!("stale");
    reject(
        &client,
        "scan_feature_library",
        json!({"document":target,"library":stale_library,"options":ScanOptions::default()}),
    )
    .await;
    let scan = call(
        &client,
        "scan_feature_library",
        json!({"document":target,"library":library,"options":ScanOptions::default()}),
    )
    .await;
    let selected = vec![scan["hits"][0]["id"].clone()];
    reject(&client,"apply_feature_library",json!({"document":target,"library":library,"scan":scan["scan"],"selected":["unknown proposal"]})).await;
    let mut stale_scan = scan["scan"].clone();
    stale_scan["revision"] = json!("stale");
    reject(
        &client,
        "apply_feature_library",
        json!({"document":target,"library":library,"scan":stale_scan,"selected":selected}),
    )
    .await;
    let other = import(&client, ">other_reference_same_bases\nTTTCGATACGGG\n").await;
    reject(
        &client,
        "apply_feature_library",
        json!({"document":other,"library":library,"scan":scan["scan"],"selected":selected}),
    )
    .await;
    let mut bounded = ScanOptions::default();
    bounded.max_comparisons = 1;
    let limited = call(
        &client,
        "scan_feature_library",
        json!({"document":target,"library":library,"options":bounded}),
    )
    .await;
    assert_eq!(limited["complete"], false);
    assert!(limited["comparisons"].as_u64().unwrap() <= 1);
    assert_eq!(call(&client, "list_documents", json!({})).await["total"], 3);
    assert_eq!(
        call(&client, "list_feature_libraries", json!({})).await["libraries"]
            .as_array()
            .unwrap()
            .len(),
        1
    );
    drop(client);
    server.abort();
}
#[tokio::test]
async fn gff_dataset_original_evidence_and_revision_bound_apply_export_are_explicit() {
    let (client, server) = setup().await;
    let target = import(&client, ">target\nTTTCGATACGGG\n").await;
    let text = "##gff-version 3\nchr1\tfixture\tmisc_feature\t3\t6\t.\t-\t.\tID=f1;Name=source%3Blabel\nchr2\tfixture\tmisc_feature\t1\t2\t.\t+\t.\tID=other\n";
    let imported = call(
        &client,
        "import_annotations",
        json!({"document":target,"text":text,"format":"Gff3","seqid":"chr1"}),
    )
    .await;
    assert_eq!(imported["applied"], false);
    assert_eq!(imported["proposal_count"], 1);
    assert_eq!(imported["skipped_seqids"], json!(["chr2"]));
    let dataset = imported["dataset"].clone();
    let exports = call(
        &client,
        "export_annotation_dataset",
        json!({"dataset":dataset}),
    )
    .await;
    assert_eq!(
        artifact(&client, &exports["original_artifact_id"]).await,
        text.as_bytes()
    );
    let full: AnnotationImportReport =
        serde_json::from_slice(&artifact(&client, &exports["dataset_artifact_id"]).await).unwrap();
    let original = molecule(&client, &target).await;
    full.validate(&original).unwrap();
    assert_eq!(
        full.proposals[0].feature.location.spans(),
        vec![(2, 6, true)]
    );
    reject(
        &client,
        "apply_annotations",
        json!({"document":target,"dataset":dataset,"selected":[0,0]}),
    )
    .await;
    assert_eq!(call(&client, "list_documents", json!({})).await["total"], 1);
    let applied = call(
        &client,
        "apply_annotations",
        json!({"document":target,"dataset":dataset,"selected":[0]}),
    )
    .await;
    let new_ref = reference(&applied["documents"][0]);
    let annotated = molecule(&client, &new_ref).await;
    assert_eq!(annotated.features[0].location.spans(), vec![(2, 6, true)]);
    assert_eq!(annotated.sequence, original.sequence);
    assert!(molecule(&client, &target).await.features.is_empty());
    reject(
        &client,
        "apply_annotations",
        json!({"document":new_ref,"dataset":dataset,"selected":[0]}),
    )
    .await;
    let exported = call(
        &client,
        "export_annotations",
        json!({"document":new_ref,"format":"Gff3","seqid":"chr1"}),
    )
    .await;
    let text = String::from_utf8(artifact(&client, &exported["artifact_id"]).await).unwrap();
    let reopened = collagene_core::annotation_formats::import_annotations(
        &original,
        &text,
        AnnotationFormat::Gff3,
        "chr1",
    )
    .unwrap();
    assert_eq!(
        reopened.proposals[0].feature.location.spans(),
        vec![(2, 6, true)]
    );
    assert_eq!(call(&client, "list_documents", json!({})).await["total"], 2);
    drop(client);
    server.abort();
}
#[tokio::test]
async fn bed_blocks_and_invalid_mapping_or_coordinates_do_not_create_hidden_documents() {
    let (client, server) = setup().await;
    let target = import(&client, ">reference\nAAAACCCCGGGGTTTT\n").await;
    let text = "chr1\t0\t8\tblocks\t20\t-\t0\t8\t0\t2\t2,2,\t0,6,\n";
    let imported = call(
        &client,
        "import_annotations",
        json!({"document":target,"text":text,"format":"Bed","seqid":"chr1"}),
    )
    .await;
    assert_eq!(imported["proposal_count"], 1);
    let mut stale = target.clone();
    stale["revision"] = json!("stale");
    reject(
        &client,
        "import_annotations",
        json!({"document":stale,"text":text,"format":"Bed","seqid":"chr1"}),
    )
    .await;
    reject(
        &client,
        "import_annotations",
        json!({"document":target,"text":"chr1\t0\t17\n","format":"Bed","seqid":"chr1"}),
    )
    .await;
    reject(
        &client,
        "import_annotations",
        json!({"document":target,"text":text,"format":"Bed","seqid":""}),
    )
    .await;
    let applied = call(
        &client,
        "apply_annotations",
        json!({"document":target,"dataset":imported["dataset"],"selected":[0]}),
    )
    .await;
    let annotated = molecule(&client, &reference(&applied["documents"][0])).await;
    assert_eq!(
        annotated.features[0].location.spans(),
        vec![(0, 2, true), (6, 8, true)]
    );
    assert_eq!(call(&client, "list_documents", json!({})).await["total"], 2);
    drop(client);
    server.abort();
}

#[tokio::test]
async fn feature_library_item_capacity_is_atomic_and_existing_refs_still_work() {
    let (client, server) = setup().await;
    let source = import(&client, SOURCE).await;
    let mut first = Value::Null;
    for index in 0..128 {
        let created=call(&client,"create_feature_library",json!({"name":format!("Own collection {index}"),"version":"1","features":[{"document":source,"feature_index":0}]})).await;
        if index == 0 {
            first = created["library"].clone();
        }
    }
    reject(&client,"create_feature_library",json!({"name":"Over limit","version":"1","features":[{"document":source,"feature_index":0}]})).await;
    assert_eq!(
        call(&client, "list_feature_libraries", json!({})).await["libraries"]
            .as_array()
            .unwrap()
            .len(),
        128
    );
    assert_eq!(call(&client, "list_documents", json!({})).await["total"], 1);
    let existing = call(&client, "export_feature_library", json!({"library":first})).await;
    let library = collagene_core::feature_library::load_library(
        &artifact(&client, &existing["artifact_id"]).await,
    )
    .unwrap();
    assert_eq!(library.name, "Own collection 0");
    assert_eq!(library.entries[0].sequence.as_str(), "CGATAC");
    drop(client);
    server.abort();
}
