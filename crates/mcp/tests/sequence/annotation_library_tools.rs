//! Real SDK contracts for immutable annotation/library evidence and explicit new-document application.
use base64::{Engine, engine::general_purpose::STANDARD};
use collagene_core::{
    formats::annotations::{AnnotationFormat, AnnotationImportReport},
    protocol,
    sequence::annotations::library::{FeatureLibrary, ScanOptions, ScanReport},
};
use collagene_mcp::server::CollageneServer;
use collagene_model::Molecule;
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
    let restored: FeatureLibrary =
        collagene_core::sequence::annotations::library::load_library(&bytes).unwrap();
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
    let preview = call(
        &client,
        "preview_annotation_merge",
        json!({
            "document":target,
            "source":{"Library":{"library":library,"scan":result["scan"]}},
            "selections":[{"proposal_id":full.hits[0].id,"decision":"Append"}]
        }),
    )
    .await;
    assert_eq!(preview["applied"], false);
    assert_eq!(preview["complete"], true);
    assert_eq!(preview["candidate_count"], 1);
    assert_eq!(call(&client, "list_documents", json!({})).await["total"], 2);
    let annotated=call(&client,"apply_annotation_merge",json!({"document":target,"source":{"Library":{"library":library,"scan":result["scan"]}},"selections":[{"proposal_id":full.hits[0].id,"decision":"Append"}]})).await;
    let new_ref = reference(&annotated["documents"][0]);
    assert_ne!(new_ref, target);
    let changed = molecule(&client, &new_ref).await;
    let original = molecule(&client, &target).await;
    assert_eq!(changed.sequence.as_str(), "TTTCGATACGGG");
    assert_eq!(changed.features[0].location.spans(), vec![(3, 9, false)]);
    assert!(changed.features[0].needs_review);
    assert!(original.features.is_empty());
    let provenance = call(
        &client,
        "inspect_annotation_provenance",
        json!({"document":new_ref,"feature_index":0}),
    )
    .await;
    let ledger: collagene_core::sequence::annotations::merge::Ledger =
        serde_json::from_value(provenance["ledger"].clone()).unwrap();
    assert_eq!(ledger.contributions.len(), 1);
    let collagene_core::sequence::annotations::merge::SourcePin::Library {
        library_revision,
        proposal_id,
        source_revision,
        ..
    } = &ledger.contributions[0].source
    else {
        panic!("not library evidence")
    };
    assert_eq!(library_revision, library["revision"].as_str().unwrap());
    assert_eq!(proposal_id, &full.hits[0].id);
    assert_eq!(source_revision, source["revision"].as_str().unwrap());
    assert_eq!(annotated["native_history_created"], false);
    assert_eq!(annotated["physical_cloning_verified"], false);

    let fresh_scan = call(
        &client,
        "scan_feature_library",
        json!({"document":new_ref,"library":library,"options":ScanOptions::default()}),
    )
    .await;
    assert_eq!(fresh_scan["hit_count"], 1);
    let fresh_id = fresh_scan["hits"][0]["id"].clone();
    assert_ne!(fresh_id, json!(full.hits[0].id));
    let repeated = call(
        &client,
        "apply_annotation_merge",
        json!({
            "document":new_ref,
            "source":{"Library":{"library":library,"scan":fresh_scan["scan"]}},
            "selections":[{
                "proposal_id":fresh_id,
                "decision":{"AddProvenance":{
                    "target_index":0,
                    "annotation_revision":collagene_core::construct_qc::annotation_revision(&changed.features[0]).unwrap()
                }}
            }]
        }),
    ).await;
    assert_eq!(
        repeated["receipt"]["actions"][0]["contribution_added"],
        false
    );
    assert_eq!(reference(&repeated["documents"][0]), new_ref);
    let retained = call(
        &client,
        "inspect_annotation_provenance",
        json!({"document":new_ref,"feature_index":0}),
    )
    .await;
    let retained: collagene_core::sequence::annotations::merge::Ledger =
        serde_json::from_value(retained["ledger"].clone()).unwrap();
    assert_eq!(retained.contributions.len(), 1);
    assert_eq!(
        retained, ledger,
        "the first actual proposal context remains stored unchanged"
    );
    assert_eq!(molecule(&client, &new_ref).await, changed);

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
    let selected = vec![json!({"proposal_id":scan["hits"][0]["id"],"decision":"Append"})];
    reject(&client,"apply_annotation_merge",json!({"document":target,"source":{"Library":{"library":library,"scan":scan["scan"]}},"selections":[{"proposal_id":"unknown proposal","decision":"Append"}]})).await;
    let mut stale_scan = scan["scan"].clone();
    stale_scan["revision"] = json!("stale");
    reject(
        &client,
        "apply_annotation_merge",
        json!({"document":target,"source":{"Library":{"library":library,"scan":stale_scan}},"selections":selected}),
    )
    .await;
    let other = import(&client, ">other_reference_same_bases\nTTTCGATACGGG\n").await;
    reject(
        &client,
        "apply_annotation_merge",
        json!({"document":other,"source":{"Library":{"library":library,"scan":scan["scan"]}},"selections":selected}),
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
        "apply_annotation_merge",
        json!({"document":target,"source":{"Dataset":{"dataset":dataset}},"selections":[{"proposal_id":imported["proposals"][0]["id"],"decision":"Append"},{"proposal_id":imported["proposals"][0]["id"],"decision":"Append"}]}),
    )
    .await;
    assert_eq!(call(&client, "list_documents", json!({})).await["total"], 1);
    let applied = call(
        &client,
        "apply_annotation_merge",
        json!({"document":target,"source":{"Dataset":{"dataset":dataset}},"selections":[{"proposal_id":imported["proposals"][0]["id"],"decision":"Append"}]}),
    )
    .await;
    let new_ref = reference(&applied["documents"][0]);
    let annotated = molecule(&client, &new_ref).await;
    assert_eq!(annotated.features[0].location.spans(), vec![(2, 6, true)]);
    assert_eq!(annotated.sequence, original.sequence);
    assert!(molecule(&client, &target).await.features.is_empty());
    reject(
        &client,
        "apply_annotation_merge",
        json!({"document":new_ref,"source":{"Dataset":{"dataset":dataset}},"selections":[{"proposal_id":imported["proposals"][0]["id"],"decision":"Append"}]}),
    )
    .await;
    let exported = call(
        &client,
        "export_annotations",
        json!({"document":new_ref,"format":"Gff3","seqid":"chr1"}),
    )
    .await;
    let text = String::from_utf8(artifact(&client, &exported["artifact_id"]).await).unwrap();
    let reopened = collagene_core::formats::annotations::import_annotations(
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
        "apply_annotation_merge",
        json!({"document":target,"source":{"Dataset":{"dataset":imported["dataset"]}},"selections":[{"proposal_id":imported["proposals"][0]["id"],"decision":"Append"}]}),
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
    let library = collagene_core::sequence::annotations::library::load_library(
        &artifact(&client, &existing["artifact_id"]).await,
    )
    .unwrap();
    assert_eq!(library.name, "Own collection 0");
    assert_eq!(library.entries[0].sequence.as_str(), "CGATAC");
    drop(client);
    server.abort();
}

#[tokio::test]
async fn gtf_quoted_evidence_exact_mapping_and_new_document_roundtrip() {
    let (client, server) = setup().await;
    let target = import(&client, ">gtf_reference\nAAAACCCCGGGGTTTT\n").await;
    let original = molecule(&client, &target).await;
    let text = "# independently supplied GTF fixture\nchr1\tfixture\texon\t3\t6\t.\t-\t.\tgene_id \"g1\"; transcript_id \"t1\"; exon_number \"1\";\nchr2\tfixture\texon\t1\t2\t.\t+\t.\tgene_id \"g2\"; transcript_id \"t2\";\n";
    let imported = call(
        &client,
        "import_annotations",
        json!({"document":target,"text":text,"format":"Gtf","seqid":"chr1"}),
    )
    .await;
    assert_eq!(imported["format"], "Gtf");
    assert_eq!(imported["applied"], false);
    assert_eq!(imported["proposal_count"], 1);
    assert_eq!(imported["skipped_seqids"], json!(["chr2"]));
    let dataset = imported["dataset"].clone();
    let evidence = call(
        &client,
        "export_annotation_dataset",
        json!({"dataset":dataset}),
    )
    .await;
    assert_eq!(
        artifact(&client, &evidence["original_artifact_id"]).await,
        text.as_bytes()
    );
    let original_chunk = call(
        &client,
        "read_artifact",
        json!({"artifact_id":evidence["original_artifact_id"],"start":0,"max_bytes":16384}),
    )
    .await;
    assert_eq!(original_chunk["filename"], "original-annotations.gtf");
    let full: AnnotationImportReport =
        serde_json::from_slice(&artifact(&client, &evidence["dataset_artifact_id"]).await).unwrap();
    full.validate(&original).unwrap();
    assert_eq!(full.format, AnnotationFormat::Gtf);
    assert_eq!(full.proposals[0].source_rows, vec![2]);
    // GTF 3..6 inclusive is forward-reference [2,6), reverse oriented.
    assert_eq!(
        full.proposals[0].feature.location.spans(),
        vec![(2, 6, true)]
    );
    reject(
        &client,
        "apply_annotation_merge",
        json!({"document":target,"source":{"Dataset":{"dataset":dataset}},"selections":[{"proposal_id":imported["proposals"][0]["id"],"decision":"Append"},{"proposal_id":imported["proposals"][0]["id"],"decision":"Append"}]}),
    )
    .await;
    let mut stale = target.clone();
    stale["revision"] = json!("stale");
    reject(
        &client,
        "import_annotations",
        json!({"document":stale,"text":text,"format":"Gtf","seqid":"chr1"}),
    )
    .await;
    reject(&client, "import_annotations", json!({"document":target,"text":"chr1\tfixture\texon\t1\t17\t.\t+\t.\tgene_id \"g1\"; transcript_id \"t1\";\n","format":"Gtf","seqid":"chr1"})).await;
    reject(&client, "import_annotations", json!({"document":target,"text":"chr1\tfixture\texon\t1\t2\t.\t+\t.\tgene_id=g1;transcript_id=t1\n","format":"Gtf","seqid":"chr1"})).await;
    assert_eq!(call(&client, "list_documents", json!({})).await["total"], 1);
    let applied = call(
        &client,
        "apply_annotation_merge",
        json!({"document":target,"source":{"Dataset":{"dataset":dataset}},"selections":[{"proposal_id":imported["proposals"][0]["id"],"decision":"Append"}]}),
    )
    .await;
    let annotated_ref = reference(&applied["documents"][0]);
    let annotated = molecule(&client, &annotated_ref).await;
    assert_eq!(annotated.sequence, original.sequence);
    assert_eq!(annotated.features.len(), 1);
    assert_eq!(annotated.features[0].location.spans(), vec![(2, 6, true)]);
    assert!(molecule(&client, &target).await.features.is_empty());
    reject(
        &client,
        "apply_annotation_merge",
        json!({"document":annotated_ref,"source":{"Dataset":{"dataset":dataset}},"selections":[{"proposal_id":imported["proposals"][0]["id"],"decision":"Append"}]}),
    )
    .await;
    let exported = call(
        &client,
        "export_annotations",
        json!({"document":annotated_ref,"format":"Gtf","seqid":"chr1"}),
    )
    .await;
    let chunk = call(
        &client,
        "read_artifact",
        json!({"artifact_id":exported["artifact_id"],"start":0,"max_bytes":16384}),
    )
    .await;
    assert_eq!(chunk["filename"], "annotations.gtf");
    let output = String::from_utf8(artifact(&client, &exported["artifact_id"]).await).unwrap();
    let reopened = collagene_core::formats::annotations::import_annotations(
        &original,
        &output,
        AnnotationFormat::Gtf,
        "chr1",
    )
    .unwrap();
    assert_eq!(reopened.proposals.len(), 1);
    assert_eq!(
        reopened.proposals[0].feature.location.spans(),
        vec![(2, 6, true)]
    );
    assert_eq!(call(&client, "list_documents", json!({})).await["total"], 2);
    drop(client);
    server.abort();
}

#[tokio::test]
async fn donor_merge_resolves_exact_workspace_events_and_declared_reverse_rotation_without_mutating_sources()
 {
    use collagene_core::{
        documents::workspace::{self, Workspace},
        material::lineage::EventBinding,
        sequence::annotations::merge,
    };
    use collagene_model::{Feature, Location, Topology};
    let (client, server) = setup().await;
    let mut donor =
        Molecule::dna("Actual independent donor", "AACCGGTA", Topology::Circular).unwrap();
    donor.features.push(Feature {
        kind: "misc_feature".into(),
        location: Location::Complement(Box::new(Location::Join(vec![
            Location::range(6, 8),
            Location::range(0, 2),
        ]))),
        qualifiers: vec![(
            "label".into(),
            Some("Independently declared donor role".into()),
        )],
        needs_review: false,
    });
    let mut other = donor.clone();
    other.name = "Different actual event with identical calls".into();
    let frame = Workspace::from_documents(
        "Actual annotation donor event owners",
        vec![donor.clone(), other],
        0,
    )
    .unwrap();
    let binding = |index: usize| EventBinding {
        document_id: frame.documents[index].id.clone(),
        event_id: frame.documents[index].head_event.clone().unwrap(),
        revision: frame.documents[index].head.clone(),
    };
    let donor_event = binding(0);
    let other_event = binding(1);
    let imported = call(
        &client,
        "import_workspace_archive",
        json!({"data_base64":STANDARD.encode(workspace::save_workspace(&frame).unwrap())}),
    )
    .await;
    let workspace_ref = reference(&imported["workspace"]);
    let target_model = Molecule::dna(
        "Actual declared donor coordinate destination",
        "CCGGTTTA",
        Topology::Circular,
    )
    .unwrap();
    let wrong_model =
        Molecule::dna("Incompatible full DNA", "ACGGTTTA", Topology::Circular).unwrap();
    let target_text = collagene_core::formats::sequence::export_with_report(
        &target_model,
        collagene_core::formats::sequence::TextFormat::GenBank,
    )
    .unwrap()
    .text;
    let wrong_text = collagene_core::formats::sequence::export_with_report(
        &wrong_model,
        collagene_core::formats::sequence::TextFormat::GenBank,
    )
    .unwrap()
    .text;
    let target = import(&client, &target_text).await;
    let wrong = import(&client, &wrong_text).await;
    let arguments = json!({
        "document":target,
        "source":{"Donor":{"workspace":workspace_ref,"event":donor_event,"transform":{"origin":2,"reverse":true},"indices":[0]}},
        "selections":[{"proposal_id":"donor-0","decision":"Append"}]
    });
    let before = call(&client, "list_documents", json!({})).await;
    let preview = call(&client, "preview_annotation_merge", arguments.clone()).await;
    assert_eq!(preview["complete"], true);
    assert_eq!(call(&client, "list_documents", json!({})).await, before);
    let applied = call(&client, "apply_annotation_merge", arguments.clone()).await;
    let result_ref = reference(&applied["documents"][0]);
    let result = molecule(&client, &result_ref).await;
    assert_eq!(result.sequence.as_str(), "CCGGTTTA");
    assert_eq!(result.features.len(), 1);
    assert_eq!(
        result.features[0]
            .location
            .extract(&result.sequence)
            .unwrap(),
        "TTTA"
    );
    assert!(result.features[0].needs_review);
    assert!(molecule(&client, &target).await.features.is_empty());
    let provenance = call(
        &client,
        "inspect_annotation_provenance",
        json!({"document":result_ref,"feature_index":0}),
    )
    .await;
    let ledger: merge::Ledger = serde_json::from_value(provenance["ledger"].clone()).unwrap();
    assert_eq!(ledger.contributions.len(), 1);
    let merge::SourcePin::Donor {
        event,
        feature_index,
        transform,
        ..
    } = &ledger.contributions[0].source
    else {
        panic!("not donor evidence")
    };
    assert_eq!(event, &donor_event);
    assert_eq!(*feature_index, 0);
    assert!(transform.reverse);
    assert_eq!(transform.origin, 2);
    let before_refusals = call(&client, "list_documents", json!({})).await;
    for field in ["document_id", "event_id", "revision"] {
        let mut invalid = arguments.clone();
        invalid["source"]["Donor"]["event"][field] =
            serde_json::to_value(&other_event).unwrap()[field].clone();
        reject(&client, "apply_annotation_merge", invalid).await;
    }
    let mut invalid = arguments.clone();
    invalid["source"]["Donor"]["transform"]["reverse"] = json!(false);
    reject(&client, "apply_annotation_merge", invalid).await;
    let mut invalid = arguments;
    invalid["document"] = wrong;
    reject(&client, "apply_annotation_merge", invalid).await;
    assert_eq!(
        call(&client, "list_documents", json!({})).await,
        before_refusals
    );
    let exported = call(&client, "export_workspace_archive", workspace_ref).await;
    let reopened =
        workspace::load_workspace(&artifact(&client, &exported["artifact_id"]).await).unwrap();
    assert_eq!(reopened.revisions[&donor_event.revision], donor);
    drop(client);
    server.abort();
}
