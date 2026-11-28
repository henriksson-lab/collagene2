//! Native editing and annotation-only curation through the actual MCP SDK.
use super::annotated_sequence_exports::{Client, bytes, call, document, reject, setup, text};
use collagene_model::{Feature, Location, Molecule};
use serde_json::{Value, json};

const SOURCE: &str = "LOCUS       native_edit              15 bp    DNA     linear   SYN 01-JAN-2026\nFEATURES             Location/Qualifiers\n     CDS             4..12\n                     /label=\"Literal biological source\"\n                     /transl_except=\"(pos:7..9,aa:Sec)\"\n                     /anticodon=\"(pos:7..9,aa:Ser,seq:tga)\"\n                     /translation=\"MU\"\n                     /codon_start=1\nORIGIN\n        1 cccatgtgat aaccc\n//\n";
async fn load(c: &Client, d: Value) -> Molecule {
    collagene_core::protocol::load_document(&bytes(c, &call(c, "export_document", d).await).await)
        .unwrap()
}
fn values<'a>(feature: &'a Feature, key: &str) -> Vec<&'a str> {
    feature
        .qualifiers
        .iter()
        .filter(|(k, _)| k.eq_ignore_ascii_case(key))
        .filter_map(|(_, v)| v.as_deref())
        .collect()
}
fn inactive(feature: &Feature) {
    assert!(feature.needs_review);
    for key in ["transl_except", "anticodon", "translation", "codon_start"] {
        assert!(values(feature, key).is_empty(), "active /{key}");
    }
    for raw in ["(pos:7..9,aa:Sec)", "(pos:7..9,aa:Ser,seq:tga)", "MU"] {
        assert!(
            values(feature, "note").iter().any(|n| n.contains(raw)),
            "missing source evidence {raw}"
        );
    }
}
#[tokio::test]
async fn native_sequence_and_annotation_strand_operations_publish_only_immutable_reference_results()
{
    let (client, server) = setup().await;
    let source = document(&call(&client, "import_document", json!({"text":SOURCE})).await);
    let original = bytes(
        &client,
        &call(&client, "export_document", source.clone()).await,
    )
    .await;
    let changed=call(&client,"edit_reference_document",json!({"document":source,"operation":{"ReplaceSequence":{"start":0,"end":0,"replacement":"GG"}}})).await;
    assert_eq!(changed["reference_only"], true);
    assert_eq!(changed["review_required"], false);
    assert_eq!(changed["annotations_requiring_review"], 0);
    assert_eq!(changed["source"], source);
    assert_eq!(changed["persistence"], "session only");
    let product_ref = document(&changed);
    let product = load(&client, product_ref.clone()).await;
    assert_eq!(product.sequence.as_str(), "GGCCCATGTGATAACCC");
    assert_eq!(product.features[0].location, Location::range(5, 14));
    assert_eq!(
        values(&product.features[0], "transl_except"),
        vec!["(pos:9..11,aa:Sec)"]
    );
    assert_eq!(
        values(&product.features[0], "anticodon"),
        vec!["(pos:9..11,aa:Ser,seq:tga)"]
    );
    assert_eq!(values(&product.features[0], "translation"), vec!["MU"]);
    assert!(!product.features[0].needs_review);
    let exported = call(
        &client,
        "export_annotated_sequence",
        json!({"document":product_ref,"format":"GenBank","accept_omissions":true}),
    )
    .await;
    let restored = document(
        &call(
            &client,
            "import_document",
            json!({"text":text(&client,&exported).await}),
        )
        .await,
    );
    assert_eq!(load(&client, restored).await.features, product.features);
    let touched=call(&client,"edit_reference_document",json!({"document":source,"operation":{"ReplaceSequence":{"start":6,"end":6,"replacement":"A"}}})).await;
    assert_eq!(touched["review_required"], true);
    assert_eq!(touched["annotations_requiring_review"], 1);
    let touched = load(&client, document(&touched)).await;
    assert_eq!(touched.sequence.as_str(), "CCCATGATGATAACCC");
    assert_eq!(touched.features[0].location, Location::range(3, 13));
    inactive(&touched.features[0]);
    let flipped = call(
        &client,
        "edit_reference_document",
        json!({"document":source,"operation":{"FlipAnnotationStrands":{"indices":[0]}}}),
    )
    .await;
    assert_eq!(flipped["reference_only"], true);
    assert_eq!(flipped["review_required"], true);
    assert_eq!(flipped["annotations_requiring_review"], 1);
    let flipped = load(&client, document(&flipped)).await;
    assert_eq!(flipped.sequence.as_str(), "CCCATGTGATAACCC");
    assert_eq!(
        flipped.features[0].location,
        Location::Complement(Box::new(Location::range(3, 12)))
    );
    inactive(&flipped.features[0]);
    assert_eq!(
        bytes(&client, &call(&client, "export_document", source).await).await,
        original
    );
    drop(client);
    server.abort();
}
#[tokio::test]
async fn invalid_native_edit_requests_and_active_unmappable_claims_publish_nothing() {
    let (client, server) = setup().await;
    let source = document(&call(&client, "import_document", json!({"text":SOURCE})).await);
    let original = bytes(
        &client,
        &call(&client, "export_document", source.clone()).await,
    )
    .await;
    let count = call(&client, "list_documents", json!({})).await["total"].clone();
    for operation in [
        json!({"ReplaceSequence":{"start":16,"end":16,"replacement":"GG"}}),
        json!({"ReplaceSequence":{"start":5,"end":4,"replacement":"GG"}}),
        json!({"ReplaceSequence":{"start":0,"end":0,"replacement":"!"}}),
        json!({"ReplaceSequence":{"start":0,"end":0,"replacement":"GG","ignored":true}}),
        json!({"FlipAnnotationStrands":{"indices":[99]}}),
        json!({"FlipAnnotationStrands":{"indices":[0,0]}}),
        json!({"FlipAnnotationStrands":{"indices":[]}}),
    ] {
        reject(
            &client,
            "edit_reference_document",
            json!({"document":source,"operation":operation}),
        )
        .await;
    }
    let mut stale = source.clone();
    stale["revision"] = json!("0".repeat(64));
    reject(&client,"edit_reference_document",json!({"document":stale,"operation":{"ReplaceSequence":{"start":0,"end":0,"replacement":"GG"}}})).await;
    assert_eq!(
        call(&client, "list_documents", json!({})).await["total"],
        count
    );
    assert_eq!(
        bytes(&client, &call(&client, "export_document", source).await).await,
        original
    );
    for replacement in [
        "/transl_except=\"(pos:7..9,aa:Sec)junk\"",
        "/transl_except=\"(pos:NC_REMOTE:7..9,aa:Sec)\"",
        "/transl_except=\"(pos:5..7,aa:Sec)\"",
        "/coded_by=\"NC_REMOTE:1..9\"",
        "/rpt_unit_range=\"7..9\"",
    ] {
        let text = SOURCE.replace("/transl_except=\"(pos:7..9,aa:Sec)\"", replacement);
        let source = document(&call(&client, "import_document", json!({"text":text})).await);
        let original = bytes(
            &client,
            &call(&client, "export_document", source.clone()).await,
        )
        .await;
        let count = call(&client, "list_documents", json!({})).await["total"].clone();
        reject(&client,"edit_reference_document",json!({"document":source,"operation":{"ReplaceSequence":{"start":0,"end":0,"replacement":"GG"}}})).await;
        assert_eq!(
            call(&client, "list_documents", json!({})).await["total"],
            count
        );
        assert_eq!(
            bytes(&client, &call(&client, "export_document", source).await).await,
            original
        );
    }
    drop(client);
    server.abort();
}
