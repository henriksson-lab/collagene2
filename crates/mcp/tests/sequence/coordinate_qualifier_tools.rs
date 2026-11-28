//! Real immutable reference edits, coordinate-qualifier projection and loss-aware
//! external/native export through the public MCP SDK. Expectations are literal.
use super::annotated_sequence_exports::{bytes, call, document, reject, setup, text};
use collagene_core::sequence::annotations::library::{FeatureLibrary, ScanOptions, ScanReport};
use serde_json::json;

// Qualifier positions are authored literally in the original eighteen-base frame.
const SOURCE: &str = "LOCUS       qualifier_source          18 bp    DNA     linear   SYN 01-JAN-2026\nFEATURES             Location/Qualifiers\n     CDS             5..10\n                     /label=\"Literal coordinate assertions\"\n                     /transl_except=\"(pos:5..7,aa:Sec)\"\n                     /anticodon=\"(pos:8..10,aa:Leu,seq:tac)\"\n                     /note=\"First retained note\"\n                     /note=\"Second retained note\"\n                     /translation=\"XX\"\n                     /codon_start=1\nORIGIN\n        1 aaaacgatac ggggtttt\n//\n";
fn expected_qualifiers() -> Vec<(String, Option<String>)> {
    [
        ("label", "Literal coordinate assertions"),
        ("transl_except", "(pos:5..7,aa:Sec)"),
        ("anticodon", "(pos:8..10,aa:Leu,seq:tac)"),
        ("note", "First retained note"),
        ("note", "Second retained note"),
        ("translation", "XX"),
        ("codon_start", "1"),
    ]
    .into_iter()
    .map(|(key, value)| (key.into(), Some(value.into())))
    .collect()
}

#[tokio::test]
async fn coordinate_qualifier_values_and_order_survive_native_and_genbank_sdk_roundtrips() {
    let (client, server) = setup().await;
    let imported = call(&client, "import_document", json!({"text":SOURCE})).await;
    let source = document(&imported);
    let exported = call(&client, "export_document", source.clone()).await;
    let before = bytes(&client, &exported).await;
    let native = collagene_core::protocol::load_document(&before).unwrap();
    assert_eq!(native.sequence.as_str(), "AAAACGATACGGGGTTTT");
    assert_eq!(native.features[0].location.spans(), vec![(4, 10, false)]);
    assert_eq!(native.features[0].qualifiers, expected_qualifiers());

    let output = call(
        &client,
        "export_annotated_sequence",
        json!({"document":source,"format":"GenBank","accept_omissions":true}),
    )
    .await;
    assert_eq!(output["state"], "prepared");
    let genbank = text(&client, &output).await;
    assert!(genbank.contains("(pos:5..7,aa:Sec)"));
    assert!(genbank.contains("(pos:8..10,aa:Leu,seq:tac)"));
    let imported = call(&client, "import_document", json!({"text":genbank})).await;
    let again = call(&client, "export_document", document(&imported)).await;
    let restored = collagene_core::protocol::load_document(&bytes(&client, &again).await).unwrap();
    assert_eq!(restored.sequence.as_str(), "AAAACGATACGGGGTTTT");
    assert_eq!(restored.features[0].qualifiers, expected_qualifiers());
    assert_eq!(restored.features[0].location.spans(), vec![(4, 10, false)]);
    assert_eq!(
        before,
        bytes(
            &client,
            &call(&client, "export_document", source.clone()).await
        )
        .await
    );
    let mut stale = source;
    stale["revision"] = json!("0".repeat(64));
    reject(&client, "export_document", stale.clone()).await;
    reject(
        &client,
        "export_annotated_sequence",
        json!({"document":stale,"format":"GenBank","accept_omissions":true}),
    )
    .await;
    drop(client);
    server.abort();
}

#[tokio::test]
async fn feature_library_preserves_coordinate_source_evidence_but_never_copies_it_to_new_frames() {
    let (client, server) = setup().await;
    let source = document(&call(&client, "import_document", json!({"text":SOURCE})).await);
    let before = bytes(
        &client,
        &call(&client, "export_document", source.clone()).await,
    )
    .await;
    let created = call(
        &client,
        "create_feature_library",
        json!({"name":"Literal local coordinate provenance","version":"1","features":[{"document":source,"feature_index":0}]}),
    ).await;
    let library = created["library"].clone();
    let output = call(
        &client,
        "export_feature_library",
        json!({"library":library}),
    )
    .await;
    let raw = text(&client, &output).await;
    let retained: FeatureLibrary = serde_json::from_str(&raw).unwrap();
    assert_eq!(retained.entries.len(), 1);
    assert_eq!(retained.entries[0].sequence.as_str(), "CGATAC");
    assert_eq!(
        retained.entries[0].source.location.spans(),
        vec![(4, 10, false)]
    );
    assert_eq!(retained.entries[0].qualifiers, expected_qualifiers());

    for (target_text, sequence, reverse) in [
        (">forward\nTTTCGATACGGG\n", "TTTCGATACGGG", false),
        (">reverse\nTTTGTATCGGGG\n", "TTTGTATCGGGG", true),
    ] {
        let target = document(&call(&client, "import_document", json!({"text":target_text})).await);
        let original = bytes(
            &client,
            &call(&client, "export_document", target.clone()).await,
        )
        .await;
        let scan = call(
            &client,
            "scan_feature_library",
            json!({"document":target,"library":library,"options":ScanOptions::default()}),
        )
        .await;
        assert_eq!(scan["complete"], true);
        assert_eq!(scan["hit_count"], 1);
        let report: ScanReport = serde_json::from_str(&text(&client, &scan).await).unwrap();
        let hit = &report.hits[0];
        assert_eq!((hit.start, hit.length, hit.reverse), (3, 6, reverse));
        let warnings = hit.warnings.join(" ");
        for key in ["transl_except", "anticodon", "translation", "codon_start"] {
            assert!(
                warnings.contains(key),
                "missing source-only loss evidence for {key}"
            );
        }
        let mut stale = target.clone();
        stale["revision"] = json!("0".repeat(64));
        reject(
            &client,
            "apply_feature_library",
            json!({"document":stale,"library":library,"scan":scan["scan"],"selected":[hit.id]}),
        )
        .await;
        let applied = call(
            &client,
            "apply_feature_library",
            json!({"document":target,"library":library,"scan":scan["scan"],"selected":[hit.id]}),
        )
        .await;
        let new_ref = document(&applied);
        let product = collagene_core::protocol::load_document(
            &bytes(
                &client,
                &call(&client, "export_document", new_ref.clone()).await,
            )
            .await,
        )
        .unwrap();
        assert_eq!(product.sequence.as_str(), sequence);
        let feature = &product.features[0];
        assert_eq!(feature.location.spans(), vec![(3, 9, reverse)]);
        assert!(feature.needs_review);
        assert!(feature.qualifiers.iter().all(|(key, _)| {
            !["transl_except", "anticodon", "translation", "codon_start"].contains(&key.as_str())
        }));
        assert_eq!(
            feature
                .qualifiers
                .iter()
                .filter(|(key, _)| key == "note")
                .map(|(_, value)| value.as_deref().unwrap())
                .collect::<Vec<_>>(),
            vec!["First retained note", "Second retained note"]
        );
        let exported = call(
            &client,
            "export_annotated_sequence",
            json!({"document":new_ref,"format":"GenBank","accept_omissions":true}),
        )
        .await;
        let genbank = text(&client, &exported).await;
        assert!(!genbank.contains("/transl_except="));
        assert!(!genbank.contains("/anticodon="));
        assert_eq!(
            original,
            bytes(&client, &call(&client, "export_document", target).await).await
        );
    }
    assert_eq!(
        before,
        bytes(&client, &call(&client, "export_document", source).await).await
    );
    drop(client);
    server.abort();
}

#[tokio::test]
async fn reference_transform_tools_project_literal_complete_local_coordinates_and_payloads() {
    let (client, server) = setup().await;
    let tools = client.list_all_tools().await.unwrap();
    let tool = tools
        .iter()
        .find(|tool| tool.name == "edit_reference_document")
        .unwrap();
    assert_eq!(tool.input_schema["additionalProperties"], false);
    let source = document(&call(&client, "import_document", json!({"text":SOURCE})).await);
    let before = bytes(
        &client,
        &call(&client, "export_document", source.clone()).await,
    )
    .await;
    let circular = document(
        &call(
            &client,
            "edit_reference_document",
            json!({"document":source,"operation":{"SetTopology":{"topology":"Circular"}}}),
        )
        .await,
    );
    for (input, operation, sequence, feature_location, transl_except, anticodon) in [
        (
            source.clone(),
            json!("ReverseComplement"),
            "AAAACCCCGTATCGTTTT",
            "complement(9..14)",
            "(pos:complement(12..14),aa:Sec)",
            "(pos:complement(9..11),aa:Leu,seq:tac)",
        ),
        (
            source.clone(),
            json!({"ExtractRegion":{"start":4,"end":10}}),
            "CGATAC",
            "1..6",
            "(pos:1..3,aa:Sec)",
            "(pos:4..6,aa:Leu,seq:tac)",
        ),
        (
            circular.clone(),
            json!({"RotateOrigin":{"new_origin":6}}),
            "ATACGGGGTTTTAAAACG",
            "join(17..18,1..4)",
            "(pos:join(17..18,1),aa:Sec)",
            "(pos:2..4,aa:Leu,seq:tac)",
        ),
    ] {
        let changed = call(
            &client,
            "edit_reference_document",
            json!({"document":input,"operation":operation}),
        )
        .await;
        assert_eq!(changed["reference_only"], true);
        assert_eq!(changed["source"], input);
        assert_eq!(changed["persistence"], "session only");
        let product_ref = document(&changed);
        let product = collagene_core::protocol::load_document(
            &bytes(
                &client,
                &call(&client, "export_document", product_ref.clone()).await,
            )
            .await,
        )
        .unwrap();
        assert_eq!(product.sequence.as_str(), sequence);
        let feature = &product.features[0];
        assert_eq!(
            collagene_model::annotations::insdc::location_text(&feature.location).unwrap(),
            feature_location
        );
        let positions = feature
            .qualifiers
            .iter()
            .filter(|(key, _)| key == "transl_except" || key == "anticodon")
            .map(|(key, value)| (key.as_str(), value.as_deref().unwrap()))
            .collect::<Vec<_>>();
        assert_eq!(
            positions,
            vec![("transl_except", transl_except), ("anticodon", anticodon)]
        );
        let exported = call(
            &client,
            "export_annotated_sequence",
            json!({"document":product_ref,"format":"GenBank","accept_omissions":true}),
        )
        .await;
        let genbank = text(&client, &exported).await;
        assert!(genbank.contains(transl_except));
        assert!(genbank.contains(anticodon));
        let roundtrip = document(&call(&client, "import_document", json!({"text":genbank})).await);
        let reloaded = collagene_core::protocol::load_document(
            &bytes(&client, &call(&client, "export_document", roundtrip).await).await,
        )
        .unwrap();
        assert_eq!(reloaded.sequence.as_str(), sequence);
        assert_eq!(reloaded.features[0].qualifiers, feature.qualifiers);
    }
    let renamed = document(
        &call(
            &client,
            "edit_reference_document",
            json!({"document":source,"operation":{"Rename":{"name":"Explicit reference name"}}}),
        )
        .await,
    );
    let renamed = collagene_core::protocol::load_document(
        &bytes(&client, &call(&client, "export_document", renamed).await).await,
    )
    .unwrap();
    assert_eq!(renamed.name, "Explicit reference name");
    assert_eq!(renamed.features[0].qualifiers, expected_qualifiers());
    let deleted = document(
        &call(
            &client,
            "edit_reference_document",
            json!({"document":source,"operation":{"DeleteFeature":{"index":0}}}),
        )
        .await,
    );
    let deleted = collagene_core::protocol::load_document(
        &bytes(&client, &call(&client, "export_document", deleted).await).await,
    )
    .unwrap();
    assert!(deleted.features.is_empty());
    assert_eq!(deleted.sequence.as_str(), "AAAACGATACGGGGTTTT");
    assert_eq!(
        before,
        bytes(
            &client,
            &call(&client, "export_document", source.clone()).await
        )
        .await
    );
    let published_count = call(&client, "list_documents", json!({})).await["total"].clone();
    reject(
        &client,
        "edit_reference_document",
        json!({"document":source,"operation":{"RotateOrigin":{"new_origin":1}}}),
    )
    .await;
    reject(
        &client,
        "edit_reference_document",
        json!({"document":circular,"operation":{"RotateOrigin":{"new_origin":6,"ignored":true}}}),
    )
    .await;
    let mut stale = source;
    stale["revision"] = json!("0".repeat(64));
    reject(
        &client,
        "edit_reference_document",
        json!({"document":stale,"operation":"ReverseComplement"}),
    )
    .await;
    assert_eq!(
        call(&client, "list_documents", json!({})).await["total"],
        published_count
    );
    drop(client);
    server.abort();
}

#[tokio::test]
async fn malformed_or_unsupported_coordinate_wrappers_fail_without_publishing_reference_products() {
    let (client, server) = setup().await;
    for wrapper in [
        "(pos:5..7,aa:Sec)junk",
        "(pos:J00194.1:5..7,aa:Sec)",
        "(pos:5..7,aa:Sec,unexpected:payload)",
        "(pos:5..8,aa:Sec)",
    ] {
        let source_text = SOURCE.replace("(pos:5..7,aa:Sec)", wrapper);
        let source = document(&call(&client, "import_document", json!({"text":source_text})).await);
        let before = bytes(
            &client,
            &call(&client, "export_document", source.clone()).await,
        )
        .await;
        let count = call(&client, "list_documents", json!({})).await["total"].clone();
        reject(
            &client,
            "edit_reference_document",
            json!({"document":source,"operation":"ReverseComplement"}),
        )
        .await;
        assert_eq!(
            call(&client, "list_documents", json!({})).await["total"],
            count
        );
        assert_eq!(
            before,
            bytes(&client, &call(&client, "export_document", source).await).await
        );
    }
    drop(client);
    server.abort();
}

#[tokio::test]
async fn manual_feature_edit_refuses_coordinate_movement_and_drops_noncoordinate_cds_claims() {
    use collagene_model::Location;
    let (client, server) = setup().await;
    let source = document(&call(&client, "import_document", json!({"text":SOURCE})).await);
    let before = bytes(
        &client,
        &call(&client, "export_document", source.clone()).await,
    )
    .await;
    let native = collagene_core::protocol::load_document(&before).unwrap();
    let mut same = native.features[0].clone();
    same.qualifiers[0].1 = Some("Explicit same-location label".into());
    let same_result = document(
        &call(
            &client,
            "edit_reference_document",
            json!({"document":source,"operation":{"SetFeature":{"index":0,"feature":same}}}),
        )
        .await,
    );
    let unchanged_geometry = collagene_core::protocol::load_document(
        &bytes(
            &client,
            &call(&client, "export_document", same_result).await,
        )
        .await,
    )
    .unwrap();
    assert_eq!(
        unchanged_geometry.features[0].qualifiers[1..],
        expected_qualifiers()[1..]
    );
    let mut moved = native.features[0].clone();
    moved.location = Location::Range {
        start: 3,
        end: 10,
        before: false,
        after: false,
    };
    let before_refusal = call(&client, "list_documents", json!({})).await["total"].clone();
    for remove_coordinate_keys in [false, true] {
        let mut proposal = moved.clone();
        if remove_coordinate_keys {
            proposal
                .qualifiers
                .retain(|(key, _)| key != "transl_except" && key != "anticodon");
        }
        reject(
            &client,
            "edit_reference_document",
            json!({"document":source,"operation":{"SetFeature":{"index":0,"feature":proposal}}}),
        )
        .await;
    }
    assert_eq!(
        call(&client, "list_documents", json!({})).await["total"],
        before_refusal
    );
    let plain_text = SOURCE
        .replace(
            "                     /transl_except=\"(pos:5..7,aa:Sec)\"\n",
            "",
        )
        .replace(
            "                     /anticodon=\"(pos:8..10,aa:Leu,seq:tac)\"\n",
            "",
        )
        .replace("/translation=", "/TrAnSlAtIoN=")
        .replace("/codon_start=", "/CoDoN_StArT=");
    let plain = document(&call(&client, "import_document", json!({"text":plain_text})).await);
    let plain_before = bytes(
        &client,
        &call(&client, "export_document", plain.clone()).await,
    )
    .await;
    let mut feature = collagene_core::protocol::load_document(&plain_before)
        .unwrap()
        .features
        .remove(0);
    feature.location = Location::Range {
        start: 3,
        end: 10,
        before: false,
        after: false,
    };
    let mut with_new_coordinates = feature.clone();
    with_new_coordinates.qualifiers.push((
        "AnTiCoDoN".into(),
        Some("(pos:8..10,aa:Leu,seq:tac)".into()),
    ));
    reject(&client,"edit_reference_document",json!({"document":plain,"operation":{"SetFeature":{"index":0,"feature":with_new_coordinates}}})).await;
    let changed = document(
        &call(
            &client,
            "edit_reference_document",
            json!({"document":plain,"operation":{"SetFeature":{"index":0,"feature":feature}}}),
        )
        .await,
    );
    let changed = collagene_core::protocol::load_document(
        &bytes(&client, &call(&client, "export_document", changed).await).await,
    )
    .unwrap();
    let feature = &changed.features[0];
    assert_eq!(feature.location.spans(), vec![(3, 10, false)]);
    assert!(feature.needs_review);
    assert!(
        feature
            .qualifiers
            .iter()
            .all(|(key, _)| !key.eq_ignore_ascii_case("translation")
                && !key.eq_ignore_ascii_case("codon_start"))
    );
    let notes = feature
        .qualifiers
        .iter()
        .filter(|(key, _)| key == "note")
        .filter_map(|(_, v)| v.as_deref())
        .collect::<Vec<_>>();
    assert!(
        notes
            .iter()
            .any(|v| v.contains("TrAnSlAtIoN") && v.contains("XX"))
    );
    assert!(
        notes
            .iter()
            .any(|v| v.contains("CoDoN_StArT") && v.contains('1'))
    );
    assert_eq!(
        before,
        bytes(&client, &call(&client, "export_document", source).await).await
    );
    assert_eq!(
        plain_before,
        bytes(&client, &call(&client, "export_document", plain).await).await
    );
    drop(client);
    server.abort();
}
