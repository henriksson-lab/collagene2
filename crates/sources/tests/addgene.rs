//! Independently authored provider fixtures use the public OpenAPI field definitions, not live
//! account data or another cloning application's output. No authenticated network is required.
use collagene_core::protocol;
use collagene_model::Topology;
use collagene_sources::*;
use serde_json::{Value, json};

const NOW: &str = "2026-11-23T18:00:00Z";
const GB: &str = "LOCUS       ADDGENE_TEST              16 bp    DNA     circular SYN 01-JAN-2026\nDEFINITION  Independently authored Addgene adapter fixture.\nFEATURES             Location/Qualifiers\n     CDS             complement(join(13..16,1..4))\n                     /label=\"origin feature\"\nORIGIN\n        1 aaaaccccgg ggtttt\n//\n";
fn sequence(id: u64) -> Value {
    json!({"sequence_id":id,"sequence_description":format!("Independent sequence {id}"),"length":16,"sequence":"aaaaccccggggtttt","genbank_api_url":format!("https://files.example.org/{id}.gb?signature=TRANSIENT_SECRET"),"genbank_url":format!("https://files.example.org/{id}.gb?signature=TRANSIENT_SECRET")})
}
fn fixture() -> Value {
    json!({"id":12345,"name":"Independent Addgene plasmid","sequences":{
        "public_user_full_sequences":[sequence(41)],"public_addgene_full_sequences":[sequence(42)],
        "public_user_partial_sequences":[sequence(43)],"public_addgene_partial_sequences":[sequence(44)]
    },"description":"Unknown catalog metadata is not interpreted as scientific sequence evidence."})
}
fn parse(value: &Value) -> Result<AddgeneCatalog, SourceError> {
    parse_addgene_catalog(12345, &serde_json::to_vec(value).unwrap())
}
#[test]
fn documented_four_categories_stable_ids_and_transient_downloads_are_separate() {
    let raw = fixture();
    let catalog = parse(&raw).unwrap();
    assert_eq!(catalog.plasmid_id(), 12345);
    assert_eq!(catalog.choices().len(), 4);
    assert_eq!(
        catalog.choices()[0].origin,
        AddgeneSequenceOrigin::Depositor
    );
    assert_eq!(catalog.choices()[1].origin, AddgeneSequenceOrigin::Addgene);
    assert_eq!(catalog.choices()[0].completeness, AddgeneCompleteness::Full);
    assert_eq!(
        catalog.choices()[2].completeness,
        AddgeneCompleteness::Partial
    );
    assert_eq!(
        catalog
            .choices()
            .iter()
            .map(|c| c.sequence_id)
            .collect::<Vec<_>>(),
        [41, 42, 43, 44]
    );
    assert!(
        catalog
            .download_url(42)
            .unwrap()
            .contains("TRANSIENT_SECRET")
    );
    let safe = format!("{catalog:?}{}", serde_json::to_string(&catalog).unwrap());
    assert!(!safe.contains("TRANSIENT_SECRET"));
    assert!(!safe.contains("files.example.org"));
    assert_eq!(
        catalog.response_blake3(),
        blake3::hash(&serde_json::to_vec(&raw).unwrap())
            .to_hex()
            .as_str()
    );
    assert_eq!(
        catalog.download_url(99).unwrap_err(),
        SourceError::AddgeneUnknownSequence
    );
    assert_eq!(
        addgene_catalog_url(12345).unwrap(),
        "https://api.developers.addgene.org/catalog/plasmid-with-sequences/12345/"
    );
    assert!(addgene_catalog_url(0).is_err());
}
#[test]
fn selected_genbank_keeps_features_original_evidence_and_secret_free_provider_provenance() {
    let catalog = parse(&fixture()).unwrap();
    let report = catalog.import_download(42, GB.as_bytes(), NOW).unwrap();
    let molecule = &report.molecules[0];
    assert_eq!(molecule.sequence.as_str(), "AAAACCCCGGGGTTTT");
    assert_eq!(molecule.topology, Topology::Circular);
    assert_eq!(
        molecule.features[0]
            .location
            .extract(&molecule.sequence)
            .unwrap(),
        "TTTTAAAA"
    );
    assert_eq!(molecule.original_genbank.as_deref(), Some(GB));
    assert_eq!(molecule.metadata["source:provider"], "Addgene");
    assert_eq!(molecule.metadata["source:identifier"], "12345");
    assert_eq!(molecule.metadata["source:sequence_id"], "42");
    assert_eq!(molecule.metadata["source:addgene_origin"], "addgene");
    assert_eq!(molecule.metadata["source:addgene_completeness"], "full");
    assert_eq!(
        molecule.metadata["source:canonical_url"],
        "https://www.addgene.org/12345/"
    );
    assert!(
        !serde_json::to_string(molecule)
            .unwrap()
            .contains("TRANSIENT_SECRET")
    );
    assert_eq!(
        protocol::load_document(&protocol::save_document(molecule).unwrap())
            .unwrap()
            .revision(),
        molecule.revision()
    );
    let partial = catalog
        .import_download(43, GB.replace("circular", "linear  ").as_bytes(), NOW)
        .unwrap();
    assert!(
        partial
            .warnings
            .iter()
            .any(|w| w.contains("explicitly partial"))
    );
    assert_eq!(
        partial.molecules[0].metadata["source:addgene_completeness"],
        "partial"
    );
    assert_eq!(
        catalog.import_download(43, GB.as_bytes(), NOW).unwrap_err(),
        SourceError::AddgeneSequenceMismatch
    );
}
#[test]
fn schema_changes_duplicate_ids_and_wrong_lengths_or_bases_fail_explicitly() {
    let mut value = fixture();
    value["id"] = json!(9);
    assert_eq!(parse(&value).unwrap_err(), SourceError::IdentifierMismatch);
    let mut value = fixture();
    value["sequences"]["future_sequence_kind"] = json!([]);
    assert_eq!(parse(&value).unwrap_err(), SourceError::AddgeneSchema);
    let mut value = fixture();
    value["sequences"]["public_user_full_sequences"] = json!([sequence(42)]);
    assert_eq!(parse(&value).unwrap_err(), SourceError::AddgeneSchema);
    let mut value = fixture();
    value["sequences"]["public_user_full_sequences"][0]["length"] = json!(17);
    assert_eq!(parse(&value).unwrap_err(), SourceError::Truncated);
    let catalog = parse(&fixture()).unwrap();
    assert_eq!(
        catalog
            .import_download(
                41,
                GB.replace("aaaaccccgg ggtttt", "aaaaccccgg ggttta")
                    .as_bytes(),
                NOW
            )
            .unwrap_err(),
        SourceError::AddgeneSequenceMismatch
    );
    assert_eq!(
        catalog
            .import_download(41, b">not_genbank\nACGT\n", NOW)
            .unwrap_err(),
        SourceError::UnexpectedContent
    );
    assert_eq!(
        catalog
            .import_download(41, format!("{GB}{GB}").as_bytes(), NOW)
            .unwrap_err(),
        SourceError::AmbiguousResponse { records: 2 }
    );
}
#[test]
fn deprecated_links_null_lengths_and_empty_sequence_evidence_are_not_invented() {
    let mut value = fixture();
    let first = &mut value["sequences"]["public_user_full_sequences"][0];
    first.as_object_mut().unwrap().remove("genbank_api_url");
    first["length"] = Value::Null;
    first["sequence"] = json!("");
    let catalog = parse(&value).unwrap();
    assert_eq!(catalog.choices()[0].length, None);
    assert_eq!(catalog.choices()[0].catalog_sequence_blake3, None);
    let report = catalog.import_download(41, GB.as_bytes(), NOW).unwrap();
    assert!(
        report
            .warnings
            .iter()
            .any(|w| w.contains("could not be cross-checked"))
    );
    let mut empty = fixture();
    for list in empty["sequences"].as_object_mut().unwrap().values_mut() {
        *list = json!([]);
    }
    assert!(parse(&empty).unwrap().choices().is_empty());
}
#[test]
fn provider_link_validation_and_private_token_configuration_are_strict() {
    for url in [
        "http://files.example.org/41.gb",
        "https://127.0.0.1/41.gb",
        "https://files.example.org:8080/41.gb",
        "https://user:password@files.example.org/41.gb",
        "https://api.developers.addgene.org/auth/test/",
    ] {
        let mut value = fixture();
        value["sequences"]["public_user_full_sequences"][0]["genbank_api_url"] = json!(url);
        assert!(parse(&value).is_err(), "{url}");
    }
    for token in ["", "contains space", "header\r\ninjection", "é"] {
        assert!(Config::default().with_addgene_api_token(token).is_err());
    }
    assert!(
        Config::default()
            .with_addgene_api_token("INDEPENDENT_TEST_TOKEN")
            .is_ok()
    );
}

#[test]
fn immutable_snapshot_keeps_the_exact_safe_summary_json_shape() {
    let raw = fixture();
    let bytes = serde_json::to_vec(&raw).unwrap();
    let catalog = parse_addgene_catalog(12345, &bytes).unwrap();
    assert_eq!(catalog.name(), "Independent Addgene plasmid");
    let sequence_hash = blake3::hash(b"AAAACCCCGGGGTTTT").to_hex().to_string();
    let response_hash = blake3::hash(&bytes).to_hex().to_string();
    let expected = json!({
        "plasmid_id":12345,
        "name":"Independent Addgene plasmid",
        "response_blake3":response_hash,
        "choices":[
            {"sequence_id":41,"description":"Independent sequence 41","origin":"depositor","completeness":"full","length":16,"catalog_sequence_blake3":sequence_hash},
            {"sequence_id":42,"description":"Independent sequence 42","origin":"addgene","completeness":"full","length":16,"catalog_sequence_blake3":sequence_hash},
            {"sequence_id":43,"description":"Independent sequence 43","origin":"depositor","completeness":"partial","length":16,"catalog_sequence_blake3":sequence_hash},
            {"sequence_id":44,"description":"Independent sequence 44","origin":"addgene","completeness":"partial","length":16,"catalog_sequence_blake3":sequence_hash}
        ]
    });
    let summary = serde_json::to_value(&catalog).unwrap();
    assert_eq!(summary, expected);
    assert_eq!(
        catalog.response_blake3(),
        expected["response_blake3"].as_str().unwrap()
    );
    let encoded = serde_json::to_string(&catalog).unwrap();
    assert!(!encoded.contains("TRANSIENT_SECRET"));
    assert!(!encoded.contains("download_url"));
    assert!(!encoded.contains("records"));
    assert!(!encoded.contains("files.example.org"));
}

#[test]
fn detached_display_choices_and_summary_cannot_relabel_the_private_import_binding() {
    let raw = fixture();
    let catalog = parse(&raw).unwrap();
    let original_summary = serde_json::to_value(&catalog).unwrap();
    let cloned_catalog = catalog.clone();
    let mut display_choices = catalog.choices().to_vec();
    display_choices[0].sequence_id = 999;
    display_choices[0].description = "Unrelated display label".into();
    display_choices[0].origin = AddgeneSequenceOrigin::Addgene;
    display_choices[0].completeness = AddgeneCompleteness::Partial;
    display_choices[0].length = Some(1);
    display_choices[0].catalog_sequence_blake3 = Some("0".repeat(64));
    display_choices.reverse();
    display_choices.clear();
    let mut editable_summary = original_summary.clone();
    editable_summary["plasmid_id"] = json!(999);
    editable_summary["name"] = json!("Unrelated plasmid");
    editable_summary["response_blake3"] = json!("0".repeat(64));
    editable_summary["choices"] = json!([]);
    assert_ne!(editable_summary, original_summary);
    assert_eq!(serde_json::to_value(&catalog).unwrap(), original_summary);
    assert_eq!(
        serde_json::to_value(&cloned_catalog).unwrap(),
        original_summary
    );
    assert_eq!(catalog.plasmid_id(), 12345);
    assert_eq!(catalog.choices().len(), 4);
    assert_eq!(
        catalog.download_url(999),
        Err(SourceError::AddgeneUnknownSequence)
    );
    assert_eq!(
        catalog
            .import_download(999, GB.as_bytes(), NOW)
            .unwrap_err(),
        SourceError::AddgeneUnknownSequence
    );

    let report = cloned_catalog
        .import_download(41, GB.as_bytes(), NOW)
        .unwrap();
    let molecule = &report.molecules[0];
    assert_eq!(molecule.metadata["source:identifier"], "12345");
    assert_eq!(
        molecule.metadata["source:canonical_url"],
        "https://www.addgene.org/12345/"
    );
    assert_eq!(molecule.metadata["source:sequence_id"], "41");
    assert_eq!(molecule.metadata["source:addgene_origin"], "depositor");
    assert_eq!(molecule.metadata["source:addgene_completeness"], "full");
    assert_eq!(
        molecule.metadata["source:addgene_catalog_blake3"],
        catalog.response_blake3()
    );
    let evidence: Value =
        serde_json::from_str(&molecule.metadata["source:addgene_catalog_evidence"]).unwrap();
    assert_eq!(evidence["plasmid_id"], 12345);
    assert_eq!(evidence["name"], "Independent Addgene plasmid");
    assert_eq!(evidence["choice"], original_summary["choices"][0]);
    assert_eq!(
        molecule.metadata["source:content_blake3"],
        blake3::hash(GB.as_bytes()).to_hex().as_str()
    );
    assert_eq!(molecule.original_genbank.as_deref(), Some(GB));
    assert_eq!(
        catalog.import_download(43, GB.as_bytes(), NOW).unwrap_err(),
        SourceError::AddgeneSequenceMismatch
    );
    assert_eq!(
        catalog
            .import_download(
                41,
                GB.replace("aaaaccccgg ggtttt", "aaaaccccgg ggttta")
                    .as_bytes(),
                NOW
            )
            .unwrap_err(),
        SourceError::AddgeneSequenceMismatch
    );
}
