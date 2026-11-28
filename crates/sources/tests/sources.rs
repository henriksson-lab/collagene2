use base64::{Engine as _, engine::general_purpose::STANDARD};
use collagene_core::formats::sequence as formats;
use collagene_model::Topology;
use collagene_sources::*;
use url::Url;

const NOW: &str = "2026-11-23T18:00:00Z";
const GB: &str = "LOCUS       TEST000001                16 bp    DNA     circular SYN 01-JAN-2026\nDEFINITION  Independently authored source-import fixture.\nACCESSION   TEST000001\nVERSION     TEST000001.2\nFEATURES             Location/Qualifiers\n     CDS             complement(join(13..16,1..4))\n                     /label=\"origin feature\"\nORIGIN\n        1 aaaaccccgg ggtttt\n//\n";

#[test]
fn resolver_accepts_specific_ncbi_records_and_keeps_public_fetch_urls_credential_free() {
    for input in [
        "test000001.2",
        "ncbi:TEST000001.2",
        "https://www.ncbi.nlm.nih.gov/nuccore/TEST000001.2/",
        "https://www.ncbi.nlm.nih.gov/nucleotide/TEST000001.2?report=genbank",
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=TEST000001.2&api_key=DO_NOT_KEEP",
    ] {
        let source = parse_source(input).unwrap();
        assert_eq!(
            source,
            Source::Ncbi {
                accession: "TEST000001.2".into()
            }
        );
        assert_eq!(
            source.canonical_url(),
            "https://www.ncbi.nlm.nih.gov/nuccore/TEST000001.2"
        );
        let download = Url::parse(&source.download_url().unwrap()).unwrap();
        assert_eq!(download.host_str(), Some("eutils.ncbi.nlm.nih.gov"));
        let query: std::collections::BTreeMap<_, _> = download.query_pairs().collect();
        assert_eq!(query["db"], "nuccore");
        assert_eq!(query["rettype"], "gbwithparts");
        assert_eq!(query["retmode"], "text");
        assert_eq!(query["tool"], "collagene2");
        assert!(!query.contains_key("api_key"));
    }
    for input in [
        "https://www.ncbi.nlm.nih.gov/nuccore/TEST000001.2?from=2&to=4",
        "https://www.ncbi.nlm.nih.gov/nuccore/?term=plasmid",
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&id=TEST000001.2",
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=TEST000001.2,TEST000002.1",
    ] {
        assert_eq!(
            parse_source(input).unwrap_err(),
            SourceError::UnsupportedSelection
        );
    }
}

#[test]
fn addgene_pages_are_typed_access_requirements_and_signed_links_remain_transient() {
    for input in [
        "12345",
        "#12345",
        "addgene:12345",
        "https://www.addgene.org/12345/",
        "https://www.addgene.org/12345/sequences/",
    ] {
        let source = parse_source(input).unwrap();
        assert_eq!(source, Source::Addgene { plasmid_id: 12345 });
        assert_eq!(
            source.download_url().unwrap_err(),
            SourceError::AddgeneAccessRequired { plasmid_id: 12345 }
        );
    }
    let source=parse_source("https://media.addgene.org/sequences/plasmid.gb?X-Amz-Signature=SECRET&Expires=123#fragment").unwrap();
    assert!(source.download_url().unwrap().contains("SECRET"));
    assert_eq!(
        source.canonical_url(),
        "https://media.addgene.org/sequences/plasmid.gb"
    );
    assert!(!format!("{source:?}").contains("SECRET"));
}

#[test]
fn local_private_mapped_addresses_credentials_and_non_https_urls_are_rejected() {
    for input in [
        "http://example.org/file.gb",
        "file:///etc/passwd",
        "https://user:SECRET@example.org/file.gb",
        "https://localhost/file.fa",
        "https://local/file.fa",
        "https://server.internal/file.gb",
        "https://127.1/file.gb",
        "https://2130706433/file.gb",
        "https://10.0.0.1/file.gb",
        "https://169.254.169.254/file.fa",
        "https://100.64.1.2/file.fa",
        "https://198.18.1.1/file.fa",
        "https://[::1]/file.gb",
        "https://[::ffff:127.0.0.1]/file.gb",
        "https://[2002:7f00:1::]/file.gb",
        "https://[2001:db8::1]/file.gb",
        "https://example.org:8443/file.gb",
    ] {
        assert!(parse_source(input).is_err(), "accepted {input}");
    }
    for input in [
        "https://8.8.8.8/file.fa",
        "https://[2001:4860:4860::8888]/file.gb",
        "https://sequence.example.org/file.gb",
    ] {
        validate_public_url(&Url::parse(input).unwrap()).unwrap();
    }
}

#[test]
fn ncbi_provenance_pins_actual_version_and_preserves_original_annotations() {
    let source =
        parse_source("https://www.ncbi.nlm.nih.gov/nuccore/TEST000001.2?api_key=SECRET").unwrap();
    let report = import_download(&source, GB.as_bytes(), NOW).unwrap();
    let m = &report.molecules[0];
    assert_eq!(m.sequence.as_str(), "AAAACCCCGGGGTTTT");
    assert_eq!(m.topology, Topology::Circular);
    assert_eq!(
        m.features[0].location.extract(&m.sequence).unwrap(),
        "TTTTAAAA"
    );
    assert_eq!(m.metadata["source:provider"], "NCBI nuccore");
    assert_eq!(m.metadata["source:identifier"], "TEST000001.2");
    assert_eq!(m.metadata["source:accession_version"], "TEST000001.2");
    assert_eq!(
        m.metadata["source:canonical_url"],
        "https://www.ncbi.nlm.nih.gov/nuccore/TEST000001.2"
    );
    assert_eq!(m.metadata["source:retrieved_at"], NOW);
    assert_eq!(m.metadata["source:format"], "GenBank");
    assert_eq!(
        m.metadata["source:content_blake3"],
        blake3::hash(GB.as_bytes()).to_hex().as_str()
    );
    assert!(
        m.original_genbank
            .as_ref()
            .unwrap()
            .contains("VERSION     TEST000001.2")
    );
    let saved = collagene_core::protocol::save_document(m).unwrap();
    let restored = collagene_core::protocol::load_document(&saved).unwrap();
    assert_eq!(restored.metadata, m.metadata);
    assert!(!serde_json::to_string(&restored).unwrap().contains("SECRET"));
    assert_eq!(
        import_download(&parse_source("TEST000001.3").unwrap(), GB.as_bytes(), NOW).unwrap_err(),
        SourceError::IdentifierMismatch
    );
    let current =
        import_download(&parse_source("TEST000001").unwrap(), GB.as_bytes(), NOW).unwrap();
    assert_eq!(
        current.molecules[0].metadata["source:accession_version"],
        "TEST000001.2"
    );
}

#[test]
fn html_missing_versions_partial_records_and_ambiguous_ncbi_responses_fail_explicitly() {
    let source = parse_source("TEST000001.2").unwrap();
    for html in [
        b"<!DOCTYPE html><html>Login required</html>".as_slice(),
        b"Error: API key invalid",
        b"ACGT",
    ] {
        assert_eq!(
            import_download(&source, html, NOW).unwrap_err(),
            SourceError::UnexpectedContent
        );
    }
    let no_version = GB.replace("VERSION     TEST000001.2\n", "");
    assert_eq!(
        import_download(&source, no_version.as_bytes(), NOW).unwrap_err(),
        SourceError::MissingAccessionVersion
    );
    let partial = GB.replace("aaaaccccgg ggtttt", "aaaaccccgg");
    assert!(matches!(
        import_download(&source, partial.as_bytes(), NOW),
        Err(SourceError::Truncated | SourceError::InvalidSequenceData)
    ));
    let doubled = format!("{GB}{GB}");
    assert_eq!(
        import_download(&source, doubled.as_bytes(), NOW).unwrap_err(),
        SourceError::AmbiguousResponse { records: 2 }
    );
    assert!(import_download(&source, GB.as_bytes(), "not-a-time").is_err());
}

#[test]
fn direct_fasta_preserves_each_original_record_and_redacts_signed_query_provenance() {
    let source = parse_source("https://files.example.org/sequence.fa?signature=SECRET").unwrap();
    let raw = b"\xef\xbb\xbf>first\nacgt\n>second\nGGCC\n";
    let report = import_download(&source, raw, NOW).unwrap();
    assert_eq!(report.molecules.len(), 2);
    assert_eq!(
        STANDARD
            .decode(&report.molecules[0].metadata["source:original_fasta_base64"])
            .unwrap(),
        b"\xef\xbb\xbf>first\nacgt\n"
    );
    assert_eq!(
        STANDARD
            .decode(&report.molecules[1].metadata["source:original_fasta_base64"])
            .unwrap(),
        b">second\nGGCC\n"
    );
    for m in &report.molecules {
        assert_eq!(
            m.metadata["source:canonical_url"],
            "https://files.example.org/sequence.fa"
        );
        assert!(!serde_json::to_string(m).unwrap().contains("SECRET"));
    }
    assert!(
        report
            .warnings
            .iter()
            .any(|w| w.contains("biological completeness"))
    );
    let mut parsed = formats::import_bytes(raw).unwrap();
    annotate_import(&mut parsed, &source, raw, "2026-11-23T19:00:00+01:00").unwrap();
    assert_eq!(parsed.molecules[0].metadata["source:retrieved_at"], NOW);
}

#[test]
fn url_and_local_fasta_share_literal_records_despite_identical_normalized_calls() {
    let first = b"\xef\xbb\xbf\r\n  >first documentary header\r\nacgt\r\n\r\n";
    let second = b" >second\nACGT";
    let raw = [first.as_slice(), second.as_slice()].concat();
    let local = formats::import_bytes(&raw).unwrap();
    let remote = import_download(
        &parse_source("https://files.example.org/two.fa").unwrap(),
        &raw,
        NOW,
    )
    .unwrap();
    assert_eq!(local.molecules.len(), 2);
    assert_eq!(local.molecules[0].sequence, local.molecules[1].sequence);
    for (index, literal) in [first.as_slice(), second.as_slice()]
        .into_iter()
        .enumerate()
    {
        for molecule in [&local.molecules[index], &remote.molecules[index]] {
            assert_eq!(
                STANDARD
                    .decode(&molecule.metadata["source:original_fasta_base64"])
                    .unwrap(),
                literal
            );
            assert_eq!(
                molecule.metadata["source:original_fasta_blake3"],
                blake3::hash(literal).to_hex().as_str()
            );
        }
        assert_eq!(
            remote.molecules[index].metadata["source:record_index"],
            index.to_string()
        );
        assert_eq!(
            remote.molecules[index].metadata["source:content_blake3"],
            blake3::hash(&raw).to_hex().as_str()
        );
    }
}

#[test]
fn public_annotation_refuses_changed_reports_atomically_for_fasta_and_ncbi() {
    let raw = b">first\nACGT\n>second\nGGCC\n";
    let source = parse_source("https://files.example.org/two.fa").unwrap();
    let canonical = formats::import_bytes(raw).unwrap();
    let mut cases = Vec::new();
    let mut changed = canonical.clone();
    changed.molecules[0].name = "edited".into();
    cases.push(changed);
    let mut changed = canonical.clone();
    changed.molecules[0].sequence =
        collagene_model::Molecule::dna("calls", "TTTT", Topology::Linear)
            .unwrap()
            .sequence;
    cases.push(changed);
    let mut changed = canonical.clone();
    changed.molecules.swap(0, 1);
    cases.push(changed);
    let mut changed = canonical.clone();
    changed.molecules[0]
        .metadata
        .insert("source:original_fasta_blake3".into(), "a".repeat(64));
    cases.push(changed);
    let mut changed = canonical;
    changed.warnings.push("Unrelated caller assertion".into());
    cases.push(changed);
    for mut changed in cases {
        let before = serde_json::to_value(&changed).unwrap();
        assert_eq!(
            annotate_import(&mut changed, &source, raw, NOW),
            Err(SourceError::InvalidSequenceData)
        );
        assert_eq!(serde_json::to_value(changed).unwrap(), before);
    }
    let mut genbank = formats::import_bytes(GB.as_bytes()).unwrap();
    genbank.molecules[0].name = "Edited imported record".into();
    let before = serde_json::to_value(&genbank).unwrap();
    assert_eq!(
        annotate_import(
            &mut genbank,
            &parse_source("TEST000001.2").unwrap(),
            GB.as_bytes(),
            NOW
        ),
        Err(SourceError::InvalidSequenceData)
    );
    assert_eq!(serde_json::to_value(genbank).unwrap(), before);
}

#[test]
fn provider_failure_and_aggregate_decoration_refusal_leave_caller_report_unchanged() {
    let mut report = formats::import_bytes(GB.as_bytes()).unwrap();
    let before = serde_json::to_value(&report).unwrap();
    assert_eq!(
        annotate_import(
            &mut report,
            &parse_source("TEST000001.3").unwrap(),
            GB.as_bytes(),
            NOW
        ),
        Err(SourceError::IdentifierMismatch)
    );
    assert_eq!(serde_json::to_value(report).unwrap(), before);
    let raw = (0..8)
        .map(|index| format!(">record{index}\nACGT\n"))
        .collect::<String>();
    let mut report = formats::import_bytes(raw.as_bytes()).unwrap();
    let before = serde_json::to_value(&report).unwrap();
    let source = Source::Direct {
        url: format!("https://files.example.org/{}", "a".repeat(1024 * 1024)),
    };
    assert_eq!(
        annotate_import(&mut report, &source, raw.as_bytes(), NOW),
        Err(SourceError::TooLarge)
    );
    assert_eq!(serde_json::to_value(report).unwrap(), before);
}
