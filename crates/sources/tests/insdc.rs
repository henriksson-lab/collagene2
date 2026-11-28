//! Independent fixtures; no provider network access or other application code.
use collagene_core::{formats::sequence as formats, protocol};
use collagene_sources::*;
const NOW: &str = "2026-11-24T12:00:00Z";
const GB: &str = "LOCUS       AB123456                  16 bp    DNA     circular SYN 01-JAN-2026\nACCESSION   AB123456\nVERSION     AB123456.2\nFEATURES             Location/Qualifiers\n     misc_feature    complement(join(13..16,1..4))\n                     /label=\"origin\"\nORIGIN\n        1 aaaaccccgg ggtttt\n//\n";
const EMBL: &str = "ID   AB123456; SV 2; circular; genomic DNA; STD; SYN; 16 BP.\nAC   AB123456;\nDE   Independently authored fixture.\nFH   Key             Location/Qualifiers\nFH\nFT   misc_feature    complement(join(13..16,1..4))\nFT                   /label=\"origin\"\nSQ   Sequence 16 BP; 4 A; 4 C; 4 G; 4 T; 0 other;\n     aaaaccccgg ggtttt       16\n//\n";
#[test]
fn exact_documented_routes_and_explicit_provider_identity() {
    for input in [
        "ena:ab123456.2",
        "ENA:AB123456.2",
        "https://www.ebi.ac.uk/ena/browser/view/AB123456.2",
        "https://www.ebi.ac.uk/ena/browser/api/fasta/AB123456.2?download=true",
    ] {
        let source = parse_source(input).unwrap();
        assert_eq!(
            source,
            Source::Ena {
                accession: "AB123456.2".into()
            }
        );
        assert_eq!(
            source.download_url().unwrap(),
            "https://www.ebi.ac.uk/ena/browser/api/embl/AB123456.2"
        );
    }
    for input in [
        "ddbj:ab123456.2",
        "DDBJ:AB123456.2",
        "https://getentry.ddbj.nig.ac.jp/getentry/na/AB123456.2/",
        "https://getentry.ddbj.nig.ac.jp/getentry?database=ddbj&accession_number=AB123456.2&format=fasta&filetype=html&limit=1",
    ] {
        let source = parse_source(input).unwrap();
        assert_eq!(
            source,
            Source::Ddbj {
                accession: "AB123456.2".into()
            }
        );
        assert_eq!(
            source.download_url().unwrap(),
            "https://getentry.ddbj.nig.ac.jp/getentry/na/AB123456.2?format=flatfile&filetype=text"
        );
    }
    assert!(matches!(
        parse_source("AB123456.2").unwrap(),
        Source::Ncbi { .. }
    ));
    assert!(matches!(
        parse_source("12345").unwrap(),
        Source::Addgene { .. }
    ));
}
#[test]
fn selectors_collections_encoded_ids_duplicate_keys_and_aliases_fail_closed() {
    for input in [
        "ena:PRJEB12345",
        "ddbj:GCA_000001.1",
        "ena:ERR123456",
        "ena:AAAA01000000",
        "ena:12345",
        "https://www.ebi.ac.uk/ena/browser/view/AB123456.2?range=1-4",
        "https://www.ebi.ac.uk/ena/browser/view/AB123456.2?download=true&download=false",
        "https://www.ebi.ac.uk/ena/browser/api/embl/AB123456.2:1-4",
        "https://www.ebi.ac.uk/ena/browser/api/embl/AB123456%2E2",
        "https://www.ebi.ac.uk/ena/browser/api/embl/search?query=plasmid",
        "https://www.ebi.ac.uk/ena/browser/view/AB123456.2,AB123457.2",
        "https://getentry.ddbj.nig.ac.jp/getentry/na/AB123456.2?revision=2020-01-01",
        "https://getentry.ddbj.nig.ac.jp/getentry/na/AB123456.2/1",
        "https://getentry.ddbj.nig.ac.jp/getentry/na/AB123456-AB123457",
        "https://getentry.ddbj.nig.ac.jp/getentry?database=aa&accession_number=AB123456.2",
        "https://getentry.ddbj.nig.ac.jp/getentry?accession_number=AB123456.2&accession_number=AB123457.2",
        "https://getentry.ddbj.nig.ac.jp/getentry/na/AB123456.2?strand=2",
        "https://getentry.ddbj.nig.ac.jp/getentry/na/AB123456.2?format=trans",
    ] {
        assert!(parse_source(input).is_err(), "accepted {input}");
    }
}
#[test]
fn both_providers_pin_actual_versions_with_exact_annotated_evidence() {
    for (prefix, raw, provider, format, evidence) in [
        ("ena", EMBL, "ENA", "EMBL", "embl.original_record"),
        ("ddbj", GB, "DDBJ", "GenBank", ""),
    ] {
        let s = parse_source(&format!("{prefix}:AB123456")).unwrap();
        let r = import_download(&s, raw.as_bytes(), NOW).unwrap();
        let m = &r.molecules[0];
        assert_eq!(m.metadata["source:provider"], provider);
        assert_eq!(m.metadata["source:format"], format);
        assert_eq!(m.metadata["source:identifier"], "AB123456");
        assert_eq!(m.metadata["source:accession_version"], "AB123456.2");
        assert!(m.metadata["source:canonical_url"].ends_with("AB123456.2"));
        assert_eq!(
            m.metadata["source:content_blake3"],
            blake3::hash(raw.as_bytes()).to_hex().as_str()
        );
        assert_eq!(
            m.features[0].location.extract(&m.sequence).unwrap(),
            "TTTTAAAA"
        );
        if evidence.is_empty() {
            assert_eq!(m.original_genbank.as_deref(), Some(GB));
        } else {
            assert_eq!(m.metadata[evidence], raw);
        }
        let saved = protocol::save_document(m).unwrap();
        assert_eq!(&protocol::load_document(&saved).unwrap(), m);
        for wrong in ["AB123456.1", "AB123457.2"] {
            assert_eq!(
                import_download(
                    &parse_source(&format!("{prefix}:{wrong}")).unwrap(),
                    raw.as_bytes(),
                    NOW
                )
                .unwrap_err(),
                SourceError::IdentifierMismatch
            );
        }
        assert!(import_download(&s, format!("{raw}{raw}").as_bytes(), NOW).is_err());
        assert!(
            import_download(&s, raw.trim_end().trim_end_matches("//").as_bytes(), NOW).is_err()
        );
        assert!(import_download(&s, b"<html>login</html>", NOW).is_err());
        assert!(import_download(&s, b">AB123456.2\nACGT\n", NOW).is_err());
    }
}
#[test]
fn legacy_locus_names_require_real_sv_and_identity_contradictions_are_rejected() {
    let source = parse_source("ena:AB123456.2").unwrap();
    let legacy = EMBL
        .replace(
            "AB123456; SV 2; circular; genomic DNA; STD; SYN; 16 BP.",
            "LOCUS_ALIAS standard; circular genomic DNA; SYN; 16 BP.",
        )
        .replace("AC   AB123456;", "AC   AB123456;\nSV   AB123456.2");
    assert_eq!(
        import_download(&source, legacy.as_bytes(), NOW)
            .unwrap()
            .molecules[0]
            .metadata["source:accession_version"],
        "AB123456.2"
    );
    for raw in [
        legacy.replace("SV   AB123456.2\n", ""),
        legacy.replace("SV   AB123456.2", "SV   AB123457.2"),
        EMBL.replace("AC   AB123456;", "AC   AB123457;"),
        EMBL.replace("AC   AB123456;", "AC   AB123456;\nSV   AB123457.2"),
    ] {
        assert!(import_download(&source, raw.as_bytes(), NOW).is_err());
    }
}
#[test]
fn portable_annotation_cannot_bind_edited_or_stale_documents_to_old_bytes() {
    for (s, raw) in [("ena:AB123456.2", EMBL), ("ddbj:AB123456.2", GB)] {
        let mut r = formats::import_bytes(raw.as_bytes()).unwrap();
        r.molecules[0].name = "edited".into();
        assert_eq!(
            annotate_import(&mut r, &parse_source(s).unwrap(), raw.as_bytes(), NOW),
            Err(SourceError::InvalidSequenceData)
        );
    }
}
