use collagene_core::formats::sequence as formats;
use collagene_model::{Molecule, Topology};
#[test]
fn annotated_exports_retain_full_supported_codon_evidence_without_imported_authority() {
    let mut m = Molecule::dna("designed", "ATGAAGTGGTAA", Topology::Linear).unwrap();
    for (k, v) in [
        ("engine", collagene_core::sequence::codon_design::CODON_ENGINE_VERSION.into()),
        ("source_revision", "a".repeat(64)),
        ("options_revision", "b".repeat(64)),
        ("required_protein", "MKW*".into()),
        (
            "options",
            "{\"genetic_code\":1,\"preferences\":{\"Explicit\":{\"weights\":[{\"codon\":\"AAG\",\"weight\":100}]}}}".into(),
        ),
    ] {
        m.metadata.insert(format!("codon_design:{k}"), v);
    }
    let comments = formats::provenance_comments(&m);
    assert!(
        comments
            .iter()
            .any(|c| c.contains("unverified producer evidence"))
    );
    assert!(comments.iter().any(|c| c.contains("required_protein=MKW*")));
    for format in [formats::TextFormat::GenBank, formats::TextFormat::Embl] {
        let report = formats::export_with_report(&m, format).unwrap();
        assert!(report.text.contains("required_protein=MKW*"));
        assert!(
            report
                .warnings
                .iter()
                .any(|w| w.contains("does not restore"))
        );
        let loaded = formats::import_text(&report.text)
            .unwrap()
            .molecules
            .remove(0);
        assert_eq!(loaded.sequence, m.sequence);
        assert!(!loaded.metadata.contains_key("codon_design:engine"));
    }
    m.metadata
        .insert("codon_design:required_protein".into(), "W".repeat(10_000));
    m.metadata
        .insert("codon_design:options".into(), "A".repeat(32 * 1024));
    let comments = formats::provenance_comments(&m);
    assert!(comments.iter().any(|c| c.ends_with(&"W".repeat(10_000))));
    assert!(comments.iter().any(|c| c.ends_with(&"A".repeat(32 * 1024))));
    assert!(!comments.iter().any(|c| c.contains("[truncated")));
}
#[test]
fn oversized_manual_metadata_is_explicitly_truncated_and_cannot_inject_control_lines() {
    let mut m = Molecule::dna("untrusted metadata", "ACGT", Topology::Linear).unwrap();
    m.metadata.insert(
        "codon_design:engine".into(),
        format!("{}\n\r\0injected", "X".repeat(129)),
    );
    let comments = formats::provenance_comments(&m);
    let c = comments
        .iter()
        .find(|c| c.contains("codon_design:engine="))
        .unwrap();
    assert!(c.contains("[truncated; retain native evidence]"));
    assert!(!c.chars().any(char::is_control));
    assert!(c.len() < 256);
}
