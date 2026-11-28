use collagene_core::formats::sequence as formats;
use collagene_model::{Molecule, Topology};
#[test]
fn annotated_contig_exports_keep_bounded_producer_claims_without_read_evidence_authority() {
    let mut m = Molecule::dna("inferred consensus", "ACGTNCGT", Topology::Linear).unwrap();
    for (k, v) in [
        ("engine", "contig-1".into()),
        ("batch_revision", "a".repeat(64)),
        ("options_revision", "b".repeat(64)),
        ("evidence_hash", "c".repeat(64)),
        ("options", "{\"min_overlap\":20}".into()),
        ("reads", "[[\"read_1\",\"source-revision\"]]".into()),
        ("id", "0".into()),
        ("complete", "false".into()),
    ] {
        m.metadata.insert(format!("contig:{k}"), v);
    }
    for format in [formats::TextFormat::GenBank, formats::TextFormat::Embl] {
        let report = formats::export_with_report(&m, format).unwrap();
        assert!(report.text.contains("contig:engine=contig-1"));
        assert!(report.text.contains("unverified producer evidence"));
        assert!(report.text.contains("contig:complete=false"));
        assert!(
            report
                .warnings
                .iter()
                .any(|v| v.contains("not a lossless evidence interchange"))
        );
        let loaded = formats::import_text(&report.text)
            .unwrap()
            .molecules
            .remove(0);
        assert_eq!(loaded.sequence, m.sequence);
        assert!(!loaded.metadata.contains_key("contig:complete"));
    }
}
#[test]
fn forged_oversized_contig_comments_are_explicitly_truncated_and_control_safe() {
    let mut m = Molecule::dna("untrusted evidence", "ACGT", Topology::Linear).unwrap();
    m.metadata.insert(
        "contig:reads".into(),
        format!("{}\n\r\0ignored", "X".repeat(32 * 1024 + 1)),
    );
    let comments = formats::provenance_comments(&m);
    let row = comments
        .iter()
        .find(|x| x.contains("contig:reads="))
        .unwrap();
    assert!(row.contains("[truncated; retain native evidence]"));
    assert!(!row.chars().any(char::is_control));
    assert!(row.len() < 33 * 1024);
}
