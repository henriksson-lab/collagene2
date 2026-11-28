//! Literal INSDC T-encoded RNA fixtures, separate from native U sequence state.
use collagene_core::{formats::sequence::*, protocol};
use collagene_model::{Alphabet, Molecule, Sequence, Topology};
const RNA: &str = "LOCUS       rna_fixture               12 bp    mRNA    linear   SYN 26-NOV-2026\nDEFINITION  Independently constructed RNA encoding fixture.\nACCESSION   ZZ000001\nVERSION     ZZ000001.1\nREFERENCE   1  (bases 3 to 8; 10 to 12)\n  AUTHORS   Example,A.\n  TITLE     Synthetic fixture citation\n  JOURNAL   Unpublished\nFEATURES             Location/Qualifiers\n     misc_feature    complement(join(1..3,7..9))\n                     /label=\"literal reverse RNA\"\nORIGIN\n        1 atggcttaan ry\n//\n";
fn import(text: &str) -> Molecule {
    import_genbank(text).unwrap().molecules.remove(0)
}
fn record(text: &str) -> gb_io::seq::Seq {
    gb_io::reader::SeqReader::new(text.as_bytes())
        .next()
        .unwrap()
        .unwrap()
}
fn calls(text: &str) -> Vec<u8> {
    record(text)
        .seq
        .iter()
        .map(|b| b.to_ascii_uppercase())
        .collect()
}
#[test]
fn standard_t_encoded_rna_keeps_original_bytes_and_native_uracil_with_exact_feature_order() {
    let report = import_genbank(RNA).unwrap();
    let m = &report.molecules[0];
    assert_eq!(m.sequence.alphabet(), Alphabet::Rna);
    assert_eq!(m.sequence.as_str(), "AUGGCUUAANRY");
    assert_eq!(m.metadata["genbank.molecule_type"], "mRNA");
    assert_eq!(m.original_genbank.as_deref(), Some(RNA));
    assert_eq!(
        m.features[0].location.extract(&m.sequence).unwrap(),
        "UUACAU"
    );
    assert!(report.warnings.iter().any(|w| w.contains("T to U")));
    let native = protocol::save_document(m).unwrap();
    assert_eq!(protocol::load_document(&native).unwrap(), *m);
    let output = export_genbank(m).unwrap();
    let rec = record(&output);
    assert_eq!(calls(&output), b"ATGGCTTAANRY");
    assert_eq!(rec.molecule_type.as_deref(), Some("mRNA"));
    assert_eq!(rec.accession.as_deref(), Some("ZZ000001"));
    assert_eq!(rec.version.as_deref(), Some("ZZ000001.1"));
    assert_eq!(rec.references[0].description, "1  (bases 3 to 8; 10 to 12)");
    assert!(rec.references[0].remark.is_none());
    let round = import(&output);
    assert_eq!(round.sequence, m.sequence);
    assert_eq!(round.features[0].location, m.features[0].location);
    assert_eq!(m.original_genbank.as_deref(), Some(RNA));
}
#[test]
fn synthetic_rna_export_uses_insdc_thymine_without_mutating_the_internal_sequence() {
    let mut m = Molecule::dna("synthetic_rna", "N", Topology::Linear).unwrap();
    m.sequence = Sequence::new(Alphabet::Rna, "AUGGCUUAANRY").unwrap();
    let before = m.clone();
    let report = export_with_report(&m, TextFormat::GenBank).unwrap();
    assert_eq!(calls(&report.text), b"ATGGCTTAANRY");
    assert_eq!(record(&report.text).molecule_type.as_deref(), Some("RNA"));
    assert!(
        report
            .warnings
            .iter()
            .any(|w| w.contains("INSDC T encoding"))
    );
    assert_eq!(import(&report.text).sequence, m.sequence);
    assert_eq!(m, before);
}
#[test]
fn rna_edit_replay_uses_normalized_baseline_and_remaps_citations_without_false_source_identity() {
    let m = import(RNA);
    let edited = m.replace(0, 0, "UU").unwrap();
    let output = export_genbank(&edited).unwrap();
    let rec = record(&output);
    assert_eq!(edited.sequence.as_str(), "UUAUGGCUUAANRY");
    assert_eq!(calls(&output), b"TTATGGCTTAANRY");
    assert_eq!(
        rec.references[0].description,
        "1  (bases 5 to 10; 12 to 14)"
    );
    let remark = rec.references[0]
        .remark
        .as_deref()
        .unwrap()
        .split_whitespace()
        .collect::<Vec<_>>()
        .join(" ");
    assert!(remark.contains("citation support for edited RNA is unverified"));
    assert!(rec.accession.is_none() && rec.version.is_none());
    assert_eq!(edited.original_genbank.as_deref(), Some(RNA));
    assert_eq!(import(&output).sequence, edited.sequence);
    assert_eq!(m, import(RNA));
    let mut stale = edited.clone();
    stale.original_genbank = Some(RNA.replace("atggcttaan", "atggcttgan"));
    assert!(export_genbank(&stale).is_err());
    let mut bypass = m.clone();
    bypass.sequence = Sequence::new(Alphabet::Rna, "UUAUGGCUUAANRY").unwrap();
    assert!(export_genbank(&bypass).is_err());
}
#[test]
fn only_rna_declared_records_normalize_and_multi_record_sources_keep_their_own_originals() {
    let dna = RNA
        .replace("mRNA    ", "DNA     ")
        .replace("rna_fixture", "dna_fixture");
    let joined = format!("{RNA}{dna}");
    let report = import_genbank(&joined).unwrap();
    assert_eq!(report.molecules.len(), 2);
    assert_eq!(report.molecules[0].sequence.as_str(), "AUGGCUUAANRY");
    assert_eq!(report.molecules[1].sequence.as_str(), "ATGGCTTAANRY");
    assert_eq!(report.molecules[0].original_genbank.as_deref(), Some(RNA));
    assert_eq!(
        report.molecules[1].original_genbank.as_deref(),
        Some(dna.as_str())
    );
    assert_eq!(
        record(&export_genbank(&report.molecules[1]).unwrap())
            .accession
            .as_deref(),
        Some("ZZ000001")
    );
    assert!(import_genbank(&dna.replace("atggcttaan", "auggcuuaan")).is_err());
}
#[test]
fn legacy_u_encoded_rna_remains_readable_and_canonical_export_does_not_count_encoding_as_an_edit() {
    let legacy = RNA.replace("atggcttaan", "auggcuuaan");
    let m = import(&legacy);
    let output = export_genbank(&m).unwrap();
    assert_eq!(m.original_genbank.as_deref(), Some(legacy.as_str()));
    assert_eq!(calls(&output), b"ATGGCTTAANRY");
    assert_eq!(record(&output).accession.as_deref(), Some("ZZ000001"));
    assert!(record(&output).references[0].remark.is_none());
    assert_eq!(import(&output).sequence, m.sequence);
}
