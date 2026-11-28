use collagene_core::formats::*;
use collagene_model::Alphabet;

#[test]
fn format_sniffing_preserves_text_record_terminators_and_provenance_exports() {
    let text = "ID   source_fixture; SV 1; linear; genomic DNA; STD; SYN; 4 BP.\nXX\nSQ   Sequence 4 BP; 1 A; 1 C; 1 G; 1 T; 0 other;\n     acgt       4\n//\n";
    let mol = import_bytes(text.as_bytes()).unwrap().molecules.remove(0);
    assert_eq!(collagene_core::embl::export_embl(&mol).unwrap(), text);
    let mut current = mol;
    current.metadata.insert(
        "source:canonical_url".into(),
        "https://www.ncbi.nlm.nih.gov/nuccore/SYN000001.1".into(),
    );
    current
        .metadata
        .insert("source:retrieved_at".into(), "2026-11-23T12:00:00Z".into());
    let gb = export_with_report(&current, TextFormat::GenBank).unwrap();
    assert!(!gb.warnings.is_empty());
    assert!(gb.text.contains("source:canonical_url="));
    assert!(gb.text.contains("source:retrieved_at="));
    let parsed = import_genbank(&gb.text).unwrap().molecules.remove(0);
    let exported = export_genbank(&parsed).unwrap();
    assert!(exported.contains("source:canonical_url="));
    let embl = export_with_report(&current, TextFormat::Embl).unwrap();
    assert!(
        embl.text
            .contains("CC   Collagene provenance source:canonical_url=")
    );
    assert!(
        !export_with_report(&current, TextFormat::Fasta)
            .unwrap()
            .warnings
            .is_empty()
    );
    assert!(import_bytes(b"ABIF").is_err());
    assert!(import_bytes(&[0xff, 0, 0xff]).is_err());
    assert_eq!(
        import_bytes(b">dna\nACGT\n").unwrap().molecules[0]
            .sequence
            .as_str(),
        "ACGT"
    );
    assert_eq!(current.sequence.alphabet(), Alphabet::Dna);
}

#[test]
fn explicit_single_stranded_genbank_is_retained_and_not_invented_as_a_duplex() {
    let text = "LOCUS       ss_fixture                 8 bp    ss-DNA  circular SYN 01-JAN-2026\nFEATURES             Location/Qualifiers\nORIGIN\n        1 gaattcaa\n//\n";
    let molecule = import_text(text).unwrap().molecules.remove(0);
    assert_eq!(molecule.metadata["genbank.molecule_type"], "ss-DNA");
    assert!(collagene_core::chemistry::Duplex::from_molecule(&molecule).is_err());
    let exported = export_genbank(&molecule).unwrap();
    assert!(exported.lines().next().unwrap().contains("ss-DNA"));
    let restored = import_text(&exported).unwrap().molecules.remove(0);
    assert!(collagene_core::chemistry::has_unmapped_physical_evidence(
        &restored
    ));
}
