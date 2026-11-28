use collagene_core::formats::{self, TextFormat};
use collagene_model::{Alphabet, Molecule, Sequence, Topology};
fn molecule(alphabet: Alphabet, name: &str, calls: &str) -> Molecule {
    let mut m = Molecule::dna(name, "N", Topology::Linear).unwrap();
    m.sequence = Sequence::new(alphabet, calls).unwrap();
    m
}
#[test]
fn conflicting_literal_name_tokens_cannot_change_the_actual_alphabet_after_export() {
    for (alphabet, name, calls, marker) in [
        (
            Alphabet::Dna,
            "DNA named [molecule=protein]",
            "ACGT",
            "[molecule=DNA]",
        ),
        (
            Alphabet::Rna,
            "RNA named [molecule=protein] [molecule=DNA]",
            "AUGC",
            "[molecule=RNA]",
        ),
        (
            Alphabet::Protein,
            "Protein named [molecule=RNA]",
            "ACGT",
            "[molecule=protein]",
        ),
    ] {
        let source = molecule(alphabet, name, calls);
        let before = source.revision();
        let report = formats::export_with_report(&source, TextFormat::Fasta).unwrap();
        assert!(
            report
                .warnings
                .iter()
                .any(|w| w.contains("Conflicting alphabet markers"))
        );
        assert!(report.text.contains("; label]"));
        assert_eq!(report.text.matches(marker).count(), 1);
        let imported = formats::import_fasta(&report.text)
            .unwrap()
            .molecules
            .remove(0);
        assert_eq!(imported.sequence, source.sequence);
        assert_eq!(source.name, name);
        assert_eq!(source.revision(), before);
        let mut current = report.text;
        for _ in 0..4 {
            let next = formats::import_fasta(&current).unwrap().molecules.remove(0);
            let exported = formats::export_fasta(&next);
            assert_eq!(exported, current);
            current = exported;
        }
    }
}
#[test]
fn correct_canonical_headers_and_established_unmarked_dna_remain_unchanged() {
    for (alphabet, name, calls) in [
        (Alphabet::Dna, "DNA", "ACGT"),
        (Alphabet::Dna, "DNA [molecule=DNA]", "ACGT"),
        (Alphabet::Rna, "RNA [molecule=RNA]", "AUGC"),
        (Alphabet::Protein, "Protein [molecule=protein]", "ACGT"),
    ] {
        let source = molecule(alphabet, name, calls);
        let report = formats::export_with_report(&source, TextFormat::Fasta).unwrap();
        assert_eq!(report.text, format!(">{name}\n{calls}\n"));
        assert!(
            !report
                .warnings
                .iter()
                .any(|w| w.contains("Conflicting alphabet markers"))
        );
        let imported = formats::import_fasta(&report.text)
            .unwrap()
            .molecules
            .remove(0);
        assert_eq!(imported.name, name);
        assert_eq!(imported.sequence, source.sequence);
    }
}
#[test]
fn contradictory_header_declarations_and_profile_mismatches_refuse_atomically() {
    for header in [
        "x [molecule=DNA] [molecule=protein]",
        "x [molecule=RNA] [molecule=protein]",
        "x [molecule=RNA] [molecule=DNA]",
    ] {
        let text = format!(">{header}\nACGN\n");
        assert!(formats::import_fasta(&text).is_err());
        assert!(formats::import_fasta_with_alphabet(&text, Alphabet::Protein).is_err());
    }
    assert!(
        formats::import_fasta_with_alphabet(">DNA [molecule=DNA]\nACGT\n", Alphabet::Protein)
            .is_err()
    );
    assert!(
        formats::import_fasta_with_alphabet(">protein [molecule=protein]\nACGT\n", Alphabet::Dna)
            .is_err()
    );
    let mixed = ">DNA [molecule=DNA]\nACGT\n>protein [molecule=protein]\nACGT\n";
    let report = formats::import_fasta(mixed).unwrap();
    assert_eq!(report.molecules[0].sequence.alphabet(), Alphabet::Dna);
    assert_eq!(report.molecules[1].sequence.alphabet(), Alphabet::Protein);
    assert!(formats::import_fasta_with_alphabet(mixed, Alphabet::Protein).is_err());
}
#[test]
fn explicit_unmarked_protein_and_rna_profiles_do_not_guess_from_letters_or_names() {
    let text = ">sp|EXAMPLE|PROTEIN An explicit protein source\nACGT\n";
    assert_eq!(
        formats::import_fasta(text).unwrap().molecules[0]
            .sequence
            .alphabet(),
        Alphabet::Dna
    );
    let report = formats::import_fasta_with_alphabet(text, Alphabet::Protein).unwrap();
    assert_eq!(report.molecules[0].sequence.alphabet(), Alphabet::Protein);
    assert_eq!(
        report.molecules[0].name,
        "sp|EXAMPLE|PROTEIN An explicit protein source"
    );
    assert!(
        report
            .warnings
            .iter()
            .any(|w| w.contains("explicitly selected as Protein"))
    );
    let text = ">RNA without a private tag\nAUGCN\n";
    assert!(formats::import_fasta(text).is_err());
    let rna = formats::import_fasta_with_alphabet(text, Alphabet::Rna).unwrap();
    assert_eq!(rna.molecules[0].sequence.as_str(), "AUGCN");
    assert_eq!(rna.molecules[0].sequence.alphabet(), Alphabet::Rna);
    let exported = formats::export_fasta(&rna.molecules[0]);
    assert_eq!(
        formats::import_fasta(&exported).unwrap().molecules[0].sequence,
        rna.molecules[0].sequence
    );
}
#[test]
fn escaped_label_text_is_not_reinterpreted_and_multirecord_profiles_keep_bounds() {
    let text = ">literal [molecule=protein; label] [molecule=DNA]\nACGT\n";
    let r = formats::import_fasta(text).unwrap();
    assert_eq!(r.molecules[0].sequence.alphabet(), Alphabet::Dna);
    assert_eq!(formats::export_fasta(&r.molecules[0]), text);
    assert!(
        formats::import_fasta_with_alphabet(">empty\n>next\nACGT\n", Alphabet::Protein).is_err()
    );
    assert!(formats::import_fasta_with_alphabet("ACGT\n>x\nACGT\n", Alphabet::Protein).is_err());
}
