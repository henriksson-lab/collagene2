use collagene_antibody_numbering::references::human_v_references;
use collagene_core::{
    Alphabet, Sequence,
    translation::{TranslationOptions, translate},
};
#[test]
fn every_reference_translation_and_cdr_boundary_agrees_with_the_existing_scientific_core() {
    let corpus = human_v_references().unwrap();
    let mut codons = 0;
    for e in corpus.entries() {
        let sequence =
            Sequence::new(Alphabet::Dna, &e.aligned_nucleotides.replace('.', "")).unwrap();
        let report = translate(&sequence, &TranslationOptions::default()).unwrap();
        assert_eq!(
            report.protein,
            e.aligned_protein.replace('-', ""),
            "{}",
            e.allele_ref
        );
        assert_eq!(report.trailing_bases, 0);
        assert_eq!(report.ambiguous_codons, 0);
        assert!(!report.protein.contains('*'));
        codons += report.translated_codons;
        let numbered = e
            .aligned_protein
            .bytes()
            .zip(&e.alignment_labels)
            .filter(|(aa, _)| *aa != b'-')
            .map(|(_, label)| label.parse::<usize>().unwrap())
            .collect::<Vec<_>>();
        for cdr in &e.cdrs {
            let (low, high) = if cdr.kind == 1 { (27, 38) } else { (56, 65) };
            assert!(
                numbered[cdr.start..cdr.end]
                    .iter()
                    .all(|n| (low..=high).contains(n))
            );
            assert!(numbered[..cdr.start].iter().all(|n| *n < low));
            assert!(numbered[cdr.end..].iter().all(|n| *n > high));
        }
    }
    assert!(codons > 30_000);
}
