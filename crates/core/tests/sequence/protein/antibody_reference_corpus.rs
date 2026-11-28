use collagene_antibody_numbering::references::{Species, load};
use collagene_core::sequence::translation::{TranslationOptions, translate};
use collagene_model::{Alphabet, Sequence};
#[test]
fn every_reference_translation_and_cdr_boundary_agrees_with_the_existing_scientific_core() {
    let corpus = load(&Species::Human).unwrap();
    let mut codons = 0;
    for e in corpus.v_entries() {
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

#[test]
fn every_j_fr4_translation_matches_the_core_and_preserves_the_untranslated_terminal_base() {
    let corpus = load(&Species::Human).unwrap();
    for entry in corpus.j_entries() {
        let dna = Sequence::new(Alphabet::Dna, &entry.fr4_nucleotides).unwrap();
        let report = translate(&dna, &TranslationOptions::default()).unwrap();
        assert_eq!(report.protein, entry.fr4_protein, "{}", entry.allele_ref);
        assert_eq!(report.trailing_bases, 0);
        let full_tail = Sequence::new(
            Alphabet::Dna,
            &entry.coding_nucleotides[entry.anchor_offset..],
        )
        .unwrap();
        let full = translate(&full_tail, &TranslationOptions::default()).unwrap();
        assert_eq!(full.protein, entry.fr4_protein);
        assert_eq!(full.trailing_bases, 1);
        assert_eq!(entry.fr4_positions.len(), full.translated_codons);
    }
}

#[test]
fn mouse_captured_v_and_authored_fr4_segments_translate_with_the_existing_core() {
    let corpus = load(&Species::MouseC57Bl6).unwrap();
    assert_eq!(corpus.species(), &Species::MouseC57Bl6);
    assert_eq!(corpus.v_entries().len(), 95);
    assert_eq!(corpus.j_entries().len(), 4);
    for entry in corpus.v_entries() {
        let dna = Sequence::dna(&entry.aligned_nucleotides.replace('.', "")).unwrap();
        let translated = translate(&dna, &TranslationOptions::default()).unwrap();
        assert_eq!(
            translated.protein,
            entry.aligned_protein.replace('-', ""),
            "{}",
            entry.allele_ref
        );
        assert_eq!(translated.trailing_bases, 0);
        assert_eq!(translated.ambiguous_codons, 0);
    }
    for entry in corpus.j_entries() {
        let dna = Sequence::dna(&entry.fr4_nucleotides).unwrap();
        let translated = translate(&dna, &TranslationOptions::default()).unwrap();
        assert_eq!(
            translated.protein, entry.fr4_protein,
            "{}",
            entry.allele_ref
        );
        assert_eq!(translated.trailing_bases, 0);
        // These are explicitly authored rule labels; raw provider frame metadata
        // is not reinterpreted as provider-supplied IMGT J numbering.
        assert_eq!(
            entry.fr4_positions,
            (118..=128).map(|n| n.to_string()).collect::<Vec<_>>()
        );
    }
}
