use collagene_antibody_numbering::{
    Insertion,
    imgt_loops::{self, Kind},
    references::human_references,
};
#[test]
fn published_literal_short_and_full_loop_placements_are_not_cdr3_gap_reuse() {
    for (length, expected) in [
        (5, vec![27, 28, 29, 37, 38]),
        (6, vec![27, 28, 29, 36, 37, 38]),
        (7, vec![27, 28, 29, 30, 36, 37, 38]),
        (8, vec![27, 28, 29, 30, 35, 36, 37, 38]),
        (9, vec![27, 28, 29, 30, 31, 35, 36, 37, 38]),
        (10, vec![27, 28, 29, 30, 31, 34, 35, 36, 37, 38]),
        (11, vec![27, 28, 29, 30, 31, 32, 34, 35, 36, 37, 38]),
        (12, (27..=38).collect()),
    ] {
        assert_eq!(
            imgt_loops::positions(Kind::Cdr1, length)
                .unwrap()
                .iter()
                .map(|p| p.number)
                .collect::<Vec<_>>(),
            expected
        );
    }
    for (length, expected) in [
        (0, vec![]),
        (1, vec![56]),
        (2, vec![56, 65]),
        (3, vec![56, 57, 65]),
        (4, vec![56, 57, 64, 65]),
        (5, vec![56, 57, 58, 64, 65]),
        (6, vec![56, 57, 58, 63, 64, 65]),
        (7, vec![56, 57, 58, 59, 63, 64, 65]),
        (8, vec![56, 57, 58, 59, 62, 63, 64, 65]),
        (9, vec![56, 57, 58, 59, 60, 62, 63, 64, 65]),
        (10, (56..=65).collect()),
    ] {
        assert_eq!(
            imgt_loops::positions(Kind::Cdr2, length)
                .unwrap()
                .iter()
                .map(|p| p.number)
                .collect::<Vec<_>>(),
            expected
        );
    }
}
#[test]
fn every_provider_v_reference_agrees_with_the_independent_loop_allocator() {
    let corpus = human_references().unwrap();
    for reference in corpus.v_entries() {
        for (kind, low, high) in [(Kind::Cdr1, 27, 38), (Kind::Cdr2, 56, 65)] {
            let observed = reference
                .alignment_labels
                .iter()
                .zip(reference.aligned_protein.bytes())
                .filter_map(|(label, aa)| {
                    let n = label.parse::<u16>().unwrap();
                    (aa != b'-' && (low..=high).contains(&n)).then_some(n)
                })
                .collect::<Vec<_>>();
            let allocated = imgt_loops::positions(kind, observed.len()).unwrap();
            assert!(allocated.iter().all(|p| p.insertion == Insertion::None));
            assert_eq!(
                allocated.iter().map(|p| p.number).collect::<Vec<_>>(),
                observed,
                "{}",
                reference.label
            );
        }
    }
}
#[test]
fn unsupported_lengths_are_refused_instead_of_inventing_extra_loop_positions() {
    for length in [0, 1, 4, 13, usize::MAX] {
        assert!(imgt_loops::positions(Kind::Cdr1, length).is_err());
    }
    for length in [11, usize::MAX] {
        assert!(imgt_loops::positions(Kind::Cdr2, length).is_err());
    }
    assert_eq!(
        imgt_loops::catalog_revision(),
        imgt_loops::catalog_revision()
    );
}
