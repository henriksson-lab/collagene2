//! Independent bounded-work regressions from review.
use collagene_core::sequence::search::*;
use collagene_model::{Feature, Location, Molecule, Topology};
#[test]
fn empty_qualifier_rows_consume_work_and_cannot_hide_unbounded_scanning() {
    let mut molecule = Molecule::dna("bounded metadata", "ACGT", Topology::Linear).unwrap();
    molecule.features.push(Feature {
        kind: "gene".into(),
        location: Location::range(0, 4),
        qualifiers: vec![(String::new(), None); 10_000],
        needs_review: false,
    });
    let options = SearchOptions {
        query: "absent".into(),
        mode: SearchMode::Annotations,
        max_work: 32,
        ..Default::default()
    };
    let report = search_document(&molecule, &options).unwrap();
    assert!(
        !report.complete,
        "empty annotation entries still cost work; a low work budget cannot prove absence across them"
    );
    assert!(report.work <= options.max_work);
    assert!(report.hits.is_empty());
    report.validate().unwrap();
}
#[test]
fn bounded_metadata_exhaustion_is_incomplete_even_when_later_qualifier_matches() {
    let mut molecule = Molecule::dna("bounded metadata", "ACGT", Topology::Linear).unwrap();
    let mut qualifiers = vec![(String::new(), None); 100];
    qualifiers.push(("label".into(), Some("needle".into())));
    molecule.features.push(Feature {
        kind: "gene".into(),
        location: Location::range(0, 4),
        qualifiers,
        needs_review: false,
    });
    let options = SearchOptions {
        query: "needle".into(),
        mode: SearchMode::Annotations,
        max_work: 32,
        ..Default::default()
    };
    let report = search_document(&molecule, &options).unwrap();
    assert!(!report.complete);
    assert!(report.hits.is_empty());
    let full = search_document(
        &molecule,
        &SearchOptions {
            max_work: 10_000,
            ..options
        },
    )
    .unwrap();
    assert!(full.complete);
    assert_eq!(full.hits.len(), 1);
}

#[test]
fn shared_exact_kernel_keeps_overlaps_circular_origins_rna_and_literal_protein() {
    use collagene_core::sequence::analysis;
    use collagene_model::{Alphabet, Sequence};
    for (alphabet, text, query, circular, expected) in [
        (Alphabet::Dna, "AAAA", "AAA", false, vec![0, 1]),
        (Alphabet::Dna, "AAAA", "AAA", true, vec![0, 1, 2, 3]),
        (Alphabet::Dna, "ATGC", "CAT", true, vec![3]),
        (Alphabet::Rna, "UAUAU", "UAU", false, vec![0, 2]),
        (Alphabet::Protein, "MXXMX", "MX", false, vec![0, 3]),
    ] {
        let mut molecule = Molecule::dna(
            "overlaps",
            "A",
            if circular {
                Topology::Circular
            } else {
                Topology::Linear
            },
        )
        .unwrap();
        molecule.sequence = Sequence::new(alphabet, text).unwrap();
        let report = search_document(
            &molecule,
            &SearchOptions {
                query: query.into(),
                mode: SearchMode::Sequence,
                both_strands: false,
                ..Default::default()
            },
        )
        .unwrap();
        assert!(report.complete);
        assert_eq!(
            report
                .hits
                .iter()
                .map(|h| h.parts[0].start)
                .collect::<Vec<_>>(),
            expected
        );
        assert!(report.hits.iter().all(|h| h.certain));
        report.validate().unwrap();
        if alphabet != Alphabet::Protein {
            assert_eq!(
                analysis::find_matches(&molecule, query, false)
                    .unwrap()
                    .iter()
                    .map(|h| h.start)
                    .collect::<Vec<_>>(),
                expected
            );
        }
    }
}

#[test]
fn shared_compatible_kernel_preserves_subset_certainty_and_reports_work_exhaustion() {
    use collagene_core::sequence::analysis;
    let molecule = Molecule::dna("compatible", "ANGRA", Topology::Circular).unwrap();
    let options = SearchOptions {
        query: "RA".into(),
        mode: SearchMode::Sequence,
        both_strands: false,
        ..Default::default()
    };
    let report = search_document(&molecule, &options).unwrap();
    // R versus R is certain; N and R versus the terminal A remain uncertain.
    // The last A crosses the stored origin.
    assert_eq!(
        report
            .hits
            .iter()
            .map(|h| (h.parts[0].start, h.certain))
            .collect::<Vec<_>>(),
        vec![(0, false), (2, false), (3, true), (4, true)]
    );
    assert_eq!(
        analysis::find_matches(&molecule, "RA", false)
            .unwrap()
            .iter()
            .map(|h| (h.start, h.certain))
            .collect::<Vec<_>>(),
        vec![(0, false), (2, false), (3, true), (4, true)]
    );
    let partial = search_document(
        &molecule,
        &SearchOptions {
            max_work: 3,
            ..options
        },
    )
    .unwrap();
    assert!(!partial.complete);
    assert!(partial.work <= 3);
    partial.validate().unwrap();
}

#[test]
fn search_report_rejects_a_superseded_engine_identity() {
    let molecule = Molecule::dna("identity", "ACGT", Topology::Linear).unwrap();
    let mut report = search_document(
        &molecule,
        &SearchOptions {
            query: "AC".into(),
            mode: SearchMode::Sequence,
            ..Default::default()
        },
    )
    .unwrap();
    assert_eq!(report.engine_version, SEARCH_ENGINE_VERSION);
    report.engine_version = "science-35".into();
    assert!(report.validate().is_err());
}
