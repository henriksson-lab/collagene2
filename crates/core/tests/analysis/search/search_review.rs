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
