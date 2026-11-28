use collagene_core::search::*;
use collagene_model::{Alphabet, Feature, Location, Molecule, Sequence, Topology};
fn options(query: &str, mode: SearchMode) -> SearchOptions {
    SearchOptions {
        query: query.into(),
        mode,
        ..Default::default()
    }
}
fn dna(sequence: &str, circular: bool) -> Molecule {
    Molecule::dna(
        "source",
        sequence,
        if circular {
            Topology::Circular
        } else {
            Topology::Linear
        },
    )
    .unwrap()
}
#[test]
fn circular_reverse_spans_follow_biological_order_and_ids_bind_options() {
    let m = dna("TCAA", true);
    let mut o = options("ATT", SearchMode::Sequence);
    let r = search_document(&m, &o).unwrap();
    assert!(r.complete);
    assert_eq!(r.hits.len(), 1);
    assert_eq!(
        r.hits[0].parts,
        vec![
            SearchSpan {
                start: 0,
                end: 1,
                reverse: true
            },
            SearchSpan {
                start: 2,
                end: 4,
                reverse: true
            }
        ]
    );
    assert_eq!(r, search_document(&m, &o).unwrap());
    o.max_hits = 2;
    assert_ne!(r.hits[0].id, search_document(&m, &o).unwrap().hits[0].id);
}
#[test]
fn exhaustive_small_dna_matches_an_independent_literal_oracle() {
    let alphabet = b"ACGT";
    for number in 0..256 {
        let mut number = number;
        let mut text = String::new();
        for _ in 0..4 {
            text.push(alphabet[number % 4] as char);
            number /= 4;
        }
        for circular in [false, true] {
            let m = dna(&text, circular);
            for q in ["AA", "CG", "TT", "ATG"] {
                let mut o = options(q, SearchMode::Sequence);
                o.both_strands = false;
                let r = search_document(&m, &o).unwrap();
                let expected: Vec<usize> = (0..if circular { 4 } else { 5 - q.len() })
                    .filter(|&start| {
                        (0..q.len()).all(|j| text.as_bytes()[(start + j) % 4] == q.as_bytes()[j])
                    })
                    .collect();
                assert_eq!(
                    r.hits.iter().map(|h| h.parts[0].start).collect::<Vec<_>>(),
                    expected,
                    "{text} {q} {circular}"
                );
            }
        }
    }
}
#[test]
fn ambiguous_certainty_rna_and_protein_are_typed() {
    let m = dna("ANAC", false);
    let mut o = options("AA", SearchMode::Sequence);
    o.both_strands = false;
    let r = search_document(&m, &o).unwrap();
    assert_eq!(r.hits.len(), 2);
    assert!(r.hits.iter().all(|h| !h.certain));
    let mut rna = dna("AAAA", false);
    rna.sequence = Sequence::new(Alphabet::Rna, "AUGU").unwrap();
    assert_eq!(
        search_document(&rna, &options("ACA", SearchMode::Sequence))
            .unwrap()
            .hits[0]
            .parts[0],
        SearchSpan {
            start: 1,
            end: 4,
            reverse: true
        }
    );
    assert!(search_document(&rna, &options("AT", SearchMode::Sequence)).is_err());
    let mut protein = dna("AAAA", false);
    protein.sequence = Sequence::new(Alphabet::Protein, "MANA").unwrap();
    let r = search_document(&protein, &options("AN", SearchMode::Sequence)).unwrap();
    assert_eq!(r.hits.len(), 1);
    assert_eq!(r.hits[0].parts[0].start, 1);
    assert!(!r.hits[0].parts[0].reverse);
}
#[test]
fn annotation_and_primer_modes_keep_compound_order_and_text_policy() {
    let mut m = dna("AAAAAAAAAAAA", false);
    m.features = vec![Feature {
        kind: "primer_bind".into(),
        location: Location::Complement(Box::new(Location::Join(vec![
            Location::range(1, 3),
            Location::range(6, 9),
        ]))),
        qualifiers: vec![
            ("label".into(), Some("Gene Alpha".into())),
            ("note".into(), Some("alpha_beta".into())),
        ],
        needs_review: false,
    }];
    let mut o = options("alpha", SearchMode::Primers);
    o.whole_word = true;
    let r = search_document(&m, &o).unwrap();
    assert_eq!(r.hits[0].kind, SearchKind::Primer);
    assert_eq!(r.hits[0].feature_index, Some(0));
    assert_eq!(
        r.hits[0].parts.iter().map(|p| p.start).collect::<Vec<_>>(),
        vec![6, 1]
    );
    o.case_sensitive = true;
    assert!(search_document(&m, &o).unwrap().hits.is_empty());
    o.whole_word = false;
    assert_eq!(search_document(&m, &o).unwrap().hits.len(), 1);
}
#[test]
fn enzymes_retain_real_type_iis_cuts_and_linear_missing_context() {
    let m = dna("GGTCTCAAAAAAAA", false);
    let r = search_document(&m, &options("BsaI", SearchMode::Enzymes)).unwrap();
    assert_eq!(r.hits.len(), 1);
    assert_eq!(r.hits[0].top_cut, Some(7));
    assert_eq!(r.hits[0].bottom_cut, Some(11));
    assert_eq!(r.hits[0].enzyme.as_deref(), Some("BsaI"));
    let r = search_document(&dna("GGTCTC", false), &options("BsaI", SearchMode::Enzymes)).unwrap();
    assert_eq!(r.hits[0].top_cut, None);
    assert_eq!(r.hits[0].bottom_cut, None);
}
#[test]
fn each_bound_is_explicit_and_all_mode_keeps_non_sequence_sources() {
    let mut m = dna("AAAAAAAAAAAA", false);
    m.features = (0..3)
        .map(|_| Feature {
            kind: "gene".into(),
            location: Location::range(0, 2),
            qualifiers: vec![],
            needs_review: false,
        })
        .collect();
    let mut o = options("A", SearchMode::Sequence);
    o.max_hits = 2;
    let r = search_document(&m, &o).unwrap();
    assert_eq!(r.hits.len(), 2);
    assert!(!r.complete);
    o.max_hits = 100;
    o.max_work = 2;
    let r = search_document(&m, &o).unwrap();
    assert_eq!(r.work, 2);
    assert!(!r.complete);
    o = options("gene", SearchMode::All);
    o.max_rows = 2;
    let r = search_document(&m, &o).unwrap();
    assert!(!r.complete);
    assert_eq!(r.rows, 2);
    assert_eq!(r.hits.len(), 2);
    assert!(!r.warnings.is_empty());
    o.max_rows = 100;
    assert!(search_document(&m, &o).unwrap().complete);
}
#[test]
fn every_iupac_symbol_pair_matches_a_concrete_set_oracle() {
    let symbols = [
        ('A', "A"),
        ('C', "C"),
        ('G', "G"),
        ('T', "T"),
        ('R', "AG"),
        ('Y', "CT"),
        ('S', "CG"),
        ('W', "AT"),
        ('K', "GT"),
        ('M', "AC"),
        ('B', "CGT"),
        ('D', "AGT"),
        ('H', "ACT"),
        ('V', "ACG"),
        ('N', "ACGT"),
    ];
    for (a, aset) in symbols {
        for (b, bset) in symbols {
            let m = dna(&a.to_string(), false);
            let mut o = options(&b.to_string(), SearchMode::Sequence);
            o.both_strands = false;
            let r = search_document(&m, &o).unwrap();
            let possible = aset.chars().any(|c| bset.contains(c));
            assert_eq!(!r.hits.is_empty(), possible, "{a} {b}");
            if possible {
                assert_eq!(
                    r.hits[0].certain,
                    aset.chars().all(|c| bset.contains(c)),
                    "{a} {b}"
                );
            }
        }
    }
}
#[test]
fn unresolved_and_fuzzy_annotation_hits_never_invent_exact_evidence() {
    let mut m = dna("AAAA", false);
    m.features = vec![
        Feature {
            kind: "gene".into(),
            location: Location::Unresolved("remote".into()),
            qualifiers: vec![],
            needs_review: false,
        },
        Feature {
            kind: "gene".into(),
            location: Location::Complement(Box::new(Location::Range {
                start: 1,
                end: 3,
                before: true,
                after: false,
            })),
            qualifiers: vec![],
            needs_review: false,
        },
    ];
    let r = search_document(&m, &options("gene", SearchMode::Annotations)).unwrap();
    assert!(r.hits[0].parts.is_empty());
    assert!(r.hits.iter().all(|h| !h.certain));
    r.validate().unwrap();
}
#[test]
fn huge_compound_and_skipped_primer_rows_respect_limits_without_partial_coordinates() {
    let mut m = dna("AAAA", false);
    m.features = vec![Feature {
        kind: "gene".into(),
        location: Location::Join(vec![Location::range(0, 1); 4097]),
        qualifiers: vec![],
        needs_review: false,
    }];
    let r = search_document(&m, &options("gene", SearchMode::Annotations)).unwrap();
    assert!(!r.complete);
    assert!(r.hits.is_empty());
    assert!(r.warnings.iter().any(|w| w.contains("4096")));
    r.validate().unwrap();
    m.features = vec![
        Feature {
            kind: "gene".into(),
            location: Location::range(0, 1),
            qualifiers: vec![],
            needs_review: false,
        },
        Feature {
            kind: "primer_bind".into(),
            location: Location::range(0, 1),
            qualifiers: vec![],
            needs_review: false,
        },
    ];
    let mut o = options("primer", SearchMode::Primers);
    o.max_rows = 1;
    let r = search_document(&m, &o).unwrap();
    assert!(!r.complete);
    assert!(r.hits.is_empty());
    assert_eq!(r.rows, 1);
}
