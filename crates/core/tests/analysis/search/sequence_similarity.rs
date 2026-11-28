use collagene_core::sequence_similarity::{self as search, Limit, Method, Options, Source, Span};
use collagene_model::{Alphabet, Molecule, Sequence, Topology};
fn source(key: &str, sequence: &str, topology: Topology) -> Source {
    Source {
        key: key.into(),
        molecule: Molecule::dna(key, sequence, topology).unwrap(),
    }
}
fn run(query: &Source, subjects: &[Source], options: &Options) -> search::Report {
    search::search(query, subjects, options, &mut || false).unwrap()
}
fn reverse(s: &str) -> String {
    s.bytes()
        .rev()
        .map(|b| match b {
            b'A' => 'T',
            b'T' => 'A',
            b'C' => 'G',
            b'G' => 'C',
            b'N' => 'N',
            _ => panic!("oracle alphabet"),
        })
        .collect()
}
// Independent tiny affine local-score oracle: no call to the application's aligner.
fn local_score(reference: &str, query: &str) -> i32 {
    let n = reference.len();
    let m = query.len();
    let mut h = vec![vec![0; m + 1]; n + 1];
    let mut deletion = vec![vec![-100000; m + 1]; n + 1];
    let mut insertion = deletion.clone();
    for i in 1..=n {
        for j in 1..=m {
            deletion[i][j] = (h[i - 1][j] - 4).max(deletion[i - 1][j] - 1);
            insertion[i][j] = (h[i][j - 1] - 4).max(insertion[i][j - 1] - 1);
            let a = reference.as_bytes()[i - 1];
            let b = query.as_bytes()[j - 1];
            let score = if a == b && a != b'N' {
                2
            } else if a == b'N' || b == b'N' {
                0
            } else {
                -3
            };
            h[i][j] = 0
                .max(h[i - 1][j - 1] + score)
                .max(deletion[i][j])
                .max(insertion[i][j]);
        }
    }
    h.into_iter().flatten().max().unwrap()
}
fn oracle(reference: &str, query: &str, circular: bool, both: bool) -> i32 {
    let mut best = 0;
    for origin in 0..if circular { reference.len() } else { 1 } {
        let rotated = format!("{}{}", &reference[origin..], &reference[..origin]);
        best = best.max(local_score(&rotated, query));
        if both {
            best = best.max(local_score(&rotated, &reverse(query)));
        }
    }
    best
}
#[test]
fn canonical_exact_matches_prove_the_score_bound_without_matrix_work() {
    let query = source("query", "CGAA", Topology::Linear);
    let reference = source("ring", "AATTCCG", Topology::Circular);
    let options = Options {
        max_cells_per_alignment: 1,
        max_origins: 1,
        ..Default::default()
    };
    let report = run(&query, &[reference], &options);
    assert!(report.complete);
    let row = &report.subjects[0];
    assert!(row.score_upper_bound_reached);
    assert_eq!(row.alignments_evaluated, 0);
    let hit = &report.hits[0];
    assert_eq!(hit.method, Method::ExactCanonicalSubstring);
    assert_eq!(hit.score, 8);
    assert_eq!(
        hit.reference_spans,
        vec![Span { start: 5, end: 7 }, Span { start: 0, end: 2 }]
    );
    assert_eq!(hit.query_span, Span { start: 0, end: 4 });
    assert!(!hit.query_reversed);
    assert_eq!(hit.statistics.certain_matches, 4);
    assert_eq!(hit.statistics.query_coverage, 1.0);
    assert_eq!(report.work_used, 14);
    let mut custom = options;
    custom.scoring.match_score = 7;
    let scored = run(
        &query,
        &[source("ring", "AATTCCG", Topology::Circular)],
        &custom,
    );
    assert_eq!(scored.hits[0].score, 28);
    assert!(scored.complete);
    assert_ne!(scored.context_revision, report.context_revision);
}
#[test]
fn reverse_complement_and_partial_reverse_coordinates_refer_to_original_query() {
    let query = source("query", "AGTC", Topology::Linear);
    let report = run(
        &query,
        &[source("reference", "GGGACTTT", Topology::Linear)],
        &Options::default(),
    );
    let hit = &report.hits[0];
    assert!(hit.query_reversed);
    assert_eq!(hit.reference_spans, vec![Span { start: 2, end: 6 }]);
    assert_eq!(hit.query_span, Span { start: 0, end: 4 });
    let query = source("query", "AAGTC", Topology::Linear);
    let report = run(
        &query,
        &[source("reference", "GAC", Topology::Linear)],
        &Options::default(),
    );
    let hit = &report.hits[0];
    assert!(hit.query_reversed);
    assert_eq!(hit.query_span, Span { start: 2, end: 5 });
    assert_eq!(hit.score, 6);
    assert_eq!(hit.statistics.paired_bases, 3);
    assert!(report.complete);
}
#[test]
fn affine_gaps_and_ambiguity_keep_exact_statistics_and_uncertainty() {
    let options = Options {
        both_query_strands: false,
        ..Default::default()
    };
    let report = run(
        &source("query", "ACGTTGAC", Topology::Linear),
        &[source("subject", "ACGTGAC", Topology::Linear)],
        &options,
    );
    let hit = &report.hits[0];
    assert_eq!(hit.score, 10);
    assert_eq!(hit.statistics.inserted_query_bases, 1);
    assert_eq!(hit.statistics.certain_matches, 7);
    assert_eq!(hit.query_span, Span { start: 0, end: 8 });
    assert_eq!(hit.reference_spans, vec![Span { start: 0, end: 7 }]);
    assert!(report.complete);
    let report = run(
        &source("query", "ANA", Topology::Linear),
        &[source("subject", "AAA", Topology::Linear)],
        &options,
    );
    let hit = &report.hits[0];
    assert_eq!(hit.method, Method::AffineLocal);
    assert_eq!(hit.score, 4);
    assert_eq!(hit.statistics.certain_matches, 2);
    assert_eq!(hit.statistics.ambiguous_compatible, 1);
    assert_eq!(hit.statistics.certain_identity, Some(2.0 / 3.0));
    assert_eq!(hit.statistics.compatible_identity, Some(1.0));
}
#[test]
fn tiny_linear_and_every_circular_rotation_agree_with_independent_score_oracle() {
    let options = Options {
        max_origins: 64,
        ..Default::default()
    };
    for query in ["ACG", "AGTC", "AAGTC", "ANA", "ACGTTGAC"] {
        for subject in ["ACG", "GGGACTTT", "GAC", "AAA", "ACGTGAC", "TGCAT"] {
            for circular in [false, true] {
                let q = source("q", query, Topology::Linear);
                let s = source(
                    "s",
                    subject,
                    if circular {
                        Topology::Circular
                    } else {
                        Topology::Linear
                    },
                );
                let report = run(&q, &[s], &options);
                let expected = oracle(subject, query, circular, true);
                assert!(report.complete, "{query} {subject} {circular}");
                assert_eq!(
                    report.subjects[0]
                        .best
                        .as_ref()
                        .map(|h| h.score)
                        .unwrap_or(0),
                    expected,
                    "{query} {subject} {circular}"
                );
            }
        }
    }
}
#[test]
fn incomplete_origin_work_cell_and_cancel_limits_never_become_absence_proof() {
    let query = source("q", "NACG", Topology::Linear);
    let subject = source("s", "TGCATCG", Topology::Circular);
    let partial = run(
        &query,
        &[subject.clone()],
        &Options {
            max_origins: 1,
            ..Default::default()
        },
    );
    assert!(!partial.complete);
    assert!(partial.subjects[0].limits.contains(&Limit::CircularOrigins));
    let full = run(
        &query,
        &[subject],
        &Options {
            max_origins: 64,
            ..Default::default()
        },
    );
    assert!(full.complete);
    let report = run(
        &query,
        &[source("s", "TTTT", Topology::Linear)],
        &Options {
            max_cells_per_alignment: 1,
            ..Default::default()
        },
    );
    assert!(!report.complete);
    assert!(report.subjects[0].limits.contains(&Limit::AlignmentCells));
    let report = run(
        &source("q", "ACG", Topology::Linear),
        &[source("s", "TTTT", Topology::Linear)],
        &Options {
            max_work: 1,
            ..Default::default()
        },
    );
    assert!(!report.complete);
    assert_eq!(report.work_used, 0);
    assert!(report.subjects[0].limits.contains(&Limit::Work));
    // This multiplication overflows usize on WASM32. Both targets must return
    // a per-subject cell limit, rather than a platform-dependent whole-job error.
    let large = run(
        &source("q", &"A".repeat(65536), Topology::Linear),
        &[source("s", &"C".repeat(65536), Topology::Linear)],
        &Options::default(),
    );
    assert!(!large.complete);
    assert!(large.subjects[0].limits.contains(&Limit::AlignmentCells));
    let subjects = vec![
        source("s1", "ACG", Topology::Linear),
        source("s2", "ACG", Topology::Linear),
    ];
    let report = search::search(&query, &subjects, &Options::default(), &mut || true).unwrap();
    assert!(!report.complete);
    assert_eq!(report.work_used, 0);
    assert!(
        report
            .subjects
            .iter()
            .all(|r| r.limits == vec![Limit::Cancelled])
    );
}
#[test]
fn ranking_is_stable_truncation_is_explicit_and_self_search_is_allowed() {
    let query = source("q", "ACG", Topology::Linear);
    let subjects = vec![
        source("z", "ACG", Topology::Linear),
        source("a", "ACG", Topology::Linear),
        query.clone(),
    ];
    let report = run(
        &query,
        &subjects,
        &Options {
            max_hits: 1,
            ..Default::default()
        },
    );
    assert!(report.complete);
    assert!(report.hits_truncated);
    assert_eq!(report.hits[0].subject_key, "a");
    assert_eq!(report.subjects.len(), 3);
    let mut renamed = query.clone();
    renamed.key = "other-query".into();
    assert_ne!(
        search::context_revision(&query, &subjects, &Options::default()).unwrap(),
        search::context_revision(&renamed, &subjects, &Options::default()).unwrap()
    );
    assert!(
        search::preflight(
            &query,
            &[subjects[0].clone(), subjects[0].clone()],
            &Options::default()
        )
        .is_err()
    );
    let report = run(
        &query,
        &subjects,
        &Options {
            minimum_score: 7,
            ..Default::default()
        },
    );
    assert!(report.complete);
    assert!(report.hits.is_empty());
    assert!(report.subjects.iter().all(|s| s.best.is_some()));
}
#[test]
fn a_circular_subject_is_never_repeated_twice_to_fake_a_full_query_match() {
    let report = run(
        &source("q", "AAAA", Topology::Linear),
        &[source("s", "AA", Topology::Circular)],
        &Options {
            max_origins: 64,
            ..Default::default()
        },
    );
    assert!(report.complete);
    assert_eq!(report.hits[0].score, 4);
    assert!(!report.subjects[0].score_upper_bound_reached);
    assert_eq!(report.hits[0].statistics.paired_bases, 2);
}
#[test]
fn rna_has_its_own_alphabet_and_incompatible_subjects_remain_visible() {
    let mut q = source("rna-query", "ATGC", Topology::Linear);
    q.molecule.sequence = Sequence::new(Alphabet::Rna, "AUGC").unwrap();
    let mut r = source("rna-reference", "CCGCAT", Topology::Linear);
    r.molecule.sequence = Sequence::new(Alphabet::Rna, "CCGCAU").unwrap();
    let report = run(
        &q,
        &[r, source("dna", "ATGC", Topology::Linear)],
        &Options::default(),
    );
    assert!(!report.complete);
    assert_eq!(report.hits[0].score, 8);
    assert!(report.hits[0].query_reversed);
    assert_eq!(report.subjects[1].limits, vec![Limit::IncompatibleAlphabet]);
    let mut protein = q;
    protein.molecule.sequence = Sequence::new(Alphabet::Protein, "MPEP").unwrap();
    assert!(
        search::search(
            &protein,
            &[source("dna", "ATGC", Topology::Linear)],
            &Options::default(),
            &mut || false
        )
        .is_err()
    );
}
