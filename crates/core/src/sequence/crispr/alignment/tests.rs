use super::*;
use std::collections::BTreeSet;

fn limits() -> AlignmentLimits {
    AlignmentLimits {
        max_mismatches: 0,
        max_dna_bulge_bases: 1,
        max_rna_bulge_bases: 1,
        max_total_bulge_bases: 1,
        max_alignments: 4096,
    }
}
fn run(q: &str, t: &str, l: AlignmentLimits) -> AlignmentEvidence {
    align(
        &Sequence::dna(q).unwrap(),
        &Sequence::dna(t).unwrap(),
        l,
        &mut AlignmentWork::new(1_000_000, 1_000_000).unwrap(),
    )
    .unwrap()
}
fn word(a: &GuideAlignment) -> String {
    a.columns
        .iter()
        .map(|c| match c {
            AlignmentColumn::Paired { .. } => 'P',
            AlignmentColumn::DnaBulge { .. } => 'D',
            AlignmentColumn::RnaBulge { .. } => 'R',
        })
        .collect()
}

#[test]
fn all_repeated_base_gap_placements_are_retained_in_both_directions() {
    let dna = run("AAA", "AAAA", limits());
    assert!(dna.complete);
    assert_eq!(
        dna.alignments.iter().map(word).collect::<Vec<_>>(),
        ["PPPD", "PPDP", "PDPP", "DPPP"]
    );
    assert!(
        dna.alignments
            .iter()
            .all(|a| a.dna_bulge_bases == 1 && a.rna_bulge_bases == 0 && a.definite)
    );
    let rna = run("AAAA", "AAA", limits());
    assert!(rna.complete);
    assert_eq!(
        rna.alignments.iter().map(word).collect::<Vec<_>>(),
        ["PPPR", "PPRP", "PRPP", "RPPP"]
    );
}

#[test]
fn ambiguity_including_identical_n_calls_never_becomes_false_certainty() {
    let l = AlignmentLimits {
        max_dna_bulge_bases: 0,
        max_rna_bulge_bases: 0,
        max_total_bulge_bases: 0,
        ..limits()
    };
    let possible = run(
        "ANR",
        "ANY",
        AlignmentLimits {
            max_mismatches: 1,
            ..l
        },
    );
    let a = &possible.alignments[0];
    assert_eq!(
        (a.min_mismatches, a.max_mismatches, a.definite),
        (1, 2, false)
    );
    assert!(run("ANR", "ANY", l).alignments.is_empty());
    let rna = Sequence::new(Alphabet::Rna, "AUG").unwrap();
    let dna = Sequence::dna("ATG").unwrap();
    let matched = align(
        &rna,
        &dna,
        l,
        &mut AlignmentWork::new(100, 1_000_000).unwrap(),
    )
    .unwrap();
    assert_eq!(
        (
            matched.alignments[0].min_mismatches,
            matched.alignments[0].max_mismatches
        ),
        (0, 0)
    );
}

#[test]
fn every_iupac_pair_matches_an_independent_concrete_realization_oracle() {
    let codes = [
        ("A", "A"),
        ("C", "C"),
        ("G", "G"),
        ("T", "T"),
        ("R", "AG"),
        ("Y", "CT"),
        ("S", "CG"),
        ("W", "AT"),
        ("K", "GT"),
        ("M", "AC"),
        ("B", "CGT"),
        ("D", "AGT"),
        ("H", "ACT"),
        ("V", "ACG"),
        ("N", "ACGT"),
    ];
    let l = AlignmentLimits {
        max_mismatches: 1,
        max_dna_bulge_bases: 0,
        max_rna_bulge_bases: 0,
        max_total_bulge_bases: 0,
        ..limits()
    };
    for (q, concrete_q) in codes {
        for (t, concrete_t) in codes {
            let mut counts = Vec::new();
            for qb in concrete_q.bytes() {
                for tb in concrete_t.bytes() {
                    counts.push(usize::from(qb != tb));
                }
            }
            let e = run(q, t, l);
            assert!(e.complete);
            assert_eq!(e.alignments.len(), 1);
            let a = &e.alignments[0];
            assert_eq!(
                (a.min_mismatches, a.max_mismatches),
                (*counts.iter().min().unwrap(), *counts.iter().max().unwrap()),
                "{q}/{t}"
            );
        }
    }
}

#[test]
fn exact_output_cap_can_be_complete_and_truncation_is_explicit() {
    let full = run(
        "AAA",
        "AAAA",
        AlignmentLimits {
            max_alignments: 4,
            ..limits()
        },
    );
    assert!(full.complete);
    let capped = run(
        "AAA",
        "AAAA",
        AlignmentLimits {
            max_alignments: 3,
            ..limits()
        },
    );
    assert!(!capped.complete);
    assert_eq!(
        capped.limit_reached,
        Some(AlignmentLimitReached::Alignments)
    );
    assert_eq!(capped.alignments, full.alignments[..3]);
}

#[test]
fn shared_work_exhaustion_does_not_return_false_absence() {
    let q = Sequence::dna("AAA").unwrap();
    let t = Sequence::dna("AAA").unwrap();
    let mut work = AlignmentWork::new(13, 1_000_000).unwrap();
    let first = align(&q, &t, limits(), &mut work).unwrap();
    assert!(!first.complete); // one definite hit is retained, remaining paths unexamined
    assert_eq!(first.alignments.len(), 1);
    assert_eq!(work.used(), 13);
    let second = align(&q, &t, limits(), &mut work).unwrap();
    assert!(!second.complete);
    assert!(second.alignments.is_empty());
    assert_eq!(second.limit_reached, Some(AlignmentLimitReached::Work));
    assert_eq!(second.work, 0);
}

/// Independent oracle: enumerate fixed-length words over P/D/R, not a pruned
/// alignment recursion. Reconstruct consumption and literal mismatches afresh.
fn oracle(q: &[u8], t: &[u8], l: AlignmentLimits) -> BTreeSet<String> {
    let mut found = BTreeSet::new();
    for width in q.len().max(t.len())..=q.len() + t.len() {
        for mut encoding in 0..3usize.pow(width as u32) {
            let (mut qi, mut ti, mut dna, mut rna, mut mismatches) = (0, 0, 0, 0, 0);
            let mut operations = String::new();
            let mut valid = true;
            for _ in 0..width {
                let operation = encoding % 3;
                encoding /= 3;
                match operation {
                    0 if qi < q.len() && ti < t.len() => {
                        mismatches += usize::from(q[qi] != t[ti]);
                        qi += 1;
                        ti += 1;
                        operations.push('P');
                    }
                    1 if ti < t.len() => {
                        ti += 1;
                        dna += 1;
                        operations.push('D');
                    }
                    2 if qi < q.len() => {
                        qi += 1;
                        rna += 1;
                        operations.push('R');
                    }
                    _ => {
                        valid = false;
                        break;
                    }
                }
            }
            if valid
                && qi == q.len()
                && ti == t.len()
                && dna <= l.max_dna_bulge_bases
                && rna <= l.max_rna_bulge_bases
                && dna + rna <= l.max_total_bulge_bases
                && mismatches <= l.max_mismatches
            {
                found.insert(operations);
            }
        }
    }
    found
}

#[test]
fn bounded_kernel_matches_independent_exhaustive_literal_oracle() {
    let sequences = ["A", "C", "AA", "AC", "CA", "CC", "AAA", "ACA", "CAC", "CCC"];
    for q in sequences {
        for t in sequences {
            for total in 0..=2 {
                let l = AlignmentLimits {
                    max_mismatches: 1,
                    max_dna_bulge_bases: 2,
                    max_rna_bulge_bases: 2,
                    max_total_bulge_bases: total,
                    ..limits()
                };
                let evidence = run(q, t, l);
                assert!(evidence.complete);
                let expected = oracle(q.as_bytes(), t.as_bytes(), l);
                let actual = evidence
                    .alignments
                    .iter()
                    .map(word)
                    .collect::<BTreeSet<_>>();
                assert_eq!(actual, expected, "guide={q}, target={t}, total={total}");
                assert_eq!(actual.len(), evidence.alignments.len());
                assert!(evidence.alignments.iter().all(|a| a.definite));
            }
        }
    }
}

#[test]
fn retained_column_budget_is_shared_and_checked_before_path_cloning() {
    let q = Sequence::dna("AAA").unwrap();
    let t = Sequence::dna("AAAA").unwrap();
    let mut work = AlignmentWork::new(10000, 4).unwrap();
    let first = align(&q, &t, limits(), &mut work).unwrap();
    assert_eq!(first.alignments.len(), 1);
    assert!(!first.complete);
    assert_eq!(
        first.limit_reached,
        Some(AlignmentLimitReached::RetainedColumns)
    );
    let second = align(&q, &t, limits(), &mut work).unwrap();
    assert!(second.alignments.is_empty());
    assert!(!second.complete);
    assert_eq!(
        second.limit_reached,
        Some(AlignmentLimitReached::RetainedColumns)
    );
}
