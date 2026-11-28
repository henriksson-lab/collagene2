//! Independent source conservation, vote and exhaustive affine-boundary oracles.
use collagene_core::{contig::*, oligo_formats};
use collagene_model::{Molecule, Topology};
use std::collections::{BTreeMap, BTreeSet};
const REFERENCE: &str = "AACCGTAGTCTGACGATGCCTAGGTCATCA";
fn plain(id: &str, s: &str) -> ReadInput {
    ReadInput {
        id: id.into(),
        molecule: Molecule::dna(id, s, Topology::Linear).unwrap(),
    }
}
fn options() -> ContigOptions {
    let mut o = ContigOptions::default();
    o.min_overlap = 6;
    o.trim.min_retained_bases = 1;
    o.both_strands = false;
    o.max_error_bp = 2000;
    o
}
fn fastq(id: &str, s: &str, low: Option<usize>) -> ReadInput {
    let mut q = vec![b'I'; s.len()];
    if let Some(i) = low {
        q[i] = b'!';
    }
    let text = format!("@{id}\n{s}\n+\n{}\n", String::from_utf8(q).unwrap());
    ReadInput {
        id: id.into(),
        molecule: oligo_formats::import_fastq(
            text.as_bytes(),
            oligo_formats::QualityEncoding::Phred33,
        )
        .unwrap()
        .molecules
        .remove(0),
    }
}
fn complement(b: u8) -> u8 {
    match b {
        b'A' => b'T',
        b'C' => b'G',
        b'G' => b'C',
        b'T' => b'A',
        b'N' => b'N',
        _ => panic!("independent fixture is concrete"),
    }
}
fn reverse(s: &str) -> String {
    s.bytes().rev().map(|b| complement(b) as char).collect()
}
fn source_audit(r: &ContigReport, inputs: &[ReadInput], o: &ContigOptions) {
    for c in &r.contigs {
        let mut calls = BTreeSet::new();
        for p in &c.placements {
            let source = &inputs[p.read_index].molecule;
            let rr = &r.reads[p.read_index];
            assert_eq!(p.calls.len(), rr.trim_end - rr.trim_start);
            for call in &p.calls {
                assert!(call.call_index >= rr.trim_start && call.call_index < rr.trim_end);
                assert!(call.column < c.columns.len());
                assert!(calls.insert((p.read_index, call.call_index)));
                let expected = source.sequence.as_str().as_bytes()[call.call_index];
                let expected = if p.reverse {
                    complement(expected)
                } else {
                    expected
                };
                let obs = c.columns[call.column]
                    .observations
                    .iter()
                    .filter(|v| v.read_index == p.read_index && v.call_index == call.call_index)
                    .collect::<Vec<_>>();
                assert_eq!(obs.len(), 1);
                assert_eq!(obs[0].base, expected as char);
                let stored = if rr.source_reverse {
                    rr.retained_end - 1 - call.call_index
                } else {
                    rr.retained_start + call.call_index
                };
                assert_eq!(call.stored_call_index, stored);
                assert_eq!(obs[0].stored_call_index, stored);
                if let Ok(q) = oligo_formats::quality_evidence(source) {
                    assert_eq!(obs[0].raw_score, Some(q.qualities[call.call_index]));
                    assert_eq!(
                        obs[0].eligible,
                        q.qualities[call.call_index] >= o.trim.raw_score_threshold
                    );
                } else {
                    assert!(obs[0].raw_score.is_none());
                    assert_eq!(obs[0].eligible, o.trim.allow_unscored);
                }
            }
        }
        assert_eq!(
            calls.len(),
            c.columns
                .iter()
                .map(|x| x.observations.len())
                .sum::<usize>()
        );
        let mut sequence = String::new();
        for (col_index, col) in c.columns.iter().enumerate() {
            assert_eq!(col.layout_index, col_index);
            let mut counted = BTreeMap::<String, u8>::new();
            for obs in &col.observations {
                if obs.eligible && b"ACGT".contains(&(obs.base as u8)) {
                    let key = r.reads[obs.read_index]
                        .duplicate_group
                        .clone()
                        .unwrap_or_else(|| format!("read:{}", obs.read_index));
                    if let Some(previous) = counted.insert(key, obs.base as u8) {
                        assert_eq!(
                            previous, obs.base as u8,
                            "duplicate evidence disagrees within a column"
                        );
                    }
                }
            }
            let mut counts = [0; 4];
            for base in counted.values() {
                counts[b"ACGT".iter().position(|b| b == base).unwrap()] += 1;
            }
            assert_eq!(col.counts, counts);
            for read in &col.eligible_gap_reads {
                assert!(col.gap_reads.contains(read));
                let before = c.columns[..col_index]
                    .iter()
                    .rev()
                    .find_map(|v| v.observations.iter().find(|v| v.read_index == *read));
                let after = c.columns[col_index + 1..]
                    .iter()
                    .find_map(|v| v.observations.iter().find(|v| v.read_index == *read));
                assert!(before.is_some_and(|v| v.eligible) && after.is_some_and(|v| v.eligible));
            }
            if let Some(base) = col.base {
                assert_eq!(col.consensus_position, Some(sequence.len()));
                sequence.push(base);
                if base != 'N' {
                    assert!(
                        col.observations
                            .iter()
                            .any(|v| v.eligible && v.base == base),
                        "consensus invented a base without eligible source evidence"
                    );
                    assert!(
                        col.eligible_gap_reads.is_empty(),
                        "concrete base suppressed contradictory eligible gaps"
                    );
                }
            } else {
                assert!(col.consensus_position.is_none());
            }
        }
        assert_eq!(sequence, c.sequence);
        if let Some(m) = &c.molecule {
            assert_eq!(m.sequence.as_str(), sequence);
        }
    }
}
#[test]
fn three_tiled_reads_reconstruct_the_literal_reference_without_reference_input() {
    let inputs = vec![
        plain("a", &REFERENCE[..18]),
        plain("b", &REFERENCE[9..27]),
        plain("c", &REFERENCE[18..]),
    ];
    let o = options();
    let r = assemble(&inputs, &o).unwrap();
    source_audit(&r, &inputs, &o);
    assert!(r.complete);
    assert!(
        r.contigs
            .iter()
            .any(|c| c.sequence == REFERENCE && c.placements.len() == 3)
    );
}
#[test]
fn a_high_quality_substitution_conflicts_but_a_low_quality_call_cannot_outvote_supported_bases() {
    let mut wrong = REFERENCE[9..].as_bytes().to_vec();
    wrong[5] = if wrong[5] == b'A' { b'C' } else { b'A' };
    let wrong = String::from_utf8(wrong).unwrap();
    let o = options();
    let inputs = vec![fastq("a", &REFERENCE[..20], None), fastq("b", &wrong, None)];
    let r = assemble(&inputs, &o).unwrap();
    source_audit(&r, &inputs, &o);
    let mut expected = REFERENCE.as_bytes().to_vec();
    expected[14] = b'N';
    let expected = String::from_utf8(expected).unwrap();
    assert!(r.contigs.iter().any(|c| {
        c.sequence == expected
            && c.columns
                .iter()
                .any(|x| x.state == ConsensusState::ConflictingCalls)
    }));
    let inputs = vec![
        fastq("a", &REFERENCE[..20], None),
        fastq("b", &wrong, Some(5)),
    ];
    let r = assemble(&inputs, &o).unwrap();
    source_audit(&r, &inputs, &o);
    assert!(r.contigs.iter().any(|c| c.sequence == REFERENCE));
}
#[test]
fn display_reversal_preserves_original_strand_support_and_original_call_indices() {
    let a = fastq("a", &REFERENCE[..20], None);
    let mut b = fastq("b", &REFERENCE[9..], None);
    b.molecule = collagene_core::annotations::reverse_complement(&b.molecule).unwrap();
    let mut o = options();
    o.both_strands = true;
    o.require_both_strands = true;
    let inputs = vec![a, b];
    let r = assemble(&inputs, &o).unwrap();
    source_audit(&r, &inputs, &o);
    assert!(r.reads[1].source_reverse);
    assert!(
        r.contigs
            .iter()
            .all(|c| c.sequence.bytes().all(|b| b == b'N')),
        "display reversal invented opposite original-strand support"
    );
}
#[test]
fn actual_reverse_reads_support_only_their_observed_common_interval() {
    let inputs = vec![
        plain("a", &REFERENCE[..20]),
        plain("b", &reverse(&REFERENCE[9..])),
    ];
    let mut o = options();
    o.both_strands = true;
    o.require_both_strands = true;
    o.min_support = 2;
    let r = assemble(&inputs, &o).unwrap();
    source_audit(&r, &inputs, &o);
    let expected = format!(
        "{}{}{}",
        "N".repeat(9),
        &REFERENCE[9..20],
        "N".repeat(REFERENCE.len() - 20)
    );
    assert!(
        r.contigs
            .iter()
            .any(|c| c.sequence == expected || c.sequence == reverse(&expected)),
        "{r:?}"
    );
}
#[test]
fn duplicate_plain_read_evidence_never_creates_two_record_support() {
    let inputs = vec![plain("a", REFERENCE), plain("duplicate", REFERENCE)];
    let mut o = options();
    o.min_support = 2;
    let r = assemble(&inputs, &o).unwrap();
    source_audit(&r, &inputs, &o);
    assert!(r.reads.iter().all(|x| x.duplicate_group.is_some()));
    assert!(
        r.contigs
            .iter()
            .all(|c| c.sequence.bytes().all(|b| b == b'N'))
    );
}
#[test]
fn indel_overlaps_retain_real_gap_evidence_and_never_turn_conflict_into_a_certain_base() {
    let reference = "AGTCACGATGCTAGCTACGTTGACCTGATCGTACAGGTCA";
    let a = &reference[..25];
    let mut b = reference[10..35].to_string();
    b.insert(8, 'A');
    let inputs = vec![plain("a", a), plain("b", &b), plain("c", &reference[15..])];
    let o = options();
    let r = assemble(&inputs, &o).unwrap();
    source_audit(&r, &inputs, &o);
    assert!(r.overlaps.iter().any(
        |e| e.accepted && (e.reference_aligned.contains('-') || e.query_aligned.contains('-'))
    ));
    assert!(
        r.contigs
            .iter()
            .any(|c| c.placements.len() >= 2 && c.columns.iter().any(|x| !x.gap_reads.is_empty()))
    );
    for col in r.contigs.iter().flat_map(|c| &c.columns) {
        if !col.eligible_gap_reads.is_empty() && col.observations.iter().any(|v| v.eligible) {
            assert_eq!(col.base, Some('N'));
            assert_eq!(col.state, ConsensusState::GapConflict);
        }
    }
}
fn exhaustive_score(a: &[u8], b: &[u8], contained: bool) -> Option<(i32, u8)> {
    fn walk(
        a: &[u8],
        b: &[u8],
        i: usize,
        j: usize,
        last: u8,
        score: i32,
        contained: bool,
        best: &mut i32,
        count: &mut u8,
    ) {
        let endpoint = if contained {
            j == b.len() && last == 0
        } else {
            i == a.len() && j < b.len()
        };
        if endpoint {
            if score > *best {
                *best = score;
                *count = 1;
            } else if score == *best {
                *count = (*count + 1).min(2);
            }
        }
        if i < a.len() && j < b.len() {
            walk(
                a,
                b,
                i + 1,
                j + 1,
                0,
                score + if a[i] == b[j] { 2 } else { -3 },
                contained,
                best,
                count,
            );
        }
        if j < b.len() {
            walk(
                a,
                b,
                i,
                j + 1,
                1,
                score + if last == 1 { -1 } else { -4 },
                contained,
                best,
                count,
            );
        }
        if i < a.len() {
            walk(
                a,
                b,
                i + 1,
                j,
                2,
                score + if last == 2 { -1 } else { -4 },
                contained,
                best,
                count,
            );
        }
    }
    let (mut best, mut count) = (i32::MIN, 0);
    for start in 0..a.len() {
        walk(
            a,
            b,
            start + 1,
            1,
            0,
            if a[start] == b[0] { 2 } else { -3 },
            contained,
            &mut best,
            &mut count,
        );
    }
    if count == 0 {
        None
    } else {
        Some((best, count))
    }
}
#[test]
fn every_short_affine_boundary_score_and_tie_count_matches_an_independent_path_enumerator() {
    let mut strings = vec![];
    for n in 1..=4 {
        for bits in 0..1usize << n {
            strings.push(
                (0..n)
                    .map(|i| if bits & (1 << i) == 0 { 'A' } else { 'C' })
                    .collect::<String>(),
            );
        }
    }
    let mut o = options();
    o.min_overlap = 1;
    o.max_error_bp = 10_000;
    for a in &strings {
        for b in &strings {
            let inputs = vec![plain("a", a), plain("b", b)];
            let r = assemble(&inputs, &o).unwrap();
            for contained in [false, true] {
                let expected = exhaustive_score(a.as_bytes(), b.as_bytes(), contained);
                let edge = r.overlaps.iter().find(|e| {
                    e.from == 0
                        && e.to == 1
                        && !e.from_reverse
                        && !e.to_reverse
                        && (e.kind != OverlapKind::Extension) == contained
                });
                match expected {
                    Some((score, count)) if score > 0 => {
                        let edge = edge.unwrap_or_else(|| {
                            panic!(
                                "missing positive boundary class {a}/{b} containment={contained}"
                            )
                        });
                        assert_eq!(
                            (edge.score, edge.optimal_paths_saturated),
                            (score, count),
                            "{a}/{b} containment={contained}: {edge:?}"
                        );
                    }
                    _ => assert!(
                        !edge.is_some_and(|e| e.accepted),
                        "ineligible boundary was accepted: {a}/{b}"
                    ),
                }
            }
        }
    }
}
