//! Independent literal-grid oracle. Does not use the implementation's index,
//! matching, span, reverse-complement or coordinate helpers.
use collagene_core::{Alphabet, Molecule, Sequence, Topology, dotplot::*};
fn molecule(alphabet: Alphabet, text: &str, circular: bool) -> Molecule {
    Molecule {
        name: "independent coordinate fixture".into(),
        sequence: Sequence::new(alphabet, text).unwrap(),
        topology: if circular {
            Topology::Circular
        } else {
            Topology::Linear
        },
        features: vec![],
        metadata: Default::default(),
        original_genbank: None,
    }
}
fn mask(b: u8) -> u32 {
    match b {
        b'A' => 1,
        b'C' => 2,
        b'G' => 4,
        b'T' | b'U' => 8,
        b'R' => 5,
        b'Y' => 10,
        b'S' => 6,
        b'W' => 9,
        b'K' => 12,
        b'M' => 3,
        b'B' => 14,
        b'D' => 13,
        b'H' => 11,
        b'V' => 7,
        b'N' => 15,
        _ => panic!("unknown nucleotide"),
    }
}
fn complement(b: u8) -> u8 {
    match b {
        b'A' => b'T',
        b'C' => b'G',
        b'G' => b'C',
        b'T' | b'U' => b'A',
        b'R' => b'Y',
        b'Y' => b'R',
        b'S' => b'S',
        b'W' => b'W',
        b'K' => b'M',
        b'M' => b'K',
        b'B' => b'V',
        b'V' => b'B',
        b'D' => b'H',
        b'H' => b'D',
        b'N' => b'N',
        _ => panic!("unknown nucleotide"),
    }
}
type Key = (usize, usize, bool, usize, usize, bool);
fn oracle(a: &str, b: &str, ac: bool, bc: bool, o: &DotPlotOptions) -> Vec<Key> {
    let mut result = vec![];
    let an = if ac && o.wrap_circular {
        a.len()
    } else {
        a.len() - o.window + 1
    };
    let bn = if bc && o.wrap_circular {
        b.len()
    } else {
        b.len() - o.window + 1
    };
    for reverse in [false, true]
        .into_iter()
        .take(if o.both_strands { 2 } else { 1 })
    {
        for x in (0..an).step_by(o.stride) {
            for y in (0..bn).step_by(o.stride) {
                let mut possible = 0;
                let mut sure = 0;
                for j in 0..o.window {
                    let ai = (x + j) % a.len();
                    let bi = if reverse {
                        b.len() - 1 - (y + j) % b.len()
                    } else {
                        (y + j) % b.len()
                    };
                    let am = mask(a.as_bytes()[ai]);
                    let bm = mask(if reverse {
                        complement(b.as_bytes()[bi])
                    } else {
                        b.as_bytes()[bi]
                    });
                    possible += usize::from(am & bm != 0);
                    sure += usize::from(am == bm && am.count_ones() == 1);
                }
                if possible >= o.min_matches {
                    result.push((
                        x,
                        if reverse { b.len() - 1 - y } else { y },
                        reverse,
                        possible,
                        sure,
                        sure >= o.min_matches,
                    ));
                }
            }
        }
    }
    result.sort();
    result
}
fn compare(a: &str, b: &str, ac: bool, bc: bool, o: &DotPlotOptions) {
    let reference = molecule(Alphabet::Dna, a, ac);
    let query = molecule(Alphabet::Dna, b, bc);
    let report = analyze(&reference, &query, o).unwrap();
    report.validate_bindings_for(&reference, &query, o).unwrap();
    assert!(
        report.grid_complete,
        "bounded small oracle should complete {report:?}"
    );
    assert_eq!(report.sampled, o.stride > 1);
    assert_eq!(report.complete, o.stride == 1);
    let mut actual = report
        .hits
        .iter()
        .map(|h| {
            (
                h.reference_start,
                h.query_start,
                h.reverse,
                h.compatible_matches,
                h.certain_matches,
                h.certain,
            )
        })
        .collect::<Vec<_>>();
    actual.sort();
    assert_eq!(
        actual,
        oracle(a, b, ac, bc, o),
        "{a} {b} circular={ac}/{bc} options={o:?}"
    );
    for hit in &report.hits {
        let mut ai = vec![];
        for s in &hit.reference_spans {
            assert!(s.start < s.end && s.end <= a.len());
            ai.extend(s.start..s.end);
        }
        let mut bi = vec![];
        for s in &hit.query_spans {
            assert!(s.start < s.end && s.end <= b.len());
            if hit.reverse {
                bi.extend((s.start..s.end).rev());
            } else {
                bi.extend(s.start..s.end);
            }
        }
        assert_eq!(
            ai,
            (0..o.window)
                .map(|j| (hit.reference_start + j) % a.len())
                .collect::<Vec<_>>()
        );
        assert_eq!(
            bi,
            (0..o.window)
                .map(|j| if hit.reverse {
                    (hit.query_start + b.len() - j % b.len()) % b.len()
                } else {
                    (hit.query_start + j) % b.len()
                })
                .collect::<Vec<_>>()
        );
    }
}
#[test]
fn every_short_ambiguous_pair_matches_literal_grid_in_all_topologies() {
    let mut texts = vec![];
    for a in b"ACN" {
        for b in b"ACN" {
            for c in b"ACN" {
                texts.push(String::from_utf8(vec![*a, *b, *c]).unwrap());
            }
        }
    }
    for a in &texts {
        for b in &texts {
            for ac in [false, true] {
                for bc in [false, true] {
                    for window in 1..=3 {
                        for min_matches in 1..=window {
                            compare(
                                a,
                                b,
                                ac,
                                bc,
                                &DotPlotOptions {
                                    window,
                                    min_matches,
                                    ..Default::default()
                                },
                            );
                        }
                    }
                }
            }
        }
    }
}
#[test]
fn nonpalindromic_inversions_strides_and_origin_windows_have_exact_original_coordinates() {
    for (a, b) in [("ACGTCGAT", "ATCGACGT"), ("AGCT", "AGCT"), ("ARYN", "NYRT")] {
        for ac in [false, true] {
            for bc in [false, true] {
                for stride in [1, 2, 3] {
                    for window in 1..=4 {
                        compare(
                            a,
                            b,
                            ac,
                            bc,
                            &DotPlotOptions {
                                window,
                                min_matches: window,
                                stride,
                                ..Default::default()
                            },
                        );
                    }
                }
            }
        }
    }
}
#[test]
fn literal_rna_reverse_complement_uses_uracil_without_relabeling_the_source() {
    let r = molecule(Alphabet::Rna, "ACGU", true);
    let q = molecule(Alphabet::Rna, "ACGU", false);
    let o = DotPlotOptions {
        window: 3,
        min_matches: 3,
        ..Default::default()
    };
    let report = analyze(&r, &q, &o).unwrap();
    let mut actual = report
        .hits
        .iter()
        .map(|h| {
            (
                h.reference_start,
                h.query_start,
                h.reverse,
                h.compatible_matches,
                h.certain_matches,
                h.certain,
            )
        })
        .collect::<Vec<_>>();
    actual.sort();
    assert_eq!(actual, oracle("ACGU", "ACGU", true, false, &o));
    report.validate_bindings_for(&r, &q, &o).unwrap();
}
#[test]
fn report_binding_rejects_source_query_and_options_substitution() {
    let a = molecule(Alphabet::Dna, "ACGT", false);
    let b = molecule(Alphabet::Dna, "ACGT", false);
    let o = DotPlotOptions {
        window: 2,
        min_matches: 2,
        ..Default::default()
    };
    let r = analyze(&a, &b, &o).unwrap();
    let changed = molecule(Alphabet::Dna, "AAAA", false);
    assert!(r.validate_bindings_for(&changed, &b, &o).is_err());
    assert!(r.validate_bindings_for(&a, &changed, &o).is_err());
    let mut altered = o.clone();
    altered.stride = 2;
    assert!(r.validate_bindings_for(&a, &b, &altered).is_err());
}
#[test]
fn reported_incomplete_limits_do_not_claim_absence_or_complete_sampling() {
    let a = molecule(Alphabet::Dna, "AAAAAAAAAA", true);
    let mut o = DotPlotOptions {
        window: 2,
        min_matches: 2,
        max_hits: 3,
        ..Default::default()
    };
    let r = analyze(&a, &a, &o).unwrap();
    assert_eq!(r.hits.len(), 3);
    assert!(!r.complete && !r.grid_complete);
    assert!(!r.limits.is_empty());
    o.max_hits = 1000;
    o.stride = 2;
    let r = analyze(&a, &a, &o).unwrap();
    assert!(r.grid_complete && r.sampled && !r.complete);
    assert_eq!(r.hits.len(), 25);
    o.stride = 1;
    o.max_work = 1;
    let r = analyze(&a, &a, &o).unwrap();
    assert!(!r.complete && !r.grid_complete);
    assert!(!r.limits.is_empty());
}
fn amino_mask(b: u8) -> u32 {
    let concrete = b"ACDEFGHIKLMNPQRSTVWYUO*";
    let symbols: &[u8] = match b {
        b'B' => b"DN",
        b'Z' => b"EQ",
        b'J' => b"IL",
        b'X' => b"ACDEFGHIKLMNPQRSTVWYUO",
        _ => std::slice::from_ref(&b),
    };
    symbols.iter().fold(0, |m, s| {
        m | (1 << concrete.iter().position(|x| x == s).unwrap())
    })
}
#[test]
fn all_protein_symbols_use_independent_sets_and_never_promote_unknown_identity() {
    let text = "ACDEFGHIKLMNPQRSTVWYBXZJUO*";
    let p = molecule(Alphabet::Protein, text, false);
    let o = DotPlotOptions {
        window: 1,
        min_matches: 1,
        both_strands: false,
        ..Default::default()
    };
    let r = analyze(&p, &p, &o).unwrap();
    assert!(r.complete);
    let mut expected = vec![];
    for (i, a) in text.bytes().enumerate() {
        for (j, b) in text.bytes().enumerate() {
            let ma = amino_mask(a);
            let mb = amino_mask(b);
            if ma & mb != 0 {
                let sure = ma == mb && ma.count_ones() == 1;
                expected.push((i, j, false, 1, usize::from(sure), sure));
            }
        }
    }
    let mut actual = r
        .hits
        .iter()
        .map(|h| {
            (
                h.reference_start,
                h.query_start,
                h.reverse,
                h.compatible_matches,
                h.certain_matches,
                h.certain,
            )
        })
        .collect::<Vec<_>>();
    actual.sort();
    expected.sort();
    assert_eq!(actual, expected);
    let x = text.find('X').unwrap();
    assert!(
        !r.hits
            .iter()
            .find(|h| h.reference_start == x && h.query_start == x)
            .unwrap()
            .certain
    );
}
#[test]
fn long_exact_sequence_uses_an_index_without_omitting_proven_nonmatches() {
    let mut state = 27u32;
    let text = (0..20_000)
        .map(|_| {
            state = state.wrapping_mul(1664525).wrapping_add(1013904223);
            b"ACGT"[((state >> 24) & 3) as usize] as char
        })
        .collect::<String>();
    let p = molecule(Alphabet::Dna, &text, false);
    let o = DotPlotOptions {
        window: 25,
        min_matches: 25,
        both_strands: false,
        ..Default::default()
    };
    let r = analyze(&p, &p, &o).unwrap();
    assert!(r.complete);
    assert_eq!(r.examined_windows, (text.len() - 24).pow(2) as u64);
    assert!(r.work_used < o.max_work);
    // Diagonal matches are mandatory even if the generated source also repeats.
    for i in 0..text.len() - 24 {
        assert!(
            r.hits
                .binary_search_by_key(&(i, i, false), |h| (
                    h.reference_start,
                    h.query_start,
                    h.reverse
                ))
                .is_ok()
        );
    }
}

#[test]
fn standalone_svg_retains_the_exact_source_options_and_engine_envelope() {
    let a = molecule(Alphabet::Dna, "ACGTACGT", true);
    let o = DotPlotOptions {
        window: 3,
        min_matches: 2,
        stride: 2,
        ..Default::default()
    };
    let r = analyze(&a, &a, &o).unwrap();
    let svg = export_svg(&r).unwrap();
    let xml = roxmltree::Document::parse(&svg).unwrap();
    let text = xml
        .descendants()
        .find(|n| n.has_tag_name("metadata") && n.attribute("id") == Some("collagene-dotplot"))
        .unwrap()
        .text()
        .unwrap();
    let metadata: serde_json::Value = serde_json::from_str(text).unwrap();
    assert_eq!(metadata["options"], serde_json::to_value(&o).unwrap());
    assert_eq!(metadata["options_revision"], r.options_revision);
    assert_eq!(metadata["engine_version"], r.engine_version);
    assert_eq!(metadata["reference_revision"], a.revision());
    assert_eq!(metadata["query_revision"], a.revision());
    assert_eq!(metadata["sampled"], true);
    assert_eq!(metadata["complete"], false);
    assert!(
        metadata["authority"]
            .as_str()
            .unwrap()
            .contains("Unverified")
    );
}
