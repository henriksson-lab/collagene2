//! Literal full-window source-coordinate and IUPAC oracle; no engine helpers.
use collagene_core::sequence::crispr::*;
use collagene_model::{Molecule, Topology};
fn molecule(s: &str, circle: bool) -> Molecule {
    Molecule::dna(
        "oracle reference",
        s,
        if circle {
            Topology::Circular
        } else {
            Topology::Linear
        },
    )
    .unwrap()
}
fn profile(side: PamSide) -> NucleaseProfile {
    NucleaseProfile {
        id: "declared-oracle".into(),
        version: "literal-1".into(),
        guide_length: 12,
        pam: "NG".into(),
        pam_side: side,
        nominal_top_offset: 5,
        nominal_bottom_offset: 8,
        custom: true,
        sources: vec![],
    }
}
fn mask(c: u8) -> u8 {
    match c {
        b'A' => 1,
        b'C' => 2,
        b'G' => 4,
        b'T' => 8,
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
        _ => panic!(),
    }
}
fn complement(c: u8) -> u8 {
    match c {
        b'A' => b'T',
        b'C' => b'G',
        b'G' => b'C',
        b'T' => b'A',
        b'R' => b'Y',
        b'Y' => b'R',
        b'S' => b'S',
        b'W' => b'W',
        b'K' => b'M',
        b'M' => b'K',
        b'B' => b'V',
        b'D' => b'H',
        b'H' => b'D',
        b'V' => b'B',
        b'N' => b'N',
        _ => panic!(),
    }
}
type Key = (
    usize,
    bool,
    Vec<usize>,
    Vec<usize>,
    usize,
    usize,
    bool,
    bool,
    Option<usize>,
    Option<usize>,
);
fn oracle(q: &str, refs: &[ReferenceInput], o: &OffTargetOptions) -> Vec<Key> {
    let mut out = vec![];
    let p = &o.profile;
    for (ri, r) in refs.iter().enumerate() {
        let m = &r.molecule;
        let s = m.sequence.as_str().as_bytes();
        let n = s.len();
        let circle = m.topology == Topology::Circular;
        let total = p.guide_length + p.pam.len();
        if n < total {
            continue;
        }
        for reverse in [false, true] {
            for start in 0..if circle { n } else { n - total + 1 } {
                let coord = |j: usize| {
                    if reverse {
                        n - 1 - (start + j) % n
                    } else {
                        (start + j) % n
                    }
                };
                let base = |j: usize| {
                    if reverse {
                        complement(s[coord(j)])
                    } else {
                        s[coord(j)]
                    }
                };
                let (go, po) = if p.pam_side == PamSide::ThreePrime {
                    (0, p.guide_length)
                } else {
                    (p.pam.len(), 0)
                };
                let mut possible = true;
                let mut certain = true;
                for (j, b) in p.pam.bytes().enumerate() {
                    let actual = mask(base(po + j));
                    let expected = mask(b);
                    possible &= actual & expected != 0;
                    certain &= actual & !expected == 0;
                }
                if !possible {
                    continue;
                }
                let mut min = 0;
                let mut max = 0;
                for (j, b) in q.bytes().enumerate() {
                    let a = mask(base(go + j));
                    let b = mask(b);
                    min += usize::from(a & b == 0);
                    max += usize::from(a != b);
                }
                if min > o.max_mismatches {
                    continue;
                }
                let boundary = |offset: i32| {
                    let at = start as i64 + go as i64 + i64::from(offset);
                    let actual = if reverse { n as i64 - at } else { at };
                    if circle {
                        Some(actual.rem_euclid(n as i64) as usize)
                    } else {
                        (0..=n as i64).contains(&actual).then_some(actual as usize)
                    }
                };
                let (top, bottom) = if reverse {
                    (
                        boundary(p.nominal_bottom_offset),
                        boundary(p.nominal_top_offset),
                    )
                } else {
                    (
                        boundary(p.nominal_top_offset),
                        boundary(p.nominal_bottom_offset),
                    )
                };
                out.push((
                    ri,
                    reverse,
                    (0..p.guide_length).map(|j| coord(go + j)).collect(),
                    (0..p.pam.len()).map(|j| coord(po + j)).collect(),
                    min,
                    max,
                    certain,
                    certain && max <= o.max_mismatches,
                    top,
                    bottom,
                ));
            }
        }
    }
    out.sort();
    out
}
fn actual(report: &OffTargetReport) -> Vec<Key> {
    let mut keys = report
        .hits
        .iter()
        .map(|h| {
            (
                h.reference_index,
                h.reverse,
                h.guide_positions.clone(),
                h.pam_positions.clone(),
                h.min_mismatches,
                h.max_mismatches,
                h.pam_certain,
                h.definite,
                h.nominal_top_cut,
                h.nominal_bottom_cut,
            )
        })
        .collect::<Vec<_>>();
    keys.sort();
    keys
}
#[test]
fn indexed_and_ambiguous_fallback_search_equal_full_iupac_window_enumeration() {
    let query = GuideQuery {
        id: "independent".into(),
        sequence: "ACGTACGTACGT".into(),
    };
    for side in [PamSide::ThreePrime, PamSide::FivePrime] {
        for case in 0..20 {
            let mut spacer = query.sequence.as_bytes().to_vec();
            if case % 3 != 0 {
                let at = case % 12;
                spacer[at] = if case % 4 == 0 {
                    b'N'
                } else if spacer[at] == b'A' {
                    b'C'
                } else {
                    b'A'
                };
            }
            let spacer = String::from_utf8(spacer).unwrap();
            let pam = if case % 5 == 0 { "RG" } else { "AG" };
            let context = if side == PamSide::ThreePrime {
                format!("TTAC{spacer}{pam}CGATG")
            } else {
                format!("TTAC{pam}{spacer}CGATG")
            };
            for circle in [false, true] {
                for origin in (0..if circle { context.len() } else { 1 }).step_by(7) {
                    let rotated = format!("{}{}", &context[origin..], &context[..origin]);
                    let refs = vec![ReferenceInput {
                        id: "a".into(),
                        molecule: molecule(&rotated, circle),
                        build: Some("synthetic-selected-set".into()),
                        accession: Some("fixture-v1".into()),
                    }];
                    for distance in 0..=6 {
                        let o = OffTargetOptions {
                            profile: profile(side),
                            max_mismatches: distance,
                            ..Default::default()
                        };
                        let r = search_off_targets(&query, &refs, &o).unwrap();
                        assert!(
                            r.complete && r.exhaustive_for_selected_records,
                            "{side:?} case{case} circle{circle} origin{origin} d{distance}: {:?}",
                            r.limits
                        );
                        assert_eq!(
                            actual(&r),
                            oracle(&query.sequence, &refs, &o),
                            "{side:?} case{case} circle{circle} origin{origin} d{distance}"
                        );
                    }
                }
            }
        }
    }
}
#[test]
fn exact_intended_role_does_not_hide_equal_dna_at_other_records() {
    let text = "TTTTACGTACGTACGTAGCCCC";
    let refs = vec![
        ReferenceInput {
            id: "a".into(),
            molecule: molecule(text, false),
            build: Some("named-build".into()),
            accession: None,
        },
        ReferenceInput {
            id: "b".into(),
            molecule: molecule(text, false),
            build: Some("named-build".into()),
            accession: None,
        },
    ];
    let query = GuideQuery {
        id: "guide".into(),
        sequence: "ACGTACGTACGT".into(),
    };
    let o = OffTargetOptions {
        profile: profile(PamSide::ThreePrime),
        max_mismatches: 0,
        intended: Some(IntendedTarget {
            reference_id: "a".into(),
            reference_revision: refs[0].molecule.revision(),
            guide_positions: (4..16).collect(),
            pam_positions: vec![16, 17],
            reverse: false,
        }),
        ..Default::default()
    };
    let r = search_off_targets(&query, &refs, &o).unwrap();
    assert!(r.complete);
    assert_eq!(r.hits.iter().filter(|h| h.is_intended).count(), 1);
    assert!(r.hits.iter().any(|h| h.reference_index == 1
        && h.guide_positions == (4..16).collect::<Vec<_>>()
        && !h.is_intended));
    let mut changed = refs.clone();
    changed[0].build = Some("a different declared build".into());
    assert!(r.validate_bindings_for(&query, &changed, &o).is_err());
}
#[test]
fn guide_gc_filter_uses_exact_rational_bounds_not_the_display_floor() {
    let m = molecule("CAAAAAAAAAAAAG", false);
    let o = GuideOptions {
        profile: profile(PamSide::ThreePrime),
        min_gc_bp: 0,
        max_gc_bp: 833,
        ..Default::default()
    };
    let r = design_guides(&m, &o).unwrap();
    let c = r
        .candidates
        .iter()
        .find(|c| c.guide_positions == (0..12).collect::<Vec<_>>() && !c.reverse)
        .unwrap();
    assert_eq!(c.gc_bp, Some(833));
    assert!(!c.passes_filters);
    let mut relaxed = o.clone();
    relaxed.max_gc_bp = 834;
    let r = design_guides(&m, &relaxed).unwrap();
    assert!(
        r.candidates
            .iter()
            .any(|c| c.guide_positions == (0..12).collect::<Vec<_>>()
                && !c.reverse
                && c.passes_filters)
    );
}
#[test]
fn every_resource_cap_invalidates_exhaustive_absence_claim() {
    let query = GuideQuery {
        id: "repetitive".into(),
        sequence: "AAAAAAAAAAAA".into(),
    };
    let refs = vec![ReferenceInput {
        id: "repeat".into(),
        molecule: molecule("AAAAAAAAAAAAAG".repeat(10).as_str(), true),
        build: None,
        accession: None,
    }];
    for kind in 0..3 {
        let mut o = OffTargetOptions {
            profile: profile(PamSide::ThreePrime),
            max_mismatches: 2,
            ..Default::default()
        };
        match kind {
            0 => o.max_comparisons = 1,
            1 => o.max_index_entries = 1,
            _ => o.max_hits = 1,
        }
        let r = search_off_targets(&query, &refs, &o).unwrap();
        assert!(!r.complete && !r.exhaustive_for_selected_records && !r.limits.is_empty());
    }
}
