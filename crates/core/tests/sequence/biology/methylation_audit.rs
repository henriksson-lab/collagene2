//! Independent literal context/coordinate checks; no engine matching or span helpers.
use collagene_core::methylation::*;
use collagene_model::{Molecule, Sequence, Topology};
fn mol(text: &str, circular: bool) -> Molecule {
    Molecule {
        name: "methylation oracle".into(),
        sequence: Sequence::dna(text).unwrap(),
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
fn options(enzyme: &str) -> MethylationOptions {
    MethylationOptions {
        enzymes: vec![enzyme.into()],
        dam: MethylationState::Unmodified,
        dcm: MethylationState::Unmodified,
        cpg: MethylationState::Unmodified,
        ..Default::default()
    }
}
fn rc(text: &str) -> String {
    text.bytes()
        .rev()
        .map(|b| match b {
            b'A' => 'T',
            b'C' => 'G',
            b'G' => 'C',
            b'T' => 'A',
            _ => panic!(),
        })
        .collect()
}
fn literal_matches(text: &str, pattern: &str, circular: bool) -> Vec<usize> {
    let b = text.as_bytes();
    (0..if circular {
        b.len()
    } else {
        b.len().saturating_sub(pattern.len()) + 1
    })
        .filter(|&i| (0..pattern.len()).all(|j| b[(i + j) % b.len()] == pattern.as_bytes()[j]))
        .collect()
}
fn overlapping(
    text: &str,
    site: usize,
    len: usize,
    pattern: &str,
    circular: bool,
) -> Vec<Vec<usize>> {
    let n = text.len();
    let positions = (0..len).map(|j| (site + j) % n).collect::<Vec<_>>();
    literal_matches(text, pattern, circular)
        .into_iter()
        .filter_map(|i| {
            let p = (0..pattern.len()).map(|j| (i + j) % n).collect::<Vec<_>>();
            p.iter().any(|x| positions.contains(x)).then_some(p)
        })
        .collect()
}
fn check_xba(text: &str, circular: bool) {
    let m = mol(text, circular);
    let mut o = options("XbaI");
    o.dam = MethylationState::FullyModified;
    let r = assess(&m, &o).unwrap();
    assert!(r.complete);
    let sites = literal_matches(text, "TCTAGA", circular);
    assert_eq!(r.sites.len(), sites.len());
    for (s, start) in r.sites.iter().zip(sites) {
        assert_eq!(s.recognition_start, start);
        assert!(s.certain);
        let expected = overlapping(text, start, 6, "GATC", circular);
        let effect = s
            .effects
            .iter()
            .find(|e| e.category == Category::Dam)
            .unwrap();
        let mut actual = effect
            .overlaps
            .iter()
            .map(|x| x.positions.clone())
            .collect::<Vec<_>>();
        actual.sort();
        let mut wanted = expected.clone();
        wanted.sort();
        assert_eq!(actual, wanted, "{text} site{start}");
        assert_eq!(
            effect.outcome,
            if expected.is_empty() {
                Outcome::PermittedUnderProfile
            } else {
                Outcome::Blocked
            }
        );
        assert_eq!(s.outcome, effect.outcome);
        for motif in &effect.overlaps {
            assert!(motif.certain);
            assert_eq!(motif.motif, "GATC");
            let mut bases = motif
                .modified_bases
                .iter()
                .map(|b| (b.position, b.reverse, b.base))
                .collect::<Vec<_>>();
            bases.sort();
            let mut expected = vec![
                (motif.positions[1], false, 'A'),
                (motif.positions[2], true, 'A'),
            ];
            expected.sort();
            assert_eq!(bases, expected);
        }
        assert_eq!(
            s.top_cut,
            Some(if circular {
                (start + 1) % text.len()
            } else {
                start + 1
            })
        );
        assert_eq!(
            s.bottom_cut,
            Some(if circular {
                (start + 5) % text.len()
            } else {
                start + 5
            })
        );
    }
}
#[test]
fn every_concrete_xbai_flank_and_origin_has_literal_dam_context_and_strand_coordinates() {
    for a in b"ACGT" {
        for b in b"ACGT" {
            for c in b"ACGT" {
                for d in b"ACGT" {
                    let text = format!(
                        "{}{}TCTAGA{}{}",
                        *a as char, *b as char, *c as char, *d as char
                    );
                    check_xba(&text, false);
                    check_xba(&rc(&text), false);
                    for origin in 0..text.len() {
                        let rotated = format!("{}{}", &text[origin..], &text[..origin]);
                        check_xba(&rotated, true);
                    }
                }
            }
        }
    }
}
#[test]
fn fully_cpg_modified_internal_sites_retain_exact_top_and_bottom_cytosines() {
    for (enzyme, text, pattern) in [
        ("NotI", "AAGCGGCCGCAA", "GCGGCCGC"),
        ("SmaI", "AACCCGGGAA", "CCCGGG"),
    ] {
        for reverse in [false, true] {
            let source = if reverse { rc(text) } else { text.into() };
            for origin in 0..source.len() {
                let rotated = format!("{}{}", &source[origin..], &source[..origin]);
                let m = mol(&rotated, true);
                let mut o = options(enzyme);
                o.cpg = MethylationState::FullyModified;
                let r = assess(&m, &o).unwrap();
                assert!(r.complete);
                let starts = literal_matches(&rotated, pattern, true);
                assert_eq!(r.sites.len(), starts.len());
                for (s, start) in r.sites.iter().zip(starts) {
                    assert_eq!(s.recognition_start, start);
                    let e = s
                        .effects
                        .iter()
                        .find(|e| e.category == Category::CpG)
                        .unwrap();
                    assert_eq!(e.outcome, Outcome::Blocked);
                    let mut wanted = overlapping(&rotated, start, pattern.len(), "CG", true);
                    wanted.sort();
                    let mut actual = e
                        .overlaps
                        .iter()
                        .map(|x| x.positions.clone())
                        .collect::<Vec<_>>();
                    actual.sort();
                    assert_eq!(actual, wanted);
                    for motif in &e.overlaps {
                        let mut bases = motif
                            .modified_bases
                            .iter()
                            .map(|b| (b.position, b.reverse, b.base))
                            .collect::<Vec<_>>();
                        bases.sort();
                        let mut expected = vec![
                            (motif.positions[0], false, 'C'),
                            (motif.positions[1], true, 'C'),
                        ];
                        expected.sort();
                        assert_eq!(bases, expected);
                    }
                }
            }
        }
    }
}
#[test]
fn unknown_and_ambiguous_dam_contexts_never_become_permitted() {
    for text in ["NNTCTAGACC", "AATCTAGANN", "GATCTAGACC", "AATCTAGATC"] {
        let m = mol(text, false);
        for state in [
            MethylationState::FullyModified,
            MethylationState::Unknown,
            MethylationState::HemiOrMixed,
        ] {
            let mut o = options("XbaI");
            o.dam = state;
            let r = assess(&m, &o).unwrap();
            assert_eq!(r.sites.len(), 1);
            assert_ne!(
                r.sites[0].outcome,
                Outcome::PermittedUnderProfile,
                "{text} {state:?}"
            );
        }
    }
}
#[test]
fn bamhi_internal_dam_motif_does_not_imply_inhibition() {
    let m = mol("AAGGATCCAA", false);
    for state in [MethylationState::FullyModified, MethylationState::Unknown] {
        let mut o = options("BamHI");
        o.dam = state;
        o.dcm = state;
        o.cpg = state;
        let r = assess(&m, &o).unwrap();
        assert!(r.complete);
        assert_eq!(r.sites.len(), 1);
        assert_eq!(r.sites[0].outcome, Outcome::PermittedUnderProfile);
    }
}
#[test]
fn type_iis_orientations_and_circular_cut_boundaries_match_literal_arithmetic() {
    for (name, pattern, top, bottom) in [
        ("BsaI", "GGTCTC", 7i64, 11i64),
        ("BsmBI", "CGTCTC", 7, 11),
        ("BbsI", "GAAGAC", 8, 12),
        ("SapI", "GCTCTTC", 8, 11),
    ] {
        let forward = format!("AAAAAAAAAAAA{pattern}AAAAAAAAAAAA");
        for reverse in [false, true] {
            let text = if reverse {
                rc(&forward)
            } else {
                forward.clone()
            };
            let recognition = if reverse { rc(pattern) } else { pattern.into() };
            for circular in [false, true] {
                for origin in 0..if circular { text.len() } else { 1 } {
                    let rotated = format!("{}{}", &text[origin..], &text[..origin]);
                    let r = assess(&mol(&rotated, circular), &options(name)).unwrap();
                    assert!(r.complete);
                    let starts = literal_matches(&rotated, &recognition, circular);
                    assert_eq!(starts.len(), 1);
                    assert_eq!(r.sites.len(), 1);
                    let s = &r.sites[0];
                    assert_eq!((s.recognition_start, s.reverse), (starts[0], reverse));
                    let at = starts[0] as i64;
                    let n = text.len() as i64;
                    let (a, b) = if reverse {
                        (
                            at + pattern.len() as i64 - bottom,
                            at + pattern.len() as i64 - top,
                        )
                    } else {
                        (at + top, at + bottom)
                    };
                    let cv = |x: i64| {
                        if circular {
                            Some(x.rem_euclid(n) as usize)
                        } else {
                            (0..=n).contains(&x).then_some(x as usize)
                        }
                    };
                    assert_eq!((s.top_cut, s.bottom_cut), (cv(a), cv(b)));
                    assert_eq!(s.outcome, Outcome::PermittedUnderProfile);
                }
            }
        }
    }
}
#[test]
fn dcm_both_strand_targets_are_retained_even_when_supplier_rule_is_unresolved() {
    for middle in ['A', 'T'] {
        let text = format!("AAAAGAGACC{middle}GGAAAA");
        for reverse in [false, true] {
            let source = if reverse { rc(&text) } else { text.clone() };
            for origin in 0..source.len() {
                let rotated = format!("{}{}", &source[origin..], &source[..origin]);
                let r = assess(&mol(&rotated, true), &options("BsaI")).unwrap();
                assert!(r.complete);
                assert_eq!(r.sites.len(), 1);
                let s = &r.sites[0];
                let effect = s
                    .effects
                    .iter()
                    .find(|e| e.category == Category::Dcm)
                    .unwrap();
                let mut wanted = overlapping(&rotated, s.recognition_start, 6, "CCAGG", true);
                wanted.extend(overlapping(&rotated, s.recognition_start, 6, "CCTGG", true));
                wanted.sort();
                let mut actual = effect
                    .overlaps
                    .iter()
                    .map(|x| x.positions.clone())
                    .collect::<Vec<_>>();
                actual.sort();
                assert_eq!(actual, wanted);
                assert!(!actual.is_empty());
                for motif in &effect.overlaps {
                    let mut bases = motif
                        .modified_bases
                        .iter()
                        .map(|b| (b.position, b.reverse, b.base))
                        .collect::<Vec<_>>();
                    bases.sort();
                    let mut wanted = vec![
                        (motif.positions[1], false, 'C'),
                        (motif.positions[3], true, 'C'),
                    ];
                    wanted.sort();
                    assert_eq!(bases, wanted);
                }
            }
        }
    }
}
