//! Literal geometry and independent exhaustive Hamming fixtures.
use collagene_core::{crispr::*, *};
const GUIDE: &str = "ACGTCAGTACGATCGATCGA";
fn m(s: &str, t: Topology) -> Molecule {
    Molecule::dna("source", s, t).unwrap()
}
fn rc(s: &str) -> String {
    s.bytes()
        .rev()
        .map(|b| match b {
            b'A' => 'T',
            b'C' => 'G',
            b'G' => 'C',
            b'T' => 'A',
            _ => 'N',
        })
        .collect()
}
fn reference(id: &str, s: &str, t: Topology) -> ReferenceInput {
    ReferenceInput {
        id: id.into(),
        molecule: m(s, t),
        build: Some("explicit-fixture-build".into()),
        accession: Some(format!("fixture:{id}")),
    }
}
#[test]
fn spcas9_coordinates_and_nominal_boundaries_follow_both_actual_orientations_and_origins() {
    assert_eq!(GUIDE.len(), 20);
    let text = format!("AAA{GUIDE}TGGCCC");
    let source = m(&text, Topology::Linear);
    let o = GuideOptions::default();
    let report = design_guides(&source, &o).unwrap();
    let guide = report
        .candidates
        .iter()
        .find(|c| !c.reverse && c.guide_positions[0] == 3)
        .unwrap();
    assert_eq!(guide.guide_sequence, GUIDE);
    assert_eq!(guide.pam_positions, vec![23, 24, 25]);
    assert_eq!(
        (guide.nominal_top_cut, guide.nominal_bottom_cut),
        (Some(20), Some(20))
    );
    validate_candidate(&source, &o, guide).unwrap();
    let reverse = m(&rc(&text), Topology::Linear);
    let report = design_guides(&reverse, &o).unwrap();
    let guide = report
        .candidates
        .iter()
        .find(|c| c.reverse && c.guide_sequence == GUIDE)
        .unwrap();
    assert_eq!(guide.guide_positions, (6..26).rev().collect::<Vec<_>>());
    assert_eq!(
        (guide.nominal_top_cut, guide.nominal_bottom_cut),
        (Some(9), Some(9))
    );
    for start in 0..text.len() {
        let rotated = format!("{}{}", &text[start..], &text[..start]);
        let report = design_guides(&m(&rotated, Topology::Circular), &o).unwrap();
        let g = report
            .candidates
            .iter()
            .find(|c| !c.reverse && c.guide_positions[0] == (3 + text.len() - start) % text.len())
            .unwrap();
        assert_eq!(g.guide_sequence, GUIDE);
        assert_eq!(
            g.nominal_top_cut,
            Some((20 + text.len() - start) % text.len())
        );
    }
}
#[test]
fn cas12a_declares_upstream_pam_and_swaps_nominal_strand_boundaries_on_reverse() {
    let guide = format!("{GUIDE}TCA");
    let text = format!("GGTTTA{guide}CC");
    let mut o = GuideOptions::default();
    o.profile = cas12a();
    let r = design_guides(&m(&text, Topology::Linear), &o).unwrap();
    let c = r
        .candidates
        .iter()
        .find(|c| !c.reverse && c.guide_sequence == guide)
        .unwrap();
    assert_eq!(c.guide_positions[0], 6);
    assert_eq!(c.pam_positions, vec![2, 3, 4, 5]);
    assert_eq!(
        (c.nominal_top_cut, c.nominal_bottom_cut),
        (Some(24), Some(29))
    );
    let r = design_guides(&m(&rc(&text), Topology::Linear), &o).unwrap();
    let c = r
        .candidates
        .iter()
        .find(|c| c.reverse && c.guide_sequence == guide)
        .unwrap();
    assert_eq!(
        (c.nominal_top_cut, c.nominal_bottom_cut),
        (Some(2), Some(7))
    );
    assert!(r.warnings.iter().any(|s| s.contains("heterogeneity")));
}
#[test]
fn ambiguous_calls_are_diagnostics_and_filter_changes_do_not_change_candidate_identity() {
    let source = m(&format!("AAA{GUIDE}TGGCCC"), Topology::Linear);
    let o = GuideOptions::default();
    let r = design_guides(&source, &o).unwrap();
    let c = r
        .candidates
        .iter()
        .find(|c| !c.reverse && c.guide_sequence == GUIDE)
        .unwrap();
    let mut strict = o.clone();
    strict.min_gc_bp = 10_000;
    let r = design_guides(&source, &strict).unwrap();
    let bad = r.candidates.iter().find(|x| x.id == c.id).unwrap();
    assert!(!bad.passes_filters);
    assert!(validate_candidate(&source, &strict, bad).is_err());
    let text = format!("AAA{}NTGGCCC", &GUIDE[..19]);
    let r = design_guides(&m(&text, Topology::Linear), &o).unwrap();
    assert!(
        r.candidates
            .iter()
            .any(|c| !c.concrete && c.gc_bp.is_none() && !c.passes_filters)
    );
    let text = format!("AAA{GUIDE}TGRCCC");
    let r = design_guides(&m(&text, Topology::Linear), &o).unwrap();
    let c = r
        .candidates
        .iter()
        .find(|c| !c.reverse && c.guide_sequence == GUIDE)
        .unwrap();
    assert!(!c.pam_certain && !c.passes_filters);
    let mut exact = GuideOptions::default();
    exact.profile = NucleaseProfile {
        id: "gc-boundary-fixture".into(),
        version: "literal-v1".into(),
        guide_length: 12,
        pam: "N".into(),
        pam_side: PamSide::ThreePrime,
        nominal_top_offset: 8,
        nominal_bottom_offset: 8,
        custom: true,
        sources: vec![],
    };
    exact.max_gc_bp = 833;
    let report = design_guides(&m("CAAAAAAAAAAAA", Topology::Linear), &exact).unwrap();
    let guide = report.candidates.iter().find(|c| !c.reverse).unwrap();
    assert_eq!(guide.gc_bp, Some(833));
    assert!(
        !guide.passes_filters,
        "1/12 is greater than833/10000 despite display flooring"
    );
}
fn oracle(reference: &ReferenceInput, query: &str, d: usize) -> Vec<(bool, Vec<usize>, usize)> {
    let mut out = vec![];
    let n = reference.molecule.sequence.len();
    let circular = reference.molecule.topology == Topology::Circular;
    for reverse in [false, true] {
        let oriented = if reverse {
            rc(reference.molecule.sequence.as_str())
        } else {
            reference.molecule.sequence.as_str().to_owned()
        };
        let count = if query.len() + 3 > n {
            0
        } else if circular {
            n
        } else {
            n - query.len() - 3 + 1
        };
        for start in 0..count {
            let base = |i: usize| oriented.as_bytes()[i % n];
            if base(start + query.len() + 1) != b'G' || base(start + query.len() + 2) != b'G' {
                continue;
            }
            let mismatches = query
                .bytes()
                .enumerate()
                .filter(|(j, b)| base(start + j) != *b)
                .count();
            if mismatches <= d {
                let positions = (0..query.len())
                    .map(|j| {
                        let i = (start + j) % n;
                        if reverse { n - 1 - i } else { i }
                    })
                    .collect();
                out.push((reverse, positions, mismatches));
            }
        }
    }
    out.sort();
    out
}
#[test]
fn indexed_pigeonhole_search_matches_an_independent_full_hamming_scan() {
    let mut mutated = GUIDE.as_bytes().to_vec();
    mutated[2] = b'T';
    mutated[11] = b'A';
    let text = format!(
        "CCAAA{GUIDE}TGGTT{}AGGAA{}",
        String::from_utf8(mutated).unwrap(),
        rc(&format!("{GUIDE}CGG"))
    );
    let q = GuideQuery {
        id: "guide".into(),
        sequence: GUIDE.into(),
    };
    for t in [Topology::Linear, Topology::Circular] {
        for shift in 0..text.len() {
            let rotated = format!("{}{}", &text[shift..], &text[..shift]);
            let refs = vec![reference("ref", &rotated, t)];
            for mismatches in 0..=3 {
                let o = OffTargetOptions {
                    max_mismatches: mismatches,
                    ..Default::default()
                };
                let r = search_off_targets(&q, &refs, &o).unwrap();
                assert!(
                    r.complete && r.exhaustive_for_selected_records,
                    "{:?}",
                    r.limits
                );
                assert!(r.index_used);
                let mut actual = r
                    .hits
                    .iter()
                    .map(|h| (h.reverse, h.guide_positions.clone(), h.min_mismatches))
                    .collect::<Vec<_>>();
                actual.sort();
                assert_eq!(
                    actual,
                    oracle(&refs[0], GUIDE, mismatches),
                    "shift{shift} mismatch{mismatches}"
                );
                assert!(
                    r.hits
                        .iter()
                        .all(|h| h.min_mismatches == h.max_mismatches && h.definite)
                );
            }
        }
    }
}
#[test]
fn ambiguity_fallback_retains_mismatch_ranges_and_exact_intended_locus_only() {
    let text = format!("{GUIDE}AGG{GUIDE}AGG");
    let refs = vec![
        reference("first", &text, Topology::Linear),
        reference(
            "same-DNA-elsewhere",
            &format!("{GUIDE}AGG"),
            Topology::Linear,
        ),
    ];
    let q = GuideQuery {
        id: "guide".into(),
        sequence: GUIDE.into(),
    };
    let intended = IntendedTarget {
        reference_id: "first".into(),
        reference_revision: refs[0].molecule.revision(),
        guide_positions: (0..20).collect(),
        pam_positions: (20..23).collect(),
        reverse: false,
    };
    let o = OffTargetOptions {
        max_mismatches: 0,
        intended: Some(intended),
        ..Default::default()
    };
    let r = search_off_targets(&q, &refs, &o).unwrap();
    assert_eq!(r.hits.iter().filter(|h| h.is_intended).count(), 1);
    assert_eq!(
        r.hits
            .iter()
            .filter(|h| !h.reverse && h.min_mismatches == 0)
            .count(),
        3
    );
    let mut text = format!("{GUIDE}AGG");
    text.replace_range(5..6, "N");
    let refs = vec![reference("ambiguous", &text, Topology::Linear)];
    let o = OffTargetOptions {
        max_mismatches: 0,
        ..Default::default()
    };
    let r = search_off_targets(&q, &refs, &o).unwrap();
    assert!(r.complete && !r.index_used);
    let h = r.hits.iter().find(|h| !h.reverse).unwrap();
    assert_eq!((h.min_mismatches, h.max_mismatches), (0, 1));
    assert!(!h.definite);
}
#[test]
fn bounds_and_source_provenance_never_turn_a_partial_search_into_absence() {
    let q = GuideQuery {
        id: "repeat".into(),
        sequence: "A".repeat(20),
    };
    let refs = vec![reference(
        "repeat",
        &format!("{}AGG{}AGG", "A".repeat(20), "A".repeat(20)),
        Topology::Linear,
    )];
    let mut o = OffTargetOptions::default();
    o.max_comparisons = 1;
    let r = search_off_targets(&q, &refs, &o).unwrap();
    assert!(!r.complete && !r.exhaustive_for_selected_records);
    let mut o = OffTargetOptions::default();
    o.max_index_entries = 1;
    let r = search_off_targets(&q, &refs, &o).unwrap();
    assert!(!r.complete);
    let mut o = OffTargetOptions::default();
    o.max_hits = 1;
    let r = search_off_targets(&q, &refs, &o).unwrap();
    assert!(!r.complete);
    assert_eq!(r.hits.len(), 1);
    let mut rich = m("ACGT", Topology::Linear);
    rich.original_genbank = Some("x".repeat(16 * 1024 * 1024));
    assert!(matches!(
        guide_revision(&rich, &GuideOptions::default()),
        Err(Error::TooLarge)
    ));
    let mut duplicate = refs.clone();
    duplicate.push(refs[0].clone());
    assert!(off_target_revision(&q, &duplicate, &OffTargetOptions::default()).is_err());
}
#[test]
fn source_options_profile_and_plausible_imported_evidence_need_fresh_replay() {
    let source = m(&format!("AAA{GUIDE}TGGCCC"), Topology::Linear);
    let o = GuideOptions::default();
    let r = design_guides(&source, &o).unwrap();
    replay_guides(&source, &o, &r).unwrap();
    let mut forged = r.clone();
    forged.candidates[0].nominal_top_cut = Some(0);
    forged.validate_bindings_for(&source, &o).unwrap();
    assert!(replay_guides(&source, &o, &forged).is_err());
    assert!(validate_candidate(&source, &o, &forged.candidates[0]).is_err());
    let mut old = r.clone();
    old.engine_version = "historical-crispr".into();
    old.validate().unwrap();
    assert!(old.validate_bindings_for(&source, &o).is_err());
    let mut fake = o.clone();
    fake.profile.pam = "NAG".into();
    assert!(design_guides(&source, &fake).is_err());
    let refs = vec![reference("ref", source.sequence.as_str(), Topology::Linear)];
    let query = GuideQuery {
        id: "guide".into(),
        sequence: GUIDE.into(),
    };
    let o = OffTargetOptions::default();
    let report = search_off_targets(&query, &refs, &o).unwrap();
    replay_off_targets(&query, &refs, &o, &report).unwrap();
    let mut forged = report.clone();
    forged.hits[0].nominal_top_cut = Some(0);
    forged.validate_bindings_for(&query, &refs, &o).unwrap();
    assert!(replay_off_targets(&query, &refs, &o, &forged).is_err());
    let mut value = serde_json::to_value(o).unwrap();
    value["assume_genome_complete"] = serde_json::json!(true);
    assert!(serde_json::from_value::<OffTargetOptions>(value).is_err());
}
#[test]
fn declared_five_prime_profiles_index_all_supported_hamming_thresholds() {
    let query = "ACGTCAGTACGA";
    let mut changed = query.as_bytes().to_vec();
    changed[1] = b'T';
    changed[5] = b'A';
    let text = format!(
        "GCTA{query}AAATG{}A{}",
        String::from_utf8(changed).unwrap(),
        rc(&format!("TA{query}"))
    );
    let profile = NucleaseProfile {
        id: "declared-fixture".into(),
        version: "literal-v1".into(),
        guide_length: 12,
        pam: "TR".into(),
        pam_side: PamSide::FivePrime,
        nominal_top_offset: 8,
        nominal_bottom_offset: 10,
        custom: true,
        sources: vec![],
    };
    let q = GuideQuery {
        id: "declared".into(),
        sequence: query.into(),
    };
    for topology in [Topology::Linear, Topology::Circular] {
        for shift in 0..text.len() {
            let text = format!("{}{}", &text[shift..], &text[..shift]);
            let refs = vec![reference("ref", &text, topology)];
            for mismatches in 0..=6 {
                let o = OffTargetOptions {
                    profile: profile.clone(),
                    max_mismatches: mismatches,
                    ..Default::default()
                };
                let r = search_off_targets(&q, &refs, &o).unwrap();
                assert!(r.complete, "{:?}", r.limits);
                let mut expected = vec![];
                let n = text.len();
                for reverse in [false, true] {
                    let s = if reverse { rc(&text) } else { text.clone() };
                    let count = if topology == Topology::Circular {
                        n
                    } else {
                        n - 14 + 1
                    };
                    for start in 0..count {
                        let b = |j: usize| s.as_bytes()[j % n];
                        if b(start) != b'T' || !b"AG".contains(&b(start + 1)) {
                            continue;
                        }
                        let d = query
                            .bytes()
                            .enumerate()
                            .filter(|(j, c)| b(start + 2 + j) != *c)
                            .count();
                        if d <= mismatches {
                            let positions = (0..12)
                                .map(|j| {
                                    let at = (start + 2 + j) % n;
                                    if reverse { n - 1 - at } else { at }
                                })
                                .collect::<Vec<_>>();
                            expected.push((reverse, positions, d));
                        }
                    }
                }
                let mut actual = r
                    .hits
                    .iter()
                    .map(|h| (h.reverse, h.guide_positions.clone(), h.min_mismatches))
                    .collect::<Vec<_>>();
                actual.sort();
                expected.sort();
                assert_eq!(actual, expected);
            }
        }
    }
}
