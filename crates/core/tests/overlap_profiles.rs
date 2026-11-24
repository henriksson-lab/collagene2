use collagene_core::{
    Alphabet, Feature, Location, Molecule, Sequence, Topology, overlap_profiles::*,
};
use std::collections::BTreeSet;
fn dna(name: &str, text: &str) -> Molecule {
    Molecule::dna(name, text, Topology::Linear).unwrap()
}
fn literal_rc(text: &str) -> String {
    text.chars()
        .rev()
        .map(|c| match c {
            'A' => 'T',
            'T' => 'A',
            'C' => 'G',
            'G' => 'C',
            _ => panic!("literal DNA only"),
        })
        .collect()
}
fn options(method: OverlapMethod) -> ProfileOptions {
    ProfileOptions {
        method,
        ..Default::default()
    }
}
fn simple(h: &str) -> Vec<Molecule> {
    vec![
        dna("left", &format!("TTGCA{h}")),
        dna("right", &format!("{h}GGTAC")),
    ]
}
#[test]
fn literal_fifteen_base_junction_has_checked_geometry_and_advisory_guidance() {
    let h = "ACGTCGATGCTAGCA";
    let m = simple(h);
    let r = check_fragments(&m, &[15], false, &options(OverlapMethod::HiFi)).unwrap();
    assert!(r.geometry_valid && r.complete);
    assert_eq!(r.recommendation_compliant, Some(true));
    assert!(!r.exhaustive_alternative_search);
    let j = &r.junctions[0];
    assert_eq!(j.homology, h);
    assert_eq!(j.downstream_homology, h);
    assert_eq!(j.gc_fraction, Some(8. / 15.));
    assert_eq!((j.recommended_min, j.recommended_max), (Some(15), Some(20)));
    assert_eq!(
        j.upstream_span,
        Some(ReferenceSpan {
            start: 5,
            end: 20,
            reverse: false
        })
    );
    assert_eq!(
        j.downstream_span,
        Some(ReferenceSpan {
            start: 0,
            end: 15,
            reverse: false
        })
    );
    assert_eq!(
        j.occurrences.iter().filter(|o| o.assigned_junction).count(),
        2
    );
    assert_eq!(
        r.source_revisions,
        m.iter().map(Molecule::revision).collect::<Vec<_>>()
    );
}
#[test]
fn reversing_inputs_preserves_original_identity_and_original_coordinate_spans() {
    let h = "ACGTCGATGCTAGCA";
    let base = simple(h);
    for a in [false, true] {
        for b in [false, true] {
            let flags = [a, b];
            let sources = base
                .iter()
                .zip(flags)
                .map(|(m, reverse)| {
                    let text = if reverse {
                        literal_rc(m.sequence.as_str())
                    } else {
                        m.sequence.as_str().to_owned()
                    };
                    dna(&m.name, &text)
                })
                .collect::<Vec<_>>();
            let r = check_oriented_fragments(
                &sources,
                &flags,
                &[15],
                false,
                &ProfileOptions::default(),
            )
            .unwrap();
            assert!(r.geometry_valid);
            assert_eq!(r.orientations, flags);
            assert_eq!(
                r.source_revisions,
                sources.iter().map(Molecule::revision).collect::<Vec<_>>()
            );
            let j = &r.junctions[0];
            assert_eq!(j.homology, h);
            assert_eq!(
                j.upstream_span,
                Some(if a {
                    ReferenceSpan {
                        start: 0,
                        end: 15,
                        reverse: true,
                    }
                } else {
                    ReferenceSpan {
                        start: 5,
                        end: 20,
                        reverse: false,
                    }
                })
            );
            assert_eq!(
                j.downstream_span,
                Some(if b {
                    ReferenceSpan {
                        start: 5,
                        end: 20,
                        reverse: true,
                    }
                } else {
                    ReferenceSpan {
                        start: 0,
                        end: 15,
                        reverse: false,
                    }
                })
            );
        }
    }
    // Reversing a declared source is a choice, not an automatically chosen best match.
    assert!(
        !check_oriented_fragments(
            &base,
            &[false, true],
            &[15],
            false,
            &ProfileOptions::default()
        )
        .unwrap()
        .geometry_valid
    );
}
#[test]
fn method_length_guidance_never_redefines_exact_sequence_compatibility() {
    let h = "ACGTACGTACGTACGTACGTACGTA";
    assert_eq!(h.len(), 25);
    let m = simple(h);
    let hi = check_fragments(&m, &[25], false, &options(OverlapMethod::HiFi)).unwrap();
    let gib = check_fragments(&m, &[25], false, &options(OverlapMethod::Gibson)).unwrap();
    let inf = check_fragments(&m, &[25], false, &options(OverlapMethod::InFusion)).unwrap();
    assert!(hi.geometry_valid && gib.geometry_valid && inf.geometry_valid);
    assert_eq!(hi.recommendation_compliant, Some(false));
    assert_eq!(gib.recommendation_compliant, Some(true));
    assert_eq!(inf.recommendation_compliant, Some(true));
    assert_eq!(inf.junctions[0].recommended_max, None);
    assert!(hi.findings.iter().any(|s| s.contains("advisory")));
    let short = check_fragments(&simple("ACGT"), &[4], false, &ProfileOptions::default()).unwrap();
    assert!(short.geometry_valid);
    assert_eq!(short.recommendation_compliant, Some(false));
}
#[test]
fn circular_three_fragment_closure_is_an_actual_last_to_first_junction() {
    let h = ["ACGTCGATGCTAGCA", "TTCAGGATCCGATGC", "GATCTAGCGTACGTT"];
    let m = vec![
        dna("a", &format!("{}AAAA{}", h[2], h[0])),
        dna("b", &format!("{}CCCC{}", h[0], h[1])),
        dna("c", &format!("{}TTTT{}", h[1], h[2])),
    ];
    let r = check_fragments(&m, &[15, 15, 15], true, &ProfileOptions::default()).unwrap();
    assert!(r.geometry_valid && r.complete);
    assert_eq!(
        r.junctions
            .iter()
            .map(|j| j.homology.as_str())
            .collect::<Vec<_>>(),
        h
    );
    let j = &r.junctions[2];
    assert_eq!((j.upstream, j.downstream), (2, 0));
    assert_eq!(
        j.upstream_span,
        Some(ReferenceSpan {
            start: 19,
            end: 34,
            reverse: false
        })
    );
    assert_eq!(
        j.downstream_span,
        Some(ReferenceSpan {
            start: 0,
            end: 15,
            reverse: false
        })
    );
    assert_eq!(r.source_bases - r.overlaps.iter().sum::<usize>(), 57);
    let mut wrong = m.clone();
    wrong[0].sequence = Sequence::dna(&format!("{}AAAA{}", h[1], h[0])).unwrap();
    assert!(
        !check_fragments(&wrong, &[15, 15, 15], true, &ProfileOptions::default())
            .unwrap()
            .geometry_valid
    );
}
#[test]
fn four_fragment_ranges_and_unsupported_table_counts_remain_explicit() {
    let h = "ACGTCGATGCTAGCATCGAT";
    assert_eq!(h.len(), 20);
    let m = (0..4)
        .map(|i| dna(&format!("f{i}"), &format!("{h}AA{h}")))
        .collect::<Vec<_>>();
    for (method, min, max) in [
        (OverlapMethod::HiFi, 20, Some(30)),
        (OverlapMethod::Gibson, 20, Some(80)),
        (OverlapMethod::InFusion, 20, None),
    ] {
        let r = check_fragments(&m, &[20; 4], true, &options(method)).unwrap();
        assert!(r.geometry_valid);
        assert_eq!(r.junctions[0].recommended_min, Some(min));
        assert_eq!(r.junctions[0].recommended_max, max);
        assert!(r.findings.iter().any(|f| f.contains("additional concrete")));
    }
    let seven = vec![dna("copy", &format!("{h}AA{h}")); 7];
    let r = check_fragments(&seven, &[20; 7], true, &ProfileOptions::default()).unwrap();
    assert_eq!(r.recommendation_compliant, None);
    assert!(r.warnings.iter().any(|w| w.contains("outside")));
    let single = dna("one", "ACGT");
    let r = check_fragments(&[single], &[], false, &ProfileOptions::default()).unwrap();
    assert!(r.geometry_valid && r.complete);
    assert_eq!(r.recommendation_compliant, None);
    assert!(r.junctions.is_empty());
}
#[test]
fn mismatch_zero_and_exhaustive_overlaps_do_not_become_valid_homology() {
    let mut m = simple("ACGTCGATGCTAGCA");
    m[1].sequence = Sequence::dna("TCGTCGATGCTAGCAGGTAC").unwrap();
    let o = ProfileOptions {
        scan_repeated_overlaps: false,
        ..Default::default()
    };
    let r = check_fragments(&m, &[15], false, &o).unwrap();
    assert!(!r.geometry_valid);
    assert_ne!(r.junctions[0].homology, r.junctions[0].downstream_homology);
    assert!(r.complete);
    assert!(!r.junctions[0].repeat_search_complete);
    for k in [0, 20, 21] {
        let r = check_fragments(&m, &[k], false, &o).unwrap();
        assert!(!r.geometry_valid);
        assert_eq!(r.junctions[0].gc_fraction, None);
        assert!(r.junctions[0].homology.is_empty());
    }
    let overlap_arms = vec![dna("x", "AAAAAA"), dna("y", "AAAAAA")];
    assert!(matches!(
        check_fragments(&overlap_arms, &[4, 4], true, &o),
        Err(collagene_core::Error::Unsupported(_))
    ));
}
#[test]
fn repeat_scans_count_extra_contexts_without_claiming_alternative_products() {
    let m = vec![dna("left", "CAAAAAA"), dna("right", "AAAAGGG")];
    let r = check_fragments(&m, &[4], false, &ProfileOptions::default()).unwrap();
    let j = &r.junctions[0];
    assert!(r.complete && j.repeat_search_complete);
    assert_eq!(j.occurrences.len(), 4);
    assert_eq!(
        j.occurrences.iter().filter(|o| o.assigned_junction).count(),
        2
    );
    assert!(r.findings.iter().any(|f| f.contains("2 additional")));
    assert!(!r.exhaustive_alternative_search);
}
#[test]
fn exhaustive_small_reference_words_have_an_independent_both_strand_occurrence_oracle() {
    let alphabet = b"ACGT";
    for code in 0..256 {
        let mut number = code;
        let word = (0..4)
            .map(|_| {
                let c = alphabet[number % 4] as char;
                number /= 4;
                c
            })
            .collect::<String>();
        let h = &word[2..];
        for flags in [[false, false], [false, true], [true, false], [true, true]] {
            let right = format!("{h}G");
            let raw = [
                if flags[0] {
                    literal_rc(&word)
                } else {
                    word.clone()
                },
                if flags[1] {
                    literal_rc(&right)
                } else {
                    right.clone()
                },
            ];
            let sources = vec![dna("a", &raw[0]), dna("b", &raw[1])];
            let report =
                check_oriented_fragments(&sources, &flags, &[2], false, &ProfileOptions::default())
                    .unwrap();
            assert!(report.geometry_valid && report.complete);
            let mut expected = BTreeSet::new();
            for (f, text) in raw.iter().enumerate() {
                let backwards = literal_rc(h);
                for (query, reverse) in [(h, false), (backwards.as_str(), true)] {
                    if reverse && backwards == h {
                        continue;
                    }
                    for start in 0..=text.len() - 2 {
                        if &text[start..start + 2] == query {
                            expected.insert((f, start, start + 2, reverse));
                        }
                    }
                }
            }
            let actual = report.junctions[0]
                .occurrences
                .iter()
                .map(|o| (o.fragment, o.start, o.end, o.reverse))
                .collect::<BTreeSet<_>>();
            assert_eq!(actual, expected, "{word} {flags:?}");
        }
    }
}
#[test]
fn budgets_are_errors_for_unestablished_geometry_and_incomplete_for_repeat_truncation() {
    let m = simple("ACGTCGATGCTAGCA");
    let o = ProfileOptions {
        max_comparisons: 14,
        ..Default::default()
    };
    assert!(check_fragments(&m, &[15], false, &o).is_err());
    let o = ProfileOptions {
        max_comparisons: 15,
        ..Default::default()
    };
    let r = check_fragments(&m, &[15], false, &o).unwrap();
    assert_eq!(r.comparisons, 15);
    assert!(!r.complete);
    assert!(!r.junctions[0].repeat_search_complete);
    assert!(r.junctions[0].occurrences.is_empty());
    let o = ProfileOptions {
        max_occurrences: 1,
        ..Default::default()
    };
    let r = check_fragments(&m, &[15], false, &o).unwrap();
    assert!(!r.complete);
    assert_eq!(r.junctions[0].occurrences.len(), 1);
    assert!(r.warnings.iter().any(|w| w.contains("stopped")));
}
#[test]
fn input_profiles_are_strict_and_source_work_is_bounded_before_hashing() {
    let m = simple("ACGT");
    assert!(check_fragments(&[], &[], false, &ProfileOptions::default()).is_err());
    assert!(
        check_oriented_fragments(&m, &[false], &[4], false, &ProfileOptions::default()).is_err()
    );
    assert!(check_fragments(&m, &[], false, &ProfileOptions::default()).is_err());
    assert!(check_fragments(&m, &[10001], false, &ProfileOptions::default()).is_err());
    for sequence in [
        Sequence::dna("ACGN").unwrap(),
        Sequence::new(Alphabet::Rna, "ACGU").unwrap(),
        Sequence::new(Alphabet::Protein, "ACGT").unwrap(),
    ] {
        let mut bad = m.clone();
        bad[0].sequence = sequence;
        assert!(check_fragments(&bad, &[2], false, &ProfileOptions::default()).is_err());
    }
    let mut circle = m.clone();
    circle[0].topology = Topology::Circular;
    assert!(check_fragments(&circle, &[4], false, &ProfileOptions::default()).is_err());
    let mut deep = Location::range(0, 1);
    for _ in 0..66 {
        deep = Location::Complement(Box::new(deep));
    }
    let mut bad = m.clone();
    bad[0].features.push(Feature {
        kind: "test".into(),
        location: deep,
        qualifiers: vec![],
        needs_review: false,
    });
    assert!(matches!(
        check_fragments(&bad, &[4], false, &ProfileOptions::default()),
        Err(collagene_core::Error::TooLarge)
    ));
    let mut huge = m;
    huge[0]
        .metadata
        .insert("giant".into(), "x".repeat(33 * 1024 * 1024));
    assert!(matches!(
        check_fragments(&huge, &[4], false, &ProfileOptions::default()),
        Err(collagene_core::Error::TooLarge)
    ));
}
#[test]
fn serialized_report_integrity_is_not_a_new_scientific_authority() {
    let r = check_fragments(
        &simple("ACGTCGATGCTAGCA"),
        &[15],
        false,
        &ProfileOptions::default(),
    )
    .unwrap();
    let bytes = serde_json::to_vec(&r).unwrap();
    let loaded: ProfileReport = serde_json::from_slice(&bytes).unwrap();
    loaded.validate().unwrap();
    assert_eq!(r, loaded);
    for mutation in 0..8 {
        let mut bad = r.clone();
        match mutation {
            0 => bad.geometry_valid = false,
            1 => bad.junctions[0].exact_match = false,
            2 => bad.junctions[0].gc_fraction = Some(f64::NAN),
            3 => bad.engine_version = "other".into(),
            4 => bad.junctions[0].occurrences[0].assigned_junction = false,
            5 => bad.comparisons = 0,
            6 => bad.orientations[0] = true,
            _ => bad.exhaustive_alternative_search = true,
        }
        assert!(bad.validate().is_err(), "mutation {mutation}");
    }
    let mut value = serde_json::to_value(&r).unwrap();
    value["options"]["unknown"] = "ignored?".into();
    assert!(serde_json::from_value::<ProfileReport>(value).is_err());
}
