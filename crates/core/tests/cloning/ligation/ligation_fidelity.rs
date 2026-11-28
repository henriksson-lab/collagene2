//! Independent literal count fixtures; no published experimental data bundled.
use collagene_core::cloning::ligation::fidelity::{
    self, Cell, Conditions, Dataset, Estimate, Prepared, Provenance, UnknownReason,
};
fn dataset(k: usize, counts: &[(&str, &str, u64)]) -> Dataset {
    let mut data = Dataset {
        schema_version: 1,
        overhang_length: k,
        content_revision: String::new(),
        conditions: Conditions {
            ligase: "Synthetic fixture ligase; not empirical".into(),
            buffer: "Declared fixture buffer".into(),
            temperature_protocol: "Constant fixture condition".into(),
            incubation: "Fixture observation period".into(),
            restriction_enzymes: vec![],
            substrate_context: "Equally represented directed synthetic observations".into(),
        },
        provenance: Provenance {
            source_url: "fixture://independent-literal-counts".into(),
            citation: "Independent test arithmetic, no experimental dataset".into(),
            license: "MIT synthetic fixture".into(),
            attribution: "Collagene test authors".into(),
            source_blake3: "a".repeat(64),
            normalization: "Literal directed counts; no transpose or pseudocount".into(),
        },
        cells: counts
            .iter()
            .map(|(a, b, count)| Cell {
                end_5to3: (*a).into(),
                partner_5to3: (*b).into(),
                count: *count,
            })
            .collect(),
    };
    data.content_revision = fidelity::content_revision(&data).unwrap();
    data
}
fn ends(v: &[&str]) -> Vec<String> {
    v.iter().map(|s| (*s).into()).collect()
}
fn competing() -> Dataset {
    dataset(
        3,
        &[
            ("AAA", "AAA", 10),
            ("AAA", "TTT", 90),
            ("AAA", "CCC", 0),
            ("AAA", "GGG", 0),
            ("TTT", "AAA", 30),
            ("TTT", "TTT", 10),
            ("TTT", "CCC", 0),
            ("TTT", "GGG", 0),
            ("CCC", "AAA", 10),
            ("CCC", "TTT", 0),
            ("CCC", "CCC", 10),
            ("CCC", "GGG", 40),
            ("GGG", "AAA", 20),
            ("GGG", "TTT", 0),
            ("GGG", "CCC", 20),
            ("GGG", "GGG", 0),
        ],
    )
}
#[test]
fn actual_directed_rows_condition_on_present_physical_partners_once_per_rc_pair() {
    let data = competing();
    let requested = ends(&["AAA", "TTT", "CCC", "AAA"]);
    let report = fidelity::score(&data, &requested).unwrap();
    assert_eq!(report.physical_ends, ends(&["AAA", "CCC", "GGG", "TTT"]));
    assert_eq!(report.pairs.len(), 2);
    assert_eq!(report.pairs[0].occurrences, vec![0, 1, 3]);
    let a = &report.pairs[0];
    assert_eq!(a.canonical_end, "AAA");
    assert_eq!(a.reverse_complement, "TTT");
    assert_eq!(a.orientations[0].total_count, Some(100));
    assert_eq!(a.orientations[1].total_count, Some(40));
    assert_eq!(
        a.estimate,
        Estimate::Known {
            correct_count: 120,
            total_count: 140,
            fidelity: 6.0 / 7.0
        }
    );
    assert_eq!(
        report.pairs[1].estimate,
        Estimate::Known {
            correct_count: 60,
            total_count: 100,
            fidelity: 0.6
        }
    );
    assert!((report.overall_fidelity.unwrap() - 18.0 / 35.0).abs() < 1e-15);
    assert!((report.log_overall_fidelity.unwrap() - (18.0f64 / 35.0).ln()).abs() < 1e-15);
    report
        .validate_binding(&fidelity::revision(&data, &requested).unwrap())
        .unwrap();
    let reversed = Prepared::new(&data)
        .unwrap()
        .score(&ends(&["GGG", "TTT"]))
        .unwrap();
    assert_eq!(reversed.overall_fidelity, report.overall_fidelity);
    assert_ne!(
        reversed.context_revision, report.context_revision,
        "request identity retains order and duplicate declarations even when score deduplicates pairs"
    );
    // Removing a competing physical pair changes the conditional scope, not the dataset.
    let one = fidelity::score(&data, &ends(&["CCC"])).unwrap();
    assert_eq!(
        one.pairs[0].estimate,
        Estimate::Known {
            correct_count: 60,
            total_count: 70,
            fidelity: 6.0 / 7.0
        }
    );
}
#[test]
fn partial_counts_missing_diagonals_or_reverse_rows_are_unknown_never_zero_imputed() {
    let base = dataset(
        3,
        &[
            ("AAA", "AAA", 10),
            ("AAA", "TTT", 90),
            ("TTT", "AAA", 30),
            ("TTT", "TTT", 10),
        ],
    );
    assert_eq!(
        fidelity::score(&base, &ends(&["AAA"]))
            .unwrap()
            .overall_fidelity,
        Some(6.0 / 7.0)
    );
    for removed in [("AAA", "AAA"), ("TTT", "AAA")] {
        let mut missing = base.clone();
        missing
            .cells
            .retain(|c| (c.end_5to3.as_str(), c.partner_5to3.as_str()) != removed);
        missing.content_revision = fidelity::content_revision(&missing).unwrap();
        let report = fidelity::score(&missing, &ends(&["AAA"])).unwrap();
        assert_eq!(report.overall_fidelity, None);
        assert_eq!(report.log_overall_fidelity, None);
        assert_eq!(
            report.pairs[0].estimate,
            Estimate::Unknown {
                reason: UnknownReason::MissingCounts
            }
        );
        let row = report.pairs[0]
            .orientations
            .iter()
            .find(|r| r.end_5to3 == removed.0)
            .unwrap();
        assert_eq!(row.total_count, None);
        assert_eq!(row.missing_partners, vec![removed.1]);
    }
    // All row cells are required even when the omitted known fixture count was0.
    let mut missing_zero = competing();
    missing_zero
        .cells
        .retain(|c| !(c.end_5to3 == "AAA" && c.partner_5to3 == "CCC"));
    missing_zero.content_revision = fidelity::content_revision(&missing_zero).unwrap();
    let r = fidelity::score(&missing_zero, &ends(&["AAA", "CCC"])).unwrap();
    assert_eq!(
        r.pairs[0].estimate,
        Estimate::Unknown {
            reason: UnknownReason::MissingCounts
        }
    );
    assert!(matches!(r.pairs[1].estimate, Estimate::Known { .. }));
    assert_eq!(r.overall_fidelity, None);
}
#[test]
fn observed_zero_no_observations_palindromes_and_u64_extremes_are_distinct() {
    let zero = dataset(
        3,
        &[
            ("AAA", "AAA", 1),
            ("AAA", "TTT", 0),
            ("TTT", "AAA", 0),
            ("TTT", "TTT", 1),
        ],
    );
    let report = fidelity::score(&zero, &ends(&["AAA"])).unwrap();
    assert_eq!(report.overall_fidelity, Some(0.0));
    assert_eq!(report.log_overall_fidelity, None);
    assert_eq!(
        report.pairs[0].estimate,
        Estimate::Known {
            correct_count: 0,
            total_count: 2,
            fidelity: 0.0
        }
    );
    let empty = dataset(
        3,
        &[
            ("AAA", "AAA", 0),
            ("AAA", "TTT", 0),
            ("TTT", "AAA", 0),
            ("TTT", "TTT", 0),
        ],
    );
    assert_eq!(
        fidelity::score(&empty, &ends(&["AAA"])).unwrap().pairs[0].estimate,
        Estimate::Unknown {
            reason: UnknownReason::ZeroTotal
        }
    );
    let palindrome = dataset(4, &[("ATAT", "ATAT", 100)]);
    let p = fidelity::score(&palindrome, &ends(&["ATAT", "ATAT"])).unwrap();
    assert_eq!(p.physical_ends, vec!["ATAT"]);
    assert_eq!(p.pairs[0].orientations.len(), 1);
    assert_eq!(
        p.pairs[0].estimate,
        Estimate::Unknown {
            reason: UnknownReason::PalindromicEnd
        }
    );
    assert_eq!(p.overall_fidelity, None);
    let max = dataset(
        4,
        &[
            ("AAAA", "AAAA", u64::MAX),
            ("AAAA", "TTTT", u64::MAX),
            ("TTTT", "AAAA", u64::MAX),
            ("TTTT", "TTTT", u64::MAX),
        ],
    );
    assert_eq!(
        fidelity::score(&max, &ends(&["AAAA"])).unwrap().pairs[0].estimate,
        Estimate::Known {
            correct_count: u128::from(u64::MAX) * 2,
            total_count: u128::from(u64::MAX) * 4,
            fidelity: 0.5
        }
    );
}
#[test]
fn dataset_revision_conditions_provenance_and_full_matrix_bounds_are_explicit() {
    let mut data = competing();
    let original = data.content_revision.clone();
    data.cells.reverse();
    assert_eq!(fidelity::content_revision(&data).unwrap(), original);
    data.validate().unwrap();
    assert_eq!(
        fidelity::decode(&serde_json::to_vec(&data).unwrap()).unwrap(),
        data
    );
    let mut changed = data.clone();
    changed.conditions.temperature_protocol = "Another declared condition".into();
    assert!(changed.validate().is_err());
    changed.content_revision = fidelity::content_revision(&changed).unwrap();
    assert_ne!(
        fidelity::revision(&changed, &ends(&["AAA"])).unwrap(),
        fidelity::revision(&data, &ends(&["AAA"])).unwrap()
    );
    let mut duplicate = data.clone();
    duplicate.cells.push(duplicate.cells[0].clone());
    assert!(fidelity::content_revision(&duplicate).is_err());
    let mut wrong = data.clone();
    wrong.provenance.source_blake3 = "Z".repeat(64);
    assert!(fidelity::content_revision(&wrong).is_err());
    wrong = data.clone();
    wrong.conditions.substrate_context.clear();
    assert!(fidelity::content_revision(&wrong).is_err());
    wrong = data.clone();
    wrong.cells[0].end_5to3 = "aaa".into();
    assert!(fidelity::content_revision(&wrong).is_err());
    assert!(fidelity::score(&data, &ends(&["AAN"])).is_err());
    assert!(fidelity::score(&data, &vec!["AAA".into(); 257]).is_err());
    assert!(fidelity::score(&data, &[]).is_err());
    // Every3/4nt cell is accepted under explicit bounds, without genotype search.
    for k in [3usize, 4] {
        let strings = (0..4usize.pow(k as u32))
            .map(|mut i| {
                let mut s = vec![b'A'; k];
                for b in s.iter_mut().rev() {
                    *b = b"ACGT"[i % 4];
                    i /= 4;
                }
                String::from_utf8(s).unwrap()
            })
            .collect::<Vec<_>>();
        let mut full = data.clone();
        full.overhang_length = k;
        full.cells = strings
            .iter()
            .flat_map(|a| {
                strings.iter().map(move |b| Cell {
                    end_5to3: a.clone(),
                    partner_5to3: b.clone(),
                    count: 1,
                })
            })
            .collect();
        full.content_revision = fidelity::content_revision(&full).unwrap();
        let e = if k == 3 {
            ends(&["AAA", "CCC"])
        } else {
            ends(&["AAAA", "CCCC"])
        };
        assert_eq!(
            fidelity::score(&full, &e).unwrap().overall_fidelity,
            Some(1.0 / 16.0)
        );
        if k == 4 {
            assert_eq!(full.cells.len(), fidelity::MAX_CELLS);
            full.cells.push(full.cells[0].clone());
            assert!(fidelity::content_revision(&full).is_err());
        }
    }
}
#[test]
fn score_receipt_tampering_and_binding_changes_do_not_validate() {
    let data = competing();
    let report = fidelity::score(&data, &ends(&["AAA", "CCC"])).unwrap();
    let token = report.context_revision.clone();
    let mut changed = report.clone();
    changed.pairs[0].orientations[0].total_count = Some(1);
    assert!(changed.validate_binding(&token).is_err());
    changed = report.clone();
    changed.pairs[0].occurrences = vec![1];
    assert!(changed.validate_binding(&token).is_err());
    changed = report.clone();
    changed.conditions.ligase = "Relabelled ligase".into();
    assert!(changed.validate_binding(&token).is_err());
    changed = report.clone();
    changed.overall_fidelity = Some(1.0);
    assert!(changed.validate_binding(&token).is_err());
    changed = report.clone();
    changed.pairs[0].orientations[0].missing_partners = vec!["AAA".into()];
    assert!(changed.validate_binding(&token).is_err());
    assert!(report.validate_binding(&"0".repeat(64)).is_err());
    assert_eq!(
        serde_json::from_slice::<fidelity::Report>(&serde_json::to_vec(&report).unwrap()).unwrap(),
        report
    );
}

#[test]
fn tiny_positive_products_preserve_log_ranking_instead_of_claiming_observed_zero() {
    let strings = (0..256usize)
        .map(|n| {
            [n / 64, (n / 16) % 4, (n / 4) % 4, n % 4]
                .into_iter()
                .map(|i| ["A", "C", "G", "T"][i])
                .collect::<String>()
        })
        .collect::<Vec<_>>();
    let complement = |s: &str| {
        s.chars()
            .rev()
            .map(|c| match c {
                'A' => 'T',
                'C' => 'G',
                'G' => 'C',
                _ => 'A',
            })
            .collect::<String>()
    };
    let requested = strings
        .iter()
        .filter(|s| s.as_str() < complement(s).as_str())
        .cloned()
        .collect::<Vec<_>>();
    assert_eq!(requested.len(), 120);
    let mut data = competing();
    data.overhang_length = 4;
    data.cells = strings
        .iter()
        .flat_map(|a| {
            strings.iter().map(move |b| Cell {
                end_5to3: a.clone(),
                partner_5to3: b.clone(),
                count: if a == b {
                    u64::MAX
                } else if *b == complement(a) {
                    1
                } else {
                    0
                },
            })
        })
        .collect();
    data.content_revision = fidelity::content_revision(&data).unwrap();
    let report = fidelity::score(&data, &requested).unwrap();
    assert_eq!(
        report.overall_fidelity,
        Some(0.0),
        "floating product underflows"
    );
    assert!(report.pairs.iter().all(|p| matches!(
        p.estimate,
        Estimate::Known {
            correct_count: 2,
            ..
        }
    )));
    let wanted = 120.0 * (1.0 / (u64::MAX as f64 + 1.0)).ln();
    assert!((report.log_overall_fidelity.unwrap() - wanted).abs() < 1e-9);
    report.validate_binding(&report.context_revision).unwrap();
}

#[test]
fn wide_count_json_is_exact_decimal_while_postcard_keeps_native_u128() {
    let data = dataset(
        4,
        &[
            ("AAAA", "AAAA", u64::MAX),
            ("AAAA", "TTTT", u64::MAX),
            ("TTTT", "AAAA", u64::MAX),
            ("TTTT", "TTTT", u64::MAX),
        ],
    );
    let report = fidelity::score(&data, &ends(&["AAAA"])).unwrap();
    let value = serde_json::to_value(&report).unwrap();
    assert_eq!(
        value["pairs"][0]["orientations"][0]["total_count"],
        "36893488147419103230"
    );
    assert_eq!(
        value["pairs"][0]["estimate"]["Known"]["correct_count"],
        "36893488147419103230"
    );
    assert_eq!(
        value["pairs"][0]["estimate"]["Known"]["total_count"],
        "73786976294838206460"
    );
    assert_eq!(
        serde_json::from_value::<fidelity::Report>(value.clone()).unwrap(),
        report
    );
    assert_eq!(
        serde_json::from_slice::<fidelity::Report>(&serde_json::to_vec(&report).unwrap()).unwrap(),
        report
    );
    let binary = postcard::to_allocvec(&report).unwrap();
    assert_eq!(
        postcard::from_bytes::<fidelity::Report>(&binary).unwrap(),
        report
    );
    // Standalone field shapes demonstrate the binary representation remains
    // rawu128, independent of the full report's positional layout.
    let known = Estimate::Known {
        correct_count: u128::from(u64::MAX) * 2,
        total_count: u128::from(u64::MAX) * 4,
        fidelity: 0.5,
    };
    let binary = postcard::to_allocvec(&known).unwrap();
    let decoded: (u32, u128, u128, f64) = postcard::from_bytes(&binary).unwrap();
    assert_eq!(
        decoded,
        (0, u128::from(u64::MAX) * 2, u128::from(u64::MAX) * 4, 0.5)
    );
    for invalid in [
        serde_json::json!(1),
        serde_json::json!(""),
        serde_json::json!("+1"),
        serde_json::json!("-1"),
        serde_json::json!("01"),
        serde_json::json!("1.0"),
        serde_json::json!(" 1"),
        serde_json::json!("1e2"),
        serde_json::json!("340282366920938463463374607431768211456"),
        serde_json::json!("9".repeat(40)),
    ] {
        let mut changed = value.clone();
        changed["pairs"][0]["estimate"]["Known"]["total_count"] = invalid.clone();
        assert!(serde_json::from_value::<fidelity::Report>(changed).is_err());
        let mut changed = value.clone();
        changed["pairs"][0]["orientations"][0]["total_count"] = invalid;
        assert!(serde_json::from_value::<fidelity::Report>(changed).is_err());
    }
    let partial = dataset(3, &[("AAA", "TTT", 1)]);
    let partial = fidelity::score(&partial, &ends(&["AAA"])).unwrap();
    let json = serde_json::to_value(&partial).unwrap();
    assert!(json["pairs"][0]["orientations"][0]["total_count"].is_null());
    assert_eq!(
        serde_json::from_value::<fidelity::Report>(json).unwrap(),
        partial
    );
    assert_eq!(
        postcard::from_bytes::<fidelity::Report>(&postcard::to_allocvec(&partial).unwrap())
            .unwrap(),
        partial
    );
}

#[test]
fn literal_competitors_change_denominators_without_inventing_complements_or_intended_joins() {
    let data = dataset(
        3,
        &[
            ("AAA", "AAA", 10),
            ("AAA", "TTT", 90),
            ("AAA", "CCC", 20),
            ("AAA", "GGG", 500),
            ("TTT", "AAA", 30),
            ("TTT", "TTT", 10),
            ("TTT", "CCC", 10),
            ("TTT", "GGG", 300),
        ],
    );
    let intended = ends(&["AAA"]);
    let extras = ends(&["CCC"]);
    let scored = fidelity::score_competition(&data, &intended, &extras).unwrap();
    assert_eq!(scored.requested_ends, intended);
    assert_eq!(scored.physical_competitors, extras);
    assert_eq!(scored.physical_ends, ends(&["AAA", "CCC", "TTT"]));
    assert_eq!(
        scored.pairs.len(),
        1,
        "free CCC does not become an intended CCC/GGG join"
    );
    assert_eq!(
        scored.pairs[0].estimate,
        Estimate::Known {
            correct_count: 120,
            total_count: 170,
            fidelity: 12.0 / 17.0
        }
    );
    assert_eq!(scored.overall_fidelity, Some(12.0 / 17.0));
    assert_eq!(
        fidelity::score(&data, &intended).unwrap().overall_fidelity,
        Some(6.0 / 7.0)
    );
    let old_incorrect_scope = fidelity::score(&data, &ends(&["AAA", "CCC"])).unwrap();
    assert_eq!(
        old_incorrect_scope.overall_fidelity, None,
        "fabricating paired CCC/GGG would require unrelated missing rows"
    );
    scored
        .validate_binding(&fidelity::competition_revision(&data, &intended, &extras).unwrap())
        .unwrap();
    assert_eq!(
        serde_json::from_value::<fidelity::CompetitionReport>(
            serde_json::to_value(&scored).unwrap()
        )
        .unwrap(),
        scored
    );
    assert_eq!(
        postcard::from_bytes::<fidelity::CompetitionReport>(
            &postcard::to_allocvec(&scored).unwrap()
        )
        .unwrap(),
        scored
    );
    let mut absent = data.clone();
    absent
        .cells
        .retain(|c| !(c.end_5to3 == "AAA" && c.partner_5to3 == "CCC"));
    absent.content_revision = fidelity::content_revision(&absent).unwrap();
    assert_eq!(
        fidelity::score_competition(&absent, &intended, &extras)
            .unwrap()
            .pairs[0]
            .estimate,
        Estimate::Unknown {
            reason: UnknownReason::MissingCounts
        }
    );
    let duplicated = fidelity::score_competition(&data, &intended, &ends(&["CCC", "CCC"])).unwrap();
    assert_eq!(duplicated.overall_fidelity, scored.overall_fidelity);
    assert_ne!(duplicated.context_revision, scored.context_revision);
    let mut forged = scored.clone();
    forged.physical_ends.push("GGG".into());
    assert!(forged.validate_binding(&scored.context_revision).is_err());
    let mut forged = scored.clone();
    forged.physical_competitors = ends(&["GGG"]);
    assert!(forged.validate_binding(&scored.context_revision).is_err());
}
#[test]
fn unpaired_palindrome_is_only_a_competitor_and_no_intended_join_has_no_overall_estimate() {
    let data = dataset(
        4,
        &[
            ("AAAA", "AAAA", 10),
            ("AAAA", "TTTT", 90),
            ("AAAA", "ATAT", 5),
            ("TTTT", "AAAA", 30),
            ("TTTT", "TTTT", 10),
            ("TTTT", "ATAT", 15),
        ],
    );
    let scored = fidelity::score_competition(&data, &ends(&["AAAA"]), &ends(&["ATAT"])).unwrap();
    assert_eq!(scored.pairs.len(), 1);
    assert_eq!(scored.overall_fidelity, Some(0.75));
    assert!(matches!(
        scored.pairs[0].estimate,
        Estimate::Known {
            correct_count: 120,
            total_count: 160,
            ..
        }
    ));
    for extras in [ends(&["ATAT"]), vec![]] {
        let empty = fidelity::score_competition(&data, &[], &extras).unwrap();
        assert_eq!(empty.overall_fidelity, None);
        assert_eq!(empty.log_overall_fidelity, None);
        assert!(empty.pairs.is_empty());
        assert_eq!(empty.physical_ends, extras);
        empty
            .validate_binding(&fidelity::competition_revision(&data, &[], &extras).unwrap())
            .unwrap();
    }
    assert!(
        fidelity::score(&data, &[]).is_err(),
        "old API nonempty contract is preserved"
    );
}
