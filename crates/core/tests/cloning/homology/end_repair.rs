use collagene_core::chemistry::*;
use collagene_model::{Alphabet, Error, Molecule, Sequence, Topology};
fn duplex(top: &str, bottom: &str, offset: i64) -> Duplex {
    let d = Duplex {
        top: Sequence::dna(top).unwrap(),
        bottom: Sequence::dna(bottom).unwrap(),
        bottom_start: offset,
        topology: Topology::Linear,
        left: Terminus::CLEAVED,
        right: Terminus::CLEAVED,
        nicks: vec![],
    };
    d.validate().unwrap();
    d
}
fn manual_complement(base: u8) -> u8 {
    match base {
        b'A' => b'T',
        b'T' => b'A',
        b'C' => b'G',
        b'G' => b'C',
        _ => panic!("oracle concrete bases only"),
    }
}
fn manual_rc(bases: &str) -> String {
    bases
        .bytes()
        .rev()
        .map(|b| manual_complement(b) as char)
        .collect()
}
fn check_ancestry(report: &EndRepairReport) {
    let mut top = vec![0usize; report.duplex.top.len()];
    let mut bottom = vec![0usize; report.duplex.bottom.len()];
    for mapping in &report.mappings {
        let from = match mapping.source_strand {
            Strand::Top => report.source.top.as_str(),
            Strand::Bottom => report.source.bottom.as_str(),
        };
        let to = match mapping.product_strand {
            Strand::Top => report.duplex.top.as_str(),
            Strand::Bottom => report.duplex.bottom.as_str(),
        };
        let source = &from[mapping.source_start..mapping.source_end];
        let expected = if mapping.reverse {
            manual_rc(source)
        } else {
            source.into()
        };
        assert_eq!(expected, &to[mapping.product_start..mapping.product_end]);
        assert_eq!(
            mapping.synthesized,
            mapping.source_strand != mapping.product_strand
        );
        let coverage = match mapping.product_strand {
            Strand::Top => &mut top,
            Strand::Bottom => &mut bottom,
        };
        for c in &mut coverage[mapping.product_start..mapping.product_end] {
            *c += 1;
        }
    }
    assert!(top.iter().chain(&bottom).all(|n| *n == 1));
}
#[test]
fn literal_five_prime_fill_extends_only_recessed_three_prime_strands() {
    // Full reference union is AATGCCGT; neither original strand spans that union.
    let source = duplex("AATGCC", "ACGGCA", 2);
    let report = repair_ends(&source, EndRepairMethod::FillFivePrimeOverhangs).unwrap();
    assert_eq!(report.duplex.top.as_str(), "AATGCCGT");
    assert_eq!(report.duplex.bottom.as_str(), "ACGGCATT");
    assert_eq!(report.duplex.bottom_start, 0);
    assert_eq!(report.top_added_3prime, "GT");
    assert_eq!(report.bottom_added_3prime, "TT");
    assert!(report.top_removed_3prime.is_empty());
    assert!(report.bottom_removed_3prime.is_empty());
    let (left, right) = report.duplex.ends().unwrap();
    assert_eq!(left.polarity, Polarity::Blunt);
    assert_eq!(right.polarity, Polarity::Blunt);
    assert_eq!(report.duplex.left, source.left);
    assert_eq!(report.duplex.right, source.right);
    check_ancestry(&report);
}
#[test]
fn literal_three_prime_removal_preserves_five_prime_ends_and_all_retained_bases() {
    let source = duplex("TGCCGT", "GGCATT", -2);
    let report = repair_ends(&source, EndRepairMethod::RemoveThreePrimeOverhangs).unwrap();
    assert_eq!(report.duplex.top.as_str(), "TGCC");
    assert_eq!(report.duplex.bottom.as_str(), "GGCA");
    assert_eq!(report.duplex.bottom_start, 0);
    assert_eq!(report.top_removed_3prime, "GT");
    assert_eq!(report.bottom_removed_3prime, "TT");
    assert!(report.top_added_3prime.is_empty());
    assert!(report.bottom_added_3prime.is_empty());
    check_ancestry(&report);
}
#[test]
fn mixed_end_polarities_and_method_specific_unchanged_ends_are_explicit() {
    let source = duplex("AATGCCGT", "GGCA", 2);
    let filled = repair_ends(&source, EndRepairMethod::FillFivePrimeOverhangs).unwrap();
    assert_eq!(filled.duplex.top.as_str(), "AATGCCGT");
    assert_eq!(filled.duplex.bottom.as_str(), "GGCATT");
    assert_eq!(
        filled.duplex.ends().unwrap().1.polarity,
        Polarity::ThreePrime
    );
    let trimmed = repair_ends(&source, EndRepairMethod::RemoveThreePrimeOverhangs).unwrap();
    assert_eq!(trimmed.duplex.top.as_str(), "AATGCC");
    assert_eq!(trimmed.duplex.bottom.as_str(), "GGCA");
    assert_eq!(
        trimmed.duplex.ends().unwrap().0.polarity,
        Polarity::FivePrime
    );
    let blunt = repair_ends(&source, EndRepairMethod::Blunt).unwrap();
    assert_eq!(blunt.duplex.top.as_str(), "AATGCC");
    assert_eq!(blunt.duplex.bottom.as_str(), "GGCATT");
    assert_eq!(blunt.top_removed_3prime, "GT");
    assert_eq!(blunt.bottom_added_3prime, "TT");
    let opposite = repair_ends(&duplex("TGCC", "ACGGCATT", -2), EndRepairMethod::Blunt).unwrap();
    assert_eq!(opposite.duplex.top.as_str(), "TGCCGT");
    assert_eq!(opposite.duplex.bottom.as_str(), "ACGGCA");
    for report in [&filled, &trimmed, &blunt, &opposite] {
        check_ancestry(report);
    }
}
#[test]
fn repair_commutes_with_reversal_for_each_geometry_and_preserves_internal_nicks() {
    for source in [
        duplex("AATGCC", "ACGGCA", 2),
        duplex("TGCCGT", "GGCATT", -2),
        duplex("AATGCCGT", "GGCA", 2),
        duplex("TGCC", "ACGGCATT", -2),
    ] {
        for method in [
            EndRepairMethod::FillFivePrimeOverhangs,
            EndRepairMethod::RemoveThreePrimeOverhangs,
            EndRepairMethod::Blunt,
        ] {
            let report = repair_ends(&source, method).unwrap();
            let reversed = repair_ends(&source.reverse().unwrap(), method).unwrap();
            assert_eq!(report.duplex.reverse().unwrap(), reversed.duplex);
            check_ancestry(&report);
            check_ancestry(&reversed);
        }
    }
    let mut source = duplex("AATGCC", "ACGGCA", 2);
    source.nicks = vec![
        Nick {
            strand: Strand::Top,
            position: 3,
        },
        Nick {
            strand: Strand::Bottom,
            position: 4,
        },
    ];
    let repaired = repair_ends(&source, EndRepairMethod::Blunt).unwrap();
    assert_eq!(repaired.duplex.nicks, source.nicks);
    assert_eq!(
        repaired.duplex.reverse().unwrap(),
        repair_ends(&source.reverse().unwrap(), EndRepairMethod::Blunt)
            .unwrap()
            .duplex
    );
}
#[test]
fn blocked_unknown_termini_unpaired_nicks_and_circles_never_become_repaired_by_assumption() {
    for state in [ThreePrime::Blocked, ThreePrime::Unknown] {
        let mut source = duplex("AATGCC", "ACGGCA", 2);
        source.left.three_prime = state;
        assert!(repair_ends(&source, EndRepairMethod::Blunt).is_err());
        source.left.three_prime = ThreePrime::Hydroxyl;
        source.right.three_prime = state;
        assert!(repair_ends(&source, EndRepairMethod::FillFivePrimeOverhangs).is_err());
        // An inapplicable operation does not invent removal of a blocking group.
        assert_eq!(
            repair_ends(&source, EndRepairMethod::RemoveThreePrimeOverhangs)
                .unwrap()
                .duplex,
            source
        );
    }
    for nick in [
        Nick {
            strand: Strand::Top,
            position: 1,
        },
        Nick {
            strand: Strand::Top,
            position: 2,
        },
        Nick {
            strand: Strand::Bottom,
            position: 7,
        },
    ] {
        let mut source = duplex("AATGCC", "ACGGCA", 2);
        source.nicks.push(nick);
        source.validate().unwrap();
        assert!(repair_ends(&source, EndRepairMethod::Blunt).is_err());
    }
    let mut circle = Duplex::exact_pcr(&Sequence::dna("ACGT").unwrap()).unwrap();
    circle.topology = Topology::Circular;
    assert!(repair_ends(&circle, EndRepairMethod::Blunt).is_err());
    let unchanged =
        Duplex::from_molecule(&Molecule::dna("unknown", "ACGT", Topology::Linear).unwrap())
            .unwrap();
    assert_eq!(
        repair_ends(&unchanged, EndRepairMethod::Blunt)
            .unwrap()
            .duplex,
        unchanged
    );
}
#[test]
fn blunting_does_not_phosphorylate_and_ligation_still_checks_real_bonds() {
    let mut source = duplex("AATGCC", "ACGGCA", 2);
    source.left.five_prime = FivePrime::Hydroxyl;
    source.right.five_prime = FivePrime::Hydroxyl;
    let product = repair_ends(&source, EndRepairMethod::Blunt).unwrap().duplex;
    assert_eq!(product.left.five_prime, FivePrime::Hydroxyl);
    assert_eq!(product.right.five_prime, FivePrime::Hydroxyl);
    assert!(ligate(&[product.clone(), product.clone()], false, false).is_err());
    let phosphorylated = product.treat(EndTreatment::Phosphorylate).unwrap();
    let joined = ligate(&[phosphorylated.clone(), phosphorylated], false, false).unwrap();
    assert_eq!(joined.top.as_str(), "AATGCCGTAATGCCGT");
    assert!(joined.nicks.is_empty());
}
fn anneal_options(minimum: usize) -> AnnealOptions {
    AnnealOptions {
        min_paired_bases: minimum,
        ..AnnealOptions::default()
    }
}
#[test]
fn literal_annealing_geometry_and_independent_strand_chemistry_match_free_ends() {
    let mut options = anneal_options(4);
    options.top_chemistry = Terminus {
        five_prime: FivePrime::Hydroxyl,
        three_prime: ThreePrime::Blocked,
    };
    options.bottom_chemistry = Terminus::CLEAVED;
    let top = Sequence::dna("AATGCC").unwrap();
    let bottom = Sequence::dna("ACGGCA").unwrap();
    let report = anneal_oligos(&top, &bottom, &options).unwrap();
    assert!(report.complete && report.unique);
    assert_eq!(report.candidates.len(), 1);
    let candidate = &report.candidates[0];
    assert_eq!(candidate.duplex.bottom_start, 2);
    assert_eq!(
        (
            candidate.paired_start,
            candidate.paired_end,
            candidate.paired_bases
        ),
        (2, 6, 4)
    );
    let (left, right) = candidate.duplex.ends().unwrap();
    assert_eq!(
        (left.polarity, left.overhang.as_str()),
        (Polarity::FivePrime, "AA")
    );
    assert_eq!(
        (right.polarity, right.overhang.as_str()),
        (Polarity::FivePrime, "AC")
    );
    assert_eq!(left.chemistry.five_prime, FivePrime::Hydroxyl);
    assert_eq!(left.chemistry.three_prime, ThreePrime::Hydroxyl);
    assert_eq!(right.chemistry.five_prime, FivePrime::Phosphate);
    assert_eq!(right.chemistry.three_prime, ThreePrime::Blocked);
    assert!(repair_ends(&candidate.duplex, EndRepairMethod::Blunt).is_err());
    let recovered: AnnealReport =
        serde_json::from_slice(&serde_json::to_vec(&report).unwrap()).unwrap();
    assert_eq!(recovered, report);
}
fn words(length: usize) -> Vec<String> {
    let mut values = vec![String::new()];
    for _ in 0..length {
        values = values
            .into_iter()
            .flat_map(|prefix| {
                b"ACGT"
                    .iter()
                    .map(move |base| format!("{prefix}{}", *base as char))
            })
            .collect();
    }
    values
}
fn independent_offsets(top: &str, bottom: &str, minimum: usize) -> Vec<(i64, usize)> {
    let mut result = vec![];
    for offset in -(bottom.len() as i64) + 1..top.len() as i64 {
        let pairs = bottom
            .bytes()
            .enumerate()
            .filter_map(|(index, base)| {
                let coordinate = offset + (bottom.len() - 1 - index) as i64;
                (coordinate >= 0 && (coordinate as usize) < top.len())
                    .then_some((coordinate as usize, base))
            })
            .collect::<Vec<_>>();
        if pairs.len() >= minimum
            && pairs
                .iter()
                .all(|(coordinate, base)| top.as_bytes()[*coordinate] == manual_complement(*base))
        {
            result.push((offset, pairs.len()));
        }
    }
    result
}
#[test]
fn exhaustive_12288_independent_three_base_pair_and_threshold_oracle_cases() {
    let values = words(3);
    for top in &values {
        for bottom in &values {
            for minimum in 1..=3 {
                let report = anneal_oligos(
                    &Sequence::dna(top).unwrap(),
                    &Sequence::dna(bottom).unwrap(),
                    &anneal_options(minimum),
                )
                .unwrap();
                let expected = independent_offsets(top, bottom, minimum);
                let actual = report
                    .candidates
                    .iter()
                    .map(|candidate| (candidate.duplex.bottom_start, candidate.paired_bases))
                    .collect::<Vec<_>>();
                assert_eq!(actual, expected, "top={top} bottom={bottom} min={minimum}");
                assert!(report.complete);
                assert_eq!(report.unique, expected.len() == 1);
            }
        }
    }
}
#[test]
fn repeated_annealing_offsets_and_budget_limits_never_create_a_longest_or_unique_claim() {
    let top = Sequence::dna("AAAA").unwrap();
    let bottom = Sequence::dna("TTTT").unwrap();
    let report = anneal_oligos(&top, &bottom, &anneal_options(2)).unwrap();
    assert!(report.complete);
    assert!(!report.unique);
    assert_eq!(
        report
            .candidates
            .iter()
            .map(|candidate| candidate.paired_bases)
            .collect::<Vec<_>>(),
        vec![2, 3, 4, 3, 2]
    );
    let mut options = anneal_options(2);
    options.max_candidates = 1;
    let capped = anneal_oligos(&top, &bottom, &options).unwrap();
    assert!(!capped.complete && !capped.unique);
    assert_eq!(capped.candidates[0].duplex.bottom_start, -2);
    options = anneal_options(2);
    options.max_comparisons = 1;
    let partial = anneal_oligos(&top, &bottom, &options).unwrap();
    assert!(partial.candidates.is_empty());
    assert!(!partial.complete && !partial.unique);
    options.max_comparisons = 2;
    let first = anneal_oligos(&top, &bottom, &options).unwrap();
    assert_eq!(first.candidates.len(), 1);
    assert!(!first.complete && !first.unique);
    assert_eq!(first.candidates[0].paired_bases, 2);
    let mut exact = anneal_options(4);
    exact.max_candidates = 1;
    let full = anneal_oligos(&top, &bottom, &exact).unwrap();
    assert!(full.complete && full.unique);
}
#[test]
fn swapped_oligos_reverse_the_entire_physical_candidate_set() {
    let top = Sequence::dna("AATGCC").unwrap();
    let bottom = Sequence::dna("ACGGCA").unwrap();
    let mut options = anneal_options(1);
    options.top_chemistry = Terminus::PCR;
    options.bottom_chemistry = Terminus::CLEAVED;
    let original = anneal_oligos(&top, &bottom, &options).unwrap();
    std::mem::swap(&mut options.top_chemistry, &mut options.bottom_chemistry);
    let swapped = anneal_oligos(&bottom, &top, &options).unwrap();
    assert!(original.complete && swapped.complete);
    for candidate in original.candidates {
        let reverse = candidate.duplex.reverse().unwrap();
        assert!(
            swapped
                .candidates
                .iter()
                .any(|other| other.duplex == reverse)
        );
    }
}
#[test]
fn unsupported_alphabets_ambiguity_and_impossible_limits_do_not_invent_physical_states() {
    for invalid in [
        Sequence::dna("ACGN").unwrap(),
        Sequence::new(Alphabet::Rna, "ACGU").unwrap(),
    ] {
        assert!(
            anneal_oligos(
                &invalid,
                &Sequence::dna("ACGT").unwrap(),
                &anneal_options(1)
            )
            .is_err()
        );
    }
    let top = Sequence::dna("AAAA").unwrap();
    let bottom = Sequence::dna("AAAA").unwrap();
    let no_match = anneal_oligos(&top, &bottom, &anneal_options(1)).unwrap();
    assert!(no_match.complete && !no_match.unique && no_match.candidates.is_empty());
    let too_short = anneal_oligos(&top, &bottom, &anneal_options(5)).unwrap();
    assert!(too_short.complete && too_short.candidates.is_empty());
    let mut options = anneal_options(0);
    assert!(anneal_oligos(&top, &bottom, &options).is_err());
    options = anneal_options(1);
    options.max_candidates = MAX_ANNEAL_CANDIDATES + 1;
    assert!(anneal_oligos(&top, &bottom, &options).is_err());
    let long = Sequence::dna(&"A".repeat(MAX_ANNEAL_OLIGO_BASES + 1)).unwrap();
    assert!(anneal_oligos(&long, &bottom, &anneal_options(1)).is_err());
}

#[test]
fn every_interval_pair_on_independent_reference_slices_has_expected_repair_geometry() {
    for reference in ["ACGTTA", "AACGCC", "GTTACA"] {
        for ta in 0..reference.len() {
            for tb in ta + 1..=reference.len() {
                for ba in 0..reference.len() {
                    for bb in ba + 1..=reference.len() {
                        if ta.max(ba) >= tb.min(bb) {
                            continue;
                        }
                        let source = duplex(
                            &reference[ta..tb],
                            &manual_rc(&reference[ba..bb]),
                            ba as i64 - ta as i64,
                        );
                        for method in [
                            EndRepairMethod::FillFivePrimeOverhangs,
                            EndRepairMethod::RemoveThreePrimeOverhangs,
                            EndRepairMethod::Blunt,
                        ] {
                            let (new_tb, new_ba) = match method {
                                EndRepairMethod::FillFivePrimeOverhangs => (tb.max(bb), ba.min(ta)),
                                EndRepairMethod::RemoveThreePrimeOverhangs => {
                                    (tb.min(bb), ba.max(ta))
                                }
                                EndRepairMethod::Blunt => (bb, ta),
                            };
                            let repaired = repair_ends(&source, method).unwrap();
                            assert_eq!(repaired.duplex.top.as_str(), &reference[ta..new_tb]);
                            assert_eq!(
                                repaired.duplex.bottom.as_str(),
                                manual_rc(&reference[new_ba..bb])
                            );
                            assert_eq!(repaired.duplex.bottom_start, new_ba as i64 - ta as i64);
                            check_ancestry(&repaired);
                        }
                    }
                }
            }
        }
    }
}
#[test]
fn oversized_end_repair_report_is_rejected_before_output_strands_are_duplicated() {
    let source = duplex(&"A".repeat(9_000_000), &"T".repeat(9_000_000), 0);
    assert!(matches!(
        repair_ends(&source, EndRepairMethod::Blunt),
        Err(Error::TooLarge)
    ));
}

#[test]
fn nested_unknown_chemical_semantics_are_not_discarded_from_annealing_options() {
    let mut value = serde_json::to_value(AnnealOptions::default()).unwrap();
    value["top_chemistry"]["five_prime_protected"] = serde_json::json!(true);
    assert!(serde_json::from_value::<AnnealOptions>(value).is_err());
}
