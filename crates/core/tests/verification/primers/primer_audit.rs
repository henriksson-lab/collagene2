//! Independent review of physical orientation, mass-action units and bounded endpoint ranking.
//! The production diagonal screens are checked against dynamic common-substring oracles;
//! binding evidence is reconstructed in each oligo's own biological orientation.
use collagene_core::{
    analysis::{Primer, simulate_pcr},
    primers::*,
};
use collagene_model::{Error, Molecule, Topology};

fn primer(text: &str) -> Primer {
    Primer::new("audit", text, "").unwrap()
}
fn near(actual: f64, expected: f64) {
    assert!((actual - expected).abs() < 1e-9, "{actual} != {expected}");
}
fn concrete(b: u8) -> Vec<u8> {
    match b {
        b'A' => vec![b'A'],
        b'C' => vec![b'C'],
        b'G' => vec![b'G'],
        b'T' => vec![b'T'],
        b'R' => vec![b'A', b'G'],
        b'Y' => vec![b'C', b'T'],
        b'N' => b"ACGT".to_vec(),
        _ => panic!("audit alphabet"),
    }
}
fn complement(b: u8) -> u8 {
    match b {
        b'A' => b'T',
        b'C' => b'G',
        b'G' => b'C',
        b'T' => b'A',
        b'R' => b'Y',
        b'Y' => b'R',
        b'N' => b'N',
        _ => panic!(),
    }
}
fn rc(text: &str) -> String {
    text.bytes().rev().map(|b| complement(b) as char).collect()
}
fn compatible(a: u8, b: u8) -> bool {
    concrete(a).iter().any(|a| concrete(b).contains(a))
}
fn words(alphabet: &[u8], length: usize) -> Vec<String> {
    let mut all = vec![String::new()];
    for _ in 0..length {
        all = all
            .into_iter()
            .flat_map(|prefix| {
                alphabet
                    .iter()
                    .map(move |b| format!("{prefix}{}", *b as char))
            })
            .collect();
    }
    all
}

#[test]
fn manual_stack_sums_satisfy_independent_half_melt_mass_action_for_both_stoichiometries() {
    // SantaLucia1998 Table2: AT/TG/GC/CG/GA plus terminal A/A;
    // AG/GC/CT plus terminal A/T and C2 symmetry, respectively.
    // These factual sums are independent of the production dinucleotide lookup.
    for (text, h, s, self_complementary) in [
        ("ATGCGA", -39.7, -108.7, false),
        ("AGCT", -20.8, -59.6, true),
    ] {
        for sodium in [0.01, 0.05, 1.] {
            for total in [1e-9, 500e-9, 1e-4] {
                let conditions = ThermoConditions {
                    sodium_molar: sodium,
                    total_strand_molar: total,
                };
                let report = analyze_primer(&primer(text), &conditions)
                    .unwrap()
                    .tm
                    .unwrap();
                near(report.delta_h_kcal_per_mol, h);
                near(report.delta_s_cal_per_mol_kelvin_at_1m, s);
                assert_eq!(report.self_complementary, self_complementary);
                let entropy = s + 0.368 * (text.len() - 1) as f64 * sodium.ln();
                let kelvin = report.tm_c + 273.15;
                let association = (-(1000. * h - kelvin * entropy) / (1.987 * kelvin)).exp();
                // At half melting, one quarter of total strand molarity is duplex molecules.
                // Self pairing has [free S]=Ct/2; equal A/B has [A]=[B]=Ct/4.
                let free = if self_complementary {
                    total / 2.
                } else {
                    total / 4.
                };
                let inferred_duplex = association * free * free;
                assert!(
                    (inferred_duplex / (total / 4.) - 1.).abs() < 1e-10,
                    "{text}: mass action failed"
                );
            }
        }
    }
}

// Longest common substring of A and reverse-complement(B), an algorithm independent
// of the production placement/run scanner. j_start==0 means B's physical 3′ end.
fn common_substring(a: &str, b: &str) -> (usize, usize, usize) {
    let reverse = rc(b);
    let (a, b) = (a.as_bytes(), reverse.as_bytes());
    let mut dp = vec![vec![0usize; b.len() + 1]; a.len() + 1];
    let (mut best, mut either_terminal, mut b_terminal) = (0, 0, 0);
    for i in 0..a.len() {
        for j in 0..b.len() {
            if compatible(a[i], b[j]) {
                let length = dp[i][j] + 1;
                dp[i + 1][j + 1] = length;
                best = best.max(length);
                let starts_at_b_three_prime = j + 1 == length;
                if i + 1 == a.len() || starts_at_b_three_prime {
                    either_terminal = either_terminal.max(length);
                }
                if starts_at_b_three_prime {
                    b_terminal = b_terminal.max(length);
                }
            }
        }
    }
    (best, either_terminal, b_terminal)
}
#[test]
fn independent_dynamic_oracle_covers_ambiguous_antiparallel_and_physical_three_prime_runs() {
    let all = words(b"ACGTRN", 3);
    for a in &all {
        for b in &all {
            let report = heterodimer(&primer(a), &primer(b)).unwrap();
            let (best, terminal, _) = common_substring(a, b);
            assert_eq!(
                (report.max_run, report.max_three_prime_run),
                (best, terminal),
                "{a} / {b}"
            );
            if let Some(run) = report.strongest {
                assert_eq!(run.a_end - run.a_start, run.length);
                assert_eq!(run.b_end - run.b_start, run.length);
                assert_eq!(run.length, best);
                let left = &a[run.a_start..run.a_end];
                let right = rc(&b[run.b_start..run.b_end]);
                assert!(
                    left.bytes()
                        .zip(right.bytes())
                        .all(|(a, b)| compatible(a, b))
                );
                assert_eq!(
                    run.certain,
                    left.bytes()
                        .chain(right.bytes())
                        .all(|b| concrete(b).len() == 1)
                );
                assert_eq!(run.a_three_prime, run.a_end == a.len());
                assert_eq!(run.b_three_prime, run.b_end == b.len());
            } else {
                assert_eq!(best, 0);
            }
        }
    }
}
#[test]
fn hairpins_match_all_legal_loop_partitions_without_counting_a_virtual_prefix_terminus() {
    for text in words(b"ATN", 8).into_iter().chain([
        "GCGAAACGC".into(),
        "GCGAAACGCTTT".into(),
        "TTTGCGAAACGC".into(),
    ]) {
        let (mut best, mut terminal) = (0, 0);
        // Any legal stem has a partition separating two arms by at least three bases.
        // The virtual end of a prefix is not a physical 3′ end of this molecule.
        for cut in 1..text.len() - 3 {
            let (run, _, right_terminal) = common_substring(&text[..cut], &text[cut + 3..]);
            best = best.max(run);
            terminal = terminal.max(right_terminal);
        }
        let report = analyze_primer(&primer(&text), &ThermoConditions::default())
            .unwrap()
            .hairpin;
        assert_eq!(
            (report.max_run, report.max_three_prime_run),
            (best, terminal),
            "{text}"
        );
        if let Some(run) = report.strongest {
            assert!(run.a_end + 3 <= run.b_start, "{text}: arms overlap loop");
            assert!(!run.a_three_prime);
            assert_eq!(run.b_three_prime, run.b_end == text.len());
        }
    }
}

#[derive(Debug, PartialEq, Eq)]
struct ExpectedHit {
    start: usize,
    reverse: bool,
    mismatches: Vec<usize>,
    ambiguous: Vec<usize>,
    terminal: usize,
    extracted: String,
}
fn binding_oracle(
    reference: &str,
    oligo: &str,
    circular: bool,
    allowed: usize,
    required: usize,
) -> Vec<ExpectedHit> {
    let mut hits = vec![];
    let limit = if circular {
        reference.len()
    } else {
        reference.len() - oligo.len() + 1
    };
    for start in 0..limit {
        for reverse in [false, true] {
            let physical: String = (0..oligo.len())
                .map(|i| reference.as_bytes()[(start + i) % reference.len()] as char)
                .collect();
            let extracted = if reverse { rc(&physical) } else { physical };
            let (mut mismatches, mut ambiguous, mut certain) = (vec![], vec![], vec![]);
            for (i, (a, b)) in oligo.bytes().zip(extracted.bytes()).enumerate() {
                if !compatible(a, b) {
                    mismatches.push(i);
                    certain.push(false);
                } else if concrete(a).len() != 1 || concrete(b).len() != 1 {
                    ambiguous.push(i);
                    certain.push(false);
                } else {
                    certain.push(true);
                }
            }
            let terminal = certain.into_iter().rev().take_while(|v| *v).count();
            if mismatches.len() <= allowed && terminal >= required {
                hits.push(ExpectedHit {
                    start,
                    reverse,
                    mismatches,
                    ambiguous,
                    terminal,
                    extracted,
                });
            }
        }
    }
    hits
}
#[test]
fn exhaustive_small_binding_oracle_checks_reverse_circular_mismatch_and_ambiguity_coordinates() {
    let queries = words(b"ACGTR", 2);
    for reference in words(b"ACGTN", 3) {
        for circular in [false, true] {
            let molecule = Molecule::dna(
                "audit",
                &reference,
                if circular {
                    Topology::Circular
                } else {
                    Topology::Linear
                },
            )
            .unwrap();
            for query in &queries {
                for (allowed, required) in [(0, 0), (1, 1), (2, 2)] {
                    let options = BindingOptions {
                        max_mismatches: allowed,
                        min_three_prime_exact: required,
                        ..BindingOptions::default()
                    };
                    let report = find_binding_sites(&molecule, &primer(query), &options).unwrap();
                    assert!(!report.truncated);
                    let actual = report
                        .hits
                        .iter()
                        .map(|h| {
                            assert_eq!(h.crosses_origin, h.start + query.len() > reference.len());
                            assert_eq!(
                                h.certain,
                                h.mismatches.is_empty() && h.ambiguous_positions.is_empty()
                            );
                            ExpectedHit {
                                start: h.start,
                                reverse: h.reverse,
                                mismatches: h.mismatches.clone(),
                                ambiguous: h.ambiguous_positions.clone(),
                                terminal: h.three_prime_exact,
                                extracted: h.location.extract(&molecule.sequence).unwrap(),
                            }
                        })
                        .collect::<Vec<_>>();
                    assert_eq!(
                        actual,
                        binding_oracle(&reference, query, circular, allowed, required),
                        "{reference}/{query} circular={circular}"
                    );
                }
            }
        }
    }
}

fn open_design(start: usize, end: usize) -> DesignOptions {
    DesignOptions {
        target_start: start,
        target_end: end,
        min_length: 8,
        max_length: 11,
        max_pairs: 4,
        min_gc_fraction: 0.,
        max_gc_fraction: 1.,
        max_homopolymer: 200,
        max_three_prime_complementarity: 200,
        max_tm_difference_c: 100.,
        forward_tail: " g g\t".into(),
        reverse_tail: " t c".into(),
        ..DesignOptions::default()
    }
}
#[test]
fn all_circular_origins_have_ranked_exact_endpoint_geometry_and_actual_tailed_pcr_products() {
    let reference = "ACGTAGTCCATGACCTGAGTCACTGACCTAGT";
    let molecule = Molecule::dna("audit", reference, Topology::Circular).unwrap();
    for start in 0..reference.len() {
        for span in [16, 24, 32] {
            let target: String = (0..span)
                .map(|i| reference.as_bytes()[(start + i) % reference.len()] as char)
                .collect();
            let options = open_design(start, start + span);
            let mut expected = vec![];
            for f in 8..=11 {
                for r in 8..=11 {
                    if f + r > span {
                        continue;
                    }
                    let ftm = analyze_primer(&primer(&target[..f]), &options.conditions)
                        .unwrap()
                        .tm
                        .unwrap()
                        .tm_c;
                    let rtm =
                        analyze_primer(&primer(&rc(&target[span - r..])), &options.conditions)
                            .unwrap()
                            .tm
                            .unwrap()
                            .tm_c;
                    expected.push((
                        (ftm - options.target_tm_c).abs()
                            + (rtm - options.target_tm_c).abs()
                            + (ftm - rtm).abs(),
                        f,
                        r,
                    ));
                }
            }
            expected.sort_by(|a, b| a.0.total_cmp(&b.0));
            expected.truncate(options.max_pairs);
            let report = design_primers(&molecule, &options).unwrap();
            assert_eq!(report.evaluated_pairs, 16);
            assert_eq!(report.pairs.len(), expected.len());
            for (pair, (score, f, r)) in report.pairs.iter().zip(expected) {
                near(pair.score, score);
                assert_eq!(pair.forward.primer.annealing.as_str(), &target[..f]);
                assert_eq!(
                    pair.reverse.primer.annealing.as_str(),
                    rc(&target[span - r..])
                );
                assert_eq!(pair.forward.primer.tail, "GG");
                assert_eq!(pair.reverse.primer.tail, "TC");
                assert_eq!(pair.product_length, span + 4);
                for candidate in [&pair.forward, &pair.reverse] {
                    assert_eq!(
                        candidate.location.extract(&molecule.sequence).unwrap(),
                        candidate.primer.annealing.as_str()
                    );
                }
                let expected_product = format!("GG{target}GA");
                let products =
                    simulate_pcr(&molecule, &pair.forward.primer, &pair.reverse.primer).unwrap();
                assert!(
                    products
                        .iter()
                        .any(|p| p.sequence.as_str() == expected_product),
                    "source origin={start}, target span={span}, f={f},r={r}"
                );
            }
        }
    }
}

#[test]
fn exact_hit_limit_is_not_mislabeled_truncated_when_no_additional_site_exists() {
    let molecule = Molecule::dna("audit", "AACG", Topology::Linear).unwrap();
    let options = BindingOptions {
        max_hits: 1,
        max_work: 16,
        ..BindingOptions::default()
    };
    let report = find_binding_sites(&molecule, &primer("AAC"), &options).unwrap();
    assert_eq!(report.hits.len(), 1);
    assert_eq!(report.placements_examined, 4);
    assert!(!report.truncated);
    let report = find_binding_sites(
        &molecule,
        &primer("AAC"),
        &BindingOptions {
            max_work: 3,
            ..options
        },
    )
    .unwrap();
    assert_eq!(report.hits.len(), 1);
    assert_eq!(report.placements_examined, 1);
    assert!(report.truncated);
}

#[test]
fn pair_work_guard_rejects_a_valid_but_excessive_cartesian_search_before_pair_materialization() {
    let molecule = Molecule::dna(
        "audit",
        &"ACGTAGTCCATGACCTGAGTCACTGACCTAGT".repeat(10),
        Topology::Linear,
    )
    .unwrap();
    let options = DesignOptions {
        target_start: 0,
        target_end: molecule.sequence.len(),
        min_length: 8,
        max_length: 71,
        forward_tail: "ACGT".repeat(30),
        reverse_tail: "TGCA".repeat(30),
        min_gc_fraction: 0.,
        max_gc_fraction: 1.,
        max_homopolymer: 200,
        max_three_prime_complementarity: 200,
        max_tm_difference_c: 100.,
        ..DesignOptions::default()
    };
    // All 64 lengths at each endpoint are individually legal, but even the shortest
    // 128-base ordered pair requires more than 134 million comparisons over that grid.
    let error = design_primers(&molecule, &options).unwrap_err();
    assert!(matches!(error, Error::Unsupported(ref message) if message.contains("20 million")));
}
