//! Manually calculated parameters and independent exhaustive interval oracles.
use collagene_core::{
    analysis::{Primer, simulate_pcr},
    primers::*,
};
use collagene_model::{Alphabet, Feature, Location, Molecule, Sequence, Topology};
fn oligo(s: &str) -> Primer {
    Primer::new("fixture", s, "").unwrap()
}
fn molecule(s: &str, t: Topology) -> Molecule {
    Molecule::dna("primer fixture", s, t).unwrap()
}
fn near(a: f64, b: f64) {
    assert!((a - b).abs() < 1e-9, "{a} != {b}");
}
#[test]
fn manual_nonself_nearest_neighbor_fixture_and_sodium_entropy_shift() {
    // CG, GT, TT, TG, GA stacks; terminal G and A initiation.
    let p = oligo("CGTTGA");
    let mut c = ThermoConditions {
        sodium_molar: 1.,
        total_strand_molar: 500e-9,
    };
    let a = analyze_primer(&p, &c).unwrap();
    let tm = a.tm.unwrap();
    near(tm.delta_h_kcal_per_mol, -41.2);
    near(tm.delta_s_cal_per_mol_kelvin_at_1m, -115.4);
    near(
        tm.tm_c,
        -41200. / (-115.4 + 1.987 * (125e-9f64).ln()) - 273.15,
    );
    assert!(!tm.self_complementary);
    c.sodium_molar = 0.05;
    let low = analyze_primer(&p, &c).unwrap().tm.unwrap();
    near(
        low.salt_corrected_delta_s,
        -115.4 + 0.368 * 5. * 0.05f64.ln(),
    );
    assert!(low.tm_c < tm.tm_c);
}
#[test]
fn self_complementary_symmetry_and_concentration_factor_are_distinct() {
    let c = ThermoConditions {
        sodium_molar: 1.,
        total_strand_molar: 500e-9,
    };
    let t = analyze_primer(&oligo("CGCG"), &c).unwrap().tm.unwrap();
    near(t.delta_h_kcal_per_mol, -30.8);
    near(t.delta_s_cal_per_mol_kelvin_at_1m, -85.8);
    near(t.tm_c, -4.456389990619925);
    assert!(t.self_complementary);
    let reverse = Sequence::dna("CGTTGA")
        .unwrap()
        .reverse_complement()
        .unwrap();
    near(
        analyze_primer(&oligo(reverse.as_str()), &c)
            .unwrap()
            .tm
            .unwrap()
            .tm_c,
        analyze_primer(&oligo("CGTTGA"), &c)
            .unwrap()
            .tm
            .unwrap()
            .tm_c,
    );
}
#[test]
fn additional_manual_stacks_cover_ta_ct_gc_and_gg_parameters() {
    let conditions = ThermoConditions {
        sodium_molar: 1.,
        ..ThermoConditions::default()
    };
    // TA+AT+TA+AG+GC with terminal T/C initiation.
    let t = analyze_primer(&oligo("TATAGC"), &conditions)
        .unwrap()
        .tm
        .unwrap();
    near(t.delta_h_kcal_per_mol, -36.8);
    near(t.delta_s_cal_per_mol_kelvin_at_1m, -107.1);
    // CC+CG+GG, terminal C/G initiation, and self-complementarity entropy.
    let t = analyze_primer(&oligo("CCGG"), &conditions)
        .unwrap()
        .tm
        .unwrap();
    near(t.delta_h_kcal_per_mol, -26.4);
    near(t.delta_s_cal_per_mol_kelvin_at_1m, -74.0);
}
#[test]
fn ambiguity_produces_gc_bounds_and_possible_structure_without_mixture_tm() {
    let p = Primer::new("degenerate", "GSN", "AT").unwrap();
    let a = analyze_primer(&p, &ThermoConditions::default()).unwrap();
    near(a.annealing_gc.minimum_fraction, 2. / 3.);
    near(a.annealing_gc.maximum_fraction, 1.);
    near(a.ordered_gc.minimum_fraction, 0.4);
    near(a.ordered_gc.maximum_fraction, 0.6);
    assert!(a.tm.is_none());
    assert!(a.tm_unavailable.unwrap().contains("unambiguous"));
    assert_eq!(a.ordered_sequence, "ATGSN");
    let possible = heterodimer(&oligo("NNN"), &oligo("AAA")).unwrap();
    assert_eq!(possible.max_run, 3);
    assert!(!possible.strongest.unwrap().certain);
}
#[test]
fn tail_is_excluded_from_initial_tm_but_included_in_structure_and_homopolymer() {
    let c = ThermoConditions::default();
    let simple = analyze_primer(&oligo("ACGTAGTC"), &c).unwrap();
    let tailed = analyze_primer(&Primer::new("tailed", "ACGTAGTC", "CCCCCC").unwrap(), &c).unwrap();
    near(simple.tm.unwrap().tm_c, tailed.tm.unwrap().tm_c);
    assert_eq!(tailed.homopolymer_run, 6);
    let pair = heterodimer(&Primer::new("a", "AAAA", "GGGG").unwrap(), &oligo("CCCC")).unwrap();
    assert_eq!(pair.max_run, 4);
    assert_eq!(pair.max_three_prime_run, 4);
}
fn complement(b: u8) -> u8 {
    match b {
        b'A' => b'T',
        b'T' => b'A',
        b'C' => b'G',
        b'G' => b'C',
        _ => panic!(),
    }
}
fn oracle(a: &str, b: &str, hairpin: bool) -> (usize, usize) {
    let (mut longest, mut terminal) = (0, 0);
    for ast in 0..a.len() {
        for bst in 0..b.len() {
            for len in 1..=a.len().min(b.len()) {
                if ast + len > a.len() || bst + len > b.len() {
                    continue;
                }
                if hairpin && ast + len + 3 > bst {
                    continue;
                }
                if (0..len)
                    .all(|j| a.as_bytes()[ast + j] == complement(b.as_bytes()[bst + len - 1 - j]))
                {
                    longest = longest.max(len);
                    if ast + len == a.len() || bst + len == b.len() {
                        terminal = terminal.max(len);
                    }
                }
            }
        }
    }
    (longest, terminal)
}
#[test]
fn contiguous_structure_screens_match_exhaustively_enumerated_antiparallel_intervals() {
    let words = (0..64)
        .map(|mut value| {
            let mut s = String::new();
            for _ in 0..3 {
                s.push(b"ACGT"[value % 4] as char);
                value /= 4;
            }
            s
        })
        .collect::<Vec<_>>();
    for a in &words {
        for b in &words {
            let result = heterodimer(&oligo(a), &oligo(b)).unwrap();
            assert_eq!(
                (result.max_run, result.max_three_prime_run),
                oracle(a, b, false),
                "{a} {b}"
            );
        }
    }
    for s in ["GCGAAACGC", "AAATTT", "GCGTTTTCGC", "ACGTACGT", "NNN"] {
        if s.contains('N') {
            continue;
        }
        let result = analyze_primer(&oligo(s), &ThermoConditions::default())
            .unwrap()
            .hairpin;
        assert_eq!(
            (result.max_run, result.max_three_prime_run),
            oracle(s, s, true),
            "{s}"
        );
    }
}
#[test]
fn binding_locations_retain_primer_orientation_and_circular_origin() {
    let m = molecule("AACCGGTT", Topology::Circular);
    let f = find_binding_sites(&m, &oligo("TTAA"), &BindingOptions::default()).unwrap();
    assert!(!f.truncated);
    assert_eq!(f.hits.len(), 2); // self-complementary oligo has two physical orientations
    for h in f.hits {
        assert_eq!(h.start, 6);
        assert!(h.crosses_origin);
        assert_eq!(h.location.extract(&m.sequence).unwrap(), "TTAA");
        assert_eq!(h.three_prime_exact, 4);
    }
    let rev = find_binding_sites(&m, &oligo("CCGGTT"), &BindingOptions::default()).unwrap();
    assert!(rev.hits.iter().any(|h| h.reverse && h.start == 0));
    for h in rev.hits {
        assert_eq!(h.location.extract(&m.sequence).unwrap(), "CCGGTT");
    }
}
#[test]
fn mismatch_positions_and_three_prime_filter_use_oligo_not_reference_orientation() {
    let m = molecule("AACCGGTT", Topology::Linear);
    let opts = BindingOptions {
        max_mismatches: 1,
        min_three_prime_exact: 3,
        ..BindingOptions::default()
    };
    // reverse ideal CCGGTT from source AACCGG. Changed primer5′ base C->A.
    let result = find_binding_sites(&m, &oligo("ACGGTT"), &opts).unwrap();
    let h = result
        .hits
        .iter()
        .find(|h| h.reverse && h.start == 0)
        .unwrap();
    assert_eq!(h.mismatches, vec![0]);
    assert_eq!(h.three_prime_exact, 5);
    assert!(!h.certain);
    let changed_three_prime = find_binding_sites(&m, &oligo("CCGGTA"), &opts).unwrap();
    assert!(
        !changed_three_prime
            .hits
            .iter()
            .any(|h| h.reverse && h.start == 0)
    );
    let ambiguous = find_binding_sites(
        &molecule("ACNT", Topology::Linear),
        &oligo("ACGT"),
        &BindingOptions::default(),
    )
    .unwrap();
    assert_eq!(ambiguous.hits[0].ambiguous_positions, vec![2]);
    assert!(!ambiguous.hits[0].certain);
}
#[test]
fn hit_and_work_limits_are_visible_not_negative_specificity_claims() {
    let m = molecule("AAAAAAAA", Topology::Linear);
    let mut opts = BindingOptions {
        max_hits: 1,
        ..BindingOptions::default()
    };
    let r = find_binding_sites(&m, &oligo("AA"), &opts).unwrap();
    assert!(r.truncated);
    assert_eq!(r.hits.len(), 1);
    opts.max_hits = 100;
    opts.max_work = 2;
    let r = find_binding_sites(&m, &oligo("AA"), &opts).unwrap();
    assert!(r.truncated);
    assert_eq!(r.placements_examined, 1);
    let r = find_binding_sites(&m, &oligo("AAAAAAAAA"), &BindingOptions::default()).unwrap();
    assert!(!r.truncated);
    assert!(r.hits.is_empty());
}
fn design_options(start: usize, end: usize) -> DesignOptions {
    DesignOptions {
        target_start: start,
        target_end: end,
        min_length: 8,
        max_length: 8,
        min_gc_fraction: 0.,
        max_gc_fraction: 1.,
        max_tm_difference_c: 100.,
        max_homopolymer: 200,
        max_three_prime_complementarity: 200,
        forward_tail: "GG".into(),
        reverse_tail: "TT".into(),
        ..DesignOptions::default()
    }
}
#[test]
fn anchored_design_has_exact_locations_tails_and_independently_reconstructed_pcr_product() {
    let text = "ACGTAGTCCATGACCTGAGTCACTGACCTAGT";
    let m = molecule(text, Topology::Linear);
    let report = design_primers(&m, &design_options(0, text.len())).unwrap();
    assert_eq!(report.pairs.len(), 1);
    let pair = &report.pairs[0];
    assert_eq!(pair.forward.primer.annealing.as_str(), "ACGTAGTC");
    assert_eq!(pair.reverse.primer.annealing.as_str(), "ACTAGGTC");
    for c in [&pair.forward, &pair.reverse] {
        assert_eq!(
            c.location.extract(&m.sequence).unwrap(),
            c.primer.annealing.as_str()
        );
    }
    assert_eq!(pair.product_length, text.len() + 4);
    let products = simulate_pcr(&m, &pair.forward.primer, &pair.reverse.primer).unwrap();
    assert!(
        products
            .iter()
            .any(|p| p.sequence.as_str() == format!("GG{text}AA"))
    );
    assert!(report.warnings.iter().any(|w| w.contains("uniqueness")));
}
#[test]
fn circular_design_uses_unwrapped_interval_and_full_complement_join_location() {
    let text = "ACGTAGTCCATGACCTGAGTCACTGACCTAGT";
    let m = molecule(text, Topology::Circular);
    let n = text.len();
    let result = design_primers(&m, &design_options(n - 4, 2 * n - 4)).unwrap();
    let pair = &result.pairs[0];
    assert!(matches!(pair.forward.location, Location::Join(_)));
    let expected = format!("{}{}", &text[n - 4..], &text[..4]);
    assert_eq!(pair.forward.primer.annealing.as_str(), expected);
    for candidate in [&pair.forward, &pair.reverse] {
        assert_eq!(
            candidate.location.extract(&m.sequence).unwrap(),
            candidate.primer.annealing.as_str()
        );
    }
}
#[test]
fn design_ranking_and_output_limits_are_reproducible_and_do_not_invent_candidates() {
    let m = molecule("ACGTAGTCCATGACCTGAGTCACTGACCTAGT", Topology::Linear);
    let mut opts = design_options(0, m.sequence.len());
    opts.max_length = 12;
    opts.max_pairs = 3;
    let r = design_primers(&m, &opts).unwrap();
    assert_eq!(r.pairs.len(), 3);
    assert_eq!(r.evaluated_pairs, 25);
    assert!(r.pairs.windows(2).all(|p| p[0].score <= p[1].score));
    assert_eq!(
        serde_json::to_value(&r).unwrap(),
        serde_json::to_value(design_primers(&m, &opts).unwrap()).unwrap()
    );
    opts.min_gc_fraction = 1.;
    opts.max_gc_fraction = 1.;
    assert!(design_primers(&m, &opts).unwrap().pairs.is_empty());
}
#[test]
fn malformed_options_alphabets_and_oversized_oligos_fail_explicitly() {
    let mut c = ThermoConditions::default();
    c.sodium_molar = f64::NAN;
    assert!(analyze_primer(&oligo("ACGT"), &c).is_err());
    c = ThermoConditions::default();
    c.total_strand_molar = 0.;
    assert!(analyze_primer(&oligo("ACGT"), &c).is_err());
    assert!(analyze_primer(&oligo(&"A".repeat(201)), &ThermoConditions::default()).is_err());
    let mut p = oligo("ACGT");
    p.annealing = Sequence::new(Alphabet::Rna, "ACGU").unwrap();
    assert!(analyze_primer(&p, &ThermoConditions::default()).is_err());
    let m = molecule("ACGTACGTACGTACGT", Topology::Linear);
    let mut opts = design_options(0, 100);
    assert!(design_primers(&m, &opts).is_err());
    opts.target_end = 16;
    opts.min_length = 0;
    assert!(design_primers(&m, &opts).is_err());
    opts.min_length = 8;
    opts.max_gc_fraction = f64::INFINITY;
    assert!(design_primers(&m, &opts).is_err());
    let mut bindings = BindingOptions::default();
    bindings.max_work = usize::MAX;
    assert!(find_binding_sites(&m, &oligo("ACGT"), &bindings).is_err());
}
#[test]
fn saved_primer_binding_fields_survive_annotated_genbank_exchange() {
    let mut m = molecule("ACGTAGTCCATGACCTGAGTCACTGACCTAGT", Topology::Linear);
    m.features.push(Feature {
        kind: "primer_bind".into(),
        location: Location::range(0, 8),
        qualifiers: vec![
            ("label".into(), Some("tailed forward".into())),
            ("oligo_sequence".into(), Some("GGACGTAGTC".into())),
            ("annealing_sequence".into(), Some("ACGTAGTC".into())),
            ("5prime_tail".into(), Some("GG".into())),
            ("source_revision".into(), Some(m.revision())),
        ],
        needs_review: false,
    });
    let text = collagene_core::formats::export_with_report(
        &m,
        collagene_core::formats::TextFormat::GenBank,
    )
    .unwrap()
    .text;
    let imported = collagene_core::formats::import_text(&text).unwrap();
    assert_eq!(
        imported.molecules[0].features[0].qualifiers,
        m.features[0].qualifiers
    );
    assert_eq!(
        imported.molecules[0].features[0].location,
        m.features[0].location
    );
}
