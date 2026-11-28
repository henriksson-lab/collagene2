//! Independent oligo constraints and explicit inverse-PCR strand diagrams.
use collagene_core::{
    cloning::pcr::{mutagenesis::*, primers::BindingOptions},
    sequence::analysis::{self, Primer},
    verification::evidence::Status,
};
use collagene_model::{Molecule, Topology};
use std::collections::BTreeSet;
fn dna(s: &str, t: Topology) -> Molecule {
    Molecule::dna("mutation fixture", s, t).unwrap()
}
fn p(s: &str, tail: &str) -> Primer {
    Primer::new("oligo", s, tail).unwrap()
}
fn opts(mismatches: usize, exact: usize) -> PcrOptions {
    PcrOptions {
        binding: BindingOptions {
            max_mismatches: mismatches,
            min_three_prime_exact: exact,
            ..Default::default()
        },
        ..Default::default()
    }
}
#[test]
fn primer_encoded_substitutions_replace_both_template_ends_and_tails_stay_distinct() {
    let source = dna("ACGACCCCGGGGTAAC", Topology::Linear);
    let report = simulate_pcr(&source, &p("ATGA", "CGT"), &p("GTCA", "GAC"), &opts(1, 1)).unwrap();
    let product = report
        .products
        .iter()
        .find(|p| {
            p.forward_start == 0
                && p.reverse_start == 12
                && p.forward_primer_index == 0
                && p.reverse_primer_index == 1
        })
        .unwrap();
    assert_eq!(product.sequence.as_str(), "CGTATGACCCCGGGGTGACGTC");
    assert_eq!(product.template_length, 16);
    assert_eq!(product.forward_tail_length, 3);
    assert_eq!(product.reverse_tail_length, 3);
    assert_eq!(
        product
            .substitutions
            .iter()
            .map(|c| (
                c.source_position,
                c.product_position,
                c.template_base,
                c.product_base
            ))
            .collect::<Vec<_>>(),
        vec![(1, 4, 'C', 'T'), (13, 16, 'A', 'G')]
    );
    assert_eq!(product.forward_binding.mismatches, vec![1]);
    assert_eq!(product.reverse_binding.mismatches, vec![2]);
    assert_eq!(product.reverse_binding.three_prime_exact, 1);
    assert!(product.certain);
    let strict = simulate_pcr(&source, &p("ATGA", "CGT"), &p("GTCA", "GAC"), &opts(1, 2)).unwrap();
    assert!(!strict.products.iter().any(|p| p.forward_start == 0
        && p.reverse_start == 12
        && p.forward_primer_index == 0
        && p.reverse_primer_index == 1));
}
#[test]
fn consistent_overlap_is_preview_evidence_but_conflicting_bases_are_not_repaired() {
    let source = dna("ACGTAA", Topology::Linear);
    let report = simulate_pcr(&source, &p("ATGT", ""), &p("TTAC", ""), &opts(1, 1)).unwrap();
    let product = report
        .products
        .iter()
        .find(|p| {
            p.forward_start == 0
                && p.reverse_start == 2
                && p.forward_primer_index == 0
                && p.reverse_primer_index == 1
        })
        .unwrap();
    assert_eq!(product.sequence.as_str(), "ATGTAA");
    assert!(product.overlapping_domains);
    assert!(report.search_complete);
    assert!(!report.enumeration_complete);
    assert!(report.unsupported_geometries > 0);
    let conflicting = simulate_pcr(&source, &p("ATGT", ""), &p("TTTC", ""), &opts(1, 1)).unwrap();
    assert!(!conflicting.products.iter().any(|p| p.forward_start == 0
        && p.reverse_start == 2
        && p.forward_primer_index == 0
        && p.reverse_primer_index == 1));
    assert!(conflicting.conflicting_overlaps > 0);
    assert!(!conflicting.enumeration_complete);
}
#[test]
fn circular_origin_geometry_and_mutations_are_recorded_in_source_and_product_coordinates() {
    let source = dna("GGTAACACGACCCCGG", Topology::Circular);
    let report = simulate_pcr(&source, &p("ATGA", "TT"), &p("GTCA", "CG"), &opts(1, 1)).unwrap();
    let product = report
        .products
        .iter()
        .find(|p| {
            p.forward_start == 6
                && p.reverse_start == 2
                && p.forward_primer_index == 0
                && p.reverse_primer_index == 1
        })
        .unwrap();
    assert_eq!(product.sequence.as_str(), "TTATGACCCCGGGGTGACCG");
    assert!(product.crosses_origin);
    assert!(
        product
            .substitutions
            .iter()
            .any(|c| c.source_position == 7 && c.product_position == 3 && c.product_base == 'T')
    );
    assert!(
        product
            .substitutions
            .iter()
            .any(|c| c.source_position == 3 && c.product_base == 'G')
    );
}
#[test]
fn ambiguity_and_all_budget_exhaustions_remain_explicit() {
    let source = dna("ACGNCCCCGGGGTAAC", Topology::Linear);
    let report = simulate_pcr(&source, &p("ACGA", ""), &p("GTTA", ""), &opts(0, 0)).unwrap();
    assert!(report.ambiguous_bindings > 0);
    assert!(!report.enumeration_complete);
    assert!(report.products.iter().any(|p| !p.certain));
    let repetitive = dna("ATATATATAT", Topology::Circular);
    for options in [
        PcrOptions {
            max_pair_checks: 1,
            ..opts(0, 0)
        },
        PcrOptions {
            max_products: 1,
            ..opts(0, 0)
        },
        PcrOptions {
            max_product_bases: 1,
            ..opts(0, 0)
        },
        PcrOptions {
            binding: BindingOptions {
                max_work: 2,
                ..Default::default()
            },
            ..opts(0, 0)
        },
    ] {
        let r = simulate_pcr(&repetitive, &p("AT", ""), &p("TA", ""), &options).unwrap();
        assert!(!r.search_complete);
        assert!(!r.enumeration_complete);
        assert!(r.warnings.iter().any(|w| w.contains("truncated")));
    }
}
fn complement(b: u8) -> u8 {
    match b {
        b'A' => b'T',
        b'T' => b'A',
        b'C' => b'G',
        b'G' => b'C',
        _ => panic!("concrete oracle"),
    }
}
fn reverse(s: &str) -> String {
    s.bytes().rev().map(|b| complement(b) as char).collect()
}
fn eligible(
    reference: &str,
    start: usize,
    primer: &str,
    reverse_role: bool,
    circular: bool,
    max_mm: usize,
    min_end: usize,
) -> bool {
    let n = reference.len();
    let m = primer.len();
    if !circular && start + m > n {
        return false;
    }
    let oriented: Vec<u8> = (0..m)
        .map(|j| {
            if reverse_role {
                complement(reference.as_bytes()[(start + m - 1 - j) % n])
            } else {
                reference.as_bytes()[(start + j) % n]
            }
        })
        .collect();
    let eq: Vec<_> = oriented
        .iter()
        .zip(primer.bytes())
        .map(|(a, b)| *a == b)
        .collect();
    eq.iter().filter(|v| !**v).count() <= max_mm
        && eq.iter().rev().take_while(|v| **v).count() >= min_end
}
fn oracle(
    source: &Molecule,
    oligos: &[Primer; 2],
    mm: usize,
    end: usize,
) -> BTreeSet<(usize, usize, usize, String, String, String, String, String)> {
    let n = source.sequence.len();
    let circular = source.topology == Topology::Circular;
    let s = source.sequence.as_str();
    let mut out = BTreeSet::new();
    // Enumerate all possible reference intervals, not binding-hit products.
    for start in 0..n {
        for length in 1..=n {
            if !circular && start + length > n {
                continue;
            }
            for f in oligos {
                for r in oligos {
                    let (fl, rl) = (f.annealing.len(), r.annealing.len());
                    if length < fl.max(rl) {
                        continue;
                    }
                    let rs = (start + length - rl) % n;
                    if !eligible(s, start, f.annealing.as_str(), false, circular, mm, end)
                        || !eligible(s, rs, r.annealing.as_str(), true, circular, mm, end)
                    {
                        continue;
                    }
                    let rt = reverse(r.annealing.as_str());
                    let mut bases = String::new();
                    let mut conflict = false;
                    for offset in 0..length {
                        let a = (offset < fl).then(|| f.annealing.as_str().as_bytes()[offset]);
                        let b =
                            (offset >= length - rl).then(|| rt.as_bytes()[offset - (length - rl)]);
                        match (a, b) {
                            (Some(a), Some(b)) if a != b => {
                                conflict = true;
                                break;
                            }
                            (Some(a), _) => bases.push(a as char),
                            (_, Some(b)) => bases.push(b as char),
                            _ => bases.push(s.as_bytes()[(start + offset) % n] as char),
                        }
                    }
                    if !conflict {
                        out.insert((
                            start,
                            rs,
                            length,
                            f.annealing.as_str().into(),
                            f.tail.clone(),
                            r.annealing.as_str().into(),
                            r.tail.clone(),
                            format!("{}{bases}{}", f.tail, reverse(&r.tail)),
                        ));
                    }
                }
            }
        }
    }
    out
}
#[test]
fn all_short_binary_templates_and_oligos_match_independent_interval_constraints() {
    let strings = |len: usize| {
        (0..1usize << len)
            .map(|bits| {
                (0..len)
                    .map(|i| if bits & (1 << i) == 0 { 'A' } else { 'T' })
                    .collect::<String>()
            })
            .collect::<Vec<_>>()
    };
    for template in strings(5) {
        for topology in [Topology::Linear, Topology::Circular] {
            let source = dna(&template, topology);
            for f in strings(2) {
                for r in strings(2) {
                    let oligos = [p(&f, "G"), p(&r, "C")];
                    for (mm, exact) in [(0, 0), (1, 0), (1, 1)] {
                        let report =
                            simulate_pcr(&source, &oligos[0], &oligos[1], &opts(mm, exact))
                                .unwrap();
                        let observed: BTreeSet<_> = report
                            .products
                            .iter()
                            .map(|product| {
                                let f = &oligos[product.forward_primer_index];
                                let r = &oligos[product.reverse_primer_index];
                                (
                                    product.forward_start,
                                    product.reverse_start,
                                    product.template_length,
                                    f.annealing.as_str().to_string(),
                                    f.tail.clone(),
                                    r.annealing.as_str().to_string(),
                                    r.tail.clone(),
                                    product.sequence.as_str().to_string(),
                                )
                            })
                            .collect();
                        assert!(report.search_complete);
                        assert_eq!(
                            observed,
                            oracle(&source, &oligos, mm, exact),
                            "{template} {topology:?} {f} {r} {mm} {exact}"
                        );
                    }
                }
            }
        }
    }
}
#[test]
fn exact_nonoverlapping_products_agree_with_the_existing_exact_pcr_definition() {
    let source = dna("AGTCCGATGCTACGTAGGACTTCAGACCTGTA", Topology::Linear);
    let f = p("AGTCCGAT", "AAC");
    let r = p("TACAGGTC", "CGG");
    let strict = simulate_pcr(&source, &f, &r, &opts(0, 1)).unwrap();
    let exact = analysis::simulate_pcr_geometry(&source, &f, &r).unwrap();
    assert!(strict.enumeration_complete);
    assert_eq!(strict.products.len(), 1);
    assert_eq!(exact.len(), 1);
    assert_eq!(strict.products[0].sequence, exact[0].product.sequence);
    assert!(strict.products[0].substitutions.is_empty());
}

#[test]
fn concrete_mismatch_product_can_be_complete_and_single_oligo_roles_are_not_ignored() {
    let source = dna("AGTCCGATGCTACGTAGGACTTCAGACCTGTA", Topology::Linear);
    let report =
        simulate_pcr(&source, &p("AGACCGAT", ""), &p("TCCAGGTC", ""), &opts(1, 2)).unwrap();
    assert!(report.enumeration_complete);
    assert_eq!(report.products.len(), 1);
    assert_eq!(
        report.products[0].sequence.as_str(),
        "AGACCGATGCTACGTAGGACTTCAGACCTGGA"
    );
    assert_eq!(
        report.products[0]
            .substitutions
            .iter()
            .map(|c| c.source_position)
            .collect::<Vec<_>>(),
        vec![2, 30]
    );
    let source = dna("ACGACCCCTCGT", Topology::Linear);
    let report = simulate_pcr(&source, &p("ACGA", "GAC"), &p("AAAA", ""), &opts(0, 1)).unwrap();
    assert!(report.enumeration_complete);
    assert_eq!(report.products.len(), 1);
    let product = &report.products[0];
    assert_eq!(
        (product.forward_primer_index, product.reverse_primer_index),
        (0, 0)
    );
    assert_eq!(product.sequence.as_str(), "GACACGACCCCTCGTGTC");
}
fn source_circle() -> Molecule {
    dna("AGTCCGATGCTACGTAGGACTTCAGACCTGTA", Topology::Circular)
}
fn desired(source: &Molecule, start: usize, end: usize, inserted: &str) -> Molecule {
    let s = source.sequence.as_str();
    let n = s.len();
    let retained: String = (0..n - (end - start))
        .map(|i| s.as_bytes()[(end + i) % n] as char)
        .collect();
    dna(&format!("{inserted}{retained}"), Topology::Circular)
}
#[test]
fn inverse_replacement_substitution_insertion_deletion_and_tail_splits_verify_independently() {
    let source = source_circle();
    for (start, end, inserted) in [
        (3, 5, "GGA"),
        (10, 10, "AACG"),
        (9, 13, ""),
        (28, 35, "CTG"),
    ] {
        let target = desired(&source, start, end, inserted);
        for split in 0..=inserted.len() {
            let report = design_inverse_pcr(
                &source,
                &target,
                None,
                &InversePcrOptions {
                    edit_start: start,
                    edit_end: end,
                    inserted: inserted.into(),
                    forward_length: 8,
                    reverse_length: 8,
                    reverse_tail_bases: split,
                },
            )
            .unwrap();
            assert_eq!(
                report.verification.status,
                Status::ConsistentUnderAssumptions,
                "{start} {end} {inserted} {split}: {:?}",
                report.verification.findings
            );
            assert_eq!(report.source_revision, source.revision());
            assert_eq!(report.target_revision, target.revision());
            assert_eq!(
                report.forward.annealing.as_str(),
                report.forward_location.extract(&source.sequence).unwrap()
            );
            assert_eq!(
                report.reverse.annealing.as_str(),
                report.reverse_location.extract(&source.sequence).unwrap()
            );
            assert_eq!(report.forward.tail, inserted[split..]);
            assert_eq!(report.reverse.tail, reverse(&inserted[..split]));
            let verification = &report.verification;
            assert_eq!(verification.reconstructed_count(), 3);
            assert_eq!(
                verification
                    .product(0)
                    .expect("PCR step must reconstruct")
                    .topology,
                Topology::Linear
            );
            assert_eq!(
                verification
                    .product(2)
                    .expect("final step must reconstruct")
                    .topology,
                Topology::Circular
            );
            assert!(
                analysis::rotation_offset(
                    verification
                        .product(2)
                        .expect("final step must reconstruct")
                        .sequence
                        .as_str(),
                    target.sequence.as_str()
                )
                .is_some()
            );
        }
    }
}
#[test]
fn inverse_design_binds_the_separate_target_and_retains_failed_or_ambiguous_verification() {
    let source = source_circle();
    let options = InversePcrOptions {
        edit_start: 3,
        edit_end: 5,
        inserted: "GGA".into(),
        forward_length: 8,
        reverse_length: 8,
        reverse_tail_bases: 0,
    };
    let wrong = dna("ACGTACGTACGTACGT", Topology::Circular);
    let report = design_inverse_pcr(&source, &wrong, None, &options).unwrap();
    assert_eq!(report.verification.status, Status::Invalid);
    assert_eq!(report.verification.reconstructed_count(), 3);
    let repeated = dna("ACGTACGTACGTACGTACGTACGTACGTACGT", Topology::Circular);
    let target = desired(&repeated, 3, 5, "GGA");
    let report = design_inverse_pcr(&repeated, &target, None, &options).unwrap();
    assert_eq!(report.verification.status, Status::Incomplete);
    assert!(
        report
            .verification
            .findings
            .iter()
            .any(|f| f.message.contains("require disambiguation"))
    );
}
#[test]
fn unsupported_inputs_limits_and_replacement_geometry_fail_explicitly() {
    let source = source_circle();
    let target = desired(&source, 3, 5, "GGA");
    let options = InversePcrOptions {
        edit_start: 3,
        edit_end: 5,
        inserted: "GGA".into(),
        forward_length: 8,
        reverse_length: 8,
        reverse_tail_bases: 0,
    };
    for changed in [
        InversePcrOptions {
            edit_start: 32,
            ..options.clone()
        },
        InversePcrOptions {
            edit_end: 36,
            edit_start: 3,
            ..options.clone()
        },
        InversePcrOptions {
            reverse_tail_bases: 4,
            ..options.clone()
        },
        InversePcrOptions {
            inserted: "N".into(),
            ..options.clone()
        },
        InversePcrOptions {
            forward_length: 80,
            ..options.clone()
        },
    ] {
        assert!(design_inverse_pcr(&source, &target, None, &changed).is_err());
    }
    let mut physical = source.clone();
    physical
        .metadata
        .insert("duplex_model".into(), "untrusted".into());
    assert!(simulate_pcr(&physical, &p("ACGT", ""), &p("ACGT", ""), &opts(0, 0)).is_err());
    assert!(design_inverse_pcr(&physical, &target, None, &options).is_err());
    assert!(simulate_pcr(&source, &p("ACGN", ""), &p("ACGT", ""), &opts(0, 0)).is_err());
    assert!(
        simulate_pcr(
            &source,
            &p("ACGT", ""),
            &p("ACGT", ""),
            &PcrOptions {
                max_pair_checks: 0,
                ..opts(0, 0)
            }
        )
        .is_err()
    );
}
