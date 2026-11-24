use collagene_model::{Molecule, Sequence, Topology};
use collagene_recombination::*;
fn dna(s: &str) -> String {
    s.split_whitespace().collect::<String>().to_uppercase()
}
fn rc(s: &str) -> String {
    Sequence::dna(s)
        .unwrap()
        .reverse_complement()
        .unwrap()
        .as_str()
        .into()
}
fn input(role: &str, s: &str, t: Topology) -> Input {
    Input {
        role: role.into(),
        molecule: Molecule::dna(role, s, t).unwrap(),
    }
}
fn profile() -> Profile {
    manufacturer_context_profile().unwrap()
}
fn site(p: &Profile, f: Family, n: u8) -> String {
    p.variants
        .iter()
        .find(|v| v.family == f && v.subtype == n)
        .unwrap()
        .sequence
        .clone()
}
fn donor() -> String {
    dna(include_str!("fixtures/pdonr221_seq.txt"))
}
fn entry() -> String {
    dna(include_str!("fixtures/pentr_dtopo_seq.txt"))
}
fn dest() -> String {
    dna(include_str!("fixtures/pdest14_seq.txt"))
}
fn bp_inputs() -> Vec<Input> {
    let p = profile();
    vec![
        input(
            "PCR",
            &format!(
                "{}GATCGTATCGAC{}",
                site(&p, Family::B, 1),
                rc(&site(&p, Family::B, 2))
            ),
            Topology::Linear,
        ),
        input("pDONR221", &donor(), Topology::Circular),
    ]
}
fn rotate(s: &str, k: usize) -> String {
    format!("{}{}", &s[k..], &s[..k])
}
fn circle_eq(a: &str, b: &str) -> bool {
    a.len() == b.len() && format!("{a}{a}").contains(b)
}
fn compare_products(report: &Report, expected: &[(String, Topology)]) {
    assert_eq!(report.candidates.len(), 1, "{:?}", report.diagnostics);
    let products = &report.candidates[0].products;
    for (s, t) in expected {
        assert!(
            products.iter().any(|p| p.molecule.topology == *t
                && if *t == Topology::Circular {
                    circle_eq(p.molecule.sequence.as_str(), s)
                        || circle_eq(p.molecule.sequence.as_str(), &rc(s))
                } else {
                    p.molecule.sequence.as_str() == s
                }),
            "missing oracle product: {} bases {:?}",
            s.len(),
            t
        )
    }
}
#[test]
fn full_profiles_are_literal_sourced_contexts() {
    let p = profile();
    assert_eq!(p.variants.len(), 8);
    assert_eq!(p.sources.len(), 8);
    assert_eq!(
        p.variants
            .iter()
            .map(|v| v.sequence.len())
            .collect::<Vec<_>>(),
        vec![25, 25, 232, 232, 125, 125, 100, 100]
    );
    for v in &p.variants {
        assert_eq!(
            &v.sequence[v.overlap_start..v.overlap_start + 7],
            if v.subtype == 1 { "GTACAAA" } else { "GTACAAG" }
        );
    }
    assert!(p.limitations.iter().any(|x| x.contains("no universal")));
}
#[test]
fn bp_literal_both_products_and_global_length() {
    let inputs = bp_inputs();
    let p = profile();
    let report = recombine(&inputs, &p, Reaction::BP, &Limits::default()).unwrap();
    let pcr = inputs[0].molecule.sequence.as_str();
    let donor = donor();
    let expected_entry = format!(
        "{}{}{}",
        &pcr[8..pcr.len() - 8],
        &donor[2902..],
        &donor[..651]
    );
    let expected_byproduct = format!(
        "{}{}{}",
        &pcr[..8],
        &donor[651..2902],
        &pcr[pcr.len() - 8..]
    );
    compare_products(
        &report,
        &[
            (expected_entry, Topology::Circular),
            (expected_byproduct, Topology::Linear),
        ],
    );
    assert!(report.complete);
    assert!(report.unique_under_declared_profile());
    for candidate in &report.candidates {
        audit_conservation(&inputs, candidate).unwrap();
        for product in &candidate.products {
            audit_product(&inputs, product).unwrap();
            assert_eq!(product.junctions.len(), 2);
            assert!(product.molecule.features.is_empty());
        }
    }
    assert_eq!(
        report.candidates[0]
            .products
            .iter()
            .map(|p| p.molecule.sequence.len())
            .sum::<usize>(),
        4823
    );
}
#[test]
fn lr_literal_source_arms_never_restore_p232() {
    let p = profile();
    let a = entry();
    let b = dest();
    let inputs = vec![
        input("entry", &a, Topology::Circular),
        input("destination", &b, Topology::Circular),
    ];
    let report = recombine(&inputs, &p, Reaction::LR, &Limits::default()).unwrap();
    let oracle0 = format!("{}{}{}", &a[651..721], &b[1904..], &b[..74]);
    let oracle1 = format!("{}{}{}", &b[74..1904], &a[721..], &a[..651]);
    compare_products(
        &report,
        &[(oracle0, Topology::Circular), (oracle1, Topology::Circular)],
    );
    let products = &report.candidates[0].products;
    assert_eq!(
        products
            .iter()
            .map(|p| p.molecule.sequence.len())
            .sum::<usize>(),
        9002
    );
    let p_junctions = products
        .iter()
        .flat_map(|p| &p.junctions)
        .filter(|j| j.family == Family::P)
        .collect::<Vec<_>>();
    assert_eq!(p_junctions.len(), 2);
    assert!(
        p_junctions
            .iter()
            .all(|j| j.matched_output_variants.is_empty())
    );
    assert!(
        report
            .diagnostics
            .iter()
            .any(|d| d.contains("not been restored"))
    );
    assert!(
        products
            .iter()
            .flat_map(|p| &p.junctions)
            .filter(|j| j.family == Family::B)
            .all(|j| !j.matched_output_variants.is_empty())
    );
}
#[test]
fn circular_origin_crossing_every_site_offset() {
    let original = bp_inputs();
    let p = profile();
    let baseline = recombine(&original, &p, Reaction::BP, &Limits::default()).unwrap();
    let expected = baseline.candidates[0]
        .products
        .iter()
        .map(|p| (p.molecule.sequence.as_str().to_owned(), p.molecule.topology))
        .collect::<Vec<_>>();
    for k in [
        0, 1, 569, 570, 650, 651, 700, 800, 2752, 2800, 2901, 2902, 2983, 4760,
    ] {
        let mut changed = original.clone();
        changed[1].molecule.sequence = Sequence::dna(&rotate(&donor(), k)).unwrap();
        let r = recombine(&changed, &p, Reaction::BP, &Limits::default()).unwrap();
        compare_products(&r, &expected);
        for c in &r.candidates {
            audit_conservation(&changed, c).unwrap();
        }
    }
}
#[test]
fn reverse_complement_sources_and_role_order() {
    let original = bp_inputs();
    let p = profile();
    let baseline = recombine(&original, &p, Reaction::BP, &Limits::default()).unwrap();
    for flipped in 0..4 {
        let mut changed = original.clone();
        for (i, input) in changed.iter_mut().enumerate() {
            if flipped & (1 << i) != 0 {
                input.molecule.sequence = input.molecule.sequence.reverse_complement().unwrap();
            }
        }
        let report = recombine(&changed, &p, Reaction::BP, &Limits::default()).unwrap();
        assert_eq!(report.candidates.len(), 1);
        for original_product in &baseline.candidates[0].products {
            let target = &original_product.molecule;
            let checked = check_target(
                &changed,
                &p,
                Reaction::BP,
                &Limits::default(),
                target,
                TargetPolicy {
                    allow_reverse_complement: true,
                },
            )
            .unwrap();
            assert!(!checked.matches.is_empty(), "flipped {flipped}");
        }
    }
    let mut swapped = original.clone();
    swapped.reverse();
    let report = recombine(&swapped, &p, Reaction::BP, &Limits::default()).unwrap();
    assert_eq!(report.candidates.len(), 1);
    for target in &baseline.candidates[0].products {
        assert!(
            !check_target(
                &swapped,
                &p,
                Reaction::BP,
                &Limits::default(),
                &target.molecule,
                TargetPolicy {
                    allow_reverse_complement: true
                }
            )
            .unwrap()
            .matches
            .is_empty()
        );
    }
}
#[test]
fn both_linear_and_circular_linear_topology() {
    let p = profile();
    let a = format!(
        "GG{}CAGTCG{}TT",
        site(&p, Family::B, 1),
        rc(&site(&p, Family::B, 2))
    );
    let b = format!(
        "AA{}GGACCA{}CC",
        site(&p, Family::P, 1),
        rc(&site(&p, Family::P, 2))
    );
    for ta in [Topology::Linear, Topology::Circular] {
        for tb in [Topology::Linear, Topology::Circular] {
            let inputs = vec![input("a", &a, ta), input("b", &b, tb)];
            let r = recombine(&inputs, &p, Reaction::BP, &Limits::default()).unwrap();
            assert_eq!(r.candidates.len(), 1);
            assert_eq!(
                r.candidates[0]
                    .products
                    .iter()
                    .filter(|p| p.molecule.topology == Topology::Circular)
                    .count(),
                usize::from(ta == Topology::Circular) + usize::from(tb == Topology::Circular)
            );
            audit_conservation(&inputs, &r.candidates[0]).unwrap();
        }
    }
}
#[test]
fn reversed_att1_pair_changes_hybrid_families() {
    let p = profile();
    let a = format!(
        "GG{}CAGTCG{}TT",
        rc(&site(&p, Family::B, 1)),
        rc(&site(&p, Family::B, 2))
    );
    let b = format!(
        "AA{}GGACCA{}CC",
        rc(&site(&p, Family::P, 1)),
        rc(&site(&p, Family::P, 2))
    );
    let inputs = vec![
        input("a", &a, Topology::Circular),
        input("b", &b, Topology::Circular),
    ];
    let r = recombine(&inputs, &p, Reaction::BP, &Limits::default()).unwrap();
    assert_eq!(r.candidates.len(), 1);
    let mut families = r.candidates[0]
        .products
        .iter()
        .map(|p| {
            p.junctions
                .iter()
                .map(|j| (j.subtype, j.family))
                .collect::<Vec<_>>()
        })
        .collect::<Vec<_>>();
    for f in &mut families {
        f.sort_by_key(|x| x.0)
    }
    assert!(families.contains(&vec![(1, Family::R), (2, Family::L)]));
    assert!(families.contains(&vec![(1, Family::L), (2, Family::R)]));
}
#[test]
fn incompatible_mutated_and_short_sites_conservative() {
    let p = profile();
    let mut inputs = bp_inputs();
    let original = inputs[0].molecule.sequence.as_str().to_owned();
    let replaced = original.replacen(&site(&p, Family::B, 2), &site(&p, Family::B, 1), 1);
    assert_eq!(replaced, original);
    let mut bytes = original.into_bytes();
    bytes[0] = b'T';
    inputs[0].molecule.sequence = Sequence::dna(std::str::from_utf8(&bytes).unwrap()).unwrap();
    let r = recombine(&inputs, &p, Reaction::BP, &Limits::default()).unwrap();
    assert!(r.candidates.is_empty());
    assert!(r.complete);
    assert!(
        r.diagnostics
            .iter()
            .any(|d| d.contains("does not establish"))
    );
    inputs[0] = input(
        "PCR",
        &format!(
            "{}ATG{}",
            site(&p, Family::B, 1),
            rc(&site(&p, Family::B, 1))
        ),
        Topology::Linear,
    );
    assert!(
        recombine(&inputs, &p, Reaction::BP, &Limits::default())
            .unwrap()
            .candidates
            .is_empty()
    );
    inputs[1] = input("donor-core-only", "GTACAAAGTACAAG", Topology::Circular);
    assert!(
        recombine(&inputs, &p, Reaction::BP, &Limits::default())
            .unwrap()
            .candidates
            .is_empty()
    );
}
#[test]
fn uncertain_calls_withhold_candidate_and_uniqueness() {
    let p = profile();
    let mut inputs = bp_inputs();
    let mut seq = inputs[0].molecule.sequence.as_str().as_bytes().to_vec();
    seq[3] = b'N';
    inputs[0].molecule.sequence = Sequence::dna(std::str::from_utf8(&seq).unwrap()).unwrap();
    let r = recombine(&inputs, &p, Reaction::BP, &Limits::default()).unwrap();
    assert!(r.sites.iter().any(|s| !s.certain));
    assert!(r.candidates.is_empty());
    assert!(!r.unique_under_declared_profile());
}
#[test]
fn duplicate_alternatives_and_caps() {
    let p = profile();
    let mut inputs = bp_inputs();
    let a = inputs[0].molecule.sequence.as_str();
    inputs[0] = input("PCR", &format!("{a}GCGCGC{a}"), Topology::Linear);
    let r = recombine(&inputs, &p, Reaction::BP, &Limits::default()).unwrap();
    assert!(r.candidates.len() > 1);
    assert!(!r.unique_under_declared_profile());
    let mut limits = Limits::default();
    limits.max_candidates = 1;
    let r = recombine(&inputs, &p, Reaction::BP, &limits).unwrap();
    assert!(!r.complete);
    assert_eq!(r.candidates.len(), 1);
    limits = Limits::default();
    limits.max_comparisons = 20;
    let r = recombine(&inputs, &p, Reaction::BP, &limits).unwrap();
    assert!(!r.complete);
    assert!(r.candidates.is_empty());
    limits = Limits::default();
    limits.max_sites = 1;
    assert!(
        !recombine(&inputs, &p, Reaction::BP, &limits)
            .unwrap()
            .complete
    );
    limits = Limits::default();
    limits.max_combinations = 1;
    assert!(
        !recombine(&inputs, &p, Reaction::BP, &limits)
            .unwrap()
            .complete
    );
    limits = Limits::default();
    limits.max_product_bases = 10;
    assert!(
        !recombine(&inputs, &p, Reaction::BP, &limits)
            .unwrap()
            .complete
    );
    limits = Limits::default();
    limits.max_report_bytes = 1024;
    assert!(recombine(&inputs, &p, Reaction::BP, &limits).is_err());
}
#[test]
fn independent_target_mutation_topology_origin_and_policy() {
    let p = profile();
    let inputs = bp_inputs();
    let r = recombine(&inputs, &p, Reaction::BP, &Limits::default()).unwrap();
    let original = &r.candidates[0]
        .products
        .iter()
        .find(|p| p.molecule.topology == Topology::Circular)
        .unwrap()
        .molecule;
    let mut target = original.clone();
    target.sequence = Sequence::dna(&rotate(original.sequence.as_str(), 17)).unwrap();
    let checked = check_target(
        &inputs,
        &p,
        Reaction::BP,
        &Limits::default(),
        &target,
        TargetPolicy {
            allow_reverse_complement: false,
        },
    )
    .unwrap();
    assert_eq!(checked.matches.len(), 1);
    assert_eq!(checked.matches[0].target_origin, 17);
    target.topology = Topology::Linear;
    assert!(
        check_target(
            &inputs,
            &p,
            Reaction::BP,
            &Limits::default(),
            &target,
            TargetPolicy {
                allow_reverse_complement: true
            }
        )
        .unwrap()
        .matches
        .is_empty()
    );
    target.topology = Topology::Circular;
    target.sequence = target.sequence.reverse_complement().unwrap();
    assert!(
        check_target(
            &inputs,
            &p,
            Reaction::BP,
            &Limits::default(),
            &target,
            TargetPolicy {
                allow_reverse_complement: false
            }
        )
        .unwrap()
        .matches
        .is_empty()
    );
    assert!(
        !check_target(
            &inputs,
            &p,
            Reaction::BP,
            &Limits::default(),
            &target,
            TargetPolicy {
                allow_reverse_complement: true
            }
        )
        .unwrap()
        .matches
        .is_empty()
    );
    let mut seq = target.sequence.as_str().as_bytes().to_vec();
    seq[14] = if seq[14] == b'A' { b'C' } else { b'A' };
    target.sequence = Sequence::dna(std::str::from_utf8(&seq).unwrap()).unwrap();
    assert!(
        check_target(
            &inputs,
            &p,
            Reaction::BP,
            &Limits::default(),
            &target,
            TargetPolicy {
                allow_reverse_complement: true
            }
        )
        .unwrap()
        .matches
        .is_empty()
    );
}
#[test]
fn imported_report_replay_and_exact_source_bindings() {
    let p = profile();
    let inputs = bp_inputs();
    let limits = Limits::default();
    let r = recombine(&inputs, &p, Reaction::BP, &limits).unwrap();
    replay(&inputs, &p, Reaction::BP, &limits, &r).unwrap();
    let mut changed = inputs.clone();
    changed[0].role = "different".into();
    assert!(
        r.validate_bindings_for(&changed, &p, Reaction::BP, &limits)
            .is_err()
    );
    let mut forged = r.clone();
    forged.candidates[0].products[0].molecule.name = "forged".into();
    forged
        .validate_bindings_for(&inputs, &p, Reaction::BP, &limits)
        .unwrap();
    assert!(replay(&inputs, &p, Reaction::BP, &limits, &forged).is_err());
    let mut forged = r.clone();
    forged.candidates[0].products[0].segments[0].source_end -= 1;
    assert!(audit_product(&inputs, &forged.candidates[0].products[0]).is_err());
}

#[test]
fn exhaustive_orientation_topology_scope_and_reference_maps() {
    let p = profile();
    for flags in 0..16 {
        let f = |family, subtype, bit| {
            let s = site(&p, family, subtype);
            if flags & (1 << bit) != 0 { rc(&s) } else { s }
        };
        let a = format!("CGG{}GATACAG{}TCC", f(Family::B, 1, 0), f(Family::B, 2, 1));
        let b = format!("CAA{}TAGGATC{}GGC", f(Family::P, 1, 2), f(Family::P, 2, 3));
        let coherent =
            ((flags & 1 != 0) != (flags & 4 != 0)) == ((flags & 2 != 0) != (flags & 8 != 0));
        for ta in [Topology::Linear, Topology::Circular] {
            for tb in [Topology::Linear, Topology::Circular] {
                let inputs = vec![input("a", &a, ta), input("b", &b, tb)];
                let r = recombine(&inputs, &p, Reaction::BP, &Limits::default()).unwrap();
                assert_eq!(
                    r.candidates.len(),
                    usize::from(coherent),
                    "flags {flags}, {ta:?}/{tb:?}"
                );
                assert!(r.complete);
                if coherent {
                    for c in &r.candidates {
                        audit_conservation(&inputs, c).unwrap();
                        for product in &c.products {
                            audit_product(&inputs, product).unwrap();
                            assert!(!product.junctions.is_empty());
                        }
                    }
                } else {
                    assert!(
                        r.diagnostics
                            .iter()
                            .any(|x| x.starts_with("Unsupported orientation"))
                    );
                    assert!(!r.unique_under_declared_profile());
                }
            }
        }
    }
}
#[test]
fn overlapping_sites_duplicate_roles_bad_profiles_and_resource_inputs_rejected() {
    let p = profile();
    let inputs = bp_inputs();
    let limits = Limits::default();
    let mut bad = inputs.clone();
    bad[1].role = bad[0].role.clone();
    assert!(recombine(&bad, &p, Reaction::BP, &limits).is_err());
    let mut bad_profile = p.clone();
    bad_profile.variants[0].overlap_start = usize::MAX;
    assert!(recombine(&inputs, &bad_profile, Reaction::BP, &limits).is_err());
    let mut bad_profile = p.clone();
    bad_profile.variants[0].sequence = "N".repeat(25);
    assert!(recombine(&inputs, &bad_profile, Reaction::BP, &limits).is_err());
    let mut bad_profile = p.clone();
    bad_profile.variants[0].subtype = 3;
    assert!(bad_profile.validate().is_err());
    let mut bad = inputs.clone();
    bad[0]
        .molecule
        .metadata
        .insert("oversized".into(), "x".repeat(17 * 1024 * 1024));
    assert!(batch_revision(&bad, &p, Reaction::BP, &limits).is_err());
    let mut bad_profile = p.clone();
    bad_profile.limitations.push("x".repeat(140000));
    assert!(bad_profile.validate().is_err());
    let mut bad_limits = limits.clone();
    bad_limits.max_candidates = usize::MAX;
    assert!(recombine(&inputs, &p, Reaction::BP, &bad_limits).is_err());
    let r = recombine(&inputs, &p, Reaction::BP, &limits).unwrap();
    let mut malformed = r.candidates[0].products[0].clone();
    malformed.segments[0].output_end = 0;
    malformed.segments[0].output_start = 2;
    assert!(audit_product(&inputs, &malformed).is_err());
}

#[test]
fn both_bp_hybrid_contexts_match_and_lr_origins_are_independent() {
    let p = profile();
    let bp = recombine(&bp_inputs(), &p, Reaction::BP, &Limits::default()).unwrap();
    for product in &bp.candidates[0].products {
        assert!(
            product
                .junctions
                .iter()
                .all(|j| !j.matched_output_variants.is_empty())
        );
        let families = product
            .junctions
            .iter()
            .map(|j| j.family)
            .collect::<Vec<_>>();
        assert!(families == vec![Family::L, Family::L] || families == vec![Family::R, Family::R]);
    }
    let a = entry();
    let b = dest();
    let original = vec![
        input("entry", &a, Topology::Circular),
        input("destination", &b, Topology::Circular),
    ];
    let baseline = recombine(&original, &p, Reaction::LR, &Limits::default()).unwrap();
    let expected = baseline.candidates[0]
        .products
        .iter()
        .map(|p| (p.molecule.sequence.as_str().to_owned(), p.molecule.topology))
        .collect::<Vec<_>>();
    for (ka, kb) in [
        (568, 66),
        (600, 100),
        (651, 74),
        (720, 1903),
        (721, 1904),
        (803, 1911),
        (2579, 6421),
    ] {
        let inputs = vec![
            input("entry", &rotate(&a, ka), Topology::Circular),
            input("destination", &rotate(&b, kb), Topology::Circular),
        ];
        let r = recombine(&inputs, &p, Reaction::LR, &Limits::default()).unwrap();
        compare_products(&r, &expected);
        r.validate().unwrap();
        for c in &r.candidates {
            audit_conservation(&inputs, c).unwrap();
        }
    }
}

#[test]
fn three_components_literal_linear_end_products_and_circle_targets() {
    let p = profile();
    let a = format!(
        "CCC{}GATACAG{}GGG",
        rc(&site(&p, Family::B, 1)),
        rc(&site(&p, Family::B, 2))
    );
    let b = format!(
        "TTA{}CGGATTC{}AAC",
        site(&p, Family::P, 1),
        site(&p, Family::P, 2)
    );
    let inputs = vec![
        input("a", &a, Topology::Linear),
        input("b", &b, Topology::Linear),
    ];
    // Two antiparallel exchanges preserve four free termini: two end joins are
    // linear products, while the two internal intervals form their own circle.
    // Independent literal boundaries: A reverse B1/B2 at3+17 and35+17;
    // B forward P1/P2 at3+82 and242+82.
    let a1 = 20;
    let a2 = 52;
    let b1 = 85;
    let b2 = 324;
    let expected_left = format!("{}{}", &a[..a1], rc(&b[..b1]));
    let expected_right = format!("{}{}", rc(&a[a2..]), &b[b2..]);
    let expected_circle = format!("{}{}", &a[a1..a2], rc(&b[b1..b2]));
    let report = recombine(&inputs, &p, Reaction::BP, &Limits::default()).unwrap();
    assert!(report.complete);
    assert_eq!(report.candidates.len(), 1);
    let c = &report.candidates[0];
    assert_eq!(c.products.len(), 3);
    audit_conservation(&inputs, c).unwrap();
    let expected = [
        (expected_left, Topology::Linear),
        (expected_right, Topology::Linear),
        (expected_circle, Topology::Circular),
    ];
    for (s, t) in expected {
        assert!(
            c.products
                .iter()
                .any(|product| product.molecule.topology == t
                    && if t == Topology::Circular {
                        circle_eq(product.molecule.sequence.as_str(), &s)
                            || circle_eq(product.molecule.sequence.as_str(), &rc(&s))
                    } else {
                        product.molecule.sequence.as_str() == s
                            || product.molecule.sequence.as_str() == rc(&s)
                    })
        );
        let target = Molecule::dna("independent required component", &s, t).unwrap();
        assert!(
            !check_target(
                &inputs,
                &p,
                Reaction::BP,
                &Limits::default(),
                &target,
                TargetPolicy {
                    allow_reverse_complement: true
                }
            )
            .unwrap()
            .matches
            .is_empty()
        );
    }
    assert_eq!(
        c.products.iter().map(|p| p.junctions.len()).sum::<usize>(),
        4
    );
    let mut limits = Limits::default();
    limits.max_components_per_candidate = 2;
    let bounded = recombine(&inputs, &p, Reaction::BP, &limits).unwrap();
    assert!(!bounded.complete);
    assert!(bounded.candidates.is_empty());
    assert!(
        bounded
            .diagnostics
            .iter()
            .any(|d| d.starts_with("Component limit"))
    );
    report.validate().unwrap();
}

#[test]
fn component_planning_precedes_dna_and_withheld_attempts_stay_charged() {
    let p = profile();
    let a = format!(
        "CCC{}GATACAG{}GGG",
        rc(&site(&p, Family::B, 1)).repeat(3),
        rc(&site(&p, Family::B, 2)).repeat(3)
    );
    let b = format!(
        "TTA{}CGGATTC{}AAC",
        site(&p, Family::P, 1).repeat(3),
        site(&p, Family::P, 2).repeat(3)
    );
    let inputs = vec![
        input("a", &a, Topology::Linear),
        input("b", &b, Topology::Linear),
    ];
    let mut limits = Limits::default();
    limits.max_components_per_candidate = 2;
    limits.max_combinations = 11;
    limits.max_product_bases = 1;
    let report = recombine(&inputs, &p, Reaction::BP, &limits).unwrap();
    assert!(!report.complete);
    assert_eq!(report.combinations, 11);
    assert_eq!(report.graph_plans, 11);
    assert_eq!(report.components_withheld, 11);
    assert_eq!(report.materializations, 0);
    assert_eq!(report.materialized_bases, 0);
    assert!(report.candidates.is_empty());
    report.validate().unwrap();
    limits = Limits::default();
    let total = a.len() + b.len();
    limits.max_product_bases = total * 2;
    limits.max_candidates = 64;
    let report = recombine(&inputs, &p, Reaction::BP, &limits).unwrap();
    assert!(!report.complete);
    assert_eq!(report.materializations, 2);
    assert_eq!(report.materialized_bases, total * 2);
    assert_eq!(report.candidates.len(), 2);
    // A candidate rejected later by serialized-report capacity still consumes its
    // full attempted-copy budget. This counter is independent of installed products.
    let simple = bp_inputs();
    let full = recombine(&simple, &p, Reaction::BP, &Limits::default()).unwrap();
    let candidate_bytes = serde_json::to_vec(&full.candidates[0]).unwrap().len();
    let full_bytes = serde_json::to_vec(&full).unwrap().len();
    let mut small = Limits::default();
    small.max_report_bytes = full_bytes - candidate_bytes / 2;
    let rejected = recombine(&simple, &p, Reaction::BP, &small).unwrap();
    assert!(!rejected.complete);
    assert!(rejected.candidates.is_empty());
    assert_eq!(rejected.materializations, 1);
    assert_eq!(
        rejected.materialized_bases,
        simple
            .iter()
            .map(|i| i.molecule.sequence.len())
            .sum::<usize>()
    );
    assert!(
        rejected
            .diagnostics
            .iter()
            .any(|d| d.starts_with("Report-byte"))
    );
    rejected.validate().unwrap();
}
#[test]
fn borrowed_native_models_reject_wide_deep_and_aggregate_trees_before_validation() {
    use collagene_model::{Error, Feature, Location};
    let p = profile();
    let inputs = bp_inputs();
    let limits = Limits::default();
    let feature = |location| Feature {
        kind: "misc_feature".into(),
        location,
        qualifiers: vec![],
        needs_review: false,
    };
    let mut wide = inputs.clone();
    wide[0]
        .molecule
        .features
        .push(feature(Location::Join(vec![Location::range(0, 1); 16385])));
    assert!(matches!(
        recombine(&wide, &p, Reaction::BP, &limits),
        Err(Error::TooLarge)
    ));
    let mut deep = inputs.clone();
    let mut tree = Location::range(0, 1);
    for _ in 0..64 {
        tree = Location::Complement(Box::new(tree));
    }
    deep[0].molecule.features.push(feature(tree));
    assert!(matches!(
        batch_revision(&deep, &p, Reaction::BP, &limits),
        Err(Error::TooLarge)
    ));
    let mut target = Molecule::dna("requirement", "ACGT", Topology::Linear).unwrap();
    target.features = deep[0].molecule.features.clone();
    assert!(matches!(
        check_target(
            &inputs,
            &p,
            Reaction::BP,
            &limits,
            &target,
            TargetPolicy {
                allow_reverse_complement: false
            }
        ),
        Err(Error::TooLarge)
    ));
    target.features = vec![feature(Location::Order(vec![Location::range(0, 1); 16385]))];
    assert!(matches!(
        check_target(
            &inputs,
            &p,
            Reaction::BP,
            &limits,
            &target,
            TargetPolicy {
                allow_reverse_complement: false
            }
        ),
        Err(Error::TooLarge)
    ));
    let mut aggregate = inputs.clone();
    for i in &mut aggregate {
        i.molecule.features = vec![feature(Location::range(0, 1)); 2000];
    }
    target.features = vec![feature(Location::range(0, 1)); 1000];
    assert!(matches!(
        check_target(
            &aggregate,
            &p,
            Reaction::BP,
            &limits,
            &target,
            TargetPolicy {
                allow_reverse_complement: false
            }
        ),
        Err(Error::TooLarge)
    ));
    // The invalid coordinate would fail recursive model.validate with Coordinates;
    // the excess aggregate node budget must be rejected first without descending it.
    let mut aggregate = inputs.clone();
    aggregate[0].molecule.features = vec![
        feature(Location::Join(vec![Location::range(0, 1); 8192])),
        feature(Location::Order(vec![Location::range(0, 1); 8192])),
        feature(Location::range(99999, 100000)),
    ];
    assert!(matches!(
        recombine(&aggregate, &p, Reaction::BP, &limits),
        Err(Error::TooLarge)
    ));
    let mut qualifiers = inputs.clone();
    qualifiers[0].molecule.features.push(Feature {
        kind: "misc_feature".into(),
        location: Location::range(0, 1),
        qualifiers: vec![("note".into(), Some("x".into())); 16385],
        needs_review: false,
    });
    assert!(matches!(
        recombine(&qualifiers, &p, Reaction::BP, &limits),
        Err(Error::TooLarge)
    ));
    let mut strings = inputs.clone();
    strings[0].role = " ".repeat(1_000_000);
    assert!(batch_revision(&strings, &p, Reaction::BP, &limits).is_err());
    let mut imported = recombine(&inputs, &p, Reaction::BP, &limits).unwrap();
    imported.candidates[0].products[0].molecule.features = deep[0].molecule.features.clone();
    assert!(matches!(imported.validate(), Err(Error::TooLarge)));
}
#[test]
fn mixed_family_parents_have_literal_reciprocal_components() {
    let p = profile();
    let a = format!(
        "{}GATACAG{}",
        site(&p, Family::B, 1),
        rc(&site(&p, Family::P, 2))
    );
    let b = format!(
        "{}CGGATTC{}",
        site(&p, Family::P, 1),
        rc(&site(&p, Family::B, 2))
    );
    let inputs = vec![
        input("mixed B1 P2", &a, Topology::Linear),
        input("mixed P1 B2", &b, Topology::Linear),
    ];
    // Independent source-arm arithmetic, without invoking the planned graph.
    let expected0 = format!("{}{}{}", &a[..8], &b[82..256], &a[182..]);
    let expected1 = format!("{}{}{}", &b[..82], &a[8..182], &b[256..]);
    let report = recombine(&inputs, &p, Reaction::BP, &Limits::default()).unwrap();
    compare_products(
        &report,
        &[(expected0, Topology::Linear), (expected1, Topology::Linear)],
    );
    assert!(report.complete);
    assert_eq!(report.unsupported_combinations, 0);
    audit_conservation(&inputs, &report.candidates[0]).unwrap();
    for product in &report.candidates[0].products {
        let families = product
            .junctions
            .iter()
            .map(|j| j.family)
            .collect::<Vec<_>>();
        assert!(families.contains(&Family::L) && families.contains(&Family::R));
    }
}
#[test]
fn target_completeness_preserves_possible_sites_unmodeled_geometry_and_limits() {
    let p = profile();
    let original = bp_inputs();
    let limits = Limits::default();
    let target =
        Molecule::dna("independent requirement", "ATGCGCTAGCTAG", Topology::Linear).unwrap();
    let policy = TargetPolicy {
        allow_reverse_complement: false,
    };
    let mut possible = original.clone();
    let mut sequence = possible[0].molecule.sequence.as_str().as_bytes().to_vec();
    sequence[3] = b'N';
    possible[0].molecule.sequence = Sequence::dna(std::str::from_utf8(&sequence).unwrap()).unwrap();
    let checked = check_target(&possible, &p, Reaction::BP, &limits, &target, policy).unwrap();
    assert!(checked.matches.is_empty());
    assert!(checked.search_complete);
    assert!(checked.possible_sites);
    assert!(checked.ambiguous_source_dna);
    assert_eq!(checked.applicability, TargetApplicability::Unresolved);
    assert!(!checked.no_match_proven_under_model);
    assert!(
        checked
            .diagnostics
            .iter()
            .any(|d| d.contains("possible sites"))
    );
    let a = format!(
        "{}GATACAG{}",
        site(&p, Family::B, 1),
        rc(&site(&p, Family::B, 2))
    );
    let b = format!(
        "{}GATACAG{}",
        rc(&site(&p, Family::P, 1)),
        rc(&site(&p, Family::P, 2))
    );
    let unsupported = vec![
        input("a", &a, Topology::Circular),
        input("b", &b, Topology::Circular),
    ];
    let checked = check_target(&unsupported, &p, Reaction::BP, &limits, &target, policy).unwrap();
    assert!(checked.matches.is_empty());
    assert!(checked.search_complete);
    assert!(checked.unsupported_geometry);
    assert!(!checked.no_match_proven_under_model);
    assert!(
        checked
            .diagnostics
            .iter()
            .any(|d| d.starts_with("Unsupported orientation"))
    );
    let mut bounded = limits.clone();
    bounded.max_comparisons = 1;
    let checked = check_target(&original, &p, Reaction::BP, &bounded, &target, policy).unwrap();
    assert!(!checked.search_complete);
    assert!(!checked.no_match_proven_under_model);
    let checked = check_target(&original, &p, Reaction::BP, &limits, &target, policy).unwrap();
    assert!(checked.search_complete);
    assert!(checked.no_match_proven_under_model);
    assert_eq!(
        checked.applicability,
        TargetApplicability::DeclaredFullSiteProfileOnly
    );
}
#[test]
fn literal_ambiguous_payload_target_match_is_not_fully_specified_dna() {
    let p = profile();
    let mut inputs = bp_inputs();
    let s = inputs[0]
        .molecule
        .sequence
        .as_str()
        .replacen("GATCGTATCGAC", "GATCGNATCGAC", 1);
    inputs[0].molecule.sequence = Sequence::dna(&s).unwrap();
    let report = recombine(&inputs, &p, Reaction::BP, &Limits::default()).unwrap();
    assert_eq!(report.candidates.len(), 1);
    let target = &report.candidates[0]
        .products
        .iter()
        .find(|p| p.molecule.sequence.as_str().contains('N'))
        .unwrap()
        .molecule;
    let check = check_target(
        &inputs,
        &p,
        Reaction::BP,
        &Limits::default(),
        target,
        TargetPolicy {
            allow_reverse_complement: false,
        },
    )
    .unwrap();
    assert!(!check.matches.is_empty());
    assert!(check.matches.iter().all(|m| !m.fully_specified_dna));
    assert!(check.ambiguous_source_dna && check.ambiguous_target_dna);
    assert!(!check.no_match_proven_under_model);
    assert!(check.warnings.iter().any(|w| w.contains("fully specified")));
    let mut only_target_ambiguous = target.clone();
    only_target_ambiguous.sequence = Sequence::dna("NNNNNN").unwrap();
    let original = bp_inputs();
    let check = check_target(
        &original,
        &p,
        Reaction::BP,
        &Limits::default(),
        &only_target_ambiguous,
        TargetPolicy {
            allow_reverse_complement: false,
        },
    )
    .unwrap();
    assert!(check.matches.is_empty());
    assert!(check.ambiguous_target_dna);
    assert!(!check.ambiguous_source_dna);
    assert!(!check.no_match_proven_under_model);
}

#[test]
fn imported_report_shape_is_bounded_before_any_json_recursion() {
    use collagene_model::{Error, Feature, Location};
    let p = profile();
    let inputs = bp_inputs();
    let r = recombine(&inputs, &p, Reaction::BP, &Limits::default()).unwrap();
    let mut forged = r.clone();
    let mut tree = Location::range(0, 1);
    for _ in 0..128 {
        tree = Location::Complement(Box::new(tree));
    }
    forged.candidates[0].products[0]
        .molecule
        .features
        .push(Feature {
            kind: "untrusted".into(),
            location: tree,
            qualifiers: vec![],
            needs_review: false,
        });
    assert!(matches!(forged.validate(), Err(Error::TooLarge)));
    let mut forged = r.clone();
    forged.candidates[0].products[0]
        .molecule
        .metadata
        .insert("claim".into(), "x".repeat(1024 * 1024));
    assert!(matches!(forged.validate(), Err(Error::TooLarge)));
    let mut forged = r.clone();
    forged.candidates[0].products[0].molecule.original_genbank = Some("x".repeat(1024 * 1024));
    assert!(matches!(forged.validate(), Err(Error::TooLarge)));
    let mut forged = r.clone();
    forged.diagnostics = vec!["x".into(); 65569];
    assert!(matches!(forged.validate(), Err(Error::TooLarge)));
    let mut forged = r.clone();
    forged.assumptions.push(" ".repeat(1024 * 1024));
    assert!(matches!(forged.validate(), Err(Error::TooLarge)));
    let mut forged = r.clone();
    forged.candidates[0].products[0].junctions[0].matched_output_variants = vec!["x".into(); 65537];
    assert!(matches!(forged.validate(), Err(Error::TooLarge)));
    let mut forged = r.clone();
    forged.candidates[0].products[0].segments =
        vec![forged.candidates[0].products[0].segments[0].clone(); 8193];
    assert!(forged.validate().is_err());
}
