//! Literal genetic-code, edit-cost, protected-site and boundary-motif oracles.
use collagene_core::sequence::codon_design::{
    self as codon, CodingRange, CodonDomain, CodonWeight, ConstrainedPolicy, DesignOptions,
    DesignStatus, ScopedMotif,
};
use collagene_model::{Feature, Location, Molecule, Topology};
fn dna_source(dna: &str) -> Molecule {
    Molecule::dna("original coding DNA", dna, Topology::Linear).unwrap()
}
fn base_options() -> DesignOptions {
    DesignOptions {
        beam_width: 4096,
        ..Default::default()
    }
}
fn candidate(
    source: &Molecule,
    options: &DesignOptions,
    policy: &ConstrainedPolicy,
) -> codon::ConstrainedDesignReport {
    let report = codon::design_constrained(source, options, policy).unwrap();
    report
        .validate_bindings_for_source(source, options, policy)
        .unwrap();
    assert_eq!(report.design.status, DesignStatus::CandidateFound);
    assert_eq!(report.check.as_ref().unwrap().satisfies, Some(true));
    report
}
#[test]
fn positional_edit_distance_precedes_global_preference_without_changing_legacy_design() {
    let source = dna_source("CTGCTG");
    let mut options = base_options();
    options.preferences = collagene_core::sequence::codon_design::CodonPreferences::Explicit {
        weights: vec![CodonWeight {
            codon: "TTA".into(),
            weight: 100,
        }],
    };
    let legacy = codon::design(&source, &options).unwrap();
    legacy.validate().unwrap();
    assert_eq!(legacy.candidate.unwrap().sequence.as_str(), "TTATTA");
    let report = candidate(&source, &options, &ConstrainedPolicy::default());
    assert_eq!(
        report.design.candidate.as_ref().unwrap().sequence.as_str(),
        "CTGCTG"
    );
    assert_eq!(report.changed_bases, Some(0));
    assert!(report.minimum_changes_proven);
    assert!(report.objective_optimality_proven);
    assert_eq!(report.design.target_protein, "LL");
}
#[test]
fn sparse_domains_and_partial_codon_locks_are_inside_the_search() {
    let source = dna_source("CTGCTG");
    let mut options = base_options();
    options.preferences = collagene_core::sequence::codon_design::CodonPreferences::Explicit {
        weights: vec![
            CodonWeight {
                codon: "TTG".into(),
                weight: 20,
            },
            CodonWeight {
                codon: "CTT".into(),
                weight: 5,
            },
            CodonWeight {
                codon: "TTA".into(),
                weight: 100,
            },
        ],
    };
    let mut policy = ConstrainedPolicy {
        allowed_codons: vec![CodonDomain {
            position: 0,
            codons: vec!["CTT".into(), "TTG".into()],
        }],
        ..Default::default()
    };
    let unrestricted = candidate(&source, &options, &policy);
    assert_eq!(
        unrestricted.design.candidate.unwrap().sequence.as_str(),
        "TTGCTG"
    );
    policy.protected = vec![CodingRange { start: 0, end: 1 }];
    let locked = candidate(&source, &options, &policy);
    assert_eq!(locked.design.candidate.unwrap().sequence.as_str(), "CTTCTG");
    assert_eq!(locked.changed_bases, Some(1));
    let bad = dna_source("TTGCTG");
    let checked = codon::check_constrained(&source, &bad, &options, &policy).unwrap();
    assert_eq!(checked.base.satisfies, Some(true));
    assert_eq!(checked.satisfies, Some(false));
    assert_eq!(checked.protected_violations, vec![0]);
}
#[test]
fn scoped_removal_preserves_a_wanted_identical_site_elsewhere_in_the_same_cds() {
    let source = dna_source("GGTCTCGGTCTC");
    let policy = ConstrainedPolicy {
        protected: vec![CodingRange { start: 6, end: 12 }],
        scoped_forbidden: vec![ScopedMotif {
            motif: "GGTCTC".into(),
            ranges: vec![CodingRange { start: 0, end: 6 }],
        }],
        ..Default::default()
    };
    let report = candidate(&source, &base_options(), &policy);
    assert_eq!(
        report.design.candidate.as_ref().unwrap().sequence.as_str(),
        "GGACTCGGTCTC"
    );
    assert_eq!(report.changed_bases, Some(1));
    assert!(report.minimum_changes_proven);
    let checked = codon::check_constrained(&source, &source, &base_options(), &policy).unwrap();
    assert_eq!(checked.scoped_violations, vec![0]);
    assert_eq!(checked.satisfies, Some(false));
    assert_eq!(
        &report.design.candidate.unwrap().sequence.as_str()[6..],
        "GGTCTC"
    );
}
#[test]
fn both_strands_and_both_fixed_context_boundaries_keep_true_coding_offsets() {
    for (dna, left, right, expected, hit) in [
        ("CTC", "GGT", "", "CTA", -3),
        ("GGT", "", "CTC", "GGA", 0),
        ("GAG", "", "ACC", "GAA", 0),
        ("ACC", "GAG", "", "ACA", -3),
    ] {
        let source = dna_source(dna);
        let options = DesignOptions {
            left_context: left.into(),
            right_context: right.into(),
            ..base_options()
        };
        let policy = ConstrainedPolicy {
            scoped_forbidden: vec![ScopedMotif {
                motif: "GGTCTC".into(),
                ranges: vec![CodingRange { start: 0, end: 3 }],
            }],
            ..Default::default()
        };
        let rejected = codon::check_constrained(&source, &source, &options, &policy).unwrap();
        assert!(rejected.scoped_violations.contains(&hit));
        let report = candidate(&source, &options, &policy);
        assert_eq!(report.design.candidate.unwrap().sequence.as_str(), expected);
        assert_eq!(report.changed_bases, Some(1));
    }
}
#[test]
fn independent_small_enumeration_checks_minimum_changes_weight_and_lexical_ties() {
    let original = "GGTCTC";
    let source = dna_source(original);
    let options = DesignOptions {
        forbidden_motifs: vec!["GGTCTC".into()],
        preferences: collagene_core::sequence::codon_design::CodonPreferences::Explicit {
            weights: vec![
                CodonWeight {
                    codon: "GGG".into(),
                    weight: 10,
                },
                CodonWeight {
                    codon: "TTG".into(),
                    weight: 100,
                },
            ],
        },
        ..base_options()
    };
    let policy = ConstrainedPolicy {
        allowed_codons: vec![CodonDomain {
            position: 1,
            codons: vec!["CTC".into(), "TTG".into()],
        }],
        protected: vec![CodingRange { start: 3, end: 4 }],
        ..Default::default()
    };
    // Hand-specified Gly/Leu synonymous sets, direct source masks and string motifs.
    let mut oracle = Vec::new();
    for gly in ["GGA", "GGC", "GGG", "GGT"] {
        for leu in ["CTC", "TTG"] {
            let dna = format!("{gly}{leu}");
            if dna.as_bytes()[3] != original.as_bytes()[3]
                || dna.contains("GGTCTC")
                || dna.contains("GAGACC")
            {
                continue;
            }
            let changes = dna
                .bytes()
                .zip(original.bytes())
                .filter(|(a, b)| a != b)
                .count();
            let weight = if gly == "GGG" { 10 } else { 1 } + if leu == "TTG" { 100 } else { 1 };
            oracle.push((changes, std::cmp::Reverse(weight), dna));
        }
    }
    oracle.sort();
    let expected = &oracle[0];
    let report = candidate(&source, &options, &policy);
    assert_eq!(report.changed_bases, Some(expected.0));
    assert_eq!(
        report.design.candidate.unwrap().sequence.as_str(),
        expected.2
    );
    assert_eq!(expected.2, "GGGCTC");
    assert!(report.objective_optimality_proven);
}
#[test]
fn empty_domains_are_infeasible_but_pruned_or_exhausted_searches_never_claim_minimality() {
    let source = dna_source("CTGCTG");
    let empty = ConstrainedPolicy {
        allowed_codons: vec![CodonDomain {
            position: 0,
            codons: vec![],
        }],
        ..Default::default()
    };
    let report = codon::design_constrained(&source, &base_options(), &empty).unwrap();
    assert_eq!(report.design.status, DesignStatus::Unsatisfiable);
    assert!(report.design.search_complete);
    assert!(report.design.candidate.is_none());
    let limited = DesignOptions {
        beam_width: 1,
        ..base_options()
    };
    let report = candidate(&source, &limited, &ConstrainedPolicy::default());
    assert!(!report.design.search_complete);
    assert!(!report.minimum_changes_proven);
    assert!(!report.objective_optimality_proven);
    let tiny = DesignOptions {
        max_work: 1,
        ..base_options()
    };
    let report = codon::design_constrained(&source, &tiny, &ConstrainedPolicy::default()).unwrap();
    assert_eq!(report.design.status, DesignStatus::SearchIncomplete);
    assert!(!report.design.search_complete);
    assert!(report.design.candidate.is_none());
    let checked =
        codon::check_constrained(&source, &source, &tiny, &ConstrainedPolicy::default()).unwrap();
    assert!(!checked.complete);
    assert_eq!(checked.satisfies, None);
    assert_eq!(checked.changed_bases, None);
}
#[test]
fn invalid_source_domains_exceptions_and_stale_policy_context_are_rejected() {
    let source = dna_source("GGTCTC");
    let foreign = ConstrainedPolicy {
        allowed_codons: vec![CodonDomain {
            position: 0,
            codons: vec!["ATG".into()],
        }],
        ..Default::default()
    };
    assert!(codon::constrained_source_preflight(&source, &base_options(), &foreign).is_err());
    let beyond = ConstrainedPolicy {
        protected: vec![CodingRange { start: 0, end: 7 }],
        ..Default::default()
    };
    assert!(codon::design_constrained(&source, &base_options(), &beyond).is_err());
    let mut exceptional = source.clone();
    exceptional.features.push(Feature {
        kind: "CDS".into(),
        location: Location::range(0, 6),
        qualifiers: vec![("transl_except".into(), Some("(pos:1..3,aa:Sec)".into()))],
        needs_review: false,
    });
    assert!(
        codon::design_constrained(&exceptional, &base_options(), &ConstrainedPolicy::default())
            .is_err()
    );
    let report = candidate(&source, &base_options(), &ConstrainedPolicy::default());
    let locked = ConstrainedPolicy {
        protected: vec![CodingRange { start: 0, end: 3 }],
        ..Default::default()
    };
    assert!(
        report
            .validate_bindings_for_source(&source, &base_options(), &locked)
            .is_err()
    );
    let mut tampered = report.clone();
    tampered.changed_bases = Some(1);
    assert!(tampered.validate_binding(&report.context_revision).is_err());
}
#[test]
fn genetic_code_is_explicit_and_initiator_interpretation_is_not_silently_added() {
    let source = dna_source("TGA");
    let options = DesignOptions {
        genetic_code: 2,
        forbidden_motifs: vec!["TGA".into()],
        ..base_options()
    };
    let report = candidate(&source, &options, &ConstrainedPolicy::default());
    assert_eq!(report.design.target_protein, "W");
    assert_eq!(report.design.candidate.unwrap().sequence.as_str(), "TGG");
    let source = dna_source("TTG");
    let report = candidate(
        &source,
        &DesignOptions {
            genetic_code: 11,
            ..base_options()
        },
        &ConstrainedPolicy::default(),
    );
    assert_eq!(report.design.target_protein, "L");
    assert_eq!(report.design.candidate.unwrap().sequence.as_str(), "TTG");
}
#[test]
fn successful_total_accounts_for_search_and_the_single_independent_checker() {
    let source = dna_source("ATGCTG");
    let options = base_options();
    let policy = ConstrainedPolicy::default();
    let report = candidate(&source, &options, &policy);
    let checked = report.check.as_ref().unwrap();
    assert!(checked.work_used > 0);
    assert_eq!(
        report.total_work_used,
        report.design.work_used + checked.work_used
    );
    assert!(report.total_work_used <= options.max_work * 2);
    let mut forged = report.clone();
    forged.total_work_used = forged.design.work_used;
    assert!(forged.validate_binding(&report.context_revision).is_err());
}
#[test]
fn failed_checker_cost_is_retained_when_its_candidate_and_receipt_are_withheld() {
    // Context-only scoped windows cannot match the coding range. The complete
    // checker still visits both full strands, while incremental search only
    // checks newly completed coding windows; its independent cost is larger.
    let source = dna_source("ATG");
    let options = DesignOptions {
        left_context: "AC".repeat(1024),
        max_work: 8000,
        ..base_options()
    };
    let policy = ConstrainedPolicy {
        scoped_forbidden: vec![ScopedMotif {
            motif: "AAAAAA".into(),
            ranges: vec![CodingRange { start: 0, end: 3 }],
        }],
        ..Default::default()
    };
    let report = codon::design_constrained(&source, &options, &policy).unwrap();
    report.validate_binding(&report.context_revision).unwrap();
    assert_eq!(report.design.status, DesignStatus::SearchIncomplete);
    assert!(report.design.candidate.is_none());
    assert!(report.check.is_none());
    assert!(
        report.total_work_used > report.design.work_used,
        "independent checker really ran even though its receipt was withheld"
    );
    assert!(report.total_work_used <= options.max_work * 2);
    let mut forged = report.clone();
    forged.total_work_used = options.max_work * 2 + 1;
    assert!(forged.validate_binding(&report.context_revision).is_err());
}
