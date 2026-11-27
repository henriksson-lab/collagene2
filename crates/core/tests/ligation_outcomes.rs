//! Independent intended/alternative sequences. No intended target drives search.
use collagene_core::{
    Sequence, Topology, analysis,
    chemistry::{self, AnnealOptions, EndTreatment},
    ligation_outcomes::{self as l, Limit, Operand, Options, Placement, Request},
    modified_duplex::{self as tracked, TrackedDuplex},
    tailing::TailingBase,
};
fn blunt(text: &str, phosphate: bool) -> TrackedDuplex {
    let p = TrackedDuplex::exact_pcr(&Sequence::dna(text).unwrap()).unwrap();
    if phosphate {
        p.treat(EndTreatment::Phosphorylate).unwrap()
    } else {
        p
    }
}
fn request(circular: bool) -> Request {
    Request {
        intended: vec![
            Placement {
                key: "vector".into(),
                reverse: false,
            },
            Placement {
                key: "insert".into(),
                reverse: false,
            },
        ],
        circular,
        allow_nicked: false,
        options: Options::default(),
    }
}
fn operands<'a>(a: &'a TrackedDuplex, b: &'a TrackedDuplex) -> [Operand<'a>; 2] {
    [
        Operand {
            key: "vector",
            material: a,
        },
        Operand {
            key: "insert",
            material: b,
        },
    ]
}
fn check(report: &l::Report, forward: &str, reverse: &str) {
    assert_eq!(report.intended.sequence.as_str(), forward);
    assert_eq!(report.alternatives.len(), 1);
    assert_eq!(report.alternatives[0].sequence.as_str(), reverse);
    assert!(!report.intended.path[1].reverse);
    assert!(report.alternatives[0].path[1].reverse);
    assert_eq!(report.paths_checked, 2);
    assert!(report.exhaustive_within_scope);
    report.validate_binding(&report.context_revision).unwrap();
}
#[test]
fn ta_tail_ligation_keeps_both_insert_orientations_with_the_two_real_junction_bases() {
    let vector = blunt("AAACCC", false)
        .tail(TailingBase::T)
        .unwrap()
        .treat(EndTreatment::Phosphorylate)
        .unwrap();
    let insert = blunt("GACTGA", false)
        .tail(TailingBase::A)
        .unwrap()
        .treat(EndTreatment::Phosphorylate)
        .unwrap();
    let report = l::analyze(&operands(&vector, &insert), &request(true)).unwrap();
    check(&report, "AAACCCTGACTGAA", "AAACCCTTCAGTCA");
    assert!(report.intended.duplex.nicks.is_empty());
    assert_eq!(report.intended.topology, Topology::Circular);
}
#[test]
fn blunt_ligation_retains_intended_and_reversed_insert_without_probabilities() {
    let vector = blunt("AGTCCG", true);
    let insert = blunt("GACTGA", true);
    let before = (vector.clone(), insert.clone());
    let report = l::analyze(&operands(&vector, &insert), &request(true)).unwrap();
    check(&report, "AGTCCGGACTGA", "AGTCCGTCAGTC");
    assert_eq!((vector, insert), before);
    assert!(
        report
            .warnings
            .iter()
            .any(|s| s.contains("Other copy numbers"))
    );
}
fn digested_circle(text: &str, names: &[&str]) -> Vec<TrackedDuplex> {
    let source = chemistry::Duplex::from_molecule(
        &collagene_core::Molecule::dna("Actual original circle", text, Topology::Circular).unwrap(),
    )
    .unwrap();
    let source = TrackedDuplex::from_unknown_duplex(&source).unwrap();
    let all = analysis::starter_enzymes();
    let enzymes = names
        .iter()
        .map(|n| all.iter().find(|e| e.name == *n).unwrap().clone())
        .collect::<Vec<_>>();
    tracked::digest(&source, &enzymes)
        .unwrap()
        .into_iter()
        .map(|f| f.duplex)
        .collect()
}
#[test]
fn identical_ecori_ends_are_nondirectional_but_distinct_ecori_hindiii_ends_are_directional() {
    let vector = digested_circle("GAATTCAAACCC", &["EcoRI"]).remove(0);
    let insert = digested_circle("GAATTCGGGTTT", &["EcoRI"]).remove(0);
    let report = l::analyze(&operands(&vector, &insert), &request(true)).unwrap();
    check(
        &report,
        "AATTCAAACCCGAATTCGGGTTTG",
        "AATTCAAACCCGAATTCAAACCCG",
    );
    let mut fragments = digested_circle("GAATTCAAACCCAAGCTTGGGTTT", &["EcoRI", "HindIII"]);
    assert_eq!(fragments.len(), 2);
    let insert = fragments.remove(1);
    let vector = fragments.remove(0);
    let report = l::analyze(&operands(&vector, &insert), &request(true)).unwrap();
    assert_eq!(
        report.intended.sequence.as_str(),
        "AATTCAAACCCAAGCTTGGGTTTG"
    );
    assert!(report.alternatives.is_empty());
    assert_eq!(report.incompatible_paths, 1);
    assert!(report.exhaustive_within_scope);
    report.validate_binding(&report.context_revision).unwrap();
}
#[test]
fn all_hydroxyl_blunt_ends_never_create_disconnected_phantom_ligation_even_allowing_nicks() {
    let vector = blunt("AGTCCG", false);
    let insert = blunt("GACTGA", false);
    let mut r = request(true);
    r.allow_nicked = true;
    let native = tracked::ligate(&[vector.clone(), insert.clone()], true, true)
        .unwrap_err()
        .to_string();
    assert!(native.contains("neither strand"), "{native}");
    assert!(l::analyze(&operands(&vector, &insert), &r).is_err());
}
fn annealed(top: &str, bottom: &str, min: usize) -> TrackedDuplex {
    let candidates = chemistry::anneal_oligos(
        &Sequence::dna(top).unwrap(),
        &Sequence::dna(bottom).unwrap(),
        &AnnealOptions {
            min_paired_bases: min,
            ..Default::default()
        },
    )
    .unwrap();
    assert_eq!(candidates.candidates.len(), 1);
    TrackedDuplex::from_unknown_duplex(&candidates.candidates[0].duplex)
        .unwrap()
        .treat(EndTreatment::Phosphorylate)
        .unwrap()
}
#[test]
fn reversing_actual_unequal_thirteen_seventeen_strands_never_invents_or_discards_bases() {
    let vector = annealed("TGCGA", "CTACTCGCA", 5);
    let insert = annealed("AGTCCGATGCTAC", "GTAGGTAGCATCGGACT", 12);
    assert_eq!(insert.duplex.top.len(), 13);
    assert_eq!(insert.duplex.bottom.len(), 17);
    let unknown = tracked::summary(&vector).unwrap().unknown_eligible
        + tracked::summary(&insert).unwrap().unknown_eligible;
    let before = insert.clone();
    let mut r = request(false);
    r.intended[1].reverse = true;
    let report = l::analyze(&operands(&vector, &insert), &r).unwrap();
    assert_eq!(report.intended.sequence.as_str(), "TGCGAGTAGGTAGCATCGGACT");
    assert_eq!(
        report.intended.duplex.bottom.as_str(),
        "AGTCCGATGCTACCTACTCGCA"
    );
    assert_eq!(report.intended.duplex.top.len(), 22);
    assert_eq!(report.intended.duplex.bottom.len(), 22);
    assert_eq!(report.intended.modifications.unknown_eligible, unknown);
    assert_eq!(insert, before);
    report.validate_binding(&report.context_revision).unwrap();
}
#[test]
fn capped_alternatives_keep_the_compatible_intended_product_extractable_and_bound_to_its_path() {
    let vector = blunt("AGTCCG", true);
    let insert = blunt("GACTGA", true);
    let mut r = request(true);
    r.options.max_path_checks = 1;
    let report = l::analyze(&operands(&vector, &insert), &r).unwrap();
    assert_eq!(report.intended.sequence.as_str(), "AGTCCGGACTGA");
    assert!(report.alternatives.is_empty());
    assert!(!report.enumeration_complete);
    assert!(report.limits.contains(&Limit::PathChecks));
    report.validate_binding(&report.context_revision).unwrap();
    let mut tampered = report.clone();
    tampered.intended.path[1].reverse = true;
    assert!(tampered.validate_binding(&report.context_revision).is_err());
    let mut r = request(true);
    r.options.max_products = 1;
    let report = l::analyze(&operands(&vector, &insert), &r).unwrap();
    assert!(report.limits.contains(&Limit::Products));
    assert_eq!(report.intended.sequence.as_str(), "AGTCCGGACTGA");
    let mut r = request(true);
    r.options.max_work = 240; // 4*24 source-base units + 6*24 intended units.
    let report = l::analyze(&operands(&vector, &insert), &r).unwrap();
    assert!(report.limits.contains(&Limit::Work));
    assert_eq!(report.work, 240);
    assert_eq!(report.intended.sequence.as_str(), "AGTCCGGACTGA");
    report.validate_binding(&report.context_revision).unwrap();
}
#[test]
fn unknown_alternative_junction_chemistry_is_not_reported_as_proven_absence() {
    let vector = blunt("AGTCCG", true);
    let mut insert = blunt("GACTGA", true);
    insert.duplex.right = chemistry::Terminus::UNKNOWN;
    insert.validate().unwrap();
    let report = l::analyze(&operands(&vector, &insert), &request(false)).unwrap();
    assert_eq!(report.intended.sequence.as_str(), "AGTCCGGACTGA");
    assert!(report.enumeration_complete);
    assert!(report.unresolved_paths > 0);
    assert!(!report.exhaustive_within_scope);
    report.validate_binding(&report.context_revision).unwrap();
}

#[test]
fn above_enumeration_limit_retains_supported_ordered_material_without_attempting_alternatives() {
    assert_eq!(l::MAX_FRAGMENTS, 1000);
    assert_eq!(l::MAX_ENUMERATION_FRAGMENTS, 8);
    for n in [9, 32, 129, 1000] {
        let material = blunt("AG", true);
        let keys = (0..n).map(|i| format!("fragment-{i}")).collect::<Vec<_>>();
        let operands = keys
            .iter()
            .map(|k| Operand {
                key: k,
                material: &material,
            })
            .collect::<Vec<_>>();
        let request = Request {
            intended: keys
                .iter()
                .map(|k| Placement {
                    key: k.clone(),
                    reverse: false,
                })
                .collect(),
            circular: true,
            allow_nicked: false,
            options: Options::default(),
        };
        let report = l::analyze(&operands, &request).unwrap();
        assert_eq!(report.intended.sequence.as_str(), "AG".repeat(n));
        assert_eq!(report.paths_checked, 1);
        assert!(report.alternatives.is_empty());
        assert!(report.limits.contains(&Limit::EnumerationFragments));
        assert!(!report.enumeration_complete && !report.exhaustive_within_scope);
        report.validate_binding(&report.context_revision).unwrap();
    }
}
