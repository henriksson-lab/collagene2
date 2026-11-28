use super::*;
use crate::{Column, Orientation, Pairing, ScoreScope};
use std::collections::BTreeSet;

fn assumptions() -> Assumptions {
    Assumptions {
        sodium_molar: 0.05,
        total_strand_molar: 500e-9,
        concentration: ConcentrationConvention::EqualConcentrationTotalStrands,
        salt: SaltConvention::SodiumEntropy1998,
        ph: PhConvention::Ph7,
    }
}
fn stack(key: &str) -> Stack {
    let b = key.as_bytes();
    Stack {
        top: [b[0] as char, b[1] as char],
        bottom: [b[3] as char, b[4] as char],
    }
}

#[test]
fn primary_literal_gt_windows_replace_rather_than_add_penalties() {
    // Primary author's supplement Table 2, AG/TT and TT/AG: top AAGAA
    // against bottom TTTTT has one internal G.T pair. Independent literal facts.
    let Outcome::Replacements(plan) = evaluate(b"AAGAA", b"TTTTT") else {
        panic!("single internal GT must have a complete replacement plan");
    };
    assert_eq!(plan.mismatch_column, 2);
    assert!(plan.heteroduplex);
    let [left, right] = plan.replacements;
    assert_eq!(left.actual, stack("AG/TT"));
    assert_eq!(right.actual, stack("GA/TT"));
    assert_eq!(left.canonical, stack("AG/TC"));
    assert_eq!(right.canonical, stack("GA/CT"));
    assert_eq!(
        (left.delta_h_tenths_kcal, left.delta_s_tenths_cal_per_kelvin),
        (10, 9)
    );
    assert_eq!(
        (
            right.delta_h_tenths_kcal,
            right.delta_s_tenths_cal_per_kelvin
        ),
        (-13, -53)
    );
    // The two absolute empirical terms sum to -0.3 kcal/mol, -4.4 cal/mol/K.
    // The canonical terms being removed sum to -16.0, -43.2 (1998 table).
    assert_eq!(
        left.delta_h_tenths_kcal + right.delta_h_tenths_kcal - (-160),
        157
    );
    assert_eq!(
        left.delta_s_tenths_cal_per_kelvin + right.delta_s_tenths_cal_per_kelvin - (-432),
        388
    );
    let Outcome::Replacements(rotated) = evaluate(b"TTTTT", b"AAGAA") else {
        panic!("physically rotated duplex must remain supported");
    };
    assert_eq!(rotated.replacements[0].delta_h_tenths_kcal, -13);
    assert_eq!(rotated.replacements[1].delta_h_tenths_kcal, 10);
}

#[test]
fn every_single_mismatch_step_has_one_physical_parameter_class() {
    // Exhaust all literal A/C/G/T paired-step possibilities independently of
    // the filtered table's source ordering or the binding scanner recurrence.
    let wc = [('A', 'T'), ('C', 'G'), ('G', 'C'), ('T', 'A')];
    let mut physical = BTreeSet::new();
    let mut oriented = 0;
    for (a, b) in wc {
        for top in ['A', 'C', 'G', 'T'] {
            for bottom in ['A', 'C', 'G', 'T'] {
                if wc.contains(&(top, bottom)) {
                    continue;
                }
                for step in [
                    Stack {
                        top: [a, top],
                        bottom: [b, bottom],
                    },
                    Stack {
                        top: [top, a],
                        bottom: [bottom, b],
                    },
                ] {
                    assert!(parameter(step).is_some(), "{}", step.key());
                    assert_eq!(parameter(step), parameter(step.rotated()));
                    physical.insert(step.key().min(step.rotated().key()));
                    oriented += 1;
                }
            }
        }
    }
    assert_eq!(oriented, 96);
    assert_eq!(physical.len(), 48);
    assert_eq!(parameters::PARAMETERS.len(), 48);
    assert_eq!(
        parameters::PARAMETERS
            .iter()
            .map(|r| r.0)
            .collect::<BTreeSet<_>>()
            .len(),
        48
    );
    assert_eq!(parameter(stack("AA/TC")), Some((23, 46)));
    assert_eq!(parameter(stack("GC/CC")), Some((36, 89)));
    assert_eq!(parameter(stack("CG/GG")), Some((-49, -153)));
    assert_eq!(parameter(stack("GG/TT")), None); // tandem GT explicitly excluded
}

#[test]
fn scientific_refusals_never_masquerade_as_perfect_duplexes() {
    for (top, bottom, expected) in [
        ("AAGAA", "TTNTT", Unsupported::AmbiguousBases),
        ("A", "T", Unsupported::TooShort),
        ("AAGGA", "TTTTT", Unsupported::MultipleMismatches),
        ("GAAAA", "TTTTT", Unsupported::MismatchNearEnd),
        ("AGAAA", "TTTTT", Unsupported::MismatchNearEnd),
        ("AAAGA", "TTTTT", Unsupported::MismatchNearEnd),
        ("AAAAG", "TTTTT", Unsupported::MismatchNearEnd),
        ("GACAGTC", "CTGACAG", Unsupported::SelfAssociation),
    ] {
        assert_eq!(
            evaluate(top.as_bytes(), bottom.as_bytes()),
            Outcome::Unsupported(expected)
        );
    }
    assert_eq!(evaluate(b"ACGT", b"TGCA"), Outcome::PerfectDuplex);
    // Top ATGCAT is self-complementary, but this actual partner is distinct.
    // The existing perfect-baseline symmetry term must therefore be undone.
    assert!(matches!(
        evaluate(b"ATGCAT", b"TATGTA"),
        Outcome::Replacements(_)
    ));
    for invalid in [f64::NAN, f64::INFINITY, -1.0, 0.0, 1.01] {
        assert!(
            Assumptions {
                sodium_molar: invalid,
                ..assumptions()
            }
            .validate()
            .is_err()
        );
    }
}

/// Literal full-footprint fixture, constructed without calling the scanner.
fn fixture(reverse: bool, circular: bool) -> (Sequence, Oligo, Options, Hit) {
    let orientation = if reverse {
        Orientation::Reverse
    } else {
        Orientation::Forward
    };
    // The source is homogeneous, but coordinates deliberately wrap circularly.
    let source = Sequence::dna(if reverse { "TTTTT" } else { "AAAAA" }).unwrap();
    let oligo = Oligo {
        id: "gt".into(),
        name: "declared tail and extra prefix".into(),
        ordered: Sequence::dna("CCAAGAA").unwrap(),
        declared_tail_length: 1,
    };
    let options = Options {
        minimum_annealing_length: 5,
        maximum_extra_unpaired: 1,
        maximum_mismatches: 1,
        minimum_three_prime_pairing: 2,
        match_score: 2,
        mismatch_score: -3,
        minimum_score: 0,
    };
    let start = if circular { 3 } else { 0 };
    let positions: Vec<_> = (0..5)
        .map(|i| {
            if reverse {
                (start + 4 - i) % 5
            } else {
                (start + i) % 5
            }
        })
        .collect();
    let columns = positions
        .iter()
        .enumerate()
        .map(|(i, &position)| Column {
            oligo_index: i + 2,
            source_position: position,
            oligo_base: "AAGAA".as_bytes()[i] as char,
            source_base: if reverse { 'T' } else { 'A' },
            template_base: 'T',
            pairing: if i == 2 {
                Pairing::Mismatch
            } else {
                Pairing::Certain
            },
        })
        .collect();
    let hit = Hit {
        oligo_index: 0,
        oligo_id: "gt".into(),
        orientation,
        three_prime_source_position: positions[4],
        source_start: start,
        length: 5,
        crosses_origin: circular,
        declared_tail_length: 1,
        extra_unpaired_length: 1,
        annealing_start: 2,
        minimum_mismatches: 1,
        maximum_mismatches: 1,
        minimum_score: 5,
        maximum_score: 5,
        three_prime_certain: 2,
        three_prime_possible: 2,
        certainly_admissible: true,
        score_scope: ScoreScope::ExactConcreteOptimum,
        columns,
    };
    (source, oligo, options, hit)
}

#[test]
fn original_coordinate_tail_and_condition_bindings_are_retained_and_checked() {
    let mut identities = BTreeSet::new();
    for reverse in [false, true] {
        for circular in [false, true] {
            let (source, oligo, options, hit) = fixture(reverse, circular);
            let topology = if circular {
                Topology::Circular
            } else {
                Topology::Linear
            };
            let context = DuplexContext::new(&source, topology).unwrap();
            let assessment = context
                .assess(&oligo, &options, &hit, assumptions())
                .unwrap();
            assessment.validate().unwrap();
            assert_eq!(assessment.top_5_to_3, "AAGAA");
            assert_eq!(assessment.bottom_3_to_5, "TTTTT");
            assert_eq!(assessment.oligo_start, 2);
            assert_eq!(
                assessment.source_positions,
                hit.columns
                    .iter()
                    .map(|c| c.source_position)
                    .collect::<Vec<_>>()
            );
            assert!(matches!(assessment.outcome, Outcome::Replacements(_)));
            assessment
                .validate_for(&context, &oligo, &options, &hit)
                .unwrap();
            identities.insert(assessment.input_revision.clone());
            let changed = context
                .assess(
                    &oligo,
                    &options,
                    &hit,
                    Assumptions {
                        sodium_molar: 0.1,
                        ..assumptions()
                    },
                )
                .unwrap();
            assert_ne!(assessment.input_revision, changed.input_revision);
            let mut forged = assessment.clone();
            forged.source_positions[0] = 99;
            assert!(
                forged
                    .validate_for(&context, &oligo, &options, &hit)
                    .is_err()
            );
            let mut forged_hit = hit.clone();
            forged_hit.columns[0].template_base = 'G';
            assert!(
                context
                    .assess(&oligo, &options, &forged_hit, assumptions())
                    .is_err()
            );
        }
    }
    assert_eq!(identities.len(), 4);
}

#[test]
fn profile_identity_rejects_tampering_and_binds_license_and_policy() {
    let mut profile = Profile::current();
    profile.validate().unwrap();
    assert_eq!(profile.parameter_count, 48);
    profile.numerical_source_sha256.replace_range(0..1, "0");
    assert!(profile.validate().is_err());
    assert_eq!(SOURCE_SHA256.len(), 64);
    assert_eq!(LICENSE_SHA256.len(), 64);
}

#[test]
fn documentary_validation_rejects_fabricated_parameters_before_source_replay() {
    let (source, oligo, options, hit) = fixture(false, false);
    let context = DuplexContext::new(&source, Topology::Linear).unwrap();
    let original = context
        .assess(&oligo, &options, &hit, assumptions())
        .unwrap();
    let mut changed = original.clone();
    let Outcome::Replacements(plan) = &mut changed.outcome else {
        panic!("literal mismatch fixture");
    };
    plan.replacements[0].delta_h_tenths_kcal += 1;
    assert!(changed.validate().is_err());
    let mut changed = original.clone();
    changed.outcome = Outcome::PerfectDuplex;
    assert!(changed.validate().is_err());
    let mut changed = original.clone();
    changed.source_positions = vec![0, 1, 4, 3, 2];
    assert!(changed.validate().is_err());
    let mut changed = original.clone();
    changed.source_positions = vec![0, 1, 1, 2, 3];
    assert!(changed.validate().is_err());
    let mut changed = original.clone();
    changed.bottom_3_to_5 = "éTTT".into();
    assert!(changed.validate().is_err());
    let mut changed = original;
    changed.profile_revision = "0".repeat(64);
    assert!(changed.validate().is_err());
}
