//! Conservation of actual retained strands through reciprocal extension.
use collagene_core::{
    chemistry::{
        FivePrime, Nick, Strand, ThreePrime,
        tracked_duplex::{
            CallState, Chemical, GenerationAssumptions, GenerationRecord, ModificationSpan,
            StrandModifications, TrackedDuplex,
        },
    },
    cloning::overlap::extension::Options,
};
use collagene_model::{Sequence, Topology};

const LEFT: &str = "AGTCCGATGCTACGTATTTTCAGACC";
const RIGHT: &str = "GCTACGTATTTTCAGACCTGTAGATCGTCA";
const PRODUCT: &str = "AGTCCGATGCTACGTATTTTCAGACCTGTAGATCGTCA";

fn options() -> Options {
    Options {
        overlap: 18,
        min_paired_bases: 12,
        ..Options::default()
    }
}
fn uniform(state: CallState) -> StrandModifications {
    StrandModifications {
        n6a_default: state,
        c5c_default: state,
        spans: vec![],
    }
}
fn record(n: usize) -> GenerationRecord {
    GenerationRecord {
        source_revision: format!("{n:064x}"),
        assumptions: GenerationAssumptions::default(),
    }
}
fn sources() -> (TrackedDuplex, TrackedDuplex) {
    let mut a = TrackedDuplex::exact_pcr(&Sequence::dna(LEFT).unwrap()).unwrap();
    let mut b = TrackedDuplex::exact_pcr(&Sequence::dna(RIGHT).unwrap()).unwrap();
    a.generation = vec![record(1), record(2)];
    b.generation = vec![record(2), record(3)];
    (a, b)
}
fn mark_endpoints(mods: &mut StrandModifications, sequence: &str, chemical: Chemical) {
    let base = match chemical {
        Chemical::N6A => 'A',
        Chemical::C5C => 'C',
    };
    let first = sequence.find(base).unwrap();
    let last = sequence.rfind(base).unwrap();
    for start in [Some(first), (last != first).then_some(last)]
        .into_iter()
        .flatten()
    {
        mods.spans.push(ModificationSpan {
            chemical,
            start,
            end: start + 1,
            state: CallState::Modified,
        });
    }
}

#[test]
fn retained_strand_states_are_copied_and_new_complements_are_unmodified() {
    let (mut a, mut b) = sources();
    a.top_modifications = uniform(CallState::Unknown);
    b.bottom_modifications = uniform(CallState::Unknown);
    // These discarded strands deliberately carry the opposite global state.
    a.bottom_modifications = uniform(CallState::Modified);
    b.top_modifications = uniform(CallState::Modified);
    for chemical in [Chemical::N6A, Chemical::C5C] {
        mark_endpoints(&mut a.top_modifications, a.duplex.top.as_str(), chemical);
        mark_endpoints(
            &mut b.bottom_modifications,
            b.duplex.bottom.as_str(),
            chemical,
        );
    }
    let before = (a.clone(), b.clone());
    let (product, geometry) = a.overlap_extend(&b, &options()).unwrap();
    assert_eq!(product.top.as_str(), PRODUCT);
    assert_eq!(product.duplex, geometry.fill.duplex);
    assert_eq!(geometry.annealed.bottom_start, 8);
    assert_eq!(geometry.fill.top_added_3prime.len(), 12);
    assert_eq!(geometry.fill.bottom_added_3prime.len(), 8);
    assert_eq!(product.generation, vec![record(1), record(2), record(3)]);
    for (strand, original, sequence) in [
        (Strand::Top, &a, product.top.as_str()),
        (Strand::Bottom, &b, product.bottom.as_str()),
    ] {
        let retained = match strand {
            Strand::Top => a.top.len(),
            Strand::Bottom => b.bottom.len(),
        };
        for (index, byte) in sequence.bytes().enumerate() {
            let chemical = match byte {
                b'A' => Chemical::N6A,
                b'C' => Chemical::C5C,
                _ => continue,
            };
            let expected = if index < retained {
                original.state_at(strand, chemical, index).unwrap()
            } else {
                Some(CallState::Unmodified)
            };
            assert_eq!(product.state_at(strand, chemical, index).unwrap(), expected);
        }
    }
    assert_eq!((a, b), before, "extension must leave both sources intact");
}

#[test]
fn unused_strand_marks_cannot_change_product_modification_state() {
    let (mut a, mut b) = sources();
    let (baseline, _) = a.overlap_extend(&b, &options()).unwrap();
    a.bottom_modifications = uniform(CallState::Modified);
    b.top_modifications = uniform(CallState::Unknown);
    let (product, _) = a.overlap_extend(&b, &options()).unwrap();
    assert_eq!(product, baseline);
}

#[test]
fn terminal_groups_are_inherited_without_inventing_phosphates() {
    for left_five in [
        FivePrime::Hydroxyl,
        FivePrime::Phosphate,
        FivePrime::Unknown,
    ] {
        for right_five in [
            FivePrime::Hydroxyl,
            FivePrime::Phosphate,
            FivePrime::Unknown,
        ] {
            let (mut a, mut b) = sources();
            a.duplex.left.five_prime = left_five;
            b.duplex.right.five_prime = right_five;
            // These 3′ groups belong only to discarded source strands.
            a.duplex.left.three_prime = ThreePrime::Blocked;
            b.duplex.right.three_prime = ThreePrime::Unknown;
            let (product, _) = a.overlap_extend(&b, &options()).unwrap();
            assert_eq!(product.left.five_prime, left_five);
            assert_eq!(product.right.five_prime, right_five);
            assert_eq!(product.left.three_prime, ThreePrime::Hydroxyl);
            assert_eq!(product.right.three_prime, ThreePrime::Hydroxyl);
        }
    }
    for blocked in [ThreePrime::Unknown, ThreePrime::Blocked] {
        let (mut a, b) = sources();
        a.duplex.right.three_prime = blocked;
        assert!(a.overlap_extend(&b, &options()).is_err());
        let (a, mut b) = sources();
        b.duplex.left.three_prime = blocked;
        assert!(a.overlap_extend(&b, &options()).is_err());
    }
}

#[test]
fn source_geometry_and_tracking_cannot_be_bypassed() {
    let (a, b) = sources();
    for source in 0..2 {
        for defect in 0..5 {
            let (mut left, mut right) = (a.clone(), b.clone());
            let operand = if source == 0 { &mut left } else { &mut right };
            match defect {
                0 => operand.duplex.nicks.push(Nick {
                    strand: Strand::Top,
                    position: 3,
                }),
                1 => operand.duplex.bottom_start = 1,
                2 => operand.duplex.topology = Topology::Circular,
                3 => operand.duplex.bottom_start = i64::MAX,
                _ => operand.top_modifications.spans.push(ModificationSpan {
                    chemical: Chemical::N6A,
                    start: 0,
                    end: usize::MAX,
                    state: CallState::Unknown,
                }),
            }
            assert!(
                left.overlap_extend(&right, &options()).is_err(),
                "source {source}, defect {defect}"
            );
        }
    }
    for overlap in [0, LEFT.len(), usize::MAX] {
        let invalid = Options {
            overlap,
            ..options()
        };
        assert!(a.overlap_extend(&b, &invalid).is_err());
    }
    let invalid = Options {
        max_comparisons: usize::MAX,
        ..options()
    };
    assert!(a.overlap_extend(&b, &invalid).is_err());
    let invalid = Options {
        max_candidates: usize::MAX,
        ..options()
    };
    assert!(a.overlap_extend(&b, &invalid).is_err());
}

#[test]
fn union_of_valid_source_generations_is_bounded_and_deduplicated() {
    let (mut a, mut b) = sources();
    a.generation = (0..600).map(record).collect();
    b.generation = (500..1000).map(record).collect();
    let (product, _) = a.overlap_extend(&b, &options()).unwrap();
    assert_eq!(
        product.generation,
        (0..1000).map(record).collect::<Vec<_>>()
    );
    b.generation.push(record(1000));
    // Both operands separately validate, but their 1001-record union must fail.
    a.validate().unwrap();
    b.validate().unwrap();
    assert!(a.overlap_extend(&b, &options()).is_err());
}

#[test]
fn reversal_is_performed_on_tracked_operands_before_extension() {
    let (a, b) = sources();
    let (forward, _) = a.overlap_extend(&b, &options()).unwrap();
    let (reverse, _) = b
        .reverse()
        .unwrap()
        .overlap_extend(&a.reverse().unwrap(), &options())
        .unwrap();
    let expected = forward.reverse().unwrap();
    assert_eq!(reverse.duplex, expected.duplex);
    assert_eq!(reverse.top_modifications, expected.top_modifications);
    assert_eq!(reverse.bottom_modifications, expected.bottom_modifications);
    assert_eq!(reverse.generation, vec![record(2), record(3), record(1)]);
}

#[test]
fn documentary_geometry_rejects_tampered_shape_without_restoring_material() {
    let (a, b) = sources();
    let (_, geometry) = a.overlap_extend(&b, &options()).unwrap();
    geometry.validate_documentary().unwrap();
    for tamper in 0..20 {
        let mut changed = geometry.clone();
        match tamper {
            0 => changed.engine_version.push('x'),
            1 => changed.source_geometry_revisions[0] = "z".repeat(64),
            2 => changed.options.overlap = usize::MAX,
            3 => changed.options.max_candidates = usize::MAX,
            4 => changed.options.max_comparisons = usize::MAX,
            5 => changed.homology.push('A'),
            6 => changed.annealed.bottom_start = i64::MAX,
            7 => changed.fill.source.left.five_prime = FivePrime::Phosphate,
            8 => changed.fill.top_added_3prime.push('A'),
            9 => changed.fill.duplex.right.five_prime = FivePrime::Phosphate,
            10 => changed.mappings[0].source = 2,
            11 => changed.mappings[0].mapping.product_end = usize::MAX,
            12 => changed.fill.mappings[0].source_end = usize::MAX,
            13 => changed.terminal_pairings.overlaps.clear(),
            14 => changed
                .terminal_pairings
                .overlaps
                .push(changed.options.overlap),
            15 => changed.terminal_pairings.first_alternative = Some(changed.options.overlap),
            16 => changed.terminal_pairings.comparisons = usize::MAX,
            17 => changed.terminal_pairings.unique = false,
            18 => changed.terminal_pairings.limits.push("Unsearched".into()),
            _ => changed.warnings.push("x".repeat(4097)),
        }
        assert!(changed.validate_documentary().is_err(), "tamper {tamper}");
    }
    // Incomplete search remains valid documentary evidence, never a positive
    // uniqueness claim or reconstructed source authority.
    let partial_options = Options {
        max_comparisons: options().overlap,
        ..options()
    };
    let (_, partial) = a.overlap_extend(&b, &partial_options).unwrap();
    assert!(!partial.terminal_pairings.complete);
    partial.validate_documentary().unwrap();
    let repetitive = TrackedDuplex::exact_pcr(&Sequence::dna("AAAAAA").unwrap()).unwrap();
    let competing_options = Options {
        overlap: 2,
        min_paired_bases: 1,
        max_candidates: 1,
        ..Options::default()
    };
    let (_, competing) = repetitive
        .overlap_extend(&repetitive, &competing_options)
        .unwrap();
    assert!(!competing.terminal_pairings.complete);
    assert!(competing.terminal_pairings.first_alternative.is_some());
    competing.validate_documentary().unwrap();
}
