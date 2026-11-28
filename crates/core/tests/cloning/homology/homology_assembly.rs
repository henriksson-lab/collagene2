use collagene_core::{
    chemistry::{Duplex, FivePrime, Nick, Strand, ThreePrime},
    homology_assembly::{self, EndCuts, Options, OuterEndPolicy, Resection},
    modified_duplex::{self, CallState, Chemical, ModificationSpan, TrackedDuplex},
};
use collagene_model::{Sequence, Topology};

const A: &str = "AGTCCGATGCTA";
const B: &str = "GCTATTCGAACC";
const CIRCLE_B: &str = "GCTATTCGAACCAGTC";
const EXPECTED: &str = "AGTCCGATGCTATTCGAACC";
fn pcr(s: &str) -> Duplex {
    Duplex::exact_pcr(&Sequence::dna(s).unwrap()).unwrap()
}
fn options(circular: bool) -> Options {
    Options {
        overlaps: if circular { vec![4, 4] } else { vec![4] },
        topology: if circular {
            Topology::Circular
        } else {
            Topology::Linear
        },
        min_paired_bases: 4,
        ..Default::default()
    }
}

#[test]
fn pcr_hydroxyl_inputs_make_independently_tiled_sealed_circle() {
    let sources = [pcr(A), pcr(CIRCLE_B)];
    let r = homology_assembly::geometry(&sources, &options(true)).unwrap();
    let p = r.product.as_ref().unwrap();
    assert_eq!(p.top.as_str(), EXPECTED);
    assert_eq!(p.bottom.as_str(), "GGTTCGAATAGCATCGGACT");
    assert_eq!(p.topology, Topology::Circular);
    assert_eq!(p.bottom_start, 0);
    assert!(p.nicks.is_empty());
    assert_eq!(r.reference_top_start, 0);
    assert_eq!(r.reference_top_end, 20);
    assert!(r.junctions.iter().all(|j| j.top_sealed
        && j.bottom_sealed
        && j.top_five_prime == FivePrime::Phosphate
        && j.bottom_five_prime == FivePrime::Phosphate));
    assert!(
        sources
            .iter()
            .all(|s| s.left.five_prime == FivePrime::Hydroxyl
                && s.right.five_prime == FivePrime::Hydroxyl)
    );
    assert!(r.mappings.iter().all(|m| !m.mapping.synthesized));
    assert_eq!(
        r.resected.iter().map(|s| s.end - s.start).sum::<usize>(),
        16
    );
    r.validate_documentary().unwrap();
}

#[test]
fn asymmetric_extra_resection_uses_actual_opposite_templates_and_standard_new_dna() {
    let mut o = options(false);
    o.resection = Resection::Explicit {
        cuts: vec![EndCuts { left: 0, right: 6 }, EndCuts { left: 7, right: 0 }],
    };
    let mut a = TrackedDuplex::from_unknown_duplex(&pcr(A)).unwrap();
    let mut b = TrackedDuplex::from_unknown_duplex(&pcr(B)).unwrap();
    // Top A[6]=A survives. Top B[4..7] is removed; its methylated A at8 survives.
    a.top_modifications.spans.push(ModificationSpan {
        chemical: Chemical::N6A,
        start: 6,
        end: 7,
        state: CallState::Modified,
    });
    b.top_modifications.spans.push(ModificationSpan {
        chemical: Chemical::N6A,
        start: 3,
        end: 4,
        state: CallState::Modified,
    });
    b.top_modifications.spans.push(ModificationSpan {
        chemical: Chemical::N6A,
        start: 8,
        end: 9,
        state: CallState::Modified,
    });
    // Bottom B[6]=A remains on the physical bottom strand. Its modification
    // never transfers through a complementary synthesis mapping.
    b.bottom_modifications.spans.push(ModificationSpan {
        chemical: Chemical::N6A,
        start: 6,
        end: 7,
        state: CallState::Modified,
    });
    let (tracked, r) = modified_duplex::homology_assemble(&[&a, &b], &o).unwrap();
    let p = tracked.unwrap();
    assert_eq!(p.top.as_str(), EXPECTED);
    assert_eq!(p.bottom.as_str(), "GGTTCGAATAGCATCGGACT");
    assert_eq!(r.junctions[0].top_gap, 3);
    assert_eq!(r.junctions[0].bottom_gap, 2);
    assert_eq!(
        p.state_at(Strand::Top, Chemical::N6A, 6).unwrap(),
        Some(CallState::Modified)
    );
    assert_eq!(
        p.state_at(Strand::Top, Chemical::N6A, 16).unwrap(),
        Some(CallState::Modified)
    );
    assert_eq!(
        p.state_at(Strand::Top, Chemical::N6A, 11).unwrap(),
        Some(CallState::Unknown)
    );
    assert_eq!(
        p.state_at(Strand::Bottom, Chemical::N6A, 6).unwrap(),
        Some(CallState::Modified)
    );
    assert_eq!(
        p.state_at(Strand::Top, Chemical::C5C, 14).unwrap(),
        Some(CallState::Unmodified)
    );
    assert_eq!(
        p.state_at(Strand::Bottom, Chemical::N6A, 12).unwrap(),
        Some(CallState::Unmodified)
    );
    let synth: Vec<_> = r
        .mappings
        .iter()
        .filter(|m| m.mapping.synthesized)
        .collect();
    assert_eq!(
        synth
            .iter()
            .map(|m| m.mapping.product_end - m.mapping.product_start)
            .sum::<usize>(),
        5
    );
    assert!(
        synth
            .iter()
            .all(|m| m.mapping.reverse && m.mapping.source_strand != m.mapping.product_strand)
    );
    assert!(p.nicks.is_empty());
}

#[test]
fn linear_outer_policy_preserves_zero_cut_groups_and_real_processed_three_prime_tails() {
    let sources = [pcr(A), pcr(B)];
    let r = homology_assembly::geometry(&sources, &options(false)).unwrap();
    let p = r.product.unwrap();
    assert_eq!(p.top.as_str(), EXPECTED);
    assert_eq!(p.left.five_prime, FivePrime::Hydroxyl);
    assert_eq!(p.right.five_prime, FivePrime::Hydroxyl);
    assert_eq!(p.bottom_start, 0);
    assert_eq!(p.top.len(), p.bottom.len());
    let mut o = options(false);
    o.outer_end_policy = OuterEndPolicy::AllAccessibleEnds { left: 2, right: 3 };
    let r = homology_assembly::geometry(&sources, &o).unwrap();
    let p = r.product.as_ref().unwrap();
    assert_eq!(p.top.as_str(), "TCCGATGCTATTCGAACC");
    assert_eq!(p.bottom.as_str(), "TCGAATAGCATCGGACT");
    assert_eq!(p.bottom_start, -2);
    assert_eq!(r.reference_top_start, 2);
    assert_eq!(r.reference_top_end, 20);
    assert_eq!(p.left.five_prime, FivePrime::Phosphate);
    assert_eq!(p.right.five_prime, FivePrime::Phosphate);
    assert_eq!(r.junctions[0].homology_start, 6);
    assert_eq!(r.junctions[0].homology_end, 10);
    assert!(r.mappings.iter().all(|m| !m.mapping.synthesized));
    r.validate_documentary().unwrap();
    o.resection = Resection::Explicit {
        cuts: vec![EndCuts { left: 0, right: 4 }, EndCuts { left: 4, right: 0 }],
    };
    assert!(homology_assembly::geometry(&sources, &o).is_err());
}

#[test]
fn unknown_or_blocked_three_prime_groups_retain_nicks_and_cannot_fill_gaps() {
    let mut sources = [pcr(A), pcr(B)];
    sources[0].right.three_prime = ThreePrime::Blocked;
    let r = homology_assembly::geometry(&sources, &options(false)).unwrap();
    assert!(r.incomplete);
    assert!(!r.junctions[0].top_sealed);
    assert_eq!(
        r.product.as_ref().unwrap().nicks,
        vec![Nick {
            strand: Strand::Top,
            position: 12
        }]
    );
    let mut tampered = r.clone();
    tampered.product.as_mut().unwrap().nicks.clear();
    assert!(tampered.validate_documentary().is_err());
    let tracked: Vec<_> = sources
        .iter()
        .map(|s| TrackedDuplex::from_unknown_duplex(s).unwrap())
        .collect();
    let (private, diagnostic) =
        modified_duplex::homology_assemble(&tracked.iter().collect::<Vec<_>>(), &options(false))
            .unwrap();
    assert!(private.is_none());
    assert!(diagnostic.product.is_some());
    assert!(diagnostic.incomplete);
    assert_eq!(diagnostic.junctions[0].top_three_prime, ThreePrime::Blocked);
    assert!(
        diagnostic
            .limits
            .iter()
            .any(|s| s.contains("nick representation"))
    );
    let mut o = options(false);
    o.resection = Resection::Explicit {
        cuts: vec![EndCuts { left: 0, right: 4 }, EndCuts { left: 5, right: 0 }],
    };
    let r = homology_assembly::geometry(&sources, &o).unwrap();
    assert!(r.product.is_none());
    assert!(r.incomplete);
    assert!(r.mappings.is_empty());
    sources[0].right.three_prime = ThreePrime::Unknown;
    assert!(
        homology_assembly::geometry(&sources, &o)
            .unwrap()
            .product
            .is_none()
    );
}

#[test]
fn reversing_actual_inputs_preserves_independent_stored_strands_and_chemistry() {
    let sources = [pcr(B).reverse().unwrap(), pcr(A).reverse().unwrap()];
    let r = homology_assembly::geometry(&sources, &options(false)).unwrap();
    let p = r.product.unwrap();
    assert_eq!(p.top.as_str(), "GGTTCGAATAGCATCGGACT");
    assert_eq!(p.bottom.as_str(), EXPECTED);
    assert!(p.nicks.is_empty());
    assert_eq!(p.left.five_prime, FivePrime::Hydroxyl);
}

#[test]
fn bounds_competitors_mismatches_and_documentary_tampering_remain_honest() {
    let sources = [pcr(A), pcr(B)];
    let o = options(false);
    let mut mismatch = sources.clone();
    mismatch[1] = pcr("ACTATTCGAACC");
    assert!(homology_assembly::geometry(&mismatch, &o).is_err());
    let mut exhausted = o.clone();
    exhausted.resection = Resection::Explicit {
        cuts: vec![
            EndCuts { left: 0, right: 12 },
            EndCuts { left: 4, right: 0 },
        ],
    };
    assert!(homology_assembly::geometry(&sources, &exhausted).is_err());
    let mut too_small = o.clone();
    too_small.max_comparisons = 1;
    assert!(homology_assembly::preflight(&sources, &too_small).is_err());
    let mut capped = o.clone();
    capped.max_comparisons = 3 * 24 + 2 * 20 + 4;
    let r = homology_assembly::geometry(&sources, &capped).unwrap();
    assert!(r.incomplete);
    assert!(!r.terminal_search_complete);
    assert!(r.product.is_some());
    let mut false_complete = r.clone();
    false_complete.incomplete = false;
    assert!(false_complete.validate_documentary().is_err());
    let repeats = [pcr("AAAAAAAAAAAA"), pcr("AAAAAAAAAAAA")];
    let r = homology_assembly::geometry(&repeats, &o).unwrap();
    assert!(r.incomplete);
    assert!(!r.junctions[0].alternatives.is_empty());
    assert!(r.product.is_some());
    let mut tampered = r.clone();
    tampered
        .mappings
        .iter_mut()
        .find(|m| m.mapping.source_start > 0)
        .unwrap()
        .mapping
        .source_start = 0;
    assert!(tampered.validate_documentary().is_err());
    let mut nicked = sources.clone();
    nicked[0].nicks.push(Nick {
        strand: Strand::Top,
        position: 2,
    });
    assert!(homology_assembly::geometry(&nicked, &o).is_err());
    let mut wrong_shape = o;
    wrong_shape.overlaps.push(4);
    assert!(wrong_shape.validate_shape(2).is_err());
}

#[test]
fn single_fragment_self_circle_and_single_linear_outer_processing_are_real_geometry() {
    let source = pcr("AGTCCGATGCTAAGTC");
    let o = Options {
        overlaps: vec![4],
        min_paired_bases: 4,
        ..Default::default()
    };
    let r = homology_assembly::geometry(&[source], &o).unwrap();
    let p = r.product.unwrap();
    assert_eq!(p.top.as_str(), A);
    assert_eq!(p.bottom.as_str(), "TAGCATCGGACT");
    assert!(p.nicks.is_empty());
    let o = Options {
        overlaps: vec![],
        topology: Topology::Linear,
        outer_end_policy: OuterEndPolicy::AllAccessibleEnds { left: 2, right: 3 },
        min_paired_bases: 4,
        ..Default::default()
    };
    let r = homology_assembly::geometry(&[pcr(A)], &o).unwrap();
    let p = r.product.unwrap();
    assert_eq!(p.top.as_str(), "TCCGATGCTA");
    assert_eq!(p.bottom.as_str(), "CATCGGACT");
    assert_eq!(p.bottom_start, -2);
    assert!(r.junctions.is_empty());
}
