use collagene_core::{
    Sequence, Topology,
    chemistry::{Duplex, FivePrime, Polarity, Strand, ThreePrime},
    homology_annealing::{self, Options},
    modified_duplex::{self, CallState, Chemical, ModificationSpan, TrackedDuplex},
    three_prime_resection::{self, Cuts},
};
const A: &str = "AGTCCGATGCTA";
const B: &str = "GCTATTCGAACC";
const CIRCLE_B: &str = "GCTATTCGAACCAGTC";
const JOINED: &str = "AGTCCGATGCTATTCGAACC";
fn pcr(s: &str) -> Duplex {
    Duplex::exact_pcr(&Sequence::dna(s).unwrap()).unwrap()
}
fn options(circle: bool) -> Options {
    Options {
        overlaps: if circle { vec![4, 4] } else { vec![4] },
        topology: if circle {
            Topology::Circular
        } else {
            Topology::Linear
        },
        cuts: if circle {
            vec![Cuts { left: 4, right: 4 }, Cuts { left: 4, right: 4 }]
        } else {
            vec![Cuts { left: 0, right: 4 }, Cuts { left: 4, right: 0 }]
        },
        min_paired_bases: 4,
        max_comparisons: 100_000,
        max_candidates: 100,
    }
}
#[test]
fn asymmetric_three_prime_removal_preserves_five_prime_groups_and_literal_strands() {
    let mut source = pcr(A);
    source.left.five_prime = FivePrime::Unknown;
    source.right.five_prime = FivePrime::Phosphate;
    let report = three_prime_resection::geometry(&source, &Cuts { left: 2, right: 3 }).unwrap();
    assert_eq!(report.product.top.as_str(), "AGTCCGATG");
    assert_eq!(report.product.bottom.as_str(), "TAGCATCGGA");
    assert_eq!(report.product.bottom_start, 2);
    assert_eq!(report.product.left.five_prime, FivePrime::Unknown);
    assert_eq!(report.product.right.five_prime, FivePrime::Phosphate);
    assert_eq!(report.product.left.three_prime, ThreePrime::Hydroxyl);
    let (left, right) = report.product.ends().unwrap();
    assert_eq!(left.polarity, Polarity::FivePrime);
    assert_eq!(left.overhang, "AG");
    assert_eq!(right.polarity, Polarity::FivePrime);
    assert_eq!(right.overhang, "TAG");
    assert_eq!(report.removed[0].start, 9);
    assert_eq!(report.removed[0].end, 12);
    assert_eq!(report.removed[1].start, 10);
    assert_eq!(report.removed[1].end, 12);
    assert!(
        report
            .mappings
            .iter()
            .all(|m| !m.synthesized && !m.reverse && m.source_start == 0)
    );
    three_prime_resection::replay(&source, &report.cuts, &report).unwrap();
    assert_eq!(source.top.as_str(), A);
}
#[test]
fn zero_cuts_keep_blocked_unknown_groups_and_positive_cuts_cannot_unblock_or_exhaust() {
    let mut source = pcr(A);
    source.left.three_prime = ThreePrime::Blocked;
    source.right.three_prime = ThreePrime::Unknown;
    assert_eq!(
        three_prime_resection::geometry(&source, &Cuts { left: 0, right: 0 })
            .unwrap()
            .product,
        source
    );
    assert!(three_prime_resection::geometry(&source, &Cuts { left: 1, right: 0 }).is_err());
    assert!(three_prime_resection::geometry(&source, &Cuts { left: 0, right: 1 }).is_err());
    let source = pcr(A);
    for cuts in [
        Cuts { left: 12, right: 0 },
        Cuts { left: 0, right: 12 },
        Cuts { left: 6, right: 6 },
    ] {
        assert!(three_prime_resection::geometry(&source, &cuts).is_err());
    }
    let mut circle = source.clone();
    circle.topology = Topology::Circular;
    assert!(three_prime_resection::geometry(&circle, &Cuts { left: 0, right: 0 }).is_err());
    let mut nicked = source;
    nicked.nicks.push(collagene_core::chemistry::Nick {
        strand: Strand::Top,
        position: 5,
    });
    assert!(three_prime_resection::geometry(&nicked, &Cuts { left: 1, right: 1 }).is_err());
}
#[test]
fn existing_staggered_free_ends_resect_without_moving_their_original_five_prime_origins() {
    let source = three_prime_resection::geometry(&pcr(A), &Cuts { left: 2, right: 3 })
        .unwrap()
        .product;
    let next = three_prime_resection::geometry(&source, &Cuts { left: 1, right: 1 })
        .unwrap()
        .product;
    assert_eq!(next.top.as_str(), "AGTCCGAT");
    assert_eq!(next.bottom.as_str(), "TAGCATCGG");
    assert_eq!(next.bottom_start, 3);
    assert_eq!(next.left.five_prime, FivePrime::Hydroxyl);
    assert_eq!(next.right.five_prime, FivePrime::Hydroxyl);
}
#[test]
fn exact_homology_circle_is_paired_but_every_covalent_join_remains_a_nick() {
    let sources = [pcr(A), pcr(CIRCLE_B)];
    let o = options(true);
    let report = homology_annealing::geometry(&sources, &o).unwrap();
    assert_eq!(report.reference.as_str(), JOINED);
    let diagnostic = report.paired_diagnostic.as_ref().unwrap();
    assert_eq!(diagnostic.top.as_str(), JOINED);
    assert_eq!(diagnostic.bottom.as_str(), "GGTTCGAATAGCATCGGACT");
    assert_eq!(diagnostic.topology, Topology::Circular);
    assert_eq!(diagnostic.nicks.len(), 4);
    assert!(report.incomplete);
    assert!(report.terminal_search_complete);
    assert!(report.junctions.iter().all(|j| j.top.start == j.top.end
        && j.bottom.start == j.bottom.end
        && j.top.five_prime == FivePrime::Hydroxyl
        && j.bottom.five_prime == FivePrime::Hydroxyl));
    assert_eq!(
        (
            report.junctions[0].top.start,
            report.junctions[0].bottom.start
        ),
        (8, 12)
    );
    assert_eq!(
        (
            report.junctions[1].top.start,
            report.junctions[1].bottom.start
        ),
        (20, 24)
    );
    assert_eq!(report.segments.len(), 4);
    assert_eq!(
        report.segments[3].sequence.as_str(),
        "GACTGGTT CGAA".replace(' ', "")
    );
    assert_eq!(
        (
            report.segments[3].reference_start,
            report.segments[3].reference_end
        ),
        (12, 24)
    );
    homology_annealing::replay(&sources, &o, &report).unwrap();
}
#[test]
fn excess_resection_leaves_real_gaps_without_synthesis_and_outer_ends_retain_five_prime_tails() {
    let mut o = options(false);
    o.cuts = vec![Cuts { left: 2, right: 6 }, Cuts { left: 7, right: 3 }];
    let r = homology_annealing::geometry(&[pcr(A), pcr(B)], &o).unwrap();
    assert_eq!(r.reference.as_str(), JOINED);
    assert!(r.paired_diagnostic.is_none());
    assert_eq!((r.junctions[0].top.start, r.junctions[0].top.end), (6, 8));
    assert_eq!(
        (r.junctions[0].bottom.start, r.junctions[0].bottom.end),
        (12, 15)
    );
    assert_eq!(r.segments[0].sequence.as_str(), "AGTCCG");
    assert_eq!(r.segments[2].sequence.as_str(), "GCTATTCGA");
    assert_eq!(
        r.preparations[0].product.left.five_prime,
        FivePrime::Hydroxyl
    );
    let mut o = options(false);
    o.cuts[0].left = 2;
    o.cuts[1].right = 3;
    let r = homology_annealing::geometry(&[pcr(A), pcr(B)], &o).unwrap();
    let p = r.paired_diagnostic.as_ref().unwrap();
    assert_eq!(p.top.as_str(), "AGTCCGATGCTATTCGA");
    assert_eq!(p.bottom.as_str(), "GGTTCGAATAGCATCGGA");
    assert_eq!(p.bottom_start, 2);
    assert_eq!(p.nicks.len(), 2);
    assert_eq!(p.left.five_prime, FivePrime::Hydroxyl);
    assert_eq!(p.right.five_prime, FivePrime::Hydroxyl);
}
#[test]
fn selected_source_homology_and_exposure_are_required_and_work_caps_do_not_claim_uniqueness() {
    let sources = [pcr("AAAACCCC"), pcr("CCCCGGGG")];
    let mut o = options(false);
    o.min_paired_bases = 2;
    let r = homology_annealing::geometry(&sources, &o).unwrap();
    assert_eq!(r.junctions[0].alternatives, vec![2, 3]);
    o.max_candidates = 1;
    let r = homology_annealing::geometry(&sources, &o).unwrap();
    assert!(!r.terminal_search_complete);
    assert_eq!(r.junctions[0].alternatives, vec![2]);
    o.max_comparisons = 4;
    let r = homology_annealing::geometry(&sources, &o).unwrap();
    assert!(!r.terminal_search_complete);
    assert_eq!(r.comparisons, 4);
    o.max_comparisons = 3;
    assert!(homology_annealing::geometry(&sources, &o).is_err());
    o = options(false);
    o.cuts[1].left = 3;
    assert!(homology_annealing::geometry(&sources, &o).is_err());
    o = options(false);
    assert!(homology_annealing::geometry(&[pcr(A), pcr("ACGTTTCGAACC")], &o).is_err());
    let original = homology_annealing::geometry(&[pcr(A), pcr(B)], &o).unwrap();
    let mut forged = original.clone();
    forged.segments[0].five_prime = FivePrime::Phosphate;
    assert!(homology_annealing::replay(&[pcr(A), pcr(B)], &o, &forged).is_err());
}
#[test]
fn one_fragment_linear_processing_and_self_circularization_keep_distinct_contracts() {
    let o = Options {
        overlaps: vec![],
        topology: Topology::Linear,
        cuts: vec![Cuts { left: 2, right: 3 }],
        min_paired_bases: 4,
        max_comparisons: 100,
        max_candidates: 10,
    };
    let r = homology_annealing::geometry(&[pcr(A)], &o).unwrap();
    assert!(r.junctions.is_empty());
    assert!(!r.incomplete);
    assert_eq!(r.paired_diagnostic.unwrap().top.as_str(), "AGTCCGATG");
    let o = Options {
        overlaps: vec![4],
        topology: Topology::Circular,
        cuts: vec![Cuts { left: 4, right: 4 }],
        ..o
    };
    let r = homology_annealing::geometry(&[pcr("AGTCCGATGCTAAGTC")], &o).unwrap();
    assert_eq!(r.reference.as_str(), A);
    assert_eq!(r.paired_diagnostic.unwrap().nicks.len(), 2);
    assert!(r.incomplete);
}
#[test]
fn tracked_resection_keeps_only_surviving_original_modifications_and_reversed_coordinates() {
    let mut source = TrackedDuplex::from_unknown_duplex(&pcr(A)).unwrap();
    source.top_modifications.spans = vec![ModificationSpan {
        chemical: Chemical::C5C,
        start: 3,
        end: 5,
        state: CallState::Modified,
    }];
    source.bottom_modifications.spans = vec![ModificationSpan {
        chemical: Chemical::C5C,
        start: 10,
        end: 11,
        state: CallState::Modified,
    }];
    let (processed, r) = source
        .resect_three_prime(&Cuts { left: 2, right: 4 })
        .unwrap();
    assert_eq!(processed.top.as_str(), "AGTCCGAT");
    assert_eq!(
        processed.state_at(Strand::Top, Chemical::C5C, 3).unwrap(),
        Some(CallState::Modified)
    );
    assert!(
        processed
            .bottom_modifications
            .spans
            .iter()
            .all(|s| s.end <= 10)
    );
    assert_eq!(processed.generation, source.generation);
    assert!(r.mappings.iter().all(|m| !m.synthesized));
    let reversed = source.reverse().unwrap();
    let (processed, _) = reversed
        .resect_three_prime(&Cuts { left: 4, right: 2 })
        .unwrap();
    assert_eq!(processed.bottom.as_str(), "AGTCCGAT");
    assert_eq!(
        processed
            .state_at(Strand::Bottom, Chemical::C5C, 3)
            .unwrap(),
        Some(CallState::Modified)
    );
    let a = TrackedDuplex::from_unknown_duplex(&pcr(A)).unwrap();
    let b = TrackedDuplex::from_unknown_duplex(&pcr(B)).unwrap();
    let (operands, report) = modified_duplex::homology_anneal(&[&a, &b], &options(false)).unwrap();
    assert_eq!(operands.len(), 2);
    assert_eq!(operands[0].top.as_str(), "AGTCCGAT");
    assert_eq!(operands[1].bottom.as_str(), "GGTTCGAA");
    assert!(report.incomplete);
}
