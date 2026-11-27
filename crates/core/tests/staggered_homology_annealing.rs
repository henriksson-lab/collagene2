use collagene_core::{
    Sequence, Topology,
    chemistry::{Duplex, FivePrime, Strand, ThreePrime},
    homology_annealing::{self as annealing, Options, PlannedRescue},
    overlap_profiles::OverlapMethod,
    three_prime_resection::Cuts,
};
fn pcr(text: &str) -> Duplex {
    Duplex::exact_pcr(&Sequence::dna(text).unwrap()).unwrap()
}
fn linear() -> Options {
    Options {
        overlaps: vec![4],
        topology: Topology::Linear,
        cuts: vec![Cuts { left: 0, right: 2 }, Cuts { left: 4, right: 0 }],
        min_paired_bases: 4,
        max_comparisons: 100_000,
        max_candidates: 100,
        planned_rescue: None,
    }
}
fn five_prime_source() -> Duplex {
    // Actual top has 12 calls; actual bottom has 12, offset by 2.
    // The last two contour calls CC exist only on bottom, not on top.
    let mut source = pcr("AGTCCGATGCTA");
    source.bottom = Sequence::dna("GGTAGCATCGGA").unwrap();
    source.bottom_start = 2;
    source.left.five_prime = FivePrime::Unknown;
    source.right.five_prime = FivePrime::Phosphate;
    source
}
#[test]
fn five_prime_extensions_are_real_strand_calls_not_top_sequence_padding() {
    let sources = [five_prime_source(), pcr("TACCGGAATTCG")];
    let original = sources.clone();
    let r = annealing::geometry(&sources, &linear()).unwrap();
    r.validate_documentary().unwrap();
    assert_eq!(r.reference.as_str(), "AGTCCGATGCTACCGGAATTCG");
    assert_eq!(r.segments[0].sequence.as_str(), "AGTCCGATGC");
    assert_eq!(r.segments[1].sequence.as_str(), "GGTAGCATCGGA");
    assert_eq!(
        (r.segments[1].reference_start, r.segments[1].reference_end),
        (2, 14)
    );
    assert_eq!(r.junctions[0].homology, "TACC");
    assert_eq!(
        (r.junctions[0].top.start, r.junctions[0].bottom.start),
        (10, 14)
    );
    let d = r.paired_diagnostic.as_ref().unwrap();
    assert_eq!(d.top.as_str(), "AGTCCGATGCTACCGGAATTCG");
    assert_eq!(d.bottom_start, 2);
    assert_eq!(d.nicks.len(), 2);
    assert_eq!(r.segments[0].five_prime, FivePrime::Unknown);
    assert_eq!(r.segments[1].five_prime, FivePrime::Phosphate);
    let placement = &r.source_placements().unwrap()[0];
    assert_eq!(placement.contour_length, 14);
    assert_eq!(placement.source_top_end, 12);
    assert_eq!(sources, original);
    assert!(r.incomplete);
    annealing::replay(&sources, &linear(), &r).unwrap();
    let mut wrong = r.clone();
    wrong.segments[1].reference_end = 12;
    assert!(wrong.validate_documentary().is_err());
}
#[test]
fn three_prime_extensions_require_explicit_removal_and_do_not_enter_the_contour() {
    // Only the middle 50 calls pair. Both extra four-call 3′ extensions
    // are explicit substrate, not missing reference calls to be filled.
    let mut vector = pcr("ACGTTGCACTGATCGAGTCCGATGCTACGTAGGACGATCCAGTACGCTTATGCA");
    vector.bottom =
        Sequence::dna("TAAGCGTACTGGATCGTCCTACGTAGCATCGGACTCGATCAGTGCAACGTACGT").unwrap();
    vector.bottom_start = -4;
    let insert = pcr("GATCCAGTACGCTTATTCAGACCTGTAGATCGTCAACGTTGCACTGATCG");
    let sources = [vector, insert];
    let options = Options {
        overlaps: vec![15, 15],
        topology: Topology::Circular,
        cuts: vec![
            Cuts {
                left: 19,
                right: 19,
            },
            Cuts {
                left: 15,
                right: 15,
            },
        ],
        min_paired_bases: 15,
        max_comparisons: 100_000,
        max_candidates: 100,
        planned_rescue: Some(PlannedRescue::EcoliTransformation),
    };
    let r = annealing::geometry(&sources, &options).unwrap();
    r.validate_documentary().unwrap();
    assert_eq!(
        r.reference.as_str(),
        "ACGTTGCACTGATCGAGTCCGATGCTACGTAGGACGATCCAGTACGCTTATTCAGACCTGTAGATCGTCA"
    );
    assert_eq!(
        r.segments.iter().map(|s| s.sequence.len()).sum::<usize>(),
        140
    );
    assert_eq!(r.paired_diagnostic.as_ref().unwrap().nicks.len(), 4);
    assert!(r.segments.iter().all(|s| s.sequence.len() == 35));
    assert_eq!(r.source_placements().unwrap()[0].source_top_end, 50);
    assert_eq!(r.preparations[0].removed.len(), 2);
    assert!(
        r.junctions
            .iter()
            .all(|j| j.top.start == j.top.end && j.bottom.start == j.bottom.end)
    );
    assert!(r.incomplete);
    let mut undercut = options.clone();
    undercut.cuts[0] = Cuts {
        left: 15,
        right: 15,
    };
    assert!(annealing::geometry(&sources, &undercut).is_err());
    let mut retained_extension = options.clone();
    retained_extension.cuts[0] = Cuts { left: 3, right: 3 };
    assert!(annealing::geometry(&sources, &retained_extension).is_err());
    for group in [ThreePrime::Blocked, ThreePrime::Unknown] {
        let mut blocked = sources.clone();
        blocked[0].right.three_prime = group;
        assert!(annealing::geometry(&blocked, &options).is_err());
        assert_eq!(blocked[0].right.three_prime, group);
    }
}
#[test]
fn existing_five_prime_repeats_expose_alternatives_without_claiming_uniqueness() {
    let mut a = pcr("AAAACCCC");
    a.bottom = Sequence::dna("GGGGGGTT").unwrap(); // plus AACCCCCC, starts at 2
    a.bottom_start = 2;
    let sources = [a, pcr("CCCCGGGG")];
    let mut o = linear();
    o.min_paired_bases = 2;
    let r = annealing::geometry(&sources, &o).unwrap();
    assert_eq!(r.junctions[0].alternatives, vec![2, 3]);
    assert!(r.terminal_search_complete);
    o.max_candidates = 1;
    let r = annealing::geometry(&sources, &o).unwrap();
    assert_eq!(r.junctions[0].alternatives, vec![2]);
    assert!(!r.terminal_search_complete);
    assert!(r.incomplete);
}
#[test]
fn planned_rescue_is_required_for_infusion_but_never_seals_or_changes_material() {
    let sources = [five_prime_source(), pcr("TACCGGAATTCG")];
    let plain = linear();
    assert!(plain.validate_profile(OverlapMethod::InFusion).is_err());
    let before = annealing::geometry(&sources, &plain).unwrap();
    let mut declared = plain.clone();
    declared.planned_rescue = Some(PlannedRescue::EcoliTransformation);
    declared.validate_profile(OverlapMethod::InFusion).unwrap();
    let after = annealing::geometry(&sources, &declared).unwrap();
    assert_eq!(before.reference, after.reference);
    assert_eq!(before.preparations, after.preparations);
    assert_eq!(before.segments, after.segments);
    assert_eq!(before.junctions, after.junctions);
    assert_eq!(before.paired_diagnostic, after.paired_diagnostic);
    assert_eq!(before.incomplete, after.incomplete);
    assert!(
        after
            .paired_diagnostic
            .unwrap()
            .nicks
            .iter()
            .any(|n| n.strand == Strand::Top)
    );
    let mut json = serde_json::to_value(plain).unwrap();
    json.as_object_mut().unwrap().remove("planned_rescue");
    assert!(serde_json::from_value::<Options>(json).is_err());
}

#[test]
fn existing_exposed_tails_do_not_require_or_invent_unblocking() {
    let mut a = pcr("AGTCCGAT");
    a.bottom = Sequence::dna("TAGCATCG").unwrap();
    a.bottom_start = 4;
    a.right.three_prime = ThreePrime::Blocked;
    a.left.three_prime = ThreePrime::Unknown;
    let sources = [a, pcr("GCTATTCGAACC")];
    let mut o = linear();
    o.cuts[0] = Cuts { left: 0, right: 0 };
    let r = annealing::geometry(&sources, &o).unwrap();
    r.validate_documentary().unwrap();
    assert_eq!(r.reference.as_str(), "AGTCCGATGCTATTCGAACC");
    assert_eq!(r.junctions[0].top.three_prime, ThreePrime::Blocked);
    assert_eq!(
        r.preparations[0].product.left.three_prime,
        ThreePrime::Unknown
    );
    assert!(r.preparations[0].removed.is_empty());
    assert!(r.incomplete);
    o.cuts[0].right = 1;
    assert!(annealing::geometry(&sources, &o).is_err());
}
