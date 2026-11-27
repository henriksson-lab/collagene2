//! Independently specified strand-state conservation and synthesis fixtures.
use collagene_core::{
    analysis,
    chemistry::{self, EndRepairMethod, EndTreatment, Strand},
    methylation::{MethylationState, Outcome},
    modified_duplex::*,
    tailing::TailingBase,
    *,
};
fn full() -> GenerationAssumptions {
    GenerationAssumptions {
        dam: MethylationState::FullyModified,
        dcm: MethylationState::FullyModified,
        cpg: MethylationState::FullyModified,
    }
}
fn tracked(s: &str, t: Topology) -> TrackedDuplex {
    TrackedDuplex::from_input(&Molecule::dna("source", s, t).unwrap(), &full()).unwrap()
}
fn enzyme(name: &str) -> analysis::Enzyme {
    analysis::starter_enzymes()
        .into_iter()
        .find(|e| e.name == name)
        .unwrap()
}
#[test]
fn source_generation_actual_strand_indices_and_unknown_defaults_are_distinct() {
    let source = tracked("GATCCCAGGCG", Topology::Linear);
    assert_eq!(
        source.state_at(Strand::Top, Chemical::N6A, 1).unwrap(),
        Some(CallState::Modified)
    );
    assert_eq!(
        source.state_at(Strand::Bottom, Chemical::N6A, 8).unwrap(),
        Some(CallState::Modified)
    );
    assert_eq!(
        source.state_at(Strand::Top, Chemical::C5C, 5).unwrap(),
        Some(CallState::Modified)
    );
    assert_eq!(
        source.state_at(Strand::Top, Chemical::C5C, 9).unwrap(),
        Some(CallState::Modified)
    );
    let unknown = TrackedDuplex::from_unknown_duplex(&source.duplex).unwrap();
    assert_eq!(
        unknown.state_at(Strand::Top, Chemical::N6A, 6).unwrap(),
        Some(CallState::Unknown)
    );
    assert_eq!(
        source.state_at(Strand::Top, Chemical::N6A, 6).unwrap(),
        Some(CallState::Unmodified)
    );
    assert_eq!(source.reverse().unwrap().reverse().unwrap(), source);
}
#[test]
fn all_origins_and_both_strands_keep_chemical_state_through_staggered_digest_and_circular_ligation()
{
    let text = "AACGATCCAGGTCTAGAGCGGCCGCAAGGATCCTTGAATTCAAGCTTCCCGGGGATATCGGTCTCAAAGCTGCAGCGTCTCAAGAAGACTTGCTCTTCAAAAAAAA";
    for at in 0..text.len() {
        let s = format!("{}{}", &text[at..], &text[..at]);
        for reverse in [false, true] {
            let source = tracked(&s, Topology::Circular);
            let source = if reverse {
                source.reverse().unwrap()
            } else {
                source
            };
            for name in ["EcoRI", "BamHI", "PstI", "NotI", "SmaI"] {
                let fragments = digest(&source, &[enzyme(name)]).unwrap();
                let origin = fragments[0].top_start;
                let product = ligate(
                    &fragments
                        .iter()
                        .map(|f| f.duplex.clone())
                        .collect::<Vec<_>>(),
                    true,
                    false,
                )
                .unwrap();
                let n = source.duplex.top.len();
                for j in 0..n {
                    for c in [Chemical::N6A, Chemical::C5C] {
                        assert_eq!(
                            product.state_at(Strand::Top, c, j).unwrap(),
                            source.state_at(Strand::Top, c, (origin + j) % n).unwrap(),
                            "top {name} origin{at}"
                        );
                        let reference = (origin + n - 1 - j) % n;
                        assert_eq!(
                            product.state_at(Strand::Bottom, c, j).unwrap(),
                            source
                                .state_at(Strand::Bottom, c, n - 1 - reference)
                                .unwrap(),
                            "bottom {name} origin{at}"
                        );
                    }
                }
            }
        }
    }
}
#[test]
fn fill_in_and_tailing_preserve_existing_calls_but_do_not_copy_template_modifications() {
    let source = tracked("AAGAATTCGATCCAGGCGTT", Topology::Circular);
    let pieces = digest(&source, &[enzyme("EcoRI")]).unwrap();
    for piece in pieces {
        let repaired = piece
            .duplex
            .repair(EndRepairMethod::FillFivePrimeOverhangs)
            .unwrap();
        for s in [Strand::Top, Strand::Bottom] {
            let before = if s == Strand::Top {
                &piece.duplex.duplex.top
            } else {
                &piece.duplex.duplex.bottom
            };
            let after = if s == Strand::Top {
                &repaired.duplex.top
            } else {
                &repaired.duplex.bottom
            };
            for i in 0..after.len() {
                for c in [Chemical::N6A, Chemical::C5C] {
                    let expected = if i < before.len() {
                        piece.duplex.state_at(s, c, i).unwrap()
                    } else {
                        let b = after.as_str().as_bytes()[i];
                        if b == if c == Chemical::N6A { b'A' } else { b'C' } {
                            Some(CallState::Unmodified)
                        } else {
                            None
                        }
                    };
                    assert_eq!(repaired.state_at(s, c, i).unwrap(), expected);
                }
            }
        }
        let tailed = repaired.tail(TailingBase::A).unwrap();
        assert_eq!(
            tailed
                .state_at(Strand::Top, Chemical::N6A, repaired.duplex.top.len())
                .unwrap(),
            Some(CallState::Unmodified)
        );
        assert_eq!(
            tailed
                .state_at(Strand::Bottom, Chemical::N6A, repaired.duplex.bottom.len())
                .unwrap(),
            Some(CallState::Unmodified)
        );
        let pcr = TrackedDuplex::exact_pcr(&tailed.duplex.top).unwrap();
        assert_eq!(summary(&pcr).unwrap().known_modified, 0);
        assert_eq!(summary(&pcr).unwrap().unknown_eligible, 0);
    }
}
#[test]
fn chemical_compatibility_conflicting_spans_and_foreign_metadata_are_rejected() {
    let source = tracked("GATCCG", Topology::Linear);
    let mut bad = source.clone();
    bad.top_modifications.spans.push(ModificationSpan {
        chemical: Chemical::N6A,
        start: 0,
        end: 1,
        state: CallState::Modified,
    });
    assert!(bad.validate().is_err());
    let mut bad = source.clone();
    bad.top_modifications.spans.push(ModificationSpan {
        chemical: Chemical::N6A,
        start: 1,
        end: 2,
        state: CallState::Unknown,
    });
    assert!(bad.validate().is_err());
    let exported = source.duplex.molecule("exported").unwrap();
    assert!(TrackedDuplex::from_input(&exported, &full()).is_err());
    let treated = source.treat(EndTreatment::Phosphorylate).unwrap();
    assert_eq!(source.top_modifications, treated.top_modifications);
    assert_eq!(source.bottom_modifications, treated.bottom_modifications);
}
#[test]
fn checked_full_cpg_blocked_circle_and_mixed_or_unknown_sites_withhold_products() {
    let source = tracked("AAGCGGCCGCAA", Topology::Circular);
    let o = CheckedDigestOptions {
        selection: checked_selection("NotI").unwrap(),
        ..Default::default()
    };
    let r = checked_digest(&source, &o).unwrap();
    assert!(r.complete && r.products_determined);
    assert_eq!(r.sites[0].outcome, Outcome::Blocked);
    assert_eq!(r.fragments.len(), 1);
    assert_eq!(r.fragments[0].duplex.duplex.topology, Topology::Circular);
    replay_checked(&source, &o, &r).unwrap();
    let unknown = TrackedDuplex::from_unknown_duplex(&source.duplex).unwrap();
    let r = checked_digest(&unknown, &o).unwrap();
    assert!(r.complete && !r.products_determined && r.fragments.is_empty());
    let mut partial = source.clone();
    partial.bottom_modifications = StrandModifications {
        n6a_default: CallState::Unmodified,
        c5c_default: CallState::Unmodified,
        spans: vec![],
    };
    assert_eq!(
        checked_digest(&partial, &o).unwrap().sites[0].outcome,
        Outcome::Unknown
    );
}
#[test]
fn checked_comparison_bounds_and_stale_catalog_do_not_produce_a_complete_product_set() {
    let source = tracked("AAGAATTCCGAATTCGCGGCCGC", Topology::Circular);
    let mut o = CheckedDigestOptions {
        selection: checked_selection("EcoRI").unwrap(),
        ..Default::default()
    };
    o.max_comparisons = 1;
    let r = checked_digest(&source, &o).unwrap();
    assert!(!r.complete && !r.products_determined && r.fragments.is_empty());
    o.selection.catalog_content_hash = "0".repeat(64);
    assert!(checked_digest(&source, &o).is_err());
}
#[test]
fn unknown_legacy_wrapper_retains_large_original_capacity_without_dense_state_arrays() {
    let d = chemistry::Duplex::from_molecule(
        &Molecule::dna("large", &"A".repeat(1_000_001), Topology::Linear).unwrap(),
    )
    .unwrap();
    let unknown = TrackedDuplex::from_unknown_duplex(&d).unwrap();
    assert!(
        unknown.top_modifications.spans.is_empty() && unknown.bottom_modifications.spans.is_empty()
    );
    assert_eq!(summary(&unknown).unwrap().unknown_eligible, 1_000_001);
    assert!(matches!(
        checked_digest(&unknown, &CheckedDigestOptions::default()),
        Err(Error::TooLarge)
    ));
}
#[test]
fn inherited_methyladenine_survives_a_broken_motif_and_single_stranded_sites_cannot_permit() {
    let source = tracked("AAGATCTAGACC", Topology::Linear);
    let cutter = analysis::Enzyme {
        name: "fixture blunt boundary".into(),
        recognition: "AAGA".into(),
        top_cut: 4,
        bottom_cut: 4,
    };
    let fragments = digest(&source, &[cutter]).unwrap();
    let product = &fragments[1].duplex;
    assert_eq!(product.duplex.top.as_str(), "TCTAGACC");
    assert!(!product.duplex.top.as_str().contains("GATC"));
    assert_eq!(
        product.state_at(Strand::Bottom, Chemical::N6A, 7).unwrap(),
        Some(CallState::Modified)
    );
    let o = CheckedDigestOptions {
        selection: checked_selection("XbaI").unwrap(),
        ..Default::default()
    };
    let report = checked_digest(product, &o).unwrap();
    assert_eq!(report.sites[0].outcome, Outcome::Unknown);
    assert!(report.fragments.is_empty() && !report.products_determined);
    let mut malformed = chemistry::DigestFragment {
        duplex: product.duplex.clone(),
        top_start: 0,
        bottom_start: i64::MIN,
    };
    assert!(map_digest_fragment(product, &malformed).is_err());
    malformed.bottom_start = i64::MAX;
    assert!(map_digest_fragment(product, &malformed).is_err());
    let d = chemistry::Duplex {
        top: Sequence::dna("TCTAGACC").unwrap(),
        bottom: Sequence::dna("GGTC").unwrap(),
        bottom_start: 4,
        topology: Topology::Linear,
        left: chemistry::Terminus::PCR,
        right: chemistry::Terminus::PCR,
        nicks: vec![],
    };
    let unmodified = TrackedDuplex::annealed_unmodified(&d).unwrap();
    let report = checked_digest(&unmodified, &o).unwrap();
    assert_eq!(report.sites[0].outcome, Outcome::Unknown);
    assert!(report.fragments.is_empty());
}
