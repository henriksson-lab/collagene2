//! Independent review regressions: ambiguity/partial certainty and circular coordinate oracles.
use collagene_core::view_tracks::*;
use collagene_model::{Feature, Location, Molecule, Topology};

fn options() -> TrackOptions {
    TrackOptions {
        features: false,
        primers: false,
        enzyme_names: vec![],
        frames: vec![0],
        orientation: TrackOrientation::Forward,
        ..Default::default()
    }
}
fn report(m: &Molecule, window: TrackWindow, options: TrackOptions) -> TrackReport {
    let report = view_tracks(m, window, &options).unwrap();
    report.validate().unwrap();
    report
}
#[test]
fn unresolved_frame_codon_is_not_certain_and_possible_stop_is_visible() {
    // TAN can be Tyr or stop. A frame display must not upgrade X to certainty.
    let m = Molecule::dna("ambiguous", "ATGTANAAA", Topology::Linear).unwrap();
    let r = report(
        &m,
        TrackWindow { start: 0, end: 9 },
        TrackOptions {
            translations: true,
            ..options()
        },
    );
    let t = r
        .items
        .iter()
        .find(|t| t.kind == TrackKind::Translation)
        .unwrap();
    assert_eq!(t.protein.as_deref(), Some("MXK"));
    assert!(!t.certain);
    assert!(
        r.warnings.iter().any(|w| w.contains("stop")),
        "possible-stop warning must survive viewport projection"
    );
}
#[test]
fn partial_and_ambiguous_orfs_are_not_certain() {
    for sequence in ["ATGAANTAA", "ATGGCNTAA"] {
        let m = Molecule::dna("candidate", sequence, Topology::Linear).unwrap();
        let r = report(
            &m,
            TrackWindow {
                start: 0,
                end: sequence.len(),
            },
            TrackOptions {
                orfs: true,
                min_orf_amino_acids: 1,
                ..options()
            },
        );
        assert!(
            !r.items.is_empty(),
            "{sequence} should retain model candidate"
        );
        assert!(
            r.items
                .iter()
                .filter(|t| t.kind == TrackKind::Orf)
                .all(|t| !t.certain),
            "{sequence}: {:#?}",
            r.items
        );
    }
    for sequence in ["ATGAAA", "ATGAAATAN"] {
        let m = Molecule::dna("partial", sequence, Topology::Linear).unwrap();
        let r = report(
            &m,
            TrackWindow {
                start: 0,
                end: sequence.len(),
            },
            TrackOptions {
                orfs: true,
                min_orf_amino_acids: 1,
                ..options()
            },
        );
        assert!(
            r.items.iter().all(|t| t.kind != TrackKind::Orf),
            "CompleteOnly must not upgrade partial/uncertain-stop candidates"
        );
    }
    let m = Molecule::dna("complete", "ATGAAATAA", Topology::Linear).unwrap();
    let r = report(
        &m,
        TrackWindow { start: 0, end: 9 },
        TrackOptions {
            orfs: true,
            min_orf_amino_acids: 1,
            ..options()
        },
    );
    assert!(
        r.items
            .iter()
            .any(|t| t.kind == TrackKind::Orf && t.certain)
    );
}
#[test]
fn viewport_display_labels_cannot_clone_unbounded_source_qualifiers() {
    let mut m = Molecule::dna("label", "ACGT", Topology::Linear).unwrap();
    m.features.push(Feature {
        kind: "gene".into(),
        location: Location::range(0, 4),
        qualifiers: vec![("label".into(), Some("a".repeat(100_000)))],
        needs_review: false,
    });
    let r = report(
        &m,
        TrackWindow { start: 0, end: 4 },
        TrackOptions {
            features: true,
            ..options()
        },
    );
    assert_eq!(r.items.len(), 1);
    assert!(
        r.items[0].label.len() <= 4096,
        "display label exceeded budget"
    );
    assert_eq!(
        m.features[0].label().len(),
        100_000,
        "source label must remain intact"
    );
}
#[test]
fn all_six_circular_frames_have_exact_source_codons_in_each_single_base_window() {
    let m = Molecule::dna("circular", "ATGCAAA", Topology::Circular).unwrap();
    // Manually decoded from ATGCAAA and its reverse complement TTTGCAT;
    // final codons wrap one source origin but never establish a global circular phase.
    let amino = [
        ["M", "Q", "N"],
        ["C", "K", ""],
        ["A", "K", ""],
        ["F", "A", "F"],
        ["L", "H", ""],
        ["C", "I", ""],
    ];
    let n = 7;
    for reverse in [false, true] {
        for frame in 0..3 {
            for visible in 0..n {
                let r = report(
                    &m,
                    TrackWindow {
                        start: visible,
                        end: visible + 1,
                    },
                    TrackOptions {
                        translations: true,
                        frames: vec![frame as u8],
                        orientation: if reverse {
                            TrackOrientation::Reverse
                        } else {
                            TrackOrientation::Forward
                        },
                        ..options()
                    },
                );
                let expected: Vec<_> = (frame..n)
                    .step_by(3)
                    .enumerate()
                    .filter_map(|(index, start)| {
                        let positions: Vec<_> = (start..start + 3)
                            .map(|p| if reverse { n - 1 - p % n } else { p % n })
                            .collect();
                        positions.contains(&visible).then_some((
                            index,
                            positions,
                            amino[usize::from(reverse) * 3 + frame][index],
                        ))
                    })
                    .collect();
                let actual: Vec<_> = r.items.iter().flat_map(|t| t.residues.iter()).collect();
                assert_eq!(
                    actual.len(),
                    expected.len(),
                    "reverse={reverse} frame={frame} window={visible}"
                );
                for (actual, (index, positions, aa)) in actual.into_iter().zip(expected) {
                    assert_eq!(actual.residue_index, index);
                    assert_eq!(actual.amino_acid, aa);
                    let decoded: Vec<_> = actual
                        .codon_parts
                        .iter()
                        .flat_map(|p| {
                            if reverse {
                                (p.start..p.end).rev().collect::<Vec<_>>()
                            } else {
                                (p.start..p.end).collect()
                            }
                        })
                        .collect();
                    assert_eq!(
                        decoded, positions,
                        "reverse={reverse} frame={frame} window={visible}"
                    );
                }
            }
        }
    }
}

#[test]
fn palindromic_recognition_is_retained_once_under_reverse_orientation() {
    let m = Molecule::dna("palindrome", "CCGAATTCCC", Topology::Linear).unwrap();
    for orientation in [
        TrackOrientation::Forward,
        TrackOrientation::Reverse,
        TrackOrientation::Both,
    ] {
        let r = report(
            &m,
            TrackWindow { start: 0, end: 10 },
            TrackOptions {
                enzyme_recognition: true,
                enzyme_cuts: true,
                enzyme_names: vec!["EcoRI".into()],
                orientation,
                ..options()
            },
        );
        let sites: Vec<_> = r
            .items
            .iter()
            .filter(|t| t.kind == TrackKind::Enzyme)
            .collect();
        assert_eq!(
            sites.len(),
            1,
            "palindromic site must not disappear or duplicate under {orientation:?}"
        );
        assert_eq!(sites[0].parts, vec![ReferencePart { start: 2, end: 8 }]);
        assert_eq!(
            sites[0].cuts.as_ref().unwrap(),
            &CutPositions {
                top: Some(3),
                bottom: Some(7)
            }
        );
    }
}
