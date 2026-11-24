use collagene_core::view_tracks::*;
use collagene_core::{Feature, Location, Molecule, Topology};
use std::sync::Arc;
fn options() -> TrackOptions {
    TrackOptions {
        features: false,
        primers: false,
        enzyme_names: vec![],
        ..Default::default()
    }
}
fn feature(kind: &str, location: Location) -> Feature {
    Feature {
        kind: kind.into(),
        location,
        qualifiers: vec![],
        needs_review: false,
    }
}
fn report(m: &Molecule, start: usize, end: usize, o: &TrackOptions) -> TrackReport {
    let r = view_tracks(m, TrackWindow { start, end }, o).unwrap();
    r.validate().unwrap();
    r
}
#[test]
fn three_frames_keep_reference_phase_across_line_boundaries() {
    // ATG AAA CCC GGG TAA -> M K P G *; offset-one TGA AAC CCG GGT -> * N P G.
    let m = Molecule::dna("phase", "ATGAAACCCGGGTAA", Topology::Linear).unwrap();
    let o = TrackOptions {
        translations: true,
        orientation: TrackOrientation::Forward,
        ..options()
    };
    let r = report(&m, 4, 11, &o);
    let f = &r.items.iter().find(|i| i.frame == Some(0)).unwrap();
    assert_eq!(f.protein.as_deref(), Some("KPG"));
    assert_eq!(
        f.residues
            .iter()
            .map(|r| r.residue_index)
            .collect::<Vec<_>>(),
        [1, 2, 3]
    );
    assert_eq!(
        f.residues[0].codon_parts,
        [ReferencePart { start: 3, end: 6 }]
    );
    assert_eq!(
        f.residues[2].codon_parts,
        [ReferencePart { start: 9, end: 12 }]
    );
    let one = r.items.iter().find(|i| i.frame == Some(1)).unwrap();
    assert_eq!(one.protein.as_deref(), Some("NPG"));
    let left = report(&m, 0, 5, &o);
    let right = report(&m, 5, 10, &o);
    let l = left.items.iter().find(|i| i.frame == Some(0)).unwrap();
    let r = right.items.iter().find(|i| i.frame == Some(0)).unwrap();
    assert_eq!(l.residues[1], r.residues[0]);
}
#[test]
fn reverse_frame_uses_full_reference_end_and_reverse_codon_mapping() {
    // Reverse complement is ATG AAA CCC GGG TAA.
    let m = Molecule::dna("reverse", "TTACCCGGGTTTCAT", Topology::Linear).unwrap();
    let o = TrackOptions {
        translations: true,
        frames: vec![0],
        orientation: TrackOrientation::Reverse,
        ..options()
    };
    let r = report(&m, 8, 13, &o);
    assert_eq!(r.items.len(), 1);
    let t = &r.items[0];
    assert_eq!(t.protein.as_deref(), Some("MKP"));
    assert_eq!(
        t.residues[0].codon_parts,
        [ReferencePart { start: 12, end: 15 }]
    );
    assert_eq!(
        t.residues[1].codon_parts,
        [ReferencePart { start: 9, end: 12 }]
    );
    assert_eq!(
        t.residues[2].codon_parts,
        [ReferencePart { start: 6, end: 9 }]
    );
}
#[test]
fn circular_origin_codon_is_same_in_both_viewports_and_has_split_parts() {
    // Frame zero codons: ATG / CAA / AAT (last codon wraps from reference 6).
    let m = Molecule::dna("origin", "ATGCAAA", Topology::Circular).unwrap();
    let o = TrackOptions {
        translations: true,
        frames: vec![0],
        orientation: TrackOrientation::Forward,
        ..options()
    };
    let r = report(&m, 6, 9, &o);
    let t = &r.items[0];
    assert_eq!(t.protein.as_deref(), Some("MN"));
    let last = &t.residues[1];
    assert_eq!(last.residue_index, 2);
    assert_eq!(last.amino_acid, "N");
    assert_eq!(
        last.codon_parts,
        [
            ReferencePart { start: 6, end: 7 },
            ReferencePart { start: 0, end: 2 }
        ]
    );
    assert_eq!(last, &report(&m, 0, 1, &o).items[0].residues[1]);
}
#[test]
fn six_frames_and_rna_use_embedded_translation() {
    let mut m = Molecule::dna("RNA", "ATGAAATAA", Topology::Linear).unwrap();
    m.sequence = collagene_core::Sequence::new(collagene_core::Alphabet::Rna, "AUGAAAUAA").unwrap();
    let o = TrackOptions {
        translations: true,
        ..options()
    };
    let r = report(&m, 0, 9, &o);
    assert_eq!(r.items.len(), 6);
    assert_eq!(r.items[0].protein.as_deref(), Some("MK*"));
    assert_eq!(
        r.items.iter().filter(|i| i.reverse == Some(true)).count(),
        3
    );
}
#[test]
fn primer_bind_is_independently_toggleable_and_reviewed_features_have_no_arrow() {
    let mut m = Molecule::dna("annotations", "ATGAAATAA", Topology::Linear).unwrap();
    m.features = vec![
        feature("gene", Location::range(1, 8)),
        feature("primer_bind", Location::range(2, 6)),
        feature(
            "misc_feature",
            Location::Order(vec![Location::range(0, 3), Location::range(6, 9)]),
        ),
    ];
    let r = report(
        &m,
        0,
        9,
        &TrackOptions {
            primers: true,
            ..options()
        },
    );
    assert_eq!(r.items.len(), 1);
    assert_eq!(r.items[0].kind, TrackKind::Primer);
    assert_eq!(r.items[0].source_feature, Some(1));
    let r = report(
        &m,
        0,
        9,
        &TrackOptions {
            features: true,
            ..options()
        },
    );
    assert_eq!(r.items.len(), 2);
    assert_eq!(r.items[1].reverse, None);
}
#[test]
fn type_iis_recognition_and_actual_strand_cuts_are_separate() {
    let m = Molecule::dna("BsaI", "AAGGTCTCAAAAAAAA", Topology::Linear).unwrap();
    let o = TrackOptions {
        enzyme_recognition: true,
        enzyme_cuts: true,
        enzyme_names: vec!["BsaI".into()],
        ..options()
    };
    let r = report(&m, 0, 16, &o);
    assert_eq!(r.items.len(), 1);
    assert_eq!(r.items[0].parts, [ReferencePart { start: 2, end: 8 }]);
    assert_eq!(
        r.items[0].cuts,
        Some(CutPositions {
            top: Some(9),
            bottom: Some(13)
        })
    );
    let cuts = report(
        &m,
        12,
        14,
        &TrackOptions {
            enzyme_recognition: false,
            ..o.clone()
        },
    );
    assert_eq!(cuts.items.len(), 1);
    let recognition = report(
        &m,
        12,
        14,
        &TrackOptions {
            enzyme_cuts: false,
            ..o
        },
    );
    assert!(recognition.items.is_empty());
}
#[test]
fn reverse_type_iis_and_circular_cut_coordinates_have_literal_oracles() {
    let m = Molecule::dna("reverse BsaI", "AAAAAAAAAAAAGAGACC", Topology::Linear).unwrap();
    let o = TrackOptions {
        enzyme_recognition: true,
        enzyme_cuts: true,
        enzyme_names: vec!["BsaI".into()],
        ..options()
    };
    let r = report(&m, 0, 18, &o);
    assert_eq!(r.items[0].reverse, Some(true));
    assert_eq!(
        r.items[0].cuts,
        Some(CutPositions {
            top: Some(7),
            bottom: Some(11)
        })
    );
    let m = Molecule::dna("circular BsaI", "GGTCTCAA", Topology::Circular).unwrap();
    let r = report(&m, 0, 8, &o);
    assert_eq!(
        r.items[0].cuts,
        Some(CutPositions {
            top: Some(7),
            bottom: Some(3)
        })
    );
}
#[test]
fn known_cds_reverse_join_and_codon_start_retain_reference_mapping() {
    // Biological CDS is ATG AAA TAA; reference pieces are reverse-complemented.
    let mut m = Molecule::dna("joined CDS", "TTATTTGGGCAT", Topology::Linear).unwrap();
    let mut cds = feature(
        "CDS",
        Location::Complement(Box::new(Location::Join(vec![
            Location::range(0, 6),
            Location::range(9, 12),
        ]))),
    );
    cds.qualifiers = vec![("label".into(), Some("joined".into()))];
    m.features.push(cds);
    let o = TrackOptions {
        known_cds: KnownCdsPolicy::Validated,
        ..options()
    };
    let r = report(&m, 3, 11, &o);
    assert_eq!(r.items[0].protein.as_deref(), Some("MK"));
    assert_eq!(
        r.items[0].residues[0].codon_parts,
        [ReferencePart { start: 9, end: 12 }]
    );
    assert_eq!(
        r.items[0].residues[1].codon_parts,
        [ReferencePart { start: 3, end: 6 }]
    );
    let mut m = Molecule::dna("phase CDS", "AATGAAATAA", Topology::Linear).unwrap();
    let mut cds = feature("CDS", Location::range(0, 10));
    cds.qualifiers
        .push(("codon_start".into(), Some("2".into())));
    m.features.push(cds);
    let r = report(&m, 3, 5, &o);
    assert_eq!(r.items[0].protein.as_deref(), Some("MK"));
    assert_eq!(
        r.items[0].residues[0].codon_parts,
        [ReferencePart { start: 1, end: 4 }]
    );
    assert_eq!(
        r.items[0].residues[1].codon_parts,
        [ReferencePart { start: 4, end: 7 }]
    );
}
#[test]
fn fuzzy_unordered_and_review_required_cds_never_invent_a_frame() {
    let mut m = Molecule::dna("unsupported CDS", "ATGAAATAA", Topology::Linear).unwrap();
    let mut reviewed = feature("CDS", Location::range(0, 9));
    reviewed.needs_review = true;
    m.features = vec![
        reviewed,
        feature("CDS", Location::Order(vec![Location::range(0, 9)])),
        feature(
            "CDS",
            Location::Range {
                start: 0,
                end: 9,
                before: true,
                after: false,
            },
        ),
    ];
    let r = report(
        &m,
        0,
        9,
        &TrackOptions {
            known_cds: KnownCdsPolicy::Validated,
            ..options()
        },
    );
    assert!(!r.complete);
    assert!(r.items.is_empty());
    assert!(r.warnings.iter().any(|w| w.contains("requires review")));
}
#[test]
fn orf_prediction_uses_whole_source_and_caps_before_filtering() {
    let m = Molecule::dna("ORFs", "ATGAAATAAATGCCCTAA", Topology::Linear).unwrap();
    let o = TrackOptions {
        orfs: true,
        min_orf_amino_acids: 1,
        frames: vec![0],
        orientation: TrackOrientation::Forward,
        ..options()
    };
    let r = report(&m, 4, 5, &o);
    assert_eq!(r.items.len(), 1);
    assert_eq!(r.items[0].kind, TrackKind::Orf);
    assert_eq!(r.items[0].protein.as_deref(), Some("K"));
    assert_eq!(r.items[0].residues[0].residue_index, 1);
    assert_eq!(r.items[0].parts, [ReferencePart { start: 0, end: 9 }]);
    let r = report(
        &m,
        0,
        18,
        &TrackOptions {
            max_items: 1,
            ..o.clone()
        },
    );
    assert!(!r.complete);
    assert_eq!(r.items.len(), 1);
    let r = report(
        &m,
        0,
        18,
        &TrackOptions {
            max_work_units: 10,
            ..o
        },
    );
    assert!(!r.complete);
    assert!(r.items.is_empty());
    assert!(r.work.used_units <= 10);
}
#[test]
fn caps_work_items_residues_and_immutable_revision_are_enforced() {
    let mut m = Molecule::dna("caps", "ATGAAATAA", Topology::Linear).unwrap();
    m.features = (0..10)
        .map(|_| feature("gene", Location::range(0, 9)))
        .collect();
    let r = report(
        &m,
        0,
        9,
        &TrackOptions {
            features: true,
            max_items: 2,
            ..options()
        },
    );
    assert_eq!(r.items.len(), 2);
    assert!(!r.complete);
    let r = report(
        &m,
        0,
        9,
        &TrackOptions {
            translations: true,
            max_residues: 2,
            ..options()
        },
    );
    assert_eq!(r.work.emitted_residues, 2);
    assert!(!r.complete);
    let r = report(
        &m,
        0,
        9,
        &TrackOptions {
            translations: true,
            max_items: 1,
            ..options()
        },
    );
    assert_eq!(r.items.len(), 1);
    assert!(!r.complete);
    let prepared = PreparedTracks::new(Arc::new(m.clone())).unwrap();
    let revision = prepared.source_revision().to_string();
    m.features.clear();
    assert_ne!(m.revision(), revision);
    assert!(
        prepared
            .report_for_revision(&m.revision(), TrackWindow { start: 0, end: 9 }, &options())
            .is_err()
    );
    assert_eq!(
        prepared
            .report(TrackWindow { start: 0, end: 9 }, &options())
            .unwrap()
            .source_revision,
        revision
    );
}
#[test]
fn invalid_configuration_corrupt_source_and_forged_report_are_rejected() {
    let m = Molecule::dna("validation", "ATGAAATAA", Topology::Linear).unwrap();
    for o in [
        TrackOptions {
            frames: vec![3],
            ..options()
        },
        TrackOptions {
            frames: vec![0, 0],
            ..options()
        },
        TrackOptions {
            genetic_code: 99,
            ..options()
        },
        TrackOptions {
            enzyme_names: vec!["unknown".into()],
            ..options()
        },
        TrackOptions {
            max_items: 1001,
            ..options()
        },
        TrackOptions {
            max_residues: 10001,
            ..options()
        },
        TrackOptions {
            max_work_units: 0,
            ..options()
        },
    ] {
        assert!(view_tracks(&m, TrackWindow { start: 0, end: 9 }, &o).is_err());
    }
    assert!(view_tracks(&m, TrackWindow { start: 8, end: 10 }, &options()).is_err());
    let mut corrupt = m.clone();
    corrupt
        .features
        .push(feature("gene", Location::range(8, 99)));
    assert!(PreparedTracks::new(Arc::new(corrupt)).is_err());
    let mut r = report(
        &m,
        0,
        9,
        &TrackOptions {
            translations: true,
            ..options()
        },
    );
    r.items[0].residues[0].codon_parts[0].end = 99;
    assert!(r.validate().is_err());
}
#[test]
fn long_sources_skip_whole_source_orf_work_but_translate_only_visible_bases() {
    let m = Molecule::dna("large", &"ATG".repeat(100000), Topology::Linear).unwrap();
    let o = TrackOptions {
        translations: true,
        orfs: true,
        max_work_units: 1000,
        frames: vec![0],
        orientation: TrackOrientation::Forward,
        ..options()
    };
    let r = report(&m, 150000, 150030, &o);
    assert_eq!(r.items.len(), 1);
    assert_eq!(r.items[0].protein.as_deref(), Some("MMMMMMMMMM"));
    assert_eq!(r.items[0].residues[0].residue_index, 50000);
    assert!(r.work.used_units <= 1000);
    assert!(!r.complete);
}
#[test]
fn reverse_circular_origin_and_join_split_codon_have_literal_parts() {
    // Reverse complement is ATGCAAA; final frame-zero codon AAT maps reference 0,6,5.
    let m = Molecule::dna("reverse origin", "TTTGCAT", Topology::Circular).unwrap();
    let o = TrackOptions {
        translations: true,
        frames: vec![0],
        orientation: TrackOrientation::Reverse,
        ..options()
    };
    let r = report(&m, 5, 8, &o);
    assert_eq!(r.items[0].protein.as_deref(), Some("MN"));
    assert_eq!(
        r.items[0].residues[1].codon_parts,
        [
            ReferencePart { start: 0, end: 1 },
            ReferencePart { start: 5, end: 7 }
        ]
    );
    let mut m = Molecule::dna("split codon", "ATCCCGAAATAA", Topology::Linear).unwrap();
    m.features.push(feature(
        "CDS",
        Location::Join(vec![Location::range(0, 2), Location::range(5, 12)]),
    ));
    let r = report(
        &m,
        0,
        2,
        &TrackOptions {
            known_cds: KnownCdsPolicy::Validated,
            ..options()
        },
    );
    assert_eq!(r.items[0].protein.as_deref(), Some("M"));
    assert_eq!(
        r.items[0].residues[0].codon_parts,
        [
            ReferencePart { start: 0, end: 2 },
            ReferencePart { start: 5, end: 6 }
        ]
    );
}
#[test]
fn many_compound_parts_remain_bounded_and_window_lookup_uses_cached_metadata() {
    let mut m = Molecule::dna("many parts", &"AT".repeat(10002), Topology::Linear).unwrap();
    m.features.push(feature(
        "misc_feature",
        Location::Join(
            (0..10001)
                .map(|i| Location::range(2 * i, 2 * i + 1))
                .collect(),
        ),
    ));
    let r = report(
        &m,
        0,
        20002,
        &TrackOptions {
            features: true,
            max_window_bases: 30000,
            max_work_units: 50000,
            ..options()
        },
    );
    assert!(!r.complete);
    assert_eq!(r.items[0].parts.len(), MAX_TRACK_PARTS);
    let r = report(
        &m,
        18000,
        18001,
        &TrackOptions {
            features: true,
            max_work_units: 2,
            ..options()
        },
    );
    assert_eq!(
        r.items[0].parts,
        [ReferencePart {
            start: 18000,
            end: 18001
        }]
    );
    assert_eq!(r.work.used_units, 1);
}
#[test]
fn unplaced_annotations_are_disclosed_only_for_affected_enabled_tracks() {
    let mut m = Molecule::dna("unplaced", "ATGAAATAA", Topology::Linear).unwrap();
    m.features
        .push(feature("CDS", Location::Unresolved("remote:1..9".into())));
    let r = report(
        &m,
        0,
        9,
        &TrackOptions {
            known_cds: KnownCdsPolicy::Validated,
            ..options()
        },
    );
    assert!(!r.complete);
    assert!(r.items.is_empty());
    assert!(
        r.warnings
            .iter()
            .any(|w| w.contains("unresolved locations"))
    );
    let r = report(&m, 0, 9, &options());
    assert!(r.complete);
    assert!(r.items.is_empty());
}
