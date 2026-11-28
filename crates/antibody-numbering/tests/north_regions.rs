//! Literal source slices for the North2011 region definition. These authored
//! rearrangements exercise conditional reference placement, not natural-domain
//! numbering ground truth. Expected loops do not use the region implementation.
use collagene_antibody_numbering::{
    Insertion, Position, domain_candidates as d,
    domain_tracks::*,
    north, profile,
    v_mapping::{self as v, Locus},
};
use collagene_model::{Alphabet, Molecule, Sequence, Topology};

const HEAVY: &str = "EVQLLESGGGLVQPGGSLRLSCAASGFTFSSYAMSWVRQAPGKGLEWVSAISGSGGSTYYADSVKGRFTISRDNSKNTLYLQMNSLRAEDTAVYYC";
const KAPPA: &str =
    "DIQMTQSPSSLSASVGDRVTITCRASQSISSYLNWYQQKPGKAPKLLIYAASSLQSGVPSRFSGSGSGTDFTLTISSLQPEDFATYYC";
const LAMBDA: &str =
    "QSALTQPASVSGSPGQSITISCTGTSSDVGGYNYVSWYQQHPGKAPKLMIYEVSNRPSGVSNRFSGSKSGNTASLTISGLQAEDEADYYC";
const H3: &str = "ARDYGMDVAAADA";
const HJ: &str = "WGQGTLVTVSS";

fn source(sequence: &str) -> Molecule {
    Molecule {
        name: "Authored North2011 region fixture".into(),
        sequence: Sequence::new(Alphabet::Protein, sequence).unwrap(),
        topology: Topology::Linear,
        features: vec![],
        metadata: Default::default(),
        original_genbank: None,
    }
}
fn options() -> d::Options {
    d::Options {
        v_scoring: v::Scoring {
            observed_reward: 5,
            frequency_bonus: 8,
            unseen_penalty: -6,
            gap_open: -12,
            gap_extend: -2,
            minimum_score: 300,
        },
        min_cdr3_residues: 2,
        max_cdr3_residues: 60,
        max_fr4_mismatches: 1,
    }
}
fn limits() -> d::Limits {
    d::Limits {
        v_mapping: v::Limits {
            max_cells: profile::MAX_CELLS,
            max_endpoints_per_locus: profile::MAX_ENDPOINTS,
            max_trace_states: profile::MAX_TRACE_STATES,
        },
        max_scan_windows: d::MAX_SCAN_WINDOWS,
        max_candidates: 32,
    }
}
fn selection(numbering_scheme: NumberingScheme) -> Selection {
    Selection {
        numbering_scheme,
        region_definition: RegionDefinition::North2011,
    }
}
fn row(track: &Track, index: usize) -> &Residue {
    track
        .residues
        .iter()
        .find(|r| r.source_index == index)
        .unwrap()
}
fn assert_loop(track: &Track, kind: Region, start: usize, end: usize, expected: &str) {
    let cdr = track.cdrs.iter().find(|c| c.kind == kind).unwrap();
    assert_eq!(cdr.status, CdrStatus::Mapped { start, end });
    let residues = track
        .residues
        .iter()
        .filter(|r| r.source_index >= start && r.source_index < end)
        .collect::<Vec<_>>();
    assert_eq!(residues.len(), end - start);
    assert!(residues.iter().all(|r| r.region == Some(kind.clone())));
    assert_eq!(
        residues.iter().map(|r| r.residue).collect::<String>(),
        expected
    );
}

#[test]
fn north_membership_uses_established_imgt_labels_and_includes_the_extra_light_second_loop_residue()
{
    // Primary North boundaries expressed through the IMGT2003 correspondence
    // table: first24..40, second55..66(H)/69(L), third105..117, inclusive.
    for locus in [Locus::Heavy, Locus::Kappa, Locus::Lambda] {
        let second_end = if locus == Locus::Heavy { 66 } else { 69 };
        for (range, expected) in [
            (1..=23, Region::Fr1),
            (24..=40, Region::Cdr1),
            (41..=54, Region::Fr2),
            (55..=second_end, Region::Cdr2),
            (second_end + 1..=104, Region::Fr3),
            (105..=117, Region::Cdr3),
            (118..=128, Region::Fr4),
        ] {
            for number in range {
                assert_eq!(
                    north::region(&locus, &number.to_string().parse().unwrap()).unwrap(),
                    expected,
                    "{locus:?} {number}"
                );
            }
        }
        for text in ["111.1", "112.1", "111.10", "112.10"] {
            let label: Position = text.parse().unwrap();
            assert_eq!(north::region(&locus, &label).unwrap(), Region::Cdr3);
            assert_eq!(label.to_string(), text);
        }
        for text in ["129", "104.1", "111A", "118.1"] {
            assert!(north::region(&locus, &text.parse().unwrap()).is_err());
        }
    }
    assert!(serde_json::from_str::<RegionDefinition>("\"North\"").is_err());
}

#[test]
fn heavy_north_loops_have_literal_source_intervals_under_each_independent_numbering_choice() {
    let molecule = source(&format!("{HEAVY}{H3}{HJ}"));
    let original = molecule.clone();
    for scheme in [
        NumberingScheme::Imgt,
        NumberingScheme::Kabat1991,
        NumberingScheme::Chothia1997,
    ] {
        let choice = selection(scheme.clone());
        let report = analyze(&molecule, &options(), &limits(), &choice, &mut || false).unwrap();
        let track = report
            .tracks
            .iter()
            .find(|t| t.locus == Locus::Heavy)
            .unwrap();
        assert_loop(track, Region::Cdr1, 22, 35, "AASGFTFSSYAMS");
        assert_loop(track, Region::Cdr2, 49, 59, "AISGSGGSTY");
        assert_loop(track, Region::Cdr3, 96, 109, "ARDYGMDVAAADA");
        for (index, expected) in [
            (21, Region::Fr1),
            (35, Region::Fr2),
            (48, Region::Fr2),
            (59, Region::Fr3),
            (95, Region::Fr3),
            (109, Region::Fr4),
        ] {
            assert_eq!(row(track, index).region, Some(expected));
        }
        let baseline = analyze(
            &molecule,
            &options(),
            &limits(),
            &Selection {
                numbering_scheme: scheme,
                region_definition: RegionDefinition::Imgt,
            },
            &mut || false,
        )
        .unwrap();
        assert_ne!(report.context_revision, baseline.context_revision);
        assert_eq!(report.evidence, baseline.evidence);
        assert_eq!(
            track
                .residues
                .iter()
                .map(|r| &r.assignment)
                .collect::<Vec<_>>(),
            baseline.tracks[0]
                .residues
                .iter()
                .map(|r| &r.assignment)
                .collect::<Vec<_>>()
        );
        report
            .validate_for(&molecule, &options(), &limits(), &choice)
            .unwrap();
        report.validate_binding(&report.context_revision).unwrap();
    }
    assert_eq!(molecule, original);
}

#[test]
fn light_north_loops_keep_literal_extra_y_and_locus_specific_source_offsets() {
    for (prefix, tail, locus, first, second, third) in [
        (
            KAPPA,
            "FGQGTKVEIK",
            Locus::Kappa,
            (23, 34, "RASQSISSYLN"),
            (48, 56, "YAASSLQS"),
            (88, 97),
        ),
        (
            LAMBDA,
            "FGGGTKLTVL",
            Locus::Lambda,
            (22, 36, "TGTSSDVGGYNYVS"),
            (50, 58, "YEVSNRPS"),
            (90, 99),
        ),
    ] {
        let molecule = source(&format!("{prefix}QQSYSTPPT{tail}"));
        for scheme in [
            NumberingScheme::Imgt,
            NumberingScheme::Kabat1991,
            NumberingScheme::Chothia1997,
        ] {
            let choice = selection(scheme);
            let report = analyze(&molecule, &options(), &limits(), &choice, &mut || false).unwrap();
            let track = report.tracks.iter().find(|t| t.locus == locus).unwrap();
            assert_loop(track, Region::Cdr1, first.0, first.1, first.2);
            assert_loop(track, Region::Cdr2, second.0, second.1, second.2);
            assert_loop(track, Region::Cdr3, third.0, third.1, "QQSYSTPPT");
            assert_eq!(row(track, second.0 - 1).region, Some(Region::Fr2));
            assert_eq!(row(track, second.0).residue, 'Y');
            assert_eq!(row(track, second.1).region, Some(Region::Fr3));
            assert_eq!(track.unassigned_regions, 0);
            report
                .validate_for(&molecule, &options(), &limits(), &choice)
                .unwrap();
        }
    }
}

#[test]
fn proven_outer_flanks_keep_north_regions_without_inventing_unsupported_interior_labels() {
    // Eight-residue H3 has established IMGT labels but no supported alternate
    // numbering window. North's whole physical loop does not depend on it.
    let short = source(&format!("{HEAVY}ARDYGMDV{HJ}"));
    for scheme in [
        NumberingScheme::Imgt,
        NumberingScheme::Kabat1991,
        NumberingScheme::Chothia1997,
    ] {
        let choice = selection(scheme.clone());
        let report = analyze(&short, &options(), &limits(), &choice, &mut || false).unwrap();
        let track = report
            .tracks
            .iter()
            .find(|t| t.locus == Locus::Heavy)
            .unwrap();
        assert_loop(track, Region::Cdr3, 96, 104, "ARDYGMDV");
        assert_eq!(track.unassigned_regions, 0);
        if scheme != NumberingScheme::Imgt {
            assert_eq!(track.unassigned_residues, 8);
            assert!((96..104).all(|i| row(track, i).assignment.is_none()));
        }
        report
            .validate_for(&short, &options(), &limits(), &choice)
            .unwrap();
    }
    let long = source(&format!(
        "{}{}{}{H3}{HJ}",
        &HEAVY[..27],
        "O".repeat(20),
        &HEAVY[27..]
    ));
    let choice = selection(NumberingScheme::Imgt);
    let report = analyze(&long, &options(), &limits(), &choice, &mut || false).unwrap();
    let track = report
        .tracks
        .iter()
        .find(|t| t.locus == Locus::Heavy)
        .unwrap();
    assert_loop(
        track,
        Region::Cdr1,
        22,
        55,
        &format!("AASGF{}TFSSYAMS", "O".repeat(20)),
    );
    assert_loop(track, Region::Cdr2, 69, 79, "AISGSGGSTY");
    assert_loop(track, Region::Cdr3, 116, 129, "ARDYGMDVAAADA");
    assert!(
        track
            .residues
            .iter()
            .any(|r| r.region == Some(Region::Cdr1) && r.assignment.is_none())
    );
    assert_eq!(track.unassigned_regions, 0);
    report
        .validate_for(&long, &options(), &limits(), &choice)
        .unwrap();
    let inserted = source(&format!("{HEAVY}ARDYGMDVAAADAAAA{HJ}"));
    let report = analyze(&inserted, &options(), &limits(), &choice, &mut || false).unwrap();
    let track = report
        .tracks
        .iter()
        .find(|t| t.locus == Locus::Heavy)
        .unwrap();
    assert_loop(track, Region::Cdr3, 96, 112, "ARDYGMDVAAADAAAA");
    assert!(track.residues.iter().any(|r| {
        r.region == Some(Region::Cdr3)
            && r.assignment
                .as_ref()
                .is_some_and(|a| matches!(a.position.insertion, Insertion::Decimal(_)))
    }));
}

#[test]
fn north_source_selection_membership_intervals_and_current_engine_remain_bound() {
    let molecule = source(&format!("{HEAVY}{H3}{HJ}"));
    let choice = selection(NumberingScheme::Imgt);
    let report = analyze(&molecule, &options(), &limits(), &choice, &mut || false).unwrap();
    for mutation in 0..5 {
        let mut forged = report.clone();
        match mutation {
            0 => forged.selection.region_definition = RegionDefinition::Imgt,
            1 => {
                forged.tracks[0]
                    .residues
                    .iter_mut()
                    .find(|r| r.source_index == 22)
                    .unwrap()
                    .region = Some(Region::Fr1)
            }
            2 => forged.tracks[0].cdrs[0].status = CdrStatus::Mapped { start: 23, end: 35 },
            3 => forged.engine = "antibody-domain-tracks-4".into(),
            4 => forged.rule_revision = "f".repeat(64),
            _ => unreachable!(),
        }
        assert!(
            forged.validate_binding(&report.context_revision).is_err(),
            "mutation{mutation}"
        );
    }
    let mut changed = molecule.clone();
    changed.name.push_str(" changed");
    assert!(
        report
            .validate_for(&changed, &options(), &limits(), &choice)
            .is_err()
    );
    let other = Selection {
        numbering_scheme: NumberingScheme::Imgt,
        region_definition: RegionDefinition::Kabat1991,
    };
    assert!(
        report
            .validate_for(&molecule, &options(), &limits(), &other)
            .is_err()
    );
}
