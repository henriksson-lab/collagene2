use collagene_antibody_numbering::{
    domain_candidates as d,
    domain_tracks::*,
    profile,
    references::{Species, load},
    v_mapping::{self as v, Locus},
};
use collagene_model::{Alphabet, Molecule, Sequence, Topology};
const V: &str = "EVQLLESGGGLVQPGGSLRLSCAASGFTFSSYAMSWVRQAPGKGLEWVSAISGSGGSTYYADSVKGRFTISRDNSKNTLYLQMNSLRAEDTAVYYC";
const J: &str = "WGQGTLVTVSS";
const LOOP: &str = "ARDYGMDVAAADA";
fn source(s: &str) -> Molecule {
    Molecule {
        name: "Authored alternate-track fixture".into(),
        sequence: Sequence::new(Alphabet::Protein, s).unwrap(),
        topology: Topology::Linear,
        features: vec![],
        metadata: Default::default(),
        original_genbank: None,
    }
}
fn options() -> d::Options {
    d::Options {
        references: collagene_antibody_numbering::v_mapping::ReferenceSelection::human_all(),
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
fn choice(s: NumberingScheme, d: RegionDefinition) -> Selection {
    Selection {
        numbering_scheme: s,
        region_definition: d,
    }
}
fn row(t: &Track, index: usize) -> &Residue {
    t.residues.iter().find(|r| r.source_index == index).unwrap()
}
fn sequence(t: &Track, kind: Region) -> Option<String> {
    match t.cdrs.iter().find(|c| c.kind == kind).unwrap().status {
        CdrStatus::Unresolved => None,
        CdrStatus::Mapped { start, end } => Some(
            t.residues
                .iter()
                .filter(|r| r.source_index >= start && r.source_index < end)
                .map(|r| r.residue)
                .collect(),
        ),
    }
}
#[test]
fn heavy_framework_insertions_windows_and_j_tail_retain_exact_source_coordinates() {
    let m = source(&format!("{V}{LOOP}{J}"));
    for scheme in [NumberingScheme::Kabat1991, NumberingScheme::Chothia1997] {
        let c = choice(scheme, RegionDefinition::Imgt);
        let r = analyze(&m, &options(), &limits(), &c, &mut || false).unwrap();
        assert_eq!(
            r.context_revision,
            request_revision(&m, &options(), &limits(), &c).unwrap()
        );
        let t = &r.tracks[0];
        for (index, label) in [
            (21, "22"),
            (35, "36"),
            (51, "52"),
            (52, "52A"),
            (83, "82A"),
            (84, "82B"),
            (85, "82C"),
            (V.len() - 1, "92"),
            (V.len(), "93"),
            (V.len() + 13, "103"),
            (m.sequence.len() - 1, "113"),
        ] {
            assert_eq!(
                row(t, index)
                    .assignment
                    .as_ref()
                    .unwrap()
                    .position
                    .to_string(),
                label
            );
        }
        assert_eq!(t.unassigned_residues, 0);
        assert_eq!(t.unassigned_regions, 0);
        assert!(t.issues.is_empty());
        r.validate_binding(&r.context_revision).unwrap();
        r.validate_for(&m, &options(), &limits(), &c).unwrap();
    }
}
#[test]
fn definition_changes_physical_cdrs_without_changing_numbering_or_scientific_candidate_identity() {
    let m = source(&format!("{V}{LOOP}{J}"));
    for scheme in [
        NumberingScheme::Imgt,
        NumberingScheme::Kabat1991,
        NumberingScheme::Chothia1997,
    ] {
        let a = analyze(
            &m,
            &options(),
            &limits(),
            &choice(scheme.clone(), RegionDefinition::Imgt),
            &mut || false,
        )
        .unwrap();
        let b = analyze(
            &m,
            &options(),
            &limits(),
            &choice(scheme, RegionDefinition::Kabat1991),
            &mut || false,
        )
        .unwrap();
        assert_ne!(a.context_revision, b.context_revision);
        assert_eq!(a.evidence, b.evidence);
        let x = &a.tracks[0];
        let y = &b.tracks[0];
        assert_eq!(
            x.residues.iter().map(|r| &r.assignment).collect::<Vec<_>>(),
            y.residues.iter().map(|r| &r.assignment).collect::<Vec<_>>()
        );
        assert_eq!(sequence(x, Region::Cdr1).as_deref(), Some("GFTFSSYA"));
        assert_eq!(sequence(y, Region::Cdr1).as_deref(), Some("SYAMS"));
        assert_eq!(sequence(x, Region::Cdr3).as_deref(), Some(LOOP));
        assert_eq!(sequence(y, Region::Cdr3).as_deref(), Some(&LOOP[2..]));
        assert_eq!(row(x, 25).region, Some(Region::Cdr1));
        assert_eq!(row(y, 25).region, Some(Region::Fr1));
        assert_eq!(row(x, 33).region, Some(Region::Fr2));
        assert_eq!(row(y, 33).region, Some(Region::Cdr1));
        b.validate_binding(&b.context_revision).unwrap();
    }
}
#[test]
fn unsupported_numbering_window_can_keep_a_known_region_without_inventing_numbers() {
    let m = source(&format!("{V}ARDYGMDV{J}"));
    let r = analyze(
        &m,
        &options(),
        &limits(),
        &choice(NumberingScheme::Kabat1991, RegionDefinition::Imgt),
        &mut || false,
    )
    .unwrap();
    let t = &r.tracks[0];
    assert_eq!(t.unassigned_residues, 8);
    assert_eq!(t.unassigned_regions, 0);
    for index in V.len()..V.len() + 8 {
        assert!(row(t, index).assignment.is_none());
        assert_eq!(row(t, index).region, Some(Region::Cdr3));
    }
    assert_eq!(sequence(t, Region::Cdr3).as_deref(), Some("ARDYGMDV"));
    assert!(!t.issues.is_empty());
    let r = analyze(
        &m,
        &options(),
        &limits(),
        &choice(NumberingScheme::Kabat1991, RegionDefinition::Kabat1991),
        &mut || false,
    )
    .unwrap();
    assert_eq!(r.tracks[0].unassigned_regions, 8);
    assert_eq!(sequence(&r.tracks[0], Region::Cdr3), None);
    r.validate_binding(&r.context_revision).unwrap();
}
#[test]
fn unsupported_first_window_does_not_discard_other_windows_or_make_a_partial_cdr_complete() {
    let m = source(&format!(
        "{}{}{}{LOOP}{J}",
        &V[..27],
        "O".repeat(20),
        &V[27..]
    ));
    let c = choice(NumberingScheme::Kabat1991, RegionDefinition::Kabat1991);
    let r = analyze(&m, &options(), &limits(), &c, &mut || false).unwrap();
    assert_eq!(sequence(&r.tracks[0], Region::Cdr1), None);
    assert_eq!(
        sequence(&r.tracks[0], Region::Cdr3).as_deref(),
        Some(&LOOP[2..])
    );
    assert!(
        r.tracks[0]
            .issues
            .iter()
            .any(|i| i.window == collagene_antibody_numbering::numbering_windows::Window::CysToTrp)
    );
    r.validate_for(&m, &options(), &limits(), &c).unwrap();
}
#[test]
fn light_locus_boundaries_and_complete_definition_intervals_are_source_bound() {
    let corpus = load(&Species::Human).unwrap();
    for (gene, tail, locus) in [
        ("IGKV1-39*01", "FGQGTKVEIK", Locus::Kappa),
        ("IGLV2-14*01", "FGGGTKLTVL", Locus::Lambda),
    ] {
        let prefix = corpus
            .v_entries()
            .iter()
            .find(|v| v.label == gene)
            .unwrap()
            .aligned_protein
            .replace('-', "");
        let m = source(&format!("{prefix}QQSYSTPPT{tail}"));
        let c = choice(NumberingScheme::Kabat1991, RegionDefinition::Kabat1991);
        let r = analyze(&m, &options(), &limits(), &c, &mut || false).unwrap();
        let t = r.tracks.iter().find(|t| t.locus == locus).unwrap();
        assert_eq!(
            row(t, prefix.len() - 1)
                .assignment
                .as_ref()
                .unwrap()
                .position
                .to_string(),
            "88"
        );
        assert_eq!(
            row(t, prefix.len() + 9)
                .assignment
                .as_ref()
                .unwrap()
                .position
                .to_string(),
            "98"
        );
        assert_eq!(
            row(t, m.sequence.len() - 1)
                .assignment
                .as_ref()
                .unwrap()
                .position
                .to_string(),
            "107"
        );
        assert_eq!(sequence(t, Region::Cdr3).as_deref(), Some("QQSYSTPPT"));
        r.validate_binding(&r.context_revision).unwrap();
    }
}
#[test]
fn selected_rules_membership_and_complete_cdr_intervals_reconstruct_in_current_reports() {
    let m = source(&format!("{V}{LOOP}{J}"));
    let c = choice(NumberingScheme::Kabat1991, RegionDefinition::Kabat1991);
    let r = analyze(&m, &options(), &limits(), &c, &mut || false).unwrap();
    let mut bad = r.clone();
    bad.tracks[0].residues[0].region = Some(Region::Cdr1);
    assert!(bad.validate_binding(&r.context_revision).is_err());
    let mut bad = r.clone();
    bad.tracks[0].cdrs[0].status = CdrStatus::Mapped { start: 0, end: 1 };
    assert!(bad.validate_binding(&r.context_revision).is_err());
    let mut bad = r.clone();
    bad.selection.region_definition = RegionDefinition::Imgt;
    assert!(bad.validate_binding(&r.context_revision).is_err());
    let mut bad = r.clone();
    bad.engine = "antibody-domain-tracks-2".into();
    assert!(bad.validate_binding(&r.context_revision).is_err());
    assert!(
        r.validate_for(
            &m,
            &options(),
            &limits(),
            &choice(NumberingScheme::Imgt, RegionDefinition::Kabat1991)
        )
        .is_err()
    );
}

#[test]
fn canonical_chart_boundaries_keep_locus_and_plain_52_distinct_from_its_insertions() {
    for (locus, checks) in [
        (
            Locus::Heavy,
            vec![
                ("25", Region::Fr1),
                ("26", Region::Cdr1),
                ("32", Region::Cdr1),
                ("33", Region::Fr2),
                ("52", Region::Fr2),
                ("52A", Region::Cdr2),
                ("55", Region::Cdr2),
                ("56", Region::Fr3),
                ("95", Region::Fr3),
                ("96", Region::Cdr3),
                ("100C", Region::Cdr3),
                ("101", Region::Cdr3),
                ("102", Region::Fr4),
                ("113", Region::Fr4),
            ],
        ),
        (
            Locus::Kappa,
            vec![
                ("25", Region::Fr1),
                ("26", Region::Cdr1),
                ("30E", Region::Cdr1),
                ("32", Region::Cdr1),
                ("33", Region::Fr2),
                ("49", Region::Fr2),
                ("50", Region::Cdr2),
                ("52", Region::Cdr2),
                ("53", Region::Fr3),
                ("90", Region::Fr3),
                ("91", Region::Cdr3),
                ("96", Region::Cdr3),
                ("97", Region::Fr4),
            ],
        ),
        (
            Locus::Lambda,
            vec![
                ("24", Region::Fr1),
                ("25", Region::Cdr1),
                ("32", Region::Cdr1),
                ("33", Region::Fr2),
                ("50", Region::Cdr2),
                ("52", Region::Cdr2),
                ("53", Region::Fr3),
                ("90", Region::Fr3),
                ("91", Region::Cdr3),
                ("95B", Region::Cdr3),
                ("96", Region::Cdr3),
                ("97", Region::Fr4),
            ],
        ),
    ] {
        for (label, expected) in checks {
            assert_eq!(
                chothia_canonical_region(&locus, &label.parse().unwrap()).unwrap(),
                expected,
                "{locus:?} {label}"
            );
        }
        for bad in ["114", "111.1"] {
            assert!(chothia_canonical_region(&locus, &bad.parse().unwrap()).is_err());
        }
    }
}
#[test]
fn canonical_definition_changes_physical_loops_under_each_independent_numbering_choice() {
    let m = source(&format!("{V}{LOOP}{J}"));
    for scheme in [
        NumberingScheme::Imgt,
        NumberingScheme::Kabat1991,
        NumberingScheme::Chothia1997,
    ] {
        let baseline = analyze(
            &m,
            &options(),
            &limits(),
            &choice(scheme.clone(), RegionDefinition::Imgt),
            &mut || false,
        )
        .unwrap();
        let c = choice(scheme, RegionDefinition::ChothiaCanonical);
        let r = analyze(&m, &options(), &limits(), &c, &mut || false).unwrap();
        let t = &r.tracks[0];
        assert_eq!(baseline.evidence, r.evidence);
        assert_eq!(
            baseline.tracks[0]
                .residues
                .iter()
                .map(|r| &r.assignment)
                .collect::<Vec<_>>(),
            t.residues.iter().map(|r| &r.assignment).collect::<Vec<_>>()
        );
        assert_eq!(sequence(t, Region::Cdr1).as_deref(), Some("GFTFSSY"));
        assert_eq!(sequence(t, Region::Cdr2).as_deref(), Some("GSGG"));
        assert_eq!(sequence(t, Region::Cdr3).as_deref(), Some("YGMDVAAAD"));
        assert_eq!(row(t, 51).region, Some(Region::Fr2));
        assert_eq!(row(t, 52).region, Some(Region::Cdr2));
        assert_eq!(row(t, V.len() + 2).region, Some(Region::Fr3));
        assert_eq!(row(t, V.len() + 3).region, Some(Region::Cdr3));
        assert_eq!(row(t, V.len() + 12).region, Some(Region::Fr4));
        assert!(t.issues.is_empty());
        assert_eq!(t.unassigned_regions, 0);
        r.validate_for(&m, &options(), &limits(), &c).unwrap();
        r.validate_binding(&r.context_revision).unwrap();
    }
    // A supported IMGT junction can remain unresolved under the finite Chothia
    // window profile without losing its established first/middle loops.
    let short = source(&format!("{V}ARDYGMDV{J}"));
    let c = choice(NumberingScheme::Imgt, RegionDefinition::ChothiaCanonical);
    let r = analyze(&short, &options(), &limits(), &c, &mut || false).unwrap();
    assert_eq!(sequence(&r.tracks[0], Region::Cdr3), None);
    assert_eq!(
        sequence(&r.tracks[0], Region::Cdr1).as_deref(),
        Some("GFTFSSY")
    );
    assert!(row(&r.tracks[0], V.len() + 3).assignment.is_some());
    assert!(row(&r.tracks[0], V.len() + 3).region.is_none());
    assert!(
        r.tracks[0]
            .issues
            .iter()
            .any(|i| i.window == collagene_antibody_numbering::numbering_windows::Window::CysToJ)
    );
}
#[test]
fn canonical_light_cdrs_are_complete_source_intervals_with_distinct_locus_start_positions() {
    let corpus = load(&Species::Human).unwrap();
    for (gene, tail, locus, start) in [
        ("IGKV1-39*01", "FGQGTKVEIK", Locus::Kappa, 26),
        ("IGLV2-14*01", "FGGGTKLTVL", Locus::Lambda, 25),
    ] {
        let prefix = corpus
            .v_entries()
            .iter()
            .find(|v| v.label == gene)
            .unwrap()
            .aligned_protein
            .replace('-', "");
        let m = source(&format!("{prefix}QQSYSTPPT{tail}"));
        let c = choice(
            NumberingScheme::Chothia1997,
            RegionDefinition::ChothiaCanonical,
        );
        let r = analyze(&m, &options(), &limits(), &c, &mut || false).unwrap();
        let t = r.tracks.iter().find(|t| t.locus == locus).unwrap();
        assert_eq!(sequence(t, Region::Cdr3).as_deref(), Some("SYSTPP"));
        let CdrStatus::Mapped { start: begin, end } = t.cdrs[0].status else {
            panic!("complete first loop")
        };
        assert_eq!(
            row(t, begin)
                .assignment
                .as_ref()
                .unwrap()
                .position
                .to_string(),
            start.to_string()
        );
        assert_eq!(
            row(t, end - 1)
                .assignment
                .as_ref()
                .unwrap()
                .position
                .to_string(),
            "32"
        );
        let mut forged = r.clone();
        forged.tracks[0].region_definition = RegionDefinition::Kabat1991;
        assert!(forged.validate_binding(&r.context_revision).is_err());
        r.validate_for(&m, &options(), &limits(), &c).unwrap();
    }
}
