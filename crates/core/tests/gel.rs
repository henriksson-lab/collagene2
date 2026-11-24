//! Independent lengths, log-coordinate geometry, mass accounting and diagnostic cases.
use collagene_core::{analysis::Primer, gel::*, *};
fn dna(s: &str, t: Topology) -> Molecule {
    Molecule::dna("gel fixture", s, t).unwrap()
}
fn wide() -> Calibration {
    Calibration {
        label: "Independent illustrative span".into(),
        points: vec![
            CalibrationPoint {
                size_bp: 1.,
                migration: 0.9,
            },
            CalibrationPoint {
                size_bp: 1000.,
                migration: 0.1,
            },
        ],
        resolution: 0.01,
    }
}
fn input(molecule: Molecule, preparation: Preparation) -> LaneInput {
    LaneInput {
        label: "Fixture".into(),
        molecule,
        preparation,
        mass_ng: None,
    }
}
fn digest(molecule: Molecule, names: &[&str]) -> LaneInput {
    input(
        molecule,
        Preparation::Digest {
            enzymes: names.iter().map(|s| s.to_string()).collect(),
        },
    )
}
fn near(a: f64, b: f64) {
    assert!((a - b).abs() < 1e-10, "{a} != {b}");
}
#[test]
fn calibration_is_piecewise_log_size_and_never_extrapolates() {
    let c = Calibration::default();
    near(c.position(1000.).unwrap().unwrap(), 0.5);
    near(
        c.position(100f64.sqrt() * 1000f64.sqrt()).unwrap().unwrap(),
        0.7,
    );
    assert!(c.position(99.).unwrap().is_none());
    assert!(c.position(10001.).unwrap().is_none());
    assert!(c.position(f64::NAN).is_err());
    let custom = Calibration {
        label: "Hand coordinate pair".into(),
        points: vec![
            CalibrationPoint {
                size_bp: 100.,
                migration: 0.8,
            },
            CalibrationPoint {
                size_bp: 400.,
                migration: 0.6,
            },
            CalibrationPoint {
                size_bp: 1600.,
                migration: 0.1,
            },
        ],
        resolution: 0.01,
    };
    near(custom.position(200.).unwrap().unwrap(), 0.7);
    near(custom.position(800.).unwrap().unwrap(), 0.35);
    let a = 10_000_000. - 0.000000004;
    let b = 10_000_000.;
    let close = Calibration {
        label: "Nearby floating point coordinates".into(),
        points: vec![
            CalibrationPoint {
                size_bp: a,
                migration: 0.9,
            },
            CalibrationPoint {
                size_bp: b,
                migration: 0.1,
            },
        ],
        resolution: 0.01,
    };
    near(close.position(b).unwrap().unwrap(), 0.1);
    near(close.position((a + b) / 2.).unwrap().unwrap(), 0.5);
}
#[test]
fn linear_blunt_digest_preserves_all_fragments_and_mass_under_comigration() {
    // 4+6+4+6+4 =24; EcoRV top/bottom boundaries 7 and17 ->7,10,7.
    let mut lane = digest(
        dna("AAAAGATATCCCCCGATATCGGGG", Topology::Linear),
        &["EcoRV"],
    );
    lane.mass_ng = Some(240.);
    let report = predict_gel(&[lane], &wide()).unwrap();
    let r = &report.lanes[0];
    assert_eq!(
        r.fragments.iter().map(|f| f.top_nt).collect::<Vec<_>>(),
        vec![7, 10, 7]
    );
    assert_eq!(
        r.fragments.iter().map(|f| f.bottom_nt).collect::<Vec<_>>(),
        vec![7, 10, 7]
    );
    assert_eq!(r.bands.len(), 2);
    assert!(r.complete_under_model);
    near(r.fragments.iter().map(|f| f.amount_ng.unwrap()).sum(), 240.);
    let repeated = r
        .bands
        .iter()
        .find(|b| b.fragment_indexes.len() == 2)
        .unwrap();
    near(repeated.min_size_bp, 7.);
    near(repeated.max_size_bp, 7.);
    near(repeated.amount_ng.unwrap(), 140.);
}
#[test]
fn circle_is_not_given_a_linear_mobility_until_cleaved() {
    let m = dna("AAAAGATATCCCCCGATATCGGGG", Topology::Circular);
    let report = predict_gel(
        &[input(m.clone(), Preparation::Uncut), digest(m, &["EcoRV"])],
        &wide(),
    )
    .unwrap();
    assert_eq!(
        report.lanes[0].fragments[0].state,
        MobilityState::CircularConformationUnknown
    );
    assert!(report.lanes[0].bands.is_empty());
    assert_eq!(report.lanes[0].unpositioned_fragments, vec![0]);
    let mut sizes = report.lanes[1]
        .fragments
        .iter()
        .map(|f| f.top_nt)
        .collect::<Vec<_>>();
    sizes.sort();
    assert_eq!(sizes, vec![10, 14]);
    assert_eq!(report.lanes[1].bands.len(), 2);
}
#[test]
fn sticky_end_sizes_use_both_strands_and_are_explicitly_approximate() {
    // PstI CTGCA|G / G|ACGTC at prefix4 -> opposite 9/5 and5/9 fragments.
    let mut lane = digest(dna("AAAACTGCAGTTTT", Topology::Linear), &["PstI"]);
    lane.mass_ng = Some(40.);
    let report = predict_gel(&[lane], &wide()).unwrap();
    let r = &report.lanes[0];
    assert_eq!(
        r.fragments
            .iter()
            .map(|f| (f.top_nt, f.bottom_nt))
            .collect::<Vec<_>>(),
        vec![(9, 5), (5, 9)]
    );
    assert!(r.fragments.iter().all(|f| f.size_bp == 7.));
    assert_eq!(r.bands.len(), 1);
    near(r.bands[0].amount_ng.unwrap(), 40.);
    assert!(r.warnings.iter().any(|w| w.contains("mean strand length")));
}
#[test]
fn resolution_grouping_does_not_bridge_a_chain_of_adjacent_sizes() {
    // Precisely positioned blunt cuts100,210 in330bp ->100,110,120.
    let sequence = format!(
        "{}GATATC{}GATATC{}",
        "A".repeat(97),
        "A".repeat(104),
        "A".repeat(117)
    );
    let c = Calibration {
        label: "Artificial three-point resolution fixture".into(),
        points: vec![
            CalibrationPoint {
                size_bp: 100.,
                migration: 0.6,
            },
            CalibrationPoint {
                size_bp: 110.,
                migration: 0.58,
            },
            CalibrationPoint {
                size_bp: 120.,
                migration: 0.56,
            },
        ],
        resolution: 0.025,
    };
    let r = predict_gel(&[digest(dna(&sequence, Topology::Linear), &["EcoRV"])], &c).unwrap();
    assert_eq!(
        r.lanes[0]
            .fragments
            .iter()
            .map(|f| f.top_nt)
            .collect::<Vec<_>>(),
        vec![100, 110, 120]
    );
    assert_eq!(r.lanes[0].bands.len(), 2);
    assert!(
        r.lanes[0]
            .bands
            .iter()
            .all(|b| b.migration_end - b.migration_start <= c.resolution)
    );
}
#[test]
fn tailed_pcr_band_is_from_reconstructed_product_and_unknown_alternatives_have_no_mass_guess() {
    let mut lane = input(
        dna("ACGACCCCGGGGTAAC", Topology::Linear),
        Preparation::Pcr {
            forward: Primer::new("F", "ACGA", "CGT").unwrap(),
            reverse: Primer::new("R", "GTTA", "GAC").unwrap(),
        },
    );
    lane.mass_ng = Some(50.);
    let r = predict_gel(&[lane], &wide()).unwrap();
    let lane = &r.lanes[0];
    assert_eq!(lane.fragments.len(), 1);
    near(lane.fragments[0].size_bp, 22.);
    near(lane.bands[0].amount_ng.unwrap(), 50.);
    assert!(lane.complete_under_model);
    let mut ambiguous = input(
        dna("AAAACCCCTTTTAAAACCCCTTTT", Topology::Linear),
        Preparation::Pcr {
            forward: Primer::new("F", "AAAA", "").unwrap(),
            reverse: Primer::new("R", "AAAA", "").unwrap(),
        },
    );
    ambiguous.mass_ng = Some(10.);
    let r = predict_gel(&[ambiguous], &wide()).unwrap();
    assert!(r.lanes[0].fragments.len() > 1);
    assert!(!r.lanes[0].complete_under_model);
    assert!(r.lanes[0].fragments.iter().all(|f| f.amount_ng.is_none()));
}
#[test]
fn unsupported_sources_and_out_of_ladder_fragments_are_not_silent_negative_bands() {
    let mut physical = dna("ACGTACGT", Topology::Linear);
    physical
        .metadata
        .insert("duplex_model".into(), "untrusted".into());
    let lanes = vec![
        digest(dna("AAAANGATATCCCC", Topology::Linear), &["EcoRV"]),
        input(physical, Preparation::Uncut),
        input(dna("ACGT", Topology::Linear), Preparation::Uncut),
    ];
    let r = predict_gel(&lanes, &Calibration::default()).unwrap();
    assert!(r.lanes[0].error.is_some());
    assert!(r.lanes[1].error.is_some());
    assert!(!r.lanes[0].complete_under_model);
    assert!(r.lanes[2].error.is_none());
    assert_eq!(r.lanes[2].unpositioned_fragments, vec![0]);
    assert!(r.lanes[2].bands.is_empty());
}
fn options() -> DiagnosticOptions {
    DiagnosticOptions {
        enzymes: vec!["BamHI".into(), "EcoRV".into()],
        include_double_digests: true,
        calibration: wide(),
        ..Default::default()
    }
}
#[test]
fn diagnostic_digest_ranks_real_discrimination_above_unknown_uncut_circles() {
    let a = dna("AAAAGATATCCCCCGATATCGGGG", Topology::Circular);
    let b = dna("AAAAGATATCCCCCGTATACGGGG", Topology::Circular);
    assert_eq!(a.sequence.len(), 24);
    assert_eq!(b.sequence.len(), 24);
    let report = diagnostic_digests(&[a.clone(), b.clone()], &options()).unwrap();
    assert_eq!(report.source_revisions, vec![a.revision(), b.revision()]);
    assert!(report.exhaustive_search);
    assert_eq!(report.evaluated_digests, 3);
    assert_eq!(report.candidates[0].enzymes, vec!["EcoRV"]);
    assert_eq!(report.candidates[0].distinguished_pairs, vec![(0, 1)]);
    let uncut = report
        .candidates
        .iter()
        .find(|c| c.enzymes == ["BamHI"])
        .unwrap();
    assert_eq!(uncut.unknown_pairs, vec![(0, 1)]);
    assert!(uncut.distinguished_pairs.is_empty());
}
#[test]
fn identical_resolved_bands_do_not_become_diagnostic_from_copy_count_or_mass() {
    let a = dna("GATATCAAAAAAGATATCAAAAAA", Topology::Circular);
    let b = dna("GATATCAAAAAA", Topology::Circular);
    let opts = DiagnosticOptions {
        enzymes: vec!["EcoRV".into()],
        include_double_digests: false,
        ..options()
    };
    let report = diagnostic_digests(&[a, b], &opts).unwrap();
    assert_eq!(report.candidates[0].profiles[0].fragments.len(), 2);
    assert_eq!(report.candidates[0].profiles[1].fragments.len(), 1);
    assert_eq!(report.candidates[0].indistinguishable_pairs, vec![(0, 1)]);
}
#[test]
fn work_and_enumeration_budgets_never_claim_exhaustive_best_digest() {
    let molecules = vec![
        dna("AAAAGATATCCCCCGATATCGGGG", Topology::Circular),
        dna("AAAAGATATCCCCCGTATACGGGG", Topology::Circular),
    ];
    let mut opts = options();
    opts.max_evaluations = 1;
    let report = diagnostic_digests(&molecules, &opts).unwrap();
    assert!(!report.exhaustive_search);
    assert_eq!(report.evaluated_digests, 1);
    assert_eq!(report.candidates[0].enzymes, vec!["BamHI"]);
    opts.max_total_bases = 1;
    let report = diagnostic_digests(&molecules, &opts).unwrap();
    assert!(!report.exhaustive_search);
    assert!(report.candidates.is_empty());
    assert_eq!(report.work_bases, 0);
}
#[test]
fn calibration_and_lane_bounds_reject_invalid_values_before_rendering() {
    let mut c = wide();
    c.points[1].migration = c.points[0].migration;
    assert!(c.validate().is_err());
    c = wide();
    c.resolution = f64::NAN;
    assert!(c.validate().is_err());
    let mut lane = input(dna("ACGT", Topology::Linear), Preparation::Uncut);
    lane.mass_ng = Some(-1.);
    assert!(predict_gel(&[lane], &wide()).is_err());
    assert!(predict_gel(&[], &wide()).is_err());
    let r = predict_gel(
        &[input(dna("ACGT", Topology::Linear), Preparation::Uncut)],
        &wide(),
    )
    .unwrap();
    let mut malformed = r.clone();
    malformed.lanes[0].bands[0].migration = f64::NAN;
    assert!(render_svg(&malformed).is_err());
    let mut value = serde_json::to_value(wide()).unwrap();
    value["gel_percent"] = serde_json::json!(1.5);
    assert!(serde_json::from_value::<Calibration>(value).is_err());
    let mut prep = serde_json::to_value(Preparation::Digest {
        enzymes: vec!["EcoRI".into()],
    })
    .unwrap();
    prep["Digest"]["assume_complete_chemical_validation"] = serde_json::json!(true);
    assert!(serde_json::from_value::<Preparation>(prep).is_err());
}
#[test]
fn svg_escapes_source_names_retains_unknown_callouts_and_is_deterministic() {
    let mut l = input(dna("ACGTACGT", Topology::Circular), Preparation::Uncut);
    l.label = "<script>&\" hostile".into();
    let report = predict_gel(&[l], &wide()).unwrap();
    let svg = render_svg(&report).unwrap();
    assert!(svg.contains("&lt;script&gt;&amp;&quot;"));
    assert!(!svg.contains("<script>"));
    assert!(svg.contains("1 unplaced"));
    assert!(svg.contains("role=\"img\""));
    assert_eq!(svg, render_svg(&report).unwrap());
    assert!(svg.contains("Ladder"));
    assert!(!svg.contains("NaN"));
    assert!(svg.contains("collagene-gel-evidence"));
    assert!(svg.contains(&report.lanes[0].source_revision));
}

#[test]
fn serialized_source_evidence_is_bounded_before_materializing_revision_buffers() {
    let mut m = dna("ACGT", Topology::Linear);
    m.metadata.insert(
        "large attachment".into(),
        "X".repeat(collagene_core::formats::MAX_IMPORT_BYTES),
    );
    let lane = input(m, Preparation::Uncut);
    assert!(matches!(
        predict_gel(&[lane], &wide()),
        Err(Error::TooLarge)
    ));
}
