use collagene_core::protein_properties::*;
use collagene_model::{Alphabet, Feature, Location, Molecule, Sequence, Topology};
fn protein(sequence: &str) -> Molecule {
    Molecule {
        name: "independent protein fixture".into(),
        sequence: Sequence::new(Alphabet::Protein, sequence).unwrap(),
        topology: Topology::Linear,
        features: vec![],
        metadata: Default::default(),
        original_genbank: None,
    }
}
fn report(sequence: &str) -> ProteinReport {
    analyze_protein(&protein(sequence), &ProteinOptions::default()).unwrap()
}
fn close(a: f64, b: f64, tolerance: f64) {
    assert!((a - b).abs() <= tolerance, "{a} != {b}");
}
#[test]
fn every_canonical_free_amino_acid_matches_published_residue_mass_plus_water() {
    // Independently transcribed UWPR table; formula code never consults this fixture.
    let rows = [
        ('G', 57.021463735, 57.05132),
        ('A', 71.037113805, 71.0779),
        ('S', 87.032028435, 87.0773),
        ('P', 97.052763875, 97.11518),
        ('V', 99.068413945, 99.13106),
        ('T', 101.047678505, 101.10388),
        ('C', 103.009184505, 103.1429),
        ('L', 113.084064015, 113.15764),
        ('I', 113.084064015, 113.15764),
        ('N', 114.042927470, 114.10264),
        ('D', 115.026943065, 115.0874),
        ('Q', 128.058577540, 128.12922),
        ('K', 128.094963050, 128.17228),
        ('E', 129.042593135, 129.11398),
        ('M', 131.040484645, 131.19606),
        ('H', 137.058911875, 137.13928),
        ('F', 147.068413945, 147.17386),
        ('R', 156.101111050, 156.18568),
        ('Y', 163.063328575, 163.17326),
        ('W', 186.079312980, 186.2099),
    ];
    for (symbol, mono, average) in rows {
        let r = report(&symbol.to_string());
        r.validate().unwrap();
        assert!(r.complete);
        let mass = r.mass.unwrap();
        close(mass.monoisotopic_da.minimum, mono + 18.010564700, 1e-9);
        close(mass.average_da.minimum, average + 18.01528, 1e-8);
        assert_eq!(mass.average_da.minimum, mass.average_da.maximum);
    }
}
#[test]
fn peptide_formation_removes_one_water_per_new_peptide_bond() {
    let r = report("AGW").mass.unwrap();
    // Three free residue masses minus two neutral waters. Not an amino-mass lookup sum.
    let free = [89.09318, 75.06660, 204.22518];
    close(
        r.average_da.minimum,
        free.iter().sum::<f64>() - 2. * 18.01528,
        1e-8,
    );
    close(
        r.monoisotopic_da.minimum,
        71.037113805 + 57.021463735 + 186.079312980 + 18.010564700,
        1e-9,
    );
}
#[test]
fn independent_high_precision_newton_roots_agree_with_bisection() {
    // Decimal50-digit Newton solve with an analytic derivative, computed independently
    // from the documented group equations (the implementation uses64-step bisection).
    for (sequence, expected) in [
        ("A", 6.1),
        ("D", 3.7499873618442117),
        ("K", 9.700002530129093),
        ("ACDEFGHIKLMNPQRSTVWY", 7.356542538059439),
        ("WWYYCC", 5.822115235298151),
        ("DDDKKK", 6.519314446400460),
        ("RRRRRRRRRR", 13.454249259279828),
        ("EEEEEEEEEE", 3.031728356330726),
    ] {
        close(report(sequence).isoelectric_point.unwrap(), expected, 2e-11);
    }
    // Increasing a neutral chain's length never creates additional terminal groups.
    close(report("AAAA").isoelectric_point.unwrap(), 6.1, 1e-11);
}
#[test]
fn zero_charge_unbracketed_in_supported_ph_range_never_becomes_fake_endpoint_pi() {
    let r = report(&"R".repeat(10000));
    assert!(r.isoelectric_point.is_none());
    assert!(r.net_charge_at_ph.unwrap() > 9999.);
    assert!(!r.complete);
    assert!(r.unavailable.iter().any(|s| s.contains("bracketed")));
}
#[test]
fn disulfides_change_mass_charge_and_extinction_consistently() {
    let m = protein("WWYYCCC");
    let reduced = analyze_protein(&m, &ProteinOptions::default()).unwrap();
    let options = ProteinOptions {
        disulfides: DisulfideState::MaximumIntrachain,
        ..Default::default()
    };
    let oxidized = analyze_protein(&m, &options).unwrap();
    assert_eq!(
        oxidized.disulfide_bonds,
        NumericInterval {
            minimum: 1.,
            maximum: 1.
        }
    );
    close(
        reduced.mass.as_ref().unwrap().average_da.minimum
            - oxidized.mass.as_ref().unwrap().average_da.minimum,
        2.01588,
        1e-10,
    );
    close(
        reduced.mass.as_ref().unwrap().monoisotopic_da.minimum
            - oxidized.mass.as_ref().unwrap().monoisotopic_da.minimum,
        2.015650070,
        1e-10,
    );
    assert!(oxidized.net_charge_at_ph > reduced.net_charge_at_ph);
    assert!(oxidized.isoelectric_point > reduced.isoelectric_point);
    let epsilon = oxidized.extinction_280.as_ref().unwrap();
    assert_eq!(epsilon.epsilon_molar_inverse_cm.minimum, 14105.);
    assert_eq!(
        reduced
            .extinction_280
            .unwrap()
            .epsilon_molar_inverse_cm
            .minimum,
        13980.
    );
    close(
        epsilon.absorbance_1mg_ml_1cm.as_ref().unwrap().minimum,
        14105. / oxidized.mass.unwrap().average_da.minimum,
        1e-12,
    );
    assert!(oxidized.warnings.iter().any(|s| s.contains("unmatched")));
    assert!(
        analyze_protein(
            &m,
            &ProteinOptions {
                disulfides: DisulfideState::Count(2),
                ..Default::default()
            }
        )
        .is_err()
    );
    assert_eq!(
        analyze_protein(
            &m,
            &ProteinOptions {
                disulfides: DisulfideState::Count(1),
                ..Default::default()
            }
        )
        .unwrap()
        .disulfide_bonds
        .maximum,
        1.
    );
}
#[test]
fn terminal_stop_policy_is_explicit_and_internal_stops_never_join_chains() {
    let a = report("MK");
    let b = report("MK*");
    assert_eq!(a.mass, b.mass);
    assert!(b.terminal_stop_removed);
    assert_ne!(a.source_revision, b.source_revision);
    assert_ne!(a.source_sequence_blake3, b.source_sequence_blake3);
    assert_eq!(b.residue_count, 2);
    assert_eq!(b.source_length, 3);
    for seq in ["*", "M*K", "MK**", "*MK"] {
        assert!(analyze_protein(&protein(seq), &ProteinOptions::default()).is_err());
    }
    assert!(
        analyze_protein(
            &protein("MK*"),
            &ProteinOptions {
                stop_policy: StopPolicy::Reject,
                ..Default::default()
            }
        )
        .is_err()
    );
}
#[test]
fn ambiguous_composition_is_retained_without_invented_frequencies() {
    let r = report("ABZJXUO");
    assert_eq!(r.composition.len(), 7);
    assert!(r.composition.iter().all(|c| c.count == 1));
    assert!(r.mass.is_none() && r.isoelectric_point.is_none() && r.extinction_280.is_none());
    assert!(!r.complete);
    r.validate().unwrap();
    let j = report("J");
    let il = report("I");
    assert_eq!(
        j.mass.as_ref().unwrap().average_da,
        il.mass.unwrap().average_da
    );
    assert!(!j.mass.unwrap().exact_identity);
    assert!(!j.complete);
    assert!(j.isoelectric_point.is_some());
    let b = report("B");
    assert!(b.mass.is_none() && b.isoelectric_point.is_none());
    assert_eq!(
        b.extinction_280.unwrap().epsilon_molar_inverse_cm.minimum,
        0.
    );
}
#[test]
fn canonical_bounds_enclose_all_literal_resolutions_including_disulfide_uncertainty() {
    let o = ProteinOptions {
        ambiguity: AmbiguityPolicy::CanonicalBounds,
        disulfides: DisulfideState::MaximumIntrachain,
        ..Default::default()
    };
    let source = analyze_protein(&protein("BCZX"), &o).unwrap();
    source.validate().unwrap();
    assert!(!source.complete);
    assert!(source.isoelectric_point.is_none());
    let mass = source.mass.unwrap();
    let e = source.extinction_280.unwrap();
    for b in ['D', 'N'] {
        for z in ['E', 'Q'] {
            for x in "ACDEFGHIKLMNPQRSTVWY".chars() {
                let s = format!("{b}C{z}{x}");
                let r = analyze_protein(&protein(&s), &o).unwrap();
                let m = r.mass.unwrap();
                assert!(
                    m.average_da.minimum >= mass.average_da.minimum - 1e-9
                        && m.average_da.maximum <= mass.average_da.maximum + 1e-9
                );
                assert!(
                    m.monoisotopic_da.minimum >= mass.monoisotopic_da.minimum - 1e-9
                        && m.monoisotopic_da.maximum <= mass.monoisotopic_da.maximum + 1e-9
                );
                let resolved = r.extinction_280.unwrap();
                assert!(
                    resolved.epsilon_molar_inverse_cm.minimum >= e.epsilon_molar_inverse_cm.minimum
                        && resolved.epsilon_molar_inverse_cm.maximum
                            <= e.epsilon_molar_inverse_cm.maximum
                );
                let absorbance = resolved.absorbance_1mg_ml_1cm.unwrap();
                let bounds = e.absorbance_1mg_ml_1cm.as_ref().unwrap();
                assert!(
                    absorbance.minimum >= bounds.minimum - 1e-12
                        && absorbance.maximum <= bounds.maximum + 1e-12
                );
            }
        }
    }
    assert!(
        analyze_protein(
            &protein("XX"),
            &ProteinOptions {
                ambiguity: AmbiguityPolicy::CanonicalBounds,
                disulfides: DisulfideState::Count(1),
                ..Default::default()
            }
        )
        .is_err()
    );
}
#[test]
fn known_modifications_are_not_silently_treated_as_unmodified() {
    let mut m = protein("AKW");
    m.features.push(Feature {
        kind: "MOD_RES".into(),
        location: Location::Range {
            start: 1,
            end: 2,
            before: false,
            after: false,
        },
        qualifiers: vec![("note".into(), Some("acetyllysine".into()))],
        needs_review: false,
    });
    let r = analyze_protein(&m, &Default::default()).unwrap();
    assert!(r.mass.is_none() && r.isoelectric_point.is_none() && r.extinction_280.is_none());
    assert_eq!(r.composition.len(), 3);
    assert!(!r.complete);
    r.validate().unwrap();
    let r = analyze_protein(
        &protein("AKW"),
        &ProteinOptions {
            modification_state: ModificationState::KnownOrUnknownModifications,
            ..Default::default()
        },
    )
    .unwrap();
    assert!(r.mass.is_none());
}
#[test]
fn quantity_bounds_and_source_alphabet_are_validated() {
    assert!(
        analyze_protein(
            &Molecule::dna("dna", "ATG", Topology::Linear).unwrap(),
            &Default::default()
        )
        .is_err()
    );
    for p in [f64::NAN, f64::INFINITY, -0.1, 14.1] {
        assert!(
            analyze_protein(
                &protein("A"),
                &ProteinOptions {
                    charge_ph: p,
                    ..Default::default()
                }
            )
            .is_err()
        );
    }
    assert!(
        analyze_protein(
            &protein("AA"),
            &ProteinOptions {
                max_residues: 1,
                ..Default::default()
            }
        )
        .is_err()
    );
    for maximum in [0, 10_000_001] {
        assert!(
            analyze_protein(
                &protein("A"),
                &ProteinOptions {
                    max_residues: maximum,
                    ..Default::default()
                }
            )
            .is_err()
        );
    }
    let mut circular = protein("A");
    circular.topology = Topology::Circular;
    assert!(analyze_protein(&circular, &Default::default()).is_err());
}
#[test]
fn source_and_model_bindings_and_strict_report_shapes_survive_json() {
    let a = report("WWYYCC");
    let bytes = serde_json::to_vec(&a).unwrap();
    let b: ProteinReport = serde_json::from_slice(&bytes).unwrap();
    assert_eq!(a, b);
    b.validate().unwrap();
    let mut changed = protein("WWYYCC");
    changed.name = "renamed immutable source".into();
    let c = analyze_protein(&changed, &Default::default()).unwrap();
    assert_ne!(a.source_revision, c.source_revision);
    assert_eq!(a.source_sequence_blake3, c.source_sequence_blake3);
    let mut invalid = a.clone();
    invalid.composition[0].count += 1;
    assert!(invalid.validate().is_err());
    let mut invalid = a.clone();
    invalid.mass.as_mut().unwrap().average_da.minimum = f64::NAN;
    assert!(invalid.validate().is_err());
    let mut invalid = a.clone();
    invalid.engine_version = "unrecognized".into();
    assert!(invalid.validate().is_err());
    let mut invalid = a.clone();
    invalid.disulfide_bonds.minimum = 0.5;
    assert!(invalid.validate().is_err());
    let mut value = serde_json::to_value(&a).unwrap();
    value["unrecognized"] = true.into();
    assert!(serde_json::from_value::<ProteinReport>(value).is_err());
}
#[test]
fn source_evidence_budget_and_deep_annotations_fail_before_source_serialization() {
    let mut m = protein("A");
    m.metadata
        .insert("large".into(), "x".repeat(32 * 1024 * 1024));
    assert!(analyze_protein(&m, &Default::default()).is_err());
    let mut m = protein("A");
    let mut location = Location::Range {
        start: 0,
        end: 1,
        before: false,
        after: false,
    };
    for _ in 0..70 {
        location = Location::Complement(Box::new(location));
    }
    m.features.push(Feature {
        kind: "region".into(),
        location,
        qualifiers: vec![],
        needs_review: false,
    });
    assert!(analyze_protein(&m, &Default::default()).is_err());
}
