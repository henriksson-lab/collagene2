//! Independent scientific review of declared free-chain/oxidation/ambiguity semantics.
use collagene_core::protein::properties::*;
use collagene_model::{Alphabet, Feature, Location, Molecule, Sequence, Topology};
fn protein(sequence: &str) -> Molecule {
    Molecule {
        name: "Independent review".into(),
        sequence: Sequence::new(Alphabet::Protein, sequence).unwrap(),
        topology: Topology::Linear,
        features: vec![],
        metadata: Default::default(),
        original_genbank: None,
    }
}
fn close(a: f64, b: f64) {
    assert!((a - b).abs() < 1e-10, "{a} != {b}");
}
#[test]
fn each_cysteine_is_half_ionized_at_its_declared_pka_and_cystine_removes_two_thiols() {
    let molecule = protein("CC");
    let options = ProteinOptions {
        charge_ph: 8.5,
        ..Default::default()
    };
    let reduced = analyze_protein(&molecule, &options).unwrap();
    let oxidized = analyze_protein(
        &molecule,
        &ProteinOptions {
            disulfides: DisulfideState::Count(1),
            ..options
        },
    )
    .unwrap();
    // pH=pKa makes each reduced Cys contribute exactly -1/2; oxidizing both
    // removes those two groups. Free termini are unchanged and cancel in the difference.
    close(
        oxidized.net_charge_at_ph.unwrap() - reduced.net_charge_at_ph.unwrap(),
        1.,
    );
    assert_eq!(
        oxidized
            .extinction_280
            .unwrap()
            .epsilon_molar_inverse_cm
            .minimum,
        125.
    );
    close(
        reduced.mass.unwrap().average_da.minimum - oxidized.mass.unwrap().average_da.minimum,
        2.01588,
    );
}
#[test]
fn unknown_pair_mass_and_absorbance_bounds_enclose_every_canonical_covalent_resolution() {
    let options = ProteinOptions {
        ambiguity: AmbiguityPolicy::CanonicalBounds,
        disulfides: DisulfideState::MaximumIntrachain,
        ..Default::default()
    };
    let source = analyze_protein(&protein("XX"), &options).unwrap();
    assert!(!source.complete);
    assert!(source.net_charge_at_ph.is_none());
    let bound = source.mass.unwrap();
    let extinction = source.extinction_280.unwrap();
    assert_eq!(
        source.disulfide_bonds,
        NumericInterval {
            minimum: 0.,
            maximum: 1.
        }
    );
    // These are separately enumerated concrete molecules, including CC (one
    // modeled disulfide) and WW (largest aromatic contribution), never averaged X.
    for a in "ACDEFGHIKLMNPQRSTVWY".chars() {
        for b in "ACDEFGHIKLMNPQRSTVWY".chars() {
            let r = analyze_protein(&protein(&format!("{a}{b}")), &options).unwrap();
            let mass = r.mass.unwrap();
            let absorbance = r.extinction_280.unwrap();
            assert!(
                mass.average_da.minimum >= bound.average_da.minimum - 1e-10
                    && mass.average_da.maximum <= bound.average_da.maximum + 1e-10
            );
            assert!(
                mass.monoisotopic_da.minimum >= bound.monoisotopic_da.minimum - 1e-10
                    && mass.monoisotopic_da.maximum <= bound.monoisotopic_da.maximum + 1e-10
            );
            assert!(
                absorbance.epsilon_molar_inverse_cm.minimum
                    >= extinction.epsilon_molar_inverse_cm.minimum
                    && absorbance.epsilon_molar_inverse_cm.maximum
                        <= extinction.epsilon_molar_inverse_cm.maximum
            );
            let a = absorbance.absorbance_1mg_ml_1cm.unwrap();
            let enclosure = extinction.absorbance_1mg_ml_1cm.as_ref().unwrap();
            assert!(
                a.minimum >= enclosure.minimum - 1e-10 && a.maximum <= enclosure.maximum + 1e-10
            );
        }
    }
}
#[test]
fn retained_chemical_evidence_always_blocks_free_unmodified_chain_claims() {
    for kind in [
        "MOD_RES",
        "modified_residue",
        "cross-link",
        "cross_link",
        "bond",
        "LIPID",
        "CARBOHYD",
        "DISULFID",
    ] {
        let mut molecule = protein("CWK");
        molecule.features.push(Feature {
            kind: kind.into(),
            location: Location::range(0, 1),
            qualifiers: vec![],
            needs_review: true,
        });
        let report = analyze_protein(&molecule, &Default::default()).unwrap();
        assert!(
            report.mass.is_none()
                && report.net_charge_at_ph.is_none()
                && report.isoelectric_point.is_none()
                && report.extinction_280.is_none(),
            "{kind}"
        );
        assert!(!report.complete);
        assert_eq!(report.residue_count, 3);
        report.validate().unwrap();
    }
    for key in ["protein:modifications", "PROTEIN:MODIFIED", "modifications"] {
        let mut molecule = protein("CWK");
        molecule
            .metadata
            .insert(key.into(), "explicit unresolved chemistry".into());
        let report = analyze_protein(&molecule, &Default::default()).unwrap();
        assert!(
            report.mass.is_none()
                && report.net_charge_at_ph.is_none()
                && report.extinction_280.is_none(),
            "{key}"
        );
    }
}
