#![cfg(feature = "serde")]
use collagene_protocol_units::*;
#[test]
fn wire_rationals_are_decimal_strings_and_unit_tags_are_enforced() {
    let r = Rational::new(u64::MAX, 7).unwrap();
    let value = serde_json::to_value(r).unwrap();
    assert!(value["numerator"].is_string());
    assert_eq!(serde_json::from_value::<Rational>(value).unwrap(), r);
    for malformed in [
        r#"{"numerator":"1","denominator":"0"}"#,
        r#"{"numerator":"-1","denominator":"1"}"#,
        r#"{"numerator":1,"denominator":"1"}"#,
        r#"{"numerator":"18446744073709551616","denominator":"1"}"#,
        r#"{"numerator":"1","denominator":"1","unexpected":true}"#,
    ] {
        assert!(serde_json::from_str::<Rational>(malformed).is_err());
    }
    let target = DnaRatioTarget {
        role: "large".into(),
        length_bp: u64::MAX,
        molar_ratio: Rational::ONE,
        stock: None,
    };
    let wire = serde_json::to_value(&target).unwrap();
    assert_eq!(wire["length_bp"], u64::MAX.to_string());
    assert_eq!(
        serde_json::from_value::<DnaRatioTarget>(wire).unwrap(),
        target
    );
    let volume = Volume::from_unit(Rational::new(1, 3).unwrap(), VolumeUnit::Microliter).unwrap();
    let value = serde_json::to_value(volume).unwrap();
    assert_eq!(value["unit"], "nL");
    assert_eq!(
        serde_json::from_value::<Volume>(value.clone()).unwrap(),
        volume
    );
    assert!(serde_json::from_value::<Mass>(value).is_err());
}
#[test]
fn recipe_and_results_roundtrip_but_unsupported_versions_or_forged_results_do_not_validate() {
    let recipe = DilutionRecipe {
        schema_version: 1,
        stock: Some(
            MolarConcentration::from_unit(Rational::integer(1000), MolarUnit::Nanomolar).unwrap(),
        ),
        target: Some(
            MolarConcentration::from_unit(Rational::integer(15), MolarUnit::Nanomolar).unwrap(),
        ),
        final_volume: Some(
            Volume::from_unit(Rational::integer(10), VolumeUnit::Microliter).unwrap(),
        ),
    };
    let plan = quantized_dilution(
        &recipe,
        Volume::from_unit(Rational::new(1, 10).unwrap(), VolumeUnit::Microliter).unwrap(),
        RoundingMode::NearestTiesUp,
    )
    .unwrap();
    let wire = serde_json::to_string(&plan).unwrap();
    let loaded: QuantizedDilution<MolarConcentration> = serde_json::from_str(&wire).unwrap();
    loaded.validate().unwrap();
    assert_eq!(loaded, plan);
    let mut forged = serde_json::to_value(&plan).unwrap();
    forged["target_reached"] = serde_json::json!(true);
    let forged: QuantizedDilution<MolarConcentration> = serde_json::from_value(forged).unwrap();
    assert_eq!(forged.validate(), Err(Error::InvalidResult));
    let mut bad = serde_json::to_value(&recipe).unwrap();
    bad["schema_version"] = serde_json::json!(2);
    let bad: DilutionRecipe<MolarConcentration> = serde_json::from_value(bad).unwrap();
    assert_eq!(dilute(&bad).unwrap_err(), Error::UnsupportedVersion);
    let mut unknown = serde_json::to_value(&recipe).unwrap();
    unknown["inventory_stock"] = serde_json::json!("outside scope");
    assert!(serde_json::from_value::<DilutionRecipe<MolarConcentration>>(unknown).is_err());
}
#[test]
fn dna_recipe_keeps_declared_model_roles_and_all_input_assumptions() {
    let recipe = DnaRatioRecipe {
        schema_version: 1,
        model: DuplexDnaMassModel {
            schema_version: 1,
            model_id: "explicit650".into(),
            grams_per_mole_per_base_pair: Rational::integer(650),
            assumptions: "Declared duplex molecule approximation; no end correction.".into(),
            source: Some(
                "https://www.neb.com/tools-and-resources/usage-guidelines/nucleic-acid-data".into(),
            ),
        },
        anchor: DnaRatioAnchor {
            role: "vector".into(),
            length_bp: 3000,
            mass: Some(Mass::from_unit(Rational::integer(50), MassUnit::Nanogram).unwrap()),
            stock: None,
        },
        targets: vec![DnaRatioTarget {
            role: "insert".into(),
            length_bp: 1000,
            molar_ratio: Rational::integer(3),
            stock: None,
        }],
        require_stock_volumes: false,
    };
    let plan = plan_dna_ratios(&recipe).unwrap();
    let wire = serde_json::to_string(&plan).unwrap();
    let loaded: DnaRatioPlan = serde_json::from_str(&wire).unwrap();
    loaded.validate().unwrap();
    assert_eq!(loaded.recipe.model.model_id, "explicit650");
    assert_eq!(
        loaded.recipe.model.grams_per_mole_per_base_pair,
        Rational::integer(650)
    );
    assert_eq!(loaded, plan);
}
