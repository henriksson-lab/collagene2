use collagene_protocol_units::*;
fn r(n: u64, d: u64) -> Rational {
    Rational::new(n, d).unwrap()
}
fn ul(n: u64, d: u64) -> Volume {
    Volume::from_unit(r(n, d), VolumeUnit::Microliter).unwrap()
}
fn nl(n: u64, d: u64) -> Volume {
    Volume::from_unit(r(n, d), VolumeUnit::Nanoliter).unwrap()
}
fn ng(n: u64, d: u64) -> Mass {
    Mass::from_unit(r(n, d), MassUnit::Nanogram).unwrap()
}
fn nm(n: u64, d: u64) -> MolarConcentration {
    MolarConcentration::from_unit(r(n, d), MolarUnit::Nanomolar).unwrap()
}
fn mass_c(n: u64, d: u64) -> MassConcentration {
    MassConcentration::from_units(r(n, d), MassUnit::Nanogram, VolumeUnit::Microliter).unwrap()
}
fn raw_c(value: Rational) -> MolarConcentration {
    MolarConcentration::from_amount_volume(
        Amount::from_unit(value, AmountUnit::Picomole).unwrap(),
        nl(1, 1),
    )
    .unwrap()
}
fn model() -> DuplexDnaMassModel {
    DuplexDnaMassModel{schema_version:SCHEMA_VERSION,model_id:"explicit-660-fixture".into(),grams_per_mole_per_base_pair:Rational::integer(660),assumptions:"Declared approximation for duplex molecules, no chemical modifications or end correction.".into(),source:None}
}
#[test]
fn rational_arithmetic_matches_an_independent_small_integer_oracle() {
    for a in 0u64..=12 {
        for b in 1u64..=12 {
            for c in 0u64..=12 {
                for d in 1u64..=12 {
                    let x = r(a, b);
                    let y = r(c, d);
                    let sum = x.checked_add(y).unwrap();
                    assert_eq!(
                        u128::from(sum.numerator()) * u128::from(b * d),
                        u128::from(sum.denominator()) * u128::from(a * d + c * b)
                    );
                    let product = x.checked_mul(y).unwrap();
                    assert_eq!(
                        u128::from(product.numerator()) * u128::from(b * d),
                        u128::from(product.denominator()) * u128::from(a * c)
                    );
                    assert_eq!(x.cmp(&y), (a * d).cmp(&(c * b)));
                    if a * d >= c * b {
                        assert_eq!(x.checked_sub(y).unwrap().checked_add(y).unwrap(), x);
                    } else {
                        assert_eq!(x.checked_sub(y), Err(Error::NegativeQuantity));
                    }
                    if c != 0 {
                        assert_eq!(x.checked_div(y).unwrap().checked_mul(y).unwrap(), x);
                    }
                }
            }
        }
    }
}
#[test]
fn rational_boundary_and_cross_cancellation_do_not_wrap_or_reject_reducible_results() {
    let max = u64::MAX;
    assert_eq!(
        r(max, 2).checked_add(r(max, 2)).unwrap(),
        Rational::integer(max)
    );
    assert_eq!(
        Rational::integer(max).checked_add(Rational::ONE),
        Err(Error::Overflow)
    );
    assert_eq!(
        r(max, max - 1).checked_mul(r(max - 1, max)).unwrap(),
        Rational::ONE
    );
    assert_eq!(
        Rational::checked_product_ratio(
            &[Rational::integer(max), Rational::integer(max)],
            &[Rational::integer(max)]
        )
        .unwrap(),
        Rational::integer(max)
    );
    assert_eq!(
        Rational::checked_product_ratio(&[], &[Rational::ZERO]),
        Err(Error::DivisionByZero)
    );
    assert_eq!(
        Rational::checked_product_ratio(&[Rational::ZERO], &[Rational::ZERO]),
        Err(Error::DivisionByZero)
    );
    assert_eq!(
        Rational::checked_product_ratio(&[Rational::ONE; 17], &[]),
        Err(Error::TooManyFactors)
    );
    assert_eq!(Rational::new(1, 0), Err(Error::DivisionByZero));
}
#[test]
fn finite_decimal_parse_is_exact_and_bounded() {
    assert_eq!(Rational::parse_decimal("0.00125").unwrap(), r(1, 800));
    assert_eq!(
        Rational::parse_decimal("18446744073709551615.0000000000000000000").unwrap(),
        Rational::integer(u64::MAX)
    );
    assert_eq!(
        Rational::parse_decimal("0.0000000000000000001").unwrap(),
        r(1, 10_000_000_000_000_000_000)
    );
    for invalid in ["NaN", "inf", "1e-3", ".2", "1.2.3", "+1", ""] {
        assert_eq!(Rational::parse_decimal(invalid), Err(Error::InvalidNumber));
    }
    assert_eq!(
        Rational::parse_decimal("-0.1"),
        Err(Error::NegativeQuantity)
    );
    assert_eq!(
        Rational::parse_decimal("0.00000000000000000001"),
        Err(Error::DecimalPrecisionLimit)
    );
}
#[test]
fn declared_units_convert_exactly_and_dimension_spelling_is_strict() {
    assert_eq!(ul(1, 1), nl(1000, 1));
    assert_eq!(
        Volume::from_unit(Rational::integer(1), VolumeUnit::Milliliter).unwrap(),
        ul(1000, 1)
    );
    assert_eq!(
        Mass::from_unit(Rational::integer(1), MassUnit::Microgram).unwrap(),
        ng(1000, 1)
    );
    assert_eq!(
        Amount::from_unit(Rational::integer(1), AmountUnit::Nanomole)
            .unwrap()
            .in_unit(AmountUnit::Picomole)
            .unwrap(),
        Rational::integer(1000)
    );
    assert_eq!(
        MolarConcentration::from_unit(Rational::ONE, MolarUnit::Micromolar).unwrap(),
        nm(1000, 1)
    );
    assert_eq!(VolumeUnit::parse("µL"), Ok(VolumeUnit::Microliter));
    assert_eq!(VolumeUnit::parse("μL"), Ok(VolumeUnit::Microliter));
    assert_eq!(MolarUnit::parse("nm"), Err(Error::UnsupportedUnit));
    assert_eq!(MassUnit::parse("nM"), Err(Error::UnsupportedUnit));
    let c = MassConcentration::from_units(
        Rational::integer(u64::MAX),
        MassUnit::Microgram,
        VolumeUnit::Milliliter,
    )
    .unwrap();
    assert_eq!(
        c.in_units(MassUnit::Microgram, VolumeUnit::Milliliter)
            .unwrap(),
        Rational::integer(u64::MAX)
    );
}
#[test]
fn ordinary_dilutions_and_concentration_amount_relationships_have_literal_answers() {
    let recipe = DilutionRecipe {
        schema_version: 1,
        stock: Some(nm(10_000, 1)),
        target: Some(nm(100, 1)),
        final_volume: Some(ul(20, 1)),
    };
    let p = dilute(&recipe).unwrap();
    assert_eq!(p.stock_volume, ul(1, 5));
    assert_eq!(p.diluent_volume, ul(99, 5));
    p.validate().unwrap();
    assert_eq!(
        nm(100, 1)
            .amount_in(ul(20, 1))
            .unwrap()
            .in_unit(AmountUnit::Picomole)
            .unwrap(),
        Rational::integer(2)
    );
    let mass = dilute(&DilutionRecipe {
        schema_version: 1,
        stock: Some(mass_c(50, 1)),
        target: Some(mass_c(5, 1)),
        final_volume: Some(ul(20, 1)),
    })
    .unwrap();
    assert_eq!(mass.stock_volume, ul(2, 1));
    assert_eq!(mass.diluent_volume, ul(18, 1));
}
#[test]
fn dilution_cancels_the_whole_formula_before_a_nonrepresentable_intermediate_ratio() {
    // x=1/2 nL, y=(M-3)/M nL; V=x+y=3(M-2)/(2M).
    // Reduced V fits u64; (x/V)'s denominator does not. Whole C2*V/C1 still yields x.
    let m = u64::MAX - 1;
    let v = r(3 * ((m - 2) / 4), m / 2);
    let stock = raw_c(v);
    let target = raw_c(r(1, 2));
    assert_eq!(
        target.canonical().checked_div(stock.canonical()),
        Err(Error::Overflow)
    );
    let p = dilute(&DilutionRecipe {
        schema_version: 1,
        stock: Some(stock),
        target: Some(target),
        final_volume: Some(Volume::from_unit(v, VolumeUnit::Nanoliter).unwrap()),
    })
    .unwrap();
    assert_eq!(p.stock_volume, nl(1, 2));
    assert_eq!(p.diluent_volume, nl(m - 3, m));
}
#[test]
fn dilution_missing_zero_and_impossible_inputs_are_distinct() {
    let mut recipe = DilutionRecipe {
        schema_version: 1,
        stock: Some(nm(100, 1)),
        target: Some(nm(10, 1)),
        final_volume: Some(ul(10, 1)),
    };
    recipe.stock = None;
    assert_eq!(
        dilute(&recipe).unwrap_err(),
        Error::MissingInput {
            field: InputField::StockConcentration,
            index: None
        }
    );
    recipe.stock = Some(nm(0, 1));
    assert_eq!(dilute(&recipe).unwrap_err(), Error::ZeroStockConcentration);
    recipe.stock = Some(nm(5, 1));
    assert_eq!(dilute(&recipe).unwrap_err(), Error::ImpossibleDilution);
    recipe.stock = Some(nm(100, 1));
    recipe.final_volume = Some(Volume::ZERO);
    assert_eq!(dilute(&recipe).unwrap_err(), Error::ZeroFinalVolume);
    recipe.schema_version = 2;
    assert_eq!(dilute(&recipe).unwrap_err(), Error::UnsupportedVersion);
}
#[test]
fn dilution_conservation_holds_for_independent_integer_concentration_grid() {
    for stock in 1..=30 {
        for target in 0..=stock {
            for volume in [1, 7, 23] {
                let p = dilute(&DilutionRecipe {
                    schema_version: 1,
                    stock: Some(nm(stock, 1)),
                    target: Some(nm(target, 1)),
                    final_volume: Some(ul(volume, 1)),
                })
                .unwrap();
                assert_eq!(
                    p.stock_volume.checked_add(p.diluent_volume).unwrap(),
                    ul(volume, 1)
                );
                assert_eq!(
                    nm(stock, 1).amount_in(p.stock_volume).unwrap(),
                    nm(target, 1).amount_in(ul(volume, 1)).unwrap()
                );
            }
        }
    }
}
#[test]
fn mixture_balance_is_exact_order_independent_and_rejects_missing_or_overfilled_components() {
    let mut recipe = MixtureRecipe {
        schema_version: 1,
        final_volume: ul(10, 1),
        components: vec![
            MixtureComponent {
                role: "buffer".into(),
                volume: Some(ul(1, 1)),
            },
            MixtureComponent {
                role: "dna".into(),
                volume: Some(ul(3, 2)),
            },
        ],
    };
    let p = balance_mixture(&recipe).unwrap();
    assert_eq!(p.specified_volume, ul(5, 2));
    assert_eq!(p.remaining_diluent, ul(15, 2));
    p.validate().unwrap();
    recipe.components.reverse();
    assert_eq!(
        balance_mixture(&recipe).unwrap().remaining_diluent,
        p.remaining_diluent
    );
    recipe.components[0].volume = None;
    assert_eq!(
        balance_mixture(&recipe).unwrap_err(),
        Error::MissingInput {
            field: InputField::ComponentVolume,
            index: Some(0)
        }
    );
    recipe.components[0].volume = Some(ul(10, 1));
    assert_eq!(
        balance_mixture(&recipe).unwrap_err(),
        Error::OverfilledMixture
    );
    recipe.components[0].role = "buffer".into();
    assert_eq!(balance_mixture(&recipe).unwrap_err(), Error::DuplicateRole);
    recipe.components = vec![
        MixtureComponent {
            role: "x".into(),
            volume: Some(Volume::ZERO)
        };
        129
    ];
    assert_eq!(
        balance_mixture(&recipe).unwrap_err(),
        Error::TooManyComponents
    );
}
#[test]
fn rounding_modes_are_explicit_and_zero_increment_or_unrepresentable_rounding_do_not_panic() {
    assert_eq!(
        round_volume(ul(3, 20), ul(1, 10), RoundingMode::ExactOnly).unwrap_err(),
        Error::RoundingRequired
    );
    let down = round_volume(ul(3, 20), ul(1, 10), RoundingMode::Floor).unwrap();
    assert_eq!(down.rounded, ul(1, 10));
    assert_eq!(down.difference, ul(1, 20));
    assert_eq!(down.direction, RoundingDirection::Down);
    down.validate().unwrap();
    let up = round_volume(ul(3, 20), ul(1, 10), RoundingMode::NearestTiesUp).unwrap();
    assert_eq!(up.rounded, ul(1, 5));
    assert_eq!(up.direction, RoundingDirection::Up);
    assert_eq!(
        round_volume(ul(1, 1), Volume::ZERO, RoundingMode::Floor).unwrap_err(),
        Error::ZeroIncrement
    );
    let huge = nl(u64::MAX, 1);
    let tiny = nl(1, u64::MAX);
    assert_eq!(
        round_volume(huge, tiny, RoundingMode::ExactOnly)
            .unwrap()
            .rounded,
        huge
    );
    assert!(round_volume(nl(u64::MAX, 1), nl(2, 1), RoundingMode::Ceiling).is_err());
}
#[test]
fn rounding_property_grid_matches_integer_tie_oracle() {
    for value in 0u64..100 {
        for step in 1u64..20 {
            let p = round_volume(nl(value, 1), nl(step, 1), RoundingMode::NearestTiesUp).unwrap();
            let q = value / step + u64::from((value % step) * 2 >= step);
            assert_eq!(p.rounded, nl(q * step, 1));
            assert!(p.difference <= nl(step, 2));
            p.validate().unwrap();
        }
    }
}
#[test]
fn rounded_dilution_reports_actual_achieved_concentration_and_cannot_overfill() {
    let recipe = DilutionRecipe {
        schema_version: 1,
        stock: Some(nm(1000, 1)),
        target: Some(nm(15, 1)),
        final_volume: Some(ul(10, 1)),
    };
    let p = quantized_dilution(&recipe, ul(1, 10), RoundingMode::NearestTiesUp).unwrap();
    assert_eq!(p.exact.stock_volume, ul(3, 20));
    assert_eq!(p.stock.rounded, ul(1, 5));
    assert_eq!(p.achieved, nm(20, 1));
    assert!(!p.target_reached);
    assert_eq!(p.diluent_volume, ul(49, 5));
    p.validate().unwrap();
    let recipe = DilutionRecipe {
        schema_version: 1,
        stock: Some(nm(1, 1)),
        target: Some(nm(1, 1)),
        final_volume: Some(nl(1, 1)),
    };
    assert_eq!(
        quantized_dilution(&recipe, nl(2, 1), RoundingMode::Ceiling).unwrap_err(),
        Error::RoundedStockExceedsFinal
    );
}
#[test]
fn dna_mass_amount_and_concentration_fixtures_follow_only_the_declared_coefficient() {
    let model = model();
    let amount = dna_amount(&model, 1000, ng(50, 1)).unwrap();
    assert_eq!(amount.in_unit(AmountUnit::Picomole).unwrap(), r(5, 66));
    assert_eq!(dna_mass(&model, 1000, amount).unwrap(), ng(50, 1));
    assert_eq!(
        dna_molar_concentration(&model, 1000, mass_c(66, 1)).unwrap(),
        nm(100, 1)
    );
    let mut other = model.clone();
    other.grams_per_mole_per_base_pair = Rational::integer(650);
    assert_eq!(
        dna_amount(&other, 1000, ng(50, 1))
            .unwrap()
            .in_unit(AmountUnit::Picomole)
            .unwrap(),
        r(1, 13)
    );
    assert_eq!(
        dna_amount(&model, 0, ng(50, 1)).unwrap_err(),
        Error::ZeroDnaLength
    );
    other.grams_per_mole_per_base_pair = Rational::ZERO;
    assert_eq!(
        dna_amount(&other, 1000, ng(50, 1)).unwrap_err(),
        Error::InvalidDnaMassModel
    );
}
#[test]
fn ratio_plan_is_a_molecule_ratio_not_a_mass_or_end_ratio_and_reports_missing_stock_policy() {
    let recipe = DnaRatioRecipe {
        schema_version: 1,
        model: model(),
        anchor: DnaRatioAnchor {
            role: "vector".into(),
            length_bp: 3000,
            mass: Some(ng(50, 1)),
            stock: Some(mass_c(25, 1)),
        },
        targets: vec![DnaRatioTarget {
            role: "insert".into(),
            length_bp: 1000,
            molar_ratio: Rational::integer(3),
            stock: Some(mass_c(10, 1)),
        }],
        require_stock_volumes: true,
    };
    let p = plan_dna_ratios(&recipe).unwrap();
    assert_eq!(p.anchor.volume, Some(ul(2, 1)));
    assert_eq!(p.targets[0].mass, ng(50, 1));
    assert_eq!(p.targets[0].volume, Some(ul(5, 1)));
    assert_eq!(
        p.targets[0].amount.canonical(),
        p.anchor
            .amount
            .canonical()
            .checked_mul(Rational::integer(3))
            .unwrap()
    );
    p.validate().unwrap();
    assert_eq!(
        dna_molar_ratio(&recipe.model, 1000, ng(50, 1), 3000, ng(50, 1)).unwrap(),
        Rational::integer(3)
    );
    let mut absent = recipe.clone();
    absent.targets[0].stock = None;
    assert_eq!(
        plan_dna_ratios(&absent).unwrap_err(),
        Error::MissingInput {
            field: InputField::StockConcentration,
            index: Some(1)
        }
    );
    absent.require_stock_volumes = false;
    assert_eq!(plan_dna_ratios(&absent).unwrap().targets[0].volume, None);
    absent.anchor.mass = None;
    assert_eq!(
        plan_dna_ratios(&absent).unwrap_err(),
        Error::MissingInput {
            field: InputField::AnchorMass,
            index: Some(0)
        }
    );
}
#[test]
fn exact_dna_inverse_and_ratio_property_grid_has_no_float_tolerance() {
    let model = model();
    for length in [1, 2, 17, 1000, 3000, 1_000_000] {
        for mass in [1, 7, 50, 1000] {
            let original = ng(mass, 1);
            let amount = dna_amount(&model, length, original).unwrap();
            assert_eq!(dna_mass(&model, length, amount).unwrap(), original);
        }
    }
    for a in [100, 1000, 3000] {
        for b in [200, 1000, 5000] {
            let ratio = dna_molar_ratio(&model, a, ng(7, 1), b, ng(11, 1)).unwrap();
            assert_eq!(ratio, r(7 * b, 11 * a));
        }
    }
}
#[test]
fn stored_results_are_not_authority_and_validation_recomputes() {
    let recipe = DilutionRecipe {
        schema_version: 1,
        stock: Some(nm(100, 1)),
        target: Some(nm(10, 1)),
        final_volume: Some(ul(10, 1)),
    };
    let mut plan = dilute(&recipe).unwrap();
    plan.stock_volume = ul(5, 1);
    assert_eq!(plan.validate(), Err(Error::InvalidResult));
    let mut p = round_volume(ul(3, 20), ul(1, 10), RoundingMode::Floor).unwrap();
    p.difference = Volume::ZERO;
    assert_eq!(p.validate(), Err(Error::InvalidResult));
}

#[test]
fn mixture_allows_a_wide_intermediate_prefix_when_the_final_exact_sum_fits() {
    let big = u64::MAX / 2 + 1;
    let components = vec![
        MixtureComponent {
            role: "half1".into(),
            volume: Some(nl(1, 2)),
        },
        MixtureComponent {
            role: "large".into(),
            volume: Some(nl(big, 1)),
        },
        MixtureComponent {
            role: "half2".into(),
            volume: Some(nl(1, 2)),
        },
    ];
    assert_eq!(nl(1, 2).checked_add(nl(big, 1)), Err(Error::Overflow));
    let mut recipe = MixtureRecipe {
        schema_version: 1,
        final_volume: nl(big + 2, 1),
        components,
    };
    let p = balance_mixture(&recipe).unwrap();
    assert_eq!(p.specified_volume, nl(big + 1, 1));
    assert_eq!(p.remaining_diluent, nl(1, 1));
    recipe.components.rotate_left(1);
    assert_eq!(
        balance_mixture(&recipe).unwrap().specified_volume,
        p.specified_volume
    );
}

#[test]
fn declared_dna_model_identity_and_assumptions_cannot_be_whitespace_only() {
    for blank in ["", " ", "\u{2003}"] {
        let mut declared = model();
        declared.model_id = blank.into();
        assert_eq!(declared.validate(), Err(Error::InvalidDnaMassModel));
    }
    for blank in ["", " ", "\n ", "\u{2003}\n"] {
        let mut declared = model();
        declared.assumptions = blank.into();
        assert_eq!(declared.validate(), Err(Error::InvalidDnaMassModel));
        assert_eq!(
            dna_amount(&declared, 1000, ng(50, 1)),
            Err(Error::InvalidDnaMassModel)
        );
    }
    let mut declared = model();
    declared.model_id = " explicit-660-fixture ".into();
    declared.assumptions = "\n Declared duplex molecule approximation. \n".into();
    declared.validate().unwrap();
    assert_eq!(declared.model_id, " explicit-660-fixture ");
    assert_eq!(
        declared.assumptions,
        "\n Declared duplex molecule approximation. \n"
    );
}

#[test]
fn component_role_ids_must_be_nonblank_but_are_not_normalized() {
    for blank in ["", " ", "\u{2003}"] {
        let recipe = MixtureRecipe {
            schema_version: 1,
            final_volume: nl(1, 1),
            components: vec![MixtureComponent {
                role: blank.into(),
                volume: Some(nl(1, 4)),
            }],
        };
        assert_eq!(balance_mixture(&recipe), Err(Error::InvalidRole));
    }
    let recipe = MixtureRecipe {
        schema_version: 1,
        final_volume: nl(1, 1),
        components: vec![
            MixtureComponent {
                role: "insert".into(),
                volume: Some(nl(1, 4)),
            },
            MixtureComponent {
                role: " insert ".into(),
                volume: Some(nl(1, 4)),
            },
        ],
    };
    let plan = balance_mixture(&recipe).unwrap();
    assert_eq!(plan.recipe.components[0].role, "insert");
    assert_eq!(plan.recipe.components[1].role, " insert ");
    assert_eq!(plan.specified_volume, nl(1, 2));
}
