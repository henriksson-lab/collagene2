use crate::*;
macro_rules! dto {($($item:tt)*)=>{#[derive(Debug,Clone,PartialEq,Eq)]#[cfg_attr(feature="serde",derive(serde::Serialize,serde::Deserialize))]#[cfg_attr(feature="serde",serde(deny_unknown_fields))]$($item)*};}
// A caller-declared scalar approximation for duplex DNA molecule molar mass.
// No default coefficient, sequence composition, chemical modification or end model is inferred.
dto! {pub struct DuplexDnaMassModel {pub schema_version:u16,pub model_id:String,pub grams_per_mole_per_base_pair:Rational,pub assumptions:String,pub source:Option<String>}}
impl DuplexDnaMassModel {
    pub fn validate(&self) -> Result<(), Error> {
        version(self.schema_version)?;
        if self.model_id.trim().is_empty()
            || self.model_id.len() > 128
            || self.model_id.chars().any(char::is_control)
            || self.grams_per_mole_per_base_pair.is_zero()
            || self.assumptions.trim().is_empty()
            || self.assumptions.len() > 2048
            || self
                .assumptions
                .chars()
                .any(|c| c.is_control() && c != '\n')
            || self
                .source
                .as_ref()
                .is_some_and(|s| s.len() > 2048 || s.chars().any(char::is_control))
        {
            return Err(Error::InvalidDnaMassModel);
        }
        Ok(())
    }
}
fn length(n: u64) -> Result<Rational, Error> {
    if n == 0 {
        Err(Error::ZeroDnaLength)
    } else {
        Ok(Rational::integer(n))
    }
}
/// pmol of duplex molecules, not pmol of ends. mass[ng]*1000/(length[bp]*MW[g/mol/bp]).
pub fn dna_amount(model: &DuplexDnaMassModel, length_bp: u64, mass: Mass) -> Result<Amount, Error> {
    model.validate()?;
    let n = length(length_bp)?;
    Amount::from_unit(
        Rational::checked_product_ratio(
            &[mass.canonical(), Rational::integer(1000)],
            &[n, model.grams_per_mole_per_base_pair],
        )?,
        AmountUnit::Picomole,
    )
}
/// Inverse of `dna_amount` under exactly the same explicitly declared mass model.
pub fn dna_mass(model: &DuplexDnaMassModel, length_bp: u64, amount: Amount) -> Result<Mass, Error> {
    model.validate()?;
    let n = length(length_bp)?;
    Mass::from_unit(
        Rational::checked_product_ratio(
            &[amount.canonical(), n, model.grams_per_mole_per_base_pair],
            &[Rational::integer(1000)],
        )?,
        MassUnit::Nanogram,
    )
}
pub fn dna_molar_concentration(
    model: &DuplexDnaMassModel,
    length_bp: u64,
    concentration: MassConcentration,
) -> Result<MolarConcentration, Error> {
    model.validate()?;
    let n = length(length_bp)?;
    let per_nl = Rational::checked_product_ratio(
        &[concentration.canonical(), Rational::integer(1000)],
        &[n, model.grams_per_mole_per_base_pair],
    )?;
    MolarConcentration::from_amount_volume(
        Amount::from_unit(per_nl, AmountUnit::Picomole)?,
        Volume::from_unit(Rational::ONE, VolumeUnit::Nanoliter)?,
    )
}
/// Molecule ratio of numerator DNA to denominator DNA under one shared duplex mass model.
pub fn dna_molar_ratio(
    model: &DuplexDnaMassModel,
    numerator_length: u64,
    numerator_mass: Mass,
    denominator_length: u64,
    denominator_mass: Mass,
) -> Result<Rational, Error> {
    model.validate()?;
    let a = length(numerator_length)?;
    let b = length(denominator_length)?;
    if denominator_mass.is_zero() {
        return Err(Error::ZeroAnchorMass);
    }
    Rational::checked_product_ratio(
        &[numerator_mass.canonical(), b],
        &[denominator_mass.canonical(), a],
    )
}
dto! {pub struct DnaRatioAnchor {pub role:String,#[cfg_attr(feature="serde",serde(with="crate::u64_string"))] pub length_bp:u64,pub mass:Option<Mass>,pub stock:Option<MassConcentration>}}
dto! {pub struct DnaRatioTarget {pub role:String,#[cfg_attr(feature="serde",serde(with="crate::u64_string"))] pub length_bp:u64,pub molar_ratio:Rational,pub stock:Option<MassConcentration>}}
dto! {pub struct DnaRatioRecipe {pub schema_version:u16,pub model:DuplexDnaMassModel,pub anchor:DnaRatioAnchor,pub targets:Vec<DnaRatioTarget>,pub require_stock_volumes:bool}}
dto! {pub struct DnaRatioComponent {pub role:String,#[cfg_attr(feature="serde",serde(with="crate::u64_string"))] pub length_bp:u64,pub molar_ratio:Rational,pub mass:Mass,pub amount:Amount,pub volume:Option<Volume>}}
dto! {pub struct DnaRatioPlan {pub schema_version:u16,pub engine_version:String,pub recipe:DnaRatioRecipe,pub anchor:DnaRatioComponent,pub targets:Vec<DnaRatioComponent>,pub assumptions:Vec<String>}}
fn stock_volume(
    mass: Mass,
    stock: Option<MassConcentration>,
    required: bool,
    index: usize,
) -> Result<Option<Volume>, Error> {
    match stock {
        Some(c) => {
            if c.is_zero() {
                return Err(Error::ZeroStockConcentration);
            }
            Ok(Some(Volume::from_unit(
                mass.canonical().checked_div(c.canonical())?,
                VolumeUnit::Nanoliter,
            )?))
        }
        None if required => Err(Error::MissingInput {
            field: InputField::StockConcentration,
            index: Some(index),
        }),
        None => Ok(None),
    }
}
/// Plan target molecule ratios against one explicitly weighed anchor. No recommended ratio is implied.
pub fn plan_dna_ratios(recipe: &DnaRatioRecipe) -> Result<DnaRatioPlan, Error> {
    version(recipe.schema_version)?;
    recipe.model.validate()?;
    if recipe.targets.len() >= MAX_COMPONENTS {
        return Err(Error::TooManyComponents);
    }
    roles(
        std::iter::once(recipe.anchor.role.as_str())
            .chain(recipe.targets.iter().map(|c| c.role.as_str())),
    )?;
    let mass = recipe.anchor.mass.ok_or(Error::MissingInput {
        field: InputField::AnchorMass,
        index: Some(0),
    })?;
    if mass.is_zero() {
        return Err(Error::ZeroAnchorMass);
    }
    let amount = dna_amount(&recipe.model, recipe.anchor.length_bp, mass)?;
    let anchor = DnaRatioComponent {
        role: recipe.anchor.role.clone(),
        length_bp: recipe.anchor.length_bp,
        molar_ratio: Rational::ONE,
        mass,
        amount,
        volume: stock_volume(mass, recipe.anchor.stock, recipe.require_stock_volumes, 0)?,
    };
    let mut targets = Vec::with_capacity(recipe.targets.len());
    for (i, t) in recipe.targets.iter().enumerate() {
        if t.molar_ratio.is_zero() {
            return Err(Error::ZeroMolarRatio);
        }
        let n = length(t.length_bp)?;
        let target_mass = Mass::from_unit(
            Rational::checked_product_ratio(
                &[mass.canonical(), t.molar_ratio, n],
                &[length(recipe.anchor.length_bp)?],
            )?,
            MassUnit::Nanogram,
        )?;
        let amount = dna_amount(&recipe.model, t.length_bp, target_mass)?;
        targets.push(DnaRatioComponent {
            role: t.role.clone(),
            length_bp: t.length_bp,
            molar_ratio: t.molar_ratio,
            mass: target_mass,
            amount,
            volume: stock_volume(target_mass, t.stock, recipe.require_stock_volumes, i + 1)?,
        });
    }
    Ok(DnaRatioPlan{schema_version:SCHEMA_VERSION,engine_version:ENGINE_VERSION.into(),recipe:recipe.clone(),anchor,targets,assumptions:vec!["Volumes are additive; optional stock concentrations are declared inputs, not measured or inventory-verified quantities.".into(),"Amounts count duplex DNA molecules, not free ends. Length-only molar mass is the explicitly declared approximation; sequence composition, chemical modifications and overhang/end states are not modeled. No reaction success is predicted.".into()]})
}
impl DnaRatioPlan {
    pub fn validate(&self) -> Result<(), Error> {
        version(self.schema_version)?;
        if self.engine_version != ENGINE_VERSION || *self != plan_dna_ratios(&self.recipe)? {
            return Err(Error::InvalidResult);
        }
        Ok(())
    }
}
