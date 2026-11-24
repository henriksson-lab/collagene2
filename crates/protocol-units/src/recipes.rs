use crate::*;
macro_rules! dto {($($item:tt)*)=>{#[derive(Debug,Clone,PartialEq,Eq)]#[cfg_attr(feature="serde",derive(serde::Serialize,serde::Deserialize))]#[cfg_attr(feature="serde",serde(deny_unknown_fields))]$($item)*};}
dto! {pub struct DilutionRecipe<C> {pub schema_version:u16,pub stock:Option<C>,pub target:Option<C>,pub final_volume:Option<Volume>}}
dto! {pub struct DilutionPlan<C> {pub schema_version:u16,pub engine_version:String,pub recipe:DilutionRecipe<C>,pub stock_volume:Volume,pub diluent_volume:Volume}}
/// C1*V1=C2*V2, additive volumes, one stock and zero-concentration diluent.
/// This does not assert solvent compatibility, volume contraction or reaction success.
pub fn dilute<C: Concentration>(recipe: &DilutionRecipe<C>) -> Result<DilutionPlan<C>, Error> {
    version(recipe.schema_version)?;
    let missing = |field| Error::MissingInput { field, index: None };
    let stock = recipe
        .stock
        .ok_or_else(|| missing(InputField::StockConcentration))?;
    let target = recipe
        .target
        .ok_or_else(|| missing(InputField::TargetConcentration))?;
    let final_volume = recipe
        .final_volume
        .ok_or_else(|| missing(InputField::FinalVolume))?;
    if stock.exact().is_zero() {
        return Err(Error::ZeroStockConcentration);
    }
    if final_volume.is_zero() {
        return Err(Error::ZeroFinalVolume);
    }
    if target.exact() > stock.exact() {
        return Err(Error::ImpossibleDilution);
    }
    let stock_volume = Volume::from_unit(
        Rational::checked_product_ratio(
            &[target.exact(), final_volume.canonical()],
            &[stock.exact()],
        )?,
        VolumeUnit::Nanoliter,
    )?;
    Ok(DilutionPlan {
        schema_version: SCHEMA_VERSION,
        engine_version: ENGINE_VERSION.into(),
        recipe: recipe.clone(),
        stock_volume,
        diluent_volume: final_volume.checked_sub(stock_volume)?,
    })
}
impl<C: Concentration> DilutionPlan<C> {
    pub fn validate(&self) -> Result<(), Error> {
        version(self.schema_version)?;
        if self.engine_version != ENGINE_VERSION || *self != dilute(&self.recipe)? {
            return Err(Error::InvalidResult);
        }
        Ok(())
    }
}
dto! {pub struct MixtureComponent {pub role:String,pub volume:Option<Volume>}}
dto! {pub struct MixtureRecipe {pub schema_version:u16,pub final_volume:Volume,pub components:Vec<MixtureComponent>}}
dto! {pub struct MixturePlan {pub schema_version:u16,pub engine_version:String,pub recipe:MixtureRecipe,pub specified_volume:Volume,pub remaining_diluent:Volume}}
/// Exact volume balance for at most 128 explicitly named components plus one fill diluent.
pub fn balance_mixture(recipe: &MixtureRecipe) -> Result<MixturePlan, Error> {
    version(recipe.schema_version)?;
    if recipe.final_volume.is_zero() {
        return Err(Error::ZeroFinalVolume);
    }
    if recipe.components.len() > MAX_COMPONENTS {
        return Err(Error::TooManyComponents);
    }
    roles(recipe.components.iter().map(|c| c.role.as_str()))?;
    let values: Vec<Rational> = recipe
        .components
        .iter()
        .enumerate()
        .map(|(i, c)| {
            c.volume.map(Volume::canonical).ok_or(Error::MissingInput {
                field: InputField::ComponentVolume,
                index: Some(i),
            })
        })
        .collect::<Result<_, _>>()?;
    let total = Volume::from_unit(Rational::checked_sum(&values)?, VolumeUnit::Nanoliter)?;
    if total > recipe.final_volume {
        return Err(Error::OverfilledMixture);
    }
    Ok(MixturePlan {
        schema_version: SCHEMA_VERSION,
        engine_version: ENGINE_VERSION.into(),
        recipe: recipe.clone(),
        specified_volume: total,
        remaining_diluent: recipe.final_volume.checked_sub(total)?,
    })
}
impl MixturePlan {
    pub fn validate(&self) -> Result<(), Error> {
        version(self.schema_version)?;
        if self.engine_version != ENGINE_VERSION || *self != balance_mixture(&self.recipe)? {
            return Err(Error::InvalidResult);
        }
        Ok(())
    }
}
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub enum RoundingMode {
    ExactOnly,
    Floor,
    Ceiling,
    NearestTiesUp,
}
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub enum RoundingDirection {
    Exact,
    Down,
    Up,
}
dto! {pub struct RoundedVolume {pub schema_version:u16,pub engine_version:String,pub requested:Volume,pub increment:Volume,pub mode:RoundingMode,pub rounded:Volume,pub direction:RoundingDirection,pub difference:Volume}}
/// Quantize explicitly; difference is the nonnegative absolute deviation, never a confidence interval.
pub fn round_volume(
    value: Volume,
    increment: Volume,
    mode: RoundingMode,
) -> Result<RoundedVolume, Error> {
    if increment.is_zero() {
        return Err(Error::ZeroIncrement);
    }
    let v = value.canonical();
    let step = increment.canonical();
    let n = u128::from(v.numerator()) * u128::from(step.denominator());
    let d = u128::from(v.denominator()) * u128::from(step.numerator());
    let floor = n / d;
    let rem = n % d;
    let up = match mode {
        RoundingMode::ExactOnly if rem != 0 => return Err(Error::RoundingRequired),
        RoundingMode::ExactOnly | RoundingMode::Floor => false,
        RoundingMode::Ceiling => rem != 0,
        RoundingMode::NearestTiesUp => rem != 0 && rem >= d - rem,
    };
    let ticks = floor.checked_add(u128::from(up)).ok_or(Error::Overflow)?;
    let count = Rational::wide(
        ticks
            .checked_mul(u128::from(step.numerator()))
            .ok_or(Error::Overflow)?,
        u128::from(step.denominator()),
    )?;
    let rounded = Volume::from_unit(count, VolumeUnit::Nanoliter)?;
    let (direction, difference) = if rounded == value {
        (RoundingDirection::Exact, Volume::ZERO)
    } else if rounded < value {
        (RoundingDirection::Down, value.checked_sub(rounded)?)
    } else {
        (RoundingDirection::Up, rounded.checked_sub(value)?)
    };
    Ok(RoundedVolume {
        schema_version: SCHEMA_VERSION,
        engine_version: ENGINE_VERSION.into(),
        requested: value,
        increment,
        mode,
        rounded,
        direction,
        difference,
    })
}
impl RoundedVolume {
    pub fn validate(&self) -> Result<(), Error> {
        version(self.schema_version)?;
        if self.engine_version != ENGINE_VERSION
            || *self != round_volume(self.requested, self.increment, self.mode)?
        {
            return Err(Error::InvalidResult);
        }
        Ok(())
    }
}
dto! {pub struct QuantizedDilution<C> {pub schema_version:u16,pub engine_version:String,pub exact:DilutionPlan<C>,pub stock:RoundedVolume,pub diluent_volume:Volume,pub achieved:C,pub target_reached:bool}}
/// Round stock volume, refill to the original total, and calculate the actual achieved concentration.
pub fn quantized_dilution<C: Concentration>(
    recipe: &DilutionRecipe<C>,
    increment: Volume,
    mode: RoundingMode,
) -> Result<QuantizedDilution<C>, Error> {
    let exact = dilute(recipe)?;
    let stock = round_volume(exact.stock_volume, increment, mode)?;
    let final_volume = recipe.final_volume.ok_or(Error::MissingInput {
        field: InputField::FinalVolume,
        index: None,
    })?;
    if stock.rounded > final_volume {
        return Err(Error::RoundedStockExceedsFinal);
    }
    let input = recipe.stock.ok_or(Error::MissingInput {
        field: InputField::StockConcentration,
        index: None,
    })?;
    let achieved = C::with_exact(Rational::checked_product_ratio(
        &[input.exact(), stock.rounded.canonical()],
        &[final_volume.canonical()],
    )?);
    let target = recipe.target.ok_or(Error::MissingInput {
        field: InputField::TargetConcentration,
        index: None,
    })?;
    let diluent_volume = final_volume.checked_sub(stock.rounded)?;
    Ok(QuantizedDilution {
        schema_version: SCHEMA_VERSION,
        engine_version: ENGINE_VERSION.into(),
        exact,
        stock,
        diluent_volume,
        achieved,
        target_reached: achieved == target,
    })
}
impl<C: Concentration> QuantizedDilution<C> {
    pub fn validate(&self) -> Result<(), Error> {
        version(self.schema_version)?;
        if self.engine_version != ENGINE_VERSION
            || *self
                != quantized_dilution(&self.exact.recipe, self.stock.increment, self.stock.mode)?
        {
            return Err(Error::InvalidResult);
        }
        Ok(())
    }
}
