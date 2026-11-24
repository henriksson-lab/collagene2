//! Decimal-input MCP adapters; all arithmetic belongs to protocol-units.
//! These calculations neither inspect stock inventory nor verify reaction chemistry.
use collagene_protocol_units as units;
use rmcp::ErrorData;
use schemars::JsonSchema;
use serde::Deserialize;
use serde_json::{Value, json};

#[derive(Debug, Clone, Copy, Deserialize, JsonSchema)]
pub enum VolumeUnit {
    #[serde(rename = "nL")]
    Nanoliter,
    #[serde(rename = "uL")]
    Microliter,
    #[serde(rename = "mL")]
    Milliliter,
}
impl From<VolumeUnit> for units::VolumeUnit {
    fn from(value: VolumeUnit) -> Self {
        match value {
            VolumeUnit::Nanoliter => Self::Nanoliter,
            VolumeUnit::Microliter => Self::Microliter,
            VolumeUnit::Milliliter => Self::Milliliter,
        }
    }
}
#[derive(Debug, Clone, Copy, Deserialize, JsonSchema)]
pub enum MassUnit {
    #[serde(rename = "ng")]
    Nanogram,
    #[serde(rename = "ug")]
    Microgram,
}
impl From<MassUnit> for units::MassUnit {
    fn from(value: MassUnit) -> Self {
        match value {
            MassUnit::Nanogram => Self::Nanogram,
            MassUnit::Microgram => Self::Microgram,
        }
    }
}
#[derive(Debug, Clone, Copy, Deserialize, JsonSchema)]
pub enum MolarUnit {
    #[serde(rename = "nM")]
    Nanomolar,
    #[serde(rename = "uM")]
    Micromolar,
}
impl From<MolarUnit> for units::MolarUnit {
    fn from(value: MolarUnit) -> Self {
        match value {
            MolarUnit::Nanomolar => Self::Nanomolar,
            MolarUnit::Micromolar => Self::Micromolar,
        }
    }
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct VolumeInput {
    /// Exact nonnegative decimal string, at most 64 bytes/19 fractional digits; no exponent notation.
    pub value: String,
    pub unit: VolumeUnit,
}
impl VolumeInput {
    fn parse(self) -> Result<units::Volume, units::Error> {
        units::Volume::from_unit(
            units::Rational::parse_decimal(&self.value)?,
            self.unit.into(),
        )
    }
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct MassInput {
    /// Exact nonnegative decimal string, at most 64 bytes/19 fractional digits; no exponent notation.
    pub value: String,
    pub unit: MassUnit,
}
impl MassInput {
    fn parse(self) -> Result<units::Mass, units::Error> {
        units::Mass::from_unit(
            units::Rational::parse_decimal(&self.value)?,
            self.unit.into(),
        )
    }
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct MassConcentrationInput {
    /// Exact nonnegative decimal string for mass_unit/volume_unit.
    pub value: String,
    pub mass_unit: MassUnit,
    pub volume_unit: VolumeUnit,
}
impl MassConcentrationInput {
    fn parse(self) -> Result<units::MassConcentration, units::Error> {
        units::MassConcentration::from_units(
            units::Rational::parse_decimal(&self.value)?,
            self.mass_unit.into(),
            self.volume_unit.into(),
        )
    }
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(tag = "kind", deny_unknown_fields)]
pub enum ConcentrationInput {
    #[serde(rename = "molar")]
    Molar { value: String, unit: MolarUnit },
    #[serde(rename = "mass")]
    Mass {
        value: String,
        mass_unit: MassUnit,
        volume_unit: VolumeUnit,
    },
}
#[derive(Debug, Clone, Copy, Deserialize, JsonSchema)]
pub enum RoundingMode {
    ExactOnly,
    Floor,
    Ceiling,
    NearestTiesUp,
}
impl From<RoundingMode> for units::RoundingMode {
    fn from(value: RoundingMode) -> Self {
        match value {
            RoundingMode::ExactOnly => Self::ExactOnly,
            RoundingMode::Floor => Self::Floor,
            RoundingMode::Ceiling => Self::Ceiling,
            RoundingMode::NearestTiesUp => Self::NearestTiesUp,
        }
    }
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct QuantizationInput {
    pub increment: VolumeInput,
    /// Explicit quantization; omitted quantization never rounds. Difference is absolute volume deviation.
    pub mode: RoundingMode,
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct DilutionQuery {
    /// Declared stock; mass and molar concentration dimensions cannot be mixed.
    pub stock: Option<ConcentrationInput>,
    pub target: Option<ConcentrationInput>,
    pub final_volume: Option<VolumeInput>,
    pub quantization: Option<QuantizationInput>,
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct MixtureComponentInput {
    /// Unique descriptive identifier (1–128 bytes); not an inventory reference.
    pub role: String,
    pub volume: Option<VolumeInput>,
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct MixtureQuery {
    pub final_volume: VolumeInput,
    /// At most 128 named component volumes; the remaining volume is one unspecified fill diluent.
    pub components: Vec<MixtureComponentInput>,
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct DnaMassModelInput {
    /// Caller-defined model identity, 1–128 bytes. No default coefficient is supplied.
    pub model_id: String,
    /// Explicit positive exact decimal coefficient in g/mol/base pair for duplex molecules.
    pub grams_per_mole_per_base_pair: String,
    /// Required statement of the model approximation, at most 2048bytes.
    pub assumptions: String,
    /// Optional descriptive citation, at most 2048bytes; never fetched or authenticated.
    pub source: Option<String>,
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct DnaAnchorInput {
    pub role: String,
    /// Exact positive whole base-pair count as a string, at most 20 digits.
    pub length_bp: String,
    pub mass: Option<MassInput>,
    pub stock: Option<MassConcentrationInput>,
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct DnaTargetInput {
    pub role: String,
    pub length_bp: String,
    /// Explicit positive decimal molecule ratio relative to the anchor; no recommended ratio is inferred.
    pub molar_ratio: String,
    pub stock: Option<MassConcentrationInput>,
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct DnaRatioQuery {
    /// Required caller declaration; no implicit 660 g/mol/bp or empirical mass model.
    pub model: DnaMassModelInput,
    pub anchor: DnaAnchorInput,
    /// At most 127 targets plus the anchor, all with unique descriptive roles.
    pub targets: Vec<DnaTargetInput>,
    /// If true, absent stock concentrations are errors rather than invented volumes.
    pub require_stock_volumes: bool,
}

fn failure(error: units::Error) -> ErrorData {
    let details = match &error {
        units::Error::MissingInput { field, index } => {
            json!({"domain":"protocol_units","code":"MissingInput","field":format!("{field:?}"),"component_index":index})
        }
        _ => json!({"domain":"protocol_units","code":format!("{error:?}")}),
    };
    ErrorData::invalid_params(error.to_string(), Some(details))
}
fn dilution<C: units::Concentration + serde::Serialize>(
    recipe: units::DilutionRecipe<C>,
    quantization: Option<QuantizationInput>,
) -> Result<Value, units::Error> {
    let exact = units::dilute(&recipe)?;
    let quantized = quantization
        .map(|q| units::quantized_dilution(&recipe, q.increment.parse()?, q.mode.into()))
        .transpose()?;
    let achieved = quantized.as_ref().map(|q| q.achieved).or(recipe.target);
    let target_reached = quantized.as_ref().is_none_or(|q| q.target_reached);
    Ok(
        json!({"engine_version":units::ENGINE_VERSION,"schema_version":units::SCHEMA_VERSION,
        "exact":exact,"quantized":quantized,"achieved_concentration":achieved,"target_reached":target_reached,
        "assumptions":["C1*V1=C2*V2 with additive volumes and one zero-concentration diluent.","Only explicitly requested stock-volume quantization is applied; diluent refills to the declared final volume. Inspect achieved concentration and absolute rounding deviation.","Declared concentrations are inputs, not measured or inventory-verified quantities. Solvent compatibility, volume contraction and reaction success are not assessed."]}),
    )
}
fn molar(
    value: Option<ConcentrationInput>,
) -> Result<Option<units::MolarConcentration>, units::Error> {
    value
        .map(|value| match value {
            ConcentrationInput::Molar { value, unit } => units::MolarConcentration::from_unit(
                units::Rational::parse_decimal(&value)?,
                unit.into(),
            ),
            _ => Err(units::Error::UnsupportedUnit),
        })
        .transpose()
}
fn mass(
    value: Option<ConcentrationInput>,
) -> Result<Option<units::MassConcentration>, units::Error> {
    value
        .map(|value| match value {
            ConcentrationInput::Mass {
                value,
                mass_unit,
                volume_unit,
            } => MassConcentrationInput {
                value,
                mass_unit,
                volume_unit,
            }
            .parse(),
            _ => Err(units::Error::UnsupportedUnit),
        })
        .transpose()
}
pub fn calculate_dilution(args: DilutionQuery) -> Result<Value, ErrorData> {
    (|| {
        let final_volume = args.final_volume.map(VolumeInput::parse).transpose()?;
        let is_mass = matches!(
            args.stock.as_ref().or(args.target.as_ref()),
            Some(ConcentrationInput::Mass { .. })
        );
        if is_mass {
            dilution(
                units::DilutionRecipe {
                    schema_version: units::SCHEMA_VERSION,
                    stock: mass(args.stock)?,
                    target: mass(args.target)?,
                    final_volume,
                },
                args.quantization,
            )
        } else {
            dilution(
                units::DilutionRecipe {
                    schema_version: units::SCHEMA_VERSION,
                    stock: molar(args.stock)?,
                    target: molar(args.target)?,
                    final_volume,
                },
                args.quantization,
            )
        }
    })()
    .map_err(failure)
}
pub fn balance_mixture(args: MixtureQuery) -> Result<Value, ErrorData> {
    (|| {
        if args.components.len()>units::MAX_COMPONENTS{return Err(units::Error::TooManyComponents)}
        let final_volume=args.final_volume.parse()?;
        let components=args.components.into_iter().map(|c|Ok(units::MixtureComponent{role:c.role,volume:c.volume.map(VolumeInput::parse).transpose()?})).collect::<Result<Vec<_>,units::Error>>()?;
        let plan=units::balance_mixture(&units::MixtureRecipe{schema_version:units::SCHEMA_VERSION,final_volume,components})?;
        Ok(json!({"plan":plan,"assumptions":["Exact additive volume balance for named declared components and one remaining fill diluent.","No concentration or chemical-compatibility calculation, inventory lookup or reaction verification is implied."]}))
    })().map_err(failure)
}
fn length(s: &str) -> Result<u64, units::Error> {
    if s.is_empty() || s.len() > 20 || !s.bytes().all(|b| b.is_ascii_digit()) {
        return Err(units::Error::InvalidNumber);
    }
    s.parse().map_err(|_| units::Error::Overflow)
}
pub fn calculate_dna_molar_ratios(args: DnaRatioQuery) -> Result<Value, ErrorData> {
    (|| {
        if args.targets.len()>=units::MAX_COMPONENTS{return Err(units::Error::TooManyComponents)}
        let m=args.model;
        let model=units::DuplexDnaMassModel{schema_version:units::SCHEMA_VERSION,model_id:m.model_id,grams_per_mole_per_base_pair:units::Rational::parse_decimal(&m.grams_per_mole_per_base_pair)?,assumptions:m.assumptions,source:m.source};
        let a=args.anchor;
        let anchor=units::DnaRatioAnchor{role:a.role,length_bp:length(&a.length_bp)?,mass:a.mass.map(MassInput::parse).transpose()?,stock:a.stock.map(MassConcentrationInput::parse).transpose()?};
        let targets=args.targets.into_iter().map(|t|Ok(units::DnaRatioTarget{role:t.role,length_bp:length(&t.length_bp)?,molar_ratio:units::Rational::parse_decimal(&t.molar_ratio)?,stock:t.stock.map(MassConcentrationInput::parse).transpose()?})).collect::<Result<Vec<_>,units::Error>>()?;
        let plan=units::plan_dna_ratios(&units::DnaRatioRecipe{schema_version:units::SCHEMA_VERSION,model,anchor,targets,require_stock_volumes:args.require_stock_volumes})?;
        Ok(json!({"plan":plan,"authority":"Exact arithmetic under the caller-declared duplex mass approximation; no inferred optimal ratio, stock inventory, end chemistry or biological reaction success."}))
    })().map_err(failure)
}
