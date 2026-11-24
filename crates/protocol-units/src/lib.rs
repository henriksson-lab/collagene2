//! Exact, bounded unit-aware arithmetic for cloning protocol calculations.
//!
//! This crate has no sequence model, inventory, network, filesystem or subprocess API.
//! Results describe arithmetic under an explicitly declared model, not reaction success.
//! No calculation silently rounds: choose a rounding mode and inspect the achieved result.
//!
//! ```
//! use collagene_protocol_units::*;
//! let recipe = DilutionRecipe {
//!     schema_version: SCHEMA_VERSION,
//!     stock: Some(MolarConcentration::from_unit(Rational::integer(10), MolarUnit::Micromolar)?),
//!     target: Some(MolarConcentration::from_unit(Rational::integer(100), MolarUnit::Nanomolar)?),
//!     final_volume: Some(Volume::from_unit(Rational::integer(20), VolumeUnit::Microliter)?),
//! };
//! let plan = dilute(&recipe)?;
//! assert_eq!(plan.stock_volume.in_unit(VolumeUnit::Microliter)?, Rational::new(1,5)?);
//! # Ok::<(), Error>(())
//! ```
//! Mass concentration and molar concentration cannot be mixed:
//!
//! ```compile_fail
//! use collagene_protocol_units::*;
//! let recipe = DilutionRecipe {
//!     schema_version: SCHEMA_VERSION,
//!     stock: Some(MassConcentration::from_units(Rational::ONE,
//!         MassUnit::Nanogram, VolumeUnit::Microliter).unwrap()),
//!     target: Some(MolarConcentration::from_unit(Rational::ONE,
//!         MolarUnit::Nanomolar).unwrap()),
//!     final_volume: Some(Volume::from_unit(Rational::ONE,
//!         VolumeUnit::Microliter).unwrap()),
//! };
//! let _ = dilute(&recipe);
//! ```
#![forbid(unsafe_code)]
mod dna;
mod quantities;
mod rational;
mod recipes;
pub use dna::*;
pub use quantities::*;
pub use rational::Rational;
pub use recipes::*;

pub const ENGINE_VERSION: &str = "protocol-units-1";
pub const SCHEMA_VERSION: u16 = 1;
pub const MAX_COMPONENTS: usize = 128;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum InputField {
    StockConcentration,
    TargetConcentration,
    FinalVolume,
    ComponentVolume,
    AnchorMass,
}
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Error {
    Overflow,
    DivisionByZero,
    NegativeQuantity,
    InvalidNumber,
    DecimalPrecisionLimit,
    UnsupportedUnit,
    UnsupportedVersion,
    MissingInput {
        field: InputField,
        index: Option<usize>,
    },
    InvalidRole,
    DuplicateRole,
    TooManyComponents,
    TooManyFactors,
    ZeroStockConcentration,
    ZeroFinalVolume,
    ImpossibleDilution,
    OverfilledMixture,
    RoundingRequired,
    ZeroIncrement,
    RoundedStockExceedsFinal,
    ZeroDnaLength,
    ZeroAnchorMass,
    ZeroMolarRatio,
    InvalidDnaMassModel,
    InvalidResult,
}
impl core::fmt::Display for Error {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::MissingInput { field, index } => {
                write!(f, "missing input {field:?} at component {index:?}")
            }
            other => f.write_str(match other {
                Self::Overflow => "exact result exceeds the bounded rational representation",
                Self::DivisionByZero => "division by zero",
                Self::NegativeQuantity => "quantity would be negative",
                Self::InvalidNumber => "expected a bounded nonnegative decimal",
                Self::DecimalPrecisionLimit => "decimal input exceeds 19 fractional digits",
                Self::UnsupportedUnit => "unsupported or dimensionally incorrect unit",
                Self::UnsupportedVersion => "unsupported recipe or result version",
                Self::InvalidRole => "role must be nonblank, contain 1–128 bytes and have no control characters",
                Self::DuplicateRole => "component role identifiers must be unique",
                Self::TooManyComponents => "at most 128 components are supported",
                Self::TooManyFactors => "at most 16 exact product/division factors are supported",
                Self::ZeroStockConcentration => "stock concentration must be positive",
                Self::ZeroFinalVolume => "final volume must be positive",
                Self::ImpossibleDilution => "dilution target exceeds stock concentration",
                Self::OverfilledMixture => "component volumes exceed the final mixture volume",
                Self::RoundingRequired => {
                    "requested volume is not an exact multiple of the increment"
                }
                Self::ZeroIncrement => "rounding increment must be positive",
                Self::RoundedStockExceedsFinal => "rounded stock volume exceeds the final volume",
                Self::ZeroDnaLength => "DNA length must be positive",
                Self::ZeroAnchorMass => "DNA ratio anchor mass must be positive",
                Self::ZeroMolarRatio => "requested DNA molar ratio must be positive",
                Self::InvalidDnaMassModel => {
                    "declare a positive per-base-pair molar mass, nonblank model identity and assumptions"
                }
                Self::InvalidResult => "stored result disagrees with fresh exact calculation",
                Self::MissingInput { .. } => "missing input",
            }),
        }
    }
}
impl std::error::Error for Error {}
pub(crate) fn version(version: u16) -> Result<(), Error> {
    if version == SCHEMA_VERSION {
        Ok(())
    } else {
        Err(Error::UnsupportedVersion)
    }
}
pub(crate) fn role(role: &str) -> Result<(), Error> {
    if role.trim().is_empty() || role.len() > 128 || role.chars().any(char::is_control) {
        Err(Error::InvalidRole)
    } else {
        Ok(())
    }
}
pub(crate) fn roles<'a>(roles: impl Iterator<Item = &'a str>) -> Result<(), Error> {
    let mut seen = std::collections::BTreeSet::new();
    for r in roles {
        role(r)?;
        if !seen.insert(r) {
            return Err(Error::DuplicateRole);
        }
    }
    Ok(())
}

#[cfg(feature = "serde")]
pub(crate) mod u64_string {
    pub fn serialize<S: serde::Serializer>(value: &u64, s: S) -> Result<S::Ok, S::Error> {
        s.serialize_str(&value.to_string())
    }
    pub fn deserialize<'de, D: serde::Deserializer<'de>>(d: D) -> Result<u64, D::Error> {
        let s = <String as serde::Deserialize>::deserialize(d)?;
        if s.is_empty() || s.len() > 20 || !s.bytes().all(|b| b.is_ascii_digit()) {
            return Err(serde::de::Error::custom("invalid exact integer string"));
        }
        s.parse().map_err(serde::de::Error::custom)
    }
}
