use crate::{Error, Rational};
macro_rules! unit {($name:ident {$($variant:ident),+})=>{#[derive(Debug,Clone,Copy,PartialEq,Eq)]#[cfg_attr(feature="serde",derive(serde::Serialize,serde::Deserialize))]pub enum $name {$($variant),+}};}
unit!(VolumeUnit {
    Nanoliter,
    Microliter,
    Milliliter
});
unit!(MassUnit {
    Nanogram,
    Microgram
});
unit!(AmountUnit { Picomole, Nanomole });
unit!(MolarUnit {
    Nanomolar,
    Micromolar
});
impl VolumeUnit {
    pub fn parse(s: &str) -> Result<Self, Error> {
        match s {
            "nL" => Ok(Self::Nanoliter),
            "uL" | "µL" | "μL" => Ok(Self::Microliter),
            "mL" => Ok(Self::Milliliter),
            _ => Err(Error::UnsupportedUnit),
        }
    }
    fn scale(self) -> Rational {
        Rational::integer(match self {
            Self::Nanoliter => 1,
            Self::Microliter => 1000,
            Self::Milliliter => 1_000_000,
        })
    }
}
impl MassUnit {
    pub fn parse(s: &str) -> Result<Self, Error> {
        match s {
            "ng" => Ok(Self::Nanogram),
            "ug" | "µg" | "μg" => Ok(Self::Microgram),
            _ => Err(Error::UnsupportedUnit),
        }
    }
    fn scale(self) -> Rational {
        Rational::integer(match self {
            Self::Nanogram => 1,
            Self::Microgram => 1000,
        })
    }
}
impl AmountUnit {
    pub fn parse(s: &str) -> Result<Self, Error> {
        match s {
            "pmol" => Ok(Self::Picomole),
            "nmol" => Ok(Self::Nanomole),
            _ => Err(Error::UnsupportedUnit),
        }
    }
    fn scale(self) -> Rational {
        Rational::integer(match self {
            Self::Picomole => 1,
            Self::Nanomole => 1000,
        })
    }
}
impl MolarUnit {
    pub fn parse(s: &str) -> Result<Self, Error> {
        match s {
            "nM" => Ok(Self::Nanomolar),
            "uM" | "µM" | "μM" => Ok(Self::Micromolar),
            _ => Err(Error::UnsupportedUnit),
        }
    }
    fn inverse_scale(self) -> u64 {
        match self {
            Self::Nanomolar => 1_000_000,
            Self::Micromolar => 1000,
        }
    }
}
macro_rules! quantity {
    ($name:ident,$tag:literal) => {
        #[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
        pub struct $name(Rational);
        impl $name {
            pub const ZERO: Self = Self(Rational::ZERO);
            pub const fn canonical(self) -> Rational {
                self.0
            }
            pub fn checked_add(self, rhs: Self) -> Result<Self, Error> {
                Ok(Self(self.0.checked_add(rhs.0)?))
            }
            pub fn checked_sub(self, rhs: Self) -> Result<Self, Error> {
                Ok(Self(self.0.checked_sub(rhs.0)?))
            }
            pub fn scale(self, factor: Rational) -> Result<Self, Error> {
                Ok(Self(self.0.checked_mul(factor)?))
            }
            pub const fn is_zero(self) -> bool {
                self.0.is_zero()
            }
        }
        #[cfg(feature = "serde")]
        impl serde::Serialize for $name {
            fn serialize<S: serde::Serializer>(&self, s: S) -> Result<S::Ok, S::Error> {
                use serde::ser::SerializeStruct;
                let mut state = s.serialize_struct(stringify!($name), 2)?;
                state.serialize_field("unit", $tag)?;
                state.serialize_field("value", &self.0)?;
                state.end()
            }
        }
        #[cfg(feature = "serde")]
        impl<'de> serde::Deserialize<'de> for $name {
            fn deserialize<D: serde::Deserializer<'de>>(d: D) -> Result<Self, D::Error> {
                #[derive(serde::Deserialize)]
                #[serde(deny_unknown_fields)]
                struct Wire {
                    unit: String,
                    value: Rational,
                }
                let w = Wire::deserialize(d)?;
                if w.unit != $tag {
                    return Err(serde::de::Error::custom(
                        "dimension or canonical unit mismatch",
                    ));
                }
                Ok(Self(w.value))
            }
        }
    };
}
quantity!(Volume, "nL");
quantity!(Mass, "ng");
quantity!(Amount, "pmol");
quantity!(MassConcentration, "ng/nL");
quantity!(MolarConcentration, "pmol/nL");
impl Volume {
    pub fn from_unit(value: Rational, unit: VolumeUnit) -> Result<Self, Error> {
        Ok(Self(value.checked_mul(unit.scale())?))
    }
    pub fn in_unit(self, unit: VolumeUnit) -> Result<Rational, Error> {
        self.0.checked_div(unit.scale())
    }
}
impl Mass {
    pub fn from_unit(value: Rational, unit: MassUnit) -> Result<Self, Error> {
        Ok(Self(value.checked_mul(unit.scale())?))
    }
    pub fn in_unit(self, unit: MassUnit) -> Result<Rational, Error> {
        self.0.checked_div(unit.scale())
    }
}
impl Amount {
    pub fn from_unit(value: Rational, unit: AmountUnit) -> Result<Self, Error> {
        Ok(Self(value.checked_mul(unit.scale())?))
    }
    pub fn in_unit(self, unit: AmountUnit) -> Result<Rational, Error> {
        self.0.checked_div(unit.scale())
    }
}
impl MassConcentration {
    pub fn from_units(value: Rational, mass: MassUnit, volume: VolumeUnit) -> Result<Self, Error> {
        Ok(Self(Rational::checked_product_ratio(
            &[value, mass.scale()],
            &[volume.scale()],
        )?))
    }
    pub fn in_units(self, mass: MassUnit, volume: VolumeUnit) -> Result<Rational, Error> {
        Rational::checked_product_ratio(&[self.0, volume.scale()], &[mass.scale()])
    }
    pub fn from_mass_volume(mass: Mass, volume: Volume) -> Result<Self, Error> {
        if volume.is_zero() {
            return Err(Error::ZeroFinalVolume);
        }
        Ok(Self(mass.0.checked_div(volume.0)?))
    }
    pub fn mass_in(self, volume: Volume) -> Result<Mass, Error> {
        Ok(Mass(self.0.checked_mul(volume.0)?))
    }
}
impl MolarConcentration {
    pub fn from_unit(value: Rational, unit: MolarUnit) -> Result<Self, Error> {
        Ok(Self(
            value.checked_div(Rational::integer(unit.inverse_scale()))?,
        ))
    }
    pub fn in_unit(self, unit: MolarUnit) -> Result<Rational, Error> {
        self.0.checked_mul(Rational::integer(unit.inverse_scale()))
    }
    pub fn from_amount_volume(amount: Amount, volume: Volume) -> Result<Self, Error> {
        if volume.is_zero() {
            return Err(Error::ZeroFinalVolume);
        }
        Ok(Self(amount.0.checked_div(volume.0)?))
    }
    pub fn amount_in(self, volume: Volume) -> Result<Amount, Error> {
        Ok(Amount(self.0.checked_mul(volume.0)?))
    }
}
mod sealed {
    pub trait Sealed {}
    impl Sealed for super::MassConcentration {}
    impl Sealed for super::MolarConcentration {}
}
/// Sealed to mass/volume and amount/volume concentrations. A dilution cannot mix dimensions.
pub trait Concentration: sealed::Sealed + Copy + Eq + core::fmt::Debug {
    #[doc(hidden)]
    fn exact(self) -> Rational;
    #[doc(hidden)]
    fn with_exact(value: Rational) -> Self;
}
impl Concentration for MassConcentration {
    fn exact(self) -> Rational {
        self.0
    }
    fn with_exact(value: Rational) -> Self {
        Self(value)
    }
}
impl Concentration for MolarConcentration {
    fn exact(self) -> Rational {
        self.0
    }
    fn with_exact(value: Rational) -> Self {
        Self(value)
    }
}
