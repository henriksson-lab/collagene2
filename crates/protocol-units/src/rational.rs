use crate::Error;
/// A canonical nonnegative rational with u64 numerator and nonzero u64 denominator.
/// Exact u128 intermediates are reduced before range checks. No floating point is used.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct Rational {
    numerator: u64,
    denominator: u64,
}
fn gcd(mut a: u128, mut b: u128) -> u128 {
    while b != 0 {
        let r = a % b;
        a = b;
        b = r;
    }
    a
}
impl Rational {
    pub const ZERO: Self = Self {
        numerator: 0,
        denominator: 1,
    };
    pub const ONE: Self = Self {
        numerator: 1,
        denominator: 1,
    };
    pub const fn integer(n: u64) -> Self {
        Self {
            numerator: n,
            denominator: 1,
        }
    }
    pub fn new(numerator: u64, denominator: u64) -> Result<Self, Error> {
        Self::wide(numerator.into(), denominator.into())
    }
    pub(crate) fn wide(n: u128, d: u128) -> Result<Self, Error> {
        if d == 0 {
            return Err(Error::DivisionByZero);
        }
        if n == 0 {
            return Ok(Self::ZERO);
        }
        let g = gcd(n, d);
        Ok(Self {
            numerator: (n / g).try_into().map_err(|_| Error::Overflow)?,
            denominator: (d / g).try_into().map_err(|_| Error::Overflow)?,
        })
    }
    pub const fn numerator(self) -> u64 {
        self.numerator
    }
    pub const fn denominator(self) -> u64 {
        self.denominator
    }
    pub const fn is_zero(self) -> bool {
        self.numerator == 0
    }
    pub fn checked_add(self, rhs: Self) -> Result<Self, Error> {
        let g = gcd(self.denominator.into(), rhs.denominator.into());
        let a = u128::from(self.numerator) * (u128::from(rhs.denominator) / g);
        let b = u128::from(rhs.numerator) * (u128::from(self.denominator) / g);
        Self::wide(
            a.checked_add(b).ok_or(Error::Overflow)?,
            u128::from(self.denominator) / g * u128::from(rhs.denominator),
        )
    }
    pub fn checked_sub(self, rhs: Self) -> Result<Self, Error> {
        let g = gcd(self.denominator.into(), rhs.denominator.into());
        let a = u128::from(self.numerator) * (u128::from(rhs.denominator) / g);
        let b = u128::from(rhs.numerator) * (u128::from(self.denominator) / g);
        Self::wide(
            a.checked_sub(b).ok_or(Error::NegativeQuantity)?,
            u128::from(self.denominator) / g * u128::from(rhs.denominator),
        )
    }
    pub fn checked_mul(self, rhs: Self) -> Result<Self, Error> {
        Self::wide(
            u128::from(self.numerator) * u128::from(rhs.numerator),
            u128::from(self.denominator) * u128::from(rhs.denominator),
        )
    }
    pub fn checked_div(self, rhs: Self) -> Result<Self, Error> {
        Self::wide(
            u128::from(self.numerator) * u128::from(rhs.denominator),
            u128::from(self.denominator) * u128::from(rhs.numerator),
        )
    }
    /// Compute a bounded product/quotient with cross-cancellation across all factors.
    /// This avoids an overflowing intermediate when the final reduced rational fits.
    pub fn checked_product_ratio(
        numerators: &[Self],
        denominators: &[Self],
    ) -> Result<Self, Error> {
        if numerators.len().saturating_add(denominators.len()) > 16 {
            return Err(Error::TooManyFactors);
        }
        if denominators.iter().any(|x| x.is_zero()) {
            return Err(Error::DivisionByZero);
        }
        if numerators.iter().any(|x| x.is_zero()) {
            return Ok(Self::ZERO);
        }
        let mut ns: Vec<u64> = numerators
            .iter()
            .map(|x| x.numerator)
            .chain(denominators.iter().map(|x| x.denominator))
            .collect();
        let mut ds: Vec<u64> = numerators
            .iter()
            .map(|x| x.denominator)
            .chain(denominators.iter().map(|x| x.numerator))
            .collect();
        for n in &mut ns {
            for d in &mut ds {
                let g = gcd((*n).into(), (*d).into()) as u64;
                *n /= g;
                *d /= g;
            }
        }
        let n = ns
            .into_iter()
            .try_fold(1u64, |a, b| a.checked_mul(b).ok_or(Error::Overflow))?;
        let d = ds
            .into_iter()
            .try_fold(1u64, |a, b| a.checked_mul(b).ok_or(Error::Overflow))?;
        Self::new(n, d)
    }
    /// Sum at most 128 terms using reduced checked u128 intermediates.
    /// The final rational must fit u64/u64; excessive intermediate denominators return Overflow.
    pub fn checked_sum(values: &[Self]) -> Result<Self, Error> {
        if values.len() > crate::MAX_COMPONENTS {
            return Err(Error::TooManyComponents);
        }
        let (mut n, mut d) = (0u128, 1u128);
        for value in values {
            let g = gcd(d, u128::from(value.denominator));
            let left = n
                .checked_mul(u128::from(value.denominator) / g)
                .ok_or(Error::Overflow)?;
            let right = u128::from(value.numerator)
                .checked_mul(d / g)
                .ok_or(Error::Overflow)?;
            let next_n = left.checked_add(right).ok_or(Error::Overflow)?;
            let next_d = (d / g)
                .checked_mul(u128::from(value.denominator))
                .ok_or(Error::Overflow)?;
            let reduce = gcd(next_n, next_d);
            n = next_n / reduce;
            d = next_d / reduce;
        }
        Self::wide(n, d)
    }
    /// Exact finite decimal input, at most 64 bytes / 19 fractional digits. No exponent notation.
    pub fn parse_decimal(value: &str) -> Result<Self, Error> {
        if value.len() > 64 {
            return Err(Error::InvalidNumber);
        }
        let value = value.trim();
        if value.starts_with('-') {
            return Err(Error::NegativeQuantity);
        }
        let mut parts = value.split('.');
        let whole = parts.next().ok_or(Error::InvalidNumber)?;
        let frac = parts.next().unwrap_or("");
        if parts.next().is_some()
            || whole.is_empty()
            || !whole.bytes().all(|b| b.is_ascii_digit())
            || !frac.bytes().all(|b| b.is_ascii_digit())
        {
            return Err(Error::InvalidNumber);
        }
        if frac.len() > 19 {
            return Err(Error::DecimalPrecisionLimit);
        }
        let n = format!("{whole}{frac}")
            .parse::<u128>()
            .map_err(|_| Error::Overflow)?;
        let d = 10u128
            .checked_pow(frac.len() as u32)
            .ok_or(Error::Overflow)?;
        Self::wide(n, d)
    }
}
impl Ord for Rational {
    fn cmp(&self, rhs: &Self) -> core::cmp::Ordering {
        (u128::from(self.numerator) * u128::from(rhs.denominator))
            .cmp(&(u128::from(rhs.numerator) * u128::from(self.denominator)))
    }
}
impl PartialOrd for Rational {
    fn partial_cmp(&self, rhs: &Self) -> Option<core::cmp::Ordering> {
        Some(self.cmp(rhs))
    }
}
impl core::fmt::Display for Rational {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        if self.denominator == 1 {
            self.numerator.fmt(f)
        } else {
            write!(f, "{}/{}", self.numerator, self.denominator)
        }
    }
}
#[cfg(feature = "serde")]
impl serde::Serialize for Rational {
    fn serialize<S: serde::Serializer>(&self, s: S) -> Result<S::Ok, S::Error> {
        use serde::ser::SerializeStruct;
        let mut state = s.serialize_struct("Rational", 2)?;
        state.serialize_field("numerator", &self.numerator.to_string())?;
        state.serialize_field("denominator", &self.denominator.to_string())?;
        state.end()
    }
}
#[cfg(feature = "serde")]
impl<'de> serde::Deserialize<'de> for Rational {
    fn deserialize<D: serde::Deserializer<'de>>(d: D) -> Result<Self, D::Error> {
        #[derive(serde::Deserialize)]
        #[serde(deny_unknown_fields)]
        struct Wire {
            numerator: String,
            denominator: String,
        }
        let w = Wire::deserialize(d)?;
        fn parse(s: &str) -> Option<u64> {
            if s.is_empty() || s.len() > 20 || !s.bytes().all(|b| b.is_ascii_digit()) {
                None
            } else {
                s.parse().ok()
            }
        }
        let n = parse(&w.numerator)
            .ok_or_else(|| serde::de::Error::custom("invalid rational numerator"))?;
        let den = parse(&w.denominator)
            .ok_or_else(|| serde::de::Error::custom("invalid rational denominator"))?;
        Self::new(n, den).map_err(serde::de::Error::custom)
    }
}
