use crate::{Error, Result};
use serde::{Deserialize, Serialize};
/// Neutral atom counts, not an ion/adduct or charge-state assignment.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct Formula {
    pub c: u64,
    pub h: u64,
    pub n: u64,
    pub o: u64,
    pub p: u64,
    pub s: u64,
    pub f: u64,
}
impl Formula {
    pub(crate) fn array(self) -> [u64; 7] {
        [self.c, self.h, self.n, self.o, self.p, self.s, self.f]
    }
    pub fn hill(self) -> String {
        let mut out = String::new();
        for (name, count) in [
            ("C", self.c),
            ("H", self.h),
            ("F", self.f),
            ("N", self.n),
            ("O", self.o),
            ("P", self.p),
            ("S", self.s),
        ] {
            if count > 0 {
                out.push_str(name);
                if count > 1 {
                    out.push_str(&count.to_string());
                }
            }
        }
        out
    }
}
pub(crate) const MONO: [f64; 7] = [
    12.0,
    1.00782503223,
    14.00307400443,
    15.99491461957,
    30.97376199842,
    31.9720711744,
    18.99840316273,
];
// Declared conventional averages, not universal exact natural-isotope masses.
pub(crate) const AVERAGE: [f64; 7] = [
    12.011,
    1.008,
    14.007,
    15.999,
    30.973761998,
    32.06,
    18.998403163,
];
pub const MASS_MODEL: &str = "neutral-acid-CHNOPS-F-NIST-isotopes-conventional-average-1";
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Mass {
    pub model: String,
    pub monoisotopic_da: f64,
    pub average_da: f64,
}
impl Mass {
    pub(crate) fn for_formula(f: Formula) -> Self {
        let sum = |weights: &[f64; 7]| {
            f.array()
                .into_iter()
                .zip(weights)
                .map(|(n, w)| n as f64 * w)
                .sum()
        };
        Self {
            model: MASS_MODEL.into(),
            monoisotopic_da: sum(&MONO),
            average_da: sum(&AVERAGE),
        }
    }
}
#[derive(Default)]
pub(crate) struct Sum([i64; 7]);
impl Sum {
    pub fn add(&mut self, delta: [i64; 7], count: usize) -> Result<()> {
        let count = i64::try_from(count).map_err(|_| Error::Limit("formula count"))?;
        for (value, d) in self.0.iter_mut().zip(delta) {
            *value = value
                .checked_add(
                    d.checked_mul(count)
                        .ok_or(Error::Limit("formula arithmetic"))?,
                )
                .ok_or(Error::Limit("formula arithmetic"))?;
        }
        Ok(())
    }
    pub fn finish(self) -> Result<Formula> {
        if self.0.iter().any(|&n| !(0..=100_000_000).contains(&n)) {
            return Err(Error::Invalid("invalid final atom counts"));
        }
        let [c, h, n, o, p, s, f] = self.0.map(|n| n as u64);
        Ok(Formula {
            c,
            h,
            n,
            o,
            p,
            s,
            f,
        })
    }
}
/// Net addition per 3'–5' phosphodiester: H3PO4 minus two waters.
pub(crate) const BOND: [i64; 7] = [0, -1, 0, 2, 1, 0, 0];
/// Terminal phosphomonoester: H3PO4 minus one water.
pub(crate) const PHOSPHATE: [i64; 7] = [0, 1, 0, 3, 1, 0, 0];
pub(crate) const TRIPHOSPHATE: [i64; 7] = [0, 3, 0, 9, 3, 0, 0];
pub(crate) const PHOSPHOROTHIOATE: [i64; 7] = [0, 0, 0, -1, 0, 1, 0];
pub(crate) const NUCLEOSIDES: [(u8, [i64; 7]); 4] = [
    (b'A', [10, 13, 5, 4, 0, 0, 0]),
    (b'C', [9, 13, 3, 5, 0, 0, 0]),
    (b'G', [10, 13, 5, 5, 0, 0, 0]),
    (b'U', [9, 12, 2, 6, 0, 0, 0]),
];
pub(crate) fn nucleoside(base: u8) -> Option<[i64; 7]> {
    NUCLEOSIDES
        .iter()
        .find(|(symbol, _)| *symbol == base)
        .map(|(_, atoms)| *atoms)
}
