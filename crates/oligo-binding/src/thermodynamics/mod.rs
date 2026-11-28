//! Isolated internal DNA mismatch parameter replacements, not a second Tm solver.
//!
//! Top bases run 5′→3′; bottom bases run 3′→5′. The caller replaces the two
//! canonical nearest-neighbor terms, retaining its established terminal terms.
//! A mismatched duplex is a heteroduplex: the caller must remove any hypothetical
//! perfect-duplex symmetry correction and use the declared equal-strand convention.
mod parameters;

use crate::{Hit, Oligo, Options, invalid};
use collagene_model::{Error, Sequence, Topology};
use serde::{Deserialize, Serialize};

pub const ENGINE: &str = "isolated-internal-dna-mismatch-replacements-1";
pub const SOURCE_SHA256: &str = "4f3cbc585ffd68d54386d40f5edb113a7c2eae4aa4388c6b496bceecdec06500";
pub const LICENSE_SHA256: &str = "d4bbb0c4e5345e6f2508c3b93fb247fcb2db95783fb8a584c23cd88e3c03df6c";
pub const SOURCE_URL: &str = "https://raw.githubusercontent.com/biopython/biopython/biopython-186/Bio/SeqUtils/MeltingTemp.py";
pub const SCOPE: &str = "concrete DNA; one internal mismatch; two canonical pairs at each end; distinct physical strands; pH7; sodium only; equal concentrations; SantaLucia1998 canonical baseline";
pub const PRIMARY_DOIS: [&str; 5] = [
    "10.1021/bi962590c",
    "10.1021/bi9803729",
    "10.1021/bi9724873",
    "10.1093/nar/26.11.2694",
    "10.1021/bi9825091",
];

#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Profile {
    #[serde(deserialize_with = "crate::decode::identifier")]
    pub engine: String,
    #[serde(deserialize_with = "crate::decode::identifier")]
    pub revision: String,
    #[serde(deserialize_with = "crate::decode::name")]
    pub scope: String,
    #[serde(deserialize_with = "crate::decode::name")]
    pub numerical_source_url: String,
    #[serde(deserialize_with = "crate::decode::identifier")]
    pub numerical_source_sha256: String,
    #[serde(deserialize_with = "crate::decode::identifier")]
    pub license_sha256: String,
    pub parameter_count: usize,
}
impl Profile {
    pub fn current() -> Self {
        Self {
            engine: ENGINE.into(),
            revision: profile_revision(),
            scope: SCOPE.into(),
            numerical_source_url: SOURCE_URL.into(),
            numerical_source_sha256: SOURCE_SHA256.into(),
            license_sha256: LICENSE_SHA256.into(),
            parameter_count: parameters::PARAMETERS.len(),
        }
    }
    pub fn validate(&self) -> Result<(), Error> {
        if self != &Self::current() {
            return Err(invalid("unsupported or altered mismatch parameter profile"));
        }
        Ok(())
    }
}

#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub enum ConcentrationConvention {
    /// Two distinct strands at equal concentrations; C_total/4 in the log term.
    EqualConcentrationTotalStrands,
}
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub enum SaltConvention {
    /// Entropy correction 0.368*(n-1)*ln[Na+]; no Mg/dNTP equivalence.
    SodiumEntropy1998,
}
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub enum PhConvention {
    /// The empirical parameter convention, not a measured buffer-pH assertion.
    Ph7,
}
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Assumptions {
    pub sodium_molar: f64,
    pub total_strand_molar: f64,
    pub concentration: ConcentrationConvention,
    pub salt: SaltConvention,
    pub ph: PhConvention,
}
impl Assumptions {
    pub fn validate(&self) -> Result<(), Error> {
        if !self.sodium_molar.is_finite()
            || !(0.01..=1.0).contains(&self.sodium_molar)
            || !self.total_strand_molar.is_finite()
            || !(1e-12..=0.01).contains(&self.total_strand_molar)
        {
            return Err(invalid(
                "mismatch model sodium/concentration outside declared bounds",
            ));
        }
        Ok(())
    }
}

/// A physical duplex step: top 5′XY3′, bottom 3′ZW5′.
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Stack {
    pub top: [char; 2],
    pub bottom: [char; 2],
}
impl Stack {
    /// Rotation by 180°, preserving strand chemical directions (not top-only RC).
    pub fn rotated(self) -> Self {
        Self {
            top: [self.bottom[1], self.bottom[0]],
            bottom: [self.top[1], self.top[0]],
        }
    }
    fn key(self) -> String {
        format!(
            "{}{}/{}{}",
            self.top[0], self.top[1], self.bottom[0], self.bottom[1]
        )
    }
}

#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Replacement {
    /// Index into this footprint, not the original source or ordered oligo.
    pub column_start: usize,
    pub actual: Stack,
    /// The canonical step to subtract using the caller's existing NN owner.
    pub canonical: Stack,
    /// Absolute mismatch-step parameter, not a penalty added to a canonical term.
    pub delta_h_tenths_kcal: i16,
    pub delta_s_tenths_cal_per_kelvin: i16,
}
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ReplacementPlan {
    pub mismatch_column: usize,
    /// The two affected steps, each replaced exactly once.
    pub replacements: [Replacement; 2],
    /// True only after excluding identical physical strand sequences. A
    /// self-complementary top can still have a distinct mismatched partner.
    pub heteroduplex: bool,
}
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub enum Unsupported {
    AmbiguousBases,
    TooShort,
    MultipleMismatches,
    /// Identical physical strands require a separate symmetry/concentration model.
    SelfAssociation,
    /// Conservative scope policy: two canonical pairs required at either end.
    MismatchNearEnd,
    MissingParameter,
}
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub enum Outcome {
    /// Route to the existing perfect-duplex owner, including its symmetry rules.
    PerfectDuplex,
    Replacements(ReplacementPlan),
    Unsupported(Unsupported),
}

#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Assessment {
    #[serde(deserialize_with = "crate::decode::identifier")]
    pub engine: String,
    #[serde(deserialize_with = "crate::decode::identifier")]
    pub profile_revision: String,
    #[serde(deserialize_with = "crate::decode::identifier")]
    pub input_revision: String,
    pub assumptions: Assumptions,
    #[serde(deserialize_with = "crate::decode::name")]
    pub top_5_to_3: String,
    #[serde(deserialize_with = "crate::decode::name")]
    pub bottom_3_to_5: String,
    /// Original plus-reference positions in annealing-column order. Reverse
    /// footprints descend and circular footprints retain their actual wrap.
    #[serde(deserialize_with = "collagene_model::decode::bounded_vec::<_, usize, 200>")]
    pub source_positions: Vec<usize>,
    pub oligo_start: usize,
    pub outcome: Outcome,
}
impl Assessment {
    /// Documentary integrity only. It does not authenticate source identity,
    /// coordinate authority or a complete/optimal search claim.
    pub fn validate(&self) -> Result<(), Error> {
        let length = self.top_5_to_3.len();
        if length == 0
            || length > crate::MAX_OLIGO_BASES
            || self.bottom_3_to_5.len() != length
            || self.source_positions.len() != length
            || self.oligo_start > crate::MAX_OLIGO_BASES - length
        {
            return Err(invalid("invalid mismatch footprint bounds"));
        }
        self.assumptions.validate()?;
        let dna = |value: &str| value.bytes().all(|b| b"ACGTRYSWKMBDHVN".contains(&b));
        if self.engine != ENGINE
            || self.profile_revision != profile_revision()
            || self.input_revision.len() != 64
            || !self
                .input_revision
                .bytes()
                .all(|b| b.is_ascii_digit() || (b'a'..=b'f').contains(&b))
            || !dna(&self.top_5_to_3)
            || !dna(&self.bottom_3_to_5)
            || self
                .source_positions
                .iter()
                .any(|&p| p >= collagene_model::MAX_SEQUENCE_LEN)
        {
            return Err(invalid(
                "unsupported mismatch profile or invalid footprint facts",
            ));
        }
        let unique: std::collections::BTreeSet<_> = self.source_positions.iter().collect();
        let ordered = |ascending: bool| {
            let mut wraps = 0;
            self.source_positions.windows(2).all(|w| {
                let (a, b) = if ascending {
                    (w[0], w[1])
                } else {
                    (w[1], w[0])
                };
                if a + 1 == b {
                    true
                } else if a > b {
                    wraps += 1;
                    wraps <= 1
                } else {
                    false
                }
            })
        };
        if unique.len() != length
            || (!ordered(true) && !ordered(false))
            || self.outcome != evaluate(self.top_5_to_3.as_bytes(), self.bottom_3_to_5.as_bytes())
        {
            return Err(invalid("altered mismatch ordering or parameter evidence"));
        }
        Ok(())
    }

    /// Recompute the bounded parameter facts before accepting retained evidence.
    pub fn validate_for(
        &self,
        context: &DuplexContext<'_>,
        oligo: &Oligo,
        options: &Options,
        hit: &Hit,
    ) -> Result<(), Error> {
        self.validate()?;
        if self != &context.assess(oligo, options, hit, self.assumptions)? {
            return Err(invalid("mismatch parameter evidence is stale or altered"));
        }
        Ok(())
    }
}

/// Identity includes exact source captures, the selected numerical table, and
/// this implementation's deliberately narrower admission policy.
pub fn profile_revision() -> String {
    let mut h = blake3::Hasher::new();
    h.update(ENGINE.as_bytes());
    h.update(SCOPE.as_bytes());
    h.update(SOURCE_SHA256.as_bytes());
    h.update(LICENSE_SHA256.as_bytes());
    for (key, dh, ds) in parameters::PARAMETERS {
        h.update(key.as_bytes());
        h.update(&dh.to_le_bytes());
        h.update(&ds.to_le_bytes());
    }
    h.finalize().to_hex().to_string()
}

/// Borrowed collection-assessment context; hashes the source once for all hits.
/// It owns no second copy of a potentially large source sequence.
pub struct DuplexContext<'a> {
    source: &'a Sequence,
    topology: Topology,
    source_revision: String,
}
impl<'a> DuplexContext<'a> {
    pub fn new(source: &'a Sequence, topology: Topology) -> Result<Self, Error> {
        crate::check_source(source)?;
        Ok(Self {
            source,
            topology,
            source_revision: crate::identity::source(source, topology),
        })
    }

    /// Validate the original hit and inspect its bounded annealing footprint.
    /// Declared 5′ tails and extra unpaired prefixes never enter these terms.
    pub fn assess(
        &self,
        oligo: &Oligo,
        options: &Options,
        hit: &Hit,
        assumptions: Assumptions,
    ) -> Result<Assessment, Error> {
        assumptions.validate()?;
        hit.validate_columns(self.source, self.topology, oligo, options)?;
        let top: String = hit.columns.iter().map(|c| c.oligo_base).collect();
        let bottom: String = hit.columns.iter().map(|c| c.template_base).collect();
        let outcome = evaluate(top.as_bytes(), bottom.as_bytes());
        let profile = profile_revision();
        let options_revision = options.content_hash()?;
        let mut hash = blake3::Hasher::new();
        for value in [
            ENGINE,
            &profile,
            &self.source_revision,
            &options_revision,
            &oligo.id,
            &oligo.name,
            oligo.ordered.as_str(),
            &top,
            &bottom,
        ] {
            hash.update(&(value.len() as u64).to_le_bytes());
            hash.update(value.as_bytes());
        }
        for n in [
            oligo.declared_tail_length,
            hit.annealing_start,
            hit.oligo_index,
        ] {
            hash.update(&(n as u64).to_le_bytes());
        }
        hash.update(&[match hit.orientation {
            crate::Orientation::Forward => 0,
            crate::Orientation::Reverse => 1,
        }]);
        for c in &hit.columns {
            hash.update(&(c.source_position as u64).to_le_bytes());
        }
        hash.update(&assumptions.sodium_molar.to_bits().to_le_bytes());
        hash.update(&assumptions.total_strand_molar.to_bits().to_le_bytes());
        Ok(Assessment {
            engine: ENGINE.to_owned(),
            profile_revision: profile,
            input_revision: hash.finalize().to_hex().to_string(),
            assumptions,
            top_5_to_3: top,
            bottom_3_to_5: bottom,
            source_positions: hit.columns.iter().map(|c| c.source_position).collect(),
            oligo_start: hit.annealing_start,
            outcome,
        })
    }
}

fn complement(base: u8) -> Option<u8> {
    match base {
        b'A' => Some(b'T'),
        b'C' => Some(b'G'),
        b'G' => Some(b'C'),
        b'T' => Some(b'A'),
        _ => None,
    }
}
fn parameter(stack: Stack) -> Option<(i16, i16)> {
    let key = stack.key();
    let rotated = stack.rotated().key();
    parameters::PARAMETERS
        .iter()
        .find_map(|&(k, h, s)| (k == key || k == rotated).then_some((h, s)))
}
fn evaluate(top: &[u8], bottom: &[u8]) -> Outcome {
    if top.iter().chain(bottom).any(|&b| complement(b).is_none()) {
        return Outcome::Unsupported(Unsupported::AmbiguousBases);
    }
    if top.len() < 2 {
        return Outcome::Unsupported(Unsupported::TooShort);
    }
    let mismatches: Vec<_> = top
        .iter()
        .zip(bottom)
        .enumerate()
        .filter_map(|(i, (&a, &b))| (complement(a) != Some(b)).then_some(i))
        .collect();
    if mismatches.is_empty() {
        return Outcome::PerfectDuplex;
    }
    if mismatches.len() != 1 {
        return Outcome::Unsupported(Unsupported::MultipleMismatches);
    }
    if top.iter().copied().eq(bottom.iter().rev().copied()) {
        return Outcome::Unsupported(Unsupported::SelfAssociation);
    }
    let mismatch = mismatches[0];
    if mismatch < 2 || mismatch + 2 >= top.len() {
        return Outcome::Unsupported(Unsupported::MismatchNearEnd);
    }
    let replacement = |start: usize| -> Option<Replacement> {
        let actual = Stack {
            top: [top[start] as char, top[start + 1] as char],
            bottom: [bottom[start] as char, bottom[start + 1] as char],
        };
        let canonical = Stack {
            top: actual.top,
            bottom: [
                complement(top[start])? as char,
                complement(top[start + 1])? as char,
            ],
        };
        let (h, s) = parameter(actual)?;
        Some(Replacement {
            column_start: start,
            actual,
            canonical,
            delta_h_tenths_kcal: h,
            delta_s_tenths_cal_per_kelvin: s,
        })
    };
    match (replacement(mismatch - 1), replacement(mismatch)) {
        (Some(left), Some(right)) => Outcome::Replacements(ReplacementPlan {
            mismatch_column: mismatch,
            replacements: [left, right],
            heteroduplex: true,
        }),
        _ => Outcome::Unsupported(Unsupported::MissingParameter),
    }
}

#[cfg(test)]
mod tests;
