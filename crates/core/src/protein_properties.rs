//! Explicit unmodified-chain protein calculations, independently implemented.
//! These are sequence/model estimates, not folding, activity, or experimental measurements.
use crate::{Alphabet, Error, Molecule};
use serde::{Deserialize, Serialize};

pub const PROTEIN_ENGINE_VERSION: &str = "protein-properties-1";
pub const MASS_MODEL: &str = "UWPR-residue-elements-2026-11-24";
pub const CHARGE_MODEL: &str = "independent-Henderson-Hasselbalch-Epk-6.6";
pub const EXTINCTION_MODEL: &str = "Pace-1995-280nm-water";
const SYMBOLS: &[u8] = b"ACDEFGHIKLMNPQRSTVWYBJXZUO";
const STANDARD: &[u8] = b"ACDEFGHIKLMNPQRSTVWY";
const MAX_SOURCE_BYTES: usize = 32 * 1024 * 1024;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
pub enum StopPolicy {
    Reject,
    AllowSingleTerminal,
}
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
pub enum AmbiguityPolicy {
    Unavailable,
    CanonicalBounds,
}
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
pub enum DisulfideState {
    Reduced,
    MaximumIntrachain,
    Count(usize),
}
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
pub enum ModificationState {
    AssumeUnmodified,
    KnownOrUnknownModifications,
}
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct ProteinOptions {
    pub stop_policy: StopPolicy,
    pub ambiguity: AmbiguityPolicy,
    pub disulfides: DisulfideState,
    pub modification_state: ModificationState,
    pub charge_ph: f64,
    pub max_residues: usize,
}
impl Default for ProteinOptions {
    fn default() -> Self {
        Self {
            stop_policy: StopPolicy::AllowSingleTerminal,
            ambiguity: AmbiguityPolicy::Unavailable,
            disulfides: DisulfideState::Reduced,
            modification_state: ModificationState::AssumeUnmodified,
            charge_ph: 7.0,
            max_residues: 1_000_000,
        }
    }
}
impl ProteinOptions {
    pub fn validate(&self) -> Result<(), Error> {
        if !self.charge_ph.is_finite()
            || !(0.0..=14.0).contains(&self.charge_ph)
            || self.max_residues == 0
            || self.max_residues > crate::model::MAX_SEQUENCE_LEN
        {
            return Err(Error::Unsupported(
                "protein pH must be finite in0..14 and residue limit in1..10,000,000".into(),
            ));
        }
        Ok(())
    }
}
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct NumericInterval {
    pub minimum: f64,
    pub maximum: f64,
}
impl NumericInterval {
    fn new(minimum: f64, maximum: f64) -> Self {
        Self { minimum, maximum }
    }
    fn point(value: f64) -> Self {
        Self::new(value, value)
    }
    fn valid(&self) -> bool {
        self.minimum.is_finite() && self.maximum.is_finite() && self.minimum <= self.maximum
    }
}
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct ResidueCount {
    pub symbol: String,
    pub count: usize,
    pub fraction: f64,
}
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct MassEstimate {
    pub model: String,
    pub average_da: NumericInterval,
    pub monoisotopic_da: NumericInterval,
    /// False for unresolved residue identity, even if isobaric J gives point mass.
    pub exact_identity: bool,
}
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct ExtinctionEstimate {
    pub model: String,
    pub wavelength_nm: u16,
    pub epsilon_molar_inverse_cm: NumericInterval,
    /// Beer–Lambert A280 for1mg/mL,1cm path. Uses average chain mass.
    pub absorbance_1mg_ml_1cm: Option<NumericInterval>,
}
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct ProteinReport {
    pub source_revision: String,
    pub source_sequence_blake3: String,
    pub engine_version: String,
    pub options: ProteinOptions,
    pub source_length: usize,
    pub residue_count: usize,
    pub terminal_stop_removed: bool,
    pub composition: Vec<ResidueCount>,
    pub disulfide_bonds: NumericInterval,
    pub mass: Option<MassEstimate>,
    pub charge_model: String,
    pub isoelectric_point: Option<f64>,
    pub net_charge_at_ph: Option<f64>,
    pub extinction_280: Option<ExtinctionEstimate>,
    pub unavailable: Vec<String>,
    pub warnings: Vec<String>,
    /// All requested quantities have resolved sequence identities under the declared model.
    /// It never means experimentally verified or complete biological semantics.
    pub complete: bool,
}

// Residue elemental counts C,H,N,O,S. H2O is added once for free linear termini.
fn atoms(symbol: u8) -> Option<[u32; 5]> {
    Some(match symbol {
        b'A' => [3, 5, 1, 1, 0],
        b'C' => [3, 5, 1, 1, 1],
        b'D' => [4, 5, 1, 3, 0],
        b'E' => [5, 7, 1, 3, 0],
        b'F' => [9, 9, 1, 1, 0],
        b'G' => [2, 3, 1, 1, 0],
        b'H' => [6, 7, 3, 1, 0],
        b'I' | b'L' | b'J' => [6, 11, 1, 1, 0],
        b'K' => [6, 12, 2, 1, 0],
        b'M' => [5, 9, 1, 1, 1],
        b'N' => [4, 6, 2, 2, 0],
        b'P' => [5, 7, 1, 1, 0],
        b'Q' => [5, 8, 2, 2, 0],
        b'R' => [6, 12, 4, 1, 0],
        b'S' => [3, 5, 1, 2, 0],
        b'T' => [4, 7, 1, 2, 0],
        b'V' => [5, 9, 1, 1, 0],
        b'W' => [11, 10, 2, 1, 0],
        b'Y' => [9, 9, 1, 2, 0],
        _ => return None,
    })
}
const AVERAGE: [f64; 5] = [12.0107, 1.00794, 14.0067, 15.9994, 32.065];
const MONO: [f64; 5] = [12.0, 1.007825035, 14.003074, 15.99491463, 31.9720707];
fn residue_mass(symbol: u8, weights: &[f64; 5]) -> Option<f64> {
    let counts = atoms(symbol)?;
    Some(
        counts
            .iter()
            .zip(weights)
            .map(|(n, w)| f64::from(*n) * w)
            .sum(),
    )
}
fn possibilities(symbol: u8) -> &'static [u8] {
    match symbol {
        b'B' => b"DN",
        b'Z' => b"EQ",
        b'J' => b"IL",
        b'X' => STANDARD,
        _ => &[],
    }
}
fn mass_bounds(
    counts: &[usize; 26],
    weights: &[f64; 5],
    bonds: &NumericInterval,
) -> NumericInterval {
    let water = 2.0 * weights[1] + weights[3];
    let mut low = water;
    let mut high = water;
    for (symbol, count) in SYMBOLS.iter().copied().zip(counts.iter().copied()) {
        if count == 0 {
            continue;
        }
        let (minimum, maximum) = if let Some(value) = residue_mass(symbol, weights) {
            (value, value)
        } else {
            let values = possibilities(symbol)
                .iter()
                .filter_map(|s| residue_mass(*s, weights));
            values.fold((f64::INFINITY, f64::NEG_INFINITY), |(lo, hi), value| {
                (lo.min(value), hi.max(value))
            })
        };
        low += count as f64 * minimum;
        high += count as f64 * maximum;
    }
    NumericInterval::new(
        low - 2.0 * weights[1] * bonds.maximum,
        high - 2.0 * weights[1] * bonds.minimum,
    )
}
fn count(counts: &[usize; 26], symbol: u8) -> usize {
    SYMBOLS
        .iter()
        .position(|s| *s == symbol)
        .map(|i| counts[i])
        .unwrap_or(0)
}
// Independent group titration: constant published Epk pKas, no residue-context corrections.
fn charge(counts: &[usize; 26], ph: f64, bonds: usize) -> f64 {
    let positive = |pka: f64| 1.0 / (1.0 + 10f64.powf(ph - pka));
    let negative = |pka: f64| 1.0 / (1.0 + 10f64.powf(pka - ph));
    positive(8.6) - negative(3.6)
        + count(counts, b'H') as f64 * positive(6.5)
        + count(counts, b'K') as f64 * positive(10.8)
        + count(counts, b'R') as f64 * positive(12.5)
        - count(counts, b'D') as f64 * negative(3.9)
        - count(counts, b'E') as f64 * negative(4.1)
        - (count(counts, b'C') - 2 * bonds) as f64 * negative(8.5)
        - count(counts, b'Y') as f64 * negative(10.1)
}
fn pi(counts: &[usize; 26], bonds: usize) -> Option<f64> {
    if charge(counts, 0., bonds) <= 0. || charge(counts, 14., bonds) >= 0. {
        return None;
    }
    let (mut lo, mut hi) = (0., 14.);
    for _ in 0..64 {
        let mid = (lo + hi) / 2.;
        if charge(counts, mid, bonds) > 0. {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    Some((lo + hi) / 2.)
}

fn source_preflight(m: &Molecule) -> Result<(), Error> {
    fn reserve(bytes: &mut usize, items: &mut usize, size: usize) -> Result<(), Error> {
        *bytes = bytes.saturating_add(size).saturating_add(64);
        *items = items.saturating_add(1);
        if *bytes > MAX_SOURCE_BYTES || *items > 1_000_000 {
            return Err(Error::TooLarge);
        }
        Ok(())
    }
    fn location(
        l: &crate::Location,
        depth: usize,
        bytes: &mut usize,
        items: &mut usize,
    ) -> Result<(), Error> {
        if depth > 64 {
            return Err(Error::TooLarge);
        }
        reserve(bytes, items, 96)?;
        match l {
            crate::Location::Complement(inner) => location(inner, depth + 1, bytes, items)?,
            crate::Location::Join(parts) | crate::Location::Order(parts) => {
                for part in parts {
                    location(part, depth + 1, bytes, items)?;
                }
            }
            crate::Location::Unresolved(text) => reserve(bytes, items, text.len())?,
            crate::Location::Range { .. } => {}
        }
        Ok(())
    }
    let mut bytes = 0usize;
    let mut items = 0usize;
    reserve(
        &mut bytes,
        &mut items,
        m.name.len().saturating_add(m.sequence.len()),
    )?;
    for (k, v) in &m.metadata {
        reserve(&mut bytes, &mut items, k.len().saturating_add(v.len()))?;
    }
    if let Some(text) = &m.original_genbank {
        reserve(&mut bytes, &mut items, text.len())?;
    }
    for f in &m.features {
        reserve(&mut bytes, &mut items, f.kind.len())?;
        location(&f.location, 0, &mut bytes, &mut items)?;
        for (k, v) in &f.qualifiers {
            reserve(
                &mut bytes,
                &mut items,
                k.len().saturating_add(v.as_ref().map_or(0, String::len)),
            )?;
        }
    }
    Ok(())
}
fn has_modification_evidence(m: &Molecule) -> bool {
    m.metadata.keys().any(|key| {
        key.eq_ignore_ascii_case("protein:modifications")
            || key.eq_ignore_ascii_case("protein:modified")
            || key.eq_ignore_ascii_case("modifications")
    }) || m.features.iter().any(|f| {
        [
            "MOD_RES",
            "modified_residue",
            "cross-link",
            "cross_link",
            "bond",
            "LIPID",
            "CARBOHYD",
            "DISULFID",
        ]
        .iter()
        .any(|kind| f.kind.eq_ignore_ascii_case(kind))
    })
}

pub fn analyze_protein(
    molecule: &Molecule,
    options: &ProteinOptions,
) -> Result<ProteinReport, Error> {
    options.validate()?;
    source_preflight(molecule)?;
    molecule.validate()?;
    if molecule.sequence.alphabet() != Alphabet::Protein {
        return Err(Error::Unsupported("protein properties require an explicit Protein document; translate nucleotide sources first".into()));
    }
    let original = molecule.sequence.as_str().as_bytes();
    let terminal_stop_removed =
        original.last() == Some(&b'*') && options.stop_policy == StopPolicy::AllowSingleTerminal;
    let residues = if terminal_stop_removed {
        &original[..original.len() - 1]
    } else {
        original
    };
    if residues.is_empty() {
        return Err(Error::EmptySequence);
    }
    if residues.len() > options.max_residues {
        return Err(Error::TooLarge);
    }
    if residues.contains(&b'*') {
        return Err(Error::Unsupported("internal/repeated stop symbols do not describe one covalent protein chain; select a chain explicitly".into()));
    }
    let mut counts = [0usize; 26];
    for &symbol in residues {
        let Some(i) = SYMBOLS.iter().position(|s| *s == symbol) else {
            return Err(Error::Unsupported("unsupported protein residue".into()));
        };
        counts[i] += 1;
    }
    let cys = count(&counts, b'C');
    let x = count(&counts, b'X');
    let ambiguous = count(&counts, b'B') + count(&counts, b'Z') + x;
    let unresolved = ambiguous + count(&counts, b'J');
    let nonstandard = count(&counts, b'U') + count(&counts, b'O');
    let bonds = match options.disulfides {
        DisulfideState::Reduced => NumericInterval::point(0.),
        DisulfideState::MaximumIntrachain => {
            NumericInterval::new((cys / 2) as f64, ((cys + x) / 2) as f64)
        }
        DisulfideState::Count(n) => {
            if n > cys / 2 {
                return Err(Error::Unsupported("declared intrachain disulfides require two explicit C residues per bond; unknown residues cannot authorize bonds".into()));
            }
            NumericInterval::point(n as f64)
        }
    };
    let mut report=ProteinReport{source_revision:molecule.revision(),source_sequence_blake3:blake3::hash(original).to_hex().to_string(),engine_version:PROTEIN_ENGINE_VERSION.into(),options:options.clone(),source_length:original.len(),residue_count:residues.len(),terminal_stop_removed,composition:SYMBOLS.iter().copied().zip(counts).filter(|(_,n)|*n>0).map(|(symbol,n)|ResidueCount{symbol:(symbol as char).to_string(),count:n,fraction:n as f64/residues.len()as f64}).collect(),disulfide_bonds:bonds.clone(),mass:None,charge_model:CHARGE_MODEL.into(),isoelectric_point:None,net_charge_at_ph:None,extinction_280:None,unavailable:vec![],warnings:vec!["Single linear chain with free, unmodified termini; no folding, processing, multimerization, other chromophores or empirical activity is inferred.".into()],complete:false};
    if terminal_stop_removed {
        report.warnings.push("One terminal stop symbol was excluded explicitly; no internal stop or downstream chain was silently removed.".into());
    }
    if unresolved > 0 {
        report.warnings.push("Ambiguous residue counts remain separate; composition does not invent fractional identities. Point-valued model estimates do not establish resolved sequence identity.".into());
    }
    if options.ambiguity == AmbiguityPolicy::CanonicalBounds && x > 0 {
        report.warnings.push("X bounds assume one of the20 canonical residues, excluding U/O and modifications; intervals are conservative enclosures, not probabilities or jointly attainable extrema.".into());
    }
    if options.modification_state == ModificationState::KnownOrUnknownModifications
        || has_modification_evidence(molecule)
    {
        report.unavailable.push("Quantitative properties unavailable: declared or retained modification/cross-link evidence requires an explicit supported chemical model.".into());
        report.validate()?;
        return Ok(report);
    }
    if nonstandard > 0 {
        report.unavailable.push("U/O are retained in composition; this profile has no audited mass/ionization/absorbance model for them.".into());
        report.validate()?;
        return Ok(report);
    }
    let allow_bounds = options.ambiguity == AmbiguityPolicy::CanonicalBounds;
    if ambiguous == 0 || allow_bounds {
        report.mass = Some(MassEstimate {
            model: MASS_MODEL.into(),
            average_da: mass_bounds(&counts, &AVERAGE, &bonds),
            monoisotopic_da: mass_bounds(&counts, &MONO, &bonds),
            exact_identity: unresolved == 0,
        });
    } else {
        report.unavailable.push(
            "Mass unavailable for B/Z/X unless explicit canonical-residue bounds are selected."
                .into(),
        );
    }
    if ambiguous == 0 {
        let n = bonds.minimum as usize;
        report.isoelectric_point = pi(&counts, n);
        if report.isoelectric_point.is_none() {
            report.unavailable.push("No bracketed zero-charge root in supported pH0..14; no endpoint is reported as a pI.".into());
        }
        report.net_charge_at_ph = Some(charge(&counts, options.charge_ph, n));
    } else {
        report.unavailable.push("pI/net charge unavailable: B/Z/X leave ionizable group identities unresolved; no arbitrary average pKa is substituted.".into());
    }
    if x == 0 || allow_bounds {
        let base = 1490.0 * count(&counts, b'Y') as f64 + 5500.0 * count(&counts, b'W') as f64;
        // All unknowns as W maximize epsilon; new cystines contribute at most125 per pair.
        let epsilon = NumericInterval::new(
            base + 125.0 * bonds.minimum,
            base + 5500.0 * x as f64 + 125.0 * bonds.maximum,
        );
        let absorbance = report.mass.as_ref().map(|mass| {
            NumericInterval::new(
                epsilon.minimum / mass.average_da.maximum,
                epsilon.maximum / mass.average_da.minimum,
            )
        });
        report.extinction_280 = Some(ExtinctionEstimate {
            model: EXTINCTION_MODEL.into(),
            wavelength_nm: 280,
            epsilon_molar_inverse_cm: epsilon,
            absorbance_1mg_ml_1cm: absorbance,
        });
    } else {
        report.unavailable.push("A280 extinction unavailable for X without an explicit canonical-residue bound assumption.".into());
    }
    if matches!(options.disulfides, DisulfideState::MaximumIntrachain) && cys % 2 != 0 {
        report.warnings.push("Maximum intrachain pairing leaves an unmatched explicit Cys when no unknown residue supplies its partner; interchain disulfides are not assumed.".into());
    }
    report.warnings.push("pI/net charge use independent published Epk pKas and ideal free termini, without structural, ionic-strength, temperature or residue-context corrections. A280 is a composition estimate in water, not a measured spectrum.".into());
    if count(&counts, b'W') == 0 {
        report.warnings.push("No explicit Trp is present; the composition-based A280 estimate can have substantial relative error.".into());
    }
    report.complete = unresolved == 0 && report.unavailable.is_empty();
    report.validate()?;
    Ok(report)
}

impl ProteinReport {
    /// Structural integrity only. Imported reports are not scientific authority;
    /// consumers needing correctness must recompute analyze_protein on the source/options.
    pub fn validate(&self) -> Result<(), Error> {
        self.options.validate()?;
        let valid_hash = |s: &str| {
            s.len() == 64
                && s.bytes()
                    .all(|b| b.is_ascii_hexdigit() && !b.is_ascii_uppercase())
        };
        let mut total = 0usize;
        let mut seen = std::collections::BTreeSet::new();
        for c in &self.composition {
            if c.symbol.len() != 1
                || !SYMBOLS.contains(&c.symbol.as_bytes()[0])
                || !seen.insert(&c.symbol)
                || c.count == 0
                || !c.fraction.is_finite()
                || c.fraction != c.count as f64 / self.residue_count.max(1) as f64
            {
                return Err(Error::Format("invalid protein composition".into()));
            }
            total = total.checked_add(c.count).ok_or(Error::TooLarge)?;
        }
        if !valid_hash(&self.source_revision)
            || !valid_hash(&self.source_sequence_blake3)
            || self.engine_version != PROTEIN_ENGINE_VERSION
            || self.charge_model != CHARGE_MODEL
            || self.residue_count == 0
            || self.residue_count > self.options.max_residues
            || self.source_length != self.residue_count + usize::from(self.terminal_stop_removed)
            || total != self.residue_count
            || !self.disulfide_bonds.valid()
            || self.disulfide_bonds.minimum < 0.
            || self.disulfide_bonds.minimum.fract() != 0.
            || self.disulfide_bonds.maximum.fract() != 0.
            || self.disulfide_bonds.maximum > self.residue_count as f64 / 2.
            || self.warnings.len() > 32
            || self.unavailable.len() > 8
            || self
                .warnings
                .iter()
                .chain(&self.unavailable)
                .any(|s| s.len() > 2048)
            || self
                .isoelectric_point
                .is_some_and(|v| !v.is_finite() || !(0.0..=14.0).contains(&v))
            || self
                .net_charge_at_ph
                .is_some_and(|v| !v.is_finite() || v.abs() > self.residue_count as f64 + 2.)
        {
            return Err(Error::Format(
                "invalid protein report bounds/bindings".into(),
            ));
        }
        if let Some(m) = &self.mass {
            if m.model != MASS_MODEL
                || !m.average_da.valid()
                || !m.monoisotopic_da.valid()
                || m.average_da.minimum <= 0.
                || m.monoisotopic_da.minimum <= 0.
            {
                return Err(Error::Format("invalid protein mass interval".into()));
            }
        }
        if let Some(e) = &self.extinction_280 {
            if e.model != EXTINCTION_MODEL
                || e.wavelength_nm != 280
                || !e.epsilon_molar_inverse_cm.valid()
                || e.epsilon_molar_inverse_cm.minimum < 0.
                || e.absorbance_1mg_ml_1cm
                    .as_ref()
                    .is_some_and(|v| !v.valid() || v.minimum < 0.)
            {
                return Err(Error::Format("invalid protein extinction interval".into()));
            }
        }
        if self.complete
            && (!self.unavailable.is_empty()
                || self.mass.as_ref().is_none_or(|m| !m.exact_identity)
                || self.isoelectric_point.is_none()
                || self.net_charge_at_ph.is_none()
                || self.extinction_280.is_none())
        {
            return Err(Error::Format(
                "unresolved protein report cannot claim completeness".into(),
            ));
        }
        Ok(())
    }
}
