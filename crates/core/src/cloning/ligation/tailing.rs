//! Explicit ideal single-nucleotide 3′ tailing of an established blunt duplex.
//! This is endpoint geometry, not polymerase kinetics or a commercial kit model.
use crate::chemistry::{Duplex, Strand, StrandMapping, ThreePrime};
use collagene_model::{Error, Sequence, Topology};
use serde::{Deserialize, Serialize};

pub const TAILING_ENGINE_VERSION: &str = concat!(env!("CARGO_PKG_VERSION"), "/tailing-1");
pub const MAX_TAILING_REPORT_BYTES: usize = 32 * 1024 * 1024;
const REPORT_OVERHEAD: usize = 4096;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
pub enum TailingBase {
    A,
    C,
    G,
    T,
}
impl TailingBase {
    pub const fn as_str(self) -> &'static str {
        match self {
            Self::A => "A",
            Self::C => "C",
            Self::G => "G",
            Self::T => "T",
        }
    }
    pub const fn complement(self) -> Self {
        match self {
            Self::A => Self::T,
            Self::T => Self::A,
            Self::C => Self::G,
            Self::G => Self::C,
        }
    }
}

/// A new base without template ancestry. Coordinates address the actual product
/// strand stored 5′→3′, not a shared top-reference range or an input sequence.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct SynthesizedTail {
    pub strand: Strand,
    pub product_start: usize,
    pub product_end: usize,
    pub sequence: String,
}
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct TailingReport {
    pub engine_version: String,
    pub source_duplex_blake3: String,
    pub base: TailingBase,
    pub source: Duplex,
    pub duplex: Duplex,
    /// Retained ancestry only. Newly added residues are described by additions.
    pub mappings: Vec<StrandMapping>,
    pub additions: Vec<SynthesizedTail>,
    pub warnings: Vec<String>,
}
struct Count {
    bytes: usize,
    limit: usize,
}
impl std::io::Write for Count {
    fn write(&mut self, bytes: &[u8]) -> std::io::Result<usize> {
        self.bytes = self
            .bytes
            .checked_add(bytes.len())
            .filter(|n| *n <= self.limit)
            .ok_or_else(|| std::io::Error::other("tailing report size limit"))?;
        Ok(bytes.len())
    }
    fn flush(&mut self) -> std::io::Result<()> {
        Ok(())
    }
}
fn json_size(value: &impl Serialize) -> Result<usize, Error> {
    let mut count = Count {
        bytes: 0,
        limit: MAX_TAILING_REPORT_BYTES,
    };
    serde_json::to_writer(&mut count, value).map_err(|_| Error::TooLarge)?;
    Ok(count.bytes)
}
struct Hash(blake3::Hasher);
impl std::io::Write for Hash {
    fn write(&mut self, bytes: &[u8]) -> std::io::Result<usize> {
        self.0.update(bytes);
        Ok(bytes.len())
    }
    fn flush(&mut self) -> std::io::Result<()> {
        Ok(())
    }
}
fn source_hash(source: &Duplex) -> Result<String, Error> {
    let mut hash = Hash(blake3::Hasher::new());
    serde_json::to_writer(&mut hash, source).map_err(|e| Error::Format(e.to_string()))?;
    Ok(hash.0.finalize().to_hex().to_string())
}
fn retained(strand: Strand, length: usize) -> StrandMapping {
    StrandMapping {
        source_strand: strand,
        source_start: 0,
        source_end: length,
        product_strand: strand,
        product_start: 0,
        product_end: length,
        reverse: false,
        synthesized: false,
    }
}

/// Append exactly one declared base to both stored 3′ termini. Accessible 3′ OH
/// is required at both ends; no 5′ state, pre-existing nick or base is changed.
/// Blunt input is required, so this does not silently fill/trim a sticky end.
pub fn tail_blunt_duplex(source: &Duplex, base: TailingBase) -> Result<TailingReport, Error> {
    // All retained source/product copies are precharged before any report cloning.
    if source.nicks.len() > 10_000 {
        return Err(Error::TooLarge);
    }
    let new_length = source
        .top
        .len()
        .checked_add(1)
        .filter(|n| *n <= collagene_model::MAX_SEQUENCE_LEN)
        .ok_or(Error::TooLarge)?;
    let precharged = json_size(source)?
        .checked_mul(2)
        .and_then(|n| n.checked_add(REPORT_OVERHEAD))
        .ok_or(Error::TooLarge)?;
    if precharged > MAX_TAILING_REPORT_BYTES {
        return Err(Error::TooLarge);
    }
    source.validate()?;
    if source.topology != Topology::Linear {
        return Err(Error::Unsupported(
            "3′ tailing requires free ends on a linear duplex".into(),
        ));
    }
    if source.bottom_start != 0 || source.bottom.len() != source.top.len() {
        return Err(Error::Unsupported("single-nucleotide 3′ tailing requires a fully paired blunt duplex; repair staggered ends explicitly first".into()));
    }
    if source.left.three_prime != ThreePrime::Hydroxyl
        || source.right.three_prime != ThreePrime::Hydroxyl
    {
        return Err(Error::Unsupported("3′ tailing requires known accessible hydroxyls at both affected ends; blocked/unknown ends are not repaired by assumption".into()));
    }
    let original_length = source.top.len();
    let mut top = String::with_capacity(new_length);
    top.push_str(source.top.as_str());
    top.push_str(base.as_str());
    let mut bottom = String::with_capacity(new_length);
    bottom.push_str(source.bottom.as_str());
    bottom.push_str(base.as_str());
    let mut duplex = source.clone();
    duplex.top = Sequence::dna(&top)?;
    duplex.bottom = Sequence::dna(&bottom)?;
    // The appended bottom 3′ residue extends to the left in top-reference space.
    // Its 5′ boundary at original_length remains fixed; the top 3′ extends right.
    duplex.bottom_start = -1;
    duplex.validate()?;
    let warnings=vec![
        "Ideal endpoint transform: exactly one declared untemplated nucleotide is appended to each strand's 3′ end. Polymerase/TdT activity, tail heterogeneity, yield and experimental completion are not predicted.".into(),
        "Five-prime chemistry and pre-existing internal nicks are retained. Tailing neither phosphorylates nor seals a bond at another strand's nick, and does not establish kit-specific end chemistry.".into(),
        "Symmetric single-base tails do not establish insertion orientation or unique products. Ordered ligation must check both actual strands and explicit phosphorylation/nick policy against an independent target.".into(),
    ];
    let report = TailingReport {
        engine_version: TAILING_ENGINE_VERSION.into(),
        source_duplex_blake3: source_hash(source)?,
        base,
        source: source.clone(),
        duplex,
        mappings: vec![
            retained(Strand::Top, original_length),
            retained(Strand::Bottom, original_length),
        ],
        additions: vec![
            SynthesizedTail {
                strand: Strand::Top,
                product_start: original_length,
                product_end: new_length,
                sequence: base.as_str().into(),
            },
            SynthesizedTail {
                strand: Strand::Bottom,
                product_start: original_length,
                product_end: new_length,
                sequence: base.as_str().into(),
            },
        ],
        warnings,
    };
    json_size(&report)?;
    Ok(report)
}
impl TailingReport {
    /// Imported previews have no authority. This comparison replays the complete
    /// transform and its source binding rather than trusting derived fields.
    pub fn validate(&self) -> Result<(), Error> {
        if self.engine_version != TAILING_ENGINE_VERSION
            || self.source_duplex_blake3.len() != 64
            || self.mappings.len() != 2
            || self.additions.len() != 2
            || self.warnings.len() > 8
            || self.warnings.iter().map(String::len).sum::<usize>() > REPORT_OVERHEAD
            || self.duplex.nicks.len() > 10_000
            || self.additions.iter().any(|a| a.sequence.len() != 1)
        {
            return Err(Error::Format(
                "invalid tailing report model or bounds".into(),
            ));
        }
        json_size(self)?;
        let expected = tail_blunt_duplex(&self.source, self.base)?;
        if *self != expected {
            return Err(Error::Format(
                "tailing report does not reconstruct from its bound source/base".into(),
            ));
        }
        Ok(())
    }
}

/// Strict JSON evidence loading, including nested field preservation. Generic
/// serde decoding alone is not a source-preservation contract for nested Duplex
/// types shared with older APIs; callers importing JSON evidence use this path.
pub fn load_tailing_report(bytes: &[u8]) -> Result<TailingReport, Error> {
    if bytes.len() > MAX_TAILING_REPORT_BYTES {
        return Err(Error::TooLarge);
    }
    let value: serde_json::Value =
        serde_json::from_slice(bytes).map_err(|e| Error::Format(e.to_string()))?;
    let report: TailingReport =
        serde_json::from_slice(bytes).map_err(|e| Error::Format(e.to_string()))?;
    if serde_json::to_value(&report).map_err(|e| Error::Format(e.to_string()))? != value {
        return Err(Error::Format(
            "tailing report contains ignored or noncanonical nested fields".into(),
        ));
    }
    report.validate()?;
    Ok(report)
}
