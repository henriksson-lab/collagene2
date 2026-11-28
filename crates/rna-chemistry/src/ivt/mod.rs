//! Declared DNA template checks. This owner neither synthesizes RNA nor grants
//! physical cleavage, capping, yield, stability or experimental authority.
pub mod profiles;

use crate::{Error, Provenance, Result};
use collagene_model::{Alphabet, Molecule, Topology};
use serde::{Deserialize, Serialize};
use std::io::Write;

pub use assessment::{
    EndGeometry, EndKind, Finding, Observation, PolyAEncoding, PolyARun, Report, Rule, RunKind,
    Verdict, assess,
};
pub use profiles::{Profile, Source, profiles};
pub const ENGINE: &str = "ivt-template-1";
pub const MAX_TEMPLATE_BASES: usize = 1_000_000;
pub const MAX_REQUEST_BYTES: usize = 16 * 1024;
pub const MAX_REPORT_BYTES: usize = 1024 * 1024;
pub const MAX_RETAINED_RUNS: usize = 1024;

macro_rules! data {
    ($($item:tt)*) => {
        #[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
        #[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
        #[serde(deny_unknown_fields)]
        $($item)*
    };
}

mod assessment;

data! {
    pub struct ProfileReference {
        pub id: String,
        pub identity: String,
    }
}
data! {
    #[derive(Copy)]
    pub enum CodingStrand {
        Forward,
        Reverse,
    }
}
data! {
    pub struct Transcript {
        /// Original zero-based index of the first coding-strand base (+1).
        /// Reverse starts at this base and proceeds towards decreasing indices.
        pub start: usize,
        pub length: usize,
        /// RNA has this strand's base order (T to U); the complementary strand
        /// is the strand actually read by the polymerase.
        pub coding_strand: CodingStrand,
    }
}
data! {
    #[derive(Copy)]
    pub enum CapChoice {
        NotSpecified,
        Uncapped,
        EnzymaticPostTranscriptional,
        /// Original TriLink N-7113 only. No transfer to AG3-primeOMe/M6/AU.
        CleanCapAgN7113,
    }
}
data! {
    #[derive(Copy)]
    pub enum InitiatorRequirement {
        NoAdditionalRequirement,
        G,
        Gg,
        Ggg,
    }
}
impl InitiatorRequirement {
    pub fn required_dna(self) -> &'static str {
        match self {
            Self::NoAdditionalRequirement => "",
            Self::G => "G",
            Self::Gg => "GG",
            Self::Ggg => "GGG",
        }
    }
}
data! {
    pub struct Linearization {
        /// Both cuts are bond offsets along the coding direction from +1.
        /// offset 0 is just before +1; offset length is after the last base.
        /// These are supplied geometry, never an assertion that cutting occurred.
        pub coding_cut_offset: usize,
        pub template_cut_offset: usize,
        pub provenance: Provenance,
    }
}
data! {
    pub struct PolyARegion {
        /// Zero-based half-open offsets within the declared transcript.
        pub start_offset: usize,
        pub end_offset: usize,
    }
}
data! {
    pub struct Request {
        pub source_revision: String,
        pub profile: ProfileReference,
        pub transcript: Transcript,
        pub cap: CapChoice,
        pub initiator: InitiatorRequirement,
        #[serde(deserialize_with = "required_option")]
        #[cfg_attr(feature = "json-schema", schemars(required))]
        pub linearization: Option<Linearization>,
        #[serde(deserialize_with = "required_option")]
        #[cfg_attr(feature = "json-schema", schemars(required))]
        pub poly_a: Option<PolyARegion>,
        pub provenance: Provenance,
    }
}
fn required_option<'de, D, T>(deserializer: D) -> std::result::Result<Option<T>, D::Error>
where
    D: serde::Deserializer<'de>,
    T: Deserialize<'de>,
{
    Option::deserialize(deserializer)
}

fn provenance(value: &Provenance) -> Result<()> {
    if value.statement.trim().is_empty()
        || value.statement.len() > 2048
        || value.statement.chars().any(char::is_control)
        || value.reference.as_ref().is_some_and(|reference| {
            reference.is_empty()
                || reference.len() > 2048
                || reference.chars().any(char::is_control)
        })
    {
        return Err(Error::Invalid(
            "bounded explicit IVT declaration evidence required",
        ));
    }
    Ok(())
}
/// Borrowed DNA/model/byte admission before source revision hashing or cloning.
/// This establishes only supported source shape, not transcription/material facts.
pub fn source_preflight(source: &Molecule) -> Result<()> {
    if source.sequence.len() > MAX_TEMPLATE_BASES {
        return Err(Error::Limit("IVT template bases"));
    }
    collagene_model::bounds::model_preflight(&[source])
        .map_err(|error| Error::Source(error.to_string()))?;
    bounded(source, 32 * 1024 * 1024)?;
    source
        .validate()
        .map_err(|error| Error::Source(error.to_string()))?;
    if source.sequence.alphabet() != Alphabet::Dna {
        return Err(Error::Invalid("IVT checks require declared DNA template"));
    }
    Ok(())
}
impl Request {
    fn validate_shape(&self) -> Result<()> {
        bounded(self, MAX_REQUEST_BYTES)?;
        if self.source_revision.len() != 64
            || !self
                .source_revision
                .bytes()
                .all(|b| b.is_ascii_hexdigit() && !b.is_ascii_uppercase())
            || self.profile != profiles::current()?.reference
        {
            return Err(Error::Binding);
        }
        provenance(&self.provenance)?;
        if let Some(linearization) = &self.linearization {
            provenance(&linearization.provenance)?;
        }
        if self.transcript.length == 0 || self.transcript.length > MAX_TEMPLATE_BASES {
            return Err(Error::Invalid("nonempty bounded transcript required"));
        }
        if self.poly_a.as_ref().is_some_and(|region| {
            region.start_offset >= region.end_offset || region.end_offset > self.transcript.length
        }) {
            return Err(Error::Invalid(
                "poly(A) region is outside declared transcript",
            ));
        }
        Ok(())
    }
    pub fn validate_for(&self, source: &Molecule) -> Result<()> {
        self.validate_shape()?;
        source_preflight(source)?;
        if self.source_revision != source.revision() {
            return Err(Error::Binding);
        }
        let transcript = &self.transcript;
        let n = source.sequence.len();
        if transcript.start >= n
            || transcript.length > n
            || (source.topology == Topology::Linear
                && match transcript.coding_strand {
                    CodingStrand::Forward => transcript.length > n - transcript.start,
                    CodingStrand::Reverse => transcript.length > transcript.start + 1,
                })
        {
            return Err(Error::Invalid(
                "transcript exceeds one original source interval",
            ));
        }
        if let Some(geometry) = &self.linearization {
            for cut in [geometry.coding_cut_offset, geometry.template_cut_offset] {
                if cut > n
                    || (source.topology == Topology::Linear
                        && match transcript.coding_strand {
                            CodingStrand::Forward => cut > n - transcript.start,
                            CodingStrand::Reverse => cut > transcript.start + 1,
                        })
                {
                    return Err(Error::Invalid("supplied cut is outside original template"));
                }
            }
        }
        Ok(())
    }
}

data! {
    pub struct MappedSpan {
        /// Original forward-source half-open interval. Reverse spans are listed
        /// in coding-strand traversal order and individually reverse-oriented.
        pub start: usize,
        pub end: usize,
        pub reverse: bool,
    }
}
impl Transcript {
    pub(super) fn position(&self, source: &Molecule, offset: isize) -> Option<usize> {
        let n = source.sequence.len() as isize;
        let position = match self.coding_strand {
            CodingStrand::Forward => self.start as isize + offset,
            CodingStrand::Reverse => self.start as isize - offset,
        };
        if source.topology == Topology::Circular {
            Some(position.rem_euclid(n) as usize)
        } else {
            (0..n).contains(&position).then_some(position as usize)
        }
    }
    pub(super) fn base(&self, source: &Molecule, offset: isize) -> Option<u8> {
        let base = *source
            .sequence
            .as_str()
            .as_bytes()
            .get(self.position(source, offset)?)?;
        Some(if self.coding_strand == CodingStrand::Reverse {
            collagene_model::complement(base, Alphabet::Dna)
        } else {
            base
        })
    }
    /// Mapped source coordinates without a per-residue position allocation.
    pub fn map_interval(
        &self,
        source: &Molecule,
        start: usize,
        end: usize,
    ) -> Result<Vec<MappedSpan>> {
        if start >= end
            || end > self.length
            || self.length > source.sequence.len()
            || self.start >= source.sequence.len()
        {
            return Err(Error::Invalid(
                "IVT interval is outside declared transcript",
            ));
        }
        let first = self
            .position(source, start as isize)
            .ok_or(Error::Invalid("missing interval start"))?;
        let last = self
            .position(source, (end - 1) as isize)
            .ok_or(Error::Invalid("missing interval end"))?;
        let n = source.sequence.len();
        Ok(match self.coding_strand {
            CodingStrand::Forward if first <= last => vec![MappedSpan {
                start: first,
                end: last + 1,
                reverse: false,
            }],
            CodingStrand::Forward => vec![
                MappedSpan {
                    start: first,
                    end: n,
                    reverse: false,
                },
                MappedSpan {
                    start: 0,
                    end: last + 1,
                    reverse: false,
                },
            ],
            CodingStrand::Reverse if first >= last => vec![MappedSpan {
                start: last,
                end: first + 1,
                reverse: true,
            }],
            CodingStrand::Reverse => vec![
                MappedSpan {
                    start: 0,
                    end: first + 1,
                    reverse: true,
                },
                MappedSpan {
                    start: last,
                    end: n,
                    reverse: true,
                },
            ],
        })
    }
    pub(super) fn bond(&self, source: &Molecule, offset: usize) -> usize {
        let bond = match self.coding_strand {
            CodingStrand::Forward => self.start as isize + offset as isize,
            CodingStrand::Reverse => self.start as isize + 1 - offset as isize,
        };
        if source.topology == Topology::Circular {
            bond.rem_euclid(source.sequence.len() as isize) as usize
        } else {
            bond as usize
        }
    }
}
pub(super) fn mask(base: u8) -> u8 {
    match base {
        b'A' => 1,
        b'C' => 2,
        b'G' => 4,
        b'T' => 8,
        b'R' => 5,
        b'Y' => 10,
        b'S' => 6,
        b'W' => 9,
        b'K' => 12,
        b'M' => 3,
        b'B' => 14,
        b'D' => 13,
        b'H' => 11,
        b'V' => 7,
        b'N' => 15,
        _ => 0,
    }
}
struct BoundedWriter(usize);
impl Write for BoundedWriter {
    fn write(&mut self, bytes: &[u8]) -> std::io::Result<usize> {
        self.0 = self
            .0
            .checked_sub(bytes.len())
            .ok_or_else(|| std::io::Error::other("IVT byte bound"))?;
        Ok(bytes.len())
    }
    fn flush(&mut self) -> std::io::Result<()> {
        Ok(())
    }
}
pub(super) fn bounded(value: &impl Serialize, max: usize) -> Result<()> {
    serde_json::to_writer(BoundedWriter(max), value)
        .map_err(|_| Error::Limit("IVT serialized bytes"))
}
pub(super) fn hash(value: &impl Serialize) -> Result<String> {
    bounded(value, MAX_REPORT_BYTES)?;
    let bytes = serde_json::to_vec(value).map_err(|_| Error::Limit("IVT serialization"))?;
    Ok(blake3::hash(&bytes).to_hex().to_string())
}
pub fn revision(source: &Molecule, request: &Request) -> Result<String> {
    request.validate_for(source)?;
    hash(&(ENGINE, source.revision(), request))
}
/// Exact installed profile/rule/source inventory; no runtime source fetching.
pub fn catalogue_revision() -> Result<String> {
    hash(&(ENGINE, "profiles", profiles()?))
}
