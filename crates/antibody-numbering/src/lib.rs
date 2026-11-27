//! Antibody reference analysis, numbering primitives and conditional IMGT rules.
//! Supplied anchors are not an automatic antibody/variable-domain detector.
use collagene_model::{Alphabet, Molecule};
use serde::{Deserialize, Serialize};
use std::{fmt, str::FromStr};
pub const ENGINE: &str = "antibody-numbering-1";
pub const MAX_SOURCE_RESIDUES: usize = 2000;
pub const MAX_CDR3_RESIDUES: usize = 2048;
const REFERENCE: &str = "https://www.imgt.org/IMGTScientificChart/Numbering/CDR3-IMGTgaps.html";
const GAP_ORDER: [u16; 11] = [111, 112, 110, 113, 109, 114, 108, 115, 107, 116, 106];
#[derive(Debug, thiserror::Error, PartialEq, Eq)]
pub enum Error {
    #[error("invalid antibody numbering input: {0}")]
    Invalid(&'static str),
    #[error("antibody numbering resource bound: {0}")]
    Limit(&'static str),
    #[error("antibody numbering source/context/report binding changed")]
    Binding,
    #[error("invalid source model: {0}")]
    Source(String),
    #[error("antibody analysis cancelled")]
    Cancelled,
}
pub type Result<T> = std::result::Result<T, Error>;
macro_rules! data {($($t:tt)*)=>{#[derive(Clone,Debug,PartialEq,Eq,Serialize,Deserialize)] #[cfg_attr(feature="json-schema",derive(schemars::JsonSchema))] #[serde(deny_unknown_fields)] $($t)*};}
pub mod domain_candidates;
pub mod domain_tracks;
pub mod imgt_loops;
pub mod numbering_windows;
pub mod profile;
pub mod references;
pub mod v_mapping;

data! {pub enum Insertion {None,Letters(String),Decimal(u16)}}
data! {
/// A label, not a floating-point quantity. No generic Ord is provided: e.g.
/// IMGT 112 insertions run in reverse insertion order along the source sequence.
pub struct Position {pub number:u16,pub insertion:Insertion}}
impl Position {
    /// Syntax only; this cannot establish membership in a biological scheme.
    pub fn validate_syntax(&self) -> Result<()> {
        if self.number == 0 || self.number > 9999 {
            return Err(Error::Invalid("position number"));
        }
        match &self.insertion {
            Insertion::None => (),
            Insertion::Letters(s)
                if !s.is_empty() && s.len() <= 4 && s.bytes().all(|b| b.is_ascii_uppercase()) =>
            {
                ()
            }
            Insertion::Decimal(n) if *n > 0 && *n <= MAX_CDR3_RESIDUES as u16 => (),
            _ => return Err(Error::Invalid("insertion label")),
        }
        Ok(())
    }
}
impl fmt::Display for Position {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.number)?;
        match &self.insertion {
            Insertion::None => Ok(()),
            Insertion::Letters(s) => write!(f, "{s}"),
            Insertion::Decimal(n) => write!(f, ".{n}"),
        }
    }
}
impl FromStr for Position {
    type Err = Error;
    fn from_str(s: &str) -> Result<Self> {
        if s.is_empty() || s.len() > 10 || !s.is_ascii() {
            return Err(Error::Invalid("position label"));
        }
        let n = s.bytes().take_while(u8::is_ascii_digit).count();
        let digits = &s[..n];
        let number = digits
            .parse::<u16>()
            .map_err(|_| Error::Invalid("position number"))?;
        if number.to_string() != digits {
            return Err(Error::Invalid("noncanonical position number"));
        }
        let suffix = &s[n..];
        let insertion = if suffix.is_empty() {
            Insertion::None
        } else if let Some(decimal) = suffix.strip_prefix('.') {
            let value = decimal
                .parse::<u16>()
                .map_err(|_| Error::Invalid("decimal insertion"))?;
            if value.to_string() != decimal {
                return Err(Error::Invalid("noncanonical insertion"));
            }
            Insertion::Decimal(value)
        } else {
            Insertion::Letters(suffix.into())
        };
        let p = Self { number, insertion };
        p.validate_syntax()?;
        Ok(p)
    }
}
fn position(number: u16) -> Position {
    Position {
        number,
        insertion: Insertion::None,
    }
}
/// IMGT CDR3 only, given an independently established length. Never sorts labels
/// lexically or numerically as floats, nor converts Kabat/Chothia labels by offset.
pub fn imgt_cdr3_positions(length: usize) -> Result<Vec<Position>> {
    if length < 2 {
        return Err(Error::Invalid(
            "CDR3 profile requires at least two residues",
        ));
    }
    if length > MAX_CDR3_RESIDUES {
        return Err(Error::Limit("CDR3 positions"));
    }
    if length <= 13 {
        let gaps = &GAP_ORDER[..13 - length];
        return Ok((105..=117)
            .filter(|n| !gaps.contains(n))
            .map(position)
            .collect());
    }
    let extra = length - 13;
    let left = extra / 2;
    let right = extra - left;
    let mut out = Vec::with_capacity(length);
    out.extend((105..=111).map(position));
    out.extend((1..=left).map(|i| Position {
        number: 111,
        insertion: Insertion::Decimal(i as u16),
    }));
    out.extend((1..=right).rev().map(|i| Position {
        number: 112,
        insertion: Insertion::Decimal(i as u16),
    }));
    out.extend((112..=117).map(position));
    Ok(out)
}
data! {pub struct Anchors {
    pub source_revision:String,
    /// Zero-based source index of the independently supplied IMGT 104 anchor.
    pub cysteine_index:usize,
    /// Zero-based source index of the independently supplied IMGT 118 anchor.
    pub j_anchor_index:usize,
    pub provenance:String,
}}
data! {pub struct NumberedResidue {pub source_index:usize,pub residue:char,pub position:Position}}
data! {pub struct ConditionalCdr3 {
    pub engine:String,pub catalog_revision:String,pub context_revision:String,pub source_revision:String,
    pub source_length:usize,pub anchors:Anchors,pub start:usize,pub end:usize,
    pub residues:Vec<NumberedResidue>,
}}
fn hash<T: Serialize>(domain: &[u8], value: &T, limit: usize) -> Result<String> {
    struct Sink {
        h: blake3::Hasher,
        n: usize,
        limit: usize,
    }
    impl std::io::Write for Sink {
        fn write(&mut self, b: &[u8]) -> std::io::Result<usize> {
            self.n = self
                .n
                .checked_add(b.len())
                .filter(|n| *n <= self.limit)
                .ok_or_else(|| std::io::Error::other("numbering JSON bound"))?;
            self.h.update(b);
            Ok(b.len())
        }
        fn flush(&mut self) -> std::io::Result<()> {
            Ok(())
        }
    }
    let mut s = Sink {
        h: blake3::Hasher::new(),
        n: 0,
        limit,
    };
    s.h.update(domain);
    serde_json::to_writer(&mut s, value).map_err(|_| Error::Limit("serialized source/context"))?;
    Ok(s.h.finalize().to_hex().to_string())
}
pub fn catalog_revision() -> String {
    hash(b"antibody-numbering-catalog-1\0",&(ENGINE,REFERENCE,GAP_ORDER,104,105,117,118,"C; F/W-G-any-residue-G; provided anchors only; add right112 then left111 alternately; source order left ascending/right descending",MAX_CDR3_RESIDUES),65536).expect("fixed numbering catalog")
}
pub fn source_revision(source: &Molecule) -> Result<String> {
    if source.sequence.alphabet() != Alphabet::Protein {
        return Err(Error::Invalid("explicit protein source required"));
    }
    if source.sequence.len() > MAX_SOURCE_RESIDUES {
        return Err(Error::Limit("source residues"));
    }
    collagene_model::bounds::ModelBudget::default()
        .molecule(source)
        .map_err(|e| Error::Source(e.to_string()))?;
    source
        .validate()
        .map_err(|e| Error::Source(e.to_string()))?;
    hash(b"", source, 32 << 20)
}
/// Conditional allocation under supplied anchors. Literal anchor checks cannot
/// establish their biological identity or antibody/domain/species classification.
pub fn number_anchored_imgt_cdr3(source: &Molecule, anchors: &Anchors) -> Result<ConditionalCdr3> {
    let revision = source_revision(source)?;
    number_validated_imgt_cdr3(source, anchors, &revision)
}
// The caller has validated this immutable source and computed its exact revision.
// Shared allocation avoids hashing the full source again for every retained pair.
fn number_validated_imgt_cdr3(
    source: &Molecule,
    anchors: &Anchors,
    revision: &str,
) -> Result<ConditionalCdr3> {
    if revision != anchors.source_revision {
        return Err(Error::Binding);
    }
    if anchors.provenance.trim().is_empty()
        || anchors.provenance.len() > 1024
        || anchors.provenance.chars().any(char::is_control)
    {
        return Err(Error::Invalid("anchor provenance"));
    }
    let b = source.sequence.as_str().as_bytes();
    let start = anchors
        .cysteine_index
        .checked_add(1)
        .ok_or(Error::Invalid("anchor index"))?;
    let end = anchors.j_anchor_index;
    let junction_end = end.checked_add(4).ok_or(Error::Invalid("J anchor index"))?;
    if start > end || junction_end > b.len() || b.get(anchors.cysteine_index) != Some(&b'C') {
        return Err(Error::Invalid(
            "anchor ordering/bounds or cysteine identity",
        ));
    }
    let j = &b[end..junction_end];
    if !matches!(j[0], b'F' | b'W')
        || j[1] != b'G'
        || j[2] == b'*'
        || j[3] != b'G'
        || b[start..end].contains(&b'*')
    {
        return Err(Error::Invalid(
            "unsupported or unresolved anchor/junction letters",
        ));
    }
    let labels = imgt_cdr3_positions(end - start)?;
    let catalog_revision = catalog_revision();
    let context_revision = hash(
        b"antibody-numbering-context-1\0",
        &(ENGINE, &revision, &catalog_revision, anchors),
        65536,
    )?;
    let residues = labels
        .into_iter()
        .enumerate()
        .map(|(i, position)| NumberedResidue {
            source_index: start + i,
            residue: b[start + i] as char,
            position,
        })
        .collect();
    Ok(ConditionalCdr3 {
        engine: ENGINE.into(),
        catalog_revision,
        context_revision,
        source_revision: revision.into(),
        source_length: b.len(),
        anchors: anchors.clone(),
        start,
        end,
        residues,
    })
}
impl ConditionalCdr3 {
    pub fn validate_for(&self, source: &Molecule, anchors: &Anchors) -> Result<()> {
        if self.engine != ENGINE
            || self.residues.len() > MAX_CDR3_RESIDUES
            || number_anchored_imgt_cdr3(source, anchors)? != *self
        {
            return Err(Error::Binding);
        }
        Ok(())
    }
}
