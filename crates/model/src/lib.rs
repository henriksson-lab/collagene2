pub mod annotations;
pub mod bounds;
pub mod decode;
pub mod genbank_edit;

use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;

pub const MAX_SEQUENCE_LEN: usize = 10_000_000;

/// Canonical column order for embedded protein scoring/matching, including O/U.
pub const PROTEIN_RESIDUES: &[u8] = b"ACDEFGHIKLMNPQRSTVWYOU";
/// IUPAC alternatives in the stored uppercase protein alphabet. A stop has no
/// amino-acid alternatives; unsupported symbols are distinct from a stop.
pub fn protein_choices(symbol: u8) -> Option<&'static [u8]> {
    match symbol {
        b'B' => Some(b"DN"),
        b'Z' => Some(b"EQ"),
        b'J' => Some(b"IL"),
        b'X' => Some(PROTEIN_RESIDUES),
        b'*' => Some(b""),
        _ => PROTEIN_RESIDUES
            .iter()
            .position(|&b| b == symbol)
            .map(|i| &PROTEIN_RESIDUES[i..i + 1]),
    }
}

#[derive(Debug, thiserror::Error)]
pub enum Error {
    #[error("Invalid {alphabet:?} character {character:?} at input position {position}")]
    InvalidBase {
        alphabet: Alphabet,
        character: char,
        position: usize,
    },
    #[error("Sequence must contain at least one base or residue")]
    EmptySequence,
    #[error("Input or result exceeds this operation's supported size limit")]
    TooLarge,
    #[error("Invalid coordinates: {0}")]
    Coordinates(String),
    #[error("Unsupported operation: {0}")]
    Unsupported(String),
    #[error("Invalid document: {0}")]
    Format(String),
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
pub enum Alphabet {
    Dna,
    Rna,
    Protein,
}

impl Alphabet {
    /// Coordinate-count label for maps and document lists; not a chemistry claim.
    pub fn map_unit(self) -> &'static str {
        match self {
            Self::Dna => "bases",
            Self::Rna => "nt",
            Self::Protein => "aa",
        }
    }
    /// Compact coordinate-count label for document summaries.
    pub fn short_unit(self) -> &'static str {
        match self {
            Self::Dna => "bp",
            Self::Rna => "nt",
            Self::Protein => "aa",
        }
    }
    pub fn selected_symbols(self) -> &'static str {
        match self {
            Self::Dna | Self::Rna => "bases",
            Self::Protein => "residues",
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
pub enum Topology {
    Linear,
    Circular,
}

/// Validated, normalized ASCII. Deserialization always goes through validation.
#[derive(Debug, Clone, PartialEq, Eq, Serialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
pub struct Sequence {
    alphabet: Alphabet,
    bases: String,
}

impl<'de> Deserialize<'de> for Sequence {
    fn deserialize<D: serde::Deserializer<'de>>(deserializer: D) -> Result<Self, D::Error> {
        #[derive(Deserialize)]
        #[serde(deny_unknown_fields)]
        struct Raw {
            alphabet: Alphabet,
            bases: String,
        }
        let raw = Raw::deserialize(deserializer)?;
        Self::new(raw.alphabet, &raw.bases).map_err(serde::de::Error::custom)
    }
}

impl Sequence {
    pub fn new(alphabet: Alphabet, input: &str) -> Result<Self, Error> {
        if input.len() > MAX_SEQUENCE_LEN.saturating_mul(2) {
            return Err(Error::TooLarge);
        }
        let valid = match alphabet {
            Alphabet::Dna => "ACGTRYSWKMBDHVN",
            Alphabet::Rna => "ACGURYSWKMBDHVN",
            Alphabet::Protein => "ACDEFGHIKLMNPQRSTVWYBXZJUO*",
        };
        let mut bases = String::with_capacity(input.len().min(MAX_SEQUENCE_LEN));
        for (i, c) in input.chars().enumerate() {
            if c.is_ascii_whitespace() {
                continue;
            }
            let c = c.to_ascii_uppercase();
            if !c.is_ascii() || !valid.contains(c) {
                return Err(Error::InvalidBase {
                    alphabet,
                    character: c,
                    position: i + 1,
                });
            }
            if bases.len() == MAX_SEQUENCE_LEN {
                return Err(Error::TooLarge);
            }
            bases.push(c);
        }
        if bases.is_empty() {
            return Err(Error::EmptySequence);
        }
        Ok(Self { alphabet, bases })
    }
    pub fn dna(input: &str) -> Result<Self, Error> {
        Self::new(Alphabet::Dna, input)
    }
    pub fn as_str(&self) -> &str {
        &self.bases
    }
    pub fn len(&self) -> usize {
        self.bases.len()
    }
    pub fn is_empty(&self) -> bool {
        self.bases.is_empty()
    }
    pub fn alphabet(&self) -> Alphabet {
        self.alphabet
    }
    pub fn reverse_complement(&self) -> Result<Self, Error> {
        if self.alphabet == Alphabet::Protein {
            return Err(Error::Unsupported("proteins have no complement".into()));
        }
        let bases = self
            .bases
            .bytes()
            .rev()
            .map(|c| complement(c, self.alphabet) as char)
            .collect();
        Ok(Self {
            alphabet: self.alphabet,
            bases,
        })
    }
    /// Fraction among unambiguous A/C/G/T/U only; unknown bases are not assumed AT.
    pub fn gc_fraction(&self) -> Option<f64> {
        if self.alphabet == Alphabet::Protein {
            return None;
        }
        let known = self.bases.bytes().filter(|c| b"ACGTU".contains(c)).count();
        let gc = self.bases.bytes().filter(|c| b"GC".contains(c)).count();
        (known > 0).then(|| gc as f64 / known as f64)
    }
}

pub fn complement(base: u8, alphabet: Alphabet) -> u8 {
    match base {
        b'A' => {
            if alphabet == Alphabet::Rna {
                b'U'
            } else {
                b'T'
            }
        }
        b'T' | b'U' => b'A',
        b'C' => b'G',
        b'G' => b'C',
        b'R' => b'Y',
        b'Y' => b'R',
        b'S' => b'S',
        b'W' => b'W',
        b'K' => b'M',
        b'M' => b'K',
        b'B' => b'V',
        b'V' => b'B',
        b'D' => b'H',
        b'H' => b'D',
        _ => b'N',
    }
}

/// Location trees preserve INSDC semantics; unresolved variants are never extracted.
#[derive(Debug, Clone, PartialEq, Eq, Serialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
pub enum Location {
    Range {
        start: usize,
        end: usize,
        before: bool,
        after: bool,
    },
    Complement(Box<Location>),
    Join(Vec<Location>),
    Order(Vec<Location>),
    /// Uninterpreted INSDC location, retained for loss-aware round trips.
    Unresolved(String),
}

impl<'de> Deserialize<'de> for Location {
    fn deserialize<D: serde::Deserializer<'de>>(deserializer: D) -> Result<Self, D::Error> {
        decode::location(deserializer)
    }
}

impl Location {
    /// A native edit which touches no range cannot break this one-piece path.
    /// Compound paths with multiple parts require the separate continuity proof.
    pub(crate) fn is_single_range(&self) -> bool {
        match self {
            Self::Range { .. } => true,
            Self::Complement(inner) => inner.is_single_range(),
            Self::Join(parts) if parts.len() == 1 => parts[0].is_single_range(),
            _ => false,
        }
    }
    pub fn range(start: usize, end: usize) -> Self {
        Self::Range {
            start,
            end,
            before: false,
            after: false,
        }
    }
    pub fn spans(&self) -> Vec<(usize, usize, bool)> {
        fn visit(l: &Location, reverse: bool, out: &mut Vec<(usize, usize, bool)>) {
            match l {
                Location::Range { start, end, .. } => out.push((*start, *end, reverse)),
                Location::Complement(inner) => visit(inner, !reverse, out),
                Location::Join(parts) | Location::Order(parts) => {
                    for p in parts {
                        visit(p, reverse, out);
                    }
                }
                Location::Unresolved(_) => {}
            }
        }
        let mut out = vec![];
        visit(self, false, &mut out);
        out
    }
    pub fn extract(&self, seq: &Sequence) -> Result<String, Error> {
        self.validate(seq.len())?;
        fn size(location: &Location, alphabet: Alphabet) -> Result<usize, Error> {
            match location {
                Location::Range {
                    start,
                    end,
                    before: false,
                    after: false,
                } => Ok(end - start),
                Location::Complement(inner) => {
                    if alphabet == Alphabet::Protein {
                        return Err(Error::Unsupported("proteins have no complement".into()));
                    }
                    size(inner, alphabet)
                }
                Location::Join(parts) => parts.iter().try_fold(0usize, |total, part| {
                    let n = size(part, alphabet)?;
                    total
                        .checked_add(n)
                        .filter(|n| *n <= MAX_SEQUENCE_LEN)
                        .ok_or(Error::TooLarge)
                }),
                _ => Err(Error::Unsupported(
                    "location does not define an exact extractable sequence".into(),
                )),
            }
        }
        let length = size(self, seq.alphabet())?;
        if length > MAX_SEQUENCE_LEN {
            return Err(Error::TooLarge);
        }
        fn write(location: &Location, seq: &Sequence, reverse: bool, out: &mut String) {
            match location {
                Location::Range { start, end, .. } => {
                    let bases = &seq.as_str()[*start..*end];
                    if reverse {
                        for b in bases.bytes().rev() {
                            out.push(complement(b, seq.alphabet()) as char);
                        }
                    } else {
                        out.push_str(bases);
                    }
                }
                Location::Complement(inner) => write(inner, seq, !reverse, out),
                Location::Join(parts) => {
                    if reverse {
                        for part in parts.iter().rev() {
                            write(part, seq, reverse, out);
                        }
                    } else {
                        for part in parts {
                            write(part, seq, reverse, out);
                        }
                    }
                }
                _ => unreachable!("exact extraction preflight"),
            }
        }
        let mut result = String::with_capacity(length);
        write(self, seq, false, &mut result);
        Ok(result)
    }
    pub fn validate(&self, len: usize) -> Result<(), Error> {
        fn check(l: &Location, len: usize, depth: usize) -> Result<(), Error> {
            if depth >= 32 {
                return Err(Error::Format("location nesting exceeds 32 levels".into()));
            }
            match l {
                Location::Range { start, end, .. } if *start >= *end || *end > len => Err(
                    Error::Coordinates(format!("feature {start}..{end} outside {len}")),
                ),
                Location::Complement(inner) => check(inner, len, depth + 1),
                Location::Join(parts) | Location::Order(parts) => {
                    if parts.is_empty() {
                        return Err(Error::Coordinates("empty compound location".into()));
                    }
                    for p in parts {
                        check(p, len, depth + 1)?;
                    }
                    Ok(())
                }
                _ => Ok(()),
            }
        }
        check(self, len, 0)
    }
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
pub struct Feature {
    pub kind: String,
    pub location: Location,
    pub qualifiers: Vec<(String, Option<String>)>,
    pub needs_review: bool,
}
impl<'de> Deserialize<'de> for Feature {
    fn deserialize<D: serde::Deserializer<'de>>(deserializer: D) -> Result<Self, D::Error> {
        decode::feature(deserializer)
    }
}

impl Feature {
    pub fn label(&self) -> &str {
        for key in ["label", "gene", "product", "note"] {
            if let Some((_, Some(value))) = self.qualifiers.iter().find(|(k, _)| k == key) {
                return value;
            }
        }
        &self.kind
    }
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
pub struct Molecule {
    pub name: String,
    pub sequence: Sequence,
    pub topology: Topology,
    pub features: Vec<Feature>,
    pub metadata: BTreeMap<String, String>,
    /// Original evidence remains immutable; it is not an instruction or current sequence.
    pub original_genbank: Option<String>,
}

impl<'de> Deserialize<'de> for Molecule {
    fn deserialize<D: serde::Deserializer<'de>>(deserializer: D) -> Result<Self, D::Error> {
        decode::molecule(deserializer)
    }
}

impl Molecule {
    pub fn dna(name: impl Into<String>, bases: &str, topology: Topology) -> Result<Self, Error> {
        Ok(Self {
            name: name.into(),
            sequence: Sequence::dna(bases)?,
            topology,
            features: vec![],
            metadata: BTreeMap::new(),
            original_genbank: None,
        })
    }
    pub fn validate(&self) -> Result<(), Error> {
        if self.name.trim().is_empty() {
            return Err(Error::Format("document name is empty".into()));
        }
        if self.sequence.alphabet() == Alphabet::Protein && self.topology == Topology::Circular {
            return Err(Error::Unsupported("circular protein documents".into()));
        }
        if self.features.len() > 100_000 {
            return Err(Error::TooLarge);
        }
        for f in &self.features {
            f.location.validate(self.sequence.len())?;
        }
        Ok(())
    }
    pub fn revision(&self) -> String {
        // Struct fields and BTreeMap have a deterministic serialization order.
        blake3::hash(&serde_json::to_vec(self).expect("serializable model"))
            .to_hex()
            .to_string()
    }
    /// Apply an edit to an explicit half-open, non-wrapping interval.
    /// Intersected annotation envelopes are retained for review; stale biological
    /// assertions become bounded source evidence rather than active destination claims.
    pub fn replace(&self, start: usize, end: usize, replacement: &str) -> Result<Self, Error> {
        if replacement.len() > MAX_SEQUENCE_LEN.saturating_mul(2) {
            return Err(Error::TooLarge);
        }
        bounds::annotation_edit_preflight(&[self])?;
        self.validate()?;
        let len = self.sequence.len();
        if start > end || end > len {
            return Err(Error::Coordinates(format!(
                "edit {start}..{end} outside {len}"
            )));
        }
        let replacement = if replacement.bytes().all(|byte| byte.is_ascii_whitespace()) {
            String::new()
        } else {
            Sequence::new(self.sequence.alphabet(), replacement)?.bases
        };
        let new_len = len - (end - start) + replacement.len();
        if new_len > MAX_SEQUENCE_LEN {
            return Err(Error::TooLarge);
        }
        // A normalized no-op preserves opaque source assertions and creates no edit recipe.
        if self.sequence.as_str()[start..end] == replacement {
            return Ok(self.clone());
        }
        let features = annotations::editing::Edit::new(start, end, replacement.len())
            .features(self, new_len)?;
        let mut bases = self.sequence.bases.clone();
        bases.replace_range(start..end, &replacement);
        let sequence = Sequence::new(self.sequence.alphabet(), &bases)?;
        let pending_recipe = if self
            .original_genbank
            .as_ref()
            .is_some_and(|raw| raw.lines().any(|line| line.starts_with("REFERENCE")))
            || self.metadata.contains_key(genbank_edit::KEY)
        {
            let recipe = genbank_edit::EditLog::append(self, start, end, &replacement, &sequence)?;
            let encoded =
                serde_json::to_string(&recipe).map_err(|e| Error::Format(e.to_string()))?;
            if encoded.len() > genbank_edit::MAX_BYTES {
                return Err(Error::TooLarge);
            }
            Some(encoded)
        } else {
            None
        };
        let mut result = self.clone();
        result.sequence = sequence;
        result.features = features;
        if let Some(raw) = pending_recipe {
            result.metadata.insert(genbank_edit::KEY.into(), raw);
            result.metadata.insert("genbank:citation_scope".into(), "Original citations require source-bound edit replay and coordinate projection before export; they do not establish support for edited sequence.".into());
        }
        bounds::annotation_edit_preflight(&[&result])?;
        result.validate()?;
        Ok(result)
    }
}
