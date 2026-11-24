use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;

pub const MAX_SEQUENCE_LEN: usize = 10_000_000;

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

std::thread_local! { static LOCATION_DECODE_DEPTH: std::cell::Cell<usize> = const { std::cell::Cell::new(0) }; }

impl<'de> Deserialize<'de> for Location {
    fn deserialize<D: serde::Deserializer<'de>>(deserializer: D) -> Result<Self, D::Error> {
        // The recursive children use Location's bounded decoder, not the unbounded raw enum.
        #[derive(Deserialize)]
        enum Raw {
            Range {
                start: usize,
                end: usize,
                before: bool,
                after: bool,
            },
            Complement(Box<Location>),
            Join(Vec<Location>),
            Order(Vec<Location>),
            Unresolved(String),
        }
        struct Guard;
        impl Drop for Guard {
            fn drop(&mut self) {
                LOCATION_DECODE_DEPTH.with(|d| d.set(d.get() - 1));
            }
        }
        LOCATION_DECODE_DEPTH.with(|depth| {
            if depth.get() >= 32 {
                return Err(serde::de::Error::custom(
                    "location nesting exceeds 32 levels",
                ));
            }
            depth.set(depth.get() + 1);
            let _guard = Guard;
            Ok(match Raw::deserialize(deserializer)? {
                Raw::Range {
                    start,
                    end,
                    before,
                    after,
                } => Location::Range {
                    start,
                    end,
                    before,
                    after,
                },
                Raw::Complement(v) => Location::Complement(v),
                Raw::Join(v) => Location::Join(v),
                Raw::Order(v) => Location::Order(v),
                Raw::Unresolved(v) => Location::Unresolved(v),
            })
        })
    }
}

impl Location {
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

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
pub struct Feature {
    pub kind: String,
    pub location: Location,
    pub qualifiers: Vec<(String, Option<String>)>,
    pub needs_review: bool,
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

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
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
    /// Intersected annotations are retained with a review flag and stale translations removed.
    pub fn replace(&self, start: usize, end: usize, replacement: &str) -> Result<Self, Error> {
        if self.features.iter().any(|f| {
            f.qualifiers.iter().any(|(key, _)| {
                matches!(
                    key.as_str(),
                    "transl_except" | "anticodon" | "rpt_unit_range" | "coded_by"
                )
            })
        }) {
            return Err(Error::Unsupported("editing coordinate-bearing feature qualifiers is not implemented; the original document remains intact".into()));
        }
        if self
            .original_genbank
            .as_ref()
            .is_some_and(|text| text.lines().any(|line| line.starts_with("REFERENCE")))
        {
            return Err(Error::Unsupported("editing GenBank reference ranges is not implemented; keep this document as an input and save designed products separately".into()));
        }
        let len = self.sequence.len();
        if start > end || end > len {
            return Err(Error::Coordinates(format!(
                "edit {start}..{end} outside {len}"
            )));
        }
        let replacement = if replacement.trim().is_empty() {
            String::new()
        } else {
            Sequence::new(self.sequence.alphabet(), replacement)?.bases
        };
        let new_len = len - (end - start) + replacement.len();
        if new_len > MAX_SEQUENCE_LEN {
            return Err(Error::TooLarge);
        }
        // Reject editing unresolved locations: preserving their stale coordinates would be misleading.
        fn map(
            l: &Location,
            a: usize,
            b: usize,
            added: usize,
        ) -> Result<(Option<Location>, bool), Error> {
            let shift = |x: usize| x - (b - a) + added;
            match l {
                Location::Range {
                    start,
                    end,
                    before,
                    after,
                } => {
                    if *end <= a {
                        return Ok((Some(l.clone()), false));
                    }
                    if *start >= b {
                        return Ok((
                            Some(Location::Range {
                                start: shift(*start),
                                end: shift(*end),
                                before: *before,
                                after: *after,
                            }),
                            false,
                        ));
                    }
                    let s = if *start < a { *start } else { a };
                    let e = if *end > b { shift(*end) } else { a + added };
                    Ok((
                        (s < e).then_some(Location::Range {
                            start: s,
                            end: e,
                            before: *before,
                            after: *after,
                        }),
                        true,
                    ))
                }
                Location::Complement(inner) => {
                    let (v, changed) = map(inner, a, b, added)?;
                    Ok((v.map(|l| Location::Complement(Box::new(l))), changed))
                }
                Location::Join(parts) | Location::Order(parts) => {
                    let mut out = vec![];
                    let mut changed = false;
                    for p in parts {
                        let (v, c) = map(p, a, b, added)?;
                        changed |= c;
                        if let Some(v) = v {
                            out.push(v);
                        }
                    }
                    Ok((
                        if out.is_empty() {
                            None
                        } else if matches!(l, Location::Join(_)) {
                            Some(Location::Join(out))
                        } else {
                            Some(Location::Order(out))
                        },
                        changed,
                    ))
                }
                Location::Unresolved(_) => Err(Error::Unsupported(
                    "editing documents with unresolved feature locations".into(),
                )),
            }
        }
        let mut result = self.clone();
        let mut bases = self.sequence.bases.clone();
        bases.replace_range(start..end, &replacement);
        result.sequence = Sequence::new(self.sequence.alphabet(), &bases)?;
        result.features.clear();
        for feature in &self.features {
            let (location, touched) = map(&feature.location, start, end, replacement.len())?;
            if let Some(location) = location {
                let mut feature = feature.clone();
                feature.location = location;
                if touched {
                    feature.needs_review = true;
                    feature.qualifiers.retain(|(key, _)| key != "translation");
                }
                result.features.push(feature);
            }
        }
        result.validate()?;
        Ok(result)
    }
}
