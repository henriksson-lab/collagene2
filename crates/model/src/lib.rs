pub mod bounds;

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

/// Unverified source-bound edit log. It cannot authorize citation coordinates
/// without replay against the actual selected original GenBank record.
pub const GENBANK_EDIT_RECIPE_KEY: &str = "genbank:sequence_edit_recipe";
pub const MAX_GENBANK_EDIT_RECIPE_BYTES: usize = 1024 * 1024;
const MAX_GENBANK_EDIT_STEPS: usize = 1024;
const MAX_GENBANK_EDIT_BASES: usize = 512 * 1024;
const MAX_GENBANK_REPLAY_WORK: usize = 64_000_000;

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct GenbankSequenceEdit {
    pub start: usize,
    pub end: usize,
    pub replacement: String,
}
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct GenbankEditRecipe {
    pub version: u8,
    pub original_record_blake3: String,
    pub source_record_index: usize,
    pub baseline_sequence_blake3: String,
    pub baseline_length: usize,
    pub edits: Vec<GenbankSequenceEdit>,
    pub result_sequence_blake3: String,
}
fn gb_edit_error(message: &str) -> Error {
    Error::Unsupported(format!("GenBank citation editing: {message}"))
}
fn gb_sequence_hash(sequence: &Sequence) -> String {
    blake3::hash(sequence.as_str().as_bytes())
        .to_hex()
        .to_string()
}
fn gb_record_index(molecule: &Molecule) -> Result<usize, Error> {
    let index = molecule
        .metadata
        .get("source_record_index")
        .map(|value| {
            value
                .parse::<usize>()
                .map_err(|_| gb_edit_error("invalid selected source record"))
        })
        .transpose()?
        .unwrap_or(0);
    if index > 100_000 {
        return Err(Error::TooLarge);
    }
    Ok(index)
}
fn gb_recipe_shape(recipe: &GenbankEditRecipe, alphabet: Alphabet) -> Result<(), Error> {
    let hex = |value: &str| {
        value.len() == 64
            && value
                .bytes()
                .all(|b| b.is_ascii_digit() || (b'a'..=b'f').contains(&b))
    };
    if recipe.version != 1
        || recipe.baseline_length == 0
        || recipe.baseline_length > MAX_SEQUENCE_LEN
        || recipe.source_record_index > 100_000
        || recipe.edits.len() > MAX_GENBANK_EDIT_STEPS
        || !hex(&recipe.original_record_blake3)
        || !hex(&recipe.baseline_sequence_blake3)
        || !hex(&recipe.result_sequence_blake3)
    {
        return Err(gb_edit_error("unsupported or invalid edit recipe"));
    }
    let mut length = recipe.baseline_length;
    let mut added = 0usize;
    let mut work = 0usize;
    for edit in &recipe.edits {
        added = added
            .checked_add(edit.replacement.len())
            .ok_or(Error::TooLarge)?;
        work = work.checked_add(length).ok_or(Error::TooLarge)?;
        if edit.start > edit.end
            || edit.end > length
            || added > MAX_GENBANK_EDIT_BASES
            || work > MAX_GENBANK_REPLAY_WORK
        {
            return Err(gb_edit_error(
                "edit log exceeds replay bounds or has invalid intervals",
            ));
        }
        if !edit.replacement.is_empty()
            && Sequence::new(alphabet, &edit.replacement)?.as_str() != edit.replacement
        {
            return Err(gb_edit_error("edit replacement is not canonical sequence"));
        }
        length = length
            .checked_sub(edit.end - edit.start)
            .and_then(|n| n.checked_add(edit.replacement.len()))
            .ok_or(Error::TooLarge)?;
        if length == 0 || length > MAX_SEQUENCE_LEN {
            return Err(Error::TooLarge);
        }
    }
    Ok(())
}
fn gb_read_recipe(molecule: &Molecule) -> Result<Option<GenbankEditRecipe>, Error> {
    let Some(raw) = molecule.metadata.get(GENBANK_EDIT_RECIPE_KEY) else {
        return Ok(None);
    };
    if raw.len() > MAX_GENBANK_EDIT_RECIPE_BYTES {
        return Err(Error::TooLarge);
    }
    let recipe: GenbankEditRecipe =
        serde_json::from_str(raw).map_err(|_| gb_edit_error("invalid edit recipe JSON"))?;
    gb_recipe_shape(&recipe, molecule.sequence.alphabet())?;
    Ok(Some(recipe))
}
fn append_genbank_edit(
    molecule: &Molecule,
    start: usize,
    end: usize,
    replacement: &str,
    result: &Sequence,
) -> Result<GenbankEditRecipe, Error> {
    let raw = molecule
        .original_genbank
        .as_ref()
        .ok_or_else(|| gb_edit_error("original record is missing"))?;
    if raw.len() > 32 * 1024 * 1024 {
        return Err(Error::TooLarge);
    }
    if !raw.starts_with("LOCUS")
        || raw.lines().filter(|line| line.starts_with("LOCUS")).count() != 1
        || !raw.lines().any(|line| line.trim_end() == "ORIGIN")
        || !raw.trim_end().ends_with("//")
    {
        return Err(gb_edit_error(
            "ordinary single-record LOCUS/ORIGIN source is required",
        ));
    }
    // The model can accept only ordinary bounded reference headers. Structured
    // parsing and independent replay still happen at GenBank export time.
    let mut references = 0usize;
    for line in raw.lines().filter(|line| line.starts_with("REFERENCE")) {
        references += 1;
        if references > 1024 || line.len() > 4096 {
            return Err(Error::TooLarge);
        }
        let description = line.strip_prefix("REFERENCE").unwrap().trim();
        parse_genbank_reference(description, MAX_SEQUENCE_LEN)?;
    }
    let raw_hash = blake3::hash(raw.as_bytes()).to_hex().to_string();
    let record_index = gb_record_index(molecule)?;
    let current_hash = gb_sequence_hash(&molecule.sequence);
    let mut recipe = if let Some(recipe) = gb_read_recipe(molecule)? {
        if recipe.original_record_blake3 != raw_hash
            || recipe.source_record_index != record_index
            || recipe.result_sequence_blake3 != current_hash
        {
            return Err(gb_edit_error("stale source or sequence binding"));
        }
        recipe
    } else {
        GenbankEditRecipe {
            version: 1,
            original_record_blake3: raw_hash,
            source_record_index: record_index,
            baseline_sequence_blake3: current_hash,
            baseline_length: molecule.sequence.len(),
            edits: vec![],
            result_sequence_blake3: String::new(),
        }
    };
    // Normal consecutive typing/backspace inside the newest replacement is one
    // bounded edit, rather than a document-length replay for every key press.
    if let Some(last) = recipe.edits.last_mut().filter(|last| {
        last.start == last.end && start >= last.start && end <= last.start + last.replacement.len()
    }) {
        last.replacement
            .replace_range(start - last.start..end - last.start, replacement);
        if last.start == last.end && last.replacement.is_empty() {
            recipe.edits.pop();
        }
    } else {
        recipe.edits.push(GenbankSequenceEdit {
            start,
            end,
            replacement: replacement.into(),
        });
    }
    recipe.result_sequence_blake3 = gb_sequence_hash(result);
    gb_recipe_shape(&recipe, molecule.sequence.alphabet())?;
    let mut scope_work = 0usize;
    for line in raw.lines().filter(|line| line.starts_with("REFERENCE")) {
        let description = line.strip_prefix("REFERENCE").unwrap().trim();
        scope_work = scope_work
            .checked_add(description.bytes().filter(|&b| b == b';').count() + 1)
            .ok_or(Error::TooLarge)?;
        if scope_work > 4096 || scope_work.saturating_mul(recipe.edits.len()) > 1_000_000 {
            return Err(Error::TooLarge);
        }
        // A deleted citation scope is not silently assigned to unrelated bases.
        remap_genbank_reference(description, &recipe)?;
    }
    Ok(recipe)
}

/// Verify an untrusted recipe by replaying against independently parsed original
/// source DNA. Hashes alone never waive the original-record/current-DNA check.
pub fn verified_genbank_edit_recipe(
    molecule: &Molecule,
    baseline: &Sequence,
    record_index: usize,
) -> Result<Option<GenbankEditRecipe>, Error> {
    let Some(recipe) = gb_read_recipe(molecule)? else {
        return Ok(None);
    };
    let raw = molecule
        .original_genbank
        .as_ref()
        .ok_or_else(|| gb_edit_error("missing original record"))?;
    if raw.len() > 32 * 1024 * 1024 {
        return Err(Error::TooLarge);
    }
    if recipe.original_record_blake3 != blake3::hash(raw.as_bytes()).to_hex().as_str()
        || recipe.source_record_index != record_index
        || gb_record_index(molecule)? != record_index
        || recipe.baseline_length != baseline.len()
        || baseline.alphabet() != molecule.sequence.alphabet()
        || recipe.baseline_sequence_blake3 != gb_sequence_hash(baseline)
        || recipe.result_sequence_blake3 != gb_sequence_hash(&molecule.sequence)
    {
        return Err(gb_edit_error(
            "source record, baseline or result binding disagrees",
        ));
    }
    let mut replay = baseline.as_str().to_owned();
    for edit in &recipe.edits {
        replay.replace_range(edit.start..edit.end, &edit.replacement);
    }
    if replay != molecule.sequence.as_str() {
        return Err(gb_edit_error(
            "edit replay does not reconstruct current DNA",
        ));
    }
    Ok(Some(recipe))
}

fn parse_genbank_reference(
    description: &str,
    length: usize,
) -> Result<(&str, Vec<(usize, usize)>), Error> {
    if description.len() > 4096 {
        return Err(Error::TooLarge);
    }
    let (number, ranges) = description.split_once("(bases ").ok_or_else(|| {
        gb_edit_error("only literal REFERENCE n (bases a to b; c to d) is supported")
    })?;
    let index = number
        .trim()
        .parse::<usize>()
        .map_err(|_| gb_edit_error("invalid reference number"))?;
    if index == 0 {
        return Err(gb_edit_error("invalid reference number"));
    }
    let ranges = ranges
        .trim()
        .strip_suffix(')')
        .ok_or_else(|| gb_edit_error("unsupported reference scope"))?;
    let mut spans = Vec::new();
    for range in ranges.split(';') {
        if spans.len() == 256 {
            return Err(Error::TooLarge);
        }
        let (start, end) = range
            .trim()
            .split_once(" to ")
            .ok_or_else(|| gb_edit_error("unsupported reference interval"))?;
        let start = start
            .trim()
            .parse::<usize>()
            .map_err(|_| gb_edit_error("invalid reference start"))?;
        let end = end
            .trim()
            .parse::<usize>()
            .map_err(|_| gb_edit_error("invalid reference end"))?;
        if start == 0 || start > end || end > length {
            return Err(gb_edit_error("reference interval outside source"));
        }
        spans.push((start - 1, end));
    }
    Ok((number, spans))
}

/// Project ordinary citation interval envelopes with the same endpoint policy
/// as feature editing. This is correspondence, never support for inserted DNA.
pub fn remap_genbank_reference(
    description: &str,
    recipe: &GenbankEditRecipe,
) -> Result<String, Error> {
    // This coordinate utility also rejects unsafe native caller inputs. Only
    // verified_genbank_edit_recipe establishes source/current sequence binding.
    if recipe.baseline_length == 0
        || recipe.baseline_length > MAX_SEQUENCE_LEN
        || recipe.edits.len() > MAX_GENBANK_EDIT_STEPS
    {
        return Err(Error::TooLarge);
    }
    let mut length = recipe.baseline_length;
    for edit in &recipe.edits {
        if edit.start > edit.end
            || edit.end > length
            || edit.replacement.len() > MAX_GENBANK_EDIT_BASES
        {
            return Err(gb_edit_error("invalid edit interval"));
        }
        length = length
            .checked_sub(edit.end - edit.start)
            .and_then(|n| n.checked_add(edit.replacement.len()))
            .filter(|&n| n > 0 && n <= MAX_SEQUENCE_LEN)
            .ok_or(Error::TooLarge)?;
    }
    let (prefix, mut spans) = parse_genbank_reference(description, recipe.baseline_length)?;
    if spans.len().saturating_mul(recipe.edits.len()) > 1_000_000 {
        return Err(Error::TooLarge);
    }
    let original = spans.clone();
    for edit in &recipe.edits {
        let shift =
            |coordinate: usize| coordinate - (edit.end - edit.start) + edit.replacement.len();
        for span in &mut spans {
            let (start, end) = *span;
            if end <= edit.start {
                continue;
            }
            if start >= edit.end {
                *span = (shift(start), shift(end));
                continue;
            }
            *span = (
                if start < edit.start {
                    start
                } else {
                    edit.start
                },
                if end > edit.end {
                    shift(end)
                } else {
                    edit.start + edit.replacement.len()
                },
            );
        }
        spans.retain(|&(start, end)| start < end);
        if spans.is_empty() {
            return Err(gb_edit_error(
                "edit removes the entire cited interval; citation remains in the original source, but automatic export cannot assign a new scope",
            ));
        }
    }
    if spans == original {
        return Ok(description.into());
    }
    Ok(format!(
        "{prefix}(bases {})",
        spans
            .iter()
            .map(|&(start, end)| format!("{} to {}", start + 1, end))
            .collect::<Vec<_>>()
            .join("; ")
    ))
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
        // Untrusted source/recipe attachments must be bounded before rich model
        // cloning, not merely when the optional citation log is decoded later.
        if self
            .metadata
            .get(GENBANK_EDIT_RECIPE_KEY)
            .is_some_and(|raw| raw.len() > MAX_GENBANK_EDIT_RECIPE_BYTES)
            || self
                .original_genbank
                .as_ref()
                .is_some_and(|raw| raw.len() > 32 * 1024 * 1024)
        {
            return Err(Error::TooLarge);
        }
        if self.features.iter().any(|f| {
            f.qualifiers.iter().any(|(key, _)| {
                ["transl_except", "anticodon", "rpt_unit_range", "coded_by"]
                    .iter()
                    .any(|known| key.eq_ignore_ascii_case(known))
            })
        }) {
            return Err(Error::Unsupported("editing coordinate-bearing feature qualifiers is not implemented; the original document remains intact".into()));
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
        let mut bases = self.sequence.bases.clone();
        bases.replace_range(start..end, &replacement);
        let sequence = Sequence::new(self.sequence.alphabet(), &bases)?;
        let pending_recipe = if self
            .original_genbank
            .as_ref()
            .is_some_and(|raw| raw.lines().any(|line| line.starts_with("REFERENCE")))
            || self.metadata.contains_key(GENBANK_EDIT_RECIPE_KEY)
        {
            let recipe = append_genbank_edit(self, start, end, &replacement, &sequence)?;
            let encoded =
                serde_json::to_string(&recipe).map_err(|e| Error::Format(e.to_string()))?;
            if encoded.len() > MAX_GENBANK_EDIT_RECIPE_BYTES {
                return Err(Error::TooLarge);
            }
            Some(encoded)
        } else {
            None
        };
        let mut result = self.clone();
        result.sequence = sequence;
        result.features.clear();
        for feature in &self.features {
            let (location, touched) = map(&feature.location, start, end, replacement.len())?;
            if let Some(location) = location {
                let mut feature = feature.clone();
                feature.location = location;
                if touched {
                    feature.needs_review = true;
                    feature
                        .qualifiers
                        .retain(|(key, _)| !key.eq_ignore_ascii_case("translation"));
                }
                result.features.push(feature);
            }
        }
        if let Some(raw) = pending_recipe {
            result.metadata.insert(GENBANK_EDIT_RECIPE_KEY.into(), raw);
            result.metadata.insert("genbank:citation_scope".into(), "Original citations require source-bound edit replay and coordinate projection before export; they do not establish support for edited sequence.".into());
        }
        result.validate()?;
        Ok(result)
    }
}
