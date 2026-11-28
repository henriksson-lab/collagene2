//! Local coordinate-bearing qualifier projection. This parses the qualifier
//! wrapper only; location syntax and rendering use the canonical INSDC codec.
use super::{Atom, MappingSegment, atoms, budget, project_location, unsupported};
use collagene_model::{Alphabet, Error, Feature, Location, Molecule};

const MAX_WRAPPER_BYTES: usize = 65_536;
const MAX_NOTE_BYTES: usize = 1024 * 1024;
const NOTE_PREFIX: &str = "Collagene source qualifier /";
const NOTE_SUFFIX: &str = " retained only as source evidence; recalculate after annotation review.";

type Qualifier = (String, Option<String>);

/// Account for source text before feature/evidence cloning or qualifier scanning.
/// Each work unit admits at most 64 bytes plus one qualifier/feature operation.
pub(super) fn preflight(source: &Molecule, work: &mut usize) -> Result<(), Error> {
    if source.features.len() > 100_000 {
        return Err(Error::TooLarge);
    }
    for feature in &source.features {
        charge(work, 1 + feature.kind.len().div_ceil(64))?;
        for (key, value) in &feature.qualifiers {
            let bytes = key
                .len()
                .checked_add(value.as_ref().map_or(0, String::len))
                .ok_or(Error::TooLarge)?;
            charge(work, 1 + bytes.div_ceil(64))?;
        }
    }
    Ok(())
}
fn charge(work: &mut usize, units: usize) -> Result<(), Error> {
    *work = work
        .checked_add(units)
        .filter(|n| *n <= super::MAX_WORK)
        .ok_or_else(|| {
            unsupported("annotation qualifier projection exceeds its shared work budget")
        })?;
    Ok(())
}
fn note(key: &str, value: &Option<String>) -> Result<Qualifier, Error> {
    let raw = value.as_deref().unwrap_or("(flag)");
    let length = NOTE_PREFIX
        .len()
        .checked_add(key.len())
        .and_then(|n| n.checked_add(1))
        .and_then(|n| n.checked_add(raw.len()))
        .and_then(|n| n.checked_add(NOTE_SUFFIX.len()))
        .filter(|n| *n <= MAX_NOTE_BYTES)
        .ok_or(Error::TooLarge)?;
    let mut text = String::with_capacity(length);
    text.push_str(NOTE_PREFIX);
    text.push_str(key);
    text.push('=');
    text.push_str(raw);
    text.push_str(NOTE_SUFFIX);
    Ok(("note".into(), Some(text)))
}
#[derive(Clone, Copy, PartialEq, Eq)]
enum Kind {
    TranslationException,
    Anticodon,
}
impl Kind {
    fn key(key: &str) -> Option<Self> {
        if key.eq_ignore_ascii_case("transl_except") {
            Some(Self::TranslationException)
        } else if key.eq_ignore_ascii_case("anticodon") {
            Some(Self::Anticodon)
        } else {
            None
        }
    }
}

pub(super) fn project(
    original: &Feature,
    source: &Molecule,
    segments: &[MappingSegment],
    destination_length: usize,
    identity: bool,
    review: bool,
    work: &mut usize,
) -> Result<Vec<Qualifier>, Error> {
    let mut result = Vec::with_capacity(original.qualifiers.len());
    let mut evidence = Vec::new();
    let mut parent: Option<Vec<Atom>> = None;
    let mut phase = None;
    for (key, value) in &original.qualifiers {
        budget(work)?;
        if !identity
            && (key.eq_ignore_ascii_case("rpt_unit_range") || key.eq_ignore_ascii_case("coded_by"))
        {
            return Err(unsupported(
                "retained rpt_unit_range/coded_by coordinates cannot be remapped",
            ));
        }
        let kind = Kind::key(key);
        let invalidate = review
            && (key.eq_ignore_ascii_case("translation")
                || key.eq_ignore_ascii_case("codon_start")
                || (!identity && kind.is_some()));
        if invalidate {
            // Inactive source evidence is never parsed into a destination claim.
            evidence.push(note(key, value)?);
        } else if !identity && let Some(kind) = kind {
            if source.sequence.alphabet() == Alphabet::Protein {
                return Err(unsupported(
                    "codon/anticodon positions require a nucleotide reference",
                ));
            }
            let raw = value.as_deref().ok_or_else(|| {
                unsupported("coordinate qualifier requires an explicit wrapper value")
            })?;
            let wrapper = Wrapper::parse(kind, raw, source.sequence.len(), work)?;
            if kind == Kind::TranslationException && phase.is_none() {
                phase = Some(codon_phase(original, work)?);
            }
            if parent.is_none() {
                parent = Some(atoms(&original.location, segments, work)?);
            }
            let mapped = wrapper.project(
                phase,
                parent.as_deref().ok_or_else(|| {
                    unsupported("parent biological coordinate proof is unavailable")
                })?,
                segments,
                source.sequence.len(),
                destination_length,
                work,
            )?;
            result.push((key.clone(), Some(mapped)));
        } else {
            result.push((key.clone(), value.clone()));
        }
    }
    result.extend(evidence);
    Ok(result)
}

struct Wrapper<'a> {
    kind: Kind,
    location: Location,
    before: &'a str,
    after: &'a str,
}
impl<'a> Wrapper<'a> {
    fn parse(kind: Kind, raw: &'a str, length: usize, work: &mut usize) -> Result<Self, Error> {
        if raw.len() > MAX_WRAPPER_BYTES {
            return Err(Error::TooLarge);
        }
        charge(work, 1 + raw.len().div_ceil(64))?;
        if !raw.is_ascii()
            || raw
                .bytes()
                .any(|b| b.is_ascii_control() && !matches!(b, b'\t' | b'\r' | b'\n'))
        {
            return Err(unsupported(
                "coordinate wrapper contains unsupported text/control bytes",
            ));
        }
        let text = raw.trim();
        if !text.starts_with('(') || !text.ends_with(')') {
            return Err(unsupported(
                "coordinate qualifier requires a complete parenthesized wrapper",
            ));
        }
        let base = raw.len() - raw.trim_start().len();
        let mut depth = 0usize;
        let mut start = base + 1;
        let mut fields = Vec::new();
        for (offset, byte) in text.bytes().enumerate() {
            match byte {
                b'(' => {
                    depth += 1;
                    if depth >= 32 {
                        return Err(unsupported("coordinate wrapper nesting exceeds 32 levels"));
                    }
                }
                b')' => {
                    depth = depth
                        .checked_sub(1)
                        .ok_or_else(|| unsupported("unbalanced coordinate wrapper"))?;
                    if depth == 0 {
                        if offset + 1 != text.len() {
                            return Err(unsupported("trailing coordinate wrapper text"));
                        }
                        fields.push((start, base + offset));
                    }
                }
                b',' if depth == 1 => {
                    if fields.len() == 3 {
                        return Err(unsupported("too many coordinate wrapper fields"));
                    }
                    fields.push((start, base + offset));
                    start = base + offset + 1;
                }
                _ => {}
            }
        }
        if depth != 0 || !(2..=3).contains(&fields.len()) {
            return Err(unsupported("invalid coordinate qualifier field structure"));
        }
        let mut position = None;
        let mut amino_acid = None;
        let mut sequence = None;
        for (start, end) in fields {
            budget(work)?;
            let field = &raw[start..end];
            let colon = field
                .find(':')
                .ok_or_else(|| unsupported("coordinate field requires a name and value"))?;
            let name = field[..colon].trim();
            let value = &field[colon + 1..];
            if value.trim().is_empty() {
                return Err(unsupported("coordinate field requires a nonempty value"));
            }
            let value_start = start + colon + 1 + value.len() - value.trim_start().len();
            let value_end = end - (value.len() - value.trim_end().len());
            if name.eq_ignore_ascii_case("pos") {
                if position.replace((value_start, value_end)).is_some() {
                    return Err(unsupported("duplicate pos field"));
                }
            } else if name.eq_ignore_ascii_case("aa") {
                if amino_acid.replace(value.trim()).is_some() {
                    return Err(unsupported("duplicate aa field"));
                }
            } else if name.eq_ignore_ascii_case("seq") && kind == Kind::Anticodon {
                if sequence.replace(value.trim()).is_some() {
                    return Err(unsupported("duplicate seq field"));
                }
            } else {
                return Err(unsupported("unknown coordinate qualifier field"));
            }
        }
        let (start, end) =
            position.ok_or_else(|| unsupported("coordinate qualifier has no pos field"))?;
        let amino_acid = amino_acid
            .filter(|value| token(value))
            .ok_or_else(|| unsupported("coordinate qualifier has no bounded aa token"))?;
        if sequence.is_some_and(|value| !token(value)) {
            return Err(unsupported("invalid anticodon seq token"));
        }
        // Normalize syntax whitespace to spaces, never concatenate separate tokens.
        let position: String = raw[start..end]
            .bytes()
            .map(|b| {
                if b.is_ascii_whitespace() {
                    ' '
                } else {
                    char::from(b)
                }
            })
            .collect();
        let location = crate::formats::sequence::parse_location(&position, length)?;
        let count = exact_size(&location, work)?;
        if count != 3
            && !(kind == Kind::TranslationException
                && amino_acid.eq_ignore_ascii_case("TERM")
                && (1..=2).contains(&count))
        {
            return Err(unsupported(
                "local codon/anticodon must span three bases; only TERM permits one or two",
            ));
        }
        Ok(Self {
            kind,
            location,
            before: &raw[..start],
            after: &raw[end..],
        })
    }
    fn project(
        &self,
        phase: Option<usize>,
        parent: &[Atom],
        segments: &[MappingSegment],
        source_length: usize,
        destination_length: usize,
        work: &mut usize,
    ) -> Result<String, Error> {
        let source_calls = calls(&self.location, source_length, work)?;
        let membership = membership(&source_calls, parent, work)?;
        if self.kind == Kind::TranslationException {
            let phase = phase
                .ok_or_else(|| unsupported("translation exception parent phase is unavailable"))?;
            if membership.start < phase || (membership.start - phase) % 3 != 0 {
                return Err(unsupported(
                    "translation exception is not aligned to the parent's CDS phase",
                ));
            }
            if source_calls.len() < 3
                && membership.start.checked_add(source_calls.len())
                    != Some(membership.parent_length)
            {
                return Err(unsupported(
                    "partial TERM exception must finish the parent's biological sequence",
                ));
            }
        }
        let (mapped, partial) = project_location(&self.location, segments, work)?;
        let mapped = mapped
            .filter(|_| !partial)
            .ok_or_else(|| unsupported("embedded codon is not fully retained"))?;
        let mapped_calls = calls(&mapped, destination_length, work)?;
        if mapped_calls != membership.calls {
            return Err(unsupported(
                "embedded codon biological order/strand changed",
            ));
        }
        let position = crate::formats::sequence::location_text(&mapped)?;
        let length = self
            .before
            .len()
            .checked_add(position.len())
            .and_then(|n| n.checked_add(self.after.len()))
            .filter(|n| *n <= MAX_WRAPPER_BYTES)
            .ok_or(Error::TooLarge)?;
        let mut result = String::with_capacity(length);
        result.push_str(self.before);
        result.push_str(&position);
        result.push_str(self.after);
        Ok(result)
    }
}
fn token(value: &str) -> bool {
    !value.is_empty()
        && value.len() <= 64
        && value.as_bytes()[0].is_ascii_alphabetic()
        && value
            .bytes()
            .all(|b| b.is_ascii_alphanumeric() || b == b'_')
}
fn exact_size(location: &Location, work: &mut usize) -> Result<usize, Error> {
    budget(work)?;
    let count = match location {
        Location::Range {
            start,
            end,
            before: false,
            after: false,
        } => end - start,
        Location::Complement(inner) => exact_size(inner, work)?,
        Location::Join(parts) => {
            let mut count = 0usize;
            for part in parts {
                count = count
                    .checked_add(exact_size(part, work)?)
                    .ok_or(Error::TooLarge)?;
                if count > 3 {
                    break;
                }
            }
            count
        }
        _ => {
            return Err(unsupported(
                "only exact local range/complement/join codon positions can be projected",
            ));
        }
    };
    if count == 0 || count > 3 {
        return Err(unsupported(
            "embedded codon exceeds three exact source calls",
        ));
    }
    Ok(count)
}
#[derive(Clone, Copy, PartialEq, Eq)]
struct Call {
    position: usize,
    reverse: bool,
}
fn calls(location: &Location, length: usize, work: &mut usize) -> Result<Vec<Call>, Error> {
    exact_size(location, work)?;
    let identity = [MappingSegment {
        source_start: 0,
        source_end: length,
        destination_start: 0,
        reverse: false,
    }];
    let mut calls = Vec::new();
    for atom in atoms(location, &identity, work)? {
        for offset in 0..atom.source_end - atom.source_start {
            budget(work)?;
            calls.push(Call {
                position: if atom.source_reverse {
                    atom.source_end - 1 - offset
                } else {
                    atom.source_start + offset
                },
                reverse: atom.source_reverse,
            });
        }
    }
    Ok(calls)
}
struct Membership {
    calls: Vec<Call>,
    start: usize,
    parent_length: usize,
}
fn codon_phase(original: &Feature, work: &mut usize) -> Result<usize, Error> {
    let mut phase = None;
    for (key, value) in &original.qualifiers {
        budget(work)?;
        if key.eq_ignore_ascii_case("codon_start") {
            let current = match value.as_deref().map(str::trim) {
                Some("1") => 0,
                Some("2") => 1,
                Some("3") => 2,
                _ => {
                    return Err(unsupported(
                        "translation exception requires a valid codon_start phase",
                    ));
                }
            };
            if phase
                .replace(current)
                .is_some_and(|previous| previous != current)
            {
                return Err(unsupported(
                    "translation exception parent has conflicting codon_start phases",
                ));
            }
        }
    }
    Ok(phase.unwrap_or(0))
}
fn membership(calls: &[Call], parent: &[Atom], work: &mut usize) -> Result<Membership, Error> {
    let mut result = Vec::with_capacity(calls.len());
    let mut previous = None;
    let mut first = None;
    let mut parent_length = 0;
    for call in calls {
        let mut ordinal = None;
        let mut prefix = 0usize;
        for atom in parent {
            budget(work)?;
            let length = atom.source_end - atom.source_start;
            if atom.source_reverse == call.reverse
                && (atom.source_start..atom.source_end).contains(&call.position)
            {
                let offset = if call.reverse {
                    atom.source_end - 1 - call.position
                } else {
                    call.position - atom.source_start
                };
                let index = prefix.checked_add(offset).ok_or(Error::TooLarge)?;
                if ordinal.replace(index).is_some() {
                    return Err(unsupported(
                        "embedded source call occurs repeatedly in the parent feature",
                    ));
                }
                result.push(Call {
                    position: if atom.destination_reverse {
                        atom.destination_end - 1 - offset
                    } else {
                        atom.destination_start + offset
                    },
                    reverse: atom.destination_reverse,
                });
            }
            prefix = prefix.checked_add(length).ok_or(Error::TooLarge)?;
        }
        let ordinal = ordinal.ok_or_else(|| {
            unsupported("embedded codon is outside the parent's biological strand")
        })?;
        if previous.is_some_and(|p: usize| p.checked_add(1) != Some(ordinal)) {
            return Err(unsupported(
                "embedded codon is not consecutive in parent biological order",
            ));
        }
        first.get_or_insert(ordinal);
        previous = Some(ordinal);
        parent_length = prefix;
    }
    Ok(Membership {
        calls: result,
        start: first.ok_or_else(|| unsupported("embedded codon has no source calls"))?,
        parent_length,
    })
}
