//! One source-bound sequence-edit ledger and fresh citation correspondence.
//! Bibliography without an exact destination scope remains ancestral evidence.
use crate::{Alphabet, Error, MAX_SEQUENCE_LEN, Molecule, Sequence};
use serde::{Deserialize, Serialize};
use std::collections::BTreeSet;

pub const KEY: &str = "genbank:sequence_edit_recipe";
pub const MAX_BYTES: usize = 1024 * 1024;
const MAX_STEPS: usize = 1024;
const MAX_INSERTED: usize = 512 * 1024;
const MAX_REPLAY_WORK: usize = 64_000_000;
const MAX_ORIGINAL: usize = 32 * 1024 * 1024;
const MAX_REFERENCES: usize = 1024;
const MAX_DESCRIPTION: usize = 4096;
const MAX_REFERENCE_SPANS: usize = 256;
const MAX_TOTAL_SPANS: usize = 4096;
const MAX_SCOPE_WORK: usize = 1_000_000;
// Map construction is at most 1+...+1024 segment visits. Existing input span ×
// edit bounds also bound intersections; this separate guard counts actual work.
const MAX_MAPPING_WORK: usize = 2_000_000;

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct SequenceEdit {
    pub start: usize,
    pub end: usize,
    pub replacement: String,
}

/// Untrusted serialized data. Only VerifiedLog establishes source/current-DNA
/// binding; the current ledger fields and version do not encode citation authority.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct EditLog {
    pub version: u8,
    pub original_record_blake3: String,
    pub source_record_index: usize,
    pub baseline_sequence_blake3: String,
    pub baseline_length: usize,
    pub edits: Vec<SequenceEdit>,
    pub result_sequence_blake3: String,
}

fn error(message: &str) -> Error {
    Error::Unsupported(format!("GenBank citation editing: {message}"))
}
fn sequence_hash(sequence: &Sequence) -> String {
    blake3::hash(sequence.as_str().as_bytes())
        .to_hex()
        .to_string()
}
fn record_index(molecule: &Molecule) -> Result<usize, Error> {
    let index = molecule
        .metadata
        .get("source_record_index")
        .map(|value| {
            value
                .parse::<usize>()
                .map_err(|_| error("invalid selected source record"))
        })
        .transpose()?
        .unwrap_or(0);
    if index > 100_000 {
        return Err(Error::TooLarge);
    }
    Ok(index)
}

impl EditLog {
    fn validate(&self, alphabet: Alphabet) -> Result<(), Error> {
        let hex = |value: &str| {
            value.len() == 64
                && value
                    .bytes()
                    .all(|b| b.is_ascii_digit() || (b'a'..=b'f').contains(&b))
        };
        if self.version != 1
            || self.baseline_length == 0
            || self.baseline_length > MAX_SEQUENCE_LEN
            || self.source_record_index > 100_000
            || self.edits.len() > MAX_STEPS
            || !hex(&self.original_record_blake3)
            || !hex(&self.baseline_sequence_blake3)
            || !hex(&self.result_sequence_blake3)
        {
            return Err(error("unsupported or invalid edit recipe"));
        }
        let mut length = self.baseline_length;
        let mut added = 0usize;
        let mut work = 0usize;
        for edit in &self.edits {
            added = added
                .checked_add(edit.replacement.len())
                .ok_or(Error::TooLarge)?;
            work = work.checked_add(length).ok_or(Error::TooLarge)?;
            if edit.start > edit.end
                || edit.end > length
                || added > MAX_INSERTED
                || work > MAX_REPLAY_WORK
            {
                return Err(error(
                    "edit log exceeds replay bounds or has invalid intervals",
                ));
            }
            if !edit.replacement.is_empty()
                && Sequence::new(alphabet, &edit.replacement)?.as_str() != edit.replacement
            {
                return Err(error("edit replacement is not canonical sequence"));
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

    /// Decode and bound an untrusted ledger without granting citation authority.
    pub fn read(molecule: &Molecule) -> Result<Option<Self>, Error> {
        let Some(raw) = molecule.metadata.get(KEY) else {
            return Ok(None);
        };
        if raw.len() > MAX_BYTES {
            return Err(Error::TooLarge);
        }
        let log: Self = serde_json::from_str(raw).map_err(|_| error("invalid edit recipe JSON"))?;
        log.validate(molecule.sequence.alphabet())?;
        Ok(Some(log))
    }

    pub(crate) fn append(
        molecule: &Molecule,
        start: usize,
        end: usize,
        replacement: &str,
        result: &Sequence,
    ) -> Result<Self, Error> {
        let raw = molecule
            .original_genbank
            .as_deref()
            .ok_or_else(|| error("original record is missing"))?;
        let descriptions = reference_descriptions(raw)?;
        let raw_hash = blake3::hash(raw.as_bytes()).to_hex().to_string();
        let index = record_index(molecule)?;
        if index != 0 {
            return Err(error("selected record is outside the single-record source"));
        }
        let current_hash = sequence_hash(&molecule.sequence);
        if start > end
            || end > molecule.sequence.len()
            || result.alphabet() != molecule.sequence.alphabet()
            || molecule.sequence.len() - (end - start) + replacement.len() != result.len()
        {
            return Err(error("invalid edit interval or result length"));
        }
        let mut log = if let Some(log) = Self::read(molecule)? {
            if log.original_record_blake3 != raw_hash
                || log.source_record_index != index
                || log.result_sequence_blake3 != current_hash
            {
                return Err(error("stale source or sequence binding"));
            }
            log
        } else {
            Self {
                version: 1,
                original_record_blake3: raw_hash,
                source_record_index: index,
                baseline_sequence_blake3: current_hash,
                baseline_length: molecule.sequence.len(),
                edits: vec![],
                result_sequence_blake3: String::new(),
            }
        };
        // Consecutive typing/backspace within the newest insertion is one bounded
        // edit. Other boundary insertions remain separate citation operations.
        if let Some(last) = log.edits.last_mut().filter(|last| {
            last.start == last.end
                && start >= last.start
                && end <= last.start + last.replacement.len()
        }) {
            last.replacement
                .replace_range(start - last.start..end - last.start, replacement);
            if last.start == last.end && last.replacement.is_empty() {
                log.edits.pop();
            }
        } else {
            log.edits.push(SequenceEdit {
                start,
                end,
                replacement: replacement.into(),
            });
        }
        log.result_sequence_blake3 = sequence_hash(result);
        log.validate(molecule.sequence.alphabet())?;
        // This admission checks framing, numeric claims and work bounds. It is
        // not export authority: independently parsed baseline replay remains required.
        log.project_references(&descriptions)?;
        Ok(log)
    }

    fn project_references<I, S>(&self, descriptions: I) -> Result<Vec<CitationProjection>, Error>
    where
        I: IntoIterator<Item = S>,
        S: AsRef<str>,
    {
        let mut result = Vec::new();
        let mut spans = 0usize;
        let mut numbers = BTreeSet::new();
        let mut retained = None;
        let mut work = MappingWork::default();
        for description in descriptions {
            if result.len() == MAX_REFERENCES {
                return Err(Error::TooLarge);
            }
            let reference = Reference::parse(description.as_ref(), self.baseline_length)?;
            if !numbers.insert(reference.number) {
                return Err(error("duplicate reference number"));
            }
            if let Scope::Bases(ref ranges) = reference.scope {
                spans = spans.checked_add(ranges.len()).ok_or(Error::TooLarge)?;
                if spans > MAX_TOTAL_SPANS
                    || spans.checked_mul(self.edits.len()).ok_or(Error::TooLarge)? > MAX_SCOPE_WORK
                {
                    return Err(Error::TooLarge);
                }
            }
            result.push(reference.project(self, &mut retained, &mut work)?);
        }
        Ok(result)
    }
}

/// An independently replayed ledger. No serializer or public field can restore
/// this authority from a caller-supplied Boolean or cached projection.
#[derive(Debug)]
pub struct VerifiedLog {
    log: EditLog,
}
impl VerifiedLog {
    pub fn verify(
        molecule: &Molecule,
        baseline: &Sequence,
        index: usize,
    ) -> Result<Option<Self>, Error> {
        let Some(log) = EditLog::read(molecule)? else {
            return Ok(None);
        };
        let raw = molecule
            .original_genbank
            .as_deref()
            .ok_or_else(|| error("missing original record"))?;
        let descriptions = reference_descriptions(raw)?;
        if index != 0
            || log.original_record_blake3 != blake3::hash(raw.as_bytes()).to_hex().as_str()
            || log.source_record_index != index
            || record_index(molecule)? != index
            || log.baseline_length != baseline.len()
            || baseline.alphabet() != molecule.sequence.alphabet()
            || log.baseline_sequence_blake3 != sequence_hash(baseline)
            || log.result_sequence_blake3 != sequence_hash(&molecule.sequence)
        {
            return Err(error("source record, baseline or result binding disagrees"));
        }
        let mut replay = baseline.as_str().to_owned();
        for edit in &log.edits {
            replay.replace_range(edit.start..edit.end, &edit.replacement);
        }
        if replay != molecule.sequence.as_str() {
            return Err(error("edit replay does not reconstruct current DNA"));
        }
        // Validate the complete source reference set even when the caller exports
        // only a subset; source-only citations never bypass ledger validation.
        log.project_references(&descriptions)?;
        Ok(Some(Self { log }))
    }
    pub fn edit_count(&self) -> usize {
        self.log.edits.len()
    }
    pub fn edits(&self) -> &[SequenceEdit] {
        &self.log.edits
    }
    /// Project descriptions from the independently parsed selected source record.
    /// This operation owns all reference/span/work bounds and preserves input order.
    pub fn project_references<I, S>(
        &self,
        descriptions: I,
    ) -> Result<Vec<CitationProjection>, Error>
    where
        I: IntoIterator<Item = S>,
        S: AsRef<str>,
    {
        self.log.project_references(descriptions)
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub enum CitationProjection {
    Mapped {
        number: u64,
        description: String,
    },
    SourceOnly {
        number: u64,
        original_description: String,
        reason: SourceOnlyReason,
    },
}
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SourceOnlyReason {
    Sites,
    Unspecified,
    UnsupportedScope,
    Deleted,
    RepresentationLimit,
}
impl SourceOnlyReason {
    pub fn description(self) -> &'static str {
        match self {
            Self::Sites => "original annotation-site citation",
            Self::Unspecified => "original citation has no positional scope",
            Self::UnsupportedScope => {
                "original citation scope is not an exact supported base interval"
            }
            Self::Deleted => "all original cited intervals were deleted",
            Self::RepresentationLimit => {
                "exact retained-source coverage exceeds the supported reference description limits"
            }
        }
    }
}

enum Scope {
    Bases(Vec<(usize, usize)>),
    SourceOnly(SourceOnlyReason),
}
struct Reference<'a> {
    number: u64,
    description: &'a str,
    scope: Scope,
}
impl<'a> Reference<'a> {
    fn parse(description: &'a str, length: usize) -> Result<Self, Error> {
        if description.len() > MAX_DESCRIPTION {
            return Err(Error::TooLarge);
        }
        if description
            .chars()
            .any(|c| c.is_control() && !c.is_whitespace())
        {
            return Err(error("invalid control character in reference description"));
        }
        let leading = description.len() - description.trim_start().len();
        let text = &description[leading..];
        let digits = text.bytes().take_while(u8::is_ascii_digit).count();
        if digits == 0 {
            return Err(error("invalid reference number"));
        }
        let number = text[..digits]
            .parse::<u64>()
            .map_err(|_| error("invalid reference number"))?;
        if number == 0 {
            return Err(error("invalid reference number"));
        }
        let after = &text[digits..];
        if after
            .chars()
            .next()
            .is_some_and(|c| !c.is_whitespace() && c != '(')
        {
            return Err(error("invalid reference number framing"));
        }
        let suffix = after.trim();
        let scope = if suffix.is_empty() {
            Scope::SourceOnly(SourceOnlyReason::Unspecified)
        } else if suffix.eq_ignore_ascii_case("(sites)") {
            Scope::SourceOnly(SourceOnlyReason::Sites)
        } else if let Some(ranges) = parse_base_scope(suffix, length)? {
            Scope::Bases(ranges)
        } else {
            Scope::SourceOnly(SourceOnlyReason::UnsupportedScope)
        };
        Ok(Self {
            number,
            description,
            scope,
        })
    }
    fn project(
        self,
        log: &EditLog,
        retained: &mut Option<RetainedMap>,
        work: &mut MappingWork,
    ) -> Result<CitationProjection, Error> {
        let mut ranges = match self.scope {
            Scope::Bases(ranges) => ranges,
            Scope::SourceOnly(reason) => {
                return Ok(CitationProjection::SourceOnly {
                    number: self.number,
                    original_description: self.description.into(),
                    reason,
                });
            }
        };
        // Source ranges must be normalized before walking the monotone map. This
        // preserves set coverage while removing duplicate/overlapping assertions.
        work.charge(ranges.len().checked_mul(9).ok_or(Error::TooLarge)?)?;
        ranges.sort_unstable();
        let mut source_ranges: Vec<(usize, usize)> = Vec::with_capacity(ranges.len());
        for (start, end) in ranges {
            if let Some(last) = source_ranges.last_mut().filter(|last| start <= last.1) {
                last.1 = last.1.max(end);
            } else {
                source_ranges.push((start, end));
            }
        }
        if retained.is_none() {
            *retained = Some(RetainedMap::build(log, work)?);
        }
        let map = retained
            .as_ref()
            .ok_or_else(|| error("missing retained-source map"))?;
        let Some(coverage) = map.coverage(&source_ranges, work)? else {
            return Ok(CitationProjection::SourceOnly {
                number: self.number,
                original_description: self.description.into(),
                reason: SourceOnlyReason::RepresentationLimit,
            });
        };
        if coverage.is_empty() {
            return Ok(CitationProjection::SourceOnly {
                number: self.number,
                original_description: self.description.into(),
                reason: SourceOnlyReason::Deleted,
            });
        }
        if coverage == source_ranges {
            return Ok(CitationProjection::Mapped {
                number: self.number,
                description: self.description.into(),
            });
        }
        // Canonical spacing avoids declaring a representation limit merely
        // because the source used excessive padding around its numeric identity.
        let description = format!(
            "{}  (bases {})",
            self.number,
            coverage
                .iter()
                .map(|&(start, end)| format!("{} to {}", start + 1, end))
                .collect::<Vec<_>>()
                .join("; ")
        );
        if description.len() > MAX_DESCRIPTION {
            return Ok(CitationProjection::SourceOnly {
                number: self.number,
                original_description: self.description.into(),
                reason: SourceOnlyReason::RepresentationLimit,
            });
        }
        Ok(CitationProjection::Mapped {
            number: self.number,
            description,
        })
    }
}

#[derive(Default)]
struct MappingWork(usize);
impl MappingWork {
    fn charge(&mut self, units: usize) -> Result<(), Error> {
        self.0 = self
            .0
            .checked_add(units)
            .filter(|&n| n <= MAX_MAPPING_WORK)
            .ok_or(Error::TooLarge)?;
        Ok(())
    }
}

#[derive(Clone, Copy)]
struct RetainedSegment {
    source_start: usize,
    source_end: usize,
    destination_start: usize,
}
impl RetainedSegment {
    fn destination_end(self) -> usize {
        self.destination_start + self.source_end - self.source_start
    }
}
struct RetainedMap {
    segments: Vec<RetainedSegment>,
}
impl RetainedMap {
    fn push(segments: &mut Vec<RetainedSegment>, segment: RetainedSegment) {
        if segment.source_start == segment.source_end {
            return;
        }
        if let Some(last) = segments.last_mut().filter(|last| {
            last.source_end == segment.source_start
                && last.destination_end() == segment.destination_start
        }) {
            last.source_end = segment.source_end;
        } else {
            segments.push(segment);
        }
    }
    /// Every edit is in the current sequence's coordinates. New bases never
    /// enter this map, even when they happen to equal deleted source letters.
    fn build(log: &EditLog, work: &mut MappingWork) -> Result<Self, Error> {
        let mut segments = vec![RetainedSegment {
            source_start: 0,
            source_end: log.baseline_length,
            destination_start: 0,
        }];
        for edit in &log.edits {
            work.charge(1)?;
            let mut next = Vec::with_capacity((segments.len() + 1).min(MAX_STEPS + 1));
            for segment in &segments {
                work.charge(1)?;
                let start = segment.destination_start;
                let end = segment.destination_end();
                if start < edit.start {
                    let left_end = end.min(edit.start);
                    Self::push(
                        &mut next,
                        RetainedSegment {
                            source_start: segment.source_start,
                            source_end: segment.source_start + left_end - start,
                            destination_start: start,
                        },
                    );
                }
                if end > edit.end {
                    let right_start = start.max(edit.end);
                    Self::push(
                        &mut next,
                        RetainedSegment {
                            source_start: segment.source_start + right_start - start,
                            source_end: segment.source_end,
                            destination_start: right_start - (edit.end - edit.start)
                                + edit.replacement.len(),
                        },
                    );
                }
            }
            if next.len() > MAX_STEPS + 1 {
                return Err(Error::TooLarge);
            }
            segments = next;
        }
        Ok(Self { segments })
    }
    /// Source and destination orders are monotone. After merging the next
    /// adjacent piece, a later piece cannot fill an earlier uncovered gap.
    fn coverage(
        &self,
        ranges: &[(usize, usize)],
        work: &mut MappingWork,
    ) -> Result<Option<Vec<(usize, usize)>>, Error> {
        let mut result: Vec<(usize, usize)> = Vec::new();
        for &(start, end) in ranges {
            let (mut low, mut high) = (0, self.segments.len());
            while low < high {
                work.charge(1)?;
                let middle = low + (high - low) / 2;
                if self.segments[middle].source_end <= start {
                    low = middle + 1;
                } else {
                    high = middle;
                }
            }
            let mut index = low;
            while let Some(segment) = self.segments.get(index).filter(|s| s.source_start < end) {
                work.charge(1)?;
                let a = start.max(segment.source_start);
                let b = end.min(segment.source_end);
                if a < b {
                    let mapped_start = segment.destination_start + a - segment.source_start;
                    let mapped_end = segment.destination_start + b - segment.source_start;
                    if let Some(last) = result.last_mut().filter(|last| mapped_start <= last.1) {
                        last.1 = last.1.max(mapped_end);
                    } else {
                        result.push((mapped_start, mapped_end));
                        if result.len() > MAX_REFERENCE_SPANS {
                            return Ok(None);
                        }
                    }
                }
                index += 1;
            }
        }
        Ok(Some(result))
    }
}

/// A recognized numeric grammar is fully consumed before its coordinate bounds
/// are interpreted. Other printable prose is source-only, never a parsed prefix.
fn parse_base_scope(scope: &str, length: usize) -> Result<Option<Vec<(usize, usize)>>, Error> {
    let Some(body) = scope.strip_prefix('(').and_then(|s| s.strip_suffix(')')) else {
        return Ok(None);
    };
    let Some(body) = body
        .strip_prefix("bases")
        .filter(|s| s.chars().next().is_some_and(char::is_whitespace))
    else {
        return Ok(None);
    };
    let mut tokens_by_span = Vec::new();
    for part in body.split(';') {
        let part = part.trim();
        let part = part
            .strip_prefix("bases")
            .filter(|s| s.chars().next().is_some_and(char::is_whitespace))
            .unwrap_or(part);
        let mut tokens = part.split_whitespace();
        let (Some(start), Some("to"), Some(end), None) =
            (tokens.next(), tokens.next(), tokens.next(), tokens.next())
        else {
            return Ok(None);
        };
        if !start.bytes().all(|b| b.is_ascii_digit()) || !end.bytes().all(|b| b.is_ascii_digit()) {
            return Ok(None);
        }
        tokens_by_span.push((start, end));
    }
    if tokens_by_span.len() > MAX_REFERENCE_SPANS {
        return Err(Error::TooLarge);
    }
    let mut spans = Vec::with_capacity(tokens_by_span.len());
    for (start, end) in tokens_by_span {
        let start = start
            .parse::<u64>()
            .map_err(|_| error("invalid reference start"))?;
        let end = end
            .parse::<u64>()
            .map_err(|_| error("invalid reference end"))?;
        if start == 0 || start > end || end > length as u64 {
            return Err(error("reference interval outside source"));
        }
        spans.push(((start - 1) as usize, end as usize));
    }
    Ok(Some(spans))
}

/// Admit source framing and reference count/description sizes before a rich
/// header parser allocates structured bibliography. This is resource admission
/// only; it establishes neither sequence replay nor citation correspondence.
pub fn preflight_record(raw: &str) -> Result<(), Error> {
    reference_descriptions(raw).map(|_| ())
}

/// Collect one bounded logical REFERENCE description. Twelve-column continuation
/// whitespace belongs to it; AUTHORS/TITLE/JOURNAL and other keyed fields do not.
fn reference_descriptions(raw: &str) -> Result<Vec<String>, Error> {
    if raw.len() > MAX_ORIGINAL {
        return Err(Error::TooLarge);
    }
    if !raw.starts_with("LOCUS")
        || raw.lines().filter(|line| line.starts_with("LOCUS")).count() != 1
        || raw
            .lines()
            .filter(|line| line.trim_end() == "ORIGIN")
            .count()
            != 1
        || raw.lines().filter(|line| line.trim_end() == "//").count() != 1
        || !raw.trim_end().ends_with("//")
        || raw.chars().any(|c| c.is_control() && !c.is_whitespace())
    {
        return Err(error(
            "ordinary single-record LOCUS/ORIGIN source is required",
        ));
    }
    let mut descriptions = Vec::new();
    let mut pending: Option<String> = None;
    let mut header = true;
    for line in raw.lines() {
        if line.starts_with("REFERENCE") {
            if !header {
                return Err(error("reference outside source header"));
            }
            if let Some(value) = pending.take() {
                descriptions.push(value);
            }
            if descriptions.len() == MAX_REFERENCES {
                return Err(Error::TooLarge);
            }
            let tail = line.strip_prefix("REFERENCE").expect("checked prefix");
            if tail.chars().next().is_some_and(|c| !c.is_whitespace()) {
                return Err(error("invalid reference header framing"));
            }
            if tail.trim().len() > MAX_DESCRIPTION {
                return Err(Error::TooLarge);
            }
            pending = Some(tail.trim().to_owned());
        } else if pending.is_some()
            && line.len() >= 12
            && line.as_bytes()[..12].iter().all(u8::is_ascii_whitespace)
        {
            let continuation = line
                .get(12..)
                .ok_or_else(|| error("invalid reference continuation"))?
                .trim();
            let value = pending.as_mut().expect("checked pending reference");
            if value
                .len()
                .checked_add(1)
                .and_then(|n| n.checked_add(continuation.len()))
                .ok_or(Error::TooLarge)?
                > MAX_DESCRIPTION
            {
                return Err(Error::TooLarge);
            }
            value.push(' ');
            value.push_str(continuation);
        } else {
            if let Some(value) = pending.take() {
                descriptions.push(value);
            }
            if line.starts_with("FEATURES") || line.trim_end() == "ORIGIN" {
                header = false;
            }
        }
    }
    if let Some(value) = pending {
        descriptions.push(value);
    }
    Ok(descriptions)
}
