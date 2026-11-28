//! Independently implemented whole-template feature matching and reviewed annotation transfer.
//! Source revisions identify provenance claims, not authenticated experimental truth.
use collagene_model::{Alphabet, Error, Feature, Location, Molecule, Sequence, Topology};
use serde::{Deserialize, Serialize};
use std::{collections::BTreeSet, io};

pub const MAX_LIBRARY_BYTES: usize = 8 * 1024 * 1024;
pub const MAX_ENTRIES: usize = 1000;
pub const MAX_TEMPLATE_BASES: usize = 1_000_000;
pub const MAX_TOTAL_TEMPLATE_BASES: usize = 10_000_000;
pub const MAX_SCAN_COMPARISONS: usize = 200_000_000;
pub const MAX_SCAN_HITS: usize = 10_000;
pub const MAX_ENTRY_QUALIFIERS: usize = 500;
pub const MAX_QUALIFIER_VALUE_BYTES: usize = 1_048_576;
const MAX_REPORT_BYTES: usize = 32 * 1024 * 1024;
pub const ENGINE_VERSION: &str = concat!(env!("CARGO_PKG_VERSION"), "/feature-library-1");

fn invalid(message: impl Into<String>) -> Error {
    Error::Format(message.into())
}
fn text(value: &str, max: usize) -> bool {
    !value.trim().is_empty() && value.len() <= max && !value.chars().any(char::is_control)
}
fn hex(value: &str) -> bool {
    value.len() == 64
        && value
            .bytes()
            .all(|b| b.is_ascii_digit() || (b'a'..=b'f').contains(&b))
}
fn kind(value: &str) -> bool {
    !value.is_empty()
        && value.len() <= 64
        && value
            .bytes()
            .all(|b| b.is_ascii_alphanumeric() || b"_'-".contains(&b))
}

struct HashWriter {
    count: usize,
    limit: usize,
    hasher: blake3::Hasher,
}
impl io::Write for HashWriter {
    fn write(&mut self, bytes: &[u8]) -> io::Result<usize> {
        self.count = self
            .count
            .checked_add(bytes.len())
            .filter(|n| *n <= self.limit)
            .ok_or_else(|| io::Error::other("feature library serialization limit exceeded"))?;
        self.hasher.update(bytes);
        Ok(bytes.len())
    }
    fn flush(&mut self) -> io::Result<()> {
        Ok(())
    }
}
fn digest(value: &impl Serialize, domain: &[u8], limit: usize) -> Result<String, Error> {
    let mut writer = HashWriter {
        count: 0,
        limit,
        hasher: blake3::Hasher::new(),
    };
    writer.hasher.update(domain);
    serde_json::to_writer(&mut writer, value).map_err(|_| Error::TooLarge)?;
    Ok(writer.hasher.finalize().to_hex().to_string())
}
fn json_size(value: &impl Serialize, limit: usize) -> Result<usize, Error> {
    struct Counter {
        count: usize,
        limit: usize,
    }
    impl io::Write for Counter {
        fn write(&mut self, bytes: &[u8]) -> io::Result<usize> {
            self.count = self
                .count
                .checked_add(bytes.len())
                .filter(|n| *n <= self.limit)
                .ok_or_else(|| io::Error::other("feature report serialization limit"))?;
            Ok(bytes.len())
        }
        fn flush(&mut self) -> io::Result<()> {
            Ok(())
        }
    }
    let mut writer = Counter { count: 0, limit };
    serde_json::to_writer(&mut writer, value).map_err(|_| Error::TooLarge)?;
    Ok(writer.count)
}
fn exact_size(location: &Location) -> Result<usize, Error> {
    match location {
        Location::Range{start,end,before:false,after:false}=>end.checked_sub(*start).filter(|n|*n>0).ok_or_else(||invalid("empty/reversed library source location")),
        Location::Complement(inner)=>exact_size(inner),
        Location::Join(parts)=>parts.iter().try_fold(0usize,|sum,part|sum.checked_add(exact_size(part)?).filter(|n|*n<=MAX_TEMPLATE_BASES).ok_or(Error::TooLarge)),
        _=>Err(Error::Unsupported("library source requires an exact extractable feature location; order/fuzzy/unresolved locations are not templates".into())),
    }
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct LibrarySource {
    pub document_name: String,
    pub document_revision: String,
    pub feature_index: usize,
    pub location: Location,
    pub source_needs_review: bool,
}
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct LibraryEntry {
    pub id: String,
    pub revision: String,
    pub name: String,
    pub sequence: Sequence,
    pub kind: String,
    /// Exact original qualifiers remain evidence. Transfer uses a conservative allowlist.
    pub qualifiers: Vec<(String, Option<String>)>,
    pub source: LibrarySource,
}
impl LibraryEntry {
    fn content_revision(&self) -> Result<String, Error> {
        digest(
            &(
                1u16,
                &self.name,
                &self.sequence,
                &self.kind,
                &self.qualifiers,
                &self.source,
            ),
            b"collagene-feature-entry-v1\0",
            MAX_LIBRARY_BYTES,
        )
    }
    fn validate(&self) -> Result<(), Error> {
        if !text(&self.name, 1024)
            || !kind(&self.kind)
            || self.sequence.alphabet() != Alphabet::Dna
            || self.sequence.len() > MAX_TEMPLATE_BASES
            || !text(&self.source.document_name, 1024)
            || !hex(&self.source.document_revision)
            || self.source.feature_index >= 100_000
            || self.qualifiers.len() > MAX_ENTRY_QUALIFIERS
        {
            return Err(invalid("invalid feature-library entry fields or limits"));
        }
        self.source
            .location
            .validate(collagene_model::MAX_SEQUENCE_LEN)?;
        if exact_size(&self.source.location)? != self.sequence.len() {
            return Err(invalid(
                "library template length disagrees with its claimed source location",
            ));
        }
        for (key, value) in &self.qualifiers {
            if key.is_empty()
                || key.len() > 128
                || !key
                    .bytes()
                    .all(|b| b.is_ascii_alphanumeric() || b"_:-".contains(&b))
                || value.as_ref().is_some_and(|value| {
                    value.len() > MAX_QUALIFIER_VALUE_BYTES
                        || value
                            .chars()
                            .any(|c| c.is_control() && !matches!(c, '\n' | '\r' | '\t'))
                })
            {
                return Err(invalid(
                    "invalid original feature-library qualifier or qualifier limit",
                ));
            }
        }
        let revision = self.content_revision()?;
        if self.revision != revision || self.id != format!("entry-{revision}") {
            return Err(invalid("feature-library entry content hash or ID mismatch"));
        }
        Ok(())
    }
}
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct FeatureLibrary {
    pub schema_version: u16,
    pub name: String,
    pub version: String,
    pub entries: Vec<LibraryEntry>,
}
impl FeatureLibrary {
    pub fn new(name: impl Into<String>, version: impl Into<String>) -> Result<Self, Error> {
        let library = Self {
            schema_version: 1,
            name: name.into(),
            version: version.into(),
            entries: vec![],
        };
        library.validate()?;
        Ok(library)
    }
    /// Content identity; call validate before using an imported library.
    pub fn revision(&self) -> String {
        digest(self, b"collagene-feature-library-v1\0", usize::MAX)
            .expect("feature library is JSON serializable")
    }
    pub fn validate(&self) -> Result<(), Error> {
        if self.schema_version != 1 {
            return Err(Error::Unsupported(
                "unknown feature-library schema version".into(),
            ));
        }
        if !text(&self.name, 1024) || !text(&self.version, 128) || self.entries.len() > MAX_ENTRIES
        {
            return Err(invalid("invalid library name/version or entry limit"));
        }
        let mut ids = BTreeSet::new();
        let mut bases = 0usize;
        for entry in &self.entries {
            entry.validate()?;
            if !ids.insert(&entry.id) {
                return Err(invalid("duplicate feature-library entry ID"));
            }
            bases = bases
                .checked_add(entry.sequence.len())
                .filter(|n| *n <= MAX_TOTAL_TEMPLATE_BASES)
                .ok_or(Error::TooLarge)?;
        }
        digest(self, b"collagene-feature-library-size\0", MAX_LIBRARY_BYTES)?;
        Ok(())
    }
    /// Biological extraction respects reverse, joined and origin-spanning exact locations.
    /// Re-adding the identical source feature is idempotent. No whole source document is copied.
    pub fn add_feature(&mut self, source: &Molecule, index: usize) -> Result<String, Error> {
        self.validate()?;
        source.validate()?;
        crate::protocol::save_document(source)?;
        if source.sequence.alphabet() != Alphabet::Dna {
            return Err(Error::Unsupported(
                "feature-library templates currently require DNA".into(),
            ));
        }
        let feature = source
            .features
            .get(index)
            .ok_or_else(|| Error::Coordinates("source annotation no longer exists".into()))?;
        if exact_size(&feature.location)? > MAX_TEMPLATE_BASES {
            return Err(Error::TooLarge);
        }
        let mut entry = LibraryEntry {
            id: String::new(),
            revision: String::new(),
            name: feature.label().into(),
            sequence: Sequence::dna(&feature.location.extract(&source.sequence)?)?,
            kind: feature.kind.clone(),
            qualifiers: feature.qualifiers.clone(),
            source: LibrarySource {
                document_name: source.name.clone(),
                document_revision: source.revision(),
                feature_index: index,
                location: feature.location.clone(),
                source_needs_review: feature.needs_review,
            },
        };
        entry.revision = entry.content_revision()?;
        entry.id = format!("entry-{}", entry.revision);
        entry.validate()?;
        if self.entries.iter().any(|existing| existing.id == entry.id) {
            return Ok(entry.id);
        }
        let id = entry.id.clone();
        let mut proposed = self.clone();
        proposed.entries.push(entry);
        proposed.validate()?;
        *self = proposed;
        Ok(id)
    }
}
/// Strict, bounded JSON: future/unknown fields are rejected, including nested model fields.
pub fn load_library(bytes: &[u8]) -> Result<FeatureLibrary, Error> {
    if bytes.len() > MAX_LIBRARY_BYTES {
        return Err(Error::TooLarge);
    }
    let library: FeatureLibrary = serde_json::from_slice(bytes)
        .map_err(|e| invalid(format!("invalid feature-library JSON: {e}")))?;
    library.validate()?;
    let raw: serde_json::Value =
        serde_json::from_slice(bytes).map_err(|e| invalid(e.to_string()))?;
    if serde_json::to_value(&library).map_err(|e| invalid(e.to_string()))? != raw {
        return Err(invalid(
            "feature-library JSON contains unsupported fields or noncanonical sequence data",
        ));
    }
    Ok(library)
}
pub fn save_library(library: &FeatureLibrary) -> Result<Vec<u8>, Error> {
    library.validate()?;
    serde_json::to_vec(library).map_err(|e| invalid(e.to_string()))
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct ScanOptions {
    pub both_strands: bool,
    /// Definitely incompatible columns; ambiguous-compatible columns remain uncertain.
    pub max_mismatches: usize,
    pub max_comparisons: usize,
    pub max_hits: usize,
}
impl Default for ScanOptions {
    fn default() -> Self {
        Self {
            both_strands: true,
            max_mismatches: 0,
            max_comparisons: 20_000_000,
            max_hits: 1000,
        }
    }
}
impl ScanOptions {
    fn validate(&self) -> Result<(), Error> {
        if self.max_mismatches > 1000
            || self.max_comparisons == 0
            || self.max_comparisons > MAX_SCAN_COMPARISONS
            || self.max_hits == 0
            || self.max_hits > MAX_SCAN_HITS
        {
            Err(invalid("invalid feature-library scan limits"))
        } else {
            Ok(())
        }
    }
}
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct AnnotationProposal {
    pub id: String,
    pub entry_id: String,
    pub entry_revision: String,
    /// Start on the forward reference, zero-based; an origin-spanning hit may wrap.
    pub start: usize,
    pub length: usize,
    pub reverse: bool,
    pub location: Location,
    pub concrete_matches: usize,
    pub compatible_ambiguous: usize,
    pub mismatches: usize,
    pub identity_lower: f64,
    pub identity_upper: f64,
    pub coverage: f64,
    pub already_annotated: bool,
    pub feature: Feature,
    pub warnings: Vec<String>,
}
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct ScanReport {
    pub source_revision: String,
    pub library_revision: String,
    pub options: ScanOptions,
    pub engine_version: String,
    pub complete: bool,
    pub comparisons: usize,
    pub hits: Vec<AnnotationProposal>,
    pub warnings: Vec<String>,
}
fn mask(base: u8) -> u8 {
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
fn location(start: usize, length: usize, n: usize, reverse: bool) -> Location {
    let forward = if start + length <= n {
        Location::range(start, start + length)
    } else {
        Location::Join(vec![
            Location::range(start, n),
            Location::range(0, start + length - n),
        ])
    };
    if reverse {
        Location::Complement(Box::new(forward))
    } else {
        forward
    }
}
fn signature(location: &Location) -> Option<Vec<(usize, usize, bool)>> {
    exact_size(location).ok()?;
    let mut spans = location.spans();
    spans.sort_unstable();
    let mut merged = Vec::<(usize, usize, bool)>::new();
    for (start, end, reverse) in spans {
        if let Some(previous) = merged.last_mut() {
            if previous.1 == start && previous.2 == reverse {
                previous.1 = end;
                continue;
            }
        }
        merged.push((start, end, reverse));
    }
    Some(merged)
}
type FeatureKey = (String, String, Vec<(usize, usize, bool)>);
fn feature_key(feature: &Feature) -> Option<FeatureKey> {
    Some((
        feature.kind.clone(),
        feature.label().to_owned(),
        signature(&feature.location)?,
    ))
}
fn descriptive(key: &str) -> bool {
    matches!(
        key,
        "label"
            | "gene"
            | "product"
            | "note"
            | "standard_name"
            | "regulatory_class"
            | "operon"
            | "allele"
            | "gene_synonym"
            | "function"
            | "color"
            | "ApEinfo_fwdcolor"
            | "ApEinfo_revcolor"
            | "ApEinfo_graphicformat"
    )
}
fn proposal(
    molecule: &Molecule,
    existing: &BTreeSet<FeatureKey>,
    library: &FeatureLibrary,
    entry: &LibraryEntry,
    palindromic: bool,
    report: &ScanReport,
    start: usize,
    reverse: bool,
    concrete: usize,
    ambiguous: usize,
    mismatches: usize,
) -> Result<AnnotationProposal, Error> {
    let length = entry.sequence.len();
    let location = location(start, length, molecule.sequence.len(), reverse);
    let id = format!(
        "annotation-{}",
        digest(
            &(
                ENGINE_VERSION,
                &report.source_revision,
                &report.library_revision,
                &entry.revision,
                &report.options,
                start,
                length,
                reverse,
                concrete,
                ambiguous,
                mismatches
            ),
            b"collagene-annotation-proposal-v1\0",
            MAX_REPORT_BYTES
        )?
    );
    let mut warnings=vec!["Whole-template, ungapped reference match only; source provenance and biological function are not independently authenticated. All transferred annotations require review.".into()];
    if palindromic {
        warnings.push("The template is its own reverse complement: sequence alone cannot resolve annotation orientation. Both orientations remain separate proposals when both strands are searched.".into());
    }
    if ambiguous > 0 {
        warnings.push("IUPAC-compatible columns do not establish concrete base identity. Identity is reported as a lower/upper bound.".into());
    }
    if mismatches > 0 {
        warnings.push("The proposed feature differs from the library template; sequence-specific source assertions are not transferred.".into());
    }
    if entry.source.source_needs_review {
        warnings.push("The source annotation was already marked for review.".into());
    }
    let removed = entry
        .qualifiers
        .iter()
        .filter(|(key, _)| !descriptive(key))
        .map(|(key, _)| key.clone())
        .collect::<BTreeSet<_>>();
    if !removed.is_empty() {
        warnings.push(format!(
            "Source-only qualifiers were retained in the library instead of transferred: {}.",
            removed.into_iter().collect::<Vec<_>>().join(", ")
        ));
    }
    let mut qualifiers = entry
        .qualifiers
        .iter()
        .filter(|(key, _)| descriptive(key))
        .cloned()
        .collect::<Vec<_>>();
    if !qualifiers.iter().any(|(key, _)| key == "label") {
        qualifiers.insert(0, ("label".into(), Some(entry.name.clone())));
    }
    for (key, value) in [
        ("library_entry_id", entry.id.clone()),
        ("library_entry_revision", entry.revision.clone()),
        ("library_revision", report.library_revision.clone()),
        ("library_name", library.name.clone()),
        ("library_version", library.version.clone()),
        (
            "library_source_revision",
            entry.source.document_revision.clone(),
        ),
        (
            "library_source_feature",
            entry.source.feature_index.to_string(),
        ),
        ("library_proposal_id", id.clone()),
        ("library_concrete_matches", concrete.to_string()),
        ("library_ambiguous_columns", ambiguous.to_string()),
        ("library_mismatches", mismatches.to_string()),
        ("library_template_length", length.to_string()),
        (
            "library_match_model",
            "whole-template ungapped; function unverified".into(),
        ),
    ] {
        qualifiers.push((key.into(), Some(value)));
    }
    let feature = Feature {
        kind: entry.kind.clone(),
        location: location.clone(),
        qualifiers,
        needs_review: true,
    };
    let already_annotated = feature_key(&feature).is_some_and(|key| existing.contains(&key));
    Ok(AnnotationProposal {
        id,
        entry_id: entry.id.clone(),
        entry_revision: entry.revision.clone(),
        start,
        length,
        reverse,
        location,
        concrete_matches: concrete,
        compatible_ambiguous: ambiguous,
        mismatches,
        identity_lower: concrete as f64 / length as f64,
        identity_upper: (concrete + ambiguous) as f64 / length as f64,
        coverage: 1.,
        already_annotated,
        feature,
        warnings,
    })
}
/// Exhaustive fixed-length placements in entry/start/forward-then-reverse order until a
/// declared bound is reached. No partial alignment or unsearched position becomes evidence.
pub fn scan(
    molecule: &Molecule,
    library: &FeatureLibrary,
    options: &ScanOptions,
) -> Result<ScanReport, Error> {
    molecule.validate()?;
    library.validate()?;
    options.validate()?;
    if molecule.sequence.alphabet() != Alphabet::Dna {
        return Err(Error::Unsupported(
            "feature-library matching currently requires a DNA reference".into(),
        ));
    }
    digest(
        molecule,
        b"collagene-feature-scan-document-size\0",
        MAX_REPORT_BYTES,
    )?;
    let mut report=ScanReport{source_revision:molecule.revision(),library_revision:library.revision(),options:options.clone(),engine_version:ENGINE_VERSION.into(),complete:true,comparisons:0,hits:vec![],warnings:vec!["Only whole-template, ungapped matches are searched. Partial/homology alignments, indels and functional validation are outside this matching model. Imported source revisions are provenance claims, not authentication.".into()]};
    let existing = molecule
        .features
        .iter()
        .filter_map(feature_key)
        .collect::<BTreeSet<_>>();
    let mut hit_bytes = 0usize;
    let bases = molecule.sequence.as_str().as_bytes();
    let n = bases.len();
    'entries: for entry in &library.entries {
        let length = entry.sequence.len();
        if length > n {
            report.warnings.push(format!("Entry '{}' is longer than the reference; no single-span/one-turn placement is possible.",entry.name));
            continue;
        }
        let reverse_sequence = entry.sequence.reverse_complement()?;
        let palindromic = reverse_sequence == entry.sequence;
        let reverse = if options.both_strands {
            Some(reverse_sequence)
        } else {
            None
        };
        let starts = if molecule.topology == Topology::Circular {
            n
        } else {
            n - length + 1
        };
        for start in 0..starts {
            for (reverse, template) in std::iter::once((false, &entry.sequence))
                .chain(reverse.as_ref().map(|sequence| (true, sequence)))
            {
                let (mut concrete, mut ambiguous, mut mismatches) = (0, 0, 0);
                for (offset, query) in template.as_str().bytes().enumerate() {
                    if report.comparisons == options.max_comparisons {
                        report.complete = false;
                        report.warnings.push("Base-comparison budget exhausted. Unfinished placements and all later work were not reported as matches.".into());
                        break 'entries;
                    }
                    report.comparisons += 1;
                    let reference = bases[(start + offset) % n];
                    let a = mask(query);
                    let b = mask(reference);
                    if a & b == 0 {
                        mismatches += 1;
                    } else if a.count_ones() == 1 && b.count_ones() == 1 {
                        concrete += 1;
                    } else {
                        ambiguous += 1;
                    }
                    if mismatches > options.max_mismatches {
                        break;
                    }
                }
                if mismatches <= options.max_mismatches {
                    if report.hits.len() == options.max_hits {
                        report.complete = false;
                        report.warnings.push("A further qualifying placement exceeded the hit limit; additional hits are unreported.".into());
                        break 'entries;
                    }
                    let hit = proposal(
                        molecule,
                        &existing,
                        library,
                        entry,
                        palindromic,
                        &report,
                        start,
                        reverse,
                        concrete,
                        ambiguous,
                        mismatches,
                    )?;
                    let size = json_size(&hit, MAX_REPORT_BYTES)?;
                    if hit_bytes.saturating_add(size).saturating_add(1) > MAX_REPORT_BYTES {
                        report.complete = false;
                        report.warnings.push("Annotation evidence reached the 32 MiB report limit; the next and later hits are unreported.".into());
                        break 'entries;
                    }
                    hit_bytes += size + 1;
                    report.hits.push(hit);
                }
            }
        }
    }
    let mut counts = std::collections::BTreeMap::<String, usize>::new();
    for hit in &report.hits {
        *counts.entry(hit.entry_id.clone()).or_default() += 1;
    }
    for hit in &mut report.hits {
        if counts[&hit.entry_id] > 1 {
            hit.warnings.push("Multiple locations/orientations were found for this entry. Each proposal is independent; no unique functional placement is asserted.".into());
        }
    }
    let repeated = counts.values().filter(|count| **count > 1).count();
    if repeated > 0 {
        report.warnings.push(format!("Multiple locations/orientations were found for {repeated} library entries; incomplete searches may hide additional placements."));
    }
    if library.entries.is_empty() {
        report.warnings.push("The library has no entries.".into());
    }
    // Measure once after adding warnings. Tail trimming serializes only removed
    // hits, so evidence accounting stays linear rather than hashing every prefix.
    let mut size = json_size(&report, usize::MAX)?;
    if size > MAX_REPORT_BYTES {
        report.complete = false;
        report
            .warnings
            .push("Report evidence limit omitted tail proposals; output is incomplete.".into());
        size = json_size(&report, usize::MAX)?;
        while size > MAX_REPORT_BYTES {
            let hit = report.hits.pop().ok_or(Error::TooLarge)?;
            size -= json_size(&hit, usize::MAX)? + usize::from(!report.hits.is_empty());
        }
    }
    debug_assert!(json_size(&report, MAX_REPORT_BYTES).is_ok());
    Ok(report)
}
/// Apply selected review-required annotations atomically, with full reconstruction of
/// match evidence. Imported reports/proposals cannot authorize arbitrary annotations.
pub fn apply_selected(
    molecule: &Molecule,
    library: &FeatureLibrary,
    report: &ScanReport,
    selected: &[String],
) -> Result<Molecule, Error> {
    molecule.validate()?;
    library.validate()?;
    report.options.validate()?;
    if report.hits.len() > MAX_SCAN_HITS {
        return Err(Error::TooLarge);
    }
    for hit in &report.hits {
        hit.location.validate(molecule.sequence.len())?;
        hit.feature.location.validate(molecule.sequence.len())?;
    }
    if selected.is_empty() || selected.len() > MAX_SCAN_HITS {
        return Err(invalid("select one or more bounded annotation proposals"));
    }
    if report.source_revision != molecule.revision()
        || report.library_revision != library.revision()
        || report.engine_version != ENGINE_VERSION
    {
        return Err(invalid(
            "annotation proposal source/library/engine is stale",
        ));
    }
    digest(report, b"collagene-feature-report-size\0", MAX_REPORT_BYTES)?;
    let rebuilt = scan(molecule, library, &report.options)?;
    if *report != rebuilt {
        return Err(invalid(
            "annotation report disagrees with reconstructed matches; imported status or edited proposals are not authority",
        ));
    }
    let proposals = rebuilt
        .hits
        .iter()
        .map(|hit| (hit.id.as_str(), hit))
        .collect::<std::collections::BTreeMap<_, _>>();
    let mut seen = BTreeSet::new();
    let mut existing = molecule
        .features
        .iter()
        .filter_map(feature_key)
        .collect::<BTreeSet<_>>();
    let mut result = molecule.clone();
    for id in selected {
        if !seen.insert(id) {
            return Err(invalid("duplicate selected annotation proposal"));
        }
        let proposal = proposals
            .get(id.as_str())
            .copied()
            .ok_or_else(|| invalid("selected annotation proposal is unavailable"))?;
        if !existing.insert(
            feature_key(&proposal.feature)
                .ok_or_else(|| invalid("proposal has no exact feature geometry"))?,
        ) {
            return Err(invalid(
                "selected feature already exists at this location/orientation; existing annotations were not overwritten",
            ));
        }
        result.features.push(proposal.feature.clone());
    }
    result.metadata.insert(
        "feature_library:source_revision".into(),
        report.source_revision.clone(),
    );
    result.metadata.insert(
        "feature_library:library_revision".into(),
        report.library_revision.clone(),
    );
    result.metadata.insert(
        "feature_library:scan_revision".into(),
        digest(report, b"collagene-feature-report-v1\0", MAX_REPORT_BYTES)?,
    );
    result.metadata.insert(
        "feature_library:proposal_ids".into(),
        serde_json::to_string(selected).map_err(|e| invalid(e.to_string()))?,
    );
    result.validate()?;
    crate::protocol::save_document(&result)?;
    Ok(result)
}
