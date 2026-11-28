#![doc = include_str!("../README.md")]
//! Immutable, versioned restriction enzyme geometry and documentary evidence.
//! This crate neither matches DNA nor validates reactions. Documentary facts are
//! formulation- and assay-context-specific; missing facts never imply compatibility.
pub mod hosts;

use serde::{Deserialize, Serialize};
use std::collections::{BTreeMap, BTreeSet};
use std::io::Write;

pub const SCHEMA_VERSION: u32 = 3;
pub const MAX_JSON_BYTES: usize = 2 * 1024 * 1024;
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Error(pub String);
impl std::fmt::Display for Error {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        self.0.fmt(f)
    }
}
impl std::error::Error for Error {}
fn fail<T>(s: &str) -> Result<T, Error> {
    Err(Error(s.into()))
}
macro_rules! dto {
    ($name:ident { $($(#[$attr:meta])* $field:ident : $ty:ty),* $(,)? }) => {
        #[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
        #[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
        #[serde(deny_unknown_fields)]
        pub struct $name {
            $($(#[$attr])* pub $field: $ty),*
        }
    };
}
dto!(Citation {
    source: String,
    locator: String
});
dto!(Snapshot {
    sha256: String,
    bytes: u64
});
dto!(Source {
    id: String,
    publisher: String,
    url: String,
    retrieved_on: String,
    edition: String,
    snapshot: Option<Snapshot>,
    snapshot_note: String,
    facts_blake3: String
});
dto!(Alias {
    name: String,
    shared: bool,
    citations: Vec<Citation>
});
// Boundaries are between bases, relative to the forward recognition start.
// The bottom boundary uses the same forward coordinate axis; Type IIS may lie outside it.
dto!(Geometry {
    recognition: String,
    top_cut: i32,
    bottom_cut: i32
});
dto!(Enzyme {
    id: String,
    name: String,
    geometry: Geometry,
    aliases: Vec<Alias>,
    citations: Vec<Citation>
});
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(rename_all = "snake_case", deny_unknown_fields)]
pub enum Condition {
    IncubationTemperature,
    BufferActivity {
        buffer: String,
        temperature_milli_c: i32,
    },
    Methylation {
        system: MethylationSystem,
    },
}
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(rename_all = "snake_case")]
pub enum MethylationSystem {
    Dam,
    Dcm,
    CpG,
}
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(rename_all = "snake_case")]
pub enum Sensitivity {
    NotSensitive,
    Blocked,
    Impaired,
    BlockedOverlapping,
    ImpairedOverlapping,
    BlockedSomeOverlappingCombinations,
    ImpairedSomeOverlappingCombinations,
}
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(rename_all = "snake_case", deny_unknown_fields)]
pub enum FactValue {
    TemperatureMilliC(i32),
    Percent(u8),
    LessThanPercent(u8),
    Methylation(Sensitivity),
    Unknown(String),
}
dto!(Evidence {
    condition: Condition,
    context: String,
    value: FactValue,
    qualifiers: Vec<String>,
    citations: Vec<Citation>
});
dto!(Formulation {
    id: String,
    name: String,
    enzyme: String,
    vendor: String,
    product_code: String,
    formulation_scope: String,
    aliases: Vec<Alias>,
    evidence: Vec<Evidence>,
    citations: Vec<Citation>
});
dto!(Edition {
    schema_version: u32,
    edition: String,
    curated_on: String,
    sources: Vec<Source>,
    enzymes: Vec<Enzyme>,
    formulations: Vec<Formulation>,
    #[serde(deserialize_with = "hosts::deserialize_hosts")]
    hosts: Vec<hosts::HostPreset>
});
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Lookup<'a, T> {
    Missing,
    Unique(&'a T),
    Ambiguous(Vec<&'a T>),
}
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Relation {
    SameIdentity,
    SameCutIsoschizomers,
    Neoschizomers,
    DifferentRecognition,
}
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Resolution<'a> {
    Unknown(Vec<&'a Evidence>),
    ContextRequired(Vec<&'a str>),
    Known {
        value: &'a FactValue,
        evidence: Vec<&'a Evidence>,
    },
    Conflict(Vec<&'a Evidence>),
}
/// A validated snapshot. No API mutates an installed edition.
#[derive(Debug, Clone)]
pub struct Catalog {
    edition: Edition,
    content_hash: String,
}
fn norm(s: &str) -> String {
    s.trim().to_ascii_lowercase()
}
fn bounded(s: &str, n: usize) -> bool {
    !s.trim().is_empty() && s.len() <= n && !s.chars().any(char::is_control)
}
fn hash(v: &impl Serialize) -> Result<String, Error> {
    Ok(
        blake3::hash(&serde_json::to_vec(v).map_err(|e| Error(e.to_string()))?)
            .to_hex()
            .to_string(),
    )
}
fn hex(s: &str) -> bool {
    s.len() == 64
        && s.bytes()
            .all(|c| c.is_ascii_hexdigit() && !c.is_ascii_uppercase())
}
struct ByteBudget(usize);
impl Write for ByteBudget {
    fn write(&mut self, bytes: &[u8]) -> std::io::Result<usize> {
        self.0 = self
            .0
            .checked_sub(bytes.len())
            .ok_or_else(|| std::io::Error::other("edition exceeds serialized byte limit"))?;
        Ok(bytes.len())
    }
    fn flush(&mut self) -> std::io::Result<()> {
        Ok(())
    }
}
/// Measure concrete proposal serialization without allocating an output buffer.
/// Native DTO callers receive the same byte bound as imported JSON before any
/// normalization, canonicalization, curated-assertion hashing or edition clone.
fn serialized_preflight(e: &Edition) -> Result<(), Error> {
    serde_json::to_writer(ByteBudget(MAX_JSON_BYTES), e)
        .map_err(|_| Error("edition exceeds serialized byte limit".into()))
}
fn canonicalize(e: &mut Edition) {
    e.sources.sort_by(|a, b| a.id.cmp(&b.id));
    e.enzymes.sort_by(|a, b| a.id.cmp(&b.id));
    e.formulations.sort_by(|a, b| a.id.cmp(&b.id));
    e.hosts.sort_by(|a, b| a.id.cmp(&b.id));
    for host in &mut e.hosts {
        host.citations
            .sort_by_key(|c| (c.source.clone(), c.locator.clone()));
    }
    for enzyme in &mut e.enzymes {
        enzyme.aliases.sort_by_key(|a| norm(&a.name));
        enzyme
            .citations
            .sort_by_key(|c| (c.source.clone(), c.locator.clone()));
    }
    for f in &mut e.formulations {
        f.aliases.sort_by_key(|a| norm(&a.name));
        f.evidence
            .sort_by_key(|x| serde_json::to_string(x).unwrap_or_default());
        f.citations
            .sort_by_key(|c| (c.source.clone(), c.locator.clone()));
    }
}
fn source_facts(e: &Edition, source: &str) -> Result<String, Error> {
    let mut assertions = Vec::new();
    for en in &e.enzymes {
        if en.citations.iter().any(|c| c.source == source) {
            assertions.push(serde_json::json!({"enzyme":en.id,"name":en.name,"geometry":en.geometry,"citations":en.citations.iter().filter(|c|c.source==source).collect::<Vec<_>>()}));
        }
        for a in &en.aliases {
            if a.citations.iter().any(|c| c.source == source) {
                assertions.push(serde_json::json!({"enzyme":en.id,"alias":a}));
            }
        }
    }
    for f in &e.formulations {
        if f.citations.iter().any(|c| c.source == source) {
            assertions.push(serde_json::json!({"formulation":f.id,"name":f.name,"enzyme":f.enzyme,"vendor":f.vendor,"product_code":f.product_code,"scope":f.formulation_scope,"citations":f.citations.iter().filter(|c|c.source==source).collect::<Vec<_>>()}));
        }
        for a in &f.aliases {
            if a.citations.iter().any(|c| c.source == source) {
                assertions.push(serde_json::json!({"formulation":f.id,"alias":a}));
            }
        }
        for v in &f.evidence {
            if v.citations.iter().any(|c| c.source == source) {
                assertions.push(serde_json::json!({"formulation":f.id,"evidence":v}));
            }
        }
    }
    for host in &e.hosts {
        if host.citations.iter().any(|c| c.source == source) {
            assertions.push(serde_json::json!({"host":host.id,"record":host}));
        }
    }
    assertions.sort_by_key(|v| v.to_string());
    hash(&assertions)
}
/// Seal independently curated proposals. This computes checksums, not authenticity.
/// Untrusted imported editions must use `Catalog::from_json`, which verifies them.
pub fn seal(mut edition: Edition) -> Result<Edition, Error> {
    capacity(&edition)?;
    serialized_preflight(&edition)?;
    for source in &mut edition.sources {
        source.facts_blake3 = "0".repeat(64)
    }
    // Empty authoring checksums can grow during normalization. Bound their final
    // fixed-width representation before hashing, too.
    serialized_preflight(&edition)?;
    validate(&edition)?;
    canonicalize(&mut edition);
    for i in 0..edition.sources.len() {
        edition.sources[i].facts_blake3 = source_facts(&edition, &edition.sources[i].id)?;
    }
    Catalog::from_edition(edition.clone())?;
    Ok(edition)
}
impl Catalog {
    pub fn builtin() -> Result<Self, Error> {
        Self::from_json(include_bytes!("../data/catalog.json"))
    }
    pub fn from_json(bytes: &[u8]) -> Result<Self, Error> {
        if bytes.len() > MAX_JSON_BYTES {
            return fail("JSON exceeds 2 MiB limit");
        };
        Self::from_edition(serde_json::from_slice(bytes).map_err(|e| Error(e.to_string()))?)
    }
    pub fn from_edition(mut e: Edition) -> Result<Self, Error> {
        capacity(&e)?;
        serialized_preflight(&e)?;
        validate(&e)?;
        canonicalize(&mut e);
        for s in &e.sources {
            if source_facts(&e, &s.id)? != s.facts_blake3 {
                return fail("source curated-facts checksum mismatch");
            }
        }
        let bytes = serde_json::to_vec(&e).map_err(|x| Error(x.to_string()))?;
        if bytes.len() > MAX_JSON_BYTES {
            return fail("edition exceeds serialized limit");
        };
        let content_hash = blake3::hash(&bytes).to_hex().to_string();
        Ok(Self {
            edition: e,
            content_hash,
        })
    }
    pub fn edition(&self) -> &Edition {
        &self.edition
    }
    pub fn content_hash(&self) -> &str {
        &self.content_hash
    }
    pub fn to_json(&self) -> Result<Vec<u8>, Error> {
        serde_json::to_vec(&self.edition).map_err(|e| Error(e.to_string()))
    }
    pub fn lookup_enzyme(&self, name: &str) -> Lookup<'_, Enzyme> {
        lookup(&self.edition.enzymes, name, |e| {
            (&e.id, &e.name, &e.aliases)
        })
    }
    pub fn lookup_formulation(&self, name: &str) -> Lookup<'_, Formulation> {
        lookup(&self.edition.formulations, name, |f| {
            (&f.id, &f.name, &f.aliases)
        })
    }
    pub fn search(&self, query: &str, limit: usize) -> Result<Vec<&Enzyme>, Error> {
        if query.len() > 128 || limit > 1000 {
            return fail("search limit exceeded");
        };
        let q = norm(query);
        Ok(self
            .edition
            .enzymes
            .iter()
            .filter(|e| {
                norm(&e.name).contains(&q)
                    || norm(&e.id).contains(&q)
                    || e.aliases.iter().any(|a| norm(&a.name).contains(&q))
            })
            .take(limit)
            .collect())
    }
    /// Compare only exact records installed in this sealed catalog. Value-identical
    /// clones are accepted; an ID alone never establishes identity or geometry.
    pub fn relation(&self, a: &Enzyme, b: &Enzyme) -> Result<Relation, Error> {
        for record in [a, b] {
            if !self
                .edition
                .enzymes
                .iter()
                .any(|installed| installed == record)
            {
                return fail("enzyme record does not belong to this sealed catalog");
            }
        }
        if a.id == b.id {
            return Ok(Relation::SameIdentity);
        }
        let ga = canonical_geometry(&a.geometry);
        let gb = canonical_geometry(&b.geometry);
        Ok(if ga.0 != gb.0 {
            Relation::DifferentRecognition
        } else if ga == gb {
            Relation::SameCutIsoschizomers
        } else {
            Relation::Neoschizomers
        })
    }
    pub fn resolve<'a>(
        &'a self,
        formulation: &str,
        condition: &Condition,
        context: Option<&str>,
    ) -> Result<Resolution<'a>, Error> {
        let f = match self.lookup_formulation(formulation) {
            Lookup::Unique(f) => f,
            Lookup::Missing => return fail("unknown formulation"),
            Lookup::Ambiguous(_) => return fail("ambiguous formulation"),
        };
        let mut evidence: Vec<_> = f
            .evidence
            .iter()
            .filter(|e| &e.condition == condition)
            .collect();
        let contexts: BTreeSet<_> = evidence.iter().map(|e| e.context.as_str()).collect();
        if context.is_none() && contexts.len() > 1 {
            return Ok(Resolution::ContextRequired(contexts.into_iter().collect()));
        };
        if let Some(c) = context {
            evidence.retain(|e| e.context == c)
        }
        let known: Vec<_> = evidence
            .iter()
            .filter(|e| !matches!(e.value, FactValue::Unknown(_)))
            .collect();
        if known.is_empty() {
            return Ok(Resolution::Unknown(evidence));
        };
        let value = &known[0].value;
        if known.iter().any(|e| &e.value != value) {
            Ok(Resolution::Conflict(evidence))
        } else {
            Ok(Resolution::Known { value, evidence })
        }
    }
}
fn lookup<'a, T>(
    items: &'a [T],
    query: &str,
    fields: impl Fn(&'a T) -> (&'a str, &'a str, &'a [Alias]),
) -> Lookup<'a, T> {
    if query.len() > 128 {
        return Lookup::Missing;
    };
    let q = norm(query);
    let found: Vec<_> = items
        .iter()
        .filter(|x| {
            let (id, name, aliases) = fields(x);
            norm(id) == q || norm(name) == q || aliases.iter().any(|a| norm(&a.name) == q)
        })
        .collect();
    match found.len() {
        0 => Lookup::Missing,
        1 => Lookup::Unique(found[0]),
        _ => Lookup::Ambiguous(found),
    }
}
fn canonical_geometry(g: &Geometry) -> (String, i32, i32) {
    let rev: String = g
        .recognition
        .chars()
        .rev()
        .map(|c| match c {
            'A' => 'T',
            'T' => 'A',
            'C' => 'G',
            'G' => 'C',
            'R' => 'Y',
            'Y' => 'R',
            'S' => 'S',
            'W' => 'W',
            'K' => 'M',
            'M' => 'K',
            'B' => 'V',
            'V' => 'B',
            'D' => 'H',
            'H' => 'D',
            _ => 'N',
        })
        .collect();
    let n = g.recognition.len() as i32;
    let f = (g.recognition.clone(), g.top_cut, g.bottom_cut);
    let r = (rev, n - g.bottom_cut, n - g.top_cut);
    std::cmp::min(f, r)
}
fn citations(cs: &[Citation], sources: &BTreeSet<&str>) -> Result<(), Error> {
    if cs.is_empty() || cs.len() > 16 {
        return fail("missing or excessive citations");
    };
    for c in cs {
        if !sources.contains(c.source.as_str()) || !bounded(&c.locator, 512) {
            return fail("invalid citation");
        }
    }
    Ok(())
}
fn names<'a>(
    rows: impl Iterator<Item = (&'a str, &'a str, &'a [Alias])>,
    sources: &BTreeSet<&str>,
) -> Result<(), Error> {
    let rows: Vec<_> = rows.collect();
    let mut canonical = BTreeMap::new();
    let mut alias_map: BTreeMap<String, Vec<(&str, bool)>> = BTreeMap::new();
    for (id, name, _) in &rows {
        if !bounded(id, 96) || !bounded(name, 128) {
            return fail("invalid identity");
        };
        for key in [norm(id), norm(name)] {
            if let Some(old) = canonical.insert(key, *id) {
                if old != *id {
                    return fail("colliding identity");
                }
            }
        }
    }
    for (id, _, aliases) in rows {
        if aliases.len() > 16 {
            return fail("too many aliases");
        };
        let mut seen = BTreeSet::new();
        for a in aliases {
            if !bounded(&a.name, 128) {
                return fail("invalid alias name");
            }
            let key = norm(&a.name);
            if !seen.insert(key.clone()) || canonical.contains_key(&key) {
                return fail("duplicate or shadowing alias");
            };
            citations(&a.citations, sources)?;
            alias_map.entry(key).or_default().push((id, a.shared));
        }
    }
    for aliases in alias_map.values() {
        if aliases.len() > 1 && aliases.iter().any(|a| !a.1) {
            return fail("unacknowledged alias collision");
        }
    }
    Ok(())
}
fn validate(e: &Edition) -> Result<(), Error> {
    if e.schema_version != SCHEMA_VERSION || !bounded(&e.edition, 128) || !date(&e.curated_on) {
        return fail("invalid edition/schema");
    };
    if e.sources.is_empty()
        || e.sources.len() > 128
        || e.enzymes.is_empty()
        || e.enzymes.len() > 1000
        || e.formulations.len() > 2000
        || e.hosts.len() > hosts::MAX_HOSTS
    {
        return fail("edition capacity exceeded");
    }
    let mut source_ids = BTreeSet::new();
    for s in &e.sources {
        if !bounded(&s.id, 96)
            || !source_ids.insert(s.id.as_str())
            || !bounded(&s.publisher, 128)
            || !s.url.starts_with("https://")
            || !bounded(&s.url, 1024)
            || s.url.contains('@')
            || !date(&s.retrieved_on)
            || !bounded(&s.edition, 256)
            || !hex(&s.facts_blake3)
            || !bounded(&s.snapshot_note, 512)
        {
            return fail("invalid source");
        };
        if let Some(snapshot) = &s.snapshot {
            if !hex(&snapshot.sha256) || snapshot.bytes == 0 {
                return fail("invalid snapshot");
            }
        }
    }
    names(
        e.enzymes
            .iter()
            .map(|x| (x.id.as_str(), x.name.as_str(), x.aliases.as_slice())),
        &source_ids,
    )?;
    names(
        e.formulations
            .iter()
            .map(|x| (x.id.as_str(), x.name.as_str(), x.aliases.as_slice())),
        &source_ids,
    )?;
    let mut host_ids = BTreeSet::new();
    for host in &e.hosts {
        host.validate()?;
        if !host_ids.insert(host.id.as_str()) {
            return fail("duplicate source-host identity");
        }
        citations(&host.citations, &source_ids)?;
    }
    let mut ids = BTreeSet::new();
    for en in &e.enzymes {
        if !ids.insert(en.id.as_str()) {
            return fail("duplicate enzyme id");
        };
        let g = &en.geometry;
        if g.recognition.is_empty()
            || g.recognition.len() > 64
            || !g
                .recognition
                .bytes()
                .all(|c| b"ACGTRYSWKMBDHVN".contains(&c))
            || !(-1024..=1024).contains(&g.top_cut)
            || !(-1024..=1024).contains(&g.bottom_cut)
        {
            return fail("invalid recognition geometry");
        };
        citations(&en.citations, &source_ids)?;
    }
    let mut formids = BTreeSet::new();
    for f in &e.formulations {
        if !formids.insert(&f.id)
            || !ids.contains(f.enzyme.as_str())
            || !bounded(&f.vendor, 128)
            || !bounded(&f.product_code, 96)
            || !bounded(&f.formulation_scope, 1024)
            || f.evidence.len() > 64
        {
            return fail("invalid formulation");
        };
        citations(&f.citations, &source_ids)?;
        for v in &f.evidence {
            if !bounded(&v.context, 256)
                || v.qualifiers.len() > 16
                || v.qualifiers.iter().any(|s| !bounded(s, 1024))
            {
                return fail("invalid evidence context");
            };
            citations(&v.citations, &source_ids)?;
            if let Condition::BufferActivity {
                buffer,
                temperature_milli_c,
            } = &v.condition
            {
                if !bounded(buffer, 128) || !(0..=100000).contains(temperature_milli_c) {
                    return fail("invalid buffer condition");
                }
            };
            let valid = match (&v.condition, &v.value) {
                (_, FactValue::Unknown(s)) => bounded(s, 1024),
                (Condition::IncubationTemperature, FactValue::TemperatureMilliC(t)) => {
                    (0..=100000).contains(t)
                }
                (
                    Condition::BufferActivity {
                        buffer,
                        temperature_milli_c,
                    },
                    FactValue::Percent(p),
                ) => {
                    bounded(buffer, 128) && (0..=100000).contains(temperature_milli_c) && *p <= 100
                }
                (
                    Condition::BufferActivity {
                        buffer,
                        temperature_milli_c,
                    },
                    FactValue::LessThanPercent(p),
                ) => {
                    bounded(buffer, 128)
                        && (0..=100000).contains(temperature_milli_c)
                        && (1..=100).contains(p)
                }
                (Condition::Methylation { .. }, FactValue::Methylation(_)) => true,
                _ => false,
            };
            if !valid {
                return fail("condition/value mismatch");
            }
        }
    }
    Ok(())
}
fn date(s: &str) -> bool {
    let b = s.as_bytes();
    if b.len() != 10
        || b[4] != b'-'
        || b[7] != b'-'
        || !b
            .iter()
            .enumerate()
            .all(|(i, c)| i == 4 || i == 7 || c.is_ascii_digit())
    {
        return false;
    }
    // Exact ASCII YYYY-MM-DD, proleptic Gregorian years 0001–9999.
    // No year zero, BCE, reduced-width dates or timezone/calendar inference.
    let Ok(year) = s[..4].parse::<u16>() else {
        return false;
    };
    let Ok(month) = s[5..7].parse::<u8>() else {
        return false;
    };
    let Ok(day) = s[8..].parse::<u8>() else {
        return false;
    };
    if year == 0 {
        return false;
    }
    let leap = year % 400 == 0 || year % 4 == 0 && year % 100 != 0;
    let days = match month {
        1 | 3 | 5 | 7 | 8 | 10 | 12 => 31,
        4 | 6 | 9 | 11 => 30,
        2 if leap => 29,
        2 => 28,
        _ => return false,
    };
    (1..=days).contains(&day)
}

fn capacity(e: &Edition) -> Result<(), Error> {
    if e.sources.len() > 128
        || e.enzymes.len() > 1000
        || e.formulations.len() > 2000
        || e.hosts.len() > hosts::MAX_HOSTS
        || e.enzymes
            .iter()
            .any(|x| x.aliases.len() > 16 || x.citations.len() > 16)
        || e.formulations
            .iter()
            .any(|x| x.aliases.len() > 16 || x.evidence.len() > 64 || x.citations.len() > 16)
        || e.hosts
            .iter()
            .any(|h| h.qualifiers.len() > 16 || h.citations.len() > 16)
    {
        return fail("edition capacity exceeded");
    };
    Ok(())
}
