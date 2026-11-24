//! Reuse declared primer sequences from actual current documents. This grants
//! neither physical oligo state nor fresh annealing/amplification evidence.
use crate::{Alphabet, Error, Feature, Location, Molecule, Sequence, Topology, analysis::Primer};
use serde::{Deserialize, Serialize};
use std::io::Write;
pub const MAX_ORDERED_BASES: usize = 200;
pub const MAX_LIST: usize = 1000;
const MAX_BYTES: usize = 32 << 20;
macro_rules! dto{($($t:tt)*)=>{#[derive(Clone,Debug,Serialize,Deserialize)]#[cfg_attr(feature="json-schema",derive(schemars::JsonSchema))]#[serde(deny_unknown_fields)]$($t)*};}
dto! {pub enum SequenceOrigin {StoredDomains,OrderedWithDeclaredTail,CurrentLocation}}
dto! {pub struct Resolved {pub source_revision:String,pub feature_index:usize,pub primer:Primer,pub sequence_origin:SequenceOrigin,pub warnings:Vec<String>}}
dto! {#[derive(PartialEq,Eq)] pub struct Summary {pub feature_index:usize,pub name:String,pub reusable:bool,pub diagnostic:Option<String>,pub annealing_bases:Option<usize>,pub tail_bases:Option<usize>,pub ordered_bases:Option<usize>,pub needs_review:bool}}
dto! {pub struct List {pub source_revision:String,pub total:usize,pub omitted:usize,pub primers:Vec<Summary>}}
fn bad(s: &str) -> Error {
    Error::Unsupported(format!("primer reuse: {s}"))
}
struct Hash {
    hash: blake3::Hasher,
    bytes: usize,
}
impl Write for Hash {
    fn write(&mut self, b: &[u8]) -> std::io::Result<usize> {
        self.bytes = self
            .bytes
            .checked_add(b.len())
            .filter(|n| *n <= MAX_BYTES)
            .ok_or_else(|| std::io::Error::other("bounded primer source"))?;
        self.hash.update(b);
        Ok(b.len())
    }
    fn flush(&mut self) -> std::io::Result<()> {
        Ok(())
    }
}
fn revision(m: &Molecule) -> Result<String, Error> {
    let mut h = Hash {
        hash: blake3::Hasher::new(),
        bytes: 0,
    };
    serde_json::to_writer(&mut h, m).map_err(|_| Error::TooLarge)?;
    Ok(h.hash.finalize().to_hex().to_string())
}
fn hex(s: &str) -> bool {
    s.len() == 64
        && s.bytes()
            .all(|b| b.is_ascii_digit() || (b'a'..=b'f').contains(&b))
}
/// Borrowed aggregate metadata/location guards precede hashing or extraction.
pub fn preflight(m: &Molecule) -> Result<(), Error> {
    crate::model::bounds::model_preflight(&[m])?;
    if m.sequence.alphabet() != Alphabet::Dna {
        return Err(bad(
            "primer reuse requires a DNA document; RNA/protein conversion is not inferred",
        ));
    }
    m.validate()
}
fn value<'a>(f: &'a Feature, key: &str) -> Result<Option<&'a str>, Error> {
    let mut fields = f
        .qualifiers
        .iter()
        .filter(|(k, _)| k.eq_ignore_ascii_case(key));
    let first = fields.next();
    if fields.next().is_some() {
        return Err(bad(
            "duplicate primer sequence/name/source qualifier; resolve ambiguity explicitly",
        ));
    }
    match first {
        None => Ok(None),
        Some((_, Some(value))) if value.len() <= 8192 => Ok(Some(value)),
        Some((_, Some(_))) => Err(Error::TooLarge),
        _ => Err(bad("primer qualifier lacks an explicit value")),
    }
}
fn name(f: &Feature) -> Result<String, Error> {
    let label = value(f, "label")?;
    let explicit = value(f, "primer_name")?;
    if label.zip(explicit).is_some_and(|(a, b)| a != b) {
        return Err(bad("conflicting explicit primer names"));
    }
    let name = label.or(explicit).unwrap_or_else(|| f.label());
    if name.trim().is_empty() || name.len() > 256 || name.chars().any(char::is_control) {
        return Err(bad("primer name must be nonempty bounded text"));
    }
    Ok(name.into())
}
fn dna(raw: &str, empty: bool) -> Result<String, Error> {
    if raw.len() > 8192 {
        return Err(Error::TooLarge);
    }
    if empty && raw.trim().is_empty() {
        return Ok(String::new());
    }
    let sequence = Sequence::dna(raw)?;
    if sequence.len() > MAX_ORDERED_BASES {
        return Err(bad(
            "stored sequence exceeds the 200 ordered-base reuse profile",
        ));
    }
    Ok(sequence.as_str().into())
}
/// Require one contiguous forward-reference interval, possibly with one explicit
/// circular origin crossing, under a uniform outer complement. No gaps, fuzzy
/// boundaries, repeated bases, mixed strands or unresolved/order sites are guessed.
fn current_domain(f: &Feature, m: &Molecule) -> Result<String, Error> {
    if f.needs_review {
        return Err(bad(
            "review-marked coordinates cannot supply missing primer sequence fields",
        ));
    }
    let mut root = &f.location;
    while let Location::Complement(next) = root {
        root = next;
    }
    let mut stack = vec![root];
    let mut ranges = Vec::new();
    let mut total = 0usize;
    while let Some(location) = stack.pop() {
        match location {
            Location::Range {
                start,
                end,
                before: false,
                after: false,
            } => {
                if start >= end || *end > m.sequence.len() {
                    return Err(bad("primer location is not an exact nonempty interval"));
                }
                total = total
                    .checked_add(end - start)
                    .filter(|n| *n <= MAX_ORDERED_BASES && *n <= m.sequence.len())
                    .ok_or_else(|| bad("primer location exceeds one turn or 200 bases"))?;
                ranges.push((*start, *end));
            }
            Location::Join(parts) => {
                for part in parts.iter().rev() {
                    stack.push(part)
                }
            }
            _ => {
                return Err(bad(
                    "missing sequence fields require exact single-strand contiguous primer geometry",
                ));
            }
        }
    }
    if ranges.is_empty() {
        return Err(bad("empty primer location"));
    }
    let mut wraps = 0;
    for pair in ranges.windows(2) {
        if pair[0].1 == pair[1].0 {
            continue;
        }
        if m.topology == Topology::Circular && pair[0].1 == m.sequence.len() && pair[1].0 == 0 {
            wraps += 1;
        } else {
            return Err(bad(
                "gapped, overlapping or disordered primer coordinates cannot supply missing sequence",
            ));
        }
    }
    if wraps > 1 {
        return Err(bad("primer location crosses the origin more than once"));
    }
    f.location.extract(&m.sequence)
}
fn fields(f: &Feature, m: &Molecule) -> Result<(Primer, SequenceOrigin, Vec<String>), Error> {
    if !f.kind.eq_ignore_ascii_case("primer_bind") {
        return Err(bad(
            "select an explicit primer_bind feature; bare DNA does not define oligo domains",
        ));
    }
    let name = name(f)?;
    let raw_annealing = value(f, "annealing_sequence")?;
    let raw_ordered = value(f, "oligo_sequence")?;
    let raw_tail = value(f, "5prime_tail")?;
    let historical = value(f, "source_revision")?;
    if historical.is_some_and(|r| !hex(r)) {
        return Err(bad("malformed historical primer source revision"));
    }
    let tail = raw_tail
        .map(|s| dna(s, true))
        .transpose()?
        .unwrap_or_default();
    let ordered = raw_ordered.map(|s| dna(s, false)).transpose()?;
    let mut warnings=vec!["Sequence reuse only: fresh binding, specificity, physical oligo chemistry and amplification are not established.".into()];
    if historical.is_some() {
        warnings.push("Stored source_revision is historical binding evidence, not current coordinate or primer authority. The actual current document revision is bound separately.".into());
    }
    if f.needs_review {
        warnings.push("The source feature is marked for review; stored sequence fields are reused without asserting that its location is a current binding site.".into());
    }
    let (annealing, origin) = if let Some(raw) = raw_annealing {
        (dna(raw, false)?, SequenceOrigin::StoredDomains)
    } else if let (Some(ordered), Some(_)) = (&ordered, raw_tail) {
        let annealing = ordered.strip_prefix(&tail).ok_or_else(|| {
            bad("stored ordered oligo does not begin with its declared 5prime tail")
        })?;
        if annealing.is_empty() {
            return Err(bad(
                "a declared tail cannot consume the entire ordered oligo",
            ));
        }
        (annealing.into(), SequenceOrigin::OrderedWithDeclaredTail)
    } else {
        let domain = current_domain(f, m)?;
        warnings.push("Missing annealing sequence was extracted from exact current coordinates; no independent binding or historical oligo identity is inferred.".into());
        (domain, SequenceOrigin::CurrentLocation)
    };
    let primer = Primer::new(name, &annealing, &tail)?;
    if primer.annealing.len() + primer.tail.len() > MAX_ORDERED_BASES {
        return Err(bad("primer exceeds the 200 ordered-base reuse profile"));
    }
    if ordered
        .as_ref()
        .is_some_and(|s| s != &primer.ordered_sequence())
    {
        return Err(bad(
            "stored annealing/tail fields disagree with the ordered oligo; a missing tail split is not invented",
        ));
    }
    if matches!(
        origin,
        SequenceOrigin::StoredDomains | SequenceOrigin::OrderedWithDeclaredTail
    ) {
        // Compare only when current geometry is defensible, but never replace a stored
        // mismatch primer with the template's calls. Failure is visible evidence.
        match current_domain(f,m){Ok(domain)if domain!=primer.annealing.as_str()=>warnings.push("Stored annealing sequence differs from the current annotated reference calls. It is retained as declared (possibly intentional mismatch/stale evidence); recheck binding before design use.".into()),Err(_)=>warnings.push("Current feature geometry cannot independently recover this primer domain. Explicit stored fields were used; coordinate/binding validity remains unresolved.".into()),_=>{}}
    }
    if primer
        .ordered_sequence()
        .bytes()
        .any(|b| !b"ACGT".contains(&b))
    {
        warnings.push("IUPAC ambiguity is preserved as sequence notation. Exact-PCR consumers may reject it; no concrete binding or mixture is inferred.".into());
    }
    Ok((primer, origin, warnings))
}
pub fn resolve(m: &Molecule, feature_index: usize) -> Result<Resolved, Error> {
    preflight(m)?;
    let feature = m
        .features
        .get(feature_index)
        .ok_or_else(|| bad("feature index is unavailable"))?;
    let (primer, sequence_origin, warnings) = fields(feature, m)?;
    Ok(Resolved {
        source_revision: revision(m)?,
        feature_index,
        primer,
        sequence_origin,
        warnings,
    })
}
pub fn resolve_at_revision(
    m: &Molecule,
    feature_index: usize,
    expected_revision: &str,
) -> Result<Resolved, Error> {
    preflight(m)?;
    if !hex(expected_revision) {
        return Err(bad("invalid expected document revision"));
    }
    let source_revision = revision(m)?;
    if source_revision != expected_revision {
        return Err(bad("document revision changed; select the primer again"));
    }
    let feature = m
        .features
        .get(feature_index)
        .ok_or_else(|| bad("feature index is unavailable"))?;
    let (primer, sequence_origin, warnings) = fields(feature, m)?;
    Ok(Resolved {
        source_revision,
        feature_index,
        primer,
        sequence_origin,
        warnings,
    })
}
/// Bounded summaries retain actual feature indices, including unusable entries.
/// One source revision is computed per call, never once per feature.
pub fn list(m: &Molecule, limit: usize) -> Result<List, Error> {
    preflight(m)?;
    if limit == 0 || limit > MAX_LIST {
        return Err(Error::TooLarge);
    }
    let source_revision = revision(m)?;
    let total = m
        .features
        .iter()
        .filter(|f| f.kind.eq_ignore_ascii_case("primer_bind"))
        .count();
    let mut primers = Vec::with_capacity(limit.min(total));
    for (feature_index, f) in m
        .features
        .iter()
        .enumerate()
        .filter(|(_, f)| f.kind.eq_ignore_ascii_case("primer_bind"))
        .take(limit)
    {
        let result = fields(f, m);
        let name = name(f).unwrap_or_else(|_| format!("primer feature {}", feature_index + 1));
        let (reusable, diagnostic, annealing_bases, tail_bases, ordered_bases) = match result {
            Ok((p, _, _)) => (
                true,
                None,
                Some(p.annealing.len()),
                Some(p.tail.len()),
                Some(p.annealing.len() + p.tail.len()),
            ),
            Err(error) => (false, Some(error.to_string()), None, None, None),
        };
        primers.push(Summary {
            feature_index,
            name,
            reusable,
            diagnostic,
            annealing_bases,
            tail_bases,
            ordered_bases,
            needs_review: f.needs_review,
        });
    }
    Ok(List {
        source_revision,
        total,
        omitted: total - primers.len(),
        primers,
    })
}
