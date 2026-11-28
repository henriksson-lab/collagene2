//! Explicit reference edit intent, independent of every molecular preparation engine.
use crate::{annotations, assembly_recipe};
use collagene_model::{Alphabet, Error, Molecule, Sequence, Topology};
use serde::{Deserialize, Serialize};
use std::io::Write;
pub const ENGINE_VERSION: &str = "mutation-intent-1";
const MAX_SOURCE: usize = 8 << 20;
const MAX_INSERTED: usize = 100_000;
const REFERENCE_MODEL: &str =
    "Reference-only mutation design intent; physical material state unavailable";
macro_rules! dto {($($t:tt)*)=>{#[derive(Clone,Debug,Serialize,Deserialize)]#[cfg_attr(feature="json-schema",derive(schemars::JsonSchema))]#[serde(deny_unknown_fields)]$($t)*};}
dto! {pub struct Request {pub name:String,pub edit_start:usize,pub edit_end:usize,pub inserted:String}}
dto! {pub struct Report {pub context_revision:String,pub source_revision:String,pub request:Request,pub target:Molecule,pub deleted:String,pub origin_shift:usize,pub warnings:Vec<String>}}
fn bad(s: &str) -> Error {
    Error::Unsupported(format!("mutation intent: {s}"))
}
fn hex(s: &str) -> bool {
    s.len() == 64
        && s.bytes()
            .all(|b| b.is_ascii_digit() || (b'a'..=b'f').contains(&b))
}
fn request_guard(r: &Request) -> Result<(), Error> {
    if r.name.trim().is_empty() || r.name.len() > 1024 || r.name.chars().any(char::is_control) {
        return Err(bad("provide a bounded explicit name"));
    }
    if r.inserted.len() > MAX_INSERTED * 2 {
        return Err(Error::TooLarge);
    }
    if !r.inserted.is_ascii() {
        return Err(bad(
            "inserted sequence must be DNA with optional ASCII whitespace",
        ));
    }
    Ok(())
}
fn inserted(r: &Request) -> Result<String, Error> {
    request_guard(r)?;
    if r.inserted.bytes().all(|b| b.is_ascii_whitespace()) {
        return Ok(String::new());
    }
    let sequence = Sequence::dna(&r.inserted)?;
    if sequence.len() > MAX_INSERTED {
        return Err(Error::TooLarge);
    }
    Ok(sequence.as_str().into())
}
fn interval(n: usize, topology: Topology, r: &Request) -> Result<usize, Error> {
    let deleted = r
        .edit_end
        .checked_sub(r.edit_start)
        .ok_or_else(|| bad("edit end must not precede start"))?;
    let valid = match topology {
        Topology::Linear => r.edit_start <= n && r.edit_end <= n,
        Topology::Circular => {
            (r.edit_start < n && deleted <= n) || (r.edit_start == n && r.edit_end == n)
        }
    };
    if !valid {
        return Err(bad("edit interval exceeds the source or one circular turn"));
    }
    Ok(deleted)
}
/// Borrowed guards before rich source copying or canonical hashing. Readable
/// physical/reference markers do not turn an intended edit into material state.
pub fn source_preflight(source: &Molecule, request: &Request) -> Result<(), Error> {
    request_guard(request)?;
    if source.sequence.len() > MAX_SOURCE {
        return Err(Error::TooLarge);
    }
    collagene_model::bounds::model_preflight(&[source])?;
    source.validate()?;
    if source.sequence.alphabet() != Alphabet::Dna {
        return Err(bad("mutation reference intent requires DNA"));
    }
    let removed = interval(source.sequence.len(), source.topology, request)?;
    let added = inserted(request)?;
    if source.sequence.len() - removed + added.len() == 0 {
        return Err(bad(
            "an intended reference must retain or introduce at least one base",
        ));
    }
    if source.features.iter().any(|f| {
        f.qualifiers.iter().any(|(k, _)| {
            matches!(
                k.as_str(),
                "transl_except" | "anticodon" | "rpt_unit_range" | "coded_by"
            )
        })
    }) {
        return Err(bad(
            "coordinate-bearing qualifiers cannot yet be safely edited; original source remains intact",
        ));
    }
    assembly_recipe::hash(&(source, request))?;
    Ok(())
}
pub fn intent_revision(source: &Molecule, request: &Request) -> Result<String, Error> {
    source_preflight(source, request)?;
    assembly_recipe::hash(&(ENGINE_VERSION, source.revision(), request))
}
struct EvidenceWriter {
    bytes: Vec<u8>,
}
impl Write for EvidenceWriter {
    fn write(&mut self, b: &[u8]) -> std::io::Result<usize> {
        if self
            .bytes
            .len()
            .checked_add(b.len())
            .is_none_or(|n| n > 16 << 20)
        {
            return Err(std::io::Error::other("source feature evidence cap"));
        }
        self.bytes.extend_from_slice(b);
        Ok(b.len())
    }
    fn flush(&mut self) -> std::io::Result<()> {
        Ok(())
    }
}
fn encode(v: &impl Serialize) -> Result<String, Error> {
    let mut w = EvidenceWriter { bytes: Vec::new() };
    serde_json::to_writer(&mut w, v).map_err(|_| Error::TooLarge)?;
    String::from_utf8(w.bytes).map_err(|_| bad("invalid evidence encoding"))
}
fn retain(target: &mut Molecule, key: String, value: String) -> Result<(), Error> {
    if target.metadata.get(&key).is_some_and(|old| old != &value) {
        return Err(bad("source evidence identity collision"));
    }
    target.metadata.insert(key, value);
    Ok(())
}
pub fn build(source: &Molecule, request: &Request) -> Result<Report, Error> {
    let context_revision = intent_revision(source, request)?;
    let source_revision = source.revision();
    let n = source.sequence.len();
    let count = interval(n, source.topology, request)?;
    let replacement = inserted(request)?;
    let wrap = request.edit_end > n;
    let origin_shift = if wrap { request.edit_start } else { 0 };
    let deleted = if wrap {
        format!(
            "{}{}",
            &source.sequence.as_str()[request.edit_start..],
            &source.sequence.as_str()[..request.edit_end - n]
        )
    } else {
        source.sequence.as_str()[request.edit_start..request.edit_end].into()
    };
    let mut working = source.clone();
    // Transform a readable reference only. Preserve original chemical descriptors
    // as source evidence; the final target always receives a restrictive guard.
    if let Some(model) = working.metadata.remove("duplex_model") {
        retain(
            &mut working,
            format!("mutation_intent:source_duplex_model:{source_revision}"),
            model,
        )?;
    }
    let snapgene = working
        .metadata
        .remove("source:snapgene:chemistry_unmapped");
    if wrap {
        working = annotations::rotate_origin(&working, request.edit_start)?;
    }
    let mut target = if wrap {
        working.replace(0, count, &replacement)?
    } else {
        working.replace(request.edit_start, request.edit_end, &replacement)?
    };
    if let Some(value) = snapgene {
        target
            .metadata
            .insert("source:snapgene:chemistry_unmapped".into(), value);
    }
    if !source.features.is_empty() {
        retain(
            &mut target,
            format!("mutation_intent:source_features:{source_revision}"),
            encode(&source.features)?,
        )?;
    }
    target.name = request.name.clone();
    for (key, value) in [
        ("duplex_model", REFERENCE_MODEL.to_owned()),
        ("mutation:intent:reference_only", "true".into()),
        ("mutation_intent:engine", ENGINE_VERSION.into()),
        ("mutation_intent:source_revision", source_revision.clone()),
        ("mutation_intent:context_revision", context_revision.clone()),
        ("mutation_intent:request", encode(request)?),
        ("mutation_intent:source_length", n.to_string()),
        ("mutation_intent:source_topology", encode(&source.topology)?),
        ("mutation_intent:origin_shift", origin_shift.to_string()),
        (
            "mutation_intent:sequence_blake3",
            blake3::hash(target.sequence.as_str().as_bytes())
                .to_hex()
                .to_string(),
        ),
        (
            "mutation_intent:deleted_blake3",
            blake3::hash(deleted.as_bytes()).to_hex().to_string(),
        ),
    ] {
        target.metadata.insert(key.into(), value);
    }
    let mut warnings=vec!["Intended reference from an explicit source edit only; no PCR, synthesis, repair, cleavage, circularization or experimental success is established.".into(),"Retained source annotations and chemical metadata are source evidence; affected annotations require review.".into()];
    if wrap {
        warnings.push(format!("Circular reference origin moved to source boundary {}; output is inserted bases followed by the retained source path.",request.edit_start));
    }
    if count == n {
        warnings.push("This replaces the complete source turn; an inverse-PCR preparation may not support this intent.".into());
    }
    let report = Report {
        context_revision: context_revision.clone(),
        source_revision,
        request: request.clone(),
        target,
        deleted,
        origin_shift,
        warnings,
    };
    report.validate_binding(&context_revision)?;
    Ok(report)
}
impl Report {
    /// Readonly consistency only. No decoded report can authorize preparation.
    pub fn validate_binding(&self, expected: &str) -> Result<(), Error> {
        if !hex(expected) || self.context_revision != expected || !hex(&self.source_revision) {
            return Err(bad("source/request context changed"));
        }
        request_guard(&self.request)?;
        if assembly_recipe::hash(&(ENGINE_VERSION, &self.source_revision, &self.request))?
            != self.context_revision
        {
            return Err(bad("request does not match its context revision"));
        }
        if self.deleted.len() > MAX_SOURCE
            || self.target.sequence.len() > MAX_SOURCE + MAX_INSERTED
            || self.warnings.len() > 8
            || self.warnings.iter().any(|w| w.len() > 4096)
        {
            return Err(Error::TooLarge);
        }
        collagene_model::bounds::model_preflight(&[&self.target])?;
        self.target.validate()?;
        if !self.deleted.is_empty() && Sequence::dna(&self.deleted)?.as_str() != self.deleted {
            return Err(bad("deleted calls changed"));
        }
        let n = self
            .target
            .metadata
            .get("mutation_intent:source_length")
            .and_then(|s| s.parse::<usize>().ok())
            .filter(|n| *n > 0 && *n <= MAX_SOURCE)
            .ok_or_else(|| bad("missing source length"))?;
        if self.target.metadata.get("mutation_intent:source_topology")
            != Some(&encode(&self.target.topology)?)
        {
            return Err(bad("source topology changed"));
        }
        let count = interval(n, self.target.topology, &self.request)?;
        let added = inserted(&self.request)?;
        let shift = if self.request.edit_end > n {
            self.request.edit_start
        } else {
            0
        };
        if self.target.sequence.alphabet() != Alphabet::Dna
            || self.target.name != self.request.name
            || self.deleted.len() != count
            || self.target.sequence.len() != n - count + added.len()
            || self.origin_shift != shift
        {
            return Err(bad("intended reference dimensions changed"));
        }
        for (key, value) in [
            ("duplex_model", REFERENCE_MODEL.to_owned()),
            ("mutation:intent:reference_only", "true".into()),
            ("mutation_intent:engine", ENGINE_VERSION.into()),
            (
                "mutation_intent:source_revision",
                self.source_revision.clone(),
            ),
            (
                "mutation_intent:context_revision",
                self.context_revision.clone(),
            ),
            ("mutation_intent:request", encode(&self.request)?),
            ("mutation_intent:origin_shift", shift.to_string()),
            (
                "mutation_intent:sequence_blake3",
                blake3::hash(self.target.sequence.as_str().as_bytes())
                    .to_hex()
                    .to_string(),
            ),
            (
                "mutation_intent:deleted_blake3",
                blake3::hash(self.deleted.as_bytes()).to_hex().to_string(),
            ),
        ] {
            if self.target.metadata.get(key) != Some(&value) {
                return Err(bad("intended reference evidence binding changed"));
            }
        }
        assembly_recipe::hash(self)?;
        Ok(())
    }
}
