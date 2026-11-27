//! Source-backed, explicitly requested ligase sealing after fresh TOPO transfers.
//! This additive graph step does not promote a contour to tracked duplex material.
use crate::{Error, assembly_recipe, topo};
use serde::{Deserialize, Serialize};
use std::collections::BTreeSet;
use topo::physical::{self, nick_sealing as nick};
pub const ENGINE_VERSION: &str = "topo-source-nick-sealing-1";
pub const MAX_REPORT_BYTES: usize = 32 * 1024 * 1024;
macro_rules! dto {($($t:tt)*)=>{#[derive(Clone,Debug,Serialize,Deserialize)]#[cfg_attr(feature="json-schema",derive(schemars::JsonSchema))]#[serde(deny_unknown_fields)]$($t)*};}
dto! {pub struct Request {pub preparation:topo::Request,pub sealing:nick::Request}}
dto! {pub struct Report {
    pub engine_version:String,pub context_revision:String,pub preparation_revision:String,
    pub request:Request,pub sources:Vec<topo::SourceIdentity>,
    pub preparation:topo::PreparationEvidence,pub sealing:nick::Report,
}}
fn bad(message: &str) -> Error {
    Error::Unsupported(format!("Source-bound TOPO nick sealing: {message}"))
}
fn physical_error(error: impl std::fmt::Display) -> Error {
    bad(&error.to_string())
}
fn digest(s: &str) -> bool {
    s.len() == 64
        && s.bytes()
            .all(|b| b.is_ascii_digit() || (b'a'..=b'f').contains(&b))
}
fn treatment_shape(r: &Request) -> Result<(), Error> {
    nick::request_preflight(&r.sealing, &r.preparation.engine_limits).map_err(physical_error)?;
    r.preparation
        .transfers
        .len()
        .checked_add(r.sealing.nicks.len())
        .filter(|n| *n <= 16 && *n <= r.preparation.engine_limits.max_steps)
        .ok_or(Error::TooLarge)?;
    Ok(())
}
/// Borrowed source/preparation/request bounds before molecule hashing, cloning,
/// PCR enumeration or graph construction. Actual terminal eligibility is fresh.
pub fn preflight(inputs: &[topo::Input], r: &Request) -> Result<(), Error> {
    treatment_shape(r)?;
    topo::source_preflight(inputs, &r.preparation)?;
    Ok(())
}
fn binding(
    preparation_revision: &str,
    sources: &[topo::SourceIdentity],
    r: &Request,
) -> Result<String, Error> {
    assembly_recipe::hash(&(
        ENGINE_VERSION,
        nick::ENGINE,
        preparation_revision,
        sources,
        r,
    ))
}
pub fn revision(inputs: &[topo::Input], r: &Request) -> Result<String, Error> {
    preflight(inputs, r)?;
    binding(
        &topo::batch_revision(inputs, &r.preparation)?,
        &topo::documentary_sources(inputs, &r.preparation)?,
        r,
    )
}
/// One preparation (including any exact PCR) and one fresh transfer reconstruction.
/// No caller-provided graph or previous positive report is accepted as input.
pub fn inspect(inputs: &[topo::Input], r: &Request) -> Result<Report, Error> {
    preflight(inputs, r)?;
    let preparation_revision = topo::batch_revision(inputs, &r.preparation)?;
    let prepared = topo::prepare(inputs, &r.preparation)?;
    let sources = topo::identities(&prepared.inputs);
    let context_revision = binding(&preparation_revision, &sources, r)?;
    let fresh = topo::execute(&prepared, &r.preparation.engine_limits)?;
    let sealing = nick::seal_fresh(
        &prepared.evidence.substrates,
        &prepared.effective_transfers,
        &r.preparation.engine_limits,
        &r.sealing,
        fresh,
    )
    .map_err(physical_error)?;
    let preparation = prepared.evidence;
    let result = Report {
        engine_version: ENGINE_VERSION.into(),
        context_revision: context_revision.clone(),
        preparation_revision,
        request: r.clone(),
        sources,
        preparation,
        sealing,
    };
    result.validate_binding(&context_revision)?;
    Ok(result)
}
impl Report {
    /// Read-only structural/source binding, not a chemical replay. Historical or
    /// imported evidence never authorizes physical-material promotion or outputs.
    pub fn validate_binding(&self, expected: &str) -> Result<(), Error> {
        treatment_shape(&self.request)?;
        if self.request.preparation.recipe.schema_version != 1
            || self.request.preparation.recipe.materials.is_empty()
            || self.request.preparation.recipe.materials.len() > 16
            || self.request.preparation.transfers.len() > 16
        {
            return Err(Error::TooLarge);
        }
        bounded_with_limit(&self.request, MAX_REPORT_BYTES - 65536)?;
        if self.engine_version != ENGINE_VERSION
            || !digest(expected)
            || self.context_revision != expected
            || !digest(&self.preparation_revision)
            || self.sources.is_empty()
            || self.sources.len() > 64
            || self.preparation.maps.len() > 4096
            || self.preparation.not_represented.len() > 4096
            || self.preparation.assumptions.len() > 128
        {
            return Err(bad("report source/request shape changed"));
        }
        let mut roles = BTreeSet::new();
        for s in &self.sources {
            if s.role.trim().is_empty()
                || s.role.len() > 256
                || !roles.insert(s.role.as_str())
                || !digest(&s.revision)
            {
                return Err(bad("report source identities changed"));
            }
        }
        let original_context = assembly_recipe::hash(&(
            topo::ENGINE_VERSION,
            physical::ENGINE,
            physical::PRIMARY_MODEL,
            physical::projection::ENGINE,
            &self.sources,
            &self.request.preparation,
        ))?;
        if original_context != self.preparation_revision
            || binding(&self.preparation_revision, &self.sources, &self.request)? != expected
        {
            return Err(bad("source/preparation/request context changed"));
        }
        if self.preparation.substrates.strands.len() > 64
            || self.preparation.substrates.pairing.len() > 1024
        {
            return Err(Error::TooLarge);
        }
        let mut prepared_bases = 0usize;
        for strand in &self.preparation.substrates.strands {
            prepared_bases = prepared_bases
                .checked_add(strand.sequence.len())
                .filter(|n| *n <= 100_000)
                .ok_or(Error::TooLarge)?;
        }
        physical::preflight(
            &self.preparation.substrates,
            &topo::effective_transfers(&self.preparation, &self.request.preparation)?,
            &self.request.preparation.engine_limits,
        )
        .map_err(physical_error)?;
        for map in &self.preparation.maps {
            let strand = self
                .preparation
                .substrates
                .strands
                .get(map.strand)
                .ok_or_else(|| bad("preparation mapping strand is absent"))?;
            if map.length == 0
                || map
                    .strand_start
                    .checked_add(map.length)
                    .is_none_or(|e| e > strand.sequence.len())
            {
                return Err(bad("preparation mapping coordinates changed"));
            }
            match &map.origin {
                topo::Origin::ReferenceCorrespondence(interval) => {
                    if !roles.contains(interval.source_role.as_str())
                        || interval.length != map.length
                    {
                        return Err(bad("preparation mapping source changed"));
                    }
                }
                topo::Origin::ExactPcrCopy { template_role, .. } => {
                    if !roles.contains(template_role.as_str()) {
                        return Err(bad("PCR mapping source changed"));
                    }
                }
                topo::Origin::SynthesizedTail { .. } => {}
            }
        }
        for gap in &self.preparation.not_represented {
            if !roles.contains(gap.source_role.as_str())
                || gap.start >= gap.end
                || gap.end > 100_000
            {
                return Err(bad("unrepresented reference interval changed"));
            }
        }
        if self.preparation.assumptions.iter().any(|s| s.len() > 4096) {
            return Err(Error::TooLarge);
        }
        nick::validate_report(
            &self.preparation.substrates,
            &topo::effective_transfers(&self.preparation, &self.request.preparation)?,
            &self.request.preparation.engine_limits,
            &self.request.sealing,
            &self.sealing,
        )
        .map_err(physical_error)?;
        bounded(self)?;
        Ok(())
    }
}
struct Counter {
    bytes: usize,
    limit: usize,
}
impl std::io::Write for Counter {
    fn write(&mut self, b: &[u8]) -> std::io::Result<usize> {
        self.bytes = self
            .bytes
            .checked_add(b.len())
            .filter(|n| *n <= self.limit)
            .ok_or_else(|| std::io::Error::other("bounded TOPO sealing report"))?;
        Ok(b.len())
    }
    fn flush(&mut self) -> std::io::Result<()> {
        Ok(())
    }
}
fn bounded(value: &impl Serialize) -> Result<(), Error> {
    bounded_with_limit(value, MAX_REPORT_BYTES)
}
fn bounded_with_limit(value: &impl Serialize, limit: usize) -> Result<(), Error> {
    serde_json::to_writer(Counter { bytes: 0, limit }, value).map_err(|_| Error::TooLarge)
}
/// Fresh source preparation and graph chemistry replay; never an import receipt.
pub fn replay(inputs: &[topo::Input], r: &Request, report: &Report) -> Result<(), Error> {
    let expected = revision(inputs, r)?;
    report.validate_binding(&expected)?;
    let current = inspect(inputs, r)?;
    if assembly_recipe::hash(&current)? != assembly_recipe::hash(report)? {
        return Err(bad("evidence differs from fresh source reconstruction"));
    }
    Ok(())
}
