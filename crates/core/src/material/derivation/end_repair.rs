//! End-repair preview, revision binding and documentary validation.
use super::ENGINE_VERSION;
use super::SourceOrigin;
use super::bad;
use super::encoded;
use super::fresh_source;
use super::limitations;
use super::preview_receipt_shape;
use super::preview_source_shape;
use super::source_origin_preflight;
use crate::chemistry;
use crate::chemistry::tracked_duplex::ModificationSummary;
use crate::verification;
use collagene_model::Error;
use collagene_model::Molecule;
use serde::Deserialize;
use serde::Serialize;

dto! {pub struct RepairPreview{pub engine_version:String,pub context_revision:String,pub source_revision:String,pub molecule:Molecule,
#[cfg_attr(feature="json-schema",schemars(with="serde_json::Value"))]pub repair:chemistry::EndRepairReport,
pub source_duplex:chemistry::DuplexSummary,pub source_modifications:ModificationSummary,
pub incomplete:bool,pub findings:Vec<String>,pub findings_total:usize,pub findings_truncated:bool}}
pub fn repair_preview_revision(
    source: &Molecule,
    origin: &SourceOrigin,
    method: chemistry::EndRepairMethod,
) -> Result<String, Error> {
    source_origin_preflight(source, origin)?;
    crate::encoding::bounded_json::hash(&(
        ENGINE_VERSION,
        "end-repair-preview",
        verification::ENGINE_VERSION,
        source,
        origin,
        method,
    ))
}
pub fn preview_repair(
    source: &Molecule,
    origin: &SourceOrigin,
    method: chemistry::EndRepairMethod,
) -> Result<RepairPreview, Error> {
    let context_revision = repair_preview_revision(source, origin, method)?;
    let fresh = fresh_source(source, origin)?;
    let repair = chemistry::repair_ends(&fresh.tracked.duplex, method)?;
    let mut local = repair.warnings.clone();
    let mut molecule = verification::terminal::repaired_document(
        source,
        &repair,
        format!("{}_end_repair", source.name),
        &mut local,
    )?
    .molecule;
    molecule.metadata.insert(
        "source:material_derivation:preview_context".into(),
        context_revision.clone(),
    );
    local.push("Reference preview of explicit template-supported end repair. It restores no material authority, infers no phosphorylation or nick sealing, and checks no intended target.".into());
    let total = fresh
        .findings_total
        .checked_add(local.len())
        .ok_or(Error::TooLarge)?;
    let (findings, findings_total, findings_truncated) = limitations(
        local
            .into_iter()
            .chain(fresh.findings.iter().cloned())
            .collect(),
        total,
        fresh.findings_truncated,
    );
    let result = RepairPreview {
        engine_version: ENGINE_VERSION.into(),
        context_revision,
        source_revision: source.revision(),
        source_duplex: fresh.tracked.duplex.summary()?,
        source_modifications: fresh.tracked.summary()?,
        molecule,
        repair,
        incomplete: fresh.incomplete,
        findings,
        findings_total,
        findings_truncated,
    };
    result.validate_binding(&result.context_revision)?;
    Ok(result)
}
impl RepairPreview {
    pub fn validate_binding(&self, expected: &str) -> Result<(), Error> {
        preview_receipt_shape(
            &self.engine_version,
            &self.context_revision,
            expected,
            &self.source_revision,
            &self.findings,
            self.findings_total,
            self.findings_truncated,
        )?;
        preview_source_shape(&self.source_duplex, &self.source_modifications)?;
        self.repair.source.validate()?;
        self.repair.duplex.validate()?;
        self.molecule.validate()?;
        if self.repair.source.summary()?.top_length != self.source_duplex.top_length
            || crate::encoding::bounded_json::hash(&self.repair.source.summary()?)?
                != crate::encoding::bounded_json::hash(&self.source_duplex)?
            || self.molecule.sequence != self.repair.duplex.top
            || self.molecule.topology != self.repair.duplex.topology
            || self
                .molecule
                .metadata
                .get("source:material_derivation:preview_context")
                .map(String::as_str)
                != Some(expected)
        {
            return Err(bad(
                "repair preview native geometry/reference correspondence differs",
            ));
        }
        encoded(self, 32 * 1024 * 1024)?;
        Ok(())
    }
}
