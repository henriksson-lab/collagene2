//! Documentary conditions attached to the exact material entering a plan step.
//! Retained hashes establish integrity; only fresh reassessment establishes the
//! current catalogue association. Neither establishes experimental compatibility.
use super::{ConditionsReport, EvidenceStatus, ReactionConditions, catalog};
use crate::chemistry::tracked_duplex::TrackedDuplex;
use collagene_model::{Error, Molecule};
use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct SourceReport {
    pub source_revision: String,
    pub material_revision: String,
    pub context_revision: String,
    pub report: ConditionsReport,
}
impl SourceReport {
    fn context(source: &str, material: &str, report: &ConditionsReport) -> String {
        let mut h = blake3::Hasher::new();
        h.update(b"plan-documentary-conditions-1\0");
        for value in [
            source,
            material,
            &report.request_revision,
            &report.batch_revision,
        ] {
            h.update(&(value.len() as u64).to_le_bytes());
            h.update(value.as_bytes());
        }
        h.finalize().to_hex().to_string()
    }
    pub fn assess(
        catalog: &catalog::Catalog,
        source: &Molecule,
        material: &TrackedDuplex,
        selection: &super::super::catalog::CatalogSelection,
        options: &ReactionConditions,
    ) -> Result<Self, Error> {
        collagene_model::bounds::model_preflight(&[source])?;
        if material.duplex.top != source.sequence || material.duplex.topology != source.topology {
            return Err(Error::Format(
                "documentary source differs from actual material top reference".into(),
            ));
        }
        let report = super::assess(catalog, selection, options)?;
        let source_revision = source.revision();
        let material_revision = material.summary()?.state_revision;
        let context_revision = Self::context(&source_revision, &material_revision, &report);
        Ok(Self {
            source_revision,
            material_revision,
            context_revision,
            report,
        })
    }
    /// Historical integrity only; does not confer current catalogue authority.
    pub fn validate_historical(&self) -> Result<(), Error> {
        self.report.validate()?;
        let valid_hash = |s: &str| s.len() == 64 && s.bytes().all(|b| b.is_ascii_hexdigit());
        if !valid_hash(&self.source_revision)
            || !valid_hash(&self.material_revision)
            || self.context_revision
                != Self::context(&self.source_revision, &self.material_revision, &self.report)
        {
            return Err(Error::Format(
                "documentary reaction-condition source binding differs".into(),
            ));
        }
        Ok(())
    }
    pub fn validate_for(
        &self,
        catalog: &catalog::Catalog,
        source: &Molecule,
        material: &TrackedDuplex,
        selection: &super::super::catalog::CatalogSelection,
        options: &ReactionConditions,
    ) -> Result<(), Error> {
        self.validate_historical()?;
        if self != &Self::assess(catalog, source, material, selection, options)? {
            return Err(Error::Format(
                "documentary reaction-condition evidence is stale".into(),
            ));
        }
        Ok(())
    }
    pub fn satisfies_request(&self) -> bool {
        self.report.documentary_status == EvidenceStatus::Known
    }
    /// Plain text: HTML callers must escape it; PDF callers use their text encoder.
    pub fn protocol_lines(&self) -> Result<Vec<String>, Error> {
        self.validate_historical()?;
        let mut lines = vec![
            format!("Documentary reaction-condition evidence: {:?}. Source {}; material {}; context {}.", self.report.documentary_status, self.source_revision, self.material_revision, self.context_revision),
            format!("Catalogue edition {}; content {}; request {}; report {}.", self.report.catalog_edition, self.report.catalog_content_hash, self.report.request_revision, self.report.batch_revision),
            format!("Requested buffer {}; temperature {}.{:03} C; incubation {} minutes; minimum documentary activity {} percent; activity context {:?}; temperature context {:?}.", self.report.options.buffer, self.report.options.temperature_milli_c / 1000, self.report.options.temperature_milli_c % 1000, self.report.options.incubation_minutes, self.report.options.minimum_activity_percent, self.report.options.activity_context, self.report.options.temperature_context),
            format!("Explicit reaction choices: {}", serde_json::to_string(&self.report.options).map_err(|e| Error::Format(e.to_string()))?),
            "Known satisfies only the exact requested documentary activity/temperature constraint. Star activity and incubation-duration compatibility remain unassessed; cleavage, PCR and experimental success are not established.".into(),
        ];
        for formulation in &self.report.formulations {
            lines.push(format!(
                "Formulation documentary record: {}",
                serde_json::to_string(formulation).map_err(|e| Error::Format(e.to_string()))?
            ));
        }
        for source in &self.report.sources {
            lines.push(format!(
                "Documentary source provenance: {}",
                serde_json::to_string(source).map_err(|e| Error::Format(e.to_string()))?
            ));
        }
        Ok(lines)
    }
}
