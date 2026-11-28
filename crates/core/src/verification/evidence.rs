//! Verification report and evidence owners with their validation behavior.
use crate::{
    chemistry::{self, tracked_duplex as modified_duplex},
    planning::plan::{OrientedSource, Source},
};
use collagene_model::{Error, Molecule};
use serde::{Deserialize, Serialize};
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
pub enum Status {
    Invalid,
    ConsistentUnderAssumptions,
    Incomplete,
}
#[derive(Debug, Clone, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct VerificationReport {
    pub plan_revision: String,
    pub target_revision: String,
    pub engine_version: String,
    pub status: Status,
    pub findings: Vec<String>,
    pub assumptions: Vec<String>,
    pub products: Vec<Molecule>,
    pub exhaustive_alternative_search: bool,
    /// Full per-step scientific evidence, including withheld TOPO contours.
    pub step_reports: Vec<StepReport>,
    /// Full processed graph evidence, including withheld material.
    pub physical_processing_reports: Vec<PhysicalProcessingStepReport>,
    /// Typed reference-design evidence, distinct from actual material state.
    pub reference_step_reports: Vec<ReferenceStepReport>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct StepReport {
    pub step_index: usize,
    #[serde(deserialize_with = "crate::encoding::required_option")]
    #[cfg_attr(feature = "json-schema", schemars(required))]
    pub topo: Option<crate::cloning::topo::Product>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
pub enum PhysicalProcessingEvidence {
    TopoActions(crate::cloning::topo::actions::Report),
    RunProduct(crate::material::capsule::ReplayEvidence),
    OverlapExtension(crate::cloning::overlap::extension::Geometry),
    HomologyAssembly(crate::cloning::overlap::assembly::Report),
    ThreePrimeResection(crate::cloning::overlap::resection::Report),
    HomologyAnnealing(crate::cloning::overlap::annealing::Report),
    DerivedOrigin(DerivedOriginEvidence),
    LigationOutcomes(LigationOutcomeEvidence),
}
#[derive(Debug, Clone, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct LigationOutcomeBinding {
    pub key: String,
    pub source: Source,
    pub source_revision: String,
    pub material_revision: String,
}
#[derive(Debug, Clone, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct LigationOutcomeEvidence {
    pub bindings: Vec<LigationOutcomeBinding>,
    pub report: crate::cloning::ligation::outcomes::Report,
}
impl LigationOutcomeEvidence {
    pub fn validate_documentary(&self) -> Result<(), Error> {
        self.report
            .validate_binding(&self.report.context_revision)?;
        if self.bindings.len() != self.report.sources.len() {
            return Err(Error::Format(
                "Ligation outcome source count changed".into(),
            ));
        }
        for (i, binding) in self.bindings.iter().enumerate() {
            if binding.key != format!("fragment-{}", i + 1)
                || !self
                    .report
                    .sources
                    .iter()
                    .any(|s| s.key == binding.key && s.state_revision == binding.material_revision)
                || [&binding.source_revision, &binding.material_revision]
                    .iter()
                    .any(|s| s.len() != 64 || !s.bytes().all(|b| b.is_ascii_hexdigit()))
                || matches!(&binding.source,Source::Input(id) if id.is_empty() || id.len()>128 || id.chars().any(char::is_control))
            {
                return Err(Error::Format(
                    "Ligation outcome exact source/material binding changed".into(),
                ));
            }
        }
        evidence_bytes(self, 32 * 1024 * 1024)?;
        Ok(())
    }
    pub fn validate_for_operation(
        &self,
        fragments: &[OrientedSource],
        circular: bool,
        allow_nicked: bool,
        options: &crate::cloning::ligation::outcomes::Options,
    ) -> Result<(), Error> {
        self.validate_documentary()?;
        if fragments.len() != self.bindings.len()
            || self.report.request.circular != circular
            || self.report.request.allow_nicked != allow_nicked
            || crate::encoding::bounded_json::hash(options)?
                != crate::encoding::bounded_json::hash(&self.report.request.options)?
        {
            return Err(Error::Format(
                "Ligation outcomes declaration changed".into(),
            ));
        }
        for ((fragment, binding), placement) in fragments
            .iter()
            .zip(&self.bindings)
            .zip(&self.report.request.intended)
        {
            if crate::encoding::bounded_json::hash(&fragment.source)?
                != crate::encoding::bounded_json::hash(&binding.source)?
                || placement.key != binding.key
                || placement.reverse != fragment.reverse
            {
                return Err(Error::Format(
                    "Ligation intended path/source correspondence changed".into(),
                ));
            }
        }
        Ok(())
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct MaterialStateEvidence {
    pub material_revision: String,
    pub duplex: chemistry::DuplexSummary,
    pub modifications: modified_duplex::ModificationSummary,
    pub origin_incomplete: bool,
    pub origin_findings: Vec<String>,
    pub origin_findings_total: usize,
    pub origin_findings_truncated: bool,
}
impl MaterialStateEvidence {
    /// Documentary shape and exact recipe binding, never imported material authority.
    pub fn validate_against(
        &self,
        expected_state: &str,
        expected_incomplete: bool,
    ) -> Result<(), Error> {
        let hex = |value: &str| {
            value.len() == 64
                && value
                    .bytes()
                    .all(|b| b.is_ascii_digit() || (b'a'..=b'f').contains(&b))
        };
        let count = self
            .modifications
            .known_modified
            .checked_add(self.modifications.unknown_eligible)
            .and_then(|n| n.checked_add(self.modifications.known_unmodified))
            .ok_or(Error::TooLarge)?;
        let total = self
            .duplex
            .top_length
            .checked_add(self.duplex.bottom_length)
            .ok_or(Error::TooLarge)?;
        if !hex(&self.material_revision)
            || self.material_revision != self.modifications.state_revision
            || self.material_revision != expected_state
            || self.origin_incomplete != expected_incomplete
            || self.modifications.engine_version != modified_duplex::ENGINE_VERSION
            || self.duplex.top_length == 0
            || self.duplex.bottom_length == 0
            || self.duplex.top_length > collagene_model::MAX_SEQUENCE_LEN
            || self.duplex.bottom_length > collagene_model::MAX_SEQUENCE_LEN
            || self.modifications.top_length != self.duplex.top_length
            || self.modifications.bottom_length != self.duplex.bottom_length
            || count > total
            || self.origin_findings.len() > 128
            || self
                .origin_findings
                .iter()
                .any(|s| s.len() > 8192 || s.chars().count() > 2048)
            || self.origin_findings_total < self.origin_findings.len()
            || (!self.origin_findings_truncated
                && self.origin_findings_total > self.origin_findings.len())
        {
            return Err(Error::Format(
                "Derived material evidence changed its bounded source/material binding".into(),
            ));
        }
        evidence_bytes(self, 32 * 1024 * 1024)?;
        Ok(())
    }
}
#[derive(Debug, Clone, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct DerivedOriginEvidence {
    pub origin: crate::material::derivation::Capsule,
    pub material: MaterialStateEvidence,
}
impl DerivedOriginEvidence {
    pub fn validate_documentary(&self) -> Result<(), Error> {
        self.origin.validate_documentary()?;
        let record = crate::material::derivation::load_record_shallow(&self.origin.record_bytes)?;
        self.material
            .validate_against(&record.product_state_revision, record.incomplete)?;
        evidence_bytes(self, 32 * 1024 * 1024)?;
        Ok(())
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct PhysicalProcessingStepReport {
    pub step_index: usize,
    pub evidence: PhysicalProcessingEvidence,
}
#[derive(Debug, Clone, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
pub enum ReferenceStepEvidence {
    MultiSiteRecombination {
        selection: crate::cloning::recombination::Selection,
        batch: crate::cloning::recombination::multisite::Batch,
    },
}
#[derive(Debug, Clone, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct ReferenceStepReport {
    pub step_index: usize,
    pub evidence: ReferenceStepEvidence,
}

struct EvidenceCounter {
    bytes: usize,
    max: usize,
}
impl std::io::Write for EvidenceCounter {
    fn write(&mut self, bytes: &[u8]) -> std::io::Result<usize> {
        self.bytes = self
            .bytes
            .checked_add(bytes.len())
            .filter(|&n| n <= self.max)
            .ok_or_else(|| std::io::Error::other("bounded plan step evidence"))?;
        Ok(bytes.len())
    }
    fn flush(&mut self) -> std::io::Result<()> {
        Ok(())
    }
}
pub(super) fn evidence_bytes(value: &impl Serialize, max: usize) -> Result<usize, Error> {
    let mut counter = EvidenceCounter { bytes: 0, max };
    serde_json::to_writer(&mut counter, value).map_err(|_| Error::TooLarge)?;
    Ok(counter.bytes)
}
