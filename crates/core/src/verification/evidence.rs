//! Verification report and evidence owners with their validation behavior.
use super::checks::{
    CheckFinding, CheckScope, CheckVerdict, MAX_FINDING_BYTES, MAX_FINDINGS, RuleId,
};
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
    #[serde(
        deserialize_with = "collagene_model::decode::bounded_vec::<_, CheckFinding, MAX_FINDINGS>"
    )]
    pub findings: Vec<CheckFinding>,
    pub assumptions: Vec<String>,
    pub planned_steps: usize,
    #[serde(deserialize_with = "collagene_model::decode::bounded_vec::<_, StepOutcome, 1000>")]
    pub steps: Vec<StepOutcome>,
    pub exhaustive_alternative_search: bool,
    /// Full per-step scientific evidence, including withheld TOPO contours.
    pub step_reports: Vec<StepReport>,
    /// Full processed graph evidence, including withheld material.
    pub physical_processing_reports: Vec<PhysicalProcessingStepReport>,
    /// Typed reference-design evidence, distinct from actual material state.
    pub reference_step_reports: Vec<ReferenceStepReport>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
pub enum StepState {
    Reconstructed,
    Rejected,
    Withheld,
    Blocked,
    NotEvaluated,
}
#[derive(Debug, Clone, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct StepOutcome {
    pub step_index: usize,
    pub state: StepState,
    #[serde(deserialize_with = "crate::encoding::required_option")]
    #[cfg_attr(feature = "json-schema", schemars(required))]
    pub product: Option<Molecule>,
    #[serde(deserialize_with = "collagene_model::decode::bounded_vec::<_, Source, 1000>")]
    pub blocked_by: Vec<Source>,
}
impl VerificationReport {
    pub fn product(&self, step_index: usize) -> Option<&Molecule> {
        self.steps
            .get(step_index)
            .filter(|step| step.step_index == step_index)?
            .product
            .as_ref()
    }
    /// The designated last step, never the last successful independent branch.
    pub fn final_product(&self) -> Option<&Molecule> {
        self.planned_steps
            .checked_sub(1)
            .and_then(|index| self.product(index))
    }
    pub fn reconstructed_count(&self) -> usize {
        self.steps
            .iter()
            .filter(|step| step.state == StepState::Reconstructed && step.product.is_some())
            .count()
    }
    /// Current documentary shape and integrity only. This neither replays the
    /// plan nor grants private chemical/source authority to imported reports.
    pub fn validate_documentary(&self) -> Result<(), Error> {
        use std::collections::BTreeSet;
        if self.planned_steps > 1000
            || self.steps.len() != self.planned_steps
            || self.findings.len() > MAX_FINDINGS
            || self.assumptions.len() > 64
            || self.assumptions.iter().any(|text| text.len() > 8192)
            || self.step_reports.len() > 1000
            || self.physical_processing_reports.len() > 64
            || self.reference_step_reports.len() > 64
        {
            return Err(Error::TooLarge);
        }
        // Native owned callers can bypass serde's location-depth admission.
        // Walk borrowed trees before recursive JSON counting or molecule checks,
        // preserving native-frame per-product and aggregate annotation limits.
        collagene_model::bounds::native_frame_preflight(
            self.steps.iter().filter_map(|step| step.product.as_ref()),
            64 * 1024 * 1024,
        )?;
        for finding in &self.findings {
            finding.validate()?;
        }
        evidence_bytes(self, 64 * 1024 * 1024)?;
        evidence_bytes(&self.findings, MAX_FINDING_BYTES)?;
        let hex = |text: &str| {
            text.len() == 64
                && text
                    .bytes()
                    .all(|b| b.is_ascii_digit() || (b'a'..=b'f').contains(&b))
        };
        let rejected_admission = self.planned_steps == 0
            && self.status != Status::ConsistentUnderAssumptions
            && self.findings.iter().any(|finding| {
                finding.rule == RuleId::Admission && finding.scope == CheckScope::Plan
            });
        if self.engine_version != super::ENGINE_VERSION
            || (!rejected_admission && (!hex(&self.plan_revision) || !hex(&self.target_revision)))
        {
            return Err(Error::Format(
                "verification report engine/binding shape changed".into(),
            ));
        }
        let mut bases = 0usize;
        for (index, step) in self.steps.iter().enumerate() {
            if step.step_index != index
                || (step.state == StepState::Reconstructed) != step.product.is_some()
                || (step.state == StepState::Blocked) != !step.blocked_by.is_empty()
            {
                return Err(Error::Format(
                    "verification outcome index/state/product contradiction".into(),
                ));
            }
            let mut blockers = BTreeSet::new();
            for source in &step.blocked_by {
                let key = match source {
                    Source::Input(id) if id.len() <= 128 && !id.chars().any(char::is_control) => {
                        format!("input:{id}")
                    }
                    Source::Step(producer)
                        if *producer < index
                            && self.steps[*producer].state != StepState::Reconstructed =>
                    {
                        format!("step:{producer}")
                    }
                    _ => return Err(Error::Format("invalid direct verification blocker".into())),
                };
                if !blockers.insert(key) {
                    return Err(Error::Format(
                        "duplicate direct verification blocker".into(),
                    ));
                }
            }
            if let Some(product) = &step.product {
                product.validate()?;
                bases = bases
                    .checked_add(product.sequence.len())
                    .ok_or(Error::TooLarge)?;
                if bases > 32_000_000 {
                    return Err(Error::TooLarge);
                }
            }
        }
        let mut state_checks = vec![[false; 4]; self.planned_steps];
        for finding in &self.findings {
            if let CheckScope::Step { step_index } = finding.scope {
                if let Some(checks) = state_checks.get_mut(step_index) {
                    checks[0] |= finding.rule == RuleId::StepCompletion
                        && finding.verdict == CheckVerdict::Pass;
                    checks[1] |= finding.verdict == CheckVerdict::Fail;
                    checks[2] |= finding.verdict == CheckVerdict::Indeterminate;
                    checks[3] |= finding.rule == RuleId::DependencyBlocked
                        && finding.verdict == CheckVerdict::Blocked;
                }
            }
            if matches!(finding.scope, CheckScope::Step{step_index}|CheckScope::Target{step_index} if step_index>=self.planned_steps)
            {
                return Err(Error::Format(
                    "verification check names an absent step".into(),
                ));
            }
        }
        for (step, checks) in self.steps.iter().zip(&state_checks) {
            let supported = match step.state {
                StepState::Reconstructed => checks[0],
                StepState::Rejected => checks[1],
                StepState::Withheld | StepState::NotEvaluated => checks[2],
                StepState::Blocked => checks[3],
            };
            if !supported {
                return Err(Error::Format(
                    "verification outcome lacks its named state predicate".into(),
                ));
            }
        }
        if self
            .findings
            .iter()
            .any(|f| f.verdict == CheckVerdict::Fail)
            && self.status != Status::Invalid
        {
            return Err(Error::Format(
                "failed verification predicate cannot claim non-invalid status".into(),
            ));
        }
        if self.status == Status::ConsistentUnderAssumptions
            && (self.planned_steps == 0
                || self.reconstructed_count() != self.planned_steps
                || self.findings.iter().any(|f| {
                    matches!(
                        f.verdict,
                        CheckVerdict::Fail | CheckVerdict::Blocked | CheckVerdict::Indeterminate
                    )
                })
                || !self
                    .findings
                    .iter()
                    .any(|f| f.rule == RuleId::TargetAgreement && f.verdict == CheckVerdict::Pass))
        {
            return Err(Error::Format(
                "positive verification status lacks complete outcomes/checks".into(),
            ));
        }
        let mut seen = BTreeSet::new();
        for step in &self.step_reports {
            if step.step_index >= self.planned_steps || !seen.insert(step.step_index) {
                return Err(Error::Format("invalid retained step evidence index".into()));
            }
            step.validate_documentary()?;
        }
        if !self.reference_step_reports.is_empty() {
            if self.reference_step_reports.len() > 64 {
                return Err(Error::TooLarge);
            }
            let mut seen = std::collections::BTreeSet::new();
            for step in &self.reference_step_reports {
                if step.step_index >= self.planned_steps || !seen.insert(step.step_index) {
                    return Err(Error::Format("Invalid reference step index".into()));
                }
                match &step.evidence {
                crate::verification::evidence::ReferenceStepEvidence::MultiSiteRecombination {
                    selection,
                    batch,
                } => {
                    batch.validate_binding(&selection.batch_revision)?;
                    if !batch.documents.iter().any(|p| p.selection == *selection) {
                        return Err(Error::Format(
                            "Reference selection is missing from its batch".into(),
                        ));
                    }
                }
            }
            }
        }
        if !self.physical_processing_reports.is_empty() {
            if self.physical_processing_reports.len() > 64 {
                return Err(Error::TooLarge);
            }
            let mut seen = std::collections::BTreeSet::new();
            for step in &self.physical_processing_reports {
                if step.step_index >= self.planned_steps || !seen.insert(step.step_index) {
                    return Err(Error::Format(
                        "Invalid physical processing step index".into(),
                    ));
                }
                match &step.evidence {
                crate::verification::evidence::PhysicalProcessingEvidence::LigationOutcomes(r) => {
                    r.validate_documentary()?;
                }
                crate::verification::evidence::PhysicalProcessingEvidence::DerivedOrigin(r) => {
                    r.validate_documentary()?;
                }
                crate::verification::evidence::PhysicalProcessingEvidence::RunProduct(r) => {
                    r.validate_documentary()?;
                }
                crate::verification::evidence::PhysicalProcessingEvidence::TopoActions(r) => {
                    r.validate_binding(&r.context_revision)?
                }
                crate::verification::evidence::PhysicalProcessingEvidence::OverlapExtension(r) => {
                    r.validate_documentary()?;
                }
                crate::verification::evidence::PhysicalProcessingEvidence::ThreePrimeResection(
                    report,
                ) => {
                    report.validate_documentary()?;
                }
                crate::verification::evidence::PhysicalProcessingEvidence::HomologyAnnealing(
                    report,
                ) => {
                    report.validate_documentary()?;
                }
                crate::verification::evidence::PhysicalProcessingEvidence::HomologyAssembly(
                    report,
                ) => {
                    report.validate_documentary()?;
                }
            }
            }
        }

        Ok(())
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct StepReport {
    pub step_index: usize,
    #[serde(deserialize_with = "crate::encoding::required_option")]
    #[cfg_attr(feature = "json-schema", schemars(required))]
    pub topo: Option<crate::cloning::topo::Product>,
    #[serde(deserialize_with = "crate::encoding::required_option")]
    #[cfg_attr(feature = "json-schema", schemars(required))]
    pub reaction_conditions: Option<crate::cloning::restriction::conditions::source::SourceReport>,
}

impl StepReport {
    /// Retained evidence integrity only. Fresh plan reconstruction remains required.
    pub fn validate_documentary(&self) -> Result<(), Error> {
        if self.step_index >= 1000 || self.topo.is_some() == self.reaction_conditions.is_some() {
            return Err(Error::Format("invalid or ambiguous step evidence".into()));
        }
        if let Some(product) = &self.topo {
            crate::cloning::topo::worker::validate_projection(product, &product.batch_revision)?;
        }
        if let Some(conditions) = &self.reaction_conditions {
            conditions.validate_historical()?;
        }
        Ok(())
    }
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
    #[serde(
        deserialize_with = "collagene_model::decode::bounded_vec::<_, crate::cloning::restriction::conditions::source::SourceReport, 32>"
    )]
    pub reaction_conditions: Vec<crate::cloning::restriction::conditions::source::SourceReport>,
    pub origin: crate::material::derivation::Capsule,
    pub material: MaterialStateEvidence,
}
impl DerivedOriginEvidence {
    pub fn validate_documentary(&self) -> Result<(), Error> {
        crate::material::derivation::reaction_conditions::validate(&self.reaction_conditions)?;
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
pub(crate) fn evidence_bytes(value: &impl Serialize, max: usize) -> Result<usize, Error> {
    let mut counter = EvidenceCounter { bytes: 0, max };
    serde_json::to_writer(&mut counter, value).map_err(|_| Error::TooLarge)?;
    Ok(counter.bytes)
}
