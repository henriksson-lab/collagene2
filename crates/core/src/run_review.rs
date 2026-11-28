//! Prepare an exact retained run for editing/reverification, never physical replay.
use crate::{
    run_lineage::{self, EventBinding, InputEventBinding, RecordedAssessment},
    transfer::Candidate,
    verification::{self, Operation},
    workspace::Workspace,
};
use collagene_model::{Error, Molecule};
use serde::{Deserialize, Serialize};
use std::io::Write;

pub const ENGINE_VERSION: &str = "run-review-1";
pub const MAX_REVIEW_BYTES: usize = 32 * 1024 * 1024;

/// A documentary receipt. The candidate has no adopted verification status.
#[derive(Debug, Clone, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct PreparedReview {
    pub engine_version: String,
    pub context_revision: String,
    pub source_workspace_fingerprint: String,
    pub run_id: String,
    pub run_attachment_hash: String,
    #[cfg_attr(feature = "json-schema", schemars(with = "serde_json::Value"))]
    pub candidate: Candidate,
    pub target: Molecule,
    pub target_event: EventBinding,
    pub inputs: Vec<InputEventBinding>,
    pub recorded_assessment: RecordedAssessment,
    pub current_verification_engine: String,
    pub recorded_engine_matches: bool,
    pub has_saved_material_origins: bool,
    pub warnings: Vec<String>,
}
struct HashWriter {
    hash: blake3::Hasher,
    remaining: usize,
}
impl Write for HashWriter {
    fn write(&mut self, b: &[u8]) -> std::io::Result<usize> {
        if b.len() > self.remaining {
            return Err(std::io::Error::other("review exceeds byte bound"));
        }
        self.remaining -= b.len();
        self.hash.update(b);
        Ok(b.len())
    }
    fn flush(&mut self) -> std::io::Result<()> {
        Ok(())
    }
}
fn hash<T: Serialize>(value: &T) -> Result<String, Error> {
    let mut writer = HashWriter {
        hash: blake3::Hasher::new(),
        remaining: MAX_REVIEW_BYTES,
    };
    writer.hash.update(b"collagene-run-review-v1\0");
    serde_json::to_writer(&mut writer, value).map_err(|_| Error::TooLarge)?;
    Ok(writer.hash.finalize().to_hex().to_string())
}
/// Request identity can be computed without preparing/cloning a candidate.
pub fn request_revision(workspace: &Workspace, run_id: &str) -> Result<String, Error> {
    request_key(&run_lineage::workspace_fingerprint(workspace)?, run_id)
}
fn request_key(workspace_fingerprint: &str, run_id: &str) -> Result<String, Error> {
    let digest = |s: &str| {
        s.len() == 64
            && s.bytes()
                .all(|c| c.is_ascii_hexdigit() && !c.is_ascii_uppercase())
    };
    if !digest(workspace_fingerprint) || !digest(run_id) {
        return Err(bad("invalid workspace/run identity"));
    }
    hash(&(
        ENGINE_VERSION,
        "request",
        verification::ENGINE_VERSION,
        workspace_fingerprint,
        run_id,
    ))
}
fn bad(message: &str) -> Error {
    Error::Format(format!("saved run review: {message}"))
}
fn resolve<'a>(workspace: &'a Workspace, binding: &EventBinding) -> Result<&'a Molecule, Error> {
    let event = workspace
        .history
        .iter()
        .find(|event| event.id == binding.event_id)
        .ok_or_else(|| bad("exact retained event is missing"))?;
    if event.document_id != binding.document_id || event.revision != binding.revision {
        return Err(bad("retained document/event/revision binding differs"));
    }
    workspace
        .revisions
        .get(&binding.revision)
        .ok_or_else(|| bad("retained revision is missing"))
}

/// Resolve the saved target and inputs by exact native events, not current heads or DNA aliases.
/// No verification, material reconstruction, workspace mutation or status adoption occurs.
pub fn prepare(workspace: &Workspace, run_id: &str) -> Result<PreparedReview, Error> {
    let source_workspace_fingerprint = run_lineage::workspace_fingerprint(workspace)?;
    let record = run_lineage::load_run(workspace, run_id)?;
    let mut attachment_hash = None;
    for (key, attachment) in &workspace.attachments {
        if attachment.kind == run_lineage::WORKSPACE_ATTACHMENT_KIND
            && run_lineage::load_run_record(&attachment.bytes)?.run_id == run_id
        {
            if attachment_hash.replace(key.clone()).is_some() {
                return Err(bad("duplicate retained run"));
            }
        }
    }
    let run_attachment_hash = attachment_hash.ok_or_else(|| bad("run attachment missing"))?;
    let target = resolve(workspace, &record.target)?;
    if target.revision() != record.plan.target_revision {
        return Err(bad("independent target binding differs"));
    }
    for (input, binding) in record.plan.inputs.iter().zip(&record.inputs) {
        let source = resolve(workspace, &binding.event)?;
        if input.id != binding.input_id
            || source.revision() != input.revision
            || hash(source)? != hash(&input.molecule)?
        {
            return Err(bad("exact original input snapshot differs"));
        }
    }
    // Structural package validation only; saved material is deliberately not replayed.
    let candidate = Candidate::new(record.plan.clone());
    candidate.validate()?;
    if candidate.plan_revision != record.plan_revision {
        return Err(bad("saved plan changed during preparation"));
    }
    let has_saved_material_origins = candidate.plan.steps.iter().any(|step| {
        matches!(
            step,
            Operation::ReconstructRunProduct { .. } | Operation::ReconstructMaterial { .. }
        )
    });
    let recorded_engine_matches = record.assessment.engine_version == verification::ENGINE_VERSION;
    let mut warnings = vec!["Recorded status is historical evidence. Review and independently reverify the candidate against its retained intended target before accepting current conclusions.".into()];
    if !recorded_engine_matches {
        warnings.push("The saved scientific engine differs from the current engine. No historical result or physical product was adopted.".into());
    }
    if has_saved_material_origins {
        warnings.push("The plan retains exact saved-material recipes. Current verification may refuse nested historical material models; this preparation neither migrates nor bypasses them.".into());
    }
    let mut result = PreparedReview {
        engine_version: ENGINE_VERSION.into(),
        context_revision: String::new(),
        source_workspace_fingerprint,
        run_id: record.run_id,
        run_attachment_hash,
        candidate,
        target: target.clone(),
        target_event: record.target,
        inputs: record.inputs,
        recorded_assessment: record.assessment,
        current_verification_engine: verification::ENGINE_VERSION.into(),
        recorded_engine_matches,
        has_saved_material_origins,
        warnings,
    };
    result.context_revision = result.receipt_revision()?;
    Ok(result)
}
impl PreparedReview {
    fn receipt_revision(&self) -> Result<String, Error> {
        // Exhaustive destructuring makes a newly added receipt field require an
        // explicit decision here; no candidate/sequence cloning is needed.
        let Self {
            engine_version,
            context_revision: _,
            source_workspace_fingerprint,
            run_id,
            run_attachment_hash,
            candidate,
            target,
            target_event,
            inputs,
            recorded_assessment,
            current_verification_engine,
            recorded_engine_matches,
            has_saved_material_origins,
            warnings,
        } = self;
        hash(&(
            engine_version,
            source_workspace_fingerprint,
            run_id,
            run_attachment_hash,
            candidate,
            target,
            target_event,
            inputs,
            recorded_assessment,
            current_verification_engine,
            recorded_engine_matches,
            has_saved_material_origins,
            warnings,
        ))
    }
    /// Validate a worker response without adopting its historical assessment.
    pub fn validate_request_binding(&self, revision: &str) -> Result<(), Error> {
        if self.engine_version != ENGINE_VERSION
            || self.current_verification_engine != verification::ENGINE_VERSION
            || request_key(&self.source_workspace_fingerprint, &self.run_id)? != revision
            || self.context_revision != self.receipt_revision()?
            || self.target_event.revision != self.target.revision()
            || self.target_event.revision != self.candidate.plan.target_revision
            || self.inputs.len() != self.candidate.plan.inputs.len()
            || self
                .inputs
                .iter()
                .zip(&self.candidate.plan.inputs)
                .any(|(b, i)| b.input_id != i.id || b.event.revision != i.revision)
        {
            return Err(bad("review request or receipt changed"));
        }
        self.candidate.validate()?;
        self.target.validate()?;
        Ok(())
    }
    /// Rebinds the complete receipt to the exact current workspace context. No scientific replay.
    pub fn validate_for_workspace(&self, workspace: &Workspace) -> Result<(), Error> {
        let expected = prepare(workspace, &self.run_id)?;
        if hash(self)? != hash(&expected)? {
            return Err(bad("review is stale or has been modified"));
        }
        Ok(())
    }
}
