//! Portable, versioned revision snapshots at the browser/MCP boundary.
//! Packages carry data only. Verification status is always recomputed by the recipient.
use crate::{planning::plan::DesignPlan, protocol};
use collagene_model::{Error, Molecule};
use serde::{Deserialize, Serialize};
use std::collections::{BTreeMap, BTreeSet};

const MAX_PLANNING_DOCUMENTS: usize = 1000;

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PlanningWorkspace {
    pub target_revision: String,
    #[serde(
        deserialize_with = "collagene_model::decode::bounded_vec::<_, Molecule, MAX_PLANNING_DOCUMENTS>"
    )]
    pub documents: Vec<Molecule>,
}
impl PlanningWorkspace {
    pub fn validate(&self) -> Result<(), Error> {
        if self.documents.is_empty() || self.documents.len() > MAX_PLANNING_DOCUMENTS {
            return Err(Error::Format(
                "planning snapshot needs 1–1000 documents".into(),
            ));
        }
        collagene_model::bounds::native_frame_preflight(
            &self.documents,
            crate::formats::sequence::MAX_IMPORT_BYTES,
        )?;
        let mut revisions = BTreeSet::new();
        for m in &self.documents {
            m.validate()?;
            if !revisions.insert(m.revision()) {
                return Err(Error::Format(
                    "duplicate document revision in planning snapshot".into(),
                ));
            }
        }
        if !revisions.contains(&self.target_revision) {
            return Err(Error::Format(
                "planning target is not one of the snapshot's documents".into(),
            ));
        }
        Ok(())
    }
}
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Candidate {
    pub plan_revision: String,
    pub plan: DesignPlan,
}
impl Candidate {
    pub fn new(plan: DesignPlan) -> Self {
        Self {
            plan_revision: if collagene_model::bounds::native_frame_preflight(
                plan.inputs.iter().map(|input| &input.molecule),
                crate::formats::sequence::MAX_IMPORT_BYTES,
            )
            .is_ok()
                && crate::material::capsule::plan_preflight(&plan).is_ok()
                && crate::verification::sources::validate_reference_steps(&plan).is_ok()
            {
                crate::planning::plan::DesignPlan::revision(&plan)
            } else {
                String::new()
            },
            plan,
        }
    }
    pub fn validate(&self) -> Result<(), Error> {
        let plan = &self.plan;
        collagene_model::bounds::native_frame_preflight(
            plan.inputs.iter().map(|input| &input.molecule),
            crate::formats::sequence::MAX_IMPORT_BYTES,
        )?;
        crate::material::capsule::plan_preflight(plan)?;
        crate::verification::sources::validate_reference_steps(plan)?;
        if plan.name.trim().is_empty()
            || plan.name.len() > 1024
            || plan.steps.is_empty()
            || plan.steps.len() > 1000
            || plan.inputs.len() > 1000
        {
            return Err(Error::Format(
                "candidate needs a name, 1–1000 steps and at most 1000 inputs".into(),
            ));
        }
        if plan.target_revision.len() != 64
            || !plan.target_revision.bytes().all(|b| b.is_ascii_hexdigit())
        {
            return Err(Error::Format("invalid candidate target revision".into()));
        }
        if self.plan_revision != crate::planning::plan::DesignPlan::revision(plan) {
            return Err(Error::Format("candidate plan revision mismatch".into()));
        }
        let mut ids = BTreeSet::new();
        for input in &plan.inputs {
            if input.id.trim().is_empty() || input.id.len() > 128 || !ids.insert(&input.id) {
                return Err(Error::Format(
                    "empty, oversized or duplicate candidate input ID".into(),
                ));
            }
            input.molecule.validate()?;
            if input.revision != input.molecule.revision() {
                return Err(Error::Format(format!(
                    "candidate input revision mismatch: {}",
                    input.id
                )));
            }
        }
        Ok(())
    }
    /// Use only explicitly open current documents, never silently replace them with package data.
    pub fn bind(&self, target: &Molecule, documents: &[Molecule]) -> Result<DesignPlan, Error> {
        self.validate()?;
        target.validate()?;
        if target.revision() != self.plan.target_revision {
            return Err(Error::Format("The active target differs from this candidate's bound revision. Select the original target or export a new planning snapshot.".into()));
        }
        let current: BTreeMap<_, _> = documents.iter().map(|m| (m.revision(), m)).collect();
        let mut plan = self.plan.clone();
        for input in &mut plan.inputs {
            let source=current.get(&input.revision).ok_or_else(||Error::Format(format!("Input '{}' is missing or has changed. Open its original revision or export a new planning snapshot.",input.id)))?;
            source.validate()?;
            input.molecule = (*source).clone();
        }
        Ok(plan)
    }
}

pub fn save_workspace(workspace: &PlanningWorkspace) -> Result<Vec<u8>, Error> {
    workspace.validate()?;
    protocol::envelope::encode(b"CGD\0", workspace)
}
pub fn load_workspace(bytes: &[u8]) -> Result<PlanningWorkspace, Error> {
    let workspace: PlanningWorkspace = protocol::envelope::decode(b"CGD\0", bytes)?;
    workspace.validate()?;
    Ok(workspace)
}
pub fn save_candidate(candidate: &Candidate) -> Result<Vec<u8>, Error> {
    candidate.validate()?;
    protocol::envelope::encode(b"CGP\0", candidate)
}
pub fn load_candidate(bytes: &[u8]) -> Result<Candidate, Error> {
    let candidate = decode_candidate_shallow(bytes)?;
    candidate.validate()?;
    Ok(candidate)
}
/// Bounded current-envelope decode only. The mixed-origin ledger charges and
/// validates each plan without recursively reentering Candidate::validate.
pub(crate) fn decode_candidate_shallow(bytes: &[u8]) -> Result<Candidate, Error> {
    protocol::envelope::decode(b"CGP\0", bytes)
}
