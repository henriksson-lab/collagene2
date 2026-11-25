//! Portable, versioned revision snapshots at the browser/MCP boundary.
//! Packages carry data only. Verification status is always recomputed by the recipient.
use crate::{
    Error, Molecule, protocol,
    verification::{DesignPlan, plan_revision},
};
use serde::{Deserialize, Serialize};
use std::collections::{BTreeMap, BTreeSet};

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PlanningWorkspace {
    pub target_revision: String,
    pub documents: Vec<Molecule>,
}
impl PlanningWorkspace {
    pub fn validate(&self) -> Result<(), Error> {
        if self.documents.is_empty() || self.documents.len() > 1000 {
            return Err(Error::Format(
                "planning snapshot needs 1–1000 documents".into(),
            ));
        }
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
            plan_revision: if crate::verification::validate_reference_steps(&plan).is_ok() {
                plan_revision(&plan)
            } else {
                String::new()
            },
            plan,
        }
    }
    pub fn validate(&self) -> Result<(), Error> {
        let plan = &self.plan;
        crate::verification::validate_reference_steps(plan)?;
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
        if self.plan_revision != plan_revision(plan) {
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
    protocol::encode(b"CGD\0", workspace)
}
pub fn load_workspace(bytes: &[u8]) -> Result<PlanningWorkspace, Error> {
    let workspace: PlanningWorkspace = protocol::decode(b"CGD\0", bytes)?;
    workspace.validate()?;
    Ok(workspace)
}
pub fn save_candidate(candidate: &Candidate) -> Result<Vec<u8>, Error> {
    candidate.validate()?;
    protocol::encode(b"CGP\0", candidate)
}
pub fn load_candidate(bytes: &[u8]) -> Result<Candidate, Error> {
    let candidate: Candidate = protocol::decode(b"CGP\0", bytes)?;
    if u16::from_le_bytes([bytes[4], bytes[5]]) < 14
        && candidate
            .plan
            .steps
            .iter()
            .any(|s| matches!(s, crate::verification::Operation::TopoActions { .. }))
    {
        return Err(Error::Unsupported(
            "Accounted chemical actions require candidate schema14".into(),
        ));
    }
    if u16::from_le_bytes([bytes[4], bytes[5]]) < 13
        && candidate
            .plan
            .steps
            .iter()
            .any(|step| matches!(step, crate::verification::Operation::ProcessedTopo { .. }))
    {
        return Err(Error::Unsupported(
            "Processed TOPO material requires candidate schema 13".into(),
        ));
    }
    if u16::from_le_bytes([bytes[4], bytes[5]]) < 12
        && candidate.plan.steps.iter().any(|step| {
            matches!(
                step,
                crate::verification::Operation::SelectedGoldenGate { .. }
            )
        })
    {
        return Err(Error::Unsupported(
            "Key-bound Golden Gate selection requires candidate schema 12".into(),
        ));
    }
    if u16::from_le_bytes([bytes[4], bytes[5]]) < 11
        && candidate
            .plan
            .steps
            .iter()
            .any(|step| matches!(step, crate::verification::Operation::CatalogDigest { .. }))
    {
        return Err(Error::Unsupported(
            "Catalog-pinned digest requires candidate schema 11".into(),
        ));
    }
    if u16::from_le_bytes([bytes[4], bytes[5]]) < 10
        && candidate.plan.steps.iter().any(|step| {
            matches!(
                step,
                crate::verification::Operation::DeclareSourceTermini { .. }
            )
        })
    {
        return Err(Error::Unsupported(
            "Declared source termini require candidate schema 10".into(),
        ));
    }
    if u16::from_le_bytes([bytes[4], bytes[5]]) < 9
        && candidate
            .plan
            .steps
            .iter()
            .any(|step| matches!(step, crate::verification::Operation::CheckedTopo { .. }))
    {
        return Err(Error::Unsupported(
            "Explicit TOPO target acceptance requires candidate schema 9".into(),
        ));
    }
    if u16::from_le_bytes([bytes[4], bytes[5]]) < 8
        && candidate
            .plan
            .steps
            .iter()
            .any(|step| matches!(step, crate::verification::Operation::Topo { .. }))
    {
        return Err(Error::Unsupported(
            "TOPO reconstruction requires candidate schema 8".into(),
        ));
    }
    if u16::from_le_bytes([bytes[4], bytes[5]]) < 7
        && candidate.plan.steps.iter().any(|step| {
            matches!(
                step,
                crate::verification::Operation::RecombineReference { .. }
            )
        })
    {
        return Err(Error::Unsupported(
            "reference recombination requires candidate schema 7".into(),
        ));
    }
    if u16::from_le_bytes([bytes[4], bytes[5]]) < 6
        && candidate.plan.steps.iter().any(|step| {
            matches!(
                step,
                crate::verification::Operation::DeclareSourceModifications { .. }
                    | crate::verification::Operation::CheckedDigest { .. }
            )
        })
    {
        return Err(Error::Unsupported(
            "tracked modifications and checked digest require candidate schema 6".into(),
        ));
    }
    if u16::from_le_bytes([bytes[4], bytes[5]]) < 5
        && candidate.plan.steps.iter().any(|step| {
            matches!(
                step,
                crate::verification::Operation::Tail3Prime { .. }
                    | crate::verification::Operation::ProfiledOverlap { .. }
            )
        })
    {
        return Err(Error::Unsupported(
            "tailing and profiled-overlap operations require candidate schema5".into(),
        ));
    }
    if u16::from_le_bytes([bytes[4], bytes[5]]) < 4
        && candidate.plan.steps.iter().any(|step| {
            matches!(
                step,
                crate::verification::Operation::RepairEnds { .. }
                    | crate::verification::Operation::AnnealOligos { .. }
            )
        })
    {
        return Err(Error::Unsupported(
            "end-repair and oligo annealing candidate operations require binary schema 4 or later"
                .into(),
        ));
    }
    if u16::from_le_bytes([bytes[4], bytes[5]]) < 3
        && candidate.plan.steps.iter().any(|step| {
            matches!(
                step,
                crate::verification::Operation::PcrWithMismatches { .. }
                    | crate::verification::Operation::GoldenGate { .. }
            )
        })
    {
        return Err(Error::Unsupported(
            "mismatch PCR and Golden Gate candidate operations require binary schema 3 or later"
                .into(),
        ));
    }
    if u16::from_le_bytes([bytes[4], bytes[5]]) == 1
        && candidate.plan.steps.iter().any(|step| {
            matches!(
                step,
                crate::verification::Operation::Digest { .. }
                    | crate::verification::Operation::TreatEnds { .. }
                    | crate::verification::Operation::Ligate { .. }
            )
        })
    {
        return Err(Error::Unsupported(
            "physical candidate operations require binary schema 2 or later".into(),
        ));
    }
    candidate.validate()?;
    Ok(candidate)
}
