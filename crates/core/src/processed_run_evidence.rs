//! Immutable historical processed graphs, separately bound to exact run events.
//! Reading this attachment never restores fresh scientific/material authority.
use crate::{
    Error, assembly_recipe,
    run_lineage::{self, EventBinding, InputEventBinding, ProductBinding, RunRecord},
    verification::{Operation, ProcessedTopoStepReport, Source},
    workspace::{Workspace, WorkspaceAttachment},
};
use serde::{Deserialize, Serialize};
use std::{
    collections::{BTreeMap, BTreeSet},
    io::Write,
};
pub const KIND: &str = "processed-topo-run-evidence-v1";
pub const MAX_BYTES: usize = 32 * 1024 * 1024;
#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Evidence {
    pub schema_version: u32,
    pub run_id: String,
    pub plan_revision: String,
    pub engine_version: String,
    pub target: EventBinding,
    pub inputs: Vec<InputEventBinding>,
    pub products: Vec<ProductBinding>,
    pub steps: Vec<ProcessedTopoStepReport>,
}
fn bad(s: &str) -> Error {
    Error::Format(format!("Historical processed run evidence: {s}"))
}
struct Bounded(Vec<u8>);
impl Write for Bounded {
    fn write(&mut self, b: &[u8]) -> std::io::Result<usize> {
        if self
            .0
            .len()
            .checked_add(b.len())
            .is_none_or(|n| n > MAX_BYTES)
        {
            return Err(std::io::Error::other("processed run evidence exceeds32MiB"));
        }
        self.0
            .try_reserve(b.len())
            .map_err(|_| std::io::Error::other("processed evidence allocation"))?;
        self.0.extend_from_slice(b);
        Ok(b.len())
    }
    fn flush(&mut self) -> std::io::Result<()> {
        Ok(())
    }
}
fn encode(e: &Evidence) -> Result<Vec<u8>, Error> {
    let mut out = Bounded(vec![]);
    serde_json::to_writer(&mut out, e).map_err(|_| Error::TooLarge)?;
    Ok(out.0)
}
fn decode(bytes: &[u8]) -> Result<Evidence, Error> {
    if bytes.len() > MAX_BYTES {
        return Err(Error::TooLarge);
    };
    serde_json::from_slice(bytes).map_err(|e| bad(&e.to_string()))
}
fn validate(e: &Evidence, r: &RunRecord) -> Result<(), Error> {
    if e.schema_version != 1
        || e.run_id != r.run_id
        || e.plan_revision != r.plan_revision
        || e.engine_version != r.assessment.engine_version
        || e.steps.is_empty()
        || e.steps.len() > 64
        || e.inputs.len() > 1000
        || e.products.len() > 1000
        || assembly_recipe::hash(&(&e.target, &e.inputs, &e.products))?
            != assembly_recipe::hash(&(&r.target, &r.inputs, &r.products))?
    {
        return Err(bad("run identity, engine or exact event bindings disagree"));
    }
    let expected: Vec<_> = r
        .plan
        .steps
        .iter()
        .enumerate()
        .filter_map(|(i, s)| matches!(s, Operation::ProcessedTopo { .. }).then_some(i))
        .collect();
    if expected.len() != e.steps.len() {
        return Err(bad("processed step evidence is missing or duplicated"));
    }
    for (index, step) in expected.into_iter().zip(&e.steps) {
        if step.step_index != index {
            return Err(bad("processed step order changed"));
        }
        let Operation::ProcessedTopo { sources, request } = &r.plan.steps[index] else {
            unreachable!()
        };
        if assembly_recipe::hash(request)? != assembly_recipe::hash(&step.evidence.request)?
            || sources.len() != step.evidence.sources.len()
        {
            return Err(bad("processed declaration changed"));
        }
        for (binding, identity) in sources.iter().zip(&step.evidence.sources) {
            let Source::Input(id) = &binding.source else {
                return Err(bad("unsupported historical processed step source"));
            };
            let input = r
                .plan
                .inputs
                .iter()
                .find(|i| &i.id == id)
                .ok_or_else(|| bad("missing original input"))?;
            if binding.role != identity.role || input.molecule.revision() != identity.revision {
                return Err(bad("original source role/revision changed"));
            }
        }
        // Historical snapshots must remain inspectable after engine/model changes.
        // Binding and bounded documentary shape are checked here, not current
        // contour/chemical consistency. Fresh replay alone grants current state.
        let report = &step.evidence;
        if report.engine_version.is_empty()
            || report.engine_version.len() > 128
            || [&report.context_revision, &report.preparation_revision]
                .iter()
                .any(|v| v.len() != 64 || !v.bytes().all(|b| b.is_ascii_hexdigit()))
            || report.bridge_notes.len() > 64
            || report.bridge_notes.iter().any(|n| n.len() > 4096)
            || report.preparation.substrates.strands.len() > 32
            || report.preparation.maps.len() > 4096
        {
            return Err(bad("historical graph documentary shape changed"));
        }
    }
    Ok(())
}
pub(crate) fn attachment(
    r: &RunRecord,
    steps: Vec<ProcessedTopoStepReport>,
) -> Result<Option<WorkspaceAttachment>, Error> {
    if steps.is_empty() {
        return Ok(None);
    }
    let e = Evidence {
        schema_version: 1,
        run_id: r.run_id.clone(),
        plan_revision: r.plan_revision.clone(),
        engine_version: r.assessment.engine_version.clone(),
        target: r.target.clone(),
        inputs: r.inputs.clone(),
        products: r.products.clone(),
        steps,
    };
    validate(&e, r)?;
    Ok(Some(WorkspaceAttachment {
        kind: KIND.into(),
        name: format!("{}.processed-graphs.json", r.run_id),
        media_type: "application/json".into(),
        bytes: encode(&e)?,
    }))
}
/// Structural archive validation only. No enzyme/PCR/contour calculation runs.
pub(crate) fn validate_workspace(w: &Workspace) -> Result<(), Error> {
    let mut records = BTreeMap::new();
    for a in w
        .attachments
        .values()
        .filter(|a| a.kind == run_lineage::WORKSPACE_ATTACHMENT_KIND)
    {
        let r = run_lineage::load_run_record(&a.bytes)?;
        records.insert(r.run_id.clone(), r);
    }
    let mut found = BTreeSet::new();
    for a in w.attachments.values().filter(|a| a.kind == KIND) {
        if a.media_type != "application/json" {
            return Err(bad("attachment media type changed"));
        }
        let e = decode(&a.bytes)?;
        if !found.insert(e.run_id.clone()) {
            return Err(bad("duplicate processed graph attachment"));
        }
        let r = records
            .get(&e.run_id)
            .ok_or_else(|| bad("graph attachment has no owning run"))?;
        validate(&e, r)?;
    }
    for r in records.values() {
        if r.plan
            .steps
            .iter()
            .any(|s| matches!(s, Operation::ProcessedTopo { .. }))
            && !found.contains(&r.run_id)
        {
            return Err(bad(
                "processed run is missing its historical graph attachment",
            ));
        }
    }
    Ok(())
}
/// Historical document only; callers must replay the literal plan for fresh state.
pub fn load(w: &Workspace, run_id: &str) -> Result<Option<Evidence>, Error> {
    w.validate()?;
    for a in w.attachments.values().filter(|a| a.kind == KIND) {
        let e = decode(&a.bytes)?;
        if e.run_id == run_id {
            return Ok(Some(e));
        }
    }
    Ok(None)
}

/// Worker download is historical raw JSON; the GUI need not parse/serialize it.
#[derive(Clone, Debug, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct Download {
    pub run_id: String,
    pub recorded_engine: String,
    pub step_count: usize,
    pub bytes: Option<Vec<u8>>,
}
pub fn download(w: &Workspace, run_id: &str) -> Result<Download, Error> {
    let r = run_lineage::load_run(w, run_id)?;
    for a in w.attachments.values().filter(|a| a.kind == KIND) {
        let e = decode(&a.bytes)?;
        if e.run_id == run_id {
            validate(&e, &r)?;
            if a.bytes.len() > MAX_BYTES - 4096 {
                return Err(bad(
                    "historical graph exceeds worker download cap; export the full native workspace instead",
                ));
            }
            return Ok(Download {
                run_id: run_id.into(),
                recorded_engine: r.assessment.engine_version,
                step_count: e.steps.len(),
                bytes: Some(a.bytes.clone()),
            });
        }
    }
    Ok(Download {
        run_id: run_id.into(),
        recorded_engine: r.assessment.engine_version,
        step_count: 0,
        bytes: None,
    })
}
