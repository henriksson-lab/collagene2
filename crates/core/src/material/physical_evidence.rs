//! Immutable historical processed graphs, separately bound to exact run events.
//! Reading this attachment never restores fresh scientific/material authority.
use crate::{
    documents::workspace::{Workspace, WorkspaceAttachment},
    material::lineage::{
        self as run_lineage, EventBinding, InputEventBinding, ProductBinding, RunRecord,
    },
    planning::plan::{Operation, Source},
    verification::evidence::PhysicalProcessingStepReport,
};
use collagene_model::Error;
use serde::{Deserialize, Serialize};
use std::{
    collections::{BTreeMap, BTreeSet},
    io::Write,
};
pub const KIND: &str = "physical-processing-run-evidence-v1";
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
    pub steps: Vec<PhysicalProcessingStepReport>,
}
fn bad(s: &str) -> Error {
    Error::Format(format!("Historical physical processing evidence: {s}"))
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
        || crate::encoding::bounded_json::hash(&(&e.target, &e.inputs, &e.products))?
            != crate::encoding::bounded_json::hash(&(&r.target, &r.inputs, &r.products))?
    {
        return Err(bad("run identity, engine or exact event bindings disagree"));
    }
    let expected: Vec<_> = r
        .plan
        .steps
        .iter()
        .enumerate()
        .filter_map(|(i, s)| {
            matches!(
                s,
                Operation::TopoActions { .. }
                    | Operation::ReconstructRunProduct { .. }
                    | Operation::ReconstructMaterial { .. }
                    | Operation::Ligate { .. }
                    | Operation::OverlapExtension { .. }
                    | Operation::HomologyAssembly { .. }
                    | Operation::ResectThreePrime { .. }
                    | Operation::HomologyAnnealing { .. }
            )
            .then_some(i)
        })
        .collect();
    if expected.len() != e.steps.len() {
        return Err(bad("processed step evidence is missing or duplicated"));
    }
    for (index, step) in expected.into_iter().zip(&e.steps) {
        if step.step_index != index {
            return Err(bad("processed step order changed"));
        }
        let report = match (&r.plan.steps[index], &step.evidence) {
            (
                Operation::ResectThreePrime { cuts, .. },
                crate::verification::evidence::PhysicalProcessingEvidence::ThreePrimeResection(
                    geometry,
                ),
            ) => {
                geometry.validate_documentary()?;
                if cuts != &geometry.cuts {
                    return Err(bad("3′ resection cuts changed"));
                }
                continue;
            }
            (
                Operation::HomologyAnnealing { policy, .. },
                crate::verification::evidence::PhysicalProcessingEvidence::HomologyAnnealing(
                    geometry,
                ),
            ) => {
                geometry.validate_documentary()?;
                if crate::encoding::bounded_json::hash(policy)?
                    != crate::encoding::bounded_json::hash(&geometry.options)?
                {
                    return Err(bad("3′ annealing policy changed"));
                }
                continue;
            }
            (
                Operation::HomologyAssembly { policy, .. },
                crate::verification::evidence::PhysicalProcessingEvidence::HomologyAssembly(
                    geometry,
                ),
            ) => {
                geometry.validate_documentary()?;
                if crate::encoding::bounded_json::hash(policy)?
                    != crate::encoding::bounded_json::hash(&geometry.options)?
                {
                    return Err(bad("homology assembly policy changed"));
                }
                continue;
            }
            (
                Operation::OverlapExtension { options, .. },
                crate::verification::evidence::PhysicalProcessingEvidence::OverlapExtension(
                    geometry,
                ),
            ) => {
                geometry.validate_documentary()?;
                if options != &geometry.options {
                    return Err(bad("overlap-extension declaration changed"));
                }
                // Run/event bindings retain the exact source and product identities.
                // This documentary geometry never restores physical material;
                // fresh plan replay against those sources is required for reuse.
                continue;
            }
            (
                Operation::ReconstructRunProduct { template, origin },
                crate::verification::evidence::PhysicalProcessingEvidence::RunProduct(replay),
            ) => {
                origin.validate_documentary()?;
                replay.validate_documentary()?;
                if crate::encoding::bounded_json::hash(origin)?
                    != crate::encoding::bounded_json::hash(&replay.origin)?
                {
                    return Err(bad("saved-material origin declaration changed"));
                }
                let revision = match template {
                    Source::Input(id) => r
                        .plan
                        .inputs
                        .iter()
                        .find(|i| &i.id == id)
                        .map(|i| i.revision.as_str()),
                    Source::Step(prior) if *prior < index => r
                        .products
                        .iter()
                        .find(|p| p.step_index == *prior)
                        .map(|p| p.revision.as_str()),
                    Source::Step(_) => None,
                }
                .ok_or_else(|| bad("saved-material source is unavailable or forward"))?;
                if revision != origin.selection.event.revision {
                    return Err(bad("saved-material exact source revision changed"));
                }
                continue;
            }
            (
                Operation::ReconstructMaterial { template, origin },
                crate::verification::evidence::PhysicalProcessingEvidence::DerivedOrigin(replay),
            ) => {
                origin.validate_documentary()?;
                replay.validate_documentary()?;
                if crate::encoding::bounded_json::hash(origin)?
                    != crate::encoding::bounded_json::hash(&replay.origin)?
                {
                    return Err(bad("saved derivation origin declaration changed"));
                }
                let Source::Input(id) = template else {
                    return Err(bad(
                        "saved derivation reconstruction requires an exact original input",
                    ));
                };
                let revision = r
                    .plan
                    .inputs
                    .iter()
                    .find(|i| &i.id == id)
                    .map(|i| i.revision.as_str())
                    .ok_or_else(|| bad("saved derivation source missing"))?;
                let binding = r
                    .inputs
                    .iter()
                    .find(|i| &i.input_id == id)
                    .ok_or_else(|| bad("saved derivation source event missing"))?;
                if revision != origin.product.revision || binding.event != origin.product {
                    return Err(bad("saved derivation exact input event changed"));
                }
                continue;
            }
            (
                Operation::Ligate {
                    fragments,
                    circular,
                    allow_nicked,
                    options,
                },
                crate::verification::evidence::PhysicalProcessingEvidence::LigationOutcomes(report),
            ) => {
                report.validate_for_operation(fragments, *circular, *allow_nicked, options)?;
                if report.report.scientific_engine != r.assessment.engine_version {
                    return Err(bad(
                        "ligation outcome scientific engine differs from owning run",
                    ));
                }
                for binding in &report.bindings {
                    let revision = match &binding.source {
                        Source::Input(id) => r
                            .plan
                            .inputs
                            .iter()
                            .find(|i| &i.id == id)
                            .map(|i| i.revision.as_str()),
                        Source::Step(prior) if *prior < index => r
                            .products
                            .iter()
                            .find(|p| p.step_index == *prior)
                            .map(|p| p.revision.as_str()),
                        Source::Step(_) => None,
                    }
                    .ok_or_else(|| bad("ligation outcome source is unavailable or forward"))?;
                    if revision != binding.source_revision {
                        return Err(bad("ligation outcome exact source revision changed"));
                    }
                }
                continue;
            }
            (
                Operation::TopoActions { .. },
                crate::verification::evidence::PhysicalProcessingEvidence::TopoActions(report),
            ) => report,
            _ => {
                return Err(bad(
                    "historical processing evidence variant does not match its operation",
                ));
            }
        };
        let Operation::TopoActions { sources, request } = &r.plan.steps[index] else {
            return Err(bad("historical processing declaration changed"));
        };
        if crate::encoding::bounded_json::hash(request)?
            != crate::encoding::bounded_json::hash(&report.request)?
        {
            return Err(bad("processed declaration changed"));
        }
        let mut external = vec![];
        for binding in sources {
            let revision = match &binding.source {
                Source::Input(id) => r
                    .plan
                    .inputs
                    .iter()
                    .find(|input| &input.id == id)
                    .map(|input| input.molecule.revision())
                    .ok_or_else(|| bad("missing original input"))?,
                Source::Step(prior) if *prior < index => r
                    .products
                    .iter()
                    .find(|product| product.step_index == *prior)
                    .map(|product| product.revision.clone())
                    .ok_or_else(|| bad("missing exact upstream product event"))?,
                Source::Step(_) => return Err(bad("cyclic or forward upstream product")),
            };
            external.push(crate::cloning::topo::SourceIdentity {
                role: binding.role.clone(),
                revision,
            });
        }
        if crate::encoding::bounded_json::hash(&external)?
            != crate::encoding::bounded_json::hash(&report.input_sources)?
        {
            return Err(bad("owning run external source identities changed"));
        }
        let actual =
            crate::cloning::topo::expanded_source_identities(&external, &request.preparation)?;
        if crate::encoding::bounded_json::hash(&actual)?
            != crate::encoding::bounded_json::hash(&report.sources)?
        {
            return Err(bad(
                "original, upstream or closed graph source identity changed",
            ));
        }
        // Historical snapshots must remain inspectable after engine/model changes.
        // Binding and bounded documentary shape are checked here, not current
        // contour/chemical consistency. Fresh replay alone grants current state.

        if report.engine_version.is_empty()
            || report.engine_version.len() > 128
            || [&report.context_revision, &report.preparation_revision]
                .iter()
                .any(|v| v.len() != 64 || !v.bytes().all(|b| b.is_ascii_hexdigit()))
            || report.bridge_notes.len() > 64
            || report.bridge_notes.iter().any(|n| n.len() > 4096)
            || report.preparation.substrates.strands.len() > 64
            || report.preparation.maps.len() > 4096
        {
            return Err(bad("historical graph documentary shape changed"));
        }
    }
    Ok(())
}
pub(crate) fn attachment(
    r: &RunRecord,
    steps: Vec<PhysicalProcessingStepReport>,
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
        name: format!("{}.physical-processing-graphs.json", r.run_id),
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
        for step in &e.steps {
            if let crate::verification::evidence::PhysicalProcessingEvidence::LigationOutcomes(
                evidence,
            ) = &step.evidence
            {
                let product = r
                    .products
                    .iter()
                    .find(|p| p.step_index == step.step_index)
                    .ok_or_else(|| bad("ligation outcome product event missing"))?;
                let molecule = w
                    .revisions
                    .get(&product.revision)
                    .ok_or_else(|| bad("ligation outcome product revision missing"))?;
                if molecule
                    .metadata
                    .get("source:ligation_outcomes:evidence_hash")
                    != Some(&crate::encoding::bounded_json::hash(evidence)?)
                    || molecule.metadata.get("source:ligation_outcomes:context")
                        != Some(&evidence.report.context_revision)
                {
                    return Err(bad(
                        "ligation outcomes full evidence differs from its recorded product receipt",
                    ));
                }
                let intended = &evidence.report.intended;
                let stored: crate::chemistry::tracked_duplex::ModificationSummary =
                    serde_json::from_str(
                        molecule
                            .metadata
                            .get("tracked_modifications")
                            .ok_or_else(|| bad("ligation product lost its modification receipt"))?,
                    )
                    .map_err(|_| bad("invalid ligation modification receipt"))?;
                let duplex: crate::chemistry::DuplexSummary = serde_json::from_str(
                    molecule
                        .metadata
                        .get("duplex_model")
                        .ok_or_else(|| bad("ligation product lost its strand receipt"))?,
                )
                .map_err(|_| bad("invalid ligation duplex receipt"))?;
                let actual = intended.duplex.summary()?;
                if molecule.sequence != intended.sequence
                    || molecule.topology != intended.topology
                    || crate::encoding::bounded_json::hash(&stored)?
                        != crate::encoding::bounded_json::hash(&intended.modifications)?
                    || crate::encoding::bounded_json::hash(&(
                        duplex.top_length,
                        duplex.bottom_length,
                        duplex.bottom_start,
                        duplex.topology,
                        &duplex.ends,
                        &duplex.nicks,
                        duplex.covalently_closed_in_model,
                    ))? != crate::encoding::bounded_json::hash(&(
                        actual.top_length,
                        actual.bottom_length,
                        actual.bottom_start,
                        actual.topology,
                        &actual.ends,
                        &actual.nicks,
                        actual.covalently_closed_in_model,
                    ))?
                {
                    return Err(bad(
                        "ligation outcomes intended product differs from owning native event",
                    ));
                }
            }
        }
    }
    for r in records.values() {
        if r.plan.steps.iter().any(|s| {
            matches!(
                s,
                Operation::TopoActions { .. }
                    | Operation::ReconstructRunProduct { .. }
                    | Operation::ReconstructMaterial { .. }
                    | Operation::Ligate { .. }
                    | Operation::OverlapExtension { .. }
                    | Operation::HomologyAssembly { .. }
                    | Operation::ResectThreePrime { .. }
                    | Operation::HomologyAnnealing { .. }
            )
        }) && !found.contains(&r.run_id)
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
    find(w, run_id)
}
pub(crate) fn find(w: &Workspace, run_id: &str) -> Result<Option<Evidence>, Error> {
    for a in w.attachments.values().filter(|a| a.kind == KIND) {
        let e = decode(&a.bytes)?;
        if e.run_id == run_id {
            return Ok(Some(e));
        }
    }
    Ok(None)
}
