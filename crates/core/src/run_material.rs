//! Reconstruct saved cloning products from original run declarations and events.
//!
//! Historical product metadata is compared for identity, never deserialized into
//! material. The only physical state returned here comes from a new verifier pass.
use crate::{
    Error, Molecule, cloning_outputs,
    modified_duplex::TrackedDuplex,
    run_lineage::{self, EventBinding},
    verification::{self, Status, VerificationReport},
    workspace::Workspace,
};
use serde::{Deserialize, Serialize};
use std::collections::{BTreeMap, BTreeSet};

pub const ENGINE_VERSION: &str = "run-material-1";

/// An exact saved event and its explicitly selected originating run step.
/// This is a request for reconstruction, not proof of material identity.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct Selection {
    pub run_id: String,
    pub step_index: usize,
    pub event: EventBinding,
}

/// Fresh in-process reconstruction receipt. Intentionally not serializable or
/// constructible by callers. Returning a documentary summary across a worker/MCP
/// boundary does not recreate this receipt on the receiving side.
#[derive(Debug)]
pub struct FreshProduct {
    selection: Selection,
    context_revision: String,
    reference: Molecule,
    tracked: TrackedDuplex,
    report: VerificationReport,
}
impl FreshProduct {
    pub fn selection(&self) -> &Selection {
        &self.selection
    }
    pub fn context_revision(&self) -> &str {
        &self.context_revision
    }
    /// The selected exact saved reference, including its presentation name.
    /// All non-name fields have been checked against the freshly replayed output.
    pub fn reference(&self) -> &Molecule {
        &self.reference
    }
    pub fn tracked(&self) -> &TrackedDuplex {
        &self.tracked
    }
    /// Retain the entire replay, including assumptions and any incomplete-state
    /// findings. A selected component never implies isolation or reaction success.
    pub fn report(&self) -> &VerificationReport {
        &self.report
    }
}

fn bad(message: &str) -> Error {
    Error::Unsupported(format!("saved cloning material: {message}"))
}
fn digest(value: &str) -> bool {
    value.len() == 64
        && value
            .bytes()
            .all(|b| b.is_ascii_digit() || (b'a'..=b'f').contains(&b))
}
fn shape(selection: &Selection) -> Result<(), Error> {
    if !digest(&selection.run_id) || selection.step_index >= 1000 {
        return Err(bad("invalid run, step or exact event binding"));
    }
    event_shape(&selection.event)
}
fn event_shape(event: &EventBinding) -> Result<(), Error> {
    if event.document_id.is_empty()
        || event.document_id.len() > 128
        || !digest(&event.revision)
        || !event.event_id.strip_prefix("event-").is_some_and(digest)
    {
        return Err(bad("invalid exact event binding"));
    }
    Ok(())
}

fn resolve<'a>(workspace: &'a Workspace, event: &EventBinding) -> Result<&'a Molecule, Error> {
    let actual = workspace
        .history
        .iter()
        .find(|entry| entry.id == event.event_id)
        .ok_or_else(|| bad("selected event is unavailable"))?;
    if actual.document_id != event.document_id || actual.revision != event.revision {
        return Err(bad("selected event binding changed"));
    }
    workspace
        .revisions
        .get(&event.revision)
        .ok_or_else(|| bad("selected event revision is unavailable"))
}

/// Follow only an unambiguous name-only lineage. Matching DNA at the endpoints
/// does not erase an intervening edit, metadata change or ambiguous merged parent.
fn name_only_ancestors<'a>(
    workspace: &'a Workspace,
    selected: &'a str,
    reference: &Molecule,
) -> Result<BTreeSet<&'a str>, Error> {
    let by_id: BTreeMap<_, _> = workspace
        .history
        .iter()
        .map(|e| (e.id.as_str(), e))
        .collect();
    let mut seen = BTreeSet::new();
    let mut current = selected;
    loop {
        if !seen.insert(current) {
            return Err(bad("cyclic product lineage"));
        }
        let event = by_id
            .get(current)
            .ok_or_else(|| bad("product ancestor is unavailable"))?;
        let molecule = workspace
            .revisions
            .get(&event.revision)
            .ok_or_else(|| bad("product ancestor revision is unavailable"))?;
        if !cloning_outputs::same_reference(molecule, reference) {
            seen.remove(current);
            return Ok(seen);
        }
        if event.parents.len() != 1 {
            return Ok(seen);
        }
        current = &event.parents[0];
    }
}

/// Find documentary origin choices for an exact selected event without running
/// science. Zero or multiple choices must remain explicit to callers. Even one
/// match is not reusable material until `reconstruct` succeeds; reference-only
/// graph products intentionally appear here but will not gain physical state.
pub fn originating_products(
    workspace: &Workspace,
    event: &EventBinding,
) -> Result<Vec<Selection>, Error> {
    // Event field checks precede potentially large workspace validation/hashing.
    event_shape(event)?;
    let runs = run_lineage::list_runs(workspace)?;
    let reference = resolve(workspace, event)?;
    let ancestors = name_only_ancestors(workspace, &event.event_id, reference)?;
    let mut choices = Vec::new();
    for run in runs {
        for product in run.products {
            if ancestors.contains(product.event_id.as_str()) {
                choices.push(Selection {
                    run_id: run.run_id.clone(),
                    step_index: product.step_index,
                    event: event.clone(),
                });
            }
        }
    }
    Ok(choices)
}

/// Reconstruct the complete originating run against its original independent
/// target. All original exact input events are checked by `load_run`; every
/// recorded output must match this pass. Neither recorded status nor imported
/// duplex/modification summaries grant material authority.
pub fn reconstruct(workspace: &Workspace, selection: &Selection) -> Result<FreshProduct, Error> {
    shape(selection)?;
    let fingerprint = run_lineage::workspace_fingerprint(workspace)?;
    let run = run_lineage::load_run(workspace, &selection.run_id)?;
    let binding = run
        .products
        .get(selection.step_index)
        .ok_or_else(|| bad("selected run step is unavailable"))?;
    let reference = resolve(workspace, &selection.event)?;
    let original = workspace
        .revisions
        .get(&binding.revision)
        .ok_or_else(|| bad("original run product revision is unavailable"))?;
    if !cloning_outputs::same_reference(reference, original)
        || !name_only_ancestors(workspace, &selection.event.event_id, reference)?
            .contains(binding.event_id.as_str())
    {
        return Err(bad(
            "selected event is not an unambiguous name-only descendant of this run product; intervening edits cannot restore material authority",
        ));
    }
    let target = resolve(workspace, &run.target)?;
    let (report, mut states) = verification::verify_with_tracked(&run.plan, target);
    if report.status == Status::Invalid
        || report.products.len() != run.products.len()
        || !report
            .products
            .iter()
            .zip(&run.products)
            .all(|(p, b)| p.revision() == b.revision)
    {
        return Err(bad(
            "fresh reconstruction no longer matches the complete recorded run; replay and review it before reuse",
        ));
    }
    let fresh_reference = report
        .products
        .get(selection.step_index)
        .ok_or_else(|| bad("selected product was not reconstructed"))?;
    if !cloning_outputs::same_reference(reference, fresh_reference) {
        return Err(bad(
            "selected saved reference differs from the freshly reconstructed product",
        ));
    }
    let tracked = states.get_mut(selection.step_index).and_then(Option::take)
        .ok_or_else(|| bad("this run step has no established reusable strand/modification state; a reference projection cannot supply physical material"))?;
    tracked.validate()?;
    let context_revision = crate::assembly_recipe::hash(&(
        ENGINE_VERSION,
        verification::ENGINE_VERSION,
        fingerprint,
        selection,
        &run.plan_revision,
    ))?;
    Ok(FreshProduct {
        selection: selection.clone(),
        context_revision,
        reference: reference.clone(),
        tracked,
        report,
    })
}
