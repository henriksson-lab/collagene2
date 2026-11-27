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

pub use crate::run_material_capsule::{
    Capsule, MAX_AGGREGATE_ORIGIN_BYTES, MaterialOriginRef, ReplayEvidence, capsule,
    capsules_preflight, material_capsules_preflight, material_origins_preflight, plan_preflight,
    validate_input_origins, validate_plan_origins, visit_candidate_plans,
};

pub const ENGINE_VERSION: &str = "run-material-1";
pub const INSPECTION_ENGINE_VERSION: &str = "run-material-inspection-1";

/// Documentary origin discovery, or a freshly built executable recipe after an
/// explicit selection. This serializable reply never carries live material.
#[derive(Debug, Clone, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct Inspection {
    pub engine_version: String,
    pub context_revision: String,
    pub event: EventBinding,
    pub origins: Vec<Selection>,
    pub selected: Option<Capsule>,
    pub findings: Vec<String>,
}

pub fn inspection_revision(
    workspace: &Workspace,
    event: &EventBinding,
    selection: &Option<Selection>,
) -> Result<String, Error> {
    event_shape(event)?;
    if let Some(selection) = selection {
        shape(selection)?;
        if &selection.event != event {
            return Err(bad("origin selection names a different saved event"));
        }
    }
    let fingerprint = run_lineage::workspace_fingerprint(workspace)?;
    resolve(workspace, event)?;
    crate::assembly_recipe::hash(&(INSPECTION_ENGINE_VERSION, fingerprint, event, selection))
}

pub fn inspect(
    workspace: &Workspace,
    event: &EventBinding,
    selection: &Option<Selection>,
) -> Result<Inspection, Error> {
    let context_revision = inspection_revision(workspace, event, selection)?;
    let origins = originating_products(workspace, event)?;
    let mut findings = Vec::new();
    let selected = if let Some(selection) = selection {
        if !origins.contains(selection) {
            return Err(bad(
                "explicit origin is not in this exact event's documentary lineage",
            ));
        }
        match capsule(workspace, selection) {
            Ok(origin) => {
                findings.push("The originating recipe was freshly reconstructed. It remains a declaration under its recorded assumptions; reuse does not establish isolation or experimental success. Every downstream preview reconstructs it again.".into());
                Some(origin)
            }
            Err(error) => {
                findings.push(error.to_string().chars().take(4096).collect());
                None
            }
        }
    } else {
        findings.push(if origins.is_empty() {
            "This exact event has no retained name-only cloning-run origin. Imported or edited reference metadata cannot restore physical material.".into()
        } else {
            "Choose an originating run explicitly to reconstruct its saved material. Documentary lineage alone is not a current material assessment.".into()
        });
        None
    };
    let result = Inspection {
        engine_version: INSPECTION_ENGINE_VERSION.into(),
        context_revision,
        event: event.clone(),
        origins,
        selected,
        findings,
    };
    result.validate_binding(&result.context_revision)?;
    Ok(result)
}
impl Inspection {
    /// Readonly shape/binding validation; never supplies physical authority.
    pub fn validate_binding(&self, expected: &str) -> Result<(), Error> {
        if !digest(expected)
            || self.context_revision != expected
            || self.engine_version != INSPECTION_ENGINE_VERSION
            || self.origins.len() > 1000
            || self.findings.len() > 128
            || self.findings.iter().any(|s| s.len() > 16384)
        {
            return Err(bad("origin inspection context or bounds changed"));
        }
        event_shape(&self.event)?;
        let mut seen = BTreeSet::new();
        for origin in &self.origins {
            shape(origin)?;
            if origin.event != self.event
                || !seen.insert((origin.run_id.as_str(), origin.step_index))
            {
                return Err(bad("origin inspection event or choices changed"));
            }
        }
        if let Some(origin) = &self.selected {
            origin.validate_documentary()?;
            if origin.selection.event != self.event || !self.origins.contains(&origin.selection) {
                return Err(bad("prepared origin is not an offered exact-event choice"));
            }
        }
        Ok(())
    }
}

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

/// Append-only choice of saved native material; no serialized choice is physical authority.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub enum MaterialSelection {
    Run(Selection),
    Derived(EventBinding),
}
impl MaterialSelection {
    pub fn event(&self) -> &EventBinding {
        match self {
            Self::Run(s) => &s.event,
            Self::Derived(e) => e,
        }
    }
    pub fn validate(&self) -> Result<(), Error> {
        match self {
            Self::Run(s) => shape(s),
            Self::Derived(e) => event_shape(e),
        }
    }
}
#[derive(Debug, Clone, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub enum MaterialOrigin {
    Run(Capsule),
    Derived(crate::material_derivation::Capsule),
}
impl MaterialOrigin {
    pub fn event(&self) -> &EventBinding {
        match self {
            Self::Run(o) => &o.selection.event,

            Self::Derived(o) => &o.product,
        }
    }
    pub fn selection(&self) -> MaterialSelection {
        match self {
            Self::Run(o) => MaterialSelection::Run(o.selection.clone()),

            Self::Derived(o) => MaterialSelection::Derived(o.product.clone()),
        }
    }
    pub fn validate_documentary(&self) -> Result<(), Error> {
        match self {
            Self::Run(o) => o.validate_documentary(),

            Self::Derived(o) => o.validate_documentary(),
        }
    }
    /// Exact retained archive attestation. No marker, alias or documentary graph
    /// replaces the fresh replay which every downstream verifier performs.
    pub fn validate_for_workspace(
        &self,
        workspace: &Workspace,
        event: &EventBinding,
    ) -> Result<(), Error> {
        event_shape(event)?;
        if self.event() != event {
            return Err(bad("material origin names another exact input event"));
        }
        workspace.validate()?;
        resolve(workspace, event)?;
        match self {
            Self::Run(origin) => {
                crate::run_material_capsule::validate_run_capsule_origins(workspace, origin)
            }
            Self::Derived(origin) => {
                crate::material_derivation::validate_capsule_for_workspace(workspace, origin, event)
            }
        }
    }
}
/// Fresh private state from one just-executed native origin. Documentary
/// capsules still require separate exact-workspace attestation at host boundaries.
#[derive(Debug)]
pub struct FreshMaterial {
    reference: Molecule,
    tracked: TrackedDuplex,
    incomplete: bool,
    findings: Vec<String>,
    findings_total: usize,
    findings_truncated: bool,
}
impl FreshMaterial {
    pub fn reference(&self) -> &Molecule {
        &self.reference
    }
    pub fn tracked(&self) -> &TrackedDuplex {
        &self.tracked
    }
    pub fn incomplete(&self) -> bool {
        self.incomplete
    }
    pub fn findings(&self) -> &[String] {
        &self.findings
    }
    pub fn findings_total(&self) -> usize {
        self.findings_total
    }
    pub fn findings_truncated(&self) -> bool {
        self.findings_truncated
    }
    pub fn into_parts(self) -> (Molecule, TrackedDuplex, bool, Vec<String>, usize, bool) {
        (
            self.reference,
            self.tracked,
            self.incomplete,
            self.findings,
            self.findings_total,
            self.findings_truncated,
        )
    }
}
impl MaterialOrigin {
    /// Execute the unchanged closed native recipe; imported state is never
    /// deserialized into material. This does not attest an unrelated workspace.
    pub fn replay(&self, source: &Molecule) -> Result<FreshMaterial, Error> {
        let (reference, tracked, incomplete, findings, total, prior_truncated) = match self {
            Self::Run(origin) => {
                let fresh = crate::run_material_capsule::replay(origin, source)?;
                let report = fresh.evidence.decode_report()?;
                let mut findings = report.findings;
                findings.extend(
                    report
                        .assumptions
                        .into_iter()
                        .map(|s| format!("Origin assumption: {s}")),
                );
                let total = findings.len();
                (
                    fresh.reference,
                    fresh.tracked,
                    fresh.incomplete,
                    findings,
                    total,
                    false,
                )
            }
            Self::Derived(origin) => {
                crate::material_derivation::replay(origin, source)?.into_parts()
            }
        };
        let mut truncated = prior_truncated || findings.len() > 128;
        let findings = findings
            .into_iter()
            .take(128)
            .map(|s| {
                if s.chars().count() > 2048 {
                    truncated = true;
                    s.chars().take(2048).collect()
                } else {
                    s
                }
            })
            .collect();
        Ok(FreshMaterial {
            reference,
            tracked,
            incomplete,
            findings,
            findings_total: total,
            findings_truncated: truncated,
        })
    }
}
pub const MATERIAL_INSPECTION_ENGINE_VERSION: &str = "material-origin-inspection-1";
#[derive(Debug, Clone, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct MaterialInspection {
    pub engine_version: String,
    pub context_revision: String,
    pub event: EventBinding,
    pub origins: Vec<MaterialSelection>,
    pub selected: Option<MaterialOrigin>,
    pub findings: Vec<String>,
}
pub fn material_inspection_revision(
    workspace: &Workspace,
    event: &EventBinding,
    selection: &Option<MaterialSelection>,
) -> Result<String, Error> {
    event_shape(event)?;
    if let Some(selection) = selection {
        selection.validate()?;
        if selection.event() != event {
            return Err(bad("material selection names a different event"));
        }
    }
    // workspace_fingerprint validates native source/model/attachment bounds first.
    let fingerprint = run_lineage::workspace_fingerprint(workspace)?;
    resolve(workspace, event)?;
    crate::assembly_recipe::hash(&(
        MATERIAL_INSPECTION_ENGINE_VERSION,
        fingerprint,
        event,
        selection,
    ))
}
pub fn inspect_material(
    workspace: &Workspace,
    event: &EventBinding,
    selection: &Option<MaterialSelection>,
) -> Result<MaterialInspection, Error> {
    let context_revision = material_inspection_revision(workspace, event, selection)?;
    let mut origins = originating_products(workspace, event)?
        .into_iter()
        .map(MaterialSelection::Run)
        .collect::<Vec<_>>();

    for asset in workspace
        .attachments
        .values()
        .filter(|a| a.kind == crate::material_derivation::WORKSPACE_ATTACHMENT_KIND)
    {
        let record = crate::material_derivation::load_record(&asset.bytes)?;
        if record.product == *event {
            origins.push(MaterialSelection::Derived(event.clone()));
        }
    }
    if origins.len() > 1000 {
        return Err(Error::TooLarge);
    }
    let mut findings = Vec::new();
    let selected = if let Some(selection) = selection {
        if !origins.contains(selection) {
            return Err(bad(
                "material selection is not an exact retained origin choice",
            ));
        }
        let result = match selection {
            MaterialSelection::Run(s) => capsule(workspace, s).map(MaterialOrigin::Run),
            MaterialSelection::Derived(e) => {
                crate::material_derivation::capsule(workspace, e).map(MaterialOrigin::Derived)
            }
        };
        match result {
            Ok(origin) => {
                findings.push("Selected material was freshly reconstructed from its exact saved native recipe. Origin limitations remain; reuse never proves isolation or experimental success. Downstream verification reconstructs it again.".into());
                Some(origin)
            }
            Err(error) => {
                findings.push(error.to_string().chars().take(4096).collect());
                None
            }
        }
    } else {
        findings.push("Choose a retained native run product, exact saved digest fragment or target-free material derivation. Documentary origin discovery does not restore physical material. Renamed digest/derivation descendants are not inferred.".into());
        None
    };
    let result = MaterialInspection {
        engine_version: MATERIAL_INSPECTION_ENGINE_VERSION.into(),
        context_revision,
        event: event.clone(),
        origins,
        selected,
        findings,
    };
    result.validate_binding(&result.context_revision)?;
    Ok(result)
}
impl MaterialInspection {
    /// Readonly bounded documentary validation; never replay science in a decoder.
    pub fn validate_binding(&self, expected: &str) -> Result<(), Error> {
        if !digest(expected)
            || self.context_revision != expected
            || self.engine_version != MATERIAL_INSPECTION_ENGINE_VERSION
            || self.origins.len() > 1000
            || self.findings.len() > 128
            || self.findings.iter().any(|s| s.len() > 16384)
        {
            return Err(bad("material inspection context or bounds changed"));
        }
        event_shape(&self.event)?;
        for (i, choice) in self.origins.iter().enumerate() {
            choice.validate()?;
            if choice.event() != &self.event || self.origins[..i].contains(choice) {
                return Err(bad("material choices duplicate or name another event"));
            }
        }
        if let Some(origin) = &self.selected {
            origin.validate_documentary()?;
            if origin.event() != &self.event || !self.origins.contains(&origin.selection()) {
                return Err(bad(
                    "prepared material is not an offered exact-event choice",
                ));
            }
        }
        Ok(())
    }
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
pub(crate) fn shape(selection: &Selection) -> Result<(), Error> {
    if !digest(&selection.run_id) || selection.step_index >= 1000 {
        return Err(bad("invalid run, step or exact event binding"));
    }
    event_shape(&selection.event)
}
pub(crate) fn event_shape(event: &EventBinding) -> Result<(), Error> {
    if event.document_id.is_empty()
        || event.document_id.len() > 128
        || !digest(&event.revision)
        || !event.event_id.strip_prefix("event-").is_some_and(digest)
    {
        return Err(bad("invalid exact event binding"));
    }
    Ok(())
}

pub(crate) fn resolve<'a>(
    workspace: &'a Workspace,
    event: &EventBinding,
) -> Result<&'a Molecule, Error> {
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
pub(crate) fn name_only_ancestors<'a>(
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
