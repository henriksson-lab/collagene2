//! Bounded existing typed-plan packages preserve an origin's reconstruction scope.
//! Embedded bytes are recipes, never serialized material or verifier authority.
use crate::{
    Error, Molecule, cloning_outputs,
    modified_duplex::{self, TrackedDuplex},
    run_lineage::{self, EventBinding, InputEventBinding, ProductBinding},
    run_material::{self, Selection},
    transfer::{self, Candidate},
    verification::{self, DesignPlan, Operation, Source, Status, VerificationReport},
    workspace::Workspace,
};
use serde::{Deserialize, Serialize};
use std::{
    collections::{BTreeMap, BTreeSet},
    io::Write,
};

pub const MAX_NESTING: usize = 8;
pub const MAX_CAPSULE_BYTES: usize = 8 * 1024 * 1024;
pub const MAX_REPORT_BYTES: usize = 32 * 1024 * 1024;
const MAX_STEPS: usize = 1000;
const MAX_INPUT_BASES: usize = 8 * 1024 * 1024;

// Human-readable archives/MCP use base64 instead of JSON integer arrays. The
// worker/candidate binary representation remains a compact length + raw bytes.
// This keeps repeated cloning generations from multiplying byte-array text.
mod compact_bytes {
    use base64::{Engine as _, engine::general_purpose::STANDARD};
    use serde::{
        Deserializer, Serializer,
        de::{self, Visitor},
    };
    use std::fmt;

    pub fn serialize<S: Serializer>(
        bytes: &[u8],
        serializer: S,
        limit: usize,
    ) -> Result<S::Ok, S::Error> {
        if bytes.len() > limit {
            return Err(serde::ser::Error::custom("bounded origin package"));
        }
        if serializer.is_human_readable() {
            serializer.serialize_str(&STANDARD.encode(bytes))
        } else {
            serializer.serialize_bytes(bytes)
        }
    }
    pub fn deserialize<'de, D: Deserializer<'de>>(
        deserializer: D,
        limit: usize,
    ) -> Result<Vec<u8>, D::Error> {
        struct Bytes {
            limit: usize,
            human: bool,
        }
        impl<'de> Visitor<'de> for Bytes {
            type Value = Vec<u8>;
            fn expecting(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
                formatter.write_str("a bounded base64 string or binary byte package")
            }
            fn visit_str<E: de::Error>(self, s: &str) -> Result<Vec<u8>, E> {
                if !self.human || s.len() > self.limit.div_ceil(3) * 4 {
                    return Err(E::custom("oversized origin package"));
                }
                let bytes = STANDARD.decode(s).map_err(E::custom)?;
                if bytes.len() > self.limit {
                    return Err(E::custom("oversized origin package"));
                }
                Ok(bytes)
            }
            fn visit_bytes<E: de::Error>(self, bytes: &[u8]) -> Result<Vec<u8>, E> {
                if self.human || bytes.len() > self.limit {
                    return Err(E::custom("oversized origin package"));
                }
                Ok(bytes.to_vec())
            }
            fn visit_borrowed_bytes<E: de::Error>(self, bytes: &'de [u8]) -> Result<Vec<u8>, E> {
                self.visit_bytes(bytes)
            }
            fn visit_byte_buf<E: de::Error>(self, bytes: Vec<u8>) -> Result<Vec<u8>, E> {
                if self.human || bytes.len() > self.limit {
                    return Err(E::custom("oversized origin package"));
                }
                Ok(bytes)
            }
        }
        let human = deserializer.is_human_readable();
        let visitor = Bytes { limit, human };
        if human {
            deserializer.deserialize_str(visitor)
        } else {
            deserializer.deserialize_bytes(visitor)
        }
    }
}
macro_rules! bounded_bytes {
    ($name:ident,$limit:expr) => {
        mod $name {
            pub fn serialize<S: serde::Serializer>(
                bytes: &[u8],
                serializer: S,
            ) -> Result<S::Ok, S::Error> {
                super::compact_bytes::serialize(bytes, serializer, $limit)
            }
            pub fn deserialize<'de, D: serde::Deserializer<'de>>(
                deserializer: D,
            ) -> Result<Vec<u8>, D::Error> {
                super::compact_bytes::deserialize(deserializer, $limit)
            }
        }
    };
}
bounded_bytes!(candidate_bytes, super::MAX_CAPSULE_BYTES);
bounded_bytes!(report_bytes, super::MAX_REPORT_BYTES);

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct Capsule {
    pub schema_version: u16,
    pub selection: Selection,
    /// Existing CGP candidate envelope, including its schema, checksum and plan.
    #[serde(with = "candidate_bytes")]
    #[cfg_attr(feature = "json-schema", schemars(with = "String"))]
    pub candidate: Vec<u8>,
    pub target: Molecule,
    pub target_event: EventBinding,
    pub inputs: Vec<InputEventBinding>,
    pub products: Vec<ProductBinding>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct ReplayEvidence {
    pub origin: Capsule,
    /// Complete typed VerificationReport encoded as bounded JSON. Child evidence
    /// also uses byte packages, avoiding an unbounded recursive report DTO.
    #[serde(with = "report_bytes")]
    #[cfg_attr(feature = "json-schema", schemars(with = "String"))]
    pub current_report: Vec<u8>,
    pub material_revision: String,
}

fn bad(message: &str) -> Error {
    Error::Unsupported(format!("saved material recipe: {message}"))
}
fn add(total: &mut usize, count: usize, cap: usize) -> Result<(), Error> {
    *total = total
        .checked_add(count)
        .filter(|n| *n <= cap)
        .ok_or(Error::TooLarge)?;
    Ok(())
}
fn digest(s: &str) -> bool {
    s.len() == 64
        && s.bytes()
            .all(|b| b.is_ascii_digit() || (b'a'..=b'f').contains(&b))
}
#[derive(Default)]
struct Budget {
    bytes: usize,
    steps: usize,
    bases: usize,
    inputs: usize,
    models: crate::model::bounds::ModelBudget,
}
impl Budget {
    fn model(&mut self, molecule: &Molecule) -> Result<(), Error> {
        add(&mut self.bases, molecule.sequence.len(), MAX_INPUT_BASES)?;
        self.models.molecule(molecule)?;
        molecule.validate()
    }
    fn plan(&mut self, plan: &DesignPlan) -> Result<(), Error> {
        add(&mut self.steps, plan.steps.len(), MAX_STEPS)?;
        add(&mut self.inputs, plan.inputs.len(), 1000)?;
        if plan.steps.is_empty()
            || plan.name.trim().is_empty()
            || plan.name.len() > 1024
            || !digest(&plan.target_revision)
        {
            return Err(bad("invalid bounded origin plan"));
        }
        let mut ids = BTreeSet::new();
        for input in &plan.inputs {
            if input.id.trim().is_empty() || input.id.len() > 128 || !ids.insert(input.id.as_str())
            {
                return Err(bad("invalid or duplicate origin input role"));
            }
            self.model(&input.molecule)?;
            if input.revision != input.molecule.revision() {
                return Err(bad("origin input revision mismatch"));
            }
        }
        for (index, op) in plan.steps.iter().enumerate() {
            if matches!(
                op,
                Operation::ReconstructRunProduct {
                    template: Source::Step(_),
                    ..
                }
            ) {
                return Err(bad(
                    "saved material reconstruction requires a named saved Input event",
                ));
            }
            for source in verification::operation_sources(op) {
                match source {
                    Source::Input(id) if ids.contains(id.as_str()) => (),
                    Source::Step(prior) if *prior < index => (),
                    _ => return Err(bad("origin has an unavailable or forward step dependency")),
                }
            }
        }
        Ok(())
    }
}
enum Pending<'a> {
    Borrowed(&'a Capsule),
    Owned(Capsule),
}
impl Pending<'_> {
    fn get(&self) -> &Capsule {
        match self {
            Self::Borrowed(c) => c,
            Self::Owned(c) => c,
        }
    }
}
fn origins(plan: &DesignPlan) -> impl Iterator<Item = &Capsule> {
    plan.steps.iter().filter_map(|op| match op {
        Operation::ReconstructRunProduct { origin, .. } => Some(origin),
        _ => None,
    })
}

fn shallow_binding(origin: &Capsule, candidate: &Candidate) -> Result<(), Error> {
    if origin.schema_version != 1 || origin.inputs.len() > 1000 || origin.products.len() > 1000 {
        return Err(bad("invalid origin recipe schema or binding counts"));
    }
    run_material::shape(&origin.selection)?;
    run_material::event_shape(&origin.target_event)?;
    let plan = &candidate.plan;
    if candidate.plan_revision != crate::assembly_recipe::hash(plan)?
        || plan.target_revision != origin.target.revision()
        || origin.target_event.revision != plan.target_revision
        || origin.inputs.len() != plan.inputs.len()
        || origin.products.len() != plan.steps.len()
        || origin.selection.step_index >= origin.products.len()
    {
        return Err(bad("origin plan, target or product binding changed"));
    }
    for (binding, input) in origin.inputs.iter().zip(&plan.inputs) {
        run_material::event_shape(&binding.event)?;
        if binding.input_id != input.id || binding.event.revision != input.revision {
            return Err(bad("exact origin input binding changed"));
        }
    }
    for (index, product) in origin.products.iter().enumerate() {
        run_material::event_shape(&EventBinding {
            document_id: product.document_id.clone(),
            event_id: product.event_id.clone(),
            revision: product.revision.clone(),
        })?;
        if product.step_index != index {
            return Err(bad("origin product order changed"));
        }
    }
    validate_input_origins(plan, &origin.inputs)?;
    Ok(())
}

/// Bind reconstruction declarations to the actual source event chosen by the
/// host. An unrelated equal-content alias is not the selected material parent.
pub fn validate_input_origins(
    plan: &DesignPlan,
    bindings: &[InputEventBinding],
) -> Result<(), Error> {
    for op in &plan.steps {
        if let Operation::ReconstructRunProduct { template, origin } = op {
            let Source::Input(id) = template else {
                return Err(bad("saved material requires a saved Input event"));
            };
            let mut matching = bindings.iter().filter(|b| &b.input_id == id);
            let actual = matching
                .next()
                .ok_or_else(|| bad("saved material source event binding is absent"))?;
            if matching.next().is_some() || actual.event != origin.selection.event {
                return Err(bad(
                    "saved material source event differs from the explicitly selected saved event",
                ));
            }
        }
    }
    Ok(())
}

/// The tree is decoded iteratively, with aggregate limits charged BEFORE each
/// child envelope is decoded. No recursive Candidate::validate call is used here.
fn walk<'a>(
    initial: impl Iterator<Item = &'a Capsule>,
    mut budget: Budget,
    mut visit: impl FnMut(&Capsule, &Candidate) -> Result<(), Error>,
) -> Result<(), Error> {
    let mut stack = initial
        .map(|c| (1, Pending::Borrowed(c)))
        .collect::<Vec<_>>();
    while let Some((depth, node)) = stack.pop() {
        if depth > MAX_NESTING {
            return Err(bad("origin recipe nesting exceeds eight generations"));
        }
        let origin = node.get();
        add(&mut budget.bytes, origin.candidate.len(), MAX_CAPSULE_BYTES)?;
        if origin.schema_version != 1 || origin.inputs.len() > 1000 || origin.products.len() > 1000
        {
            return Err(bad("invalid origin recipe schema or binding counts"));
        }
        budget.model(&origin.target)?;
        let candidate = transfer::decode_candidate_shallow(&origin.candidate)?;
        budget.plan(&candidate.plan)?;
        shallow_binding(origin, &candidate)?;
        visit(origin, &candidate)?;
        for step in candidate.plan.steps {
            if let Operation::ReconstructRunProduct { origin, .. } = step {
                stack.push((depth + 1, Pending::Owned(origin)));
            }
        }
    }
    Ok(())
}

impl Capsule {
    /// Documentary consistency only. This never reconstructs a physical product.
    pub fn validate_documentary(&self) -> Result<(), Error> {
        walk(std::iter::once(self), Budget::default(), |_, _| Ok(()))
    }
}

/// Guard all nested existing plans before outer hashing, cloning or science.
pub fn plan_preflight(plan: &DesignPlan) -> Result<(), Error> {
    if origins(plan).next().is_none() {
        return Ok(());
    }
    let mut budget = Budget::default();
    budget.plan(plan)?;
    walk(origins(plan), budget, |_, _| Ok(()))
}

/// Aggregate borrowed recipe origins before a row compiler clones or executes
/// them. Distinct uses are charged separately because each will be reconstructed.
pub fn capsules_preflight<'a>(origins: impl IntoIterator<Item = &'a Capsule>) -> Result<(), Error> {
    walk(origins.into_iter(), Budget::default(), |_, _| Ok(()))
}

/// Build a durable declaration only after the actual saved event has been freshly
/// reconstructed. Returned bytes carry instructions, never physical authority.
pub fn capsule(workspace: &Workspace, selection: &Selection) -> Result<Capsule, Error> {
    let _fresh = run_material::reconstruct(workspace, selection)?;
    let run = run_lineage::load_run(workspace, &selection.run_id)?;
    let origin = Capsule {
        schema_version: 1,
        selection: selection.clone(),
        candidate: transfer::save_candidate(&Candidate::new(run.plan))?,
        target: run_material::resolve(workspace, &run.target)?.clone(),
        target_event: run.target,
        inputs: run.inputs,
        products: run.products,
    };
    origin.validate_documentary()?;
    Ok(origin)
}

/// Workspace ownership only, safe inside workspace/run structural validation.
/// Does not call Workspace::validate/load_run (which would recurse back here),
/// replay science, deserialize material, or confer a current verdict.
pub fn validate_plan_origins(workspace: &Workspace, plan: &DesignPlan) -> Result<(), Error> {
    if origins(plan).next().is_none() {
        return Ok(());
    }
    let mut records = BTreeMap::new();
    for attachment in workspace
        .attachments
        .values()
        .filter(|a| a.kind == run_lineage::WORKSPACE_ATTACHMENT_KIND)
    {
        let run = run_lineage::load_run_record(&attachment.bytes)?;
        records.insert(run.run_id.clone(), run);
    }
    let mut budget = Budget::default();
    budget.plan(plan)?;
    walk(origins(plan), budget, |origin, candidate| {
        let run = records
            .get(&origin.selection.run_id)
            .ok_or_else(|| bad("origin run is not retained in this workspace"))?;
        if run.plan_revision != candidate.plan_revision
            || crate::assembly_recipe::hash(&run.plan)?
                != crate::assembly_recipe::hash(&candidate.plan)?
            || run.target != origin.target_event
            || crate::assembly_recipe::hash(&run.inputs)?
                != crate::assembly_recipe::hash(&origin.inputs)?
            || crate::assembly_recipe::hash(&run.products)?
                != crate::assembly_recipe::hash(&origin.products)?
            || run_material::resolve(workspace, &run.target)? != &origin.target
        {
            return Err(bad("origin declaration differs from retained run/events"));
        }
        for (binding, input) in origin.inputs.iter().zip(&candidate.plan.inputs) {
            if run_material::resolve(workspace, &binding.event)? != &input.molecule {
                return Err(bad(
                    "original input event no longer owns the declared snapshot",
                ));
            }
        }
        for product in &origin.products {
            run_material::resolve(
                workspace,
                &EventBinding {
                    document_id: product.document_id.clone(),
                    event_id: product.event_id.clone(),
                    revision: product.revision.clone(),
                },
            )?;
        }
        let selected = run_material::resolve(workspace, &origin.selection.event)?;
        let product = &origin.products[origin.selection.step_index];
        let original = workspace
            .revisions
            .get(&product.revision)
            .ok_or_else(|| bad("origin product revision is unavailable"))?;
        if !cloning_outputs::same_reference(selected, original)
            || !run_material::name_only_ancestors(
                workspace,
                &origin.selection.event.event_id,
                selected,
            )?
            .contains(product.event_id.as_str())
        {
            return Err(bad(
                "selected saved event is not a name-only descendant of the origin product",
            ));
        }
        Ok(())
    })
}

struct Bounded(Vec<u8>);
impl Write for Bounded {
    fn write(&mut self, bytes: &[u8]) -> std::io::Result<usize> {
        if bytes.len() > MAX_REPORT_BYTES.saturating_sub(self.0.len()) {
            return Err(std::io::Error::other("saved material report exceeds32MiB"));
        }
        self.0.extend_from_slice(bytes);
        Ok(bytes.len())
    }
    fn flush(&mut self) -> std::io::Result<()> {
        Ok(())
    }
}

impl ReplayEvidence {
    /// Decodes one bounded documentary report. It is never a material receipt.
    pub fn decode_report(&self) -> Result<VerificationReport, Error> {
        if self.current_report.len() > MAX_REPORT_BYTES {
            return Err(Error::TooLarge);
        }
        serde_json::from_slice(&self.current_report)
            .map_err(|_| bad("invalid documentary origin report"))
    }
    pub fn validate_documentary(&self) -> Result<(), Error> {
        self.origin.validate_documentary()?;
        if !digest(&self.material_revision) {
            return Err(bad("invalid material revision summary"));
        }
        let report = self.decode_report()?;
        let mut bases = 0;
        for product in &report.products {
            add(&mut bases, product.sequence.len(), 32_000_000)?;
            crate::model::bounds::model_preflight(&[product])?;
            product.validate()?;
        }
        let candidate = transfer::decode_candidate_shallow(&self.origin.candidate)?;
        if report.plan_revision != candidate.plan_revision
            || report.target_revision != self.origin.target_event.revision
            || report.status == Status::Invalid
            || report.products.len() != self.origin.products.len()
            || !report
                .products
                .iter()
                .zip(&self.origin.products)
                .all(|(p, b)| p.revision() == b.revision)
            || report.engine_version.is_empty()
            || report.engine_version.len() > 128
        {
            return Err(bad("documentary origin report binding changed"));
        }
        Ok(())
    }
}

pub(crate) struct Replayed {
    pub reference: Molecule,
    pub tracked: TrackedDuplex,
    pub evidence: ReplayEvidence,
    pub incomplete: bool,
}
/// Native operation implementation: every use performs a new scoped evaluation.
pub(crate) fn replay(origin: &Capsule, source: &Molecule) -> Result<Replayed, Error> {
    origin.validate_documentary()?;
    if source.revision() != origin.selection.event.revision {
        return Err(bad("selected outer source revision changed"));
    }
    let candidate = transfer::decode_candidate_shallow(&origin.candidate)?;
    let (report, mut states) = verification::verify_with_tracked(&candidate.plan, &origin.target);
    if report.status == Status::Invalid
        || report.products.len() != origin.products.len()
        || !report
            .products
            .iter()
            .zip(&origin.products)
            .all(|(p, b)| p.revision() == b.revision)
    {
        return Err(bad(
            "origin no longer reconstructs its complete recorded product set against its own target",
        ));
    }
    let product = &report.products[origin.selection.step_index];
    if !cloning_outputs::same_reference(source, product) {
        return Err(bad(
            "outer source differs from freshly reconstructed origin material",
        ));
    }
    let tracked = states
        .get_mut(origin.selection.step_index)
        .and_then(Option::take)
        .ok_or_else(|| bad("selected origin has no reusable physical strand/modification state"))?;
    tracked.validate()?;
    let material_revision = modified_duplex::summary(&tracked)?.state_revision;
    let incomplete = report.status == Status::Incomplete;
    let mut writer = Bounded(Vec::new());
    serde_json::to_writer(&mut writer, &report).map_err(|_| Error::TooLarge)?;
    Ok(Replayed {
        reference: source.clone(),
        tracked,
        incomplete,
        evidence: ReplayEvidence {
            origin: origin.clone(),
            current_report: writer.0,
            material_revision,
        },
    })
}
