//! Bounded existing typed-plan packages preserve an origin's reconstruction scope.
//! Embedded bytes are recipes, never serialized material or verifier authority.
use crate::{
    chemistry::tracked_duplex::{self as modified_duplex, TrackedDuplex},
    documents::{
        transfer::{self, Candidate},
        workspace::Workspace,
    },
    material::{
        lineage::{self as run_lineage, EventBinding, InputEventBinding, ProductBinding},
        outputs as cloning_outputs,
        run::{self as run_material, Selection},
    },
    planning::plan::{DesignPlan, Operation, Source},
    verification::{
        self,
        evidence::{Status, VerificationReport},
    },
};
use collagene_model::{Error, Molecule};
use serde::{Deserialize, Serialize};
use std::{
    collections::{BTreeMap, BTreeSet},
    io::Write,
};

pub const MAX_NESTING: usize = 8;
pub const MAX_CAPSULE_BYTES: usize = 8 * 1024 * 1024;
/// Worker91 policy: one ledger for all root/descendant native origin wrappers.
pub const MAX_AGGREGATE_ORIGIN_BYTES: usize = 32 * 1024 * 1024;
/// One library may declare distinct saved sources for all bounded alternatives.
/// Operation-specific input limits and the shared ancestry/byte ledger still apply.
pub const MAX_DOCUMENTARY_ORIGINS: usize = 256;
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
    models: collagene_model::bounds::ModelBudget,
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
                } | Operation::ReconstructMaterial {
                    template: Source::Step(_),
                    ..
                }
            ) {
                return Err(bad(
                    "saved material reconstruction requires a named saved Input event",
                ));
            }
            for source in crate::planning::plan::Operation::sources(op) {
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
    Run(std::borrow::Cow<'a, Capsule>),
    Derived(std::borrow::Cow<'a, crate::material::derivation::Capsule>),
    DerivedRecord(std::borrow::Cow<'a, crate::material::derivation::Record>),
    TopoGraph(std::borrow::Cow<'a, crate::cloning::topo::graph_history::Capsule>),
    TopoGraphRecord(std::borrow::Cow<'a, crate::cloning::topo::graph_history::Record>),
}
fn material_node(origin: &run_material::MaterialOrigin) -> Pending<'_> {
    match origin {
        run_material::MaterialOrigin::Run(o) => Pending::Run(std::borrow::Cow::Borrowed(o)),

        run_material::MaterialOrigin::Derived(o) => Pending::Derived(std::borrow::Cow::Borrowed(o)),
    }
}
fn material_owned(origin: run_material::MaterialOrigin) -> Pending<'static> {
    match origin {
        run_material::MaterialOrigin::Run(o) => Pending::Run(std::borrow::Cow::Owned(o)),

        run_material::MaterialOrigin::Derived(o) => Pending::Derived(std::borrow::Cow::Owned(o)),
    }
}
fn graph_nodes(request: &crate::cloning::topo::Request) -> impl Iterator<Item = Pending<'_>> {
    request.recipe.materials.iter().filter_map(|m| match m {
        crate::cloning::topo::MaterialRecipe::SavedGraph { origin, .. } => {
            Some(Pending::TopoGraph(std::borrow::Cow::Borrowed(origin)))
        }
        _ => None,
    })
}
fn owned_graph_nodes(
    request: crate::cloning::topo::Request,
) -> impl Iterator<Item = Pending<'static>> {
    request
        .recipe
        .materials
        .into_iter()
        .filter_map(|m| match m {
            crate::cloning::topo::MaterialRecipe::SavedGraph { origin, .. } => {
                Some(Pending::TopoGraph(std::borrow::Cow::Owned(origin)))
            }
            _ => None,
        })
}
fn origins(plan: &DesignPlan) -> Vec<(usize, Pending<'_>)> {
    let mut result = vec![];
    for op in &plan.steps {
        match op {
            Operation::ReconstructRunProduct { origin, .. } => {
                result.push((1, Pending::Run(std::borrow::Cow::Borrowed(origin))))
            }
            Operation::ReconstructMaterial { origin, .. } => {
                result.push((1, Pending::Derived(std::borrow::Cow::Borrowed(origin))))
            }
            Operation::Topo { request, .. } | Operation::CheckedTopo { request, .. } => {
                result.extend(graph_nodes(request).map(|n| (1, n)))
            }
            Operation::TopoActions { request, .. } => {
                result.extend(graph_nodes(&request.preparation).map(|n| (1, n)))
            }
            _ => {}
        }
    }
    result
}
// Count bounded serialized structure without allocating a duplicate capsule. In
// particular an opaque derived record is charged before its bytes are decoded.
fn charge_serialized(budget: &mut Budget, value: &impl Serialize) -> Result<(), Error> {
    charge_serialized_bounded(budget, value, MAX_CAPSULE_BYTES)
}
fn charge_serialized_bounded(
    budget: &mut Budget,
    value: &impl Serialize,
    limit: usize,
) -> Result<(), Error> {
    struct Count {
        bytes: usize,
        limit: usize,
    }
    impl Write for Count {
        fn write(&mut self, bytes: &[u8]) -> std::io::Result<usize> {
            self.bytes = self
                .bytes
                .checked_add(bytes.len())
                .filter(|n| *n <= self.limit)
                .ok_or_else(|| std::io::Error::other("individual origin package exceeds bound"))?;
            Ok(bytes.len())
        }
        fn flush(&mut self) -> std::io::Result<()> {
            Ok(())
        }
    }
    let mut count = Count { bytes: 0, limit };
    serde_json::to_writer(&mut count, value).map_err(|_| Error::TooLarge)?;
    add(&mut budget.bytes, count.bytes, MAX_AGGREGATE_ORIGIN_BYTES)
}

fn shallow_binding(origin: &Capsule, candidate: &Candidate) -> Result<(), Error> {
    if origin.schema_version != 1 || origin.inputs.len() > 1000 || origin.products.len() > 1000 {
        return Err(bad("invalid origin recipe schema or binding counts"));
    }
    run_material::shape(&origin.selection)?;
    run_material::event_shape(&origin.target_event)?;
    let plan = &candidate.plan;
    if candidate.plan_revision != crate::cloning::assembly::recipe::hash(plan)?
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
        let (template, selected) = match op {
            Operation::ReconstructRunProduct { template, origin } => {
                (template, &origin.selection.event)
            }
            Operation::ReconstructMaterial { template, origin } => (template, &origin.product),
            _ => continue,
        };
        let Source::Input(id) = template else {
            return Err(bad("saved material requires a saved Input event"));
        };
        let mut matching = bindings.iter().filter(|b| &b.input_id == id);
        let actual = matching
            .next()
            .ok_or_else(|| bad("saved material source event binding is absent"))?;
        if matching.next().is_some() || &actual.event != selected {
            return Err(bad(
                "saved material source event differs from the explicitly selected saved event",
            ));
        }
    }
    Ok(())
}

/// The tree is decoded iteratively, with aggregate limits charged BEFORE each
/// child envelope is decoded. No recursive Candidate::validate call is used here.

fn walk<'a>(
    roots: Vec<(usize, Pending<'a>)>,
    mut budget: Budget,
    mut visit: impl FnMut(&Capsule, &Candidate) -> Result<(), Error>,
    mut derived_visit: impl FnMut(
        &crate::material::derivation::Record,
        Option<&[crate::material::derivation::SourceSnapshot]>,
    ) -> Result<(), Error>,
    mut graph_visit: impl FnMut(
        &crate::cloning::topo::graph_history::Record,
        Option<&[crate::cloning::topo::Input]>,
    ) -> Result<(), Error>,
) -> Result<(), Error> {
    let mut stack = roots;
    while let Some((depth, node)) = stack.pop() {
        if depth > MAX_NESTING {
            return Err(bad("mixed origin recipe nesting exceeds eight generations"));
        }
        match node {
            Pending::Run(origin) => {
                if origin.candidate.len() > MAX_CAPSULE_BYTES {
                    return Err(Error::TooLarge);
                }
                charge_serialized_bounded(
                    &mut budget,
                    origin.as_ref(),
                    MAX_AGGREGATE_ORIGIN_BYTES,
                )?;
                if origin.schema_version != 1
                    || origin.inputs.len() > 1000
                    || origin.products.len() > 1000
                {
                    return Err(bad("invalid origin recipe schema or binding counts"));
                }
                budget.model(&origin.target)?;
                let candidate = transfer::decode_candidate_shallow(&origin.candidate)?;

                budget.plan(&candidate.plan)?;
                shallow_binding(&origin, &candidate)?;
                visit(&origin, &candidate)?;
                for step in candidate.plan.steps {
                    match step {
                        Operation::ReconstructRunProduct { origin, .. } => {
                            stack.push((depth + 1, Pending::Run(std::borrow::Cow::Owned(origin))))
                        }
                        Operation::ReconstructMaterial { origin, .. } => stack
                            .push((depth + 1, Pending::Derived(std::borrow::Cow::Owned(origin)))),
                        Operation::Topo { request, .. }
                        | Operation::CheckedTopo { request, .. } => {
                            stack.extend(owned_graph_nodes(request).map(|n| (depth + 1, n)))
                        }
                        Operation::TopoActions { request, .. } => stack
                            .extend(owned_graph_nodes(request.preparation).map(|n| (depth + 1, n))),
                        _ => {}
                    }
                }
            }
            Pending::Derived(origin) => {
                // Account the opaque bytes and full wrapper BEFORE decoding it.
                charge_serialized_bounded(
                    &mut budget,
                    origin.as_ref(),
                    crate::material::derivation::MAX_CAPSULE_BYTES,
                )?;
                if origin.schema_version != 2 {
                    return Err(bad("unsupported material capsule schema"));
                }
                run_material::event_shape(&origin.product)?;
                let record =
                    crate::material::derivation::load_record_shallow(&origin.record_bytes)?;
                record.validate_shallow()?;
                crate::material::derivation::validate_source_snapshots(
                    &record.request,
                    &origin.source_snapshots,
                )?;
                if origin.product != record.product {
                    return Err(bad("material capsule source/product snapshot changed"));
                }
                for source in &origin.source_snapshots {
                    budget.model(&source.molecule)?;
                }
                add(&mut budget.steps, 1, MAX_STEPS)?;
                derived_visit(&record, Some(&origin.source_snapshots))?;
                for source in record.request.sources {
                    if let crate::material::derivation::SourceOrigin::Saved(saved) = source.origin {
                        stack.push((depth + 1, material_owned(saved)));
                    }
                }
            }
            Pending::TopoGraph(origin) => {
                charge_serialized_bounded(
                    &mut budget,
                    origin.as_ref(),
                    crate::cloning::topo::graph_history::MAX_CAPSULE_BYTES,
                )?;
                for source in &origin.source_snapshots {
                    budget.model(&source.molecule)?;
                }
                let record = origin.validate_shallow()?;
                add(
                    &mut budget.steps,
                    record
                        .request
                        .request
                        .steps
                        .len()
                        .checked_add(1)
                        .ok_or(Error::TooLarge)?,
                    MAX_STEPS,
                )?;
                graph_visit(&record, Some(&origin.source_snapshots))?;
                for b in record.request.bindings {
                    if let Some(o) = b.origin {
                        stack.push((depth + 1, material_owned(o)));
                    }
                }
                stack.extend(owned_graph_nodes(record.request.request).map(|n| (depth + 1, n)));
            }
            Pending::TopoGraphRecord(record) => {
                charge_serialized_bounded(
                    &mut budget,
                    record.as_ref(),
                    crate::cloning::topo::graph_history::MAX_RECORD_BYTES,
                )?;
                record.validate_documentary()?;
                add(
                    &mut budget.steps,
                    record
                        .request
                        .request
                        .steps
                        .len()
                        .checked_add(1)
                        .ok_or(Error::TooLarge)?,
                    MAX_STEPS,
                )?;
                graph_visit(&record, None)?;
                match record {
                    std::borrow::Cow::Borrowed(record) => {
                        for b in &record.request.bindings {
                            if let Some(o) = &b.origin {
                                stack.push((depth + 1, material_node(o)));
                            }
                        }
                        stack.extend(graph_nodes(&record.request.request).map(|n| (depth + 1, n)));
                    }
                    std::borrow::Cow::Owned(record) => {
                        for b in record.request.bindings {
                            if let Some(o) = b.origin {
                                stack.push((depth + 1, material_owned(o)));
                            }
                        }
                        stack.extend(
                            owned_graph_nodes(record.request.request).map(|n| (depth + 1, n)),
                        );
                    }
                }
            }
            Pending::DerivedRecord(record) => {
                charge_serialized(&mut budget, record.as_ref())?;
                record.validate_shallow()?;
                add(&mut budget.steps, 1, MAX_STEPS)?;
                derived_visit(&record, None)?;
                match record {
                    std::borrow::Cow::Borrowed(record) => {
                        for source in &record.request.sources {
                            if let crate::material::derivation::SourceOrigin::Saved(saved) =
                                &source.origin
                            {
                                stack.push((depth + 1, material_node(saved)));
                            }
                        }
                    }
                    std::borrow::Cow::Owned(record) => {
                        for source in record.request.sources {
                            if let crate::material::derivation::SourceOrigin::Saved(saved) =
                                source.origin
                            {
                                stack.push((depth + 1, material_owned(saved)));
                            }
                        }
                    }
                }
            }
        }
    }
    Ok(())
}
impl Capsule {
    /// Documentary consistency only. This never reconstructs a physical product.
    pub fn validate_documentary(&self) -> Result<(), Error> {
        walk(
            vec![(1, Pending::Run(std::borrow::Cow::Borrowed(self)))],
            Budget::default(),
            |_, _| Ok(()),
            |_, _| Ok(()),
            |_, _| Ok(()),
        )
    }
}
pub fn plan_preflight(plan: &DesignPlan) -> Result<(), Error> {
    if origins(plan).is_empty() {
        return Ok(());
    }
    let mut budget = Budget::default();
    budget.plan(plan)?;
    let nodes = origins(plan);
    walk(nodes, budget, |_, _| Ok(()), |_, _| Ok(()), |_, _| Ok(()))
}
pub fn capsules_preflight<'a>(origins: impl IntoIterator<Item = &'a Capsule>) -> Result<(), Error> {
    walk(
        origins
            .into_iter()
            .map(|o| (1, Pending::Run(std::borrow::Cow::Borrowed(o))))
            .collect(),
        Budget::default(),
        |_, _| Ok(()),
        |_, _| Ok(()),
        |_, _| Ok(()),
    )
}
/// Borrowed roots let compilers reserve their aggregate budget before cloning.
pub enum MaterialOriginRef<'a> {
    Run(&'a Capsule),
    Derived(&'a crate::material::derivation::Capsule),
}
pub fn material_capsules_preflight<'a>(
    origins: impl IntoIterator<Item = MaterialOriginRef<'a>>,
) -> Result<(), Error> {
    walk(
        origins
            .into_iter()
            .map(|origin| {
                (
                    1,
                    match origin {
                        MaterialOriginRef::Run(o) => Pending::Run(std::borrow::Cow::Borrowed(o)),
                        MaterialOriginRef::Derived(o) => {
                            Pending::Derived(std::borrow::Cow::Borrowed(o))
                        }
                    },
                )
            })
            .collect(),
        Budget::default(),
        |_, _| Ok(()),
        |_, _| Ok(()),
        |_, _| Ok(()),
    )
}
/// Scope references only; graph origins never become ordinary duplex MaterialOrigin.
pub enum NativeOriginRef<'a> {
    Run(&'a Capsule),
    Derived(&'a crate::material::derivation::Capsule),
    TopoGraph(&'a crate::cloning::topo::graph_history::Capsule),
}
fn native_nodes<'a>(
    roots: impl IntoIterator<Item = NativeOriginRef<'a>>,
) -> Vec<(usize, Pending<'a>)> {
    roots
        .into_iter()
        .map(|o| {
            (
                1,
                match o {
                    NativeOriginRef::Run(o) => Pending::Run(std::borrow::Cow::Borrowed(o)),
                    NativeOriginRef::Derived(o) => Pending::Derived(std::borrow::Cow::Borrowed(o)),
                    NativeOriginRef::TopoGraph(o) => {
                        Pending::TopoGraph(std::borrow::Cow::Borrowed(o))
                    }
                },
            )
        })
        .collect()
}
pub fn native_origins_preflight<'a>(
    roots: impl IntoIterator<Item = NativeOriginRef<'a>>,
) -> Result<(), Error> {
    walk(
        native_nodes(roots),
        Budget::default(),
        |_, _| Ok(()),
        |_, _| Ok(()),
        |_, _| Ok(()),
    )
}
/// Exact native ownership only. Safe inside Workspace validation; no recursive
/// Workspace::validate, material replay or deserialized-state adoption.
pub fn validate_native_origins<'a>(
    workspace: &Workspace,
    roots: impl IntoIterator<Item = NativeOriginRef<'a>>,
) -> Result<(), Error> {
    workspace_walk(workspace, native_nodes(roots), Budget::default())
}
/// Aggregate both origin types before recipe compilers copy them.
pub fn material_origins_preflight<'a>(
    origins: impl IntoIterator<Item = &'a run_material::MaterialOrigin>,
) -> Result<(), Error> {
    walk(
        origins
            .into_iter()
            .map(|origin| {
                (
                    1,
                    match origin {
                        run_material::MaterialOrigin::Run(o) => {
                            Pending::Run(std::borrow::Cow::Borrowed(o))
                        }

                        run_material::MaterialOrigin::Derived(o) => {
                            Pending::Derived(std::borrow::Cow::Borrowed(o))
                        }
                    },
                )
            })
            .collect(),
        Budget::default(),
        |_, _| Ok(()),
        |_, _| Ok(()),
        |_, _| Ok(()),
    )
}

/// One mixed ledger, including opaque bytes before any derived-record decode.
pub(crate) fn derived_capsule_preflight(
    origin: &crate::material::derivation::Capsule,
) -> Result<(), Error> {
    walk(
        vec![(1, Pending::Derived(std::borrow::Cow::Borrowed(origin)))],
        Budget::default(),
        |_, _| Ok(()),
        |_, _| Ok(()),
        |_, _| Ok(()),
    )
}
pub(crate) fn derived_record_preflight(
    record: &crate::material::derivation::Record,
) -> Result<(), Error> {
    walk(
        vec![(
            1,
            Pending::DerivedRecord(std::borrow::Cow::Borrowed(record)),
        )],
        Budget::default(),
        |_, _| Ok(()),
        |_, _| Ok(()),
        |_, _| Ok(()),
    )
}
/// Documentary dependency enumeration across all native material families.
pub(crate) fn derived_ancestry_walk(
    record: &crate::material::derivation::Record,
    visit: impl FnMut(&Capsule, &Candidate) -> Result<(), Error>,

    derived_visit: impl FnMut(
        &crate::material::derivation::Record,
        Option<&[crate::material::derivation::SourceSnapshot]>,
    ) -> Result<(), Error>,
    graph_visit: impl FnMut(
        &crate::cloning::topo::graph_history::Record,
        Option<&[crate::cloning::topo::Input]>,
    ) -> Result<(), Error>,
) -> Result<(), Error> {
    walk(
        vec![(
            1,
            Pending::DerivedRecord(std::borrow::Cow::Borrowed(record)),
        )],
        Budget::default(),
        visit,
        derived_visit,
        graph_visit,
    )
}

/// Typed documentary dependencies use the SAME mixed-family ledger, not a
/// separate graph parser or per-family ancestry allowance.
pub(crate) fn graph_ancestry_walk(
    record: &crate::cloning::topo::graph_history::Record,
    run_visit: impl FnMut(&Capsule, &Candidate) -> Result<(), Error>,
    derived_visit: impl FnMut(
        &crate::material::derivation::Record,
        Option<&[crate::material::derivation::SourceSnapshot]>,
    ) -> Result<(), Error>,
    graph_visit: impl FnMut(
        &crate::cloning::topo::graph_history::Record,
        Option<&[crate::cloning::topo::Input]>,
    ) -> Result<(), Error>,
) -> Result<(), Error> {
    walk(
        vec![(
            1,
            Pending::TopoGraphRecord(std::borrow::Cow::Borrowed(record)),
        )],
        Budget::default(),
        run_visit,
        derived_visit,
        graph_visit,
    )
}

/// Build a durable declaration only after the actual saved event has been freshly
/// reconstructed. Returned bytes carry instructions, never physical authority.
pub fn capsule(workspace: &Workspace, selection: &Selection) -> Result<Capsule, Error> {
    let mut roots = documentary_origins(
        workspace,
        &[run_material::MaterialSelection::Run(selection.clone())],
        &[],
    )?
    .materials;
    let run_material::MaterialOrigin::Run(origin) = roots.remove(0) else {
        return Err(bad("run capsule kind differs"));
    };
    let _fresh = run_material::reconstruct(workspace, selection)?;
    Ok(origin)
}

/// Build closed documentary recipes for explicitly selected native events.
/// This performs ownership/current record validation and ONE mixed-root budget,
/// before copying candidate bytes, source models or closed child capsules. It
/// never verifies a plan, calls reconstruct/replay, or creates live material.
/// The eventual native consumer must freshly replay the returned origins.
pub struct DocumentaryOrigins {
    pub materials: Vec<run_material::MaterialOrigin>,
    pub graphs: Vec<crate::cloning::topo::graph_history::Capsule>,
}
pub fn documentary_origins(
    workspace: &Workspace,
    selections: &[run_material::MaterialSelection],
    graph_selections: &[crate::cloning::topo::graph_history::Selection],
) -> Result<DocumentaryOrigins, Error> {
    if selections
        .len()
        .checked_add(graph_selections.len())
        .is_none_or(|n| n > MAX_DOCUMENTARY_ORIGINS)
    {
        return Err(Error::TooLarge);
    }
    for selection in selections {
        selection.validate()?;
    }
    workspace.validate()?;
    enum Root {
        Run {
            index: usize,
            record: run_lineage::RunRecord,
        },
        Derived {
            index: usize,
            hash: String,
            record: crate::material::derivation::Record,
        },
    }
    let mut roots = Vec::with_capacity(selections.len());
    let mut loaded = Budget::default();
    // Current native loaders own all record decoding, never a partial/custom reader.
    for (index, selection) in selections.iter().enumerate() {
        match selection {
            run_material::MaterialSelection::Run(selection) => {
                let record = run_lineage::load_run(workspace, &selection.run_id)?;
                // Bound retained decoded roots as they enter the collection too;
                // repeated selection of a large run cannot accumulate uncharged copies
                // before the mixed ancestry phase checks its aggregate allowance.
                charge_serialized_bounded(&mut loaded, &record, MAX_AGGREGATE_ORIGIN_BYTES)?;
                run_material::selected_reference(workspace, selection, &record)?;
                roots.push(Root::Run { index, record });
            }
            run_material::MaterialSelection::Derived(event) => {
                let (hash, record) =
                    crate::material::derivation::documentary_capsule_parts(workspace, event)?;
                charge_serialized_bounded(&mut loaded, &record, MAX_AGGREGATE_ORIGIN_BYTES)?;
                roots.push(Root::Derived {
                    index,
                    hash,
                    record,
                });
            }
        }
    }
    let mut graph_records = Vec::with_capacity(graph_selections.len());
    for selection in graph_selections {
        let record =
            crate::cloning::topo::graph_history::documentary_capsule_parts(workspace, selection)?;
        charge_serialized_bounded(&mut loaded, &record, MAX_AGGREGATE_ORIGIN_BYTES)?;
        graph_records.push(record);
    }
    let mut budget = Budget::default();
    let mut nodes = vec![];
    for (record, selection) in graph_records.iter().zip(graph_selections) {
        let bytes = &workspace.attachments[&selection.record_hash].bytes;
        add(
            &mut budget.bytes,
            bytes.len().div_ceil(3) * 4,
            MAX_AGGREGATE_ORIGIN_BYTES,
        )?;
        add(&mut budget.bytes, 2048, MAX_AGGREGATE_ORIGIN_BYTES)?;
        for binding in &record.request.bindings {
            let source = run_material::resolve(workspace, &binding.event)?;
            charge_serialized_bounded(
                &mut budget,
                &(&binding.role, source),
                MAX_AGGREGATE_ORIGIN_BYTES,
            )?;
            budget.model(source)?;
            if let Some(origin) = &binding.origin {
                nodes.push((2, material_node(origin)));
            }
        }
        add(
            &mut budget.steps,
            record
                .request
                .request
                .steps
                .len()
                .checked_add(1)
                .ok_or(Error::TooLarge)?,
            MAX_STEPS,
        )?;
        nodes.extend(graph_nodes(&record.request.request).map(|n| (2, n)));
    }
    for root in &roots {
        match root {
            Root::Run { index, record } => {
                let target = run_material::resolve(workspace, &record.target)?;
                // Count the actual compact candidate shape without allocating
                // its envelope or cloning a rich source. Conservative reserves
                // cover the current CGP header and outer capsule field names.
                #[derive(Serialize)]
                struct CandidateDeclaration<'a> {
                    plan_revision: &'a str,
                    plan: &'a DesignPlan,
                }
                let candidate_bytes =
                    postcard::experimental::serialized_size(&CandidateDeclaration {
                        plan_revision: &record.plan_revision,
                        plan: &record.plan,
                    })
                    .map_err(|_| Error::TooLarge)?
                    .checked_add(64)
                    .filter(|n| *n <= MAX_CAPSULE_BYTES)
                    .ok_or(Error::TooLarge)?;
                add(
                    &mut budget.bytes,
                    candidate_bytes.div_ceil(3) * 4,
                    MAX_AGGREGATE_ORIGIN_BYTES,
                )?;
                add(&mut budget.bytes, 1024, MAX_AGGREGATE_ORIGIN_BYTES)?;
                charge_serialized_bounded(
                    &mut budget,
                    &(
                        1u16,
                        selections.get(*index),
                        &record.target,
                        &record.inputs,
                        &record.products,
                        target,
                    ),
                    MAX_AGGREGATE_ORIGIN_BYTES,
                )?;
                budget.plan(&record.plan)?;
                budget.model(target)?;
                nodes.extend(
                    origins(&record.plan)
                        .into_iter()
                        .map(|(depth, node)| (depth + 1, node)),
                );
            }
            Root::Derived { record, hash, .. } => {
                let bytes = &workspace.attachments[hash].bytes;
                // Charge every role-bound snapshot and descendant under the one
                // enclosing ledger BEFORE cloning any source or closed child.
                add(
                    &mut budget.bytes,
                    bytes.len().div_ceil(3) * 4,
                    MAX_AGGREGATE_ORIGIN_BYTES,
                )?;
                charge_serialized_bounded(
                    &mut budget,
                    &(2u16, &record.product),
                    MAX_AGGREGATE_ORIGIN_BYTES,
                )?;
                add(&mut budget.bytes, 1024, MAX_AGGREGATE_ORIGIN_BYTES)?;
                add(&mut budget.steps, 1, MAX_STEPS)?;
                for binding in &record.request.sources {
                    let source = run_material::resolve(workspace, &binding.event)?;
                    charge_serialized_bounded(
                        &mut budget,
                        &(binding.role, source),
                        MAX_AGGREGATE_ORIGIN_BYTES,
                    )?;
                    budget.model(source)?;
                    if let crate::material::derivation::SourceOrigin::Saved(origin) =
                        &binding.origin
                    {
                        nodes.push((2, material_node(origin)));
                    }
                }
            }
        }
    }
    workspace_walk(workspace, nodes, budget)?;
    let mut result = Vec::with_capacity(roots.len());
    for root in roots {
        let origin = match root {
            Root::Run { index, record } => {
                let run_material::MaterialSelection::Run(selection) = &selections[index] else {
                    return Err(bad("root kind differs"));
                };
                let target = run_material::resolve(workspace, &record.target)?.clone();
                run_material::MaterialOrigin::Run(Capsule {
                    schema_version: 1,
                    selection: selection.clone(),
                    candidate: transfer::save_candidate(&Candidate::new(record.plan))?,
                    target,
                    target_event: record.target,
                    inputs: record.inputs,
                    products: record.products,
                })
            }
            Root::Derived {
                index,
                hash,
                record,
            } => {
                let run_material::MaterialSelection::Derived(product) = &selections[index] else {
                    return Err(bad("root kind differs"));
                };
                run_material::MaterialOrigin::Derived(crate::material::derivation::Capsule {
                    schema_version: 2,
                    product: product.clone(),
                    record_bytes: workspace.attachments[&hash].bytes.clone(),
                    source_snapshots: record
                        .request
                        .sources
                        .iter()
                        .map(|source| {
                            Ok(crate::material::derivation::SourceSnapshot {
                                role: source.role,
                                molecule: run_material::resolve(workspace, &source.event)?.clone(),
                            })
                        })
                        .collect::<Result<Vec<_>, Error>>()?,
                })
            }
        };
        result.push(origin);
    }
    // Check actual current serialized wrappers too; this cannot expand the
    // pre-copy allowance or grant another ancestry/replay budget.
    let mut graphs = Vec::with_capacity(graph_records.len());
    for (record, selection) in graph_records.into_iter().zip(graph_selections) {
        let mut snapshots = Vec::with_capacity(record.request.bindings.len());
        for b in record.request.bindings {
            snapshots.push(crate::cloning::topo::Input {
                role: b.role,
                molecule: run_material::resolve(workspace, &b.event)?.clone(),
            });
        }
        graphs.push(crate::cloning::topo::graph_history::Capsule {
            schema_version: 2,
            selection: selection.clone(),
            record_bytes: workspace.attachments[&selection.record_hash].bytes.clone(),
            source_snapshots: snapshots,
        });
    }
    native_origins_preflight(
        result
            .iter()
            .map(|o| match o {
                run_material::MaterialOrigin::Run(o) => NativeOriginRef::Run(o),
                run_material::MaterialOrigin::Derived(o) => NativeOriginRef::Derived(o),
            })
            .chain(graphs.iter().map(NativeOriginRef::TopoGraph)),
    )?;
    for origin in &result {
        origin.validate_for_workspace(workspace, origin.event())?;
    }
    validate_native_origins(workspace, graphs.iter().map(NativeOriginRef::TopoGraph))?;
    Ok(DocumentaryOrigins {
        materials: result,
        graphs,
    })
}

/// Workspace ownership only, safe inside workspace/run structural validation.
/// Does not call Workspace::validate/load_run (which would recurse back here),
/// replay science, deserialize material, or confer a current verdict.
fn attest_run(
    workspace: &Workspace,
    records: &BTreeMap<String, run_lineage::RunRecord>,
    origin: &Capsule,
    candidate: &Candidate,
) -> Result<(), Error> {
    let run = records
        .get(&origin.selection.run_id)
        .ok_or_else(|| bad("origin run is not retained in this workspace"))?;
    if run.plan_revision != candidate.plan_revision
        || crate::cloning::assembly::recipe::hash(&run.plan)?
            != crate::cloning::assembly::recipe::hash(&candidate.plan)?
        || run.target != origin.target_event
        || crate::cloning::assembly::recipe::hash(&run.inputs)?
            != crate::cloning::assembly::recipe::hash(&origin.inputs)?
        || crate::cloning::assembly::recipe::hash(&run.products)?
            != crate::cloning::assembly::recipe::hash(&origin.products)?
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
}
fn workspace_walk(
    workspace: &Workspace,
    nodes: Vec<(usize, Pending<'_>)>,
    budget: Budget,
) -> Result<(), Error> {
    // Reject the complete mixed declaration under one ledger before loading or
    // hashing retained owner records. The second pass only performs attestation;
    // it cannot introduce children or grant another ancestry allowance.
    let mut checked = Budget::default();
    checked.bytes = budget.bytes;
    checked.steps = budget.steps;
    checked.bases = budget.bases;
    checked.inputs = budget.inputs;
    let borrowed = nodes
        .iter()
        .map(|(depth, node)| {
            (
                *depth,
                match node {
                    Pending::Run(origin) => {
                        Pending::Run(std::borrow::Cow::Borrowed(origin.as_ref()))
                    }
                    Pending::Derived(origin) => {
                        Pending::Derived(std::borrow::Cow::Borrowed(origin.as_ref()))
                    }
                    Pending::DerivedRecord(record) => {
                        Pending::DerivedRecord(std::borrow::Cow::Borrowed(record.as_ref()))
                    }
                    Pending::TopoGraph(origin) => {
                        Pending::TopoGraph(std::borrow::Cow::Borrowed(origin.as_ref()))
                    }
                    Pending::TopoGraphRecord(record) => {
                        Pending::TopoGraphRecord(std::borrow::Cow::Borrowed(record.as_ref()))
                    }
                },
            )
        })
        .collect();
    walk(
        borrowed,
        checked,
        |_, _| Ok(()),
        |_, _| Ok(()),
        |_, _| Ok(()),
    )?;
    let mut records = BTreeMap::new();
    for attachment in workspace
        .attachments
        .values()
        .filter(|a| a.kind == run_lineage::WORKSPACE_ATTACHMENT_KIND)
    {
        let run = run_lineage::load_run_record(&attachment.bytes)?;
        if records.insert(run.run_id.clone(), run).is_some() {
            return Err(bad("duplicate retained origin run"));
        }
    }
    walk(
        nodes,
        budget,
        |origin, candidate| attest_run(workspace, &records, origin, candidate),
        |record, source| {
            crate::material::derivation::attest_record_shallow(workspace, record, source)
        },
        |record, snapshots| {
            crate::cloning::topo::graph_history::attest_record_shallow(workspace, record, snapshots)
        },
    )
}
pub fn validate_plan_origins(workspace: &Workspace, plan: &DesignPlan) -> Result<(), Error> {
    if origins(plan).is_empty() {
        return Ok(());
    }
    let mut budget = Budget::default();
    budget.plan(plan)?;
    let nodes = origins(plan);
    workspace_walk(workspace, nodes, budget)
}
pub(crate) fn validate_run_capsule_origins(
    workspace: &Workspace,
    origin: &Capsule,
) -> Result<(), Error> {
    workspace_walk(
        workspace,
        vec![(1, Pending::Run(std::borrow::Cow::Borrowed(origin)))],
        Budget::default(),
    )
}

pub(crate) fn validate_derived_request_origins(
    workspace: &Workspace,
    request: &crate::material::derivation::Request,
) -> Result<(), Error> {
    request.validate_shallow()?;
    let mut budget = Budget::default();
    charge_serialized(&mut budget, request)?;
    let mut nodes = vec![];
    for binding in &request.sources {
        let source = run_material::resolve(workspace, &binding.event)?;
        budget.model(source)?;
        if let crate::material::derivation::SourceOrigin::Saved(origin) = &binding.origin {
            if origin.event() != &binding.event {
                return Err(bad(
                    "derived request origin names another exact source event",
                ));
            }
            nodes.push((2, material_node(origin)));
        }
    }
    workspace_walk(workspace, nodes, budget)
}
pub(crate) fn validate_derived_capsule_origins(
    workspace: &Workspace,
    origin: &crate::material::derivation::Capsule,
) -> Result<(), Error> {
    workspace_walk(
        workspace,
        vec![(1, Pending::Derived(std::borrow::Cow::Borrowed(origin)))],
        Budget::default(),
    )
}
/// Attestation can visit current nested plans without another recursive decoder or
/// per-generation budget. The initial bytes are an origin's candidate envelope.
pub fn visit_candidate_plans(
    first: &[u8],
    mut visit: impl FnMut(&Candidate) -> Result<(), Error>,
) -> Result<(), Error> {
    let mut budget = Budget::default();
    if first.len() > MAX_CAPSULE_BYTES {
        return Err(Error::TooLarge);
    }
    add(&mut budget.bytes, first.len(), MAX_AGGREGATE_ORIGIN_BYTES)?;
    let candidate = transfer::decode_candidate_shallow(first)?;
    budget.plan(&candidate.plan)?;
    visit(&candidate)?;
    let nodes = origins(&candidate.plan)
        .into_iter()
        .map(|(depth, node)| (depth + 1, node))
        .collect();
    walk(
        nodes,
        budget,
        |_, candidate| visit(candidate),
        |_, _| Ok(()),
        |_, _| Ok(()),
    )
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
            collagene_model::bounds::model_preflight(&[product])?;
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
    let (report, mut states) =
        verification::engine::verify_with_tracked(&candidate.plan, &origin.target);
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
