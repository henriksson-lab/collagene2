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
/// Worker91 policy: one ledger for all root/descendant native origin wrappers.
pub const MAX_AGGREGATE_ORIGIN_BYTES: usize = 32 * 1024 * 1024;
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
                } | Operation::ReconstructDigestFragment {
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
    Run(std::borrow::Cow<'a, Capsule>),
    Digest(std::borrow::Cow<'a, crate::run_digest_history::Capsule>),
    Record(&'a crate::run_digest_history::Record),
    Derived(std::borrow::Cow<'a, crate::material_derivation::Capsule>),
    DerivedRecord(std::borrow::Cow<'a, crate::material_derivation::Record>),
}
fn material_node(origin: &run_material::MaterialOrigin) -> Pending<'_> {
    match origin {
        run_material::MaterialOrigin::Run(o) => Pending::Run(std::borrow::Cow::Borrowed(o)),
        run_material::MaterialOrigin::Digest(o) => Pending::Digest(std::borrow::Cow::Borrowed(o)),
        run_material::MaterialOrigin::Derived(o) => Pending::Derived(std::borrow::Cow::Borrowed(o)),
    }
}
fn material_owned(origin: run_material::MaterialOrigin) -> Pending<'static> {
    match origin {
        run_material::MaterialOrigin::Run(o) => Pending::Run(std::borrow::Cow::Owned(o)),
        run_material::MaterialOrigin::Digest(o) => Pending::Digest(std::borrow::Cow::Owned(o)),
        run_material::MaterialOrigin::Derived(o) => Pending::Derived(std::borrow::Cow::Owned(o)),
    }
}
fn origins(plan: &DesignPlan) -> Vec<(usize, Pending<'_>)> {
    plan.steps
        .iter()
        .filter_map(|op| match op {
            Operation::ReconstructRunProduct { origin, .. } => {
                Some((1, Pending::Run(std::borrow::Cow::Borrowed(origin))))
            }
            Operation::ReconstructDigestFragment { origin, .. } => {
                Some((1, Pending::Digest(std::borrow::Cow::Borrowed(origin))))
            }
            Operation::ReconstructMaterial { origin, .. } => {
                Some((1, Pending::Derived(std::borrow::Cow::Borrowed(origin))))
            }
            _ => None,
        })
        .collect()
}
// Count bounded serialized structure without allocating a duplicate capsule. In
// particular a digest layer is charged before its inner candidate is decoded.
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
        let (template, selected) = match op {
            Operation::ReconstructRunProduct { template, origin } => {
                (template, &origin.selection.event)
            }
            Operation::ReconstructDigestFragment { template, origin } => {
                (template, &origin.record.product)
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
    stack: Vec<(usize, Pending<'a>)>,
    budget: Budget,
    visit: impl FnMut(&Capsule, &Candidate) -> Result<(), Error>,
    digest_visit: impl FnMut(&crate::run_digest_history::Record, Option<&Molecule>) -> Result<(), Error>,
) -> Result<(), Error> {
    walk_enclosed(stack, None, budget, visit, digest_visit)
}
fn walk_mixed<'a>(
    nodes: Vec<(usize, Pending<'a>)>,
    budget: Budget,
    visit: impl FnMut(&Capsule, &Candidate) -> Result<(), Error>,
    digest_visit: impl FnMut(&crate::run_digest_history::Record, Option<&Molecule>) -> Result<(), Error>,
    derived_visit: impl FnMut(
        &crate::material_derivation::Record,
        Option<&Molecule>,
    ) -> Result<(), Error>,
) -> Result<(), Error> {
    walk_enclosed_mixed(nodes, None, budget, visit, digest_visit, derived_visit)
}
fn walk_enclosed<'a>(
    nodes: Vec<(usize, Pending<'a>)>,
    enclosing: Option<u16>,
    budget: Budget,
    visit: impl FnMut(&Capsule, &Candidate) -> Result<(), Error>,
    digest_visit: impl FnMut(&crate::run_digest_history::Record, Option<&Molecule>) -> Result<(), Error>,
) -> Result<(), Error> {
    walk_enclosed_mixed(nodes, enclosing, budget, visit, digest_visit, |_, _| Ok(()))
}
fn walk_enclosed_mixed<'a>(
    roots: Vec<(usize, Pending<'a>)>,
    enclosing_version: Option<u16>,
    mut budget: Budget,
    mut visit: impl FnMut(&Capsule, &Candidate) -> Result<(), Error>,
    mut digest_visit: impl FnMut(
        &crate::run_digest_history::Record,
        Option<&Molecule>,
    ) -> Result<(), Error>,
    mut derived_visit: impl FnMut(
        &crate::material_derivation::Record,
        Option<&Molecule>,
    ) -> Result<(), Error>,
) -> Result<(), Error> {
    let mut stack = roots
        .into_iter()
        .map(|(depth, node)| (depth, enclosing_version, node))
        .collect::<Vec<_>>();
    while let Some((depth, enclosing, node)) = stack.pop() {
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
                let version = u16::from_le_bytes([origin.candidate[4], origin.candidate[5]]);
                let enclosing = Some(enclosing.map_or(version, |prior| prior.min(version)));
                if enclosing.is_some_and(|version| version < 23)
                    && candidate
                        .plan
                        .steps
                        .iter()
                        .any(|step| matches!(step, Operation::LigateWithOutcomes { .. }))
                {
                    return Err(bad(
                        "ligation outcomes require candidate schema23 in every enclosing native capsule",
                    ));
                }
                if enclosing.is_some_and(|version| version < 21)
                    && candidate
                        .plan
                        .steps
                        .iter()
                        .any(|step| matches!(step, Operation::ReconstructDigestFragment { .. }))
                {
                    return Err(bad(
                        "digest reconstruction requires candidate schema21 in every enclosing native capsule",
                    ));
                }
                budget.plan(&candidate.plan)?;
                shallow_binding(&origin, &candidate)?;
                visit(&origin, &candidate)?;
                for step in candidate.plan.steps {
                    match step {
                        Operation::ReconstructRunProduct { origin, .. } => stack.push((
                            depth + 1,
                            enclosing,
                            Pending::Run(std::borrow::Cow::Owned(origin)),
                        )),
                        Operation::ReconstructDigestFragment { origin, .. } => stack.push((
                            depth + 1,
                            enclosing,
                            Pending::Digest(std::borrow::Cow::Owned(origin)),
                        )),
                        Operation::ReconstructMaterial { origin, .. } => stack.push((
                            depth + 1,
                            enclosing,
                            Pending::Derived(std::borrow::Cow::Owned(origin)),
                        )),
                        _ => {}
                    }
                }
            }
            Pending::Digest(origin) => {
                add(&mut budget.steps, 1, MAX_STEPS)?;
                charge_serialized(&mut budget, origin.as_ref())?;
                if origin.schema_version != 1 {
                    return Err(bad("unsupported derived capsule schema"));
                }
                budget.model(&origin.source)?;
                origin.record.validate_shallow()?;
                if origin.source.revision() != origin.record.request.source.revision {
                    return Err(bad("digest capsule source snapshot differs from its event"));
                }
                digest_visit(&origin.record, Some(&origin.source))?;
                match origin {
                    std::borrow::Cow::Borrowed(origin) => stack.push((
                        depth + 1,
                        enclosing,
                        Pending::Run(std::borrow::Cow::Borrowed(&origin.record.request.origin)),
                    )),
                    std::borrow::Cow::Owned(origin) => stack.push((
                        depth + 1,
                        enclosing,
                        Pending::Run(std::borrow::Cow::Owned(origin.record.request.origin)),
                    )),
                }
            }
            Pending::Record(record) => {
                add(&mut budget.steps, 1, MAX_STEPS)?;
                charge_serialized(&mut budget, record)?;
                record.validate_shallow()?;
                digest_visit(record, None)?;
                stack.push((
                    depth + 1,
                    enclosing,
                    Pending::Run(std::borrow::Cow::Borrowed(&record.request.origin)),
                ));
            }
            Pending::Derived(origin) => {
                if enclosing.is_some_and(|version| version < 22) {
                    return Err(bad(
                        "material derivation requires candidate schema22 in every enclosing native capsule",
                    ));
                }
                // Account the opaque bytes and full wrapper BEFORE decoding it.
                charge_serialized_bounded(
                    &mut budget,
                    origin.as_ref(),
                    crate::material_derivation::MAX_CAPSULE_BYTES,
                )?;
                if origin.schema_version != 1 {
                    return Err(bad("unsupported material capsule schema"));
                }
                run_material::event_shape(&origin.product)?;
                budget.model(&origin.source_snapshot)?;
                let record = crate::material_derivation::load_record_shallow(&origin.record_bytes)?;
                record.validate_shallow()?;
                if origin.product != record.product
                    || origin.source_snapshot.revision() != record.request.source.revision
                {
                    return Err(bad("material capsule source/product snapshot changed"));
                }
                add(&mut budget.steps, 1, MAX_STEPS)?;
                derived_visit(&record, Some(&origin.source_snapshot))?;
                if let crate::material_derivation::SourceOrigin::Saved(saved) =
                    record.request.origin
                {
                    stack.push((depth + 1, enclosing, material_owned(saved)));
                }
            }
            Pending::DerivedRecord(record) => {
                if enclosing.is_some_and(|version| version < 22) {
                    return Err(bad(
                        "material derivation requires candidate schema22 in every enclosing native capsule",
                    ));
                }
                charge_serialized(&mut budget, record.as_ref())?;
                record.validate_shallow()?;
                add(&mut budget.steps, 1, MAX_STEPS)?;
                derived_visit(&record, None)?;
                match record {
                    std::borrow::Cow::Borrowed(record) => {
                        if let crate::material_derivation::SourceOrigin::Saved(saved) =
                            &record.request.origin
                        {
                            stack.push((depth + 1, enclosing, material_node(saved)));
                        }
                    }
                    std::borrow::Cow::Owned(record) => {
                        if let crate::material_derivation::SourceOrigin::Saved(saved) =
                            record.request.origin
                        {
                            stack.push((depth + 1, enclosing, material_owned(saved)));
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
        )
    }
}
pub fn plan_preflight(plan: &DesignPlan) -> Result<(), Error> {
    if !plan.steps.iter().any(|step| {
        matches!(
            step,
            Operation::ReconstructRunProduct { .. }
                | Operation::ReconstructDigestFragment { .. }
                | Operation::ReconstructMaterial { .. }
        )
    }) {
        return Ok(());
    }
    let mut budget = Budget::default();
    budget.plan(plan)?;
    let nodes = origins(plan);
    walk(nodes, budget, |_, _| Ok(()), |_, _| Ok(()))
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
    )
}
/// Borrowed roots let compilers reserve their aggregate budget before cloning.
pub enum MaterialOriginRef<'a> {
    Run(&'a Capsule),
    Digest(&'a crate::run_digest_history::Capsule),
    Derived(&'a crate::material_derivation::Capsule),
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
                        MaterialOriginRef::Digest(o) => {
                            Pending::Digest(std::borrow::Cow::Borrowed(o))
                        }
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
    )
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
                        run_material::MaterialOrigin::Digest(o) => {
                            Pending::Digest(std::borrow::Cow::Borrowed(o))
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
    )
}
pub(crate) fn digest_capsule_preflight(
    origin: &crate::run_digest_history::Capsule,
) -> Result<(), Error> {
    walk(
        vec![(1, Pending::Digest(std::borrow::Cow::Borrowed(origin)))],
        Budget::default(),
        |_, _| Ok(()),
        |_, _| Ok(()),
    )
}
pub(crate) fn digest_record_preflight(
    record: &crate::run_digest_history::Record,
) -> Result<(), Error> {
    walk(
        vec![(1, Pending::Record(record))],
        Budget::default(),
        |_, _| Ok(()),
        |_, _| Ok(()),
    )
}
/// A new digest occupies one ancestry layer above its saved native source.
pub(crate) fn digest_origin_walk(
    origin: &Capsule,
    visit: impl FnMut(&Capsule, &Candidate) -> Result<(), Error>,
    digest_visit: impl FnMut(&crate::run_digest_history::Record, Option<&Molecule>) -> Result<(), Error>,
) -> Result<(), Error> {
    walk(
        vec![(2, Pending::Run(std::borrow::Cow::Borrowed(origin)))],
        Budget::default(),
        visit,
        digest_visit,
    )
}

/// One mixed ledger, including opaque bytes before any derived-record decode.
pub(crate) fn derived_capsule_preflight(
    origin: &crate::material_derivation::Capsule,
) -> Result<(), Error> {
    walk(
        vec![(1, Pending::Derived(std::borrow::Cow::Borrowed(origin)))],
        Budget::default(),
        |_, _| Ok(()),
        |_, _| Ok(()),
    )
}
pub(crate) fn derived_record_preflight(
    record: &crate::material_derivation::Record,
) -> Result<(), Error> {
    walk(
        vec![(
            1,
            Pending::DerivedRecord(std::borrow::Cow::Borrowed(record)),
        )],
        Budget::default(),
        |_, _| Ok(()),
        |_, _| Ok(()),
    )
}
/// Documentary dependency enumeration across all native material families.
pub(crate) fn derived_ancestry_walk(
    record: &crate::material_derivation::Record,
    visit: impl FnMut(&Capsule, &Candidate) -> Result<(), Error>,
    digest_visit: impl FnMut(&crate::run_digest_history::Record, Option<&Molecule>) -> Result<(), Error>,
    derived_visit: impl FnMut(
        &crate::material_derivation::Record,
        Option<&Molecule>,
    ) -> Result<(), Error>,
) -> Result<(), Error> {
    walk_mixed(
        vec![(
            1,
            Pending::DerivedRecord(std::borrow::Cow::Borrowed(record)),
        )],
        Budget::default(),
        visit,
        digest_visit,
        derived_visit,
    )
}
/// Existing digest dependencies can opt into explicit derived-record callbacks
/// without starting a separate recursively-reset ledger.
pub(crate) fn digest_origin_walk_mixed(
    origin: &Capsule,
    visit: impl FnMut(&Capsule, &Candidate) -> Result<(), Error>,
    digest_visit: impl FnMut(&crate::run_digest_history::Record, Option<&Molecule>) -> Result<(), Error>,
    derived_visit: impl FnMut(
        &crate::material_derivation::Record,
        Option<&Molecule>,
    ) -> Result<(), Error>,
) -> Result<(), Error> {
    walk_mixed(
        vec![(2, Pending::Run(std::borrow::Cow::Borrowed(origin)))],
        Budget::default(),
        visit,
        digest_visit,
        derived_visit,
    )
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
                    Pending::Digest(origin) => {
                        Pending::Digest(std::borrow::Cow::Borrowed(origin.as_ref()))
                    }
                    Pending::Record(record) => Pending::Record(record),
                    Pending::Derived(origin) => {
                        Pending::Derived(std::borrow::Cow::Borrowed(origin.as_ref()))
                    }
                    Pending::DerivedRecord(record) => {
                        Pending::DerivedRecord(std::borrow::Cow::Borrowed(record.as_ref()))
                    }
                },
            )
        })
        .collect();
    walk_mixed(
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
    walk_mixed(
        nodes,
        budget,
        |origin, candidate| attest_run(workspace, &records, origin, candidate),
        |record, source| {
            crate::run_digest_history::attest_record_shallow(workspace, record, source)
        },
        |record, source| {
            crate::material_derivation::attest_record_shallow(workspace, record, source)
        },
    )
}
pub fn validate_plan_origins(workspace: &Workspace, plan: &DesignPlan) -> Result<(), Error> {
    if !plan.steps.iter().any(|step| {
        matches!(
            step,
            Operation::ReconstructRunProduct { .. }
                | Operation::ReconstructDigestFragment { .. }
                | Operation::ReconstructMaterial { .. }
        )
    }) {
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
pub(crate) fn validate_digest_request_origins(
    workspace: &Workspace,
    origin: &Capsule,
) -> Result<(), Error> {
    workspace_walk(
        workspace,
        vec![(2, Pending::Run(std::borrow::Cow::Borrowed(origin)))],
        Budget::default(),
    )
}
pub(crate) fn validate_digest_capsule_origins(
    workspace: &Workspace,
    origin: &crate::run_digest_history::Capsule,
) -> Result<(), Error> {
    workspace_walk(
        workspace,
        vec![(1, Pending::Digest(std::borrow::Cow::Borrowed(origin)))],
        Budget::default(),
    )
}
pub(crate) fn validate_derived_request_origins(
    workspace: &Workspace,
    request: &crate::material_derivation::Request,
) -> Result<(), Error> {
    let source = run_material::resolve(workspace, &request.source)?;
    let mut budget = Budget::default();
    charge_serialized(&mut budget, request)?;
    budget.model(source)?;
    let nodes = match &request.origin {
        crate::material_derivation::SourceOrigin::OriginalReference => vec![],
        crate::material_derivation::SourceOrigin::Saved(origin) => {
            if origin.event() != &request.source {
                return Err(bad(
                    "derived request origin names another exact source event",
                ));
            }
            vec![(2, material_node(origin))]
        }
    };
    workspace_walk(workspace, nodes, budget)
}
pub(crate) fn validate_derived_capsule_origins(
    workspace: &Workspace,
    origin: &crate::material_derivation::Capsule,
) -> Result<(), Error> {
    workspace_walk(
        workspace,
        vec![(1, Pending::Derived(std::borrow::Cow::Borrowed(origin)))],
        Budget::default(),
    )
}
/// Capability checks can visit nested plans without another recursive decoder or
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
    walk_enclosed(
        nodes,
        Some(u16::from_le_bytes([first[4], first[5]])),
        budget,
        |_, candidate| visit(candidate),
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
