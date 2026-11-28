//! Private receipts of coordinate projections actually performed by native code.
//! These structures are neither executable plans nor deserializable authority.
//! Sequence correspondence never asserts physical atom/strand inheritance.
use crate::planning::plan::Source;
use collagene_model::{
    Error, Molecule,
    annotations::{MappingSegment, Projection},
};
use std::{
    collections::{BTreeMap, BTreeSet},
    sync::{
        Arc, Mutex,
        atomic::{AtomicUsize, Ordering},
    },
};

static NEXT_TRACE_ID: AtomicUsize = AtomicUsize::new(1);

pub(crate) mod native;
#[cfg(test)]
mod tests;

#[derive(Clone, Debug)]
pub(crate) struct Limits {
    pub max_snapshots: usize,
    pub max_projections: usize,
    pub max_segments: usize,
    pub max_feature_ordinals: usize,
    pub max_scanned_bytes: usize,
    pub max_lookup_work: usize,
}
impl Default for Limits {
    fn default() -> Self {
        Self {
            max_snapshots: 4096,
            max_projections: 4096,
            max_segments: 4096,
            max_feature_ordinals: 32768,
            max_scanned_bytes: 64 << 20,
            max_lookup_work: 20_000_000,
        }
    }
}
impl Limits {
    fn validate(&self) -> Result<(), Error> {
        if !(1..=4096).contains(&self.max_snapshots)
            || !(1..=4096).contains(&self.max_projections)
            || !(1..=4096).contains(&self.max_segments)
            || !(1..=32768).contains(&self.max_feature_ordinals)
            || !(1..=(64 << 20)).contains(&self.max_scanned_bytes)
            || !(1..=50_000_000).contains(&self.max_lookup_work)
        {
            return Err(Error::TooLarge);
        }
        Ok(())
    }
}
fn pin(value: &str) -> Result<(), Error> {
    if value.len() != 64
        || !value
            .bytes()
            .all(|byte| byte.is_ascii_digit() || (b'a'..=b'f').contains(&byte))
    {
        return Err(bad("exact current native context revision required"));
    }
    Ok(())
}
fn reference_event_bytes(source: &crate::material::lineage::EventBinding) -> Result<usize, Error> {
    pin(&source.revision)?;
    if source.document_id.is_empty()
        || source.event_id.is_empty()
        || source.document_id.len() > 4096
        || source.event_id.len() > 4096
        || source
            .document_id
            .chars()
            .chain(source.event_id.chars())
            .any(char::is_control)
    {
        return Err(bad("exact bounded reference source event required"));
    }
    Ok(source.document_id.len() + source.event_id.len() + source.revision.len())
}
fn bad(message: &str) -> Error {
    Error::Unsupported(format!("native correspondence trace: {message}"))
}
fn add(value: &mut usize, amount: usize, limit: usize) -> Result<(), Error> {
    *value = value
        .checked_add(amount)
        .filter(|value| *value <= limit)
        .ok_or(Error::TooLarge)?;
    Ok(())
}
#[derive(Default, Debug)]
struct Ledger {
    planned_slots: usize,
    snapshots: usize,
    projections: usize,
    segments: usize,
    feature_ordinals: usize,
    namespace_links: usize,
    scanned_bytes: usize,
}
/// One enclosing native replay job, including all retained-origin replays.
/// Cloning this capability shares counters; it never resets construction limits.
#[derive(Clone, Debug)]
pub(crate) struct TraceBudget {
    ledger: Arc<Mutex<Ledger>>,
    limits: Limits,
}
#[derive(Clone, Copy)]
enum Resource {
    PlannedSlots,
    Snapshots,
    Projections,
    Segments,
    FeatureOrdinals,
    NamespaceLinks,
    ScannedBytes,
}
impl TraceBudget {
    pub fn new(limits: Limits) -> Result<Self, Error> {
        limits.validate()?;
        Ok(Self {
            ledger: Arc::new(Mutex::new(Ledger::default())),
            limits,
        })
    }
    /// Refuse an unrecordable new replay before invoking scientific engines.
    /// This reads the same monotonic counters; the actual builder charges once.
    pub fn preflight_native_job(&self, planned_steps: usize) -> Result<(), Error> {
        let ledger = self
            .ledger
            .lock()
            .map_err(|_| bad("shared native trace ledger is unavailable"))?;
        if planned_steps == 0
            || ledger
                .planned_slots
                .checked_add(planned_steps)
                .is_none_or(|count| count > 1000)
            || ledger.snapshots >= self.limits.max_snapshots
            || ledger.projections >= self.limits.max_projections
            || ledger.segments >= self.limits.max_segments
            || ledger.feature_ordinals >= self.limits.max_feature_ordinals
            || ledger.namespace_links >= 32768
            || ledger.scanned_bytes >= self.limits.max_scanned_bytes
        {
            return Err(Error::TooLarge);
        }
        Ok(())
    }
    /// Reserve one bounded native-model evidence copy before retaining it.
    /// All nested origin owners share this same monotonic ledger.
    pub fn retain_model(&self, molecule: &Molecule) -> Result<(), Error> {
        let remaining = self.remaining_scan()?;
        let bytes = match molecule_bytes(molecule, remaining) {
            Ok(bytes) => bytes,
            Err(error) => {
                self.exhaust_scan();
                return Err(error);
            }
        };
        self.charge(Resource::ScannedBytes, bytes)
    }
    fn charge(&self, resource: Resource, amount: usize) -> Result<(), Error> {
        let mut ledger = self
            .ledger
            .lock()
            .map_err(|_| bad("shared native trace ledger is unavailable"))?;
        let (counter, maximum) = match resource {
            Resource::PlannedSlots => (&mut ledger.planned_slots, 1000),
            Resource::Snapshots => (&mut ledger.snapshots, self.limits.max_snapshots),
            Resource::Projections => (&mut ledger.projections, self.limits.max_projections),
            Resource::Segments => (&mut ledger.segments, self.limits.max_segments),
            Resource::FeatureOrdinals => (
                &mut ledger.feature_ordinals,
                self.limits.max_feature_ordinals,
            ),
            Resource::NamespaceLinks => (&mut ledger.namespace_links, 32768),
            Resource::ScannedBytes => (&mut ledger.scanned_bytes, self.limits.max_scanned_bytes),
        };
        add(counter, amount, maximum)
    }
    fn remaining_scan(&self) -> Result<usize, Error> {
        let ledger = self
            .ledger
            .lock()
            .map_err(|_| bad("shared native trace ledger is unavailable"))?;
        Ok(self
            .limits
            .max_scanned_bytes
            .saturating_sub(ledger.scanned_bytes))
    }
    fn exhaust_scan(&self) {
        if let Ok(mut ledger) = self.ledger.lock() {
            ledger.scanned_bytes = self.limits.max_scanned_bytes;
        }
    }
}
/// A token names one actual native snapshot, including a serial that is never
/// reused after a rejected staged operation. Equal bytes never select a token.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub(crate) struct SnapshotToken {
    trace_id: usize,
    index: usize,
    serial: usize,
}
#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) enum ProjectionKind {
    ReferenceCorrespondence,
    TemplateSynthesis,
    EditedTemplate,
    HomologyCollapse,
}
#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) enum IntroductionKind {
    PrimerTail,
    DeclaredLiteral,
    NewlySynthesizedUnmappedCalls,
}
#[derive(Clone, Debug)]
pub(crate) struct MappedChange {
    pub source_position: usize,
    pub destination_position: usize,
    /// Source-forward call before copying/orientation.
    pub before: char,
    /// Actual destination-forward call.
    pub after: char,
}
#[derive(Clone, Debug)]
pub(crate) struct FeatureProjection {
    pub source_feature_index: usize,
    pub destination_feature_index: usize,
    pub source_annotation_revision: String,
    pub projected_annotation_revision: String,
}
#[derive(Clone, Debug)]
pub(crate) struct ProjectionReceipt {
    source: SnapshotToken,
    pub occurrence: u32,
    pub kind: ProjectionKind,
    pub segments: Vec<MappingSegment>,
    pub features: Vec<FeatureProjection>,
    pub unmapped_source_feature_indices: Vec<usize>,
    pub dropped_source_feature_indices: Vec<usize>,
    pub changes: Vec<MappedChange>,
}
#[derive(Clone, Debug)]
pub(crate) struct Introduction {
    pub start: usize,
    pub end: usize,
    pub kind: IntroductionKind,
}
/// Exact native-validated retained recipe binding. This is documentary identity,
/// not an executable nested plan or a deserializable source of authority.
#[derive(Clone, Debug)]
pub(crate) enum OriginBinding {
    Run {
        selection: crate::material::run::Selection,
        plan_revision: String,
        target_event: crate::material::lineage::EventBinding,
        inputs: Vec<crate::material::lineage::InputEventBinding>,
        products: Vec<crate::material::lineage::ProductBinding>,
    },
    TopoGraph {
        selection: crate::cloning::topo::graph_history::Selection,
        context_revision: String,
        sources: Vec<TopoRootBinding>,
    },
    Material {
        product: crate::material::lineage::EventBinding,
        record_hash: String,
        context_revision: String,
        sources: Vec<crate::material::derivation::SourceBinding>,
    },
}
#[derive(Clone, Debug)]
pub(crate) struct TopoRootBinding {
    pub role: String,
    /// A declared immutable root may lack retained event-level provenance.
    pub event: Option<crate::material::lineage::EventBinding>,
}
#[derive(Clone, Debug)]
pub(crate) struct OriginOccurrence {
    pub step_index: usize,
    pub occurrence: u32,
    pub binding: Arc<OriginBinding>,
}
impl OriginBinding {
    fn admitted_bytes(&self) -> Result<usize, Error> {
        fn text(bytes: &mut usize, value: &str) -> Result<(), Error> {
            if value.is_empty() || value.len() > 4096 || value.chars().any(char::is_control) {
                return Err(bad(
                    "retained origin identity requires bounded nonempty native bindings",
                ));
            }
            add(bytes, value.len(), 1 << 20)
        }
        fn event(
            bytes: &mut usize,
            value: &crate::material::lineage::EventBinding,
        ) -> Result<(), Error> {
            text(bytes, &value.document_id)?;
            text(bytes, &value.event_id)?;
            text(bytes, &value.revision)
        }
        let mut bytes = 0;
        match self {
            Self::Run {
                selection,
                plan_revision,
                target_event,
                inputs,
                products,
            } => {
                if inputs.len() > 1000 || products.len() > 1000 {
                    return Err(Error::TooLarge);
                }
                event(&mut bytes, &selection.event)?;
                text(&mut bytes, &selection.run_id)?;
                text(&mut bytes, plan_revision)?;
                event(&mut bytes, target_event)?;
                for input in inputs {
                    text(&mut bytes, &input.input_id)?;
                    event(&mut bytes, &input.event)?;
                }
                for product in products {
                    if product.step_index >= 1000 {
                        return Err(Error::TooLarge);
                    }
                    text(&mut bytes, &product.document_id)?;
                    text(&mut bytes, &product.event_id)?;
                    text(&mut bytes, &product.revision)?;
                }
            }
            Self::TopoGraph {
                selection,
                context_revision,
                sources,
            } => {
                if sources.len() > 64 {
                    return Err(Error::TooLarge);
                }
                text(&mut bytes, &selection.record_hash)?;
                text(&mut bytes, &selection.state_revision)?;
                text(&mut bytes, context_revision)?;
                for source in sources {
                    text(&mut bytes, &source.role)?;
                    if let Some(binding) = &source.event {
                        event(&mut bytes, binding)?;
                    }
                }
            }
            Self::Material {
                product,
                record_hash,
                context_revision,
                sources,
            } => {
                if sources.len() > 1000 {
                    return Err(Error::TooLarge);
                }
                event(&mut bytes, product)?;
                text(&mut bytes, record_hash)?;
                text(&mut bytes, context_revision)?;
                for source in sources {
                    event(&mut bytes, &source.event)?;
                }
            }
        }
        Ok(bytes)
    }
}
#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) enum SourceHandle {
    Plan(Source),
    Material(crate::material::derivation::SourceRole),
    MaterialProduct,
    ReferenceEvent(crate::material::lineage::EventBinding),
    ReferenceProduct { action_index: usize },
    TopoRole(String),
}
#[derive(Clone, Copy)]
enum TraceKind {
    Plan,
    Material,
    ReferenceEdit,
}
#[derive(Clone, Debug)]
pub(crate) struct SourceIdentity {
    pub scope: Vec<OriginOccurrence>,
    pub source: SourceHandle,
    pub revision: String,
}
#[derive(Debug)]
struct Snapshot {
    token: SnapshotToken,
    scope: Vec<OriginOccurrence>,
    stage: Option<usize>,
    step: Option<usize>,
    identity: Option<SourceIdentity>,
    revision: Option<String>,
    dna_revision: String,
    length: usize,
    feature_count: usize,
    projections: Vec<ProjectionReceipt>,
    introductions: Vec<Introduction>,
}
#[derive(Clone, Debug)]
pub(crate) struct StepTrace {
    pub step_index: usize,
    pub product_revision: String,
    pub product_length: usize,
    /// Direct final-snapshot span coverage only; transitive lookup may be incomplete.
    pub local_coverage_complete: bool,
    final_snapshot: SnapshotToken,
}
#[derive(Clone, Debug)]
pub(crate) enum TraceContext {
    UnboundNative,
    Plan {
        revision: String,
    },
    Material {
        context_revision: String,
    },
    ReferenceEdit {
        source: crate::material::lineage::EventBinding,
        actions_revision: String,
    },
}
#[derive(Debug)]
pub(crate) struct ReconstructionTrace {
    context: TraceContext,
    steps: Vec<Option<StepTrace>>,
    snapshots: Vec<Snapshot>,
    limits: Limits,
    budget: TraceBudget,
}
pub(crate) struct TraceBuilder {
    trace: ReconstructionTrace,
    budget: TraceBudget,
    trace_id: usize,
    kind: TraceKind,
    next_serial: usize,
    next_step: usize,
    input_ids: BTreeSet<String>,
}
/// One actual native projection. The destination is an explicit draft token;
/// sources must already be sealed. Feature ordinals come from the canonical
/// Projection owner, never a label, location guess or sequence search.
pub(crate) struct ProjectionCapture<'a> {
    pub source: SnapshotToken,
    pub source_molecule: &'a Molecule,
    pub destination: SnapshotToken,
    pub destination_molecule: &'a Molecule,
    pub occurrence: u32,
    pub kind: ProjectionKind,
    pub segments: &'a [MappingSegment],
    pub projection: &'a Projection,
    pub destination_feature_indices: &'a [usize],
    pub changes: &'a [MappedChange],
}
pub(crate) struct StepBuilder<'a> {
    builder: &'a mut TraceBuilder,
    step: usize,
    snapshot_start: usize,
    committed: bool,
}
fn molecule_bytes(molecule: &Molecule, maximum: usize) -> Result<usize, Error> {
    collagene_model::bounds::model_preflight(&[molecule])?;
    molecule.validate()?;
    crate::verification::evidence::evidence_bytes(molecule, maximum)
}
fn dna_revision(molecule: &Molecule) -> Result<String, Error> {
    crate::encoding::bounded_json::hash(&(
        molecule.sequence.alphabet(),
        molecule.topology,
        &molecule.sequence,
    ))
}
fn feature_revision(
    feature: &collagene_model::Feature,
    budget: &TraceBudget,
) -> Result<String, Error> {
    collagene_model::bounds::feature_preflight(&[std::slice::from_ref(feature)])?;
    let remaining = budget.remaining_scan()?.min(256 << 10);
    let bytes = match crate::verification::evidence::evidence_bytes(feature, remaining) {
        Ok(bytes) => bytes,
        Err(error) => {
            budget.exhaust_scan();
            return Err(error);
        }
    };
    budget.charge(Resource::ScannedBytes, bytes)?;
    crate::encoding::bounded_json::hash(feature)
}
impl TraceBuilder {
    pub fn new(planned_steps: usize, limits: Limits) -> Result<Self, Error> {
        Self::with_budget(planned_steps, TraceBudget::new(limits)?)
    }
    pub fn with_budget(planned_steps: usize, budget: TraceBudget) -> Result<Self, Error> {
        budget.charge(Resource::PlannedSlots, planned_steps)?;
        let limits = budget.limits.clone();
        Ok(Self {
            trace: ReconstructionTrace {
                context: TraceContext::UnboundNative,
                steps: vec![None; planned_steps],
                snapshots: vec![],
                limits,
                budget: budget.clone(),
            },
            budget,
            trace_id: NEXT_TRACE_ID
                .fetch_update(Ordering::Relaxed, Ordering::Relaxed, |id| id.checked_add(1))
                .map_err(|_| Error::TooLarge)?,
            kind: TraceKind::Plan,
            next_serial: 0,
            next_step: 0,
            input_ids: BTreeSet::new(),
        })
    }
    fn snapshot(&self, token: SnapshotToken) -> Result<&Snapshot, Error> {
        self.trace
            .snapshots
            .get(token.index)
            .filter(|snapshot| snapshot.token == token)
            .ok_or_else(|| bad("snapshot token is stale, missing or from a discarded operation"))
    }
    fn charge_model(&mut self, molecule: &Molecule) -> Result<(), Error> {
        self.budget.retain_model(molecule)
    }
    fn allocate(
        &mut self,
        molecule: &Molecule,
        step: Option<usize>,
    ) -> Result<SnapshotToken, Error> {
        self.charge_model(molecule)?;
        self.budget.charge(Resource::Snapshots, 1)?;
        self.next_serial = self.next_serial.checked_add(1).ok_or(Error::TooLarge)?;
        let token = SnapshotToken {
            trace_id: self.trace_id,
            index: self.trace.snapshots.len(),
            serial: self.next_serial,
        };
        self.trace.snapshots.push(Snapshot {
            token,
            scope: vec![],
            stage: step,
            step,
            identity: None,
            revision: None,
            dna_revision: dna_revision(molecule)?,
            length: molecule.sequence.len(),
            feature_count: 0,
            projections: vec![],
            introductions: vec![],
        });
        Ok(token)
    }
    pub fn register_input(
        &mut self,
        id: &str,
        molecule: &Molecule,
    ) -> Result<SnapshotToken, Error> {
        if !matches!(self.kind, TraceKind::Plan)
            || self.next_step != 0
            || id.trim().is_empty()
            || id.len() > 128
            || id.chars().any(char::is_control)
            || self.input_ids.contains(id)
        {
            return Err(bad(
                "distinct bounded original input role required before execution",
            ));
        }
        let token = self.allocate(molecule, None)?;
        let revision = molecule.revision();
        let snapshot = &mut self.trace.snapshots[token.index];
        snapshot.feature_count = molecule.features.len();
        snapshot.revision = Some(revision.clone());
        snapshot.identity = Some(SourceIdentity {
            scope: vec![],
            source: SourceHandle::Plan(Source::Input(id.into())),
            revision,
        });
        self.input_ids.insert(id.into());
        Ok(token)
    }
    /// One actual canonical material derivation, without constructing a plan.
    pub fn for_material(budget: TraceBudget, context_revision: &str) -> Result<Self, Error> {
        pin(context_revision)?;
        let mut builder = Self::with_budget(1, budget)?;
        builder.kind = TraceKind::Material;
        builder.trace.context = TraceContext::Material {
            context_revision: context_revision.into(),
        };
        Ok(builder)
    }
    pub fn bind_plan(&mut self, revision: &str) -> Result<(), Error> {
        pin(revision)?;
        if self.next_step != 0
            || !matches!(self.kind, TraceKind::Plan)
            || !matches!(self.trace.context, TraceContext::UnboundNative)
        {
            return Err(bad(
                "plan identity must be pinned once before actual native execution",
            ));
        }
        self.trace.context = TraceContext::Plan {
            revision: revision.into(),
        };
        Ok(())
    }
    pub fn for_reference_edits(
        budget: TraceBudget,
        source: &crate::material::lineage::EventBinding,
        actions_revision: &str,
        planned_steps: usize,
    ) -> Result<Self, Error> {
        pin(actions_revision)?;
        budget.charge(
            Resource::ScannedBytes,
            reference_event_bytes(source)? + actions_revision.len(),
        )?;
        let mut builder = Self::with_budget(planned_steps, budget)?;
        builder.kind = TraceKind::ReferenceEdit;
        builder.trace.context = TraceContext::ReferenceEdit {
            source: source.clone(),
            actions_revision: actions_revision.into(),
        };
        Ok(builder)
    }
    pub fn register_reference_source(
        &mut self,
        binding: &crate::material::lineage::EventBinding,
        molecule: &Molecule,
    ) -> Result<SnapshotToken, Error> {
        if self.next_step != 0
            || !matches!(&self.trace.context, TraceContext::ReferenceEdit { source, .. } if source == binding)
            || !self.trace.snapshots.is_empty()
        {
            return Err(bad(
                "the exact bound native reference source can be registered once before editing",
            ));
        }
        self.budget
            .charge(Resource::ScannedBytes, reference_event_bytes(binding)?)?;
        self.charge_model(molecule)?;
        let revision = molecule.revision();
        if revision != binding.revision {
            return Err(bad(
                "native reference event revision differs from its actual source",
            ));
        }
        let token = self.allocate(molecule, None)?;
        let snapshot = &mut self.trace.snapshots[token.index];
        snapshot.feature_count = molecule.features.len();
        snapshot.revision = Some(revision.clone());
        snapshot.identity = Some(SourceIdentity {
            scope: vec![],
            source: SourceHandle::ReferenceEvent(binding.clone()),
            revision,
        });
        Ok(token)
    }
    /// Register actual immutable donor events before action execution. Equal
    /// sequence never aliases events; exact repeated event bindings reuse the
    /// same original source, while copy occurrences stay on projection edges.
    pub fn register_reference_donor(
        &mut self,
        binding: &crate::material::lineage::EventBinding,
        molecule: &Molecule,
    ) -> Result<SnapshotToken, Error> {
        if self.next_step != 0
            || !matches!(self.kind, TraceKind::ReferenceEdit)
            || self.trace.snapshots.is_empty()
        {
            return Err(bad(
                "reference donors require the initial bound source before action execution",
            ));
        }
        let bytes = reference_event_bytes(binding)?;
        self.charge_model(molecule)?;
        let revision = molecule.revision();
        if revision != binding.revision {
            return Err(bad(
                "native donor event revision differs from its actual immutable source",
            ));
        }
        if let Some(snapshot) = self.trace.snapshots.iter().find(|snapshot| snapshot.identity.as_ref().is_some_and(|identity| matches!(&identity.source, SourceHandle::ReferenceEvent(event) if event == binding))) {
            return Ok(snapshot.token);
        }
        self.budget.charge(Resource::ScannedBytes, bytes)?;
        let token = self.allocate(molecule, None)?;
        let snapshot = &mut self.trace.snapshots[token.index];
        snapshot.feature_count = molecule.features.len();
        snapshot.revision = Some(revision.clone());
        snapshot.identity = Some(SourceIdentity {
            scope: vec![],
            source: SourceHandle::ReferenceEvent(binding.clone()),
            revision,
        });
        Ok(token)
    }
    pub fn register_material_source(
        &mut self,
        role: crate::material::derivation::SourceRole,
        molecule: &Molecule,
    ) -> Result<SnapshotToken, Error> {
        if !matches!(self.kind, TraceKind::Material)
            || self.next_step != 0
            || self.trace.snapshots.iter().any(|snapshot| {
                snapshot
                    .identity
                    .as_ref()
                    .is_some_and(|identity| identity.source == SourceHandle::Material(role))
            })
        {
            return Err(bad(
                "actual distinct material role required before native derivation capture",
            ));
        }
        let token = self.allocate(molecule, None)?;
        let revision = molecule.revision();
        let snapshot = &mut self.trace.snapshots[token.index];
        snapshot.feature_count = molecule.features.len();
        snapshot.revision = Some(revision.clone());
        snapshot.identity = Some(SourceIdentity {
            scope: vec![],
            source: SourceHandle::Material(role),
            revision,
        });
        Ok(token)
    }
    pub fn begin_step(&mut self, step: usize) -> Result<StepBuilder<'_>, Error> {
        if step >= self.trace.steps.len() || step < self.next_step {
            return Err(bad(
                "step trace must use its original forward execution index once",
            ));
        }
        self.next_step = step.checked_add(1).ok_or(Error::TooLarge)?;
        let snapshot_start = self.trace.snapshots.len();
        Ok(StepBuilder {
            builder: self,
            step,
            snapshot_start,
            committed: false,
        })
    }
    pub fn discard_step(&mut self, index: usize) -> Result<(), Error> {
        if index.checked_add(1) != Some(self.next_step) {
            return Err(bad("only the latest attempted step can be discarded"));
        }
        let first = self
            .trace
            .snapshots
            .iter()
            .position(|snapshot| snapshot.stage == Some(index));
        if let Some(first) = first {
            if self.trace.snapshots[first..]
                .iter()
                .any(|snapshot| snapshot.stage != Some(index))
            {
                return Err(bad("discard would invalidate a later committed snapshot"));
            }
            self.trace.snapshots.truncate(first);
        }
        if let Some(slot) = self.trace.steps.get_mut(index) {
            *slot = None;
        }
        // Charges and serials survive rollback; surviving indices never move.
        Ok(())
    }
    pub fn finish(self) -> ReconstructionTrace {
        self.trace
    }
}
fn segment_end(segment: &MappingSegment) -> Result<usize, Error> {
    segment
        .destination_start
        .checked_add(
            segment
                .source_end
                .checked_sub(segment.source_start)
                .ok_or(Error::TooLarge)?,
        )
        .ok_or(Error::TooLarge)
}
impl StepBuilder<'_> {
    pub fn destination(&mut self, molecule: &Molecule) -> Result<SnapshotToken, Error> {
        self.builder.allocate(molecule, Some(self.step))
    }
    pub fn project(&mut self, capture: ProjectionCapture<'_>) -> Result<(), Error> {
        let ProjectionCapture {
            source,
            source_molecule,
            destination,
            destination_molecule,
            occurrence,
            kind,
            segments,
            projection,
            destination_feature_indices,
            changes,
        } = capture;
        if occurrence >= 4096 || segments.len() > 4096 || changes.len() > 32768 {
            return Err(Error::TooLarge);
        }
        self.builder.snapshot(source)?;
        self.builder.snapshot(destination)?;
        self.builder.charge_model(source_molecule)?;
        self.builder.charge_model(destination_molecule)?;
        let from = self.builder.snapshot(source)?;
        let to = self.builder.snapshot(destination)?;
        if from.revision.is_none()
            || to.revision.is_some()
            || to.step != Some(self.step)
            || source == destination
            || from.scope.is_empty() && from.step.is_some_and(|step| step > self.step)
            || destination_feature_indices.len() != projection.features.len()
            || projection.source_feature_indices.len() != projection.features.len()
            || projection.unmapped_source_feature_indices.len()
                != projection.unmapped_features.len()
            || projection.dropped_source_feature_indices.len() != projection.dropped_features
        {
            return Err(bad(
                "actual sealed source, current draft destination and aligned feature ordinals required",
            ));
        }
        if from.revision.as_deref() != Some(source_molecule.revision().as_str())
            || dna_revision(destination_molecule)? != to.dna_revision
            || from.length != source_molecule.sequence.len()
            || to.length != destination_molecule.sequence.len()
            || source_molecule.sequence.alphabet() != destination_molecule.sequence.alphabet()
        {
            return Err(bad(
                "source/destination token does not bind the actual native molecules",
            ));
        }
        if to
            .projections
            .iter()
            .any(|projection| projection.source == source && projection.occurrence == occurrence)
        {
            return Err(bad(
                "source occurrence must be captured once with all of its actual spans",
            ));
        }
        let source_length = from.length;
        let source_feature_count = from.feature_count;
        let destination_length = to.length;
        self.builder.budget.charge(Resource::Projections, 1)?;
        self.builder
            .budget
            .charge(Resource::Segments, segments.len())?;
        let feature_count = projection
            .source_feature_indices
            .len()
            .checked_add(projection.unmapped_source_feature_indices.len())
            .and_then(|count| count.checked_add(projection.dropped_source_feature_indices.len()))
            .and_then(|count| count.checked_add(changes.len()))
            .ok_or(Error::TooLarge)?;
        self.builder
            .budget
            .charge(Resource::FeatureOrdinals, feature_count)?;
        let mut change_positions = BTreeMap::new();
        for change in changes {
            if kind != ProjectionKind::EditedTemplate
                || !change.before.is_ascii()
                || !change.after.is_ascii()
                || change_positions
                    .insert(change.destination_position, change)
                    .is_some()
            {
                return Err(bad("distinct native edited-template changes required"));
            }
        }
        let mut matched_changes = BTreeSet::new();
        for segment in segments {
            let end = segment_end(segment)?;
            if segment.source_start >= segment.source_end
                || segment.source_end > source_length
                || end > destination_length
            {
                return Err(bad("projection span exceeds actual native references"));
            }
            self.builder.budget.charge(
                Resource::ScannedBytes,
                segment.source_end - segment.source_start,
            )?;
            for destination_position in segment.destination_start..end {
                let offset = destination_position - segment.destination_start;
                let source_position = if segment.reverse {
                    segment.source_end - 1 - offset
                } else {
                    segment.source_start + offset
                };
                let original = source_molecule.sequence.as_str().as_bytes()[source_position];
                let expected = if segment.reverse {
                    collagene_model::complement(original, source_molecule.sequence.alphabet())
                } else {
                    original
                };
                let actual =
                    destination_molecule.sequence.as_str().as_bytes()[destination_position];
                if let Some(change) = change_positions.get(&destination_position) {
                    matched_changes.insert(destination_position);
                    if change.source_position != source_position
                        || change.before != original as char
                        || change.after != actual as char
                        || expected == actual
                    {
                        return Err(bad(
                            "edited-template change differs from actual mapped calls",
                        ));
                    }
                } else if actual != expected {
                    return Err(bad(
                        "copied correspondence differs from actual native calls",
                    ));
                }
            }
        }
        if matched_changes.len() != changes.len() {
            return Err(bad(
                "edited-template change lies outside its actual projection",
            ));
        }
        let all_indices = projection
            .source_feature_indices
            .iter()
            .chain(&projection.unmapped_source_feature_indices)
            .chain(&projection.dropped_source_feature_indices);
        if all_indices
            .clone()
            .any(|index| *index >= source_feature_count)
        {
            return Err(bad(
                "source feature ordinal is outside the sealed native snapshot",
            ));
        }
        let mapped = projection
            .source_feature_indices
            .iter()
            .copied()
            .collect::<BTreeSet<_>>();
        let mut absent = BTreeSet::new();
        for index in projection
            .unmapped_source_feature_indices
            .iter()
            .chain(&projection.dropped_source_feature_indices)
        {
            if mapped.contains(index) || !absent.insert(*index) {
                return Err(bad(
                    "source feature disposition is duplicated or contradictory",
                ));
            }
        }
        if mapped.len().checked_add(absent.len()) != Some(source_feature_count) {
            return Err(bad(
                "canonical projection did not retain every original feature disposition",
            ));
        }
        let mut features = Vec::with_capacity(projection.features.len());
        for (offset, (feature, source_index)) in projection
            .features
            .iter()
            .zip(&projection.source_feature_indices)
            .enumerate()
        {
            let destination_feature_index = destination_feature_indices[offset];
            features.push(FeatureProjection {
                source_feature_index: *source_index,
                destination_feature_index,
                source_annotation_revision: feature_revision(
                    &source_molecule.features[*source_index],
                    &self.builder.budget,
                )?,
                projected_annotation_revision: feature_revision(feature, &self.builder.budget)?,
            });
        }
        self.builder.trace.snapshots[destination.index]
            .projections
            .push(ProjectionReceipt {
                source,
                occurrence,
                kind,
                segments: segments.to_vec(),
                features,
                unmapped_source_feature_indices: projection.unmapped_source_feature_indices.clone(),
                dropped_source_feature_indices: projection.dropped_source_feature_indices.clone(),
                changes: changes.to_vec(),
            });
        Ok(())
    }
    fn topo_scope(
        &self,
        role: &str,
        scopes: &[crate::verification::projection::TopoScope],
    ) -> Result<Vec<OriginOccurrence>, Error> {
        if role.is_empty()
            || role.len() > 4096
            || role.chars().any(char::is_control)
            || scopes.is_empty()
            || scopes.len() > 8
        {
            return Err(bad(
                "actual bounded retained graph root role/scope required",
            ));
        }
        self.builder
            .budget
            .charge(Resource::NamespaceLinks, scopes.len())?;
        self.builder
            .budget
            .charge(Resource::ScannedBytes, role.len())?;
        for scope in scopes {
            if scope.occurrence > 128 || !matches!(&*scope.binding, OriginBinding::TopoGraph { .. })
            {
                return Err(bad(
                    "native graph scope must retain its actual graph binding/node occurrence",
                ));
            }
            self.builder
                .budget
                .charge(Resource::ScannedBytes, scope.binding.admitted_bytes()?)?;
        }
        let scope = scopes
            .iter()
            .map(|scope| OriginOccurrence {
                step_index: self.step,
                occurrence: scope.occurrence,
                binding: scope.binding.clone(),
            })
            .collect::<Vec<_>>();
        Ok(scope)
    }
    fn register_topo_source(
        &mut self,
        role: &str,
        molecule: &Molecule,
        scopes: &[crate::verification::projection::TopoScope],
    ) -> Result<SnapshotToken, Error> {
        let scope = self.topo_scope(role, scopes)?;
        let token = self.builder.allocate(molecule, None)?;
        let revision = molecule.revision();
        let snapshot = &mut self.builder.trace.snapshots[token.index];
        snapshot.stage = Some(self.step);
        snapshot.scope = scope.clone();
        snapshot.revision = Some(revision.clone());
        snapshot.feature_count = molecule.features.len();
        snapshot.identity = Some(SourceIdentity {
            scope,
            source: SourceHandle::TopoRole(role.into()),
            revision,
        });
        Ok(token)
    }
    fn bind_topo_source(
        &mut self,
        token: SnapshotToken,
        role: &str,
        molecule: &Molecule,
        scopes: &[crate::verification::projection::TopoScope],
    ) -> Result<(), Error> {
        let scope = self.topo_scope(role, scopes)?;
        self.builder.charge_model(molecule)?;
        let snapshot = self.builder.snapshot(token)?;
        if snapshot.stage != Some(self.step)
            || snapshot.identity.is_some()
            || snapshot.revision.as_deref() != Some(molecule.revision().as_str())
        {
            return Err(bad(
                "retained graph role must label its exact sealed native presentation copy",
            ));
        }
        let snapshot = &mut self.builder.trace.snapshots[token.index];
        snapshot.scope = scope.clone();
        snapshot.identity = Some(SourceIdentity {
            scope,
            source: SourceHandle::TopoRole(role.into()),
            revision: snapshot.revision.clone().expect("sealed source"),
        });
        Ok(())
    }
    fn bind_reference_event(
        &mut self,
        token: SnapshotToken,
        event: &crate::material::lineage::EventBinding,
        molecule: &Molecule,
    ) -> Result<(), Error> {
        self.builder
            .budget
            .charge(Resource::ScannedBytes, reference_event_bytes(event)?)?;
        self.builder.charge_model(molecule)?;
        let snapshot = self.builder.snapshot(token)?;
        if snapshot.stage != Some(self.step)
            || snapshot.identity.is_some()
            || snapshot.revision.as_deref() != Some(event.revision.as_str())
            || molecule.revision() != event.revision
        {
            return Err(bad(
                "native event label must name its exact sealed retrieved reference",
            ));
        }
        let snapshot = &mut self.builder.trace.snapshots[token.index];
        snapshot.identity = Some(SourceIdentity {
            scope: snapshot.scope.clone(),
            source: SourceHandle::ReferenceEvent(event.clone()),
            revision: event.revision.clone(),
        });
        Ok(())
    }
    fn import(
        &mut self,
        binding: OriginBinding,
        mut trace: ReconstructionTrace,
        selected_step: usize,
        occurrence: u32,
    ) -> Result<SnapshotToken, Error> {
        if occurrence >= 4096 || !Arc::ptr_eq(&self.builder.budget.ledger, &trace.budget.ledger) {
            return Err(bad(
                "nested trace must come from this same enclosing native job budget",
            ));
        }
        let bound = match (&binding, &trace.context) {
            (
                OriginBinding::Run {
                    selection,
                    plan_revision,
                    ..
                },
                TraceContext::Plan { revision },
            ) => selection.step_index == selected_step && plan_revision == revision,
            (
                OriginBinding::Material {
                    context_revision, ..
                },
                TraceContext::Material {
                    context_revision: actual,
                },
            ) => selected_step == 0 && context_revision == actual,
            _ => false,
        };
        if !bound {
            return Err(bad(
                "retained origin binding does not name this actual native trace context/selected product",
            ));
        }
        let selected = trace
            .step(selected_step)
            .ok_or_else(|| bad("nested selected step has no actual committed trace"))?
            .final_snapshot;
        let links = trace.snapshots.iter().try_fold(0usize, |total, snapshot| {
            if snapshot.scope.len() >= 8 {
                return Err(Error::TooLarge);
            }
            total
                .checked_add(snapshot.scope.len() + 1)
                .ok_or(Error::TooLarge)
        })?;
        self.builder
            .budget
            .charge(Resource::NamespaceLinks, links)?;
        self.builder
            .budget
            .charge(Resource::ScannedBytes, binding.admitted_bytes()?)?;
        let binding = Arc::new(binding);
        let prefix = OriginOccurrence {
            step_index: self.step,
            occurrence,
            binding,
        };
        let mut translated = Vec::with_capacity(trace.snapshots.len());
        let base = self.builder.trace.snapshots.len();
        for (offset, _) in trace.snapshots.iter().enumerate() {
            self.builder.next_serial = self
                .builder
                .next_serial
                .checked_add(1)
                .ok_or(Error::TooLarge)?;
            translated.push(SnapshotToken {
                trace_id: self.builder.trace_id,
                index: base.checked_add(offset).ok_or(Error::TooLarge)?,
                serial: self.builder.next_serial,
            });
        }
        let rewrite = |token: SnapshotToken| -> Result<SnapshotToken, Error> {
            if trace
                .snapshots
                .get(token.index)
                .is_none_or(|snapshot| snapshot.token != token)
            {
                return Err(bad("nested trace contains a foreign or stale native token"));
            }
            translated
                .get(token.index)
                .copied()
                .ok_or_else(|| bad("nested native token translation is incomplete"))
        };
        let selected = rewrite(selected)?;
        // Validate every edge before moving any snapshot into the enclosing job.
        for snapshot in &trace.snapshots {
            for projection in &snapshot.projections {
                rewrite(projection.source)?;
            }
        }
        for (offset, mut snapshot) in trace.snapshots.drain(..).enumerate() {
            snapshot.token = translated[offset];
            snapshot.stage = Some(self.step);
            snapshot.scope.insert(0, prefix.clone());
            if let Some(identity) = &mut snapshot.identity {
                identity.scope.insert(0, prefix.clone());
            }
            for projection in &mut snapshot.projections {
                projection.source = translated[projection.source.index];
            }
            self.builder.trace.snapshots.push(snapshot);
        }
        Ok(selected)
    }
    /// Exact native dedup retain-loop remapping, including earlier receipts.
    pub fn remap_features(
        &mut self,
        destination: SnapshotToken,
        indices: &[usize],
    ) -> Result<(), Error> {
        let snapshot = self.builder.snapshot(destination)?;
        if snapshot.step != Some(self.step) || snapshot.revision.is_some() || indices.len() > 4096 {
            return Err(bad("feature remapping requires the actual current draft"));
        }
        for projection in &snapshot.projections {
            if projection
                .features
                .iter()
                .any(|feature| feature.destination_feature_index >= indices.len())
            {
                return Err(bad("native dedup map omits a previously captured feature"));
            }
        }
        let snapshot = &mut self.builder.trace.snapshots[destination.index];
        for projection in &mut snapshot.projections {
            for feature in &mut projection.features {
                feature.destination_feature_index = indices[feature.destination_feature_index];
            }
        }
        Ok(())
    }
    pub fn introduce(
        &mut self,
        destination: SnapshotToken,
        start: usize,
        end: usize,
        kind: IntroductionKind,
    ) -> Result<(), Error> {
        let snapshot = self.builder.snapshot(destination)?;
        if snapshot.step != Some(self.step)
            || snapshot.revision.is_some()
            || start >= end
            || end > snapshot.length
        {
            return Err(bad("bounded actual draft introduction range required"));
        }
        self.builder.budget.charge(Resource::Segments, 1)?;
        self.builder.trace.snapshots[destination.index]
            .introductions
            .push(Introduction { start, end, kind });
        Ok(())
    }
    pub fn seal(&mut self, destination: SnapshotToken, molecule: &Molecule) -> Result<(), Error> {
        self.builder.snapshot(destination)?;
        self.builder.charge_model(molecule)?;
        let snapshot = self.builder.snapshot(destination)?;
        if snapshot.step != Some(self.step)
            || snapshot.revision.is_some()
            || dna_revision(molecule)? != snapshot.dna_revision
        {
            return Err(bad(
                "seal must bind the actual unchanged-reference native draft once",
            ));
        }
        let snapshot = &self.builder.trace.snapshots[destination.index];
        for projection in &snapshot.projections {
            for feature in &projection.features {
                let actual = molecule
                    .features
                    .get(feature.destination_feature_index)
                    .ok_or_else(|| {
                        bad("projected feature was not appended at its actual ordinal")
                    })?;
                if feature_revision(actual, &self.builder.budget)?
                    != feature.projected_annotation_revision
                {
                    return Err(bad(
                        "actual appended feature differs from its captured projection",
                    ));
                }
            }
        }
        let snapshot = &mut self.builder.trace.snapshots[destination.index];
        snapshot.revision = Some(molecule.revision());
        snapshot.feature_count = molecule.features.len();
        Ok(())
    }
    pub fn commit(mut self, final_snapshot: SnapshotToken) -> Result<SnapshotToken, Error> {
        let snapshot = self.builder.snapshot(final_snapshot)?;
        if snapshot.step != Some(self.step) {
            return Err(bad("final native snapshot belongs to another operation"));
        }
        let revision = snapshot
            .revision
            .clone()
            .ok_or_else(|| bad("final native snapshot is not sealed"))?;
        let product_length = snapshot.length;
        let local_coverage_complete = coverage(snapshot)?;
        self.builder.trace.steps[self.step] = Some(StepTrace {
            step_index: self.step,
            product_revision: revision.clone(),
            product_length,
            local_coverage_complete,
            final_snapshot,
        });
        self.builder.trace.snapshots[final_snapshot.index].identity = Some(SourceIdentity {
            scope: vec![],
            source: match self.builder.kind {
                TraceKind::Plan => SourceHandle::Plan(Source::Step(self.step)),
                TraceKind::Material => SourceHandle::MaterialProduct,
                TraceKind::ReferenceEdit => SourceHandle::ReferenceProduct {
                    action_index: self.step,
                },
            },
            revision,
        });
        self.committed = true;
        Ok(final_snapshot)
    }
}
impl Drop for StepBuilder<'_> {
    fn drop(&mut self) {
        if !self.committed {
            self.builder.trace.snapshots.truncate(self.snapshot_start);
        }
        // Ledger charges and serials deliberately survive failed/discarded work.
    }
}
fn coverage(snapshot: &Snapshot) -> Result<bool, Error> {
    let mut spans = snapshot
        .introductions
        .iter()
        .map(|span| (span.start, span.end))
        .collect::<Vec<_>>();
    for projection in &snapshot.projections {
        for segment in &projection.segments {
            spans.push((segment.destination_start, segment_end(segment)?));
        }
    }
    spans.sort_unstable();
    let mut end = 0usize;
    for (start, stop) in spans {
        if start > end {
            return Ok(false);
        }
        end = end.max(stop);
    }
    Ok(end == snapshot.length)
}

#[derive(Clone, Debug)]
pub(crate) enum OccurrenceMapping {
    Base {
        destination_position: usize,
        source_position: usize,
        reverse: bool,
    },
    Feature {
        destination_feature_index: usize,
        source_feature_index: usize,
    },
}
#[derive(Clone, Debug)]
pub(crate) struct OccurrenceKey {
    pub source_snapshot: SnapshotToken,
    pub mapping: OccurrenceMapping,
    pub scope: Vec<OriginOccurrence>,
    pub step_index: usize,
    pub destination_snapshot: SnapshotToken,
    pub occurrence: u32,
    pub kind: ProjectionKind,
}
#[derive(Clone, Debug)]
pub(crate) struct ChangeAtCall<'a> {
    pub change: &'a MappedChange,
    pub occurrence: OccurrenceKey,
    /// Orientation from this intermediate's forward reference to requested output.
    pub output_reverse: bool,
}
#[derive(Clone, Debug)]
pub(crate) struct OriginCall<'a> {
    pub source: &'a SourceIdentity,
    pub position: usize,
    pub reverse: bool,
    pub path: Vec<OccurrenceKey>,
    pub changes: Vec<ChangeAtCall<'a>>,
}
#[derive(Debug)]
pub(crate) struct IntroducedCall<'a> {
    pub snapshot: SnapshotToken,
    pub position: usize,
    pub reverse: bool,
    pub kind: &'a IntroductionKind,
    pub path: Vec<OccurrenceKey>,
}
#[derive(Debug)]
pub(crate) struct LookupEvidence<'a> {
    pub origins: Vec<OriginCall<'a>>,
    pub introductions: Vec<IntroducedCall<'a>>,
    pub complete: bool,
    pub reason: Option<&'static str>,
}
/// An exact baseline identity selected from the actual native receipt graph.
#[derive(Clone, Copy, Debug)]
pub(crate) enum LookupStop<'a> {
    PlanSource(&'a Source),
    Snapshot(SnapshotToken),
    Targets(&'a LookupTargets),
}
#[derive(Debug)]
pub(crate) struct LookupTargets {
    start: SnapshotToken,
    mask: Vec<u64>,
}
/// One session carries actual receipt-walking work across a complete QC job;
/// callers must not create a new session for every codon/base comparison.
pub(crate) struct LookupSession<'a> {
    trace: &'a ReconstructionTrace,
    remaining: usize,
}
impl ProjectionReceipt {
    pub fn source_snapshot(&self) -> SnapshotToken {
        self.source
    }
}
impl StepTrace {
    pub fn final_snapshot(&self) -> SnapshotToken {
        self.final_snapshot
    }
}
impl ReconstructionTrace {
    pub fn lookup_targets(
        &self,
        start: SnapshotToken,
        tokens: &[SnapshotToken],
    ) -> Result<LookupTargets, Error> {
        checked_snapshot(self, start)?;
        if tokens.len() > 4096 {
            return Err(Error::TooLarge);
        }
        let mut mask = vec![0u64; self.snapshots.len().div_ceil(64)];
        for &token in tokens {
            validate_snapshot_stop(self, start, token)?;
            mask[token.index / 64] |= 1u64 << (token.index % 64);
        }
        Ok(LookupTargets { start, mask })
    }
    pub fn snapshot_scope(&self, token: SnapshotToken) -> Option<&[OriginOccurrence]> {
        self.snapshots
            .get(token.index)
            .filter(|snapshot| snapshot.token == token)
            .map(|snapshot| snapshot.scope.as_slice())
    }
    pub fn snapshot_step(&self, token: SnapshotToken) -> Option<usize> {
        self.snapshots
            .get(token.index)
            .filter(|snapshot| snapshot.token == token)
            .and_then(|snapshot| snapshot.step)
    }
    pub fn sources(&self) -> impl Iterator<Item = (SnapshotToken, &SourceIdentity)> {
        self.snapshots.iter().filter_map(|snapshot| {
            snapshot
                .identity
                .as_ref()
                .map(|source| (snapshot.token, source))
        })
    }
    pub fn context(&self) -> &TraceContext {
        &self.context
    }
    pub fn step(&self, step_index: usize) -> Option<&StepTrace> {
        self.steps.get(step_index).and_then(Option::as_ref)
    }
    pub fn projections(&self, step_index: usize) -> Option<&[ProjectionReceipt]> {
        let step = self.step(step_index)?;
        Some(&self.snapshots.get(step.final_snapshot.index)?.projections)
    }
    pub fn snapshot_source(&self, token: SnapshotToken) -> Option<&SourceIdentity> {
        self.snapshots
            .get(token.index)
            .filter(|snapshot| snapshot.token == token)?
            .identity
            .as_ref()
    }
    pub fn snapshot_projections(&self, token: SnapshotToken) -> Option<&[ProjectionReceipt]> {
        Some(
            &self
                .snapshots
                .get(token.index)
                .filter(|snapshot| snapshot.token == token)?
                .projections,
        )
    }
    pub fn lookup_session(&self) -> LookupSession<'_> {
        LookupSession {
            trace: self,
            remaining: self.limits.max_lookup_work,
        }
    }
}
#[derive(Clone)]
struct PendingCall<'a> {
    snapshot: SnapshotToken,
    position: usize,
    reverse: bool,
    path: Vec<OccurrenceKey>,
    changes: Vec<ChangeAtCall<'a>>,
}
fn checked_snapshot(trace: &ReconstructionTrace, token: SnapshotToken) -> Result<&Snapshot, Error> {
    trace
        .snapshots
        .get(token.index)
        .filter(|snapshot| snapshot.token == token && snapshot.revision.is_some())
        .ok_or_else(|| bad("native snapshot is foreign, stale, missing or unsealed"))
}
fn validate_snapshot_stop(
    trace: &ReconstructionTrace,
    start: SnapshotToken,
    token: SnapshotToken,
) -> Result<(), Error> {
    let root = checked_snapshot(trace, start)?;
    let snapshot = checked_snapshot(trace, token)?;
    let Some(identity) = &snapshot.identity else {
        return Err(bad("comparison snapshot needs its actual source identity"));
    };
    if start == token {
        return Err(bad(
            "comparison source cannot be the requested product itself",
        ));
    }
    if root.scope.is_empty() && snapshot.scope.is_empty() {
        let source_index = match &identity.source {
            SourceHandle::Plan(Source::Step(index)) => Some(*index),
            SourceHandle::ReferenceProduct { action_index } => Some(*action_index),
            _ => None,
        };
        if source_index
            .zip(root.step)
            .is_some_and(|(source, product)| source >= product)
        {
            return Err(bad(
                "comparison source step must precede the requested actual product",
            ));
        }
    }
    Ok(())
}
fn validate_stop(
    trace: &ReconstructionTrace,
    start: SnapshotToken,
    stop: Option<LookupStop<'_>>,
) -> Result<(), Error> {
    match stop {
        Some(LookupStop::Snapshot(token)) => validate_snapshot_stop(trace, start, token),
        Some(LookupStop::Targets(targets)) => {
            if targets.start != start {
                return Err(bad(
                    "lookup targets belong to a different actual starting snapshot",
                ));
            }
            Ok(())
        }
        Some(LookupStop::PlanSource(Source::Step(index))) => {
            let root = checked_snapshot(trace, start)?;
            if root.scope.is_empty() && root.step.is_some_and(|step| *index >= step) {
                return Err(bad(
                    "comparison source step must precede the requested product",
                ));
            }
            Ok(())
        }
        _ => Ok(()),
    }
}
fn matches_stop(snapshot: &Snapshot, stop: Option<LookupStop<'_>>) -> bool {
    match stop {
        Some(LookupStop::Snapshot(token)) => snapshot.token == token,
        Some(LookupStop::Targets(targets)) => targets.mask.get(snapshot.token.index / 64)
            .is_some_and(|word| word & (1u64 << (snapshot.token.index % 64)) != 0),
        Some(LookupStop::PlanSource(expected)) => snapshot.scope.is_empty()
            && snapshot.identity.as_ref().is_some_and(|identity| matches!(&identity.source, SourceHandle::Plan(source) if source == expected)),
        None => false,
    }
}
fn path_cost(path: &[OccurrenceKey]) -> usize {
    path.iter().map(|edge| 1 + edge.scope.len()).sum()
}
fn changes_cost(changes: &[ChangeAtCall<'_>]) -> usize {
    changes
        .iter()
        .map(|change| 1 + change.occurrence.scope.len())
        .sum()
}
impl<'a> LookupSession<'a> {
    pub fn remaining(&self) -> usize {
        self.remaining
    }
    pub fn restrict_remaining(&mut self, maximum: usize) {
        self.remaining = self.remaining.min(maximum);
    }
    fn charge(&mut self, amount: usize) -> bool {
        if amount > self.remaining {
            self.remaining = 0;
            false
        } else {
            self.remaining -= amount;
            true
        }
    }
    /// Composes only actual captured segments. All original input roles remain
    /// distinct. stop_at can select an earlier actual step as a comparison source;
    /// unrelated input origins are retained too, never silently selected away.
    pub fn at(
        &mut self,
        step_index: usize,
        position: usize,
        stop_at: Option<LookupStop<'_>>,
    ) -> Result<LookupEvidence<'a>, Error> {
        let start = self
            .trace
            .step(step_index)
            .ok_or_else(|| bad("selected native step has no committed trace"))?
            .final_snapshot;
        self.at_snapshot(start, position, stop_at)
    }
    pub fn at_snapshot(
        &mut self,
        start: SnapshotToken,
        position: usize,
        stop_at: Option<LookupStop<'_>>,
    ) -> Result<LookupEvidence<'a>, Error> {
        let trace = self.trace;
        let root = checked_snapshot(trace, start)?;
        validate_stop(trace, start, stop_at)?;
        if position >= root.length {
            return Err(bad("lookup call is outside the actual selected snapshot"));
        }
        let mut result = LookupEvidence {
            origins: vec![],
            introductions: vec![],
            complete: true,
            reason: None,
        };
        let mut pending = vec![PendingCall {
            snapshot: start,
            position,
            reverse: false,
            path: vec![],
            changes: vec![],
        }];
        while let Some(call) = pending.pop() {
            if !self.charge(1) {
                result.complete = false;
                result.reason = Some("Shared receipt lookup work was exhausted.");
                break;
            }
            let snapshot = trace
                .snapshots
                .get(call.snapshot.index)
                .filter(|snapshot| snapshot.token == call.snapshot)
                .ok_or_else(|| bad("receipt references a missing/stale native snapshot"))?;
            if snapshot.revision.is_none() {
                return Err(bad("unsealed snapshot cannot supply causal correspondence"));
            }
            if let Some(source) = &snapshot.identity {
                if matches!(
                    source.source,
                    SourceHandle::Plan(Source::Input(_)) | SourceHandle::Material(_)
                ) || matches!(
                    source.source,
                    SourceHandle::TopoRole(_) | SourceHandle::ReferenceEvent(_)
                ) && snapshot.projections.is_empty()
                    || matches_stop(snapshot, stop_at)
                {
                    if result.origins.len() + result.introductions.len() >= 64 {
                        result.complete = false;
                        result.reason =
                            Some("Native origin alternatives exceeded the 64-result lookup bound.");
                        break;
                    }
                    result.origins.push(OriginCall {
                        source,
                        position: call.position,
                        reverse: call.reverse,
                        path: call.path,
                        changes: call.changes,
                    });
                    continue;
                }
            }
            if call.path.len() >= 64 {
                result.complete = false;
                result.reason = Some("Native correspondence path exceeds the 64-edge bound.");
                continue;
            }
            let mut matched = false;
            for introduction in &snapshot.introductions {
                if !self.charge(1) {
                    result.complete = false;
                    result.reason = Some("Shared receipt lookup work was exhausted.");
                    break;
                }
                if introduction.start <= call.position && call.position < introduction.end {
                    matched = true;
                    if result.origins.len() + result.introductions.len() >= 64 {
                        result.complete = false;
                        result.reason =
                            Some("Native origin alternatives exceeded the 64-result lookup bound.");
                        break;
                    }
                    if !self.charge(path_cost(&call.path)) {
                        result.complete = false;
                        result.reason = Some("Shared receipt lookup work was exhausted.");
                        break;
                    }
                    result.introductions.push(IntroducedCall {
                        snapshot: snapshot.token,
                        position: call.position,
                        reverse: call.reverse,
                        kind: &introduction.kind,
                        path: call.path.clone(),
                    });
                }
            }
            for projection in &snapshot.projections {
                for segment in &projection.segments {
                    if !self.charge(1) {
                        result.complete = false;
                        result.reason = Some("Shared receipt lookup work was exhausted.");
                        break;
                    }
                    let end = segment_end(segment)?;
                    if call.position < segment.destination_start || call.position >= end {
                        continue;
                    }
                    matched = true;
                    if pending.len() >= 64
                        || !self.charge(
                            path_cost(&call.path)
                                + changes_cost(&call.changes)
                                + snapshot.scope.len()
                                + 1,
                        )
                    {
                        result.complete = false;
                        result.reason =
                            Some("Native origin breadth or shared lookup work bound was reached.");
                        break;
                    }
                    let offset = call.position - segment.destination_start;
                    let position = if segment.reverse {
                        segment.source_end - 1 - offset
                    } else {
                        segment.source_start + offset
                    };
                    let occurrence = OccurrenceKey {
                        source_snapshot: projection.source,
                        mapping: OccurrenceMapping::Base {
                            destination_position: call.position,
                            source_position: position,
                            reverse: segment.reverse,
                        },
                        scope: snapshot.scope.clone(),
                        step_index: snapshot.step.ok_or_else(|| {
                            bad("projected native snapshot lacks its producing step")
                        })?,
                        destination_snapshot: snapshot.token,
                        occurrence: projection.occurrence,
                        kind: projection.kind.clone(),
                    };
                    let mut path = call.path.clone();
                    path.push(occurrence.clone());
                    let mut changes = call.changes.clone();
                    for change in &projection.changes {
                        if !self.charge(1) {
                            result.complete = false;
                            result.reason = Some("Shared receipt lookup work was exhausted.");
                            break;
                        }
                        if change.destination_position == call.position {
                            changes.push(ChangeAtCall {
                                change,
                                occurrence: occurrence.clone(),
                                output_reverse: call.reverse,
                            });
                        }
                    }
                    pending.push(PendingCall {
                        snapshot: projection.source,
                        position,
                        reverse: call.reverse ^ segment.reverse,
                        path,
                        changes,
                    });
                }
                if self.remaining == 0 {
                    break;
                }
            }
            if !matched {
                result.complete = false;
                result.reason = Some(
                    "This actual product call has no captured source/introduction correspondence.",
                );
            }
            if self.remaining == 0 {
                result.complete = false;
                result.reason = Some("Shared receipt lookup work was exhausted.");
                break;
            }
        }
        Ok(result)
    }
}

#[derive(Debug)]
pub(crate) struct FeatureOrigin<'a> {
    pub source: &'a SourceIdentity,
    pub feature_index: usize,
    pub annotation_revision: &'a str,
    pub path: Vec<OccurrenceKey>,
}
#[derive(Debug)]
pub(crate) struct FeatureLookupEvidence<'a> {
    pub origins: Vec<FeatureOrigin<'a>>,
    pub complete: bool,
    pub reason: Option<&'static str>,
}
impl<'a> LookupSession<'a> {
    pub fn feature_at(
        &mut self,
        step_index: usize,
        feature_index: usize,
        stop_at: Option<LookupStop<'_>>,
    ) -> Result<FeatureLookupEvidence<'a>, Error> {
        let start = self
            .trace
            .step(step_index)
            .ok_or_else(|| bad("selected native step has no committed trace"))?
            .final_snapshot;
        self.feature_at_snapshot(start, feature_index, stop_at)
    }
    pub fn feature_at_snapshot(
        &mut self,
        start: SnapshotToken,
        feature_index: usize,
        stop_at: Option<LookupStop<'_>>,
    ) -> Result<FeatureLookupEvidence<'a>, Error> {
        let trace = self.trace;
        let root = checked_snapshot(trace, start)?;
        validate_stop(trace, start, stop_at)?;
        if feature_index >= root.feature_count {
            return Err(bad("feature ordinal is outside actual selected snapshot"));
        }
        let mut result = FeatureLookupEvidence {
            origins: vec![],
            complete: true,
            reason: None,
        };
        let mut pending = vec![(
            start,
            feature_index,
            None::<&'a str>,
            Vec::<OccurrenceKey>::new(),
        )];
        while let Some((token, index, expected, path)) = pending.pop() {
            if !self.charge(1) {
                result.complete = false;
                result.reason = Some("Shared receipt lookup work was exhausted.");
                break;
            }
            let snapshot = trace
                .snapshots
                .get(token.index)
                .filter(|snapshot| snapshot.token == token)
                .ok_or_else(|| bad("feature receipt references a missing/stale snapshot"))?;
            if snapshot.revision.is_none() {
                return Err(bad("unsealed feature receipt cannot supply correspondence"));
            }
            if let Some(source) = &snapshot.identity {
                if matches!(
                    source.source,
                    SourceHandle::Plan(Source::Input(_)) | SourceHandle::Material(_)
                ) || matches!(
                    source.source,
                    SourceHandle::TopoRole(_) | SourceHandle::ReferenceEvent(_)
                ) && snapshot.projections.is_empty()
                    || matches_stop(snapshot, stop_at)
                {
                    let Some(annotation_revision) = expected else {
                        return Err(bad(
                            "source feature correspondence lacks its actual captured annotation pin",
                        ));
                    };
                    if result.origins.len() >= 64 {
                        result.complete = false;
                        result.reason =
                            Some("Feature origin alternatives exceed the 64-result bound.");
                        break;
                    }
                    result.origins.push(FeatureOrigin {
                        source,
                        feature_index: index,
                        annotation_revision,
                        path,
                    });
                    continue;
                }
            }
            if path.len() >= 64 {
                result.complete = false;
                result.reason = Some("Feature correspondence path exceeds the 64-edge bound.");
                continue;
            }
            let mut matched = false;
            for projection in &snapshot.projections {
                for feature in &projection.features {
                    if !self.charge(1) {
                        result.complete = false;
                        result.reason = Some("Shared receipt lookup work was exhausted.");
                        break;
                    }
                    if feature.destination_feature_index != index {
                        continue;
                    }
                    matched = true;
                    if expected.is_some_and(|pin| pin != feature.projected_annotation_revision) {
                        result.complete = false;
                        result.reason = Some(
                            "The actual annotation changed beyond this captured feature correspondence.",
                        );
                        continue;
                    }
                    if pending.len() >= 64
                        || !self.charge(path_cost(&path) + snapshot.scope.len() + 1)
                    {
                        result.complete = false;
                        result.reason =
                            Some("Feature origin breadth or lookup work bound was reached.");
                        break;
                    }
                    let mut next = path.clone();
                    next.push(OccurrenceKey {
                        source_snapshot: projection.source,
                        mapping: OccurrenceMapping::Feature {
                            destination_feature_index: index,
                            source_feature_index: feature.source_feature_index,
                        },
                        scope: snapshot.scope.clone(),
                        step_index: snapshot.step.ok_or_else(|| {
                            bad("feature projection lacks its actual producing step")
                        })?,
                        destination_snapshot: snapshot.token,
                        occurrence: projection.occurrence,
                        kind: projection.kind.clone(),
                    });
                    pending.push((
                        projection.source,
                        feature.source_feature_index,
                        Some(feature.source_annotation_revision.as_str()),
                        next,
                    ));
                }
                if self.remaining == 0 {
                    break;
                }
            }
            if !matched {
                result.complete = false;
                result.reason = Some(
                    "This destination feature has no captured ordinal correspondence; labels/equal content do not identify an ancestor.",
                );
            }
            if self.remaining == 0 {
                result.complete = false;
                result.reason = Some("Shared receipt lookup work was exhausted.");
                break;
            }
        }
        Ok(result)
    }
}
