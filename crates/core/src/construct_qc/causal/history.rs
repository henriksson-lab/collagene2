//! Cached actual native event replay and exact event-bound coordinate traversal.
//! Labels, equal sequences and annotation names never select a predecessor.
use super::{Request, TransitionPin, Work};
use crate::{
    chemistry::tracked_duplex::CleavageBudget,
    documents::{
        history::transition::{AnalysisProducer, AnalysisRecord, Transition},
        reference_edit::Operation,
        workspace::Workspace,
    },
    material::{
        self,
        lineage::{self, EventBinding, InputEventBinding, ProductBinding},
    },
    planning::plan::Source,
    verification::{
        self,
        ancestry::{
            LookupSession, LookupStop, LookupTargets, OccurrenceKey, OriginBinding,
            ReconstructionTrace, SnapshotToken, SourceHandle, SourceIdentity, TraceBudget,
        },
        replay_budget::ReplayBudget,
    },
};
use collagene_model::{Error, annotations::biological_path::BiologicalCall};
use std::collections::{BTreeMap, BTreeSet, HashMap};

pub(super) struct CallMapping {
    pub parent_calls: Vec<BiologicalCall>,
    pub introduced: bool,
    pub foreign: bool,
    pub complete: bool,
    pub changed: bool,
    pub causes: Vec<TransitionPin>,
    pub reason: Option<String>,
}
pub(super) struct FeatureMapping {
    pub parent_feature_indices: Vec<usize>,
    pub complete: bool,
    pub reason: Option<String>,
}
pub(super) struct HistoryEvidence {
    pub available: bool,
    pub reason: String,
}
#[derive(Clone)]
struct RunBindings {
    inputs: Vec<InputEventBinding>,
    products: Vec<ProductBinding>,
}
enum Root {
    Run(RunBindings),
    Material(OriginBinding),
    Reference { event: EventBinding, last: usize },
}
struct NativeTrace {
    trace: ReconstructionTrace,
    bindings: HashMap<SnapshotToken, EventBinding>,
}
#[derive(Clone)]
enum Node {
    Root,
    Fork {
        source: EventBinding,
        pin: TransitionPin,
    },
    Native {
        trace: usize,
        start: SnapshotToken,
        pin: TransitionPin,
    },
    Unknown(String),
}
pub(super) struct History {
    events: BTreeMap<String, EventBinding>,
    nodes: BTreeMap<String, Node>,
    traces: Vec<NativeTrace>,
    history_events: usize,
}
struct Builder<'a, 'budget> {
    replay: ReplayBudget<'budget>,
    workspace: &'a Workspace,
    history: History,
    budget: TraceBudget,
    maximum: usize,
    runs: BTreeMap<String, String>,
    records: BTreeMap<String, lineage::RunRecord>,
    materials: BTreeMap<String, (String, String)>,
    run_traces: BTreeMap<String, Option<usize>>,
    admitted: BTreeSet<String>,
}
fn reason(error: impl std::fmt::Display) -> String {
    error.to_string().chars().take(1024).collect()
}
fn product_binding(product: &ProductBinding) -> EventBinding {
    EventBinding {
        document_id: product.document_id.clone(),
        event_id: product.event_id.clone(),
        revision: product.revision.clone(),
    }
}
impl History {
    pub fn prepare(
        workspace: &Workspace,
        request: &Request,
        trace_budget: TraceBudget,
    ) -> Result<Self, Error> {
        request.validate(workspace)?;
        // The workspace was admitted above. Reuse its canonical record parser
        // once per retained current-kind attachment, never a full-workspace
        // validation for every snapshot/role in the same assessment.
        let records = workspace
            .attachments
            .values()
            .filter(|asset| asset.kind == lineage::WORKSPACE_ATTACHMENT_KIND)
            .map(|asset| {
                let record = lineage::load_run_record(&asset.bytes)?;
                Ok((record.run_id.clone(), record))
            })
            .collect::<Result<BTreeMap<_, _>, Error>>()?;
        let runs = records
            .iter()
            .map(|(id, record)| (id.clone(), record.plan_revision.clone()))
            .collect();
        let materials = workspace
            .attachments
            .iter()
            .filter(|(_, asset)| asset.kind == material::derivation::WORKSPACE_ATTACHMENT_KIND)
            .map(|(hash, asset)| {
                let record = material::derivation::load_record(&asset.bytes)?;
                Ok((
                    record.product.event_id,
                    (record.context_revision, hash.clone()),
                ))
            })
            .collect::<Result<BTreeMap<_, _>, Error>>()?;
        let mut cleavage = CleavageBudget::new(100_000_000)?;
        let mut builder = Builder {
            replay: ReplayBudget::qc_job(&mut cleavage),
            workspace,
            budget: trace_budget,
            maximum: request.limits.max_history_events,
            runs,
            records,
            materials,
            run_traces: BTreeMap::new(),
            admitted: BTreeSet::new(),
            history: Self {
                events: workspace
                    .history
                    .iter()
                    .map(|event| {
                        (
                            event.id.clone(),
                            EventBinding {
                                document_id: event.document_id.clone(),
                                event_id: event.id.clone(),
                                revision: event.revision.clone(),
                            },
                        )
                    })
                    .collect(),
                nodes: BTreeMap::new(),
                traces: vec![],
                history_events: 0,
            },
        };
        for pair in &request.pairs {
            builder.prepare_event(&pair.product);
        }
        Ok(builder.history)
    }
    pub fn history_events(&self) -> usize {
        self.history_events
    }
    pub fn session(&self) -> Session<'_> {
        Session {
            history: self,
            lookups: self
                .traces
                .iter()
                .map(|native| native.trace.lookup_session())
                .collect(),
            targets: HashMap::new(),
        }
    }
    fn valid(&self, event: &EventBinding) -> bool {
        self.events.get(&event.event_id) == Some(event)
    }
}
impl Builder<'_, '_> {
    fn prepare_event(&mut self, binding: &EventBinding) {
        // Native owners resolve their immutable sources directly. There is no
        // postorder dependency, so an explicit bounded agenda handles deep DAGs.
        let mut agenda = vec![binding.clone()];
        let mut queued = BTreeSet::from([binding.event_id.clone()]);
        while let Some(binding) = agenda.pop() {
            queued.remove(&binding.event_id);
            if !self.history.valid(&binding)
                || self.history.nodes.contains_key(&binding.event_id)
                || !self.admit(&binding.event_id)
            {
                continue;
            }
            self.history.nodes.insert(
                binding.event_id.clone(),
                Node::Unknown("Native predecessor is unavailable or cyclic".into()),
            );
            match self.replay(&binding) {
                Ok((node, parents)) => {
                    self.history.nodes.insert(binding.event_id.clone(), node);
                    for parent in parents {
                        if !self.history.nodes.contains_key(&parent.event_id)
                            && agenda.len() < self.maximum
                            && queued.insert(parent.event_id.clone())
                        {
                            agenda.push(parent);
                        }
                    }
                }
                Err(error) => {
                    self.history
                        .nodes
                        .insert(binding.event_id.clone(), Node::Unknown(reason(error)));
                }
            }
        }
    }
    fn admit(&mut self, id: &str) -> bool {
        if self.admitted.contains(id) {
            return true;
        }
        if self.history.history_events >= self.maximum {
            return false;
        }
        self.history.history_events += 1;
        self.admitted.insert(id.into());
        true
    }
    fn run_plan_revision(&self, id: &str) -> Result<String, Error> {
        self.runs
            .get(id)
            .cloned()
            .ok_or_else(|| Error::Format("Exact retained run metadata is unavailable".into()))
    }
    fn resolve_identity(&mut self, identity: &SourceIdentity, root: &Root) -> Option<EventBinding> {
        let scoped_run = identity
            .scope
            .iter()
            .rev()
            .find_map(|scope| match &*scope.binding {
                OriginBinding::Run {
                    inputs, products, ..
                } => Some((inputs.as_slice(), products.as_slice())),
                _ => None,
            });
        let scoped_material = identity
            .scope
            .iter()
            .rev()
            .find_map(|scope| match &*scope.binding {
                binding @ OriginBinding::Material { .. } => Some(binding),
                _ => None,
            });
        let found = match &identity.source {
            SourceHandle::ReferenceEvent(event) => Some(event.clone()),
            SourceHandle::ReferenceProduct { action_index } => match root {
                Root::Reference { event, last } if action_index == last => Some(event.clone()),
                _ => None,
            },
            SourceHandle::Plan(source) => {
                let (inputs, products) = scoped_run.or_else(|| match root {
                    Root::Run(bindings) => {
                        Some((bindings.inputs.as_slice(), bindings.products.as_slice()))
                    }
                    _ => None,
                })?;
                match source {
                    Source::Input(id) => inputs
                        .iter()
                        .find(|input| input.input_id == *id)
                        .map(|input| input.event.clone()),
                    Source::Step(index) => products
                        .get(*index)
                        .filter(|product| product.step_index == *index)
                        .map(product_binding),
                }
            }
            SourceHandle::Material(role) => {
                let binding = scoped_material.or_else(|| match root {
                    Root::Material(binding) => Some(binding),
                    _ => None,
                })?;
                match binding {
                    OriginBinding::Material { sources, .. } => sources
                        .iter()
                        .find(|source| source.role == *role)
                        .map(|source| source.event.clone()),
                    _ => None,
                }
            }
            SourceHandle::MaterialProduct => {
                let binding = scoped_material.or_else(|| match root {
                    Root::Material(binding) => Some(binding),
                    _ => None,
                })?;
                match binding {
                    OriginBinding::Material { product, .. } => Some(product.clone()),
                    _ => None,
                }
            }
            SourceHandle::TopoRole(role) => {
                identity
                    .scope
                    .iter()
                    .rev()
                    .find_map(|scope| match &*scope.binding {
                        OriginBinding::TopoGraph { sources, .. } => sources
                            .iter()
                            .find(|source| source.role == *role)
                            .and_then(|source| source.event.clone()),
                        _ => None,
                    })
            }
        }?;
        (found.revision == identity.revision && self.history.valid(&found)).then_some(found)
    }
    fn event_pin(&mut self, event: &EventBinding) -> Option<TransitionPin> {
        let native = self
            .workspace
            .history
            .iter()
            .find(|native| native.id == event.event_id)?;
        match &native.transition {
            Transition::ReferenceEdit {
                source, actions, ..
            } => Some(TransitionPin::ReferenceEdit {
                event: event.clone(),
                source: source.clone(),
                actions_revision: Operation::actions_revision(source, actions).ok()?,
            }),
            Transition::Fork { source } => Some(TransitionPin::Fork {
                event: event.clone(),
                source: source.clone(),
            }),
            Transition::CloningStep { run_id, step_index } => Some(TransitionPin::CloningStep {
                event: event.clone(),
                run_id: run_id.clone(),
                step_index: *step_index,
                plan_revision: self.run_plan_revision(run_id).ok()?,
            }),
            Transition::Analysis { origin }
                if origin.producer == AnalysisProducer::MaterialDerivation =>
            {
                match &origin.record {
                    AnalysisRecord::Context { revision } => {
                        let (context, hash) = self.materials.get(&event.event_id)?;
                        if context != revision {
                            return None;
                        }
                        Some(TransitionPin::MaterialDerivation {
                            event: event.clone(),
                            context_revision: revision.clone(),
                            record_hash: hash.clone(),
                        })
                    }
                    _ => None,
                }
            }
            _ => None,
        }
    }
    fn retain_trace(
        &mut self,
        trace: ReconstructionTrace,
        root: Root,
    ) -> Result<(usize, Vec<EventBinding>), Error> {
        let mut bindings = HashMap::new();
        let mut aliases = Vec::new();
        let mut parents = BTreeMap::new();
        for (token, identity) in trace.sources() {
            let Some(event) = self.resolve_identity(identity, &root) else {
                continue;
            };
            if !self.admit(&event.event_id) {
                continue;
            }
            let derived = matches!(
                &identity.source,
                SourceHandle::Plan(Source::Step(_))
                    | SourceHandle::MaterialProduct
                    | SourceHandle::ReferenceProduct { .. }
            ) || matches!(&identity.source, SourceHandle::ReferenceEvent(_))
                && trace
                    .snapshot_projections(token)
                    .is_some_and(|projections| !projections.is_empty());
            if derived {
                if let Some(pin) = self.event_pin(&event) {
                    aliases.push((event.event_id.clone(), token, pin));
                }
            } else {
                parents.insert(event.event_id.clone(), event.clone());
            }
            bindings.insert(token, event);
        }
        let index = self.history.traces.len();
        self.history.traces.push(NativeTrace { trace, bindings });
        for (id, start, pin) in aliases {
            // A genuine derived snapshot from this same execution can satisfy a
            // later request without replaying that native recipe a second time.
            self.history.nodes.entry(id).or_insert(Node::Native {
                trace: index,
                start,
                pin,
            });
        }
        Ok((index, parents.into_values().collect()))
    }
    fn replay(&mut self, binding: &EventBinding) -> Result<(Node, Vec<EventBinding>), Error> {
        if !self.history.valid(binding) {
            return Err(Error::Format(
                "Exact retained document/event/revision is missing or stale".into(),
            ));
        }
        let workspace = self.workspace;
        let native = workspace
            .history
            .iter()
            .find(|event| event.id == binding.event_id)
            .unwrap();
        let expected = material::run::resolve(workspace, binding)?;
        match &native.transition {
            Transition::IndependentReference => Ok((Node::Root, vec![])),
            Transition::Fork { source } => {
                let actual = material::run::resolve(workspace, source)?;
                self.budget.retain_model(actual)?;
                if actual != expected {
                    return Err(Error::Format("Native fork source differs from retained product".into()));
                }
                Ok((Node::Fork {
                    source: source.clone(),
                    pin: TransitionPin::Fork { event: binding.clone(), source: source.clone() },
                }, vec![source.clone()]))
            }
            Transition::ReferenceEdit { source, actions, .. } => {
                let actual = material::run::resolve(workspace, source)?;
                self.budget.preflight_native_job(actions.len())?;
                let (product, trace) = Operation::apply_all_with_ancestry(
                    actual, source, actions, self.budget.clone(),
                )?;
                if &product != expected {
                    return Err(Error::Format("Actual native reference edit differs from retained event".into()));
                }
                let trace = trace?;
                let last = actions.len() - 1;
                let start = trace.step(last)
                    .ok_or_else(|| Error::Format("Native edit lacks its exact final action receipt".into()))?
                    .final_snapshot();
                let pin = TransitionPin::ReferenceEdit {
                    event: binding.clone(), source: source.clone(),
                    actions_revision: Operation::actions_revision(source, actions)?,
                };
                let (index, parents) = self.retain_trace(trace, Root::Reference {
                    event: binding.clone(), last,
                })?;
                Ok((Node::Native { trace: index, start, pin }, parents))
            }
            Transition::CloningStep { run_id, step_index } => {
                if let Some(index) = self.run_traces.get(run_id).copied().flatten() {
                    let start = self.history.traces[index].trace.step(*step_index)
                        .ok_or_else(|| Error::Format("Selected native run outcome lacks a trace".into()))?
                        .final_snapshot();
                    let pin = self.event_pin(binding)
                        .ok_or_else(|| Error::Format("Exact native run binding is unavailable".into()))?;
                    return Ok((Node::Native { trace: index, start, pin }, vec![]));
                }
                if self.run_traces.contains_key(run_id) {
                    return Err(Error::Unsupported("Native run receipt retention was refused".into()));
                }
                self.run_traces.insert(run_id.clone(), None);
                let record=self.records.remove(run_id).ok_or_else(||Error::Format("Exact retained native run is unavailable".into()))?;
                let roots = RunBindings {
                    inputs: record.inputs.clone(), products: record.products.clone(),
                };
                self.runs.insert(run_id.clone(), record.plan_revision.clone());
                let target = material::run::resolve(workspace, &record.target)?;
                self.budget.preflight_native_job(record.plan.steps.len())?;
                if self.replay.exhausted() {
                    return Err(Error::Unsupported("Shared QC native scientific work is exhausted".into()));
                }
                let (report, _, trace) = verification::engine::verify_with_tracked_ancestry_budget(
                    &record.plan, target, &mut self.replay, self.budget.clone(),
                );
                if self.replay.exhausted() {
                    return Err(Error::Unsupported("Shared QC native scientific work is exhausted; complete source correspondence is unavailable".into()));
                }
                for product in &record.products {
                    let actual = material::run::resolve(workspace, &product_binding(product))?;
                    if report.product(product.step_index) != Some(actual) {
                        return Err(Error::Format("Fresh native run products differ from exact retained products".into()));
                    }
                }
                let trace = trace?;
                let start = trace.step(*step_index)
                    .ok_or_else(|| Error::Format("Selected native run product has no captured receipt".into()))?
                    .final_snapshot();
                let pin = TransitionPin::CloningStep {
                    event: binding.clone(), run_id: run_id.clone(), step_index: *step_index,
                    plan_revision: record.plan_revision,
                };
                let (index, parents) = self.retain_trace(trace, Root::Run(roots))?;
                self.run_traces.insert(run_id.clone(), Some(index));
                Ok((Node::Native { trace: index, start, pin }, parents))
            }
            Transition::Analysis { origin } => {
                match (&origin.producer, &origin.record) {
                    (AnalysisProducer::MaterialDerivation, AnalysisRecord::Context { revision }) => {
                        let mut origins = material::capsule::documentary_origins(
                            workspace, &[material::run::MaterialSelection::Derived(binding.clone())], &[],
                        )?.materials;
                        let material::run::MaterialOrigin::Derived(capsule) = origins.remove(0) else {
                            return Err(Error::Format("Retained material recipe kind differs".into()));
                        };
                        self.budget.preflight_native_job(1)?;
                        if self.replay.exhausted() {
                            return Err(Error::Unsupported("Shared QC native scientific work is exhausted".into()));
                        }
                        let (_, root, trace) = material::derivation::replay_scoped(
                            &capsule, expected, &mut self.replay, Some(self.budget.clone()),
                        )?;
                        let trace = trace.ok_or_else(|| Error::Format("Requested native material trace is unavailable".into()))?;
                        let OriginBinding::Material { context_revision, record_hash, .. } = &root else {
                            return Err(Error::Format("Actual material replay origin kind differs".into()));
                        };
                        if context_revision != revision {
                            return Err(Error::Format("Typed material context differs from native replay".into()));
                        }
                        let pin = TransitionPin::MaterialDerivation {
                            event: binding.clone(), context_revision: context_revision.clone(),
                            record_hash: record_hash.clone(),
                        };
                        let trace = trace?;
                        let start = trace.step(0)
                            .ok_or_else(|| Error::Format("Material replay has no native output trace".into()))?
                            .final_snapshot();
                        let (index, parents) = self.retain_trace(trace, Root::Material(root))?;
                        Ok((Node::Native { trace: index, start, pin }, parents))
                    }
                    (_, AnalysisRecord::UnretainedPreview { .. }) => Ok((Node::Unknown(
                        "An unretained computational preview supplies no executable native recipe".into(),
                    ), vec![])),
                    _ => Ok((Node::Unknown(
                        "This typed analysis producer has no supported native causal replay for this request".into(),
                    ), vec![])),
                }
            }
        }
    }
}

pub(super) struct Session<'a> {
    history: &'a History,
    lookups: Vec<LookupSession<'a>>,
    targets: HashMap<(usize, SnapshotToken), BTreeMap<String, LookupTargets>>,
}
impl Session<'_> {
    pub fn history(&self, product: &EventBinding) -> HistoryEvidence {
        if !self.history.valid(product) {
            return HistoryEvidence {
                available: false,
                reason: "Exact product document/event/revision is stale or unavailable".into(),
            };
        }
        match self.history.nodes.get(&product.event_id) {
            Some(Node::Unknown(reason)) => HistoryEvidence {
                available: false,
                reason: reason.clone(),
            },
            Some(_) => HistoryEvidence {
                available: true,
                reason: "The selected event has a current typed native producer or reference root; full source correspondence is checked separately"
                    .into(),
            },
            None => HistoryEvidence {
                available: false,
                reason: "Requested native history-event bound was reached".into(),
            },
        }
    }
    fn prepare_targets(
        &mut self,
        index: usize,
        start: SnapshotToken,
        baseline: &EventBinding,
        work: &mut Work,
    ) -> Result<(), String> {
        let key = (index, start);
        if self
            .targets
            .get(&key)
            .is_some_and(|targets| targets.contains_key(&baseline.event_id))
        {
            return Ok(());
        }
        let native = &self.history.traces[index];
        let mut tokens = Vec::new();
        for (token, event) in &native.bindings {
            if !work.charge(1) {
                return Err("Shared causal lookup work was exhausted".into());
            }
            if event == baseline && *token != start {
                // The canonical owner rejects future/nonancestor comparison
                // tokens. Such a token cannot be a reached causal predecessor.
                if native.trace.lookup_targets(start, &[*token]).is_ok() {
                    tokens.push(*token);
                }
            }
        }
        let targets = native
            .trace
            .lookup_targets(start, &tokens)
            .map_err(reason)?;
        self.targets
            .entry(key)
            .or_default()
            .insert(baseline.event_id.clone(), targets);
        Ok(())
    }
    pub fn map_call(
        &mut self,
        product: &EventBinding,
        baseline: &EventBinding,
        call: BiologicalCall,
        work: &mut Work,
    ) -> CallMapping {
        let mut result = CallMapping {
            parent_calls: vec![],
            introduced: false,
            foreign: false,
            complete: true,
            changed: false,
            causes: vec![],
            reason: None,
        };
        if !self.history.valid(product) || !self.history.valid(baseline) {
            result.complete = false;
            result.reason =
                Some("Exact source/product event binding is unavailable or stale".into());
            return result;
        }
        let mut pending = vec![(product.clone(), call, 0usize)];
        let mut visited = BTreeSet::new();
        while let Some((event, call, depth)) = pending.pop() {
            if !work.charge(1) || depth > 64 {
                result.complete = false;
                result.reason = Some("Shared causal history traversal bound was reached".into());
                break;
            }
            if !visited.insert((event.event_id.clone(), call.position, call.reverse)) {
                continue;
            }
            if event == *baseline {
                if !result.parent_calls.contains(&call) {
                    result.parent_calls.push(call);
                }
                continue;
            }
            match self.history.nodes.get(&event.event_id).cloned() {
                Some(Node::Root) => result.foreign = true,
                Some(Node::Fork { source, pin }) => {
                    if !push_cause(&mut result.causes, pin, work) {
                        result.complete = false;
                        result.reason =
                            Some("Native cause retention or work bound was reached".into());
                        break;
                    }
                    pending.push((source, call, depth + 1));
                }
                Some(Node::Unknown(message)) => {
                    result.complete = false;
                    result.reason = Some(message);
                }
                None => {
                    result.complete = false;
                    result.reason =
                        Some("Native history-event bound prevents predecessor traversal".into());
                }
                Some(Node::Native {
                    trace: index,
                    start,
                    pin,
                }) => {
                    if !push_cause(&mut result.causes, pin, work) {
                        result.complete = false;
                        result.reason =
                            Some("Native cause retention or work bound was reached".into());
                        break;
                    }
                    if let Err(message) = self.prepare_targets(index, start, baseline, work) {
                        result.complete = false;
                        result.reason = Some(message);
                        break;
                    }
                    let native = &self.history.traces[index];
                    let targets = &self.targets[&(index, start)][&baseline.event_id];
                    let lookup = &mut self.lookups[index];
                    lookup.restrict_remaining(work.remaining());
                    let before = lookup.remaining();
                    let evidence = lookup.at_snapshot(
                        start,
                        call.position,
                        Some(LookupStop::Targets(targets)),
                    );
                    work.charge(before - lookup.remaining());
                    match evidence {
                        Err(error) => {
                            result.complete = false;
                            result.reason = Some(reason(error));
                        }
                        Ok(evidence) => {
                            result.complete &= evidence.complete;
                            if let Some(message) = evidence.reason {
                                result.reason = Some(message.into());
                            }
                            for origin in evidence.origins {
                                if !path_causes(&mut result.causes, &origin.path, work) {
                                    result.complete = false;
                                    result.reason = Some(
                                        "Native cause retention or work bound was reached".into(),
                                    );
                                    break;
                                }
                                result.changed |= !origin.changes.is_empty();
                                let token = origin.path.last().map(|edge| edge.source_snapshot);
                                let binding = token.and_then(|token| native.bindings.get(&token));
                                // A stopped root may be the initial snapshot
                                // (zero-edge identity); bind that exact start.
                                let binding = binding.or_else(|| {
                                    origin
                                        .path
                                        .is_empty()
                                        .then(|| native.bindings.get(&start))
                                        .flatten()
                                });
                                if let Some(binding) = binding {
                                    if pending.len() >= 64 {
                                        result.complete = false;
                                        result.reason = Some(
                                            "Native origin branching exceeds the 64-state bound"
                                                .into(),
                                        );
                                        break;
                                    }
                                    pending.push((
                                        binding.clone(),
                                        BiologicalCall {
                                            position: origin.position,
                                            reverse: call.reverse ^ origin.reverse,
                                        },
                                        depth + 1,
                                    ));
                                } else {
                                    result.complete = false;
                                    result.reason=Some("Native source snapshot has no exact retained event binding".into());
                                }
                            }
                            for introduced in evidence.introductions {
                                result.introduced = true;
                                if !path_causes(&mut result.causes, &introduced.path, work) {
                                    result.complete = false;
                                    result.reason = Some(
                                        "Native cause retention or work bound was reached".into(),
                                    );
                                    break;
                                }
                            }
                        }
                    }
                }
            }
        }
        result
    }
    pub fn map_feature(
        &mut self,
        product: &EventBinding,
        baseline: &EventBinding,
        feature_index: usize,
        work: &mut Work,
    ) -> FeatureMapping {
        let mut result = FeatureMapping {
            parent_feature_indices: vec![],
            complete: true,
            reason: None,
        };
        if !self.history.valid(product) || !self.history.valid(baseline) {
            result.complete = false;
            result.reason =
                Some("Exact source/product event binding is unavailable or stale".into());
            return result;
        }
        let mut pending = vec![(product.clone(), feature_index, 0usize)];
        let mut visited = BTreeSet::new();
        while let Some((event, index, depth)) = pending.pop() {
            if !work.charge(1) || depth > 64 {
                result.complete = false;
                result.reason = Some("Shared feature-history traversal bound was reached".into());
                break;
            }
            if !visited.insert((event.event_id.clone(), index)) {
                continue;
            }
            if event == *baseline {
                if !result.parent_feature_indices.contains(&index) {
                    result.parent_feature_indices.push(index);
                }
                continue;
            }
            match self.history.nodes.get(&event.event_id).cloned() {
                Some(Node::Root) => {}
                Some(Node::Fork { source, .. }) => pending.push((source, index, depth + 1)),
                Some(Node::Unknown(message)) => {
                    result.complete = false;
                    result.reason = Some(message);
                }
                None => {
                    result.complete = false;
                    result.reason =
                        Some("Native history-event bound prevents feature traversal".into());
                }
                Some(Node::Native {
                    trace: trace_index,
                    start,
                    ..
                }) => {
                    if let Err(message) = self.prepare_targets(trace_index, start, baseline, work) {
                        result.complete = false;
                        result.reason = Some(message);
                        break;
                    }
                    let native = &self.history.traces[trace_index];
                    let targets = &self.targets[&(trace_index, start)][&baseline.event_id];
                    let lookup = &mut self.lookups[trace_index];
                    lookup.restrict_remaining(work.remaining());
                    let before = lookup.remaining();
                    let evidence = lookup.feature_at_snapshot(
                        start,
                        index,
                        Some(LookupStop::Targets(targets)),
                    );
                    work.charge(before - lookup.remaining());
                    match evidence {
                        Err(error) => {
                            result.complete = false;
                            result.reason = Some(reason(error));
                        }
                        Ok(evidence) => {
                            result.complete &= evidence.complete;
                            if let Some(message) = evidence.reason {
                                result.reason = Some(message.into());
                            }
                            for origin in evidence.origins {
                                let token = origin
                                    .path
                                    .last()
                                    .map(|edge| edge.source_snapshot)
                                    .unwrap_or(start);
                                if let Some(binding) = native.bindings.get(&token) {
                                    if pending.len() >= 64 {
                                        result.complete = false;
                                        result.reason = Some(
                                            "Feature origin branching exceeds the 64-state bound"
                                                .into(),
                                        );
                                        break;
                                    }
                                    pending.push((
                                        binding.clone(),
                                        origin.feature_index,
                                        depth + 1,
                                    ));
                                } else {
                                    result.complete = false;
                                    result.reason = Some(
                                        "Native feature source has no exact retained event binding"
                                            .into(),
                                    );
                                }
                            }
                        }
                    }
                }
            }
        }
        result
    }
}
fn push_cause(causes: &mut Vec<TransitionPin>, pin: TransitionPin, work: &mut Work) -> bool {
    if !work.charge(1) {
        return false;
    }
    if causes.contains(&pin) {
        return true;
    }
    if causes.len() >= 64 {
        return false;
    }
    causes.push(pin);
    true
}
fn path_causes(causes: &mut Vec<TransitionPin>, path: &[OccurrenceKey], work: &mut Work) -> bool {
    for edge in path {
        for scope in &edge.scope {
            let pin = match &*scope.binding {
                OriginBinding::Run {
                    selection,
                    plan_revision,
                    ..
                } => TransitionPin::CloningStep {
                    event: selection.event.clone(),
                    run_id: selection.run_id.clone(),
                    step_index: selection.step_index,
                    plan_revision: plan_revision.clone(),
                },
                OriginBinding::Material {
                    product,
                    context_revision,
                    record_hash,
                    ..
                } => TransitionPin::MaterialDerivation {
                    event: product.clone(),
                    context_revision: context_revision.clone(),
                    record_hash: record_hash.clone(),
                },
                OriginBinding::TopoGraph {
                    selection,
                    context_revision,
                    ..
                } => TransitionPin::TopoGraph {
                    record_hash: selection.record_hash.clone(),
                    state_revision: selection.state_revision.clone(),
                    context_revision: context_revision.clone(),
                },
            };
            if !push_cause(causes, pin, work) {
                return false;
            }
        }
    }
    true
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{
        construct_qc::{AnnotationAssumption, annotation_revision},
        documents::history::{DocumentHistory, transition::AnalysisOrigin},
        verification::ancestry::Limits,
    };
    use collagene_model::{Feature, Location, Molecule};
    fn event(history: &DocumentHistory) -> EventBinding {
        EventBinding {
            document_id: history.id.clone(),
            event_id: history.head_event().into(),
            revision: history.head_revision().into(),
        }
    }
    fn declared(molecule: &Molecule) -> AnnotationAssumption {
        AnnotationAssumption {
            feature_index: 0,
            source_revision: molecule.revision(),
            annotation_revision: annotation_revision(&molecule.features[0]).unwrap(),
            rationale: "Declared synthetic CDS coordinates; not biological identity".into(),
        }
    }
    fn fixture(preview_baseline: bool) -> (Workspace, Request, EventBinding) {
        let mut source = Molecule::dna(
            "identical labels",
            "ATGCAACCCTAA",
            collagene_model::Topology::Linear,
        )
        .unwrap();
        source.features.push(Feature {
            kind: "CDS".into(),
            location: Location::range(0, 12),
            qualifiers: vec![],
            needs_review: false,
        });
        let unrelated = DocumentHistory::new("independent-other-owner", source.clone()).unwrap();
        let foreign = event(&unrelated);
        let mut history = DocumentHistory::new("actual-owner", source.clone()).unwrap();
        if preview_baseline {
            let mut result = history.molecule.clone();
            result.name = "Unretained design".into();
            let origin =
                AnalysisOrigin::preview(AnalysisProducer::AssemblyIntent, result.revision());
            history
                .commit_analysis(result, "Actual unretained computed fixture", origin, &[])
                .unwrap();
        }
        let baseline = event(&history);
        let source_cds = declared(&history.molecule);
        history
            .commit(
                vec![Operation::Rename {
                    name: "Actual renamed product".into(),
                }],
                "Actual native rename",
            )
            .unwrap();
        let product = event(&history);
        let request = Request {
            pairs: vec![super::super::CdsPair {
                id: "pair".into(),
                source: baseline,
                product,
                source_cds,
                product_cds: declared(&history.molecule),
                genetic_code: 1,
            }],
            ..Default::default()
        };
        let mut workspace =
            Workspace::from_documents("Native causal literal fixture", vec![source], 0).unwrap();
        workspace.documents.clear();
        workspace.history.clear();
        workspace.revisions.clear();
        history.contribute(&mut workspace).unwrap();
        unrelated.contribute(&mut workspace).unwrap();
        workspace.active_document = history.id.clone();
        workspace.validate().unwrap();
        (workspace, request, foreign)
    }
    #[test]
    fn selected_baseline_stops_before_its_unretained_predecessor_and_shares_one_native_trace() {
        let (workspace, mut request, _) = fixture(true);
        let mut second = request.pairs[0].clone();
        second.id = "second-declared-pair".into();
        request.pairs.push(second);
        let history = History::prepare(
            &workspace,
            &request,
            TraceBudget::new(Limits::default()).unwrap(),
        )
        .unwrap();
        assert_eq!(history.traces.len(), 1);
        let pair = &request.pairs[0];
        let mut session = history.session();
        let mut work = Work::new(10000);
        assert!(!session.history(&pair.source).available);
        let mapping = session.map_call(
            &pair.product,
            &pair.source,
            BiologicalCall {
                position: 4,
                reverse: false,
            },
            &mut work,
        );
        assert!(mapping.complete && !mapping.foreign && !mapping.introduced);
        assert_eq!(
            mapping.parent_calls,
            vec![BiologicalCall {
                position: 4,
                reverse: false
            }]
        );
        let feature = session.map_feature(&pair.product, &pair.source, 0, &mut work);
        assert!(feature.complete);
        assert_eq!(feature.parent_feature_indices, vec![0]);
    }
    #[test]
    fn equal_dna_from_another_event_owner_never_supplies_selected_ancestry() {
        let (workspace, request, foreign) = fixture(false);
        let history = History::prepare(
            &workspace,
            &request,
            TraceBudget::new(Limits::default()).unwrap(),
        )
        .unwrap();
        let mut session = history.session();
        let mut work = Work::new(10000);
        let mapping = session.map_call(
            &request.pairs[0].product,
            &foreign,
            BiologicalCall {
                position: 4,
                reverse: true,
            },
            &mut work,
        );
        assert!(mapping.complete && mapping.foreign);
        assert!(mapping.parent_calls.is_empty());
    }
    #[test]
    fn history_and_shared_lookup_caps_remain_unknown_instead_of_absence() {
        let (workspace, mut request, _) = fixture(false);
        request.limits.max_history_events = 1;
        let history = History::prepare(
            &workspace,
            &request,
            TraceBudget::new(Limits::default()).unwrap(),
        )
        .unwrap();
        assert_eq!(history.history_events(), 1);
        let mut session = history.session();
        let pair = &request.pairs[0];
        let mut work = Work::new(10000);
        let limited = session.map_call(
            &pair.product,
            &pair.source,
            BiologicalCall {
                position: 4,
                reverse: false,
            },
            &mut work,
        );
        assert!(!limited.complete && limited.reason.is_some());
        let mut empty_work = Work::new(1);
        let limited = session.map_call(
            &pair.product,
            &pair.source,
            BiologicalCall {
                position: 4,
                reverse: false,
            },
            &mut empty_work,
        );
        assert!(!limited.complete && limited.reason.is_some());
    }
    #[test]
    fn a_thousand_native_event_preparations_use_an_explicit_agenda() {
        let (mut workspace, mut request, _) = fixture(false);
        let mut histories = DocumentHistory::restore_all(&workspace).unwrap();
        let history = histories
            .iter_mut()
            .find(|history| history.id == request.pairs[0].product.document_id)
            .unwrap();
        // The fixture already has one real rename; add999more distinct events.
        for index in 0..999 {
            history
                .commit(
                    vec![Operation::Rename {
                        name: format!("native-event-{index}"),
                    }],
                    "Actual native name edit",
                )
                .unwrap();
        }
        request.pairs[0].product = event(history);
        request.pairs[0].product_cds = declared(&history.molecule);
        request.limits.max_history_events = 1000;
        workspace.documents.clear();
        workspace.history.clear();
        workspace.revisions.clear();
        for history in histories {
            history.contribute(&mut workspace).unwrap();
        }
        let prepared = History::prepare(
            &workspace,
            &request,
            TraceBudget::new(Limits::default()).unwrap(),
        )
        .unwrap();
        assert_eq!(prepared.history_events(), 1000);
        assert_eq!(prepared.traces.len(), 1000);
    }
}
