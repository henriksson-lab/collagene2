//! Retained cloning runs are historical evidence, never cached verification authority.
use crate::{
    documents::{
        history::transition::Transition,
        transfer::Candidate,
        workspace::{RevisionEvent, Workspace, WorkspaceAttachment, WorkspaceDocument},
    },
    planning::plan::{DesignPlan, Operation, Source},
    verification::{
        self,
        evidence::{Status, VerificationReport},
    },
};
use collagene_model::{Error, Molecule};
use serde::{Deserialize, Serialize};
use std::{collections::BTreeSet, io::Write};
pub const WORKSPACE_ATTACHMENT_KIND: &str = "cloning-run-v1";
pub const MAX_RUN_BYTES: usize = 32 * 1024 * 1024;
pub const MAX_RUNS: usize = 128;
macro_rules! dto { ($($i:item)*) => {$($i)*}; }
dto! {
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature="json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct EventBinding { pub document_id:String, pub event_id:String, pub revision:String }
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature="json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct InputEventBinding { pub input_id:String, pub event:EventBinding }
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature="json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct ProductBinding { pub step_index:usize, pub document_id:String, pub event_id:String, pub revision:String }
#[derive(Debug, Clone, Serialize, Deserialize)]
#[cfg_attr(feature="json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct RecordedAssessment {
    #[serde(deserialize_with="crate::encoding::required_option")]
    #[cfg_attr(feature="json-schema",schemars(required))]
    pub standard_review: Option<crate::cloning::assembly::standards::declaration::DeclaredReview>,
    pub engine_version: String,
    pub status: Status,
    pub findings: Vec<crate::verification::checks::CheckFinding>,
    pub assumptions: Vec<String>,
    pub exhaustive_alternative_search: bool,
}
#[derive(Debug, Clone, Serialize, Deserialize)]
#[cfg_attr(feature="json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct RunRecord { pub schema_version:u32, pub plan_schema:u32, pub run_id:String, pub plan_revision:String, pub plan:DesignPlan, pub target:EventBinding, pub inputs:Vec<InputEventBinding>, pub products:Vec<ProductBinding>, pub assessment:RecordedAssessment }
#[derive(Debug, Clone, Serialize, Deserialize)]
#[cfg_attr(feature="json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct RunSummary { pub run_id:String, pub attachment_hash:String, pub plan_name:String, pub target:EventBinding, pub products:Vec<ProductBinding>, pub recorded_status:Status, pub recorded_engine:String, pub input_count:usize }
#[derive(Debug, Clone, Serialize, Deserialize)]
#[cfg_attr(feature="json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct RunReceipt { pub run_id:String, pub attachment_hash:String, pub products:Vec<ProductBinding>, pub recorded_status:Status, pub recorded_engine:String, pub already_recorded:bool, pub source_workspace_fingerprint:String }
#[derive(Debug, Clone, Serialize, Deserialize)]
#[cfg_attr(feature="json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct RunCommit { #[cfg_attr(feature="json-schema", schemars(with="serde_json::Value"))] pub workspace:Workspace, pub receipt:RunReceipt }
#[derive(Debug, Clone, Serialize, Deserialize)]
#[cfg_attr(feature="json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct ReplayReport { pub run_id:String, pub source_workspace_fingerprint:String, pub recorded_status:Status, pub recorded_engine:String, pub current_report:VerificationReport, pub engine_matches:bool, pub products_match:bool, pub status_matches:bool }
}
fn bad(s: &str) -> Error {
    Error::Format(format!("cloning run: {s}"))
}
fn digest(s: &str) -> bool {
    s.len() == 64
        && s.bytes()
            .all(|c| c.is_ascii_digit() || (b'a'..=b'f').contains(&c))
}
struct LimitedWriter<W> {
    inner: W,
    remaining: usize,
}
impl<W: Write> Write for LimitedWriter<W> {
    fn write(&mut self, b: &[u8]) -> std::io::Result<usize> {
        if b.len() > self.remaining {
            return Err(std::io::Error::other("size limit"));
        }
        self.inner.write_all(b)?;
        self.remaining -= b.len();
        Ok(b.len())
    }
    fn flush(&mut self) -> std::io::Result<()> {
        self.inner.flush()
    }
}
fn encoded<T: Serialize>(v: &T) -> Result<Vec<u8>, Error> {
    let mut out = LimitedWriter {
        inner: Vec::new(),
        remaining: MAX_RUN_BYTES,
    };
    serde_json::to_writer(&mut out, v).map_err(|_| Error::TooLarge)?;
    Ok(out.inner)
}
fn hash<T: Serialize>(v: &T) -> Result<String, Error> {
    let mut out = LimitedWriter {
        inner: blake3::Hasher::new(),
        remaining: crate::documents::workspace::MAX_EXPANDED_BYTES,
    };
    serde_json::to_writer(&mut out, v).map_err(|_| Error::TooLarge)?;
    Ok(out.inner.finalize().to_hex().to_string())
}
/// Hash validated content-addressed descriptors without expanding attachment bytes into JSON integers.
pub fn workspace_fingerprint(w: &Workspace) -> Result<String, Error> {
    w.validate()?;
    hash(&(
        "cloning-workspace-fingerprint-v1",
        &w.name,
        &w.active_document,
        &w.documents,
        w.revisions.keys().collect::<Vec<_>>(),
        &w.history,
        w.candidates
            .iter()
            .map(|c| &c.plan_revision)
            .collect::<Vec<_>>(),
        &w.active_candidate,
        w.attachments
            .iter()
            .map(|(k, a)| (k, &a.kind, &a.name, &a.media_type))
            .collect::<Vec<_>>(),
    ))
}
fn sources(op: &Operation) -> Vec<&Source> {
    crate::planning::plan::Operation::sources(op)
}
fn operation_label(run_id: &str, index: usize) -> String {
    format!("Cloning run {run_id} step {}", index + 1)
}
fn resolve<'a>(w: &'a Workspace, b: &EventBinding) -> Result<&'a Molecule, Error> {
    let e = w
        .history
        .iter()
        .find(|e| e.id == b.event_id)
        .ok_or_else(|| bad("event unavailable"))?;
    if e.document_id != b.document_id || e.revision != b.revision {
        return Err(bad("event binding mismatch"));
    }
    w.revisions
        .get(&b.revision)
        .ok_or_else(|| bad("revision unavailable"))
}
fn validate_plan(plan: &DesignPlan) -> Result<(), Error> {
    crate::material::capsule::plan_preflight(plan)?;
    verification::sources::validate_reference_steps(plan)?;
    encoded(plan)?;
    Candidate::new(plan.clone()).validate()?;
    for (i, op) in plan.steps.iter().enumerate() {
        if matches!(
            op,
            Operation::DeclareSourceModifications {
                template: Source::Step(_),
                ..
            } | Operation::DeclareSourceTermini {
                template: Source::Step(_),
                ..
            }
        ) {
            return Err(bad(
                "source generation declarations require an explicit input source",
            ));
        }
        let refs = sources(op);
        if refs.len()
            > if matches!(op, Operation::Ligate { .. }) {
                1000
            } else {
                128
            }
        {
            return Err(Error::TooLarge);
        }
        for s in refs {
            match s {
                Source::Input(id) if plan.inputs.iter().any(|v| &v.id == id) => (),
                Source::Step(n) if *n < i => (),
                _ => return Err(bad("unavailable or forward operation source")),
            }
        }
    }
    Ok(())
}
fn identity(r: &RunRecord) -> Result<String, Error> {
    hash(&(
        "cloning-run-v1",
        r.schema_version,
        r.plan_schema,
        &r.plan_revision,
        &r.plan,
        &r.target,
        &r.inputs,
        r.products
            .iter()
            .map(|p| (&p.step_index, &p.revision))
            .collect::<Vec<_>>(),
        &r.assessment,
    ))
}
fn validate_texts(values: &[String], limit: usize) -> Result<(), Error> {
    if values.len() > limit || values.iter().any(|v| v.len() > 65_536) {
        return Err(Error::TooLarge);
    }
    Ok(())
}
fn validate_record(r: &RunRecord) -> Result<(), Error> {
    if r.schema_version != 2
        || r.plan_schema != crate::documents::workspace::PLAN_SCHEMA_VERSION
        || !digest(&r.run_id)
        || r.assessment.engine_version.is_empty()
        || r.assessment.engine_version.len() > 128
        || r.assessment.status == Status::Invalid
    {
        return Err(bad("invalid schema or assessment"));
    }
    if r.assessment.findings.len() > 10_000 {
        return Err(Error::TooLarge);
    }
    match (&r.plan.standard, &r.assessment.standard_review) {
        (Some(declaration), Some(review)) => {
            review.validate_binding(&r.plan_revision, declaration)?
        }
        (None, None) => {}
        _ => return Err(bad("declared standard and retained advisory review differ")),
    }
    for finding in &r.assessment.findings {
        finding.validate()?;
    }
    validate_texts(&r.assessment.assumptions, 1_000)?;
    if !r
        .assessment
        .engine_version
        .bytes()
        .all(|b| b.is_ascii_graphic())
    {
        return Err(bad("invalid recorded engine"));
    }
    validate_plan(&r.plan)?;
    if r.plan_revision != crate::planning::plan::DesignPlan::revision(&r.plan)
        || r.target.revision != r.plan.target_revision
        || r.inputs.len() != r.plan.inputs.len()
        || r.products.len() != r.plan.steps.len()
    {
        return Err(bad("plan or binding count mismatch"));
    }
    for (input, b) in r.plan.inputs.iter().zip(&r.inputs) {
        if input.id != b.input_id || input.revision != b.event.revision {
            return Err(bad("input bindings must follow plan input order"));
        }
    }
    for (i, p) in r.products.iter().enumerate() {
        if i != p.step_index || !digest(&p.revision) {
            return Err(bad("invalid product binding"));
        }
    }
    if identity(r)? != r.run_id {
        return Err(bad("run identity mismatch"));
    }
    Ok(())
}
pub fn save_run_record(r: &RunRecord) -> Result<Vec<u8>, Error> {
    validate_record(r)?;
    encoded(r)
}
pub fn load_run_record(bytes: &[u8]) -> Result<RunRecord, Error> {
    if bytes.len() > MAX_RUN_BYTES {
        return Err(Error::TooLarge);
    }
    let r: RunRecord = serde_json::from_slice(bytes).map_err(|e| bad(&e.to_string()))?;
    let raw: serde_json::Value = serde_json::from_slice(bytes).map_err(|e| bad(&e.to_string()))?;
    if raw != serde_json::to_value(&r).map_err(|e| bad(&e.to_string()))? {
        return Err(bad("unknown or ignored record fields"));
    }
    validate_record(&r)?;
    Ok(r)
}
fn parents(r: &RunRecord, i: usize) -> Result<Vec<String>, Error> {
    let mut seen = BTreeSet::new();
    let mut out = Vec::new();
    for s in sources(&r.plan.steps[i]) {
        let id = match s {
            Source::Input(id) => {
                &r.inputs
                    .iter()
                    .find(|b| &b.input_id == id)
                    .ok_or_else(|| bad("missing input"))?
                    .event
                    .event_id
            }
            Source::Step(n) => &r.products[*n].event_id,
        };
        if seen.insert(id.clone()) {
            out.push(id.clone())
        }
    }
    // Saved graph recipes can be the entire native input, with no current
    // PlanInput operands. Preserve their exact retained event closure rather
    // than pretending that a graph with no explicit Source has no parents.
    let graph_request = match &r.plan.steps[i] {
        Operation::Topo { request, .. } | Operation::CheckedTopo { request, .. } => Some(request),
        Operation::TopoActions { request, .. } => Some(&request.preparation),
        _ => None,
    };
    if let Some(request) = graph_request {
        for material in &request.recipe.materials {
            if let crate::cloning::topo::MaterialRecipe::SavedGraph { origin, .. } = material {
                let record = origin.validate_shallow()?;
                for binding in crate::cloning::topo::graph_history::source_bindings(&record)? {
                    if seen.insert(binding.event_id.clone()) {
                        if out.len() == 128 {
                            return Err(Error::TooLarge);
                        }
                        out.push(binding.event_id);
                    }
                }
            }
        }
    }
    Ok(out)
}
fn crossrefs(w: &Workspace, r: &RunRecord) -> Result<(), Error> {
    crate::material::capsule::validate_plan_origins(w, &r.plan)?;
    crate::material::capsule::validate_input_origins(&r.plan, &r.inputs)?;
    resolve(w, &r.target)?;
    for (b, input) in r.inputs.iter().zip(&r.plan.inputs) {
        let m = resolve(w, &b.event)?;
        if serde_json::to_value(m).map_err(|e| bad(&e.to_string()))?
            != serde_json::to_value(&input.molecule).map_err(|e| bad(&e.to_string()))?
        {
            return Err(bad("input snapshot mismatch"));
        }
    }
    let mut docs = BTreeSet::new();
    for (i, p) in r.products.iter().enumerate() {
        if !docs.insert(&p.document_id) {
            return Err(bad("products require separate documents"));
        }
        resolve(
            w,
            &EventBinding {
                document_id: p.document_id.clone(),
                event_id: p.event_id.clone(),
                revision: p.revision.clone(),
            },
        )?;
        let e = w.history.iter().find(|e| e.id == p.event_id).unwrap();
        if e.parents != parents(r, i)?
            || e.operation != operation_label(&r.run_id, i)
            || e.transition
                != (Transition::CloningStep {
                    run_id: r.run_id.clone(),
                    step_index: i,
                })
        {
            return Err(bad("product event provenance mismatch"));
        }
    }
    Ok(())
}
/// Structural archive validation only; deliberately does not replay scientific calculations.
pub fn validate_workspace_runs(w: &Workspace) -> Result<(), Error> {
    let mut ids = BTreeSet::new();
    for a in w
        .attachments
        .values()
        .filter(|a| a.kind == WORKSPACE_ATTACHMENT_KIND)
    {
        if ids.len() >= MAX_RUNS {
            return Err(Error::TooLarge);
        }
        if a.media_type != "application/json" {
            return Err(bad("run attachment media type"));
        }
        let r = load_run_record(&a.bytes)?;
        if !ids.insert(r.run_id.clone()) {
            return Err(bad("duplicate run ID"));
        }
        crossrefs(w, &r)?
    }
    crate::material::physical_evidence::validate_workspace(w)?;
    Ok(())
}
pub fn list_runs(w: &Workspace) -> Result<Vec<RunSummary>, Error> {
    w.validate()?;
    w.attachments
        .iter()
        .filter(|(_, a)| a.kind == WORKSPACE_ATTACHMENT_KIND)
        .map(|(h, a)| {
            let r = load_run_record(&a.bytes)?;
            Ok(RunSummary {
                run_id: r.run_id,
                attachment_hash: h.clone(),
                plan_name: r.plan.name,
                target: r.target,
                products: r.products,
                recorded_status: r.assessment.status,
                recorded_engine: r.assessment.engine_version,
                input_count: r.inputs.len(),
            })
        })
        .collect()
}
pub fn load_run(w: &Workspace, id: &str) -> Result<RunRecord, Error> {
    w.validate()?;
    for a in w
        .attachments
        .values()
        .filter(|a| a.kind == WORKSPACE_ATTACHMENT_KIND)
    {
        let r = load_run_record(&a.bytes)?;
        if r.run_id == id {
            return Ok(r);
        }
    }
    Err(bad("run unavailable"))
}
fn receipt(r: &RunRecord, h: String, fingerprint: String, already_recorded: bool) -> RunReceipt {
    RunReceipt {
        run_id: r.run_id.clone(),
        attachment_hash: h,
        products: r.products.clone(),
        recorded_status: r.assessment.status,
        recorded_engine: r.assessment.engine_version.clone(),
        already_recorded,
        source_workspace_fingerprint: fingerprint,
    }
}
pub fn record_run(
    w: &Workspace,
    plan: &DesignPlan,
    target_binding: EventBinding,
    input_bindings: Vec<InputEventBinding>,
    allow_incomplete: bool,
) -> Result<RunCommit, Error> {
    crate::material::capsule::plan_preflight(plan)?;
    crate::material::capsule::validate_plan_origins(w, plan)?;
    crate::material::capsule::validate_input_origins(plan, &input_bindings)?;
    let fingerprint = workspace_fingerprint(w)?;
    validate_plan(plan)?;
    let target = resolve(w, &target_binding)?;
    if target.revision() != plan.target_revision {
        return Err(bad("target revision mismatch"));
    }
    if input_bindings.len() != plan.inputs.len() {
        return Err(bad("input binding count mismatch"));
    }
    let mut inputs = Vec::new();
    for input in &plan.inputs {
        let matches: Vec<_> = input_bindings
            .iter()
            .filter(|b| b.input_id == input.id)
            .collect();
        if matches.len() != 1 {
            return Err(bad("missing or duplicate input binding"));
        }
        let b = matches[0];
        if resolve(w, &b.event)?.revision() != input.revision {
            return Err(bad("source revision mismatch"));
        }
        inputs.push(b.clone());
    }
    let mut report = verification::engine::verify(plan, target);
    report.validate_documentary()?;
    report.validate_standard_for_plan(plan)?;
    let physical_graphs = std::mem::take(&mut report.physical_processing_reports);
    let reference_graphs = std::mem::take(&mut report.reference_step_reports);
    if report.steps.len() != plan.steps.len()
        || report.steps.iter().any(|step| {
            step.state != crate::verification::evidence::StepState::Reconstructed
                || step.product.is_none()
        })
        || !report
            .final_product()
            .is_some_and(|p| verification::engine::plan_target_matches(plan, p, target))
        || report.status == Status::Invalid
        || (report.status == Status::Incomplete && !allow_incomplete)
    {
        return Err(bad(
            "capture requires a complete reconstruction matching the independent target",
        ));
    }
    let mut r = RunRecord {
        schema_version: 2,
        plan_schema: crate::documents::workspace::PLAN_SCHEMA_VERSION,
        run_id: String::new(),
        plan_revision: report.plan_revision.clone(),
        plan: plan.clone(),
        target: target_binding,
        inputs,
        products: report
            .steps
            .iter()
            .map(|step| {
                let molecule = report
                    .product(step.step_index)
                    .ok_or_else(|| bad("complete run lost a reconstructed step product"))?;
                Ok(ProductBinding {
                    step_index: step.step_index,
                    document_id: String::new(),
                    event_id: String::new(),
                    revision: molecule.revision(),
                })
            })
            .collect::<Result<Vec<_>, Error>>()?,
        assessment: RecordedAssessment {
            standard_review: report.standard_review,
            engine_version: report.engine_version,
            status: report.status,
            findings: report.findings,
            assumptions: report.assumptions,
            exhaustive_alternative_search: report.exhaustive_alternative_search,
        },
    };
    r.run_id = identity(&r)?;
    for (h, a) in &w.attachments {
        if a.kind == WORKSPACE_ATTACHMENT_KIND {
            let existing = load_run_record(&a.bytes)?;
            if existing.run_id == r.run_id {
                return Ok(RunCommit {
                    workspace: w.clone(),
                    receipt: receipt(&existing, h.clone(), fingerprint, true),
                });
            }
        }
    }
    let mut next = w.clone();
    for step in report.steps {
        let i = step.step_index;
        let m = step
            .product
            .ok_or_else(|| bad("complete run lost a reconstructed step product"))?;
        let base = format!("run-{}-{}", r.run_id, i + 1);
        let mut id = base.clone();
        let mut n = 0;
        while next.documents.iter().any(|d| d.id == id) {
            n += 1;
            id = format!("{base}-{n}")
        }
        r.products[i].document_id = id.clone();
        let event = RevisionEvent::new(
            id.clone(),
            m.revision(),
            parents(&r, i)?,
            operation_label(&r.run_id, i),
            Transition::CloningStep {
                run_id: r.run_id.clone(),
                step_index: i,
            },
        )?;
        r.products[i].event_id = event.id.clone();
        next.documents.push(WorkspaceDocument {
            id,
            head: m.revision(),
            head_event: Some(event.id.clone()),
            undo: vec![],
            redo: vec![],
            undo_events: vec![],
            redo_events: vec![],
        });
        next.revisions.insert(m.revision(), m);
        next.history.push(event);
    }
    if let Some(attachment) = crate::material::reference_evidence::attachment(&r, reference_graphs)?
    {
        let hash = blake3::hash(&attachment.bytes).to_hex().to_string();
        next.attachments.insert(hash, attachment);
    }
    if let Some(attachment) = crate::material::physical_evidence::attachment(&r, physical_graphs)? {
        let hash = blake3::hash(&attachment.bytes).to_hex().to_string();
        next.attachments.insert(hash, attachment);
    }

    let bytes = save_run_record(&r)?;
    let h = blake3::hash(&bytes).to_hex().to_string();
    next.attachments.insert(
        h.clone(),
        WorkspaceAttachment {
            kind: WORKSPACE_ATTACHMENT_KIND.into(),
            name: format!("{}.cloning-run.json", r.run_id),
            media_type: "application/json".into(),
            bytes,
        },
    );
    next.validate()?;
    Ok(RunCommit {
        workspace: next,
        receipt: receipt(&r, h, fingerprint, false),
    })
}
pub fn replay_run(w: &Workspace, id: &str) -> Result<ReplayReport, Error> {
    let fingerprint = workspace_fingerprint(w)?;
    let r = load_run(w, id)?;
    let report = verification::engine::verify(&r.plan, resolve(w, &r.target)?);
    report.validate_documentary()?;
    report.validate_standard_for_plan(&r.plan)?;
    let products_match = report.steps.len() == r.products.len()
        && report.steps.iter().zip(&r.products).all(|(step, binding)| {
            step.state == crate::verification::evidence::StepState::Reconstructed
                && step
                    .product
                    .as_ref()
                    .is_some_and(|product| product.revision() == binding.revision)
        });
    Ok(ReplayReport {
        run_id: r.run_id,
        source_workspace_fingerprint: fingerprint,
        recorded_status: r.assessment.status,
        status_matches: report.status == r.assessment.status,
        engine_matches: report.engine_version == r.assessment.engine_version,
        recorded_engine: r.assessment.engine_version,
        products_match,
        current_report: report,
    })
}
impl RunCommit {
    pub fn validate(&self) -> Result<(), Error> {
        self.workspace.validate()?;
        let r = load_run(&self.workspace, &self.receipt.run_id)?;
        let expected = receipt(
            &r,
            self.receipt.attachment_hash.clone(),
            self.receipt.source_workspace_fingerprint.clone(),
            self.receipt.already_recorded,
        );
        if !digest(&self.receipt.source_workspace_fingerprint)
            || self
                .workspace
                .attachments
                .get(&self.receipt.attachment_hash)
                .is_none_or(|a| {
                    a.kind != WORKSPACE_ATTACHMENT_KIND
                        || load_run_record(&a.bytes).map_or(true, |v| v.run_id != r.run_id)
                })
            || encoded(&expected)? != encoded(&self.receipt)?
        {
            return Err(bad("invalid receipt"));
        }
        Ok(())
    }
}
impl ReplayReport {
    /// Bind an untrusted transport response to retained evidence without rerunning science.
    pub fn validate_for_workspace(&self, workspace: &Workspace) -> Result<(), Error> {
        self.validate()?;
        if self.current_report.engine_version != verification::ENGINE_VERSION {
            return Err(bad(
                "fresh replay does not match the current scientific engine",
            ));
        }
        if self.source_workspace_fingerprint != workspace_fingerprint(workspace)? {
            return Err(bad("replay workspace changed"));
        }
        let record = load_run(workspace, &self.run_id)?;
        self.current_report
            .validate_standard_for_plan(&record.plan)?;
        let products_match = self.current_report.steps.len() == record.products.len()
            && self
                .current_report
                .steps
                .iter()
                .zip(&record.products)
                .all(|(step, binding)| {
                    step.state == crate::verification::evidence::StepState::Reconstructed
                        && step
                            .product
                            .as_ref()
                            .is_some_and(|product| product.revision() == binding.revision)
                });
        if self.recorded_status != record.assessment.status
            || self.recorded_engine != record.assessment.engine_version
            || self.products_match != products_match
            || self.current_report.plan_revision != record.plan_revision
            || self.current_report.target_revision != record.target.revision
        {
            return Err(bad("replay evidence binding mismatch"));
        }
        Ok(())
    }
    pub fn validate(&self) -> Result<(), Error> {
        if !digest(&self.run_id)
            || !digest(&self.source_workspace_fingerprint)
            || self.recorded_engine.is_empty()
            || self.engine_matches != (self.recorded_engine == self.current_report.engine_version)
            || self.status_matches != (self.recorded_status == self.current_report.status)
            || self.recorded_status == Status::Invalid
        {
            return Err(bad("invalid replay report"));
        }
        let report = &self.current_report;
        report.validate_documentary()?;
        if !digest(&report.plan_revision)
            || !digest(&report.target_revision)
            || report.engine_version.is_empty()
            || report.engine_version.len() > 128
            || report.steps.len() > 1000
            || self.recorded_engine.len() > 128
            || (self.products_match && report.steps.is_empty())
        {
            return Err(bad("invalid current replay assessment"));
        }

        validate_texts(&report.assumptions, 1_000)?;
        for product in report.steps.iter().filter_map(|step| step.product.as_ref()) {
            product.validate()?;
        }
        hash(self)?;
        Ok(())
    }
}
