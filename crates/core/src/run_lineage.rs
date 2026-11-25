//! Retained cloning runs are historical evidence, never cached verification authority.
use crate::{
    Error, Molecule,
    transfer::Candidate,
    verification::{self, DesignPlan, Operation, Source, Status, VerificationReport},
    workspace::{RevisionEvent, Workspace, WorkspaceAttachment, WorkspaceDocument},
};
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
#[derive(Debug, Clone, Serialize, Deserialize)]
#[cfg_attr(feature="json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct InputEventBinding { pub input_id:String, pub event:EventBinding }
#[derive(Debug, Clone, Serialize, Deserialize)]
#[cfg_attr(feature="json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct ProductBinding { pub step_index:usize, pub document_id:String, pub event_id:String, pub revision:String }
#[derive(Debug, Clone, Serialize, Deserialize)]
#[cfg_attr(feature="json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct RecordedAssessment { pub engine_version:String, pub status:Status, pub findings:Vec<String>, pub assumptions:Vec<String>, pub exhaustive_alternative_search:bool }
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
        remaining: crate::workspace::MAX_EXPANDED_BYTES,
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
    verification::operation_sources(op)
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
    verification::validate_reference_steps(plan)?;
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
        if refs.len() > 128 {
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
    if r.schema_version != 1
        || !matches!(r.plan_schema, 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10)
        || !digest(&r.run_id)
        || r.assessment.engine_version.is_empty()
        || r.assessment.engine_version.len() > 128
        || r.assessment.status == Status::Invalid
    {
        return Err(bad("invalid schema or assessment"));
    }
    if r.plan_schema < 10
        && (r.assessment.engine_version.ends_with("/science-15")
            || r.plan
                .steps
                .iter()
                .any(|op| matches!(op, Operation::SelectedGoldenGate { .. })))
    {
        return Err(bad(
            "Key-bound Golden Gate requires workspace plan schema 10",
        ));
    }
    if r.plan_schema < 9
        && (r.assessment.engine_version.ends_with("/science-14")
            || r.plan
                .steps
                .iter()
                .any(|op| matches!(op, Operation::CatalogDigest { .. })))
    {
        return Err(bad(
            "Catalog-pinned digest requires workspace plan schema 9",
        ));
    }
    if r.plan_schema < 8
        && (r.assessment.engine_version.ends_with("/science-13")
            || r.plan
                .steps
                .iter()
                .any(|op| matches!(op, Operation::DeclareSourceTermini { .. })))
    {
        return Err(bad(
            "Declared source termini require workspace plan schema 8",
        ));
    }
    if r.plan_schema < 7
        && (r.assessment.engine_version.ends_with("/science-12")
            || r.plan
                .steps
                .iter()
                .any(|op| matches!(op, Operation::CheckedTopo { .. })))
    {
        return Err(bad(
            "Explicit TOPO target acceptance requires workspace plan schema 7",
        ));
    }
    if r.plan_schema < 6
        && ((r.assessment.engine_version.ends_with("/science-10")
            || r.assessment.engine_version.ends_with("/science-11"))
            || r.plan
                .steps
                .iter()
                .any(|op| matches!(op, Operation::Topo { .. })))
    {
        return Err(bad(
            "TOPO-aware reconstruction requires workspace plan schema 6",
        ));
    }
    if r.plan_schema < 5
        && (r.assessment.engine_version.ends_with("/science-9")
            || r.plan
                .steps
                .iter()
                .any(|op| matches!(op, Operation::RecombineReference { .. })))
    {
        return Err(bad(
            "reference recombination requires workspace plan schema 5",
        ));
    }
    if r.plan_schema < 4
        && (r.assessment.engine_version.ends_with("/science-8")
            || r.plan.steps.iter().any(|op| {
                matches!(
                    op,
                    Operation::DeclareSourceModifications { .. } | Operation::CheckedDigest { .. }
                )
            }))
    {
        return Err(bad(
            "tracked modification semantics require workspace plan schema 4",
        ));
    }
    validate_texts(&r.assessment.findings, 10_000)?;
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
    if r.plan_revision != verification::plan_revision(&r.plan)
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
    Ok(out)
}
fn crossrefs(w: &Workspace, r: &RunRecord) -> Result<(), Error> {
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
        if e.parents != parents(r, i)? || e.operation != operation_label(&r.run_id, i) {
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
    let report = verification::verify(plan, target);
    if report.products.len() != plan.steps.len()
        || !report
            .products
            .last()
            .is_some_and(|p| verification::plan_target_matches(plan, p, target))
        || report.status == Status::Invalid
        || (report.status == Status::Incomplete && !allow_incomplete)
    {
        return Err(bad(
            "capture requires a complete reconstruction matching the independent target",
        ));
    }
    let mut r = RunRecord {
        schema_version: 1,
        plan_schema: crate::workspace::PLAN_SCHEMA_VERSION,
        run_id: String::new(),
        plan_revision: report.plan_revision.clone(),
        plan: plan.clone(),
        target: target_binding,
        inputs,
        products: report
            .products
            .iter()
            .enumerate()
            .map(|(i, m)| ProductBinding {
                step_index: i,
                document_id: String::new(),
                event_id: String::new(),
                revision: m.revision(),
            })
            .collect(),
        assessment: RecordedAssessment {
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
    for (i, m) in report.products.into_iter().enumerate() {
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
    let report = verification::verify(&r.plan, resolve(w, &r.target)?);
    let products_match = report.products.len() == r.products.len()
        && report
            .products
            .iter()
            .zip(&r.products)
            .all(|(m, b)| m.revision() == b.revision);
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
        let products_match = self.current_report.products.len() == record.products.len()
            && self
                .current_report
                .products
                .iter()
                .zip(&record.products)
                .all(|(m, p)| m.revision() == p.revision);
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
        if !digest(&report.plan_revision)
            || !digest(&report.target_revision)
            || report.engine_version.is_empty()
            || report.engine_version.len() > 128
            || report.products.len() > 1000
            || self.recorded_engine.len() > 128
            || (self.products_match && report.products.is_empty())
        {
            return Err(bad("invalid current replay assessment"));
        }
        validate_texts(&report.findings, 10_000)?;
        validate_texts(&report.assumptions, 1_000)?;
        for product in &report.products {
            product.validate()?;
        }
        hash(self)?;
        Ok(())
    }
}
