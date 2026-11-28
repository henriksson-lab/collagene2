//! Private sequential batch computation and immutable historical summaries.
//! No decoded record/session flag can authorize fresh science or product release.
use crate::{
    documents::workspace::{self, Workspace, WorkspaceAttachment},
    material::lineage::{self as run_lineage, EventBinding, InputEventBinding},
    planning::combinatorial::{
        self as design, ChoiceSelection, Design, SelectedConstruct, SelectedMatrix,
    },
    verification::{self, evidence::Status},
};
use collagene_model::Error;
use serde::{Deserialize, Serialize};
use std::{collections::BTreeSet, io};
pub const KIND: &str = "combinatorial-batch-record-v3";
pub const ENGINE_VERSION: &str = verification::ENGINE_VERSION;
const CONTRACT: &str = "combinatorial-batch-4";
pub const MAX_RECORD_BYTES: usize = 32 * 1024 * 1024;
pub const MAX_RECORDS: usize = 128;
/// Aggregate declared DNA-search/repeat work, not a measured CPU-time claim.
pub const MAX_BATCH_WORK: usize = 2_000_000_000;
const MAX_ROW_BYTES: usize = 128 * 1024;
macro_rules! dto {($($t:tt)*)=>{#[derive(Clone,Debug,Serialize,Deserialize)]#[cfg_attr(feature="json-schema",derive(schemars::JsonSchema))]#[serde(deny_unknown_fields)]$($t)*};}
dto! {pub enum RowRequest {
    Ready(SelectedConstruct),
    /// Caller-reported editable projection error, explicitly not evaluated.
    DraftError{id:String,name:String,choices:Vec<ChoiceSelection>,message:String},
}}
dto! {pub struct Request {pub design:Option<Design>,pub rows:Vec<RowRequest>}}
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
pub enum RowState {
    DraftError,
    Error,
    Evaluated,
}
dto! {pub struct TailingProfileEvidence {pub profile:crate::cloning::ligation::tailing_recipe::Profile,pub status:crate::cloning::ligation::tailing_recipe::ProfileStatus,pub findings:Vec<String>}}
dto! {pub struct RowEvidence {
    pub index:usize,pub construct_id:String,pub name:String,pub state:RowState,
    pub status:Option<Status>,pub error:Option<String>,pub combination_id:Option<String>,pub declaration_revision:Option<String>,
    pub input_bindings:Vec<InputEventBinding>,pub target_binding:Option<EventBinding>,
    pub plan_revision:Option<String>,pub target_revision:Option<String>,pub final_product_revision:Option<String>,
    pub completed_steps:usize,pub planned_steps:usize,pub preparation_complete:bool,
    #[serde(deserialize_with="crate::encoding::required_option")]#[cfg_attr(feature="json-schema",schemars(required))]pub tailing_profile:Option<TailingProfileEvidence>,
    pub findings:Vec<String>,pub assumptions:Vec<String>,pub exhaustive_alternative_search:Option<bool>,
}}
dto! {pub struct Progress {pub context_revision:String,pub next_index:usize,pub total_rows:usize,pub finished:bool,pub row:Option<RowEvidence>}}
dto! {pub struct Record {pub schema_version:u16,pub context_revision:String,pub engine_version:String,pub request:Request,pub rows:Vec<RowEvidence>}}
dto! {pub struct Commit {#[cfg_attr(feature="json-schema",schemars(with="serde_json::Value"))]pub workspace:Workspace,pub context_revision:String,pub capture_revision:String,pub source_workspace_fingerprint:String,pub attachment_hash:String}}
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
pub struct Summary {
    pub attachment_hash: String,
    pub context_revision: String,
    pub design_name: String,
    pub total_rows: usize,
    pub evaluated_rows: usize,
    pub error_rows: usize,
    pub recorded_engine: String,
}
dto! {pub struct Replayed {pub record_hash:String,pub row_index:usize,pub recorded:RowEvidence,pub current:RowEvidence,pub current_engine:String}}

/// Not serializable or constructible from an imported receipt. One owned frame,
/// one request and bounded summaries; full per-row products/plans are dropped.
pub struct Session {
    workspace: Workspace,
    request: Request,
    context_revision: String,
    rows: Vec<RowEvidence>,
    work_used: usize,
    work_exhausted: bool,
}
fn bad(message: &str) -> Error {
    Error::Format(format!("Combinatorial batch: {message}"))
}
fn attachment_versions(workspace: &Workspace) -> Result<(), Error> {
    if workspace
        .attachments
        .values()
        .any(|asset| asset.kind.starts_with("combinatorial-batch-record-") && asset.kind != KIND)
    {
        return Err(bad("unsupported batch attachment version"));
    }
    Ok(())
}
fn digest(value: &str) -> bool {
    value.len() == 64
        && value
            .bytes()
            .all(|b| b.is_ascii_digit() || (b'a'..=b'f').contains(&b))
}
fn identifier(value: &str) -> bool {
    !value.is_empty()
        && value.len() <= 48
        && value
            .bytes()
            .all(|b| b.is_ascii_alphanumeric() || matches!(b, b'_' | b'-' | b'.'))
}
fn label(value: &str) -> bool {
    !value.trim().is_empty() && value.len() <= 256 && !value.chars().any(char::is_control)
}
struct Counter {
    bytes: usize,
    limit: usize,
}
impl io::Write for Counter {
    fn write(&mut self, bytes: &[u8]) -> io::Result<usize> {
        self.bytes = self
            .bytes
            .checked_add(bytes.len())
            .filter(|n| *n <= self.limit)
            .ok_or_else(|| io::Error::other("batch summary byte limit"))?;
        Ok(bytes.len())
    }
    fn flush(&mut self) -> io::Result<()> {
        Ok(())
    }
}
fn count(value: &impl Serialize, limit: usize) -> Result<usize, Error> {
    let mut counter = Counter { bytes: 0, limit };
    serde_json::to_writer(&mut counter, value).map_err(|_| Error::TooLarge)?;
    Ok(counter.bytes)
}
fn row_identity(row: &RowRequest) -> (&str, &str) {
    match row {
        RowRequest::Ready(row) => (&row.id, &row.name),
        RowRequest::DraftError { id, name, .. } => (id, name),
    }
}
impl Request {
    pub fn from_matrix(matrix: SelectedMatrix) -> Self {
        Self {
            design: Some(matrix.design),
            rows: matrix
                .constructs
                .into_iter()
                .map(RowRequest::Ready)
                .collect(),
        }
    }
    pub fn validate(&self) -> Result<(), Error> {
        if self.rows.is_empty() || self.rows.len() > design::MAX_SELECTED {
            return Err(bad("select 1–256 batch rows"));
        }
        if let Some(design) = &self.design {
            design.validate()?;
        }
        let mut ids = BTreeSet::new();
        for row in &self.rows {
            let (id, name) = row_identity(row);
            if !identifier(id) || !label(name) || !ids.insert(id) {
                return Err(bad("invalid or duplicate batch row identity"));
            }
            let choices = match row {
                RowRequest::Ready(construct) => {
                    if self.design.is_none() {
                        return Err(bad("a ready row requires its typed source design"));
                    }
                    // Detailed semantic row validation belongs to advance. A bad
                    // preparation/profile in one ready row must not hide others.
                    if let Some(overrides) = &construct.overrides {
                        if overrides.preparations.len() > design::MAX_BINS {
                            return Err(Error::TooLarge);
                        }
                        if matches!(&overrides.method, Some(design::Method::Tailing { actions, .. }) if actions.len() > design::MAX_BINS)
                        {
                            return Err(Error::TooLarge);
                        }
                        if let Some(design::Method::Overlap {
                            overlaps,
                            finishing,
                            ..
                        }) = &overrides.method
                        {
                            use crate::cloning::assembly::recipe::OverlapFinishing;
                            let excessive_policy = match finishing {
                                OverlapFinishing::ReferenceOnly => false,
                                OverlapFinishing::HomologyAssembly { policy } => {
                                    policy.overlaps.len() > design::MAX_BINS
                                        || matches!(&policy.resection,crate::cloning::overlap::assembly::Resection::Explicit{cuts} if cuts.len()>design::MAX_BINS)
                                }
                                OverlapFinishing::HomologyAnnealing { policy } => {
                                    policy.overlaps.len() > design::MAX_BINS
                                        || policy.cuts.len() > design::MAX_BINS
                                }
                            };
                            if overlaps.len() > design::MAX_BINS || excessive_policy {
                                return Err(Error::TooLarge);
                            }
                        }
                    }
                    design::declaration_bounds(&construct.target)?;
                    &construct.choices
                }
                RowRequest::DraftError {
                    choices, message, ..
                } => {
                    if message.is_empty() || message.len() > 4096 {
                        return Err(Error::TooLarge);
                    }
                    choices
                }
            };
            if choices.len() > design::MAX_BINS
                || choices
                    .iter()
                    .any(|c| !identifier(&c.bin_id) || !identifier(&c.choice_id))
            {
                return Err(bad("invalid bounded batch choice descriptors"));
            }
        }
        // Reserve enough room for bounded local-error summaries for every row.
        count(self, MAX_RECORD_BYTES - self.rows.len() * 8192 - 4096)?;
        Ok(())
    }
}
pub fn request_revision(request: &Request) -> Result<String, Error> {
    request.validate()?;
    crate::encoding::bounded_json::hash(&(CONTRACT, ENGINE_VERSION, "request", request))
}
pub fn preflight(workspace: &Workspace, request: &Request) -> Result<(), Error> {
    request.validate()?;
    design::workspace_guard(workspace)
}
fn begin_token(fingerprint: &str, request: &Request) -> Result<String, Error> {
    crate::encoding::bounded_json::hash(&(
        CONTRACT,
        ENGINE_VERSION,
        "session",
        fingerprint,
        request,
    ))
}
pub fn begin_revision(workspace: &Workspace, request: &Request) -> Result<String, Error> {
    preflight(workspace, request)?;
    begin_token(&run_lineage::workspace_fingerprint(workspace)?, request)
}
fn capture_token(fingerprint: &str, context: &str) -> Result<String, Error> {
    if !digest(fingerprint) || !digest(context) {
        return Err(bad("invalid capture context"));
    }
    crate::encoding::bounded_json::hash(&(
        CONTRACT,
        ENGINE_VERSION,
        "capture",
        fingerprint,
        context,
    ))
}
pub fn capture_revision(workspace: &Workspace, context: &str) -> Result<String, Error> {
    design::workspace_guard(workspace)?;
    capture_token(&run_lineage::workspace_fingerprint(workspace)?, context)
}
fn abbreviated(value: &str) -> String {
    let mut text = value
        .chars()
        .take(1000)
        .map(|c| if c.is_control() { ' ' } else { c })
        .collect::<String>();
    if value.chars().count() > 1000 {
        text.push_str("… (abbreviated diagnostic; inspect the retained row separately)");
    }
    text
}
fn error_row(
    index: usize,
    row: &RowRequest,
    state: RowState,
    error: String,
    combination_id: Option<String>,
) -> RowEvidence {
    let (id, name) = row_identity(row);
    RowEvidence {
        index,
        construct_id: id.into(),
        name: name.into(),
        state,
        status: None,
        error: Some(abbreviated(&error)),
        combination_id,
        declaration_revision: None,
        input_bindings: vec![],
        target_binding: None,
        plan_revision: None,
        target_revision: None,
        final_product_revision: None,
        completed_steps: 0,
        planned_steps: 0,
        preparation_complete: false,
        tailing_profile: None,
        findings: vec![],
        assumptions: vec![],
        exhaustive_alternative_search: None,
    }
}
fn evaluate_row(
    index: usize,
    row: &RowRequest,
    design: Option<&Design>,
    resolved: Result<design::Resolved, Error>,
) -> RowEvidence {
    let construct = match row {
        RowRequest::DraftError { message, .. } => {
            return error_row(index, row, RowState::DraftError, message.clone(), None);
        }
        RowRequest::Ready(construct) => construct,
    };
    let design = match design {
        Some(design) => design,
        None => {
            return error_row(
                index,
                row,
                RowState::Error,
                "Typed design is unavailable.".into(),
                None,
            );
        }
    };
    let combination = design::combination_key(design, &construct.choices).ok();
    let result = (|| -> Result<RowEvidence, Error> {
        let evaluated = resolved?.evaluate()?;
        let report = evaluated.verification();
        let plan = evaluated.plan();
        if report.findings.len() > 4096
            || report.assumptions.len() > 4096
            || report
                .findings
                .iter()
                .chain(&report.assumptions)
                .any(|s| s.len() > 65_536)
        {
            return Err(Error::TooLarge);
        }
        count(
            &(&report.findings, &report.assumptions),
            MAX_ROW_BYTES - 16 * 1024,
        )?;
        if report.engine_version != verification::ENGINE_VERSION
            || report.target_revision != evaluated.resolved.target.revision()
            || report.plan_revision != crate::planning::plan::DesignPlan::revision(plan)
            || report.products.len() > plan.steps.len()
        {
            return Err(bad("native row report identity is inconsistent"));
        }
        let preparation_complete = evaluated.preparation_complete();
        let final_product_revision =
            if !plan.steps.is_empty() && report.products.len() == plan.steps.len() {
                report.products.last().map(|m| m.revision())
            } else {
                None
            };
        let evidence = RowEvidence {
            index,
            construct_id: construct.id.clone(),
            name: construct.name.clone(),
            state: RowState::Evaluated,
            status: Some(report.status),
            error: None,
            combination_id: Some(evaluated.resolved.combination_id.clone()),
            declaration_revision: Some(evaluated.resolved.declaration_revision.clone()),
            input_bindings: evaluated.resolved.bindings.clone(),
            target_binding: evaluated.resolved.target_binding.clone(),
            plan_revision: Some(report.plan_revision.clone()),
            target_revision: Some(report.target_revision.clone()),
            final_product_revision,
            completed_steps: report.products.len(),
            planned_steps: plan.steps.len(),
            preparation_complete,
            tailing_profile: match &evaluated.preview {
                design::Preview::Tailing(p) => Some(TailingProfileEvidence {
                    profile: p.recipe.profile,
                    status: p.profile_status,
                    findings: p.profile_findings.clone(),
                }),
                _ => None,
            },
            findings: report.findings.clone(),
            assumptions: report.assumptions.clone(),
            exhaustive_alternative_search: Some(report.exhaustive_alternative_search),
        };
        evidence.validate()?;
        Ok(evidence)
    })();
    match result {
        Ok(row) => row,
        Err(error) => error_row(index, row, RowState::Error, error.to_string(), combination),
    }
}
impl RowEvidence {
    pub fn validate(&self) -> Result<(), Error> {
        if self.index >= design::MAX_SELECTED
            || !identifier(&self.construct_id)
            || !label(&self.name)
            || self.input_bindings.len() > design::MAX_INPUTS
            || self.completed_steps > self.planned_steps
            || self.planned_steps > 128
            || self.findings.len() > 4096
            || self.assumptions.len() > 4096
        {
            return Err(bad("invalid row summary dimensions"));
        }
        if self.tailing_profile.as_ref().is_some_and(|profile| {
            profile.findings.len() > 64 || profile.findings.iter().any(|s| s.len() > 4096)
        }) {
            return Err(Error::TooLarge);
        }
        if self
            .error
            .as_ref()
            .is_some_and(|s| s.is_empty() || s.len() > 4096)
            || self
                .findings
                .iter()
                .chain(&self.assumptions)
                .any(|s| s.len() > 65_536)
        {
            return Err(Error::TooLarge);
        }
        if self
            .combination_id
            .as_ref()
            .is_some_and(|s| !s.strip_prefix("combination-").is_some_and(digest))
            || [
                &self.declaration_revision,
                &self.plan_revision,
                &self.target_revision,
                &self.final_product_revision,
            ]
            .into_iter()
            .any(|s| s.as_ref().is_some_and(|s| !digest(s)))
        {
            return Err(bad("invalid row revision identity"));
        }
        let mut ids = BTreeSet::new();
        for input in &self.input_bindings {
            if input.input_id.is_empty()
                || input.input_id.len() > 128
                || !ids.insert(&input.input_id)
            {
                return Err(bad("invalid row source identifiers"));
            }
            crate::cloning::restriction::history::source_role(&input.event)?;
        }
        if let Some(target) = &self.target_binding {
            crate::cloning::restriction::history::source_role(target)?;
        }
        match self.state {
            RowState::Evaluated => {
                if self.status.is_none()
                    || self.error.is_some()
                    || self.combination_id.is_none()
                    || self.declaration_revision.is_none()
                    || self.plan_revision.is_none()
                    || self.target_revision.is_none()
                    || self.exhaustive_alternative_search.is_none()
                    || self.input_bindings.is_empty()
                {
                    return Err(bad("evaluated row lacks actual report identity"));
                }
                if self.final_product_revision.is_some()
                    != (!self.planned_steps.eq(&0) && self.completed_steps == self.planned_steps)
                {
                    return Err(bad("final product is not the actual completed final step"));
                }
                if self.status == Some(Status::ConsistentUnderAssumptions)
                    && (!self.preparation_complete || self.final_product_revision.is_none())
                {
                    return Err(bad(
                        "positive row summary contradicts incomplete reconstruction",
                    ));
                }
            }
            RowState::DraftError | RowState::Error => {
                if self.status.is_some()
                    || self.error.is_none()
                    || self.declaration_revision.is_some()
                    || !self.input_bindings.is_empty()
                    || self.target_binding.is_some()
                    || self.plan_revision.is_some()
                    || self.target_revision.is_some()
                    || self.final_product_revision.is_some()
                    || self.completed_steps != 0
                    || self.planned_steps != 0
                    || self.preparation_complete
                    || self.tailing_profile.is_some()
                    || !self.findings.is_empty()
                    || !self.assumptions.is_empty()
                    || self.exhaustive_alternative_search.is_some()
                {
                    return Err(bad(
                        "uncomputed/error row cannot retain a scientific result",
                    ));
                }
            }
        }
        count(self, MAX_ROW_BYTES)?;
        Ok(())
    }
}
impl Progress {
    pub fn validate(&self) -> Result<(), Error> {
        if !digest(&self.context_revision)
            || self.total_rows == 0
            || self.total_rows > design::MAX_SELECTED
            || self.next_index > self.total_rows
            || self.finished != (self.next_index == self.total_rows)
        {
            return Err(bad("invalid batch progress"));
        }
        if let Some(row) = &self.row {
            row.validate()?;
            if row.index >= self.next_index {
                return Err(bad("progress row is not completed"));
            }
        } else if self.next_index != 0 {
            return Err(bad("completed progress lacks its row summary"));
        }
        Ok(())
    }
}
impl Session {
    pub fn begin(workspace: Workspace, matrix: Request) -> Result<Self, Error> {
        let context_revision = begin_revision(&workspace, &matrix)?;
        Ok(Self {
            workspace,
            request: matrix,
            context_revision,
            rows: vec![],
            work_used: 0,
            work_exhausted: false,
        })
    }
    pub fn summary(&self) -> Progress {
        Progress {
            context_revision: self.context_revision.clone(),
            next_index: self.rows.len(),
            total_rows: self.request.rows.len(),
            finished: self.rows.len() == self.request.rows.len(),
            row: self.rows.last().cloned(),
        }
    }
    pub fn advance(&mut self, context: &str, index: usize) -> Result<Progress, Error> {
        if context != self.context_revision {
            return Err(bad("batch context changed"));
        }
        if index < self.rows.len() {
            let mut progress = self.summary();
            progress.row = Some(self.rows[index].clone());
            return Ok(progress);
        }
        if index != self.rows.len() || index >= self.request.rows.len() {
            return Err(bad("advance exactly the next selected construct"));
        }
        let request = &self.request.rows[index];
        let mut row = match request {
            RowRequest::DraftError { message, .. } => {
                error_row(index, request, RowState::DraftError, message.clone(), None)
            }
            RowRequest::Ready(construct) => {
                let withheld = |index, request| {
                    error_row(
                        index,
                        request,
                        RowState::Error,
                        format!(
                            "Not computed: aggregate batch work budget ({MAX_BATCH_WORK} declared units) was reached. Remaining ready rows require a smaller selected batch."
                        ),
                        None,
                    )
                };
                if self.work_exhausted {
                    withheld(index, request)
                } else {
                    let resolved = self
                        .request
                        .design
                        .as_ref()
                        .ok_or_else(|| bad("typed source design unavailable"))
                        .and_then(|design| {
                            design::resolve_validated(design, construct, &self.workspace)
                        });
                    match resolved {
                        Ok(resolved) => {
                            let charge = resolved.declared_work();
                            if self
                                .work_used
                                .checked_add(charge)
                                .is_none_or(|n| n > MAX_BATCH_WORK)
                            {
                                self.work_exhausted = true;
                                withheld(index, request)
                            } else {
                                self.work_used += charge;
                                evaluate_row(
                                    index,
                                    request,
                                    self.request.design.as_ref(),
                                    Ok(resolved),
                                )
                            }
                        }
                        Err(error) => {
                            evaluate_row(index, request, self.request.design.as_ref(), Err(error))
                        }
                    }
                }
            }
        };
        let remaining = self.request.rows.len() - index - 1;
        let limit = MAX_RECORD_BYTES - 4096 - remaining * 8192;
        if count(&(&self.request, &self.rows, &row), limit).is_err() {
            row = error_row(
                index,
                &self.request.rows[index],
                match &self.request.rows[index] {
                    RowRequest::DraftError { .. } => RowState::DraftError,
                    _ => RowState::Error,
                },
                "Batch summary budget reached; inspect this construct separately.".into(),
                None,
            );
        }
        count(&(&self.request, &self.rows, &row), limit)?;
        self.rows.push(row);
        let progress = self.summary();
        progress.validate()?;
        Ok(progress)
    }
    pub fn finish_record(&self) -> Result<Record, Error> {
        if self.rows.len() != self.request.rows.len() {
            return Err(bad("finish every selected row before capturing a batch"));
        }
        let context_revision = record_revision(ENGINE_VERSION, &self.request, &self.rows)?;
        let record = Record {
            schema_version: 3,
            context_revision,
            engine_version: ENGINE_VERSION.into(),
            request: self.request.clone(),
            rows: self.rows.clone(),
        };
        record.validate()?;
        Ok(record)
    }
    pub fn capture_to(&self, current: &Workspace, context: &str) -> Result<Commit, Error> {
        if context != self.context_revision {
            return Err(bad("batch context changed"));
        }
        let record = self.finish_record()?;
        design::workspace_guard(current)?;
        source_guard(current, &record)?;
        let fingerprint = run_lineage::workspace_fingerprint(current)?;
        install(current.clone(), record, context.to_string(), fingerprint)
    }
    pub fn into_capture(self, context: &str) -> Result<Commit, Error> {
        if context != self.context_revision {
            return Err(bad("batch context changed"));
        }
        let record = self.finish_record()?;
        source_guard(&self.workspace, &record)?;
        let fingerprint = run_lineage::workspace_fingerprint(&self.workspace)?;
        install(self.workspace, record, context.to_string(), fingerprint)
    }
}
fn record_revision(engine: &str, request: &Request, rows: &[RowEvidence]) -> Result<String, Error> {
    count(&(engine, request, rows), MAX_RECORD_BYTES - 1024)?;
    crate::encoding::bounded_json::hash(&("combinatorial-batch-record-3", engine, request, rows))
}

impl Record {
    /// Explicit indivisible documentary publication must attest every embedded
    /// closed source, including unselected choices and diagnostic/error rows.
    /// No historical verdict or literal edit is recomputed or upgraded here.
    pub fn validate_source_closure(&self, workspace: &Workspace) -> Result<(), Error> {
        self.validate()?;
        let Some(design) = &self.request.design else {
            return Ok(());
        };
        design::source_preflight(design, workspace)?;
        for row in &self.request.rows {
            let RowRequest::Ready(row) = row else {
                continue;
            };
            if let Some(overrides) = &row.overrides {
                for p in &overrides.preparations {
                    let Some(preparation) = &p.preparation else {
                        continue;
                    };
                    if !matches!(
                        preparation,
                        crate::cloning::assembly::recipe::Preparation::FromRun { .. }
                            | crate::cloning::assembly::recipe::Preparation::FromDerived { .. }
                    ) {
                        continue;
                    }
                    let bin = design
                        .bins
                        .iter()
                        .find(|b| b.id == p.bin_id)
                        .ok_or_else(|| bad("embedded saved source override has no exact bin"))?;
                    let selected = row
                        .choices
                        .iter()
                        .find(|c| c.bin_id == bin.id)
                        .ok_or_else(|| bad("embedded saved source override lacks exact choice"))?;
                    let choice = bin
                        .choices
                        .iter()
                        .find(|c| c.id == selected.choice_id)
                        .ok_or_else(|| bad("embedded saved source override choice absent"))?;
                    design::preparation_binding(design, &choice.source, preparation)?;
                }
            }
        }
        Ok(())
    }
    pub fn validate(&self) -> Result<(), Error> {
        if self.schema_version != 3
            || !digest(&self.context_revision)
            || self.engine_version.is_empty()
            || self.engine_version.len() > 256
            || self.rows.len() != self.request.rows.len()
        {
            return Err(bad("invalid historical batch record"));
        }
        self.request.validate()?;
        for (index, (request, row)) in self.request.rows.iter().zip(&self.rows).enumerate() {
            row.validate()?;
            let (id, name) = row_identity(request);
            if row.index != index || row.construct_id != id || row.name != name {
                return Err(bad("recorded row identity differs from its request"));
            }
            match (request, row.state) {
                (RowRequest::DraftError { .. }, RowState::DraftError) => {}
                (RowRequest::Ready(construct), RowState::Evaluated) => {
                    let design = self
                        .request
                        .design
                        .as_ref()
                        .ok_or_else(|| bad("missing recorded design"))?;
                    if row.combination_id.as_ref()
                        != Some(&design::combination_key(design, &construct.choices)?)
                    {
                        return Err(bad("recorded combination differs from selected choices"));
                    }
                    let method = construct
                        .overrides
                        .as_ref()
                        .and_then(|o| o.method.as_ref())
                        .unwrap_or(&design.method);
                    match (method, &row.tailing_profile) {
                        (design::Method::Tailing { profile, .. }, Some(evidence))
                            if *profile == evidence.profile => {}
                        (design::Method::Tailing { .. }, _) | (_, Some(_)) => {
                            return Err(bad(
                                "recorded tailing profile differs from the selected method",
                            ));
                        }
                        (_, None) => {}
                    }
                    let expected = design.selected_input_bindings(construct)?;
                    if row.input_bindings.len() != expected.len()
                        || row
                            .input_bindings
                            .iter()
                            .zip(&expected)
                            .any(|(a, b)| a.input_id != b.input_id || a.event != b.event)
                    {
                        return Err(bad(
                            "recorded source events differ from the selected alternatives",
                        ));
                    }
                    match &construct.target {
                        design::TargetDeclaration::Existing(target) => {
                            if row.target_binding.as_ref() != Some(target)
                                || row.target_revision.as_ref() != Some(&target.revision)
                            {
                                return Err(bad(
                                    "recorded target differs from the explicitly selected event",
                                ));
                            }
                        }
                        design::TargetDeclaration::OriginalRegions { .. } => {
                            if row.target_binding.is_some() {
                                return Err(bad(
                                    "literal target has an invented existing target event",
                                ));
                            }
                        }
                    }
                }
                (RowRequest::Ready(_), RowState::Error) => {}
                _ => return Err(bad("recorded row state contradicts its request")),
            }
        }
        if self.context_revision
            != record_revision(&self.engine_version, &self.request, &self.rows)?
        {
            return Err(bad("historical batch content revision mismatch"));
        }
        count(self, MAX_RECORD_BYTES)?;
        Ok(())
    }
}
fn source_guard(workspace: &Workspace, record: &Record) -> Result<(), Error> {
    attachment_versions(workspace)?;
    record.validate_source_closure(workspace)?;
    for row in &record.rows {
        for input in &row.input_bindings {
            design::resolve_event(workspace, &input.event)?;
        }
        if let Some(target) = &row.target_binding {
            design::resolve_event(workspace, target)?;
        }
    }
    Ok(())
}
struct Bytes(Vec<u8>);
impl io::Write for Bytes {
    fn write(&mut self, bytes: &[u8]) -> io::Result<usize> {
        if self.0.len().saturating_add(bytes.len()) > MAX_RECORD_BYTES {
            return Err(io::Error::other("batch record limit"));
        }
        self.0.extend_from_slice(bytes);
        Ok(bytes.len())
    }
    fn flush(&mut self) -> io::Result<()> {
        Ok(())
    }
}
fn install(
    mut workspace: Workspace,
    record: Record,
    context_revision: String,
    source_workspace_fingerprint: String,
) -> Result<Commit, Error> {
    let mut encoded = Bytes(vec![]);
    serde_json::to_writer(&mut encoded, &record).map_err(|_| Error::TooLarge)?;
    let attachment_hash = blake3::hash(&encoded.0).to_hex().to_string();
    let name = record
        .request
        .design
        .as_ref()
        .map(|d| d.name.as_str())
        .unwrap_or("Unprepared combinatorial batch");
    if let Some(existing) = workspace.attachments.get(&attachment_hash) {
        if existing.kind != KIND
            || existing.name != name
            || existing.media_type != "application/json"
            || existing.bytes != encoded.0
        {
            return Err(bad(
                "identical record bytes have inconsistent attachment metadata",
            ));
        }
    } else {
        if workspace.attachments.len() >= workspace::MAX_ATTACHMENTS
            || workspace
                .attachments
                .values()
                .filter(|a| a.kind == KIND)
                .count()
                >= MAX_RECORDS
        {
            return Err(Error::TooLarge);
        }
        workspace.attachments.insert(
            attachment_hash.clone(),
            WorkspaceAttachment {
                kind: KIND.into(),
                name: name.into(),
                media_type: "application/json".into(),
                bytes: encoded.0,
            },
        );
    }
    workspace.validate()?;
    let capture_revision = capture_token(&source_workspace_fingerprint, &context_revision)?;
    Ok(Commit {
        workspace,
        context_revision,
        capture_revision,
        source_workspace_fingerprint,
        attachment_hash,
    })
}
pub fn load(workspace: &Workspace, hash: &str) -> Result<Record, Error> {
    attachment_versions(workspace)?;
    design::workspace_guard(workspace)?;
    load_validated(workspace, hash)
}
fn load_validated(workspace: &Workspace, hash: &str) -> Result<Record, Error> {
    if !digest(hash) {
        return Err(bad("invalid batch attachment hash"));
    }
    let asset = workspace
        .attachments
        .get(hash)
        .ok_or_else(|| bad("saved batch is unavailable"))?;
    if asset.kind != KIND
        || asset.media_type != "application/json"
        || asset.bytes.len() > MAX_RECORD_BYTES
        || blake3::hash(&asset.bytes).to_hex().as_str() != hash
    {
        return Err(bad("invalid batch attachment type, hash or size"));
    }
    let record: Record = serde_json::from_slice(&asset.bytes).map_err(|e| bad(&e.to_string()))?;
    source_guard(workspace, &record)?;
    Ok(record)
}
pub fn list(workspace: &Workspace) -> Result<Vec<Summary>, Error> {
    attachment_versions(workspace)?;
    design::workspace_guard(workspace)?;
    if workspace
        .attachments
        .values()
        .filter(|a| a.kind == KIND)
        .count()
        > MAX_RECORDS
    {
        return Err(Error::TooLarge);
    }
    workspace
        .attachments
        .iter()
        .filter(|(_, asset)| asset.kind == KIND)
        .map(|(hash, _)| {
            let record = load_validated(workspace, hash)?;
            let evaluated_rows = record
                .rows
                .iter()
                .filter(|r| r.state == RowState::Evaluated)
                .count();
            Ok(Summary {
                attachment_hash: hash.clone(),
                context_revision: record.context_revision,
                design_name: record
                    .request
                    .design
                    .as_ref()
                    .map(|d| d.name.clone())
                    .unwrap_or_else(|| "Unprepared combinatorial batch".into()),
                total_rows: record.rows.len(),
                evaluated_rows,
                error_rows: record.rows.len() - evaluated_rows,
                recorded_engine: record.engine_version,
            })
        })
        .collect()
}
pub fn replay(workspace: &Workspace, hash: &str, index: usize) -> Result<Replayed, Error> {
    let record = load(workspace, hash)?;
    let row = record
        .request
        .rows
        .get(index)
        .ok_or_else(|| bad("recorded row is unavailable"))?;
    let resolved = match (row, record.request.design.as_ref()) {
        (RowRequest::Ready(construct), Some(design)) => {
            design::resolve_validated(design, construct, workspace)
        }
        _ => Err(bad("no computable declaration")),
    };
    let current = evaluate_row(index, row, record.request.design.as_ref(), resolved);
    let replayed = Replayed {
        record_hash: hash.into(),
        row_index: index,
        recorded: record.rows[index].clone(),
        current,
        current_engine: ENGINE_VERSION.into(),
    };
    replayed.validate()?;
    Ok(replayed)
}
impl Replayed {
    pub fn validate(&self) -> Result<(), Error> {
        if !digest(&self.record_hash)
            || self.current_engine != ENGINE_VERSION
            || self.row_index != self.recorded.index
            || self.row_index != self.current.index
            || self.recorded.construct_id != self.current.construct_id
        {
            return Err(bad("invalid replay identity"));
        }
        self.recorded.validate()?;
        self.current.validate()
    }
}
impl Commit {
    pub fn validate(&self) -> Result<(), Error> {
        if !digest(&self.context_revision)
            || !digest(&self.capture_revision)
            || !digest(&self.source_workspace_fingerprint)
            || !digest(&self.attachment_hash)
            || self.capture_revision
                != capture_token(&self.source_workspace_fingerprint, &self.context_revision)?
        {
            return Err(bad("invalid batch installation receipt"));
        }
        load(&self.workspace, &self.attachment_hash)?;
        Ok(())
    }
    pub fn validate_binding(&self, revision: &str) -> Result<(), Error> {
        if self.capture_revision != revision {
            return Err(bad("batch installation context changed"));
        }
        self.validate()
    }
}
