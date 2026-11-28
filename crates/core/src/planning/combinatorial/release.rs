//! Atomic selected-row release over existing native finalization.
//! Historical batch summaries never authorize outputs. Each row is compiled afresh
//! and cloning_outputs::finalize performs the only scientific verifier pass here.
use crate::{
    cloning::assembly::recipe as assembly_recipe,
    documents::workspace::{
        self, RevisionEvent, Workspace, WorkspaceAttachment, WorkspaceDocument,
    },
    material::{
        lineage::{self as run_lineage, EventBinding, InputEventBinding},
        outputs as cloning_outputs,
    },
    planning::combinatorial::{self as design, batch},
};
use collagene_model::{Error, Molecule};
use serde::{Deserialize, Serialize};
use std::{collections::BTreeSet, io};
pub const KIND: &str = "combinatorial-release-v1";
pub const CONTRACT: &str = "combinatorial-release-1";
pub const MAX_MANIFEST_BYTES: usize = 32 * 1024 * 1024;
pub const MAX_INSPECTION_BYTES: usize = 32 * 1024 * 1024;
const MAX_ROW_BYTES: usize = 256 * 1024;
macro_rules! dto {($($t:tt)*)=>{#[derive(Clone,Debug,Serialize,Deserialize)]#[cfg_attr(feature="json-schema",derive(schemars::JsonSchema))]#[serde(deny_unknown_fields)]$($t)*};}
dto! {pub struct Selection{pub construct_id:String,pub options:cloning_outputs::Options}}
dto! {pub struct Request{pub batch_hash:String,#[serde(deserialize_with="decode_selections")]pub rows:Vec<Selection>}}
dto! {pub struct ReleasedRow{
    pub ordinal:usize,pub construct_id:String,pub combination_id:String,pub declaration_revision:String,
    #[serde(deserialize_with="decode_inputs")]pub inputs:Vec<InputEventBinding>,pub target:EventBinding,pub options:cloning_outputs::Options,
    pub output:cloning_outputs::Receipt,
}}
dto! {pub struct Progress{pub context_revision:String,pub next_index:usize,pub total_rows:usize,pub finished:bool,pub failed:Option<String>,pub row:Option<ReleasedRow>}}
dto! {pub struct Manifest{
    pub schema_version:u16,pub context_revision:String,pub source_workspace_fingerprint:String,
    pub batch_hash:String,pub batch_context:String,pub engine_version:String,
    pub request:Request,#[serde(deserialize_with="decode_rows")]pub rows:Vec<ReleasedRow>,
}}
dto! {pub struct Commit{
    #[cfg_attr(feature="json-schema",schemars(with="serde_json::Value"))]pub workspace:Workspace,
    pub context_revision:String,pub source_workspace_fingerprint:String,pub manifest_hash:String,
}}
dto! {pub struct Inspection{pub context_revision:String,pub index:usize,pub row:ReleasedRow,pub product:Molecule,pub target:Molecule}}
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
pub struct Summary {
    pub manifest_hash: String,
    pub context_revision: String,
    pub batch_hash: String,
    pub row_count: usize,
    pub recorded_engine: String,
}
fn bounded<'de, D, T, const N: usize>(deserializer: D) -> Result<Vec<T>, D::Error>
where
    D: serde::Deserializer<'de>,
    T: Deserialize<'de>,
{
    struct Visitor<T, const N: usize>(std::marker::PhantomData<T>);
    impl<'de, T: Deserialize<'de>, const N: usize> serde::de::Visitor<'de> for Visitor<T, N> {
        type Value = Vec<T>;
        fn expecting(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
            write!(f, "at most {N} release values")
        }
        fn visit_seq<A: serde::de::SeqAccess<'de>>(self, mut seq: A) -> Result<Vec<T>, A::Error> {
            if seq.size_hint().is_some_and(|n| n > N) {
                return Err(serde::de::Error::custom("release value limit"));
            }
            let mut values = vec![];
            while values.len() < N {
                match seq.next_element()? {
                    Some(value) => values.push(value),
                    None => return Ok(values),
                }
            }
            if seq.next_element::<serde::de::IgnoredAny>()?.is_some() {
                return Err(serde::de::Error::custom("release value limit"));
            }
            Ok(values)
        }
    }
    deserializer.deserialize_seq(Visitor::<T, N>(std::marker::PhantomData))
}
fn decode_selections<'de, D: serde::Deserializer<'de>>(d: D) -> Result<Vec<Selection>, D::Error> {
    bounded::<D, Selection, 256>(d)
}
fn decode_rows<'de, D: serde::Deserializer<'de>>(d: D) -> Result<Vec<ReleasedRow>, D::Error> {
    bounded::<D, ReleasedRow, 256>(d)
}
fn decode_inputs<'de, D: serde::Deserializer<'de>>(
    d: D,
) -> Result<Vec<InputEventBinding>, D::Error> {
    bounded::<D, InputEventBinding, 12>(d)
}
/// Private, non-serializable all-or-nothing state. A failed row destroys the
/// private working frame; no caller can extract a successful prefix workspace.
pub struct ReleaseSession {
    working: Option<Workspace>,
    record: batch::Record,
    request: Request,
    source_workspace_fingerprint: String,
    context_revision: String,
    rows: Vec<ReleasedRow>,
    failure: Option<String>,
    work_used: usize,
}
fn bad(message: &str) -> Error {
    Error::Format(format!("Combinatorial release: {message}"))
}
fn digest(value: &str) -> bool {
    value.len() == 64
        && value
            .bytes()
            .all(|b| b.is_ascii_digit() || (b'a'..=b'f').contains(&b))
}
fn identifier(value: &str, max: usize) -> bool {
    !value.is_empty()
        && value.len() <= max
        && value
            .bytes()
            .all(|b| b.is_ascii_alphanumeric() || matches!(b, b'_' | b'-' | b'.' | b':'))
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
            .ok_or_else(|| io::Error::other("bulk release byte limit"))?;
        Ok(bytes.len())
    }
    fn flush(&mut self) -> io::Result<()> {
        Ok(())
    }
}
fn count(value: &impl Serialize, limit: usize) -> Result<(), Error> {
    serde_json::to_writer(&mut Counter { bytes: 0, limit }, value).map_err(|_| Error::TooLarge)
}
impl Request {
    pub fn validate(&self) -> Result<(), Error> {
        if !digest(&self.batch_hash)
            || self.rows.is_empty()
            || self.rows.len() > design::MAX_SELECTED
        {
            return Err(bad("select a saved batch and 1–256 explicit rows"));
        }
        let mut ids = BTreeSet::new();
        for row in &self.rows {
            if row.construct_id.is_empty()
                || row.construct_id.len() > 48
                || !row
                    .construct_id
                    .bytes()
                    .all(|b| b.is_ascii_alphanumeric() || matches!(b, b'_' | b'-' | b'.'))
                || !ids.insert(&row.construct_id)
            {
                return Err(bad("invalid or duplicate selected construct ID"));
            }
            row.options.validate()?;
        }
        count(self, 1024 * 1024)
    }
}
fn selected<'a>(
    record: &'a batch::Record,
    id: &str,
) -> Result<&'a design::SelectedConstruct, Error> {
    record
        .request
        .rows
        .iter()
        .find_map(|row| match row {
            batch::RowRequest::Ready(row) if row.id == id => Some(row),
            _ => None,
        })
        .ok_or_else(|| bad("release requires an explicitly selected ready-row declaration"))
}
fn attest_original_sources(
    workspace: &Workspace,
    record: &batch::Record,
    request: &Request,
) -> Result<(), Error> {
    let design = record
        .request
        .design
        .as_ref()
        .ok_or_else(|| bad("saved batch has no typed source design"))?;
    for selection in &request.rows {
        let row = selected(record, &selection.construct_id)?;
        design::validate_construct(design, row)?;
        for bin in &design.bins {
            let chosen = row
                .choices
                .iter()
                .find(|choice| choice.bin_id == bin.id)
                .ok_or_else(|| bad("missing selected source bin"))?;
            let choice = bin
                .choices
                .iter()
                .find(|choice| choice.id == chosen.choice_id)
                .ok_or_else(|| bad("missing selected source alternative"))?;
            design::resolve_event(workspace, &choice.source)?;
        }
        if let design::TargetDeclaration::Existing(target) = &row.target {
            design::resolve_event(workspace, target)?;
        }
    }
    Ok(())
}
pub fn request_revision(request: &Request) -> Result<String, Error> {
    request.validate()?;
    assembly_recipe::hash(&(CONTRACT, "request", request))
}
pub fn preflight(workspace: &Workspace, request: &Request) -> Result<(), Error> {
    request.validate()?;
    let record = batch::load(workspace, &request.batch_hash)?;
    // Freeze membership against the original frame. A later row must not acquire
    // a guessed source/target event from an earlier private generated product.
    attest_original_sources(workspace, &record, request)
}
pub fn revision(workspace: &Workspace, request: &Request) -> Result<String, Error> {
    preflight(workspace, request)?;
    assembly_recipe::hash(&(
        CONTRACT,
        crate::verification::ENGINE_VERSION,
        cloning_outputs::ENGINE_VERSION,
        run_lineage::workspace_fingerprint(workspace)?,
        request,
    ))
}
fn target_event(
    workspace: &mut Workspace,
    resolved: &design::Resolved,
) -> Result<EventBinding, Error> {
    if let Some(target) = &resolved.target_binding {
        design::resolve_event(workspace, target)?;
        return Ok(target.clone());
    }
    if resolved.intent.is_none() {
        return Err(bad("missing independent target declaration"));
    }
    let expected_revision = resolved.target.revision();
    let parents = resolved
        .bindings
        .iter()
        .map(|b| b.event.event_id.clone())
        .collect::<BTreeSet<_>>()
        .into_iter()
        .collect::<Vec<_>>();
    for input in &resolved.bindings {
        design::resolve_event(workspace, &input.event)?;
    }
    let operation = format!(
        "Declared combinatorial target for {}; original source expectation only, not a physical product or experimental confirmation",
        resolved.construct_id
    );
    // Historical exact targets may be reused without rewinding their edited heads.
    if let Some(event) = workspace.history.iter().find(|event| {
        event.revision == expected_revision
            && event.parents == parents
            && event.operation == operation
    }) {
        return Ok(EventBinding {
            document_id: event.document_id.clone(),
            event_id: event.id.clone(),
            revision: event.revision.clone(),
        });
    }
    if workspace.documents.len() >= workspace::MAX_DOCUMENTS
        || workspace.history.len() >= workspace::MAX_HISTORY_EVENTS
    {
        return Err(Error::TooLarge);
    }
    let key = assembly_recipe::hash(&(
        CONTRACT,
        "intent",
        &resolved.construct_id,
        &resolved.declaration_revision,
        &parents,
    ))?;
    let base = format!("bulk-intent-{key}");
    let mut id = base.clone();
    let mut suffix = 0usize;
    while workspace.documents.iter().any(|d| d.id == id) {
        suffix = suffix.checked_add(1).ok_or(Error::TooLarge)?;
        if suffix > workspace::MAX_DOCUMENTS {
            return Err(Error::TooLarge);
        }
        id = format!("{base}-{suffix}");
    }
    let event = RevisionEvent::new(id.clone(), expected_revision.clone(), parents, operation)?;
    let target = EventBinding {
        document_id: id.clone(),
        event_id: event.id.clone(),
        revision: expected_revision.clone(),
    };
    workspace.documents.push(WorkspaceDocument {
        id,
        head: expected_revision.clone(),
        head_event: Some(event.id.clone()),
        undo: vec![],
        redo: vec![],
        undo_events: vec![],
        redo_events: vec![],
    });
    workspace
        .revisions
        .insert(expected_revision, resolved.target.clone());
    workspace.history.push(event);
    workspace.validate()?;
    Ok(target)
}
impl ReleaseSession {
    pub fn begin(workspace: Workspace, request: Request) -> Result<Self, Error> {
        let context_revision = revision(&workspace, &request)?;
        let record = batch::load(&workspace, &request.batch_hash)?;
        let source_workspace_fingerprint = run_lineage::workspace_fingerprint(&workspace)?;
        Ok(Self {
            working: Some(workspace),
            record,
            request,
            source_workspace_fingerprint,
            context_revision,
            rows: vec![],
            failure: None,
            work_used: 0,
        })
    }
    pub fn summary(&self) -> Progress {
        Progress {
            context_revision: self.context_revision.clone(),
            next_index: self.rows.len(),
            total_rows: self.request.rows.len(),
            finished: self.failure.is_none() && self.rows.len() == self.request.rows.len(),
            failed: self.failure.clone(),
            row: self.rows.last().cloned(),
        }
    }
    pub fn advance(&mut self, context: &str, index: usize) -> Result<Progress, Error> {
        if context != self.context_revision {
            return Err(bad("release context changed"));
        }
        if self.failure.is_some() {
            return Err(bad("release failed; begin a new private release"));
        }
        if index < self.rows.len() {
            let mut progress = self.summary();
            progress.row = Some(self.rows[index].clone());
            return Ok(progress);
        }
        if index != self.rows.len() || index >= self.request.rows.len() {
            return Err(bad("release exactly the next explicitly selected row"));
        }
        let outcome = (|| -> Result<ReleasedRow, Error> {
            let selection = &self.request.rows[index];
            let construct = selected(&self.record, &selection.construct_id)?;
            let design = self
                .record
                .request
                .design
                .as_ref()
                .ok_or_else(|| bad("missing typed source design"))?;
            let working = self
                .working
                .as_mut()
                .ok_or_else(|| bad("private release state is unavailable"))?;
            let resolved = design::resolve_validated(design, construct, working)?;
            let charge = resolved.declared_work();
            self.work_used = self
                .work_used
                .checked_add(charge)
                .filter(|n| *n <= batch::MAX_BATCH_WORK)
                .ok_or_else(|| bad("aggregate release work limit reached; no outputs published"))?;
            // Compile + native finalization. Do NOT call evaluate here: finalization
            // freshly verifies once and records exact operation products itself.
            let plan = design::compile_resolved(&resolved)?;
            let compiled = design::Compiled { resolved, plan };
            let target = target_event(working, &compiled.resolved)?;
            if target.revision != compiled.plan.target_revision {
                return Err(bad("independent target binding changed during release"));
            }
            // Every preparation/product/target remains retained by the established
            // output layer. Capacity failure must roll back the complete release.
            if working
                .documents
                .len()
                .checked_add(compiled.plan.steps.len())
                .is_none_or(|n| n > workspace::MAX_DOCUMENTS)
            {
                return Err(Error::TooLarge);
            }
            let output = cloning_outputs::finalize(
                working,
                &compiled.plan,
                target.clone(),
                compiled.resolved.bindings.clone(),
                &selection.options,
            )?;
            let row = ReleasedRow {
                ordinal: index,
                construct_id: construct.id.clone(),
                combination_id: compiled.resolved.combination_id,
                declaration_revision: compiled.resolved.declaration_revision,
                inputs: compiled.resolved.bindings,
                target,
                options: selection.options.clone(),
                output: output.receipt,
            };
            row.validate()?;
            count(
                &(&self.request, &self.rows, &row),
                MAX_MANIFEST_BYTES - 4096,
            )?;
            *working = output.workspace;
            Ok(row)
        })();
        match outcome {
            Ok(row) => {
                self.rows.push(row);
                Ok(self.summary())
            }
            Err(error) => {
                let full = error.to_string();
                let mut message = full.chars().take(1000).collect::<String>();
                if full.chars().count() > 1000 {
                    message
                        .push_str("… (abbreviated release failure; original sources are retained)");
                }
                self.failure = Some(message);
                self.working = None;
                Err(error)
            }
        }
    }
    /// Read-only inspection of an actually retained private result. It does not
    /// rerun science or create a successful report from stored molecular text.
    pub fn inspect(&self, index: usize) -> Result<Inspection, Error> {
        if self.failure.is_some() {
            return Err(bad("failed release has discarded every private result"));
        }
        let row = self
            .rows
            .get(index)
            .ok_or_else(|| bad("this selected row has not been released privately"))?;
        let workspace = self
            .working
            .as_ref()
            .ok_or_else(|| bad("private release state is unavailable"))?;
        let product = design::resolve_event(workspace, &row.output.construct)?;
        let target = design::resolve_event(workspace, &row.target)?;
        product.validate()?;
        target.validate()?;
        count(&(product, target, row), MAX_INSPECTION_BYTES - 4096)?;
        let inspection = Inspection {
            context_revision: self.context_revision.clone(),
            index,
            row: row.clone(),
            product: product.clone(),
            target: target.clone(),
        };
        inspection.validate()?;
        Ok(inspection)
    }
    /// Browser integration verifies the actual live frame before one IndexedDB
    /// transaction; a changed original frame never receives a partial rebase.
    pub fn finish_to(self, current: &Workspace, context: &str) -> Result<Commit, Error> {
        if run_lineage::workspace_fingerprint(current)? != self.source_workspace_fingerprint {
            return Err(bad("source workspace changed; release was not installed"));
        }
        self.into_commit(context)
    }
    /// Immutable WorkspaceRef/MCP path. Host still verifies the returned source
    /// fingerprint against the explicitly imported original workspace reference.
    pub fn into_commit(mut self, context: &str) -> Result<Commit, Error> {
        if context != self.context_revision
            || self.failure.is_some()
            || self.rows.len() != self.request.rows.len()
        {
            return Err(bad(
                "release must complete every selected row before publication",
            ));
        }
        let mut workspace = self
            .working
            .take()
            .ok_or_else(|| bad("private release state was discarded"))?;
        let manifest = Manifest {
            schema_version: 1,
            context_revision: self.context_revision.clone(),
            source_workspace_fingerprint: self.source_workspace_fingerprint.clone(),
            batch_hash: self.request.batch_hash.clone(),
            batch_context: self.record.context_revision,
            engine_version: crate::verification::ENGINE_VERSION.into(),
            request: self.request,
            rows: self.rows,
        };
        count(&manifest, MAX_MANIFEST_BYTES)?;
        let bytes = serde_json::to_vec(&manifest).map_err(|e| bad(&e.to_string()))?;
        let manifest_hash = blake3::hash(&bytes).to_hex().to_string();
        if let Some(existing) = workspace.attachments.get(&manifest_hash) {
            if existing.kind != KIND
                || existing.media_type != "application/json"
                || existing.bytes != bytes
            {
                return Err(bad("manifest attachment identity conflict"));
            }
        } else {
            if workspace.attachments.len() >= workspace::MAX_ATTACHMENTS {
                return Err(Error::TooLarge);
            }
            workspace.attachments.insert(
                manifest_hash.clone(),
                WorkspaceAttachment {
                    kind: KIND.into(),
                    name: "Combinatorial output release".into(),
                    media_type: "application/json".into(),
                    bytes,
                },
            );
        }
        workspace.validate()?;
        let commit = Commit {
            workspace,
            context_revision: self.context_revision,
            source_workspace_fingerprint: self.source_workspace_fingerprint,
            manifest_hash,
        };
        commit.validate()?;
        Ok(commit)
    }
}
impl ReleasedRow {
    pub fn validate(&self) -> Result<(), Error> {
        self.options.validate()?;
        if self.ordinal >= design::MAX_SELECTED
            || !identifier(&self.construct_id, 48)
            || !self
                .combination_id
                .strip_prefix("combination-")
                .is_some_and(digest)
            || !digest(&self.declaration_revision)
            || self.inputs.is_empty()
            || self.inputs.len() > design::MAX_BINS
        {
            return Err(bad("invalid released row identity or dimensions"));
        }
        let mut ids = BTreeSet::new();
        for input in &self.inputs {
            if !input.input_id.strip_prefix("source-").is_some_and(digest)
                || !ids.insert(&input.input_id)
            {
                return Err(bad(
                    "released row source identities are invalid or duplicate",
                ));
            }
            crate::cloning::restriction::history::source_role(&input.event)?;
        }
        crate::cloning::restriction::history::source_role(&self.target)?;
        crate::cloning::restriction::history::source_role(&self.output.construct)?;
        let output = &self.output;
        let run = &output.run;
        if !digest(&output.context_revision)
            || !digest(&output.source_workspace_fingerprint)
            || output.source_workspace_fingerprint != run.source_workspace_fingerprint
            || !digest(&run.run_id)
            || !digest(&run.attachment_hash)
            || run.products.is_empty()
            || run.products.len() > 128
            || run.recorded_engine.is_empty()
            || run.recorded_engine.len() > 128
            || !run.recorded_engine.bytes().all(|b| b.is_ascii_graphic())
            || run.recorded_status == crate::verification::evidence::Status::Invalid
            || (run.recorded_status == crate::verification::evidence::Status::Incomplete
                && !self.options.allow_incomplete)
        {
            return Err(bad("invalid native output/run receipt"));
        }
        let mut products = BTreeSet::new();
        for (index, product) in run.products.iter().enumerate() {
            if product.step_index != index || !products.insert(product.document_id.clone()) {
                return Err(bad("released operation product identity is invalid"));
            }
            crate::cloning::restriction::history::source_role(&EventBinding {
                document_id: product.document_id.clone(),
                event_id: product.event_id.clone(),
                revision: product.revision.clone(),
            })?;
        }
        if output.visible_products.len() > 129
            || output.hidden_products.len() > 128
            || self.options.save_primers != output.primer_attachment.is_some()
            || output
                .primer_attachment
                .as_ref()
                .is_some_and(|hash| !digest(hash))
        {
            return Err(bad("invalid released output options or counts"));
        }
        let mut seen = BTreeSet::new();
        for id in output
            .visible_products
            .iter()
            .chain(&output.hidden_products)
        {
            if !identifier(id, 128)
                || !seen.insert(id)
                || (!products.contains(id) && *id != output.construct.document_id)
            {
                return Err(bad(
                    "released document visibility differs from its products",
                ));
            }
        }
        let extra = usize::from(!products.contains(&output.construct.document_id));
        if seen.len() != run.products.len() + extra
            || !output
                .visible_products
                .contains(&output.construct.document_id)
        {
            return Err(bad("released output presentation is incomplete"));
        }
        count(self, MAX_ROW_BYTES)
    }
}
impl Progress {
    pub fn validate(&self) -> Result<(), Error> {
        if !digest(&self.context_revision)
            || self.total_rows == 0
            || self.total_rows > design::MAX_SELECTED
            || self.next_index > self.total_rows
            || self.finished != (self.failed.is_none() && self.next_index == self.total_rows)
            || self
                .failed
                .as_ref()
                .is_some_and(|s| s.is_empty() || s.len() > 4096)
        {
            return Err(bad("invalid private release progress"));
        }
        if let Some(row) = &self.row {
            row.validate()?;
            if row.ordinal >= self.next_index {
                return Err(bad("private released row is not completed"));
            }
        } else if self.next_index != 0 {
            return Err(bad("private progress lacks its completed row"));
        }
        Ok(())
    }
}
impl Inspection {
    pub fn validate(&self) -> Result<(), Error> {
        if !digest(&self.context_revision) || self.index != self.row.ordinal {
            return Err(bad("private inspection context differs from released row"));
        }
        self.row.validate()?;
        self.product.validate()?;
        self.target.validate()?;
        if self.product.sequence.alphabet() != collagene_model::Alphabet::Dna
            || self.target.sequence.alphabet() != collagene_model::Alphabet::Dna
            || self.product.revision() != self.row.output.construct.revision
            || self.target.revision() != self.row.target.revision
            || self.product.name != self.row.options.name
        {
            return Err(bad(
                "private product/target models differ from their exact retained events",
            ));
        }
        count(self, MAX_INSPECTION_BYTES)
    }
}
impl Manifest {
    /// Reconstruct the current literal independent target and declared recipe,
    /// reusing fresh native source views without executing the new requested
    /// preparations. Returns their actual declared costs from that same context.
    /// Fresh export/release reuse must use this
    /// before verification: verifying a forged plan against a forged target alone
    /// would otherwise turn an imported self-target into apparent agreement.
    pub fn validate_current_declarations(
        &self,
        workspace: &Workspace,
    ) -> Result<std::collections::BTreeMap<String, usize>, Error> {
        self.validate(workspace)?;
        let batch = batch::load(workspace, &self.batch_hash)?;
        let source_design = batch
            .request
            .design
            .as_ref()
            .ok_or_else(|| bad("missing source design"))?;
        let mut costs = std::collections::BTreeMap::new();
        for row in &self.rows {
            let declared = selected(&batch, &row.construct_id)?;
            let resolved = design::resolve_validated(source_design, declared, workspace)?;
            if resolved.declaration_revision != row.declaration_revision
                || resolved.target.revision() != row.target.revision
            {
                return Err(bad(
                    "released target differs from its independently declared original-source intent",
                ));
            }
            if let Some(target) = &resolved.target_binding {
                if target != &row.target {
                    return Err(bad(
                        "released target differs from its explicit original event",
                    ));
                }
            }
            costs.insert(row.construct_id.clone(), resolved.declared_work());
            let plan = design::compile_resolved(&resolved)?;
            let asset = workspace
                .attachments
                .get(&row.output.run.attachment_hash)
                .ok_or_else(|| bad("released run unavailable"))?;
            let run = run_lineage::load_run_record(&asset.bytes)?;
            if run.plan_revision != crate::planning::plan::DesignPlan::revision(&plan) {
                return Err(bad(
                    "released plan differs from the explicitly declared per-construct preparation/method",
                ));
            }
        }
        Ok(costs)
    }
    /// Read-only linkage validation. Imported manifests never restore a live
    /// ReleaseSession or grant product/experimental authority.
    pub fn validate(&self, workspace: &Workspace) -> Result<(), Error> {
        self.request.validate()?;
        if self.schema_version != 1
            || !digest(&self.context_revision)
            || !digest(&self.source_workspace_fingerprint)
            || self.batch_hash != self.request.batch_hash
            || self.rows.len() != self.request.rows.len()
            || self.engine_version.is_empty()
            || self.engine_version.len() > 256
        {
            return Err(bad("invalid historical release manifest"));
        }
        let batch = batch::load(workspace, &self.batch_hash)?;
        if batch.context_revision != self.batch_context {
            return Err(bad("release points to a different historical batch"));
        }
        for (index, (row, selection)) in self.rows.iter().zip(&self.request.rows).enumerate() {
            row.validate()?;
            if row.ordinal != index
                || row.construct_id != selection.construct_id
                || serde_json::to_vec(&row.options).map_err(|e| bad(&e.to_string()))?
                    != serde_json::to_vec(&selection.options).map_err(|e| bad(&e.to_string()))?
            {
                return Err(bad(
                    "released row identity/options differ from explicit selection",
                ));
            }
            let declared = selected(&batch, &row.construct_id)?;
            let design = batch
                .request
                .design
                .as_ref()
                .ok_or_else(|| bad("release source design unavailable"))?;
            if row.combination_id != design::combination_id(design, &declared.choices)?
                || !digest(&row.declaration_revision)
            {
                return Err(bad("release source-choice identity changed"));
            }
            if row.inputs.len() != design.bins.len() {
                return Err(bad("released input count differs from source bins"));
            }
            for (bin, input) in design.bins.iter().zip(&row.inputs) {
                let selected = declared
                    .choices
                    .iter()
                    .find(|choice| choice.bin_id == bin.id)
                    .ok_or_else(|| bad("missing declared choice"))?;
                let choice = bin
                    .choices
                    .iter()
                    .find(|choice| choice.id == selected.choice_id)
                    .ok_or_else(|| bad("missing declared source alternative"))?;
                let (input_id, _) =
                    design::source_role_ids(&design.id, &declared.id, &bin.id, &choice.id)?;
                if input.input_id != input_id || input.event != choice.source {
                    return Err(bad(
                        "released plan is not bound to the selected source alternatives",
                    ));
                }
            }
            for input in &row.inputs {
                design::resolve_event(workspace, &input.event)?;
            }
            design::resolve_event(workspace, &row.target)?;
            design::resolve_event(workspace, &row.output.construct)?;
            let asset = workspace
                .attachments
                .get(&row.output.run.attachment_hash)
                .ok_or_else(|| bad("released run artifact unavailable"))?;
            if asset.kind != run_lineage::WORKSPACE_ATTACHMENT_KIND
                || blake3::hash(&asset.bytes).to_hex().as_str() != row.output.run.attachment_hash
            {
                return Err(bad("released run attachment differs"));
            }
            let run = run_lineage::load_run_record(&asset.bytes)?;
            if run.run_id != row.output.run.run_id
                || run.target != row.target
                || run.assessment.engine_version != self.engine_version
                || run.assessment.engine_version != row.output.run.recorded_engine
                || run.assessment.status != row.output.run.recorded_status
                || run.inputs.len() != row.inputs.len()
                || run
                    .inputs
                    .iter()
                    .zip(&row.inputs)
                    .any(|(a, b)| a.input_id != b.input_id || a.event != b.event)
            {
                return Err(bad("release target/input/run linkage differs"));
            }
            if run.products.len() != row.output.run.products.len()
                || run
                    .products
                    .iter()
                    .zip(&row.output.run.products)
                    .any(|(a, b)| {
                        a.step_index != b.step_index
                            || a.document_id != b.document_id
                            || a.event_id != b.event_id
                            || a.revision != b.revision
                    })
            {
                return Err(bad("release operation product events differ"));
            }
            for product in &run.products {
                design::resolve_event(
                    workspace,
                    &EventBinding {
                        document_id: product.document_id.clone(),
                        event_id: product.event_id.clone(),
                        revision: product.revision.clone(),
                    },
                )?;
            }
            let last = run
                .products
                .last()
                .ok_or_else(|| bad("released run has no final product"))?;
            let named = design::resolve_event(workspace, &row.output.construct)?;
            let original = workspace
                .revisions
                .get(&last.revision)
                .ok_or_else(|| bad("released final source unavailable"))?;
            if !cloning_outputs::ancestry(workspace, &row.output.construct.event_id, &last.event_id)
                || !cloning_outputs::same_reference(named, original)
            {
                return Err(bad(
                    "released named construct is not a reference-preserving descendant of its recorded final product",
                ));
            }
            if row.options.save_primers != row.output.primer_attachment.is_some() {
                return Err(bad("released primer output differs from selected options"));
            }
            if let Some(hash) = &row.output.primer_attachment {
                let asset = workspace
                    .attachments
                    .get(hash)
                    .ok_or_else(|| bad("released primer artifact unavailable"))?;
                if asset.kind != "oligo-table-v1"
                    || blake3::hash(&asset.bytes).to_hex().as_str() != hash
                {
                    return Err(bad("released primer collection differs"));
                }
                let table: crate::formats::oligos::OligoTable =
                    serde_json::from_slice(&asset.bytes).map_err(|e| bad(&e.to_string()))?;
                table.validate()?;
                if table.records.is_empty() {
                    return Err(bad("released primer collection is empty"));
                }
            }
        }
        count(self, MAX_MANIFEST_BYTES)
    }
}
impl Commit {
    pub fn validate_binding(&self, context: &str, source_fingerprint: &str) -> Result<(), Error> {
        if self.context_revision != context
            || self.source_workspace_fingerprint != source_fingerprint
        {
            return Err(bad("release install context changed"));
        }
        self.validate()
    }
    pub fn validate(&self) -> Result<(), Error> {
        if !digest(&self.context_revision)
            || !digest(&self.source_workspace_fingerprint)
            || !digest(&self.manifest_hash)
        {
            return Err(bad("invalid release commit identity"));
        }
        self.workspace.validate()?;
        let asset = self
            .workspace
            .attachments
            .get(&self.manifest_hash)
            .ok_or_else(|| bad("release manifest unavailable"))?;
        if asset.kind != KIND
            || asset.media_type != "application/json"
            || asset.bytes.len() > MAX_MANIFEST_BYTES
            || blake3::hash(&asset.bytes).to_hex().as_str() != self.manifest_hash
        {
            return Err(bad("release manifest hash or type differs"));
        }
        let manifest: Manifest =
            serde_json::from_slice(&asset.bytes).map_err(|e| bad(&e.to_string()))?;
        if manifest.context_revision != self.context_revision
            || manifest.source_workspace_fingerprint != self.source_workspace_fingerprint
        {
            return Err(bad("release installation receipt differs from manifest"));
        }
        manifest.validate(&self.workspace)
    }
}
/// Historical manifests only. Loading never reconstructs a private release session.
pub fn load(workspace: &Workspace, hash: &str) -> Result<Manifest, Error> {
    design::workspace_guard(workspace)?;
    if !digest(hash) {
        return Err(bad("invalid release manifest hash"));
    }
    let asset = workspace
        .attachments
        .get(hash)
        .ok_or_else(|| bad("release manifest unavailable"))?;
    if asset.kind != KIND
        || asset.media_type != "application/json"
        || asset.bytes.len() > MAX_MANIFEST_BYTES
        || blake3::hash(&asset.bytes).to_hex().as_str() != hash
    {
        return Err(bad("release manifest hash/type/size differs"));
    }
    let manifest: Manifest =
        serde_json::from_slice(&asset.bytes).map_err(|e| bad(&e.to_string()))?;
    manifest.validate(workspace)?;
    Ok(manifest)
}
pub fn list(workspace: &Workspace) -> Result<Vec<Summary>, Error> {
    design::workspace_guard(workspace)?;
    workspace
        .attachments
        .iter()
        .filter(|(_, asset)| asset.kind == KIND)
        .map(|(hash, _)| {
            let manifest = load(workspace, hash)?;
            Ok(Summary {
                manifest_hash: hash.clone(),
                context_revision: manifest.context_revision,
                batch_hash: manifest.batch_hash,
                row_count: manifest.rows.len(),
                recorded_engine: manifest.engine_version,
            })
        })
        .collect()
}
/// Optional later archive seam: never join a caller-controlled construct ID,
/// even legal native IDs such as '.' or '..', into an archive path.
pub fn archive_entry(ordinal: usize, construct_id: &str, extension: &str) -> Result<String, Error> {
    if ordinal >= design::MAX_SELECTED
        || construct_id.is_empty()
        || construct_id.len() > 48
        || !matches!(
            extension,
            "gb" | "fasta" | "html" | "pdf" | "csv" | "tsv" | "json"
        )
    {
        return Err(bad("unsupported bounded archive entry"));
    }
    let hash = blake3::hash(construct_id.as_bytes()).to_hex().to_string();
    Ok(format!(
        "constructs/{:04}-{}.{}",
        ordinal + 1,
        &hash[..24],
        extension
    ))
}
