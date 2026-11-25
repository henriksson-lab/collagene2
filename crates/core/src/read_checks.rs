//! Immutable construct-owned read checks. Loading is historical evidence;
//! capture and replay alone invoke the existing fresh scientific validator.
use crate::{
    Error, Molecule, assembly_recipe,
    run_lineage::{self, EventBinding},
    trace_validation::{
        self, SiteState, TraceValidationOptions, TraceValidationReport, ValidationStatus,
    },
    workspace::{self, Workspace, WorkspaceAttachment},
};
use serde::{Deserialize, Serialize};
use std::{
    collections::{BTreeMap, BTreeSet},
    io::Write,
};
pub const ATTACHMENT_KIND: &str = "read-check-v1";
pub const MAX_SAVED_BYTES: usize = 17 << 20;
pub const MAX_CHECKS: usize = 128;
const MAX_READS: usize = 32;
const MAX_COLUMNS: usize = 200_000;
macro_rules! dto {($($t:tt)*)=>{#[derive(Clone,Debug,Serialize,Deserialize)]#[cfg_attr(feature="json-schema",derive(schemars::JsonSchema))]#[serde(deny_unknown_fields)]$($t)*};}
dto! {pub struct SavedCheck {pub schema_version:u32,pub context_revision:String,pub reference:EventBinding,pub reads:Vec<EventBinding>,#[cfg_attr(feature="json-schema",schemars(with="serde_json::Value"))] pub report:TraceValidationReport}}
dto! {pub struct Commit {#[cfg_attr(feature="json-schema",schemars(with="serde_json::Value"))] pub workspace:Workspace,pub context_revision:String,pub attachment_hash:String}}
dto! {pub struct Summary {pub attachment_hash:String,pub reference:EventBinding,pub read_count:usize,pub recorded_engine:String,#[cfg_attr(feature="json-schema",schemars(with="String"))] pub recorded_status:ValidationStatus}}
dto! {pub struct Checked {pub context_revision:String,pub saved:SavedCheck,#[cfg_attr(feature="json-schema",schemars(with="serde_json::Value"))] pub report:TraceValidationReport}}
fn bad(s: &str) -> Error {
    Error::Format(format!("read checks: {s}"))
}
fn hex(s: &str) -> bool {
    s.len() == 64
        && s.bytes()
            .all(|b| b.is_ascii_digit() || (b'a'..=b'f').contains(&b))
}
fn binding_guard(b: &EventBinding) -> Result<(), Error> {
    if b.document_id.is_empty()
        || b.document_id.len() > 128
        || !b
            .document_id
            .bytes()
            .all(|v| v.is_ascii_alphanumeric() || matches!(v, b'_' | b'-' | b'.' | b':'))
        || !hex(&b.revision)
        || !b.event_id.strip_prefix("event-").is_some_and(hex)
    {
        return Err(bad("invalid exact event binding"));
    }
    Ok(())
}
pub(crate) fn resolve<'a>(w: &'a Workspace, b: &EventBinding) -> Result<&'a Molecule, Error> {
    binding_guard(b)?;
    let event = w
        .history
        .iter()
        .find(|e| e.id == b.event_id)
        .ok_or_else(|| bad("missing retained event"))?;
    if event.document_id != b.document_id || event.revision != b.revision {
        return Err(bad("event/document/revision mismatch"));
    }
    w.revisions
        .get(&b.revision)
        .ok_or_else(|| bad("missing retained revision"))
}
pub(crate) fn workspace_guard(w: &Workspace) -> Result<(), Error> {
    if w.documents.len() > workspace::MAX_DOCUMENTS
        || w.revisions.len() > workspace::MAX_REVISIONS
        || w.history.len() > workspace::MAX_HISTORY_EVENTS
        || w.attachments.len() > workspace::MAX_ATTACHMENTS
        || w.candidates.len() > workspace::MAX_CANDIDATES
    {
        return Err(Error::TooLarge);
    }
    for m in w.revisions.values() {
        crate::model::bounds::model_preflight(&[m])?;
    }
    for c in &w.candidates {
        if c.plan.inputs.len() > 1000 || c.plan.steps.len() > 1000 {
            return Err(Error::TooLarge);
        }
        for i in &c.plan.inputs {
            crate::model::bounds::model_preflight(&[&i.molecule])?;
        }
    }
    if w.attachments
        .values()
        .filter(|a| a.kind == ATTACHMENT_KIND)
        .count()
        > MAX_CHECKS
    {
        return Err(Error::TooLarge);
    }
    w.validate()
}
fn bindings_guard(reference: &EventBinding, reads: &[EventBinding]) -> Result<(), Error> {
    binding_guard(reference)?;
    if reads.is_empty() || reads.len() > MAX_READS {
        return Err(bad("native read validation requires1..32 ordered reads"));
    }
    for b in reads {
        binding_guard(b)?;
    }
    Ok(())
}
fn inputs_guard(
    w: &Workspace,
    reference: &EventBinding,
    reads: &[EventBinding],
    options: &TraceValidationOptions,
) -> Result<(), Error> {
    trace_validation::validate_options(options)?;
    bindings_guard(reference, reads)?;
    workspace_guard(w)?;
    let mut models = Vec::with_capacity(reads.len() + 1);
    models.push(resolve(w, reference)?);
    for b in reads {
        models.push(resolve(w, b)?);
    }
    crate::model::bounds::model_preflight(&models)?;
    assembly_recipe::hash(&(models, options))?;
    Ok(())
}
fn saved_revision(
    reference: &EventBinding,
    reads: &[EventBinding],
    report: &TraceValidationReport,
) -> Result<String, Error> {
    assembly_recipe::hash(&(
        "read-check-evidence-1",
        reference,
        reads,
        &report.options,
        &report.engine_version,
        &report.model,
    ))
}
/// Installation token includes the entire current workspace. Saved evidence has
/// a separate stable event/options/engine identity so retries can deduplicate.
pub fn capture_revision(
    w: &Workspace,
    reference: &EventBinding,
    reads: &[EventBinding],
    options: &TraceValidationOptions,
) -> Result<String, Error> {
    inputs_guard(w, reference, reads, options)?;
    assembly_recipe::hash(&(
        "capture-read-check-1",
        run_lineage::workspace_fingerprint(w)?,
        reference,
        reads,
        options,
    ))
}
struct Bytes {
    bytes: Vec<u8>,
}
impl Write for Bytes {
    fn write(&mut self, b: &[u8]) -> std::io::Result<usize> {
        if self
            .bytes
            .len()
            .checked_add(b.len())
            .is_none_or(|n| n > MAX_SAVED_BYTES)
        {
            return Err(std::io::Error::other("saved read-check cap"));
        }
        self.bytes.extend_from_slice(b);
        Ok(b.len())
    }
    fn flush(&mut self) -> std::io::Result<()> {
        Ok(())
    }
}
fn encode(v: &impl Serialize) -> Result<Vec<u8>, Error> {
    let mut w = Bytes { bytes: Vec::new() };
    serde_json::to_writer(&mut w, v).map_err(|_| Error::TooLarge)?;
    Ok(w.bytes)
}
fn texts(values: &[String]) -> Result<(), Error> {
    if values.len() > 200_000 || values.iter().any(|v| v.len() > 65_536) {
        return Err(Error::TooLarge);
    }
    Ok(())
}
fn fraction(v: f64) -> bool {
    v.is_finite() && (0.0..=1.0).contains(&v)
}
pub(crate) fn report_guard(
    report: &TraceValidationReport,
    reference: &EventBinding,
    reads: &[EventBinding],
) -> Result<(), Error> {
    bindings_guard(reference, reads)?;
    trace_validation::validate_options(&report.options)?;
    if report.reference_revision != reference.revision
        || report.read_revisions.len() != reads.len()
        || report.reads.len() != reads.len()
        || report
            .read_revisions
            .iter()
            .zip(reads)
            .any(|(r, b)| r != &b.revision)
        || report.engine_version.is_empty()
        || report.engine_version.len() > 128
        || report.model.is_empty()
        || report.model.len() > 1024
        || report.consensus.is_empty()
        || report.consensus.len() > crate::model::MAX_SEQUENCE_LEN
        || report.sites.len() > report.options.max_consensus_sites
        || report.insertions.len() > MAX_COLUMNS
        || report.uncovered_ranges.len() > MAX_COLUMNS + 1
        || !fraction(report.reference_coverage_fraction)
        || report.cells_charged > report.options.max_total_cells
        || report.exact_comparisons > report.options.max_exact_comparisons
    {
        return Err(bad("invalid report identity, limits or dimensions"));
    }
    if !report
        .consensus
        .bytes()
        .all(|b| b"ACGTRYSWKMBDHVN".contains(&b))
    {
        return Err(bad("invalid consensus calls"));
    }
    let n = report.consensus.len();
    let mut columns = 0usize;
    for (read, binding) in report.reads.iter().zip(reads) {
        if read.source_revision != binding.revision
            || read.name.len() > 1024
            || !fraction(read.paired_reference_fraction)
            || read.columns.len() > MAX_COLUMNS
            || read.source_hash.as_ref().is_some_and(|s| !hex(s))
            || read.error.as_ref().is_some_and(|s| s.len() > 65_536)
        {
            return Err(bad("invalid read identity or coordinates"));
        }
        texts(&read.warnings)?;
        columns = columns
            .checked_add(read.columns.len())
            .filter(|n| *n <= MAX_COLUMNS)
            .ok_or(Error::TooLarge)?;
        if read.reference_window_start > n || read.reference_window_length > n.saturating_mul(2) {
            return Err(bad("invalid reference window"));
        }
        let calls = read.source_call_length.unwrap_or(0);
        if calls > crate::traces::MAX_TRACE_BASES
            || read
                .source_retained_start
                .zip(read.source_retained_end)
                .is_some_and(|(a, b)| a > b || b > calls)
        {
            return Err(bad("invalid retained call coordinates"));
        }
        if let Some(trim) = &read.trim {
            if trim.original_calls > calls
                || trim.start > trim.end
                || trim.end > trim.original_calls
                || trim.low_or_unknown_calls.len() > calls
                || trim.low_or_unknown_calls.iter().any(|i| *i >= calls)
            {
                return Err(bad("invalid trim coordinates"));
            }
            texts(&trim.warnings)?;
        }
        for c in &read.columns {
            if c.reference_position.is_some_and(|i| i >= n)
                || c.reference_boundary > n
                || c.call_index.is_some_and(|i| i >= calls)
                || c.stored_call_index.is_some_and(|i| i >= calls)
                || c.reference_base
                    .is_some_and(|b| !b.is_ascii() || crate::analysis::base_mask(b as u8) == 0)
                || c.observed_base
                    .is_some_and(|b| !b.is_ascii() || crate::analysis::base_mask(b as u8) == 0)
            {
                return Err(bad("invalid alignment column"));
            }
        }
        if let Some(a) = &read.alignment {
            if a.reference_aligned.len() > MAX_COLUMNS
                || a.query_aligned.len() != a.reference_aligned.len()
                || a.operations.len() > MAX_COLUMNS
                || a.differences.len() > MAX_COLUMNS
                || a.cigar.len() > MAX_COLUMNS * 12
                || a.algorithm.len() > 1024
                || a.reference_start > a.reference_end
                || a.reference_end > read.reference_window_length
                || a.query_start > a.query_end
                || a.query_end > calls
                || !hex(&a.reference_hash)
                || !hex(&a.query_hash)
            {
                return Err(bad("invalid retained alignment shape"));
            }
            if a.stats.columns != a.reference_aligned.len()
                || !fraction(a.stats.reference_coverage)
                || !fraction(a.stats.query_coverage)
                || a.stats.certain_identity.is_some_and(|v| !fraction(v))
                || a.stats.compatible_identity.is_some_and(|v| !fraction(v))
                || !a
                    .reference_aligned
                    .bytes()
                    .chain(a.query_aligned.bytes())
                    .all(|b| b == b'-' || b"ACGTRYSWKMBDHVN".contains(&b))
            {
                return Err(bad("invalid alignment calls or statistics"));
            }
            let consumed = a.operations.iter().try_fold(0usize, |sum, run| {
                if run.length == 0 {
                    None
                } else {
                    sum.checked_add(run.length)
                }
            });
            if consumed != Some(a.stats.columns) {
                return Err(bad("invalid alignment run lengths"));
            }
            for d in &a.differences {
                if d.reference_start > d.reference_end
                    || d.reference_end > read.reference_window_length
                    || d.query_start > d.query_end
                    || d.query_end > calls
                    || d.reference_bases.len() != d.reference_end - d.reference_start
                    || d.query_bases.len() != d.query_end - d.query_start
                {
                    return Err(bad("invalid alignment difference coordinates"));
                }
            }
            texts(&a.warnings)?;
        }
    }
    let mut last = None;
    let mut coverage = 0usize;
    let mut discrepancy = false;
    for s in &report.sites {
        if s.position >= n
            || last.is_some_and(|i| s.position <= i)
            || !s.reference_base.is_ascii()
            || crate::analysis::base_mask(s.reference_base as u8) == 0
            || !s.consensus_base.is_ascii()
            || report.consensus.as_bytes()[s.position] != s.consensus_base as u8
            || s.eligible_read_indices.len() > reads.len()
            || s.eligible_read_indices.iter().any(|i| *i >= reads.len())
            || s.eligible_read_indices
                .iter()
                .collect::<BTreeSet<_>>()
                .len()
                != s.eligible_read_indices.len()
            || s.counts.iter().any(|c| *c > reads.len())
            || s.forward_reads > reads.len()
            || s.reverse_reads > reads.len()
        {
            return Err(bad("invalid consensus site"));
        }
        last = Some(s.position);
        if matches!(
            s.state,
            SiteState::SupportedMatch
                | SiteState::SupportedDifference
                | SiteState::ReferenceAmbiguous
        ) {
            coverage += 1;
        }
        discrepancy |= matches!(
            s.state,
            SiteState::SupportedDifference
                | SiteState::ConflictingCalls
                | SiteState::DeletionEvidence
        );
    }
    let mut insertion_groups: BTreeMap<(usize, &str), BTreeSet<usize>> = BTreeMap::new();
    let mut insertion_calls = 0usize;
    for i in &report.insertions {
        insertion_calls = insertion_calls
            .checked_add(i.sequence.len())
            .filter(|n| *n <= MAX_COLUMNS)
            .ok_or(Error::TooLarge)?;
        if i.reference_boundary > n
            || i.read_index >= reads.len()
            || i.sequence.is_empty()
            || i.sequence.len() != i.call_indices.len()
            || !i.sequence.bytes().all(|b| b"ACGTRYSWKMBDHVN".contains(&b))
            || i.call_indices
                .iter()
                .any(|c| *c >= report.reads[i.read_index].source_call_length.unwrap_or(0))
        {
            return Err(bad("invalid insertion evidence"));
        }
        if i.eligible {
            insertion_groups
                .entry((i.reference_boundary, &i.sequence))
                .or_default()
                .insert(i.read_index);
        }
    }
    discrepancy |= insertion_groups.values().any(|v| {
        v.len() >= report.options.min_support
            && (!report.options.require_both_strands
                || (v
                    .iter()
                    .any(|i| report.reads[*i].original_read_reverse_to_reference)
                    && v.iter()
                        .any(|i| !report.reads[*i].original_read_reverse_to_reference)))
    });
    let mut end = 0;
    let mut uncovered = 0usize;
    for &(a, b) in &report.uncovered_ranges {
        if a >= b || b > n || a < end {
            return Err(bad("invalid uncovered interval"));
        }
        end = b;
        uncovered += b - a;
    }
    if report.eligible_coverage != coverage
        || uncovered != n - coverage
        || (report.reference_coverage_fraction - coverage as f64 / n as f64).abs() > 1e-12
    {
        return Err(bad("coverage receipt disagrees with retained sites"));
    }
    let expected = if discrepancy {
        ValidationStatus::Discrepancies
    } else if coverage == n
        && report
            .sites
            .iter()
            .all(|s| s.state == SiteState::SupportedMatch)
        && report
            .reads
            .iter()
            .all(|r| r.error.is_none() && r.placement_unique_proven && r.search_complete)
    {
        ValidationStatus::ConsistentUnderRecordedThresholds
    } else {
        ValidationStatus::Incomplete
    };
    if report.status != expected {
        return Err(bad(
            "reported scientific status disagrees with retained evidence",
        ));
    }
    texts(&report.warnings)?;
    encode(report)?;
    Ok(())
}
impl SavedCheck {
    pub(crate) fn validate(&self) -> Result<(), Error> {
        if self.schema_version != 1 || !hex(&self.context_revision) {
            return Err(bad("invalid saved check schema"));
        }
        report_guard(&self.report, &self.reference, &self.reads)?;
        if self.context_revision != saved_revision(&self.reference, &self.reads, &self.report)? {
            return Err(bad("saved context mismatch"));
        }
        encode(self)?;
        Ok(())
    }
}
/// Load bounded historical evidence. No current scientific status is restored.
pub fn load(w: &Workspace, hash: &str) -> Result<SavedCheck, Error> {
    if !hex(hash) {
        return Err(bad("invalid check attachment hash"));
    }
    workspace_guard(w)?;
    load_bounded(w, hash)
}
fn load_bounded(w: &Workspace, hash: &str) -> Result<SavedCheck, Error> {
    let a = w
        .attachments
        .get(hash)
        .ok_or_else(|| bad("missing read check"))?;
    if a.kind != ATTACHMENT_KIND
        || a.media_type != "application/json"
        || a.bytes.len() > MAX_SAVED_BYTES
        || blake3::hash(&a.bytes).to_hex().as_str() != hash
    {
        return Err(bad("invalid check attachment"));
    }
    let saved: SavedCheck =
        serde_json::from_slice(&a.bytes).map_err(|_| bad("malformed saved check"))?;
    saved.validate()?;
    let reference = resolve(w, &saved.reference)?;
    if reference.sequence.len() != saved.report.consensus.len() {
        return Err(bad("retained reference length changed"));
    }
    for s in &saved.report.sites {
        if reference.sequence.as_str().as_bytes()[s.position] != s.reference_base as u8 {
            return Err(bad("retained reference calls changed"));
        }
    }
    for b in &saved.reads {
        resolve(w, b)?;
    }
    Ok(saved)
}
pub fn list(w: &Workspace) -> Result<Vec<Summary>, Error> {
    workspace_guard(w)?;
    w.attachments
        .iter()
        .filter(|(_, a)| a.kind == ATTACHMENT_KIND)
        .map(|(hash, _)| {
            let s = load_bounded(w, hash)?;
            Ok(Summary {
                attachment_hash: hash.clone(),
                reference: s.reference,
                read_count: s.reads.len(),
                recorded_engine: s.report.engine_version,
                recorded_status: s.report.status,
            })
        })
        .collect()
}
pub fn capture(
    w: &Workspace,
    reference: EventBinding,
    reads: Vec<EventBinding>,
    options: &TraceValidationOptions,
) -> Result<Commit, Error> {
    let context_revision = capture_revision(w, &reference, &reads, options)?;
    let actual = resolve(w, &reference)?;
    let molecules = reads
        .iter()
        .map(|b| resolve(w, b).cloned())
        .collect::<Result<Vec<_>, _>>()?;
    let report = trace_validation::validate_reads(actual, &molecules, options)?;
    let saved = SavedCheck {
        schema_version: 1,
        context_revision: saved_revision(&reference, &reads, &report)?,
        reference,
        reads,
        report,
    };
    saved.validate()?;
    let bytes = encode(&saved)?;
    let hash = blake3::hash(&bytes).to_hex().to_string();
    if !w.attachments.contains_key(&hash)
        && w.attachments
            .values()
            .filter(|a| a.kind == ATTACHMENT_KIND)
            .count()
            >= MAX_CHECKS
    {
        return Err(Error::TooLarge);
    }
    let mut next = w.clone();
    if let Some(a) = next.attachments.get(&hash) {
        if a.kind != ATTACHMENT_KIND || a.bytes != bytes {
            return Err(bad("attachment identity conflict"));
        }
    } else {
        next.attachments.insert(
            hash.clone(),
            WorkspaceAttachment {
                kind: ATTACHMENT_KIND.into(),
                name: format!("Read check {}", &saved.context_revision[..12]),
                media_type: "application/json".into(),
                bytes,
            },
        );
    }
    let commit = Commit {
        workspace: next,
        context_revision: context_revision.clone(),
        attachment_hash: hash,
    };
    commit.validate_binding(&context_revision)?;
    Ok(commit)
}
impl Commit {
    pub fn validate_binding(&self, expected: &str) -> Result<(), Error> {
        if !hex(expected) || self.context_revision != expected {
            return Err(bad("capture request changed"));
        }
        load(&self.workspace, &self.attachment_hash)?;
        Ok(())
    }
}
pub fn replay_revision(w: &Workspace, hash: &str) -> Result<String, Error> {
    let saved = load(w, hash)?;
    inputs_guard(w, &saved.reference, &saved.reads, &saved.report.options)?;
    assembly_recipe::hash(&(
        "replay-read-check-1",
        run_lineage::workspace_fingerprint(w)?,
        hash,
        crate::verification::ENGINE_VERSION,
        trace_validation::MODEL,
    ))
}
pub fn replay(w: &Workspace, hash: &str) -> Result<Checked, Error> {
    let context_revision = replay_revision(w, hash)?;
    let saved = load(w, hash)?;
    let reads = saved
        .reads
        .iter()
        .map(|b| resolve(w, b).cloned())
        .collect::<Result<Vec<_>, _>>()?;
    let report = trace_validation::validate_reads(
        resolve(w, &saved.reference)?,
        &reads,
        &saved.report.options,
    )?;
    let checked = Checked {
        context_revision: context_revision.clone(),
        saved,
        report,
    };
    checked.validate_binding(&context_revision)?;
    Ok(checked)
}
impl Checked {
    pub fn validate_binding(&self, expected: &str) -> Result<(), Error> {
        if !hex(expected) || self.context_revision != expected {
            return Err(bad("replay request changed"));
        }
        self.saved.validate()?;
        report_guard(&self.report, &self.saved.reference, &self.saved.reads)?;
        if self.report.options != self.saved.report.options
            || self.report.engine_version != crate::verification::ENGINE_VERSION
            || self.report.model != trace_validation::MODEL
            || self.report.consensus.len() != self.saved.report.consensus.len()
        {
            return Err(bad("fresh replay engine/options/reference mismatch"));
        }
        assembly_recipe::hash(self)?;
        Ok(())
    }
}
