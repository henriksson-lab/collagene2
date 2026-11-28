//! Immutable read-assembly evidence. Loading checks structure and exact source
//! events, never restores a scientific result as current verification authority.
use crate::{
    contig::{ContigReport, ReadInput},
    run_lineage::InputEventBinding,
    workspace::{RevisionEvent, Workspace, WorkspaceAttachment},
};
use collagene_model::Error;
use serde::{Deserialize, Serialize};
use std::{
    collections::{BTreeMap, BTreeSet},
    io::Write,
};
pub const WORKSPACE_ATTACHMENT_KIND: &str = "read-assembly-v1";
pub const MAX_SAVED_BYTES: usize = 17 * 1024 * 1024;
pub const MAX_SAVED_ASSEMBLIES: usize = 128;
#[derive(Debug, Clone, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct SavedAssembly {
    pub schema_version: u32,
    pub report: ContigReport,
    pub inputs: Vec<InputEventBinding>,
}
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct AssemblySummary {
    pub attachment_hash: String,
    pub batch_revision: String,
    pub engine_version: String,
    pub read_count: usize,
    pub contig_count: usize,
    pub recorded_complete: bool,
}
fn bad(text: &str) -> Error {
    Error::Format(format!("saved read assembly: {text}"))
}
fn digest(text: &str) -> bool {
    text.len() == 64
        && text
            .bytes()
            .all(|b| b.is_ascii_digit() || (b'a'..=b'f').contains(&b))
}
fn identifier(text: &str) -> bool {
    !text.is_empty()
        && text.len() <= 128
        && text
            .bytes()
            .all(|b| b.is_ascii_alphanumeric() || matches!(b, b'_' | b'-' | b'.' | b':'))
}
struct Limited(Vec<u8>);
impl Write for Limited {
    fn write(&mut self, b: &[u8]) -> std::io::Result<usize> {
        if b.len() > MAX_SAVED_BYTES.saturating_sub(self.0.len()) {
            return Err(std::io::Error::other("saved assembly byte limit"));
        }
        self.0.extend_from_slice(b);
        Ok(b.len())
    }
    fn flush(&mut self) -> std::io::Result<()> {
        Ok(())
    }
}
impl SavedAssembly {
    pub fn new(report: ContigReport, inputs: Vec<InputEventBinding>) -> Result<Self, Error> {
        let value = Self {
            schema_version: 1,
            report,
            inputs,
        };
        value.validate()?;
        Ok(value)
    }
    /// Historical structural validation only, including older producer engines.
    pub fn validate(&self) -> Result<(), Error> {
        if self.schema_version != 1
            || self.inputs.is_empty()
            || self.inputs.len() > 32
            || self.inputs.len() != self.report.reads.len()
        {
            return Err(bad("schema or read binding count"));
        }
        self.report.validate()?;
        let mut ids = BTreeSet::new();
        for (binding, read) in self.inputs.iter().zip(&self.report.reads) {
            if binding.input_id != read.id
                || binding.event.revision != read.source_revision
                || binding.input_id.is_empty()
                || binding.input_id.len() > 256
                || binding.input_id.chars().any(char::is_control)
                || !ids.insert(&binding.input_id)
                || !identifier(&binding.event.document_id)
                || !binding
                    .event
                    .event_id
                    .strip_prefix("event-")
                    .is_some_and(digest)
                || !digest(&binding.event.revision)
            {
                return Err(bad(
                    "input roles must bind the exact ordered read revisions",
                ));
            }
        }
        Ok(())
    }
}
pub fn make_attachment(record: &SavedAssembly) -> Result<(String, WorkspaceAttachment), Error> {
    record.validate()?;
    let mut out = Limited(Vec::new());
    serde_json::to_writer(&mut out, record).map_err(|_| Error::TooLarge)?;
    let hash = blake3::hash(&out.0).to_hex().to_string();
    Ok((
        hash.clone(),
        WorkspaceAttachment {
            kind: WORKSPACE_ATTACHMENT_KIND.into(),
            name: format!("read-assembly-{hash}.json"),
            media_type: "application/json".into(),
            bytes: out.0,
        },
    ))
}
pub fn load_record(bytes: &[u8]) -> Result<SavedAssembly, Error> {
    if bytes.len() > MAX_SAVED_BYTES {
        return Err(Error::TooLarge);
    }
    let record: SavedAssembly = serde_json::from_slice(bytes).map_err(|e| bad(&e.to_string()))?;
    let raw: serde_json::Value = serde_json::from_slice(bytes).map_err(|e| bad(&e.to_string()))?;
    if raw != serde_json::to_value(&record).map_err(|_| Error::TooLarge)? {
        return Err(bad("unknown or ignored record fields"));
    }
    record.validate()?;
    Ok(record)
}
fn event_index(w: &Workspace) -> BTreeMap<&str, &RevisionEvent> {
    w.history.iter().map(|e| (e.id.as_str(), e)).collect()
}
fn crossrefs(
    w: &Workspace,
    record: &SavedAssembly,
    events: &BTreeMap<&str, &RevisionEvent>,
) -> Result<(), Error> {
    for (binding, read) in record.inputs.iter().zip(&record.report.reads) {
        let e = events
            .get(binding.event.event_id.as_str())
            .ok_or_else(|| bad("bound source event is unavailable"))?;
        if e.document_id != binding.event.document_id || e.revision != binding.event.revision {
            return Err(bad("source event/document/revision mismatch"));
        }
        let source = w
            .revisions
            .get(&e.revision)
            .ok_or_else(|| bad("source molecular revision is unavailable"))?;
        if source.sequence.len() != read.source_length {
            return Err(bad("source length does not match retained evidence"));
        }
    }
    Ok(())
}
/// Called by Workspace::validate after generic content hashes and event validation.
/// Does not call Workspace::validate recursively or rerun scientific assembly.
pub fn validate_workspace_assemblies(w: &Workspace) -> Result<(), Error> {
    let attachments = w
        .attachments
        .values()
        .filter(|a| a.kind == WORKSPACE_ATTACHMENT_KIND)
        .collect::<Vec<_>>();
    if attachments.len() > MAX_SAVED_ASSEMBLIES {
        return Err(Error::TooLarge);
    }
    if attachments.is_empty() {
        return Ok(());
    }
    let events = event_index(w);
    for a in attachments {
        if a.media_type != "application/json" {
            return Err(bad("assembly attachment must use application/json"));
        }
        let record = load_record(&a.bytes)?;
        crossrefs(w, &record, &events)?;
    }
    Ok(())
}
pub fn list(w: &Workspace) -> Result<Vec<AssemblySummary>, Error> {
    w.validate()?;
    w.attachments
        .iter()
        .filter(|(_, a)| a.kind == WORKSPACE_ATTACHMENT_KIND)
        .map(|(hash, a)| {
            let record = load_record(&a.bytes)?;
            Ok(AssemblySummary {
                attachment_hash: hash.clone(),
                batch_revision: record.report.batch_revision,
                engine_version: record.report.engine_version,
                read_count: record.report.reads.len(),
                contig_count: record.report.contigs.len(),
                recorded_complete: record.report.complete,
            })
        })
        .collect()
}
/// Resolve the original retained events, independently of current document heads.
/// Call the current assembler explicitly on the returned snapshots before use.
pub fn resolve_inputs(w: &Workspace, record: &SavedAssembly) -> Result<Vec<ReadInput>, Error> {
    w.validate()?;
    record.validate()?;
    crossrefs(w, record, &event_index(w))?;
    record
        .inputs
        .iter()
        .map(|b| {
            Ok(ReadInput {
                id: b.input_id.clone(),
                molecule: w
                    .revisions
                    .get(&b.event.revision)
                    .ok_or_else(|| bad("source revision missing"))?
                    .clone(),
            })
        })
        .collect()
}
