//! Immutable event-bound pairwise comparison evidence. Saving retains a privately
//! computed report; loading never replays it or establishes a current assessment.
use crate::{
    alignment::review::Report,
    documents::workspace::{self, Workspace, WorkspaceAttachment},
    material::lineage::EventBinding,
};
use collagene_model::{Error, Molecule};
use serde::{Deserialize, Serialize};
use std::io::{self, Write};

pub const KIND: &str = "pairwise-review-v1";
pub const MAX_SAVED_BYTES: usize = workspace::MAX_ENTRY_BYTES;
const VERSION: u16 = 1;

#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Record {
    pub schema_version: u16,
    pub reference: EventBinding,
    pub query: EventBinding,
    pub report: Report,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Summary {
    pub attachment_hash: String,
    pub reference: EventBinding,
    pub query: EventBinding,
    pub context_revision: String,
    pub reference_name: String,
    pub query_name: String,
}

fn bad(message: &str) -> Error {
    Error::Format(format!("pairwise history: {message}"))
}
fn digest(value: &str) -> bool {
    value.len() == 64
        && value
            .bytes()
            .all(|b| b.is_ascii_digit() || (b'a'..=b'f').contains(&b))
}
fn binding_guard(binding: &EventBinding) -> Result<(), Error> {
    if binding.document_id.is_empty()
        || binding.document_id.len() > 128
        || !binding
            .document_id
            .bytes()
            .all(|b| b.is_ascii_alphanumeric() || matches!(b, b'_' | b'-' | b'.' | b':'))
        || !digest(&binding.revision)
        || !binding.event_id.strip_prefix("event-").is_some_and(digest)
    {
        return Err(bad("invalid exact event binding"));
    }
    Ok(())
}
fn workspace_guard(frame: &Workspace) -> Result<(), Error> {
    if frame.documents.len() > workspace::MAX_DOCUMENTS
        || frame.revisions.len() > workspace::MAX_REVISIONS
        || frame.history.len() > workspace::MAX_HISTORY_EVENTS
        || frame.candidates.len() > workspace::MAX_CANDIDATES
        || frame.attachments.len() > workspace::MAX_ATTACHMENTS
        || frame
            .attachments
            .values()
            .any(|asset| asset.bytes.len() > workspace::MAX_ENTRY_BYTES)
    {
        return Err(Error::TooLarge);
    }
    let mut minimum_bytes = 0usize;
    for bytes in frame
        .revisions
        .values()
        .map(|m| m.sequence.len())
        .chain(frame.attachments.values().map(|asset| asset.bytes.len()))
    {
        minimum_bytes = minimum_bytes
            .checked_add(bytes)
            .filter(|total| *total <= workspace::MAX_EXPANDED_BYTES)
            .ok_or(Error::TooLarge)?;
    }
    for molecule in frame.revisions.values() {
        collagene_model::bounds::model_preflight(&[molecule])?;
    }
    for candidate in &frame.candidates {
        if candidate.plan.inputs.len() > 1000 || candidate.plan.steps.len() > 1000 {
            return Err(Error::TooLarge);
        }
        for input in &candidate.plan.inputs {
            collagene_model::bounds::model_preflight(&[&input.molecule])?;
        }
    }
    frame.validate()
}
fn resolve<'a>(frame: &'a Workspace, binding: &EventBinding) -> Result<&'a Molecule, Error> {
    binding_guard(binding)?;
    if !frame
        .documents
        .iter()
        .any(|doc| doc.id == binding.document_id)
    {
        return Err(bad("the original document is not retained"));
    }
    let event = frame
        .history
        .iter()
        .find(|event| event.id == binding.event_id)
        .ok_or_else(|| bad("the exact source event is not retained"))?;
    if event.document_id != binding.document_id || event.revision != binding.revision {
        return Err(bad("source document/event/revision mismatch"));
    }
    frame
        .revisions
        .get(&binding.revision)
        .ok_or_else(|| bad("the source revision is not retained"))
}
fn validate_record(frame: &Workspace, record: &Record) -> Result<(), Error> {
    if record.schema_version != VERSION {
        return Err(Error::Unsupported(
            "unknown pairwise review archive version".into(),
        ));
    }
    let reference = resolve(frame, &record.reference)?;
    let query = resolve(frame, &record.query)?;
    record.report.validate_sources(reference, query)
}

struct Limited(Vec<u8>);
impl Write for Limited {
    fn write(&mut self, bytes: &[u8]) -> io::Result<usize> {
        if self.0.len().saturating_add(bytes.len()) > MAX_SAVED_BYTES {
            return Err(io::Error::other("pairwise review archive byte limit"));
        }
        self.0.extend_from_slice(bytes);
        Ok(bytes.len())
    }
    fn flush(&mut self) -> io::Result<()> {
        Ok(())
    }
}

/// Retain an already privately computed, event-bound comparison without running
/// alignment again. The host must check its live draft/epoch before installation.
/// This function changes only the attachment map in its returned workspace.
pub fn capture(
    frame: &Workspace,
    reference: &EventBinding,
    query: &EventBinding,
    report: &Report,
) -> Result<(Workspace, String), Error> {
    workspace_guard(frame)?;
    let source_reference = resolve(frame, reference)?;
    let source_query = resolve(frame, query)?;
    report.validate_sources(source_reference, source_query)?;
    // Borrow while serializing; malformed large reports are never cloned before
    // validation or the existing attachment-size boundary.
    #[derive(Serialize)]
    struct Borrowed<'a> {
        schema_version: u16,
        reference: &'a EventBinding,
        query: &'a EventBinding,
        report: &'a Report,
    }
    let mut encoded = Limited(Vec::new());
    serde_json::to_writer(
        &mut encoded,
        &Borrowed {
            schema_version: VERSION,
            reference,
            query,
            report,
        },
    )
    .map_err(|_| Error::TooLarge)?;
    let hash = blake3::hash(&encoded.0).to_hex().to_string();
    if let Some(existing) = frame.attachments.get(&hash) {
        if existing.kind != KIND
            || existing.media_type != "application/json"
            || existing.bytes != encoded.0
        {
            return Err(bad("same bytes are retained under another attachment type"));
        }
        return Ok((frame.clone(), hash));
    }
    if frame.attachments.len() >= workspace::MAX_ATTACHMENTS {
        return Err(Error::TooLarge);
    }
    let mut staged = frame.clone();
    staged.attachments.insert(
        hash.clone(),
        WorkspaceAttachment {
            kind: KIND.into(),
            name: format!("Pairwise review {}", &hash[..12]),
            media_type: "application/json".into(),
            bytes: encoded.0,
        },
    );
    staged.validate()?;
    Ok((staged, hash))
}

fn load_guarded(frame: &Workspace, hash: &str) -> Result<Record, Error> {
    if !digest(hash) {
        return Err(bad("invalid attachment identity"));
    }
    let asset = frame
        .attachments
        .get(hash)
        .ok_or_else(|| bad("saved comparison is unavailable"))?;
    if asset.kind != KIND
        || asset.media_type != "application/json"
        || asset.bytes.len() > MAX_SAVED_BYTES
    {
        return Err(bad(
            "invalid comparison attachment kind, media type or size",
        ));
    }
    if blake3::hash(&asset.bytes).to_hex().as_str() != hash {
        return Err(bad("comparison attachment hash mismatch"));
    }
    let record: Record = serde_json::from_slice(&asset.bytes).map_err(|e| bad(&e.to_string()))?;
    validate_record(frame, &record)?;
    Ok(record)
}

/// Open historical recorded evidence only. No current/fresh/verified flag is
/// persisted or restored; explicit worker recomputation is a separate action.
pub fn load(frame: &Workspace, hash: &str) -> Result<Record, Error> {
    workspace_guard(frame)?;
    load_guarded(frame, hash)
}

/// List every retained pairwise review within the existing workspace budget.
/// Records with invalid bindings fail explicitly instead of being silently hidden.
pub fn list(frame: &Workspace) -> Result<Vec<Summary>, Error> {
    workspace_guard(frame)?;
    let mut summaries = Vec::new();
    for (hash, asset) in &frame.attachments {
        if asset.kind != KIND {
            continue;
        }
        let record = load_guarded(frame, hash)?;
        summaries.push(Summary {
            attachment_hash: hash.clone(),
            reference_name: resolve(frame, &record.reference)?.name.clone(),
            query_name: resolve(frame, &record.query)?.name.clone(),
            context_revision: record.report.context_revision,
            reference: record.reference,
            query: record.query,
        });
    }
    Ok(summaries)
}
