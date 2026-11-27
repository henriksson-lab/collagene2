//! Exact-event RNA hypotheses. Archives store declarations, never cached quantities
//! or a material/experimental certificate. Reopening requires explicit analysis.
use crate::{
    Error, Molecule, rna_chemistry as rna,
    run_lineage::EventBinding,
    workspace::{Workspace, WorkspaceAttachment},
};
use serde::{Deserialize, Serialize};
use std::io::{self, Write};
pub const KIND: &str = "rna-declaration-v1";
pub const VERSION: u16 = 1;
pub const MAX_BYTES: usize = rna::MAX_DECLARATION_BYTES + 4096;
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Record {
    pub version: u16,
    pub name: String,
    pub source: EventBinding,
    pub declaration: rna::Declaration,
}
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Summary {
    pub hash: String,
    pub name: String,
    pub source: EventBinding,
    pub source_name: String,
}
fn bad(s: impl std::fmt::Display) -> Error {
    Error::Format(format!("RNA declaration: {s}"))
}
fn digest(s: &str) -> bool {
    s.len() == 64
        && s.bytes()
            .all(|b| b.is_ascii_digit() || (b'a'..=b'f').contains(&b))
}
pub fn source<'a>(w: &'a Workspace, b: &EventBinding) -> Result<&'a Molecule, Error> {
    crate::restriction_history::source_role(b)?;
    if !w.documents.iter().any(|d| d.id == b.document_id) {
        return Err(bad("source document is unavailable"));
    }
    let e = w
        .history
        .iter()
        .find(|e| e.id == b.event_id)
        .ok_or_else(|| bad("exact source event is unavailable"))?;
    if e.document_id != b.document_id || e.revision != b.revision {
        return Err(bad("source event/document/revision differs"));
    }
    w.revisions
        .get(&b.revision)
        .ok_or_else(|| bad("source revision is unavailable"))
}
fn validate(w: &Workspace, r: &Record) -> Result<(), Error> {
    if r.version != VERSION {
        return Err(Error::Unsupported("RNA declaration record version".into()));
    }
    if r.name.trim().is_empty() || r.name.len() > 256 || r.name.chars().any(char::is_control) {
        return Err(bad("a name of at most 256 bytes is required"));
    }
    rna::declaration_revision(source(w, &r.source)?, &r.declaration).map_err(bad)?;
    Ok(())
}
struct Limited(Vec<u8>);
impl Write for Limited {
    fn write(&mut self, b: &[u8]) -> io::Result<usize> {
        if b.len() > MAX_BYTES.saturating_sub(self.0.len()) {
            return Err(io::Error::other("RNA record size"));
        }
        self.0.extend_from_slice(b);
        Ok(b.len())
    }
    fn flush(&mut self) -> io::Result<()> {
        Ok(())
    }
}
fn decode(w: &Workspace, hash: &str) -> Result<Record, Error> {
    if !digest(hash) {
        return Err(bad("invalid attachment hash"));
    }
    let a = w
        .attachments
        .get(hash)
        .ok_or_else(|| bad("saved declaration is unavailable"))?;
    if a.kind != KIND
        || a.media_type != "application/json"
        || a.bytes.len() > MAX_BYTES
        || blake3::hash(&a.bytes).to_hex().as_str() != hash
    {
        return Err(bad("attachment role/size/hash differs"));
    }
    serde_json::from_slice(&a.bytes).map_err(bad)
}
/// Called by native archive validation; never calls Workspace::validate recursively.
pub(crate) fn validate_workspace_records(w: &Workspace) -> Result<(), Error> {
    for (hash, a) in &w.attachments {
        if a.kind == KIND {
            validate(w, &decode(w, hash)?)?;
        }
    }
    Ok(())
}
pub fn capture(w: &Workspace, r: &Record) -> Result<(Workspace, String), Error> {
    w.validate()?;
    validate(w, r)?;
    let mut bytes = Limited(Vec::new());
    serde_json::to_writer(&mut bytes, r).map_err(|_| Error::TooLarge)?;
    let hash = blake3::hash(&bytes.0).to_hex().to_string();
    if let Some(a) = w.attachments.get(&hash) {
        if a.kind != KIND || a.media_type != "application/json" || a.bytes != bytes.0 {
            return Err(bad("same bytes have another attachment role"));
        }
        return Ok((w.clone(), hash));
    }
    if w.attachments.len() >= crate::workspace::MAX_ATTACHMENTS {
        return Err(Error::TooLarge);
    }
    let mut next = w.clone();
    next.attachments.insert(
        hash.clone(),
        WorkspaceAttachment {
            kind: KIND.into(),
            name: r.name.clone(),
            media_type: "application/json".into(),
            bytes: bytes.0,
        },
    );
    next.validate()?;
    Ok((next, hash))
}
pub fn load(w: &Workspace, hash: &str) -> Result<Record, Error> {
    w.validate()?;
    let r = decode(w, hash)?;
    validate(w, &r)?;
    Ok(r)
}
pub fn list(w: &Workspace) -> Result<Vec<Summary>, Error> {
    w.validate()?;
    let mut out = vec![];
    for (hash, a) in &w.attachments {
        if a.kind == KIND {
            let r = decode(w, hash)?;
            out.push(Summary {
                hash: hash.clone(),
                name: r.name,
                source_name: source(w, &r.source)?.name.clone(),
                source: r.source,
            });
        }
    }
    Ok(out)
}
pub fn delete(w: &Workspace, hash: &str) -> Result<Workspace, Error> {
    load(w, hash)?;
    let mut next = w.clone();
    next.attachments.remove(hash);
    next.validate()?;
    Ok(next)
}
