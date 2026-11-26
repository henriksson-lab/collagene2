//! Local synchronization bookkeeping. Never upload this attachment as part of
//! a plasmid's publication closure, and never interpret its server as a URL to fetch.
use collagene_library_protocol::{Digest, DocumentId, LabId, RevisionId};
use serde::{Deserialize, Serialize};

pub const ATTACHMENT_KIND: &str = "lab-library-bindings-v1";
const MAX_BYTES: usize = 1024 * 1024;
const MAX_BINDINGS: usize = 4096;

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Remote {
    /// Canonical origin supplied by the connected transport, without credentials.
    pub server: String,
    pub lab: LabId,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Binding {
    pub remote: Remote,
    pub local_document: String,
    /// The exact local content last received from or successfully sent to this remote.
    pub local_revision: Digest,
    pub local_head_event: Option<String>,
    pub remote_document: DocumentId,
    /// CAS base belongs to the remote history, not the native molecule hash.
    pub remote_revision: RevisionId,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Bindings {
    version: u16,
    records: Vec<Binding>,
}

impl Default for Bindings {
    fn default() -> Self {
        Self {
            version: 1,
            records: Vec::new(),
        }
    }
}

impl Bindings {
    pub fn get(&self, remote: &Remote, local_document: &str) -> Option<&Binding> {
        self.records
            .iter()
            .find(|binding| &binding.remote == remote && binding.local_document == local_document)
    }

    /// Only record a validated server receipt. A different document requires an
    /// explicit unlink/save-copy action, never an accidental replacement of a CAS base.
    pub fn record(&mut self, binding: Binding) -> Result<(), String> {
        validate_binding(&binding)?;
        if let Some(previous) = self.records.iter_mut().find(|previous| {
            previous.remote == binding.remote && previous.local_document == binding.local_document
        }) {
            if previous.remote_document != binding.remote_document {
                return Err("This local plasmid is linked to another remote document; save a separate copy or unlink it first.".into());
            }
            *previous = binding;
        } else {
            if self.records.len() >= MAX_BINDINGS {
                return Err("The workspace has reached its remote-link limit.".into());
            }
            self.records.push(binding);
        }
        Ok(())
    }

    pub fn unlink(&mut self, remote: &Remote, local_document: &str) {
        self.records.retain(|binding| {
            &binding.remote != remote || binding.local_document != local_document
        });
    }

    pub fn encode(&self) -> Result<Vec<u8>, String> {
        self.validate()?;
        let bytes = serde_json::to_vec(self).map_err(|error| error.to_string())?;
        if bytes.len() > MAX_BYTES {
            return Err("Workspace remote links exceed their size limit.".into());
        }
        Ok(bytes)
    }

    pub fn decode(bytes: &[u8]) -> Result<Self, String> {
        if bytes.len() > MAX_BYTES {
            return Err("Workspace remote links exceed their size limit.".into());
        }
        let value: Self = serde_json::from_slice(bytes).map_err(|error| error.to_string())?;
        value.validate()?;
        Ok(value)
    }

    fn validate(&self) -> Result<(), String> {
        if self.version != 1 || self.records.len() > MAX_BINDINGS {
            return Err("Unsupported or oversized workspace remote links.".into());
        }
        let mut identities = std::collections::BTreeSet::new();
        for binding in &self.records {
            validate_binding(binding)?;
            if !identities.insert((
                &binding.remote.server,
                &binding.remote.lab,
                &binding.local_document,
            )) {
                return Err("Workspace contains duplicate remote links.".into());
            }
        }
        Ok(())
    }
}

fn validate_binding(binding: &Binding) -> Result<(), String> {
    fn native_id(value: &str) -> bool {
        !value.is_empty()
            && value.len() <= 128
            && value
                .bytes()
                .all(|byte| byte.is_ascii_alphanumeric() || matches!(byte, b'_' | b'-'))
    }
    let server = &binding.remote.server;
    if server.is_empty()
        || server.len() > 2048
        || server.chars().any(char::is_control)
        || server.contains(['@', '?', '#'])
        || !native_id(&binding.local_document)
        || binding
            .local_head_event
            .as_deref()
            .is_some_and(|event| !native_id(event))
    {
        return Err("Invalid local library link metadata.".into());
    }
    Ok(())
}
