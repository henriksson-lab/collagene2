//! Operator-selected, saved native documents in an isolated in-process MCP.
#![cfg(not(target_arch = "wasm32"))]
use collagene_core::{protocol, workspace};
use collagene_desktop_host::open_local_repository;
use collagene_lab_library::{DocumentId, PayloadKind, RevisionId};
use collagene_library_sharing::{self as sharing, Inventory};
use collagene_mcp::CollageneServer;
use std::{collections::BTreeSet, path::Path, str::FromStr};

pub const MAX_SELECTIONS: usize = 32;
const MAX_PRELOAD_BYTES: usize = 64 * 1024 * 1024;
#[derive(Debug, thiserror::Error)]
pub enum Error {
    #[error("select a saved local document as doc-32hex@revision64hex")]
    InvalidSelection,
    #[error("at most 32 unique saved document/revision selections are allowed")]
    SelectionLimit,
    #[error("native local repository is unavailable")]
    Repository,
    #[error("the exact selected saved revision is unavailable or invalid")]
    InvalidSavedRevision,
    #[error("the selected document dependency closure cannot be exported")]
    Closure,
    #[error("selected session data exceeds memory limits")]
    Capacity,
    #[error("managed MCP session storage configuration is invalid")]
    Storage,
    #[error("selected native data could not be preloaded")]
    Preload,
}
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct NativeSelection {
    pub document: DocumentId,
    pub revision: RevisionId,
}
impl NativeSelection {
    pub fn parse(text: &str) -> Result<Self, Error> {
        let (document, revision) = text.split_once('@').ok_or(Error::InvalidSelection)?;
        Ok(Self {
            document: DocumentId::parse(document).map_err(|_| Error::InvalidSelection)?,
            revision: RevisionId::parse(revision).map_err(|_| Error::InvalidSelection)?,
        })
    }
}
impl FromStr for NativeSelection {
    type Err = Error;
    fn from_str(s: &str) -> Result<Self, Error> {
        Self::parse(s)
    }
}
/// Secret-free description of exactly what entered the memory-only session.
/// The inventory is documentary scope, never inherited scientific authority.
pub struct Preloaded {
    pub selection: NativeSelection,
    pub inventory: Option<Inventory>,
    pub receipt: serde_json::Value,
}
pub struct Session {
    server: CollageneServer,
    pub preloaded: Vec<Preloaded>,
}
impl Session {
    /// Bare --mcp: no repository or managed scratch directory is opened.
    pub fn empty() -> Self {
        Self {
            server: CollageneServer::default(),
            preloaded: vec![],
        }
    }
    pub fn server(&self) -> &CollageneServer {
        &self.server
    }
    pub fn into_server(self) -> CollageneServer {
        self.server
    }
}
/// `root` is operator/shell configuration, never an MCP tool argument. Opening
/// the repository-only native factory validates its private root. No working
/// database is opened. Empty sessions do not access the configured root at all.
pub fn create(root: impl AsRef<Path>, selections: &[NativeSelection]) -> Result<Session, Error> {
    create_with_storage(root, selections, None)
}
/// Reuse the existing managed-session storage policy when the operator supplies
/// --storage-root. This path is never accepted from an MCP request.
pub fn create_with_storage(
    root: impl AsRef<Path>,
    selections: &[NativeSelection],
    managed_storage_root: Option<&Path>,
) -> Result<Session, Error> {
    if selections.len() > MAX_SELECTIONS {
        return Err(Error::SelectionLimit);
    }
    let mut seen = BTreeSet::new();
    for pin in selections {
        if !seen.insert((pin.document.as_str(), pin.revision.as_str())) {
            return Err(Error::SelectionLimit);
        }
    }
    let server = match managed_storage_root {
        Some(path) => CollageneServer::with_storage_root(path).map_err(|_| Error::Storage)?,
        None => CollageneServer::default(),
    };
    if selections.is_empty() {
        return Ok(Session {
            server,
            preloaded: vec![],
        });
    }
    let repository = open_local_repository(root).map_err(|_| Error::Repository)?;
    let mut preloaded = Vec::new();
    let mut bytes = 0usize;
    for pin in selections {
        let payload = repository
            .open_document(&pin.document, Some(&pin.revision))
            .map_err(|_| Error::InvalidSavedRevision)?;
        payload
            .validate()
            .map_err(|_| Error::InvalidSavedRevision)?;
        let (receipt, inventory) = match payload.kind {
            PayloadKind::WorkspaceArchive => {
                let selected = payload
                    .selection
                    .as_ref()
                    .ok_or(Error::InvalidSavedRevision)?;
                let frame = workspace::load_workspace(&payload.bytes)
                    .map_err(|_| Error::InvalidSavedRevision)?;
                let exported = sharing::document_closure(
                    &frame,
                    &sharing::Selection {
                        document_id: selected.document_id.clone(),
                        revision: selected.revision.as_str().into(),
                        head_event: selected.head_event.clone(),
                    },
                )
                .map_err(|_| Error::Closure)?;
                bytes = bytes
                    .checked_add(exported.archive.len())
                    .ok_or(Error::Capacity)?;
                if bytes > MAX_PRELOAD_BYTES {
                    return Err(Error::Capacity);
                }
                (
                    server
                        .preload_archive_bytes(&exported.archive)
                        .map_err(|_| Error::Preload)?,
                    Some(exported.inventory),
                )
            }
            PayloadKind::NativeDocument => {
                protocol::load_document(&payload.bytes).map_err(|_| Error::InvalidSavedRevision)?;
                bytes = bytes
                    .checked_add(payload.bytes.len())
                    .ok_or(Error::Capacity)?;
                if bytes > MAX_PRELOAD_BYTES {
                    return Err(Error::Capacity);
                }
                (
                    server
                        .preload_document_bytes(&payload.bytes)
                        .map_err(|_| Error::Preload)?,
                    None,
                )
            }
        };
        preloaded.push(Preloaded {
            selection: pin.clone(),
            inventory,
            receipt: receipt.structured_content.ok_or(Error::Preload)?,
        });
    }
    Ok(Session { server, preloaded })
}
pub async fn serve(session: Session) -> Result<(), Box<dyn std::error::Error + Send + Sync>> {
    collagene_mcp::serve_stdio(session.into_server()).await
}
