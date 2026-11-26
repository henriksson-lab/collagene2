use crate::{Fault, Result};
use collagene_desktop_sync_protocol as wire;
use collagene_library_sync as sync;
use collagene_repository::{self as repository, folders};
use serde::{Serialize, de::DeserializeOwned};
use std::io::Write;
struct Bounded(Vec<u8>);
impl Write for Bounded {
    fn write(&mut self, b: &[u8]) -> std::io::Result<usize> {
        if self
            .0
            .len()
            .checked_add(b.len())
            .is_none_or(|n| n > 1 << 20)
        {
            return Err(std::io::Error::other("desktop sync metadata bound"));
        }
        self.0.extend_from_slice(b);
        Ok(b.len())
    }
    fn flush(&mut self) -> std::io::Result<()> {
        Ok(())
    }
}
pub(crate) fn convert<T: Serialize, U: DeserializeOwned>(value: &T) -> Result<U> {
    let mut out = Bounded(vec![]);
    serde_json::to_writer(&mut out, value).map_err(|_| Fault::too_large())?;
    serde_json::from_slice(&out.0)
        .map_err(|_| Fault::invalid("Invalid native metadata projection."))
}
pub(crate) fn remote(r: &repository::RemoteConfig) -> Result<wire::Remote> {
    Ok(wire::Remote {
        id: wire::RemoteId::new(&r.id)
            .map_err(|_| Fault::invalid("Invalid configured remote identity."))?,
        name: wire::library::Text::new(&r.name).map_err(|_| Fault::too_large())?,
        endpoint: wire::library::Text::new(&r.endpoint).map_err(|_| Fault::too_large())?,
        expected_lab: r.lab.clone(),
    })
}
pub(crate) fn rejection(r: &sync::Rejection) -> Result<wire::Rejection> {
    let code = match r.code {
        sync::RejectionCode::Invalid => wire::library::ErrorCode::Invalid,
        sync::RejectionCode::NotFound => wire::library::ErrorCode::NotFound,
        sync::RejectionCode::AlreadyExists => wire::library::ErrorCode::AlreadyExists,
        sync::RejectionCode::NotEmpty => wire::library::ErrorCode::NotEmpty,
        sync::RejectionCode::TooLarge => wire::library::ErrorCode::TooLarge,
        sync::RejectionCode::Conflict => wire::library::ErrorCode::Conflict,
    };
    Ok(wire::Rejection {
        code,
        message: r.message.clone(),
        conflict: r.conflict.as_ref().map(convert).transpose()?,
    })
}
pub(crate) fn pending(p: &repository::OutboxMetadata) -> Result<wire::PendingPush> {
    let target = match &p.target {
        repository::OutboxTarget::Document(id) => wire::LocalTarget::Document(convert(id)?),
        repository::OutboxTarget::Folder(id) => wire::LocalTarget::Folder(convert(id)?),
        repository::OutboxTarget::DeletedFolder(id) => {
            wire::LocalTarget::DeletedFolder(convert(id)?)
        }
    };
    Ok(wire::PendingPush {
        job: wire::PushJobId::new(&p.id)
            .map_err(|_| Fault::invalid("Invalid pending identity."))?,
        remote: wire::RemoteId::new(&p.remote_id)
            .map_err(|_| Fault::invalid("Invalid pending remote identity."))?,
        target,
        local_revision: p.local_revision.as_ref().map(convert).transpose()?,
        attempted: p.attempted,
        rejection: match &p.rejection {
            None => None,
            Some(receipt) => match &receipt.outcome {
                sync::MutationOutcome::Rejected(rejected) => {
                    if !p.attempted {
                        return Err(Fault::invalid(
                            "An unattempted outbox entry carries a rejection.",
                        ));
                    }
                    Some(rejection(rejected)?)
                }
                sync::MutationOutcome::Applied(_) => {
                    return Err(Fault::invalid("Invalid retained pending outcome."));
                }
            },
        },
    })
}
pub(crate) fn receipt(r: &sync::Receipt) -> Result<wire::PushReceipt> {
    Ok(wire::PushReceipt {
        job: wire::PushJobId::new(r.id.as_str())
            .map_err(|_| Fault::invalid("Invalid receipt identity."))?,
        intent: r.intent.clone(),
        catalog_generation: r.catalog_generation,
        outcome: match &r.outcome {
            sync::MutationOutcome::Applied(_) => wire::PushOutcome::Applied,
            sync::MutationOutcome::Rejected(r) => wire::PushOutcome::Rejected(rejection(r)?),
        },
    })
}
pub(crate) fn tree(j: &folders::TreeJob) -> Result<wire::TreeStatus> {
    let p = j.progress();
    let phase = match p.phase {
        folders::TreePhase::Discovering => wire::TreePhase::Discovering,
        folders::TreePhase::Scanning => wire::TreePhase::Scanning,
        folders::TreePhase::Importing => wire::TreePhase::Importing,
        folders::TreePhase::Complete => wire::TreePhase::Complete,
        folders::TreePhase::Invalidated => wire::TreePhase::Invalidated,
        folders::TreePhase::Blocked => wire::TreePhase::Blocked,
    };
    let root = match &j.root {
        folders::TreeRoot::LabRoot => wire::TreeRoot::LabRoot,
        folders::TreeRoot::Folder { id, parent } => wire::TreeRoot::Folder {
            id: id.clone(),
            parent: parent.clone(),
        },
    };
    let issue = p
        .issue
        .as_ref()
        .map(|issue| -> Result<_> {
            Ok(match issue {
                folders::TreeIssue::RemoteUnavailable => wire::TreeIssue::RemoteUnavailable,
                folders::TreeIssue::CatalogChanged { actual } => {
                    wire::TreeIssue::CatalogChanged { actual: *actual }
                }
                folders::TreeIssue::RemoteRejected { code } => {
                    wire::TreeIssue::RemoteRejected { code: *code }
                }
                folders::TreeIssue::InvalidGraph => wire::TreeIssue::InvalidGraph,
                folders::TreeIssue::Limit => wire::TreeIssue::Limit,
                folders::TreeIssue::AlreadyBound => wire::TreeIssue::AlreadyBound,
                folders::TreeIssue::LocalParentChanged { folder } => {
                    wire::TreeIssue::LocalParentChanged {
                        folder: convert(folder)?,
                    }
                }
                folders::TreeIssue::LocalMutationRejected => wire::TreeIssue::LocalMutationRejected,
                folders::TreeIssue::LocalStateChanged => wire::TreeIssue::LocalStateChanged,
            })
        })
        .transpose()?;
    Ok(wire::TreeStatus {
        job: wire::TreeJobId::new(&j.id)
            .map_err(|_| Fault::invalid("Invalid tree job identity."))?,
        remote: wire::RemoteId::new(&j.remote.id)
            .map_err(|_| Fault::invalid("Invalid tree remote identity."))?,
        root,
        destination: j.destination.as_ref().map(convert).transpose()?,
        phase,
        catalog_generation: p.catalog_generation,
        scanned_nodes: u16::try_from(p.scanned_nodes).map_err(|_| Fault::too_large())?,
        folders_created: u16::try_from(p.folders_created).map_err(|_| Fault::too_large())?,
        documents_created: u16::try_from(p.documents_created).map_err(|_| Fault::too_large())?,
        issue,
    })
}
pub(crate) fn page<T>(
    items: Vec<T>,
    request: &wire::LocalPageRequest,
) -> Result<(Vec<T>, Option<u16>)> {
    let start = usize::from(request.offset);
    let end = start
        .saturating_add(usize::from(request.limit))
        .min(items.len());
    let next = if end < items.len() {
        Some(u16::try_from(end).map_err(|_| Fault::too_large())?)
    } else {
        None
    };
    Ok((
        items
            .into_iter()
            .skip(start)
            .take(usize::from(request.limit))
            .collect(),
        next,
    ))
}
