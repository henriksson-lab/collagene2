//! One guarded durable current/previous transaction; callers own cancellation and
//! restore only after this function reaches the IndexedDB commit point.
use crate::worker::CachedWorker;
use collagene_core::protocol::{AnalysisJob, AnalysisResult, Request};
use collagene_core::workspace::Workspace;
use leptos::prelude::*;
use std::{cell::RefCell, rc::Rc};

pub async fn persist(
    database: &crate::storage::Database,
    next: Workspace,
    previous: Workspace,
    validate: Callback<Workspace, Result<(), String>>,
    encoder: Rc<CachedWorker>,
    guard: impl Fn() -> Result<(), String>,
    transaction: Rc<RefCell<Option<crate::storage::Transaction>>>,
) -> Result<Workspace, String> {
    guard()?;
    validate.run(next.clone())?;
    let current = encode(&encoder, next.clone(), u32::MAX - 560).await?;
    guard()?;
    let previous = encode(&encoder, previous, u32::MAX - 561).await?;
    guard()?;
    crate::storage::replace_workspace_cancellable(database, &current, &previous, transaction)
        .await?;
    Ok(next)
}

async fn encode(encoder: &CachedWorker, frame: Workspace, id: u32) -> Result<Vec<u8>, String> {
    let revision = frame
        .head(&frame.active_document)
        .map_err(|e| e.to_string())?
        .revision();
    let response = encoder
        .ask(Request {
            id,
            revision: revision.clone(),
            job: AnalysisJob::SaveWorkspace { workspace: frame },
        })
        .await?;
    if response.id != id || response.revision != revision {
        return Err("Workspace encoding response differs from its request.".into());
    }
    match response.result {
        AnalysisResult::WorkspaceSaved(bytes) => Ok(bytes),
        AnalysisResult::Error(error) => Err(error),
        _ => Err("Unexpected workspace encoding response.".into()),
    }
}
