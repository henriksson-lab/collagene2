//! Trusted native adapter. Caller authenticates membership for every invocation.
//! Only ordinary public mutations are dispatched; repository checkpoint callbacks
//! and local scope administration are deliberately absent from this interface.
use crate::*;
use collagene_lab_library::{self as store, mutations as durable};
use collagene_library_host as host;
use serde::{Serialize, de::DeserializeOwned};
use std::io::Write;

#[derive(Debug, thiserror::Error)]
pub enum DispatchError {
    #[error("invalid sync request/reply")]
    Wire,
    #[error("request expected a different authenticated lab")]
    WrongLab,
    #[error("reliable library service unavailable")]
    Unavailable,
}
struct Bounded(Vec<u8>);
impl Write for Bounded {
    fn write(&mut self, b: &[u8]) -> std::io::Result<usize> {
        if b.len() > (128usize << 10).saturating_sub(self.0.len()) {
            return Err(std::io::Error::other("sync receipt limit"));
        }
        self.0.extend_from_slice(b);
        Ok(b.len())
    }
    fn flush(&mut self) -> std::io::Result<()> {
        Ok(())
    }
}
fn convert<T: Serialize, U: DeserializeOwned>(v: &T) -> Result<U, DispatchError> {
    let mut b = Bounded(Vec::new());
    serde_json::to_writer(&mut b, v).map_err(|_| DispatchError::Unavailable)?;
    serde_json::from_slice(&b.0).map_err(|_| DispatchError::Unavailable)
}
fn folder(id: &FolderId) -> Result<store::FolderId, DispatchError> {
    store::FolderId::parse(id.as_str()).map_err(|_| DispatchError::Wire)
}
fn document(id: &DocumentId) -> Result<store::DocumentId, DispatchError> {
    store::DocumentId::parse(id.as_str()).map_err(|_| DispatchError::Wire)
}
fn revision(id: &RevisionId) -> Result<store::RevisionId, DispatchError> {
    store::RevisionId::parse(id.as_str()).map_err(|_| DispatchError::Wire)
}
fn payload(p: &Payload) -> Result<store::Payload, DispatchError> {
    Ok(store::Payload {
        kind: convert(&p.kind)?,
        bytes: p.bytes.as_slice().to_vec(),
        archive_hash: store::Digest::parse(p.archive_hash.as_str())
            .map_err(|_| DispatchError::Wire)?,
        selection: p.selection.as_ref().map(convert).transpose()?,
        sharing: p.sharing.as_ref().map(convert).transpose()?,
    })
}
/// Shape/integrity translation only; never grants membership, publication scope
/// or native material authority. Dispatch performs native review before first use.
pub fn to_storage_mutation(m: &Mutation) -> Result<durable::Mutation, DispatchError> {
    m.validate().map_err(|_| DispatchError::Wire)?;
    use MutationOperation as W;
    use durable::MutationOperation as S;
    let operation = match &m.operation {
        W::CreateFolder { parent, name } => S::CreateFolder {
            parent: parent.as_ref().map(folder).transpose()?,
            name: name.to_string(),
        },
        W::RenameFolder {
            id,
            expected_generation,
            name,
        } => S::RenameFolder {
            id: folder(id)?,
            expected_generation: *expected_generation,
            name: name.to_string(),
        },
        W::MoveFolder {
            id,
            expected_generation,
            parent,
        } => S::MoveFolder {
            id: folder(id)?,
            expected_generation: *expected_generation,
            parent: parent.as_ref().map(folder).transpose()?,
        },
        W::DeleteEmptyFolder {
            id,
            expected_generation,
        } => S::DeleteEmptyFolder {
            id: folder(id)?,
            expected_generation: *expected_generation,
        },
        W::ImportDocument {
            folder: destination,
            name,
            payload: p,
        } => S::ImportDocument {
            folder: destination.as_ref().map(folder).transpose()?,
            name: name.to_string(),
            payload: payload(p)?,
        },
        W::SaveRevision {
            id,
            base_revision,
            payload: p,
        } => S::SaveRevision {
            id: document(id)?,
            base_revision: revision(base_revision)?,
            payload: payload(p)?,
        },
        W::RenameDocument {
            id,
            expected_generation,
            name,
        } => S::RenameDocument {
            id: document(id)?,
            expected_generation: *expected_generation,
            name: name.to_string(),
        },
        W::MoveDocument {
            id,
            expected_generation,
            folder: destination,
        } => S::MoveDocument {
            id: document(id)?,
            expected_generation: *expected_generation,
            folder: destination.as_ref().map(folder).transpose()?,
        },
    };
    Ok(durable::Mutation {
        id: durable::MutationId::parse(m.id.as_str()).map_err(|_| DispatchError::Wire)?,
        operation,
    })
}
fn failure(code: library::ErrorCode, text: &str) -> ResponseOutcome {
    ResponseOutcome::Failure(library::Failure {
        code,
        message: Text::new(text).expect("bounded static message"),
        conflict: None,
    })
}
fn storage_failure(id: &MutationId, e: durable::MutationError) -> ResponseOutcome {
    match e {
        durable::MutationError::IntentMismatch => ResponseOutcome::IntentMismatch {
            mutation_id: id.clone(),
        },
        durable::MutationError::Invalid(_) => {
            failure(library::ErrorCode::Invalid, "Invalid reliable mutation.")
        }
        durable::MutationError::Limit(_) => failure(
            library::ErrorCode::TooLarge,
            "Reliable mutation storage limit reached; retain the pending operation.",
        ),
        _ => failure(
            library::ErrorCode::Unavailable,
            "Reliable library storage is unavailable; retain the same pending mutation ID.",
        ),
    }
}
fn execute(
    db: &mut store::Library,
    access: &store::LabAccess,
    request: &Request,
) -> Result<ResponseOutcome, DispatchError> {
    match request {
        Request::Lookup {
            mutation_id,
            expected_intent,
        } => {
            let id = durable::MutationId::parse(mutation_id.as_str())
                .map_err(|_| DispatchError::Wire)?;
            let intent =
                store::Digest::parse(expected_intent.as_str()).map_err(|_| DispatchError::Wire)?;
            Ok(match db.lookup_reliable(access, &id, &intent) {
                Ok(Some(receipt)) => ResponseOutcome::Receipt(convert(&receipt)?),
                Ok(None) => ResponseOutcome::Unknown {
                    mutation_id: mutation_id.clone(),
                    intent: expected_intent.clone(),
                },
                Err(e) => storage_failure(mutation_id, e),
            })
        }
        Request::Apply { mutation } => {
            let id = durable::MutationId::parse(mutation.id.as_str())
                .map_err(|_| DispatchError::Wire)?;
            let digest = mutation.intent_digest().map_err(|_| DispatchError::Wire)?;
            let intent = store::Digest::parse(digest.as_str()).map_err(|_| DispatchError::Wire)?;
            // An already committed operation is historical receipt lookup, not
            // another write or a fresh scientific verdict. Raw payload integrity
            // was checked by request.validate before reaching this branch.
            match db.lookup_reliable(access, &id, &intent) {
                Ok(Some(receipt)) => return Ok(ResponseOutcome::Receipt(convert(&receipt)?)),
                Ok(None) => {}
                Err(e) => return Ok(storage_failure(&mutation.id, e)),
            }
            if let Some(payload) = mutation.operation.payload() {
                // Same automatic native/closure policy as ordinary server saves.
                // A request's WholeWorkspace enum alone cannot grant approval.
                if let Err(error) = host::document_publication_scope(access, payload) {
                    return Ok(match error {
                        host::Error::PublicationScope => failure(
                            library::ErrorCode::Forbidden,
                            "Native workspace is not the independently validated selected-document publication closure.",
                        ),
                        host::Error::Wire | host::Error::Native => failure(
                            library::ErrorCode::Invalid,
                            "Invalid native document/archive or exact selected event.",
                        ),
                        host::Error::Unavailable => failure(
                            library::ErrorCode::Unavailable,
                            "Native publication validation is unavailable.",
                        ),
                    });
                }
            }
            let native = to_storage_mutation(mutation)?;
            // This assertion binds the independently implemented WASM and native
            // intent contracts before any transactional effect is attempted.
            if native
                .intent_digest()
                .map_err(|_| DispatchError::Unavailable)?
                .as_str()
                != digest.as_str()
            {
                return Err(DispatchError::Unavailable);
            }
            Ok(match db.apply_reliable(access, &native) {
                Ok(result) => ResponseOutcome::Receipt(convert(&result.receipt)?),
                Err(e) => storage_failure(&mutation.id, e),
            })
        }
    }
}
/// Membership and lab selection are trusted host inputs, never wire fields.
/// Native publication-validation failures have no receipt/effect; committed
/// domain results (including CAS conflicts) are durably replayable exact receipts.
pub fn dispatch(
    db: &mut store::Library,
    access: &store::LabAccess,
    request: &RequestEnvelope,
) -> Result<ResponseEnvelope, DispatchError> {
    if request.expected_lab.as_str() != access.lab_id().as_str() {
        return Err(DispatchError::WrongLab);
    }
    request.validate().map_err(|_| DispatchError::Wire)?;
    let lab = db.lab(access).map_err(|_| DispatchError::Unavailable)?;
    let context = LabContext {
        lab_id: library::LabId::new(lab.id.as_str()).map_err(|_| DispatchError::Unavailable)?,
        name: Text::new(lab.name).map_err(|_| DispatchError::Unavailable)?,
    };
    let response = ResponseEnvelope {
        request_id: request.request_id,
        context: Some(context),
        outcome: execute(db, access, &request.request)?,
    };
    response
        .validate_for(request)
        .map_err(|_| DispatchError::Unavailable)?;
    Ok(response)
}
