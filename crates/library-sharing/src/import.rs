//! Atomic native graph import. Never rewrites events, flattens saved products,
//! treats a remote locator as authority, or replaces unrelated local documents.
use crate::{Selection, inspect_document_closure};
use collagene_core::documents::{
    transfer::Candidate,
    workspace::{self, RevisionEvent, Workspace, WorkspaceAttachment, WorkspaceDocument},
};
use collagene_model::Molecule;
use serde::{Deserialize, Serialize};
use std::{
    collections::{BTreeMap, BTreeSet},
    io::Write,
};

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
#[serde(deny_unknown_fields)]
pub enum SelectedPolicy {
    /// A differing selected head/event is a conflict, including same DNA under
    /// a different event. No automatic sequence-based branch resolution.
    Strict,
    /// Explicitly requested update with exact local pin and proven event ancestry.
    FastForward { expected_local: Selection },
    /// Explicit branch switch. The returned displaced bookmark must be durably
    /// retained by the caller BEFORE committing the resulting workspace.
    Checkout { expected_local: Selection },
}
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
#[serde(deny_unknown_fields)]
pub struct ImportOptions {
    pub selected_policy: SelectedPolicy,
    pub activate_selected: bool,
}
impl Default for ImportOptions {
    fn default() -> Self {
        Self {
            selected_policy: SelectedPolicy::Strict,
            activate_selected: true,
        }
    }
}
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct HeadBookmark {
    /// The complete old navigation branch, not just a sequence revision.
    pub document: WorkspaceDocument,
    pub imported_selection: Selection,
}
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ImportReport {
    /// A proposed replacement state. Explicit checkout/fast-forward bookmarks
    /// must be persisted before the caller commits this returned value.
    pub workspace: Workspace,
    pub selected: Selection,
    pub inserted_documents: Vec<String>,
    pub preserved_dependency_heads: Vec<Selection>,
    pub displaced_selected: Option<HeadBookmark>,
}
#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq)]
pub enum IdentityKind {
    Revision,
    Event,
    Candidate,
    Attachment,
}
#[derive(Debug, thiserror::Error)]
pub enum ImportError {
    #[error("the local workspace is invalid")]
    InvalidLocal,
    #[error("the fetched native workspace is invalid")]
    InvalidFetched,
    #[error("the fetched workspace is not the exact supported selected-document closure")]
    InvalidClosure,
    #[error("immutable {kind:?} identity collision: {id}")]
    Collision { kind: IdentityKind, id: String },
    #[error("the selected native document has a different local head/event")]
    SelectedConflict {
        local: Selection,
        incoming: Selection,
    },
    #[error("the local selected head/event changed since the explicit import choice")]
    StaleLocal {
        expected: Selection,
        actual: Option<Selection>,
    },
    #[error("fast-forward requires exact event ancestry; sequence equality is insufficient")]
    NotFastForward {
        local: Selection,
        incoming: Selection,
    },
    #[error("the combined native workspace exceeds its bounded dimensions or size")]
    TooLarge,
    #[error("the combined native workspace is inconsistent; no state was changed")]
    InvalidMerged,
}
fn pin(d: &WorkspaceDocument) -> Selection {
    Selection {
        document_id: d.id.clone(),
        revision: d.head.clone(),
        head_event: d.head_event.clone(),
    }
}
fn event_equal(a: &RevisionEvent, b: &RevisionEvent) -> bool {
    a.id == b.id
        && a.document_id == b.document_id
        && a.revision == b.revision
        && a.parents == b.parents
        && a.operation == b.operation
}
fn attachment_equal(a: &WorkspaceAttachment, b: &WorkspaceAttachment) -> bool {
    a.kind == b.kind && a.name == b.name && a.media_type == b.media_type && a.bytes == b.bytes
}
struct Count {
    used: usize,
    cap: usize,
}
impl Write for Count {
    fn write(&mut self, b: &[u8]) -> std::io::Result<usize> {
        self.used = self
            .used
            .checked_add(b.len())
            .filter(|n| *n <= self.cap)
            .ok_or_else(|| std::io::Error::other("bounded native import"))?;
        Ok(b.len())
    }
    fn flush(&mut self) -> std::io::Result<()> {
        Ok(())
    }
}
fn size(v: &impl Serialize, cap: usize) -> Result<usize, ImportError> {
    let mut out = Count { used: 0, cap };
    serde_json::to_writer(&mut out, v).map_err(|_| ImportError::TooLarge)?;
    Ok(out.used)
}
struct Bytes(Vec<u8>);
impl Write for Bytes {
    fn write(&mut self, b: &[u8]) -> std::io::Result<usize> {
        if self
            .0
            .len()
            .checked_add(b.len())
            .is_none_or(|n| n > workspace::MAX_ENTRY_BYTES)
        {
            return Err(std::io::Error::other("bounded import candidate"));
        }
        self.0.extend_from_slice(b);
        Ok(b.len())
    }
    fn flush(&mut self) -> std::io::Result<()> {
        Ok(())
    }
}
struct Compare<'a> {
    bytes: &'a [u8],
    offset: usize,
}
impl Write for Compare<'_> {
    fn write(&mut self, b: &[u8]) -> std::io::Result<usize> {
        let end = self
            .offset
            .checked_add(b.len())
            .ok_or_else(|| std::io::Error::other("candidate length"))?;
        if self.bytes.get(self.offset..end) != Some(b) {
            return Err(std::io::Error::other("candidate differs"));
        }
        self.offset = end;
        Ok(b.len())
    }
    fn flush(&mut self) -> std::io::Result<()> {
        Ok(())
    }
}
fn candidate_equal(a: &Candidate, b: &Candidate) -> Result<bool, ImportError> {
    let mut first = Bytes(Vec::new());
    serde_json::to_writer(&mut first, a).map_err(|_| ImportError::TooLarge)?;
    let mut other = Compare {
        bytes: &first.0,
        offset: 0,
    };
    Ok(serde_json::to_writer(&mut other, b).is_ok() && other.offset == first.0.len())
}
fn exact_closure(
    fetched: &Workspace,
    selected: &Selection,
    private_local: bool,
) -> Result<(), ImportError> {
    let inventory = if private_local {
        crate::inspect_private_document_closure(fetched, selected)
    } else {
        inspect_document_closure(fetched, selected)
    }
    .map_err(|_| ImportError::InvalidClosure)?;
    if fetched.active_document != selected.document_id
        || fetched.active_candidate.is_some()
        || inventory.documents.len() != fetched.documents.len()
        || inventory.revisions.len() != fetched.revisions.len()
        || inventory.events.len() != fetched.history.len()
        || inventory.candidates.len() != fetched.candidates.len()
        || inventory.attachments.len() != fetched.attachments.len()
        || fetched.documents.iter().any(|d| {
            inventory
                .document_heads
                .get(&d.id)
                .is_none_or(|h| h.revision != d.head || h.head_event != d.head_event)
                || !d.undo.is_empty()
                || !d.redo.is_empty()
                || !d.undo_events.is_empty()
                || !d.redo_events.is_empty()
        })
        || fetched
            .revisions
            .keys()
            .any(|id| !inventory.revisions.contains(id))
        || fetched
            .history
            .iter()
            .any(|e| !inventory.events.contains(&e.id))
        || fetched
            .candidates
            .iter()
            .any(|c| !inventory.candidates.contains(&c.plan_revision))
        || fetched
            .attachments
            .keys()
            .any(|id| !inventory.attachments.contains(id))
    {
        return Err(ImportError::InvalidClosure);
    }
    Ok(())
}
fn ancestor(events: &BTreeMap<&str, &RevisionEvent>, old: &Selection, new: &Selection) -> bool {
    let (Some(old), Some(new)) = (old.head_event.as_deref(), new.head_event.as_deref()) else {
        return false;
    };
    let mut visited = BTreeSet::new();
    let mut pending = vec![new];
    while let Some(id) = pending.pop() {
        if id == old {
            return true;
        }
        if !visited.insert(id) {
            continue;
        }
        let Some(event) = events.get(id) else {
            return false;
        };
        pending.extend(event.parents.iter().map(String::as_str));
    }
    false
}
/// Merge one fetched native publication closure into a validated local workspace.
/// The two inputs remain untouched on success and every error. This is a native
/// graph union, not remote sync ancestry negotiation or scientific verification.
pub fn merge_document_closure(
    local: &Workspace,
    fetched: &Workspace,
    selected: &Selection,
    options: &ImportOptions,
) -> Result<ImportReport, ImportError> {
    merge_inner(local, fetched, selected, options, false)
}
/// Trusted native-local counterpart. Scope validation allows exact private SBOL
/// originals, but every immutable collision, head CAS and bookmark check is the
/// same as public import. Never call this with an HTTP-origin document.
pub fn merge_private_document_closure(
    local: &Workspace,
    fetched: &Workspace,
    selected: &Selection,
    options: &ImportOptions,
) -> Result<ImportReport, ImportError> {
    merge_inner(local, fetched, selected, options, true)
}
fn merge_inner(
    local: &Workspace,
    fetched: &Workspace,
    selected: &Selection,
    options: &ImportOptions,
    private_local: bool,
) -> Result<ImportReport, ImportError> {
    local.validate().map_err(|_| ImportError::InvalidLocal)?;
    fetched
        .validate()
        .map_err(|_| ImportError::InvalidFetched)?;
    exact_closure(fetched, selected, private_local)?;
    let incoming = fetched
        .documents
        .iter()
        .find(|d| d.id == selected.document_id)
        .ok_or(ImportError::InvalidClosure)?;
    if pin(incoming) != *selected {
        return Err(ImportError::InvalidClosure);
    }
    // Borrowed identity collision checks precede every potentially large clone.
    let mut revisions: BTreeMap<&str, &Molecule> = local
        .revisions
        .iter()
        .map(|(id, m)| (id.as_str(), m))
        .collect();
    for (id, m) in &fetched.revisions {
        if let Some(old) = revisions.get(id.as_str()) {
            if *old != m {
                return Err(ImportError::Collision {
                    kind: IdentityKind::Revision,
                    id: id.clone(),
                });
            }
        } else {
            revisions.insert(id, m);
        }
    }
    let mut events: BTreeMap<&str, &RevisionEvent> =
        local.history.iter().map(|e| (e.id.as_str(), e)).collect();
    for event in &fetched.history {
        if let Some(old) = events.get(event.id.as_str()) {
            if !event_equal(old, event) {
                return Err(ImportError::Collision {
                    kind: IdentityKind::Event,
                    id: event.id.clone(),
                });
            }
        } else {
            events.insert(&event.id, event);
        }
    }
    let mut candidates: BTreeMap<&str, &Candidate> = local
        .candidates
        .iter()
        .map(|c| (c.plan_revision.as_str(), c))
        .collect();
    for candidate in &fetched.candidates {
        if let Some(old) = candidates.get(candidate.plan_revision.as_str()) {
            if !candidate_equal(old, candidate)? {
                return Err(ImportError::Collision {
                    kind: IdentityKind::Candidate,
                    id: candidate.plan_revision.clone(),
                });
            }
        } else {
            candidates.insert(&candidate.plan_revision, candidate);
        }
    }
    let mut attachments: BTreeMap<&str, &WorkspaceAttachment> = local
        .attachments
        .iter()
        .map(|(id, a)| (id.as_str(), a))
        .collect();
    for (id, a) in &fetched.attachments {
        if let Some(old) = attachments.get(id.as_str()) {
            if !attachment_equal(old, a) {
                return Err(ImportError::Collision {
                    kind: IdentityKind::Attachment,
                    id: id.clone(),
                });
            }
        } else {
            attachments.insert(id, a);
        }
    }
    let existing = local
        .documents
        .iter()
        .find(|d| d.id == selected.document_id);
    if let SelectedPolicy::FastForward { expected_local }
    | SelectedPolicy::Checkout { expected_local } = &options.selected_policy
    {
        if existing.map(pin).as_ref() != Some(expected_local)
            || expected_local.document_id != selected.document_id
        {
            return Err(ImportError::StaleLocal {
                expected: expected_local.clone(),
                actual: existing.map(pin),
            });
        }
    }
    let mut displaced_selected = None;
    if let Some(old) = existing {
        let old_pin = pin(old);
        if old_pin != *selected {
            match &options.selected_policy {
                SelectedPolicy::Strict => {
                    return Err(ImportError::SelectedConflict {
                        local: old_pin,
                        incoming: selected.clone(),
                    });
                }
                SelectedPolicy::FastForward { .. } => {
                    if !ancestor(&events, &old_pin, selected) {
                        return Err(ImportError::NotFastForward {
                            local: old_pin,
                            incoming: selected.clone(),
                        });
                    }
                }
                SelectedPolicy::Checkout { .. } => {}
            }
            displaced_selected = Some(HeadBookmark {
                document: old.clone(),
                imported_selection: selected.clone(),
            });
        }
    }
    let local_ids = local
        .documents
        .iter()
        .map(|d| d.id.as_str())
        .collect::<BTreeSet<_>>();
    let inserted_documents = fetched
        .documents
        .iter()
        .filter(|d| !local_ids.contains(d.id.as_str()))
        .map(|d| d.id.clone())
        .collect::<Vec<_>>();
    if local
        .documents
        .len()
        .checked_add(inserted_documents.len())
        .is_none_or(|n| n > workspace::MAX_DOCUMENTS)
        || revisions.len() > workspace::MAX_REVISIONS
        || events.len() > workspace::MAX_HISTORY_EVENTS
        || candidates.len() > workspace::MAX_CANDIDATES
        || attachments.len() > workspace::MAX_ATTACHMENTS
    {
        return Err(ImportError::TooLarge);
    }
    // Conservative aggregate accounting before cloning large models. Native
    // final validation remains authoritative for all model/history/run bounds.
    let mut charged = 0usize;
    for m in revisions.values() {
        charged = charged
            .checked_add(size(m, workspace::MAX_ENTRY_BYTES)?)
            .ok_or(ImportError::TooLarge)?;
        if charged > workspace::MAX_EXPANDED_BYTES {
            return Err(ImportError::TooLarge);
        }
    }
    for c in candidates.values() {
        charged = charged
            .checked_add(size(c, workspace::MAX_ENTRY_BYTES)?)
            .ok_or(ImportError::TooLarge)?;
        if charged > workspace::MAX_EXPANDED_BYTES {
            return Err(ImportError::TooLarge);
        }
    }
    for a in attachments.values() {
        charged = charged
            .checked_add(a.bytes.len())
            .ok_or(ImportError::TooLarge)?;
        if charged > workspace::MAX_EXPANDED_BYTES {
            return Err(ImportError::TooLarge);
        }
    }
    let mut documents = local.documents.clone();
    let mut preserved_dependency_heads = Vec::new();
    for imported in &fetched.documents {
        if let Some(at) = documents.iter().position(|d| d.id == imported.id) {
            if imported.id == selected.document_id {
                if displaced_selected.is_some() {
                    let d = &mut documents[at];
                    d.head = selected.revision.clone();
                    d.head_event = selected.head_event.clone();
                    d.undo.clear();
                    d.redo.clear();
                    d.undo_events.clear();
                    d.redo_events.clear();
                }
            } else {
                preserved_dependency_heads.push(pin(&documents[at]));
            }
        } else {
            documents.push(imported.clone());
        }
    }
    let active_document = if options.activate_selected {
        selected.document_id.clone()
    } else {
        local.active_document.clone()
    };
    size(
        &(
            &local.name,
            &active_document,
            &documents,
            events.values().collect::<Vec<_>>(),
            &local.active_candidate,
        ),
        workspace::MAX_MANIFEST_BYTES,
    )?;
    let local_event_ids = local
        .history
        .iter()
        .map(|e| e.id.as_str())
        .collect::<BTreeSet<_>>();
    let local_candidate_ids = local
        .candidates
        .iter()
        .map(|c| c.plan_revision.as_str())
        .collect::<BTreeSet<_>>();
    let result = Workspace {
        name: local.name.clone(),
        active_document,
        documents,
        revisions: revisions
            .into_iter()
            .map(|(id, m)| (id.into(), m.clone()))
            .collect(),
        // Preserve local ordering; fetched nonduplicates follow their archive order.
        history: local
            .history
            .iter()
            .chain(
                fetched
                    .history
                    .iter()
                    .filter(|e| !local_event_ids.contains(e.id.as_str())),
            )
            .cloned()
            .collect(),
        candidates: local
            .candidates
            .iter()
            .chain(
                fetched
                    .candidates
                    .iter()
                    .filter(|c| !local_candidate_ids.contains(c.plan_revision.as_str())),
            )
            .cloned()
            .collect(),
        active_candidate: local.active_candidate.clone(),
        attachments: attachments
            .into_iter()
            .map(|(id, a)| (id.into(), a.clone()))
            .collect(),
    };
    result.validate().map_err(|_| ImportError::InvalidMerged)?;
    Ok(ImportReport {
        workspace: result,
        selected: selected.clone(),
        inserted_documents,
        preserved_dependency_heads,
        displaced_selected,
    })
}
