//! Immutable revision events with reversible navigation and retained branches.
//! Navigation never establishes scientific verification or executes event descriptions.
use crate::{
    Error, Molecule,
    workspace::{RevisionEvent, Workspace, WorkspaceDocument},
};
use std::{
    collections::{BTreeMap, BTreeSet},
    io,
    sync::Arc,
};
pub const MAX_UNDO_STATES: usize = 100;
pub const MAX_HISTORY_EVENTS: usize = 10_000;
pub const MAX_HISTORY_BYTES: usize = 64 * 1024 * 1024;
#[derive(Clone)]
struct Node {
    event: Arc<RevisionEvent>,
    molecule: Arc<Molecule>,
    bytes: usize,
}
#[derive(Debug, Clone)]
pub struct HistorySummary {
    pub event_id: String,
    pub revision: String,
    pub operation: String,
    pub parents: Vec<String>,
    pub current: bool,
    pub sequence_length: usize,
}
#[derive(Clone)]
pub struct DocumentHistory {
    pub id: String,
    pub molecule: Molecule,
    nodes: Vec<Node>,
    head: String,
    undo: Vec<String>,
    redo: Vec<String>,
}
struct Count(usize);
impl io::Write for Count {
    fn write(&mut self, bytes: &[u8]) -> io::Result<usize> {
        self.0 = self
            .0
            .checked_add(bytes.len())
            .filter(|n| *n <= MAX_HISTORY_BYTES)
            .ok_or_else(|| io::Error::other("history capacity"))?;
        Ok(bytes.len())
    }
    fn flush(&mut self) -> io::Result<()> {
        Ok(())
    }
}
fn size(m: &Molecule) -> Result<usize, Error> {
    let mut count = Count(0);
    serde_json::to_writer(&mut count, m).map_err(|_| Error::TooLarge)?;
    Ok(count.0)
}
impl DocumentHistory {
    pub fn new(id: impl Into<String>, molecule: Molecule) -> Result<Self, Error> {
        molecule.validate()?;
        let id = id.into();
        let revision = molecule.revision();
        let bytes = size(&molecule)?;
        let event = RevisionEvent::new(id.clone(), revision, vec![], "Opened document")?;
        Ok(Self {
            id,
            molecule: molecule.clone(),
            head: event.id.clone(),
            nodes: vec![Node {
                event: Arc::new(event),
                molecule: Arc::new(molecule),
                bytes,
            }],
            undo: vec![],
            redo: vec![],
        })
    }
    pub fn can_undo(&self) -> bool {
        !self.undo.is_empty()
    }
    pub fn can_redo(&self) -> bool {
        !self.redo.is_empty()
    }
    pub fn event_count(&self) -> usize {
        self.nodes.len()
    }
    pub fn head_event(&self) -> &str {
        &self.head
    }
    /// Exact reversible navigation event identities, oldest entry first and next
    /// operation last. These stacks describe navigation, not biological ancestry.
    pub fn navigation_event_ids(&self) -> (Vec<String>, Vec<String>) {
        (self.undo.clone(), self.redo.clone())
    }
    /// Open an immutable retained event as a new document while keeping its exact
    /// source-event parent. The caller must retain the source event in the workspace.
    pub fn fork_event(&self, event: &str, id: impl Into<String>) -> Result<Self, Error> {
        let source = self.node(event)?;
        let id = id.into();
        if id == self.id {
            return Err(Error::Format(
                "a fork requires a distinct document identity".into(),
            ));
        }
        let fork = RevisionEvent::new(
            id.clone(),
            source.event.revision.clone(),
            vec![source.event.id.clone()],
            "Opened retained revision as new document",
        )?;
        Ok(Self {
            id,
            molecule: source.molecule.as_ref().clone(),
            head: fork.id.clone(),
            nodes: vec![Node {
                event: Arc::new(fork),
                molecule: source.molecule.clone(),
                bytes: source.bytes,
            }],
            undo: vec![],
            redo: vec![],
        })
    }
    /// Retain a transformed document with the exact selected source event as its
    /// parent. This records computational lineage only; operation text and parent
    /// links never establish physical chemistry or scientific verification.
    /// The caller must retain this source event in the containing workspace.
    pub fn derive_event(
        &self,
        event: &str,
        id: impl Into<String>,
        molecule: Molecule,
        operation: impl Into<String>,
    ) -> Result<Self, Error> {
        let source = self.node(event)?.event.id.clone();
        let id = id.into();
        if id == self.id {
            return Err(Error::Format(
                "a derived document requires a distinct document identity".into(),
            ));
        }
        let mut derived = Self::new(id, molecule)?;
        let node = RevisionEvent::new(
            derived.id.clone(),
            derived.molecule.revision(),
            vec![source],
            operation,
        )?;
        derived.head = node.id.clone();
        derived.nodes[0].event = Arc::new(node);
        Ok(derived)
    }
    /// Retain a computational product with every explicitly bound source event.
    /// Sources remain immutable and must be retained in the same workspace.
    /// Parent edges and operation text never establish physical or verification authority.
    pub fn derive_from_events(
        sources: &[(&Self, &str)],
        id: impl Into<String>,
        molecule: Molecule,
        operation: impl Into<String>,
    ) -> Result<Self, Error> {
        if sources.is_empty() {
            return Err(Error::Format(
                "a derived document requires retained source events".into(),
            ));
        }
        if sources.len() > 128 {
            return Err(Error::TooLarge);
        }
        let id = id.into();
        let mut seen = BTreeSet::new();
        let mut parents = Vec::new();
        for (source, event) in sources {
            if source.id == id {
                return Err(Error::Format(
                    "a derived document requires a distinct identity from every source".into(),
                ));
            }
            let node = source.node(event)?;
            if seen.insert(node.event.id.clone()) {
                parents.push(node.event.id.clone());
            }
        }
        let mut derived = Self::new(id, molecule)?;
        let event = RevisionEvent::new(
            derived.id.clone(),
            derived.molecule.revision(),
            parents,
            operation,
        )?;
        derived.head = event.id.clone();
        derived.nodes[0].event = Arc::new(event);
        Ok(derived)
    }
    pub fn revision_molecule(&self, event: &str) -> Result<Molecule, Error> {
        Ok(self.node(event)?.molecule.as_ref().clone())
    }
    fn node(&self, id: &str) -> Result<&Node, Error> {
        self.nodes
            .iter()
            .find(|n| n.event.id == id)
            .ok_or_else(|| Error::Format("history event is missing".into()))
    }
    pub fn revision_sizes(&self) -> BTreeMap<String, usize> {
        self.nodes
            .iter()
            .map(|n| (n.event.revision.clone(), n.bytes))
            .collect()
    }
    pub fn retained_bytes(&self) -> usize {
        self.revision_sizes().values().sum()
    }
    /// Atomic commit: stale current state, invalid data or capacity errors leave every
    /// navigation stack and event unchanged. Redo branches stay in immutable history.
    pub fn commit(&mut self, next: Molecule, operation: impl Into<String>) -> Result<(), Error> {
        next.validate()?;
        if self.node(&self.head)?.event.revision != self.molecule.revision() {
            return Err(Error::Format(
                "current state was changed outside its history controller".into(),
            ));
        }
        let revision = next.revision();
        if revision == self.node(&self.head)?.event.revision {
            return Ok(());
        }
        let event = RevisionEvent::new(
            self.id.clone(),
            revision.clone(),
            vec![self.head.clone()],
            operation,
        )?;
        let existing = self.nodes.iter().find(|n| n.event.id == event.id).cloned();
        let node = if let Some(existing) = existing {
            existing
        } else {
            if self.nodes.len() >= MAX_HISTORY_EVENTS {
                return Err(Error::Unsupported("Document history capacity reached; export the workspace or open a revision as a new document before continuing".into()));
            }
            let shared = self.nodes.iter().find(|n| n.event.revision == revision);
            let (molecule, bytes) = if let Some(shared) = shared {
                (shared.molecule.clone(), shared.bytes)
            } else {
                let bytes = size(&next)?;
                if self.retained_bytes().saturating_add(bytes) > MAX_HISTORY_BYTES {
                    return Err(Error::Unsupported("Retained document history exceeds 64 MiB; export the workspace or fork a new document before continuing".into()));
                }
                (Arc::new(next.clone()), bytes)
            };
            Node {
                event: Arc::new(event),
                molecule,
                bytes,
            }
        };
        if !self.nodes.iter().any(|n| n.event.id == node.event.id) {
            self.nodes.push(node.clone());
        }
        self.undo.push(self.head.clone());
        if self.undo.len() > MAX_UNDO_STATES {
            self.undo.drain(..self.undo.len() - MAX_UNDO_STATES);
        }
        self.redo.clear();
        self.head = node.event.id.clone();
        self.molecule = next;
        Ok(())
    }
    pub fn undo(&mut self) -> Result<bool, Error> {
        let Some(id) = self.undo.last().cloned() else {
            return Ok(false);
        };
        let molecule = self.node(&id)?.molecule.as_ref().clone();
        self.undo.pop();
        self.redo.push(self.head.clone());
        self.head = id;
        self.molecule = molecule;
        Ok(true)
    }
    pub fn redo(&mut self) -> Result<bool, Error> {
        let Some(id) = self.redo.last().cloned() else {
            return Ok(false);
        };
        let molecule = self.node(&id)?.molecule.as_ref().clone();
        self.redo.pop();
        self.undo.push(self.head.clone());
        self.head = id;
        self.molecule = molecule;
        Ok(true)
    }
    pub fn checkout(&mut self, event: &str) -> Result<(), Error> {
        let molecule = self.node(event)?.molecule.as_ref().clone();
        if event == self.head {
            return Ok(());
        }
        self.undo.push(self.head.clone());
        if self.undo.len() > MAX_UNDO_STATES {
            self.undo.drain(..self.undo.len() - MAX_UNDO_STATES);
        }
        self.redo.clear();
        self.head = event.into();
        self.molecule = molecule;
        Ok(())
    }
    pub fn summaries(&self) -> Vec<HistorySummary> {
        self.nodes
            .iter()
            .map(|n| HistorySummary {
                event_id: n.event.id.clone(),
                revision: n.event.revision.clone(),
                operation: n.event.operation.clone(),
                parents: n.event.parents.clone(),
                current: n.event.id == self.head,
                sequence_length: n.molecule.sequence.len(),
            })
            .collect()
    }
    pub fn contribute(&self, workspace: &mut Workspace) -> Result<(), Error> {
        if workspace.documents.iter().any(|d| d.id == self.id) {
            return Err(Error::Format("duplicate document identity".into()));
        }
        if self.node(&self.head)?.event.revision != self.molecule.revision() {
            return Err(Error::Format(
                "current state does not match history head".into(),
            ));
        }
        for n in &self.nodes {
            workspace
                .revisions
                .entry(n.event.revision.clone())
                .or_insert_with(|| n.molecule.as_ref().clone());
            workspace.history.push(n.event.as_ref().clone());
        }
        workspace.documents.push(WorkspaceDocument {
            id: self.id.clone(),
            head: self.molecule.revision(),
            head_event: Some(self.head.clone()),
            undo: self
                .undo
                .iter()
                .map(|id| self.node(id).map(|n| n.event.revision.clone()))
                .collect::<Result<_, _>>()?,
            redo: self
                .redo
                .iter()
                .map(|id| self.node(id).map(|n| n.event.revision.clone()))
                .collect::<Result<_, _>>()?,
            undo_events: self.undo.clone(),
            redo_events: self.redo.clone(),
        });
        Ok(())
    }
    /// Preserve all events, including non-head branches and cross-document parents.
    /// Legacy stacks without event references acquire explicit recovery roots rather
    /// than guessed edit ancestry. Original recorded events remain retained.
    pub fn restore_all(workspace: &Workspace) -> Result<Vec<Self>, Error> {
        workspace.validate()?;
        let shared = workspace
            .revisions
            .iter()
            .map(|(id, m)| Ok((id.clone(), (Arc::new(m.clone()), size(m)?))))
            .collect::<Result<BTreeMap<_, _>, Error>>()?;
        let mut output = Vec::new();
        for doc in &workspace.documents {
            let mut nodes = workspace
                .history
                .iter()
                .filter(|e| e.document_id == doc.id)
                .map(|event| {
                    let (molecule, bytes) = &shared[&event.revision];
                    Node {
                        event: Arc::new(event.clone()),
                        molecule: molecule.clone(),
                        bytes: *bytes,
                    }
                })
                .collect::<Vec<_>>();
            let mut recovered = |revision: &str| -> Result<String, Error> {
                let event = RevisionEvent::new(
                    doc.id.clone(),
                    revision.to_owned(),
                    vec![],
                    "Recovered navigation snapshot; original event was not recorded",
                )?;
                let id = event.id.clone();
                if !nodes.iter().any(|n| n.event.id == id) {
                    let (molecule, bytes) = &shared[revision];
                    nodes.push(Node {
                        event: Arc::new(event),
                        molecule: molecule.clone(),
                        bytes: *bytes,
                    });
                }
                Ok(id)
            };
            let head = if let Some(head) = &doc.head_event {
                head.clone()
            } else {
                recovered(&doc.head)?
            };
            let undo = if !doc.undo_events.is_empty() {
                doc.undo_events.clone()
            } else {
                doc.undo
                    .iter()
                    .map(|r| recovered(r))
                    .collect::<Result<_, _>>()?
            };
            let redo = if !doc.redo_events.is_empty() {
                doc.redo_events.clone()
            } else {
                doc.redo
                    .iter()
                    .map(|r| recovered(r))
                    .collect::<Result<_, _>>()?
            };
            output.push(Self {
                id: doc.id.clone(),
                molecule: shared[&doc.head].0.as_ref().clone(),
                nodes,
                head,
                undo,
                redo,
            });
        }
        // An archive may intentionally hold more history than the current edit cap;
        // restoring it preserves evidence instead of silently trimming it.
        Ok(output)
    }
    pub fn retained_revision_ids(&self) -> BTreeSet<String> {
        self.nodes
            .iter()
            .map(|n| n.event.revision.clone())
            .collect()
    }
}
