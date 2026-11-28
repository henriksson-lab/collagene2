//! Immutable revision events with reversible navigation and retained branches.
//! Navigation never establishes scientific verification or executes event descriptions.
use crate::documents::workspace::{RevisionEvent, Workspace, WorkspaceDocument};
use collagene_model::{Error, Molecule};
use serde::Serialize;
use std::{
    collections::{BTreeMap, BTreeSet},
    io,
    sync::Arc,
};
pub mod transition;
use crate::{
    documents::reference_edit::Operation as ReferenceEdit, material::lineage::EventBinding,
};
use transition::{AnalysisOrigin, Transition};
pub const MAX_UNDO_STATES: usize = 100;
pub const MAX_HISTORY_EVENTS: usize = 10_000;
pub const MAX_HISTORY_BYTES: usize = 64 * 1024 * 1024;
#[derive(Clone)]
struct Node {
    event: Arc<RevisionEvent>,
    molecule: Arc<Molecule>,
    bytes: usize,
    event_bytes: usize,
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
fn size(m: &impl Serialize) -> Result<usize, Error> {
    let mut count = Count(0);
    serde_json::to_writer(&mut count, m).map_err(|_| Error::TooLarge)?;
    Ok(count.0)
}
impl DocumentHistory {
    pub fn new(id: impl Into<String>, molecule: Molecule) -> Result<Self, Error> {
        collagene_model::bounds::native_frame_preflight(
            std::iter::once(&molecule),
            crate::documents::reference_edit::MAX_EDIT_BYTES,
        )?;
        molecule.validate()?;
        let id = id.into();
        let revision = molecule.revision();
        let bytes = size(&molecule)?;
        let event = RevisionEvent::new(
            id.clone(),
            revision,
            vec![],
            "Opened document",
            Transition::IndependentReference,
        )?;
        let event_bytes = size(&event)?;
        Ok(Self {
            id,
            molecule: molecule.clone(),
            head: event.id.clone(),
            nodes: vec![Node {
                event: Arc::new(event),
                molecule: Arc::new(molecule),
                bytes,
                event_bytes,
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
    /// Revision of the immutable current event, without cloning history summaries.
    pub fn head_revision(&self) -> &str {
        &self
            .node(&self.head)
            .expect("history head is a retained event")
            .event
            .revision
    }
    pub fn binding(&self) -> EventBinding {
        EventBinding {
            document_id: self.id.clone(),
            event_id: self.head.clone(),
            revision: self.head_revision().into(),
        }
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
            Transition::Fork {
                source: EventBinding {
                    document_id: self.id.clone(),
                    event_id: source.event.id.clone(),
                    revision: source.event.revision.clone(),
                },
            },
        )?;
        let event_bytes = size(&fork)?;
        Ok(Self {
            id,
            molecule: source.molecule.as_ref().clone(),
            head: fork.id.clone(),
            nodes: vec![Node {
                event: Arc::new(fork),
                molecule: source.molecule.clone(),
                bytes: source.bytes,
                event_bytes,
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
        origin: AnalysisOrigin,
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
            Transition::Analysis { origin },
        )?;
        derived.nodes[0].event_bytes = size(&node)?;
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
        origin: AnalysisOrigin,
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
            Transition::Analysis { origin },
        )?;
        derived.nodes[0].event_bytes = size(&event)?;
        derived.head = event.id.clone();
        derived.nodes[0].event = Arc::new(event);
        Ok(derived)
    }
    pub fn revision_molecule(&self, event: &str) -> Result<Molecule, Error> {
        Ok(self.node(event)?.molecule.as_ref().clone())
    }
    /// Borrow immutable retained sequence data without cloning it. Reading a
    /// history event does not restore a scientific assessment or material state.
    pub fn revision_molecule_ref(&self, event: &str) -> Result<&Molecule, Error> {
        Ok(self.node(event)?.molecule.as_ref())
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
        self.revision_sizes()
            .values()
            .copied()
            .chain(self.nodes.iter().map(|node| node.event_bytes))
            .fold(0, usize::saturating_add)
    }
    /// Atomic commit: stale current state, invalid data or capacity errors leave every
    /// navigation stack and event unchanged. Redo branches stay in immutable history.
    pub fn commit(
        &mut self,
        actions: Vec<ReferenceEdit>,
        operation: impl Into<String>,
    ) -> Result<(), Error> {
        self.commit_with_source_events(actions, operation, &[])
    }
    /// Commit a reversible reference edit with explicitly retained source events.
    /// The host must validate cross-document event ownership in its staged workspace.
    /// Parent edges record computational ancestry, never physical verification.
    /// Unlike an ordinary no-op edit, explicit source ancestry records an event
    /// even when the resulting molecular content is unchanged.
    pub fn commit_with_source_events(
        &mut self,
        actions: Vec<ReferenceEdit>,
        operation: impl Into<String>,
        source_events: &[String],
    ) -> Result<(), Error> {
        if source_events.len() > 128 {
            return Err(Error::TooLarge);
        }
        if source_events
            .iter()
            .any(|event| !crate::documents::workspace::event_id(event))
        {
            return Err(Error::Format(
                "invalid supplemental source event identity".into(),
            ));
        }
        self.ensure_current()?;
        let source = self.binding();
        let next = ReferenceEdit::apply_all(&self.molecule, &actions)?;
        let supplemental_sources =
            source_events
                .iter()
                .filter(|id| **id != self.head)
                .fold(Vec::new(), |mut ids, id| {
                    if !ids.contains(id) {
                        ids.push(id.clone());
                    }
                    ids
                });
        let transition = Transition::ReferenceEdit {
            source,
            actions,
            supplemental_sources,
        };
        self.commit_result(next, operation, source_events, transition)
    }
    /// Retain an actual named computational producer, without claiming a
    /// reference edit or treating its descriptive label as a replay recipe.
    pub fn commit_analysis(
        &mut self,
        next: Molecule,
        operation: impl Into<String>,
        origin: AnalysisOrigin,
        source_events: &[String],
    ) -> Result<(), Error> {
        self.commit_result(
            next,
            operation,
            source_events,
            Transition::Analysis { origin },
        )
    }
    fn ensure_current(&self) -> Result<(), Error> {
        collagene_model::bounds::native_frame_preflight(
            std::iter::once(&self.molecule),
            crate::documents::reference_edit::MAX_EDIT_BYTES,
        )?;
        if self.node(&self.head)?.event.revision != self.molecule.revision() {
            return Err(Error::Format(
                "current state was changed outside its history controller".into(),
            ));
        }
        Ok(())
    }
    fn commit_result(
        &mut self,
        next: Molecule,
        operation: impl Into<String>,
        source_events: &[String],
        transition: Transition,
    ) -> Result<(), Error> {
        if source_events.len() > 128 {
            return Err(Error::TooLarge);
        }
        if source_events
            .iter()
            .any(|event| !crate::documents::workspace::event_id(event))
        {
            return Err(Error::Format(
                "invalid supplemental source event identity".into(),
            ));
        }
        self.ensure_current()?;
        collagene_model::bounds::native_frame_preflight(
            std::iter::once(&next),
            crate::documents::reference_edit::MAX_EDIT_BYTES,
        )?;
        next.validate()?;
        let revision = next.revision();
        if revision == self.node(&self.head)?.event.revision && source_events.is_empty() {
            return Ok(());
        }
        let mut parents = vec![self.head.clone()];
        for source in source_events {
            if !parents.contains(source) {
                parents.push(source.clone());
            }
        }
        let event = RevisionEvent::new(
            self.id.clone(),
            revision.clone(),
            parents,
            operation,
            transition,
        )?;
        let event_bytes = size(&event)?;
        let existing = self.nodes.iter().find(|n| n.event.id == event.id).cloned();
        let node = if let Some(existing) = existing {
            existing
        } else {
            if self.nodes.len() >= MAX_HISTORY_EVENTS {
                return Err(Error::Unsupported("Document history capacity reached; export the workspace or open a revision as a new document before continuing".into()));
            }
            let shared = self.nodes.iter().find(|n| n.event.revision == revision);
            let additional_molecule_bytes = if shared.is_some() { 0 } else { size(&next)? };
            if self
                .retained_bytes()
                .saturating_add(additional_molecule_bytes)
                .saturating_add(event_bytes)
                > MAX_HISTORY_BYTES
            {
                return Err(Error::Unsupported("Retained document history and typed actions exceed64MiB; export the workspace or fork a new document before continuing".into()));
            }
            let (molecule, bytes) = if let Some(shared) = shared {
                (shared.molecule.clone(), shared.bytes)
            } else {
                (Arc::new(next.clone()), additional_molecule_bytes)
            };
            Node {
                event: Arc::new(event),
                molecule,
                bytes,
                event_bytes,
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
    /// Every navigation position must name an existing exact event. Restoration
    /// never synthesizes events or infers ancestry from molecular content.
    pub fn restore_all(workspace: &Workspace) -> Result<Vec<Self>, Error> {
        workspace.validate()?;
        let shared = workspace
            .revisions
            .iter()
            .map(|(id, m)| Ok((id.clone(), (Arc::new(m.clone()), size(m)?))))
            .collect::<Result<BTreeMap<_, _>, Error>>()?;
        let mut output = Vec::new();
        for doc in &workspace.documents {
            let nodes = workspace
                .history
                .iter()
                .filter(|e| e.document_id == doc.id)
                .map(|event| {
                    let (molecule, bytes) = &shared[&event.revision];
                    Ok(Node {
                        event: Arc::new(event.clone()),
                        molecule: molecule.clone(),
                        bytes: *bytes,
                        event_bytes: size(event)?,
                    })
                })
                .collect::<Result<Vec<_>, Error>>()?;
            let head = doc
                .head_event
                .as_ref()
                .ok_or_else(|| Error::Format("document history head event is missing".into()))?
                .clone();
            let undo = doc.undo_events.clone();
            let redo = doc.redo_events.clone();
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
