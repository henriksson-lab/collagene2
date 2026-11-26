//! Nonmutating selected-document publication, never a whole-workspace fallback.
//! Native immutable events, run declarations and origin capsules remain exact.
//! This module computes documentary dependency scope, not scientific authority.
pub mod import;
mod sbol_source;
pub use sbol_source::KIND as SBOL_SOURCE_ATTACHMENT_KIND;

use collagene_core::{
    run_lineage::{self, EventBinding, ProductBinding},
    verification::Operation,
    workspace::{self, Workspace, WorkspaceDocument},
};
use serde::{Deserialize, Serialize};
use std::collections::{BTreeMap, BTreeSet, VecDeque};
pub const ENGINE_VERSION: &str = "library-document-closure-1";
const MAX_REFERENCE_EDGES: usize = 1_000_000;

#[derive(Debug, thiserror::Error)]
pub enum Error {
    #[error("native workspace validation: {0}")]
    Native(#[from] collagene_core::Error),
    #[error("selected native document/head/event does not match")]
    Selection,
    #[error("unsupported attachment reference contract: {0}")]
    UnsupportedAttachment(String),
    #[error("selected SBOL source cannot be safely published: {0}")]
    SbolSource(&'static str),
    #[error("unresolved native publication dependency: {0}")]
    Unresolved(&'static str),
    #[error("publication dependency limits exceeded")]
    TooLarge,
}
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Selection {
    pub document_id: String,
    pub revision: String,
    pub head_event: Option<String>,
}
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct HeadPin {
    pub revision: String,
    pub head_event: Option<String>,
}
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Omitted {
    pub documents: usize,
    pub revisions: usize,
    pub events: usize,
    pub candidates: usize,
    pub attachments: usize,
}
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Inventory {
    pub engine_version: String,
    pub selection: Selection,
    pub documents: BTreeSet<String>,
    pub revisions: BTreeSet<String>,
    pub events: BTreeSet<String>,
    pub candidates: BTreeSet<String>,
    pub attachments: BTreeSet<String>,
    pub document_heads: BTreeMap<String, HeadPin>,
    pub omitted: Omitted,
}
pub struct Export {
    pub workspace: Workspace,
    pub archive: Vec<u8>,
    pub inventory: Inventory,
}
/// Encoding can be offloaded to an existing browser/native workspace worker.
pub struct Projection {
    pub workspace: Workspace,
    pub inventory: Inventory,
}

/// The supported run schema has self-contained plans. Standalone candidate
/// proposals and workspace UI state are not required foreign keys of a plasmid.
fn presentation(kind: &str) -> bool {
    matches!(
        kind,
        "cloning-output-visibility-v1"
            | "oligo-table-v1"
            | "restriction-enzyme-sets-v1"
            | "named-digest-setup-v1"
            | "assembly-library-history-v1"
            // Exact versions read by App's draft restoration/autosave. Do not
            // prefix-match future schemas that might acquire different roles.
            | "assembly-draft-v1"
            | "assembly-draft-v2"
            | "assembly-draft-v3"
            | "assembly-draft-v4"
            | "assembly-draft-v5"
            | "assembly-draft-v6"
            | "assembly-draft-v7"
            | "assembly-draft-v8"
            | "assembly-draft-v9"
            | "assembly-draft-v10"
            // Named working designs, unfinished primer text and visibility
            // preferences produce no immutable plasmid/run event themselves.
            | "assembly-design-v1"
            | "primer-pair-draft-v1"
            | "primer-pair-draft-v2"
            | "annotation-display-v1"
            | "lab-library-bindings-v1"
            | "lab-library-branches-v1"
    )
}
fn run_evidence(kind: &str) -> bool {
    matches!(
        kind,
        collagene_core::reference_run_evidence::KIND
            | collagene_core::physical_run_evidence::KIND
            | collagene_core::processed_run_evidence::KIND
    )
}
// Explicit known-schema index projections. The original Workspace::validate
// checks complete run/evidence schemas before these borrowed-size indexes exist.
#[derive(Deserialize)]
struct RunIndex {
    run_id: String,
    products: Vec<ProductBinding>,
}
#[derive(Deserialize)]
struct EvidenceOwner {
    run_id: String,
}

struct Graph<'a> {
    source: &'a Workspace,
    events: BTreeMap<&'a str, &'a workspace::RevisionEvent>,
    run_headers: BTreeSet<String>,
    run_assets: BTreeMap<String, String>,
    owners: BTreeMap<String, String>,
    companions: BTreeMap<String, Vec<String>>,
    kept_events: BTreeSet<String>,
    kept_revisions: BTreeSet<String>,
    kept_runs: BTreeSet<String>,
    kept_assets: BTreeSet<String>,
    pending_events: VecDeque<String>,
    pending_runs: VecDeque<String>,
    edges: usize,
}
impl<'a> Graph<'a> {
    fn new(source: &'a Workspace) -> Result<Self, Error> {
        source.validate()?;
        let mut g = Self {
            source,
            events: source.history.iter().map(|e| (e.id.as_str(), e)).collect(),
            run_headers: BTreeSet::new(),
            run_assets: BTreeMap::new(),
            owners: BTreeMap::new(),
            companions: BTreeMap::new(),
            kept_events: BTreeSet::new(),
            kept_revisions: BTreeSet::new(),
            kept_runs: BTreeSet::new(),
            kept_assets: BTreeSet::new(),
            pending_events: VecDeque::new(),
            pending_runs: VecDeque::new(),
            edges: 0,
        };
        for (hash, a) in &source.attachments {
            if a.kind == run_lineage::WORKSPACE_ATTACHMENT_KIND {
                let header: RunIndex =
                    serde_json::from_slice(&a.bytes).map_err(|_| Error::Unresolved("run index"))?;
                if g.run_headers.contains(&header.run_id) {
                    return Err(Error::Unresolved("duplicate run"));
                }
                for p in &header.products {
                    g.edge()?;
                    if g.owners
                        .insert(p.event_id.clone(), header.run_id.clone())
                        .is_some()
                    {
                        return Err(Error::Unresolved("multiple run owners for a product event"));
                    }
                }
                g.run_assets.insert(header.run_id.clone(), hash.clone());
                g.run_headers.insert(header.run_id);
            } else if run_evidence(&a.kind) {
                let owner: EvidenceOwner = serde_json::from_slice(&a.bytes)
                    .map_err(|_| Error::Unresolved("run evidence owner"))?;
                g.companions
                    .entry(owner.run_id)
                    .or_default()
                    .push(hash.clone());
            } else if a.kind == sbol_source::KIND {
                // Ownership is derived later from retained native revisions.
                // An unrelated documentary original must not block publication.
            } else if !presentation(&a.kind) {
                // The attachment extension point has no generic foreign-key
                // contract. Do not guess ownership by matching arbitrary text,
                // silently discard unknown ancestry, or publish unrelated data.
                return Err(Error::UnsupportedAttachment(a.kind.clone()));
            }
        }
        Ok(g)
    }
    fn edge(&mut self) -> Result<(), Error> {
        self.edges = self
            .edges
            .checked_add(1)
            .filter(|n| *n <= MAX_REFERENCE_EDGES)
            .ok_or(Error::TooLarge)?;
        Ok(())
    }
    fn event(&mut self, id: &str) -> Result<(), Error> {
        self.edge()?;
        if !self.events.contains_key(id) {
            return Err(Error::Unresolved("event"));
        }
        if self.kept_events.insert(id.into()) {
            self.pending_events.push_back(id.into());
        }
        Ok(())
    }
    fn binding(&mut self, b: &EventBinding) -> Result<(), Error> {
        let event = self
            .events
            .get(b.event_id.as_str())
            .ok_or(Error::Unresolved("bound event"))?;
        if event.document_id != b.document_id || event.revision != b.revision {
            return Err(Error::Unresolved("exact event pin"));
        }
        self.event(&b.event_id)
    }
    fn run(&mut self, id: &str) -> Result<(), Error> {
        self.edge()?;
        if !self.run_headers.contains(id) {
            return Err(Error::Unresolved("owning origin run"));
        }
        if self.kept_runs.insert(id.into()) {
            self.pending_runs.push_back(id.into());
        }
        Ok(())
    }
    fn resolve(&mut self) -> Result<(), Error> {
        while !self.pending_events.is_empty() || !self.pending_runs.is_empty() {
            while let Some(id) = self.pending_events.pop_front() {
                let event = *self
                    .events
                    .get(id.as_str())
                    .ok_or(Error::Unresolved("event"))?;
                if !self.source.revisions.contains_key(&event.revision) {
                    return Err(Error::Unresolved("event revision"));
                }
                self.kept_revisions.insert(event.revision.clone());
                for parent in &event.parents {
                    self.event(parent)?;
                }
                if let Some(owner) = self.owners.get(&id).cloned() {
                    self.run(&owner)?;
                }
                // This is a reserved native producer label, not a free-form
                // user description contract. An orphan cannot be published as
                // if its original reconstruction record were present.
                else if event.operation.starts_with("Cloning run ") {
                    return Err(Error::Unresolved("orphan cloning-run product event"));
                }
            }
            while let Some(id) = self.pending_runs.pop_front() {
                let hash = self
                    .run_assets
                    .get(&id)
                    .cloned()
                    .ok_or(Error::Unresolved("run attachment"))?;
                self.kept_assets.insert(hash.clone());
                if let Some(assets) = self.companions.get(&id) {
                    self.kept_assets.extend(assets.iter().cloned());
                }
                let asset = self
                    .source
                    .attachments
                    .get(&hash)
                    .ok_or(Error::Unresolved("run bytes"))?;
                let run = run_lineage::load_run_record(&asset.bytes)?;
                self.binding(&run.target)?;
                for b in &run.inputs {
                    self.binding(&b.event)?;
                }
                for p in &run.products {
                    self.binding(&EventBinding {
                        document_id: p.document_id.clone(),
                        event_id: p.event_id.clone(),
                        revision: p.revision.clone(),
                    })?;
                }
                // Native origin validation already proves each capsule equals
                // its retained run declaration, including nested candidate bytes.
                // Following owner plans avoids independently interpreting opaque
                // serialized candidate packages or inventing material state.
                for op in &run.plan.steps {
                    if let Operation::ReconstructRunProduct { origin, .. } = op {
                        self.run(&origin.selection.run_id)?;
                        self.binding(&origin.selection.event)?;
                        self.binding(&origin.target_event)?;
                        for b in &origin.inputs {
                            self.binding(&b.event)?;
                        }
                        for p in &origin.products {
                            self.binding(&EventBinding {
                                document_id: p.document_id.clone(),
                                event_id: p.event_id.clone(),
                                revision: p.revision.clone(),
                            })?;
                        }
                    }
                }
            }
        }
        Ok(())
    }
    fn inventory(mut self, selection: &Selection) -> Result<Inventory, Error> {
        let selected = self
            .source
            .documents
            .iter()
            .find(|d| d.id == selection.document_id)
            .ok_or(Error::Selection)?;
        if selected.head != selection.revision || selected.head_event != selection.head_event {
            return Err(Error::Selection);
        }
        self.kept_revisions.insert(selection.revision.clone());
        if let Some(event) = &selection.head_event {
            self.event(event)?;
        }
        self.resolve()?;
        let (sbol_assets, used) = sbol_source::required_assets(
            self.source,
            &self.kept_revisions,
            MAX_REFERENCE_EDGES.saturating_sub(self.edges),
        )?;
        self.edges = self
            .edges
            .checked_add(used)
            .filter(|n| *n <= MAX_REFERENCE_EDGES)
            .ok_or(Error::TooLarge)?;
        self.kept_assets.extend(sbol_assets);
        let mut document_heads: BTreeMap<String, HeadPin> = BTreeMap::new();
        // Choose a stable retained head, rather than a dependency document's
        // possibly unrelated current/private branch. Keep all required events.
        // Lexical event-ID order is deterministic and grants no new ancestry.
        for id in &self.kept_events {
            let event = self
                .events
                .get(id.as_str())
                .ok_or(Error::Unresolved("retained event"))?;
            document_heads
                .entry(event.document_id.clone())
                .or_insert_with(|| HeadPin {
                    revision: event.revision.clone(),
                    head_event: Some(event.id.clone()),
                });
        }
        document_heads.insert(
            selection.document_id.clone(),
            HeadPin {
                revision: selection.revision.clone(),
                head_event: selection.head_event.clone(),
            },
        );
        for id in document_heads.keys() {
            if !self.source.documents.iter().any(|d| &d.id == id) {
                return Err(Error::Unresolved("source document identity"));
            }
        }
        let documents = document_heads.keys().cloned().collect::<BTreeSet<_>>();
        let omitted = Omitted {
            documents: self.source.documents.len() - documents.len(),
            revisions: self.source.revisions.len() - self.kept_revisions.len(),
            events: self.source.history.len() - self.kept_events.len(),
            candidates: self.source.candidates.len(),
            attachments: self.source.attachments.len() - self.kept_assets.len(),
        };
        Ok(Inventory {
            engine_version: ENGINE_VERSION.into(),
            selection: selection.clone(),
            documents,
            revisions: self.kept_revisions,
            events: self.kept_events,
            candidates: BTreeSet::new(),
            attachments: self.kept_assets,
            document_heads,
            omitted,
        })
    }
}
/// Derive permitted native scope independently from the source graph. No clone
/// of unneeded molecules/attachments, scientific replay or filesystem I/O occurs.
pub fn inspect_document_closure(
    source: &Workspace,
    selection: &Selection,
) -> Result<Inventory, Error> {
    Graph::new(source)?.inventory(selection)
}

/// Export a fresh view; the caller's local library, history and unpublished
/// branches remain untouched. Whole-workspace sharing is a separate host action.
pub fn project_document_closure(
    source: &Workspace,
    selection: &Selection,
) -> Result<Projection, Error> {
    let inventory = inspect_document_closure(source, selection)?;
    let documents = inventory
        .document_heads
        .iter()
        .map(|(id, pin)| WorkspaceDocument {
            id: id.clone(),
            head: pin.revision.clone(),
            head_event: pin.head_event.clone(),
            undo: vec![],
            redo: vec![],
            undo_events: vec![],
            redo_events: vec![],
        })
        .collect();
    let frame = Workspace {
        name: "Plasmid publication".into(),
        active_document: selection.document_id.clone(),
        documents,
        revisions: source
            .revisions
            .iter()
            .filter(|(id, _)| inventory.revisions.contains(*id))
            .map(|(id, m)| (id.clone(), m.clone()))
            .collect(),
        history: source
            .history
            .iter()
            .filter(|e| inventory.events.contains(&e.id))
            .cloned()
            .collect(),
        candidates: vec![],
        active_candidate: None,
        attachments: source
            .attachments
            .iter()
            .filter(|(id, _)| inventory.attachments.contains(*id))
            .map(|(id, a)| (id.clone(), a.clone()))
            .collect(),
    };
    frame.validate()?;
    Ok(Projection {
        workspace: frame,
        inventory,
    })
}
/// Synchronous convenience for native hosts; browser hosts can encode Projection
/// through their existing bounded workspace worker.
pub fn document_closure(source: &Workspace, selection: &Selection) -> Result<Export, Error> {
    let projected = project_document_closure(source, selection)?;
    let archive = workspace::save_workspace(&projected.workspace)?;
    Ok(Export {
        workspace: projected.workspace,
        archive,
        inventory: projected.inventory,
    })
}
