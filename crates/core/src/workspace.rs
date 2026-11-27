//! Durable revision/workspace archive. ZIP is an in-memory container, never a filesystem action.
//! Content hashes bind evidence, not permission, trusted chemistry or scientific verification.
use crate::{
    Alphabet, Error, Feature, Molecule, Sequence, Topology,
    transfer::Candidate,
    verification::{DesignPlan, Operation, PlanInput},
};
use base64::{Engine, engine::general_purpose::STANDARD};
use serde::{
    Deserialize, Serialize,
    de::{DeserializeOwned, MapAccess, SeqAccess, Visitor},
};
use std::{
    collections::{BTreeMap, BTreeSet, VecDeque},
    fmt,
    io::{self, Cursor, Read, Write},
};
use zip::{CompressionMethod, ZipArchive, ZipWriter, write::SimpleFileOptions};

pub const MAX_ARCHIVE_BYTES: usize = 64 * 1024 * 1024;
pub const MAX_EXPANDED_BYTES: usize = 128 * 1024 * 1024;
pub const MAX_ENTRY_BYTES: usize = 32 * 1024 * 1024;
pub const MAX_MANIFEST_BYTES: usize = 16 * 1024 * 1024;
pub const MAX_DOCUMENTS: usize = 1000;
pub const MAX_REVISIONS: usize = 10_000;
pub const MAX_HISTORY_EVENTS: usize = 20_000;
pub const MAX_CANDIDATES: usize = 128;
pub const MAX_ATTACHMENTS: usize = 1000;
const MAX_ENTRIES: usize = 50_000;
const MAX_STACK: usize = 1000;
const FORMAT: &str = "collagene-workspace";
const VERSION: u32 = 1;
pub const PLAN_SCHEMA_VERSION: u32 = 23;

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct WorkspaceDocument {
    pub id: String,
    pub head: String,
    pub undo: Vec<String>,
    pub redo: Vec<String>,
    pub head_event: Option<String>,
    pub undo_events: Vec<String>,
    pub redo_events: Vec<String>,
}
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct RevisionEvent {
    pub id: String,
    pub document_id: String,
    pub revision: String,
    /// Ordered parent event IDs. Cross-document parents can describe a derived design.
    pub parents: Vec<String>,
    /// A descriptive operation label, not executable replay code or a grant of authority.
    pub operation: String,
}
impl RevisionEvent {
    pub fn new(
        document_id: impl Into<String>,
        revision: impl Into<String>,
        parents: Vec<String>,
        operation: impl Into<String>,
    ) -> Result<Self, Error> {
        let mut event = Self {
            id: String::new(),
            document_id: document_id.into(),
            revision: revision.into(),
            parents,
            operation: operation.into(),
        };
        event.check_fields()?;
        event.id = event.expected_id()?;
        Ok(event)
    }
    fn expected_id(&self) -> Result<String, Error> {
        let mut hash = blake3::Hasher::new();
        hash.update(b"collagene-history-event-v1\0");
        hash.update(
            &serde_json::to_vec(&(
                &self.document_id,
                &self.revision,
                &self.parents,
                &self.operation,
            ))
            .map_err(format_error)?,
        );
        Ok(format!("event-{}", hash.finalize().to_hex()))
    }
    fn check_fields(&self) -> Result<(), Error> {
        if !identifier(&self.document_id)
            || !digest(&self.revision)
            || self.parents.len() > 128
            || self.operation.len() > 4096
            || self
                .operation
                .chars()
                .any(|c| c.is_control() && !matches!(c, '\n' | '\r' | '\t'))
            || self.parents.iter().any(|id| !event_id(id))
            || self.parents.iter().collect::<BTreeSet<_>>().len() != self.parents.len()
        {
            return Err(bad("invalid history event fields"));
        }
        Ok(())
    }
}
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct WorkspaceAttachment {
    pub kind: String,
    pub name: String,
    pub media_type: String,
    pub bytes: Vec<u8>,
}
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Workspace {
    pub name: String,
    pub active_document: String,
    pub documents: Vec<WorkspaceDocument>,
    pub revisions: BTreeMap<String, Molecule>,
    pub history: Vec<RevisionEvent>,
    pub candidates: Vec<Candidate>,
    pub active_candidate: Option<String>,
    pub attachments: BTreeMap<String, WorkspaceAttachment>,
}
impl Workspace {
    pub fn from_documents(
        name: impl Into<String>,
        documents: Vec<Molecule>,
        active: usize,
    ) -> Result<Self, Error> {
        if documents.is_empty() || documents.len() > MAX_DOCUMENTS || active >= documents.len() {
            return Err(bad("workspace needs an available active document"));
        }
        let mut workspace = Self {
            name: name.into(),
            active_document: format!("document-{}", active + 1),
            documents: vec![],
            revisions: BTreeMap::new(),
            history: vec![],
            candidates: vec![],
            active_candidate: None,
            attachments: BTreeMap::new(),
        };
        for (i, molecule) in documents.into_iter().enumerate() {
            molecule_size(&molecule)?;
            molecule.validate()?;
            let revision = molecule.revision();
            let id = format!("document-{}", i + 1);
            let event = RevisionEvent::new(&id, &revision, vec![], "Imported document")?;
            workspace.documents.push(WorkspaceDocument {
                id,
                head: revision.clone(),
                undo: vec![],
                redo: vec![],
                head_event: Some(event.id.clone()),
                undo_events: vec![],
                redo_events: vec![],
            });
            workspace.history.push(event);
            workspace.revisions.insert(revision, molecule);
        }
        workspace.validate()?;
        Ok(workspace)
    }
    pub fn head(&self, document_id: &str) -> Result<&Molecule, Error> {
        self.documents
            .iter()
            .find(|d| d.id == document_id)
            .and_then(|d| self.revisions.get(&d.head))
            .ok_or_else(|| bad("workspace document/head is unavailable"))
    }
    pub fn unresolved_candidate_targets(&self) -> Vec<String> {
        self.candidates
            .iter()
            .map(|c| &c.plan.target_revision)
            .filter(|r| !self.revisions.contains_key(*r))
            .cloned()
            .collect::<BTreeSet<_>>()
            .into_iter()
            .collect()
    }
    pub fn validate(&self) -> Result<(), Error> {
        if self.name.trim().is_empty()
            || self.name.len() > 1024
            || self.name.chars().any(char::is_control)
            || self.documents.is_empty()
            || self.documents.len() > MAX_DOCUMENTS
            || self.revisions.is_empty()
            || self.revisions.len() > MAX_REVISIONS
            || self.history.len() > MAX_HISTORY_EVENTS
            || self.candidates.len() > MAX_CANDIDATES
            || self.attachments.len() > MAX_ATTACHMENTS
        {
            return Err(bad("invalid workspace dimensions or name"));
        }
        json_size(
            &(
                &self.name,
                &self.active_document,
                &self.documents,
                &self.history,
                &self.active_candidate,
            ),
            MAX_MANIFEST_BYTES,
        )?;
        let mut expanded = 0usize;
        let mut metadata_fields = 0usize;
        for (revision, molecule) in &self.revisions {
            let size = molecule_size(molecule)?;
            charge(&mut expanded, size, MAX_EXPANDED_BYTES)?;
            metadata_fields = metadata_fields
                .checked_add(molecule.metadata.len())
                .ok_or(Error::TooLarge)?;
            if metadata_fields > 200_000 || molecule.metadata.len() > 10_000 {
                return Err(Error::TooLarge);
            }
            molecule.validate()?;
            if !digest(revision) || molecule.revision() != *revision {
                return Err(bad("workspace revision hash does not match its molecule"));
            }
        }
        let mut documents = BTreeSet::new();
        for document in &self.documents {
            if !identifier(&document.id)
                || !documents.insert(document.id.as_str())
                || document.undo.len() > MAX_STACK
                || document.redo.len() > MAX_STACK
            {
                return Err(bad(
                    "invalid or duplicate workspace document identity/history stack",
                ));
            }
            for revision in std::iter::once(&document.head)
                .chain(&document.undo)
                .chain(&document.redo)
            {
                if !self.revisions.contains_key(revision) {
                    return Err(bad(
                        "document history references a missing molecule revision",
                    ));
                }
            }
        }
        if !documents.contains(self.active_document.as_str()) {
            return Err(bad("active workspace document is unavailable"));
        }
        let mut events = BTreeMap::new();
        for event in &self.history {
            event.check_fields()?;
            if event.id != event.expected_id()?
                || events.insert(event.id.as_str(), event).is_some()
                || !documents.contains(event.document_id.as_str())
                || !self.revisions.contains_key(&event.revision)
            {
                return Err(bad("invalid, duplicate or dangling history event"));
            }
        }
        let mut indegrees = BTreeMap::new();
        let mut children: BTreeMap<&str, Vec<&str>> = BTreeMap::new();
        for event in &self.history {
            for parent in &event.parents {
                if !events.contains_key(parent.as_str()) {
                    return Err(bad("history event references a missing parent"));
                }
                children.entry(parent).or_default().push(&event.id);
            }
            indegrees.insert(event.id.as_str(), event.parents.len());
        }
        let mut queue = indegrees
            .iter()
            .filter(|(_, n)| **n == 0)
            .map(|(id, _)| *id)
            .collect::<VecDeque<_>>();
        let mut visited = 0;
        while let Some(id) = queue.pop_front() {
            visited += 1;
            if let Some(nodes) = children.get(id) {
                for child in nodes {
                    let n = indegrees
                        .get_mut(child)
                        .ok_or_else(|| bad("invalid history topology"))?;
                    *n -= 1;
                    if *n == 0 {
                        queue.push_back(child);
                    }
                }
            }
        }
        if visited != self.history.len() {
            return Err(bad("workspace history contains a cycle"));
        }
        for document in &self.documents {
            let check = |id: &str, revision: &str| -> Result<(), Error> {
                let event = events
                    .get(id)
                    .ok_or_else(|| bad("document navigation references a missing history event"))?;
                if event.document_id != document.id || event.revision != revision {
                    return Err(bad(
                        "document navigation event does not match its revision/identity",
                    ));
                }
                Ok(())
            };
            if let Some(id) = &document.head_event {
                check(id, &document.head)?;
            }
            for (stack, revisions) in [
                (&document.undo_events, &document.undo),
                (&document.redo_events, &document.redo),
            ] {
                if !stack.is_empty() {
                    if stack.len() != revisions.len() {
                        return Err(bad(
                            "document navigation event/revision stack lengths disagree",
                        ));
                    }
                    for (id, revision) in stack.iter().zip(revisions) {
                        check(id, revision)?;
                    }
                }
            }
        }
        let mut plans = BTreeSet::new();
        for candidate in &self.candidates {
            crate::run_material::plan_preflight(&candidate.plan)?;
            charge(
                &mut expanded,
                json_size(&candidate.plan, MAX_ENTRY_BYTES)?,
                MAX_EXPANDED_BYTES,
            )?;
            candidate.validate()?;
            crate::run_material::validate_plan_origins(self, &candidate.plan)?;
            if !plans.insert(candidate.plan_revision.as_str()) {
                return Err(bad("duplicate workspace candidate revision"));
            }
            for input in &candidate.plan.inputs {
                if self.revisions.get(&input.revision) != Some(&input.molecule) {
                    return Err(bad(
                        "candidate input is missing or differs from its workspace revision",
                    ));
                }
            }
        }
        if self
            .active_candidate
            .as_ref()
            .is_some_and(|id| !plans.contains(id.as_str()))
        {
            return Err(bad("active candidate is unavailable"));
        }
        for (hash, attachment) in &self.attachments {
            if !digest(hash)
                || blake3::hash(&attachment.bytes).to_hex().as_str() != hash
                || attachment.bytes.len() > MAX_ENTRY_BYTES
                || attachment.kind.is_empty()
                || attachment.kind.len() > 64
                || !attachment
                    .kind
                    .bytes()
                    .all(|b| b.is_ascii_alphanumeric() || b == b'_' || b == b'-')
                || attachment.name.len() > 1024
                || attachment.name.chars().any(char::is_control)
                || attachment.media_type.is_empty()
                || attachment.media_type.len() > 128
                || !attachment.media_type.bytes().all(|b| b.is_ascii_graphic())
            {
                return Err(bad("invalid workspace attachment metadata/hash"));
            }
            if attachment.kind == crate::alignment_formats::WORKSPACE_ATTACHMENT_KIND {
                if attachment.media_type != "application/json" {
                    return Err(bad("alignment attachment must use application/json"));
                }
                crate::alignment_formats::load_alignment(&attachment.bytes)?;
            }
            if attachment.kind == crate::genome::WORKSPACE_ATTACHMENT_KIND {
                if attachment.media_type != "application/json" {
                    return Err(bad("genome reference link must use application/json"));
                }
                crate::genome::load_link(&attachment.bytes)?;
            }
            if matches!(
                attachment.kind.as_str(),
                crate::read_alignment_archive::SAM_SOURCE_KIND
                    | crate::read_alignment_archive::BAM_SOURCE_KIND
            ) {
                // Source metadata is validated here; parsed collections and checked
                // mappings must be recreated from the exact bytes on reopening.
                crate::read_alignment_archive::inspect(attachment)?;
            }
            if matches!(
                attachment.kind.as_str(),
                crate::variants_archive::VCF_SOURCE_KIND | crate::variants_archive::BCF_SOURCE_KIND
            ) {
                // Archive membership restores exact source bytes, never parsed calls.
                crate::variants_archive::inspect(attachment)?;
            }
            charge(&mut expanded, attachment.bytes.len(), MAX_EXPANDED_BYTES)?;
        }
        crate::run_lineage::validate_workspace_runs(self)?;
        crate::material_derivation::validate_workspace_records(self)?;
        crate::reference_run_evidence::validate_workspace(self)?;
        crate::contig_history::validate_workspace_assemblies(self)?;
        Ok(())
    }
}
fn identifier(value: &str) -> bool {
    !value.is_empty()
        && value.len() <= 128
        && value
            .bytes()
            .all(|b| b.is_ascii_alphanumeric() || matches!(b, b'_' | b'-' | b'.' | b':'))
}
fn digest(value: &str) -> bool {
    value.len() == 64
        && value
            .bytes()
            .all(|b| b.is_ascii_digit() || (b'a'..=b'f').contains(&b))
}
fn event_id(value: &str) -> bool {
    value.strip_prefix("event-").is_some_and(digest)
}
fn bad(message: &str) -> Error {
    Error::Format(message.into())
}
fn format_error(_: impl fmt::Display) -> Error {
    bad("workspace archive data is invalid")
}
fn charge(total: &mut usize, amount: usize, maximum: usize) -> Result<(), Error> {
    *total = total
        .checked_add(amount)
        .filter(|n| *n <= maximum)
        .ok_or(Error::TooLarge)?;
    Ok(())
}
struct Counter {
    bytes: usize,
    maximum: usize,
}
impl Write for Counter {
    fn write(&mut self, bytes: &[u8]) -> io::Result<usize> {
        self.bytes = self
            .bytes
            .checked_add(bytes.len())
            .filter(|n| *n <= self.maximum)
            .ok_or_else(|| io::Error::other("workspace size limit"))?;
        Ok(bytes.len())
    }
    fn flush(&mut self) -> io::Result<()> {
        Ok(())
    }
}
fn json_size(value: &impl Serialize, maximum: usize) -> Result<usize, Error> {
    let mut out = Counter { bytes: 0, maximum };
    serde_json::to_writer(&mut out, value).map_err(|_| Error::TooLarge)?;
    Ok(out.bytes)
}
fn molecule_size(value: &Molecule) -> Result<usize, Error> {
    json_size(value, MAX_ENTRY_BYTES)
}

#[derive(Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
struct BlobInfo {
    length: usize,
}
#[derive(Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
struct Manifest {
    format: String,
    version: u32,
    molecule_schema: u32,
    plan_schema: u32,
    name: String,
    active_document: String,
    documents: Vec<WorkspaceDocument>,
    history: Vec<RevisionEvent>,
    revisions: BTreeMap<String, String>,
    candidates: Vec<String>,
    active_candidate: Option<String>,
    unresolved_candidate_targets: Vec<String>,
    attachments: BTreeMap<String, AttachmentRecord>,
    blobs: BTreeMap<String, BlobInfo>,
}
#[derive(Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
struct MoleculeRecord {
    name: String,
    sequence: String,
    alphabet: Alphabet,
    topology: Topology,
    features: String,
    metadata: BTreeMap<String, ValueRecord>,
    original_genbank: Option<String>,
}
#[derive(Debug, Serialize, Deserialize)]
#[serde(
    tag = "encoding",
    content = "blob",
    rename_all = "snake_case",
    deny_unknown_fields
)]
enum ValueRecord {
    Utf8(String),
    Base64(String),
}
#[derive(Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
struct PlanRecord {
    plan_revision: String,
    name: String,
    target_revision: String,
    inputs: Vec<InputRecord>,
    steps: String,
}
#[derive(Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
struct InputRecord {
    id: String,
    revision: String,
}
#[derive(Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
struct AttachmentRecord {
    kind: String,
    name: String,
    media_type: String,
    blob: String,
}
#[derive(Default)]
struct Blobs {
    values: BTreeMap<String, Vec<u8>>,
    total: usize,
}
impl Blobs {
    fn put(&mut self, bytes: &[u8]) -> Result<String, Error> {
        if bytes.len() > MAX_ENTRY_BYTES {
            return Err(Error::TooLarge);
        }
        let hash = blake3::hash(bytes).to_hex().to_string();
        if !self.values.contains_key(&hash) {
            if self.values.len() + 1 >= MAX_ENTRIES {
                return Err(Error::TooLarge);
            }
            charge(&mut self.total, bytes.len(), MAX_EXPANDED_BYTES)?;
            self.values.insert(hash.clone(), bytes.to_vec());
        }
        Ok(hash)
    }
    fn json(&mut self, value: &impl Serialize) -> Result<String, Error> {
        self.put(&serde_json::to_vec(value).map_err(format_error)?)
    }
}
/// Deterministic, uncompressed ZIP writing: deduplication avoids repeated source evidence,
/// while Stored entries keep save cost independent of compression heuristics/CPU backends.
pub fn save_workspace(workspace: &Workspace) -> Result<Vec<u8>, Error> {
    workspace.validate()?;
    let mut blobs = Blobs::default();
    let mut revisions = BTreeMap::new();
    for (revision, molecule) in &workspace.revisions {
        let mut metadata = BTreeMap::new();
        for (key, value) in &molecule.metadata {
            let stored = if key.ends_with("base64") {
                STANDARD
                    .decode(value)
                    .ok()
                    .filter(|bytes| STANDARD.encode(bytes) == *value)
                    .map(|bytes| blobs.put(&bytes).map(ValueRecord::Base64))
                    .transpose()?
            } else {
                None
            };
            metadata.insert(
                key.clone(),
                match stored {
                    Some(value) => value,
                    None => ValueRecord::Utf8(blobs.put(value.as_bytes())?),
                },
            );
        }
        let record = MoleculeRecord {
            name: blobs.put(molecule.name.as_bytes())?,
            sequence: blobs.put(molecule.sequence.as_str().as_bytes())?,
            alphabet: molecule.sequence.alphabet(),
            topology: molecule.topology,
            features: blobs.json(&molecule.features)?,
            metadata,
            original_genbank: molecule
                .original_genbank
                .as_ref()
                .map(|text| blobs.put(text.as_bytes()))
                .transpose()?,
        };
        revisions.insert(revision.clone(), blobs.json(&record)?);
    }
    let mut candidates = Vec::new();
    for candidate in &workspace.candidates {
        let plan = &candidate.plan;
        let record = PlanRecord {
            plan_revision: candidate.plan_revision.clone(),
            name: plan.name.clone(),
            target_revision: plan.target_revision.clone(),
            inputs: plan
                .inputs
                .iter()
                .map(|i| InputRecord {
                    id: i.id.clone(),
                    revision: i.revision.clone(),
                })
                .collect(),
            steps: blobs.json(&plan.steps)?,
        };
        candidates.push(blobs.json(&record)?);
    }
    let mut attachments = BTreeMap::new();
    for (hash, attachment) in &workspace.attachments {
        attachments.insert(
            hash.clone(),
            AttachmentRecord {
                kind: attachment.kind.clone(),
                name: attachment.name.clone(),
                media_type: attachment.media_type.clone(),
                blob: blobs.put(&attachment.bytes)?,
            },
        );
    }
    let manifest = Manifest {
        format: FORMAT.into(),
        version: VERSION,
        molecule_schema: 1,
        plan_schema: PLAN_SCHEMA_VERSION,
        name: workspace.name.clone(),
        active_document: workspace.active_document.clone(),
        documents: workspace.documents.clone(),
        history: workspace.history.clone(),
        revisions,
        candidates,
        active_candidate: workspace.active_candidate.clone(),
        unresolved_candidate_targets: workspace.unresolved_candidate_targets(),
        attachments,
        blobs: blobs
            .values
            .iter()
            .map(|(hash, value)| {
                (
                    hash.clone(),
                    BlobInfo {
                        length: value.len(),
                    },
                )
            })
            .collect(),
    };
    json_size(&manifest, MAX_MANIFEST_BYTES)?;
    let manifest = serde_json::to_vec(&manifest).map_err(format_error)?;
    if manifest.len().saturating_add(blobs.total) > MAX_EXPANDED_BYTES {
        return Err(Error::TooLarge);
    }
    let mut predicted = manifest
        .len()
        .checked_add(blobs.total)
        .ok_or(Error::TooLarge)?;
    charge(
        &mut predicted,
        22 + 76 + 2 * "manifest.json".len(),
        MAX_ARCHIVE_BYTES,
    )?;
    for hash in blobs.values.keys() {
        charge(&mut predicted, 76 + 2 * (6 + hash.len()), MAX_ARCHIVE_BYTES)?;
    }
    let mut zip = ZipWriter::new(Cursor::new(Vec::new()));
    let options = SimpleFileOptions::default()
        .compression_method(CompressionMethod::Stored)
        .last_modified_time(
            zip::DateTime::from_date_and_time(1980, 1, 1, 0, 0, 0).map_err(format_error)?,
        )
        .unix_permissions(0o600);
    zip.start_file("manifest.json", options)
        .map_err(format_error)?;
    zip.write_all(&manifest).map_err(format_error)?;
    for (hash, value) in &blobs.values {
        zip.start_file(format!("blobs/{hash}"), options)
            .map_err(format_error)?;
        zip.write_all(value).map_err(format_error)?;
    }
    let bytes = zip.finish().map_err(format_error)?.into_inner();
    if bytes.len() > MAX_ARCHIVE_BYTES {
        return Err(Error::TooLarge);
    }
    Ok(bytes)
}

// serde_json normally accepts duplicate object keys. Archive maps cannot do so: which
// revision/manifest value wins must not depend on a parser's duplicate-key policy.
struct UniqueJson(serde_json::Value);
impl<'de> Deserialize<'de> for UniqueJson {
    fn deserialize<D: serde::Deserializer<'de>>(deserializer: D) -> Result<Self, D::Error> {
        struct V;
        impl<'de> Visitor<'de> for V {
            type Value = UniqueJson;
            fn expecting(&self, f: &mut fmt::Formatter) -> fmt::Result {
                f.write_str("JSON without duplicate object keys")
            }
            fn visit_bool<E: serde::de::Error>(self, v: bool) -> Result<Self::Value, E> {
                Ok(UniqueJson(v.into()))
            }
            fn visit_i64<E: serde::de::Error>(self, v: i64) -> Result<Self::Value, E> {
                Ok(UniqueJson(v.into()))
            }
            fn visit_u64<E: serde::de::Error>(self, v: u64) -> Result<Self::Value, E> {
                Ok(UniqueJson(v.into()))
            }
            fn visit_f64<E: serde::de::Error>(self, v: f64) -> Result<Self::Value, E> {
                serde_json::Number::from_f64(v)
                    .map(|n| UniqueJson(n.into()))
                    .ok_or_else(|| E::custom("non-finite JSON"))
            }
            fn visit_str<E: serde::de::Error>(self, v: &str) -> Result<Self::Value, E> {
                Ok(UniqueJson(v.into()))
            }
            fn visit_string<E: serde::de::Error>(self, v: String) -> Result<Self::Value, E> {
                Ok(UniqueJson(v.into()))
            }
            fn visit_unit<E: serde::de::Error>(self) -> Result<Self::Value, E> {
                Ok(UniqueJson(serde_json::Value::Null))
            }
            fn visit_none<E: serde::de::Error>(self) -> Result<Self::Value, E> {
                self.visit_unit()
            }
            fn visit_seq<A: SeqAccess<'de>>(self, mut seq: A) -> Result<Self::Value, A::Error> {
                let mut values = Vec::new();
                while let Some(value) = seq.next_element::<UniqueJson>()? {
                    values.push(value.0);
                }
                Ok(UniqueJson(values.into()))
            }
            fn visit_map<A: MapAccess<'de>>(self, mut map: A) -> Result<Self::Value, A::Error> {
                let mut values = serde_json::Map::new();
                while let Some((key, value)) = map.next_entry::<String, UniqueJson>()? {
                    if values.insert(key, value.0).is_some() {
                        return Err(serde::de::Error::custom("duplicate JSON object key"));
                    }
                }
                Ok(UniqueJson(values.into()))
            }
        }
        deserializer.deserialize_any(V)
    }
}
fn json<T: DeserializeOwned>(bytes: &[u8]) -> Result<T, Error> {
    let mut parser = serde_json::Deserializer::from_slice(bytes);
    let value = UniqueJson::deserialize(&mut parser).map_err(format_error)?;
    parser.end().map_err(format_error)?;
    serde_json::from_value(value.0).map_err(format_error)
}
fn json_exact<T: DeserializeOwned + Serialize>(bytes: &[u8]) -> Result<T, Error> {
    let value: serde_json::Value = json(bytes)?;
    let parsed: T = serde_json::from_value(value.clone()).map_err(format_error)?;
    if serde_json::to_value(&parsed).map_err(format_error)? != value {
        return Err(bad(
            "workspace typed data contains unsupported or noncanonical fields",
        ));
    }
    Ok(parsed)
}
fn text(bytes: &[u8]) -> Result<String, Error> {
    String::from_utf8(bytes.to_vec()).map_err(format_error)
}
fn get<'a>(blobs: &'a BTreeMap<String, Vec<u8>>, hash: &str) -> Result<&'a [u8], Error> {
    if !digest(hash) {
        return Err(bad("invalid workspace blob reference"));
    }
    blobs
        .get(hash)
        .map(Vec::as_slice)
        .ok_or_else(|| bad("missing workspace blob"))
}
fn u16_at(bytes: &[u8], offset: usize) -> Result<u16, Error> {
    Ok(u16::from_le_bytes(
        bytes
            .get(offset..offset + 2)
            .ok_or_else(|| bad("truncated ZIP structure"))?
            .try_into()
            .map_err(format_error)?,
    ))
}
fn u32_at(bytes: &[u8], offset: usize) -> Result<u32, Error> {
    Ok(u32::from_le_bytes(
        bytes
            .get(offset..offset + 4)
            .ok_or_else(|| bad("truncated ZIP structure"))?
            .try_into()
            .map_err(format_error)?,
    ))
}
#[derive(Debug)]
struct Entry {
    name: String,
    uncompressed: usize,
}
/// Strict small ZIP profile, checked before zip-rs builds metadata maps (which otherwise
/// deduplicate repeated names). No prepended executables, ZIP64, links or hidden entries.
fn zip_directory(bytes: &[u8]) -> Result<Vec<Entry>, Error> {
    if bytes.len() > MAX_ARCHIVE_BYTES {
        return Err(Error::TooLarge);
    }
    if bytes.len() < 22 || !bytes.starts_with(b"PK\x03\x04") {
        return Err(bad("not a workspace ZIP archive"));
    }
    let start = bytes.len().saturating_sub(65_557);
    let eocd = (start..=bytes.len() - 22)
        .rev()
        .find(|&i| {
            &bytes[i..i + 4] == b"PK\x05\x06"
                && u16_at(bytes, i + 20)
                    .ok()
                    .is_some_and(|n| i + 22 + n as usize == bytes.len())
        })
        .ok_or_else(|| bad("workspace ZIP footer is unavailable"))?;
    let count = u16_at(bytes, eocd + 10)? as usize;
    let size = u32_at(bytes, eocd + 12)? as usize;
    let directory = u32_at(bytes, eocd + 16)? as usize;
    if u16_at(bytes, eocd + 4)? != 0
        || u16_at(bytes, eocd + 6)? != 0
        || u16_at(bytes, eocd + 8)? as usize != count
        || count == 0
        || count > MAX_ENTRIES
        || directory.checked_add(size) != Some(eocd)
        || u16_at(bytes, eocd + 20)? != 0
    {
        return Err(bad("unsupported workspace ZIP footer/profile"));
    }
    let mut entries = Vec::new();
    let mut names = BTreeSet::new();
    let mut ranges = Vec::new();
    let mut position = directory;
    let mut expanded = 0usize;
    for _ in 0..count {
        if bytes.get(position..position + 4) != Some(b"PK\x01\x02") {
            return Err(bad("invalid workspace ZIP directory"));
        }
        let flags = u16_at(bytes, position + 8)?;
        let method = u16_at(bytes, position + 10)?;
        let crc = u32_at(bytes, position + 16)?;
        let compressed = u32_at(bytes, position + 20)? as usize;
        let uncompressed = u32_at(bytes, position + 24)? as usize;
        let name_len = u16_at(bytes, position + 28)? as usize;
        let extra = u16_at(bytes, position + 30)? as usize;
        let comment = u16_at(bytes, position + 32)? as usize;
        let offset = u32_at(bytes, position + 42)? as usize;
        let name_start = position.checked_add(46).ok_or(Error::TooLarge)?;
        let next = name_start
            .checked_add(name_len)
            .and_then(|n| n.checked_add(extra))
            .and_then(|n| n.checked_add(comment))
            .filter(|n| *n <= eocd)
            .ok_or_else(|| bad("truncated workspace ZIP directory"))?;
        let name = std::str::from_utf8(
            bytes
                .get(name_start..name_start + name_len)
                .ok_or_else(|| bad("invalid workspace entry name"))?,
        )
        .map_err(format_error)?;
        if !(name == "manifest.json" || name.strip_prefix("blobs/").is_some_and(digest))
            || !names.insert(name.to_owned())
            || flags & !0x080e != 0
            || !matches!(method, 0 | 8)
            || u16_at(bytes, position + 34)? != 0
            || (u32_at(bytes, position + 38)? >> 16) & 0o170000 == 0o120000
        {
            return Err(bad("duplicate or unsupported workspace ZIP entry"));
        }
        let cap = if name == "manifest.json" {
            MAX_MANIFEST_BYTES
        } else {
            MAX_ENTRY_BYTES
        };
        if uncompressed > cap {
            return Err(Error::TooLarge);
        }
        charge(&mut expanded, uncompressed, MAX_EXPANDED_BYTES)?;
        if offset >= directory {
            return Err(bad("workspace ZIP local offset is invalid"));
        }
        if bytes.get(offset..offset + 4) != Some(b"PK\x03\x04")
            || u16_at(bytes, offset + 6)? != flags
            || u16_at(bytes, offset + 8)? != method
        {
            return Err(bad("workspace ZIP local/directory headers disagree"));
        }
        let local_name = u16_at(bytes, offset + 26)? as usize;
        let local_extra = u16_at(bytes, offset + 28)? as usize;
        let name_begin = offset.checked_add(30).ok_or(Error::TooLarge)?;
        if local_name != name_len
            || bytes.get(name_begin..name_begin + local_name) != Some(name.as_bytes())
        {
            return Err(bad("workspace ZIP local name differs from its directory"));
        }
        let data = name_begin
            .checked_add(local_name)
            .and_then(|n| n.checked_add(local_extra))
            .ok_or(Error::TooLarge)?;
        let mut end = data
            .checked_add(compressed)
            .filter(|n| *n <= directory)
            .ok_or_else(|| bad("workspace ZIP file range is invalid"))?;
        if flags & 8 == 0 {
            if u32_at(bytes, offset + 14)? != crc
                || u32_at(bytes, offset + 18)? as usize != compressed
                || u32_at(bytes, offset + 22)? as usize != uncompressed
            {
                return Err(bad("workspace ZIP sizes or CRC disagree"));
            }
        } else {
            let local_crc = u32_at(bytes, offset + 14)?;
            let local_compressed = u32_at(bytes, offset + 18)? as usize;
            let local_uncompressed = u32_at(bytes, offset + 22)? as usize;
            if (local_crc != 0 && local_crc != crc)
                || (local_compressed != 0 && local_compressed != compressed)
                || (local_uncompressed != 0 && local_uncompressed != uncompressed)
            {
                return Err(bad("workspace ZIP local descriptor placeholders disagree"));
            }
            let descriptor = if bytes.get(end..end + 4) == Some(b"PK\x07\x08") {
                end + 4
            } else {
                end
            };
            if u32_at(bytes, descriptor)? != crc
                || u32_at(bytes, descriptor + 4)? as usize != compressed
                || u32_at(bytes, descriptor + 8)? as usize != uncompressed
            {
                return Err(bad("invalid workspace ZIP data descriptor"));
            }
            end = descriptor
                .checked_add(12)
                .filter(|n| *n <= directory)
                .ok_or_else(|| bad("workspace ZIP descriptor exceeds its entry"))?;
        }
        ranges.push((offset, end));
        entries.push(Entry {
            name: name.into(),
            uncompressed,
        });
        position = next;
    }
    if position != eocd || !names.contains("manifest.json") {
        return Err(bad("workspace ZIP directory/manifest is invalid"));
    }
    ranges.sort_unstable();
    let mut expected = 0;
    for (begin, end) in ranges {
        if begin != expected || end <= begin {
            return Err(bad(
                "workspace ZIP contains overlapping or hidden local data",
            ));
        }
        expected = end;
    }
    if expected != directory {
        return Err(bad("workspace ZIP has unreferenced data"));
    }
    Ok(entries)
}
/// Recognizes a bounded archive's explicit format marker, including future versions so
/// callers can report an unsupported version instead of misrouting it as sequence text.
pub fn is_workspace_archive(bytes: &[u8]) -> bool {
    let Ok(entries) = zip_directory(bytes) else {
        return false;
    };
    let Some(entry) = entries.iter().find(|e| e.name == "manifest.json") else {
        return false;
    };
    if entry.uncompressed > MAX_MANIFEST_BYTES {
        return false;
    }
    let Ok(mut archive) = ZipArchive::new(Cursor::new(bytes)) else {
        return false;
    };
    let Ok(mut file) = archive.by_name("manifest.json") else {
        return false;
    };
    let mut content = Vec::new();
    if file
        .by_ref()
        .take(MAX_MANIFEST_BYTES as u64 + 1)
        .read_to_end(&mut content)
        .is_err()
    {
        return false;
    }
    json::<serde_json::Value>(&content)
        .ok()
        .and_then(|v| {
            v.get("format")
                .and_then(|v| v.as_str())
                .map(|f| f == FORMAT)
        })
        .unwrap_or(false)
}
/// Load all data and validate graph/integrity before exposing any mutable application state.
pub fn load_workspace(bytes: &[u8]) -> Result<Workspace, Error> {
    let entries = zip_directory(bytes)?;
    let mut archive = ZipArchive::new(Cursor::new(bytes)).map_err(format_error)?;
    if archive.len() != entries.len()
        || archive.offset() != 0
        || archive.has_overlapping_files().map_err(format_error)?
    {
        return Err(bad(
            "workspace ZIP directory differs from its validated structure",
        ));
    }
    let mut files = BTreeMap::new();
    for (i, entry) in entries.iter().enumerate() {
        let mut file = archive.by_index(i).map_err(format_error)?;
        if file.name() != entry.name
            || file.size() != entry.uncompressed as u64
            || file.is_dir()
            || file.is_symlink()
            || file.encrypted()
            || !matches!(
                file.compression(),
                CompressionMethod::Stored | CompressionMethod::Deflated
            )
        {
            return Err(bad("unsupported workspace ZIP entry metadata"));
        }
        let mut data = Vec::new();
        file.by_ref()
            .take(entry.uncompressed as u64 + 1)
            .read_to_end(&mut data)
            .map_err(format_error)?;
        if data.len() != entry.uncompressed {
            return Err(bad("workspace entry decompression length differs"));
        }
        files.insert(entry.name.clone(), data);
    }
    let manifest_bytes = files
        .remove("manifest.json")
        .ok_or_else(|| bad("workspace manifest is missing"))?;
    let header: serde_json::Value = json(&manifest_bytes)?;
    if header.get("format").and_then(|v| v.as_str()) != Some(FORMAT) {
        return Err(bad("not a Collagene workspace manifest"));
    }
    if header.get("version").and_then(|v| v.as_u64()) != Some(VERSION as u64)
        || header.get("molecule_schema").and_then(|v| v.as_u64()) != Some(1)
        || header.get("plan_schema").and_then(|v| v.as_u64()) != Some(PLAN_SCHEMA_VERSION as u64)
    {
        return Err(Error::Unsupported(
            "workspace archive/model schema version".into(),
        ));
    }
    let manifest: Manifest = serde_json::from_value(header).map_err(format_error)?;
    if manifest.revisions.len() > MAX_REVISIONS
        || manifest.documents.len() > MAX_DOCUMENTS
        || manifest.history.len() > MAX_HISTORY_EVENTS
        || manifest.candidates.len() > MAX_CANDIDATES
        || manifest.attachments.len() > MAX_ATTACHMENTS
    {
        return Err(Error::TooLarge);
    }
    if manifest.blobs.len() != files.len() {
        return Err(bad("workspace blob manifest and archive entries disagree"));
    }
    let mut blobs = BTreeMap::new();
    for (hash, info) in &manifest.blobs {
        if !digest(hash) {
            return Err(bad("invalid workspace blob digest"));
        }
        let value = files
            .remove(&format!("blobs/{hash}"))
            .ok_or_else(|| bad("listed workspace blob is missing"))?;
        if value.len() != info.length || blake3::hash(&value).to_hex().as_str() != hash {
            return Err(bad("workspace blob hash or length differs"));
        }
        blobs.insert(hash.clone(), value);
    }
    if !files.is_empty() {
        return Err(bad("workspace archive contains unlisted entries"));
    }
    let mut revisions = BTreeMap::new();
    let mut sizes = BTreeMap::new();
    let mut expanded = 0usize;
    let mut used = BTreeSet::new();
    for (revision, blob) in &manifest.revisions {
        if !digest(revision) {
            return Err(bad("invalid workspace molecule revision"));
        }
        used.insert(blob.clone());
        let record: MoleculeRecord = json(get(&blobs, blob)?)?;
        used.extend([
            record.name.clone(),
            record.sequence.clone(),
            record.features.clone(),
        ]);
        used.extend(record.original_genbank.iter().cloned());
        used.extend(record.metadata.values().map(|value| match value {
            ValueRecord::Utf8(hash) | ValueRecord::Base64(hash) => hash.clone(),
        }));
        if record.metadata.len() > 10_000 {
            return Err(Error::TooLarge);
        }
        // Charge a conservative lower bound before duplicating large shared evidence.
        let mut lower = 0usize;
        for hash in [&record.name, &record.sequence, &record.features]
            .into_iter()
            .chain(record.original_genbank.iter())
        {
            charge(&mut lower, get(&blobs, hash)?.len(), MAX_ENTRY_BYTES)?;
        }
        for value in record.metadata.values() {
            let (hash, base64) = match value {
                ValueRecord::Utf8(h) => (h, false),
                ValueRecord::Base64(h) => (h, true),
            };
            let n = get(&blobs, hash)?.len();
            charge(
                &mut lower,
                if base64 {
                    n.checked_add(2)
                        .and_then(|n| n.checked_div(3))
                        .and_then(|n| n.checked_mul(4))
                        .ok_or(Error::TooLarge)?
                } else {
                    n
                },
                MAX_ENTRY_BYTES,
            )?;
        }
        if expanded
            .checked_add(lower)
            .is_none_or(|n| n > MAX_EXPANDED_BYTES)
        {
            return Err(Error::TooLarge);
        }
        let raw = get(&blobs, &record.sequence)?;
        let sequence = Sequence::new(
            record.alphabet,
            std::str::from_utf8(raw).map_err(format_error)?,
        )?;
        if sequence.as_str().as_bytes() != raw {
            return Err(bad("noncanonical workspace sequence chunk"));
        }
        let features: Vec<Feature> = json_exact(get(&blobs, &record.features)?)?;
        let mut metadata = BTreeMap::new();
        for (key, value) in record.metadata {
            let value = match value {
                ValueRecord::Utf8(hash) => text(get(&blobs, &hash)?)?,
                ValueRecord::Base64(hash) => STANDARD.encode(get(&blobs, &hash)?),
            };
            metadata.insert(key, value);
        }
        let molecule = Molecule {
            name: text(get(&blobs, &record.name)?)?,
            sequence,
            topology: record.topology,
            features,
            metadata,
            original_genbank: record
                .original_genbank
                .map(|hash| text(get(&blobs, &hash)?))
                .transpose()?,
        };
        molecule.validate()?;
        let size = molecule_size(&molecule)?;
        charge(&mut expanded, size, MAX_EXPANDED_BYTES)?;
        if molecule.revision() != *revision {
            return Err(bad(
                "workspace molecule revision differs from reconstructed content",
            ));
        }
        sizes.insert(revision.clone(), size);
        revisions.insert(revision.clone(), molecule);
    }
    let mut candidates = Vec::new();
    for blob in &manifest.candidates {
        used.insert(blob.clone());
        let record: PlanRecord = json(get(&blobs, blob)?)?;
        used.insert(record.steps.clone());
        if record.inputs.len() > 1000 || record.name.len() > 1024 {
            return Err(Error::TooLarge);
        }
        let before_candidate = expanded;
        charge(
            &mut expanded,
            get(&blobs, &record.steps)?.len(),
            MAX_EXPANDED_BYTES,
        )?;
        let mut inputs = Vec::new();
        for input in record.inputs {
            charge(
                &mut expanded,
                *sizes
                    .get(&input.revision)
                    .ok_or_else(|| bad("candidate source revision is missing"))?,
                MAX_EXPANDED_BYTES,
            )?;
            inputs.push(PlanInput {
                id: input.id,
                revision: input.revision.clone(),
                molecule: revisions
                    .get(&input.revision)
                    .ok_or_else(|| bad("candidate source revision is missing"))?
                    .clone(),
            });
        }
        let steps: Vec<Operation> = json_exact(get(&blobs, &record.steps)?)?;
        let candidate = Candidate {
            plan_revision: record.plan_revision,
            plan: DesignPlan {
                name: record.name,
                target_revision: record.target_revision,
                inputs,
                steps,
            },
        };
        // Bound complete reconstructed plans before revision hashing can allocate JSON.
        let plan_size = json_size(&candidate.plan, MAX_ENTRY_BYTES)?;
        expanded = before_candidate;
        charge(&mut expanded, plan_size, MAX_EXPANDED_BYTES)?;
        candidate.validate()?;
        candidates.push(candidate);
    }
    let mut attachments = BTreeMap::new();
    for (hash, record) in manifest.attachments {
        used.insert(record.blob.clone());
        let bytes = get(&blobs, &record.blob)?;
        charge(&mut expanded, bytes.len(), MAX_EXPANDED_BYTES)?;
        if record.kind == crate::run_lineage::WORKSPACE_ATTACHMENT_KIND {
            let run = crate::run_lineage::load_run_record(bytes)?;
            if manifest.plan_schema != run.plan_schema {
                return Err(Error::Unsupported(
                    "retained cloning run has an unsupported current plan schema".into(),
                ));
            }
        }
        attachments.insert(
            hash,
            WorkspaceAttachment {
                kind: record.kind,
                name: record.name,
                media_type: record.media_type,
                bytes: bytes.to_vec(),
            },
        );
    }
    if used.len() != blobs.len() {
        return Err(bad("workspace archive contains unreachable blobs"));
    }
    let workspace = Workspace {
        name: manifest.name,
        active_document: manifest.active_document,
        documents: manifest.documents,
        revisions,
        history: manifest.history,
        candidates,
        active_candidate: manifest.active_candidate,
        attachments,
    };
    workspace.validate()?;
    if workspace.unresolved_candidate_targets() != manifest.unresolved_candidate_targets {
        return Err(bad("unresolved candidate target manifest differs"));
    }
    Ok(workspace)
}
