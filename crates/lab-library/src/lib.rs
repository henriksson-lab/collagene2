//! Embedded shared-lab storage. The trusted server owns this handle and authenticates
//! membership before granting a [`LabAccess`]. No scientific or identity claim in a
//! payload grants permission or material authority.
use rusqlite::{Connection, OptionalExtension, Transaction, TransactionBehavior, params};
use serde::{Deserialize, Serialize};
use std::{path::Path, time::Duration};

pub const MAX_WORKSPACE_BYTES: usize = 64 * 1024 * 1024;
pub const MAX_NATIVE_DOCUMENT_BYTES: usize = 32 * 1024 * 1024;
pub const MAX_PAGE_SIZE: u16 = 200;
pub const MAX_FOLDER_DEPTH: usize = 128;
pub const MAX_LAB_ENTITIES: u64 = 100_000;
pub const MAX_DOCUMENT_REVISIONS: u64 = 10_000;
pub const DEFAULT_MAX_LAB_ARCHIVE_BYTES: u64 = 1024 * 1024 * 1024;
#[derive(Debug, Clone, Copy)]
pub struct StorageOptions {
    pub max_lab_archive_bytes: u64,
}
impl Default for StorageOptions {
    fn default() -> Self {
        Self {
            max_lab_archive_bytes: DEFAULT_MAX_LAB_ARCHIVE_BYTES,
        }
    }
}
impl StorageOptions {
    fn validate(self) -> Result<()> {
        if self.max_lab_archive_bytes == 0 || self.max_lab_archive_bytes > (1u64 << 46) {
            return Err(Error::Invalid("lab byte quota"));
        }
        Ok(())
    }
}
const APPLICATION_ID: i64 = 0x43474c42;
const SCHEMA_VERSION: i64 = 1;

type Result<T> = std::result::Result<T, Error>;
#[derive(Debug, thiserror::Error)]
pub enum Error {
    #[error("invalid {0}")]
    Invalid(&'static str),
    #[error("item is not present in this lab")]
    NotFound,
    #[error("an item with that name already exists in this folder")]
    AlreadyExists,
    #[error("revision conflict: expected {expected}, current {actual}")]
    Conflict {
        expected: RevisionId,
        actual: RevisionId,
    },
    #[error("metadata conflict: expected generation {expected}, current {actual}")]
    MetadataConflict { expected: u64, actual: u64 },
    #[error("catalog changed: expected generation {expected}, current {actual}")]
    CatalogConflict { expected: u64, actual: u64 },
    #[error("folder is not empty")]
    NotEmpty,
    #[error("storage limit exceeded: {0}")]
    Limit(&'static str),
    #[error("stored data failed integrity checks: {0}")]
    Corrupt(&'static str),
    #[error("unsupported or unrelated database")]
    UnsupportedDatabase,
    #[error(transparent)]
    Database(#[from] rusqlite::Error),
}
macro_rules! identifier {
    ($name:ident, $prefix:literal, $digits:expr) => {
        #[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize)]
        #[serde(try_from = "String", into = "String")]
        pub struct $name(String);
        impl $name {
            pub fn parse(value: impl Into<String>) -> Result<Self> {
                let value = value.into();
                let suffix = value
                    .strip_prefix($prefix)
                    .ok_or(Error::Invalid(stringify!($name)))?;
                if suffix.len() != $digits
                    || !suffix
                        .bytes()
                        .all(|b| b.is_ascii_digit() || (b'a'..=b'f').contains(&b))
                {
                    return Err(Error::Invalid(stringify!($name)));
                }
                Ok(Self(value))
            }
            pub fn as_str(&self) -> &str {
                &self.0
            }
        }
        impl std::fmt::Display for $name {
            fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
                self.0.fmt(f)
            }
        }
        impl TryFrom<String> for $name {
            type Error = Error;
            fn try_from(value: String) -> Result<Self> {
                Self::parse(value)
            }
        }
        impl From<$name> for String {
            fn from(value: $name) -> String {
                value.0
            }
        }
    };
}
identifier!(LabId, "lab-", 32);
identifier!(FolderId, "fld-", 32);
identifier!(DocumentId, "doc-", 32);
identifier!(RevisionId, "", 64);
identifier!(Digest, "", 64);

/// A server-side scope, deliberately neither serializable nor client-constructible.
#[derive(Debug, Clone)]
pub struct LabAccess {
    lab: LabId,
}
impl LabAccess {
    pub fn lab_id(&self) -> &LabId {
        &self.lab
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum PayloadKind {
    NativeDocument,
    WorkspaceArchive,
}
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum WorkspaceSharing {
    DocumentClosure,
    ExplicitWholeWorkspace,
}
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct NativeSelection {
    pub document_id: String,
    pub revision: Digest,
    pub head_event: Option<String>,
}
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Payload {
    pub kind: PayloadKind,
    pub bytes: Vec<u8>,
    pub archive_hash: Digest,
    pub selection: Option<NativeSelection>,
    pub sharing: Option<WorkspaceSharing>,
}
impl Payload {
    /// Shape/integrity only. The host MUST decode/validate native data and sharing
    /// boundaries before passing a payload to a write operation.
    pub fn validate(&self) -> Result<()> {
        let maximum = match self.kind {
            PayloadKind::NativeDocument => MAX_NATIVE_DOCUMENT_BYTES,
            PayloadKind::WorkspaceArchive => MAX_WORKSPACE_BYTES,
        };
        if self.bytes.is_empty() || self.bytes.len() > maximum {
            return Err(Error::Limit("payload bytes"));
        }
        if blake3::hash(&self.bytes).to_hex().as_str() != self.archive_hash.as_str() {
            return Err(Error::Invalid("archive hash"));
        }
        match self.kind {
            PayloadKind::NativeDocument if self.selection.is_none() && self.sharing.is_none() => (),
            PayloadKind::WorkspaceArchive if self.selection.is_some() && self.sharing.is_some() => {
                let selection = self.selection.as_ref().unwrap();
                native_id(&selection.document_id)?;
                if let Some(event) = &selection.head_event {
                    native_id(event)?;
                }
            }
            _ => return Err(Error::Invalid("payload selection/sharing")),
        }
        Ok(())
    }
}
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct Lab {
    pub id: LabId,
    pub name: String,
    pub catalog_generation: u64,
}
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct Folder {
    pub id: FolderId,
    pub parent: Option<FolderId>,
    pub name: String,
    pub generation: u64,
}
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct RevisionSummary {
    pub id: RevisionId,
    pub parent: Option<RevisionId>,
    pub number: u64,
    pub archive_hash: Digest,
    pub archive_size: u64,
    pub kind: PayloadKind,
    pub selection: Option<NativeSelection>,
    pub sharing: Option<WorkspaceSharing>,
    pub created_at: u64,
}
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct Document {
    pub id: DocumentId,
    pub folder: Option<FolderId>,
    pub name: String,
    pub generation: u64,
    pub head: RevisionSummary,
}
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub enum Entry {
    Folder(Folder),
    Document(Document),
}
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct Page {
    pub catalog_generation: u64,
    pub entries: Vec<Entry>,
    pub next_offset: Option<u32>,
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub enum SearchScope {
    AllFolders,
    Folder {
        folder: Option<FolderId>,
        recursive: bool,
    },
}
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct RevisionPage {
    pub document: DocumentId,
    pub catalog_generation: u64,
    pub revisions: Vec<RevisionSummary>,
    pub next_offset: Option<u32>,
}

/// One native SQLite connection; open additional handles for concurrent requests.
/// SQLite serializes writes; CAS and tree validation run under BEGIN IMMEDIATE.
pub struct Library {
    connection: Connection,
    options: StorageOptions,
}
impl Library {
    pub fn open(path: impl AsRef<Path>) -> Result<Self> {
        Self::open_with_options(path, StorageOptions::default())
    }
    pub fn open_in_memory() -> Result<Self> {
        Self::open_in_memory_with_options(StorageOptions::default())
    }
    pub fn open_with_options(path: impl AsRef<Path>, options: StorageOptions) -> Result<Self> {
        options.validate()?;
        Self::configure(Connection::open(path)?, options)
    }
    pub fn open_in_memory_with_options(options: StorageOptions) -> Result<Self> {
        options.validate()?;
        Self::configure(Connection::open_in_memory()?, options)
    }
    fn configure(mut connection: Connection, options: StorageOptions) -> Result<Self> {
        connection.busy_timeout(Duration::from_secs(5))?;
        connection.execute_batch("PRAGMA foreign_keys=ON; PRAGMA trusted_schema=OFF;")?;
        let tx = connection.transaction_with_behavior(TransactionBehavior::Immediate)?;
        let version: i64 = tx.query_row("PRAGMA user_version", [], |r| r.get(0))?;
        let application: i64 = tx.query_row("PRAGMA application_id", [], |r| r.get(0))?;
        if version == 0 && application == 0 {
            let tables: i64 = tx.query_row(
                "SELECT count(*) FROM sqlite_master WHERE name NOT LIKE 'sqlite_%'",
                [],
                |r| r.get(0),
            )?;
            if tables != 0 {
                return Err(Error::UnsupportedDatabase);
            }
            tx.execute_batch(SCHEMA)?;
            tx.pragma_update(None, "application_id", APPLICATION_ID)?;
            tx.pragma_update(None, "user_version", SCHEMA_VERSION)?;
        } else if version != SCHEMA_VERSION || application != APPLICATION_ID {
            return Err(Error::UnsupportedDatabase);
        }
        tx.commit()?;
        connection.execute_batch("PRAGMA journal_mode=WAL; PRAGMA synchronous=FULL;")?;
        Ok(Self {
            connection,
            options,
        })
    }
    /// Host administration only; do not expose this as unauthenticated lab creation.
    pub fn create_lab(&mut self, name: &str) -> Result<Lab> {
        valid_name(name)?;
        let tx = self
            .connection
            .transaction_with_behavior(TransactionBehavior::Immediate)?;
        let id = LabId::parse(random_id(&tx, "lab-")?)?;
        tx.execute(
            "INSERT INTO labs(id,name) VALUES(?1,?2)",
            params![id.as_str(), name],
        )?;
        tx.commit()?;
        Ok(Lab {
            id,
            name: name.into(),
            catalog_generation: 0,
        })
    }
    /// TRUSTED HOST ONLY: invoke after authenticating current membership in this lab.
    /// A request-supplied lab/member identifier must never be treated as authorization.
    pub fn trusted_host_access(&self, lab: &LabId) -> Result<LabAccess> {
        ensure_lab(&self.connection, lab)?;
        Ok(LabAccess { lab: lab.clone() })
    }
    pub fn lab(&self, access: &LabAccess) -> Result<Lab> {
        let found = self
            .connection
            .query_row(
                "SELECT name,catalog_generation FROM labs WHERE id=?1",
                [access.lab.as_str()],
                |r| Ok((r.get(0)?, r.get(1)?)),
            )
            .optional()?;
        let (name, generation): (String, u64) = found.ok_or(Error::NotFound)?;
        Ok(Lab {
            id: access.lab.clone(),
            name,
            catalog_generation: generation,
        })
    }
    pub fn folder(&self, access: &LabAccess, id: &FolderId) -> Result<Folder> {
        get_folder(&self.connection, &access.lab, id)
    }
    pub fn document(&self, access: &LabAccess, id: &DocumentId) -> Result<Document> {
        get_document(&self.connection, &access.lab, id)
    }
    pub fn create_folder(
        &mut self,
        access: &LabAccess,
        parent: Option<&FolderId>,
        name: &str,
    ) -> Result<Folder> {
        valid_name(name)?;
        let tx = self
            .connection
            .transaction_with_behavior(TransactionBehavior::Immediate)?;
        ensure_lab(&tx, &access.lab)?;
        if parent_depth(&tx, &access.lab, parent, None)? + 1 > MAX_FOLDER_DEPTH {
            return Err(Error::Limit("folder depth"));
        }
        entity_capacity(&tx, &access.lab)?;
        unique_name(&tx, &access.lab, parent, name, None, None)?;
        let id = FolderId::parse(random_id(&tx, "fld-")?)?;
        tx.execute(
            "INSERT INTO folders(id,lab_id,parent_id,name) VALUES(?1,?2,?3,?4)",
            params![
                id.as_str(),
                access.lab.as_str(),
                parent.map(FolderId::as_str),
                name
            ],
        )?;
        bump_catalog(&tx, &access.lab)?;
        let folder = get_folder(&tx, &access.lab, &id)?;
        tx.commit()?;
        Ok(folder)
    }
    pub fn rename_folder(
        &mut self,
        access: &LabAccess,
        id: &FolderId,
        expected_generation: u64,
        name: &str,
    ) -> Result<Folder> {
        valid_name(name)?;
        let tx = self
            .connection
            .transaction_with_behavior(TransactionBehavior::Immediate)?;
        let folder = get_folder(&tx, &access.lab, id)?;
        metadata_cas(expected_generation, folder.generation)?;
        unique_name(
            &tx,
            &access.lab,
            folder.parent.as_ref(),
            name,
            Some(id),
            None,
        )?;
        tx.execute(
            "UPDATE folders SET name=?1,generation=generation+1 WHERE id=?2 AND lab_id=?3",
            params![name, id.as_str(), access.lab.as_str()],
        )?;
        bump_catalog(&tx, &access.lab)?;
        let folder = get_folder(&tx, &access.lab, id)?;
        tx.commit()?;
        Ok(folder)
    }
    pub fn move_folder(
        &mut self,
        access: &LabAccess,
        id: &FolderId,
        expected_generation: u64,
        parent: Option<&FolderId>,
    ) -> Result<Folder> {
        let tx = self
            .connection
            .transaction_with_behavior(TransactionBehavior::Immediate)?;
        let folder = get_folder(&tx, &access.lab, id)?;
        metadata_cas(expected_generation, folder.generation)?;
        let depth = parent_depth(&tx, &access.lab, parent, Some(id))?;
        let height:usize=tx.query_row("WITH RECURSIVE branch(id,depth) AS (SELECT id,1 FROM folders WHERE id=?1 AND lab_id=?2 UNION ALL SELECT f.id,b.depth+1 FROM folders f JOIN branch b ON f.parent_id=b.id WHERE f.lab_id=?2 AND b.depth<=?3) SELECT max(depth) FROM branch",params![id.as_str(),access.lab.as_str(),MAX_FOLDER_DEPTH as u64],|r|r.get(0))?;
        if depth + height > MAX_FOLDER_DEPTH {
            return Err(Error::Limit("folder depth"));
        }
        unique_name(&tx, &access.lab, parent, &folder.name, Some(id), None)?;
        tx.execute(
            "UPDATE folders SET parent_id=?1,generation=generation+1 WHERE id=?2 AND lab_id=?3",
            params![
                parent.map(FolderId::as_str),
                id.as_str(),
                access.lab.as_str()
            ],
        )?;
        bump_catalog(&tx, &access.lab)?;
        let folder = get_folder(&tx, &access.lab, id)?;
        tx.commit()?;
        Ok(folder)
    }
    /// There is deliberately no recursive deletion or document/revision deletion API.
    pub fn delete_empty_folder(
        &mut self,
        access: &LabAccess,
        id: &FolderId,
        expected_generation: u64,
    ) -> Result<()> {
        let tx = self
            .connection
            .transaction_with_behavior(TransactionBehavior::Immediate)?;
        let folder = get_folder(&tx, &access.lab, id)?;
        metadata_cas(expected_generation, folder.generation)?;
        let children:bool=tx.query_row("SELECT EXISTS(SELECT 1 FROM folders WHERE lab_id=?1 AND parent_id=?2 UNION ALL SELECT 1 FROM documents WHERE lab_id=?1 AND folder_id=?2)",params![access.lab.as_str(),id.as_str()],|r|r.get(0))?;
        if children {
            return Err(Error::NotEmpty);
        }
        tx.execute(
            "DELETE FROM folders WHERE lab_id=?1 AND id=?2",
            params![access.lab.as_str(), id.as_str()],
        )?;
        bump_catalog(&tx, &access.lab)?;
        tx.commit()?;
        Ok(())
    }
    pub fn create_document(
        &mut self,
        access: &LabAccess,
        folder: Option<&FolderId>,
        name: &str,
        payload: &Payload,
    ) -> Result<Document> {
        valid_name(name)?;
        payload.validate()?;
        let tx = self
            .connection
            .transaction_with_behavior(TransactionBehavior::Immediate)?;
        ensure_lab(&tx, &access.lab)?;
        parent_depth(&tx, &access.lab, folder, None)?;
        entity_capacity(&tx, &access.lab)?;
        unique_name(&tx, &access.lab, folder, name, None, None)?;
        let id = DocumentId::parse(random_id(&tx, "doc-")?)?;
        charge_archive(
            &tx,
            &access.lab,
            payload.bytes.len(),
            self.options.max_lab_archive_bytes,
        )?;
        let revision = revision_id(&id, None, payload)?;
        tx.execute(
            "INSERT INTO documents(id,lab_id,folder_id,name,head_revision) VALUES(?1,?2,?3,?4,?5)",
            params![
                id.as_str(),
                access.lab.as_str(),
                folder.map(FolderId::as_str),
                name,
                revision.as_str()
            ],
        )?;
        insert_revision(&tx, &id, &revision, None, 1, payload)?;
        bump_catalog(&tx, &access.lab)?;
        let document = get_document(&tx, &access.lab, &id)?;
        tx.commit()?;
        Ok(document)
    }
    /// Content-only save. Stale writers receive the current head; no bytes/history
    /// are changed on conflict. Metadata moves/renames use independent generation CAS.
    pub fn save(
        &mut self,
        access: &LabAccess,
        id: &DocumentId,
        base_revision: &RevisionId,
        payload: &Payload,
    ) -> Result<Document> {
        payload.validate()?;
        let tx = self
            .connection
            .transaction_with_behavior(TransactionBehavior::Immediate)?;
        let document = get_document(&tx, &access.lab, id)?;
        if &document.head.id != base_revision {
            return Err(Error::Conflict {
                expected: base_revision.clone(),
                actual: document.head.id,
            });
        }
        if document.head.number >= MAX_DOCUMENT_REVISIONS {
            return Err(Error::Limit("document revisions"));
        }
        charge_archive(
            &tx,
            &access.lab,
            payload.bytes.len(),
            self.options.max_lab_archive_bytes,
        )?;
        let revision = revision_id(id, Some(base_revision), payload)?;
        insert_revision(
            &tx,
            id,
            &revision,
            Some(base_revision),
            document.head.number + 1,
            payload,
        )?;
        let changed=tx.execute("UPDATE documents SET head_revision=?1,generation=generation+1 WHERE id=?2 AND lab_id=?3 AND head_revision=?4",params![revision.as_str(),id.as_str(),access.lab.as_str(),base_revision.as_str()])?;
        if changed != 1 {
            return Err(Error::Corrupt("CAS update"));
        }
        bump_catalog(&tx, &access.lab)?;
        let document = get_document(&tx, &access.lab, id)?;
        tx.commit()?;
        Ok(document)
    }
    pub fn rename_document(
        &mut self,
        access: &LabAccess,
        id: &DocumentId,
        expected_generation: u64,
        name: &str,
    ) -> Result<Document> {
        valid_name(name)?;
        self.change_document_metadata(access, id, expected_generation, Some(name), None)
    }
    pub fn move_document(
        &mut self,
        access: &LabAccess,
        id: &DocumentId,
        expected_generation: u64,
        folder: Option<&FolderId>,
    ) -> Result<Document> {
        self.change_document_metadata(access, id, expected_generation, None, Some(folder))
    }
    fn change_document_metadata(
        &mut self,
        access: &LabAccess,
        id: &DocumentId,
        expected_generation: u64,
        name: Option<&str>,
        folder: Option<Option<&FolderId>>,
    ) -> Result<Document> {
        let tx = self
            .connection
            .transaction_with_behavior(TransactionBehavior::Immediate)?;
        let document = get_document(&tx, &access.lab, id)?;
        metadata_cas(expected_generation, document.generation)?;
        let name = name.unwrap_or(&document.name);
        let folder = folder.unwrap_or(document.folder.as_ref());
        parent_depth(&tx, &access.lab, folder, None)?;
        unique_name(&tx, &access.lab, folder, name, None, Some(id))?;
        tx.execute("UPDATE documents SET name=?1,folder_id=?2,generation=generation+1 WHERE id=?3 AND lab_id=?4",params![name,folder.map(FolderId::as_str),id.as_str(),access.lab.as_str()])?;
        bump_catalog(&tx, &access.lab)?;
        let document = get_document(&tx, &access.lab, id)?;
        tx.commit()?;
        Ok(document)
    }
    /// Snapshot-stable pagination: require the returned catalog generation on
    /// subsequent pages. Any intervening mutation produces an explicit conflict.
    pub fn browse(
        &mut self,
        access: &LabAccess,
        parent: Option<&FolderId>,
        offset: u32,
        limit: u16,
        expected_catalog_generation: Option<u64>,
    ) -> Result<Page> {
        pagination(offset, limit)?;
        continuation(offset, expected_catalog_generation)?;
        let tx = self.connection.transaction()?;
        parent_depth(&tx, &access.lab, parent, None)?;
        let generation = catalog_generation(&tx, &access.lab)?;
        if let Some(expected) = expected_catalog_generation {
            if expected != generation {
                return Err(Error::CatalogConflict {
                    expected,
                    actual: generation,
                });
            }
        }
        let mut statement=tx.prepare("SELECT kind,id FROM (SELECT 0 AS kind,id,name FROM folders WHERE lab_id=?1 AND parent_id IS ?2 UNION ALL SELECT 1 AS kind,id,name FROM documents WHERE lab_id=?1 AND folder_id IS ?2) ORDER BY name COLLATE BINARY,id LIMIT ?3 OFFSET ?4")?;
        let ids = statement
            .query_map(
                params![
                    access.lab.as_str(),
                    parent.map(FolderId::as_str),
                    u64::from(limit) + 1,
                    offset
                ],
                |r| Ok((r.get::<_, u8>(0)?, r.get::<_, String>(1)?)),
            )?
            .collect::<std::result::Result<Vec<_>, _>>()?;
        let more = ids.len() > usize::from(limit);
        let mut entries = Vec::new();
        for (kind, id) in ids.into_iter().take(usize::from(limit)) {
            entries.push(if kind == 0 {
                Entry::Folder(get_folder(&tx, &access.lab, &FolderId::parse(id)?)?)
            } else {
                Entry::Document(get_document(&tx, &access.lab, &DocumentId::parse(id)?)?)
            });
        }
        drop(statement);
        tx.commit()?;
        Ok(Page {
            catalog_generation: generation,
            entries,
            next_offset: more.then_some(offset + u32::from(limit)),
        })
    }
    /// Name substring search, with ASCII case folding from SQLite lower(). Unicode
    /// text remains intact; non-ASCII case equivalence is not invented.
    pub fn search(
        &mut self,
        access: &LabAccess,
        query: &str,
        scope: &SearchScope,
        offset: u32,
        limit: u16,
        expected_catalog_generation: Option<u64>,
    ) -> Result<Page> {
        if query.trim().is_empty() || query.len() > 256 || query.chars().any(char::is_control) {
            return Err(Error::Invalid("search query"));
        }
        pagination(offset, limit)?;
        continuation(offset, expected_catalog_generation)?;
        let tx = self.connection.transaction()?;
        let (mode, parent) = match scope {
            SearchScope::AllFolders => (0, None),
            SearchScope::Folder { folder, recursive } => {
                (if *recursive { 2 } else { 1 }, folder.as_ref())
            }
        };
        parent_depth(&tx, &access.lab, parent, None)?;
        let generation = catalog_generation(&tx, &access.lab)?;
        snapshot(expected_catalog_generation, generation)?;
        let mut statement = tx.prepare("WITH RECURSIVE scope(id,depth) AS (SELECT coalesce(?2,''),0 UNION ALL SELECT f.id,s.depth+1 FROM folders f JOIN scope s ON coalesce(f.parent_id,'')=s.id WHERE f.lab_id=?1 AND s.depth<?7), entries(kind,id,name,parent) AS (SELECT 0,id,name,parent_id FROM folders WHERE lab_id=?1 UNION ALL SELECT 1,id,name,folder_id FROM documents WHERE lab_id=?1) SELECT kind,id FROM entries WHERE instr(lower(name),lower(?3))>0 AND (?4=0 OR (?4=1 AND parent IS ?2) OR (?4=2 AND coalesce(parent,'') IN (SELECT id FROM scope))) ORDER BY name COLLATE BINARY,id LIMIT ?5 OFFSET ?6")?;
        let ids = statement
            .query_map(
                params![
                    access.lab.as_str(),
                    parent.map(FolderId::as_str),
                    query,
                    mode,
                    u64::from(limit) + 1,
                    offset,
                    MAX_FOLDER_DEPTH as u64
                ],
                |r| Ok((r.get::<_, u8>(0)?, r.get::<_, String>(1)?)),
            )?
            .collect::<std::result::Result<Vec<_>, _>>()?;
        let more = ids.len() > usize::from(limit);
        let mut entries = Vec::new();
        for (kind, id) in ids.into_iter().take(usize::from(limit)) {
            entries.push(if kind == 0 {
                Entry::Folder(get_folder(&tx, &access.lab, &FolderId::parse(id)?)?)
            } else {
                Entry::Document(get_document(&tx, &access.lab, &DocumentId::parse(id)?)?)
            });
        }
        drop(statement);
        tx.commit()?;
        Ok(Page {
            catalog_generation: generation,
            entries,
            next_offset: more.then_some(offset + u32::from(limit)),
        })
    }
    pub fn revision_page(
        &mut self,
        access: &LabAccess,
        id: &DocumentId,
        offset: u32,
        limit: u16,
        expected_catalog_generation: Option<u64>,
    ) -> Result<RevisionPage> {
        pagination(offset, limit)?;
        continuation(offset, expected_catalog_generation)?;
        let tx = self.connection.transaction()?;
        get_document(&tx, &access.lab, id)?;
        let generation = catalog_generation(&tx, &access.lab)?;
        snapshot(expected_catalog_generation, generation)?;
        let mut statement = tx.prepare(
            "SELECT id FROM revisions WHERE document_id=?1 ORDER BY number DESC LIMIT ?2 OFFSET ?3",
        )?;
        let ids = statement
            .query_map(params![id.as_str(), u64::from(limit) + 1, offset], |r| {
                r.get::<_, String>(0)
            })?
            .collect::<std::result::Result<Vec<_>, _>>()?;
        let more = ids.len() > usize::from(limit);
        let revisions = ids
            .into_iter()
            .take(usize::from(limit))
            .map(|revision| get_revision(&tx, id, &RevisionId::parse(revision)?))
            .collect::<Result<Vec<_>>>()?;
        drop(statement);
        tx.commit()?;
        Ok(RevisionPage {
            document: id.clone(),
            catalog_generation: generation,
            revisions,
            next_offset: more.then_some(offset + u32::from(limit)),
        })
    }
    pub fn revisions(
        &self,
        access: &LabAccess,
        id: &DocumentId,
        offset: u32,
        limit: u16,
    ) -> Result<Vec<RevisionSummary>> {
        pagination(offset, limit)?;
        get_document(&self.connection, &access.lab, id)?;
        let mut statement = self.connection.prepare(
            "SELECT id FROM revisions WHERE document_id=?1 ORDER BY number DESC LIMIT ?2 OFFSET ?3",
        )?;
        let ids = statement
            .query_map(params![id.as_str(), limit, offset], |r| {
                r.get::<_, String>(0)
            })?
            .collect::<std::result::Result<Vec<_>, _>>()?;
        ids.into_iter()
            .map(|revision| get_revision(&self.connection, id, &RevisionId::parse(revision)?))
            .collect()
    }
    pub fn revision(
        &self,
        access: &LabAccess,
        id: &DocumentId,
        revision: &RevisionId,
    ) -> Result<RevisionSummary> {
        get_document(&self.connection, &access.lab, id)?;
        get_revision(&self.connection, id, revision)
    }
    pub fn payload(
        &self,
        access: &LabAccess,
        id: &DocumentId,
        revision: &RevisionId,
    ) -> Result<Payload> {
        let metadata = self.revision(access, id, revision)?;
        let bytes:Option<Vec<u8>>=self.connection.query_row("SELECT CASE WHEN length(archive)<=?3 THEN archive ELSE NULL END FROM revisions WHERE document_id=?1 AND id=?2",params![id.as_str(),revision.as_str(),MAX_WORKSPACE_BYTES as u64],|r|r.get(0))?;
        let payload = Payload {
            kind: metadata.kind,
            bytes: bytes.ok_or(Error::Corrupt("oversized archive"))?,
            archive_hash: metadata.archive_hash,
            selection: metadata.selection,
            sharing: metadata.sharing,
        };
        payload
            .validate()
            .map_err(|_| Error::Corrupt("archive integrity"))?;
        if payload.bytes.len() as u64 != metadata.archive_size {
            return Err(Error::Corrupt("archive size"));
        }
        if revision_id(id, metadata.parent.as_ref(), &payload)? != metadata.id {
            return Err(Error::Corrupt("revision binding"));
        }
        Ok(payload)
    }
}

#[derive(Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
struct PayloadMetadata {
    kind: PayloadKind,
    selection: Option<NativeSelection>,
    sharing: Option<WorkspaceSharing>,
}
fn payload_metadata(payload: &Payload) -> PayloadMetadata {
    PayloadMetadata {
        kind: payload.kind,
        selection: payload.selection.clone(),
        sharing: payload.sharing,
    }
}
fn metadata_bytes(payload: &Payload) -> Result<Vec<u8>> {
    serde_json::to_vec(&payload_metadata(payload)).map_err(|_| Error::Invalid("payload metadata"))
}
fn valid_name(name: &str) -> Result<()> {
    if name.trim().is_empty()
        || name.len() > 1024
        || name.chars().count() > 256
        || name.chars().any(char::is_control)
    {
        return Err(Error::Invalid("name"));
    }
    Ok(())
}
fn native_id(id: &str) -> Result<()> {
    if id.is_empty() || id.len() > 128 || id.chars().any(|c| c.is_control() || c.is_whitespace()) {
        return Err(Error::Invalid("native selection ID"));
    }
    Ok(())
}
fn pagination(offset: u32, limit: u16) -> Result<()> {
    if limit == 0 || limit > MAX_PAGE_SIZE || u64::from(offset) > MAX_LAB_ENTITIES {
        return Err(Error::Invalid("pagination"));
    }
    Ok(())
}
fn continuation(offset: u32, generation: Option<u64>) -> Result<()> {
    if offset > 0 && generation.is_none() {
        return Err(Error::Invalid("continuation snapshot"));
    }
    Ok(())
}
fn snapshot(expected: Option<u64>, actual: u64) -> Result<()> {
    if let Some(expected) = expected {
        if expected != actual {
            return Err(Error::CatalogConflict { expected, actual });
        }
    }
    Ok(())
}
fn metadata_cas(expected: u64, actual: u64) -> Result<()> {
    if expected != actual {
        return Err(Error::MetadataConflict { expected, actual });
    }
    Ok(())
}
fn ensure_lab(connection: &Connection, lab: &LabId) -> Result<()> {
    if !connection.query_row(
        "SELECT EXISTS(SELECT 1 FROM labs WHERE id=?1)",
        [lab.as_str()],
        |r| r.get::<_, bool>(0),
    )? {
        return Err(Error::NotFound);
    }
    Ok(())
}
fn catalog_generation(connection: &Connection, lab: &LabId) -> Result<u64> {
    connection
        .query_row(
            "SELECT catalog_generation FROM labs WHERE id=?1",
            [lab.as_str()],
            |r| r.get(0),
        )
        .optional()?
        .ok_or(Error::NotFound)
}
fn bump_catalog(tx: &Transaction<'_>, lab: &LabId) -> Result<()> {
    tx.execute(
        "UPDATE labs SET catalog_generation=catalog_generation+1 WHERE id=?1",
        [lab.as_str()],
    )?;
    Ok(())
}
fn charge_archive(tx: &Transaction<'_>, lab: &LabId, bytes: usize, quota: u64) -> Result<()> {
    let current: u64 = tx
        .query_row(
            "SELECT archive_bytes FROM labs WHERE id=?1",
            [lab.as_str()],
            |r| r.get(0),
        )
        .optional()?
        .ok_or(Error::NotFound)?;
    let total = current
        .checked_add(bytes as u64)
        .ok_or(Error::Limit("lab archive bytes"))?;
    if total > quota {
        return Err(Error::Limit("lab archive bytes"));
    }
    tx.execute(
        "UPDATE labs SET archive_bytes=?1 WHERE id=?2",
        params![total, lab.as_str()],
    )?;
    Ok(())
}
fn random_id(tx: &Transaction<'_>, prefix: &str) -> Result<String> {
    let suffix: String = tx.query_row("SELECT lower(hex(randomblob(16)))", [], |r| r.get(0))?;
    Ok(format!("{prefix}{suffix}"))
}
fn entity_capacity(connection: &Connection, lab: &LabId) -> Result<()> {
    let count:u64=connection.query_row("SELECT (SELECT count(*) FROM folders WHERE lab_id=?1)+(SELECT count(*) FROM documents WHERE lab_id=?1)",[lab.as_str()],|r|r.get(0))?;
    if count >= MAX_LAB_ENTITIES {
        return Err(Error::Limit("lab entities"));
    }
    Ok(())
}
fn parent_depth(
    connection: &Connection,
    lab: &LabId,
    parent: Option<&FolderId>,
    exclude: Option<&FolderId>,
) -> Result<usize> {
    let mut current = parent.cloned();
    let mut depth = 0;
    while let Some(id) = current {
        if exclude == Some(&id) {
            return Err(Error::Invalid("folder cycle"));
        }
        depth += 1;
        if depth > MAX_FOLDER_DEPTH {
            return Err(Error::Limit("folder depth"));
        }
        current = get_folder(connection, lab, &id)?.parent;
    }
    Ok(depth)
}
fn unique_name(
    connection: &Connection,
    lab: &LabId,
    parent: Option<&FolderId>,
    name: &str,
    exclude_folder: Option<&FolderId>,
    exclude_document: Option<&DocumentId>,
) -> Result<()> {
    let exists:bool=connection.query_row("SELECT EXISTS(SELECT 1 FROM folders WHERE lab_id=?1 AND parent_id IS ?2 AND name=?3 AND id IS NOT ?4 UNION ALL SELECT 1 FROM documents WHERE lab_id=?1 AND folder_id IS ?2 AND name=?3 AND id IS NOT ?5)",params![lab.as_str(),parent.map(FolderId::as_str),name,exclude_folder.map(FolderId::as_str),exclude_document.map(DocumentId::as_str)],|r|r.get(0))?;
    if exists {
        return Err(Error::AlreadyExists);
    }
    Ok(())
}
fn get_folder(connection: &Connection, lab: &LabId, id: &FolderId) -> Result<Folder> {
    let row:Option<(Option<String>,String,u64)>=connection.query_row("SELECT parent_id,substr(name,1,1025),generation FROM folders WHERE lab_id=?1 AND id=?2",params![lab.as_str(),id.as_str()],|r|Ok((r.get(0)?,r.get(1)?,r.get(2)?))).optional()?;
    let (parent, name, generation) = row.ok_or(Error::NotFound)?;
    valid_name(&name).map_err(|_| Error::Corrupt("folder name"))?;
    Ok(Folder {
        id: id.clone(),
        parent: parent.map(FolderId::parse).transpose()?,
        name,
        generation,
    })
}
fn get_document(connection: &Connection, lab: &LabId, id: &DocumentId) -> Result<Document> {
    let row:Option<(Option<String>,String,u64,String)>=connection.query_row("SELECT folder_id,substr(name,1,1025),generation,head_revision FROM documents WHERE lab_id=?1 AND id=?2",params![lab.as_str(),id.as_str()],|r|Ok((r.get(0)?,r.get(1)?,r.get(2)?,r.get(3)?))).optional()?;
    let (folder, name, generation, head) = row.ok_or(Error::NotFound)?;
    valid_name(&name).map_err(|_| Error::Corrupt("document name"))?;
    Ok(Document {
        id: id.clone(),
        folder: folder.map(FolderId::parse).transpose()?,
        name,
        generation,
        head: get_revision(connection, id, &RevisionId::parse(head)?)?,
    })
}
fn get_revision(
    connection: &Connection,
    id: &DocumentId,
    revision: &RevisionId,
) -> Result<RevisionSummary> {
    type Row = (Option<String>, u64, String, u64, Option<String>, u64);
    let row:Option<Row>=connection.query_row("SELECT parent_id,number,archive_hash,archive_size,CASE WHEN length(metadata)<=4096 THEN metadata ELSE NULL END,created_at FROM revisions WHERE document_id=?1 AND id=?2",params![id.as_str(),revision.as_str()],|r|Ok((r.get(0)?,r.get(1)?,r.get(2)?,r.get(3)?,r.get(4)?,r.get(5)?))).optional()?;
    let (parent, number, hash, size, metadata, created_at) = row.ok_or(Error::NotFound)?;
    let metadata: PayloadMetadata =
        serde_json::from_str(&metadata.ok_or(Error::Corrupt("oversized metadata"))?)
            .map_err(|_| Error::Corrupt("payload metadata"))?;
    let max = match metadata.kind {
        PayloadKind::NativeDocument => MAX_NATIVE_DOCUMENT_BYTES,
        PayloadKind::WorkspaceArchive => MAX_WORKSPACE_BYTES,
    };
    if number == 0 || number > MAX_DOCUMENT_REVISIONS || size == 0 || size > max as u64 {
        return Err(Error::Corrupt("revision bounds"));
    }
    Ok(RevisionSummary {
        id: revision.clone(),
        parent: parent.map(RevisionId::parse).transpose()?,
        number,
        archive_hash: Digest::parse(hash)?,
        archive_size: size,
        kind: metadata.kind,
        selection: metadata.selection,
        sharing: metadata.sharing,
        created_at,
    })
}
fn revision_id(
    id: &DocumentId,
    parent: Option<&RevisionId>,
    payload: &Payload,
) -> Result<RevisionId> {
    let mut hash = blake3::Hasher::new();
    hash.update(b"collagene-lab-revision-v1\0");
    for part in [
        id.as_str().as_bytes(),
        parent.map(RevisionId::as_str).unwrap_or("").as_bytes(),
        payload.archive_hash.as_str().as_bytes(),
        &metadata_bytes(payload)?,
    ] {
        hash.update(&(part.len() as u64).to_le_bytes());
        hash.update(part);
    }
    RevisionId::parse(hash.finalize().to_hex().to_string())
}
fn insert_revision(
    tx: &Transaction<'_>,
    id: &DocumentId,
    revision: &RevisionId,
    parent: Option<&RevisionId>,
    number: u64,
    payload: &Payload,
) -> Result<()> {
    let metadata = String::from_utf8(metadata_bytes(payload)?)
        .map_err(|_| Error::Invalid("metadata UTF-8"))?;
    tx.execute("INSERT INTO revisions(id,document_id,parent_id,number,archive,archive_hash,archive_size,metadata) VALUES(?1,?2,?3,?4,?5,?6,?7,?8)",params![revision.as_str(),id.as_str(),parent.map(RevisionId::as_str),number,&payload.bytes,payload.archive_hash.as_str(),payload.bytes.len() as u64,metadata])?;
    Ok(())
}
const SCHEMA: &str = r#"
CREATE TABLE labs(
 id TEXT PRIMARY KEY NOT NULL,
 name TEXT NOT NULL CHECK(length(name) BETWEEN 1 AND 1024),
 catalog_generation INTEGER NOT NULL DEFAULT 0 CHECK(typeof(catalog_generation)='integer' AND catalog_generation>=0),
 archive_bytes INTEGER NOT NULL DEFAULT 0 CHECK(typeof(archive_bytes)='integer' AND archive_bytes>=0)
);
CREATE TABLE folders(
 id TEXT PRIMARY KEY NOT NULL,
 lab_id TEXT NOT NULL REFERENCES labs(id),
 parent_id TEXT,
 name TEXT NOT NULL CHECK(length(name) BETWEEN 1 AND 1024),
 generation INTEGER NOT NULL DEFAULT 1 CHECK(typeof(generation)='integer' AND generation>=1),
 UNIQUE(lab_id,id),
 CHECK(parent_id IS NULL OR parent_id!=id),
 FOREIGN KEY(lab_id,parent_id) REFERENCES folders(lab_id,id)
);
CREATE UNIQUE INDEX folder_names ON folders(lab_id,coalesce(parent_id,''),name);
CREATE INDEX folder_parents ON folders(lab_id,parent_id);
CREATE TABLE documents(
 id TEXT PRIMARY KEY NOT NULL,
 lab_id TEXT NOT NULL REFERENCES labs(id),
 folder_id TEXT,
 name TEXT NOT NULL CHECK(length(name) BETWEEN 1 AND 1024),
 generation INTEGER NOT NULL DEFAULT 1 CHECK(typeof(generation)='integer' AND generation>=1),
 head_revision TEXT NOT NULL,
 FOREIGN KEY(lab_id,folder_id) REFERENCES folders(lab_id,id),
 FOREIGN KEY(id,head_revision) REFERENCES revisions(document_id,id) DEFERRABLE INITIALLY DEFERRED
);
CREATE UNIQUE INDEX document_names ON documents(lab_id,coalesce(folder_id,''),name);
CREATE INDEX document_folders ON documents(lab_id,folder_id);
CREATE TABLE revisions(
 id TEXT PRIMARY KEY NOT NULL CHECK(length(id)=64),
 document_id TEXT NOT NULL REFERENCES documents(id),
 parent_id TEXT,
 number INTEGER NOT NULL CHECK(number BETWEEN 1 AND 10000),
 archive BLOB NOT NULL CHECK(length(archive) BETWEEN 1 AND 67108864),
 archive_hash TEXT NOT NULL CHECK(length(archive_hash)=64),
 archive_size INTEGER NOT NULL CHECK(archive_size=length(archive)),
 metadata TEXT NOT NULL CHECK(length(metadata)<=4096),
 created_at INTEGER NOT NULL DEFAULT (unixepoch()) CHECK(created_at>=0),
 UNIQUE(document_id,id),
 UNIQUE(document_id,number),
 FOREIGN KEY(document_id,parent_id) REFERENCES revisions(document_id,id)
);
CREATE TRIGGER immutable_revision_update BEFORE UPDATE ON revisions BEGIN SELECT RAISE(ABORT,'immutable revision'); END;
CREATE TRIGGER immutable_revision_delete BEFORE DELETE ON revisions BEGIN SELECT RAISE(ABORT,'immutable revision'); END;
"#;
