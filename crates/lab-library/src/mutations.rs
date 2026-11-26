//! Idempotent native mutations and private repository checkpoints. The host still
//! authenticates LabAccess and validates native payloads before calling writes.
//! Public intents contain no actor, SQL, remote URL, credential or checkpoint.
use super::*;

pub const MAX_CHECKPOINT_BYTES: usize = 64 * 1024 * 1024;
pub const MAX_CHECKPOINT_TOTAL_BYTES: u64 = 128 * 1024 * 1024;
pub const MAX_CHECKPOINT_WRITES: usize = 8;
pub const MAX_CHECKPOINTS_PER_LAB: u64 = 4096;
pub const MAX_RECEIPT_BYTES: usize = 128 * 1024;
pub const MAX_RECEIPT_TOTAL_BYTES: u64 = 64 * 1024 * 1024;
pub const MAX_RECEIPTS_PER_LAB: u64 = 100_000;
const EXTENSION_VERSION: u32 = 1;

pub type Result<T> = std::result::Result<T, MutationError>;
#[derive(Debug, thiserror::Error)]
pub enum MutationError {
    #[error(transparent)]
    Storage(#[from] super::Error),
    #[error("mutation ID was already used for a different intent")]
    IntentMismatch,
    #[error("private checkpoint changed")]
    CheckpointConflict {
        key: CheckpointKey,
        current: Option<Digest>,
    },
    #[error("invalid reliable mutation field: {0}")]
    Invalid(&'static str),
    #[error("reliable storage limit exceeded: {0}")]
    Limit(&'static str),
    #[error("reliable storage integrity failed: {0}")]
    Corrupt(&'static str),
    #[error("unsupported reliable storage version {0}")]
    Version(u32),
    #[error("native repository transaction aborted: {0}")]
    Aborted(String),
}
impl From<rusqlite::Error> for MutationError {
    fn from(value: rusqlite::Error) -> Self {
        Self::Storage(super::Error::Database(value))
    }
}

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize)]
#[serde(try_from = "String", into = "String")]
pub struct MutationId(String);
impl MutationId {
    pub fn parse(value: impl Into<String>) -> Result<Self> {
        let value = value.into();
        let digits = value
            .strip_prefix("mut-")
            .ok_or(MutationError::Invalid("mutation ID"))?;
        if digits.len() != 32
            || !digits
                .bytes()
                .all(|b| b.is_ascii_digit() || (b'a'..=b'f').contains(&b))
        {
            return Err(MutationError::Invalid("mutation ID"));
        }
        Ok(Self(value))
    }
    pub fn as_str(&self) -> &str {
        &self.0
    }
}
impl TryFrom<String> for MutationId {
    type Error = MutationError;
    fn try_from(value: String) -> Result<Self> {
        Self::parse(value)
    }
}
impl From<MutationId> for String {
    fn from(value: MutationId) -> Self {
        value.0
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Mutation {
    pub id: MutationId,
    pub operation: MutationOperation,
}
#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum MutationOperation {
    CreateFolder {
        parent: Option<FolderId>,
        name: String,
    },
    RenameFolder {
        id: FolderId,
        expected_generation: u64,
        name: String,
    },
    MoveFolder {
        id: FolderId,
        expected_generation: u64,
        parent: Option<FolderId>,
    },
    DeleteEmptyFolder {
        id: FolderId,
        expected_generation: u64,
    },
    ImportDocument {
        folder: Option<FolderId>,
        name: String,
        payload: Payload,
    },
    SaveRevision {
        id: DocumentId,
        base_revision: RevisionId,
        payload: Payload,
    },
    RenameDocument {
        id: DocumentId,
        expected_generation: u64,
        name: String,
    },
    MoveDocument {
        id: DocumentId,
        expected_generation: u64,
        folder: Option<FolderId>,
    },
}
#[derive(Serialize)]
struct Descriptor<'a> {
    kind: PayloadKind,
    archive_hash: &'a Digest,
    archive_size: u64,
    selection: &'a Option<NativeSelection>,
    sharing: Option<WorkspaceSharing>,
}
impl<'a> From<&'a Payload> for Descriptor<'a> {
    fn from(p: &'a Payload) -> Self {
        Self {
            kind: p.kind,
            archive_hash: &p.archive_hash,
            archive_size: p.bytes.len() as u64,
            selection: &p.selection,
            sharing: p.sharing,
        }
    }
}
#[derive(Serialize)]
enum IntentOperation<'a> {
    CreateFolder {
        parent: &'a Option<FolderId>,
        name: &'a str,
    },
    RenameFolder {
        id: &'a FolderId,
        expected_generation: u64,
        name: &'a str,
    },
    MoveFolder {
        id: &'a FolderId,
        expected_generation: u64,
        parent: &'a Option<FolderId>,
    },
    DeleteEmptyFolder {
        id: &'a FolderId,
        expected_generation: u64,
    },
    ImportDocument {
        folder: &'a Option<FolderId>,
        name: &'a str,
        payload: Descriptor<'a>,
    },
    SaveRevision {
        id: &'a DocumentId,
        base_revision: &'a RevisionId,
        payload: Descriptor<'a>,
    },
    RenameDocument {
        id: &'a DocumentId,
        expected_generation: u64,
        name: &'a str,
    },
    MoveDocument {
        id: &'a DocumentId,
        expected_generation: u64,
        folder: &'a Option<FolderId>,
    },
}
impl Mutation {
    /// Canonical wire/native algorithm: domain bytes followed directly by JSON
    /// of the ordered externally tagged operation and bounded payload descriptor.
    pub fn intent_digest(&self) -> Result<Digest> {
        use MutationOperation as O;
        let operation = match &self.operation {
            O::CreateFolder { parent, name } => IntentOperation::CreateFolder { parent, name },
            O::RenameFolder {
                id,
                expected_generation,
                name,
            } => IntentOperation::RenameFolder {
                id,
                expected_generation: *expected_generation,
                name,
            },
            O::MoveFolder {
                id,
                expected_generation,
                parent,
            } => IntentOperation::MoveFolder {
                id,
                expected_generation: *expected_generation,
                parent,
            },
            O::DeleteEmptyFolder {
                id,
                expected_generation,
            } => IntentOperation::DeleteEmptyFolder {
                id,
                expected_generation: *expected_generation,
            },
            O::ImportDocument {
                folder,
                name,
                payload,
            } => IntentOperation::ImportDocument {
                folder,
                name,
                payload: payload.into(),
            },
            O::SaveRevision {
                id,
                base_revision,
                payload,
            } => IntentOperation::SaveRevision {
                id,
                base_revision,
                payload: payload.into(),
            },
            O::RenameDocument {
                id,
                expected_generation,
                name,
            } => IntentOperation::RenameDocument {
                id,
                expected_generation: *expected_generation,
                name,
            },
            O::MoveDocument {
                id,
                expected_generation,
                folder,
            } => IntentOperation::MoveDocument {
                id,
                expected_generation: *expected_generation,
                folder,
            },
        };
        preflight(&self.operation)?;
        let json = serde_json::to_vec(&operation).map_err(|_| MutationError::Invalid("intent"))?;
        let mut hash = blake3::Hasher::new();
        hash.update(b"collagene-library-mutation-v1\0");
        hash.update(&json);
        Ok(Digest::parse(hash.finalize().to_hex().to_string())?)
    }
}
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub enum Applied {
    CreatedFolder(Folder),
    RenamedFolder(Folder),
    MovedFolder(Folder),
    DeletedEmptyFolder(FolderId),
    Imported(Document),
    Saved {
        document: Document,
        revision: RevisionSummary,
    },
    RenamedDocument(Document),
    MovedDocument(Document),
}
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum RejectionCode {
    Invalid,
    NotFound,
    AlreadyExists,
    NotEmpty,
    TooLarge,
    Conflict,
}
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub enum MutationConflict {
    Folder { current: Folder },
    Document { current: Document },
}
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct Rejection {
    pub code: RejectionCode,
    pub message: String,
    pub conflict: Option<MutationConflict>,
}
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub enum MutationOutcome {
    Applied(Applied),
    Rejected(Rejection),
}
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct Receipt {
    pub id: MutationId,
    pub intent: Digest,
    pub catalog_generation: u64,
    pub outcome: MutationOutcome,
}
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ApplyResult {
    pub receipt: Receipt,
    pub replayed: bool,
}

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
pub struct CheckpointKey {
    pub namespace: String,
    pub key: String,
}
impl CheckpointKey {
    pub fn repository(key: impl Into<String>) -> Self {
        Self {
            namespace: "repository".into(),
            key: key.into(),
        }
    }
    fn validate(&self) -> Result<()> {
        if self.namespace != "repository"
            || self.key.is_empty()
            || self.key.len() > 128
            || !self
                .key
                .bytes()
                .all(|b| b.is_ascii_alphanumeric() || matches!(b, b'-' | b'_'))
        {
            return Err(MutationError::Invalid("checkpoint key"));
        }
        Ok(())
    }
}
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CheckpointValue {
    pub version: u32,
    pub bytes: Vec<u8>,
}
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Checkpoint {
    pub key: CheckpointKey,
    pub value: CheckpointValue,
    pub digest: Digest,
}
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CheckpointWrite {
    pub key: CheckpointKey,
    pub expected: Option<Digest>,
    pub value: Option<CheckpointValue>,
}

impl Library {
    /// Native ID generation only; not an authorization credential.
    pub fn new_mutation_id(&self) -> Result<MutationId> {
        let id: String =
            self.connection
                .query_row("SELECT 'mut-'||lower(hex(randomblob(16)))", [], |r| {
                    r.get(0)
                })?;
        MutationId::parse(id)
    }
    /// Owned-file local repository scope, never a network authentication API.
    pub fn local_repository_access(&mut self, namespace: &str, name: &str) -> Result<LabAccess> {
        if namespace.is_empty()
            || namespace.len() > 64
            || !namespace
                .bytes()
                .all(|b| b.is_ascii_alphanumeric() || matches!(b, b'-' | b'_'))
        {
            return Err(MutationError::Invalid("local namespace"));
        }
        valid_name(name)?;
        if self.connection.is_readonly(rusqlite::DatabaseName::Main)? {
            let version: u32 = self.connection.query_row(
                "SELECT version FROM cg_reliable_schema WHERE singleton=1",
                [],
                |r| r.get(0),
            )?;
            if version != EXTENSION_VERSION {
                return Err(MutationError::Version(version));
            }
            let id: String = self.connection.query_row(
                "SELECT lab_id FROM cg_reliable_local_scopes WHERE namespace=?1",
                [namespace],
                |r| r.get(0),
            )?;
            let lab = LabId::parse(id)?;
            ensure_lab(&self.connection, &lab)?;
            return Ok(LabAccess { lab });
        }
        let tx = self
            .connection
            .transaction_with_behavior(TransactionBehavior::Immediate)?;
        extension(&tx)?;
        let existing: Option<String> = tx
            .query_row(
                "SELECT lab_id FROM cg_reliable_local_scopes WHERE namespace=?1",
                [namespace],
                |r| r.get(0),
            )
            .optional()?;
        let lab = if let Some(id) = existing {
            let lab = LabId::parse(id)?;
            ensure_lab(&tx, &lab)?;
            lab
        } else {
            let count: u64 =
                tx.query_row("SELECT count(*) FROM cg_reliable_local_scopes", [], |r| {
                    r.get(0)
                })?;
            if count >= 128 {
                return Err(MutationError::Limit("local repository scopes"));
            }
            let lab = LabId::parse(random_id(&tx, "lab-")?)?;
            tx.execute(
                "INSERT INTO labs(id,name) VALUES(?1,?2)",
                params![lab.as_str(), name],
            )?;
            tx.execute(
                "INSERT INTO cg_reliable_local_scopes(namespace,lab_id) VALUES(?1,?2)",
                params![namespace, lab.as_str()],
            )?;
            lab
        };
        tx.commit()?;
        Ok(LabAccess { lab })
    }
    pub fn apply_reliable(
        &mut self,
        access: &LabAccess,
        mutation: &Mutation,
    ) -> Result<ApplyResult> {
        self.apply_reliable_with(access, mutation, |_| Ok(Vec::new()))
    }
    /// Finalize is trusted native code, not serialized client behavior. It runs
    /// once for a new receipt, including a domain rejection, inside the same tx.
    /// A replay never reruns it because its prior checkpoint effect already committed.
    pub fn apply_reliable_with<F>(
        &mut self,
        access: &LabAccess,
        mutation: &Mutation,
        finalize: F,
    ) -> Result<ApplyResult>
    where
        F: FnOnce(&Receipt) -> Result<Vec<CheckpointWrite>>,
    {
        let intent = mutation.intent_digest()?;
        let tx = self
            .connection
            .transaction_with_behavior(TransactionBehavior::Immediate)?;
        ensure_lab(&tx, &access.lab)?;
        extension(&tx)?;
        usage(&tx, &access.lab)?;
        if let Some(receipt) = lookup(&tx, &access.lab, &mutation.id, &intent)? {
            tx.commit()?;
            return Ok(ApplyResult {
                receipt,
                replayed: true,
            });
        }
        let (count, bytes): (u64, u64) = tx.query_row(
            "SELECT receipt_count,receipt_bytes FROM cg_reliable_usage WHERE lab_id=?1",
            [access.lab.as_str()],
            |r| Ok((r.get(0)?, r.get(1)?)),
        )?;
        if count >= MAX_RECEIPTS_PER_LAB || bytes >= MAX_RECEIPT_TOTAL_BYTES {
            return Err(MutationError::Limit("mutation receipts"));
        }
        tx.execute_batch("SAVEPOINT reliable_effect")?;
        let outcome = match execute(
            &tx,
            &access.lab,
            &mutation.operation,
            self.options.max_lab_archive_bytes,
        ) {
            Ok(applied) => {
                tx.execute_batch("RELEASE reliable_effect")?;
                MutationOutcome::Applied(applied)
            }
            Err(error) => {
                tx.execute_batch("ROLLBACK TO reliable_effect; RELEASE reliable_effect")?;
                MutationOutcome::Rejected(rejection(&tx, &access.lab, &mutation.operation, error)?)
            }
        };
        let receipt = Receipt {
            id: mutation.id.clone(),
            intent,
            catalog_generation: catalog_generation(&tx, &access.lab)?,
            outcome,
        };
        let encoded =
            serde_json::to_vec(&receipt).map_err(|_| MutationError::Corrupt("receipt encoding"))?;
        if encoded.len() > MAX_RECEIPT_BYTES
            || bytes + encoded.len() as u64 > MAX_RECEIPT_TOTAL_BYTES
        {
            return Err(MutationError::Limit("mutation receipts"));
        }
        let writes = finalize(&receipt)?;
        checkpoints(&tx, &access.lab, &writes)?;
        tx.execute("INSERT INTO cg_reliable_mutations(lab_id,id,intent,receipt,receipt_hash) VALUES(?1,?2,?3,?4,?5)",params![access.lab.as_str(),mutation.id.as_str(),receipt.intent.as_str(),&encoded,blake3::hash(&encoded).to_hex().as_str()])?;
        tx.execute("UPDATE cg_reliable_usage SET receipt_count=receipt_count+1,receipt_bytes=receipt_bytes+?2 WHERE lab_id=?1",params![access.lab.as_str(),encoded.len() as u64])?;
        tx.commit()?;
        Ok(ApplyResult {
            receipt,
            replayed: false,
        })
    }
    pub fn lookup_reliable(
        &mut self,
        access: &LabAccess,
        id: &MutationId,
        expected_intent: &Digest,
    ) -> Result<Option<Receipt>> {
        let tx = self
            .connection
            .transaction_with_behavior(TransactionBehavior::Immediate)?;
        ensure_lab(&tx, &access.lab)?;
        extension(&tx)?;
        let receipt = lookup(&tx, &access.lab, id, expected_intent)?;
        tx.commit()?;
        Ok(receipt)
    }
    pub fn read_checkpoint(
        &mut self,
        access: &LabAccess,
        key: &CheckpointKey,
    ) -> Result<Option<Checkpoint>> {
        key.validate()?;
        let behavior = if self.connection.is_readonly(rusqlite::DatabaseName::Main)? {
            TransactionBehavior::Deferred
        } else {
            TransactionBehavior::Immediate
        };
        let tx = self.connection.transaction_with_behavior(behavior)?;
        ensure_lab(&tx, &access.lab)?;
        extension(&tx)?;
        let value = read_checkpoint(&tx, &access.lab, key)?;
        tx.commit()?;
        Ok(value)
    }
    /// Native metadata retirement transaction: at most 128 checkpoint deletions
    /// and three replacement records. Ordinary callers retain the eight-write
    /// limit. Every key, byte quota and expected digest uses the same validator.
    pub fn commit_checkpoint_retirements(
        &mut self,
        access: &LabAccess,
        writes: &[CheckpointWrite],
    ) -> Result<()> {
        let replacements = writes.iter().filter(|write| write.value.is_some()).count();
        if replacements > 3 || writes.len().saturating_sub(replacements) > 128 {
            return Err(MutationError::Limit("checkpoint retirement batch"));
        }
        let tx = self
            .connection
            .transaction_with_behavior(TransactionBehavior::Immediate)?;
        ensure_lab(&tx, &access.lab)?;
        extension(&tx)?;
        usage(&tx, &access.lab)?;
        checkpoints_bounded(&tx, &access.lab, writes, 131)?;
        tx.commit()?;
        Ok(())
    }
    /// Atomic native-only acknowledgments/remote bookkeeping, without a DNA mutation.
    pub fn commit_checkpoints(
        &mut self,
        access: &LabAccess,
        writes: &[CheckpointWrite],
    ) -> Result<()> {
        let tx = self
            .connection
            .transaction_with_behavior(TransactionBehavior::Immediate)?;
        ensure_lab(&tx, &access.lab)?;
        extension(&tx)?;
        usage(&tx, &access.lab)?;
        checkpoints(&tx, &access.lab, writes)?;
        tx.commit()?;
        Ok(())
    }
}

fn preflight(operation: &MutationOperation) -> Result<()> {
    use MutationOperation as O;
    let (name, payload) = match operation {
        O::CreateFolder { name, .. }
        | O::RenameFolder { name, .. }
        | O::RenameDocument { name, .. } => (Some(name), None),
        O::ImportDocument { name, payload, .. } => (Some(name), Some(payload)),
        O::SaveRevision { payload, .. } => (None, Some(payload)),
        _ => (None, None),
    };
    if name.is_some_and(|name| name.len() > 1024) {
        return Err(MutationError::Limit("mutation name"));
    }
    if let Some(payload) = payload {
        payload.validate()?;
    }
    Ok(())
}
fn execute(
    tx: &Transaction<'_>,
    lab: &LabId,
    operation: &MutationOperation,
    quota: u64,
) -> super::Result<Applied> {
    use MutationOperation as O;
    Ok(match operation {
        O::CreateFolder { parent, name } => {
            valid_name(name)?;
            if parent_depth(tx, lab, parent.as_ref(), None)? + 1 > MAX_FOLDER_DEPTH {
                return Err(Error::Limit("folder depth"));
            }
            entity_capacity(tx, lab)?;
            unique_name(tx, lab, parent.as_ref(), name, None, None)?;
            let id = FolderId::parse(random_id(tx, "fld-")?)?;
            tx.execute(
                "INSERT INTO folders(id,lab_id,parent_id,name) VALUES(?1,?2,?3,?4)",
                params![
                    id.as_str(),
                    lab.as_str(),
                    parent.as_ref().map(FolderId::as_str),
                    name
                ],
            )?;
            bump_catalog(tx, lab)?;
            Applied::CreatedFolder(get_folder(tx, lab, &id)?)
        }
        O::RenameFolder {
            id,
            expected_generation,
            name,
        } => {
            valid_name(name)?;
            let current = get_folder(tx, lab, id)?;
            metadata_cas(*expected_generation, current.generation)?;
            unique_name(tx, lab, current.parent.as_ref(), name, Some(id), None)?;
            tx.execute(
                "UPDATE folders SET name=?1,generation=generation+1 WHERE id=?2 AND lab_id=?3",
                params![name, id.as_str(), lab.as_str()],
            )?;
            bump_catalog(tx, lab)?;
            Applied::RenamedFolder(get_folder(tx, lab, id)?)
        }
        O::MoveFolder {
            id,
            expected_generation,
            parent,
        } => {
            let current = get_folder(tx, lab, id)?;
            metadata_cas(*expected_generation, current.generation)?;
            let depth = parent_depth(tx, lab, parent.as_ref(), Some(id))?;
            let height:usize=tx.query_row("WITH RECURSIVE branch(id,depth) AS (SELECT id,1 FROM folders WHERE id=?1 AND lab_id=?2 UNION ALL SELECT f.id,b.depth+1 FROM folders f JOIN branch b ON f.parent_id=b.id WHERE f.lab_id=?2 AND b.depth<=?3) SELECT max(depth) FROM branch",params![id.as_str(),lab.as_str(),MAX_FOLDER_DEPTH as u64],|r|r.get(0))?;
            if depth + height > MAX_FOLDER_DEPTH {
                return Err(Error::Limit("folder depth"));
            }
            unique_name(tx, lab, parent.as_ref(), &current.name, Some(id), None)?;
            tx.execute(
                "UPDATE folders SET parent_id=?1,generation=generation+1 WHERE id=?2 AND lab_id=?3",
                params![
                    parent.as_ref().map(FolderId::as_str),
                    id.as_str(),
                    lab.as_str()
                ],
            )?;
            bump_catalog(tx, lab)?;
            Applied::MovedFolder(get_folder(tx, lab, id)?)
        }
        O::DeleteEmptyFolder {
            id,
            expected_generation,
        } => {
            let current = get_folder(tx, lab, id)?;
            metadata_cas(*expected_generation, current.generation)?;
            let children:bool=tx.query_row("SELECT EXISTS(SELECT 1 FROM folders WHERE lab_id=?1 AND parent_id=?2 UNION ALL SELECT 1 FROM documents WHERE lab_id=?1 AND folder_id=?2)",params![lab.as_str(),id.as_str()],|r|r.get(0))?;
            if children {
                return Err(Error::NotEmpty);
            }
            tx.execute(
                "DELETE FROM folders WHERE lab_id=?1 AND id=?2",
                params![lab.as_str(), id.as_str()],
            )?;
            bump_catalog(tx, lab)?;
            Applied::DeletedEmptyFolder(id.clone())
        }
        O::ImportDocument {
            folder,
            name,
            payload,
        } => {
            valid_name(name)?;
            parent_depth(tx, lab, folder.as_ref(), None)?;
            entity_capacity(tx, lab)?;
            unique_name(tx, lab, folder.as_ref(), name, None, None)?;
            let id = DocumentId::parse(random_id(tx, "doc-")?)?;
            charge_archive(tx, lab, payload.bytes.len(), quota)?;
            let revision = revision_id(&id, None, payload)?;
            tx.execute("INSERT INTO documents(id,lab_id,folder_id,name,head_revision) VALUES(?1,?2,?3,?4,?5)",params![id.as_str(),lab.as_str(),folder.as_ref().map(FolderId::as_str),name,revision.as_str()])?;
            insert_revision(tx, &id, &revision, None, 1, payload)?;
            bump_catalog(tx, lab)?;
            Applied::Imported(get_document(tx, lab, &id)?)
        }
        O::SaveRevision {
            id,
            base_revision,
            payload,
        } => {
            let current = get_document(tx, lab, id)?;
            if &current.head.id != base_revision {
                return Err(Error::Conflict {
                    expected: base_revision.clone(),
                    actual: current.head.id,
                });
            }
            if current.head.number >= MAX_DOCUMENT_REVISIONS {
                return Err(Error::Limit("document revisions"));
            }
            charge_archive(tx, lab, payload.bytes.len(), quota)?;
            let revision = revision_id(id, Some(base_revision), payload)?;
            insert_revision(
                tx,
                id,
                &revision,
                Some(base_revision),
                current.head.number + 1,
                payload,
            )?;
            let changed=tx.execute("UPDATE documents SET head_revision=?1,generation=generation+1 WHERE id=?2 AND lab_id=?3 AND head_revision=?4",params![revision.as_str(),id.as_str(),lab.as_str(),base_revision.as_str()])?;
            if changed != 1 {
                return Err(Error::Corrupt("CAS update"));
            }
            bump_catalog(tx, lab)?;
            let document = get_document(tx, lab, id)?;
            Applied::Saved {
                revision: document.head.clone(),
                document,
            }
        }
        O::RenameDocument {
            id,
            expected_generation,
            name,
        } => {
            valid_name(name)?;
            let current = get_document(tx, lab, id)?;
            metadata_cas(*expected_generation, current.generation)?;
            unique_name(tx, lab, current.folder.as_ref(), name, None, Some(id))?;
            tx.execute(
                "UPDATE documents SET name=?1,generation=generation+1 WHERE id=?2 AND lab_id=?3",
                params![name, id.as_str(), lab.as_str()],
            )?;
            bump_catalog(tx, lab)?;
            Applied::RenamedDocument(get_document(tx, lab, id)?)
        }
        O::MoveDocument {
            id,
            expected_generation,
            folder,
        } => {
            let current = get_document(tx, lab, id)?;
            metadata_cas(*expected_generation, current.generation)?;
            parent_depth(tx, lab, folder.as_ref(), None)?;
            unique_name(tx, lab, folder.as_ref(), &current.name, None, Some(id))?;
            tx.execute("UPDATE documents SET folder_id=?1,generation=generation+1 WHERE id=?2 AND lab_id=?3",params![folder.as_ref().map(FolderId::as_str),id.as_str(),lab.as_str()])?;
            bump_catalog(tx, lab)?;
            Applied::MovedDocument(get_document(tx, lab, id)?)
        }
    })
}
fn rejection(
    connection: &Connection,
    lab: &LabId,
    operation: &MutationOperation,
    error: Error,
) -> Result<Rejection> {
    use MutationOperation as O;
    let (code, message) = match error {
        Error::Invalid(_) => (
            RejectionCode::Invalid,
            "The mutation fields or folder geometry are invalid.",
        ),
        Error::NotFound => (
            RejectionCode::NotFound,
            "The item is not present in this lab.",
        ),
        Error::AlreadyExists => (
            RejectionCode::AlreadyExists,
            "An item with that name already exists in this folder.",
        ),
        Error::NotEmpty => (
            RejectionCode::NotEmpty,
            "The folder contains retained items and was not deleted.",
        ),
        Error::Limit(_) => (
            RejectionCode::TooLarge,
            "The operation exceeds the library's storage limits.",
        ),
        Error::Conflict { .. } | Error::MetadataConflict { .. } => (
            RejectionCode::Conflict,
            "The item changed; retain your edits and review the recorded current state.",
        ),
        other => return Err(other.into()),
    };
    let conflict = if code == RejectionCode::Conflict {
        match operation {
            O::RenameFolder { id, .. }
            | O::MoveFolder { id, .. }
            | O::DeleteEmptyFolder { id, .. } => Some(MutationConflict::Folder {
                current: get_folder(connection, lab, id)?,
            }),
            O::SaveRevision { id, .. }
            | O::RenameDocument { id, .. }
            | O::MoveDocument { id, .. } => Some(MutationConflict::Document {
                current: get_document(connection, lab, id)?,
            }),
            _ => return Err(MutationError::Corrupt("conflict operation")),
        }
    } else {
        None
    };
    Ok(Rejection {
        code,
        message: message.into(),
        conflict,
    })
}
fn lookup(
    connection: &Connection,
    lab: &LabId,
    id: &MutationId,
    expected: &Digest,
) -> Result<Option<Receipt>> {
    let row:Option<(String,Option<Vec<u8>>,String)>=connection.query_row("SELECT intent,CASE WHEN length(receipt)<=?3 THEN receipt ELSE NULL END,receipt_hash FROM cg_reliable_mutations WHERE lab_id=?1 AND id=?2",params![lab.as_str(),id.as_str(),MAX_RECEIPT_BYTES as u64],|r|Ok((r.get(0)?,r.get(1)?,r.get(2)?))).optional()?;
    let Some((intent, bytes, hash)) = row else {
        return Ok(None);
    };
    if intent != expected.as_str() {
        return Err(MutationError::IntentMismatch);
    }
    let bytes = bytes.ok_or(MutationError::Corrupt("receipt bounds"))?;
    if blake3::hash(&bytes).to_hex().as_str() != hash {
        return Err(MutationError::Corrupt("receipt checksum"));
    }
    let receipt: Receipt =
        serde_json::from_slice(&bytes).map_err(|_| MutationError::Corrupt("receipt fields"))?;
    if &receipt.id != id
        || &receipt.intent != expected
        || receipt.catalog_generation > catalog_generation(connection, lab)?
    {
        return Err(MutationError::Corrupt("receipt binding"));
    }
    Ok(Some(receipt))
}
fn checkpoint_digest(lab: &LabId, key: &CheckpointKey, value: &CheckpointValue) -> Result<Digest> {
    let mut hash = blake3::Hasher::new();
    hash.update(b"collagene-private-checkpoint-v1\0");
    for part in [
        lab.as_str().as_bytes(),
        key.namespace.as_bytes(),
        key.key.as_bytes(),
        &value.version.to_le_bytes(),
        &value.bytes,
    ] {
        hash.update(&(part.len() as u64).to_le_bytes());
        hash.update(part);
    }
    Ok(Digest::parse(hash.finalize().to_hex().to_string())?)
}
fn read_checkpoint(
    connection: &Connection,
    lab: &LabId,
    key: &CheckpointKey,
) -> Result<Option<Checkpoint>> {
    let row:Option<(u32,Option<Vec<u8>>,String)>=connection.query_row("SELECT version,CASE WHEN length(value)<=?4 THEN value ELSE NULL END,digest FROM cg_reliable_checkpoints WHERE lab_id=?1 AND namespace=?2 AND key=?3",params![lab.as_str(),key.namespace,key.key,MAX_CHECKPOINT_BYTES as u64],|r|Ok((r.get(0)?,r.get(1)?,r.get(2)?))).optional()?;
    let Some((version, bytes, digest)) = row else {
        return Ok(None);
    };
    let value = CheckpointValue {
        version,
        bytes: bytes.ok_or(MutationError::Corrupt("checkpoint size"))?,
    };
    if version == 0 {
        return Err(MutationError::Corrupt("checkpoint version"));
    }
    let digest = Digest::parse(digest)?;
    if checkpoint_digest(lab, key, &value)? != digest {
        return Err(MutationError::Corrupt("checkpoint binding"));
    }
    Ok(Some(Checkpoint {
        key: key.clone(),
        value,
        digest,
    }))
}
fn checkpoints(tx: &Transaction<'_>, lab: &LabId, writes: &[CheckpointWrite]) -> Result<()> {
    checkpoints_bounded(tx, lab, writes, MAX_CHECKPOINT_WRITES)
}
fn checkpoints_bounded(
    tx: &Transaction<'_>,
    lab: &LabId,
    writes: &[CheckpointWrite],
    limit: usize,
) -> Result<()> {
    if writes.len() > limit {
        return Err(MutationError::Limit("checkpoint batch"));
    }
    let mut seen = std::collections::BTreeSet::new();
    let mut charge = 0u64;
    let mut old_bytes = 0u64;
    let mut created = 0u64;
    let mut deleted = 0u64;
    for write in writes {
        write.key.validate()?;
        if !seen.insert(&write.key) {
            return Err(MutationError::Invalid("duplicate checkpoint key"));
        }
        if let Some(value) = &write.value {
            if value.version == 0 {
                return Err(MutationError::Invalid("checkpoint version"));
            }
            if value.bytes.len() > MAX_CHECKPOINT_BYTES {
                return Err(MutationError::Limit("checkpoint bytes"));
            }
            charge += value.bytes.len() as u64;
        }
        let current = read_checkpoint(tx, lab, &write.key)?;
        if current.as_ref().map(|value| &value.digest) != write.expected.as_ref() {
            return Err(MutationError::CheckpointConflict {
                key: write.key.clone(),
                current: current.map(|value| value.digest),
            });
        }
        if current.is_none() && write.value.is_some() {
            created += 1;
        }
        if current.is_some() && write.value.is_none() {
            deleted += 1;
        }
        old_bytes += current
            .map(|value| value.value.bytes.len() as u64)
            .unwrap_or(0);
    }
    let count: u64 = tx.query_row(
        "SELECT count(*) FROM cg_reliable_checkpoints WHERE lab_id=?1",
        [lab.as_str()],
        |r| r.get(0),
    )?;
    if count + created - deleted > MAX_CHECKPOINTS_PER_LAB {
        return Err(MutationError::Limit("checkpoint count"));
    }
    let current: u64 = tx.query_row(
        "SELECT checkpoint_bytes FROM cg_reliable_usage WHERE lab_id=?1",
        [lab.as_str()],
        |r| r.get(0),
    )?;
    let next = current
        .checked_sub(old_bytes)
        .and_then(|n| n.checked_add(charge))
        .ok_or(MutationError::Corrupt("checkpoint accounting"))?;
    if charge > MAX_CHECKPOINT_TOTAL_BYTES || next > MAX_CHECKPOINT_TOTAL_BYTES {
        return Err(MutationError::Limit("checkpoint total bytes"));
    }
    for write in writes {
        if let Some(value) = &write.value {
            let digest = checkpoint_digest(lab, &write.key, value)?;
            tx.execute("INSERT INTO cg_reliable_checkpoints(lab_id,namespace,key,version,value,digest) VALUES(?1,?2,?3,?4,?5,?6) ON CONFLICT(lab_id,namespace,key) DO UPDATE SET version=excluded.version,value=excluded.value,digest=excluded.digest",params![lab.as_str(),write.key.namespace,write.key.key,value.version,&value.bytes,digest.as_str()])?;
        } else {
            tx.execute(
                "DELETE FROM cg_reliable_checkpoints WHERE lab_id=?1 AND namespace=?2 AND key=?3",
                params![lab.as_str(), write.key.namespace, write.key.key],
            )?;
        }
    }
    tx.execute(
        "UPDATE cg_reliable_usage SET checkpoint_bytes=?2 WHERE lab_id=?1",
        params![lab.as_str(), next],
    )?;
    Ok(())
}
fn usage(tx: &Transaction<'_>, lab: &LabId) -> Result<()> {
    tx.execute(
        "INSERT OR IGNORE INTO cg_reliable_usage(lab_id) VALUES(?1)",
        [lab.as_str()],
    )?;
    Ok(())
}
fn extension(tx: &Transaction<'_>) -> Result<()> {
    let exists:bool=tx.query_row("SELECT EXISTS(SELECT 1 FROM sqlite_master WHERE type='table' AND name='cg_reliable_schema')",[],|r|r.get(0))?;
    if exists {
        let version: u32 = tx.query_row(
            "SELECT version FROM cg_reliable_schema WHERE singleton=1",
            [],
            |r| r.get(0),
        )?;
        if version != EXTENSION_VERSION {
            return Err(MutationError::Version(version));
        }
    } else {
        let count: u64 = tx.query_row(
            "SELECT count(*) FROM sqlite_master WHERE name GLOB 'cg_reliable_*'",
            [],
            |r| r.get(0),
        )?;
        if count != 0 {
            return Err(MutationError::Corrupt("partial extension schema"));
        }
        tx.execute_batch(EXTENSION_SCHEMA)?;
    }
    Ok(())
}
const EXTENSION_SCHEMA: &str = r#"
CREATE TABLE cg_reliable_schema(singleton INTEGER PRIMARY KEY CHECK(singleton=1),version INTEGER NOT NULL);
INSERT INTO cg_reliable_schema VALUES(1,1);
CREATE TABLE cg_reliable_usage(
 lab_id TEXT PRIMARY KEY REFERENCES labs(id),
 receipt_count INTEGER NOT NULL DEFAULT 0 CHECK(receipt_count BETWEEN 0 AND 100000),
 receipt_bytes INTEGER NOT NULL DEFAULT 0 CHECK(receipt_bytes BETWEEN 0 AND 67108864),
 checkpoint_bytes INTEGER NOT NULL DEFAULT 0 CHECK(checkpoint_bytes BETWEEN 0 AND 134217728)
);
CREATE TABLE cg_reliable_mutations(
 lab_id TEXT NOT NULL REFERENCES labs(id),id TEXT NOT NULL CHECK(length(id)=36),
 intent TEXT NOT NULL CHECK(length(intent)=64),receipt BLOB NOT NULL CHECK(length(receipt)<=131072),
 receipt_hash TEXT NOT NULL CHECK(length(receipt_hash)=64),PRIMARY KEY(lab_id,id)
);
CREATE TRIGGER cg_reliable_receipt_update BEFORE UPDATE ON cg_reliable_mutations BEGIN SELECT RAISE(ABORT,'immutable receipt');END;
CREATE TRIGGER cg_reliable_receipt_delete BEFORE DELETE ON cg_reliable_mutations BEGIN SELECT RAISE(ABORT,'immutable receipt');END;
CREATE TABLE cg_reliable_checkpoints(
 lab_id TEXT NOT NULL REFERENCES labs(id),namespace TEXT NOT NULL CHECK(namespace='repository'),
 key TEXT NOT NULL CHECK(length(key) BETWEEN 1 AND 128),version INTEGER NOT NULL CHECK(version>0),
 value BLOB NOT NULL CHECK(length(value)<=67108864),digest TEXT NOT NULL CHECK(length(digest)=64),
 PRIMARY KEY(lab_id,namespace,key)
);
CREATE TABLE cg_reliable_local_scopes(namespace TEXT PRIMARY KEY,lab_id TEXT NOT NULL UNIQUE REFERENCES labs(id));
"#;
