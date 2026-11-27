//! Managed saved-catalog backup. No working editor file, external credential,
//! frontend path, remote transport or biological archive IPC is involved.
use crate::{HostControl, HostError, HostResult, JobStart, WindowAccess, permissions};
use collagene_desktop_protocol::OperationId;
use collagene_lab_library::{self as store, backup};
use collagene_library_protocol as base;
use rusqlite::OptionalExtension;
use serde::{Deserialize, Serialize};
use std::{
    fs::{self, File, OpenOptions},
    io::{Read, Write},
    path::{Path, PathBuf},
};
#[path = "catalog_activation.rs"]
mod activation;
pub use activation::catalog_activation_pending;
pub(crate) use activation::{
    CatalogLease, catalog_generation, catalog_session_status, check_catalog_connection,
    current_catalog_path, startup_catalog_path,
};
const MAX_RECORD: usize = 64 * 1024;
const MAX_OPERATIONS: usize = 4096;
const MAX_PAGE: u16 = 200;
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Receipt {
    pub id: String,
    pub basename: String,
    pub archive_hash: base::Digest,
    pub archive_bytes: u64,
    pub created_unix_seconds: u64,
    pub database_hash: base::Digest,
    pub database_bytes: u64,
    pub schema_hash: base::Digest,
    pub schema_version: u32,
    pub lab: base::LabId,
    pub lab_count: u32,
    pub documents: u64,
    pub folders: u64,
    pub revisions: u64,
    pub mutation_receipts: u64,
    /// Observation before the initial snapshot request, not a claim that no writer
    /// committed while SQLite obtained its consistent snapshot.
    pub request_checkpoint: Option<base::Digest>,
}
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Entry {
    pub id: String,
    pub basename: String,
    pub receipt: Option<Receipt>,
}
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Page {
    pub entries: Vec<Entry>,
    pub next_offset: Option<u32>,
    pub catalog_revision: base::Digest,
}
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct RecoveryReceipt {
    pub operation_id: String,
    pub backup: Receipt,
    pub directory_basename: String,
    pub database_basename: String,
}
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
enum Action {
    Create {
        request_checkpoint: Option<base::Digest>,
    },
    Restore {
        backup_id: String,
        archive_hash: base::Digest,
    },
}
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
struct Intent {
    version: u32,
    operation_id: String,
    lab: base::LabId,
    action: Action,
}
#[derive(Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
struct Record<T> {
    version: u32,
    hash: base::Digest,
    value: T,
}
pub struct CatalogBackupService {
    root: PathBuf,
    backups: PathBuf,
    operations: PathBuf,
    recoveries: PathBuf,
    lab: store::LabId,
    library: store::Library,
    database: PathBuf,
}
fn invalid(message: &'static str) -> HostError {
    HostError::Invalid(message)
}
fn io<T>(_: T) -> HostError {
    HostError::Unavailable
}
fn convert_digest(d: &store::Digest) -> HostResult<base::Digest> {
    base::Digest::new(d.as_str()).map_err(io)
}
fn parse_operation(value: &str) -> HostResult<OperationId> {
    OperationId::new(value).map_err(|_| invalid("Unknown or invalid managed backup ID."))
}
fn private(path: &Path, directory: bool) -> HostResult<()> {
    let m = fs::symlink_metadata(path).map_err(io)?;
    if m.file_type().is_symlink() || if directory { !m.is_dir() } else { !m.is_file() } {
        return Err(invalid(
            "Managed backup paths must remain private regular files/directories.",
        ));
    }
    #[cfg(unix)]
    {
        use std::os::unix::fs::{MetadataExt, PermissionsExt};
        if m.uid() != unsafe { libc::geteuid() }
            || m.permissions().mode() & 0o777 != if directory { 0o700 } else { 0o600 }
        {
            return Err(invalid(
                "Managed backups require current-user-only storage permissions.",
            ));
        }
    }
    Ok(())
}
fn mkdir(path: &Path) -> HostResult<()> {
    let mut b = fs::DirBuilder::new();
    #[cfg(unix)]
    {
        use std::os::unix::fs::DirBuilderExt;
        b.mode(0o700);
    }
    match b.create(path) {
        Ok(()) => sync(path.parent().ok_or(HostError::Unavailable)?)?,
        Err(e) if e.kind() == std::io::ErrorKind::AlreadyExists => (),
        Err(e) => return Err(io(e)),
    }
    private(path, true)
}
fn sync(path: &Path) -> HostResult<()> {
    #[cfg(unix)]
    File::open(path).map_err(io)?.sync_all().map_err(io)?;
    let _ = path;
    Ok(())
}
fn exists(path: &Path) -> HostResult<bool> {
    match fs::symlink_metadata(path) {
        Ok(_) => Ok(true),
        Err(e) if e.kind() == std::io::ErrorKind::NotFound => Ok(false),
        Err(e) => Err(io(e)),
    }
}
fn read<T: for<'a> Deserialize<'a> + Serialize>(path: &Path) -> HostResult<T> {
    private(path, false)?;
    let file = File::open(path).map_err(io)?;
    if file.metadata().map_err(io)?.len() > MAX_RECORD as u64 {
        return Err(HostError::Corrupt);
    }
    let mut bytes = vec![];
    file.take(MAX_RECORD as u64 + 1)
        .read_to_end(&mut bytes)
        .map_err(io)?;
    if bytes.len() > MAX_RECORD {
        return Err(HostError::Corrupt);
    }
    let record: Record<T> = serde_json::from_slice(&bytes).map_err(|_| HostError::Corrupt)?;
    let value = serde_json::to_vec(&record.value).map_err(|_| HostError::Corrupt)?;
    if record.version != 1 || base::Digest::of(&value) != record.hash {
        return Err(HostError::Corrupt);
    }
    Ok(record.value)
}
fn write<T: Serialize>(path: &Path, value: &T) -> HostResult<()> {
    use std::sync::atomic::{AtomicU64, Ordering};
    static NEXT: AtomicU64 = AtomicU64::new(0);
    let value_bytes = serde_json::to_vec(value).map_err(io)?;
    let bytes = serde_json::to_vec(&Record {
        version: 1,
        hash: base::Digest::of(&value_bytes),
        value,
    })
    .map_err(io)?;
    if bytes.len() > MAX_RECORD {
        return Err(invalid("Managed backup record exceeds its bound."));
    }
    let parent = path.parent().ok_or(HostError::Unavailable)?;
    let staging = parent.join("staging");
    private(&staging, true)?;
    if fs::read_dir(&staging).map_err(io)?.take(65).count() > 64 {
        return Err(invalid(
            "Review abandoned backup record staging before continuing.",
        ));
    }
    struct Pending(PathBuf);
    impl Drop for Pending {
        fn drop(&mut self) {
            let _ = fs::remove_file(&self.0);
        }
    }
    let mut allocated = None;
    for _ in 0..100 {
        let p = staging.join(format!(
            "record-{}-{}",
            std::process::id(),
            NEXT.fetch_add(1, Ordering::Relaxed)
        ));
        let mut options = OpenOptions::new();
        options.write(true).create_new(true);
        #[cfg(unix)]
        {
            use std::os::unix::fs::OpenOptionsExt;
            options.mode(0o600).custom_flags(libc::O_NOFOLLOW);
        }
        match options.open(&p) {
            Ok(file) => {
                allocated = Some((Pending(p), file));
                break;
            }
            Err(e) if e.kind() == std::io::ErrorKind::AlreadyExists => (),
            Err(e) => return Err(io(e)),
        }
    }
    let (pending, mut file) = allocated.ok_or(HostError::Unavailable)?;
    file.write_all(&bytes).map_err(io)?;
    file.sync_all().map_err(io)?;
    drop(file);
    fs::hard_link(&pending.0, path).map_err(io)?;
    sync(parent)?;
    drop(pending);
    sync(&staging)?;
    Ok(())
}
fn file_hash(path: &Path, max: u64) -> HostResult<(base::Digest, u64)> {
    private(path, false)?;
    let mut f = File::open(path).map_err(io)?;
    let len = f.metadata().map_err(io)?.len();
    if len == 0 || len > max {
        return Err(HostError::Corrupt);
    }
    let mut h = blake3::Hasher::new();
    let mut b = [0u8; 64 * 1024];
    let mut total = 0u64;
    loop {
        let n = f.read(&mut b).map_err(io)?;
        if n == 0 {
            break;
        }
        total = total
            .checked_add(n as u64)
            .filter(|n| *n <= max)
            .ok_or(HostError::Corrupt)?;
        h.update(&b[..n]);
    }
    if total != len {
        return Err(HostError::Corrupt);
    }
    Ok((
        base::Digest::new(h.finalize().to_hex().as_str()).map_err(io)?,
        len,
    ))
}
impl CatalogBackupService {
    pub fn open(root: impl AsRef<Path>, local_lab: &store::LabId) -> HostResult<Self> {
        let root = permissions::root(root.as_ref())?;
        let lease = CatalogLease::acquire(&root)?;
        Self::open_with_lease(&root, local_lab, lease)
    }
    pub(crate) fn open_with_lease(
        root: impl AsRef<Path>,
        local_lab: &store::LabId,
        lease: std::sync::Arc<CatalogLease>,
    ) -> HostResult<Self> {
        let root = permissions::root(root.as_ref())?;
        let database = current_catalog_path(&root)?;
        private(&database, false)?;
        if fs::metadata(&database).map_err(io)?.len() == 0 {
            return Err(invalid(
                "Save a local catalog before creating its backup service.",
            ));
        }

        {
            let c = rusqlite::Connection::open_with_flags(
                &database,
                rusqlite::OpenFlags::SQLITE_OPEN_READ_ONLY,
            )
            .map_err(io)?;
            c.execute_batch("PRAGMA trusted_schema=OFF; PRAGMA query_only=ON")
                .map_err(io)?;
            store::validate_database_header(
                c.query_row("PRAGMA application_id", [], |r| r.get::<_, i64>(0))
                    .map_err(io)?,
                c.query_row("PRAGMA user_version", [], |r| r.get::<_, i64>(0))
                    .map_err(io)?,
            )
            .map_err(|_| HostError::Corrupt)?;
            let mut q = c
                .prepare("SELECT id FROM labs ORDER BY id LIMIT 2")
                .map_err(io)?;
            let labs = q
                .query_map([], |r| r.get::<_, String>(0))
                .map_err(io)?
                .collect::<std::result::Result<Vec<_>, _>>()
                .map_err(io)?;
            if labs.len() != 1 || labs[0] != local_lab.as_str() {
                return Err(invalid(
                    "The saved catalog has a different complete local lab scope.",
                ));
            }
        }
        let mut library = store::Library::open(&database).map_err(io)?;
        library.retain_native_owner(lease).map_err(io)?;
        library.trusted_host_access(local_lab).map_err(io)?;
        let backups = root.join("backups");
        mkdir(&backups)?;
        let operations = backups.join("operations");
        mkdir(&operations)?;
        mkdir(&operations.join("staging"))?;
        let recoveries = root.join("recoveries");
        mkdir(&recoveries)?;
        let service = Self {
            root,
            database,
            backups,
            operations,
            recoveries,
            lab: local_lab.clone(),
            library,
        };
        service.intents()?;
        Ok(service)
    }
    fn scope(&self) -> HostResult<()> {
        private(&self.root, true)?;
        private(&self.backups, true)?;
        private(&self.operations, true)?;
        private(&self.recoveries, true)?;
        private(&self.database, false)
    }
    fn policy(&self) -> backup::Policy {
        backup::Policy::single_lab(self.lab.clone())
    }
    fn wire_lab(&self) -> HostResult<base::LabId> {
        base::LabId::new(self.lab.as_str()).map_err(io)
    }
    fn artifact(&self, id: &str) -> HostResult<PathBuf> {
        id::check(id)?;
        Ok(self.backups.join(format!("{id}.cglbak")))
    }
    fn intent_path(&self, op: &str) -> PathBuf {
        self.operations.join(format!("{op}.intent.json"))
    }
    fn receipt_path(&self, op: &str) -> PathBuf {
        self.operations.join(format!("{op}.receipt.json"))
    }
    fn validate_intent(&self, intent: &Intent, op: &str) -> HostResult<()> {
        id::check(op)?;
        if intent.version != 1 || intent.operation_id != op || intent.lab != self.wire_lab()? {
            return Err(HostError::Corrupt);
        }
        if let Action::Restore { backup_id, .. } = &intent.action {
            id::check(backup_id)?;
        }
        Ok(())
    }
    fn intents(&self) -> HostResult<Vec<Intent>> {
        self.scope()?;
        let mut intents = vec![];
        let mut count = 0usize;
        let mut receipt_ids = vec![];
        for e in fs::read_dir(&self.operations).map_err(io)? {
            count += 1;
            if count > MAX_OPERATIONS * 2 + 1 {
                return Err(invalid("Managed backup operation bound reached."));
            }
            let e = e.map_err(io)?;
            let name = e
                .file_name()
                .into_string()
                .map_err(|_| HostError::Corrupt)?;
            if name == "staging" {
                private(&e.path(), true)?;
                continue;
            }
            private(&e.path(), false)?;
            if let Some(op) = name.strip_suffix(".intent.json") {
                let intent: Intent = read(&e.path())?;
                self.validate_intent(&intent, op)?;
                intents.push(intent)
            } else if let Some(op) = name.strip_suffix(".receipt.json") {
                id::check(op)?;
                receipt_ids.push(op.to_string());
            } else {
                return Err(invalid(
                    "Unrecognized file in managed backup operation storage.",
                ));
            }
        }
        if intents.len() > MAX_OPERATIONS {
            return Err(invalid("Managed backup operation bound reached."));
        }
        intents.sort_by(|a, b| a.operation_id.cmp(&b.operation_id));
        for op in receipt_ids {
            let intent = intents
                .iter()
                .find(|i| i.operation_id == op)
                .ok_or_else(|| invalid("Managed backup receipt has no owned durable intent."))?;
            match &intent.action {
                Action::Create { request_checkpoint } => {
                    let r: Receipt = read(&self.receipt_path(&op))?;
                    if r.id != op
                        || r.basename != format!("{op}.cglbak")
                        || r.lab != self.wire_lab()?
                        || &r.request_checkpoint != request_checkpoint
                        || r.lab_count != 1
                    {
                        return Err(HostError::Corrupt);
                    }
                }
                Action::Restore {
                    backup_id,
                    archive_hash,
                } => {
                    let r: RecoveryReceipt = read(&self.receipt_path(&op))?;
                    if r.operation_id != op
                        || r.directory_basename != op
                        || r.database_basename != "library.sqlite"
                        || &r.backup.id != backup_id
                        || &r.backup.archive_hash != archive_hash
                        || r.backup.lab != self.wire_lab()?
                        || r.backup.lab_count != 1
                    {
                        return Err(HostError::Corrupt);
                    }
                }
            }
        }
        let mut count = 0usize;
        for e in fs::read_dir(&self.recoveries).map_err(io)? {
            count += 1;
            if count > MAX_OPERATIONS + 64 {
                return Err(invalid("Managed recovery directory bound reached."));
            }
            let e = e.map_err(io)?;
            let op = e
                .file_name()
                .into_string()
                .map_err(|_| HostError::Corrupt)?;
            if op.starts_with(".collagene-backup-") {
                private(&e.path(), true)?;
                continue;
            }
            id::check(&op)?;
            private(&e.path(), true)?;
            if !intents
                .iter()
                .any(|i| i.operation_id == op && matches!(i.action, Action::Restore { .. }))
            {
                return Err(invalid("Recovery directory has no owned durable intent."));
            }
        }
        let mut count = 0usize;
        for e in fs::read_dir(&self.backups).map_err(io)? {
            count += 1;
            if count > MAX_OPERATIONS + 65 {
                return Err(invalid("Managed backup directory bound reached."));
            }
            let e = e.map_err(io)?;
            let name = e
                .file_name()
                .into_string()
                .map_err(|_| HostError::Corrupt)?;
            if name == "operations" {
                private(&e.path(), true)?;
                continue;
            }
            if name.starts_with(".collagene-backup-") {
                private(&e.path(), true)?;
                continue;
            }
            let op = name
                .strip_suffix(".cglbak")
                .ok_or_else(|| invalid("Unrecognized file in managed backup storage."))?;
            id::check(op)?;
            private(&e.path(), false)?;
            if !intents
                .iter()
                .any(|i| i.operation_id == op && matches!(i.action, Action::Create { .. }))
            {
                return Err(invalid("Backup file has no owned durable creation intent."));
            }
        }
        Ok(intents)
    }
    fn request_checkpoint(&self) -> HostResult<Option<base::Digest>> {
        let c = rusqlite::Connection::open_with_flags(
            &self.database,
            rusqlite::OpenFlags::SQLITE_OPEN_READ_ONLY,
        )
        .map_err(io)?;
        let exists:bool=c.query_row("SELECT EXISTS(SELECT 1 FROM sqlite_master WHERE type='table' AND name='cg_reliable_checkpoints')",[],|r|r.get(0)).map_err(io)?;
        if !exists {
            return Ok(None);
        }
        let value:Option<String>=c.query_row("SELECT digest FROM cg_reliable_checkpoints WHERE lab_id=?1 AND namespace='repository' AND key='state'",[self.lab.as_str()],|r|r.get(0)).optional().map_err(io)?;
        value
            .as_ref()
            .map(base::Digest::new)
            .transpose()
            .map_err(io)
    }
    fn new_intent(&self, intent: &Intent) -> HostResult<()> {
        if activation::activation_operation_in_use(&self.root, &intent.operation_id)? {
            return Err(invalid("Operation ID already belongs to catalog adoption."));
        }
        let all = self.intents()?;
        if all.len() >= MAX_OPERATIONS {
            return Err(invalid("Managed backup operation bound reached."));
        }
        write(&self.intent_path(&intent.operation_id), intent)
    }
    fn verified(&self, op: &str, intent: &Intent) -> HostResult<Receipt> {
        self.validate_intent(intent, op)?;
        let Action::Create { request_checkpoint } = &intent.action else {
            return Err(invalid(
                "This operation ID belongs to recovery, not backup creation.",
            ));
        };
        let path = self.artifact(op)?;
        let manifest = backup::verify(&path, &self.policy()).map_err(|_| HostError::Corrupt)?;
        let (archive_hash, archive_bytes) =
            file_hash(&path, backup::DEFAULT_MAX_DATABASE_BYTES + 65536 + 44)?;
        Ok(Receipt {
            id: op.into(),
            basename: format!("{op}.cglbak"),
            archive_hash,
            archive_bytes,
            created_unix_seconds: manifest.created_unix_seconds,
            database_hash: convert_digest(&manifest.database_hash)?,
            database_bytes: manifest.database_bytes,
            schema_hash: convert_digest(&manifest.schema_hash)?,
            schema_version: manifest.schema_version.try_into().map_err(io)?,
            lab: self.wire_lab()?,
            lab_count: manifest.labs.len().try_into().map_err(io)?,
            documents: *manifest
                .table_rows
                .get("documents")
                .ok_or(HostError::Corrupt)?,
            folders: *manifest
                .table_rows
                .get("folders")
                .ok_or(HostError::Corrupt)?,
            revisions: *manifest
                .table_rows
                .get("revisions")
                .ok_or(HostError::Corrupt)?,
            mutation_receipts: manifest
                .table_rows
                .get("cg_reliable_mutations")
                .copied()
                .unwrap_or(0),
            request_checkpoint: request_checkpoint.clone(),
        })
    }
    /// An existing creation intent never authorizes taking a later replacement
    /// snapshot. Only its already-published immutable artifact can complete retry.
    pub fn create(&mut self, op: &OperationId) -> HostResult<Receipt> {
        check_catalog_connection(&self.root, &self.database)?;
        self.scope()?;
        let op = op.as_str();
        let intent_path = self.intent_path(op);
        let artifact = self.artifact(op)?;
        let intent = if exists(&intent_path)? {
            let i: Intent = read(&intent_path)?;
            self.validate_intent(&i, op)?;
            if !matches!(i.action, Action::Create { .. }) {
                return Err(invalid(
                    "Operation ID already belongs to a different backup action.",
                ));
            }
            if !exists(&artifact)? {
                return Err(invalid(
                    "This interrupted backup has no published snapshot. Use a new operation ID; retry cannot silently take a newer snapshot.",
                ));
            }
            i
        } else {
            if exists(&artifact)? || exists(&self.receipt_path(op))? {
                return Err(invalid("Unowned pre-existing backup or receipt."));
            }
            let i = Intent {
                version: 1,
                operation_id: op.into(),
                lab: self.wire_lab()?,
                action: Action::Create {
                    request_checkpoint: self.request_checkpoint()?,
                },
            };
            self.new_intent(&i)?;
            backup::create(&self.library, &artifact, &self.policy()).map_err(io)?;
            i
        };
        let receipt = self.verified(op, &intent)?;
        let path = self.receipt_path(op);
        if exists(&path)? {
            let old: Receipt = read(&path)?;
            if old != receipt {
                return Err(HostError::Corrupt);
            }
        } else {
            write(&path, &receipt)?;
        }
        Ok(receipt)
    }
    pub fn verify(&self, backup_id: &str) -> HostResult<Receipt> {
        self.verify_pinned(backup_id, None)
    }
    fn verify_pinned(
        &self,
        backup_id: &str,
        expected: Option<&base::Digest>,
    ) -> HostResult<Receipt> {
        self.scope()?;
        id::check(backup_id)?;
        let intent: Intent = read(&self.intent_path(backup_id))?;
        let receipt = self.verified(backup_id, &intent)?;
        if expected.is_some_and(|h| h != &receipt.archive_hash) {
            return Err(invalid(
                "The selected backup changed; verify and choose it again.",
            ));
        }
        let path = self.receipt_path(backup_id);
        if exists(&path)? {
            let recorded: Receipt = read(&path)?;
            if recorded != receipt {
                return Err(HostError::Corrupt);
            }
        } else {
            write(&path, &receipt)?;
        }
        Ok(receipt)
    }
    pub fn list(
        &self,
        offset: u32,
        limit: u16,
        expected_catalog: Option<&base::Digest>,
    ) -> HostResult<Page> {
        if limit == 0 || limit > MAX_PAGE || offset as usize > MAX_OPERATIONS {
            return Err(invalid("Invalid managed backup page."));
        }
        let mut entries = vec![];
        for intent in self.intents()? {
            if !matches!(intent.action, Action::Create { .. }) {
                continue;
            }
            let op = &intent.operation_id;
            let receipt = if exists(&self.receipt_path(op))? {
                let r: Receipt = read(&self.receipt_path(op))?;
                if r.id != *op
                    || r.basename != format!("{op}.cglbak")
                    || r.lab != self.wire_lab()?
                    || r.request_checkpoint
                        != match intent.action {
                            Action::Create {
                                ref request_checkpoint,
                            } => request_checkpoint.clone(),
                            _ => None,
                        }
                {
                    return Err(HostError::Corrupt);
                }
                let artifact = self.artifact(op)?;
                private(&artifact, false)?;
                if fs::metadata(artifact).map_err(io)?.len() != r.archive_bytes {
                    return Err(HostError::Corrupt);
                }
                Some(r)
            } else {
                None
            };
            entries.push(Entry {
                id: op.clone(),
                basename: format!("{op}.cglbak"),
                receipt,
            });
        }
        let catalog_revision = base::Digest::of(&serde_json::to_vec(&entries).map_err(io)?);
        if offset > 0 && expected_catalog != Some(&catalog_revision) {
            return Err(invalid(
                "Managed backup list changed; refresh its first page.",
            ));
        }
        let start = offset as usize;
        if start > entries.len() {
            return Err(invalid("Managed backup page exceeds the catalog."));
        }
        let end = (start + usize::from(limit)).min(entries.len());
        Ok(Page {
            entries: entries[start..end].to_vec(),
            next_offset: (end < entries.len()).then_some(end as u32),
            catalog_revision,
        })
    }
    pub fn restore_copy(
        &mut self,
        op: &OperationId,
        backup_id: &str,
        expected_archive_hash: &base::Digest,
    ) -> HostResult<RecoveryReceipt> {
        check_catalog_connection(&self.root, &self.database)?;
        self.scope()?;
        id::check(backup_id)?;
        let source = self.verify_pinned(backup_id, Some(expected_archive_hash))?;
        let op = op.as_str();
        let path = self.intent_path(op);
        let action = Action::Restore {
            backup_id: backup_id.into(),
            archive_hash: source.archive_hash.clone(),
        };
        let directory = self.recoveries.join(op);
        if exists(&path)? {
            let intent: Intent = read(&path)?;
            self.validate_intent(&intent, op)?;
            if intent.action != action {
                return Err(invalid(
                    "Operation ID already belongs to a different backup recovery.",
                ));
            }
        } else {
            if exists(&directory)? || exists(&self.receipt_path(op))? {
                return Err(invalid("Recovery destination or receipt already exists."));
            }
            self.new_intent(&Intent {
                version: 1,
                operation_id: op.into(),
                lab: self.wire_lab()?,
                action,
            })?;
        }
        if !exists(&directory)? {
            backup::restore_new(&self.artifact(backup_id)?, &directory, &self.policy())
                .map_err(io)?;
        }
        private(&directory, true)?;
        let mut count = 0;
        for e in fs::read_dir(&directory).map_err(io)? {
            let e = e.map_err(io)?;
            if e.file_name() != "library.sqlite" {
                return Err(invalid(
                    "Recovery copy was opened or changed; no overwrite/cutover is permitted.",
                ));
            }
            count += 1;
        }
        if count != 1 {
            return Err(HostError::Corrupt);
        }
        let (hash, bytes) = file_hash(
            &directory.join("library.sqlite"),
            backup::DEFAULT_MAX_DATABASE_BYTES,
        )?;
        if hash != source.database_hash || bytes != source.database_bytes {
            return Err(HostError::Corrupt);
        }
        let receipt = RecoveryReceipt {
            operation_id: op.into(),
            backup: source,
            directory_basename: op.into(),
            database_basename: "library.sqlite".into(),
        };
        let path = self.receipt_path(op);
        if exists(&path)? {
            let old: RecoveryReceipt = read(&path)?;
            if old != receipt {
                return Err(HostError::Corrupt);
            }
        } else {
            write(&path, &receipt)?;
        }
        Ok(receipt)
    }
}
/// Dedicated native raw-IPC lane, independent of the editor/working host mutex.
pub struct DesktopCatalogBackups {
    service: std::sync::Mutex<Option<CatalogBackupService>>,
    root: PathBuf,
    generation: Option<u64>,
    control: HostControl,
    _lease: std::sync::Arc<CatalogLease>,
}
fn backup_failure(error: HostError) -> collagene_desktop_protocol::Failure {
    crate::failure(match error {
        HostError::Corrupt => invalid(
            "Saved-catalog backup or recovery failed native integrity checks; preserve the existing files.",
        ),
        other => other,
    })
}
impl CatalogBackupService {
    fn verify_recovery(
        &self,
        operation: &OperationId,
        backup_id: &str,
        expected: &base::Digest,
    ) -> HostResult<RecoveryReceipt> {
        let source = self.verify_pinned(backup_id, Some(expected))?;
        let op = operation.as_str();
        let intent: Intent = read(&self.intent_path(op))?;
        self.validate_intent(&intent, op)?;
        if intent.action
            != (Action::Restore {
                backup_id: backup_id.into(),
                archive_hash: expected.clone(),
            })
        {
            return Err(HostError::Corrupt);
        }
        let directory = self.recoveries.join(op);
        private(&directory, true)?;
        let mut entries = fs::read_dir(&directory).map_err(io)?;
        let only = entries.next().ok_or(HostError::Corrupt)?.map_err(io)?;
        if only.file_name() != "library.sqlite" || entries.next().is_some() {
            return Err(HostError::Corrupt);
        }
        let (hash, bytes) = file_hash(&only.path(), backup::DEFAULT_MAX_DATABASE_BYTES)?;
        if hash != source.database_hash || bytes != source.database_bytes {
            return Err(HostError::Corrupt);
        }
        let expected = RecoveryReceipt {
            operation_id: op.into(),
            backup: source,
            directory_basename: op.into(),
            database_basename: "library.sqlite".into(),
        };
        let recorded: RecoveryReceipt = read(&self.receipt_path(op))?;
        if recorded != expected {
            return Err(HostError::Corrupt);
        }
        Ok(recorded)
    }
}
impl DesktopCatalogBackups {
    /// Fresh durable admission for this exact process generation. Historical
    /// retries and rejected requests must not themselves close a healthy session.
    pub fn catalog_session_ready(&self) -> bool {
        matches!(
            catalog_session_status(&self.root, self.generation),
            collagene_desktop_protocol::CatalogStatus::Ready { .. }
        )
    }
    pub fn open(
        root: impl AsRef<Path>,
        lab: &store::LabId,
        control: HostControl,
    ) -> HostResult<Self> {
        let root = permissions::root(root.as_ref())?;
        let lease = CatalogLease::acquire(&root)?;
        Self::open_with_lease(root, lab, control, lease)
    }
    pub(crate) fn open_with_lease(
        root: impl AsRef<Path>,
        lab: &store::LabId,
        control: HostControl,
        lease: std::sync::Arc<CatalogLease>,
    ) -> HostResult<Self> {
        let root = permissions::root(root.as_ref())?;
        let generation = Some(catalog_generation(&root)?);
        let service = CatalogBackupService::open_with_lease(&root, lab, lease.clone())?;
        Ok(Self {
            service: std::sync::Mutex::new(Some(service)),
            root,
            generation,
            control,
            _lease: lease,
        })
    }
    /// Root-only control/status/cancellation remains available even when no catalog
    /// connection or biological scope can be trusted. All data-bearing commands refuse.
    pub(crate) fn open_recovery_with_lease(
        root: impl AsRef<Path>,
        lab: Option<&store::LabId>,
        control: HostControl,
        lease: std::sync::Arc<CatalogLease>,
        generation: Option<u64>,
    ) -> HostResult<Self> {
        let root = permissions::root(root.as_ref())?;
        let service = lab
            .and_then(|lab| CatalogBackupService::open_with_lease(&root, lab, lease.clone()).ok());
        Ok(Self {
            service: std::sync::Mutex::new(service),
            root,
            generation,
            control,
            _lease: lease,
        })
    }
    fn execute(
        &self,
        operation: &OperationId,
        command: &collagene_desktop_protocol::CatalogBackupCommand,
    ) -> HostResult<collagene_desktop_protocol::CatalogBackupSuccess> {
        use collagene_desktop_protocol::{CatalogBackupCommand as C, CatalogBackupSuccess as S};
        match command {
            C::CancelActivation {
                accepted_operation,
                review_revision,
                activation_generation,
            } => {
                return Ok(S::ActivationCancelled(activation::cancel_activation(
                    &self.root,
                    operation,
                    accepted_operation,
                    review_revision,
                    *activation_generation,
                )?));
            }
            C::AdoptRecovery { review }
                if !matches!(
                    catalog_session_status(&self.root, self.generation),
                    collagene_desktop_protocol::CatalogStatus::Ready { .. }
                ) =>
            {
                let prior=activation::replay_adoption(&self.root,operation,review)?.ok_or_else(||invalid("Restart before accepting another catalog selection; only an exact durable retry is allowed in this session."))?;
                return Ok(S::AdoptionPending(prior));
            }
            C::ActivationStatus => {
                let status = activation::activation_status(&self.root)?;
                if matches!(
                    status,
                    collagene_desktop_protocol::CatalogActivationStatus::Cancelled(_)
                ) {
                    return Err(invalid(
                        "Use recovery bootstrap to inspect the cancelled activation.",
                    ));
                }
                return Ok(S::ActivationStatus(status));
            }
            _ => {}
        }
        let mut holder = self.service.lock().map_err(io)?;
        let service=holder.as_mut().ok_or_else(||invalid("Saved catalog is unavailable; only recovery status and exact pending cancellation are permitted."))?;
        if !matches!(command, C::AdoptRecovery { .. } | C::ListRecoveries { .. })
            && !matches!(
                catalog_session_status(&self.root, self.generation),
                collagene_desktop_protocol::CatalogStatus::Ready { .. }
            )
        {
            return Err(invalid(
                "Catalog session is paused. Restart after resolving activation before using catalog data.",
            ));
        }
        if !matches!(
            command,
            C::AdoptRecovery { .. } | C::ActivationStatus | C::ListRecoveries { .. }
        ) {
            check_catalog_connection(&service.root, &service.database)?;
        }
        Ok(match command {
            C::ListRecoveries {
                offset,
                limit,
                expected_catalog,
            } => S::RecoveriesListed(service.list_recoveries(
                *offset,
                *limit,
                expected_catalog.as_ref(),
            )?),
            C::ReviewOriginal => S::Reviewed(service.review_original()?),
            C::ReviewRecovery { id, expected_hash } => {
                S::Reviewed(service.review_recovery(id, expected_hash)?)
            }
            C::AdoptRecovery { review } => {
                S::AdoptionPending(service.adopt_recovery(operation, review)?)
            }
            C::ActivationStatus | C::CancelActivation { .. } => {
                unreachable!("handled root-only commands")
            }
            C::Create => S::Created(crate::convert(&service.create(operation)?)?),
            C::List {
                offset,
                limit,
                expected_catalog,
            } => S::Listed(crate::convert(&service.list(
                *offset,
                *limit,
                expected_catalog.as_ref(),
            )?)?),
            C::Verify { id, expected_hash } => S::Verified(crate::convert(
                &service.verify_pinned(id.as_str(), expected_hash.as_ref())?,
            )?),
            C::RestoreCopy { id, expected_hash } => S::Restored(crate::convert(
                &service.restore_copy(operation, id.as_str(), expected_hash)?,
            )?),
        })
    }
    fn revalidate_cached(
        &self,
        request: &collagene_desktop_protocol::RequestEnvelope,
        reply: &[u8],
    ) -> HostResult<()> {
        use collagene_desktop_protocol::{
            self as wire, CatalogBackupCommand as C, CatalogBackupSuccess as S,
        };
        let response = wire::decode_response_for(reply, request).map_err(|_| HostError::Wire)?;
        let wire::Outcome::Success(wire::Success::CatalogBackup(recorded)) = &response.outcome
        else {
            return Ok(());
        };
        let wire::Command::CatalogBackup(command) = &request.command else {
            return Err(HostError::Wire);
        };
        if let (
            C::CancelActivation {
                accepted_operation,
                review_revision,
                activation_generation,
            },
            S::ActivationCancelled(_),
        ) = (command, recorded)
        {
            let fresh = S::ActivationCancelled(activation::cancel_activation(
                &self.root,
                &request.operation_id,
                accepted_operation,
                review_revision,
                *activation_generation,
            )?);
            return if &fresh == recorded {
                Ok(())
            } else {
                Err(HostError::Corrupt)
            };
        }
        if let (C::AdoptRecovery { review }, S::AdoptionPending(r)) = (command, recorded) {
            let prior = activation::replay_adoption(&self.root, &request.operation_id, review)?
                .ok_or(HostError::Corrupt)?;
            return if &prior == r {
                Ok(())
            } else {
                Err(HostError::Corrupt)
            };
        }
        let holder = self.service.lock().map_err(io)?;
        let service = holder
            .as_ref()
            .ok_or_else(|| invalid("Saved catalog remains paused; use recovery bootstrap."))?;
        let fresh = match (command, recorded) {
            (C::AdoptRecovery { review }, S::AdoptionPending(_)) => {
                S::AdoptionPending(service.adopt_recovery(&request.operation_id, review)?)
            }
            (C::ListRecoveries { .. }, S::RecoveriesListed(_)) => return Ok(()),
            (C::ActivationStatus, S::ActivationStatus(_)) => {
                S::ActivationStatus(service.activation_status()?)
            }
            (C::ReviewOriginal, S::Reviewed(_)) => S::Reviewed(service.review_original()?),
            (C::ReviewRecovery { id, expected_hash }, S::Reviewed(_)) => {
                S::Reviewed(service.review_recovery(id, expected_hash)?)
            }
            (C::Create, S::Created(_)) => S::Created(crate::convert(
                &service.verify(request.operation_id.as_str())?,
            )?),
            (C::Verify { id, expected_hash }, S::Verified(_)) => S::Verified(crate::convert(
                &service.verify_pinned(id.as_str(), expected_hash.as_ref())?,
            )?),
            (C::RestoreCopy { id, expected_hash }, S::Restored(_)) => S::Restored(crate::convert(
                &service.verify_recovery(&request.operation_id, id.as_str(), expected_hash)?,
            )?),
            (C::List { .. }, S::Listed(_)) => return Ok(()),
            _ => return Err(HostError::Wire),
        };
        if &fresh != recorded {
            return Err(HostError::Corrupt);
        }
        Ok(())
    }
    pub fn dispatch_raw(&self, access: &WindowAccess, bytes: &[u8]) -> HostResult<Vec<u8>> {
        use collagene_desktop_protocol as wire;
        self.control.check_window(access)?;
        if bytes.len() > wire::MAX_CATALOG_BACKUP_REQUEST_BYTES {
            return Err(HostError::Wire);
        }
        let request = wire::decode_request(bytes).map_err(|_| HostError::Wire)?;
        let wire::Command::CatalogBackup(command) = &request.command else {
            return Err(HostError::Wire);
        };
        let job =
            match self
                .control
                .begin_job(access, &request.operation_id, base::Digest::of(bytes))
            {
                Ok(JobStart::Run(job)) => job,
                Ok(JobStart::Cached(reply)) => {
                    if reply.len() > wire::MAX_CATALOG_BACKUP_RESPONSE_BYTES {
                        return Err(HostError::Wire);
                    }
                    match self.revalidate_cached(&request, &reply) {
                        Ok(()) => return Ok(reply),
                        Err(error) => {
                            let response = wire::ResponseEnvelope {
                                request_id: request.request_id,
                                operation_id: request.operation_id.clone(),
                                outcome: wire::Outcome::Failure(backup_failure(error)),
                            };
                            response
                                .validate_for(&request)
                                .map_err(|_| HostError::Wire)?;
                            return wire::encode_response_for_frame(&response, bytes)
                                .map_err(|_| HostError::Wire);
                        }
                    }
                }
                Err(error) => {
                    let response = wire::ResponseEnvelope {
                        request_id: request.request_id,
                        operation_id: request.operation_id.clone(),
                        outcome: wire::Outcome::Failure(backup_failure(error)),
                    };
                    response
                        .validate_for(&request)
                        .map_err(|_| HostError::Wire)?;
                    return wire::encode_response_for_frame(&response, bytes)
                        .map_err(|_| HostError::Wire);
                }
            };
        let result = if matches!(
            command,
            wire::CatalogBackupCommand::ListRecoveries { .. }
                | wire::CatalogBackupCommand::List { .. }
                | wire::CatalogBackupCommand::ActivationStatus
                | wire::CatalogBackupCommand::ReviewOriginal
                | wire::CatalogBackupCommand::ReviewRecovery { .. }
        ) {
            self.execute(&request.operation_id, command)
        } else {
            job.accept()
                .and_then(|()| self.execute(&request.operation_id, command))
        };
        let outcome = if job.cancelled_before_acceptance() {
            wire::Outcome::Failure(wire::Failure::Cancelled)
        } else {
            match result {
                Ok(success) => wire::Outcome::Success(wire::Success::CatalogBackup(success)),
                Err(error) => wire::Outcome::Failure(backup_failure(error)),
            }
        };
        let terminal = match &outcome {
            wire::Outcome::Success(_) => wire::OperationState::Completed,
            wire::Outcome::Failure(wire::Failure::Cancelled) => wire::OperationState::Cancelled,
            _ => wire::OperationState::Failed,
        };
        let response = wire::ResponseEnvelope {
            request_id: request.request_id,
            operation_id: request.operation_id.clone(),
            outcome,
        };
        response
            .validate_for(&request)
            .map_err(|_| HostError::Wire)?;
        let reply =
            wire::encode_response_for_frame(&response, bytes).map_err(|_| HostError::Wire)?;
        if reply.len() > wire::MAX_CATALOG_BACKUP_RESPONSE_BYTES {
            return Err(HostError::Wire);
        }
        // Accepted file publication is authoritative even if the caller closes its
        // window before in-memory status/cache bookkeeping can be retained.
        let _ = job.finish(terminal, &reply);
        Ok(reply)
    }
}

mod id {
    use super::*;
    pub fn check(value: &str) -> HostResult<()> {
        parse_operation(value).map(|_| ())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use collagene_core::{
        Molecule, Topology,
        analysis::Primer,
        cloning_outputs,
        run_lineage::{EventBinding, InputEventBinding},
        run_material,
        verification::{DesignPlan, Operation, PlanInput, Source},
        workspace::{self, Workspace},
    };
    use std::sync::atomic::{AtomicU64, Ordering};
    use store::mutations::{Applied, MutationOperation};
    struct Temp(PathBuf);
    impl Temp {
        fn new() -> Self {
            static NEXT: AtomicU64 = AtomicU64::new(0);
            Self(std::env::temp_dir().join(format!(
                    "collagene-catalog-backup-{}-{}-{}",
                    std::process::id(),
                    std::time::SystemTime::now()
                        .duration_since(std::time::UNIX_EPOCH)
                        .unwrap()
                        .as_nanos(),
                    NEXT.fetch_add(1, Ordering::Relaxed)
                )))
        }
    }
    impl Drop for Temp {
        fn drop(&mut self) {
            let _ = fs::remove_dir_all(&self.0);
        }
    }
    fn op(n: u64) -> OperationId {
        OperationId::new(format!("op-{n:032x}")).unwrap()
    }
    fn event(w: &Workspace, n: usize) -> EventBinding {
        let d = &w.documents[n];
        EventBinding {
            document_id: d.id.clone(),
            event_id: d.head_event.clone().unwrap(),
            revision: d.head.clone(),
        }
    }
    fn native() -> (store::Payload, run_material::Selection) {
        let source = Molecule::dna(
            "Original physical PCR template",
            "ACGACCCCGGGGTAAC",
            Topology::Linear,
        )
        .unwrap();
        let target = Molecule::dna(
            "Independent literal expected PCR",
            "ACGACCCCGGGGTAAC",
            Topology::Linear,
        )
        .unwrap();
        let plan = DesignPlan {
            name: "Real saved native PCR".into(),
            target_revision: target.revision(),
            inputs: vec![PlanInput {
                id: "original".into(),
                revision: source.revision(),
                molecule: source.clone(),
            }],
            steps: vec![Operation::Pcr {
                template: Source::Input("original".into()),
                forward: Primer::new("F", "ACGACCCC", "").unwrap(),
                reverse: Primer::new("R", "GTTACCCC", "").unwrap(),
                template_strand: collagene_core::chemistry::Strand::Top,
            }],
        };
        let w = Workspace::from_documents(
            "Complete private saved catalog archive",
            vec![source, target],
            1,
        )
        .unwrap();
        let commit = cloning_outputs::finalize(
            &w,
            &plan,
            event(&w, 1),
            vec![InputEventBinding {
                input_id: "original".into(),
                event: event(&w, 0),
            }],
            &cloning_outputs::Options {
                name: "Native saved PCR product".into(),
                show_prepared: true,
                save_primers: true,
                allow_incomplete: false,
            },
        )
        .unwrap();
        let selection = run_material::Selection {
            run_id: commit.receipt.run.run_id.clone(),
            step_index: 0,
            event: commit.receipt.construct.clone(),
        };
        let mut frame = commit.workspace;
        let private = b"private unfinished saved-catalog design fixture".to_vec();
        frame.attachments.insert(
            blake3::hash(&private).to_hex().to_string(),
            workspace::WorkspaceAttachment {
                kind: "private-catalog-fixture-v1".into(),
                name: "Unpublished saved design state".into(),
                media_type: "application/octet-stream".into(),
                bytes: private,
            },
        );
        let bytes = workspace::save_workspace(&frame).unwrap();
        let payload = store::Payload {
            kind: store::PayloadKind::WorkspaceArchive,
            archive_hash: store::Digest::parse(blake3::hash(&bytes).to_hex().to_string()).unwrap(),
            bytes,
            selection: Some(store::NativeSelection {
                document_id: selection.event.document_id.clone(),
                revision: store::Digest::parse(selection.event.revision.clone()).unwrap(),
                head_event: Some(selection.event.event_id.clone()),
            }),
            sharing: Some(store::WorkspaceSharing::ExplicitWholeWorkspace),
        };
        (payload, selection)
    }
    fn setup(t: &Temp) -> (collagene_repository::Repository, store::LabId) {
        let root = permissions::root(&t.0).unwrap();
        let path = permissions::database(&root, "repository.sqlite").unwrap();
        let repo = collagene_repository::Repository::open(path).unwrap();
        let lab = store::LabId::parse(repo.local_context().unwrap().lab_id.as_str()).unwrap();
        (repo, lab)
    }
    #[test]
    fn real_saved_catalog_history_private_state_restart_and_recovery_copy_do_not_touch_working_autosave()
     {
        let t = Temp::new();
        let (mut repo, lab) = setup(&t);
        let (payload, selection) = native();
        let before = payload.bytes.clone();
        let Applied::Imported(doc) = repo
            .mutate_local(MutationOperation::ImportDocument {
                folder: None,
                name: "Private plasmid archive".into(),
                payload,
            })
            .unwrap()
        else {
            panic!("import")
        };
        let sentinel = t.0.join("working.sqlite");
        let mut options = OpenOptions::new();
        options.write(true).create_new(true);
        #[cfg(unix)]
        {
            use std::os::unix::fs::OpenOptionsExt;
            options.mode(0o600);
        }
        options
            .open(&sentinel)
            .unwrap()
            .write_all(b"not a SQLite database; do not open editor storage")
            .unwrap();
        let mut service = CatalogBackupService::open(&t.0, &lab).unwrap();
        let receipt = service.create(&op(1)).unwrap();
        assert_eq!(receipt.lab_count, 1);
        assert_eq!(receipt.documents, 1);
        assert_eq!(receipt.revisions, 1);
        assert!(receipt.request_checkpoint.is_some());
        assert!(receipt.created_unix_seconds > 0);
        assert!(
            !serde_json::to_string(&receipt)
                .unwrap()
                .contains(t.0.to_str().unwrap())
        );
        repo.mutate_local(MutationOperation::CreateFolder {
            parent: None,
            name: "Newer private live branch".into(),
        })
        .unwrap();
        assert_eq!(
            service.create(&op(1)).unwrap(),
            receipt,
            "retry acknowledges original bytes, not a newer snapshot"
        );
        drop(service);
        let mut restarted = CatalogBackupService::open(&t.0, &lab).unwrap();
        assert_eq!(restarted.verify(&receipt.id).unwrap(), receipt);
        let recovery = restarted
            .restore_copy(&op(2), &receipt.id, &receipt.archive_hash)
            .unwrap();
        assert_eq!(
            restarted
                .restore_copy(&op(2), &receipt.id, &receipt.archive_hash)
                .unwrap(),
            recovery
        );
        assert_eq!(recovery.directory_basename, op(2).as_str());
        assert_eq!(recovery.database_basename, "library.sqlite");
        assert_eq!(
            fs::read(&sentinel).unwrap(),
            b"not a SQLite database; do not open editor storage"
        );
        let recovered = store::Library::open(
            t.0.join("recoveries")
                .join(op(2).as_str())
                .join("library.sqlite"),
        )
        .unwrap();
        let access = recovered.trusted_host_access(&lab).unwrap();
        let bytes = recovered
            .payload(&access, &doc.id, &doc.head.id)
            .unwrap()
            .bytes;
        assert_eq!(bytes, before);
        let frame = workspace::load_workspace(&bytes).unwrap();
        assert!(
            frame
                .attachments
                .values()
                .any(|a| a.kind == "private-catalog-fixture-v1")
        );
        assert_eq!(
            run_material::reconstruct(&frame, &selection)
                .unwrap()
                .reference()
                .sequence
                .as_str(),
            "ACGACCCCGGGGTAAC"
        );
        assert_eq!(repo.document(&doc.id).unwrap(), doc);
        assert_eq!(
            fs::read(&sentinel).unwrap(),
            b"not a SQLite database; do not open editor storage"
        );
    }
    #[test]
    fn interrupted_creation_and_lost_receipt_never_take_a_replacement_snapshot_and_recovery_hash_is_pinned()
     {
        let t = Temp::new();
        let (mut repo, lab) = setup(&t);
        let mut service = CatalogBackupService::open(&t.0, &lab).unwrap();
        let first = service.create(&op(10)).unwrap();
        fs::remove_file(service.receipt_path(op(10).as_str())).unwrap();
        repo.mutate_local(MutationOperation::CreateFolder {
            parent: None,
            name: "After snapshot".into(),
        })
        .unwrap();
        drop(service);
        let mut service = CatalogBackupService::open(&t.0, &lab).unwrap();
        assert_eq!(service.create(&op(10)).unwrap(), first);
        assert_eq!(service.verify(&first.id).unwrap(), first);
        let wrong = base::Digest::of(b"wrong selected backup bytes");
        assert!(service.restore_copy(&op(11), &first.id, &wrong).is_err());
        assert!(!service.intent_path(op(11).as_str()).exists());
        assert!(!service.recoveries.join(op(11).as_str()).exists());
        assert!(
            service
                .restore_copy(&op(10), &first.id, &first.archive_hash)
                .is_err()
        );
        let interrupted = Intent {
            version: 1,
            operation_id: op(12).as_str().into(),
            lab: service.wire_lab().unwrap(),
            action: Action::Create {
                request_checkpoint: service.request_checkpoint().unwrap(),
            },
        };
        service.new_intent(&interrupted).unwrap();
        assert!(service.create(&op(12)).is_err());
        assert!(!service.artifact(op(12).as_str()).unwrap().exists());
        let page = service.list(0, 200, None).unwrap();
        assert!(
            page.entries
                .iter()
                .any(|e| e.id == op(12).as_str() && e.receipt.is_none())
        );
        let newer = service.create(&op(13)).unwrap();
        assert_eq!(newer.folders, 1);
        assert_ne!(newer.archive_hash, first.archive_hash);
    }
    #[test]
    fn complete_checkpoint_pagination_exposes_the_201st_owned_pending_intent_and_rejects_stale_pages()
     {
        let t = Temp::new();
        let (_repo, lab) = setup(&t);
        let service = CatalogBackupService::open(&t.0, &lab).unwrap();
        for n in 1..=201 {
            service
                .new_intent(&Intent {
                    version: 1,
                    operation_id: op(n).as_str().into(),
                    lab: service.wire_lab().unwrap(),
                    action: Action::Create {
                        request_checkpoint: None,
                    },
                })
                .unwrap();
        }
        let a = service.list(0, 200, None).unwrap();
        assert_eq!(a.entries.len(), 200);
        assert_eq!(a.next_offset, Some(200));
        let b = service.list(200, 200, Some(&a.catalog_revision)).unwrap();
        assert_eq!(b.entries.len(), 1);
        assert_eq!(b.entries[0].id, op(201).as_str());
        assert!(b.next_offset.is_none());
        assert!(service.list(200, 200, None).is_err());
        service
            .new_intent(&Intent {
                version: 1,
                operation_id: op(202).as_str().into(),
                lab: service.wire_lab().unwrap(),
                action: Action::Create {
                    request_checkpoint: None,
                },
            })
            .unwrap();
        assert!(service.list(200, 200, Some(&a.catalog_revision)).is_err());
        assert!(service.list(0, 0, None).is_err());
        assert!(service.list(0, 201, None).is_err());
    }
    #[test]
    fn unknown_ids_foreign_names_corrupt_archives_and_mismatched_local_lab_fail_closed() {
        let t = Temp::new();
        let (_repo, lab) = setup(&t);
        let mut service = CatalogBackupService::open(&t.0, &lab).unwrap();
        let receipt = service.create(&op(1)).unwrap();
        for value in [
            "../repository.sqlite",
            "op-00000000000000000000000000000001.cglbak",
            "/tmp/a",
            "wrong",
        ] {
            assert!(service.verify(value).is_err());
        }
        let other = store::LabId::parse("lab-11111111111111111111111111111111").unwrap();
        assert!(CatalogBackupService::open(&t.0, &other).is_err());
        let path = service.artifact(&receipt.id).unwrap();
        let mut bytes = fs::read(&path).unwrap();
        let n = bytes.len() - 1;
        bytes[n] ^= 1;
        fs::write(&path, bytes).unwrap();
        assert!(service.verify(&receipt.id).is_err());
        assert!(service.create(&op(1)).is_err());
        assert!(
            service
                .restore_copy(&op(2), &receipt.id, &receipt.archive_hash)
                .is_err()
        );
        assert!(!service.recoveries.join(op(2).as_str()).exists());
        let foreign = service.backups.join("unexpected.bin");
        let mut options = OpenOptions::new();
        options.write(true).create_new(true);
        #[cfg(unix)]
        {
            use std::os::unix::fs::OpenOptionsExt;
            options.mode(0o600);
        }
        options
            .open(&foreign)
            .unwrap()
            .write_all(b"foreign file must not become a backup")
            .unwrap();
        assert!(service.list(0, 200, None).is_err());
        assert!(CatalogBackupService::open(&t.0, &lab).is_err());
    }
    #[cfg(unix)]
    #[test]
    fn managed_directory_symlinks_are_refused_without_reading_or_writing_their_target() {
        use std::os::unix::fs::symlink;
        let t = Temp::new();
        let (_repo, lab) = setup(&t);
        let outside = Temp::new();
        permissions::root(&outside.0).unwrap();
        symlink(&outside.0, t.0.join("backups")).unwrap();
        assert!(CatalogBackupService::open(&t.0, &lab).is_err());
        assert_eq!(fs::read_dir(&outside.0).unwrap().count(), 0);
    }
    #[test]
    fn recovery_adoption_is_exact_idempotent_restart_bound_and_preserves_private_editor_and_original()
     {
        let t = Temp::new();
        let (mut repo, lab) = setup(&t);
        let (payload, selection) = native();
        let original_archive = payload.bytes.clone();
        let Applied::Imported(doc) = repo
            .mutate_local(MutationOperation::ImportDocument {
                folder: None,
                name: "Saved native PCR".into(),
                payload,
            })
            .unwrap()
        else {
            panic!("import")
        };
        let mut service = CatalogBackupService::open(&t.0, &lab).unwrap();
        let backup = service.create(&op(300)).unwrap();
        repo.mutate_local(MutationOperation::CreateFolder {
            parent: None,
            name: "Newer private branch folder".into(),
        })
        .unwrap();
        service
            .restore_copy(&op(301), &backup.id, &backup.archive_hash)
            .unwrap();
        let page = service.list_recoveries(0, 1, None).unwrap();
        assert_eq!(page.entries.len(), 1);
        let review = service
            .review_recovery(&op(301), &backup.archive_hash)
            .unwrap();
        let accepted = service.adopt_recovery(&op(302), &review).unwrap();
        assert_eq!(service.adopt_recovery(&op(302), &review).unwrap(), accepted);
        assert!(service.adopt_recovery(&op(303), &review).is_err());
        assert!(catalog_activation_pending(&t.0).unwrap());
        assert!(check_catalog_connection(&t.0, &t.0.join("repository.sqlite")).is_err());
        assert!(matches!(
            service.activation_status().unwrap(),
            collagene_desktop_protocol::CatalogActivationStatus::PendingRestart(_)
        ));
        // Editor state is a separate namespace; even opaque unsaved bytes are untouched.
        let working = t.0.join("working.sqlite");
        let mut options = OpenOptions::new();
        options.write(true).create_new(true);
        #[cfg(unix)]
        {
            use std::os::unix::fs::OpenOptionsExt;
            options.mode(0o600);
        }
        options
            .open(&working)
            .unwrap()
            .write_all(b"independent editor recovery")
            .unwrap();
        drop(repo);
        drop(service);
        let lease = CatalogLease::acquire(&t.0).unwrap();
        let selected = startup_catalog_path(&t.0, &lease).unwrap();
        assert_eq!(
            selected,
            t.0.join("recoveries")
                .join(op(301).as_str())
                .join("library.sqlite")
        );
        assert_eq!(fs::read(&working).unwrap(), b"independent editor recovery");
        assert!(t.0.join("repository.sqlite").is_file());
        assert!(!catalog_activation_pending(&t.0).unwrap());
        assert!(check_catalog_connection(&t.0, &t.0.join("repository.sqlite")).is_err());
        assert!(check_catalog_connection(&t.0, &selected).is_ok());
        let active = CatalogBackupService::open_with_lease(&t.0, &lab, lease.clone()).unwrap();
        assert!(matches!(
            active.activation_status().unwrap(),
            collagene_desktop_protocol::CatalogActivationStatus::ActiveAfterRestart(_)
        ));
        let access = active.library.trusted_host_access(&lab).unwrap();
        let bytes = active
            .library
            .payload(&access, &doc.id, &doc.head.id)
            .unwrap()
            .bytes;
        assert_eq!(
            bytes, original_archive,
            "native IDs, private archives and run evidence remain exact"
        );
        let frame = workspace::load_workspace(&bytes).unwrap();
        assert_eq!(
            run_material::reconstruct(&frame, &selection)
                .unwrap()
                .reference()
                .sequence
                .as_str(),
            "ACGACCCCGGGGTAAC"
        );
        // Legitimate writes to selected catalog survive a second restart: do not require
        // immutable pre-activation bytes after it has become an ordinary writable catalog.
        let mut restored = collagene_repository::Repository::open(&selected).unwrap();
        restored
            .mutate_local(MutationOperation::CreateFolder {
                parent: None,
                name: "Post-activation folder".into(),
            })
            .unwrap();
        drop(restored);
        drop(active);
        drop(lease);
        let lease = CatalogLease::acquire(&t.0).unwrap();
        assert_eq!(startup_catalog_path(&t.0, &lease).unwrap(), selected);
    }
    #[test]
    fn adoption_refuses_stale_review_and_changed_recovery_and_startup_rechecks_source() {
        let t = Temp::new();
        let (mut repo, lab) = setup(&t);
        let mut service = CatalogBackupService::open(&t.0, &lab).unwrap();
        let backup = service.create(&op(310)).unwrap();
        service
            .restore_copy(&op(311), &backup.id, &backup.archive_hash)
            .unwrap();
        let stale = service
            .review_recovery(&op(311), &backup.archive_hash)
            .unwrap();
        repo.mutate_local(MutationOperation::CreateFolder {
            parent: None,
            name: "After review".into(),
        })
        .unwrap();
        assert!(service.adopt_recovery(&op(312), &stale).is_err());
        assert!(!catalog_activation_pending(&t.0).unwrap());
        let fresh = service
            .review_recovery(&op(311), &backup.archive_hash)
            .unwrap();
        service.adopt_recovery(&op(312), &fresh).unwrap();
        // Simulate an unsupported old writer after acceptance. Startup must refuse rather
        // than silently abandoning a newly committed source branch.
        repo.mutate_local(MutationOperation::CreateFolder {
            parent: None,
            name: "Late old writer".into(),
        })
        .unwrap();
        drop(repo);
        drop(service);
        let lease = CatalogLease::acquire(&t.0).unwrap();
        assert!(startup_catalog_path(&t.0, &lease).is_err());
        assert!(catalog_activation_pending(&t.0).unwrap());
        assert!(t.0.join("repository.sqlite").is_file());
    }
    #[test]
    fn current_catalog_frames_reject_every_noncurrent_header() {
        use collagene_desktop_protocol as w;
        let q = w::RequestEnvelope {
            request_id: 991,
            operation_id: op(320),
            command: w::Command::CatalogBackup(w::CatalogBackupCommand::ActivationStatus),
        };
        let request = w::encode_request(&q).unwrap();
        assert_eq!(
            u16::from_le_bytes([request[4], request[5]]),
            w::SCHEMA_VERSION
        );
        let reply = w::ResponseEnvelope {
            request_id: q.request_id,
            operation_id: q.operation_id.clone(),
            outcome: w::Outcome::Success(w::Success::CatalogBackup(
                w::CatalogBackupSuccess::ActivationStatus(w::CatalogActivationStatus::None),
            )),
        };
        let response = w::encode_response_for_request(&reply, &q).unwrap();
        w::decode_response_for(&response, &q).unwrap();
        for unsupported in [0, 1, w::SCHEMA_VERSION - 1, w::SCHEMA_VERSION + 1] {
            let mut bytes = request.clone();
            bytes[4..6].copy_from_slice(&unsupported.to_le_bytes());
            assert!(w::decode_request(&bytes).is_err());
            let mut bytes = response.clone();
            bytes[4..6].copy_from_slice(&unsupported.to_le_bytes());
            assert!(w::decode_response_for(&bytes, &q).is_err());
        }
    }

    #[test]
    fn adopted_selection_never_falls_back_after_lost_intent_and_process_lease_is_exclusive() {
        let t = Temp::new();
        let (_repo, lab) = setup(&t);
        let lease = CatalogLease::acquire(&t.0).unwrap();
        assert!(CatalogLease::acquire(&t.0).is_err());
        drop(lease);
        let mut service = CatalogBackupService::open(&t.0, &lab).unwrap();
        let backup = service.create(&op(330)).unwrap();
        service
            .restore_copy(&op(331), &backup.id, &backup.archive_hash)
            .unwrap();
        let review = service
            .review_recovery(&op(331), &backup.archive_hash)
            .unwrap();
        let mut forged = review.clone();
        forged.current_lab = base::LabId::new("lab-11111111111111111111111111111111").unwrap();
        forged.review_revision = forged.revision().unwrap();
        assert!(service.adopt_recovery(&op(332), &forged).is_err());
        service.adopt_recovery(&op(332), &review).unwrap();
        drop(service);
        let lease = CatalogLease::acquire(&t.0).unwrap();
        startup_catalog_path(&t.0, &lease).unwrap();
        fs::remove_file(t.0.join("catalog-activation/00000000000000000001.accepted.json")).unwrap();
        assert!(
            startup_catalog_path(&t.0, &lease).is_err(),
            "orphan applied receipt must not silently select original source"
        );
        assert!(t.0.join("repository.sqlite").is_file());
    }
    #[test]
    fn escaped_repository_and_backup_connections_retain_the_catalog_session_lease() {
        let t = Temp::new();
        let host = crate::DesktopHost::open(&t.0).unwrap();
        let mut repository = host.open_repository_connection().unwrap();
        let backups = host.open_catalog_backup().unwrap();
        // Both secondary handles share the original lock; construction must not
        // try acquiring a second exclusive OS lease during ordinary startup.
        drop(host);
        assert!(crate::DesktopHost::open(&t.0).is_err());
        let readonly = crate::open_local_repository(&t.0).unwrap();
        assert_eq!(
            readonly.local_context().unwrap(),
            repository.local_context().unwrap()
        );
        drop(readonly);
        repository
            .mutate_local(MutationOperation::CreateFolder {
                parent: None,
                name: "Escaped live connection".into(),
            })
            .unwrap();
        drop(repository);
        assert!(
            crate::DesktopHost::open(&t.0).is_err(),
            "backup connection still owns lease"
        );
        drop(backups);
        let reopened = crate::DesktopHost::open(&t.0).unwrap();
        let second = reopened.open_repository_connection().unwrap();
        drop(reopened);
        assert!(
            crate::DesktopHost::open(&t.0).is_err(),
            "repository alone still owns lease"
        );
        drop(second);
        assert!(crate::DesktopHost::open(&t.0).is_ok());
    }

    #[test]
    fn memory_only_catalog_reads_coexist_with_desktop_but_cannot_mutate_or_initialize() {
        let t = Temp::new();
        let mut host = crate::DesktopHost::open(&t.0).unwrap();
        let (payload, _) = native();
        let Applied::Imported(doc) = host
            .catalog
            .as_mut()
            .unwrap()
            .repository
            .mutate_local(MutationOperation::ImportDocument {
                folder: None,
                name: "MCP readable saved PCR".into(),
                payload,
            })
            .unwrap()
        else {
            panic!("import")
        };
        let mut readonly = crate::open_local_repository(&t.0).unwrap();
        assert_eq!(
            readonly
                .open_document(&doc.id, Some(&doc.head.id))
                .unwrap()
                .bytes,
            host.catalog
                .as_mut()
                .unwrap()
                .repository
                .open_document(&doc.id, Some(&doc.head.id))
                .unwrap()
                .bytes
        );
        assert!(
            readonly
                .mutate_local(MutationOperation::CreateFolder {
                    parent: None,
                    name: "Forbidden MCP catalog mutation".into()
                })
                .is_err()
        );
        assert_eq!(
            readonly.local_context().unwrap(),
            host.catalog
                .as_mut()
                .unwrap()
                .repository
                .local_context()
                .unwrap()
        );
        let missing = Temp::new();
        assert!(crate::open_local_repository(&missing.0).is_err());
        assert!(
            !missing.0.join("repository.sqlite").exists(),
            "read-only factory cannot create a database"
        );
    }
    #[test]
    fn repeated_adoption_and_rollback_keep_each_catalog_branch_and_generation() {
        let t = Temp::new();
        let (mut repo, lab) = setup(&t);
        let mut service = CatalogBackupService::open(&t.0, &lab).unwrap();
        let backup = service.create(&op(350)).unwrap();
        repo.mutate_local(MutationOperation::CreateFolder {
            parent: None,
            name: "Original branch preserved".into(),
        })
        .unwrap();
        service
            .restore_copy(&op(351), &backup.id, &backup.archive_hash)
            .unwrap();
        let first = service
            .review_recovery(&op(351), &backup.archive_hash)
            .unwrap();
        assert_eq!(first.activation_generation, 0);
        let receipt = service.adopt_recovery(&op(352), &first).unwrap();
        drop(service);
        drop(repo);
        let lease = CatalogLease::acquire(&t.0).unwrap();
        let recovered = startup_catalog_path(&t.0, &lease).unwrap();
        let mut branch = collagene_repository::Repository::open(&recovered).unwrap();
        branch
            .mutate_local(MutationOperation::CreateFolder {
                parent: None,
                name: "Recovery private branch preserved".into(),
            })
            .unwrap();
        drop(branch);
        drop(lease);
        let service = CatalogBackupService::open(&t.0, &lab).unwrap();
        let rollback = service.review_original().unwrap();
        assert_eq!(rollback.activation_generation, 1);
        assert_eq!(rollback.target_snapshot.folders, 1);
        service.adopt_recovery(&op(353), &rollback).unwrap();
        drop(service);
        let lease = CatalogLease::acquire(&t.0).unwrap();
        assert_eq!(
            startup_catalog_path(&t.0, &lease).unwrap(),
            t.0.join("repository.sqlite")
        );
        drop(lease);
        let service = CatalogBackupService::open(&t.0, &lab).unwrap();
        assert_eq!(
            service.adopt_recovery(&op(352), &first).unwrap(),
            receipt,
            "historical exact retry is not a new transition"
        );
        let third = service
            .review_recovery(&op(351), &backup.archive_hash)
            .unwrap();
        assert_eq!(third.activation_generation, 2);
        assert_eq!(
            third.target_snapshot.folders, 1,
            "previous recovery writes were retained"
        );
        assert!(
            service.adopt_recovery(&op(354), &first).is_err(),
            "old generation cannot adopt a new operation"
        );
        service.adopt_recovery(&op(354), &third).unwrap();
        drop(service);
        let lease = CatalogLease::acquire(&t.0).unwrap();
        assert_eq!(startup_catalog_path(&t.0, &lease).unwrap(), recovered);
        assert!(t.0.join("repository.sqlite").is_file());
    }
}

#[cfg(test)]
mod raw_tests {
    use super::*;
    use collagene_desktop_protocol as wire;
    #[test]
    fn scoped_native_frames_cache_exact_receipts_reverify_artifacts_and_cancel_before_acceptance() {
        let root = std::env::temp_dir().join(format!(
            "collagene-backup-raw-{}-{}",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        struct Cleanup(PathBuf);
        impl Drop for Cleanup {
            fn drop(&mut self) {
                let _ = fs::remove_dir_all(&self.0);
            }
        }
        let _cleanup = Cleanup(root.clone());
        let mut host = crate::DesktopHost::open(&root).unwrap();
        let access = host.authorize_window().unwrap();
        let backups = host.open_catalog_backup().unwrap();
        let op = |n: u64| OperationId::new(format!("op-{n:032x}")).unwrap();
        let create = wire::RequestEnvelope {
            request_id: 1,
            operation_id: op(1),
            command: wire::Command::CatalogBackup(wire::CatalogBackupCommand::Create),
        };
        let bytes = wire::encode_request(&create).unwrap();
        let reply = backups.dispatch_raw(&access, &bytes).unwrap();
        assert_eq!(backups.dispatch_raw(&access, &bytes).unwrap(), reply);
        assert!(reply.len() < wire::MAX_CATALOG_BACKUP_RESPONSE_BYTES);
        let response = wire::decode_response_for(&reply, &create).unwrap();
        let wire::Outcome::Success(wire::Success::CatalogBackup(
            wire::CatalogBackupSuccess::Created(receipt),
        )) = response.outcome
        else {
            panic!("created backup")
        };
        assert_eq!(receipt.id, op(1));
        let verify = wire::RequestEnvelope {
            request_id: 2,
            operation_id: op(2),
            command: wire::Command::CatalogBackup(wire::CatalogBackupCommand::Verify {
                id: op(1),
                expected_hash: Some(base::Digest::of(b"different explicit selected hash")),
            }),
        };
        let reply = backups
            .dispatch_raw(&access, &wire::encode_request(&verify).unwrap())
            .unwrap();
        assert!(matches!(
            wire::decode_response_for(&reply, &verify).unwrap().outcome,
            wire::Outcome::Failure(_)
        ));
        let cancelled = op(3);
        let cancel = wire::RequestEnvelope {
            request_id: 4,
            operation_id: op(4),
            command: wire::Command::Cancel {
                operation_id: cancelled.clone(),
            },
        };
        let control = host.control();
        assert!(
            control
                .dispatch_raw(&access, &wire::encode_request(&cancel).unwrap())
                .unwrap()
                .is_some()
        );
        let attempt = wire::RequestEnvelope {
            request_id: 3,
            operation_id: cancelled.clone(),
            command: wire::Command::CatalogBackup(wire::CatalogBackupCommand::Create),
        };
        let reply = backups
            .dispatch_raw(&access, &wire::encode_request(&attempt).unwrap())
            .unwrap();
        assert!(matches!(
            wire::decode_response_for(&reply, &attempt).unwrap().outcome,
            wire::Outcome::Failure(wire::Failure::Cancelled)
        ));
        assert!(
            !root
                .join("backups")
                .join(format!("{}.cglbak", cancelled.as_str()))
                .exists()
        );
        // An in-memory cached ACK never bypasses actual artifact verification.
        let artifact = root.join("backups").join(receipt.basename.as_str());
        let mut file = fs::read(&artifact).unwrap();
        let last = file.len() - 1;
        file[last] ^= 1;
        fs::write(&artifact, file).unwrap();
        let reply = backups.dispatch_raw(&access, &bytes).unwrap();
        assert!(matches!(
            wire::decode_response_for(&reply, &create).unwrap().outcome,
            wire::Outcome::Failure(_)
        ));
        assert!(
            backups
                .dispatch_raw(
                    &access,
                    &vec![0; wire::MAX_CATALOG_BACKUP_REQUEST_BYTES + 1]
                )
                .is_err()
        );
        let wrong = wire::RequestEnvelope {
            request_id: 5,
            operation_id: op(5),
            command: wire::Command::Bootstrap,
        };
        assert!(
            backups
                .dispatch_raw(&access, &wire::encode_request(&wrong).unwrap())
                .is_err()
        );
        let other_root = root.join("other-host");
        let mut other = crate::DesktopHost::open(other_root).unwrap();
        let foreign = other.authorize_window().unwrap();
        assert!(matches!(
            backups.dispatch_raw(&foreign, &bytes),
            Err(HostError::Window)
        ));
    }
}
