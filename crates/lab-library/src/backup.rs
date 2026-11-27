//! Native, whole-database snapshots. This is an operator backup, never a public
//! publication or a source of scientific/physical authority.
use super::*;
use rusqlite::{
    OpenFlags,
    backup::{Backup, StepResult},
};
use std::{
    collections::BTreeMap,
    fs::{self, File, OpenOptions},
    io::{Read, Write},
    path::{Path, PathBuf},
    sync::atomic::{AtomicU64, Ordering},
    time::{Instant, SystemTime, UNIX_EPOCH},
};
const MAGIC: &[u8; 8] = b"CGLBAK\0\x01";
const MAX_MANIFEST: usize = 64 * 1024;
const HEADER_BYTES: usize = 44;
pub const DEFAULT_MAX_DATABASE_BYTES: u64 = 8 * 1024 * 1024 * 1024;
pub const MAX_DATABASE_BYTES: u64 = 1 << 46;
#[derive(Debug, thiserror::Error)]
pub enum BackupError {
    #[error("invalid or unsupported backup: {0}")]
    Invalid(&'static str),
    #[error("backup storage I/O failed: {0}")]
    Io(#[from] std::io::Error),
    #[error("backup SQLite operation failed: {0}")]
    Sqlite(#[from] rusqlite::Error),
    #[error("library integrity failed: {0}")]
    Storage(#[from] super::Error),
    #[error("reliable storage integrity failed: {0}")]
    Reliable(#[from] mutations::MutationError),
}
pub type Result<T> = std::result::Result<T, BackupError>;
#[derive(Debug, Clone)]
pub struct Policy {
    pub expected_labs: Vec<LabId>,
    pub max_database_bytes: u64,
}
impl Policy {
    pub fn single_lab(lab: LabId) -> Self {
        Self {
            expected_labs: vec![lab],
            max_database_bytes: DEFAULT_MAX_DATABASE_BYTES,
        }
    }
    fn validate(&self) -> Result<()> {
        if self.expected_labs.is_empty()
            || self.expected_labs.len() > 128
            || self.max_database_bytes == 0
            || self.max_database_bytes > MAX_DATABASE_BYTES
        {
            return Err(invalid("policy limits"));
        }
        let mut ids = self.expected_labs.clone();
        ids.sort();
        ids.dedup();
        if ids.len() != self.expected_labs.len() {
            return Err(invalid("duplicate expected lab"));
        }
        Ok(())
    }
}
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct LabPin {
    pub id: LabId,
    pub catalog_generation: u64,
}
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Manifest {
    pub format: String,
    pub version: u32,
    pub application_id: i64,
    pub schema_version: i64,
    pub schema_hash: Digest,
    pub database_hash: Digest,
    pub database_bytes: u64,
    pub created_unix_seconds: u64,
    pub max_lab_archive_bytes: u64,
    pub labs: Vec<LabPin>,
    pub table_rows: BTreeMap<String, u64>,
}
fn invalid(s: &'static str) -> BackupError {
    BackupError::Invalid(s)
}
fn digest(bytes: &[u8]) -> Digest {
    Digest(blake3::hash(bytes).to_hex().to_string())
}
fn private_dir(path: &Path) -> Result<()> {
    let m = fs::symlink_metadata(path)?;
    if !m.is_dir() || m.file_type().is_symlink() {
        return Err(invalid("directory must be real, not a symlink"));
    }
    #[cfg(unix)]
    {
        use std::os::unix::fs::MetadataExt;
        if m.uid() != unsafe { libc::geteuid() } || m.mode() & 0o077 != 0 {
            return Err(invalid("directory must be owner-only"));
        }
    }
    Ok(())
}
fn private_file(path: &Path) -> Result<()> {
    let m = fs::symlink_metadata(path)?;
    if !m.is_file() || m.file_type().is_symlink() {
        return Err(invalid("file must be regular, not a symlink"));
    }
    #[cfg(unix)]
    {
        use std::os::unix::fs::MetadataExt;
        if m.uid() != unsafe { libc::geteuid() } || m.mode() & 0o077 != 0 {
            return Err(invalid("file must be owner-only"));
        }
    }
    Ok(())
}
fn mkdir(path: &Path) -> Result<()> {
    let mut d = fs::DirBuilder::new();
    #[cfg(unix)]
    {
        use std::os::unix::fs::DirBuilderExt;
        d.mode(0o700);
    }
    d.create(path)?;
    Ok(())
}
fn new_file(path: &Path) -> Result<File> {
    let mut o = OpenOptions::new();
    o.write(true).read(true).create_new(true);
    #[cfg(unix)]
    {
        use std::os::unix::fs::OpenOptionsExt;
        o.mode(0o600);
    }
    Ok(o.open(path)?)
}
fn sync_dir(path: &Path) -> Result<()> {
    #[cfg(unix)]
    File::open(path)?.sync_all()?;
    let _ = path;
    Ok(())
}
fn target(path: &Path) -> Result<PathBuf> {
    let parent = path
        .parent()
        .ok_or_else(|| invalid("missing parent directory"))?;
    private_dir(parent)?;
    let name = path
        .file_name()
        .ok_or_else(|| invalid("missing filename"))?;
    let p = parent.canonicalize()?.join(name);
    if fs::symlink_metadata(&p).is_ok() {
        return Err(invalid("destination already exists"));
    }
    Ok(p)
}
struct Stage(PathBuf);
impl Stage {
    fn new(parent: &Path) -> Result<Self> {
        static NEXT: AtomicU64 = AtomicU64::new(0);
        private_dir(parent)?;
        for _ in 0..100 {
            let p = parent.join(format!(
                ".collagene-backup-{}-{}",
                std::process::id(),
                NEXT.fetch_add(1, Ordering::Relaxed)
            ));
            match mkdir(&p) {
                Ok(()) => return Ok(Self(p)),
                Err(BackupError::Io(e)) if e.kind() == std::io::ErrorKind::AlreadyExists => (),
                Err(e) => return Err(e),
            }
        }
        Err(invalid("cannot allocate private staging directory"))
    }
}
impl Drop for Stage {
    fn drop(&mut self) {
        let _ = fs::remove_dir_all(&self.0);
    }
}
fn schema(c: &Connection) -> Result<Vec<(String, String, String, Option<String>)>> {
    let mut s=c.prepare("SELECT type,name,tbl_name,sql FROM sqlite_master WHERE name NOT LIKE 'sqlite_%' ORDER BY type,name")?;
    let rows = s.query_map([], |r| Ok((r.get(0)?, r.get(1)?, r.get(2)?, r.get(3)?)))?;
    let mut out = vec![];
    for row in rows {
        if out.len() >= 128 {
            return Err(invalid("schema object count"));
        }
        out.push(row?)
    }
    Ok(out)
}
fn schema_check(c: &Connection) -> Result<Digest> {
    validate_database_header(
        c.query_row("PRAGMA application_id", [], |r| r.get::<_, i64>(0))?,
        c.query_row("PRAGMA user_version", [], |r| r.get::<_, i64>(0))?,
    )
    .map_err(|_| invalid("foreign or unsupported database"))?;
    let actual = schema(c)?;
    let mut expected = Library::open_in_memory()?;
    if actual.iter().any(|(_, n, _, _)| n == "cg_reliable_schema") {
        expected.local_repository_access("backup-schema-probe", "Schema probe")?;
        let v: i64 = c.query_row(
            "SELECT version FROM cg_reliable_schema WHERE singleton=1",
            [],
            |r| r.get(0),
        )?;
        if v != 1
            || c.query_row("SELECT count(*) FROM cg_reliable_schema", [], |r| {
                r.get::<_, i64>(0)
            })? != 1
        {
            return Err(invalid("reliable schema version"));
        }
    }
    if actual != schema(&expected.connection)? {
        return Err(invalid("unknown or modified schema"));
    }
    Ok(digest(
        &serde_json::to_vec(&actual).map_err(|_| invalid("schema encoding"))?,
    ))
}
fn inspect_database(
    path: &Path,
    policy: &Policy,
    options: StorageOptions,
) -> Result<(Digest, Vec<LabPin>, BTreeMap<String, u64>)> {
    private_file(path)?;
    if fs::metadata(path)?.len() > policy.max_database_bytes {
        return Err(invalid("database size bound"));
    }
    let c = Connection::open_with_flags(
        path,
        OpenFlags::SQLITE_OPEN_READ_ONLY | OpenFlags::SQLITE_OPEN_NO_MUTEX,
    )?;
    c.execute_batch("PRAGMA trusted_schema=OFF; PRAGMA foreign_keys=ON;")?;
    let schema_hash = schema_check(&c)?;
    if c.query_row("PRAGMA integrity_check", [], |r| r.get::<_, String>(0))? != "ok" {
        return Err(invalid("SQLite integrity check"));
    }
    if c.prepare("PRAGMA foreign_key_check")?
        .query([])?
        .next()?
        .is_some()
    {
        return Err(invalid("foreign key violation"));
    }
    let mut labs = vec![];
    let mut st = c.prepare("SELECT id,catalog_generation,archive_bytes FROM labs ORDER BY id")?;
    for row in st.query_map([], |r| {
        Ok((
            r.get::<_, String>(0)?,
            r.get::<_, u64>(1)?,
            r.get::<_, u64>(2)?,
        ))
    })? {
        let (id, generation, charged) = row?;
        if labs.len() >= 128 {
            return Err(invalid("lab count"));
        }
        let id = LabId::parse(id)?;
        let actual:u64=c.query_row("SELECT coalesce(sum(r.archive_size),0) FROM revisions r JOIN documents d ON d.id=r.document_id WHERE d.lab_id=?1",[id.as_str()],|r|r.get(0))?;
        if actual != charged {
            return Err(invalid("lab archive accounting"));
        }
        labs.push(LabPin {
            id,
            catalog_generation: generation,
        });
    }
    let mut expected = policy.expected_labs.clone();
    expected.sort();
    if labs.iter().map(|l| &l.id).collect::<Vec<_>>() != expected.iter().collect::<Vec<_>>() {
        return Err(invalid("expected complete lab-ID set differs"));
    }
    let mut counts = BTreeMap::new();
    for (ty, name, _, _) in schema(&c)? {
        if ty == "table" {
            let n: u64 = c.query_row(&format!("SELECT count(*) FROM \"{name}\""), [], |r| {
                r.get(0)
            })?;
            if n > 12_800_000 {
                return Err(invalid("table row bound"));
            }
            counts.insert(name, n);
        }
    }
    // All archives are SQLite BLOBs. Check one bounded archive at a time and bind
    // each immutable revision to its exact document, parent and payload metadata.
    drop(st);
    let library = Library {
        connection: c,
        options,
        native_owner: None,
    };
    let mut st=library.connection.prepare("SELECT r.id,r.document_id,d.lab_id,r.number,r.parent_id FROM revisions r JOIN documents d ON d.id=r.document_id ORDER BY r.document_id,r.number")?;
    let mut previous: Option<(String, String, u64)> = None;
    for row in st.query_map([], |r| {
        Ok((
            r.get::<_, String>(0)?,
            r.get::<_, String>(1)?,
            r.get::<_, String>(2)?,
            r.get::<_, u64>(3)?,
            r.get::<_, Option<String>>(4)?,
        ))
    })? {
        let (id, doc, lab, number, parent) = row?;
        if let Some((pdoc, pid, pnum)) = &previous {
            if pdoc == &doc {
                if number != pnum + 1 || parent.as_ref() != Some(pid) {
                    return Err(invalid("revision chain"));
                }
            } else if number != 1 || parent.is_some() {
                return Err(invalid("revision root"));
            }
        } else if number != 1 || parent.is_some() {
            return Err(invalid("revision root"));
        }
        let access = library.trusted_host_access(&LabId::parse(lab)?)?;
        library.payload(
            &access,
            &DocumentId::parse(doc.clone())?,
            &RevisionId::parse(id.clone())?,
        )?;
        previous = Some((doc, id, number));
    }
    collaboration::check_integrity(&library.connection)?;
    if counts.contains_key("cg_reliable_schema") {
        check_reliable(&library.connection, &labs)?;
    }
    // Detect cycles/depth violations and malformed current metadata without changing it.
    for lab in &labs {
        let access = library.trusted_host_access(&lab.id)?;
        let mut q = library
            .connection
            .prepare("SELECT id FROM folders WHERE lab_id=?1")?;
        for row in q.query_map([lab.id.as_str()], |r| r.get::<_, String>(0))? {
            let id = FolderId::parse(row?)?;
            let folder = library.folder(&access, &id)?;
            let mut seen = std::collections::BTreeSet::new();
            seen.insert(id);
            let mut parent = folder.parent;
            while let Some(p) = parent {
                if seen.len() >= MAX_FOLDER_DEPTH || !seen.insert(p.clone()) {
                    return Err(invalid("folder ancestry"));
                }
                parent = library.folder(&access, &p)?.parent;
            }
        }
        let mut q = library
            .connection
            .prepare("SELECT id FROM documents WHERE lab_id=?1")?;
        for row in q.query_map([lab.id.as_str()], |r| r.get::<_, String>(0))? {
            let id = DocumentId::parse(row?)?;
            let doc = get_document(&library.connection, &lab.id, &id)?;
            let latest: u64 = library.connection.query_row(
                "SELECT max(number) FROM revisions WHERE document_id=?1",
                [id.as_str()],
                |r| r.get(0),
            )?;
            if doc.head.number != latest {
                return Err(invalid("document head is not latest immutable revision"));
            }
        }
    }
    Ok((schema_hash, labs, counts))
}
fn check_reliable(c: &Connection, labs: &[LabPin]) -> Result<()> {
    for lab in labs {
        let mut receipts = 0u64;
        let mut receipt_bytes = 0u64;
        let mut q=c.prepare("SELECT id,intent,CASE WHEN length(receipt)<=131072 THEN receipt ELSE NULL END,receipt_hash FROM cg_reliable_mutations WHERE lab_id=?1")?;
        for row in q.query_map([lab.id.as_str()], |r| {
            Ok((
                r.get::<_, String>(0)?,
                r.get::<_, String>(1)?,
                r.get::<_, Option<Vec<u8>>>(2)?,
                r.get::<_, String>(3)?,
            ))
        })? {
            let (id, intent, bytes, hash) = row?;
            let bytes = bytes.ok_or_else(|| invalid("receipt bound"))?;
            if digest(&bytes).as_str() != hash {
                return Err(invalid("receipt hash"));
            }
            let receipt: mutations::Receipt =
                serde_json::from_slice(&bytes).map_err(|_| invalid("receipt data"))?;
            if receipt.id.as_str() != id
                || receipt.intent.as_str() != intent
                || receipt.catalog_generation > lab.catalog_generation
            {
                return Err(invalid("receipt binding"));
            }
            receipts += 1;
            receipt_bytes += bytes.len() as u64;
        }
        let mut checkpoint_bytes = 0u64;
        let mut checkpoints = 0u64;
        let mut q=c.prepare("SELECT namespace,key,version,CASE WHEN length(value)<=67108864 THEN value ELSE NULL END,digest FROM cg_reliable_checkpoints WHERE lab_id=?1")?;
        for row in q.query_map([lab.id.as_str()], |r| {
            Ok((
                r.get::<_, String>(0)?,
                r.get::<_, String>(1)?,
                r.get::<_, u32>(2)?,
                r.get::<_, Option<Vec<u8>>>(3)?,
                r.get::<_, String>(4)?,
            ))
        })? {
            let (ns, key, version, bytes, stored) = row?;
            let bytes = bytes.ok_or_else(|| invalid("checkpoint bound"))?;
            if ns != "repository"
                || key.is_empty()
                || key.len() > 128
                || !key
                    .bytes()
                    .all(|b| b.is_ascii_alphanumeric() || matches!(b, b'-' | b'_'))
                || version == 0
            {
                return Err(invalid("checkpoint key/version"));
            }
            let mut h = blake3::Hasher::new();
            h.update(b"collagene-private-checkpoint-v1\0");
            for part in [
                lab.id.as_str().as_bytes(),
                ns.as_bytes(),
                key.as_bytes(),
                &version.to_le_bytes(),
                &bytes,
            ] {
                h.update(&(part.len() as u64).to_le_bytes());
                h.update(part);
            }
            if h.finalize().to_hex().as_str() != stored {
                return Err(invalid("checkpoint hash"));
            }
            checkpoint_bytes += bytes.len() as u64;
            checkpoints += 1;
        }
        if receipts > 100_000
            || receipt_bytes > 64 * 1024 * 1024
            || checkpoint_bytes > 128 * 1024 * 1024
            || checkpoints > 4096
        {
            return Err(invalid("reliable quota"));
        }
        let usage:Option<(u64,u64,u64)>=c.query_row("SELECT receipt_count,receipt_bytes,checkpoint_bytes FROM cg_reliable_usage WHERE lab_id=?1",[lab.id.as_str()],|r|Ok((r.get(0)?,r.get(1)?,r.get(2)?))).optional()?;
        if usage.unwrap_or((0, 0, 0)) != (receipts, receipt_bytes, checkpoint_bytes) {
            return Err(invalid("reliable accounting"));
        }
    }
    Ok(())
}
fn stream_hash(file: &mut File, bytes: u64, mut output: Option<&mut File>) -> Result<Digest> {
    let mut h = blake3::Hasher::new();
    let mut left = bytes;
    let mut b = [0u8; 64 * 1024];
    while left > 0 {
        let n = (left.min(b.len() as u64)) as usize;
        file.read_exact(&mut b[..n])?;
        h.update(&b[..n]);
        if let Some(out) = output.as_deref_mut() {
            out.write_all(&b[..n])?;
        }
        left -= n as u64;
    }
    Ok(Digest(h.finalize().to_hex().to_string()))
}
fn encode_manifest(m: &Manifest) -> Result<Vec<u8>> {
    let b = serde_json::to_vec(m).map_err(|_| invalid("manifest encoding"))?;
    if b.len() > MAX_MANIFEST {
        return Err(invalid("manifest size"));
    }
    Ok(b)
}
/// Online consistent snapshot, including live WAL contents. The final file is
/// published atomically without overwriting any existing directory entry.
pub fn create(library: &Library, backup_path: &Path, policy: &Policy) -> Result<Manifest> {
    policy.validate()?;
    let final_path = target(backup_path)?;
    let parent = final_path
        .parent()
        .ok_or_else(|| invalid("backup parent"))?;
    let stage = Stage::new(parent)?;
    let db = stage.0.join("database.sqlite");
    new_file(&db)?;
    {
        let mut destination = Connection::open(&db)?;
        let page_size: u64 = library
            .connection
            .query_row("PRAGMA page_size", [], |r| r.get(0))?;
        if !(512..=65536).contains(&page_size) || !page_size.is_power_of_two() {
            return Err(invalid("SQLite page size"));
        }
        let backup = Backup::new(&library.connection, &mut destination)?;
        let deadline = Instant::now() + Duration::from_secs(120);
        loop {
            if Instant::now() >= deadline {
                return Err(invalid("online snapshot deadline"));
            }
            let step = backup.step(256)?;
            let progress = backup.progress();
            if progress.pagecount > 0
                && (progress.pagecount as u64) * page_size > policy.max_database_bytes
            {
                return Err(invalid("snapshot page bound"));
            }
            if fs::metadata(&db)?.len() > policy.max_database_bytes {
                return Err(invalid("snapshot byte bound"));
            }
            match step {
                StepResult::Done => break,
                StepResult::Busy | StepResult::Locked => {
                    std::thread::sleep(Duration::from_millis(5))
                }
                _ => (),
            }
        }
        drop(backup);
        destination.execute_batch("PRAGMA journal_mode=DELETE; PRAGMA synchronous=FULL;")?;
    }
    OpenOptions::new()
        .read(true)
        .write(true)
        .open(&db)?
        .sync_all()?;
    let (schema_hash, labs, table_rows) = inspect_database(&db, policy, library.options)?;
    let size = fs::metadata(&db)?.len();
    let hash = stream_hash(&mut File::open(&db)?, size, None)?;
    let manifest = Manifest {
        format: "collagene-lab-backup".into(),
        version: 1,
        application_id: APPLICATION_ID,
        schema_version: SCHEMA_VERSION,
        schema_hash,
        database_hash: hash,
        database_bytes: size,
        created_unix_seconds: SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .map_err(|_| invalid("system clock"))?
            .as_secs(),
        max_lab_archive_bytes: library.options.max_lab_archive_bytes,
        labs,
        table_rows,
    };
    let bytes = encode_manifest(&manifest)?;
    let packed = stage.0.join("backup.pending");
    let mut out = new_file(&packed)?;
    out.write_all(MAGIC)?;
    out.write_all(&(bytes.len() as u32).to_le_bytes())?;
    out.write_all(blake3::hash(&bytes).as_bytes())?;
    out.write_all(&bytes)?;
    stream_hash(&mut File::open(&db)?, size, Some(&mut out))?;
    out.sync_all()?;
    drop(out);
    fs::hard_link(&packed, &final_path)?;
    sync_dir(parent)?;
    Ok(manifest)
}
fn extract(backup_path: &Path, policy: &Policy, stage: &Stage) -> Result<(Manifest, PathBuf)> {
    policy.validate()?;
    private_file(backup_path)?;
    let mut source = File::open(backup_path)?;
    let mut header = [0u8; HEADER_BYTES];
    source.read_exact(&mut header)?;
    if &header[..8] != MAGIC {
        return Err(invalid("magic/version"));
    }
    let n = u32::from_le_bytes(header[8..12].try_into().map_err(|_| invalid("header"))?) as usize;
    if n == 0 || n > MAX_MANIFEST {
        return Err(invalid("manifest size"));
    }
    let mut bytes = vec![0; n];
    source.read_exact(&mut bytes)?;
    if blake3::hash(&bytes).as_bytes() != &header[12..44] {
        return Err(invalid("manifest content hash"));
    }
    let manifest: Manifest =
        serde_json::from_slice(&bytes).map_err(|_| invalid("manifest shape"))?;
    if encode_manifest(&manifest)? != bytes
        || manifest.format != "collagene-lab-backup"
        || manifest.version != 1
        || manifest.application_id != APPLICATION_ID
        || manifest.schema_version != SCHEMA_VERSION
        || manifest.database_bytes == 0
        || manifest.database_bytes > policy.max_database_bytes
        || source.metadata()?.len() != HEADER_BYTES as u64 + n as u64 + manifest.database_bytes
    {
        return Err(invalid("manifest bounds/version/size"));
    }
    let options = StorageOptions {
        max_lab_archive_bytes: manifest.max_lab_archive_bytes,
    };
    options.validate()?;
    let db = stage.0.join("database.sqlite");
    let mut out = new_file(&db)?;
    let hash = stream_hash(&mut source, manifest.database_bytes, Some(&mut out))?;
    out.sync_all()?;
    drop(out);
    if hash != manifest.database_hash {
        return Err(invalid("database content hash"));
    }
    let (schema_hash, labs, counts) = inspect_database(&db, policy, options)?;
    if schema_hash != manifest.schema_hash || labs != manifest.labs || counts != manifest.table_rows
    {
        return Err(invalid("manifest inventory binding"));
    }
    Ok((manifest, db))
}
/// Read and verify all database bytes without modifying the source or live DB.
pub fn verify(backup_path: &Path, policy: &Policy) -> Result<Manifest> {
    let parent = backup_path
        .parent()
        .ok_or_else(|| invalid("backup parent"))?;
    let stage = Stage::new(parent)?;
    Ok(extract(backup_path, policy, &stage)?.0)
}
/// Restore exclusively to a previously absent directory. No live or existing
/// destination is opened, replaced or renamed. Keep the backup for rollback.
pub fn restore_new(
    backup_path: &Path,
    destination_root: &Path,
    policy: &Policy,
) -> Result<Manifest> {
    let destination = target(destination_root)?;
    let parent = destination
        .parent()
        .ok_or_else(|| invalid("restore parent"))?;
    let stage = Stage::new(parent)?;
    let (manifest, db) = extract(backup_path, policy, &stage)?;
    mkdir(&destination)?;
    let published = destination.join("library.sqlite");
    if let Err(error) = fs::hard_link(&db, &published) {
        let _ = fs::remove_dir(&destination);
        return Err(error.into());
    }
    sync_dir(&destination)?;
    sync_dir(parent)?;
    Ok(manifest)
}
