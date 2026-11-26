//! Append-only, restart-bound catalog selection. Never overwrite a database branch.
use super::*;
use collagene_desktop_protocol::{
    CatalogActivationReceipt as Activation, CatalogActivationStatus as Status,
    CatalogCancellationReceipt as Cancellation, CatalogLocation as Location,
    CatalogRecoveryReview as Review, CatalogSnapshot as Snapshot,
};
const MAX_ACTIVATIONS: u64 = 4096;
fn directory(root: &Path) -> PathBuf {
    root.join("catalog-activation")
}
fn record(root: &Path, g: u64, kind: &str) -> PathBuf {
    directory(root).join(format!("{g:020}.{kind}.json"))
}
struct Journal {
    generation: u64,
    current: Location,
    last: Option<Activation>,
    pending: Option<Activation>,
    cancelled: Option<Cancellation>,
    visited: Vec<Location>,
}
fn journal(root: &Path) -> HostResult<Journal> {
    let mut j = Journal {
        generation: 0,
        current: Location::Original,
        last: None,
        pending: None,
        cancelled: None,
        visited: vec![Location::Original],
    };
    let dir = directory(root);
    if !exists(&dir)? {
        return Ok(j);
    }
    private(&dir, true)?;
    let mut generations = std::collections::BTreeSet::new();
    let mut count = 0;
    for e in fs::read_dir(&dir).map_err(io)? {
        count += 1;
        if count > MAX_ACTIVATIONS * 3 + 1 {
            return Err(HostError::Corrupt);
        }
        let e = e.map_err(io)?;
        let name = e.file_name().into_string().map_err(io)?;
        if name == "staging" {
            private(&e.path(), true)?;
            continue;
        }
        private(&e.path(), false)?;
        let number = name
            .strip_suffix(".accepted.json")
            .or_else(|| name.strip_suffix(".applied.json"))
            .or_else(|| name.strip_suffix(".cancelled.json"))
            .ok_or(HostError::Corrupt)?;
        let n = number.parse::<u64>().map_err(io)?;
        if n == 0 || n > MAX_ACTIVATIONS || number != format!("{n:020}") {
            return Err(HostError::Corrupt);
        }
        generations.insert(n);
    }
    let mut operations = std::collections::BTreeSet::new();
    for n in generations {
        if n != j.generation + 1 || j.pending.is_some() {
            return Err(HostError::Corrupt);
        }
        let a: Activation = read(&record(root, n, "accepted"))?;
        a.review.validate().map_err(|_| HostError::Corrupt)?;
        if a.review.activation_generation != j.generation
            || a.review.current_catalog != j.current
            || !operations.insert(a.operation_id.clone())
        {
            return Err(HostError::Corrupt);
        }
        let p = record(root, n, "applied");
        let cancelled_path = record(root, n, "cancelled");
        if exists(&p)? && exists(&cancelled_path)? {
            return Err(HostError::Corrupt);
        }
        if exists(&cancelled_path)? {
            let cancelled: Cancellation = read(&cancelled_path)?;
            cancelled.validate().map_err(|_| HostError::Corrupt)?;
            if cancelled.accepted != a
                || cancelled.generation != n
                || !operations.insert(cancelled.operation_id.clone())
            {
                return Err(HostError::Corrupt);
            }
            j.generation = n;
            j.cancelled = Some(cancelled);
        } else if exists(&p)? {
            let saved: Activation = read(&p)?;
            if saved != a {
                return Err(HostError::Corrupt);
            }
            j.generation = n;
            j.current = a.review.target.clone();
            j.visited.push(j.current.clone());
            j.last = Some(a);
            j.cancelled = None;
        } else {
            j.pending = Some(a)
        }
    }
    Ok(j)
}
pub(super) fn activation_operation_in_use(root: &Path, op: &str) -> HostResult<bool> {
    let j = journal(root)?;
    for n in 1..=j.generation + u64::from(j.pending.is_some()) {
        let a: Activation = read(&record(root, n, "accepted"))?;
        if a.operation_id.as_str() == op {
            return Ok(true);
        }
        let p = record(root, n, "cancelled");
        if exists(&p)? {
            let c: Cancellation = read(&p)?;
            if c.operation_id.as_str() == op {
                return Ok(true);
            }
        }
    }
    Ok(false)
}
fn location_path(root: &Path, location: &Location) -> HostResult<PathBuf> {
    match location {
        Location::Original => {
            let p = root.join("repository.sqlite");
            private(&p, false)?;
            permissions::database(root, "repository.sqlite")
        }
        Location::Recovery { id } => {
            let parent = root.join("recoveries");
            private(&parent, true)?;
            let dir = parent.join(id.as_str());
            private(&dir, true)?;
            private(&dir.join("library.sqlite"), false)?;
            permissions::database(&dir, "library.sqlite")
        }
    }
}
fn scope(path: &Path, lab: &base::LabId) -> HostResult<()> {
    private(path, false)?;
    let c = rusqlite::Connection::open_with_flags(path, rusqlite::OpenFlags::SQLITE_OPEN_READ_ONLY)
        .map_err(io)?;
    c.execute_batch("PRAGMA trusted_schema=OFF; PRAGMA query_only=ON")
        .map_err(io)?;
    let app: i64 = c
        .query_row("PRAGMA application_id", [], |r| r.get(0))
        .map_err(io)?;
    let schema: i64 = c
        .query_row("PRAGMA user_version", [], |r| r.get(0))
        .map_err(io)?;
    let mut q = c
        .prepare("SELECT id FROM labs ORDER BY id LIMIT 2")
        .map_err(io)?;
    let labs = q
        .query_map([], |r| r.get::<_, String>(0))
        .map_err(io)?
        .collect::<Result<Vec<_>, _>>()
        .map_err(io)?;
    if app != 0x43474c42 || schema != 1 || labs != [lab.as_str()] {
        return Err(HostError::Corrupt);
    }
    Ok(())
}
pub fn catalog_activation_pending(root: impl AsRef<Path>) -> HostResult<bool> {
    let root = permissions::root(root.as_ref())?;
    Ok(journal(&root)?.pending.is_some())
}
pub(crate) fn current_catalog_path(root: &Path) -> HostResult<PathBuf> {
    let j = journal(root)?;
    let p = location_path(root, &j.current)?;
    if let Some(cancelled) = j.cancelled {
        scope(&p, &cancelled.accepted.review.current_lab)?;
    } else if let Some(last) = j.last {
        scope(&p, &last.review.current_lab)?;
    }
    Ok(p)
}
pub(crate) fn check_catalog_connection(root: &Path, path: &Path) -> HostResult<()> {
    let j = journal(root)?;
    if j.pending.is_some() || location_path(root, &j.current)? != path {
        return Err(invalid(
            "Catalog adoption was accepted. Restart before using the saved catalog; editor recovery remains available.",
        ));
    }
    Ok(())
}
pub(crate) fn startup_catalog_path(
    root: &Path,
    lease: &std::sync::Arc<CatalogLease>,
) -> HostResult<PathBuf> {
    let j = journal(root)?;
    if j.generation == 0 && j.pending.is_none() {
        return permissions::database(root, "repository.sqlite");
    }
    if let Some(a) = j.pending {
        let lab = store::LabId::parse(a.review.current_lab.as_str()).map_err(io)?;
        let service = CatalogBackupService::open_with_lease(root, &lab, lease.clone())?;
        if service.review_location(&a.review.target)? != a.review {
            return Err(invalid(
                "Catalog changed after review. Both branches were preserved; activation was refused.",
            ));
        }
        let p = location_path(root, &a.review.target)?;
        scope(&p, &a.review.current_lab)?;
        write(&record(root, j.generation + 1, "applied"), &a)?;
    }
    current_catalog_path(root)
}
fn checkpoint(path: &Path, lab: &store::LabId) -> HostResult<Option<base::Digest>> {
    let c = rusqlite::Connection::open_with_flags(path, rusqlite::OpenFlags::SQLITE_OPEN_READ_ONLY)
        .map_err(io)?;
    let exists:bool=c.query_row("SELECT EXISTS(SELECT 1 FROM sqlite_master WHERE type='table' AND name='cg_reliable_checkpoints')",[],|r|r.get(0)).map_err(io)?;
    if !exists {
        return Ok(None);
    }
    let value:Option<String>=c.query_row("SELECT digest FROM cg_reliable_checkpoints WHERE lab_id=?1 AND namespace='repository' AND key='state'",[lab.as_str()],|r|r.get(0)).optional().map_err(io)?;
    value
        .as_ref()
        .map(base::Digest::new)
        .transpose()
        .map_err(io)
}
impl CatalogBackupService {
    fn snapshot(&self, path: &Path) -> HostResult<Snapshot> {
        scope(path, &self.wire_lab()?)?;
        let dir = directory(&self.root);
        mkdir(&dir)?;
        let staging = dir.join("staging");
        mkdir(&staging)?;
        if fs::read_dir(&staging).map_err(io)?.take(65).count() > 64 {
            return Err(invalid(
                "Inspect abandoned catalog review staging before continuing.",
            ));
        }
        use std::sync::atomic::{AtomicU64, Ordering};
        static NEXT: AtomicU64 = AtomicU64::new(0);
        let artifact = staging.join(format!(
            "snapshot-{}-{}.cglbak",
            std::process::id(),
            NEXT.fetch_add(1, Ordering::Relaxed)
        ));
        struct Cleanup(PathBuf);
        impl Drop for Cleanup {
            fn drop(&mut self) {
                let _ = fs::remove_file(&self.0);
            }
        }
        if exists(&artifact)? {
            return Err(HostError::Unavailable);
        }
        let _cleanup = Cleanup(artifact.clone());
        let library = store::Library::open_read_only(path).map_err(io)?;
        let m = backup::create(&library, &artifact, &self.policy()).map_err(io)?;
        let generation = m.labs.first().ok_or(HostError::Corrupt)?.catalog_generation;
        let rows = |name: &str| m.table_rows.get(name).copied().unwrap_or(0);
        Ok(Snapshot {
            database_hash: convert_digest(&m.database_hash)?,
            database_bytes: m.database_bytes,
            schema_hash: convert_digest(&m.schema_hash)?,
            schema_version: u32::try_from(m.schema_version).map_err(io)?,
            documents: rows("documents"),
            folders: rows("folders"),
            revisions: rows("revisions"),
            mutation_receipts: rows("cg_reliable_mutations"),
            catalog_generation: generation,
            checkpoint: checkpoint(path, &self.lab)?,
        })
    }
    pub fn list_recoveries(
        &self,
        offset: u32,
        limit: u16,
        expected: Option<&base::Digest>,
    ) -> HostResult<collagene_desktop_protocol::CatalogRecoveryPage> {
        if limit == 0 || limit > 200 || offset > 4096 || (offset > 0 && expected.is_none()) {
            return Err(HostError::Wire);
        }
        let mut entries = Vec::new();
        for intent in self.intents()? {
            if let Action::Restore { .. } = intent.action {
                let p = self.receipt_path(&intent.operation_id);
                if exists(&p)? {
                    let r: RecoveryReceipt = read(&p)?;
                    entries.push(crate::convert::<
                        _,
                        collagene_desktop_protocol::CatalogRecoveryReceipt,
                    >(&r)?);
                }
            }
        }
        let bytes = serde_json::to_vec(&entries).map_err(io)?;
        let revision = base::Digest::of(&bytes);
        if expected.is_some_and(|h| h != &revision) {
            return Err(invalid("Recovery list changed; refresh the first page."));
        }
        let end = (offset as usize + usize::from(limit)).min(entries.len());
        let start = (offset as usize).min(entries.len());
        let next = (end < entries.len()).then_some(end as u32);
        Ok(collagene_desktop_protocol::CatalogRecoveryPage {
            entries: base::Rows::new(entries.drain(start..end).collect()).map_err(io)?,
            next_offset: next,
            catalog_revision: revision,
        })
    }
    fn review_location(&self, target: &Location) -> HostResult<Review> {
        let j = journal(&self.root)?;
        if *target == j.current {
            return Err(invalid("That catalog is already active."));
        }
        let recovery = match target {
            Location::Original => None,
            Location::Recovery { id } => {
                let receipt: RecoveryReceipt = read(&self.receipt_path(id.as_str()))?;
                let fresh = if !j.visited.contains(target) {
                    self.verify_recovery(id, &receipt.backup.id, &receipt.backup.archive_hash)?
                } else {
                    let intent: Intent = read(&self.intent_path(id.as_str()))?;
                    self.validate_intent(&intent, id.as_str())?;
                    if intent.action
                        != (Action::Restore {
                            backup_id: receipt.backup.id.clone(),
                            archive_hash: receipt.backup.archive_hash.clone(),
                        })
                        || receipt.operation_id != id.as_str()
                        || receipt.directory_basename != id.as_str()
                        || receipt.database_basename != "library.sqlite"
                    {
                        return Err(HostError::Corrupt);
                    }
                    if self.verify_pinned(&receipt.backup.id, Some(&receipt.backup.archive_hash))?
                        != receipt.backup
                    {
                        return Err(HostError::Corrupt);
                    }
                    receipt
                };
                Some(crate::convert(&fresh)?)
            }
        };
        let current = self.snapshot(&location_path(&self.root, &j.current)?)?;
        let target_snapshot = self.snapshot(&location_path(&self.root, target)?)?;
        let mut r = Review {
            target: target.clone(),
            target_recovery: recovery,
            current_catalog: j.current,
            activation_generation: j.generation,
            current_lab: self.wire_lab()?,
            current_generation: current.catalog_generation,
            current_checkpoint: current.checkpoint,
            current_database_hash: current.database_hash,
            target_snapshot,
            review_revision: base::Digest::of(b"pending"),
        };
        r.review_revision = r.revision().map_err(io)?;
        Ok(r)
    }
    pub fn review_original(&self) -> HostResult<Review> {
        if journal(&self.root)?.pending.is_some() {
            return Err(invalid(
                "Restart to complete the accepted catalog activation first.",
            ));
        }
        self.review_location(&Location::Original)
    }
    pub fn review_recovery(&self, id: &OperationId, hash: &base::Digest) -> HostResult<Review> {
        if journal(&self.root)?.pending.is_some() {
            return Err(invalid(
                "Restart to complete the accepted catalog activation first.",
            ));
        }
        let r = self.review_location(&Location::Recovery { id: id.clone() })?;
        if r.target_recovery
            .as_ref()
            .is_none_or(|v| &v.backup.archive_hash != hash)
        {
            return Err(HostError::Corrupt);
        }
        Ok(r)
    }
    pub fn adopt_recovery(
        &self,
        operation: &OperationId,
        review: &Review,
    ) -> HostResult<Activation> {
        review.validate().map_err(|_| HostError::Wire)?;
        let wanted = Activation {
            operation_id: operation.clone(),
            review: review.clone(),
        };
        let j = journal(&self.root)?;
        // Durable retries are operation+intent bound even after later activations.
        for n in 1..=j.generation + u64::from(j.pending.is_some()) {
            let prior: Activation = read(&record(&self.root, n, "accepted"))?;
            if prior.operation_id == *operation {
                if prior != wanted {
                    return Err(invalid(
                        "Operation ID has a different accepted activation intent.",
                    ));
                }
                return Ok(prior);
            }
        }
        if j.pending.is_some() {
            return Err(invalid("A catalog activation is awaiting restart."));
        }
        if activation_operation_in_use(&self.root, operation.as_str())? {
            return Err(invalid(
                "Operation ID already belongs to a cancelled activation.",
            ));
        }
        // Refuse before review snapshots or journal writes. Exact historical
        // retries above still return their durable outcome at the capacity limit.
        let next_generation = next_generation(j.generation)?;
        if exists(&self.intent_path(operation.as_str()))?
            || exists(&self.receipt_path(operation.as_str()))?
        {
            return Err(invalid("Operation ID belongs to another backup action."));
        }
        if self.review_location(&review.target)? != *review {
            return Err(invalid(
                "Source, target or activation generation changed; review again before adoption.",
            ));
        }
        let dir = directory(&self.root);
        mkdir(&dir)?;
        mkdir(&dir.join("staging"))?;
        write(&record(&self.root, next_generation, "accepted"), &wanted)?;
        Ok(wanted)
    }
    pub fn activation_status(&self) -> HostResult<Status> {
        activation_status(&self.root)
    }
}
pub(crate) fn catalog_generation(root: &Path) -> HostResult<u64> {
    Ok(journal(root)?.generation)
}
pub(crate) fn activation_status(root: &Path) -> HostResult<Status> {
    let j = journal(root)?;
    if let Some(p) = j.pending {
        Ok(Status::PendingRestart(p))
    } else if let Some(c) = j.cancelled {
        Ok(Status::Cancelled(c))
    } else if let Some(a) = j.last {
        Ok(Status::ActiveAfterRestart(a))
    } else {
        Ok(Status::None)
    }
}
pub(crate) fn catalog_session_status(
    root: &Path,
    start: Option<u64>,
) -> collagene_desktop_protocol::CatalogStatus {
    use collagene_desktop_protocol::{CatalogStatus as C, library::Text};
    let activation = activation_status(root).ok();
    let current = catalog_generation(root).ok();
    if start.is_some()
        && start == current
        && !matches!(activation, Some(Status::PendingRestart(_)))
        && activation.is_some()
    {
        return C::Ready {
            activation: activation.expect("checked present"),
        };
    }
    C::Paused{message:Text::new(if activation.is_none(){"Saved catalog activation history failed integrity checks. Editor recovery remains available; catalog access and activation cancellation are blocked."}else{"Saved catalog could not be opened safely or this session is paused. Editor recovery remains available. If an activation is pending, cancel that exact choice; restart after resolving the catalog."}).expect("bounded static catalog status"),activation}
}
pub(super) fn replay_adoption(
    root: &Path,
    operation: &OperationId,
    review: &Review,
) -> HostResult<Option<Activation>> {
    let j = journal(root)?;
    for n in 1..=j.generation + u64::from(j.pending.is_some()) {
        let prior: Activation = read(&record(root, n, "accepted"))?;
        if prior.operation_id == *operation {
            if &prior.review != review {
                return Err(invalid("Activation operation has another accepted intent."));
            }
            return Ok(Some(prior));
        }
    }
    Ok(None)
}

pub(crate) fn cancel_activation(
    root: &Path,
    operation: &OperationId,
    accepted_operation: &OperationId,
    review_revision: &base::Digest,
    generation: u64,
) -> HostResult<Cancellation> {
    let j = journal(root)?;
    // Durable exact retries remain bound even when a later startup has resumed catalog work.
    for n in 1..=j.generation {
        let p = record(root, n, "cancelled");
        if exists(&p)? {
            let c: Cancellation = read(&p)?;
            if c.operation_id == *operation {
                if &c.accepted.operation_id != accepted_operation
                    || &c.accepted.review.review_revision != review_revision
                    || c.accepted.review.activation_generation != generation
                {
                    return Err(invalid(
                        "Cancellation operation has a different accepted intent.",
                    ));
                }
                return Ok(c);
            }
        }
    }
    if activation_operation_in_use(root, operation.as_str())? {
        return Err(invalid(
            "Operation ID already belongs to activation history.",
        ));
    }
    for suffix in ["intent.json", "receipt.json"] {
        if exists(
            &root
                .join("backups/operations")
                .join(format!("{}.{suffix}", operation.as_str())),
        )? {
            return Err(invalid(
                "Operation ID already belongs to another backup action.",
            ));
        }
    }
    let accepted = j
        .pending
        .ok_or_else(|| invalid("There is no matching pending activation to cancel."))?;
    if accepted.operation_id != *accepted_operation
        || accepted.review.review_revision != *review_revision
        || accepted.review.activation_generation != generation
    {
        return Err(invalid(
            "Pending activation changed; refresh its exact receipt before cancellation.",
        ));
    }
    let next = next_generation(j.generation)?;
    let receipt = Cancellation {
        operation_id: operation.clone(),
        accepted,
        generation: next,
    };
    receipt.validate().map_err(|_| HostError::Corrupt)?;
    // No catalog database is opened or rewritten here. Cancellation consumes the
    // accepted generation and retains the current location plus both branches.
    write(&record(root, next, "cancelled"), &receipt)?;
    Ok(receipt)
}

fn next_generation(current: u64) -> HostResult<u64> {
    if current >= MAX_ACTIVATIONS {
        return Err(invalid(
            "Catalog activation history is full; no catalog was changed.",
        ));
    }
    Ok(current + 1)
}

/// Held until all writable host handles close; SQLite readers may coexist.
pub(crate) struct CatalogLease {
    _file: File,
}
impl CatalogLease {
    pub(crate) fn acquire(root: &Path) -> HostResult<std::sync::Arc<Self>> {
        let p = root.join("catalog-session.lock");
        let mut o = OpenOptions::new();
        o.read(true).write(true).create(true);
        #[cfg(unix)]
        {
            use std::os::unix::fs::OpenOptionsExt;
            o.mode(0o600).custom_flags(libc::O_NOFOLLOW);
        }
        #[cfg(windows)]
        {
            use std::os::windows::fs::OpenOptionsExt;
            o.share_mode(0);
        }
        let file = o.open(&p).map_err(|_| {
            invalid("Another desktop catalog session is open; close it before activation.")
        })?;
        private(&p, false)?;
        #[cfg(unix)]
        {
            use std::os::fd::AsRawFd;
            if unsafe { libc::flock(file.as_raw_fd(), libc::LOCK_EX | libc::LOCK_NB) } != 0 {
                return Err(invalid(
                    "Another desktop catalog session is open; close it before activation.",
                ));
            }
        }
        Ok(std::sync::Arc::new(Self { _file: file }))
    }
}

#[cfg(test)]
mod admission_tests {
    use super::*;
    #[test]
    fn activation_capacity_accepts_last_slot_then_refuses_without_overflow() {
        assert_eq!(next_generation(0).unwrap(), 1);
        assert_eq!(
            next_generation(MAX_ACTIVATIONS - 1).unwrap(),
            MAX_ACTIVATIONS
        );
        for generation in [MAX_ACTIVATIONS, MAX_ACTIVATIONS + 1, u64::MAX] {
            assert!(matches!(
                next_generation(generation),
                Err(HostError::Invalid(_))
            ));
        }
    }
}
