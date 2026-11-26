# Native saved-catalog backup service

`CatalogBackupService` is a framework-independent native adapter over the checked
SQLite [backup format](LAB_LIBRARY_BACKUP.md). The shell supplies its configured
private root and exact local lab ID. Requests supply only typed operation/backup
IDs; they never choose filesystem paths, source databases, credentials or a remote.
No external program, browser, HTTP request or key provider participates.

## Scope and ownership

The service opens only the existing, nonempty `ROOT/repository.sqlite`, validates
its native database identity and exact single local-lab scope, and keeps a separate
SQLite handle. It never opens, reads, initializes or overwrites `working.sqlite`.
This is a backup of the **saved local catalog**: folders, document revisions,
complete stored private archives, native history/evidence and reliable repository
receipts/checkpoints/outbox state. Unsaved editor content and current/previous
working autosaves are not claimed to be included.

The shell gives backup work its own admission lane/mutex, independent of working
file persistence. SQLite's online backup API provides a consistent committed
snapshot while ordinary local writers continue. The shell supplies cancellation
and operation-status controls; once mutation acceptance occurs it must not report
that accepted filesystem work was rolled back. List replies remain read-only and cancellable before acceptance. Create, Verify
and RestoreCopy accept before native work because Verify can recover a missing
owned receipt after checking an already-published artifact.

Native API:

```rust
CatalogBackupService::open(root, &local_lab)
service.create(&operation_id)
service.list(offset, limit, expected_catalog_revision)
service.verify(backup_id)
service.restore_copy(&operation_id, backup_id, &expected_archive_hash)
```

`backup_id` is exactly the original `op-` plus32 lowercase hexadecimal digits,
validated again inside the service. The root protocol can compare an optional
expected hash after `verify`; `restore_copy` itself requires and checks the exact
chosen archive hash before publishing a recovery intent or copy.

## Managed storage and receipts

Fixed locations are:

- `ROOT/backups/OPERATION_ID.cglbak`: the immutable checked backup container.
- `ROOT/backups/operations/OPERATION_ID.intent.json`: durable exact action/scope.
- `ROOT/backups/operations/OPERATION_ID.receipt.json`: checksummed publication reply.
- `ROOT/backups/operations/staging/`: private, bounded atomic record staging.
- `ROOT/recoveries/OPERATION_ID/library.sqlite`: a new checked recovery copy.

No live cutover, source replacement or automatic application of recovered state is
performed. Recovery is deliberately a new managed copy; it is not a new active
Tauri root. In particular, the service does not rename that file to
`repository.sqlite` or change startup configuration. A later explicit cutover policy
can review and adopt a copy while retaining the old source for rollback.

Receipts expose safe basenames, IDs, file/image hashes and byte sizes, creation
time from the actual backup manifest, schema identity, exact local lab/counts and
a **request-time checkpoint observation**. That checkpoint is observed before
snapshot creation; a writer may commit while SQLite obtains its consistent
snapshot. It must never be labelled the exact snapshot checkpoint or current
editor state. Native IDs/evidence in stored archives remain unchanged; backup does
not renew scientific verification or turn imported summaries into physical stock.

## Exact retries and interruption

An action's versioned, local-lab-bound intent is synced and atomically published
before its snapshot/recovery work. Record files use private staging, fsync and an
exclusive same-filesystem hard-link publication, so a partial JSON record is not
installed as a completed intent or receipt. Backup artifacts use the same
no-overwrite checked native backup publication. There is no overwrite fallback.

A repeated creation operation behaves as follows:

- Published artifact and receipt: freshly verify both before returning the original
  receipt, even if the live source has since changed.
- Published artifact with a lost receipt: freshly verify the owned artifact and
  recover its original metadata/receipt.
- Intent without an artifact: refuse automatic re-snapshot. The original snapshot
  may never have completed, and a later source must not silently replace it. The
  user can explicitly start a new operation/backup ID.
- Wrong action, different lab, corrupt record, changed artifact or orphaned file:
  refuse; never reinterpret it as permission to create or overwrite data.

A recovery intent binds the selected backup ID and **whole-container hash**. Retry
can reconstruct the same immutable source into an absent managed child or verify
an already-published copy before returning its original receipt. A recovery copy
must have exactly one private `library.sqlite` file whose full image hash/length
match the freshly verified pinned backup. Unexpected files, WAL/SHM, partial roots,
changed copies and cross-action ID reuse refuse. Nothing is deleted to force a
retry. Opening a recovery with a writable SQLite adapter may change its header or
create sidecars; that is then a changed copy, not an identical idempotent retry.

An error after publication may leave a valid artifact/copy. Retain the operation
ID and inspect/retry it; do not infer rollback or start a silent replacement.
Abandoned native `.collagene-backup-*` staging directories remain unpublished and
are recognized only as private staging, never as complete entries. Normal drops
remove them; filesystem administration may review crash leftovers. Record staging
is bounded to64 abandoned files before further writes refuse. Total managed
operation intents are bounded to4096.

## Listing and integrity checks

Listing returns all owned creation intents, including unfinished entries with no
receipt. It pages at most200 entries, in stable operation-ID order, and binds every
continuation to the hash of the complete current entry set. Changed lists reject
stale continuation; the201st item is not silently hidden. An unfinished entry is
not a completed backup. Resuming its creation uses the **same entry ID as the
operation ID**, rather than assigning a new ID under an old name.

Listed receipts are bounded historical publication metadata with current
regular-file/length/scope checks. Listing does not rehash every multi-GiB backup;
the UI must not call the list a fresh integrity verification. Explicit `verify`
and every recovery run full native manifest/schema/SQLite/FK/revision/archive/
receipt/checkpoint checks and stream the complete container hash before returning
a freshly checked result. No archive bytes cross IPC. Foreign filenames, unknown
IDs, orphaned receipts/artifacts, wrong labs and symlinks refuse.

Unix checks effective ownership and exact0700/0600 private modes and syncs directory
entries. Windows uses the existing shell-selected user-private ACL policy; this
module does not claim independent DACL verification or power-loss directory-entry
durability beyond the documented native backup boundary. Filesystem privilege
against the owning account is outside this native/window trust boundary.

## Focused source fixtures

Six embedded native fixture groups exercise:

- A real saved PCR/FromRun workspace with private saved design data, backup/restart,
  unchanged original bytes, fresh native replay from a recovery and an unreadable-as-
  SQLite `working.sqlite` sentinel that remains untouched.
- Lost receipt, source advance, original-byte acknowledgement, unpublished intent
  refusal, action-ID reuse and wrong selected recovery hash before durable writes.
- Complete201-entry pending-intent pagination and stale continuation rejection.
- Invalid IDs, foreign filenames/labs and archive corruption without replacement.
- Unix managed-directory symlink refusal without accessing its target.
- Real scoped raw IPC, exact frame/cache binding, fresh artifact checks even for
  cached success replies, explicit hash mismatch, cancellation tombstones, foreign
  windows, wrong commands and the16KiB request limit.

Fixtures use temporary native storage only. No compiler, browser, application HTTP,
GUI, live user database or credential-vault checks are performed by this source
worker; the parent runs the coordinated native checks after registration/freeze.


## Native raw-IPC facade

`DesktopCatalogBackups::open(root, &local_lab, control)` owns the service in its own
mutex and shares only the host's branded window/operation controls. Its
`dispatch_raw(&WindowAccess, bytes)` accepts only a complete validated DTQ
`CatalogBackup` command, at most16KiB; replies are at most256KiB. The full validated
frame digest binds the operation registry/cache, including its namespace, command,
request ID and expected hashes. Responses are converted to bounded typed metadata,
validated against the request and encoded at the request's supported protocol
floor. No large native/archive byte vector crosses the command boundary.

Cached creation/verification/recovery success is rechecked against actual owned
artifacts or copy hashes before acknowledgement; an old in-memory reply cannot
hide subsequently corrupted files. Cached list replies remain historical metadata
for that exact read request. Optional Verify hash and mandatory recovery hash are
checked before any missing receipt is recovered. Wrong-hash verification cannot
write a replacement receipt. An accepted completed publication remains the outcome
even if window disposal prevents later cache/status bookkeeping; no successful
write is relabelled as cancelled. The root shell routes small cancellation/status
controls outside this service mutex and keeps backup admission separate from
working-file saves.
