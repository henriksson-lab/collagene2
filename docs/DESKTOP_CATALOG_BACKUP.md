# Native saved-catalog backup service

`CatalogBackupService` is a framework-independent native adapter over the checked
SQLite [backup format](LAB_LIBRARY_BACKUP.md). The shell supplies its configured
private root and exact local lab ID. Requests supply only typed operation/backup
IDs; they never choose filesystem paths, source databases, credentials or a remote.
No external program, browser, HTTP request or key provider participates.

## Scope and ownership

The service opens only the existing, nonempty selected managed catalog (initially
`ROOT/repository.sqlite`), validates
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

RestoreCopy only creates a checked managed copy. A separate explicit review/adoption
can select that copy for the next native startup; there is no live SQLite-handle
swap or rename over `repository.sqlite`. The private root stays the same, preserving
editor recovery paths and OS/session credential scope. The original catalog remains
at its original path as a rollback source.

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


## Reviewed restart-bound adoption and rollback (DTQ/DTR3)

The appended commands are `ListRecoveries`, `ReviewRecovery`, `ReviewOriginal`,
`AdoptRecovery` and `ActivationStatus`. Existing ordinary desktop commands retain
v1 and backup/copy commands retain v2 floors. Header downgrade cannot carry new
variants, including a failure response to a v3 request. Inputs contain typed
managed IDs, hashes and a checked review, never paths or arbitrary database bytes.

`ListRecoveries` pages every completed owned restore receipt, with complete-set
hash continuation binding and at most200 rows. These are historical receipts,
not fresh integrity claims. Selectors are only `Original` and `Recovery { id }`.
`ReviewOriginal` explicitly reviews rollback to the preserved original catalog;
`ReviewRecovery` also supports a known copy that was active earlier and has since
accumulated private edits. No branch is overwritten or flattened.

Review freshly checks native single-lab/schema scope and original restore intent,
receipt and backup identity. A never-activated recovery must still match its exact
immutable backup image. A previously-active branch is permitted to contain its
subsequent legitimate writes. Both source and target receive checked **online
SQLite snapshots**, including WAL, through the native backup API. The temporary
snapshot containers are removed after extracting database/schema hashes, byte
counts and saved-catalog counts. The descriptor additionally binds source/target
catalog generations, private repository checkpoint pins and current activation
generation. It never presents an ordinary primary-file hash as a complete WAL
snapshot. Native backup size/time limits also bound this review work.

The client supplies the reviewed descriptor verbatim. Acceptance independently
reconstructs that descriptor before publishing an immutable checksummed
`catalog-activation/GENERATION.accepted.json`, via private staging, fsync and
exclusive no-overwrite publication. Generations are contiguous, at most4096;
unknown files, gaps, orphan applied records, reused operation IDs and changed
intents refuse. Same-operation/exact-intent retry returns its historical original
receipt even after later selections; it never creates another transition. A
historical receipt is not proof that its target is currently pending or active:
clients query fresh ActivationStatus.

The shell holds exclusive catalog admission while accepting; running catalog or
sync jobs make it busy. Accepted selection is **pending restart**, not live. Shell
admission rechecks durable pending state even if the response is lost, and the
native library dispatcher refuses old-catalog work. Working-editor reads/saves,
current/previous recovery, preferences, exports and operation controls remain
available. Exact acceptance is discoverable on reload without a reply receipt.

At startup, before opening writable catalog handles, an exclusive native
`catalog-session.lock` lease prevents a second writable desktop process using this
root. Startup recomputes source and target snapshots and checks the accepted
activation generation before publishing the identical `GENERATION.applied.json`.
Late old-catalog writes, changed target bytes, wrong lab/schema, symlinks or damaged
records refuse cutover and retain both branches. A crash before applied publication
retries the same choice; after publication startup follows that selected managed
path. Subsequent writes there are normal catalog writes, so later startups validate
native identity/lab/private paths rather than enforcing the old pre-activation
image hash. There is no silent fallback to another branch if the selection journal
or its selected database is damaged.

Repeated recovery adoptions and rollback-to-original use the same review/CAS flow;
all prior catalog files remain in place, together with the immutable selection
journal. The root, working.sqlite/current/previous autosaves, preferences and OS
credential namespaces do not change. No in-process SQLite-handle swap is claimed.

MCP preload uses `Repository::open_read_only`: SQLite READ_ONLY and query_only,
existing schema/scope/state validation, no initialization/migration/default inserts.
It can read the active saved catalog alongside the GUI, but writes through that
handle fail in SQLite. New sessions refuse while a selection awaits restart.
Already-preloaded MCP sessions remain independent memory snapshots and cannot save
AI results into either branch. Restoring old outbox state does not send network
requests automatically; publication remains explicit and receipt/CAS guarded.

Additional authored fixtures cover preserved saved-PCR bytes and replay after
activation, independent editor sentinels, exact retry/restart, active branch writes,
repeated adoption and rollback, stale source/generation refusal, lost-intent refusal,
exclusive desktop leases, read-only MCP coexistence/mutation refusal, and v3 versus
legacy frame floors. The parent executes coordinated tests; this source worker ran
no compiler, GUI, browser, network or live-data checks.

## Editor-only recovery startup and exact cancellation (DTQ/DTR4)

A failed catalog activation does **not** prevent opening the native working editor.
The host opens its independent working store before attempting catalog startup and
keeps an optional catalog connection. If a pending source/target changed, a catalog
is missing/unsafe, or the activation journal is corrupt, the window can still load
current/previous working archives, preserve drafts/history, save working data and
preferences, export files and use operation controls. Library/sync factories refuse
without trusted catalog availability. No empty replacement catalog or fictitious
lab context stands in for the missing connection.

The old Bootstrap DTO/prefix remains unchanged. New `BootstrapRecovery` request and
result require v4 and carry real working metadata/preferences, an optional genuine
library context, and `CatalogStatus::Ready` or `Paused`. Paused always has no catalog
context. Its bounded public diagnostic and optional activation record distinguish
an exact known pending choice from a journal that could not be trusted. Cached
recovery bootstrap reads refresh their observational status so an earlier pending
or available reply cannot conceal later acceptance/cancellation. The small control
fast path explicitly recognizes the new ordinary bootstrap without admitting
unknown selectors or downgraded headers.

`CancelActivation` supplies the exact accepted operation ID, review revision and
baseline activation generation. With a valid journal and that choice still pending,
the native root-only service appends `GENERATION.cancelled.json`. The terminal
receipt embeds the exact accepted record and new cancellation operation. It consumes
the pending generation, retains the current catalog location, and never opens,
repairs, merges, deletes or rewrites either database. A changed-source or corrupt-
target failure therefore does not prevent cancellation. A damaged/ambiguous journal
refuses cancellation; editor recovery is still available. An already-applied choice
cannot be undone by cancellation—use normal reviewed rollback for that case.

The same operation and exact cancellation intent replay their durable receipt after
lost replies or restart. Different pins/nonce reuse refuse. Applied and cancelled
terminal records cannot coexist. A cancelled generation can be followed by another
reviewed activation without renumbering or rewriting history. The current process
remains catalog-paused after cancellation, even if its previous catalog location is
now selectable; a fresh restart establishes new handles/context. Native generation
checks and shell admission prevent old catalog/sync writes during that interval.
Fresh startup still validates the retained catalog's exact lab scope, including
when cancellation was the first terminal record. Working editor and credential
storage roots never change.

The backup facade can be constructed using only the trusted private root, shared
lease and window control when no catalog can be opened. It exposes exact pending
cancellation (and safe diagnostics through recovery bootstrap), while all catalog-
data operations refuse. Recovery metadata requests never read credential stores or
transmit biological data.

`desktop-host/tests/catalog_recovery_startup.rs` contains five authored temporary
native fixture groups: changed-source editor startup with exact cancel/cache/restart
and preserved private branches; corrupt-target cancellation without touching target
bytes; corrupt-journal editor access with cancellation refusal; v4/control-prefix
and legacy-bootstrap floors; and refusal to cancel an already-applied selection.
Current and caller-supplied previous private workspace archives/preferences remain
byte-exact. These are parent-run fixtures; the source worker ran no compiler,
browser, application HTTP, external program or live-user-data checks.
