# Shared lab database backup and new-destination recovery

These are trusted native operator commands. They use SQLite's online backup API
inside Rust; no external SQLite command, server HTTP request, browser, credential
provider or scientific executable is needed. Offline local editing and existing
remote publication/sync remain independent.

## Commands

Use the exact configured lab ID and an ASCII backup name (letters, digits,
hyphens or underscores,1–64 characters):

```sh
collagene-server backup --storage-root /private/lab --lab-id lab-0123456789abcdef0123456789abcdef --backup-id nightly_001
collagene-server verify-backup --storage-root /private/lab --lab-id lab-0123456789abcdef0123456789abcdef --backup-id nightly_001
collagene-server restore-new --storage-root /private/lab --lab-id lab-0123456789abcdef0123456789abcdef --backup-id nightly_001 --destination-root /private/recovered-lab
```

Backup paths are fixed to
`STORAGE_ROOT/backups/backup-BACKUP_ID.cglbak`; the name cannot escape that directory.
Backup refuses an existing filename. Verify/restore do not open the original live
library database or read the access-key file. Restore refuses **any existing
destination**, including an empty directory. Its parent must already exist and be
private. It creates `NEW_DIRECTORY/library.sqlite`; it never replaces a live
SQLite database, its WAL or an existing recovery copy.

Check the recovered database and configure the server to use that new root before
any later service cutover. The original source database and immutable backup remain
available for rollback. Existing-database replacement, automatic live cutover and
offline rollback orchestration are not implemented by this first restore policy.
Do not copy a backup's SQLite image over a running database by hand.

`--max-backup-bytes N` controls the SQLite image limit (default8GiB, maximum64TiB).
This excludes the at-most64KiB manifest and44-byte container header. Snapshot copy
uses256-page chunks (at most16MiB temporary incremental work); a failed oversize
snapshot is not published. A120-second copy-loop deadline bounds repeated
busy/restart retries; an individual SQLite step can finish after the last deadline
check, and later integrity/hash checks are separate work. This is not a hard total
request timeout. Large databases may require a quieter backup window; no partial
snapshot is represented as complete.

Backup accepts the same optional `--max-lab-bytes` runtime setting as serving and
records it in the manifest. Runtime storage options are not themselves a database
schema or permission. Verify does not reject a valid database merely because its
operator has subsequently lowered the quota below already-stored bytes. Recovery
prints the recorded quota; choose the intended `--max-lab-bytes` when serving the
new root. Other runtime settings, key files, OS credentials and static assets are
not embedded or silently restored.

## What is preserved

The snapshot is the **whole database**, including committed WAL contents:

- Every lab, folder, document, immutable revision and archive BLOB.
- Exact revision parents, metadata, native document/event selection pins and hashes.
- Durable applied/rejected mutation receipts and their original conflict snapshots.
- Private repository checkpoints, outbox/bindings and trusted local-scope mappings
  when present in this database.

The native API policy names the exact complete expected lab-ID set. The CLI
intentionally supports one exact lab ID and refuses a database containing additional
labs. An operator embedding the API can explicitly supply the complete multi-lab
set; backup never silently exports other labs under one lab's label.

This operator backup may contain private local state and unpublished biological
history. It is not a selected publication closure, a sequence export, a shared
workspace attachment or a remotely callable library command. Stored native/archive
BLOBs are checked byte-for-byte and retain their existing scientific provenance;
backup does not parse them into new molecules, invent repaired material or renew
verification/physical authority.

## Format and integrity

`CGLBAK\0\x01` is followed by a little-endian u32 manifest length,32-byte BLAKE3
of the manifest, canonical bounded JSON, and the SQLite image. Manifest schema1
records application/schema versions, exact schema fingerprint, image size/hash,
creation time, configured archive quota, sorted lab/catalog pins and table counts.
Unexpected fields, future versions, noncanonical JSON, truncation, appended bytes,
wrong expected labs and oversized images refuse.

Verification independently compares SQLite schema definitions with the current
native schema and exact supported reliable-storage extension, rather than trusting
a manifest-provided schema alone. It checks SQLite integrity and foreign keys,
revision chains and current heads, each archive's size/hash/native storage revision
binding, folder ancestry, and receipt/checkpoint hashes and accounting. Unsupported
or modified tables/triggers/extensions refuse. Per-table row counts are capped at
12,800,000, lab count at128, and existing per-record/storage limits also apply.

These checks detect corruption and accidental foreign-database restoration; they
are **not a digital signature or authentication of a malicious backup author**.
Use backups from the trusted operator-controlled source. Hashes alone cannot prove
that an attacker did not replace both a valid database and its manifest.

## Publication and filesystem boundary

Staging directories and files are private, source/target symlinks are refused, and
new files are created exclusively. After all checks, a same-filesystem hard link
publishes the completed single backup file atomically and fails if its destination
already exists. There is no overwrite/copy fallback on filesystems without hard-link
support. Restore likewise links a fully checked, synced image into its newly
reserved private directory. Temporary staging is removed on normal failure/drop;
process death may leave a `.collagene-backup-*` directory that is not a completed
backup and can be reviewed/removed by the operator.

Unix checks current effective ownership and owner-only permissions; staging uses
0700 directories/0600 files. File data is synced before publication and parent
directories are synced on Unix. Windows inherits the operator-selected directory's
ACL; this module has no independent Windows DACL inspector. File data is synced,
but equivalent directory-entry power-loss durability is not claimed by the standard
filesystem API. Keep the configured storage/backup/recovery parents private using
native OS permissions. The storage root is trusted native configuration, never a
path accepted from a web/MCP/library request.

An I/O error after publication may leave a valid final backup or recovered database.
Do not delete or overwrite it automatically: verify the backup/recovery outcome
before retrying with another new destination. Local filesystem privilege against
the owning account is outside this boundary.

## Source validation

Five native SQLite fixture groups cover live-WAL/second-connection commits,
folders/revisions/receipts/checkpoints after restart, immutable backups and independent
recovery, corruption/truncation/manifest hashes/future versions/foreign labs,
resource bounds, unsupported schema/BLOB corruption and Unix permissions/symlinks.
Two CLI helper groups exercise the actual backup/verify/new-restore commands and
path confinement using temporary fixture directories only. No live user database,
access-key contents, application HTTP, browser or GUI checks are involved. Compiler
and fixture execution are performed by the parent batch after source freeze.
