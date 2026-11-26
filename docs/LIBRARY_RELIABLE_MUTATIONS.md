# Reliable native library mutations

`lab-library::mutations` adds an optional, namespaced SQLite extension to the
existing lab library. Base schema/API version 1 remains usable. The new module
owns its extension tables and initializes them transactionally; it does not
change base tables, archive formats, scientific behavior or authentication.
The module is registered in the native library and exercised by the server,
offline repository and desktop host.

## Idempotency and exact receipts

A trusted host submits `Mutation { id, operation }`. IDs are client-generated
128-bit random nonces encoded as `mut-` plus 32 lowercase hex digits;
`Library::new_mutation_id()` provides native generation. They are correlation
identifiers, not authorization credentials. `LabAccess` still comes from trusted
host authentication. All eight base write operations are supported: create,
rename, move and delete-empty folder; import, save, rename and move document.

For a new `(lab, mutation_id)`, `apply_reliable` executes the operation and stores
its receipt in the **same SQLite transaction**. It returns
`ApplyResult { receipt, replayed }`. If a reply is lost, resubmit the same ID and
intent or call `lookup_reliable` with its intent digest. An exact retry returns
the original receipt, including its generated IDs, revision, catalog generation
and original conflict snapshot. It does not recompute an outcome from later
state or rerun a rejected operation after conditions change. Another intent
under the same ID returns `IntentMismatch`; another lab has a separate namespace.

The canonical digest algorithm is shared with the pure binary sync protocol:
BLAKE3 of `b"collagene-library-mutation-v1\0"` followed directly by the UTF-8
`serde_json` encoding of the externally tagged operation. Variant and field order
are declared in `IntentOperation`. A payload becomes a descriptor in this exact
field order: `kind`, `archive_hash`, `archive_size`, `selection`, `sharing`.
Selection uses the existing `document_id`, `revision`, `head_event` order. The
operation's names, IDs, CAS bases, size, byte hash, native selection and publication
scope are bound. The nonce and private checkpoints are excluded. Payload shape,
size and actual byte hash are checked before opening the write transaction;
large binary data are never expanded into JSON integer arrays for this digest.

Receipts contain either `Applied` with the exact operation result or `Rejected`
with a bounded domain failure and, when relevant, the original current folder or
document snapshot. Content/metadata conflicts, missing items, duplicate names,
nonempty deletion, invalid folder geometry and storage limits are durable
rejections. Database, integrity and private transaction failures abort completely
without pretending a usable receipt exists. To deliberately retry a rejected
operation against changed state, create a **new mutation ID** and updated intent.

Each operation runs under an internal savepoint. Domain rejection rolls its
partial effects back before the receipt is written; an outer failure rolls back
both operation and receipt. CAS, lab isolation, sibling names, cycle/depth checks,
immutable revision history and archive quota accounting match the base API. A
concurrent duplicate creation generates one document and one receipt, while the
other caller replays that exact identity.

The ledger is append-only, with immutable update/delete triggers. Receipts are
bounded to 128 KiB each, 100,000 per lab and 64 MiB total receipt bytes per lab.
Hitting a limit refuses a new mutation atomically; existing receipt lookup/replay
remains available. No automatic pruning forgets old IDs and silently makes their
retries execute again. Any future compaction must preserve the explicit
idempotency boundary and be designed as a versioned retention policy.

## Private native repository checkpoints

`CheckpointKey { namespace: "repository", key }` addresses native-only state.
Keys are at most 128 ASCII letters/digits, `_` or `-`. Values contain an explicit
positive schema version and opaque bounded binary bytes. These types have no
wire serialization implementation. No client SQL, path, actor, remote endpoint
or credential operation is accepted. Remote credentials/access keys must remain
outside checkpoints and scientific archives; this API is not a credential vault.
The native repository owns decoding and validating each checkpoint schema.

`read_checkpoint` returns its value and digest. `commit_checkpoints` applies a
batch atomically, comparing each expected digest (`None` means absent) before
changing any row. `None` as the new value deletes a checkpoint. Digests bind the
lab, namespace, key, version and exact bytes. A stale expected value returns
`CheckpointConflict`; earlier writes in the batch do not become visible.
Checkpoint changes are private bookkeeping and do not increment the public
folder catalog generation.

For an offline edit that must queue a push or update a binding alongside its
local native revision, use:

```text
apply_reliable_with(access, mutation, |receipt| -> checkpoint writes)
```

This trusted, synchronous finalizer runs inside the same outer transaction after
the exact outcome is known. It can therefore bind an outbox/state checkpoint to
a newly generated document/revision ID. Returning an error, stale checkpoint,
over-limit value or duplicate key rolls back mutation, quota charge, receipt and
all checkpoint changes. No network call or external process belongs inside this
callback. An exact replay does **not** rerun the callback: its original private
updates already committed, and newer repository state must not be overwritten.
An acknowledgment-only transition can use `commit_checkpoints` without changing
DNA. The durable repository still owns remote push/pull state, divergence and
last-common revision policy; this layer does not infer synchronization decisions.

Limits are eight checkpoint writes per transaction, 64 MiB per value, 128 MiB
aggregate values per lab and 4,096 checkpoint keys per lab. Zero-length values
still consume a key. These are ceilings, not preallocations. Application-level
bounded decoding remains necessary before constructing a checkpoint value.

`local_repository_access(namespace, name)` is an explicit **native owned-file**
factory. It atomically creates or rediscovers a stable namespace-to-lab mapping
for an offline repository across restarts. It is never a network authentication
method and must not be exposed by the shared server. Namespace names are bounded
ASCII identifiers, with at most 128 local scopes per database. Existing scopes
retain their original lab and display name on restart.

## Validation

The new `tests/reliable.rs` fixtures cover lost-reply replay across restart and
later edits, canonical intent/nonce mismatch, lab isolation, exact original
conflict snapshots, rollback on finalizer or checkpoint CAS failure, replay
without callback rerun, private checkpoint batch/restart behavior, simultaneous
duplicate creation and preservation of base tree/CAS/quota invariants. Tests use
opaque storage payloads; host native loaders/publication checks still validate
real scientific bytes before any remote mutation is submitted.
