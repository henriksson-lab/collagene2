# Native desktop host

`crates/desktop-host` is a framework-independent native host over the real offline
`Repository` and a private working-state SQLite database. It exposes one bounded
binary interface suitable for the Tauri command `collagene_rpc`. No listener,
HTTP operation, login, external executable, frontend filesystem path or frontend
identity assertion is part of this interface.

The separate [desktop protocol](DESKTOP_PROTOCOL.md) describes its portable raw
byte frames and DTOs. The [offline repository](OFFLINE_REPOSITORY.md) remains the
local plasmid/folder/history model; scientific native archive formats and saved
`FromRun` replay rules remain unchanged.

## Shell-owned authority

The native shell chooses an application-data root and calls:

```rust,ignore
let mut host = DesktopHost::open(configured_private_root)?;
let access = host.authorize_window()?;
let control = host.control();
let reply_bytes = host.dispatch_raw(&access, request_bytes)?;
```

`WindowAccess` is an opaque, cloneable native capability branded to that host
instance. It has no constructor or serialization implementation available to the
frontend. A capability from another host instance, or a revoked capability, is
rejected. Authorization is a **shell** decision: issuing this capability does not
authenticate arbitrary web content. The shell must restrict the IPC command to
its approved local application window/origin and normal native application ACL.
A request body never chooses a window, lab membership, storage root or path.

`DesktopHost` is Send and intended to run behind the shell's native mutex/job
admission gate; `WindowAccess` and `HostControl` are Send + Sync. The shell should
run ordinary work outside the UI thread. These APIs are deliberately independent
of Tauri so the same persistence/validation rules are native-testable.

The configured root is created/canonicalized natively. Only fixed
`repository.sqlite` and `working.sqlite` names, plus SQLite's corresponding
sidecars, are used. On Unix a new final root is created0700 and both primary
files0600 before SQLite opens them; existing root/database/WAL/SHM records must be
owned by the effective OS user and have owner-only modes (0700 directory,0600
files). A shared0755 root, readable0644 database or symlink is rejected rather than
silently adopted or chmodded. Checks use the OS effective user, not environment
identity claims. The working database also uses SQLite's no-follow open flag.

Windows/non-Unix builds retain regular-file/no-symlink checks, but Rust std has no
portable DACL inspector. The native shell must choose and enforce a user-private
AppData ACL there; the host does not falsely claim Unix-mode checks validate
Windows access control. Protect the configured root against other filesystem
writers on every platform. No filename, SQL, root, remote endpoint or credential
from a raw IPC command reaches a path join or file open.

Bootstrap returns the native local library context, working-slot metadata and
small view preferences; it does not disclose private filesystem paths. A native
window capability is local application authority, not a remote lab login. Local
plasmid analysis/storage does not require a network account.

## Local plasmids and folders

`Command::Library` uses the existing bounded library v1 request shape. Reads use
the same scoped local SQLite library. Writes are routed through
`Repository::mutate_local_outcome_with_gate`, so native data, repository state,
bindings and relevant private checkpoints retain the repository's atomic mutation
semantics. The narrow structured-outcome APIs preserve the existing
`mutate_local` wrapper while keeping NotFound, AlreadyExists, NotEmpty and CAS
conflicts typed instead of parsing error strings.

For local v1 mutations, the operation ID reaches the same native receipt
transaction as a stable mutation nonce. Its exact full intent and local repository
scope are checked before replay. A matching retry returns its original typed
outcome without rerunning mutation, checkpoint or cancellation callback, even
after the transient reply cache is pruned or the desktop process restarts. Reusing
the nonce for another intent is rejected. Frontend operation IDs therefore must
be globally unique random128-bit IDs and must be retained when retrying the same
request; they are not reusable sequential window counters. Working-slot and
preference writes retain their separate event-version CAS behavior.

Private local workspace writes use the ordinary native loader and exact selected
native document/head/event pins. They retain the **complete private workspace**,
including unpublished branches, unfinished drafts, UI state and run evidence.
Publication projection is not used for local persistence. Names are virtual
labels; a slash in a plasmid/folder display name is not a filesystem path.

Ordinary generation/base-revision CAS still applies. Library rejections return
bounded public typed errors; database, native parser and filesystem internals are
not returned over IPC. Remote synchronization remains a separate explicit native
repository operation and is not exposed by these desktop commands.

## Full private working state

The fixed working database has an explicit application/schema identity and
refuses unrelated databases. It uses SQLite foreign keys, WAL, `synchronous=FULL`
and bounded busy waits. Current and previous each hold at most **64 MiB** of
native workspace archive bytes. Their combined logical archive budget is
**128 MiB**, without preallocating that space.

Working slots, the independent **1,024-byte** preference table and repository
checkpoints/outbox occupy separate fixed storage scopes. Preferences and remote
bookkeeping do not consume the two archive slots or require a quota increase.
SQLite indexes/WAL and filesystem overhead are separate from logical archive
bounds; supported database backups must account for WAL sidecars.

A working version is a digest of the version domain, monotonic generation, prior
version and both current/previous archive descriptors. Saving identical bytes
creates a new event version. Thus an older native window cannot overwrite a later
save merely because the archive's content hash happens to match again. There is
one CAS version for the **pair**, and a separate CAS version for preferences.

The storage API matches the browser adapter's semantics:

- `SaveCurrent { expected_version, archive }` changes **current only**, preserving
  the existing previous slot.
- `ReplaceWorking { expected_version, current, previous }` atomically installs
  **both caller-provided archives**. The provided backup may contain unsaved
  in-memory work; the host does not substitute the on-disk current archive.
- `ReadWorking { slot, expected_version }` returns one exact pinned slot, never
  both 64 MiB archives in one reply. Metadata/bootstrap does not read giant BLOBs.
- `RecoverPrevious { expected_version }` independently validates previous,
  installs it as a new current event, and leaves the known-good previous backup
  intact. A corrupted current cannot be copied into and destroy that backup.
- Preferences use independent `ReadPreferences`/`SavePreferences` and their own
  expected version. An empty preference value can explicitly clear it.

Every archive write is validated with `workspace::load_workspace` before it can
commit. Both replacement archives must pass validation; one invalid backup leaves
the entire pair unchanged. Read checks stored archive length/hash and native
loading as well. No checksum, imported duplex summary or retained historical report
is turned into live scientific material by storage. Fresh scientific verification
and exact originating-run reconstruction are still required.

The pair replacement and its version update occur under one `BEGIN IMMEDIATE`
transaction. Reopening after a crash observes an old or new complete pair. Slot
queries guard actual BLOB length before allocating; another process cannot create
an unbounded read merely by corrupting a size field. Metadata and slot body reads
share a SQLite read transaction, so a concurrent save cannot mix generations.

## Cancellation and operation state

A cloneable control handle is intentionally separate from the host's long-lived
work mutex:

```rust,ignore
if let Some(reply) = control.dispatch_raw(&access, request_bytes)? {
    return reply; // small Cancel / OperationStatus handled without the work lock
}
// None is not validation or authorization of an ordinary request.
// Schedule host.dispatch_raw under the normal bounded native work gate.
```

For frames larger than the 1 KiB control cap this fast path returns None without
hashing, decoding or copying large archive data. Small control frames are fully
checked. Large ordinary requests undergo full checked decoding only on the work
thread. A replacement may contain two 64 MiB archives, so raw request admission
is bounded by the protocol's 129 MiB profile; replies by 65 MiB.

Operation IDs are scoped to the actual native window capability. Cancelling an
unknown/not-yet-registered job creates a cancellation tombstone, so an ordinary
job queued behind the native mutex can still be stopped before it runs.
Cancelling Running requests a stop before the commit gate. Read-only results can
be discarded after validation if cancellation arrived meanwhile.

Immediately before durable commit, the host atomically changes Running to
**Committing**. In repository writes this gate runs inside the same transaction
before its prepared checkpoint/receipt commits; rejecting the gate rolls back the
entire native mutation. Working-state writes use the same linearization rule.
Once Committing is accepted, cancellation reports Committing or a completed
terminal state and **does not undo a committed write**. A worker-local acceptance
marker survives window revocation/removal of transient registry entries, so a
committed success or storage failure is not relabeled Cancelled after revocation.
A precommit revocation still cancels and rolls back the operation. Committing does
not promise success: a later storage error can still fail normally.

The registry is bounded to 2,048 operation entries and 4 MiB of recent small reply
cache. It preserves uncaught cancellation tombstones, prunes completed entries
when necessary, and does not cache large archive replies. A recent identical
operation frame can receive its cached reply; reusing an ID for different data is
rejected. This cache is bounded in-process RPC bookkeeping. Local library writes
also use the durable scoped mutation receipt described above; cache loss cannot
turn a repeated Import into a second record. Persisted pair/preference CAS remains
the protection for working storage against stale windows and restart retries. `revoke_window` disables that native
capability and removes its transient records.

## Deliberate initial scope

This host exposes native persistence/local library operations and cancellation.
It does not expose arbitrary SQL, arbitrary filesystem commands, remote sync
credentials, shell execution, or a catch-all scientific job executor. The current
browser scientific worker remains available; a native data-only analysis bridge
can later add an explicit job/size contract without giving its output storage or
material authority.

The shell still owns window/origin ACL, raw command registration, admission limits,
application-data path selection and presentation of conflicts/recovery. The WASM
adapter owns raw Uint8Array framing and retains the pair/preference versions from
bootstrap or successful replies. It must surface a conflict instead of silently
retrying with another window's newer version.

## Focused native fixtures

The host's integration/unit fixtures exercise real temporary SQLite repositories
and private native archives:

- Caller-provided unsaved backup versus current-only save, independent preference
  CAS, identical-byte event-version changes and full restart persistence.
- Atomic rejection when either proposed replacement archive is invalid; corrupted
  current detection and recovery without damaging known-good previous data.
- Window-brand/revocation isolation, cancel-before-dispatch tombstones and the
  bounded control lane.
- Typed local folder/library conflicts and complete private workspace roundtrip.
- Real saved PCR history through working save/reopen and fresh native material
  reconstruction.
- A cancellation inside a real repository transaction rolling back its local
  effect, and too-late cancellation after Committing leaving a successful durable
  working save intact.
- Unix root/file owner-only creation and rejection of insecure existing modes.
- Import replay after rename/restart, when its old display name is free again;
  exact receipt replay and intent-mismatch rejection prevent duplicate imports.
- Commit acceptance surviving window revocation, plus cancelled receipt rollback
  and replay never rerunning its callback.
- Compile-time Send/Sync requirements for native shell scheduling.

No GUI, browser, HTTP or application smoke checks are used. Source is formatted
and frozen for the parent's focused native host/protocol checks.

## Scoped external native jobs

The same branded window authority and bounded operation registry can supervise a
separate validated native protocol, such as DSQ synchronization, without sharing
the local repository/working-file mutex:

```rust,ignore
match control.begin_job(&access, &operation_id, validated_frame_digest)? {
    JobStart::Cached(checked_reply) => return Ok(checked_reply),
    JobStart::Run(job) => {
        // Reversible planning/validation may happen before acceptance.
        // The service calls job.accept() at its indivisible acceptance boundary.
        // It must retain durable pending/outcome state if remote effects are unknown.
        let checked_secret_free_reply = run_with_gate(|| job.accept().is_ok())?;
        job.finish(OperationState::Completed, &checked_secret_free_reply)?;
    }
}
```

`JobStart` and `ScopedJob` are exported native types, not wire authority. The
caller must validate its protocol framing, pins and request first, and bind the
supplied digest to the complete frame **including its protocol namespace**.
A different frame under an already used operation ID cannot replace the original
entry or cached reply. `ScopedJob` has no public constructor, Clone, Debug or
serialization implementation. Its only authority is one existing branded
window, operation ID and validated intent digest.

`accept(&self)` linearizes with cancellation and window revocation under the
registry lock and records acceptance. Repeated acceptance of the same guard is
idempotent. Pre-accept cancellation/revocation prevents acceptance;
`cancelled_before_acceptance()` lets a service discard an unaccepted proposal.
After acceptance, a late cancellation/revocation cannot relabel a committed or
possibly remote effect as rolled back. `accepted()` records that distinction.

`finish(self, terminal_state, checked_secret_free_reply)` permits only Completed,
Failed or pre-accept Cancelled. The caller has already checked its response
protocol and excluded credentials; the generic registry does not parse DSQ or
hold a password/URL/path service. An accepted job cannot finish as Cancelled.
Dropping an unfinished guard terminalizes known pre-accept cancellation as
Cancelled and other abandoned jobs as Failed. Failed after acceptance describes
an unconfirmed outcome; it does not claim that an accepted effect was undone.
Native durable receipts/outboxes remain authoritative across restart.

The registry stays bounded to the existing 2,048 entries and 4 MiB total cached
reply bytes; only replies up to 32 KiB are cached. Oversized checked replies still
get a terminal state but no retained data copy. Completed records can be evicted,
so this in-memory facility does not replace durable mutation receipts. Revoked
terminal records and unseen cancellation tombstones are discarded; in-flight
accepted records remain until their owner completes, then normal eviction can
reclaim them. DTQ Cancel/OperationStatus addresses DTQ and DSQ jobs uniformly,
within the same window scope.

The shell uses separate ordinary local-DTQ and remote-DSQ admission lanes. A slow
remote request must not prevent separate `working.sqlite` autosaves. The external
service owns a separate Repository connection; SQL revision/checkpoint CAS still
rejects intervening edits rather than permitting stale catalog overwrites. No
additional frontend command, filesystem grant or OS-vault access follows from
this job guard.

Focused native unit fixtures exercise cross-protocol cancellation, pre-accept
revocation, post-accept cancellation/revocation, foreign-window denial,
changed-intent replay protection, abandoned guards and the existing cache bounds.
They need no HTTP, GUI, browser or live vault.
