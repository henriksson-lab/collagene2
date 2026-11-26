# Explicit native remote removal

Removing a configured database is local bookkeeping, not deleting remote folders
or plasmids. `prepare_remote_removal(remote_id)` returns a `RemovalReview` with
read-only exact remote configuration, trusted local lab, captured repository
checkpoint, and document-binding/folder-binding/completed-clone counts. Its
private digest also binds complete affected binding snapshots, the tree-job index
and every exact completed job checkpoint. Opaque approval is not deserializable.

Any pending mutation for that remote blocks removal, including unattempted,
rejected and uncertain uploads. Resolve it through existing explicit cancel,
receipt lookup/resume or acknowledged-rejection handling first. Any unfinished,
blocked or invalidated tree clone for that remote also blocks removal. No outcome
is silently discarded. Completed selected-remote clone bookkeeping can be retired;
other remotes' jobs and all local cloned folders/documents remain intact.

The host presents the review and obtains explicit acknowledgement before calling
`review.approve()` and `remove_remote_reviewed`. The latter freshly derives the
same scope. It atomically removes only the chosen remote configuration, its
binding map entries and owned completed tree-job checkpoints/index entries. It
writes a secret-free durable `RemovalReceipt` in the SAME SQLite transaction.
Local document/folder rows, all immutable biological archive revisions, native
history, branch bookmarks, private working/recovery files and other remotes are
never deleted. The method performs no remote API request and no vault access.

`remote_removal_receipt(remote_id)` reads a prior outcome. A repeated matching
opaque approval returns that exact receipt after restart without another write.
The receipt retains the full old remote scope and review digest for host lifecycle
recovery; it contains no biological document/event IDs or credentials. It proves
local metadata retirement, not scientific verification, remote deletion or a
vault operation. Native hosts must order credential cleanup explicitly. SQLite
and an OS vault do not share an atomic transaction: a key-forget failure must not
be reported as success, and a later SQLite CAS failure may leave configuration
and data present after its key was forgotten. The frontend must describe that
truthfully. Credentials never enter a repository/archive/receipt.

The storage-only `commit_checkpoint_retirements` API permits at most 128 deletion
records plus three replacement records (131 total), preserving ordinary
`commit_checkpoints`' eight-write limit. It uses the same trusted LabAccess,
expected-digest CAS, duplicate-key validation, namespace/version/byte bounds,
aggregate quotas and single SQLite transaction. Repository code supplies only
its exact state/index/receipt records and owned completed job deletions; no SQL
is accepted from callers. A late conflict rolls back the complete batch.

Starting a new tree clone asserts its captured repository state checkpoint in
the same transaction as job/index creation. Thus a paused source validation
cannot create a new job for a concurrently retired remote. If clone creation
wins first, removal's index/job checks invalidate its earlier review instead.

Focused fixtures are lab-library `checkpoint_retirements` (full 131-record
success, ordinary-limit preservation, late-CAS rollback, count limits and lab
isolation) and repository `remote_removal` (real native binary remote bindings,
completed clone cleanup, all local archives/history preservation, independent
remote preservation, restart/replay, pending/unfinished refusal and concurrent
review invalidation). No HTTP, GUI, user vault or live data is touched. The
implementation worker formats/checks diffs only; parent runs compilation/tests.
Host/DSQ/UI lifecycle integration is parent-owned separate work.

## Native service and user controls

The DSQ/DSR schema-3 host now implements this lifecycle. A maximum of four opaque
reviews is retained in native memory; oldest expiry bounds abandoned reviews.
Fresh native repository scope is checked before accepting and forgetting a key.
The exact old configured scope is used; no endpoint is accepted in the approval.
Already-missing vault records succeed idempotently. A refused forget retains
configuration and bindings, while reporting possible session clearing. Held
publication reviews for that remote expire BEFORE the credential change, even
when its vault refuses. A successful unlink consumes its review and clears other
held reviews for the same configuration. Other remotes' publication reviews stay.

The native Databases accordion includes a compact removal disclosure. It prepares
a source-bound counts/configuration review and requires an explicit acknowledgement
before “Forget key and unlink database”. Changing the selection invalidates the
review. A last-outcome check uses only the selected or previously captured known
remote; there is no raw ID/path field. Interrupted responses retain that captured
remote for checking the secret-free durable receipt. Both review and outcome
validation bind the complete configuration; receipt validation additionally checks
its exact reviewed digest/counts/local context when available. Parent list state
removes only the acknowledged remote ID. Local designs remain untouched.

Focused integration targets: sync-protocol `removal_frames` (new minimum floor 3,
old 1/2 floors, exact response binding and bounded/versioned receipts), sync-host
`removal_controls` (owner acceptance, exact scoped key-forget, durable restart
outcome, locked vault, pre-vault stale source, post-forget CAS change, foreign/old/
session-reset review IDs, Missing idempotency and failed-forget publication-review
expiry), and pure web `remote_removal_model` (complete configuration/selection and
receipt binding). Injected adapters perform no real vault, HTTP or GUI calls.
These source additions are formatted/diff-checked; parent batch validation is
required before claiming runtime acceptance.
