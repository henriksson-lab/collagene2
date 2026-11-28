# Desktop remote-library protocol

`collagene-desktop-sync-protocol` defines a portable Rust/WASM metadata protocol,
separate from the DTQ/DTR native storage protocol. It provides current DSQ/DSR
schema **5** frames for the native
`collagene_rpc` dispatcher. This crate contains no networking, filesystem access,
repository execution, native credentials or scientific algorithms. A native
service must implement those operations over the trusted window's repository.

The 42-byte header contains magic `DSQ\0` or `DSR\0`, little-endian schema u16,
postcard payload length u32 and a BLAKE3 payload checksum. Request frames are at
most 64 KiB and response frames at most 3 MiB, **including the header**. Checksums
detect corruption; they are not authentication. Decoding rejects trailing bytes,
unknown versions, invalid identifiers and oversized declared collections before
their allocation. Serialization uses a capped buffer. All listed pages have at
most 200 rows. There are no biological archive bytes in either direction.

Every request has `request_id`, a desktop `OperationId`, and a typed command.
Every reply repeats both pins. `decode_response_for` validates the reply shape,
operation and command-specific identity/page pins. `validate_in_remote` additionally
compares the entire remote config with the UI's known config. Remote read success
requires a context containing the actual configured remote and returned lab.
The expected lab is an assertion checked by the native HTTP/sync boundary before
effects, never a client grant of membership. Native host authority supplies the
local lab, window, repository, credential account and all transport permissions.

## Commands and results

* `ListBookmarks` is a fully offline, read-only local repository query. It returns
  bounded snapshot references (local catalog document/revision) and native branch
  references (native document, old revision/event, undo/redo counts and incoming
  revision/event). A native graph bookmark has no stored catalog association;
  clients must not infer one from names, sequence equality or the active editor.
  These references convey documentary history, not physical material authority.
  The reply has no remote context and includes the exact repository checkpoint.
  Continuation requests require that checkpoint; concurrent changes return
  `BookmarkChanged` so a client refreshes instead of mixing pages. Reopening
  retained archives uses the existing local immutable revision-history workflow.

* `Remotes` / `AddRemote` list or create nonsecret configurations: generated
  `rem-` ID, bounded name, HTTP(S) origin and expected lab. Native execution must
  parse and canonicalize origins, enforce HTTPS or its explicit literal-loopback
  exception, and apply the existing no-redirect transport policy. No credentials,
  query parameters, paths or fragments belong in an endpoint.
* `SetCredential`, `DeleteCredential`, `CredentialStatus` and `ClearSession`
  operate on that configured remote's host-derived credential scope. Retention is
  explicitly `Session` or `Remember`. Status is `Missing`, `Session` or
  `Remembered`; replies never contain a token. `ClearSession(None)` forgets only
  in-memory session values, not persisted vault values.
* `Browse`, `Search` and `ListRevisions` reuse library snapshot paging and its
  parent, revision, catalog generation and continuation validation. These are the
  only general remote requests accepted by this protocol; an arbitrary server
  mutation or uploaded archive cannot be submitted.
* `CloneDocument` imports the fetched document into a declared local folder and
  returns the actual source document/revision and new local catalog document.
  The service must validate/preserve the native closure and reject an already
  bound source instead of replacing a private local head. The ordinary single
  clone follows the fetched current head; recursive clones use immutable scanned
  revisions as described below.
* `QueuePush` names a saved local catalog document, its exact local library
  revision and a configured remote. The host rechecks that revision before
  deriving publication closure from saved native data. A queued reply pins that
  same revision. No browser working archive, arbitrary target or metadata marker
  grants material authority. The service captures the freshly loaded repository
  checkpoint internally and uses native transaction CAS; it must not substitute
  the latest local head after a stale request.
* `ResumePush` reuses the durable job's exact mutation ID, bound intent and
  immutable publication. Its summary includes intent digest, original catalog
  generation and applied/rejected outcome. A lost response resumes by receipt
  lookup; it cannot create a new mutation. `Pending` reports bounded job summaries;
  `CancelUnattempted` and `DiscardRejected` retain their existing native guards.
  Pending/status list offsets are best-effort local views, not a remote catalog
  snapshot and not permission to clear unseen jobs.
* `Pull` requires the exact saved local library revision. `CleanFastForward`
  maps to the native guarded fast-forward policy; `ExplicitCheckout` is an
  explicit selection change that preserves the entire displaced local branch
  and ancestry. Neither mode invents a local history event, flattens a saved run,
  or silently replaces a native root identity. The service returns the actual
  saved local document after native merge/validation and transactional CAS.
* `BeginTreeClone`, `ResumeTreeClone`, `RetryTreeClone` and `TreeJobs` expose
  bounded existing recursive jobs. Roots distinguish a virtual lab root from
  a selected folder with its declared parent. Local mount destination is a
  separate field; it never implies a remote move. Scans are snapshot-stable;
  imports fetch the exact scanned immutable revisions even if the remote advances.
  Responses expose progress, invalidation and blocking issues without sending
  whole trees or archives. Partial jobs remain partial across restarts.
* `QueueTreeFolderMetadata` and `QueueTreeDocumentMetadata` invoke the native
  mount-aware helpers. They preserve local placement versus remote parentage,
  generation conflicts and the durable outbox. No raw metadata mutation or private
  checkpoint is supplied by the browser.

DTQ cancellation/status controls remain unchanged and address these operations'
window-scoped `OperationId`s. DSQ has no second cancellation state machine. Native
execution must not claim cancellation after durable commit linearization; durable
sync-job cancellation is separately restricted to an unattempted outbox entry.

## Secrets and errors

`Secret` accepts at most 2048 printable nonspace ASCII bytes, matching the native
vault/token bound. Its custom serde validates borrowed text before allocation;
its custom `Debug` prints only `[REDACTED]`, including inside request debug output.
It has no `Display` or implicit string conversion. Only native credential handling
should call `expose_secret`. This portable DTO does not promise allocator
zeroization. Hosts must consume promptly, avoid logging frames/serialized requests,
and transfer into the native zeroizing secret representation. No response type
contains a secret. Native vault access errors must remain explicit and must not
fall back to plaintext persistence.

Errors are typed safe codes, a guarded current local document on revision conflict,
or the existing bounded remote failure. They cannot return arbitrary native errors,
filesystem paths or token content. Transport errors are not evidence that a remote
mutation failed; durable receipt lookup determines that outcome.

Portable fixtures cover malformed/corrupt/oversized frames, secret decode bounds
and redaction, request/operation/lab/origin/page binding, saved-revision queue pins,
and mount/progress semantics. Native adapter and real repository integration tests
belong to the service crate. Adding this wire crate alone does not enable remote
management in the current desktop or public web preview.

## Reviewed native remote removal

Append-only PrepareRemoteRemoval, ConfirmRemoteRemoval, DiscardRemoteRemovalReview
and CheckRemoteRemovalOutcome commands require frame version 3, as do their new
results and credential-forget/post-forget-conflict failure codes. Existing
publication commands retain minimum version 2; original commands retain version
1. Legacy request-prefix decoding and response-for-frame encoding remain usable.
No credential, arbitrary path, biological archive or remotely supplied authority
is added to these DTOs. Reviews contain full configured remote scope, trusted
local lab, checkpoint and exact binding/completed-job counts. Receipt schema 1 is
a distinct local record format, carried inside the current schema-5 transport.

The native host retains at most four opaque unpredictable removal IDs, expiring
the oldest. A review is re-derived before approval acceptance and key-forget;
IDs from another native instance, expired/discarded reviews and changed scope
refuse. A missing vault entry is idempotent success. Other vault failures retain
configuration/bindings and truthfully report that the session may already have
been cleared. After successful key-forget, SQLite unlink remains independently
CAS guarded. A conflict reports key-forgotten/configuration-retained; an uncertain
I/O outcome is checked against its durable receipt before declaring doubt.
No remote server delete is performed. CheckOutcome is read-only and binds the
full old remote configuration to its native local-context receipt.


## Remote human discussion

Current desktop sync IPC schema **5** carries
`Discussion { remote, request: library::discussion::Request }` and a checked
`Success::Discussion` response. Earlier IPC schemas are rejected without a
migration or fallback decoder. Discussion uses its own current schema 2 on the
native HTTP transport; molecular document and reliable-sync frames retain their
separate current versions.

Requests contain only bounded human text, exact remote document/revision pins,
page/generation controls and an expected authenticated actor pin for appends.
No role claim, bearer, scientific archive, native path or arbitrary URL is carried.
The native host resolves the saved remote and uses its separately managed
credential. Response checks bind the outer desktop operation ID, request ID,
complete remote context and expected lab, then the exact nested discussion request,
trusted member, revision, action and generation/page sequence.

The 3 MiB response frame cap admits a maximum 200-event discussion page with
8 KiB text bodies and remote context. Discussion requests retain the ordinary
64 KiB IPC frame cap and the tighter 16 KiB service request bound. The dedicated
HTTP route is fixed by the `Discussion` transport channel; no endpoint is inferred
from incoming bytes. Invalid channel/version/framing is refused before requesting
a credential.

Discussion is explicit online work, not a durable plasmid push. There is no
outbox entry, replay, clone/pull side effect or scientific history attachment.
After a lost reply, refresh the exact discussion; the unchanged expected generation
refuses a duplicate compare-and-append. A changed actor similarly refuses rather
than impersonating the previously reviewed author.

## Explicit folder-clone recovery

`RestartTreeClone { job }` restarts an invalidated or blocked scan against the
current remote catalog only while no folder or document has been imported.
An active scan, completed job or partially imported snapshot cannot be reset.
The reply is `TreeUpdated` with empty scan state; no remote context or native
biological mutation is claimed.

`InspectTreeParentChange { job }` returns `TreeParentReview`: the exact blocked
`TreeStatus` and the current native `Folder`, including name, parent and generation.
The folder must be the clone's reported `LocalParentChanged` identity. This is
local metadata review, not a remote server assertion or a reusable authority token.
The UI presents those values and requires explicit acknowledgment.

`AcceptTreeParentChange { review }` rereads and compares the full job status and
folder metadata before accepting the native operation. Repository authority also
checks blocked/completed-scan state, exact folder membership and current generation
when committing the checkpoint. A stale review is refused. Acceptance updates the
acknowledged generation and clears the resolved issue while leaving the clone
blocked. It does not import or overwrite any node, resume the job, or publish a
folder move remotely. `RetryTreeClone` remains a separate explicit action.

Recovery survives host restart because the job checkpoint owns its acknowledged
local generations. Cancellation before acceptance does not change it; interrupted
response delivery can be resolved by refreshing the persisted clone state. Normal
native operation IDs and receipt handling remain unchanged. All preceding desktop
sync encodings are rejected without readers or migrations.
