# Desktop remote-library protocol

`collagene-desktop-sync-protocol` defines a portable Rust/WASM metadata protocol,
separate from DTQ/DTR v1. It provides DSQ/DSR schema 1 frames for a future native
`collagene_rpc` dispatcher. This crate contains no networking, filesystem access,
repository execution, native credentials or scientific algorithms. A native
service must implement those operations over the trusted window's repository.

The 42-byte header contains magic `DSQ\0` or `DSR\0`, little-endian schema u16,
postcard payload length u32 and a BLAKE3 payload checksum. Request frames are at
most 64 KiB and response frames at most 1 MiB, **including the header**. Checksums
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
