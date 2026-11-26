# Native shared-library server

`collagene-server` hosts an immutable GUI asset snapshot and a same-origin binary
`POST /api/library` endpoint plus the retry-safe `POST /api/library/sync`
namespace. It dispatches through `collagene-library-host` into
the embedded SQLite lab library. It does not start external scientific tools,
database executables or git processes. Static assets remain public; library
operations require the server's shared access key.

This is an initial deployment mode with one trusted startup lab and one shared
member key. Individual accounts, membership administration and review remain
future work. The host never accepts a client member identity or lab selector as
authorization. Holders of the current key share access to this lab's folders and
plasmids; restarting with a replacement key revokes the old key for subsequent
requests. Do not hand out a key to someone who should not have lab-wide access.

Local/offline documents remain a separate deployment concern. This server is not
a mandatory remote, nor a complete synchronization engine. Future Tauri local
libraries and user-selected remotes need explicit bindings, retained local forks,
last-common revisions and deliberate pull/push conflict handling. CAS revision
history is groundwork for that behavior, not a claim that it exists already.
Credentials and access keys never belong inside scientific archives. Credential-free
remote bindings may be retained in a local workspace backup, but are excluded
from per-plasmid publication and biological document exports.

## Operator configuration

The binary has three explicit actions:

```text
collagene-server keygen --key-file /private/lab/access.key
collagene-server init --storage-root /private/lab/database --lab-name "Example lab"
collagene-server serve --storage-root /private/lab/database \
  --lab-id lab-ID-PRINTED-BY-INIT --key-file /private/lab/access.key \
  --assets /srv/collagene/static-snapshot
```

The parent directory must exist. `init` creates the selected database directory
if absent, and prints the generated stable lab ID. Serving requires an existing
initialized database and the exact lab ID from trusted operator configuration.
The database filename is fixed to `library.sqlite` under that root. Root and
SQLite sidecar symlinks/non-regular files are rejected. There is no endpoint for
choosing paths, listing directories, or uploading arbitrary filesystem files.

The default listener is `127.0.0.1:8084`, with origin
`http://127.0.0.1:8084`. A different explicit origin may be needed if the browser
uses a hostname instead of the bind address. Non-loopback binding requires both
an explicit `--bind IP:PORT` and `--origin SCHEME://HOST[:PORT]`. For example, behind
a TLS reverse proxy the operator would configure the public HTTPS origin and have
the proxy forward the original `Host` header. Forwarded identity/origin headers
are ignored. No deployment, firewall change or access-key publication happens
as a side effect of building this crate.

The server itself serves HTTP. Use a TLS reverse proxy for any non-loopback
access: a bearer key sent over unencrypted HTTP can be read in transit. Keep the
backend on a protected interface and configure normal proxy connection/header
limits. The canonical public origin is operator configuration; it is never
learned from a request.

`keygen` uses the operating system random source for 32 random bytes, writes 64
hex characters plus a newline, and creates a new file exclusively. On Unix it
creates mode `0600`, and key loading refuses group/other permissions. The key
is never accepted on the command line, printed or included in request logs. Key
loading accepts exactly 64 hex characters with an optional final newline; an
imported key must also come from a cryptographic random generator. Length/encoding
checks alone cannot prove that a manually supplied secret has entropy.

On Unix the database directory must be owner-only (`0700` is appropriate). On
Windows protect the key and database directories using administrator-configured
ACLs; the process inherits those ACLs rather than invoking platform permission
tools. Run under a restricted service account. The key file must be outside the
static asset root, and asset/database roots must be separate trees. Only trusted
operators may modify these startup directories.

## Binary API and authentication

Clients send `Authorization: Bearer <64 hex characters>` and content type
`application/vnd.collagene.library+postcard`. Requests and responses use the
versioned/checksummed library protocol, not JSON. Native document/workspace data
remain binary. No cookies or CORS permissions are provided.

Authentication uses constant-time comparison of the decoded 32-byte secret before
reading or decoding the request body. The API then checks `Host` against the
configured origin and, when present, requires `Origin` to match. `Origin: null`,
foreign origins, multiple origin/authorization values, wrong hosts, unsupported
content encoding/type and query-string API requests are rejected. Non-browser
clients may omit `Origin` but still need the key and correct host.

Every known lab context is also pinned with the
`X-Collagene-Expected-Lab: lab-…` request header. The header is an assertion about
the already selected startup database, never a namespace selector or permission.
After authentication and origin checks, a present malformed/duplicate pin returns
400 and a different lab returns 412 **before the body is polled**. Reusing a
shared access key while reconfiguring the server cannot redirect a pinned write
into another lab. The sync endpoint requires this header on every request
(missing means 428). V1 read/bootstrap requests may omit it to discover the
trusted lab context, but all eight v1 mutation operations require it before
dispatch; missing mutation pins return 428 without a database effect. V1 binary
frames and media types remain unchanged.

HTTP transport failures have generic bounded messages and meaningful statuses:
401 authentication, 403 origin, 413 body size, 415 media type, 408 body timeout,
429 concurrent capacity, and 400 malformed framing. Successful transport returns
a typed binary envelope, including typed operation failures/conflicts and the
trusted lab context. Responses containing library information use `no-store`.
The service emits no request URLs, document names, payloads or bearer tokens to
logs. Startup messages contain only listener/lab/snapshot information.

Native validation and database operations run on Tokio's blocking pool under the
bounded library permit. A disconnected request does not release its permit while
its native task is still running. A mutex serializes this instance's database
handle, and SQLite provides transactional CAS protection across other handles.
Native processing is bounded by the protocol/core/archive limits and never calls
external software. There is no forced cancellation halfway through an atomic
write. If a client loses the connection while saving, retain local edits and
check the current remote head before retrying; do not assume failure means no
commit. In particular, v1 creation requests are not an idempotent job API; use the separate
retry-safe mutation namespace below for durable mutation IDs and receipt lookup.

The default maximum is two concurrent library requests, at most 65 MiB plus the
42-byte frame header per body, and 30 seconds to receive the body. Admission is
immediate: saturation returns 429 rather than building an unbounded queue.
`--max-concurrent N` permits 1 through 16; raising it increases peak native/archive
memory use. `--max-lab-bytes N` sets the storage layer's retained-archive quota
(default 1 GiB, not preallocated). See `LAB_LIBRARY_STORAGE.md` for other database
bounds, durability and backup requirements.

## Retry-safe mutation endpoint

`POST /api/library/sync` accepts only the separate
`application/vnd.collagene.library-sync+postcard` namespace and its checksummed
`CSQ`/`CSR` frames. The existing `/api/library` media type, frames and behavior
remain unchanged. Sending one namespace's frame or media type to the other
endpoint fails; neither route guesses a format.

Each sync request additionally carries `expected_lab` inside its checksummed
envelope. The server checks the borrowed frame prefix before decoding an archive,
and the native adapter independently checks it before any lookup, native review
or mutation. Header/body disagreement or a changed startup lab rejects the
operation without writing; response-context checks remain a final consistency
check rather than the first database-choice guard.

The sync route exposes only `Apply` and `Lookup`. An authenticated lab member
supplies a random mutation ID and an explicitly typed folder/document mutation.
A committed mutation and its exact receipt share one SQLite transaction. If its
response is lost, lookup or resubmission of the same ID/intent returns that
receipt without creating another document/revision. Reusing an ID with another
intent is rejected. CAS conflicts are durable rejections: retry returns the
original conflict snapshot even after the current head changes. Changed intent
or a deliberate retry after resolving a conflict requires a new mutation ID.

Receipt lookup happens before reconsidering a committed historical operation.
For new payload mutations the native adapter validates actual document bytes,
workspace pins and independently approved publication closure before executing
any database write. Malformed/unapproved native payloads return bounded failures
without storing a mutation receipt. It never accepts private checkpoints, local
repository scope creation, SQL, credentials or administration operations over
this endpoint. Both routes use the same trusted startup lab, header/origin
checks, body limits/timeouts and admission pool.

The shared admission lease follows both actual blocking work and the response's
`Bytes` backing allocation. Polling/dropping a response body does not free capacity
while a transport still holds one of its data frames, clones or slices. A slow
reader therefore cannot accumulate arbitrarily many large replies outside the
configured concurrent-request ceiling.

This endpoint provides a reliable remote mutation boundary. Full local/offline
repository synchronization, user-selected remotes and user-facing divergence
resolution remain separate layers; existing browser v1 saves do not silently
become a durable sync outbox.

## Archive publication and history

The native host validates `.cg2` document bytes and native workspace archives,
including the exact selected native document/revision/event. Workspace writes
must declare `DocumentClosure` and pass the server-owned, independently derived
publication closure check. It retains actual scientific ancestry and dependent
run/capsule evidence while rejecting unrelated documents, private later heads,
or unsupported graph-bearing sidecars. The server does not approve whatever
inventory a client sends, and never flattens history to a sequence-only document.

`ExplicitWholeWorkspace` on the wire does not itself authorize sharing everything;
this endpoint rejects it. A future explicit authenticated review interaction can
use the native host's separate trusted scope API. Current automatic publication
requires the exact approved document/dependency closure. Old native archives and
physical status labels do not create molecular authority: downstream reuse still
requires native origin reconstruction and verification.

## Immutable public assets

Startup snapshots only the explicitly selected compiled asset directory through
confined directory capabilities. File/subdirectory symlinks and special files
are rejected; hidden entries are not published. Limits are 4,096 files, 128 MiB
per file, 512 MiB total and 32 directory levels. These are bounds, not allocations.
An `index.html` is required. HTTP reads only in-memory bytes after startup; it
cannot read newly created files, directories, database records or arbitrary host
paths. Unknown, traversal/hidden and backslash paths return 404. Directory
listing and SPA fallback to unrelated paths are absent. Static methods are GET
and HEAD, with BLAKE3 ETags, no-cache revalidation, correct WASM/JS types and
`nosniff`/frame/referrer headers.

Preserve the original hashed paths of workers and their WASM assets in the
snapshot. Include any prior versions still needed by open client tabs before
restarting the server with a new snapshot. The server does not rewrite worker
URLs or inspect current disk files on requests.

## Validation

`tests/router.rs` uses the in-process Axum router only, without binding a listener
or making application HTTP/browser checks. Fixtures cover authentication before
body polling, origin/host/content/body rejection, trusted lab binary dispatch and
CAS, precise archive closure versus unrelated/whole-workspace publication,
body-timeout/concurrency limits, immutable hashed assets and traversal handling,
exclusive/private key generation and symlink refusal. Native host/storage/sharing
fixtures separately exercise actual archived cloning evidence and restart history.
Additional sync router fixtures cover early authentication, namespace/media
separation, shared response admission, lost-reply lookup/replay, nonce mismatch,
exact CAS rejection replay, native payload rejection without a receipt and
history-preserving closure publication without unrelated-document exposure.

The expected-lab router fixture covers wrong/malformed/duplicate pre-body pins,
missing sync and every v1 mutation pin, unpinned read/bootstrap compatibility,
contradictory sync envelope/header pins, absence of a receipt or data effect,
and successful writes under the same correct explicit lab pin.
