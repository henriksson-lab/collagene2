# Durable offline plasmid repository

`crates/repository` provides an owned-file native Rust repository over the embedded
SQLite lab library. Local plasmids, folders, complete private design workspaces
and immutable local history work without any remote connection or authentication.
The crate starts no listener, launches no external program and depends on no
remote service for cloning analysis. Tauri IPC, asynchronous UI orchestration and
OS credential-store integration are separate adapters.

The repository is registered in the native workspace and used by the desktop
host. The Tauri shell and native remote-control service use separate repository
connections, keeping slow remote requests independent from editor recovery saves.
See [reliable mutations](LIBRARY_RELIABLE_MUTATIONS.md),
[sync protocol](LIBRARY_SYNC_PROTOCOL.md), and
[native closure import](LAB_LIBRARY_IMPORT.md) for their separate contracts.

## Owned local persistence

`Repository::open(path)` opens the real SQLite library and obtains its stable
`offline-plasmids` local scope through `local_repository_access`. Reopening the
same file rediscovers that scope and its versioned checkpoint; a remote login,
remote lab locator or external configuration file is not needed. This is an
owned-file API, never a network authentication mechanism. The caller must protect
the file using normal local OS permissions.

`mutate_local` accepts the existing typed native folder/document operations:
create, rename, move and empty-folder delete; native import and immutable content
save with exact base-revision CAS. `browse`, `folder`, `document`, and
`open_document` use local storage only. Every payload is checked by the existing
native document/workspace loader and its exact selected document/head/event pins.
Local full workspaces deliberately retain private/unpublished branches, UI state,
scientific events, run attachments and saved `FromRun` evidence.

Every successful mutation and the repository checkpoint commit in the **same**
SQLite transaction using `apply_reliable_with`. Generated local IDs and heads are
available to the trusted checkpoint callback before commit. A checkpoint conflict
rolls back the DNA/folder mutation as well as its receipt; it cannot leave a
half-created local record or lose another connection's binding changes. After a
stale result the caller uses `reload` and reconsiders its operation, rather than
automatically overwriting a newer head.

Local data and repository bookkeeping occupy separate native tables. Remote
bindings, queued network frames and branch bookmarks never become biological
workspace attachments or publication dependencies. A local workspace payload's
`ExplicitWholeWorkspace` storage descriptor describes its private archive scope;
it is **not** authority to publish that archive.

## Independent selected remotes

`add_remote` creates a stable local `rem-…` identity with a display name, endpoint
locator and expected remote lab ID. There may be independent bindings to several
remotes for the same local plasmid or folder. Endpoint URLs must be HTTP(S), with
no user information, passwords, query strings or fragments. Credentials are not
fields of these DTOs. The injected `RemoteTransport` supplies authenticated HTTP,
IPC or another bounded transport without persisting secrets in repository state.

Document bindings retain the exact remote document record, remote immutable head,
last common native molecule/event selection, and the corresponding immutable local
library revision. Folder bindings retain remote folder IDs and expected metadata
generations. Native IDs/events, local database revisions and remote database
revisions remain separate identities.

`RemoteTransport::exchange(remote, channel, bytes)` is invoked only by explicit
clone, pull or resume operations. Replies must match the binary request and the
selected remote lab context. Transport errors contain no credential-bearing raw
HTTP messages. Ordinary local reads, saves and analysis do not invoke transport.
These are blocking native methods: a UI adapter should run them outside its render
thread and respect local snapshot/CAS guards at completion.

## Durable push and interrupted replies

`queue_push` derives a new exact selected-document/dependency closure from the
private native workspace **before** queueing. Unrelated private plasmids, UI
sidecars, credentials and repository checkpoints are excluded. A native `.cg2`
snapshot remains that format; no absent ancestry is fabricated. The queued sync
operation is either an import or a save against the binding's exact last remote
base revision.

The complete bounded binary Apply frame, stable mutation ID, selected native pin,
immutable local source revision and destination configuration are persisted before
any network call. Frames are stored in at most two 64 MiB private checkpoint
chunks, with a full-frame digest. State plus chunks commit atomically. The native
128 MiB aggregate checkpoint quota bounds retained outbox data; it is a ceiling,
not a preallocation. Too many/large queued archives fail before a partial queue is
committed.

`resume_push` first performs a durable receipt lookup:

- A matching receipt is validated against the **entire** queued mutation, not
  merely its ID, before updating any binding.
- Unknown resends the same mutation ID and immutable intent/body. A server that
  committed an import before its response was lost returns its original generated
  document ID; it does not create a duplicate document.
- A CAS/domain rejection is retained as the exact durable rejected receipt. The
  outbox remains visible and the common base is not advanced. `discard_rejected`
  removes only an acknowledged no-effect operation so a user can explicitly pull
  or resolve a branch before queueing a new intent.
- An interrupted/ambiguous outcome remains queued. It cannot be discarded and
  recreated under a new ID. `cancel_unattempted` is available only before any
  transport attempt; cancellation of an unknown in-flight outcome would require
  a separate server cancellation contract.

Accepted receipt, binding changes and outbox-chunk removal commit atomically in
one local checkpoint transaction. If that commit is interrupted or conflicts,
the persisted queue still contains the same ID and a later lookup recovers it.
A later private edit does not alter queued bytes. Acknowledging the old queued
commit updates the last-common pin while leaving the newer private working head
untouched; another explicit push can then publish that later work.

No automatic last-writer-wins retry, blind base advancement or sequence-based
claim of branch equivalence occurs.

## Clone and pull

`clone_from` validates the pinned remote payload, imports it locally, and binds the
generated local document to its remote record in the same transaction. A repeated
clone of an already bound remote document is refused and directs the caller to
pull; it cannot silently create duplicate local copies after restart.

`pull` fetches the current remote record and narrows its workspace to the exact
selected native dependency closure. It uses the native immutable graph merge
against the complete current private local workspace. Unrelated local documents,
UI state and newer private dependency heads remain. A different native root ID
is a separate identity and must be cloned separately. Default strict selected
head/event conflict requires a deliberate decision; explicit expected-pin
fast-forward or checkout follows the existing import policies.

For checkout/fast-forward, the complete displaced local navigation bookmark is
persisted in the repository checkpoint **in the same transaction as the new local
content head**. It is not merely returned as a transient UI signal. Both immutable
branches and exact run evidence remain. Restarting can recover the bookmark and
open old immutable local revisions.

Bare `.cg2` snapshots have no native event graph to prove ancestry. They use the
separate explicit `pull_snapshot(expected_local_revision)` operation: the old
immutable local body is durably bookmarked before the head changes. An archive
must use graph-preserving `pull`; the snapshot path cannot flatten a workspace.
Neither path restores physical material from detached historical summaries.

## Folders and metadata

`clone_folder` clones only the selected folder metadata after locating it in a
validated, generation-bound remote parent listing. Children/documents are chosen
and cloned separately; this supports selective import without silently downloading
an entire remote tree.

`queue_folder` explicitly creates a remote folder or pushes one rename/move at a
time using its bound expected metadata generation. If both name and parent changed,
rename is queued first; after its receipt updates the generation a subsequent
explicit queue can move it. Parent folders must already be bound/pushed for that
remote. An offline local empty-folder deletion can be followed by
`queue_deleted_folder`; remote nonempty/conflicting state remains a typed rejected
receipt. Nothing recursively deletes documents or revisions.

`queue_document_metadata` similarly sends explicit name/location changes with the
expected remote generation. Content saves do not silently rename or move remote
records. A content push to an existing remote preserves that remote location even
if a local folder move remains private. New remote imports require the selected
local parent to be bound, or the local root.

This is not complete recursive folder synchronization. Automatic tree discovery,
remote folder pull/merge, move-cycle conflict presentation, remote removal/editing,
deletions/tombstones beyond empty-folder operations, remote-to-remote transfer and
multi-device ancestor negotiation remain explicit later work. The present APIs
provide real offline persistence, selected clones, immutable document pull/push,
retry-safe mutation recovery and generation-bound folder metadata operations.

## Bounds, storage and remaining adapters

Repository metadata is capped at 4 MiB, 128 configured remotes, 1,000 outbox entries
and 1,000 graph/snapshot bookmarks each; storage additionally bounds all native
entities, revisions, checkpoint bytes and reliable receipts. Native archives retain
their existing 32/64 MiB payload and 128 MiB expansion profiles. Large binary
payloads stay binary in checkpoint frames, never JSON integer arrays. Private
SQLite/WAL files need normal supported backup and disk monitoring.

The transport/UI adapter still needs authenticated remote sessions and OS-managed
credentials, Tauri filesystem/IPC permissions, progress/cancellation presentation,
branch/bookmark review and explicit remote selection controls. Local analysis does
not wait for those adapters or require a remote account. Do not describe these
selected sync primitives as a completed git clone/merge implementation.

## Focused real fixtures

`tests/repository.rs` uses durable local SQLite files, the actual native host/sync
dispatchers and real remote SQLite stores. It covers offline CRUD/full private
workspace persistence across reopen; a committed remote creation with a lost
response followed by restart/lookup recovery without duplicate creation; later
private edits during an old queued push; independent remotes and offline folder
rename/move generations; divergent concurrent pushes and durable explicit pull
bookmarks; real saved PCR evidence through push/clone/reopen/fresh native replay;
and local checkpoint conflicts rolling back a folder mutation entirely.

No browser, network HTTP or compilation checks were run during the coordinated
GUI publication freeze. Source is formatted/frozen for the parent's focused native
batch after the new storage/sync modules are registered together.

## Native HTTP transport

`collagene_repository::http::HttpRemoteTransport` is the native HTTP implementation
of the same injected transport interface. It uses the existing locked/cached
`reqwest` **0.12.28** with blocking Rustls TLS. Construct it with a trusted
`CredentialProvider` and bounded `HttpOptions`, then pass it to explicit repository
clone, pull or resume methods. It is not constructed or called by ordinary local
CRUD. A Tauri/UI caller should run these blocking operations outside its render
thread.

A credential provider returns a `BearerToken` for the selected `RemoteConfig`.
The token and optional `MemoryCredentials` provider deliberately implement neither
Debug nor serialization. Memory credentials are bound to local remote ID, canonical origin and expected
lab; changing a locator or lab cannot reuse a token. They can
be removed/cleared or replaced using `credentials_mut`; they are not written to
repository checkpoints, URLs, workspaces or archives. An OS credential store can
implement the same provider later. Tokens are bounded visible ASCII without
whitespace or control characters, and the Authorization header is marked sensitive
so standard reqwest request diagnostics redact it. Token strings are not promised
to be locked or zeroized memory; deployment should use its appropriate OS/session
credential controls.

The configured endpoint must be a **trusted origin**, such as
`https://lab.example:9443`, with only `/` as its path and no URL user information,
query or fragment. The adapter constructs exactly:

- ordinary v1 POST `/api/library` with the library postcard media type;
- reliable sync POST `/api/library/sync` with the separate sync postcard media type.

Both Content-Type and Accept are fixed to that channel. Every request also sends
`X-Collagene-Expected-Lab` from the mandatory selected `RemoteConfig.lab`. The HTTP
host checks that pin against its authenticated database scope before polling the
body or allowing any mutation. Response-only checking would be too late if an
origin/key had been retargeted to another lab.

The borrowed request frame must have its correct magic/version, bounded exact
length and checksum before the credential provider is queried. Sync additionally
uses the borrowed `request_prefix` reader to compare its mandatory envelope
`expected_lab` with `RemoteConfig.lab` before credentials or execution, without
deserializing/copying a potentially 64 MiB archive. The native sync dispatcher
checks that same pin before receipt lookup or mutation. Ordinary v1 retains its
existing binary shape and is scoped by the mandatory HTTP header at that boundary.
The server still performs full typed decoding, and the repository validates reply
correlation, lab context and full mutation receipt binding before committing
anything locally.

HTTPS uses normal certificate validation. HTTP is rejected by default. An explicit
`allow_loopback_http` option permits only literal loopback IPv4/IPv6 addresses,
for example `http://127.0.0.1:8085` or `http://[::1]:8085`. A hostname such as
`localhost` does not qualify by DNS convention; use a literal address instead.

Redirects, environment/system proxies, automatic request retries and transparent
response decompression are disabled. A 3xx response is an error, never a new
credential-bearing request to its Location. HTTP failures, credential failures and
interrupted reads produce sanitized classified errors; no raw server error body,
reqwest error URL or credential is retained in them. `last_error()` exposes the
last static classification for a UI adapter. Mutation retries remain controlled
by the repository's durable Lookup/same-ID Apply logic.

Only status 200 with the exact expected binary media type is accepted. Unexpected
Content-Encoding and duplicate/malformed key headers fail before reading a body.
Ordinary replies are bounded by the existing library frame profile; sync replies
by the smaller 256 KiB receipt profile plus framing. Advertised oversized lengths
are rejected before allocation/read. A streamed or dishonest missing length can
request at most the cap plus one overflow probe; overflow bytes are not appended.
Truncated/mismatched content lengths and read failures withhold the reply.

Default connect timeout is 5 seconds and the request/read budget is 60 seconds;
configuration is bounded to a 60-second connect timeout and 300-second budget.
An overall deadline is checked before and after every body chunk to stop
indefinite trickling. Reqwest's blocking reader also bounds each individual read;
a read already in progress may take its configured timeout after the overall
budget expires. These are bounded blocking calls, not an instant-cancel promise.
No transport configuration grants a remote authority over local analysis.

Five additional offline HTTP boundary fixture groups construct requests and use
fake executors/body readers. They verify origin/loopback rules, exact channel
paths/media, sensitive credential redaction, rejection before credential access,
redirect/error/media/encoding handling without reading those bodies, database
header pins on both channels and a mismatched sync lab pin rejected before
credential access, and bounded unknown/lying/truncated body lengths and expired
deadlines. They perform **no
real HTTP request, listener, DNS, TLS handshake or browser smoke check**.


## Saved-head guards for native controls

`local_context` reads the owned database identity; it does not grant authority to
a message. `checkpoint_version` is an internal service snapshot pin. Token-free
`outbox_metadata` exposes IDs, saved source revisions, attempted status and durable
rejections without exposing upload frames or archive bytes.

`queue_push_pinned` and `pull_pinned` check the saved library revision displayed
by the native UI before selecting source payloads or starting remote requests.
The native service captures the checkpoint itself. Helpers reload and check that
internal checkpoint, then retain the existing transaction/CAS checks so a writer
that intervenes during network I/O cannot overwrite newer local work. No global
checkpoint is required from the frontend for ordinary reads.
