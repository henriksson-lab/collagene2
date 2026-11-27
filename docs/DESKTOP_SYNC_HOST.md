# Native desktop synchronization service

`crates/desktop-sync-host` is the framework-independent native executor for the
metadata-only [DSQ/DSR protocol](DESKTOP_SYNC_PROTOCOL.md). It operates on an
injected, real offline `Repository` and a native transport/credential adapter.
The shell supplies window/IPC authority separately; no DSQ request chooses a
filesystem path, local repository, window identity or credential-vault namespace.

This service does not replace local editing, native working-state storage or the
scientific worker. Local plasmids/designs continue to work without remote login.
It implements explicit remote controls over **saved local catalog records**; it
cannot publish an unsaved in-memory frontend draft or a caller-supplied archive.

## Native API and adapter

```rust,ignore
let mut service = SyncHost::production(saved_repository, http_options)?;
let reply = service.dispatch_raw(&request_bytes, &mut || {
    // Native shell marks this operation accepted/committing before returning true.
    mark_native_operation_accepted()
})?;
```

`SyncHost<A: SyncAdapter>::new(repository, adapter)` supports native injection.
Construction compares the adapter's trusted local repository lab ID with the
repository's own native `local_context`; request data cannot substitute it.
`into_parts` can return the native repository/adapter to the owning application.
The service never stores credentials or endpoints into biological archives.

The production adapter is the existing
`HttpRemoteTransport<Credentials<OsVault>>`. Its factory derives the local vault
scope from the injected repository, not the expected remote lab locator. It uses
the existing strict HTTP origin/TLS/loopback options, expected-lab headers and sync
envelope pinning, bounded binary responses, disabled redirects/proxies/retries,
and sanitized failures. See [native credentials](DESKTOP_CREDENTIALS.md) and
[offline repository](OFFLINE_REPOSITORY.md). No live vault or network action is
performed merely by constructing credentials or listing local remote metadata.

`SyncAdapter` exposes token-free credential status and consuming store/delete/
session-clear operations in addition to `RemoteTransport`. The injectable backend
is native code, not a serialized behavior or a credential supplied by a remote
server. Tokens appear only in the bounded SetCredential request and the native
credential provider. They are never success/failure fields, outbox metadata,
repository configuration, native payloads or Debug output of the secret wrapper.
Do not log raw request frames; they can contain a submitted token.

## Controls and credential behavior

All current DSQ commands are implemented:

- List/add configured remotes, with native URL/identity validation.
- Set a session token, explicitly remember it in the native vault, inspect its
  token-free status, delete the exact scoped entry, or clear one/all sessions.
- Read-only remote folder browse, name search and immutable revision listings.
- Clone a remote document into an explicit local catalog folder.
- Queue a saved local document, resume its durable push, page pending metadata,
  cancel an unattempted queue, or discard a known rejected operation.
- Pull with clean event-proven fast-forward or an explicit native checkout.
- Begin/resume/retry selected tree jobs, page their status, and queue job-scoped
  folder/document metadata changes that preserve root-mount placement.

Session mode does not invoke the OS vault. Remember is explicit and fails closed
if the vault is locked/unavailable; it never silently installs a session fallback.
A failed deletion does not claim that the persisted credential disappeared.
Credentials remain scoped to local repository, configured remote ID, canonical
origin and expected remote lab. ClearSession does not delete remembered entries.

The service emits no arbitrary base-library write and accepts no archive upload.
All remote content writes pass through the repository's durable immutable outbox.
Read-only wire requests are restricted to Browse/Search/ListRevisions, and each
response is validated against its complete request, correlation and expected lab.
Typed read failures retain their validated remote context and conflict metadata;
free error text is replaced by a bounded native public description.

Clone/pull obtain a validated remote context and use the repository's own exact
payload loading, context checks and graph-preserving import. Their replies carry
the actual remote context and local saved record, never the archive bytes. No
metadata response confers physical-material authority or replaces fresh native
saved-origin replay.

## Saved-head and checkpoint guards

QueuePush/Pull require an explicit **local library revision**, distinct from a
molecule revision or a remote history revision. The service reloads the native
repository, captures its checkpoint internally, and compares the requested local
record with its current saved head. A stale window receives LocalConflict with
the current metadata before any acceptance or upload.

`queue_push_pinned` and `pull_pinned` recheck that captured checkpoint and saved
head before payload selection. The queue records the exact saved source revision;
it never silently reads a later head between a frontend guard and construction.
The existing atomic checkpoint CAS rejects an intervening canonical writer before
committing local queue/import changes. The checkpoint pin is native implementation
state, not a blanket frontend guard or an authorization claim in DSQ.

Clean pull requires the current native source pin to equal the binding's last
common pin and then proves native event fast-forward ancestry. Private branches
are refused before requesting a remote context. Explicit checkout still requires
the exact saved local library revision and persists the displaced native branch
bookmark atomically with the new saved body. Eventless `.cg2` bodies require the
separate explicit revision-pinned snapshot checkout; no event ancestry is invented.
Unrelated local drafts, UI state, private dependency heads and `FromRun` evidence
are preserved by the native import path.

The local design must be saved to its catalog record before queue/pull. This
contract intentionally has no request field for a biological archive, generic
scientific plan, private draft bytes or a replacement remote base revision.

## Durable outcomes and cancellation

The shell owns request/window admission and operation-state tracking. The supplied
accept callback must mark a mutation irrevocably accepted before returning true.
A false return stops the operation before its durable/native-vault/network mutation
step. Read-only metadata calls do not use it. A later canceled window must not
relabel an accepted or uncertain upload result as an unattempted cancellation.

For ResumePush, acceptance precedes the repository's Lookup/Apply sequence. A lost
response can therefore return Unavailable while the durable queue remains marked
attempted. A later resume looks up the same mutation ID/full intent and recovers
its exact receipt without duplicate import. The service does not throw away an
unknown result, advance the remote base blindly, or recreate it under a new ID.

CancelUnattempted and DiscardRejected use the repository's strict existing rules:
only a queue not yet attempted can be canceled, and only an acknowledged rejected
receipt can be discarded. Rejected receipts remain inspectable. The service does
not provide an unsafe “forget unknown upload” operation.

## Selected tree jobs and bounded metadata

Tree commands use durable native job records. Resume performs at most four native
steps per request, even when the caller gives a larger maximum budget; the field
is an upper bound, not a promise to finish a whole subtree. The caller explicitly
polls/resumes later. Catalog snapshot changes, limits, local parent-generation
changes and preserved existing bindings remain in the returned status.

Metadata queue commands call `queue_tree_folder_metadata` and
`queue_tree_document_metadata`, not the ordinary unscoped move logic. A selected
remote root cloned under a private local mount does not acquire an invented remote
parent. An explicit rename can be queued while retaining the original remote
placement. Parent-placement review/acknowledgment remains the repository's separate
native control; this IPC increment does not silently accept a changed local mount.

Pending and tree pages contain at most 200 items, with an explicit next offset.
The host projects only summaries: no full tree node graph, archive frame, native
vault account or credential appears in them. Corrupt typed outbox identities fail
rather than panicking or being silently skipped. Requests are bounded to64 KiB;
complete replies to3 MiB, including framing. Unavailable operations never cause
unbounded metadata or a partial success presented as complete.

## Native fixtures and remaining integration

The service fixtures use real SQLite library/repository stores, the native host
and reliable sync dispatchers, and an injected in-memory vault. They cover:

- Offline metadata reads; session-only credentials; locked Remember failing
  without fallback; token-free responses and wrong local credential scope.
- Stale saved-head rejection before acceptance; pre-accept cancellation; a remote
  import committed before its response was lost; restart and exact receipt lookup
  without duplicate creation or unknown-outcome deletion.
- Pending pagination with a continuation cursor and full remote read correlation/
  context checks.
- Real clone and divergent private native branches, explicit pull with durable
  branch bookmarks, and preserved native saved PCR history with fresh replay.
- Selected tree cloning into a private mount, followed by job-scoped metadata
  queueing that leaves the remote root's parent unchanged.

No fixture opens a real HTTP listener, sends an application HTTP request, invokes
an OS vault, or runs a GUI/browser smoke check. The native shell still needs to
route DSQ frames through its existing trusted `collagene_rpc` boundary and preserve
accepted-operation outcomes. Frontend controls/presentation and OS credential
interaction policy belong to those adapters, not to a path or permission claim
inside a frame.


## Remote revision discussion

The native service dispatches the current bounded `Discussion` command against
one saved, fully pinned remote using its existing native credential provider.
It checks window admission before transport, exact nested request correlation,
remote/lab context, actor and revision/generation bindings. No local plasmid,
scientific history, durable outbox or credential record is changed by discussion.
Transport interruption is not automatically retried. The shared browser/Tauri
review controls handle an uncertain result by requiring refresh and review.
See [discussion workflow](LAB_COLLABORATION.md).
