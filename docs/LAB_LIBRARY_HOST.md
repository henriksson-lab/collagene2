# Native lab-library host boundary

`crates/library-host` joins the storage-independent binary library protocol,
embedded SQLite library, native scientific loaders and document-sharing closure
rules. It is native Rust; it starts no listener, launches no executables and
implements no login or identity provider. See [storage](LAB_LIBRARY_STORAGE.md)
and [wire messages](LAB_LIBRARY_PROTOCOL.md) for their independent contracts.

This is an integration boundary. The registered [native server](SERVER_LIBRARY.md)
supplies shared-key authentication before dispatch. Tauri local storage and full
remote synchronization remain separate delivery work. The browser integration
uses the same binary operations and existing native archive loaders.

## Authority and dispatch

The trusted host calls either:

```rust,ignore
let response = collagene_library_host::dispatch(
    &mut library, &access, &request, publication_scope.as_ref(),
)?;
let bytes = collagene_library_host::dispatch_frame(
    &mut library, &access, &request_bytes, publication_scope.as_ref(),
)?;
```

`access` is an existing native `collagene_lab_library::LabAccess`. The caller must
establish current lab membership before **each** dispatch; the adapter neither
constructs access from wire data nor caches membership. There is no actor, role,
lab override, filesystem path or SQL parameter in this API. Shared folders remain
one lab-wide tree, without introducing ELN, inventory or organizational reporting.

The reply's lab ID/name come from the database under that exact native scope.
All operations, including read, search, immutable revision lookup and typed
conflict recovery, use it. Cross-lab IDs yield no access. A malformed frame has
no trustworthy correlation ID, so `dispatch_frame` returns a transport error
rather than inventing a reply. Valid requests receive protocol-validated replies
bound to their request IDs. Error strings are bounded public messages; SQLite
errors, local paths and native-loader internals are not sent to clients.

Folder/document changes preserve the storage generation CAS. Content changes
preserve the immutable base-revision CAS. A stale save returns the current typed
document and makes no revision; the adapter does not retry or overwrite. After
a metadata conflict, the returned current record is a scoped fresh read and may
already include another intervening change. Clients keep their own unsaved work.
Correlation IDs are not durable idempotency keys: a lost response must not cause
an automatic import/create retry.

## Native byte validation

Every content write first checks the protocol descriptor, byte bound and archive
hash, then uses the existing core loader:

- Native `CG2\0` document: `protocol::load_document`, preserving its exact molecule
  and metadata. This format has no archive selection or sharing declaration.
- Workspace archive: `workspace::load_workspace`, retaining its bounded revisions,
  history, candidates, attachments and saved run evidence. The selected document
  ID must exist; its exact head molecule revision and optional head event must
  match the payload pin. A missing supplied event is accepted only if the loaded
  document really has no event. Display names never identify a document.

Opening stored content repeats its descriptor, byte-integrity and current native
loader/pin checks. It does not demand a new publication scope for a previously
saved immutable archive. A database or checksum alone is not scientific validity.
Native loaders check their own format/history structure; downstream cloning
still needs ordinary fresh verification and exact saved-origin reconstruction.
The adapter does not deserialize a historical duplex summary into live material,
flatten a run to obtain a convenient plasmid, or fabricate missing ancestry.

The binary payload is converted directly between protocol/storage byte containers;
it is never expanded into a JSON array. Small structural metadata uses a bounded
4 MiB conversion followed by protocol response validation. The complete request,
native loader state and storage payload can coexist temporarily. A real host
must bound concurrent requests and received bodies, not assume the per-frame
65 MiB cap is an aggregate memory guarantee.

## Ordinary Save and publication scope

A wire `DocumentClosure` enum is a declaration, not authorization to expose the
whole local workspace. Ordinary Save should derive the selected document's exact
dependency closure automatically, without repeated user confirmations:

1. Construct the archive using `collagene_library_sharing::document_closure` for
   the exact selected native document/revision/event.
2. Set the payload selection to those exact pins and its sharing descriptor to
   `DocumentClosure`.
3. On the trusted native host call
   `document_publication_scope(&access, &payload)`. This loads the archive and
   independently derives the required closure using
   `collagene_library_sharing::inspect_document_closure`.
4. The adapter compares every retained document/revision/event/candidate/attachment
   ID with that derived closure, and checks every dependency document's projected
   head revision/event. Unrelated documents, private sidecars or extra history
   cannot be admitted merely by relabeling a full archive as one plasmid.
5. Pass the resulting native `PublicationScope` to dispatch. The write checks it
   against the exact lab, archive hash, selection, sharing and entity inventory
   again before calling storage CAS.

The scope has private fields and no serialization/deserialization implementation.
It is trusted publication validation, not a user approval workflow or a scientific
receipt. Native documents return `None` from `document_publication_scope` after
native loading, since they carry no workspace data to authorize.

A lower-level `PublicationScope::document_closure` accepts a native host's
**independently permitted** `ArchiveInventory` for another trusted integration.
It must never receive an inventory simply echoed from a client, or an inspection
inventory treated as proof that all included data should be shared. The automatic
helper above supplies the concrete server-owned dependency derivation.

Whole-workspace publication remains a separate explicit action. The automatic
ordinary-Save helper rejects `ExplicitWholeWorkspace`. A future authenticated
host interaction may call `PublicationScope::explicit_whole_workspace` only after
establishing that explicit intent for the inspected exact bytes. The request enum
alone cannot grant it. Dispatch without a matching scope returns typed Forbidden;
it never silently substitutes a whole workspace or a sequence-only snapshot.

## Offline work and later synchronization

The adapter/library can be embedded by a native Tauri host; a remote connection
is not required for local plasmid analysis or scientific verification. Local
folders and native histories must remain usable offline. LabAccess is a trusted
local-host scope in that deployment, not a requirement for a remote identity.

The present immutable import/save CAS operations are groundwork for synchronization,
not a complete git-like sync protocol. Chosen remotes, ancestor negotiation,
merge/import of fetched dependency graphs and unresolved divergent heads remain
separate work. A fetched archive must not replace unrelated local documents or
unpublished heads. Native document IDs, exact event identities, parent ancestry
and required run evidence remain distinct from database document/revision IDs.
A future sync layer must preserve both branches on conflict and never flatten
scientific dependencies or silently overwrite a newer remote/local head.

## Focused validation

`tests/boundary.rs` covers real native documents over binary request/reply frames,
immutable CAS conflicts and pinned old revision opens; folder metadata conflicts,
search and lab isolation; malformed-but-correctly-hashed native records; wrong
native document/revision/event pins; scoped closure extras and cross-lab replay;
automatic selected-document closure publication versus explicit whole sharing;
and real saved PCR ancestry through shared storage and fresh native reconstruction.
The same saved molecule flattened into a bare native snapshot cannot reconstruct
its missing original run. Malformed frames and external actor claims are rejected
before storage changes.

No browser or application HTTP checks are part of this boundary increment. Native
checks are coordinated after source freeze with the existing scientific feature
graph, rather than compiling a second large core dependency graph during GUI
publication. Before public shared-library deployment, the real host still needs
authenticated sessions/transport, membership revocation behavior, bounded body and
concurrency handling, persistent deployment configuration, and the shared-folder
UI wired to these validated exact-source operations.
