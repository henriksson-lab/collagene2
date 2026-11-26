# Shared lab-library binary protocol

The registered `crates/library-protocol` crate supplies bounded binary messages
for the shared-library server and browser folder pane. It has no dependency on
scientific core, SQLite, the browser or an authentication provider. The public
viewer preview and deployment status are recorded in [the plan](PROJECT_PLAN.md);
this protocol alone does not provide a login or complete offline synchronization.

The companion storage contract is [LAB_LIBRARY_STORAGE.md](LAB_LIBRARY_STORAGE.md).
Both interfaces describe one shared folder tree per lab. All authenticated lab
members use that tree; projects and per-folder member permissions are not required.
Inventory, ELN and organizational workflow reporting remain outside this scope.

## Host authority

Requests contain no member, lab, filesystem path, SQL or authentication assertion.
The host authenticates the connection/session, verifies current lab membership,
and obtains the storage `LabAccess` before dispatch. A route or selected lab ID
is a locator; its presence never grants membership. A caller must not invoke
`trusted_host_access` merely because a request supplied a plausible lab ID.

Authorized responses carry a host-derived `LabContext { lab_id, name }`.
Unauthenticated/forbidden responses may omit it. Clients bind replies to the
request ID and active host context using `validate_for_in_context`; they cancel
in-flight operations and clear scoped caches when changing labs or signing out.
Checksums and opaque identifiers do not authenticate a host or restore scientific
material. The native server supplies a shared-key authenticated transport without ambient
cookies; this crate itself provides neither a login nor an exposed server.

## Framing and bounds

The suggested future endpoint accepts
`application/vnd.collagene.library+postcard`. Each message contains:

| Bytes | Field |
| --- | --- |
| 0–3 | `CLQ\0` request or `CLR\0` response magic |
| 4–5 | little-endian schema version, currently **1** |
| 6–9 | little-endian payload length |
| 10–41 | BLAKE3 of the complete postcard payload |
| 42 onward | postcard payload, with no trailing bytes |

The schema version is exact: future or older versions are rejected before
deserialization. This format has a separate namespace from worker messages,
native `CG2\0` documents and workspace archives. Enum order and field order are
wire semantics; changing them requires a protocol version change and fixtures.

The payload cap is 65 MiB, accommodating one 64 MiB workspace plus bounded
metadata. Native documents are capped at 32 MiB; workspace archives at 64 MiB.
Pages contain at most 200 entries. Names have at most 1,024 UTF-8 bytes and 256
Unicode scalar values; queries 256 bytes; native selection IDs 128 bytes; server
error messages 1,024 bytes. Custom string/blob/row visitors check binary lengths
before cloning or reserving owned data, including oversized announced row counts.
Encoding uses a bounded serialization flavor. Frame length, version and checksum
are checked before typed decoding, and trailing bytes are rejected.

These are per-message bounds, not an aggregate server-memory promise. The future
host must bound body reception, concurrent requests and native expanded archive
work. Buffer ownership can otherwise temporarily retain both a frame and decoded
archive. Authenticate early and never read an unbounded body before enforcing the
frame cap.

## Folder and document operations

The virtual root is `None`, not a client-selected disk path. Host-generated IDs
are `fld-` or `doc-` followed by 32 lowercase hex digits; lab context uses the
same form with `lab-`. IDs remain stable across rename/move.

Requests support:

- Browse a folder; create, rename, move or delete an **empty** folder.
- Import a new plasmid document; open its head or a pinned immutable revision.
- Save a new revision with an expected current head; rename/move document metadata.
- Search names across all folders or a declared root/subtree scope.
- List a document's immutable revisions, newest first.

Display names are virtual labels, not path components. Validation matches storage:
nonblank, no control characters, at most 256 scalar values/1,024 bytes. A displayed
slash does not traverse the filesystem. The host resolves each opaque ID inside
the authenticated lab; another lab's ID yields no access. Folder existence,
cycles/depth, duplicate names, transaction isolation and authorization are storage/
host checks, not conclusions inferred from a syntactically valid request.

Browse and search pages use `(name, id)` byte ordering, matching SQLite BINARY
collation. `PageRequest` carries `offset`, `limit` and optional
`expected_catalog_generation`. Every continuation (`offset > 0`) requires that
generation. Any mutation changes the lab catalog generation; stale pagination
returns a typed catalog conflict and the client restarts its listing. Generation
zero is valid for an empty new lab. Revision pages use descending sequential
revision numbers and retained parent identities. A response is checked against
the requested page size, continuation offset, folder/document IDs and snapshot.

## Revision and conflict semantics

Folder/document metadata has a generation starting at one. Rename/move/delete
requests carry its expected generation; stale changes return the current typed
record as a conflict. Saving content is a separate operation and does not rename
or move the document. It carries `base_revision`; storage atomically creates an
immutable child and advances the head only if that base is still current.

A revision ID is a server-derived history identity bound to the document, parent,
archive hash and payload metadata. It is distinct from the raw archive BLAKE3 and
the selected native molecule's revision. The protocol uses distinct types for
history revisions and archive/content digests. Revision metadata also carries its
number, parent, archive size, kind, sharing declaration, native selection and
creation time. No endpoint destroys document content or immutable revisions.

`Conflict` responses distinguish catalog changes, folder metadata changes and a
document's current state. `AlreadyExists` and `NotEmpty` are ordinary typed errors.
Clients preserve their unsaved work and offer reopen/review/save-copy handling;
they must not automatically overwrite the newer head. Requests have a correlation
ID, **not** a durable idempotency key: automatically retrying an import/create after
a lost response can create another entity. A future host can add a transactionally
persisted retry policy under a separately versioned contract.

## Native payloads and scientific ancestry

`NativeDocument` contains an existing native `CG2\0` document envelope and has no
workspace selection or sharing declaration. It preserves the molecule, but cannot
by itself supply absent workspace history, saved run capsules or event ancestry.

`WorkspaceArchive` contains the existing validated workspace archive, plus a
selected native document ID, molecule revision and optional native head event.
The host checks those pins against the loaded archive. If that document has a
head event, the supplied event must match; absence is only acceptable for an
actually eventless supported archive. Imported metadata cannot bootstrap tracked
DNA or alter the existing physical replay rules.

Workspace publication additionally declares **DocumentClosure** or
**ExplicitWholeWorkspace**. The host must construct/validate the dependency
closure needed for the selected document and its retained ancestry, or receive an
explicit whole-workspace sharing choice. A client-supplied enum alone is not proof
that unrelated local documents were excluded or that wider sharing was approved.
Do not silently upload an entire local workspace when a user saves one plasmid.

Before a write, the future host must:

1. Authenticate membership and resolve the destination inside that lab.
2. Bound/check the frame and typed request.
3. Load native data through core's document/workspace loaders, checking their own
   versions, hashes, archive expansion bounds, selection pins and sharing scope.
4. Pass the validated opaque bytes/descriptor to storage's atomic mutation.
5. Map only bounded, public error details into the response; do not send SQL,
   filesystem paths, session secrets or raw internal errors.

Opening reverses this boundary: verify stored bytes/metadata, deliver a bound
revision and payload, then use the ordinary native loader on the client. A saved
scientific report still needs its existing verification/replay rules before it
can supply cloning material. Transport success says nothing about whether a
cloning design is scientifically correct.

## Integration and validation status

The future Rust host adapter maps these storage-independent DTOs to
`collagene-lab-library` under a trusted `LabAccess`. A Leptos server-library pane
can browse the folder tree, open/download native records and retain `{document,
base_revision}` for explicit saves. Local workspaces remain the editing/history
model; database IDs locate shared records and do not replace native event pins.

The unregistered crate includes focused fixtures for framing/corruption,
malicious length prefixes, identity parsing, page/context binding, CAS ancestry,
conflicts and workspace sharing/payload integrity. All **six groups passed** in
an authorized isolated offline harness, using one build job and the storage
worker's shared temporary target after its tests completed. The harness points
at these source/test files and leaves the root manifest/lock untouched. No network
or browser/application check was run. Workspace registration and host/UI
integration are explicit remaining work.
