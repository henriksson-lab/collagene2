# Bounded desktop repository IPC

`collagene-desktop-protocol` defines the raw byte contract for one desktop
`collagene_rpc` command. The browser supplies a `Uint8Array`; the native host
returns a checked binary response. Scientific core data continues to use its
existing native document/workspace formats. No JSON integer-array transport,
filesystem path, actor claim, window credential or network command is introduced.

The trusted host supplies window access, the owned repository/profile root and
the local library scope out of band. Bootstrap returns the host-derived local
`LabContext`; supplying an identifier in a request does not grant access. This
crate is portable Rust/serde and contains no Tauri, SQLite, native filesystem,
HTTP or external-program dependency.

## Frame and resource bounds

Frames use a separate namespace from browser/server library and sync traffic:

| Bytes | Meaning |
| --- | --- |
| 0–3 | `DTQ\0` request or `DTR\0` response |
| 4–5 | little-endian exact schema, currently **1** |
| 6–9 | little-endian payload length |
| 10–41 | BLAKE3 of postcard payload bytes |
| 42 onward | postcard payload, with no trailing data |

The request payload limit is 129 MiB, allowing two independently supplied native
64 MiB workspaces plus bounded metadata. Responses are limited to 65 MiB, working
slots to 64 MiB each, and view preferences to 1,024 bytes. Custom bounded byte,
string and collection DTOs are reused from the library protocol. Framing,
version, direction, declared length, checksum and complete consumption are checked
before a typed request is executed. Byte integrity does not establish native
scientific correctness or material authority.

`RequestEnvelope { request_id, operation_id, command }` uses a nonzero exchange ID
and an `op-` plus 32 lowercase-hex-digit operation ID. Responses bind both fields.
The host scopes operation IDs to its actual window capability, never to a
window/actor assertion accepted from the byte body.

## Commands

- `Bootstrap`: local library context, current/previous working metadata and view
  preferences.
- `Library`: existing bounded ordinary local folder/document requests. Results
  preserve typed library successes and domain failures. They are dispatched only
  against the host's owned offline repository, not a client-selected database.
- `WorkingMetadata` and `ReadWorking`: inspect current/previous slots and read a
  slot, optionally pinned to the known working-pair version.
- `SaveCurrent`: save a native workspace with expected pair-version CAS; the
  native storage policy controls current/previous rotation.
- `ReplaceWorking`: supply **both** current and caller-selected previous native
  archives under one expected pair version. These are independent byte bodies;
  the previous body must never silently become another copy of current.
- `RecoverPrevious`: explicitly recover under an exact known pair-version pin.
- `ReadPreferences`/`SavePreferences`: independent bounded preferences version/CAS.
- `Cancel`/`OperationStatus`: window-scoped operation control; no data/profile
  selector or unrestricted host command is carried.

Working metadata is `{version, current, previous}`; each present slot includes
its exact archive hash and byte count. A present slot requires a version. An
initialized empty pair may also have a version, depending on native storage.
Preferences have their own `{version, bytes}` and do not share working-pair CAS.
An uninitialized preference value has no version and empty bytes; explicit
clearing may produce a versioned empty value.

The expected version on a write is a compare-and-swap precondition, not a request
to overwrite whatever is current. Working and preference conflicts carry their
current typed metadata separately. Replies for Save/Replace bind returned slot
hashes/sizes to the supplied bodies. Read replies bind slot and requested version.
Local library replies retain ordinary document/folder IDs, archive descriptors,
revision parents, metadata generations and page snapshots.

## Cancellation and commit linearization

States are `Unknown`, `Queued`, `Running`, `CancelRequested`, `Committing`,
`Completed`, `Cancelled` and `Failed`. The native host owns transitions and the
atomic boundary before a durable write. A cancel arriving after that boundary
cannot truthfully undo the write: its acknowledgment reports Committing or
Completed. `Success::Cancelled` is a cancellation-request acknowledgment whose
**status** must be inspected; it is not a boolean assertion that work was undone.
The protocol rejects an acknowledgment that labels an unchanged Running/Queued
operation as an accepted cancellation.

The host can create a Cancelled tombstone for an unknown operation ID so a queued
job that has not registered yet cannot subsequently begin. Completed IDs are not
permission to rerun mutations. These are in-process window-scoped controls, not
durable remote mutation receipts or scientific verification.

`decode_control` is a deliberately narrow fast path. It checks the request
namespace/version/exact declared size and borrows only the exchange ID, operation
ID and command tag. For normal jobs it returns `None` immediately, without
decoding/allocating a possible 128 MiB archive pair or hashing that entire body.
**None is not normal-request validation**: the ordinary host path must still call
`decode_request`, which verifies the full checksum and payload before execution.

For Cancel/Status tags the complete frame must be at most 1 KiB. Its checksum and
full typed command are verified before `Some(ControlRequest)` is returned. This
lets native control dispatch run before a working-archive host mutex while
preserving all validation required to perform a control action.

## Native and security boundary

Only the actual native host can bind a message to WindowAccess and its configured
profile. Native document/archive loaders still check formats, expansion limits,
exact history/events and source ancestry. Imported metadata, operation status,
hashes and IPC success do not bootstrap reusable physical material. Working data
and caller-provided backups require those same loaders before being committed.

Failures distinguish bounded public Invalid messages, working/preferences
conflicts, cancellation, service unavailability and typed local library failures.
They carry no raw database errors, paths, authentication material or arbitrary
host execution instructions. Preferences are a small opaque view-data payload,
not a replacement for native archives or a credential storage API.

## Source fixtures and integration status

Focused fixtures cover namespace/version/length/checksum rejection, caller-owned
previous bytes independently bound from current, read/version conflicts, isolated
preference CAS, cancellation at Committing, small-control limits, hostile ID
lengths, required full validation after a normal fast-path result, and unrelated
library conflicts/operation IDs failing reply binding.

The desktop host and browser bridge are separate work owned by their respective
adapters. No builds, existing source changes, root manifest or lockfile edits were
performed by this worker; the parent coordinates native/WASM validation and
registration. Adding this protocol alone does not enable a Tauri application or
claim completed offline/remote-sync integration.
