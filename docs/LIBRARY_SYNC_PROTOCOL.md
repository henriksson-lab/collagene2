# Reliable mutation channel for offline library synchronization

`collagene-library-sync` adds a separate retry-safe mutation namespace. It does
not change the existing version-1 `/api/library` requests or pretend those
ordinary create/save operations are idempotent. The new contract uses
`/api/library/sync` with `application/vnd.collagene.library-sync+postcard`.
Route registration and deployment belong to the native host; adding these files
does not enable an endpoint in the current frozen public assets/server build.

Pure wire types and intent hashing work in Rust/WASM. The `native` adapter and its
SQLite/host dependencies compile only for native targets. No external `git`,
scientific executable, filesystem path, SQL, actor claim or private checkpoint
operation is carried by this wire format.

## Authentication and namespace

The existing trusted host authenticates current membership and selects the lab
before polling/buffering a body or dispatching work. A shared lab access key or
other established authentication supplies that scope; an ID on the wire grants
no permission. Responses carry a host-derived ordinary `LabContext` so clients
can reject another lab's response. Durable mutation IDs are scoped by lab in the
storage ledger, and local clients separately retain which remote/lab a pending
operation belongs to. The envelope carries a required `expected_lab` database
pin. Native dispatch checks it against the independently authenticated `LabAccess`
before receipt lookup, native validation or mutation. The pin is a precondition,
never a membership grant. A server reconfigured at the same origin therefore
rejects an old client's Import before it can commit into the wrong lab.

The new frame keeps the familiar 42-byte structure but uses separate magic:

| Bytes | Meaning |
| --- | --- |
| 0–3 | `CSQ\0` request or `CSR\0` response |
| 4–5 | little-endian exact sync schema, **1** |
| 6–9 | little-endian payload byte length |
| 10–41 | BLAKE3 of the postcard payload |
| 42 onward | postcard bytes, with no trailing payload |

Requests are bounded at 65 MiB, retaining the ordinary native document 32 MiB /
workspace 64 MiB payload limits. Responses are separately bounded at 256 KiB;
the native durable receipt is capped at 128 KiB. Identifiers, names, native pins,
error text and byte deserializers reuse the ordinary bounded protocol DTOs. Exact
version/magic/length/checksum checks precede deserialization; unknown versions,
wrong-direction frames, excess announced lengths and trailing bytes fail.

## Mutation identity and intent

`Mutation { id, operation }` uses a generated `mut-` plus 32 lowercase hexadecimal
digits. A client generates and durably records a new ID **before** attempting a
new logical operation; a retry uses the same ID and exact operation. The envelope
`request_id` correlates an exchange and is not the durable operation identity.
The final envelope is `{request_id, expected_lab, request}`. `request_prefix`
validates the frame and borrows the exchange ID/lab pin without decoding or
allocating its archive. HTTP adapters additionally bind
`X-Collagene-Expected-Lab` before body polling or credential/executor access.

Eight public operations mirror the ordinary mutation fields:

- `CreateFolder`, `RenameFolder`, `MoveFolder`, `DeleteEmptyFolder`.
- `ImportDocument`, `SaveRevision`, `RenameDocument`, `MoveDocument`.

Rename/move/delete retain their expected metadata generation; save retains its
exact base library revision. There is no force-write or private sidecar variant.

`Mutation::intent_digest()` and the native storage implementation independently
derive the same canonical digest:

1. Serialize an externally tagged `IntentOperation` with the same variant and
   field names/order as these eight operations.
2. Replace a native payload with the ordered descriptor
   `{kind, archive_hash, archive_size, selection, sharing}`. The archive hash is
   checked against actual bytes before dispatch. Selection fields remain
   `{document_id, revision, head_event}`; absent fields serialize as `null`.
3. Hash domain bytes `collagene-library-mutation-v1\0`, then those UTF-8 JSON bytes
   directly, using BLAKE3. There are no intervening length prefixes.

The nonce, exchange ID, authentication, remote URL, local checkpoint and raw DNA
byte array are excluded from the intent. Native bytes are still cryptographically
bound by their validated hash, size, kind and selected native pins/sharing policy.
This avoids expanding a 64 MiB archive into a JSON integer array. Field order and
the domain are explicit versioned contracts; changing them requires deliberate
compatibility work and the wire/native cross-fixture.

Apply carries only the typed operation and mutation ID; its intent is computed,
not accepted as a caller's assertion. A reused lab/ID with different intent yields
`IntentMismatch` and performs no new domain effect.

## Durable receipts and lost replies

Requests are `Apply { mutation }` or
`Lookup { mutation_id, expected_intent }`.

A durable receipt is
`{id, intent, catalog_generation, outcome}`. Its outcome is either:

- `Applied`, containing the exact returned created/renamed/moved record, deleted
  empty-folder ID, imported document, or saved document+revision.
- `Rejected`, containing a bounded domain error and, for CAS conflicts, the exact
  folder/document snapshot observed by that original attempt.

The SQLite extension commits the public effect, receipt and any trusted native
repository callback/checkpoints in the same transaction. On failure before that
commit, no successful effect or receipt can become half-visible. A matching
retry returns the exact previously committed receipt, without generating another
document/folder/revision or substituting a newer conflict snapshot. A replay does
not rerun a trusted native checkpoint callback that already committed.

The native adapter's wire dispatch exposes **no** checkpoint operations or local
repository administration. Those callbacks are trusted in-process facilities for
the repository implementation, outside the remote protocol and intent.

Lookup returns the exact receipt, `Unknown` bound to the requested ID/intent, or
`IntentMismatch`. Unknown is not proof that an earlier HTTP job never started:
resending the **same** Apply remains safe because its atomic ledger check decides
which attempt commits. The client must not generate another nonce merely because
a response was lost.

Recommended recovery sequence:

1. Persist the exact pending Apply bytes, nonce, intent, local publication
   selection/source and expected remote base in the offline repository.
2. After uncertain delivery or restart, Lookup that same nonce and intent.
3. If Unknown, resend those same Apply bytes. If a receipt exists, validate it
   against the full pending mutation before acknowledging it locally.
4. Atomically update the appropriate per-remote binding and outbox state only
   after a matching applied receipt. Keep local edits on conflict or uncertainty.

`decode_response_for` binds the exchange and Lookup's ID/intent.
`Receipt::validate_for(&pending_mutation)` additionally checks the complete
operation, returned identities, destination/name, saved parent and payload
descriptor. A repository using Lookup must perform that second check before
acknowledging its outbox. `validate_for_in_context` also binds the authorized lab.

## Native data and failure boundaries

The native adapter first validates bounded fields and raw payload integrity. A
matching committed receipt can be replayed immediately as historical operation
evidence. Before a first native import/save it uses the existing host's native
document/archive loaders and independently derived selected-document closure.
Whole-workspace sharing cannot be granted by the enum alone. Native IDs, event
pins, source ancestry and all existing physical replay rules remain unchanged.

Transport/malformed-native/publication-scope validation errors occur before a
domain mutation and return bounded non-durable `Failure`, not an invented durable
success. Storage/domain results such as stale CAS, duplicate name or missing item
are durable `Rejected` receipts. A storage outage, ledger quota or aborted native
transaction must retain the pending operation for recovery; no fallback through
the non-idempotent v1 write endpoint is permitted.

Receipts do not establish current scientific correctness or reusable material.
Opening a saved archive still uses the ordinary native loaders, and downstream
cloning still freshly verifies/reconstructs its exact source/event ancestry.

## Scope and validation status

This is a reliable mutation channel, not a complete automatic merge algorithm.
Offline repository state, folders, multiple configurable remotes, unpublished
heads and explicit divergence handling remain separate typed repository work.
No remote credentials are serialized into plasmid archives or mutation receipts.

Focused source fixtures cover malformed frames/large length prefixes, namespace
separation, descriptor/nonce intent binding and native/wire digest parity. Native
workflow fixtures simulate a lost import reply, intervening metadata edits,
receipt lookup/replay without duplicate creation, lost save without duplicate
revision, durable stale-conflict replay after later writes, cross-lab isolation
and native-invalid payload rejection without a domain effect.

No builds, root registration, existing route changes or root lockfile changes were
performed while the Trunk publication snapshot was frozen. The parent coordinates
registration and combined native/WASM validation after that build completes.
