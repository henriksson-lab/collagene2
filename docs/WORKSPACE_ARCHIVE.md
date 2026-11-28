# Native workspace archive, version 2

A `.cg2workspace` backs up an open scientific workspace: document identities,
immutable molecule revisions, retained event branches, exact undo/redo navigation,
candidate designs and original evidence. The `.cg2` single-document exchange
format is separate. Loading validates data before the application changes its
state. It does not execute history descriptions or restore scientific authority.

The Rust core owns the same reader and writer for native and WASM builds. ZIP
processing is in memory and never extracts paths or invokes an external archiver.
Only the current internal versions are accepted: workspace **2**, molecule schema
**1** and plan schema **32**. Superseded internal layouts are rejected explicitly.
Intentional external sequence import formats have their own readers.

## Documents and exact events

Document IDs distinguish tabs independently of molecule content. Molecule revision
keys are lowercase BLAKE3 digests of the canonical model representation. Each
native document requires its retained head event. Undo/redo event stacks must
match the revision stacks in length, document ownership and content. Restoration
does not guess missing identities from equal DNA or manufacture recovery events.

A history event contains its document ID, molecule revision, ordered parent event
IDs, descriptive operation label and a required typed transition:

| Transition | Retained meaning |
| --- | --- |
| `IndependentReference` | Opened or created without an edit ancestry claim; no parents |
| `ReferenceEdit` | Exact source event, ordered core edit operations and explicit supplemental source parents |
| `Fork` | An exact source event opened under a distinct document ID with unchanged molecular content |
| `CloningStep` | Native run identity and original step index; fresh scientific reconstruction remains required |
| `Analysis` | Named computational producer and its attachment, native context or explicitly unretained preview |

Reference edits share one operation type between the editor, MCP and history.
They include literal sequence edits, reverse complement, origin rotation,
annotation changes and source-bound annotated paste. Paste retains the actual
source event and selection. The operation label supplies no replay authority.
An analysis context alone does not prove that its producer record is available;
an unretained preview is explicitly not a causal edit recipe.

`RevisionEvent::new` hashes the literal bytes `collagene-history-event-v2\0`
(the last byte is NUL), followed by compact JSON serialization of
`(document_id, revision, ordered_parents, operation, transition)`. Parent references
must resolve, source bindings must match their owning event, and the graph must
be acyclic. Returning to previous content can create a distinct event identity.
Inactive branches remain retained. See [history navigation](HISTORY_GRAPH.md).

Loading checks structural validity and identity. Causal QC freshly reconstructs
supported typed transitions against exact retained sources when requested;
structural archive validity alone does not certify scientific reconstruction.

## Container and content addressing

The ZIP contains one `manifest.json` and entries named
`blobs/<64 lowercase hex digits>`. Every blob has a declared length and digest and
must be reachable from a recognized record. Undeclared files, directory entries,
unsafe paths, duplicate entries and orphan blobs are rejected.

The manifest stores bounded history descriptors, including a `transition_blob`
reference. Each typed transition is a separately bounded content-addressed JSON
blob. Identical payloads deduplicate; expanded logical transition ownership is
still charged for each event when materialized. Large edit literals do not
inflate the manifest itself beyond its separate limit.

Molecule records reference sequence, name, feature JSON, metadata values and
optional original GenBank chunks. Candidate records reference existing immutable
source revisions and operations. Raw attachments preserve their exact bytes,
kind, name and media type. Content-addressed values deduplicate across records.
Canonical standard-base64 metadata values whose keys end in `base64` are stored
as decoded chunks and restored to the same strings; other metadata stays exact
UTF-8. No evidence filename grants permission to read an external file.

Candidates remain proposals. Missing target revisions are explicitly listed as
unresolved and cannot acquire a fresh verification result merely by reopening.
Saved runs, material derivations, read checks and other recognized attachments
retain their current native contracts. Their stored assessments are historical;
fresh native validation is required before using them as scientific authority.

## Resource bounds

| Resource | Maximum |
| --- | --- |
| Archive bytes | 64 MiB |
| Expanded archive entries / materialized logical data | 128 MiB each |
| One blob or serialized molecule | 32 MiB |
| Manifest | 16 MiB |
| Documents | 1,000 |
| Molecule revisions | 10,000 |
| History events | 20,000 |
| Candidate plans | 128 |
| Attachments | 1,000 |
| ZIP entries | 50,000 |
| Undo or redo stack | 1,000 |

Metadata has per-molecule and workspace aggregate limits. Native shape admission
precedes recursive serialization and hashing. Reads, decompression, declared
lengths, repeated evidence and expanded typed actions are bounded with checked
arithmetic. These are logical serialized-data limits, not a process RSS limit:
validated blobs, owned models and transport buffers can coexist.

The editor separately bounds retained document history and edit transactions;
opening a valid archive does not silently discard evidence to fit an edit window.
Binary worker save/load jobs return validated data. Server persistence, filesystem
permissions and desktop synchronization remain separate owners.

This version's source changes and fixtures are part of the current integration
batch; compiler and scientific verification remain pending until that batch is
complete.
