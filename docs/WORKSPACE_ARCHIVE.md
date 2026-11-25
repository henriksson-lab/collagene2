# Native workspace archive, version 1

A `.cg2workspace` is an explicit backup/restoration of an open scientific workspace. It preserves document identities, immutable molecule revisions, retained history branches, exact undo/redo navigation, candidate designs, and original evidence. The existing `.cg2` single-document format remains separate and readable; use that format to exchange one molecule. Loading an archive returns validated data before the application mutates its state. Loading never writes external files or executes history descriptions.

The independently implemented MIT core provides `workspace::{Workspace, WorkspaceDocument, RevisionEvent, WorkspaceAttachment, save_workspace, load_workspace, is_workspace_archive}`. It runs natively and in WASM using the same code. ZIP processing uses pinned `zip` 4.6.1 with Rust-only DEFLATE support; no external archiver is invoked.

## State and history

`Workspace::from_documents(name, molecules, active_index)` creates `document-1`, `document-2`, etc. These are document identities, independent of content hashes, so two tabs with identical content remain distinct. IDs are 1–128 ASCII letters, digits, underscores, hyphens, periods or colons. `Workspace::head(document_id)` resolves a document's current immutable molecule. All revision keys are lowercase BLAKE3 hex digests of the existing core molecule representation.

Each document stores a head revision and ordered undo/redo revision stacks. `head_event`, `undo_events` and `redo_events` identify exact history events: two different edits can return to the same molecule revision. Nonempty event stacks must match their corresponding revision stacks in length, document identity and content. Optional empty event stacks support callers which have content history only; they cannot distinguish repeated-content events by themselves.

History events have a document identity, molecule revision, ordered parent event IDs and a descriptive operation string. Parent edges can cross documents to describe derivation. All references resolve and the event graph must be acyclic. `RevisionEvent::new` constructs the validated event ID as `event-` followed by BLAKE3 of:

1. The literal bytes `collagene-history-event-v1\0` (the final byte is NUL).
2. Compact JSON serialization of the tuple `(document_id, revision, ordered_parents, operation)`.

Descriptions are reviewable history, not typed edit replay. The archive retains every supplied valid event, including inactive branches; history-window and memory policies belong to the application. It does not silently prune events.

Candidate records preserve exact plan revisions, input identities/revisions and ordered operations. All source molecules must match immutable revisions present in the workspace. An imported candidate may reference a target revision which is unavailable; those references are listed explicitly in `unresolved_candidate_targets` and remain unbound until the application resolves them. `active_candidate` is an optional plan-revision marker; retaining an inactive proposal does not reopen it automatically.

Generic attachments are keyed by BLAKE3 of raw bytes and preserve kind, name and media type. This supports original oligo-table evidence stored once, while derived row molecules reference its hash and row identity. Attachment labels and original source documents are evidence, not assertions that they describe a subsequently edited molecule.

## Container and content addressing

The ZIP contains exactly one `manifest.json` and entries named `blobs/<64 lowercase hex digits>`. Every file is declared in the manifest, hash/length checked, and reachable from a recognized record. There are no directory entries or extraction paths. The manifest has format `collagene-workspace`, version `1`, molecule schema `1` and plan schema `9` (with compatible plan schemas 1–8 still readable).

The manifest maps molecule revisions to molecule-record blobs. A molecule record references separate sequence, name, feature JSON, metadata-value and optional original-GenBank chunks. Sequences use canonical ASCII and carry explicit alphabet/topology. A plan record references existing source revisions and its operations JSON; it does not duplicate source molecules inside each plan. Hash-addressed raw attachments, metadata values and original sequence evidence are shared across revisions.

Canonical standard-base64 metadata values whose keys end in `base64` are stored as decoded binary chunks and restored to exactly the same base64 strings. Other values remain exact UTF-8 strings, including unknown metadata and whitespace. This avoids repeatedly storing large original trace evidence as base64. Deduplication applies to identical byte chunks regardless of which records reference them. The restored molecule revision must match the original content digest.

The writer emits Stored entries, a fixed 1980-01-01 timestamp, mode 0600, manifest first and sorted blob names. Identical workspace values produce identical archive bytes. Stored output gives predictable save cost without compression heuristics. The reader also accepts bounded Deflate entries within this profile. ZIP64, encryption, multidisk files, symlinks, comments, prepended executables, trailing bytes, overlapping ranges, hidden local entries and unsupported compression methods are rejected. General-purpose flags are restricted; local headers and descriptors must agree with the central directory. This is a deliberately narrow archive profile, not a general ZIP importer.

The independent structural preflight runs before the ZIP library builds its filename map, so duplicate central-directory names cannot disappear through map replacement. CRC, blob hash, declared length, graph and model checks all complete before the returned workspace can be installed. Duplicate JSON keys and unknown record fields are rejected. Feature and operation JSON must round-trip through the supported typed schema without silently dropping fields. A future molecule/plan representation requires explicit schema migration or a version bump; unknown scientific semantics are not treated as understood.

ZIP layout follows the primary [PKWARE APPNOTE specification](https://pkware.cachefly.net/webdocs/casestudies/APPNOTE.TXT). Container handling uses [zip-rs](https://github.com/zip-rs/zip2); no cloning application's implementation was borrowed.

## Bounds and transport

| Bound | Limit |
|---|---:|
| Archive bytes | 64 MiB |
| Expanded archive entries / materialized logical workspace data | 128 MiB each |
| Single blob or serialized molecule | 32 MiB |
| Manifest | 16 MiB |
| Documents | 1,000 |
| Molecule revisions | 10,000 |
| History events | 20,000 |
| Candidate plans | 128 |
| Attachments | 1,000 |
| ZIP entries | 50,000 |
| Undo or redo stack | 1,000 |

Metadata is limited to 10,000 values per molecule and 200,000 across revisions. Arithmetic is checked; declared ZIP expansion is bounded before decompression, each read is capped, and reconstructed repeated evidence is charged before large duplication. These are logical serialized-data bounds, not a promise that process RSS stays below 128 MiB: validated blobs, materialized objects, JSON parser state and transport buffers can coexist.

Worker protocol version 6 introduced `SaveWorkspace`, `LoadWorkspace` and `ValidateTraces` jobs. Save binds the request revision to the validated active document's head; load binds it to BLAKE3 of the exact archive bytes. Results are `WorkspaceSaved`, `WorkspaceLoaded` and `TraceValidation`. Workspace results are data and carry no imported verification status or permissions.

Ordinary documents and scientific worker messages retain the existing 32 MiB payload cap. Only workspace job/result variants may use 128 MiB payloads, and only with supported worker envelope versions 6–19, subject to each nested operation/report floor. A borrowed prefix preflight rejects oversized non-workspace variants and older envelopes before deserializing their payload. Serialization itself is bounded. Envelope checksums detect accidental/crafted content changes; they provide no authenticity or authorization. Reopened candidates must pass the deterministic verifier again before use as verified designs.

## Validation coverage and limits

Tests include a manually encoded ZIP/reference manifest fixture, exact event navigation with repeated content and branches, candidate preservation, binary-evidence deduplication, Deflate input, malicious paths/duplicate entries/headers/CRC, truncation, graph errors, future schemas and an actual workspace transport round trip exceeding 32 MiB. Native single-document interchange, project synchronization, server persistence, authentication, organizational records and automated external writes are outside this archive layer.

The current worker writer uses version 8, adding feature-library scan/application and annotation exchange while retaining version-6 workspace layouts. The web app stores validated immutable feature-library JSON snapshots as `feature-library-v1` attachments keyed by the BLAKE3 of their exact bytes; complete definitions survive workspace reopening.

Plan schema 2 adds end-repair and oligo-annealing operations; schema-1 archives remain readable when their plans use earlier operations. Labeling the new operations as schema 1 is rejected before candidates can be installed. Candidate binary exchange now writes version 6; imported candidates remain unverified data.


Plan schema 3 introduced tailing and profiled overlap; plan schema 4 adds tracked source-generation declarations and positional-profile checked digests. Native loading rejects those operations under plan schemas 1–3, including operations inside retained run attachments. New retained runs use plan schema 4; old plan-schema-3 records remain historical evidence and retain their original engine identity. Imported modification summaries cannot restore physical state or verification authority. Worker schema 17 gates the new semantics in both inline and archive-bearing workspace jobs/results.


A `genome-reference-link-v1` attachment is a bounded typed external-reference link: role, exact manifest/content hashes, declared build identifiers/counts and an optional opaque cache hint. Loading validates its closed schema but does not load genome DNA, open an arbitrary path, restore a semantic-verification receipt or grant coverage authority. Inline/binary worker workspace jobs/results carrying this recognized link require schema 19. Archive format 1 and plan schema 4 remain unchanged. Large-reference bytes stay outside ordinary molecule/workspace payloads.


Plan schema 8 adds explicit original-source terminal hypotheses. Loading refuses
`DeclareSourceTermini` under older plan schemas, including nested retained-run
records. Candidate 10 and worker 34 carry the same historical operation floor; that
operation was introduced with science-13. Earlier operation layouts and recorded engine
identities remain historical data, and the outer archive stays version 1.

## Assembly draft declarations

The web workspace retains its current declared assembly draft as an
`assembly-draft-v1` JSON attachment under the normal content hash. It contains
named ordered slots, exact source document/event/revision snapshots, one-turn
source intervals, orientation, topology and editable preparation/primer settings.
The web codec caps it at 32 MiB and validates source content and row consistency;
workspace restore additionally resolves each original source event in the native
history. There is one active draft attachment. Existing archive/document/plan
layouts are unchanged.

Frozen intent reports, verification status and private reconstructed chemistry
are not draft fields. Reopening retains declarations only; reviewing and checking
uses fresh worker computation. Created intended-reference documents record the
actual source parents and remain explicitly distinct from simulated material.

Catalog-pinned digest plans write plan schema 9. Loading that operation from a
manifest below 9 fails explicitly, including an archive whose manifest was
relabeled. Existing historical operation layouts remain readable. Shared draft
rows can retain the optional catalog selection and original top/bottom source
intervals; these are declared preparation inputs and geometry receipts, not
restored verification or physical intermediates.

## Finalized cloning outputs

`cloning-output-visibility-v1` records presentation-only hidden preparation document
IDs. All run documents, exact source events and causal parents remain in the
archive and available for replay. Optional saved primers use `oligo-table-v1`.
Finalization does not change archive version 1 and does not restore historical
verification authority. See [output contract](CLONING_OUTPUTS.md).
