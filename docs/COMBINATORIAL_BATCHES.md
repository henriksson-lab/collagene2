# Selected-construct batch evaluation and historical records

The existing assembly dock evaluates explicitly selected constructs with progress,
cancellation and row-local outcomes. This extends the
[source-bin matrix](COMBINATORIAL_SOURCE_DESIGN.md); it does not enumerate all
Cartesian combinations or create product documents automatically.

## Per-row declarations and native reconstruction

Every ready row retains its own preparation, terminal treatment, method, topology
and literal junction intent. Source-choice geometry remains immutable. The shared
catalog does not replace per-construct manual primers with another row's defaults.
Duplicate source combinations remain separate construct identities.

Unfinished primers/options and unassigned selected sources become explicit
`DraftError` rows. A pending working source edit affects only its active row;
the batch cannot silently compute that row's last saved source instead. Other
valid rows remain eligible. Raw invalid text and pending geometry participate in
the host declaration key, even when two invalid strings produce the same error.
Ordinary row/slot navigation leaves an unchanged matrix key stable.

The current batch route verifies declared prepared-fragment assembly. A row asking
for whole-original/whole-amplicon Golden Gate reaction-pool analysis is explicitly
uncomputed, with directions to its existing reaction inspector. Its raw choices
remain saved. Connecting those existing finite reaction models to batch evaluation
remains work; no kinetic engine is added as a prerequisite.

A private, nonserializable core session owns one validated source workspace and
advances exactly one construct at a time through existing preparation compilers
and the deterministic verifier. The original regions or explicit independent
target event declare the requirement; a simulated product never becomes the target.
Exact historical source events remain usable without rewinding current heads.

Each result distinguishes a declared draft failure, reconstruction error or native
evaluation. Native statuses, completed/planned steps, preparation completeness,
findings, assumptions, target/plan/final-product revisions and source bindings remain
separate. A complete final product can still fail the independent target. An
evaluated or consistent-under-assumptions row is not experimental confirmation.

The session retains bounded summaries, not full duplicated plans/products for every
row. Result rows and evidence notes are fully paged. A current result can open its
construct in the same dock for local correction and explicit rerun. Historical
records cannot silently replace the current editable design.

## Capture, cancellation and replay

Only a completed private session can capture a fresh computation into an immutable
`combinatorial-batch-record-v1` attachment. Capture adds evidence without creating
targets/products, moving source heads or changing history. The record retains the
typed request, ordered row evidence and recorded scientific engine.

The host checks declaration, request, epoch, generation and source-workspace
fingerprint before installation. A save validates and encodes both current and
previous recovery states, then commits them in one cancellable IndexedDB
transaction. Assembly editing controls are disabled during this commit. Cancellation
terminates the same worker used for encoding; an old cancelled completion cannot
unlock a subsequent save because commit leases have distinct operation owners.
An active transaction retains its lease until abort/commit is resolved. A durable
successful commit is restored even if its panel closes just after commit.

Native loading/listing validates record bytes, hashes, shape and retained source
events, but remains historical. Neither JSON nor a saved receipt reconstructs a
private session or restores fresh authority. Explicit replay uses the retained
request and current native engine. The web can replay a historical matrix without
replacing the current editable one. MCP row replay returns recorded and newly
computed evidence separately.

## Limits and transport

Existing limits remain12 bins,256 choices and256 selected rows. A row summary is
at most128 KiB and a record32 MiB; the native workspace retains its aggregate
budgets. Aggregate declared search/repeat work is bounded to2,000,000,000 units,
in addition to each existing scientific-operation limit. This is a conservative
declared-work bound, not a measured CPU-time estimate. Reaching it leaves remaining
ready rows explicitly uncomputed rather than omitting them or claiming completion.
Unrelated valid rich documents keep the native model profile; scientific limits
apply to the actual selected sources and target.

Worker57 appends Begin128/Started120, Advance129/Advanced121 and Capture130/Captured122.
Begin transmits the workspace once; advances carry only context/index and return
small progress/evidence. Begin/Capture and the captured native workspace use explicit
128 MiB native-workspace envelope cases. Advances and progress retain32 MiB.
Pre57 labels reject the new variable-length selectors before payload decode.
Stateless processing cannot fabricate a private batch session. Existing tags,
candidate11, plan9, science14 and native archive1 remain unchanged.

## MCP operations

All operations use explicitly imported immutable `WorkspaceRef` values. No tool
accepts an arbitrary host filesystem path or starts external scientific software.

| Tool | Purpose |
| --- | --- |
| `page_combinatorial_design` | Bounded Cartesian pages with exact source ownership checks. |
| `inspect_combinatorial_construct` | Fresh native reconstruction of one explicit construct against its independent target. |
| `evaluate_combinatorial_batch` | Evaluate selected rows and return a native workspace artifact containing historical batch evidence. |
| `list_combinatorial_batches` | Page historical batch summaries. |
| `open_combinatorial_batch` | Open historical evidence without recomputation. |
| `replay_combinatorial_batch` | Recompute a selected recorded row, separating recorded/current evidence. |

Complete evidence is available through immutable artifacts rather than a truncated
summary. MCP batch evaluation is one bounded request behind the existing compute
gate; it is not a scheduling interface or public persistent job service. These six
operations bring the current registry to142 tools.

## Remaining delivery work

Batch product release now uses a separate fresh-reconstruction transaction and
the existing run/output finalizer. Per-construct output manifests, actual staged
product/target inspection and bulk sequence/primer/HTML/PDF ZIP exports are
implemented; see [contracts](COMBINATORIAL_RELEASE.md). Batch reaction-pool
composition and library topology editing remain in scope. A batch record is not
permission to release imported products or claim a verified library.
Browser/visual/performance acceptance stays deferred until all in-scope code is
implemented.

Validation passed:19 focused core tests (batch6, design8, binary protocol5), two
real MCP SDK fixture groups and all98 native web tests, including four projection
groups and the cancelled-encoder/new-commit lease race. The final WASM compile
passed. Tests cover literal products/wrong targets, local draft failures, per-row
overrides, historical aliases, work limits, capture deduplication, variable-length
selector gates, native preservation and explicit replay. No browser checks or
preview refresh were run for this checkpoint.
