# Current code organization plan

28 November 2026. Three parallel reviewers refreshed this queue against the
current working tree: web/worker, scientific core/model, and native/MCP/repository.
This is the single current refactor queue. Earlier audit records document
completed changes; they are not additional unfinished plans. This audit changes
organization priorities, not production behavior.

## Rules and evidence

[AGENTS.md](../../AGENTS.md) records the requirements: at most **20 direct Rust
files per directory**, aiming below that ceiling; at most **5,000 lines per
file**; an ownership review at **3,000 lines** before further growth. Keep a
type and its meaningful implementations together. Prefer inherent methods for
operations that own or query its invariants. Split by genuine domain or
lifecycle boundaries, never arbitrary chunks or a desire to hit a line count.

The [current full inventory](refactor-layout-current.json) covers **1,254 Rust
files / 217 directories** with **zero violations and zero size warnings**.
The largest directory has **17** direct Rust files; the largest file is
`web/src/app/mod.rs`, **2,413 lines** after the current controller/owner source batch.
The existing `/home/mahogny/github/claude/rust_refactor` release binary's read-only
`module-size-check --include-generated` independently reports **667 compliant
production files**. It excludes non-production files and supplements the
repository-wide gate. No helper build, automatic rewrite, compiler run or
semantic scatter analysis was performed for this audit. Ownership findings are source review,
not exhaustive semantic proof or scientific verification.

The refreshed [core audit](refactor-core-current.md), [web audit](refactor-web-current.md)
and [native audit](refactor-native-current.md) record concrete source findings and
justified domain/lifecycle splits. They feed this queue, not separate plans.

## Ordered implementation batches

| Priority | Concrete change and evidence | Required boundaries |
| --- | --- | --- |
| 1 | **Own remote clone review state.** Extract tree pagination, parent review/acknowledgment, invalidation and recovery controls from `web/src/library/remote_library_view.rs` into a cohesive owner. | Retain one command dispatch, exact context pins, host CAS checks, cancellation, and separate acknowledgment/Retry actions. Keep unrelated sync lifecycles separate. |
| 2 | **Contig evidence methods.** In `web/src/reads/contigs.rs`, make `current_contig` and `needs_ack` methods of `BoundAssembly`; let `EvidenceState` own construction, reset and selection instead of callback duplication. | Preserve evidence freshness and uncertainty acknowledgment. Keep complete owners together; no physical file split is required. |
| 3 | **Download preparation constructor, with the next export edit.** Replace `prepare_download` with `ExportArtifact::prepare` and update callers. | Preserve reservation before copying, byte bounds and guarded durable completion. |
| 3 | **TreeJob encoding, with its next persistence edit.** Move `repository/src/folders.rs`'s `encoded(job)` to private `TreeJob::encode` in the same owner file. | Preserve exact encoding, validation and CAS. Keep Repository transaction methods with Repository; no schema bump for a receiver-only change. |
| 1 | **Contig identity methods, with the next contig change.** Move `options_revision:223` and `report_revision:1902` in `core/src/review/contigs/mod.rs` to `ContigOptions` and `ContigReport`. | Preserve exact hash inputs/bytes and freshness. Keep multi-input batch identity free; no generic identity abstraction. This does not advance deferred BAM/BAI work. |
| 2 | **Group growing domains before adding more siblings.** Start with the 17-file core sequence-format suite; review the 16-file core read, material-derivation, annotation and protein suites, MCP analysis/cloning and sharing tests, and web workspace/preparation/method directories with their next relevant addition. | Group actual shared concepts into registered children, preserving complete owners. No speculative empty tree. Update Cargo/module registration and compare executable test names/counts before and after moves. The ceiling is not a target. |

Independent web and core ownership batches can proceed in parallel. One integrator owns shared
callers, registration and verification. Narrow broad parent imports during
substantive edits when the dependency boundary becomes clearer; do not create
an isolated cosmetic batch for imports.

## Completed work and justified owners

The refreshed audit's trace-channel duplication is fixed in the pending delivery
batch: `TraceChannels::all` is the sole borrowed A/C/G/T accessor; core display
and web rendering use it directly. The duplicate free function is deleted. No
wire change or new test is needed for that receiver-only consolidation; the
connected trace fixtures and native/WASM checks cover the batch.


The [controller/scientific-owner batch](controller-ownership-2026-11-28.md) is
implemented in source: App callbacks now belong to existing protocol, restriction,
biology and comparison controllers; Preparation and TrackedDuplex own their
operations; renderer State owns the selected gesture/draw methods. App is 2,413
lines. Static comparison preserves the moved bodies, with the explicit request
supersession/complete PlanState invalidation fixes documented separately.
Native/WASM compiler checks and 107 focused cases pass; all 2,509 baseline case
names in 33 harnesses are preserved. Exact-source CI remains required. Do not
repeat these source moves or treat comparison/native tests as UI runtime acceptance.

Do not reopen proposal/evidence and prepared-pool extraction, Draft consolidation
and natural methods, intended-reference accessors, removed import facades,
combinatorial methods, TreeJob methods, worker domain routing, restoration owners,
MCP root/test registration or the model annotation projection/qualifier/INSDC
owner. The existing [implementation record](ownership-implementation-2026-11-28.md)
and [native editing contract](../NATIVE_ANNOTATION_EDITING.md) retain their
verification evidence. The citation batch now passes native/WASM checks and 427 focused cases, including
14 new fixtures; [discovery and execution](citation-integration.json) remain
separate evidence. Its complete edit-log owner is 733 lines and preserves the
current encoding. Its `82937c1` checkpoint passed scientific and desktop CI. A potential admission versus
projection performance cleanup needs evidence of repeated work before changing
this owner; do not create another citation representation.

Native dispatch is also complete: `NativeState` now owns catalog backup, sync
and source dispatch in the same desktop main file, with direct command callers
and no free-function forwarding entry points. The desktop compiler check passes.
Static operation-body comparison found only receiver/formatting changes. Catalog
guard lifetimes, independent operation lanes, nonce binding, cancellation and
durable receipts remain in those methods. Cross-platform scientific and desktop
CI passed at `82937c1`; current source changes retain their own CI gate.

Keep model Sequence, Location and Molecule definitions with their meaningful
implementations. Annotation editing, qualifiers and INSDC parsing have genuine
responsibility boundaries. Keep serializable Draft separate from Leptos review
lifecycles, CachedWorker with its cleanup, and the small browser worker adapter
over core AnalysisSession.

Keep domain-specific MCP handler impls beside request schemas and framework
registration while session/storage authority stays centralized. Retain Repository
transaction extensions, SyncHost publication/removal lifecycles, sharing's bounded
Graph and complete material/SBOL closure indexes, and desktop backup/activation
owners. TOPO graph execution, verifier Execution and private material replay/install
kernels warrant cohesive ownership. The review supports targeted changes rather
than a wholesale scientific or native rewrite.

## Verification and closure

Implement each coherent owner/caller/schema/fixture batch before compiling.
Delete superseded paths and update callers directly: no compatibility aliases,
reexports, delegators, fallback decoders or scaffolding. Temporary intermediate
breakage is acceptable. Preserve scientific correctness, data integrity, resource
bounds, revision/cancellation guards and atomic publication.

Run formatting, whitespace checks and `python3 scripts/check_rust_layout.py`.
For moved tests, compare executable discovery, not just source attributes. Run
necessary focused native and WASM checks for the finished batch in the background,
using the shared Cargo cache and consistent package/features. Parallel workers
must not compile their individual patches. Do not repeat unchanged passing cases
without a concrete reason. Do not let automatic release verification in the
refactor tool override this workflow; use read-only reports or reviewed edits
followed by the repository's checks.

No browser, screenshot or application HTTP checks until all in-scope code is
complete. Cross-platform CI and product acceptance remain separate delivery gates.
Passing this audit does not complete the product goal.
