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

The [current full inventory](refactor-layout-current.json) covers **1,247 Rust
files / 215 directories** with **zero violations and zero size warnings**.
The largest directory has **17** direct Rust files; the largest file is
`web/src/app/mod.rs`, **2,869 lines**. This includes uncommitted citation work.
The existing `/home/mahogny/github/claude/rust_refactor` release binary's read-only
`module-size-check --include-generated` independently reports **664 compliant
production files**. It excludes non-production files and supplements the
repository-wide gate. No helper build, automatic rewrite, compiler run or
semantic scatter analysis was performed. Ownership findings are source review,
not exhaustive semantic proof or scientific verification.

## Ordered implementation batches

| Priority | Concrete change and evidence | Required boundaries |
| --- | --- | --- |
| 1 | **Reduce App responsibilities before further growth.** Move `run_plan` at `web/src/app/mod.rs:1283` through its complete asynchronous lifecycle into existing `protocol/controller.rs` owners. Move enzyme setup persistence (`:1797`), retained digest operations (`:1866`), antibody operations (`:1973`), RNA persistence (`:2098`) and comparison persistence/admission (`:2159`) into their existing domain controllers. | Move request capture, reset, cancellation/freshness checks and result publication together. Preserve generation, workspace epoch, exact target and all source-event guards. Keep App as composition; no universal callback container or scattered App impls. Its small material-modal composition wrapper can remain. |
| 2 | **Preparation declaration methods, with the next preparation batch.** Move `preparation_catalog_selection` and `preparation_operation_count` at `core/src/cloning/assembly/recipe.rs:77` and `:101` onto `Preparation` in that same file. | Update combinatorial, ligation/tailing, web and fixture callers. Preserve counts, exact catalog pins and shallow saved-origin semantics. Keep multi-input compilation/preflight as algorithms. |
| 3 | **Tracked-material methods, with the next tracked-material batch.** Move `summary:869` and `source_preflight:1051` in `core/src/chemistry/tracked_duplex.rs` onto `TrackedDuplex`; review private `strand:304` and `index_at:1058` at the same time. | Keep the cohesive 1,651-line owner together. Update assembly, edited-primer, overlap-extension, domestication and fixture callers. Preserve validation, positional chemistry and exact hash bytes. Multioperand ligation and homology assembly remain free algorithms. |
| 4 | **Renderer State methods, with related viewer work.** Convert private `reference_position:1220`, `advance_gesture:1254`, `feature_spans:1296`, `nearest_span:1319` and `draw:1497` in `web/src/viewer/map/renderer.rs` where they naturally own/query State. | Keep State, GPU resources and cleanup together. Preserve rotation, pointer behavior, selection previews and frame scheduling. Pure geometry, tessellation and colors remain free functions. No arbitrary renderer split. |
| 5 | **Contig identity methods, with the next contig change.** Move `options_revision:223` and `report_revision:1902` in `core/src/review/contigs/mod.rs` to `ContigOptions` and `ContigReport`. | Preserve exact hash inputs/bytes and freshness. Keep multi-input batch identity free; no generic identity abstraction. This does not advance deferred BAM/BAI work. |
| 6 | **Group growing domains before adding more siblings.** Start with the 17-file core sequence-format suite; review the 16-file core read, material-derivation, annotation and protein suites, MCP analysis/cloning and sharing tests, and web workspace/preparation/method directories with their next relevant addition. | Group actual shared concepts into registered children, preserving complete owners. No speculative empty tree. Update Cargo/module registration and compare executable test names/counts before and after moves. The ceiling is not a target. |

Independent web and core ownership batches can proceed in parallel. One integrator owns shared
callers, registration and verification. Narrow broad parent imports during
substantive edits when the dependency boundary becomes clearer; do not create
an isolated cosmetic batch for imports.

## Completed work and justified owners

Do not reopen proposal/evidence and prepared-pool extraction, Draft consolidation
and natural methods, intended-reference accessors, removed import facades,
combinatorial methods, TreeJob methods, worker domain routing, restoration owners,
MCP root/test registration or the model annotation projection/qualifier/INSDC
owner. The existing [implementation record](ownership-implementation-2026-11-28.md)
and [native editing contract](../NATIVE_ANNOTATION_EDITING.md) retain their
verification evidence. The citation batch now passes native/WASM checks and 427 focused cases, including
14 new fixtures; [discovery and execution](citation-integration.json) remain
separate evidence. Its complete edit-log owner is 733 lines and preserves the
current encoding. Exact-source CI remains required. A potential admission versus
projection performance cleanup needs evidence of repeated work before changing
this owner; do not create another citation representation.

Native dispatch is also complete: `NativeState` now owns catalog backup, sync
and source dispatch in the same desktop main file, with direct command callers
and no free-function forwarding entry points. The desktop compiler check passes.
Static operation-body comparison found only receiver/formatting changes. Catalog
guard lifetimes, independent operation lanes, nonce binding, cancellation and
durable receipts remain in those methods. Cross-platform desktop CI remains
required for the current source.

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
