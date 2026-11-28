# Web ownership audit

Static review of the working tree on 28 November 2026 by a parallel reviewer.
Web has 297 Rust files, a largest file of 2,413 lines, and at most 16 direct
Rust files per directory. Repository-wide layout: 1,254 files / 217 directories,
zero size/count violations and zero 3,000-line warnings. These are layout and
source-review findings, not scientific or UI runtime verification.

## Bounded findings

1. `library/remote_library_view.rs:123` combines independent sync lifecycles.
   Extract a cohesive clone-review owner for tree pagination, parent review and
   acknowledgment state (152–156), invalidation (174–178), success admission
   (363–376), and recovery controls (665–704). Preserve single command dispatch,
   exact context pins, host CAS validation, cancellation and separate acknowledgment
   and Retry actions. Keep complete owners together rather than fragmenting callbacks.
2. In `reads/contigs.rs`, `BoundAssembly` (57) should own `current_contig` (145)
   and `needs_ack` (148). `EvidenceState` (733) should own construction (743),
   reset and selection; callbacks currently duplicate initialization around
   340–372. Preserve current-evidence uncertainty acknowledgment. No file split
   is required.
3. During the next export change, replace `io/download.rs:17`'s
   `prepare_download` with `ExportArtifact::prepare`. Preserve reservation before
   copying, byte bounds, source guards and durable receipt semantics.

## Completed work and boundaries

Cloning proposal/pool ownership, Draft methods, intended-reference accessors,
renderer State methods and App controller delegation are already implemented.
Do not reopen these findings. App is now 2,413 lines; no size-triggered extraction
is required. Retain domain controllers and worker cleanup with their lifecycle
owners. Stateless and symmetric algorithms need not become methods.

The local rust_refactor module workflow was read. Automatic release checks and
compatibility reexports/delegators conflict with repository rules and were not
used. No builds, tests, browser, screenshots or application HTTP checks ran for
this audit. Source line locations refer to the reviewed uncommitted delivery
batch. Findings belong in the single current REFACTOR_PLAN.md queue.
