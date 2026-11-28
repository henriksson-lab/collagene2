# Scientific core ownership audit

Static review of the working tree on 28 November 2026. Existing uncommitted
delivery work was inspected without modification. Findings concern ownership
and organization, not a claim of scientific incorrectness.

## Measured limits

`python3 scripts/check_rust_layout.py` passed: **1,254 Rust files in 217
directories, zero violations and zero size warnings**. Core has 214 source and
374 integration-test files; model has seven source files. Largest core files:
`review/contigs/mod.rs` 2,252 lines; `cloning/assembly/recipe.rs` 2,245;
`planning/combinatorial/mod.rs` 2,105; `verification/engine.rs` 2,071. No hard
limit currently forces a scientific-core split.

## Concrete remaining work

| Priority | Evidence | Refactor decision |
| --- | --- | --- |
| Next affected owner batch | `core/src/review/contigs/mod.rs:223` has `options_revision(&ContigOptions)`; `:1902` has `report_revision(&ContigReport)`. Both solely validate/hash their receiver and both types already have local inherent implementations. | Make these `ContigOptions::revision` and `ContigReport::revision` in the same module. Directly rewrite callers and remove free functions. Preserve validation order and revision bytes. Keep mixed-input `batch_revision(reads, options)` free. No wrapper type or new file is justified. |
| Pending trace batch or next trace edit | `core/src/review/reads/traces/mod.rs:42` has private inherent `TraceChannels::all`; `traces/display.rs:110` duplicates it as public free `channels(&TraceChannels)`. Display preparation/validation and `web/src/reads/trace_view.rs:185` use the latter. | Keep one inherent borrowed A/C/G/T accessor on `TraceChannels`, with visibility required by the web caller. Directly update consumers and remove the duplicate free function. Preserve channel order and zero-copy slices. |
| Before further sequence-format test growth | `core/tests/formats/sequences/` has 17 direct Rust files. Three SnapGene siblings and two oligo siblings form natural families. | Group complete SnapGene scenarios under registered `snapgene/`, and oligo scenarios under `oligos/`. Retain the existing Cargo harness; directly update fixture paths and compare before/after executable test discovery. Preserve assertions. |
| Before the next related test expansion | Annotation, protein-sequence and material-derivation test directories each have 16 direct Rust files. | Group only the actual growing domain. For annotations, coordinate/native/source-evidence editing is a reasonable cohesive family. Preserve explicit registration and executable inventory; avoid one directory per scenario merely to lower the count. |

Paths in the table are relative to `crates/`. Pending trace line positions may
change before commit. These findings come from inspected source, not a semantic
scatter-tool report.

## Completed findings: do not reopen

The earlier combinatorial method conversion is complete: `Design` owns its
combination identity, construct-shape validation, matrix size and paging;
`SelectedMatrix` owns validation. Stateless Cartesian algorithms remain free.

Model now owns annotation projection and the INSDC location codec under
`model/src/annotations/`, and source-bound citation/edit replay under
`model/src/genbank_edit.rs`. `Molecule` remains with its meaningful editing
implementation. The earlier dependency/ownership finding is resolved; do not
duplicate this parser in core or reopen the extraction.

Pending trace regrouping keeps the complete `Trace` owner/parser together and
places derived display preparation in a separate lifecycle module. Registered
trace-test children reduce the read-analysis parent to 13 direct Rust files.
This is a substantive domain boundary rather than scattered `impl Trace` files.

## Boundaries worth preserving

- The verifier’s private `Execution` state and implementation remain together.
  Source admission, evidence projection and checked output have real authority
  boundaries; do not replace them with a universal shared context.
- Material preview types, computation and validation are colocated in their
  scientific domain modules. Replay and atomic workspace installation remain
  with the parent owner. Avoid DTO-only extraction or forwarding facades.
- Contig assembly has real private `Prepared` and `Layout` owners. Symmetric
  affine-overlap algorithms can remain free. At 3,000 lines, review a complete
  overlap kernel or layout owner before growth; do not scatter report methods.
- Worker request/response/envelope/session/execution boundaries have concrete
  transport and lifecycle reasons. Stateless domain handlers need no catch-all
  context. Preserve one central response resource-admission authority.
- Macro-defined domain types share names such as `RepairPreview`, `Capsule`,
  `Report` and `Options`. Name-only searches falsely imply scattered ownership;
  do not merge types without resolved identity.

## Guidance and verification

Read the local rust-refactor module workflow and ROZT roadmap. Their ownership,
complete-item, resolved-caller and macro/cfg cautions apply. Their write mode’s
release checks, optional compatibility delegators and old-path reexports conflict
with repository rules and were not used. No helper binary was built or run.

Only the static layout gate ran. No Rust compiler/test, browser, screenshot or
application HTTP check ran for this audit. Future source work should change
owners/callers in one batch, preserve scientific/serialization behavior, compare
moved test discovery, then perform one focused native/WASM check batch in the
existing cache. Documentation alone does not warrant another compilation.

The trace accessor finding is now resolved in the connected delivery batch: core
and web use the one inherent `TraceChannels::all`; the duplicate helper is deleted.
See `refactor-current-integration.json` for current check and execution evidence.
