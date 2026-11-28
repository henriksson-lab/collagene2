# Code ownership follow-up

28 November 2026. This is a static follow-up on the current working tree,
including unfinished paired-oligo changes. It supplements the completed domain
refactor; it does not establish scientific or build correctness for those changes.

## Limits and measurements

The repository rules already enforce at most 20 direct Rust files per directory,
at most 5,000 lines per file, and an ownership review at 3,000 lines. Keep a type
with its meaningful implementation; use inherent methods for operations on its
state and invariants. A domain or lifecycle boundary must justify distributing
its implementation. Stateless and symmetric algorithms may remain free functions.

The current layout check covers 1,213 Rust files in 210 directories with zero
violations and zero size warnings. The largest directory has 17 direct Rust
files (`core/tests/formats/sequences`); several have 16. Group the next related
family before approaching the ceiling instead of treating 20 as a target.

| File | Lines | Review concern |
| --- | ---: | --- |
| `core/src/protocol/execution.rs` | 2,990 | One large worker dispatch function; almost at the review threshold. |
| `web/src/app/mod.rs` | 2,945 | Cross-feature composition still includes workspace restoration and preparation lifecycle callbacks. |
| `core/src/material/derivation/mod.rs` | 2,827 | Durable record/capsule authority and several method-specific preview owners share one module. |
| `web/src/cloning/preparation/cloning_methods.rs` | 2,511 | Method forms, request construction and evidence rendering share a panel module. |

Paths in this table are relative to `crates/`. Measurements are a working-tree
snapshot, not a grandfathering baseline or a permanent exception.

## Current action queue

Three parallel reviewers checked core, web, and MCP/native/sharing ownership.
Core has 193 production and 370 test Rust files; web has 290 Rust files. These
are targeted organization reviews plus size/count inventories, not full
behavioral audits. Previously resolved flat-directory and scattered Draft
implementation findings must not be reopened as though still present.

| Order | Concrete change | Required preservation |
| --- | --- | --- |
| 1 | Split `protocol/execution.rs`'s single `process_request` into complete typed domain handlers, starting with workspace/material, cloning and sequence/review. Keep a short exhaustive central router. | Central envelope admission, resource limits and revision checks remain mandatory before execution; retain enum selectors and binary bytes. No string dispatch or fallback handler. |
| 2 | Extract material preview owners under `material/derivation/`: digest, PCR, terminal processing, overlap extension and annealing. Move each declaration, revision calculation, producer and validation together. | The private fresh-material/record/capsule replay kernel remains authoritative for ancestry, chemistry and bounds. Avoid DTO-only files and impl fragments. `DigestPreview` is currently declared near line 122 and implemented near 1,590. |
| 3 | Move feature reset/install behavior out of App's `apply_restored` callback into existing narrow workspace, read, PCR, assembly, library and viewer owners. App composes restoration. | Bump the workspace epoch before resets and preserve atomic restoration, cancellation, stale-result rejection and document/event identity. Do not replace the callback with a universal state bag. |
| 4 | Separate cloning proposal/evidence owners and the prepared-pool UI from `cloning_methods.rs`. Colocate suggestions, evidence and binding/request/apply methods; move pool context, receipt, evidence, settings and component together. | Preserve all source/options/target/epoch pins and cancellation guards. Keep actual chemistry in the core. |
| 5 | Convert natural operations while moving their owners: `Record` source/run/derivation/graph queries; `Design::validate`; `Resolved::compile/evaluate`; library active/allocate/matrix/enumeration operations. Put slot geometry comparison on `Slot`. | Update callers directly. Keep explicitly Draft-owned mutation transactions on `Draft`; symmetric algorithms remain free functions when appropriate. |
| 6 | Move generic required-nullable decoding and bounded JSON identity helpers out of `cloning/assembly/recipe.rs` into a narrowly named codec owner. | Preserve exact serialized/hash bytes, bounds and validation. No old-path reexports; this changes ownership rather than scientific contracts. |
| 7 | When editing repository tree jobs, move `folders.rs`'s `add_node` and `scan_step` onto `TreeJob`. | Preserve queue ordering, node/depth/page limits, transport rejection and checkpoint atomicity. This is a small cleanup, not a delivery blocker. |

`Draft` already has a cohesive production implementation and separate tests.
Do not distribute it again. The 17-file sequence-format test directory should
gain registered semantic groups before further growth; preserve explicit Cargo
test discovery when regrouping it.

The MCP/native/sharing review found no reason for another large refactor in its
scope. MCP's largest source is 1,641 lines; desktop backup is 2,069 lines,
including about 880 test lines. Keep rmcp's domain-specific `CollageneServer`
impls: the macro requires that registration boundary while `Session`/`Store`
own mutable state and atomic publication. Repository/library impl extensions
have actual remote-job and transaction responsibilities. Do not split cohesive
backup code just for a lower line count. MCP analysis/cloning and library-sharing
test directories each have 16 files and need the same proactive grouping policy.

When next editing repository remote/folder/publication, lab-library transaction,
desktop catalog activation or material-closure owners, replace broad parent
imports with explicit dependencies where that improves clarity. These are
existing domain extensions, not compatibility shims; no separate cosmetic
rewrite or extra build batch is warranted.

Complete the active scientific batch's integration before changing these shared
interfaces. Then implement connected ownership batches in the order above,
updating production callers, wire schemas where actually necessary and fixtures
together. Items 5 and 6 accompany the relevant owner moves rather than forcing
separate mechanical build cycles. Work on independent core and UI owners may run
in parallel, with one integrator and one shared Cargo cache.

For each completed source batch, run formatting/layout checks, focused native
and WASM checks and affected scientific/state-transition fixtures. Compiler
checks run in the background while independent work continues. Compare executable
test inventories when moving tests. Do not add browser/application HTTP checks
or automatic release verification to this refactor.

## Refactoring constraints

Keep one authoritative verifier, physical-material replay owner and current
internal representation. Move definitions and meaningful implementations together,
update imports and callers directly, and remove superseded code. Do not add
forwarding facades, compatibility aliases, fallback decoders or migrations.

The helper at `/home/mahogny/github/claude/rust_refactor` documents useful
`struct-impl-scatter-stats`, `struct-module-stats` and `module-coupling-stats`
reports. Its write commands automatically run release-profile verification;
that workflow conflicts with this repository's batch/check policy. Its README
was inspected; no automatic transformation, tool build or semantic analysis was
run for this follow-up. Static findings below are review evidence rather than
semantic-tool proof.
