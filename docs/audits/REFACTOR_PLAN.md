# Current code organization plan

28 November 2026. This is the current remaining-work queue, based on three
parallel source reviews. It supersedes the unfinished-work queue in the earlier
[ownership follow-up](organization-followup-2026-11-28.md); completed refactors
remain completed. The audit includes the current uncommitted annotation work,
whose scientific verification is separate. This audit changes the plan, not
production behavior.

## Rules and evidence

The rules are enforced in [AGENTS.md](../../AGENTS.md): at most 20 direct Rust
files per directory, aim below that limit, at most 5,000 lines per file, and an
ownership review at 3,000 lines before further growth. Keep definitions and
meaningful implementations together; use methods for naturally owned operations.
A justified domain, transaction or lifecycle boundary can span files. Stateless
and symmetric algorithms can remain free functions.

The [fresh layout inventory](refactor-layout-2026-11-28.json) covers **1,238 Rust
files in 212 directories**: **zero violations, zero size warnings**, a maximum
of **17 direct Rust files**, and a largest file of **2,869 lines** (`web/app/mod.rs`).
The existing `rust_refactor` binary's read-only `module-size-check`, including
generated sources, independently reports 658 compliant production files. It
excludes 580 non-production files, so it supplements rather than replaces the
repository-wide gate. No helper build, automatic rewrite or compiler invocation
was needed. Semantic scatter scores were not run; ownership findings below come
from source inspection, not a claim of exhaustive semantic proof.

Detailed reviews:

- [Scientific core, model and auxiliary crates](refactor-core-2026-11-28.md).
- [Web application and worker ownership](refactor-web-2026-11-28.md).
- [MCP, native, sharing and repository ownership](refactor-native-2026-11-28.md).

## Ordered implementation batches

| Priority | Batch and concrete change | Completion requirements |
| --- | --- | --- |
| 1 | **Cloning preparation owners.** Separate proposal/evidence handling and prepared-pool handling from `web/src/cloning/preparation/cloning_methods.rs` (2,511 lines). Move each owner's state, request/binding logic, guarded apply behavior and views together. Keep the panel as composition. | Preserve source/options/target/epoch pins, cancellation and stale-response rejection. Choose small complete owners, not DTO-only files or distributed impl fragments. Core retains scientific authority. |
| 2 | **Combinatorial methods.** Put `combination_id`, `validate_construct`, `matrix_size` and `matrix_page` on `Design`; put `validate_matrix` on `SelectedMatrix`. Review their private helpers during the same change. | Directly update core, transport, UI, MCP and test callers. Preserve ordering, pagination, deterministic identities, matrix bounds and exact encoded bytes. Delete old free-function entry points. |
| 3 | **Finish web owner methods and remove two import facades, alongside 1.** Convert naturally owned Draft edit/request operations in place; move context/warning accessors onto `IntendedReference`. Remove the panel's `prepared_rows` reexport and draft model's `MaterialDomain` reexport; import canonical owners directly. | Preserve draft mutation transactions, exact intent/source pins and serialized bytes. Keep Draft definition and implementation together. No old-path forwarding API remains. |
| 4 | **Small native lifecycle cleanup, parallel with 1–3.** Move repository folder `add_node` and `scan_step` onto `TreeJob` in its current owner file. | Preserve queue ordering, transport/generation rejection, page/node/depth bounds and atomic checkpoint publication. Keep repository transactions on `Repository`; do not move them merely because they accept a job. |
| 5 | **Contig identity methods, with the next contig change.** Move `options_revision` to `ContigOptions` and `report_revision` to `ContigReport`. Keep the multi-input batch identity operation free unless a genuine owner emerges. | Preserve hash inputs and bytes, source order and report freshness. No new generic identity abstraction. This does not increase priority of deferred BAM/BAI work. |
| 6 | **Native annotation-editing batch.** Move the complete projection kernel, qualifier owner and canonical INSDC codec down into `model::annotations::{mod, qualifiers, insdc}`. Keep core document builders and physical-evidence checks in core. Use the lower owner for native replace, destination paste and annotation-only strand flips together. | No model-to-core cycle, duplicate parser, old-path reexport or second annotation representation. Keep Location and its meaningful implementation together. Preserve native edit envelopes: insertion inside a feature expands its envelope, unlike source-subset projection. Retain citation replay. |
| 7 | **Proactive registered test/domain grouping.** Before related additions grow existing 15–17-file directories, introduce cohesive child domains. Start with the 17-file sequence-format suite and the 16-file MCP analysis/cloning and sharing suites when their next batch lands. Apply the same rule to web workspace/method groups. | Group actual shared concepts; create no empty speculative tree. Preserve explicit Cargo/module registration and compare executable test inventories before/after any move. The ceiling is not a target. |
| 8 | **Growth review of App and renderer, with relevant feature work.** Move saved restriction setup, antibody, RNA and comparison callback groups into their existing narrow controllers before expanding them. Convert private renderer gesture/draw operations into State methods where they own its invariants. | Keep App as composition and GPU resources with their lifetime owner; preserve async epoch/source guards, cleanup and pointer behavior. Retain pure geometry as free functions. No arbitrary file splitting. |

Items 1–4 are implemented and pass native/WASM checks plus 495 focused cases; the
[implementation and comparison record](ownership-implementation-2026-11-28.md)
tracks integrated verification. The annotation codec's two forwarding functions
were also deleted, with callers using the canonical format owner directly.
Items 1–4 form one coherent ownership batch. Independent workers may own
web, core and native changes concurrently; one integrator owns shared callers,
module registration and verification. Items 5–8 are attached to the next relevant
feature change rather than creating isolated cosmetic compiler cycles. Complete
in-flight scientific work before overlapping its files with a refactor.

## Next connected native editing behavior

Item 6 is now implemented in source: the lower model owner, native edits, paste
and strand flips use one supported-qualifier policy. The former global refusal,
duplicated paste checks and index-paired phase postpass are deleted. The
[native editing contract](../NATIVE_ANNOTATION_EDITING.md) records behavior,
operation-specific bounds and integration evidence. Native/WASM checks and 1,349 focused cases pass; exact-source CI remains required.

Preserve supported assertions when the parent biological calls, order, strand and
phase remain intact. An interrupted parent/exception must require review and
retain affected coordinate/translation/phase declarations only as bounded source
evidence. Identical normalized replacement preserves claims. Annotation-only
strand flips change biological meaning without changing DNA, so deactivate these
claims rather than complementing them as though the same biological calls were
preserved. Keep arbitrary manual feature relocation guarded and unsupported
active qualifier mappings explicitly refused. Delete the duplicate classifiers
and paste-only phase postpass after one authoritative edit policy exists.

Independent fixtures must cover forward/reverse/spliced features, codon phase,
partial TERM, boundary insertion, complete deletion, identity, malformed and
unsupported qualifiers, resource bounds, unchanged sources and citation replay.
These semantics are implemented; fixture execution and exact-source CI remain
separate requirements. The common request budget is unchanged, while the bounded
native edit profile retains 100,000-feature / 32 MiB original-record support.

## Explicit decisions to retain current owners

Do not repeat the completed worker routing, material previews, App restoration,
Draft consolidation, MCP root split or original test relocation. Those findings
were resolved. App (2,869 lines) remains a growth watch: review a new responsibility
before adding it, and extract a complete controller only when its lifecycle
warrants one. Line count alone does not justify scattering App state.

Keep domain-specific rmcp handler impls beside their request schemas; the framework
registration boundary justifies this split while Session/Store retain mutation
and publication authority. Keep repository/library transaction extensions,
desktop backup/activation lifecycles, TOPO graph execution, and the private
material replay/install kernel cohesive. No wholesale scientific or native
rewrite is supported by this review. Large renderer/layout algorithms may stay
free functions when they do not naturally belong to a single input type.

## Batch verification and closure

Implement owner moves, callers, schemas and tests as one coherent source batch.
Remove superseded paths and functions directly; no old-path reexports, delegators,
compatibility decoders or scaffolding. Formatting, whitespace and the repository
layout gate are required. For test moves, compare before/after executable names
and counts, not just source attributes.

After the source batch is complete, run the affected native fixtures and web/worker
WASM check using the existing shared Cargo cache and consistent selections.
Run necessary compiler work in the background while independent implementation
continues. Preserve scientific correctness, data integrity, resource bounds,
revision/cancellation and atomicity checks. Do not rerun unchanged passing cases
without a concrete reason. Browser, screenshot and application HTTP checks remain
deferred until full in-scope code is complete. Exact-source CI and product
acceptance remain separate delivery gates; passing this audit does not complete
the product goal.

## Next citation and edit-log owner

Legal GenBank `(sites)`, unscoped/prose reference descriptions, repeated `bases`
after semicolons, and fully deleted citation coverage still prevent otherwise
valid DNA edits. This is bibliographic handling, not unknown molecular biology.
The next connected batch should move complete citation/edit-log ownership into
`model::genbank_edit`, retaining Molecule and its meaningful implementation together.
Use typed mapped/source-only citation outcomes; keep numbering and bibliography
while explicitly withholding destination coverage for ancestral-only evidence.
Always replay the source-bound edit log independently before export. Malformed,
stale or unsupported-version logs remain errors, never source-only fallbacks.
Update all callers and fixtures directly, with one current encoding and no root
forwarding exports. The [delivery audit](../DELIVERY_AUDIT.md) records the gap.
