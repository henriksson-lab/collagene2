# MCP and native ownership audit, current source

28 November 2026. Source-only follow-up covering MCP, desktop/Tauri hosts and
protocols, native credentials/storage, library sharing/synchronization, repository
and server. This audit includes the current uncommitted qualifier/reference-edit
batch. It does not certify that batch or supersede its compiler/test evidence.
No production changes, compiler/test invocations, browser/application requests,
runtime database reads, credentials or release builds were used.

## Measurements

A physical-source inventory found 257 Rust files across the 16 scoped crates.
No scoped directory exceeds 20 direct Rust files, no Rust file exceeds 5,000
lines, and no file reaches the 3,000-line ownership review threshold. The
repository-wide layout script remains the required gate; these scope counts
are supplementary evidence, not an exemption or baseline.

| Scope | Rust files | Maximum direct files/directory | Largest file, lines |
| --- | ---: | ---: | --- |
| MCP | 147 | 16 | `src/tools/analysis/mod.rs`, 1,641 |
| Desktop host | 12 | 8 | `src/catalog_backup.rs`, 2,069 |
| Desktop sync host | 8 | 5 | `tests/service.rs`, 1,396 |
| Repository | 12 | 7 | `src/folders.rs`, 1,001 |
| Lab library | 9 | 5 | `src/lib.rs`, 1,119 |
| Library sharing | 22 | 16 | `src/lib.rs`, 763 |
| Server | 7 | 5 | `tests/router.rs`, 1,939 |

Other scoped crates are desktop, desktop-credentials, desktop-mcp,
desktop-protocol, desktop-sync-protocol, library-host, library-protocol,
library-sync and local-storage. None needs size-driven splitting. Counts
include test/example/build-script sources wherever present.

The former 8,973-line MCP root and 95-file integration-test directory are
resolved. The root is now a small entry; mutable session/storage behavior lives
in `session/mod.rs` (1,250 lines). Eight explicit Cargo integration harnesses in
`mcp/Cargo.toml:32` retain discovery with `autotests = false`. Do not reopen the
historical monolith findings as current defects.

## Refactor queue

### 1. Repository tree job methods: small, concrete ownership improvement

`repository/src/folders.rs:134` defines `TreeJob`; its existing inherent impl at
line 163 contains `progress`. The free function `add_node` at line 296 operates
only on a job: it checks node/depth limits, uniqueness, updates its discovery
queue, and appends a node. `scan_step` at line 319 implements the same job's
remote scanning state machine with an explicit transport dependency.

Move these to `TreeJob::add_node` and `TreeJob::scan_step` in the same owner
file, update all callers directly, and delete the old functions. Inspect
`browse` during the move; retain it as a free transport operation if its job
argument is merely configuration rather than owned lifecycle. Do not move
repository checkpoint/transaction methods onto the job simply because they
accept it. This is a low-risk cohesive cleanup, not a large native rewrite.

Preserve queue ordering, discovery/page/generation checks, node/depth bounds,
remote rejection and atomic checkpoint publication. Existing tree-job cases
are the appropriate focused verification. No wire version change is warranted
for a receiver-only refactor.

### 2. Proactive integration-test grouping before the next related additions

Three directories have 16 direct Rust files:
`mcp/tests/analysis`, `mcp/tests/cloning`, and `library-sharing/tests`.
The MCP assembly suite has 13. These are below the ceiling but should gain
cohesive registered child groups before substantial further growth.

Suggested grouping follows current domains: MCP analysis read/alignment/variant
review versus protein/RNA analysis; cloning PCR/overlap-extension versus terminal
preparation/digest; sharing physical-material closure versus selected-document
publication and import. Choose groups from actual incoming work instead of
creating empty future subdirectories. Keep every former source module isolated
where helper names overlap. Record the pre-move executable test inventory,
explicitly register nested modules/harnesses, and compare after discovery.
Moving files alone must never reduce coverage.

### 3. Narrow imports during the next substantive owner change

Production broad parent imports remain in repository remote/folder/publication
modules, lab-library transaction/collaboration/backup modules, desktop catalog
activation, native sync publication/removal, and material derivation closure.
Examples: `repository/src/folders.rs:4`,
`desktop-host/src/catalog_activation.rs:2`, and
`library-sharing/src/material_derivation_closure.rs:4`.

Replace these with explicit dependencies as those domains are next edited.
This makes cross-owner coupling reviewable. It does not justify a separate
cosmetic batch or scattering one atomic owner across files. Test-local
`use super::*` is not the same concern.

## Reviewed splits that should stay

- `CollageneServer` has domain-specific impls across MCP tools, but these contain
  the actual rmcp registered handlers alongside their request schemas. The
  server owns compute scheduling and router composition (`server/mod.rs:75`);
  `Session`/`Store` own state, quota and atomic publication
  (`session/mod.rs:176`, `:289`). Keep this registration/lifecycle separation.
  Do not replace domain methods with root forwarding handlers or introduce a
  universal mutable context to reduce an impl-scatter score.
- `StoredWorkspace` and its meaningful `split`, `restore`, `summary` methods
  are together (`session/mod.rs:68`, `:78`). Session artifact/workspace commit
  methods share one quota/lock boundary. Splitting every map into a different
  store would weaken the design. There is no current size-driven need to split
  this owner.
- `Repository` extensions in remote, folder import/publication, publication
  review and remote removal represent distinct job/transaction lifecycles.
  `Library` extensions in mutations and collaboration similarly preserve SQL
  transaction boundaries. Their extra impl files have domain reasons. No
  generic rewrite is justified solely by counting files.
- `CatalogBackupService` is defined in `desktop-host/src/catalog_backup.rs:94`
  and has an activation-specific impl in `catalog_activation.rs:226`. Backup
  creation/recovery and append-only restart-bound activation are different
  lifecycles. The 2,069-line backup file includes approximately 880 test lines;
  retain it until a real responsibility change warrants a `catalog_backup/`
  owner subtree. Do not split it just to obtain smaller numbers.
- The native sync `publication::Reviews` and `removal::Reviews` are different
  types with similar bounded queue mechanics, not one struct scattered between
  files (`publication.rs:5`, `removal.rs:5`). Their owned review payloads and
  execution authority differ. Their duplication is small; do not merge their
  admission/consumption rules merely for code reuse. A future shared private
  bounded queue would need to preserve typed payloads and independent queues.
- Sharing's digest and material-derivation `Index` types are separate module
  owners. Stateless document-closure comparison/merge algorithms in `import.rs`
  operate symmetrically on source/destination documents; making them methods on
  one document would imply misleading ownership.

## Integration requirements

Apply the small tree-job method move with its callers as one batch. Test
regrouping is independently parallelizable after agreeing on module registration;
one integrator owns the Cargo harness inventory and compilation selection.
Keep meaningful scientific/source-authority, filesystem, revision, cancellation,
quota and atomicity tests. No compatibility exports or delegating old methods.

The read-only helper in `/home/mahogny/github/claude/rust_refactor` may assist
future semantic scatter/coupling analysis, but this audit used direct source
inspection. Its automatic write-command release verification conflicts with
this repository's check policy. No helper build or automatic rewrite was run.
