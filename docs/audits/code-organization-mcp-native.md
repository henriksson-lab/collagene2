# MCP and native organization audit

Audit date: 2026-11-28. Scope: MCP, desktop hosts/protocols, server, repository, shared-library storage/synchronization/publication. This is a source inspection, not a compiler or correctness assessment. No compiler, tests, browser, application HTTP requests, runtime databases, or credentials were used.

The working tree includes unfinished SBOL integration, including `mcp/src/sbol_tools.rs` and its registrations. Counts below include that working tree; its behavior has not been verified by this audit. Do not discard or overwrite that batch during restructuring.

Implementation update: the 95 original MCP integration targets have now been grouped into eight explicit Cargo harnesses under `tests/{sequence,analysis,cloning,assembly,methods,workspace,publication,protocol}/main.rs`, with automatic discovery disabled. All 176 source-attributed tests retain their original isolated modules; the largest test directory has 16 Rust files. Static review confirmed original file contents are unchanged except the 13 necessary fixture/include path updates, and checked each module chain/manifest entry. The native build, executable discovery and focused SDK cases pass; see the [verification record](refactor-verification-2026-11-28.md) for the exact scope. Historical counts and recommendations below describe the pre-refactor audit baseline.

## Measured baseline

The rules are at most 20 Rust files in **every** directory, at most 5,000 lines in **every** Rust file, and review/refactor pressure at 3,000 lines. These limits include tests. Counts are physical `.rs` files, not Cargo modules.

| Crate | Rust files, whole crate | Most files in one directory | Largest Rust file, lines |
| --- | ---: | ---: | --- |
| mcp | 112 | **95**, `tests/` | **8,973**, `src/lib.rs` |
| desktop | 7 | 4 | 422, `src/main.rs` |
| desktop-credentials | 6 | 5 | 301, `src/linux.rs` |
| desktop-host | 12 | 8 | 2,068, `src/catalog_backup.rs` |
| desktop-mcp | 2 | 1 | 307, `tests/session.rs` |
| desktop-protocol | 5 | 3 | 814, `src/lib.rs` |
| desktop-sync-host | 8 | 5 | 1,393, `tests/service.rs` |
| desktop-sync-protocol | 7 | 4 | 1,043, `src/lib.rs` |
| lab-library | 9 | 5 | 1,119, `src/lib.rs` |
| library-host | 4 | 2 | 636, `tests/boundary.rs` |
| library-protocol | 4 | 3 | 800, `src/lib.rs` |
| library-sharing | 21 | 16 | 744, `src/lib.rs` |
| library-sync | 3 | 2 | 729, `src/lib.rs` |
| local-storage | 2 | 1 | 284, `src/lib.rs` |
| repository | 12 | 7 | 1,001, `src/folders.rs` |
| server | 7 | 5 | 1,937, `tests/router.rs` |

Only MCP violates the hard rules in this scope. Its source directory already contains 17 Rust files: adding more flat helper modules would soon violate the directory limit as well. None of the native crates requires indiscriminate file splitting.

## Findings and ownership

### P0: MCP root combines the entire application surface

`mcp/src/lib.rs` mixes request schemas (roughly lines 43–1831), session records and quota logic (1833–2071), encoding/output helpers (2073 onward), `CollageneServer` state and scientific/storage helpers (2265–3309), `ServerHandler` resource transport (3313 onward), and the single tool router (3438–8863). This is an ownership problem beyond its line-count failure.

Put request structures, enum conversions/defaults, tool implementations, and domain response shaping together. For example, `SavedTopoTask` and its `preparation`/`native` methods currently reside in root while TOPO behavior also lives in `topo_actions_tools.rs` and `topo_graph_tools.rs`. They belong to a TOPO family, with a split between physical preparation and annotation projection only where those responsibilities warrant it. `AssemblyIntentRequest`, `SubstitutedIntentRequest`, and `ContributionEditIntentRequest` already have useful `native` methods; retain their methods beside their definitions when moving them.

Do not replace the monolith with arbitrary `types_1.rs`, `handlers_2.rs`, or a single giant `types.rs`. Avoid retaining the old root request paths with forwarding reexports. Change Rust callers to their current owning module paths in the same batch. Preserve MCP tool names and JSON shapes: these are the deliberate external product interface, not an obsolete internal Rust layout.

### P0: session mutation must retain one atomic owner

`Store` (1964) owns document revisions, plans, artifacts, uploads, workspace metadata, evidence collections and shared byte accounting. `StoredWorkspace` has sensible `split`, `restore`, and `summary` methods (1866 onward); move that complete owner rather than distributing its methods between generic utility files.

The server currently owns substantial publication/storage behavior: `text_artifact` (2335), `candidate_document_artifact_bounded` (2501), `paired_text_artifacts` (2670), `binary_artifact` (2769), `restore_workspace` (3028), and `commit_workspace_outputs_named` (3092). These should become methods of session storage, or a narrowly scoped transaction owner borrowing that storage. The server should own scheduling/transport and invoke one storage transaction after computation. This is a real ownership improvement; merely moving 15 `impl CollageneServer` blocks to 15 files would preserve the god object.

Do not split each map into independent mutexes or unrelated stores solely to reduce file size. Workspace publication currently reserves capacity across documents, metadata and artifacts before publishing. Preserve single-lock atomic commit, duplicate accounting, no partial outputs on failure, upload limits, stale-reference rejection, revision identity, and source-authority separation. Physical session submodules may cooperate through restricted visibility without inventing a second storage representation.

Keep `Artifact` with its encoding/size/access behavior. `ArtifactData::Text` and `Binary` represent genuine payload alternatives. SBOL imported-source membership is an authority distinction; moving it must not convert arbitrary stored JSON artifacts into accepted source handles.

### P1: router separation is supported by the installed dependency

The checked local `rmcp-macros` 3.11.0 documentation explicitly describes distinct named `#[tool_router(router = ..., vis = ...)]` impl blocks and adding their `ToolRouter` values. Use named domain routers and compose them once during server construction. This allows actual tool methods and their schemas to live together without maintaining duplicate root forwarding handlers.

Keep the explicit `#[tool_handler(router = self.tool_router)] impl ServerHandler` boundary and its custom resource behavior. Avoid deriving a separate server handler for every router. Check generated schema/attribute name resolution after moves; macro-generated associated methods can require restricted visibility across sibling modules. Use explicit imports, not `use crate::*` as currently found in `sbol_tools.rs`. Keep helper structs/functions private unless sibling access actually requires `pub(super)`/`pub(crate)`.

The current declared tool count is 208, including the unfinished five-tool SBOL batch. This is a source observation, not proof all 208 registrations compile. Refactoring acceptance must compare actual registered names, duplicate detection and schemas, not simply update a numerical assertion until it passes.

### P0: moving integration tests can silently remove coverage

There are 95 automatic Cargo integration targets directly in `mcp/tests`; no explicit target declarations currently compensate for nested files. Moving them into folders alone makes Cargo stop discovering them.

Use one harness entry file per semantic suite at `tests/`, each with a same-named subdirectory and conventional `mod.rs` for deeper groups. For example, `tests/cloning.rs` explicitly declares `#[path = "cloning/mod.rs"] mod suite;`, and `tests/cloning/mod.rs` declares its recipe/assembly submodules. The explicit path avoids the same-named root file colliding with directory-module lookup. This preserves ordinary Cargo discovery while putting every physical directory below 20 files. Put each former test file in a distinct Rust module to preserve helper-name isolation and recognizable filtered test names. Update commands/documentation that select former Cargo target names directly; do not retain empty compatibility targets.

Bounded in-memory SDK sessions are repeated, for example `Session` in `tests/cloning_gel_tools.rs`. A small `tests/support/mod.rs` may own generic JSON-RPC session lifetime/call helpers. Keep biological fixtures and domain assertions in their own suites. Do not introduce a test fixture framework or consolidate independent scientific expectations into implementation-derived golden output.

### P2: native crates mostly have useful existing boundaries

Desktop credentials, explicit filesystem policy, stdio MCP hosting, desktop IPC and synchronization are already distinct crates. Keep these boundaries; avoid folding credentials or trusted native paths into shared request structures.

`desktop-host/src/catalog_backup.rs` has 2,068 lines, including test modules beginning at 1189 and 1942, and an existing `catalog_backup/activation.rs`. Its cohesive service is below the warning threshold. A low-risk improvement is to make it `catalog_backup/mod.rs`, with `CatalogBackupService` and its impls together, the `DesktopCatalogBackups` shell-facing owner in `catalog_backup/desktop.rs`, existing activation beside them, and private test modules under `catalog_backup/tests/`. Extract file-record I/O into a cohesive private owner only if there is demonstrated duplication or independent lifecycle; do not scatter a service method per file.

`DesktopHost` is defined in `desktop-host/src/lib.rs` with a library-specific impl in `library.rs`. `Repository` has impls in root, `remote.rs`, `folders.rs`, `folder_publication.rs`, and `remote_removal.rs`. `Library` has its root impl plus transaction/collaboration impls. These are understandable domain extensions, but should be treated as named owner families when reviewing scatter. Future growth should use `host/`, `repository/`, or `library/` subtrees around the owner, not another set of disconnected root impl files.

There is one concrete method opportunity in `repository/src/folders.rs`: private `add_node(job, node)` (296) changes only `TreeJob` and checks that job's bounds/queue, so make it `TreeJob::add_node`. Consider `scan_step` (319) as a job method with explicit transport only after checking that its ownership remains clear. Keep symmetric conversion algorithms and shared validation free functions when no receiver naturally owns them. Do not turn every function into a method to satisfy a metric.

`lab-library/src/lib.rs` is 1,119 lines and `mutations.rs` 1,038; there is no current size-driven need to split SQL, transactions and owned `Library` methods further. Any later extraction must preserve one transaction boundary for quota, revision insertion, generation/CAS updates and mutation receipts.

`library-sharing/src` already has 16 files, approaching the cap. Group future work into semantic `publication/`, `import/`, and `evidence/` families before adding five more siblings. Existing `sbol_source.rs` and `sbol_subset.rs` have different source-closure and subset-publication responsibilities; proximity does not justify merging their authority decisions.

Server routing and repository HTTP code are below warning size. Preserve separate browser/server authentication and native credential scopes; this audit found no reason to rewrite their public transport formats.

## Proposed MCP subtree

This is a responsibility map, not a mandate for an empty file at every leaf. Each directory should target roughly 8–15 Rust files, with a hard ceiling of 20. Target coherent files below roughly 2,500 lines to leave growth room; review anything reaching 3,000.

```text
mcp/src/
  lib.rs                 module declarations/current public API only
  main.rs                native stdio entry point
  server/
    mod.rs               CollageneServer + construction/scheduling
    registry.rs          named router composition, registry metadata
    resources.rs         ServerHandler resources/tool dispatch
    response.rs          bounded serialization/ErrorData mapping
  session/
    mod.rs               Store + shared accounting/transaction owner
    references.rs        document/workspace/evidence reference contracts
    artifact.rs          Artifact + bounded artifact operations
    workspace.rs         StoredWorkspace/Candidate + restore/publication
    upload.rs            WorkspaceUpload + lifecycle
    evidence.rs          stored collection lookup/publication
    genome.rs            GenomeResources/NativeSearch ownership
    genome_store.rs      existing managed native dataset implementation
  tools/
    mod.rs
    sequence/            import/export, tracks/annotations, editing
    cloning/             digest, PCR, end preparation, intent, plans
    assembly/            ligation/homology/Golden Gate/recombination/TOPO
    analysis/            alignment, variants, reads, proteins, CRISPR
    workspace/           archive, history, runs, material origins
    publication/         combinatorial releases, assembly libraries
    exchange/            SBOL, structures, quantities/protocol outputs
    catalog/             enzyme lookup/evidence
```

Each family contains its request types and actual registered methods. Further partition `assembly/` by method family if it approaches 20 files; do not put all current helper filenames into one flat `tools/` directory. Move existing unit/quantity helpers into the relevant protocol domain rather than duplicating validation. `DocumentRef`, `WorkspaceRef`, and `EvidenceRef` can share a reference-contract module while retaining their distinct meanings.

Proposed integration harnesses: `sequence`, `cloning`, `assembly`, `analysis`, `workspace`, `publication`, `protocol`, and `session`. Each owns a bounded subtree; subdivide assembly by chemistry. A mapping of all 95 old files to new modules and of old target-specific commands to current suite filters is part of the patch, not optional follow-up work.

## Ordered implementation batches

1. **Inventory and integration ownership.** Freeze concurrent writes to `mcp/src/lib.rs` during its move. Record current tool names/request schemas, public Rust consumers, every integration target and test name. Incorporate current SBOL work directly; do not create an old/new exchange branch. Add the organization rules to repository guidance and prepare a metadata-only size/count gate.
2. **Session and transport ownership.** Move complete session records/impls, references, bounded responses, resource handling and server construction. Change callers directly. Retain atomic workspace/artifact publication and compute-gate behavior. No forwarding root APIs. Keep in-module access tests beside the owner or update their restricted imports.
3. **Domain routers and schemas.** Move request types and actual handlers together into named domain routers; compose them once. Carry helper implementations into the same concept family. Extract storage-owned methods to the session transaction boundary, avoiding behavior changes to chemistry verification.
4. **Integration suites.** Group all 95 tests into explicit discoverable harnesses, preserve independent scientific/data-integrity assertions and generic session teardown. Update smoke scripts and target selectors directly. A common SDK helper extraction can happen here after comparing actual repetition.
5. **Native housekeeping, lower priority.** Apply the `TreeJob` ownership improvement and, if useful to upcoming work, the backup owner subtree. Group sharing modules before approaching the directory cap. Leave small cohesive crates intact.
6. **Combined verification and gate.** Finish each assigned integration batch before compiling. Use shared-cache background native checks while continuing unrelated documentation/edits; collect failures and fix together. Run focused MCP integration suites plus embedded desktop MCP/session checks and relevant native host/repository data-integrity tests. Native-only modules need native validation; WASM checks belong to touched shared core/schema/protocol consumers, not unsupported native targets. Preserve existing platform CI coverage and require Linux/macOS/Windows builds for platform-sensitive path/credential changes. No browser checks are required for a source-only reorganization.

Parallel work is safe by disjoint semantic families after the root/session interfaces are agreed. One integrator owns root/router composition and shared Cargo target selection. Workers should not each compile their patch, rewrite the shared root, or update the same schema imports concurrently.

## Using rust-refactor without conflicting with repository rules

Read helper sources: `README.md`, `docs/MODULE_REFACTOR_WORKFLOW.md`, and `docs/MODULE_SIZE_CHECK.md` under `/home/mahogny/github/claude/rust_refactor`. Suitable analysis commands include `struct-impl-scatter-stats`, `struct-module-stats`, `module-coupling-stats`, `visibility-stats`, and `to-oop-stats`. Interpret macro/cfg uncertainty explicitly; the MCP proc-macro surface is not ordinary handwritten callers.

Do not blindly follow the helper workflow's temporary compatibility glob imports/reexports or `--keep-delegator`: those conflict with this repository's direct-update rule. Its README also states that **all write-command verification uses release builds**, which conflicts with the repository's prohibition on automatic release verification. Use read-only inventories and reviewable dry-run plans where they perform no builds; apply manual coherent direct edits unless the helper can honor this repository's verification policy.

The helper's `module-size-check` excludes integration tests/benches/examples, and does not implement the per-directory Rust-file count rule. It is insufficient as the sole repository gate. Add a fast source inventory that covers all tracked repository `.rs` files, excludes build output/dependencies explicitly, fails above 5,000 lines or 20 files/directory, and reports 3,000-line warnings. Do not exempt a source just because its filename suggests generated output. Ownership review remains necessary; size compliance alone cannot detect scattered structs.

## Acceptance criteria

- No scoped directory exceeds 20 `.rs` files; no scoped Rust file exceeds 5,000 lines; files at or above 3,000 have been redesigned or explicitly flagged for the next coherent batch.
- Root MCP is a small entry point; schemas and actual handlers belong to named semantic families. Struct definitions and their methods are co-located or have a documented, domain-specific owner-family reason.
- Tool names, descriptions, schemas, resources and externally meaningful response contracts match the pre-refactor baseline, except any separately authorized unfinished feature changes. Router construction detects accidental omissions/duplicates; no count-only acceptance.
- Every pre-refactor integration test remains discoverable, or its removal has a concrete obsolete-only justification. Cargo target/filter commands resolve the new suites directly.
- Store quota and stale-reference checks, atomic publication, cancellation, revision/authority validation, safe filesystem boundaries and credentials remain covered by meaningful existing tests.
- All callers use current module paths. No compatibility root reexports, aliases, fallback encodings, duplicate implementations or empty legacy harnesses are introduced.
- Focused checks for the complete batch pass before declaring implementation finished. This audit itself provides no build/test pass claim.
