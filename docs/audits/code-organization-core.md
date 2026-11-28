# Scientific core organization audit

Audit date: 2026-11-28. Scope: the current `crates/core` source and integration tests, including uncommitted work visible at audit time. This is a static organization audit; no build, scientific test, browser check, application HTTP request, or refactor-tool execution was performed. Source behavior and current SBOL integration have not been verified by this audit.

## Findings and priorities

Implementation update: the model facade is deleted, worker protocol owners now
live under `src/protocol/`, and all 320 test targets are registered as nine domain
suites. The table and source anchors below record the audit baseline. Integration
verification passes for this batch; verifier and remaining source-family moves are open. See the [verification record](refactor-verification-2026-11-28.md) for its exact scope.

| Priority | Finding | Evidence | Required change |
| --- | --- | --- | --- |
| P0 | Worker transport exceeds the absolute file limit. | `src/protocol.rs`: 5,797 lines. | Separate transport responsibilities and domain execution before adding more worker operations. |
| P0 | Source and integration tests exceed the directory limit by large margins. | `src/`: 156 direct `.rs` files, 157 recursively; `tests/`: 320 direct `.rs` files, 328 recursively. Existing child directories are `src/annotation_formats/` (one Rust file) and `tests/common/` (eight). | Introduce concept subtrees; move callers and test registration in the same batches. |
| P1 | Verification mixes plan declarations, scientific reconstruction, previews, and protocol rendering. | `src/verification.rs`: 4,528 lines; already above the 3,000-line warning threshold. | Extract coherent owners and rendering while retaining one authoritative verifier. |
| P1 | Explicit compatibility facade remains in new software. | `src/model.rs:1` describes itself as a compatibility facade and only reexports `collagene_model::*`; `src/lib.rs:16` reexports through it. | Delete this facade and update imports to the canonical model crate, including downstream consumers. |
| P2 | Native material derivation is the next likely oversized owner. | `src/material_derivation.rs`: 2,729 lines, with preview type declarations distant from their behavior. | Reorganize by owned preview/derivation concepts before further growth. |
| P2 | Module names encode relationships through repeated prefixes rather than directory boundaries. | Flat `assembly_*`, `overlap_extension_*`, `topo_*`, `run_*`, `pairwise_*`, `read_alignment_*` families. | Move families into coherent namespaces with explicit dependencies. |

The next largest source files are `contig.rs` (2,249), `assembly_recipe.rs` (2,235), `pdf.rs` (2,034), `combinatorial_design.rs` (1,939), and `topo.rs` (1,879). These do not require splitting solely for line count. Move their complete owners first; split only when responsibility or size justifies it.

## Worker transport: split by responsibility, preserve one schema

Current anchors in `src/protocol.rs`:

- Lines 9–157: current schema versions, envelope header/checksum/size handling, bounded serialization, document save/load.
- Line 159: `AnalysisJob`; line 1,121: `Request`; line 1,127: `AnalysisResult`; line 1,330: `Response`.
- Line 1,336: request validation/encoding; line 1,945: response decoding and validation.
- Line 2,789: stateless execution, revision binding, and result construction in one long dispatcher.
- Lines 5,581–5,787: `AnalysisSession` definition and stateful execution.

Proposed owner family: `worker/{mod.rs, schema.rs, envelope.rs, request.rs, response.rs, session.rs, dispatch/mod.rs, dispatch/<domain>.rs}`. Keep each request or response type and its inherent behavior together; `schema.rs` owns the two cross-domain enums and versions, not a second schema. Keep `AnalysisSession` and its implementation together in `session.rs`; its private genome/read/variant/combinatorial state remains private.

Extract complete domain handlers from the 2,789+ dispatcher. A short exhaustive central match routes typed operations; domain handlers perform the actual revision-bound operation. Avoid optional handler chains, catch-all fallback dispatch, repeated decoding, or string routing. Group handlers into cloning, sequence, evidence, tracks, and document families, introducing child directories before any family exceeds 20 Rust files. Do not introduce a file per enum variant.

`WORKSPACE_REQUEST_SELECTORS` at line 14, numeric response selectors near line 1,945, and session tag routing near line 5,596 encode serialized enum ordering. A mechanical relocation must retain variant order and the exact current wire representation. If domain work deliberately changes that representation, change the current schema version and every caller/fixture together and explicitly reject the old version. No dual decoder or old enum aliases. Preserve predecode resource bounds and checksum/version validation before routing or allocation.

## Plan and verification ownership

The current file has useful existing local ownership: `LigationOutcomeEvidence` has its own implementation at line 288, `MaterialStateEvidence` at line 366, and `DerivedOriginEvidence` at line 426. Preserve those pairings rather than moving declarations into a generic shared-types file.

Specific improvements:

1. `Operation` is defined at `verification.rs:45`, while `operation_sources` and `operation_sources_mut` are free functions at lines 1,247 and 1,291. They inspect only the operation's own variants. Move them beside the enum as `Operation::sources()` and `Operation::sources_mut()` and update callers directly. `plan_editing.rs:5` currently imports both visitors.
2. `DesignPlan` at `verification.rs:217` owns inputs, ordered steps, and target revision. `plan_revision` at line 703 is naturally `DesignPlan::revision()`. Put the struct and its direct invariant/identity behavior in `planning/plan.rs`.
3. `plan_editing.rs:108`, `149`, `185`, and `218` implement moved/removed/inserted/replaced operations on input IDs and step arrays. Review these as plan graph behavior, not generic helpers: either make atomic `DesignPlan` editing methods if callers own a complete plan, or introduce one clearly owned graph-editing value when editing partial drafts is a real requirement. Do not force callers to invent target revisions merely to call a method. Keep lower-level symmetric graph algorithms free where that is clearer.
4. `verification.rs:2,328` starts `reconstruct`, passing products, findings, step reports, physical processing reports, reference reports, tracked physical states, final Golden Gate evidence, and a cleavage budget separately. Introduce a private verifier execution owner holding these correlated states and their bounds; make run/step execution methods on it. It should consume the single `DesignPlan`/`Operation` representation, not copy them into a private alternative plan encoding. Retain separate lifetimes and immutable inputs where needed.
5. `PcrDocument`/PCR functions at lines 463–586, `MutagenesisPreview` at 587–702, `RepairPreview` at 1,714, `AnnealPreview` at 1,724, and `TailPreview` at 1,893 belong with the relevant preview execution. Move complete type/function groups into PCR and terminal-chemistry owner modules. Cross-domain annotation projection belongs in a small shared projection module, not duplicated in each handler.
6. `verification.rs:3,995–4,528` contains HTML escaping, protocol tables, `CheckedProtocol`, and HTML rendering. Move to protocol-output ownership. `CheckedProtocol` (4,086) must remain non-deserializable and constructible only from fresh verification. `pdf.rs:1,395–1,401` already obtains it before rendering PDF; keep HTML and PDF on this same checked context. Method candidates are `CheckedProtocol::new`/render operations; do not accept a caller-supplied report as authority.

Suggested verification family: `verification/{mod.rs, engine.rs, evidence.rs, projections.rs, sources.rs, operations/mod.rs, operations/pcr.rs, operations/digest.rs, operations/ligation.rs, operations/overlap.rs, operations/terminal.rs, operations/topo.rs, operations/recombination.rs}`. Evidence types and their methods stay together. Keep a bounded exhaustive operation dispatch in the engine; handlers receive only the context they need. Avoid scattering multiple `impl Operation` blocks among the scientific handlers.

## Cohesive source subtrees

This is a namespace proposal, not a requirement to split every existing file. Use explicit imports and direct updated paths; keep root `lib.rs` as a small module declaration list. Aim for 6–12 direct Rust files per directory, allowing up to 20 only where the concepts justify it.

| Destination | Current related files or concepts | Ownership notes |
| --- | --- | --- |
| `sequence/` | `analysis`, translation/reference, editing/templates, annotations/batch/library, search/similarity, codon design | Further divide annotations and design if the direct count approaches 20. Model types live in `collagene-model`, not a new core facade. |
| `alignment/` | `alignment`, pairwise context/export/history/review/search, dotplot, antibody alignment | Pairwise record behavior stays with its record owner. External alignment formats belong under format exchange. |
| `cloning/restriction/` | catalog digest, enzyme sets, restriction selection/inspection/history, digest setups/cut selection, run catalog digest | Keep physical digestion and site inspection distinguishable. |
| `cloning/pcr/` | primers, PCR source/outputs, mismatch outputs/inspection, mutagenesis, material PCR, primer reuse | Shared PCR geometry stays authoritative; move preview owners with their implementations. |
| `cloning/overlap/` | overlap profiles, homology assembly/annealing, three-prime resection, overlap extension/PCR/execution/primer design/inspection | Distinct physical protocols remain distinct algorithms. |
| `cloning/assembly/` | intent/edit/substitution intent, junctions/editing, edited primer design, diagnostics, library history, selective junction editing | Keep intent declarations distinct from reconstructed material evidence. |
| `cloning/golden_gate/` | Golden Gate/reaction, fusion design/rank, primer/junction design, prepared reactions, domestication/preparation/rounds | About a dozen current files; no need to flatten into the parent assembly directory. |
| `cloning/ligation/` | ligation profiles/fidelity/outcomes/recipe, tailing/recipe, terminal chemistry previews | Preserve alternative products and selected intended path in recipes and output. |
| `cloning/topo/` | TOPO algorithm/actions/preparation/run preparation/annotations/graph history/worker | Do not merge reference annotations with source-backed physical reconstruction. |
| `cloning/recombination/` | recombination adapters and multisite recombination | Dedicated scientific crates remain canonical; no algorithm duplication. |
| `planning/` | plan/operation definitions, editing, assembly recipes, oligos, combinatorial design/batch/export/release | Combinatorial sessions and recipes can each be child namespaces rather than pushing this directory to the limit. |
| `material/` | run material/capsules, material derivation/previews, lineage, run evidence/review, physical/reference evidence, cloning outputs/reference transfer | Use `material/derivation/` for preview owners. Avoid a large generic `types.rs`. |
| `documents/` | workspace, document history, transfer, contig/RNA history | Keep `Workspace` declaration and ordinary methods together; archive codec is a justified separate responsibility. |
| `formats/` | import dispatch, GenBank/FASTA, EMBL, SnapGene, XDNA, SCF, GCG, original Collagene XML, oligo/annotation/protein/alignment exchange, fidelity checks | Use `formats/sequence/`, `formats/annotations/`, `formats/alignment/` so format expansion cannot exceed 20 direct files. Move existing nested annotation-format code with its owner. |
| `output/` | checked protocol context/HTML, PDF, SVG/export, protocol layouts/quantities, gel | Layout and quantities each keep their request/report validation implementations together. Gel prediction may instead live in analysis if that simplifies dependencies. |
| `review/` | contigs, trace validation/traces, read checks/review, indexed BAM and read reference/store/archive/session, variant reference/archive/session | Divide into `reads/`, `variants/`, and `contigs/`; BAM/BAI feature expansion remains postponed. |
| `tracks/` | tracks/view tracks, genome/genome session/CRISPR | Preserve prepared-index sessions and resource bounds. |
| `worker/` | transport schema, envelope/codec, domain dispatch, analysis session | Transport depends on science; core scientific algorithms should not depend on worker dispatch. |

Chemistry, modified-duplex tracking, methylation, protein analysis/structures, CRISPR, and other remaining owners must be explicitly assigned in the move manifest rather than left as an unreviewed flat remainder. Chemistry should have its own subtree if terminal methods make it grow; protein analysis can be its own small subtree. The proposal intentionally leaves fewer than 20 root module declarations and only a handful of root Rust files.

### Material ownership details

`material_derivation.rs:34–45` declares macro-generated selection/request/record types, line 104 declares its capsule, and line 111 its digest preview. The digest preview implementation is at line 1,522; additional preview implementations appear at 1,827, 1,895, 2,037, 2,265, and 2,656. Keep each preview's declaration, validation, revision computation, and preview operation together when extracting PCR/digest/terminal/annealing owners. Preserve one shared capsule/record replay owner for ancestry validation and byte budgets.

`run_material.rs:15–22` reexports a substantial capsule API implemented by `run_material_capsule.rs` (capsule declaration at line 125, implementation at 577). This is a responsibility split rather than proof of an obsolete implementation. Make it an explicit nested `material::capsule` owner and update callers to that owner; do not keep old flat names as transitional aliases. Preserve private freshly reconstructed physical state and documentary evidence as separate concepts.

`legacy_formats.rs`, `collagene_xml_compat.rs`, and retained source bytes are intentional external-file support. Their names do not justify deletion. Rename them to clear external-format owners if useful, while preserving supported original Collagene/XDNA imports and fidelity checks. This is different from the internal `model.rs` compatibility facade.

## Integration-test discovery and organization

`crates/core/Cargo.toml` currently uses default Cargo integration-test discovery; it has no `autotests = false` and no `[[test]]` entries. The 320 direct test files are therefore separate test executables. Simply moving them to arbitrary nested directories would silently stop discovery. Existing helpers are explicitly included with relative `#[path]` attributes, for example `tests/actual_strand_protocol.rs:12` and `tests/protocol_quantities_transport.rs:1`; all such paths need updating when their owning test moves.

Preferred result: a small set of domain integration-test harnesses, each loading nested test modules through explicit module declarations. Keep no more than 20 direct Rust files at every level, including `tests/` and helper directories. For example:

```text
tests/
  cloning.rs                # registers cloning/{restriction,pcr,overlap,...}
  planning.rs               # registers intent/editing/combinatorial modules
  verification.rs           # registers invariants and adversarial cases
  formats.rs                # registers formats/{sequence,annotations,...}
  material.rs
  workspace.rs
  review.rs
  transport.rs
  output.rs
  common/                   # eight existing helpers, kept cohesive
  cloning/
    pcr.rs                  # declares pcr/<scenario>.rs modules
    pcr/...
    restriction.rs
    restriction/...
```

Use actual module trees, not `include!` fragments. Set `autotests = false` and explicitly register these current harnesses with `[[test]]` if that makes the intended discovery surface unambiguous. Alternatively retain default discovery with only the intended root harnesses. Choose one coherent strategy for the batch; do not retain 320 obsolete target aliases. A second valid design is explicit per-file `[[test]]` entries preserving 320 executables, but it adds manifest/compilation overhead with little ownership benefit here.

Record a static before/after inventory of `#[test]` and parameterized/generated tests, including fixture paths, before moving them. Compiler discovery and `--list` output from the final harnesses must confirm that scientific coverage survived the regrouping. Update targeted command references from former executable names to domain target plus test filter. CI's package-wide `cargo test` in `.github/workflows/check.yml:29` will discover registered new harnesses, but cannot detect tests silently omitted by bad module declarations. Preserve test semantics, assertions, bounds, and current binary-format validation; this audit found no basis to delete scientific tests.

## Ordered implementation batches

1. **Ownership map and policy gate.** Record every source/test destination and owner, with no duplicate destinations and a directory-count projection. Add a lightweight repository-wide gate: `.rs` files at most 5,000 lines; warning at 3,000; at most 20 direct Rust files per directory. Include tests, examples, benches, build scripts and generated Rust unless the user explicitly approves an exception. Keep warnings actionable rather than treating 3,000 as the hard limit.
2. **Worker and model cleanup.** Split the oversized transport, delete the model facade, and update every core/web/MCP/source/desktop/test consumer directly. Keep current binary bytes unchanged for mechanical edits. Inspect numeric selector tables and enum matches statically before compiling.
3. **Plan/verifier/output owners.** Move plan declarations and natural methods, introduce private execution state, extract coherent operation handlers and checked protocol rendering. Update all recipe, UI, MCP, export, and fixture callers together. Preserve one internal plan representation and fresh verification authority.
4. **Scientific subtree relocation.** Move independent domain families in parallel with disjoint source ownership; update central module declarations/imports in one coordinated integration batch. Co-locate material preview types and methods. Delete obsolete files and paths immediately; no compatibility facade, delegator, fallback decoder, or temporary glob import.
5. **Test harness regrouping.** Move domain tests and shared fixtures with registration, update helper paths and focused command documentation, and compare discovery coverage. This can be prepared alongside source work, but a single integrator should own Cargo and central declarations.
6. **One coherent validation batch.** After implementation and caller updates, run the focused native checks/tests and WASM check using the shared cache and existing package/feature selection. Compile selected integration harnesses together with `--keep-going` or the repository's combined no-run JSON build when the web harness is included; execute resulting harnesses directly where instructed. Collect failures, fix them as a batch, and avoid repeating unchanged passing checks. Required compiler work runs in the background while independent implementation continues. Full platform CI is a milestone check.

## Using rust_refactor without violating repository rules

Reviewed helper documents: `README.md`, `docs/MODULE_REFACTOR_WORKFLOW.md`, `docs/MODULE_SIZE_CHECK.md`, `docs/EXTRACT_ITEMS_MODULE.md`, and `docs/RELOCATE_MODULE_TREE.md` under `/home/mahogny/github/claude/rust_refactor`.

Useful audit candidates are `struct-impl-scatter-stats`, `struct-module-stats`, `to-oop-stats`, `module-coupling-stats`, `visibility-stats`, and reviewed relocation dry runs. Semantic analysis is useful for ambiguous ownership and cfg-dependent callers; raw free-function count alone is not a reason to create methods. Native/WASM profiles must match the owner being analyzed.

Important conflicts: the helper README and workflow say all tool-launched compiler verification uses the release profile. Current repository rules disallow automatic release verification and require coherent native/WASM batches. `extract-items-module` also deliberately preserves old paths and inserts `use super::*`; the workflow suggests compatibility delegators/glob phases. Do not use those write modes as supplied. Use reviewed reports/dry-run plans and direct source edits, or adjust the helper's behavior explicitly first. `relocate-module-tree` supports a manifest with compatibility policy `none`, but its automatic release verification still conflicts. No helper binary was built or run during this audit.

The helper's `module-size-check` excludes tests/examples/benches and optionally generated code. It therefore cannot enforce the user's universal file limits by itself. A separate simple all-Rust inventory must enforce both directory counts and file sizes.

## Acceptance criteria

- Every repository directory has at most 20 direct `.rs` files; core aims for fewer. Every Rust file is at most 5,000 lines, and every file at/above 3,000 has an explicit reviewed extraction plan.
- `protocol.rs` and `verification.rs` monoliths are replaced by the coherent current owners; scientific algorithms, evidence bounds and current transport semantics are preserved.
- Structs/enums and their core inherent implementations are together. Remaining cross-file implementations have a concrete ownership/trait/target justification.
- No internal model facade or superseded module aliases remain. All callers use current paths; obsolete internal formats remain rejected.
- Scientific test coverage is preserved and discoverable under the new harnesses, including intended versus alternative products, physical strand geometry, fresh verification authority, revision binding, and malformed/oversized payload rejection.
- Focused native tests and WASM checks pass after the complete batch; cross-platform CI uses the final current layout. No browser checks are needed for this structural audit.
