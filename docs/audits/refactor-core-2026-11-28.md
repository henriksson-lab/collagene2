# Current scientific ownership audit

28 November 2026. This is a fresh static review of the current working tree,
including the coordinate-qualifier batch. It covers core/model, scientific
auxiliary crates and the reusable layout kernels. Web rendering/controller
ownership and MCP/native services have separate reviewers. No production source
was changed and no compiler, test, browser or application HTTP check was run.

## Measured limits

`python3 scripts/check_rust_layout.py` reports **1,238 Rust files in 212
directories, zero violations and zero size warnings**. Core currently has 215
source and 371 integration-test Rust files. The largest core source is
`crates/core/src/review/contigs/mod.rs` at 2,252 lines; the largest auxiliary
scientific source reviewed is `crates/topo/src/lib.rs` at 2,315. Neither requires
an artificial split to pass the limits. The old flat-directory, worker monolith,
model facade and material-preview scatter findings have already been addressed;
the earlier audit tables are historical evidence, not current defects.

The highest-pressure core directories are tests: sequence-format tests have 17
direct Rust files; read-analysis, material-derivation and protein-sequence tests
have 16 each. Annotation tests now have 15. The next related additions should
first create registered semantic groups. Preserve `autotests = false` and the
explicit domain harnesses in `crates/core/Cargo.toml`; compare before/after test
inventories whenever moving these modules.

## Ranked remaining work

| Order | Current evidence | Refactor decision and scope |
| --- | --- | --- |
| 1 | `planning/combinatorial/mod.rs:987`, `:1107`, `:1248`, `:1259` expose `combination_id`, `validate_construct`, `matrix_size`, `matrix_page` as functions whose dominant owner is `Design`; `:1299` validates only `SelectedMatrix`. `Design` is declared at `:68` with its existing implementation at `:162`. | Complete the existing ownership conversion: put these operations on `Design` and `SelectedMatrix` in this same owner module. Update core/UI/MCP callers directly. Preserve mixed-radix order, bounded paging, IDs and validation bytes. Keep `cartesian_size`/`cartesian_page` free: they operate on generic counts, not design state. This is a small coherent batch, not grounds to split the 2,102-line module into many files. |
| 2 | `review/contigs/mod.rs:223` validates/hashes one `ContigOptions`; `:1902` validates/hashes one `ContigReport`. Both types and their other methods are already local (`:126`, `:1936`). | On the next contig edit, make these `ContigOptions::revision` and `ContigReport::revision` and update callers directly. Keep `batch_revision(reads, options)` free unless introducing an actual owned batch is independently justified. Do not invent a wrapper solely to remove a free function. |
| 3 | `model/src/lib.rs:426`–`:788` contains GenBank edit-log declarations, bounded replay and citation-coordinate remapping between `Molecule`'s declaration (`:413`) and implementation (`:790`). Core has the canonical GenBank location codec in `formats/sequence/location.rs` and coordinate-bearing qualifier projection in `sequence/annotations/qualifiers.rs`. | Before adding native edit support for those qualifiers, explicitly choose the lower-level annotation/coordinate owner. Move complete codec/qualifier concepts downward only if model editing needs them, and remove the old definitions/callers in the same batch. Keep `Molecule` and its meaningful edit implementation together; format-specific edit-log/citation mechanics can have a cohesive `model::genbank_edit` owner. Do not make model depend on core, duplicate the parser, or distribute `impl Molecule` across format files. The current 976-line model file itself is comfortably within limits; this is a dependency/ownership decision, not an urgent size split. |
| 4 | `core/tests/formats/sequences/mod.rs` directly registers 16 sibling scenario files, giving 17 files with the registrar. `snapgene`, `snapgene_export`, `snapgene_protein` and `oligo_formats`, `oligo_orders` are obvious families. | Before expanding those families, group them under explicitly registered `snapgene/` and `oligos/` domains. Maintain exactly one harness registration for every test. Do not move files merely to make a count lower without retaining discovery. |

Paths in the table are relative to `crates/core/src` unless explicitly prefixed
with a crate name. These are organization findings, not claims of scientific
incorrectness. No hard limit currently forces a core refactor before other work.

## Owners worth preserving

- Worker execution now has a small admission entry, an exhaustive router and
  typed domain handlers under `core/src/protocol/execution/`. Those handlers are
  stateless algorithms over concrete payloads. Turning them into implementations
  scattered around one universal context type would worsen ownership.
- `core/src/material/derivation/{digest,pcr,end_repair,terminal,overlap_extension,annealing}.rs`
  keep preview declarations, computation, revision and validation together. The
  replay/install authority remains in the 1,611-line parent. Do not reintroduce
  DTO-only files or move record authority into scientific preview owners.
- `core/src/verification/engine.rs:359` owns correlated execution state and its
  implementation; `planning/plan.rs` owns `Operation` and `DesignPlan` with their
  methods; checked protocol rendering has its own freshness-bound owner. These
  are substantive domain boundaries, not compatibility facades.
- `core/src/review/contigs/mod.rs` has real private `Prepared` and `Layout`
  owners. Its affine-overlap helpers are stateless/symmetric and need not become
  methods. If growth reaches 3,000, review extraction of the complete overlap
  kernel or layout owner with its types; do not split report methods away from
  `ContigReport` to meet the count.
- `topo/src/lib.rs` keeps private `Graph` (`:537`) and `ExecutionState` (`:1110`)
  with implementations. `topo/src/projection.rs:277` defines a different private
  `Graph`; this is not a scattered implementation of the simulation graph.
  `projection`, `processing` and `nick_sealing` are different physical proof and
  lifecycle responsibilities. Preserve that distinction.
- Recombination's `lib.rs` and `multisite.rs` define different `Profile`/`Report`
  types. Likewise, core has many distinct module-local `Report`, `Options` and
  `Request` types. Name-only searches cannot establish type scatter and must not
  be used as evidence for merging these owners.
- `enzyme-catalog/src/lib.rs:116`/`:225` colocate `Catalog` and its implementation.
  `view-layout/src/lib.rs` colocates `Interval`/`Layout` methods; `pack` is a
  stateless packing algorithm. `map-label-layout/src/lib.rs` colocates `Rect`
  and `Side` methods; radial label placement lives under `radial.rs` with its
  private `FanBlock` implementation. Geometry routines such as `layout`,
  `fan_out` and `arrow_head` are appropriate free functions. Web-specific drawing
  belongs to the web audit, not these reusable geometry kernels.
- Auxiliary read-alignment codecs place readers/writers and their implementations
  together (`sam.rs`, `bam.rs`, `bgzf.rs`); structures keeps `Structure` methods
  in `lib.rs`, CIF parser lifecycle in `cif.rs` and scene generation in `scene.rs`.
  No structural finding justifies restarting deferred BAM/BAI feature work.

## Tool and verification policy

The helper `/home/mahogny/github/claude/rust_refactor` offers useful read-only
semantic ownership reports. Its write workflow automatically invokes release
checks, which conflicts with the repository verification policy; no such mode
was used. This review uses source and call-shape evidence, not a claim that a
semantic scatter analysis was executed.

Implement owner and caller changes as a single coherent batch, retain one
current representation, remove superseded paths and use no forwarding aliases.
After source work is complete, run formatting/layout checks, the affected native
harnesses and WASM compilation with the shared cache. Compiler work should run
in the background while independent implementation continues. Documentation-only
audit changes do not require a new compiler or scientific test run.
