# Targeted ownership implementation

This records implementation of the next batch from the
[current refactor plan](REFACTOR_PLAN.md). It is separate from the preceding
annotation semantic change and does not establish full-product acceptance.

## Implemented owners

- `Design` owns combination identity, construct validation and bounded matrix
  size/paging, together with its private selection/shape helpers.
  `SelectedMatrix` owns validation. Ten original core/native bodies match after
  receiver and direct caller substitutions, including the two TreeJob operations.
- `TreeJob` owns node admission and remote scan progression in its existing module.
  Repository transaction/checkpoint authority stays on `Repository`.
- Five Draft operations now belong to its existing implementation: digest
  fragment installation, input/edit request construction, orientation, and
  material-origin assignment. All five bodies match after receiver substitutions;
  52 test names in the touched web source files remain unchanged.
- `IntendedReference` owns context and warning accessors. The draft model's
  `MaterialDomain` forwarding export is removed; digest consumers use the core
  declaration directly. Preparation consumers use the canonical preparation owner.
- The annotation editor's location parser/formatter forwarding functions are
  deleted; all source and test callers use the canonical format codec directly.

No organization-only change intentionally alters serialized bytes, deterministic
identities, scientific algorithms, current schemas or annotation engine stamps.
The separate qualifier batch advances scientific annotation stamps as documented
in its [profile](../COORDINATE_QUALIFIERS.md).

The preparation panel is now a 289-line composition owner. Its dedicated
1,472-line review owner contains proposal declarations, request construction,
current-context admission, explicit apply methods and the shared worker lifecycle.
The 711-line pool owner keeps context admission/input binding with receipts,
settings, candidate review and selection. A 110-line diagnostics module owns
presentation helpers; existing preparation compilation owns shared recipe builders.
The panel's `prepared_rows` forwarding export is deleted.

An independent source comparison accounts for 27 moved function/DTO bodies,
including complete pool execution/rendering, request builders and context types.
Proposal navigation/apply/export and evidence finalization retain their distinct
source/revision/target/options/epoch checks. The author also compared all original
string literals and the shared request bodies. Static evidence is retained in the
[comparison record](ownership-implementation-evidence.json); it supplements rather
than replaces compilation and execution. No tests or harnesses were relocated.

The layout gate passes 1,241 Rust files in 214 directories with no violations or
3,000-line warnings. The final web/worker WASM check passes after correcting one module path and a
missing panel type import. Native harness compilation passes. Discovery preserves all 2,451 baseline cases
across the same 27 harnesses and adds 12 repository cases in two selected harnesses
(2,463 total / 29 harnesses). Focused execution passes 495 unique cases: 64 core, 18 MCP SDK, 407 web model
and six repository folder cases. The preceding annotation batch passes 2,152
selected cases; these counts overlap and must not be added as distinct cases.
Formatting, whitespace and documentation links pass. Exact-source CI remains
required for this combined checkpoint. The initial native attempt exhausted disk
space; 2.3 GB of superseded Collagene test executables were removed while current
artifacts were retained, and the same build selection then passed. No browser/application
HTTP, screenshot or automatic release verification is included.
