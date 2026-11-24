# collagene-map-label-layout

MIT, dependency-free Rust placement of measured label boxes. No text measurement,
glyph slicing, sequence model, biological coordinates, rendering or runtime I/O.
The caller reuses cached text measurements while updating projected anchors during
rotation, zoom and resize. This is a display layout, not scientific validation.

```rust
use collagene_map_label_layout::*;
let labels=[Label{id:7,stable_rank:42,anchor:Point{x:80.,y:50.},
    preferred_center:Point{x:30.,y:50.},side:Side::Left,
    measured:Size{width:40.,height:14.},priority:0}];
let columns=[Column{side:Side::Left,rect:Rect{x:5.,y:5.,width:50.,height:90.}}];
let result=layout(&labels,Viewport{width:100.,height:100.,padding:5.},&columns,&Options::default())?;
assert_eq!(result.displayed,1);
assert_eq!(result.labels[0].id,7);
# Ok::<(),Error>(())
```

The API takes borrowed labels with stable unique integer IDs, projected anchors,
preferred centers/sides, a camera-independent stable rank, caller-measured full text-box dimensions and an explicit
presentation priority. Left/Right regions pack vertically (circular map sides);
Above/Below pack horizontally (linear map sides). Regions must be disjoint and
inside the supplied padded viewport. A caller may use any subset of four sides.
The caller chooses map exclusion regions; the library has no knowledge of map
arcs, tracks, titles, source positions or font baselines. Higher priority means
presentation preference only, such as selection, never biological importance.

All input IDs appear exactly once in original order. Each result is either a
placed measured box + leader with the EXACT original anchor, or explicit overflow
with a reason. Reasons distinguish display-budget LOD, missing region, too wide,
too tall, dense packing and floating-point precision failure. Available region dimensions accompany size/density
overflow for explicit caller remeasurement. Text is never silently clipped or
truncated: fitting ellipsis requires the caller to measure actual glyphs and
resubmit. Fitting does not change the measured box dimensions. A final linear boundary /
adjacent-collision check never publishes a box outside its region or overlapping
another box: finite dimensions that collapse at the projected coordinate's
precision, or cannot survive exact-fit floating-point rounding, receive explicit
Precision overflow with available bounds. No epsilon tolerance silently converts
such failures into placed labels.

## Display budget / LOD

Defaults: 32768 inputs, 128 considered/displayed labels maximum, 4 coordinate-unit
gap. Hard caps: 131072 inputs, 1024 considered labels, 4 side regions, 64MiB managed
storage charge. These are explicit display limits, not feature/biology limits.
ALL features must remain in the caller's independent scene and inspector. A
label overflow must not remove a source feature or selectable geometry.

Priorities are considered high to low. If a priority tier exceeds the remaining
budget, its labels are sorted by caller-provided stable rank and stable ID; midpoint ranks from evenly spaced strata enter the
budget. Thus dense 10k inputs are sampled across the caller-provided ranking,
not cut off at the first feature in a file. A reference/spatial rank is a caller
choice; arbitrary ranks do not prove geometric or biological coverage. Exact ties use IDs deterministically. This is not an importance
or coverage score. Every unconsidered identity has DisplayBudget overflow.
The budget-selected subset is invariant under input permutation and camera
projection changes when IDs, stable ranks and priorities stay unchanged. Stable
rank is caller data (for example whole-reference ordering), not a projection
computed by this crate. Projection changes can change positions, side capacity
failures and final placement, while preserving the budget-selected identities. Stable IDs make those changes inspectable, but no
temporal hysteresis or optimal movement is promised.

Within a side region, eligible selected labels are admitted greedily in the above
priority/stable-rank order while their total extent+gaps fits. Overflow does not trigger
an unbounded search or budget backfill. Admitted boxes are ordered along the flow
axis, swept forward near preferences then backward to fit the region. Boxes in a
region cannot overlap; disjoint regions ensure boxes on different sides cannot
overlap either. Leaders may cross other leaders or map content. No claim of
optimal coloring, maximum label count, shortest leaders or optimal displacement.

## Bounds and portability

Allocation-free preflight checks counts/options, finite values, region geometry,
viewport/padding and a checked storage charge. Positive viewport/box dimensions
are at most 1e6 units; projected anchors/preferences have absolute value at most 1e9
and may lie outside the viewport (zoomed geometry), but remain unchanged in leader
evidence. NaN/infinity/zero/negative dimensions fail explicitly. Duplicate IDs are
detected with a bounded sorted index after resource preflight.

Output/index/selected-list storage is charged before allocation; allocator
bookkeeping and caller-owned input are outside this managed-storage accounting.
Rust `sort_unstable` has documented worst-case O(N logN) and in-place allocation
behavior. All sort keys are fixed-size finite numeric values. Total time is
O(N logN+B logB), memory O(N+B), for N inputs and B explicit display budget;
four-column lookup is constant. No all-pairs collision scan, recursion, external
dependencies or per-base arrays. Native and wasm32 use the same pure code.

## Checks

Independent literal and exhaustive small-grid collision/boundary/anchor tests,
input-permutation and invalid-value tests, dense 10k and sparse 10k scenarios.

```sh
cargo test --locked --offline -p collagene-map-label-layout -j2
cargo check --locked --offline -p collagene-map-label-layout --target wasm32-unknown-unknown -j2
```
