# collagene-view-layout

MIT, dependency-free Rust library for static annotation lane layout. It does not
read sequence data, draw pixels, validate biological annotations or optimize graph
coloring. No runtime I/O or native-only dependency; native and WASM compatible.

```rust
use collagene_view_layout::{pack, Item, Interval, Limits, Placement};
let joined = Item { id: 7, spans: vec![
    Interval { start: 90, end: 100 }, Interval { start: 0, end: 10 }
] };
let gap = Item { id: 9, spans: vec![Interval { start: 10, end: 90 }] };
let layout = pack(&[joined, gap], 100, &Limits::default())?;
assert_eq!(layout.lane_count, 1);
assert!(layout.assignments.iter().all(|a| a.placement == Placement::Lane(0)));
# Ok::<(), collagene_view_layout::Error>(())
```

`pack(&[Item], source_length: u64, &Limits)` borrows the caller's data. Each item
has a unique stable `u64` ID and at least one nonempty half-open source interval.
Intervals must satisfy `start < end <= source_length`. Use
`Interval::from_start_length` when constructing from lengths to reject overflow;
already-wrapped external arithmetic cannot be reconstructed by this library.

Own overlapping/adjacent spans are sorted and merged. Gaps are preserved. Items
are processed by normalized interval-vector lexicographic order, then stable ID.
Every span of a joined or circular-split item must fit the same lane. First-fit
selects the lowest nonconflicting lane; adjacent items can share a lane. This is
deterministic across input permutations, including which items overflow, but can
be less compact than an optimal coloring. Changing the item set can change lanes.

Assignments preserve original input order, with every ID exactly once and its
normalized spans plus `Placement::Lane(index)` or `Overflow`. `overflow_ids()` is
an iterator over those retained assignments, not a second stored ID collection.
`lane_count` is actual used lanes, zero for empty input. Zero lane capacity is
valid and retains every item as overflow. Empty input permits source length zero.

## Resource policy

Limits are checked before any input clone or library heap allocation. Absolute
caps: 131072 items, 524288 original intervals, 256 lanes and 2,000,000,000
conservative work units. Defaults: 32768 items, 131072 intervals, 64 lanes and
128,000,000 work units. A fixed 256MiB managed-storage charge cap also applies.
Thus inputs below count caps may still fail work/storage budgets. Counts refer to
original input, including duplicate/overlapping spans before normalization.
Invalid limits fail explicitly rather than silently clamp. Default limits cover
10,000 sparse annotations on a 1Mb source without truncation or fallback.

`work_bound` performs allocation-free count and storage preflight. Sorting uses
explicit bottom-up mergesort with one reserved scratch buffer, avoiding quadratic
insertion sorting and implementation-dependent sorting bounds. Duplicate IDs are
checked in a sorted index. Each merge level of canonical item sorting consumes
every key once; a lexicographic span-vector comparison costs no more than the
consumed key's length plus one, so that pass costs O((N+S)logN), not O(N*S).

Each lane has a disjoint `BTreeMap<start,end>`. Only predecessor/successor intervals
can conflict with a candidate span. ALL normalized spans of a joined item are
checked before ANY are inserted. Failed attempts leave lane indexes untouched.
At most2*L*S neighbor queries and S insertions are charged before allocation, each
at 4*(ceil(log2(S+1))+1) logical tree-work units. Other charges cover validation,
merge comparisons/copies and N*L lane visitation. Work units are conservative
algorithmic accounting, not CPU instructions or elapsed-time guarantees.
The full formula is in `work_bound`; checked arithmetic rejects overflow.

`memory_charge` accounts for assignment storage, all original-sized normalized
span vectors, a largest-vector normalization scratch buffer, item ordering plus
order scratch, lane headers, and 256 bytes per possible B-tree interval entry.
Caller-owned input, allocator bookkeeping and process/runtime overhead are not
managed-storage charges; this is not an OS memory sandbox. The B-tree charge is a
conservative container allowance rather than an allocator-specific size promise.
B-tree nodes allocate incrementally within the preflighted entry budget (the
standard container has no reservation API). No per-base arrays, randomness,
recursion or unbounded search. Memory is O(N+S+L); work is
O((N+S)logN + S*log(max spans per item) + L*S*log(S+1)). Errors publish no partial
layout. Bounds validation occurs before allocation; duplicate-ID detection then
uses reserved sorted-index scratch rather than quadratic borrowed comparisons.

## Future viewer integration contract

Compute a whole-source, source-revision/options-bound plan once, with independent
category budgets if desired. Keep that plan and row height stable throughout
asynchronous viewport loading. Use actual `lane_count` to remove reserved empty
rows; expose overflow items in an inspectable list. Never re-pack only the latest
viewport response, use modulo assignment, drop overflowing identities or infer a
biological relationship from lane membership. Pixel offsets, scroll anchoring,
source revision identity and cancellation remain the host's responsibility.

## Verification

Independent literal joined/circular/adjacent/overflow/boundary fixtures, exact
resource-boundary tests and an exhaustive four-position / three-item bitset oracle
check first-fit outcomes, span unions, collisions, coverage, ID order and reversed
input invariance. Independent 10,000-item sparse and mixed/overflow fixtures
exercise realistic annotation scale; a rejected-joined-item fixture checks atomic
index updates. The per-position oracle exists only in tests.

```sh
cargo test -p collagene-view-layout --locked --offline -j2
cargo check -p collagene-view-layout --locked --offline --target wasm32-unknown-unknown -j2
```
