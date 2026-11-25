# Select an actual digested fragment by boundaries

`core::digest_cut_selection` connects a user's explicit boundary choices to the
existing catalog inspection/product engines. It does not introduce digestion
biology, a target, a new plan operation or a reference slicing shortcut.

`Boundary` is `SourceStart`, `SourceEnd` or `Cut {site_index}`. A cut index identifies
one actual site instance in the exact fresh inspection, not an enzyme name or the
first occurrence of a motif. `Choice {first,second,complement}` records the click
choices even when normalization selects the same physical fragment.

`preview(catalog,selection,source_role,source,options,expected_batch_revision,choice)`
returns `Report {context_revision,choice,inspection,product}`. Both nested objects
are the existing full `CatalogRestrictionInspection` and `CatalogDigestProduct`.
The chosen index is `report.product.product.fragment_index`; its actual global
source top/bottom geometry is `report.product.product.summary`.
`selection_revision` accepts the same arguments for worker binding.

## Freshness and physical boundaries

Borrowed native-model and existing catalog/source guards run before hashing or
copying the source. The expected outer catalog batch must equal the current batch
for the exact source role/revision, catalog edition/content/ordered selection and
inspection options. A fresh full catalog inspection must then have complete
recognition and determined physical products with no exhausted limits.
Possible recognition, missing/out-of-substrate cuts, diagnostic strand ambiguity
and absent site instances cannot supply physical endpoints. Imported physical
summaries cannot restore a source; the unchanged engines apply their normal
substrate restrictions and assumed-complete cleavage model.

Boundary semantics are:

- **Circular:** first → second follows source-forward coordinates, crossing the
  origin when needed. Complement swaps the endpoints. Choosing the same actual
  cut at both ends requests one full source turn; this succeeds only when one
  full-turn fragment actually exists. Coincident cut instances are resolved
  against the determined duplex fragment set, not arbitrarily merged by labels.
- **Linear:** SourceStart/SourceEnd mean the actual intact reference's free ends.
  The two endpoints normalize into ascending source order. Equal endpoints are
  empty and refused. Complement is refused because it would describe disconnected
  pieces rather than one fragment.

The selected top interval **and bottom interval** must match exactly one actual
linear fragment. Circular bottom boundaries use the selected site's signed
strand-cut stagger in the unrolled global source coordinate system. Global bottom
start is distinct from the fragment's local `duplex.bottom_start`; both actual
strand lengths remain in the returned summary. If an internal cut splits the
requested path, selection fails with an additional-cuts explanation. A matching
reference substring is never substituted for the missing physical fragment.

`catalog_digest::product` freshly computes the actual selected product through the
existing engine. Its source/options/geometry/inspection pins and exact summary
must agree with the selected fresh fragment before publication. Free-ended linear
pieces retain their unknown inherited terminal groups; selection does not make
those groups phosphorylated or accessible. Methylation, supplier conditions and
experimental cleavage remain the existing explicit assumptions/limitations.

## Binding and bounds

The context hash binds the cut choice to the verified outer batch, which already
binds original source identity, selection geometry and options. Changing source,
role, options, catalog or endpoints invalidates the receipt. `Report::validate`
and `validate_binding` perform readonly shape/context and nested receipt checks;
they do not rerun digestion or restore private physical state from decoded data.

Existing selection/site/fragment/work/output limits apply. Site indices are capped
before indexing, model/location bounds precede source serialization, and the
combined report is bounded by the ordinary 32 MiB canonical envelope. Inspection
and selected-product creation each use their existing fresh bounded paths; this
adapter does not trust a caller-supplied inspection instead of recomputing it.
Original molecules and source attachments are not mutated.

Literal fixtures cover two-cut circular interval/complement, single-cut full-turn
linearization, linear free-end fragments with unequal actual strand lengths,
additional-cut refusal, possible/out-of-range sites, stale batches, immutable
sources and changed choice/product receipts. For example, circular
`GAATTCAAAAGAATTCCCC` cut at its two EcoRI instances yields first selected top
`AATTCAAAAG` (top start1, bottom start5, lengths10/10), or complementary top
`AATTCCCCG` (top start11, bottom start15). These are concrete existing fragment
records, not stitched reference intervals.
