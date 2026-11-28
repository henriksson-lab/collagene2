# Named digest setups and predicted gels

A digest setup now retains its name, exact source document/event/revision, complete
ordered catalog selection, inspection options, optional supplier reaction setup
and historical inspection/condition reports in a content-hashed native-workspace
attachment (current kind/schema 3). The declared model assumes complete
ideal cleavage with methylation unassessed. Separately recorded unit-assay
conditions do not establish partial-digestion kinetics, star activity or
experimental cutting. The collapsed Reaction conditions control offers exact
formulation choices, buffer/temperature/time, an activity threshold, explicit
contexts and a manual worker assessment. No formulation is chosen implicitly.

Save requires a current source-bound inspection. Incomplete evidence can be kept
with its original status, without inventing products. Reopen displays source,
catalog/options, recognition/cut diagnostics and paged fragment evidence separately
from the current digest. Saved raw JSON remains downloadable. Loading/listing does
not rescan or restore a current-success flag. Catalog/engine changes remain visible
as historical evidence; fresh reuse must resolve the exact pinned catalog.

Use settings applies the full saved pin and limits only to the original exact
active event. It never restores enzyme names and silently selects different
geometry. Explicit enzyme changes clear the pin. An explicit source-fork action
creates a new document from the retained version, preserving original heads. The
returned fork must contain the same revision and be the actual new active event
before its pin/limits are reused. The saved record remains unchanged, no fresh
result is restored, and Run virtual digest is still required.

Both current and saved setups can be added directly to predicted gel lanes without
retyping enzymes. Current inputs are guarded against source/options changes; saved
recipes carry their original event and catalog/limits. Lanes re-resolve that event
for every request and asynchronous completion, so later source-head changes never
redirect them to a current index or equal-DNA alias. Multiple additions survive
navigation between inspectors; actual workspace restoration clears old lane and
diagnostic source selections.

Pinned gel prediction freshly inspects the actual retained molecule using the
full saved recipe. Only complete scans with determined products supply fragment
lengths. Unresolved/partial cleavage or more than 1,024 fragments yields an explicit
withheld/error lane, without truncation or silently changing saved limits. Actual
top/bottom strand lengths determine mean duplex size and mass allocation. Circular
uncut material remains unpositioned under the linear-duplex migration model.
Existing illustrative calibration, co-migration, mass and uncertainty rules remain.
Full end evidence stays in the digest inspection; gel report DTOs retain their
existing layout.

Pinned gel JSON includes the exact submitted recipes, lane indices, labels, setup
names and masses alongside the result, not a later read of mutable controls.
Legacy-only gel JSON and SVG remain supported. No observed images or OpenGel
comparison/formats are introduced. All reconstruction is embedded native/WASM Rust.

MCP `predict_gel` requires a native workspace reference on each catalog-pinned
lane. It verifies the exact retained source owner/event/revision and actual
molecule against the immutable document reference, including historical events
after head advance. A matching DNA hash or caller-supplied event alone is not
ownership evidence. Restores happen outside the document-store lock; full declared
lane recipes and references remain in analysis artifacts. Legacy tool inputs
without workspace references remain supported, and prediction grants no writes.

Implementation: `core::digest_setups`, `gel::Preparation::CatalogDigest`, append-only
worker55 job124/result116, the digest inspector and typed gel-lane bridge. Native
fixtures cover historical ownership/catalog drift, immutable archive roundtrips,
nonstarter enzymes, staggered strands/mass, circular uncut DNA, partial/oversized
withholding, binary compatibility and unchanged legacy annotation capacity.
Browser acceptance remains deferred until full in-scope implementation is complete.


## Retained reaction conditions

A saved setup with explicit reaction options also retains their freshly recomputed
catalogue assessment. Changing a condition changes its documentary identity and
invalidates displayed current evidence; it does not change assumed cleavage
geometry. Missing choice/facts, conflict, required context, below-threshold activity
and temperature mismatch remain distinct. Common-condition suggestions use only
exact recorded buffer/temperature and assay contexts shared by all chosen
formulations. Star risk and duration adequacy remain unassessed even when the
recorded activity meets the threshold. Qualifiers and citations remain visible.

Loading saved settings also restores or explicitly clears the optional conditions;
reopening the historical report never restores a current assessment. Detached
GenBank/EMBL exports retain bounded condition declarations as unverified producer
comments. Keep the native workspace for exact source and report provenance.
