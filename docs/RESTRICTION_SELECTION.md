# Restriction regions and named enzyme sets

Restriction discovery uses the existing source-bound scan. Changing a region,
relation, geometry basis or count filter derives a view of that evidence; it does
not rescan sequence or construct digest products. The same derived sites/counts
drive the enzyme picker, sequence track and circular/linear map.

Choose a nonempty exact range, joined annotation or origin-spanning selection,
then activate a region relation. This pins the selected union to the document
context and revision. Site navigation keeps the pin; **Use current selection**
replaces it explicitly. Changing the source invalidates the pin. A missing valid
pin produces unknown scoped counts and no displayed sites, rather than silently
showing whole-source results.

Recognition classification tests all recognition bases: fully covered is inside,
disjoint is outside, partial coverage crosses the boundary. Cut classification
tests both strand-cut boundaries, including selected endpoints: both inside,
both outside, or mixed. Circular boundary zero and sequence length are equivalent.
Type IIS recognition and cleavage can therefore classify differently. Missing,
outside-substrate or ambiguous cuts are undetermined. This is geometry evidence,
not an experimental cleavage claim.

Possible recognition remains possible; incomplete scans remain incomplete after
filtering. Unclassifiable cuts prevent exact inside/outside/crossing count claims
for the affected enzyme. Exact-count filters exclude uncertain counts. Site indices
remain the original scan indices so navigation retains its identity.

Named sets save ordered enzyme choices as bounded native-workspace attachments,
with catalog edition/content hash, stable catalog IDs and complete recognition/cut
geometry. Save replaces an existing name case-insensitively; empty sets are valid.
Loading refuses missing enzymes or changed geometry rather than silently selecting
same-name substitutes. A catalog revision change with identical definitions is
reported and can load. Save/delete change attachments only; they do not change
sequence history or assert that a digest was checked.

Implementation: `core::restriction_selection`, `core::enzyme_sets`, and
`web::restriction_view`. Native fixtures cover exact-location extraction, interval
unions, circular/linear geometry, uncertain evidence, catalog drift, archive
roundtrip, attachment integrity and atomic bounds. Browser acceptance remains
explicitly deferred until the full in-scope implementation is complete.
