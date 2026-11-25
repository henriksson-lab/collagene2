# PCR source overlays

The applied exact-PCR fields, captured pair domains and selected exact/mismatch
PCR product now connect to the existing sequence and map views. Overlays are
transient display projections; they do not add features, alter DNA, create
products, change an independent target or authorize a cloning plan.

Select an exact-PCR product or a mismatch-PCR geometry to publish its original
source intervals, oriented primer domains, amplified interval and applicable
mismatch/ambiguity positions. A sole result retains the existing auto-selection
policy; multiple products require an explicit choice. Summaries preserve search
and enumeration limitations. Showing one geometry does not establish uniqueness,
physical end chemistry, experimental amplification or target verification.

The baseline exact-PCR result stores role names without oligo indices. The overlay
therefore checks at most four role assignments against actual source domains,
selected product letters and tails. Names are not assumed unique, and one oligo
may occupy both roles. Compatible assignments remain visible with ambiguity
warnings; the projection does not rescan the template or silently choose a winner.
Mismatch-PCR uses its selected indexed binding geometry, including reverse-strand
5′ mismatch coordinates and circular intervals.

A weaker fallback shows matching captured pair anchors as selected-reference
comparisons. Actual original event/bases and current oligos must still agree.
These marks are explicitly labeled as captured selections, not binding searches;
they contain no inferred amplified interval. Unanchored or changed pairs have no
invented placement. Synthetic 5′ tails are always listed separately as added DNA,
never mapped onto source bases.

Sequence drawing adds two transient primer lanes while an overlay is present,
with directional arrows, amplified-source shading and diagnostic source-base
marks. The existing reflow anchor preserves the visible reference position.
Circular/linear maps share role colors and source-span navigation; annotations and
restriction labels keep their existing hit priority. Source IDs are not replaced
with fabricated feature indices. DNA hashing and primer search do not occur during
painting or scrolling.

View-template/primer actions navigate the same geometry. When its exact source is
already active, all source parts are highlighted in the main editor. A different
source opens a read-only snapshot using the shared viewer, without changing the
independent target or source heads. The snapshot keeps the selected overlay only
while its original owner, exact source event and workspace epoch remain current.

Publications carry component owners, so cleanup of an old panel cannot clear a
newer overlay. Source/event, primer or method changes invalidate the appropriate
publication. The main editor keeps a stable document identity for cursor/scroll;
a separate full event context protects overlays against equal-DNA aliases.
Map, sequence and read-only views suppress stale overlays. New workspace epochs
clear pending publications and snapshots.

Implementation: `web::pcr_overlay`, `web::pcr_overlay_draw`, selected-output panels,
shared renderer/sequence view and `web::pcr_source_view`. Seven independent native
fixtures cover aliased names, same-oligo roles, circular geometry, reverse mismatch
mapping, ambiguity/incomplete evidence, stale anchors and bounds. Browser visual
and interaction acceptance remains deferred until full in-scope code is complete.
