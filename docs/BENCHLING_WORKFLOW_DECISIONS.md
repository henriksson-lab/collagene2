# Benchling workflow decisions after parallel follow-up

Reviewed 25 November 2026 by three parallel read-only research workers against
official public Benchling documentation and the current repository. This records
interaction decisions, not browser acceptance or evidence of Benchling's scientific
correctness. No application browser checks or compiler runs were needed for this
documentation checkpoint. Earlier gap tables are historical where superseded below.

## Primer creation, existing oligos and recovery

Benchling documents directional primer creation from a selection separately from
attaching an existing primer: configure a binding search, choose a returned hit,
then attach. Pair linkage and PCR product creation are subsequent explicit actions.
Its assembly guide exposes preferred-primer failures locally and identifies which
primer was actually used. Sources: [Use primers](https://help.benchling.com/hc/en-us/articles/39072057302541-Use-primers),
[Run PCR in silico](https://help.benchling.com/hc/en-us/articles/9684267245325-Run-PCR-in-silico),
[assembly guide](https://help.benchling.com/hc/en-us/articles/39922419457165-How-to-use-the-combinatorial-cloning-and-concatenation-assembly-tools).

Keep our existing editable pair → Apply → analysis → selected geometry → named
output route. Endpoint/tail previews, circular anchors, persistence, source overlays
and order sheets already exist; do not rebuild them. Complete these connections:

1. Offer **Compare existing oligo to selection** separately from **Use selection**.
   Preserve literal name, domain and tail; require equal canonical domain length,
   explicit direction and exact source-event reference. A comparison is not binding
   evidence. Current `SelectionSeed::oligo` intentionally replaces the domain and
   must remain the creation operation. This preserve-literal action is our inferred
   usability improvement, not a documented Benchling feature.
2. Let an explicitly selected computed binding hit populate a pair slot with its
   actual placement and unchanged oligo. Selecting a hit does not establish unique
   binding, especially when the search is incomplete. Reuse the existing pair and
   source contracts rather than inventing a second attachment model.
3. Page the single-oligo binding list beyond its current 100-row preview. Full JSON
   is already available, but selecting a later placement should not require it.
4. Put relevant edit/options/alternative-selection/rerun actions beside invalid
   domains, zero hits, multiple hits and unresolved comparisons. Preserve intentional
   edits; never silently relax mismatch rules or replace primers. Benchling's
   automatic replacement is not a requirement; our explicit proposal/Apply policy
   and deferred Primer3 work remain unchanged.

Public documentation does not establish exact endpoint-edit semantics, mismatch
preservation on selection replacement, exhaustive ambiguity handling or zero-hit
recovery. Do not attribute our algorithms or uncertainty contracts to Benchling.

## Ordinary and combinatorial assembly

The documented useful pattern is source slots/preparation → junction/product
review → named construct and optional primers/fragments. Ordered alternative bins
represent different constructs. Sources: [ordinary assembly wizard](https://help.benchling.com/hc/en-us/articles/39656605989901-Create-assemblies-with-the-assembly-wizard)
and [combinatorial assembly guide](https://help.benchling.com/hc/en-us/articles/39922419457165-How-to-use-the-combinatorial-cloning-and-concatenation-assembly-tools).

Our existing dock already has source capture/orientation, digest/PCR preparation,
selective junction proposals, actual-product review and shared output/protocol
creation. Improve the continuity and local recovery of that dock. Retain explicit
Apply; prepared reaction alternatives are not interchangeable design sources.
Genuine alternative source bins and a chosen construct matrix remain planned as a
layer over the shared preparation compiler, with distinct per-combination output
identity. No second wizard, ELN/project hierarchy or external design engine is needed.
See the [detailed assembly check](BENCHLING_ASSEMBLY_WORKFLOW_CHECK.md).

## Comparison and sequencing verification

Benchling starts known-reference checking from the construct, distinguishes it from
consensus reconstruction, and provides mismatch navigation, traces/quality,
annotation/translation context and explicit realignment. Saved alignments reopen
from participating sequences; aligned FASTA/Clustal and PDF are documented outputs.
Source: [How to create alignments](https://help.benchling.com/hc/en-us/articles/39922093033869-How-to-create-alignments).

Our exact-event comparisons, complete paging/difference navigation, construct-level
**Check sequencing reads**, original-call navigation, saved read checks, historical
replay and explicit derived-reference proposals already exist. The older
[verification gap table](BENCHLING_CLONE_VERIFICATION_GAPS.md) predates those
implementations. Preserve [read checks](READ_CHECKS.md) and
[review decisions](READ_REVIEW.md) rather than creating duplicate implementations.

Remaining computed-pairwise connections are sequence search, mapped annotation and
translation context within aligned rows, and aligned FASTA/Clustal plus HTML/PDF
output. Reuse native coordinate mappings and existing exchange/report primitives;
keep gaps, original coordinates, orientation/origin policy, exact events and
historical/fresh status explicit. The imported-alignment viewer's exports do not
make these actions available from a computed pairwise report.

A direct read-check action in `AssemblyProductView` is an optional continuity
improvement. The supported existing route is to open the actual final product and
use its document-toolbar action. This is not a missing validation engine.

## Execution decision

Continue primer recovery and ordinary assembly connections first, then computed
comparison presentation/exports and genuine combinatorial source choices. Reuse
these documented interactions while retaining our simpler open-document model,
independent-target verification and exact-event evidence. No external scientific
runtime, deferred primer/MSA engine, BAM/BAI expansion or browser check is added.
