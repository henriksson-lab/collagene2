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

Implementation follow-up: the four connections above now have explicit controls
in the existing pair/single-oligo route: unchanged-oligo comparison previews,
selected-hit handoff, full returned-hit pages and local edit/policy/rerun actions.
[Contracts and limits](PRIMER_REFERENCE_RECOVERY.md). This closes these specific
connections, not every PCR/mutagenesis workflow or final GUI acceptance.

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

Implementation follow-up: the dock now has genuine source bins and explicit
selected construct rows, with per-construct raw settings and stable recipe IDs.
[Contracts and remaining batch/output work](COMBINATORIAL_SOURCE_DESIGN.md).
This closes the basic alternative-choice connection, not the full bulk workflow.

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

Implementation follow-up: these computed-pairwise connections are now present as
explicit worker-backed search, cached actual-source annotation/CDS context and
aligned FASTA/Clustal/HTML/PDF exports. Existing report/history/coordinate engines
remain the foundation. [Details and bounds](PAIRWISE_PRESENTATION.md).

A direct read-check action in `AssemblyProductView` is an optional continuity
improvement. The supported existing route is to open the actual final product and
use its document-toolbar action. This is not a missing validation engine.

## Execution decision

Continue primer recovery and ordinary assembly connections first, then computed
comparison presentation/exports and genuine combinatorial source choices. Reuse
these documented interactions while retaining our simpler open-document model,
independent-target verification and exact-event evidence. No external scientific
runtime, deferred primer/MSA engine, BAM/BAI expansion or browser check is added.


## Further setup review: preparation, topology and recovery

A further parallel review of the same official [combinatorial guide](https://help.benchling.com/hc/en-us/articles/39922419457165-How-to-use-the-combinatorial-cloning-and-concatenation-assembly-tools)
and [ordinary wizard guide](https://help.benchling.com/hc/en-us/articles/39656605989901-Create-assemblies-with-the-assembly-wizard)
clarifies several interactions rather than adding another wizard:

- Preparation belongs beside the source bin. Benchling documents Gibson PCR,
  existing-overlap and digest preparation, with digestion limited to one bin;
  Golden Gate supports existing sites or primer preparation per bin. Preserve our
  independent per-row preparations without copying that product restriction.
- Topology can change while the selected cloning method stays fixed in Benchling.
  Our topology edit must invalidate scientific results and update the existence
  of the closing junction. Per-construct topology overrides are our capability,
  not something established by these documents.
- Source contribution orientation and forward/reverse product display are distinct
  operations. Benchling's ordinary mismatch recovery explicitly suggests reversing
  orientation and checking the preview. Offer local explicit edits; never flip a
  contribution or rewrite primers silently.
- The combinatorial table distinguishes blocking errors from warnings. Keep our
  diagnostic failures inspectable and exact drafts retained, while allowing output
  only through the existing fresh native finalizer and explicit incomplete policy.

Public documentation does not resolve mixed preparation among alternatives within
one bin, per-row topology overrides or exhaustive scientific correctness. Benchling's
wizard limitations do not define our physical-duplex model. These findings guide
interaction choices, not verification authority or final application acceptance.


## Further matrix review: explicit choices, diagnostics and finalization

The official combinatorial guide documents either autopopulating all combinations
or manually choosing desired combinations. Generated names remain editable.
Construct rows show overlap lengths and warning/error messages; selected rows open
sequence/map and primer review, with navigation between constructs. Table filtering
and sorting are documented. Finalization requires readiness, creates constructs
and optional primer/fragment outputs, then locks the assembly. The saved assembly
object remains available alongside its outputs and separate history.

Our selected-row matrix, retained editable draft and deliberate product publication
follow that interaction pattern. Status filtering and faster navigation between
reviewed constructs are useful subsequent improvements. The public guide does not
establish atomic transaction behavior, failed-row subset publication, deduplication
of equal combinations or per-row output-option granularity. Keep our contracts
explicit rather than attributing those policies to Benchling. Enterprise folder,
schema and notebook steps remain outside this project's scope.


## Further mutation and recombination review

The same combinatorial guide lists Gateway under Homology when input fragments
already contain ATT regions. It does not establish dedicated BP/LR selectors,
site-family compatibility, reciprocal-component review or recombinase chemistry.
Preserve our explicit profile/source/site/component contracts and independent
intent. The documented larger/smaller-fragment preference must not silently choose
our biological target.

Focused official-source searches did not establish current dedicated inverse-PCR
setup/recovery. This is a documentation gap, not evidence that Benchling lacks it.
Its [ordinary PCR guide](https://help.benchling.com/hc/en-us/articles/9684267245325-Run-PCR-in-silico)
does not establish circular mutagenesis treatment/ligation semantics. Retain our
existing independent mutation intent and explicit preparation route.

The distinct [CRISPR ssODN guide](https://help.benchling.com/hc/en-us/articles/37748593861133-Use-CRISPR-tools)
documents editing a reference copy, confirming changes, adjusting homology arms
and inspecting silent-mutation suggestions. This is donor-reference design. If
expanded here, reuse explicit proposal/Apply and source-bound intent; neither a
suggested mutation nor a donor reference establishes physical repair success.


## Further output handoff and export review

The [DNA/RNA overview](https://help.benchling.com/hc/en-us/articles/39110680274189-DNA-and-RNA-sequence-overview)
documents selected-sequence export as Multi-GenBank, Multi-FASTA, CSV or individual
GenBank files in ZIP. The [primer guide](https://help.benchling.com/hc/en-us/articles/39072057302541-Use-primers)
uses checkbox selection and a bulk CSV action within the sequence panel. This
supports our selected-construct format chooser and attached-primer exports without
adding enterprise search/folder management. Keep unrelated workspace oligos out
of a construct's declared primer sheet.

The current combinatorial guide places producing-assembly details in an Assembly
History tab on outputs. An [official clarification](https://help.benchling.com/hc/en-us/articles/17806823063053-Why-can-I-no-longer-see-the-assembly-history-in-the-lineage-tab)
distinguishes that tab from sequence lineage. Keep our producing-run, immutable
source events and saved outputs directly navigable; they need not be collapsed
into one history concept.

None of these sources establishes that ordinary GenBank/FASTA export preserves
the complete assembly graph. Our inference is to retain native workspace/manifest
links and expose format losses, rather than describe exchange files as full history
backups. The current help also does not establish generated wet-laboratory HTML/PDF
protocols with reaction quantities and conditions. Our computed protocols remain
an explicit project requirement; do not claim that feature as documented Benchling
parity. Atomic publication, exact-event replay and independent native verification
are our own contracts, not guarantees inferred from its finalization UI.
