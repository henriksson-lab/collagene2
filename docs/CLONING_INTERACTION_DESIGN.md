# Cloning interaction design

The primary reference is Benchling's sequence-centred single-construct assembly
wizard. This is a workflow specification, not a claim that the present inspector
implements it. Its public written guides and animated walkthroughs were reviewed
on 2026-11-25. No Benchling code is used.

## User workflow

1. Start **Assemble** from an open sequence. Choose the method and circular or
   linear output. Use backbone plus insert slots for an ordinary construct; permit
   additional ordered fragments without introducing projects or organizational
   objects. Preserve a draft when changing document tabs.
2. Activate a slot, open its source in the sequence viewer, select its region or
   restriction boundaries, then **Set fragment**. A slot retains the actual source
   document/event, selected coordinates, orientation and preparation choice. Show
   source names, lengths and small directional previews. Allow a full source,
   inversion of a circular selection, reversal, replacement and reordering.
   Selecting a region states design intent; it does not prove that a molecular
   preparation produces that fragment.
3. Choose fragment preparation in context: existing material/overlap, restriction
   digest, or PCR. Show compatible enzymes and actual cut positions in the source
   viewer. For PCR, offer suggested or saved primers and keep the annealing domain
   visibly distinct from introduced sequence. Common homology and primer settings
   belong in one compact settings view; detailed filters remain optional.
4. Review the intended construct and freshly reconstructed product. Keep the
   sequence/map central, with fragment boundaries and junctions selectable.
   Display error/warning markers at the affected fragment or junction, not just a
   final report. Review oligo sequences and binding/homology lengths in a primer
   table, with an order-sheet export. Primer/junction edits should visibly show
   their effect on the construct and invalidate obsolete checks.
5. Create the chosen product, optionally prepared fragments and primers, with
   source relationships and the exact recipe retained. Provide HTML/PDF steps
   from the same plan. Saving into a simple document collection is sufficient;
   ELN registration, inventory and organizational schemas remain excluded.

For batch cloning, extend the same slots to bins of alternatives and an explicit
construct matrix. Reuse method parameters, fragment preparation and product review;
do not make single-construct users configure a combinatorial table first.

## Intent and independent verification

Normal manual design must work without importing a pre-existing complete target.
The selected source regions, order, orientation and explicit sequence edits can
specify the intended sequence independently of the preparation simulator. Freeze
that design intent before choosing primers or reconstructing the molecular steps.
Then verify the actual PCR/digest/assembly result against it. An optional separately
supplied target supports reverse design and AI checking.

Never construct the expected target by copying the verifier's product and then
call the equality check independent. A region selection is not a digest operation;
a plausible overlap is not evidence of unique amplification, reagent efficacy or
experimentally sealed DNA. These distinctions belong in specific diagnostics,
while the primary flow remains selecting, preparing and reviewing fragments.

## Current implementation and next gaps

The first dedicated assembly dock now retains named backbone/insert slots and
source-event snapshots while the main viewer switches documents. It can capture a
contiguous selection, a whole source, or the complementary circular path; reorder
and reverse slots; build an expected reference from those declarations; and show
its source contributions and bounded intended-junction windows. Source-derived
PCR suggestions, saved-primer reuse, preparation and fresh verification use the
existing methods panel against that frozen expected reference.

Draft declarations and manual preparation settings are saved in the native
workspace/autosave as an `assembly-draft-v1` attachment. Restoring the draft checks
its source-event snapshots but does not restore a positive check or frozen
analysis. Reviewing a proposed preparation creates an intended-reference document
with exact source parents and opens the existing plan review. Historical source
navigation opens a separate retained-source document without rewinding its head.

The virtual digest table now commits a freshly determined fragment into the active
slot, including the original source event, selected catalog definitions and both
strand intervals. The operation reconstructs that digest afresh. Reversal uses the
original bottom-strand interval for the intended reversed top DNA; top/bottom
lengths and positions need not agree. Explicit repair stays a separate operation.
Fresh overlap/ligation evidence also displays the actual final product in a
read-only sequence/WebGL map view and can launch fresh junction inspection.

This is the first connected workflow, not full Benchling interaction parity.
Direct cut choices now connect sequence/map restriction labels to the digest
boundary picker and a freshly reconstructed assembly fragment. The reviewed
product also exports declared primer CSV/TSV and HTML/PDF steps from its plan.
Remaining work includes integrated product/primer/junction editing and an
optional-primer/fragment output chooser. Keep these connections
central rather than adding more standalone parameter panels. Browser checks remain
deferred until implementation is complete, as instructed; no GUI acceptance claim
is made from source or compiler checks alone.

## Implementation order

1. Keep a dedicated assembly draft separate from whichever source document is
   currently visible. Switching source tabs must preserve the chosen output and
   populated slots, rather than treating each viewed source as a new target.
2. Connect an active slot to selection in the sequence/map viewer, including the
   retained part of a circular backbone. Bind exact source events and show each
   region and orientation before preparation. Build a frozen intent sequence from
   those declarations independently of reaction reconstruction.
3. Add a central intended/product preview and junction navigation. Put primer
   domains, actual cut/end problems, per-fragment warnings and preparation review
   alongside the selected region. Reuse the core planner, simulator and source
   bindings already implemented.
4. Connect output creation, optional primers/fragments and order/protocol export
   to the same reviewed draft. Preserve producing-run/source navigation in each
   product. Then extend slots to batch alternatives and method-specific controls.

Parallel reviews are complete and identify concrete connections to reuse:

| Review | Main implementation connections |
| --- | --- |
| [Restriction and Golden Gate](BENCHLING_RESTRICTION_GOLDEN_GATE.md) | Persistent slots; exact cut-pair and inverted backbone selection; existing-site versus PCR-added Type IIS preparation; connected primer/junction review. |
| [Primers and mutagenesis](BENCHLING_PRIMERS_MUTAGENESIS.md) | Directional selection into a shared oligo editor; saved-primer reuse; explicit edit intent; sequence-centred domain/junction editing; direct order sheets. |
| [Verification and history](BENCHLING_VERIFICATION_HISTORY.md) | Discrepancy-to-source chromatogram navigation; complete review queues; saved construct-owned read checks; product-to-producing-run and exact source-event navigation. |

The notes distinguish current documented behavior, historical examples, inspected
images and our proposals. Their local P0/P1 recommendations fit inside the main
cloning-first execution queue; they do not move read verification ahead of core
assembly work, add product scope or introduce new test gates. Reuse our existing
scientific engines, source guards, trace painter, history and export machinery.
Do not build another parallel set of unconnected parameter panels.

## Reference sources

- [Benchling assembly wizard](https://help.benchling.com/hc/en-us/articles/39656605989901-Create-assemblies-with-the-assembly-wizard):
  backbone/insert slots, selection in the source editor, reversing orientation,
  preview and output creation; includes public digest/ligation and Gibson animations.
- [Benchling combinatorial assembly](https://help.benchling.com/hc/en-us/articles/39922419457165-How-to-use-the-combinatorial-cloning-and-concatenation-assembly-tools):
  method parameters, ordered fragment bins, per-bin preparation, preferred primers,
  construct/primer review and optional outputs. Its registration and project
  machinery are outside this product's scope.
- [Benchling assembly wizard introduction](https://www.benchling.com/blog/assembly-wizard):
  source-region selection, editable binding/homology domains, junction edits and
  source history. This older article supplies historical interaction examples;
  current guides above are the primary reference.
