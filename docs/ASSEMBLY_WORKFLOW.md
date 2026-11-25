# Preparing ordered overlap fragments

The **Assemble** action opens the persistent assembly dock. It retains source slots,
regions and orientation while switching documents, then declares an intended
reference before preparation and verification. The standalone methods inspector
remains available for separately supplied targets. The sections below describe
its shared preparation controls. PCR, Assemble and Plan are directly accessible beside Features, RE sites and Digest. The standalone workflow uses original open documents as source roles and the active document as an independently chosen target. The dock uses its pinned sources and frozen intended reference. Each ordered row explicitly chooses one preparation:

- Whole linear DNA: use the complete source. Circular DNA requires a declared preparation; the stored origin is not treated as a cut.
- Restriction digest fragment: state enzyme names and the 1-based fragment number. The compiler records the digest operation and zero-based core index. Optional explicit end repair can fill 5′ overhangs, remove 3′ overhangs or blunt ends; the default retains the actual cohesive state, and repair is a separate typed operation. Cleavage is an explicit model assumption, not experimental activity or reagent/methylation permission.
- Exact PCR: supply forward/reverse annealing domains separately from their 5′ extensions, all written 5′→3′. Exact binding and competing products are handled by the shared core verifier, not inferred from names or imported primer metadata.

The reverse-complement option acts after preparation. Moving rows changes assembly order. Enter one overlap length per junction, including the closing junction for circular products, and select Gibson, HiFi or In-Fusion advisory guidance. Manual preparation remains available; the source-bound basic primer suggestion route below uses the existing embedded heuristic engine. No external scientific process is called.

Compute method evidence sends original sources, explicit recipe and independent target to the worker. The core reconstructs prepared fragments and checks their oriented end homology. The UI displays prepared reference DNA and actual upstream/downstream junction words; full evidence download retains the compiled preparation graph, source bindings, findings and target assessment. Preparation, sequence consistency and published length guidance are separate findings and do not establish a unique mixture product or experimental success.

Review method candidate recompiles the same source-bound recipe into a typed plan, retaining original source inputs and explicit preparation steps. It does not promote cached prepared references to source material. Any source selection/content, row option/order, target or workspace epoch change invalidates displayed evidence; cancellation closes the private worker. Native plan verification and retained source-event history use the existing candidate workflow.

The desktop form supports up to12 rows. Focused core/MCP and WASM checks run in the background during implementation. Browser checks remain deferred until the full in-scope implementation is complete. Mobile work remains deferred.

A [small PCR + EcoRV example](../fixtures/prepared-overlap-demo.md) provides independently specified fragments, a 15-base junction and correct/wrong targets.


Exact-PCR rows also support explicit saved-primer copying from any open document
with usable primer annotations, including imported oligo documents. Forward and
reverse names, annealing domains and 5′ extensions are copied together and remain
editable; no field is selected automatically. The same controls are available
for PCR preparation in TA/GC workflows. Source event/revision changes invalidate
the chooser, and changed recipe fields invalidate previous assembly evidence.


## Suggest PCR overlap primers

For a PCR-prepared overlap assembly, select the source interval and orientation
for each ordered fragment, then choose the homology length and basic annealing
length/Tm settings. **Set from sequence selection** copies the current contiguous
source selection into a row; numeric source coordinates remain editable.
**Suggest PCR primers** proposes anchored annealing domains and source-derived
5′ homology extensions. It does not alter the row's current primers.

Review the ordered oligos and their separate annealing/tail domains, then use
**Apply suggested primers** to copy the suggestion into the preparation rows.
Run the existing method/plan check to simulate the actual PCR and assembly and
compare the product with the independent target. Changed source events, document
identity, intervals, orientation, settings or target invalidate a suggestion,
including switching between distinct documents with identical sequence content.
The suggestions do not establish binding uniqueness or experimental success.

The same preparation controls also serve the dedicated sequence-centred assembly
dock described in [the interaction design](CLONING_INTERACTION_DESIGN.md). Its
expected reference is declared from selected source regions before preparing
primers or reconstructing the molecular steps.

## Choose digest boundaries and export a reviewed design

With the assembly dock open, select a restriction label in the sequence/map viewer
or use the digest's site-row boundary buttons. Run the selected enzyme reaction
and choose first/second boundaries; linear sources also offer their physical
source ends. Circular selection follows the source-forward path, with an explicit
opposite-path toggle. Choosing the same circular cut twice retains a full turn.
Commit only a fresh determined fragment; additional internal cuts must be resolved
through an actual preparation, not a merged reference slice.

After method reconstruction, the actual-product review offers Primer CSV/TSV and
Steps HTML/PDF. These actions capture the same reviewed plan/target and use private
workers. Order export preserves added tails and every use's domain/mismatch
metadata. Missing order details stop the automatic all-required download; report
completeness is not amplification or experimental verification.
