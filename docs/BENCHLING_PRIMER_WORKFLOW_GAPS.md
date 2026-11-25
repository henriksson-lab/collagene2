# Primer and mutation workflow: remaining interaction decisions

Reviewed 2026-11-25 from official public documentation and current source. This is a delta to [the earlier primer/mutagenesis research](BENCHLING_PRIMERS_MUTAGENESIS.md), not a second feature inventory. No authenticated Benchling session, local browser, compiler or tests were used. No third-party assets or code were retained.

## What the public documentation establishes

Benchling's current primer guide describes creating a primer from selected bases with an explicit forward/reverse choice; its pair editor can take the reverse domain from another selection. Tail length is editable and introduced bases have distinct presentation. Its wizard separately offers constraints for designing across junctions and forcing binding endpoints. Saved primers can be selected after a binding search, linked into pairs, and selectively exported. These documented controls do not establish our implementation's thermodynamic equivalence. [Use primers, updated October 2, 2025](https://help.benchling.com/hc/en-us/articles/39072057302541-Use-primers).

The current assembly guide documents per-fragment preferred primers, reports unsuitable choices through status messages, and can design replacements. It also separates construct and primer review. This supports contextual feedback; our recommendation remains explicit approval before replacing a manual or saved oligo. [Combinatorial cloning and concatenation, updated December 12, 2025](https://help.benchling.com/hc/en-us/articles/39922419457165-How-to-use-the-combinatorial-cloning-and-concatenation-assembly-tools).

The historical 2014 wizard article explicitly describes changing binding/homology lengths and typing substitutions or insertions at primer-covered junctions, with affected oligos updated and changes undoable. This is historical evidence for the interaction, not verification of the present assembly UI. [Assembly Wizard, March 16, 2014](https://www.benchling.com/blog/assembly-wizard). Its Golden Gate follow-up likewise describes primer-covered point edits, parental mismatch labels and undo. [Golden Gate Assembly Wizard, December 16, 2014](https://www.benchling.com/blog/golden-gate-assembly-wizard).

Current generic sequence-editing guidance describes caret insertion and selection deletion. It does not state that arbitrary reference edits automatically recompile assembly or mutagenesis recipes. Treat those as separate actions until explicitly connected. [DNA and RNA sequence overview, “Edit DNA and RNA sequences”](https://help.benchling.com/hc/en-us/articles/39110680274189-DNA-and-RNA-sequence-overview).

No dedicated current inverse-PCR procedure was established by this review. Public mutagenesis mentions and historical junction editing are insufficient to infer Benchling's physical preparation sequence, mutation placement rules, or how it handles competing PCR geometries.

## Current source: gaps that are still real

These observations are source-level, not GUI acceptance results.

| Interaction | Already available | Remaining decision or connection |
| --- | --- | --- |
| Selection → oligo | [Selection toolbar](../crates/web/src/selection_toolbar.rs) opens the primer panel; [primer analysis](../crates/web/src/primers.rs) can use selected reference bases. | Neither entry offers an explicit reverse-oriented seed. Define which strand is copied, display the reference interval, and preserve oligo 5′→3′ orientation. |
| Manual/saved pair → mutation | Ordinary PCR and assembly preparations already use the [shared picker](../crates/web/src/primer_picker.rs) and strict resolver. | [Mismatch-PCR controls](../crates/web/src/mutagenesis.rs) remain separate manual domain/tail fields. Connect the existing picker/editor rather than add another primer library. |
| Junction → actual source and oligos | [Assembly workspace](../crates/web/src/assembly_workspace.rs) has intended-junction windows and source navigation. [Fresh junction review](../crates/web/src/assembly_junctions_view.rs) shows actual DNA, homology, frames and CDS evidence. | Actual junction review has no callback to select that product interval or open its contributing oligo pair. Introduce a junction identity that distinguishes intended boundaries from freshly reconstructed homology/physical boundaries. |
| Edit → consistent draft | Source regions, preparation options and proposed overlap oligos are retained; proposal application is explicit. | There is no shared sequence-centered transaction that updates declared intent and all affected oligos together. Direct typing in a document remains a reference edit, not such a transaction. |
| Mutation intent → inverse PCR | Existing inverse-PCR form accepts source interval, replacement, domain lengths and tail split, and checks a separately supplied target. | Users must translate a viewer selection into numeric boundaries and arrange the intended target themselves. Add an explicit edit-intent handoff; do not use the simulator's answer as its own expected target. |
| Orders and final output | [Plan oligo export](../crates/core/src/plan_oligos.rs), reviewed-plan CSV/TSV, and [finalization](../crates/web/src/assembly_finalize.rs) now exist, including optional primer collections. | Do not repeat the older note's claim that plan order export or optional output selection is entirely missing. Standalone primer and inverse-PCR results could hand off the selected plan/pair to these same controls. |

## Recommended adaptations, not claims about Benchling

**Use one explicit pair draft.** Its fields should retain source document/event/revision, domain interval and strand, full ordered bases, tail/domain split and user name. Offer Forward/Reverse “Use selection”, manual entry and the existing saved picker. Changing a split changes the interpretation of the same ordered bases; changing ordered bases is a separate edit. Show both before applying. A typed oligo without a unique location can still be analyzed; do not invent coordinates merely to fill the editor.

**Keep creation and verification primers distinct.** Overlap primers introduce homology to create a construct; a junction-spanning verification primer is a different design task. The documented wizard has the latter constraint, but our current anchored pair designer and overlap planner do not expose an equivalent junction-spanning binding requirement. Do not label ordinary overlap suggestions as junction-verification primers. A future task should require an explicit junction and minimum binding on both sides, using bounded source-bound analysis; this is not a prerequisite for the editor connection.

**Treat an intentional mutation as a draft operation.** A future Replace/Insert/Delete action should show frozen source bases, requested replacement and intended DNA before compiling the existing preparation. Present proposed oligo changes as a diff. One acceptance/undo action must cover intent plus affected oligos, and invalidate prior reconstruction, junction and order evidence. Unsupported internal mismatch geometry stays explicit; an edited target or a favorable codon display cannot establish a successful preparation.

**Keep errors attached to the user's choice.** Preserve the selected saved oligo when it cannot bind, misses the requested source interval, or violates current assembly constraints. Show the failed condition and offer a separately selectable suggestion. Ambiguous or competing geometries remain visible. This deliberately avoids silently adopting the replacement behavior described in Benchling's guide.

## Practical next implementation task

Build the shared directional pair editor and connect it to the existing primer panel and mismatch-PCR form first. It should:

1. Copy a forward or reverse viewer selection into one explicit named domain, with coordinates and tail shown separately.
2. Load either slot through the existing strict saved-primer resolver, preserving names, literal bases and domain/tail decomposition.
3. Render an aligned source/domain comparison for supported exact geometry and list intentional mismatches in oligo orientation; show unresolved placement otherwise.
4. Apply the pair only on request, then rerun the existing worker analysis; edits immediately clear earlier results without erasing the draft.

This closes a concrete everyday cloning interaction without new primer-design engines, Primer3 ports, external programs, projects, ELN or inventory. Sequence-centered mutation/junction transactions should follow on the same draft contract, rather than grow another disconnected form.
