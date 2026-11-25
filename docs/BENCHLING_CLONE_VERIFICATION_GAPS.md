# Clone verification and history: Benchling workflow gaps

Research date: 2026-11-25. This follow-up checks the earlier
[BENCHLING_VERIFICATION_HISTORY.md](BENCHLING_VERIFICATION_HISTORY.md) against
public official Benchling documentation and current repository source. It is a
source review, not an executed browser walkthrough. No application code, external
assets, compiler checks or local browser checks were used. Scope is desktop
Sanger clone review and sequence lineage; ELN, inventory, scheduling, BAM and
new multiple-alignment engines are excluded.

## What the official workflow establishes

**Start on the construct, then attach reads.** Benchling's current alignment guide
(updated July 21, 2026) recommends uploading AB1 reads directly on the sequenced
plasmid. Users can instead select existing sequences, choose pairwise comparison
and explicitly assign templates. Saved alignments reopen from participating
sequences. The review offers mismatch navigation, trace/quality visibility,
trimming and explicit realignment after edits or changed inputs. Circular
alignment presentation can change origin without modifying stored DNA. This is
an interaction reference; its MAFFT/MARS implementation is outside our scope.
[How to create alignments](https://help.benchling.com/hc/en-us/articles/39922093033869-How-to-create-alignments).

**Unknown sequence reconstruction and known-reference checking are different
intentions.** The older teaching module (updated October 13, 2022) first creates a
new consensus from sequencing reads and copies its insert into a cloned vector.
A separate exercise begins on a known template, which is preselected for the
alignment, and attaches reads there for comparison. It warns that consensus
orientation may differ. This is teaching material, not evidence of the exact
current production dialog or an automated discrepancy-approval transaction.
[DNA Sequencing](https://help.benchling.com/hc/en-us/articles/9684245367821-DNA-Sequencing).

**A revised construct need not overwrite its starting version.** The current
copy/paste cloning guide (updated September 22, 2025) starts by choosing a
historical sequence version and cloning it, then renaming and editing the copy.
This supports an explicit fork-before-revision interaction, not a documented
claim that every pasted base carries chromatogram evidence.
[Clone by copy and paste](https://help.benchling.com/hc/en-us/articles/9684297597453-Clone-by-copy-and-paste).

**Design history has its own entry point.** The current assembly guide exposes
Assembly History on the resulting construct. The older July 24, 2023 correction
explains that assembly information moved out of the generic History/Lineage
sidebar into this tab. Neither establishes scientific replay, trace validation,
or experimental success merely because a construct has assembly history.
[Assembly workflow](https://help.benchling.com/hc/en-us/articles/39922419457165-How-to-use-the-combinatorial-cloning-and-concatenation-assembly-tools),
[history location correction](https://help.benchling.com/hc/en-us/articles/17806823063053-Why-can-I-no-longer-see-the-assembly-history-in-the-lineage-tab).

The reviewed pages do not establish exact event-bound evidence retention,
one-click mismatch acceptance, or immutable review decisions. These remain
**unverified Benchling behavior**, not asserted product omissions. The adaptation
below is our design recommendation.

## Current implementation and missing connections

| User task | Current source foundation | Concrete gap / suggested adaptation |
| --- | --- | --- |
| Bring in Sanger evidence | [trace_view.rs](../crates/web/src/trace_view.rs) renders original channels, calls and native quality; the app exposes trace documents and chromatogram viewing. | Add an **Add sequencing reads** route from the actual construct/product. Keep generic import available, but make the destination reference explicit before checking. Do not select it by matching filename. |
| Check a known construct | [trace_validation.rs](../crates/web/src/trace_validation.rs) selects up to 32 reads, runs bounded worker validation and displays eligible coverage, unresolved evidence and support-only consensus. | The panel receives document indices and molecular revisions, not stable document/event bindings. Upgrade its host contract to pin reference and read events; invalidate on equal-content event/document switches as well as DNA changes. Existing generation/revision checks protect changed content but do not establish event identity. |
| Inspect a discrepancy | Validation exposes read/reference mapping evidence; [contigs.rs](../crates/web/src/contigs.rs) already has event-bound `ReadNavigation` and exact current/original call and sample navigation. | Reuse that navigation contract in reference validation. One action should retain the reference coordinate and open the exact supporting chromatogram call, including reverse mapping. Today `TraceValidationPanel` has no source-call navigation callback. |
| Reach all differences | Reference validation shows bounded site/column previews and full downloadable JSON; [compare.rs](../crates/web/src/compare.rs) provides bounded pairwise differences. | Add paged Previous/Next discrepancy navigation, filters and insertion-boundary navigation. Do not require JSON inspection for differences beyond preview limits or collapse uncovered, low-quality and conflicting states into one mismatch count. |
| Reopen verification on a product | Reference-free contigs have saved, source-event-bound assembly evidence and fresh replay before consensus release. Reference validation can open a consensus carrying report/revision metadata. | Add a product-owned **Read checks** list with exact reference/read events, options and historical/fresh status. An embedded report on a consensus document is not the same as a saved verification relation on the checked construct. Reopening must not restore a current positive badge without reconstruction. |
| Review the actual cloning result | New [assembly_product_view.rs](../crates/web/src/assembly_product_view.rs) binds plan, independent target, current engine and exact final step before showing the actual product. [cloning_outputs.rs](../crates/core/src/cloning_outputs.rs) retains run history and controls prepared-document visibility. | Route **Validate this product with reads** from this actual final reference. Keep preparation products, intended reference and reconstructed final product visibly distinct. Finalization alone does not add observed read evidence. |
| Follow ancestry | [history_view.rs](../crates/web/src/history_view.rs) supports event selection, restore/fork and local parent navigation; [run_capture.rs](../crates/web/src/run_capture.rs) retains exact run sources/target and fresh replay. | Make cross-document source and independent-target event identifiers actionable, preserving current heads. Link product → producing run → original event; matching DNA content alone must not create that ancestry. |
| Make a revised sequence | Direct reference editing and document history exist; a support-only consensus can be opened separately. | Add an explicit revision proposal containing selected differences and reasons. Preview all edits, create/fork a derived reference, record exact evidence parents, then rerun checks. Never rewrite original chromatogram calls or silently use uncovered reference bases as observations. |

These are host/view/persistence connections around existing engines. They do not
require Primer3, a new basecaller, an external aligner or a generic laboratory
management layer.

## Proposed compact workflow

1. From a finalized construct, choose **Check sequencing reads**. Show the exact
   construct event and, separately, its independently specified design target.
2. Import or select AB1/SCF reads. Preserve raw bytes and call/sample mappings;
   choose trimming, orientation and support policies explicitly. Human trimming
   changes eligibility, not the original signal or measured confidence.
3. Review coverage and a navigable discrepancy queue alongside chromatograms.
   Distinguish supported matches, conflicting calls, insertions/deletions,
   ambiguous placement, unscored calls and uncovered intervals. Keep native caller
   scores distinct from calibrated error probabilities.
4. Retain reference, leave unresolved, or propose a derived revision for each
   reviewed difference. A human decision is a review record, not a newly measured
   base or chemical repair. A consensus assembled without a reference remains a
   separate route; it should not silently replace the expected construct.
5. Preview and explicitly commit revisions with exact source parents. Rerun against
   the new reference event; the previous assessment remains historical. Returning
   to an identical sequence through a different event does not reuse its authority.
6. Reopen checks from the construct and navigate back to original traces or the
   producing run. Export source-bound assessment and review decisions separately
   from the cloning protocol and intended-reference artifacts.

## Meaningful later acceptance cases

These are proposed checks, **not executed results**:

- Two reads agree over an insert while backbone bases remain uncovered: the view
  reports partial eligible support, never whole-plasmid confirmation.
- A reverse read crosses a circular origin: a discrepancy opens the correct
  original call/sample and preserves both coordinate systems.
- Duplicate read imports, tied placements and missing/native quality remain
  explicit; duplicated files are not assumed independent physical observations.
- Switching to another document/event with identical DNA invalidates the active
  review context. A late worker response cannot attach to the new selection.
- A single reviewed edit creates a new event/derived document, retains original
  reads and their bytes, and requires fresh validation. Undo, reopen and export
  preserve the distinction between historical review and current evidence.
- A failed final cloning step never offers its last successful PCR preparation as
  the sequenced final construct; a wrong-target but reconstructed final product
  remains inspectable with its mismatch status.

First implementation priority: event-bound validation capture and direct
mismatch-to-chromatogram navigation, followed by construct-owned saved checks.
Revision proposals should use those same identities rather than introduce a
second unbound sequence-edit path.
