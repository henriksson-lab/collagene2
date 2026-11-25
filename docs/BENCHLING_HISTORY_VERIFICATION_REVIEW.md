# Cloning history and sequencing review: current workflow audit

Research date: 2026-11-25. This supplements
[the initial review](BENCHLING_VERIFICATION_HISTORY.md) and
[the earlier gap list](BENCHLING_CLONE_VERIFICATION_GAPS.md). Several gaps in
those documents are now implemented; they must not be treated as an unchanged
backlog. This audit reads current source and public official documentation. It
does not report an executed local UI walkthrough. No application code, external
assets, compiler checks or local browser checks were used.

Scope is desktop cloning, construct lineage, pairwise comparison and Sanger
review. ELN, inventory, organizational functions, OpenGel and further BAM/BAI
work remain excluded or deferred. Benchling's external alignment engines are
interaction references only; they are not proposed runtime dependencies.

Subsequent implementation: saved runs now expose exact target/source/product
version actions. A shared read-only reference viewer opens the retained event;
explicit forking rechecks membership and preserves original document heads.
The remaining-source-navigation assessment below records the original audit and
is superseded for that connection by the [current run UI contract](RUN_LINEAGE_UI.md).
Pairwise comparison and integrated trace-review recommendations remain active.

## What the public documentation establishes

**Construct-owned review.** Benchling recommends attaching AB1 reads to the
sequenced plasmid. Its alignment dialog distinguishes template comparisons from
creating a new consensus. Saved alignments reopen from participating sequences.
Review provides mismatch arrows, search, row sorting, manual trim, traces,
quality and translations. Realign explicitly updates edited inputs, trims or
parameters. Circular alignment presentation can rotate/reindex sequences without
changing stored DNA. Export includes FASTA, Clustal and PDF. These documented
interactions favor a persistent review surface on the construct.
[How to create alignments](https://help.benchling.com/hc/en-us/articles/39922093033869-How-to-create-alignments),
updated July 21, 2026.

**Row order is a review preference.** Benchling documents name sorting as the
default and a position sort that saves the chosen order for an individual
alignment. Position means the first non-gap base. This is a useful usability
detail when reviewing reads along a plasmid, separate from their scientific
input order.
[Why are my reads not aligning in sequential order anymore?](https://help.benchling.com/hc/en-us/articles/24903603244557-Why-are-my-reads-not-aligning-in-sequential-order-anymore),
updated October 9, 2024.

**Assembly history belongs on the resulting construct.** The assembly guide
documents a saved assembly object that retains unfinished work and survives
finalization. Each resulting construct has an Assembly History entry point for
revisiting its design. The separate correction notice explicitly locates this
outside the generic History/Lineage sidebar. This supports a direct product-to-
design route, without adopting project administration.
[Assembly guide](https://help.benchling.com/hc/en-us/articles/39922419457165-How-to-use-the-combinatorial-cloning-and-concatenation-assembly-tools),
updated December 12, 2025;
[Assembly History location](https://help.benchling.com/hc/en-us/articles/17806823063053-Why-can-I-no-longer-see-the-assembly-history-in-the-lineage-tab),
updated July 24, 2023.

**Historical versions can seed new designs.** The copy/paste cloning tutorial
starts by selecting a sequence version and cloning it, then renaming and editing
the copy. This establishes a version-to-derived-design interaction, not an
immutable scientific evidence model.
[Clone by copy and paste](https://help.benchling.com/hc/en-us/articles/9684297597453-Clone-by-copy-and-paste),
updated September 22, 2025.

The reviewed sources do not establish deterministic replay against immutable
source events, one-click discrepancy approval, or experimentally verified sample
identity. Those behaviors remain unverified for Benchling. Our evidence and
replay requirements are independent design choices.

## Current implementation: close the stale gaps

These are source observations, not claims of browser-validated behavior.

| Previously requested connection | Current implementation |
| --- | --- |
| Event-bound reference/read checks | [trace_validation.rs](../crates/web/src/trace_validation.rs) now accepts exact reference/source information and a workspace epoch. Guarded evidence carries reference and ordered read event bindings. Equal-content documents are not substituted merely by DNA revision. |
| Reference-to-original-read navigation | The same panel has `on_read` and `on_reference` callbacks. Read-column and insertion-call controls navigate pinned reads; insertion boundaries use a reference caret. Current molecule call coordinates remain distinct from original attachment call/sample coordinates. |
| Reach all reported validation evidence | Uncovered ranges, insertions, reference sites and per-read alignment columns are paged. The old first-200/first-100-only gap does not describe this panel anymore. An option includes supported matches in the reference-site list. |
| Construct-owned saved checks | [read_checks_view.rs](../crates/web/src/read_checks_view.rs) provides capture, construct-associated listing, exact source events and explicit fresh replay. Capture recomputes evidence from retained originals and commits the workspace durably. Reopening historical checks does not restore a current positive assessment. |
| Explicit review decisions and derived references | [read_review_view.rs](../crates/web/src/read_review_view.rs), [read_review_capture.rs](../crates/web/src/read_review_capture.rs) and [read_review_history.rs](../crates/web/src/read_review_history.rs) now retain decisions, preview nominated edits and create derived references. Retain, unresolved, observed-base, deletion and insertion choices stay separate; the derived reference requires a new check. |
| Actual final-product inspection | [assembly_product_view.rs](../crates/web/src/assembly_product_view.rs) exposes the actual reconstructed final step, separately from the independent target. A failed final operation is not represented by its last successful preparation. |

Keep these implementations and their authority boundaries. Do not introduce
another persistence or approval system solely because the earlier research
listed them as missing.

## Remaining concrete connections

These recommendations adapt the documented interaction patterns. They do not
claim Benchling implements our exact evidence semantics.

### 1. Make retained run sources and target actionable

[RunDetailsPanel](../crates/web/src/run_capture.rs) still renders retained input
and independent-target document/event/revision identifiers as text. Its host
contract exposes only `on_open_product`, whose control opens a product document.
Add an exact-event navigation callback for retained inputs, target and recorded
product events. Opening a historical event should preserve the current document
head and offer an explicit fork if the user wants to edit it.

Keep the product's producing run easy to reach and group Sources, Steps,
Junctions and Read checks into a coherent review route. Content equality,
matching target DNA or generic parent labels must never manufacture a producing
run. The current recorded-status/fresh-replay distinction is already useful and
should remain visible.

### 2. Replace the separate trace/report round trip with connected review

Reference validation now reaches original calls, but its evidence view is tables
and text, while [trace_view.rs](../crates/web/src/trace_view.rs) renders the raw
chromatogram separately. A compact reference/aligned-read/trace view would retain
the selected discrepancy and native quality while moving through supporting
calls. Start with existing coordinate mappings and the trace painter.

Add Previous/Next discrepancy, coordinate jump and category filters. Existing
pagination solves reachability, but does not yet provide a focused mismatch
review queue. Keep conflicts, gaps, low-quality observations and uncovered
ranges distinct. Allow read rows to sort visually by placement or name without
changing the ordered source bindings used by checks and decisions.

### 3. Bring ordinary construct comparison into retained review

[compare.rs](../crates/web/src/compare.rs) still limits the displayed alignment
preview to 4,000 columns and differences to 200; full data is JSON-downloadable.
Its difference button is disabled for an empty reference span, so an insertion
does not have an equivalent caret-navigation action there. Unlike reference-read
validation, it has not gained complete paged evidence navigation.

Retain pairwise comparisons with exact reference/query events and scoring, reopen
them from participating constructs, and provide complete discrepancy navigation.
Reuse the retained alignment presentation where appropriate, while preserving
native-coordinate mappings and distinguishing imported alignments from fresh
computed results. Explicit rerun must follow changed inputs or options.

The panel currently compares stored origins/orientations and explicitly says it
does not search circular rotations or reverse complements. Add an explicit
comparison-local orientation/origin policy with displayed native-coordinate
mappings, ambiguity and search limits. Do not rotate the saved source merely to
get an interpretable comparison or silently infer topology from length.

### 4. Make manual trace trimming a retained review choice

The trace viewer selects calls/ranges, while reference validation supplies native-
score trimming policies. It does not expose the documented alignment-style
manual trim boundaries as a retained, per-read input to reference checking.
Offer a previewed retained interval, preserve the original attachment and
call/sample mapping, then explicitly rerun. A manually retained call is not a
new measurement of quality.

## Suggested implementation order and later acceptance

Finish ordinary cloning output paths first. For this review area, the next small
connected change is exact historical-source navigation from a saved run, followed
by complete ordinary pairwise navigation and saved comparisons. An integrated
sequencing/trace workspace can reuse the now-existing saved-check and review
machinery; it should not delay core cloning.

Later acceptance should verify product → run → exact source event → product
navigation after a source head advances; a reverse read crossing a circular
origin; differences beyond every preview/page boundary; insertion-caret
navigation; explicit rerun after manual trim or origin-policy change; and
reopening checks/decisions without promoting historical status. These are
proposed acceptance cases, not executed results.
