# Benchling verification, comparison and assembly-history workflow research

Research date: 2026-11-25. Primary public Benchling documentation was read against
current repository source. This is workflow research, not a live tenant usability
study or a claim that our unexecuted UI behaves exactly as its source suggests.
No third-party code/assets were copied. No compiler or local browser checks were
run. Scope excludes ELN/projects, inventory, NGS/BAM additions, MSA-engine ports
and Primer3 ports.

## Documented Benchling interactions

**B1 — Alignment review.** From a construct, open Alignments → Create new
alignment, upload traces or choose existing sequences, choose pairwise and assign
the template. Saved alignments reopen from the construct. Review provides mismatch
arrows, search, trace/quality and translation visibility, manual trimming and
realignment after edits/input changes. Circular comparison can rotate/reindex
alignment presentation without changing stored sequences. Alignment export includes
FASTA, Clustal and PDF. The documented DNA engine is MAFFT; copying this interaction
pattern does not require adopting that engine. These are documented interactions,
not evidence of automatic discrepancy acceptance. Source: [How to create
alignments](https://help.benchling.com/hc/en-us/articles/39922093033869-How-to-create-alignments),
updated 2026-07-21, sections Create new DNA/RNA alignment, About plasmid alignments,
Use alignment tools and Export alignment.

**B2 — Consensus versus template verification.** The official teaching workflow
imports sequencing results into a new consensus, then copies an identified insert
into a separately cloned vector. Its template-alignment exercise starts from the
known reference, which is preselected when opening alignment creation, and attaches
results there; translations can be shown to inspect changes. This documents two
separate starting intentions, not an automatic “correct reference” command.
Source: [DNA Sequencing](https://help.benchling.com/hc/en-us/articles/9684245367821-DNA-Sequencing),
updated 2022-10-13, sections Create a DNA consensus sequence alignment and Align
sequencing results to a DNA template.

**B3 — Editing and simultaneous inspection.** The sequence overview documents split
workspace views and direct sequence insertion/deletion. It recommends AB1 for trace
data and describes trace formats as carrying chromatograms/quality when imported
as alignments. This establishes supported interaction primitives, not a particular
read-evidence approval transaction. Source: [DNA and RNA sequence
overview](https://help.benchling.com/hc/en-us/articles/39110680274189-DNA-and-RNA-sequence-overview),
updated 2025-09-30, sections Supported sequence file formats, Edit and personalize
your sequence view and Edit DNA and RNA sequences.

**B4 — Product review and provenance entry point.** The assembly workflow offers a
pre-assembly results preview; finishing creates construct outputs. A construct's
Assembly History tab revisits its design. The official correction notice explains
that this information moved out of the generic History/Lineage side panel into
its own entity tab for visualization/metadata. Sources: [Combinatorial assembly
instructions](https://help.benchling.com/hc/en-us/articles/39922419457165-How-to-use-the-combinatorial-cloning-and-concatenation-assembly-tools),
updated 2025-12-12, final assembly/preview and View assembly history sections;
[Assembly History location](https://help.benchling.com/hc/en-us/articles/17806823063053-Why-can-I-no-longer-see-the-assembly-history-in-the-lineage-tab),
updated 2023-07-24.

**B5 — Saved representation.** The official alignment API describes sequence rows
with gapped calls and trim boundaries plus a reference-sequence link. That supports
a retained relation between inputs, reference and review representation; it does
not document deterministic scientific replay. Source: [Get NucleotideAlignment by
ID](https://docs.benchling.com/docs/api/benchling-api/nucleotide-alignment/get).

### Official visual references and limits

The official help pages link these assets directly:

- [Plasmid-alignment GIF](https://help.benchling.com/hc/article_attachments/39922065839885).
- [Alignment-trimming GIF](https://help.benchling.com/hc/article_attachments/39922093030797).
- [Assembly-history GIF](https://help.benchling.com/hc/article_attachments/39922491662477).

These links were resolved from the primary help articles. The research web tool
could not render their `image/gif` responses; no frame-by-frame visual claim is
made. The documented controls above are supported by article text. Assets remain
external links, not repository files.

The reviewed primary documentation does **not** establish a one-click
accept/reject-mismatch reference patch workflow, immutable evidence semantics for
such patches, or an assembly-history “rerun exact source revisions” command.
Those remain unverified, not asserted absent. A Benchling-hosted PipeBio article
also surfaced, but its product-specific Sanger feature list is not attributed to
the core Benchling UI here. No inventory or Notebook workflow is proposed.

## What our current source already connects

The following are repository observations, not Benchling claims. Links identify
the implementation inspected, and names below are visible control labels.

| Area | Existing implementation and concrete interaction |
| --- | --- |
| Raw trace | [trace_view.rs](../crates/web/src/trace_view.rs) plots channels, base calls and available native quality, supports zoom/fit/pan and keyboard call selection. [app.rs](../crates/web/src/app.rs) switches a trace document between sequence and chromatogram. Current stale-call evidence is explicitly refused rather than substituted. |
| Known-reference reads | [trace_validation.rs](../crates/web/src/trace_validation.rs), `TraceValidationPanel`: active document is reference; select up to 32 reads; choose trim/quality/orientation/support policies; Validate selected reads; inspect eligible coverage, unresolved sites, mapping columns and insertion observations. Site buttons select a reference base. Open support-only consensus and JSON/FASTA downloads exist. |
| Reference-free read assembly | [contigs.rs](../crates/web/src/contigs.rs): ordered read picker, trimming/quality options, overlap evidence and layout canvas; selected-column observations include current/original call and sample coordinates. Open source call already invokes an event-bound navigation callback. Keep assembly evidence, saved-assembly replay and explicit contig release exist. |
| Pairwise constructs | [compare.rs](../crates/web/src/compare.rs): active reference plus query picker; global/local/query-end-to-end modes; explicit affine scores/budget; differences select reference spans; CIGAR and JSON export. It explicitly uses supplied origins/orientation, without automatic circular rotation/reverse-complement search. Preview/differences stop at 4,000 columns/200 entries. |
| Imported alignment review | [alignment_view.rs](../crates/web/src/alignment_view.rs): retained alignment collection, bounded scrollable rows, keyboard cell navigation, original/native downloads, exports and ungapped row extraction. It imports existing alignments and does not perform MSA. Pairwise comparison results are not directly sent here by the inspected ComparePanel. |
| Generic document history | [history_view.rs](../crates/web/src/history_view.rs): event graph, exact event lookup, restore/fork, local-parent navigation, candidate-target content-match notice. Cross-document parents have names/revisions/IDs but only local parents receive Show parent event. Content matches are correctly distinguished from recorded ancestry. |
| Cloning run review | [run_capture.rs](../crates/web/src/run_capture.rs): Record cloning run binds exact source/target events and retained products. Review shows per-step operations/parameters, product document buttons, recorded engine/status, fresh Replay cloning run and replay JSON. Source bindings and independent target are displayed as identifiers, not direct historical-source navigation controls. HTML/PDF protocol and reconstruction exports are connected in app.rs. |

## Proposed practical UI connections

These are design recommendations inferred from the comparison. They are **not**
claims that Benchling already implements each detail or that a new scientific
engine is required.

### P0 — Reference discrepancy → exact source chromatogram

Add an event-bound source-call navigation callback to `TraceValidationPanel`,
following the existing `contigs::ReadNavigation` pattern. It currently receives
no navigation callback despite reporting original/current read call and sample
indices. A discrepancy should open its supporting read at the exact call/peak,
including reverse mapping and a clear reference coordinate. Selection must keep
the immutable source event pinned; an advanced document head must not silently
replace the evidence.

A compact split view should retain reference/aligned calls alongside the raw trace
window and native quality. Preserve the existing distinction between native ABIF
scores and calibrated probabilities. Initially reuse the existing trace painter
and stored alignment coordinates; this is a connection between existing evidence
and views, not a new basecaller or NGS feature.

### P0 — A review queue instead of a capped report dump

Provide Previous/Next discrepancy, an overview strip and filters for supported
mismatch, insertion, deletion, ambiguous/low-quality call and uncovered interval.
Keep these states distinct. Pagination or jump-to-coordinate must reach every
item currently only available beyond the first 200 sites/100 table columns in
JSON. Empty between-base insertion spans need a boundary marker, not a disabled
selection with no alternate navigation.

A review decision can be “retain reference”, “propose edited derived reference” or
“unresolved”, with a reason and exact source evidence. This is a proposed user
workflow; the primary Benchling documentation did not establish an equivalent
single-click approval mechanism. Preserve original reads and reference, preview
all proposed changes, apply only on explicit user action using current revision
pins, and create history rather than overwriting raw trace calls. Recompute read
validation after changes; old reports must remain historical. Do not treat a
human choice as a newly measured base or automatically fill uncovered calls from
the reference.

### P0 — Construct-owned saved verification

Add a concise “Read checks” list on the reference/product: exact reference event,
read source events, options, coverage/unresolved summary and date/engine. Reopen
its discrepancy/trace view from that construct. This complements the existing
saved reference-free read assemblies and consensus-document report attachments;
it is not a request to rebuild their persistence. A downloaded report alone should
not be the main route back to verification. Reopening must visibly distinguish a
historical assessment from a fresh rerun.

### P0 — Product → producing run → original source event

Expose the recorded producing run from the product's main review/history surface.
Within run review, make retained input and target identifiers actionable: open
that exact historical event without moving the source document's current head.
Reuse the existing historical-source navigation/fork behavior from saved contigs.
Add “inspect junctions”, “compare with independent target” and “validate reads for
this product” routes with explicit identities. Do not infer a producing run from
a matching DNA hash, candidate target or generic parent label.

The existing event graph, exact source bindings and replay are useful foundations.
A compact product Sources → Steps → Junctions → Read checks presentation can expose
them without adding Benchling-style project/folder administration. Recorded status
must remain distinct from fresh replay and experimental confirmation.

### P1 — Comparison-local orientation/origin policy

Offer explicit “as stored”, “reverse complement query” and circular-origin-aware
comparison views. Display the chosen orientation/offset and both native coordinate
mappings. This should affect comparison presentation, not mutate a source via
Rotate origin merely to get an interpretable result. Preserve ambiguity/ties and
bounded-search limitations rather than label all rotations a unique match. Keep
exact rotation equivalence separate from approximate alignment with indels.
Existing Rust alignment/coordinate machinery should be reused; no MAFFT/MARS port
or new alignment dependency follows from this workflow recommendation.

### P1 — Retain pairwise review, then rerun explicitly

Allow a fresh pairwise result to be retained as a source/target/event/options-bound
review object and reopened from either construct. A two-row alignment viewer can
reuse the existing retained-alignment presentation while preserving that alignment
columns are not native reference coordinates. Add search and discrepancy navigation
rather than asking users to scroll a truncated text preview. Changing sources,
trim choices or scoring should provide explicit rerun controls and invalidate the
current badge, without rewriting the historical comparison.

### P1 — Manual trace trim as derived input

Current read validation exposes automatic native-score-window policies, and FASTQ
selection trimming is separately available. Raw trace view lacks the equivalent
manual retained-range interaction shown in the documented alignment workflow.
An explicit selected-call interval/trim preview could produce a derived read while
retaining the original attachment and call/sample mapping. Reuse the trace core
rather than altering original peak data. Changing that range should rerun the
reference check and keep quality eligibility separate from a user-drawn trim.

## Suggested first practical walkthrough

1. Open a saved assembly product; follow its recorded run to inspect an exact
   source event and return without changing source heads.
2. Select reference plus two Sanger reads, run existing validation, then step
   directly from every mismatch/unresolved site to supporting chromatograms.
3. Record an unresolved decision or preview one explicit derived-reference edit;
   preserve originals and rerun validation against the new revision.
4. Reopen the saved check from the product; distinguish historical/fresh results.
5. Compare the construct with an independently selected expected sequence under
   explicit orientation/origin policy, then inspect relevant assembly junctions.
6. Export existing protocol/reconstruction/replay artifacts alongside read-check
   evidence, keeping expected cloning design distinct from observed read support.

This is a proposed acceptance walkthrough, not a report of executed tests.
