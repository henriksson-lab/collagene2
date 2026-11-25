# Search, aligned source context and portable comparison review

These connections extend the existing event-bound pairwise viewer. They neither
realign sources nor replace the native comparison engine. Fresh and historical
reports retain the existing exact source-event, independent navigation and explicit
recomputation contracts in [PAIRWISE_REVIEW.md](PAIRWISE_REVIEW.md).

## Sequence search

An explicit Rust-worker search operates over reference, query or both displayed
rows, skipping alignment gaps. A hit retains its complete alignment-column window,
ungapped oriented coordinates and original source ranges. Intervening gap columns
remain inside that window. Overlapping hits remain separate. Search does not inspect
unaligned source flanks, add another circular wrap at the ends of the displayed row,
or silently search a different query orientation.

Literal search compares normalized IUPAC symbols using a linear prefix matcher.
Compatible search uses a multiword bitset matcher over intersecting IUPAC call sets;
a compatible match is not a certain sequence identity. Both support patterns of
1–4096 calls. Default work is 32,000,000 units, with an explicit hard maximum of
160,000,000, and at most 10,000 returned hits. Work includes visited columns and
matcher comparisons/bitset words. An interrupted search reports incomplete; zero
returned hits cannot establish absence then. The result list is paged and each hit
can navigate the aligned page or its exact original source bases.

## Annotation and CDS context

The viewer retains the actual pinned reference and query molecules, not names or
matching DNA aliases. Annotation directions and optional CDS translations are
prepared explicitly in a private Rust worker, then cached across 600-column page
changes. Every mapped position keeps its original source coordinate; a gap has no
source base. Feature indices remain stable and overlapping features are paged.

CDS translation follows the complete original biological location, including join,
complement, circular spans, genetic code and codon_start. Page boundaries and gaps
never restart the frame. The residue is placed at the middle source base of its
codon, with dots at its other two source calls. Reverse CDS context follows biological
direction even when the display runs the other way. Terminal stop codons remain
visible as context, separately from the native retained CDS protein convention.
Repeated contradictory mappings use an explicit marker rather than choosing a call.

Unsupported, fuzzy or unresolved CDS geometry and unsupported translation settings
retain the annotation and an explicit withheld reason. Preparation has an aggregate
1,500,000-source-position limit and 64,000,000-base work budget. Withheld or unmapped
features remain visible in paged limitation notes. None of these tracks edits an
annotation, repairs its frame or infers gene function or physical chemistry.

## Portable outputs

Aligned FASTA and Clustal export the actual reported rows, including gaps and the
declared orientation/origin. They reuse the alignment-exchange writers through a
validated content-only path, without inventing original imported-file provenance.
Stable reference/query role identifiers avoid duplicate-name ambiguity. Standard
files cannot retain the complete native coordinate maps, source identity, scoring,
uncertainty and historical/fresh status; the export review discloses those losses.
An empty local alignment is an explicit error for row formats; rich review and JSON
can still represent it without fabricating aligned sequence.

HTML and PDF share a full review model: exact document/event/revision identifiers,
source display names, comparison policies/options, coverage and uncertainty, every
reported alignment column, original-coordinate difference ranges/carets and CIGAR.
They distinguish historical evidence from a guarded host label of fresh computation.
Export itself never recomputes alignment or supplies experimental verification.

PDF reuses the embedded font, pagination and renderer of the existing exporter, with
comparison-specific headings/footers. Base symbols use explicit fixed cell positions
so a proportional font cannot misalign the rows. This is a comparison review, not a
cloning protocol. No browser print API, external converter or scientific executable
is used.

Generation is cancellable and runs in a private worker. Format changes, source/form
changes, workspace restoration and component cleanup invalidate late results.
Download requires the same retained comparison guard and an explicit action after
reviewing limitations. Exports are bounded to 16 MiB; PDF additionally uses existing
text/page/font limits. Oversized outputs fail explicitly rather than omit columns.
Full native comparison JSON remains available.

## Binary compatibility

Worker schema56 appends export125/result117, search126/result118 and context
preparation127/result119. Existing variants remain in place and the ordinary32 MiB
wire cap applies; no workspace-size exception is added. Borrowed request guards
precede copies/serialization and exact request revisions bind all generation inputs.
Standalone response guards check shapes/bounds, not scientific truth. The host
retains exact events and checks current request/generation/epoch before installation.
Candidate, plan, scientific-verifier and native-workspace schemas are unchanged.

This closes the specific comparison connections identified in the
[parallel Benchling workflow review](BENCHLING_WORKFLOW_DECISIONS.md). It introduces
no external aligner, multiple-alignment engine or additional BAM/BAI work. Browser
and visual/performance acceptance remain deferred until all in-scope code is complete.

Validation:40 focused native fixtures passed: search/context8, export5, protocol7,
existing alignment exchange12 and existing PDF8. The combined Trunk app/worker
build passed. Checks cover literal coordinate/codon cases, gap and reverse/origin
mapping, malformed/stale inputs, incomplete searches, version/budget gates,
standard-file round trips and independently decoded PDF text/cell positions. No
browser or HTTP smoke checks were performed. The read-only8082 preview was refreshed
with the immutable schema56 build and retained older worker bundles.
