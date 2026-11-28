# Native annotation editing

Native sequence replacement, destination paste and annotation-only strand flips
use one bounded annotation policy. The projection kernel, coordinate qualifier
proofs and complete INSDC location codec live in `collagene_model::annotations`.
Core retains reference-document construction, trace/FASTQ handling and physical
source-evidence checks. Old projection and codec entry points are removed; callers
use the canonical owners directly.

## Sequence replacement and paste

`Molecule::replace` accepts an explicit zero-based half-open interval. Empty or
ASCII-whitespace replacement deletes that interval; other input follows the
sequence alphabet's normalization and validation. Non-ASCII whitespace is invalid
sequence input, not an instruction to delete bases. A normalized replacement that
is identical to the original interval preserves the complete source model,
including raw qualifiers, review flags and the citation recipe.

Native editing retains endpoint envelopes: insertion strictly inside a feature
expands its envelope, whereas insertion at its start or end is excluded. This is
different from projecting only retained source bases. Ordered exon geometry and
strand are retained; editing an intronic gap can preserve the same biological
feature sequence.

Supported local `transl_except` and `anticodon` assertions are remapped only when
the parent biological calls, order, strand, phase and original continuity remain
intact. Inserting between previously adjacent joined parts, including a circular
origin crossing, cannot silently invent a new exon gap. Genuine existing intronic
gaps remain supported. Untouched single-span annotations already have a one-piece
native mapping and do not need the compound continuity proof. Touched, newly
interrupted or uncertain parents require review; affected coordinate, translation and phase
assertions become bounded source-evidence notes. Keeping an exception's three
positions does not justify retaining it after a frameshift elsewhere in its
parent. Fully deleted features are omitted without parsing discarded claims.
Malformed retained active claims and unsupported active `rpt_unit_range`/`coded_by`
mappings refuse. Amino-acid/anticodon payloads remain declarations, not interpreted
biological truth. The [qualifier profile](COORDINATE_QUALIFIERS.md) defines the
supported grammar and phase/TERM checks.

Destination paste uses this native edit policy for existing annotations, then
projects copied annotations independently. It does not pair old/new features by
index or run a second phase-removal pass. Complete deletion of an earlier feature
cannot accidentally invalidate or associate a later feature. Source and
destination context checks, explicit paste preview and immutable failures remain.

## Annotation curation

An annotation-only strand flip changes biological meaning while leaving DNA
unchanged. Coordinate, translation and phase claims become reviewed source
evidence, including malformed or unsupported coordinate text. Their positions
are not complemented as though the same biological calls had been preserved.
Repeated declarations and their literal values remain evidence.

Arbitrary manual feature relocation still refuses coordinate-bearing qualifiers
on either the old or proposed annotation. Removing a claim and moving the feature
in one request cannot bypass this check. Review/remove the claim at the unchanged
location before explicitly curating the relocated annotation. Other location
changes use the same bounded source-evidence operation.

## Admission and evidence

Native editing has an explicit resource profile in the existing `ModelBudget`
implementation. It retains import-scale support rather than reducing native
editing to the smaller shared request profile.

| Bound | Native annotation edit | Ordinary shared model admission |
| --- | --- | --- |
| Features | 100,000 | 4,096 |
| Location nodes | 1,000,000 | 16,384 |
| Aggregate qualifiers | 1,000,000 | 16,384 |
| Qualifiers per feature | 16,384 | 16,384 |
| Location depth | 32 | 32 |
| Original GenBank records | Separate 32 MiB aggregate allowance | 4 MiB per record within the shared text allowance |
| Annotation/metadata text | 8 MiB, excluding original records | 8 MiB, including original records |

Existing per-string limits and the one-million-unit annotation proof budget also
apply; the table does not guarantee every combination at the maxima fits. Notes
are at most 1 MiB. Raw replacement input is at most twice the ten-million-symbol
sequence limit. Borrowed admission precedes recursive validation and copying;
output admission completes before returning an edited model. Failed operations
leave their sources unchanged. Source-bound GenBank citation replay remains
mandatory for real edits; identity edits do not invent an edit record.

Native MCP's existing `edit_reference_document` tool additionally exposes typed
`ReplaceSequence` and `FlipAnnotationStrands` operations. It still has the ordinary
session/request admission profile, exact immutable source references, atomic
publication and reference-only receipts. Results expose `review_required` and
`annotations_requiring_review`, including marks already present on the source. It creates no physical material, native
workspace ancestry or verified cloning result. The registry remains 209 tools;
thin WebMCP permissions do not expand.

Scientific engine `science-34` and the standalone mutation-intent, assembly-edit,
junction-edit, domestication and reference-transfer engine stamps advance because
annotation content can change. Current worker/candidate/plan layouts are unchanged;
stale contexts are rejected without fallback implementations.

Independent fixtures cover literal forward/reverse/spliced mappings, phase and
partial TERM, edit boundaries, normalized identity, malformed claims, feature
removal, source immutability, clipboard indexing, strand curation, native/GenBank
round trips and real MCP SDK requests. Positive cases include 4,097 rich features
and an original record larger than 5 MiB. A full 100,000-feature single-span
case with ordinary label/note fields remains usable under the shared work budget;
a 100,001-feature input is refused.
Compilation and execution evidence is recorded separately in the
[delivery audit](DELIVERY_AUDIT.md). No browser/application HTTP checks are part of
this implementation batch.
