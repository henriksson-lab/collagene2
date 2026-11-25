# Catalog-pinned digest preparation

`Operation::CatalogDigest {template,selection,fragment_index}` is appended as
operation tag 17, after original-source terminal declarations. The shared
`Preparation::CatalogDigest {selection,fragment_index,repair}` is appended as tag 4.
The scientific engine is science14; worker/candidate/plan version gates distinguish
these new variants without reordering older variants or changing report layouts.

`selection` is the existing `catalog_digest::CatalogSelection`. Reconstruction
loads the embedded catalog, validates its exact edition/content pin and each
ordered selected ID/name/recognition/cut geometry, then resolves those enzymes.
Tampered IDs, geometry, names or catalog pins are rejected; a name never silently
falls back to the small legacy starter list or a supplier alias. Any supported
geometry in the complete embedded catalog can be selected. Documentary supplier
conditions are not chemical or experimental authority.

The actual digest path is shared with existing `Digest`: resolve the fresh tracked
source, invoke `modified_duplex::digest`, select the explicit zero-based fragment
in existing top-reference cut ordering, project annotations with the existing
source map and preserve actual strand/end/modification state. Imported detached
physical summaries cannot restore a substrate. `Source::Step` consumes only a
fresh supported tracked intermediate, preserving restrictive reference-only and
unknown-state behavior.

This remains the same **assumed-complete cleavage model** as legacy Digest. It
does not become a positional methylation check merely because the enzyme comes
from the larger catalog. Explicit source-generation states that are known
modified or unknown retain the existing incomplete-cleavage warning. Checked
modification-aware digestion remains a separate operation/profile with its own
supported sensitivities and limits. No complete cleavage, reaction compatibility,
kinetics or experimental outcome follows from supplier observations.

The preparation compiler emits `CatalogDigest` and, only when requested, a
separate `RepairEnds` step. Both overlap and ligation use the same preparation
mapping, preserve original source documents/revisions and apply orientation to
actual prepared material. `repair: None` retains real cohesive geometry. Borrowed
selection shape/pin checks precede cloning, and the existing preparation work
reservation charges `2*N*sum(recognition_lengths)` from every resolved enzyme.
Selections retain their existing 64-enzyme/32-KiB bounds; aggregate native-model,
100-million preparation comparison, materialization and output bounds remain in
force. HTML and PDF protocols print the full exact selection and assumption scope.

## Digest-region intent and reverse orientation

For an original intact source, `FragmentSummary.top_start` is the source reference
boundary of its retained top strand. `FragmentSummary.bottom_start` is an
**unrolled global original-source top-reference boundary** of its retained bottom
span. It is distinct from `FragmentSummary.duplex.bottom_start`, which is local to
the selected fragment. The existing physical mapper checks:

`global_bottom_start - top_start == duplex.bottom_start`.

Forward intended reference contribution uses the top interval. Reversing a real
staggered duplex swaps its stored top/bottom strands: its new top is the reverse
complement of the original-source interval
`[global_bottom_start, global_bottom_start + bottom_length)`, modulo the original
length on an explicitly circular source. Reversing only the old top interval can
give the wrong reference origin/sequence. Repair can change geometry and synthesized
calls, so old cut intervals cannot silently serve as repaired intent geometry.

The literal nonstarter NheI fixture uses original circle `AAAAGCTAGCCCC` (13 nt),
recognition `GCTAGC`, top/bottom offsets 1/5. Its selected top begins at 5 and is
`CTAGCCCCAAAAG`; the global bottom span begins at 9, its local offset is 4, and
its stored bottom/new reversed top is `CTAGCTTTTGGGG`. Both ends carry actual
5′ `CTAG` cohesive sequence with cleavage phosphate/OH groups. The naive reverse
complement of the old top, `CTTTTGGGGCTAG`, is a different linear reference and is
explicitly rejected as the independent reversed target.

Focused independently written fixtures cover both strand strings/end geometry,
source-origin annotations, unchanged originals, two actual catalog-digested arms
ligating against a literal independent circle, optional explicit repair, tampered
catalog/source bindings, detached physical rejection, tracked methylation
uncertainty and reversed global-bottom coordinate intent. No external software or
new restriction/ligation algorithm is introduced.
