# Antibody numbering implementation

`collagene-antibody-numbering` is a new independent MIT Rust/WASM module for exact
position labels and the published IMGT CDR3 allocation rule. This is a numbering
primitive, not automatic antibody or variable-domain detection. No worker, UI,
MCP integration, full-domain numbering, Kabat/Chothia assignment or automatic
CDR/framework derivation is implemented yet. Those remain part of F34.

## Labels and conditional CDR3 allocation

A `Position` retains an integer number and an explicit insertion kind: none,
uppercase letters or a positive decimal insertion index. Labels such as `111.1`
and `111.10` remain distinct. Decimal suffixes are indices, never floating-point
fractions. Syntax parsing is strict and cannot certify a numbering scheme or a
biological position. It rejects leading-zero aliases, lowercase insertion letters,
zero decimal indices and malformed labels. No generic ordering is supplied:
source order must follow the actual scheme, not lexical or floating-point sorting.

`imgt_cdr3_positions` follows the published [IMGT CDR3 correspondence chart](https://www.imgt.org/IMGTScientificChart/Numbering/CDR3-IMGTgaps.html).
For shorter loops it leaves the specified central gaps. Longer loops receive
alternating extra positions on the 112 and 111 sides; source order traverses the
111 insertions upward and the 112 insertions downward. The supported length range
is 2–2,048 residues as a resource-bounded rule allocator, not a claim that every
such loop is biologically realizable. This function does not renumber other loops
or translate labels between schemes.

`number_anchored_imgt_cdr3` takes an exact immutable protein source and explicitly
provided zero-based Cys/J-anchor coordinates with provenance. It checks literal C
and F/W-G-X-G anchor context, source bounds, loop stops and length, then allocates
labels to the intervening source residues. The flanking anchors are excluded from
the CDR3 residue list. Ambiguous conserved anchor identities refuse. This is
conditional on the user's anchor interpretation: those letters alone cannot prove
that the positions are actual IMGT 104/118 or that an antibody/domain exists.

Source, anchor assumptions, rule data and results remain bound together. The
streaming source hash matches the shared molecule revision; source annotation
changes also invalidate prior contexts. `ConditionalCdr3::validate_for` reconstructs
from independently supplied source and anchors. Current engine is
`antibody-numbering-1`, with no older internal reader. The conditional source
profile is capped at 2,000 residues and uses the shared model/text/location bounds.

## Full-domain implementation still required

Numbering, domain classification and CDR definitions are separate layers. An
automatic detector must derive candidate variable-domain alignments against
explicit independently checked references, preserve ambiguous/partial outcomes,
and assign actual labels before deriving regions. It must support the accepted
scheme/chain/species scope, independently verify assignments and expose the
numbering/region tracks on proteins and alignments. Raw residue offsets and motif
matches cannot replace these requirements.

The [IMGT chart](https://www.imgt.org/IMGTScientificChart/Numbering/IMGTIGVLsuperfamily.html)
also warns that some cross-numbering correspondences cannot be converted
automatically. CDR1/2 positioning has additional context-specific rules; the
CDR3 allocator must not be reused as an unverified general loop-numbering rule.
Future reference manifests must retain source bytes/hashes, explicit numbering
scheme evidence, authorship and dataset terms. Depositor residue numbers alone
will not be assumed to be Kabat, Chothia or IMGT numbering.

The [wwPDB archive policy](https://www.rcsb.org/pages/policies) provides a CC0
source for future structure-derived fixtures. That policy applies to archive
data, not every integrated external annotation. No database dump, third-party
code, antibody sequences or trained model is bundled in this initial module.
Published numbering rules and entirely authored anchor fixtures are used here.

Four native fixture groups and the optional-schema WASM check pass. Fixtures
include literal published loop allocations, distinct long insertion labels,
a 2–2,048 length/uniqueness/boundary sweep, source/provenance/report tampering,
resolved and unresolved anchors, interrupted loops and resource/coordinate refusal.
These are not evidence of automatic detection, full-domain numbering or integrated
GUI/platform acceptance.
