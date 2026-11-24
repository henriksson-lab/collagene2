# Literal homology-donor templates

`core::edit_template` builds a donor and an intended edited-reference sequence
from an exact original DNA source, a literal replacement and explicit homology
parameters. A separately supplied intended target is a comparison requirement;
its bases and annotations are never copied into either generated molecule.
The engine is `edit-template-1`. This is sequence/homology design, not a model of
cellular HDR, nuclease activity, repair efficiency or an experimental outcome.

`EditTemplateOptions` gives an original zero-based edit boundary, number of
source bases to delete, replacement DNA, positive left/right arm lengths, donor
orientation and a maximum distance from a bound guide's nominal top-strand cut.
Replacement can be empty for deletion. Insertion uses zero deleted bases. Empty
edits are rejected; a replacement which leaves the biological reference unchanged
is explicitly not a consistent editing design. Inputs/replacement are concrete
DNA. Unresolved source bases are not silently chosen for synthesis.

Both homology arms lie outside the edit and use distinct source positions.
Linear inputs must contain both requested flanks. Circular arms may cross the
origin, but the left arm, deleted interval and right arm together traverse the
source at most once. If deletion crosses the original origin, the edited
reference starts with the replacement followed by retained source bases. Other
edits retain the original reference boundary. Circular targets can be rotationally
equivalent, with orientation preserved. Donor reversal changes the donor sequence
and mapped annotations, not the intended edit.

A `GuideBinding` includes the complete declared nuclease/filter options and an
exact eligible `GuideCandidate`. Local source/PAM/coordinate reconstruction
rejects forged or foreign guides before donor design. The nominal cut-distance
bound is an explicit design policy; it is not a measured repair limit. A missing
nominal cut or a cut outside the requested bound fails the associated design.
A template without a guide is allowed, but establishes no cleavage association.
[Guide and off-target model](CRISPR.md).

The designer concatenates source-derived arms and the literal replacement. The
independent checker inspects each supplied donor base, arm and replacement, and
walks original source coordinates to reconstruct the intended edited reference.
It does not call the designer or its arm builder. The checker separately compares
that reconstruction to the supplied target. Mutated arms, wrong replacement,
wrong length/orientation, a mismatched target and a violated guide-distance policy
cannot obtain `consistent=true`. A wrong target can produce an inspectable draft
report; it never changes the donor to satisfy that target.

Generated molecules retain projected source annotations and review-required
component features. Interrupted or unmapped source evidence remains on the
original source graph. Raw source attachments, chromatograms, physical strand
state and original GenBank bytes are not copied onto changed DNA. Donor chemistry
and strandedness remain unspecified reference-sequence properties. Source-event
history is the design parent; the target event is retained as requirement
provenance rather than a causal material parent.

Preparation bounds apply before validation/hashing: one million source bases,
up to 1.1 million intended-target bases, 8 MiB of rich source data per molecule,
4,096 annotations/metadata entries, bounded structural nodes and location depth.
Donors are at most 100,000 bases including both arms and replacement. Full report
serialization is bounded to 8 MiB. These are independent limits; rich annotations
or position vectors can reach a byte limit before a base-count limit.

`validate` checks bounded historical shape and internally consistent producer
claims. `validate_bindings_for` checks the current engine and exact source,
independent target, options, guide and candidate identities without redesigning.
Neither grants authority to imported flags. `check` computes a fresh independent
check; `replay` recomputes the complete design. UI opening/export requires fresh
worker reconstruction and an independent supplied-donor check. Reloaded producer
metadata starts unverified. Native workspaces and full JSON retain the recipe and
bindings; bounded GenBank/EMBL comments are summaries with explicit export loss.

Remaining F25 work includes calibrated activity/specificity models, validated
whole-build reference stores, bulged/variant search, additional nuclease biology,
editing mechanisms and experimental protocol models. This literal donor profile
does not replace those requirements.
