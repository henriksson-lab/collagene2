# Coordinate-bearing annotation qualifiers

Reference projection now remaps exact local `transl_except` and `anticodon`
positions through extraction, origin rotation, reverse complement and cloning
product construction. It uses the existing complete INSDC location codec and
one bounded parser for the enclosing qualifier fields.

Supported wrappers have `pos` and `aa` fields; `anticodon` may additionally carry
`seq`. Field order, repeated qualifiers and the original amino-acid/sequence
payload bytes are preserved. Only the position expression changes. Duplicate,
unknown or missing fields, malformed nesting and trailing non-whitespace text
reject an active projection. The [INSDC feature-table specification](https://www.insdc.org/submitting-standards/feature-table/)
defines the source syntax.

Positions must contain exact local ranges, complements or joins. Each retained
call must occur exactly once on the parent's biological strand and consecutively
in its biological order. Display-sorted spans cannot establish that order.
Translation exceptions must align with the parent's `codon_start` phase;
conflicting or invalid phase declarations reject the operation. Complete codons
contain three bases. A one- or two-base `TERM` exception is allowed only at the
end of the parent's biological sequence. Anticodons contain three bases.

For example, extracting bases 4–12 from a linear source moves an exception at
`(pos:7..9,aa:Sec)` to `(pos:4..6,aa:Sec)`. Reverse complementation changes its
coordinate orientation while retaining the same biological calls and the same
declared `aa`/`seq` text. This does not interpret recoding, apply translation
exceptions, validate the declared amino acid, or establish experimental biology.
The ordinary CDS translator retains its explicit unsupported-exception boundary.

Partial, interrupted or uncertain mapped features require review. Affected
coordinate, translation and phase assertions become bounded source-evidence
notes instead of active destination claims. Circular origin crossings can remain
continuous; equivalent cuts in a linear destination can interrupt the feature.
Wholly dropped features do not cause unrelated projections to fail. Identity
maps retain source coordinate text when biological continuity is unchanged.

Remote, fuzzy, ordered and unresolved embedded positions remain unsupported for
active remapping. `rpt_unit_range` and `coded_by` also retain their explicit
refusal profile. Wrappers are at most 64 KiB, source-evidence notes at most 1 MiB,
and parsing/mapping shares the existing one-million-unit projection work budget.
Input text is charged before cloning; failures leave sources unchanged.

Manual annotation relocation is a separate operation from projecting a molecule.
If the old or proposed feature contains coordinate-bearing fields, location
changes refuse until those declarations are reviewed or removed at the original
location. Other relocations deactivate `translation` and `codon_start`
case-insensitively and retain their values in source-evidence notes. New or
same-location annotations remain explicit user declarations. Native sequence
editing, destination paste and annotation-only strand flips now
use the [shared native edit policy](NATIVE_ANNOTATION_EDITING.md); arbitrary manual
relocation retains its separate source-correspondence guard.

Native MCP exposes these operations through `edit_reference_document`, which
publishes immutable session references with an explicit reference-only receipt.
It does not create physical material, workspace ancestry or cloning verification.
The projection checkpoint advanced the scientific engine to `science-33` and
standalone TOPO/recombination/reference-transfer stamps because annotations can
change. The connected native edit batch uses `science-34` and advances the affected
edit-engine stamps described in its contract, invalidating stale contexts.
Worker schema 123 and the current plan encodings retain their layouts.

Independent fixtures cover literal positions, spliced/reverse/circular cases,
partial termination, malformed and unsupported wrappers, work bounds, manual
editing, actual tailed PCR products and GenBank/native/MCP round trips. Current
execution evidence is recorded in the [delivery audit](DELIVERY_AUDIT.md).
