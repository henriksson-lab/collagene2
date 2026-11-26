# Primer proposals for an overlap assembly with internal edits

The intended construct is declared from original source intervals, orientations,
junction trims/inserts and literal internal edits. A proposal uses that declaration
to prepare every contribution: exact PCR for an unchanged core, or two source
PCRs, reciprocal extension and an outer PCR for an internally edited core.

Each outgoing junction uses its literal synthetic insert followed by a declared
length of homology from the next retained, edited contribution. This includes
edits within that next prefix. Reverse contributions place the reverse-complement
addition on the stored-source forward primer; forward contributions place it on
the reverse primer. The resulting products are oriented only when the recipe is
assembled. A circular construct also plans the last-to-first junction. A linear
construct has no outgoing tail on its last contribution.

The native algorithm ranks primer candidates using explicit conditions and
filters. Bounded geometry searches use original templates and declared endpoints;
they do not choose products by matching an intended target. The entire proposed
recipe is then reconstructed and compared with the independently built literal
reference. Neither a successful local primer score nor a matching diagnostic
intermediate establishes the whole assembly.

All original source revisions, optional originating runs, intent and design
settings belong to the request context. A shared work/material ledger covers the
roles together and reserves the final verification work. Capped searches and
originating limitations remain visible. Returned evidence is advisory; saved
product metadata cannot supply the next physical substrate.

The main workspace applies suggestions only after checking that the original
sources, intended construct and editable draft are unchanged. Only active PCR
oligos and the explicitly requested uniform overlap lengths may change. Existing
pins, origins, inactive fields and the independent target remain intact. Active
digestion, paired-oligo or declared-duplex routes must be changed explicitly;
conflicting nonempty tails cannot be silently overwritten. The ordinary complete
assembly verification remains the final workflow step after applying oligos.
Switching the active contribution is navigation; it preserves the proposal and
current selection, while actual declaration changes invalidate it.

Native validation passes29 focused cases, including unchanged single-fragment
behavior after the shared-budget refactor. Native compilation passes; the final
WASM check is running. Closed-origin replay and literal/annotation operations
retain their own bounds outside the100M comparison ledger. Worker85 and MCP
tool178 carry this advisory whole-recipe proposal; persisted operation semantics
remain science22. The public viewer has not yet been rebuilt for these controls.
Browser checks remain deferred until the complete in-scope implementation is
finished.
