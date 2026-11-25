# Annotated reference clipboard

The editor's Clipboard menu offers **Copy with annotations**, **Copy reverse with
annotations** and **Paste with annotations**. This internal clipboard belongs to
the open workspace session and works without system clipboard permission,
including the HTTP preview. Existing plain-text clipboard actions and shortcuts
remain available in secure browser contexts; external plain text does not acquire
internal source provenance or annotations.

Copy uses an exact selected reference interval, ordered displayed regions, or an
exact annotation location. Joins/complements follow biological order; a global
reverse-complement choice reverses the complete selection. Circular wrapping
retains ordered parts. Fuzzy, ordered, between-base and unresolved locations are
not flattened into invented exact selections. Repeated/overlapping source parts
are explicitly rejected by this one-to-one projection profile.

The clipboard freezes the source molecule, exact source event, copy declaration
and native context. It survives document switches and source-head advances while
that original event remains retained. Workspace replacement clears it. System
clipboard permission is not needed, and copying creates no document or source edit.

## Preview and apply

Paste requires one explicit destination interval or insertion caret. An active
compound destination selection must be resolved to an explicit interval first.
The native worker reconstructs source letters and annotations from the actual
original molecule and copy declaration; serialized clipboard letters are never
reconstruction authority. The preview pages product DNA, active annotations,
warnings and unmapped evidence and offers complete preview JSON.

Apply reconstructs again in the worker. The host checks the source/destination
events, literal destination selection, epoch and displayed clipboard identity,
then stages one reversible edit. The destination event has both its prior head
and the original copied event as parents. A compact `reference-paste-v1` attachment
records the declared transformation and bindings. Workspace validation and byte
budgets precede publishing any document or attachment. No original source head
is modified. Undo/Redo restore exact prior/next destination events, including after
native save/load; unfinished unrelated drafts are retained.

## Annotation and chemistry semantics

Existing destination annotation geometry uses the native reference edit
transform. Touched annotations require review; stale translations and unsupported
coordinate-bearing qualifiers cannot survive through case-variant names. Native
reference edits now handle translation/coordinate-qualifier case variants
consistently too. Source feature projection uses the original declared parts
directly against the **final** destination topology. A break in a linear copy
preview is not automatically a break after a complete circular replacement;
insertion of extra destination DNA can still interrupt that same gene.

Complete source/destination feature evidence and relevant raw headers are retained
under revision-bound metadata keys, with bounded chunks for large evidence.
Unmapped features remain evidence, not invented active coordinates. Source
citation qualifiers become `source_citation`, so their numbers do not falsely
address destination GenBank references. Original qualifier values remain in source
evidence. Current trace/FASTQ/physical descriptors are archived as ancestral
evidence; paste does not invent remapped chromatogram calls or prepared chemistry.

Copied fragments and pasted documents have an explicit restrictive reference-only
descriptor. Cleavage geometry, phosphorylation, synthesis, repair and ligation
must still be reconstructed through actual cloning plans. Reference paste alone
does not certify a physical cloning reaction or experimental success.

## Implementation and validation

`core::reference_transfer` provides `Selection::from_location`, `RequestCopy`,
`Fragment`, `RequestPaste` and `Report`. Copy and paste contexts bind full source
models and declarations. `Fragment.context_revision` is the expected copied
identity, not merely the fragment sequence hash. Paste always rebuilds that
fragment from the actual retained source. Read-only validators check bounded
receipt structure; fresh native execution supplies transformation evidence.

Worker53 appends jobs120/121 and results112/113, with old-schema prefix rejection
and the ordinary 32 MiB limit. Candidate11, plan9, science-14, archive1 and MCP136
remain unchanged. `DocumentHistory::commit_with_source_events` adds explicit
computational parent events to the normal undoable controller; its ordinary
`commit` behavior is preserved. Hosts validate external parent ownership in the
complete staged workspace before publication.

Focused fixtures specify literal RC/ordered/wrapping DNA, source and destination
feature coordinates, circular gene continuity versus interruption, raw evidence
and citation handling, RNA/protein boundaries, stale contexts and binary gates.
History fixtures cover cross-document parents, content-identical ancestry,
rollback and persisted Undo/Redo. Browser acceptance remains deferred.
