# Explicit mutation reference intent

`mutation_intent` builds an expected DNA reference from an actual immutable source and an explicit `Request {name, edit_start, edit_end, inserted}`. It calls no primer, PCR, ligation or assembly simulator. The intended target is therefore available before choosing or checking a preparation.

Coordinates are zero-based half-open source boundaries. Equal boundaries insert at a caret; empty inserted DNA deletes a region; both nonempty express replacement. Linear intervals must fit the source, and insertion at its terminal boundary is supported. Circular starts must be below the source length; unrolled ends can extend through at most one source turn. The exceptional circular terminal caret `start=end=length` inserts at the displayed end while preserving the original origin; deleting from a start at `length` remains invalid. A nonwrapping edit preserves the original origin. A wrapping edit first rotates the source to its selected start using the existing annotation transformation, then replaces the deleted prefix. Its output is `inserted + retained_source_path`, and `Report.origin_shift` records the original start boundary. A complete-turn replacement is valid intent if it introduces nonempty DNA, even when inverse-PCR preparation cannot implement it.

`source_preflight` validates bounded model content before copying or hashing. Sources are DNA and at most 8 MiB; inserted sequence has at most 100,000 normalized DNA calls and 200,000 raw ASCII bytes, including optional whitespace. Names are nonempty, control-free and at most 1024 bytes. Output and retained evidence must fit existing model limits and the 32 MiB report serialization profile. Ambiguous DNA is retained as reference intent, not promoted to a concrete preparation result.

Native `Molecule::replace` handles affected feature coordinates, review flags and stale translation removal. Unresolved coordinates and unsupported coordinate-bearing qualifiers fail explicitly. Original feature records are additionally retained under a source-revision-specific metadata key, including fully deleted annotations. Original raw GenBank and citation evidence follow existing source-bound native edit/rotation handling. Wrapping transformations may archive active original headers as source evidence rather than leave stale current coordinates. No source document or historical event is modified.

Readable physical or reference-only sources may supply reference letters. Their original `duplex_model` descriptor is preserved as ancestral evidence; only the temporary reference-transform copy removes that descriptor to use the coordinate engine. The output always receives `mutation:intent:reference_only=true` and a restrictive `duplex_model` explaining that physical material state is unavailable. Other retained chemical metadata remains source evidence. The intended target cannot become an intact duplex substrate merely because its letters are valid.

`intent_revision` binds the engine domain, actual source revision and literal request. `Report` retains that context, source revision, request, target, deleted DNA, origin shift and explicit limitations. `validate_binding` is bounded readonly consistency checking of those fields, target sequence receipts and metadata. It neither replays a molecular preparation nor makes imported report JSON authoritative. Hosts must supply current source/event bindings and compare the expected token before applying a response.

Six literal fixtures cover replacement, insertion, deletion, origin wrapping, full-turn replacement, affected/deleted annotations, immutable raw source evidence, readable physical sources, stale/tampered receipts and input bounds. This module does not choose mutation architecture, synthesize oligos, infer a repaired product or establish experimental success.

## Host integration and validation

The sequence toolbar and Mutagenesis use the same exact-range/caret eligibility.
Fuzzy, unresolved and compound locations are not flattened; selection is independent
of the primer length limit. Preview runs in a cancellable worker bound to source
identity, edit fields and workspace epoch. Opening recomputes the declared reference
from the exact retained source event, rejects altered target content, then stages
its child document and `mutation-intent-v1` evidence together. Repeat requests reuse
only the same target and source event; equal-content source aliases remain distinct.
Normal browser workspace recovery saves the resulting state. This is not a new
experimental product or a stored positive verification result.

Six core edit groups, two binary protocol groups, two real MCP SDK groups and two
native history-staging cases passed. Existing native web unit cases, the full
122-tool stdio smoke and final web/worker WASM compiler check also passed. Browser
checks remain deferred, so these checks do not establish GUI acceptance.

## Inverse-PCR review and outputs

The inverse-PCR candidate handoff now preserves the exact original source and
independent target events, rather than dropping source preferences. It synchronously
checks the retained result, form, generation and epoch before opening the general
plan panel. Fresh verification leads to the shared actual-product viewer, order
CSV/TSV, HTML/PDF and construct/primer/prepared-reference output chooser. Explicit
historical source events are retained without rewinding current heads. The active
target must still match its exact reviewed event. See [output finalization](CLONING_OUTPUTS.md).

An intended edit remains reference-only; this connection does not promote it to
prepared material or use simulated DNA as its own expected target. Existing
inverse-PCR PCR/treatment/ligation assumptions and finalizer assessment policy
remain explicit.

Mismatch-PCR now separately reviews every returned geometry, selected DNA and
substitution evidence in pages. Selected named reference/oligo outputs require
fresh native preview reconstruction and exact draft/event checks before atomic
installation. Search bounds and uncertainty remain explicit; this does not turn
the preview into a complete independently checked cloning plan.
[Reference output contract](PCR_OUTPUTS.md#conditional-mismatch-pcr-outputs).
