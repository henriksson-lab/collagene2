# Selected PCR reference outputs

Ordinary exact PCR now computes into an explicit review instead of immediately
opening every returned geometry. The source is the original current document
event; names, annealing domains and tails are frozen with the result. If there
are several geometries, none is selected automatically. One geometry can be
chosen, named and retained as a reference document, with an optional declared
oligo collection. The named pair can also be exported as spreadsheet-safe CSV
or TSV through the existing oligo adapter, including ordered 5′→3′ bases and
separate domains/tails. No quantities, modifications or stocks are inferred.

The review pages all products and the selected product's DNA. Actual forward
and reverse binding roles are shown: the native model searches either oligo in
either role, including single-oligo amplification. The order table retains both
declared oligos, rather than quietly discarding one according to a chosen role.
Projected annotation findings remain visible. This is reference construction,
not an independent target check, physical end reconstruction or a prediction of
experimental yield. Ambiguous annealing remains unsupported by the existing
exact model; supported IUPAC tail bases remain literal reference uncertainty.

## Reconstruction and installation

`core::pcr_outputs::Request` binds the complete source, named oligos, selected
index, original preview revision and output options. `revision` performs bounded
preflight and computes a context identity. `prepare` reruns the existing native
PCR/annotation reconstruction and requires the indexed preview's full revision
to match before naming it or building optional outputs. An index alone, imported
preview, matching DNA or user-facing oligo name cannot authorize selection.

Worker50 appends job117 `PreparePcrOutputs` and result109 `PcrOutputs`. The
ordinary 32 MiB message limit applies. Both borrowed-prefix gates reject these
variants under older schemas before decoding the new payload. Candidate11,
plan9, science-14, native archive1 and the 136-tool MCP registry are unchanged.

The web host checks the exact original document/event/revision, literal six-field
primer draft, workspace epoch, generation and displayed result identity before
and after the asynchronous worker. Changing the source or draft clears review
and cancels pending output work. Disposing the app cancels the initial PCR
worker too; both phases check their owning view before accessing signals. The complete staged workspace is validated
before publishing any document or attachment. The new reference's history parent
is the exact template event; that event is also retained in its PCR metadata.
Primer-only output leaves the current source active. Usual workspace autosave
and native export retain the installed references and collections.

An opaque `pcr-output-reference-v1` JSON attachment records the original source
event, request context, geometry coordinates/roles, original preview revision,
named oligos and chosen output event/collection attachment. It records reference
provenance, not verification authority. It deliberately omits duplicate full
product DNA. `oligo-table-v1` uses the existing collection format; oligos outside
its reusable-primer profile retain full literal order evidence and diagnostics
without becoming falsely reusable primer records. Oversized table fields fail
explicitly without imposing the table limit on product-only reconstruction.

## Validation and remaining work

Focused native fixtures cover independently specified alternative geometry,
literal product DNA and feature coordinates, source/name/index mismatches,
optional outputs, complete ordered tails, spreadsheet safety, IUPAC/long tails
and new binary gates. Browser acceptance remains deferred until the full in-scope
implementation is complete.

This closes the ordinary-PCR automatic-all-products handoff identified in the
[Benchling review](BENCHLING_PCR_OUTPUT_REVIEW.md). Inverse-PCR now forwards exact
source/target events into general candidate review and reuses the actual-product,
construct/primer/prepared-output chooser and order exports described in
[cloning outputs](CLONING_OUTPUTS.md). Mismatch-PCR and standalone
primer-analysis now reuse the same oligo table and collection contracts, as
described below.
No expected target is synthesized from a PCR result merely to enter an assembly
finalizer.

## Conditional mismatch-PCR outputs

Mismatch-PCR now retains the actual submitted named pair, source molecule/event,
viewed-reference event, raw form and applied draft with its private worker result.
All reported products are reachable in pages; selected DNA, primer substitutions
and projection findings are also paged. The selected geometry shows actual oligo
roles, including same-oligo amplification, source/product mismatch coordinates,
certainty and search completeness. An alternative is never selected by index-zero
convention when several results exist.

`core::mismatch_pcr_outputs::prepare` freshly calls the existing mismatch preview
engine and requires the selected original preview revision before naming outputs.
The result preserves the selected native bindings/substitutions and a full search
summary: options, products reported, pair checks, unsupported/conflicting/ambiguous
counts, rejected geometries, completeness and warnings. It does not duplicate the
DNA of every unselected product in the output response. Capped searches stay
incomplete; a selected preview is not evidence of uniqueness or target agreement.

Worker51 appends job118 `PrepareMismatchPcrOutputs` and result110
`MismatchPcrOutputs`, with old-schema prefix rejection and the ordinary 32 MiB
envelope limit. Original source/viewed-reference events, epoch, literal form,
applied pair and displayed result identity are checked around the worker and at
installation. The host stages the optional named reference, optional oligo
collection and compact `mismatch-pcr-output-v1` provenance attachment together,
then validates the complete workspace before publishing. The reference history
parent is the original source event. The viewed reference is context, not a
verified independent target of this output operation.

An existing exact-PCR or mismatch-PCR preview never gains physical verification
through output creation. The existing inverse-PCR graph/finalizer is the separate
route for an explicit complete preparation checked against independently declared
DNA.

## Standalone oligo outputs

`core::oligo_orders::from_primers` constructs the existing bounded oligo table for
one or more literal named primers. It preserves row order, duplicate names, full
ordered bases and domain/tail fields. It reuses the existing table's field/byte
limits and evidence-only policy for oligos outside the reusable-primer profile;
no binding, synthesis or stock state is inferred. `pcr_outputs::primer_table`
delegates to that helper, so exact and mismatch PCR use the same format.

The primer panel offers CSV/TSV and collection saving for a single analyzed oligo
or a particular ranked pair. Pair actions are collapsed behind “Order or save this
pair.” Source event, epoch, generation and literal calculation/form fields are
checked synchronously before actions; no action relies only on an eventual render
invalidation. Collection saving uses the existing atomic oligo installer; export
uses the standard spreadsheet-safe format while retained source rows keep literal
names. Changing a name/domain/tail requires analyzing or selecting again, rather
than exporting stale calculation evidence under new inputs.
