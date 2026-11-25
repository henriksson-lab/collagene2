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
[Benchling review](BENCHLING_PCR_OUTPUT_REVIEW.md). Inverse-PCR still needs exact
source-event preferences forwarded into candidate review and connection to the
existing construct/primer/prepared-output chooser. Mismatch-PCR and standalone
primer-analysis output actions should reuse these contracts where appropriate.
No expected target is synthesized from a PCR result merely to enter an assembly
finalizer.
