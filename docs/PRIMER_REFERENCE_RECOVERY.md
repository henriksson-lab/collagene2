# Existing-oligo reference comparison and binding-hit handoff

The shared pair editor distinguishes making a primer from selected DNA from
comparing an existing oligo to that DNA. **Use selection** still copies the chosen
directional domain. **Compare existing oligo to selection** instead previews a new
reference anchor while preserving the exact name, annealing text and tail. It
requires a valid bounded oligo and equal canonical domain/reference lengths;
invalid or unfinished text is retained rather than repaired automatically.

Both comparison directions are explicit, independently of the pair slot. The
preview displays actual oriented reference calls, literal domain calls, disjoint
mismatches and compatible ambiguity. Acceptance requires the same whole pair,
workspace epoch and selected source/region; it resolves the exact event and checks
actual selected DNA again before installing only the new anchor. Changing the
selection makes acceptance fail and requires another preview. No sequence alignment,
binding scan, thermodynamic calculation or amplification claim is inferred.

The single-oligo analysis view pages all returned binding hits in bounded groups.
The report's truncation and warnings remain visible; reaching the last returned
page does not establish exhaustive specificity. Zero or multiple hits offer local
guidance to edit the domain/tail or declared search options, inspect alternatives
and explicitly rerun. No rule is relaxed automatically.

Each returned binding hit can be explicitly copied into either editable pair slot.
This retains the analyzed oligo and actual directional source reference, including
intentional mismatches and circular-origin geometry. The handoff rechecks source,
epoch, calculation generation and form context, resolves the exact source, and
verifies the hit's coordinates, location and per-base facts without rescanning
the entire template. It changes only the chosen draft slot, leaves the other slot
intact and requires separate pair Apply and subsequent workflow analysis.

This handoff stores a **reference comparison**, not a persistent fresh binding
receipt. Native draft restoration retains the same weaker comparison contract;
incomplete search, alternatives or scientific authority cannot be upgraded by
copying a hit. The existing pair attachment version2 and worker schema remain
unchanged. Single-oligo annotations, primer orders and full evidence JSON remain
available with their source guards.

Binding navigation highlights every returned source part. A circular crossing is
represented by its two ordered ranges and does not enable a fabricated contiguous
base edit. Pair-slot identity is distinct from binding strand: a reverse-strand
placement can be explicitly put in either slot.

Implementation: `primer_pair_model`, `primer_reference_capture_view`,
`primer_binding_view`, the shared pair editor and document-host navigation. This
follows the [parallel workflow decisions](BENCHLING_WORKFLOW_DECISIONS.md) while
retaining explicit proposal/accept and Apply. Public Benchling documentation does
not establish the same preserve-literal capture semantics.

Browser acceptance remains deferred until all in-scope implementation is complete.

Five independent native fixture groups cover raw case/spacing and intentional-call
preservation, equal-length refusal without mutation, forward mismatch/ambiguity
facts, reverse circular-origin placement and malformed/stale evidence rejection.
All 85 native web tests and the WASM compile passed for this checkpoint. These
fixtures validate the capture contracts; they do not establish experimental PCR
success or GUI usability.
