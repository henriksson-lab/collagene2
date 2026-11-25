# Benchling design-review continuity: current source audit

Reviewed 2026-11-25. Official public help only; no local browser, screenshots, builds or production edits. Distinct from the parallel assembly-setup review. Existing Benchling reviews were read first; earlier gaps now closed are not repeated as missing implementations.

## Documented interactions and present comparison

Benchling starts manual primer creation directly from a selected sequence region, asks forward/reverse direction, and opens its primer editor. It also supports populating another pair domain from a selection. Existing-primer attachment is a separate search/choose/attach flow. Linked pairs can be unlinked while primers stay attached; primer visibility is separate again. Export uses selected primers. These are distinct interaction operations, not proof of current binding or amplification. [Use primers, updated2025-10-02](https://help.benchling.com/hc/en-us/articles/39072057302541-Use-primers).

Our shared pair editor already preserves domain/tail/name and exact selected source events; saved oligo reuse, reference comparison, actual hit handoff, full hit paging and source overlays exist. `selection_toolbar.rs` already offers **Make forward primer** and **Make reverse primer**; `app.rs` binds these to the existing pair draft and primer panel. An initial review missed the separate toolbar when inspecting `SequenceView`; a shortcut is therefore not missing and must not be duplicated. Primer creation, reuse and comparison remain separate operations; no Primer3 integration is added.

Benchling's separate PCR guide links a selected pair, then offers Create PCR Product, options, and a Copy/save step. [Run PCR in silico, updated2022-10-20](https://help.benchling.com/hc/en-us/articles/9684267245325-Run-PCR-in-silico). Our ordinary and mismatch PCR selected-output panels already supply explicit geometry, named output, optional declared primer collection, CSV/TSV and exact-source guards; `pcr_outputs.rs:219–257` and `PCR_OUTPUTS.md` close the old automatic-all-products gap. No missing PCR output engine was found.

Benchling displays selected-region coordinates, length, GC and Tm in the sequence status area. [DNA/RNA sequence overview, updated2025-09-30](https://help.benchling.com/hc/en-us/articles/39110680274189-DNA-and-RNA-sequence-overview). Our separate selection toolbar already computes unique selected-reference length and bounded GC statistics, explicitly reporting known calls and ambiguous/unresolved limits. The footer alone is not the full selection UI. Thermodynamic analysis remains in the existing condition-labelled primer route; a generic selection Tm must not be presented as a PCR annealing recommendation.

For outputs, the same Benchling overview distinguishes per-sequence format selection from bulk selected-file GenBank/FASTA/CSV/ZIP exports. Our existing document exports, SVG map and released-construct archive already cover these purposes (`MAP_EXPORT.md`, `assembly_archive_view.rs:260–281`). No new format engine is justified by this research.

## Highest-value product review connection

Benchling exposes visual Assembly History from the resulting sequence's own tab, separate from generic history/lineage. [Assembly History location](https://help.benchling.com/hc/en-us/articles/17806823063053-Why-can-I-no-longer-see-the-assembly-history-in-the-lineage-tab); [current assembly guide](https://help.benchling.com/hc/en-us/articles/39922419457165-How-to-use-the-combinatorial-cloning-and-concatenation-assembly-tools).

Our exact product→run→source/product version navigation, editable fork, historical status, current replay and retained graph downloads already exist. `RunDetailsPanel` currently presents paged input cards and ordered product-step cards with expandable textual operation parameters (`run_capture.rs:575–617`). The run panel now adds a compact source→operation/preparation→product overview over those exact bindings. It pages ten steps and six source references per step, opens the existing guarded exact-event viewers, and labels historical/reference-only ancestry. Dependencies come from the shared typed-operation source function; matching DNA, current names and adjacent step order invent no edges. This closes the identified presentation connection; GUI acceptance remains deferred.

## Explicit unknowns / non-goals

The reviewed help does not establish exact-event immutability, bounded ambiguity/exhaustiveness rules, target independence, internal chemical-state preservation or valid protocol-reagent quantities. It does not establish that GenBank/FASTA/CSV or map exports preserve an assembly graph. Product-linked Assembly History is not evidence that export recovers it. Public Notebook PDF/HTML export is also not a documented automatic verified cloning-protocol generator. Retain our native workspace/history and fresh protocol verification; no ELN, inventory, project hierarchy, external runtime or undocumented scientific parity is added.

The compact existing run provenance view is now implemented; retain the existing selection toolbar and primer analysis route. Keep current PCR/output/export implementations and their exact-event contracts. GUI smoothness/layout acceptance remains deferred by the user.
