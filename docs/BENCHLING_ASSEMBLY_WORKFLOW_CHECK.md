# Benchling assembly workflow check

Reviewed 25 November 2026 against official public documentation and current local source. This is a static workflow audit, without local browser checks, builds or runtime acceptance. It supplements [the earlier assembly review](BENCHLING_ASSEMBLY_WORKFLOW_REVIEW.md) and [PROJECT_PLAN.md](PROJECT_PLAN.md); the implementation checkpoints in the latter supersede older missing-feature tables, particularly [the Golden Gate gap note](BENCHLING_GOLDEN_GATE_WORKFLOW_GAPS.md).

## Documented behavior and limits of the evidence

**Ordinary single-construct assembly.** Benchling starts from an open sequence and a method choice. Restriction/ligation uses retained-backbone and insert slots filled from cut-bounded source selections, with selection inversion for the backbone and an orientation correction. Preview precedes naming and output creation. The same guide describes Golden Gate as this flow with Type IIS enzymes. Gibson selects a linearized backbone and insert, reviews the result and creates sequence/primer outputs. Folder choice happens at creation. Its FAQ describes sense-strand-cut handling for Gibson and directs already digested material to a separate copy/paste route; it does not establish a general physical-duplex model. [Official assembly-wizard guide, updated 30 September 2025](https://help.benchling.com/hc/en-us/articles/39656605989901-Create-assemblies-with-the-assembly-wizard).

**More complex preparation and combinations.** The newer tool saves an assembly object and configures ordered bins. Gibson production distinguishes PCR-added homology, existing homology and restriction digestion; Golden Gate distinguishes existing sites and primer-based preparation. Multiple alternatives in a bin generate alternative constructs rather than multiple contributions to one construct. Preferred existing primers are checked and may be replaced with a warning. Fragment/construct status, per-junction overlap lengths, sequence/map review and a primer table support correction before creation. Finalization offers constructs and optional primer/fragment outputs; assembly history remains accessible. Its separate Homology mode includes an ambiguity preference. Exact primer-replacement and orientation behavior should not be generalized across its different dialogs. [Official combinatorial cloning guide](https://help.benchling.com/hc/en-us/articles/39922419457165-How-to-use-the-combinatorial-cloning-and-concatenation-assembly-tools).

**Existing digested material.** The documented copy/paste route clones a chosen plasmid version, selects the insert between restriction sites, replaces a selected destination region and confirms insertion. This establishes a reference-editing interaction, not its chemistry or completeness guarantees. [Official copy/paste guide, updated 22 September 2025](https://help.benchling.com/hc/en-us/articles/9684297597453-Clone-by-copy-and-paste).

**Historical design context.** The 2014 Golden Gate announcement documents source selection, site-bearing primer generation, spacer handling, primer export and a printable summary. It supports the usefulness of keeping preparation and outputs together; its screenshots are historical, and its four-base wording must not become a universal Type IIS assumption. [Official Golden Gate announcement, 16 December 2014](https://www.benchling.com/blog/golden-gate-assembly-wizard). The earlier wizard announcement documents editable Gibson overlap/Tm parameters and restriction incompatibility/internal-cut warnings. [Official assembly announcement, 16 March 2014](https://www.benchling.com/blog/assembly-wizard).

These documents describe product behavior, not implementation correctness. They do not establish exhaustive mixture search, physical end-state tracking, source-event replay or experimental success. No undocumented Benchling algorithm is assumed here. Linked animations were not inspected frame by frame.

## What already exists here

The following are source observations. They must not be reclassified as missing engines because an older review predates their implementation.

| Workflow stage | Current source foundation | Consequence for further work |
| --- | --- | --- |
| Select backbone and insert | [assembly_workspace.rs](../crates/web/src/assembly_workspace.rs) defaults to named Backbone/Insert slots, retains captured original sources, regions, orientation and ordered contributions, and serializes a bounded draft. | Extend this dock; avoid a disconnected second wizard. |
| Restriction preparation | [digest.rs](../crates/web/src/digest.rs) and the dock connect actual reconstructed cut-bounded fragments, complementary circular paths and source highlighting. | Keep recognition ranges, strand cuts and physical fragments distinct. |
| Gibson/overlap preparation | [cloning_methods.rs](../crates/web/src/cloning_methods.rs) provides explicit PCR/digest/existing-material rows, saved oligo reuse, suggestions followed by Apply and method profiles. | Primer suggestion and ordinary overlap preparation are present. Improve the continuous interaction instead of adding another primer tool. |
| Golden Gate preparation and alternatives | The same methods/dock code supports PCR-added Type IIS preparation, existing preparations, selected prepared pools and explicit whole-original/amplicon reaction substrates. [golden_gate_reaction.rs](../crates/core/src/golden_gate_reaction.rs) and [prepared_golden_gate.rs](../crates/core/src/prepared_golden_gate.rs) retain bounded alternatives and physical preparation. | Older claims that the dock lacks Golden Gate or PCR→digest composition are stale. Finite complete-cleavage/one-copy enumeration remains a distinct model from kinetic cycling. |
| Junction editing and evidence | [selective_junction_edit.rs](../crates/core/src/selective_junction_edit.rs), [golden_gate_junction_edit.rs](../crates/core/src/golden_gate_junction_edit.rs) and [assembly_junctions.rs](../crates/core/src/assembly_junctions.rs) support selective regeneration and actual selected-path junction windows. The methods UI provides row-level recovery/source actions. | Local diagnostics and Golden Gate junction inspection already exist. Review their presentation and continuity, without claiming they are absent. |
| Product and output review | [assembly_product_view.rs](../crates/web/src/assembly_product_view.rs), [assembly_finalize.rs](../crates/web/src/assembly_finalize.rs) and [cloning_outputs.rs](../crates/core/src/cloning_outputs.rs) connect actual reconstructed product review, named construct creation, optional oligos/preparation visibility, retained source events and order/protocol exports. | Reuse these outputs for every method. HTML/PDF and primer export are present; no new export engine is implied by this comparison. |

## Implementation follow-up

The editable reaction-choice persistence gap below is now closed in source:
`MethodsState` retains an optional bounded pool form in the existing draft.
Substrate choices follow reordered/removed roles; older drafts use their topology
and existing defaults. Receipts remain transient and require explicit inspection.
[Details](GOLDEN_GATE_POOL_SETTINGS.md). This does not close combinatorial source
bins or the remaining interaction recommendations.

## Prioritized actionable work

**P0 — Keep the ordinary route visually continuous.** The recommended interaction is method → backbone/insert capture → preparation → junction/product review → named outputs. This is our adaptation of the documented pattern. Present the active slot, its source/direction and preparation status together; let an error take the user to the existing source/preparation control. Treat discoverability as a later focused user-review question, not a demonstrated runtime defect. The existing diagnostics table, junction controls and actual-product viewer should supply the information.

**P0 — Finish persistence of editable reaction choices.** `Draft` already retains `MethodsState` and explicit junction edits. However, `PreparedPool` constructs a component-local `PoolForm` containing reaction mode, substrate choices, topology and enumeration bounds; these fields are outside that serialized draft. Preserve these editable choices across draft save/reopen and navigation. Reopened scientific evidence must still be recomputed from exact inputs; persisting settings must not restore a fresh success state. This is a concrete remaining connection, not a request to rebuild draft persistence generally.

**P0 — Preserve deliberate manual work.** Our explicit preview/Apply and selected-row regeneration should remain the policy when a primer or junction changes. Explain which retained preparation becomes incompatible and which rows would be replaced. Benchling's documented replacement warning is useful interaction evidence; automatic replacement is not required for compatibility. The relevant algorithms already exist here.

**P1 — True combinatorial design after the ordinary route.** Current `Slot` stores one captured source, and the UI's bounded rows describe one ordered design. Enumerating Golden Gate reaction alternatives is different from a bin containing interchangeable promoter/CDS/backbone choices. Add alternative source choices, explicit selected combinations and per-construct review/output identity as a separate design layer over the shared preparation compiler. Avoid equating more reaction candidates with a completed combinatorial cloning workflow. Do not adopt enterprise tables or Benchling's advertised numerical limits as requirements without our own resource design.

**P1 — Fusion-end design and internal-site domestication.** Explicit Golden Gate preparations, fusion definitions and junction editing are present. Automated choice/ranking of fusion ends and intentional removal of internal sites remain separate design capabilities in the plan. Build them as inspectable proposals with explicit edits and fresh preparation/reaction reconstruction. Public help does not specify Benchling's scoring algorithm sufficiently to reproduce it, and its historical efficiency language does not supply a validated scoring dataset.

**P1/P2 — Broader physical reaction coverage.** Complete-cleavage, finite-state, one-copy bounds remain explicit in the current Golden Gate model. Partial/kinetic cycling and repeated-copy mixtures require additional scientific work if pursued. They are independent of matching Benchling's documented wizard interaction; do not block ordinary source selection and review on implementing them.

## Keep projects simple

Use open documents plus a named assembly draft as the normal entry point. Keep the existing source viewer central, ordinary two-slot cloning the default, and advanced multi-part/combination controls progressive. Create a construct with optional oligos and visible preparation documents in the current workspace; retain hidden preparation/history as the existing core requires. Export the same reviewed design to sequence, primer sheet and HTML/PDF steps.

This is a product recommendation, not a claim that Benchling has this exact simplified mode. Project folders, studies, registration schemas, worklists, ELN and inventory from its enterprise finalization flow are outside our scope. They are unnecessary to reproduce the useful source → design → review → output interaction.
