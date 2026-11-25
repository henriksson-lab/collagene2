# PCR and primer workflow check

Reviewed 25 November 2026 against official public Benchling documentation and the current Collagene source, including the changes through `e12f93f`. This is a source and documentation review; no local browser checks, builds or runtime acceptance were performed. It supersedes stale gap statements rather than reopening implemented work.

## Documented reference behavior

| Official source | Supported behavior | Limits of the evidence |
|---|---|---|
| [Use primers](https://help.benchling.com/hc/en-us/articles/39072057302541-Use-primers), updated 2 October 2025 | A template selection creates a directional primer; an editor handles named singles/pairs, another selected domain, overhangs and restriction-site additions. Existing oligos undergo a configurable binding search before attachment. Linked pairs are reusable; selected primers export as CSV. Calculation settings and structure review are exposed. Its automatic wizard uses Primer3. | The guide does not specify exact mismatch enumeration, completeness guarantees, circular-product geometry, or experimentally validated PCR success. Its structure calculations do not establish equivalence with our simpler screens. |
| [Run PCR in silico](https://help.benchling.com/hc/en-us/articles/9684267245325-Run-PCR-in-silico), updated 20 October 2022 | Select/link two primers from the map or panel, choose the pair, request a PCR product, review options, then save the resulting sequence. | This describes interaction, not a precise physical-end or mismatch model. |
| [PCR and Primer Design](https://help.benchling.com/hc/en-us/articles/9684234653837-PCR-and-Primer-Design), updated 11 September 2024 | The teaching workflow starts with selected forward/reverse domains, links them, displays expected product size and creates a separate product. Primer length, GC content and melting temperature aid review. | Its design guidelines are teaching advice, not hard acceptance criteria or an amplification guarantee. Notebook instructions are outside our scope. |
| [Golden Gate Assembly Wizard](https://www.benchling.com/blog/golden-gate-assembly-wizard), published 16 December 2014 | Historical documentation connects fragment selection, generated primers, product review, primer-covered base edits, parent mismatch labels, undo, primer export, summary printing and finalization. | This is historical evidence for connected mutation/primer interaction. It does not establish a current dedicated inverse-PCR wizard, a particular kit mechanism or terminal chemistry. |

Searches of the official help site for inverse PCR and mutagenesis did not establish a current dedicated inverse-PCR procedure. That is a documentation limitation, not evidence that Benchling lacks one. The historical Golden Gate workflow must not be generalized into undocumented support for every mutation architecture.

## What Collagene already implements

These findings come from repository source and are independent of the public preview deployment.

| Capability | Current source evidence | Assessment |
|---|---|---|
| Shared directional editing | [primer_pair_editor.rs](../crates/web/src/primer_pair_editor.rs), [primer_pair_model.rs](../crates/web/src/primer_pair_model.rs) | Either slot accepts forward/reverse-complement selection or saved oligos; names, annealing domains and tails are explicit. Captured-domain comparison marks mismatches and ambiguity. Apply remains explicit. |
| Unfinished pair persistence | [primer_pair_storage.rs](../crates/web/src/primer_pair_storage.rs), attachment load/save wiring in [app.rs](../crates/web/src/app.rs) | Workspace persistence exists. The older recommendation to add it is closed. Restored text is not restored positive binding evidence. |
| Binding and basic properties | [primers.rs](../crates/web/src/primers.rs) | Binding options include allowed mismatches and required 3′ matches. Source-bound analysis, declared Tm conditions, simple structure screens, saved annotations and selected pair transfer already exist. |
| Ordinary PCR review | [pcr_outputs.rs](../crates/web/src/pcr_outputs.rs) | Multiple returned geometries remain selectable, with paged review, source geometry, named oligos, optional product/primer outputs and CSV/TSV export. A sole result may be preselected; multiple results require a choice. |
| Mismatch PCR review | [mismatch_pcr_outputs.rs](../crates/web/src/mismatch_pcr_outputs.rs) | Review includes returned geometries, template positions, primer-derived substitutions, paged reference letters/findings, search/model completeness and optional named reference/primer outputs. Fresh reconstruction precedes installation. |
| Primer orders | [primer_orders.rs](../crates/web/src/primer_orders.rs) and the two PCR output components | Selected oligos can become collections and CSV/TSV sheets with complete ordered 5′→3′ sequence and domain/tail decomposition. No JSON reconstruction is required. |
| Mutation intent and inverse PCR | [mutagenesis.rs](../crates/web/src/mutagenesis.rs), [mutation_review.rs](../crates/web/src/mutation_review.rs), candidate handoff in [app.rs](../crates/web/src/app.rs) | Literal source-bound mutation intent creates a separate intended reference. Inverse PCR reaches candidate review with exact original source events, independent-target checking and the existing shared output/order route. Earlier claims that event handoff or outputs are absent are stale. |

## Overlay implementation follow-up

The applied pair/selected geometry connection is now implemented: captured-domain
comparisons and actual selected exact/mismatch PCR products feed transient map and
sequence overlays, with full source navigation and explicit added tails and
uncertainty. Foreign-source inspection preserves the independent target.
[Details](PCR_OVERLAYS.md). The first recommendation below remains research history;
endpoint/tail-boundary editing is closed by the next checkpoint and local recovery
remains follow-up work. Browser
acceptance has not been run.

## Endpoint-editing implementation follow-up

The existing shared editor now previews endpoint and tail-boundary changes, preserves
literal edits/ordered DNA according to the operation, captures fresh exact-event
comparison geometry and requires explicit acceptance before pair Apply. Circular
anchors have a versioned migration path. [Implementation](PRIMER_DOMAIN_EDITS.md).
The second recommendation below is now implemented as explicit coordinate/boundary
controls; pointer-drag manipulation remains an optional interaction refinement.
Local failure recovery and output-presentation recommendations remain follow-ups.
The [parallel follow-up decisions](BENCHLING_WORKFLOW_DECISIONS.md) refine recovery
into literal-preserving reference capture, selected-hit handoff and full binding-list
paging. Those controls and local edit/policy/rerun actions are now implemented
([details](PRIMER_REFERENCE_RECOVERY.md)). Recommendations1/2 below are retained as
completed research history; broader output continuity remains follow-up work.

## Actionable priorities

The following are our recommendations from the comparison, not claims about undocumented Benchling internals.

1. **P0 — Bring the applied pair, binding placements and selected product together on the sequence.** Existing sequence tracks render saved primer annotations; draft comparison and PCR product geometry are currently reviewed through separate components. Add a transient source-bound overlay for the applied pair and selected geometry, with direction, annealing domain, introduced tail/mismatches and amplified interval distinguished. Navigate from either oligo or product row to the same placement. Preserve all alternatives and incomplete-search status; showing one geometry does not establish uniqueness.
2. **P0 — Edit annealing endpoints and tail boundaries with explicit preview and Apply.** Current text and selection actions are functional; inline boundary manipulation is the remaining interaction improvement. Preview both ordered oligos, anchor consequences and resulting product changes before changing an applied pair. An edited domain length must obtain a new explicit anchor rather than silently inheriting the old comparison. Preserve names and intentional mismatches.
3. **P1 — Make binding failure recovery local to the pair.** Distinguish invalid domain text, no placement under the declared 3′/mismatch policy, multiple placements, incomplete search and independent-target mismatch. Offer the relevant domain, region or parameter edit beside that outcome, then rerun explicitly. This improves connections between existing controls and evidence; it does not require a new primer engine or automatic relaxation of criteria.
4. **P1 — Unify the output review presentation across exact PCR, mismatch PCR and inverse-PCR candidates.** The underlying output and order routes already exist. Standardize where names, complete ordered oligos, optional saved artifacts and reference/verified-product distinctions appear. Keep workflow-specific scientific contracts: an ordinary simulator-derived reference must not be relabeled an independently verified target to enter the assembly finalizer.
5. **P1 — Put interpretation of primer bases beside the picture.** Persistently distinguish tail from annealing domain, mark intentional mismatches and ambiguous comparisons, and show whether properties apply to the annealing domain or whole ordered oligo. Existing calculations and order tables provide this data. A compact visual legend can reduce the need to infer it across panels.

## Scientific boundaries and deferred work

Retain the current independent-target check for designed mutations, original event identities, explicit source-change invalidation, full ordered oligos and uncertainty reporting. A captured reference comparison is not a specificity search; template-local binding is not genome-wide specificity; a simple complementary-run screen is not a folding/ΔG calculation. Product letters alone do not establish end chemistry or successful amplification.

Primer3, advanced thermodynamic ports and external scientific executables remain deferred or excluded as specified in [PROJECT_PLAN.md](PROJECT_PLAN.md). This audit does not introduce ELN, stock/inventory, vendor-order submission, organizational project structure or automation. No additional parser, planning engine or primer editor is justified merely to reproduce the documented workflow; improve the existing connected surfaces first.
