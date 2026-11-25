# PCR and primer output handoff: focused workflow review

Reviewed 2026-11-25 using official public Benchling documentation and repository source. No local browser checks, compilation, tests, third-party source code or authenticated Benchling session were used. This updates the output-handoff gap in [the primer/validation review](BENCHLING_PRIMER_VALIDATION_WORKFLOW_REVIEW.md); it does not reopen completed shared primer-editor work. Source inspection is not GUI acceptance.

Subsequent implementation: ordinary exact-PCR now retains all returned geometries
for explicit selection, naming, optional product/oligo output and CSV/TSV. Fresh
worker reconstruction and exact original-event guards precede atomic staging.
The observations below record the pre-change audit; the ordinary-PCR handoff is
no longer wholly absent. [Current contract](PCR_OUTPUTS.md). Inverse-PCR now also
preserves exact source/target events through candidate review and reaches the
shared actual-product/output/order route; [current output contract](CLONING_OUTPUTS.md).
The comparison below records the original audit, rather than an unchanged backlog.
Mismatch-PCR and standalone primer-analysis follow-ups remain active.

## Documented Benchling behavior

**Primer selection precedes saving and export.** Benchling's primer guide describes selecting template bases and direction, editing a named single primer or pair, and saving. Wizard results let users choose primers before saving or CSV export. Attached sequence primers also have checkboxes and a bulk CSV action. These are two selection contexts: choosing a designed candidate and exporting already attached primers. A named oligo can be reused across templates; a primer pair defines the intended amplification region. The documented wizard uses Primer3, which remains outside our implementation scope. [Use primers, updated October 2, 2025](https://help.benchling.com/hc/en-us/articles/39072057302541-Use-primers).

**PCR creation is an explicit product action.** Users link a pair, select it in the primer-pairs panel, request a PCR product, review output options, then copy/save the resulting sequence. The documentation does not describe selecting among competing amplification geometries or establish the algorithm's physical-end and mismatch assumptions. We should copy the clear interaction progression without inventing equivalent scientific guarantees. [Run PCR in silico, updated October 20, 2022](https://help.benchling.com/hc/en-us/articles/9684267245325-Run-PCR-in-silico).

**Assembly finalization separates output categories.** Benchling exposes construct preview and a separate primer table before finalization. Finalization requires constructs and offers optional primers and amplified fragments. It subsequently stores assembly history. This supports a consistent output review for primer-led cloning, but the public guide does not state that standalone PCR uses the identical finalization component. Folder, registry, notebook and worklist organization are unnecessary for our scoped adaptation. [Combinatorial cloning and concatenation, updated December 12, 2025](https://help.benchling.com/hc/en-us/articles/39922419457165-How-to-use-the-combinatorial-cloning-and-concatenation-assembly-tools).

A focused official-site search did not establish a dedicated current inverse-PCR/mutagenesis procedure. Historical junction-editing articles are insufficient to infer its current physical preparation or output semantics. Our inverse-PCR route should remain explicitly PCR, declared phosphorylation and ligation, independently checked against a separately specified target.

## Current source comparison

| Route | Existing behavior | Remaining handoff |
| --- | --- | --- |
| Ordinary exact PCR | `app.rs::do_pcr` validates named domains/tails, sends `PcrDocuments`, binds the response to the source event, draft inputs, workspace epoch and generation, then atomically derives and opens **every** returned product document. Workspace limits are checked before installation. | Retain results for explicit product/geometry selection and naming before installation. Offer the chosen pair as a reusable primer collection and CSV/TSV order reference without reconstructing it from evidence JSON. Preserve existing source-event checks and atomic limits. |
| Primer analysis/design | `primers.rs` can retain a primer in either shared draft slot, save source-bound primer annotations, and transfer a selected pair to PCR. The shared pair editor already supports directional selection, manual and saved oligos. JSON evidence downloads exist. | Add selected named oligo export/save actions using existing oligo tables. Do not build another editor or claim Primer3-equivalent ranking. |
| Mismatch PCR | `mutagenesis.rs` displays bounded product geometries, opens an explicitly chosen reference preview, and retains full reported evidence in JSON. The previews explicitly withhold chemistry and experimental-success claims. | Reuse product-selection/output controls while preserving the conditional mismatch model and all search/truncation information. Do not convert a reference preview into verified cloning output. |
| Inverse PCR | `mutagenesis.rs` retains the plan, independent-target verification and full ordered oligos, and sends a candidate to the general plan panel. That panel supports fresh verification, candidate download, reconstructed-product retention, actual-junction inspection and HTML/PDF protocols. | The `MutagenesisPanel` callback in `app.rs` explicitly clears `candidate_events`; it does not forward exact source-event preferences. General plan review does not mount the shared optional-output or order-export components. Preserve the current review route, add exact original event handoff, and expose the shared final output choice from this context. |
| Assembly dock | `assembly_outputs.rs` exports plan oligos, while `AssemblyFinalizePanel` handles named construct, optional primers/prepared references and atomic output commit. Its binding capture requires ordered exact source events and the active independent target. | Reuse these scientific/export/output builders. Current `finalize_assembly` accepts an assembly `ReviewRequest` and stages its intent; standalone PCR should not invent an assembly intent or disguise an inferred PCR preview as an independent expected target to reach it. |

Relevant implementation: [app.rs](../crates/web/src/app.rs), [primers.rs](../crates/web/src/primers.rs), [mutagenesis.rs](../crates/web/src/mutagenesis.rs), [assembly_outputs.rs](../crates/web/src/assembly_outputs.rs), [assembly_finalize.rs](../crates/web/src/assembly_finalize.rs), [plan oligos](../crates/core/src/plan_oligos.rs), [cloning outputs](../crates/core/src/cloning_outputs.rs).

## Recommended connected workflow

These are design recommendations from the comparison, not undocumented Benchling claims.

1. Keep one named applied pair and exact template event. Compute into a retained result instead of immediately installing all ordinary PCR products.
2. Show reported product length, binding geometry, original template interval, annotation consequences and model/search limitations. Require an explicit geometry when several products exist; never silently choose index zero. Keep all results reachable rather than adding an unexplained first-N cutoff.
3. Let users select outputs: a named PCR reference preview, the declared primer pair as a collection, and CSV/TSV export. Export full ordered 5′→3′ bases including tails, retaining names and domain/tail decomposition in metadata. Export is an order reference, not stock management or proof of amplification.
4. Install selected documents in one preflighted transaction with event provenance. Draft, source-event, workspace or options changes invalidate result authority without deleting the editable pair. Preserve existing document-count, byte-budget and recovery validation.
5. For inverse PCR, keep the separately declared target and original source events through candidate review. Freshly reconstruct before finalization, then expose the existing construct/primer/prepared-reference choices and shared order export. A verification result, a saved candidate and an order sheet remain different artifacts.

The first task should be ordinary exact-PCR result selection plus named oligo output reuse. It closes an everyday cloning interruption without new scientific engines, while inverse-PCR can subsequently use the same presentation and provenance plumbing. Do not relax physical verification merely to fit reference previews into the assembly finalizer.
