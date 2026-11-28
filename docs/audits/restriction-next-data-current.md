# F21 next data admission: current source

Reviewed 2026-11-28. This report supersedes the **proposed**, unimplemented status of the host roster in [the earlier research](restriction-source-presets-current.md); its biological cautions still apply. This batch changes only the independent enzyme catalogue and its tests/docs. Parent-owned core, UI and MCP integration is a separate part of the shared batch. No compiler or scientific tests were run by this worker.

## Actual gap and admitted boundary

The [delivery queue](../DELIVERY_AUDIT.md#current-implementation-queue) correctly identifies missing source-host presets and sparse formulation data. Before this batch, the runtime catalogue had 40 enzyme geometries, four NEB formulations and 32 condition observations, with no host owner. The existing shared reaction editor already retains explicit formulation, buffer, temperature, incubation minutes, activity threshold and assay contexts; the source-generation controls expose independent Dam/Dcm/CpG assumptions. Those controls are present functionality, not a missing workflow.

The current `2026-11-28.2` catalogue, schema **3**, adds eight documentary host records and four product formulations: **21 sources, 40 geometries, eight formulations, 64 condition observations**. The geometry set is unchanged. A selected host is provenance about a supplier product, never a measurement of a DNA preparation or an instruction to change its positional modifications. All eight host records leave CpG unknown.

## Host facts admitted

| Exact preset ID | Supplier product | Dam / Dcm expectation | Primary evidence |
| --- | --- | --- | --- |
| `neb-5-alpha-c2987` | NEB 5-alpha C2987 | Active / Active | NEB selection-guide phenotype row; separate cloning-cell product identity table. |
| `neb-10-beta-c3019` | NEB 10-beta C3019 | Active / Active | Same primary tables; DH10B-derived product, without a universal DH10B alias. |
| `neb-stable-c3040` | NEB Stable C3040 | Active / Active | Explicit phenotype row; a stability-oriented name does not imply methylase deficiency. |
| `neb-dam-dcm-c2925` | NEB C2925 | Deficient / Deficient | Phenotype table and product genotype `dcm-6`, `dam13::Tn9`. |
| `neb-bl21-de3-c2527` | NEB BL21(DE3) C2527 | Active / Deficient | Explicit mixed phenotype; this adds source documentation, not expression workflows. |
| `agilent-jm110-200239` | Agilent JM110 200239 | Deficient / Deficient | Product specifications and mutant `dam dcm` genotype. |
| `thermo-inv110-c717103` | Invitrogen INV110 C7171-03 | Deficient / Deficient | Manufacturer manual MAN0000103, revision A.0, page 1; online code C717103. |
| `thermo-top10-c404003` | Invitrogen TOP10 C404003 | Active / **Unknown** | Manufacturer TOP10/Dam troubleshooting statement and exact product identity. Dcm is not inferred from a genotype omission. |

The NEB rows were verified in the current [online selection guide](https://www.neb.com/en/tools-and-resources/selection-charts/competent-cell-selection-guide), including the separate M.EcoKI column, with product identities from [the cloning-cell table](https://www.neb.com/en-us/products/competent-cells/cloning-competent-cell-strains). C2925's [product page](https://www.neb.com/en/products/c2925-dam-dcm-competent-e-coli) supplies its two mutant alleles and cautions against primary cloning/ligation in dam-deficient hosts. The independent supplier records come from [Agilent JM110](https://www.agilent.com/store/en_US/Prod-200239/200239), the [INV110 manufacturer manual](https://documents.thermofisher.com/TFS-Assets/LSG/manuals/oneshotinv110_man.pdf), and [TOP10 troubleshooting](https://www.thermofisher.com/us/en/home/technical-resources/technical-reference-library/cloning-technical-support-center/comp-cells/comp-cells-troubleshooting.html) plus [C404003](https://www.thermofisher.com/order/catalog/product/C404003).

M.EcoKI activity in C2925 is retained as a qualification: Dam/Dcm deficiency is not absence of every modification. Supplier genotype/phenotype is distinct from actual occupancy, propagation history, mixed material and subsequent treatment. NEB documents incomplete modification in some substrates and loss of eukaryotic CpG patterns through bacterial cloning; neither observation establishes every site's state in an imported plasmid. [NEB methylation guide](https://www.neb.com/en/tools-and-resources/selection-charts/dam-dcm-and-cpg-methylation).

## Formulation facts admitted

Each row below retains its own product identity, recombinant-albumin lot transition, four-buffer observations, recommendation and methylation categories. Percentages are approximate supplier unit-assay facts, not cleavage probabilities, reaction yield, substrate-independent activity or complete-mixture compatibility.

| Added formulation | Incubation / four-buffer activity, r1.1 · r2.1 · r3.1 · rCutSmart | Primary source |
| --- | --- | --- |
| `neb-r0136` BamHI | 37 °C; 75 · 100 · 100 · 100%; star warnings in r1.1/r2.1/rCutSmart retained | [R0136 properties](https://www.neb.com/en/products/r0136-bamhi), rAlbumin lot 10184085 onward |
| `neb-r3136` BamHI-HF | 37 °C; 100 · 50 · 10 · 100%; reduced star activity does not mean none | [R3136 properties](https://www.neb.com/en-us/products/r3136-bamhi-hf), lot 10133983 onward |
| `neb-r3733` BsaI-HFv2 | 37 °C; 100 · 100 · 100 · 100%; qualified overlapping Dcm/CpG effects retained | [R3733 properties](https://www.neb.com/en-gb/products/r3733-bsai-hf-v2), lot 10156743 onward |
| `neb-r0739` BsmBI-v2 | 55 °C unit-assay context; <10 · 50 · 100 · 25%; distinct 37 °C activity observation does not define these buffer values at 37 °C | [R0739 properties](https://www.neb.com/en-us/products/r0739-bsmbi-v2), lot 10284618 onward |

The temperature association follows each product's stated unit definition alongside its buffer table; no interpolation is made. In particular the BsmBI-v2 source also reports 10% activity at 37 °C without a full four-buffer matrix at that temperature. The catalogue does **not** manufacture a `rCutSmart / 37 °C / 25%` fact from it. Current assessment therefore returns unknown for that uncurated condition. A 37 °C Golden Gate cycle does not inherit the 55 °C single-enzyme unit-assay result. The [supplier chart's notes](https://www.neb.com/en/tools-and-resources/usage-guidelines/nebuffer-performance-chart-with-restriction-enzymes) likewise qualify buffer percentages by assay substrate and distinguish star warnings from activity values.

The [BamHI comparison FAQ](https://www.neb.com/en-us/faqs/how-does-the-level-of-star-activity-of-bamhi-hf-compare-to-bamhi) supplies a useful narrowly scoped star observation: 1× buffer, 5% glycerol, λ DNA 1 µg/50 µL, 1 hour at 37 °C; reported amounts are 32 units for BamHI/r3.1, 4 units for BamHI/rCutSmart, and **≥4,000** units for HF/rCutSmart. The qualifier preserves all those dimensions, censoring and unspecified assay lot. It does not claim that the observation was measured on a current rAlbumin lot or make it a universally safe enzyme-dose threshold.

## Current API and integrity

`hosts::HostPreset` owns exact ID, name, species, vendor/product code, product scope, Dam/Dcm/CpG `MethylaseExpectation::{Active,Deficient,Unknown}`, qualifications and citations. `Catalog::hosts()` borrows records; `host(id)` accepts only an exact installed ID; `validate_host(record)` checks bounds and **full current value membership**, including citations. An ID, URL, altered phenotype or caller-computed hash cannot replace installed facts. There is no host-to-`FullyModified` conversion and no automatic demethylation.

The core's cohesive `chemistry::source_host::SourceHost` additionally retains the catalogue edition/content pin, preset and **full cited source records together**. `from_current(id)` captures only installed facts with bounds before deep copies; `validate()` refuses old pins, altered/missing/extra/duplicate sources and noncanonical source order. One immutable cached catalogue serves repeated generation-record validation. This closes the bare-preset provenance gap: historical source metadata is not silently replaced with current metadata. It remains independent of every editable chemical generation assumption. Five dedicated core owner tests are authored and unrun.

`Edition.hosts` is required in the one current schema. Sources hash their cited complete host records as curated assertions; these are integrity hashes, not authenticated vendor snapshots. Native source sealing was performed by the parent using the existing offline authoring tool. The compact dataset is **58,951 bytes**, SHA-256 **`086c3a8b2926ca31a9e27362406b47fb91ef4d9496e0f5fda2a8120d2637bd92`**. All 21 source fact digests are populated; no zero placeholder remains. Parent owns catalogue-dependent scientific pins and cross-crate schema updates.

Dedicated unrun tests cover the eight literal phenotypes, unknown CpG/Dcm, current-value membership/tampering, source checksums, required current schema, native/decode bounds, Postcard, the four literal activity matrices, BsmBI temperature separation and exact star-assay qualification. Existing unsupported-version tests were updated directly to reject 1/2/4 around current schema 3.

## Remaining useful work and license boundary

1. Integrate the documentary host choice into exact source-generation declarations, serial material provenance, controls and exported protocols; preserve independently editable modification assumptions. Parent owns that connected work. Selecting a host must not mutate chemical state or revive an old verification result.
2. Add product-specific data for common remaining selected geometries, especially HindIII, XbaI, XhoI, NotI-HF, BbsI-HF, Esp3I and SapI, after the same source/lot/context review. “40 geometries” must not appear as “40 products with known compatible conditions.”
3. A general typed star/duration assessment still needs declared enzyme quantities, stock glycerol fractions, final volume, other glycerol sources, substrate and matched assay predicates. Current `ReactionConditions` lacks those quantities. Preserve documentary warnings and unknowns; do not implement a probability or universal >5%/long-incubation rule. Golden Gate's ligase/buffer/cycling mixture remains a distinct unassessed scope.
4. Do not fill TOP10 Dcm or generic DH5α/DH10B aliases from memory. Admit stronger product-specific primary evidence when found; preserve unknown until then.

Code and independently authored fact selection remain MIT. No vendor HTML, images, table layout, manual or scientific implementation is redistributed. No blanket open-data grant was found: [NEB Terms of Use](https://www.neb.com/en-us/terms-of-use), revised 2023-08-02, expressly restrict reuse of website content; supplier documents retain their notices. The catalogue's minimal factual encoding does not purport to license vendor documents or trademarks. Permission for copying vendor prose, bulk compilations or images is not established by this work. Raw snapshot fields remain absent for newly reviewed sources because exact payload bytes were not captured; an assertion digest is never presented as an upstream-document digest.


Parent integration checkpoint: the shared native/WASM batch and 349 unique focused
cases pass. See [source and execution evidence](cloning-extensions-current-integration.json)
for the exact scope, initial catalogue-pin corrections and final metadata-only
wording edit. No graphical or application HTTP acceptance is claimed.
