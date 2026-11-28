# Restriction source presets and condition facts

Primary-source research: 2026-11-28. Proposed bounded inputs for the next F21 implementation; no new scientific code or catalog facts are imported here. A host's documented methylase phenotype is **not a measurement of the submitted plasmid's positional modification state**.

## Small supplier-scoped preset roster

These records describe host methylase expectations, not `FullyModified`/`Unmodified` calls. CpG generation state remains **Unknown for every preset**. Supplier products/derivatives are separate identities; do not silently alias every DH5α/DH10B variant to one product.

| Proposed preset identity | Documentary Dam / Dcm expectation | Primary source and exact scope |
| --- | --- | --- |
| NEB 5-alpha, C2987 | Active / Active | NEB *Competent Cells* brochure, printed p.10, “Cloning Strain Properties” methylation-phenotype column: Dam+, Dcm+, M.EcoKI+. The same brochure identifies it as a DH5α derivative. |
| NEB 10-beta, C3019 | Active / Active | Same p.10 column: Dam+, Dcm+, M.EcoKI−; described separately as a DH10B derivative. |
| NEB dam−/dcm−, C2925 | Deficient / Deficient | Supplier product describes plasmid growth free of Dam/Dcm methylation; genotype explicitly contains `dcm-6` and `dam13::Tn9`. The brochure still lists M.EcoKI+, so “Dam/Dcm-free” must never mean “all DNA modifications absent.” |
| Agilent JM110, 200239 | Deficient / Deficient | Product specifications explicitly say preparation of plasmid/phagemid DNA free from Dam/Dcm methylation; genotype lists mutant `dam dcm`. |
| NEB BL21(DE3), C2527 | Active / Deficient | NEB brochure p.10 “Protein Expression Strain Properties” lists Dam+, Dcm−, M.EcoKI−. Useful mixed-state control, without importing protein-expression workflows. |
| Invitrogen TOP10 | Active / Unknown in this captured source set | Thermo Fisher's competent-cell troubleshooting explicitly identifies TOP10 as Dam+. It does not establish Dcm status in the reviewed passage; do not fill that field from memory or absence of a genotype mutation. |

Sources: [NEB brochure, printed p.10 / PDF page 10](https://www.neb.com/-/media/nebuk/files/brochures/competent_cells_brochure.pdf?rev=1556823e86c14ba3a2a96fe8afcc9f4f#page=10); [NEB C2925, Product Information/Genotype](https://www.neb.com/en/products/c2925-dam-dcm-competent-e-coli); [Agilent JM110 specifications](https://www.agilent.com/store/en_US/Prod-200239/200239); [Thermo Fisher TOP10 troubleshooting](https://www.thermofisher.com/us/en/home/technical-resources/technical-reference-library/cloning-technical-support-center/comp-cells/comp-cells-troubleshooting.html). Brochure values were read from its parsed table; the web PDF screenshot endpoint failed, so no visual table-inspection claim is made.

For named DH10B itself, a [primary experimental publication's strain table](https://pmc.ncbi.nlm.nih.gov/articles/PMC3312211/) identifies Invitrogen DH10B as dam+/dcm+ with ΔhsdRMS. This is useful corroboration, not proof of every laboratory derivative's genotype. Zymo's [DH5 Alpha product FAQ](https://www.zymoresearch.com/products/mix-and-go-competent-cells-dh5-alpha) gives a general “most cloning strains” Dam+/Dcm+ statement; that is weaker than a product-specific phenotype assertion. Initial implementation should use the explicit NEB product records above rather than manufacture a universal DH5α alias. `hsdR/hsdM`, `mcr/mrr`, `endA` and Dam/Dcm describe different systems; a restriction-deficient genotype is not automatically methylase-free.

## Safe generation semantics

Published facts: NEB warns that methylation extent varies even among DNA from *E. coli*, contrasting fully modified pBR322 with approximately half of λ Dam sites modified. Dam targets GATC adenine; Dcm targets the second cytosine of CCWGG. NEB also distinguishes overlapping methylase/restriction contexts and says eukaryotic CpG patterns are not retained through bacterial cloning. None of this establishes the state of a particular uploaded molecule. [NEB Dam/Dcm/CpG guide](https://www.neb.com/en/tools-and-resources/selection-charts/dam-dcm-and-cpg-methylation).

Proposed API/UI contract:

1. A pinned preset retains `id`, edition/content identity, supplier/product, genotype/phenotype assertions, source URL/locator/review date and explicit Dam/Dcm expectations (`Active`, `Deficient`, `Unknown`). It does not itself call `TrackedDuplex::from_input` or overwrite retained chemistry.
2. Selecting a host only selects documentary provenance. If the user wants a planning assumption, require an explicit **generation declaration** bound to the exact original source revision. Keep its editable Dam/Dcm/CpG states visible. “Assume all eligible Dam/Dcm sites modified after this declared propagation” is a separate user assertion, never an automatic consequence of Dam+/Dcm+.
3. Negative genotype alone likewise cannot erase an incoming modified template. A proposed absence assumption must declare that this is newly propagated plasmid material from the stated host under the stated production history, without relevant additional methylase treatment. Retained positional evidence and mixed/unknown history take precedence; do not simulate demethylation merely by changing a label.
4. CpG remains unknown unless separately declared/evidenced. A human insert, host name or bacterial methylase genotype cannot set every CpG modified or unmodified. Never clear other modification classes outside the supported model.
5. Apply assumptions once to an original source generation. Subsequent supported cuts/ligation/repair/PCR use the existing positional propagation owner. Do not rescan a product to reconstruct methylation lost from a split motif.

The existing `GenerationAssumptions` and source-bound `GenerationRecord` can retain the explicit expanded states; preset provenance needs a deliberate current representation extension, not an adapter. See [tracked owner](../../crates/core/src/chemistry/tracked_duplex.rs#L25). Host propagation/growth is not an automated workflow claim. NEB specifically discourages Dam− hosts for primary cloning/ligation because of mutation/transformation concerns; retain that documentary qualification without predicting mutation rates. [NEB generation guide](https://www.neb.com/en/tools-and-resources/usage-guidelines/making-unmethylated-dam-dcm-dna).

## Reaction/star facts with exact applicability

| Fact to curate separately | Exact conditions/qualification; permitted interpretation |
| --- | --- |
| EcoRI R0101 | Product notes flag possible star activity specifically in NEBuffer r2.1 or rCutSmart. Recommended reaction is 1× EcoRI/SspI buffer, 37 °C; stock contains 50% glycerol. Do not treat the 100% r2.1 unit-assay activity as specificity approval. [NEB R0101](https://www.neb.com/en/products/r0101-ecori). |
| EcoRI-HF R3101 | Separate formulation: recommended 1× rCutSmart, 37 °C, 50% glycerol stock; engineered for reduced star activity. “Reduced” is not zero or a transferable numeric limit. [NEB R3101](https://www.neb.com/en/products/r3101-ecori-hf). |
| BamHI/BamHI-HF comparative assay | NEB reports first-star observations in **1× buffer, 5% glycerol, λ DNA 1 µg/50 µL, 1 h at 37 °C**: BamHI-HF in rCutSmart ≥4,000 units; BamHI in r3.1 32 units; BamHI in rCutSmart 4 units. Preserve substrate, formulation, buffer, time, temperature, units and ≥ censoring together. These are assay observations, not safe-dose thresholds for arbitrary plasmids or extrapolations to other volumes/times. [Exact supplier FAQ](https://www.neb.com/en-us/faqs/how-does-the-level-of-star-activity-of-bamhi-hf-compare-to-bamhi). |
| SmaI R0141 duration qualification | NEB permits either listed temperature but recommends 25 °C for extended digestions **>1 hour**. It is a temperature/duration recommendation, not a quantitative star-activity threshold. [NEB FAQ](https://www.neb.com/en/faqs/which-incubation-temperature-should-i-use). |
| Generic NEB warning | Its guideline lists **>5% v/v glycerol**, excessive enzyme, prolonged reaction, solvents and alternative divalent ions as contributing conditions. It supplies no universal probability or exact safe time/temperature. Its separate double-digest guidance says **less than 5%**, so retain original comparison and scope rather than normalize all sources into one universal cutoff. [Star guideline](https://www.neb.com/tools-and-resources/usage-guidelines/star-activity); [double digests](https://www.neb.com/en-gb/tools-and-resources/usage-guidelines/double-digests). |

Glycerol can be calculated only from explicitly known stock fractions and all added volumes: `Σ(stock_volume × stock_glycerol_fraction) / final_volume`. Two 1 µL stocks at 50% in 20 µL contribute exactly 5%; that is arithmetic, not proof that star activity occurs or is absent. Unknown stock fraction, glycerol in another reagent or changed final volume must remain unknown. The current conditions request has no stock-volume/unit/DNA-mass fields, so it cannot evaluate those assay predicates; keep them documentary until a coherent declared-quantity extension exists.

## Provenance and redistribution

No bulk supplier database, page image, genotype table or scientific implementation is copied into code by this audit. Small independently curated factual assertions need source locators, retrieval date, original qualification and a new content-pinned edition. Supplier pages carry their own copyright/trademark terms; this review found no explicit blanket open-data license for these pages. Collagene's MIT code license must not be presented as relicensing supplier material. Record that uncertainty; curate minimal facts with citations rather than redistributing full vendor tables/HTML. Source-content hashes, if captured later, bind an actual retrieved payload; an assertion hash must not masquerade as an upstream snapshot hash.
