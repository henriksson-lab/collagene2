# Selected-reference CRISPR extensions: research inputs

Research date: 2026-11-28. This is a bounded implementation proposal, not a claim that these extensions or empirical scores are implemented. Current scientific behavior remains documented in [CRISPR.md](../CRISPR.md): SpCas9 NGG/20, declared Cas12a TTTV/23 planning geometry, custom geometry, and bounded Hamming search over explicitly selected records or a streamed manifest.

## Workflow evidence

Benchling's documented workflow selects a genome/transcript and guide parameters, highlights a sequence region, displays a sortable guide table, and saves selected guides as oligos. A guide's aggregate off-target score opens its individual candidate sites. Selected guides can enter vector assembly after choosing the vector and insertion region. Its donor workflow edits a sequence copy, reviews edits, adjusts homology arms, then considers silent target/PAM changes. Table export is available. These are useful interaction patterns without adopting its folder/registry/notebook machinery. [Benchling: Use CRISPR tools](https://help.benchling.com/hc/en-us/articles/37748593861133-Use-CRISPR-tools).

Keep aggregate specificity and individual-site activity labels distinct: Benchling documents opposite desirable directions and acknowledges that a generic “off-target score” label causes confusion. Collagene's selected-reference results must additionally expose the actual searched records, mismatch bound and truncation/completion state; they cannot imply a genome-wide assessment. [Benchling score clarification](https://help.benchling.com/hc/en-us/articles/10682547899021-Are-higher-or-lower-CRISPR-gRNA-off-target-scores-better).

## Small added nuclease roster

All sequences below are written on the **PAM-bearing, non-target DNA strand**, 5′→3′: `protospacer | PAM`. The guide matches that protospacer with RNA U replacing DNA T; it hybridizes to the opposite strand. PAM is excluded from guide length. Both reference orientations must be searched, including circular-origin windows. An expression-added 5′ G is a separate oligo/expression choice, not an extra matched genomic base.

| Proposed distinct preset | Matched guide | Design PAM and side | Nominal top/bottom offsets from oriented protospacer's 5′ boundary |
| --- | --- | --- | --- |
| SaCas9, 21-nt guide | 21 | `NNGRRT`, 3′ | `18 / 18` |
| SaCas9, 22-nt guide | 22 | `NNGRRT`, 3′ | `19 / 19` |
| SpCas9-VRER | 20 | `NGCG`, 3′ | `17 / 17` |

**Published facts:** Ran et al. tested SaCas9 guide lengths and found that NNGRRT supports its strongest genomic targeting, while other NNGRR PAMs can also be cleaved and merit off-target consideration. The authors' deposited SaCas9 manual recommends 21- or 22-nt guides and explicitly notes NNGRR(N) recognition. These two fixed-length presets are useful choices, not a claim that other lengths are inactive. [Ran et al., 2015](https://pmc.ncbi.nlm.nih.gov/articles/PMC4393360/); [author SaCas9 manual](https://media.addgene.org/data/plasmids/61/61592/61592-attachment_iAbvIKnbqNRO.pdf).

Kleinstiver et al. identify SpCas9-VRER as D1135V/G1218R/R1335E/T1337R, characterize preference for NGCG, and evaluate targeting with 20-nt spacers. Activity against other PAMs is not a universal zero. Their VQR variant recognizes a broader, uneven PAM spectrum, so it is a less clean first addition. [Kleinstiver et al., 2015](https://pmc.ncbi.nlm.nih.gov/articles/PMC4540238/).

**Planning convention/inference:** the offsets above deliberately declare a `guide_length − 3` blunt boundary on each strand. They are nominal rendering/planning inputs, not experimentally guaranteed product termini for each sequence, assay or engineered variant. The underlying SpCas9 cleavage convention is supported by [Jinek et al., 2012](https://pmc.ncbi.nlm.nih.gov/articles/PMC6286148/); applying a fixed pair of boundaries to these added presets is an explicit software approximation. Do not silently convert these marks into verified end chemistry or deterministic ligation products. A versioned source-bound preset must retain that qualification.

**Implementation consequence:** design admission PAM and selected-reference screening PAM need separate explicit policy. For SaCas9, screen at least `NNGRRN` when reporting the broader published scope; filtering only NNGRRT would miss acknowledged alternatives. The result must record which policy was searched. The current single-PAM representation cannot imply that broader screening happened. Preserve strand-aware cut mapping, immutable source fingerprints and independent candidate verification. Nickases, paired-guide rearrangements, base/prime editing and variant-specific empirical models are separate features, not consequences of adding these presets.

## Empirical data candidate and license boundary

The strongest next **on-target evaluation dataset candidate** is the authors' new validation experiment accompanying *Accounting for small variations in the tracrRNA sequence improves sgRNA activity predictions for CRISPR screening* (2022), Supplementary Data 4/Source Data. It measures guide abundance changes in triplicate A375 Cas9 dropout screens using Hsu, Chen and DeWeirdt tracrRNA variants. These are growth/guide-depletion observations, not direct per-molecule cleavage probabilities. The study shows that tracrRNA choice changes predictive performance. Its article states CC BY 4.0, subject to material-specific credit exceptions. Use only the authors' newly generated data after checking the exact supplement's notices; do not assume all collated third-party training data share that license. [Primary article, data availability and license](https://pmc.ncbi.nlm.nih.gov/articles/PMC9448816/).

No table, model weight or external implementation has been imported. Admission still requires the exact supplement file, publisher URL, content digest, per-file notices, column/coordinate dictionary and assay metadata. Retain attribution and CC BY notice separately from Collagene's MIT code license. Preserve tracrRNA, cell line, replicate, target/gene and measured endpoint; split evaluation by gene and duplicate spacer to avoid leakage. A locally implemented empirical model would need its own versioned derivation and held-out evaluation. This proposal does not promise a Rule Set 3 reproduction.

For comparison, the authors' [Azimuth repository](https://github.com/MicrosoftResearch/Azimuth) is archived and publishes a [BSD-3-Clause software license](https://raw.githubusercontent.com/MicrosoftResearch/Azimuth/master/LICENSE.txt). That establishes a software notice, not a verified dataset-specific grant for every underlying measurement or the separate CFD mismatch/PAM tables. The [2016 primary paper](https://www.nature.com/articles/nbt.3437) is relevant scientifically; no CFD tables are admitted by this audit. Resolve their exact data provenance/license before redistribution.

An optional score must declare its applicable nuclease, guide/PAM/context lengths, scaffold, assay and endpoint. Missing context or unsupported SaCas9/Cas12a/VRER inputs must produce “not applicable,” never a transplanted SpCas9 score. Do not label empirical rankings as expression, cleavage, repair or success probabilities. Geometry and Hamming counts remain independently useful without such scores; neither proves experimental safety or complete off-target coverage.

## Bounded next implementation batch

1. Add the three source-bound nominal presets and direct orientation/origin/cut-coordinate fixtures, keeping the explicit planning qualification.
2. Extend one current profile/search representation with explicit intended versus screening PAM policy; update callers directly and record searched policy in results/exports.
3. Expose the selected-region guide table, site drilldown and selected-guide-to-assembly handoff using existing immutable scientific captures. Keep donor edits explicit and independently verifiable.
4. Curate the single authors-generated validation supplement as an offline evaluation artifact only after file-specific license review. Empirical score implementation remains a separate evidence-gated batch; no external executables or runtime genome parsing are required by this proposal.
