# Beyond-competition feature research

Research baseline: 28 November 2026. Six parallel primary-source research agents
covered practitioner complaints, methods literature, community standards and
competitions, incumbent product gaps, long-read clone verification, and modern
modalities with AI-agent-native design.

Scope filter applied throughout: embedded Rust/WASM implementation only, no
external executables or cloud services, and the existing plan exclusions
(ELN, inventory, instruments, scheduling, reporting; biosecurity screening;
Primer3 ports; multiple-alignment engine ports). This document records research
findings and candidate features. It is not an implementation claim, and a cited
method is not evidence that our implementation of it is correct.

A recurring and important result: **the strongest opportunities are deterministic
checks, not predictions.** Published, enumerated, calibration-free checks are
both the biggest documented user need and the safest thing for this project to
ship. Predictive models mostly arrive with weak calibration, restrictive data
licensing, or both, and are recorded here with those limits attached.

## 1. The headline evidence: basic construct design is not a solved problem

Two independent agents converged on the same survey. Bai et al.,
*Prevalence of errors in lab-made plasmids across the globe*, Nucleic Acids
Research 53(14):gkaf697, 5 August 2025, [doi:10.1093/nar/gkaf697](https://doi.org/10.1093/nar/gkaf697),
open access (CC BY): design review identified significant errors in **384 of
2,521 plasmids**. Separate experimental cohorts had inconsistent restriction
patterns in 128/852 and sequence variations in 91/259 directly sequenced plasmids.
These denominators differ; the latter percentage must not be attributed to every
plasmid in the design survey or combined into an observed whole-cohort error rate.
Addgene's own QC reportedly finds errors in
[30% of deposits](https://blog.addgene.org/a-look-at-addgenes-qc-process).

Table 1 groups reported design-error categories. The study identifies 384
plasmids with design errors among 2,521 examined; this is not an exhaustive count
of distinct error instances or proof of 21 automatically decidable checks.
Reported examples include: one promoter driving two genes in a lentiviral vector (39),
insert too large for the vector (23), an extra stop codon upstream of a 2A peptide (31),
no stop codon upstream of an IRES (14), sense-orientation internal polyadenylation
signal in a lentiviral vector (7), Kozak absent or not adjacent to the ORF (17),
a CRISPR donor that does not mutate the PAM (9), a gRNA vector wrongly including
the PAM (7), shRNA sense/antisense flipped (7), and unstable high-GC or repeat
sequence in a high-copy origin (27).

Conflict of interest to record: all authors are affiliated with VectorBuilder,
which sells a competing design tool. The published taxonomy motivates checks, but each rule still needs explicit
inputs and scientific validation. Context-dependent biological effects may remain
indeterminate rather than become deterministic sequence verdicts.

This is the direct justification for the [F36 construct QC panel](../PROJECT_PLAN.md#f36-construct-qc-panel--user-decision-28-november-2026),
and the error table is its prioritisation guide.

## 2. Practitioner-reported gaps (written-down, attributable asks)

Reddit, SEQanswers, Biostars and Bluesky were unreachable from the research
environment, so this slice is weighted toward issue trackers of tools bench
scientists actually file against, and under-weights price and speed complaints.

| Need | Evidence | Echoes |
| --- | --- | --- |
| Declare the intended product, then check every assembly against it: frameshifts from primer placement, internal stops in fusions, **out-of-frame upstream ATG**, missing Kozak/RBS, protease sites created at junctions | [pydna #515](https://github.com/pydna-group/pydna/issues/515) (15 Jan 2026): graphical inspection in SnapGene or Benchling "is as error prone as the ad-hoc cut-and-paste cloning itself… Each of these problems might not be caught at the design level, as the result may seem correct"; [#395](https://github.com/pydna-group/pydna/issues/395) asks to "screen for conflicting start codons after 'Assembly'" | 2 threads + a named failure mode |
| Sweep an entire primer collection against a construct, mismatch-aware: every binding site, strand, footprint, mismatch count, bound-stretch Tm | [OpenCloning #83](https://github.com/OpenCloning/OpenCloning/issues/83) (8 Sep 2026), filed by a university RSE team for a lab: "given a plasmid and a box of primers, which of them bind, where, and would they work at my annealing temperature"; [#73](https://github.com/OpenCloning/OpenCloning/issues/73) | 3 independent asks |
| Persist the as-built/as-designed difference on the record rather than overwriting the map | [Bioinformatics SE 18801](https://bioinformatics.stackexchange.com/questions/18801) (NCBI's answer: `misc_difference` for lab-produced variation); pydna maintainer proposes HGVS | 3 |
| Non-destructive auto-annotation: merge with existing features and keep per-feature provenance, never erase | [OpenCloning #77](https://github.com/OpenCloning/OpenCloning/issues/77); [pLannotate #28](https://github.com/mmcguffi/pLannotate/issues/28) "keep the original annotations… just add the new ones without erasing the originals" | 5 |
| Enzymology-accurate assembly simulation: HiFi removes up to ~10 mismatched 3′ bases; mismatched primers still amplify; offer restriction-site regeneration | [pydna #686](https://github.com/pydna-group/pydna/issues/686) (28 Sep 2026); [OpenCloning #86](https://github.com/OpenCloning/OpenCloning/issues/86) | 4 |
| Screening-strategy designer: which digest or colony PCR distinguishes the intended clone from empty vector, wrong orientation, tandem multi-insert and parent carryover | [Biology SE 78144](https://biology.stackexchange.com/questions/78144) and siblings; encodes real traps (bands too close to resolve, "one primer inside your fragment and the other in the backbone… also verifies correct direction") | several |
| Score a hand-designed PCR, and re-fit the scoring on the lab's own recorded outcomes | [Lattice primers #4](https://github.com/Lattice-Automation/primers/issues/4): "want to quantify a 'good' or 'bad' design… potentially change the penalty score settings based on our actual lab data" | 1, but unique |
| Degenerate libraries as first-class objects, not 10⁵ files | [Bioinformatics SE 15662](https://bioinformatics.stackexchange.com/questions/15662): "I don't want to end up with >10⁵ different definitions to describe the sample" | 1 strong |
| Verification alignment tuned for sequencing reality: one semi-global hit with soft ends, not a fragmented gapped alignment | [OpenCloning #81](https://github.com/OpenCloning/OpenCloning/issues/81) | 2 |

Positioning context, not features: "ApE is old-school but works well. SnapGene is
too expensive, and Benchling is not very power-user friendly"
([HN, 12 Feb 2023](https://news.ycombinator.com/item?id=34759358)). Several
open-source trackers now answer feature requests with funding quotes, so the
maintained-free-tool gap is acknowledged by its own maintainers.

## 3. Incumbent gaps worth exploiting

Vendor help centres (`help.benchling.com`, `support.snapgene.com`,
`help.geneious.com`) return 403 to automated fetches; those quotes come from the
search index of the vendor's own page and are marked accordingly in the source
report. Product and release-note pages were read directly.

- **CRISPR scale.** SnapGene shipped CRISPR guide design only on 3 November 2026
  (v9.0.0), restricted to sequences **under 5,000 bp**, with efficiency and
  specificity scoring "on a cloud server" and specificity against "a selection of
  common genomes" ([release notes](https://www.snapgene.com/updates/snapgene-version-9-0-0)).
  Geneious Prime 2027.0 surfaces only "the 100 highest-scoring off-target sites"
  ([release notes](https://www.geneious.com/updates/geneious-prime-2027-0)).
  Collagene's implemented searches are offline and bounded, with explicit work,
  hit and reference limits and provenance. They must not be described as uncapped
  or as proof of genome-wide coverage beyond the actual searched inputs.
- **Verification is nobody's product.** No vendor documentation of nanopore
  read-to-construct QC was found in Geneious, SnapGene or Benchling. SnapGene's
  whole-plasmid path requires a circular reference, "aligns each read
  independently", does not assemble, and leaves pass/fail to a manual
  "Confirmed experimentally" checkbox. Geneious added circular alignment only in
  2027.0.
- **Cloud dependence.** Benchling "requires an internet connection to access and
  save data", with no offline mode, desktop version or local install. SnapGene
  requires internet for activation and licence checks. When Addgene began
  requiring login for sequence downloads in April 2026, SnapGene's Addgene import
  broke and its own docs still say it "is currently not working" seven months later.
- **Licence mechanics.** SnapGene academic $360/user/yr, corporate $1,945, with a
  $750 purchase-order minimum and a permanent licence that excludes upgrades;
  on expiry it reverts to Viewer mode where "most other functions are disabled".
  Geneious reserves its **public REST API** and **multiple insertion point
  cloning** for Enterprise. Those two are worth making headline free features.
- **Interop is lossy in quotable ways.** Benchling exports SBOL only from the
  single-sequence interface, with no bulk path. SnapGene imports `.geneious` files
  but "other file types such as alignments and trees are ignored". Geneious 2026.0
  added SnapGene **history** import, so history interop is now contested ground.
  A published round-trip conformance suite with per-field loss reports is a
  defensible differentiator and fits our existing declared-loss export machinery.
- **No collaboration primitives in SnapGene.** No commenting, review or approval
  feature appears in any release note from November 2023 to November 2026.
  Range-anchored comments and reviewer sign-off stored *in the construct file*,
  diffable over a shared folder, need no cloud tenant.

## 4. The agent-native verification layer (strongest strategic finding)

Benchling and DNAForge published **BenchBench-Cloning**
([arXiv 2610.33355](https://arxiv.org/abs/2610.33355), 23 October 2026, CC BY 4.0;
[blog](https://www.benchling.com/blog/how-good-is-ai-at-molecular-biology),
30 October 2026): 65 expert-built cloning tasks, where protocol tasks are scored by
**binary pass in a cloning simulator** and construct-design tasks by
**task-specific deterministic checks**. The best model scored **30.6%**
single-attempt; protocol pass@1 was 27% against pass@10 of 72%. Their conclusion:
"The capability is emerging, but the reliability isn't there yet."

Their enumerated failure modes are exactly what a deterministic engine catches and
a language model cannot: primers acceptable individually but failing as a pair,
annealing sites too far apart, self-compatible vector ends not dephosphorylated,
blunt or compatible ends treated as directional, ignored Dam/Dcm methylation,
missing experimental controls, over-building a factorial panel beyond what was
asked, and construction records with unresolvable sources or duplicate identifiers.

Benchling's own 2026 MCP offering is data access to the ELN and registry, not
design. SnapGene has no public API. **The design-and-verify surface is open.**

What to expose that nobody does: a structured **design-intent object** separable
from the chosen route; **deterministic replay** of a route to a byte-identical
product with every intermediate; **typed machine-checkable verdicts** carrying the
rule's citation rather than prose; **all independent blocking errors at once**, with dependent checks explicitly
blocked and resource limits retained; and
"explain why this failed" grounded in which named check fired.

Precedent for the representation, and a license-clean one: the **Construction File**
spec (Anderson lab, *PLOS ONE* 18(11):e0294469, 13 November 2023, CC BY 4.0), whose
simulator "emphasizes detecting design errors rather than only computing products".
Supporting evidence for the deterministic split: LLM judges flip verdicts at
temperature zero and under position swaps. Also instructive, a generative plasmid
design study where **16 of 1,000 generated plasmids survived filtering and 3 were
built**; its hand-rolled filters (exactly one origin at ≥99% identity, one
resistance gene at 100%, no internal repeat ≥50 bp) are precisely the deterministic
gate a generative agent needs.

Regression-suite data: LAB-Bench CloningScenarios is CC BY-SA 4.0
([arXiv 2407.10362](https://arxiv.org/abs/2407.10362)). BenchBench's graders are
stated to be at `github.com/benchling/benchbench-cloning`, but that URL returned
404 — verify before depending on it.

## 5. Long-read clone verification

What vendors return, from their own documentation:
[Plasmidsaurus](https://plasmidsaurus.com/faq/plasmid) gives one consensus
"starting at the origin", FASTA plus pLannotate-annotated GenBank, per-base
confidence (.txt/.tsv), coverage and read-length plots, a multimer analysis file,
and synthetic .ab1 chromatograms; FASTQ on request, but "reads that don't align to
the consensus aren't included". Azenta Plasmid-EZ adds **VCF plus a tab-delimited
variant file** and a per-position `base_count.tsv`. ONT's
[wf-clone-validation](https://nanoporetech.com/document/epi2me-workflows/wf-clone-validation)
emits a consensus FASTQ with per-base quality, annotations, and a pass/fail derived
from just two global numbers (95% coverage, 99% identity); its docs state
"The workflow has no way of reporting contaminants."

What users still do by hand is all of the interpretation. No vendor compares the
consensus to the **design**; Plasmidsaurus and Unveil explicitly tell users to
align it themselves. Published practice is to treat homopolymer indels as
artifacts: in *Cell Genomics* 2024, 98% of clones "showed a perfect match to the
anticipated map", and "small insertions or deletions in homopolymer runs were
treated as nanopore artifacts and considered as passing".

Error modes with numbers, for the uncertainty layer:

- Homopolymers dominate. nPoRe ([PMC10022090](https://pmc.ncbi.nlm.nih.gov/articles/PMC10022090/))
  reports a 41.8% indel error rate at homopolymer length 10; a 2026 framework finds
  accuracy "frequently falls below 50% in long homopolymer tracts (N > 10)",
  skewed to deletions.
- **Mixed-clone sensitivity has a hard floor.** On clonal plasmids with no true
  variants (*PLOS ONE*, 5 October 2026, [doi:10.1371/journal.pone.0356259](https://doi.org/10.1371/journal.pone.0356259)),
  LoFreq at 2000× averaged 8.25 false calls on R9.4.1 versus 2.00 on R10.4.1; the
  lowest reliably detected allele frequency was **10%**, 5% unreliable, and
  calibration is "unsuitable" for indels. Any mixed-clone claim below ~10% is not
  supportable.
- Structural artifacts near ITRs are real: recombination frequencies of 82.5%, 60%
  and 20% at 4, 41 and 182 bp from an intact ITR versus 0% at 543 bp
  ([PMC13566322](https://pmc.ncbi.nlm.nih.gov/articles/PMC13566322/)).

There is **no verification reporting standard**. VCF is used in practice but has no
specification or convention for circular contigs or origin-spanning variants, so an
exported VCF must pin a linearisation and say so in its header. GFA appears only as
an assembler intermediate and no vendor delivers it, so GFA support would be
speculative. Features 1–5 below need neither MSA nor BAM/BAI, so both deferrals hold.

## 6. Modality-specific checks with hard published rules

**AAV**, the best-quantified area in this research. Radukic et al., *NAR*
53(2):gkae1170, 9 December 2024,
[doi:10.1093/nar/gkae1170](https://doi.org/10.1093/nar/gkae1170): ITR deletions
arise by slipped misalignment at short direct repeats; **minimum stable ITR-to-origin
distance 438 bp**; and notably **Stbl3 does not prevent ITR deletion**, contradicting
standard vendor advice. Extended repository-wide (bioRxiv 2026.09.24.753884): of
**7,041 Addgene AAV plasmids, 4,773 (67.8%) carry an exact Δ11 motif**, and in
**95.3%** of those the affected ITR is the **origin-proximal** one; spacing above
1.4 kb makes the loss near-undetectable. GenScript's own flyer tabulates five real
ITR mutations where **the observed gel band equals the expected band** — the usual
manual check silently passes broken ITRs. Addgene concedes Sanger and NGS "might
fail to read through ITR sequences".

Computable without calibration: ITR-to-origin distance against the 438 bp and 1.4 kb
findings, origin-proximal ITR warning, exact Δ11 motif detection, ITR length
classification, 5′/3′ asymmetry, short-direct-repeat scan at junctions, and an
ITR-to-ITR budget with per-element accounting. Report the size ceiling as a graded
band with sources, not a threshold: the "4.7 kb limit" is vendor-dependent
(VectorBuilder states 4.7 and 4.2 on different pages; Addgene ~4.5) and the measured
titre drop is a **step at ~4.9 kb, not a slope**. CpG has **no published pass/fail
criterion** — count and report it, marking the ITR-resident CpGs as irreducible.

**Lentiviral**: cis-element presence and geometry (SIN U3 deletion, R-U5, Ψ, RRE,
cPPT/CTS, WPRE); measure **LTR-to-LTR, not plasmid size**; scan both strands of the
packaged genome for AATAAA/ATTAAA and splice consensus; warn on antisense internal
cassettes. Per-kb titre coefficients are folklore (Addgene's "10-fold per 2 kb"
versus others' "3-fold per kb") and there is **no published homology-length to RCL
threshold**, so a cross-plasmid homology report must be labelled geometry, not risk.

**mRNA/IVT**, unusually hard rules and untouched by incumbents:

- TriLink CleanCap AG requires the template to read `TAATACGACTCACTATAAGG` —
  **AG at +1/+2**. This **directly conflicts** with the highest-yield T7 start.
- T7 start context (Conrad et al., *Commun Biol* 3:439, 2020,
  [doi:10.1038/s42003-020-01167-x](https://doi.org/10.1038/s42003-020-01167-x),
  open access): +1..+3 GGG best, **+4..+8 motifs span a 5-fold yield range**,
  sequence-specific rather than composition-driven, full ranking in supplementary
  data; **79% of GGG-initiated transcripts carry an extra non-templated G**.
  Do not transfer T7 rules to SP6.
- Linearisation-site picker: blunt or 5′-overhang preferred, 3′ overhangs
  "generally should be avoided" (Thermo MAN0012652; NEB concurs). **This reuses our
  existing restriction catalogue and cut-geometry machinery almost for free.**
- Poly(A) recombination, measured (Trepotec et al., *RNA* 25:507, 2019,
  [PMC6426288](https://pmc.ncbi.nlm.nih.gov/articles/PMC6426288/)): continuous A120
  recombined in >50% of cases; 2×60 split by a single **G** gave **0%**; a single C
  ~50%; spacers ≥12 nt lost expression with no stability gain.
- dsRNA byproducts come from promoter-independent antisense initiation at the blunt
  DNA end; the dominant lever is Mg²⁺ 30→5 mM, so this is protocol advice and
  explicitly not a sequence fix.

**Editing beyond cut sites**, deterministic parts only: base-editor windows
(positions 4–8 for BE3/BE4max/ABE7.10, ABE8e 3–10), bystander-site enumeration with
codon consequence, CBE stop installation, PrimeDesign's published pegRNA rules
(PBS 10–15 nt, RTT 10–30, no C at extension position 1, PAM-disrupted ranked above
intact), ssODN geometry, HITI junction frame checks, Cas12a PAM tiering and the 38
validated direct-repeat variants for multiplex arrays, and library filters (5′G for
U6, Type IIS site exclusion, poly-T termination). Coverage needs a declared
sampling model. For S equally abundant variants sampled independently with
replacement, the expected observed fraction after T draws is exactly
`1 - (1 - 1/S)^T`; for S>1 and 0<p<1, the corresponding draw requirement is
`ceil(ln(1-p) / ln(1-1/S))`. The expression `-S*ln(1-p)` is a large-library
approximation, not the exact finite-library expression or the probability of
observing every variant. Nonuniform pools need their actual declared frequencies.
Software licence claims still require artifact-specific verification: **Rule Set 3 is Apache-2.0**, **PRIDICT2 is MIT
with weights**. Do **not** vendor BE-Hive (no LICENSE file); DeepPrime and
DeepBaseEditor are academic-only. Bystander genotype frequencies and pegRNA
efficiency remain ML-only, and Richardson's asymmetric ssODN design failed to
reproduce in mouse zygotes, so it belongs as a selectable template, not a default.

**Protein-design to construct compilation**, the clearest unmet need found, with a
leading lab's own admission. SAPP (Wicky et al., *Nat Commun*,
[doi:10.1038/s41467-026-76740-9](https://doi.org/10.1038/s41467-026-76740-9)) states
that "DNA fragments are typically manually assembled into plasmids" and that after
automation "at least 80% of SAPP's total cost comes from synthetic DNA". Their
in-house code is exactly the missing layer: reverse-translate, constrain for
repeats and Type IIS sites and alternative starts, **re-weight adaptively when a
sequence cannot be made**, split into orthogonal-overhang fragments, and append
vector-specific tags across 25 vectors. Baker-lab binder campaigns pad designs to
65 aa *because 195 nt fits a cheap array oligo* — DNA constraints are already
dictating protein design. OpenAntigens (*NAR*, 29 October 2026; annotations
CC BY 4.0) does the structural half and deliberately stops at FASTA, saying signal
peptides, tags and vectors "remain part of the user's experimental design".

Computable: frame-safe fusion assembly with register bookkeeping, host-specific
construct variant matrices, pool-wide k-mer uniqueness, and structural annotations
(pLDDT, burial, N-X-S/T sequons, free cysteines, domain edges) carried onto DNA
coordinates. **UniProt is CC BY 4.0** and is the right source for signal peptides,
domains and disulfides; **SignalP 6.0 is fee-licensed and cannot be bundled**.
Authors' own caveats: pLDDT/PAE boundary cutoffs are uncalibrated
("reproducible hypotheses, not expression-success predictions"), and 2A versus IRES
is locus-dependent, so both belong as annotated options rather than recommendations.

**Library and pool design** is self-contained with no licensing friction: the
22-c-trick NDT+VHG+TGG mixture at 12:9:1 (Kille et al., *ACS Synth Biol* 2013,
[doi:10.1021/sb300037w](https://doi.org/10.1021/sb300037w)); SwiftLib's dynamic
program, fully specified in its paper (*NAR* 43:e34, 2015) and therefore
reimplementable from the text, since the reference repository has no LICENSE;
sequence-Levenshtein barcodes (*BMC Bioinformatics* 14:272, 2013) and FREE
divergence (*PNAS* 115:E6217, 2018) with the ≥2m+1 bound; sequencer colour-balance
checking; and nicking mutagenesis's single-site requirement. Cannot ship: the RBS
Library Calculator (proprietary) and DeCoDe (requires Gurobi). Anderson promoter
strengths must ship as a *convention with provenance* — independent iGEM
remeasurements disagree substantially.

**Control-construct proposal, not a DoE recommender.** The published failure mode
cuts against eagerness: BenchBench's top two models **over-built** factorial panels
beyond what was asked while simultaneously **omitting required controls**. The
defensible capability is narrow and deterministic — propose the control constructs
a stated design intent implies (closed-frame, non-targeting, safe-harbour, empty
vector) and enumerate a factorial panel only on explicit request.

## 7. Standards and toolkit grammars

- **Assembly-standard grammars as data.** `Edinburgh-Genome-Foundry/genedom` (MIT)
  carries `EMMA.csv` and `YTK.csv` with the schema
  `slot_name,left_overhang,right_overhang,enzyme,left_addition,right_addition,extra_avoided_sites,description`.
  Adopt that schema rather than inventing one. `althonos/moclo` (MIT) adds MoClo,
  YTK/PTK, CIDAR, EcoFlex and plant kits as overhang signature tuples.
  **Source the overhang tables from the MIT code repositories, not the literature**:
  Patron 2015 (Common Syntax), Pollak 2019 (Loop) and Iverson 2016 (CIDAR) are all
  paywalled and not in PMC.
- **iGEM.** The 2026 Judge Handbook (p.18) states that submitted sequences "are
  automatically screened for assembly-method compatibility" and that parts "must be
  compatible with either BioBrick Standard [RFC 10] or iGEM Type IIS Standard
  [RFC 1000]", enforced in the judging form. Teams find out after submitting. A
  pre-flight panel plus a documentation-completeness checklist mirroring the rubric
  is the highest leverage per line of code in this research. **Unverified:** the
  exact RFC 1000 fusion sites and RFC 10 prefix/suffix and forbidden-site list —
  iGEM blocks automated fetches, so this needs a human browser session.
- **SEVA** is a grammar we can lint and generate: the four-position cipher, module
  flanks (`FseI`–`AscI`, `SwaI`–`PshAI`, `PacI`–`SpeI`) and the 20 sites stripped
  from all SEVA backbones. SEVA 4.0 is CC BY; the site warns some sequences may
  carry IP restrictions and offers no bulk download.
- **SBOL beyond structure**, three classes worth adding to our existing structural
  profile: **CombinatorialDerivation + VariableFeature** maps one-to-one onto a
  Golden Gate library, with `strategy`, `cardinality` and `Constraint` objects
  ruling out illegal combinations — exporting one derivation instead of N enumerated
  plasmids is the single best SBOL win and directly answers the practitioner ask in
  §2. **Implementation** links a design to a physical realisation; **PROV-O
  Activity/Plan/Agent** serialises history we already hold. Skip Model, Interaction,
  ExperimentalData and OM units: no consumer in our product.
  **Licence flag:** sbolstandard.org content is CC BY-NC-ND 4.0, so read the spec
  freely but copy no prose or tables into the product. SBOL Visual glyph licensing
  is **unresolved** and blocking for glyph rendering.
- **Addgene**: a deposit-package export is light and useful, but do **not** mirror
  their catalogue — API access requires an accepted per-scope data licence and their
  terms restrict content to informational, non-commercial use.

## 8. Predictive capabilities, with their limits

Recorded for completeness. Each would be a new predictive model and must arrive
with method, parameter licence and validation, as required of CRISPR scoring data.

- **Golden Gate set-fidelity optimisation.** The Pryor 2020 matrices are **already
  bundled** (see [published ligation profiles](../PUBLISHED_LIGATION_PROFILES.md)),
  so the remaining opportunity is the **optimizer** — repositioning junctions to
  maximise predicted set fidelity, choosing 3 bp versus 4 bp overhangs, and
  distributing a fidelity budget across hierarchical levels. Validation is
  published as predicted-versus-observed (13-fragment SapI 79% predicted / 91%
  observed; 35-fragment BsmBI-v2 65% / 71%) and the authors call their predictions
  "qualitative". Reimplement the arithmetic; NEB's own code is AGPL-3.0.
- **Genetic stability linting** (EFM-style repeat-mediated deletion and SSR hotspots).
  Implementable with a suffix automaton; the risk score is relative to an *E. coli*
  baseline from a modest measured set, so present ranked risk, never a failure
  probability.
- **Synthesis manufacturability.** 303 commercially ordered fragments, ~22% failed;
  three thresholds carry most of the signal (longest repeat >26 bp, ≥16 overlapping
  9-mers per 100 bp, GC outside 29–63%). Twist publishes a full numeric grid; ship
  vendor grids as **versioned profiles with provenance** and surface their live
  contradictions (Twist 5 kb versus 7 kb; IDT eBlocks 900 versus 1500 bp) rather
  than silently resolving them. Footgun worth catching: **Twist oligo pools accept
  only `N`**, not R/Y/S/W/K/M/B/D/H/V, so NNK and NDT schemes cannot be ordered as
  pools.
- **Secondary structure** is the enabling dependency for several of the above. A
  clean-room Zuker/McCaskill over Turner-2004 NNDB parameters is the safe route;
  **ViennaRNA permits commercial use but bars redistribution for a fee**, and the
  parameter files' own licence still needs verifying. **DegScore is MIT** and is
  ridge regression over sequence and loop type, trivially portable.
  **LinearDesign is not redistributable.**
- **Weak or uncalibrated, do not ship as scores:** Gibson/HiFi overlap rules
  (vendor heuristics, no peer-reviewed dose-response study; overlap *non-uniqueness*
  is the one well-grounded warning); RBS Calculator (parameters not fully open,
  code academic-use-restricted, roughly within 2-fold for about half of its
  validation set); terminator strength (~40% of variance, R² ≈ 0.29 on
  re-evaluation — a catalogue of *measured* terminators is more honest and more
  useful than the model); internal Shine-Dalgarno pausing (contested as a footprint
  artefact); %MinMax (its own authors state it does not predict local translation
  rate); any slippery-site frameshift score (**no supporting validation found**);
  pLDDT/PAE construct boundaries (the authors call them uncalibrated); protein
  liability motifs (~60% false-positive without burial context).
- **No public labelled dataset exists** of Gibson/HiFi or Golden Gate *reaction*
  success keyed to junction sequence. Any assembly-success probability beyond the
  ligation-fidelity model would be uncalibrated and must be presented as a feature
  report instead.

## 9. Research limits

Reddit, SEQanswers, Biostars, Bluesky, protocols.io comments and the
Benchling/SnapGene feature-request boards were unreachable, so practitioner venting
about price, speed and spreadsheet rituals is under-represented. Vendor help centres
block automated fetches. Two agents exhausted their search budgets. iGEM RFC
details, SBOL Visual glyph licensing, Addgene per-plasmid deposit fields and several
parameter-table licences remain unverified and are flagged above.

Biosecurity screening was **excluded by user decision on 28 November 2026** and was
not pursued; one agent independently reported that its related policy sub-research
was halted by model safeguards, which is consistent with that decision.


## Implementation review — 28 November 2026

The active plan is the implementation authority; this research is not a source of
automatic scientific or licensing authority. Parallel source work is assigned to
F37 assembly-standard data/grammar and F38 oligo-binding crates, plus a core-owned
F36 linter that reuses translation/annotation algorithms without a dependency
cycle. No worker runs its own compilation or application checks.

Corrections and admission limits:

- The current BenchBench-Cloning primary text describes deterministic graders
  plus rubric-based LLM judging on ten construct-design tasks. Its aggregate score
  is not an entirely deterministic benchmark. Collagene's verdicts remain native
  deterministic results; no LLM judge is adopted.
  [Primary text, section 3.2](https://arxiv.org/html/2610.33355v1).
- A search that found no vendor documentation cannot establish that no vendor
  offers a feature. The broad negative competitor statements above are bounded
  research findings, not verified universal product comparisons.
- Primary-source review separates design errors (384/2,521), restriction-pattern
  inconsistencies (128/852) and directly sequenced variations (91/259). These are
  distinct cohorts. Taxonomy categories are candidate checks, not proof every
  error is decidable from a sequence or its annotations.
  [Primary paper](https://academic.oup.com/nar/article/53/14/gkaf697/8222435).
- Existing barcode scope is unresolved for the new genetic barcode-set proposal;
  no barcode implementation starts until the user clarifies that earlier exclusion.
- Design-versus-consensus work can use current sequence imports and alignment;
  this does not reopen deferred BAM/BAI or raw long-read processing.
- Software licences for Rule Set 3/PRIDICT2 do not automatically establish licences
  for every parameter table, model weight or evaluation dataset. Exact admitted
  artifacts, scientific applicability and independent validation are still required.
- Public overhang/standard facts must be independently curated with data-specific
  provenance. No implementation is copied from another cloning program. A bounded
  RFC10 part-sequence profile is not full vector/strain/registry compliance.
