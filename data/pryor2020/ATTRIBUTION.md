# Pryor 2020 ligation-frequency data

Pryor JM, Potapov V, Kucera RB, Bilotti K, Cantor EJ, Lohman GJS (2020).
**Enabling one-pot Golden Gate assemblies of unprecedented complexity using
data-optimized assembly design.** PLOS ONE 15(9): e0238592.
https://doi.org/10.1371/journal.pone.0238592

Original CC-BY-4.0 supporting information, downloaded 2026-11-25 directly from PLOS:

- S1, BsaI-HFv2: https://doi.org/10.1371/journal.pone.0238592.s001
- S2, BsmBI-v2: https://doi.org/10.1371/journal.pone.0238592.s002
- S3, Esp3I: https://doi.org/10.1371/journal.pone.0238592.s003
- S4, BbsI-HF: https://doi.org/10.1371/journal.pone.0238592.s004
- S5, SapI: https://doi.org/10.1371/journal.pone.0238592.s005

`manifest.json` records the original file byte sizes, SHA-256 hashes and publisher
URLs. Runtime metadata additionally computes the actual original byte BLAKE3 and
the canonical converted-data revision. No runtime network request is made.

## Independent conversion and interpretation

Original XLSX bytes are unchanged. The bounded Rust reader maps literal column
and row labels to canonical 5′→3′ end strings. Columns are lexically ordered;
rows are reverse-complement ordered. The spreadsheet's visual diagonal therefore
represents Watson–Crick pairs, not identically named self-pairs. All explicit
integer entries, including zeros, are preserved. Complete label coverage, missing
cells, duplicate labels/cells, symmetry, formulas, external parts and resource
bounds are independently checked. There is no transposition fill, averaging,
row normalization, rescaling or pseudocount.

The matrices are symmetric **published ligation-frequency weights**, not a new
claim of independent raw directional read observations. Replicate pooling and
normalization before publication are not specified in the workbooks. This
interpretation is versioned separately from the original caller-supplied directed
observation model. The same independently written conditional numeric kernel
pools the two physical-end rows and multiplies once per unique RC pair. Its
qualitative prediction does not establish reaction efficiency or yield.

The exact assay formulation, buffer, cycling and substrate conditions accompany
every profile. Shared recognition/cut geometry never silently substitutes one
enzyme formulation for another. In particular Esp3I is not labelled BsmBI-v2,
and these 2020 cycling datasets are not the static 2018 dataset described by
SnapGene. Palindromes retain the scorer's explicit Unknown policy.

No NEB optimizer implementation was accessed or borrowed; the importer and
scoring code are independently written and MIT licensed. Data remains CC BY 4.0.
