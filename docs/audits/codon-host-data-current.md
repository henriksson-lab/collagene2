# Host codon data: concrete next F26 input

Research checkpoint, 2026-11-28. Recommendation: derive small, versioned **observed
reference-CDS count profiles**, starting with E. coli K-12 MG1655, then BL21(DE3)
and S. cerevisiae S288C nuclear CDS. Bundle derived counts and provenance, not
copied vendor/Kazusa optimization tables. This supplies useful preferences without
claiming calibrated expression, strain fitness or a tissue/culture-specific model.

## Existing integration surface

`sequence/codon_design.rs` already accepts up to 64 explicit `CodonWeight` values
in `DesignOptions`, each 0–1,000,000; omitted codons score 1, and 0 forbids a codon.
Its search/checker maximize an additive score and bind complete options to their
revision. `translation::CodeDescriptor` supplies the current genetic table and
recognized starts. Thus data derivation can be a small independent Rust owner;
core needs a direct profile-selection/provenance extension, not a replacement
optimizer or external executable. Update the current options/report representation
and callers together if profile identity becomes part of retained evidence.

## Sources, provenance and reuse

NCBI documents CDS FASTA as nucleotide sequences for annotated genomic CDS, with
matching annotation/report files and checksums. It explicitly warns that annotation
can change for the same assembly version: accession.version alone does not pin a
profile. Capture the exact CDS and annotation/report URLs, retrieval date, annotation
release/date, original checksums, our BLAKE3 hashes, and exclusion-policy/engine
version. Obtain only CDS plus needed metadata; do not fetch whole human genomes.
[NCBI FTP data definition and versioning](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/data-processing/policies-annotation/genomeftp/)

NCBI itself imposes no molecular-data use/distribution restrictions, but does not
transfer submitter rights or guarantee absence of third-party claims. Consequently
“everything on NCBI is public domain/MIT” would be incorrect. Keep NCBI and original
assembly/annotation acknowledgments and the molecular-data disclaimer with derived
profiles; do not relicense upstream records as MIT. The MIT application can carry
our independently derived numeric counts with this source notice. Do not copy
article text or third-party codon tables under an assumed license.
[NCBI molecular-data and copyright policy](https://www.ncbi.nlm.nih.gov/home/about/policies/)

| Initial candidate | Assembly candidate | Scope |
| --- | --- | --- |
| E. coli K-12 MG1655 | GCF_000005845.2, ASM584v2 | Chromosome CDS; code 11; first implementation |
| E. coli BL21(DE3) | GCF_000022665.1, ASM2266v1 | Separate named strain profile; code 11 |
| S. cerevisiae S288C | GCF_000146045.2, R64 | Nuclear CDS only; code 1; exclude mitochondrial code-3 CDS |

Code selection and organelle separation follow the [NCBI genetic-code catalogue](https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi?chapter=tgencodes).
These are acquisition candidates, not already captured/verified datasets. MG1655
reference identity is confirmed by its [NCBI nucleotide record](https://www.ncbi.nlm.nih.gov/nuccore/NC_000913.3).
BL21 assembly identity is supported by an [NCBI depositor record](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM4097390).
S288C identity and current propagated annotation are visible in
[NCBI Genome Data Viewer](https://www.ncbi.nlm.nih.gov/gdv/browser/genome/?id=GCF_000146045.2).
Before bundling, verify assembly metadata, source URLs/checksums and original
annotation attribution, particularly SGD propagation. The research browsing tool could not acquire FTP/API payloads. The integrator
subsequently captured MG1655 through bounded direct HTTPS, checked the NCBI MD5
manifest, and recorded SHA256 hashes for compressed and expanded CDS/GFF bytes.
The pinned assembly report names taxid511145 and chromosome NC_000913.3
(4,641,652 bases). Total CDS/GFF/annotation-hash acquisition was under2MiB.
No derived codon table is claimed validated yet.

## Bounded pure-Rust derivation

Implement a streaming CDS record counter, returning `[u64;64]` in declared A/C/G/T
lexical order plus separate start/terminal-stop counts and skip-reason totals.
Bound compressed input to 20 MiB, expanded input to 32 MiB, records to 20,000,
individual CDS to 300,000 bases and headers to 4 KiB; cap work by expanded bytes.
Reject limit exhaustion, never publish partial aggregate counts as a full profile.
Pure-Rust gzip/FASTA handling and the embedded translator work natively and in WASM;
acquisition belongs to a maintainer HTTP step, not a runtime process dependency.

Use annotation metadata to admit only complete, unambiguous, in-frame CDS with the
explicit expected code; exclude pseudogenes, partial boundaries, translation
exceptions/slippage/recoding and organelles. Validate recognized initiation and a
single terminal stop. Exclude initiation and terminal stop from internal sense-codon
counts. Count each admitted annotated locus once, not each unique protein sequence;
do not silently deduplicate genuine paralogs. Record all exclusions. Additional
splicing/isoform policies are unnecessary for the first bacterial profile.

Expose observed counts, per-thousand usage, and within-amino-acid relative usage.
For the existing additive preference API, use a separately named deterministic
mapping such as `1 + floor(999999 * count / max_synonymous_count)`; zero observations
remain allowed with score 1, rather than accidentally becoming forbidden codons.
Missing synonymous-family data must be reported; use uniform preferences there.
Keep terminal-stop weights uniform unless an explicit separate policy is selected.
This arithmetic preference is **not CAI**, an expression prediction, or an empirical
success probability. Preserve raw counts so later scoring need not reverse rounded
weights. A source hash, derivation policy and selected profile hash must accompany
the expanded weights in retained reports.

Next work is bounded derivation plus one captured bacterial dataset and UI/MCP
selection. Human/CHO/tissue panels and empirical expression calibration can follow;
they are not prerequisites for this useful feature. No source changes, builds or
tests were performed for this research checkpoint.

## Captured first source

The exact acquisition record is [mg1655-source-acquisition.json](mg1655-source-acquisition.json).
CDS expanded bytes:5,043,083; GFF expanded bytes:2,314,238. Original metadata,
compressed sources and the upstream checksum manifest are retained in the
curation scratch directory pending the separate bounded derivation batch.
No empirical expression or fitness calibration follows from this acquisition.
