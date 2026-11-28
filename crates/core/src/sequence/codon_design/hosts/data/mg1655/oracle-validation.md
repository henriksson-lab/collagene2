# Independent count oracle

The captured MG1655 source SHA256s (and upstream MD5s where supplied) were checked
before decompression. A Python standard-library maintainer oracle then built direct
GFF `(locus_tag, protein_id)` exclusion sets and separately tallied each admitted
FASTA CDS's triplets. It did not invoke the Rust parser, translator, profile code,
third-party biological algorithms or external tools. This is a curation check, not
application runtime functionality.

The independent path asserts one ordinary single-range GFF CDS per admitted
locus/product; exact FASTA coordinate/strand/length agreement; phase zero; no
partial/pseudo/translation-exception flags; concrete frame-complete DNA; explicit
code-11 start membership and ordinary TAA/TAG/TGA termination; no internal stops.
Initiation and terminal codons are counted separately. Exceptional multipart
products are excluded by their direct GFF identities, rather than replaying Rust's
general location parser and skip-priority logic. Distinct paralog loci remain
separate even with an identical protein identifier or sequence. No independent
whole-genome base extraction or expression inference is claimed.

All **192** entries matched the freshly derived Rust output: 64 internal codon
counts, 64 initiation counts and 64 terminal-stop counts. Totals were 4,318 observed
locus/product records, 4,294 admitted, 18 pseudogene exclusions, six translation
exceptions and 1,322,130 internal codons.

Observed starts: ATG 3,870; GTG 338; TTG 80; ATT 4; CTG 2. Observed stops: TAA
2,749; TAG 305; TGA 1,240. `expected-counts.json` stores all arrays in explicit
A/C/G/T lexical triplet order, raw source SHA256s and the compared Rust profile
identity. Its counts are independently tallied observations, not expression data.

The independently executable maintainer oracle is retained as
`independent_oracle.py` in this directory. Run `python3 independent_oracle.py`
from here, or pass `--data` and `--profile` explicitly. It verifies all source pins
before tallying and does not modify files by default. Optional `--output` writes
the detailed result; `--write-fixture` deliberately regenerates a numeric fixture.
The application and its scientific runtime never invoke this script.

Rust profile identity compared:
`52d52265b4f997e5f95c2d4178bcfc575e5002bef278a7a327a9efe6f8d61f43`.
