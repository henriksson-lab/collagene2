# Embedded NCBI genetic-code catalogue

The Rust translator embeds the 27 extant tables in NCBI `gc.prt` version 4.6:
**1–6, 9–16, 21–33**. Tables absent from this finite catalogue, including removed
or reserved IDs, are rejected. No runtime download or external translator is used.

The complete private `Code` owner, lookup, descriptor projection and table data
live in `core/src/sequence/translation/codes.rs`; the decoder, report types and
meaningful implementations remain in `translation/mod.rs`. The canonical public
module path is unchanged. There is no forwarding facade or alternate decoder.
`CodeDescriptor::for_id` replaces duplicated admission allowlists in reference
translation, codon design and assembly junction inspection. Existing web dropdowns
already derive their entries from this same catalogue.

## Primary data and provenance

- [Machine-readable NCBI gc.prt](https://ftp.ncbi.nih.gov/entrez/misc/data/gc.prt),
  version 4.6, captured unmodified in `core/tests/data/genetic_codes/gc.prt`.
- [NCBI genetic-code descriptions](https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi),
  displaying revision 23 September 2024; its explanations distinguish dual-use
  sense/termination codons from ordinary definite stops.
- `core/tests/data/genetic_codes/provenance.json` retains the URL, retrieval time,
  exact byte count, SHA-256 and supported IDs. The captured file is 16,444 bytes;
  SHA-256 is `2aecbdd09ecc5e0d233c105090a593737298ef93d88f6bc0c1bef26df6051a04`.

The rows are scientific table data, not code from another cloning application.
TCAG ordering is explicit. The `ncbieaa` row gives internal residue assignments;
`sncbieaa` marks accepted initiators with `M` and possible termination with `*`.
An asterisk in `sncbieaa` with a sense residue in `ncbieaa` therefore requires
context rather than replacement with a definite stop.

## Context-dependent termination

| Table | Dual-use codons | Report behavior without organism/position context |
| --- | --- | --- |
| 27, Karyorelict nuclear | TGA: tryptophan or termination | `X`, possible stop |
| 28, Condylostoma nuclear | TAA/TAG: glutamine or termination; TGA: tryptophan or termination | `X`, possible stop |
| 31, Blastocrithidia nuclear | TAA/TAG: glutamate or termination | `X`, possible stop |

DNA/RNA and IUPAC expansions use the same decoder. A concrete dual-use codon
increments `possible_stop_codons` and `unresolved_codons`, but does not become an
ambiguous nucleotide codon. A synonymous IUPAC codon that includes a dual-use
expansion remains unresolved even if all sense assignments agree. Definite
stop codons remain `*`; this does not infer programmed recoding or exceptions.

A terminal dual-use codon never supplies `has_terminal_stop` or CDS completeness.
ORFs stop at it as an explicit `AmbiguousStop` boundary; complete-only searches
exclude that unproven candidate. Partial candidates retain the existing partial
flags and termination reason. The same rule applies on circular sources.
Codon design excludes uncertain codons from canonical-residue domains and its
independent checker rejects candidates translating through them. Required coding
DNA with unresolved residues is refused. A protein requiring a terminal `*`
has no solution in a table that contains no definite stop.

## Validation scope

Fixtures contain literal reassigned-codon expectations for every added table,
RNA/reverse-oriented equivalents, explicit starts, reserved-ID rejection,
IUPAC/context-stop cases, CDS and linear/circular ORF uncertainty, unsafe codon
weighting/checker refusal, reference output and current binary boundary cases,
and junction crossing-codon evidence. A separate captured-row correspondence
case checks all 64 codons and starts for each table against the retained primary
record. That correspondence is distinct from the independently written literal
scientific expectations.

Tables 1, 2 and 11 retain their previous assignments, initiation sets, names and
output behavior. No existing worker payload layout changes for this catalogue;
`CodeDescriptor` adds explanatory stop-context metadata. No cloning-verifier
chemistry algorithm changes. The parent integration batch owns current worker
schema changes, final compiler/test execution and cross-platform CI evidence.
