# PHYLIP, MSF and NEXUS alignment exchange

This independently written MIT Rust module exchanges existing alignments. It
adds no alignment engine, subprocess, network, filesystem access or dependency.
The same implementation is intended for Windows, macOS, Linux and WASM. PIR, GDE
and Selex remain a separate subsequent batch; BAM/BAI and OpenGel are unaffected.

## API and native preservation

`alignment_exchange::import_alignment(bytes, format, alphabet, name)` returns the
existing `alignment_formats::AlignmentDocument`. Its stable row identities derive
from the exact source hash and source row index. Native `.cg2aln` save/load,
original-byte download, edited-source detection, ungapped extraction and alignment
column coordinates use the existing shared validation/provenance path. No new
native wrapper or second scientific document model is introduced.

Six variants are appended after the historical three AlignmentFormat variants:
`PhylipSequential`, `PhylipInterleaved`, `PhylipRelaxedSequential`,
`PhylipRelaxedInterleaved`, `Msf`, `Nexus`. Old binary enum tags keep their meaning.
The main `alignment_formats` import/export functions dispatch these profiles too.
Unsupported older enum decoders reject the new variants. Worker/native transport
version floors are an integration responsibility, separate from source format
parsing.

`export_alignment(document, format)` uses conservative defaults.
`export_alignment_with_options(document, format, &ExportOptions)` accepts
`truncate_phylip_identifiers` (default false) and `block_width` (default 80,
allowed 1..1000). Both return the existing AlignmentExport: text, ordered explicit
loss messages and an identifier renamings map. Canonical output is different from
exact evidence export. Original bytes, including uninterpreted blocks, survive
native save/load even after current rows are edited.

The caller always supplies DNA, RNA or protein. An alphabet-looking sequence
never selects the alphabet. Native identifiers now preserve bounded printable
ASCII labels, including internal spaces. Legacy FASTA, CLUSTAL and Stockholm
exports explicitly refuse whitespace-bearing labels because their unquoted name
syntax cannot preserve them. NEXUS quotes names; strict PHYLIP uses its fixed name
field. No import silently replaces spaces or truncates identifiers.

## PHYLIP profiles

The official PHYLIP documentation defines a two-number dimension header,
fixed ten-character names, sequential versus interleaved layouts and ignored
sequence whitespace/numeric counters. Only first-block rows carry names in the
interleaved profile. This implementation follows those layouts and requires
identical ordered block widths. Strict fields preserve internal name spaces;
relaxed profiles use unbounded-by-ten whitespace-delimited identifiers. The mode
is explicit, without ambiguous sequential/interleaved or strict/relaxed guessing.
[Official PHYLIP sequence input documentation](https://phylipweb.github.io/phylip/doc/sequence.html).

This profile accepts native-supported uppercase/lowercase molecular symbols and
hyphen gaps. It rejects dot/missing-data/deletion extensions that would require
semantic conversion; DNA/RNA X, O and '?' are not silently changed into N or gaps.
PHYLIP protein J/O/U extensions are outside its documented legacy alphabet and
are rejected on import/export. Protein stop '*' remains a residue symbol.
Options in headers, multiple concatenated data sets, repeated continuation names
and additional rows are rejected rather than guessed.

Strict export rejects names exceeding ten bytes by default. Explicit truncation
records every changed identifier. Any resulting collision or ambiguous trailing
padding is an error, even with truncation enabled. Relaxed export and MSF reject
name whitespace. Canonical output never inserts a tree or computes confidence.

## GCG MSF profile

GCG's manual specifies the AA/NA type header, alignment dimensions/checksum line,
individual name/length/checksum/weight records, separating `//`, and aligned
sequence rows. The checksums describe file contents, not biological validity.
[Original GCG Wisconsin Package user guide](https://bip.weizmann.ac.il/education/materials/gcg/using_sequences.html).
[EMBOSS official MSF format example](https://obf.io/emboss.open-bio.org/html/use/apas03.html).

This implementation requires the `!!AA_MULTIPLE_ALIGNMENT 1.0` or
`!!NA_MULTIPLE_ALIGNMENT 1.0` header, P/N type agreement with the explicit alphabet,
equal declared aligned lengths and complete ordered blocks. Numeric ruler lines
are ignored. Both `Check:` and EMBOSS `CompCheck:` aggregate headings are accepted.
Finite nonnegative weights are retained as GS/WT annotations and preserved on MSF
export. Other preamble content is retained in original bytes and summarized as
metadata; source text is not executed.

Nonzero checksums are verified over the original uppercase sequence characters,
including their actual gap glyphs, using weights cycling 1..57 and modulo 10000.
The aggregate is the sum of row checksums modulo 10000. Zero fields are accepted
as unspecified and explicitly reported, never called verified. Dot and tilde gap
cells become canonical hyphens without changing column coordinates. Export uses
dot gaps, regenerates checksums, preserves supplied WT values and defaults absent
weights to 1.00. MSF's nucleotide N type does not distinguish DNA from RNA; the
explicit caller alphabet is still required. Profiles with differing per-row Len
interpretations or reordered/incomplete sequence blocks fail explicitly.

## NEXUS molecular matrix profile

The original NEXUS paper describes modular blocks and molecular character
matrices. A file can contain information that a particular reader does not
interpret. That extensibility is preserved as original evidence, not a claim that
this module implements phylogenetics or understands arbitrary private commands.
[Maddison, Swofford and Maddison, 1997, Systematic Biology 46:590–621](https://doi.org/10.1093/sysbio/46.4.590).
[Primary paper PDF](https://www.phylo.org/sub_sections/portal2_portal2/nexus_document.pdf).

Supported: one DATA or CHARACTERS matrix; optional single TAXA/TAXLABELS block;
explicit DNA/RNA/PROTEIN datatype; checked NTAX/NCHAR; sequential and interleaved
rows; repeated labels or consistently unlabelled continuation blocks; configurable
gap symbol and explicitly declared MATCHCHAR; nested comments and quoted labels
with doubled apostrophes. Unquoted underscore label syntax becomes a space as
specified; a quoted underscore remains literal. Quoted export preserves literal
names and punctuation. MATCHCHAR copies only the actual prior first-row cell at
the same declared column; it never uses another desired sequence.

Whole unknown blocks such as TREES, SETS and private blocks are retained exactly
in original bytes with explicit not-interpreted notes. Canonical export emits one
molecular DATA matrix and reports the omitted semantics. TITLE is documentary;
LINK and unknown commands in interpreted blocks are rejected. Multiple character
matrices/taxa blocks, TRANSPOSE, mixed datatypes, EQUATE, TOKENS, RESPECTCASE,
polymorphic/set-valued cells and other FORMAT semantics are unsupported errors.
Actual missing-data cells are also rejected: '?' is not silently asserted to be
an existing residue or a gap in the current native cell model. Declaring MISSING
without using it is accepted. NEXUS trees are neither inferred nor executed.

## Bounds and independent fixtures

Shared bounds remain 8 MiB source/output, 1000 rows, 8,000,000 aligned cells,
65,536 bytes per line and 20,000 annotations. Dimensions are checked with checked
arithmetic before sized allocation. Every row append checks remaining declared
width. NEXUS tokenization additionally caps 300,000 tokens and comment nesting at
128. Annotation appends are capped during parsing. Output lines and total bytes
are checked while writing, so an excessively repeated long label cannot allocate
an unbounded canonical file. MSF weights are indexed once rather than scanned
for every output row. Errors return no truncated or partially accepted document.

The `alignment_exchange` test target contains eight independent fixture groups:
known gapped rows/coordinates and source-preserving native exchange; explicit
PHYLIP truncation collisions; malformed headers/block geometry/alphabet mismatch;
hand-computed MSF checksums and tampering; quoted NEXUS labels, literal underscore,
MATCHCHAR and opaque blocks; RNA/protein and unlabeled interleaving; unsupported
matrix semantics; input/output option bounds. No real aligner or cloning software
source code is used to define expected sequences. Tests must be run by the owning
integration batch; this source change did not run compiler or browser checks.
