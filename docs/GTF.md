# GTF annotation exchange

`AnnotationFormat::Gtf` extends the existing reference-bound annotation import,
review/application, normalized export and exact-original-source download APIs.
It is appended after GFF3 and BED; existing DTO fields and source-evidence JSON
layouts are unchanged. Import report schema remains 1. Worker transport gates
new GTF enum payloads separately.

The supported row profile follows the nine tab-separated columns documented by
[Ensembl](https://static.ensembl.org/info/website/upload/gff.html) and the column
and repeated attribute conventions described by
[GENCODE](https://www.gencodegenes.org/pages/data_format.html). These sources
specify inclusive coordinates beginning at one and quoted attribute values;
GENCODE documents CDS phase values 0, 1, 2 or missing, plus gene/transcript IDs.
This implementation is independent MIT code.

## Import and review

The caller explicitly maps a file sequence name to an actual DNA/RNA document.
Matching names are not evidence of reference identity. Coordinates become native
half-open intervals: source `3..6` becomes `[2,6)`. Every row remains one separate
reviewed feature, including single-base intervals. Gene and transcript rows,
exons and separate CDS segments are **not** automatically merged or validated as
a transcript model. Off-reference rows are syntax-checked, retained in the exact
source and listed as skipped sequence names.

Forward/reverse strands preserve biological direction. A missing strand keeps
reference coordinates in an `Order` location, preventing exact biological
extraction. A numeric CDS phase on a known strand becomes `codon_start = phase +
1` for that individual row, using its biological beginning; reverse-strand phase
is measured from the high-coordinate end. Missing CDS phase supplies no inferred
frame or translation. Non-CDS frame values remain explicit source evidence.
Phases across rows, splice sites, coding completeness and gene ontology are not
independently verified.

Attributes use `tag "value";` syntax. Repeated tags, duplicate values, their order,
spaces, semicolons and UTF-8 text within quoted values are preserved. Explicit
escaped quotes and backslashes are supported; other backslash escape dialects
are rejected. Finite unquoted numeric attributes are accepted as a compatibility
extension and normalized to quoted strings on export. GFF3 `tag=value` syntax,
unterminated quotes/semicolons, control characters, nonfinite scores, invalid
strands/phases and out-of-bounds mapped intervals fail atomically. Missing or
repeated relationship IDs produce review notices rather than invented IDs.

Comments and track/browser directives remain inert original-source evidence.
Embedded FASTA and GFF3-version declarations are rejected under this GTF profile.
Circular origin wrapping is not inferred. Unrecognized feature types and
attributes remain uninterpreted source data.

## Export and fidelity

Unchanged imported GTF features export their original row geometry, source,
score, strand, phase and ordered/repeated attributes. Quoting, whitespace and
unquoted numeric spelling can be normalized. Exact original bytes, including
comments and CRLF line endings, remain separately downloadable; normalized
export makes no byte-preservation claim. Retained row evidence is reconstructed
and compared with current feature geometry and qualifiers before use. Changed
source-bound features require review instead of replaying stale phase or IDs.

Generic current features require one explicit `gtf:attribute:gene_id`, plus one
`gtf:attribute:transcript_id` except for a gene row. Synthetic gene/transcript
relationships are not generated. Exact nonoverlapping joined intervals export
as individual rows with those supplied IDs. Generic CDS export requires one
explicit `codon_start`; subsequent row phases are derived in biological order
from actual interval lengths. Unknown/fuzzy/mixed-strand/origin-crossing
locations are refused. Labels are retained as `cg2_label`; other unrepresented
qualifiers remain native evidence and generate an export loss notice.

GTF-to-GFF3/BED conversion uses existing conservative exporters and their loss
notices; it does not infer GFF3 parent graphs or BED coding extents. Native
workspaces remain necessary for complete evidence, review flags and provenance.

## Resource profile

Existing annotation bounds apply: source at most 8 MiB, physical lines at most
64 KiB, at most 10,000 data rows and 20,000 directives, output/report at most
32 MiB. Each row has at most 1,000 attributes, a tag at most 128 bytes and a value
at most 8,192 bytes. The parser scans quoted values once and does no network,
external software invocation, transcript graph search or scientific inference.
Normalized output is reparsed before publication. Focused literal fixtures cover
reverse CDS phase, multiplicity, independent rows, unknown strand, missing phase,
strict invalid syntax, stale evidence and explicit joined-CDS export.

A small [synthetic example](../fixtures/annotations-demo.md) pairs GTF rows with the existing seven-base GenBank reference.
