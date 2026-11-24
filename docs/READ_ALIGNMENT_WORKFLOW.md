# SAM/BAM source inspection and reference checks

Collagene retains SAM/BAM as read-alignment source collections, independently of
editable molecule documents and reference-free read assemblies. An imported
mapping, CIGAR, MAPQ, tag, header name, M5 checksum or UR address is not proof
that the mapping is correct. No header URI is fetched implicitly.

The parser and writer are independently implemented in the MIT-licensed
`collagene-read-alignments` crate. The core collection owns the exact original
bytes and privately completed parsed records. It publishes a collection only
after the stream completes, including BAM/BGZF framing and integrity checks.
JSON summaries and worker receipts cannot recreate that collection.

In the browser, open a `.sam` or `.bam` file through **Open files**. A SAM header
also identifies a SAM source when the filename lacks its usual extension.
Choose a header reference and an actual DNA document explicitly before checking
a selected record. The checker binds the original record/source hashes and the
actual molecule revision, interprets CIGAR against that DNA and reports original
query/reference coordinates. Missing calls, unresolved symbols, unsupported
states and bounded output remain explicit. Quality and MAPQ are source claims;
the viewer does not infer calibrated confidence, coverage, variants or consensus.

The browser profile accepts up to 16 MiB of source bytes, charges at most 64 MiB
of retained collection buffers and limits collections to 100,000 records. Pages
contain at most 100 records; display previews are bounded independently. A source
that exceeds these bounds fails atomically rather than publishing partial
records as a complete collection. The native MCP profile has separate smaller
transport limits; consult each tool's schema and error for its active limits.

**Keep original source in workspace** retains exact bytes as a content-addressed
`sam-source-v1` or `bam-source-v1` attachment. Open it again from the history
inspector. Every opening reparses the source in a fresh worker; archive metadata
and earlier reference-check reports restore no scientific authority. Removal
removes the workspace source entry. Original download preserves source bytes;
SAM/BAM conversions can change framing or representation and carry explicit
notices. A checked report is evidence for its named source/reference/options,
not a replacement for those inputs.

MCP exposes `import_read_alignments`, `inspect_read_alignments`,
`check_read_alignment` and `export_read_alignments` against the same immutable
collection and reference checker. Binary worker schema 24 appends source
prepare/page/check/export/dispose messages; prior message tags remain stable.

This implementation is not yet browser-validated. Browser checks remain deferred
until all in-scope implementation is complete, as requested. CRAM, indexes,
large-genome range access and coverage/variant inference remain outstanding.
