# Reviewed sequence compatibility controls

The File menu adds **Export SnapGene** and **Export GCG** to the existing document
exchange controls. Both prepare an exchange file against the exact selected
source and workspace context; saving checks that context again. The source
molecule, native history and retained originals remain unchanged.

SnapGene export displays the selected-profile assessment before permitting
omissions. Accepting omissions cannot bypass the physical-material restrictions.
Only fresh supported reference packets are emitted; edited files do not receive
stale opaque packets from their original import. The profile is not a claim of
full SnapGene history or application compatibility. See [profile](SNAPGENE_EXPORT.md).

GCG exports require a record ID and documentary date, with the actual UTC date
initially shown for review. Acknowledgment covers the listed annotation,
topology, history and original-evidence omissions. The native file remains the
lossless application format. See [profile](GCG_SEQUENCE_FORMAT.md).

Opening recognized GCG text in the file picker opens a separate review with an
explicit DNA/RNA/protein choice. GCG's nucleotide type does not distinguish DNA
from RNA. The parser checks type, length, numbered positions and checksum.
This review imports the whole selected source atomically after a second context
check; multi-file batches requiring this choice are refused before publication.

FASTA sources use the same review panel with a record-declaration mode and strict
DNA/RNA/protein choices. A bounded summary shows the actual imported alphabet and
length. Explicit record declarations must agree with a forced choice; conflicting
name labels cannot override the native export alphabet. The paste dialog also
offers **Choose FASTA alphabet**. [Alphabet contract](FASTA_ALPHABETS.md).

Ordinary original Collagene XML continues through its existing importer. A project
containing lower-strand-only records opens an explicit migration review instead.
The opt-in profile reverses the actual lower calls into their 5′→3′ sequence;
it does not invent a complementary strand. All original rows, source bytes and
record indexes remain retained, and uncertain annotation mappings remain marked
for review. See [migration profile](COLLAGENE_XML_COMPAT.md).

Native/embedded MCP exposes text profiles through `import_compat_sequence` (4 MiB)
and explicit binary profiles through `import_sequence_file` (8 MiB decoded).
Binary choices are SnapGene DNA/protein, XDNA v0 with a supplied document name,
ABIF, SCF and current native `.cg2`. Wrong profiles, malformed source data and
oversized payloads refuse before any session documents are installed. Supported
original attachments and parser notices remain in the immutable native documents.

`export_snapgene`, `export_gcg` and `export_xdna` use the same core reference writers
as the browser. `export_trace_scf` regenerates SCF 3.00 from actual retained SCF
signals, original calls/peaks, four-channel accuracies and declared orientation;
it cannot turn an edited sequence or ABIF scalar scores into a chromatogram.
Precision is explicit. Missing/changed source evidence and precision loss refuse.
Exports return `review_required` without an artifact when documentary omissions
need acknowledgment, then use the existing base64 artifact reader after acceptance.
These tools accept no filesystem paths, invoke no external programs and leave
source documents unchanged. The native registry has 203 tools. Binary worker,
scientific and workspace schema versions are unchanged by this MCP-only batch.

The adapters and pure/source-binding fixtures are separate from application
acceptance. Browser checks remain deferred until all in-scope code is complete;
these controls have not been visually inspected in a browser. Current evidence
and remaining proprietary-container gaps are recorded in the
[delivery matrix](FORMAT_DELIVERY_MATRIX.md).
