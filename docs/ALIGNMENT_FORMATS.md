# Independent alignment exchange, profile 1

`collagene_core::alignment_formats` exchanges existing multiple alignments. It does not run or implement MAFFT, MUSCLE, Clustal, Primer3, or another aligner. The implementation is independent Rust, with no scientific subprocesses or new dependencies, and is portable to native Windows/macOS/Linux and WASM.

## Public entry points

- `import_alignment(bytes, AlignmentFormat, Alphabet, name)` requires an explicit DNA, RNA, or protein alphabet. Ambiguous letters never determine the alphabet.
- `export_alignment(&document, format)` returns canonical text, explicit `losses`, and a `renamings` map. Current profiles preserve identifiers and report no renamings; identifiers conflicting with Stockholm framing must be explicitly changed by the caller.
- `document.original_bytes()` returns the exact imported bytes after validating their BLAKE3 digest. It remains available after edits.
- `save_alignment` / `load_alignment` use bounded JSON attachments, schema version 1, strict unknown-field rejection, source hash validation, canonical residue validation and source-derived stable row IDs.
- `document.column_to_ungapped(row_id, zero_based_column)` returns no position for a gap, or an ungapped row index with orientation `AsStored`. `row_position_map(row_id)` builds a whole row mapping in one linear pass for viewers. Neither sequence identifiers such as `name/12-99` nor alignment columns become external reference coordinates.
- `document.extract_rows(&[row_id])` creates ordinary linear molecules with gaps removed. It copies no features, original GenBank record, alignment scores, consensus, GR/GC annotations, or alignment confidence. An all-gap row cannot become an empty molecule.

Every row has a stable source-derived ID and source row index, a separate editable identifier, an optional FASTA description, and uppercase ASCII aligned residues with `-` gaps. All current rows have the same nonzero column count. All-gap rows are supported in alignment documents. Duplicate identifiers are rejected rather than conflated.

`source_is_current()` reimports the exact original bytes using the retained import alphabet and name and compares the full document. No serialized boolean asserts that original evidence is current. Current row IDs must still derive from the retained source, so adding independently originated rows requires a new import. Subsets/reordering/editing can be stored with provenance retained. Changed row sequences/identifiers or a subset must clear the old CLUSTAL consensus and GR/GC annotations before validation, avoiding stale column evidence. GF/GS free text is retained as source annotation, not verified biological fact. A hash establishes integrity of retained bytes, not authenticity of a remote publisher.

## Accepted formats

**Aligned FASTA**: `>` header with unique first-token identifier and optional description; wrapped residue lines; empty lines permitted; `-` is the only gap glyph. Residue lines cannot contain internal spaces, digits, comments, dots, or silently discarded punctuation. Lowercase residues become uppercase while original bytes remain available. Unequal row lengths and empty rows are errors. FASTA export wraps at 80 columns and preserves descriptions; Stockholm annotations and CLUSTAL consensus cannot be represented and are reported as losses.

**CLUSTAL**: first line `CLUSTAL`, `CLUSTAL ...`, or `CLUSTALW ...`; ordered equal-width sequence blocks; optional checked cumulative ungapped residue counts; optional space-indented consensus using `*`, `:`, `.`, and spaces. Sequence residue columns must agree within a block. Consensus must line up with the sequence column and is preserved as imported symbols, not recalculated or assigned nucleotide conservation significance. Missing consensus contributes spaces. Tabs in consensus are rejected because their displayed position is ambiguous. Gaps are `-`; a dot in a sequence is rejected. Repeated blocks must contain every row in the same order. Export uses 80-column blocks without residue counters; row descriptions and Stockholm metadata are reported as losses. Long identifiers are never truncated.

**Stockholm 1.0**: exact first line `# STOCKHOLM 1.0`, final `//`, one record only, and no nonblank trailing content. Sequence gap `.` is normalized to `-` per the specification; dots inside annotations remain dots. Wrapped blocks are supported with blank boundaries or a repeated first row; all blocks must have identical ordered row sets and widths. GF and GS free-text entries, including repeated entries, remain ordered. GR and GC annotations must have exactly one ASCII non-whitespace character per alignment column, refer to existing identifiers, and occur once with identical annotation keys in every block. Duplicate GR/GC keys within a block and inconsistent annotation coverage are errors. Plain `#` comments are retained; unknown `#=` classes are explicit unsupported errors. Export retains GF/GS/GR/GC and comments, using 80-column blocks; FASTA descriptions and CLUSTAL consensus are reported as losses. Annotation placement/line wrapping may change, so original download is the exact round-trip route.

Export always reports that canonical output changes source layout/case/header and does not encode the application document name. The importer is deliberately stricter than permissive legacy tools that discard unsupported characters. Inputs outside this profile fail explicitly; they are not silently repaired.

## Bounds and validation

Maximum source size: 8 MiB; rows: 1,000; total aligned cells: 8,000,000; line size: 65,536 bytes; source annotation lines and retained annotation entries: 20,000; retained annotation text budget: 8 MiB; JSON attachment: 32 MiB. Output is subject to the source and line limits. Errors do not return truncated alignments. Input must be ASCII UTF-8 text with LF/CRLF line endings; unsupported control characters are rejected. Parsing and exports never access files, network, or executables.

## Application and workspace integration

**Alignments** opens a separate retained collection view. The main file picker routes `.aln`/`.clustal`, `.sto`/`.stockholm`, `.afa`/`.afasta` and saved `.cg2aln` attachments there. Use **Choose alignment file** or paste within Alignments for ordinary `.fa`/`.fasta` files that are intended as alignments; equal-length ordinary sequence records do not silently become an alignment. Choose the supplied alphabet explicitly. Saved native alignments retain their recorded alphabet.

The two-axis virtual Canvas2D view paints only visible cells and keeps row labels fixed. Arrow/Home/End/Page keys and clicks inspect a cell; the readable DOM view reports alignment columns, gaps, ungapped row positions and observed symbol counts. A compact immutable prefix index answers row-coordinate queries without rescanning full rows during painting. Native CSS scroll extent is capped while logical coordinates cover the entire supported column range. GR/GC and CLUSTAL symbols are retained imported evidence, never recalculated confidence. Row and annotation lists are paginated; full native/source downloads retain omitted preview text.

**Open selected ungapped rows** creates ordinary linear DNA/RNA/protein documents in stored orientation. An invalid or all-gap selected row rejects the whole extraction. Source row IDs and current/source alignment hashes survive as provenance; no external reference mapping, features or alignment confidence are invented. Canonical exports are prepared in the worker and show losses before download. Original-source download retains exact bytes. Source/options/workspace changes and unmount cancel pending worker jobs and reject obsolete replies.

Collections are stored as hashed `alignment-document-v1` workspace attachments. Recognized attachments must pass strict native/source validation before a workspace can save/load or replace browser storage. They survive autosave, recovery and full workspace exchange. This web collection is bounded to 128 retained alignments and the existing archive/transport limits. Source integrity is not publisher authenticity or alignment correctness.

Schema 10 binary workers and the session-only MCP server support alignment import/read/column inspection/export/save/extraction. MCP workspace imports expose recognized alignment attachments atomically. Scientific interfaces and selection do not grant persistent application writes. See [worker protocol](WORKER_PROTOCOL.md) and [MCP tools](MCP.md).

## Primary format references

- [Stockholm format 1.0 specification, Erik Sonnhammer](https://sonnhammer.sbc.su.se/Stockholm.html): framing, dot/hyphen gaps, GF/GS/GR/GC annotations, wrapping, repeated free-text metadata.
- [Official Clustal X help, Clustal authors](https://www.theclustal.com/download/clustalx_help/): CLUSTAL recognition, hyphen gaps, consensus glyphs, and optional residue numbering. We use an explicitly stricter input profile and require caller-supplied alphabet rather than the program's heuristic alphabet guessing.

References reviewed 2026-11-24. No implementation from another cloning application or alignment library was copied.
