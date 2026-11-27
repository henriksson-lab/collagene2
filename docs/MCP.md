# Local MCP server

`collagene-mcp` exposes the portable scientific core through the official Rust `rmcp` SDK. It uses standard newline-delimited UTF-8 JSON-RPC over stdio, separate from the application's binary worker protocol. stdout is reserved for MCP; startup errors go to stderr. The executable never starts scientific programs or AI clients.

The current registry has 190 tools. Internal exchange uses only worker 100, candidate 26, native workspace plan 24, verifier `science-27`, and assembly draft 15. Noncurrent internal encodings are rejected; they are not migrated or relabeled. External sequence-file profiles keep their own explicitly documented versions. Selected-construct design and batch operations
add exact-version matrix paging, native inspection, historical capture/list/open
and explicit fresh replay. [Batch contracts](COMBINATORIAL_BATCHES.md).
Reference recombination retains explicit profile discovery, source inspection,
independent target checking and fresh selected-component reconstruction. The
optional `plan_sources` bridge binds a returned selection to typed plan sources
without requiring a client to reproduce internal hash recipes.
[Recombination tool contracts](RECOMBINATION_MCP.md).

Sequence-file exchange provides `import_compat_sequence`, `export_snapgene` and
`export_gcg`. Imports require explicit profiles/alphabet. Lossy exports return a
review without an artifact until omissions are acknowledged, and restrictive
physical or modified-protein evidence cannot be erased through that choice.
[Selected profiles and controls](SEQUENCE_COMPATIBILITY_CONTROLS.md).

Build and configure:

```sh
cargo build -p collagene-mcp --release --locked
```

A typical MCP client's configuration is:

```json
{
  "mcpServers": {
    "collagene": {
      "command": "/absolute/path/to/collagene2/target/release/collagene-mcp",
      "args": []
    }
  }
}
```

On Windows use the absolute `.exe` path, escaping backslashes as JSON requires. The client owns process startup and shutdown. No server configuration files, database, network endpoint or authentication secret is needed for this session-only stdio mode. A client can access only documents it explicitly imports into that session, not browser storage or local files.

Optionally restrict managed native genome storage to a chosen existing folder:

```json
{"mcpServers":{"collagene":{"command":"/absolute/path/collagene-mcp","args":["--storage-root","/absolute/path/CollageneData"]}}}
```

Create and protect that folder beforehand. The server canonicalizes it, rejects a file/missing folder/final symlink or junction, and creates a unique per-session child beneath the opened directory capability. All managed genome staging, readback and export I/O uses retained file/directory handles there. Small documents, plans and artifacts remain memory-only. Default startup still uses private OS-temp session storage. `--help` and `--version` write only to stderr; invalid arguments/storage configuration return a structured JSON diagnostic on stderr and exit2, preserving stdout exclusively for MCP.

This option does **not** enable arbitrary file imports/exports, access to existing files in that folder, a public filesystem HTTP API or browser access to the native root. Upload/download chunks and opaque session handles remain the only file-backed interface. It limits managed application I/O, not all OS access by other same-user processes. Unix creation modes are0700/0600; Windows privacy follows the configured parent's ACL. Cleanup is best-effort for owned files and empty session directories on normal disposal; crashes, foreign handles or external namespace changes may leave data behind. Keep the parent folder protected. The reusable implementation is [collagene-local-storage](../crates/local-storage/README.md).

## Browser exchange

The web Plan panel provides an explicit file exchange without exposing a network bridge:

1. Open your input documents and select the desired target. Click **Export for AI**. The `.cg2design` file contains exact document snapshots, including topology, features, original evidence and revisions; duplicate revisions are included once. It is a planning snapshot, not a backup of undo history or a complete project archive.
2. Give that file to an MCP client with file-reading capability. Call `import_workspace` with its bytes encoded as base64. The response includes exact document references and the target selected in the browser. The initial stdio exchange supports at most 8 MiB of decoded snapshot data; the web exporter enforces the same limit.
3. The client proposes a plan bound to those references and may call `verify_plan`. Call `export_candidate` with the exact plan revision. The response identifies a base64-encoded binary artifact; concatenate `read_artifact` text chunks until `done`, decode the concatenated base64, and save as `.cg2plan` using the client's file facilities.
4. Back in the web workspace, use **Open candidate**. Inspect bound inputs and ordered scientific steps. Select the intended target and click **Verify candidate**. The browser binds the candidate to its current open target/input revisions and reconstructs it in a Rust/WASM worker. Missing or changed revisions block the check. The imported file never replaces open documents, applies a product or imports an asserted verification status.
5. Export the HTML or embedded PDF protocol after the check. **Save candidate** preserves a portable proposal with no trusted status. Candidate proposals also survive browser autosave/recovery; computed reports must be regenerated after reload. **Discard candidate** returns to the manual overlap editor.

A `.cg2design` file can also be opened directly in the web app to recover its explicit document snapshots. This is an intentional file-open action. Opening only a `.cg2plan` does not silently restore missing inputs or an old target; open those original document revisions explicitly, or export a fresh snapshot and redesign.

Both exchange formats are experimental versioned binary envelopes with bounded Postcard payloads, BLAKE3 integrity checks and model validation. A checksum detects corruption; it does not confer authority or scientific validity. MCP's required JSON transport uses base64 for this binary file boundary; worker communication remains binary.

## Full native workspaces and chunked import

`import_workspace_archive` accepts `{data_base64}` for a complete `.cg2workspace` archive up to 8 MiB decoded. This is separate from the original planning-only `import_workspace` tool. It preserves the exact document identities, all retained revisions, history graph and undo/redo stacks, active document/candidate, candidates and opaque attachments. Imported candidates remain unverified data. To compute one as a session plan, explicitly rebind its inputs and independent target using `propose_plan`; archive import does not register a trusted report or silently execute a plan.

For larger native archives:

1. Call `begin_workspace_upload` with the exact decoded `expected_bytes` and lower-case `blake3` of the whole archive. Native archives can contain up to 64 MiB compressed and 128 MiB expanded; the separate session limit still applies. The response returns `upload_id`, `next_offset: 0` and the maximum decoded chunk length.
2. Repeatedly call `append_workspace_upload` with that ID, the exact next decoded `offset` and one independently base64-encoded chunk of 1–262,144 bytes. Every accepted response returns the next offset. Overlaps, replayed chunks, gaps, invalid base64 and overrun are rejected without changing the upload.
3. Call `finalize_workspace_upload`. Length and whole-file hash are checked before parsing. Core validation and session size/item checks complete before any imported documents become visible. Incomplete uploads remain appendable. A complete finalize attempt consumes the upload and releases its reservation even when hashing, parsing or session storage fails; retry those cases with a new upload. `abort_workspace_upload` explicitly releases a pending reservation without importing anything.

At most four uploads can be pending. Each reserves its full expected decoded byte count against the session's shared 64 MiB budget at begin, preventing unaccounted append growth. Hashing/decompression retain that reservation until parsing completes. Native workspace structure is stored separately from shared immutable molecule revisions: reimporting a workspace deduplicates both, and history does not require another full raw archive snapshot. Metadata and expanded documents can exceed the session budget even when the archive passes native compressed/expanded limits; that fails explicitly and atomically rather than dropping history, candidates or attachments.

`list_workspaces` returns immutable workspace IDs/revisions and bounded candidate summaries. `export_workspace_archive` requires `{id, revision}` and reconstructs the complete native archive. Its binary artifact stores raw bytes; `read_artifact` encodes only the requested base64 text interval, including arbitrary non-four-character-aligned offsets. Concatenate chunks before decoding. No file path, external writer or executable is accepted. Existing and future source documents are not replaced by an archive import.

`list_documents {start?, limit?}` returns the first 128 summaries plus `total` and `next_start`. Supply `{start, limit}` for later pages, with a limit of 1–128. A session can hold up to 10,000 distinct molecule revisions, accommodating native workspaces with 1,000 document identities and retained history. Plan/artifact/workspace/oligo-collection limits remain separate.

## Caller reads and oligo collections

`import_fastq` takes `{text, encoding}` with required `Phred33` or `Phred64`; overlapping ASCII quality ranges are never guessed. `export_fastq` takes an exact `document` and a required target encoding. Stale mappings or scores unrepresentable in that encoding fail explicitly instead of fabricating or clamping confidence. Reverse/trim transformations retain the original quality source and coordinate mapping. Each source import is bounded to 8 MiB.

`validate_reads` takes `reference`, ordered `reads` document references and explicit `TraceValidationOptions`. It reconstructs SCF/ABIF/FASTQ quality, trimming, pairwise placement uniqueness, coverage and differences against the independently selected reference. The response contains bounded confidence/placement/coverage/difference previews and a full revision/options/engine-bound JSON artifact. Native ABIF/SCF confidence is not assumed calibrated Phred. Uncovered or low-confidence positions are not filled from the reference, and partial or ambiguous placement cannot establish whole-plasmid confirmation. See [read-validation model and limits](TRACE_VALIDATION.md).

`import_oligos` accepts CSV/TSV `text`, explicit `options` mapping, and optional `selection` of zero-based source data-row indices. The full table is retained once as an immutable collection, with original rows, unknown columns, modifications and quantity units. Explicit selection sets are separately retained and deduplicated; selected supported primers become ssDNA reference documents carrying only their own original row and the table hash. Unsupported modified rows remain collection evidence and cannot be silently converted to DNA. Selecting such a row rejects the entire import without partially creating documents.

`export_oligos` takes the exact `collection: {id, revision}`, an explicit delimiter and `SpreadsheetSafe` or `MachineData` mode, with optional source-row selection metadata. It returns three bounded artifacts: exact original bytes, a normalized complete CSV/TSV table, and full JSON evidence containing the table, requested selection and saved import selections. The normalized table contains every collection row; selection is preserved explicitly in JSON rather than by inventing a smaller source table. Spreadsheet-safe escaping and machine-data behavior remain explicit, and modification strings are retained. Retrieve large artifacts through `read_artifact`.

## Design workflow

1. `import_document` with GenBank, FASTA or plain DNA contents. Inspect import warnings, including linear topology defaults for FASTA/plain DNA. The response returns immutable document IDs, revisions and bounded summaries. Reimporting identical documents is idempotent. Imports are atomic.
2. `list_documents` and `read_sequence_region` inspect those revisions. Read coordinates are zero-based half-open, and reads cannot wrap across a circular origin. Split origin-spanning reads into two regions.
3. `list_enzymes`, `analyze_restriction`, `digest_fragments` and `simulate_pcr` expose the same deterministic models as the browser. No scientific result is supplied by an AI. PCR primers contain `name`, `annealing: {alphabet: "Dna", bases: "..."}` and a 5′ `tail` string. Both annealing sequences are written 5′→3′.
4. `propose_plan` binds an independent target and explicit inputs. It returns `state: "candidate"` and `verified: false`. A `verified` field supplied by the client is rejected. Duplicate proposals return the same plan revision.
5. `verify_plan` accepts that exact plan revision. The core reconstructs products, checks input hashes/dependencies and compares the final product with the separately bound target. Its status is `Invalid`, `Incomplete`, or `ConsistentUnderAssumptions`. None establishes experimental success or exhaustive product uniqueness.
6. `export_protocol` re-verifies the exact saved plan, creates a session HTML artifact and returns its ID. `read_artifact` reads bounded UTF-8 chunks; concatenate their `text` fields using `next_start` until `done` is true. Save the HTML with the MCP client's artifact/file facilities and use browser printing for PDF. For binary candidate artifacts, decode base64 after concatenation instead of saving the base64 text as a binary file. Invalid/incomplete drafts retain their findings visibly.

For the synthetic `fixtures/overlap-demo.fasta`, replace the ID/revision placeholders below with the values returned by import. Current operation schemas are generated from Rust and also exposed by `tools/list`.

```json
{
  "name": "Two-fragment example",
  "target": {"id": "TARGET_ID", "revision": "TARGET_REVISION"},
  "inputs": [
    {"id": "left", "document": {"id": "LEFT_ID", "revision": "LEFT_REVISION"}},
    {"id": "right", "document": {"id": "RIGHT_ID", "revision": "RIGHT_REVISION"}}
  ],
  "steps": [
    {"Overlap": {
      "fragments": [{"Input": "left"}, {"Input": "right"}],
      "overlaps": [4],
      "circular": false
    }}
  ]
}
```

`{"Step": 0}` refers to the first prior reconstructed step. A PCR operation is `{"Pcr": {"template": {"Input": "template"}, "forward": PRIMER, "reverse": PRIMER}}`. Forward/cyclic dependencies are rejected during verification. A syntactically valid candidate may be scientifically invalid; saving it does not promote its status.

## Translation and pairwise tools

`translate_sequence` takes an immutable `document` reference plus `TranslationOptions`; `translate_cds` takes that reference and a zero-based `feature_index`, using the annotation's compound location, `codon_start` and `transl_table`. `find_orfs` takes a document and explicit `OrfOptions`. The initial genetic codes are 1, 2 and 11. Ambiguous/partial predictions and unsupported exceptions remain explicit; a stored translation is compared as evidence rather than trusted.

`align_sequences` takes `reference` and `query` document references plus `AlignmentOptions` for global, local or query-in-reference affine alignment. `compare_molecules` uses `ComparisonOptions` to explore bounded reference-origin and reverse-query candidates. Inspect its evaluated candidate counts, exhaustive-search and optimality flags before interpreting the best result. The model never labels a budget-limited origin search exhaustive.

All five tools return bounded previews/summaries and an `artifact_id` for the full JSON report. Use `read_artifact` to retrieve it in chunks. Reports bind source revisions, options and scientific engine version; these tools do not mutate documents or confer verification status on AI-proposed cloning steps. Current parameter schemas are available through `tools/list`.

## Primer and PDF tools

`analyze_primer` takes a document reference, `primer`, explicit `conditions` and `binding_options`; an optional `partner` enables heterodimer screening. Reports distinguish annealing-domain thermodynamics from full ordered-oligo structure screens, and exact 3′ runs from mismatched/ambiguous compatibility. `design_primers` takes a document and `DesignOptions` for anchored endpoint selection with separate tails, explicit filters and reproducible ranking. Both return bounded previews and full JSON artifacts. These calculations do not establish specificity, amplification or experimental success. See [primer models](PRIMERS.md).

`export_protocol_pdf` takes the same exact `plan_revision` as `export_protocol`. It re-verifies the plan and independent target, then returns a standalone PDF artifact with `encoding: "base64"`, byte/page counts, status, bound revisions, engine version and rendering warnings. Read its text chunks with `read_artifact`, concatenate and base64-decode them into `cloning-protocol.pdf`. Generation is embedded Rust; no external printer, viewer or process is invoked. Invalid/incomplete drafts stay marked on every page. See [PDF limits and font notice](PDF_EXPORT.md).

## Restriction/ligation and URL tools

`digest_fragments` accepts `{document: DOCUMENT_REF, enzymes: ["EcoRI", "BamHI"]}` and returns complete-digest summaries in ascending top-cut order: separate strand lengths, overhang polarity and 5′→3′ protruding sequences, terminal chemistry and source boundaries. Uncut circular input remains circular. Known source chemical/strand evidence that is not modeled is incomplete, not silently assumed to be ordinary duplex DNA.

Plans also accept:

```json
[
  {"Digest": {"template": {"Input": "vector"}, "enzymes": ["EcoRI"], "fragment_index": 0}},
  {"TreatEnds": {"template": {"Step": 0}, "treatment": "Dephosphorylate"}},
  {"Digest": {"template": {"Input": "insert"}, "enzymes": ["EcoRI"], "fragment_index": 0}},
  {"Ligate": {
    "fragments": [
      {"source": {"Step": 1}, "reverse": false},
      {"source": {"Step": 2}, "reverse": true}
    ],
    "circular": true,
    "allow_nicked": true,
    "options": {
      "max_enumeration_fragments": 8, "max_path_checks": 4096,
      "max_products": 100, "max_total_product_bases": 8000000,
      "max_report_bytes": 33554432, "max_work": 20000000
    }
  }}
]
```

These steps require separately bound `vector` and `insert` inputs; step indexes refer to earlier reconstructed products. Forward/self-references are invalid. Digest fragment indexes are zero-based in the API, while the manual web form numbers them from one. `Phosphorylate` is the other terminal treatment. Ligation reconstructs both strands and accounts for junction bonds; `allow_nicked` does not repair a nick or make a non-ligatable junction valid. Verification returns end/nick summaries, and incomplete/invalid plans retain already reconstructed intermediates. The intended path and independently reconstructed alternatives retain actual DNA, topology, ends, nicks and modification summaries. Enumeration limits do not change target acceptance or imply isolation. This is a bounded one-copy path model, not probabilities or exhaustive real-mixture behavior. PCR on a known staggered/nicked physical intermediate is explicitly incomplete until a two-strand PCR template model is available.

`import_url` accepts `{url: "NC_001422.1"}` (or a supported URL) and imports bounded source records through the native Rust HTTP client. NCBI key/contact may be supplied through the server process environment; they are never stored in sequence metadata or exported URLs. For Addgene, configure `COLLAGENE_ADDGENE_API_TOKEN` with approved catalog/sequence access. Call `list_source_sequences` with the plasmid URL/ID, review its full/partial/deposited/Addgene choices, then call `import_url` with that URL and an explicit `sequence_id`. The catalog is refreshed and downloaded sequence is checked against available catalog evidence. Choice previews are capped at 100 with 512-character descriptions; the full safe public summary is a JSON artifact. No token, raw sequence or transient download URL appears in that summary.

After native URL import, `export_document` takes `{id, revision}` and returns a base64 `.cg2` artifact. Read chunks, concatenate and decode before opening it in the browser. This preserves the exact fetched revision, annotations and evidence when CORS blocks direct import. [URL adapters](URL_IMPORT.md) documents policies and unverified live Addgene interoperability.

The `.cg2plan` candidate and native workspace plan use the current internal encodings listed above. Nested plans, capsules and worker requests/replies must satisfy those same current encodings, integrity checks and aggregate resource bounds. Unsupported versions fail explicitly; no engine identity is synthesized. Single-molecule `.cg2` has its own molecule format version and does not grant retained physical history. See [worker protocol](WORKER_PROTOCOL.md).

## Golden Gate, mutagenesis and gels

`simulate_golden_gate` takes immutable original `documents` and `GoldenGateOptions`. It starts with unknown original modification state, which can withhold undetermined products. Use the prepared/whole-reaction inspection tools for explicit generation assumptions or saved native material origins. It returns bounded directed candidate paths, stable physical sequence-group counts, recut/carryover evidence, finite-search completeness and a full report artifact. To propose a specific construction, use `GoldenGate {templates, options, candidate_index}` in a separately target-bound plan. The verifier reconstructs it again from current inputs or trusted prior physical steps. A site-free candidate does not establish kinetic success or exhaustive mixture uniqueness.

`simulate_mismatch_pcr` takes `document`, actual `forward`/`reverse` primers and `PcrOptions`. Its artifact retains introduced source/product changes, both binding orientations, overlap/ambiguity evidence and completeness flags. `PcrWithMismatches {template, forward, reverse, options}` is the corresponding verified step; uncertain/overlapping or multiple-product cases cannot silently become a unique verified duplex.

`design_inverse_pcr` takes an immutable circular `document`, independent circular `target` and `InversePcrOptions`. It returns actual primers, computed diagnostics and a stored candidate `plan_revision` with explicit PCR, terminal phosphorylation and sealed circular ligation. Review/reverify that candidate before export/use. This stores a session proposal, not a project edit or experimental result. Full evidence is a JSON artifact.

`predict_gel` takes bound `lanes` (label, document, preparation, optional loaded mass) and an explicit `Calibration`. It returns full JSON and standalone SVG artifact IDs, with bounded lane summaries. Uncut circular conformation, unsupported/error profiles and out-of-calibration fragments remain visible and unpositioned. `design_diagnostic_digest` takes immutable `documents` and `DiagnosticOptions`; ranked enzyme sets distinguish supported band profiles under that calibration and resolution, while unknown profiles stay unknown. Index-pair results refer to the ordered document list/source revisions. Search caps and completeness accompany all reports. There is no observed-gel comparison or OpenGel integration.

Model details: [Golden Gate](GOLDEN_GATE.md), [mutagenesis](MUTAGENESIS.md), [gel prediction](GEL_PREDICTION.md). All five tools use the same embedded core and immutable revision checks; schemas are available through `tools/list`. Export artifacts and proposals do not acquire browser/project-write authority.

## Reusable feature libraries and annotation exchange

These tools provide session-only, immutable evidence handling:

- `create_feature_library` takes `{name, version, features: [{document, feature_index}]}`. Each source is an exact imported `DocumentRef`; biological extraction handles supported compound/reverse locations. A collection may combine features from several documents. Duplicate selections, unsupported extraction or size limits reject the entire operation. No proprietary feature database is bundled.
- `import_feature_library` takes `{data_base64}` containing validated library JSON, at most 8 MiB decoded. `list_feature_libraries` lists bounded summaries. `export_feature_library` takes `{library: {id, revision}}` and returns the complete source/version/entry JSON artifact.
- `scan_feature_library` takes `{document, library, options}`. `ScanOptions` explicitly controls strands, allowed mismatches, comparison budget and hit budget. It matches whole entries only: coverage is one, while concrete matches, ambiguity-compatible columns and definite mismatches give distinct identity bounds. Circular origins, reverse orientation, repeated sites and incomplete bounded searches remain explicit. The result provides an immutable `scan: {id, revision}`, a bounded preview and a full JSON artifact.
- `apply_feature_library` takes `{document, library, scan, selected: [proposal_id]}`. The original target/library revisions must still match. The core reconstructs the proposal set before applying selected IDs; imported reports cannot assert coordinates or successful verification. It returns a **new** immutable document, appends only selected supported features and marks inferred annotations for review. It does not modify the source document or grant physical cloning authority.
- `import_annotations` takes `{document, text, format, seqid}` with explicit `Gff3`, `Gtf` or `Bed` and source sequence ID. It returns an immutable `dataset: {id, revision}` and bounded proposal/coordinate previews. It never guesses a sequence-name mapping or installs rows automatically. Original source bytes, unrelated sequence IDs, row evidence and format warnings remain retained. Input text is limited to 8 MiB.
- `apply_annotations` takes `{document, dataset, selected: [proposal_index]}`. It validates the same reference revision and original dataset again, then returns a **new** document with selected annotations appended and source evidence retained. Existing annotations are not silently replaced.
- `export_annotations` takes `{document, format, seqid}` and returns a full GFF3/GTF/BED text artifact with a format-specific extension plus loss diagnostics. GTF uses quoted attributes and explicit gene/transcript identifiers; it does not infer transcript joins or hierarchy. `export_annotation_dataset` takes `{dataset}` and returns the exact original source bytes and complete JSON proposal evidence as separate artifacts. The strict format profile rejects unsupported phase/relationship/location semantics instead of inventing them.

Libraries, scans and annotation datasets each have a separate 128-item cap within the existing shared 64 MiB serialized session budget. Full scan reports and their text artifacts are accounted separately when both are retained. Creation/application is atomic: invalid selections, stale references or capacity failures do not partially create document/library entries. Previews include at most 100 matches/proposals, truncated labels and explicit truncation flags; use artifacts for complete evidence. These tools accept no paths, file writes or claims of verification.

Core semantics: [feature libraries](FEATURE_LIBRARY.md) and [annotation format profiles](ANNOTATION_FORMATS.md). Worker envelope version 7 introduced library scanning/application and annotation import/application/export; version 8 added physical end-repair/annealing previews and materialization, while preserving supported earlier message layouts and the ordinary 32 MiB message cap. Full native workspace jobs retain their separate bounded large-message profile.

## Resources, limits and authority

`resources/list` exposes immutable `collagene://document/REVISION`, `collagene://plan/PLAN_REVISION` and `collagene://report/PLAN_REVISION` URIs. Document/plan resources are summaries; sequence regions use the bounded read tool. Report reads recompute the same bound plan using the running engine version. No filesystem or HTTP resources are accepted. Annotations, names and original evidence are data and cannot authorize tool calls.

The server limits ordinary text imports to 4 MiB; FASTQ/oligo source text to 8 MiB; document revisions to 10,000; plans, artifacts, imported workspaces, oligo collections, saved oligo selections, feature libraries, feature scans and annotation datasets to 128 each; accounted serialized storage and pending-upload reservations together to 64 MiB; generic `propose_plan` requests to 100 steps/128 inputs/2 MiB; sequence reads to 10,000 bases; artifact reads to 16,384 text bytes; and restriction/PCR summaries to the first 100 entries with explicit truncation flags. Binary artifacts are accounted by their stored raw byte count, not an eagerly retained base64 expansion. The accounting bounds stored serialized content, not total process RSS. The SDK bounds incoming stdio lines at 16 MiB; normal scientific-core complexity/materialization limits also apply.

Scientific jobs use one blocking-compute slot per session so they do not block async transport handling or multiply materialized products concurrently. SDK cancellation can stop awaiting a job, but an already running bounded computation currently finishes; cooperative scientific cancellation/status remains to implement. Disconnecting loses all session data. Export artifacts before closing the client.

Live browser synchronization, persistent projects, scoped transactional `apply_plan`, richer searches/maps and broader diagnostic workflows, hosted Streamable HTTP with authentication/origin validation, and progress/cancellation are still pending. The manual snapshot/candidate exchange grants no implicit project-write authority and does not provide live synchronization with open browser documents. The manual exchange covers the first review/verification path; the broader hosted/persistent MCP contracts remain incomplete.

## Validation and dependency provenance

```sh
cargo test -p collagene-core -p collagene-mcp --locked
cargo build -p collagene-mcp --locked
python3 scripts/mcp_smoke.py
```

The integration tests use a real SDK client/server over an asynchronous byte stream, including 1,000-document native workspaces, exact history/candidate/attachment roundtrips, a native archive above 8 MiB transferred in chunks, strict upload offsets/hashes/reservations, stale revisions, explicit FASTQ quality and oligo modification evidence. Feature/annotation SDK tests also exercise independently supplied source features, immutable library JSON, exact/limited scans, unknown/tampered/stale proposal rejection, reverse GFF3 coordinates, BED blocks, exact original-source exports and explicit application to new documents. The separate development-only Python helper launches the actual executable and checks the standard stdio handshake, tools, immutable resources, candidate verification and HTML export. The helper's subprocess use is test infrastructure, not application runtime behavior. It also checks negotiation with the 2025-11-25 protocol version.

Collagene's application code is MIT. The [official Rust MCP SDK](https://github.com/modelcontextprotocol/rust-sdk) is an Apache-2.0 standalone dependency; it supplies protocol transport/routing, not cloning algorithms. The [MCP transport specification](https://modelcontextprotocol.io/specification/2025-11-25/basic/transports) defines the stdio boundary. No code from another cloning application was used.

## End repair and declared oligo annealing

`simulate_end_repair` takes an exact `document`, nonempty `enzymes` list, zero-based `fragment_index` and `EndRepairMethod` (`FillFivePrimeOverhangs`, `RemoveThreePrimeOverhangs` or `Blunt`). It reconstructs complete restriction cleavage from that source, then repairs the selected physical fragment. The bounded summary reports actual strand lengths, end chemistry, added/removed bases and annotation findings; a JSON artifact retains the full source-bound repair/ancestry evidence. It never accepts imported `duplex_model` metadata as a physical substrate. In a target-bound plan use `Digest` → `RepairEnds` → optional `Ligate` rather than reimporting a preview as an asserted end state.

`simulate_oligo_annealing` takes exact `top` and `bottom` document references and explicit `AnnealOptions`. Inputs are declared linear 5′→3′ ssDNA oligo designs; known double/mixed-strand declarations, circular inputs and imported physical-intermediate metadata are rejected. Concrete DNA only is supported. Options record `min_paired_bases`, `max_comparisons`, `max_candidates`, `top_chemistry` and `bottom_chemistry`. Each chemistry has `five_prime` (`Phosphate`, `Hydroxyl`, `Unknown`) and `three_prime` (`Hydroxyl`, `Blocked`, `Unknown`). Ordinary unmodified synthetic oligos use hydroxyl at both ends; caller-supplied chemistry is a plan assumption, not authentication of an experimental reagent.

The tool enumerates antiparallel placements, retaining competing offsets and completeness. Full-length/longest pairing is never silently selected. Only complete search with one candidate establishes unique placement under this model. The artifact includes both input revisions and full pair/mapping evidence. An `AnnealOligos {top, bottom, options, candidate_index}` plan operation selects one explicit offset. Alternative or unproven placement remains incomplete, even if that selected reference sequence matches the target. Annealing creates no covalent bond or implicit phosphate; `TreatEnds(Phosphorylate)` and `Ligate` remain separate operations when appropriate. A completely reconstructed product that misses its separately bound target remains invalid.

Physical preview documents carry documentary evidence only. HTML/PDF protocols include the explicit end-repair and annealing operations, actual strand/terminal state and stated assumptions. See [end repair and annealing models](END_REPAIR.md).

Original Collagene XML is accepted by `import_document` through the shared text parser, subject to the 4 MiB MCP text limit and explicit external-format profile guards. Binary XDNA can pass through the web/native workspace exchange; a direct generic binary-file upload tool is not yet provided. Original external-file attachments remain source data and cannot grant duplex chemistry.


## Unified search and browser tracks

`search_document { document, options }` runs readonly searches over an exact immutable document reference. Options explicitly select `All`, `Sequence`, `Annotations`, `Primers` or `Enzymes`, the query, nucleotide strand policy, annotation case/word matching, optional enzyme names and hit/row/work limits. DNA/RNA queries retain their alphabet and IUPAC certainty. Protein search uses literal residues without DNA coercion. Each hit binds its coordinates, feature index/enzyme, compound biological order, source revision and full options through a stable ID. Limits return explicit incomplete status. The tool previews at most 100 hits and stores the complete bounded report as a JSON artifact. See [search semantics](SEARCH.md).

`view_tracks { document, window, options }` uses the same bounded reference track engine as the web viewer. `window` has zero-based half-open `start`/`end`; circular windows can cross the origin once. Explicit options select annotations, saved primers, enzyme recognition and actual strand cuts, translations, validated CDS and ORFs. Source revision, engine/version, viewport, complete options, work/residue budgets and completeness are retained. The preview contains at most 100 items; the full report is a JSON artifact. Codon mappings preserve whole-reference frame/residue identity across viewport changes. No tool changes the selected region, installs annotations or upgrades search/track evidence to scientific verification. See [browser track semantics](VIEW_TRACKS.md).

`AnalysisSession` retains one validated immutable source/index; cached viewport requests carry their options/window and exact expected source revision. Missing, stale and replaced sources fail explicitly; replacing a prepared source is atomic. One-shot calls remain available. Resource caps and source identity apply to both routes.

## Alignment exchange and protein estimates

Six session-only alignment tools integrate the independent exchange profiles:

- `import_alignment` accepts exact `data_base64`, explicit `format` (`AlignedFasta`, `Clustal`, `Stockholm`), `alphabet` and `name`; decoded source limit 8 MiB.
- `import_alignment_document` restores a strict native JSON attachment as `data_base64`, decoded limit 32 MiB, validating retained source bytes and hashes before storage.
- `list_alignments` returns bounded immutable collection summaries. IDs/revisions hash canonical typed native JSON, including original evidence.
- `read_alignment` takes an exact `alignment: {id, revision}`, `row_start`, `max_rows` (1–50), `column_start`, and `max_columns` (1–1000). Columns and ungapped row positions are zero-based; gaps return null. These are as-stored row coordinates, not inferred external genomic coordinates. Source consensus remains imported evidence.
- `export_alignment` accepts `format` as an exchange profile, `Native` or `Original`. A UTF-8 artifact retains complete data; read all `read_artifact` chunks. Original exports preserve exact ASCII source bytes including CRLF. Native retains all supported metadata/source evidence; exchange exports disclose losses. `source_is_current` is independently derived, not an imported claim.
- `extract_alignment_rows` requires explicit stable `row_ids`, removes gaps and creates new immutable session molecule references. Unknown, duplicate or all-gap row selections fail atomically. It copies no fabricated external coordinates, confidence or feature annotations and never changes browser documents.

Alignment session evidence shares the 64 MiB accounting budget and 128-collection cap. Summary row previews cap 100 rows, names/identifiers/descriptions cap 256 characters, and larger collections remain accessible through paged reads/full exports. Importing a native workspace atomically exposes its retained alignment collections to these tools, with deduplicated identities and shared session accounting. It validates recognized `alignment-document-v1` attachments before installing any state; arbitrary future attachment kinds remain opaque. No tool runs MAFFT, MUSCLE, Clustal or any other installed aligner. Client/transport message limits still apply to direct base64 imports; the decoded file bound does not establish arbitrary JSON-RPC frame capacity.

`analyze_protein` accepts an exact immutable molecule `document` and `protein_properties::ProteinOptions`. It returns the report and a full JSON artifact: explicit stop/ambiguity/disulfide/modification assumptions, composition, available mass/charge/pI/extinction estimates, model IDs and unavailable-value diagnostics. This is a deterministic sequence-model calculation with no browser editing authority or claim of measured folding/activity.


## Explicit tailing and overlap-method checks

`simulate_tailing` accepts an immutable `document`, a strict `TailSubstrate` recipe, and `base` A/C/G/T. Digest recipes specify enzymes, a zero-based fragment index and optional explicit blunting; exact-PCR recipes supply actual primers and require a unique exact product. The bounded summary and full JSON artifact preserve both ordered strands, source revision, recipe, retained mappings and synthesized additions. This ideal one-base-per-3′-end transform never assumes phosphorylation, repairs internal nicks, or predicts Taq/TdT/kit activity. Unknown/blocked 3′ ends and unsupported geometry remain errors. The input sequence is not changed.

`check_overlap_profile` accepts immutable `documents`, equally sized `reverse` flags, ordered `overlaps`, `circular`, and `ProfileOptions`. Reports bind original source revisions/orientations, exact endpoint homology, original-coordinate spans, GC and advisory method guidance, plus bounded repeat evidence. Published Gibson/HiFi/In-Fusion recommendations do not prove reagent compatibility or empirical success and do not invalidate an otherwise exact DNA design merely because a length recommendation is unmet. The inline preview bounds junctions/occurrences/homology text; its full artifact retains the complete bounded report.

To propose a physical TA/GC design, use explicit `Digest`/`RepairEnds` or `Pcr` steps, `Tail3Prime`, optional `TreatEnds`, and `Ligate`. Both added junction bases and all 5′/3′ chemistry are reconstructed against the independently bound target; symmetric tails do not establish orientation or mixture uniqueness. `ProfiledOverlap` uses `OrientedSource` entries with explicit reverse-complement flags and the declared profile options. It performs source-bound checks and the same exact sequence/annotation reconstruction as ordinary overlaps. A known blocked 3′ physical terminus cannot be silently removed by a sequence profile. HTML/PDF protocols print these operations and their assumptions. Imported reports never confer verification authority.


## Retained cloning runs

`record_workspace_run` takes an exact immutable `workspace: {id, revision}`, a saved `plan_revision`, explicit `target_binding: {document_id, event_id, revision}`, `input_bindings: [{input_id, event: {document_id, event_id, revision}}]`, and optional `allow_incomplete` (default false). Use `list_workspace_runs` on the source workspace to inspect exact historical events; the full JSON artifact retains all events when its first100-event preview is truncated. Equal molecule hashes never select an event implicitly.

The tool freshly reconstructs every step against the separately bound target before recording. It returns a **new immutable workspace reference**, a receipt and a complete JSON artifact containing the resulting workspace. The source reference stays unchanged; nothing writes to a browser workspace or filesystem. The new workspace preserves original heads/navigation and adds separate product documents/events with actual operation parents. Incomplete evidence requires explicit opt-in and a complete reconstruction matching the target. Invalid or partially reconstructed plans cannot record products. Repeating the same recording deduplicates it.

`list_workspace_runs` returns at most100 run summaries with at most20 product bindings each, explicit truncation flags, bounded names and a full JSON artifact. `replay_workspace_run` takes `{workspace, run_id}` and recomputes the retained plan from its exact historical inputs with the current embedded verifier. It returns engine/status/product comparisons, at most100 product summaries and a full report artifact. Stored records and their historical status never establish current verification authority; replay evidence is bound to the complete supplied workspace fingerprint. Recognized `cloning-run-v1` attachments are structurally validated during native archive import. All session accounting and artifact limits above still apply.

## Codon design and independent checking

`design_codons` accepts an immutable `document` reference and strict `codon_design::DesignOptions`. The source is an explicitly supported protein or exact DNA coding sequence, with genetic code, terminal-stop semantics and limitations declared in the [codon design profile](CODON_DESIGN.md). Explicit integer codon preferences are not organism expression predictions. Returned candidates must pass a separately implemented sequence/translation/constraint checker. Bounded search reports incompleteness and only asserts optimality when proved.

The tool returns bounded source/target/candidate previews, status/work/pruning/score, a full JSON `artifact_id`, and an immutable `candidate` document reference when successful. The candidate and report artifact install atomically after shared session capacity checks. The candidate can be read with `read_sequence_region` and used as an explicit input in a separately target-bound cloning plan. This does not create persistent application writes or physical synthesis/verification authority.

`check_codon_design` accepts separate immutable `source` and `candidate` document references plus the same strict options. It derives the required protein again and checks supplied DNA directly; it does not optimize or trust a report embedded in metadata. It returns complete/unknown/pass/fail evidence, bounded protein/violation previews and a full JSON artifact. Source/candidate revisions and the options revision bind both tools. Artifacts paginate through `read_artifact`; their saved findings remain historical evidence and require a new current check when used again.


## Revision-bound window dot plots

`analyze_dotplot` takes `{ reference: DocumentRef, query: DocumentRef, options: DotPlotOptions, hit_start?: 0, max_hits?: 50 }`. Both document references are immutable and must match their exact revisions. The embedded Rust engine compares DNA/RNA/protein windows and reports original coordinates, forward/reverse orientation where supported, compound circular spans, compatible versus certain matches, and sampling/work/index/hit-limit evidence. It does not infer an alignment, biological homology or material ancestry.

The result exposes at most200 retained hits per page, stable hit IDs, `next_hit_start`/`done`, both source revisions, options revision, algorithm identity, lengths, coverage/completeness flags and limits. Repeat the same deterministic request with a new `hit_start` to inspect further retained hits. Full JSON and vector SVG are returned through `artifact_id` and `svg_artifact_id`; both artifacts install atomically after capacity checks. The standalone SVG embeds exact engine/options/revisions, alphabet/length/circularity fields and completion/sampling flags as unverified producer metadata. Original molecules remain unchanged. `read_artifact` pages the complete evidence; previews never claim to contain omitted hits. The session has no persistent application-write authority. [Model, matching semantics and limits](DOTPLOT.md).


## Reference-free read contig assembly

`assemble_contigs` takes `{reads: [{id, document}], options}`. Each `document` is an exact immutable `DocumentRef`; the role IDs must be distinct and their order is part of the batch binding. `ContigOptions` declares all trimming/raw-score/missing-quality policies, affine scores, minimum overlap/error/support rules, orientations and pair/total-cell/edge/column budgets. There is no external reference sequence and no installed assembler invocation.

The result includes a bounded read summary, contig and overlap previews, all counts/completion flags, a complete JSON artifact and immutable source-derived consensus document summaries. Optional `contig_start`/`max_contigs` (1–10) and `overlap_start`/`max_overlaps` (1–50) paginate retained evidence. The exact batch can reuse a report computed within this server session; `reused_result` discloses that reuse. Only internally computed results enter this cache, and exact source/options bindings are rechecked. Imported report bytes cannot populate it or gain verification authority.

The full report retains source-call/original-trace indices, raw quality provenance, selected and alternative overlap evidence, layout placements, per-column observations and consensus states. Conflicting or insufficient support remains `N`; supported gaps stay inspectable even when omitted from the emitted DNA. Reciprocal-best layout choices and bounded model completion do not prove one physical molecule, empirical correctness or independent biological observations. Duplicate retained evidence does not become extra independent support.

All report/cache/document capacity checks precede installation. A failed call cannot leave an orphan consensus document, artifact or cached result. The tool creates session-only immutable resources and grants no persistent application writes. A saved producer claim still requires fresh source-bound reconstruction before use. [Contig model](CONTIGS.md).

## Methylation context assessment

`assess_methylation` takes `{document, options, site_start?, max_sites?}`. The document is an exact immutable `DocumentRef`; options explicitly declare selected starter enzymes, Dam/Dcm/CpG source modification states and comparison/site limits. States are `Unmodified`, `FullyModified`, `Unknown` or `HemiOrMixed`. The default is unknown for all three categories. This is an assessment of recognition-site sequence contexts against versioned supplier profiles, with no physical modified-base propagation or digest-product simulation.

The result binds source, options, engine and supplier-profile revisions. It includes profile identities/source links, completion/work/limit warnings, site coordinates, category effects and overlapping modification motifs. Optional pages use `site_start` and `max_sites` (1–100); a complete immutable JSON artifact retains every reported site. Unsupported formulations or profile details, ambiguous recognition contexts, unknown or hemi/mixed states and limited searches remain explicit. An incomplete result cannot establish absence of affected sites. Profile permission is conditional on the declared model and does not prove empirical activity or grant cloning verification authority.

Sources are bounded before cloning or hashing. Invalid page ranges and artifact capacity failures leave no new resources. The tool creates no product documents, saved operations, workspace edits or inherited positive verification claims. [Scientific profile](METHYLATION.md).

## Tracked source generation and checked digest plans

`propose_plan` accepts `DeclareSourceModifications {template: {Input: id}, assumptions: {dam, dcm, cpg}}` and `CheckedDigest {template, selection, options, fragment_index}`. `selection` is the full exact `CatalogSelection`: edition/content hash, ordered enzyme identities and strand geometry. `options` is `restriction_inspection::InspectionOptions`, including comparison/site/fragment/report bounds. The positional sensitivity profile is fixed by the current core and its profile revision is recorded in evidence; callers cannot choose a different profile by passing an arbitrary revision. `list_enzymes` exposes the current profile identity. This is distinct from documentary supplier-context assessment.

`verify_plan`, HTML export and PDF export reconstruct original inputs and the independent target using the current verifier. Input-generation declarations state assumptions in the closed Dam/Dcm/CpG model, not sample measurements. `generation` on shared preparation `FragmentSpec`, whole-reaction `Row`, combinatorial `Choice` and `PreparationOverride` is a required nullable field. `null` means unknown source generation, never unmodified. A present declaration records explicit Dam/Dcm/CpG states before preparation. A present override's `generation: null` explicitly clears that declaration; only absence of the entire override inherits a choice. Saved Run/Derived sources inherit freshly replayed positional state and reject redeclaration. Standard newly synthesized PCR strands are unmodified under the model; parental modifications are not copied to them.

Checked cleavage preserves blocked sites as uncut and withholds undetermined physical products when unknown, mixed, impaired, unsupported or capped evidence cannot determine them. Reports expose positional site outcomes, profile/catalog pins, shared comparisons, limits and warnings. Product summaries retain tracked modification counts, state revision and generation provenance. Imported producer metadata grants no physical authority. Recorded scientific history is documentary; fresh native replay remains necessary.

## CRISPR guides, selected references and literal editing donors

These four sequence-bound CRISPR tools create immutable computational evidence, with no physical editing operation or application mutation authority:

- `design_crispr_guides` takes `{document, options?, candidate_start?, max_candidates?}`. Default options use canonical SpCas9; the result also exposes the canonical SpCas9/Cas12a profile DTOs. Explicit custom geometry, PAM side, GC/homopolymer filters and search limits remain bound to the report. Pages contain at most 100 candidates with exact stored-reference guide/PAM positions, orientation, ambiguity/filter evidence and nominal strand cuts. Full JSON retains every candidate. Nominal geometry does not establish empirical cleavage or guide efficiency.
- `search_crispr_off_targets` takes `{query: {id, sequence}, references: [{id, document, build, accession}], options?, hit_start?, max_hits?}`. It searches only the explicitly selected immutable references under a bounded Hamming/PAM model, returning at most 100 hits per page plus a complete JSON artifact. Build/accession/scope labels are declarations, not verified reference-build coverage. `is_intended` requires the exact reference role, revision, guide/PAM position vectors and orientation; every hit remains retained, including identical sequences elsewhere or under another role. Completion describes this selected set and model, not a verified whole-genome search or safety result. Bulges and biological activity are outside the model.
- `design_edit_template` takes `{source, target, options, guide?}`. The independent target is a noncausal requirement. Donor arms and the edited reference are reconstructed from the source and explicit edit/deletion/replacement options. An optional guide binds the exact source candidate/options and declared cut-distance condition. The full draft report is always retained when design returns; an immutable donor document is installed only if the independent sequence check is consistent. Artifact/document capacity checks precede installation, so failures cannot leave orphan donors.
- `check_edit_template` takes `{source, target, donor, options, guide?}`, with exact immutable document references for all three DNA inputs. It independently checks the supplied bounded linear donor sequence against source-derived arms/edit and the target requirement. The result preserves its exact candidate sequence hash, source/target/options bindings and all component booleans in a complete immutable JSON artifact. Imported donor metadata or a foreign positive report cannot substitute for this calculation.

All sources and aggregate reference sets are preflighted while borrowed, before owned copies or revision hashing. Pages are deterministic and source-bound. Full reports are session-only artifacts. A consistent donor check proves literal sequence agreement under these options; it does not predict cellular HDR, editing efficiency, repair outcome or organism-level safety. [Guide/off-target model](CRISPR.md); [editing-template model](EDIT_TEMPLATE.md).

## File-backed large references and cooperative genome searches

Large references use managed private scratchfiles, outside the 64 MiB small-document store. No tool accepts an arbitrary filesystem path. Normal disposal uses platform-specific owned-handle cleanup: Unix unlinks a matching managed entry while it is open, and Windows uses delete-on-close with delete sharing denied. The backend supports two active imports, four published datasets and four live searches, with a separate 64 GiB managed-byte reservation budget.

The eleven new tools are:

- `begin_genome_import {format, expected_bytes, expected_hash?, build?, options?}` accepts `Fasta` or `Cgref`. It pins the raw input length and optional BLAKE3 before returning an opaque import handle. Build/ingestion overrides apply only to FASTA; native stores preserve their exact manifest.
- `append_genome_import {import_id, offset, data_base64}` accepts at most 1 MiB decoded per call. Offsets must match exactly. Format/hash/write failures discard private partial state instead of publishing a dataset.
- `finish_genome_import {import_id, max_blocks?}` finalizes framing, then verifies at most 1–64 canonical blocks per call. Repeat while `done=false`. Publication requires every encoded block and full canonical contig hash to pass private readback; a caller-supplied progress object cannot bypass this pass.
- `inspect_genome {genome?, contig_start?, max_contigs?}` lists current dataset summaries or up to 100 contig summaries for one exact binding. Dataset handles include `id`, `manifest_hash` and `content_revision`; stale bindings fail. Build names/accessions and expected-component lists remain declarations.
- `read_genome_export {genome, offset?, max_bytes?}` reads up to 1 MiB of raw `.cgref` bytes as base64. Concatenate decoded chunks in order. No whole-file base64 artifact is constructed. Reimport performs full readback before publication.
- `read_genome_interval {genome, options}` reads at most 16,384 local bases with explicit global contig coordinates and source role. Required blocks are verified locally; no giant `Molecule` is created and no fresh whole-contig hash claim is made.
- `design_genome_interval_guides {genome, options, candidate_start?, max_candidates?}` checks bounded interval owner starts with the exact necessary right-halo/circular context. It returns up to 100 candidates plus immutable JSON evidence, retaining global coordinates, source role and exact manifest/content hashes. Local integrity is distinct from full-source coverage.
- `start_genome_off_target_search {genome, query, options?}` creates a private live whole-manifest cursor. The exact reference role, guide/PAM model, intended global coordinates and resource budgets remain bound.
- `advance_genome_off_target_search {job_id, max_blocks?}` processes 1–8 owner steps per call. The backend reads exact next owner/halo blocks itself. Clients cannot skip blocks or import a serialized checkpoint. Errors withhold complete coverage evidence.
- `read_genome_off_target_hits {job_id, offset?, limit?}` returns at most 100 retained hits with current private progress. Page completion is not search completeness. Every hit retains exact global coordinates/orientation, mismatch bounds and intended-role classification.
- `cancel_genome_job {id, genome?}` disposes an import/search or, with an exact dataset binding, releases that dataset and related searches. Owned scratchfiles are deleted; persistent application workspaces are untouched.

`source_integrity_verified` requires all canonical contig hashes; `exhaustive_for_manifest` requires every owned start under the declared model without a work/hit limit; expected-component completeness refers only to the explicitly declared component set. Even all three cannot authenticate a public assembly accession or establish biological safety. Native stores retain original-stream provenance such as `raw_source_hash`; when the original FASTA stream is unavailable, its provenance remains producer metadata rather than a newly measured input hash. The supplied upload hash, when present, checks the actual uploaded bytes.

The native tests stream 70 MiB of FASTA into a `.cgref` larger than 32 MiB, export/reimport it in bounded chunks and repeat full semantic readback without a whole-file RAM buffer. Actual SDK tests exercise private file lifecycle, cross-block guides, ordered search, native reimport, corruption, stale handles, cancellation, strict path/checkpoint rejection and unchanged small-document inventory. [Scientific coverage model](GENOME_CRISPR.md).


## Desktop restriction inspection

`inspect_restriction {document, source_role, enzymes, catalog_pin?, options?, site_start?, max_sites?, fragment_start?, max_fragments?}` inspects an exact immutable DNA document with 1–64 embedded catalog enzyme IDs or uniquely resolved names/aliases. Names resolve case-insensitively; aliases that name multiple identities fail rather than choosing one. The report pins the exact resolved geometry, ordered selection, source role/revision and computational limits. The exact catalog wrapper binds edition/hash and identities; the single tagged artifact contract below retains the complete evidence. Each site/fragment page contains at most 100 records (default 20); a complete immutable JSON artifact is available through `read_artifact`. Repeating the same request with different page offsets yields the same artifact. Invalid offsets, stale document revisions and extra authority fields fail before an artifact is installed.

The tool preserves possible recognition sites, exact origin-split spans, cuts outside substrate and bounded-search diagnostics. `scan_complete` describes recognition enumeration; `products_determined` separately describes whether the assumed-complete duplex product set exists. Neither empirical cleavage nor methylation/condition permission is inferred. Fragment summaries use stable source-cut-order indices, original strand boundaries and explicit end geometry; reference-strand lengths are a separate view. Limited scans and possible/edge evidence do not authorize releasing a product. The tool creates no sequence documents or persistent history. Browser product opening uses fresh source-event-bound worker reconstruction and native history installation.

`analyze_restriction` reports reference-strand lengths; `digest_fragments` reports assumed-complete duplex summaries from the starter enzyme model. Actual SDK tests cover the new tool's negotiation/schema, exact multi-enzyme evidence, immutable artifact pages, stale source/role changes, malformed pages, unchanged document inventory, possible/edge evidence and computational limits. Protocol tests exercise the same native/WASM dispatch kernel, current-schema rejection and integrity guards, ordinary message caps, fresh product derivation and exact source-event installation. A real executable stdio smoke also exercises the new tool; test results are recorded in the implementation status after execution.


## Exact protocol arithmetic

Three stateless tools use `collagene-protocol-units` directly: `calculate_dilution`, `balance_mixture` and `calculate_dna_molar_ratios`. These calculations neither inspect inventory nor verify chemistry, measured stock concentrations or reaction success. They do not create documents, plans or artifacts. No source DNA, network call or external program is required.

Quantity values are nonnegative **decimal strings**, at most 64 bytes and 19 fractional digits, without exponent notation. Volume units are the closed enum `nL`, `uL`, `mL`; mass units are `ng`, `ug`; molar concentration units are `nM`, `uM`. Mass concentrations specify separate `mass_unit` and `volume_unit`. Input dimensions cannot be interchanged. Exact outputs use canonical units (`nL`, `ng`, `pmol`, `ng/nL`, `pmol/nL`) and `{numerator: "…", denominator: "…"}` string pairs, avoiding floating-point loss. Arithmetic that cannot fit the bounded reduced rational returns an explicit overflow error.

`calculate_dilution` accepts `stock`, `target`, `final_volume` and optional `quantization`. Stock and target use either `{kind: "molar", value: "10", unit: "uM"}` or `{kind: "mass", value: "50", mass_unit: "ng", volume_unit: "uL"}`. Missing concentrations/volume remain missing-input errors; mixed mass/molar dimensions fail. The model is `C1*V1=C2*V2`, additive volumes and zero-concentration diluent. For example:

```json
{"stock":{"kind":"molar","value":"10","unit":"uM"},"target":{"kind":"molar","value":"100","unit":"nM"},"final_volume":{"value":"20","unit":"uL"}}
```

This needs exactly 200 nL stock and 19,800 nL diluent. Optional `quantization: {increment: {value: "0.3", unit: "uL"}, mode: "NearestTiesUp"}` rounds stock to 300 nL, refills diluent to 19,700 nL and reports 150 nM achieved concentration (`3/20000 pmol/nL`), a 100 nL absolute volume deviation and `target_reached: false`. Modes are `ExactOnly`, `Floor`, `Ceiling`, `NearestTiesUp`. No quantization is applied unless requested. `ExactOnly` rejects inexact increments; a rounded stock volume exceeding the final volume fails rather than returning a negative diluent.

`balance_mixture` accepts `final_volume` and up to 128 `{role, volume}` components. Roles must be unique nonblank identifiers of at most 128 bytes. Volumes are additive, with one remaining fill diluent. Missing component volumes retain their component index; duplicate roles and overfilled mixtures fail. This is volume balance alone; it does not infer dilution, compatibility or the diluent identity.

`calculate_dna_molar_ratios` requires `model`, `anchor`, `targets` and `require_stock_volumes`. The model explicitly supplies `model_id`, positive decimal `grams_per_mole_per_base_pair`, nonblank `assumptions` (at most 2048 bytes), and optional descriptive `source`. There is **no implicit 660 coefficient**. The anchor supplies `role`, positive integer-string `length_bp`, `mass` and optional `stock` mass concentration. Each of up to 127 targets supplies `role`, integer-string `length_bp`, positive decimal `molar_ratio` and optional `stock`. Ratios are duplex molecule ratios relative to the anchor, not end counts or recommended experimental ratios. Missing stocks yield absent volumes only when `require_stock_volumes` is false; true preserves indexed missing-stock errors.

For a synthetic declared coefficient of 650 g/mol/bp, a 3000 bp anchor weighing 50 ng is exactly `1/39 pmol`. A 1000 bp target at an explicit 3:1 molecule ratio needs 50 ng (`1/13 pmol`). With declared 25 ng/uL anchor stock and 10 ng/uL target stock, volumes are exactly 2000 nL and 5000 nL. The full caller model and assumptions remain in the returned recipe. Changing the coefficient changes molecular amounts; it cannot establish empirical mass, sequence-dependent chemistry or ligation success.

Library errors retain typed MCP error data `{domain: "protocol_units", code: "…"}`. `MissingInput` also retains `field` and zero-based `component_index` (null for unindexed dilution fields). Unsupported enum units, malformed schemas and undeclared extra fields are rejected by the MCP parameter parser. SDK tests use independent literal arithmetic, exercise exact/quantized results, dimensions, missing models/stocks, impossible dilution, overfill and overflow; the executable stdio smoke exercises all three tools too.


## Versioned restriction catalog and pinned digestion

`lookup_enzyme_catalog {query, mode?, offset?, limit?, catalog_pin?}` searches enzyme IDs/names/aliases or performs exact `Enzyme`/`Formulation` lookup. Default `Search` mode pages at most100 records (default20); exact alias ambiguity stays `Ambiguous`. A formulation is a separate vendor/product record referring to an enzyme identity. Records preserve stable IDs, geometry, alias provenance, qualifiers and source citations. Search does not collapse isoschizomers, neoschizomers or supplier formulations into interchangeable products. Every response includes actual catalog edition/content hash and a versioned `collagene://enzyme-catalog/v1/{edition}/{hash}` resource URI. That exact listed resource returns the complete sealed snapshot and documentary scope; unknown editions/hashes fail explicitly. No mutable `latest` URI or caller-uploaded catalog is accepted. Optional lookup/evidence `catalog_pin` contains exactly `{catalog_edition, catalog_content_hash}` and rejects a stale snapshot.

`resolve_enzyme_evidence {formulation, condition, context?, catalog_pin?}` uses an explicit formulation identity, typed condition and optional exact assay context. Conditions are `{kind: "incubation_temperature"}`, `{kind: "buffer_activity", buffer: "…", temperature_milli_c: 37000}`, or `{kind: "methylation", system: "dam"|"dcm"|"cp_g"}`. `Known`, `Unknown`, `Conflict` and `ContextRequired` remain distinct. Evidence keeps original contexts, qualifiers, citations and source metadata/URLs. A documented `less_than_percent` remains an inequality, and approximate relative unit-assay activity cannot become predicted digestion efficiency. Missing evidence and conflicting observations never become compatibility, an unmodified-source assumption or permission to cleave. Documentary source hashes identify the curated claims; they do not authenticate laboratory activity.

`list_enzymes` describes the starter enzyme geometry and separate pinned tracked-modification profile. Catalog resources provide the broader identity/geometry dataset. `analyze_restriction` and `digest_fragments` have the explicit starter-model scopes described above; supplier observations are not silently interpreted as permission to cleave.

`inspect_restriction` resolves real catalog identities, including enzymes beyond the starter12. An optional full `catalog_pin` is the exact `CatalogSelection` returned by an earlier inspection: schema1, edition/hash and ordered IDs with full geometry. Selection mismatch, forged geometry and stale editions fail before source cloning/computation. Inline fields and pages are bounded summaries. The sole `artifact_id` identifies a complete JSON `RestrictionInspectionEvidence`, tagged by `source`: `Original` contains `inspection: CatalogRestrictionInspection`; `SavedRun` contains `inspection: run_catalog_digest::Inspection`, including the full catalog report, freshly reconstructed strand/end/modification evidence and origin findings. Each artifact publishes as one bounded atomic reservation and deduplicates across pages. No separate raw report or second artifact ID is produced. Neither a catalog pin nor documentary methylation/buffer facts upgrade assumed-complete products to experimentally verified cleavage.

Fresh catalog fragment products retain exact source-event, catalog selection, engine and product-revision producer fields. GenBank COMMENT and EMBL CC export preserve them in bounded physical-producer marker chunks. Reimport stores only `source:restriction_catalog:*` claims and marks chemistry unmapped; it never promotes comments into native fresh-product keys. Invalid/truncated claims stay disclosed and cannot restore physical or scientific authority. FASTA remains sequence-only. Current source-event ancestry and a fresh embedded calculation are required for product installation.

Protocol fixtures cover current request/result selectors, exact selection and source binding, rejection of noncurrent encodings, ordinary 32 MiB versus installation-only 128 MiB budgets and fresh native event ancestry. Format fixtures cover both annotated text formats, exact claims, malformed/truncated producer data and nonrestoration of chemistry. Actual SDK fixtures cover versioned resource discovery, explicit starter scope, contextual/unknown evidence with citations, expanded geometry, complete tagged artifacts and atomic single-item capacity rejection. A separate adapter fixture preserves conflicting/context-dependent evidence without modifying the installed runtime catalog.

## Reference-only BP/LR exchange

`inspect_recombination`, `check_recombination_target` and `select_recombination_component` require two ordered immutable source roles and a mandatory pinned full-site profile/reaction/limits request. Selection replays source-derived reference DNA and independently checks the specifically selected component before atomic session installation. All reciprocal components and uncertainty remain evidence; no biochemical activity or experimental winner is inferred. Explicit profile discovery is available through `collagene://recombination/profile/manufacturer-context-v1`. See [contracts, privacy and limits](RECOMBINATION_MCP.md). The typed operation remains reference-only; imported graphs cannot restore material authority.


### Ordered overlap assembly recipes

`propose_overlap_assembly` compiles an explicit practical recipe into the existing typed candidate plan. Supply `inputs: [{id, document: {id, revision}}]`, an independently imported `target`, and `recipe: {name, fragments, overlaps, circular, options}`. Each ordered fragment declares a distinct `role`, an original `input_id`, `reverse`, and one preparation:

- `"WholeLinear"`: use the actual complete linear DNA source. A circular source cannot silently become a linear fragment.
- `{"Digest":{"enzymes":["EcoRV"],"fragment_index":1,"repair":null}}`: reconstruct complete digestion and select the explicit zero-based fragment in existing top-cut order. Original uncut source remains a plan input. Explicit `repair: null` preserves the actual cut ends; optional `FillFivePrimeOverhangs`, `RemoveThreePrimeOverhangs` or `Blunt` declares an additional embedded end-repair operation rather than silently converting cohesive ends.
- `{"ExactPcr":{"forward":Primer,"reverse":Primer}}`: reconstruct PCR using actual full ordered oligos represented by their explicit `name`, 3′ `annealing` sequence and 5′ `tail`. No tails or primers are invented.

`overlaps` declares exact lengths between adjacent prepared, oriented fragments (including the closing junction for a circle). `options` contains the existing explicit Gibson/HiFi/InFusion advisory method profile and bounded repeat-search controls. Advice is distinct from exact homology geometry, physical preparation and independent target agreement. Shared source inputs can deliberately feed multiple named preparation roles; this does not select a unique experimental assembly outcome.

The response retains `state: "candidate"` and `verified: false`, plus a fresh deterministic `verification` summary, `preparation_complete`, context and plan revisions, and an immutable `overlap-assembly-recipe.json` artifact. The full artifact contains preparation rows, original inputs, fresh prepared products, junction/profile evidence and complete verifier findings. Partial preparation, invalid geometry and wrong targets remain inspectable drafts; they do not become a proof. Candidate and evidence publication is atomic under existing 128-item and 64 MiB session limits; report JSON is limited to 32 MiB. Sources and target are bounded before copying, and sources are never modified or automatically replaced by products.

Use the returned `plan_revision` with existing `verify_plan`, `export_candidate`, `export_protocol` and `export_protocol_pdf`; verification and protocol export reconstruct again. Saved artifacts and source annotations never restore scientific authority. Successful exact agreement is conditional on declared preparation assumptions, not laboratory success or an exhaustive reaction search.

The source-authored SDK fixture specifies a literal PCR left fragment and EcoRV-prepared right fragment sharing 15 bases, with an independently supplied target. It preserves strict equality checks and covers wrong targets, explicit reverse orientation, stale sources, unchanged original documents, candidate reuse and HTML/PDF export. The stdio smoke also exercises a WholeLinear recipe through the actual tool. Execution results are recorded separately; adding these fixtures does not assert that they have run.


### Saved primer reuse

`resolve_primer` resolves one `primer_bind` feature from an exact immutable `{document: {id, revision}, feature_index}`. It returns the actual `resolved.primer` (name, annealing domain and explicit 5′ tail), current source revision, feature index, sequence origin and warnings through the same strict Rust resolver as the web primer controls. Repeated/conflicting sequence qualifiers, inconsistent ordered/domain splits, unsupported missing-field geometry and stale sources reject. An explicitly stored mismatch remains a literal oligo with a warning; it is not corrected against the template. Generic missing-field inference requires supported exact strand-oriented geometry. The 200-base oligo bound and IUPAC ambiguity are preserved.

Resolution creates no document or proof. Explicitly pass the returned `primer` into `simulate_pcr` or `analyze_primer` with the actual current template to reconstruct binding/products independently. Saved locations and historical analysis do not establish amplification, specificity, chemical preparation or thermodynamic conditions. See [primer reuse](PRIMERS.md#reuse-from-mcp).


### Ordered ligation recipes

`propose_ligation_assembly` accepts `{inputs, target, recipe}` with the same exact immutable source document references and independent target as overlap recipes. Its recipe is `{name, fragments, circular, allow_nicked}`. Each fragment is `{source: FragmentSpec, treatment}`: `source` declares actual `role`, `input_id`, `preparation` (WholeLinear, Digest with optional explicit end repair, or exact PCR) and `reverse`, plus required nullable `generation`. `treatment` is explicitly `null`, `"Phosphorylate"` or `"Dephosphorylate"`. The compiler emits each actual preparation, optional repair and terminal treatment before its final ordered ligation step. It does not infer phosphorylation, nick repair or physical duplex state from imported metadata.

Fresh verifier reconstruction establishes actual oriented strands, 5′/3′ overhangs and terminal chemistry before checking target sequence/topology. `allow_nicked: true` permits explicitly retained nicks; it never seals or hides them. Unknown or blocked termini, incompatible ends, omitted required treatment and mismatching targets remain explicit failed drafts with their available preparation evidence. The intended order remains explicit. The ligation-outcome evidence also retains alternative compatible one-copy paths under recorded work/product/enumeration limits; completion within that scope is not exhaustive real-mixture behavior or evidence of isolation.

Atomic candidate/evidence publication uses the existing session caps and 32 MiB JSON artifact limit. Returned `plan_revision` works with `verify_plan`, binary candidate export and HTML/PDF protocols. `verified: false` records candidate status; the separate fresh `verification` summary and `product_physical` summary show conditional results. Complete preparation rows, physical fragment summaries, each junction (including circular closure), product nick/end evidence and verifier findings remain in `ligation-assembly-recipe.json`. Full source documents remain unchanged in the plan; imported annotations and saved artifact claims never grant chemistry authority.

Source-authored SDK fixtures independently specify two actual PCR products with explicit phosphorylation, one omitted-treatment nick case, wrong targets, reversal and stale-source rejection, plus three separately cut circular EcoRI sources whose compatible AATT cohesive ends generate a literal 42-base circle. The stdio smoke includes the three-piece case. These are fixture definitions, not a claim of execution or laboratory validation.


Dedicated overlap/ligation recipe compilers have separate bounded generation
profiles: up to 128 ordered source roles, with at most 257 overlap-preparation
steps or 385 ligation-preparation/treatment steps. They still use existing typed
operations, fresh verification, bounded candidate storage, 32 MiB evidence
artifacts and the shared 64 MiB managed session budget. Generated plans do not
bypass those storage/resource checks or establish experimental success.

### Assembly junction inspection

`inspect_assembly_junctions` accepts an exact saved `plan_revision` and explicit bounded flank/genetic-code `options`. It reuses that plan's actual immutable source inputs and separately stored independent target; a caller cannot substitute a newly claimed product or target during inspection. Fresh verifier reconstruction remains part of inspection. Supported final ordered overlap/ligation methods expose actual reference-coordinate junction DNA and available CDS/frame context. Unsupported final operations are reported explicitly rather than converted into guessed junctions.

The response is a bounded summary with a complete immutable JSON artifact. DNA homology, boundary-crossing codons and CDS annotations describe the declared reconstruction model and current annotation evidence. They do not establish expression, fusion-protein function, a validated biological reading frame, experimental success or an exhaustive assembly search. Source annotations retain uncertainty/review status. A wrong independent target remains a failed verification even if locally reconstructed junction DNA can still be inspected.

Inspection changes no source, candidate or verification authority. Existing `verify_plan`, HTML and PDF protocol exports remain the way to recompute and review the saved candidate. Complete evidence is downloadable through `read_artifact` in bounded chunks; imported artifact contents cannot recreate scientific authority.


Junction options explicitly specify `flank_bases` (at most 150 per side), `genetic_code` (1, 2 or 11), `frames: [{reverse, frame}]` (frame 0, 1 or 2 in each selected orientation), `max_junctions` and `include_cds`. Frame offsets are relative to the displayed window or its reverse complement, not a silently inferred full-product CDS frame. Boundary-crossing codons expose the oriented-window offset; source-CDS context is separately conditional on exact location, strand, phase and review status. Circular windows retain their ordered origin-crossing coordinate parts, while shared overlap homology remains an explicit interval, not duplicated DNA.

Inline previews retain at most ten junctions, six local frames and ten CDS descriptions per junction, with explicit omissions/truncation; full translations, original verifier findings and physical evidence remain in `assembly-junctions.json` (32 MiB cap). A physical junction's nick count refers to the complete fresh product nick state; it is not an inferred bond-by-bond seal proof. Literal source-authored SDK fixtures cover a GAA codon crossing an actual phosphorylated-PCR ligation boundary, wrong target, reversed preparation, shared overlap homology and cohesive circular-origin closure. Adding these fixtures does not assert their execution.


### Explicit original-source terminal assumptions

The typed plan and overlap/ligation recipe tools expose `Operation::DeclareSourceTermini {template, left, right}` and `Preparation::DeclaredDuplex {left, right}`. Each terminus must explicitly declare `five_prime` (`Phosphate`, `Hydroxyl` or `Unknown`) and `three_prime` (`Hydroxyl`, `Blocked` or `Unknown`); there is no default terminal state. In a recipe, for example, `{"DeclaredDuplex":{"left":{"five_prime":"Phosphate","three_prime":"Hydroxyl"},"right":{"five_prime":"Phosphate","three_prime":"Hydroxyl"}}}` states an operator hypothesis about an original intact linear duplex input. It does not report measured chemistry, perform synthesis or infer state from source annotations.

The operation accepts only an original `Source::Input`. Circular, ambiguous, single/mixed-stranded, source-step, imported physical, tracked, variant or reference-only states cannot be promoted into intact duplexes. Source sequence, features and citations are retained unchanged; methylation remains unknown. Existing physical/guard metadata cannot be overridden by a declaration, and its presence never grants authority. Descendants require their originating graph. This is not nick repair, end repair, denaturation, chemical recovery or a generic materialization bypass.

An ordinary `WholeLinear` input still has unknown chemistry. Omitting the declaration does not silently enable ligation; unknown 3′ chemistry remains incomplete and blocked termini cannot form the missing bond under a sealed-junction policy. Successful ligation is conditional on the explicit source/end assumptions and fresh independent-target verification. New producer summaries are restrictive evidence, not a future import capability. HTML/PDF exports retain the exact declaration and current verifier revision.

Source-authored SDK fixtures ligate imported linear AAAA and CCCC under declared 5′-phosphate/3′-OH assumptions against an independent AAAACCCC target; they check unchanged original native exports, direct typed-plan parity, omitted/unknown/blocked ends and refusal to promote prior steps, ssDNA, circular and imported physical/reference-only sources. The stdio smoke includes a declared-source recipe plus the omitted-declaration incomplete case. Execution is recorded separately; these descriptions do not imply fixtures have run.


### Practical overlap-primer planning

`design_overlap_primers` accepts `{inputs, target, request}`. Inputs are actual immutable document bindings (`id` plus `document: {id, revision}`); the target is independently supplied. Request fields are `name`, ordered `fragments: [{role, input_id, start, end, reverse}]`, one `overlap` length, `circular`, explicit `primer_options` and `profile_options`. Coordinates are zero-based half-open source intervals, with supported circular intervals described by the core profile. Up to 32 fragment roles and overlap lengths 1–80 are supported under bounded native-source preflight. Pre-existing primer tails must be empty so synthetic homology is not silently combined with unrelated sequence.

The existing anchored endpoint heuristic selects a pair for each actual template interval under explicit length, GC, structure and basic Tm/solution conditions. One-sided homology comes from the next **source** interval in its requested orientation, never from target bases. Actual annealing domains and generated 5′ tails remain separate; complete ordered oligos, orientation and original input IDs are returned with a practical `assembly_recipe::Recipe`. This is a reproducible heuristic, not Primer3, a genome-wide specificity engine or an empirical success prediction. Reported Tm is not automatically a PCR protocol temperature.

The result explicitly has `state: "advisory"`, `verified: false`, a context revision, ordered source revisions, `recipe`, bounded design summaries, order-ready `ordered_oligos`, warnings and a complete `overlap-primer-design.json` artifact (32 MiB). The artifact retains exact source and independent target references alongside the typed report. No source, candidate or verification state is installed. A wrong target does not choose different source bases or make an advisory design valid; it remains independently bound and must fail subsequent verification when mismatched.

Pass the returned `recipe` with the same actual `inputs` and independent `target` to `propose_overlap_assembly`. That separate call reconstructs exact PCR preparations and homology, then runs the verifier; it returns a candidate and fresh conditional status for normal candidate/HTML/PDF workflows. Unsupported geometry, stale source revisions, pre-existing tails and exceeded limits reject explicitly. Source-authored SDK fixtures specify actual original 20-base templates and a literal 40-base target, assert the hand-derived next-source homology tail, feed the recipe through real preparation/verification, retain mismatched-target designs as unverified, and exercise invalid bindings/intervals/tails/limits. The stdio smoke exercises the same planner-to-compiler seam; execution is recorded separately.


### Assembly intent before molecular preparation

`build_assembly_intent` accepts `{inputs, request}` and **no target argument**. Bind actual immutable sources with `inputs: [{id, document: {id, revision}}]`. Request is `{name, fragments: [{role, input_id, start, end, reverse}], circular}`: exact source regions in explicit order and orientation define intended reference DNA. Coordinates are zero-based half-open; supported circular source intervals remain limited to one source turn. Region selection, orientation and intended topology are design decisions, not physical preparation operations.

The core constructs this intended reference without running PCR, digestion, ligation or overlap simulation. The returned `state: "intended-reference"`, `verified: false`, `preparation_simulated: false` and exact `target: {id, revision}` make that separation explicit. Source bindings, context revision, every source/output placement and a full `assembly-intent.json` artifact are retained. The target carries intent provenance and a restrictive reference-only physical-state marker; it is not measured DNA, prepared duplex material or an inherited chemical state. Building intent changes no original source snapshot and creates no candidate plan.

Report and target installation are atomic under existing managed-session limits. Borrowed source guards precede copying/hashing; full report JSON is bounded to 32 MiB, with the same shared 64 MiB session accounting and document/artifact capacities. Repeated identical requests deduplicate. The exact target reference can be passed directly to `design_overlap_primers` and `propose_overlap_assembly` with the original inputs. Those separate steps design actual oligos, reconstruct preparations and verify their output against the already frozen intent. Expected DNA is never obtained by copying the verifier's product. A separately imported desired target remains an alternative input route for reverse design and AI checking.

The source-authored SDK fixture constructs a literal 40-base intended reference from two actual 20-base sources before invoking the primer planner or simulator, then checks real PCR/overlap agreement. It also covers reversed regions and exact output placement, circular intent, unchanged original exports, idempotence, stale/invalid/missing bindings and explicit wrong-source rebinding followed by failed original preparation. The stdio smoke now follows source selection → intended target → oligo design → separate assembly verification, with no complete-target import prerequisite. Execution evidence is recorded separately.

## Catalog-pinned cloning preparation

Existing `propose_overlap_assembly`, `propose_ligation_assembly` and `propose_plan`
accept the catalog-pinned digest preparation/operation. The shared preparation is
`CatalogDigest { selection, fragment_index, repair }`, retaining the exact
`CatalogSelection` from embedded-catalog inspection and a zero-based fragment.
The compiler emits `Operation::CatalogDigest` and optional explicit `RepairEnds`.
The original immutable document is the input; a downloaded fragment's producer
metadata does not restore a physical substrate. Validate the actual plan against
its independent intended target and export HTML/PDF from that same plan.

Fresh evaluation checks the catalog edition/content/IDs and exact geometry before
using the existing tracked digest engine. The catalog pin does not predict
methylation compatibility, reaction conditions or experimental completion.
Native plans reconstruct the selected pinned catalog geometry and explicit processing against the separate target.

### Order sheets for a reviewed plan

`export_plan_oligos` accepts `{plan_revision, format: "Csv"|"Tsv", options: {deduplicate: true}}` (options default to deduplication). It uses the exact saved plan and stored independent target, extracts actual declared PCR/mismatch-PCR oligos and eligible original anneal inputs, and atomically publishes spreadsheet-safe order text plus full JSON use evidence. It does not run verification or design primers. The returned plan/target/context revisions bind the output; `verified: false` explicitly avoids a cloning claim.

Exact ordered DNA/modification dedup preserves separate names, sources, annealing/tail partitions and mismatch policies in `uses`. Same-name different-sequence oligos retain distinct IDs. A differing partition is not silently represented by one misleading domain column. No quantities or supplier modifications are invented; downstream phosphorylation is not an order request. Original anneal termini are per-use state evidence only. Opaque/unrepresentable original order data sets `complete: false` and `unavailable`; partial sheets are review evidence and must not be treated as all-required orders. Step-derived anneal material does not imply a new synthesized oligo order.

The response includes counts, at most 100 uses and explicit truncation, report/export warnings, `table_artifact_id`, `report_artifact_id` and filename. Full `plan-oligos.json` retains all uses. Existing 128-artifact/64 MiB session accounting and paired-output atomic checks apply; each core order receipt is bounded to 32 MiB. Complete order data still does not prove that the plan amplifies or matches its target. Use `verify_plan` and existing HTML/PDF exports independently. See [plan oligo output](PLAN_OLIGOS.md) for profile, bounds and current unsupported cases.


### Finalize reviewed cloning outputs

`finalize_cloning_outputs` accepts an immutable uploaded native `workspace: {id, revision}`, an exact saved `plan_revision`, explicit `target_binding: {document_id, event_id, revision}`, and `input_bindings: [{input_id, event: {document_id, event_id, revision}}]`. The bindings refer to actual retained workspace events; equal DNA, names or imported claims never choose them automatically. Options default to `{name: "Assembly", show_prepared: false, save_primers: false, allow_incomplete: false}`. A target binding differing from the saved independent target rejects.

The core reconstructs through the existing retained-run path and creates a **new immutable session workspace**. It retains every actual step product and source/target relationship. `show_prepared` changes presentation visibility only, never historical material retention. The final construct receives its chosen current name through a new history event (or a separate descendant when an edited final head must be preserved); the original recorded event and original sources are not rewritten. `save_primers` retains the existing `oligo-table-v1` attachment only when faithful complete nonempty plan order data is available. An opaque/unavailable required order, failed final step, incompatible ends or invalid bindings cannot partially publish outputs. Incomplete accepted evidence requires explicit opt-in and remains labeled as such.

The response contains a new `workspace` reference, `receipt`, full commit JSON `artifact_id`, and `native_artifact_id` for the binary `.cg2workspace` archive (`encoding: "base64"`). Receipt includes context and original-workspace fingerprint, the retained-run receipt, exact construct event, optional primer attachment hash, and visible/hidden output document IDs. The workspace, revision pool, both artifacts and capacity reservation publish atomically. Identical retry context/options deduplicate without appending events or artifacts. Existing session limits apply, including full commit JSON at most 32 MiB, native archive limits and the shared 64 MiB managed session budget.

Read binary artifacts in chunks, concatenate base64 and decode before saving. Downloading or explicitly importing that archive is separate from finalizing: this tool does not edit a browser workspace, hosted service or arbitrary filesystem path. It does not establish experimental material, laboratory success or inherited positive verification. Reopen/replay with the retained exact sources to reconstruct current evidence again; hiding prepared rows does not remove their replay history.

Source-authored SDK coverage uses real chunked native upload, explicit source-event mapping, actual PCR/overlap preparation, saved oligos, renamed final reference, complete native reload and replay, unchanged original archive and identical retries. Negative cases cover wrong/stale event bindings, physically unsealed preparation and opaque order data with no new workspace publication. A focused private capacity fixture stages two artifacts with only one slot remaining and checks that no workspace/revision/artifact is partially installed. The optional stdio `--native-workspace` flow exercises finalization and replay against its explicit historical fixture; ordinary smoke still checks tool discovery. Execution evidence is recorded separately.

## Type IIS primer proposals

`design_golden_gate_primers {inputs, target, request}` accepts exact immutable
source `BoundInput` records and an independently supplied target. Request fields
are `name`, `enzyme` (`BsaI`, `BsmBI`, `BbsI`, `SapI`), concrete `clamp` and `spacer`,
and `fragments: [{role, input_id, start, end, reverse, forward_length,
reverse_length, right_overhang}]`. Coordinates are zero-based half-open original
intervals; circular intervals may cross the origin by at most one turn. The right
fusion is explicit downstream assembly-top DNA. Spacer and fusion lengths follow
the actual enzyme geometry, including SapI's three-base fusion.

The response is advisory (`verified:false`) with context revision, requests,
original-source `PcrCatalogDigest` fragment preparations, bounded row previews and
an immutable full JSON artifact. Each oligo separates clamp/motif/spacer/fusion
from its annealing domain. Retained-top previews are limited to160 calls with
explicit truncation; full sequence stays in the artifact. Sources and candidates
are unchanged. Reports are bounded to32MiB and existing session limits apply.

Wrap each returned fragment as `{source: fragment, treatment: null}` in an ordered
`propose_ligation_assembly` recipe with the same original inputs and independent
target. Existing-site catalog digest fragments may be retained alongside proposed
PCR rows. That tool freshly reconstructs PCR→digest→optional repair/treatment and
ligation using actual chemistry; normal verification and HTML/PDF/order outputs
reuse the resulting plan. Changing the target changes the context token but never
selects different oligo bases. No source is domesticated to remove internal sites,
no external tool is invoked, and ordered consistency is not experimental yield or
an exhaustive mixture claim. [Detailed profile](GOLDEN_GATE_PRIMER_DESIGN.md).

## Mutation intent before preparation

`build_mutation_intent {document, request}` accepts an immutable original reference
and `{name, edit_start, edit_end, inserted}`. Coordinates are zero-based half-open
boundaries; equal boundaries insert, empty inserted DNA deletes, and circular edits
can cross the origin by at most one turn. A circular terminal insertion caret at
`start=end=length` is valid and preserves the original origin.

The response is `state:"intended-reference"`, `verified:false` and
`preparation_simulated:false`, with an atomically retained target `DocumentRef`,
full immutable JSON artifact, context/source revisions, deleted length, explicit
origin shift and warnings. Source DNA is unchanged. No PCR simulation generates
the expected DNA. Pass the original document and returned target to
`design_inverse_pcr`, then verify/export the resulting plan separately. Intended
references remain restricted against reuse as physical inputs. Borrowed guards
precede source cloning; source DNA is bounded to8MiB, inserted DNA to100,000 bases
and report output to32MiB, with existing session capacities. Unsupported coordinate
qualifiers fail before publication rather than receive guessed mappings.
[Worked flow and interpretation](MCP_MUTATION_INTENT_FLOW.md).

## Construct-owned read check capture and replay

`record_read_check {workspace, reference, reads, options}` takes an immutable native
`WorkspaceRef`, exact `EventBinding {document_id,event_id,revision}` reference and
ordered1–32 read bindings, plus existing `TraceValidationOptions`. It freshly checks
retained events and atomically returns a new immutable workspace, receipt, full
saved-check JSON artifact and native archive artifact. All artifact/workspace slots
and bytes are reserved before publication. The original workspace, source heads and
raw reads are unchanged. Repeated identical evidence reuses its attachment; event
aliases remain distinct. Recorded status and engine are historical receipt fields.

`list_read_checks {workspace}` returns historical summaries by attachment hash,
construct reference event, read count and recorded engine/status. It never replays.
`replay_read_check {workspace,attachment_hash}` resolves the exact saved original
events and reconstructs with the current engine/model. It returns fresh computational
status/coverage and a full immutable JSON artifact; source edits and navigation are
not rewound. Neither capture nor replay establishes physical sample identity,
experimental success or unobserved bases. Partial read support remains incomplete.
Existing quality, ambiguity, placement and duplicate-file policies are retained.

Saved evidence is bounded to17MiB and128 records per workspace; existing native
workspace/session budgets remain. Replay artifacts use the bounded32MiB scientific
JSON cap. No external aligner or basecaller is invoked.
[Native read checks](READ_CHECKS.md).

## Explicit read review decisions

`preview_read_review {workspace,request}` and `commit_read_review {workspace,request}`
accept `{attachment_hash,name,decisions:[{key,choice,reason}]}`. Keys are
`Site {position}` (zero-based reference base) or `Insertion {index}` in the fresh
report. Choices are `RetainReference`, `Unresolved`, `UseBase {read_index,call_index,
base}`, `DeleteBase {read_index}`, or `UseInsertion {read_index,boundary,sequence}`.
The nominated call index is `ReadColumn.call_index` in the bound source event;
original stored-call/sample mappings remain in the report. Every changed allele,
location and eligible observation must match fresh replay exactly. No physical
read, confidence value or source reference is rewritten.

Preview returns a bounded proposal summary and full immutable JSON artifact, with
no document installation. Commit freshly reconstructs and atomically returns a new
immutable workspace, historical receipt, full review JSON and native archive.
Explicit retain/unresolved-only requests produce no target; changed references
retain all original source events as computational parents and remain restricted
proposals requiring a new check. Conflicting observations and unreviewed positions
are never relabelled as a consensus or experimental confirmation.

`list_read_reviews {workspace}` returns historical decision/change counts, original
reference bindings and optional derived-reference bindings. Reopening or importing
these records does not restore current verification. Limits are256 decisions,
1024-byte reasons,100,000 replacement bases and128 native review records under the
existing model/workspace/session budgets. Both output artifacts and workspace slots
are reserved before publication. [Read review scope](READ_REVIEW.md).

## Junction-edited assembly intent and original-source primers

`build_edited_assembly_intent {inputs,base,edits}` takes exact immutable `BoundInput`
records and the existing base assembly-intent request. Each edit is
`{junction,remove_left,remove_right,inserted}`: trim the upstream oriented core's
end, trim the downstream oriented core's start, then place literal DNA between them.
Junctions are zero-based and circular assemblies include closing junction n−1.
Sources remain unchanged. The tool atomically retains an unverified intended target
`DocumentRef` and full `{inputs,intent}` JSON artifact before any reaction simulation.

`design_edited_assembly_primers {inputs,request}` takes
`{base,edits,overlap,primer_options,profile_options}` and returns advisory `recipe`,
`ordered_oligos`, design count, context/intent/target revisions and full
`{inputs,report}` artifact. It does not choose a target from a simulated product or
silently install oligos. Verify that `target_revision` matches the independently
built target, then use `propose_overlap_assembly` with that target, same originals
and returned recipe. Primer analysis includes actual complete synthetic tails;
existing filters/conditions and scientific uncertainty remain explicit.

This route proposes every row as exact PCR followed by profiled overlap assembly;
it is separate from selecting an existing digest fragment or Golden Gate preparation.
Reverse rows use original source intervals and the appropriate original primer tail.
Limits are32 rows/edits,100 normalized bases per insert,4096 inserted bases overall,
existing200-base ordered oligos, bounded source/annotation/design work and32MiB
artifacts. Invalid trims, stale inputs, ambiguous concrete proposals and failed
filters do not mutate originals or promote candidates. [Model](ASSEMBLY_JUNCTION_EDIT.md).


## Selective junction primer proposals

`design_selective_assembly_primers {inputs,request}` is an advisory source-bound
planner. Ordered rows explicitly retain a complete preparation or request new PCR
primers; the report records aligned optional designs and retained row indices.
Only newly proposed oligos appear in `ordered_oligos`. Full immutable JSON evidence
retains both source bindings and the entire report. Existing manual primers,
catalog pins, fragment choices and repairs are not silently replaced. The target
revision agrees with separately built `build_edited_assembly_intent` for the same
sources/base/edits. Fresh `propose_overlap_assembly` against that independent target
must still check all preparations and actual products. See
[selective SDK contract](SELECTIVE_JUNCTION_MCP.md).


## Local assembly diagnostics

Existing `propose_overlap_assembly` and `propose_ligation_assembly` responses and
full JSON artifacts add `diagnostics`, bound to the same preview/plan/target.
Rows distinguish input references, reconstructed preparation chains, failed chains
and blocked work. The first missing operation is located from the recorded dense
product prefix and compiled row ownership, never by parsing findings. Target
mismatch after all operations remains a separate stage. Exact homology/profile
locations can identify overlap problems; available ligation endpoints do not by
themselves identify an incompatible junction or establish chemistry. Full original
findings, prepared products and verification remain available. This adds no tool,
new authority or protocol version. [Detailed contract](ASSEMBLY_DIAGNOSTICS_MCP.md).


## Golden Gate junction primer proposals

`design_golden_gate_junction_primers {inputs,request}` proposes original-source
PCR→catalog-digest preparations only for selected rows of an independently
junction-edited construct. Retained names/domains/tails, source direction, catalog
pins, fragment index, repair and terminal treatment stay exact. The full bounded
artifact retains source bindings and the entire report; `ordered_oligos` includes
only selected replacements. The returned intent/target revision agrees with
separately building `build_edited_assembly_intent` for those sources/base/edits.
Pass the complete recipe and that target to `propose_ligation_assembly` for fresh
preparation and final verification. This does not simulate exhaustive reaction
cycling, optimize fusion choices or domesticate internal sites. See
[SDK contract](GOLDEN_GATE_JUNCTION_MCP.md).


## Prepared Golden Gate alternatives

`OrderedTypeIis {enzyme, allow_nicked}` is the combinatorial ordered Type IIS restriction/ligation method. It uses positional checked preparative cleavage, specified fragment order and independent-target verification; it does not model cycling or enumerate a reaction mixture. `GoldenGatePool` is the distinct finite-pool method with explicit substrates, candidate paths and recut assessment. There is no `Method::GoldenGate` alias for the ordered method.

Both `inspect_prepared_golden_gate_pool` and `inspect_golden_gate_reaction` require `origin_workspaces` in their request, including an explicit empty list when no saved origin is used. Every saved-source role must name its exact retained immutable workspace/event proof. Source and target document references, native material origins, catalog/profile pins and generation declarations are checked before reconstruction. Whole-source substrates use `WholeSource {origin}` with required nullable origin: `null` is an original source, otherwise the saved Run/Derived material is freshly replayed without flattening it to sequence.

All preparative cuts, source-pool cleavage and candidate recut checks share one comparison budget. `source_cleavage`, `profile_revision`, `products_determined` and `comparisons` disclose the source assessment. If any requested source has undetermined products, the entire pool is withheld rather than enumerating an available subset. Each candidate retains its checked recut evidence and each recognition site its permitted/blocked/impaired/unknown outcome. Candidate uncertainty is not a site-free, isolated-product or successful-reaction claim.

`inspect_prepared_golden_gate_pool {inputs,target,request,origin_workspaces}` uses fresh actual
preparation strands with the existing bounded Golden Gate pool model. It reports
all compact candidate summaries, target matches, subset/source roles, recutting
and exact search completeness; the full immutable JSON artifact retains source
bindings, independent target, ordered preparation evidence and every pool result.
Only selected prepared roles are included. Unrepresented discarded digest pieces
and original-source carryover are outside that declared pool. Missing material
withholds enumeration instead of silently narrowing the pool.

`select_prepared_golden_gate_candidate {artifact_id,selection}` requires exact
parent context, candidate index and both keys. It atomically stores an unverified
native plan and binary `candidate.cg2plan` artifact, retaining every preparation
step and original input. Only the final operation changes to existing GoldenGate
with the recorded search options/index. Its `artifact_id` and
`candidate_artifact_id` identify the same binary bytes; `export_candidate` and
`verify_plan` remain available by the new plan revision. Selection and reported
target correspondence do not grant a positive cloning verdict. See
[MCP contracts](PREPARED_GOLDEN_GATE_MCP.md).


## Golden Gate path junction inspection

The existing `inspect_assembly_junctions {plan_revision,options}` tool now handles
selected GoldenGate plans. Its `method.GoldenGate` data identifies the exact
native candidate and selected digest path, including original Input/Step revisions,
raw source cuts, actual orientation and output spans. Boundaries, physical ends,
flank windows and local codons come from that selected product, not target DNA or
recipe order. The full JSON artifact retains all parts and all reported windows;
existing compact-response truncation flags still apply.

Target mismatch does not hide actual reconstructed windows; failed final assembly
does not substitute a preparation. Finite-pool enumeration status remains separate
from geometry/window completeness. Inspection does not promote a candidate or
claim a functional fusion. Large Golden Gate reconstruction runs outside the
session store lock after borrowed source/target preflight. The current typed report binds the selected native path.


## Declared whole-reaction Golden Gate

`inspect_golden_gate_reaction` includes complete original circles/linear DNA,
complete PCR amplicons or explicitly isolated preparations according to each row.
It retains all digest pieces, separate uncut-circle carryover, backgrounds and
finite-search limits. `probe_plan_revision`/`probe_verification` are diagnostics,
not an automatically saved winner. `select_golden_gate_reaction_candidate` binds
an exact artifact/context/index/key selection to an atomically stored unverified
native plan and binary candidate; fresh `verify_plan` remains mandatory. Both
operations preserve independent targets and explicit material choices. See
[reaction MCP contract](GOLDEN_GATE_REACTION_MCP.md).

## Named digest lanes in predicted gels

`predict_gel` supports catalog-pinned lanes via
`preparation: {"CatalogDigest":{"recipe":...}}`. Such a lane also requires
`workspace: {"id":...,"revision":...}` naming an imported native workspace.
The complete recipe includes the retained source event, catalog edition/hash/IDs
and geometry, inspection limits and explicit complete-cleavage/methylation-unassessed
mode. The source is resolved historically without moving heads and must equal the
lane's immutable document. The tool rejects missing or stale workspace proofs,
wrong event owners/revisions and metadata/content aliases, then freshly inspects
source DNA. The exported analysis JSON retains all declared lane recipes/references.
Uncut, starter-digest and PCR lanes use their explicitly typed preparations.

Saved setup capture/reopen in the browser remains historical until explicitly
recomputed. See [named digest semantics](NAMED_DIGESTS.md).


## Published profiles and source-backed TOPO declarations

`list_ligation_profiles` returns five exact pinned experimental profiles with
formulation, conditions and separately licensed data provenance.
`score_published_ligation_profile` accepts `profile`, intended `ends` and optional
`physical_competitors` (literal5′→3′ singleton ends). `rank_golden_gate_fusions`
accepts original `inputs` plus `{design,profile,enzyme_formulation}`. It returns
advisory rows and a full artifact; use the original selection context and explicit
candidate ID with `materialize_golden_gate_fusions`, then freshly verify the recipe.

`score_ligation_fidelity` accepts optional `physical_competitors`; omission preserves
the old directed-count report, presence uses the separately versioned competition
report. Free physical ends are not scored as intended joins. No combined estimate
is produced when there are no joins. [Model contract](PUBLISHED_LIGATION_PROFILES.md).

`compile_topo_preparation` accepts exact historical `sources` and a shared typed
TA/directional-CACC `declaration`. It returns the source-backed native request,
strand indexes and direction flag for existing `inspect_topo_orientations` /
projection/check tools. Declared activation is an assumption; neither compilation
nor a transferred source-call graph implies a repaired or closed clone.


## Ordered physical TOPO reactions and saved graphs

`simulate_topo` reconstructs exact original sources using a single ordered
`request.steps` transcript. `Transfer` names actual source `Call` endpoints;
`Process` declares flap cleavage or selected nick sealing with explicit chemistry
and profile assumptions. Failed steps and every source call/component remain
visible. No report or projected sequence is accepted as native state.

`project_topo_actions` uses `request: {preparation, selection, options}`. Selection
pins `source_batch_revision` from the physical report, a component and original
`reference_seed: {strand,index}`. It freshly replays the complete transcript,
proves a paired contour where supported and bridges only faithful representable
material. Other components, released activation, nicks and unpaired calls remain
in its full evidence artifact. There are no separate phase executors.

`topo_graph` with its explicit capture task retains the exact source-backed transcript in a schema2 native
graph record; graph continuation reconstructs original sources and retained steps
before applying new steps. The consumed activation and effective cut-end chemistry
are not reset at a save/reopen boundary. Graph imports and downloads are documentary
until fresh replay. Unsupported internal versions are refused.

`propose_plan`/`verify_plan` use the sole `TopoActions` physical bridge. Verification
returns `scientific_step_count` plus bounded full typed evidence, including
withheld products. `record_workspace_run` retains current physical/reference
sidecars with exact run/event ownership. `inspect_workspace_run_evidence` downloads
a bounded current scientific-evidence bundle without restoring material authority.
Fresh replay remains a separate action. [Graph history](TOPO_GRAPH_HISTORY.md) ·
[Chemical steps](TOPO_PROCESSING_ACTIONS.md).


## Native upstream PCR composition

Generic `propose_plan`/`verify_plan` support a `TopoActions` role bound to a prior
step only when that role is solely an ExactPcr template. Fresh native tracked
material is required; no serialized graph restores it. Historical run inspection
retains the exact upstream product identity. `inspect_saved_topo` resolves explicit immutable workspace/event-bound Run or Derived origins and freshly replays only the supported PCR-template role; prepared vectors remain explicit original inputs and unsupported staggered/nicked templates are refused without automatic repair. [Contract and limits](TOPO_UPSTREAM_PCR.md).

## Source-bound TA/GC assembly

`propose_tailing_assembly` accepts `{inputs, target, recipe}` with exact immutable
sources and an independent target. The recipe is `{name, profile, fragments,
circular, allow_nicked}`, where `profile` is `"Ta"` or `"Gc"`. Each fragment is
`{source: FragmentSpec, action, treatment}`. The action is `"KeepExisting"` or,
for example, `{"AddTail":{"base":"A"}}`. TA permits newly added A/T and GC
permits newly added G/C. KeepExisting emits no tail step and preserves the actual
prepared material, including selected saved origins. It does not infer end
chemistry from a stock name or imported annotation.

Preparations use the shared original-source PCR, digest, end repair and native
saved-material contracts. Optional `treatment` is `null`, `"Phosphorylate"` or
`"Dephosphorylate"`, applied after the selected keep/add action. Adding a tail
requires an established blunt duplex with accessible actual 3′ ends and adds one
declared nucleotide to both strands. Primer extensions remain separate.

Fresh compilation and verification produce an immutable candidate and bounded
complete `tailing-assembly-recipe.json`, including unsuccessful designs, physical
junctions, alternative products and retained nicks. The response exposes
`profile_status` and `profile_findings` separately from `verification.status`.
Successful generic blunt ligation can have an incompatible TA/GC profile;
matching the profile cannot supply a missing phosphate or establish isolation.

The response also includes plan/context revisions, preparation completeness and
an artifact ID. Candidates remain `verified: false`; `verify_plan`,
`export_protocol` and `export_protocol_pdf` provide fresh checks and protocol
handoff. Identical requests reuse immutable outputs. Source documents remain
unchanged. Polymerase activity, factory-vector chemistry, host repair, kit
performance and probabilities are not inferred. See the
[per-row recipe contract](TAILING_RECIPE_ACTIONS.md).

The existing overlap, ligation and tailing recipe tools also accept `Preparation::AnnealedOligos {bottom_input_id, options, candidate_index}`. Supply both original ordered oligos in `inputs`; a fragment's `input_id` is the top original and `bottom_input_id` identifies its separately revision-bound bottom original. Consequently three original inputs can yield two prepared ligation fragments, such as one pinned backbone digest plus one explicitly selected oligo pair. `options` declares pairing/search bounds and each oligo's terminal chemistry; annealing supplies no phosphorylation or covalent bond. Candidate availability and completeness are reconstructed freshly, and the independently supplied target never chooses a placement. The paired preparation itself adds no MCP tool or plan operation.

## PCR-product digest fragment inspection

`inspect_pcr_digest` accepts `{document, request, fragment_start, max_fragments}`. `document` is an exact immutable original-source reference; `request` is `{source_role, forward, reverse, selection, options}`, using typed primers (annealing domains and 5′ extensions separately), a full exact catalog pin and restriction-inspection limits. Pagination defaults to the first100 fragments and permits1–100 per call.

The backend freshly reconstructs actual unique PCR geometry and then inspects the resulting catalog digest. The response includes geometry count, bounded geometry preview, amplicon length, scan/product status, actual fragment/terminus summaries, findings and an immutable complete `pcr-digest-inspection.json` artifact. The fragment page retains its source/request context, and repeated identical inspection reuses the artifact. Competing PCR geometries withhold digest products; outer PCR5′OH and cut-derived5′P termini remain distinct. This tool creates no plan, target or product document and selects no fragment. Existing `propose_ligation_assembly`/`propose_overlap_assembly` recipes can subsequently bind a deliberately selected PCR/catalog fragment to the same original source and independent target.

## Saved-run material origins

`list_run_material_origins {workspace, event}` discovers exact retained run/step choices in an explicitly supplied immutable `WorkspaceRef`. It returns documentary selections and an immutable `run-material-origins.json` inspection artifact, without automatically choosing even a sole match. `prepare_run_material {workspace, selection}` deliberately chooses one run/step/exact event and freshly builds its bounded `Capsule`; retrieve it from the complete `run-material-preparation.json` artifact's `selected` field. Both tools check the exact source workspace and event ownership. Equal sequence letters in another document or event do not substitute for the requested source.

The capsule is a reproducible documentary request, not serialized live material. Existing recipe tools accept `Preparation::FromRun {origin: Capsule, preparation: RunPreparation}` with one flat explicit WholeLinear, Digest, CatalogDigest, ExactPcr, PcrCatalogDigest or MismatchPcr downstream preparation. They freshly reconstruct the original run again before preparing its actual material and comparing the separately chosen target. Discovery/preparation tools create no documents, candidates or targets, and an imported artifact alone grants no scientific authority.

`inspect_pcr_digest` additionally accepts optional `origin: Capsule` for an explicitly chosen saved-run material source. Without it, the tool and bare `pcr-digest-inspection.json` artifact remain unchanged, including rejection of saved physical-marker documents. With it, the tool freshly reconstructs that exact saved source before PCR/digestion and returns an outer origin-bound `context_revision` plus `origin_revision`. Its `run-pcr-digest-inspection.json` artifact is `{context_revision, origin_revision, preview}`, where `preview` has the existing ordinary inspection layout. No marker is deleted, no imported summary becomes material, and no fragment is automatically chosen. Tool count remains173.

When `inspect_pcr_digest` supplies an `origin`, it must also supply `workspace: WorkspaceRef`. This exact immutable native workspace must retain the selected source event and the original run's complete candidate, target event/snapshot, input bindings and product bindings. Documentary ownership is checked before the single fresh reconstruction/inspection pass; a self-consistent capsule cannot invent an event owner or substitute an equal-DNA alias. Omitting the workspace or changing its revision rejects the saved-source request. Calls with no origin keep the existing original-source API and artifact behavior.

## Direct catalog digestion of saved material

Existing `inspect_restriction` accepts optional `origin: Capsule` and `workspace: WorkspaceRef` (the workspace is required when origin is present). It checks exact retained source/run/independent-target/input/product bindings and event lineage before one fresh saved-material inspection. These documentary checks share the same server helper as saved-run PCR/digest inspection and run no duplicate reconstruction. Tool count stays173.

Without an origin, existing flat fields and both ordinary inspection artifact layouts remain unchanged. With an origin, both `artifact_id` identify a complete `run-catalog-restriction-inspection.json` wrapper: `{engine_version, context_revision, origin_revision, origin_status, source_modifications, fragment_modifications, inspection, findings}`. `inspection` retains the existing catalog report and batch format. The response adds outer `context_revision`, `origin_revision`, `origin_status` and `origin_findings`; ordinary site/fragment pages and actual termini remain available. Outer PCR hydroxyl groups, cut-derived phosphate groups and tracked modification summaries are preserved rather than inferred from saved sequence metadata. The tool neither publishes a document nor selects/installs a digest product. Cleavage is the declared complete-cleavage model, not evidence of enzyme activity, methylation permission, isolation or experimental success.


Shared overlap, ordered-ligation and TA/GC recipe tools also accept
`Preparation::MismatchPcr {forward, reverse, options}` (and the corresponding
flat `RunPreparation` under a saved origin). `options` is the recorded native
mismatch-PCR binding/search policy. These tools still require the separate target
and freshly verify against it; they never replace the target with a simulated
mutant product. Complete-amplicon reaction substrates preserve the mismatch
policy. Existing protocol HTML/PDF export describes the original mismatch-PCR
operation and explicit oligo sequences. The same native mismatch policy is retained throughout preparation and verification.

## Independently declared internal substitutions

`build_substituted_assembly_intent {inputs, request}` accepts1–32 exact immutable `BoundInput {id, document: DocumentRef}` snapshots and the native request `{base, junction_edits, substitutions}`. Each substitution is `{role, offset, expected, replacement}` with a zero-based offset in the original contribution's construct orientation and single concrete DNA bases. Expected source bases, retained-core placement, role identities and bounds are checked before copying/building. Junction trims/inserts and substitutions define the intended reference independently; primers and simulated products never supply its DNA. Repeated contributions of the same original source retain distinct roles and may carry separate declared changes.

The tool atomically returns an unverified intended-reference session target and complete `substituted-assembly-intent.json` artifact, including exact source refs, the native report and applied source/output coordinates. `verified` and `preparation_simulated` are false. Source snapshots, including historical imported revisions, stay immutable; stale expected bases or source refs reject without publishing a partial target. As with existing standalone assembly-intent tools, bindings cover MCP document roles/revisions only. This tool creates no native workspace/event ancestry, physical material, primer proposal or scientific acceptance; explicit native adoption establishes actual event parents separately.

`inspect_pcr_digest` additionally accepts optional `pcr_options: mutagenesis::PcrOptions`. Omitting it preserves the exact-PCR tool behavior, response and artifact layouts. Supplying it selects the distinct native mismatch-PCR/catalog inspector and `mismatch-pcr-digest-inspection.json` artifact. The request's existing primer/catalog fields remain in `request`; explicit mismatch/3′-exact/search policy is separate in `pcr_options`. Primer tails, encoded substitutions, candidate binding evidence and the actual selected amplicon remain inspectable in the full typed artifact.

Mismatch mode returns `mode: "mismatch-pcr"`, actual policy, native search/enumeration completeness, candidate summaries, selected amplicon summary, paged real fragments and `incomplete`. A certain nonoverlapping conditional singleton may be inspectable under a capped search; the incomplete flag stays explicit and no uniqueness or successful reaction is inferred. A native Run/Derived `origin` requires its exact owning `workspace` and shares the existing full ownership attestation before fresh science. No independent target, plan or product document is created. Existing shared recipe tools also accept explicit `MismatchPcrCatalogDigest` preparations and their flat saved-run equivalent; an explicit fragment index is verified separately against the independently declared target.


`verify_plan` retains a complete verification artifact whenever TOPO
`step_reports` or `annotation_reports` contain evidence, including graph-only
outcomes with no DNA contour. These use the existing artifact field alongside
physical and reference evidence. A plain plan with no typed
step-evidence arrays keeps its existing response. Artifacts remain documentary;
no intermediate or graph-only result is published as a molecule by verification.

## Explicit physical overlap previews

`preview_homology_recipe` binds immutable original `inputs`, a separate `target`,
an explicit preparation `recipe`, and a `homology_assembly::Options` policy. It
reconstructs declared 5′ resection, template-supported fill and eligible sealing.
`preview_annealing_recipe` takes the same source/target/recipe bindings and a
`homology_annealing::Options` policy with explicit per-oriented-source left/right
3′ cuts, junction lengths, topology and pairing-search limits. Neither operation
infers its physical policy from a vendor profile or selects source sequence to
match the target. Every saved `FromRun` or `FromDerived` role additionally requires its exact
`origin_workspaces` binding and fresh closed-origin replay.

Both previews save no candidate or product; their bounded JSON artifacts include
complete typed plans and physical/junction evidence. A client can explicitly
propose the retained operations with `propose_plan` and verify them separately.
The 3′ annealing route retains actual separate strands, missing bases and unsealed
nicks. Its joined reference always remains incomplete and cannot be recovered as
reusable assembled physical material. Standalone `ResectThreePrime` is available
as a typed plan operation and preserves its actual shortened substrate for later
explicit processing. [Scientific contract](THREE_PRIME_HOMOLOGY.md).


## Current saved-material and review tools

`material_derivation {workspace, task}` uses an immutable native `WorkspaceRef` and exact source `EventBinding`. An original reference is explicit; a saved origin is a small `MaterialSelection`, resolved and replayed natively rather than a client-supplied physical capsule. PCR preview exposes bounded geometry indices, lengths and a preview revision. Installation names the exact index and expected preview revision; the server reconstructs the full geometry itself. Exact/mismatch PCR, pinned catalog digestion and explicit end repair use the single derivation path. New standalone products retain exact native source events and derivation evidence; inherited incompleteness and search limits remain visible. Optional primer collections publish atomically with the new immutable workspace/archive. No tool writes a native desktop catalog or invents an independent assembly target. Later assembly verification still requires its separately declared target.

`prepare_saved_run_review {workspace, run_id}` prepares the exact retained native run with its original source and independent-target events. It returns the documentary review JSON and a portable unverified candidate artifact; it does not adopt the recorded verdict, material status or workspace heads. Fresh verification is required.

`export_map_svg {id, revision}` exports canonical features of that exact immutable native document to a bounded UTF-8 SVG artifact. It returns the exact source revision, warnings and artifact identity. The MCP route requests no restriction overlay; it makes no cloning/physical-material claim and changes no source or native catalog.


## Documentary plate layouts

`protocol_layout` uses an exact immutable native workspace and a typed task.
`Preview` resolves explicit96/384-well assignments to planned batch declarations
or released run/product events. `Capture` requires the reviewed context and
returns a new workspace; existing heads and historical assessments are unchanged.
`List` and `Open` inspect saved layouts. `Export` produces CSV, HTML or embedded
PDF with exact identities and explicit missing quantities. Assignment validity
does not freshly verify a construct or establish physical execution.

Ordinary selected-plasmid publication omits plate-layout sidecars.
`ReviewPublication` computes the complete layout dependency scope;
`ProjectPublication` requires its exact scope revision before producing a native
archive. Indivisible batch/release records may retain references outside assigned
wells, disclosed by the review counts and notices. This creates an export
projection, without uploading it to a database.
[Native model](PROTOCOL_LAYOUT.md) · [Transport and sharing](PROTOCOL_LAYOUT_TRANSPORT.md).


## Saved-template inverse PCR

The existing `design_inverse_pcr` tool accepts required nullable `origin` and
`origin_workspace` fields. Null/null selects the original-reference hypothesis.
A saved Run/Derived selection is resolved in its exact immutable workspace and
compared with the queried source before design. The target remains separately
supplied. A saved template adds native reconstruction before PCR, phosphorylation
and intended circular ligation. Unsupported nicks/strand geometry do not acquire
repair or source authority from metadata.

The returned context binds source, full origin, target and raw options; actual
ordered oligos use normalized DNA. Candidate and bounded evidence are published
atomically in the session store. Historical source findings remain visible, and
acceptance/export still require fresh ordinary verification.
[Native scope](INVERSE_PCR_MATERIAL.md) · [Wire and SDK contract](INVERSE_PCR_ORIGIN_TRANSPORT.md).


## Saved-source domestication and tailing libraries

`prepare_domestication` accepts explicit `request.origin` and `origin_workspace`
(required nullable values). Null/null declares an original reference; a native
Run/Derived origin requires its exact retained workspace and event. The source
session document must equal that event snapshot. The checked genotype is fixed
before inverse PCR, sequential inverse-PCR rounds or OE planning. Returned plans
contain the native reconstruction prefix and only new primer orders; inherited
source limitations remain visible. Missing/foreign owners and marker-only source
claims are refused. Imported Prepared reports confer no material authority.

Combinatorial `Method::Tailing` carries `profile`, `allow_nicked`, and one explicit
`{bin_id, action}` per source bin. Actions use the ordinary KeepExisting/AddTail
contract; per-construct whole-method overrides preserve deliberate differences.
Every computed tailing row retains typed `tailing_profile` evidence separately
from physical verification. A matched end family cannot repair a wrong independent
target, and a valid generic blunt ligation does not establish TA/GC profile match.
Current batch records are schema2. See [tailing libraries](TAILING_SOURCE_LIBRARIES.md)
and [saved domestication](DOMESTICATION_SAVED_MATERIAL.md).
