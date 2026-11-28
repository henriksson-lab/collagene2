# Exact FASTA record retention: current source contract

28 November 2026. Native and WASM compilation, 700 focused cases and the default memory-only real MCP stdio check (229 tools) pass. [Exact sources, cases and execution limits](fasta-retention-current-integration.json). No provider/live-service or GUI acceptance is claimed.

## One parser and one existing metadata representation

[formats::sequence](../../crates/core/src/formats/sequence/mod.rs) owns local and downloaded FASTA parsing and exact record retention. Existing `import_bytes`, `import_text`, `import_fasta` and `import_fasta_with_alphabet` entry points retain each record through the same private `FastaInput` owner. Format detection uses a trimmed view, while FASTA parsing/retention receives the original input slice. Explicit alphabet/header-marker conflict rules, normalized calls, literal header names, linear topology and absent feature annotations remain unchanged. A standard byte-zero UTF-8 BOM is accepted and retained.

The existing keys are reused directly:

| Key | Meaning |
| --- | --- |
| `source:original_fasta_base64` | Base64 of the exact corresponding original record bytes. |
| `source:original_fasta_blake3` | BLAKE3 of those decoded bytes. |

There is no new raw-record encoding, internal-format compatibility path, molecule/CG2 version bump or provider claim. Normal FASTA export remains a normalized sequence presentation, not a reconstruction or replacement of the retained original.

## Exact record partition

The raw file is partitioned into disjoint borrowed record slices before encoding:

- The first record starts at byte zero, retaining any standard initial BOM, leading blank lines and whitespace.
- A later record starts at the beginning of its physical header line, including indentation before `>`.
- Blank lines and line endings before the next header belong to the preceding record.
- The last record includes all trailing bytes, including whitespace and an absent final newline.

CRLF, header/sequence casing, indentation and internal sequence whitespace remain byte-exact in the source record while the ordinary sequence parser normalizes calls. Decoded records concatenate to the original input. A multi-record file is never copied into every molecule; each record receives only its own source slice. Exact bytes and hashes are evidence of retained data integrity, not sequencing completeness, read quality, sample identity, function or physical material.

## Admission before amplification

The existing input cap is 32 MiB and the record cap is 1,000. Native molecule names are admitted at 1,024 bytes before copying. Existing per-sequence length/alphabet validation still applies.

`FastaInput::new` scans borrowed lines and bounds the complete normalized-model plus encoded-original footprint before allocating headers, sequence buffers or base64 strings. Its conservative aggregate calculation includes per-record base64 padding, header JSON escaping and structural/warning overhead, and refuses a result exceeding 32 MiB. A raw file fitting the input cap therefore does not automatically fit the retained import report. No truncated source record is returned to force an import through a bound.

Each attached molecule passes the canonical native annotation/import-scale guard. Public `ImportReport::validate_bounds()` then checks the borrowed molecule count, canonical aggregate native-frame admission and streamed JSON report size against 32 MiB. The same report method validates final URL provider decoration. It does not clone or serialize an unbounded temporary output buffer and does not establish provider/scientific authority.

[ModelBudget](../../crates/model/src/bounds.rs) now classifies the **existing** FASTA base64 key through the same known-original owner as trace, SnapGene, original XML/XDNA, FASTQ and supported annotation-original keys. It no longer treats retained FASTA bytes as ordinary metadata with a 1 MiB value cap. Classification is consistent in borrowed guards and native decoding; ordinary metadata retains its existing limit.

The distinct native resource ledgers remain:

| Admission profile | Actual original-byte accounting |
| --- | --- |
| Ordinary scientific `model_preflight` | Each known-original encoded value is bounded at 4 MiB and shares the existing aggregate 8 MiB text budget with other admitted text. These are encoded byte counts, not a promise that 4 MiB raw FASTA fits. |
| Native annotation/edit/import-scale model | Original records share a separate 32 MiB aggregate allowance; ordinary text retains its separate 8 MiB allowance. |
| Native frame | Each model first passes import-scale admission; all models then share the declared frame allowance. This FASTA report uses 32 MiB, including encoded originals and other retained text. |

No scientific or native ceiling is raised to admit arbitrary files. Multiple known originals share the appropriate ledger rather than obtaining independent unlimited allowances. Successful import/native save does not promise that every scientific operation accepts that size.

## URL imports retain the same source records atomically

[sources::import_download/annotate_import](../../crates/sources/src/lib.rs) reuse the core parser; the former duplicate FASTA splitter and encoder are removed. `import_download` parses once and decorates only that fresh result. Public `annotate_import` admits the supplied models, reparses through the canonical importer and requires exact equality of both molecules and parser warnings before adding any provider facts. Edited/unrelated calls, names, annotations, raw-record hashes or warnings cannot acquire fetched-source provenance.

Provider decoration reserves aggregate space before copying URLs/metadata across records, preserves provider-specific accession/content checks, validates native saves and runs the final report guard. The public caller's report is replaced only after the complete staged operation succeeds; a late failure leaves it unchanged. URL metadata is documentary fetched-source information, not independent biological authentication. There is no provider network or credential-based acceptance check in this audit.

## Supplied-consensus provenance and authored coverage

The existing [design-consensus provenance owner](../../crates/core/src/alignment/design_consensus/provenance.rs) can now inspect locally imported FASTA through the same keys previously used for downloaded FASTA. It checks bounded decoded bytes/digest and reparses the actual source record, distinguishing matching normalized calls, different current calls and invalid/unavailable originals. Native edits retain the original bytes; those bytes never overwrite edited current calls. This does not manufacture per-base qualities or measured evidence.

Passing fixtures cover literal BOM/CRLF/whitespace partitions, RNA/protein markers and conflicts, native CG2 reopen, 257 disjoint records, retained values larger than the former ordinary-metadata cap, aggregate amplification refusal, model-ledger classification, local/URL parity, atomic annotation refusal, MCP import retention and supplied-consensus original-versus-edited-call evidence. The model-ledger fixture was corrected before its first execution to reflect the existing 4 MiB per-value and shared 8 MiB text limits; production bounds did not change during compilation.
