# collagene-read-alignments

MIT pure Rust, native and `wasm32-unknown-unknown`, with no subprocesses, runtime network requests, filesystem paths or external scientific programs. This standalone exchange layer is not yet integrated with the Collagene workspace or viewer.

`SamReader<R: BufRead>` and `BamReader<R: Read>` return one immutable `Record` at a time. `Header` retains exact SAM header text, ordered reference names/lengths/topology, and the original complete decompressed BAM header when present. `Record` retains flags (including unknown bits), reported positions and mate fields, CIGAR, sequence or its absence, raw Phred scores or their absence, ordered auxiliary fields, original record bytes, and a BAM BGZF virtual offset. Header dictionaries are indexed once; records bind both the dictionary and exact full header context, including unverified M5/UR/build claims. A same-name/same-length foreign header cannot silently reinterpret the record. Unknown tag **names** remain opaque evidence; unsupported type encodings fail explicitly.

```rust
use collagene_read_alignments::{Limits, SamReader, SamWriter};
use std::io::Cursor;
let source = "@SQ\tSN:ref\tLN:100\nr\t0\tref\t11\t30\t2M\t*\t0\t0\tAC\t!I\n";
let mut reader = SamReader::new(Cursor::new(source), Limits::default())?;
let mut writer = SamWriter::new(Vec::new(), reader.header(), Limits::default())?;
while let Some(record) = reader.next_record()? {
    assert_eq!(record.fields().start, Some(10));
    let map = record.coordinate_map(reader.header())?;
    assert!(map.usable_reported_mapping()); // Reported geometry only.
    writer.write_record(&record)?;
}
assert_eq!(reader.completion().unwrap().records, 1);
assert_eq!(writer.finish()?, source.as_bytes());
# Ok::<(), collagene_read_alignments::Error>(())
```

Read-only accessors expose private snapshots; there is no mutable record/header accessor and no JSON deserializer that can restore a completion claim. Clone a `Fields` value for downstream edits, but it cannot replace a record's original evidence through this API.

## Scientific scope

All positions exposed in typed fields are zero-based; intervals are half-open. A `CoordinateMap` has one run per reported CIGAR operation. Insertions are anchored between reference bases, deletions and skips remain different reference-only spans, soft clipping remains present query sequence, hard clipping remains absent calls, padding consumes neither coordinate, and M/=/X preserve different imported assertions. These are not reference-checked matches or variants. The unmapped flag withholds projections even when RNAME/POS/CIGAR are populated. Secondary and supplementary records remain separate; there is no unique-read, coverage, consensus or structural-variant inference. Missing SEQ/QUAL stay missing. `reported_query_length` can come from CIGAR when `sequence_present` is false; it is not an observed call count.

Reverse FLAG16 maps the SAM/BAM stored sequence orientation to the **implied** original full read, including CIGAR hard-clip extents. For `2H3S2M1I2=1X1D4N1P2S5H`, stored query length is 11 and implied original length is 18. Reverse stored calls 0/10 imply original positions 15/5: five omitted positions precede the reversed read and two follow it. No omitted calls or bases are generated. If CIGAR is absent, original full-read coordinates are unavailable. Circular wrapping requires explicit `@SQ TP:circular`; linear overruns withhold projections. Dictionary identity is not a genome-sequence hash or an attestation of an assembly release. Header sort order, M5/UR, MAPQ, NM/MD, SA, MM/ML and other evidence remain unverified source assertions; no UR is fetched.

## Fidelity and conversion

- Same-format SAM writes preserve exact header/record bytes, including CRLF, lexical values, tag order and an unterminated final line. Appending after an unterminated line fails rather than silently adding a delimiter.
- BAM-to-BAM writes preserve exact decompressed header and record payloads, including integer wire widths, float bits, unknown tag names, odd-base padding and original CG transport. BGZF is recompressed: original compressed framing is **not** preserved. Keep the original external stream by the completion source hash if exact compressed evidence is required. The library never buffers an entire source file to attach it to a record.
- Cross-format writes explicitly report formatting, integer representation, case and CG transport changes. BAM uppercase-normalizes representable SAM bases; nonrepresentable bases (including U or `.`) fail instead of silently becoming N. Non-finite BAM float payloads cannot be represented by supported SAM syntax and fail. A single BAM Q9 cannot become SAM `*` without losing its meaning, so conversion fails instead of adjusting the score. Empty arrays are supported. `writer.notices()` accumulates notices; inspect them before publishing output.
- Standard auxiliary A/i/f/Z/H/B types and BAM c/C/s/S/i/I widths are supported. Arrays retain their element type and raw integer/float bits. BAM CG:B:I restores >65,535-operation effective CIGAR only with its canonical placeholder. CG remains in raw evidence; SAM conversion removes only that transport field and writes the effective CIGAR. Zero-length operations, deprecated B/back CIGAR and nonstandard auxiliary type extensions are explicitly unsupported. SAM auxiliary integers use the BAM-representable range.
- BAM requires a matching explicit dictionary; conversion never invents one. BAM with binary-only dictionary can be read and rewritten as BAM, but SAM export requires the original matching SQ text dictionary; adding missing SQ lines or substituting header context needs a separate explicit evidence-bearing transformation, not this writer. Optional binary header NUL padding remains in `original_bam_header()`.

CRAM, VCF/BCF and BAI/CSI are explicitly rejected as alignment records. Ordinary gzip is not BGZF. `detect` returns only a `SamCandidate`/`BgzfCandidate`/`Unknown` probe; the reader establishes the actual transport/format. BAI/CSI querying and reference resolution are future work. Generic I/O and preserved virtual offsets leave a seam for bounded indexed adapters; no serialized index or imported claim grants integrity or biological authority.

## Resource and completion contract

Default limits are 8 MiB aggregate header, 4 MiB record/line, 100,000 references, 1,000,000 reported query bases, 100,000 CIGAR operations, 1,024 auxiliary tags, 1,000,000 array elements, 10,000,000 records, 64 GiB source bytes, 256 GiB decompressed bytes, and 10,000,000 BGZF blocks. Memory-policy maxima are 16 MiB header, 32 MiB record, 1,000,000 references, 10,000,000 bases, 1,000,000 operations, 4,096 tags and 4,000,000 array elements. Caller-selected stream/record counts remain finite u64 budgets with checked arithmetic.

BGZF reads compressed frames and inflates a **single** block at a time; advertised compressed size and decoded ISIZE cannot exceed 65,536. CRC, exact DEFLATE termination/consumption, final empty EOF block, partial headers, declared counts and cumulative budgets are checked. Embedded empty blocks do not end a stream. Allocation checks precede untrusted declared lengths. Memory is proportional to bounded header, one record and its typed fields (not only its raw wire size), one operation-map vector when requested, and one BGZF frame/output block. No whole-file decompression or per-reference-base mapping allocation occurs. Caller-owned sinks can still choose an in-memory whole output; streaming sinks are recommended for large data.

A reader emits `Completion` only after a caller asks for the next record and verified EOF is reached. Failures poison the reader; no completion hash is returned after truncation, corruption, limits or parse errors. Stopping early is incomplete. The receipt binds exact original source bytes, not validated biological content. Writers enforce record and cumulative budgets, poison on errors, and require explicit `finish`; there is no successful finalization after a failed record/write. A caller-owned sink can contain partial bytes after any streaming error and must be discarded rather than published.

## Primary references and dependency decision

Reviewed 2026-11-24: [GA4GH/SAMtools SAM/BAM 1.6, printing b5341fb dated 2025-08-12](https://samtools.github.io/hts-specs/SAMv1.pdf), especially alignment fields, optional fields, BGZF and BAM encodings; [optional-tag specification a71164e dated 2024-09-09](https://samtools.github.io/hts-specs/SAMtags.pdf); [canonical official specifications](https://github.com/samtools/hts-specs). All fixtures here were independently assembled, including literal little-endian BAM and malformed framing; no cloning application code or fixture corpus was borrowed.

[noodles](https://docs.rs/noodles/latest/noodles/) is MIT-compatible and offers broader format/index APIs. The cached 0.46 BGZF framing source and 0.88 BAM source were reviewed: application-level limits and strict partial-header/ISIZE checks are still needed before those APIs can allocate from wire declarations. This crate uses narrow explicit bounds and pure-Rust flate2 `rust_backend` (MIT/Apache-2.0), crc32fast (MIT/Apache-2.0), and blake3 (CC0/Apache-2.0/Apache-2.0 WITH LLVM-exception). Optional `interop-tests` enables noodles-sam0.84/noodles-bam0.88 for an independent reader check. It is optional because enabling noodles also enables flate2's alternative pure-Rust zlib-rs backend; default tests exercise the production miniz backend independently.

```sh
cargo test -p collagene-read-alignments --locked
cargo test -p collagene-read-alignments --locked --features interop-tests
cargo build -p collagene-read-alignments --locked --release --target wasm32-unknown-unknown --lib
```

## Shared sequence boundary

`reported_symbol(index)` distinguishes IUPAC, BAM/SAM `ReferenceEqual` (`=`), and unresolved symbols (`.` or non-DNA letters). `dna_calls()` returns a borrowed `DnaCallView` only for present all-IUPAC DNA; it rejects reference equality and missing/unresolved calls without replacing them with N or a supplied reference base. It reports whether an uppercase-normalizing host parser would change lexical case. A host can then explicitly pass `view.as_str()` to `collagene_model::Sequence::dna`, retaining this record's original bytes/hash, coordinates, missing quality and evidence warnings separately. The standalone crate intentionally has no dependency on an absolute root-model path. A constructed Sequence still does not establish matching reference DNA, full-read coverage, unique alignment, physical modification state or biological validation. Resolving `=` requires a separately bound actual reference and a new evidence-bearing transformation, outside this exchange API.

`progress()` reports emitted-record count and source bytes accounted so far with explicit Reading/Complete/Failed state. This is diagnostic progress, not a resumable checkpoint or a completion hash; BGZF input can read ahead to a block boundary. A partial failing I/O operation need not have all consumed bytes reflected in that diagnostic count. Only `completion()` after verified EOF binds the exact full source.

The immutable source/evidence binding cannot be changed through the public API:

```compile_fail
use collagene_read_alignments::Record;
fn forge(record: &mut Record) { record.fields.flags = 0; }
```

```compile_fail
use collagene_read_alignments::Header;
fn replace_dictionary(header: &mut Header) { header.references = Vec::new().into(); }
```

`CoordinateMap` is also a sealed derived snapshot. Its getters return values or immutable slices; untrusted mutable report fields cannot enter projection arithmetic. No API accepts a caller-constructed map for projection.

```compile_fail
use collagene_read_alignments::CoordinateMap;
fn forge(map: &mut CoordinateMap) { map.circular = true; map.reference_length = None; }
```

Header lines use the same `max_record_bytes` line cap in constructors, readers and writers, in addition to aggregate header limits. Export under a tighter policy rechecks header lines, reference count and complete encoded BAM header overhead **before** allocation or writing. SAM cross-format rendering uses one bounded formatting/counting pass before an exactly sized String allocation. BAM encoding pre-counts fixed fields, packed calls, qualities, all auxiliary payloads and long-CIGAR transport before allocation. Both counts use the smaller of the record cap and remaining cumulative budget. A compact BAM auxiliary array cannot bypass a SAM output cap through decimal expansion. Writers retain only sealed header identity and needed shared lookup data after emitting the header, rather than cloning unused original BAM evidence.
