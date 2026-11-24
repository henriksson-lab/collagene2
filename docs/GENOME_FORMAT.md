# Immutable large-reference storage

`core::genome` stores genomic reference data independently of `Molecule` and
its ten-million-base limit. It implements streaming FASTA ingestion, packed
IUPAC blocks, immutable metadata/layout identities, portable `.cgref` files,
local block verification and full canonical sequence verification. See
[GENOME_CRISPR.md](GENOME_CRISPR.md) for interval/global-coordinate operations
and private cooperative off-target search.

## Data and identities

Canonical sequence calls are uppercase IUPAC DNA. All fifteen nonempty base
masks are preserved. Lowercase calls normalize to uppercase; whitespace is not
part of the sequence. No gap, digit, invalid byte or uncertain residue is
silently converted to a concrete base.

Contigs have distinct stable role IDs, source header IDs, optional accession
versions, descriptions, topology, `u64` lengths and canonical DNA hashes. Equal
DNA in two roles remains two loci. Expected component records can map an input
header ID to a declared stable role and accession. Duplicate source IDs or role
collisions are errors.

Three identities serve different purposes:

* `sequence_hash` is BLAKE3 over the canonical contig DNA calls in order.
* `content_revision` hashes a domain-tagged, role-sorted collection of role ID,
  accession, length, topology and sequence hash. It is independent of FASTA
  record order, descriptions, line wrapping, case, build labels and file layout.
* `manifest_hash` hashes the exact canonical serialized manifest, including
  block offsets/hashes, build metadata, expected components and source evidence.

Jobs and portable links must pin **both** manifest and content identities.
Changing a provider label or expected-component declaration can change the
manifest without changing canonical reference content.

Build labels, assembly/version strings and provenance are recorded declarations,
not independent authentication of a provider's real genome assembly. Missing
expected components are listed and prevent `expected_complete`. A present
component contradicting a supplied expected identity, length, topology or DNA
hash is an error. Without an expected component list, the result is a declared
reference collection; a name alone cannot establish a complete build.

## Streaming FASTA profile

`FastaBuilder::begin`, bounded `push` calls, and `finish` never construct an
entire contig. Headers begin with `>` at the beginning of a line, IDs are unique
bounded first tokens, and header metadata is UTF-8. LF and CRLF headers are
supported. Sequence data can be wrapped or one arbitrarily long line. ASCII
space/tab/CR/LF are ignored between calls. Empty contigs and sequence before a
header are errors. This initial profile does not accept gzip or alignment gaps.

The builder buffers at most one 65,536-call DNA block, emits packed block events
for host persistence, and retains bounded metadata. Each event supplies its
exact append-only file offset, contig role, block index, original contig start,
count and encoded-byte hash. The host checks/writes the offset and can release
the event immediately. A parser, byte, base, record or metadata-bound failure
poisons the builder; it cannot later produce a positive finished reference.
Hosts discard its staging data.

The manifest retains the raw input byte hash/count and parsed header IDs and
descriptions. **It does not retain verbatim original FASTA evidence**: original
whitespace, case and header formatting cannot be reconstructed from `.cgref`.
A reconstructed FASTA would be a derived normalized export. Preserve the
original input separately when verbatim evidence is needed.

## Portable `.cgref` version 1

The layout is append-only and independently versioned from ordinary document
and worker formats:

| Section | Encoding |
|---|---|
| Header | Eight bytes: `43 47 52 45 46 00 01 00` |
| Sequence blocks | Packed bytes in manifest contig/block order |
| Manifest | Version-1 typed Postcard payload, at most 16 MiB |
| Tail | 48 bytes: manifest byte length as little-endian `u64`, 32 raw BLAKE3 hash bytes, eight-byte `CGREND1\0` magic |

Each block holds at most 65,536 calls. The first call occupies the high nibble;
the second occupies the low nibble. The nibble is the conventional nonempty
A/C/G/T bitmask (1/2/4/8). Zero is invalid for a call; the unused low nibble of an
odd-length final block must be zero. Blocks before a contig's final block must
have the full canonical block size. Descriptors cover the contig without gaps
or overlaps and cover the packed file data contiguously after the header.

`parse_tail` validates framing and calculates bounded random-read ranges.
`decode_manifest` checks the exact manifest length/hash, runs a no-allocation
typed shape/count preflight before deserializing, rejects noncanonical/trailing
encoding, then validates layout and metadata consistency. File lengths and
offsets use `u64`; only bounded loaded blocks become local `usize` indices.

## Integrity and scientific authority

`GenomeManifest::validate` and `FinishedReference::validate` check bounded
metadata, identities and layout. They do not read all DNA or attest that its
canonical hashes are correct.

`BlockVerifier::begin` validates/caches the layout and manifest hash once.
`verify` checks a named role/block, exact encoded length/hash and nibble/padding
semantics. It returns a `VerifiedBlock` with private construction and read-only
DNA/identity accessors. Local interval operations can use these checks while
explicitly disclosing that whole-contig canonical hashes remain unverified.

`SemanticValidator` accepts every canonical block exactly once, in manifest
contig/block order. It rejects changed identities, omissions, duplicates and
reordering, and recomputes each contig's canonical DNA hash at EOF. An error
poisons the validator. Only complete successful EOF yields a receipt. Hosts
must own fresh validator state; an imported/serialized receipt or status flag
cannot restore current authority. Halo blocks used in analysis are not pushed
again as canonical coverage.

Even successful canonical hash verification authenticates data against the
pinned manifest, not the biological truth of its provider labels. Whole-build
component coverage and exhaustive scientific search are separate checks.

## Bounds, adapters and portable links

Independent limits include 10,000 contigs, 250,000 total blocks, a 16 MiB
manifest, 4 KiB headers and bounded provenance. Defaults allow up to 1 MiB per
push, 32 GiB raw source bytes, 16 billion total bases and 10 billion bases per
contig. Other count/metadata limits can reject a source before its base ceiling;
these are ceilings, not a promise that every maximal combination fits. All
source hashing is incremental and block events are bounded by input size.

The browser adapter persists staging blocks separately from ordinary documents
and publishes only after readback/semantic verification. File reads use bounded
slices, rather than the ordinary document import's whole-file buffer.

The native MCP adapter owns private temporary files outside the 64 MiB document
store. It accepts opaque epoch-bound handles and exact upload offsets, never
caller-selected paths. It limits active imports/datasets, reserves disk quota,
reads back at most 64 blocks per verification step, and publishes no dataset
before semantic EOF. Parser/hash/I/O errors remove staging files; offset errors
leave valid uploads unchanged. Unix parent directories use mode 0700; Windows
inherits the temporary-directory ACL. Handles close before owned files are
deleted, including cancellation and server shutdown.

`ReferenceLink` is a small schema-1 `genome-reference-link-v1` workspace
attachment. It preserves a logical role, manifest/content identities and build
summary. Its optional local cache hint is not a path, permission or scientific
authority. Missing stores require reattachment of the same immutable identity;
loading a link does not restore source integrity or completed search status.
The packed reference remains a separate portable file, not a multi-gigabyte
ordinary molecule or workspace payload.
