# BAI and CSI index groundwork

`collagene-genomic-index` is a pure bounded format library. It supplies candidate BGZF chunks, not a claim that an index matches any BAM/BCF/VCF source or proves absence. No GUI indexed browsing, variant calling or biological validation is introduced by this crate.

The implementation was independently authored from the official [SAM/BAM specification, BAI section 5.2](https://samtools.github.io/hts-specs/SAMv1.pdf) and [CSI specification](https://samtools.github.io/hts-specs/CSIv1.pdf), consulted 2026-11-24. The SAM PDF identifies version 1.6, revision b5341fb dated 2025-08-12. These primary specifications define packed virtual offsets and generalized hierarchical bins. No external scientific program or cloning-app code is used.

The API accepts explicit `Format::Bai` or `Format::Csi` and uncompressed index bytes. The core BAM adapter separately decodes CSI BGZF transport under compressed/decoded/block caps, CRC and exact terminal-EOF checks; arbitrary gzip is not accepted by that adapter. BAI metadata 37450 and generalized CSI metadata(bin_limit+1) are separated from actual chunks; mapped/placed-unmapped u64 counts are never interpreted as offsets. Opaque CSI auxiliary bytes remain exact.

Index files do not carry trustworthy source identity or a reference-name dictionary. Hosts must retain actual index bytes/hash and actual source bytes/hash, compare dictionary ordinal count/order, check BGZF address/block boundaries and filter decoded records against the requested interval. `Association` is an explicitly caller-declared hash/size/ordered-name tuple; optional binding checks numeric source bounds and subsequent expected hash but does not authenticate the declaration. Querying without it is explicitly ordinal candidate lookup. A complete file parse is not a complete genome or alignment search. An empty/corrupt/incomplete index must not become an absence claim; the host needs an independent retained-record audit or another verified indexing mechanism.

CSI query calculations follow the generalized hierarchical half-open bin ranges, with depth<=10 and total coordinate bits<=63 to fit checked u64/u32 representations. BAI uses its fixed 29-bit coordinate range. Per-level binary searches traverse only stored overlapping bins; candidate chunks are capped, sorted and merged without per-base arrays. The current conservative query intentionally does not prune by imported BAI linear offsets or CSI loffsets. This can return extra candidate records and cannot skip a candidate solely on those imported optimization assertions.

All retained allocation classes are counted in a first structural pass before cloning/allocating parsed data. A second pass constructs sealed state. Query chunk/work limits are checked before candidate vector allocation; duplicate bins, malformed metadata, reversed real/reference ranges, negative signed counts, truncation and unexpected tails fail atomically. The crate is registered in the workspace, shared core and native/WASM CI jobs. Focused fixtures cover literal index framing and boundary queries.


`core::indexed_bam::IndexedBam::bind` accepts only an actually parsed, completed
BAM collection. It retains the actual source hash and raw index hash, compares
reference counts, and resolves names through the actual BAM dictionary order.
The core adapter accepts raw CSI or strict BGZF-compressed CSI and retains the exact original transport hash. The initial `select` adapter checks candidate offsets against source byte bounds
and scans every retained record to independently check interval omissions. A
matching mapped record outside all returned chunks rejects the whole query.

Selection uses zero-based half-open *reported alignment spans*. CIGAR D/N consume
reference; clips and insertions do not. A zero-reference-consuming mapped record
occupies one reference position. Placed unmapped records are excluded, while
secondary/supplementary mapped records remain. All matching indices are returned
atomically under explicit audit/match limits. A deletion/skip intersection is not
an observed aligned base or coverage. No actual-reference check is inferred.

The full scan intentionally costs O(retained records) and supports the current
bounded, fully imported collection. Streamed indexed decoding, HTTP ranges,
TBI/CRAI and variant-index hosts remain future work.
This adapter does not establish full real-world index interoperability or search
performance for large sources.


`core::indexed_bam::build_bai` can also rebuild a BAI from the actual completed
BAM. It computes bins from reported spans and chunks from actual record offsets,
checks real coordinate/reference ordering (including the unplaced tail), and
includes placed unmapped records. It emits no imported metadata counts or linear
pruning claims. Construction is atomic and bounded by output and managed-memory
limits; coordinates beyond the BAI 29-bit scheme require CSI. CSI writing remains
pending. Rebuilding an index does not validate alignment biology.
