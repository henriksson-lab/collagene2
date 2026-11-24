# Variant exchange and explicit DNA checks

`collagene-variants` is MIT-licensed exchange infrastructure for declared VCF4
and BCF2.2 evidence. It uses pinned, inspected Noodles parsers with additional
allocation, framing and codec guards. It is separate from sequence documents,
variant calling, sample inventory and laboratory workflows.

An immutable `Collection` retains exact original bytes, a BLAKE3 identity and
complete parsed header/record state. Typed INFO, FORMAT, missing values, sample
column order and genotype phasing remain available through readonly record APIs.
Normalized pages have at most100 records and a separate aggregate text bound.
Original export is byte-exact; conversion has explicit normalization notices and
returns no partial artifact on errors. See the [crate's bounds and codec
limitations](../crates/variants/README.md). CRAM, BCF1, compressed VCF, indexes and
large-file random access remain pending.

Preflight bounds sparse structured-header IDX slots before the upstream parser
can expand indexed maps. BCF extended typed lengths, sample dimensions and raw
payload consumption are checked before lazy record indexing. Unsupported codec
combinations fail explicitly rather than calling unimplemented upstream paths.
BGZF input is framed, size-bounded and CRC-checked before decoded blocks are
retained. A completed parse establishes a file receipt, not allele truth.

`core::variant_reference::check_record` binds a collection hash, record index,
explicit file-contig mapping, actual molecule revision and one-based alternative
index. It compares concrete REF bases at the record's actual forward coordinates
with the selected DNA and checks a declared contig length when available.
Matching names and genotypes never select or authorize an edit. Missing positions,
symbolic/breakend/spanning alleles, mismatches, ambiguous DNA and unsupported
ranges remain explicit; circular wrapping is not inferred.

`preview` recomputes the check before applying a concrete selected allele through
the existing source-preserving sequence editor. Annotation/citation review and
source history remain subject to that editor's rules. A decoded report is never
accepted as an editing proof. The preview retains a chemistry-unmapped reference guard: a file allele cannot silently become prepared cloning material. It can still be an independently specified design target. This creates a selected-allele design preview;
independent target verification and actual genotype inference are separate tasks.

This checkpoint exposes native/WASM library APIs. Browser import and MCP session
integration are not yet wired. Native core/MCP/web fixture type checks and WASM
compiler checks passed; new fixtures are authored but unexecuted locally. Browser
checks remain deferred until in-scope implementation is complete.
