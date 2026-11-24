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
limitations](../crates/variants/README.md). VCF accepts bounded standard gzip (including concatenated members) and strict BGZF; compressed BCF remains strict BGZF. CRAM, BCF1, indexes and large-file random access remain pending.

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

The browser now opens VCF/BCF sources through **Open files** (including `.vcf.gz`
and `.bcf.gz`) or pasted VCF. A retained worker holds the complete private
collection, serves bounded record pages and freshly checks the explicit DNA
mapping and ALT choice. Opening a preview keeps the exact original file in the
workspace and derives a new document from the exact reference history event;
source DNA is unchanged. Reopen retained files from the history inspector.
Archive metadata does not restore parsed records or passing checks. Saved design metadata includes a bounded unverified record/ALT recipe. `replay_preview` requires the original private file collection and exact actual reference, reconstructs again, and compares literal DNA/topology; renaming a design is allowed, while a changed result is detected. Metadata and a matching replay do not restore earlier verification or chemical authority.

MCP exposes `import_variants`, `inspect_variants`, `check_variant_reference`,
`preview_variant_allele`, `export_variants`, `dispose_variants` and
`replay_variant_preview` (112 tools total). Collections are immutable, scoped to the session and charged before
insertion. Reference previews install atomically with bounded evidence. The
browser source profile is16MiB compressed/16MiB decoded,50,000 records and a
96MiB conservative reservation; MCP uses4MiB source/8MiB decoded,25,000 records
and48MiB reservation. These are managed budgets, not a measurement of process RSS.
The worker protocol appends schema29 jobs88–93/results80–85; older layouts stay
stable. Try the [synthetic example](../fixtures/variants-demo.md).

Native fixture type checks passed at this checkpoint. The development WASM
build passed and the immutable port8082 preview was refreshed; new fixtures are authored but unexecuted locally. Browser
checks remain deferred until in-scope implementation is complete.
