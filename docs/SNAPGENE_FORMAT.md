# SnapGene DNA import profile

`collagene_core::snapgene::import_snapgene` is an independently implemented MIT importer, not a wrapper around another cloning application. `original_snapgene_bytes` recovers an original source attachment; it does **not** serialize the current edited document as SnapGene. No external executable is invoked.

## Mapped data

The importer accepts the known uncompressed DNA packet profile with a SnapGene cookie, circular/linear topology, and explicitly double-stranded DNA. DNA alphabets are validated and case is normalized. Cookie sequence types other than DNA, compressed sequence packets, unknown property bits, single-stranded DNA, and producer/import revisions above 15 are rejected. The accepted revision range is a structural profile, not a claim that every producing SnapGene release has been tested.

Standard feature segments use one-based inclusive source ranges. Internal locations are zero-based half-open. Circular origin crossings, multiple segments and reverse complements are represented as location trees. Duplicate qualifier values and flags are retained. Supported uniform segment colors become ApE color qualifiers. Non-directional/bidirectional arrow styling, named subsegments, differing segment colors and unmapped attributes produce explicit warnings/review flags. Gaps or unknown segment kinds become unresolved locations; no joined biological sequence is asserted for them.

Mapped notes remain source-prefixed metadata, including any claim of experimental confirmation. Such claims never set verification status or authorize tools. Rich-text feature values remain literal data; this importer does not reproduce HTML formatting or execute it.

Primers/binding predictions, history, trace data, attachments, preferences and unmapped packets are not converted. Packet counts and explicit fidelity warnings report that limitation. Methylation/property bits and additional end/overhang packets set `source:snapgene:chemistry_unmapped=true`; the cloning core must not treat those source assertions as a fully modeled physical intermediate.

## Original evidence and bounds

Successful imports retain all original bytes as base64 with a BLAKE3 integrity hash in source-prefixed metadata. `.cg2` persistence preserves that evidence. Standard GenBank/FASTA exports are derived views and cannot preserve the complete SnapGene document. Source extraction still returns the original bytes after edits and is intentionally distinct from writing an edited `.dna` file. Checksums detect accidental corruption, not provenance authenticity.

Limits are 32 MiB input, 4,096 packets, 8 MiB per parsed XML packet, 100,000 XML nodes, 32 element levels, 10,000 features/segments per feature, and 100 deduplicated fidelity warnings with an explicit truncation notice. Parsed XML forbids DTD/entity declarations and external resolution. Duplicate singleton packets, truncated lengths, invalid coordinates and conflicting qualifier encodings fail explicitly. The resulting serialized molecule is also limited to 32 MiB so source retention cannot silently create an unsavable native document.

## Layout provenance and confidence

The cookie/packet framing and annotation structure were researched from [Damien Goutte-Gattat's original format investigation](https://incenp.org/dvlpt/docs/binary-sequence-formats/index.html), and compared narrowly with the [sgffp author's newer reverse-engineered layout notes](https://merv1n34k.github.io/sgffp/format-spec). These describe observations, not an official complete SnapGene specification. SnapGene's [developer guidance](https://www.snapgene.com/features/convert-file-formats) directs developers to contact the vendor for format information.

One discrepancy required a careful format-only check: the older prose assigns topology to a different flag bit. The newer layout notes and the [DNA packet interpretation in Biopython](https://github.com/biopython/biopython/blob/master/Bio/SeqIO/SnapGeneIO.py) agree that bit 0 denotes circular topology. The implementation uses that interpretation, with bit 1 declaring double-stranded DNA. Only layout/coordinate details were inspected; no application or scientific code was copied, translated or used as an algorithm correctness oracle.

Tests construct packets and XML independently, with hand-specified sequences and coordinates. They cover topology, origin-spanning reverse compound locations, duplicate qualifiers, integrity-preserved source recovery, unmapped history/primers/chemistry, malformed lengths/duplicates, invalid locations, DTDs and nesting. A representative corpus of real vendor-produced files, compressed sequence profiles, exact unknown feature semantics, primers and experimental history remain required before claiming robust migration support. Current support is the planned import spike.
