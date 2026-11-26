# SnapGene DNA/protein import profiles

`collagene_core::snapgene::import_snapgene` is an independently implemented MIT importer, not a wrapper around another cloning application. `original_snapgene_bytes` recovers an original source attachment; it does **not** serialize the current edited document as SnapGene. No external executable is invoked.

## Mapped data

The importer accepts the known uncompressed DNA profile (cookie type1, one packet0) with circular/linear topology and explicitly double-stranded DNA, and the selected linear protein profile (cookie type2, one uncompressed packet21, property byte0). The cookie determines the alphabet: protein letters A/C/G/T are not guessed to be DNA. Alphabets are validated and case is normalized. Other cookie types, compressed packet1, RNA packet32, conflicting/mixed DNA/protein packets, duplicate sequence packets, unknown property bits and producer/import revisions above15 reject. Single-stranded DNA and circular/nonzero-property proteins refuse. The accepted revision range is a structural profile, not a claim that every producing SnapGene release has been tested.

Standard feature segments use one-based inclusive source ranges. Internal locations are zero-based half-open. Circular origin crossings, multiple segments and reverse complements are represented as location trees. Duplicate qualifier values and flags are retained. Supported uniform segment colors become ApE color qualifiers. Non-directional/bidirectional arrow styling, named subsegments, differing segment colors and unmapped attributes produce explicit warnings/review flags. Gaps or unknown segment kinds become unresolved locations; no joined biological sequence is asserted for them.

Protein annotation ranges use amino-acid coordinates, with the same one-based-inclusive to half-open mapping. Exact uniformly forward domains/multipart features are documentary annotations. Reverse protein features and origin-wrapping ranges refuse; amino acids are never reverse-complemented. Protein crosslink/bond feature types retain unresolved reviewed locations and original evidence because this profile does not establish segment-to-bond endpoint semantics. It never converts a crosslink span into a contiguous peptide or proven covalent bond.

Mapped notes remain source-prefixed metadata, including any claim of experimental confirmation. Such claims never set verification status or authorize tools. Rich-text feature values remain literal data; this importer does not reproduce HTML formatting or execute it.

Primers/binding predictions, history, trace data, attachments, preferences and unmapped packets are not converted. Packet counts and explicit fidelity warnings report that limitation. DNA methylation/property bits and additional end/overhang packet8 set `source:snapgene:chemistry_unmapped=true`; the cloning core must not treat those source assertions as a fully modeled physical intermediate. Protein imports explicitly describe primary-sequence/annotation reference scope in `protein:authority`, not PTM, disulfide, terminal processing, folding, material state or experimental verification. Additional protein properties packet8 remain opaque and also set `protein:modifications`, withholding a quantitative unmodified/free-chain material interpretation. Annotated modifications remain source claims and cannot establish physical chemistry.

## Original evidence and bounds

Successful imports retain all original bytes as base64 with a BLAKE3 integrity hash in source-prefixed metadata. `.cg2` persistence preserves that evidence. Standard GenBank/FASTA exports are derived views and cannot preserve the complete SnapGene document. Source extraction still returns the original bytes after edits and is intentionally distinct from writing an edited `.dna` file. Checksums detect accidental corruption, not provenance authenticity.

Limits are 32 MiB input, 4,096 packets, 8 MiB per parsed XML packet, 100,000 XML nodes, 32 element levels, 10,000 features/segments per feature, and 100 deduplicated fidelity warnings with an explicit truncation notice. Parsed XML forbids DTD/entity declarations and external resolution. Duplicate singleton packets, truncated lengths, invalid coordinates and conflicting qualifier encodings fail explicitly. The resulting serialized molecule is also limited to 32 MiB so source retention cannot silently create an unsavable native document.

## Layout provenance and confidence

The cookie/packet framing and annotation structure were researched from [Damien Goutte-Gattat's original format investigation](https://incenp.org/dvlpt/docs/binary-sequence-formats/index.html), and compared narrowly with the [sgffp author's newer reverse-engineered layout notes](https://merv1n34k.github.io/sgffp/format-spec). These describe observations, not an official complete SnapGene specification. SnapGene's [developer guidance](https://www.snapgene.com/features/convert-file-formats) directs developers to contact the vendor for format information.

One discrepancy required a careful format-only check: the older prose assigns topology to a different flag bit. The newer layout notes and the [DNA packet interpretation in Biopython](https://github.com/biopython/biopython/blob/master/Bio/SeqIO/SnapGeneIO.py) agree that bit 0 denotes circular topology. The implementation uses that interpretation, with bit 1 declaring double-stranded DNA. Only layout/coordinate details were inspected; no application or scientific code was copied, translated or used as an algorithm correctness oracle.

Tests construct packets and XML independently, with hand-specified sequences and coordinates. They cover topology, origin-spanning reverse compound locations, duplicate qualifiers, integrity-preserved source recovery, unmapped history/primers/chemistry, malformed lengths/duplicates, invalid locations, DTDs and nesting. A representative corpus of real vendor-produced files, compressed sequence profiles, exact unknown feature semantics, primers and experimental history remain required before claiming robust migration support. Current support is the bounded DNA/linear-protein import profile and its separately documented selected writer, not complete vendor-workspace migration.


## Linear protein read/write increment

The [sgffp author's published packet21 description](https://merv1n34k.github.io/sgffp/format-spec#block-21--protein-sequence)
identifies the same uncompressed payload structure as packet0, with amino-acid
single-letter data; its cookie table identifies protein as type2. This increment
uses only those format-layout observations, not third-party parser/application
code. No extra compressed packet profile, RNA, trace/index or physical-protein
chemistry feature was added.

Five independent literal protein groups establish cookie/packet/AA-feature behavior
before a secondary round trip with the [selected writer](SNAPGENE_EXPORT.md).
They cover case/original-byte/native preservation, cookie-versus-letter alphabet
selection, duplicate/mixed/compressed/RNA and nonzero/circular flags, invalid/truncated
sequences, reverse/wrapped/out-of-bounds residue features, unresolved crosslinks,
additional-properties chemistry guards and opaque-history retention. Existing DNA
fixtures and guards remain in place, including type2-cookie plus DNA-packet refusal.
No vendor GUI, browser, compiler or scientific executable is used by this source
worker; the parent executes the coordinated fixture/compiler checks.
