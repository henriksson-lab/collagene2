# Explicit single-nucleotide 3′ tailing

`tailing::tail_blunt_duplex(source, TailingBase::{A,C,G,T})` applies an ideal declared endpoint transform to an established physical duplex. It appends exactly one chosen untemplated nucleotide to the 3′ suffix of each actual strand, both stored 5′→3′. This module does not infer that a polymerase performed the transform, predict tail heterogeneity/yield, or implement a commercial kit's chemistry.

The operation requires concrete DNA, a linear fully paired blunt duplex, and known accessible 3′ hydroxyls at both affected ends. Circular, staggered, discordant, non-DNA, unknown/blocked 3′ and invalid-nick sources fail explicitly. A staggered source must undergo its own supported repair operation first. Retained internal paired-region nicks and both 5′ chemical states remain unchanged; no phosphorylation, nick sealing or deblocking is inferred.

## Coordinates and ancestry

For an input pair of length `n` with bottom offset zero, the product has two strands of length `n+1` and bottom offset `-1`. The new bottom 3′ residue extends left in top-reference coordinates, while the new top 3′ residue extends right. Each free end therefore exposes a one-base 3′ overhang of the declared base. Reversing the duplex and then tailing gives the same physical product as reversing the tailed duplex, including transformed nick boundaries and swapped end chemistry.

`TailingReport` binds the engine, canonical serialized-source BLAKE3 digest, selected base, exact source/product strands and chemistry. Its two `StrandMapping` entries cover only truly retained input residues. `SynthesizedTail` records each new residue's actual product-strand coordinates and sequence; it has no fabricated source/template coordinate. The new top reference has length `n+1`; the unpaired left bottom residue is retained in physical evidence rather than invented as an extra top-reference base.

`TailingReport::validate()` replays the complete transform from its retained source/base and compares all derived fields, mappings, additions, warnings and digest. `load_tailing_report` additionally rejects unknown or silently normalized nested JSON fields before replay. Generic serde decoding by itself is not an original-source preservation contract. An imported report is not physical authority; the design verifier must reconstruct the originating graph independently.

## TA and GC interpretation

The [manufacturer's TA Cloning manual, MAN0000008 revision 4](https://documents.thermofisher.com/TFS-Assets/LSG/manuals/ta_man.pdf) describes single 3′ A insert tails and complementary 3′ T vector tails. The independently implemented geometry supports that declared state; it does not infer completion from an enzyme name.

[Lucigen's manufacturer-authored pSMART GC manual, MA034 revision C](https://cdn.primeweb.fi/www.oligomer.fi/products/Product-Manual-(40736-1).pdf), distributed by Oligomer, describes G-tailed inserts and C-tailed vectors. That manual also explicitly specifies dephosphorylated vector ends and phosphorylated inserts. Generic base tailing therefore preserves the declared 5′ state rather than guessing kit activation or claiming a fully sealed product.

Both junctions matter. In the literal fixture, blunt `ACGC` vector plus blunt `GGTA` insert become the circular top reference `ACGCTGGTAA` after T/A tailing and compatible ordered ligation: ten bases, not the original eight. The other insert orientation gives `ACGCTTACCA`. C/G tailing analogously gives `ACGCCGGTAG`. A dephosphorylated vector with phosphorylated insert leaves explicit nicks under the existing ligation model; a sealed product requires all necessary phosphates. Two 5′-OH-only substrates cannot form a covalent ligation junction merely because their tails pair.

These symmetric tails do not establish orientation or a unique mixture product. Repeated-copy alternation, incomplete tailing, contaminants, off-pathway products, reagent conditions and experimental success are not exhaustively predicted. TOPO protein–DNA adducts and their bond-transfer mechanism are separate unsupported chemistry, not ordinary tailing or ligation.

## Bounds and evidence

The serialized report is limited to 32 MiB. Before cloning source/product buffers, a bounded streaming serializer counts source bytes and conservatively prepays two source-sized states plus 4 KiB of fixed report overhead. Both strands must also remain within the shared ten-million-residue bound. Output serialization checks the actual final limit; oversized inputs reject rather than produce partial reports. This is an output/work bound, not a claim about peak allocator memory or elapsed laboratory time.

Twelve focused native tests include 1,360 independently checked short source/tail combinations, manual paired-coordinate and ancestry coverage, reversal with nicks, literal two-junction TA/GC and phosphate-sensitive products, invalid-source/chemistry cases, source/report tamper rejection and preallocation rejection for an 8.4 Mb source. No cloning application's implementation or output is used as a correctness oracle. References were reviewed on 2026-11-24; no external scientific program is invoked.

## Application workflow

The **Cloning methods** inspector retains its source settings while a candidate is reviewed in **Plan**. It supports a restriction fragment (optionally explicitly blunted) or a unique exact PCR product, a declared single nucleotide on each actual 3′ end, separate terminal phosphorylation/dephosphorylation, physical fragment orientation, topology and nick policy. Opening an intermediate as a reference preserves evidence but never grants physical authority to a reimported molecule. Candidate construction always replays the preparation graph from the original immutable sources and binds the independently selected active target.

Standalone tail previews and their JSON downloads are source/recipe evidence, not target verification. Current source, recipe, workspace and target changes cancel pending work and invalidate displayed evidence. HTML/PDF protocols reverify the same ordered operations, oligos and declared chemistry; the new tail bases are explicitly nontemplated. The runtime invokes no external scientific program.
