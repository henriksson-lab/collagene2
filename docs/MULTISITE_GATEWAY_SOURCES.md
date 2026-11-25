# MultiSite Gateway primary sequence acquisition

Acquired2026-11-25. Files and derived manifest are in `/tmp/collagene2-multisite-primary/`; no repository/production edits, compiler or application/browser checks. Public manufacturer HTTP reads only. Local rendering inspected a manufacturer PDF, not the application. No third-party cloning code or vector sequences were borrowed.

## Result

Exact primary sequences now supply R3/R4 (including the correct VectorII edition), P3 and the numbered3/4 crossover contexts. A different manufacturer source supplies a fully printed L4/R1 entry-vector context. **A complete independently sequenced R2/L3 entry plasmid has not been located.** Its manufacturer MCS exposes a literal right flank, but that must not be relabelled as a fully resolved family-specific L3 site by filling adjacent insert Ns.

`acquisition-manifest.json` records every raw file URL/size/SHA256 and all extracted context sequences,1-based inclusive intervals, canonical orientation and core offset. The TXT source files contain literal sequence plus whitespace; stripping whitespace yields unambiguous ACGT of the reported lengths. PDF map annotations were read independently from the source sequence extraction.

## Downloaded exact manufacturer plasmid sequences

| Plasmid | Length | Sequence and annotation map |
|---|---:|---|
|pDONR P2R-P3|4773|[sequence](https://documents.thermofisher.com/TFS-Assets/vectors/pdonrp2rp3_seq.txt),[map](https://documents.thermofisher.com/TFS-Assets/vectors/pdonrp2rp3_map.pdf)|
|pDEST R4-R3 **VectorII**|4555|[sequence](https://documents.thermofisher.com/TFS-Assets/vectors/pdestr4-r3_vector_ii_seq.txt),[map](https://documents.thermofisher.com/TFS-Assets/vectors/pdestr4-r3_vector_ii_map.pdf)|
|pDEST R4-R3 **original**|4107|[sequence](https://documents.thermofisher.com/TFS-Assets/vectors/pdestr4r3_seq.txt),[map](https://documents.thermofisher.com/TFS-Assets/vectors/pdestr4r3_map.pdf)|
|pMS/GW control|5898|[sequence](https://documents.thermofisher.com/TFS-Assets/vectors/pmsgw_seq.txt),[map](https://documents.thermofisher.com/TFS-Assets/vectors/pmsgw_map.pdf)|

Do not substitute the4107bp original vector for MAN0000299's4555bp VectorII. Both contain matching125bp R3/R4 context strings, but their complete reference DNAs differ. VectorII moves the annotated R3 interval by448bp; exact event/source identity must remain distinct.

Raw sequence SHA256:

- `pdonrp2rp3_seq.txt`: `16b796fec73d1e06b38f5feeec059aa20f7e149d503ca75e0fc605222c1b10d9`
- `pdestr4-r3_vector_ii_seq.txt`: `e91c6ae899edcfd029be1c796177f4e09466eb50e36a346edc239228858da3ea`
- `pdestr4r3_seq.txt`: `7ecd76ddc3d7c005bc43e0118b1b3fd9c146d528991f54c4d5f48d1f2ebbbf8f`
- `pmsgw_seq.txt`: `0b5248b5592a779b1aa98539f99944adb7ba4e3a91450a541112262e4e0a6b6d`

## Literal source intervals and orientation

“RC” below means reverse-complement extraction to read the canonical seven-base specificity overlap. It is determined from literal DNA, not inferred from an annotation strand or the letterr in a vector name.

| Source | Context interval,1-based inclusive | Length | Canonical extraction | Core / zero-based offset |
|---|---|---:|---|---|
|pDONR P2R-P3|P2R591–822|232|RC|GTACAAG/82|
|pDONR P2R-P3|P3 2746–2977|232|RC|GTATAAT/82|
|pDEST R4-R3 VectorII|R4 37–161|125|forward|GTATAGA/7|
|pDEST R4-R3 VectorII|R3 2064–2188|125|RC|GTATAAT/7|
|pDEST R4-R3 original|R4 37–161|125|forward|GTATAGA/7|
|pDEST R4-R3 original|R3 1616–1740|125|RC|GTATAAT/7|
|pMS/GW|B4 37–57|21|forward|GTATAGA/7|
|pMS/GW|B1 1285–1308|24|forward|GTACAAA/7|
|pMS/GW|B2 3154–3174|21|RC|GTACAAG/7|
|pMS/GW|B3 3510–3530|21|RC|GTATAAT/7|

The literal B contexts are21/24bp here, not the earlier profile's25bp primer contexts. Do not pad them or silently declare them identical variants. VectorII's map does not print(c) next to R3 even though canonical overlap orientation is RC; raw sequence determines the extraction. These contexts are versioned full *manufacturer-annotated spans*, not claimed universal minimal recognition sequences.

## Primary L4/R1 source: pENTR5′-TOPO

The [official map](https://documents.thermofisher.com/TFS-Assets/vectors/pentr5ptopo_map.pdf) identifies the2680nt pENTR5′-TOPO reference:

- attL4:592–688,97nt,forward canonical GTATAGA at offset83;
- attR1:721–845,125nt,forward canonical GTACAAA at offset8.

The [official detailed MCS](https://documents.thermofisher.com/TFS-Assets/vectors/pentr5ptopo_mcs.pdf) prints **all bases within both intervals**. Neither interval crosses its displayed PCR-insert/TOPO gap. Both literal intervals were transcribed from printed row coordinates and visually checked. File outputs `pentr5ptopo-L4-literal.txt` and `pentr5ptopo-R1-literal.txt` plus manifest contain them. Source is a published partial reference diagram; a complete2680nt raw sequence has NOT been retrieved, and none was reconstructed by filling the gap.

Useful coordinate transcription checks: displayed left rows start501 and581 with80bases each, then661; right rows start718 (73bases to790) and791 (80bases to870). This matters because spacing expands codons and homopolymers. The resulting literal L4 is97bases, not padded to100. R1 has a different arm context from the existing pDEST14 profile; retain its own variant identity.

Raw PDF SHA256:

- `pentr5ptopo_map.pdf`: `2d3b121b6a93ca4c07891c6c1c3de0c69bd56d23d440144c6d819330d02c45f3`
- `pentr5ptopo_mcs.pdf`: `485adcb080bacad812a3fd053090751d53b4a28b8dafa8da058d7204fb0f7182`

## L3 / P4 remaining gap and useful primary evidence

Downloaded [pDONR P4-P1R map](https://documents.thermofisher.com/TFS-Assets/vectors/pdonrp4p1r_map.pdf) identifies4777nt donor, P4 593–824(c), P1R2748–2979(c). Its [MCS](https://documents.thermofisher.com/TFS-Assets/vectors/pdonrp4p1r_mcs.pdf) illustrates L4/R1 entry-flank sequence with an explicit unknown insert; it is not a full donor/entry sequence file. Tested plausible sequence URLs returned403, which may indicate nonexistent assets or access behavior; it does not prove sequence confidentiality.

Downloaded [pDONR P2R-P3 MCS](https://documents.thermofisher.com/TFS-Assets/vectors/pdonrp2rp3_mcs.pdf) illustrates R2/L3 entry-flanks and explicitly shows an unknown3′element. Its post-insert flank labelled2884–2977 is94literal bases:

`ACTTTATTATACAAAGTTGGCATTATAAAAAAGCATTGCTTATCAATTTGTTGCAACGAACAGGTCACTATCAGTCAAAATAAAATCATTATTT`

This94nt suffix is **also exactly present in the full P3 donor sequence at2884–2977**. Consequently, treating it alone as a distinguishing full attL3 profile would also match the unrecombined attP3 donor. The missing upstream family-discriminating context must be resolved from a real entry sequence or an explicitly source-derived BP model with honest computed provenance. It must not be filled from insert Ns or labelled an independently published full L3 site. A raw pENTR R2/L3 plasmid source remains unresolved in this acquisition.

The manual's donor/PCR-arm diagrams can guide the expected site-family relationships, but do not authorize inventing unprinted arms or pretending a computed hybrid is a primary sequence. The existing direct L1/L2 source remains the earlier full pENTR-D-TOPO vector asset; these variants should not be replaced with presumed equivalent arms.

## Reproducible discovery path / catalog caveat

The current public product app exposes a read-only document search:

`https://www.thermofisher.com/api/store/Assets/Documents/Others/search?skus=12537023&country=us&targetSite=TF&partialLotNumber=true`

and the same endpoint forK59120. Both JSON responses are saved. Their vector `path` fields provided the exact VectorII/pENTR5′MCS filenames. **Several catalogue titles are visibly misindexed**: e.g.the VectorII path is returned under a GENEART control-insert title; pDONR P2R-P3 under pThioHis. Verify actual PDF label/length and literal DNA against the associated map, not catalogue display title. No account, write operation or private endpoint was used.

Original manuals are retained with raw hashes:

- [MAN0000299](https://documents.thermofisher.com/TFS-Assets/LSG/manuals/multisite_gateway_man.pdf): `20b9d15a7b53d94066bdc94064318ebba9ae1fd0db73ec0512ce65c4d5fd98a9`.
- [pENTR5′-TOPO manual](https://documents.thermofisher.com/TFS-Assets/LSG/manuals/pentr_5primetopo_man.pdf): `3d58e1ba75cb79040b5ba1a1ef8260e7ef3a9e4323e5e9e7ed3de954decdfcca`.

This is source acquisition only. It does not grant a new profile biological authority, establish reagent performance or license repackaging all manufacturer documents into MIT code. No data were committed or bundled. Next implementation can use the fully resolved contexts while leaving L3/full-entry recognition explicitly unresolved until its distinguishing bases and provenance are supplied.
