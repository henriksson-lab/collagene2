# Finite three-entry source-context facts

Independently curated 2026-11-25; profile edition `2026-11-25.1`, ID
`tol2kit-228-383-302-vectorii-lr-reference`. The JSON contains nine functional
attachment contexts, 1,188 DNA bases in total. It contains no whole plasmid,
GenBank annotation record, manufacturer manual text or copied implementation.
Our implementation and curation code use the repository's MIT license; source
citations do not assert manufacturer endorsement or relicense source documents.

This is a finite reference-matching profile for the cited source editions,
not a universal Gateway/Tol2kit profile, a physical strand mechanism or predicted
kit success. Exact full context and canonical seven-base overlap are retained;
missing distal arms are never invented. The known full P3 competitor is necessary:
its first96 canonical bases equal the complete L3 context. Classify all registered
families before selecting LR-reactive contexts. Unregistered variants remain
outside the declared model.

## Context coordinates and orientation

Coordinates are one-based inclusive in the original source reference. “RC” means
reverse-complement the extracted interval to obtain the canonical JSON sequence.
Overlap offsets are zero-based within that canonical sequence.

| Source edition | Context | Interval | Orientation | Bases | Overlap offset / canonical core |
| --- | --- | --- | --- | ---: | --- |
| p5E-MCS clone228,2810bp | L4 |593–688|Forward|96|82 / GTATAGA|
| Same | R1 |856–1012|Forward|157|7 / GTACAAA|
| pME-EGFP clone383,3327bp | L1 |569–668|Forward|100|83 / GTACAAA|
| Same | L2 |1452–1551|RC|100|83 / GTACAAG|
| p3E-polyA clone302,2838bp | R2 |591–747|RC|157|7 / GTACAAG|
| Same | L3 |947–1042|RC|96|82 / GTATAAT|
| pDEST R4-R3 VectorII,4555bp | R4 |37–161|Forward|125|7 / GTATAGA|
| Same | R3 |2064–2188|RC|125|7 / GTATAAT|
| pDONR P2R-P3,4773bp | P3 competitor |2746–2977|RC|232|82 / GTATAAT|

## Primary source provenance

The original depositor lab supplies [p5E-MCS](https://tol2kitkwan.genetics.utah.edu/index.php/P5E-MCS),
[pME-EGFP](https://tol2kitkwan.genetics.utah.edu/index.php/PME-EGFP) and
[p3E-polyA](https://tol2kitkwan.genetics.utah.edu/index.php/P3E-polyA).
The JSON `url`/`raw_sha256` pair identifies the literal downloaded ZIP; its
`map_url`/`map_sha256` identifies the source page linking that GenBank.
Those public pages are provenance records, not claims of enzyme behavior.

| Original archive member | Member SHA256 |
| --- | --- |
| `228. p5E-MCS.gb` | `b0fc54681b822e36ad17bc50148fbf057c2893ad9098d3e7c2e1bf54c4780889` |
| `383. pME-EGFP.gb` | `696622988621fa7246ec4c67dd8d249fc78058fff0d4ed5c7d14da5145a25e9d` |
| `302. p3E-polyA.gb` | `58627c3e0fa0c054faedc219db51fa924aa363bb4b1d96dddb7ca359ca2f81aa` |

The manufacturer supplies the [VectorII sequence](https://documents.thermofisher.com/TFS-Assets/vectors/pdestr4-r3_vector_ii_seq.txt)
and [map](https://documents.thermofisher.com/TFS-Assets/vectors/pdestr4-r3_vector_ii_map.pdf),
and the [P2R-P3 donor sequence](https://documents.thermofisher.com/TFS-Assets/vectors/pdonrp2rp3_seq.txt)
and [map](https://documents.thermofisher.com/TFS-Assets/vectors/pdonrp2rp3_map.pdf).
Their JSON hashes identify exact TXT/PDF bytes, not a normalized DNA digest.
The [three-fragment manual](https://documents.thermofisher.com/TFS-Assets/LSG/manuals/multisite_gateway_man.pdf)
documents the L4–R1 / L1–L2 / R2–L3 / R4–R3 arrangement; no manual prose or figures
are incorporated into the profile.

Canonical orientations for the lab contexts were checked from the complete source
DNA and shared cores rather than trusting unqualified annotation strand syntax.
No source methylation comment establishes inherited chemical state. The numerical
intervals, sequences, hashes and citations are factual context observations;
profile selection and reference reconstruction remain our explicit bounded model.
