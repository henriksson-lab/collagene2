# iGEM RFC1000 Level-0 factual profile

Primary source: [iGEM Foundation, DNA Assembly](https://technology.igem.org/technologies/assembly/),
“iGEM Type IIS Assembly” technical specification, Assembly Compatibility and
Basic Parts (Level 0s), Fusion Sites. Reviewed 2026-11-28. The page supplies no
immutable revision ID. Search retrieval exposed the primary text; direct opening
returned the site's JavaScript shell. No original-byte digest is claimed.

This directory contains independently authored MIT factual curation and code,
not a copy of upstream prose, images, HTML or software. No upstream prose/assets
license was verified; MIT is **not** asserted for those materials. The short DNA
motifs, role associations and declared overlap facts are independently recorded;
retain this source attribution when distributing the curation. Nothing asserts
endorsement by iGEM or rights in particular biological constructs or materials.

| Declared Level-0 role | Left fusion | Right fusion |
|---|---|---|
| Promoter | GGAG | TACT |
| 5′ UTR / RBS | TACT | AATG |
| CDS | AATG | GCTT |
| Terminator | GCTT | CGCT |

Columns follow the assembled reference strand 5′→3′. They are not independently
oriented right-end physical overhang strands. Forbidden internal recognition
motifs are BsaI `GGTCTC` and SapI `GCTCTTC`, in either orientation. Additional
restriction-enzyme rules are not implied. The documented CDS begins with `ATG`;
these three bases overlap the left fusion site's final three bases.

Input is the documented body, with no flanks and explicit interval `[0,len)`.
Both profile additions are deliberately empty. The local grammar is our declared
**contiguous subpath** admission scope, not an assertion that any partial chain is
a complete functional transcription unit. The four-role table induces three
adjacent transitions; every role may be a subpath start/end and occurs at most
once. The profile omits vector architecture, spacers/clamps, other roles and
higher assembly levels rather than inventing requirements. Catalog geometry is
bound to the separately curated installed enzyme catalog edition/content hash.
