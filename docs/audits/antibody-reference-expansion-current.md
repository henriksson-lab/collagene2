# Antibody reference expansion: concrete source audit

Research date: 2026-11-28. Antibody GUI/workflow support remains in scope. This audit supplies reference-data inputs; it does not claim broader detection, species classification, numbering certification or biological calibration is implemented. No Rust/Cargo changes, builds, tests, external scientific executables or Collagene application checks were performed.

## Current owner and smallest coherent next batch

The current [reference owner](../../crates/antibody-numbering/src/references.rs#L1) admits exactly three human OGRDB sources: G00060/v10 IGH, G00061/v5 IGK, G00062/v4 IGL. It retains 354 functional V-prefixes and 23 J/FR4 segments with a CC0 manifest. V references have explicit 104-column provider labels; J positions beginning at IMGT118 are an authored complete-codon rule, not provider-supplied J labels. [Current attribution and source URLs](../../crates/antibody-numbering/data/ogrdb-human-vj/ATTRIBUTION.md).

**Recommended next batch: mouse C57BL/6 IGH only, from two paired CC0 OGRDB v5 representations.** This adds one species/strain and reuses the Heavy locus without pretending to support every mouse strain, light chain, TCR or numbering convention. Begin with supported IMGT correspondence; retain unresolved mappings. Kabat/Chothia/North applicability needs separate chain/species evidence.

| Pinned acquisition | Observed payload and qualification |
| --- | --- |
| [C57BL/6 IGH JSON v5](https://ogrdb.airr-community.org/api_v2/germline/set/10090.C57BL.sl.6%20IGH.C57BL.sl.6/5) | `OGRDB:G00004`, `OGRDB:C57BL/6 IGH.5`, release 2024-10-15; 114 alleles: 102 V, eight D, four J. Of the V records, 95 have provider `functional: true`. These counts are not validated numbering-model counts. |
| [Same v5 provider-gapped FASTA](https://ogrdb.airr-community.org/api_v2/germline/set/10090.C57BL.sl.6%20IGH.C57BL.sl.6/5/gapped) | 35,235 bytes, 114 records. The [official API guide](https://github.com/airr-community/ogrdb/blob/master/schema/ogrdb_api_v2_guide.md) calls this format IMGT-gapped. Its labels differ from JSON, e.g. `IGHV0-24BS*00` versus `IGHV-24BS`; preserve both, not an invented alias. Independent literal comparison found every ungapped FASTA sequence matches exactly one JSON `coding_sequence` in this pinned pair. Require that checked correspondence during derivation; reject missing/ambiguous joins. |

Important: all 102 V JSON delineations have null `aligned_sequence` and `alignment_labels`. Their region coordinates alone do **not** satisfy the current human corpus's full column-label contract. The paired provider-gapped representation is additional evidence, not a reason to treat unlabeled JSON as already numbered. Validate frame, whole-codon gaps, prefix length, translation, region boundaries and anchors before admitting individual models. Keep exclusions and atypical anchors explicit.

## License findings and preserved acquisition

[OGRDB's content license](https://ogrdb.airr-community.org/render_page/licensing_statement.html) is CC0; its website code is separately EUPL, and no code is borrowed. This supplies an unambiguous data route for the proposed mouse batch.

IMGT's **current** [official terms](https://www.imgt.org/about/termsofuse.php) explicitly grant data **and metadata** under CC BY4.0 to public/private users; tools are separate. However, the exact [V-QUEST download README](https://www.imgt.org/download/V-QUEST/README.txt) still says academic research only/contact for other use, and the [GENE-DB README](https://www.imgt.org/download/GENE-DB/README.txt) mentions private-company agreements. These older file-specific notices conflict with the current general policy. This audit records the conflict without resolving it as legal advice. **No IMGT raw or derived corpus is admitted into the MIT crate here.** Retain attribution, modifications and CC BY notices separately if the conflict is subsequently assessed in favor of admission; never relabel source data MIT.

Sixteen original research payloads, format/license captures and API enumerations are retained at `/tmp/collagene2-antibody-reference-research-2026-11-28/`. Its `manifest.json` records exact requested URL, capture timestamp, byte count, SHA-256 and per-file license-policy/capture linkage. It includes both mouse JSON versions inspected, the chosen v5 gapped file, alpaca V/J, human TRAV and actual OGRDB species/set/version responses. This directory is research retention, not a repository-distributed dataset.

## Broader source roster and exact access

Live [OGRDB species enumeration](https://ogrdb.airr-community.org/api_v2/germline/species) returned human, mouse, rhesus, *Peromyscus* and rainbow trout. It did not list alpaca, rat, rabbit or pig. [Mouse sets](https://ogrdb.airr-community.org/api_v2/germline/sets/10090) contains 45 strain-specific sets; do not flatten them into universal mouse coverage. [Human sets](https://ogrdb.airr-community.org/api_v2/germline/sets/9606) lists IGHC plus the three IG V/J sets, with **no current TCR set**. [Rhesus sets](https://ogrdb.airr-community.org/api_v2/germline/sets/9544) lists `9544.IGH_VDJ`, `9544.IGK_VJ`, `9544.IGL_VJ`; [IGH versions](https://ogrdb.airr-community.org/api_v2/germline/set/9544.IGH_VDJ/versions) returned 1 and 2. Full rhesus payloads were not decoded in this audit; no exact column-label claim is made. The [provider release announcement](https://wordpress.vdjbase.org/index.php/ogrdb/rhesus-macaque-ig-germline-sets-released/) describes population-derived variation and ASC naming: this is not a one-animal complete germline or a classifier.

IMGT's [V-QUEST directory page](https://www.imgt.org/vquest/refseqh.html) currently identifies release **202643-4, 22 October 2026**. Direct file URLs are mutable; pin each downloaded response digest and captured release page instead of assuming the URL is immutable. The following existing directory listings expose actual files, subject to the license conflict above:

| Scope | Authoritative download directory and files |
| --- | --- |
| Mouse IG | [Mus_musculus/IG](https://www.imgt.org/download/V-QUEST/IMGT_V-QUEST_reference_directory/Mus_musculus/IG/): `IGHV/IGHJ`, `IGKV/IGKJ`, `IGLV/IGLJ` as `.fasta`. The generic directory can include multiple mouse taxa/strains; preserve each record's species rather than rewrite it. |
| Alpaca | [Vicugna_pacos/IG](https://www.imgt.org/download/V-QUEST/IMGT_V-QUEST_reference_directory/Vicugna_pacos/IG/): [IGHV.fasta](https://www.imgt.org/download/V-QUEST/IMGT_V-QUEST_reference_directory/Vicugna_pacos/IG/IGHV.fasta), [IGHJ.fasta](https://www.imgt.org/download/V-QUEST/IMGT_V-QUEST_reference_directory/Vicugna_pacos/IG/IGHJ.fasta). The listing has heavy V/D/J only, not alpaca kappa/lambda. |
| Human TCR | [Homo_sapiens/TR](https://www.imgt.org/download/V-QUEST/IMGT_V-QUEST_reference_directory/Homo_sapiens/TR/): `TRAV/TRAJ`, `TRBV/TRBJ`, `TRGV/TRGJ`, `TRDV/TRDJ` `.fasta`; D files also exist. |
| Mouse TCR | [Mus_musculus/TR](https://www.imgt.org/download/V-QUEST/IMGT_V-QUEST_reference_directory/Mus_musculus/TR/): the same four V/J pairs and applicable D files. |
| Rat | [Rattus_norvegicus/IG](https://www.imgt.org/download/V-QUEST/IMGT_V-QUEST_reference_directory/Rattus_norvegicus/IG/): H/K/L V/J pairs. |
| Rabbit | [Oryctolagus_cuniculus/IG](https://www.imgt.org/download/V-QUEST/IMGT_V-QUEST_reference_directory/Oryctolagus_cuniculus/IG/): H/K/L V/J pairs. |
| Rhesus | [Macaca_mulatta/IG](https://www.imgt.org/download/V-QUEST/IMGT_V-QUEST_reference_directory/Macaca_mulatta/IG/): H/K/L V/J pairs; use independently audited CC0 OGRDB if suitable. |
| Pig | [Sus_scrofa/IG](https://www.imgt.org/download/V-QUEST/IMGT_V-QUEST_reference_directory/Sus_scrofa/IG/): H/K/L V/J pairs. The indexed listing had an older per-file date than the global release page: another reason to hash actual bytes. |

The [GENE-DB bulk index](https://www.imgt.org/download/GENE-DB/) additionally exposes [gapped amino-acid references](https://www.imgt.org/download/GENE-DB/IMGTGENEDB-ReferenceSequences.fasta-AA-WithGaps-F%2BORF%2BinframeP), gapped nucleotide references and `RELEASE`. Its README documents V/C gapping and 15-field headers. These can inform constant-domain work, but a complete constant chain/exon sequence is not automatically a delineated constant-domain ground truth. Exact domain boundaries, exon phase and numbering must be independently retained and validated. Bulk payload and release value were not acquired here.

## Format checks and literal anchors

The [V-QUEST format README](https://www.imgt.org/download/V-QUEST/README.txt) explicitly identifies IMGT gaps in V-REGION and describes the 15 header fields, including functionality, region, nucleotide codon start, sequence edits, partial and reverse-complement indicators. Do not discard ORF/P, partial or corrected-source flags. `rev-compl` describes provenance; sequences are already supplied in gene coding orientation, not a request to complement them a second time. A J header's codon start does not supply CDR3/FR4 position labels.

Independently decoded literal checks, usable as future curation fixtures:

- Mouse paired reference `IGHV0-24BS*00` joins uniquely to `OGRDB:Mus musculus/C57BL/6_IGH:IGHV-24BS.1`. Gapped codon columns 23/41/104 are `TGC/TGG/TGT`, translating to C/W/C. This is a reference check, not a query-domain classifier.
- Downloaded alpaca IGHV has 84 records (73 `F`, ten `P`, one `ORF`), 35,027 bytes, SHA-256 `1008c346a4c3f8c646b843598059d7f2d9a631727014f82b3c9a710b5e403453`. `AM773729|IGHV1-1*01` declares codon start1; the same aligned columns are `TGC/TGG/TGT`. These 73 F headers are **not 73 admitted or validated models**.
- Alpaca IGHJ has seven records (six F, one ORF), 961 bytes, SHA-256 `bfbfe40904f033eb82b3efd1ee907ac311af86e9b020cfdf0f72263e9ff50145`. `AM773729|IGHJ2*01` declares codon start3. Starting at raw zero-based offset2 gives whole codons; codons7–10 are `TGG GGC CAG GGC` = W-G-Q-G. That motif is a checked candidate anchor, not provider-supplied IMGT118–121 labels or a complete rearranged CDR3.
- Human TRAV acquisition has 136 records, 56,861 bytes; headers include `F`, `(F)`, `[F]` and P. Do not flatten bracketed functionality into the human corpus's boolean without documenting its provider meaning.

## Necessary owner/algorithm changes

1. Replace the human-only corpus assumption with one current source-set selector and species/strain/locus-scoped manifest. Update [hard-coded source counts/labels](../../crates/antibody-numbering/src/references.rs#L78), [human profile construction](../../crates/antibody-numbering/src/v_mapping.rs#L50), candidate binding and UI/MCP selectors directly. Preserve exact evidence/provider labels and per-source licenses; no compatibility loader.
2. Keep provider-gapped columns distinct from an authored gap/label interpretation and from query correspondence. Audit every proposed functional V prefix/J frame, incomplete codon and anchor; retain excluded records with reasons. Do not train or certify a model merely by including its germline sequence.
3. A mouse Heavy-only batch does not need new receptor-locus variants. Later TCR work needs explicit Alpha/Beta/Gamma/Delta loci, locus-specific V/J anchors and FR4 handling, dual-use gene identity treatment, IMGT-only scheme gating, and new region/profile validation. Existing Heavy-versus-light branches in [domain tracks](../../crates/antibody-numbering/src/domain_tracks.rs#L68) cannot be reused as generic TCR rules.
4. Alpaca IGH alone does not establish VHH versus conventional heavy-chain subtype, nanobody specificity, a constant domain or whole-antibody topology. Missing light/constant evidence remains missing. Broader species scores must not be compared as calibrated likelihoods without independent evaluation.
5. Add source-specific literal/frame/gap fixtures and independent natural-domain fixtures before claims of broader numbering correctness. Reference germline facts, CDR-definition rules, query-domain detection and experimental antibody behavior remain separate evidence layers. GUI display/menu comparison remains in [the Benchling audit](benchling-antibody-gui-current.md).
