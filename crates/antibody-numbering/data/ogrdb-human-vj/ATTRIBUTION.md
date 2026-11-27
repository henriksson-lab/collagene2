# Human IMGT V-prefix and J/FR4 reference corpus

The data in `references.json` derives from **Open Germline Receptor Database
(OGRDB)** releases G00060 version 10 (human IGH), G00061 version 5 (human IGK),
and G00062 version 4 (human IGL), released 27 May 2026 and retrieved
27 November 2026. Curator: William Lees; AIRR Community germline working groups.

OGRDB publishes its content under **CC0-1.0**. This data dedication is separate
from Collagene 2's MIT software license; no OGRDB implementation code is included.
[Provider license statement](https://ogrdb.airr-community.org/render_page/licensing_statement.html)
· [CC0 legal text](https://creativecommons.org/publicdomain/zero/1.0/legalcode)
· [Provider citation information](https://ogrdb.airr-community.org/help).

This is a Collagene-derived corpus, not an unmodified provider release or a trained
numbering/classification model. Each source URL, exact original-response SHA-256,
release identity, retained count and excluded allele reference is in the manifest.
The three API responses are the original inputs; unrelated source metadata and
supporting genomic contigs are not embedded. Source responses are about 4.4 MB
in total; the derived corpus is below 0.5 MB.

Provider-annotated functional V alleles are retained: 194 IGH, 74 IGK and 86 IGL.
Functional J alleles are also retained: 7 IGH, 7 IGK and 9 IGL. All 44 other records
have an explicit exclusion ledger. For each retained allele,
the first 104 explicitly labeled IMGT codon columns, their nucleotide sequence,
translated residues, provider allele identity and CDR1/2 delineations are retained.
Gaps remain gaps. Provider-annotated functional alleles with atypical conserved
cysteines are retained with an explicit flag; no residue is repaired or relabeled.

Germline CDR3 tails beyond the provider's labeled columns are excluded. The corpus
does not supply rearranged CDR3, query-domain assignments, Kabat/Chothia mappings, nonhuman references
or a calibrated domain classifier. Reference positions do not automatically
establish the corresponding positions in a query protein.

Rebuild from the three **original** API JSON files, in this order:

```sh
python3 scripts/build_antibody_references.py \
  G00060-v10.json G00061-v5.json G00062-v4.json \
  --output crates/antibody-numbering/data/ogrdb-human-vj/references.json
```

The generator verifies all original hashes and source identities before producing
output. It checks nucleotide prefix agreement, complete-codon gaps, explicit label
columns, standard-code translation and CDR boundary agreement. It performs no
network or external scientific-tool calls. A native integration fixture compares
every retained translation with Collagene's existing Rust translation engine.
[Numbering/delineation format](https://docs.airr-community.org/en/latest/datarep/germline.html).


The current derived corpus format is version 2, replacing the V-only development
format directly. It has separate `v_entries` and `j_entries`; no old-format reader
or migration is retained. Original provider release identities are preserved as
source provenance, independently of this derived internal version.

For each J allele, original sequence, coding interval, coding DNA and raw provider
frame/anchor fields are retained. In these pinned responses the anchor locates
the conserved W/F codon within coding DNA, and its phase agrees with the provider
complete-codon start. The generator checks these identities rather than applying
a generic interpretation of AIRR coordinates. FR4 proteins contain complete
codons only. One terminal nucleotide is preserved separately and remains
untranslated; it contributes to the splice-junction codon of the constant region.

J position labels are authored from the published IMGT conserved-anchor rule:
118 followed sequentially through complete J codons. They are **derived rule
labels**, unlike the explicit provider V-column labels. The pinned human heavy
FR4 proteins have 11 residues (118–128), light FR4 proteins 10 (118–127).
[FR4/anchor definition](https://www.imgt.org/IMGTeducation/Tutorials/3Dstructure/_UK/DOMAINs/fleches.html)
· [Coding-region/splice boundary](https://imgt.org/IMGTScientificChart/Numbering/IMGTnumberingCR.html).
