# Saved primer reuse with an explicit 5′ extension

Open `prepared-overlap-demo.fa`, then import `primer-reuse-demo.csv` as an oligo
table. Map name, ordered sequence, annealing domain and tail to their separately
named columns. Select both rows to open the two oligo documents.

Make `left_template` active and choose **PCR**. Expand saved forward-primer reuse,
explicitly select the `saved-forward` source document and copy its primer. Repeat
for the reverse primer from `saved-reverse`. Review the editable names and domains:

| Primer | Annealing, 5′→3′ | 5′ extension | Complete ordered oligo, 5′→3′ |
| --- | --- | --- | --- |
| saved-forward | TTGCAATCGT | GG | GGTTGCAATCGT |
| saved-reverse | GCTAGCATCG | empty | GCTAGCATCG |

Exact PCR on `TTGCAATCGTCGATGCTAGC` yields the independently specified reference
preview `GGTTGCAATCGTCGATGCTAGC` (22 bases). The two forward extension bases are
introduced, not copied from template coordinates. The preview should retain the
exact template event as its parent; saved primer names and actual ordered
sequences remain part of its PCR evidence. This is a reference preview, not
verified physical chemistry or laboratory success.

The same explicit saved-primer picker is available in exact-PCR assembly
preparation. Copying is a snapshot of sequence fields: manual edits remain
possible, and changing preparation fields invalidates previous assembly evidence.
The `GG` extension changes the independently required target; it must not silently
verify against the unextended target in the earlier assembly example.

If a selected primer document changes before copying, select its new revision
and review it again. Missing, conflicting or unsupported sequence qualifiers stay
visible as unavailable choices. An annotation label or imported binding claim
does not establish amplification or specificity.
