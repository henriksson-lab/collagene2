# Three-piece ligation with explicit PCR chemistry

These deliberately short synthetic sequences demonstrate strand/ligation
accounting, not practical primer design. Open `ligation-demo.fa`, make
`ligation-target` active, choose **Assemble**, then the ordered-ligation workflow.
Choose a **linear** product and add three rows in the following order:

| Original source | Preparation | Forward annealing, 5′→3′ | Reverse annealing, 5′→3′ | Treatment |
| --- | --- | --- | --- | --- |
| ligation-left | Exact PCR | AAA | GGG | Phosphorylate |
| ligation-middle | Exact PCR | GGG | AAA | Phosphorylate |
| ligation-right | Exact PCR | ATG | ACG | Phosphorylate |

Leave tails empty and orientations forward. Require sealed junctions, then
compute evidence. The independently expected top strand is
`AAACCCGGGTTTATGCGT` (18 bases). Each exact-PCR fragment starts blunt with 5′-OH /
3′-OH termini under the existing declared model; phosphorylation is a separate
step. Inspect the actual termini, two internal junctions, absence of retained
nicks, and the independently supplied target assessment.

Remove phosphorylation from the left fragment: sealed-only ligation must fail.
Explicitly allowing retained nicks permits one unsealed bottom-strand bond at the
first junction while its other strand seals. The source preparation, missing
phosphorylation and nick remain evidence; no host repair is assumed. Keeping all
fragments unphosphorylated cannot be rescued merely by enabling nicks, because
neither strand can seal at a junction.

Switching to `ligation-wrong-target` must invalidate old evidence and a fresh
check must reject the changed final base. Reversing a row changes its actual
oriented strands and the expected product. Candidate review and HTML/PDF export
retain all three original source roles, PCR primers, terminal treatments, order
and declared nick policy.

An imported linear DNA reference alone has unknown termini. Selecting
whole-linear input does not declare 3′-OH or infer physical chemistry from a file
annotation. Use explicit supported preparation; phosphorylation by itself does
not resolve unknown 3′ chemistry. Results establish only the declared ordered
candidate, not uniqueness among a mixture or experimental success.
