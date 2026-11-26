# Explicit 3′ resection and homology annealing

This native scientific foundation models **declared hydrolytic 3′ shortening and
exact pairing of exposed 5′ tails**. It is separate from the existing
`homology_assembly` 5′ resection/filling/ligation kernel. No existing physical
operation, In-Fusion profile, source authority, transport schema or user workflow
changes merely because this module exists.

## Primary source boundary

[Takara's In-Fusion planning guidance](https://help.takarabio.com/36168/kb/article/109886/how-should-i-plan-my-in-fusion-experiment)
describes 3′ exonuclease-generated 5′ overhangs, complementary terminal annealing
and rescue of the circular recombinant in *E. coli*. It explicitly excludes
covalent assembly of linear molecules. This supports the mechanistic distinction
from Gibson; it does not disclose the enzyme formulation, deterministic cut
positions, blocked-end processing or cellular repair outcomes.

[NEB's Exonuclease III product description](https://www.neb.com/en-us/products/m0206-exonuclease-iii-e-coli)
describes removal from accessible 3′ hydroxyl termini. Its
[activity chart](https://www.neb.com/en-us/tools-and-resources/selection-charts/properties-of-exonucleases-and-nonspecific-endonucleases)
and [protocol](https://www.neb.com/en/protocols/protocol-for-exonuclease-iii-m0206)
distinguish substrate geometry and note cut-length dependence on conditions.
These are primary examples for a declared hydrolytic mechanism, **not a claim
that Exonuclease III is the proprietary In-Fusion enzyme**. No enzyme-specific
kinetics, sequence inhibition, phosphorothioate susceptibility, AP activity,
3′ phosphatase or cellular rescue is imported into this model. Research baseline:
26 November 2026. Implementation is independent.

## Resection contract

`three_prime_resection::geometry(&Duplex, &Cuts)` takes explicit counts in the
supplied, already oriented source geometry:

- `left` removes the **bottom strand's stored 3′ suffix** at the physical left.
- `right` removes the **top strand's stored 3′ suffix** at the physical right.

The top 5′ reference origin remains zero. The bottom reference start advances by
`left`; both original 5′ termini and their P/OH/unknown groups remain unchanged.
Positive cuts require known accessible 3′OH. A new retained terminus after actual
hydrolytic cleavage has 3′OH. Zero cuts preserve the original 3′ group, including
blocked or unknown. There is no assumed unblocking, phosphorylation, filling,
ligation or repair. Source bases are removed, never synthesized.

The generic transform accepts an intact nick-free linear duplex, including
existing staggered ends, as long as neither strand nor its last paired core is
exhausted. It does not claim a particular enzyme can process every such
substrate. Circular/nicked substrates are rejected because additional accessible
ends or internal-nick processing would require a separate declaration.

The report retains exact source geometry revision, removed stored-strand spans,
actual resulting strands and exhaustive copy mappings. `replay` freshly checks
the full report against the original source and counts. Neither a report nor a
serialized geometry restores native material authority.

## Ordered annealing contract

`homology_annealing::Options` contains `overlaps`, `topology`, one explicit `Cuts`
per source, minimum paired length and comparison/candidate bounds. There are
1–32 ordered operands and `n-1` linear or `n` circular junctions, including
single-source processing or self-circularization. Initially intact **blunt**
operands are required by this assembly adapter. Nonblunt input must not be
silently flattened; its independent resection remains supported above.

`geometry`/`geometry_borrowed` have no target argument. They validate each selected
suffix/prefix homology directly from original source strands and require both
resections to expose the whole declared overlap. Cuts are not selected from a
target or predicted automatically. Every cut is explicit, including linear outer
ends; the model never silently spares them or restores blunt ends.

The report contains:

- Every actual resected operand with its copy mappings.
- Every separate surviving covalent strand, literal sequence, source indices,
  contour placement and exact 5′/3′ groups.
- Each homology and the unsealed top/bottom discontinuities. Equal boundaries
  mean a nick; distinct boundaries delimit an actual gap.
- A source-derived **coordinate reference**, explicitly distinguished from the
  extant strand pieces and any absent bases.
- Adjacent exposed-terminal alternatives and separate bounded-search completion.

The final circular bottom piece may extend past the contour origin. Its unwrapped
interval retains the actual span; modulo coordinates are a display operation,
not permission to erase a gap or change strand identity. Reversing a source
occurs before supplying it and before interpreting its two cuts.

For exact-H resection, paired gapless coverage fits a public diagnostic `Duplex`
with every joining nick retained. Original PCR 5′OH stays 5′OH; no joining bond
is sealed even when 5′ phosphate is present. That diagnostic is never returned
as an assembled `TrackedDuplex`. The current private nick representation cannot
preserve all nick chemistry, and serialization is not a material channel.
Extra resection leaves actual top/bottom gaps: `paired_diagnostic` is absent,
while all strand pieces and gap chemistry remain available. No polymerase is
silently added. Every joined assembly remains explicitly incomplete pending an
independently declared supported finishing/host model. Single-linear processing
with no junction is a different, completely specified transform.

## Provenance bridge and integration

`TrackedDuplex::resect_three_prime` reconstructs the actual processed operand
from freshly supplied tracked input via the existing source-piece mapper. Only
surviving modification calls and generation provenance remain; no new nucleotide
or default unmodified tract is invented. `modified_duplex::homology_anneal`
returns **the separate processed tracked operands plus the diagnostic report**,
not a fictitious covalently assembled molecule. Original source objects are
unchanged.

These APIs are a native kernel foundation. A later append-only verifier/recipe
operation, versioned worker/MCP transport and compact preparation review must
resolve exact tracked sources and preserve independent intended-target checking.
The ordinary In-Fusion reference route remains valid and unchanged; enabling
this kernel cannot mean a completed proprietary In-Fusion/host-rescue model.

## Bounds and scientific fixtures

Total supplied stored-strand bases are bounded to 2,000,000 before assembly
traversal/cloning; source and evidence serialization are streamed through a
bounded hash writer. The comparison budget applies to selected homology and
alternative terminal comparisons. Bounded source-validation, hashing and literal
copy work is separately linear; it is not mislabeled as terminal comparison
work. Alternatives are only those compatible with the supplied, exposed tails,
not an exhaustive reaction mixture search. Exceeding a work/candidate cap retains
explicit incomplete evidence and never claims uniqueness.

Eight independent fixture groups specify actual expected sequences, end groups,
cut coordinates and strand discontinuities: asymmetric clipping; unchanged
zero-cut groups and blocked/exhausted/nicked inputs; repeated clipping of
staggered ends; a paired unsealed circle; asymmetric gaps and real linear outer
5′ tails; repeated homology/budget limits and forged report replay; one-linear
versus self-circularization; and retained/removed/reversed modification calls with
a tracked-operand bridge. No experimental yield, cellular repair or proprietary
activity is inferred from those fixtures.
