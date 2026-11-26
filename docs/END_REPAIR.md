# End repair and oligo annealing

This independently implemented model transforms exact duplex geometry and enumerates complementary oligo placements. It is embedded Rust and runs without an installed scientific executable. It is not an enzyme-kinetics, melting-temperature, concentration, yield or laboratory-completion prediction.

## Physical representation and provenance

Every strand is stored in its own 5′→3′ order. The top strand starts at reference coordinate zero; `bottom_start` is the left boundary of the antiparallel bottom strand in those reference coordinates. Only overlapping, concretely complementary A/C/G/T strands become `Duplex` states. Ambiguous/RNA input is rejected by this physical model rather than silently instantiated as one concrete molecule.

`StrandMapping` addresses the original and product strands in their own stored 5′→3′ coordinates. Identity blocks copy retained bases. A cross-strand synthesized block has `reverse=true`: the new strand bases are the reverse complement of the indicated existing template interval. Removed bases have no product mapping. Reports preserve both original and resulting actual strand sequences, added/removed 3′ suffixes, chemical states and nick boundaries. These are inspectable reconstruction data; imported metadata does not establish trusted physical chemistry.

## End-repair methods

- `FillFivePrimeOverhangs` extends the recessed opposite 3′ terminus across its existing 5′-overhang template. It does not fill a protruding 3′ overhang.
- `RemoveThreePrimeOverhangs` removes a protruding 3′ suffix up to the initial paired boundary. It does not remove 5′ protrusions or resect into the paired duplex.
- `Blunt` performs the applicable combination at each end.

These ideal endpoint transforms follow the synthesis/exonuclease directionality described by the [NEB T4 DNA Polymerase mechanism](https://www.neb.com/en-us/products/m0203-t4-dna-polymerase) and [NEB blunting guidance](https://www.neb.com/en/faqs/which-dna-polymerase-can-i-use-to-blunt-dna). The [polymerase directionality FAQ](https://www.neb.com/en-us/faqs/can-i-use-a-dna-polymerase-to-fill-in-3-overhangs-or-remove-5-overhangs) distinguishes unsupported 3′ fill-in and 5′ removal. These sources support the modeled directionality; selecting a software method does not assert complete conversion with a particular enzyme or reaction condition. The [manufacturer's blunting protocol](https://www.neb.com/en-us/protocols/protocol-for-blunting-ends-by-3-overhang-removal-and-fill-in-of-3-recessed-5-overhang-ends-using1) also describes undesired resection under other conditions, which is outside this ideal endpoint model.

An affected 3′ end must already have an explicitly known accessible hydroxyl. Blocked/unknown ends are unsupported; no blocking-group removal is inferred. Unaffected ends retain their exact state. No mode changes a 5′ phosphate/hydroxyl/unknown state, adds phosphorylation, or seals an internal nick. A no-op on an already blunt or nonmatching-polarity end preserves chemistry rather than asserting successful treatment. Circular inputs are rejected because they have no free ends.

Nicks strictly within the initial paired interval are retained. A nick at or outside that interval can separate an unpaired/template segment, so this model refuses to establish the repaired product. More detailed handling of disconnected overhang pieces, strand displacement, gaps, partial/excessive exonucleolysis and unknown nick-end chemistry remains outside this operation.

## Oligo annealing

`anneal_oligos(top, bottom, options)` takes two declared 5′→3′ DNA strands. It checks all antiparallel offsets with at least `min_paired_bases` overlapping bases. Every base in that overlap must be concretely complementary; mismatches and IUPAC compatibility do not become an exact physical duplex. Nonoverlapping original strand ends remain actual overhangs. No bases are synthesized or removed, and the model does not choose the longest overlap.

The top and bottom input `Terminus` values are explicit source/laboratory assumptions. Defaults are unphosphorylated 5′ hydroxyls and 3′ hydroxyls, matching ordinary unmodified oligo assumptions; callers must explicitly specify modified/unknown states. Annealing itself does not phosphorylate or ligate. The [IDT annealing guidance](https://www.idtdna.com/page/support-and-education/decoded-plus/annealing-oligonucleotides/) separately discusses making duplexes from complementary single strands and supplying 5′ phosphates when required for ligation. The [NEB cloning guide](https://www.neb.com/en/tools-and-resources/usage-guidelines/cloning-guide) describes using complementary oligos with overhangs as cloning inserts.

Candidates retain source sequence hashes, exact offset, actual strands, overlap coordinates, ancestry and stable IDs. A report sets `unique=true` only after completing the bounded supplied-pair offset search and finding exactly one candidate. Multiple candidates remain separate alternatives, and capped/unfinished searches never establish uniqueness. Search completeness refers only to this supplied pair and this model: hairpins, self-association, mismatches, secondary structure and thermodynamic/kinetic competition are not searched.

Defaults are 12 paired bases, 20 million compared bases and 100 retained candidates. Absolute bounds are 10,000 bases per oligo, 200 million compared bases, 1,000 retained candidates and 32 MiB serialized reports. Repair preflights original-plus-product strand bytes before constructing output copies and checks the same report cap. Work/serialized bounds are not an assertion that process RSS equals the serialized size.

## Independent validation

`cargo test -p collagene-core --test end_repair --locked` uses manually specified 5′/3′/mixed overhangs and resulting actual strand strings, both chemical termini, reverse-orientation equivalence, retained/rejected nick boundaries, explicit phosphorylation before ligation and complete source-to-product ancestry coverage. An independent absolute-reference interval matrix checks every overlapping strand interval pair on multiple concrete reference strings. Annealing is compared against 12,288 independent three-base pair/minimum-overlap cases using manual allele complements and coordinate overlays, without another cloning program as an oracle. Separate tests cover repeated offsets, partial/candidate limits, no matches, unknown/blocked states, RNA/ambiguity, JSON reconstruction and output bounds.

The web End repair / annealing panel computes restriction-derived repair and annealing previews in workers. Its active document remains an independently selected target; reviewing a candidate does not verify it. Optional phosphorylation and no-nick circular ligation are separate typed steps. Opening a preview creates a reference document whose metadata cannot authorize later physical operations. HTML/PDF step sheets include exact ordered oligos and declared terminal chemistry; native workspace plan schema 2 and candidate v4 distinguish the new operations.

## Annealed insert plus a digested backbone

Use **Design assembly → Digest and ligate** for a connected source-to-product
workflow. Set the backbone slot from its actual pinned virtual-digest fragment.
For the insert slot, open the ordered top oligo and choose **Use whole source**
(or explicitly declare its intended source interval). Open the ordered bottom
oligo, expand **Use an annealed oligo insert**, and capture the viewer source as
the bottom strand. Both original sequences are written5′→3′.

Review the intended construct from the independently selected backbone and top
contributions. The bottom source remains bound to the declaration but is not
concatenated into that reference as a third fragment. In the active preparation
card, set pairing/search limits and declared termini, inspect the annealing
geometries, and explicitly choose one. Phosphorylation after annealing is a
separate terminal-treatment choice. Compute the full preparation/ligation against
the frozen target before reviewing the product and creating named outputs.

The plan retains the backbone and both oligo originals; it reconstructs catalog
digestion, annealing, treatment and ligation in one verifier pass. Step sheets
include both complete ordered oligos. Opened reference previews cannot replace
these original operations. Saving/restoring a paired draft checks all three
retained source events. Paired drafts use version4; combinatorial paired-source
projection is not yet implemented. Existing plan14/candidate16 formats already
represent the physical operations; worker71 adds the shared recipe preparation.

### Editing independent ligation intent

In **Declare intended junction DNA**, select an internal junction, circle closure
or final linear suffix. Set the intended upstream/downstream reference trims and
literal DNA. The final linear suffix has no downstream trim. These declarations
retain the captured sources, primers, digest choices, end treatments and selected
annealing geometry. **Review intended construct** checks the reference regions;
subsequent recipe reconstruction verifies whether the prepared strands actually
produce that DNA. An incorrect literal target cannot finalize as verified output.

### Reversing a paired insert

Use **Swap ordered oligo roles** in the active source card. The previous bottom
becomes the new top and its whole ordered sequence becomes the independent
reference contribution. Complete top/bottom terminal chemistry follows the
captured originals. The action resets only the two neighboring junction
declarations and the old geometry choice; unrelated edits and preparations are
retained. Review the new bounds, select a fresh annealing geometry, and reconstruct
ligation again. Both exact original oligo events remain in the plan.

This explicit role swap handles unequal overhangs correctly: the physical new top
is the original bottom strand. Ordinary annealing starts with continuous declared
single strands, so arbitrary pre-existing nicked duplex material still requires
its own supported upstream reconstruction.
