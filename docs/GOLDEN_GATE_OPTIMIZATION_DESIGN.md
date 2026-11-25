# Next cloning design: fusion placement and domestication

Parallel source and official-documentation reviews, 25 November 2026. These are
design decisions for unfinished work, not implementation or application acceptance.
Retain the existing source-bound primers, selective junction proposals, reaction
discovery, explicit candidate keys, fresh verification and recovery machinery.

Implementation follow-up: [signed fusion placement](GOLDEN_GATE_FUSION_DESIGN.md)
now composes actual original-template preparations, and
[source-bound domestication](DOMESTICATION.md) adds positional synonymous search,
full-source/CDS checking and explicit inverse-PCR preparation. The missing-adapter
discussion below records the original design assessment; broader optimization,
joint coding constraints, curated data and physical routes remain in scope.

## Actual missing adapters

The existing `golden_gate_primer_design` and `golden_gate_junction_edit` modules
already generate actual PCR/digest preparations, use profile-sized fusion ends,
reject extra sites and retain unmodified manual rows. `assembly_junction_edit`
already builds independent trim/insert intent. Add source-bound search over
alternative junction placements that preserves the intended assembled DNA, rather
than another primer or cleavage simulator.

For adjacent original contributions A and B, a rightward boundary move can transfer
a literal prefix of B into A's extension while trimming that prefix from B's
contribution. Freeze expected target DNA before searching; independently check
every proposed representation against it. General placement windows, including
leftward movement and circular-origin policy, remain part of complete optimization
scope. A first implementation profile must not redefine that completion boundary.

Search mathematical fusion constraints cheaply, preserve fixed/manual ends, expose
bounded work and incomplete results, then materialize a chosen proposal through
the existing actual primer/digest compiler. Palindrome, distinct-end and Hamming
rules are design constraints, not empirical fidelity estimates. Explicit Apply
must invalidate prior pool choices and native evidence. Optimization must never
select a reaction winner using the target.

Generic synonymous motif removal already exists in `codon_design`, with forbidden
motifs, contexts and an independent checker. Domestication needs a source/CDS
adapter: preserve actual phase, strand, genetic code, exceptions and affected
feature constraints; show every proposed base change; independently check protein
and motifs before accepted intent. Do not erase unsupported qualifiers to make a
feature appear silent. Reference editing, annotation-only changes and supported
mutagenic preparation are separate outcomes. An exact PCR cannot silently enact
arbitrary internal synonymous substitutions.

## Documented interaction references

SnapGene documents junction-local overhang adjustment with reading-frame context
and affected-primer redesign. Reuse that interaction in our existing junction
review, preserving manual settings until explicit Apply.
[Official fusion-adjustment guide](https://support.snapgene.com/hc/en-us/articles/10384177073300-Golden-Gate-Cloning-of-Multiple-Fragments-into-a-Vector-with-no-Type-IIS-sites).

Its published scoring convention uses condition-specific empirical ligation data
and multiplies unique-end fidelities. Dataset conditions and missing-data cases
must be visible; a ranking is not a universal efficiency or experimental outcome.
[Official fidelity explanation](https://support.snapgene.com/hc/en-us/articles/10384163522196-How-is-Golden-Gate-Fidelity-Predicted-in-SnapGene).

The silent-removal guide distinguishes direct sequence mutation, a mutagenic PCR
product and annotation only. It previews synonymous codon changes and unchanged
translation. Those distinctions fit our existing independent mutation intent and
physical reconstruction boundaries.
[Official silent-site removal](https://support.snapgene.com/hc/en-us/articles/10383898655124-Remove-a-Restriction-Site-From-a-Coding-Sequence-Without-Changing-the-Translation).

Current Benchling public documentation does not establish an optimizer formula or
complete domestication procedure. Do not attribute our algorithm to it.
[Current assembly guide](https://help.benchling.com/hc/en-us/articles/39922419457165-How-to-use-the-combinatorial-cloning-and-concatenation-assembly-tools).

## Data and licensing boundary

Keep independently written software MIT licensed. Published ligation matrices are
a possible separately attributed data package; verify actual supplement notices,
conditions, original hashes and normalization before including them. The Pryor
2020 publisher article lists relevant supplementary matrices and a CC-BY license.
No dataset has been bundled by this review.
[Publisher source](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0238592).

Publication access does not authorize borrowing optimizer implementation code.
The protocol paper describes noncommercial code terms, and NEB's tool page has its
own restrictions. Reimplement published algorithms independently; preserve any
data attribution and avoid relabelling third-party data as MIT.
[Protocol statement](https://currentprotocols.onlinelibrary.wiley.com/doi/10.1002/cpz1.882),
[NEB tool terms](https://www.neb.com/en-us/research/nebeta-tools).

### Supplement inspection follow-up

The five original Pryor2020 supplementary XLSX files were inspected separately,
without adding a dataset to the repository. They contain complete integer
ligation-frequency matrices with explicit zeros. Row labels use reverse-complement
ordering relative to lexical columns; import must read the actual labels rather
than treat the visual diagonal as same-label self-ligation. Matrices are symmetric
after mapping by their labels. Raw directional-read preprocessing is not specified
in the files, so preserve published-frequency semantics instead of claiming
independent directional observations.

Their cycling conditions differ from the older static profile described above,
and formulation identity must remain explicit. A bounded independently written
ZIP/XML importer can use the existing Rust dependencies. Publisher-level CC-BY
evidence supports a separately attributed data package; it does not relabel the
data MIT or authorize copying an optimizer implementation.
[Paper, methods and supplementary inventory](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0238592),
[publisher licensing policy](https://journals.plos.org/plosone/s/licenses-and-copyright).
