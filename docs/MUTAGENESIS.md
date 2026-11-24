# Primer-encoded mutations and inverse PCR

`collagene_core::mutagenesis` is independently implemented Rust. It invokes no external scientific executable, Primer3, alignment engine or proprietary design service. The module separates declared primer-binding compatibility, final oligo-encoded sequence, and the physical steps required to close an inverse-PCR product.

## Conditional mismatch PCR

`simulate_pcr(source, forward, reverse, options)` takes two actual 5′→3′ oligos using the existing `analysis::Primer` representation: an annealing domain, a separate 5′ tail and a name. The domains may intentionally differ from the reference. Both oligos are searched on both strands and are evaluated in both physical roles, including same-oligo pairs. Names do not restrict an oligo to its intended forward/reverse role. Identical aliases with the same domain/tail are evaluated once; different binding origins or domain/tail decompositions remain separate geometries even if their final sequence bytes coincide.

`PcrOptions.binding` records allowed disjoint mismatches, a required consecutive certain 3′ run, binding-hit and work limits. Mismatch and ambiguity positions are in each oligo's own 5′→3′ coordinates, including reverse hits. Defaults permit one disjoint mismatch and require one certain terminal 3′ base. These defaults are an explicit enumeration policy, not a recommended experimental condition. No mismatch thermodynamic/kinetic effect or empirical probability of extension is inferred. Increasing the allowed mismatch count is not evidence that those sites amplify.

For an inward-facing non-overlapping pair, the final top reference is reconstructed as:

1. The actual forward oligo's 5′ tail.
2. The retained template interval, with its prefix replaced by the actual forward annealing-domain bases and its suffix by the reverse-complemented reverse domain.
3. The reverse complement of the reverse oligo's 5′ tail.

The product retains source positions, domain/tail lengths, both binding hits and the introduced source/product base changes. A deliberate concrete substitution does not become an uncertain product base merely because its initial binding contains a mismatch. Ambiguous source evidence remains distinct: compatible IUPAC template binding or unresolved copied bases prevents concrete complete enumeration, even when an oligo provides a concrete terminal base. Degenerate ordered-oligo mixtures are explicitly unsupported.

The model covers at most one circular template turn and bounded inward intervals on a linear reference. Binding sites crossing the circular origin retain exact locations. If the annealing domains overlap, compatible encoded bases can produce a **sequence preview**, but `unsupported_geometries` increases and `enumeration_complete` remains false. Conflicting overlap bases do not acquire an invented repaired duplex; `conflicting_overlaps` records the unresolved geometry. This is not a model of QuikChange-style nicked whole-plasmid synthesis, strand displacement, heteroduplex repair or polymerase-specific terminal additions.

`search_complete` means the configured placement/pair/product budgets were not exhausted. `enumeration_complete` additionally requires no unsupported overlap/capacity geometry, conflicting overlap or ambiguous evidence. It means complete under these declared rules, not exhaustive empirical PCR outcomes. A truncated or unsupported report cannot establish a unique product. Callers constructing a verified physical intermediate must check these flags rather than selecting a convenient matching sequence from an incomplete report.

Defaults are 10,000 checked pairs, 1,000 products and 32 million materialized product bases. Hard limits are one million pair checks and 10,000 products; the base budget cannot exceed 32 million. Individual oligos are at most 200 ordered concrete DNA bases. Binding work/hit limits retain the basic primer module's bounds. Conflict and capacity checks occur before copying full template intervals, and rejected pair geometries do not accumulate an unbounded deduplication structure.

## Inverse-PCR replacement design

`design_inverse_pcr(source, target, options)` requires a circular DNA source and a **separately supplied intended circular DNA target**. The source interval `[edit_start, edit_end)` is replaced by `inserted`. An empty interval inserts sequence; empty inserted sequence deletes the interval; a nonempty replacement expresses substitutions or combined changes. A circular end may pass the origin, but the deletion must retain some source DNA. Input coordinates are zero-based half-open; the web form uses one-based boundaries with an exclusive end.

The user specifies bounded forward/reverse annealing lengths of 8–80 bases. The domains flank the omitted interval in back-to-back orientation and may not overlap within the retained source turn. Exact biological locations are extracted independently to construct each oligo. The requested replacement is supplied by 5′ tails, not smuggled into a copied reference product.

`reverse_tail_bases` splits the inserted DNA as `A + B`: the forward tail is `B`, and the reverse tail is `reverse_complement(A)`. The expected linear product ordering is therefore `B + retained_source + A`; circular closure makes the inserted sequence contiguous as `A + B`. Tails remain subject to the 200-base total oligo limit. The original source and intended target remain unchanged.

The design returns the actual primers, locations, source/target revisions and a fully specified candidate plan:

1. Exact PCR from the immutable original source, retaining competing binding geometries.
2. Explicit modeled phosphorylation of ideal PCR termini.
3. Sealed blunt circular ligation, with no accepted nicks.

The existing independent verifier reconstructs that plan against the separately bound target. Wrong targets remain invalid; competing exact-PCR products remain incomplete; already reconstructed intermediates remain available. The returned candidate does not claim that a commercial enzyme mixture, template digestion, transformation or cellular repair has been performed. Missing laboratory conditions and experimental uncertainty remain visible.

Back-to-back and overlapping primer architectures are distinct experimental methods, described in [NEB's site-directed mutagenesis guidance](https://www.neb.com/applications/cloning-and-synthetic-biology/site-directed-mutagenesis). Its [primary primer-design FAQ](https://www.neb.com/en-us/faqs/how-do-i-design-primers-to-use-with-the-q5-site-directed-mutagenesis-kit) distinguishes insertion tails and deletion-flanking designs. Those sources inform the method boundaries; no external design service or implementation is used. The phosphorylation/ligation sequence above is Collagene's explicit chemical model, not a claim to reproduce a proprietary reagent formulation.

## Web and verification boundaries

The **Mutagenesis** panel chooses an immutable source document and runs conditional mismatch-PCR previews in a binary WASM worker. Actual primer domains, tails, mismatch and 3′ policy are editable. Product references use the verifier's annotation projection helper; altered or uncertain CDS/feature coordinates require review and stale translation assertions are removed. Reference previews do not confer physical end chemistry or target verification.

For inverse PCR, the active document is the intended target and the source is chosen separately. The panel displays the generated oligos, explicit chemical steps, current deterministic findings and full evidence download. **Review inverse-PCR candidate** transfers a proposal without an authoritative status; the Plan panel binds current source/target revisions and re-verifies it. Source, target, option, selector, request-generation and component-lifetime checks discard stale worker results. Product and plan JSON downloads retain the full returned evidence and limits.

## Validation

`cargo test -p collagene-core --test mutagenesis --locked` checks manually specified forward/reverse substitutions and tails, primer-oriented 3′ mismatch rules, circular origin source/product coordinates, compatible versus conflicting overlaps, ambiguous evidence and all truncation budgets. An independent reference-interval/oligo-constraint oracle checks 3,072 short template/oligo/topology/policy combinations rather than reusing the binding/product enumerator. Zero-mismatch non-overlapping products are also compared with the existing exact-PCR definition.

Inverse-PCR fixtures cover substitution, insertion, deletion, origin-spanning replacement and every split of each inserted tail. The existing verifier independently checks all three physical steps against the external target. Wrong targets, repeated-template alternatives, invalid splits, overlapping domains, unsupported alphabets/physical metadata and oversized inputs fail or remain incomplete explicitly. Reports are bounded to the existing 32 MiB portable serialization limit.
