# Codon design and reverse translation

The independently implemented `codon_design` engine (`codon-design-1`) designs synonymous, concrete linear DNA against a separately bound source. It supports both a canonical protein document and a complete, concrete linear DNA coding sequence in frame zero. DNA sources are first translated to derive the required protein; a candidate cannot supply or change that requirement. At least one actual amino acid is required. A single terminal `*` is preserved as a stop codon; internal/multiple stops, ambiguous residues, U/O recoding, ambiguous DNA, partial codons and RNA sources are unsupported.

Genetic codes 1, 2 and 11 use the existing audited translator and the [NCBI genetic-code tables](https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi). Ordinary codon decoding applies at every position: this design profile does not recode alternative initiators as methionine. It does not infer promoter, RBS, signal peptide, cleavage, post-translational processing, expression, host fitness, RNA folding or synthesis success. Retained translation exceptions, incompatible `/transl_table` or `/codon_start`, and known modified/cross-linked source evidence reject rather than silently changing the biological requirement.

## Explicit preferences and constraints

`DesignOptions` records:

- The selected genetic code.
- Per-codon nonnegative integer weights, 0–1,000,000. Unspecified codons have weight 1; weight 0 disallows a codon. The objective maximizes the sum of declared weights, including any terminal stop codon. These are user preferences, not a built-in host-frequency or expression model. Equal scores use lexicographically smallest coding DNA (`A<C<G<T`). Duplicate codon declarations reject.
- Forbidden uppercase IUPAC DNA motifs, checked on both strands. A concrete DNA window matches an IUPAC motif when each base is admitted by its symbol; no ambiguity expansion is needed. Duplicate/reverse-complement-equivalent patterns are checked once.
- Fixed left/right concrete DNA context, which participates in motif, homopolymer and window constraints, including context-only matches and both junctions. Context is not appended to the returned coding DNA.
- Global coding GC bounds in integer basis points. The exact inequalities are `minimum_bp * length <= 10000 * GC_count <= maximum_bp * length`. Context does not contribute to global GC; any terminal stop codon does.
- Maximum homopolymer length, across the complete left-context/coding-DNA/right-context construct.
- Optional full-window GC bounds across that same construct. Only complete windows are constrained; a construct shorter than the window has no such windows. No partial window or circular-origin rule is invented.
- Beam width and a deterministic work budget. Every declaration, including resource options and codon preference order, is bound by `options_revision`.

Motifs, codons and contexts must already be normalized uppercase strings without whitespace; invalid or ignored JSON fields reject. Empty context is allowed. An empty motif is an error. The UI can offer normalization as an explicit preprocessing step, but the core never reports a hash for different effective options.

## Bounded search and independent checking

The search enumerates synonymous codon extensions in deterministic order, rejects newly completed forbidden motifs, homopolymers and windows, and uses safe integer future-GC enclosures to prune impossible prefixes. It orders surviving prefixes by score, then exact lexical rank. It retains only the configured beam and a bounded constraint suffix. Compact shared prefix-arena nodes reconstruct the chosen full sequence; the search does not copy and sort complete growing DNA prefixes at every extension.

The independent checker translates the complete proposed DNA against the original source and scans all declared constraints from scratch. It does not call the incremental search predicate. `design` releases a candidate only when this checker completes and passes. A user or MCP client may also submit arbitrary DNA to `check` against its independent source and options.

`CandidateFound` means a candidate passes the declared sequence constraints. It does not mean that all designs or biochemical outcomes were explored. `search_complete` and `optimality_proven` become false after any beam or memory pruning. If such a search later finds no candidate, the result is `SearchIncomplete`, never `Unsatisfiable`. `Unsatisfiable` is reserved for a complete finite search with no passing design. If the checker exhausts its own budget, `satisfies` is unavailable (`None`), rather than false or true. At most 256 violation details are retained; truncation is disclosed. Lack of a candidate after incomplete search does not establish infeasibility.

Work units account for introduced bases, actual motif/run/window comparisons, suffix copies in 64-byte blocks, conservative fixed-key ranking costs, retained prefix records and final DNA reconstruction. They are an algorithmic budget, not a promise about CPU cycles or elapsed time. The search and independent checker each receive the configured budget; report their usage separately. Prefix-node allocation and fixed-size sorting remain bounded even when a work limit is reached.

Current hard bounds are 10,000 residues (including a possible terminal stop), 30,000 DNA bases, 4 MiB bounded source evidence, 4,096 features/metadata entries, location depth 32 and 100,000 location nodes with separate evidence-item/byte accounting, at most 64 distinct codon weights, 64 motifs of at most 128 bases, 2,048 bases per context/window, beam 1–4,096, and work 1–100,000,000. A million compact prefix nodes and a 16 MiB frontier budget bound retained search data. Default work is 20,000,000 and beam is 128; kiloresidue genes are exercised with independent native assertions. A local unoptimized native measurement produced checked 3,000 bp candidates from 1,000 leucines in 1.041 seconds (17,875,661 work units) and from 50 copies of the 20 canonical residues in 0.460 seconds (9,471,111 units). Both searches pruned alternatives and made no global-optimality claim; these are local development measurements, not cross-platform or WASM performance guarantees. Large motifs, windows, context or declared homopolymer lookback consume more work. Reaching resource limits is explicit and never proves a negative result or optimum.

## Evidence and transport

Each report binds source revision, options revision, engine and required protein. Checker reports also bind the exact candidate revision. Returned DNA retains the required protein, source/options hashes, complete options JSON and engine as provenance metadata. Original source annotations/evidence remain on the source; peptide or original-CDS coordinates are not copied to changed DNA. These metadata are historical declarations, never physical chemistry or imported verification authority.

`validate` checks bounded report structure. `validate_bindings_for_source` and `validate_bindings_for` additionally check exact source/options/candidate bindings without a full constraint scan, for responsive UI result handling. `validate_for_source` and `validate_for` recompute the independent checker and are useful for native transport/audit boundaries. None of these binding checks authenticates a remote producer's search-exhaustiveness or global-optimality claim; those claims require a fresh trusted search. The browser runs design and an explicit fresh checker in workers, and invalidates results when source/options change.

The original analysis has its own engine version. It adds no physical `Operation`, changes no cloning verifier semantics, and does not change candidate-plan or archive-plan schemas. Advanced host-frequency datasets, CAI-style objectives, codon-pair/RNA models, repeat/structure/supplier-specific constraints and embedded synthesis-order integrations remain further F26 work.

## Positional constraints and minimum changes

The separate versioned `ConstrainedPolicy` API reuses the same bounded search and
independent checker while preserving the original serialized options and reports.
It accepts per-position synonymous codon domains, protected coding-base ranges
and scoped forbidden motifs on both strands. Scoped recognition includes fixed
flanks and either coding boundary; it does not turn the whole motif list into a
global deletion rule. This permits preserving a wanted occurrence elsewhere in
the same CDS.

Its objective minimizes changed original nucleotides, then maximizes declared
codon preference, then chooses a deterministic lexical tie. Locks and domains
filter extensions inside the search. Independent checking recalculates positional
violations, edit distance, translation and generic constraints from full DNA.
Pruning and exhausted work still withhold global-minimum or impossibility claims.

`ConstrainedDesignReport` contains its own `ConstrainedOutcome`, rather than a
legacy `DesignReport` incorrectly claiming unconstrained weight optimality.
`total_work_used` includes both search and checker work, even when a failed checker
withholds its candidate and receipt. Composite callers can reserve a per-pass
ceiling, charge actual completed work and refund unused allowance; a successful
small first gene must not consume an entire multi-gene design budget.

This API retains ordinary internal-codon semantics. The
[source/CDS domestication adapter](DOMESTICATION.md) separately maps real feature
phase, strand, genetic code and initiation/stop policy, independently checks every
affected feature and builds the explicit physical preparation. Positional search
alone does not establish a synonymous change in a differently framed overlapping
gene, material stock or experimental success.
