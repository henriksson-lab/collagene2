# Pairwise nucleotide comparison

`collagene_core::alignment` is independently implemented MIT Rust. It does not call any executable, use another cloning application as an oracle, or include the deferred MAFFT/MUSCLE/Clustal engines. The scientific algorithm is a bounded three-state affine-gap dynamic program, with local restart and explicit global/end-free boundaries. The mathematical basis is the [Gotoh affine-gap formulation](https://pubmed.ncbi.nlm.nih.gov/7166760/) and [Smith–Waterman local alignment](https://pubmed.ncbi.nlm.nih.gov/7265238/); the implementation and synthetic tests were written for this repository rather than copied from an application.

## Portable API

```rust,no_run
use collagene_model::Sequence;
use collagene_core::alignment::*;
# fn example() -> Result<(), collagene_core::Error> {
let reference = Sequence::dna("TTTACGTGGG")?;
let query = Sequence::dna("ACGT")?;
let options = AlignmentOptions {
    mode: AlignmentMode::QueryInReference,
    ..AlignmentOptions::default()
};
let report = align(&reference, &query, &options)?;
assert_eq!((report.reference_start, report.reference_end), (3, 7));
# Ok(()) }
```

The inputs must use matching DNA or RNA alphabets. Protein alignment, mixed DNA/RNA normalization, quality-weighted alignment and trace consensus are not implemented. Sequence characters remain IUPAC values rather than being resolved to invented concrete bases.

Modes are `Global` (both sequences end-to-end), `Local` (one positive-scoring pair of substrings), and `QueryInReference` (the complete query, with free reference flanks). The latter may legitimately prefer an all-insertion query when no pairing supports a better score; identity is then unavailable and paired coverage is zero. A zero-scoring empty local result is not evidence of matching sequence.

Default scores are certain match +2, compatible ambiguity 0, disjoint mismatch −3, gap open −4 and extension −1. A length-k gap costs `open + (k−1)×extension`. Options record the actual scoring. Strictly negative bounded gap values, positive matches and bounded ambiguity/mismatch values are validated. Transitions between opposing gap runs are permitted with a new opening penalty, so unusually large mismatch penalties do not silently restrict the space of legal alignments.

All output/source ranges are zero-based half-open. `reference_aligned` and `query_aligned` have one character per emitted column, with `-` for a gap. Extended CIGAR uses `=` for a certain concrete match, `M` for compatible IUPAC ambiguity, `X` for disjoint sets, `I` for query insertion/reference gap and `D` for reference deletion/query gap. These `M` semantics are explicitly narrower than a generic unclassified SAM CIGAR `M`; this is not a SAM exporter.

Differences include maximal substitution, insertion, deletion and ambiguous-compatible runs, their actual source symbols, and both coordinate intervals. Zero-width intervals denote boundaries. Identity is calculated over paired columns, with compatible ambiguity reported separately from certain identity. Coverage uses paired columns divided by the **full** reference/query length; clipped flanks and deletion/insertion gaps are not implicitly verified coverage. Span lengths are reported separately.

One deterministic optimum is returned. State tie preferences are fixed, with gap extension preferred when tied, then fixed state order; global/local DP transposes row storage when appropriate. Endpoint ties keep the first visited optimum. Alternative equally scoring alignments are not enumerated. Scores are not statistical significance or experimental confirmation.

## Bounds

Default matrix limit is 16 million cells, configurable up to a 64 million cell hard maximum. Cells include the zero row/column. Three score states use rolling rows; traceback packs all three predecessors into one byte/cell. Matrix/trace allocation is checked against 96 MiB before allocating. A separate numeric bound prevents score overflow, and materialized alignments are limited to one million columns. Identical global inputs with nonnegative ambiguity scoring can use a proven linear shortcut, so no quadratic matrix is charged/allocated for that case; ambiguity still appears as uncertain columns. Errors explicitly report unsupported resource limits instead of switching to an undeclared heuristic. These limits are implementation bounds, not established browser performance measurements.

## Whole-molecule orientations and circular origins

`compare_molecules` accepts `ComparisonOptions`, including global `AlignmentOptions`, query reverse-complement permission, reference-rotation permission, candidate-origin limit and a total cell budget. It returns the winning `Alignment`, input molecule revisions, reference origin, query orientation, evaluated-candidate counts and separate exhaustive/optimality flags. Topology differences are recorded rather than accepted as physical equivalence.

Alignment reference positions address a reference rotated left by `reference_origin`; source boundaries map as `(origin + x) % reference_length`. A reversed query interval `[start,end)` maps to source `[length−end,length−start)`. Report consumers must apply those mappings rather than displaying oriented coordinates as if they were original source coordinates.

An identical symbolic rotation/orientation can use a linear scan when ambiguity scores are nonnegative: it attains the scoring upper bound, proving the best score without enumerating other optima. Otherwise small references can enumerate all permitted origins. Larger references use input origin plus bounded candidates from at most five spaced exact, unambiguous seeds (up to 4,096 hits per seed). Total work is bounded by a separate maximum of 128 million DP cells. Seed-limited or budget-limited searches explicitly have `optimality_proven=false`, with corresponding exhaustive flags/warnings. A better untested origin may exist. No whole-genome search completeness is claimed.

## Validation

Tests use hand-specified gap, clipping, RNA and ambiguity examples, plus exhaustive enumeration of all alignment paths for tiny A/C strings. Local results are checked against all substring-pair paths, and query-in-reference results against all reference-substring paths. Emitted global tracebacks are independently rescored. Circular/RC mapping, complete small-origin search, heuristic limits and resource rejections have explicit regression cases. No external aligner is invoked.

The comparison core is portable to native/WASM. A separate web `ComparePanel` uses the worker API and checks current source/query/options before accepting or exporting a result. Its browser rendering must be validated by the integration test rather than inferred from compilation.
