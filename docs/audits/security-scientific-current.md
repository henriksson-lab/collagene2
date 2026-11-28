# Scientific-core security review

Integrated verification: native/WASM compilation and the connected **490-case**
focused batch pass, including the relevant admission/owner/search regressions.
[Exact case list and scope](features-security-current-integration.json).
Application HTTP, browser and live-provider acceptance remain unexecuted.

Source review of the current working tree, 28 November 2026. No compiler,
scientific test, browser, application HTTP request or network attack was run.
This is a bounded review, not a complete security certificate. Native document
and candidate decoding amplification is being handled by another worker.

## Confirmed: exact sequence/PCR search lacked a work bound

`core/src/sequence/analysis.rs::find_matches_view` tested each possible start by
walking the whole compatible prefix. Its 100,000-hit cap bounds retained match
memory, but does not bound unsuccessful comparisons. `Primer::new` accepts
sequence-valid annealing strings up to the model's 10-million-base limit.
`simulate_pcr_geometry_view` called the matcher four times before its later
10,000-pair bound. The public native MCP `simulate_pcr_inner` calls this exact
engine directly, so stricter controls in the separate mismatch-PCR engine do not
protect this path. Basic worker `AnalysisJob::Search` also used this kernel.

A prior-imported circular template `A` repeated 1,000,000 times and forward
primer `A` repeated 99,999 times followed by `C` are valid concrete DNA. Every
start compares all 100,000 primer bases before the final mismatch: 100 billion
comparisons and no hit. A reverse primer of `T` repeated 100,000 times also
passes the primer constructor. The request's primer data is about 200 kB; it
fits hosted request limits. Scaling template/primer lengths to 10 million /
1 million remains within model limits and gives 10 trillion comparisons.
The zero-hit example never reaches the retained-hit cap. Running native
`spawn_blocking` work cannot be stopped simply by dropping its RPC future.

Impact: authenticated scientific-work exhaustion in the shared hosted process,
or a stalled dedicated worker/native session. This finding does not grant
cross-member data access or fabricate a verified cloning product.

### Source resolution, integrated checks pass

The duplicate exact/typed sequence loops now call `sequence::matching`: linear
KMP for concrete nucleotide or literal protein text, and a comparison-budgeted
compatible-IUPAC fallback. Overlapping/circular/oriented hits and existing subset
certainty semantics are retained. Exact/PCR APIs return an explicit error if work
cannot be completed; typed search retains its explicit incomplete result. All
four exact-PCR binding searches share one 200-million-work job budget, including
normalization and reverse-complement costs. Prefix entries are u32 (at most
40 MB for the maximum model pattern), with work admission before allocation.

Document search now has its own `sequence-search-2` engine identity and rejects
superseded report engines, because its work accounting now includes preprocessing
and KMP rather than only naïve comparisons. No wire shape or verifier chemistry
changes. Six new cases cover independent overlap/origin/RNA/protein expectations,
IUPAC uncertainty/exhaustion, linear-work late mismatch, current engine admission
and the actual native MCP exact-PCR path. No new test was executed by this worker;
execution belongs to the parent's connected integration batch.

## Reviewed guards and negative evidence

- Pairwise alignment matrix layout uses checked addition/multiplication plus
  64-million-cell / 96-MiB bounds before allocation. A separate length × maximum
  score bound of 500 million prevents the apparent affine i32 arithmetic issue,
  including long narrow matrices. No overflow finding is claimed there.
- Contig overlap matrices receive a checked `(n+1)*(m+1)` count and pair/cumulative
  cell admission before the private allocation. The visible private unchecked
  multiplication is dominated by that guard; it is not a confirmed wasm32 bug.
- Dot-plot grid counts are converted to u64 before multiplication; the apparent
  10-million × 10-million wasm32 product is not performed in usize.
- `Sequence` construction/deserialization rejects empty input and keeps its
  fields private. A proposed modulo-zero restriction/digest input therefore
  cannot be constructed through the supported public/native model path.
- Ligation profile shifts are dominated by k=3/4 admission; the Matrix's k field
  is private. No arbitrary-shift input was demonstrated.

These observations apply to inspected callers and guards, not every exported
algorithm. Scientific approximation limits and missing calibration profiles
remain product scope items, rather than security defects by themselves.
