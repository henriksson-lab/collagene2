# Endpoint-aware affine row-band review

Independent mathematical/source review and current integration contract, 28 November 2026. The current source uses registered `alignment/{mod,kernel,grid}.rs` and the existing review/history/export owners. Earlier unregistered storage work has been integrated into one native recurrence. The registered implementation passes the shared native/WASM checks and independent literal/oracle fixtures; [source and execution evidence](connected-design-current-integration.json) retains corrections and limits. Graphical acceptance is unperformed.

## Mathematical result

The certificate is sound under the admitted global affine model. Paths comprise paired steps and single-axis gap calls. Changing directly between gap directions opens a new gap and does not invalidate the per-call upper bound. For endpoint delta `n-m`, leaving the endpoint-aware corridor requires at least `abs(delta)+2*(r+1)` total gap calls. This holds on either side of the corridor and for unequal endpoints; every complete path satisfies `2P+G=n+m`.

The implementation uses `max(0,match_score,ambiguous_score,mismatch_score)`, rather than assuming ambiguity scores are positive. With strictly negative gap scores, replacing pairs by more gap calls only decreases `P*maxSub+G*maxGapUnit`. Thus the bound covers negative ambiguity, opening-versus-extension preferences, alternating opposite gaps and all actual omitted paths. Using zero when all substitution scores are negative is conservative. If minimum gaps exceed `n+m`, no complete omitted path exists. Otherwise parity is exact and integer division loses no half pair.

Strict selected score greater than the omitted-path bound excludes off-band ties as well as better paths. Equality remains unproven. Neither this certificate nor a full domain proves unique correspondence inside the domain; traceback retains the existing deterministic tie ordering. Circular origin/strand coverage and candidate-level proof are separate facts.

A comparison's observed incumbent may safely discharge another candidate when its score strictly exceeds that candidate's omitted-path upper bound. The candidate's computed in-band maximum is already no better than the incumbent; omitted paths are strictly worse. The candidate's **own** unrestricted optimum remains unproven, reflected in typed evidence. This is not a biological confidence rule.

No mathematical counterexample or concrete coordinate-storage fault was identified in this source review. That is a bounded review result, not execution evidence or a proof that all integration behavior is correct.

## Current arithmetic and resource owners

- `AlignmentOptions::validate_for_lengths` is the canonical admission owner before recurrence and boundary initialization. Absolute path scores are bounded by 500 million, increments by 10,000, and `NEG=-1 billion` remains disjoint. Unreachable-neighbor addition preserves the sentinel; admitted reachable arithmetic stays within `i32`.
- Proof geometry and score bounds use checked `i64`. Failed conversion, radius increment, multiplication or addition refuses; unsigned radii cannot wrap into narrow domains. Exact admitted-cell counting uses checked `u128`, allowing narrow-band admission on wasm32 without overflowing an intermediate full rectangle.
- `RowLayout` retains original-coordinate inclusive ranges, endpoint-aware addressing and checked offsets. Outside neighbors are three unreachable states; strip edges never create free starts or gap closure. Full domains use direct rectangular indexing without row-offset amplification.
- Trace bytes, two maximum-width score rows and offsets are admitted against the existing 96 MiB payload bound before allocation. Each attempt owns and drops its grid; widening retains only the completed fallback path/receipt. Output is bounded to one million columns and is never silently truncated.
- `CellBudget` borrows one shared actual-cell ledger. Evaluated row cells consume both per-alignment and comparison caps monotonically; retries and candidate changes do not refund or reset work. Allocation/admission failure does not claim a completed evaluation. Identity/seed scans are distinct from DP-cell accounting.

## Current recurrence, scheduler and evidence contract

There is one Pair/Insert/Delete recurrence in `kernel.rs`, with the previous opening/extension model, opposite-gap transitions and predecessor ordering. Full and corridor layouts feed that same recurrence. Local and QueryInReference keep full domains and original initialization/endpoints; adaptive-global policy with those modes refuses explicitly. Global transposition preserves source coordinates and Insert/Delete direction.

Standalone adaptive alignment widens deterministically until strict proof, full domain or a typed radius/attempt/cell/memory stop. The last completed path survives a refused wider admission with unrestricted optimum explicitly unproven. Whole comparison performs breadth-first initial-band passes over bounded requested origins/strands before widening unresolved candidates, so a wrong early orientation does not monopolize the ledger. The actual winning path and bounded scalar receipts are retained; grids are not cloned across candidates.

One `Alignment` includes actual full input lengths, options and `SearchEvidence`. Completed attempts retain domain, actual cells, score, omitted-path bound and certificate. `SearchEvidence::validate` checks exact domain cell counts, progression, cumulative limits, final score and stop reason; `Attempt::completed` alone is not treated as full documentary validation.

`Alignment::validate_search` additionally validates canonical scoring admission, current algorithm, full source lengths, mode boundaries, CIGAR/score/differences/statistics and exact-symbolic proof premises against actual complete equal nongap rows. `MoleculeComparison::validate_search` binds the selected path and cumulative candidate evidence; review validation checks candidate rosters against actual requested strand/origin policy and completeness claims. Removing an unselected orientation and adjusting totals cannot preserve a false exhaustive declaration.

Raw worker alignment requests bind both molecules and options using the existing review revision owner. Fresh worker/MCP outputs are checked against known input lengths/options; bounded response decoding validates the current native proof/shape. Historical loading/export remains documentary and never realigns or creates material authority. Current internal option/report/history forms are updated directly without old-format fallback decoding.

## Independently authored fixtures and pending acceptance

The tiny oracle independently enumerates legal diagonal/horizontal/vertical paths and sums substitution or gap-run opening/extension costs; it does not share production DP storage/recurrence. It permits opposite-gap switches. Escape assessment records actual maximum omitted-path scores rather than reproducing the bound. Both length signs, empty endpoints, short A/C strings and all tiny radii exercise endpoint geometry; ambiguous inputs and extreme/nondefault gap preferences add other score classes. A literal off-band tie exercises strict inequality.

Authored cases also compare full/adaptive scores and deterministic CIGAR, check row-index bijection/unreachable neighbors, actual cumulative retry/candidate charges, overflow/resource admission and forged receipt fields. Integration fixtures cover literal 10/20 kb variants, reverse/rotated candidate scheduling under a shared cap, preserving an unproven path after refused widening, worker source/options/length/proof binding and actual MCP summary/artifact scope. Existing local/query-in-reference, history, source-mapping and export cases remain registered. A specific roster-tampering fixture removes genuine candidate evidence and repairs totals while retaining false completeness claims; the current validator refuses it.

These fixtures are **authored, not claimed executed here**. Remaining shared-batch acceptance includes compiler/native/WASM checks and fixture execution, especially full/adaptive tie parity, source/context/archive tampering, source-mode boundary behavior and finite resource refusal. No timing or empirical accuracy result is asserted.

The role-specific supplied-consensus/design owner and path-conditional CDS consequence workflow remain subsequent work, documented separately in [design-consensus-current.md](design-consensus-current.md). This affine phase supplies a bounded computational alignment contract, not read confidence, depth, sample identity, biological function or experimental success.
