# Inspect and select a freshly prepared Golden Gate pool

`inspect_prepared_golden_gate_pool` takes immutable original `inputs`, an independent
`target` DocumentRef and `request: { recipe, options }`. The recipe is the existing
typed ligation preparation graph; options are the explicit finite Golden Gate
pool/orientation/work policy. The tool reconstructs actual preparations and their
physical duplexes before enumerating. Missing any required physical preparation
withholds the entire pool rather than silently enumerating a successful subset.

The response remains `state: "advisory"`, `verified: false`. It reports exact
context/plan/target revisions, preparation completeness, missing row indices,
enumeration completeness, work/truncation notices, roles and every retained
candidate summary. Candidate summaries retain directed/physical identity keys,
source-row indices, target-match result, length and recognition-site/recut evidence.
The full JSON artifact is `{ kind: "prepared-golden-gate-pool-v1", inputs, target,
report }` and includes the original request, preparation preview, raw pool and
candidate assessments. Target matching is separate from whether background
products exist. `exhaustive_mixture_search` remains false even when the stated
finite one-copy pool enumeration is complete.

`select_prepared_golden_gate_candidate` accepts that artifact ID and
`selection: { context_revision, candidate_index, directed_sequence_key,
physical_sequence_key }`. It returns an unverified saved candidate, parent/new
context revisions, plan revision and a binary native candidate artifact. Both
`artifact_id` and `candidate_artifact_id` name the same base64-readable
`candidate.cg2plan` artifact.

Selection preserves every original input and PCR/catalog-digest/repair/treatment
operation. It changes the final ordered ligation to the existing Golden Gate
operation using the actual prepared sources, explicit pool policy and selected
candidate. Selecting is not verification or restoration of authority from stored
JSON: call `verify_plan` on the returned plan revision so the deterministic
verifier reconstructs the preparations and selected physical product against the
same independently retained target. A selected self-ligation background can be a
real reconstructed product while still being an invalid match to that target.
No target or original source is rewritten to make it pass.

The SDK fixture in `crates/mcp/tests/prepared_golden_gate_tools.rs` imports a literal
circular 40-base target before primer design. Two 20-base originals share an
`ATGA` end convention, producing self-A, self-B and two-fragment physical groups.
It checks distinct background/target summaries, complete evidence bindings,
selected-plan prefix/input preservation, fresh target and self-product verification,
immutable deduplication and original exports. A one-state cap remains explicitly
incomplete; failed preparation withholds the pool. Stale source/context and
incorrect candidate identity keys are rejected. Validation results are recorded in [implementation status](IMPLEMENTATION_STATUS.md).
