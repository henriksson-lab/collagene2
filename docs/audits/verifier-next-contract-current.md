# F39: next verifier contract and source batch

2026-11-28, static readiness review only. No source edits, builds or tests.

## Current owners and concrete gap

- `planning/plan.rs` owns the sole `DesignPlan`, `Operation`, `Source` and typed
  `Operation::sources()` visitor. `planning/editing.rs` owns dependency editing,
  count admission and strict backward step references. Keep these owners; add no
  prose-to-plan or alternative execution representation.
- `verification/engine.rs` owns `Execution`, fresh molecular reconstruction,
  physical/tracked variants and final Golden Gate capture. `reconstruct` currently
  returns at the first failed guard or `?`; `evaluate_captured` turns that error
  into one string and stops. Earlier projection/incomplete notes are strings too.
- `verification/evidence.rs` owns `VerificationReport`: three overall statuses,
  `findings: Vec<String>`, dense `products: Vec<Molecule>` and indexed scientific
  evidence. `sources.rs`, `terminal.rs`, `pcr.rs` and `projection.rs` own source,
  physical and annotation facts. Emit named findings at these actual rule sites.
- Every current successful product is a contiguous **leading prefix** of plan
  steps. `products[i]`, private `physical[i]` and `reference_only[i]` all mean
  original step `i`. Current recipe diagnostics infer the failed step from
  `products.len()`. Continuing independent branches while compacting products
  would silently bind a different product to `Source::Step(i)`.

## One current report contract

Freeze exact DTO names before the batch; the following is the recommended shape:

- Replace string findings with `CheckFinding {rule: RuleId, verdict: CheckVerdict,
  scope: CheckScope, subjects: Vec<Source>, coordinates, citations, message}`.
  Scope identifies global/input/step/target, including the original step index.
  Verdict distinguishes pass, fail, indeterminate, informational and blocked.
  Each named rule has a versioned invariant definition and citation locator;
  physical literature/profile citations supplement an internal structural rule.
  Input binding and byte limits need honest internal-spec citations, not invented
  scientific references. Bound ids, messages, subjects, coordinates and citations.
- Replace the prefix-only product field with `steps: Vec<StepOutcome>` indexed by
  original plan step. Each outcome retains `step_index`, an explicit state,
  required nullable product and bounded direct blocking dependencies. Suggested
  states: reconstructed, rejected, withheld, blocked and not evaluated. A
  reconstructed diagnostic/reference view still carries its actual incomplete
  or reference-only authority; product presence does not imply usable chemistry.
  Remove old prefix fields and update callers directly. Do not retain parallel
  compact products as a compatibility facade.
- Keep private tracked/material slots at the same original indices. A missing
  physical state on a reconstructed reference-only step is different from a
  failed/unexecuted step. Never reconstruct tracked authority from report JSON,
  metadata, target agreement or another step's successful material.
- Preserve the existing overall `Status` meaning: any established invalidity is
  invalid; otherwise unresolved/blocked/resource-limited evaluation is incomplete;
  consistency requires the final designated step, independent-target agreement
  and all relevant required checks under the explicit assumptions. Also report
  a target match when a surviving final branch matches despite an unrelated
  invalid branch; that local pass must not erase the separate invalidity.

The same current `DesignPlan` and independent target remain the replay inputs.
The report is derived evidence, never an agent verdict or a second source of
physical authority. Do not allow product success to overwrite failed findings.

## Execution refactor and admission

1. Separate **resource/representation admission** from per-input/per-operation
   semantic checks. Keep `material::capsule::plan_preflight`, model guards and
   nested-origin admission before hashing, cloning or allocating dependency
   vectors. Unbounded or undecodable requests may stop globally with an explicit
   admission finding. A stale input or bad operation profile should block its own
   consumers while unrelated admitted inputs/branches are still checked.
2. Enumerate typed sources once after their lengths are admitted. Record unknown
   ids, stale bindings and forward/cyclic references as named failures, not a
   fabricated material. Ambiguous duplicate input ids cannot be resolved by
   choosing whichever map insertion won. Dependent steps retain direct producer
   indices; do not copy full transitive blocker paths into every descendant.
3. Extract one current per-step execution method returning a staged result.
   Existing guards construct typed failures at their callsites; domain calls use
   an explicitly named check in `map_err` or their typed native outcomes. The
   original `Error` text can be message detail. **Never parse error strings into
   rule ids**, and do not use one generic reconstruction code for all normal
   primer/end/profile/dependency failures.
4. Evaluate independent cheap checks together. Block dependent computation when
   its prerequisite is unresolved. Commit each step's validated product, tracked
   state and related evidence atomically. Current arms append evidence/findings
   before later guards can fail: retain useful failed/withheld evidence as such,
   but do not leak a ghost product, a successful trace or physical authority.
5. Continue later independent nodes in deterministic plan order. Shared cleavage,
   materialization, evidence and finding budgets apply to the whole job, including
   failed branches. Budget exhaustion means explicit not-evaluated/indeterminate
   outcomes for affected remaining work, never scientific incompatibility or
   silently completed enumeration. Do not retry a failed expensive branch or
   reset its work budget to keep collecting errors.

Retain current caps: 1,000 plan steps/inputs, per-operation source bounds,
32M accumulated product/material bases, 32 MiB step evidence, 64 MiB combined
serialized products/evidence, bounded nested origins and the shared 100M cleavage
comparison budget. Source ancestry replay retains its separately declared bounded
scope. Add a bounded total check count and serialized check budget before finding
copies; reserve capacity for step blocking/limit outcomes. Continue streaming
byte admission instead of serializing arbitrary large diagnostic strings first.

## Direct callers requiring the same batch

| Owner | Required direct change |
| --- | --- |
| `verification/output.rs`, `output/pdf.rs` | Render typed rule/verdict/citation and select products by actual step; escape all displayed text. |
| `cloning/assembly/{recipe,diagnostics,junctions}.rs`, `cloning/ligation/recipe.rs`, prepared-pool and primer/domestication owners | Replace prefix length/first-failure inference with explicit per-step outcomes. Keep every current physical/reference-only prerequisite. Diagnostics should expose multiple failed rows and blocked descendants. |
| `material/{lineage,run,capsule}.rs`, `material/derivation` | Replay and compare step-indexed products. Saved complete runs still require every intended step/event; failed independent drafts must not become reusable saved material. Preserve exact source events and current fresh replay. |
| `planning/combinatorial/{mod,batch,export,release}.rs` | Count evaluated/reconstructed/failed steps explicitly; require designated final product and complete authority for release. |
| `protocol/{schema,validation}.rs`, worker response and native candidate envelopes | Current version bump, typed schemas and bounded decode; documentary shape validation must reject duplicate/missing/incorrect step indices and state/product contradictions. Fresh request binding remains separate. |
| MCP `tools/assembly/intents.rs::verify_plan_inner`, assembly/cloning/workspace tools | Return original step ids with typed checks. Its current enumerate-products output incorrectly becomes sparse step identity if unchanged. Expose full bounded report artifact even for failed ordinary plans, not only plans with specialized scientific evidence. |
| Web `cloning/preparation/review`, assembly reviews/models and `history/saved_run_review_*` | Render failed independent rows and precise blocked descendants, typed check links, product lookup by actual step; retain current source/epoch/generation guards. |

Existing saved report fields that intentionally summarize documentary notes must
be reviewed individually. Do not broadly convert typed findings back into strings
at every boundary, preserve obsolete internal encodings or add hidden fallback
loaders. Bump the affected current envelopes and reject superseded formats.

## Practical parallel edit partition

First freeze `CheckFinding`, `StepOutcome`, state/binding rules and budget APIs in
one short source design message. Then use four mutually exclusive source owners:

1. **Execution/contract owner:** `verification/{engine,evidence,sources,mod}` and
   necessary cohesive new execution subowners. Own graph evaluation, slot identity,
   per-step staging and domain failure sites. The current 2,127-line engine merits
   a responsibility split by execution/operation domain if this batch grows it;
   never split a match into numbered fragments or introduce forwarding facades.
2. **Scientific helpers/native consumers:** `verification/{pcr,terminal,projection}`,
   assembly/ligation/primer diagnostics, material replay and combinatorial owners.
   Emit typed findings directly and update product/material consumers against the
   frozen contract. Keep `verification/output.rs` with this owner unless explicitly
   assigned elsewhere.
3. **Transport/UI/export owner:** protocol schema/validation, MCP tools, web reviews
   and PDF/HTML presentation; update current versions, schemas, bindings and artifact
   summaries together. Own all affected shared protocol files to avoid overlap.
4. **Independent fixtures/security owner:** focused current verifier and integration
   tests, bounded report admission/tampering fixtures and documentation; no edits to
   production files claimed above. Preserve Cargo discovery and test inventory.

Required literal cases: two failing independent PCR/digest branches plus a valid
third branch; a consumer of each failure explicitly blocked; a later valid branch
whose step index differs from compact success count; multiple cheap errors in one
step; reference-only success still unusable as chemical substrate; stale input
isolated to its consumers; target mismatch and unrelated failure both visible;
withheld TOPO evidence without a product; a failure after staged evidence creation;
shared-budget exhaustion without reset; report slot/blocked-index tampering; and
identical successful replay products/annotations/material under the same current
plan and source bindings. Finish all source owners before one combined native/WASM
compiler batch; fix concrete diagnostics together.
