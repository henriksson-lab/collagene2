# Current agent verification contract

2026-11-28. Current source passes native/WASM checks and 1,860 unique selected
cases across the integration batch. [Execution/source evidence](verifier-and-design-current-integration.json).
This describes the sole current representation.

## Ownership and authority

`planning/plan.rs` owns DesignPlan, Operation and Source. The existing editing
owner counts and edits source references. There is no second execution plan,
prose-to-rule classifier or superseded decoder.

`verification/engine.rs` is the entry. `engine/execution.rs` owns borrowed
admission, independent branch scheduling, atomic commits and target comparison.
`engine/operation.rs` reconstructs products and assigns rule identifiers at the
actual guards. Native source, PCR, terminal and annotation owners retain their
computations. Language models propose designs; they never supply verdicts or
chemical authority.

## Indexed outcomes and named predicates

VerificationReport retains planned_steps and exactly that many ordered
StepOutcome records. Each has its original step_index, a state, a required-nullable
product and direct unavailable dependencies. States are Reconstructed, Rejected,
Withheld, Blocked and NotEvaluated.

Products and private physical/reference-only state occupy the same original
slots and commit together only on success. Consumers resolve exact indices.
A failed branch cannot shift another product into its slot. Dependent steps
remain blocked while independent branches continue. Useful withheld evidence
can be retained without creating a product or usable material.

CheckFinding carries a RuleId, verdict, plan/input/step/target scope, source
subjects, original coordinates, citations and explanatory text. Verdicts are
Pass, Fail, Indeterminate, Informational and Blocked. Named predicates are
published at collagene://verification-rules/v1; domain citations supplement that
contract. Rules are never inferred from an error message.

An established failure remains Invalid even when an independent final product
matches the target. Missing final products never cause comparison with an earlier
product. A wrong target binding leaves independent operations available but blocks
comparison and rejects target-dependent CheckedTopo acceptance.

## Admission and shared budgets

Borrowed resource admission precedes hashing and reconstruction. Malformed or
excessive request shapes are refused before expensive work. Bounded semantic
input identity/revision, source selection and operation failures are checked
locally. A bounded wrong saved-origin selection cannot invalidate unrelated
branches globally.

Shared limits include 1,000 steps/inputs, 100,000 source references, 8,192 findings,
4 MiB finding evidence and a 64 MiB report. Optional detail cannot consume reserved
outcome evidence. Material, cleavage and native method budgets remain explicit.
Exhaustion cannot produce a truncated success.

ExactSearchBudget counts actual exact-PCR preprocessing, matching and copying
across successful and failed operations under one 200M-unit ledger. Native source
PCR uses the same kernel and retains its separate materialized-byte usage result.
Intentional-mismatch and other embedded methods retain their declared work scopes.
Portable candidate exchange preserves import-scale original records; scientific
replay retains tighter computation admission.

## Consumers and verification scope

Assembly, ligation and tailing retain nullable original-row products. Diagnostics
retain multiple failures and independent successes. Complete-run publication still
requires every intended product and fresh private authority. Saved-origin notes
remain documentary text, distinct from typed checks.

MCP always retains a complete bounded verification artifact, including ordinary
failed plans. Compact output preserves original indices, states, blockers,
nullable products and typed checks. Web and HTML/PDF use the same exact slots.
Product navigation skips missing products without changing their identity.

Current versions are candidate34, workspace plan32, editable draft23, material
derivation5, batch6/contract7, release4 and science38. Worker135 also contains
independent IVT/base-editor analysis additions. Unsupported internal versions are
rejected explicitly.

Fixtures cover simultaneous failures, later independent products, direct blocked
descendants, source/target binding, saved origins, actual shared PCR work,
documentary tampering and before-hash resource admission. Authored fixtures are
not execution evidence; results belong to the completed integration record.
No browser, application HTTP, screenshot or live-provider checks run in this
batch. Full feature and graphical acceptance remain open.
