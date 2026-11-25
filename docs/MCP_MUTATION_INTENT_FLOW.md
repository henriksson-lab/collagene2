# Source-bound mutation intent before inverse PCR

`build_mutation_intent` accepts an immutable `document: { id, revision }` and
`request: { name, edit_start, edit_end, inserted }`. Coordinates are zero-based
source boundaries, start inclusive and end exclusive. Equal boundaries insert;
empty inserted DNA deletes. The target is constructed from the operator's edit
on the original reference **before** invoking PCR, primer design or ligation.
It is an intended reference, not measured material or a prepared PCR substrate.

The response contains `state: "intended-reference"`, `verified: false`,
`preparation_simulated: false`, the exact source/request/context, an immutable
`target` DocumentRef, target summary, deleted length, source-origin shift,
warnings and a full JSON evidence artifact. Source data remains unchanged.
Repeated exact requests reuse the immutable target/artifact. Changing the name,
edit, source revision or input DNA changes context; equal output DNA does not
make two source histories interchangeable.

For an ordinary internal edit, source origin is preserved. A circular edit that
crosses the stored origin uses the core module's explicit rotated-origin rule;
the returned shift must remain visible to clients. Neither coordinate convention
selects a PCR geometry toward a favorable target.

## Practical client flow

1. Import/select an original circular DNA source and retain its exact DocumentRef.
2. Declare replacement/insertion/deletion boundaries and bases. Call
   `build_mutation_intent`; review the literal intended target and edit evidence.
3. Pass the **original source** as `document`, returned intended target as
   `target`, and explicit `InversePcrOptions` to `design_inverse_pcr`. Primer
   lengths, inserted bases and tail split are deliberate preparation choices.
4. Inspect its actual ordered primers and typed PCR → phosphorylation → ligation
   plan. Use `verify_plan` on the saved plan revision; mismatch or incomplete
   reconstruction remains visible. Intent construction alone has no plan status.
5. Existing HTML/PDF plan export and subsequent cloning-output finalization use
   that reviewed saved plan and independent target. No read support or
   experimental success follows from computational consistency.

A generated intended reference must not be reused as an original physical PCR
substrate merely because its DNA is concrete. It retains reference-only evidence;
use the original source or reconstruct a supported physical preparation instead.
There is no metadata-based promotion, silent end repair or automatic reagent/order
assumption in the intent transaction.

## SDK fixture scope

`crates/mcp/tests/mutation_intent_tools.rs` authors two actual in-process MCP SDK
flows using a 48-base circular GenBank source. Literal replacement, insertion and
deletion targets are asserted before any simulator runs; unchanged original export,
immutable deduplication, changed context and stale references are covered. A
separately invoked inverse-PCR preparation uses the frozen replacement target and
explicit 8-base domains, then checks its saved three-step plan. Another intended
edit is supplied to that same preparation to expose a wrong-target failure.
Reference-only input reuse must not produce a consistent physical plan.

These are source-authored fixtures, not a claim of executed test results. The
parent integration task coordinates compilation and execution.
