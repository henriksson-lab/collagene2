# Code organization and refactoring plan

28 November 2026. This is part of the current product plan, not an alternative
implementation or protocol. Three parallel static audits examined the scientific
core, web application, and MCP/native crates. They included visible unfinished
SBOL work; they did not verify that work or run builds/application checks.

## Measured violations

| Directory | Direct Rust files | Maximum |
| --- | ---: | ---: |
| `crates/core/tests` | 320 | 20 |
| `crates/web/src` | 228 | 20 |
| `crates/core/src` | 156 | 20 |
| `crates/mcp/tests` | 95 | 20 |

| File | Lines | Assessment |
| --- | ---: | --- |
| `crates/mcp/src/lib.rs` | 8,973 | Above absolute 5,000 limit |
| `crates/web/src/app.rs` | 7,273 | Above absolute limit |
| `crates/core/src/protocol.rs` | 5,797 | Above absolute limit |
| `crates/core/src/verification.rs` | 4,528 | Above 3,000 review threshold |
| `crates/web/src/assembly_draft_model.rs` | 3,174 | Above review threshold |

Run `python3 scripts/check_rust_layout.py` to refresh the baseline. It checks
tracked and unignored Rust files, including tests and build scripts. It currently
fails on the actual violations; there are no grandfathering exemptions. Ignored
build output and dependency caches are not repository sources. Ownership and
method cohesion require review in addition to this mechanical check.

## Execution order

Implementation checkpoint: the model facade has been deleted and callers now
depend directly on `collagene-model`. Worker transport is split into nine
cohesive files; the largest is below 3,000 lines. Request encoding and response
decoding belong to their types. Core/MCP tests have moved into nine/eight explicit
domain harnesses. Cargo metadata confirms registration; the
[test inventory](integration-test-inventory.json) records all 1,931 preserved
cases plus four new SBOL SDK cases. All 17 integration harnesses compile and
their executable test lists match all 1,935 recorded cases. Focused native
execution passed 158 tests. Draft production code is now 1,594 formatted lines;
its identity/comparison methods live with `Draft`, its workspace facade is
removed, and all 33 extracted unit cases pass. The web/worker WASM check passes.
Remaining work includes verifier owners, the flat core/web source trees, and
the oversized web/MCP roots. The exhaustive
[core move manifest](core-source-moves.json) assigns all 164 current source files
exactly once and projects a maximum of 11 Rust files per directory; it has not
yet been applied.

This refactor takes immediate priority over adding features. Preserve the full
product scope, unfinished work and scientific behavior. Each batch finishes all
callers before its focused background compilation; intermediate builds may break.

| Batch | Work and owner | Completion evidence |
| --- | --- | --- |
| R1 | Record source destinations, current tool registry/schema and test discovery. Settle session, draft and plan ownership. One integrator owns central declarations, Cargo targets and cross-crate changes. | Every existing source/test has a destination or justified deletion; projected directories have headroom, normally at most 16 files. |
| R2 | Core: replace the model compatibility facade with direct canonical model imports; split worker envelope/schema/session/domain dispatch; colocate `DesignPlan` and `Operation` behavior. Extract a private verifier execution owner and checked HTML/PDF context. | No facade or second plan/schema; mechanical moves preserve current wire bytes/enum ordering. Fresh deterministic verification remains the authority for protocol output. |
| R3 | MCP: separate server/resources, atomic session storage and semantic tool routers. Move schemas and actual handlers together, using named rmcp routers rather than root forwarding methods. | Registry names/schemas preserved, no omissions/duplicates; quota reservation and workspace/artifact publication remain atomic. |
| R4 | Web: introduce cohesive workspace/document/selection and assembly controllers; move the complete export lifecycle out of `app`; return draft identity/mutations to the draft owner. Extract shell composition and remove model/view reexport facades. | App is composition rather than a universal state container; exact source/event/epoch guards, cancellation and bounded exports remain intact. |
| R5 | Relocate remaining scientific and UI families into semantic subtrees. Regroup core/MCP integration tests into domain harnesses with explicit module registration; update target selectors, helpers and fixture paths directly. | No directory over 20 files; no file over 5,000 lines. Before/after discovery inventory accounts for every meaningful test. No obsolete root wrappers. |
| R6 | Finish focused native/WASM validation, repair concrete failures together, then enable the layout checker in existing CI. Review every remaining 3,000-line warning and owner spread. | Focused checks pass, layout gate passes without exemptions, platform CI retains Windows/macOS/Linux coverage. Browser checks remain deferred under the existing user instruction. |

R2–R4 can be prepared by separate workers once interfaces and file ownership are
agreed. R5 family moves can run in parallel on disjoint trees. Workers must not
rewrite shared roots or trigger separate builds; central integration and one
consistent Cargo package/feature selection prevent conflicting edits and rebuilds.
Native crates below the thresholds stay intact unless a concrete ownership fix
helps the current batch; this is not a blanket rewrite.

## Detailed findings and destination maps

- [Scientific core](code-organization-core.md): transport selector constraints,
  verifier/material ownership, semantic source families and integration harnesses.
- [Web](code-organization-web.md): allocation of all 228 current files, controller
  boundaries, draft cohesion and natural method candidates.
- [MCP and native](code-organization-mcp-native.md): supported router composition,
  transaction ownership, test harnesses and native boundary assessment.

The helper at `/home/mahogny/github/claude/rust_refactor` has useful owner/scatter,
coupling and method-discovery analyses. Its documented automatic write workflow
uses release verification and allows temporary compatibility imports/delegators.
Those behaviors conflict with repository rules. Use reviewed read-only reports or
safe dry runs; apply direct coherent edits unless the helper can obey our rules.
Do not split owners into numbered chunks, scatter a large struct's impls across
many arbitrary files, or use compatibility reexports to keep old imports working.
