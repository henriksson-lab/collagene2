# Code organization refactor status

28 November 2026. The domain refactor is implemented across the scientific core, web and MCP sources. Integrated native/WASM validation and executable rediscovery pass locally; cross-platform CI and desktop packaging passed at refactor commit `416b26c`. The current ownership and remaining checks are recorded in the [refactor checkpoint](domain-refactor-2026-11-28.md).

## Current execution checklist

The [working-tree ownership follow-up](organization-followup-2026-11-28.md)
records the new parallel audit and current action queue. The original domain
move is complete; remaining work concerns worker dispatch, material preview
owners, restoration/proposal lifecycles and natural methods. Current limits pass
across 1,213 Rust files, with worker execution and App composition close to the
3,000-line review threshold. Historical measurements below describe the completed
refactor checkpoint, not the latest working tree.

| Batch | Implementation status | Verification still required |
| --- | --- | --- |
| R1: inventory and ownership | Implemented. Source destinations, current tool/schema identities and test inventory are recorded; real state owners and shared integration boundaries are explicit. | Passed: all 1,935 integration and 401 web unit cases are discoverable; layout gate passes. |
| R2: core owners | Implemented. Model facade removed; canonical model imports updated. Worker envelope/schema/request/response/session owners are separate. `DesignPlan` and `Operation` own identity/source methods. Private verifier execution and checked HTML/PDF output owners replace mixed responsibilities. | Passed: native harnesses, web/worker WASM, scientific/transport and fresh-output authority cases. |
| R3: MCP owners | Implemented. Server/resources, atomic session storage and domain tool schemas/handlers are separate. Named rmcp routers compose the current 208-tool registry. | Passed: all 180 SDK cases and nine private/session/catalog cases; 208-tool registry retained. |
| R4: web owners | Implemented. App composes narrow controllers; workspace, PCR and export lifecycles have concrete owners. Draft identity/comparison behavior lives with `Draft`; old view/model facades are removed. | Passed: final WASM ownership check and all 401 web unit cases. |
| R5: domain and test relocation | Implemented. Scientific/UI families live in semantic subtrees. Core/MCP tests use nine/eight explicit domain harnesses with direct caller/fixture updates. | Passed: exact executable inventory comparison for all 2,336 integration/web cases. |
| R6: final verification and enforcement | Local verification complete. The layout checker is enabled in the CI web job. | Native/WASM, formatting, layout, independent review and 1,034 executed cases pass. Passed: exact-commit Windows/macOS/Linux CI and desktop packaging. |

The refactor preserved scientific version `science-32` and worker/candidate schemas `121`/`31`; subsequent feature batches advance their own current contracts. The current implementation uses one plan representation and one current schema. Unsupported internal versions remain rejected. No old module facades, compatibility aliases or fallback implementations are intended to survive this change.

## Layout policy and current measurements

`python3 scripts/check_rust_layout.py` checks tracked and unignored Rust files, including tests and build scripts. More than 20 direct Rust files per directory or 5,000 lines per file fails the gate; 3,000 lines triggers an ownership warning. Ignored build/dependency caches are not source. This gate has no grandfathering exceptions, and size compliance does not replace review of coherent owners and methods.

The current static checkpoint found a maximum of 16 direct Rust files per directory and no Rust file at or above 3,000 lines. Core now has only `lib.rs` directly in `src`, declaring 14 actual domain modules. Worker execution and web App composition remain close enough to the warning threshold to monitor when adding functionality. The final gate reports zero violations and zero warnings across 1,202 Rust files.

The [test inventory](integration-test-inventory.json) records all 1,935 integration cases and 401 web unit cases. Executable discovery matches every case. The focused run passed 1,034 cases: 444 scientific/protocol cases, 180 MCP SDK cases, nine MCP private/session/catalog cases, and all 401 web unit cases. Source-level web test attributes remain 402; one target-gated case is outside native discovery.

## Ownership review and audit artifacts

- [Current domain checkpoint](domain-refactor-2026-11-28.md): actual scientific, session and controller owners; framework-specific method distribution; remaining validation.
- [Original scientific-core audit](code-organization-core.md): pre-refactor evidence and reasoning, including wire selectors, verifier/material authority and test discovery.
- [Original web audit](code-organization-web.md): pre-refactor controller, draft and method-cohesion findings.
- [Original MCP/native audit](code-organization-mcp-native.md): pre-refactor router, transaction and native-boundary findings.

The executed source move plan has been deleted. The current owner trees are
defined directly by the [core](../../crates/core/src/lib.rs),
[MCP](../../crates/mcp/src/lib.rs), and [web](../../crates/web/src/main.rs)
module declarations. The
[independent verifier review](verifier-ownership-review-2026-11-28.md) accounts
for every original scientific function and type without finding a semantic
regression. Compiler and execution checks remain separate requirements.

The helper at `/home/mahogny/github/claude/rust_refactor` offers useful ownership/scatter and coupling reports. Its automatic release verification and suggested compatibility phases conflict with repository rules; use reviewed read-only analysis unless its write behavior is explicitly adjusted. Continue using coherent batches, direct caller updates and one shared compilation cache. Browser/application HTTP checks remain deferred under the user's instruction.

The refactor commit `416b26c` passed [Check](https://github.com/henriksson-lab/collagene2/actions/runs/46095578568) (web and all three native platforms) and [Desktop build](https://github.com/henriksson-lab/collagene2/actions/runs/46095578560). This evidence covers that source checkpoint; subsequent scientific feature changes require their own validation.
