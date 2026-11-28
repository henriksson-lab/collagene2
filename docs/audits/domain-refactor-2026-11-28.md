# Domain ownership refactor checkpoint

28 November 2026. Implementation and local native/WASM verification are complete. Exact-commit cross-platform CI and desktop packaging passed at `416b26c`; full-product and visual acceptance remain open.

| Area | Current ownership |
| --- | --- |
| Scientific core | `src/lib.rs` declares 14 real domains. Scientific families live under sequence, alignment, chemistry, cloning, documents, formats, material, output, planning, protein, protocol, review, tracks and verification. The model compatibility facade is deleted; callers use `collagene-model` directly. |
| Plans and verification | `planning/plan.rs` owns `DesignPlan`, `Operation`, their source visitors and plan revision. Private verifier `Execution` owns the correlated report/physical prefix and reconstruction bounds. Evidence, source guards, projections and preview owners have explicit modules. |
| Worker transport | `protocol` separates schema, bounded envelopes, request/response types, stateless execution, persistent sessions and binding checks. Encoding and decoding belong to `Request` and `Response`. |
| Protocol output | `verification/output.rs` owns non-deserializable `CheckedProtocol`, created from fresh verification with private fields. HTML and PDF use that current checked context; imported reports cannot construct it. |
| MCP | `server` owns transport/router composition and resources; `session` owns bounded storage, references and atomic publication; `tools` owns semantic request schemas and their actual handlers. The current registry remains 208 tools. |
| Web | App composes feature controllers and rendering. `WorkspaceSession` owns six document/workspace handles and guarded admission/edit/export context. `PcrState` owns PCR lifecycle/evidence handles; `PcrController` owns PCR actions. `ExportArtifact`/`ExportService` own bounded download reservation, guarded completion and lifecycle cleanup. |
| Tests | Core and MCP integration tests use nine and eight explicit domain harnesses with nested module registration. Fixture paths and current API imports have been updated directly. |

The refactor does not introduce another internal plan, migration path or compatibility decoder. The refactor checkpoint retained scientific engine `science-32`, worker schema `121` and candidate schema `31`. Mechanical moves preserve serialized enum ordering and scientific assumptions; the local checks and independent review below verify those contracts within their stated scope.

## Why some implementations remain distributed

The rmcp framework associates tool metadata, schemas and router generation with methods implemented on `CollageneServer`. Domain-specific server implementations therefore remain in the corresponding `tools` modules, with their request types and real handler bodies. Named routers are composed by the server. This is a framework boundary with a concrete registration purpose, rather than an arbitrary split of stateful methods or a root of forwarding wrappers.

Mutable MCP storage remains cohesive: `Store` and `Session` own reservation, quotas, reference validation and atomic document/workspace/artifact publication. Genome storage has its own resource owner. Moving those mutations into the distributed tool implementations would weaken the transaction boundary; tools invoke the session owner instead.

Web controllers likewise own narrow features rather than collecting every Leptos signal into one application state bag. Document/workspace, shell presentation and viewer selection contexts have explicit small scopes; PCR lifecycle state remains separate. Assembly, restriction, import, library, protocol, annotation, comparison, reference, biology and read controllers own their guarded callbacks. App retains cross-feature composition and restoration application because those coordinate the whole view. A universal `AppController` owning all state would only relocate the monolith.

## Layout enforcement

`python3 scripts/check_rust_layout.py` is enabled in the existing `web` job in `.github/workflows/check.yml`. It examines tracked and unignored Rust source files, including tests and build scripts. It rejects more than 20 direct Rust files per directory and files above 5,000 lines, and reports ownership-review warnings from 3,000 lines. There are no source grandfathering exceptions.

A static inventory at this checkpoint found a maximum of 16 direct Rust files in a directory and no file at or above 3,000 lines. The largest observed files were worker execution (2,977 lines), App composition (2,943) and material derivation (2,757). The final gate reports zero violations and zero warnings across 1,202 Rust files. File limits complement review of ownership, method cohesion and resource authority.

## Remaining acceptance

Exact-commit Windows/macOS/Linux CI and desktop packaging passed for this refactor tree. Installed
platform and full-product/visual acceptance remain separate product gates;
browser/application HTTP checks stay deferred under the user's instruction.

The original audit documents describe pre-refactor findings. The executed move
plan has been deleted; current module declarations define the implementation.
The [independent verifier review](verifier-ownership-review-2026-11-28.md)
accounts for all 49 original functions and 27 types, including execution order,
resource bounds, partial results, revision hashing and fresh protocol authority.


## Final local verification

Native harness compilation and web/worker WASM checks pass. All 1,935 integration
cases and 401 native web unit cases match executable discovery. Execution passed
1,034 cases: core verification 175, material 220, protocol output 49; MCP SDK 180;
MCP private/session/catalog 9; web 401. Assertions were retained; no browser,
application HTTP, screenshot or release-build checks were run.

Formatting, whitespace and local documentation-link checks pass. The strict
layout check covers 1,202 Rust files and reports zero violations and zero warnings.
Exact-commit Windows/macOS/Linux CI and desktop packaging passed at `416b26c`. Earlier
intermediate diagnostics above describe the integration process, not remaining
source failures.
