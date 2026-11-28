# Native and MCP ownership audit

28 November 2026. Read-only review of current source outside core/web, including
pending folder recovery. No production changes, compiler/tests, application HTTP,
browser, release build, runtime database or credential access were used.

`python3 scripts/check_rust_layout.py` reported **1,254 Rust files in 217
directories, zero violations and zero 3,000-line warnings**. Outside core/web,
the largest file is `topo/src/lib.rs` (2,315 lines), then desktop catalog backup
(2,069). No urgent native size-driven rewrite is warranted.

The referenced rust_refactor module workflow was read. Its owner-first analysis
is useful; automatic release verification and retained delegators conflict with
repository policy and were not used. Findings are direct source review, not
exhaustive semantic scatter analysis or verification of pending behavior.

## Bounded current queue

| Priority | Finding/action | Required boundary |
| --- | --- | --- |
| Before related test growth | `mcp/tests/analysis`, `mcp/tests/cloning`, and `library-sharing/tests` each have 16 direct Rust files. Group actual incoming domains: read/alignment versus protein analysis, PCR/overlap-extension versus terminal preparation, physical-material closure versus publication. | Preserve helpers in their source modules; explicitly register child modules. MCP has `autotests = false` and eight named harnesses. Compare executable inventories before/after. No empty speculative directories. |
| With next tree-job persistence change | `repository/src/folders.rs:289` has `encoded(job)` validating TreeJob's own version, queue/bounds and serialization. Make it private `TreeJob::encode` in that same owner file during the next substantive edit. Review `job_write` at 305 separately. | Keep Repository atomic checkpoint/transaction methods with Repository. Preserve exact bytes/bounds/CAS. No schema bump for a receiver-only move. This is optional until this owner changes again. |
| With further tree-control growth | `desktop-sync-host/src/lib.rs:359` contains roughly 500 lines of credential/remote/outbox/pull/tree command dispatch, in a 901-line file. Keep now; if more controls arrive, group tree commands with their status/scope/review methods (`:317–358`) into one genuine lifecycle module. | Keep reload, lab validation, cancellation admission, credential-triggered review clearing and final receipt handling central. No per-command file fragmentation or forwarding APIs. Verify recovery/cancellation/stale-review service cases. |
| During substantive owner changes | Broad production parent imports hide dependencies in repository folders, desktop activation and SBOL projection. Narrow them while editing those owners. | No standalone cosmetic rewrite; test-local parent imports are distinct. |

## Resolved findings and justified owners

The historical TreeJob move is already complete: `TreeJob::add_node` is at
`repository/src/folders.rs:175`, `TreeJob::scan_step` at 198. Do not carry that
earlier audit finding as unfinished work. MCP's former root monolith and flat
test directory are likewise resolved by current domain routers/session ownership
and explicit integration harnesses.

- Keep CollageneServer's domain handler impls beside rmcp request schemas and
  framework registration. Server construction/router/compute scheduling live in
  `server/mod.rs`; quota, lock and atomic publication stay in `session/mod.rs`.
  Impl scattering here has a concrete domain reason. Do not create a universal
  mutable context or a new root monolith to improve a scatter metric.
- Repository transaction/job extensions and Library mutation/collaboration
  extensions have real lifecycle boundaries. Their central state must remain
  cohesive with transaction dependencies.
- CatalogBackupService backup/recovery and append-only restart-bound activation
  are distinct lifecycles. Its definition and backup methods stay together;
  `catalog_activation.rs` is justified. The backup file includes substantial
  meaningful tests and is below the warning threshold.
- SyncHost publication/removal modules own separate typed review queues with
  different authority and consumption rules. Similar queue mechanics do not
  justify an untyped abstraction. Finish the pending recovery batch before
  reorganizing those callers.
- SourceClient is defined with its methods in `sources/src/native.rs`; portable
  Source/Config methods are with their types. Native transport versus portable
  URL resolution is a platform boundary. Keep lock-based rate scheduling and
  symmetric graph algorithms free where methods would imply false ownership.
- TOPO Graph, WorkLedger and ExecutionState have cohesive local impls. Preserve
  independent source-audit/replay correctness roles; no size-driven split is
  currently required.

Integrate these decisions into the single current refactor plan. This audit is
not another implementation queue or certification of pending scientific checks.
