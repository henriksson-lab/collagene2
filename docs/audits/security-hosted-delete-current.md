# Hosted MCP stale deletion

Integrated verification: native/WASM compilation and the connected **490-case**
focused batch pass, including the relevant admission/owner/search regressions.
[Exact case list and scope](features-security-current-integration.json).
Application HTTP, browser and live-provider acceptance remain unexecuted.

Source review, 28 November 2026. Confirmed lifecycle/logic race in the new hosted
MCP implementation, found before publication. Impact is loss of the same
member's replacement scientific workspace; this does not bypass authentication
or access another member's workspace or the shared plasmid database.

## Failure path

`hosted_mcp::serve` captured a quota lease for workspace A, then awaited the
request body. A concurrent DELETE could retire A and a POST could create a new
workspace B under the same member/lab/workspace label. The delayed DELETE then
called key-only `Inner::remove`, retiring B despite having been admitted for A.
The later SDK cancellation guard did not protect the earlier DELETE branch.

## Resolution

`WorkspacePool::retire` now requires the exact `Arc<Workspace>` captured by the
request lease. It compares pointer identity under the same pool mutex as removal
and returns Conflict if the current entry is a replacement. `Inner::remove` and
the HTTP caller pass the captured owner directly; there is no key-only fallback.
Normal retirement still closes compute admission and retains quota/storage until
both response leases and running native work finish.

The pure owner regression
`delayed_delete_cannot_retire_a_replacement_under_the_same_workspace_key`
constructs that interleaving and requires B to remain current and uncancelled.
It needs no browser or HTTP request. Test execution and native/WASM checks are
in the connected feature batch. This is one resolved source finding, not
an exhaustive security certificate.
