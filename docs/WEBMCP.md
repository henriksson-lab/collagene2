# Optional browser WebMCP

The open editor registers seven tools when the browser provides the current
`document.modelContext.registerTool` API. Unsupported browsers retain the ordinary
application. Tauri WebView availability is not assumed. Registration uses the
current AbortSignal lifetime API; there is no polyfill or obsolete API fallback.

| Tool | Behavior |
| --- | --- |
| `collagene_current_context` | Returns the active editor document name, alphabet, topology, reference length, selection, selected feature index, feature count and exact source pin. Does not return the entire sequence. |
| `collagene_read_range` | Reads at most 10,000 reference symbols from an exact pinned event. |
| `collagene_navigate_range` | Uses the existing Rust editor navigation callback to select and highlight an exact pinned reference range. |
| `collagene_assembly_context` | Reads the active declared assembly's source events/selectors, unfinished method settings, literal edits, draft pin and current intended-reference review state. |
| `collagene_open_assembly_workspace` | Opens the existing assembly controls without changing the draft or target. |
| `collagene_review_assembly_draft` | Invokes the existing “Review intended construct” worker action for the exact draft pin. Returns the running review token; completion is read through assembly context. |
| `collagene_cancel_assembly_review` | Invokes the same cancellation command as the visible controls, only for the exact running owner/generation token and draft pin. |

Read and navigation arguments are `{expected, start, end}`. `expected` is the
complete `pin` returned by current context:

```json
{
  "workspace_epoch": 7,
  "source": {
    "document_id": "document ID from current context",
    "event_id": "event ID from current context",
    "revision": "revision from current context"
  }
}
```

Coordinates are zero-based half-open reference sequence boundaries. The tools
refuse empty, reversed, oversized or out-of-range requests. Circular origin
crossings need separate explicit ranges; no wrapping or clamping occurs.
Reference coordinates do not represent extruded bottom-strand material calls.

The tools read the live editor, including unsaved changes. They refuse access
while workspace restoration is in progress. Every read/navigation operation
checks document, event, source revision and workspace epoch together. Equal DNA
at another event does not satisfy the pin. Selection changes do not mutate the
sequence, history or saved catalog. Responses are bounded to 64 KiB, requests to
16 KiB. Names are explicitly untrusted content; long names have an explicit
truncation flag. No cross-origin exposure is requested.

Editor disposal aborts registration. Execution checks both registration lifetime
and the browser's per-execution cancellation signal before entering Rust. Each
operation is synchronous after admission; cancellation does not claim to undo an
already applied view selection. A registration error aborts the whole tool group.

## Assembly review connection

Read `collagene_assembly_context` with `{}`. Its pin is
`{workspace_epoch, draft_revision}`. It binds the source event/revision, selector,
material selection and cached origin identity, physical preparation orientation,
raw method settings, literal edits, and active library declaration. It does not
serialize or replay source molecules or closed material capsules. The active
construct is returned, rather than every construct in an alternatives library.
Material identities are returned separately from capsule-free method settings;
those settings cannot be used to manufacture physical material.

To start review, open the workspace using `collagene_open_assembly_workspace({})`,
then read assembly context after the controls mount. Pass its exact pin as
`{expected: pin}` to `collagene_review_assembly_draft`. This invokes the same
native worker path as the normal “Review intended construct” button. The immediate
response means work has started; it is not scientific completion or acceptance.
Poll assembly context. `intended_reference_ready` supplies the independent
requirement's target identity and intent context, **not** a verified cloning plan.
When the draft changes, that old completion is not returned as a current review.
Captured historical source events stay pinned until explicitly replaced in the UI;
new document heads do not silently replace them.

Cancel using `{expected: pin, review: token}` where the token is the exact current
`{owner, generation}` from review status. Tokens cannot cancel another mounted
workspace or a later review. Closing the dock disposes the endpoint and its worker.
As in a normal button click, a successfully admitted request starts a shared UI
workflow: cancellation of the short WebMCP request after its running receipt does
not cancel that workflow. Use the explicit cancel command or normal UI control.

The reviewed intended reference exposes the existing primer proposal and apply
controls. Agents do not bypass those controls, supply a Frozen receipt, apply
primers, accept a candidate, or assert method success. Existing proposal freshness,
source attestation, explicit apply and subsequent verification requirements remain
in force. Automated proposal/application, draft mutation, scientific verification
commands, storage, exports and unattended automation remain outside this adapter.
The broader project-plan acceptance for agent-driven design proposals/apply is
therefore still open. Ordinary MCP remains a separate interface.

Pure model fixtures are in `web_mcp::tests` and `web_mcp_assembly::tests`. Browser registration and Tauri runtime
availability have not been exercised in this implementation batch.

The API follows the [current WebMCP community draft](https://webmachinelearning.github.io/webmcp/)
and [Chrome's WebMCP updates](https://developer.chrome.com/docs/ai/webmcp/updates).
The community draft is an evolving API, not a W3C Recommendation.
