# Optional browser WebMCP

The open editor registers three tools when the browser provides the current
`document.modelContext.registerTool` API. Unsupported browsers retain the ordinary
application. Tauri WebView availability is not assumed. Registration uses the
current AbortSignal lifetime API; there is no polyfill or obsolete API fallback.

| Tool | Behavior |
| --- | --- |
| `collagene_current_context` | Returns the active editor document name, alphabet, topology, reference length, selection, selected feature index, feature count and exact source pin. Does not return the entire sequence. |
| `collagene_read_range` | Reads at most 10,000 reference symbols from an exact pinned event. |
| `collagene_navigate_range` | Uses the existing Rust editor navigation callback to select and highlight an exact pinned reference range. |

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

This increment does **not** provide design proposals, reviewed apply, scientific
verification, physical-material reconstruction, draft editing, storage, exports
or unattended automation. It therefore does not complete the broader WebMCP design
integration in the project plan. Ordinary MCP remains a separate interface.

Pure model fixtures are in `web_mcp::tests`. Browser registration and Tauri runtime
availability have not been exercised in this implementation batch.

The API follows the [current WebMCP community draft](https://webmachinelearning.github.io/webmcp/)
and [Chrome's WebMCP updates](https://developer.chrome.com/docs/ai/webmcp/updates).
The community draft is an evolving API, not a W3C Recommendation.
