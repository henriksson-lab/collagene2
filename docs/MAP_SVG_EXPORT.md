# Portable map SVG export

`web/src/map_export_model.rs` and `web/src/map_export.rs` implement the
**Export map SVG** action in the existing File menu. The dialog prepares its
artifact in the Rust worker before saving it. Source integration is complete;
the parent-managed compiler/native fixture batch is still pending. No browser,
renderer screenshot or application check has been performed.

The output uses the existing Rust `svg::render_map` or
`svg::render_map_with_restriction_sites` implementation. It exports a reference
map, not an experimental cleavage prediction or material/provenance archive.
No browser canvas capture, third-party renderer, external program, raster image,
script, external font or additional dependency is introduced.

The native renderer has a distinct portable layout: canonical orientation,
whole-document extent, all source features subject to its documented caps,
32 feature labels and 64 printed restriction labels. Its warnings describe
omitted or unresolved annotation geometry. Viewport rotation, zoom, selection
shading, hidden-feature filters, search/PCR overlays and selected/hovered-site
highlighting are not reproduced. These limitations are visible before saving.
The exact displayed restriction-coordinate list is passed through when supplied;
there is no new scan or silently substituted enzyme set. Native workspace export
remains the route for retaining full annotation and biological history.

## Source and worker contract

The model is registered for native fixture tests and WASM; the dialog is
WASM-only. `MapExportDialog` receives:

- `snapshot: Callback<(), Result<map_export_model::Snapshot, String>>`;
- `render: Callback<RenderRequest, map_export::LocalResult<Rendered>>`;
- `on_saved: Callback<String>`.

App supplies the pinned snapshot and `worker_renderer()` supplies the isolated
worker adapter. The renderer's new read-only `export_view()` getter returns
actual camera state, current highlight/feature visibility, restriction evidence
and immutable overlay context; it never resets a gesture or substitutes a camera
default. An active drag refuses capture until the gesture finishes.

The snapshot callback must reject workspace restoration/unloaded state and bind
its opening exact document/head event using App's export-context guard. It then
returns the unchanged full molecule, exact `EventBinding`, workspace epoch,
reference selection, actual finite renderer rotation/zoom, selected feature,
sorted hidden-feature indices and presentation-settings context. Restriction
state is `None` when hidden, or the exact displayed list plus its source revision,
scan/filter context and selected/hovered-site index. Do not fill actual rotation
with a default. App refuses capture while source, selection, visibility or
restriction-filter
effects are still updating, and refuses stale/missing restriction scans rather
than substituting an empty complete scan.

The dialog binds its opening capture. Source/event/epoch/selection/view changes
require reopening it; a response for a different capture is refused even when DNA
is identical. Snapshot capture occurs before rendering, on render completion,
before export preparation, and again at native save admission. The native file
helper returns acknowledged completion. A source change after acknowledgement
never turns a real saved file into a claimed failure: the completion message names
the captured document/event. The parent callback should use live/try-set guards
for its own UI lifetime, but must not require that captured source is still active.
Browser completion retains the existing browser-download dispatch semantics.

Worker93 appends `RenderMapSvg { request: core::svg_export::Request }` and
`MapSvgExported(core::svg_export::Rendered)`. The request contains the unchanged
molecule, exact source revision and optional displayed restriction-site list.
The response contains its source revision, SVG bytes as UTF-8 text and bounded
warnings. The outer request/reply revision is the complete capture context;
`worker_renderer()` checks both it and the request ID before accepting the typed
result. `core::svg_export::render` is the shared native renderer adapter used by
both the worker and pure model fixtures. No rendering runs synchronously in an
App event handler. This adds no scientific operation or material authority.

Pure fixtures exercise equal-content wrong document/event refusal, workspace
change, rotation/zoom/selection changes, stale molecule event binding,
restriction-filter changes with identical empty displayed lists, crossed worker
responses, filename binding and invalid numeric view state. They await the
parent-managed combined native test batch.
