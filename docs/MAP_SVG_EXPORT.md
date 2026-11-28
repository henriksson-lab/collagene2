# Portable SVG and PNG map export

The canonical map owner is `core::output::map_export`. It accepts a molecule,
its exact revision, optional displayed restriction-site evidence, and explicit
`raster::Format::Svg` or `Png(PngOptions { width })`. The output owns its source
revision, selected image representation and bounded renderer warnings. No older
map request/result or SVG-only internal adapter remains.

Both formats share existing deterministic feature and restriction geometry.
The renderer exports the whole document in canonical orientation with bounded
feature/label layout. It does not reproduce viewport rotation, zoom, range
shading, hidden-feature filters, search/PCR overlays or selected-site highlights.
The exact supplied restriction list survives; export neither rescans nor invents
enzyme cleavage. Map display is documentary reference geometry, not evidence of
experimental cleavage, isolated material or cloning success.

## Pixels and resources

SVG is vector text capped at 2 MiB. PNG is produced entirely in Rust from the same
bounded SVG, using pinned resvg 0.47 with text support and the existing Liberation
Sans font. No CLI/browser, external program, system-font scan, uploaded SVG,
remote image, local path or raster codec is used. Image resource resolvers are
explicitly disabled. Unsupported font scalars are visibly escaped as [U+XXXX]
and reported instead of silently disappearing. Original names remain in source.

PNG width is 16–4096 pixels; height follows the diagram's aspect ratio. Each side
is at most 4096 pixels and total pixels at most 4,000,000 (16MB RGBA). Encoding uses
a 24 MiB byte budget and source revision text metadata. Response validation checks
PNG signature, bounded dimensions, RGBA encoding, chunk order/CRCs, exact source
context and encoded-content hash without pixel decompression on the UI thread.
The worker/native renderer and independent fixtures fully decode output with
bounded row buffers; no untrusted dimensions trigger full-image admission allocation. PNG is pixels rather than searchable vector text or full coordinate
metadata; retain the source document or SVG for complete evidence.

The upstream renderer is a portable Rust library available under MIT or Apache2:
[resvg source and license](https://github.com/linebender/resvg). This is general
rendering infrastructure, not borrowed cloning-software logic.

## Captured web workflow

One `MapExportDialog` offers SVG and PNG with explicit width controls. App takes
an exact source/event/epoch/view snapshot, including actual rotation, zoom,
selection, feature visibility, restriction list/filter and selected-site context.
No synchronous duplicate SVG button bypasses that dialog. Heavy PNG rendering
runs in the worker.

The context also binds image format and pixel settings. Preparation checks the
opening source/view, current context on completion, exact worker ID/revision,
source revision and image format. Save preparation and native save admission
check the complete capture again. A stale completion, changed format/width or
closed dialog cannot publish bytes. Native success is acknowledged for the
captured document/event; later source changes cannot reverse a real completed
save. Filename and MIME match the selected format.

Current worker 126 uses `RenderMap`/`MapExported` at tags 174/164, both under the
ordinary 32 MiB framing budget. Source/model validation remains authoritative;
serialization integrity never grants physical-material authority.

## Predicted gels

`core::output::gel_export` similarly binds a full documentary report/calibration,
its report revision and selected format. Worker `RenderGelImage`/`GelImageExported`
(tags 199/189) independently recomputes that context. The gel pane prepares PNG
in a worker with a displayed-report generation/workspace epoch guard, then
checks it again before download/save admission. SVG/report exports remain
intentional external formats. No observed-image/OpenGel comparison is included.

## Authored verification

Native fixtures independently check PNG backbone/feature pixels and font output,
calibrated 1000 bp gel position, unknown circular mobility, exact source/report
pins, malformed CRC/compressed-stream refusal, dimension/pixel/byte limits and
stateless/session current-format behavior. Web model fixtures bind format and
width as well as source/view. Actual MCP SDK cases download SVG and PNG artifacts,
preserve immutable sources, and reject stale sources/oversized settings. These
fixtures await the parent-managed coherent native/WASM check batch; they are not
browser, screenshot or installed-platform acceptance.
