# Interactive annotation labels

The Rust WebGL map keeps source geometry, label text and camera state separate.
`annotation_geometry` validates borrowed location bounds and prepares local parts
once per document revision. The sequence-row adapter uses the same preparation.
Original feature indices survive `source` omissions and unknown locations.

The independent MIT `collagene-map-label-layout` crate places measured text boxes
in disjoint side regions. Every result is placed or explicitly omitted, with
bounded sorting and storage. It does not measure text, alter annotations, validate
biology or replace the separate SVG/PDF export layout.

The web adapter caches full names, measurements, colors and an anchor in an actual
local part. Joined annotations use the longest known part, retaining source order
for ties; their anchor never averages across an unannotated gap. There is one
label per annotation, independently of the number of rendered parts.

A maximum of 128 candidate identities is selected across source-coordinate ranks,
with the selected annotation given presentation priority. The library selects
that budget when the source/selection or rendered geometry set changes. Rotation
and resize project and place only those cached candidates; they do not rerun
scientific analysis, clone source names or upload static WebGL geometry.

Circular labels use left/right regions; linear labels use regions above/below the
reference. Names can wrap onto two measured lines or end with a measured ellipsis.
The interface reports omitted and abbreviated labels. Hover and the annotation
inspector retain the full original name. Clicking a text box selects its original
feature index, including when two annotations have identical coordinates.

The full source and annotation inspector retain annotations omitted from the map.
Unknown positions count as unplaced. Label-cache count, aggregate name bytes
(4 MiB), individual name bytes (4 KiB), finite text metrics, a 16-pixel glyph line budget
and location limits fail
explicitly rather than silently displaying an arbitrary feature prefix. Text
boxes have opaque map-background fills, drawn after all leaders; zoomed map
geometry may be occluded by a label box. Leaders can cross. Placement is bounded
greedy layout, not an optimal solution or a guarantee of temporal stability.

**Split**, **Sequence** and **Map** change the displayed panels while keeping them
mounted. View controls and map transforms do not edit molecular data or undo
history. Enlarging the map makes more label space available. This does not
complete the planned resizable split panes, shared interactive/export scenes,
context restoration or cross-browser performance acceptance.

Development checks live in the layout crate and
`scripts/map_labels_browser_qa.py`. The browser helper imports actual GenBank
records, checks numeric hit-box bounds/collisions and exact original identity,
exercises all three views and uses a real 1 Mb/10,000-annotation record. Functional
scale checks do not establish a portable frame-rate or memory benchmark.
