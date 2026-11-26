# Direct sequence editing and shared restriction sites

Implementation checkpoint, not a browser acceptance claim. Browser verification
is deferred until the in-scope code is complete, as required by the project plan.

The sequence canvas handles pointer drag selection with timed edge scrolling on reference/complement rows,
Shift-click extension, insertion carets and Shift-arrow selection. Typing replaces
a contiguous selection or inserts at the caret; Backspace/Delete remove selected
symbols or adjacent symbols. Ctrl/Cmd+A selects the reference, Ctrl/Cmd+Z undoes,
Shift+Ctrl/Cmd+Z redoes, and Ctrl+Y redoes. Clipboard shortcuts retain their existing
handler. Escape clears the selection. The on-screen hint makes these controls
visible rather than relying on a modal edit dialog.

Exact local typing/cut/paste transitions now preserve their committed caret
through the later reactive source update. Unrelated document writes preserve
active selection; external event changes, undo/redo and restored workspaces clear
it. [Identity and normalization rules](SEQUENCE_EDIT_TRANSITIONS.md).

Edits use the existing model's alphabet validation, annotation transformation and
history commit. Invalid edits preserve the source and report the reason. Compound
or unresolved annotation selections are not silently converted into contiguous
replacements. Ordinary literal GenBank citation ranges can be edited through a bounded,
source-bound recipe. The exact original record remains unchanged; export replays
the declared edits against that source, remaps citation ranges, marks ancestral
citation support unverified and removes stale active accession/version claims.
Unsupported citation forms, deleted citation coverage and unremappable
coordinate-bearing qualifiers still fail explicitly.

The RE sites inspector and editor toolbar control one enzyme selection and site
count filter shared by the sequence and circular/linear map. The bundled catalog
currently contains 40 named enzyme geometries. Map **Select range** mode also supports linear and origin-wrapping circular drag ranges with a local live preview; **Rotate** retains the circular camera interaction. Multi-part selections remain explicit and are not silently edited as one contiguous region.

A scan runs after source changes,
with a short debounce during typing. Filtering and rotating do not rescan DNA.
Every request is bound to the document event context, revision, enzyme geometries
and options; switching sources cancels old work, including equal-content events.

`restriction_inspection::scan()` enumerates recognition evidence without digest
fragment reconstruction. Existing virtual digestion remains a separate operation.
Worker schema 23 appends `ScanRestriction` at job tag 76 and reuses the established
restriction-inspection result. Older job indices and minimum schemas are retained.

The picker exposes all/zero/one/two/three-or-more counts. Exact-count filters only
admit fully scanned enzymes with no possible IUPAC matches; uncertain/incomplete
counts have their own filter and never imply absence. Counts describe recognition
sites, not experimental cleavage. Outside-substrate and unsupported cuts remain
unavailable rather than receiving fabricated coordinates. Results are paged, and
the downloadable JSON retains the original source-bound scan reports and counts.

Recognition spans and actual top/bottom cut boundaries appear as separate marks.
The map uses purple recognition, blue top cuts and orange bottom cuts; possible
matches are dashed/marked with `?`. Labels are bounded separately from markers.
All returned marker geometry remains available, while text-density omissions are
reported. RE clicks retain original site identity and select every recognition
span, including origin-wrapping sites. Full tooltips identify recognition ranges,
orientation and both cuts. The SVG map action includes the same selected RE layer when a fresh scan is available; its bounded native renderer preserves every marker or returns a size error, and reports label omissions explicitly. The map's document-event context owns its gesture;
new analysis results for that same context do not cancel rotation.

The sequence worker's viewport report remains separate from the shared RE scan.
Its completeness flag describes that report alone. Displayed RE tracks consume
the source scan directly; no second viewport-specific enzyme scan can disagree
with the map or the enzyme picker.

Compiler checks run in the background while implementation continues. Focused native scan/wire and MCP fixtures were added; their execution and browser acceptance are not asserted by this document. The `scan_restriction` MCP tool exposes the same scan-only evidence, catalog pins, paged sites and uncertainty-aware counts without manufacturing digest products. The complete product goal remains active.
