# Benchling editor workflow review

Reviewed 2026-11-25 against official public documentation and the current Collagene source. This is a source assessment, not a runtime acceptance check. No local browser checks, builds, Benchling application code or authenticated sessions were used. Earlier gap documents describe earlier implementations; the assessment below supersedes their editor-specific claims where the code has since changed.

Subsequent implementation: the editor now has a session-local annotated reference
clipboard, explicit paste preview and native source/destination projection with
exact-event parents and Undo/Redo. Plain-text OS clipboard operations remain
separate; rich transfer does not assert digest/ligation chemistry. This supersedes
the original missing-rich-transfer recommendation below. [Current contract](REFERENCE_TRANSFER.md).

## Documented Benchling behavior

**Sequence selection and editing.** The DNA/RNA overview describes drag selection with coordinates, selected length and sequence statistics; sequence, linear and circular views; display toggles; direct base insertion/deletion; and copying regions with their annotations into another location or document. It distinguishes changing a circular sequence's index from moving a map visually. [DNA and RNA sequence overview](https://help.benchling.com/hc/en-us/articles/39110680274189-DNA-and-RNA-sequence-overview), updated September 30, 2025.

**Restriction site discovery.** The restriction guide describes enzyme-name/site-count search, navigating to an individual recognition site, and highlighting selected enzymes on the map. It supports reusable custom enzyme lists and selection-relative filters such as inside/outside the selected region. Its advanced filters affect map display separately from the enzyme list. This distinction should be made explicit in any corresponding Collagene controls. [View restriction enzymes](https://help.benchling.com/hc/en-us/articles/38609319462669-View-restriction-enzymes), updated October 2, 2025.

**Annotation management.** The feature guide documents region-to-annotation creation, coordinates including circular wrapping, strand, type, color and custom fields. Expanded feature views support selecting multiple annotations and changing type, strand or color together. Visibility is distinct from deletion. A feature library provides reusable annotation proposals. These descriptions do not specify how every compound/fuzzy GenBank location behaves during edits. [How to characterize sequences with annotations and translations](https://help.benchling.com/hc/en-us/articles/39768268305805-How-to-characterize-sequences-with-annotations-and-translations), updated September 30, 2025.

**Copy/paste cloning.** A separate guide selects restriction-flanked insert and destination regions, then confirms insertion. It says transferred annotations/translations are retained and the application may reverse fragments and handle sticky ends. This documents an interaction; it does not establish a reproducible physical verification algorithm. Collagene should keep orientation explicit and use its native digest/ligation verifier for chemical claims. [Clone by copy and paste](https://help.benchling.com/hc/en-us/articles/9684297597453-Clone-by-copy-and-paste), updated September 22, 2025.

## What current Collagene code already implements

These are implementation observations, not claims that appearance, performance or every gesture has passed final browser acceptance:

| Workflow | Current source evidence |
|---|---|
| Select and edit bases | `crates/web/src/sequence_view.rs` handles pointer selection, Shift extension, arrows, Select All, typed replacement, Delete/Backspace and Undo/Redo. Carets use boundaries; complex or unresolved selections are protected from accidental contiguous replacement. |
| Contextual selection actions | `selection_toolbar.rs` reports coordinates/length/GC and exposes annotation, forward/reverse primer, mutation and base-edit actions. It preserves selected feature parts rather than substituting a bounding interval. |
| Circular overview | `app.rs` wires select-range versus rotate modes, zoom, annotation and restriction hit navigation into the map. `restriction_map.rs` draws recognition spans and separate strand-cut evidence, including circular coordinates. Visual rotation is separate from document reindexing. |
| Restriction count filters | `restriction_view.rs` shares a source-bound scan between map, sequence and enzyme picker. Filters include 0, 1, 2, 3+ and uncertain/incomplete sites; incomplete or ambiguous scans do not establish an exact count. Its name/recognition search and count filter affect both picker and displayed overlay. |
| Physical fragment selection | `digest.rs` supports fresh digest fragment navigation, paired restriction-boundary picking and assignment of a reconstructed fragment to assembly. `assembly_workspace.rs` supports keeping the outside of a circular selection. These are stronger existing connections than the older research document's inspector-only assessment. |
| Clipboard | `clipboard.rs` handles Copy, reverse-complement Copy, Cut, replacement Paste and insertion Paste, with revision/selection guards across asynchronous browser permissions. It transfers text; paste calls the source's replacement operation. An additional internal annotated-reference clipboard now projects source features into the destination with an explicit preview; plain text remains text. Clipboard API access requires HTTPS or localhost; the existing edit dialogs remain a text-entry route on the HTTP preview. |
| Annotation/document editing | `editors.rs` edits individual feature kind, exact/compound location and qualifiers and marks changed locations for review. Document controls cover topology, reverse complement and origin changes. `feature_tools.rs` supplies user-owned feature libraries and reviewable automatic annotation proposals. |

## Implementation follow-up, 25 November 2026

The first three recommendations below have now been implemented. Annotated internal
reference transfer retains source-event ancestry and previews annotation projection
([details](REFERENCE_TRANSFER.md)). Selection-relative restriction discovery pins
exact range unions and separates recognition from strand-cut geometry; named enzyme
sets retain catalog identity and refuse changed definitions
([details](RESTRICTION_SELECTION.md)). These recommendations remain as research
history, not the current missing-feature queue. Annotation list operations and
editor-to-design continuity remain candidates for follow-up.

## Prioritized remaining connections

These are Collagene recommendations inferred from the research and source comparison, not assertions about Benchling internals.

1. **Annotation-preserving reference transfer.** Add an explicit, bounded internal sequence-fragment payload alongside ordinary plain text. Capture original source/event, ordered selected parts and intersecting annotations; preview what survives, splits or needs review; transform coordinates and strand on reverse-complement insertion; commit the destination once with Undo and provenance. Keep untrusted external plain text as sequence text. Do not claim ligation or silently choose orientation from clipboard content. The existing digest-to-assembly route remains the physical cloning path.
2. **Selection-relative enzyme discovery.** Extend the shared restriction scan with clearly labeled inside/outside selection filtering. Define whether the rule tests recognition coverage or actual cleavage boundaries, especially for Type IIS enzymes and origin-spanning selections; expose the interpretation instead of conflating the two. Preserve possible/incomplete results separately. Offer the same explicit choices to map and picker rather than reproducing Benchling's surprising filter mismatch.
3. **Reusable enzyme choices.** Add named user-owned enzyme sets with create/update/select actions and native-workspace persistence. Preserve catalog identity and missing-enzyme diagnostics. This is a local molecular-biology preference; organization/project administration is unnecessary.
4. **Annotation list operations.** Add multi-select and an expanded list/table, then bulk type/strand/color edits and show/hide actions. Distinguish display preferences from edits to the scientific record. Keep existing qualifiers and unsupported location detail intact; use one guarded history commit for a batch. The current panel primarily selects and edits one annotation at a time; the current feature editor has generic qualifier rows rather than dedicated color/strand controls.
5. **Editor-to-design continuity.** Keep selection, source identity and the receiving assembly slot visible as users move among sequence, map, restriction and primer actions. Existing selection-to-slot and digest-to-slot code should be reused. A contextual shortcut can expose those actions without introducing another coordinate or chemistry engine.

The next editor work should complete these connections rather than add a second map, keyboard editor or restriction scanner. Browser acceptance and smoothness measurements remain deferred until implementation is complete, as requested. No ELN, inventory, organizational workflow or mobile requirements are introduced by this review.
