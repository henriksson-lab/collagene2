# Viewer drawing and command hierarchy

The 2026-11-25 user screenshots established concrete deficiencies: restriction
labels were detached from their sites, feature tracks lacked clear direction,
long map names broke inside words, and exposed actions consumed the editing area.
The supplied Benchling hover/menu screenshots and additional official references
are recorded in [the visual comparison](BENCHLING_VIEWER_VISUAL_REFERENCE.md).

## Drawing contracts

Sequence restriction labels sit above the aligned strands with bounded staggered
lanes and short leaders. Hover highlights recognition bases on both strands and
shows the actual staggered interbase cleavage path in the letters. Shared hover
uses the same source-bound site index in sequence and map, without changing the
DNA selection or opening a tool panel. Possible/unknown cut evidence is not turned
into a determined cut. Label density has an explicit overflow affordance.

Features use filled directional bodies. Only the actual biological terminal gets
an arrowhead; row clipping, internal joins and circular-origin splits do not
create extra heads. Complemented joined features retain biological order.

Circular map RE names follow angular order outside the backbone with outward
elbow leaders. Co-located names form compact groups. Labels and leaders avoid
the backbone disk; focused groups are bold with reserved bold text bounds.
Feature bands sit inside the backbone. Nonoverlapping features share a lane;
compound/origin-split intervals use one stable lane and actual overlaps use extra
lanes. Excess overlaps are explicitly counted, preserving source annotations.
The same cached lanes and label boxes drive rendering and hit testing. Feature
text remains readable with deliberate abbreviation rather than arbitrary word
breaks.

Dragging a circular map rotates the view by default. Precise base ranges are
selected in Sequence; the overview's approximate range-mode buttons are removed.
Feature, site and primer-overlay clicks still navigate exact known coordinates.
Empty map-space clicks do not fabricate single-base selections.

## Fixed command surface

A narrow right-hand icon rail opens one floating tool card. Opening a tool does
not resize the sequence or move the map camera. Closing it restores the unobscured
canvas. Advanced tools are grouped in a menu rather than added as permanent tabs.

File, Create and Analyze menus group document commands. Selection actions,
including forward/reverse primer creation, are contextual dropdown entries.
Sequence view toggles, routine statuses, instructions and detailed track evidence
live in a View popover; actual failures and omitted annotations remain visible.
The restriction picker keeps enzyme search, cut counts and the enzyme list near
the top, with region rules and advanced explanations collapsed. These are display
changes over the existing scientific operations and history, not new algorithms.

## Assembly reference and validation

The [assembly wizard visual review](BENCHLING_ASSEMBLY_WIZARD_VISUAL.md) identifies
the ordinary method-picker → persistent bottom-dock workflow, distinct from the
full combinatorial assembly document. Further assembly presentation should reuse
the existing source, preparation, independent-intent, review and output engines,
showing one active fragment rather than exposing every preparation form at once.

The ordinary dock now has an ordered fragment strip with source name, retained
length and direction. Only the active fragment editor is mounted; switching slots
changes view state without replacing pinned sources or scientific evidence.
Fragment capture sits beside that editor, while name/method/topology, reordering
and library alternatives use expandable controls. Intended junction edits remain
available separately; explicit intended review opens the existing preparation
and verification controls. The initial dock height is reduced to preserve source
viewing. A new draft starts with a compact method choice; a draft with captured sources
resumes directly in its active editor. The complete compact product/finalization
presentation remains future work.

Twenty focused geometry cases plus one doctest pass: existing label placement,
radial leaders/boxes, terminal heads, clipping, origin/cut footprints and actual
overlap lanes. The full web/worker build validates the application integration.
No application browser/HTTP smoke checks were performed, as instructed. These
checks establish numerical and build contracts; visual acceptance remains with
the user's review of the published viewer.

## Preparation follows the active fragment

The ordinary draft passes a presentation-only active-row memo and row-navigation
callback to the existing method panel. Only that preparation editor is mounted.
All fragments still participate in the same reconstruction, primer proposals and
verification; changing the active row does not change scientific settings or
clear retained evidence. The standalone method panel retains its multiple-row
layout. Duplicate workflow and disabled row-management controls are omitted in
the fixed-source dock.

Diagnostic Edit preparation links first select the affected draft row, then focus
the newly mounted editor on the next animation frame. Deferred focus checks the
workspace epoch, preparation settings and current selected row. Intended-source
navigation also selects the matching row. Source capture collapses after explicit
intended review, while the intended sequence/junction inspection remains available
in a disclosure and warnings stay outside it. Primer setting grids expand beside
the existing Suggest → inspect → explicit Apply actions.

## Product review controls

The actual-product surface now groups its reconstruction status and output actions
in the header. Primer/HTML/PDF export controls expand together; verification
findings, scientific engine version and retained source/preparation identities
remain in their existing evidence disclosure. Invalid/incomplete and stale-context
messages remain visible.

A compact Sequence/Map/Split toolbar groups map zoom, exact coordinate navigation
and annotation browsing. Split view puts precise sequence at left and the map at
right. The product map uses the same rotation-first interaction as the document
map; no approximate map range or empty-space base selection is offered. Coordinate
navigation from map-only mode reveals the sequence. Actual map omission/lane
statistics link users to the complete retained annotation browser. All drawing is
over the existing actual reconstructed reference; this presentation adds no
scientific verdict or new transformation.
