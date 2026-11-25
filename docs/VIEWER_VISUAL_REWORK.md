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
