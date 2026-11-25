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

## Working drafts and primer proposal review

The Draft menu can start a new assembly or switch to the previous working draft.
Each action first shows the current draft name and number of captured fragments;
confirmation rejects a changed draft, changed previous draft, changed workspace
epoch or a busy workspace. Replacement cancels current intended-review work,
drops method-review owners and their workers, clears junction undo and positive
review state, then installs only the editable declaration. A restored draft must
review intent and reconstruct evidence again. Finalized constructs and recorded
plans are unaffected.

The previous declaration is held at application scope so hiding/reopening the dock
does not lose it. This is deliberately one temporary previous draft, not a saved
design collection: workspace replacement/updates and session termination clear it.
Only the current working draft is serialized by the existing workspace format.
The UI states this limit; users retain permanent designs through existing saved
workspaces, recorded plans and finalized outputs. Persistent multiple editable
drafts remain separate delivery work.

Overlap and Golden Gate primer proposals now use compact accessible tables.
Fragment navigation uses the same exact-snapshot guards as the proposal and opens
the active preparation editor without applying anything. Per-row disclosures
retain both full ordered own-5′→3′ oligos, their separate annealing/tail domains,
and Type IIS tail segments where applicable. Existing calculated Tm/ΔTm/product
length or actual cleaved ends remain visible in the summary; no new numerical
metric is inferred. Warnings stay outside collapsed details. A row with no pair
says so, and an over-bound report shows an explicit error rather than a partial
primer list. Existing explicit Apply and evidence download behavior is retained.

## Physical fragment ends

A pure native/WASM layout module now projects typed physical ends into small
strand diagrams. It places the protruding strand using side and5′/3′ polarity,
keeps terminal phosphate/hydroxyl/blocked/unknown groups on their actual side, and
reverses bottom-strand display order without complementing the reported letters.
No paired letters or end compatibility are invented. Interior geometry and nicks
are explicitly outside the schematic; complete physical lengths, reference
offset, end sequences and nick positions remain in a disclosure.

The selected ligation fragment shows its freshly reconstructed ends in assembly
orientation, after any declared repair/treatment. Existing fragment and junction
inspection retain the same diagrams, with upstream-right and downstream-left
ends facing each other. Type IIS primer proposals show their actual predicted
cleaved ends inside the proposed-chain review. Missing physical state remains
explicitly unavailable; a circular summary has no external linear end pair.
These are views of existing receipts/proposals, not another molecular engine.

End input shape, IUPAC letters and sequence-length limits are checked before
display. Long protrusions show10 leading and10 trailing letters separated by an
explicit elision; full length remains labeled and the complete end sequence stays
in metadata. At most24 glyphs are allocated. Diagrams scroll rather than shrink
long end letters into unreadable text. Four native fixture groups cover all four
nonpalindromic side/polarity cases, terminal groups, malformed input and bounded
long-end display. Visual acceptance remains deferred to user review.

## Retained documents and working tabs

The explorer uses the retained document catalogue; the tab strip uses a separate
ordered, unique set of open document IDs. Opening a retained document selects or
appends its tab. Source navigation and newly created products open the requested
working tab through the existing active-document routes. Preparation-history
visibility affects the explorer only, not the identity or presence of open tabs.

Closing a background tab keeps the current editor and focus. Closing the active
tab selects its right neighbor, falling back to its left neighbor. Closing the
last tab shows an empty view with file/URL import and access to Documents. It
hides document mutation/export controls, clears transient selection and closes
the assembly dock while keeping its draft and all retained documents/history.
Map interaction pauses while no viewer is open. Reopening the same last document
also reinstates its viewer before deferred range navigation applies.

Close buttons are siblings of tab buttons, not nested buttons. Delete and middle
click can close tabs; arrow/Home/End navigation remains available. Active closure
focuses the resulting tab or empty view on the next frame, guarded against a
changed tab list, selection or owner. Background closure does not steal focus.

The archive's retained working-document index remains valid when no tab is open.
Open/closed tabs are presentation state and are not yet serialized: a new session
opens the recorded working document. Frame updates prune absent IDs and keep
matching working source tabs when finalization adds a product. Tab closure deletes no document, source
event, primer, preparation or history graph and changes no scientific/archive
schema. The four pure-model fixture groups cover opening/order,
background/active/last closure and reopening, reconciliation without catalogue
mutation, and reset between workspaces. Browser acceptance remains deferred.

## Clipboard belongs to the displayed sequence

The old window-level shortcut listener matched any `sequence-canvas` element,
then used the app's active source molecule. An inline read-only product canvas
could therefore route a copy/cut/paste shortcut to a different source. That
listener is removed. Trusted native copy/cut/paste events are handled by each
SequenceView using its own molecule, feature, range and read-only/locked state.
Readonly products can copy their displayed DNA; cut/paste cannot edit them or the
active source. Ordinary form fields, dialogs and selectable DOM text keep browser
clipboard behavior. Choosing bases on a canvas clears an older DOM text selection.

The implementation follows the [W3C Clipboard API/events draft](https://www.w3.org/TR/clipboard-apis/#clipboard-event-interfaces):
clipboard events reach focused contexts even when they are not native text
editors, with event clipboardData for explicit handling. Canvas shortcuts use
those user-agent events rather than the asynchronous navigator clipboard API.
Existing menu buttons retain asynchronous access and its secure-context/browser
permission requirements. Browser policy and platform acceptance remain runtime
constraints; no browser checks were run.

Copy preserves an annotation's biological complement/join order. Other copies
require one matching reference span; multiple highlights cannot silently copy
only their first part. Cut/paste require explicit reference coordinates, reject
read-only/busy/ambiguous/stale selections and pass the displayed revision through
the host's existing edit/history callback. Cut stages text in the clipboard event
before requesting an edit; host failure leaves sequence unchanged. Undo retains
the removed sequence in history. Paste validates the input byte limit, normalizes
using the displayed DNA/RNA/protein alphabet and positions the caret using the
normalized length. Menu paste/insert use these same guards, including selected
parts in their asynchronous freshness context. Annotated-reference paste retains
its separate explicit preview/apply route. Four pure fixture groups cover scoped
source/product copy, biological joins, rejected edits and normalized edits/bounds.
