# Assembly wizard: visual workflow reference

Reviewed 2026-11-25. This supplements the semantic review in
[BENCHLING_ASSEMBLY_WORKFLOWS.md](BENCHLING_ASSEMBLY_WORKFLOWS.md). The observations
below come from actual official screenshots/GIF frames, not text alone. No
Collagene application browser, build or test was run. Images remain research
references under `/tmp`, not product assets.

## Ordinary single-construct wizard: a small launcher, then a bottom dock

The official [digest/ligation GIF](https://help.benchling.com/hc/article_attachments/39656605985165)
was inspected at frames0,20,40,60,99,120,199,250,298,350,397 of398.
The current [wizard guide](https://help.benchling.com/hc/en-us/articles/39656605989901-Create-assemblies-with-the-assembly-wizard)
still links this interface.

- Frame20: a centered, roughly850×470 strategy modal over the source viewer.
  Three radio choices—Digest and Ligate, Gibson, Golden Gate—have short
  descriptions; Start/Cancel sit at bottom right. It is not a long parameter form.
- After Start, the modal closes. A bottom dock occupies roughly the last quarter
  of the1080px-high viewer, while source tabs, sequence/map and selection stay
  usable above it.
- Bottom chevrons represent Backbone, Insert and added contributions. Each
  retains source name/length and interval or enzyme information while switching
  documents.
- The active fragment occupies the dock's left half: selection instruction or
  duplex-end preview, Set from Selection/Change, reverse and source-jump actions.
  The right half holds Overall Assembly findings, a name field and Assemble.
  Red/green dividing rules distinguish the illustrated error/ready states.
- Final frames open the created sequence in an ordinary document tab.

These visual states establish the hierarchy, not Benchling's scientific guarantees.
The accompanying instructions use cut-bounded selection, an inverted retained
backbone and explicit orientation correction. Our exact duplex/source-event
semantics remain stronger than a generic visual success tick.

## Gibson: the same persistent dock, with preparation and spacers in context

The official [Gibson GIF](https://help.benchling.com/hc/article_attachments/39656605986061)
was inspected at frames0,35,87,140,210,280,350 of351.

Frame87 shows a whole circular backbone selected at a single BamHI boundary;
the dock states its retained interval and enzyme. Frame140 keeps that backbone
chevron populated while an insert is selected from another source tab. Frame210
shows the prepared fragment's two strand ends, compact DIGEST/PRIMER labels,
Change/reverse/source-jump controls and local status. Frame280 edits a named
spacer with just a bases field and Save in the active-fragment area. Naming and
Assemble remain in the same right-hand position. The final frame opens a small
sequence/primer destination modal with Create/Cancel.

This animation does not show every primer-editing dialog or a separate complete
primer-review table. Do not claim those exact controls were visually verified.
The useful pattern is source selection → one active preparation → local review,
with the rest of the construct still visible.

## Golden Gate junction/primer review: historical evidence, clearly labeled

The [2014 official Golden Gate article](https://www.benchling.com/blog/golden-gate-assembly-wizard)
provides a [setup/error screenshot](https://images.ctfassets.net/kzeezny59h5p/2E9MdmXrBoUB9SzsN2DGAx/096b5f3fa986a13bb7fddacea2a318d0/2015-04-screenshot-2014-12-16-14-41-56.png)
and [review screenshot](https://images.ctfassets.net/kzeezny59h5p/3PlP20jcEJWvKNevhzYgkS/a184cee594d1782540a2f47e20fbd5fe/2015-04-screenshot-2014-12-16-14-55-43.png).
Both were inspected. The first pairs an active fragment's strand-end letters
with an error explicitly naming two identical internal sticky ends and their
fragments. The second pairs a primer table with an upstream/downstream junction
table, sticky-end bases and reverse complements; export/print/finalization actions
are together. This is useful visual precedent, not evidence that today's dialogs
or thermodynamic algorithms are identical. Unsupported numerical metrics must
not be invented to reproduce the table.

## New combinatorial tool: a different document-level workflow

The current [combinatorial guide](https://help.benchling.com/hc/en-us/articles/39922419457165-How-to-use-the-combinatorial-cloning-and-concatenation-assembly-tools)
shows a full assembly document, not the ordinary bottom dock. Actual visuals
inspected:

- [Add-fragment GIF](https://help.benchling.com/hc/article_attachments/39922419445901),
  frames0,27,55,82,109 of110: ordered bin cards above Fragments and Constructs
  tables; per-bin preparation dropdowns; a search modal followed by a fragment
  modal with Start/End/Orientation beside a map. Each bin contains alternatives.
- [Construct-review GIF](https://help.benchling.com/hc/article_attachments/39922491660941),
  frames0,57,114 of115: a construct overview, then a large review modal with
  construct identity at left and Sequence/Primers tabs, view/zoom/gear and Done.
- [Finalization GIF](https://help.benchling.com/hc/article_attachments/41889229351693),
  frames0,60,147,293 of294: two steps, Define outputs then Review and finalize;
  an optional amplified-fragment card; final creation; then a read-only assembly
  record. Enterprise schema/folder/notebook/worklist settings are excluded here.

Keep this visual scale for an explicitly selected library task. Do not make a
simple backbone+insert user populate a matrix or navigate enterprise output forms.

## Proposed Collagene presentation using existing models

This is our design recommendation, not another engine or a claim of current
Benchling behavior.

1. **Launch once.** Assemble opens a compact strategy choice or resumes the
   persistent draft. Keep an immediately reachable current-method selector.
   Library mode remains a deliberate secondary action.
2. **Keep source viewing central.** Present the ordinary draft as a compact
   bottom area with Backbone → Insert → additional slots. Preserve source tabs,
   selected region, exact captured event, orientation and declared interval.
   Switching source documents must not replace the active slot or lose the draft.
3. **Show one active fragment, not every form at once.** Its summary contains
   source/region/direction, preparation and both physical ends where established.
   Set from selection, Choose digest fragment, Reverse and Change are local
   actions. Detailed preparation/primer settings expand for that row only.
4. **Keep intent and findings legible.** A compact intended-contribution ribbon
   and local junction/error area remain visible beside the active row. Name the
   offending row/junction; expose source navigation and the existing correction
   controls. Keep independent intended DNA distinct from actual reconstruction.
5. **Review deliberately.** Open the existing actual-product viewer and primer/
   junction evidence after an explicit Review action. Preserve manual primers;
   show a proposal before Apply. A wrong-target actual product remains useful
   for diagnosis, with its invalid status visible.
6. **Finalize the reviewed recipe.** A compact output sheet names the construct
   and offers optional primers/prepared references, order export and HTML/PDF.
   Reuse existing fresh reconstruction and immutable lineage; hidden prepared
   documents remain retained. No folder/project hierarchy is required.

The newly supplied [22:42:42 Create screenshot](</tmp/Screenshot 2026-11-25 at 22.42.42.png>)
shows Annotation, Primer, Translation, New AA sequence, New DNA, New RNA and New
part, with directional submenus. It is not an assembly-wizard screenshot; the
wizard's documented launcher is the separate Assembly control. This resolves the
previous visual note's uncertainty about Create without conflating the menus.

New `/tmp/Screenshot*` files were checked during this pass; no additional copies
beyond22:42:42 were present at completion. The sampled frames establish visible
states, not every gesture, animation transition or runtime acceptance result.
