# Viewer visual reference and drawing decisions

Reviewed 2026-11-25. Actual images were inspected: the user's two22:34 screenshots,
all nine22:40 menu screenshots, both22:42 toolbar-menu screenshots, five additional
official Benchling PNGs, and the
first frame of its primer-visibility GIF. Public image search helped locate
references; third-party tutorial images were excluded. No Collagene browser or
build was run. This records visual observations and our implementation choices,
not undocumented Benchling algorithms or pixel-perfect copying.

## Restriction cuts and sequence tracks

The user's [Benchling screenshot](</tmp/Screenshot 2026-11-25 at 22.34.48 benchling.png>)
shows AscI with a bold label above the sequence, translucent recognition-span
highlight across both aligned strands, and a dark staircase drawn **at the actual
interbase boundaries in the letters**. The hovered map label is also darker.
The [Collagene screenshot](</tmp/Screenshot 2026-11-25 at 22.34.21 collagene.png>)
places a small cut diagram below the sequence and separates the cut description
from the bases, making the relationship harder to read.

Implement the in-sequence geometry using the existing exact top/bottom cut
coordinates: one vertical stroke at each strand's interbase boundary and a
horizontal joining stroke in the gap between strands. Keep recognition shading
separate from cuts, including Type IIS cuts outside the motif. A hover/focus state
should emphasize the same exact site in sequence and map. Preserve distinct
unknown cuts rather than drawing guessed geometry. These are our drawing rules.

The official [enzyme-detail image](https://help.benchling.com/hc/article_attachments/38609319456141)
adds another useful treatment: a two-strand motif with opposed cut markers in a
small enzyme card beside the digest panel. It supports a compact hover/focus
explanation rather than permanent explanatory paragraphs. The
[site-label image](https://help.benchling.com/hc/article_attachments/38609323397645)
shows short colored enzyme-name chips, fine leaders and ordinary monochrome DNA.

The official [sequence/linear-map screenshot](https://help.benchling.com/hc/article_attachments/10465904427533)
shows filled, labeled feature arrows below the DNA, separate translation cells,
and multiple compact restriction-label lanes above it. The user's LTR spans
several sequence rows: intermediate pieces have continuation marks; the actual
feature end has the directional point. Use real source endpoints for arrowheads,
not a fresh biological arrowhead at every wrapped row. Keep direction, feature
color and readable labels; ordinary features should not look like thin underline
rules alone. At default zoom, reduce our oversized track headers and repeated
help text before shrinking readable bases.

## Circular and linear map organization

The user's circular map places feature arrows just inside the backbone; longer
names follow their arcs. Restriction names sit outside, with short radial/elbow
leaders and collision-aware fan-out. Nearby names sometimes share a text group;
dense regions show `+N`. Hovered AscI is promoted visually. Our screenshot instead
uses long vertical label columns and wrapped external feature names.

Implement measured label boxes with angularly ordered placement, compact leaders
and separate feature/enzyme label zones. Keep the chosen/hovered label readable.
Grouping or `+N` must provide the complete site list on inspection, with exact
positions retained; it must not remove scientific hits. Feature labels that do
not fit should use a short label or inspected detail, never spill wrapped text
through neighboring tracks. The user's22:40:47 linear map similarly keeps dense
RE labels above the axis and filled directional feature lanes below it, with
shortened labels inside the arrows. Its large central cut cluster demonstrates
why label density should not be confused with missing sites.

## The right rail is an anchored tool launcher

All nine new user images show the same narrow vertical icon rail. The selected
icon gains a circular backplate and opens **one floating card immediately to its
left**. The underlying pane geometry stays unchanged. Cards are about350px wide
in the2219px screenshots and fit their contents; they are not permanently open
full-height columns. These are direct screenshot observations:

| Local screenshot22:40:ss | Visible card/tab structure |
| --- | --- |
|29|Sequence Info / Tags; compact labeled fields, topology dropdown, update action, separated export controls.|
|35|CRISPR heading; primary design action, saved-analysis empty state, separate template action.|
|38|Sequence Alignments; create action and saved-list/empty state.|
|40|Version History; compact change/user/time rows, reset action, then a separately labeled clone-version form.|
|41|Primers / Pairs; empty state, Create primers dropdown and Attach existing.|
|43|New Digest / Saved Digests; enzyme-list and map-visibility dropdowns, search, Name/Cuts and Selected/Color lists, region dropdown, compatibility checkbox, Run digest.|
|44|Annotations / Translations; visibility dropdown, Create new, total/hidden counts, bounded list of chevron + color + name/range + eye, expanded-view footer.|
|47|Linear Map is a main-pane tab; it is not another rail card.|
|52|Metadata is a main-pane tab; its enterprise registration controls are not part of our requested cloning interface.|

The same compact pattern appears in the official
[digest rail/card screenshot](https://help.benchling.com/hc/article_attachments/38609319453325),
[map-filter dropdown](https://help.benchling.com/hc/article_attachments/38609323397901)
and [open region dropdown](https://help.benchling.com/hc/article_attachments/38609319456141).
The [restriction guide](https://help.benchling.com/hc/en-us/articles/38609319462669-View-restriction-enzymes)
confirms that list search and map visibility are different controls. Preserve
that distinction when compressing our existing controls.

Our implementation choices: a roughly40–48px rail; one closable340–380px card,
clamped to available viewport height; compact tabs and grouped controls; long
lists scroll inside the card; large editors/results open the existing pane.
Retain tool state when closing a card. Use accessible icon names/tooltips,
`aria-expanded`, keyboard focus, Escape and click-outside dismissal. Icons should
have consistent meaning; do not recreate nine text tabs beside the rail.

## Top toolbar dropdowns: two additional supplied images

The user's [22:42:23 screenshot](</tmp/Screenshot 2026-11-25 at 22.42.23.png>)
shows **Analyze** open: Run Primer3, Submit to NCBI BLAST, Analyze as translation,
and Optimize codons → Forward / Reverse. This is direct evidence for a short
primary menu with a directional submenu. It is not a request to add the deferred
Primer3 engine or network analysis. Create is adjacent but its menu is not open
in these images; do not invent its exact items.

The [22:42:34 screenshot](</tmp/Screenshot 2026-11-25 at 22.42.34.png>) shows the
sequence gear menu, about220px wide rather than the350px tool cards. Every item
has its label on the left and one switch on the right:

- Annotations with count8; Primers with count0; Primer Bases.
- Cut Sites; Complement; ORFs with an inline secondary gear.
- Axis; Translations with count2 and an inline secondary gear.
- Amino Acid Indices; Show Cuts.

Counts are small badges, secondary settings stay beside their owning label,
and hover emphasizes a single row. Cut Sites is on while Show Cuts is off:
label visibility and visible cleavage geometry are distinct options. Our
implementation should retain that distinction, and can reveal precise geometry
on site hover/focus even when persistent cut glyphs are hidden. This latter
hover policy is our design choice, not a timing guarantee derived from stills.

Keep Create/Analyze and view-display controls in compact dropdowns. The official
[sequence overview](https://help.benchling.com/hc/en-us/articles/39110680274189-DNA-and-RNA-sequence-overview)
places annotation, primer, complement, axis, translation and cut display under
a gear menu. The [primer-visibility GIF](https://help.benchling.com/hc/article_attachments/39072049400717)
was inspected only at its first frame, so this note does not infer its full menu
animation or timing.

Research images remain in `/tmp/collagene2-benchling-visual-ref`; they are reference
material, not bundled product assets. The user's images are the most relevant
layout target. Preserve our data, coordinate and uncertainty semantics while
adopting their compact visual hierarchy.
