# Restriction/ligation and Golden Gate interaction research

Reviewed 2026-11-25 using official public Benchling documentation and its publicly served screenshots/GIFs. No Benchling application code, authenticated session, local application browser checks, tests or builds were used. This document separates documented behavior, visual observations and proposed Collagene connections. It does not claim feature parity.

## Documented workflows

### Restriction digest and ligation

The wizard asks for a method, then a backbone and insert. Its instructions select between two restriction sites, invert that selection to retain the backbone, and commit it to the backbone slot. The insert uses the selected interval between its sites. An orientation action addresses incompatible ends; a preview is reviewed before naming and creating the construct. The text uses Shift-click, while the older teaching guide describes Ctrl/Cmd+Shift; the gesture is therefore not a portable contract to copy verbatim. [Assembly wizard guide](https://help.benchling.com/hc/en-us/articles/39656605989901-Create-assemblies-with-the-assembly-wizard), [teaching example](https://help.benchling.com/hc/en-us/articles/9684255457805-Molecular-Cloning-Methods).

The digest interface supports enzyme-name/cut-count search, navigation from a site to sequence, and map highlighting. Digest fragment-table rows select their corresponding sequence regions. The documented map filters and enzyme-list filters have distinct behavior. [Restriction enzyme guide](https://help.benchling.com/hc/en-us/articles/38609319462669-View-restriction-enzymes).

A separate copy/paste route selects cut-flanked insert DNA and replaces a destination interval. That guide describes automatic reversal and sticky-end handling. We should not translate a sequence paste into chemical proof or silently flip an explicit user orientation. [Copy/paste cloning guide](https://help.benchling.com/hc/en-us/articles/9684297597453-Clone-by-copy-and-paste).

### Golden Gate preparation and review

The current combinatorial guide distinguishes existing Type IIS sites from PCR introduction of sites/overhangs. Settings include pre-recognition bases/length and annealing constraints. A fragment is assigned to an ordered bin, with region and orientation review. Saved preferred primers are checked for binding, required sites and configured constraints; Benchling documents fallback designs with status warnings. Constructs and primer review precede finalization. Construct, primer and amplified-fragment outputs are separate choices. [Combinatorial assembly guide](https://help.benchling.com/hc/en-us/articles/39922419457165-How-to-use-the-combinatorial-cloning-and-concatenation-assembly-tools).

The teaching example explicitly uses PCR-added directional BsaI sites when the original sources lack suitable sites. This is a preparation decision, not merely searching the starting sequences for recognition motifs. [Golden Gate worked example](https://help.benchling.com/hc/en-us/articles/9684255457805-Molecular-Cloning-Methods).

An older official walkthrough also distinguishes cut-site selection from arbitrary sequence selection requiring new primer tails. It documents spacers, review/editing and primer export. These are historical interaction examples, not evidence of a current API or current scientific guarantees. [2014 Golden Gate walkthrough](https://www.benchling.com/blog/golden-gate-assembly-wizard).

## What the official images actually show

I inspected sampled frames from the [digest GIF](https://help.benchling.com/hc/article_attachments/39656605985165) and [fragment-addition GIF](https://help.benchling.com/hc/article_attachments/39922419445901). Source tabs and the main map remain visible while a lower assembly area retains fragment choices. The fragment modal combines region controls with a map. This supports a persistent design workspace rather than repeatedly opening unrelated inspector forms.

The historical [Golden Gate setup screenshot](https://images.ctfassets.net/kzeezny59h5p/2E9MdmXrBoUB9SzsN2DGAx/096b5f3fa986a13bb7fddacea2a318d0/2015-04-screenshot-2014-12-16-14-41-56.png) shows strand-end letters beside a fragment preview and an explicit duplicate internal-overhang error naming both fragments. Its [review screenshot](https://images.ctfassets.net/kzeezny59h5p/3PlP20jcEJWvKNevhzYgkS/a184cee594d1782540a2f47e20fbd5fe/2015-04-screenshot-2014-12-16-14-55-43.png) pairs a primer table with a junction table listing upstream/downstream parts, sticky-end bases and reverse complements. Those metrics must not be reproduced with invented values.

Static research artifacts were kept under `/tmp/collagene2-benchling-research`, not added as product assets. Digest frames inspected: 0, 99, 199, 298, 397 of 398; fragment GIF: 0, 27, 55, 82, 109 of 110. Sampling confirms visible states, not every timing/gesture transition.

## Current repository assessment

This is a read-only source assessment, not a runtime claim:

- `cloning_methods.rs` now has ordered source rows, explicit digest/PCR/declared-duplex preparation, end treatment, orientation, selection-to-interval copying and proposal/application of overlap primers. It remains an inspector form; it has no persistent central backbone/insert-slot workspace or paired restriction-boundary selection contract.
- `digest.rs` and restriction map overlays already expose recognition/cut evidence and fragment navigation. They are not yet a single action that commits the exact retained digest fragment into an assembly slot.
- `golden_gate.rs` selects complete source templates, enumerates a bounded Type IIS pool and opens a candidate. It is not yet a region-based PCR-tail/site preparation wizard. Enumeration and target verification must remain distinct from designing the substrates.
- Shared recipe compilers, source-event guards, plan verification, junction review and native lineage are useful foundations. A new interface should call them rather than introduce a second coordinate/chemistry implementation.

## Required connections for this project

The following are our implementation recommendations, informed by the research rather than claims about Benchling internals:

| Connection | Concrete behavior and review boundary |
|---|---|
| Persistent assembly draft | Method, ordered slots and source-event references survive source-tab switching. Start with a backbone and insert, allow extra slots; no project/folder hierarchy required. |
| Set backbone/insert | A named active slot receives either a source interval or an exact digest fragment. Show source name, revision, length, retained direction and preparation. Committing a slot never edits the source. |
| Inverted circular selection | Provide an explicit kept-versus-removed choice. Preserve an origin-crossing range as ordered split parts/unrolled coordinates; never approximate it by the smallest bounding interval or confuse two strand cuts with two unrelated recognition spans. |
| Restriction boundaries | Select actual source-bound cut instances, including strand geometry. Commit the physical fragment through the existing digest reconstruction, not by slicing between label positions. Show both retained ends before accepting the slot. |
| Orientation and end errors | Display upstream-right against downstream-left ends in the same declared 5′→3′ convention. Reversal is explicit and triggers fresh checks. Unknown/blocked termini, mismatches and retained nicks remain separate reasons; no automatic end repair. |
| Golden Gate production choice | Each slot chooses existing sites or PCR-added sites. The latter needs source region, enzyme geometry, clamp, recognition orientation, spacer/cut offset and intended fusion overhang. Show added bases separately from annealing domains. |
| Tail correctness | Fresh PCR followed by the selected Type IIS digest must expose the requested ends. Flag internal sites, duplicate/reverse-complement or palindromic end ambiguities and alternatives within the actual bounded model. Do not assume every enzyme creates four-base ends; the catalog includes different geometries. |
| Preferred/suggested primers | Preserve saved names and domains. Show ordered oligo, tail, binding interval, calculated metrics and failed constraints. Proposed replacement requires explicit application. Unsupported thermodynamic metrics stay unavailable. |
| Central product review | Keep sources, assembly preview and junction/primer review connected. Selecting a junction navigates to its actual product bases and contributing source ranges. Distinguish missing substrate preparation, incompatible assembly and independent-target mismatch. |
| Final outputs | Choose construct, designed primers and prepared references separately; name them before commit. Reconstruct and atomically retain actual original-source parents. HTML/PDF steps and machine-readable evidence should describe the same recipe. A design completion action does not claim laboratory success. |

First delivery should connect source selection, two-ended digest fragment commitment, visible junctions and explicit primer review in one draft. Combinatorial libraries, automatic preferred-primer replacement, ELN, inventory, worklists, organizational registration and scheduling are not prerequisites. The existing algorithmic uncertainty/alternative-product evidence should remain visible even where reference screenshots use a simple green check.
