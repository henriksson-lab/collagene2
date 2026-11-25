# Assembly workflow review against official Benchling guides

Research date: 2026-11-25. This review uses public official Benchling help pages
and current repository source. It supplements
[CLONING_INTERACTION_DESIGN.md](CLONING_INTERACTION_DESIGN.md),
[BENCHLING_RESTRICTION_GOLDEN_GATE.md](BENCHLING_RESTRICTION_GOLDEN_GATE.md) and the
primer/history reviews. No Benchling code or assets were copied. No local browser,
compiler or test checks were run. Source-present functionality below is not a
claim of completed GUI acceptance. Projects, ELN, inventory, scheduling and new
external primer/alignment engines are outside scope.

## Documented interaction patterns

### Single construct: select sources inside the viewer

The official assembly-wizard guide, updated September 30, 2025, starts from an
open sequence. Digest/ligation selects a retained backbone using cut boundaries
and inverted selection, then an insert in another source. Slots receive the
selection; orientation can be reversed to resolve incompatibility. Users review
a preview, name the construct and create the output. Gibson similarly chooses a
linearized backbone and insert region, previews, resolves errors and creates
sequence/primer outputs. This is a useful compact default for ordinary cloning,
without introducing a combinatorial matrix first.

The guide also states that its Gibson implementation uses sense-strand cut
coordinates rather than overhang directionality. That is a documented modeling
limitation, not a reason to weaken our actual duplex/end reconstruction. Its
FAQ's restrictions are specific to that wizard; they should not be generalized
to every Benchling assembly mode.
[Create assemblies with the assembly wizard](https://help.benchling.com/hc/en-us/articles/39656605989901-Create-assemblies-with-the-assembly-wizard).

### Restriction inspection: connect sites, fragments and source coordinates

The current restriction guide, updated October 2, 2025, documents enzyme search
by name/cut count, selecting a site to navigate to it, and region-based map
filters. Running a digest opens fragment results; clicking a fragment row
highlights it on the sequence. Named digests can be reopened. A predicted-gel
view is also documented. For this product, the useful pattern is direct
site → actual fragment → assembly selection; predicted-versus-observed gel
comparison remains excluded.
[View restriction enzymes](https://help.benchling.com/hc/en-us/articles/38609319462669-View-restriction-enzymes).

### Multiple preparations: preserve row choices and show status in context

The December 12, 2025 combinatorial guide supports ordered fragment bins,
region/orientation preview and preparation choices per bin. Gibson distinguishes
new PCR homology, existing homology and restriction preparation. Tables expose
fragment/construct status and per-junction overlap lengths, with construct review
before output creation. Optional primer/amplified-fragment outputs and assembly
history retain the design connection.

Previously saved oligos can be preferred primers. The guide says unsuitable
preferred primers may be replaced by newly designed primers, accompanied by a
status warning. Our adaptation should preserve manual rows and request explicit
Apply before replacement. The same article's fragment-edit section restricts
orientation changes there to the backbone, despite a broader earlier reference
to region/orientation preview: exact behavior across its dialogs is uncertain.
This does not establish arbitrary per-fragment orientation editing in every
Benchling mode.
[Combinatorial cloning and concatenation](https://help.benchling.com/hc/en-us/articles/39922419457165-How-to-use-the-combinatorial-cloning-and-concatenation-assembly-tools).

### Output portability is separate from creating a construct

Benchling's sequence overview, updated September 30, 2025, documents individual
sequence export from sequence information, including format/options selection
before download. It does not establish that a GenBank file recreates the full
assembly process or experimental evidence.
[DNA and RNA sequence overview](https://help.benchling.com/hc/en-us/articles/39110680274189-DNA-and-RNA-sequence-overview).

## Current repository connections

These are our source observations, not Benchling claims or executed workflows.

| Stage | Current implementation | Remaining question / proposed adaptation |
| --- | --- | --- |
| Choose sources and order | [assembly_workspace.rs](../crates/web/src/assembly_workspace.rs) retains an assembly draft independently of the active source viewer, with named slots, source-event capture, selected/whole/complementary circular regions, ordering and reversal. | Keep the active slot, pinned source and intended orientation immediately visible while navigating source documents. Treat any discoverability issue as an interaction question to review later, not missing scientific scope. |
| Choose digest fragments | [digest.rs](../crates/web/src/digest.rs) connects map/sequence cut picks, opposite circular path, actual fragment reconstruction, source highlighting and Use in active assembly slot. It reports both strand lengths/ends and preserves the catalog selection. | Keep a one-route path from a cut label to the prepared slot. Explain internal cuts, unknown source termini and end repair at that specific slot/boundary. Never substitute a selected top-reference slice for a physical digest. |
| Prepare overlap rows | [cloning_methods.rs](../crates/web/src/cloning_methods.rs) has PCR, digest/catalog digest and existing-material preparation, saved-primer copying, method parameters, explicit suggestions and Apply. It uses original sources and separately frozen intended target. | Make preparation state legible per slot: existing digest, manually reviewed PCR, or selected new design. Show what changes before Apply and retain manual names/domains/tails when their row is not selected. |
| Edit junctions selectively | [assembly_junction_edit.rs](../crates/core/src/assembly_junction_edit.rs), [selective_junction_edit.rs](../crates/core/src/selective_junction_edit.rs) and the assembly dock expose explicit junction trims/inserts, independent edited intent and selected PCR regeneration with variable overlaps. Retained preparations remain exact. | This is now source-present; do not report it as an absent engine. The next useful review is whether a user can find a junction, understand synthetic versus source-derived bases, see affected rows and apply only chosen replacements without losing context. |
| Review actual products | [assembly_product_view.rs](../crates/web/src/assembly_product_view.rs) checks the exact final-step receipt and displays a read-only virtual sequence/map, actual product coordinates and fresh junction inspection. A wrong-target actual product stays inspectable; a failed final step cannot masquerade as its last preparation. | Keep intended reference and actual reconstruction visually distinct. Locate errors at their source row/junction; expose alternatives and incomplete enumeration without turning warnings into a generic success tick. |
| Save and export | [assembly_finalize.rs](../crates/web/src/assembly_finalize.rs), [cloning_outputs.rs](../crates/core/src/cloning_outputs.rs), [plan_oligos.rs](../crates/core/src/plan_oligos.rs) retain source/run history, name the construct, optionally save primers and control prepared-document visibility. Actual-product review offers order CSV/TSV and HTML/PDF steps. | Keep these actions on the same reviewed plan and target. Prepared fragments hidden from Documents must remain retained, and incomplete order data must not silently produce an all-required sheet. Native history and text exports have different preservation guarantees. |

The older central design note's remaining-work wording predates some new junction
and selective-regeneration source. Read it with the current implementation and
[SELECTIVE_JUNCTION_MCP.md](SELECTIVE_JUNCTION_MCP.md), rather than starting a
second disconnected primer/junction tool.

## Cloning-first priorities

**P0: One ordinary two-fragment route.** Present select backbone → select insert
→ prepare → inspect junction/product → create outputs as the primary path.
Additional ordered rows should extend it. Users should not need an imported
full target: region/order/orientation and explicit edits define independent
intent before the preparation simulator runs. An optional separately supplied
target remains useful for checking a pre-existing specification.

**P0: Preserve user choices during edits.** When a junction changes, show changed
source boundaries, literal inserted DNA, new primer domains/tails and the precise
rows proposed for regeneration. Keep retained PCR/digest rows exact. If a retained
row cannot produce the edited contribution, show its failure; do not silently
redesign it. The selective compiler and explicit Apply are already suitable
foundations. Do not infer an experimental or chemical state from oligo names,
annotations, a target hash or imported producer claims.

**P0: Useful local diagnostics and recovery.** Pair each error with an actionable
location: extra source cut, unsupported end geometry, alternative PCR binding,
wrong orientation, nonmatching overlap, invalid independent target or stale source
event. Preserve the rest of the draft while the user corrects that location.
Recompute affected preparation and the complete assembly; no cached positive
status survives changed inputs. No new simulation engine is implied by this UI
recommendation.

**P1: A compact row/junction overview.** Show source, orientation, preparation,
overlap/overhang, primer origin (retained versus newly proposed) and local status
in one reachable overview. Keep sequence/map review central and allow exact
source/junction navigation. Advanced filters and raw JSON are secondary evidence,
not the main way to understand an ordinary cloning failure.

**P1: Output continuity.** Final creation, primer sheet and HTML/PDF steps should
all refer to the same reviewed recipe. Product → producing run → exact historical
source should remain navigable; source names or matching DNA do not establish
ancestry. Read checks can follow a created construct later without delaying this
core assembly flow. Do not add project administration to solve output naming.

## What remains uncertain

- The current guides describe controls and linked animations, not every mouse
  gesture or validation edge case. No frame-by-frame animation or live Benchling
  tenant study was performed for this note.
- Current public documentation does not fully specify how manual primer edits,
  junction edits and selection persistence interact in each wizard generation.
  Historical blog examples are not evidence that the exact same editing UI exists
  today. Our explicit retain/regenerate model is a proposed adaptation.
- A preview/status label in help text is not proof of exhaustive product search,
  exact source-event replay, nick sealing, methylation state or experimental
  success. Our independent intended-reference and fresh reconstruction contract
  should remain explicit.
- No GUI acceptance result is claimed for the current repository. Later focused
  user review should assess the connected cloning route once in-scope code is
  complete, as requested; this research introduces no new test gate.

## Focused follow-up: topology, retained primers and uncertainty (2026-11-25)

Reopened the two official guides above rather than introducing another broad
workflow audit. The current combinatorial article explicitly permits linear or
circular products, per-bin production choices, and an explicitly populated
construct table instead of every combination. It documents primer status that can
differ between constructs. Errors must be resolved before assembly, while warnings
may be accepted. Sequence/map and primer review precede optional primer/fragment
output creation. These are documented interactions; the page does **not** specify
probabilistic fusion-set ranking, outer-end competition, domestication search,
exhaustive alternative-product verification, or exact source-event identity.
[Official combinatorial guide](https://help.benchling.com/hc/en-us/articles/39922419457165-How-to-use-the-combinatorial-cloning-and-concatenation-assembly-tools).

The ordinary wizard's documented recovery is to change insert orientation and
review the preview. Its Gibson FAQ explicitly uses sense-strand cuts rather than
overhang direction. This is a stated limitation of that wizard, not a requirement
to drop strand/end evidence in our implementation.
[Official single-construct guide](https://help.benchling.com/hc/en-us/articles/39656605989901-Create-assemblies-with-the-assembly-wizard).

Current-source implications, not Benchling claims:

- The earlier Golden Gate gap table predates the shared PCR/digest planner,
  whole-pool selection, fusion placement and immutable batch/release work. These
  engines and their GUI connections now exist in source; do not commission a
  replacement implementation from that historical table.
- `crates/web/src/assembly_fusion_view.rs` now carries explicit signed placement,
  retained manual rows, profile choice, and separately scoped physical outer-end
  competitors. Keep mathematical candidate feasibility, empirical conditional
  ranking and actual product reconstruction visibly distinct. The official pages
  establish no basis for presenting ranking as experimental success.
- `crates/web/src/assembly_workspace.rs` already offers explicit selected-row
  regeneration. Retain that deliberate application rather than copying the
  documented preferred-primer replacement behavior without review. Diagnostics
  should name the retained/new row and affected junction, including per-construct
  differences.
- `crates/web/src/assembly_release_view.rs` already reconstructs selected rows
  privately before publication. Keep its optional primers/prepared-document
  choices connected to those same receipts. The remaining user-facing assessment
  is discoverability and continuity, not a missing output compiler.

Unresolved public-documentation questions remain narrow: exact manual-edit
persistence through every dialog, invalid-junction focus/recovery gestures, and
physical competition/chemistry guarantees are unspecified. No additional wizard,
project administration, runtime GUI inspection or new acceptance gate follows
from this review.
