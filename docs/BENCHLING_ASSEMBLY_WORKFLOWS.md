# Assembly workflows: current source-linked decisions

Reviewed 2026-11-25 against official public Benchling documentation and current
repository source. This is a focused clarification of
[BENCHLING_ASSEMBLY_WORKFLOW_CHECK.md](BENCHLING_ASSEMBLY_WORKFLOW_CHECK.md),
not a new feature inventory. No application browser checks, builds or authenticated
Benchling session were used. Earlier gap tables are historical; source-present
features below are not claims of GUI acceptance.

## What the current documentation establishes

| Area | Documented interaction | Implication for our interface |
| --- | --- | --- |
| Ordinary restriction/ligation | Start from an open sequence, set the retained backbone and insert from cut-bounded selections, correct orientation if needed, inspect the preview, then name/create the construct. Golden Gate follows that pattern with Type IIS enzymes. [Assembly wizard](https://help.benchling.com/hc/en-us/articles/39656605989901-Create-assemblies-with-the-assembly-wizard). | Preserve the existing two-slot source-view route as the default. An alternative library is an extension, not a prerequisite. |
| Gibson limitations | The single-construct wizard uses sense-strand cuts and does not expose restriction-site retention as a choice. These are stated limitations of that wizard. [Wizard FAQ](https://help.benchling.com/hc/en-us/articles/39656605989901-Create-assemblies-with-the-assembly-wizard). | Keep our explicit strand/end reconstruction and source-bound junction edits; do not reproduce those limitations. |
| Preparation and order | The newer tool distinguishes PCR-added homology, existing overlaps and restriction preparation. Golden Gate distinguishes existing sites from PCR-introduced sites. Bins are ordered; alternatives within one bin occupy different constructs. [Combinatorial guide, preparation/bins](https://help.benchling.com/hc/en-us/articles/39922419457165-How-to-use-the-combinatorial-cloning-and-concatenation-assembly-tools). | Keep ordered contributions separate from alternatives and from a physical reaction pool. |
| Primers and review | Preferred primers are checked and may be replaced with warnings. Construct review includes sequence/map and primer views. Finalization offers optional primer/fragment outputs and retains assembly history. [Combinatorial guide, primers/review/finalization](https://help.benchling.com/hc/en-us/articles/39922419457165-How-to-use-the-combinatorial-cloning-and-concatenation-assembly-tools). | Our deliberate adaptation is proposed changes followed by explicit Apply, preserving manual rows. |
| HiFi and homology | HiFi is described under Homology; that mode also mentions Gateway with existing ATT regions and In-Fusion with already amplified homology-bearing fragments. [Combinatorial guide, assembly methods](https://help.benchling.com/hc/en-us/articles/39922419457165-How-to-use-the-combinatorial-cloning-and-concatenation-assembly-tools). | These names do not establish an enzyme-specific physical model or justify declaring our Gateway scope complete. |

The worked tutorial specifically distinguishes new Gibson homology arms from
pre-existing arms, and shows PCR introduction of BsaI sites when sources lack
them. It asks users to inspect generated primers before finalization. Its
four-base overhang wording belongs to that BsaI example; our pinned enzyme
geometry remains authoritative, including SapI's different geometry.
[Official cloning tutorial](https://help.benchling.com/hc/en-us/articles/9684255457805-Molecular-Cloning-Methods).

## Concrete corrections to the implementation plan

These are our source observations and design decisions, not assertions about
Benchling internals.

1. **Do not reopen completed preparation engines as missing work.**
   [assembly_workspace.rs](../crates/web/src/assembly_workspace.rs) already retains
   ordered exact-source slots and exposes method/per-slot preparation summaries
   before intended-reference review.
   [cloning_methods.rs](../crates/web/src/cloning_methods.rs) already connects
   existing material, PCR, pinned digest, Gibson/HiFi/In-Fusion profiles and
   Type IIS primer proposals with explicit Apply. Whole-source/whole-amplicon
   Golden Gate exploration also exists. The early Golden Gate gap table predates
   these connections.
2. **Keep preparation recovery local and preserve declarations.** Current
   method diagnostics already link upstream/downstream preparation and source
   rows. Continue that route when improving the interface; do not introduce
   another wizard. A changed source, orientation, primer or method must invalidate
   computed evidence while retaining unrelated draft settings. The official
   guides do not specify every invalid-junction recovery gesture.
3. **Keep intended DNA, reconstructed DNA and laboratory validation distinct.**
   [assembly_product_view.rs](../crates/web/src/assembly_product_view.rs) already
   exposes actual-product and fresh junction review, including wrong-target
   diagnostic products. The Benchling preview/check instructions establish a
   review sequence; they do not establish exhaustive mixture enumeration,
   positional chemistry, independent-target verification or experimental success.
   Preserve our stricter contracts in [PROJECT_PLAN.md](PROJECT_PLAN.md).
4. **Use the same reviewed recipe for outputs and ancestry.** Existing
   [assembly_finalize.rs](../crates/web/src/assembly_finalize.rs),
   [assembly_release_view.rs](../crates/web/src/assembly_release_view.rs) and
   [run_capture.rs](../crates/web/src/run_capture.rs) provide the connected
   construct/oligo/preparation choices, immutable run history and exact-event
   navigation. Further presentation work should expose these existing actions,
   not add another output compiler. HTML/PDF remain declared cloning instructions
   tied to that recipe; native history preserves richer evidence.

Public documentation still leaves manual-edit persistence across all wizard
variants and exact junction-error focus behavior uncertain. It also gives
orientation controls different descriptions in initial fragment preview versus
later fragment editing; do not infer universal restrictions from one dialog.
No new scientific requirement follows from those unknowns. Our explicit
orientation, Apply and fresh-check decisions remain deliberate.

Projects, folder administration, ELN, worklists, inventory and registration are
excluded. The useful pattern is source selection → explicit preparation →
junction/product/primer review → named outputs with ancestry. This note adds no
browser acceptance gate and does not change the user's deferred GUI-review plan.
