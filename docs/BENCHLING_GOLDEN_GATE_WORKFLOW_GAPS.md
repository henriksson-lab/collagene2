# Golden Gate preparation and review: remaining workflow gaps

Reviewed 2026-11-25. This is an official-documentation and repository-source audit, not a runtime acceptance report. No local browser, compiler or test was run. It updates the repository assessment in [the earlier restriction/Golden Gate note](BENCHLING_RESTRICTION_GOLDEN_GATE.md): several of that note's general assembly and digest-selection gaps now have source implementations.

## What Benchling documents

**Current multi-fragment workflow.** Golden Gate offers existing Type IIS cut sites or PCR-added sites/overhangs. PCR settings include pre-recognition bases/length, binding-region length/Tm and pair Tm difference. Users add source fragments to ordered bins, review regions/direction and choose production per bin. Saved preferred primers are checked against amplification, site introduction and configured parameters; unsuitable choices can be replaced with warnings. The construct review exposes sequence/map views and a primer table. Finalization offers constructs plus optional primers and fragments. These are documented product interactions, not evidence of the underlying algorithms or chemistry guarantees. [Combinatorial assembly documentation, Golden Gate, preferred primers, construct review and finalization sections](https://help.benchling.com/hc/en-us/articles/39922419457165-How-to-use-the-combinatorial-cloning-and-concatenation-assembly-tools).

**Single-construct selection.** The assembly wizard places a retained backbone and insert into named slots. The documented backbone gesture selects between cut sites, inverts the selection and sets the fragment. The insert selects its cut-bounded interval; orientation can be changed before reviewing the preview and naming the output. Its Golden Gate section references this workflow with Type IIS enzymes. The exact modifier-key gesture differs from the teaching guide, so copying a specific key combination is unnecessary. [Assembly wizard instructions](https://help.benchling.com/hc/en-us/articles/39656605989901-Create-assemblies-with-the-assembly-wizard).

**PCR-added sites are an essential branch.** The worked BsaI example starts with sources lacking the needed sites. Users retain the desired backbone region, set the insert, and review primers that introduce directional BsaI sites before finalizing. Thus, selecting the desired DNA is not equivalent to finding usable restriction sites in the original source. [Molecular Cloning Methods, Golden Gate example](https://help.benchling.com/hc/en-us/articles/9684255457805-Molecular-Cloning-Methods).

**Historical context only.** The 2014 walkthrough explicitly distinguishes selecting existing cut sites from selecting arbitrary DNA and generating site-bearing primers. It also shows spacer handling, primer review/export and a printable summary. Its images and mutation-editing behavior are historical examples, not a current specification we should promise to reproduce. [Golden Gate Assembly Wizard, including public setup/review screenshots](https://www.benchling.com/blog/golden-gate-assembly-wizard).

The current guide embeds fragment-editing and construct-review GIFs beside the corresponding instructions. They support the documented connection between source selection and assembly review; this audit does not infer undocumented interaction timing from those animations. Projects, worklists, registration tables and ELN integration are excluded from our scope.

## What our source currently connects

| Source area | Implemented foundation | Golden Gate-specific missing connection |
|---|---|---|
| `web/src/assembly_workspace.rs` | Persistent pinned original sources, named ordered slots, selection/whole/complement capture, explicit direction and intended-reference ribbon. | No Golden Gate method in this dock; selecting a source region does not create a Type IIS PCR preparation. |
| `web/src/digest.rs` and `assembly_workspace.rs::set_digest_fragment` | Two actual cut boundaries, circular complementary path and exact catalog-pinned digest-fragment handoff. Retains both strand intervals for sticky-end reversal. | Existing-site Golden Gate should reuse this path and restrict/check the intended enzyme geometry; it should not reinterpret a recognition span as a cut or silently use unrelated enzymes. |
| `web/src/cloning_methods.rs` | Exact PCR domains/tails, saved-primer reuse, explicit proposal/application, ordered preparation and ligation, physical end/junction evidence. | Automatic primer suggestion is for overlap homology. There is no Golden Gate clamp/site/spacer/fusion-overhang planner or PCR-followed-by-digest preparation row. |
| `web/src/golden_gate.rs` | BsaI/BsmBI/BbsI/SapI complete-digest pool enumeration, candidate paths, recut warnings, alternatives and bounded-search status. | Selects whole source documents rather than the dock's pinned prepared slots; no PCR-site addition, region-to-primer route or central product/primer review. |
| `web/src/assembly_product_view.rs`, `assembly_finalize.rs` | Fresh reconstructed-product view, junction inspection and explicit construct/prepared-fragment/primer output choices with source history. | The dock callback path currently serves overlap/ligation preparation evidence. Golden Gate candidates need to enter this same output path with their full preparation graph. |

These are static source findings. Their presence does not establish that a complete Golden Gate workflow has passed runtime acceptance.

## Concrete next implementation

1. **Add Golden Gate to the existing persistent draft.** Keep method, enzyme pin, per-slot production choice and exact original source events together. Retain the current whole-template pool enumerator as a separate advanced exploration tool. Changing method, source region, direction, enzyme or primer inputs must clear computed evidence.
2. **Offer two explicit production choices per slot.** “Use existing Type IIS sites” should set an actual reconstructed digest fragment. “Add sites by PCR” should retain the chosen source interval and propose a named primer pair. Manual primers remain editable; no proposal overwrites them until Apply.
3. **Implement one shared scientific preparation compiler.** Emit actual original-source PCR, then the selected enzyme digest and chosen fragment, followed by assembly. Extend the shared recipe layer rather than concatenate tail strings in the UI or open a PCR preview and treat that document's metadata as physical authority. This compiler/planner is missing; the UI cannot truthfully fill it in with existing overlap suggestions.
4. **Make the proposed oligo inspectable.** Show annealing DNA separately from clamp, recognition motif, cleavage spacer and fusion-end bases. Use the pinned enzyme's actual cut offsets and overhang length: the existing SapI profile has a three-base overhang, so a universal four-base field is incorrect. Saved primers require fresh checks for the intended amplicon and resulting digested ends.
5. **Review junctions before committing.** Show upstream/right and downstream/left strand ends, retained source DNA and added junction bases in one convention. Report internal sites, recreated recognition sites, duplicate or self-compatible ends and bounded alternative paths. Sequence compatibility and a finite search receipt must remain distinct from predicted reaction efficiency or uniqueness in a real mixture.
6. **Connect actual results to the central product view.** Keep intended DNA separate from reconstructed output. Selecting a junction should expose its real product DNA and contributing original intervals; primer review should refer to the same frozen recipe. Final output creation must freshly reconstruct that recipe and preserve all original source events, including hidden preparation documents.

These choices are our adaptations, not claims about Benchling internals. We deliberately require explicit primer application, preserve orientation decisions and expose uncertainty instead of adopting automatic replacement, automatic direction changes or a single success badge.

## Delivery acceptance examples

- Existing-site circular backbone plus insert: choose two real cuts, retain the complementary backbone path, reverse a sticky fragment, and verify against its actual bottom-strand-derived top DNA.
- Sources without the chosen sites: select regions, propose and explicitly apply directional site-bearing primers, reconstruct PCR and digestion, and inspect the resulting fusion ends before assembly.
- Mixed preparation: an existing-site backbone and PCR-prepared insert share one pinned enzyme/recipe without using the intended target as a source.
- Failure diagnosis: internal extra cut, incompatible ends, unsuitable saved primer, wrong independent target and search truncation remain distinct, visible outcomes.
- Saved/reopened draft: source pins and declared options survive; computed positive status does not. Recompute before creating outputs. Original source/primer/preparation history remains retained regardless of which documents are shown.

Combinatorial bins with multiple alternatives, spacer mutation tools and advanced oligo optimization can follow this complete single-design route. They are not prerequisites for the source-selection → preparation → end review → product/primer review → outputs connection.
