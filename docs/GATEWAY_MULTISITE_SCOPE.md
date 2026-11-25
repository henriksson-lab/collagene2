# Gateway F15 scope and next coherent deliverable

Historical baseline audit,2026-11-25. Later implementation increments are recorded below and in [the connected workflow contract](MULTISITE_CONNECTED_WORKFLOW.md); the original baseline is retained as scope history. No application/browser acceptance is claimed.

## What is actually implemented

`crates/recombination/src/lib.rs` is a bounded **two-source, two-site reciprocal reference exchange**, not full Gateway physical simulation. `Family` B/P/L/R and BP/LR are real typed distinctions. Each full-site variant is exact and source-pinned; the profile validator admits only subtype1/2 and explicitly checks their seven-base overlaps (:101–140). Source preflight requires exactly two molecules (:388). `recombine` groups one subtype1 and one subtype2 hit per source (:889–909), then builds a conservative reference port graph for those four selected sites. It retains possible sites/unsupported orientation/limits and all components, not only a target-like clone.

The embedded eight-context profile is specifically pDONR221 P232, pDEST14 R125, pENTR-D-TOPO L100 and manual attB25 primer contexts. Custom core/MCP snapshots exist, but this is not an arbitrary universal attachment-site matcher. Modified full arms, att1.1/2.1, additional numbered sites, MultiSite and intramolecular/network reactions are outside the current model. Derived BP/LR hybrid arms are not magically padded to a standard profile length; the README correctly distinguishes generated hybrids from registered source contexts.

`core::recombination` already projects conservative annotations, returns every reciprocal component and pins profile/source identities. `Operation::RecombineReference` carries literal request + explicit candidate/component selection, reconstructs freshly and remains Incomplete; its reference-only ancestry cannot acquire physical authority downstream (`verification.rs:2265–2300`, then sticky reference-only flags). The web panel already has two-source BP/LR setup, explicit profile opt-in, component selection, independent target check, Plan handoff and exact preferred source events (`web/src/recombination.rs:65–110,317–377,390–445`). There is no reason to rebuild that existing panel or call it missing Gateway support wholesale.

## What remains physically unimplemented

Current overlap source maps are co-sourced reference provenance, not the atom/strand provenance of recombinase cleavage, heteroduplex exchange, Holliday-junction resolution or ligation. They do not preserve positional modifications, real existing nicks/adducts, strand topology, enzyme complexes or complete reaction mixtures. The manufacturer describes conservative recombination without DNA synthesis and a core strand-exchange mechanism, but that does not identify every inherited strand atom or establish kit efficiency. [Gateway manual, MAN0000282](https://documents.thermofisher.com/TFS-Assets/LSG/manuals/gatewayman.pdf).

Therefore no fresh TrackedDuplex should be manufactured by reverse-complementing the chosen reference product. Nor should lack of synthesis be treated as license to arbitrarily pick one source for the seven-base co-sourced overlap. Full physical Gateway composition requires independently established strand crossover geometry and source-call conservation, rather than relabelling the existing report. Keep that work explicitly in F15.

## Recommended next product deliverable

Implement a **connected three-fragment MultiSite Gateway LR design**, taking existing entry constructs and a destination vector through ordered source selection, att-site compatibility review, explicit product choice, named outputs, exact ancestry and protocol/order exports. This produces an everyday promoter–ORF–terminator construct workflow rather than a profile-only maintenance patch. Architect the new engine for bounded2–4 inserts, while first pinning and delivering the completely documented three-fragment route as its explicit profile. Do not claim Pro/ProPlus kit equivalence from one historical scheme.

Primary evidence is explicit: the Three-Fragment manual describes L4–5′element–R1 + L1–gene–L2 + R2–3′element–L3 with destination R4–R3, producing B4–5′element–B1–gene–B2–3′element–B3 (printed quick-reference p.v; PDFp4). The same manual supplies donor identities, primer-site strings and map references. [MAN0000299 three-fragment manual](https://documents.thermofisher.com/TFS-Assets/LSG/manuals/multisite_gateway_man.pdf).

Do not conflate that route with newer Pro/ProPlus designs. Manufacturer support says the three-fragment construction kit and Pro families use different vector/site arrangements; it also distinguishes LR ClonaseII Plus from ordinary LR formulations. [Official support comparison](https://www.thermofisher.com/uk/en/home/technical-resources/technical-reference-library/cloning-technical-support-center/recombination-based-cloning/recombination-based-cloning-getting-started.html). Keep formulation as explicit documentary context, not a predicted success factor. The manufacturer separately describes current MultiSite as assembling2–4 fragments. [MultiSite technology overview](https://www.thermofisher.com/us/en/home/life-science/cloning/gateway-cloning/multisite-gateway-technology.html).

### Native algorithm/profile work

Use a NEW additive multi-site request/report/profile model, keeping old two-site structs and their native layouts stable. Subtype cannot simply become3/4 in the old validator: matching groups and graph construction hard-code two sites per each of two sources. A profile needs explicit site specificity IDs, family and canonical orientation, sourced full contexts, crossover overlap boundaries, allowed cognate pairs and documentary kit/version. Treat reversed sites as orientation of a typed specificity, never by stripping an “r” suffix heuristically. Preserve variant distinctions even when an overlap is equal.

For a three-fragment reaction there are four parent molecules and four explicit compatible pairings: R4/L4, R1/L1, L2/R2 and L3/R3. Enumerate source-bound site placements under global work/state/materialization budgets, build one whole port graph, and reconnect all selected pairs within that same graph. Discover all resulting components and preserve every original reference base exactly once. Do NOT run the old two-parent cassette engine repeatedly and silently discard intermediate reciprocal products to force the desired final construct. Do not assume every topology/orientation gives exactly two components.

Reuse the independently written source-scan, interval, orientation and segment-audit concepts from the old engine through narrow shared helpers where useful. Independently test the new general graph; the old two-site oracle is not sufficient evidence for four-source behavior. Include missing/ambiguous/mutated sites, duplicate alternatives, noncognate overlaps, reordered/reverse parents, origin-spanning contexts and atomic bounded partial search. No candidate may be chosen by matching the expected target.

Before shipping the named profile, obtain exact full-site contexts from official vector sequences/maps, with raw source hashes and attribution. Manual diagrams containing abbreviated arms or N placeholders are not full-site sequences. The retrieved Pro manual URL returned403 during this audit, so exact broader Pro schemes remain unverified here; this is why the fully retrievable MAN0000299 route is the concrete first profile. Never synthesize supposedly missing arms or disguise a source-derived hybrid as an official full variant.

### Connected design and output work

Use existing ordered source-slot/region/orientation contracts, independent intent builder, source-event capture and actual-product/junction review. An ordinary user should choose destination + ordered entry fragments, see local site compatibility failures and preview the scar-bearing intended reference before reaction reconstruction. If the user supplies an independent complete target, retain it. Otherwise freeze intent from explicitly selected source regions/order plus source-pinned literal expected junction sequences BEFORE running the reaction graph; never copy the reaction product into the target field and call that independent verification.

A future attB-primer/BP-entry preparation route should reuse the existing shared primer editor, actual PCR compiler and order exports. Manufacturer examples include extra5′primer sequence outside the attachment core; keep clamps/full attachment context/annealing domain visibly distinct and do not invent an attL-tailed PCR shortcut. A three-entry LR deliverable is already substantial without pretending those entry clones were experimentally made or synthesizing their long attachment arms.

Add one new typed plan operation selecting an exact multi-site candidate/component from the SAME fresh enumeration. Keep reference-only scope explicit until physical strand simulation exists. Reuse current run/source/target binding and immutable historical artifacts, shared HTML/PDF rendering and native archive output rather than inventing another organizational/project object. The current combinatorial design Method/Recipe types have no Gateway variant; add the method through the same compiler so selected alternative bins can later produce multiple reference constructs consistently, rather than leave MultiSite only as another detached tool panel.

### Completion criterion for this deliverable

A user opens four actual GenBank constructs, selects destination plus promoter/ORF/terminator entries, chooses the explicit supported profile, reviews every candidate/component and scar/junction, freezes/uses independent intent, then saves a named reference construct with exact producing sources and printable declared steps. Wrong targets fail; incomplete/possible evidence remains explicit; all reciprocal components remain available. Output remains a reference design, not a chemically verified Gateway clone. Physical Gateway strand-state composition, other site families/formulations,2/4-fragment profiles and broader batch/remaining methods remain on the original goal.


## Primary source acquisition follow-up

The [source acquisition record](MULTISITE_GATEWAY_SOURCES.md) now pins manufacturer
R3/R4, P3, B3/B4 and fully printed L4/R1 contexts with raw hashes and exact ranges.
The original 4107 bp and VectorII 4555 bp destinations remain distinct sources.
A published 94-base putative L3 flank also occurs verbatim in unrecombined P3,
so it cannot define family-specific L3 alone. Resolve its distinguishing context
from a real entry sequence or explicit source-derived BP provenance before using
it as a closed built-in LR profile. No missing arms or insert bases were invented.


## Whole-network reference implementation increment

The separate [whole-network model and Plan bridge](MULTISITE_REFERENCE.md) now
implements the bounded multi-parent graph, explicit candidate/component selection,
source annotation evidence, independent-target checks, binary worker/MCP and exact
historical attachments. It does not yet complete the destination/entry intent and
combinatorial interaction criteria above. Its explicit profile import is not a
bundled universally applicable Gateway kit profile.

The actual original-lab L3 context is96 bases and remains nested in the longer P3
context. Classification must therefore consider all registered full contexts
before filtering LR families. Adding an arbitrary insert base to make L3 unique
would be incorrect. The new finite model uses maximal certain physical contexts
and withholds conflicting maximal families. Complete source fixtures remain outside
the repository while redistribution terms are unresolved; synthetic fixtures
independently establish reference graph bookkeeping, not biological kit validity.


## Connected workflow and finite preset increment

The [ordinary/library compiler and UI increment](MULTISITE_CONNECTED_WORKFLOW.md)
now reuses ordered original contributions, optional junction intent, exact native
events and the same named-output/history routes. Whole parents feed one network;
source regions define a target before enumeration. Explicit component pins include
source-event and contribution/target identity; neither a short motif nor the
target chooses a family or candidate. A single construct does not require creating
a library. The exact short nine-context source-edition profile is now bundled as
functional scientific facts, with source hashes and independently written curator
code; complete plasmid records remain unbundled. The next broader work includes
physical chemistry, additional biological profile/preparation coverage and full
GUI acceptance, rather than rebuilding this source/intent/output compiler.
