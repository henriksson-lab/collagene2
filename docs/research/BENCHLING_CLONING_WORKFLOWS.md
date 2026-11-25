# Benchling TOPO/Gateway follow-up

Reviewed **25 November 2026**. This is a narrow update to [Methods follow-up](../BENCHLING_METHODS_FOLLOWUP.md), [workflow decisions](../BENCHLING_WORKFLOW_DECISIONS.md) and [clone verification gaps](../BENCHLING_CLONE_VERIFICATION_GAPS.md). Existing selection, assembly, primer/output review and read-verification findings are not repeated here. No local application/browser check or compilation was used.

## Public evidence and limits

The official [combinatorial cloning guide](https://help.benchling.com/hc/en-us/articles/39922419457165-How-to-use-the-combinatorial-cloning-and-concatenation-assembly-tools), updated **12 December 2025**, lists Gateway under Homology when source fragments already contain ATT regions. Its flow remains settings → bins → fragments → constructs → sequence/map/primer review → finalize outputs. This supplies no dedicated BP/LR site compatibility, reciprocal-product or scar-review procedure.

The [assembly wizard guide](https://help.benchling.com/hc/en-us/articles/39656605989901-Create-assemblies-with-the-assembly-wizard), updated **30 September 2025**, documents digest/ligation, Gibson and Golden Gate source-selection and preview interactions. The [cloning overview](https://www.benchling.com/blog/cloning), published **14 April 2023**, mentions TOPO and Gateway but is educational material, not an application walkthrough or mechanistic validation contract.

Additional official-domain searches for TOPO, directional TOPO, topoisomerase, CACC and Gateway BP/LR did **not** locate a dedicated current procedure. These searches add no supported method-specific interaction beyond the previous audit. This is an evidence limit, not proof that Benchling lacks the feature or tenant-only guidance.

## Current implementation implications — our decisions

- **TOPO delivery has advanced since Methods follow-up.** [Current preparation controls](../TOPO_WEB_WORKFLOW.md) now expose explicit TA/directional CACC modes, actual vector strands, literal primers and declared activation through the existing native orientation/component/source-call review. The older recommendation that directional controls are missing is superseded. Actual graph composition, accounted flap processing, a proven graph-to-duplex bridge and durable graph-only handoff remain separate work. A withheld contour must remain withheld; this increment does not complete F15.
- **Gateway already has a more explicit reference contract.** The [BP/LR workflow](../RECOMBINATION_UI.md) uses exact source events, supported profiles, reciprocal components, independent targets and advisory plans. Retain that flow rather than replacing it with generic homology because of Benchling's label. Broader profiles/MultiSite and physical/repair biology remain distinct scope.
- **Product review needs connection, not another engine.** Reuse existing sequence/primer review, named outputs, bulk release, protocols and construct-owned read checks. Later user review should assess navigation from preparation diagnostics to chosen components and outputs. Public documentation does not establish physical completeness, experimental confirmation or permission to hide unresolved graph outcomes.

These are implementation recommendations, not claims that Benchling models the same chemistry. No new gap is inferred solely from absent public documentation.
