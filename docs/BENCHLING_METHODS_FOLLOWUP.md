# Remaining cloning-method workflow review

Reviewed 25 November 2026 against official public Benchling sources and current
repository contracts. This supplements [workflow decisions](BENCHLING_WORKFLOW_DECISIONS.md)
and [cloning interaction design](CLONING_INTERACTION_DESIGN.md). No application
browser checks, compiler runs, authenticated sessions or third-party code were
used. Search results that mention a technique are distinguished from documented
application controls; failure to locate a procedure does not establish that the
product lacks the capability.

## Documented behavior and explicit unknowns

| Method | Official evidence | What remains unestablished |
| --- | --- | --- |
| Gateway/recombination | The current combinatorial guide lists Gateway under Homology when fragments already contain ATT regions. The ordinary assembly settings and bins apply. Its Homology ambiguity preference selects a larger or smaller first-bin fragment. [Combinatorial guide, updated 12 December 2025](https://help.benchling.com/hc/en-us/articles/39922419457165-How-to-use-the-combinatorial-cloning-and-concatenation-assembly-tools). | Dedicated BP/LR setup, compatible site-family rules, reciprocal-product selection and method-specific recovery are not specified. The educational [cloning overview](https://www.benchling.com/blog/cloning) is not an application procedure or a scientific correctness contract. |
| TOPO | The same educational overview discusses TOPO as a cloning technique. This review did not locate a dedicated current public help procedure for its setup and preview. [Official overview](https://www.benchling.com/blog/cloning). | Vector activation, prepared strand geometry, directional choices, transfer selection, retained nicks and acceptance policy cannot be inferred from the technique name or article. No kit assumptions are adopted. |
| TA/GC and end repair | The ordinary wizard documents source selections, orientation correction, preview and creation for digest/ligation. Copy/paste cloning uses cut-bounded insert/destination selections and documents automatic reversal and sticky-end handling. [Wizard, updated 30 September 2025](https://help.benchling.com/hc/en-us/articles/39656605989901-Create-assemblies-with-the-assembly-wizard); [copy/paste, updated 22 September 2025](https://help.benchling.com/hc/en-us/articles/9684297597453-Clone-by-copy-and-paste). | No dedicated public TA/GC tailing, fill-in, overhang removal, phosphorylation or dephosphorylation setup/recovery was established. General sequence editing and sticky-end insertion do not establish those chemical transformations. |
| Oligo concatenation | The combinatorial guide explicitly distinguishes sequence concatenation, without additional joining logic, from cloning; DNA/RNA oligos are supported inputs/outputs. [Combinatorial guide](https://help.benchling.com/hc/en-us/articles/39922419457165-How-to-use-the-combinatorial-cloning-and-concatenation-assembly-tools). | Concatenation is not documented evidence of annealing, ligation, strand synthesis or multi-oligo physical assembly. |
| Oligo duplex setup | A distinct guide creates a composite entity by selecting forward/reverse oligos and their offset, with optional components/attachments. It subsequently supports editing the displayed connections. [Duplex guide, updated 1 October 2025](https://help.benchling.com/hc/en-us/articles/34874363508365-Create-oligo-duplex-and-oligo-conjugate-entities). April 2025 release notes explicitly document a local error for invalid duplex offsets. [Release notes, updated 7 August 2025](https://help.benchling.com/hc/en-us/articles/36130712944909-Release-Notes-Volume-4-2025). | The guide does not specify complete complementary-offset enumeration, mismatch handling, thermodynamics, cloning-ready terminal chemistry or downstream ligation authority. Composite chemical connections must not be confused with a validated cloning reaction. |

These findings come from written documentation; linked animations were not
inspected frame by frame. Registry schemas, component registration, notebook
handoff and chemical-conjugate administration are outside our scope.

## Recommendations against the existing implementation

These are our adaptations, not additional claims about Benchling:

1. **Keep explicit preparation beside source selection.** Retain the current
   fill/remove/blunt and terminal-treatment choices, showing actual before/after
   strand ends and added/removed bases. [End repair/annealing](END_REPAIR.md) and
   [tailing](TAILING.md) already implement these transformations and workers;
   another physical engine is unnecessary. Errors should point to the affected
   source/end or parameter. Do not adopt undocumented automatic repair.
2. **Extend the existing source/recipe handoff where coverage is missing.** Current
   [cloning methods](../crates/web/src/cloning_methods.rs) support explicit TA/GC
   graphs, but disable that route for stable library source roles. Native
   [combinatorial design](../crates/core/src/combinatorial_design.rs) does not yet
   represent tailing or two-input annealing as a library preparation. Future
   integration should reuse existing typed operations and exact original events;
   copying an intermediate reference into a slot must not invent physical stock.
3. **Preserve the existing TOPO inspection contract.** The current
   [web workflow](TOPO_WEB_WORKFLOW.md) exposes one declared TA geometry; broader
   directional preparation controls remain a real delivery gap. Reuse the
   [general core adapter](TOPO_ADAPTER.md), explicit component/source-call choice,
   independent target and acceptance policy. Public Benchling evidence supplies
   no justification for hiding transfer outcomes or defaulting away retained nicks.
4. **Keep annealing distinct from reference concatenation.** Reuse the existing
   [annealing candidates](END_REPAIR.md#oligo-annealing), actual antiparallel
   overlap and per-strand termini. Present offset/overhang previews and local
   invalid-placement messages; preserve alternatives and explicit selection.
   Multi-oligo physical assembly needs its own supported composition rules;
   Benchling's concatenation feature does not establish them.
5. **Do not relabel reference recombination as full Gateway chemistry.** Our
   [BP/LR workflow](RECOMBINATION_UI.md) already has explicit profiles, exact
   sources, reciprocal components and independent-target checking. Broader
   profiles and physical chemistry remain separate work documented in the
   [adapter scope](RECOMBINATION_ADAPTER.md). An ambiguity preference must not
   silently choose the biological target.

Continue source → declared preparation → actual product/diagnostic review → shared
named outputs and HTML/PDF steps. These findings justify connecting existing
surfaces and exposing their limits, not copying enterprise structure or inventing
unsupported Benchling methods.
