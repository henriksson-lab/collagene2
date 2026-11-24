# Preparing ordered overlap fragments

The **Assemble** inspector action opens the ordered overlap workflow. PCR, Assemble and Plan are directly accessible beside Features, RE sites and Digest. The overlap workflow uses original open documents as source roles and the active document as an independently chosen target. Each ordered row explicitly chooses one preparation:

- Whole linear DNA: use the complete source. Circular DNA requires a declared preparation; the stored origin is not treated as a cut.
- Restriction digest fragment: state enzyme names and the 1-based fragment number. The compiler records the digest operation and zero-based core index. Optional explicit end repair can fill 5′ overhangs, remove 3′ overhangs or blunt ends; the default retains the actual cohesive state, and repair is a separate typed operation. Cleavage is an explicit model assumption, not experimental activity or reagent/methylation permission.
- Exact PCR: supply forward/reverse annealing domains separately from their 5′ extensions, all written 5′→3′. Exact binding and competing products are handled by the shared core verifier, not inferred from names or imported primer metadata.

The reverse-complement option acts after preparation. Moving rows changes assembly order. Enter one overlap length per junction, including the closing junction for circular products, and select Gibson, HiFi or In-Fusion advisory guidance. No primer-design engine or external scientific process is called.

Compute method evidence sends original sources, explicit recipe and independent target to the worker. The core reconstructs prepared fragments and checks their oriented end homology. The UI displays prepared reference DNA and actual upstream/downstream junction words; full evidence download retains the compiled preparation graph, source bindings, findings and target assessment. Preparation, sequence consistency and published length guidance are separate findings and do not establish a unique mixture product or experimental success.

Review method candidate recompiles the same source-bound recipe into a typed plan, retaining original source inputs and explicit preparation steps. It does not promote cached prepared references to source material. Any source selection/content, row option/order, target or workspace epoch change invalidates displayed evidence; cancellation closes the private worker. Native plan verification and retained source-event history use the existing candidate workflow.

The desktop form supports up to12 rows. Focused core/MCP and WASM checks run in the background during implementation. Browser checks remain deferred until the full in-scope implementation is complete. Mobile work remains deferred.

A [small PCR + EcoRV example](../fixtures/prepared-overlap-demo.md) provides independently specified fragments, a 15-base junction and correct/wrong targets.
